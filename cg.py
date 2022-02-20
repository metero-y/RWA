# -*- coding: utf-8 -*-
# Authors:   张䶮
"""
实现列生成
"""
import collections
import copy

import networkx as nx
from gurobipy import *

#还需要添加bound and price
'''

'''

# 判断任务请求在图G中是否存在路径, 若存在返回对应路径
def has_path(G, task):
    """
    :param           G: 图
    :param        task: 任务请求
    :return:  all_path: 所有最短路径
    """
    try:
        path = nx.all_shortest_paths(G, source=task[0], target=task[1])
        all_path = []
        for i in path:
            tmp = []
            for j in range(len(i) - 1):
                tmp.append((i[j], i[j + 1]))
            all_path.append(tmp)

    except nx.NetworkXNoPath:
        return False, None
    return True, all_path


def RF(arcs, V, SD, SD_idx,
       W):  # routing formulation that ignores the wavelength continuity constraints and omits integrality requirements
    try:
        m = gurobipy.Model('RF')
        fsdl = m.addVars(arcs, range(len(SD)), lb=0.0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='fsdl')
        dsd = []
        for i in SD_idx:
            tmp = m.addVar(lb=0.0, ub=SD[i], vtype=GRB.CONTINUOUS, name='dsd' + str(i))
            dsd.append(tmp)

        m.addConstrs((fsdl.sum(i, j, '*') <= len(W) for i, j in arcs), name='每条边最多只能满足一个任务')
        m.addConstrs(
            (fsdl.sum(k, '*', j) == fsdl.sum('*', k, j) for j in range(len(SD)) for k in V if
             k != SD_idx[j][0] and k != SD_idx[j][1]), "流守恒")
        m.addConstrs(
            (fsdl.sum(SD_idx[j][0], '*', j) == dsd[j] for j in range(len(SD))), "负载")
        m.addConstrs(
            (fsdl.sum('*', SD_idx[j][1], j) == dsd[j] for j in range(len(SD))), "负载")
        m.addConstrs(
            (fsdl.sum(SD_idx[j][1], '*', j) == 0 for j in range(len(SD))), "防止出现环路")
        m.addConstrs(
            (fsdl.sum('*', SD_idx[j][0], j) == 0 for j in range(len(SD))), "防止出现环路")

        expression = 0
        for i in dsd:
            expression += i
        m.setObjective(expression, GRB.MAXIMIZE)

        m.optimize()
        m.write('RF.lp')
        dic = collections.defaultdict(int)
        for v in m.getVars():
            if v.varName[:4] == 'fsdl':
                num = -1
                acr_x = ''  # 边的入结点
                acr_y = ''  # 边的出结点
                for i in range(len(v.varName)):
                    if num == 2:
                        break
                    if num == 1:
                        acr_y += v.varName[i]
                    if v.varName[i] == ',':
                        num += 1
                    if num == 0:
                        acr_x += v.varName[i]
                    if v.varName[i] == '[':
                        num += 1
                dic[(int(acr_x), int(acr_y[:-1]))] += v.x
        return dic

    except gurobipy.GurobiError as e:
        print('Errorcode ' + str(e.errno) + ": " + str(e))
        return 0

    except AttributeError:
        print('Encountered an attribute error')
        return 0


def has_K_path(G, task, strategy, SD, SD_idx, arcs, V, W):
    """
    :param          G: 图
    :param       task: 任务请求
    :param   strategy: k最短路径的选择策略
    :return:     path: 所有k最短路径
    """
    try:
        Paths = nx.all_simple_paths(G, source=task[0], target=task[1])  # 存储所有简单路径
        num = []  # 存储所有简单路径对应的路径数值
        paths = []  # 列表化Paths
        for i in Paths:
            num.append(len(i)-1)
            tmp = []
            for j in range(len(i) - 1):
                tmp.append((i[j], i[j + 1]))
            paths.append(tmp)
        order = sorted(range(len(num)), key=lambda k: num[k], reverse=False)  # 所有简单路径按长度升序排列, 存储对应的路径序号
        Path = copy.deepcopy(paths)
        nums = copy.deepcopy(num)
        for i in range(len(order)):
            Path[i] = paths[order[i]]
            nums[i] = num[order[i]]

        if strategy == 1:  # 策略1
            for i in nums:
                if i != nums[1]:
                    break
            if len(nums) - i <= 15:  # 如果除去最短路径不超过15个，就全部接受，否则接受前15个
                all_path = Path[:]
            else:
                all_path = Path[:i + 15]
        elif strategy == 2:  # 策略2
            if len(nums) <= 2 * SD[task]:  # 如果路径数不超过2倍的需求，就全部接受，否则接受前2倍的需求个
                all_path = Path[:]
            else:
                max_k = nums[2 * SD[task] - 1]
                p_l = -1
                length = 0
                for i in nums:
                    if i < max_k:
                        p_l += 1
                    elif i == max_k:
                        length += i
                    else:
                        break
                if p_l + length + 1 == 2 * SD[task]:
                    all_path = Path[:p_l + length + 1]  # 最后一组pi全部需要
                else:
                    pi = Path[p_l + 1:p_l + length + 1]  # 提取最后一组pi
                    load_acr = RF(arcs, V, SD, SD_idx, W)
                    if not load_acr:
                        return False, None
                    load_l = []
                    for i in pi:
                        num = 0
                        for j in i:
                            num += load_acr[j]
                        load_l.append(num)
                    order_l = sorted(range(len(load_l)), key=lambda k: load_l[k],
                                     reverse=False)  # 最后一组路径按负载升序排列, 存储对应的路径序号
                    Path_pi = copy.deepcopy(pi)
                    for i in range(len(order_l)):
                        Path_pi[i] = pi[order_l[i]]
                    all_path = Path[:p_l + 1] + Path_pi[:2 * SD[task] - p_l - 1]
        elif strategy == 3:  # 策略3
            if nums[0] == 1:  # 有一跳路径就只用一跳路径
                if nums[-1] == 1:
                    all_path = Path[:]
                    return True, all_path
                for i in range(len(nums)):
                    if nums[i] > 1:
                        break
                all_path = Path[:i]
                return True, all_path
            if len(nums) <= nums[0] * SD[task]:  # 如果路径数不超过2倍的需求，就全部接受，否则接受前2倍的需求个
                all_path = Path[:]
            else:
                max_k = nums[nums[0] * SD[task] - 1]
                p_l = -1
                length = 0
                for i in nums:
                    if i < max_k:
                        p_l += 1
                    elif i == max_k:
                        length += i
                    else:
                        break
                if p_l + length + 1 == nums[0] * SD[task]:
                    all_path = Path[:p_l + length + 1]  # 最后一组pi全部需要
                else:
                    pi = Path[p_l + 1:p_l + length + 1]  # 提取最后一组pi
                    load_acr = RF(arcs, V, SD, SD_idx, W)
                    if not load_acr:
                        return False, None
                    load_l = []
                    for i in pi:
                        num = 0
                        for j in i:
                            num += load_acr[j]
                        load_l.append(num)
                    order_l = sorted(range(len(load_l)), key=lambda k: load_l[k],
                                     reverse=False)  # 最后一组路径按负载升序排列, 存储对应的路径序号
                    Path_pi = copy.deepcopy(pi)
                    for i in range(len(order_l)):
                        Path_pi[i] = pi[order_l[i]]
                    all_path = Path[:p_l + 1] + Path_pi[:nums[0] * SD[task] - p_l - 1]
        else:
            print('策略有误！')
            return False, None

    except nx.NetworkXNoPath:
        return False, None
    return True, all_path


# 得到连接请求的最短路径, 并按照降序进行排列
def ALL_Shortest_Path(SD_idx, G):
    """
    :param   SD_idx: 任务请求
    :param       G: 图的结构信息
    :return:  Path: 每个任务连接请求的对应所有最短路径， 找不到路径的序对的路径为空
    """
    Path = []  # 存储每个连接请求的最短路径
    for task in SD_idx:
        exist_path, path = has_path(G, task)
        if exist_path:
            Path.append(path)
        else:  # 找不到最短路径
            Path.append(0)
    return Path


def ALL_K_Shortest_Path(SD, SD_idx, G, strategy, arcs, V, W):
    """
    :param     SD_idx: 任务请求
    :param          G: 图的结构信息
    :param   strategy: k最短路径的选择策略
    :return:     Path: 每个任务连接请求的对应的所有k最短路径， 找不到路径的序对的路径为空
    """
    Path = []  # 存储每个连接请求的最短路径
    for task in SD_idx:
        exist_path, path = has_K_path(G, task, strategy, SD, SD_idx, arcs, V, W)
        if exist_path:
            Path.append(path)
        else:  # 找不到最短路径
            Path.append(0)
    return Path


def CG(arcs, V, SD, SD_idx, W, C, a):
    while 1:
        dual = RMP(SD, SD_idx, W, C, a)
        solver = LPP(arcs, V, SD, SD_idx, C, a, dual)
        if not solver:
            break
    return MP(SD, W, C, a)


def CG_plus(arcs, V, SD, SD_idx, W, C, a, P):
    while 1:
        dual = RMP(SD, SD_idx, W, C, a)
        solver = PPP(arcs, V, SD, SD_idx, C, a, dual, P)
        if solver:
            continue
        solver = LPP(arcs, V, SD, SD_idx, C, a, dual)
        if not solver:
            break
    return MP(SD, W, C, a)


def CG_plus_H(arcs, V, SD, SD_idx, W, C, a, P):
    while 1:
        dual = RMP(SD, SD_idx, W, C, a)
        solver = PPP(arcs, V, SD, SD_idx, C, a, dual, P)
        if solver:
            break
    return MP(SD, W, C, a)


def RMP(SD, SD_idx, W, C, a):  # restrict master problem
    try:
        m = gurobipy.Model('RMP')
        zc = m.addVars(len(C), lb=0.0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='zc')
        ysd = m.addVars(len(SD), lb=0.0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='ysd')

        m.addConstr((zc.sum() <= len(W)), name='不能超过可使用波长数')

        for i in range(len(SD)):
            prod = [a[j][SD_idx[i]] for j in range(len(C))]

            m.addConstr(ysd[i] <= (zc.prod(prod)), name='不能超过' + str(SD_idx[i]) + '的总路由数')
        m.addConstrs((ysd[i] <= SD[SD_idx[i]] for i in range(len(SD))), name='不能超过任务数')

        m.setObjective(ysd.sum(), GRB.MAXIMIZE)

        m.optimize()
        m.write('RMP.lp')
        return (m.pi)


    except gurobipy.GurobiError as e:
        print('Errorcode ' + str(e.errno) + ": " + str(e))
        return 0

    except AttributeError:
        print('Encountered an attribute error')
        return 0


def LPP(arcs, V, SD, SD_idx, C, a, dual):  # price problem：link
    try:
        m = gurobipy.Model('LPP')
        asdl = m.addVars(arcs, range(len(SD)), vtype=GRB.BINARY, name='asdl')

        m.addConstrs((asdl.sum(i, j, '*') <= 1 for i, j in arcs), name='每条边最多只能满足一个任务')
        m.addConstrs(
            (asdl.sum(k, '*', j) == asdl.sum('*', k, j) for j in range(len(SD)) for k in V if
             k != SD_idx[j][0] and k != SD_idx[j][1]), "流守恒")
        m.addConstrs(
            (asdl.sum(SD_idx[j][0], '*', j) <= SD[SD_idx[j]] for j in range(len(SD))), "不能超过任务数")
        m.addConstrs(
            (asdl.sum(SD_idx[j][1], '*', j) == 0 for j in range(len(SD))), "防止出现环路")
        m.addConstrs(
            (asdl.sum('*', SD_idx[j][0], j) == 0 for j in range(len(SD))), "防止出现环路")

        expression = 0
        for i in range(len(SD)):
            port_sum = asdl.sum(SD_idx[i][0], '*', i) * dual[i + 1]
            expression += port_sum

        m.setObjective(expression - dual[0], GRB.MAXIMIZE)

        m.optimize()
        m.write('LPP.lp')
        if m.objVal > 0:  # 添加新的策略
            c = []  # 单层策略
            for v in m.getVars():
                if v.x == 1:  # 得到被路由的边
                    num = -1
                    sd_idx = ''  # 任务下标
                    start_node = ''  # 起始结点
                    for i in range(len(v.varName)):
                        if v.varName[i] == ']':
                            break
                        if num == 2:
                            sd_idx += v.varName[i]
                        if v.varName[i] == ',':
                            num += 1
                        if num == 0:
                            start_node += v.varName[i]
                        if v.varName[i] == '[':
                            num += 1
                    if SD_idx[int(sd_idx)][0] == int(start_node):
                        c.append(SD_idx[int(sd_idx)])
            C.append(c)
            dic = collections.Counter(c)
            a.append(dic)
            return 1
        else:
            return 0

    except gurobipy.GurobiError as e:
        print('Errorcode ' + str(e.errno) + ": " + str(e))
        return 0

    except AttributeError:
        print('Encountered an attribute error')
        return 0


def PPP(arcs, V, SD, SD_idx, C, a, dual, P):  # price problem：path
    try:
        m = gurobipy.Model('PPP')
        bsdp = []  # beta^sd_p
        for i in range(len(SD)):
            tmp = m.addVars(len(P[i]), vtype=GRB.BINARY, name='bsdp' + str(SD_idx[i]))
            bsdp.append(tmp)

        for i in arcs:
            expression = 0
            for j in range(len(SD)):
                prod = [1 if i in p else 0 for p in P[j]]
                expression += bsdp[j].prod(prod)
            m.addConstr(expression <= 1, name='每条边' + str(i) + '最多只能满足一个任务')

        for j in range(len(SD)):
            m.addConstr(bsdp[j].sum() <= SD[SD_idx[j]], "不能超过任务数" + str(SD_idx[j]))

        expression = 0
        for i in range(len(SD)):
            port_sum = bsdp[i].sum() * dual[i + 1]
            expression += port_sum
        m.setObjective(expression - dual[0], GRB.MAXIMIZE)

        m.optimize()
        m.write('PPP.lp')
        if m.objVal > 0:  # 添加新的策略
            c = []  # 单层策略
            for v in m.getVars():
                if v.x == 1:  # 得到被路由的路径
                    num = -1
                    sd_f = ''  # 任务起始点
                    sd_e = ''  # 任务终止点
                    for i in range(len(v.varName)):
                        if v.varName[i] == ')':
                            break
                        if num == 2:
                            sd_e += v.varName[i]
                        if v.varName[i] == ',' or v.varName[i] == ' ':
                            num += 1
                        if num == 0:
                            sd_f += v.varName[i]
                        if v.varName[i] == '(':
                            num += 1
                    c.append((int(sd_f), int(sd_e)))
            C.append(c)
            dic = collections.Counter(c)
            a.append(dic)
            return 1
        else:
            return 0

    except gurobipy.GurobiError as e:
        print('Errorcode ' + str(e.errno) + ": " + str(e))
        return 0

    except AttributeError:
        print('Encountered an attribute error')
        return 0


def MP(SD, W, C, a):
    try:
        m = gurobipy.Model('cg')
        zc = m.addVars(len(C), lb=0.0, ub=GRB.INFINITY, vtype=GRB.INTEGER, name='zc')

        m.addConstr((zc.sum() <= len(W)), name='不能超过可使用波长数')
        expression = 0
        for i in SD:
            prod = [a[j][i] for j in range(len(C))]
            expression += zc.prod(prod)
            m.addConstr((zc.prod(prod) <= SD[i]), name='不能超过' + str(i[0]) + ',' + str(i[1]) + '的任务数')

        m.setObjective(expression, GRB.MAXIMIZE)

        m.optimize()
        m.write('cg.lp')
        ans = {}
        for v in m.getVars():
            ans[v.varName] = v.x
        return (ans)

    except gurobipy.GurobiError as e:
        print('Errorcode ' + str(e.errno) + ": " + str(e))
        return 0

    except AttributeError:
        print('Encountered an attribute error')
        return 0


if __name__ == "__main__":
    L = [(0, 1), (1, 3), (0, 2), (2, 3)]  # 边
    arcs = []  # 双向边
    for i in range(len(L)):
        arcs.append(L[i])
        arcs.append((L[i][1], L[i][0]))
    V = range(4)  # 结点
    SD = {(0, 3): 2, (0, 2): 2}  # 任务
    SD_idx = []
    for rq in SD:
        SD_idx.append(rq)

    W = range(2)  # 可使用波长
    G = nx.Graph()
    for i in V:
        G.add_node(V[i])
    for x, y in arcs:  # edges::
        G.add_edges_from([(x, y)])

    # CG
    C = [[]]  # 策略
    dit = collections.Counter(C[0])
    a = [dit]
    ans = CG(arcs, V, SD, SD_idx, W, C, a)
    print(C, a)
    print(ans)

    # CG+与CG+H
    P1 = ALL_Shortest_Path(SD_idx, G)  # CG+&CG+H
    print(P1)
    if 0 in P1:
        print('有任务没有路径')
    else:
        # CG+
        C = [[]]  # 策略
        dit = collections.Counter(C[0])
        a = [dit]
        ans = CG_plus(arcs, V, SD, SD_idx, W, C, a, P1)
        print(C, a)
        print(ans)
        # CG+H
        C = [[]]  # 策略
        dit = collections.Counter(C[0])
        a = [dit]
        ans = CG_plus_H(arcs, V, SD, SD_idx, W, C, a, P1)
        print(C, a)
        print(ans)

    # CG++与CG++H
    strategy = 3 # 1,2,3
    P2 = ALL_K_Shortest_Path(SD, SD_idx, G, strategy, arcs, V, W)  # CG+&CG+H
    print('P2',P2)
    if 0 in P2:
        print('有任务没有路径')
    else:
        # CG+
        C = [[]]  # 策略
        dit = collections.Counter(C[0])
        a = [dit]
        ans = CG_plus(arcs, V, SD, SD_idx, W, C, a, P2)
        print(C, a)
        print(ans)
        # CG+H
        C = [[]]  # 策略
        dit = collections.Counter(C[0])
        a = [dit]
        ans = CG_plus_H(arcs, V, SD, SD_idx, W, C, a, P2)
        print(C, a)
        print(ans)
