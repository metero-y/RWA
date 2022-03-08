# -*- coding: utf-8 -*-
# Authors:   张䶮
"""
启发式算法：最小阻塞
"""
from cmath import inf
import collections
import copy
import heapq
import networkx as nx
import random
import numpy as np
from gurobipy import *

from sqlalchemy import true


def solve(G, sp_len, al_simple_path, SD, W, residue_rq, bandwidth):
    order = []
    for rq in sp_len:
        if SD[rq] == 0:
            continue
        # 起始点出deg，终止点入deg，最短路径值,越小越优先
        tmp = G[0].out_degree(rq[0])*len(W) * \
            G[0].in_degree(rq[1])*len(W) * sp_len[rq][0]*al_simple_path[rq]
        if tmp != 0:  # 为0的时候自动忽略，因为这本身就不可行
            order.append((tmp, rq))
        else:
            SD[rq] = 0
    heapq.heapify(order)

    C = [[] for _ in W]  # 总策略
    C_len = [[] for _ in W]  # 总策略对应路径长度
    while order:  # 边赋权还没实现
        require = heapq.heappop(order)
        lightpath = sp_len[require[1]][0]
        wavelength = sp_len[require[1]][1]
        C[wavelength].append(require[1])
        C_len[wavelength].append(lightpath)
        for i in range(len(lightpath)-1):  # 更新图
            G[wavelength].remove_edges_from(
                [(lightpath[i], lightpath[i+1])])
            bandwidth[(lightpath[i], lightpath[i+1])] -= 1  # 更新带宽
        for arc in arcs:
            if arc[0] == require[1][0] or arc[1] == require[1][1]:
                residue_rq[arc] -= 1  # 更新该边的起点或终点是一个未完成任务的个数
        for w in W:  # 更新权值
            for x, y in list(G[w].edges):
                # 权重：该边的起点或终点是一个未完成任务的个数/该边的空余带宽
                G[w][x][y]['weight'] = residue_rq[(x, y)]/bandwidth[(x, y)]

        SD[require[1]] = SD[require[1]]-1

        order = []  # 更新order
        for rq in sp_len:
            if SD[rq] == 0:
                continue
            # 起始点出deg，终止点入deg，最短路径值,越小越优先
            outdeg = 0
            indeg = 0
            min_val = float(inf)
            al_simple_path[rq] = 0
            for w in range(len(W)):
                outdeg += G[w].out_degree(rq[0])
                indeg += G[w].in_degree(rq[1])
                try:
                    path = nx.shortest_path(
                        G[w], rq[0], rq[1], weight='weight')
                    paths = nx.all_simple_paths(
                        G[w], rq[0], rq[1], len(path))
                    for _ in paths:
                        al_simple_path[rq] += 1
                except:
                    continue
                now_val = len(path)-1
                if now_val < min_val:
                    min_val = now_val
                    wavelength = w
                    lightpath = path
            if min_val < float(inf):
                sp_len[rq] = [lightpath, wavelength]
            else:
                SD[rq] = 0
                continue
            tmp = outdeg*indeg*(len(sp_len[rq][0])-1)*al_simple_path[rq]
            if tmp != 0:  # 为0的时候自动忽略，因为这本身就不可行
                order.append((tmp, rq))
            else:
                SD[rq] = 0
        heapq.heapify(order)

    return C, C_len


def opt(W, K, arcs, V, assignment, C_len):
    try:
        m = gurobipy.Model('opt')
        xwke = m.addVars(W, range(len(K)), arcs, vtype=GRB.BINARY, name='xkew')
        ywk = m.addVars(W, range(len(K)), vtype=GRB.BINARY, name='ywk')

        m.addConstrs((ywk.sum('*', k) <= 1 for k in range(len(K))),
                     name='每个任务最多只分配一个波长')
        m.addConstrs((xwke.sum(w, '*', i, j) <= 1 for w in W for i,
                     j in arcs), name='同一波长链路只能有一个信息')
        m.addConstrs((xwke[w, k, i, j] <= ywk[w, k] for w in W for k in range(len(K)) for i, j in arcs),
                     name='不在任务波长下的任务链路为0')
        m.addConstrs(
            (xwke.sum(w, k, '*', j) == xwke.sum(w, k, j, '*') for w in W for k in range(len(K)) for j in V if
             j != K[k][0] and j != K[k][1]), "流守恒")
        m.addConstrs((xwke.sum(w, k, K[k][0], '*') - xwke.sum(w, k, '*', K[k][0]) == ywk[w, k] for w in W for k in
                      range(len(K))), "发源地")
        m.addConstrs((xwke.sum(w, k, K[k][1], '*') - xwke.sum(w, k, '*', K[k][1]) == -ywk[w, k] for w in W for k in
                      range(len(K))), "接受地")

        m.setObjective(ywk.sum(), GRB.MAXIMIZE)

        for w in W:
            for k in range(len(K)):
                ywk[w, k].Start = 0
                for x, y in arcs:
                    xwke[w, k, x, y].Start = 0

        for w in W:
            for k in range(len(assignment[w])):
                ywk[w, assignment[w][k]].Start = 1
                for i in range(len(C_len[w][k])-1):
                    xwke[w, assignment[w][k], C_len[w][k]
                         [i], C_len[w][k][i+1]].Start = 1

        m.optimize()
        print(m.objVal)
    except gurobipy.GurobiError as e:
        print('Errorcode ' + str(e.errno) + ": " + str(e))
        return 0

    except AttributeError:
        print('Encountered an attribute error')
        return


if __name__ == "__main__":
    L = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 7), (2, 5), (3, 4), (3, 8), (4, 5), (4, 6), (5, 10), (5, 12), (6, 7),
         (7, 9), (8, 11), (8, 13), (9, 10), (9, 11), (9, 13), (11, 12), (12, 13)]  # 边
    arcs = []  # 双向边
    bandwidth = collections.defaultdict(int)  # 每一条边的空余带宽
    for i in range(len(L)):
        arcs.append(L[i])
        arcs.append((L[i][1], L[i][0]))
        bandwidth[arcs[-1]] = 30
        bandwidth[arcs[-2]] = 30
    V = range(14)  # 结点

    random.seed(1015)  # 设置随机种子random.seed(1015)
    K = []
    residue_rq = collections.defaultdict(int)  # 该边的起点或终点是一个任务的个数
    for _ in range(436):
        task = random.sample(range(14), 2)  # 结点不可重复
        K.append((task[0], task[1]))
        for arc in arcs:
            if arc[0] == task[0] or arc[1] == task[1]:
                residue_rq[arc] += 1

    G = []
    W = range(30)
    for _ in W:  # 构建分层图
        g = nx.DiGraph((x, y, {'weight': residue_rq[(x, y)]/bandwidth[(x, y)]})
                       for (x, y) in arcs)  # 权重：该边的起点或终点是一个未完成任务的个数/该边的空余带宽
        G.append(g)
    # g.out_degree(*)记录出*的度
    # g.in_degree(*)记录入*的度

    SD = collections.Counter(K)  # 任务
    print('平均值', np.mean(list(SD.values())))
    print('方差', np.var(list(SD.values())))

    sp_len = collections.defaultdict(list)  # 每种任务的最短路径长度以及波长
    al_simple_path = collections.defaultdict(int)
    for rq in SD:
        try:
            sp_len[rq] = [nx.shortest_path(
                G[0], rq[0], rq[1], weight='weight'), 0]
            paths = nx.all_simple_paths(G[0], rq[0], rq[1], len(sp_len[rq][0]))
            for _ in paths:
                al_simple_path[rq] += 1
        except:
            SD[rq] = 0

    C, C_len = solve(G, sp_len, al_simple_path, SD, W, residue_rq, bandwidth)

    # assignment = [[] for _ in W]  # 验证结果可靠性
    # K_copy = copy.deepcopy(K)
    # for w in range(len(W)):
    #     for rq in C[w]:
    #         assignment[w].append(K_copy.index(rq))
    #         K_copy[assignment[w][-1]] = 0
    # opt(W, K, arcs, V, assignment, C_len)

    num = 0
    for i in C:
        num += len(i)
    print(num)
    num = 0
    for i in C_len:
        for j in i:
            num += len(j)-1
    print(num)
