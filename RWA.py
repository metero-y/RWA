# -*- coding: utf-8 -*-
import numpy as np
import networkx as nx
import pickle as cp
import random
import ctypes
import os
import sys
import time
import glob
import re
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from gurobipy import *
## 表示无穷大
INF_val = 9999999


class Dijkstra_Path():
    def __init__(self, node_map):
        self.node_map = node_map
        self.node_length = len(node_map)
        self.used_node_list = []
        self.collected_node_dict = {}

    def __call__(self, from_node, to_node):
        self.from_node = from_node
        self.to_node = to_node
        self._init_dijkstra()
        return self._format_path()

    def _init_dijkstra(self):
        ## Add from_node to used_node_list
        self.used_node_list.append(self.from_node)
        for index1 in range(self.node_length):
            self.collected_node_dict[index1] = [INF_val, -1]

        self.collected_node_dict[self.from_node] = [0, -1]  # from_node don't have pre_node
        for index1, weight_val in enumerate(self.node_map[self.from_node]):
            if weight_val:
                self.collected_node_dict[index1] = [weight_val, self.from_node]  # [weight_val, pre_node]
        self._foreach_dijkstra()

    def _foreach_dijkstra(self):
        while (len(self.used_node_list) < self.node_length - 1):
            min_key = -1
            min_val = INF_val
            for key, val in self.collected_node_dict.items():  # 遍历已有权值节点
                if val[0] < min_val and key not in self.used_node_list:
                    min_key = key
                    min_val = val[0]

                    ## 把最小的值加入到used_node_list
            if min_key != -1:
                self.used_node_list.append(min_key)
            else:
                return [(self.from_node, INF_val)]

            for index1, weight_val in enumerate(self.node_map[min_key]):
                ## 对刚加入到used_node_list中的节点的相邻点进行遍历比较
                if weight_val > 0 and self.collected_node_dict[index1][0] > weight_val + min_val:
                    self.collected_node_dict[index1][0] = weight_val + min_val  # update weight_val
                    self.collected_node_dict[index1][1] = min_key


    def _format_path(self):
        node_list = []
        temp_node = self.to_node
        node_list.append((temp_node, self.collected_node_dict[temp_node][0]))
        while self.collected_node_dict[temp_node][1] != -1:
            temp_node = self.collected_node_dict[temp_node][1]
            node_list.append((temp_node, self.collected_node_dict[temp_node][0]))
        node_list.reverse()
        return node_list


def set_node_map(node_map, node, node_list):
    for x, y in node_list:
        node_map[node.index(x)][node.index(y)] = node_map[node.index(y)][node.index(x)] = 1

def gen_graph(max_n,min_n,density,newVal):
    test_max_n = max_n
    test_min_n = min_n

    cur_n = np.random.randint(test_max_n - test_min_n  + 1) + test_min_n

    m = density
    g = nx.barabasi_albert_graph(n = cur_n, m = m, seed=newVal)

    return g

def solve(dna,K,W,edges,nodes,val):
    LightPath=[]
    ans=0
    need=[]
    for jjd in range(len(dna)):
        k = K[dna[jjd]]
        from_node = k[0]
        to_node = k[1]
        min_val=val
        for i in range(len(W)):
            w_edges = edges[:]
            for j in range(1, len(LightPath) + 1, 2):
                if LightPath[j] == i:
                    temp = LightPath[j - 1]
                    for apk in range(len(temp) - 1):
                        if temp[apk][0] < temp[apk + 1][0]:
                            edge = (temp[apk][0], temp[apk + 1][0])
                        else:
                            edge = (temp[apk + 1][0], temp[apk][0])
                        w_edges.remove(edge)

            node_map = [[0 for val in range(len(nodes))] for val in range(len(nodes))]
            set_node_map(node_map, nodes, w_edges)
            dijkstrapath = Dijkstra_Path(node_map)
            path = dijkstrapath(from_node, to_node)
            now_val = path[-1][1]
            if now_val < min_val:
                min_val = now_val
                wavelength = i
                lightpath = path
        if min_val < val:
            LightPath.append(lightpath)
            LightPath.append(wavelength)
            ans += min_val
        else:
            need.append(jjd)
    if need != []:
        return need,LightPath
    else:
        return ans,LightPath

def decode(path,W):
    graph=[[] for i in range(len(W))]
    for w in range(len(path)):
        LightPath=path[w]
        for i in range(len(W)):
            for j in range(1, len(LightPath) + 1, 2):
                if LightPath[j] == i:
                    temp = LightPath[j - 1]
                    for apk in range(len(temp) - 1):
                        if temp[apk][0] < temp[apk + 1][0]:
                            edge = (temp[apk][0], temp[apk + 1][0])
                        else:
                            edge = (temp[apk + 1][0], temp[apk][0])
                        if edge in graph[i]:
                            pass
                        else:
                            graph[i].append(edge)
    return graph

seed=20
np.random.seed(seed)
random.seed((seed))
val = 9999999
if __name__ == "__main__":
    percentage=[]
    answer=[]
    answer_simple=[]
    objVal=[]
    renwushu=200
    pop_size=10
    pc=0.8#交叉概率
    pm=0.01#变异概率
    train=0
    while train <10:
        train+=1

        np.random.seed(train)
        print('train',train)
        edges = [(0, 1), (0, 2), (0, 7), (1, 2), (1, 3), (2, 5), (3, 4), (3, 10), (4, 5), (4, 6), (5, 9), (5, 13),
                 (6, 7), (7, 8), (8, 9), (8, 11), (8, 12), (10, 11), (10, 12), (11, 13), (12, 13)]
        nodes = range(14)
        K=[]
        for vad in range(renwushu):
            res = random.sample(range(len(nodes)), 2)
            K.append((res[0], res[1]))
        if train==6 or train ==7:
            answer.append('无解')
            continue
        W = range(30)
        pop=[]
        meaningless=[]
        fitness = []
        times=0
        graph=[]
        time_1=time.time()
        while len(pop)!=pop_size and times<30:
            if (times+1)%31==0:
                print('30遍了',len(pop))
            dna=list(np.random.choice(range(renwushu), renwushu, replace=False))
            if dna in meaningless:
                print('重复啦')
                continue
            times+=1
            dna_is_false=0
            www=[]
            # while dna_is_false < 50:
            #     sol, LightPath = solve(dna, K, W, edges, nodes, val)
            #     if type(sol) != list:
            #         fitness.append(sol)
            #         pop.append(dna)
            #         graph.append(LightPath)
            #         if dna_is_false != 0:
            #             print('很好奇有几个小坏蛋', dna_is_false)
            #         break
            #     else:
            #         meaningless.append(dna)
            #         frag = dna[:sol[0]]
            #         frag2 = dna[sol[0]:]
            #         dna = frag2 + frag
            #         if dna in www:
            #             print('循环了', dna_is_false)
            #         else:
            #             www.append(dna)
            #         dna_is_false += 1
            #         if dna_is_false % 300 == 0:
            #             print('天哪')
            # while dna_is_false < 1000:
            #     sol, LightPath = solve(dna, K, W, edges, nodes, val)
            #     if type(sol) != list:
            #         fitness.append(sol)
            #         pop.append(dna)
            #         graph.append(LightPath)
            #         if dna_is_false != 0:
            #             print('很好奇有几个小坏蛋', dna_is_false)
            #         break
            #     else:
            #         meaningless.append(dna)
            #         frag = dna[:sol[0]]
            #         for i in range(1, len(sol)):
            #             frag += dna[sol[i - 1] + 1:sol[i]]
            #         if sol[-1] != renwushu - 1:
            #             frag += dna[sol[-1] + 1:]
            #         frag2 = [dna[sol[0]]]
            #         for i in range(1, len(sol)):
            #             frag2 += [dna[sol[i]]]
            #         dna = frag2 + frag
            #         if dna in www:
            #             print('循环了', dna_is_false)
            #         else:
            #             www.append(dna)
            #         dna_is_false += 1
            #         if dna_is_false % 300 == 0:
            #             print('天哪')
        print(len(pop), len(meaningless))
        if len(pop)==0:
            train-=1
            print(len(pop), len(meaningless))
            print('出错了，再来一遍')
            continue
        while len(pop)!=pop_size:
            dna=list(np.random.choice(range(renwushu), renwushu, replace=False))
            if dna in meaningless:
                print('重复啦')
                continue
            sol,LightPath=solve(dna,K,W,edges,nodes,val)
            if type(sol)!=list:
                fitness.append(sol)
                pop.append(dna)
                graph.append(LightPath)
            else:
                meaningless.append(dna)
                fitness.append(val)
                pop.append(dna)
                graph.append('')
        print(len(graph))
        print(fitness, 'fitness')
        iterate=0
        improvement=0
        regist=[]
        regist_graph=[]
        temp_graph=graph[:]
        regist_graph.append(temp_graph)
        regist.append(min(fitness))
        while iterate<=100:
            now_sol=min(fitness)
            new_pop = pop[:]
            #交叉遗传
            for i in range(pop_size):
                if (np.random.rand() < pc):
                    i_ = np.random.randint(0, pop_size)
                    # 交叉位置
                    y = np.random.randint(0, renwushu)
                    # 记录交叉项
                    fragment1 = pop[i][y:]
                    fragment2 = pop[i_][y:]
                    a1_2 = []
                    a2_2 = []
                    for j in pop[i][:y]:
                        while j in fragment2:
                            j = fragment1[fragment2.index(j)]
                        a1_2.append(j)
                    for j in pop[i_][:y]:
                        while j in fragment1:
                            j = fragment2[fragment1.index(j)]
                        a2_2.append(j)
                    new_pop[i] = a1_2 + fragment2
                    new_pop[i_] = a2_2 + fragment1
            #变异
            for i in range(pop_size):
                if (np.random.rand() < pm):
                    mp_1,mp_2 = np.random.choice(range(renwushu), 2, replace=False)
                    new_pop[i][mp_1],new_pop[i][mp_2]=new_pop[i][mp_2],new_pop[i][mp_1]
            # 适者生存
            new_fitness=[]
            new_graph=[]
            for i in range(pop_size):
                sol,LightPath = solve(new_pop[i], K, W, edges, nodes,val)
                if type(sol)!=list:
                    new_fitness.append(sol)
                    new_graph.append(LightPath)
                else:
                    new_fitness.append(val)
                    new_graph.append('')
            newpop = pop + new_pop
            newfitness=fitness+new_fitness

            newgraph=graph+new_graph


            minval=min(newfitness)
            index=0
            for i in range(len(newfitness)):
                if newfitness[i]==minval:
                    pop[index]=newpop[i]
                    fitness[index]=newfitness[i]
                    graph[index]=newgraph[i]
                    index+=1
                if index==pop_size//3:
                    break


            total = 0
            for i in range(pop_size*2):
                total += 1 / newfitness[i]
            p = []
            for i in range(pop_size*2):
                p.append(1 / newfitness[i] / total)  # 将所有个体的适应度归一化列表
            normal_fitness = []
            for i in range(pop_size*2):  # （计算累计概率）
                total = 0
                j = 0
                while (j <= i):
                    total += p[j]
                    j += 1
                normal_fitness.append(total)
            MS = []
            for i in range(index,pop_size):
                MS.append(np.random.rand())
            for i in range(pop_size-index):
                if (MS[i] <= normal_fitness[0]):
                    pop[i+index] = newpop[0]
                    fitness[i+index]=newfitness[0]
                    graph[i+index]=newgraph[0]
                else:
                    for j in range(1, pop_size*2):
                        if MS[i] <= normal_fitness[j] and MS[i] > normal_fitness[j - 1]:
                            pop[i+index] = newpop[j]
                            fitness[i+index] = newfitness[j]
                            graph[i+index] = newgraph[j]        # 对原种群，将存活个体进行保存到新种群

            #更新fitness
            aft_sol=min(fitness)
            print(fitness,'fitness')
            difference=now_sol-aft_sol
            regist.append(aft_sol)
            temp_graph=graph[:]
            regist_graph.append(temp_graph)
            iterate += 1
            if difference==0 and aft_sol!=val:
                improvement+=1
            else:
                improvement=0
            if improvement==3:
                break
        answer.append(min(regist))
        time_2=time.time()
        print(time_2-time_1)


        # time_3=time.time()
        # best_graph = [[] for i in range(len(W))]
        # try:
        #     arcs = []
        #     node = []
        #     for i in range(len(edges)):
        #         arcs.append(edges[i])
        #         arcs.append((edges[i][1], edges[i][0]))
        #         node.append(edges[i][0])
        #         node.append(edges[i][1])
        #
        #     from_node = []
        #     end_node = []
        #     for demand in range(len(K)):
        #         from_node.append(K[demand][0])
        #         end_node.append((K[demand][1]))
        #
        #
        #     m = gurobipy.Model('demo')
        #     xwke = m.addVars(W, range(len(K)), arcs, vtype=GRB.BINARY, name='xkew')
        #     ywk = m.addVars(W, range(len(K)), vtype=GRB.BINARY, name='ywk')
        #
        #     m.addConstrs((ywk.sum('*', k) == 1 for k in range(len(K))), name='每个任务只分配一个波长')
        #     m.addConstrs((xwke.sum(w, '*', i, j) <= 1 for w in W for i, j in arcs), name='同一波长链路只能有一个信息')
        #     m.addConstrs((xwke[w, k, i, j] <= ywk[w, k] for w in W for k in range(len(K)) for i, j in arcs),
        #                  name='不在任务波长下的任务链路为0')
        #
        #     m.addConstrs((xwke[w, k, i, j] + xwke[w, k, j, i] <= 1 for w in W for k in range(len(K)) for i, j in arcs),
        #                  "正反向只有一个通")
        #     m.addConstrs(
        #         (xwke.sum(w, k, '*', j) == xwke.sum(w, k, j, '*') for w in W for k in range(len(K)) for j in node if
        #          j != K[k][0] and j != K[k][1]), "流守恒")
        #     m.addConstrs((xwke.sum(w, k, K[k][0], '*') - xwke.sum(w, k, '*', K[k][0]) == ywk[w, k] for w in W for k in
        #                   range(len(K))), "发源地")
        #     m.addConstrs((xwke.sum(w, k, K[k][1], '*') - xwke.sum(w, k, '*', K[k][1]) == -ywk[w, k] for w in W for k in
        #                   range(len(K))), "接受地")
        #
        #     m.setObjective(xwke.sum(), GRB.MINIMIZE)
        #
        #     m.optimize()
        #     for v in m.getVars():
        #         if v.x == 1 and v.varName[0] == 'x':
        #             num = -1
        #             a = ''
        #             b=''
        #             for i in range(len(v.varName)):
        #                 if v.varName[i] == ']':
        #                     break
        #                 if num == 2 or num == 3:
        #                     a += v.varName[i]
        #                 if v.varName[i] == ',':
        #                     num += 1
        #                 if num==0:
        #                     b+=v.varName[i]
        #                 if v.varName[i]=='[':
        #                     num+=1
        #             for i in range(len(a)):
        #                 if a[i]==',':
        #                     break
        #             num1=int(a[:i])
        #             num2=int(a[i+1:])
        #             if num1<num2:
        #                 best_graph[int(b)].append((num1,num2))
        #             else:
        #                 best_graph[int(b)].append((num2, num1))
        #     objVal.append(int(m.objVal))
        # except gurobipy.GurobiError as e:
        #     print('Errorcode ' + str(e.errno) + ": " + str(e))
        #     objVal.append(None)
        # except AttributeError:
        #     print('Encountered an attribute error')
        #     objVal.append(None)
        # time_4=time.time()
        # print(time_4-time_3)
        # ssd=[]
        # a=regist_graph[-1]
        # for i in range(len(a)):
        #     temp=decode([a[i]],W)
        #     if temp not in ssd:
        #         ssd.append(temp)
        # print(len(ssd))
        #
        # graph = decode(regist_graph[regist.index(min(regist))], W)#-1
        # time_5 = time.time()
        # try:
        #     from_node = []
        #     end_node = []
        #     for demand in range(len(K)):
        #         from_node.append(K[demand][0])
        #         end_node.append((K[demand][1]))
        #     all_arcs = []
        #     all_node = []
        #     m = gurobipy.Model('demo')
        #     ywk = []
        #     xwke = []
        #     for jj in range(len(W)):
        #         arcs = []
        #         node = []
        #         for i in range(len(graph[jj])):
        #             arcs.append(graph[jj][i])
        #             arcs.append((graph[jj][i][1], graph[jj][i][0]))
        #             if graph[jj][i][0] in node:
        #                 pass
        #             else:
        #                 node.append(graph[jj][i][0])
        #             if graph[jj][i][1] in node:
        #                 pass
        #             else:
        #                 node.append(graph[jj][i][1])
        #         all_arcs.append(arcs)
        #         all_node.append(node)
        #         xwke.append(m.addVars([jj], range(len(K)), arcs, vtype=GRB.BINARY, name='xkew'))
        #         ywk.append(m.addVars([jj], range(len(K)), vtype=GRB.BINARY, name='ywk'))
        #     for w in W:
        #         m.addConstrs((xwke[w].sum(w, '*', i, j) <= 1 for i, j in all_arcs[w]), name='同一波长链路只能有一个信息')
        #         m.addConstrs((xwke[w][w, k, i, j] <= ywk[w][w, k] for k in range(len(K)) for i, j in all_arcs[w]),
        #                      name='不在任务波长下的任务链路为0')
        #
        #         m.addConstrs(
        #             (xwke[w][w, k, i, j] + xwke[w][w, k, j, i] <= 1 for k in range(len(K)) for i, j in all_arcs[w]),
        #             "正反向只有一个通")
        #         m.addConstrs(
        #             (xwke[w].sum(w, k, '*', j) == xwke[w].sum(w, k, j, '*') for k in range(len(K)) for j in all_node[w]
        #              if
        #              (j != K[k][0] and j != K[k][1] and K[k][0] in all_node[w] and K[k][1] in all_node[w]) or (
        #                          K[k][0] not in all_node[w] or K[k][1] not in all_node[w])), "流守恒")
        #         m.addConstrs((xwke[w].sum(w, k, K[k][0], '*') - xwke[w].sum(w, k, '*', K[k][0]) == ywk[w][w, k] for k in
        #                       range(len(K)) if K[k][0] in all_node[w] and K[k][1] in all_node[w]), "发源地")
        #         m.addConstrs(
        #             (xwke[w].sum(w, k, K[k][1], '*') - xwke[w].sum(w, k, '*', K[k][1]) == -ywk[w][w, k] for k in
        #              range(len(K)) if K[k][1] in all_node[w] and K[k][0] in all_node[w]), "接受地")
        #
        #     for k in range(len(K)):
        #         obj = LinExpr()
        #         for i in range(len(W)):
        #             if K[k][1] in all_node[i] and K[k][0] in all_node[i]:
        #                 obj += ywk[i][i, k]
        #
        #         m.addConstr(obj == 1, name='每个任务只分配一个波长{}'.format(k))
        #     obj = LinExpr()
        #     for i in range(len(W)):
        #         obj += xwke[i].sum()
        #     m.setObjective(obj, GRB.MINIMIZE)
        #
        #     m.optimize()
        #     print(m.objVal)
        # except gurobipy.GurobiError as e:
        #     print('Errorcode ' + str(e.errno) + ": " + str(e))
        # except AttributeError:
        #     print('Encountered an attribute error')
        # time_6 = time.time()
        # print(time_6 - time_5)
        #
        #
        #
        #


        # count = 0
        # all_count = 0
        # co=0
        # for i in range(len(W)):
        #     for j in range(len(best_graph[i])):
        #         if best_graph[i][j] in graph[i]:
        #             count += 1
        #     if best_graph[i] != []:
        #         all_count += len(best_graph[i])
        #     if graph[i]!=[]:
        #         co+=len(graph[i])
        # print(count,all_count,co)
    #     percentage.append((count,all_count,co))
    #     sol,LightPath=solve(range(renwushu),K,W,edges,nodes,val)
    #     if type(sol)!=list:
    #         answer_simple.append(sol)
    #     else:
    #        answer_simple.append('无解')
    # print(percentage)
    print(answer)
    # print(objVal)
    # print(answer_simple)
