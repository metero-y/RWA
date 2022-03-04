# -*- coding: utf-8 -*-
# Authors:   张䶮
"""
启发式算法：最小阻塞
"""
import collections
from copy import deepcopy
from math import copysign
import networkx as nx
import random


def solve(dna, K, W, edges, nodes, val):
    LightPath = []
    ans = 0
    need = []
    for jjd in range(len(dna)):
        k = K[dna[jjd]]
        from_node = k[0]
        to_node = k[1]
        min_val = val
        for i in range(len(W)):
            pass  # 寻找最短加权路径
        if min_val < val:
            LightPath.append(lightpath)
            LightPath.append(wavelength)
            ans += min_val
        else:
            need.append(jjd)
    if need != []:
        return need, LightPath
    else:
        return ans, LightPath


if __name__ == "__main__":
    L = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 7), (2, 5), (3, 4), (3, 8), (4, 5), (4, 6), (5, 10), (5, 12), (6, 7),
         (7, 9), (8, 11), (8, 13), (9, 10), (9, 11), (9, 13), (11, 12), (12, 13)]  # 边
    arcs = []  # 双向边
    for i in range(len(L)):
        arcs.append(L[i])
        arcs.append((L[i][1], L[i][0]))
    V = range(14)  # 结点

    random.seed(1015)  # 设置随机种子random.seed(2)
    K = []
    for _ in range(436):
        task = random.sample(range(14), 2)  # 结点不可重复
        K.append((task[0], task[1]))
    SD = collections.Counter(K)  # 任务
    SD_idx = []
    for rq in SD:
        SD_idx.append(rq)

    G=[]
    W = range(30)
    for _ in W:  # 构建分层图
        g = nx.DiGraph()
        g.add_nodes_from(V)
        g.add_edges_from(arcs)
        print(g)
        G.append(g)
    deg=nx.degree(G[0])
    
