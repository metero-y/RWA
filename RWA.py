# -*- coding: utf-8 -*-
# Authors:   张䶮
"""
启发式算法：最小阻塞
"""
from cmath import inf
import collections
import heapq
import networkx as nx
import random
import numpy as np

from sqlalchemy import true


def solve(G, sp_len, SD, SD_idx, W):
    order = []
    for rq in sp_len:
        if SD[rq] == 0:
            continue
        # 起始点出deg，终止点入deg，最短路径值,越小越优先
        tmp = G[0].out_degree(rq[0])*len(W) * \
            G[0].in_degree(rq[1])*len(W) * sp_len[rq]
        if tmp != 0:  # 为0的时候自动忽略，因为这本身就不可行
            order.append((tmp, rq))
        else:
            SD[rq] = 0
    heapq.heapify(order)

    C = [[] for _ in W] # 总策略
    C_len = [[] for _ in W]  # 总策略对应路径长度
    while order:  # 边赋权还没实现
        require = heapq.heappop(order)
        min_val = float(inf)
        for w in range(len(W)):
            try:
                path = nx.shortest_path(G[w], require[1][0], require[1][1])
            except:
                continue
            now_val = len(path)-1
            if now_val < min_val:
                min_val = now_val
                wavelength = w
        if min_val < float(inf):
            C[wavelength].append(require[1])
            C_len[wavelength].append(min_val)
            for i in range(len(path)-1):  # 更新图
                G[wavelength].remove_edges_from([(path[i], path[i+1])])
            SD[require[1]] = SD[require[1]]-1

            order = []  # 更新order
            for rq in sp_len:
                if SD[rq] == 0:
                    continue
                # 起始点出deg，终止点入deg，最短路径值,越小越优先
                outdeg = 0
                indeg = 0
                for j in range(len(W)):
                    outdeg += G[j].out_degree(rq[0])
                    indeg += G[j].in_degree(rq[1])
                tmp = outdeg*indeg*sp_len[rq]
                if tmp != 0:  # 为0的时候自动忽略，因为这本身就不可行
                    order.append((tmp, rq))
                else:
                    SD[rq] = 0
            heapq.heapify(order)
        else:
            SD[require[1]] = 0
    return C,C_len


if __name__ == "__main__":
    L = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 7), (2, 5), (3, 4), (3, 8), (4, 5), (4, 6), (5, 10), (5, 12), (6, 7),
         (7, 9), (8, 11), (8, 13), (9, 10), (9, 11), (9, 13), (11, 12), (12, 13)]  # 边
    arcs = []  # 双向边
    for i in range(len(L)):
        arcs.append(L[i])
        arcs.append((L[i][1], L[i][0]))
    V = range(14)  # 结点

    G = []
    W = range(30)
    for _ in W:  # 构建分层图
        g = nx.DiGraph((x, y, {'weight': len(W)})
                       for (x, y) in arcs)  # 这里的权重还没改！！！！！
        G.append(g)
    # g.out_degree(*)记录出*的度
    # g.in_degree(*)记录入*的度

    random.seed(1015)  # 设置随机种子random.seed(2)
    K = []
    for _ in range(436):
        task = random.sample(range(14), 2)  # 结点不可重复
        K.append((task[0], task[1]))
    SD = collections.Counter(K)  # 任务
    print('平均值', np.mean(list(SD.values())))
    print('方差', np.var(list(SD.values())))
    
    SD_idx = []
    sp_len = collections.defaultdict(int)  # 每种任务的最短路径长度
    for rq in SD:
        try:
            sp_len[rq] = nx.shortest_path_length(G[0], rq[0], rq[1])
        except:
            sp_len[rq] = 0  # 没有最短路径时，设置为0，在求解时会自动忽略
            SD[rq] = 0
        SD_idx.append(rq)
        
    C,C_len = solve(G, sp_len, SD, SD_idx, W)

    num=0
    for i in C:
        num+=len(i)
    print(num)
