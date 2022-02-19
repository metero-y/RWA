# -*- coding: utf-8 -*-
import collections
import copy

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
import math
import matplotlib

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

## 表示无穷大
INF_val = float("inf")


# 生成任务
def generator_task(task_num, node_nums):
    """
    :param task_num  : 需要生成的任务数
    :param node_nums : 图中结点数量
    :return: tasks   : 生成的任务连接请求
    """
    tasks = []
    task_index = range(node_nums)
    for _ in range(task_num):
        task = random.sample(task_index, 2)
        tasks.append((task[0], task[1]))
    return tasks


# 判断任务请求在图G中是否存在路径, 若存在返回对应路径
def has_path(G, task):
    """
    :param G:     图
    :param task: 任务请求
    :return:  path:路径
    """
    try:
        path = nx.shortest_path(G, source=task[0], target=task[1])

    except nx.NetworkXNoPath:
        return False, None
    return True, path


# 得到连接请求的最短路径, 并按照降序进行排列
def Paths_order(tasks, G):
    """
    :param   tasks: 任务请求
    :param       G: 图的结构信息
    :return: order: 任务连接请求按照最短路径的降序排列
             Paths: 每个任务连接请求的对应最短路径， 找不到路径的序对的路径为空
    """
    num = []  # 存储每个连接请求的最短路径数值
    Paths = []  # 存储每个连接请求的最短路径
    for task in tasks:
        exist_path, path = has_path(G, task)
        if exist_path:
            num.append(len(path) - 1)
            Paths.append(path)
        else:  # 找不到最短路径
            num.append(0)
            Paths.append([0])
    order = sorted(range(len(num)), key=lambda k: num[k], reverse=True)  # 每个连接请求的最短路径按长度降序排列, 存储对应的任务请求序号
    Path = copy.deepcopy(Paths)
    for i in range(len(order)):
        Path[i] = Paths[order[i]]

    return order, Path


# 判断加权图中是否存在路径
def has_weight_path(G, task):
    """
       :param G:     图
       :param task: 任务请求
       :return:  exist: 判断是否存在路径
    """

    path = nx.shortest_path(G, source=task[0], target=task[1], weight='weight')
    path_1 = nx.shortest_path_length(G, source=task[0], target=task[1], weight='weight')
    if path_1 < INF_val:
        return True, path
    else:
        return False, None


def LFFP(tasks, edges, nodes):
    """
    :param          tasks: 任务请求序列
    :param          edges: 边序列
    :param          nodes: 结点序列
    :return: Occupathions: 每个波长层中边的占用情况
                LightPath: 存储每条路径所及其所对应的波长
                    order: 任务请求的降序排列顺序
                  require: 存储每个波长层中的连接请求
    """
    can_route = True  # 标识路径在该波长层中是否存在
    w = 0  # 标识波长层序号
    LightPath = []  # 存储每个连接请求的路径及其对应的波长层

    Graphs = []
    Graph = nx.Graph()
    for i in range(len(nodes)):
        Graph.add_node(nodes[i])
    for x, y in edges:  # edges:
        Graph.add_edges_from([(x, y)])
    Graphs.append(Graph.copy())

    order, Path = Paths_order(tasks, Graph)  # 得到连接请求的最短路径, 并按照降序进行排列

    Occupathions = []  # 存储每个波长层中边的占用情况
    occupathion = {}
    for e in edges[:]:
        occupathion[e] = False
    Occupathions.append(occupathion.copy())
    require = [[]]

    for i in range(len(Path)):
        path = Path[i]
        Edge = []  # 存储路径中所设计的边
        for apk in range(len(path) - 1):
            if path[apk] < path[apk + 1]:
                edge = (path[apk], path[apk + 1])
            else:
                edge = (path[apk + 1], path[apk])
            Edge.append(edge)
        for w in range(len(Occupathions)):
            can_route = True
            # 判断该连接请求的最短路径是否存在当前的波长层中
            for e in Edge:
                if Occupathions[w][e]:
                    can_route = False
                    break
            if can_route == True:  # 分配波长
                LightPath.append((path, w))  # 记录每条路径对应的波长层
                require[w].append(order[i])  # 记录每个波长层所分配的任务请求
                for e in Edge:
                    Occupathions[w][e] = True
                break
        if can_route == False:  # 如果当前所有波长层都无法路由，就添加新的波长层
            Occupathions.append(occupathion.copy())
            LightPath.append((path, w + 1))
            require.append([order[i]])
            for e in Edge:
                Occupathions[-1][e] = True

    return Occupathions, LightPath, order, require


def HPLD(tasks, edges, nodes):
    """
        :param          tasks: 任务请求序列
        :param          edges: 边序列
        :param          nodes: 结点序列
        :return: Occupathions: 每个波长层中边的占用情况
                    LightPath: 存储每条路径所及其所对应的波长
                      require: 存储每个波长层中的连接请求
    """
    Occupathions, LightPath, order, require = LFFP(tasks, edges, nodes)  # 利用LFFP分配路径
    task_num = len(tasks)  # 任务数量
    edge_num = len(edges)  # 边数量
    reduce_wave = []  # 存储被简化的波长层
    sigma = 1
    edges_cost = {}
    edge_load = {}
    for e in edges:
        edge_load[e] = 0
    for i in range(task_num):
        path = LightPath[i][0]
        for apk in range(len(path) - 1):
            if path[apk] < path[apk + 1]:
                edge = (path[apk], path[apk + 1])
            else:
                edge = (path[apk + 1], path[apk])
            edge_load[edge] += 1
    leiji = [] #防止出现循环
    while LightPath not in leiji:
        leiji.append(LightPath[:])
        # 负载偏差量的确定
        mean_load = 0
        max_load = 0
        num_load = 0
        for e in edge_load:
            if edge_load[e] > 0:
                mean_load += edge_load[e]
                num_load += 1
                if max_load < edge_load[e]:
                    max_load = edge_load[e]
        mean_load /= edge_num  # 分母应该为总边数
        deta = int(sigma * (max_load - mean_load))
        if deta == 0:
            return Occupathions, LightPath, require, reduce_wave
        for e in edges:
            if max_load - 1 - edge_load[e] > 0:
                edges_cost[e] = 1 / (max_load - 1 - edge_load[e])
            else:
                edges_cost[e] = INF_val
        Max_Load_Path_index = []  # 记录最大负载路径的索引
        for i in range(task_num):
            path = LightPath[i][0]
            for apk in range(len(path) - 1):
                if path[apk] < path[apk + 1]:
                    e = (path[apk], path[apk + 1])
                else:
                    e = (path[apk + 1], path[apk])
                if edge_load[e] == max_load:
                    Max_Load_Path_index.append(i)
                    break
        Graph = nx.Graph()
        for i in range(len(nodes)):
            Graph.add_node(nodes[i])
        for e in edges:  # edges:
            x, y = e[0], e[1]
            Graph.add_weighted_edges_from([(x, y, edges_cost[e])])

        Exist_Path = []
        for path_idx in Max_Load_Path_index:
            exist_path, haspath = has_weight_path(Graph, tasks[order[path_idx]])
            if exist_path:
                Exist_Path.append((path_idx, haspath))
        if not Exist_Path:
            return Occupathions, LightPath, require, reduce_wave
        else:
            not_change = True
            while Exist_Path and not_change:
                random_index_index = np.random.choice(range(len(Exist_Path)))
                random_index = Exist_Path[random_index_index]
                Exist_Path.remove(random_index)
                # print(Exist_Path)
                p = LightPath[random_index[0]][0]  # 得到原路径
                ori = LightPath[random_index[0]][1]  # 得到原路径所在波长层
                ori_edges = []  # 记录原路径的边
                # 移除原路径
                for apk in range(len(p) - 1):
                    if p[apk] < p[apk + 1]:
                        e = (p[apk], p[apk + 1])
                    else:
                        e = (p[apk + 1], p[apk])
                    ori_edges.append(e)
                    edge_load[e] -= 1
                    Occupathions[ori][e] = False
                    # print(require,order,random_index)
                require[ori].remove(order[random_index[0]])
                if len(require[ori]) == 0:
                    reduce_wave.append(ori)
                # 计算新路径
                path = random_index[1]
                Edge = []  # 存储新的路径
                for apk in range(len(path) - 1):
                    if path[apk] < path[apk + 1]:
                        edge = (path[apk], path[apk + 1])
                    else:
                        edge = (path[apk + 1], path[apk])
                    Edge.append(edge)
                can_route = True
                for w in range(len(Occupathions)):
                    if w in reduce_wave:
                        continue
                    can_route = True
                    for e in Edge:
                        if Occupathions[w][e]:
                            can_route = False
                            break
                    if can_route:
                        require[w].append(order[random_index[0]])
                        LightPath[random_index[0]] = (path, w)
                        for e in Edge:
                            Occupathions[w][e] = True
                            edge_load[e] += 1
                        not_change = False
                        break
                if can_route == False:
                    for e in ori_edges:
                        edge_load[e] += 1
                        Occupathions[ori][e] = True
                    require[ori].append(order[random_index[0]])
                    if len(require[ori]) == 1:
                        reduce_wave.remove(ori)
            if not_change:
                return Occupathions, LightPath, require, reduce_wave
    return Occupathions, LightPath, require, reduce_wave


if __name__ == "__main__":
    random.seed(2)  # 设置随机种子
    wave_link_length = []  # 波长链路数
    times = []  # 时间
    wave_num = []  # 波长数
    node_nums = 14  # 结点数量
    nodes = np.arange(node_nums)  # 结点序列
    edges = [(0, 1), (0, 2), (0, 7), (1, 2), (1, 3), (2, 5), (3, 4), (3, 10), (4, 5), (4, 6), (5, 9), (5, 13),
             (6, 7), (7, 8), (8, 9), (8, 11), (8, 12), (10, 11), (10, 13), (11, 13), (12, 13)]
    # 边连接关系
    Max_num_tasks = 200

    for task_num in range(1, Max_num_tasks + 1):

        tasks = generator_task(task_num, node_nums)  # 生成任务

        start_time = time.time()  # 求解起始时间

        Occupathions, LightPath, require, reduce_wave = HPLD(tasks, edges, nodes)

        end_time = time.time()  # 求解终止时间

        times.append(end_time - start_time)
        ns = 0
        for i in range(len(LightPath)):
            ns = len(LightPath[i][0]) + ns - 1
        wave_link_length.append(ns)

        wave_num.append(len(Occupathions) - len(reduce_wave))
    print(wave_link_length)
    print(sum(wave_link_length))
    print(sum(times))
    print(wave_num)
    print(sum(wave_num))
