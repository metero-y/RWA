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
    for _ in range(task_num):
        task = random.sample(range(node_nums), 2)
        tasks.append((task[0], task[1]))
    tasks_array = np.array(tasks)
    return tasks_array

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


#连接请求按照最短路径的递减顺序排列，用FF-RWA方法进行路由。
def FFD_RWA(order, tasks, edges, nodes):
    """
    :param   order: 任务连接请求序号
    :param   tasks: 任务连接请求
    :param   edges: 边序列
    :param   nodes: 结点序列
    :return:Graphs: 每一个波长层的剩余边连接关系
         LightPath: 存储每条路径所及其所对应的波长
              load: 每个波长层的负载数
           require: 每个波长层中任务请求
    """
    LightPath = []               # 存储每个连接请求的路径及其对应的波长层
    Graph = nx.DiGraph()
    for i in range(len(nodes)):
        Graph.add_node(nodes[i])
    for x, y in edges:  # edges::
        Graph.add_edges_from([(x, y)])
        Graph.add_edges_from([(y, x)])
    Graphs = []
    Graphs.append(Graph.copy())
    require = [[]]
    load = [0]
    for id in order:
        task = tasks[id]
        can_rounte = False
        w = 0
        for w in range(len(Graphs)):
            G = Graphs[w]
            exist_path, path = has_path(G, task)
            if exist_path:
              LightPath.append((path, w))
              load[w] = load[w] + len(path) - 1
              require[w].append(id)
              can_rounte = True
              for apk in range(len(path) - 1):
                  e = (path[apk + 1], path[apk])
                  Graphs[w].remove_edge(e[0], e[1])
                  Graphs[w].remove_edge(e[1], e[0])
              break
        if can_rounte == False:
                Graphs.append(Graph.copy())
                G = Graphs[-1]
                path = nx.shortest_path(G, source=task[0], target=task[1])
                LightPath.append((path, w+1))
                require.append([])
                load.append(0)
                load[-1] = load[-1] + len(path) - 1
                require[-1].append(id)
                for apk in range(len(path) - 1):
                    e = (path[apk], path[apk + 1])
                    Graphs[-1].remove_edge(e[0], e[1])
                    Graphs[-1].remove_edge(e[1], e[0])
    return Graphs, LightPath, load, require


# 局部搜索算法: 对于每一个连接请求, 都试图将其重新路由到使用链路更多的波长层中
def Shift_paths(order, tasks, Graphs, require, load, LightPath):
    """
    :param     order: 任务连接请求序号
    :param     tasks: 任务连接请求
    :param    Graphs: 每一个波长层的剩余边连接关系
    :param   require: 每个波长层中任务请求
    :param      load: 每个波长层的负载数
    :param LightPath: 存储每条路径所及其所对应的波长
    :return:  Graphs: 每一个波长层的剩余边连接关系
           LightPath: 存储每条路径所及其所对应的波长
                load: 每个波长层的负载数
             require: 每个波长层中任务请求
         reduce_wave: 可简化的波长
    """
    nums_task = len(tasks)           # 任务连接请求数量
    reduce_wave = []
    for id in range(nums_task):
        can_route_W = []             # 存储连接请求的可路由波长
        can_route   = []             # 存储连接请求的路径
        task        = tasks[order[id]]
        now_w       = LightPath[id][1]
        for w in range(len(Graphs)):
            if w in reduce_wave:
                continue
            G = Graphs[w]
            exist_path, path = has_path(G, task)
            if exist_path:
                can_route_W.append(w)
                can_route.append(path)

        if can_route_W != []:
            new_w = can_route_W[0]   # 存储负载最大的波长层
            for i in can_route_W:
                if load[i] > load[new_w]:
                    new_w = i

            if load[new_w] > load[now_w]:
                orign_path  = LightPath[id][0]
                load[now_w] = load[now_w] - len(orign_path) + 1
                for apk in range(len(orign_path) - 1):
                    e = (orign_path[apk + 1], orign_path[apk])
                    Graphs[now_w].add_edges_from([(e[0], e[1])])
                    Graphs[now_w].add_edges_from([(e[1], e[0])])
                if load[now_w] == 0:
                    reduce_wave.append(now_w)
                require[now_w].remove(order[id])
                path = can_route[can_route_W.index(new_w)]
                LightPath[id] = (path, new_w)

                for apk in range(len(path) - 1):
                    e = (path[apk + 1], path[apk])
                    Graphs[new_w].remove_edge(e[0], e[1])
                    Graphs[new_w].remove_edge(e[1], e[0])
                load[new_w] = load[new_w] + len(path) - 1
                require[new_w].append(order[id])


    return Graphs, LightPath, load, require, reduce_wave

# 随机选择两个波长层，将使用链路较少的波长层中的连接请求的路径移到使用链路较多的波长层中，并将发生冲突的链路所对应的连接请求重新路由。
def Mutate(order, tasks, Graphs, require, load, reduce_wave, LightPath, nodes, edges):
    """
    :param        order: 任务连接请求序号
    :param        tasks: 任务连接请求
    :param       Graphs: 每一个波长层的剩余边连接关系
    :param      require: 每个波长层中任务请求
    :param         load: 每个波长层的负载数
    :param  reduce_wave: 可简化的波长
    :param    LightPath: 存储每条路径所及其所对应的波长
    :param        nodes: 结点序列
    :param        edges: 边序列
    :return:     Graphs: 每一个波长层的剩余边连接关系
              LightPath: 存储每条路径所及其所对应的波长
                   load: 每个波长层的负载数
                require: 每个波长层中任务请求
            reduce_wave: 可简化的波长
    """
    waves = list(range(len(Graphs)))
    for wave in reduce_wave:
        waves.remove(wave)
    if len(waves) < 2:
        return Graphs, LightPath, load, require, reduce_wave
    w_1, w_2 = np.random.choice(waves, 2, replace=False)  # 随机选择两个波长层
    if load[w_1] < load[w_2]:
        w_1, w_2 = w_2, w_1
    r_2 = require[w_2]  # 链路较少的波长层中的连接请求
    out_path_index = np.random.choice(r_2)  # 在链路较少的波长层中随机选择一个连接请求
    out_index = order.index(out_path_index)
    out_path = LightPath[out_index][0]     # 所选择连接请求的路径
    out_edge = {}  # 存储所选择连接请求的路径中的边
    for apk in range(len(out_path) - 1):
        if out_path[apk] < out_path[apk + 1]:
            e = (out_path[apk], out_path[apk + 1])
        else:
            e = (out_path[apk + 1], out_path[apk])
        out_edge[e] = True
        Graphs[w_2].add_edges_from([(e[0], e[1])])
        Graphs[w_2].add_edges_from([(e[1], e[0])])

    load[w_2] = load[w_2] - len(out_path) + 1
    require[w_2].remove(out_path_index)

    r_1 = require[w_1].copy() # 链路较多的波长层中的连接请求

    R   = []
    for i in range(len(r_1)):
        r = r_1[i]
        conflict = False
        index = order.index(r)
        path  = LightPath[index][0]  # 链路较多的波长层中的连接请求对应路径
        for apk in range(len(path) - 1):
            if path[apk] < path[apk + 1]:
                e = (path[apk], path[apk + 1])
            else:
                e = (path[apk + 1], path[apk])
            if e in out_edge:
                conflict = True
                R.append(r)
                require[w_1].remove(r)
                load[w_1] = load[w_1] - len(path) + 1
                break
        if conflict:
            for apk in range(len(path) - 1):
                e = (path[apk], path[apk + 1])
                Graphs[w_1].add_edges_from([(e[0], e[1])])
                Graphs[w_1].add_edges_from([(e[1], e[0])])

    LightPath[out_index] = (out_path, w_1)
    load[w_1] = load[w_1] + len(out_path) - 1
    require[w_1].append(out_path_index)
    for apk in range(len(out_path) - 1):
        e = (out_path[apk + 1], out_path[apk])
        Graphs[w_1].remove_edge(e[0], e[1])
        Graphs[w_1].remove_edge(e[1], e[0])

    for r in R:
        index = order.index(r)
        task  = tasks[index]   # 该连接请求所对应的任务
        can_route = False
        for w in waves:
            if w == w_1:
                continue
            G = Graphs[w]
            exist_path, path = has_path(G, task)
            if exist_path:
                load[w] = load[w] + len(path) - 1
                can_route = True
                LightPath[index] = (path, w)
                for apk in range(len(path) - 1):
                    e = (path[apk + 1], path[apk])
                    Graphs[w].remove_edge(e[0], e[1])
                    Graphs[w].remove_edge(e[1], e[0])
                break
        if can_route == False:
            Graph = nx.DiGraph()
            for i in range(len(nodes)):
                Graph.add_node(nodes[i])
            for x, y in edges:  # edges
                Graph.add_edges_from([(x, y)])
                Graph.add_edges_from([(y, x)])
            Graphs.append(Graph.copy())
            G = Graphs[-1]
            path = nx.shortest_path(G, source=task[0], target=task[1])
            LightPath[index] = (path, len(Graphs) - 1)
            require.append([])
            load.append(0)
            load[-1] = load[-1] + len(path) - 1
            require[-1].append(index)
            for apk in range(len(path) - 1):
                e = (path[apk], path[apk + 1])
                Graphs[-1].remove_edge(e[0], e[1])
                Graphs[-1].remove_edge(e[1], e[0])
    return Graphs, LightPath, load, require, reduce_wave


if __name__ == "__main__":
    # 设置随机种子
    random.seed(2)
    np.random.seed(2)
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
        order = list(range(len(tasks)))

        Graphs, LightPath, load, require = FFD_RWA(order, tasks, edges, nodes)

        Graphs, LightPath, load, require, reduce_wave = Shift_paths(order, tasks, Graphs, require, load, LightPath)
        if len(Graphs) >= 2:
            Graphs, LightPath, load, require, reduce_wave = Mutate(order, tasks, Graphs, require, load, reduce_wave, LightPath, nodes, edges)
        end_time = time.time()  # 求解终止时间
        # print(end_time, start_time)
        times.append(end_time - start_time)
        ns = 0
        for i in range(len(LightPath)):
            ns = len(LightPath[i][0]) + ns - 1
        wave_link_length.append(ns)

        wave_num.append(len(Graphs) - len(reduce_wave))
    print(wave_link_length)
    print(sum(wave_link_length))
    print(sum(times))
    print(wave_num)
    print(sum(wave_num))