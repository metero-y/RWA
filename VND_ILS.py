# -*- coding: utf-8 -*-
import heapq

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
import copy

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

def shortestpath(order, tasks, edges, nodes):
    """
    :param order:  任务连接请求序号
    :param tasks:  任务连接请求
    :param edges:  边序列
    :param nodes:  结点序列
    :return:  sp:  任务请求的最短路径
    """
    sp = []
    Graph = nx.DiGraph()
    for i in range(len(nodes)):
        Graph.add_node(nodes[i])
    for x, y in edges:  # edges::
        Graph.add_edges_from([(x, y)])
        Graph.add_edges_from([(y, x)])
    for id in order:
        task = tasks[id]
        path = nx.shortest_path(Graph, source=task[0], target=task[1])
        sp.append(path)
    return sp

#连接请求按照最短路径的递减顺序排列，用FF-RWA方法进行路由。
def FFD_RWA(order, tasks, edges, nodes):
    """
    :param   order: 任务连接请求序号
    :param   tasks: 任务连接请求
    :param   edges: 边序列
    :param   nodes: 结点序列
    :return:Graphs: 每一个波长层的剩余边连接关系
         LightPath: 存储每条路径所及其所对应的波长
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

#执行微扰
def ils(orders, taskss, Graphss, loads, nodess, LightPaths, requires, lambda_t):

    R = []
    order, tasks, Graphs, load, nodes, LightPath, require = copy.deepcopy(orders), copy.deepcopy(taskss), copy.deepcopy(
        Graphss), copy.deepcopy(loads), copy.deepcopy(nodess), copy.deepcopy(LightPaths), copy.deepcopy(requires)
    for w in range(len(Graphs)):
        if w == lambda_t:
            continue
        # 随机选一个任务连接请求
        r = list(np.random.choice(require[w], 1, replace=False))[0]
        path = LightPath[order.index(r)][0]
        R.append((r, w))
        for apk in range(len(path) - 1):
            e = (path[apk + 1], path[apk])
            Graphs[w].add_edges_from([(e[0], e[1])])
            Graphs[w].add_edges_from([(e[1], e[0])])
        require[w].remove(r)
        load[w] = load[w] - len(path) + 1
    R_index = list(np.random.choice(range(len(R)), len(R), replace=False))
    change = False
    for item in R_index:
        id, wave = R[item]
        for w in range(len(Graphs)):
            if w == lambda_t:
                continue
            G = Graphs[w]
            task = tasks[id]
            exist_path, path = has_path(G, task)
            if exist_path:
                if LightPath[order.index(id)][0] != path or wave!=w:#同一层中不同路径也应该是进行了改变
                    change =True
                load[w] = load[w] + len(path) - 1
                LightPath[order.index(id)] = (path, w)
                for apk in range(len(path) - 1):
                    e = (path[apk + 1], path[apk])
                    Graphs[w].remove_edge(e[0], e[1])
                    Graphs[w].remove_edge(e[1], e[0])
                require[w].append(id)
                break

        if id not in require[w]:
            return [], [], [], [], False
    return Graphs, LightPath, load, require, change

def VND(order, tasks, Graphs, load, nodes, LightPath, require, sp):
    """
    :param order:      任务连接请求顺序
    :param tasks:      任务连接请求
    :param Graphs:     每个波长层中的剩余图
    :param load:       每个波长层中被占用的波长链路数
    :param nodes:      结点序列
    :param LightPath:  每个连接请求对应的光路径
    :param require:    每个波长层中中被分配链路路径的任务请求
    :param sp:         每个任务请求对应的最短路径
    :return:
    """
    require_num = []              # 存储每个波长层中的任务请求数及其波长层序号
    for i in range(len(require)):
        require_num.append((len(require[i]), i))
    heapq.heapify(require_num)    # 根据每个波长层中的任务请求数进行推排序
    lambda_t = require_num[0][1]  # 找到利用链路最小的层

    can_list = []
    #print(require[lambda_t])
    for i in require[lambda_t]:
        can_list.append((-len(sp[order.index(i)]), i)) # 得到对应波长层中连接请求以最短路径长度递减顺序的排列
    #print(len(tasks))
    heapq.heapify(can_list)  # 最短路径降序排列
    while can_list != []:
        temple = heapq.heappop(can_list) # 得到要转移的第一个连接请求
        r_t = temple[1]  # 任务请求序号
        #print(can_list, r_t)
        transfer = False
        for w in range(len(Graphs)):  # 第一邻域

            if w == lambda_t:  # 不同于要被转移的波长层
                continue
            G = Graphs[w]
            task = tasks[r_t]
            exist_path, path = has_path(G, task)
            if exist_path:
                #print(require[lambda_t],1)
                transfer = True
                r_t_path = LightPath[order.index(r_t)][0] # 转移路径在原来波长层中的路径
                for apk in range(len(r_t_path) - 1): # 恢复转移路径在原来波长层中的边
                    e = (r_t_path[apk + 1], r_t_path[apk])
                    Graphs[lambda_t].add_edges_from([(e[0], e[1])])
                    Graphs[lambda_t].add_edges_from([(e[1], e[0])])
                load[lambda_t] = load[lambda_t] - len(r_t_path) + 1
                load[w] = load[w] + len(path) - 1
                LightPath[order.index(r_t)] = (path, w)
                #print(require[lambda_t], r_t)
                require[lambda_t].remove(r_t)
                require[w].append(r_t)
                for apk in range(len(path) - 1):
                    e = (path[apk + 1], path[apk])
                    Graphs[w].remove_edge(e[0], e[1])
                    Graphs[w].remove_edge(e[1], e[0])
                break

        if transfer == False:  # 第一领域没遇到可行波长层, 转移到第二邻域, 即对于每一个不同于lambda_t的波长层, 尝试将其每个连接请求重新分配给除lambda_t之外的其他波长层。
            for w_1 in range(len(Graphs)):  # 该层不再还原

                if w_1 == lambda_t: # 不同于要被转移的波长层
                    continue
                for id in require[w_1]:  # 任务连接请求下标
                    for w_2 in range(len(Graphs)):  # 重新路由
                        if w_2 == lambda_t or w_2 == w_1:
                            continue
                        G = Graphs[w_2]
                        task = tasks[id]
                        exist_path, path = has_path(G, task)
                        if exist_path:
                            for apk in range(len(path) - 1):
                                e = (path[apk + 1], path[apk])
                                Graphs[w_2].remove_edge(e[0], e[1])
                                Graphs[w_2].remove_edge(e[1], e[0])
                            ori_path = LightPath[order.index(id)][0]
                            for apk in range(len(ori_path) - 1):  # 恢复路径在原来波长层中的边
                                e = (ori_path[apk + 1], ori_path[apk])
                                Graphs[w_1].add_edges_from([(e[0], e[1])])
                                Graphs[w_1].add_edges_from([(e[1], e[0])])  # 添加回来原来的层，然后看看r_t能不能路由进来
                            load[w_1] = load[w_1] - len(ori_path) + 1
                            load[w_2] = load[w_2] + len(path) - 1
                            LightPath[order.index(id)] = (path, w_2)
                            require[w_1].remove(id)
                            require[w_2].append(id)
                            break
                    G = Graphs[w_1]
                    task = tasks[r_t]
                    exist_path, path = has_path(G, task)
                    if exist_path:
                        #print(require[lambda_t], 2)
                        transfer = True
                        r_t_path = LightPath[order.index(r_t)][0]  # 转移路径在原来波长层中的路径
                        for apk in range(len(r_t_path) - 1):  # 恢复转移路径在原来波长层中的边
                            e = (r_t_path[apk + 1], r_t_path[apk])
                            Graphs[lambda_t].add_edges_from([(e[0], e[1])])
                            Graphs[lambda_t].add_edges_from([(e[1], e[0])])
                        load[lambda_t] = load[lambda_t] - len(r_t_path) + 1
                        load[w_1] = load[w_1] + len(path) - 1
                        LightPath[order.index(r_t)] = (path, w_1)
                        require[lambda_t].remove(r_t)
                        require[w_1].append(r_t)
                        for apk in range(len(path) - 1):
                            e = (path[apk + 1], path[apk])
                            Graphs[w_1].remove_edge(e[0], e[1])
                            Graphs[w_1].remove_edge(e[1], e[0])
                        break
                if transfer:
                    break
        if transfer == False:  # 第三邻域， 寻找最短路径长度小于被转移路径长度的其他层中的连接请求, 与被转移路径进行交换
            for i in range(len(sp)):

                if i == order.index(r_t) or sp[i] >= sp[order.index(r_t)] or LightPath[i][1] == lambda_t:
                    continue
                w = LightPath[i][1]
                r_t_path = LightPath[order.index(r_t)][0]  # 转移路径在原来波长层中的路径
                for apk in range(len(r_t_path) - 1):  # 恢复转移路径在原来波长层中的边
                    e = (r_t_path[apk + 1], r_t_path[apk])
                    Graphs[lambda_t].add_edges_from([(e[0], e[1])])
                    Graphs[lambda_t].add_edges_from([(e[1], e[0])])

                ori_path = LightPath[i][0]
                for apk in range(len(ori_path) - 1):  # 恢复路径在原来波长层中的边
                    e = (ori_path[apk + 1], ori_path[apk])
                    Graphs[w].add_edges_from([(e[0], e[1])])
                    Graphs[w].add_edges_from([(e[1], e[0])])  # 添加回来原来的层，然后看看r_t能不能路由进来

                # 看交换是否可行
                G_w = Graphs[w]
                task_r_t = tasks[r_t]
                exist_r_t_path, new_r_t_path = has_path(G_w, task_r_t)
                if exist_r_t_path:
                    G_r_t = Graphs[lambda_t]
                    task_i = tasks[order[i]]
                    exist_i_path, i_path = has_path(G_r_t, task_i)
                    if exist_i_path: # 交换路径成功
                        transfer = True
                        #print(require[lambda_t], 3)
                        load[lambda_t] = load[lambda_t] - len(r_t_path) + 1
                        load[w] = load[w] + len(new_r_t_path) - 1
                        LightPath[order.index(r_t)] = (new_r_t_path, w)
                        require[lambda_t].remove(r_t)
                        require[w].append(r_t)
                        for apk in range(len(new_r_t_path) - 1):
                            e = (new_r_t_path[apk + 1], new_r_t_path[apk])
                            Graphs[w].remove_edge(e[0], e[1])
                            Graphs[w].remove_edge(e[1], e[0])
                        load[lambda_t] = load[lambda_t] + len(i_path) - 1
                        load[w] = load[w] - len(ori_path) + 1
                        LightPath[i] = (i_path, lambda_t)
                        require[lambda_t].append(order[i])
                        require[w].remove(order[i])

                        for apk in range(len(i_path) - 1):
                            e = (i_path[apk + 1],  i_path[apk])
                            Graphs[lambda_t].remove_edge(e[0], e[1])
                            Graphs[lambda_t].remove_edge(e[1], e[0])
                        heapq.heappush(can_list, (-len(sp[i]), order[i]))  # 最短路径降序排列
                        break
                for apk in range(len(r_t_path) - 1):  # 恢复转移路径在原来波长层中的边
                    e = (r_t_path[apk + 1], r_t_path[apk])
                    Graphs[lambda_t].remove_edge(e[0], e[1])
                    Graphs[lambda_t].remove_edge(e[1], e[0])

                for apk in range(len(ori_path) - 1):  # 恢复路径在原来波长层中的边
                    e = (ori_path[apk + 1], ori_path[apk])
                    Graphs[w].remove_edge(e[0], e[1])
                    Graphs[w].remove_edge(e[1], e[0])

        if transfer == False:  # 第三邻域也不行, 执行微扰
            threshold_value = 3
            times = 0
            while times <= threshold_value:
                times += 1
                Graphs_1, LightPath_1, load_1, require_1, change = ils(order, tasks, Graphs, load, nodes, LightPath, require, lambda_t)

                if Graphs_1 and change:
                    Graphs, LightPath, load, require = Graphs_1, LightPath_1, load_1, require_1

                    heapq.heappush(can_list, (-len(sp[order.index(r_t)]), r_t))
                    break
                if times == 4:
                    return Graphs, LightPath, load, require

    return Graphs, LightPath, load, require



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
        sp = shortestpath(order, tasks, edges, nodes) # 求解各个任务连接请求的最短路径
        Graphs, LightPath, load, require = FFD_RWA(order, tasks, edges, nodes)  # 利用FFD方法得到初始解
        Occupathions, LightPath, load,require = VND(order, tasks, Graphs, load, nodes, LightPath, require, sp)

        end_time = time.time()  # 求解终止时间
        # print(end_time, start_time)
        times.append(end_time - start_time)
        ns = 0
        for i in range(len(LightPath)):
            ns = len(LightPath[i][0]) + ns - 1
        wave_link_length.append(ns)
        wave_num.append(len(Occupathions))
    print(wave_link_length)
    print(sum(wave_link_length))
    print(sum(times))
    print(wave_num)
    print(sum(wave_num))