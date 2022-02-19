# -*- coding: utf-8 -*-
import collections

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


# 连接请求根据结点排序得到连接请求顺序
def ordering(nodes, tasks):
    """
    :param   nodes: 结点排列顺序
    :param   tasks: 任务连接请求
    :return: order: 任务连接请求排列顺序
    """
    order = []
    task_nums = len(tasks)  # 任务数
    task_index = list(range(task_nums))  # 任务连接请求索引
    node_nums = len(nodes)  # 结点数量
    node_index = np.arange(node_nums)  # 存储每个结点在结点排列中的位置
    for index in range(node_nums):
        node_index[nodes[index]] = index
    while len(order) != task_num:
        for i in nodes:
            same_i = []  # 存储起点相同的任务请求
            task_end_index = []  # 存储起点相同的任务请求的终点在结点排列中的位置
            j = 0
            while j < len(task_index):
                index = task_index[j]
                task = tasks[index]
                if task[0] == i and index not in order:
                    same_i.append(index)
                    task_end_index.append(node_index[task[1]])
                    task_index.pop(j)
                else:
                    j += 1
            # 起点相同的任务请求根据终点在结点排列中的位置进行排序
            for j in nodes:  # 应该按照nodes中的顺序排序，而不应该按照自然顺序
                for k in task_end_index:
                    if k == j:
                        order.append(same_i[task_end_index.index(k)])
                        break
    return order


# 得到任务请求顺序对应的评价函数
def Evaluation_function(order, tasks, edges, nodes):
    """
    :param order: 任务请求顺序
    :param tasks: 任务请求
    :param edges: 边序列
    :param nodes: 结点序列
    :return:
    """
    LightPath = []  # 存储每个连接请求的路径及其对应的波长层
    Graph = nx.Graph()
    for i in range(len(nodes)):
        Graph.add_node(nodes[i])
    for x, y in edges:  # edges::
        Graph.add_edges_from([(x, y)])

    Graphs = []
    Graphs.append(Graph.copy())
    load = [0]
    for id in order:
        can_rounte = False
        task = tasks[id]
        w = 0
        for w in range(len(Graphs)):
            G = Graphs[w]
            exist_path, path = has_path(G, task)
            if exist_path:
                load[w] = load[w] + len(path) - 1
                LightPath.append((path, w))
                can_rounte = True
                for apk in range(len(path) - 1):
                    e = (path[apk + 1], path[apk])
                    Graphs[w].remove_edge(e[0], e[1])
                break

        if can_rounte == False:
            Graphs.append(Graph.copy())
            G = Graphs[-1]
            path = nx.shortest_path(G, source=task[0], target=task[1])
            LightPath.append((path, w + 1))
            load.append(0)
            load[-1] = load[-1] + len(path) - 1
            for apk in range(len(path) - 1):
                e = (path[apk], path[apk + 1])
                Graphs[-1].remove_edge(e[0], e[1])

    return Graphs, LightPath, load


# 禁忌搜索算法
def RWA_Tabu(tasks, edges, nodes):
    """
    :param          tasks: 任务请求序列
    :param          edges: 边序列
    :param          nodes: 结点序列
    :return:
    """
    H = []  # 禁忌表
    opt_load = [INF_val]
    node_nums = len(nodes)  # 结点数量
    cur_solution = list(nodes)  # 当前解向量
    Neighbor = []  # 存储当前解向量的邻域解向量的索引
    for i in range(node_nums - 1):
        for j in range(i + 1, node_nums):
            Neighbor.append((i, j))
    while True:
        while True:
            loads_search = []
            Graphs_search = []
            LightPath_search = []
            loads_sum_search = []
            Neighbor_search = []
            # 在当前解向量不在禁忌表中的领域解向量进行随机搜索
            cur_Neighbor = Neighbor.copy()  # 存储当前解向量的未检查邻域解向量的索引
            for k in range(10):
                find = False
                if len(cur_Neighbor) == 0:
                    break
                item_index = np.random.choice(range(len(cur_Neighbor)))  # 随机选择一个邻居结点，使用random会改变后面的问题实例
                temp_solution = cur_solution.copy()
                item = cur_Neighbor[item_index]  # 修改item
                i,j=item
                cur_Neighbor.remove(item)
                temp_solution[i], temp_solution[j] = temp_solution[j], temp_solution[i]
                while True:  # 该邻居解在禁忌表中
                    if temp_solution not in H:
                        find = True
                        break
                    temp_solution[i], temp_solution[j] = temp_solution[j], temp_solution[i]  # 应该回滚为当前排序，而不是继续在修改后的排序中搜索
                    if len(cur_Neighbor) == 0:
                        break
                    item_index = np.random.choice(range(len(cur_Neighbor)))  # 随机选择一个邻居结点，使用random会改变后面的问题实例
                    temp_solution = cur_solution.copy()
                    item = cur_Neighbor[item_index]  # 修改item
                    i, j = item
                    cur_Neighbor.remove(item)
                    temp_solution[i], temp_solution[j] = temp_solution[j], temp_solution[i]

                if find:
                    H.append(temp_solution)
                    order = ordering(temp_solution, tasks)
                    Graphs, LightPath, load = Evaluation_function(order, tasks, edges, nodes)
                    loads_search.append(load)
                    Graphs_search.append(Graphs)
                    LightPath_search.append(LightPath)
                    loads_sum_search.append(sum(load))
                    Neighbor_search.append(item)
            if len(loads_sum_search) == 0:
                break
            optimal = loads_sum_search.index(min(loads_sum_search))  # 局部搜索最佳解
            load = loads_search[optimal]
            if sum(load) < sum(opt_load):  # 找到新的最优解
                opt_G = Graphs_search[optimal]
                opt_LightPath = LightPath_search[optimal]
                opt_load = loads_search[optimal]
                neighbor = Neighbor_search[optimal]
                i, j = neighbor
                cur_solution[i], cur_solution[j] = cur_solution[j], cur_solution[i]

        local_opt = True
        for item in Neighbor:
            temp_solution = cur_solution.copy()
            i, j = item
            temp_solution[i], temp_solution[j] = temp_solution[j], temp_solution[i]
            order = ordering(temp_solution, tasks)
            Graphs, LightPath, load = Evaluation_function(order, tasks, edges, nodes)
            if sum(load) < sum(opt_load):
                opt_G = Graphs
                opt_LightPath = LightPath
                opt_load = load
                cur_solution = temp_solution
                H.remove(cur_solution)
                local_opt = False
                break
        if local_opt:
            break

    return opt_G, opt_LightPath, opt_load


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

        Occupathions, LightPath, opt_load = RWA_Tabu(tasks, edges, nodes)

        end_time = time.time()  # 求解终止时间
        # print(end_time, start_time)
        times.append(end_time - start_time)
        ns = 0
        for i in range(len(LightPath)):
            ns = len(LightPath[i][0]) + ns - 1
        wave_link_length.append(ns)
        print(ns)
        wave_num.append(len(Occupathions))
    print(wave_link_length)
    print(sum(wave_link_length))
    print(sum(times))
    print(wave_num)
    print(sum(wave_num))
