# -*- coding: utf-8 -*-
# Authors:   张, 刘剑锋
"""
实现FFD-RWA：连接请求按照最短路径的递减顺序排列，用FF-RWA方法进行路由。
"""
import numpy as np
import copy
import networkx as nx
import random
import time
import matplotlib
matplotlib.use('TkAgg')

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
    num = []      # 存储每个连接请求的最短路径数值
    Paths = []    # 存储每个连接请求的最短路径
    for task in tasks:
        exist_path, path = has_path(G, task)
        if exist_path:
            num.append(len(path)-1)
            Paths.append(path)
        else:  #找不到最短路径
            num.append(0)
            Paths.append([0])
    order = sorted(range(len(num)), key=lambda k: num[k], reverse=True)  #每个连接请求的最短路径按长度降序排列, 存储对应的任务请求序号
    Path = copy.deepcopy(Paths)
    for i in range(len(order)):
        Path[i] = Paths[order[i]]

    return order, Path



#连接请求按照最短路径的递减顺序排列，用FF-RWA方法进行路由。
def FFD_RWA(tasks, edges, nodes):
    """
    :param   tasks: 任务连接请求
    :param   edges: 边序列
    :param   nodes: 结点序列
    :return:Graphs: 每一个波长层的剩余边连接关系
         LightPath: 存储每条路径所及其所对应的波长
    """

    LightPath = []               # 存储每个连接请求的路径及其对应的波长层
    Graph = nx.Graph()
    for i in range(len(nodes)):
        Graph.add_node(nodes[i])
    for x, y in edges:  # edges::
        Graph.add_edges_from([(x, y)])
    Graphs = []
    Graphs.append(Graph.copy())
    order, Path = Paths_order(tasks, Graph)  # 得到连接请求的最短路径, 并按照降序进行排列
    require = [[]]

    for id in order:
        task = tasks[id]
        can_rounte = False
        for w in range(len(Graphs)):
            exist_path, path = has_path(Graphs[w], task)
            if exist_path:
                LightPath.append((path, w))
                require[w].append(id)  # 记录每个波长层所分配的任务请求
                can_rounte = True
                for apk in range(len(path) - 1):
                    e = (path[apk + 1], path[apk])
                    Graphs[w].remove_edge(e[0], e[1])
                break
        if can_rounte == False:
            Graphs.append(Graph.copy())
            path = nx.shortest_path(Graphs[-1], source=task[0], target=task[1])
            LightPath.append((path, w + 1))
            require.append([id])
            for apk in range(len(path) - 1):
                e = (path[apk], path[apk + 1])
                Graphs[-1].remove_edge(e[0], e[1])
    return Graphs, LightPath, order, require


if __name__ == "__main__":

    random.seed(2)  # 设置随机种子
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

        Graphs, LightPath, order, require = FFD_RWA(tasks, edges, nodes)

        end_time = time.time()  # 求解终止时间

        times.append(end_time - start_time)
        ns = 0
        for i in range(len(LightPath)):
            ns = len(LightPath[i][0]) + ns - 1
        wave_link_length.append(ns)

        wave_num.append(len(Graphs))
    print(wave_link_length)
    print(sum(wave_link_length))
    print(sum(times))
    print(wave_num)
    print(sum(wave_num))




