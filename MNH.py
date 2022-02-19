# -*- coding: utf-8 -*-
# Authors:   张, 刘剑锋
"""
实现MNH：在路由问题中进行局部搜索（利用最短路径不唯一这个特点）
     路由问题：
	首先为所有的连接请求都分配一个最短路径。
	然后依次为连接请求重新分配最短路径，如果新的最短路径上链路的最大负载小于	之前的最短路径。
	循环往复，直到寻找不到新的最短路径为止。
    波长分配问题：
	将路径按照递减顺序排列，依次为其分配波长。

"""
import collections
import numpy as np
import random
import time
import copy


## 表示无穷大
INF_val = float("inf")

class node:
    """
    结点类
    """
    def __init__(self,val, path, dis):
        """
        结点类的初始化
        :param val:  结点序号
        :param path: 存储该结点到起点的当前最短路径
        :param dis:  存储该结点到起点的当前最短距离
        """
        self.val = val
        self.path = path
        self.dis = dis


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

# 得到每个任务请求的所有最短路径
def ALL_Shortest_Path(tasks, edges):
    """
     :param                tasks: 任务连接请求
     :param                edges: 边序列
    :return:  All_Shortest_paths: 每个任务连接请求所对应的所有最短路径
    """
    adjs = {}
    for e in edges:
        if e[0] not in adjs:
            adjs[e[0]] = []
        adjs[e[0]].append(e[1])
        if e[1] not in adjs:
            adjs[e[1]] = []
        adjs[e[1]].append(e[0])
    node_num = len(adjs)    #图中结点个数
    All_Shortest_paths = []
    for task in tasks:
        start_node  = task[0]                            # 起点
        end_node    = task[1]                            # 终点
        distance    = [float('inf')] * node_num          # 存储每个结点到起点的当前最短距离
        queue       = []                                 # 存储广度优先遍历的结点序列队列
        min_dis     = float('inf')                       # 记录所找到的所有起点到终点路径的最小距离
        queue.append(node(start_node, [start_node], 0))  # 广度优先遍历的起点结点
        all_min_path = []                                # 存储起点到终点的所有最短路径
        all_possible_paths = []                          # 存储起点到终点的所有可能最短路径
        while queue:
            s = queue.pop(0) # 出队列
            if s.val == end_node and s.dis <= min_dis:   # 更新并记录当前所找到的所有可能最短路径
                min_dis = s.dis
                all_possible_paths.append(s)
            elif s.dis >= min_dis:    # 剪枝不可能扩展为最短路径的局部路径
                continue
            else:
                for next_node in adjs[s.val]:
                    if distance[next_node] >= s.dis + 1: # 只广度优先搜索有可能找到最短路径的局部路径
                        queue.append(node(next_node, s.path + [next_node], s.dis + 1))
                        distance[next_node] = s.dis + 1

        for p in all_possible_paths:  # 找到所有最短路径
            if p.dis == min_dis:
                all_min_path.append(p.path)
        All_Shortest_paths.append(all_min_path)

    return All_Shortest_paths

# 得到每个任务请求的最小最大负载对应的路径
def Min_Load_Paths(tasks, edges):
    """
         :param           tasks: 任务连接请求
         :param           edges: 边序列
         :return: min_load_path: 每个任务请求所对应的最小最大负载路径
    """
    tasks_num          = len(tasks)                         # 任务请求数量
    weight_edges       = []
    All_Shortest_paths = ALL_Shortest_Path(tasks, edges)    # 得到每个任务连接请求的所有最短路径,
    min_load_path      = []
    min_indexs          = []                                # 记录任务请求所分配路径的索引
    for paths in All_Shortest_paths:                        # 每个任务请求的当前所分配路径初始化为第一个最短路径, 记录当前路径分配下每条边的负载数
        path = paths[0]
        min_load_path.append(path)
        min_indexs.append(0)
        for apk in range(len(path) - 1):
            if path[apk] < path[apk + 1]:
                edge = (path[apk], path[apk + 1])
            else:
                edge = (path[apk + 1], path[apk])
            weight_edges.append(edge)
    weight_edges = collections.Counter(weight_edges)  # 记录每条边的负载数
    improvement = 0
    while True:
        old_improvement = improvement
        for i in range(tasks_num):
            if len(All_Shortest_paths[i]) == 1:  # 只有一条最短路径, 无需重新分配
                continue
            min_load  = 1  # 记录所有最短路径中的最小最大负载值
            min_index = min_indexs[i]   #记录最小最大负载的对应路径索引
            path      = All_Shortest_paths[i][min_index]
            for apk in range(len(path) - 1):                # 计算原先所分配最短路径上的最大负载数
                if path[apk] < path[apk + 1]:
                    edge = (path[apk], path[apk + 1])
                else:
                    edge = (path[apk + 1], path[apk])
                if min_load < weight_edges[edge]:
                    min_load = weight_edges[edge]
                weight_edges[edge] -= 1
            for j in range(len(All_Shortest_paths[i])):  # 计算其余最短路径上的最大负载数, 并找到最小最大负载的对应路径
                if j != min_indexs[i]:
                  path = All_Shortest_paths[i][j]
                  load = 0
                  for apk in range(len(path) - 1):
                     if path[apk] < path[apk + 1]:
                        edge = (path[apk], path[apk + 1])
                     else:
                        edge = (path[apk + 1], path[apk])
                     if load < weight_edges[edge] + 1:
                        load = weight_edges[edge] + 1
                  if load < min_load:
                     min_load  = load
                     min_index = j
            if min_index != min_indexs[i]: #在循环之后再判断是否进行了修改，这样应该可以进一步提升运行速度？
                min_indexs[i] = min_index
                improvement += 1
            path = All_Shortest_paths[i][min_index]
            for apk in range(len(path) - 1):  # 计算更新路径上每条边的负载数
                if path[apk] < path[apk + 1]:
                    edge = (path[apk], path[apk + 1])
                else:
                    edge = (path[apk + 1], path[apk])
                weight_edges[edge] += 1
        if old_improvement == improvement:
            break
    for i in range(tasks_num):
        index = min_indexs[i]
        if index != 0:
            min_load_path[i] = All_Shortest_paths[i][index]
    return min_load_path



# 连接请求按照最短路径降序排列，依次分配波长，如果当前所有波长层都无法路由，就添加新的波长层
def LFFP(order, Path, edges):
    """
    :param          order: 任务请求的序号
    :param           Path: 任务请求的最短路径
    :param          edges: 边序列
    :return: Occupathions: 每个波长层中边的占用情况
                LightPath: 存储每条路径所及其所对应的波长
                  require: 存储每个波长层中的连接请求
    """
    can_route = True  # 标识路径在该波长层中是否存在
    w = 0             # 标识波长层序号
    LightPath = []    # 存储每个连接请求的路径及其对应的波长层
    Occupathions = [] # 存储每个波长层中边的占用情况
    occupathion  = {}
    for e in edges[:]:
        occupathion[e] = False
    Occupathions.append(occupathion.copy())
    require = [[]]
    for i in range(len(Path)):
        path = Path[i]
        Edge = [] # 存储路径中所设计的边
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
            if can_route == True: # 分配波长
                LightPath.append((path, w)) # 记录每条路径对应的波长层
                require[w].append(order[i]) # 记录每个波长层所分配的任务请求
                for e in Edge:
                   Occupathions[w][e] = True
                break
        if can_route == False: # 如果当前所有波长层都无法路由，就添加新的波长层
            Occupathions.append(occupathion.copy())
            LightPath.append((path, w+1))
            require.append([order[i]])
            for e in Edge:
                Occupathions[-1][e] = True

    return Occupathions, LightPath, require

"""
  MNH：在路由问题中进行局部搜索（利用最短路径不唯一这个特点）
    路由问题：
	首先为所有的连接请求都分配一个最短路径。
	然后依次为连接请求重新分配最短路径，如果新的最短路径上链路的最大负载小于	之前的最短路径。
	循环往复，直到寻找不到新的最短路径为止。
     波长分配问题：
	将路径按照递减顺序排列，依次为其分配波长。
"""

def MNH(tasks, edges):
    """
       :param          tasks: 任务连接请求
       :param          edges: 边序列
       :return: Occupathions: 每个波长层中边的占用情况
                   LightPath: 存储每条路径所及其所对应的波长
                     require: 存储每个波长层中的连接请求
    """
    min_load_path = Min_Load_Paths(tasks, edges) # 得到每个任务连接请求所对应的所有最短路径

    num = [len(min_load_path[i]) -1 for i in range(len(min_load_path))]

    order = sorted(range(len(num)), key=lambda k: num[k], reverse=True)  # 每个任务连接请求的最短路径按长度降序排列, 存储对应的任务请求序号
    Path = [[] for _ in range(len(num))]
    for i in range(len(order)):
        Path[i] = min_load_path[order[i]]
    Occupathions, LightPath, require = LFFP(order, Path, edges)

    return Occupathions, LightPath, require

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

    for task_num in range(1, Max_num_tasks+1 ):

        tasks = generator_task(task_num, node_nums)  # 生成任务

        start_time = time.time()  # 求解起始时间

        Occupathions, LightPath, require = MNH(tasks, edges)

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
