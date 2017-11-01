# -*- coding: utf-8 -*-
# Utilities
import sys, getopt
# Libraries
from mpi4py import MPI
import numpy as np
from scipy.stats import bernoulli
import graph_json_io
import collections

#################################
# -s # seed_num s
# -t # simulate_time t
# -f # file_path
#################################

"""
#规约函数，取平均值
def dictSum(dict1, dict2, datatype):
    for item in dict2:
        if item in dict1:
            dict1[item] += dict2[item]
        else:
            dict1[item] = dict2[item]
    return dict1

dictSumOp = MPI.Op.Create(dictSum, commute=True)
"""


class Calculator:
    
    def __init__(self, comm, seed_num, simulate_time, file_path, isWeight=True):
        
        self.comm = comm
        self.comm_rank = comm.Get_rank()
        self.comm_size = comm.Get_size()
        
        self.network = None
        self.simulate_data = {}
        self.simulate_time_per_node = simulate_time/(self.comm_size-1)+1
        self.seed_nodes = []
        self.nodes_influence = 0
        
        if self.comm_rank == 0:
            self.seed_num = seed_num
            self.simulate_time = self.simulate_time_per_node*(self.comm_size-1)
            self.network = self.read_graph(file_path,isWeight)
            
        # master rank
        self.master = 0
        
        self.network = comm.bcast(self.network if self.comm_rank == 0 else None, root=0)
        print self.comm_rank,self.network.number_of_nodes(),self.simulate_time_per_node
        #if self.comm_rank > 0:
        #    print self.influence_simulate_IC([],[0],self.simulate_time_per_node)
            
    def read_graph(self,file_path,isWeight):
        
        if True:    #暂时只处理有概率权重图
            return graph_json_io.read_json_file(file_path)


        #单一节点的数据模拟
    def influence_simulate_IC(self,A,B,s_time):
        """
        A:上轮激活节点
        B:本轮应激活节点
        s_time:本轮应模拟次数
        """
        simulate_process = {}   
        #已激活的节点
        C = A + B
        
        if s_time == 1:
            #迭代更新产生激活节点
            iB = B[:]
            #迭代更新本次激活的节点
            iD = []
            while True:
                for n in iB:
                    #print n
                    for nbr in self.network.successors(n):
                        edge_weight = self.network.edge[n][nbr]['weight']
                        if nbr not in C and edge_weight>0:
                            if (n,nbr) not in self.simulate_data:
                                self.simulate_data[(n,nbr)] = list(bernoulli.rvs(edge_weight, loc=0, size=50,random_state=None))     
                            elif len(self.simulate_data[(n,nbr)])<s_time:
                                self.simulate_data[(n,nbr)] += list(bernoulli.rvs(edge_weight, loc=0, size=50,random_state=None))     
                                
                            isActivate = self.simulate_data[(n,nbr)][0:s_time][0]
                            self.simulate_data[(n,nbr)] = self.simulate_data[(n,nbr)][s_time:]
                            if isActivate:    
                                iD.append(nbr)
                                C.append(nbr)
                if len(iD) == 0:
                    break
                else:
                    iB = iD[:]
                    #print iD
                    del iD[:]       
            #返回本次模拟的长度
            return len(C)
        
        for n in B:
            for nbr in self.network.successors(n):
                edge_weight = self.network.edge[n][nbr]['weight']
                if nbr not in C and edge_weight >0:
                    if (n,nbr) not in self.simulate_data:
                        self.simulate_data[(n,nbr)] = list(bernoulli.rvs(edge_weight, loc=0, size=2000,random_state=None))     
                    elif len(self.simulate_data[(n,nbr)])<s_time:
                        self.simulate_data[(n,nbr)] += list(bernoulli.rvs(edge_weight, loc=0, size=2000,random_state=None))     
                            
                    simulate_list = self.simulate_data[(n,nbr)][0:s_time]
                    if sum(simulate_list)>0:
                        if nbr in simulate_process:
                            simulate_process[nbr][simulate_process[nbr]+simulate_list>0]=1
                        else:
                            simulate_process[nbr] = np.array(simulate_list)
      
                    self.simulate_data[(n,nbr)] = self.simulate_data[(n,nbr)][s_time:]
        
        #本次激活的节点集合
        D = simulate_process.keys()
        if len(D)==0:
            return len(C)*s_time
        #print D
        simulate_matrix = [0]*s_time
        for node_simulate in simulate_process.values():
            simulate_matrix = np.column_stack((simulate_matrix,node_simulate))
        #print simulate_matrix
        count_lst=[]
        for i in range(s_time):
            activate_nodes = filter(lambda (x,y):x>0, zip(simulate_matrix[i,1:],D))
            activate_nodes = [y for (x,y) in activate_nodes]
            count_lst.append(tuple(activate_nodes))
        count_dict = collections.Counter(count_lst)
        #print count_dict
        activate_node_num = 0
        for k in count_dict:
            if len(k)==0:
                activate_node_num += len(C)*count_dict[k]
            else:
                activate_node_num += self.influence_simulate_IC(C,list(k),count_dict[k])
        return activate_node_num

    def CELF(self):
        
        if self.comm_rank > 0:
            while True:
                B = comm.bcast(self.seed_nodes if self.comm_rank == 0 else None, root=0)
                self.nodes_influence = self.influence_simulate_IC([],B,self.simulate_time_per_node)
                self.comm.reduce(self.nodes_influence, root=0,op=MPI.SUM)
            
        if self.comm_rank == 0:
            SEED=[]
            candidate_list = []
            # init candidate_list
            for n in self.network.nodes():
                self.seed_nodes = [n]
                self.comm.bcast(self.seed_nodes if self.comm_rank == 0 else None, root=0)
                increase = self.comm.reduce(self.nodes_influence, root=0,op=MPI.SUM)/float(self.simulate_time)
                candidate_list.append([n,increase,0])
            candidate_list.sort(key=lambda x:x[1],reverse=True)
            
            #先使用自带排序方法，后续考虑是否优化
            seed_value = 0
            for i in range(seed_num):
                while candidate_list[0][2]<i:
                    self.seed_nodes = SEED+[candidate_list[0][0]]
                    self.comm.bcast(self.seed_nodes if self.comm_rank == 0 else None, root=0)
                    candidate_list[0][1] = self.comm.reduce(self.nodes_influence, root=0,op=MPI.SUM)/float(self.simulate_time)-seed_value
                    candidate_list[0][2]=i
                    candidate_list.sort(key=lambda x: x[1], reverse=True)
            
                SEED.append(candidate_list[0][0])
                seed_value += candidate_list[0][1]
                del candidate_list[0]
                print SEED,seed_value 
                
        self.comm.Abort()
    
if __name__ == "__main__":
    
    opts, _ = getopt.getopt(sys.argv[1:], "s:t:f:")
    seed_num=5
    simulate_time=100
    file_path = '../data/weight_graph.g'
    
    for op, value in opts:
        if op == "-s":
            seed_num = int(value)
        elif op == "-s":
            simulate_time = int(value)
        elif op == "-f":
            file_path = str(value)
            
    comm = MPI.COMM_WORLD
    node = Calculator(comm,seed_num,simulate_time,file_path)
    node.CELF()