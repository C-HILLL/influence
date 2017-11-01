# -*- coding: utf-8 -*-
# Utilities
import sys, getopt
# Libraries
import numpy as np
from mpi4py import MPI

import graph_json_io

#################################
# -a # 阻尼系数,即α
# -i # 最大迭代次数
# -e # 确定迭代是否结束的参数,即ϵ
# -f # file_path
#################################

def dictSum(dict1, dict2, datatype):
    for item in dict2:
        if item in dict1:
            dict1[item] += dict2[item]
        else:
            dict1[item] = dict2[item]
    return dict1

dictSumOp = MPI.Op.Create(dictSum, commute=True)



class Calculator:
    
    def __init__(self, comm, a_factor, max_iterate, end_factor,file_path,isWeight=False):
        
        self.comm = comm
        self.comm_rank = comm.Get_rank()
        self.comm_size = comm.Get_size()
        # which rows&columns each node holds
        #行矩阵
        self.rows = None
        #列矩阵
        self.columns = None
        # master rank
        self.master = 0
        # height of the matrices
        
        self.a_factor = a_factor
        self.max_iterate = max_iterate
        self.end_factor = end_factor
        
        self.result = {}
        self.data_send_list = []
        self.columns = None
        self.rows = None
        self.indices = None
        if self.comm_rank == 0:
            
            self.full_columns = np.mat(self.read_graph(file_path,isWeight))
            self.height = self.full_columns.shape[0]
            self.full_rows = np.mat(np.ones((self.height,1)))/float(self.height)
            self.full_indices = range(self.height)
            if self.comm_size<2:
                print "can't run on master model!"
                self.comm.Abort()
            elif self.comm_size-1>self.height:
                self.length = 1
                self.step = self.height
                self.rest_length = 0
                self.rest_step = self.comm_size-1-self.height
            else:
                self.step = self.height%(self.comm_size-1)
                self.length = self.height/(self.comm_size-1)+1
                self.rest_length = self.length-1
                self.rest_step = self.comm_size-1-self.step  
            self.assign('colums')
        self.columns,self.indices = comm.scatter(self.data_send_list, root=0)

    def read_graph(self,filename,isWeight):
        simulate_matrix = None
        if True:    #暂时只处理无权重的图
            network = graph_json_io.read_json_file(file_path)
            nodes =  network.nodes()
            nodes.sort()
            simulate_matrix = [0]*len(nodes)
            for i in nodes:
                nbrs = [0]*len(nodes)
                nbrs_list = network.successors(i)
                
                nbrs_num = len(nbrs_list)
                for nbr in nbrs_list:
                    nbrs[nbr] = 1.0/nbrs_num
                simulate_matrix = np.column_stack((simulate_matrix,nbrs))

        return simulate_matrix[:,1:]
        
    def assign(self,cr):
        del self.data_send_list[:]
        if cr=='colums':
            self.data_send_list.append([None,None])
            for i in range(self.step):
                data_send=[self.full_columns[:,i*self.length:(i+1)*self.length],self.full_indices[i*self.length:(i+1)*self.length]]
                self.data_send_list.append(data_send)
            base = self.step*self.length
            for i in range(self.rest_step):
                data_send=[self.full_columns[:,base+i*self.rest_length:base+(i+1)*self.rest_length],self.full_indices[base+i*self.rest_length:base+(i+1)*self.rest_length]]
                self.data_send_list.append(data_send)
        elif cr=='rows':
            self.data_send_list.append([None,None])
            for i in range(self.step):
                data_send=[self.full_rows[i*self.length:(i+1)*self.length,:],self.full_indices[i*self.length:(i+1)*self.length]]
                self.data_send_list.append(data_send)
            base = self.step*self.length
            for i in range(self.rest_step):
                data_send=[self.full_rows[base+i*self.rest_length:base+(i+1)*self.rest_length,:],self.full_indices[base+i*self.rest_length:base+(i+1)*self.rest_length]]
                self.data_send_list.append(data_send)
         
    def calculate(self):
        
        if self.comm_rank == 0:
            preivous_rows = self.full_rows
        step = self.max_iterate
        while step>0:
            if self.comm_rank == 0:
                self.assign('rows')
            self.rows,self.indices = comm.scatter(self.data_send_list, root=0)
            if self.comm_rank > 0:
                result_dic = {}
                for i,index in enumerate(self.indices):
                    for r_i,r in enumerate(self.rows[i,:].getA1()):
                        for c_i,c in enumerate(self.columns[:,i].getA1()):
                            if (c_i,r_i) in result_dic:
                                result_dic[(c_i,r_i)] += c*r
                            else:
                                result_dic[(c_i,r_i)] = c*r
                #print self.comm_rank,result_dic
                self.result=result_dic  
            sumup = self.collect()
            
            if self.comm_rank == 0:
                self.result={}
                matrix = [(i,sumup[i]) for i in sumup]
                matrix.sort(key = lambda x:x[0][0]*self.height+x[0][1])
                #print matrix
                #print self.full_columns*self.full_rows
                self.full_rows = (1-self.a_factor)*np.mat(np.array([i[1] for i in matrix]).reshape((self.height,1)))+\
                            self.a_factor/self.height*np.mat(np.ones((self.height,1)))
                #print self.full_rows
                if sum(abs(self.full_rows-preivous_rows))>self.end_factor:
                    preivous_rows = self.full_rows
                else:
                    break
            step-=1
        print 'step:',step
        print self.full_rows
        self.comm.Abort()
        
    def collect(self):
        sumup = self.comm.reduce(self.result, op=dictSumOp)
        return sumup

if __name__ == "__main__":
    
    opts, _ = getopt.getopt(sys.argv[1:], "a:i:e:f:")
    a_factor=0
    max_iterate=100
    end_factor=0.00001
    file_path = '../data/weight_graph.g'
    
    for op, value in opts:
        if op == "-a":
            a_factor = float(value)
        elif op == "-i":
            max_iterate = int(value)
        elif op == "-e":
            end_factor = float(value)
        elif op == "-f":
            file_path = str(value)
            
    comm = MPI.COMM_WORLD
    node = Calculator(comm,a_factor,max_iterate,end_factor,file_path)
    node.calculate()