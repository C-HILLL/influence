# -*- coding: utf-8 -*-
# Utilities
import sys, getopt
# Libraries
import numpy as np
#import scipy.sparse as sparse
from mpi4py import MPI

import graph_json_io

#################################
# -a # 阻尼系数,即α
# -i # 最大迭代次数
# -e # 确定迭代是否结束的参数,即ϵ
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
    
    def __init__(self, a_factor, max_iterate, end_factor,filename='weight_graph.g',isWeight=False):
        
        self.a_factor = a_factor
        self.max_iterate = max_iterate
        self.end_factor = end_factor
        
        self.columns = self.read_graph(filename,isWeight)
        
        self.height = self.columns.shape[0]
        
        self.rows = np.mat(np.ones((self.height,1)))/float(self.height)
        #print self.columns
        #print self.rows
        #print self.height
        
    def read_graph(self,filename,isWeight):
        simulate_matrix = None
        if True:    #暂时只处理无权重的图
            network = graph_json_io.read_json_file('data/'+filename)
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
     
    def calculate(self):
        preivous_rows = self.rows
        step = self.max_iterate
        while step>0:
            self.rows = (1-self.a_factor)*self.columns*self.rows+\
                            self.a_factor/self.height*np.mat(np.ones((self.height,1)))
            
            if sum(abs(self.rows-preivous_rows))>self.end_factor:
                preivous_rows = self.rows
            else:
                break
            step-=1
        print step
        print self.rows
if __name__ == "__main__":
    
    opts, _ = getopt.getopt(sys.argv[1:], "a:i:e:")
    a_factor=0
    max_iterate=100
    end_factor=0.00001
    
    for op, value in opts:
        if op == "-a":
            a_factor = float(value)
        elif op == "-i":
            max_iterate = int(value)
        elif op == "-e":
            end_factor = float(value)
            
    node = Calculator(a_factor,max_iterate,end_factor)
    node.calculate()