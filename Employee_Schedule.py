# -*- coding: utf-8 -*-
"""
Created on Sun May 23 22:19:51 2021

@author: KUNAL
"""

#@author: KUNAL
# No days 7 total slot 7*3=21 no of employee E=7, Max munber of Employee per slot =2 

from dwave_qbsolv import QBSolv
import networkx as nx
from collections import defaultdict
from dwave.system import DWaveSampler,EmbeddingComposite
import dwave.inspector
from dwave.system import FixedEmbeddingComposite
from dwave.system import LeapHybridSampler
import dimod
from dimod import BinaryQuadraticModel

n_E=7
n_D=7
n_S=21
size= n_E*n_S
chainstrength=100
#G=nx.Graph()
#G.add_edges_from([(0,4),(0,5),(1,2),(1,6),(2,4),(3,7),(5,6),(6,7)])
Q=defaultdict(int)
lagrage= 5
def get_index(employee_index,shift_index):

    return employee_index *  n_S + shift_index

# Inverse of get_index - given a composite index in a 1D list, return the

# employee_index and day_index

def get_employee_and_day_and_shift(index):

    employee_index, shift_index = divmod(index, n_S)
    return employee_index, shift_index


#Contraint
for i in range(n_E):
    for s in range(n_S):
        ind=get_index(i,s)
        Q[(ind,ind)] += -3*lagrage
        for j in range(i+1,n_E):
            indj=get_index(j,s)
            Q[(ind,indj)] += 2*lagrage

#Objective 
for i in range(n_E):
    for s in range(n_S):
        ind=get_index(i,s)
        for j in range(i+1,n_E):
            indj=get_index(j,s)
            Q[(ind,indj)] += -1
            
sampler =EmbeddingComposite(DWaveSampler())
sampleset =sampler.sample_qubo(Q, chain_strength=chainstrength, num_reads=10000)
print(Q)
print(sampleset) 
dwave.inspector.show(sampleset)

isng = dimod.qubo_to_ising(Q,0)
print(isng[0])
print(isng[1])              

            
            
sampler =EmbeddingComposite(DWaveSampler())
response =sampler.sample_qubo(Q, chain_strength=chainstrength, num_reads=100)
bqm = BinaryQuadraticModel.from_qubo(Q, 9*lagrage)
sampler = LeapHybridSampler()
response = sampler.sample(bqm)
print(response)
samples = response.first.sample
energy = response.first.energy



print("Size ", size)
print("Energy ", energy)
print(samples)

# Graphics
for x in range(n_S):
    print() 
    sched = [get_employee_and_day_and_shift(j) for j in range(size) if samples[j] == 1]
    str_header_for_output1 = "Day: "
    str_header_for_output1 += " " * 6
    str_header_for_output1 += "             ".join(map(str, range(n_D)))
    
    print(str_header_for_output1)
    
    str_header_for_output = "Shift: "
    str_header_for_output += " " * 4
    str_header_for_output += "  ".join(map(str, range(n_S)))
    print(str_header_for_output)
    for n in range(n_E):
        str_row = ""
        for d in range(n_S):
            outcome = "X" if (n, d) in sched else " "
            if d > 9:
                outcome += " "
            str_row += "  " + outcome
        print("Employee ", n, str_row)
