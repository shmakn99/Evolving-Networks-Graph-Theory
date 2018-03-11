# -*- coding: utf-8 -*-
"""
Created on Sat Mar 10 15:03:00 2018

@author: Manik Sharma
"""

import networkx as nx
import random as rd
import re

class net:
    
    def __init__(self):
        self.graph=nx.Graph()
        self.m_gene=['','']
        self.rtn_gene=['','']
    
    def contruct_graph(self):
        routine=get_routine()
        N=get_magnitude()/len(routine)
        
        for tp in routine:
            if tp[0]=='EY':
                erdren(tp[1],N)
            else:
                balbert(tp[1],N,tp[2])
        
    def get_routine(self):
        
        lst=[]
        i=0
        while i<len(self.rtn_gene):
            if self.rtn_gene[i]=='A':
                for j in range(i,len(self.rtn_gene)):
                    if self.rtn_gene[j]=='Z':
                        lst.append((i,j))
                        break
            i+=1
        
        #print (gene)
            
        ret=[]
        
        for pt in lst:
            cpt=self.rtn_gene[pt[0]:pt[1]+1]
            b=cpt.count('B')
            e=cpt.count('E')
            x=cpt.count('X')
            y=cpt.count('Y')
            
            if ((b>0 and x>0) or (e>0 and y>0)):
                if (b>e):
                    if (b-x==0):
                        continue
                    else:
                        if max(c,d,f)==c:
                            ret.append(('BX',abs(b-x),1))
                        elif max(d,f)==d:
                            ret.append(('BX',abs(b-x),2))
                        else:
                            ret.append(('BX',abs(b-x),3))
                else:
                    if (e==0 or y==0):
                        continue
                    else:
                        ret.append(('EY',min(e/y,y/e)))
            else:
                continue
                
        return(ret)

    def erdren(self,p,n):
        if nx.number_of_nodes(self.graph)==0:
            for i in range(n):
                for j in range(i+1,n):
                    if rd.uniform(0,1)<p:
                        self.graph.add_edge(i,j)
        
        else:
            o_n=nx.number_of_nodes(self.graph)
            
            for i in range(n+o_n):
                for j in range(i,o_n):
                    self.graph.add_edge(i,j)
        
        
    def balbert(self,m,n,var):
        if var==1:
            dic=nx.degree_centrality(self.graph)
        if var==2:
            dic=nx.betweenness_centrality(self.graph)
        if var==3:
            dic==nx.clustering(self.graph)
        
        s=sum(dic.keys())
        
        for k in dic.keys():
            dic[k]=dic[k]/s
        
        if nx.number_of_nodes(self.graph)==0:
            i=0,j=1
            while (nx.number_of_nodes(self.graph)<m):
                self.graph.add_edge(i,j)
                i=j
                j+=1
            
            while (nx.number_of_nodes(self.graph)<n):
                prob=re.uniform(0,1)
                c_p=0
                for k in dic.keys():
                    c_p+=dic[k]
                    if prob<c_p:
                        self.graph.add_edge(k,nx.number_of_nodes(self.graph))
            
        else:   
            t_n=nx.number_of_nodes(self.graph)+n
            while (nx.number_of_nodes(self.graph)<n):
                prob=re.uniform(0,1)
                c_p=0
                for k in dic.keys():
                    c_p+=dic[k]
                    if prob<c_p:
                        self.graph.add_edge(k,nx.number_of_nodes(self.graph))
            
        
    def get_magnitude(self):
        lst=[]
        i=0
        while i<len(gene):
            if gene[i]=='A':
                for j in range(i,len(gene)):
                    if gene[j]=='Z':
                        lst.append((i,j))
                        break
            i+=1
        
        #print (gene)
        
        n= (len(lst))
        
        ret=[]
        
        for pt in lst:
            cpt=gene[pt[0]:pt[1]+1]
            b=cpt.count('B')
            e=cpt.count('E')
            x=cpt.count('X')
            y=cpt.count('Y')
            
            if ((b>0 and x>0) or (e>0 and y>0)):
                if (b>e):
                    if (b-x==0):
                        continue
                    else:
                        ret.append(('BX',len(cpt)/(b+x)))
                else:
                    if (e==0 or y==0):
                        continue
                    else:
                        ret.append(('EY',len(cpt)/(e+y)))
            else:
                continue
        
        f1=np.mean([tp[1] for tp in ret if tp[0]=='EY'])
        f2=np.mean([tp[1] for tp in ret if tp[0]=='BX'])
        #print (f1,f2,n)
        
        return (int(n*f1*f2))
    
def genesis(ntgsm,glm,glrtn):
    