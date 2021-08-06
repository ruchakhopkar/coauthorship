# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 14:31:08 2021

@author: rucha
"""

from collections import Counter
from itertools import chain, combinations
import pandas as pd
import networkx as nx
from matplotlib.pyplot import figure

df=pd.read_csv("E://PTRC//coauthorship//2020_cell.csv")

#creating full names of authors
df['author_0']=(df['author_0_first_name'])+' '+(df['author_0_last_name'])
df['author_1']=(df['author_1_first_name'])+' '+(df['author_1_last_name'])

lst_authors=[]
for i in range(len(df)):
    authors=[]
    authors.append(df.iloc[i,-2])
    authors.append(df.iloc[i,-1])
    lst_authors.append(set(authors))
    
#get pairs of authors who worked together
meets = Counter(chain.from_iterable(combinations(line, 2) for line in lst_authors))

#convert into dataframe
adjacency_matrix= pd.DataFrame.from_dict(meets, orient='index').reset_index()

#rename column names
adjacency_matrix=adjacency_matrix.rename(columns={"index": "Authors", 0: "Count"})

#creating separate columns for authors
author0=[]
author1=[]
for i in range(len(adjacency_matrix)):
    a0,a1=adjacency_matrix.iloc[i,0]
    author0.append(a0)
    author1.append(a1)

#adding the separate authors and deleting the combined authors
adjacency_matrix['author0']=author0
adjacency_matrix['author1']=author1
adjacency_matrix.drop(columns=['Authors'], axis=1, inplace=True)

#creating a graph
G = nx.Graph()
G = nx.from_pandas_edgelist(adjacency_matrix, 'author0', 'author1',edge_attr='Count')
figure(figsize=(10, 8))
nx.draw_shell(G, with_labels=True)