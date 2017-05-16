# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

class pairedGene(object):
    def __init__(self):
        self.stats = stats
        self.distance = distance
        
import random
from scipy import stats
import pandas as pd
df = pd.read_excel("SampleData.xls", "Sheet1")

print(df)

#get random pair numbers from user
random_pairs = int(input('Please enter number of random pairs desired: '))

#create a dictionary to maintain info on what genes have been used
genesUsed = {}

i = 0
while i < random_pairs:
    gene1 = df.loc[randint(0,len(df)-1)]
    print(gene1)
    gene2 = df.loc[randint(0,len(df)-1)]
    print(gene2)
    if gene1[0]+ gene2[0] not in genesUsed:
        genesUsed = {gene1[0]+ gene2[0] : i}
        pairedGene.stats = stats.spearmanr(gene1[1:len(gene1)], gene2[1:len(gene2)])
        print(pairedGene.stats)
        i += 1

print(genesUsed)