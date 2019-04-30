import math
import random
import bisect
import collections
import numpy as np
from statistics import *

from util import *

""" SELECTION FUNCTION LIBRARY """

def selectionTournament(chromosome:list, k:int, P:float):
    """ K-Way Tournament Selection: select K random individuals from population, P probability keep best """
    assert k > 0, "selectionTournament requires a k value greater than 0"
    outlist = list()
    tempC = chromosome
    for j in range(0, len(tempC)):
       random.shuffle(tempC)
       genekeep = list()
       tempfit = float("inf")
       new = tempC[0:k]
       #tempC[0:k] = []
       for i in range(0, len(new)):
            if getGeneFitness(new[i]) < tempfit:
                if random.uniform(0,1) >= P:
                    genekeep = new[i]
                    tempfit = getGeneFitness(new[i])
                else:
                    genekeep = tempC[i]
       outlist.append(genekeep)
    return outlist

def selectionRankingLinear(chromosome:list, max:int):
    """ Linear Ranking Selection: from max create a ranking offset to sample """
    # printChromosome(chromosome)
    lSize = len(chromosome)
    rankList = list()
    ranks = list()
    
    for x in range(lSize):
        rankList.append(x) 
    # gauss algorithm
    rankSum = len(rankList) * ((rankList[0]+rankList[len(rankList)-1]) / 2)
    # normalize to fitness sum
    for x in range(lSize):
        rankList[x] = ( rankList[x] ) / rankSum
        ranks.append(x)
   
    # verify
    assert len(rankList) * ((rankList[0]+rankList[len(rankList)-1]) / 2) == 1 
    counts = collections.defaultdict(int)
    for i in range(len(chromosome)):
        counts[choice(ranks, rankList)] += 1
      
    # build new chromosome
    chrm01 = list()
    for entry in counts:
        print(f'Rank:{entry} Freq:{counts[entry]}')
        for x in range(counts[entry]):
            chrm01.append(chromosome[entry])
    chrm01 = selectionSimple(chrm01, lSize)



    return chrm01

def selectionRankingFitness(chromosome:list):
    """ Fitness Ranking Selection: y = 0.5x linear line ranking by fitness """
    lSize = len(chromosome)
    rankList = list()
    ranks = list()  
    for x in range(lSize):
        rankList.append(x+1)       # +1 allows r0 to have a probability != 0
    rankSum = len(rankList) * ((rankList[0]+rankList[len(rankList)-1]) / 2)
    
    for x in range(lSize):
        rankList[x] = ( rankList[x] ) / rankSum
        ranks.append(x)

    assert len(rankList) * ((rankList[0]+rankList[len(rankList)-1]) / 2) == 1
    counts = collections.defaultdict(int)
    for i in range(len(chromosome)):
        counts[choice(ranks, rankList)] += 1
      
    chrm01 = list()
    for entry in counts:
        for x in range(counts[entry]):
            chrm01.append(chromosome[entry])
    chrm01 = selectionSimple(chrm01, lSize)
    chrm01 = selectionSimple(chrm01, lSize)
    return chrm01

def selectionProportional(chromosome:list, size):
    """ Proportional Selection: Better fitness gets greater selection chance """
    assert len(chromosome) == size, "Error: selectionProportional size parameter incorrect"
    flist = list()
    for gene in chromosome:
        flist.append(getGeneFitness(gene))
    # gauss algorithm
    flistSum = len(flist) * ((flist[0]+flist[len(flist)-1]) / 2)
    problist = list()
    for gene in chromosome:
        fit = getGeneFitness(gene)
        problist.append(fit/flistSum)  
    idxlist = list()
    for idx in range(0, len(chromosome)):
        idxlist.append(idx)
    newChromo = list()
    newlist = np.random.choice(idxlist, len(chromosome), problist)
    outlist = list()
    for idx in range(0, len(chromosome)):
        outlist.append(chromosome[newlist[idx]])
    return outlist

def selectionTruncate(chromosome:list, percentile:int):
    """ Truncation Selection: sample from the best of a % of the population """
    chromosome.sort(key=sum)
    partitionIdx = int(math.ceil((percentile * len(chromosome)) / 100))
    newlist = list()
    outlist = list()
    for idx in range(0, partitionIdx):
        newlist.append(chromosome[idx])
    for idx in range(0, len(chromosome)):
        outlist.append(random.choice(newlist))
    return outlist

def selectionSimple(chromosome:list, size:int):
    """ Simple selection algorithm, simply shuffles the genes of a given chromosome 1:1 """
    assert size == len(chromosome), "invalid selectionSimple() size paramater"
    return random.sample(chromosome, size)
