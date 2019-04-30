import math
import random
import bisect
import collections
from statistics import *

""" UTILITY FUNCTION LIBRARY """

def initGene(size: int):
    """Create single vector gene of size size"""
    gene = []
    for i in range(size):
        gene.append(getRandomNum(-1, 5, False))
    return gene

def initChromosome(population: int, geneLength: int):
    """Create list of genes at size population"""
    chromosome = []
    for i in range(population):
        chromosome.append(initGene(geneLength))
    return chromosome

def getGeneFitness(gene:list):
    """ Return the value of gene after fitness """
    value = 0
    for i in range(0, len(gene)):
        value += ((gene[0]**2) + (gene[1]**2) + (gene[2]**2))
    return value

def geneFitnessList(chromosome:list):
    """Return a list of gene fitness results"""
    fitResults = list()
    for gene in chromosome:
        fitResults.append(getGeneFitness(gene))
    return fitResults

def getFitness(lst:list):
    """Get list of fitness results for genes of chromosome"""
    fitResults:list = geneFitnessList(lst)
    bestFit, bestId = min((bestFit, bestId) for (bestId, bestFit) in enumerate(fitResults))
    worstFit, worstId = max((worstFit, worstId) for (worstId, worstFit) in enumerate(fitResults))
    averageFitness = mean(fitResults)
    std = stdev(fitResults)
    #print(f'\t\tBest Fitness:\t\t{bestFit} @{bestId}: {lst[bestId]}\n\t\tWorst Fitness:\t\t{worstFit} @{worstId}: {lst[worstId]}\n\t\tFitness Average:\t{averageFitness}\n')
    return bestFit, worstFit, averageFitness, std

def getRandomNum(minVal: float, maxVal: float, maximize: bool):
    """Create standard variation random number of range, max or min"""
    val1 = random.uniform(minVal, maxVal)
    val2 = random.uniform(minVal, maxVal)
    return round(max(val1,val2),4) if maximize else round(min(val1, val2),4)

def printChromosome(chromosome: list):
    """Print chromosome to stdout with formatting"""
    print(f'\t\n'.join(['\t'.join(['\t\t'+str(format(str(cell), '.6')) for cell in row]) for row in chromosome]))
    print(f'\t\t--------------------')
    getFitness(chromosome)
    return

def choice(population, weights):
    assert len(population) == len(weights)
    cdf_vals = cdf(weights)
    x = random.random()
    idx = bisect.bisect(cdf_vals, x)
    return population[idx]

def cdf(weights):
    total = sum(weights)
    result = []
    cumsum = 0
    for w in weights:
        cumsum += w
        result.append(cumsum / total)
    return result