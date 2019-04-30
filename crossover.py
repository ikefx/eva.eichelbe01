import math
import random
from statistics import *

""" CROSSOVER FUNCTION LIBRARY """

def crossover(list1:list, list2:list):
    """ One point crossover: swap genes at random partition """
    assert len(list1) == len(list2), "genes do not have equal length for crossover()"
    lSize = len(list1)
    cutpoint = random.randint(1,lSize-1)
    child1 = []
    child2 = []
    for i in range(0,cutpoint):
        child1.append(list1[i])
        child2.append(list2[i])
    for j in range(cutpoint, lSize):
        child1.append(list2[j])
        child2.append(list1[j])
    
    return child1, child2

def crossoverShuffle(list1:list, list2:list):
    """ Improved One point crossover: shuffle allele before cross """
    assert len(list1) == len(list2), "genes do not have equal length for crossover()"
    lSize = len(list1)
    list1new = list1
    list2new = list2
    random.shuffle(list1new)
    random.shuffle(list2new)
    cutpoint = random.randint(1,lSize-1)
    child1 = []
    child2 = []
    for i in range(0,cutpoint):
        child1.append(list1new[i])
        child2.append(list2new[i])
    for j in range(cutpoint, lSize):
        child1.append(list2new[j])
        child2.append(list1new[j])
    
    return child1, child2
   
def crossoverUniform(parent1:list, parent2:list):
    """ Uniform Crossover: uniform random if each allele is swapped """
    assert len(parent1) == len(parent2), "genes do not have equal length for crossoverUniform()"
    child1 = list()
    child2 = list()
    for idx in range(0, len(parent1)):
        tempAllele1 = parent1[idx] if random.uniform(1,100) >= 50 else parent2[idx]
        tempAllele2 = parent1[idx] if random.uniform(1,100) >= 50 else parent2[idx]
        child1.append(tempAllele1)
        child2.append(tempAllele2)

    return child1, child2


def crossoverArithmetic(gene1:list, gene2:list, alpha:int):
    """ Arithmetic Crossover: Create children weighted from averages of x[i] """
    assert len(gene1) == len(gene2), "genes do not have equal length for crossoverArithmetic()"
    child1 = list()
    child2 = list()
    for idx in range(0, len(gene1)):
        child1.append( (alpha * gene1[idx]) + ((1 - alpha) * gene2[idx]) )
        child2.append( (1-alpha) * gene1[idx] + (alpha * gene2[idx]) )

    return child1, child2

  
