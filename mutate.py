import math
import random
from statistics import *

""" MUTATE FUNCTION LIBRARY """

def mutateGaussian(list1:list, probM: int, sigma:float):
    """ Mutates allele by random offset, maximum of sigma strength """
    for index in range(0, len(list1)):
        if(random.uniform(0,10) <= (probM * 10)):
            temp = round(random.uniform(-1*sigma, sigma), 4)
            list1[index] += temp
    return list1

def mutateRandom(list1:list, probM: int):
    """ Reset Mutation: allele is replace by random val of initialize range: [-1,5] """
    for index in range(0, len(list1)):
        if(random.uniform(0,10) <= (probM * 10)):
            list1[index] = round(random.uniform(-1, 5), 4)
    return list1

