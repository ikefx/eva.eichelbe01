import math
import random
from prettytable import PrettyTable
from statistics import *
from selection import *
from mutate import *
from crossover import *
from util import *

evals:int = 0   # global value for number of evals of program

def main():
    """ PROJECT 02 """
    N = 10              # population size
    Pc = 0.8            # probability crossover
    Pm = 0.1            # probability mutation
    sigma = 1           # strength of mutation
    maxGenerations = 30 # number of generations

    alpha = 0.75        # crossoverArithmetic()         Alg Parameter
    percentile = 50     # selectionTruncate()           Alg Parameter
    k = 2               # selectionTournament() k       Alg Parameter
    P = 0.40            # selectionTournament() P       Alg Parameter
    
    listBOR = list()
    listAVG = list()
    listSTD = list()
    #1
    # init first generation
    chrm01 = initChromosome(N, 3)
    chrm01 = selectionSimple(chrm01, 10)    # 1:1 resample
    for gene in range(0, len(chrm01)):
        if(gene % 2 == 0):
            if(random.uniform(1, 10) <= (Pc * 10)):
                tempTuple = crossover(chrm01[gene], chrm01[gene+1])
                chrm01[gene] = mutateRandom(tempTuple[0], sigma)
                chrm01[gene+1] = mutateRandom(tempTuple[1], sigma)
    print(f'\tSimulation 1:')
    printChromosome(chrm01)

    fitness = list()
    avgList = list()
    for i in range(0, maxGenerations):
        "Do 30 runs"    
        selectionTournament(chrm01, k, P)
        for gene in range(0, len(chrm01)):
            if(gene % 2 == 0):
                if(random.uniform(1, 10) <= (Pc * 10)):
                    tempTuple = crossoverShuffle(chrm01[gene], chrm01[gene+1])
                    chrm01[gene] = mutateRandom(tempTuple[0], Pm)
                    chrm01[gene+1] = mutateRandom(tempTuple[1], Pm)
        fitness.append(getFitness(chrm01))
        minFit = float("inf")
        maxFit = 0
        for stat in fitness:
            if abs(stat[0]) < minFit:
                minFit = stat[0]
                if abs(stat[1]) > maxFit:
                    maxFit = stat[1]
                avgList.append(stat[2])
    print("\tResults after 30 generations with Tournament Selection, Shuffle Crossover, Random Mutation:")
    print(f'\tFitness \'Best in run\' : {minFit}\n\tFitness \'Worst in Run\' : {maxFit}\n\tFitness Average (overall) : {mean(avgList)}\n')
    listBOR.append(minFit)
    listAVG.append(mean(avgList))
    listSTD.append(getFitness(chrm01)[3])
    #2
    # init first generation
    chrm01 = initChromosome(N, 3)
    chrm01 = selectionSimple(chrm01, 10)    # 1:1 resample
    for gene in range(0, len(chrm01)):
        if(gene % 2 == 0):
            if(random.uniform(1, 10) <= (Pc * 10)):
                tempTuple = crossover(chrm01[gene], chrm01[gene+1])
                chrm01[gene] = mutateRandom(tempTuple[0], Pm)
                chrm01[gene+1] = mutateRandom(tempTuple[1], Pm)
    print(f'\tSimulation 2:')
    printChromosome(chrm01)

    fitness = list()
    avgList = list()
    for i in range(0, maxGenerations):
        "Do 30 runs"    
        selectionTournament(chrm01, k, P)
        for gene in range(0, len(chrm01)):
            if(gene % 2 == 0):
                if(random.uniform(1, 10) <= (Pc * 10)):
                    tempTuple = crossoverUniform(chrm01[gene], chrm01[gene+1])
                    chrm01[gene] = mutateGaussian(tempTuple[0], Pm, sigma)
                    chrm01[gene+1] = mutateGaussian(tempTuple[1], Pm, sigma)
        fitness.append(getFitness(chrm01))
        minFit = 1000
        maxFit = 0
        for stat in fitness:
            if abs(stat[0]) < minFit:
                minFit = stat[0]
                if abs(stat[1]) > maxFit:
                    maxFit = stat[1]
                avgList.append(stat[2])
    print("\tResults after 30 generations with Tournament Selection, 1-Point Crossover, Gaussian Mutation:")
    print(f'\tFitness \'Best in run\' : {minFit}\n\tFitness \'Worst in Run\' : {maxFit}\n\tFitness Average (overall) : {mean(avgList)}\n')
    listBOR.append(minFit)
    listAVG.append(mean(avgList))
    listSTD.append(getFitness(chrm01)[3])

    #3
    # init first generation
    chrm01 = initChromosome(N, 3)
    chrm01 = selectionTournament(chrm01, k, P)    # 1:1 resample
    for gene in range(0, len(chrm01)):
        if(gene % 2 == 0):
            if(random.uniform(1, 10) <= (Pc * 10)):
                tempTuple = crossover(chrm01[gene], chrm01[gene+1])
                chrm01[gene] = mutateRandom(tempTuple[0], Pm)
                chrm01[gene+1] = mutateRandom(tempTuple[1], Pm)
    print(f'\tSimulation 3:')
    printChromosome(chrm01)
    fitness = list()
    avgList = list()
    for i in range(0, maxGenerations):
        "Do 30 runs"    
        selectionProportional(chrm01, N)
        for gene in range(0, len(chrm01)):
            if(gene % 2 == 0):
                if(random.uniform(1, 10) <= (Pc * 10)):
                    tempTuple = crossover(chrm01[gene], chrm01[gene+1])
                    chrm01[gene] = mutateGaussian(tempTuple[0], Pm, sigma)
                    chrm01[gene+1] = mutateGaussian(tempTuple[1], Pm, sigma)
        fitness.append(getFitness(chrm01))
        minFit = 1000
        maxFit = 0
        for stat in fitness:
            if abs(stat[0]) < minFit:
                minFit = stat[0]
                if abs(stat[1]) > maxFit:
                    maxFit = stat[1]
                avgList.append(stat[2])
    print("\tResults after 30 generations with Tournament Selection, 1-Point Crossover, Gaussian Mutation:")
    print(f'\tFitness \'Best in run\' : {minFit}\n\tFitness \'Worst in Run\' : {maxFit}\n\tFitness Average (overall) : {mean(avgList)}\n')
    listBOR.append(minFit)
    listAVG.append(mean(avgList))
    listSTD.append(getFitness(chrm01)[3])
   
    
    #4
    # init first generation
    chrm01 = initChromosome(N, 3)
    chrm01 = selectionSimple(chrm01, 10)    # 1:1 resample
    for gene in range(0, len(chrm01)):
        if(gene % 2 == 0):
            if(random.uniform(1, 10) <= (Pc * 10)):
                tempTuple = crossover(chrm01[gene], chrm01[gene+1])
                chrm01[gene] = mutateRandom(tempTuple[0], Pm)
                chrm01[gene+1] = mutateRandom(tempTuple[1], Pm)
    print(f'\tSimulation 4:')
    printChromosome(chrm01)
    fitness = list()
    avgList = list()
    for i in range(0, maxGenerations):
        "Do 30 runs"    
        selectionProportional(chrm01, N)
        for gene in range(0, len(chrm01)):
            if(gene % 2 == 0):
                if(random.uniform(1, 10) <= (Pc * 10)):
                    tempTuple = crossover(chrm01[gene], chrm01[gene+1])
                    chrm01[gene] = mutateGaussian(tempTuple[0], Pm, sigma)
                    chrm01[gene+1] = mutateGaussian(tempTuple[1], Pm, sigma)
        fitness.append(getFitness(chrm01))
        minFit = float("inf")
        maxFit = 0
        for stat in fitness:
            if abs(stat[0]) < minFit:
                minFit = stat[0]
                if abs(stat[1]) > maxFit:
                    maxFit = stat[1]
                avgList.append(stat[2])
    print("\tResults after 30 generations with Proportional Selection, 1-Point Crossover, Gaussian Mutation:")
    print(f'\tFitness \'Best in run\' : {minFit}\n\tFitness \'Worst in Run\' : {maxFit}\n\tFitness Average (overall) : {mean(avgList)}\n')
    listBOR.append(minFit)
    listAVG.append(mean(avgList))
    listSTD.append(getFitness(chrm01)[3])

    #5
    # init first generation
    chrm01 = initChromosome(N, 3)
    chrm01 = selectionSimple(chrm01, 10)    # 1:1 resample
    for gene in range(0, len(chrm01)):
        if(gene % 2 == 0):
            if(random.uniform(1, 10) <= (Pc * 10)):
                tempTuple = crossover(chrm01[gene], chrm01[gene+1])
                chrm01[gene] = mutateRandom(tempTuple[0], Pm)
                chrm01[gene+1] = mutateRandom(tempTuple[1], Pm)
    print(f'\tSimulation 5:')
    printChromosome(chrm01)
    fitness = list()
    avgList = list()
    for i in range(0, maxGenerations):
        "Do 30 runs"    
        selectionTruncate(chrm01, percentile)
        for gene in range(0, len(chrm01)):
            if(gene % 2 == 0):
                if(random.uniform(1, 10) <= (Pc * 10)):
                    tempTuple = crossoverArithmetic(chrm01[gene], chrm01[gene+1], alpha)
                    chrm01[gene] = mutateGaussian(tempTuple[0], Pm, sigma)
                    chrm01[gene+1] = mutateGaussian(tempTuple[1], Pm, sigma)
        fitness.append(getFitness(chrm01))
        minFit = float("inf")
        maxFit = 0
        for stat in fitness:
            if abs(stat[0]) < minFit:
                minFit = stat[0]
                if abs(stat[1]) > maxFit:
                    maxFit = stat[1]
                avgList.append(stat[2])
    print("\tResults after 30 generations with Truncate Selection, Arithmetic Crossover, Gaussian Mutation:")
    print(f'\tFitness \'Best in run\' : {minFit}\n\tFitness \'Worst in Run\' : {maxFit}\n\tFitness Average (overall) : {mean(avgList)}\n')
    listBOR.append(minFit)
    listAVG.append(mean(avgList))
    listSTD.append(getFitness(chrm01)[3])

    #6
    # init first generation
    chrm01 = initChromosome(N, 3)
    chrm01 = selectionSimple(chrm01, 10)    # 1:1 resample
    for gene in range(0, len(chrm01)):
        if(gene % 2 == 0):
            if(random.uniform(1, 10) <= (Pc * 10)):
                tempTuple = crossover(chrm01[gene], chrm01[gene+1])
                chrm01[gene] = mutateRandom(tempTuple[0], Pm)
                chrm01[gene+1] = mutateRandom(tempTuple[1], Pm)
    print(f'\tSimulation 6:')
    printChromosome(chrm01)
    fitness = list()
    avgList = list()
    for i in range(0, maxGenerations):
        "Do 30 runs"    
        selectionRankingFitness(chrm01)
        for gene in range(0, len(chrm01)):
            if(gene % 2 == 0):
                if(random.uniform(1, 10) <= (Pc * 10)):
                    tempTuple = crossoverShuffle(chrm01[gene], chrm01[gene+1])
                    chrm01[gene] = mutateRandom(tempTuple[0], .10)
                    chrm01[gene+1] = mutateRandom(tempTuple[1], .10)
        fitness.append(getFitness(chrm01))
        minFit = 1000
        maxFit = 0
        for stat in fitness:
            if abs(stat[0]) < minFit:
                minFit = stat[0]
                if abs(stat[1]) > maxFit:
                    maxFit = stat[1]
                avgList.append(stat[2])
    print("\tResults after 30 generations with Ranking Selection, Shuffle Crossover, Random Mutation:")
    print(f'\tFitness \'Best in run\' : {minFit}\n\tFitness \'Worst in Run\' : {maxFit}\n\tFitness Average (overall) : {mean(avgList)}\n')
    listBOR.append(minFit)
    listAVG.append(mean(avgList))
    listSTD.append(getFitness(chrm01)[3])

    #7
    # init first generation
    chrm01 = initChromosome(N, 3)
    chrm01 = selectionSimple(chrm01, 10)    # 1:1 resample
    for gene in range(0, len(chrm01)):
        if(gene % 2 == 0):
            if(random.uniform(1, 10) <= (Pc * 10)):
                tempTuple = crossover(chrm01[gene], chrm01[gene+1])
                chrm01[gene] = mutateRandom(tempTuple[0], Pm)
                chrm01[gene+1] = mutateRandom(tempTuple[1], Pm)
    print(f'\tSimulation 7:')
    printChromosome(chrm01)
    fitness = list()
    avgList = list()
    for i in range(0, maxGenerations):
        "Do 30 runs"    
        selectionRankingFitness(chrm01)
        for gene in range(0, len(chrm01)):
            if(gene % 2 == 0):
                if(random.uniform(1, 10) <= (Pc * 10)):
                    tempTuple = crossover(chrm01[gene], chrm01[gene+1])
                    chrm01[gene] = mutateGaussian(tempTuple[0], .10, sigma)
                    chrm01[gene+1] = mutateGaussian(tempTuple[1], .10, sigma)
        fitness.append(getFitness(chrm01))
        minFit = float("inf")
        maxFit = 0
        for stat in fitness:
            if abs(stat[0]) < minFit:
                minFit = stat[0]
                if abs(stat[1]) > maxFit:
                    maxFit = stat[1]
                avgList.append(stat[2])
    print("\tResults after 30 generations with Ranking Selection, 1-Point Crossover, Gaussian Mutation:")
    print(f'\tFitness \'Best in run\' : {minFit}\n\tFitness \'Worst in Run\' : {maxFit}\n\tFitness Average (overall) : {mean(avgList)}\n')
    listBOR.append(minFit)
    listAVG.append(mean(avgList))
    listSTD.append(getFitness(chrm01)[3])

    #8
    # init first generation
    chrm01 = initChromosome(N, 3)
    chrm01 = selectionSimple(chrm01, 10)    # 1:1 resample
    for gene in range(0, len(chrm01)):
        if(gene % 2 == 0):
            if(random.uniform(1, 10) <= (Pc * 10)):
                tempTuple = crossover(chrm01[gene], chrm01[gene+1])
                chrm01[gene] = mutateRandom(tempTuple[0], Pm)
                chrm01[gene+1] = mutateRandom(tempTuple[1], Pm)
    print(f'\tSimulation 8:')
    printChromosome(chrm01)
    fitness = list()
    avgList = list()
    for i in range(0, maxGenerations):
        "Do 30 runs"    
        selectionTournament(chrm01, 3, P)
        for gene in range(0, len(chrm01)):
            if(gene % 2 == 0):
                if(random.uniform(1, 10) <= (Pc * 10)):
                    tempTuple = crossoverShuffle(chrm01[gene], chrm01[gene+1])
                    chrm01[gene] = mutateGaussian(tempTuple[0], .10, sigma)
                    chrm01[gene+1] = mutateGaussian(tempTuple[1], .10, sigma)
        fitness.append(getFitness(chrm01))
        minFit = 1000
        maxFit = 0
        for stat in fitness:
            if abs(stat[0]) < minFit:
                minFit = stat[0]
                if abs(stat[1]) > maxFit:
                    maxFit = stat[1]
                avgList.append(stat[2])
    print("\tResults after 30 generations with Tournament(3) Selection, Shuffle Crossover, Gaussian Mutation:")
    print(f'\tFitness \'Best in run\' : {minFit}\n\tFitness \'Worst in Run\' : {maxFit}\n\tFitness Average (overall) : {mean(avgList)}\n')
    listBOR.append(minFit)
    listAVG.append(mean(avgList))
    listSTD.append(getFitness(chrm01)[3])

     #9
    # init first generation
    chrm01 = initChromosome(N, 3)
    chrm01 = selectionSimple(chrm01, 10)    # 1:1 resample
    for gene in range(0, len(chrm01)):
        if(gene % 2 == 0):
            if(random.uniform(1, 10) <= (Pc * 10)):
                tempTuple = crossover(chrm01[gene], chrm01[gene+1])
                chrm01[gene] = mutateRandom(tempTuple[0], Pm)
                chrm01[gene+1] = mutateRandom(tempTuple[1], Pm)
    print(f'\tSimulation 9:')
    printChromosome(chrm01)
    fitness = list()
    avgList = list()
    for i in range(0, maxGenerations):
        "Do 30 runs"    
        selectionTournament(chrm01, 5, P)
        for gene in range(0, len(chrm01)):
            if(gene % 2 == 0):
                if(random.uniform(1, 10) <= (Pc * 10)):
                    tempTuple = crossoverArithmetic(chrm01[gene], chrm01[gene+1], alpha)
                    chrm01[gene] = mutateGaussian(tempTuple[0], .10, sigma)
                    chrm01[gene+1] = mutateGaussian(tempTuple[1], .10, sigma)
        fitness.append(getFitness(chrm01))
        minFit = float("inf")
        maxFit = 0
        for stat in fitness:
            if abs(stat[0]) < minFit:
                minFit = stat[0]
                if abs(stat[1]) > maxFit:
                    maxFit = stat[1]
                avgList.append(stat[2])
    print("\tResults after 30 generations with Tournament(5) Selection, Arithmetic Crossover, Gaussian Mutation:")
    print(f'\tFitness \'Best in run\' : {minFit}\n\tFitness \'Worst in Run\' : {maxFit}\n\tFitness Average (overall) : {mean(avgList)}\n')
    listBOR.append(minFit)
    listAVG.append(mean(avgList))
    listSTD.append(getFitness(chrm01)[3])

     #10
    # init first generation
    chrm01 = initChromosome(N, 3)
    chrm01 = selectionSimple(chrm01, 10)    # 1:1 resample
    for gene in range(0, len(chrm01)):
        if(gene % 2 == 0):
            if(random.uniform(1, 10) <= (Pc * 10)):
                tempTuple = crossover(chrm01[gene], chrm01[gene+1])
                chrm01[gene] = mutateRandom(tempTuple[0], .5)
                chrm01[gene+1] = mutateRandom(tempTuple[1], .5)
    print(f'\tSimulation 10:')
    printChromosome(chrm01)
    fitness = list()
    avgList = list()
    for i in range(0, maxGenerations):
        "Do 30 runs"    
        selectionTruncate(chrm01, percentile)
        for gene in range(0, len(chrm01)):
            if(gene % 2 == 0):
                if(random.uniform(1, 10) <= (Pc * 10)):
                    tempTuple = crossoverUniform(chrm01[gene], chrm01[gene+1])
                    chrm01[gene] = mutateGaussian(tempTuple[0], Pm, .3)
                    chrm01[gene+1] = mutateGaussian(tempTuple[1], Pm, .3)
        fitness.append(getFitness(chrm01))
        minFit = float("inf")
        maxFit = 0
        for stat in fitness:
            if abs(stat[0]) < minFit:
                minFit = stat[0]
                if abs(stat[1]) > maxFit:
                    maxFit = stat[1]
                avgList.append(stat[2])
    print("\tResults after 30 generations with Truncate Selection, Uniform Crossover, Gaussian Mutation:")
    print(f'\tFitness \'Best in run\' : {minFit}\n\tFitness \'Worst in Run\' : {maxFit}\n\tFitness Average (overall) : {mean(avgList)}\n')
    listBOR.append(minFit)
    listAVG.append(mean(avgList))
    listSTD.append(getFitness(chrm01)[3])

    print("---> After 10 simulations of different configurations (each sim with 30 generations)")
    print(f'\tBest of all runs: {min(listBOR, key=abs)} | stdev: {stdev(listBOR)}\n\tAverage of all best of runs: {mean(listBOR)} | stdev: {stdev(listBOR)}\n\tAverage of all run averages: {mean(listAVG)} | stdev: {stdev(listAVG)}\n')
    
    t = PrettyTable(['Simulation #', 'Best of 30 Runs', 'Average of 30 BORs', 'Standard Deviation of 30 BORs'])
    for i in range(0,10):
        t.add_row([f'{i+1}', f'{listBOR[i]}', f'{listAVG[i]}', f'{listSTD[i]}'])
    print(t)


        #for runs in range (1, 31):
    #    if(runs % 10 == 0):
    #        print(f'-----Run #{runs} ----->')
    #        fitness = []
    #        avgList = []

    #    chrm01 = initChromosome(N, 3)
    #    for i in range(1,maxGenerations+1): # should just be an init case
    #        chrm01 = random.sample(chrm01, 10) # this is 1:1 sampling, fine for init
    #        for index in range(0,len(chrm01)):
    #            if(index % 2 == 0):
    #                if(random.uniform(1,10) <= (Pc * 10)):
    #                    tempTuple = crossover(chrm01[index], chrm01[index+1])
    #                    chrm01[index] = mutate(tempTuple[0], Pm, sigma)
    #                    chrm01[index+1] = mutate(tempTuple[1], Pm, sigma)
    #        if i %  10 == 0 and runs % 10 == 0:
    #            print(f'\tNew generation: {i}') 
    #            printChromosome(chrm01)
    #            fitness.append(getFitness(chrm01))
 
    #    if(runs % 10 == 0):
    #        minFit = fitness[0][0]
    #        maxFit = fitness[0][1]
    #        for stat in fitness:
    #            if abs(stat[0]) < minFit:
    #                minFit = stat[0]
    #            if abs(stat[1]) > maxFit:
    #                maxFit = stat[1]
    #            avgList.append(stat[2])
    #        print(f'\tFitness \'Best in run\' : {minFit}\n\tFitness \'Worst in Run\' : {maxFit}\n\tFitness Average (overall) : {mean(avgList)}\n')
    #        listBOR.append(minFit)
    #        listAVG.append(mean(avgList))
    #print(f'Best of all runs: {min(listBOR, key=abs)} | stdev: {stdev(listBOR)}\nAverage of all best of runs: {mean(listBOR)} | stdev: {stdev(listBOR)}\nAverage of all run averages: {mean(listAVG)} | stdev: {stdev(listAVG)}\n')










    # init first generation
    #chrm01 = initChromosome(N, 3)
    #chrm01 = selectionSimple(chrm01, 10)    # 1:1 resample
    #for gene in range(0, len(chrm01)):
    #    if(gene % 2 == 0):
    #        if(random.uniform(1, 10) <= (Pc * 10)):
    #            tempTuple = crossover(chrm01[gene], chrm01[gene+1])
    #            chrm01[gene] = mutateRandom(tempTuple[0], sigma)
    #            chrm01[gene+1] = mutateRandom(tempTuple[1], sigma)
    #print(f'\tInitial Chromosome:')
    #printChromosome(chrm01)

    #nextGen = selectionRankingFitness(chrm01)
    #for gene in range(0, len(nextGen)):
    #    if(gene % 2 == 0):
    #        if(random.uniform(1, 10) <= (Pc * 10)):
    #            tempTuple = crossover(nextGen[gene], nextGen[gene+1])
    #            nextGen[gene] = mutateGaussian(tempTuple[0], Pm, sigma)
    #            nextGen[gene+1] = mutateGaussian(tempTuple[1], Pm, sigma)
    #print(f'\tNext Generation')
    #printChromosome(nextGen)

    #nextGen2 = selectionSimple(nextGen, 10)
    #for gene in range(0, len(nextGen2)):
    #    if(gene % 2 == 0):
    #        if(random.uniform(1, 10) <= (Pc * 1000)):
    #            tempTuple = crossoverShuffle(nextGen[gene], nextGen2[gene+1])
    #            nextGen2[gene] = mutateGaussian(tempTuple[0], Pm, sigma)
    #            nextGen2[gene+1] = mutateGaussian(tempTuple[1], Pm, sigma)
    #print(f'\tNext Generation')
    #printChromosome(nextGen2)

    #nextGen3 = selectionTournament(nextGen2, k, P)
    #printChromosome(nextGen3)
    #nextGen4 = selectionProportional(nextGen3, len(nextGen3))
    #printChromosome(nextGen4)




    #project 01
    #    listBOR = []
    #    listAVG = []
    #
    #""" do 30 runs """
    #for runs in range (1, 31):
    #    if(runs % 10 == 0):
    #        print(f'-----Run #{runs} ----->')
    #        fitness = []
    #        avgList = []

    #    chrm01 = initChromosome(N, 3)
    #    for i in range(1,maxGenerations+1): # should just be an init case
    #        chrm01 = random.sample(chrm01, 10) # this is 1:1 sampling, fine for init
    #        for index in range(0,len(chrm01)):
    #            if(index % 2 == 0):
    #                if(random.uniform(1,10) <= (Pc * 10)):
    #                    tempTuple = crossover(chrm01[index], chrm01[index+1])
    #                    chrm01[index] = mutate(tempTuple[0], Pm, sigma)
    #                    chrm01[index+1] = mutate(tempTuple[1], Pm, sigma)
    #        if i %  10 == 0 and runs % 10 == 0:
    #            print(f'\tNew generation: {i}') 
    #            printChromosome(chrm01)
    #            fitness.append(getFitness(chrm01))
 
    #    if(runs % 10 == 0):
    #        minFit = fitness[0][0]
    #        maxFit = fitness[0][1]
    #        for stat in fitness:
    #            if abs(stat[0]) < minFit:
    #                minFit = stat[0]
    #            if abs(stat[1]) > maxFit:
    #                maxFit = stat[1]
    #            avgList.append(stat[2])
    #        print(f'\tFitness \'Best in run\' : {minFit}\n\tFitness \'Worst in Run\' : {maxFit}\n\tFitness Average (overall) : {mean(avgList)}\n')
    #        listBOR.append(minFit)
    #        listAVG.append(mean(avgList))
    #print(f'Best of all runs: {min(listBOR, key=abs)} | stdev: {stdev(listBOR)}\nAverage of all best of runs: {mean(listBOR)} | stdev: {stdev(listBOR)}\nAverage of all run averages: {mean(listAVG)} | stdev: {stdev(listAVG)}\n')

if __name__ == "__main__":
    main() 