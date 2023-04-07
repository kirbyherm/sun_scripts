#!/mnt/simulations/secarml/soft/anaconda3/bin/python

# make sure above path points to the version of python where you have pygmo installed 
# nscl servers
#!/mnt/misc/sw/x86_64/all/anaconda/python3.7/bin/python
# hpcc servers
#!/mnt/home/herman67/anaconda3/envs/pygmo/bin/python

#import commands
import sys, math
import os, shutil, signal
import subprocess as commands
import re
import random
import numpy as np
import matplotlib.pyplot as plt
import time
import itertools
import timeit

import pygmo as pg
from problem import optimizeRes
import pandas as pd
from utils import load_all, multi_minimize, chisq, run, load_configs
from pandas import read_csv
import json


config = load_configs()
energies = config['levels']
WORK_DIR = './'
SCRATCH_DIR = '/scratch/hermanse/{}'.format(config['isotope'])
names = []
for h in energies:
	names.append('h'+str(h))

# take batch_no as input from command line (this specifies outputFile)
script, batch_no = sys.argv

# set a specific seed for the algorithm
#   note that each island uses a random init pop
#   thus none of the islands are identical
#   (unless init pop is identical)
seed = 56448189

# MOEAD hyperparameters
#   default parameters have worked well
generations = 20000
cr_p = 1.0 # crossover parameter, 1.0 by default
f_p = 0.5 # diff evolution operator parameter, 0.5 by default
eta_m = 20 # distribution index used by the polynomial mutation, 20 by default
realb = 0.9 # chance that the neighborhood is considered at each generation, rather than the whole population, 0.9 by default
neighbors = 50 # size of the weight's neighborhood, 20 by default
limit = 2 # max number of copies reinserted in the population if preserve_diversity=True, 2 by default
preserve_diversity=True # activates diversity preservation mechanisms

# specify number of magnets to tune
magnet_dim = len(names)
# specify output file
fileName = "output_4f_moead_FP2_FP3_{}_{}.csv".format(generations, batch_no)
outputFile = "{}{}".format(SCRATCH_DIR, fileName)
n_obj = 2

objs = ['tas','tas_m1','branch_sum']
objs = ['tas','branch_sum']
columns = names
for i in range(len(objs)):
    columns.append(objs[i])

# write population to h5 file
def save_pop(pop, db_out):

    xs = pop.get_x()
    fs = pop.get_f()
    df = pd.DataFrame(np.hstack((xs, fs)), columns=columns)            

    # Check whether the specified path exists or not
    isExist = os.path.exists(db_out)
    df_orig = None    
    if isExist:
        df_orig = pd.read_hdf(db_out)      
        df = pd.concat([df_orig, df], axis=0)
    df.to_hdf(db_out,key='df')
    
    return    

# main function
def main(pop_init=None):

    # start timer
    startTime = timeit.default_timer()

    # there should be a way to run batch_fitness evaluations, haven't gotten it to work on NSCL
    #b = pg.bfe()
    #alg.set_bfe(b)
    #alg.set_verbosity(1)
    
    # initialize problem
    df_tas, df_tas_sim, df_ss, df_ss_sim, df_mult, df_mult_sim, df_ss_m1, df_ss_m1_sim = load_all()

    p_optimizeRes = pg.problem(optimizeRes(magnet_dim,df_tas, df_tas_sim, df_ss, df_ss_sim, df_mult, df_mult_sim, df_ss_m1, df_ss_m1_sim, outputFile))

    # relic from old evolutions, which launched all islands internally
    #   this can allow for interconnectivity between islands
    #   now each run of this script calls its own island
    n_islands = 1 
    # if more than one island and want to exchange info need to define a topology to connect them    
    #top = pg.topology(pg.fully_connected(n_islands,1.0))

    # when running 5 objectives, pop needed to be 70    
    pop_n = 969
    if p_optimizeRes.get_nobj() == 3:
        # with 4 objs need pop=84
        pop_n = 990 
        pop_n = 1770 

    # check if we are using an input population or need to initialize one
    pop_new = None
    if (pop_init==None):    
        # randomly create a new population of size pop_n
        pop_new = pg.population(p_optimizeRes, size=pop_n) 
        print("initialize pop")
    else:
        pop_new = pop_init
        print("provided pop")

    print(p_optimizeRes, pop_n)

    save_freq = 200
    db_out = SCRATCH_DIR + "fit_{}d_db_{}s_{}.h5".format(n_obj, batch_no, generations)
    perc_prog = 0
    i = 0
    while i < int(generations):
        if perc_prog > 50:
            save_freq = 100
        if perc_prog > 70:
            save_freq = 20
        # initialize algorithm with hyperparameters
        alg_constructor = pg.moead(gen=save_freq,neighbours=neighbors,seed=seed)
        alg = pg.algorithm(alg_constructor)
        #alg.set_verbosity(1)
    
        pop_new = alg.evolve(pop_new)
        save_pop(pop_new, db_out)
        if int(i / generations *100 ) > perc_prog:
#            print ( i * save_freq )
            perc_prog = int(i/generations * 100)
        i += save_freq
    # check total time    
    print ('Running time (sec): %f' % (timeit.default_timer() - startTime))

    # when I tried to use multiprocessing I needed this
    #pg.mp_island.shutdown_pool()
    #pg.mp_bfe.shutdown_pool()
    pop_final = pop_new

    return pop_final

# if want to initialize a pop with a single nominal gene and random others
#   the evolution will quickly "forget" the nominal though as it prefers to find new solutions
def init_pop(dim,pop_size):
    qNom = np.array([0.        , 0.        , 0.00242834, 0.24587247, 0.26753184,0.03340626, 0.02989332, 0.16553433, 0.01720398, 0.00967029,0.        , 0.00282008, 0.0465006 , 0.04302853, 0.        ,0.01044767, 0.01423525, 0.02762151, 0.        , 0.01414714,0.01830467, 0.01514356, 0.              ])
    df_tas = read_csv("ti57_tas_obs.csv",names=["bin","content","error"],dtype=(float,float))
    df_ss = read_csv("ti57_ss_obs.csv",names=["bin","content","error"],dtype=(float,float))
    df_mult = read_csv("ti57_mult_obs.csv",names=["bin","content","error"],dtype=(float,float))
    scale_factor = 1.0
    scale_factor = np.sum(df_tas['content'])
    df_ss_sim = load_ss_sims(scale_factor)
    df_tas_sim = load_tas_sims(scale_factor,int(10000/df_tas.shape[0]))
    scale_factor = np.sum(df_ss['content'])
    df_ss_sim = load_ss_sims(scale_factor)
    scale_factor = np.sum(df_mult['content'])
    df_mult_sim = load_mult_sims(scale_factor)
#    print(np.sum(qNom),total_minimize(qNom,df_tas, df_tas_sim, df_ss, df_ss_sim, df_mult, df_mult_sim))
    #print(df_tas_sim, df_ss_sim, df_mult_sim)
    p_optimizeRes = pg.problem(optimizeRes(magnet_dim,df_tas, df_tas_sim, df_ss, df_ss_sim, df_mult, df_mult_sim, outputFile))
    pop = pg.population(prob=optimizeRes(dim,df_tas, df_tas_sim, df_ss, df_ss_sim, df_mult, df_mult_sim,out=outputFile))
    pop.push_back(x=qNom,f=np.zeros(p_optimizeRes.get_nobj()))#[0:dim])
#    print(pop)
#    print(pop.problem.get_fevals())
    pop.set_x(0,qNom)
#    print(pop)
    for i in range(pop_size-1):
        pop.push_back(x=qNom,f=np.zeros(p_optimizeRes.get_nobj()))#[0:dim])
    return pop


# if want to read in a specific population (as from a previous run)
#   read in from csv file, with header x0-x10, f0-f3
#     for a 11 magnet, 4 objective case
def read_pop(filename):
    df = pd.read_csv(filename)
    p_optimizeRes = pg.problem(optimizeRes(magnet_dim,out=outputFile))
    pop = pg.population(p_optimizeRes)
    nrow, ncol = df.shape
    print(df)
    for i in range(nrow):
        xs = []
        for j in range(magnet_dim):
            xs.append(df["x"+str(j)][i]) 
        xs = np.asarray(xs)
        fs = []
        for j in range(ncol-magnet_dim):
            fs.append(df["f"+str(j)][i])
        pop.push_back(xs,f=fs)
    return pop    


# when called from console (as the batch script will)
if __name__=='__main__':


    # Check whether the specified path exists or not
    isExist = os.path.exists(SCRATCH_DIR)
    
    if not isExist:
      
      # Create a new directory because it does not exist 
      os.makedirs(SCRATCH_DIR)
    # if using a preexisting or want to pass an identical population
    #popi = init_pop(magnet_dim,990)
    #popi = read_pop("init_pop.csv")

    # else just initialize and evolve random pop
    #print(run())
    pop2 = main()
    db_scratch = SCRATCH_DIR + "fit_{}d_db_{}s_{}.h5".format(n_obj, batch_no, generations)
    db_out = WORK_DIR + "fit_{}d_db_{}s_{}.h5".format(n_obj, batch_no, generations)
#    save_pop(pop2, db_out)
    shutil.move(db_scratch, db_out)
#    os.rmdir(SCRATCH_DIR)
