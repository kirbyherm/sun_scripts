#!/mnt/simulations/secarml/soft/anaconda3/bin/python
#!/mnt/misc/sw/x86_64/all/anaconda/python3.7/bin/python

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
import matplotlib
import time
import itertools
import timeit
#from cosy import cosyrun, write_fox
import pygmo as pg
#from problem import optimizeRes
import pandas as pd

from utils import load_configs, get_best_point

from scipy.stats.distributions import chi2

# specify Tex details for pretty plots
#os.environ['PATH'] = os.environ['PATH'] + ':/mnt/misc/sw/indep/all/texlive/2013/bin/x86_64-linux/latex'
##os.environ['PATH'] = os.environ['PATH'] + ':/usr/bin/tex'
#plt.rcParams.update({
#    "text.usetex": True,
#})

script, batch, generations = sys.argv
optimized_params = 3
filename = './fit_{}d_db_{}s_{}.h5'.format(optimized_params, batch, generations)
fNom = np.zeros(optimized_params)+1
fNames = ["TAS","SumSeg","Mult","sum1"]
fNames = fNames[:optimized_params]
#print(len(fNames))
config = load_configs()
energies = config['levels']
magnet_dim = len(energies)
names = []
for h in energies:
	names.append('h'+str(h))

best_point = get_best_point("tas")

# read pop from h5 file (i.e. after running view_db.py)
def read_pop_df(filename, pop=None):
    df = pd.read_hdf(filename)
    return df

def plot_hists(df, df_reduce, filename):

    write_qnames = {'q1':'q1','q2':'q2','q3':'q3','q4':'q4','q5':'q5','q6':'q6','q7':'q7','q8':'q8','q9':'q9','q10':'q10','q11':'q11','q12':'q12','q13':'q13','q14':'q14','q15':'q15','q16':'h1','q17':'h2','q18':'h3','q19':'o1'}
    df = df.iloc[:,:magnet_dim]
#    df_clos = df_reduce.loc[df_reduce['closest']==True].sort_values(by=['mult']).reset_index(drop=True)
    df_best = df_reduce.copy()
    df_reduce = df_reduce.iloc[:,:magnet_dim]
#    df = df.rename(columns=write_qnames)
#    fig = plt.figure(1)
    print(magnet_dim, math.ceil(magnet_dim / 8))
    fig, axs = plt.subplots(math.ceil(magnet_dim / 8),8)
#    print(axs)
    if magnet_dim % 8 > 0:
        for i in range(magnet_dim % 8,8):
            fig.delaxes(axs[math.floor(magnet_dim/8)][i])
    axs = (fig.axes)
    bins = np.linspace(0,0.5,100)
    hists = df.hist(bins=bins,log=True,ax=axs,grid=False)
#    print(hists[0].bins)
#    fig.yscale('log')
#    fig.tight_layout()
    fig.set_figheight(10)
    fig.set_figwidth(19)
#    ax = plt.gca()
#    print(ax.get_ylim())
    axs2 = axs
    for i in range(len(axs)):
#    for j in range(len(axs[i])):
#        print(axs[i].get_ylim())
        axs[i].axvline(x=best_point[i], color='red')
#        axs[i].set_ylim([100,1000000])
        axs[i].set_ylabel('all points')
        if i not in [2,3]:
            axs[i].set_xlim([-0.02,0.2])
        else:
            axs[i].set_xlim([0,0.5])
#        axs2[i] = axs[i].twinx()
#        axs2[i].set_ylim([1,100000])
#        axs2[i].set_ylabel('optimized points',rotation=-90)
#    number_of_clusters = np.max(df_best['kcluster']+1)
    colors = list(plt.get_cmap('tab20').colors)
    write_qnames = names[:magnet_dim]
#    print(df_best,axs2)
#    df_clos['yplot'] = 0
#    for i in range(number_of_clusters):
##                ax = df.loc[(df['kcluster']==i)].plot(x='FP4_BeamSpot',y=obj,style='o',color=colors[i],label=df_clos.loc[df_clos['kcluster']==i].index[0]+1,markersize=3.0)
#        hist = df_best.loc[df_best['kcluster']==i][write_qnames].hist(bins=bins,log=True,ax=axs2,color=np.array([colors[i+1]]),label=df_clos.loc[(df_clos['kcluster']==i)].index[0]+1,grid=False)
#        for j in range(len(axs)):
#            y_min, y_max = axs2[j].get_ylim()
#            df_clos['yplot'][i] = df_clos.index*(y_max-y_min)+y_min
#            print((df_clos['kcluster'][i])*(np.logspace(1.0,3.0,num=10)[i]),(df_clos['kcluster'][i]),(np.logspace(0.0,2.5,num=10)[i]))
#            axs2[j].text(df_clos["q{0}".format(j+1)][i],(np.logspace(0.0,2.8,num=10)[df_clos['kcluster'][i]]),str(np.max(df_clos.loc[(df_clos['kcluster']==i)]['kcluster'])+1),color='black')
#    for i in range(len(axs)):
#        axs2[i].set_title("")

#    print(hist)
    plt.savefig(filename+'full_hist.png')
    plt.cla()
    
    return

def plot_chi2_func(df, metric, filename):

#    df = df.iloc[:,:magnet_dim]
#    df_clos = df_reduce.loc[df_reduce['closest']==True].sort_values(by=['mult']).reset_index(drop=True)
#    df = df.rename(columns=write_qnames)
#    fig = plt.figure(1)
#    print(axs)
    fig, axs = plt.subplots(math.ceil(magnet_dim / 8),8)
#    print(axs)
    if magnet_dim % 8 > 0:
        for i in range(magnet_dim % 8,8):
            fig.delaxes(axs[math.floor(magnet_dim/8)][i])
    axs = (fig.axes)
#    bins = np.linspace(0,0.5,100)
#    hists = df.hist(bins=bins,log=True,ax=axs,grid=False)
#    print(hists[0].bins)
#    fig.yscale('log')
    fig.set_figheight(10)
    fig.set_figwidth(19)
#    ax = plt.gca()
#    print(ax.get_ylim())
    axs2 = axs
    metric_chi2 = df[metric].min()
    dof = 51
    if metric == 'mult':
        dof = 1
    elif metric == 'total':
        dof = 148+87
    df_cut = df.loc[df[metric] < metric_chi2 + chi2.ppf(0.683,dof)]
    continuum = 0
    continuum_m = 0
    continuum_p = 0
    for i in range(len(axs)):
#    for j in range(len(axs[i])):
#        print(axs[i].get_ylim())
        axs[i].set_title(df.columns[i])
        axs[i].plot(df.iloc[:,i], df[metric],'.',linestyle='none' )
        axs[i].axvline(x=best_point[i], color='red')
        axs[i].axvline(x=df_cut.iloc[:,i].min(), color='red', linestyle='dashed')
        axs[i].axvline(x=df_cut.iloc[:,i].max(), color='red', linestyle='dashed')
        axs[i].axhline(y=metric_chi2+chi2.ppf(0.683,dof), color='black')
#        axs[i].set_ylim([0,(metric_chi2+chi2.ppf(0.683,dof))*2])
        if i % 8 == 0:
            axs[i].set_ylabel('chi2')
        if names[i] not in ['h2000','h1750','h2475']:
            axs[i].set_xlim([-0.01,0.2])
        else:
            axs[i].set_xlim([0,0.4])
        if i < 7:
            print("{0}, {1:.2f}, {2:.2f}, {3:.2f}".format(names[i], best_point[i]*100, abs(df_cut.iloc[:,i].min()-best_point[i])*100, abs(df_cut.iloc[:,i].max()-best_point[i])*100))
        else:
            continuum += best_point[i]
            continuum_m += abs(df_cut.iloc[:,i].min()-best_point[i])**(2)
            continuum_p += abs(df_cut.iloc[:,i].max()-best_point[i])**(2)
#        axs[i].set_yscale('log')
        axs[i].set_ylim([metric_chi2-51,metric_chi2+chi2.ppf(0.683,dof)+100])
#        axs2[i] = axs[i].twinx()
#        axs2[i].set_ylabel('optimized points',rotation=-90)
#    number_of_clusters = np.max(df_best['kcluster']+1)
    print("{0}, {1:.2f}, {2:.2f}, {3:.2f}".format("continuum", continuum*100, (continuum_m)**(0.5)*100, (continuum_p)**(0.5)*100))
    plt.tight_layout()
    colors = list(plt.get_cmap('tab20').colors)
    write_qnames = names[:magnet_dim]
#    print(df_best,axs2)
#    df_clos['yplot'] = 0
#    for i in range(number_of_clusters):
##                ax = df.loc[(df['kcluster']==i)].plot(x='FP4_BeamSpot',y=obj,style='o',color=colors[i],label=df_clos.loc[df_clos['kcluster']==i].index[0]+1,markersize=3.0)
#        hist = df_best.loc[df_best['kcluster']==i][write_qnames].hist(bins=bins,log=True,ax=axs2,color=np.array([colors[i+1]]),label=df_clos.loc[(df_clos['kcluster']==i)].index[0]+1,grid=False)
#        for j in range(len(axs)):
#            y_min, y_max = axs2[j].get_ylim()
#            df_clos['yplot'][i] = df_clos.index*(y_max-y_min)+y_min
#            print((df_clos['kcluster'][i])*(np.logspace(1.0,3.0,num=10)[i]),(df_clos['kcluster'][i]),(np.logspace(0.0,2.5,num=10)[i]))
#            axs2[j].text(df_clos["q{0}".format(j+1)][i],(np.logspace(0.0,2.8,num=10)[df_clos['kcluster'][i]]),str(np.max(df_clos.loc[(df_clos['kcluster']==i)]['kcluster'])+1),color='black')
#    for i in range(len(axs)):
#        axs2[i].set_title("")

    print(np.min(df[metric]))
    plt.savefig(filename+metric+'chi2func.png')
    plt.cla()
    
    return

def plot_chi2_hists(df, df_reduce, filename):

    df = df.iloc[:,magnet_dim:-2]/np.array([87-magnet_dim,87-magnet_dim,70,87+87+10-magnet_dim])
#    df_clos = df_reduce.loc[df_reduce['closest']==True].sort_values(by=['mult']).reset_index(drop=True)
    df_best = df_reduce.copy()
    df_reduce = df_reduce.iloc[:,magnet_dim:-2]
#    df = df.rename(columns=write_qnames)
#    fig = plt.figure(1)
    fig, axs = plt.subplots(2,2)
#    print(axs)
#    fig.delaxes(axs[5][4])
#    fig.delaxes(axs[3][5])
#    fig.delaxes(axs[5][2])
    axs = (fig.axes)
    bins = np.linspace(1.0,5,100)
    hists = df.hist(bins=bins,log=True,ax=axs,grid=False)
#    print(hists[0].bins)
#    fig.yscale('log')
#    fig.tight_layout()
    fig.set_figheight(10)
    fig.set_figwidth(19)
#    ax = plt.gca()
#    print(ax.get_ylim())
    axs2 = axs
#    for i in range(len(axs)):
##    for j in range(len(axs[i])):
##        print(axs[i].get_ylim())
#        axs[i].axvline(x=best_point[i], color='red')
##        axs[i].set_ylim([100,1000000])
#        axs[i].set_ylabel('all points')
#        if i not in [2,3]:
#            axs[i].set_xlim([-0.02,0.2])
#        else:
#            axs[i].set_xlim([0,0.5])
#        axs2[i] = axs[i].twinx()
#        axs2[i].set_ylim([1,100000])
#        axs2[i].set_ylabel('optimized points',rotation=-90)
#    number_of_clusters = np.max(df_best['kcluster']+1)
    colors = list(plt.get_cmap('tab20').colors)
    plt.savefig(filename+'chi2_hist.png')
    plt.cla()
    
    return


def main(filename,batch):

    file_extension = os.path.splitext(filename)[-1]
    print(os.path.split(filename))
    popi = None
    print("reading {} file".format(file_extension))
    df = None
    if file_extension == ".h5":
        df = read_pop_df(filename)
    df_full = pd.read_hdf(filename)
    df = df.loc[df['branch_sum'] < 1.5e3]
    plot_hists(df, df_full, filename)
#    plot_chi2_hists(df, df_full, filename)
#    objs = ['tas','ss','mult','tas_m1','total']
    objs = ['tas']
    for i in range(len(objs)):
        if i > 0 and i!=1:
            continue 
        plot_chi2_func(df, objs[i], filename)

if __name__=='__main__':
    main(filename,batch)



