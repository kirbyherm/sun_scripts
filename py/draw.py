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
import pygmo as pg
import pandas as pd
from utils import run
from scipy.stats import chi2

from utils import load_configs
# specify Tex details for pretty plots
#os.environ['PATH'] = os.environ['PATH'] + ':/mnt/misc/sw/indep/all/texlive/2013/bin/x86_64-linux/latex'
##os.environ['PATH'] = os.environ['PATH'] + ':/usr/bin/tex'
#plt.rcParams.update({
#    "text.usetex": True,
#})
config = load_configs()
energies = config['levels']
magnet_dim = len(energies)
names = []
for h in energies:
	names.append('h'+str(h))

inputs = sys.argv
optimized_params = 5
fNom = np.zeros(optimized_params)+1
fNames = [r"{FP1-res}${}^{-1}$",r"{FP2-res}${}^{-1}$",r"{FP3-res}${}^{-1}$",r"MaxBeamWidth",r"BeamSpotSize"]
fNames = fNames[:optimized_params]
#print(len(fNames))
#magnet_dim = 41
n_sims = int(inputs[1])
dof = 137-magnet_dim

# read pop from h5 file (i.e. after running view_db.py)
def read_pop_df(filename, pop=None):
    df = pd.read_hdf(filename)
    return df

def main(n_sims):

    filename = 'correct_chi2mc_{}.h5'.format(n_sims)
    file_extension = os.path.splitext(filename)[-1]
    print(os.path.split(filename))
    print("reading {} file".format(file_extension))
    df = None
    if file_extension == ".h5":
        df = read_pop_df(filename)
#    df = df.sort_values(by='ssobjs',ignore_index=True,ascending=False)
    fig, axs = plt.subplots(6,8)
    for i in range(3,8):
        fig.delaxes(axs[5][i])
    fig.set_figheight(10)
    fig.set_figwidth(19)
#    fig.delaxes(ax[3][5])
    ax = (fig.axes)
    for i in range(len(ax)):
        print(i, magnet_dim)
        if i >= magnet_dim:
            break
#        ax[i].plot(df.iloc[:,i], df['chi2'], '.', linestyle='none')
        ax[i].hist(df.iloc[:,i], log=True)
#        ax[i].plot(df.loc[df["true_value"]==True].iloc[:,i], df.loc[df["true_value"]==True]['chi2'], 'rx', linestyle='none')
        ax[i].axvline(x=np.sum(df.loc[df["true_value"]==True].iloc[:,i]), color='red')
#        ax[i].set_yscale('log')
#        ax[i].set_xscale('log')
        ax[i].set_title(df.columns[i])
#        ax[i].set_ylim(100,200)
#        xdiff = max(abs(df.iloc[:,i].mean()-df.iloc[:,i].max()),abs(df.iloc[:,i].mean()-df.iloc[:,i].min()))
#        xmin = max(df.iloc[:,i].mean()-2*xdiff,0)
#        xmax = df.iloc[:,i].mean()+2*xdiff
#        ax[i].set_xlim(xmin, xmax)
        print(len(names), len(ax), magnet_dim, names[i])
        if names[i] not in ['h2000','h1750','h2475']:
            ax[i].set_xlim([-0.02,0.1])
        else:
            ax[i].set_xlim([0,0.5])
#        ax[i].axes.yaxis.set_visible(False) 
#        if i < 48:
#            ax[i].axes.xaxis.set_ticklabels([]) 
#        if i % 10 == 0:
#            ax[i].set_ylabel('chi2')
    plt.tight_layout()
#    plt.show()
    plt.savefig('branch_chi2.png')
    fig, ax = plt.subplots()
    dof = 1
    df['chi2'].apply(lambda x: x/(dof)).hist(bins=30,log=False)
    plt.axvline(df.loc[df["true_value"] == True].at[n_sims,'chi2']/(dof), 0, n_sims, color='r', label='true chi2 = {:.2f}'.format(df.loc[df["true_value"] == True].at[n_sims,'chi2']/(dof)))
    ax.set_xlabel('chi2 distribution')
#    x = np.linspace(chi2.ppf(0.01, dof),chi2.ppf(0.99, dof), 100)
#    ax.plot(x/dof, chi2.pdf(x, dof))
    plt.legend()
    plt.savefig('chi2dist.png')
    plt.cla()
    fig, ax = plt.subplots()
    dof = 1
    df['ss_m1_chi2'].apply(lambda x: x/(dof)).hist(bins=30,log=False)
    plt.axvline(df.loc[df["true_value"] == True].at[n_sims,'ss_m1_chi2']/(dof), 0, n_sims, color='r', label='true chi2 = {:.2f}'.format(df.loc[df["true_value"] == True].at[n_sims,'ss_m1_chi2']/(dof)))
    ax.set_xlabel('chi2 distribution')
#    x = np.linspace(chi2.ppf(0.01, dof),chi2.ppf(0.99, dof), 100)
#    ax.plot(x/dof, chi2.pdf(x, dof))
    plt.legend()
    plt.savefig('ss_m1_chi2dist.png')
    plt.cla()
    df['ss_chi2'].apply(lambda x: x/(dof)).hist(bins=30,log=False)
    plt.axvline(df.loc[df["true_value"] == True].at[n_sims,'ss_chi2']/(dof), 0, n_sims, color='r', label='true chi2 = {:.2f}'.format(df.loc[df["true_value"] == True].at[n_sims,'ss_chi2']/(dof)))
    ax.set_xlabel('chi2 distribution')
#    x = np.linspace(chi2.ppf(0.01, dof),chi2.ppf(0.99, dof), 100)
#    ax.plot(x/dof, chi2.pdf(x, dof))
    plt.legend()
    plt.savefig('ss_chi2dist.png')
    plt.cla()
    df['mult_chi2'].apply(lambda x: x/(dof)).hist(bins=30,log=False)
    plt.axvline(df.loc[df["true_value"] == True].at[n_sims,'mult_chi2']/(dof), 0, n_sims, color='r', label='true chi2 = {:.2f}'.format(df.loc[df["true_value"] == True].at[n_sims,'mult_chi2']/(dof)))
    ax.set_xlabel('chi2 distribution')
#    x = np.linspace(chi2.ppf(0.01, dof),chi2.ppf(0.99, dof), 100)
#    ax.plot(x/dof, chi2.pdf(x, dof))
    plt.legend()
    plt.savefig('mult_chi2dist.png')
    print(df.loc[df["chi2"]>df.at[n_sims,'chi2']].shape[0]/df.shape[0])
    print(df.at[n_sims,'chi2']/dof)

#    plt.clf()
#    df['1000bin'].loc[(df['1000bin']<121) & (df['1000bin']>59)].hist(bins=60)
#    plt.axvline(df.loc[df["true_value"] == True].at[n_sims,'1000bin'], 0, n_sims, color='r', label='true count = {:}'.format(df.loc[df["true_value"] == True].at[n_sims,'1000bin']))
#    plt.savefig('1000bindist.png')
#    plt.clf()
#    df['170bin'].loc[(df['170bin']<601) & (df['170bin']>449)].hist(bins=150)
#    plt.axvline(df.loc[df["true_value"] == True].at[n_sims,'170bin'], 0, n_sims, color='r', label='true count = {:}'.format(df.loc[df["true_value"] == True].at[n_sims,'170bin']))
#    plt.savefig('170bindist.png')



if __name__=='__main__':
    main(n_sims)



