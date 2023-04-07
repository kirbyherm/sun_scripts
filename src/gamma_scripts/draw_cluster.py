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
#from problem import optimizeRes
import pandas as pd

# specify Tex details for pretty plots
os.environ['PATH'] = os.environ['PATH'] + ':/mnt/misc/sw/indep/all/texlive/2013/bin/x86_64-linux/latex'
#os.environ['PATH'] = os.environ['PATH'] + ':/usr/bin/tex'
#plt.rcParams.update({
#    "text.usetex": True,
#})

script, filename = sys.argv
optimized_params = 4
fNom = np.zeros(optimized_params)+1
fNames = [r"{FP1-res}${}^{-1}$",r"{FP2-res}${}^{-1}$",r"{FP3-res}${}^{-1}$",r"MaxBeamWidth",r"BeamSpotSize"]
fNames = fNames[:optimized_params]
#print(len(fNames))
magnet_dim = 19

# read pop from h5 file (i.e. after running view_db.py)
def read_pop_df(filename, pop=None):
    df = pd.read_hdf(filename)

    magnet_dim = 19
    p_optimizeRes = pg.problem(optimizeRes(magnet_dim))
    nobj = p_optimizeRes.get_nobj()
    if pop == None:
        pop = pg.population(p_optimizeRes)
    nrow, ncol = df.shape
    for i in df.index:
        append=True
        xs = []
        for j in range(1,magnet_dim+1):
            xs.append(df["q"+str(j)][i]) 
        xs = np.asarray(xs)
        fs = []
        for j in range(magnet_dim,magnet_dim+nobj):
#            if i == 0:
#                print(df.iloc[i,magnet_dim:magnet_dim+p_optimizeRes.get_nobj()])
            fs.append(df.iloc[i,j])
        if append:
            pop.push_back(xs,f=fs)
    return pop, df


def plot_4d(popi,filename,df):

    sort_param = 4
    good_results=0
    magnet_dim = len(popi.get_x()[0])
#    hv = pg.hypervolume(popi)
#    ref_point = np.zeros(optimized_params)+1e10 
    ndf, dl, dc, ndl = pg.fast_non_dominated_sorting(popi.get_f())
    plot_x, plot_y = 0,0
    fig, axs = plt.subplots(optimized_params-1,sharex=True)
    fig.suptitle('Pareto Fronts of each parameter vs. BeamSpotSize at FP4')
    axs[optimized_params-2].set_xlabel(fNames[sort_param])
    reduced_ndf = []
    first = True
    df_closest = df
    cmap = matplotlib.colors.ListedColormap(matplotlib.cm.get_cmap("Pastel1").colors[:3])

    for j in range(0,optimized_params-1):
        ndf_champ = []
        axs[plot_y].axvline(x=fNom[0],linestyle="dashed",color="red")
        axs[plot_y].axhline(y=1.0,linestyle="dashed",color="red")
        for i in ndf[0]:
            check_val=1e9
            if np.all(np.array(popi.get_f()[i]) < check_val) == True:
                good_results+=1
#                print(filename[-6:-4], good_results, popi.get_f()[i])
                ndf_champ.append(popi.get_f()[i])
                reduced_ndf.append(i)
        try:
            pg.plot_non_dominated_fronts(ndf_champ,comp=[sort_param,j],axes=axs[plot_y])
#            if j == 1:
#                print(filename[-6:-4], len(ndf_champ))
        except:
            print(filename[-6:-4], "no better than nominal solutions")
            return
        if "closest" in df.columns:
            df_closest = df.loc[df['closest']==True]
            df_closest = df_closest.reset_index(drop=True)
#            print(df_closest.iloc[:,15:19])
            for i_closest in df_closest.index:
                axs[plot_y].text(df_closest.iloc[:,magnet_dim+sort_param][i_closest],df_closest.iloc[:,magnet_dim+j][i_closest],str(i_closest+1),color='red')
        axs[plot_y].set_ylabel(fNames[j])
        axs[plot_y].set_yscale('log')
#        print(axs[plot_y].get_xlim(),axs[plot_y].get_ylim())
        max_y = max(abs(np.log10(axs[plot_y].get_ylim())))
        axs[plot_y].set_ylim(np.power(10,-max_y),np.power(10,max_y))
        axs[plot_y].set_xscale('log')
        max_x = max(abs(np.log10(axs[plot_y].get_xlim())))
        axs[plot_y].set_xlim(np.power(10,-max_x),np.power(10,max_x))
#        print(math.ceil(max(axs[plot_y].get_ylim())),math.ceil(max(axs[plot_y].get_ylim())))
        colors = np.zeros((int(math.ceil(max(axs[plot_y].get_ylim()))),int(math.ceil(max(axs[plot_y].get_xlim()))))) + 1
        for colori in range(1,colors.shape[0]):
            for colorj in range(1, colors.shape[1]):
                colors[colori,colorj] = 0
        colors[0,0] = 2
#        print(colors)
        axs[plot_y].pcolormesh(colors, cmap=cmap)
        plot_y += 1

    fig.tight_layout()
    fig.savefig(filename+"_paretos.png")
    plt.cla()
    fig2, axs2 = plt.subplots(5,4)
    plot_x, plot_y = 0,0
    reduced_qs = np.array(popi.get_x())[reduced_ndf]
    ycolumns = []
    for i in range(magnet_dim):
        ycolumns.append('y{}'.format(i))
    df = pd.DataFrame(reduced_qs, columns = ycolumns)
    qNom = np.zeros(magnet_dim)
    write_qnames = ['q1','q2','q3','q4','q5','q6','q7','q8','q9','q10','q11','q12','q13','q14','q15','h1','h2','h3','o1']
    for i in range(magnet_dim):
        if plot_x > 3:
            plot_x, plot_y = 0, plot_y+1 
        axs2[plot_y,plot_x] = df['y{0}'.format(i)].plot.hist(ax=axs2[plot_y,plot_x],bins=100,range=(-3,3))
#        axs2[plot_y,plot_x].axvline( x = qNom[i], ymin=0,ymax=20,color='red',linestyle='dashed')
#        axs[plot_y,plot_x].axvline( x = max_y[i], ymin=0,ymax=20,color='green',linestyle='dashed')
        axs2[plot_y,plot_x].axes.yaxis.set_visible(False)
#        axs2[plot_y,plot_x].axes.set_yscale("log")
#        axs2[plot_y,plot_x].axes.set_ylim(0.1,20)
        xlower, xupper = popi.problem.get_bounds()
        xlower, xupper = np.min(xlower), np.max(xupper)
        axs2[plot_y,plot_x].axes.set_xlim(xlower,xupper)
        axs2[plot_y,plot_x].set_title("{0}".format(write_qnames[i]))
        y_min, y_max = axs2[plot_y,plot_x].get_ylim()
        df_closest['yplot'] = pd.Series(df_closest.index).apply(lambda x: x/len(df_closest.index)*(y_max-y_min)+y_min)
#        print(df_closest.iloc[:,:15])
#        if "closest" in df.columns:
        for i_closest in df_closest.index:
            axs2[plot_y,plot_x].text(df_closest["q{0}".format(i+1)][i_closest],df_closest['yplot'][i_closest],str(i_closest+1),color='red')
        
        plot_x += 1
    
    axs2[plot_y,plot_x].axis('off')
    axs2[plot_y,plot_x].text(0.3, 0.5, 'x-axis is in log2', horizontalalignment='center', verticalalignment='center', transform=axs2[plot_y,plot_x].transAxes)
    fig2.tight_layout()
    plt.savefig(filename + "_magnet_hists.png")
    return

def main(filename):

    file_extension = os.path.splitext(filename)[-1]
    print(os.path.split(filename))
    popi = None
    print("reading {} file".format(file_extension))
#    if file_extension == ".h5":
#        popi, df_list = read_pop_df(filename)
#    for df in df_list:
#        plot_4d(popi,filename,df)

    df = pd.read_hdf(filename)
    number_of_clusters = np.max(df['kcluster']+1)
    colors = list(plt.get_cmap('tab20').colors)
    for obj in ['tas','ss']:
        ax = None
        for i in range(number_of_clusters):
            print(df.loc[df['closest']==True])
            df_clos = df.loc[df['closest']==True].sort_values(by=['mult']).reset_index(drop=True)
            print(df_clos.loc[df_clos['kcluster']==i].index)
            if i == 0:
                ax = df.loc[(df['kcluster']==i)].plot(x='mult',y=obj,style='o',color=np.array(colors[i+1]),label=df_clos.loc[df_clos['kcluster']==i].index[0]+1,markersize=3.0)
#                ax.text(df_clos.loc[df_clos['kcluster']==i]['FP4_BeamSpot'],df_clos.loc[df_clos['kcluster']==i][obj],str(df_clos.loc[df_clos['kcluster']==i].index[0]),color=colors[i],backgroundcolor='w')
#                ax.plot([ax.axes.get_xlim()[0],df_clos.loc[df_clos['kcluster']==i]['FP4_BeamSpot']],[ax.axes.get_ylim()[0],df_clos.loc[df_clos['kcluster']==i][obj]],color=colors[i])
            else:
                ax = df.loc[(df['kcluster']==i)].plot(x='mult',y=obj,style='o',color=np.array(colors[i+1]),ax=ax,label=df_clos.loc[df_clos['kcluster']==i].index[0]+1,markersize=3.0)
#                ax.plot([ax.axes.get_xlim()[0],df_clos.loc[df_clos['kcluster']==i]['FP4_BeamSpot']],[ax.axes.get_ylim()[0],df_clos.loc[df_clos['kcluster']==i][obj]],color=colors[i])
#                ax.text(df_clos.loc[df_clos['kcluster']==i]['FP4_BeamSpot'],df_clos.loc[df_clos['kcluster']==i][obj],str(df_clos.loc[df_clos['kcluster']==i].index[0]),color=colors[i],backgroundcolor='w')
                ax.axes.set_ylabel(obj)
                ax.axes.set_xlabel('mult')
                ax.axes.set_yscale('log')
                ax.axes.set_xscale('log')
        for i in range(number_of_clusters):
            plt.plot(df.loc[(df['kcluster']==i)].mean()['mult'],df.loc[(df['kcluster']==i)].mean()[obj],marker='s',fillstyle='none',color='b',markersize=20.0,label='_nolegend_')
            ax = df.loc[(df['closest']==True) & (df['kcluster']==i)].plot(x='mult',y=obj,style='o',fillstyle='none',color='b',ax=ax,markersize=20.0,label='_nolegend_')
            plt.plot(df.loc[(df['kcluster']==i)].mean()['mult'],df.loc[(df['kcluster']==i)].mean()[obj],marker='x',color=np.array(colors[i+1]),markersize=20.0,label='_nolegend_', markeredgewidth=5.0)
            ax = df.loc[(df['closest']==True) & (df['kcluster']==i)].plot(x='mult',y=obj,style='x',color=np.array(colors[i+1]),ax=ax,markersize=20.0,label='_nolegend_', markeredgewidth=5.0)
        plt.tight_layout()
        plt.savefig(filename+'_'+obj+'.png')


if __name__=='__main__':
    main(filename)



