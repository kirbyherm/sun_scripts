#!/mnt/misc/sw/x86_64/all/anaconda/python3.7/bin/python
# Author: Narine Kokhlikyan <narine@slice.com>
# License: BSD

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

from matplotlib.ticker import NullFormatter
from sklearn import manifold, datasets
from time import time

n_samples = 150
n_components = 2
names = []
hist = np.loadtxt('levels.txt')
for h in hist:
	names.append('h'+str(int(h)))


def plot_tsne(filename):

    (fig, subplots) = plt.subplots(3, 3, figsize=(8, 8))
    perplexity = 30
    perplexities = [perplexity]
#    perplexities = [100,1000,10000,100000]
    
    X, y = datasets.make_circles(
        n_samples=n_samples, factor=0.5, noise=0.05, random_state=0
    )
    
    df = pd.read_hdf(filename)
#    df = df.loc[df['kcluster']==3]
    print(df)
    plot_combos = [[0,1],[2,3],[3,4],[2,4],[4,6],[6,9],[9,14],[12,14]]

    X = np.array(df[names])
    print(X)
    y = np.array(df['kcluster'])
    z = np.array(df['closest'])
    number_of_clusters = np.max(y+1)
    colors = list(plt.get_cmap('tab20').colors)
    clusters = []
    for i in range(number_of_clusters):
        clusters.append(y == i)
    
    close = z == True
    
    
#    print (X,y[close])
#
    plot_x, plot_y = 0, 0
#    for pair in plot_combos:
#    
#        ax = subplots[plot_y][plot_x]
#        for i in range(number_of_clusters):
#            if i < 3:
#                continue
#            ax.scatter(X[clusters[i],pair[0]], X[clusters[i],pair[1]], c=np.array(colors[i+1]).reshape(1,-1))
#        for i, txt in enumerate(y[close]):
#            if i < 3:
#                continue
#            ax.annotate(txt+1, (X[close,pair[0]][txt], X[close,pair[1]][txt]))
##        ax.xaxis.set_major_formatter(NullFormatter())
##        ax.yaxis.set_major_formatter(NullFormatter())
#        ax.set_title("{1}:{0}".format(magnets[pair[0]], magnets[pair[1]]))
#        plt.axis("tight")
#        plot_x += 1
#        if plot_x > 2:
#            plot_y+=1
#            plot_x=0    
#    
    for i, perplexity in enumerate(perplexities):
    
        ax = subplots[plot_y][plot_x]
    
        t0 = time()
        tsne = manifold.TSNE(
            n_components=n_components,
            init="random",
            random_state=0,
            perplexity=perplexity,
            learning_rate=10.0,
            n_iter=5000,
        )
        Y = tsne.fit_transform(X)
        t1 = time()
        print("circles, perplexity=%d in %.2g sec" % (perplexity, t1 - t0))
        ax.set_title("Perplexity=%d" % perplexity)
        for j in range(number_of_clusters):
#            if j < 3:
#                continue
            ax.scatter(Y[clusters[j],0], Y[clusters[j],1], c=np.array(colors[j+1]).reshape(1,-1))
        for k, txt in enumerate(y[close]):
#            if k < 3:
#                continue
            print(k,txt)
            ax.annotate(txt+1, (Y[close,0][txt], Y[close,1][txt]))
        ax.xaxis.set_major_formatter(NullFormatter())
        ax.yaxis.set_major_formatter(NullFormatter())
        ax.axis("tight")
    plt.savefig(filename+'_tsne1.png')
    return


if __name__=='__main__':
    inputs = sys.argv
    plot_tsne(inputs[1])


