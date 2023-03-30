import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

script, column = sys.argv
names = ['chi2']
hist = np.loadtxt('levels.txt')
for h in hist:
	names.append('h'+str(int(h)))
df = pd.read_csv('tas_chi2.csv',names=names)
ymin = (np.min(df['chi2']))
df.plot.scatter(int(column),'chi2')
#xmin, xmax, ymin, ymax = plt.axis()
plt.ylim(ymin-100,ymin+1000)
#plt.yscale('log')
plt.show()
#df.plot.scatter(2,'chi2')
#plt.yscale('log')
#plt.show()
#df.plot.scatter(3,'chi2')
#plt.yscale('log')
#plt.show()
#df.plot.scatter(4,'chi2')
#plt.yscale('log')
#plt.show()
#df.plot.scatter(5,'chi2')
#plt.yscale('log')
#plt.show()
#df.iloc[4184,:]
#np.sum(df.iloc[4184,:])
print(np.sum(df.iloc[-1,1:]))
