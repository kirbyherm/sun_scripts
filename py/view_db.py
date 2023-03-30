#!/mnt/misc/sw/x86_64/all/anaconda/python3.7/bin/python

# make sure above path points to the version of python where you have pygmo installed 
# nscl servers
#!/mnt/misc/sw/x86_64/all/anaconda/python3.7/bin/python
# hpcc servers
#!/mnt/home/herman67/anaconda3/envs/pygmo/bin/python


#import commands
import pandas as pd
import numpy as np
import os,sys
from kmeans_cluster import run_kmeans
#from load_data import run


# set pandas view options to print everything
pd.set_option("max_rows", None)
pd.set_option("max_columns", None)

# set important directories and files
PYGMO_DIR = "../"
OUTPUT_DIR = PYGMO_DIR

Qnom = np.array([-0.40033,0.219852,0.2552369,-0.246677876,0.11087109,0.175336731,-0.0268214976,-0.14859,0.2855,-0.0335,0.149432825,-0.182,0.1910,0.12900,-0.1380,0,0,0,0])
# Q1H:=0.003703;
# H1:=0.0103564;
# H2:=0.0052735{*0.5};
# H3:=-0.008774463{*1.5};
# O1:=0.031283{*2.0};

# function i found which constructs a pareto front from an input set of points
def is_pareto_efficient_simple(costs):
    """
    Find the pareto-efficient points
    :param costs: An (n_points, n_costs) array
    :return: A (n_points, ) boolean array, indicating whether each point is Pareto efficient
    """
    is_efficient = np.ones(costs.shape[0], dtype = bool)
    for i, c in enumerate(costs):
#        print(i,c,is_efficient[i])
        if is_efficient[i]:
            is_efficient[is_efficient] = np.any(costs[is_efficient]<c, axis=1)  # Keep any point with a lower cost
            is_efficient[i] = True  # And keep self
    return is_efficient

# only show best [#] since we get a lot of points
show_best = 10
batch = 210
kclusters = 5

def main(start_i=batch):
    # specify database for input
    db_out = OUTPUT_DIR + "secar_4d_db_{}s.h5".format(start_i)
    # initialize empty df
    df = None
    # list the objective names
    objectives = ['tas','ss','mult','tas_m1','total','branch_sum']
    # read df
    if os.path.exists(db_out):
        print('opening existing db')
        df = pd.read_hdf(db_out)    
    
    # drop na if any
    df = df.dropna()

    df = df.loc[df['branch_sum'] < 1]    
    # restrict the df to only the points that fit the problem constraints
    #   (can also change this to any value, e.g. 1 to show only better than nominal)
    magnet_dim = 36
    reduc_factors = np.array([137-magnet_dim,137-magnet_dim,70,137-magnet_dim,137+137+10-magnet_dim,1])
    max_obj = np.multiply(np.array([2.5,4.0,1.5,1.0,5.0,1e5]), reduc_factors)
    max_obj = np.multiply(np.array([2.5,10.0,10.0,1.0,10.0,1e5]), reduc_factors)
    max_obj = np.zeros(6)+1e3
    max_obj[4] = 2500
#    df = df.loc[(df['FP1_res'] < max_obj) & (df['FP2_res'] < max_obj) & (df['MaxBeamWidth'] < max_obj) & (df['FP3_res'] < max_obj) & (df['FP4_BeamSpot'] < max_obj)]
#    df = df.loc[(df['FP2_res'] < max_obj) & (df['MaxBeamWidth'] < max_obj) & (df['FP3_res'] < max_obj) & (df['FP4_BeamSpot'] < max_obj)]
    query_txt = '' 
    for i in range(len(objectives)):
        query_txt += objectives[i] + "<{}".format(max_obj[i])
        if i < len(objectives)-1:
            query_txt+="&"
    print(query_txt)
    df = df.query(query_txt)
    
    # get costs and pass to pareto function
#    costs = df[['FP1_res','FP2_res','FP3_res','MaxBeamWidth','FP4_BeamSpot']]
    costs = df[objectives]
#    print(df[objectives])
    costs = np.array(costs)
#    pareto = is_pareto_efficient_simple(costs)
#    # add pareto column to df
#    df['pareto'] = pareto
#    print("pareto front points: {}".format(np.count_nonzero(pareto) ))
#    # restrict df to only those points on the pareto front
#    df = (df.loc[(df['pareto']==True)])
#    print(df.iloc[:100,-5:-1])
    # additional code which provides an alternate way of sorting/filtering df, based on sum-squares of objs
    df_objs = df[objectives]
    df_objs['ssobjs'] = (df_objs**2).sum(axis=1)
    df['ssobjs'] = df_objs['ssobjs']
##    print(df_objs)
#    df['ssobjs'] = np.sqrt(df['FP1_res']**2+df['FP2_res']**2+df['FP3_res']**2+df['MaxBeamWidth']**2+df['FP4_BeamSpot']**2)
#    df = df.sort_values(by='ssobjs',ignore_index=True)
##    df = df.loc[df['ssobjs'] < df['ssobjs'][show_best]]
##    print(df)
#    quads = df.columns
##    print(quads)
#    for q in range(len(quads)):
#        print(q, quads[q])
#        if "q" in quads[q]:
#            magnet_dim = q+1
    df = df.reset_index(drop=True)
#    df = run_kmeans(df, magnet_dim, kclusters)
#    
#    # print objective values for [show_best] number of points, sorted by FP4_BeamSpot
#    print_columns = objectives
#    print_columns.append('closest')
#    print_columns.append('ssobjs')
#    print_columns.append('kcluster')
##    print(df.loc[:show_best,['FP1_res','FP2_res','FP3_res','MaxBeamWidth','FP4_BeamSpot','closest','ssobjs','kcluster']].sort_values(by=['ssobjs']))
#    print(df.loc[:show_best,print_columns].sort_values(by=['ssobjs']))
#    # sort df by FP4_BeamSpot values, and reindex
    df = df.sort_values(by='ssobjs',ignore_index=True,ascending=False)
    df = df.reset_index(drop=True)
#    # print the magnet scale factors for the best FP4_BeamSpot points
#    write_qnames = ['q1','q2','q3','q4','q5','q6','q7','q8','q9','q10','q11','q12','q13','q14','q15','h1','h2','h3','o1']
#    (np.power(2,df.loc[df['closest']==True].iloc[:,:19])).round(5).to_csv('magnet_factors.csv',header=write_qnames,index=False)
#    (df.iloc[:,:].sort_values(by='ssobjs')).round(5).to_csv('all_good_points.csv',header=df.columns,index=False)
#    # write only the magnet values and objective values to df
##    print(df.columns)
#    df = df.drop('pareto',1)
    #df = df.iloc[:,:19]
    df.to_hdf('best{}.h5'.format(start_i),key='df')
    print(df.iloc[-100:,:44])
    print(df.iloc[-100:,-6:])
#    run(df.iloc[-1,:44],"1")
    return

if __name__=='__main__':
    inputs = sys.argv
    print(inputs)
    if len(inputs) > 1:
        batch = int(inputs[1])
    main(batch)
