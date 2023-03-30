#!/mnt/misc/sw/x86_64/all/anaconda/python3.7/bin/python

# make sure above path points to the version of python where you have pygmo installed 
# nscl servers
#!/mnt/misc/sw/x86_64/all/anaconda/python3.7/bin/python
# hpcc servers
#!/mnt/home/herman67/anaconda3/envs/pygmo/bin/python

# import commands
import numpy as np
import pandas as pd
import os, sys

# specify important directories and names
PYGMO_DIR = "../"
OUTPUT_DIR = PYGMO_DIR
generations = 750
population_size = 969
batch = 230
# specify number of magnets
magnet_dim = 36
names = []
hist = np.loadtxt(PYGMO_DIR + 'levels_new.txt')
for h in hist:
	names.append('h'+str(int(h)))

def main(gens=generations, batch_id=210):

    OUTPUT_PREFIX = OUTPUT_DIR + 'py/output_4f_moead_FP2_FP3_{}_'.format(gens)
    # set up columns for dataframe
#    quads = []
#    for i in range(magnet_dim):
#        quads.append("h{}".format(i+1))
    columns = names 
    objs = ['tas','ss','mult','tas_m1','total','branch_sum']
    for i in range(len(objs)):
        columns.append("{}".format(objs[i]))
    
    # i run batches in 10s, so i specify the first id of the batch
    start_i = batch_id 
    end_i = start_i + 1
    # name the output
    db_out = OUTPUT_DIR + "secar_4d_db_{}s.h5".format(start_i)
    
    # check if we want to append to an existing df
    df = None
    #if os.path.exists(db_out):
    #    print('appending')
    #    df = pd.read_hdf(db_out)    
    #    print(df)
    
    # load results from all islands
    for i in range(start_i, end_i):
        df_new = pd.read_csv('{}{}.csv'.format(OUTPUT_PREFIX,i),names=columns)
        max_obj = 1
        # append to existing df
        if i > start_i: # or os.path.exists(db_out):
            df = df.append(df_new,ignore_index=True)
        # initialize df
        else:
            df = df_new
    
    print("percent done: {:.2f}%".format(len(df.index)/((generations+1)*population_size)*10))
    # write df to h5
    max_obj = 1 
#    df['chi2'] = df['tas'] + df['ss']
    obj_val = [592.3877639215518, 1914.312304257408, 83.14181256225928, 634092.0660473696]
    obj_val = np.zeros(4)+1

    # check for solutions strictly better than nominal (all objs < 1)
    query_txt = '' 
    for i in range(len(objs)):
#        df[objs[i]] = df[objs[i]].apply(lambda x: x/obj_val[i])
        query_txt += objs[i] + "<{}".format(max_obj)
        if i < len(objs)-1:
            query_txt+="&"
        else:
            query_txt = objs[i] + "<{}".format(max_obj)
    print(query_txt)
    print(df)
#    df = df.query(query_txt)
    print(df)
    df.to_hdf(db_out,key='df')
    print(df)
    print(np.sum(df.iloc[-1,:44]))
    return

if __name__=='__main__':
    inputs = sys.argv
    print(inputs)
    if len(inputs) > 1:
        generations = int(inputs[1])
    if len(inputs) > 2:
        batch = int(inputs[2])
    main(generations,batch)
