#!/mnt/simulations/secarml/soft/anaconda3/bin/python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys, os
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import normalize
from scipy.optimize import lsq_linear, minimize, Bounds, NonlinearConstraint
from multiprocessing import Pool, Process, Queue, current_process, freeze_support
import json

from . import utils as utils
#from utils import load_all, total_minimize, chisq, running_chisq, monte_carlo, plot_all, normalize_sims

# script, first, second = sys.argv
WORK_DIR = '../'
energies = []
def load_configs():

    f = open(WORK_DIR+'config.json')
    config = json.load(f)
    print(config)
    f.close()

    return config

config = load_configs()
GEANT_DIR = "/mnt/analysis/e17028/SuNGEANT/geant/sunNewDetector_{}/rootFiles/".format(config['isotope_daughter'])
NGEANT_DIR = "/mnt/analysis/e17028/SuNGEANT/geant/sunNewDetector_{}/rootFiles/".format(config['isotope_ndaughter'])

true_fit = np.array([1.93467728e-02, 4.64154590e-02, 2.04564634e-01, 2.51136937e-01,
       4.23914069e-02, 3.00546166e-02, 1.32635819e-01, 0.00000000e+00,
       7.12683989e-03, 0.00000000e+00, 3.20977082e-02, 6.85470172e-03,
       1.96794809e-02, 2.18927770e-14, 3.39806745e-02, 0.00000000e+00,
       2.50895912e-02, 1.68366184e-03, 3.45207888e-02, 5.22420991e-03,
       2.80327745e-02, 7.73878312e-03, 1.66384095e-02, 1.29155549e-02,
       8.36737551e-03, 1.23299237e-02, 5.90454122e-03, 8.35875226e-03,
       4.74507340e-04, 2.72550054e-03, 3.00468485e-15, 2.87576222e-03,
       8.34812635e-04, 0.00000000e+00, 3.70935804e-14, 9.72214920e-14])
if len(true_fit) != len(energies):
    if len(true_fit) < len(energies):
        for i in range(len(energies)-len(true_fit)):
            true_fit.append(0)
    else:
        true_fit = true_fit[:len(energies)]
if len(true_fit) != len(energies):
    print("still not equal")

def parallel_eval(i):

    # load data
    df_tas, df_tas_sim, df_ss, df_ss_sim, df_mult, df_mult_sim, df_ss_m1, df_ss_m1_sim = utils.load_all()

    total_n_tas = np.sum(df_tas['content'])
    total_n_ss = np.sum(df_ss['content'])
    total_n_mult = np.sum(df_mult['content'])
    total_n_ss_m1 = np.sum(df_ss_m1['content'])
#    df_tas = df_tas.iloc[:400,:]
#    df_ss = df_ss.iloc[:750,:]

    # make variables for optimization
    test = np.zeros(df_ss_sim.shape[1])
    ensdf_lvl = np.array([23.8, 5, 1.1, 26, 26.6,0,0,0, 17.3])/100.0
    for test_i in range(len(ensdf_lvl)):
        test[test_i] += ensdf_lvl[test_i]
    bounds = Bounds(0,1) 
    con = lambda x: np.sum(x)-1
    nlc = NonlinearConstraint(con, 0.999, 1.001)
    con1 = ({'type': 'eq', 'fun': con})

    # sample data from observation
    np.random.seed()
#    print(df_tas['content'].iloc[8])
#    print(df_tas.loc[df_tas['content']<4])
#    print(df_ss.loc[df_ss['content']<4])
#    print(df_mult.loc[df_mult['content']<4])
    tas_rand = np.random.poisson(df_tas['content'])
    ss_rand = np.random.poisson(df_ss['content'])
    mult_rand = np.random.poisson(df_mult['content'])
    ss_m1_rand = np.random.poisson(df_ss_m1['content'])
#    tas_rand = np.random.normal(df_tas['content'], df_tas['error'],size=df_tas.shape[0])
#    ss_rand = np.random.normal(df_ss['content'], df_ss['error'],size=df_ss.shape[0])
#    mult_rand = np.random.normal(df_mult['content'], df_mult['error'],size=df_mult.shape[0])
    if i > -1:
        df_tas = utils.monte_carlo(df_tas,tas_rand,i)
        df_ss = utils.monte_carlo(df_ss,ss_rand,i)
        df_ss_m1 = utils.monte_carlo(df_ss_m1,ss_m1_rand,i)
        df_mult = utils.monte_carlo(df_mult,mult_rand,i)
    else:
        print("not randomizing")

#    print(np.sum(df_tas['content']), total_n_tas)
    total_i_tas = np.sum(df_tas['content'])
    total_i_ss = np.sum(df_ss['content'])
    total_i_ss_m1 = np.sum(df_ss_m1['content'])
    total_i_mult = np.sum(df_mult['content'])
    df_tas_sim = utils.normalize_sims(df_tas_sim, total_i_tas)
    df_ss_sim = utils.normalize_sims(df_ss_sim, total_i_ss)
    df_ss_m1_sim = utils.normalize_sims(df_ss_m1_sim, total_i_ss_m1)
    df_mult_sim = utils.normalize_sims(df_mult_sim, total_i_mult)
#    df_tas['error'] = df_tas['content'].apply(lambda x: np.sqrt(x))
#    df_ss['error'] = df_ss['content'].apply(lambda x: np.sqrt(x))
#    df_mult['error'] = df_mult['content'].apply(lambda x: np.sqrt(x))
#    print(np.sum(df_tas['content']), total_n_tas)
#    print(np.sum(df_ss['content']), total_n_ss)
#    print(np.sum(df_mult['content']), total_n_mult)
#    print(i,type(i))
    chi2 = 0
#    if i == 0:
#        print("chi2: {}, reduced_chi2: {}".format(chi2, chi2/reduce_factor))
    # run minimization
    popt_all = minimize(utils.total_minimize, test, args=(df_tas, df_tas_sim, df_ss, df_ss_sim, df_mult, df_mult_sim, df_ss_m1, df_ss_m1_sim, 'all'), bounds=bounds, constraints=con1)
    chi2 += utils.chisq(popt_all.x, df_tas_sim, df_tas)
    reduce_factor = max(df_tas.shape[0]-test.shape[0],1)
    ss_chi2 = utils.chisq(popt_all.x, df_ss_sim, df_ss)
    ss_m1_chi2 = utils.chisq(popt_all.x, df_ss_m1_sim, df_ss_m1)
#    reduce_factor = max(df_ss.shape[0]+df_tas.shape[0]-test.shape[0],1)
    mult_chi2 = utils.chisq(popt_all.x, df_mult_sim, df_mult)
#    reduce_factor = max(df_ss.shape[0]+df_tas.shape[0]+df_mult.shape[0]-test.shape[0],1)
    return_array = popt_all.x
    return_array = np.append(return_array, [df_tas.at[int(170/25),'content']])
#    print([df_tas.at[int(170/25),'content']])
    return_array = np.append(return_array, [df_tas.at[int(1000/25),'content']])
    return_array = np.append(return_array, [chi2])
    return_array = np.append(return_array, [ss_chi2])
    return_array = np.append(return_array, [mult_chi2])
    return_array = np.append(return_array, [ss_m1_chi2])

    data_coef=popt_all.x

    return return_array 


def main():

    # load data
    df_tas, df_tas_sim, df_ss, df_ss_sim, df_mult, df_mult_sim, df_ss_m1, df_ss_m1_sim = utils.load_all()

#    df_tas = df_tas.iloc[:400,:]
#    df_ss = df_ss.iloc[:750,:]

    # make variables for optimization
    test = np.zeros(df_ss_sim.shape[1])
    ensdf_lvl = np.array([23.8, 5, 1.1, 26, 26.6,0,0,0, 17.3])/100.0
    for i in range(len(ensdf_lvl)):
        test[i] += ensdf_lvl[i]
    bounds = Bounds(0,1) 
    con = lambda x: np.sum(x)-1
    nlc = NonlinearConstraint(con, 0.999, 1.001)
    con1 = ({'type': 'eq', 'fun': con})

    # run minimization
#    df_tas['error'] = df_tas['content'].apply(lambda x: np.sqrt(x))
#    df_ss['error'] = df_ss['content'].apply(lambda x: np.sqrt(x))
#    df_mult['error'] = df_mult['content'].apply(lambda x: np.sqrt(x))
    popt_all = minimize(utils.total_minimize, test, args=(df_tas, df_tas_sim, df_ss, df_ss_sim, df_mult, df_mult_sim, df_ss_m1, df_ss_m1_sim, 'all'), bounds=bounds, constraints=con1)
    data_coef=popt_all.x.reshape((1,df_ss_sim.shape[1]))
#    data_coef=np.append(popt_all.x, popt_all.fun).reshape(1,df_ss_sim.shape[1]+1)
#    data_coef=np.array([data_coef[-1][-1]])
#    print(popt_all.x)
    data_chi = popt_all.fun * (df_tas.shape[0]-test.shape[0]) 
#    data_chi = popt_all.fun * (df_tas.shape[0]+df_ss.shape[0]-test.shape[0]) 
            
    # set up and run parallel minimization for n_sims samples
    NUMBER_OF_PROCESSES = 10
    n_sims = int(5e4)
    pool = Pool(processes = NUMBER_OF_PROCESSES)
    r = pool.map_async(parallel_eval, range(n_sims))
    r.wait()
    r_return = np.array(r.get())
#    print(r_return.shape)
    sim_coef = r_return[:,:-4].reshape((n_sims,df_ss_sim.shape[1]+2))
    print(sim_coef, sim_coef.shape)
    sim_chi = r_return[:,-4]   
    sim_ss_chi = r_return[:,-3]   
    sim_mult_chi = r_return[:,-2]   
    sim_ss_m1_chi = r_return[:,-1]   
    print(sim_chi, sim_chi.shape )

    sim_columns = []
    for i in range(len(df_tas_sim.columns)):
        sim_columns.append(df_tas_sim.columns[i])
    sim_columns.append('170bin')
    sim_columns.append('1000bin')
    df_sims = pd.DataFrame(sim_coef, columns=sim_columns)
    df_sims['chi2'] = sim_chi
    df_sims['ss_chi2'] = sim_ss_chi
    df_sims['mult_chi2'] = sim_mult_chi
    df_sims['ss_m1_chi2'] = sim_ss_m1_chi
    ss_chi = utils.chisq(popt_all.x, df_ss_sim, df_ss)
    ss_m1_chi = utils.chisq(popt_all.x, df_ss_m1_sim, df_ss_m1)
    mult_chi = utils.chisq(popt_all.x, df_mult_sim, df_mult)
    data_coef = np.append(data_coef, [df_tas.at[int(170/25),'content']])
    data_coef = np.append(data_coef, [df_tas.at[int(1000/25),'content']])
    print(data_coef)
    data_coef=data_coef.reshape((1,df_ss_sim.shape[1]+2))
    df_true = pd.DataFrame(data_coef, columns=sim_columns)
    df_true['chi2'] = data_chi
    df_true['ss_chi2'] = ss_chi
    df_true['ss_m1_chi2'] = ss_m1_chi
    df_true['mult_chi2'] = mult_chi
    df_true['true_value'] = True
##    df_sims['chi2'] = sim_chi
    df_sims['true_value'] = False
    df_sims = pd.concat([df_sims,df_true], ignore_index=True)
    df_sims.to_hdf('correct_chi2mc_{}.h5'.format(n_sims), key='df')
    print(df_sims)
#    bins = np.linspace(0,0.3,30)
#    fig, axs = plt.subplots(4,6)
#    fig.delaxes(axs[3][5])
#    axs = (fig.axes)
##    print(len(axs))
#    hists = df_sims.iloc[:,:23].hist(bins=bins,grid=False,ax=axs)
#    axs2 = hists
##    fig.tight_layout()
#    fig.set_figheight(10)
#    fig.set_figwidth(19)
##    axs2 = axs
#    for i in range(len(axs)):
#        axs[i].axvline(x=data_coef[i], color='red')
#        if i not in [3,4,7]:
#            axs[i].set_xlim([-0.02,0.1])
#    plt.savefig("mc_{}.png".format(n_sims))
    return



if __name__ == "__main__":
    main()
