#!/mnt/simulations/secarml/soft/anaconda3/bin/python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys, os
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import normalize
from scipy.optimize import lsq_linear, minimize, Bounds, NonlinearConstraint
from datetime import datetime

import shutil
from . import utils as utils

import json

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

# script, first, second = sys.argv
WORK_DIR = '../'
energies = []
#energies = np.loadtxt(WORK_DIR+'levels_orig.txt',dtype='str')[:46]
f = open(WORK_DIR+'config.json')
bin_widths = json.load(f)
tas_bin_width = bin_widths['tas_bin_width']
tas_lower_cut = bin_widths['tas_lower_limit']
tas_upper_cut = bin_widths['tas_upper_limit']
ss_bin_width = bin_widths['ss_bin_width']
ss_lower_cut = bin_widths['ss_lower_limit']
ss_upper_cut = bin_widths['ss_upper_limit']
energies = bin_widths['levels']
isotope = bin_widths['isotope']
print(bin_widths)
f.close()
config=bin_widths


def main(samp=0):

    # datetime object containing current date and time
    now = datetime.now()
    dt_string = now.strftime("%d%m%Y-%H.%M.%S")
    os.makedirs(dt_string)
    df_tas, df_tas_sim, df_ss, df_ss_sim, df_mult, df_mult_sim, df_ss_m1, df_ss_m1_sim = utils.load_all(samp)
#    df_tas, df_tas_sim, df_ss, df_ss_sim, df_mult, df_mult_sim = utils.load_all()

#    df_tas = df_tas.iloc[:150,:]
#    df_ss = df_ss.iloc[:700,:]
#    df_tas = rebin_tas(df_tas)
#    df_tas_sim = rebin_tas_sim(df_tas_sim)
#    df_ss = rebin_ss(df_ss)
#    df_ss_sim = rebin_ss_sim(df_ss_sim)
#    print(df_tas.sum(axis=0), df_tas_og.sum(axis=0))
#    plt.plot(df_tas['keV'], df_tas['content'], 'x')
#    plt.show()


#    return
            
    df_all = pd.concat([df_tas,df_ss,df_mult])
    df_all_sim = pd.concat([df_tas_sim,df_ss_sim,df_mult_sim])
    n_sims = int(1e5)
    sim_coef = []
    test = np.zeros(len(energies))
    ensdf_lvl = np.array([23.8, 5])/100.0
    for i in range(len(ensdf_lvl)):
        test[i] += ensdf_lvl[i]
    bounds = Bounds(0,1) 
    con = lambda x: np.sum(x)-1
    nlc = NonlinearConstraint(con, 0.999, 1.001)
    con1 = ({'type': 'eq', 'fun': con})
    fit_param = 'm1'
    popt_all = minimize(utils.total_minimize, test, args=(df_tas, df_tas_sim, df_ss, df_ss_sim, df_mult, df_mult_sim, df_ss_m1, df_ss_m1_sim, fit_param), bounds=bounds, constraints=con1)
#    popt_all = minimize(utils.total_minimize, test, args=(df_tas, df_tas_sim, df_ss, df_ss_sim, df_mult, df_mult_sim), bounds=bounds, constraints=con1)
#    popt_all = minimize(utils.total_minimize, test, args=(df_tas, df_tas_sim, df_ss, df_ss_sim, df_mult, df_mult_sim))
    data_coef=popt_all.x
#    df_tas['error'] = df_tas['error'].apply(lambda x: x*10)
    df_tas['chi2'] = utils.running_chisq(data_coef, df_tas_sim, df_tas)#/int(df_tas.shape[0]+df_ss.shape[0]-len(data_coef))
    running_chi2 = []
    for i in range((df_tas.shape[0])):
        if i > 0:
            running_chi2.append(running_chi2[i-1]+df_tas['chi2'][i])
        else:
            running_chi2.append(df_tas['chi2'][i])
    df_tas['sumchi2'] = running_chi2
#    df_ss['error'] = df_ss['error'].apply(lambda x: x*10)
    df_ss['chi2'] = utils.running_chisq(data_coef, df_ss_sim, df_ss)#/int(df_tas.shape[0]+df_ss.shape[0]-len(data_coef))
    running_chi2 = []
    for i in range((df_ss.shape[0])):
        if i > 0:
            running_chi2.append(running_chi2[i-1]+df_ss['chi2'][i])
        else:
            running_chi2.append(df_ss['chi2'][i])
    df_ss['sumchi2'] = running_chi2
    df_mult['chi2'] = utils.running_chisq(data_coef, df_mult_sim, df_mult)#/int(df_tas.shape[0]+df_ss.shape[0]-len(data_coef))
    running_chi2 = []
    for i in range((df_mult.shape[0])):
        if i > 0:
            running_chi2.append(running_chi2[i-1]+df_mult['chi2'][i])
        else:
            running_chi2.append(df_mult['chi2'][i])
    df_mult['sumchi2'] = running_chi2
    df_tas['fit'] = np.sum(np.multiply(df_tas_sim,popt_all.x),axis=1)
    df_ss['fit'] = np.sum(np.multiply(df_ss_sim,popt_all.x),axis=1)
    df_mult['fit'] = np.sum(np.multiply(df_mult_sim,popt_all.x),axis=1)
    df_ss_m1['fit'] = np.sum(np.multiply(df_ss_m1_sim,popt_all.x),axis=1)
    utils.plot_all(df_tas, df_ss, df_mult, df_ss_m1)
    f = open("{}/{}".format(dt_string, "output.txt"),'w')
    f.write("{}".format(popt_all))
    f.write("\n")

    continuum = 0
    for i in range(len(energies)):
        print("{0} , {1:.2f}".format(energies[i],popt_all.x[i]*100))
        f.write("{0} , {1:.2f}\n".format(energies[i],popt_all.x[i]*100))
        if 'n' not in energies[i]:
            if 'm' in energies[i]:
                continue
            elif int(energies[i]) > config['cutoff']:
                continuum += popt_all.x[i]*100
    f.write("{0:.2f}\n".format(continuum))
    f.write(df_tas.to_string())
    f.write(df_ss.to_string())
    f.write(df_ss_m1.to_string())
    f.write(df_mult.to_string())
    print(df_ss.shape[0] + df_tas.shape[0] - len(energies))
    print(df_ss.shape[0] + df_tas.shape[0] - len(energies))
    reduce_factor = max(df_tas.shape[0]-test.shape[0],1)
    total_chi2 = np.sum(df_tas['chi2'])/reduce_factor
    print("tas chi2: {}, reduced chi2: {}".format(total_chi2*reduce_factor, total_chi2))
    f.write("tas chi2: {}, reduced chi2: {}\n".format(total_chi2*reduce_factor, total_chi2))
    total_chi2_tas = total_chi2
    reduce_factor = max(df_ss.shape[0]-test.shape[0],1)
    total_chi2_ss = np.sum(df_ss['chi2'])/reduce_factor
    f.write("ss chi2: {}, reduced chi2: {}\n".format(total_chi2_ss*reduce_factor, total_chi2_ss))
    print("ss chi2: {}, reduced chi2: {}".format(total_chi2_ss*reduce_factor, total_chi2_ss))
    reduce_factor = max(df_mult.shape[0]-test.shape[0],1)
    total_chi2_mult = np.sum(df_mult['chi2'])/reduce_factor
    f.write("mult chi2: {}, reduced chi2: {}\n".format(total_chi2_mult*reduce_factor, total_chi2_mult))
    print("mult chi2: {}, reduced chi2: {}".format(total_chi2_mult*reduce_factor, total_chi2_mult))
    reduce_factor = max(df_ss.shape[0]+df_tas.shape[0]-test.shape[0],1)
    total_chi2_tasss = (np.sum(df_tas['chi2'])+np.sum(df_ss['chi2']))/reduce_factor
    print("tas+ss chi2: {}, reduced chi2: {}".format(total_chi2_tasss*reduce_factor, total_chi2_tasss))
    f.write("tas+ss chi2: {}, reduced chi2: {}\n".format(total_chi2_tasss*reduce_factor, total_chi2_tasss))

    for plot in ['tas_fit_all_{}bins.png'.format(tas_bin_width),'ss_fit_hi_{}bins.png'.format(tas_bin_width),'ss_fit_lo_{}bins.png'.format(tas_bin_width),'ss_m1_fit_all_{}bins.png'.format(tas_bin_width),'mult_fit_{}bins.png'.format(tas_bin_width)]:
        shutil.move(plot, "{}/{}".format(dt_string, plot))
    shutil.copy('/mnt/analysis/e17028/RAINIER/{}/settings.h'.format(config['isotope_daughter']), "{}/{}".format(dt_string, "settings.h"))
    shutil.copy('/mnt/analysis/e17028/RAINIER/{}/z023_edit.datcopy'.format(config['isotope_daughter']), "{}/{}".format(dt_string, "z023_edit.datcopy"))
    f.close()
    return





if __name__ == "__main__":
    main()
