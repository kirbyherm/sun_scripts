#!/mnt/simulations/secarml/soft/anaconda3/bin/python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys, os
from scipy.optimize import lsq_linear, minimize, Bounds, NonlinearConstraint
import json

# script, first, second = sys.argv
WORK_DIR = "../"
energies = []

def load_configs():

    f = open(WORK_DIR+'config.json')
    config = json.load(f)
#    print(config)
    f.close()

    return config

config = load_configs()
GEANT_DIR = "/mnt/analysis/e17028/SuNGEANT/geant/sunNewDetector_{}/rootFiles/".format(config['isotope_daughter'])
NGEANT_DIR = "/mnt/analysis/e17028/SuNGEANT/geant/sunNewDetector_{}/rootFiles/".format(config['isotope_ndaughter'])

def monte_carlo(df,rands,i):
    df_rand = df.copy(deep=True)
    df_rand['content'] = rands
    return df_rand

def running_chisq(coef, obs, exp):
    sim = np.sum(np.multiply(obs, coef),axis=1)
    chi2 = (np.sum((sim - exp['content']) ** 2 / np.maximum(np.zeros(exp.shape[0])+1,(exp['error']) ** 2)))
    return (sim - exp['content']) ** 2 / (np.maximum(np.zeros(exp.shape[0])+1,exp['error']) ** 2)

def chisq(coef, obs, exp):
    sim = np.sum(np.multiply(obs, coef),axis=1)
    chi2 = (np.sum((sim - exp['content']) ** 2 / (np.maximum(np.zeros(exp.shape[0])+1,exp['error'])**2)))
    return chi2 


def multi_minimize(coef, df_tas, df_tas_sim, df_ss, df_ss_sim, df_mult, df_mult_sim, df_ss_m1, df_ss_m1_sim):
    test = coef
    df_all = pd.concat([df_tas,df_ss,df_mult])
    df_all_sim = pd.concat([df_tas_sim,df_ss_sim,df_mult_sim])
    popt_tas = chisq( test, df_tas_sim, df_tas)
    popt_ss = chisq(test, df_ss_sim, df_ss)
    popt_mult = chisq(test, df_mult_sim, df_mult)
    popt_ss_m1 = chisq(test, df_ss_m1_sim, df_ss_m1)
    popt_all = np.sum([popt_tas, popt_ss, popt_mult, popt_ss_m1])
    reduce_factor = max(df_tas.shape[0]-test.shape[0],1)
    total_chi2 = np.sum([popt_tas])/reduce_factor
    total_chi2_tas = total_chi2
    reduce_factor = max(df_ss.shape[0]-test.shape[0],1)
    total_chi2 = np.sum([popt_ss])/reduce_factor
    reduce_factor = max(df_ss.shape[0]+df_tas.shape[0]-test.shape[0],1)
    total_chi2 = np.sum([popt_tas, popt_ss])/reduce_factor
    reduce_factor = max(df_ss.shape[0]+df_tas.shape[0]+df_mult.shape[0]-test.shape[0],1)
    total_chi2 = np.sum([popt_tas, popt_ss, popt_mult])/reduce_factor
    return [popt_tas, popt_ss, popt_mult, popt_ss_m1, popt_all]
     

     
def total_minimize(coef, df_tas, df_tas_sim, df_ss, df_ss_sim, df_mult, df_mult_sim, df_ss_m1, df_ss_m1_sim, fit_param):
    test = coef 
    popt_tas = chisq( test, df_tas_sim, df_tas)
    popt_ss = chisq(test, df_ss_sim, df_ss)
    popt_mult = chisq(test, df_mult_sim, df_mult)
    popt_ss_m1 = chisq(test, df_ss_m1_sim, df_ss_m1)
    reduce_factor = max(df_tas.shape[0]-test.shape[0],1)
    total_chi2 = np.sum([popt_tas])/reduce_factor
#    print("tas chi2: {}, reduced chi2: {}".format(total_chi2*reduce_factor, total_chi2))
    total_chi2_tas = total_chi2
    reduce_factor = max(df_ss.shape[0]-test.shape[0],1)
    total_chi2_ss = np.sum([popt_ss])/reduce_factor
#    print("ss chi2: {}, reduced chi2: {}".format(total_chi2_ss*reduce_factor, total_chi2_ss))
    reduce_factor = max(df_mult.shape[0]-test.shape[0],1)
    total_chi2_mult = np.sum([popt_mult])/reduce_factor
#    print("mult chi2: {}, reduced chi2: {}".format(total_chi2_mult*reduce_factor, total_chi2_mult))
    reduce_factor = max(df_ss.shape[0]+df_tas.shape[0]-test.shape[0],1)
    total_chi2_tasss = np.sum([popt_tas, popt_ss])/reduce_factor
#    print("tas+ss chi2: {}, reduced chi2: {}".format(total_chi2_tasss*reduce_factor, total_chi2_tasss))
    reduce_factor = max(df_ss.shape[0]+df_tas.shape[0]+df_mult.shape[0]-test.shape[0],1)
    total_chi2 = np.sum([popt_tas, popt_ss, popt_mult])/reduce_factor
    reduce_factor = max(df_ss.shape[0]+df_tas.shape[0]+df_mult.shape[0]+df_ss_m1.shape[0]-test.shape[0],1)
    total_chi2_m1 = np.sum([popt_tas, popt_ss, popt_mult,popt_ss_m1])/reduce_factor
#    print("total chi2: {}, reduced chi2: {}".format(total_chi2*reduce_factor, total_chi2))
    if fit_param == 'tas':
        total_chi2 = total_chi2_tas
    elif fit_param == 'ss':
        total_chi2 = total_chi2_ss
    elif fit_param == 'tasss':
        total_chi2 = total_chi2_tasss
    elif fit_param == 'm1':
        total_chi2 = total_chi2_m1
    return total_chi2

def plot_all(df_tas, df_ss, df_mult, df_ss_m1, figname="", config=load_configs()):

    tas_bin_width = config['tas_bin_width']
    tas_lower_cut = config['tas_lower_limit']
    tas_upper_cut = config['tas_upper_limit']
    ss_bin_width = config['ss_bin_width']
    ss_lower_cut = config['ss_lower_limit']
    ss_upper_cut = config['ss_upper_limit']
    energies = config['levels']
    isotope = config['isotope']

    print(df_tas, tas_bin_width, tas_lower_cut, tas_upper_cut)

    # Draw tas
    lo_plot, hi_plot = int(tas_lower_cut/tas_bin_width), int(tas_upper_cut/tas_bin_width)
    fig, ax = plt.subplots()
#    plt.errorbar((np.arange(df_tas.shape[0])*tas_bin_width)[lo_plot:hi_plot],df_tas['content'][lo_plot:hi_plot],yerr=df_tas['error'][lo_plot:hi_plot],label='obs',capsize=2,elinewidth=0.5,barsabove=True,errorevery=1)
#    plt.plot((np.arange(df_tas.shape[0])*tas_bin_width)[lo_plot:hi_plot],df_tas['fit'][lo_plot:hi_plot],label='fit')
    plt.errorbar((df_tas['bin'])[lo_plot:hi_plot],df_tas['content'][lo_plot:hi_plot],yerr=df_tas['error'][lo_plot:hi_plot],label='obs',capsize=2,elinewidth=0.5,barsabove=True,errorevery=1)
    plt.plot((df_tas['bin'])[lo_plot:hi_plot],df_tas['fit'][lo_plot:hi_plot],label='fit')
    plt.xlabel("keV")
    plt.ylabel("counts/{}keV".format(tas_bin_width))
    xmin, xmax, ymin, ymax = plt.axis()
    fig.set_size_inches(15,5)
    ax2 = ax.twinx()
    ax2.plot((df_tas['bin'])[lo_plot:hi_plot],df_tas['sumchi2'][lo_plot:hi_plot], 'r', linestyle='dashed', label='red chi2')
    ax2.set_ylabel('cumulative chi2')
    fig.legend()
    plt.savefig('tas_fit_all_{}bins{}.png'.format(tas_bin_width,figname),dpi=200)
    plt.cla()
#    ax2 = ax.twinx()
#    ax2.plot((np.arange(df_tas.shape[0])*tas_bin_width)[:lo_plot],df_tas['sumchi2'][:lo_plot], 'r', linestyle='dashed', label='red chi2')
#    ax2.set_ylabel('cumulative chi2')

    # Draw ss
    plt.cla()
    ss_lower_cut = 0
    ss_upper_cut = 2100
    lo_plot, hi_plot = int(ss_lower_cut/ss_bin_width), int(ss_upper_cut/ss_bin_width)
    print(lo_plot, hi_plot, df_ss)
    fig, ax = plt.subplots()
#    plt.errorbar((np.arange(df_ss.shape[0])*ss_bin_width)[lo_plot:hi_plot],df_ss['content'][lo_plot:hi_plot],yerr=df_ss['error'][lo_plot:hi_plot],label='obs',capsize=2,elinewidth=0.5,barsabove=True,errorevery=1)
#    plt.plot((np.arange(df_ss.shape[0])*ss_bin_width)[lo_plot:hi_plot],df_ss['fit'][lo_plot:hi_plot],label='fit')
    plt.errorbar(df_ss['bin'][lo_plot:hi_plot],df_ss['content'][lo_plot:hi_plot],yerr=df_ss['error'][lo_plot:hi_plot],label='obs',capsize=2,elinewidth=0.5,barsabove=True,errorevery=1)
#    plt.plot((np.arange(df_ss.shape[0])*ss_bin_width)[lo_plot:hi_plot],df_ss['fit'][lo_plot:hi_plot],label='fit')
    plt.plot(df_ss['bin'][lo_plot:hi_plot],df_ss['fit'][lo_plot:hi_plot],label='fit')
    plt.legend()
    xmin, xmax, ymin, ymax = plt.axis()
    plt.xlabel("keV")
    plt.ylabel("counts/{}keV".format(ss_bin_width))
    plt.yscale('log')
#    plt.ylim([1,ymax])
    ax2 = ax.twinx()
    ax2.plot(df_ss['bin'][lo_plot:hi_plot],df_ss['sumchi2'][lo_plot:hi_plot], 'r', linestyle='dashed', label='red chi2')
    ax2.set_ylabel('cumulative chi2')
    fig.legend()
    plt.savefig('ss_fit_lo_{}bins{}.png'.format(tas_bin_width,figname),dpi=200)

    # Draw ss
    plt.cla()
    ss_lower_cut = 2100
    ss_upper_cut = 12800
    lo_plot, hi_plot = int(ss_lower_cut/ss_bin_width), int(ss_upper_cut/ss_bin_width)
    print(lo_plot, hi_plot, df_ss)
    fig, ax = plt.subplots()
#    plt.errorbar((np.arange(df_ss.shape[0])*ss_bin_width)[lo_plot:hi_plot],df_ss['content'][lo_plot:hi_plot],yerr=df_ss['error'][lo_plot:hi_plot],label='obs',capsize=2,elinewidth=0.5,barsabove=True,errorevery=1)
#    plt.plot((np.arange(df_ss.shape[0])*ss_bin_width)[lo_plot:hi_plot],df_ss['fit'][lo_plot:hi_plot],label='fit')
    plt.errorbar(df_ss['bin'][lo_plot:hi_plot],df_ss['content'][lo_plot:hi_plot],yerr=df_ss['error'][lo_plot:hi_plot],label='obs',capsize=2,elinewidth=0.5,barsabove=True,errorevery=1)
#    plt.plot((np.arange(df_ss.shape[0])*ss_bin_width)[lo_plot:hi_plot],df_ss['fit'][lo_plot:hi_plot],label='fit')
    plt.plot(df_ss['bin'][lo_plot:hi_plot],df_ss['fit'][lo_plot:hi_plot],label='fit')
    plt.legend()
    xmin, xmax, ymin, ymax = plt.axis()
    plt.xlabel("keV")
    plt.ylabel("counts/{}keV".format(ss_bin_width))
    plt.yscale('log')
    plt.ylim([1,ymax])
    ax2 = ax.twinx()
    ax2.plot(df_ss['bin'][lo_plot:hi_plot],df_ss['sumchi2'][lo_plot:hi_plot], 'r', linestyle='dashed', label='red chi2')
    ax2.set_ylabel('cumulative chi2')
    fig.legend()
    plt.savefig('ss_fit_hi_{}bins{}.png'.format(tas_bin_width,figname),dpi=200)


    # Draw mult
    plt.clf()
#    print(df_mult, np.arange(10))
    print(df_mult)
    fig, ax = plt.subplots()
    plt.errorbar((np.arange(df_mult.shape[0])),df_mult['content'],yerr=df_mult['error'],label='obs',capsize=2,elinewidth=0.5,barsabove=True,errorevery=1)
    plt.plot((np.arange(df_mult.shape[0])),df_mult['fit'],label='fit')
    plt.legend()
    plt.xlabel("multiplicity")
    plt.ylabel("counts/mult")
    ax2 = ax.twinx()
    ax2.plot((np.arange(df_mult.shape[0])),df_mult['sumchi2'], 'r', linestyle='dashed', label='red chi2')
    ax2.set_ylabel('cumulative chi2')
    plt.savefig('mult_fit_{}bins{}.png'.format(tas_bin_width,figname),dpi=200)

    # Draw ss_m1
    plt.cla()
    ss_lower_cut = config['ss_lower_limit']
    ss_upper_cut = config['ss_upper_limit']
    lo_plot, hi_plot = int(ss_lower_cut/ss_bin_width), int(ss_upper_cut/ss_bin_width)
    fig, ax = plt.subplots()
#    plt.errorbar((np.arange(df_ss.shape[0])*ss_bin_width)[lo_plot:hi_plot],df_ss['content'][lo_plot:hi_plot],yerr=df_ss['error'][lo_plot:hi_plot],label='obs',capsize=2,elinewidth=0.5,barsabove=True,errorevery=1)
    plt.errorbar(df_ss_m1['bin'][lo_plot:hi_plot],df_ss_m1['content'][lo_plot:hi_plot],yerr=df_ss_m1['error'][lo_plot:hi_plot],label='obs',capsize=2,elinewidth=0.5,barsabove=True,errorevery=1)
#    plt.plot((np.arange(df_ss_m1.shape[0])*ss_m1_bin_width)[lo_plot:hi_plot],df_ss_m1['fit'][lo_plot:hi_plot],label='fit')
    plt.plot(df_ss_m1['bin'][lo_plot:hi_plot],df_ss_m1['fit'][lo_plot:hi_plot],label='fit')
    plt.legend()
    xmin, xmax, ymin, ymax = plt.axis()
    plt.xlabel("keV")
    plt.ylabel("counts/{}keV".format(ss_bin_width))
    plt.yscale('log')
    plt.ylim(1,4000)
    plt.xlim(0,10000)
#    ax2 = ax.twinx()
#    ax2.plot((np.arange(df_ss.shape[0])*ss_bin_width)[lo_plot:hi_plot],df_ss['sumchi2'][lo_plot:hi_plot], 'r', linestyle='dashed', label='red chi2')
#    ax2.set_ylabel('cumulative chi2')
    plt.legend()
    fig.set_size_inches(15,5)
    plt.savefig('ss_m1_fit_all_{}bins{}.png'.format(tas_bin_width,figname),dpi=200)
    plt.cla()

    return

def rebin_data(df, config=load_configs()):
    tas_bin_width = config['tas_bin_width']
    tas_lower_cut = config['tas_lower_limit']
    tas_upper_cut = config['tas_upper_limit']
    ss_bin_width = config['ss_bin_width']
    ss_lower_cut = config['ss_lower_limit']
    ss_upper_cut = config['ss_upper_limit']
    energies = config['levels']
    isotope = config['isotope']
    cutoff = config['cutoff']

    energies_int = []
    cont_bins = 0
    n_levels = 0
    for i in energies:
        if "n" in i:
            n_levels+=1
            continue
        if cutoff > tas_upper_cut:
            continue
        energies_int.append(int(i))
        if int(i) >= cutoff:
            cont_bins += 1
    discrete_bin_width = tas_bin_width
    cutoff_bin = int(min(cutoff,tas_upper_cut)/discrete_bin_width)
    running_bin_values = [0] 
    bin_widths = np.zeros(cutoff_bin + cont_bins) + discrete_bin_width 
    for i in range(len(bin_widths)):
        if i > cutoff_bin and i < len(bin_widths)-1:
            bin_widths[i] = energies_int[i-cutoff_bin+len(energies)-cont_bins-n_levels]-energies_int[i-cutoff_bin+len(energies)-cont_bins-1-n_levels]
        elif i > cutoff_bin:
            bin_widths[i] = 400
        if i > 0:
            running_bin_values.append(bin_widths[i]+running_bin_values[i-1])
#    print(running_bin_values)
    df_new = pd.DataFrame()
    bins = np.zeros(len(bin_widths),dtype=int) 
    contents = np.zeros(len(bin_widths),dtype=float) 
    errors = np.zeros(len(bin_widths),dtype=float) 
    bin_width_counter = 0
    bin_counter = 0
#    print(df.at[0,'content'])
    for i in range(df.shape[0]):
        if bin_width_counter >= bin_widths[bin_counter]:
            errors[bin_counter] = np.sqrt(errors[bin_counter])
            bin_counter += 1
            if bin_counter >= len(bin_widths):
                break
            bins[bin_counter] = bin_width_counter 
            if bin_counter>1:
                bins[bin_counter] += bins[bin_counter-1]
            bin_width_counter = 0
        contents[bin_counter] += df.at[i,'content']*discrete_bin_width/bin_widths[bin_counter]
        errors[bin_counter] += df.at[i,'error'] * df.at[i,'error']
        bin_width_counter += discrete_bin_width
    df_new['content'] = contents
    df_new['bin'] = running_bin_values
    df_new['error'] = errors
    return df_new

def rebin_sims(df, config):
    tas_bin_width = config['tas_bin_width']
    tas_lower_cut = config['tas_lower_limit']
    tas_upper_cut = config['tas_upper_limit']
    ss_bin_width = config['ss_bin_width']
    ss_lower_cut = config['ss_lower_limit']
    ss_upper_cut = config['ss_upper_limit']
    energies = config['levels']
    isotope = config['isotope']
    cutoff = config['cutoff']

    energies_int = []
    cont_bins = 0
    n_levels = 0
    for i in energies:
        if "n" in i:
            n_levels+=1
            continue
        if cutoff > tas_upper_cut:
            continue
        energies_int.append(int(i))
        if int(i) >= cutoff:
            cont_bins += 1
    discrete_bin_width = tas_bin_width
    running_bin_values = [0] 
    cutoff_bin = int(min(cutoff,tas_upper_cut)/discrete_bin_width)
    bin_widths = np.zeros(cutoff_bin + cont_bins) + discrete_bin_width 
    for i in range(len(bin_widths)):
        if i > cutoff_bin and i < len(bin_widths)-1:
            bin_widths[i] = energies_int[i-cutoff_bin+len(energies)-cont_bins-n_levels]-energies_int[i-cutoff_bin+len(energies)-cont_bins-1-n_levels]
        elif i > cutoff_bin:
            bin_widths[i] = 400
        if i > 0:
            running_bin_values.append(bin_widths[i]+running_bin_values[i-1])
#    bin_widths = np.zeros(cont_bins+cutoff_bin)
#    running_bin_values = [0]
#    for i in range(len(bin_widths)):
#        bin_widths[i] = discrete_bin_width 
#        if i >cutoff_bin:
#            bin_widths[i] = 100
#        if i >cutoff_bin+3:
#            bin_widths[i] = 150
#        if i >cutoff_bin+9:
#            bin_widths[i] = 200
#        if i >cutoff_bin+25:
#            bin_widths[i] = 300
#        if i >cutoff_bin+36:
#            bin_widths[i] = 400
#        if i > 0:
#            running_bin_values.append(bin_widths[i]+running_bin_values[i-1])
#    print(running_bin_values)
    df_new = pd.DataFrame()
    for j in range(df.shape[1]):
        contents = np.zeros(len(bin_widths),dtype=float) 
        bin_width_counter = 0
        bin_counter = 0
        for i in range(df.shape[0]):
            if bin_width_counter >= bin_widths[bin_counter]:
                bin_counter += 1
                bin_width_counter = 0
            if bin_counter >= len(bin_widths):
                break
            contents[bin_counter] += df.at[i,'h{}'.format(energies[j])]*discrete_bin_width/bin_widths[bin_counter]
            bin_width_counter += discrete_bin_width
        df_new['h{}'.format(energies[j])] = contents
    return df_new


def load_tas_sims(scale_factor, config=load_configs(), samp=0):

    tas_bin_width = config['tas_bin_width']
    energies = config['levels']
    upper_cut = int(config['tas_upper_limit']/tas_bin_width)

    df = pd.DataFrame()
    for i in energies:
        data = None
        if samp > 0:
            if 'n' in i:
                data = np.loadtxt(NGEANT_DIR+'{}/tas_new/output_{}.csv'.format(samp,i),delimiter=',')
            else:
                data = np.loadtxt(GEANT_DIR+'{}/tas_new/output_{}.csv'.format(samp,i),delimiter=',')
        if 'n' in i:
            data = np.loadtxt(NGEANT_DIR+'tas_new/output_{}.csv'.format(i),delimiter=',')
        else:
            data = np.loadtxt(GEANT_DIR+'tas_new/output_{}.csv'.format(i),delimiter=',')
        content = data[:,1]/np.sum(data[:,1])*scale_factor
        df['h{}'.format(i)] = content
    df = df.groupby(df.index // tas_bin_width).sum()
    return df.iloc[:upper_cut,:]

def load_mult_sims(scale_factor, config=load_configs(), samp=0):
    energies=config['levels']
    df = pd.DataFrame()
    for i in energies:
        data = None
        if samp > 0:
            if 'n' in i:
                data = np.loadtxt(NGEANT_DIR+'{}/mult_new/output_{}.csv'.format(samp,i),delimiter=',')
            else:
                data = np.loadtxt(GEANT_DIR+'{}/mult_new/output_{}.csv'.format(samp,i),delimiter=',')
        if 'n' in i:
            data = np.loadtxt(NGEANT_DIR+'mult_new/output_{}.csv'.format(i),delimiter=',')
        else:
            data = np.loadtxt(GEANT_DIR+'mult_new/output_{}.csv'.format(i),delimiter=',')
        content = data[:,1]/np.sum(data[:,1])*scale_factor
        content = content[:9]
#        print(data)
        df['h{}'.format(i)] = content
    return df.iloc[:9,:]

def load_ss_sims(scale_factor, config=load_configs(), samp=0):
    df = pd.DataFrame()

    ss_bin_width = config['ss_bin_width']
    energies = config['levels']
    upper_cut = int(config['ss_upper_limit']/ss_bin_width)
    for i in energies:
        data = None
        if samp > 0:
            if 'n' in i:
                data = np.loadtxt(NGEANT_DIR+'{}/ss_new/Enew{}.csv'.format(samp,i),delimiter=',')
            else:
                data = np.loadtxt(GEANT_DIR+'{}/ss_new/Enew{}.csv'.format(samp,i),delimiter=',')
        if 'n' in i:
            data = np.loadtxt(NGEANT_DIR+'ss_new/Enew{}.csv'.format(i),delimiter=',')
        else:
            data = np.loadtxt(GEANT_DIR+'ss_new/Enew{}.csv'.format(i),delimiter=',')
        content = data[:,1]/np.sum(data[:,1])*scale_factor
        df['h{}'.format(i)] = content
    df = df.groupby(df.index // ss_bin_width).sum()
    return df.iloc[:upper_cut,:]

def load_ss_m1_sims(scale_factor, config=load_configs(), samp=0):
    df = pd.DataFrame()
    ss_bin_width = config['ss_bin_width']
    energies = config['levels']
    upper_cut = int(config['ss_upper_limit']/ss_bin_width)
    for i in energies:
        data = None
        if samp > 0:
            if 'n' in i:
                data = np.loadtxt(NGEANT_DIR+'{}/ss_new_m1/Enew{}.csv'.format(samp,i),delimiter=',')
            else:
                data = np.loadtxt(GEANT_DIR+'{}/ss_new_m1/Enew{}.csv'.format(samp,i),delimiter=',')
        if 'n' in i:
            data = np.loadtxt(NGEANT_DIR+'ss_new_m1/Enew{}.csv'.format(i),delimiter=',')
        else:
            data = np.loadtxt(GEANT_DIR+'ss_new_m1/Enew{}.csv'.format(i),delimiter=',')
#        if i < 2476:
#            data = np.loadtxt(WORK_DIR+'ss_new/Enew{}.csv'.format(i),delimiter=',')
#        else:
#            data = np.loadtxt(WORK_DIR+'ss_csv/Enew{}.csv'.format(i),delimiter=',')
        content = data[:,1]/np.sum(data[:,1])*scale_factor
        df['h{}'.format(i)] = content
    df = df.groupby(df.index // ss_bin_width).sum()
    return df.iloc[:upper_cut,:]


def normalize_sims(df, scale_factor, config=load_configs()):
    energies = config['levels']
    new_df = pd.DataFrame()
    for i in energies:
        content = df['h{}'.format(i)]/np.sum(df['h{}'.format(i)])*scale_factor
        df['h{}'.format(i)] = content
    return df

def combine_continuum(df, cont):
    df["continuum"] = 0
    cont_sum = 0
    for i in cont:
        df["continuum"] += df[i] * cont[i]
        cont_sum += cont[i]
        df = df.drop(i,axis=1)
    df['continuum'] = df['continuum']/cont_sum
    print(df)
    return df

def load_all(samp=0, temp=""):

    config = load_configs()
    tas_bin_width = config['tas_bin_width']
    tas_lower_cut = config['tas_lower_limit']
    tas_upper_cut = config['tas_upper_limit']
    ss_bin_width = config['ss_bin_width']
    ss_lower_cut = config['ss_lower_limit']
    ss_upper_cut = config['ss_upper_limit']
    energies = config['levels']
    isotope = config['isotope']
    
    df_tas = pd.read_csv(WORK_DIR+"{}_tas_obs{}.csv".format(isotope, temp),names=["bin","content","error"],dtype=(float,float))
    df_tas = rebin_data(df_tas)
    tas_lower_cut = int(tas_lower_cut/tas_bin_width)
    tas_upper_cut = int(tas_upper_cut/tas_bin_width)
    df_tas = df_tas.iloc[tas_lower_cut:tas_upper_cut,:]
#    df_tas['error'] = df_tas['error'].apply(lambda x: x*0.10+5)
    df_tas = df_tas.reset_index(drop=True)
    df_ss = pd.read_csv(WORK_DIR+"{}_ss_obs{}.csv".format(isotope, temp),names=["bin","content","error"],dtype=(float,float))
    df_ss_m1 = pd.read_csv(WORK_DIR+"mult_1/{}_ss_obs{}.csv".format(isotope,temp),names=["bin","content","error"],dtype=(float,float))
    ss_lower_cut = int(ss_lower_cut/ss_bin_width)
    ss_upper_cut = int(ss_upper_cut/ss_bin_width)
    df_ss = df_ss.iloc[ss_lower_cut:ss_upper_cut,:]
    df_ss = rebin_data(df_ss)
    df_ss = df_ss.reset_index(drop=True)
    df_ss_m1 = rebin_data(df_ss_m1)
    df_ss_m1 = df_ss_m1.reset_index(drop=True)
#    print(ss_upper_cut)
#    df_ss['error'] = df_ss['error'].apply(lambda x: x*0.50+20)
    df_mult = pd.read_csv(WORK_DIR+"{}_mult_obs{}.csv".format(isotope, temp),names=["bin","content","error"],dtype=(float,float))
    df_mult = df_mult.iloc[:9,:]
#    df_mult['error'] = df_mult['error'].apply(lambda x: x*10)
    scale_factor = np.sum(df_tas['content'])
    df_tas_sim = load_tas_sims(scale_factor, config, samp)
    df_tas_sim = rebin_sims(df_tas_sim, config)
    df_tas_sim = df_tas_sim.iloc[tas_lower_cut:tas_upper_cut,:]
    df_tas_sim = normalize_sims(df_tas_sim, scale_factor, config)
    df_tas_sim = df_tas_sim.reset_index(drop=True)
    scale_factor = np.sum(df_ss['content'])
    df_ss_sim = load_ss_sims(scale_factor, config, samp)
    df_ss_sim = rebin_sims(df_ss_sim, config)
#    df_ss_sim = df_ss_sim.iloc[ss_lower_cut:ss_upper_cut,:]
    df_ss_sim = normalize_sims(df_ss_sim, scale_factor, config)
    df_ss_sim = df_ss_sim.reset_index(drop=True)
    scale_factor = np.sum(df_ss_m1['content'])
    df_ss_m1_sim = load_ss_m1_sims(scale_factor, config ,samp)
    df_ss_m1_sim = rebin_sims(df_ss_m1_sim, config)
#    df_ss_sim = df_ss_sim.iloc[ss_lower_cut:ss_upper_cut,:]
    df_ss_m1_sim = normalize_sims(df_ss_m1_sim, scale_factor, config)
    df_ss_m1_sim = df_ss_m1_sim.reset_index(drop=True)
    scale_factor = np.sum(df_mult['content'])
    df_mult_sim = load_mult_sims(scale_factor, config, samp)

    return df_tas, df_tas_sim, df_ss, df_ss_sim, df_mult, df_mult_sim, df_ss_m1, df_ss_m1_sim


def run(coef=np.zeros(44), nameout="test"):

    config = load_configs()
    tas_bin_width = config['tas_bin_width']
    tas_lower_cut = config['tas_lower_limit']
    tas_upper_cut = config['tas_upper_limit']
    ss_bin_width = config['ss_bin_width']
    ss_lower_cut = config['ss_lower_limit']
    ss_upper_cut = config['ss_upper_limit']
    energies = config['levels']
    isotope = config['isotope']
    

    df_tas = pd.read_csv(WORK_DIR+"{}_tas_obs.csv".format(isotope),names=["bin","content","error"],dtype=(float,float))
#    df_tas['content'].mask(df_tas['content'] < 0,0,inplace=True)
    df_ss = pd.read_csv(WORK_DIR+"{}_ss_obs.csv".format(isotope),names=["bin","content","error"],dtype=(float,float))
#    df_ss['content'].mask(df_ss['content'] < 0,0,inplace=True)
    df_mult = pd.read_csv(WORK_DIR+"{}_mult_obs.csv".format(isotope),names=["bin","content","error"],dtype=(float,float))
#    df_mult['content'].mask(df_mult['content'] < 0,0,inplace=True)
    df_mult = df_mult.iloc[:10,:]
#    print(np.sum(df_tas),np.sum(df_ss), np.sum(df_mult))
    scale_factor = np.sum(df_tas['content'])
    df_tas_sim = load_tas_sims(scale_factor,config)
    scale_factor = np.sum(df_ss['content'])
    df_ss_sim = load_ss_sims(scale_factor, config)
    scale_factor = np.sum(df_mult['content'])
    df_mult_sim = load_mult_sims(scale_factor, config)

    test = [0.        , 0.        , 0.00242834, 0.24587247, 0.26753184,
       0.03340626, 0.02989332, 0.16553433, 0.01720398, 0.00967029,
       0.        , 0.00282008, 0.0465006 , 0.04302853, 0.        ,
       0.01044767, 0.01423525, 0.02762151, 0.        , 0.01414714,
       0.01830467, 0.01514356, 0.        , 0.        , 0.00389088,
       0.00410051, 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        ]
#    test = np.zeros(44)

    #resol = total_minimize(test,df_tas,df_tas_sim,df_ss, df_ss_sim, df_mult, df_mult_sim)
#    X = np.array([1]) 
#    for column in df_ss_sim.columns:
#        if len(X) == 1 :
#            X = np.array(np.vstack((np.vstack((np.array(df_tas_sim[column])[:1000].reshape(-1,1),np.array(df_ss_sim[column])[:2500].reshape(-1,1))),np.array(df_mult_sim[column]).reshape(-1,1)))).flatten().reshape(-1,1)
#
#            print(column, X.shape)
#        else:
##            temp_vec = np.array(np.vstack((np.array(df_tas_sim[column])[:1000].reshape(-1,1),np.array(df_ss_sim[column])[:2500].reshape(-1,1)))).flatten().reshape(-1,1)
#            temp_vec = np.array(np.vstack((np.vstack((np.array(df_tas_sim[column])[:1000].reshape(-1,1),np.array(df_ss_sim[column])[:2500].reshape(-1,1))),np.array(df_mult_sim[column]).reshape(-1,1)))).flatten().reshape(-1,1)
#            print(column, X.shape, temp_vec.shape)
#            X = np.column_stack((X,temp_vec))        
#            print(column, X.shape)
##    X = np.array([df_tas_sim.iloc[:900,:],df_ss_sim.iloc[:900,:]])
##    y = pd.DataFrame({"tas_obs":df_tas['content'],"ss_obs":df_ss['content']}).iloc[:900,:]
#    y = np.vstack((np.vstack((np.array(df_tas['content'],dtype=np.float)[:1000].reshape(-1,1),np.array(df_ss['content'],dtype=np.float)[:2500].reshape(-1,1))),np.array(df_mult['content'],dtype=np.float).reshape(-1,1))).flatten()

#    df_all = df_all.append(df_ss)
#    df_all = df_all.append(df_mult)
#    print(df_all)
#    df_all = pd.concat([df_tas,df_ss,df_mult])
#    df_all_sim = pd.concat([df_tas_sim,df_ss_sim,df_mult_sim])
##    print(df_all_sim)
#    test = np.zeros(df_ss_sim.shape[1])
#    ensdf_lvl = np.array([23.8, 5, 1.1, 26, 26.6,0,0,0, 17.3])/100.0
#    for i in range(len(ensdf_lvl)):
#        test[i] += ensdf_lvl[i]
#    bounds = Bounds(0,1) 
#    con = lambda x: np.sum(x)-1
#    nlc = NonlinearConstraint(con, 0.999, 1.001)
#    con1 = ({'type': 'eq', 'fun': con})
#    popt_all = minimize(total_minimize, test, args=(df_tas, df_tas_sim, df_ss, df_ss_sim, df_mult, df_mult_sim), bounds=bounds, constraints=con1)
#    
#    print(popt_all)
##    upper_bound = np.zeros(df_ss_sim.shape[1])+np.inf
##    print(chisq(test, df_tas_sim, df_tas))
##    popt_tas = minimize(chisq, test, args=(df_tas_sim, df_tas), bounds=bounds)
##    print(popt_tas)
##    print(np.multiply(df_tas_sim,popt_tas.x))
##    popt_ss = minimize(chisq, test, args=(df_ss_sim, df_ss), bounds=bounds)
##    print(popt_ss)
##    popt_mult = minimize(chisq, test, args=(df_mult_sim, df_mult), bounds=bounds)
##    print(popt_mult)
##    popt_all = minimize(chisq, test, args=(df_all_sim, df_all), bounds=bounds)
##    print(popt_all)
##    print(chisq(popt_all.x, df_tas_sim, df_tas))
    
    df_tas['fit'] = np.sum(np.multiply(df_tas_sim,coef),axis=1)
    df_ss['fit'] = np.sum(np.multiply(df_ss_sim,coef),axis=1)
    df_mult['fit'] = np.sum(np.multiply(df_mult_sim,coef),axis=1)
    plot_all(df_tas, df_ss, df_mult, nameout)

    return coef

def get_best_point(fit_param='all'):

    config = load_configs()
    energies = config['levels']
    df_tas, df_tas_sim, df_ss, df_ss_sim, df_mult, df_mult_sim, df_ss_m1, df_ss_m1_sim = load_all()
    test = np.zeros(len(energies))
    bounds = Bounds(0,1) 
    con = lambda x: np.sum(x)-1
    nlc = NonlinearConstraint(con, 0.999, 1.001)
    con1 = ({'type': 'eq', 'fun': con})
    popt_all = minimize(total_minimize, test, args=(df_tas, df_tas_sim, df_ss, df_ss_sim, df_mult, df_mult_sim, df_ss_m1, df_ss_m1_sim, fit_param), bounds=bounds, constraints=con1)
    data_coef=popt_all.x
    return data_coef


if __name__ == "__main__":
    run()
