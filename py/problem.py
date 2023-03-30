#!/mnt/misc/sw/x86_64/all/anaconda/python3.7/bin/python

# make sure above path points to the version of python where you have pygmo installed 
# nscl servers
#!/mnt/misc/sw/x86_64/all/anaconda/python3.7/bin/python
# hpcc servers
#!/mnt/home/herman67/anaconda3/envs/pygmo/bin/python

# import commands
from utils import multi_minimize
from numpy import array, zeros, multiply, power, sum, append
from pandas import read_csv

# set important directories
PYGMO_DIR = './'
FOX_DIR = PYGMO_DIR 
OUTPUT_DIR = PYGMO_DIR 
SCRATCH_DIR = '/scratch/hermanse'

# make pygmo problem 
class optimizeRes:
    # define the output file to write solutions
    #   this is how we see all the attempted solutions
    #   and check how the algorithm performs
    # output file is defined when initiating the problem in optimize.py
    def __init__(self, dim,df_tas, df_tas_sim, df_ss, df_ss_sim, df_mult, df_mult_sim,df_ss_m1, df_ss_m1_sim, out="out"):
        self.dim = dim
        self.out = out
        self.df_tas = df_tas
        self.df_tas_sim = df_tas_sim
        self.df_ss = df_ss
        self.df_ss_sim = df_ss_sim
        self.df_mult = df_mult
        self.df_mult_sim = df_mult_sim
        self.df_ss_m1 = df_ss_m1
        self.df_ss_m1_sim = df_ss_m1_sim



    # fitness function that is called by each iteration
    #   x is between -1 and 1, but cosyrun expects a factor
    #   between 2^-1 and 2^1
    def fitness(self, x):
        # convert to scale factor
        pass_x = x #power(zeros(self.dim)+2.0,x)
        # run cosy
#        print(self.df_tas, self.df_tas_sim)
        resol = multi_minimize(pass_x, self.df_tas, self.df_tas_sim, self.df_ss, self.df_ss_sim, self.df_mult, self.df_mult_sim, self.df_ss_m1, self.df_ss_m1_sim)
        con = power((sum(pass_x)-1)*100,4)
        if con < 1:
            con=0
        resol = append(resol, con)

#        # write output to file
#        f = open(self.out,"a")
#        # write magnet values (as power of 2)
#        for i in range(len(x)):
#            f.write("{0},".format(pass_x[i]))
#        # write objective values (as ratio to nom)
#        for i in range(len(resol)):
#            if i == len(resol)-1:
#                f.write("{0}\n".format(resol[i]))
#            else: 
#                f.write("{0},".format(resol[i]))
#        f.close()
        # return objective values to evolve
        return [resol[0],resol[5]]

    # define number of objectives
    def get_nobj(self):
        return 2

    # define constraint (np.sum(coef) ==1)
#    def get_nec(self):
#        return 1

    # define bounds of x values
    #   i.e. from -1 to 1 (powers of 2)
    def get_bounds(self):
        qLower =zeros(self.dim)
        qUpper =zeros(self.dim) + 1.0
#        qLower[9] *= 1.5 
#        qUpper[9] *= 1.5 
#        return ([0]*41,[1]*41)
        return (qLower, qUpper)

    # define problem name
    def get_name(self):
        return "resolution optimizer"

    # print dimensions 
    def get_extra_info(self):
        return "\tDimensions: " + str(self.dim)




