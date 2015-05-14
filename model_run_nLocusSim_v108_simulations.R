###########################################################
# Run model simulations for melanistic model              #
# 20 Sept. 2014                                           #
# A. A. Comeault                                          #
# simulation program and model name: "nLocusSim_v108"     #
# model writen by S. Flaxman                              #
###########################################################
#
#
# note that for this script to work properly,
# the nLocusSim_v108 program must be in your path.
#
#

rm(list=ls())

############################################################
# Function to run the "nLocusSum_v108" simulation program: #
############################################################
run_nLocusSim <- function(base_dir, prefix, a, A, b, B, c, C, d, D, f, F, g, G, h, H, L, m, M, N, o, p, P, r, R, s, S, t, T, w, Z)
{ 
  # define the directory to write ouput to:
  dir <- paste(base_dir, "/", prefix, sep="")
  
  # make and move into the appropriate directory for writing files:
  system(paste("rm -rf", dir))
  system(paste("mkdir", dir))
  setwd(dir)
  
  # run program:
  system(paste("nLocusSim_v108", 
               "-a", a, "-A", A, "-b", b, "-B", B, "-c", c, "-C", C, "-d", d, "-D", D, "-f", f, "-F", F, "-g", g, "-G", G, "-h", h, "-H", H, "-L", L, "-m", m, "-M", M, "-N", N, "-o", o, "-p", p, "-P", P, "-r", r, "-R", R, "-s", s, "-S", S, "-t", t, "-T", T, "-w", w, "-Z", Z))
  
  setwd(base_dir)  
}


#############################
# Set up default parameters #
#############################
# note that is important that the output for runs run with w = 0 differs from those
# run with w = 1. w = 1 will only output data on ending states.

a_default <- 0.50        # host type proportion ADENOSTOMA
A_default <- 0           # number of generations in allopatry before 2ndary contact
b_default <- 0           # hold mel allele freq constant (1 == yes, 0 == no)
B_default <- 0.50        # biased niche mating ratio
c_default <- 0           # cost of migration prob(0 to 1)
C_default <- 3           # number of COLUMNS of PATCHES
d_default <- 0.2         # prop_host_dark_niche prob(0 to 1)
D_default <- 0           # read in random number seed
f_default <- 0.300       # BASE_MATING_FITNESS
F_default <- 0.120       # CHC_MATCHING_ADVANTAGE
g_default <- 0.0572         # MELANISTIC_MATING_ADVANTAGE_M
G_default <- 0.0800         # MELANISTIC_MATING_ADVANTAGE_F
h_default <- 1           # print headers (1 == yes, 0 == no)
H_default <- 1           # DOMINANCE_COEFFICIENT of the UNSTRIPED_ALLELE (0.5 to 1)
L_default <- 1           # layout of hosts (0 = A clumped with A; 1 random patches)
m_default <- 0.06        # SD of MIGRATION distance MELANISTICS
M_default <- 0.04        # SD of MIGRATION distance GREENS
N_default <- 4000        # sum total N of all patches
o_default <- 1           # location of offspring within deme 1 = random location; 0 = offspring born at location of one parent
p_default <- 0.35        # initial FREQ of MELANISTIC ALLELE in AD
P_default <- 0.35        # initial FREQ of MELANISTIC ALLELE in CEA
r_default <- 0.5         # RECOMBINATION rate prob(0 to 1)
R_default <- 6           # number of ROWS of PATCHES
s_default <- 0.30        # S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE
S_default <- 0.2225      # S_COEFF_MELANISTIC_IN_GREEN_NICHE
t_default <- 0.1000      # S_COEFF_GREEN_IN_DARK_NICHE
T_default <- 1000        # max NUMBER of GENERATIONS to run the simulation
w_default <- 0           # 0 = single run; 1 = multiple runs.
Z_default <- 100         # frequency of data recording (in generations)


# Set-up base directory information.
# output files from the simulations will be written to this location:
base_dir <- "~/Documents/Research/Timema/modeling_dark_morph/MelanisticModel/Source"
setwd(base_dir)


##################################################
# single runs with pre-defined parameter values. #
##################################################
# set up parameters that will be varied across runs:

ps <- c(0, 0.35)  # frequency of melanistic allele
rs <- c(0.1, 0.5) # recombination rate

ms <- c(0.0266)   # sd of melanistic migration rate 
Ms <- c(0.0200)   # sd of green migration rate

Ss <- c(0.2225)    # s against melanistic in green niche
ts <- c(0.1000)    # s against green in dark niche

for (runs in seq_along(ms)) {
  # make a directory to write outputs for simulations using different migration rates and selection.
  system(paste("mkdir", paste("108_singleRun", ms[runs], Ms[runs], sep="_")))
  
  for (i in 1:200) { # define the number of independent simulations to run (here we run 200).
    for (p in ps) {
      for (r in rs) {
        # name output files with names describing parameters that were varried (here parameters "p" and "r"):
        prefix <- paste( paste( "108_singleRun", ms[runs], Ms[runs], sep="_"), '/', i, '_p', p, '_r', r, sep="" )
                        
        # call the run_nLocusSim function:
        run_nLocusSim(base_dir, prefix, a_default, A_default, b_default, B_default, c_default, C_default, d_default,
                      D_default, f_default, F_default, g_default, G_default, h_default, H_default, L_default, 
                      ms[runs], Ms[runs],   # migration parameters (m and M) defined on lines 89 and 90.
                      N_default, o_default,
                      p, p, r,              # freq of melanistic (p, P) and recombination rate (r) parameters defined on lines 86 and 87
                      R_default, s_default, 
                      Ss[runs], ts[runs],   # selection parameters defined on lines 92 and 93
                      T_default, w_default, Z_default)
      }
    }
  }
         
}

