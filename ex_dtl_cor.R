
rm(list = ls())
options(error = recover)
require(devtools)

# for local package
# load_all("dtlcor")
# document("dtlcor")

# for online package
# install_github("ashwz/dtlcor")

require(dtlcor)
ls(getNamespace("dtlcor"))

# description of parameters
describe_dtlcor()

# set design parameters
nsim      = 1000  

delta     = 0.05  
n         = 80    
N         = 152   
alpha     = 0.025

q_seq     = seq(0.19, 0.32, 0.01) 
gamma_seq = seq(0.14, 0.34, 0.01) 

D           = 162
mPFS        = c(180, 276, 300)
q           = c(0.2, 0.4, 0.5)
gamma       = 0.15
drop_rate   = 0.05
follow_time = 10^6
enroll      = 240 / 365.25
t           = c(0.5, 1)
fix_rho     = NULL

# --------------------basic function-----------------
# theoretical FWER
dtl_tier_the(delta, n, t = 1, rho = 0.4, q = 0.3, alpha_s = alpha)

# get alpha_s given alpha
dtl_get_alpha_s(delta, n, t = 1, rho = 0.4, q = 0.3, alpha)

# get upper bound of correlation
dtl_cor_the_PH_upper_bound(tau_k = n/N, pi_ar = 0.5, q = 0.3, gamma = 0.2)

# --------------------real application-----------------
# get alpha_t (minimum of alpha_s)
# values of fix_rho: NULL (use upper bound of rho) or real values from 0 to 1.
dtl_app_get_alpha_t(N, n, q_seq, gamma_seq, alpha, fix_rho = NULL)
dtl_app_get_alpha_t(N, n, q_seq, gamma_seq, alpha, fix_rho = 1)
dtl_app_get_alpha_t(N, n, q_seq, gamma_seq, alpha, fix_rho = 0)

# simulate a single trial
dtl_single = dtl_app_sim_single(D, N, n, mPFS, q, gamma, delta, drop_rate, follow_time, enroll, t)

# simulation results
alpha_t     = dtl_app_get_alpha_t(N, n, q_seq, gamma_seq, alpha, fix_rho = fix_rho)$alpha_t[1]
dtl_summary = dtl_app_sim(nsim, alpha_t, D, N, n, mPFS, q, gamma, delta, drop_rate, follow_time, enroll, t)

