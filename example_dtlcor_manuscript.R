
# install dtlcor package
rm(list = ls())
options(error = recover)
require(devtools)
install_github("ashwz/dtlcor")
require(dtlcor)

# description of parameters
# describe_dtlcor()

# set design parameters
nsim      = 1000  
q_seq     = seq(0.19, 0.32, 0.01) 
gamma_seq = seq(0.14, 0.34, 0.01) 

delta       = 0.05  
n           = 80    
D           = 162
N           = 152   
alpha       = 0.025
mPFS        = c(180, 276, 300)
q           = c(0.2, 0.4, 0.5)
gamma       = 0.15
drop_rate   = 0.05
follow_time = 10^6
enroll      = 240 / 365.25
t           = c(0.5, 1)
fix_rho     = NULL

# simulate a single trial
set.seed(1000)
dtl_single = dtl_app_sim_single(D, N, n, mPFS, q, gamma, delta, 
                                drop_rate, follow_time, enroll, t)

# get alpha_t (minimum of alpha_s)
# values of fix_rho: NULL (use upper bound of rho) or real values from 0 to 1.
alpha_t = dtl_app_get_alpha_t(N, n, q_seq, gamma_seq, 
                              alpha, fix_rho = fix_rho)$alpha_t[1]

# simulation results
set.seed(1000)
dtl_summary = dtl_app_sim(nsim, alpha_t, D, N, n, mPFS, q, gamma, delta, 
                          drop_rate, follow_time, enroll, t)
dtl_summary %>% 
    select(alpha_t, c.1, c.2, cen_rate, dur, prob_21, rej_12, rej_1, rej_2)
