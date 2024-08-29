
#' get survival function of non-responder
#'
#' @export
#'
get_S0 <- function(mPFS, q, gamma, t){
    solve_S0 = function(x){
        S = 1 - pexp(t, log(2) / mPFS)
        q*x^gamma + (1-q)*x - S
    }
    uniroot(solve_S0, interval = c(0, 1), tol = 10^(-6))$root
}

#' get density function of non-responder
#'
#' @export
#'
get_f0 <- function(mPFS, q, gamma, t){
    S0 = get_S0(mPFS, q, gamma, t)
    f0 = dexp(t, log(2) / mPFS) / (q*gamma*S0^(gamma-1) + (1-q))
    f0
}

#' get alpha_t (minimum of alpha_s) given alpha and ranges of q and gamma
#'
#' @export
#'
dtl_app_get_alpha_t = function(n, N, delta, q_seq, gamma_seq, alpha, fix_rho = NULL){
    
    # sample size calculation (naive and ignoring correlation)
    seq_all = expand.grid(q_seq, gamma_seq)
    
    alpha_s_all = apply(seq_all, 1, function(x){
        q_s     = x[1]
        gamma_s = x[2]
        
        # calculate rho
        if (!is.null(fix_rho)){
            rho = fix_rho
        } else{
            rho = dtl_cor_the_PH_upper_bound(tau_k = n / N,
                                             pi_ar = 0.5,
                                             q     = q_s, 
                                             gamma = gamma_s)
        }
        alpha_s  = dtl_get_alpha_s(delta, n, 1, rho, q_s, alpha) # no interim t = 1
        data.frame(q_s, gamma_s, rho, alpha_s)
    })
    rst_alpha_s = data.frame(t(sapply(1:length(alpha_s_all), function(x){as.numeric(alpha_s_all[[x]])})))
    colnames(rst_alpha_s) = c("q", "gamma", "rho", "alpha_s")
    
    rst_alpha_t = unique(rst_alpha_s %>%
                             filter(alpha_s == min(rst_alpha_s$alpha_s)) %>%
                             rename_with(~"alpha_t", "alpha_s"))
    
    list(rst_alpha_t = rst_alpha_t, rst_alpha_s = rst_alpha_s)
    
}


#' simulate a single trial
#'
#' @export
#'
dtl_app_sim_single <- function(D, N, n, mPFS, q, gamma, delta, drop_rate, follow_time, enroll, t){

    accr_time = N / enroll # day per arm

    q_rep    = rep(q, each = N)
    mPFS_rep = rep(mPFS, each = N)
    lambda_C = -log(1-drop_rate) / 365.25 # censor survival at 1 year equals 1 - drop_rate

    ID       = 1:(3*N)
    arm      = rep(0:2, each = N)

    Eve_Time = rexp(3*N, log(2) / mPFS_rep)

    f0       = apply(cbind(mPFS_rep, q_rep, Eve_Time), 1, function(x){
        get_f0(x[1], x[2], gamma, x[3])
    })
    q_t      = 1 - (1-q_rep)*f0 / dexp(Eve_Time, log(2)/mPFS_rep)
    X        = rbinom(3*N, 1, q_t)

    Cen_Time = rexp(3*N, lambda_C)
    tt_accr  = runif(3*N, 0, accr_time)
    tt_eve   = tt_accr + Eve_Time
    tt_cen   = tt_accr + Cen_Time

    dat_final_temp = tibble(ID, arm, X, Eve_Time, Cen_Time, tt_accr, tt_eve, tt_cen) %>%
        mutate(censor = case_when(tt_eve < tt_cen & tt_eve < tt_accr + follow_time   ~ 0,
                                  tt_cen <= tt_eve & tt_cen < tt_accr + follow_time  ~ 1,
                                  tt_accr + follow_time <= tt_cen & tt_accr + follow_time <= tt_eve ~ 2),
               tt     = case_when(censor == 0 ~ tt_eve,
                                  censor == 1 ~ tt_cen,
                                  censor == 2 ~ tt_accr + follow_time)) %>%
        arrange(tt)

    # DTL stage
    dat_DTL = dat_final_temp %>% arrange(tt_accr) %>% filter(arm !=0) %>% group_by(arm) %>% slice(1:n)
    rst_W   = dat_DTL %>% group_by(arm) %>% summarise(W = mean(X))

    W_names   = c("W_1", "W_2")
    W_dat_all = data.frame(rbind(rst_W$W))
    colnames(W_dat_all) = W_names

    if (W_dat_all$W_2 - W_dat_all$W_1 - delta <= 0){
        dat_final_temp_2 = dat_final_temp %>% filter(arm != 2)
    } else {
        dat_final_temp_2 = dat_final_temp %>% filter(arm != 1)
    }

    # Final stage
    dat_final_temp_3 = dat_final_temp_2 %>%
        mutate(D_cumsum = cumsum(censor==0))

    t_length  = length(t)
    dat_final = list()
    Z_all     = NULL
    for (k in 1:t_length){

        D_k = ceiling(D*t[k])

        if (max(dat_final_temp_3$D_cumsum) < D_k){
            tt_end = max(dat_final_temp_3$tt)
        } else{
            tt_end = min(dat_final_temp_3$tt[dat_final_temp_3$D_cumsum == D_k])
        }

        dat_final[[k]] = dat_final_temp_3 %>%
            mutate(tt_end = tt_end) %>%
            mutate(censor = if_else(tt > tt_end, 3, censor),
                   tt     = if_else(censor != 3, tt, tt_end),
                   Time   = tt - tt_accr,
                   Delta  = if_else(censor == 0, 1, 0))

        rst_test = logrank_test(Surv(Time, Delta) ~ factor(arm), data = dat_final[[k]])
        Z        = -rst_test@statistic@teststatistic

        Z_all = c(Z_all,
                  c(if_else(1 %in% dat_final[[k]]$arm, Z, NA),
                    if_else(2 %in% dat_final[[k]]$arm, Z, NA)))

    }

    Z_names = apply(expand.grid(1:2, 1:t_length), 1, function(x){
        paste0("Z_", x[1], x[2])
    })

    Z_dat_all           = data.frame(rbind(Z_all))
    rownames(Z_dat_all) = NULL
    colnames(Z_dat_all) = Z_names

    dat_WZ = data.frame(W_dat_all, Z_dat_all)

    return(list(dat_WZ = dat_WZ, dat_final = dat_final))

}

#' analysis data
#'
#' @export
#'
dtl_app_ana <- function(dat_all, t, c){

    t_length     = length(t)
    dat_WZ       = dat_all$dat_WZ
    dat_final    = dat_all$dat_final[[length(dat_all$dat_final)]]

    stop_interim = if_else(1 %in% dat_final$arm,
                           which(dat_WZ[1, 1 + 2*(1:t_length)] > c)[1],
                           which(dat_WZ[1, 2 + 2*(1:t_length)] > c)[1])
    stop_interim = if_else(is.na(stop_interim), length(dat_all$dat_final), stop_interim)


    I_21  = 2 %in% dat_final$arm
    rej   = case_when(1 %in% dat_final$arm & !all(dat_WZ[1, 1 + 2*(1:t_length)] <= c) ~ 1,
                      2 %in% dat_final$arm & !all(dat_WZ[1, 2 + 2*(1:t_length)] <= c) ~ 2,
                      .default = 0)

    cen_rate  = mean(dat_final$censor!=0)
    cen_1     = mean(dat_final$censor==1)
    cen_2     = mean(dat_final$censor==2)
    cen_3     = mean(dat_final$censor==3)
    dur       = dat_all$dat_final[[stop_interim]]$tt_end[1] - min(dat_final$tt_accr)

    rst = data.frame(t = rbind(t), c = rbind(c), I_21, rej, cen_rate, cen_1, cen_2, cen_3, dur)
    rownames(rst) = NULL

    return(rst)

}

#' simulation function
#'
#' @export
#'
dtl_app_sim <- function(nsim, alpha_t,
                       D, N, n, mPFS, q, gamma, delta,
                       drop_rate, follow_time, enroll, t){

    t_length  = length(t)

    OF_Design = gsDesign(k = length(t), test.type=1, sfu="OF", alpha = alpha_t, timing = t)
    c         = OF_Design$upper$bound

    rst_all = NULL
    for (i in 1:nsim){
        # simulate data
        dat_all  = dtl_app_sim_single(D, N, n, mPFS, q, gamma, delta, drop_rate, follow_time, enroll, t)

        # analysis data
        rst     = dtl_app_ana(dat_all, t, c)
        rst_all = rbind(rst_all, data.frame(rep = i, rst))
    }

    rst_all_wide = rst_all %>%
        mutate(prob_21 = I_21,
               rej_12  = (rej == 1 | rej ==2),
               rej_1   = (rej == 1),
               rej_2   = (rej == 2)) %>%
        dplyr::select(!c(rej, I_21))

    rst_dtl = data.frame(mPFS_0      = mPFS[1],
                         mPFS_1      = mPFS[2],
                         mPFS_2      = mPFS[3],
                         q_0         = q[1],
                         q_1         = q[2],
                         q_2         = q[3],
                         gamma       = gamma,
                         delta       = delta,
                         drop_rate   = drop_rate,
                         follow_time = follow_time,
                         enroll      = enroll,
                         D           = D,
                         N           = N,
                         n           = n,
                         alpha_t     = alpha_t,
                         rbind(apply(rst_all_wide[, -1], 2, mean)))

    return(rst_dtl)
}

