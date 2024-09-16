
#' description of parameters
#'
#' @export
describe_dtlcor <- function(){
    cat("nsim:         number of replicates \n")
    cat("q_seq:        range of response rate to get alpha_t \n")
    cat("              (can be 100(1-2*alpha)% CI of q) \n")
    cat("gamma_seq:    range of hazards ratio of responder and non-responder \n")
    cat("              to get alpha_t (can be 100(1-2*alpha)% CI  of gamma) \n")
    cat("delta:        least difference to decide superiority of high dose \n")
    cat("n:            sample size per arm at DTL look \n")
    cat("D:            number of events \n")
    cat("N:            sample size per arm at final analysis \n")
    cat("alpha:        targeted FWER \n")
    cat("mPFS:         median progression-free survival for standard of care, \n")
    cat("              low dose and high dose \n")
    cat("q:            responser rate for standard of care, \n")
    cat("              low dose and high dose \n")
    cat("gamma:        hazards ratio of responder and non-responder \n")
    cat("drop_rate:    drop-out rate \n")
    cat("enroll:       enrollment rate \n")
    cat("interim_t:    information fraction at interim and final \n")
    cat("fix_rho:      NULL (use upper bound of rho) or values from 0 to 1 \n")
}

#' get the theoretical FWER given fixed correlation coefficient
#'
#' @export
dtl_tier_the <- function(delta, n, t, rho, q, alpha_s){

    if (length(t) == 1){

        c = qnorm(1-alpha_s)

        int_prob = function(z){
            sigma_11 = sqrt(2*q*(1-q)/n)
            d_1 = (sqrt(2)*delta/sigma_11 + z*rho)/sqrt(2-rho^2)
            d_2 = (sqrt(2)*delta/sigma_11 - z*rho)/sqrt(2-rho^2)
            int_prob = (pnorm(d_1) - pnorm(d_2))*dnorm(z)
            return(int_prob)
        }

        Type_I_Error =
            alpha_s +
            integrate(int_prob, lower = c, upper = Inf)$value

    } else if (length(t) > 1){

        OF_Design = gsDesign(k = length(t), test.type=1, sfu="OF", alpha = alpha_s, timing = t)
        c         = OF_Design$upper$bound

        int_prob = function(z){
            sigma_11 = sqrt(2*q*(1-q)/n)
            Sigma_22 = diag(length(t))
            for (i in 1:length(t)){
                for (j in i:length(t)){
                    Sigma_22[i,j] = sqrt(t[i]/t[j])
                    Sigma_22[j,i] = Sigma_22[i,j]
                }
            }
            d_1 = (sqrt(2)*delta/sigma_11 + t(rho)%*%solve(Sigma_22)%*%z)/sqrt(2-t(rho)%*%solve(Sigma_22)%*%rho)
            d_2 = (sqrt(2)*delta/sigma_11 - t(rho)%*%solve(Sigma_22)%*%z)/sqrt(2-t(rho)%*%solve(Sigma_22)%*%rho)
            int_prob = (pnorm(d_1) - pnorm(d_2))*dmvnorm(z, mean = rep(0,dim(Sigma_22)[1]), sigma = Sigma_22)
            return(int_prob)
        }

        Type_I_Error = alpha_s
        for (i in 1:length(t)){
            lower            = rep(-Inf, length(t))
            upper            = rep(Inf, length(t))
            lower[i]         = c[i]
            if (i != 1) {
                upper[1:(i-1)] = c[1:(i-1)]
            }
            Type_I_Error = Type_I_Error +
                cubature::cubintegrate(int_prob, lower = lower, upper = upper)$integral
        }

    }

    return(Type_I_Error)

}

#' get the alpha_s given significance level alpha
#'
#' @export
#'
dtl_get_alpha_s <- function(delta, n, t, rho, q, alpha){
    alpha_s = uniroot(
        function(x,alpha,delta,n,t,rho,q){
            dtl_tier_the(delta,n,t,rho,q,x) - alpha
        },
        alpha = alpha,
        delta = delta,
        n     = n,
        t     = t,
        rho   = rho,
        q     = q,
        lower = 0.001,
        upper = 0.5)$root

    return(alpha_s)
}

#' get theoretical upper bound of correlation coefficient
#'
#' @export
#'
dtl_cor_the_PH_upper_bound <- function(tau_k, pi_ar = 0.5, q, gamma){

    integrand = function(u) {
        return(1 / ( (1+q/(1-q)*u^(1-1/gamma))^2 * (q+(1-q)/gamma*u^(1/gamma-1)) ))
    }

    int = integrate(integrand, lower = 0, upper = 1)$value
    the = sqrt((1-pi_ar) * q/(1-q) * tau_k) * (1/gamma-1) * sqrt(int)
    return(the)

}
