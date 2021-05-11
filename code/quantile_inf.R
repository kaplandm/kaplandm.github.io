# Feedback: KaplanDM@Missouri.edu
# Author (original): David M. Kaplan
# 29 March 2013 (and subsequent updates)
# Confidence intervals (etc.) for various quantile-related objects.
# Examples: http://faculty.missouri.edu/kaplandm/code/quantile_inf_examples_qte.R
# 
# References (chronological):
# 'Calculating nonparametric confidence intervals for quantiles using fractional order statistics' by Alan D. Hutson, 1999, https://doi.org/10.1080/02664769922458
# 'Improved quantile inference via fixed-smoothing asymptotics and Edgeworth expansion' by David M. Kaplan, 2015, https://doi.org/10.1016/j.jeconom.2014.08.011
# 'Fractional order statistic approximation for nonparametric conditional quantile inference' by Matt Goldman and David M. Kaplan, 2017, https://doi.org/10.1016/j.jeconom.2016.09.015
# 'Non-parametric inference on conditional quantile differences and linear combinations, using L-statistics' by Matt Goldman and David M. Kaplan, 2018, https://doi.org/10.1111/ectj.12095
# 
# See also http://faculty.missouri.edu/kaplandm/code/quantile_inf_np.R for related nonparametric inference on conditional quantile-related objects


###################
#                 #
# Main function   #
#                 #
###################
# Function
# # Given sample data and quantile(s) of interest, this function can produce a confidence interval (CI) for a single quantile of a population distribution, for multiple quantiles jointly, for a linear combination of multiple quantiles, for quantile treatment effects (given treatment and control samples), or for treatment effects on linear combinations of quantiles.  See below (and examples in quantile_inf_examples.R) for details.
# Arguments
# # X: numeric vector of data, or for METHOD.TYPE 'qte', X$t is the numeric vector of treatment group data and X$c is the numeric vector of control group data.
# # p: quantile(s) of interest, a numeric vector (or scalar) with values strictly between zero and one; e.g., the 0.5-quantile is the median.  Should be scalar for METHOD.TYPE 'single'; should be vector with length(p)>1 for METHOD.TYPE 'joint' or 'lincom'; can be scalar or vector for METHOD.TYPE 'qte'.
# # ALPHA: the nominal coverage probability of the returned confidence interval (CI) is 1-ALPHA; e.g., ALPHA=0.05 produces a 95% CI.
# # ONESIDED: ONESIDED=0 produces a two-sided CI; ONESIDED=1 produces an upper one-sided CI; ONESIDED=-1 produces a lower one-sided CI.
# # METHOD.TYPE: 'single' produces a CI for a single quantile using one sample of data; 'joint' produces a multidimensional CI (the Cartesian product of the separate intervals returned) for multiple quantiles with 1-ALPHA nominal joint coverage probability; 'lincom' produces a CI for a linear combination of quantiles using the specified vector of weights PSI (corresponding to specified vector of quantiles p), e.g. p=c(0.75,0.25) and PSI=c(1,-1) for the interquartile range; 'qte' produces a CI for the quantile treatment effect (or, under fewer assumptions, simply the difference between the p-quantiles of two populations) if length(p)==1, or a CI for the treatment effect on a linear combination of quantiles if length(p)>1 and length(PSI)==length(p).  The two numeric vectors of data from the two (i.e., treatment and control) populations should be passed in X$t and X$c, respectively.
# # PSI: numeric vector specifying weights for linear combinations, where first element should be normalized to be 1 and length(PSI)==length(p)>1, for METHOD.TYPE 'lincom' or 'qte'.  For METHOD.TYPE 'single' or 'joint', or 'qte' with length(p)==1, PSI should be NULL.
# # NORM.APPROX.CAL: if true, then (where implemented) use analytic normal approximation for alpha calibration; if false, then use joint beta simulation for calibration
# # NORM.APPROX.Q: if true, then (where implemented) use analytic normal approximation for fractional order statistic index determination
# # SPACING.FLAG: if true, then for nuisance parameter estimation, use quantile spacing density estimator with Goldman-Kaplan zero-bias bandwidth; if false, then use standard kernel density estimator
# # BETA.BLOCK.SIZE, BETA.BLOCKS: for simulating joint beta distributions of fractional order statistics, take BETA.BLOCK.SIZE draws and repeat BETA.BLOCKS times, so that the total number of draws is BETA.BLOCK.SIZE*BETA.BLOCKS
# # GAMMA allows the PDF ratio nuisance parameter value(s) to be specified by the user rather than estimated.  In general, they should be estimated, so the default of NULL is recommended.  If non-NULL, for QTE, must have length(GAMMA$t)==length(GAMMA$c)==length(p), or length(GAMMA)==length(p) for lincom method
# # NUM.INT: for QTE with length(p)==1, use numerical integration rather than simulation to greatly improve speed (with same accuracy).
# # QTE.GK.ONLY: do not use Kaplan or bootstrap for QTE when length(p)==1 if Goldman-Kaplan cannot be used--just return NA.  Note: primary purpose is for simulation studies, not practice.
# # QTE.NO.GK: do not use Goldman-Kaplan for QTE when length(p)==1, but rather Kaplan (or bootstrap).  Note: primary purpose is for nonparametric quantile marginal effect inference when local sample is small.
# # QTE.NO.BS: only use Goldman-Kaplan or Kaplan methods for QTE with length(p)==1; do not use bootstrap.  Note: even if FALSE, bootstrap is only used if 1) GK is not computable, 2) treatment and control samples are very different sizes.  Set TRUE if computation time seems slow (esp. if used as subroutine of nonparametric QME inference).
# SINGLE.CALIB: if TRUE then use the calibrated Hutson CI from (11) in Goldman & Kaplan; only applies to METHOD.TYPE 'single'
#
# Return value: a list with the following components.
# # methname: 'Hutson' for the Hutson (1999) single-quantile CI, 'GoldmanKaplan' for any of the Goldman and Kaplan (2017/X) methods, 'Kaplan' for the Kaplan (2015) 'single' or 'qte' CI, 'bootstrap-perc-t' for a bootstrap percentile-t with bootstrap-estimated variance used in Studentization.
# # CI: data.frame with 'lo' and 'hi' the lower and upper endpoint(s), respectively, of confidence interval(s).  For METHOD.TYPE 'single', 'lincom', or 'qte', dim(CI)[1]==length(CI$lo)==length(CI$hi)==1.  For METHOD.TYPE 'joint', the overall joint CI is the Cartesian product of the individual CIs specified, and dim(CI)[1]==length(CI$lo)==length(CI$hi)==length(p)
# # ALPHAtilde: calibrated alpha used for individual CIs for METHOD.TYPE 'qte', 'joint', or 'lincom'
#
quantile.inf <- function(X,p,ALPHA=0.05,ONESIDED=0,METHOD.TYPE='single',PSI=NULL,NORM.APPROX.CAL=FALSE, NORM.APPROX.Q=FALSE, SPACING.FLAG=TRUE, BETA.BLOCK.SIZE=10^4, BETA.BLOCKS=5, GAMMA=NULL, NUM.INT=TRUE, QTE.GK.ONLY=FALSE, QTE.NO.BS=FALSE, QTE.NO.GK=FALSE, SINGLE.CALIB=FALSE) {

  # Check arguments: errors
  if (missing(X)) stop("Argument X is missing.")
  if (missing(p)) stop("Argument p is missing.")
  if (!is.logical(NUM.INT)) stop("Argument NUM.INT (numerical integration instead of simulation for QTE method on single quantile) must be TRUE or FALSE.")
  if (is.na(NUM.INT)) stop("Argument NUM.INT must be TRUE or FALSE (not NA).")
  quantile.inf.validate.inputs(X=X,p=p,ALPHA=ALPHA,ONESIDED=ONESIDED,METHOD.TYPE=METHOD.TYPE,PSI=PSI,NORM.APPROX.CAL=NORM.APPROX.CAL,NORM.APPROX.Q=NORM.APPROX.Q,GAMMA=GAMMA)
  for (varname in c('QTE.GK.ONLY','QTE.NO.BS','QTE.NO.GK')) {
    if (!is.logical(get(varname)) || is.na(get(varname))) stop(sprintf("%s must be TRUE or FALSE.",varname))
  }
  if (QTE.GK.ONLY && QTE.NO.GK) stop("Contradiction: QTE.GK.ONLY and QTE.NO.GK are both TRUE; one (or both) must be FALSE.")
  
  # Check arguments: warnings
  if (!is.null(ALPHA) && ALPHA>0.2) warning('ALPHA is usually 0.01, 0.05, or 0.10; e.g., 0.05 gives a 95% confidence interval.')
  if (length(p)==1 && METHOD.TYPE %in% c('joint','lincom')) { warning("Setting METHOD.TYPE to 'single' since length(p)==1"); METHOD.TYPE <- 'single' }
  if (!is.null(PSI) && METHOD.TYPE %in% c('single','joint')) warning(sprintf("Ignoring argument PSI since METHOD.TYPE=='%s'",METHOD.TYPE))
  if (METHOD.TYPE %in% c('lincom','qte') && length(p)==1) PSI <- 1

  # Constants
  PDFMIN<-1e-4 #avoid dividing by zero for 2-sample Hutson
  UNIROOT.TOL <- 0.0001
  UNIROOT.INT <- c(0.000000001,1-0.000000001)
  KERNEL.TYPE <- 'gaussian'
  BW.TYPE <- 'SJ'

  # Prep data
  if (METHOD.TYPE=='qte') { X$t <- sort(X$t);  X$c <- sort(X$c) } else { X <- sort(X) } #X[i] is now order statistic i
  vec1s <- rep.int(ONESIDED,length(p))
  if (METHOD.TYPE %in% c('qte','lincom') && length(p)>1) vec1s <- ONESIDED*sign(PSI)
  
  # Return if METHOD.TYPE 'single'
  if (METHOD.TYPE=='single') {
    if (SINGLE.CALIB) {
      ret <- tryCatch(list(methname='Hutson',CI=quantile.inf.single.GK(X=X,p=p,ALPHA=ALPHA,ONESIDED=ONESIDED,NORM.APPROX=NORM.APPROX.Q,CALIB=SINGLE.CALIB)), 
                      error=function(CIE) NULL)
      if (!is.null(ret)) return(ret)
    }
    ret <- tryCatch(list(methname='Hutson',CI=quantile.inf.single.GK(X=X,p=p,ALPHA=ALPHA,ONESIDED=ONESIDED,NORM.APPROX=NORM.APPROX.Q)), 
                    error=function(CIE) list(methname='Kaplan',CI=quantile.inf.single.K(X=X,p=p,ALPHA=ALPHA,ONESIDED=ONESIDED)))
    return(ret)
  }

  # Calibrate alpha
  ALPHAtilde <- ALPHA
  if (METHOD.TYPE=='single') stop("Should have returned already for METHOD.TYPE 'single'")
  ALPHAtilde <- quantile.inf.calibrate(X=X,p=p,ALPHA=ALPHA, ONESIDED=ONESIDED, METHOD.TYPE=METHOD.TYPE, PSI=PSI, NORM.APPROX=NORM.APPROX.CAL, KERNEL.TYPE=KERNEL.TYPE, BW.TYPE=BW.TYPE, SPACING.FLAG=SPACING.FLAG, BETA.BLOCK.SIZE=BETA.BLOCK.SIZE, BETA.BLOCKS=BETA.BLOCKS, GAMMA=GAMMA, NUM.INT=NUM.INT)
  
  # Calculate individual CIs at each quantile
  if (METHOD.TYPE=='qte') {
    CI.ind <- list(t=data.frame(lo=rep.int(-Inf,length(p)), hi=rep.int(Inf,length(p))), c=data.frame(lo=rep.int(-Inf,length(p)), hi=rep.int(Inf,length(p))))
  } else {
    CI.ind <- data.frame(lo=rep.int(-Inf,length(p)), hi=rep.int(Inf,length(p)))
  }
  for (i in 1:length(p)) {
    if (length(p)==1) {
      tmp.fn <- function(CIE) c(Inf,-Inf)
    } else {
      tmp.fn <- function(CIE) stop("Could not compute with GK approach; an alternative method and possibly extreme value theory is needed for this combination of sample size, quantiles, and nominal level.  Original error message: ",CIE$message,call.=FALSE)
    }
    if (METHOD.TYPE=='qte') {
      if (QTE.NO.GK && length(p)==1) {
        CI.ind$t[i,] <- CI.ind$c[i,] <- tmp.fn(0)
      } else {
        CI.ind$t[i,] <- tryCatch(quantile.inf.single.GK(X=X$t,p=p[i],ALPHA=ALPHAtilde,ONESIDED=vec1s[i],NORM.APPROX=NORM.APPROX.Q), error=tmp.fn)
        CI.ind$c[i,] <- tryCatch(quantile.inf.single.GK(X=X$c,p=p[i],ALPHA=ALPHAtilde,ONESIDED=-vec1s[i],NORM.APPROX=NORM.APPROX.Q), error=tmp.fn)
      }
      if (all(CI.ind$t[i,]==c(Inf,-Inf)) || all(CI.ind$c[i,]==c(Inf,-Inf))) {
        if (QTE.GK.ONLY) {
          return(list(methname='GoldmanKaplan',CI=list(lo=NA,hi=NA),ALPHAtilde=NA))
        }
        if (length(p)==1 && (abs((length(X$t)-length(X$c))/max(length(X$t),length(X$c)))<=0.5 || QTE.NO.BS)) { #If sample sizes similar or QTE.NO.BS==TRUE, use Kaplan 2015
          return(list(methname='Kaplan',CI=quantile.inf.qte.K(X=X,p=p,ALPHA=ALPHA,ONESIDED=ONESIDED,PDFratio=GAMMA$c)))
        } else if (length(p)==1 && !QTE.NO.BS) {
          return(list(methname='bootstrap-perc-t',CI=quantile.inf.qte.BSt(X=X,p=p,ALPHA=ALPHA,ONESIDED=ONESIDED)))
        } else { stop("Unexpected error; please report to KaplanDM@Missouri.edu") }
      }
    } else {
      CI.ind[i,] <- tryCatch(quantile.inf.single.GK(X=X,p=p[i],ALPHA=ALPHAtilde,ONESIDED=vec1s[i],NORM.APPROX=NORM.APPROX.Q), error=tmp.fn)
    }
  }
  
  # Combine individuals CIs into overall CI
  CI <- data.frame(lo=-Inf,hi=Inf)
  if (METHOD.TYPE=='joint') {
    return(list(methname='GoldmanKaplan',CI=CI.ind,ALPHAtilde=ALPHAtilde))
  } else if (METHOD.TYPE=='lincom') {
    tmphi <- (sign(PSI)+1)/2 + 1
    tmplo <- 3-tmphi
    if (ONESIDED==-1) {
      CI$hi <- sum(PSI*CI.ind[cbind(1:length(p),tmphi)])
    } else if (ONESIDED==1) {
      CI$lo <- sum(PSI*CI.ind[cbind(1:length(p),tmplo)])
    } else { #two-sided
      CI$lo <- sum(PSI*CI.ind[cbind(1:length(p),tmplo)])
      CI$hi <- sum(PSI*CI.ind[cbind(1:length(p),tmphi)])
    }
    return(list(methname='GoldmanKaplan',CI=CI,ALPHAtilde=ALPHAtilde))
  } else { # qte
    tmphi <- (sign(PSI)+1)/2 + 1
    tmplo <- 3-tmphi
    if (ONESIDED==-1) {
      CI$hi <- sum(PSI*CI.ind$t[cbind(1:length(p),tmphi)]) - sum(PSI*CI.ind$c[cbind(1:length(p),tmplo)])
    } else if (ONESIDED==1) {
      CI$lo <- sum(PSI*CI.ind$t[cbind(1:length(p),tmplo)]) - sum(PSI*CI.ind$c[cbind(1:length(p),tmphi)])
    } else { #two-sided
      CI$lo <- sum(PSI*CI.ind$t[cbind(1:length(p),tmplo)]) - sum(PSI*CI.ind$c[cbind(1:length(p),tmphi)])
      CI$hi <- sum(PSI*CI.ind$t[cbind(1:length(p),tmphi)]) - sum(PSI*CI.ind$c[cbind(1:length(p),tmplo)])
    }
    return(list(methname='GoldmanKaplan',CI=CI,ALPHAtilde=ALPHAtilde))
  }

} #end of function quantile.inf()



###############################################################
# Pre-calculated simulated critical values from Kaplan (2015) #
###############################################################
sub.simcv <- function() {
simalphas <- c(0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,
             0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,
             0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,
             0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,
             0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,
             0.49,0.5,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.6,
             0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.7,0.71,0.72,
             0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.8,0.81,0.82,0.83,0.84,
             0.85,0.86,0.87,0.88,0.89,0.9,0.91,0.92,0.93,0.94,0.95,0.96,
             0.97,0.98,0.99)
simCVs.bi <- matrix(0,2,length(simalphas))
simCVs.bi[1,] <-
  c(10.9435,8.8905,7.8981,7.2232,6.7377,6.351,6.0229,5.7584,
    5.5418,5.3417,4.1981,3.6015,3.2189,2.9354,2.714,2.535,2.3849,2.2609,
    2.149,2.0503,1.9611,1.8804,1.8057,1.7389,1.6765,1.6177,1.5631,
    1.5123,1.4647,1.4196,1.3779,1.3379,1.2998,1.2629,1.2281,1.1953,
    1.1624,1.1308,1.1012,1.0723,1.0442,1.0169,0.99113,0.96603,0.94104,
    0.91779,0.8943,0.8723,0.85036,0.82907,0.80833,0.7879,0.76856,
    0.74943,0.73068,0.7125,0.69448,0.67636,0.65898,0.64175,0.62518,
    0.60846,0.59187,0.57614,0.56043,0.54512,0.53009,0.51467,0.49973,
    0.48494,0.47025,0.45579,0.44124,0.42726,0.4134,0.3995,0.38613,
    0.37253,0.35938,0.34629,0.33318,0.32009,0.30725,0.29457,0.28195,
    0.26964,0.25698,0.24456,0.23229,0.22018,0.20801,0.19615,0.18394,
    0.17234,0.16065,0.14903,0.13742,0.12555,0.11394,0.10242,0.090718,
    0.079071,0.067377,0.055898,0.04451,0.033114,0.021501,0.010147)
simCVs.bi[2,] <-
  c(5.9696,5.2541,4.8539,4.5537,4.3365,4.1601,4.0172,3.8932,3.7904,
    3.6945,3.1109,2.7768,2.5495,2.3795,2.2387,2.1231,2.0241,1.9387,
    1.8618,1.7919,1.7286,1.6696,1.6155,1.5656,1.5192,1.4756,1.4325,
    1.3925,1.3547,1.3189,1.2847,1.2522,1.2219,1.1917,1.1633,1.1353,
    1.1082,1.0819,1.0566,1.032,1.0086,0.98543,0.96225,0.94023,0.91818,
    0.89684,0.87671,0.85666,0.83658,0.81765,0.79923,0.78084,0.76276,
    0.74445,0.72695,0.70947,0.69251,0.67591,0.65949,0.64369,0.62773,
    0.6121,0.59668,0.58137,0.56623,0.55174,0.53661,0.52165,0.50672,
    0.49258,0.47876,0.46438,0.45022,0.43646,0.42261,0.40875,0.39521,
    0.38141,0.36817,0.35498,0.34185,0.3288,0.31576,0.30296,0.29003,
    0.27721,0.26461,0.25188,0.23945,0.2274,0.21519,0.20266,0.19044,
    0.17843,0.16643,0.15419,0.14254,0.13036,0.11814,0.10607,0.094184,
    0.08209,0.069939,0.057994,0.046069,0.034169,0.0224,0.010453)
simCVs.uni <- matrix(0,2,length(simalphas))
simCVs.uni[1,] <-
c(44.8943,31.1802,25.2111,21.6883,19.2587,17.4917,16.1646,14.8876,
 13.9785,13.1943,8.939,7.0955,5.9917,5.226,4.6791,4.2355,3.8938,3.6012,
 3.3534,3.1435,2.9647,2.8014,2.6569,2.5254,2.4059,2.2988,2.1986,2.1045,
 2.0205,1.9421,1.8691,1.8022,1.7373,1.6756,1.617,1.5641,1.5122,1.4644,
 1.4184,1.3736,1.3299,1.289,1.249,1.2112,1.1749,1.1398,1.107,1.0753,
 1.044,1.0146,0.98538,0.95693,0.92908,0.90212,0.87608,0.85127,0.82678,
 0.80332,0.78023,0.75758,0.7353,0.71397,0.69282,0.67244,0.65224,0.63245,
 0.61313,0.5938,0.57494,0.55698,0.5391,0.52137,0.50417,0.4875,0.47138,
 0.45488,0.43891,0.42298,0.40778,0.39269,0.37756,0.3624,0.34747,0.33319,
 0.3183,0.30394,0.2894,0.27537,0.26169,0.24795,0.23425,0.22045,0.20687,
 0.19345,0.18045,0.16714,0.15436,0.14132,0.12855,0.11578,0.1028,
 0.090287,0.077914,0.064961,0.052481,0.040066,0.027524,0.014898)
simCVs.uni[2,] <-
c(11.8392,9.6368,8.5465,7.7846,7.2393,6.8183,6.502,6.2144,5.9611,5.7361,
4.531,3.8996,3.4738,3.1766,2.942,2.7496,2.5918,2.4522,2.3315,2.2232,
2.1273,2.0418,1.9634,1.8916,1.8254,1.7632,1.7051,1.6499,1.5986,1.5497,
1.5024,1.4595,1.4175,1.3791,1.3411,1.3046,1.27,1.2356,1.2038,1.1726,
1.1425,1.1137,1.0863,1.059,1.0322,1.0064,0.98103,0.9573,0.93325,0.91049,
0.88797,0.86565,0.84382,0.82256,0.80216,0.78183,0.76214,0.74232,0.72384,
0.705,0.6862,0.66792,0.64977,0.63206,0.61481,0.598,0.58103,0.56456,
0.54839,0.53288,0.51723,0.50117,0.4856,0.47051,0.45574,0.44051,0.42602,
0.41153,0.39682,0.38232,0.36861,0.35478,0.34084,0.32665,0.313,0.29897,
0.28536,0.27191,0.25866,0.24548,0.23193,0.21882,0.20566,0.19261,0.17974,
0.16659,0.15379,0.14096,0.12814,0.11529,0.10286,0.090379,0.077713,
0.064897,0.052336,0.040043,0.027465,0.014858)
return(list(simalphas=simalphas,simCVs.bi=simCVs.bi,simCVs.uni=simCVs.uni))
}


#########################
#                       #
# Check input args      #
#                       #
#########################
# Check for valid input data
quantile.inf.validate.inputs <- function(X,p,ALPHA,ONESIDED,METHOD.TYPE,PSI,NORM.APPROX.CAL,NORM.APPROX.Q,GAMMA) {
  for (varname in c('X','p','ALPHA','ONESIDED','METHOD.TYPE')) {
    if (is.null(get(varname))) stop(sprintf("%s cannot be NULL.",varname))
  }
  if (!is.character(METHOD.TYPE)) { stop("Need is.character(METHOD.TYPE)=TRUE") }
  if (!(METHOD.TYPE %in% c('single','qte','joint','lincom'))) stop("Argument METHOD.TYPE must be one of: 'single', 'qte', 'joint', 'lincom', which correspond respectively to one-sample inference on a single quantile, two-sample quantile treatment effect inference (including treatment effects on linear combinations of quantiles), one-sample joint inference on multiple quantiles, and one-sample inference on a linear combination of quantiles (such as the interquartile range).")
  if (is.list(X)) {
    if (METHOD.TYPE!='qte') stop(sprintf("X should be a numeric vector for METHOD.TYPE '%s'",METHOD.TYPE))
    chk.num.vec.fn(X$t);  chk.num.vec.fn(X$c)
  } else {
    if (METHOD.TYPE=='qte') stop("X should be a list with X$t and X$c containing the two data samples for METHOD.TYPE 'qte'")
    chk.num.vec.fn(X)
  }
  if (missing(p)) stop("Need to supply quantile of interest, e.g. p=0.5 for median")
  chk.num.vec.fn(p)
  if (min(p)<0 || max(p)>1) stop('Argument p must be between 0 and 1; e.g., 0.5 for the median.')
  if (METHOD.TYPE=='single') {
    if (length(p)!=1) stop('For METHOD.TYPE "single", p should be a scalar.')
  }
  if (METHOD.TYPE %in% c('qte','lincom') && length(p)>1) {
    if (is.null(PSI)) stop('For two-sample quantile linear combination treatment effect inference, need to provide weight vector PSI')
    if (length(PSI)!=length(p)) stop('Need to have length(PSI)==length(p)')
    if (PSI[1]!=1) stop('PSI should be normalized to have PSI[1]=1')
    if (any(PSI==0)) stop(sprint('PSI should not have any zero entries, but element #%d is zero',which(PSI==0)))
  }
  tmp <- length(X)
  if (METHOD.TYPE=='qte') tmp <- min(length(X$t),length(X$c))
  if (min(p)<1/tmp || max(p)>=1-1/tmp) stop(paste0('For this function, 1/n<=p<1-1/n is required.  ',
           'For values nearer to 0 or 1, ',
           'use methods based on extreme value theory.'))
  if (!(ONESIDED %in% c(-1,0,1))) stop(paste0('ONESIDED must be one of three values: 0 for two-sided inference;',
           ' -1 for one-sided with alternative hypothesis H1:xi<xi0',
           ' where confidence interval is of form (-Inf,b];',
           ' and 1 for one-sided with alternative hypothesis H1:xi>xi0',
           ' where confidence interval is of form [a,Inf).'))
  if (ONESIDED!=0 && !is.null(ALPHA) && ALPHA>=0.5) stop('ALPHA must be less than 0.5 for one-sided inference.')
  if (!is.logical(NORM.APPROX.CAL) || is.na(NORM.APPROX.CAL)) stop("NORM.APPROX.CAL must be either TRUE or FALSE")
  if (!is.logical(NORM.APPROX.Q) || is.na(NORM.APPROX.Q)) stop("NORM.APPROX.Q must be either TRUE or FALSE")
  if (!is.null(GAMMA)) {
    if (METHOD.TYPE=='qte') {
      if (length(GAMMA$t)!=length(GAMMA$c) || length(GAMMA$c)!=length(p)) {
        stop("For QTE method, GAMMA must be either NULL or a list with length(GAMMA$t)=length(GAMMA$c)=length(p).")
      }
    } else if (METHOD.TYPE %in% c('single','joint')) {
      warning("non-NULL GAMMA is being ignored since METHOD.TYPE is single or joint.")
    } else { #lincom
      if (length(GAMMA)!=length(p)) {
        stop("For linear combination method, GAMMA must be either NULL or numeric with length(GAMMA)=length(p).")
      }
    }
  }
}


#########################
#                       #
# Check numeric vector  #
#                       #
#########################
chk.num.vec.fn <- function(X) {
    varname <- deparse(substitute(X))
    if (!is.vector(X) || is.list(X)) stop(sprintf("%1$s must be a vector (and not a list).  If %1$s is a one-dimensional array or matrix, you can use as.vector(%1$s) in order to cast %1$s as a vector.",varname))
    if (!is.numeric(X)) stop(sprintf("%1$s must be numeric; i.e., is.numeric(%1$s) must be TRUE.",varname))
}


#########################
#                       #
# list[a,b,...] <- fn() #
#                       #
#########################
# Cool list[a,b,c] <- multivalued.fn(...) functionality
list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
   args <- as.list(match.call())
   args <- args[-c(1:2,length(args))]
   length(value) <- length(args)
   for(i in seq(along=args)) {
     a <- args[[i]]
     if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
   }
   x
}


#########################
#                       #
# GK single CI          #
#                       #
#########################
# Pre-condition: X is numeric vector sorted in ascending order
quantile.inf.single.GK <- function(X,p,ALPHA,ONESIDED,NORM.APPROX=FALSE,CALIB=FALSE) {
  a <- ifelse(ONESIDED!=0,ALPHA,ALPHA/2)
  lo <- -Inf;  hi <- Inf
  n <- length(X)
  if (ONESIDED>=0) {
    ul <- quantile.inf.CIul(p=p,n=n,a=a,NORM.APPROX=NORM.APPROX)
    if (CALIB) {
      z <- qnorm(1-a)
      epsl <- (n+1)*ul - floor((n+1)*ul)
      a.calib <- a + (epsl*(1-epsl)*z*exp(-z^2/2))/(sqrt(2*pi)*p*(1-p)*n)
      uc <- quantile.inf.CIul(p=p,n=n,a=a.calib,NORM.APPROX=NORM.APPROX)
      kc <- floor((n+1)*uc)
      if (kc>floor((n+1)*ul)) ul <- kc/(n+1) else ul <- uc
    }
    lo <- quantile.inf.interp(X=X,u=ul) 
  }
  if (ONESIDED<=0) {
    uh <- quantile.inf.CIuh(p=p,n=n,a=a,NORM.APPROX=NORM.APPROX)
    if (CALIB) {
      z <- qnorm(1-a)
      epsh <- (n+1)*uh - floor((n+1)*uh)
      a.calib <- a + (epsh*(1-epsh)*z*exp(-z^2/2))/(sqrt(2*pi)*p*(1-p)*n)
      uc <- quantile.inf.CIuh(p=p,n=n,a=a.calib,NORM.APPROX=NORM.APPROX)
      kc <- floor((n+1)*uc)
      if (kc<floor((n+1)*uh)) uh <- floor((n+1)*uh)/(n+1) else uh <- uc
    }
    hi <- quantile.inf.interp(X=X,u=uh) 
  }
  return(data.frame(lo=lo, hi=hi))
}


#########################
#                       #
# GK CI indices         #
#                       #
#########################
quantile.inf.CIul <- function(p,n,a,UNIROOT.INT=0:1+.Machine$double.eps*c(1,-1),UNIROOT.TOL=0.0001, APPROX=FALSE, NORM.APPROX=FALSE) {
  if (NORM.APPROX) { #takes precedent over APPROX; use high-order approx from Goldman and Kaplan lemma
    ul <- p - qnorm(1-a)*sqrt(p*(1-p)/n) - (2*p-1)*(qnorm(1-a)^2+2)/(6*n)
  } else if (APPROX) {  #APPROX (Hutson 1999, eqn 9): 
    ul <- qbeta(a,(n+1)*p,(n+1)*(1-p))
  } else {
    fun7 <- function(u1) (pbeta(p,(n+1)*u1,(n+1)*(1-u1)) - (1-a))
    ul <- uniroot(fun7,interval=UNIROOT.INT,tol=UNIROOT.TOL)$root
  }
  return(ul)
}
quantile.inf.CIuh <- function(p,n,a,UNIROOT.INT=0:1+.Machine$double.eps*c(1,-1),UNIROOT.TOL=0.0001, APPROX=FALSE, NORM.APPROX=FALSE) {
  if (NORM.APPROX) { #takes precedent over APPROX; use high-order approx from Goldman and Kaplan lemma
    uh <- p + qnorm(1-a)*sqrt(p*(1-p)/n) - (2*p-1)*(qnorm(1-a)^2+2)/(6*n)
  } else if (APPROX) {  #APPROX (Hutson 1999, eqn 9): 
    uh <- qbeta(1-a,(n+1)*p,(n+1)*(1-p))
  } else {
    fun8 <- function(u2) (pbeta(p,(n+1)*u2,(n+1)*(1-u2)) - a)
    uh <- uniroot(fun8,interval=UNIROOT.INT,tol=UNIROOT.TOL)$root
  }
  return(uh)
}


##########################
#                        #
# Quantile interpolation #
#                        #
##########################
# Pre-condition: X is numeric vector sorted in ascending order, u is between 0 and 1
quantile.inf.interp <- function(X,u) {
  n <- length(X)
  ind.u <- (n+1)*u
  ind.int <- floor(ind.u)
  eps <- ind.u - ind.int
  if (ind.int==0) stop("Order statistic does not exist: index 0")
  if (ind.int==n && eps>0) stop(sprintf("Order statistic does not exist: index %d=length(X)+1",n+1))
  return((1-eps)*X[ind.int]+eps*X[ind.int+1])
}


#########################
#                       #
# Kaplan single CI      #
#                       #
#########################
quantile.inf.single.K <- function(X,p,ALPHA,ONESIDED) {
  #Load critical values simulated ahead of time for small values of m.
  simcvs <- sub.simcv()
  simalphas <- simcvs$simalphas
  simCVs.uni <- simcvs$simCVs.uni
  #Initialize depending on ONESIDED.
  if (ONESIDED==-1) {
    CI.lo <- -Inf;  zz <- qnorm(1-ALPHA,0,1)
  } else if (ONESIDED==1) {
    CI.hi <- Inf;  zz <- qnorm(ALPHA,0,1)
  } else if (ONESIDED==0) {
    zz <- qnorm(1-ALPHA/2,0,1)
  } else {
    stop('Argument ONESIDED can only be -1, 0, or 1.')
  }
  nx <- length(X);  sqrtnx <- sqrt(nx);  sqrtp1p <- sqrt(p*(1-p))
  rx <- floor(nx*p)+1  #order statistic index for the sample quantile
  Xnr <- X[rx]  #order statistic, and estimator of quantile
  mmax <- min(rx-1,nx-rx)      #smoothing parameter m can't be larger than this
  if (mmax<1) stop(sprintf("Can't have p as close to zero or one as p=%4.2f is for n=%d",p,nx))
  #
  localpwr <- 0.5 #1st-order power of alternative against which to max power
  #OPTIMIZE IS FASTER THAN OPTIM, AND MORE ROBUST THAN UNIROOT
  CC <- optimize(f=function(x) (pnorm(x+zz)-pnorm(x-zz)-(1-localpwr))^2,
                 interval=c(.1,6), maximum=FALSE, tol=0.00001)$minimum
  m <- nx^(2/3)*(CC*zz)^(1/3)*(3/4)^(1/3) *
       (((dnorm(qnorm(p)))^2)/(2*(qnorm(p))^2+1))^(1/3) *
       ((dnorm(zz-CC)-dnorm(zz+CC))/(dnorm(zz-CC)+dnorm(zz+CC)))^(1/3)
  m <- max(1,min(mmax,floor(m))) #make sure not <1 or >mmax
  Smn <- (nx/(2*m)) * (X[rx+m]-X[rx-m])  #SBG sparsity estimator
  ccv <- abs(zz + (zz^3)/(4*m))  #+ (zz^5+8*zz^3)/(96*m^2); #2nd-order works fine, no need for 3rd-order
  #use absolute value; eqn holds for 1- or 2-sided (ONESIDED=-1,0,1)
  if (m<=dim(simCVs.uni)[1]) { #For small m, ccv approx is inaccurate=>use simulated ccv
    if (ONESIDED!=0) ALPHA <- 2*ALPHA
    if (sum(simalphas==ALPHA)>0) {
      ccv <- simCVs.uni[m,simalphas==ALPHA]
    } else {
      #INTERPOLATE BTWN CONSECUTIVE ALPHAS
      tmp <- which(simalphas>=ALPHA)[1]
      if (ALPHA>max(simalphas)) {
        ccv <- simCVs.uni[m,end]
      } else if (ALPHA<min(simalphas)) {
        ccv <- simCVs.uni[m,1]
      } else {
        tmp2 <- (simalphas[tmp]-ALPHA)/(simalphas[tmp]-simalphas[tmp-1])
        ccv <- tmp2*simCVs.uni[m,tmp-1] + (1-tmp2)*simCVs.uni[m,tmp]
      }
    }
  }
  #Usual calculation of CI endpoints given critical value & std error.
  if (ONESIDED>=0) CI.lo <- Xnr - ccv*Smn*sqrtp1p/sqrtnx
  if (ONESIDED<=0) CI.hi <- Xnr + ccv*Smn*sqrtp1p/sqrtnx
  return(data.frame(lo=CI.lo,hi=CI.hi))
}


#########################
#                       #
# Kaplan QTE CI         #
#                       #
#########################
#Pre-condition: X$t and X$c already sorted (ascending)
#SIM.CV=FALSE recommended since doesn't account for theta as well as ccv
quantile.inf.qte.K <- function(X,p,ALPHA,ONESIDED,PDFMIN=0.0001,KERNEL.TYPE='epanechnikov',BW.TYPE='SJ',PDFratio=NULL,SIM.CV=FALSE) {
  Y <- X$c;  X <- X$t
  #Load simulated critical values.
  if (SIM.CV) {
    simCVs.uni <- simCVs.bi <- NULL
    tryCatch({
      simcvs <- sub.simcv()
      simalphas<- simcvs$simalphas
      simCVs.uni <- simcvs$simCVs.uni
      simCVs.bi  <- simcvs$simCVs.bi
    }, error=function(ERR) {}) #warning("Could not load simulated fixed-m critical values inside function quantile.inf.qte.K(); using approximated CVs instead.")
  }
  nx <- length(X);  ny <- length(Y)
  rx <- floor(nx*p)+1; ry <- floor(ny*p)+1
  Xnr <- X[rx];  Ynr <- Y[ry] #sample quantiles for X, Y
  sqrtp1p <- sqrt(p*(1-p));  sqrtnx <- sqrt(nx)  #sqrtny=sqrt(ny);
  zz <- qnorm(1-ALPHA/2,0,1)  #standard normal critical value, 2-sided
  if (ONESIDED!=0) zz <- qnorm(1-ALPHA,0,1)
  localpwr <- 0.5 #1st-order power of alternative against which to max power
  CC <- optimize(f=function(x) (pnorm(x+zz)-pnorm(x-zz)-(1-localpwr))^2,
                 interval=c(.1,6), maximum=FALSE, tol=0.00001)$minimum
  #
  #Pilot estimate of theta to plug into m
  xeval <- X[rx];  yeval<-Y[ry] #where to evaluate PDFs when calculating theta
  if (is.null(PDFratio)) {
    fxest<-max(PDFMIN,
               density(X,from=xeval,to=xeval,n=1,kernel=KERNEL.TYPE,bw=BW.TYPE)$y)
    fyest<-max(PDFMIN,
               density(Y,from=yeval,to=yeval,n=1,kernel=KERNEL.TYPE,bw=BW.TYPE)$y)
    delta <- fxest/fyest
  } else {
    delta <- PDFratio
  }
#     theta.est <- (fxest^(-2)+(nx/ny)*fyest^(-2))^(-2) * (fxest^(-4)+(nx^2/ny^2)*fyest^(-4))
  theta.est <- (1+(nx/ny)*delta^2)^(-2) * (1+(nx^2/ny^2)*delta^4)
  #
  mmax <- min(min(rx-1,nx-rx),min(ry-1,ny-ry))
  m <- max(nx,ny)^(2/3)*(3/4)^(1/3) *
       ((1-theta.est)+theta.est*(CC*zz)*(dnorm(zz-CC)-dnorm(zz+CC)) /
                       (dnorm(zz-CC)+dnorm(zz+CC)))^(1/3) *
       (((dnorm(qnorm(p)))^2)/(2*(qnorm(p))^2+1))^(1/3)
  m <- max(1,min(mmax,floor(m))) #make sure not <1 or >mmax
  Sx <- (nx/(2*m)*(X[rx+m]-X[rx-m]))
  Sy <- (ny/(2*m)*(Y[ry+m]-Y[ry-m]))
  Smn <- sqrt(Sx^2+(nx/ny)*Sy^2)
  if (is.null(PDFratio)) {
    delta <- Sy/Sx #rem: S is sparsity, not density
  } else {
    delta <- PDFratio
  }
  theta.est <- (1+(nx/ny)*delta^2)^(-2) * (1+(nx^2/ny^2)*delta^4)
#         ((nx/(2*m)*(X[rx+m]-X[rx-m]))^4+(nx^2/ny^2)*(ny/(2*m)*(Y[ry+m]-Y[ry-m]))^4) /  (Smn^4)  #shrink a little toward 1 by adding (.1+...)/(.1+...)?
  c1 <- (1/4)*(theta.est*zz^3-zz*(1-theta.est))
  ccv1 <- zz + c1/m
  c2 <- (1/32)*(zz*(17*theta.est^2-30*theta.est+13) +zz^3*(16*theta.est^2-20*theta.est+(20/3)) +zz^5*(7*theta.est^2-10*theta.est+(10/3)))
  ccv2 <- zz + c1/m + c2/m^2
  ccv <- ccv2 #for 2-sample, sometimes ccv1>ccv2, like for theta=3/4; unlike 1s, can be >simCV, too
  #Using simulated CVs is somewhat ad hoc when nx!=ny
  if (SIM.CV && m<=dim(simCVs.uni)[1]) {
    if (ONESIDED!=0) ALPHA <- 2*ALPHA
    tmpw1 <- 2*theta.est-1; tmpw2 <- 2-2*theta.est #weight btwn uni and "bi" sim'd ccvs
    #rem: theta is between 1/2 and 1, hence not just tmpw1=theta
    if (sum(simalphas==ALPHA)>0)
       ccv <- tmpw1*simCVs.uni[m,simalphas==ALPHA] +
              tmpw2*simCVs.bi[m,simalphas==ALPHA]  #now vector, nreplic x 1
    else {
       #INTERPOLATE BTWN CONSECUTIVE ALPHAS
       tmp <- which(simalphas>=ALPHA)[1]
       if (ALPHA>max(simalphas))
           ccv <- tmpw1*simCVs.uni[m,end] +
                  tmpw2*simCVs.bi[m,end]
       else if (ALPHA<min(simalphas))
           ccv <- tmpw1*simCVs.uni[m,1] +
                  tmpw2*simCVs.bi[m,1]
       else {
           tmp2<-(simalphas[tmp]-ALPHA)/(simalphas[tmp]-simalphas[tmp-1])
           ccv <- tmpw1*(tmp2*simCVs.uni[m,tmp-1] + (1-tmp2)*simCVs.uni[m,tmp]) +
                  tmpw2*(tmp2*simCVs.bi[m,tmp-1] + (1-tmp2)*simCVs.bi[m,tmp])
       }
    }
  }
  #Usual CI endpoint calculation based on critical value, std error.
  CI.lo <- -Inf; CI.hi <- Inf
  if (ONESIDED>=0) CI.lo <- Xnr-Ynr - ccv*Smn*sqrtp1p/sqrtnx
  if (ONESIDED<=0) CI.hi <- Xnr-Ynr + ccv*Smn*sqrtp1p/sqrtnx
  return(data.frame(lo=CI.lo,hi=CI.hi))
}


#########################
#                       #
# Bootstrap QTE CI      #
#                       #
#########################
# Pre-condition: X is a list containing numeric vectors X$t and X$c, both of which are sorted in ascending order.
quantile.inf.qte.BSt <- function(X,p,ALPHA,ONESIDED) {
    BREP1<-99;  BREP2<-100 #Kaplan (2015) found these best in simulations.
    #BREP1 is outer replications, BREP2 is for calculating std deviation.
    nt <- length(X$t); nc <- length(X$c); sqrtnt <- sqrt(nt) #sqrtp1p=sqrt(p*(1-p));
    rx <- floor(nt*p)+1; ry <- floor(nc*p)+1
    Y <- X$c;  X <- X$t  #or rewrite following code to use X$t and X$c instead of X and Y
    Xnr <- X[rx]; Ynr <- Y[ry] #Sample quantiles for X and Y. Could also use quantile(X,p).
    #
    T0stub <- sqrtnt*(Xnr-Ynr)  #Test statistic "stub"
    #Use bootstrap to calculate standard error, to Studentize test
    #  statistic with.
    diffstars2 <- rep.int(0,BREP2) #matrix(0,BREP2,1)
    for (iBS2 in 1:BREP2) {
        starindx2 <- sample(1:nt,replace=TRUE) #from uniform integers 1:n, nx1 mtx
        starindy2 <- sample(1:nc,replace=TRUE) #from uniform integers 1:n, nx1 mtx
        Xstar2 <- X[starindx2];        Ystar2 <- Y[starindy2]
        Xstarsort2 <- sort(Xstar2);    Ystarsort2 <- sort(Ystar2)
        Xnrstar2 <- Xstarsort2[rx];     Ynrstar2 <- Ystarsort2[ry]
        diffstars2[iBS2] <- Xnrstar2-Ynrstar2
    }
    stdstar <- sd(sqrtnt*diffstars2)
    if (BREP2==0)  stdstar <- 1 #not Studentized
    BST0 <- T0stub / stdstar    #This is the test statistic.
    #
    #Now bootstrap to get dist'n of (pivotal) test statistic.
    Tstars <- rep.int(0,BREP1) #matrix(0,BREP1,1)
    for (iBS in 1:BREP1) {
        starindx <- sample(1:nt,replace=TRUE) #from uniform integers 1:n, nx1 mtx
        starindy <- sample(1:nc,replace=TRUE) #from uniform integers 1:n, nx1 mtx
        Xstar <- X[starindx];        Ystar <- Y[starindy]
        Xstarsort <- sort(Xstar);    Ystarsort <- sort(Ystar)
        Xnrstar <- Xstarsort[rx];     Ynrstar <- Ystarsort[ry]
        diffstars2 <- rep.int(0,BREP2) #matrix(0,BREP2,1)
        for (iBS2 in 1:BREP2) {
            starindx2 <- sample(1:nt,replace=TRUE) #from uniform integers 1:n, nx1 mtx
            starindy2 <- sample(1:nc,replace=TRUE) #from uniform integers 1:n, nx1 mtx
            Xstar2 <- Xstar[starindx2];       Ystar2 <- Ystar[starindy2]
            Xstarsort2 <- sort(Xstar2);       Ystarsort2 <- sort(Ystar2)
            Xnrsort2 <- Xstarsort2[rx];        Ynrsort2 <- Ystarsort2[ry]
            diffstars2[iBS2] <- Xnrsort2-Ynrsort2
        }
        stdstar<-sd(sqrtnt*diffstars2)
        if (BREP2==0)   stdstar<-1  #not Studentized
        Tstars[iBS] <- sqrtnt*((Xnrstar-Ynrstar)-(Xnr-Ynr))/(stdstar)
    }
    #
    #Take critical value from bootstrap distribution.
    Tstars.abssort <- sort(abs(Tstars)); Tstars.sort <- sort(Tstars)
    if (ONESIDED==0) cv<-Tstars.abssort[1+floor((1-ALPHA)*BREP1)] #symmetric performs best for 2-sided
    else if (ONESIDED==-1) cv<-Tstars.sort[1+floor(ALPHA*BREP1)]
    else if (ONESIDED== 1) cv<-Tstars.sort[1+floor((1-ALPHA)*BREP1)]
    #
    #Usual CI endpoint calculation given critical value and std error.
    CI.lo <- -Inf; CI.hi<-Inf
    if (ONESIDED>=0) CI.lo <- Xnr-Ynr - cv*stdstar/sqrtnt
    if (ONESIDED<=0) CI.hi <- Xnr-Ynr + abs(cv)*stdstar/sqrtnt
    return(data.frame(lo=CI.lo,hi=CI.hi))
}


#########################
#                       #
# Alpha calibration     #
#                       #
#########################
# Pre-condition: X (or X$t and X$c) is sorted in ascending order
quantile.inf.calibrate <- function(X,p,ALPHA,ONESIDED,METHOD.TYPE,PSI,NORM.APPROX,BETA.BLOCK.SIZE=10^5,BETA.BLOCKS=5,KERNEL.TYPE='gaussian',BW.TYPE='SJ',PDFMIN=0.0001,SPACING.FLAG=TRUE, GAMMA=NULL, NUM.INT=TRUE) {
  len.p <- length(p)
  # BETA.REPS <- BETA.BLOCKS*BETA.BLOCK.SIZE
  # Estimate PDF values
  if (METHOD.TYPE=='joint') { #no PDF estimation needed
    n <- length(X)
  } else if (METHOD.TYPE=='lincom') {
    n <- length(X)
  } else if (METHOD.TYPE=='qte') {
    n <- data.frame(t=length(X$t), c=length(X$c))
  } else { stop(sprintf("Function quantile.inf.calibrate() called with METHOD.TYPE=='%s' instead of 'qte', 'lincom', or 'joint'",METHOD.TYPE)) }
  
  # Normal approximation
  if (NORM.APPROX) {
    if (METHOD.TYPE!='qte' || len.p>1) { 
      warning(sprintf("Normal approximation currently only available for METHOD.TYPE 'qte' with length(p)==1; called with METHOD.TYPE=='%s' and length(p)==%d",METHOD.TYPE,len.p))
    } else {
      if (is.null(GAMMA)) {
        delta.est <- (sqrt(n$t)*quantile.inf.density(X=X$t,p=p,n=n$t, SPACING.FLAG=SPACING.FLAG))/(sqrt(n$c)*quantile.inf.density(X=X$c,p=p,n=n$c, SPACING.FLAG=SPACING.FLAG))
      } else {
        delta.est <- sqrt(n$t) / sqrt(n$c) * GAMMA$c
      }
      thetastar.est <- (1+delta.est)/sqrt(1+delta.est^2)
      tmp <- 2-abs(ONESIDED) #2 for 2-sided, 1 for 1-sided
      return(tmp*pnorm(qnorm(ALPHA/tmp)/thetastar.est))
    }
  }
  #
  if (METHOD.TYPE=='qte') { #search over interval (ALPHA,2*pnorm(qnorm(ALPHA/2)/sqrt(2*sum(abs(PSI))))) if ONESIDED==0, else w/o two 2's
    # Estimate \gamma_j=f(Q(u_1))/f(Q(u_j))
    if (!is.null(GAMMA)) {
      gammat <- GAMMA$t;  gammac <- GAMMA$c
    } else {
      gammat <- gammac <- numeric(len.p)
      f1 <- quantile.inf.density(X=X$t,p=p[1],n=n$t, SPACING.FLAG=SPACING.FLAG)
      for (j in 1:len.p) {
        gammat[j] <- f1 / quantile.inf.density(X=X$t,p=p[j],n=n$t, SPACING.FLAG=SPACING.FLAG)
        gammac[j] <- f1 / quantile.inf.density(X=X$c,p=p[j],n=n$c, SPACING.FLAG=SPACING.FLAG)
      }
    }
    # Generate calibration function
    comp.val <- sum(PSI*p*(gammat-gammac))
    if (ONESIDED!=0) {
      if (NUM.INT && len.p==1) {
        cal.fn <- function(a) {
          RELTOL <- 0.0001 #0.001: if ALPHA=0.10, within 0.0001. but not the most reliable request, so made smaller. (default is .Machine$double.eps^0.25, around 0.00012 on mine)
          #eqn 6 refers to Hutson (1999) equation (6)
          if (ONESIDED<0) {
            uth <- quantile.inf.CIuh(p=p,n=n$t,a=a)
            ucl <- quantile.inf.CIul(p=p,n=n$c,a=a)
            if (gammac==1) {
              eqn6 <- integrate(function(y){dbeta(y,(n$c+1)*ucl,(n$c+1)*(1-ucl))*
                                            pbeta(y,(n$t+1)*uth,(n$t+1)*(1-uth))},
                                max(0,ucl-3/sqrt(n$c)), min(1,ucl+3/sqrt(n$c)),
                                rel.tol=RELTOL)$value
            } else {
              eqn6 <- integrate(function(y){dbeta(y,(n$c+1)*ucl,(n$c+1)*(1-ucl))*
                                            pbeta(y*gammac+p*(1-gammac),(n$t+1)*uth,(n$t+1)*(1-uth))},
                                max(0,ucl-3/sqrt(n$c)), min(1,ucl+3/sqrt(n$c)),
                                rel.tol=RELTOL)$value
            }
          } else {
            utl <- quantile.inf.CIul(p=p,n=n$t,a=a)
            uch <- quantile.inf.CIuh(p=p,n=n$c,a=a)
            if (gammac==1) {
              eqn6 <- integrate(function(x){dbeta(x,(n$t+1)*utl,(n$t+1)*(1-utl))*
                                            pbeta(x,(n$c+1)*uch,(n$c+1)*(1-uch))},
                                max(0,utl-3/sqrt(n$t)), min(1,utl+3/sqrt(n$t)),
                                rel.tol=RELTOL)$value
            } else {
              eqn6 <- integrate(function(x){dbeta(x,(n$t+1)*utl,(n$t+1)*(1-utl))*
                                            pbeta(x/gammac+p*(1-1/gammac),(n$c+1)*uch,(n$c+1)*(1-uch))},
                                max(0,utl-3/sqrt(n$t)), min(1,utl+3/sqrt(n$t)),
                                rel.tol=RELTOL)$value
            }
          }
          return((1 - eqn6) - (1-ALPHA))
        }
      } else {
        cal.fn <- function(a) {
          ust <- usc <- rep.int(NA,len.p)
          for (i in 1:len.p) {
            if (ONESIDED*PSI[i]==-1) {
              ust[i] <- quantile.inf.CIuh(p=p[i],n=n$t,a=a)
              usc[i] <- quantile.inf.CIul(p=p[i],n=n$c,a=a)
            } else {
              ust[i] <- quantile.inf.CIul(p=p[i],n=n$t,a=a)
              usc[i] <- quantile.inf.CIuh(p=p[i],n=n$c,a=a)
            }
          }
          blk.cps <- rep.int(NA,BETA.BLOCKS)
          for (i in 1:BETA.BLOCKS) {
            Bst <- quantile.inf.betas(u=sort(ust),n=n$t, reps=BETA.BLOCK.SIZE)
            Bst[,order(ust)] <- Bst
            Bsc <- quantile.inf.betas(u=sort(usc),n=n$c, reps=BETA.BLOCK.SIZE)
            Bsc[,order(usc)] <- Bsc
            tmpind <- 1:len.p
            tmpt <- PSI[1]*gammat[1]*Bst[,1]
            tmpc <- PSI[1]*gammac[1]*Bsc[,1]
            for (j in tmpind[-1]) {
              tmpt <- tmpt + PSI[j]*gammat[j]*Bst[,j]
              tmpc <- tmpc + PSI[j]*gammac[j]*Bsc[,j]
            }
            blk.cps[i] <- mean(ONESIDED*(tmpt-tmpc)<ONESIDED*comp.val)
          }
          return(mean(blk.cps) - (1-ALPHA))
        }
      }
    } else { #two-sided
      if (NUM.INT && len.p==1) {
        cal.fn <- function(a) {
          utl <- quantile.inf.CIul(p=p,n=n$t,a=a/2)
          uth <- quantile.inf.CIuh(p=p,n=n$t,a=a/2)
          ucl <- quantile.inf.CIul(p=p,n=n$c,a=a/2)
          uch <- quantile.inf.CIuh(p=p,n=n$c,a=a/2)
          RELTOL <- 0.0001 #0.001: if ALPHA=0.10, within 0.0001. but not the most reliable request, so made smaller. (default is .Machine$double.eps^0.25, around 0.00012 on mine)
          #eqn 6 refers to Hutson (1999) equation (6)
          if (gammac==1) {
            eqn6_tlch <- integrate(function(x){dbeta(x,(n$t+1)*utl,(n$t+1)*(1-utl))*
                                               pbeta(x,(n$c+1)*uch,(n$c+1)*(1-uch))},
                                   max(0,utl-3/sqrt(n$t)), min(1,utl+3/sqrt(n$t)),
                                   rel.tol=RELTOL)$value
          } else {
            eqn6_tlch <- integrate(function(x){dbeta(x,(n$t+1)*utl,(n$t+1)*(1-utl))*
                                               pbeta(x/gammac+p*(1-1/gammac),(n$c+1)*uch,(n$c+1)*(1-uch))},
                                   max(0,utl-3/sqrt(n$t)), min(1,utl+3/sqrt(n$t)),
                                   rel.tol=RELTOL)$value
          }
          if (gammac==1) {
            eqn6_clth <- integrate(function(y){dbeta(y,(n$c+1)*ucl,(n$c+1)*(1-ucl))*
                                               pbeta(y,(n$t+1)*uth,(n$t+1)*(1-uth))},
                                   max(0,ucl-3/sqrt(n$c)), min(1,ucl+3/sqrt(n$c)),
                                   rel.tol=RELTOL)$value
          } else {
            eqn6_clth <- integrate(function(y){dbeta(y,(n$c+1)*ucl,(n$c+1)*(1-ucl))*
                                               pbeta(y*gammac+p*(1-gammac),(n$t+1)*uth,(n$t+1)*(1-uth))},
                                   max(0,ucl-3/sqrt(n$c)), min(1,ucl+3/sqrt(n$c)),
                                   rel.tol=RELTOL)$value
          }
          return((1 - eqn6_tlch - eqn6_clth) - (1-ALPHA))
        }
      } else {
        cal.fn <- function(a) {
          ust <- usc <- rep.int(NA,2*len.p) #lower endpoint first
          for (i in 1:len.p) {
            ust[i] <- quantile.inf.CIul(p=p[i],n=n$t,a=a/2)
            ust[len.p+i] <- quantile.inf.CIuh(p=p[i],n=n$t,a=a/2)
            usc[i] <- quantile.inf.CIul(p=p[i],n=n$c,a=a/2)
            usc[len.p+i] <- quantile.inf.CIuh(p=p[i],n=n$c,a=a/2)
            if (PSI[i]<0) {
              ust[c(i,len.p+i)] <- ust[c(len.p+i,i)]
              usc[c(i,len.p+i)] <- usc[c(len.p+i,i)]
            }
          }
          blk.cps <- rep.int(NA,BETA.BLOCKS)
          for (i in 1:BETA.BLOCKS) {
            Bst <- quantile.inf.betas(u=sort(ust),n=n$t, reps=BETA.BLOCK.SIZE)
            Bst[,order(ust)] <- Bst
            Bsc <- quantile.inf.betas(u=sort(usc),n=n$c, reps=BETA.BLOCK.SIZE)
            Bsc[,order(usc)] <- Bsc
            tmpind <- 1:len.p
            tmpLt <- PSI[1]*gammat[1]*Bst[,1]
            tmpHt <- PSI[1]*gammat[1]*Bst[,1+len.p]
            tmpLc <- PSI[1]*gammac[1]*Bsc[,1]
            tmpHc <- PSI[1]*gammac[1]*Bsc[,1+len.p]
            for (j in tmpind[-1]) {
              tmpLt <- tmpLt + PSI[j]*gammat[j]*Bst[,j]
              tmpHt <- tmpHt + PSI[j]*gammat[j]*Bst[,j+len.p]
              tmpLc <- tmpLc + PSI[j]*gammac[j]*Bsc[,j]
              tmpHc <- tmpHc + PSI[j]*gammac[j]*Bsc[,j+len.p]
            }
            blk.cps[i] <- mean((tmpLt-tmpHc<comp.val) & (tmpHt-tmpLc>comp.val))
          }
          return(mean(blk.cps) - (1-ALPHA))
        }
      }
    }
    tmp <- cal.fn(0.99)
    if (tmp>0) {
      ALPHAtilde <- 0.99
    } else {
      ALPHAtilde <- uniroot(f=cal.fn, interval=c(ALPHA,0.99), tol=max(0.00001,min(0.001,min(n$t,n$c)^(-3/2))), f.upper=tmp)$root
    }
  } else if (METHOD.TYPE=='lincom') { #search over interval (ALPHA,1)
    # Estimate \gamma_j=f(Q(u_1))/f(Q(u_j))
    if (!is.null(GAMMA)) {
      gamma <- GAMMA
    } else {
      gamma <- numeric(len.p);  gamma[1] <- 1
      f1 <- quantile.inf.density(X=X,p=p[1],n=n, SPACING.FLAG=SPACING.FLAG)
      for (j in 2:len.p) {
        gamma[j] <- f1 / quantile.inf.density(X=X,p=p[j],n=n, SPACING.FLAG=SPACING.FLAG)
      }
    }
    # Generate calibration function
    comp.val <- sum(PSI*gamma*p)
    if (ONESIDED!=0) {
      cal.fn <- function(a) {
        us <- rep.int(NA,len.p)
        for (i in 1:len.p) {
          if (ONESIDED*sign(PSI[i])==-1) { us[i] <- quantile.inf.CIuh(p=p[i],n=n,a=a) } else { us[i] <- quantile.inf.CIul(p=p[i],n=n,a=a) }
        }
        blk.cps <- rep.int(NA,BETA.BLOCKS)
        for (i in 1:BETA.BLOCKS) {
          Bs <- quantile.inf.betas(u=sort(us),n=n,reps=BETA.BLOCK.SIZE)
          Bs[,order(us)] <- Bs
          tmp <- PSI[1]*gamma[1]*Bs[,1]
          tmp2 <- 1:len.p
          for (j in tmp2[-1]) { tmp <- tmp + PSI[j]*gamma[j]*Bs[,j] }
          blk.cps[i] <- mean(ONESIDED*tmp<ONESIDED*comp.val)
        }
        return(mean(blk.cps) - (1-ALPHA))
      }
    } else { #two-sided
      cal.fn <- function(a) {
        us <- rep.int(NA,2*len.p) #lower endpoint first
        for (i in 1:len.p) {
          us[i] <- quantile.inf.CIul(p=p[i],n=n,a=a/2)
          us[len.p+i] <- quantile.inf.CIuh(p=p[i],n=n,a=a/2)
          if (PSI[i]<0) us[c(i,len.p+i)] <- us[c(len.p+i,i)]
        }
        blk.cps <- rep.int(NA,BETA.BLOCKS)
        for (i in 1:BETA.BLOCKS) {
          Bs <- quantile.inf.betas(u=sort(us),n=n,reps=BETA.BLOCK.SIZE)
          Bs[,order(us)] <- Bs
          tmpind <- 1:len.p
          tmpL <- PSI[1]*gamma[1]*Bs[,1]
          for (j in tmpind[-1]) { tmpL <- tmpL + PSI[j]*gamma[j]*Bs[,j] }
          tmpH <- PSI[1]*gamma[1]*Bs[,len.p+1]
          for (j in tmpind[-1]) { tmpH <- tmpH + PSI[j]*gamma[j]*Bs[,len.p+j] }
          blk.cps[i] <- mean((tmpL<comp.val) & (tmpH>comp.val))
        }
        return(mean(blk.cps) - (1-ALPHA))
      }
    }
    tmp <- cal.fn(0.99)
    if (tmp>0) {
      ALPHAtilde <- 0.99
    } else {
      ALPHAtilde <- uniroot(f=cal.fn, interval=c(ALPHA,0.99), tol=max(0.00001,min(0.001,n^(-3/2))), f.upper=tmp)$root
    }
  } else if (METHOD.TYPE=='joint') { #search over interval (0,ALPHA)
  	comp.mat <- NULL;  for (i in 1:len.p) { comp.mat <- cbind(comp.mat,rep.int(p[i],BETA.BLOCK.SIZE)) }
  	if (ONESIDED!=0) {
      cal.fn <- function(a) {
        us <- rep.int(NA,len.p)
        for (i in 1:len.p) { 
          if (ONESIDED==-1) { us[i] <- quantile.inf.CIuh(p=p[i],n=n,a=a) } else { us[i] <- quantile.inf.CIul(p=p[i],n=n,a=a) }
        }
        blk.cps <- rep.int(NA,BETA.BLOCKS)
        for (i in 1:BETA.BLOCKS) {
          Bs <- quantile.inf.betas(u=sort(us),n=n,reps=BETA.BLOCK.SIZE)
          Bs[,order(us)] <- Bs
          if (ONESIDED==-1) {
            blk.cps[i] <- mean(prod.matrix(Bs>comp.mat)) #percentage of rows/draws where *all* Bs>us
          } else {
            blk.cps[i] <- mean(prod.matrix(Bs<comp.mat)) #want lower endpoint below true quantiles
          }
        }
        return(mean(blk.cps) - (1-ALPHA))
      }
    } else { #two-sided CI: both upper and lower endpoints
      cal.fn <- function(a) {
        uls <- uhs <- rep.int(NA,len.p)
        for (i in 1:len.p) { 
          uls[i] <- quantile.inf.CIul(p=p[i],n=n,a=a/2)
          uhs[i] <- quantile.inf.CIuh(p=p[i],n=n,a=a/2)
        }
        us <- c(uls,uhs)
        blk.cps <- rep.int(NA,BETA.BLOCKS)
        for (i in 1:BETA.BLOCKS) {
          Bs <- quantile.inf.betas(u=sort(us),n=n,reps=BETA.BLOCK.SIZE)
          Bs[,order(us)] <- Bs
          blk.cps[i] <- mean(prod.matrix(Bs[,1:len.p]<comp.mat)*prod.matrix(Bs[,-(1:len.p)]>comp.mat)) #percentage of rows/draws where *all* Bls<uls and Bhs>uhs
        }
        return(mean(blk.cps) - (1-ALPHA))
      }
    }
    ALPHAtilde <- uniroot(f=cal.fn, interval=c(0,ALPHA), f.lower=ALPHA, tol=max(0.00001,min(0.001,n^(-3/2))))$root
  }
  #
  return(ALPHAtilde)
}


#########################
#                       #
# Beta draws            #
#                       #
#########################
quantile.inf.betas <- function(u,n,reps) {
  len.u <- length(u)
  ret <- matrix(nrow=reps,ncol=len.u)
  #
  # rbeta alternative
  ret[,1] <- rbeta(n=reps,shape1=(n+1)*u[1],shape2=(n+1)*(1-u[1]))
  tmp <- 1:len.u
  for (i in tmp[-1]) {
    ret[,i] <- ret[,i-1] + (1-ret[,i-1])*rbeta(n=reps,shape1=(n+1)*(u[i]-u[i-1]),shape2=(n+1)*(1-u[i]))
  }
  return(ret)
  # #
  # # bayesm alternative: matches, but much slower
  # library(bayesm)
  # dir.params <- (n+1)*u[1]
  # tmp <- 1:len.u
  # for (i in tmp[-1]) { dir.params <- c(dir.params,(n+1)*(u[i]-u[i-1])) }
  # dir.params <- c(dir.params,(n+1)*(1-u[len.u]))
  # for (i in 1:reps) {
    # ret[i,] <- cumsum(rdirichlet(dir.params))[-(len.u+1)]
  # }
  # return(ret)
}


#########################
#                       #
# Density estimation    #
#                       #
#########################
# Pre-condition: X is a numeric vector sorted in ascending order, p is scalar between 0 and 1, n=length(X)
quantile.inf.bias0.bandwidth <- function(u,n)  min(u-1/(n+1),  min(n/(n+1)-u,  n^(2/3) *((3*dnorm(qnorm(u))^2)/(2+4*qnorm(u)^2))^(1/3) /(n+1) ) )
quantile.inf.density <- function(X,p,n,KERNEL.TYPE='gaussian',BW.TYPE='SJ',PDFMIN=0.0001,SPACING.FLAG=FALSE) {
  if (SPACING.FLAG) { #quantile spacing estimator
    # Compute bandwidth
    tmp <- quantile.inf.bias0.bandwidth(p,n)
    bias0.band <- min(tmp, p-1/(n+1)-2*.Machine$double.eps, n/(n+1)-p-2*.Machine$double.eps)
    # Compute spacing
    sparsity.est <- (1/(2*bias0.band))*(quantile.inf.interp(X,p+bias0.band)-quantile.inf.interp(X,p-bias0.band))
    return(1/sparsity.est)
  } else { #kernel estimator
    xeval <- quantile.inf.interp(X,p)
    fxest <- density(X, from=xeval, to=xeval, n=1, kernel=KERNEL.TYPE, bw=BW.TYPE)$y
    fxest <- max(fxest,PDFMIN)
  }
}


#########################
#                       #
# Row product           #
#                       #
#########################
prod.matrix <- function(x) {
  y <- x[,1]
  tmp <- 1:dim(x)[2]
  for (i in tmp[-1]) y <- y*x[,i]
  return(y)
}


## EOF
