#' @title Converts trial by trial experimental information into response counts.
#'
#' @description Given data from an experiment where an observer discriminates between two stimulus alternatives on every trial and provides confidence ratings, converts trial by trial experimental information for N trials into response counts.
#' @param stimID 1xN vector. stimID[i] = 0 --> stimulus on i'th trial was S1; stimID[i] = 1 --> stimulus on i'th trial was S2.
#' @param response 1xN vector. response[i] = 0 --> response on i'th trial was "S1"; response[i] = 1 --> response on i'th trial was "S2".
#' @param rating 1xN vector. rating(i) = X --> rating on i'th trial was X; X must be in the range 1 <= X <= nRatings.
#' @param nRatings A number. Total number of available subjective ratings available for the subject. e.g. if subject can rate confidence on a scale of 1-4, then nRatings = 4.
#' @param padCells If set to TRUE, each response count in the output has the value of padAmount added to it. Padding cells is desirable if trial counts of 0 interfere with model fitting. If set to FALSE, trial counts are not manipulated and 0s may be present in the response count output.
#' @param padAmount The value to add to each response count if padCells == TRUE. Default value is 1/(2*nRatings).
#' @return A list contains two vectors, nR_S1 and nR_S2. These are vectors containing the total number of responses in each response category, conditional on presentation of S1 and S2.
#'
#'     If nR_S1 = c(100,50,20,10,5,1), then when stimulus S1 was presented, the subject had the following response counts:
#'
#'     responded S1, rating=3 : 100 times
#'
#'     responded S1, rating=2 : 50 times
#'
#'     responded S1, rating=1 : 20 times
#'
#'     responded S2, rating=1 : 10 times
#'
#'     responded S2, rating=2 : 5 times
#'
#'     responded S2, rating=3 : 1 time
#' @export

trials2counts <- function(stimID, response, rating, nRatings, padCells = FALSE, padAmount = 1 / (2*nRatings)){

  # trials2counts(stimID, response, rating, nRatings, padCells, padAmount)
  #
  # Given data from an experiment where an observer discriminates between two
  # stimulus alternatives on every trial and provides confidence ratings,
  # converts trial by trial experimental information for N trials into response
  # counts.
  #
  # INPUTS
  # stimID:   1xN vector. stimID(i) = 0 --> stimulus on i'th trial was S1.
  #                       stimID(i) = 1 --> stimulus on i'th trial was S2.
  #
  # response: 1xN vector. response(i) = 0 --> response on i'th trial was "S1".
  #                       response(i) = 1 --> response on i'th trial was "S2".
  #
  # rating:   1xN vector. rating(i) = X --> rating on i'th trial was X.
  #                       X must be in the range 1 <= X <= nRatings.
  #
  # N.B. all trials where stimID is not 0 or 1, response is not 0 or 1, or
  # rating is not in the range [1, nRatings], are omitted from the response
  # count.
  #
  # nRatings: total # of available subjective ratings available for the
  #           subject. e.g. if subject can rate confidence on a scale of 1-4,
  #           then nRatings = 4
  #
  # optional inputs
  #
  # padCells: if set to 1, each response count in the output has the value of
  #           padAmount added to it. Padding cells is desirable if trial counts
  #           of 0 interfere with model fitting.
  #           if set to 0, trial counts are not manipulated and 0s may be
  #           present in the response count output.
  #           default value for padCells is 0.
  #
  # padAmount: the value to add to each response count if padCells is set to 1.
  #            default value is 1/(2*nRatings)
  #
  #
  # OUTPUTS
  # list(nR_S1, nR_S2)
  # these are vectors containing the total number of responses in
  # each response category, conditional on presentation of S1 and S2.
  #
  # e.g. if nR_S1 = c(100 50 20 10 5 1), then when stimulus S1 was
  # presented, the subject had the following response counts:
  # responded S1, rating=3 : 100 times
  # responded S1, rating=2 : 50 times
  # responded S1, rating=1 : 20 times
  # responded S2, rating=1 : 10 times
  # responded S2, rating=2 : 5 times
  # responded S2, rating=3 : 1 time
  #
  # The ordering of response / rating counts for S2 should be the same as it
  # is for S1. e.g. if nR_S2 = c(3 7 8 12 27 89), then when stimulus S2 was
  # presented, the subject had the following response counts:
  # responded S1, rating=3 : 3 times
  # responded S1, rating=2 : 7 times
  # responded S1, rating=1 : 8 times
  # responded S2, rating=1 : 12 times
  # responded S2, rating=2 : 27 times
  # responded S2, rating=3 : 89 times


  ## sort inputs

  # check for valid inputs
  if (!( length(stimID) == length(response) && length(stimID) == length(rating) )){
    print("stimID, response, and rating input vectors must have the same lengths")
  }

  # filter bad trials
  f <- (stimID == 0 | stimID == 1) & (response == 0 | response == 1) & (rating >=1 & rating <= nRatings);
  stimID   <- stimID[f];
  response <- response[f];
  rating   <- rating[f];

  ## compute response counts

  nR_S1 <- array(dim = 2 * nRatings);
  nR_S2 <- array(dim = 2 * nRatings);

  # S1 responses
  for (r in nRatings : 1){
    nR_S1[nRatings - r + 1] <- sum(stimID==0 & response==0 & rating==r);
    nR_S2[nRatings - r + 1] <- sum(stimID==1 & response==0 & rating==r);
  }

  # S2 responses
  for (r in 1 : nRatings){
    nR_S1[nRatings + r] <- sum(stimID==0 & response==1 & rating==r);
    nR_S2[nRatings + r] <- sum(stimID==1 & response==1 & rating==r);
  }


  # pad response counts to avoid zeros
  if(padCells){
    nR_S1 <- nR_S1 + padAmount;
    nR_S2 <- nR_S2 + padAmount;
  }

  return(list("nR_S1" = nR_S1, "nR_S2" = nR_S2))
}

#' @title Calculate meta-d' using Maximum likelihood estimation.
#'
#' @description Given data from an experiment where an observer discriminates between two stimulus alternatives on every trial and provides confidence ratings, provides a type 2 SDT analysis of the data.
#' @param nR_S1 A vector. Contains the total number of responses in the first response category, conditional on presentation of S1 and S2.
#'
#' If nR_S1 = c(100,50,20,10,5,1), then when stimulus S1 was presented, the subject had the following response counts:
#'
#' responded S1, rating=3 : 100 times
#'
#' responded S1, rating=2 : 50 times
#'
#' responded S1, rating=1 : 20 times
#'
#' responded S2, rating=1 : 10 times
#'
#' responded S2, rating=2 : 5 times
#'
#' responded S2, rating=3 : 1 time
#' @param nR_S2 A vector. Contains the total number of responses in the second response category, conditional on presentation of S1 and S2.
#'
#' The ordering of response / rating counts for S2 should be the same as it is for S1.
#'
#' If nR_S2 = c(3,7,8,12,27,89), then when stimulus S2 was presented, the subject had the following response counts:
#'
#' responded S1, rating=3 : 3 times
#'
#' responded S1, rating=2 : 7 times
#'
#' responded S1, rating=1 : 8 times
#'
#' responded S2, rating=1 : 12 times
#'
#' responded S2, rating=2 : 27 times
#'
#' responded S2, rating=3 : 89 times
#' @param s A number. The ratio of standard deviations for type 1 distributions, i.e. s = sd(S1) / sd(S2). If not specified, s is set to a default value of 1.
#' @param fncdf A function for the CDF of the type 1 distribution. If not specified, fncdf defaults to \code{pnorm} (i.e. CDF for normal distribution).
#' @param fninv A function for the inverse CDF of the type 1 distribution. If not specified, fninv defaults to \code{qnorm} (i.e. inverse CDF for normal distribution).
#' @return A list. In the following, let S1 and S2 represent the distributions of evidence
#' generated by stimulus classes S1 and S2. Then the fields is a list are as follows:
#'
#' da           = mean(S2) - mean(S1), in root-mean-square(sd(S1),sd(S2)) units
#'
#' s            = sd(S1) / sd(S2)
#'
#' meta_da      = meta-d' in RMS units
#'
#' M_diff       = meta_da - da
#'
#' M_ratio      = meta_da / da
#'
#' meta_ca      = type 1 criterion for meta-d' fit, RMS units
#'
#' t2ca_rS1     = type 2 criteria of "S1" responses for meta-d' fit, RMS units
#'
#' t2ca_rS2     = type 2 criteria of "S2" responses for meta-d' fit, RMS units
#'
#' S1units      = contains same parameters in sd(S1) units.
#'
#' logL         = log likelihood of the data fit
#'
#' est_HR2_rS1  = estimated (from meta-d' fit) type 2 hit rates for S1 responses
#'
#' obs_HR2_rS1  = actual type 2 hit rates for S1 responses
#'
#' est_FAR2_rS1 = estimated type 2 false alarm rates for S1 responses
#'
#' obs_FAR2_rS1 = actual type 2 false alarm rates for S1 responses
#'
#' est_HR2_rS2  = estimated type 2 hit rates for S2 responses
#'
#' obs_HR2_rS2  = actual type 2 hit rates for S2 responses
#'
#' est_FAR2_rS2 = estimated type 2 false alarm rates for S2 responses
#'
#' obs_FAR2_rS2 = actual type 2 false alarm rates for S2 responses
#'
#' If there are N ratings, then there will be N-1 type 2 hit rates and falsealarm rates.
#' @export

fit_meta_d_MLE <- function(nR_S1, nR_S2, s = 1, fncdf = pnorm, fninv = qnorm){

#
# fit_meta_d_MLE(nR_S1, nR_S2, s, fncdf, fninv)
#
# Given data from an experiment where an observer discriminates between two
# stimulus alternatives on every trial and provides confidence ratings,
# provides a type 2 SDT analysis of the data.
#
# INPUTS
#
# * nR_S1, nR_S2
# these are vectors containing the total number of responses in
# each response category, conditional on presentation of S1 and S2.
#
# e.g. if nR_S1 = c(100 50 20 10 5 1), then when stimulus S1 was
# presented, the subject had the following response counts:
  # responded S1, rating=3 : 100 times
  # responded S1, rating=2 : 50 times
  # responded S1, rating=1 : 20 times
  # responded S2, rating=1 : 10 times
  # responded S2, rating=2 : 5 times
  # responded S2, rating=3 : 1 time
  #
  # The ordering of response / rating counts for S2 should be the same as it
  # is for S1. e.g. if nR_S2 = c(3 7 8 12 27 89), then when stimulus S2 was
  # presented, the subject had the following response counts:
    # responded S1, rating=3 : 3 times
    # responded S1, rating=2 : 7 times
    # responded S1, rating=1 : 8 times
    # responded S2, rating=1 : 12 times
    # responded S2, rating=2 : 27 times
    # responded S2, rating=3 : 89 times
    #
    # N.B. if nR_S1 or nR_S2 contain zeros, this may interfere with estimation of
    # meta-d'.
    #
    # Some options for dealing with response cell counts containing zeros are:
    #
    # (1) Add a small adjustment factor, e.g. adj_f = 1/(length(nR_S1), to each
    # input vector:
    #
    # adj_f = 1/length(nR_S1);
    # nR_S1_adj = nR_S1 + adj_f;
    # nR_S2_adj = nR_S2 + adj_f;
    #
    # This is a generalization of the correction for similar estimation issues of
    # type 1 d' as recommended in
    #
    # Hautus, M. J. (1995). Corrections for extreme proportions and their biasing
    #     effects on estimated values of d'. Behavior Research Methods, Instruments,
    #     & Computers, 27, 46-51.
    #
    # When using this correction method, it is recommended to add the adjustment
    # factor to ALL data for all subjects, even for those subjects whose data is
    # not in need of such correction, in order to avoid biases in the analysis
    # (cf Snodgrass & Corwin, 1988).
    #
    # (2) Collapse across rating categories.
    #
    # e.g. if your data set has 4 possible confidence ratings such that length(nR_S1)==8,
    # defining new input vectors
    #
    # nR_S1_new = c(sum(nR_S1(1:2)), sum(nR_S1(3:4)), sum(nR_S1(5:6)), sum(nR_S1(7:8)));
    # nR_S2_new = c(sum(nR_S2(1:2)), sum(nR_S2(3:4)), sum(nR_S2(5:6)), sum(nR_S2(7:8)));
    #
    # might be sufficient to eliminate zeros from the input without using an adjustment.
    #
    # * s
    # this is the ratio of standard deviations for type 1 distributions, i.e.
    #
    # s = sd(S1) / sd(S2)
    #
    # if not specified, s is set to a default value of 1.
    # For most purposes, we recommend setting s = 1.
    # See http://www.columbia.edu/~bsm2105/type2sdt for further discussion.
    #
    # * fncdf
    # a function handle for the CDF of the type 1 distribution.
    # if not specified, fncdf defaults to @normcdf (i.e. CDF for normal
    # distribution)
    #
    # * fninv
    # a function handle for the inverse CDF of the type 1 distribution.
    # if not specified, fninv defaults to @norminv
    #
    # OUTPUT
    #
    # Output is packaged in the struct "fit."
    # In the following, let S1 and S2 represent the distributions of evidence
    # generated by stimulus classes S1 and S2.
    # Then the fields is a list are as follows:
    #
    # da           = mean(S2) - mean(S1), in room-mean-square(sd(S1),sd(S2)) units
    # s            = sd(S1) / sd(S2)
    # meta_da      = meta-d' in RMS units
    # M_diff       = meta_da - da
    # M_ratio      = meta_da / da
    # meta_ca      = type 1 criterion for meta-d' fit, RMS units
    # t2ca_rS1     = type 2 criteria of "S1" responses for meta-d' fit, RMS units
    # t2ca_rS2     = type 2 criteria of "S2" responses for meta-d' fit, RMS units
    #
    # S1units      = contains same parameters in sd(S1) units.
    #                 these may be of use since the data-fitting is conducted
    #                 using parameters specified in sd(S1) units.
    #
    # logL         = log likelihood of the data fit
    #
    # est_HR2_rS1  = estimated (from meta-d' fit) type 2 hit rates for S1 responses
    # obs_HR2_rS1  = actual type 2 hit rates for S1 responses
    # est_FAR2_rS1 = estimated type 2 false alarm rates for S1 responses
    # obs_FAR2_rS1 = actual type 2 false alarm rates for S1 responses
    #
    # est_HR2_rS2  = estimated type 2 hit rates for S2 responses
    # obs_HR2_rS2  = actual type 2 hit rates for S2 responses
    # est_FAR2_rS2 = estimated type 2 false alarm rates for S2 responses
    # obs_FAR2_rS2 = actual type 2 false alarm rates for S2 responses
    #
    # If there are N ratings, then there will be N-1 type 2 hit rates and false
    # alarm rates.

# 2017/06/22 - created R version
# 2015/07/23 - fixed bug for output fit.meta_ca and fit.S1units.meta_c1.
#            - added comments to help section as well as a warning output
#              for nR_S1 or nR_S2 inputs containing zeros
# 2014/10/14 - updated discussion of "s" input in the help section above.
# 2010/09/07 - created


## parse inputs

# check inputs
if (length(nR_S1)!=length(nR_S2)){
  print('input arrays must have the same number of elements')
  return()
}

if (min(nR_S1) == 0 | min(nR_S2) == 0){
  print(' ')
  print('WARNING!!')
  print('---------')
  print("Your inputs")
  cat("nR_S1 = [", nR_S1,"]")
  cat(" ")
  cat("nR_S2 = [", nR_S2,"]")
  print("contain zeros! This may interfere with proper estimation of meta-d'.")
  print("See 'help fit_meta_d_MLE' for more information.")
}

nRatings <- length(nR_S1) / 2;
nCriteria <- 2*nRatings - 1;

## set up constraints for MLE estimation

# parameters
# meta-d' - 1
# t2c     - nCriteria-1

A <- matrix(data = 0, nrow = nCriteria - 2, ncol = nCriteria);
b <- numeric(nCriteria - 2);

# constrain type 2 criteria values,
# such that t2c(i) is always <= t2c(i+1)
# want t2c(i)   <= t2c(i+1)
# -->  t2c(i+1) >= c(i) + 1e-5 (i.e. very small deviation from equality)
# -->  t2c(i) - t2c(i+1) <= -1e-5
for (i in 2 : (nCriteria - 1)){
  A[i-1,i:(i+1)] <- c(1,-1);
  b[i-1] <- -1e-5;
}

# lower bounds on parameters
LB <- c(-10, rep(-20,(nCriteria-1)/2), rep(0,(nCriteria-1)/2));
  ## meta-d'; criteria lower than t1c; criteria higher than t1c, respectively

# upper bounds on parameters
UB <- c(10, rep(0,(nCriteria-1)/2), rep(20,(nCriteria-1)/2))
  ## meta-d'; criteria lower than t1c; criteria higher than t1c, respectively



## select constant criterion type

constant_criterion <- expression(meta_d1 * (t1c1 / d1)); # relative criterion


## set up initial guess at parameter values

ratingHR  <- numeric(nRatings*2 - 1);
ratingFAR <- numeric(nRatings*2 - 1);
for(c in 2:(nRatings*2)){
  ratingHR[c-1] <- sum(nR_S2[c:length(nR_S2)]) / sum(nR_S2);
  ratingFAR[c-1] <- sum(nR_S1[c:length(nR_S2)]) / sum(nR_S1);
}

t1_index <- nRatings;
t2_index <- setdiff(1:(2*nRatings-1), t1_index);

d1 <- (1/s) * fninv( ratingHR[t1_index] ) - fninv( ratingFAR[t1_index] );
meta_d1 <- d1;

c1 <- (-1/(1+s)) * ( fninv( ratingHR ) + fninv( ratingFAR ) );
t1c1 <- c1[t1_index];
t2c1 <- c1[t2_index];

guess <- c(meta_d1, t2c1 - eval(constant_criterion));

#print(guess)

## find the best fit for type 2 hits and FAs

### save nR_S1 nR_S2 nRatings d1 t1c1 s constant_criterion fncdf fninv to global list "metadSettings"
metadSettings <<- list("fncdf" = fncdf, "fninv" = fninv, "nR_S1" = nR_S1, "nR_S2" = nR_S2,
                       "nRatings" = nRatings, "d1" = d1, "t1c1" = t1c1, "s" = s, "LB" = LB, "UB" = UB,
                       "constant_criterion" = constant_criterion, "A" = A, "b" = b)

#op = optimset(@fmincon);
#op = optimset(op,'MaxFunEvals',100000);

resNlopt <- nloptr::cobyla(x0 = guess, fn = fit_meta_d_logL, lower = LB, upper = UB,
                           hin = constraintFunc, control = list("maxeval" = 100000))

meta_d1  <- resNlopt$par[1];
t2c1     <- resNlopt$par[2:length(resNlopt$par)] + eval(constant_criterion);
logL     <- -resNlopt$value;


## data is fit, now to package it...

## find observed t2FAR and t2HR

# I_nR and C_nR are rating trial counts for incorrect and correct trials
# element i corresponds to # (in)correct w/ rating i
I_nR_rS2 = nR_S1[(nRatings+1):length(nR_S1)];
I_nR_rS1 = nR_S2[nRatings:1];

C_nR_rS2 = nR_S2[(nRatings+1):length(nR_S1)];
C_nR_rS1 = nR_S1[nRatings:1];

obs_FAR2_rS2 <- numeric(nRatings - 1)
obs_HR2_rS2  <- numeric(nRatings - 1)
obs_FAR2_rS1 <- numeric(nRatings - 1)
obs_HR2_rS1  <- numeric(nRatings - 1)

for (i in 2:nRatings){
  obs_FAR2_rS2[i-1] = sum( I_nR_rS2[i:length(I_nR_rS2)] ) / sum(I_nR_rS2);
  obs_HR2_rS2[i-1]  = sum( C_nR_rS2[i:length(C_nR_rS2)] ) / sum(C_nR_rS2);

  obs_FAR2_rS1[i-1] = sum( I_nR_rS1[i:length(I_nR_rS1)] ) / sum(I_nR_rS1);
  obs_HR2_rS1[i-1]  = sum( C_nR_rS1[i:length(C_nR_rS1)] ) / sum(C_nR_rS1);
}


## find estimated t2FAR and t2HR

S1mu <- -meta_d1/2;
S1sd <- 1;
S2mu <-  meta_d1/2;
S2sd <- S1sd/s;

mt1c1 <- eval(constant_criterion);

C_area_rS2 <- 1-fncdf(mt1c1,S2mu,S2sd);
I_area_rS2 <- 1-fncdf(mt1c1,S1mu,S1sd);

C_area_rS1 <- fncdf(mt1c1,S1mu,S1sd);
I_area_rS1 <- fncdf(mt1c1,S2mu,S2sd);

est_FAR2_rS2 <- numeric(nRatings-1)
est_HR2_rS2 <- numeric(nRatings-1)
est_FAR2_rS1 <- numeric(nRatings-1)
est_HR2_rS1 <- numeric(nRatings-1)

for (i in 1:(nRatings-1)){

  t2c1_lower <- t2c1[nRatings-i];
  t2c1_upper <- t2c1[nRatings-1+i];

  I_FAR_area_rS2 <- 1-fncdf(t2c1_upper,S1mu,S1sd);
  C_HR_area_rS2  <- 1-fncdf(t2c1_upper,S2mu,S2sd);

  I_FAR_area_rS1 <- fncdf(t2c1_lower,S2mu,S2sd);
  C_HR_area_rS1  <- fncdf(t2c1_lower,S1mu,S1sd);

  est_FAR2_rS2[i] <- I_FAR_area_rS2 / I_area_rS2;
  est_HR2_rS2[i]  <- C_HR_area_rS2 / C_area_rS2;

  est_FAR2_rS1[i] <- I_FAR_area_rS1 / I_area_rS1;
  est_HR2_rS1[i]  <- C_HR_area_rS1 / C_area_rS1;

}


## package output
mt1c1 <- eval(constant_criterion);
t2ca  <- ( sqrt(2)*s / sqrt(1+s^2) ) * t2c1;
meta_ca <- ( sqrt(2)*s / sqrt(1+s^2) ) * mt1c1
fit_meta_da <- sqrt(2/(1+s^2)) * s * meta_d1;
fit_da <- sqrt(2/(1+s^2)) * s * d1;


fitList <- list(
  "da"            = fit_da,
  "s"             = s,
  "nR_S1"         = nR_S1,
  "nR_S2"         = nR_S2,
  "meta_da"       = fit_meta_da,
  "M_diff"        = fit_meta_da - fit_da,
  "M_ratio"       = fit_meta_da / fit_da,
  "mt1c1"         = eval(constant_criterion),
  "meta_ca"       = meta_ca,
  "t2ca_rS1"      = t2ca[1:(nRatings-1)],
  "t2ca_rS2"      = t2ca[nRatings:length(t2ca)],
  "S1units"       = list("d1" = d1, "meta_d1" = meta_d1, "s" = s, "meta_c1" = mt1c1,
                         "t2c1_rS1" = t2c1[1:(nRatings-1)], "t2c1_rS2" = t2c1[nRatings:length(t2c1)]),
  "logL"          = logL,
  "est_HR2_rS1"   = est_HR2_rS1,
  "obs_HR2_rS1"   = obs_HR2_rS1,
  "est_FAR2_rS1"  = est_FAR2_rS1,
  "obs_FAR2_rS1"  = obs_FAR2_rS1,
  "est_HR2_rS2"   = est_HR2_rS2,
  "obs_HR2_rS2"   = obs_HR2_rS2,
  "est_FAR2_rS2"  = est_FAR2_rS2,
  "obs_FAR2_rS2"  = obs_FAR2_rS2
)

## clean up
rm(metadSettings)

return(fitList)
}


#' @title Calculate response-specific meta-d' using Maximum likelihood estimation.
#'
#' @description Given data from an experiment where an observer discriminates between two stimulus
#' alternatives on every trial and provides confidence ratings, provides a RESPONSE-SPECIFIC type 2
#' SDT analysis of the data. i.e. meta-d' is computed separately for each response type.
#' @param nR_S1 A vector. Contains the total number of responses in the first response category, conditional on presentation of S1 and S2.
#'
#' If nR_S1 = c(100,50,20,10,5,1), then when stimulus S1 was presented, the subject had the following response counts:
#'
#' responded S1, rating=3 : 100 times
#'
#' responded S1, rating=2 : 50 times
#'
#' responded S1, rating=1 : 20 times
#'
#' responded S2, rating=1 : 10 times
#'
#' responded S2, rating=2 : 5 times
#'
#' responded S2, rating=3 : 1 time
#' @param nR_S2 A vector. Contains the total number of responses in the second response category, conditional on presentation of S1 and S2.
#'
#' The ordering of response / rating counts for S2 should be the same as it is for S1.
#'
#' If nR_S2 = c(3,7,8,12,27,89), then when stimulus S2 was presented, the subject had the following response counts:
#'
#' responded S1, rating=3 : 3 times
#'
#' responded S1, rating=2 : 7 times
#'
#' responded S1, rating=1 : 8 times
#'
#' responded S2, rating=1 : 12 times
#'
#' responded S2, rating=2 : 27 times
#'
#' responded S2, rating=3 : 89 times
#' @param s A number. The ratio of standard deviations for type 1 distributions, i.e. s = sd(S1) / sd(S2). If not specified, s is set to a default value of 1.
#' @param fncdf A function for the CDF of the type 1 distribution. If not specified, fncdf defaults to \code{pnorm} (i.e. CDF for normal distribution).
#' @param fninv A function for the inverse CDF of the type 1 distribution. If not specified, fninv defaults to \code{qnorm} (i.e. inverse CDF for normal distribution).
#' @return A list. In the following, let S1 and S2 represent the distributions of evidence
#' generated by stimulus classes S1 and S2. Then the fields is a list are as follows:
#' da        = mean(S2) - mean(S1), in room-mean-square(sd(S1),sd(S2)) units
#'
#' t1ca      = type 1 criterion for overall data, RMS units
#'
#' s         = sd(S1) / sd(S2)
#'
#' meta_da_rS1   = meta-d' for S1 responses, RMS units
#'
#' t1ca_rS1      = type 1 criteron for meta-d' fit for S1 responses, RMS units
#'
#' t2ca_rS1      = type 2 criteria for meta-d' fit for S1 responses, RMS units
#'
#' M_ratio_rS1   = meta_da_rS1 / da
#'
#' M_diff_rS1    = meta_da_rS1 - da
#'
#' logL_rS1      = log likelihood of meta-d' fit for S1 responses
#'
#' obs_HR2_rS1   = actual type 2 hit rates for S1 responses
#'
#' est_HR2_rS1   = estimated type 2 hit rates for S1 responses
#'
#' obs_FAR2_rS1  = actual type 2 false alarm rates for S1 responses
#'
#' est_FAR2_rS1  = estimated type 2 false alarm rates for S1 responses
#'
#'
#' meta_da_rS2   = meta-d' for S2 responses, RMS units
#'
#' t1ca_rS2      = type 1 criteron for meta-d' fit for S2 responses, RMS units
#'
#' t2ca_rS2      = type 2 criteria for meta-d' fit for S2 responses, RMS units
#'
#' M_ratio_rS2   = meta_da_rS2 / da
#'
#' M_diff_rS2    = meta_da_rS2 - da
#'
#' logL_rS2      = log likelihood of meta-d' fit for S2 responses
#'
#' obs_HR2_rS2   = actual type 2 hit rates for S2 responses
#'
#' est_HR2_rS2   = estimated type 2 hit rates for S2 responses
#'
#' obs_FAR2_rS2  = actual type 2 false alarm rates for S2 responses
#'
#' est_FAR2_rS2  = estimated type 2 false alarm rates for S2 responses
#'
#' S1units contains same parameters in sd(S1) units.
#' @export

fit_rs_meta_d_MLE <- function(nR_S1, nR_S2, s = 1, fncdf = pnorm, fninv = qnorm){
  # check inputs
  if (length(nR_S1)%%2 !=0){
   warning('input arrays must have an even number of elements')
   return();
  }
  if (length(nR_S1)!=length(nR_S2)){
    warning('input arrays must have the same number of elements')
  }

  if (min(nR_S1) == 0 | min(nR_S2) == 0){
    print(' ')
    print('WARNING!!')
    print('---------')
    print("Your inputs")
    cat("nR_S1 = [", nR_S1,"]")
    cat(" ")
    cat("nR_S2 = [", nR_S2,"]")
    print("contain zeros! This may interfere with proper estimation of meta-d'.")
    print("See 'help fit_meta_d_MLE' for more information.")
  }

  nRatings <- length(nR_S1) / 2;
  nCriteria <- 2*nRatings - 1;


  ## find actual type 2 FAR and HR (data to be fit)

  # I_nR and C_nR are rating trial counts for incorrect and correct trials
  # element i corresponds to # (in)correct w/ rating i
  I_nR_rS2 <- nR_S1[(nRatings+1):length(nR_S1)];
  I_nR_rS1 <- nR_S2[nRatings:1];

  C_nR_rS2 <- nR_S2[(nRatings+1):length(nR_S2)];
  C_nR_rS1 <- nR_S1[nRatings:1];

  t2FAR_rS2 <- numeric(nRatings - 1)
  t2HR_rS2 <- numeric(nRatings - 1)
  t2FAR_rS1 <- numeric(nRatings - 1)
  t2HR_rS1 <- numeric(nRatings - 1)

  for (i in 2:nRatings){
    t2FAR_rS2[i-1] <- sum( I_nR_rS2[i:length(I_nR_rS2)] ) / sum(I_nR_rS2);
    t2HR_rS2[i-1]  <- sum( C_nR_rS2[i:length(C_nR_rS2)] ) / sum(C_nR_rS2);

    t2FAR_rS1[i-1] <- sum( I_nR_rS1[i:length(I_nR_rS1)] ) / sum(I_nR_rS1);
    t2HR_rS1[i-1]  <- sum( C_nR_rS1[i:length(C_nR_rS1)] ) / sum(C_nR_rS1);
  }

  t2FAR_rS1 <- t2FAR_rS1[length(t2FAR_rS1):1];
  t2HR_rS1  <- t2HR_rS1[length(t2HR_rS1):1];


  ## set up constraints

  # parameters
  # meta-d' - 1
  # t2c     - (nCriteria-1)/2

  nCritPerResp <- (nCriteria-1)/2;
  offset <- 1;

  A = matrix(0, nrow = nCritPerResp - 1, ncol = offset + nCritPerResp);
  b = numeric(nCritPerResp - 1);

  for (crit in 1:(nCritPerResp - 1)){
    A[crit,(offset+crit):(offset+crit+1)] <- c(1,-1);
    b[crit] <- -.001;
  }

    LB_rS1 <- c(-10, rep(-20,nCritPerResp));

    UB_rS1 <- c(10, rep(20,nCritPerResp));

    LB_rS2 <- c(-10, rep(-20,nCritPerResp));

    UB_rS2 <- c(10, rep(20,nCritPerResp));


  ## select constant criterion type

  constant_criterion_rS1 <- expression(meta_d1_rS1 * (t1c1 / d1)); # relative criterion
  constant_criterion_rS2 <- expression(meta_d1_rS2 * (t1c1 / d1)); # relative criterion


  ## set up guess

  ratingHR  <- numeric(nRatings*2 - 1);
  ratingFAR <- numeric(nRatings*2 - 1);
  for(c in 2:(nRatings*2)){
    ratingHR[c-1] <- sum(nR_S2[c:length(nR_S2)]) / sum(nR_S2);
    ratingFAR[c-1] <- sum(nR_S1[c:length(nR_S2)]) / sum(nR_S1);
  }

  t1_index <- nRatings;
  t2_index <- setdiff(1:(2*nRatings-1), t1_index);

  d1 <- (1/s) * fninv( ratingHR[t1_index] ) - fninv( ratingFAR[t1_index] );
  meta_d1_rS1 <- d1;
  meta_d1_rS2 <- d1;

  c1 <- (-1/(1+s)) * ( fninv( ratingHR ) + fninv( ratingFAR ) );
  t1c1 <- c1[t1_index];
  t2c1 <- c1[t2_index];

  guess_rS1 <- c(meta_d1_rS1, t2c1[1:nCritPerResp] - eval(constant_criterion_rS1));
  guess_rS2 <- c(meta_d1_rS2, t2c1[(nCritPerResp+1):length(t2c1)] - eval(constant_criterion_rS1));

  # save fitM_data.mat nR_S1 nR_S2 t2FAR_rS2 t2HR_rS2 t2FAR_rS1 t2HR_rS1 nRatings
  #      t1c1 s d1 fncdf constant_criterion_rS1 constant_criterion_rS2

  metadSettings_rs <<- list("fncdf" = fncdf, "fninv" = fninv, "nR_S1" = nR_S1, "nR_S2" = nR_S2,
                         "nRatings" = nRatings, "d1" = d1, "t1c1" = t1c1, "s" = s,
                         "A" = A, "b" = b, "t2FAR_rS2" = t2FAR_rS2, "t2HR_rS2" = t2HR_rS2,
                         "t2FAR_rS1" = t2FAR_rS1, "t2HR_rS1" = t2HR_rS1,
                         "constant_criterion_rS1" = constant_criterion_rS1,
                         "constant_criterion_rS2" = constant_criterion_rS2)


  ## fit for S1 responses

  ## find the best fit for type 2 hits and FAs
  resNlopt_rS1 <- nloptr::cobyla(x0 = guess_rS1, fn = fitM_rS1_logL, lower = LB_rS1, upper = UB_rS1,
                                 hin = constraintFunc_rs, control = list("maxeval" = 100000))

  meta_d1_rS1 <- resNlopt_rS1$par[1];
  meta_c1_rS1 <- eval(constant_criterion_rS1);

  t2c1_rS1    <- resNlopt_rS1$par[2:length(resNlopt_rS1$par)] + eval(constant_criterion_rS1);

  logL_rS1    <- -resNlopt_rS1$value;


  ## find model-estimated t2FAR and t2HR for S1 responses

  # find the estimated type 2 FAR and HR
  S1mu <- -meta_d1_rS1/2;
  S1sd <- 1;
  S2mu <-  meta_d1_rS1/2;
  S2sd <- S1sd/s;

  # adjust so that everything is centered on t1c1=0
  h <- 1-pnorm(0,S2mu,S2sd);
  f <- 1-pnorm(0,S1mu,S1sd);

  shift_c1 <- (-1/(1+s)) * (qnorm(h)+qnorm(f)); # this is the value of c1 midway b/t S1 and S2

  S1mu <- S1mu + shift_c1; # shift S1 and S2mu so that they lie on an axis for 0 --> c1=0
  S2mu <- S2mu + shift_c1;


  C_area_rS1 <- fncdf(meta_c1_rS1,S1mu,S1sd);
  I_area_rS1 <- fncdf(meta_c1_rS1,S2mu,S2sd);

  est_t2FAR_rS1 <- numeric(length(t2c1_rS1))
  est_t2HR_rS1  <- numeric(length(t2c1_rS1))

  for (i in 1:length(t2c1_rS1)){

    t2c1_lower <- t2c1_rS1[i];

    I_FAR_area_rS1 <- fncdf(t2c1_lower,S2mu,S2sd);
    C_HR_area_rS1  <- fncdf(t2c1_lower,S1mu,S1sd);

    est_t2FAR_rS1[i] <- I_FAR_area_rS1 / I_area_rS1;
    est_t2HR_rS1[i]  <- C_HR_area_rS1 / C_area_rS1;

  }



  ## fit for S2 responses

  ## find the best fit for type 2 hits and FAs

  resNlopt_rS2 <- nloptr::cobyla(x0 = guess_rS2, fn = fitM_rS2_logL, lower = LB_rS2, upper = UB_rS2,
                                 hin = constraintFunc_rs, control = list("maxeval" = 100000))

  meta_d1_rS2 <- resNlopt_rS2$par[1];
  meta_c1_rS2 <- eval(constant_criterion_rS2);

  t2c1_rS2    <- resNlopt_rS2$par[2:length(resNlopt_rS2$par)] + eval(constant_criterion_rS2);

  logL_rS2    <- -resNlopt_rS2$value;

  ## find estimated t2FAR and t2HR

  # find the estimated type 2 FAR and HR
  S1mu <- -meta_d1_rS2/2;
  S1sd <- 1;
  S2mu <-  meta_d1_rS2/2;
  S2sd <- S1sd/s;

  # adjust so that everything is centered on t1c1=0
  h <- 1-pnorm(0,S2mu,S2sd);
  f <- 1-pnorm(0,S1mu,S1sd);

  shift_c1 <- (-1/(1+s)) * (qnorm(h)+qnorm(f)); # this is the value of c1 midway b/t S1 and S2

  S1mu <- S1mu + shift_c1; # shift S1 and S2mu so that they lie on an axis for 0 --> c1=0
  S2mu <- S2mu + shift_c1;

  C_area_rS2 <- 1-fncdf(meta_c1_rS2,S2mu,S2sd);
  I_area_rS2 <- 1-fncdf(meta_c1_rS2,S1mu,S1sd);

  est_t2FAR_rS2 <- numeric(length(t2c1_rS2))
  est_t2HR_rS2  <- numeric(length(t2c1_rS2))

  for (i in 1:length(t2c1_rS2)){

    t2c1_upper <- t2c1_rS2[i];

    I_FAR_area_rS2 <- 1-fncdf(t2c1_upper,S1mu,S1sd);
    C_HR_area_rS2  <- 1-fncdf(t2c1_upper,S2mu,S2sd);

    est_t2FAR_rS2[i] <- I_FAR_area_rS2 / I_area_rS2;
    est_t2HR_rS2[i]  <- C_HR_area_rS2 / C_area_rS2;

  }


  ## package output

  # type 1 params
  meta_da_rS1 <- SDT_s_convert(meta_d1_rS1, s,'d1', expression(da))
  meta_da_rS2 <- SDT_s_convert(meta_d1_rS2, s,'d1', expression(da))
  da <- SDT_s_convert(d1,   s,'d1', expression(da))

  fitList <- list(da         = da,
                  t1ca       = SDT_s_convert(t1c1, s,'c1', expression(ca)),
                  s          = s,


                  # type 2 fits for rS1
                  meta_da_rS1   = meta_da_rS1,
                  t1ca_rS1      = SDT_s_convert(meta_c1_rS1, s,'c1', expression(ca)),
                  t2ca_rS1      = SDT_s_convert(t2c1_rS1,    s,'c1', expression(ca)),

                  M_ratio_rS1   = meta_da_rS1 / da,
                  M_diff_rS1    = meta_da_rS1 - da,

                  logL_rS1      = logL_rS1,

                  obs_HR2_rS1   = t2HR_rS1,
                  est_HR2_rS1   = est_t2HR_rS1,
                  obs_FAR2_rS1  = t2FAR_rS1,
                  est_FAR2_rS1  = est_t2FAR_rS1,


                  # type 2 fits for rS2
                  meta_da_rS2   = meta_da_rS2,
                  t1ca_rS2      = SDT_s_convert(meta_c1_rS2, s, 'c1', expression(ca)),
                  t2ca_rS2      = SDT_s_convert(t2c1_rS2,    s, 'c1', expression(ca)),

                  M_ratio_rS2   = meta_da_rS2 / da,
                  M_diff_rS2    = meta_da_rS2 - da,

                  logL_rS2      = logL_rS2,

                  obs_HR2_rS2   = t2HR_rS2,
                  est_HR2_rS2   = est_t2HR_rS2,
                  obs_FAR2_rS2  = t2FAR_rS2,
                  est_FAR2_rS2  = est_t2FAR_rS2,


                  # S1 units
                  S1units.d1   = d1,
                  S1units.t1c1 = t1c1,

                  S1units.meta_d1_rS1 = meta_d1_rS1,
                  S1units.t1c1_rS1    = meta_c1_rS1,
                  S1units.t2c1_rS1    = t2c1_rS1,

                  S1units.meta_d1_rS2 = meta_d1_rS2,
                  S1units.t1c1_rS2    = meta_c1_rS2,
                  S1units.t2c1_rS2    = t2c1_rS2 )


  return (fitList)
}



## function to find the likelihood of parameter values, given observed data
fit_meta_d_logL <- function(parameters){

  # set up parameters
  meta_d1  <- parameters[1]
  t2c1     <- parameters[2:length(parameters)]

  # loads:
  # nR_S1 nR_S2 nRatings d1 t1c1 s constant_criterion fncdf fninv
  nR_S1 <- metadSettings$nR_S1
  nR_S2 <- metadSettings$nR_S2
  t1c1 <- metadSettings$t1c1
  d1 <- metadSettings$d1
  s <- metadSettings$s
  nRatings <- metadSettings$nRatings
  fncdf <- metadSettings$fncdf
  fninv <- metadSettings$fninv

  # define mean and SD of S1 and S2 distributions
  S1mu <- -meta_d1/2
  S1sd <- 1
  S2mu <-  meta_d1/2
  S2sd <- S1sd/metadSettings$s;


  # adjust so that the type 1 criterion is set at 0
  # (this is just to work with optimization toolbox constraints...
  #  to simplify defining the upper and lower bounds of type 2 criteria)
  S1mu <- S1mu - eval(metadSettings$constant_criterion);
  S2mu <- S2mu - eval(metadSettings$constant_criterion);

  t1c1 <- 0;



  ### set up MLE analysis

  # get type 2 response counts
  nC_rS1 <- numeric(nRatings)
  nI_rS1 <- numeric(nRatings)
  nC_rS2 <- numeric(nRatings)
  nI_rS2 <- numeric(nRatings)

  for (i in 1:nRatings){

    # S1 responses
    nC_rS1[i] <- nR_S1[i];
    nI_rS1[i] <- nR_S2[i];

    # S2 responses
    nC_rS2[i] <- nR_S2[nRatings+i];
    nI_rS2[i] <- nR_S1[nRatings+i];

  }

  # get type 2 probabilities
  C_area_rS1 <- fncdf(t1c1,S1mu,S1sd);
  I_area_rS1 <- fncdf(t1c1,S2mu,S2sd);

  C_area_rS2 <- 1-fncdf(t1c1,S2mu,S2sd);
  I_area_rS2 <- 1-fncdf(t1c1,S1mu,S1sd);

  t2c1x <- c(-Inf, t2c1[1:(nRatings-1)], t1c1, t2c1[nRatings:length(t2c1)], Inf);

  prC_rS1 <- numeric(nRatings)
  prI_rS1 <- numeric(nRatings)
  prC_rS2 <- numeric(nRatings)
  prI_rS2 <- numeric(nRatings)

  for (i in 1:nRatings){
    prC_rS1[i] = ( fncdf(t2c1x[i+1],S1mu,S1sd) - fncdf(t2c1x[i],S1mu,S1sd) ) / C_area_rS1;
    prI_rS1[i] = ( fncdf(t2c1x[i+1],S2mu,S2sd) - fncdf(t2c1x[i],S2mu,S2sd) ) / I_area_rS1;

    prC_rS2[i] = ( (1-fncdf(t2c1x[nRatings+i],S2mu,S2sd)) - (1-fncdf(t2c1x[nRatings+i+1],S2mu,S2sd)) ) / C_area_rS2;
    prI_rS2[i] = ( (1-fncdf(t2c1x[nRatings+i],S1mu,S1sd)) - (1-fncdf(t2c1x[nRatings+i+1],S1mu,S1sd)) ) / I_area_rS2;
  }


  # calculate logL
  logL = 0;
  for (i in 1:nRatings){
    logL = logL + nC_rS1[i]*log(prC_rS1[i]) + nI_rS1[i]*log(prI_rS1[i]) +
      nC_rS2[i]*log(prC_rS2[i]) + nI_rS2[i]*log(prI_rS2[i]);
  }

  if (is.na(logL)){
    logL=-Inf
  }

  logL = -logL;

  return(logL)
}

## function to find the SSE for type 2 HR and FAR for S1 responses
fitM_rS1_logL <- function(input){

  # set up parameters
  meta_d1_rS1 <- input[1];
  t2c1        <- input[2:length(input)];

  # load the behavioral data
  # t2FAR_rS2 t2HR_rS2 t2FAR_rS1 t2HR_rS1 nRatings t1c1 s d1 fncdf constant_criterion_rS1 constant_criterion_rS2
  nR_S1     <- metadSettings_rs$nR_S1
  nR_S2     <- metadSettings_rs$nR_S2
  t2FAR_rS2 <- metadSettings_rs$t2FAR_rS2
  t2HR_rS2  <- metadSettings_rs$t2HR_rS2
  t2FAR_rS1 <- metadSettings_rs$t2FAR_rS1
  t2HR_rS1  <- metadSettings_rs$t2HR_rS1
  nRatings  <- metadSettings_rs$nRatings
  t1c1      <- metadSettings_rs$t1c1
  s         <- metadSettings_rs$s
  d1        <- metadSettings_rs$d1
  fncdf     <- metadSettings_rs$fncdf
  constant_criterion_rS1 <- metadSettings_rs$constant_criterion_rS1
  constant_criterion_rS2 <- metadSettings_rs$constant_criterion_rS2

  # find the estimated type 2 FAR and HR
  S1mu <- -meta_d1_rS1/2;
  S1sd <- 1;
  S2mu <-  meta_d1_rS1/2;
  S2sd <- S1sd/s;


  # adjust so that everything is centered on t1c1=0
  h <- 1-pnorm(0,S2mu,S2sd);
  f <- 1-pnorm(0,S1mu,S1sd);

  shift_c1 <- (-1/(1+s)) * (qnorm(h)+qnorm(f)); # this is the value of c1 midway b/t S1 and S2

  S1mu <- S1mu + shift_c1; # shift S1 and S2mu so that they lie on an axis for 0 --> c1=0
  S2mu <- S2mu + shift_c1;


  # adjust so that t1c1 = 0
  S1mu <- S1mu - eval(constant_criterion_rS1);
  S2mu <- S2mu - eval(constant_criterion_rS1);
  t1c1 <- 0;



  ### set up MLE analysis

  # get type 2 response counts
  nC_rS1 <- numeric(nRatings)
  nI_rS1 <- numeric(nRatings)

  for (ic in 1:nRatings){

    # S1 responses
    nC_rS1[ic] = nR_S1[ic];
    nI_rS1[ic] = nR_S2[ic];

  }

  # get type 2 probabilities
  C_area_rS1 <- fncdf(t1c1,S1mu,S1sd);
  I_area_rS1 <- fncdf(t1c1,S2mu,S2sd);


  t2c1x <- c(-Inf, t2c1, t1c1);

  # for ic = 1:4
  prC_rS1 <- numeric(nRatings)
  prI_rS1 <- numeric(nRatings)

  for (ic in 1:nRatings){
    prC_rS1[ic] <- ( fncdf(t2c1x[ic+1],S1mu,S1sd) - fncdf(t2c1x[ic],S1mu,S1sd) ) / C_area_rS1;
    prI_rS1[ic] <- ( fncdf(t2c1x[ic+1],S2mu,S2sd) - fncdf(t2c1x[ic],S2mu,S2sd) ) / I_area_rS1;
  }

  logL <- 0;
  for (ic in 1:nRatings){
    logL <- logL + nC_rS1[ic]*log(prC_rS1[ic]) + nI_rS1[ic]*log(prI_rS1[ic]);
  }

  if (is.na(logL)){
    logL=-Inf
  }

  logL <- -logL;

  return(logL)
}

## function to find the SSE for type 2 HR and FAR for S2 responses
fitM_rS2_logL <- function(input){

  # set up parameters
  meta_d1_rS2 <- input[1];
  t2c1        <- input[2:length(input)];

  # load the behavioral data
  nR_S1     <- metadSettings_rs$nR_S1
  nR_S2     <- metadSettings_rs$nR_S2
  t2FAR_rS2 <- metadSettings_rs$t2FAR_rS2
  t2HR_rS2  <- metadSettings_rs$t2HR_rS2
  t2FAR_rS1 <- metadSettings_rs$t2FAR_rS1
  t2HR_rS1  <- metadSettings_rs$t2HR_rS1
  nRatings  <- metadSettings_rs$nRatings
  t1c1      <- metadSettings_rs$t1c1
  s         <- metadSettings_rs$s
  d1        <- metadSettings_rs$d1
  fncdf     <- metadSettings_rs$fncdf
  constant_criterion_rS1 <- metadSettings_rs$constant_criterion_rS1
  constant_criterion_rS2 <- metadSettings_rs$constant_criterion_rS2



  # find the estimated type 2 FAR and HR
  S1mu <- -meta_d1_rS2/2;
  S1sd <- 1;
  S2mu <-  meta_d1_rS2/2;
  S2sd <- S1sd/s;


  # adjust so that everything is centered on t1c1=0
  h <- 1-pnorm(0,S2mu,S2sd);
  f <- 1-pnorm(0,S1mu,S1sd);

  shift_c1 <- (-1/(1+s)) * (qnorm(h)+qnorm(f)); # this is the value of c1 midway b/t S1 and S2

  S1mu <- S1mu + shift_c1; # shift S1 and S2mu so that they lie on an axis for 0 --> c1=0
  S2mu <- S2mu + shift_c1;


  # adjust so that t1c1 = 0
  S1mu <- S1mu - eval(constant_criterion_rS2);
  S2mu <- S2mu - eval(constant_criterion_rS2);
  t1c1 <- 0;



  ### set up MLE analysis

  # get type 2 response counts
  nC_rS2 <- numeric(nRatings)
  nI_rS2 <- numeric(nRatings)

  for (ic in 1:nRatings){

    # S2 responses
    nC_rS2[ic] <- nR_S2[ic+nRatings];
    nI_rS2[ic] <- nR_S1[ic+nRatings];

  }

  # get type 2 probabilities
  C_area_rS2 <- 1-fncdf(t1c1,S2mu,S2sd);
  I_area_rS2 <- 1-fncdf(t1c1,S1mu,S1sd);

  t2c1x = c(t1c1, t2c1, Inf);

  prC_rS2 <- numeric(nRatings)
  prI_rS2 <- numeric(nRatings)

  for (ic in 1:nRatings){
    prC_rS2[ic] <- ( (1-fncdf(t2c1x[ic],S2mu,S2sd)) - (1-fncdf(t2c1x[ic+1],S2mu,S2sd)) ) / C_area_rS2;
    prI_rS2[ic] <- ( (1-fncdf(t2c1x[ic],S1mu,S1sd)) - (1-fncdf(t2c1x[ic+1],S1mu,S1sd)) ) / I_area_rS2;
  }


  logL <- 0;
  for (ic in 1:nRatings){
    logL <- logL + nC_rS2[ic]*log(prC_rS2[ic]) + nI_rS2[ic]*log(prI_rS2[ic]);
  }

  if (is.na(logL)){
    logL <- -Inf
  }

  logL <- -logL;

  return(logL)
}

SDT_s_convert <- function(in1,s,cont1,cont2){

  # out = SDT_s_convert(in1,s,cont1,cont2)
  #
  # in1 - input var to be converted
  # s - sd(S1)/sd(S2)
  # cont1 - identity of in1
  # cont2 - identity of out
  #
  # cont1 and cont2 can be these tokens
  # 'da','d1','d2','ca','c1',c2'

  # convert d'
  # s = d2 / d1
  # da = sqrt(2./(1+s.^2)) .* d2
  if (cont1 == "da"){
    da = in1;
    d2 = da / sqrt(2/(1+s^2));
    d1 = d2 / s;
  }
  if (cont1 == "d1"){
    d1 = in1;
    d2 = d1 * s;
    da = sqrt(2/(1+s^2)) * d2;
  }

  if (cont1 == "d2"){
    d2 = in1;
    d1 = d2 / s;
    da = sqrt(2/(1+s^2)) * d2;
  }

  # convert c
  # s = c2 / c1
  # ca = ( sqrt(2).*s ./ sqrt(1+s.^2) ) .* c1;
  if (cont1 == "ca"){
    ca = in1;
    c1 = ( sqrt(1+s^2) / sqrt(2)*s ) * ca;
    c2 = c1 * s;
  }

  if (cont1 == "c1"){
    c1 = in1;
    c2 = c1 * s;
    ca = ( sqrt(2)*s / sqrt(1+s^2) ) * c1;
  }

  if (cont1 == "c2"){
    c2 = in1;
    c1 = c2 / s;
    ca = ( sqrt(2)*s / sqrt(1+s^2) ) * c1;
  }

  return (eval(cont2));

}

## function of constrains
constraintFunc <- function(parameters){
  A <- metadSettings$A
  b <- metadSettings$b
  return(as.numeric(t(t(b))-A%*%parameters))
}

constraintFunc_rs <- function(parameters){
  A <- metadSettings_rs$A
  b <- metadSettings_rs$b
  return(as.numeric(t(t(b))-A%*%parameters))
}
