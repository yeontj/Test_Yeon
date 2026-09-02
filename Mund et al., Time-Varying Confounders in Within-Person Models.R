#################################################################
## THIS FILE IS LICENSED UNDER THE TERMS AND CONDITIONS OF THE ##
## CC-BY 4.0 INTERNATIONAL                                     ##
## for details on the license, see                             ##
## http://creativecommons.org/licenses/by/4.0/                 ##
#################################################################

# (C) 2021 Marcus Mund, Matthew D. Johnson, Steffen Nestler

#### +++++++++++++++++++++++++++++++++++++++++++++++++
####   PACKAGES                                   ####
#### +++++++++++++++++++++++++++++++++++++++++++++++++

library(rio)
library(data.table)
library(lavaan)
library(xtable)
library(ggplot2)
library(MBESS)
library(corx)

#### +++++++++++++++++++++++++++++++++++++++++++++++++
####   UNCONDITIONAL MODELS                       ####
#### +++++++++++++++++++++++++++++++++++++++++++++++++

### all model syntax adapted from Mund & Nestler, 2019, doi:10.1016/j.alcr.2018.10.002


#### Random-Intercept Cross-Lagged Panel Model ####

riclpm.mod <- '
# Define intercept factors
ix =~ 1*sat6.1  + 1*sat6.2  + 1*sat6.3  + 1*sat6.4
iy =~ 1*log_inc.1 + 1*log_inc.2 + 1*log_inc.3 + 1*log_inc.4

# Define phantom latent variables
etax1 =~ 1*sat6.1
etax2 =~ 1*sat6.2
etax3 =~ 1*sat6.3
etax4 =~ 1*sat6.4

etay1 =~ 1*log_inc.1
etay2 =~ 1*log_inc.2
etay3 =~ 1*log_inc.3
etay4 =~ 1*log_inc.4

# Autoregressive effects
etax2 ~ a1*etax1
etax3 ~ a1*etax2
etax4 ~ a1*etax3

etay2 ~ a2*etay1
etay3 ~ a2*etay2
etay4 ~ a2*etay3

# Crosslagged effects
etay2 ~ c1*etax1
etay3 ~ c1*etax2
etay4 ~ c1*etax3

etax2 ~ c2*etay1
etax3 ~ c2*etay2
etax4 ~ c2*etay3

# Some further constraints on the variance structure
# 1. Set error variances of the observed variables to zero
sat6.1 ~~ 0*sat6.1
sat6.2 ~~ 0*sat6.2
sat6.3 ~~ 0*sat6.3
sat6.4 ~~ 0*sat6.4
     
log_inc.1 ~~ 0*log_inc.1
log_inc.2 ~~ 0*log_inc.2
log_inc.3 ~~ 0*log_inc.3
log_inc.4 ~~ 0*log_inc.4

# 2. Let lavaan estimate the variance of the latent variables
etax1 ~~ varx1*etax1
etax2 ~~ varx2*etax2
etax3 ~~ varx3*etax3
etax4 ~~ varx4*etax4
     
etay1 ~~ vary1*etay1
etay2 ~~ vary2*etay2
etay3 ~~ vary3*etay3
etay4 ~~ vary4*etay4

# 3. We also want estimates of the intercept factor variances and an
#    estimate of their covariance
ix ~~ varix*ix
iy ~~ variy*iy
ix ~~ covi*iy

# 4. We have to define that the covariance between the intercepts and
#    the latents of the first time point are zero
etax1 ~~ 0*ix
etay1 ~~ 0*ix
etax1 ~~ 0*iy
etay1 ~~ 0*iy

# 5. Finally, we estimate the covariance between the latents of x and y
#    of the first time point, the second time-point and so on. note that
#    for the second to fourth time point the correlation is constrained to
#    the same value
etax1 ~~ cov1*etay1
etax2 ~~ cove*etay2
etax3 ~~ cove*etay3
etax4 ~~ cove*etay4

# The model also contains a mean structure and we have to define some
# constraints for this part of the model. the assumption is that we
# want estimates of the manifest variables. the means of the intercept
# factors will be constrained to 0
sat6.1 ~ 1
sat6.2 ~ 1
sat6.3 ~ 1
sat6.4 ~ 1

log_inc.1 ~ 1
log_inc.2 ~ 1
log_inc.3 ~ 1
log_inc.4 ~ 1
     
etax1 ~ 0*1
etax2 ~ 0*1
etax3 ~ 0*1
etax4 ~ 0*1
     
etay1 ~ 0*1
etay2 ~ 0*1
etay3 ~ 0*1
etay4 ~ 0*1
     
ix ~ 0*1
iy ~ 0*1

# define correlations for CI
cori  := covi / (sqrt(varix) * sqrt(variy))
cort1 := cov1 / (sqrt(varx1) * sqrt(vary1))
cort2 := cove / (sqrt(varx2) * sqrt(vary2))
cort3 := cove / (sqrt(varx3) * sqrt(vary3))
cort4 := cove / (sqrt(varx4) * sqrt(vary4))
'

riclpm.fit <- sem(riclpm.mod, missing = "fiml", data = pf, verbose = T)
capture.output(summary(riclpm.fit, fit.measures = T, std = T, ci = T), file = "./results/unconditional/RI-CLPM.txt")
# inspect(riclpm.fit, what = "cor.lv")




#### +++++++++++++++++++++++++++++++++++++++++++++++++
####   MODELS WITH TIC                            ####
#### +++++++++++++++++++++++++++++++++++++++++++++++++

#### Cross-Lagged Panel Model (unconstrained tic) ####

clpm_tic_con.mod <- '
# Define phantom latent variables
etax1 =~ 1*sat6.1
etax2 =~ 1*sat6.2
etax3 =~ 1*sat6.3
etax4 =~ 1*sat6.4

etay1 =~ 1*log_inc.1
etay2 =~ 1*log_inc.2
etay3 =~ 1*log_inc.3
etay4 =~ 1*log_inc.4

# Autoregressive effects
etax2 ~ a1*etax1
etax3 ~ a1*etax2
etax4 ~ a1*etax3

etay2 ~ a2*etay1
etay3 ~ a2*etay2
etay4 ~ a2*etay3

# Crosslagged effects
etay2 ~ c1*etax1
etay3 ~ c1*etax2
etay4 ~ c1*etax3

etax2 ~ c2*etay1
etax3 ~ c2*etay2
etax4 ~ c2*etay3

# Some further constraints on the variance structure
# 1. Set error variances of the observed variables to the same value

sat6.1 ~~ sat6.1
sat6.2 ~~ sat6.2
sat6.3 ~~ sat6.3
sat6.4 ~~ sat6.4

log_inc.1 ~~ log_inc.1
log_inc.2 ~~ log_inc.2
log_inc.3 ~~ log_inc.3
log_inc.4 ~~ log_inc.4

# 2. Set the variance of the latent variables, except the first one, to the same value
etax2 ~~ varx*etax2
etax3 ~~ varx*etax3
etax4 ~~ varx*etax4

etay2 ~~ vary*etay2
etay3 ~~ vary*etay3
etay4 ~~ vary*etay4

# 3. Define the variance and the covariance for the first latent variable

etax1 ~~ varx1*etax1
etay1 ~~ vary1*etay1

# 4. Finally, we estimate the covariance between the latents of x and y
#    of the first time point, the second time-point and so on. note that
#    for the second to fourth time point the correlation is constrained to
#    the same value
etax1 ~~ cov1*etay1
etax2 ~~ cove*etay2
etax3 ~~ cove*etay3
etax4 ~~ cove*etay4

# 5. The model also contains a mean structure and we have to define some
#    constraints for this part of the model. We set the intercepts of the
#    observed indicators to zero and estimate the means of the latent
#    variables:

sat6.1 ~ 0*1
sat6.2 ~ 0*1
sat6.3 ~ 0*1
sat6.4 ~ 0*1

log_inc.1 ~ 0*1
log_inc.2 ~ 0*1
log_inc.3 ~ 0*1
log_inc.4 ~ 0*1

etax1 ~ 1
etax2 ~ 1
etax3 ~ 1
etax4 ~ 1
etay1 ~ 1
etay2 ~ 1
etay3 ~ 1
etay4 ~ 1

# incorporating tic
etax1 + etax2 + etax3 + etax4 ~ conx*sex_gen
etay1 + etay2 + etay3 + etay4 ~ cony*sex_gen

sex_gen ~~ sex_gen

# define correlations for CI
cort1 := cov1 / (sqrt(varx1) * sqrt(vary1))
cort2 := cove / (sqrt(varx) * sqrt(vary))
cort3 := cove / (sqrt(varx) * sqrt(vary))
cort4 := cove / (sqrt(varx) * sqrt(vary))
' 

clpm_tic_con.fit <- sem(clpm_tic_con.mod, missing = "fiml", data = pf, verbose = T)
capture.output(summary(clpm_tic_con.fit, fit.measures = T, std = T, ci = T), file = "./results/TIC/CLPM_constrained.txt")

anova(clpm_tic.fit, clpm_tic_con.fit)


#### Random-Intercept Cross-Lagged Panel Model (tic as constant effect) ####

riclpm_tic.mod <- '
# Define intercept factors
ix =~ 1*sat6.1  + 1*sat6.2  + 1*sat6.3  + 1*sat6.4
iy =~ 1*log_inc.1 + 1*log_inc.2 + 1*log_inc.3 + 1*log_inc.4

# Define phantom latent variables
etax1 =~ 1*sat6.1
etax2 =~ 1*sat6.2
etax3 =~ 1*sat6.3
etax4 =~ 1*sat6.4

etay1 =~ 1*log_inc.1
etay2 =~ 1*log_inc.2
etay3 =~ 1*log_inc.3
etay4 =~ 1*log_inc.4

# Autoregressive effects
etax2 ~ a1*etax1
etax3 ~ a1*etax2
etax4 ~ a1*etax3

etay2 ~ a2*etay1
etay3 ~ a2*etay2
etay4 ~ a2*etay3

# Crosslagged effects
etay2 ~ c1*etax1
etay3 ~ c1*etax2
etay4 ~ c1*etax3

etax2 ~ c2*etay1
etax3 ~ c2*etay2
etax4 ~ c2*etay3

# Some further constraints on the variance structure
# 1. Set error variances of the observed variables to zero
sat6.1 ~~ 0*sat6.1
sat6.2 ~~ 0*sat6.2
sat6.3 ~~ 0*sat6.3
sat6.4 ~~ 0*sat6.4

log_inc.1 ~~ 0*log_inc.1
log_inc.2 ~~ 0*log_inc.2
log_inc.3 ~~ 0*log_inc.3
log_inc.4 ~~ 0*log_inc.4

# 2. Let lavaan estimate the variance of the latent variables
etax1 ~~ varx1*etax1
etax2 ~~ varx2*etax2
etax3 ~~ varx3*etax3
etax4 ~~ varx4*etax4

etay1 ~~ vary1*etay1
etay2 ~~ vary2*etay2
etay3 ~~ vary3*etay3
etay4 ~~ vary4*etay4

# 3. We also want estimates of the intercept factor variances and an
#    estimate of their covariance
ix ~~ varix*ix
iy ~~ variy*iy
ix ~~ covi*iy

# 4. We have to define that the covariance between the intercepts and
#    the latents of the first time point are zero
etax1 ~~ 0*ix
etay1 ~~ 0*ix
etax1 ~~ 0*iy
etay1 ~~ 0*iy

# 5. Finally, we estimate the covariance between the latents of x and y
#    of the first time point, the second time-point and so on. note that
#    for the second to fourth time point the correlation is constrained to
#    the same value
etax1 ~~ cov1*etay1
etax2 ~~ cove*etay2
etax3 ~~ cove*etay3
etax4 ~~ cove*etay4

# The model also contains a mean structure and we have to define some
# constraints for this part of the model. the assumption is that we
# want estimates of the manifest variables. the means of the intercept
# factors will be constrained to 0
sat6.1 ~ 1
sat6.2 ~ 1
sat6.3 ~ 1
sat6.4 ~ 1

log_inc.1 ~ 1
log_inc.2 ~ 1
log_inc.3 ~ 1
log_inc.4 ~ 1

etax1 ~ 0*1
etax2 ~ 0*1
etax3 ~ 0*1
etax4 ~ 0*1

etay1 ~ 0*1
etay2 ~ 0*1
etay3 ~ 0*1
etay4 ~ 0*1

ix ~ 0*1
iy ~ 0*1


# incorporating tic
sat6.1    + sat6.2    + sat6.3    + sat6.4    ~ ticx*sex_gen
log_inc.1 + log_inc.2 + log_inc.3 + log_inc.4 ~ ticy*sex_gen

sex_gen ~~ sex_gen

# define correlations for CI
cori  := covi / (sqrt(varix) * sqrt(variy))
cort1 := cov1 / (sqrt(varx1) * sqrt(vary1))
cort2 := cove / (sqrt(varx2) * sqrt(vary2))
cort3 := cove / (sqrt(varx3) * sqrt(vary3))
cort4 := cove / (sqrt(varx4) * sqrt(vary4))
'

riclpm_tic.fit <- sem(riclpm_tic.mod, missing = "fiml", data = pf, verbose = T)
capture.output(summary(riclpm_tic.fit, fit.measures = T, std = T, ci = T), file = "./results/TIC/RI-CLPM.txt")


#### Random-Intercept Cross-Lagged Panel Model (tic with time-varying effect) ####

riclpm_tic_var.mod <- '
# Define intercept factors
ix =~ 1*sat6.1  + 1*sat6.2  + 1*sat6.3  + 1*sat6.4
iy =~ 1*log_inc.1 + 1*log_inc.2 + 1*log_inc.3 + 1*log_inc.4

# Define phantom latent variables
etax1 =~ 1*sat6.1
etax2 =~ 1*sat6.2
etax3 =~ 1*sat6.3
etax4 =~ 1*sat6.4

etay1 =~ 1*log_inc.1
etay2 =~ 1*log_inc.2
etay3 =~ 1*log_inc.3
etay4 =~ 1*log_inc.4

# Autoregressive effects
etax2 ~ a1*etax1
etax3 ~ a1*etax2
etax4 ~ a1*etax3

etay2 ~ a2*etay1
etay3 ~ a2*etay2
etay4 ~ a2*etay3

# Crosslagged effects
etay2 ~ c1*etax1
etay3 ~ c1*etax2
etay4 ~ c1*etax3

etax2 ~ c2*etay1
etax3 ~ c2*etay2
etax4 ~ c2*etay3

# Some further constraints on the variance structure
# 1. Set error variances of the observed variables to zero
sat6.1 ~~ 0*sat6.1
sat6.2 ~~ 0*sat6.2
sat6.3 ~~ 0*sat6.3
sat6.4 ~~ 0*sat6.4

log_inc.1 ~~ 0*log_inc.1
log_inc.2 ~~ 0*log_inc.2
log_inc.3 ~~ 0*log_inc.3
log_inc.4 ~~ 0*log_inc.4

# 2. Let lavaan estimate the variance of the latent variables
etax1 ~~ varx1*etax1
etax2 ~~ varx2*etax2
etax3 ~~ varx3*etax3
etax4 ~~ varx4*etax4

etay1 ~~ vary1*etay1
etay2 ~~ vary2*etay2
etay3 ~~ vary3*etay3
etay4 ~~ vary4*etay4

# 3. We also want estimates of the intercept factor variances and an
#    estimate of their covariance
ix ~~ varix*ix
iy ~~ variy*iy
ix ~~ covi*iy

# 4. We have to define that the covariance between the intercepts and
#    the latents of the first time point are zero
etax1 ~~ 0*ix
etay1 ~~ 0*ix
etax1 ~~ 0*iy
etay1 ~~ 0*iy

# 5. Finally, we estimate the covariance between the latents of x and y
#    of the first time point, the second time-point and so on. note that
#    for the second to fourth time point the correlation is constrained to
#    the same value
etax1 ~~ cov1*etay1
etax2 ~~ cove*etay2
etax3 ~~ cove*etay3
etax4 ~~ cove*etay4

# The model also contains a mean structure and we have to define some
# constraints for this part of the model. the assumption is that we
# want estimates of the manifest variables. the means of the intercept
# factors will be constrained to 0
sat6.1 ~ 1
sat6.2 ~ 1
sat6.3 ~ 1
sat6.4 ~ 1

log_inc.1 ~ 1
log_inc.2 ~ 1
log_inc.3 ~ 1
log_inc.4 ~ 1

etax1 ~ 0*1
etax2 ~ 0*1
etax3 ~ 0*1
etax4 ~ 0*1

etay1 ~ 0*1
etay2 ~ 0*1
etay3 ~ 0*1
etay4 ~ 0*1

ix ~ 0*1
iy ~ 0*1


# incorporating tic
sat6.1    + sat6.2    + sat6.3    + sat6.4    ~ sex_gen
log_inc.1 + log_inc.2 + log_inc.3 + log_inc.4 ~ sex_gen

sex_gen ~~ sex_gen

# define correlations for CI
cori  := covi / (sqrt(varix) * sqrt(variy))
cort1 := cov1 / (sqrt(varx1) * sqrt(vary1))
cort2 := cove / (sqrt(varx2) * sqrt(vary2))
cort3 := cove / (sqrt(varx3) * sqrt(vary3))
cort4 := cove / (sqrt(varx4) * sqrt(vary4))
'

riclpm_tic_var.fit <- sem(riclpm_tic_var.mod, missing = "fiml", data = pf, verbose = T)
capture.output(summary(riclpm_tic_var.fit, fit.measures = T, std = T, ci = T), file = "./results/TIC/RI-CLPM_varyingTIC.txt")



#### +++++++++++++++++++++++++++++++++++++++++++++++++
####   MODELS WITH TVC                            ####
#### +++++++++++++++++++++++++++++++++++++++++++++++++


#### Random-Intercept Cross-Lagged Panel Model ####

riclpm_tvc.mod <- '
# Define intercept factors
ix =~ 1*sat6.1  + 1*sat6.2  + 1*sat6.3  + 1*sat6.4
iy =~ 1*log_inc.1 + 1*log_inc.2 + 1*log_inc.3 + 1*log_inc.4

# Define phantom latent variables
etax1 =~ 1*sat6.1
etax2 =~ 1*sat6.2
etax3 =~ 1*sat6.3
etax4 =~ 1*sat6.4

etay1 =~ 1*log_inc.1
etay2 =~ 1*log_inc.2
etay3 =~ 1*log_inc.3
etay4 =~ 1*log_inc.4

# Autoregressive effects
etax2 ~ a1*etax1
etax3 ~ a1*etax2
etax4 ~ a1*etax3

etay2 ~ a2*etay1
etay3 ~ a2*etay2
etay4 ~ a2*etay3

# Crosslagged effects
etay2 ~ c1*etax1
etay3 ~ c1*etax2
etay4 ~ c1*etax3

etax2 ~ c2*etay1
etax3 ~ c2*etay2
etax4 ~ c2*etay3

# Some further constraints on the variance structure
# 1. Set error variances of the observed variables to zero
sat6.1 ~~ 0*sat6.1
sat6.2 ~~ 0*sat6.2
sat6.3 ~~ 0*sat6.3
sat6.4 ~~ 0*sat6.4

log_inc.1 ~~ 0*log_inc.1
log_inc.2 ~~ 0*log_inc.2
log_inc.3 ~~ 0*log_inc.3
log_inc.4 ~~ 0*log_inc.4

# 2. Let lavaan estimate the variance of the latent variables
etax1 ~~ varx1*etax1
etax2 ~~ varx2*etax2
etax3 ~~ varx3*etax3
etax4 ~~ varx4*etax4

etay1 ~~ vary1*etay1
etay2 ~~ vary2*etay2
etay3 ~~ vary3*etay3
etay4 ~~ vary4*etay4

# 3. We also want estimates of the intercept factor variances and an
#    estimate of their covariance
ix ~~ varix*ix
iy ~~ variy*iy
ix ~~ covi*iy

# 4. We have to define that the covariance between the intercepts and
#    the latents of the first time point are zero
etax1 ~~ 0*ix
etay1 ~~ 0*ix
etax1 ~~ 0*iy
etay1 ~~ 0*iy

# 5. Finally, we estimate the covariance between the latents of x and y
#    of the first time point, the second time-point and so on. note that
#    for the second to fourth time point the correlation is constrained to
#    the same value
etax1 ~~ cov1*etay1
etax2 ~~ cove*etay2
etax3 ~~ cove*etay3
etax4 ~~ cove*etay4

# The model also contains a mean structure and we have to define some
# constraints for this part of the model. the assumption is that we
# want estimates of the manifest variables. the means of the intercept
# factors will be constrained to 0
sat6.1 ~ 1
sat6.2 ~ 1
sat6.3 ~ 1
sat6.4 ~ 1

log_inc.1 ~ 1
log_inc.2 ~ 1
log_inc.3 ~ 1
log_inc.4 ~ 1

etax1 ~ 0*1
etax2 ~ 0*1
etax3 ~ 0*1
etax4 ~ 0*1

etay1 ~ 0*1
etay2 ~ 0*1
etay3 ~ 0*1
etay4 ~ 0*1

ix ~ 0*1
iy ~ 0*1

# incorporate tvc
etax1 + etay1 ~ se.1
etax2 + etay2 ~ se.2
etax3 + etay3 ~ se.3
etax4 + etay4 ~ se.4

se.1 ~~ se.1
se.2 ~~ se.2
se.3 ~~ se.3
se.4 ~~ se.4

se.1 ~~ se.2 + se.3 + se.4
se.2 ~~ se.3 + se.4
se.3 ~~ se.4


# define correlations for CI
cori  := covi / (sqrt(varix) * sqrt(variy))
cort1 := cov1 / (sqrt(varx1) * sqrt(vary1))
cort2 := cove / (sqrt(varx2) * sqrt(vary2))
cort3 := cove / (sqrt(varx3) * sqrt(vary3))
cort4 := cove / (sqrt(varx4) * sqrt(vary4))
'

riclpm_tvc.fit <- sem(riclpm_tvc.mod, missing = "fiml", data = pf, verbose = T)
capture.output(summary(riclpm_tvc.fit, fit.measures = T, std = T, ci = T), file = "./results/TVC/RI-CLPM.txt")



#### PREPARE RESULTS ####

#### collect outputs ####

## unconditional ##

riclpm.unc <- parameterEstimates(riclpm.fit, ci = T)



riclpm.unc$model <- "RI-CLPM"


riclpm.unc$type <- "Unconditional"




riclpm.unc$param <- NA
riclpm.unc$param <- ifelse(riclpm.unc$label == "a1", "Path a1 \n (Satisfaction --> Satisfaction)", riclpm.unc$param)
riclpm.unc$param <- ifelse(riclpm.unc$label == "a2", "Path a2 \n (Income --> Income)", riclpm.unc$param)
riclpm.unc$param <- ifelse(riclpm.unc$label == "c1", "Path c1 \n (Satisfaction --> Income)", riclpm.unc$param)
riclpm.unc$param <- ifelse(riclpm.unc$label == "c2", "Path c2 \n (Income --> Satisfaction)", riclpm.unc$param)
riclpm.unc$param <- ifelse(riclpm.unc$label == "cort1", "Initial Correlation \n (Satisfaction <--> Income)", riclpm.unc$param)
riclpm.unc$param <- ifelse(riclpm.unc$label == "cort2", "Within-Time Correlation (T2)", riclpm.unc$param)
riclpm.unc$param <- ifelse(riclpm.unc$label == "cort3", "Within-Time Correlation (T3)", riclpm.unc$param)
riclpm.unc$param <- ifelse(riclpm.unc$label == "cort4", "Within-Time Correlation (T4)", riclpm.unc$param)
riclpm.unc$param <- ifelse(riclpm.unc$label == "cori", "Correlation Between Latent Intercepts", riclpm.unc$param) 
riclpm.unc <- subset(riclpm.unc, !is.na(riclpm.unc$param))

## tic ##

riclpm.tic <- parameterEstimates(riclpm_tic.fit, ci = T)

riclpm.tic$model <- "RI-CLPM"

riclpm.tic$type <- "Time-Invariant Covariate (constant effect)"

riclpm.tic$param <- NA
riclpm.tic$param <- ifelse(riclpm.tic$label == "a1", "Path a1 \n (Satisfaction --> Satisfaction)", riclpm.tic$param)
riclpm.tic$param <- ifelse(riclpm.tic$label == "a2", "Path a2 \n (Income --> Income)", riclpm.tic$param)
riclpm.tic$param <- ifelse(riclpm.tic$label == "c1", "Path c1 \n (Satisfaction --> Income)", riclpm.tic$param)
riclpm.tic$param <- ifelse(riclpm.tic$label == "c2", "Path c2 \n (Income --> Satisfaction)", riclpm.tic$param)
riclpm.tic$param <- ifelse(riclpm.tic$label == "cort1", "Initial Correlation \n (Satisfaction <--> Income)", riclpm.tic$param)
riclpm.tic$param <- ifelse(riclpm.tic$label == "cort2", "Within-Time Correlation (T2)", riclpm.tic$param)
riclpm.tic$param <- ifelse(riclpm.tic$label == "cort3", "Within-Time Correlation (T3)", riclpm.tic$param)
riclpm.tic$param <- ifelse(riclpm.tic$label == "cort4", "Within-Time Correlation (T4)", riclpm.tic$param)
riclpm.tic$param <- ifelse(riclpm.tic$label == "cori", "Correlation Between Latent Intercepts", riclpm.tic$param) 
riclpm.tic <- subset(riclpm.tic, !is.na(riclpm.tic$param))

## tic with varying effects ##
riclpm.tic_var <- parameterEstimates(riclpm_tic_var.fit, ci = T)

riclpm.tic_var$model <- "RI-CLPM"

riclpm.tic_var$type <- "Time-Invariant Covariate (varying effect)"

riclpm.tic_var$param <- NA
riclpm.tic_var$param <- ifelse(riclpm.tic_var$label == "a1", "Path a1 \n (Satisfaction --> Satisfaction)", riclpm.tic_var$param)
riclpm.tic_var$param <- ifelse(riclpm.tic_var$label == "a2", "Path a2 \n (Income --> Income)", riclpm.tic_var$param)
riclpm.tic_var$param <- ifelse(riclpm.tic_var$label == "c1", "Path c1 \n (Satisfaction --> Income)", riclpm.tic_var$param)
riclpm.tic_var$param <- ifelse(riclpm.tic_var$label == "c2", "Path c2 \n (Income --> Satisfaction)", riclpm.tic_var$param)
riclpm.tic_var$param <- ifelse(riclpm.tic_var$label == "cort1", "Initial Correlation \n (Satisfaction <--> Income)", riclpm.tic_var$param)
riclpm.tic_var$param <- ifelse(riclpm.tic_var$label == "cort2", "Within-Time Correlation (T2)", riclpm.tic_var$param)
riclpm.tic_var$param <- ifelse(riclpm.tic_var$label == "cort3", "Within-Time Correlation (T3)", riclpm.tic_var$param)
riclpm.tic_var$param <- ifelse(riclpm.tic_var$label == "cort4", "Within-Time Correlation (T4)", riclpm.tic_var$param)
riclpm.tic_var$param <- ifelse(riclpm.tic_var$label == "cori", "Correlation Between Latent Intercepts", riclpm.tic_var$param) 
riclpm.tic_var <- subset(riclpm.tic_var, !is.na(riclpm.tic_var$param))


## tvc ##
riclpm.tvc <- parameterEstimates(riclpm_tvc.fit, ci = T)

riclpm.tvc$model <- "RI-CLPM"

riclpm.tvc$type <- "Time-Varying Covariate"

riclpm.tvc$param <- NA
riclpm.tvc$param <- ifelse(riclpm.tvc$label == "a1", "Path a1 \n (Satisfaction --> Satisfaction)", riclpm.tvc$param)
riclpm.tvc$param <- ifelse(riclpm.tvc$label == "a2", "Path a2 \n (Income --> Income)", riclpm.tvc$param)
riclpm.tvc$param <- ifelse(riclpm.tvc$label == "c1", "Path c1 \n (Satisfaction --> Income)", riclpm.tvc$param)
riclpm.tvc$param <- ifelse(riclpm.tvc$label == "c2", "Path c2 \n (Income --> Satisfaction)", riclpm.tvc$param)
riclpm.tvc$param <- ifelse(riclpm.tvc$label == "cort1", "Initial Correlation \n (Satisfaction <--> Income)", riclpm.tvc$param)
riclpm.tvc$param <- ifelse(riclpm.tvc$label == "cort2", "Within-Time Correlation (T2)", riclpm.tvc$param)
riclpm.tvc$param <- ifelse(riclpm.tvc$label == "cort3", "Within-Time Correlation (T3)", riclpm.tvc$param)
riclpm.tvc$param <- ifelse(riclpm.tvc$label == "cort4", "Within-Time Correlation (T4)", riclpm.tvc$param)
riclpm.tvc$param <- ifelse(riclpm.tvc$label == "cori", "Correlation Between Latent Intercepts", riclpm.tvc$param) 
riclpm.tvc <- subset(riclpm.tvc, !is.na(riclpm.tvc$param))


#### combine ####
riclpm.res <- rbind(riclpm.unc, riclpm.tic, riclpm.tic_var, riclpm.tvc)

results <- riclpm.res



#### tables ####

library(xtable)


## ri-clpm ##
riclpm.tab <- data.frame(riclpm.unc[, c("label", "est", "ci.lower", "ci.upper", "pvalue")]
                         , riclpm.tic[, c("est", "ci.lower", "ci.upper", "pvalue")]
                         , riclpm.tic_var[, c("est", "ci.lower", "ci.upper", "pvalue")]
                         , riclpm.tvc[, c("est", "ci.lower", "ci.upper", "pvalue")]
)

print(xtable(riclpm.tab, digits = c(1, 1, 2, 2, 2, 3, 2, 2, 2, 3, 2, 2, 2, 3, 2, 2, 2, 3))
      , include.rownames = F
      , include.colnames = F
      , booktabs = T
      , only.contents = T
)




#### plot ####

library(ggplot2)

## paths ##
paths <- subset(results, label == "a1" | label == "a2" | label == "c1" | label == "c2")

g <- ggplot(data = paths, aes(x = factor(model, level = c("CLPM", "RI-CLPM", "LCM-SR"))
                              , y = est
                              , group = factor(type, level = c("Time-Varying Covariate", "Time-Invariant Covariate (varying effect)", "Time-Invariant Covariate (constant effect)", "Unconditional" )))) + xlab("Model") + ylab("Unstandardized Estimate") + theme_bw()
g <- g + geom_hline(aes(yintercept = 0), lty = "dotted", lwd = 0.2)
g <- g + geom_point(aes(pch = type), position = position_dodge(width = -0.5), size = 3) + 
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), width = 0.2, position = position_dodge(width = -0.5))
g <- g + facet_wrap(~param, scales = "fixed")
# g <- g + guides(shape = guide_legend(reverse = T))
g <- g + theme(legend.position = "top"
               , legend.title = element_blank()
               , legend.text = element_text(size = 16)
               , axis.text = element_text(size = 16, color = "black")
               , axis.title = element_text(size = 18, face = "bold")
               , axis.title.x = element_blank()
               , strip.text = element_text(size = 16, face = "bold")
               , panel.spacing = unit(0.5, "cm")
)
g <- g + guides(shape = guide_legend(nrow = 2, byrow = TRUE, reverse = T))
g
ggsave(g, file = "./../pics/paths.pdf", width = 15, height = 20)
ggsave(g, file = "./../pics/paths.jpg", width = 15, height = 20, dpi = 300)


## correlations ##
cors <- subset(results, label == "cori" | label == "cors" | label == "cort1" | label == "cort2" | label == "cort3" | label == "cort4" | label == "core" | label == "cort" | label == "corixsy" | label == "coriysx")


g <- ggplot(data = cors, aes(x = factor(model, level = c("CLPM", "RI-CLPM", "LCM-SR")), y = est, group = factor(type, level = c("Time-Varying Covariate", "Time-Invariant Covariate (varying effect)", "Time-Invariant Covariate (constant effect)", "Unconditional" )))) + xlab("Model") + ylab("Correlation") + theme_bw()
g <- g + geom_hline(aes(yintercept = 0), lty = "dotted", lwd = 0.2)
g <- g + geom_point(aes(pch = type), position = position_dodge(width = -0.5), size = 3) + 
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), width = 0.2, position = position_dodge(width = -0.5))
g <- g + facet_wrap(~factor(param, level = c("Initial Correlation \n (Satisfaction <--> Income)", "Within-Time Correlation (T2)", "Within-Time Correlation (T3)", "Within-Time Correlation (T4)", "Correlation Between Latent Intercepts", "Correlation Between Latent Slopes", "Correlation Between Intercept Satisfaction \n and Slope Income", "Correlation Between Intercept Income \n and Slope Satisfaction")), scales = "fixed")
g <- g + theme(legend.position = "top"
               , legend.title = element_blank()
               , legend.text = element_text(size = 16)
               , axis.text = element_text(size = 16, color = "black")
               , axis.title = element_text(size = 18, face = "bold")
               , axis.title.x = element_blank()
               , strip.text = element_text(size = 16, face = "bold")
               , panel.spacing = unit(0.5, "cm")
)
g <- g + guides(shape = guide_legend(nrow = 2, byrow = TRUE, reverse = T))
g
ggsave(g, file = "./../pics/cors.pdf", width = 15, height = 20)
ggsave(g, file = "./../pics/cors.jpg", width = 15, height = 20, dpi = 300)