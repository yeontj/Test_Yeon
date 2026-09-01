##########################################
#RI-CLPM##################################
#FRIEND SUPPORT X LONELINESS##############
#MASTER'S THESIS #########################
#START DATE: 08.28.2026. #################
#JUNGYEON SUH ############################
##########################################

library('lavaan')

#read data
alive = readRDS('thesisdata_alive.rds')

alive

RICLPM_3wave <- '
  # Create between components (random intercepts) 
  # Not average but latent variables
  
  RI_lon =~ 1*loneliness11.1 +
            1*loneliness11.2 +
            1*loneliness11.3

  RI_friend =~ 1*frsupport1 +
               1*frsupport2 +
               1*frsupport3
               
  ##########################
  # Within-person components
  # Single-indicator latent variable from observed
  
  wlon1 =~ 1*loneliness11.1
  wlon2 =~ 1*loneliness11.2
  wlon3 =~ 1*loneliness11.3

  wfriend1 =~ 1*frsupport1
  wfriend2 =~ 1*frsupport2
  wfriend3 =~ 1*frsupport3

  ##########################
  # Freely estimated lagged effects
  # Including auto-regression & cross-lagged effect
  
  wlon2 + wfriend2 ~ wlon1 + wfriend1
  wlon3 + wfriend3 ~ wlon2 + wfriend2
  
  ##########################
  # Initial within-person covariance
  wlon1 ~~ wfriend1

  ##########################
  # Innovation covariances / correlated change
  # after first wave
  # They are residual covariance as they are predicted from previous waves
  wlon2 ~~ wfriend2
  wlon3 ~~ wfriend3

  ##########################
  # Random-intercept variances and covariance
  
  RI_lon ~~ RI_lon
  RI_friend ~~ RI_friend
  RI_lon ~~ RI_friend

  # Within-person variances
  wlon1 ~~ wlon1
  wfriend1 ~~ wfriend1

  # Residual within-person variances
  wlon2 ~~ wlon2
  wfriend2 ~~ wfriend2
  wlon3 ~~ wlon3
  wfriend3 ~~ wfriend3
'

fit_free <- lavaan(
  RICLPM_3wave,
  data = alive,
  missing = "fiml", #full-information maximum likelihood
  estimator = "MLR", 
  meanstructure = TRUE, #estimate observed-variable mean/intercept by wave
  int.ov.free = TRUE #free estimation of observed-variable intercept
)

summary(fit_free, standardized = F)
summary(fit_free, standardized = T)


