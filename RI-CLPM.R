##########################################
#RI-CLPM##################################
#FRIEND SUPPORT X LONELINESS##############
#MASTER'S THESIS #########################
#START DATE: 08.28.2026. #################
#JUNGYEON SUH 


library('lavaan')

#read data
alive = read.csv('thesis_analysis.csv')

str(alive)

################################################
#BASE MODEL#####################################

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
  missing = "ML", #maximum likelihood
  meanstructure = TRUE, #estimate observed-variable mean/intercept by wave
  int.ov.free = TRUE #free estimation of observed-variable intercept
)

summary(fit_free, standardized = T,fit.measures=T)
##################################################
#IMPOSING CONSTRAINTS OVER TIME###################

RICLPM_constrained <- '

  # Between-person components: random intercepts
  RI_lon =~ 1*loneliness11.1 +
            1*loneliness11.2 +
            1*loneliness11.3

  RI_friend =~ 1*frsupport1 +
               1*frsupport2 +
               1*frsupport3


  # Within-person components
  wlon1 =~ 1*loneliness11.1
  wlon2 =~ 1*loneliness11.2
  wlon3 =~ 1*loneliness11.3

  wfriend1 =~ 1*frsupport1
  wfriend2 =~ 1*frsupport2
  wfriend3 =~ 1*frsupport3


  # Constrained autoregressive and cross-lagged paths
  wlon2 ~ a*wlon1 + b*wfriend1
  wfriend2 ~ c*wlon1 + d*wfriend1

  wlon3 ~ a*wlon2 + b*wfriend2
  wfriend3 ~ c*wlon2 + d*wfriend2


  # Constrained innovation covariance
  wlon2 ~~ cov*wfriend2
  wlon3 ~~ cov*wfriend3


  # Initial within-person covariance
  wlon1 ~~ wfriend1


  # Random-intercept variances and covariance
  RI_lon ~~ RI_lon
  RI_friend ~~ RI_friend
  RI_lon ~~ RI_friend


  # Initial within-person variances
  wlon1 ~~ wlon1
  wfriend1 ~~ wfriend1


  # Constrained residual within-person variances
  wlon2 ~~ vx*wlon2
  wlon3 ~~ vx*wlon3

  wfriend2 ~~ vy*wfriend2
  wfriend3 ~~ vy*wfriend3


  # Constrained grand means
  loneliness11.1 +
  loneliness11.2 +
  loneliness11.3 ~ mx*1

  frsupport1 +
  frsupport2 +
  frsupport3 ~ my*1
'

fit_constrained = lavaan(RICLPM_constrained,
                         data = alive,
                         missing = 'ML',
                         meanstructure = T,
                         int.ov.free = T
                         )
summary(fit_constrained, standardized = T, fit.measures = T)

library(restriktor)

H1 <- "abs(b2) < abs(c2); abs(b3) < abs(c3); 
       abs(b4) < abs(c4); abs(b5) < abs(c5)" 

GORICA_free <- goric(
  RICLPM_3wave, 
  standardized = TRUE,
  hypotheses = list(H1ws.l = H1), 
  comparison = "complement", # Test informative hypothesis versus its complement
  type = "gorica"
)

GORICA_free


