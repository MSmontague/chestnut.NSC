######################################################################
######    M Montague   ###   Chestnut NSC   ###   July 2020   ########
######################################################################
##############################  modified July 7, 2021  ###############

setwd("C:/Users/madel/OneDrive - purdue.edu/publishing/R_scripts_data")

####### Packages, functions & colors                       ############

library(ggplot2)
library(ggpubr)
library(gridExtra)
library(plyr)
library(dplyr)
library(tidyr)
library(nlme)
library(lme4)
library(car)
library(emmeans)
library(sjstats)

# if having trouble---try install.packages("XXX", type = "binary")

# functions for transformations 
asinTransform <- function(p) { asin(sqrt(p)) }
in.asin = function(x){(sin(x)**2)}

# colors for figures
colorset = c("black","orangered","blue")  #  prev. "darkred"
colorset2 = c("blacK","red","chocolate","blue") 

#### 1. Seasonal analysis of NSC concentration data        ############

##  A. Data  

NSC <- read.csv("P1_NSC_data.csv", header=TRUE)
  # N_fraction and C_fraction are from COSCO.  C_fraction is the total C (NSC+structural)
  NSC$tree <- factor(NSC$tree)
  NSC$tissue <- factor(NSC$tissue, levels = c("Leaves", "Twigs", "Branches", "Inner bark", "Xylem", "Root collar", "Coarse roots", "Fine roots"))
  NSC$collection2 <- revalue(NSC$collection, c("1-leaf out"="1-leaf out", "7-leaf out II"="1-leaf out"))
  NSC$collection3 <- revalue(NSC$collection2, c("1-leaf out"="Leaf-out","2-shoot expansion"="Shoot expansion", "3-bud set"="Bud set", "4-leaf fall"="Leaf fall","5-dormancy"="Dormancy" , "6-sap flow" ="Sap flow" ))
  NSC$starch_fraction2 <- ifelse(NSC$starch_fraction <= 0, 0.0000001, NSC$starch_fraction)
  NSC$NSC_Arcsin <- asin(sqrt(NSC$tot_NSC_fraction))
  NSC$starch_Arcsin <- asin(sqrt(NSC$starch_fraction2))
  NSC$sugar_Arcsin <- asin(sqrt(NSC$sugar_fraction))

  
## B. compare leaf-out 2019 (1-leaf out) vs. leaf out 2020 (7-leaf out II) 
  
  selected<-c("1-leaf out","7-leaf out II")
  LO <- NSC[ NSC$collection %in% selected,]
  table(LO$tissue, LO$collection) # one inner bark sample was lost from leaf out II (Table S1)
  
  # compare models differing in their inclusion of collection (LO 2019 vs. LO 2020)
  
  # starch
  lmer.starch.LO = lmer(starch_Arcsin ~  tissue * collection + (1|tree), data=LO)
  lmer.starch.LO_n = lmer(starch_Arcsin ~  tissue + (1|tree), data=LO)
  anova(lmer.starch.LO,lmer.starch.LO_n, test="LRT")
  # sugar
  lmer.sugar.LO = lmer(sugar_Arcsin ~  tissue * collection + (1|tree), data=LO)
  lmer.sugar.LO_n = lmer(sugar_Arcsin ~  tissue + (1|tree), data=LO)
  anova(lmer.sugar.LO,lmer.sugar.LO_n, test="LRT")
  # total NSC
  lmer.NSC.LO = lmer(NSC_Arcsin ~  tissue * collection + (1|tree), data=LO)
  lmer.NSC.LO_n = lmer(NSC_Arcsin ~  tissue + (1|tree), data=LO)
  anova(lmer.NSC.LO,lmer.NSC.LO_n, test="LRT")
  

## C. seasonal analysis
  
 # I. starch model  

  starch.lme <- lmer(starch_Arcsin ~  tissue * collection2 + (1|tree), data=NSC)  
  # these models will be rank-deficient because there is no NSC data for leaves in the winter months.  
  # expect to drop 6 coefficients.
  plot(starch.lme)
  qqnorm(resid(starch.lme))
  qqline(resid(starch.lme))
  summary(starch.lme)
  Anova(starch.lme, type="III")
  starch.lme.null <- lmer(starch_Arcsin ~  1 + (1|tree), data=NSC)
  anova(starch.lme, starch.lme.null, test="F")
  performance::r2(starch.lme) # : Nakagawa's R2: The marginal R2 considers only the variance of the fixed effects, while the conditional R2 takes both the fixed and random effects into account.
  # https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.2041-210x.2012.00261.x
  
  # visualize model fit
  NSC_starch <- NSC 
  NSC_starch$starch_fit <- in.asin(predict(starch.lme, type="response"))
  ggplot(data=NSC_starch) + 
    geom_point(aes(x = collection, y = starch_fraction)) +   # observed data
    geom_point(aes(x=collection2, y = starch_fit),col='red', size=2) +  # fitted values
    theme_bw() +
    labs(x = "Collection", y = "Starch %") +
    theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1)) +
    facet_wrap(~tissue) 

 # II. sugar model  
  
  sugar.lme <- lmer(sugar_Arcsin ~  tissue * collection2 + (1|tree), data=NSC)  
  plot(sugar.lme)
  qqnorm(resid(sugar.lme))
  qqline(resid(sugar.lme))
  summary(sugar.lme)
  Anova(sugar.lme, type="III")
  sugar.lme.null <- lmer(sugar_Arcsin ~  1 + (1|tree), data=NSC)
  anova(sugar.lme, sugar.lme.null, test="F")
  performance::r2(sugar.lme) 
  
  # visualize model fit
  NSC_sugar <- NSC 
  NSC_sugar$sugar_fit <- in.asin(predict(sugar.lme, type="response"))
  ggplot(data=NSC_sugar) + 
    geom_point(aes(x = collection, y = sugar_fraction)) +   # observed data
    geom_point(aes(x=collection2, y = sugar_fit),col='red', size=2) +  # fitted values
    theme_bw() +
    labs(x = "Collection", y = "sugar %") +
    theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1)) +
    facet_wrap(~tissue)
  
 # III. NSC model  
  
  NSC.lme <- lmer(NSC_Arcsin ~  tissue * collection2 + (1|tree) , data=NSC)  
  # these models will be rank-deficient because there is no NSC data for leaves in the winter months
  plot(NSC.lme)
  qqnorm(resid(NSC.lme))
  qqline(resid(NSC.lme))
  summary(NSC.lme)
  Anova(NSC.lme, type="III")
  NSC.lme.null <- lmer(NSC_Arcsin ~  1 + (1|tree), data=NSC)
  anova(NSC.lme, NSC.lme.null, test="F")
  performance::r2(NSC.lme) 
  
  # visualize model fit
  NSC_NSC <- NSC 
  NSC_NSC$NSC_fit <- in.asin(predict(NSC.lme, type="response"))
  ggplot(data=NSC_NSC) + 
    geom_point(aes(x = collection, y = tot_NSC_fraction)) +   # observed data
    geom_point(aes(x=collection2, y = NSC_fit),col='red', size=2) +  # fitted values
    theme_bw() +
    labs(x = "Collection", y = "NSC %") +
    theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1)) +
    facet_wrap(~tissue) 
  
  
## D. summarize mean & SD for starch, sugar & total NSC concentrations
  detach(package:plyr)
  concentrations_summary <- NSC %>% 
    group_by(collection2, tissue) %>% 
    summarize(mean_NSC = mean(tot_NSC_fraction),
              mean_starch = mean(starch_fraction),
              mean_sugar = mean(sugar_fraction),
              st_su = mean(starch_fraction/sugar_fraction),
              sd_NSC = sd(tot_NSC_fraction),
              sd_starch = sd(starch_fraction),
              sd_sugar = sd(sugar_fraction))
  concentrations_summary  ## Table S2  ##
  write.csv(concentrations_summary,"C:/Users/madel/OneDrive - purdue.edu/publishing/R_scripts_data//concentrations_summary.csv", row.names = FALSE)

  
#### 2. Chestnut allometry                                 ####

## A. Biomass data 
biomass <-na.omit(read.csv("biomass_clean_long_full.csv"))       # no fractions left out
  biomass$tissue <- factor(biomass$tissue, levels = c("Leaves", "Twigs", "Branches", "Inner bark", "Xylem", "Root collar", "Coarse roots", "Fine roots", "Total biomass"))

## B. Calculate root:shoot ratio
bb <- spread(biomass, tissue, biomass)
  bb$shoot = bb$Leaves  + bb$Twigs + bb$Branches + bb$'Inner bark' + bb$Xylem
  bb$root = bb$'Coarse roots' + bb$'Fine roots' + bb$'Root collar'
  bb$root_shoot = bb$root/bb$shoot
  
bb%>% 
  # group_by(tree) %>%  # remove comment to examine RSR & DBH for each tree
  summarize(mean_DBH = mean(DBH),
            sd_DBH   = sd(DBH),
            mean_RSR = mean(root_shoot),
            sd_RSR   = sd(root_shoot),
            se_RSR   = sd(root_shoot, na.rm=T)/sqrt(sum(!is.na(root_shoot))))
# visualize
ggplot() + theme_classic() +
  labs(y="Root:shoot ratio") +
  geom_smooth(data=bb, aes(x=DBH, y=root_shoot), method="lm", se=F) +
  geom_point(data=bb, aes(x=DBH, y=root_shoot)) +
  scale_y_continuous(limits=c(0,0.5))
# model
rs <- glm(data=bb, root_shoot~DBH)
  plot(rs)
  summary(rs)
  car::Anova(rs,type="III")

## C. Fit allometric glm (Table 2)
  
glm_gamma <- glmer(biomass ~ DBH  * tissue + Ht + (1|tree) , data=biomass, family=Gamma(link="log"),
                     glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
  plot(glm_gamma)
  qqnorm(resid(glm_gamma))
  qqline(resid(glm_gamma))
  summary(glm_gamma)
  performance::r2(glm_gamma) # : Nakagawa's R2: The marginal R2 considers only the variance of the fixed effects, while the conditional R2 takes both the fixed and random effects into account.
  # https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.2041-210x.2012.00261.x
  glm_gamma_null <- glmer(biomass ~ 1 + (1|tree), data=bio, family=Gamma(link="log"))
  anova(glm_gamma, glm_gamma_null, test="Chisq")
  Anova(glm_gamma, type="III")
  
## D. Predict organ biomass for all trees in the plantation

# Census DBH and height data for our chestnut stand
# file has space for predicted organ biomass and predicted NSC concentration throughout the year
C_2m <- read.csv("C_2m_biomass - with total biomass.csv")
  C_2m$tag <- factor(C_2m$tag)
  levels(C_2m$tag)   
  
  bio <- biomass
  
  Base <- ggplot() + 
    geom_point(aes(x = DBH, y = biomass, color=tissue), data=biomass) + 
    scale_color_discrete() + 
    facet_wrap(~tissue) 
  Base
  
  # fit a glm (allometric model)
  glm_gamma <- glmer(biomass ~ DBH  * tissue + Ht + (1|tree), data=bio, family=Gamma(link="log"),
                     glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
  # get fitted values
  bio$bio_fit <- predict(glm_gamma, type="response")
 
  #plot fitted values
  Base + geom_point(aes(x=DBH, y = bio_fit), data=bio) + facet_wrap(~tissue)
  
  # introduce new data 
  rm(nd)
  nd <- C_2m
  nd$bio_pred <- predict(glm_gamma, newdata=nd, type="response", re.form=NA)   # re.form = do not include random effects.  Because these are different trees!
  
  # plot
  Base + 
    geom_point(aes(x=DBH, y = bio_fit), data=bio) + 
    geom_point(aes(x=DBH, y = bio_pred), data = nd, size = 0.8, col="brown", alpha = 0.2) +
    facet_wrap(~tissue)
  
  # function to predict confidence intervals. 
  # from Ben Bolker: https://bbolker.github.io/mixedmodels-misc/ecostats_chap.html
  easyPredCI <- function(model,newdata=NULL,alpha=0.05) {
    ## baseline prediction, on the linear predictor (logit) scale:
    pred0 <- predict(model,re.form=NA,newdata=newdata)
    ## fixed-effects model matrix for new data
    X <- model.matrix(formula(model,fixed.only=TRUE)[-2],newdata)
    beta <- fixef(model) ## fixed-effects coefficients
    V <- vcov(model)     ## variance-covariance matrix of beta
    pred.se <- sqrt(diag(X %*% V %*% t(X))) ## std errors of predictions
    ## inverse-link function
    linkinv <- family(model)$linkinv
    ## construct 95% Normal CIs on the link scale and
    ##  transform back to the response (probability) scale:
    crit <- -qnorm(alpha/2)
    linkinv(cbind(resp_lwr=pred0-crit*pred.se,
                  resp_upr=pred0+crit*pred.se,
                  fit_resp=pred0))
  }
  
  pred.CI <- easyPredCI(glm_gamma,nd)
  
  nd <- cbind(nd,pred.CI)
  
  # plot! confidence intervals get v. large >16cm DBH
  Base + 
    geom_point(aes(x=DBH, y = bio_fit), data=bio) + 
    geom_point(aes(x=DBH, y = bio_pred), data = nd, size = 0.8, col="brown", alpha = 0.2) +
    geom_point(aes(x=DBH, y = fit_resp), data = nd, size = 1, col="darkgoldenrod", alpha = 0.5) +
    geom_errorbar(aes(x = DBH, ymax = resp_upr, ymin = resp_lwr), data=nd, col="brown") +
    facet_wrap(~tissue, scales="free", nrow=2) + 
    theme(legend.position = "none") + ylab("Biomass (kg)")
  
  # We approximated the standard error by taking the mean of the upper and lower CI length /1.96
  se_upr <- (nd$resp_upr - nd$fit_resp)/1.96 
  se_lwr <- (nd$fit_resp - nd$resp_lwr)/1.96 
  se <- (se_upr + se_lwr)/2
  
  cbind(nd$fit_resp,se_upr, se_lwr,se)
  
  nd$tissue_biomass <- nd$fit_resp
  nd$tissue_biomass_lci <- nd$resp_lwr
  nd$tissue_biomass_uci <- nd$resp_upr
  nd$tissue_biomass_se <- se   
  nd
  
  C_2m$tissue_biomass      <- nd$tissue_biomass
  C_2m$tissue_biomass_lci  <- nd$tissue_biomass_lci
  C_2m$tissue_biomass_uci  <- nd$tissue_biomass_uci
  C_2m$tissue_biomass_se   <- nd$tissue_biomass_se
  
  C_2m
  
  write.csv(C_2m,"C:/Users/madel/OneDrive - purdue.edu/publishing/R_scripts_data//C_2m_out.csv", row.names = FALSE)
  

#### 3. NSC pool size calculation                          ####  

C_2m <- C_2m[-which(C_2m$tissue=="Total biomass"),]
  
## A. Estimate NSC concentrations                         
# for all trees in the plantation
# I. Starch 
  
NSC_starch <- NSC
  NSC_starch$collection <- NSC_starch$collection2

starchBase <- ggplot(data=NSC_starch) + 
  geom_point(aes(x = collection, y = starch_fraction)) +
  scale_color_discrete() + 
  theme_bw() +
  labs(x = "Collection", y = "Starch %") +
  theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1)) +
  facet_wrap(~tissue) 
starchBase

# fit a glm
starch.glm <- lmer(starch_Arcsin ~  tissue * collection + (1|tree), data=NSC_starch)
  
# fit values
NSC_starch$starch_fit <- in.asin(predict(starch.glm, type="response"))
  
# visualize
starchBase + 
  geom_point(aes(x=collection, y = starch_fit),col="blue", size=3, data=NSC_starch) + # fitted 
  
# introduce new data: C_2m (every tree in the plantation)
  rm(nd_st)
  nd_st <- C_2m
  nd_st$starch_pred <- in.asin(predict(starch.glm, newdata=nd_st, type="response", re.form=NA))
  
# plot
starchBase + 
  geom_point(aes(x=collection, y = starch_fit),col="blue", size=3, data=NSC_starch, alpha = 0.5) +     # interaction
  geom_point(aes(x=collection, y = starch_pred), data = nd_st, size = 1.5, col="green") +
  facet_wrap(~tissue)

### bootMer takes a fitted merMod object then the function that takes a fitted merMod object and returns the thing we want as a vector. 
  # 1000 simulations takes a minute
  bigBoot_starch <- bootMer(starch.glm, nsim=1000, 
                            FUN=function(x)predict(x, newdata=nd_st, re.form=NA, type="response"),
                            parallel = "multicore", ncpus = 4)
  bigBoot_starch

starch_boot_ci <- t(apply(bigBoot_starch$t,2,quantile,c(0.025,0.975),na.rm=TRUE))
  nd_st$starch_lci <- in.asin(starch_boot_ci[,1])
  nd_st$starch_uci <- in.asin(starch_boot_ci[,2])
  
# PLOT! 
starchBase + 
    geom_point(aes(x=collection, y = starch_fit), data=NSC_starch) + 
    geom_point(aes(x=collection, y = starch_pred), data = nd_st, size = 1.3, col="red", alpha = 0.5) +
    geom_errorbar(aes(x = collection, ymax = starch_lci, ymin = starch_uci), data=nd_st, col="brown") +
    facet_wrap(~tissue)
  
  se_upr <- (nd_st$starch_uci - nd_st$starch_pred)/1.96 
  se_lwr <- (nd_st$starch_pred - nd_st$starch_lci)/1.96 
  se <- (se_upr + se_lwr)/2
  sq <- (se/nd_st$fit_resp)^2
  
  nd_st$starch_fraction <- nd_st$starch_pred
  nd_st$starch_fraction_lci <- nd_st$starch_lci
  nd_st$starch_fraction_uci <- nd_st$starch_uci
  nd_st$starch_fraction_se <- se   # this se is really the average of the two se BACKTRANSFORMED.  b/c the transformation stretched out the data!
  
  
# II. Sugar 
  
  NSC_sugar <- NSC
  NSC_sugar$collection <- NSC_sugar$collection2
  
  sugarBase <- ggplot(data=NSC_sugar) + 
    geom_point(aes(x = collection, y = sugar_fraction)) +
    scale_color_discrete() + 
    theme_bw() +
    labs(x = "Collection", y = "sugar %") +
    theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1)) +
    facet_wrap(~tissue) 
  sugarBase
  
  # fit a glm
  sugar.glm <- lmer(sugar_Arcsin ~  tissue * collection + (1|tree), data=NSC_sugar)
  
  # fit values
  NSC_sugar$sugar_fit <- in.asin(predict(sugar.glm, type="response"))
  
  # visualize
  sugarBase + 
    geom_point(aes(x=collection, y = sugar_fit),col="blue", size=3, data=NSC_sugar) + # fitted 
    
  # introduce new data: C_2m (every tree in the plantation)
  rm(nd_su)
  nd_su <- C_2m
  nd_su$sugar_pred <- in.asin(predict(sugar.glm, newdata=nd_su, type="response", re.form=NA))
  
  # plot
  sugarBase + 
    geom_point(aes(x=collection, y = sugar_fit),col="blue", size=3, data=NSC_sugar, alpha = 0.5) +     # interaction
    geom_point(aes(x=collection, y = sugar_pred), data = nd_su, size = 1.5, col="green") +
    facet_wrap(~tissue)
  
  bigBoot_sugar <- bootMer(sugar.glm, nsim=1000, 
                            FUN=function(x)predict(x, newdata=nd_su, re.form=NA, type="response"),
                            parallel = "multicore", ncpus = 4)
  bigBoot_sugar
  
  sugar_boot_ci <- t(apply(bigBoot_sugar$t,2,quantile,c(0.025,0.975),na.rm=TRUE))
  nd_su$sugar_lci <- in.asin(sugar_boot_ci[,1])
  nd_su$sugar_uci <- in.asin(sugar_boot_ci[,2])
  
  # PLOT! 
  sugarBase + 
    geom_point(aes(x=collection, y = sugar_fit), data=NSC_sugar) + 
    geom_point(aes(x=collection, y = sugar_pred), data = nd_su, size = 1.3, col="red", alpha = 0.5) +
    geom_errorbar(aes(x = collection, ymax = sugar_lci, ymin = sugar_uci), data=nd_su, col="brown") +
    facet_wrap(~tissue)
  
  se_upr <- (nd_su$sugar_uci - nd_su$sugar_pred)/1.96 
  se_lwr <- (nd_su$sugar_pred - nd_su$sugar_lci)/1.96 
  se <- (se_upr + se_lwr)/2
  sq <- (se/nd_su$fit_resp)^2
  
  nd_su$sugar_fraction <- nd_su$sugar_pred
  nd_su$sugar_fraction_lci <- nd_su$sugar_lci
  nd_su$sugar_fraction_uci <- nd_su$sugar_uci
  nd_su$sugar_fraction_se <- se   # this se is really the average of the two se BACKTRANSFORMED.  b/c the transformation stretched out the data!
  
  
# III. NSC  
  
  NSC_NSC <- NSC
  NSC_NSC$collection <- NSC_NSC$collection2
  
  NSCBase <- ggplot(data=NSC_NSC) + 
    geom_point(aes(x = collection, y = tot_NSC_fraction)) +
    scale_color_discrete() + 
    theme_bw() +
    labs(x = "Collection", y = "NSC %") +
    theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1)) +
    facet_wrap(~tissue) 
  NSCBase
  
  # fit a glm
  NSC.glm <- lmer(NSC_Arcsin ~  tissue * collection + (1|tree), data=NSC_NSC)
  
  # fit values
  NSC_NSC$NSC_fit <- in.asin(predict(NSC.glm, type="response"))
  
  # visualize
  NSCBase + 
    geom_point(aes(x=collection, y = NSC_fit),col="blue", size=3, data=NSC_NSC) + # fitted 
    
    # introduce new data: C_2m (every tree in the plantation)
    rm(nd_NSC)
  nd_NSC <- C_2m
  nd_NSC$NSC_pred <- in.asin(predict(NSC.glm, newdata=nd_NSC, type="response", re.form=NA))
  
  # plot
  NSCBase + 
    geom_point(aes(x=collection, y = NSC_fit),col="blue", size=3, data=NSC_NSC, alpha = 0.5) +     # interaction
    geom_point(aes(x=collection, y = NSC_pred), data = nd_NSC, size = 1.5, col="green") +
    facet_wrap(~tissue)
  
  bigBoot_NSC <- bootMer(NSC.glm, nsim=1000, 
                            FUN=function(x)predict(x, newdata=nd_NSC, re.form=NA, type="response"),
                            parallel = "multicore", ncpus = 4)
  bigBoot_NSC
  
  NSC_boot_ci <- t(apply(bigBoot_NSC$t,2,quantile,c(0.025,0.975),na.rm=TRUE))
  nd_NSC$NSC_lci <- in.asin(NSC_boot_ci[,1])
  nd_NSC$NSC_uci <- in.asin(NSC_boot_ci[,2])
  
  # PLOT! 
  NSCBase + 
    geom_point(aes(x=collection, y = NSC_fit), data=NSC_NSC) + 
    geom_point(aes(x=collection, y = NSC_pred), data = nd_NSC, size = 1.3, col="red", alpha = 0.5) +
    geom_errorbar(aes(x = collection, ymax = NSC_lci, ymin = NSC_uci), data=nd_NSC, col="brown") +
    facet_wrap(~tissue)
  
  se_upr <- (nd_NSC$NSC_uci - nd_NSC$NSC_pred)/1.96 
  se_lwr <- (nd_NSC$NSC_pred - nd_NSC$NSC_lci)/1.96 
  se <- (se_upr + se_lwr)/2
  sq <- (se/nd_NSC$fit_resp)^2
  
  nd_NSC$NSC_fraction <- nd_NSC$NSC_pred
  nd_NSC$NSC_fraction_lci <- nd_NSC$NSC_lci
  nd_NSC$NSC_fraction_uci <- nd_NSC$NSC_uci
  nd_NSC$NSC_fraction_se <- se   # this se is really the average of the two se BACKTRANSFORMED.  b/c the transformation stretched out the data!
  
## B. Create pool sizes 
  
  # Multiply est. concentrations x biomass. Propagation of error for pool SE 
  
  # Starch
  C_2m$starch_fraction     <- nd_st$starch_fraction
  C_2m$starch_fraction_se  <- nd_st$starch_fraction_se 
  C_2m$starch_fraction_lci <- nd_st$starch_fraction_lci
  C_2m$starch_fraction_uci <- nd_st$starch_fraction_uci
  C_2m$starch_fraction_lci2 <- ifelse( C_2m$starch_fraction_lci < 0, 0.0000001,  C_2m$starch_fraction_lci)

  C_2m$starch_pool     <- C_2m$starch_fraction*C_2m$tissue_biomass
  C_2m$starch_pool_se  <- (sqrt((C_2m$starch_fraction_se/C_2m$starch_fraction)^2 + (C_2m$tissue_biomass_se/C_2m$tissue_biomass)^2))*C_2m$starch_pool 
  C_2m$starch_pool_lci <- C_2m$starch_pool - 1.96*C_2m$starch_pool_se
  C_2m$starch_pool_uci <- C_2m$starch_pool + 1.96*C_2m$starch_pool_se
  C_2m$starch_pool_lci2 <- ifelse( C_2m$starch_pool_lci < 0, 0.0000001,  C_2m$starch_pool_lci)
  
  # Sugar
  C_2m$sugar_fraction     <- nd_su$sugar_fraction
  C_2m$sugar_fraction_se  <- nd_su$sugar_fraction_se 
  C_2m$sugar_fraction_lci <- nd_su$sugar_fraction_lci
  C_2m$sugar_fraction_uci <- nd_su$sugar_fraction_uci
  
  C_2m$sugar_pool     <- C_2m$sugar_fraction*C_2m$tissue_biomass
  C_2m$sugar_pool_se  <- (sqrt((C_2m$sugar_fraction_se/C_2m$sugar_fraction)^2 + (C_2m$tissue_biomass_se/C_2m$tissue_biomass)^2))*C_2m$sugar_pool
  C_2m$sugar_pool_lci <- C_2m$sugar_pool - 1.96*C_2m$sugar_pool_se
  C_2m$sugar_pool_uci <- C_2m$sugar_pool + 1.96*C_2m$sugar_pool_se
  
  #NSC  
  C_2m$NSC_fraction     <- nd_NSC$NSC_fraction
  C_2m$NSC_fraction_se  <- nd_NSC$NSC_fraction_se 
  C_2m$NSC_fraction_lci <- nd_NSC$NSC_fraction_lci
  C_2m$NSC_fraction_uci <- nd_NSC$NSC_fraction_uci 
  
  C_2m$NSC_pool     <- C_2m$NSC_fraction*C_2m$tissue_biomass
  C_2m$NSC_pool_se  <- (sqrt((C_2m$NSC_fraction_se/C_2m$NSC_fraction)^2 + (C_2m$tissue_biomass_se/C_2m$tissue_biomass)^2))*C_2m$NSC_pool 
  C_2m$NSC_pool_lci <- C_2m$NSC_pool - 1.96*C_2m$NSC_pool_se
  C_2m$NSC_pool_uci <- C_2m$NSC_pool + 1.96*C_2m$NSC_pool_se
  
write.csv(C_2m,"C:/Users/madel/OneDrive - purdue.edu/publishing/R_scripts_data//C_2m_out_NSCfull.csv", row.names = FALSE)
 
#### 4. Identify median tree pools                         ####       

# Use the pools from median-sized tree to illustrate 
 
detach(package:plyr)

# want a list of each tree (only once)
C_2m <- read.csv("C_2m_out_NSCfull.csv")
C_2m_trees <- C_2m[which(C_2m$collection=="4-leaf fall" & C_2m$tissue=="Coarse roots"),] 

summary(C_2m_trees$Ht)
summary(C_2m_trees$DBH)
# median-sized tree has a DBH of 10.98cm and a height of 11m
# tree #795 is closest: DBH = 10.6cm and Height = 10.3m

C_2m_med1 <- C_2m[which(C_2m$tag == "795"),]

C_2m_med <- C_2m_med1 %>% 
  group_by(tissue, collection, Julian.date) %>% 
  summarize(NSC_pool = NSC_pool,
            starch_pool = starch_pool,
            sugar_pool = sugar_pool,
            NSC_pool_se = NSC_pool_se,
            starch_pool_se = starch_pool_se,
            sugar_pool_se = sugar_pool_se,
            tissue_biomass = tissue_biomass,
            tissue_biomass_se = tissue_biomass_se)
C_2m_med 
med.pools.summary <- as.data.frame(C_2m_med)
write.csv(med.pools.summary,"C:/Users/madel/OneDrive - purdue.edu/publishing/R_scripts_data//med.pools.summary.csv", row.names = FALSE)

# create whole-tree NSC pools for the median tree
med_wt <- C_2m_med %>% 
  group_by(collection, Julian.date) %>% 
  summarize(NSC_pool = sum(NSC_pool),
            starch_pool = sum(starch_pool),
            sugar_pool = sum(sugar_pool),
            tissue_biomass=sum(tissue_biomass),
            NSC_pool_se = sum(sqrt(NSC_pool_se^2)),
            starch_pool_se = sum(sqrt(starch_pool_se^2)),
            sugar_pool_se = sum(sqrt(sugar_pool_se^2)),
            tissue_biomass_se=sum(sqrt(sugar_pool_se^2)))
med_wt
med_wt$tissue = "Whole tree"
selected<-c("2-shoot expansion","4-leaf fall")
med_wt <- med_wt[med_wt$collection %in% selected,]

med_wt <- rbind(C_2m_med, med_wt)
med_wt
med_wt$tissue <- factor(med_wt$tissue, levels =  c("Leaves", "Twigs", "Branches", "Inner bark", "Xylem","Root collar", "Coarse roots", "Fine roots", "Whole tree"))
med_wt


#### 5. Delta NSC                                          ####

# calculation
detach(package:plyr)
remove(t)    

# create C_2m_med in Step 4 above
t <- C_2m_med[-which(C_2m_med$tissue=="Leaves"),] %>%
  group_by(tissue) %>% 
  summarise(NSC_flux = max(NSC_pool) - min(NSC_pool),
            starch_flux = max(starch_pool) - min(starch_pool),
            sugar_flux = max(sugar_pool) - min(sugar_pool),
            NSC_flux_se = sqrt(max(NSC_pool_se)^2 + min(NSC_pool_se)^2),
            starch_flux_se = sqrt(max(starch_pool_se)^2 + min(starch_pool_se)^2),
            sugar_flux_se = sqrt(max(sugar_pool_se)^2 + min(sugar_pool_se)^2))

# modelling 

####### Fig. 2.  Concentration                             ####
  
  detach(package:plyr)
  
  nu <- NSC %>% 
    group_by(collection, tissue, Julian.date) %>% 
    summarize(mean_NSC = mean(tot_NSC_fraction),
              mean_starch = mean(starch_fraction),
              mean_sugar = mean(sugar_fraction),
              st_su = mean(starch_fraction/sugar_fraction),
              sd_NSC = sd(tot_NSC_fraction),
              sd_starch = sd(starch_fraction),
              sd_sugar = sd(sugar_fraction),
              se_NSC = sd(tot_NSC_fraction)/sqrt(n()),
              se_starch = sd(starch_fraction)/sqrt(n()),
              se_sugar = sd(sugar_fraction)/sqrt(n()))
  nu
  nu <- nu[-which(nu$Julian.date==111),]
  
  ggplot() + theme_classic() + ylab("Concentration (%)") + 
    geom_point(aes(x=Julian.date, y=sugar_fraction, col="Sugar"), 
               data=NSC, alpha= 0.7, size=1.3, shape=16, position="jitter", show.legend = FALSE) +
    geom_point(aes(x=Julian.date, y=starch_fraction, col="Starch"), 
               data=NSC, alpha= 0.7, size=1.3, shape=16, position="jitter",show.legend = FALSE) + 
    geom_point(aes(x=Julian.date, y=tot_NSC_fraction, col="Total NSC"), 
               data=NSC, alpha= 0.7, size=1.3, shape=16, position="jitter", show.legend = FALSE) +
    geom_line(aes(x=Julian.date, y=mean_NSC, col= "Total NSC"), 
              data=nu, size=1.2, show.legend = FALSE) +
    geom_line(aes(x=Julian.date, y=mean_sugar, col= "Sugar"), 
              data=nu, size=1.2, show.legend = FALSE) +
    geom_line(aes(x=Julian.date, y=mean_starch, col= "Starch"), 
              data=nu, size=1.2, show.legend = FALSE) +
    geom_ribbon(aes(x=Julian.date, ymin=mean_NSC-1.96*se_NSC, ymax=mean_NSC+1.96*se_NSC, fill= "Total NSC"), 
                data=nu, alpha=0.3) +
    geom_ribbon(aes(x=Julian.date, ymin=mean_starch-1.96*se_starch, ymax=mean_starch+1.96*se_starch, fill= "Starch"), 
                data=nu,  alpha=0.4) +
    geom_ribbon(aes(x=Julian.date,ymin=mean_sugar-1.96*se_sugar, ymax=mean_sugar+1.96*se_sugar, fill= "Sugar"), 
                data=nu,  alpha=0.4) +
    scale_color_manual(values = c("Total NSC" = "purple","Starch" = "red","Sugar" = "blue")) +
    scale_fill_manual(values = c("Total NSC" = "purple","Starch" = "red","Sugar" = "blue")) +
    guides(fill = guide_legend(reverse = TRUE)) +
    theme(legend.title = element_blank(),
          legend.position = c(0.82, 0.15),        # (L --> R, B --> T))
          legend.background = element_blank(),    # remove grey rectangle from behind color dots on legend
          legend.margin = margin(0, 8, 6, 8),   # change padding around legend (top, right, bottom, left)
          legend.key.width = unit(2.5, "line"), 
          legend.text=element_text(size=16),
          strip.text =element_text(size=18),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 13), #, hjust=0.8
          axis.title.y = element_text(size=18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))  +
    coord_cartesian(ylim=c(0,0.30)) +
    scale_y_continuous(limits=c(-0.5,0.30),
                       labels=c("0", "10","20","30")) +
    scale_x_continuous(breaks=c(25,55,111,156,233,293),
                       labels=c("Jan. \nD","  Mar.\nSF","Apr.\nLO","Jun.\nSE","Aug.\nBS","Oct.\nLF")) +
    facet_wrap(vars(tissue), nrow=3, as.table = T) 
  
  ggsave("Fig_2.jpeg", width = 9.5, height = 8, units="in", dpi=330)
  
  
####### Fig. 3.  Whole-tree biomass                        #####
  
  rm(bioT)
  bioT <- biomass[which(biomass$tissue=="Total biomass"),]
  bioT
  c <- read.csv("C_2m_out.csv")
  c <- c[which(c$tissue=="Total biomass"),]
  
  fitlme <- glmer(biomass ~ DBH  * tissue + Ht + (1|tree) , data=biomass, family=Gamma(link="log"))
  newdat.lme = data.frame(DBH = bioT$DBH,
                          tissue=bioT$tissue,
                          Ht = median(bioT$Ht))   # 
  head(newdat.lme)
  bioT$predlme = predict(fitlme, newdata = newdat.lme, re.form=NA, type="response")
  bioT
  pred.CI <- easyPredCI(fitlme, newdat.lme)
  bioT <- cbind(bioT,pred.CI)
  
  newdat.lme
  
  newdat.lme <- cbind(newdat.lme, pred.CI)
  newdat.lme <- newdat.lme[which(newdat.lme$tissue=="Total biomass"),]
  newdat.lme <- newdat.lme[which(newdat.lme$DBH>2.5),]
  newdat.lme <- newdat.lme[which(newdat.lme$DB<16),]
  min(newdat.lme$DBH)
  
  ggplot() + theme_bw() + 
    geom_line(data = bioT, aes(x=DBH,y = predlme), size = 1, col="grey") +
    geom_ribbon(data = bioT, aes(x=DBH,ymin = resp_lwr, ymax=resp_upr), size = 1, fill="grey", alpha=0.3) +
    geom_point(data=bioT, aes(x=DBH, y=biomass), size=3, shape=16) +
    labs(y="Oven-dry biomass (kg)", 
         x = "DBH") +
    theme_classic() + 
    theme(title=element_text(size=17),
          legend.title = element_blank(),
          legend.text = element_text(size=16),
          legend.background=element_blank(),     # remove grey rectangle from behind color dots on legend
          legend.margin = margin(0, 10, 6, 10),     # change padding around legend (top, right, bottom, left)
          legend.position="none",
          axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size=16),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size=16),
          strip.text = element_text(size=16),
          strip.background = element_blank()) 
  
  ggsave("fig.3.jpeg", width = 4, height = 3.5, units="in", dpi=330)
  
####### Fig. 4.  Organ biomass                             #####
  
  rm(bioT)
  bioT <- biomass
  bioT
  c <- read.csv("C_2m_out.csv")
  
  fitlme <- glmer(biomass ~ DBH  * tissue + Ht + (1|tree) , data=biomass, family=Gamma(link="log"))
  newdat.lme = data.frame(DBH = bioT$DBH,
                          tissue=bioT$tissue,
                          Ht = median(bioT$Ht))   # 
  head(newdat.lme)
  bioT$predlme = predict(fitlme, newdata = newdat.lme, re.form=NA, type="response")
  bioT
  pred.CI <- easyPredCI(fitlme, newdat.lme)
  bioT <- cbind(bioT,pred.CI)
  
  bbT <- spread(bioT[,c(1:6,8)], tissue, predlme)
  bbT 
  
  bbT$t <- bbT$'Fine roots' + bbT$'Coarse roots' + bbT$'Root collar' + bbT$Xylem + bbT$'Inner bark' + bbT$Branches + bbT$Twigs + bbT$Leaves
  bbT$fr_p <- bbT$'Fine roots'/bbT$t
  bbT$cr_p <- (bbT$'Fine roots' + bbT$'Coarse roots')/bbT$t
  bbT$rc_p <- (bbT$'Fine roots' + bbT$'Coarse roots' + bbT$'Root collar')/bbT$t
  bbT$x_p <- (bbT$'Fine roots' + bbT$'Coarse roots' + bbT$'Root collar' + bbT$Xylem)/bbT$t
  bbT$p_p <- (bbT$'Fine roots' + bbT$'Coarse roots' + bbT$'Root collar' + bbT$Xylem + bbT$'Inner bark')/bbT$t
  bbT$b_p <- (bbT$'Fine roots' + bbT$'Coarse roots' + bbT$'Root collar' + bbT$Xylem + bbT$'Inner bark' + bbT$Branches)/bbT$t
  bbT$t_p <- (bbT$'Fine roots' + bbT$'Coarse roots' + bbT$'Root collar' + bbT$Xylem + bbT$'Inner bark' + bbT$Branches + bbT$Twigs)/bbT$t
  bbT$l_p <- (bbT$'Fine roots' + bbT$'Coarse roots' + bbT$'Root collar' + bbT$Xylem + bbT$'Inner bark' + bbT$Branches + bbT$Twigs + bbT$Leaves)/bbT$t
  bbT
  
  ggplot(data=bbT) + theme_bw() + 
    geom_ribbon(aes(x=DBH,ymin = 0, ymax=fr_p, fill="Fine roots")) +
    geom_ribbon(aes(x=DBH,ymin = fr_p, ymax=cr_p, fill="Coarse roots")) +
    geom_ribbon(aes(x=DBH,ymin = cr_p, ymax=rc_p, fill="Root collar")) +
    geom_ribbon(aes(x=DBH,ymin = rc_p, ymax=x_p, fill="Xylem")) +
    geom_ribbon(aes(x=DBH,ymin = x_p, ymax=p_p,  fill="Phloem + bark")) +
    geom_ribbon(aes(x=DBH,ymin = p_p, ymax=b_p,  fill="Branches")) +
    geom_ribbon(aes(x=DBH,ymin = b_p, ymax=t_p,  fill="Twigs")) +
    geom_ribbon(aes(x=DBH,ymin = t_p, ymax=l_p, fill="Leaves")) +
    scale_fill_manual(values = c("Leaves" = "green4","Twigs" = "red","Branches" = "blue", "Phloem + bark" = "tan4",
                                 "Xylem" = "orange","Root collar" = "orchid2","Coarse roots" = "indianred4","Fine roots" = "black"),
                      breaks = c("Leaves","Twigs","Branches","Phloem + bark","Xylem","Root collar","Coarse roots","Fine roots")) +
    labs(y="Proportion of total biomass", 
         x = "DBH") +
    theme_classic() + 
    theme(title=element_text(size=17),
          legend.title = element_blank(),
          legend.text = element_text(size=14),
          legend.background=element_blank(),     # remove grey rectangle from behind color dots on legend
          legend.margin = margin(0, 10, 6, 10),     # change padding around legend (top, right, bottom, left)
          #legend.position="none",
          axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size=16),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size=16),
          strip.text = element_text(size=16),
          strip.background = element_blank()) 
  
  ggsave("Fig.4.jpeg", width = 6.5, height = 4, units="in", dpi=330)
  
####### Fig. 5.  Median-tree pools                         ####
  
  ggplot(data=med_wt) + theme_minimal() + 
    geom_point(aes(x=Julian.date, y=NSC_pool, col="Total NSC")) +
    geom_point(aes(x=Julian.date, y=starch_pool, col="Starch")) +
    geom_point(aes(x=Julian.date, y=sugar_pool, col="Sugar")) +
    geom_line(aes(x=Julian.date, y=NSC_pool, col= "Total NSC"), size=1.2) +
    geom_line(aes(x=Julian.date, y=starch_pool, col= "Starch"), size=1.2) +
    geom_line(aes(x=Julian.date, y=sugar_pool, col= "Sugar"), size=1.2) +
    geom_ribbon(aes(x=Julian.date, ymin=NSC_pool-1.96*NSC_pool_se, ymax=NSC_pool+1.96*NSC_pool_se, fill="Total NSC"), 
                alpha=0.3, show.legend=F) + 
    geom_ribbon(aes(x=Julian.date, ymin=starch_pool-1.96*starch_pool_se, ymax=starch_pool+1.96*starch_pool_se, fill="Starch"), 
                alpha=0.3, show.legend=F) +
    geom_ribbon(aes(x=Julian.date, ymin=sugar_pool-1.96*sugar_pool_se, ymax=sugar_pool+1.96*sugar_pool_se, fill="Sugar"), 
                alpha=0.4, show.legend=F) +
    scale_color_manual(values = c("Total NSC" = "purple","Starch" = "red","Sugar" = "blue")) +
    scale_fill_manual(values = c("Total NSC" = "purple","Starch" = "red","Sugar" = "blue")) +
    ylab("Carbohydrate pool (kg)") + 
    guides(color = guide_legend(reverse = TRUE,override.aes=list(fill=NA))) + 
    scale_x_continuous(breaks=c(25,55,111,156,233,293),
                       labels=c("Jan. \nD","  Mar.\nSF","Apr.\nLO","Jun.\nSE","Aug.\nBS","Oct.\nLF")) +
    theme(legend.position = "bottom",                        
          legend.background = element_blank(),             # remove grey rectangle from behind color dots on legend
          legend.margin = margin(0, 8, 6, 8),              # change padding around legend (top, right, bottom, left)
          legend.key.width = unit(2.8, "line"), 
          legend.text=element_text(size=18),
          legend.title=element_blank(),
          strip.text =element_text(size=16, face="bold"),
          axis.title.y = element_text(size=18, face="bold"),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 13),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) +
    facet_wrap(vars(tissue), scales="free_y", nrow=3) 
  
   ggsave("Fig_5.jpeg", width = 9.5, height = 8, units="in", dpi=330)
   
####### Fig. 6.  Delta-NSC                                 ####

   t  # from Delta NSC calculation above (Steps 4 & 5)
   t$tissue <- factor(t$tissue, levels = c("Twigs", "Branches", "Inner bark", "Xylem","Root collar", "Coarse roots", "Fine roots"))
   t
   # make seperate starch, sugar, NSC plots then stitch together using arrangeGrob()
   
   t_St <- ggplot(t,aes(x=tissue, y=starch_flux))  + theme_classic() +
     geom_bar(stat = "identity", fill="red", alpha=0.7) + 
     geom_errorbar(aes(x=tissue, ymin=starch_flux-starch_flux_se*1.96, ymax=starch_flux+starch_flux_se*1.96), size=0.9, width=0.7) +
     labs(y=expression(""*Delta*" starch (kg)")) +
     xlab("") +                                     
     theme(axis.title.y=element_text(size=15),
           axis.text.y=element_text(size=14),
           axis.title.x=element_text(size=15),
           axis.text.x=element_text(size=14),
           axis.ticks = element_blank(),
           plot.margin = unit(c(1,1,0,1), "lines")) +
     scale_x_discrete(labels=c("Twigs","Branches", "Inner\nbark","Xylem*",
                               "Root\ncollar","Coarse\nroots", "Fine\nroots", "Whole*\ntree")) +
     scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits=c(0,1.55))             
   t_St
   
   t_Su <- ggplot(t,aes(x=tissue, y=sugar_flux))  + theme_classic() +
     geom_bar(stat = "identity", fill="blue", alpha=0.7) + 
     geom_errorbar(aes(x=tissue, ymin=ifelse(sugar_flux-sugar_flux_se*1.96 < 0, 0, sugar_flux-sugar_flux_se*1.96), ymax=sugar_flux+sugar_flux_se*1.96), size=0.9, width=0.7) +
     labs(y=expression(""*Delta*" sugar (kg)")) +
     xlab("") +                                                  
     theme(axis.title.y=element_text(size=15),
           axis.text.y=element_text(size=14),
           axis.title.x= element_blank(),
           axis.text.x= element_blank(),
           axis.ticks = element_blank(),
           plot.margin = unit(c(1,1,0,1), "lines")) +
     scale_x_discrete(labels=c("Twigs","Branches", "Inner\nbark","Xylem*",
                               "Root\ncollar","Coarse\nroots", "Fine\nroots", "Whole*\ntree")) +
     scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits=c(0,1.55)) 
   t_Su
   
   t_NSC <- ggplot(t,aes(x=tissue, y=NSC_flux))  + theme_classic() +
     geom_bar(stat = "identity", fill="purple", alpha=0.7) + 
     geom_errorbar(aes(x=tissue, ymin=NSC_flux-NSC_flux_se*1.96, ymax=NSC_flux+NSC_flux_se*1.96), size=0.9, width=0.7) +
     labs(y=expression(""*Delta*" NSC (kg)")) +
     xlab("") +                                     
     theme(axis.title.y=element_text(size=15),
           axis.text.y=element_text(size=14),
           axis.title.x= element_blank(),
           axis.text.x= element_blank(),
           axis.ticks = element_blank(),
           plot.margin = unit(c(1,1,0,1), "lines")) +
     scale_x_discrete(labels=c("Twigs","Branches", "Inner\nbark","Xylem*",
                               "Root\ncollar","Coarse\nroots", "Fine\nroots", "Whole*\ntree")) +
     scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits=c(0,1.55))
   t_NSC
   
   t_bp <- ggarrange(t_NSC, t_Su,t_St, nrow=3, heights=c(2.4,2.4,3.3)) # check what it looks like
   t_bp
   t_bp <- arrangeGrob(t_NSC, t_Su,t_St, nrow=3, heights=c(2.4,2.4,3.3)) # can be exported using ggsave
   ggsave("fig.6.jpeg", t_bp, width = 7, height = 6, units="in", dpi=330)
   
####### Fig. S1. DBH x concentration                       ####
  
  # make seperate figures for starch, sugar, & total NSC, then stich them together with ggarange()
  
  st <- ggplot(data=NSC) + theme_classic() + 
    ylab("Starch (%)") +
    geom_point(aes(x=DBH, y=starch_fraction),
               alpha= 0.6, size=1.3, shape=16, position="jitter") +
    geom_smooth(aes(x=DBH, y=starch_fraction),
                method="lm",se=T, size=1, fill="blue", alpha= 0.3) +
    scale_y_continuous(limits=c(0,0.3), labels=c("0", "10","20","30"))+
    theme(legend.background = element_blank(), # remove grey rectangle from behind color dots on legend  
          legend.margin = margin(0, 0,0,0),  # change padding around legend (top, right, bottom, left)
          legend.key.width = unit(2.5, "line"),
          legend.text=element_text(size=16),
          legend.title=element_blank(),
          strip.text =element_text(size=14),
          axis.text.x = element_blank(),
          axis.title.x =element_blank(),
          axis.text.y =element_text(size=14),
          axis.title.y =element_text(size=16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) +
    facet_wrap(~collection3,nrow=1)
  st
  
  su <- ggplot(data=NSC) + theme_classic() + 
    ylab("Sugar (%)") +
    geom_point(aes(x=DBH, y=sugar_fraction),
               alpha= 0.6, size=1.3, shape=16, position="jitter") +
    geom_smooth(aes(x=DBH, y=sugar_fraction),
                method="lm",se=T, size=1, fill="blue", alpha= 0.3) +
    scale_y_continuous(limits=c(0,0.3), labels=c("0", "10","20","30"))+
    theme(legend.background = element_blank(), # remove grey rectangle from behind color dots on legend  
          legend.margin = margin(0, 8, 6, 8),  # change padding around legend (top, right, bottom, left)
          legend.key.width = unit(2.5, "line"),
          legend.text=element_text(size=18),
          legend.title=element_blank(),
          strip.text =element_blank(),
          axis.text.x = element_blank(),
          axis.title.x =element_blank(),
          axis.text.y =element_text(size=14),
          axis.title.y =element_text(size=16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) +
    facet_wrap(~collection3,nrow=1)
  su
  
  NN <- ggplot(data=NSC) + theme_classic() + 
    ylab("NSC (%)") +
    geom_point(aes(x=DBH, y=tot_NSC_fraction),
               alpha= 0.6, size=1.3, shape=16, position="jitter") +
    geom_smooth(aes(x=DBH, y=tot_NSC_fraction),
                method="lm",se=T, size=1, fill="blue", alpha= 0.3) +
    scale_y_continuous(limits=c(0,0.3), labels=c("0", "10","20","30"))+
    theme(legend.background = element_blank(), # remove grey rectangle from behind color dots on legend  
          legend.margin = margin(0, 8, 6, 8),  # change padding around legend (top, right, bottom, left)
          legend.key.width = unit(2.5, "line"),
          legend.text=element_text(size=18),
          legend.title=element_blank(),
          strip.text =element_blank(),
          axis.text.x = element_text(angle = 0, size = 12, color="black"),
          axis.title.x=element_text(size=16),
          axis.text.y =element_text(size=14),
          axis.title.y =element_text(size=16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) +
    facet_wrap(~collection3,nrow=1)
  NN
  
  grid.arrange(arrangeGrob(st + theme(legend.position="none"),
                           su + theme(legend.position="none"),
                           NN + theme(legend.position="none"),
                           nrow=3, heights=c(2.2,2,2.2)))  # ,mylegend, ncol=2,heights=c(6, 1),  widths = c(4, 0.7)
  
  g <- arrangeGrob(st, su, NN, nrow=3, heights=c(2.2,2,2.4))
  ggsave("g.jpeg", g, width = 9.7, height = 6, units="in", dpi=330)
  
####### Fig. S2. Median biomass                            ######
  
  med_wt_bio <- med_wt[which(med_wt$collection=="2-shoot expansion"),]
  med_wt_bio$tissue2 <- factor(med_wt_bio$tissue, levels =  c("Whole tree","Leaves", "Twigs", "Branches", "Inner bark", "Xylem","Root collar", "Coarse roots", "Fine roots"))
  med_wt_bio
  
  ggplot(med_wt_bio,aes(x=tissue2, y=tissue_biomass))  + theme_classic() +
    geom_bar(stat = "identity", alpha=0.7) + 
    geom_errorbar(aes(x=tissue2, ymin=tissue_biomass-tissue_biomass_se*1.96, ymax=tissue_biomass+tissue_biomass_se*1.96), size=0.9, width=0.7) +
    labs(y=expression("Biomass (kg)")) +
    xlab("") +                                      # coord_trans(y = "log10") +                # ylim(0, 3.0) +  used previously
    theme(axis.title.y=element_text(size=16),
          axis.text.y=element_text(size=14),
          axis.text.x=element_text(size=14),
          axis.ticks = element_blank()) +
    scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                       breaks=c(0, 1, 10, 20,30,40,50)) +
    scale_x_discrete(labels=c("Whole\ntree","Leaves","Twigs","Branches", "Inner\nbark","Xylem",
                              "Root\ncollar","Coarse\nroots", "Fine\nroots")) 
  
  ggsave("fig.S2.jpeg", width = 8, height = 6, units="in", dpi=330)
  
