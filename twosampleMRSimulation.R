sink('fixed_delta.txt')
#Fn===========
library(ggpubr)

onesample_twostep = function(X,Y,G){
  #SNP-X associations by regressing X on all G
  BXG_temp = summary(lm(X[,1] ~ .,data = data.frame(G)))$coefficients[-1,1]
  BXG=seBXG=pv_df=(matrix(nrow = length(BXG_temp),ncol = ncol(X)))
  for (k in 1:ncol(X)){
    BXG[,k] = summary(lm(X[,k] ~ .,data = data.frame(G)))$coefficients[-1,1]
    seBXG[,k] =summary(lm(X[,k] ~ .,data = data.frame(G)))$coefficients[-1,2]
    pv_df[,k] =summary(lm(X[,k] ~ .,data = data.frame(G)))$coefficients[-1,4]
  }
  #print('Converted to summary statistics')
  colnames(BXG)=paste('BXG',1:ncol(X),sep = '_')
  seBetaYG=summary(lm(Y~.,data = data.frame(G)))$coefficients[-1,2]
  
  #Fit and induce Collider bias
  YXGdata = data.frame(Y,X,G)
  #
  #Logistic Regression _that conditions on X_. Here, the collider bias is induced
  #
  FIT4             = summary(lm(YXGdata)) #estimated associations
  alphahatstar     = FIT4$coef[-seq(1,ncol(X)+1), 1]  #obtain the estimates for the genetic variant coefficient,
  se.alphahatstar  = FIT4$coef[-seq(1,ncol(X)+1), 2]
  betastar         = FIT4$coef[seq(2,ncol(X)+1),1]
  betastar_SE      = FIT4$coef[seq(2,ncol(X)+1),2]
  # print('Collider bias induced')
  
  ##IVW but alphahatstar instead of the outcome
  
  # #Grapple: newer version data input
  # #library(GRAPPLE)
  # data_grapl = data.frame(SNP=1:nrow(BXG),BXG,seBXG,alphahatstar,se.alphahatstar)
  # colnames(data_grapl)[2:(ncol(BXG)+1)]=paste('gamma_exp',1:ncol(BXG),sep='')
  # colnames(data_grapl)[(ncol(BXG)+2):(2*ncol(BXG)+1)]=paste('se_exp',1:ncol(BXG),sep='')
  # colnames(data_grapl)[(2*ncol(BXG)+2):ncol(data_grapl)] = c('gamma_out1','se_out1')
  # 
  # grapl1=grappleRobustEst(data = data_grapl)
  # bGRAPPLE=grapl1$beta.hat
  # MRE4      = betastar + bGRAPPLE
  # SE_MRE4 = sqrt(betastar_SE^2 + diag(grapl1$beta.var))
  # 
  # #You can fit any two-sample MR method (median, mode, lasso, etc.) with this
  # #Here I'm just using GRAPPLE
  # #
  # #The general formula would be ahat* ~ X1 + X2+ ... 
  # #And then correct the estimates by beta_X1_c = beta_X1 + betastar
  return(list(BXG=BXG,
              seBXG=seBXG,
              beta_star = betastar, 
              beta_star_SE = betastar_SE,
              alphahatstar=alphahatstar,
              se.alphahatstar=se.alphahatstar,
              pv_df=pv_df))
}
#Sim===========
library(mr.raps)
N_SIM = 400; bias_df = empSE_df = data.frame(matrix(nrow = N_SIM,ncol = 2))
RMSE_df = data.frame(matrix(nrow = N_SIM,ncol = 2))
N_vec = round(seq(2e2,5e3,length.out = 8))
Delta_S_X_vec = 2

RMSE_fin = bias_fin = empSE_fin = data.frame(matrix(nrow = length(N_vec)*length(Delta_S_X_vec),
                                                    ncol = 3))
f_strat_r = f_strat_temp=f_strat_sd = fish_z_sd=NULL
fis_z_temp = fis_z =fis_z_sd= NULL
rownames(RMSE_fin)=rownames(bias_fin)=rownames(empSE_fin)=
  apply(expand.grid(paste(N_vec,'SNPs'), paste('Delta',Delta_S_X_vec)), 
        1, paste, collapse="_")

colnames(RMSE_fin)=colnames(bias_fin)=colnames(empSE_fin)=
  c('RAPS_Strat','RAPS_unstrat','fish_z')

outlier_list=list()

for(di in 1:length(Delta_S_X_vec)){
  for (np in 1:length(N_vec)){
    n_temp = N_vec[np]
    for (nsim in 1:N_SIM){
      #sex-stratified
      n_ind = n_temp; p = 35
      S = rbinom(n = n_ind,size = 1,prob = 0.5)
      G = matrix(rbinom(n = n_ind*p,size = 1,prob = 0.2),
                 nrow = n_ind,ncol = p)
      U = rnorm(n = n_ind,0,3); epsilon_X = rnorm(n = n_ind,mean = 0,sd = 3)
      
      alpha_G_Y = matrix(rnorm(p,0.02,1), nrow = p, ncol = 1)
      
      gamma_SNP_X = matrix(rnorm(p,0.3,0.1), nrow = p, ncol = 1) + alpha_G_Y
      Delta_S_X = matrix(rnorm(p,Delta_S_X_vec[di],0.01), 
                         nrow = p, ncol = 1)
      
      X = scale(G %*% gamma_SNP_X) + U + epsilon_X
      X[which(S==1)] = X[which(S==1)] + S[which(S==1)] * (G[which(S==1),] %*% Delta_S_X)
      #X=scale(X)
      beta_X_Y = matrix(0.3,nrow = 1,ncol = 1)
      
      beta_2 = 1
      epsilon_Y = rnorm(n = n_ind,mean = 0,sd = 8)
      
      #SECOND INDEPENDENT SAMPLE FOR OUTCOME DATA
      S_2 = rbinom(n = n_ind,size = 1,prob = 0.5)
      G_2 = matrix(rbinom(n = n_ind*p,size = 1,prob = 0.2),
                 nrow = n_ind,ncol = p)
      U_2 = rnorm(n = n_ind,0,5); epsilon_X_2 = rnorm(n = n_ind,mean = 0,sd = 5)
      
      #Gamma and Delta same as above
      #gamma_SNP_X = matrix(0.4, nrow = p, ncol = 1)
      #Delta_S_X = matrix(Delta_S_X_vec[di], nrow = p, ncol = 1)
      
      #Using same S and G as in the exposures as the S_2 and G_2 are biased
      #Those should be different , G_2, S_2
      X_2 = G_2 %*% gamma_SNP_X + S_2 * (G_2 %*% Delta_S_X) + U_2 + epsilon_X_2
      
      #beta_X_Y = matrix(0.3,nrow = 1,ncol = 1)
      #alpha_G_Y = matrix(0.3, nrow = p, ncol = 1)
      #beta_2 = 0.2
      epsilon_Y_2 = rnorm(n = n_ind,mean = 0,sd = 5)
      
      Y = X_2 %*% beta_X_Y + S_2*beta_2 + G_2 %*% alpha_G_Y + U_2 + epsilon_Y_2
     # Y=scale(Y)
      
      coll_sam_s0=onesample_twostep(X = data.frame(X[which(S==0),]),
                                    G = data.frame(G[which(S==0),]),
                                    Y = as.numeric(Y[which(S==0),]))
      coll_sam_s1=onesample_twostep(X = data.frame(X[which(S==1),]),
                                    G = data.frame(G[which(S==1),]),
                                    Y = as.numeric(Y[which(S==1),]))
      #plot(coll_sam_s0$BXG,coll_sam_s1$BXG)
      
      #G-X
      #1. Collider Corrction
      gx_s1=coll_sam_s1$BXG; gx_s1_SE=coll_sam_s1$seBXG
      gx_s0=coll_sam_s0$BXG; gx_s0_SE=coll_sam_s0$seBXG
      
      #2. Summary Estimates, no correction
      gy_s1=summary(lm(Y[S==1,] ~ G_2[S==1,]))$coef[-1,1]
      gy_s1_SE=summary(lm(X[S==1,] ~ G[S==1,]))$coef[-1,2]
      
      gy_s0=summary(lm(Y[S==0,] ~ G_2[S==0,]))$coef[-1,1]
      gy_s0_SE=summary(lm(Y[S==0,] ~ G_2[S==0,]))$coef[-1,2]
      
      
      #M1: RAPS+stratification
      data_rap = data.frame(beta.exposure = gx_s1-gx_s0, 
                            beta.outcome = gy_s1-gy_s0,
                            se.outcome = sqrt(gy_s1_SE^2+gy_s0_SE^2),
                            se.exposure = sqrt(gx_s1_SE^2+gx_s0_SE^2))
      colnames(data_rap) = c('beta.exposure','beta.outcome',
                             'se.exposure','se.outcome')
      #raps_strat=mr.raps(data_rap,diagnostics = F,
      #                  over.dispersion = F,
      #                   loss.function = 'huber')   
      raps_strat=mr.raps(b_exp=data_rap$beta.exposure,b_out=data_rap$beta.outcome,
                         se_exp=data_rap$se.exposure,se_out=data_rap$se.outcome,
                         loss.function = 'huber')
      
      #M2: RAPS+no stratification
      coll_unstrat = onesample_twostep(X = data.frame(X),
                                       G = data.frame(G),
                                       Y = as.numeric(Y))
      gy_unstr=summary(lm(Y ~ G_2))$coef[-1,1]
      gy_SE_unstr=summary(lm(Y ~ G_2))$coef[-1,2]
      
      
      data_rap_unstr = data.frame(beta.exposure = coll_unstrat$BXG, 
                                  beta.outcome = gy_unstr,
                                  se.exposure = coll_unstrat$seBXG,
                                  se.outcome = gy_SE_unstr)
      colnames(data_rap_unstr) = c('beta.exposure','beta.outcome',
                                   'se.exposure','se.outcome')
      #raps_unstrat=mr.raps(data_rap_unstr,diagnostics = F,
      #                     over.dispersion = F,
      #                     loss.function = 'huber')
      raps_unstrat=mr.raps(b_exp=data_rap_unstr$beta.exposure,b_out=data_rap_unstr$beta.outcome,
                           se_exp=data_rap_unstr$se.exposure,se_out=data_rap_unstr$se.outcome,
                           loss.function = 'huber')
      
      #=======================================================================================================
      #=======================================================================================================
      #=======================================================================================================
      bias_df[nsim,1]=raps_strat$beta.hat - beta_X_Y
      bias_df[nsim,2]=raps_unstrat$beta.hat - beta_X_Y
      
      empSE_df[nsim,1]=raps_strat$beta.se
      empSE_df[nsim,2]=raps_unstrat$beta.se
      
      RMSE_df[nsim,1]=(raps_strat$beta.hat - beta_X_Y)^2
      RMSE_df[nsim,2]=(raps_unstrat$beta.hat - beta_X_Y)^2
      
      #fisher's z
      fis_z_temp[nsim]=mean((gx_s1 - gx_s0)^2 / (gx_s1_SE^2+gx_s0_SE^2))
      
      
      #fstat
      f_strat_temp[nsim]=mean((gx_s1 / gx_s1_SE)^2)
    }
    print(np/length(N_vec)) ###############
    outlier_thresh=100
    b_outl=unlist(apply(bias_df,2,function(x) which(abs(x)>outlier_thresh)))
    r_outl=unlist(apply(RMSE_df,2,function(x) which(abs(x)>outlier_thresh)))
    e_outl=unlist(apply(empSE_df,2,function(x) which(abs(x)>outlier_thresh)))
    
    outlier_list[[np]]=unique(c(b_outl,r_outl,e_outl))
    if (length(outlier_list[[np]]) !=0){
    bias_df=bias_df[-outlier_list[[np]],]
    RMSE_df=RMSE_df[-outlier_list[[np]],]
    empSE_df=empSE_df[-outlier_list[[np]],]
    }
    
    fis_z[np]=mean(fis_z_temp);fis_z_sd[np]= sd(fis_z_temp)
    f_strat_r[np] = mean(f_strat_temp); f_strat_sd[np] = sd(f_strat_temp)
    #print(paste(paste(np,'np'), paste('np + np*(di-1)',np + (np-1)*(di-1))))
    RMSE_fin[np + length(N_vec ) * (di-1),1:2] = sqrt(apply(RMSE_df,2,mean)/N_SIM)
    bias_fin[np + length(N_vec ) * (di-1),1:2] = apply(bias_df,2,mean)
    empSE_fin[np + length(N_vec ) * (di-1),1:2] = apply(empSE_df,2,mean)
    RMSE_fin[np + length(N_vec ) * (di-1),3] = bias_fin[np + length(N_vec ) * (di-1),3] = fis_z[np]}
  }

outl_excl = unique(c(which(complete.cases(bias_fin$RAPS_Strat)), 
                     which(complete.cases(RMSE_fin$RAPS_unstrat))))
bias_fin = bias_fin[outl_excl,]
RMSE_fin = RMSE_fin[outl_excl,]
empSE_fin = empSE_fin[outl_excl,]

#Plot----------
library(ggplot2);library(stringr);library(tidyverse)
df1 = data.frame(SNP = as.numeric(gsub("([0-9]+).*$", "\\1", rownames(bias_fin))),
                 delta=unlist(str_split(string = rownames(bias_fin),pattern = ' '))[seq(3,length(unlist(str_split(string = rownames(RMSE_fin),pattern = ' '))),3)],
                 RAPS_strat = bias_fin$RAPS_Strat,
                 RAPS_unstrat = bias_fin$RAPS_unstrat ,
                 fish_z = RMSE_fin$fish_z,
                 fish_z_sd=fis_z_sd[outl_excl],
                 fstat = f_strat_r[outl_excl],
                 fstat_sd=f_strat_sd[outl_excl])
df1$SNP[which(df1$SNP == 1)] = 1e5
df1_t=df1 %>% 
  pivot_longer(c(RAPS_strat,RAPS_unstrat ), names_to = "method", values_to = "bias")
#df1_t

s1=ggplot(df1, aes(x=SNP)) +geom_line(aes(y = fish_z), color = "darkred") + geom_point(aes(y = fstat)) + 
  geom_errorbar(aes(ymin=fstat-fstat_sd, ymax=fstat+fstat_sd), width=.2,position=position_dodge(0.05))+
  geom_line(aes(y = fstat), color="steelblue", linetype="twodash") + 
  theme_bw()+ylab('F-statistic,Fishers z')+ geom_point(aes(y = fish_z)) + 
  geom_errorbar(aes(ymin=fish_z-fish_z_sd, ymax=fish_z+fish_z_sd), width=.2,position=position_dodge(0.05))
bias_pl=df1_t %>% ggplot( aes(x=SNP, y=bias, group=method, color=method)) +
  geom_line() +  theme_bw()+ theme(legend.position = c(0.8, max(df1_t$bias)/9))+
  ylab('Bias') + geom_hline(yintercept = 0) #+ ylim(0,2)

library(cowplot)
library(ggplot2);library(stringr);library(tidyverse)
################################################################################################

df1 = data.frame(SNP = as.numeric(gsub("([0-9]+).*$", "\\1", rownames(empSE_fin))),
                 delta=rep(Delta_S_X_vec,length(N_vec)),
                 RAPS_strat = empSE_fin$RAPS_Strat,
                 RAPS_unstrat = empSE_fin$RAPS_unstrat ,
                 fish_z = RMSE_fin$fish_z,
                 fish_z_sd=fis_z_sd[outl_excl],
                 fstat = f_strat_r[outl_excl],
                 fstat_sd=f_strat_sd[outl_excl])
df1$SNP[which(df1$SNP == 1)] = 1e5
df1_t=df1 %>% 
  pivot_longer(c(RAPS_strat,RAPS_unstrat ), names_to = "method", values_to = "bias")
#df1_t

s1=ggplot(df1, aes(x=SNP)) +geom_line(aes(y = fish_z), color = "darkred") + geom_point(aes(y = fstat)) + 
  geom_errorbar(aes(ymin=fstat-fstat_sd, ymax=fstat+fstat_sd), width=.2,position=position_dodge(0.05))+
  geom_line(aes(y = fstat), color="steelblue", linetype="twodash") + 
  theme_bw()+ylab('F-statistic,Fishers z')+ geom_point(aes(y = fish_z)) + 
  geom_errorbar(aes(ymin=fish_z-fish_z_sd, ymax=fish_z+fish_z_sd), width=.2,position=position_dodge(0.05))
empSE_pl=df1_t %>% ggplot( aes(x=SNP, y=bias, group=method, color=method)) +
  geom_line() +  theme_bw()+ theme(legend.position = c(0.8, max(df1_t$bias)/9))+
  ylab('EmpSE') + geom_hline(yintercept = 0) #+ ylim(0,2)

################################################################################################

#stack df1 and df2 with rbind(), add an extra column 'type_of_plot' with bias, RMSE
################################################################################################

df1 = data.frame(SNP = as.numeric(gsub("([0-9]+).*$", "\\1", rownames(RMSE_fin))),
                 delta=unlist(str_split(string = rownames(RMSE_fin),pattern = ' '))[seq(3,length(unlist(str_split(string = rownames(RMSE_fin),pattern = ' '))),3)],
                 RAPS_strat = RMSE_fin$RAPS_Strat,
                 RAPS_unstrat = RMSE_fin$RAPS_unstrat ,
                 fish_z = RMSE_fin$fish_z,
                 fish_z_sd=fis_z_sd[outl_excl],
                 fstat = f_strat_r[outl_excl],
                 fstat_sd=f_strat_sd[outl_excl])
df1$SNP[which(df1$SNP == 1)] = 1e5
df1_t=df1 %>% 
  pivot_longer(c(RAPS_strat,RAPS_unstrat ), names_to = "method", values_to = "bias")
#df1_t

s1=ggplot(df1, aes(x=SNP)) +geom_line(aes(y = fish_z), color = "darkred") + geom_point(aes(y = fstat)) + 
  geom_errorbar(aes(ymin=fstat-fstat_sd, ymax=fstat+fstat_sd), width=.2,position=position_dodge(0.05))+
  geom_line(aes(y = fstat), color="steelblue", linetype="twodash") + 
  theme_bw()+ylab('F-statistic,Fishers z')+ geom_point(aes(y = fish_z)) + 
  geom_errorbar(aes(ymin=fish_z-fish_z_sd, ymax=fish_z+fish_z_sd), width=.2,position=position_dodge(0.05))
rmse_pl=df1_t %>% ggplot( aes(x=SNP, y=bias, group=method, color=method)) +
  geom_line() +  theme_bw()+ theme(legend.position = c(0.8, max(df1_t$bias)/9))+
  ylab('RMSE') + geom_hline(yintercept = 0) #+ ylim(0,2)

################################################################################################

#facet_wrap to combine bias & rmse plots

# 1. join together, specify bias & rmse, 2. facet wrap to get the bias rmse, 3. facet grid

bias_pl2=plot_grid(s1, bias_pl,rmse_pl,empSE_pl,ncol=1, align = "h", axis = "bt")
#Save---------
#bias_pl2
ggsave(filename = 'bias_pl_delta.pdf',plot = bias_pl2)

sink()
