
library(plyr)
library(tidyverse)
library(nimble)
library(MCMCvis)
library(abind)
library(bayesplot)

source("R/misc_functions.r")

##---- load data ----

load("Data/deerCE.Rdata")

##---- construct data for Nimble ----

K<- dim(y)[1]# number of sites
T <- dim(y)[3] # number of primary periods
J<- dim(y)[2] # secondary periods per primary period

n<- apply(y, c(1,3), sum)  # marginal total removals per site and primary period


constants<- list(K=K, J=J, T=T) 

data <- list(y=y, n=n, eff=eff, garea=garea, elev=elevation, season=season)


##---- Nimble model code ----


code<- nimbleCode({
  for (i in 1:K) {
    for(t in 1:T){
      pi[i,1,t]<- p[i,1,t]
      for(j in 2:J){
        pi[i,j,t]<- prod(1-p[i,1:(j-1),t]) * p[i,j,t]
      }
      pcap[i,t]<- sum(pi[i,1:J,t])
      for(j in 1:J){
        cloglog(p[i,j,t])<- alpha + log(eff[i,j,t]+1e-10) + kappa[i]*log(j)
        pic[i,j,t]<- pi[i,j,t]/pcap[i,t]
      }
      y[i,1:J,t] ~ dmulti(pic[i,1:J,t], n[i,t])
      n[i,t] ~ dbin(pcap[i,t], N[i,t])
      nrep[i,t] ~ dbin(pcap[i,t], N[i,t])
    }
    N[i,1] ~ dpois(lambda[i,1])
    Nr[i,1]<- N[i,1] - n[i,1]
    log(lambda[i,1]) <- beta[i] + log(garea[i])
    for(t in 2:T) {
      lambda[i,t]<- exp(r[i,t-1]) * Nr[i,t-1]
      N[i,t] ~ dpois(lambda[i,t])
      Nr[i,t]<- N[i,t] - n[i,t]
    }
  }
  
  for(i in 1:K) {
    for(j in 1:(T-1)){
      r[i,j] <- eta[1] + eta[2]*elev[i] + eta[3]*season[j] + eta[4]*elev[i]*season[j]  
    }
  }
  
  for(j in 1:K){
    beta[j] ~ dnorm(0, sd = 5)
    kappa[j] ~ dnorm(muk, sd = sdk)
  }
  for(j in 1:4) {
    eta[j] ~ dnorm(0, 1)
  }
  alpha ~ dnorm(0, sd = 5)
  muk ~ dnorm(0, sd = 1)
  sdk ~ T(dnorm(0, sd=1),0,)
  
})


##---- Initial values ----

pic.init<- array(runif(K*J*T,0.05,0.15), c(K,J,T))

Nin<- n+200
rinit<- matrix(runif(K*(T-1), -0.5, 0.5),K, T-1)


inits <- function(){list(alpha=-2, p=pic.init, pi=pic.init,pic=pic.init, r=rinit,eta=rnorm(4), 
                         beta=rnorm(K), lambda=Nin, N=Nin, kappa=rnorm(K),muk=rnorm(1),sdk=runif(1))}

##---- Create the Nimble model

Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits())
Rmcmc<- compileNimble(Rmodel, showCompilerOutput = F)
ModSpec <- configureMCMC(Rmodel, onlySlice = TRUE) # Slice sampling
ModSpec$resetMonitors()
ModSpec$addMonitors(c("alpha","N","Nr","r","beta","kappa","eta","n","muk","sdk","nrep"))
Cmcmc <- buildMCMC(ModSpec)
Cmodel <- compileNimble(Cmcmc, project = Rmodel, resetFunctions = TRUE)

## ----- Run MCMC ----

ni<- 20000
nb<- 10000
nc<- 3
nt<- 1


samp<- runMCMC(Cmodel, niter = ni, nburnin = nb, nchains = nc, thin = nt, inits=inits,  
               samplesAsCodaMCMC = TRUE)


MCMCsummary(samp, params=c("beta","alpha","r","kappa","muk","sdk"))


# Operational area IDs
zones<- c("Alpine NP - Bogong High Plains","Alpine NP - Eastern Alps", "Burrowa-Pine Mt NP",
           "Croajingolong NP","Coopracambra NP", "Errinundra NP", "Mt Buffalo NP", 
           "Mt Mitta Mitta RP","Snowy River NP","Wabba Wilderness Park")


pars<- MCMCpstr(samp, params=c("N","Nr"), type="chains")

calc_N(pars, zones)

calc_N(pars, zones, A=garea)

win.graph(15,10)
mcmc_trace(samp, regex_pars = c("beta","alpha"))

win.graph(15,10)
mcmc_dens_overlay(samp, regex_pars = c("beta","alpha"))


##---- Density Plots with date---------------------------------------
session<- 1:T
period<- c("start","end")
pdates<- expand_grid(session,period)
pdates<- pdates %>% mutate(Date = as.Date(c("2020-02-10","2020-05-8","2020-06-11","2020-10-30","2021-03-01",
                         "2021-05-27","2021-09-13","2021-12-03","2022-03-07","2022-05-27")))

pars<- MCMCpstr(samp, params=c("N","Nr"), type="chains")

dout<- calc_N_period(pars, zones, A=garea)
Nhat<- dout$Nhat
Nhat<- left_join(Nhat, pdates) %>% mutate(period=factor(period, levels = c("start","end"),
                                                      labels = c("Initial","Residual")))

brks<- pdates %>% group_by(session) %>% summarise(Date=median(Date))
brk_nms<- c("Mar 2020","Aug 2020","Apr 2021","Oct 2021","Apr 2022")

win.graph(11,6)
Nhat %>% ggplot(aes(Date, N, color=period, group=1)) +
  geom_line(position = position_dodge(width=0.1)) +
  geom_pointrange(aes(ymin=nlcl, ymax=nucl), position=position_dodge(width=0.1)) +
  facet_wrap(~site, nrow=2) +
  scale_color_manual(values = c("cornflowerblue","brown1")) +
  scale_x_continuous(breaks=brks$Date, labels = brk_nms) +
  scale_y_continuous(limits=c(0, 3.1)) +
  labs(x="Period",y=expression(paste("Density (deer/k",m^2,")")), color="Estimate") +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        axis.text.x = element_text(size=8, angle = 45, hjust=1),
        axis.text.y = element_text(size=12),
        legend.position = "bottom",
        legend.title = element_text(face="bold", size=12),
        legend.text = element_text(size=12))
