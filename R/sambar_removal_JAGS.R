
library(tidyverse)
library(jagsUI)
library(MCMCvis)
library(abind)

source("R/misc_functions.r")

##---- load data ----

load("Data/deerCE.Rdata")

##---- construct data for JAGS ----

K<- dim(y)[1]# number of sites
T <- dim(y)[3] # number of primary periods
J<- dim(y)[2] # secondary periods per primary period

n<- apply(y, c(1,3), sum)  # marginal total removals per site and primary period

jags.data <- list(K=K, J=J, T=T, y=y, n=n, eff=eff, garea=garea, elev=elevation, season=season)

##---- JAGS model ---------------------------------------

sink("removal_model.jags")

cat("
model{
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
    beta[j] ~ dnorm(0, 0.2)
    kappa[j] ~ dnorm(muk, sprec)
  }
  for(j in 1:4) {
    eta[j] ~ dnorm(0, 0.2)
  }
  alpha ~ dnorm(0, 0.2)
  muk ~ dnorm(0, 0.2)
  sprec <- pow(sdk, -2)
  sdk ~ dunif(0, 5)
  
}
", fill=TRUE)
sink()

##----- Inits


pic.init<- array(runif(K*J*T,0.05,0.15), c(K,J,T))

Nin<- n+200
rinit<- matrix(runif(K*(T-1), -0.5, 0.5),K, T-1)


inits <- function(){list(alpha=-2,eta=rnorm(4),beta=rnorm(K), N=Nin, kappa=rnorm(K),muk=rnorm(1),sdk=runif(1))}

parameters<- c("alpha","N","Nr","r","beta","kappa","eta","muk","sdk")

ni<- 200000
nb<- 100000
nc<- 3
nt<- 10

samp<- jags(jags.data, inits, parameters, "removal_model.jags", n.chains = nc, n.thin = nt, 
            n.iter=ni, n.burnin=nb, parallel = TRUE)


MCMCsummary(samp, params=c("beta","alpha","r","kappa","muk","sdk"))


# Sambar Deer operational zone IDs
zones<- c("Alpine NP - Bogong High Plains","Alpine NP - Eastern Alps", "Burrowa-Pine Mt NP",
          "Croajingolong NP","Coopracambra NP", "Errinundra NP", "Mt Buffalo NP", 
          "Mt Mitta Mitta RP","Snowy River NP","Wabba Wilderness Park")


pars<- MCMCpstr(samp, params=c("N","Nr"), type="chains")

calc_N(pars, zones)

calc_N(pars, zones, A=garea)

win.graph(15,10)
jagsUI::traceplot(samp, parameters = c("beta","alpha"))


