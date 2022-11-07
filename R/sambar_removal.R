
library(tidyverse)
library(nimble)
library(MCMCvis)
library(abind)

source("R/misc_functions.r")

##---- load data ----

load("Data/deerCE.Rdata")

##---- construct data for Nimble ----

K<- dim(y)[1]# number of sites
T <- dim(y)[3] # number of primary periods
J<- dim(y)[2] # secondary periods per primary period

n<- apply(y, c(1,3), sum)  # marginal total removals per site and primary period


constants<- list(K=K, J=J, T=T) 

data <- list(y=y, n=n, eff=eff, garea=garea)


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
      #nrep[i,t] ~ dbin(pcap[i,t], N[i,t])
    }
    N[i,1] ~ dpois(lambda[i,1])
    Nr[i,1]<- N[i,1] - n[i,1]
    log(lambda[i,1]) <- beta[i] + log(garea[i])
    for(t in 2:T) {
      lambda[i,t]<- exp(r[t-1]) * Nr[i,t-1]
      N[i,t] ~ dpois(lambda[i,t])
      Nr[i,t]<- N[i,t] - n[i,t]
    }
  }
  
  for(j in 1:2){
    r[j] ~ dnorm(0, sd=1)
  }
  for(j in 1:K){
    beta[j] ~ dnorm(0, sd = 5)
    kappa[j] ~ dnorm(muk, sd = sdk)
  }
  
  alpha ~ dnorm(0, sd = 5)
  muk ~ dnorm(0, sd = 5)
  sdk ~ T(dnorm(0, sd=1),0,)

})

##---- Initial values ----

pic.init<- array(0.05, c(K,J,T))

Nin<- n+200

inits1 <- list(alpha=rnorm(1), p=pic.init, pi=pic.init,pic=pic.init, r=rnorm(2), 
               beta=rnorm(K), lambda=Nin, N=Nin, kappa=rnorm(K),muk=rnorm(1),sdk=runif(1))

##---- Create the Nimble model

Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits1)
Rmcmc<- compileNimble(Rmodel, showCompilerOutput = F)

ModSpec <- configureMCMC(Rmodel) # add multivariate slice sampler for beta
ModSpec$removeSamplers(c('beta'), print=FALSE)
ModSpec$addSampler(target=c('beta'), type='AF_slice')
ModSpec$resetMonitors()
ModSpec$addMonitors(c("alpha","N","Nr","r","beta","kappa","muk","sdk"))

Cmcmc <- buildMCMC(ModSpec)

Cmodel <- compileNimble(Cmcmc, project = Rmodel, resetFunctions = TRUE)

## ----- Run MCMC ----

ni<- 20000
nb<- 10000
nc<- 3
nt<- 1

# Use the following for serious inference
#ni<- 200000
#nb<- 100000
#nc<- 3
#nt<- 10

inits = function(){inits1}

samp<- runMCMC(Cmodel, niter = ni, nburnin = nb, nchains = nc, thin = nt, inits=inits,
               samplesAsCodaMCMC = TRUE)


MCMCsummary(samp, params=c("beta","alpha","r","kappa","muk","sdk"))


# Sambar Deer operational zone IDs
zones<- c("Alpine NP - Bogong High Plains","Alpine NP - Eastern Alps", "Burrowa-Pine Mt NP",
           "Croajingolong NP","Coopracambra NP", "Errinundra NP", "Mt Buffalo NP", 
           "Mt Mitta Mitta RP","Snowy River NP","Wabba Wilderness Park")


pars<- MCMCpstr(samp, params=c("N","Nr"), type="chains")

calc_N(pars, zones)

calc_N(pars, zones, A=garea)

win.graph(12,12)
diagPlot(samp, params = c("beta","alpha"))



##---- Plots---------------------------------------

dout<- calc_N(pars, zones, A=garea)
Nhat<- dout$Nhat

Nb<- Nhat %>% select(site,season,N,nlcl,nucl)
Na<- Nhat %>% select(site,season,Nr,rlcl,rucl)
names(Na)[3:5]<- c("N","nlcl","nucl")
Nb<- Nb %>% mutate(status="Initial")
Na<- Na %>% mutate(status="Residual")

Nout<- bind_rows(Nb,Na)

win.graph(11,6)
Nout %>% ggplot(aes(season, N, color=status)) +
  geom_line(position = position_dodge(width=0.1)) +
  geom_pointrange(aes(ymin=nlcl, ymax=nucl), position=position_dodge(width=0.1)) +
  facet_wrap(~site, nrow=2) +
  scale_color_manual(values = c("cornflowerblue","brown1")) +
  scale_x_continuous(breaks=c(1,2,3)) +
  labs(x="Period",y=expression(paste("Density (deer/k",m^2,")")), color="Estimate") +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        axis.text = element_text(size=12),
        legend.position = "bottom",
        legend.title = element_text(face="bold", size=12),
        legend.text = element_text(size=12))

#--- Cumulative catch curves----------------------

rmat<- abind(y[,,1],y[,,2],y[,,3])
emat<- abind(eff[,,1],eff[,,2],eff[,,3])

cumc<- t(apply(rmat,1,cumsum))
cpue<- rmat/emat

colnames(cumc)<-  paste0("P_",1:15)
colnames(cpue)<- paste0("P_",1:15)

cumc<- as_tibble(cumc)
cpue<- as_tibble(cpue)

cumc<- cumc %>% mutate(site=zones)
cpue<- cpue %>% mutate(site=zones)

cumc<- cumc %>% pivot_longer(!site, names_to = "Period", values_to = "Cumcatch")
cpue<- cpue %>% pivot_longer(!site, names_to = "Period", values_to = "CPUE")

catch<- left_join(cumc, cpue, by=c("site","Period"))

win.graph(14,6)
catch %>% filter(!(site %in% c("Brodribb area","Nunniong area"))) %>% 
  ggplot(aes(Cumcatch, CPUE)) +
  geom_point() +
  geom_smooth(method="lm") +
  facet_wrap(~site, scales="free", nrow=2) +
  labs(x ="Cumulative removals") +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        axis.text = element_text(size=12),
        legend.position = "bottom")

