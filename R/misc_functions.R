


calc_N<- function(samples, zones, A=NULL) {
  N<- samples$N  # initial
  Nr<- samples$Nr # residual
  nsites<- length(zones)
  if(!is.null(A) && length(A) != nsites) stop("A is wrong length")
  results<- list()
  results2<- list()
  for(i in 1:nsites) {
    if(!is.null(A)) {
      Nin<- N[i,,]/A[i]
      Ninr<- Nr[i,,]/A[i]
    }
    else {
      Nin<- N[i,,]
      Ninr<- Nr[i,,]
    }
    T<- dim(Nin)[1]
    Nhat<- apply(Nin,1,mean)
    Ncl<- apply(Nin, 1, quantile, c(0.05,0.95))
    Nhatr<- apply(Ninr,1,mean)
    Nclr<- apply(Ninr, 1, quantile, c(0.05,0.95))
    dec<- round(mean(1-Ninr[T,]/Nin[1,]),3)
    deccl<- round(quantile(1-Ninr[T,]/Nin[1,], c(0.05,0.95)),3)
    df<- data.frame(site=zones[i],season=seq_len(T),N=Nhat,nlcl=Ncl[1,],nucl=Ncl[2,],
                    Nr=Nhatr,rlcl=Nclr[1,],rucl=Nclr[2,])
    df2<- data.frame(site=zones[i],dec=dec*100,lcl=deccl[1]*100,ucl=deccl[2]*100)
    results[[i]]<- df
    results2[[i]]<- df2
  }
  results<- do.call(rbind, results)
  results2<- do.call(rbind, results2)
  list(Nhat=results,dec=results2)
}


calc_N_period<- function(samples, zones, A=NULL) {
  
  N<- samples$N  # initial
  Nr<- samples$Nr # residual
  nsites<- length(zones)
  if(!is.null(A) && length(A) != nsites) stop("A is wrong length")
  results<- list()
  results2<- list()
  for(i in 1:nsites) {
    if(!is.null(A)) {
      Nin<- N[i,,]/A[i]
      Ninr<- Nr[i,,]/A[i]
    }
    else {
      Nin<- N[i,,]
      Ninr<- Nr[i,,]
    }
    T<- dim(Nin)[1]
    Nhat<- apply(Nin,1,mean)
    Ncl<- apply(Nin, 1, quantile, c(0.05,0.95))
    Nhatr<- apply(Ninr,1,mean)
    Nclr<- apply(Ninr, 1, quantile, c(0.05,0.95))
    dec<- round(mean(1-Ninr[T,]/Nin[1,]),3)
    deccl<- round(quantile(1-Ninr[T,]/Nin[1,], c(0.05,0.95)),3)
    dfs<- data.frame(site=zones[i],session=seq_len(T),period="start",N=Nhat,nlcl=Ncl[1,],nucl=Ncl[2,])
    dfe<- data.frame(site=zones[i],session=seq_len(T),period="end",N=Nhatr,nlcl=Nclr[1,],nucl=Nclr[2,])
    df<- bind_rows(dfs,dfe) %>% mutate(period = factor(period, levels = c("start","end"))) %>% 
      arrange(session,period)
    df2<- data.frame(site=zones[i],dec=dec*100,lcl=deccl[1]*100,ucl=deccl[2]*100)
    results[[i]]<- df
    results2[[i]]<- df2
  }
  results<- do.call(rbind, results)
  results2<- do.call(rbind, results2)
  list(Nhat=results,dec=results2)
}
