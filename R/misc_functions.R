

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
    Nhat<- apply(Nin,1,median)
    Ncl<- apply(Nin, 1, quantile, c(0.05,0.95))
    Nhatr<- apply(Ninr,1,median)
    Nclr<- apply(Ninr, 1, quantile, c(0.05,0.95))
    dec<- round(mean(1-Ninr[3,]/Nin[1,]),3)
    deccl<- round(quantile(1-Ninr[3,]/Nin[1,], c(0.05,0.95)),3)
    df<- data.frame(site=zones[i],season=c(1,2,3),N=Nhat,nlcl=Ncl[1,],nucl=Ncl[2,],
                    Nr=Nhatr,rlcl=Nclr[1,],rucl=Nclr[2,])
    df2<- data.frame(site=zones[i],dec=dec*100,lcl=deccl[1]*100,ucl=deccl[2]*100)
    results[[i]]<- df
    results2[[i]]<- df2
  }
  results<- do.call(rbind, results)
  results2<- do.call(rbind, results2)
  list(Nhat=results,dec=results2)
}

