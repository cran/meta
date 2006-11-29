summary.meta <- function(object,
                         byvar,
                         bylab,
                         bystud=FALSE,
                         level.comb=0.95,
                         ...){
  
  if (!inherits(object, "meta"))
    stop("Argument 'object' must be an object of class \"meta\"")

  if (inherits(object, "metainf"))
    stop("Summary method not defined for objects of class  \"metainf\"")
  ##
  if (inherits(object, "metacum"))
    stop("Summary method not defined for objects of class \"metacum\"")

  
  sm <- object$sm
  k  <- object$k
  Q  <- object$Q
  
  
  ##
  ## Higgins & Thompson (2002), Statistics in Medicine, 21, 1539-58
  ##
  H <- sqrt(Q/(k-1))
  ##
  selogH <- ifelse(k>2,
                   ifelse(Q<=k,
                          sqrt(1/(2*(k-2))*(1-1/(3*(k-2)^2))),
                          0.5*(log(Q)-log(k-1))/(sqrt(2*Q)-sqrt(2*k-3))),
                   NA)
  ##
  tres <- ci(log(H), selogH, level.comb)
  ##
  ci.H <- list(TE=max(exp(tres$TE),1),
               lower=max(exp(tres$lower),1),
               upper=max(exp(tres$upper),1))

  
  func.t <- function(x) (x^2-1)/x^2
  ##
  ci.I2 <- list(TE=func.t(ci.H$TE),
                lower=func.t(ci.H$lower),
                upper=func.t(ci.H$upper))


  ci.lab <- paste(round(100*level.comb, 1), "%-CI", sep="")
  ##
  ci.f <- ci(object$TE.fixed , object$seTE.fixed , level.comb)
  ci.r <- ci(object$TE.random, object$seTE.random, level.comb)
  
  
  if (!missing(byvar)){

    if (any(is.na(byvar))) stop("Missing values in 'byvar'")
    
    by.levs <- unique(byvar)
    
    if (missing(bylab)) bylab <- deparse(substitute(byvar))
    

    res.w <- matrix(NA, ncol=4, nrow=length(by.levs))
    j <- 0
    ##
    for ( i in by.levs){
      j <- j+1
      sel <- byvar == i
      ##
      if (inherits(object, "metabin")){
        meta1 <- metabin(object$event.e[sel], object$n.e[sel],
                         object$event.c[sel], object$n.c[sel],
                         studlab=object$studlab[sel],
                         method="Inverse",
                         sm=object$sm,
                         incr=object$incr,
                         allincr=object$allincr,
                         addincr=object$addincr,
                         allstudies=object$allstudies,
                         MH.exact=object$MH.exact,
                         RR.cochrane=object$RR.cochrane,
                         warn=object$warn)
      }
      ##
      if (inherits(object, "metacont")){
        meta1 <- metacont(object$n.e[sel], object$mean.e[sel], object$sd.e[sel],
                          object$n.c[sel], object$mean.c[sel], object$sd.c[sel],
                          sm=sm,
                          studlab=object$studlab[sel])
      }
      ##
      if (inherits(object, "metagen")){
        meta1 <- metagen(object$TE[sel], object$seTE[sel], sm=sm,
                         studlab=object$studlab[sel])
      }
      ##
      if (bystud){
        bylab2 <- paste(bylab, " = ", i, sep="")
        lab <- paste(rep("-", nchar(bylab2)), collapse="")
        ##
        cat(lab, "\n", sep="")
        cat(bylab2)
        cat("\n", lab, "\n", sep="")
        ##
        print(meta1, details=FALSE, ma=FALSE)
      }
      ##
      res.w[j,] <- c(meta1$TE.fixed, meta1$seTE.fixed,
                     meta1$Q, meta1$k)
    }
    ##
    TE.fixed.w <- res.w[,1]
    seTE.fixed.w <- res.w[,2]
    Q.w <- res.w[,3]
    k.w <- res.w[,4]
    ##
    ci.w <- ci(TE.fixed.w, seTE.fixed.w, level.comb)
    ##
    Q.b <- sum(1/seTE.fixed.w^2*(TE.fixed.w - object$TE.fixed)^2,
               na.rm=TRUE)
    ##
    if ((round(Q-sum(Q.w, na.rm=TRUE),2) - round(Q.b,2)) != 0)
      warning(paste("Q-sum(Q.w) != Q.b\nQ.b =", round(Q.b,2)))
    ##
    if (bystud) cat("\n")
  }


  res <- list(fixed=ci.f, random=ci.r,
              k=k, Q=Q, tau=object$tau, H=ci.H, I2=ci.I2,
              k.all=length(object$TE),
              Q.CMH=object$Q.CMH,
              sm=sm, method=object$method,
              call=match.call(),
              ci.lab=ci.lab)

  
  if (!missing(byvar)){
    res$within  <- ci.w
    res$k.w     <- k.w
    res$Q.w     <- Q.w
    res$bylab   <- bylab
    res$by.levs <- by.levs
  }
  
  

  if (inherits(object, "trimfill")){
    res$object <- object
    class(res) <- c("summary.meta", "trimfill")
  }
  else
    class(res) <- c("summary.meta")

  res
}
