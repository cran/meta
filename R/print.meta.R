print.meta <- function(x,
                       sortvar,
                       level=0.95,
                       level.comb=level,
                       details=FALSE,
                       ma=TRUE,
                       digits=max(0, .Options$digits - 3),
                       ...){
    
  if (!inherits(x, "meta"))
    stop("Argument 'x' must be an object of class \"meta\"")
  
  
  format.TE <- function(TE, na=FALSE){
    TE <- rmSpace(TE)
    if (na) res <- format(TE)
    else res <- ifelse(is.na(TE), "", format(TE))
    res
  }


  k.all <- length(x$TE)
  
  if (missing(sortvar)) sortvar <- 1:k.all
  
  if (length(sortvar) != k.all)
    stop("'x' and 'sortvar' have different length")


  ci.lab <- paste(round(100*level, 1), "%-CI", sep="")
  
  
  if (details){

    if (inherits(x, "metabin")){
      res <- cbind(event.e=x$event.e, n.e=x$n.e,
                   event.c=x$event.c, n.c=x$n.c,
                   p.e=round(x$event.e/x$n.e, digits),
                   p.c=round(x$event.c/x$n.c, digits))
    }
    
    else if (inherits(x, "metacont")){
      res <- cbind(n.e=x$n.e, mean.e=x$mean.e, sd.e=x$sd.e,
                   n.c=x$n.c, mean.c=x$mean.c, sd.c=x$sd.c)
    }
    
    else{
      res <- cbind(TE=x$TE, seTE=x$seTE)
    }
    
    dimnames(res)[[1]] <- x$studlab
    
    prmatrix(res[order(sortvar),])
    cat("\n\n")
  }
  
  
  if (k.all == 1){
    cat("Trial-ID:", x$studlab, "\n\n")
    if (inherits(x, "metabin"))
      print(summary(metabin(x$event.e, x$n.e,
                            x$event.c, x$n.c,
                            studlab=x$studlab,
                            method="Inverse",
                            sm=x$sm,
                            incr=x$incr,
                            allincr=x$allincr,
                            addincr=x$addincr,
                            allstudies=x$allstudies,
                            MH.exact=x$MH.exact,
                            RR.cochrane=x$RR.cochrane,
                            warn=x$warn), level.comb=level.comb), digits=digits)
    else
      print(summary(x, level.comb=level.comb), digits=digits)
  }
  else{
    tres <- ci(x$TE, x$seTE, level)
    ##
    if (x$sm == "RR" | x$sm == "OR" | x$sm == "HR"){
      TE <- round(exp(tres$TE), digits)
      lowTE <- round(exp(tres$lower), digits)
      uppTE <- round(exp(tres$upper), digits)
    }
    else{
      TE <- round(tres$TE, digits)
      lowTE <- round(tres$lower, digits)
      uppTE <- round(tres$upper, digits)
    }
    
    if (sum(x$w.fixed)>0)
     w.fixed.p <- 100*round(x$w.fixed/sum(x$w.fixed, na.rm=TRUE), 4)
    else w.fixed.p <- x$w.fixed

    if (sum(x$w.random)>0)
     w.random.p <- 100*round(x$w.random/sum(x$w.random, na.rm=TRUE), 4)
    else w.random.p <- x$w.random
    
    if (inherits(x, "metainf")|inherits(x, "metacum")){
      I2 <- 100*x$I2
      sel1 <- is.na(I2)
      I2 <- ifelse(sel1, "", format(round(I2, 1)))
      ##
      sel2 <- is.na(x$p.value)
      p.value <- format.p(x$p.value)
      p.value <- ifelse(sel2, "", p.value)
      ##
      res <- cbind(format.TE(TE), p.ci(format(lowTE), format(uppTE)),
                   p.value,
                   paste("    ", I2, ifelse(sel1, "", "%"), sep=""))
      dimnames(res) <- list(paste(x$studlab, "  ", sep=""),
                            c(x$sm, ci.lab, "p.value", "I^2"))
      
      if (inherits(x, "metainf")){
        if (x$pooled=="fixed")
          cat("\nInfluential analysis (Fixed effect model)\n")
        else
          cat("\nInfluential analysis (Random effects model)\n")
      }
      
      if (inherits(x, "metacum")){
        if (x$pooled=="fixed")
          cat("\nCumulative meta-analysis (Fixed effect model)\n")
        else
          cat("\nCumulative meta-analysis (Random effects model)\n")
      }
      cat("\n")
      prmatrix(res, quote=FALSE, right=TRUE)
      
      method <- ifelse(x$method=="MH",
                       "Mantel-Haenszel method",
                       ifelse(x$method=="Peto", "Peto method",
                              ifelse(x$method=="Inverse",
                                     "Inverse variance method",
                                     x$method)))
      
      cat(paste("\nMethod:", method, "\n"))
    }
    else{
      res <- cbind(format.TE(TE, na=TRUE), p.ci(format(lowTE), format(uppTE)),
                   format(w.fixed.p), format(w.random.p))
      dimnames(res) <-
        list(x$studlab, c(x$sm, ci.lab, "%W(fixed)", "%W(random)"))
      prmatrix(res[order(sortvar),], quote=FALSE, right=TRUE)
      cat("\n")
    }
    
    
    if (ma&!(inherits(x, "metainf")|inherits(x, "metacum")))
      print(summary(x, level.comb=level.comb), digits=digits)
    
  }
  
  
  invisible(NULL)
}
