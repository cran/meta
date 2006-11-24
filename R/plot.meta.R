plot.meta <- function(x,
                      byvar, bylab, sortvar,
                      studlab=TRUE,
                      level=0.95, level.comb=level,
                      comb.f=FALSE, comb.r=FALSE,
                      text.f="Fixed effect model",
                      text.r="Random effects model",
                      lty.f=2, lty.r=3,
                      xlab=NULL,
                      xlim, ylim, lwd=1, cex=1,
                      cex.comb=1.2*cex,
                      cex.axis=cex,
                      cex.lab=cex,
		      log=ifelse(x$sm=="RR"|x$sm=="OR"|x$sm=="HR", "x", ""),
                      axes=TRUE,
                      allstudies=TRUE,
                      weight="fixed",
                      scale.diamond=1,
                      scale.square=1,
                      col.i="black",
                      ...){

  
  if (!inherits(x, "meta"))
    stop("Argument 'x' must be an object of class \"meta\"")

  
  oldpar <- par(las=1, mar=c(5.1, 2.1, 4.1, 2.1))
  on.exit(par(oldpar))


  if (inherits(x, "metainf")|inherits(x, "metacum")){
    x$TE.fixed    <- rev(x$TE)[1]
    x$seTE.fixed  <- rev(x$seTE)[1]
    x$TE.random   <- rev(x$TE)[1]
    x$seTE.random <- rev(x$seTE)[1]
    ##
    x$TE <- rev(rev(x$TE)[-(1:2)])
    x$seTE <- rev(rev(x$seTE)[-(1:2)])
    x$studlab <- rev(rev(x$studlab)[-(1:2)])
    ##
    if (x$pooled=="fixed") text.f <- "Fixed effect model"
    else  text.f <- "Random effects model"
    comb.f <- TRUE
  }
  
  sm <- x$sm

  
  if (is.null(xlab)){
    if      (sm=="OR" ) xlab <- "Odds Ratio"
    else if (sm=="RD" ) xlab <- "Risk Difference"
    else if (sm=="RR" ) xlab <- "Relative Risk"
    else if (sm=="SMD") xlab <- "Standardised mean difference"
    else if (sm=="WMD") xlab <- "Weighted mean difference"
    else if (sm=="HR" ) xlab <- "Hazard Ratio"
    else if (sm=="AS" ) xlab <- "Arcus Sinus Transformation"
    else xlab <- sm
  }
  ##  if (inherits(x, "meta")){
  ##    if (sm == "OR"  & missing(xlab)) xlab <- "Odds Ratio"
  ##    else if (sm == "RD"  & missing(xlab)) xlab <- "Risk Difference"
  ##    else if (sm == "RR"  & missing(xlab)) xlab <- "Relative Risk"
  ##    else if (sm == "SMD" & missing(xlab)) xlab <- "Standardised mean difference"
  ##    else if (sm == "WMD" & missing(xlab)) xlab <- "Weighted mean difference"
  ##    else if (sm == "HR"  & missing(xlab)) xlab <- "Hazard Ratio"
  ##    else if (sm == "AS"  & missing(xlab)) xlab <- "Arcus Sinus Transformation"
  ##    else xlab <- sm
  ##  }

  
  ##
  ## total number of trials to plot (*not* number of trials combined)
  ##
  k.all <- length(x$TE)


  if (allstudies) n.stud <- k.all
  else n.stud <- x$k # number of trials combined in meta-analysis

  
  by <- !missing(byvar)
  sort <- !missing(sortvar)

  if (!by) byvar <- rep(1, k.all)
  if (!sort) sortvar <- rep(1, k.all)

  if (sort & length(sortvar) != k.all)
    stop("'x' and 'sortvar' have different length")
  if (by & length(byvar) != k.all)
    stop("'x' and 'byvar' have different length")
  if (by & any(is.na(byvar)))
    stop("Missing values in 'byvar'")


  if (length(studlab) == 1 & is.logical(studlab))
    if (studlab == FALSE) studlab <- rep("", k.all)
    else studlab <- x$studlab
  ##
  if (length(studlab) != k.all)
    stop("'x' and 'studlab' have different length")
  
  if (allstudies)
    sel <- 1:length(x$TE)
  else
    sel <- !is.na(x$TE)
  ##
  if (n.stud != sum(sel>0)) warning("n.stud != sum(sel)")
  ##
  n.e <- x$n.e[sel]
  n.c <- x$n.c[sel]
  ##
  event.e <- x$event.e[sel]
  event.c <- x$event.c[sel]
  ##
  mean.e <- x$mean.e[sel]
  mean.c <- x$mean.c[sel]
  ##
  sd.e <- x$sd.e[sel]
  sd.c <- x$sd.c[sel]
  ##
  TE <- x$TE[sel]
  seTE <- x$seTE[sel]
  ##
  w.fixed <- x$w.fixed[sel]
  w.random <- x$w.random[sel]
  studlab <- studlab[sel]
  ##
  byvar <- byvar[sel]
  sortvar <- sortvar[sel]

  
  if (sort | by){
    ##
    o <- order(byvar, sortvar)
    ##
    n.e <- n.e[o]
    n.c <- n.c[o]
    ##
    event.e <- event.e[o]
    event.c <- event.c[o]
    ##
    mean.e <- mean.e[o]
    mean.c <- mean.c[o]
    ##
    sd.e <- sd.e[o]
    sd.c <- sd.c[o]
    ##
    TE <- TE[o]
    seTE <- seTE[o]
    ##
    w.fixed <- w.fixed[o]
    w.random <- w.random[o]
    studlab <- studlab[o]
    ##
    byvar <- byvar[o]
    sortvar <- sortvar[o]
  }

  
  by.levs <- unique(byvar)

  
  tres <- ci(TE, seTE, level)
  TE <- tres$TE
  lowTE <- tres$lower
  uppTE <- tres$upper
  ##
  tres <- ci(x$TE.fixed, x$seTE.fixed, level.comb)
  TE.fixed <- tres$TE
  lowTE.fixed <- tres$lower
  uppTE.fixed <- tres$upper
  ##
  tres <- ci(x$TE.random, x$seTE.random, level.comb)
  TE.random <- tres$TE
  lowTE.random <- tres$lower
  uppTE.random <- tres$upper

  
  if (by){

    res.w <- matrix(NA, ncol=7, nrow=length(by.levs))
    j <- 0
    ##
    for ( i in by.levs){
      j <- j+1
      sel <- byvar == i
      ##
      if (inherits(x, "metabin")){
        meta1 <- metabin(event.e[sel], n.e[sel], event.c[sel], n.c[sel],
                         method=x$method, sm=sm,
                         incr=x$incr, allincr=x$allincr,
                         addincr=x$addincr, allstudies=x$allstudies,
                         MH.exact=x$MH.exact, RR.cochrane=x$RR.cochrane,
                         warn=x$warn)
      }
      ##
      if (inherits(x, "metacont")){
        meta1 <- metacont(n.e[sel], mean.e[sel], sd.e[sel],
                          n.c[sel], mean.c[sel], sd.c[sel],
                          sm=sm)
      }
      ##
      if (inherits(x, "metagen")){
        meta1 <- metagen(TE[sel], seTE[sel], sm=sm)
      }
      ##
      res.w[j,] <- c(meta1$TE.fixed, meta1$seTE.fixed,
                     meta1$TE.random, meta1$seTE.random,
                     meta1$Q, meta1$k, length(meta1$TE))
      ##
    }
    ##
    k.w     <- res.w[,6]
    k.all.w <- res.w[,7]
    ##
    if (allstudies) sel <- 1:length(k.w)
    else sel <- k.w>0
    ##
    TE.fixed.w    <- res.w[sel,1]
    seTE.fixed.w  <- res.w[sel,2]
    TE.random.w   <- res.w[sel,3]
    seTE.random.w <- res.w[sel,4]
    Q.w           <- res.w[sel,5]
    ##
    by.levs <- by.levs[sel]
    k.all.w <- k.all.w[sel]
    ##
    if (comb.f | !(comb.f|comb.r)){
      if (comb.r) warning("FE estimate used in groups defined by 'byvar'")
      ##
      tres <- ci(TE.fixed.w, seTE.fixed.w, level)
    }
    if (!comb.f & comb.r)
      tres <- ci(TE.random.w, seTE.random.w, level)
    ##
    TE.w <- tres$TE
    lowTE.w <- tres$lower
    uppTE.w <- tres$upper
  }
  
  
  if (by) n.by <- length(by.levs)
  else n.by <- 0

  
  if (comb.f & comb.r){
    dy.comb <- 3.0
    yTE.fixed <- -0.5
    yTE.random <- -2.0
  }
  if (comb.f & !comb.r){
    dy.comb <- 1.5
    yTE.fixed <- -0.5
  }
  if (!comb.f & comb.r){
    dy.comb <- 1.5
    yTE.random <- -0.5
  }
  if (!comb.f & !comb.r){
    dy.comb <- 0
  }
  
  ##
  ## Plot symbol
  ##
  if (weight=="same"|inherits(x, "metainf")|inherits(x, "metacum"))
    cex.i <- rep(cex, n.stud)
  else
    if (weight=="fixed")
      cex.i <- 3*w.fixed/max(w.fixed)*scale.square
    else
      cex.i <- 3*w.random/max(w.random)*scale.square


  ##
  ## y-axis:
  ##
  rows.w <- 3
  N <- n.stud + rows.w*n.by
  ##
  if (!by){
    yTE <- N:1
    yTE.w <- NA
  }
  else{
    ##
    ## Don't ask me why I did it this way, but it works.
    ##
    sel <- rep(rep(c(FALSE,TRUE), n.by),
               as.vector(rbind(rep(rows.w,n.by),
                               rev(k.all.w))))
    ##
    yTE <- rev((1:N)[sel])
    yTE.w <- rev((1:N)[!sel][2+cumsum(c(0, rep(rows.w, n.by-1)))])
    ##
    ##
    if (!comb.f & !comb.r){
      yTE <- yTE - 1
      yTE.w <- yTE.w - 1
      N <- N-1
    }
  }
  ##
  if (missing(ylim))
    ylim <- c(0.5-dy.comb, N+0.5)

  dy <- 0.25*scale.diamond

  
  if (sm == "RR" | sm == "OR" | sm == "HR"){
    TE.null <- 1
    ##
    TE <- exp(TE)
    lowTE <- exp(lowTE)
    uppTE <- exp(uppTE)
    ##
    TE.fixed <- exp(TE.fixed)
    lowTE.fixed <- exp(lowTE.fixed)
    uppTE.fixed <- exp(uppTE.fixed)
    ##
    TE.random <- exp(TE.random)
    lowTE.random <- exp(lowTE.random)
    uppTE.random <- exp(uppTE.random)
    ##
    if (by){
      TE.w <- exp(TE.w)
      lowTE.w <- exp(lowTE.w)
      uppTE.w <- exp(uppTE.w)
    }
  }
  else TE.null <- 0

  
  ##
  ## x-axis:
  ##
  if (missing(xlim))
    xlim <- c(min(lowTE, na.rm=TRUE), max(uppTE, na.rm=TRUE))

  
  if (missing(bylab)) bylab <- paste("byvar", by.levs, sep=" = ")
  else bylab <- paste(bylab, by.levs, sep=" = ")

  
  ##
  ## Forest plot
  ##
  plot(NA, NA, xlim=xlim, ylim=ylim,
       xlab="", ylab="", type="n", log=log,
       axes=FALSE, ...)
  
  ##
  ## Add axis
  ##
  if (axes){
    axis(1, cex.axis=cex.axis)
    mtext(xlab, 1, line=2, cex=cex.lab)
  }

  ##
  ## Add no effect line
  ##
  ##abline(v=TE.null, lty=2, lwd=2, col="dimgray")
  abline(v=TE.null, lty=1, lwd=lwd)

  
  ##
  ## Add line at pooled effect
  ##
  if (comb.f)
    segments(TE.fixed, yTE.fixed+dy, TE.fixed, N,
             lty=lty.f, lwd=2, col="dimgray")
  if (comb.r)
    segments(TE.random, yTE.random+dy, TE.random, N,
             lty=lty.r, lwd=2, col="dimgray")
  ##  if (inherits(x, "metainf")|inherits(x, "metacum"))
  ##    abline(v=TE.fixed, lty=3, lwd=2, col="dimgray")

  ##
  ## Add single trial results
  ##
  if (length(col.i) == 1) col.i <- rep(col.i, length(TE))
  for ( i in seq(along=TE)){
    if (cex.i[i] < 0.3)
      points(TE[i], yTE[i], pch=3, col=col.i[i], lwd=lwd)
    else
      points(TE[i], yTE[i], cex=as.numeric(cex.i[i]), pch=15, col=col.i[i])
  }
  ##
  ## Plot confidence intervals
  ##
  segments(lowTE, yTE, uppTE, yTE, lwd=lwd, col=col.i)
  ##
  ## Add trial labels
  ##
  mtext(studlab, side=2, line=-0.2, outer=FALSE,
        at=yTE+0.25, adj=0, cex=0.7*cex)

  ##
  ## Add group results
  ##
  if (by){
    ##
    for ( i in 1:length(lowTE.w)){
      polygon(c(lowTE.w[i], TE.w[i], uppTE.w[i], TE.w[i]),
              c(yTE.w[i], yTE.w[i]-dy, yTE.w[i], yTE.w[i]+dy),
              col="lightgray")
    }
    ##
    mtext(bylab,
          side=2, line=-0.2, outer=FALSE,
          at=yTE.w+0.25,
          adj=0, cex=0.8*cex)
  }

  ##
  ## Add pooled estimates
  ##
  if (comb.f){
    ##
    ## FE estimate
    ##
    polygon(c(lowTE.fixed, TE.fixed, uppTE.fixed, TE.fixed),
            c(yTE.fixed, yTE.fixed-dy, yTE.fixed, yTE.fixed+dy),
            col="black")

    mtext(text.f,
          side=2, line=-0.2, outer=FALSE,
          at=yTE.fixed+0.25, adj=0, cex=cex.comb)
  }
  ##
  if (comb.r){
    ##
    ## RE estimate
    ##
    polygon(c(lowTE.random, TE.random, uppTE.random, TE.random),
            c(yTE.random, yTE.random-dy, yTE.random, yTE.random+dy),
            col="darkgray")
    ##
    mtext(text.r, 
          side=2, line=-0.2, outer=FALSE,
          at=yTE.random+0.25, adj=0, cex=cex.comb)
  }
  
  invisible(NULL)
}
