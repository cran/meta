funnel <- function(x, y,
                   xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL,
                   comb.fixed=FALSE, axes=TRUE,
                   pch=1, text=NULL, cex=1, col=NULL,
                   log="", yaxis="se", sm=NULL,
                   level=NULL, ...){
  
  if (inherits(x, "meta")){
    TE <- x$TE
    seTE <- x$seTE
    TE.fixed <- x$TE.fixed
    seTE.fixed <- x$seTE.fixed
    sm <- x$sm
  }
  else{
    TE <- x
    seTE <- y
    TE.fixed <- sum(TE*(1/seTE^2), na.rm=TRUE)/sum(1/seTE^2, na.rm=TRUE)
    if (is.null(sm)) sm <- "other"
  }
  
  if(length(TE) != length(seTE))
    stop("length of argument TE and seTE must be equal")
  
  if (!is.null(level) && (level<=0|level>=1))
    stop("no valid level for confidence interval")
  
  iyaxis <- charmatch(yaxis,
                     c("se", "size", "invvar", "invse"),
                     nomatch = NA)
  if(is.na(iyaxis) | iyaxis==0)
    stop("yaxis should be \"se\", \"size\", \"invvar\", or \"invse\"")
  ##
  yaxis <- c("se", "size", "invvar", "invse")[iyaxis]
  
  
  seTE[is.infinite(seTE)] <- NA
  
  if ( yaxis == "se" )
    seTE.min <- 0
  else
    seTE.min <- min(seTE, na.rm=TRUE)
  ##
  seTE.max <- max(seTE, na.rm=TRUE)
  ##
  if (!is.null(level)){
    ##
    seTE.seq <- seq(seTE.min, seTE.max, length.out=100)
    ##
    ciTE <- ci(TE.fixed, seTE.seq, level)
    ##
    TE.xlim <- 1.025*c(min(c(TE, ciTE$lower), na.rm=TRUE),
                       max(c(TE, ciTE$upper), na.rm=TRUE))
  }
  ##
  if (match(sm, c("OR", "RR", "HR"), nomatch=0)>0){
    TE <- exp(TE)
    TE.fixed <- exp(TE.fixed)
    if (!is.null(level)){
      ciTE$lower <- exp(ciTE$lower)
      ciTE$upper <- exp(ciTE$upper)
      TE.xlim <- exp(TE.xlim)
    }
    ##
    if (log=="") log <- "x"
  }
  
  
  ##
  ## y-value: weight
  ##
  if (yaxis=="invvar") weight <- 1/seTE^2
  if (yaxis=="invse")  weight <- 1/seTE
  if (yaxis=="se") weight <- seTE
  if (yaxis=="size")
    if (inherits(x, "metabin") || inherits(x, "metacont"))
      weight <- floor(x$n.e)+floor(x$n.c)
    else
      weight <- y
  
  
  ##
  ## x-axis: labels / xlim
  ##
  if (is.null(xlab)){
    if      (sm=="OR" ) xlab <- "Odds Ratio"
    else if (sm=="RD" ) xlab <- "Risk Difference"
    else if (sm=="RR" ) xlab <- "Relative Risk"
    else if (sm=="SMD") xlab <- "Standardised mean difference"
    else if (sm=="WMD"|sm=="MD") xlab <- "Mean difference"
    else if (sm=="HR" ) xlab <- "Hazard Ratio"
    else if (sm=="AS" ) xlab <- "Arcus Sinus Transformation"
    else if (sm=="proportion" ){
      if (inherits(x, "metaprop"))
        if (!x$freeman.tukey)
          xlab <- "Proportion (arcsine transformation)"
        else
          xlab <- "Proportion (Freeman-Tukey double arcsine transformation)"
    }
    else xlab <- sm
  }
  ##
  if (is.null(xlim) & !is.null(level) &
      (yaxis == "se" |
       yaxis == "invse" |
       yaxis == "invvar"))
    xlim <- TE.xlim
  else if (is.null(xlim))
    xlim <- range(TE, na.rm=TRUE)
  
  
  ##
  ## y-axis: labels / ylim
  ##
  if (yaxis=="se"   & is.null(ylab)) ylab <- "Standard error"
  if (yaxis=="size" & is.null(ylab)) ylab <- "Study size"
  if (yaxis=="invvar" & is.null(ylab)) ylab <- "Inverse of variance"
  if (yaxis=="invse" & is.null(ylab)) ylab <- "Inverse of standard error"
  ##
  if (is.null(ylim) & yaxis=="se") ylim <- c(max(weight, na.rm=TRUE), 0)
  if (is.null(ylim)              ) ylim <- range(weight, na.rm=TRUE)
  
  ##
  ## plot
  ##
  plot(TE, weight, type="n",
       xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab,
       axes=axes, log=log, ...)
  
  if (is.null(text)){
    if (is.null(col))
      points(TE, weight, pch=pch, cex=cex)
    else
      points(TE, weight, pch=pch, cex=cex, col=col)
  }
  else{
    if (is.null(col))
      text(TE, weight, labels=text, cex=cex)
    else
      text(TE, weight, labels=text, cex=cex, col=col)
  }
  
  if (comb.fixed)
    abline(v=TE.fixed, lty=2)
  
  if (!is.null(level) & yaxis=="se"){
    points(ciTE$lower, seTE.seq, type="l", lty=2)
    points(ciTE$upper, seTE.seq, type="l", lty=2)
  }
  
  if (!is.null(level) & yaxis=="invvar"){
    points(ciTE$lower, 1/seTE.seq^2, type="l", lty=2)
    points(ciTE$upper, 1/seTE.seq^2, type="l", lty=2)
  }
  
  if (!is.null(level) & yaxis=="invse"){
    points(ciTE$lower, 1/seTE.seq, type="l", lty=2)
    points(ciTE$upper, 1/seTE.seq, type="l", lty=2)
  }
  
  invisible(NULL)
}
