funnel <- function(x, y,
                   xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL,
                   comb.f=FALSE, axes=TRUE,
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
  
  seTE.max <- max(seTE, na.rm=TRUE)
  ##
  if (!is.null(level)){
    ##
    ciTE.low <- TE.fixed-qnorm(1-(1-level)/2)*seTE.max
    ciTE.upp <- TE.fixed+qnorm(1-(1-level)/2)*seTE.max
    ##
    TE.xlim <- 1.025*c(min(c(TE, ciTE.low), na.rm=TRUE),
                       max(c(TE, ciTE.upp), na.rm=TRUE))
  }
  ##
  if (match(sm, c("OR", "RR", "HR"), nomatch=0)>0){
    TE <- exp(TE)
    TE.fixed <- exp(TE.fixed)
    if (!is.null(level)){
      ciTE.low <- exp(ciTE.low)
      ciTE.upp <- exp(ciTE.upp)
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
    else if (sm=="WMD") xlab <- "Weighted mean difference"
    else if (sm=="HR" ) xlab <- "Hazard Ratio"
    else if (sm=="AS" ) xlab <- "Arcus Sinus Transformation"
    else xlab <- sm
  }
  ##
  if (is.null(xlim) & !is.null(level) & yaxis == "se")
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
  
  if (comb.f) abline(v=TE.fixed, lty=2)

  if (!is.null(level)){
    lines(c(TE.fixed, TE.fixed), c(seTE.max, 0), lty=2)
    lines(c(ciTE.upp, TE.fixed), c(seTE.max, 0), lty=2)
    lines(c(ciTE.low, TE.fixed), c(seTE.max, 0), lty=2)
  }
  
  invisible(NULL)
}
