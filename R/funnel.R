funnel <- function(x, y,
                   xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL,
                   comb.f=FALSE, axes=TRUE, labels=NULL, cex.lab=0.8,
                   log="", yaxis="se", sm=NULL, ...){
  
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
  
  if (match(sm, c("OR", "RR", "HR"), nomatch=0)>0){
    TE <- exp(TE)
    TE.fixed <- exp(TE.fixed)
    ##
    if (log=="") log <- "x"
  }  

  ##
  ## y-value: weight
  ##
  if (yaxis=="inv") weight <- 1/seTE^2
  if (yaxis=="se") weight <- seTE
  if (yaxis=="size")
    if (inherits(x, "metabin") || inherits(x, "metacont"))
      weight <- floor(x$n.e)+floor(x$n.c)
    else
      weight <- y
  
  
  ##
  ## x-axis: labels / xlim
  ##
  if (sm=="OR"  & is.null(xlab)) xlab <- "Odds Ratio"
  if (sm=="RD"  & is.null(xlab)) xlab <- "Risk Difference"
  if (sm=="RR"  & is.null(xlab)) xlab <- "Relative Risk"
  if (sm=="SMD" & is.null(xlab)) xlab <- "Standardised mean difference"
  if (sm=="WMD" & is.null(xlab)) xlab <- "Weighted mean difference"
  if (sm=="HR"  & is.null(xlab)) xlab <- "Hazard Ratio"
  ##
  if (is.null(xlim)) xlim <- range(TE, na.rm=TRUE)
  
  ##
  ## y-axis: labels / ylim
  ##
  if (yaxis=="se"   & is.null(ylab)) ylab  <- "Standard error"
  if (yaxis=="size" & is.null(ylab)) ylab  <- "Study size"
  if (yaxis=="inv"  & is.null(ylab)) ylab  <- "Inverse of variance"
  ##
  if (is.null(ylim) & yaxis=="se") ylim <- c(max(weight, na.rm=TRUE), 0)
  if (is.null(ylim)              ) ylim <- range(weight, na.rm=TRUE)
  
  ##
  ## plot
  ##
  plot(TE, weight, type="n",
       xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab,
       axes=axes, log=log, ...)
  
  if (is.null(labels)) points(TE, weight, pch=1)
  else text(TE, weight, labels=labels, cex=cex.lab)
  
  if (comb.f) abline(v=TE.fixed, lty=2)
   
  invisible(NULL)
}
