radial <- function(x, y,
                   xlim=NULL, ylim=NULL,
                   xlab="Inverse of standard error",
		   ylab="Standardised treatment effect (z-score)",
		   comb.f=TRUE, axes=TRUE, labels=NULL, cex.lab=0.8,
                   level=NULL, ...){

  
  if (inherits(x, "meta")){
    TE <- x$TE
    seTE <- x$seTE
    TE.f <- x$TE.f
    seTE.f <- x$seTE.f
  }
  else{
    TE <- x
    seTE <- y
  }


  if(length(TE) != length(seTE))
    stop("length of argument TE and seTE must be equal")

  if (!is.null(level) && (level<=0|level>=1))
    stop("no valid level for confidence interval")


  zscore <- TE/seTE

  if (is.null(xlim)) xlim <- c(0, max(1/seTE, na.rm=TRUE))
  if (is.null(ylim)) ylim <- c(min(c(-2, zscore), na.rm=TRUE),
                               max(c(2, zscore), na.rm=TRUE))

  plot(1/seTE, zscore,
       xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
       axes=axes, type="n", ...)

  if (is.null(labels)) points(1/seTE, zscore, pch=1)
  else text(1/seTE, zscore, labels=labels, cex=cex.lab)

  if (comb.f){
    lmcomb <- lm(zscore ~ I(1/seTE) - 1)
    abline(lmcomb, lty=4)
    if (!is.null(level)){
      alpha <- 1-level
      abline(qnorm(1-alpha/2), coef(lmcomb), lty=2)
      abline(qnorm(alpha/2), coef(lmcomb), lty=2)
    }
  }
  
  invisible(NULL)
}
