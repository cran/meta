metainf <- function(x, pooled="fixed", sortvar){

  
  if (!inherits(x, "meta"))
    stop("Argument 'x' must be an object of class \"meta\"")


  imeth <- charmatch(tolower(pooled), c("fixed", "random"), nomatch = NA)
  ##
  if (is.na(imeth)) 
        stop("'pooled' should be \"fixed\" or \"random\"")
  ##
  pooled <- c("fixed", "random")[imeth]


  k.all <- length(x$TE)
  sort <- !missing(sortvar)

  if (!sort) sortvar <- rep(1, k.all)

  if (sort & length(sortvar) != k.all)
    stop("'x' and 'sortvar' have different length")

  
  n.e <- x$n.e
  n.c <- x$n.c
  ##
  event.e <- x$event.e
  event.c <- x$event.c
  ##
  mean.e <- x$mean.e
  mean.c <- x$mean.c
  ##
  sd.e <- x$sd.e
  sd.c <- x$sd.c
  ##
  TE <- x$TE
  seTE <- x$seTE
  ##
  studlab <- x$studlab
  sortvar <- sortvar
  
  
  if (sort){
    ##
    o <- order(sortvar)
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
    studlab <- studlab[o]
    sortvar <- sortvar[o]
  }

  
  res.i <- matrix(NA, ncol=2, nrow=k.all)
  ##
  for (i in 1:k.all){
    sel <- -i
    ##
    if (inherits(x, "metabin")){
      m <- metabin(event.e[sel], n.e[sel], event.c[sel], n.c[sel],
                   method=x$method, sm=x$sm,
                   incr=x$incr, allincr=x$allincr, addincr=x$addincr,
                   allstudies=x$allstudies, MH.exact=x$MH.exact,
                   RR.cochrane=x$RR.cochrane,
                   warn=x$warn)
    }
    ##
    if (inherits(x, "metacont")){
      m <- metacont(n.e[sel], mean.e[sel], sd.e[sel],
                    n.c[sel], mean.c[sel], sd.c[sel], sm=x$sm)
    }
    ##
    if (inherits(x, "metagen")){
      m <- metagen(TE[sel], seTE[sel], sm=x$sm)
    }
    ##
    if (pooled == "fixed"){
      res.i[i,] <- c(m$TE.fixed, m$seTE.fixed)
      TE.sum <- x$TE.fixed
      seTE.sum <- x$seTE.fixed
    }
    if (pooled == "random"){
      res.i[i,] <- c(m$TE.random, m$seTE.random)
      TE.sum <- x$TE.random
      seTE.sum <- x$seTE.random
    }
  }
  
  
  slab <- c(paste("Omitting", studlab), "Pooled estimate")
  
  res <- list(TE=c(res.i[,1], NA, TE.sum),
              seTE=c(res.i[,2], NA, seTE.sum),
              studlab=c(rev(rev(slab)[-1]), " ", rev(slab)[1]),
              sm=x$sm, method=x$method, k=x$k,
              pooled=pooled)
  
  class(res) <- c("metainf", "meta")
  
  res
}
