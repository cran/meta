metabin <- function(event.e, n.e, event.c, n.c, studlab,
                    data=NULL, subset=NULL, method="MH",
                    sm=
                    ifelse(!is.na(charmatch(method, c("Peto", "peto"),
                                            nomatch = NA)),
                           "OR", "RR"),
                    incr=0.5, allincr=FALSE, addincr=FALSE, allstudies=FALSE,
                    MH.exact=FALSE, RR.cochrane=FALSE,
                    warn=TRUE){
  
  
  if (is.null(data)) data <- sys.frame(sys.parent())
  ##
  mf <- match.call()
  mf$data <- mf$subset <- mf$method <- mf$sm <- NULL
  mf$incr <- mf$allincr <- mf$allstudies <- mf$MH.exact <- NULL
  mf[[1]] <- as.name("data.frame")
  mf <- eval(mf, data)
  ##
  if (!is.null(subset)) mf <- mf[subset,]
  ##
  event.e <- mf$event.e
  n.e <- mf$n.e
  event.c <- mf$event.c
  n.c <- mf$n.c
  ##
  if (!missing(studlab)) studlab <- as.character(mf$studlab)

  
  k.all <- length(event.e)
  ##
  if (k.all == 0) stop("event.e = numeric(0)")
  ##
  if (missing(studlab)) studlab <- as.character(1:k.all)


  if (match(sm, c("OR", "RD", "RR"), nomatch=0) == 0)
    stop("possible summary measures are \"OR\", \"RD\" and \"RR\"")
  ##
  if (!(is.numeric(event.e) & is.numeric(n.e) &
        is.numeric(event.c) & is.numeric(n.c)))
    stop("Non-numeric value for event.e, n.e, event.c or n.c")
  ##
  if (any(n.e <= 0 | n.c <= 0)) stop("n.e and n.c must be positive")
  ##
  if (any(event.e < 0 | event.c < 0))
    stop("event.e and event.c must be larger equal zero")
  ##
  if (any(event.e > n.e)) stop("event.e > n.e")
  if (any(event.c > n.c)) stop("event.c > n.c")
  ##
  if (length(studlab) != k.all)
    stop("Number of studies and labels are different")


  imeth <- charmatch(tolower(method),
                     c("inverse", "mh", "peto"), nomatch = NA)
  ##
  if(is.na(imeth))
    stop("method should be \"Inverse\", \"MH\", or \"Peto\"")
  ##
  method <- c("Inverse", "MH", "Peto")[imeth]
  ##
  if (method == "Peto" & sm != "OR")
    stop("Peto's method only possible with 'sm=OR'")
  ##
  if (k.all == 1 & method == "MH"){
    warning("For single trials, inverse variance method used instead of Mantel Haenszel method.")
    method <- "Inverse"
  }

  
  ##
  ## Include non-informative trials?
  ##
  if (sm == "RD") incl <- rep(1, k.all)
  else{
    if (allstudies) incl <- rep(1, k.all)
    else
      incl <- ifelse((event.c==0 & event.e==0) |
                     (event.c==n.c & event.e==n.e), NA, 1)
  }
  ##
  ## k: effective number of trials
  ##
  k <- sum(!is.na(incl))


  ##
  ## Sparse computation
  ##
  sel <- switch(sm,
                OR=((n.e - event.e) == 0 | event.e == 0 |
                    (n.c - event.c) == 0 | event.c == 0),
                RD=((n.e - event.e) == 0 | event.e == 0 |
                    (n.c - event.c) == 0 | event.c == 0),
                RR=(event.e == 0 | event.c == 0))
  ##
  sel[is.na(incl)] <- FALSE
  ##
  sparse <- any(sel)
  ##
  if (addincr){
    i <- rep(incr, k.all)
    if (warn)
      warning(paste(incr, "added to each cell frequency of all studies"))
  }
  else{
    if (sparse){
      if (allincr){
        i <- rep(incr, k.all)
        if (warn)
          warning(paste(incr, "added to each cell frequency of all studies"))
      }
      else{
        ##
        ## Bradburn, Deeks, Altman, Stata-procedure "metan":
        ## & SAS PROC FREQ (for method="Inverse")
        ##
        i <- incr*sel
        if (warn)
          warning(paste(incr, "added to each cell in 2x2 tables with zero cell frequencies"))
      }
    }
    else i <- rep(0, k.all)
  }
  
  
  n11 <- event.e*incl
  n21 <- event.c*incl
  n1. <- n.e*incl
  n2. <- n.c*incl
  ##
  n.. <- n1. + n2.
  n12 <- n1. - n11
  n22 <- n2. - n21
  n.1 <- n11 + n21
  n.2 <- n12 + n22

  
  Q.CMH <- (sum(n11 - n1.*n.1/n.., na.rm=TRUE)^2/
            sum(n1.*n2.*n.1*n.2/n..^3, na.rm=TRUE))


  if (sm == "OR"){
    if (method == "MH" || method == "Inverse"){
      ## 
      ## Cooper & Hedges (1994), p. 251-2
      ## 
      TE <- log(((n11+i)*(n22+i)) / ((n12+i)*(n21+i)))
      varTE <- 1/(n11+i) + 1/(n12+i) + 1/(n21+i) + 1/(n22+i)
    }
    else if (method == "Peto"){
      ## 
      ## Cooper & Hedges (1994), p. 252
      ## 
      O <- n11
      E <- n1.*n.1/n..
      V <- n1.*n2.*n.1*n.2/((n..-1)*n..^2)
      ##
      TE <- (O-E)/V
      varTE <- 1/V
    }
  }
  else if (sm == "RR"){
    ## 
    ## Cooper & Hedges (1994), p. 247-8
    ## 
    if (!RR.cochrane){
      TE <- log(((n11+i)/(n1.+i))/((n21+i)/(n2.+i)))
      varTE <- 1/(n11+i) - 1/(n1.+i) + 1/(n21+i) - 1/(n2.+i)
      }
    else{
      TE <- log(((n11+i)/(n1.+2*i))/((n21+i)/(n2.+2*i)))
      varTE <- 1/(n11+i) - 1/(n1.+2*i) + 1/(n21+i) - 1/(n2.+2*i)
    }
  }
  else if (sm == "RD"){
    ## 
    ## Cooper & Hedges (1994), p. 246-7
    ## 
    ##TE <- (n11+i)/(n1.+2*i) - (n21+i)/(n2.+2*i)
    TE <- n11/n1. - n21/n2.
    varTE <- (n11+i)*(n12+i)/(n1.+2*i)^3 +
      (n21+i)*(n22+i)/(n2.+2*i)^3
  }

  
  ##
  ## Calculate random effects estimate:
  ##
  m <- metagen(TE, sqrt(varTE))
  ##
  TE.random <- m$TE.random
  seTE.random <- m$seTE.random
  w.random <- m$w.random

  
  if (method == "Inverse" || method=="Peto"){
    w.fixed <- m$w.fixed
    TE.fixed <- m$TE.fixed
    varTE.fixed <- m$seTE.fixed^2
  }
  else if (method == "MH"){
    if (!is.logical(MH.exact))
      stop("MH.exact must be of type 'logical'")
    ##
    i <- i*(!MH.exact)
    ##
    if (sm == "OR"){
      ## 
      ## Cooper & Hedges (1994), p. 253-5 (MH.exact==TRUE)
      ##
      ## Bradburn, Deeks, Altman, Stata-procedure "metan"
      ## und RevMan 3.1 (MH.exact==FALSE)
      ## 
      A <- (n11+i)*(n22+i)/(n..+4*i)
      B <- (n11+i + n22+i)/(n..+4*i)
      C <- (n12+i)*(n21+i)/(n..+4*i)
      D <- (n12+i + n21+i)/(n..+4*i)
      ##
      ## Cooper & Hedges (1994), p. 265-6
      ##
      w.fixed <- C
      TE.fixed <- log(sum(A, na.rm=TRUE)/sum(C, na.rm=TRUE))
      varTE.fixed <- (1/(2*sum(A, na.rm=TRUE)^2) *
                      (sum(A*B, na.rm=TRUE) +
                       exp(TE.fixed)*(sum(B*C, na.rm=TRUE)+
                                      sum(A*D, na.rm=TRUE)) +
                       exp(TE.fixed)^2*sum(C*D, na.rm=TRUE)))
    }
    else if (sm =="RR"){
      ##
      ## Greenland, Robins (1985) (MH.exact==TRUE)
      ##
      ## Bradburn, Deeks, Altman, Stata-procedure "metan"
      ## (MH.exact==FALSE)
      ## 
      D <- ((n1.+2*i)*(n2.+2*i)*(n.1+2*i) -
            (n11+i)*(n21+i)*(n..+4*i))/(n..+4*i)^2
      R <- (n11+i)*(n2.+2*i)/(n..+4*i)
      S <- (n21+i)*(n1.+2*i)/(n..+4*i)
      ##
      w.fixed <- S
      TE.fixed <- log(sum(R, na.rm=TRUE)/sum(S, na.rm=TRUE))
      varTE.fixed <- sum(D, na.rm=TRUE)/(sum(R, na.rm=TRUE)*sum(S, na.rm=TRUE))
    }
    else if (sm == "RD"){
      ##
      ## Jon Deeks (1999) (MH.exact==TRUE)
      ##
      ## Bradburn, Deeks, Altman, Stata-procedure "metan"
      ## und RevMan 3.1 (MH.exact==FALSE)
      ## 
      R <- ((n11+i)*(n12+i)*(n2.+2*i)^3 +
            (n21+i)*(n22+i)*(n1.+2*i)^3)/
              ((n1.+2*i)*(n2.+2*i)*(n..+4*i)^2)
      S <- n1.*n2./n..
      ##
      w.fixed <- S
      TE.fixed <- weighted.mean(TE, w.fixed, na.rm=TRUE)
      varTE.fixed <- sum(R, na.rm=TRUE)/sum(S, na.rm=TRUE)^2
    }
  }


  ##
  ## Modify fixed effects estimate:
  ##
  if (is.nan(TE.fixed)){
    TE.fixed <- NA
    varTE.fixed <- NA
  }
  ##
  w.fixed[is.nan(w.fixed)] <- 0
  w.fixed[is.na(w.fixed)] <- 0

  
  ##
  ## Cooper & Hedges (1994), p. 274-5
  ##
  if (!is.na(TE.fixed)){
    Q <- sum(1/varTE*(TE-TE.fixed)^2, na.rm=TRUE)
    ##
    if (Q<=(k-1)) tau2 <- 0
    else
      tau2 <- (Q-(k-1))/(sum(1/varTE  , na.rm=TRUE) -
                         sum(1/varTE^2, na.rm=TRUE)/
                         sum(1/varTE  , na.rm=TRUE))
  }
  else{
    Q <- NA
    tau2 <- NA
  }

  
  res <- list(event.e=event.e, n.e=n.e,
              event.c=event.c, n.c=n.c,
              studlab=studlab,
              TE=TE, seTE=sqrt(varTE),
              w.fixed=w.fixed, w.random=w.random,
              TE.fixed=TE.fixed, seTE.fixed=sqrt(varTE.fixed),
              TE.random=TE.random, seTE.random=seTE.random,
              k=k, Q=Q, tau=sqrt(tau2),
              Q.CMH=Q.CMH,
              sm=sm, method=method,
              sparse=sparse,
              incr=incr,
              allincr=allincr,
              addincr=addincr,
              allstudies=allstudies,
              MH.exact=MH.exact,
              RR.cochrane=RR.cochrane,
              call=match.call(),
              warn=warn)
  
  class(res) <- c("metabin", "meta")
  
  res
}
