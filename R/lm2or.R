lm2or <- function(coef, se.coef, studlab,
                  data = NULL, subset = NULL, exclude = NULL,
                  ##
                  pval, tval, lower, upper, level.ci = 0.95,
                  n, m = 1,
                  ##
                  backtransf = gs("backtransf"), ...) {
  
  
  ##
  ##
  ## (1) Read data
  ##
  ##
  nulldata <- is.null(data)
  ##
  if (nulldata)
    data <- sys.frame(sys.parent())
  ##
  mf <- match.call()
  ##
  ## Check whether essential variables are available
  ##
  if (missing(n))
    stop("Argument 'n' providing sample size(s) is mandatory.", call. = FALSE)
  ##
  missing.coef <- missing(coef)
  missing.se.coef <- missing(se.coef)
  missing.tval <- missing(tval)
  missing.pval <- missing(pval)
  missing.lower <- missing(lower)
  missing.upper <- missing(upper)
  ##
  if (missing.tval) {
    if (missing.coef)
      stop("Argument 'coef' mandatory if argument 'tval' is missing.",
           call. = FALSE)
    ##
    if (!(!missing.se.coef | !missing.pval | !(missing.lower | missing.upper)))
      stop("Either standard errors, confidence limits, or p-values ",
           "must be provided.", call. = FALSE)
  }
  ##
  ## Catch 'coef', 'se.coef', 'tval', 'pval', 'lower', 'upper',
  ## 'level.ci', and 'n' from data:
  ##
  coef <- eval(mf[[match("coef", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  chknull(coef)
  ##
  se.coef <- eval(mf[[match("se.coef", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  ##
  tval <- eval(mf[[match("tval", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  ##
  pval <- eval(mf[[match("pval", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  ##
  lower <- eval(mf[[match("lower", names(mf))]],
                data, enclos = sys.frame(sys.parent()))
  ##
  upper <- eval(mf[[match("upper", names(mf))]],
                data, enclos = sys.frame(sys.parent()))
  ##
  if (!missing(level.ci))
    level.ci <- eval(mf[[match("level.ci", names(mf))]],
                     data, enclos = sys.frame(sys.parent()))
  ##
  n <- eval(mf[[match("n", names(mf))]],
            data, enclos = sys.frame(sys.parent()))
  k.All <- length(n)
  ##
  m <- eval(mf[[match("m", names(mf))]],
            data, enclos = sys.frame(sys.parent()))
  ##
  ## Catch 'studlab', 'subset', and 'exclude' from data:
  ##
  studlab <- eval(mf[[match("studlab", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  studlab <- setstudlab(studlab, k.All)
  ##
  missing.subset <- missing(subset)
  subset <- eval(mf[[match("subset", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  ##
  missing.exclude <- missing(exclude)
  exclude <- eval(mf[[match("exclude", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  
  
  ##
  ##
  ## (2) Check length of essential variables
  ##
  ##
  arg <- "n"
  ##
  if (!missing.se.coef)
    chklength(se.coef, k.All, arg)
  if (!missing.tval)
    chklength(tval, k.All, arg)
  if (!missing.pval)
    chklength(pval, k.All, arg)
  if (!missing.lower)
    chklength(lower, k.All, arg)
  if (!missing.upper)
    chklength(upper, k.All, arg)
  if (length(level.ci) == 1)
    level.ci <- rep_len(level.ci, k.All)
  else
    chklength(level.ci, k.All, arg)
  if (length(m) == 1)
    m <- rep_len(m, k.All)
  else
    chklength(m, k.All, arg)
  ##
  ## Subset and exclude studies
  ##
  if (!missing.subset)
    if ((is.logical(subset) & (sum(subset) > k.All)) ||
        (length(subset) > k.All))
      stop("Length of argument 'subset' is larger than number of studies.")
  ##
  if (!missing.exclude)
    if ((is.logical(exclude) & (sum(exclude) > k.All)) ||
        (length(exclude) > k.All))
      stop("Length of argument 'exclude' is larger than number of studies.")
  
  
  ##
  ##
  ## (3) Check formats
  ##
  ##
  if (!missing.coef)
    chknumeric(coef)
  if (!missing.se.coef)
    chknumeric(se.coef, min = 0)
  if (!missing.tval)
    chknumeric(tval)
  if (!missing.pval)
    chklevel(pval, length = 0)
  if (!missing.lower)
    chknumeric(lower)
  if (!missing.upper)
    chknumeric(upper)
  chklevel(level.ci, length = 0)
  ##
  chknumeric(n, min = 1)
  if (!all(is.wholenumber(n)))
    stop("Sample sizes in argument 'n' must be whole numbers.", call. = FALSE)
  ##
  chknumeric(m, min = 1)
  if (!all(is.wholenumber(m[!is.na(m)])))
    stop("Argument 'm' must contain whole numbers.", call. = FALSE)
  ##
  chklogical(backtransf)
  
  
  ##
  ##
  ## (4) Calculate t-value(s)
  ##
  ##
  if (missing.tval)
    tval <- rep_len(NA, k.All)
  ##
  method.tval <- rep_len("", k.All)
  method.tval[!is.na(tval)] <- "tval"
  ##
  ## (a) Use standard errors
  ##
  sel.NA <- is.na(tval)
  if (any(sel.NA) & !missing.se.coef) {
    j <- sel.NA & !is.na(se.coef)
    method.tval[j] <- "se"
    tval[j] <- coef[j] / se.coef[j]
  }
  ##
  ## (b) Use p-values
  ##
  sel.NA <- is.na(tval)
  if (any(sel.NA) & !missing.pval) {
    j <- sel.NA & !is.na(coef) & !is.na(pval)
    method.tval[j] <- "pval"
    se.coef.j <- seTE.pval(coef[j], pval[j])$seTE
    tval[j] <- coef[j] / se.coef.j
  }
  ##
  ## (c) Use confidence limits
  ##
  sel.NA <- is.na(tval)
  if (any(sel.NA) & !missing.lower & !missing.upper) {
    j <- sel.NA & !is.na(lower) & !is.na(upper)
    method.tval[j] <- "ci"
    se.coef.j <- TE.seTE.ci(lower[j], upper[j], level.ci[j])$seTE
    tval[j] <- coef[j] / se.coef.j
  }
  
  
  ##
  ##
  ## (5) Calculate correlation(s) and log odds ratio(s)
  ##
  ##
  r <- tval / sqrt(tval^2 + (n - m - 1))
  ##
  lnOR <- log(((pi - acos(r)) / acos(r))^2)
  ##varlnOR <- (1 - r^2)^2 / (n - 1) *
  ##  (2 * pi / (sqrt(1 - r^2) * (pi * acos(r) - acos(r)^2)))^2
  varlnOR <- 4 * pi^2 * (1 - r^2) / (n - 1) / (acos(r) * (pi - acos(r)))^2
  selnOR <- sqrt(varlnOR)
  
  
  dat <- data.frame(lnOR, selnOR, OR = exp(lnOR),
                    studlab,
                    tval = tval,
                    n = n,
                    m = m,
                    method.tval = method.tval,
                    stringsAsFactors = FALSE)
  ##
  if (!missing.coef)
    dat$coef <- coef
  if (!missing.se.coef)
    dat$se.coef <- se.coef
  if (!missing.pval)
    dat$pval <- pval
  if (!missing.lower)
    dat$lower <- lower
  if (!missing.upper)
    dat$upper <- upper
  dat$level.ci <- level.ci
  ##
  dat$subset <- subset
  dat$exclude <- exclude
  ##
  res <- metagen(lnOR, selnOR, studlab = studlab,
                 data = dat, subset = subset, exclude = exclude,
                 sm = "OR", backtransf = backtransf, ...)
  
  res
}
