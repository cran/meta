#' L'Abbé plot for meta-analysis with binary outcomes
#' 
#' @description
#' Draw a L'Abbé plot for meta-analysis with binary outcomes.
#' 
#' @aliases labbe labbe.metabin labbe.default
#' 
#' @param x An object of class \code{metabin}. Alternatively, the x
#'   coordinates of points of the L'Abbé plot.
#' @param y The y coordinates of the L'Abbé plot, if argument \code{x}
#'   is not an object of class \code{metabin}.
#' @param xlim The x limits (min, max) of the plot.
#' @param ylim The y limits (min, max) of the plot.
#' @param xlab A label for the x-axis.
#' @param ylab A label for the y-axis.
#' @param TE.fixed A numeric or vector specifying combined fixed
#'   effect estimate(s).
#' @param TE.random A numeric or vector specifying combined random
#'   effects estimate(s).
#' @param comb.fixed A logical indicating whether the pooled fixed
#'   effect estimate should be plotted.
#' @param comb.random A logical indicating whether the pooled random
#'   effects estimate should be plotted.
#' @param backtransf A logical indicating which values should be
#'   printed on x- and y-axis (see Details).
#' @param axes A logical indicating whether axes should be drawn on
#'   the plot.
#' @param pch The plotting symbol used for individual studies.
#' @param text A character vector specifying the text to be used
#'   instead of plotting symbol.
#' @param cex The magnification to be used for plotting symbol.
#' @param col A vector with colour of plotting symbols.
#' @param bg A vector with background colour of plotting symbols (only
#'   used if \code{pch} in \code{21:25}).
#' @param lwd The line width.
#' @param lwd.fixed The line width(s) for fixed effect estimate(s) (if
#'   \code{comb.fixed} is not \code{NULL} or \code{FALSE}).
#' @param lwd.random The line width(s) for random effects estimate(s)
#'   (if \code{comb.random} is not \code{NULL} or \code{FALSE}).
#' @param lty.fixed Line type(s) for fixed effect estimate(s).
#' @param lty.random Line type(s) for random effects estimate(s).
#' @param col.fixed Colour of line(s) for fixed effect estimate(s).
#' @param col.random Colour of line(s) for random effects estimate(s).
#' @param nulleffect A logical indicating whether line for null effect
#'   should be added to the plot..
#' @param lwd.nulleffect Width of line for null effect.
#' @param col.nulleffect Colour of line for null effect.
#' @param sm A character string indicating underlying summary measure,
#'   i.e., \code{"RD"}, \code{"RR"}, \code{"OR"}, or \code{"ASD"}.
#' @param weight Either a numeric vector specifying relative sizes of
#'   plotting symbols or a character string indicating which type of
#'   plotting symbols is to be used for individual treatment
#'   estimates. One of missing (see Details), \code{"same"},
#'   \code{"fixed"}, or \code{"random"}, can be abbreviated. Plot
#'   symbols have the same size for all studies or represent study
#'   weights from fixed effect or random effects model.
#' @param studlab A logical indicating whether study labels should be
#'   printed in the graph. A vector with study labels can also be
#'   provided (must be of same length as \code{x$event.e} then).
#' @param cex.studlab Size of study labels.
#' @param pos.studlab Position of study labels, see argument
#'   \code{pos} in \code{\link{text}}.
#' @param label.e Label for experimental group.
#' @param label.c Label for control group.
#' @param \dots Graphical arguments as in \code{par} may also be
#'   passed as arguments.
#'
#' @details
#' A L'Abbé plot is a scatter plot with the risk in the control group
#' on the x-axis and the risk in the experimental group on the y-axis
#' (L'Abbé et al., 1987). It can be used to evaluate heterogeneity in
#' meta-analysis.  Furthermore, this plot can aid to choose a summary
#' measure (odds ratio, risk ratio, risk difference) that will result
#' in more consistent results (Jiménez et al., 1997; Deeks, 2002).
#' 
#' If argument \code{backtransf} is TRUE (default), event
#' probabilities will be printed on x- and y-axis. Otherwise,
#' transformed event probabilities will be printed as defined by the
#' summary measure, i.e., log odds of probabilities for odds ratio as
#' summary measure (\code{sm = "OR"}), log probabilities for \code{sm
#' = "RR"}, and arcsine-transformed probabilities for \code{sm =
#' "ASD"}.
#' 
#' If \code{comb.fixed} is TRUE, the pooled estimate of the fixed
#' effect model is plotted as a line. If \code{comb.random} is TRUE,
#' the pooled estimate of the random effects model is plotted as a
#' line.
#' 
#' Information from object \code{x} is utilised if argument
#' \code{weight} is missing. Weights from the fixed effect model are
#' used (\code{weight = "fixed"}) if argument \code{x$comb.fixed} is
#' \code{TRUE}; weights from the random effects model are used
#' (\code{weight = "random"}) if argument \code{x$comb.random} is
#' \code{TRUE} and \code{x$comb.fixed} is \code{FALSE}.
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{metabin}}
#' 
#' @references
#' Deeks JJ (2002):
#' Issues in the selection of a summary statistic for meta-analysis of
#' clinical trials with binary outcomes.
#' \emph{Statistics in Medicine},
#' \bold{21}, 1575--600
#'
#' Jiménez FJ, Guallar E, Martín-Moreno JM (1997):
#' A graphical display useful for meta-analysis.
#' \emph{European Journal of Public Health},
#' \bold{1}, 101--5
#' 
#' L'Abbé KA, Detsky AS, O'Rourke K (1987):
#' Meta-analysis in clinical research.
#' \emph{Annals of Internal Medicine},
#' \bold{107}, 224--33
#' 
#' @keywords hplot
#' 
#' @examples
#' data(Olkin1995)
#' m1 <- metabin(ev.exp, n.exp, ev.cont, n.cont,
#'               data = Olkin1995,
#'               studlab = paste(author, year),
#'               sm = "RR", method = "I")
#' 
#' # L'Abbe plot for risk ratio
#' #
#' labbe(m1)
#' 
#' # L'Abbe plot for odds ratio
#' #
#' labbe(m1, sm = "OR")
#' # same plot
#' labbe(update(m1, sm = "OR"))
#' 
#' # L'Abbe plot for risk difference
#' #
#' labbe(m1, sm = "RD")
#' 
#' # L'Abbe plot on log odds scale
#' #
#' labbe(m1, sm = "OR", backtransf = FALSE)
#' 
#' # L'Abbe plot for odds ratio with coloured lines for various
#' # treatment effects (defined as log odds ratios)
#' #
#' mycols <- c("blue", "yellow", "green", "red",
#'             "green", "yellow", "blue")
#' labbe(m1, sm = "OR",
#'       comb.random = FALSE,
#'       TE.fixed = log(c(1 / 10, 1 / 5, 1 / 2, 1, 2, 5, 10)),
#'       col.fixed = mycols, lwd.fixed = 2)
#' 
#' # L'Abbe plot on log odds scale with coloured lines for various
#' # treatment effects (defined as log odds ratios)
#' #
#' labbe(m1, sm = "OR",
#'       comb.random = FALSE,
#'       TE.fixed = log(c(1 / 10, 1 / 5, 1 / 2, 1, 2, 5, 10)),
#'       col.fixed = mycols, lwd.fixed = 2,
#'       backtransf = FALSE)
#' 
#' @rdname labbe
#' @method labbe metabin
#' @export
#' @export labbe.metabin


labbe.metabin <- function(x,
                          xlim, ylim,
                          xlab = NULL, ylab = NULL,
                          TE.fixed = x$TE.fixed,
                          TE.random = x$TE.random,
                          comb.fixed = x$comb.fixed,
                          comb.random = x$comb.random,
                          backtransf = x$backtransf,
                          axes = TRUE,
                          pch = 21, text = NULL, cex = 1,
                          col = "black", bg = "lightgray",
                          lwd = 1, lwd.fixed = lwd, lwd.random = lwd,
                          lty.fixed = 2, lty.random = 9,
                          col.fixed = col, col.random = col,
                          nulleffect = TRUE,
                          lwd.nulleffect = lwd, col.nulleffect = "lightgray",
                          sm = x$sm, weight,
                          studlab = FALSE, cex.studlab = 0.8, pos.studlab = 2,
                          label.e = x$label.e, label.c = x$label.c,
                          ...) {
  
  
  ##
  ##
  ## (1) Check for meta object and upgrade older meta objects
  ##
  ##
  chkclass(x, "metabin")
  x <- updateversion(x)
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  chklogical(comb.fixed)
  chklogical(comb.random)
  chklogical(backtransf)
  chklogical(axes)
  chknumeric(cex)
  chknumeric(lwd)
  chknumeric(lwd.fixed)
  chknumeric(lwd.random)
  chknumeric(lty.fixed)
  chknumeric(lty.random)
  chklogical(nulleffect)
  chknumeric(lwd.nulleffect)
  chknull(sm)
  sm <- setchar(sm, .settings$sm4bin)
  chknumeric(cex.studlab)
  pos.studlab <- as.numeric(setchar(pos.studlab, as.character(1:4)))
  
  
  if (!backtransf) {
    xpos <- (x$event.c + x$incr.c) / (x$n.c + 2 * x$incr.c)
    ypos <- (x$event.e + x$incr.e) / (x$n.e + 2 * x$incr.e)
  }
  else {
    xpos <- x$event.c / x$n.c
    ypos <- x$event.e / x$n.e
  }
  ##
  if(length(xpos) != length(ypos))
    stop("event rates must be of same length")
  
  
  if (sm != x$sm) {
    m <- update(x, sm = sm)
    if (missing(TE.fixed))
      TE.fixed <- m$TE.fixed
    if (missing(TE.random))
      TE.random <- m$TE.random
    w.fixed <- m$w.fixed
    w.random <- m$w.random
  }
  else {
    w.fixed <- x$w.fixed
    w.random <- x$w.random
  }
  ##
  ## Exclude studies from meta-analysis
  ##
  if (!is.null(x$exclude)) {
    xpos <- xpos[!x$exclude]
    ypos <- ypos[!x$exclude]
    w.fixed <- w.fixed[!x$exclude]
    w.random <- w.random[!x$exclude]
  }
  
  
  if (!backtransf) {
    if (sm %in% c("OR", "DOR")) {
      xpos <- log(xpos / (1 - xpos))
      ypos <- log(ypos / (1 - ypos))
    }
    else if (sm == "RR") {
      xpos <- log(xpos)
      ypos <- log(ypos)
    }
    else if (sm == "ASD") {
      xpos <- asin(sqrt(xpos))
      ypos <- asin(sqrt(ypos))
    }
  }
  
  
  if (missing(weight))
    weight <- ifelse(comb.random & !comb.fixed, "random", "fixed")
  ##
  iweight <- charmatch(tolower(weight),
                       c("same", "fixed", "random"), nomatch = NA)
  ##
  if(is.na(iweight))
    stop("weight should be \"same\", \"fixed\", or \"random\"")
  ##
  weight <- c("same", "fixed", "random")[iweight]
  ##
  if (weight == "fixed")
    cex.i <- 4 * cex * sqrt(w.fixed) / sqrt(max(w.fixed))
  else if (weight == "random")
    cex.i <- 4 * cex * sqrt(w.random) / sqrt(max(w.random))
  else if (weight == "same")
    cex.i <- cex
  ##
  if (min(cex.i) < 0.5)
    cex.i <- cex.i + (0.5 - min(cex.i))
  
  
  if (is.logical(studlab) && studlab)
    studlab <- x$studlab
  
  
  if (length(comb.fixed) == 0)
    comb.fixed <- TRUE
  if (length(comb.random) == 0)
    comb.random <- TRUE
  
  
  if (!backtransf)
    minval <- min(c(xpos, ypos), na.rm = TRUE)
  else
    minval <- 0
  ##
  if (missing(xlim) & missing(ylim)) {
    xlim <- c(minval, max(c(xpos, ypos), na.rm = TRUE))
    ylim <- xlim
  }
  if (missing(xlim))
    xlim <- c(minval, max(c(xpos, ypos), na.rm = TRUE))
  if (missing(ylim))
    ylim <- xlim
  
  
  oldpar <- par(pty = "s")
  on.exit(par(oldpar))
  
  
  if (is.null(xlab)) {
    if (length(label.c) > 0) {
      xlab <- paste0("Event rate (", label.c, ")")
      if (!backtransf)
        if (sm %in% c("OR", "DOR"))
          xlab <- paste0("ln (odds) ", label.c)
        else if (sm == "RR")
          xlab <- paste0("ln (event rate) ", label.c)
        else if (sm == "ASD")
          xlab <- paste0("Arcsin-transformed event rate (", label.c, ")")
    }
    else {
      xlab <- "Event rate (Control)"
      if (!backtransf)
        if (sm %in% c("OR", "DOR"))
          xlab <- "ln (odds) Control"
        else if (sm == "RR")
          xlab <- "ln (event rate) Control"
        else if (sm == "ASD")
          xlab <- "Arcsin-transformed event rate (Control)"
    }
  }
  ##
  if (is.null(ylab)) {
    if (length(label.e) > 0) {
      ylab <- paste0("Event rate (", label.e, ")")
      if (!backtransf)
        if (sm %in% c("OR", "DOR"))
          ylab <- paste0("ln (odds) ", label.e)
        else if (sm == "RR")
          ylab <- paste0("ln (event rate) ", label.e)
        else if (sm == "ASD")
          ylab <- paste0("Arcsin-transformed event rate (", label.e, ")")
    }
    else {
      ylab <- "Event rate (Experimental)"
      if (!backtransf)
        if (sm %in% c("OR", "DOR"))
          ylab <- "ln (odds) Experimental"
        else if (sm == "RR")
          ylab <- "ln (event rate) Experimental"
        else if (sm == "ASD")
          ylab <- "Arcsin-transformed event rate (Experimental)"
    }
  }
  
  
  if (comb.fixed && length(lty.fixed) == 1 & length(TE.fixed) > 1)
    lty.fixed <- rep(lty.fixed, length(TE.fixed))
  ##
  if (comb.fixed && length(lwd.fixed) == 1 & length(TE.fixed) > 1)
    lwd.fixed <- rep(lwd.fixed, length(TE.fixed))
  ##
  if (comb.fixed && length(col.fixed) == 1 & length(TE.fixed) > 1)
    col.fixed <- rep(col.fixed, length(TE.fixed))
  
  if (comb.random && length(lty.random) == 1 & length(TE.random) > 1)
    lty.random <- rep(lty.random, length(TE.random))
  ##
  if (comb.random && length(lwd.random) == 1 & length(TE.random) > 1)
    lwd.random <- rep(lwd.random, length(TE.random))
  ##
  if (comb.random && length(col.random) == 1 & length(TE.random) > 1)
    col.random <- rep(col.random, length(TE.random))
  
  
  ##
  ## Generate L'Abbe plot
  ##
  plot(xpos, ypos, type = "n", 
       xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
       axes = axes, ...)
  ##
  if (nulleffect)
    abline(0, 1, lwd = lwd.nulleffect, col = col.nulleffect)
  ##
  points(xpos, ypos, pch = pch, cex = cex.i, col = col, bg = bg, lwd = lwd)
  
  
  ##
  ## Auxillary function
  ##
  addlines <- function(x.line, y.line, ylim, lty, lwd, col) {
    sel <- min(ylim) <= y.line & y.line <= max(ylim)
    if (sum(sel) > 1)
      lines(x.line[sel], y.line[sel],
            lty = lty, lwd = lwd, col = col)
  }
  
  
  ##
  ## Add results for common effect model
  ##
  if (comb.fixed & length(TE.fixed) > 0) {
    x.line <- seq(min(xlim), max(xlim), len = 100)
    ##
    if (!backtransf)
      for (i in 1:length(TE.fixed))
        abline(TE.fixed[i], 1,
               lty = lty.fixed[i], lwd = lwd.fixed[i],
               col = col.fixed[i])
    else {
      if (sm == "RR") {
        for (i in 1:length(TE.fixed)) {
          y.line <- x.line * exp(TE.fixed[i])
          addlines(x.line, y.line, ylim,
                   lty.fixed[i], lwd.fixed[i], col.fixed[i])
        }
      }
      else if (sm == "RD") {
        for (i in 1:length(TE.fixed)) {
          y.line <- x.line + TE.fixed[i]
          addlines(x.line, y.line, ylim,
                   lty.fixed[i], lwd.fixed[i], col.fixed[i])
        }
      }
      else if (sm %in% c("OR", "DOR")) {
        for (i in 1:length(TE.fixed)) {
          y.line <- exp(TE.fixed[i]) * (x.line / (1 - x.line)) /
            (1 + exp(TE.fixed[i]) * x.line / (1 - x.line))
          addlines(x.line, y.line, ylim,
                   lty.fixed[i], lwd.fixed[i], col.fixed[i])
        }
      }
      else if (sm == "ASD" & length(TE.fixed) > 0) {
        for (i in 1:length(TE.fixed)) {
          y.line <- sin(asin(sqrt(x.line)) + TE.fixed[i])^2
          addlines(x.line, y.line, ylim,
                   lty.fixed[i], lwd.fixed[i], col.fixed[i])
        }
      }
    }
  }
  
  
  ##
  ## Add results for random effects model
  ##
  if (comb.random & length(TE.random) > 0) {
    x.line <- seq(min(xlim), max(xlim), len = 100)
    ##
    if (!backtransf)
      for (i in 1:length(TE.random))
        abline(TE.random[i], 1,
               lty = lty.random[i], lwd = lwd.random[i],
               col = col.random[i])
    else {
      if (sm == "RR") {
        for (i in 1:length(TE.random)) {
          y.line <- x.line * exp(TE.random[i])
          addlines(x.line, y.line, ylim,
                   lty.random[i], lwd.random[i], col.random[i])
        }
      }
      else if (sm == "RD") {
        for (i in 1:length(TE.random)) {
          y.line <- x.line + TE.random[i]
          addlines(x.line, y.line, ylim,
                   lty.random[i], lwd.random[i], col.random[i])
        }
      }
      else if (sm %in% c("OR", "DOR")) {
        for (i in 1:length(TE.random)) {
          y.line <- exp(TE.random[i]) * (x.line / (1 - x.line)) /
            (1 + exp(TE.random[i]) * x.line / (1 - x.line))
          addlines(x.line, y.line, ylim,
                   lty.random[i], lwd.random[i], col.random[i])
        }
      }
      else if (sm == "ASD" & length(TE.random) > 0) {
        for (i in 1:length(TE.random)) {
          y.line <- sin(asin(sqrt(x.line)) + TE.random[i])^2
          addlines(x.line, y.line, ylim,
                   lty.random[i], lwd.random[i], col.random[i])
        }
      }
    }
  }


  ##
  ## Add study labels
  ##
  if (!is.logical(studlab) && length(studlab) > 0)
    text(xpos, ypos, labels = studlab, pos = pos.studlab, cex = cex.studlab)
  
  
  invisible(NULL)
}





#' @rdname labbe
#' @method labbe default
#' @export
#' @export labbe.default


labbe.default <- function(x, y,
                          xlim, ylim,
                          xlab = NULL, ylab = NULL,
                          TE.fixed = NULL, TE.random = NULL,
                          comb.fixed = !is.null(TE.fixed),
                          comb.random = !is.null(TE.random),
                          backtransf = TRUE,
                          axes = TRUE,
                          pch = 21, text = NULL, cex = 1,
                          col = "black", bg = "lightgray",
                          lwd = 1, lwd.fixed = lwd, lwd.random = lwd,
                          lty.fixed = 2, lty.random = 9,
                          col.fixed = col, col.random = col,
                          nulleffect = TRUE,
                          lwd.nulleffect = lwd, col.nulleffect = "lightgray",
                          sm = "", weight,
                          studlab = FALSE, cex.studlab = 0.8, pos.studlab = 2,
                          label.e = NULL, label.c = NULL,
                          ...) {
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  xpos <- x
  ypos <- y
  ##
  if(length(xpos) != length(ypos))
    stop("arguments 'x' and 'y' must be of same length")
  ##
  chklogical(comb.fixed)
  chklogical(comb.random)
  chklogical(backtransf)
  chklogical(axes)
  chknumeric(cex)
  chknumeric(lwd)
  chknumeric(lwd.fixed)
  chknumeric(lwd.random)
  chknumeric(lty.fixed)
  chknumeric(lty.random)
  chklogical(nulleffect)
  chknumeric(lwd.nulleffect)
  chknull(sm)
  sm <- setchar(sm, .settings$sm4bin)
  chknumeric(cex.studlab)
  pos.studlab <- as.numeric(setchar(pos.studlab, as.character(1:4)))
  
  
  if (!backtransf) {
    if (sm %in% c("OR", "DOR")) {
      xpos <- log(xpos / (1 - xpos))
      ypos <- log(ypos / (1 - ypos))
    }
    else if (sm == "RR") {
      xpos <- log(xpos)
      ypos <- log(ypos)
    }
    else if (sm == "ASD") {
      xpos <- asin(sqrt(xpos))
      ypos <- asin(sqrt(ypos))
    }
  }
  
  
  if (!missing(weight))
    cex.i <- 4 * cex * sqrt(weight) / sqrt(max(weight))
  else
    cex.i <- rep(cex, length(xpos))
  ##
  if (min(cex.i) < 0.5)
    cex.i <- cex.i + (0.5 - min(cex.i))
  
  
  if (backtransf)
    minval <- 0
  else
    minval <- min(c(xpos, ypos), na.rm = TRUE)
  ##
  if (missing(xlim) & missing(ylim)) {
    xlim <- c(minval, max(c(xpos, ypos), na.rm = TRUE))
    ylim <- xlim
  }
  if (missing(xlim))
    xlim <- c(minval, max(c(xpos, ypos), na.rm = TRUE))
  if (missing(ylim))
    ylim <- xlim
  
  
  oldpar <- par(pty = "s")
  on.exit(par(oldpar))
  
  
  if (is.null(xlab)) {
    if (length(label.c) > 0) {
      xlab <- paste0("Event rate (", label.c, ")")
      if (!backtransf)
        if (sm %in% c("OR", "DOR"))
          xlab <- paste0("ln (odds) ", label.c)
        else if (sm == "RR")
          xlab <- paste0("ln (event rate) ", label.c)
        else if (sm == "ASD")
          xlab <- paste0("Arcsin-transformed event rate (", label.c, ")")
    }
    else {
      xlab <- "Event rate (Control)"
      if (!backtransf)
        if (sm %in% c("OR", "DOR"))
          xlab <- "ln (odds) Control"
        else if (sm == "RR")
          xlab <- "ln (event rate) Control"
        else if (sm == "ASD")
          xlab <- "Arcsin-transformed event rate (Control)"
    }
  }
  ##
  if (is.null(ylab)) {
    if (length(label.e) > 0) {
      ylab <- paste0("Event rate (", label.e, ")")
      if (!backtransf)
        if (sm %in% c("OR", "DOR"))
          ylab <- paste0("ln (odds) ", label.e)
        else if (sm == "RR")
          ylab <- paste0("ln (event rate) ", label.e)
        else if (sm == "ASD")
          ylab <- paste0("Arcsin-transformed event rate (", label.e, ")")
    }
    else {
      ylab <- "Event rate (Experimental)"
      if (!backtransf)
        if (sm %in% c("OR", "DOR"))
          ylab <- "ln (odds) Experimental"
        else if (sm == "RR")
          ylab <- "ln (event rate) Experimental"
        else if (sm == "ASD")
          ylab <- "Arcsin-transformed event rate (Experimental)"
    }
  }
  
  
  if (comb.fixed && length(lty.fixed) == 1 & length(TE.fixed) > 1)
    lty.fixed <- rep(lty.fixed, length(TE.fixed))
  ##
  if (comb.fixed && length(lwd.fixed) == 1 & length(TE.fixed) > 1)
    lwd.fixed <- rep(lwd.fixed, length(TE.fixed))
  ##
  if (comb.fixed && length(col.fixed) == 1 & length(TE.fixed) > 1)
    col.fixed <- rep(col.fixed, length(TE.fixed))
  
  if (comb.random && length(lty.random) == 1 & length(TE.random) > 1)
    lty.random <- rep(lty.random, length(TE.random))
  ##
  if (comb.random && length(lwd.random) == 1 & length(TE.random) > 1)
    lwd.random <- rep(lwd.random, length(TE.random))
  ##
  if (comb.random && length(col.random) == 1 & length(TE.random) > 1)
    col.random <- rep(col.random, length(TE.random))
  
  
  ##
  ## Generate L'Abbe plot
  ##
  plot(xpos, ypos, type = "n", 
       xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
       axes = axes, ...)
  ##
  if (nulleffect)
    abline(0, 1, lwd = lwd.nulleffect, col = col.nulleffect)
  ##
  points(xpos, ypos, pch = pch, cex = cex.i, col = col, bg = bg, lwd = lwd)
  
  
  ##
  ## Auxillary function
  ##
  addlines <- function(x.line, y.line, ylim, lty, lwd, col) {
    sel <- min(ylim) <= y.line & y.line <= max(ylim)
    if (sum(sel) > 1)
      lines(x.line[sel], y.line[sel],
            lty = lty, lwd = lwd, col = col)
  }
  
  
  ##
  ## Add results for common effect model
  ##
  if (comb.fixed & length(TE.fixed) > 0) {
    x.line <- seq(min(xlim), max(xlim), len = 100)
    ##
    if (!backtransf)
      for (i in 1:length(TE.fixed))
        abline(TE.fixed[i], 1,
               lty = lty.fixed[i], lwd = lwd.fixed[i],
               col = col.fixed[i])
    else {
      if (sm == "RR") {
        for (i in 1:length(TE.fixed)) {
          y.line <- x.line * exp(TE.fixed[i])
          addlines(x.line, y.line, ylim,
                   lty.fixed[i], lwd.fixed[i], col.fixed[i])
        }
      }
      else if (sm == "RD") {
        for (i in 1:length(TE.fixed)) {
          y.line <- x.line + TE.fixed[i]
          addlines(x.line, y.line, ylim,
                   lty.fixed[i], lwd.fixed[i], col.fixed[i])
        }
      }
      else if (sm %in% c("OR", "DOR")) {
        for (i in 1:length(TE.fixed)) {
          y.line <- exp(TE.fixed[i]) * (x.line / (1 - x.line)) /
            (1 + exp(TE.fixed[i]) * x.line / (1 - x.line))
          addlines(x.line, y.line, ylim,
                   lty.fixed[i], lwd.fixed[i], col.fixed[i])
        }
      }
      else if (sm == "ASD" & length(TE.fixed) > 0) {
        for (i in 1:length(TE.fixed)) {
          y.line <- sin(asin(sqrt(x.line)) + TE.fixed[i])^2
          addlines(x.line, y.line, ylim,
                   lty.fixed[i], lwd.fixed[i], col.fixed[i])
        }
      }
    }
  }
  
  
  ##
  ## Add results for random effects model
  ##
  if (comb.random & length(TE.random) > 0) {
    x.line <- seq(min(xlim), max(xlim), len = 100)
    ##
    if (!backtransf)
      for (i in 1:length(TE.random))
        abline(TE.random[i], 1,
               lty = lty.random[i], lwd = lwd.random[i],
               col = col.random[i])
    else {
      if (sm == "RR") {
        for (i in 1:length(TE.random)) {
          y.line <- x.line * exp(TE.random[i])
          addlines(x.line, y.line, ylim,
                   lty.random[i], lwd.random[i], col.random[i])
        }
      }
      else if (sm == "RD") {
        for (i in 1:length(TE.random)) {
          y.line <- x.line + TE.random[i]
          addlines(x.line, y.line, ylim,
                   lty.random[i], lwd.random[i], col.random[i])
        }
      }
      else if (sm %in% c("OR", "DOR")) {
        for (i in 1:length(TE.random)) {
          y.line <- exp(TE.random[i]) * (x.line / (1 - x.line)) /
            (1 + exp(TE.random[i]) * x.line / (1 - x.line))
          addlines(x.line, y.line, ylim,
                   lty.random[i], lwd.random[i], col.random[i])
        }
      }
      else if (sm == "ASD" & length(TE.random) > 0) {
        for (i in 1:length(TE.random)) {
          y.line <- sin(asin(sqrt(x.line)) + TE.random[i])^2
          addlines(x.line, y.line, ylim,
                   lty.random[i], lwd.random[i], col.random[i])
        }
      }
    }
  }


  ##
  ## Add study labels
  ##
  if (!is.logical(studlab) && length(studlab) > 0)
    text(xpos, ypos, labels = studlab, pos = pos.studlab, cex = cex.studlab)
  
  
  invisible(NULL)
}
