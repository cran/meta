## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = TRUE,
  warning = FALSE
  )
options(knitr.kable.NA = ".")

load("mmiss_limit.rda")

## ----eval = FALSE-------------------------------------------------------------
#  install.packages(c("meta", "metasens"))

## -----------------------------------------------------------------------------
library(meta)

## ----eval = FALSE-------------------------------------------------------------
#  library(metasens)

## -----------------------------------------------------------------------------
settings.meta(digits = 2, method.tau = "PM")

## -----------------------------------------------------------------------------
joy = read.csv("Joy2006.txt")
# Add new variable: miss
joy$miss = ifelse((joy$drop.h + joy$drop.p) == 0, 
  "Without missing data", "With missing data")
head(joy)
str(joy)

## -----------------------------------------------------------------------------
m.publ = metabin(resp.h, resp.h + fail.h, resp.p, resp.p + fail.p,
  data = joy, studlab = paste0(author, " (", year, ")"),
  label.e = "Haloperidol", label.c = "Placebo",
  label.left = "Favours placebo", label.right = "Favours haloperidol")

## -----------------------------------------------------------------------------
summary(m.publ)

## ----eval = FALSE-------------------------------------------------------------
#  print(summary(m.publ))

## ----eval = FALSE-------------------------------------------------------------
#  forest(m.publ, sortvar = year, prediction = TRUE,
#    file = "figure2.pdf", width = 10)

## ----echo = FALSE, out.width = "95%"------------------------------------------
knitr::include_graphics("figure2.pdf")

## -----------------------------------------------------------------------------
m.publ.sub = update(m.publ, subgroup = miss, print.subgroup.name = FALSE)
m.publ.sub

## ----eval = FALSE-------------------------------------------------------------
#  forest(m.publ.sub, sortvar = year,
#    xlim = c(0.1, 100), at = c(0.1, 0.3, 1, 3, 10, 30, 100),
#    test.subgroup.common = FALSE,
#    label.test.subgroup.random = "Test for subgroup differences:",
#    file = "figure3.pdf", width = 10)

## ----echo = FALSE, out.width = "95%"------------------------------------------
knitr::include_graphics("figure3.pdf")

## ----eval = FALSE-------------------------------------------------------------
#  # Impute as no events (ICA-0) - default
#  mmiss.0 = metamiss(m.publ, drop.h, drop.p)
#  # Impute as events (ICA-1)
#  mmiss.1 = metamiss(m.publ, drop.h, drop.p, method = "1")
#  # Observed risk in control group (ICA-pc)
#  mmiss.pc = metamiss(m.publ, drop.h, drop.p, method = "pc")
#  # Observed risk in experimental group (ICA-pe)
#  mmiss.pe = metamiss(m.publ, drop.h, drop.p, method = "pe")
#  # Observed group-specific risks (ICA-p)
#  mmiss.p = metamiss(m.publ, drop.h, drop.p, method = "p")
#  # Best-case scenario (ICA-b)
#  mmiss.b = metamiss(m.publ, drop.h, drop.p, method = "b", small.values = "bad")
#  # Worst-case scenario (ICA-w)
#  mmiss.w = metamiss(m.publ, drop.h, drop.p, method = "w", small.values = "bad")
#  # Gamble-Hollis method
#  mmiss.gh = metamiss(m.publ, drop.h, drop.p, method = "GH")
#  # IMOR.e = 2 and IMOR.c = 2 (same as available case analysis)
#  mmiss.imor2 = metamiss(m.publ, drop.h, drop.p, method = "IMOR", IMOR.e = 2)
#  # IMOR.e = 0.5 and IMOR.c = 0.5
#  mmiss.imor0.5 = metamiss(m.publ, drop.h, drop.p, method = "IMOR", IMOR.e = 0.5)

## -----------------------------------------------------------------------------
meths = c("Available case analysis (ACA)",
  "Impute no events (ICA-0)", "Impute events (ICA-1)",
  "Observed risk in control group (ICA-pc)",
  "Observed risk in experimental group (ICA-pe)",
  "Observed group-specific risks (ICA-p)",
  "Best-case scenario (ICA-b)", "Worst-case scenario (ICA-w)",
  "Gamble-Hollis analysis",
  "IMOR.e = 2, IMOR.c = 2", "IMOR.e = 0.5, IMOR.c = 0.5")
# Use inverse-variance method for pooling (which is used for
# imputation methods)
m.publ.iv = update(m.publ, method = "Inverse")
# Combine results (random effects)
mbr = metabind(m.publ.iv,
  mmiss.0, mmiss.1,
  mmiss.pc, mmiss.pe, mmiss.p,
  mmiss.b, mmiss.w, mmiss.gh,
  mmiss.imor2, mmiss.imor0.5,
  name = meths, pooled = "random")

## ----eval = FALSE-------------------------------------------------------------
#  forest(mbr, xlim = c(0.5, 4),
#    leftcols = c("studlab", "I2.w", "tau2.w", "Q.w", "pval.Q.w"),
#    leftlab = c("Meta-Analysis Method", "I2", "Tau2", "Q", "P-value"),
#    type.study = "diamond",
#    digits.addcols = c(4, 2, 2, 2), just.addcols = "right",
#    file = "figure4.pdf", width = 10)

## ----echo = FALSE, out.width = "95%"------------------------------------------
knitr::include_graphics("figure4.pdf")

## ----eval = FALSE-------------------------------------------------------------
#  funnel(m.publ)

## -----------------------------------------------------------------------------
metabias(m.publ, method.bias = "score")

## -----------------------------------------------------------------------------
tf.publ = trimfill(m.publ)
tf.publ

## -----------------------------------------------------------------------------
summary(tf.publ)

## ----eval = FALSE-------------------------------------------------------------
#  funnel(tf.publ)

## ----eval = FALSE-------------------------------------------------------------
#  l1.publ = limitmeta(m.publ)

## ----eval = FALSE-------------------------------------------------------------
#  l1.publ

## ----eval = FALSE-------------------------------------------------------------
#  pdf("figure5.pdf", width = 10, height = 10)
#  #
#  par(mfrow = c(2, 2), pty = "s",
#      oma = c(0, 0, 0, 0), mar = c(4.1, 3.1, 2.1, 1.1))
#  #
#  funnel(m.publ, xlim = c(0.05, 50), axes = FALSE)
#  axis(1, at = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 50))
#  axis(2, at = c(0, 0.5, 1, 1.5))
#  box()
#  title(main = "Panel A: Funnel plot", adj = 0)
#  #
#  funnel(m.publ, xlim = c(0.05, 50), axes = FALSE,
#         contour.levels = c(0.9, 0.95, 0.99),
#         col.contour = c("darkgray", "gray", "lightgray"))
#  legend("topright",
#         c("p < 1%", "1% < p < 5%", "5% < p < 10%", "p > 10%"),
#         fill = c("lightgray", "gray", "darkgray", "white"),
#         border = "white", bg = "white")
#  axis(1, at = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 50))
#  axis(2, at = c(0, 0.5, 1, 1.5))
#  box()
#  title(main = "Panel B: Contour-enhanced funnel plot", adj = 0)
#  #
#  funnel(tf.publ, xlim = c(0.05, 50), axes = FALSE)
#  axis(1, at = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 50))
#  axis(2, at = c(0, 0.5, 1, 1.5))
#  box()
#  title(main = "Panel C: Trim-and-fill method", adj = 0)
#  #
#  funnel(l1.publ, xlim = c(0.05, 50), axes = FALSE,
#         col.line = 8, lwd.line = 3)
#  axis(1, at = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 50))
#  axis(2, at = c(0, 0.5, 1, 1.5))
#  box()
#  title(main = "Panel D: Limit meta-analysis", adj = 0)
#  #
#  dev.off()

## ----echo = FALSE, out.width = "95%"------------------------------------------
knitr::include_graphics("figure5.pdf")

