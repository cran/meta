print.trimfill <- function(object){
  
  if (!inherits(object, "trimfill"))
    stop("Argument 'object' must be an object of class \"trimfill\"")
  
  
  cat("Trim and fill:\n")
  ##
  cat(paste("- number of filled trials: ",
            sum(object$trimfill),
            "\n", sep=""))
  cat(paste(##"- Parameters:",
            ##"\n",
            "- underlying model: ",
            ifelse(object$ma.fixed,
                   "fixed effect model",
                   "random effects model"),
            "\n",
            "- missing trials on ",
            ifelse(object$left, "left ", "right "), "side",
            "\n",
            "- type=", object$type,
            "\n",
            "- number of iterations: ", object$n.iter,
            " (n.iter.max=",
            object$n.iter.max, ")\n", sep=""))
  
  invisible(NULL)
}
