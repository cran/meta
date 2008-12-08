addSpace <- function(x, z="", end=TRUE){
  ##
  n.char <- max(c(nchar(x), nchar(z))) - nchar(x)
  ##
  label <- vector()
  for ( i in 1:length(n.char) ){
    label[i] <- paste(rep(" ", n.char[i]), collapse="")
  }
  ##
  if (end)
    res <- paste(x, label, sep="")
  else
    res <- paste(label, x, sep="")
  ##
  res
}
