#".First.lib" <- function(lib,pkg) {   
#   library.dynam("meta",pkg,lib)
#}

.onLoad <- function(libname,pkgname)
{
   cat("load meta: ", libname, "...\n")
   library.dynam("meta",pkgname,libname)
}

.onUnload <- function(libpath)
{   
   library.dynam.unload("meta",libpath)
}
