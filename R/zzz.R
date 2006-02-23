.First.lib <- function(lib, pkg) library.dynam("kohonen",pkg,lib)
.Last.lib <- function(libpath) library.dynam.unload("kohonen", libpath)
