.onAttach <- function(libname, pkgname) {
  msg <- "\nWelcome to the FindCurve package for computing solution curves to systems of non-linear equations\n"
  msg <- c(msg, "Explore the demos (shown by demo(package=\"FindCurve\") to get an overview\n")
  msg <- c(msg, "Also check out the help pages of the exported functions:\n\n")
  msg <- c(msg, "FindCurve, CleanFindCurve, and showsystem\n\n")

  packageStartupMessage(msg)
}
