#'
#' Shows the system definition file of one of the examples provided with FindCurve
#'
#' \code{showsystem} displays the file contents of one of the systems of equations
#' that is provided as an example.
#'
#'   showsystem(modelname = NULL)
#'
#' @param  modelname  (string)
#' \preformatted{}
#'               Name of the example system to be displayed.
#'
#' @examples
#' showsystem("CjaP.h")
#'
#' @export
showsystem <- function(modelname = NULL) {
  oldwd <- getwd()
  modeldir <- system.file("Systems", package="FindCurve")
  setwd(modeldir)

  if ((!length(modelname)) || (!nchar(modelname)) || (!file.exists(modelname))) {
    cat("\nAvailable example models:\n\n")
    allmodels <- list.files(".", ".h")
    cat(" ")
    cat(paste0(allmodels, sep="\n"))
    cat("\nYou have to specify one of the above model names\n\n")
  }
  else {
    file.show(modelname)
  }

  setwd(oldwd)
}
