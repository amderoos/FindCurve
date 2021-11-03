#' Deletes on request all files produced by the FindCurve package.
#'
#' \code{CleanFindCurve} deletes all FindCurve result files (default) and/or
#' all executables (hit 'F') in the current directory.
#'
#' @param str Character (optional). Only valid argument is 'F'. If not or wrongly
#' specified the user will be asked whether to do a full clean up or whether to
#' quit the clean up.
#'
#' @return None.
#'
#' @examples
#' CleanFindCurve()
#'
#' CleanFindCurve("F")
#'
#' @export
CleanFindCurve <- function(str = NULL) {

  if ((!length(str)) || ((str != "f") && (str != "F"))) {
    cat("\nDelete all FindCurve result files (default) and/or all executables (hit 'F') in the current directory?\n\n")
    str = readline("Press 'Q' to abort, 'F' for a full clean up, any other key to only clean the result files: ")
  }

  if ((str == "q") || (str == "Q")) {
    cat("\nDelete operation aborted\n\n")
  }
  else {
    allsrcs = dir(".", "*\\.h")

    for (i in allsrcs) {
      if ((str == "f") || (str == "F")) {
        cat("\nDeleting all executables in current and temporary directory")
        fname = sub("\\.h", "equi.so", i)
        if (file.exists(fname)) file.remove(fname)

        fname = sub("\\.h", "equi.dll", i)
        if (file.exists(fname)) file.remove(fname)
      }

      fname = sub("\\.h", "", i)
      allres = dir(".", paste0(fname, "-.*-[0-9]*\\.bif"))
      for (j in allres) file.remove(j)
      allres = dir(".", paste0(fname, "-.*-[0-9]*\\.err"))
      for (j in allres) file.remove(j)
      allres = dir(".", paste0(fname, "-.*-[0-9]*\\.out"))
      for (j in allres) file.remove(j)
      allres = dir(".", paste0(fname, "-.*-[0-9]*\\.mat"))
      for (j in allres) file.remove(j)
    }
  }
}
