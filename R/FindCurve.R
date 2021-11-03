#'
#' Computes a solution curve determined by a system of non-linear equations
#'
#' \code{FindCurve} computes a solution curve that is determined by a system
#' of non-linear equations. When computing solution curves as a function of one
#' parameter, \code{FindCurve} can detect transcritical bifurcation points
#' (branching points, BP), limit points (LP) in the solution curve and
#' points where an arbitrary function of the solution reaches an extreme value (EXT).
#' The location of these bifurcation points can subsequently be computed as a
#' function of second parameter.
#'
#'   output <- FindCurve(modelname = NULL, curvetype = NULL, startpoint = NULL,
#'                       stepsize = NULL, bounds = NULL, parameters = NULL,
#'                       options = NULL, cflags = NULL,
#'                       clean = FALSE, force = FALSE, debug = FALSE, silent = FALSE)
#'
#' @param   modelname  (string, required)
#' \preformatted{}
#'                Basename of the file with the specification of the system of
#'                equations. The file should have extension '.h'. For example, if
#'                the system of equations is specified in "Cja.h", the name should
#'                be given as "Cja.h" or "Cja".
#'
#' @param   curvetype  (string, required)
#' \preformatted{}
#'                Type of curve to compute: BP, EQ, LP or EXT
#'
#' @param   startpoint (row vector, required)
#' \preformatted{}
#'                The initial point from which to start the continuation of the curve
#'
#' @param   stepsize   (double value, required)
#' \preformatted{}
#'                The target step size along the computed solution curve
#'
#' @param   bounds     (row vector of length 2, 4 or twice the number of unknowns, required)
#' \preformatted{}
#'                The bounds to the region to which to restrict the computation
#'                of the curve. To only put restrictions on the parameter
#'                values, the vector should be of length 2 or 4 in case
#'                of a one or two parameter continuation, respectively.
#'                Otherwise, the vector should be of a length twice the
#'                number of unknowns in the problem.
#'                The format is always [min1 max1 min2 max2....]
#'
#' @param   parameters (row vector, optional, can be left equal to its default NULL)
#' \preformatted{}
#'                Vector of length PARAMETER_NR (set in the model program
#'                file), specifying the values for the model parameters to use in
#'                the computation. Vectors of other lengths, including an empty
#'                vector will be ignored.
#'
#' @param   options    (row vector of strings, optional, can be left equal to its default NULL)
#' \preformatted{}
#'               Vector with pairs of strings, consisting of an option name and a value (for
#'               example c("par1", "1")) or single options (i.e. c("test")).
#'               Possible option names and their values are:
#' \preformatted{}
#'               \verb{"par1", "<index>"}:    Index of the first bifurcation parameter
#' \preformatted{}
#'               \verb{"par2", "<index>"}:    Index of the second bifurcation parameter
#' \preformatted{}
#'               \verb{"EXTfun", "<index>"}:  Index of the element in the ExtraOutput[]
#'                                            vector, for which to test for maximum and
#'                                            minimum values along the curve
#' \preformatted{}
#'               \verb{"EXTpar", "<index>"}:  Index of the parameter, with respect
#'                                            to which to test for maximum and minimum
#'                                            values along the curve
#' \preformatted{}
#'               \verb{"report", "<value>"}:  Interval between consecutive output of
#'                                            computed points to the console ( >= 1).
#'                                            Minimum value of 1 implies output of every
#'                                            point
#' \preformatted{}
#'               \verb{"negative"}:           Allow negative values in the solution point
#' \preformatted{}
#'               \verb{"noBP"}:               Do not check for branching points while
#'                                            computing equilibrium curves
#' \preformatted{}
#'               \verb{"noEXT"}:              Do not check for extremum points while
#'                                            computing equilibrium curves
#' \preformatted{}
#'               \verb{"noLP"}:               Do not check for limit points while
#'                                            computing equilibrium curves
#' \preformatted{}
#'               \verb{"silent"}:             Do not report any error or diagnostic
#'                                            messages to the console
#' \preformatted{}
#'               \verb{"single"}:             Only compute the first point of the
#'                                            solution curve, do not continue the curve
#' \preformatted{}
#'               \verb{"test"}:               Perform only a single integration over
#'                                            reporting the values of the unknowns and
#'                                            the value of the resulting systems of
#'                                            non-liear equations
#'
#' @param   cflags    (string, optional, can be left equal to its default value NULL)
#' \preformatted{}
#'               String containing possible arguments to pass to the compiler
#'               when compiling the C code. This string can be used for conditional
#'               compilation purposes (i.e. cflags = "-D<macro name>=1") to switch
#'               between different versions of equations to solve for. If cflags = NULL,
#'               the compile script looks whether an environment variable CFLAGS is defined
#'               and uses its value instead.
#'
#' @param   clean      (Boolean, optional argument)
#' \preformatted{}
#'               Specify clean = TRUE as argument to remove all the result files
#'               of previous computations before proceeding
#'
#' @param   force      (Boolean, optional argument)
#' \preformatted{}
#'               Specify force = TRUE as argument to force a rebuilding of the model
#'               before the computation
#'
#' @param   debug      (Boolean, optional argument)
#' \preformatted{}
#'               Specify debug = TRUE as argument to compile the model in verbose
#'               mode and with debugging flag set
#'
#' @param   silent     (Boolean, optional argument)
#' \preformatted{}
#'               Specify silent = TRUE as argument to suppress reporting of compilation
#'               commands and results on the console
#'
#' @return  The output is a list containing the following elements:
#' \preformatted{}
#'   \verb{curvepoints}: Matrix with output for all computed points along the curve
#' \preformatted{}
#'   \verb{curvedesc}:   Column vector with strings, summarizing the numerical details
#'                of the computed curve (i.e., initial point, parameter values,
#'                numerical settings used)
#' \preformatted{}
#'   \verb{bifpoints}:   Matrix with the located bifurcation points along the curve
#' \preformatted{}
#'   \verb{biftypes}:    Column vector of strings, containing a description of the
#'                type of bifurcation for each of the located bifurcation points
#'
#' @examples
#' \dontrun{
#' parameters <- c(0.1, 30.0, 1.0, 0.1, 0.5, 0.1, 0.015, 0.5, 3.0, 0.0, 0.0, 0.0, 0.0,
#'                 10.0, 0.0, 0.32, 0.032, 0.5, 0.005, 0.0, 3.0, 0.0, 1.0, 0.0);
#' output1 <- FindCurve("CjaP", "EQ", c(0.0, 2.12019, 0.892170, 11.3152, 0.0), 0.5, c(0.0, 1.0),
#'                      parameters, c("par1", "10", "EXTfun", "1", "EXTpar", "10"),
#'                      clean = TRUE, force = TRUE);
#' }
#'
#' @import pkgbuild
#' @import utils
#' @export
FindCurve <- function(modelname = NULL, curvetype = NULL, startpoint = NULL, stepsize = NULL, bounds = NULL, parameters = NULL, options = NULL, cflags = NULL, clean = FALSE, force = FALSE, debug = FALSE, silent = FALSE) {

  if ((!length(modelname)) || (!nchar(modelname))) stop("You have to specify a basename of the file")

  pkgsrcdir <- system.file("C", package="FindCurve")
  pkgsystemsdir <- system.file("Systems", package="FindCurve")

  if (regexpr("\\.h", modelname) != (nchar(modelname)-1)) modelname <- paste(modelname, ".h", sep='')
  modelfound <- 0
  if (file.exists(paste0("./", modelname))) {
      hfile.fullname <- paste0("./", modelname)
      modelfound <- 1
    }
  else {
      str <- readline(paste0("\nModel file not found in directory ", getwd(), "\nContinue searching for model file in package directory? [y(es)/q(uit)] : "))
      if ((str == '') || (str == 'y') || (str == 'Y')) {
          if (file.exists(paste0(pkgsystemsdir, "/", modelname))) {
              hfile.fullname <- paste0(pkgsystemsdir, "/", modelname)
              modelfound <- 1
            }
        }
    }

  if (modelfound != 1) {
    if (regexpr("\\.h", modelname) == (nchar(modelname)-1))
      stop(paste('No model source file named', modelname, 'can be found', sep=' '))
    else
      stop(paste('No model source file named', paste(modelname, ".h", sep=''), 'can be found', sep=' '))
  }

  if ((length(curvetype) != 1) || (!is.character(curvetype))) stop('Curve type should be one of the strings BP, EQ, LP or EXT')
  if ((!length(startpoint)) || (!is.double(startpoint))) stop('Starting values should be a vector with double values')
  if ((length(stepsize) != 1) || (!is.double(stepsize))) stop('Step size argument should be a single double value')
  if (((length(bounds) %% 2) != 0) || (!is.double(bounds))) stop('Bounds argument should be a vector with double values of length 2, 4 or twice the number of unknowns')
  if ((length(parameters)) && (!is.double(parameters))) stop('If specified parameter values should be a vector with double values')
  if ((length(options)) && (!is.character(options))) stop('If specified options should be an array with strings')

  FCsrcdir.name <- normalizePath(system.file("C", package="FindCurve"))
  tmpdir <- normalizePath(tempdir())

  hfile.dirname=dirname(normalizePath(hfile.fullname))
  hfile.basename=basename(normalizePath(hfile.fullname))
  hfile.abspath=paste0(hfile.dirname, '/', hfile.basename)
  model.name=sub("\\.h", "", hfile.basename)

  libfile.basename <- paste0(model.name, "equi", .Platform$dynlib.ext)
  libfile.fullname = paste0(tmpdir, '/', libfile.basename)

  if (clean) {
    outlist=list.files(pattern=paste0(model.name, "-.*-.*.", "[beo][iru][frt]"))
    if (debug) cat("\nCleaning :", outlist, "\n\n")
    for (i in outlist) file.remove(i)
  }

  oldwd = getwd()
  setwd(hfile.dirname)

  # Check whether the executable is up-to-date
  build.flag <- !file.exists(libfile.fullname);

  if (!build.flag) {
    # The comma in front of mtime is needed to get output as class POSIXct
    build.flag <- (difftime(file.info(hfile.basename)[,'mtime'], file.info(libfile.fullname)[,'mtime'], units = "secs") >= 0)
  }

  if (build.flag || force) {
    # Switch to the temporary directory for building
    setwd(tmpdir)

    srclist <- c("FindCurve.c", "biftest.c", "curve.c", "dopri5.c", "io.c")
    objlist <- c("FindCurve.o", "biftest.o", "curve.o", "dopri5.o", "io.o")

    if (file.exists(libfile.fullname)) file.remove(libfile.fullname)
    if (file.exists(libfile.basename)) file.remove(libfile.basename)
    if (file.exists("FindCurve.o")) file.remove("FindCurve.o");

    if (force) file.remove(Filter(file.exists, objlist))

    if (!silent) cat("\nBuilding executable", libfile.basename, "using sources from", FCsrcdir.name, "...\n\n")

    # Check whether the source files are to be found and copy them to the temporary directory
    for (i in srclist) {
      fname <- paste0(FCsrcdir.name, "/", i)
      if (!file.exists(fname)) {
        setwd(oldwd)
        stop(paste("\nFile:", fname, "not found! Reinstall the FindCurve package\n\n"))
      }
      else file.copy(fname, tmpdir, overwrite = TRUE)
    }

    # Copy the header file to the temporary directory
    file.copy(hfile.abspath, tmpdir, overwrite = TRUE)

    # Construct the command line
    buildargs <- c(paste0("--output=\'", libfile.basename, "\'"), srclist)

    # Define the basic compilation flags
    cppflags <- "-DR_PACKAGE"
    if (!is.null(cflags)) {
      cppflags <- paste0(cppflags, " ", cflags)
    } else {
      if (exists("CFLAGS")) cppflags <- paste0(cppflags, " ", get("CFLAGS"))
    }
    if (debug) cppflags <- paste0(cppflags, " -DDEBUG=1 -g -Wall")
    else cppflags <- paste0(cppflags, " -Wno-format-overflow -Wno-unknown-warning-option")
    cppflags <- paste0(cppflags, " -I.", " -I\"", FCsrcdir.name, "\"")

    # Define the model-specific flags
    modelflags <- paste0(" -DPROBLEMHEADER=\"", hfile.basename, "\"")

    # Define the environment variables to be set
    buildenv <- c(PKG_CFLAGS = paste(cppflags, modelflags), PKG_LIBS = "$(LAPACK_LIBS) $(BLAS_LIBS)")

    if (debug) {
      cat("Current working directory: ", getwd(), "\n")
      cat("Command                  :  R CMD SHLIB\n")
      cat(paste0("Arguments                :  ", paste(buildargs, sep=" ", collapse=" "), "\n"))
      cat("Environment variables    :  \n")
      for (i in 1:length(buildenv)) {
          cat(sprintf("\t%14s = %s\n", names(buildenv)[i], buildenv[i]))
        }
      cat("\n")
    }

    # Compilation steps using the newer R package 'pkgbuild'
    if (silent) sink(tempfile())
    result <- try(pkgbuild::rcmd_build_tools("SHLIB", cmdargs=buildargs, env=buildenv, wd=tmpdir), silent = TRUE)
    if (silent) sink()
    if (is.list(result) && ('status' %in% names(result))) {
      if (result$status != 0) {
        setwd(oldwd)
        cat(result$stdout)
        cat(result$stderr)
        stop(paste0("\nCompilation of ", libfile.basename, " failed!\n"))
      }
      else if (!silent) cat(result$stderr)
    }

    if (!silent) cat("\nCompilation of ", libfile.basename, " succeeded!\n\n")
  } else if (!silent) cat("Dynamic library file ", libfile.fullname, " is up-to-date\n\n");

  setwd(oldwd)

  dyn.load(libfile.fullname)
  if (silent) sink(tempfile())
  if (silent) options <- c(options, "silent")
  cout <- try(.Call("FindCurve", model.name, curvetype, startpoint, stepsize, bounds,
                    parameters, options, PACKAGE=paste0(model.name, "equi")),
              silent = TRUE)
  if (silent) sink()
  dyn.unload(libfile.fullname)

  suspendInterrupts(
    {
      desc = data = bifpoints = biftypes = NULL
      if (exists("cout")) {
        outfile.name = paste0(cout, ".out")
        if (file.exists(outfile.name) && (file.info(outfile.name)$size > 0)) {
          desc <- readLines(outfile.name)
          data <- as.matrix(read.table(text=desc, blank.lines.skip = TRUE, fill=TRUE))
          desc <- desc[grepl("^#", desc)]
          cat(desc, sep='\n')
        }

        biffile.name = paste0(cout, ".bif")
        if (file.exists(biffile.name) && (file.info(biffile.name)$size > 0)) {
          bifinput <- readLines(biffile.name)
          bifpoints <- as.matrix(read.table(text=bifinput, blank.lines.skip = TRUE, comment.char='*', fill=TRUE))
          biftypes = gsub("^.*\\*\\*\\*\\*\\s+|\\s+\\*\\*\\*\\*.*$", "", bifinput)
        }
      }

      setwd(oldwd)
      if (length(desc) || length(data) || length(bifpoints) || length(biftypes)) {
        if (length(bifpoints) || length(biftypes)) {
          output = list(curvedesc = desc, curvepoints = data, bifpoints = bifpoints, biftypes = biftypes)
        }
        else output = list(curvedesc = desc, curvepoints = data)
        return(output)
      } else cat("\nComputations with ", model.name, " produced no output\n")
    }
  )
}

