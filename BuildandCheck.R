oldwd <- getwd()

findcurvedir <- "/Users/andre/programs/FindCurve"
pkgdir <- paste0(findcurvedir, "/R")

setwd(pkgdir)
curver <- gsub("Version:[ ]*", "", readLines("DESCRIPTION")[grep('Version:.*',readLines("DESCRIPTION"))])
tgzfile <- paste0("/Users/andre/programs/FindCurve/FindCurve_", curver, ".tar.gz")

############################# REFRESHING THE DOCUMENTATION
# Build the man pages
pkgbuild::compile_dll()
devtools::document()

############################# BUILDING THE PACKAGE
#
devtools::build()

############################# INSTALLING THE PACKAGE
# Install the package
detach("package:FindCurve", unload=TRUE)
remove.packages("FindCurve")
install.packages("/Users/andre/programs/FindCurve/FindCurve_0.0.1.tar.gz", repos = NULL)
library(FindCurve)

############################# CHECKING THE PACKAGE
# Check the package locally
# setwd(tempdir()); system(paste("/Users/andre/bin/R CMD check --as-cran", tgzfile)); setwd(pkgdir)

devtools::check()

# Check the package via rhub
rhub::check_for_cran(tgzfile, email = "A.M.deRoos@uva.nl", check_args = "--as-cran", show_status = FALSE)

# Upload the package for checking to https://win-builder.r-project.org/upload.aspx
devtools::check_win_release()
devtools::check_win_devel()

devtools::release()

############################# EXTENSIVE LOCAL TESTING OF THE PACKAGE (TAKES LONG!!!!!!)
#
devtools::test(stop_on_failure = TRUE)

############################# Spell-checking the documentation
#
devtools::spell_check()
