# Copyright 2019 Fred Hutchinson Cancer Research Center
# See the included LICENSE file for details on the licence that is granted to the user of this software.

# Output the CXX flags. These flags are propagated to sourceCpp via the 
# inlineCxxPlugin (defined below) and to packages via a line in Makevars[.win]
# like this:
#
#  PKG_CXXFLAGS += $(shell "${R_HOME}/bin${R_ARCH_BIN}/Rscript.exe" -e "cytolib::CxxFlags()")
#
CxxFlags <- function() {
   cat(cytolibCxxFlags())
}


#' Output the LD flags for building against cytolib. These flags are propagated
#' to sourceCpp via the inlineCxxPlugin (defined below) and to packages 
#' via a line in Makevars[.win] like this:
#'
#'   PKG_LIBS += $(shell "${R_HOME}/bin${R_ARCH_BIN}/Rscript.exe" -e "cytolib::cytolib_LdFlags()")
#' @export
#' @importFrom RProtoBufLib LdFlags
#' @importFrom RcppParallel RcppParallelLibs
cytolib_LdFlags <- function() {
      libDir <- "lib/"
      if (.Platform$OS.type == "windows")
         libDir <- paste(libDir, .Platform$r_arch, "/", sep="")
      cat(asBuildPath(system.file(paste(libDir, "libcytolib.a", sep = ""), package = "cytolib")))
  }




# Helper function to ape the behavior of the R build system
# when providing paths to libraries
asBuildPath <- function(path) {
   if (.Platform$OS.type == "windows") {
      path <- normalizePath(path)
      if (grepl(' ', path, fixed=TRUE))
         path <- utils::shortPathName(path)
      path <- gsub("\\\\", "/", path)
   }
   return(path)
}
