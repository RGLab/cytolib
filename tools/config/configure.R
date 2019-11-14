# Prepare your package for installation here.
# Use 'define()' to define configuration variables.
# Use 'configure_file()' to substitute configuration values.

define(CYTOLIB_VERSION = read.dcf("DESCRIPTION")[1,][["Version"]])

configure_file("cytolibConfig.h.in", "inst/include/cytolib/cytolibConfig.h")
