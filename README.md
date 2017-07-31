# C++ library for the gated cytometry data

`cytolib` provides the c++ headers and library file for users to use and interact with the `GatingSet` (the gated cytometry data structure) at c++ level.


The **cytolib** package is installed in the normal `R` manner without the need of any user efforts.

All packages wishing to use the libraries in `cytolib` only need to:

- add `cytolib` to **LinkingTo** field in **DESCRIPTION** file so that the compiler knows where to find the headers when user package is complied
e.g.
```
LinkingTo: cytolib
```

- set **PKG_LIBS** in **src/Makevars** file so that linker can find and linked to the **libcytolib.a** file 
e.g.
```bash
PKG_LIBS =`${R_HOME}/bin/Rscript -e "cytolib:::LdFlags()"`
```

See **flowWorkspace** package for the example of using `cytolib`.
