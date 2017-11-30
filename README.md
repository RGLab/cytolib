# C++ library for the gated cytometry data

`cytolib` provides the c++ headers for users to use and interact with the `GatingSet` (the gated cytometry data structure) at c++ level.


The **cytolib** package is installed via `R CMD INSTALL ...`. 

R packages wishing to use the libraries in `cytolib` only need to:

- add `cytolib` to **LinkingTo** field in the **DESCRIPTION** file so the compiler knows where to find the headers when the user package is complied
e.g.

```
LinkingTo: cytolib
```

See **flowWorkspace** package for the example of using `cytolib`.
