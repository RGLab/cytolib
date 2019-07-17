# C++ library for the gated cytometry data

`cytolib` provides the c++ library for users to use and interact with the `GatingSet` (the gated cytometry data structure) at c++ level.

### Installation
The **cytolib** package is installed via `R CMD INSTALL ...`. 

R packages wishing to use the libraries in `cytolib` need to:

- add `cytolib` to **LinkingTo** field in the **DESCRIPTION** file so the compiler knows where to find the headers when the user package is complied
e.g.

```
LinkingTo: cytolib
```
### Usage

See **flowWorkspace** package for the example of using `cytolib`.
