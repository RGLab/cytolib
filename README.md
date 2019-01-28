# C++ library for the gated cytometry data

`cytolib` provides the c++ headers for users to use and interact with the `GatingSet` (the gated cytometry data structure) at c++ level.

### Installation
The **cytolib** package is installed via `R CMD INSTALL ...`. 

R packages wishing to use the libraries in `cytolib` only need to:

- add `cytolib` to **LinkingTo** field in the **DESCRIPTION** file so the compiler knows where to find the headers when the user package is complied
e.g.

```
LinkingTo: cytolib
```
### Usage
Make sure to call `CYTOLIB_INIT()` in user c code

See **flowWorkspace** package for the example of using `cytolib`.

### Switching between the old and new tool chains
```
#new
devtools::install_github(c("RGLab/RProtoBufLib", "RGLab/cytolib", "RGLab/flowWorkspace", "RGLab/openCyto","RGLab/CytoML"), ref = "cytoset")
#old
devtools::install_github(c("RGLab/RProtoBufLib", "RGLab/cytolib", "RGLab/flowWorkspace", "RGLab/openCyto","RGLab/CytoML"), ref = "trunk")
```
