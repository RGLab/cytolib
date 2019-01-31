switch_cytotool_version <- function(cytoset = TRUE) {
  if(cytoset)
  {
    pkg <- "flowWorkspaceData"
    if(require(pkg, quietly = T, warn.conflicts = F, character.only = T)&&packageVersion(pkg) < "3.17.1")
      devtools::install_github(file.path("RGLab", pkg), ref = "cytoset")
                               
    pkg <- "RProtoBufLib"
    if(require(pkg, quietly = T, warn.conflicts = F, character.only = T)&&packageVersion(pkg) < "2.1.7")
      devtools::install_github(file.path("RGLab", pkg), ref = "cytoset")

    pkg <- "cytolib"
    if(require(pkg, quietly = T, warn.conflicts = F, character.only = T)&&packageVersion(pkg) < "2.1.12")
      devtools::install_github(file.path("RGLab", pkg), ref = "cytoset")

    pkg <- "flowWorkspace"
    if(require(pkg, quietly = T, warn.conflicts = F, character.only = T)&&packageVersion(pkg) < "4.0.3")
      devtools::install_github(file.path("RGLab", pkg), ref = "cytoset")

    pkg <- "openCyto"
    if(require(pkg, quietly = T, warn.conflicts = F, character.only = T)&&packageVersion(pkg) < "2.21.3")
      devtools::install_github(file.path("RGLab", pkg), ref = "cytoset")

    pkg <- "ggcyto"
    if(require(pkg, quietly = T, warn.conflicts = F, character.only = T)&&packageVersion(pkg) < "1.11.4")
      devtools::install_github(file.path("RGLab", pkg), ref = "trunk")

    pkg <- "COMPASS"
    if(require(pkg, quietly = T, warn.conflicts = F, character.only = T)&&packageVersion(pkg) < "1.21.1")
      devtools::install_github(file.path("RGLab", pkg), ref = "trunk")

    pkg <- "CytoML"
    if(require(pkg, quietly = T, warn.conflicts = F, character.only = T)&&packageVersion(pkg) < "2.9.3")
      devtools::install_github(file.path("RGLab", pkg), ref = "cytoset")
    

  }else
  {
    pkg <- "flowWorkspaceData"
    if(require(pkg, quietly = T, warn.conflicts = F, character.only = T)&&packageVersion(pkg) >= "3.17.1")
      devtools::install_github(file.path("RGLab", pkg), ref = "trunk")
    
    pkg <- "RProtoBufLib"
    if(require(pkg, quietly = T, warn.conflicts = F, character.only = T)&&packageVersion(pkg) >= "2.1.7")
      devtools::install_github(file.path("RGLab", pkg), ref = "trunk")
    
    pkg <- "cytolib"
    if(require(pkg, quietly = T, warn.conflicts = F, character.only = T)&&packageVersion(pkg) >= "2.1.12")
      devtools::install_github(file.path("RGLab", pkg), ref = "trunk")
    
    pkg <- "flowWorkspace"
    if(require(pkg, quietly = T, warn.conflicts = F, character.only = T)&&packageVersion(pkg) >= "4.0.3")
      devtools::install_github(file.path("RGLab", pkg), ref = "trunk")
    
    pkg <- "openCyto"
    if(require(pkg, quietly = T, warn.conflicts = F, character.only = T)&&packageVersion(pkg) >= "2.21.3")
      devtools::install_github(file.path("RGLab", pkg), ref = "trunk")
    
    pkg <- "ggcyto"
    if(require(pkg, quietly = T, warn.conflicts = F, character.only = T)&&packageVersion(pkg) < "1.11.4")
      devtools::install_github(file.path("RGLab", pkg), ref = "trunk")
    
    pkg <- "COMPASS"
    if(require(pkg, quietly = T, warn.conflicts = F, character.only = T)&&packageVersion(pkg) < "1.21.1")
      devtools::install_github(file.path("RGLab", pkg), ref = "trunk")
    
    pkg <- "CytoML"
    if(require(pkg, quietly = T, warn.conflicts = F, character.only = T)&&packageVersion(pkg) >= "2.9.3")
      devtools::install_github(file.path("RGLab", pkg), ref = "trunk")
    
  }
  
}
