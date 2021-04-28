.onLoad <- function(libname, pkgname){

  path <- system.file("python", package = "scEC")

  reticulate::source_python(paste0(path, "/splitCells.py"))
  reticulate::source_python(paste0(path, "/Heterogeneity.py"))

  reticulate::import("scipy", delay_load = TRUE)
  reticulate::import("numpy", delay_load = TRUE)
}
