
.onLoad <- function(libname, pkgname){

  path <- system.file("python", package = "scEC")

  splitCells <- NULL
  splitCells <<- reticulate::import_from_path("splitCells", path = path)

}
