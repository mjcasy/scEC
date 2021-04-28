
.onLoad <- function(libname, pkgname){

  path <- system.file("python", package = "scEC")

  PyFunc <- NULL
  PyFunc <<- reticulate::import_from_path("PyFunc", path = path)

}
