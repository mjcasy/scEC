
.onLoad <- function(libname, pkgname){

  path <- system.file("python", package = "scEC")

  PyFunc <- reticulate::import_from_path("PyFunc", path = path)
  assign(x = "PyFunc", value = PyFunc, envir = topenv())

}

