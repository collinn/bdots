
if (requireNamespace("tinytest", quietly = TRUE)) {
  if (length(unclass(packageVersion("bdots"))[[1]]) == 4 | Sys.getenv("R_FORCE_TEST") == 'TRUE') {
    tinytest::test_package("bdots", pattern = "*\\.[rR]$")
  }
}

