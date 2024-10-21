## helper functions for constructing matlab calls used by targets ----

assign_variable <- function (var, val, force_unquote = FALSE) {
  stopifnot(length(val) == 1) # do NOT want vectorized behavior
  if (is.character(val) & length(val) == 1 & !force_unquote) val <- wrap_single_quotes(val)
  glue("{var} = {val}")
}

call_function <- function (fn, args) {
  args <- glue_collapse(args, sep = ",")
  glue("{fn}({args})")
}

call_script <- function (script) glue("run({wrap_single_quotes(script)})")

wrap_single_quotes <- function (x) glue("'{x}'")
