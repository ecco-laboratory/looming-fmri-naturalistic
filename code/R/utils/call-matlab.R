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

# wrapper around matlabr::run_matlab_code that returns the path to the output object
# and that always prepends the path to the matlab install to the R search path
# so that the system call works
run_matlab_target <- function (commands, out_path, matlab_path) {
  matlab_exit_status <- with_path(matlab_path, run_matlab_code(commands))
  stopifnot(matlab_exit_status == 0)
  
  if (dir.exists(out_path)) {
    return (list.files(out_path, full.names = TRUE))
  } else {
    return (out_path)
  }
}
