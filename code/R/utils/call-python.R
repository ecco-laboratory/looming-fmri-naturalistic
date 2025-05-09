## helper functions for constructing python calls used by targets ----

run_python_target <- function (command_args, out_path, conda_path) {
  with_path(conda_path, code = system2("python", args = command_args))
  return (out_path)
}
