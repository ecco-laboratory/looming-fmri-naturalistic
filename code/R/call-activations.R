calc_flynet_activations <- function (videos, out_path, script, weights, output_type = "activations") {
  # TODO: May still need to batch this in case there are too many videos to feed into one arg
  video_paths <- paste(videos, collapse = " ")
  
  command_args <- c(script,
                    "-l 132",
                    paste("-i", video_paths),
                    paste("-o", out_path),
                    paste("-w", weights),
                    paste("-q", output_type))
  
  # assume conda_path is an environmental variable defined in the calling script
  out <- run_python_target(command_args, out_path, conda_path)
  return (out)
}

calc_alexnet_activations <- function (videos, out_path, script) {
  # TODO: May still need to batch this in case there are too many videos to feed into one arg
  video_paths <- paste(videos, collapse = " ")
  
  command_args <- c(script,
                    paste("-i", video_paths),
                    paste("-o", out_path))
  
  out <- run_python_target(command_args, out_path, conda_path)
  return (out)
}

resample_video_fps <- function (in_path, out_path, script) {
  with_path(conda_path,
            code = system2("python",
                           args = c(script,
                                    "-i", in_path,
                                    "-o", out_path))
  )
  
  list.files(out_path,
             full.names = TRUE)
}
