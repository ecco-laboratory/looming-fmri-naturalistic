# various helper functions for working with BOLD niftis in targets pipeline ----

get_bold_gz <- function (subject, task, run) {
  inject(here::here(!!!path_here_derivatives, subject, "func",
                    paste(subject, task, run, "space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz", sep = "_")))
}

gunzip_bold <- function (bold_gz) {
  system2("gunzip",
          # do not overwrite if already exists unzipped. skips under the hood
          args = c("-k",
                   bold_gz))
  
  return (str_remove(bold_gz, ".gz"))
}
