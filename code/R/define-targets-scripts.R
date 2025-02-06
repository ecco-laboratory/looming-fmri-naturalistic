# This script is designed to be regular-sourced in the targets pipelines
# to make all of the variables available in the tar_make environment

store_flynet.looming <- "/home/mthieu/Repos/emonet-py/ignore/_targets"
ratings_ck2017 <- tar_read(ratings_ck2017, store = file.path(store_flynet.looming, "subjective"))
weights_flynet <- tar_read(weights_flynet, store = file.path(store_flynet.looming, "subjective"))
py_calc_flynet_activations <- tar_read(py_calc_flynet_activations, store = file.path(store_flynet.looming, "subjective"))
py_make_looming_video <- tar_read(py_make_looming_video, store = file.path(store_flynet.looming, "eyeblink"))

targets_scripts <- list(
  tar_target(
    name = py_get_video_metadata,
    command = here::here("code", "python", "get_video_metadata.py"),
    format = "file"
  ),
  tar_target(
    name = py_resample_video_fps,
    command = here::here("code", "python", "resample_video_fps.py"),
    format = "file"
  ),
  tar_target(
    name = py_calc_alexnet_activations,
    command = here::here("code", "python", "calc_alexnet_activations.py"),
    format = "file"
  ),
  tar_target(
    name = matlab_optimize_ga_trial_order,
    command = here::here("code", "matlab", "optimize_ga_trial_order.m"),
    format = "file"
  ),
  tar_target(
    name = matlab_spmbatch_smooth,
    command = here::here("code", "matlab", "smooth.m"),
    format = "file"
  ),
  tar_target(
    name = matlab_spmbatch_spec_est_level1,
    command = here::here("code", "matlab", "specify_estimate_level1.m"),
    format = "file"
  ),
  tar_target(
    name = matlab_spmbatch_spec_est_level2,
    command = here::here("code", "matlab", "specify_estimate_level2.m"),
    format = "file"
  ),
  tar_target(
    name = matlab_spmbatch_voi,
    command = here::here("code", "matlab", "voi.m"),
    format = "file"
  ),
  tar_target(
    name = matlab_parcellate_betas,
    command = here::here("code", "matlab", "parcellate_betas_canlabtools.m"),
    format = "file"
  ),
  tar_target(
    name = matlab_fit_pls,
    command = here::here("code", "matlab", "fit_pls_canlabtools.m"),
    format = "file"
  )
)
