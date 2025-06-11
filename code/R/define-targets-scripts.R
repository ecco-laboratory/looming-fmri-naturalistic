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
    name = matlab_select_atlas_subset,
    command = here::here("code", "matlab", "select_this_atlas_subset.m"),
    format = "file"
  ),
  tar_target(
    name = matlab_mask_fmri_data,
    command = {
      matlab_select_atlas_subset
      here::here("code", "matlab", "mask_fmri_data_canlabtools.m")
      },
    format = "file"
  ),
  tar_target(
    name = matlab_parcellate_betas,
    command = {
      matlab_select_atlas_subset
      here::here("code", "matlab", "parcellate_betas_canlabtools.m")
      },
    format = "file"
  ),
  tar_target(
    name = matlab_load_encoding_activations,
    command = here::here("code", "matlab", "load_encoding_activations_allsubs.m"),
    format = "file"
  ),
  tar_target(
    name = matlab_load_bold_for_pls,
    command = here::here("code", "matlab", "load_fmri_data_for_pls_allsubs.m"),
    format = "file"
  ),
  tar_target(
    name = matlab_fit_pls,
    command = {
      matlab_load_encoding_activations
      matlab_load_bold_for_pls
      here::here("code", "matlab", "fit_pls_canlabtools.m")
      },
    format = "file"
  ),
  tar_target(
    name = matlab_fit_pls_no_xval,
    command = {
      matlab_load_encoding_activations
      matlab_load_bold_for_pls
      here::here("code", "matlab", "fit_pls_no_xval_canlabtools.m")
    },
    format = "file"
  ),
  tar_target(
    name = matlab_pred_pls,
    command = {
      matlab_load_encoding_activations
      matlab_load_bold_for_pls
      here::here("code", "matlab", "pred_pls_canlabtools.m")
    },
    format = "file"
  ),
  tar_target(
    name = matlab_fit_model_connectivity,
    command = here::here("code", "matlab", "fit_modelbased_connectivity_canlabtools.m"),
    format = "file"
  ),
  tar_target(
    name = matlab_export_statmap,
    command = here::here("code", "matlab", "export_statmap_canlabtools.m"),
    format = "file"
  ),
  tar_target(
    name = matlab_parcellate_avg,
    command = {
      matlab_select_atlas_subset
      here::here("code", "matlab", "parcellate_avg_fmri_data_canlabtools.m")
      },
    format = "file"
  ),
  tar_target(
    name = matlab_apply_wb_signature,
    command = {
      here::here("code", "matlab", "apply_wb_signature_canlabtools.m")
    },
    format = "file"
  )
)
