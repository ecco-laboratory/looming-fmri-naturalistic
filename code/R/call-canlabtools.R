## helper functions for constructing canlabtools matlab calls used by targets ----

canlabtools_fit_model_connectivity <- function (out_path,
                                                tr_duration,
                                                trs_to_use,
                                                bolds,
                                                confounds,
                                                pred.encoding.roi,
                                                script = matlab_fit_model_connectivity) {
  
  matlab_commands = c(
    assign_variable("out_path", out_path),
    # need to pass this in for the band-pass filter preprocessing
    assign_variable("tr_duration", tr_duration),
    # and this for accessing the right volumes from the 4D niftis using spm syntax
    rvec_to_matlab(trs_to_use, matname = "trs_to_use"),
    # these expect lists with one vector for one subject's targets by run
    rvec_to_matlabcell(bolds, matname = "paths_nifti"),
    rvec_to_matlabcell(confounds, matname = "paths_confounds"),
    assign_variable("seed_yhat_path", pred.encoding.roi),
    call_script(script)
  )
  
  out <- run_matlab_target(matlab_commands, out_path, matlab_path)
  return (out)
}

canlabtools_export_statmap <- function (out_path,
                                        roi,
                                        values,
                                        script = matlab_export_statmap) {
  matlab_commands = c(
    assign_variable("out_path", out_path),
    rvec_to_matlab(values, matname = "zs"),
    call_script(script)
  )
  
  if (!is.null(roi)) {
    matlab_commands <- c(
      matlab_commands,
      # the ROI has already been constrained by the source data
      # this will fail weirdly if you don't make usre this matches up with the prereq values
      rvec_to_matlabcell(roi, matname = "region")
    )
  }
  
  matlab_commands <- c(
    matlab_commands,
    call_script(script)
  )
  
  out <- run_matlab_target(matlab_commands, out_path, matlab_path)
  return (out)
}

canlabtools_pred_encoding_pls <- function (out_path_perf,
                                           out_path_pred,
                                           tr_duration,
                                           bolds_masked,
                                           activations1 = NA,
                                           activations2 = NA,
                                           betas,
                                           script = matlab_pred_pls) {
  
  matlab_commands <- c(
    assign_variable("out_path_perf", out_path_perf),
    assign_variable("out_path_pred", out_path_pred),
    # need to pass this in for the HRF convolution
    assign_variable("tr_duration", tr_duration),
    # this modeling script operates across subjects but takes in data by run
    # so these expect lists with one list-field per subject containing vectors for runs
    # the ROI has already been selected by the choice of paths_masked
    rvec_to_matlabcell(bolds_masked, matname = "paths_masked"),
    assign_variable("path_betas", betas)
  )
  
  # stop if both sets of activations are specified. no reason for this to run if that's the case 
  # bc then it's just the regular predictions that come out of the original model fit
  stopifnot(!(all(is.na(activations1)) & all(is.na(activations2))))
  # is.na is vectorized, including for lists. lol... 
  # I guess we could just test if it's length 1 or not but this seems safe?
  if (!all(is.na(activations1))) {
    matlab_commands <- c(matlab_commands,
                         rvec_to_matlabcell(activations1, matname = "paths_activations_1"))
  } else if (!all(is.na(activations2))) {
    matlab_commands <- c(matlab_commands,
                         rvec_to_matlabcell(activations2, matname = "paths_activations_2"))
  }
  
  matlab_commands <- c(matlab_commands,
                       call_script(script))
  
  out <- run_matlab_target(matlab_commands, 
                           c(out_path_perf, out_path_pred), 
                           matlab_path)
  return (out)
}

canlabtools_fit_encoding_pls <- function (out_path_perf,
                                          out_path_pred,
                                          out_path_betas,
                                          tr_duration,
                                          bolds_masked,
                                          activations1,
                                          # must set as NA instead of NULL because the value comes in from a vector
                                          activations2 = NA,
                                          script = matlab_fit_pls) {
  
  matlab_commands <- c(
    assign_variable("out_path_perf", out_path_perf),
    assign_variable("out_path_pred", out_path_pred),
    assign_variable("out_path_betas", out_path_betas),
    # need to pass this in for the HRF convolution
    assign_variable("tr_duration", tr_duration),
    # this modeling script operates across subjects but takes in data by run
    # so these expect lists with one list-field per subject containing vectors for runs
    # the ROI has already been selected by the choice of paths_masked
    rvec_to_matlabcell(bolds_masked, matname = "paths_masked"),
    rvec_to_matlabcell(activations1, matname = "paths_activations_1")
  )
  
  # is.na is vectorized, including for lists. lol... # I guess we could just test if it's length 1 or not but this seems safe?
  if (!all(is.na(activations2))) {
    matlab_commands <- c(matlab_commands,
                         rvec_to_matlabcell(activations2, matname = "paths_activations_2"))
  }
  
  matlab_commands <- c(matlab_commands,
                       call_script(script))
  
  out <- run_matlab_target(matlab_commands, 
                           c(out_path_perf, out_path_pred, out_path_betas), 
                           matlab_path)
  return (out)
}

canlabtools_mask_fmri_data <- function (out_path,
                                        tr_duration,
                                        trs_to_use,
                                        bolds,
                                        confounds,
                                        roi = "Bstem_SC",
                                        script = matlab_mask_fmri_data) {
  
  matlab_commands = c(
    assign_variable("out_path", out_path),
    # need to pass this in for the band-pass filter preprocessing
    assign_variable("tr_duration", tr_duration),
    # and this for accessing the right volumes from the 4D niftis using spm syntax
    rvec_to_matlab(trs_to_use, matname = "trs_to_use"),
    # these expect lists with one vector for one subject's targets by run
    rvec_to_matlabcell(bolds, matname = "paths_nifti"),
    rvec_to_matlabcell(confounds, matname = "paths_confounds"),
    rvec_to_matlabcell(roi, matname = "region"),
    call_script(script)
  )
  
  out <- run_matlab_target(matlab_commands, out_path, matlab_path)
  return (out)
}
