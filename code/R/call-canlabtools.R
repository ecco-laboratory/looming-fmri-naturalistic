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
    # this modeling script operates across subjects but takes in data by run
    # so these expect lists with one list-field per subject containing vectors for runs
    # the ROI has already been selected by the choice of paths_masked
    rvec_to_matlabcell(roi, matname = "region"),
    call_script(script)
  )
  
  out <- run_matlab_target(matlab_commands, out_path, matlab_path)
  return (out)
}

canlabtools_fit_encoding_pls <- function (out_path_perf,
                                          out_path_pred,
                                          tr_duration,
                                          bolds_masked,
                                          activations, 
                                          script = matlab_fit_pls) {
  
  matlab_commands = c(
    assign_variable("out_path_perf", out_path_perf),
    assign_variable("out_path_pred", out_path_pred),
    # need to pass this in for the HRF convolution
    assign_variable("tr_duration", tr_duration),
    # this modeling script operates across subjects but takes in data by run
    # so these expect lists with one list-field per subject containing vectors for runs
    # the ROI has already been selected by the choice of paths_masked
    rvec_to_matlabcell(bolds_masked, matname = "paths_masked"),
    rvec_to_matlabcell(activations, matname = "paths_activations"),
    call_script(script)
  )
  
  out <- run_matlab_target(matlab_commands, c(out_path_perf, out_path_pred), matlab_path)
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
