## helper functions for constructing canlabtools matlab calls used by targets ----

canlabtools_fit_encoding_pls <- function (out_paths,
                                          tr_duration,
                                          trs_to_use,
                                          bolds,
                                          activations, 
                                          confounds,
                                          roi = "Bstem_SC",
                                          script = matlab_fit_pls) {
  
  matlab_commands = c(
    rvec_to_matlabcell(out_paths, matname = "out_paths"),
    # need to pass this in for the HRF convolution
    assign_variable("tr_duration", tr_duration),
    # and this for accessing the right volumes from the 4D niftis using spm syntax
    rvec_to_matlab(trs_to_use, matname = "trs_to_use"),
    # this modeling script operates across subjects but takes in data by run
    # so these expect lists with one list-field per subject containing vectors for runs
    rvec_to_matlabcell(bolds, matname = "paths_nifti"),
    rvec_to_matlabcell(activations, matname = "paths_activations"),
    rvec_to_matlabcell(confounds, matname = "paths_confounds"),
    rvec_to_matlabcell(roi, matname = "regions"),
    # TODO: Add confounds back in when the script takes them for lightweight canlabtools preproc
    call_script(script)
  )
  
  with_path(
    matlab_path,
    run_matlab_code(matlab_commands)
  )
  
  return (unlist(out_paths))
}
