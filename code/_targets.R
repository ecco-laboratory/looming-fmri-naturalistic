# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(tidyverse)
library(rlang)

# Set target options:
tar_option_set(
  packages = c("withr",
               "matlabr",
               "tidyverse",
               "magrittr",
               "glue",
               "rlang",
               "qualtRics"), # packages that your targets need to run
  # format = "qs", # Optionally set the default storage format. qs is fast.
  #
  # For distributed computing in tar_make(), supply a {crew} controller
  # as discussed at https://books.ropensci.org/targets/crew.html.
  # Choose a controller that suits your needs. For example, the following
  # sets a controller with 2 workers which will run as local R processes:
  #
controller = crew::crew_controller_local(workers = 4,
                                         seconds_idle = 30)
  #
  # Alternatively, if you want workers to run on a high-performance computing
  # cluster, select a controller from the {crew.cluster} package. The following
  # example is a controller for Sun Grid Engine (SGE).
  # 
  #   controller = crew.cluster::crew_controller_sge(
  #     workers = 50,
  #     # Many clusters install R as an environment module, and you can load it
  #     # with the script_lines argument. To select a specific verison of R,
  #     # you may need to include a version string, e.g. "module load R/4.3.0".
  #     # Check with your system administrator if you are unsure.
  #     script_lines = "module load R"
  #   )
  #
  # Set other options as needed.
)

# tar_make_clustermq() is an older (pre-{crew}) way to do distributed computing
# in {targets}, and its configuration for your machine is below.
options(clustermq.scheduler = "slurm")
options(clustermq.template = "clustermq.tmpl")

# tar_make_future() is an older (pre-{crew}) way to do distributed computing
# in {targets}, and its configuration for your machine is below.
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Run the R scripts in the R/ folder with your custom functions:
tar_source(c("code/R/"))
# source("other_functions.R") # Source other scripts as needed.

## other useful global vars ----
# My dear little conda env
conda_path <- "/home/mthieu/Repos/emonet-py/env/bin"
matlab_path <- "/opt/MATLAB/R2024a/bin"
# fMRI parameters that pass into multiple targets
tr_length <- 2 # seconds
nback_rating_duration <- 12 # tr units
nback_blocks_per_run <- 4
nback_block_threatenings <- c("dog", "frog", "spider", "above", "below")

# target definition

## targets: scripts (sometimes in other project folders) ----
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

## targets: looming stimuli of various kinds ----

folder_videos_naturalistic <- here::here("ignore", "stimuli", "videos", "naturalistic")

targets_stimuli <- list(
  # this will really only be used to generate placeholder activations
  # to patch into a full encoding model activation timecourse for the naturalistic task
  # I made it by running a fixation screenshot in iMovie for 1 second like a clown
  tar_target(
    video_fixation,
    command = here::here("ignore", "stimuli", "videos", "fixation.mp4"),
    format = "file"
  ),
  tar_target(
    name = images_controlled,
    command = list.files(here::here("ignore", "stimuli", "images", "use"),
                         full.names = TRUE),
    format = "file"
  ),
  tar_target(
    name = videos_controlled,
    command = with_path(conda_path, {
      
      # local ignore/stimuli/videos/controlled is symlinked to this folder as well
      # but write to the server-wide folder to share out permissions
      out_path <- here::here("ignore", "stimuli", "videos", "controlled")
      
      # generate _almost_ every crossing combination
      # but we only need reversed middle-near looms, not reversed middle-far recedes 
      conditions <- crossing(img = images_controlled,
                             angle = c("above", "below"),
                             nesting(direction = c("looming", "looming", "receding"),
                                     reverse = c("", "--reverse", "")))
      
      for (i in 1:nrow(conditions)) {
        cat("starting image:", crayon::green(conditions$img[i]), fill = TRUE)
        cat("direction:", crayon::magenta(conditions$direction[i]),
            "approach angle:", crayon::magenta(conditions$angle[i]),
            crayon::red(conditions$reverse[i]),
            fill = TRUE)
        
        system2("python", args = c(py_make_looming_video,
                                   paste("--infile", conditions$img[i]),
                                   paste("--outpath", out_path),
                                   # Adult Ps sat 40 cm away from monitor
                                   # and I think it needs to come 4 obj-widths away
                                   "--diststart 15",
                                   "--distend 0.1",
                                   "--posstart middle",
                                   paste("--direction", conditions$direction[i]),
                                   conditions$reverse[i],
                                   paste("--angle", conditions$angle[i]),
                                   "--loomtime 0.5"
        )) 
      }
      # print the path to the output videos
      list.files(out_path,
                 pattern = "mp4",
                 full.names = TRUE)
    }),
    format = "file"
  ),
  tar_target(
    # these don't really need to get shown to humans, we can show them the native fps ones
    # this is for the encoding models
    name = videos_controlled_10fps,
    command = {
      videos_controlled
      
      out_path <- here::here("ignore", "stimuli", "videos", "controlled_10fps")
      
      with_path(conda_path,
                code = system2("python",
                               args = c(py_resample_video_fps,
                                        "-i", here::here("ignore", "stimuli", "videos", "controlled"),
                                        "-o", out_path))
      )
      
      list.files(out_path,
                 full.names = TRUE)
    },
    format = "file"
  ),
  tar_target(
    name = qualtrics.ids_videos_controlled,
    command = here::here("ignore", "stimuli", "videos", "controlled_qualtrics_ids.csv"),
    format = "file"
  ),
  tar_target(
    name = videos_naturalistic,
    command = list.files(folder_videos_naturalistic, pattern = ".mp4", full.names = TRUE),
    format = "file"
  ),
  tar_target(
    # these don't really need to get shown to humans, we can show them the native fps ones
    # this is for the encoding models
    name = videos_naturalistic_10fps,
    command = {
      videos_naturalistic
      
      out_path <- here::here("ignore", "stimuli", "videos", "naturalistic_10fps")
      
      with_path(conda_path,
                code = system2("python",
                               args = c(py_resample_video_fps,
                                        "-i", folder_videos_naturalistic,
                                        "-o", out_path))
      )
      
      list.files(out_path,
                 full.names = TRUE)
    },
    format = "file"
  ),
  tar_target(
    name = metadata_videos_controlled,
    command = tibble(filename = videos_controlled) %>% 
      separate_wider_delim(filename, 
                           delim = "/", 
                           names = c(rep(NA, 9), "filename")) %>% 
      # don't drop the filename (without folders)
      separate_wider_delim(filename, 
                           delim = "_", 
                           names = c("image_num", 
                                     "direction", 
                                     NA, 
                                     "hemifield", 
                                     NA, 
                                     NA, 
                                     "reversed"), 
                           too_few = "align_start", 
                           cols_remove = FALSE) %>% 
      # extract the start location information from reversed
      # to distinguish looms, recedes from near, and recedes from middle
      mutate(start_location = if_else(!is.na(reversed), "near", "middle"), 
             direction = if_else(!is.na(reversed), "receding", direction)) %>% 
      select(-reversed) %>% 
      # need image_num temporarily bc it makes labeling animal type easier
      mutate(image_num = as.integer(image_num),
             animal_type = case_when(image_num <= 9 ~ "dog",
                                     image_num <= 29 ~ "spider",
                                     image_num <= 35 ~ "frog",
                                     TRUE ~ NA_character_)) %>% 
      nest(rows = -c(image_num, animal_type)) %>% 
      group_by(animal_type) %>% 
      # animal_num is for binding onto the generated stimlists
      mutate(animal_num = 1:n()) %>% 
      select(-image_num) %>% 
      ungroup() %>% 
      unnest(rows)
  ),
  # this one is for the ones going into the main nback runs
  # aka dropping the reverse-loom-from-near stimuli
  # that makes every other target so much easier to leave as is
  tar_target(
    metadata_videos_nback,
    command = metadata_videos_controlled %>% 
      filter(start_location == "middle") %>% 
      select(-start_location) %>% 
      # condition_num is for binding onto the trial orders from matlab
      nest(rows = -c(animal_type, hemifield, direction)) %>% 
      arrange(animal_type, hemifield, direction) %>% 
      mutate(condition_num = 1:n()) %>% 
      unnest(rows)
  ),
  tar_target(
    metadata_videos_luminance,
    command = metadata_videos_controlled %>% 
      filter(!(direction == "receding" & start_location == "middle")) %>% 
      select(-start_location) %>% 
      # condition_num is for binding onto the trial orders from matlab
      nest(rows = -c(animal_type, hemifield, direction)) %>% 
      arrange(animal_type, hemifield, direction) %>% 
      mutate(condition_num = 1:n()) %>% 
      unnest(rows)
  ),
  tar_target(
    name = metadata_videos_naturalistic,
    command = {
      out_path <- file.path(folder_videos_naturalistic, "metadata.csv")
      
      video_paths <- paste(videos_naturalistic, collapse = " ")
      
      with_path(conda_path,
                code = system2("python",
                               args = c(py_get_video_metadata,
                                        paste("-i", video_paths),
                                        paste("-o", out_path))))
      out_path
    },
    format = "file"
  ),
  tar_target(
    name = annotations_videos_naturalistic,
    command = file.path(folder_videos_naturalistic, "annotations.csv"),
    format = "file"
  ),
  tar_target(
    name = ck2017_imagenet_categories,
    command = here::here("ignore", "data", "norm", "ck2017_imagenet_categories.csv"),
    format = "file"
  )
)

## targets: fMRI stimlist materials for n-back task ----

targets_stimlists_nback <- list(
  tar_target(
    name = stims_nback_base,
    # threatening category is now handled later
    # this just expands to have sufficient stims for each of these categories
    # although! mind that attended is actually set in the psychopy task
    # this just has the effect of doubling the rows to have the desired count in each attended category
    command = expand_grid(attended = c("animal_type", "hemifield"),
                          animal_type = c("dog", "frog", "spider"),
                          hemifield = c("above", "below"))
  ),
  # do a separate stimlist dataframe for each block
  # because each block will get its own PsychoPy run
  # so the between-block counterbalancing order will be handled over there
  tar_map(
    # values are what directly go into the target names
    # so sprintf them here and extract numeric again later /shrug
    values = tidyr::expand_grid(block_num = sprintf("block%02d", 1:5),
                                order_num = sprintf("order%02d", 1:2)),
    tar_target(
      name = stims_nback,
      command = {
        block_num_numeric <- as.integer(str_sub(block_num, start = -2L))
        this_threatening <- nback_block_threatenings[block_num_numeric]
        if (this_threatening %in% c("dog", "frog", "spider")) {
          this_predictable <- "animal_type"
        } else {
          this_predictable <- "hemifield"
        }
        
        stims_nback_base %>% 
          make_stims_nback(threatening = this_threatening,
                           predictable = this_predictable)
      }
    ),
    tar_target(
      name = trial.order_nback,
      command = {
        out_path <- here::here("ignore",
                               "stimlists",
                               "nback",
                               paste0("trialorder_", block_num, "_", order_num, ".txt"))
        n_conditions <- stims_nback %>% 
          distinct(animal_type, hemifield, direction) %>% 
          nrow()
        freq_conditions <- stims_nback %>% 
          count(animal_type, hemifield, direction) %>% 
          pull(n)
        
        t_contrasts <- stims_nback %>% 
          distinct(animal_type, hemifield, direction) %>% 
          mutate(contrast_animal = if_else(animal_type == "dog", 1, -0.5), 
                 contrast_hemifield = if_else(hemifield == "above", 1, -1), 
                 contrast_direction = if_else(direction == "looming", 1, -1)) %>% 
          select(starts_with("contrast")) %>% 
          as.matrix() %>% 
          # so contrasts go along rows
          t()
        
        matlab_commands = c(
          assign_variable("nConditions", n_conditions),
          assign_variable("freqConditions", rmat_to_matlab_mat(t(freq_conditions))),
          assign_variable("isi", 2),
          assign_variable("itiPerTrial", 2),
          assign_variable("ratingRestLength", nback_rating_duration),
          assign_variable("nRatingRests", nback_blocks_per_run),
          assign_variable("tContrasts", rmat_to_matlab_mat(t_contrasts)),
          assign_variable("fitnessWeights", rmat_to_matlab_mat(t(c(0, 1, 0, 1)))),
          call_script(matlab_optimize_ga_trial_order),
          call_function("writematrix", 
                        args = c("Models",
                                 wrap_single_quotes(out_path)))
        )
        
        with_path(
          matlab_path,
          run_matlab_code(matlab_commands)
        )
        
        out_path
      },
      format = "file"
    ),
    tar_target(
      name = stimlist_nback,
      command = {
        out_path_glue <- here::here("ignore",
                                    "stimlists",
                                    "nback",
                                    paste0("stimlist_", block_num, "_", order_num,
                                           "_{sprintf('miniblock%02d', miniblock_num)}.csv"))
        
        block_num_numeric <- as.integer(str_sub(block_num, start = -2L))
        this_threatening <- nback_block_threatenings[block_num_numeric]
        if (this_threatening %in% c("dog", "frog", "spider")) {
          this_predictable <- "animal_type"
        } else {
          this_predictable <- "hemifield"
        }
        
        # output to nowhere so that the file paths coming out of walk(write_csv)
        # will become the file path outputs
        make_stimlist_nback(
          out_path_glue = out_path_glue,
          trial.order_path = trial.order_nback,
          metadata = metadata_videos_nback,
          threatening = this_threatening,
          predictable = this_predictable
        )
        
      },
      format = "file"
    )
  ),
  tar_target(
    name = stimlist_nback_miniblock_order,
    command = {
      out_path <- here::here("ignore", "stimlists", "nback", "miniblock_order.csv")
      tibble(miniblock_num = 1:nback_blocks_per_run,
             miniblock_file = sprintf("miniblock%02d.csv", miniblock_num)) %>% 
        write_csv(file = out_path)
      
      out_path
      },
    format = "file"
  ),
  tar_target(
    name = stims_luminance,
    command = stims_nback_base %>% 
      make_stims_nback()
  ),
  tar_map(
    values = tibble(order_num = sprintf("order%02d", 1:2)),
    tar_target(
      name = trial.order_luminance,
      command = {
        out_path <- here::here("ignore",
                               "stimlists",
                               "nback",
                               paste0("trialorder_luminance_", order_num, ".txt"))
        n_conditions <- stims_luminance %>% 
          distinct(animal_type, hemifield, direction) %>% 
          nrow()
        freq_conditions <- stims_luminance %>% 
          count(animal_type, hemifield, direction) %>% 
          pull(n)
        
        t_contrasts <- stims_luminance %>% 
          distinct(animal_type, hemifield, direction) %>% 
          mutate(contrast_animal = if_else(animal_type == "dog", 1, -0.5), 
                 contrast_hemifield = if_else(hemifield == "above", 1, -1), 
                 contrast_direction = if_else(direction == "looming", 1, -1)) %>% 
          select(starts_with("contrast")) %>% 
          as.matrix() %>% 
          # so contrasts go along rows
          t()
        
        matlab_commands = c(
          assign_variable("nConditions", n_conditions),
          assign_variable("freqConditions", rmat_to_matlab_mat(t(freq_conditions))),
          assign_variable("isi", 2),
          assign_variable("itiPerTrial", 2),
          assign_variable("ratingRestLength", nback_rating_duration),
          assign_variable("nRatingRests", nback_blocks_per_run),
          assign_variable("tContrasts", rmat_to_matlab_mat(t_contrasts)),
          assign_variable("fitnessWeights", rmat_to_matlab_mat(t(c(0, 1, 0, 1)))),
          call_script(matlab_optimize_ga_trial_order),
          call_function("writematrix", 
                        args = c("Models",
                                 wrap_single_quotes(out_path)))
        )
        
        with_path(
          matlab_path,
          run_matlab_code(matlab_commands)
        )
        
        out_path
      }
    ),
    tar_target(
      stimlist_luminance,
      command = {
        out_path_glue <- here::here("ignore",
                                    "stimlists",
                                    "nback",
                                    paste0("stimlist_luminance_", order_num,
                                           "_{sprintf('miniblock%02d', miniblock_num)}.csv"))
        
        # output to nowhere so that the file paths coming out of walk(write_csv)
        # will become the file path outputs
        make_stimlist_nback(
          out_path_glue = out_path_glue,
          trial.order_path = trial.order_luminance,
          metadata = metadata_videos_luminance
        )
      }
    )
  )
)

## targets: fMRI stimlist materials for naturalistic jumpy videos task ----

targets_stimlists_naturalistic <- list(
  tar_target(
    name = stims_naturalistic,
    command = read_csv(metadata_videos_naturalistic) %>% 
      select(video, height, width, duration) %>% 
      left_join(read_csv(annotations_videos_naturalistic) %>% 
                  mutate(video_id = paste0(video_id, ".mp4")) %>% 
                  rename(animal_type = animal), 
                by = c("video" = "video_id")) %>% 
      filter(animal_type %in% c("dog", "frog", "spider", "cat", "food"), is.na(ignore)) %>% 
      select(-ignore)
  ),
  tar_map(
    values = tibble(order_num = sprintf("order%02d", 1:2)),
    tar_target(
      name = stimlist_naturalistic,
      command = {
        out_path_glue <- here::here("ignore",
                                    "stimlists",
                                    "naturalistic",
                                    paste0("stimlist_{sprintf('block%02d', block_num)}_",
                                           order_num, ".csv"))
        
        allowed_isis <- 6:14 # this seems to land us at 91 stimuli when multiplied into that poisson pseudo randomizer below
        
        onsets <- stims_naturalistic %>% 
          # must not be called "filename" bc that is used by psychopy to name the data files
          rename(video_id = video) %>% 
          # generate pseudo-random ISIs to jitter between videos
          # these are pre-generated to have a fixed mean and then shuffled in order
          mutate(iti_before = sample(rep(allowed_isis,
                                  times = round(n() * dpois(allowed_isis, lambda = 7) / sum(dpois(allowed_isis, lambda = 7)))))) %>% 
          # THIS randomizes the trial order
          # Since we're basically doing an old-school slow event related design anyway,
          # There is no added utility to using optimizeGA
          slice_sample(prop = 1) %>% 
          mutate(# 2024-08-26: NO LONGER RATING INSIDE SCANNER (to save time)
                 pct_task_elapsed = cumsum(duration + iti_before) / sum(duration + iti_before),
                 block_num = case_when(pct_task_elapsed < 0.34 ~ 1L,
                                       pct_task_elapsed < 0.67 ~ 2L,
                                       TRUE ~ 3L)) %>% 
          group_by(block_num) %>% 
          mutate(trial_num = 1:n(),
                 # because these ISIs set the interval _before_ the given video,
                 # manually go in and make sure the first one is always 8
                 # so that will serve as the discarded acquisitions time
                 # the last post-stim ISI is set as a constant in the PsychoPy task. don't worry about it here
                 iti_before = if_else(trial_num == 1, 8L, iti_before)) %>% 
          ungroup()
        
        onsets %>% 
          select(-pct_task_elapsed) %>% 
          nest(trials = -block_num) %>% 
          mutate(out_file = glue::glue(out_path_glue)) %$% 
          # to write the files out
          # walk2 should invisibly return out_file for targets
          walk2(out_file, trials,
                \(x, y) {write_csv(y,
                                   file = x)})
        
      },
      format = "file"
    )
  )
)

## targets: encoding model predictions (depend only on stimuli) ----

targets_encoding.models <- list(
  tar_target(
    name = activations.flynet_fixation,
    command = {
      video_path <- video_fixation
      out_path <- here::here("ignore",
                             "data",
                             "encoding",
                             "activations.flynet_fixation.csv") 
      
      with_path(conda_path,
                code = system2("python",
                               args = c(py_calc_flynet_activations,
                                        "-l 132",
                                        paste("-i", video_path),
                                        paste("-o", out_path),
                                        paste("-w", weights_flynet),
                                        "-q activations")))
      
      out_path
    },
    format = "file"
  ),
  tar_target(
    name = activations.flynet_naturalistic,
    command = {
      video_paths <- paste(videos_naturalistic_10fps, collapse = " ")
      out_path <- here::here("ignore",
                             "data",
                             "encoding",
                             "activations.flynet_naturalistic.csv") 
      
      with_path(conda_path,
                code = system2("python",
                               args = c(py_calc_flynet_activations,
                                        "-l 132",
                                        paste("-i", video_paths),
                                        paste("-o", out_path),
                                        paste("-w", weights_flynet),
                                        "-q activations")))
      
      out_path
    },
    format = "file"
  ),
  tar_target(
    name = hitprobs.flynet_naturalistic,
    command = {
      video_paths <- paste(videos_naturalistic_10fps, collapse = " ")
      out_path <- here::here("ignore",
                             "data",
                             "encoding",
                             "hitprobs.flynet_naturalistic.csv") 
      
      with_path(conda_path,
                code = system2("python",
                               args = c(py_calc_flynet_activations,
                                        "-l 132",
                                        paste("-i", video_paths),
                                        paste("-o", out_path),
                                        paste("-w", weights_flynet),
                                        "-q hit_probs")))
      
      out_path
    },
    format = "file"
  ),
  tar_target(
    name = activations.flynet_controlled,
    command = {
      # TODO: May still need to batch this in case there are too many videos to feed into one arg
      video_paths <- paste(videos_controlled_10fps, collapse = " ")
      out_path <- here::here("ignore",
                             "data",
                             "encoding",
                             "activations.flynet_controlled.csv") 
      
      with_path(conda_path,
                code = system2("python",
                               args = c(py_calc_flynet_activations,
                                        "-l 132",
                                        paste("-i", video_paths),
                                        paste("-o", out_path),
                                        paste("-w", weights_flynet),
                                        "-q activations")))
      
      out_path
    },
    format = "file"
  )
)

# targets: maps out by task x subject x run ----

path_here_fmri <- c("ignore", "data", "fmri")
path_here_derivatives <- c(path_here_fmri, "derivatives", "fmriprep-23.1.4")
tr_duration_mb8 <- 0.492
disdaq_duration <- 8 # in seconds! CONSTANT PER PHIL!
n_trs_kept_controlled <- 1001
# for pilots sub-9902-9904, this was 1407
# but before sub-0001, removed in-scanner ratings so the run is now shorter
n_trs_kept_naturalistic <- 989

map_values <- crossing(task = c("controlled", "naturalistic"),
                       # best only to include already fmriprepped subjects
                       # EXCLUDE FOR MOTION: 5, 11
                       # any other missing numbers have just not been fmriprepped yet
                       subject = c(1:4, 6:10, 12:15, 17, 19:20),
                       run = 1:5) %>% 
  filter(!(task == "naturalistic" & run > 3)) %>% 
  mutate(task = paste("task", task, sep = "-"),
         subject = sprintf("sub-%04d", subject),
         run = sprintf("run-%02d", run))

targets_fmri <- list(
  tar_map(
    unlist = FALSE,
    values = map_values,
    tar_target(
      events_raw,
      command = list.files(here::here("ignore", "data", "beh", subject, "raw"),
                           pattern = paste(subject, task, run, sep = "_"),
                           full.names = TRUE),
      format = "file"
    ),
    tar_target(
      events,
      command = {
        # n_trs fed in here is not including any TRs during the 8-second wait period
        # it counts from the first TR whose acquisition crosses the 8-second boundary
        if (task == "task-controlled") {
          parse_events_nback(events_raw, tr_duration = tr_duration_mb8, n_trs = n_trs_kept_controlled)
        } else if (task == "task-naturalistic") {
          parse_events_naturalistic(events_raw, stims_naturalistic, tr_duration = tr_duration_mb8, n_trs = n_trs_kept_naturalistic)
        }
      },
    ),
    tar_target(
      confounds,
      command = inject(here::here(!!!path_here_derivatives, subject, "func",
                                  paste(subject, task, run, "desc-confounds_timeseries.tsv", sep = "_"))),
      format = "file"
    ),
    tar_target(
      confounds_prespm,
      command = discard_confound_start(file_path = confounds,
                                       tr_duration = tr_duration_mb8,
                                       disdaq_duration = disdaq_duration),
      format = "file"
    ),
    tar_target(
      events_prespm_boxcar,
      command = {
        matlab_info <- format_events_matlab(onsets = events,
                                            raw_path = events_raw,
                                            onset_type = "boxcar")
        with_path(
          matlab_path,
          run_matlab_code(matlab_info$matlab_commands)
        )
        matlab_info$out_path
      },
      format = "file"
    ),
    tar_target(
      events_prespm_endspike,
      command = {
        matlab_info <- format_events_matlab(onsets = events,
                                            raw_path = events_raw,
                                            onset_type = "endspike")
        with_path(
          matlab_path,
          run_matlab_code(matlab_info$matlab_commands)
        )
        matlab_info$out_path
      },
      format = "file"
    ),
    tar_target(
      bold_gz,
      command = inject(here::here(!!!path_here_derivatives, subject, "func",
                                  paste(subject, task, run, "space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz", sep = "_"))),
      format = "file"
    ),
    tar_target(
      bold,
      command = {
        system2("gunzip",
                args = c("-k",
                         bold_gz))
        
        str_remove(bold_gz, ".gz")
      },
      format = "file"
    ),
    tar_target(
      bold_smoothed.4mm,
      command = spm_smooth(bold_path = bold,
                           n_trs = nrow(read_tsv(confounds)),
                           kernel = 4,
                           script = matlab_spmbatch_smooth),
      format = "file"
    )
  ),
  # stuff that's currently only being done on one task
  # move it into the main tar_map if it does eventually get done on both
  tar_eval(
    tar_target(
      target_name,
      command = {
        out_path <- here::here("ignore", "outputs", out_filename)
        
        make_encoding_timecourse(onsets = input_name,
                                 path_stim_activations = activations.flynet_naturalistic,
                                 run_duration = n_trs_kept_naturalistic * tr_duration_mb8) %>% 
          write_csv(file = out_path,
                    # because they're going into evil MATLAB
                    col_names = FALSE)
        
        out_path
        },
      format = "file"
    ),
    values = map_values %>% 
      # we only need the full single-trial events dataframes for naturalistic
      # to prepare to set the single-trial beta niftis as targets later
      filter(endsWith(task, "naturalistic")) %>% 
      mutate(across(everything(), \(x) str_replace(x, "-", ".")),
             input_name = paste("events", task, subject, run, sep = "_"),
             target_name = paste("activations.flynet", task, subject, run, sep = "_"),
             out_filename = paste(task, subject, run,"acts-flynet.csv", sep = "_"),
             across(c(input_name, target_name), syms)) %>% 
      select(input_name, target_name, out_filename)
  )
)

# targets: maps out by task x subject (combines across run) ----

# this helper function creates combined target name & input target pairs
# for tar_eval to combine across some, but not all, targets
# should be usable for combining across run but keeping task and subject
# and then combining across subject but keeping task
make_eval_values <- function (values, summarize_fmt, target_prefix) {
  values %>% 
    mutate(input_names = map2(suffix, combine_vals, \(x, y) syms(sprintf(glue::glue("%s_%s_{summarize_fmt}"), target_prefix, x, y))),
           target_name = syms(paste(target_prefix, suffix, sep = "_"))) %>% 
    select(target_name, input_names)
}

map_values_across.run <- map_values %>% 
  # so that the n column counts the number of valid runs per this subject/task
  # use the n column later to populate the valid target names to combine
  count(task, subject) %>% 
  expand(nesting(task, subject, n), 
         model_type = c("boxcar", "endspike")) %>% 
  arrange(task, subject, model_type) %>% 
  # because some of the targets that depend on this aren't SPM models
  # so for those we'll pop model_type in before the final suffix
  unite(suffix, task, subject, remove = FALSE) %>% 
  mutate(suffix = str_replace_all(suffix, "-", "."),
         combine_vals = map(n, \(x) 1:x)) %>% 
  select(-n)

targets_fmri_level1 <- list(
  tar_eval(
    expr = tar_target(
      target_name,
      command = vctrs::vec_c(!!!input_names),
      format = "file"
    ),
    values = map_values_across.run %>% 
      distinct(suffix, task, subject, combine_vals) %>% 
      make_eval_values(summarize_fmt = "run.%02d", target_prefix = "bold_smoothed.4mm")
  ),
  tar_eval(
    expr = tar_target(
      target_name,
      command = bind_rows(!!!input_names,
                          .id = "run") %>% 
        mutate(run = as.integer(run)),
    ),
    values = map_values_across.run %>% 
      # we only need the full single-trial events dataframes for naturalistic
      # to prepare to set the single-trial beta niftis as targets later
      filter(endsWith(task, "naturalistic")) %>% 
      distinct(suffix, task, subject, combine_vals) %>% 
      make_eval_values(summarize_fmt = "run.%02d", target_prefix = "events")
  ),
  tar_eval(
    expr = tar_target(
      target_name,
      command = vctrs::vec_c(!!!input_names),
      format = "file"
    ),
    values = map_values_across.run %>% 
      unite(col = "suffix", model_type, suffix) %>% 
      make_eval_values(summarize_fmt = "run.%02d", target_prefix = "events_prespm")
  ),
  tar_eval(
    expr = tar_target(
      target_name,
      command = vctrs::vec_c(!!!input_names),
      format = "file"
    ),
    values = map_values_across.run %>% 
      distinct(suffix, task, subject, combine_vals) %>% 
      make_eval_values(summarize_fmt = "run.%02d", target_prefix = "confounds_prespm")
  ),
  tar_eval(
    expr = tar_target(
      target_name,
      command = vctrs::vec_c(!!!input_names),
      format = "file"
    ),
    values = map_values_across.run %>% 
      filter(endsWith(task, "naturalistic")) %>% 
      distinct(suffix, task, subject, combine_vals) %>% 
      make_eval_values(summarize_fmt = "run.%02d", target_prefix = "activations.flynet")
  ),
  tar_eval(
    tar_target(
      name = target_name,
      command = {
        this_end_tr <- if_else(task == "task-controlled", n_trs_kept_controlled, n_trs_kept_naturalistic)
        # the output file to be tracked is the SPM.mat file
        this_model_type <- paste0("model-", model_type)
        spm_spec_est_level1(model_path = file.path("ignore", "models", task, "acq-mb8", subject, this_model_type, "smoothed-4mm"),
                            tr_duration = tr_duration_mb8,
                            trs_to_use = 1:this_end_tr + (disdaq_duration %/% tr_duration_mb8),
                            # these three need to take vectors containing the values for each run
                            bolds = a,
                            events_prespm = b,
                            confounds_prespm = c,
                            script = matlab_spmbatch_spec_est_level1)
      },
      format = "file"
    ),
    values = map_values_across.run %>% 
      select(-combine_vals) %>% 
      mutate(a = syms(paste("bold_smoothed.4mm", suffix, sep = "_")),
             b = syms(paste("events_prespm", model_type, suffix, sep = "_")),
             c = syms(paste("confounds_prespm", suffix, sep = "_")),
             target_name = syms(paste("spm_level1_smoothed.4mm", model_type, suffix, sep = "_"))) %>% 
      select(-suffix)
  ),
  tar_eval(
    tar_target(
      name = target_name,
      command = {
        out_path <- file.path(dirname(input_name), "betas_by_parcel")
        dir.create(out_path, showWarnings = FALSE)
        
        matlab_commands = c(
          assign_variable("model_path", input_name),
          assign_variable("out_folder", out_path),
          call_script(matlab_parcellate_betas)
        )
        
        with_path(
          matlab_path,
          run_matlab_code(matlab_commands)
        )
        
        list.files(out_path, full.names = TRUE)
        },
      format = "file"
    ),
    values = map_values_across.run %>% 
      filter(task == "task-naturalistic") %>% 
      select(-combine_vals) %>% 
      mutate(input_name = syms(paste("spm_level1_smoothed.4mm", model_type, suffix, sep = "_")),
             target_name = syms(paste("betas.by.parcel_smoothed.4mm", model_type, suffix, sep = "_"))) %>% 
      select(-suffix)
  )
)

## targets: maps out by task x contrast (combines across subject) ----

map_values_across.task <- map_values_across.run %>% 
  select(task, subject, model_type) %>% 
  mutate(across(everything(), \(x) str_replace_all(x, "-", "."))) %>% 
  chop(subject) %>% 
  # now there are 2 rows for the 2 tasks
  # this is just easiest to directly make a list-column of length 2 I think
  mutate(contrast = if_else(endsWith(task, "controlled"),
                            list(c("attend.animal",
                                   "dog",
                                   "frog",
                                   "spider",
                                   "above",
                                   "looming",
                                   "looming.baseline",
                                   "stimuli",
                                   "ratings")),
                            list(c("dog",
                                   "frog",
                                   "spider",
                                   "looming",
                                   "looming.baseline",
                                   "stimuli"))),
         contrast_num = map(contrast, \(x) 1:length(x))) %>% 
  unchop(cols = c(contrast, contrast_num)) %>% 
  rename(suffix = task, combine_vals = subject) %>% 
  arrange(suffix, contrast_num, model_type)

### SPM level 2 group analyses ----

targets_fmri_level2 <- list(
  # this target factory is still split out by subject, but very much pre-combining across subject
  tar_eval(
    tar_target(
      name = target_name,
      command = dep_name %>% 
        dirname() %>% 
        file.path(sprintf("con_%04d.nii", contrast_num)),
      format = "file"
    ),
    # this has 1 row per contrast x subject
    # multiple contrasts for the same subject will have the same SPM.mat dependency
    values = map_values_across.task %>% 
      # bc this one doesn't yet combine by subject. wait till the next one!
      unchop(combine_vals) %>% 
      mutate(combine_vals = str_replace_all(combine_vals, "-", "."),
             dep_name = syms(paste("spm_level1_smoothed.4mm", model_type, suffix, combine_vals, sep = "_")),
             # THIS IS THE TARGET NAME FORMAT FOR THE DOWNSTREAM ONES! U CAN DO IT
             target_name = syms(sprintf("con.%s_smoothed.4mm_%s_%s_%s", contrast, model_type, suffix, combine_vals))) %>% 
      select(target_name, dep_name, contrast_num)
  ),
  tar_eval(
    expr = tar_target(
      target_name,
      command = vctrs::vec_c(!!!input_names),
      format = "file"
    ),
    values = map_values_across.task %>% 
      mutate(input_names = pmap(list(contrast, model_type, suffix, combine_vals), \(a, b, c, d) syms(sprintf("con.%s_smoothed.4mm_%s_%s_%s", a, b, c, d))),
             target_name = syms(sprintf("con.%s_smoothed.4mm_%s_%s", contrast, model_type, suffix))) %>% 
      select(target_name, input_names)
  ),
  tar_eval(
    tar_target(
      name = target_name,
      command = {
        # the output file to be tracked is the SPM.mat file
        this_model_type <- paste0("model-", model_type)
        spm_spec_est_level2(model_path = file.path("ignore", "models", str_replace(task, "\\.", "-"), "acq-mb8", "group", this_model_type, "smoothed-4mm", contrast),
                            # this need to take a vector containing the values for each subject
                            cons = input_names,
                            con_name = contrast,
                            script = matlab_spmbatch_spec_est_level2)
      },
      format = "file"
    ),
    values = map_values_across.task %>% 
      mutate(input_names = syms(sprintf("con.%s_smoothed.4mm_%s_%s", contrast, model_type, suffix)),
             target_name = syms(sprintf("level2.%s_smoothed.4mm_%s_%s", contrast, model_type, suffix))) %>% 
      select(target_name, input_names, task = suffix, contrast, model_type)
  )
)

### canlabtools-based group analyses so help me god ----

targets_fmri_canlabtools <- list(
  tar_eval(
    expr = tar_target(
      activations.flynet_task.naturalistic_all.subs,
      # lists cannot be format = "file"
      # yeah sorta goofy but whatevs
      command = list(!!!combine_vals)
    ),
    values = map_values_across.task %>% 
      distinct(combine_vals) %>% 
      unchop(combine_vals) %>% 
      mutate(combine_vals = syms(paste0("activations.flynet_task.naturalistic_", combine_vals))) %>% 
      chop(combine_vals)
  ),
  tar_eval(
    expr = tar_target(
      bold_smoothed.4mm_task.naturalistic_all.subs,
      command = list(!!!combine_vals)
    ),
    values = map_values_across.task %>% 
      distinct(combine_vals) %>% 
      unchop(combine_vals) %>% 
      mutate(combine_vals = syms(paste0("bold_smoothed.4mm_task.naturalistic_", combine_vals))) %>% 
      chop(combine_vals)
  ),
  tar_target(
    perf_encoding.flynet_task.naturalistic_region.sc,
    command = canlabtools_fit_encoding_pls(out_path = here::here("ignore",
                                                                 "outputs", 
                                                                 "naturalistic_perf.flynet_sc.csv"),
                                           tr_duration = tr_duration_mb8,
                                           trs_to_use = 1:n_trs_kept_naturalistic + (disdaq_duration %/% tr_duration_mb8),
                                           # these three need to take vectors containing the values for each run
                                           bolds = bold_smoothed.4mm_task.naturalistic_all.subs,
                                           activations = activations.flynet_task.naturalistic_all.subs,
                                           script = matlab_fit_pls),
    format = "file"
  )
)

## targets: maps out by parcel ROI (for whole-brainy analyses) ----

map_values_across.run_by.parcel <- map_values_across.run %>% 
  expand(nesting(suffix, task, subject, model_type), 
         roi = tar_read(betas.by.parcel_smoothed.4mm_boxcar_task.naturalistic_sub.0001) %>% 
           basename() %>% 
           str_sub(start = 7L, end = -5L)) %>% 
  unite(col = "suffix", model_type, suffix, remove = FALSE)

map_values_across.task_by.parcel <- map_values_across.run_by.parcel %>% 
  select(task, subject, roi, model_type) %>% 
  mutate(across(everything(), \(x) str_replace_all(x, "-", "."))) %>% 
  chop(subject) %>% 
  rename(combine_vals = subject) %>% 
  unite(col = "suffix", model_type, task, roi, remove = FALSE)

targets_wholebrain <- list(
  tar_eval(
    tar_target(
      name = target_name,
      command = {
        all_betas <- input_name
        these_betas <- all_betas[grepl(roi, all_betas)]
        
        these_betas %>%
          read_csv() %>% 
          # VarX and Row are the names that come in from matlab writetable()
          pivot_longer(cols = starts_with("Var")) %>% 
          pivot_wider(names_from = Row) %>% 
          select(-name) %>%
          cor() %>% 
          as_tibble(rownames = "condition_row") %>% 
          pivot_longer(cols = -condition_row,
                       names_to = "condition_col",
                       values_to = "correlation")
      }
    ),
    values = map_values_across.run_by.parcel %>% 
      filter(task == "task-naturalistic") %>% 
      mutate(input_name = syms(sprintf("betas.by.parcel_smoothed.4mm_%s", suffix)),
             target_name = syms(sprintf("rdms_smoothed.4mm_%s_%s", suffix, roi))) %>% 
      select(target_name, input_name, roi)
  ),
  tar_eval(
    tar_target(
      name = target_name,
      command = bind_rows(!!!input_names,
                          .id = "target_name")
    ),
    values = map_values_across.task_by.parcel %>% 
      filter(task == "task-naturalistic") %>% 
      mutate(input_names = pmap(list(model_type, task, roi, combine_vals), \(a, b, c, d) syms(sprintf("rdms_smoothed.4mm_%s_%s_%s_%s", a, b, d, c))),
             target_name = syms(sprintf("rdms_smoothed.4mm_%s", suffix))) %>% 
      select(target_name, input_names)
  )
)

## targets: behavioral data ----

targets_beh <- list(
  # in-scanner behavioral data for the controlled 1-back task
  # this target is a single aggregated target across all subjects!!
  tar_combine(
    name = events_raw_task.controlled,
    targets_fmri[[1]][["events_raw"]],
    command = vctrs::vec_c(!!!.x)
  ),
  tar_target(
      name = beh_controlled,
      command = {
        # for whatever reason this was not working when it was externalized to another function...?
        # tar_eval problems?
        events_raw_task.controlled %>% 
          map(\(x) read_csv(x) %>% 
                mutate(frameRate = as.character(frameRate))) %>% 
          bind_rows() %>% 
          # Drop the wait-for-trigger screen and last row of ISI before the end-of-run screen
          filter(!is.na(miniblock_file)) %>% 
          # we want both of the RT measures to start from the beginning of the video display
          mutate(video_late_resp.rt = video_late_resp.rt + video_late_resp.started - video_resp.started,
                 rt = coalesce(video_resp.rt, video_late_resp.rt),
                 participant = as.integer(participant)) %>% 
          select(subj_num = participant,
                 run_num = run,
                 threatening,
                 attended,
                 miniblock_num,
                 animal_type,
                 hemifield,
                 direction, 
                 rt,
                 ends_with("rating")
          ) %>% 
          group_by(subj_num, run_num, miniblock_num) %>% 
          mutate(attended = attended[1],
                 across(ends_with("rating"), \(x) x[length(x)])) %>% 
          mutate(attended = case_match(attended, "LOCATION" ~ "hemifield", "ANIMAL" ~ "animal"),
                 resp_expected = if_else(attended == "hemifield",
                                         hemifield == lag(hemifield),
                                         animal_type == lag(animal_type)),
                 resp_expected = coalesce(resp_expected, FALSE)) %>% 
          ungroup() %>% 
          # remove the rows for the rating trials after the ratings are already filled in upward
          filter(!is.na(animal_type))
      }
  ),
  tar_target(
    name = events_raw_task.naturalistic,
    # note! this does not use the subjects listed out in map_values
    # but instead will return every file uploaded to beh that matches the pattern
    # so may include some subjects who are not yet fmriprepped
    command = list.files(here::here("ignore", "data", "beh"),
                         pattern = "task-naturalistic_beh", 
                         full.names = TRUE, 
                         recursive = TRUE),
    format = "file"
  ),
  tar_target(
    name = beh_naturalistic,
    command = events_raw_task.naturalistic %>% 
      map(read_csv) %>% 
      bind_rows() %>% 
      select(subj_num = participant, video_id, has_loom, animal_type, ends_with("rating")) %>% 
      mutate(subj_num = as.integer(subj_num)) %>% 
      filter(!is.na(video_id))
  ),
  # from 3 colleagues I was able to shake down for ratings in March 2024, lol
  tar_target(
    name = norms_naturalistic_raw,
    command = list.files(here::here("ignore", "data", "norm"),
                         recursive = TRUE,
                         full.names = TRUE) %>% 
      .[!grepl("debug", .)]
  ),
  tar_target(
    name = norms_naturalistic,
    command = norms_naturalistic_raw %>% 
      map(\(x) read_csv(x, na = c("", "NA", "None"))) %>% 
      bind_rows()
  ),
  tar_target(
    name = norms_naturalistic_qualtrics,
    command = {
      annotations_naturalistic <- read_csv(annotations_videos_naturalistic)
      annotations_controlled <- metadata_videos_nback %>% 
        left_join(read_csv(qualtrics.ids_videos_controlled), by = "filename") %>% 
        # to get the cols in the same format as the naturalistic annotations
        mutate(has_loom = as.numeric(direction == "looming")) %>% 
        select(video_id = filename, qualtrics_id, animal = animal_type, has_loom, hemifield)
      
      annotations <- bind_rows(naturalistic = annotations_naturalistic, 
                               controlled = annotations_controlled,
                               .id = "stim_type")
        
      
      # TODO: Make sure your qualtrics API token is reproducible yet safe from prying eyes
      # TODO also: set up another target so that targets will think the survey is updated when it is
      qualtrics_data <- fetch_survey("SV_6nZOkqUU0MEAnMa",
                                     include_display_order = FALSE)
      
      video_qualtrics_ids <- qualtrics_data %>% 
        select(contains("caption")) %>% 
        sjlabelled::get_label() 
      
      df_video_qualtrics_ids <- tibble(video_num = names(video_qualtrics_ids),
                                       qualtrics_id = video_qualtrics_ids %>% 
                                         str_split_i(pattern = " ", i = 1)) %>% 
        separate_wider_regex(cols = video_num,
                             patterns = c(video_num = "[[:digit:]]+",
                                          "_[[:alpha:]]+",
                                          block_num = "[[:digit:]]+")) %>% 
        mutate(across(ends_with("num"), as.integer))
      
      qualtrics_data %>% 
        filter(Status != "Survey Preview", participantId_check != "test") %>% 
        select(-contains("quit"), -starts_with("Recipient"), -DistributionChannel, -UserLanguage) %>% 
        pivot_longer(cols = c(contains("caption"), 
                              contains("pleasantness"), 
                              contains("arousal"), 
                              contains("fear")), 
                     names_to = c("video_num", ".value", "block_num"), 
                     names_pattern = "([[:digit:]]+)_([[:alpha:]]+)([[:digit:]]+)",
                     names_transform = list(video_num = as.integer,
                                            block_num = as.integer)) %>% 
        left_join(df_video_qualtrics_ids, by = c("video_num", "block_num")) %>% 
        left_join(annotations,
                  by = "qualtrics_id")
    }
  )
)

list(
  targets_scripts,
  targets_stimuli,
  targets_stimlists_nback,
  targets_stimlists_naturalistic,
  targets_encoding.models,
  targets_fmri,
  targets_fmri_level1,
  targets_fmri_level2,
  targets_wholebrain,
  targets_fmri_canlabtools,
  targets_beh
)
