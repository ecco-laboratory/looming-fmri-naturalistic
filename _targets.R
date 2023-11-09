# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(tibble)

# Set target options:
tar_option_set(
  packages = c("withr",
               "matlabr",
               "tidyverse",
               "magrittr",
               "glue") # packages that your targets need to run
  # format = "qs", # Optionally set the default storage format. qs is fast.
  #
  # For distributed computing in tar_make(), supply a {crew} controller
  # as discussed at https://books.ropensci.org/targets/crew.html.
  # Choose a controller that suits your needs. For example, the following
  # sets a controller with 2 workers which will run as local R processes:
  #
  #   controller = crew::crew_controller_local(workers = 2)
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
tar_source(c("R/utils/",
             "R/make-stimlist.R"))
# source("other_functions.R") # Source other scripts as needed.

## other useful global vars ----
# My dear little conda env
conda_path <- "/home/mthieu/Repos/emonet-py/env/bin"
matlab_path <- "/opt/MATLAB/2022a/bin"
# fMRI parameters that pass into multiple targets
tr_length <- 2 # seconds
nback_rating_duration <- 12 # tr units
nback_blocks_per_run <- 4
nback_block_threatenings <- c("dog", "frog", "spider", "above", "below")

# target definition

## targets: scripts in other project folders ----

targets_ext <- list(
 tar_target(
   name = py_make_looming_video,
   command = "/home/mthieu/Repos/emonet-py/python/myutils/make_looming_video.py",
   format = "file"
 ),
 tar_target(
   name = py_get_video_metadata,
   command = here::here("python", "get_video_metadata.py"),
   format = "file"
 ),
 tar_target(
   name = matlab_optimize_ga_trial_order,
   command = here::here("matlab", "optimize_ga_trial_order.m"),
   format = "file"
 ),
 tar_target(
   name = matlab_optimize_ga_trial_order_naturalistic,
   command = here::here("matlab", "optimize_ga_trial_order_naturalistic.m"),
   format = "file"
 )
)

## targets: looming stimuli of various kinds ----

targets_stimuli <- list(
  tar_target(
    name = images_controlled,
    command = list.files(here::here("ignore", "stimuli", "images", "use"),
                         full.names = TRUE),
    format = "file"
  ),
  tar_target(
    name = videos_controlled,
    command = with_path("/home/mthieu/Repos/emonet-py/env/bin", {
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
                 full.names = TRUE)
    }),
    format = "file"
  ),
  tar_target(
    name = folder_videos_naturalistic,
    command = "/home/data/eccolab/looming/stimuli/youtube/",
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
      with_path(conda_path, 
                system2("python",
                        args = c(py_get_video_metadata,
                                 paste("--path", folder_videos_naturalistic)))
      )
      paste0(folder_videos_naturalistic, "metadata.csv")
    },
    format = "file"
  ),
  tar_target(
    name = annotations_videos_naturalistic,
    command = "/home/data/eccolab/looming/stimuli/youtube/annotations.csv",
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
             filename = sprintf("miniblock%02d.csv", miniblock_num)) %>% 
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
      filter(animal_type %in% c("dog", "frog", "spider"), is.na(ignore)) %>% 
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
        
        allowed_isis <- 6:16
        
        naturalistic_rating_duration <- 4*3 + 8 # 3 ratings plus 8-second washout
        onsets <- stims_naturalistic %>% 
          rename(filename = video) %>% 
          # generate pseudo-random ISIs to jitter between video offset and rating onset
          # these are pre-generated to have a fixed mean and then shuffled in order
          mutate(isi = sample(rep(allowed_isis,
                                  times = round(n() * dpois(allowed_isis, lambda = 7) / sum(dpois(allowed_isis, lambda = 7)))))) %>% 
          # THIS randomizes the trial order
          # Since we're basically doing an old-school slow event related design anyway,
          # There is no added utility to using optimizeGA
          slice_sample(prop = 1) %>% 
          mutate(# need to add in the constant rating interval per trial
                 # it nonlinearly changes the pct elapsed bc videos vary in length
                 pct_task_elapsed = cumsum(duration + isi + naturalistic_rating_duration) / max(cumsum(duration + isi + naturalistic_rating_duration)),
                 block_num = case_when(pct_task_elapsed < 0.33 ~ 1L,
                                       pct_task_elapsed < 0.66 ~ 2L,
                                       TRUE ~ 3L)) %>% 
          group_by(block_num) %>% 
          mutate(trial_num = 1:n()) %>% 
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

list(
  targets_ext,
  targets_stimuli,
  targets_stimlists_nback,
  targets_stimlists_naturalistic
)
