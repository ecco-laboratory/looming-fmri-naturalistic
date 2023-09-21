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
  packages = c("tibble",
               "dplyr",
               "tidyr",
               "purrr") # packages that your targets need to run
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
tar_source("R/utils/")
# source("other_functions.R") # Source other scripts as needed.

# Replace the target list below with your own:

# scripts in other project folders
targets_ext <- list(
 tar_target(
   name = "py_make_looming_video",
   command = "/home/mthieu/Repos/emonet-py/python/myutils/make_looming_video.py",
   format = "file"
 ) 
)

# looming stimuli of various kinds
targets_stimuli <- list(
  tar_target(
    name = images_controlled,
    command = list.files(here::here("ignore", "stimuli", "images", "use"),
                         full.names = TRUE),
    format = "file"
  ),
  tar_target(
    name = videos_controlled,
    command = withr::with_path("/home/mthieu/Repos/emonet-py/env/bin", {
      out_path <- here::here("ignore", "stimuli", "videos", "controlled")
      for (img in images_controlled) {
        cat("starting image:", crayon::green(img), fill = TRUE)
        for (direction in c("looming", "receding")) {
          cat("direction:", crayon::magenta(direction), fill = TRUE)
          for (angle in c("above", "below")) {
            # print any of the associated stimuli that need to be tracked
            # run the python script
            cat("approach angle:", crayon::magenta(angle), fill = TRUE)
            system2("python", args = c(py_make_looming_video,
                                       paste("--infile", img),
                                       paste("--outpath", out_path),
                                       # Adult Ps sat 40 cm away from monitor
                                       # and I think it needs to come 4 obj-widths away
                                       "--diststart 15",
                                       "--distend 0.1",
                                       "--posstart middle",
                                       paste("--direction", direction),
                                       paste("--angle", angle),
                                       "--loomtime 0.5"
            )) 
          }
        }
      }
      # print the path to the output videos
      list.files(out_path,
                 full.names = TRUE)
    }),
    format = "file"
  )
)

targets_stimlists <- list(
  tar_target(
    name = metadata_videos_controlled,
    command = tibble(filename = videos_controlled) %>% 
      separate_wider_delim(filename, delim = "/", names = c(rep(NA, 9), "filename")) %>% 
      # don't drop the filename (without folders)
      separate_wider_delim(filename, delim = "_", names = c("image_num", "direction", NA, "hemifield", NA, NA), cols_remove = FALSE) %>% 
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
  tar_target(
    name = stimlist_nback_base,
    command = expand_grid(attended = c("animal_type", "hemifield"),
                          threatening = c("dog", "spider", "frog", "above", "below"),
                          animal_type = c("dog", "spider", "frog"),
                          hemifield = c("above", "below")) %>% 
      mutate(predictable = if_else(threatening %in% c("above", "below"), "hemifield", "animal_type"),
             # this is a constant--assigns animal_nums 1:6 for every row
             animal_num = map(animal_type, \(x) 1:6)) %>% 
      unchop(animal_num) %>% 
      nest(trials = -c(attended, predictable, threatening))
  ),
  tar_map(
    values = tibble(stimlist_num = 1:10),
    tar_target(
      name = stimlist_nback,
      command = stimlist_nback_base %>% 
        # shuffle between block
        slice(get_balanced_latin_square_order(nrow(stimlist_nback_base), stimlist_num)) %>% 
        mutate(block_num = 1:n()) %>% 
        unnest(trials) %>% 
        # nest into the little groups of 6 to assign animal_nums to loom vs recede
        nest(trials = -c(attended, predictable, threatening, animal_type, hemifield)) %>% 
        # this is a temp helper column for generating the stim-specific loom/recedes
        mutate(threat_status = case_when(
          (predictable == "animal_type" & animal_type == threatening) | (predictable == "hemifield" & hemifield == threatening) ~ "threatening",
          (predictable == "animal_type" & animal_type != threatening) | (predictable == "hemifield" & hemifield != threatening) ~ "safe",
          TRUE ~ NA_character_),
          trials = map2(trials, threat_status, \(x, y) {
            if (y == "threatening") {
              directions <- sample(rep(c("looming", "receding"), times = c(4, 2)))
            } else {
              directions <- sample(rep(c("looming", "receding"), times = c(2, 4)))
            }
            x$direction <- directions
            return (x)
          })) %>%
        select(-threat_status) %>% 
        # bind the stimulus file names on
        unnest(trials) %>% 
        left_join(metadata_videos_controlled,
                  by = c("animal_type", "direction", "hemifield", "animal_num")) %>% 
        # now shuffle within block
        nest(trials = -c(attended, predictable, threatening)) %>% 
        mutate(trials = map(trials, \(x) x %>% 
                              slice_sample(prop = 1) %>% 
                              mutate(trial_num = 1:nrow(.)))) %>% 
        unnest(trials)
    )
  )
)

list(
  targets_ext,
  targets_stimuli,
  targets_stimlists
)
