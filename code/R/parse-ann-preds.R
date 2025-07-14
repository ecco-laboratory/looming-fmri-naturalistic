# targets-compatible helper functions for pulling ANN-predicted categories for stimulus videos

# n_top determines top-n accuracy. it will always return top 1, but you can also return a higher top n alongside
get_alexnet_guesses <- function (path_alexnet_activations, stim_labels, path_imagenet_categories, n_top = 1) {
  imagenet_categories <- read_csv(path_imagenet_categories)
  
  guesses <- read_csv(path_alexnet_activations) %>% 
    right_join(stim_labels, 
               by = "video") %>% 
    # relabel a couple of them to match the category labels from the imagenet key
    mutate(animal_type_imagenet = map(animal_type, \(x) {
      if (x == "spider") {
        return ("arachnid")
        } else if (x == "food") {
          return (c("food", "fruit", "vegetable"))
        } else {
          return (x)
        }
      })) %>% 
    nest(.by = c(video, frame, animal_type, animal_type_imagenet, looming), .key = "activations") %>% 
    # unlisting a 1-row tibble of all numbers turns it into a (named) vector
    mutate(activations = map(activations, unlist),
           # so that the highest activations will be first to reduce my own working memory load
           indices = map(activations, \(x) sort(x, decreasing = TRUE, index.return = TRUE) %>% 
                           # keep just the sorted indices! these are R indices so they go from 1!!!
                           pluck("ix") %>% 
                           # and keep only the top n indices (to save space from here)
                           magrittr::extract(1:n_top))) %>% 
    select(-activations) %>% 
    mutate(guesses = map(indices, \(x) imagenet_categories$categories[x]),
           correct_top.1 = map2_lgl(animal_type_imagenet, guesses, \(x, y) any(x == y[1])),
           # remember, %in% vectorizes over the LEFT argument
           correct_top.n = map2_lgl(animal_type_imagenet, guesses, \(x, y) any(x %in% y))) %>% 
    select(-where(is.list))
  
  # if the correct columns are redundant, keep only 1
  if (n_top == 1) {
    guesses %<>%
      select(-correct_top.n)
  }
  return (guesses)
}

get_flynet_guesses <- function (path_flynet_activations, stim_labels) {
  # 2025-07-09: I went and checked the Numpy data direct from Baohua Zhou,
  # and for the 256-unit model's classifier layer to calculate P(hit),
  # these are the weight and bias values to multiply by after summing activations and before sigmoiding
  # I don't feel like importing them direct from Python... pls trust.
  zhou2022_classifier_weight <- 1
  zhou2022_classifier_bias <- -2.9806
  
  # For the way we extracted FlyNet information from Study 3 of the looming secondary analysis paper
  # For each video x RF unit, get the slope of activation over time
  # Then if looming is linearly decodable from RF slopes, it counts as looming?
  guesses <- read_csv(path_flynet_activations) %>% 
    pivot_longer(cols = -c(video, frame),
                 names_to = "rf",
                 values_to = "activation") %>%
    inner_join(stim_labels, 
               by = "video") %>% 
    nest(activations = c(frame, activation)) %>% 
    mutate(coefs = map(activations,
                       ~lm(activation ~ scale(frame, scale = FALSE), data = .) %>% 
                         pluck("coefficients"),
                       .progress = list(name = "RF activation slopes"))) %>% 
    select(-activations) %>% 
    unnest_wider(coefs) %>% 
    rename(intercept = "(Intercept)", slope = "scale(frame, scale = FALSE)") %>% 
    pivot_wider(names_from = rf,
                values_from = c(intercept, slope))
    
    
    nest(.by = c(video, frame, animal_type, looming), .key = "activations") %>% 
    # first, sum across units for each frame
    mutate(activation_sum = map_dbl(activations, \(x) sum(unlist(x)))) %>% 
    select(-activations) %>% 
    # then adjust by the linear transformation weight and bias
    # and convert to probability for each frame
    mutate(p_hit = plogis(zhou2022_classifier_weight*activation_sum + zhou2022_classifier_bias)) %>% 
    nest(.by = c(video, animal_type, looming), .key = "frames") %>% 
    mutate(p_hit_video = map_dbl(frames, \(x) mean(x$p_hit)))
  
  return (guesses)
}
