#' Simulates judgments of frequency and duration with a competitive learning
#' network based on the PASS-family
#'
#' @param patterns matrix with input patterns, one row is one pattern
#' @param frequency presentation frequency for each pattern in the matrix
#' @param duration presentation duration for each pattern in the matrix
#' @param lrate_onset learning rate at the onset of a stimulus
#' @param lrate_drop_time point at which the learning rate drops, must be lower
#'   than duration
#' @param lrate_drop_perc how much the learning rate drops at lrate_drop_time
#' @param n_runs: number of simulations to be run, default is 100
#' @param n_output_units number of output units, defaults to number of input
#'   units
#' @param pulses_per_second how many time steps should be simulated per second
#' @param random_order whether stimuli should be presented in random order,
#'   default is TRUE
#' @return list with following elements:
#'   output_activation_sum: the sum of the activation strengths of the output
#'   units for each input pattern
#'   weight_matrix: final weight_matrix
#'   pres_matrix: presentation matrix
#' @export
run_sim <- function(patterns, frequency, duration, lrate_onset, lrate_drop_time,
                    lrate_drop_perc, n_runs = 100,
                    n_output_units = ncol(patterns), pulses_per_second = 1,
                    random_order = TRUE){
  output_activation_sum <- NULL
  for (j in 1:n_runs) {
    n_input_units <- ncol(patterns)
    weight_matrix <- init_weight_mtrx(n_input_units, n_output_units)
    pres <- create_pres_matrix(patterns, frequency, duration, lrate_onset,
                               lrate_drop_time, lrate_drop_perc,
                               pulses_per_second = pulses_per_second,
                               random_order = random_order)
    # present input i and update weights
    for (i in 1:nrow(pres$input)) {
      weight_matrix <- updt_winner_weights(pres$input[i, ], weight_matrix,
                                           pres$lrate[i])
    }
    # output for simulation j
    output_activation_sum <- rbind(output_activation_sum,
                                   calc_output_sum(patterns, weight_matrix))
  }
  final_list <- list("output" = output_activation_sum,
                    "weight_matrix" = weight_matrix,
                    "pres_matrix" = pres)
  return(final_list)
}

#' Create presentation matrix
#' @inheritParams run_sim
#' @return list with two lists: the presentation matrix and the learning weights
#'   (both in one random order if random_order = TRUE)
#' @export
create_pres_matrix <- function(patterns, frequency, duration, lrate_onset,
                               lrate_drop_time, lrate_drop_perc,
                               pulses_per_second, random_order = TRUE){
  attention <-  get_attention(duration, lrate_onset, lrate_drop_time,
                              lrate_drop_perc, pulses_per_second)
  pres_list <- list()
  final_result <- list()
  # produce list (stimuli) of lists (presentation time) of smallest entity
  for (i in 1:nrow(patterns)) {
    duration_pulse <- duration[i] * pulses_per_second
    pres_list[[i]] <- list(cbind(matrix(rep(patterns[i, ], duration_pulse),
                                       nrow = duration_pulse, byrow = TRUE),
                                attention[[i]]))
    pres_list[[i]] <- lapply(pres_list[i], rep, frequency[i])
  }
  # unlist so that all presentations remain lists, apply unlist two times
  pres_list <- unlist(pres_list, recursive = F)
  pres_list <- unlist(pres_list, recursive = F)

  # randomize presentations
  if (random_order) {
    random_index <- sample(1:length(pres_list), length(pres_list))
  } else {
    random_index <- 1:length(pres_list)
  }
  pres_list <- pres_list[random_index]

  # create df
  pres_df <- do.call(rbind, pres_list)

  # split into two lists
  final_result[["input"]] <- pres_df[, 1:(ncol(pres_df) - 1)]
  final_result[["lrate"]] <- pres_df[, ncol(pres_df)]
  return(final_result)
}

#' Initialize weight matrix for competitive learning network. The function draws
#' from a normal distribution and then normalizes the weights, such that the
#' sum of weights is 1.
#'
#' @n_input_units: number of input units
#' @n_output_units: number of output units
#' @mean mean of normal distribution
#' @sd sd of normal distribution
#' @return mtrx with n_input_units rows and n_output_units columns, the sum of
#'   every column is 1
#' @export
init_weight_mtrx <- function(n_input_units, n_output_units, mean = 0.5,
                             sd = 0.005){
  n_weights <- n_output_units * n_input_units
  # initialize random weight mtrx
  weight_mtrx <- matrix(rnorm(n_weights, mean, sd), n_output_units,
                        byrow = TRUE)
  # weights for every output unit (row sum) must equal 1
  weight_mtrx <- prop.table(weight_mtrx, margin = 1)
  return(weight_mtrx)
}

#' Updates weights for competitive learning network algorithm
#'
#' @input: input vector
#' @weight_matrix: weight matrix of network
#' @lrate: learning rate
#' @return new weight matrix for step t + 1
#' @export
updt_winner_weights <- function(input, weight_matrix, lrate){
  output <- weight_matrix %*% input
  winner <- which(output == max(output))
  if (length(winner) > 1) {
    print ("more than one winner, random selection")
    winner <- sample(winner, 1)
  }
  delta_w <- lrate * (input / sum(input) - weight_matrix[winner, ])
  weight_matrix[winner, ] <- weight_matrix[winner, ] + delta_w
  return(weight_matrix)
}

#' Calculates sum of activation in output units for input patterns
#'
#' @inheritParams updt_winner_weights
#'
#' @return sum of activations of output units for input patterns
#' @export
calc_output_sum <- function(inputs, weight_matrix){
  if (class(inputs) == "matrix") {
    return(colSums(apply(inputs, 1, function(x) weight_matrix %*% x)))
  }
  if (is.vector(inputs) == TRUE) {
    return(sum(weight_matrix %*% inputs))
  }
}

#' Calculates attention (learning rate) development over time
#'
#' @inheritParams run_sim
#'
#' @return list of attention values (learning rate) for every time pulse for
#'   every stimulus
get_attention <- function(duration, lrate_onset, lrate_drop_time,
                          lrate_drop_perc, pulses_per_second = 1){
  duration_pulses <- duration * pulses_per_second
  lrate_drop_pulse <- lrate_drop_time * pulses_per_second
  low_value <- lrate_onset * lrate_drop_perc
  y <- list()
  for (i in 1:length(duration_pulses)) {
    high <- rep(lrate_onset, lrate_drop_pulse - 1)
    n_low_pulses <- duration_pulses[i] - lrate_drop_pulse + 1
    if (n_low_pulses > 0) {
      low <- rep(low_value, n_low_pulses)
    } else {
      low <- NULL
    }
    y[[i]] <- c(high, low)
  }
  return(y)
}

#' Runs one "experiment" with several "participants", which means that a
#' simulation is run for each "participant" and correlative effect sizes are
#' calculated; stimuli are orthogonal
#'
#' @inheritParams run_sim
#' @number_of_participants corresponds with number of simulations run
#' @cor_noise_sd the amount of noise added to the final activations of the
#'   network, set to 0 if you do not want any noise
#' @export
#' @importFrom magrittr "%>%
do_exp <- function(duration, frequency, lrate_onset, lrate_drop_time,
                   lrate_drop_perc, number_of_participants, cor_noise_sd){
  inputs <- diag(length(duration))

  sim_low_a <- run_sim(patterns = inputs, lrate_onset = lrate_onset,
                       lrate_drop_time = 2, frequency = frequency,
                       duration = duration, n_runs = number_of_participants,
                       lrate_drop_perc = lrate_drop_perc,
                       n_output_units = length(duration), random_order = TRUE)

  strength <- as.data.frame(sim_low_a$output)
  # create useful data structure
  d <- cbind(id = 1:nrow(strength), strength)
  # with factor_key = TRUE we can leave the ordering of the columns
  d <- tidyr::gather(d, key = "condition", value = "dv_activation", -id,
              factor_key = TRUE)
  d <- dplyr::arrange(d, id, condition)
  d <- cbind(d, iv_f = frequency, iv_d = duration,
             iv_total_d = duration * frequency)

  # add noise to dv
  d <- dplyr::mutate(d,
              dv_activation = dv_activation +
                rnorm(length(dv_activation), mean = 0, sd = cor_noise_sd))
  # calculate effect sizes
  r_contrast <- d %>%
    dplyr::summarize(f_dv = cor(iv_f, dv_activation),
              td_dv = cor(iv_total_d, dv_activation),
              d_dv = cor(iv_d, dv_activation))
  r_contrast
}