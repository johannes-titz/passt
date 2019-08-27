library(dplyr)
library(tidyr)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(gridExtra)
source("pass-t.R")

run_exp <- function(presTime, presRate, lrate, number_of_participants,
                   attentionC, noise_sd){
  # runs a single experiment with some predefined parameters that were used
  # in Titz & Sedlmeier (submitted)
  #
  # Args:
  #   presRate: presentation rates for stimuli
  #   presTime: presentation times for stimuli
  #   lrate: determines how fast the neural network will learn, value
  #     must be between 0-1, values too low and too high will not work
  #   number_of_participants: number of different simulation runs, comparable
  #     to having several participants in an experiment.
  #   attentionC: parameter of attention function; time point at which attention
  #     drops
  #   noise_sd: SD of Gaussian noise added to final activations of stimuli
  #     (standing for the judgments)
  #
  # Returns:
  #   a one-row data frame that contains correlations (constrasts) between
  #   presentation rates, times and the final activations (judgments)
  attentionFreq <- c(rep("constant", length(presTime)))
  inputs <- diag(length(presTime))

  sim_low_a <- SimulatePassT(inputMatrix = inputs,
                             lrate = lrate,
                             presRate = presRate, presTime = presTime,
                             runs = number_of_participants,
                             attentionFunction = attentionFreq,
                             attentionC = attentionC,
                             n_output_units = length(presTime),
                             random_order = T)

  sim_low_aStrength <- as.data.frame(sim_low_a$output)

  # create useful data structure
  d <- cbind(id = 1:nrow(sim_low_aStrength), sim_low_aStrength)
  # with factor_key = T we can leave the ordering of the columns
  d <- gather(d, key = "condition", value = "dv_activation", -id,
              factor_key = T)
  d <- arrange(d, id, condition)
  d <- cbind(d, iv_f = presRate, iv_d = presTime,
             iv_total_d = presTime * presRate)

  d <- mutate(d,
              dv_activation = dv_activation +
                rnorm(length(dv_activation),
                      mean = 0,
                      sd = noise_sd))
  # calc effects (correlations but also partial correlations to study
  # methodological effects)
  cont <- d %>%
    summarize(ff = cor(iv_f, dv_activation),
              tdf = cor(iv_total_d, dv_activation),
              df = cor(iv_d, dv_activation),
              df_semi = ppcor::pcor.test(iv_d, dv_activation, iv_f)$estimate,
              ff_semi = ppcor::pcor.test(iv_f, dv_activation, iv_d)$estimate,
              tdf_semi = ppcor::pcor.test(iv_total_d, dv_activation,
                                          iv_f)$estimate)
  cont
}

run_exps <- function(number_of_exps, presTime, presRate, lrate,
                    number_of_participants, attentionC, noise_sd){
  # runs several experiments by using run_exp function with some predefined
  # parameters that were used in Titz & Sedlmeier (submitted)
  #
  # Args:
  #   number_of_exps: number of experiments to run
  #   presRate: presentation rates for stimuli
  #   presTime: presentation times for stimuli
  #   lrate: determines how fast the neural network will learn, value
  #     must be between 0-1, values too low and too high will not work
  #   number_of_participants: number of different simulation runs, comparable
  #     to having several participants in an experiment.
  #   attentionC: parameter of attention function; time point at which attention
  #     drops
  #   noise_sd: SD of Gaussian noise added to final activations of stimuli
  #     (standing for the judgments)
  #
  # Returns:
  #   a data frame that contains correlations (constrasts) between presentation
  #   rates, times and the final activations (judgments), number of rows equals
  #   number of experiments
  conts <- NULL
  for (i in seq(number_of_exps)){
    conts[[i]] <- run_exp(presTime, presRate, lrate,
                         number_of_participants, attentionC, noise_sd)
    cat("#")
  }
  plyr::ldply(conts, "data.frame")
}

# prestnation time and rate
presTime <- c(4, 2, 1, 8, 4, 2, 12, 6, 3)
presRate <- c(2, 4, 8, 2, 4, 8, 2, 4, 8)

# JOF --------------------------------------------------------------------------
# produce parameter space
pars_jof <- expand.grid(attention_drop = seq(0.02, 0.16, 0.02),
                    attention_init = 0.05,
                    cor_noise_sd = seq(0, 0.2, 0.01))

effects_jof <- Map(function(x, y, z) run_exps(1, presTime, presRate, y, 10, x,
                                             z),
            x = pars_jof$attention_drop, y = pars_jof$attention_init,
            z = pars_jof$cor_noise_sd)
effects_jof <- plyr::ldply(effects_jof, "data.frame")

res_jof <- bind_cols(effects_jof, pars_jof)
res_jof <- mutate(res_jof, "dv" = "jof")

# JOD --------------------------------------------------------------------------
# produce parameter space
pars_jod <- expand.grid(attention_drop = seq(0.16, 0.6, 0.02),
                     attention_init = 0.05,
                     cor_noise_sd = seq(0, 0.2, 0.01))

effects_jod <- Map(function(x, y, z) run_exps(1, presTime, presRate, y, 10, x,
                                             z),
             x = pars_jod$attention_drop, y = pars_jod$attention_init,
             z = pars_jod$cor_noise_sd)
effects_jod <- plyr::ldply(effects_jod, "data.frame")

res_jod <- bind_cols(effects_jod, pars_jod)
res_jod <- mutate(res_jod, "dv" = "jod")

# change legend to percent
res_jof <- mutate(res_jof, attention_drop = attention_drop * 100)
res_jod <- mutate(res_jod, attention_drop = attention_drop * 100)

res <- bind_rows(res_jof, res_jod)

# plotting----------------------------------------------------------------------
p_jof <- ggplot(res_jof, aes(ff, tdf, color = attention_drop,
                         size = cor_noise_sd)) +
  geom_point(alpha = 1) +
  guides(color = guide_colorbar(order=1)) +
  scale_size_continuous("Gaussian noise [SD]", range=c(5,1)) +
  scale_color_gradient2(name = "Attention [%]", low = "white", high = "black") +
  theme_bw() +
  theme(legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill=alpha('white', 0.3))) +
  scale_x_continuous(name = "Influence of frequency",
                     labels = c("0", ".25", ".50", ".75", "1"),
                     limits = c(0, 1)) +
  scale_y_continuous("",
                     labels = c("0", ".25", ".50", ".75", "1"),
                     limits = c(0, 1)) +
  ggtitle("Judgment of Frequency") + theme(plot.title = element_text(hjust = 0.5))

p_jod <- ggplot(res_jod, aes(ff, tdf, color = attention_drop, size = cor_noise_sd)) +
  geom_point(alpha = 1) +
  scale_size_continuous(guide = F, "Noise [SD]", range=c(5,1)) +
  scale_color_gradient2(name = "Attention [%]", low = "white", high = "black") +
  theme_bw() +
  theme(legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill=alpha('white', 0.3))) +
  scale_x_continuous(name = "Influence of frequency",
                     labels = c("0", ".25", ".50", ".75", "1"),limits = c(0, 1)) +
  scale_y_continuous("Influence of duration",
                     labels = c("0", ".25", ".50", ".75", "1"),limits = c(0, 1))+
  ggtitle("Judgment of Duration") + theme(plot.title = element_text(hjust = 0.5))

# save combined plot
pdf("sim_bw.pdf", width = 10, height = 5)
grid.arrange(p_jod, p_jof, ncol = 2)
dev.off()
