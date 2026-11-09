rm(list = ls())
gc()
graphics.off()

library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)
library(forcats)

## 1. Signal Generation Functions -------------------------------------------------

smooth_H <- function(t, t0, gamma = 0.2) {
  0.5 * (tanh(gamma * (t - t0)) + 1)
}

generate_signal_trace <- function(type, times) {
  df <- tibble(time = times, S1 = 0, S2 = 0)
  
  if (type == "alternating") {
    pulse_width <- 60
    pulse_starts <- seq(0, max(times), by = pulse_width)
    for (i in seq_along(pulse_starts)) {
      start <- pulse_starts[i]
      end <- start + pulse_width
      if (i %% 2 == 1) {
        df$S1 <- df$S1 + ifelse(times >= start & times < end, 1, 0)
      } else {
        df$S2 <- df$S2 + ifelse(times >= start & times < end, 1, 0)
      }
    }
  } else if (type == "simultaneous") {
    pulse_centers <- seq(50, max(times) - 50, by = 150)
    pulse_width <- 60
    for (center in pulse_centers) {
      active <- ifelse(times >= (center - pulse_width / 2) & times < (center + pulse_width / 2), 1, 0)
      df$S1 <- df$S1 + active
      df$S2 <- df$S2 + active
    }
  }
  
  df %>% mutate(across(c(S1, S2), ~ pmin(pmax(.x, 0), 1)))
}

## 2. Models ----------------------------------------------------------------------

competitive_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    sig <- parameters$signal_trace[which.min(abs(parameters$signal_trace$time - t)), ]
    S1 <- sig$S1
    S2 <- sig$S2
    
    free_T <- T_total - P1T1 - P2T2
    free_P1 <- P1_total - P1T1
    free_P2 <- P2_total - P2T2
    
    dP1T1 <- k_on1 * free_P1 * free_T * S1 - k_off1 * P1T1
    dP2T2 <- k_on2 * free_P2 * free_T * S2 - k_off2 * P2T2
    dOutput <- k_phos * P1T1 - k_dephos * P2T2
    
    # Conservation verification
    T_conserved <- T_total - (free_T + P1T1 + P2T2)
    if (abs(T_conserved) > 1e-6) warning("T not conserved in competitive model at t=", t)
    
    list(c(dP1T1, dP2T2, dOutput))
  })
}

non_competitive_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    sig <- parameters$signal_trace[which.min(abs(parameters$signal_trace$time - t)), ]
    S1 <- sig$S1
    S2 <- sig$S2
    
    # Two independent T pools
    free_T1 <- T1_total - P1T1 - P1T1P2T2
    free_T2 <- T2_total - P2T2 - P1T1P2T2
    
    free_P1 <- P1_total - P1T1 - P1T1P2T2
    free_P2 <- P2_total - P2T2 - P1T1P2T2
    
    # Primary binding to independent sites
    dP1T1 <- k_on1 * free_P1 * free_T1 * S1 - k_off1 * P1T1 - 
      k_on2 * free_P2 * P1T1 * S2 + k_off12 * P1T1P2T2
    
    dP2T2 <- k_on2 * free_P2 * free_T2 * S2 - k_off2 * P2T2 - 
      k_on1 * free_P1 * P2T2 * S1 + k_off12 * P1T1P2T2
    
    # Ternary complex formation
    dP1T1P2T2 <- k_on2 * free_P2 * P1T1 * S2 + k_on1 * free_P1 * P2T2 * S1 - 
      2 * k_off12 * P1T1P2T2
    
    dOutput <- k_phos * (P1T1 + P1T1P2T2) - k_dephos * (P2T2 + P1T1P2T2)
    
    # Conservation checks
    T1_conserved <- T1_total - (free_T1 + P1T1 + P1T1P2T2)
    T2_conserved <- T2_total - (free_T2 + P2T2 + P1T1P2T2)
    if (abs(T1_conserved) > 1e-6 || abs(T2_conserved) > 1e-6) {
      warning("T not conserved in non-competitive model at t=", t)
    }
    
    list(c(dP1T1, dP2T2, dP1T1P2T2, dOutput))
  })
}

somewhat_noncompetitive_improved <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    sig <- parameters$signal_trace[which.min(abs(parameters$signal_trace$time - t)), ]
    S1 <- sig$S1
    S2 <- sig$S2
    
    # Single shared T pool
    free_T <- T_total - P1T1 - P2T2 - P1T1P2T2
    free_P1 <- P1_total - P1T1 - P1T1P2T2
    free_P2 <- P2_total - P2T2 - P1T1P2T2
    
    # Primary competitive binding
    dP1T1 <- k_on1 * free_P1 * free_T * S1 - k_off1 * P1T1 - 
      k_on2 * free_P2 * P1T1 * S2 + k_off12 * P1T1P2T2
    
    dP2T2 <- k_on2 * free_P2 * free_T * S2 - k_off2 * P2T2 - 
      k_on1 * free_P1 * P2T2 * S1 + k_off12 * P1T1P2T2
    
    # Secondary binding after primary binding occurs
    dP1T1P2T2 <- k_on2 * free_P2 * P1T1 * S2 + k_on1 * free_P1 * P2T2 * S1 - 
      2 * k_off12 * P1T1P2T2
    
    # Output calculation gives partial weight to ternary complex
    dOutput <- k_phos * (P1T1 + P1T1P2T2) - k_dephos * (P2T2 + P1T1P2T2)
    
    # Conservation check
    T_conserved <- T_total - (free_T + P1T1 + P2T2 + P1T1P2T2)
    if (abs(T_conserved) > 1e-6) warning("T not conserved in somewhat model at t=", t)
    
    list(c(dP1T1, dP2T2, dP1T1P2T2, dOutput))
  })
}

## 3. Parameter Sets --------------------------------------------------------------

base_params_competitive <- list(
  k_on1 = 0.1,
  k_off1 = 0.02,
  k_on2 = 0.1,
  k_off2 = 0.02,
  T_total = 10,
  P1_total = 100,
  P2_total = 100,
  k_phos = 1,
  k_dephos = 1
)

base_params_noncompetitive <- list(
  k_on1 = 0.1,
  k_off1 = 0.02,
  k_on2 = 0.1,
  k_off2 = 0.02,
  k_off12 = 0.02,
  T1_total = 5,  # Split T pool
  T2_total = 5,  # Split T pool
  P1_total = 100,
  P2_total = 100,
  k_phos = 1,
  k_dephos = 1
)

base_params_somewhat <- list(
  k_on1 = 0.1,
  k_off1 = 0.02,
  k_on2 = 0.1,
  k_off2 = 0.02,
  k_off12 = 0.02,
  T_total = 10,  # Single T pool
  P1_total = 100,
  P2_total = 100,
  k_phos = 1,
  k_dephos = 1
)

times <- seq(0, 1000, by = 1)

## 4. Simulation Function --------------------------------------------------------

simulate_scenario <- function(signal_type, model_func, model_label, params) {
  signal_trace <- generate_signal_trace(signal_type, times)
  
  # Determine initial state based on model
  if (identical(model_func, competitive_model)) {
    state <- c(P1T1 = 0, P2T2 = 0, Output = 0)
  } else {
    state <- c(P1T1 = 0, P2T2 = 0, P1T1P2T2 = 0, Output = 0)
  }
  
  out <- ode(
    y = state,
    times = times,
    func = model_func,
    parms = c(params, list(signal_trace = signal_trace)),
    method = "lsoda"
  )
  
  # Convert to data frame and add metadata
  result <- as.data.frame(out) %>%
    mutate(Model = model_label, 
           Signal = signal_type,
           free_T = if (model_label == "Competitive") {
             params$T_total - P1T1 - P2T2
           } else if (model_label == "Non-Competitive") {
             params$T1_total + params$T2_total - P1T1 - P2T2 - 2*P1T1P2T2
           } else {
             params$T_total - P1T1 - P2T2 - P1T1P2T2
           })
  
  return(result)
}

## 5. Run All Simulations --------------------------------------------------------

model_order <- c("Non-Competitive", "Somewhat Non-Competitive", "Competitive")
signal_order <- c("simultaneous", "alternating")

results <- map_dfr(signal_order, function(sig_type) {
  bind_rows(
    simulate_scenario(sig_type, competitive_model, "Competitive", base_params_competitive),
    simulate_scenario(sig_type, somewhat_noncompetitive_improved, "Somewhat Non-Competitive", base_params_somewhat),
    simulate_scenario(sig_type, non_competitive_model, "Non-Competitive", base_params_noncompetitive)
  )
}) %>%
  mutate(Model = factor(Model, levels = model_order),
         Signal = factor(Signal, levels = signal_order))

## 6. Plotting -------------------------------------------------------------------

# Prepare signal data for plotting
signal_data <- map_dfr(signal_order, function(sig_type) {
  generate_signal_trace(sig_type, times) %>%
    pivot_longer(cols = c(S1, S2), names_to = "Variable", values_to = "Value") %>%
    mutate(Model = "Signal", Signal = sig_type)
})

# Prepare model results for plotting
plot_data <- results %>%
  pivot_longer(cols = -c(time, Model, Signal, free_T), 
               names_to = "Variable", values_to = "Value") %>%
  mutate(Model = factor(Model, levels = model_order),
         Signal = factor(Signal, levels = signal_order))

# Create the plot
ggplot() +
  # Plot signal traces first
  geom_line(data = signal_data,
            aes(x = time, y = Value, color = Variable),
            linetype = "dashed", linewidth = 0.8) +
  
  # Plot model results
  geom_line(data = plot_data,
            aes(x = time, y = Value, color = Model, linetype = Model),
            linewidth = 1.2) +
  
  # Facet by variable and signal type
  facet_grid(Variable ~ Signal, scales = "free_y") +
  
  # Theme and labels
  theme_minimal(base_size = 14) +
  labs(x = "Time", y = "Concentration", 
       title = "Comparison of Binding Models with Signal Inputs",
       subtitle = "Dashed lines show input signals (S1, S2)") +
  
  # Color and line type scales
  scale_color_manual(
    values = c("S1" = "#E41A1C", "S2" = "#377EB8",
               "Competitive" = "#4DAF4A", 
               "Somewhat Non-Competitive" = "#FF7F00",
               "Non-Competitive" = "#984EA3")) +
  
  scale_linetype_manual(
    values = c("Competitive" = "solid",
               "Somewhat Non-Competitive" = "solid",
               "Non-Competitive" = "solid")) +
  
  # Theme adjustments
  theme(legend.position = "top",
        legend.box = "horizontal",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12))


## 7. Conservation Verification -------------------------------------------------

# Verify that total T is conserved in all models
conservation_check <- results %>%
  group_by(Model, Signal) %>%
  summarise(
    T_conservation_error = max(abs(
      if (first(Model) == "Competitive") {
        base_params_competitive$T_total - free_T - P1T1 - P2T2
      } else if (first(Model) == "Non-Competitive") {
        (base_params_noncompetitive$T1_total + base_params_noncompetitive$T2_total) - 
          (free_T + P1T1 + P2T2 + 2*P1T1P2T2)
      } else {
        base_params_somewhat$T_total - free_T - P1T1 - P2T2 - P1T1P2T2
      }
    )),
    .groups = "drop"
  )

print("Conservation Verification Results:")
print(conservation_check)