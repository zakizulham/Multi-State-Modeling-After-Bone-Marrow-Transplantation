# ==============================================================================
# STOCHASTIC MODELING 1 - FINAL PROJECT
# Continuous-Time Markov Chain (CTMC) Analysis on Bone Marrow Transplant Data
# Department of Mathematics, Universitas Indonesia
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. SETUP & DEPENDENCIES
# ------------------------------------------------------------------------------
# To install dependencies, run: install.packages(c("mstate", "tidyverse", "Matrix"))

library(mstate)
library(tidyverse)
library(Matrix)

# Load the historical EBMT dataset
data("ebmt4")
raw_data <- ebmt4

# ------------------------------------------------------------------------------
# 2. TRANSITION MATRIX ARCHITECTURE (6-STATE CTMC)
# ------------------------------------------------------------------------------
# Designing the official 6-state multi-state structure from the JSS paper:
# State 1: Transplant (Tx)      -> Can transition to Rec, AE, Rel, or Death
# State 2: Recovery (Rec)       -> Can transition to Rec+AE, Relapse, or Death
# State 3: Adverse Event (AE)   -> Can transition to Rec+AE, Relapse, or Death
# State 4: Rec+AE               -> Can transition to Relapse or Death
# State 5: Relapse (Rel)        -> Absorbing State (Terminal due to data boundaries)
# State 6: Death (Absorbing)    -> Terminal State
# ------------------------------------------------------------------------------

t_mat <- transMat(
  x = list(
    c(2, 3, 5, 6), # Transitions from State 1 (Tx) to Rec(2), AE(3), Rel(5), or Death(6)
    c(4, 5, 6), # Transitions from State 2 (Rec) to Rec+AE(4), Rel(5), or Death(6)
    c(4, 5, 6), # Transitions from State 3 (AE) to Rec+AE(4), Rel(5), or Death(6)
    c(5, 6), # Transitions from State 4 (Rec+AE) to Rel(5) or Death(6)
    c(), # State  (Rel) is an Absorbing State (Terminal)
    c() # State 6 (Death) is an Absorbing State (Terminal)
  ),
  names = c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death")
)

print("Transition Architecture Matrix:")
print(t_mat)

# ------------------------------------------------------------------------------
# 3. DATA PREPARATION
# ------------------------------------------------------------------------------
ebmt_long <- msprep(
  time = c(NA, "rec", "ae", "recae", "rel", "srv"),
  status = c(NA, "rec.s", "ae.s", "recae.s", "rel.s", "srv.s"),
  data = ebmt4,
  trans = t_mat,
  keep = c("agecl", "match")
) |>
  as_tibble()

# ------------------------------------------------------------------------------
# 4. INTENSITY GENERATOR MATRIX (Q-MATRIX) ESTIMATION
# ------------------------------------------------------------------------------
# Estimating constant transition intensities (q_ij) using MLE:
# q_ij = total_events / total_time_at_risk
# ------------------------------------------------------------------------------

q_estimates <- ebmt_long |>
  group_by(from, to) |>
  summarise(
    total_events = sum(status),
    total_exposure = sum(time),
    q_intensity = total_events / total_exposure,
    .groups = "drop"
  )

print(q_estimates)

# Constructing the actual 6x6 infinitesimal generator matrix
Q_matrix <- matrix(
  0,
  nrow = 6,
  ncol = 6,
  dimnames = list(
    c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death"),
    c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death")
  )
)

# Populate off-diagonal intensities
for (i in 1:nrow(q_estimates)) {
  Q_matrix[q_estimates$from[i], q_estimates$to[i]] <- q_estimates$q_intensity[i]
}

# Compute diagonal elements: q_ii = -sum_{j \neq i} q_ij
diag(Q_matrix) <- -rowSums(Q_matrix)

print("Estimated Infinitesimal Generator Matrix (Q):")
print(Q_matrix)

# ------------------------------------------------------------------------------
# 5. EXPECTED SOJOURN TIMES ANALYSIS (WAKTU SINGGAH)
# ------------------------------------------------------------------------------
# In a time-homogeneous CTMC, the time spent in state i follows an Exponential
# distribution with rate v_i = -q_ii.
# The expected (mean) sojourn time is given by E[S_i] = 1 / -q_ii.
# ------------------------------------------------------------------------------

# Extract the transition rates out of each state (negative diagonal of Q)
rates_leaving <- -diag(Q_matrix)

# Compute the expected sojourn times (in days)
expected_sojourn <- 1 / rates_leaving

sojourn_analysis <- tibble(
  State = names(rates_leaving),
  Leaving_Rate = rates_leaving,
  Expected_Sojourn_Days = expected_sojourn
) |>
  mutate(
    Interpretation = if_else(
      is.infinite(Expected_Sojourn_Days),
      "Absorbing State (Permanent Stay)",
      paste(
        "Spends an average of",
        round(Expected_Sojourn_Days, 2),
        "days before transitioning"
      )
    )
  )

print("Expected Sojourn Times for Each Clinical State:")
print(sojourn_analysis)

# ------------------------------------------------------------------------------
# 6. CHAPMAN-KOLMOGOROV PROJECTION VIA MATRIX EXPONENTIAL
# ------------------------------------------------------------------------------
# Computing transient state probabilities P(t) = exp(Q * t) over a 2000-day horizon
# ------------------------------------------------------------------------------

time_horizon <- seq(0, 2000, by = 10)

patient_projection <- tibble(t = time_horizon) |>
  mutate(
    P_Tx = map_dbl(t, \(day) expm(Q_matrix * day)["Tx", "Tx"]),
    P_Rec = map_dbl(t, \(day) expm(Q_matrix * day)["Tx", "Rec"]),
    P_AE = map_dbl(t, \(day) expm(Q_matrix * day)["Tx", "AE"]),
    P_Rec_AE = map_dbl(t, \(day) expm(Q_matrix * day)["Tx", "Rec+AE"]),
    P_Relapse = map_dbl(t, \(day) expm(Q_matrix * day)["Tx", "Rel"]),
    P_Death = map_dbl(t, \(day) expm(Q_matrix * day)["Tx", "Death"])
  ) |>
  pivot_longer(
    cols = starts_with("P_"),
    names_to = "Status",
    values_to = "Probability"
  ) |>
  mutate(
    Status = recode_values(
      Status,
      "P_Tx" ~ "Initial State (Tx)",
      "P_Rec" ~ "Recovery Phase (Rec)",
      "P_AE" ~ "Adverse Event (AE)",
      "P_Rec_AE" ~ "Combined Phase (Rec+AE)",
      "P_Relapse" ~ "Relapse (Rel)",
      "P_Death" ~ "Deceased (Death)"
    ),
    Status = factor(
      Status,
      levels = c(
        "Initial State (Tx)",
        "Recovery Phase (Rec)",
        "Adverse Event (AE)",
        "Combined Phase (Rec+AE)",
        "Relapse (Rel)",
        "Deceased (Death)"
      )
    )
  )

# ------------------------------------------------------------------------------
# 7. EXPLICIT TRANSITION PROBABILITY MATRICES P(t) = exp(Q * t)
# ------------------------------------------------------------------------------
# Computing the Transition Probability Matrix for key clinical milestones.
# This directly applies the Chapman-Kolmogorov relation for time-homogeneous CTMC.
# ------------------------------------------------------------------------------

# Milestone 1: 100 Days Post-Transplantation (Standard Short-Term Landmark)
t_100 <- 100
P_matrix_100d <- as.matrix(expm(Q_matrix * t_100))
# Ensure row names and column names match the state space definition
dimnames(P_matrix_100d) <- dimnames(Q_matrix)

print("Transition Probability Matrix P(t) at t = 100 Days:")
print(round(P_matrix_100d, 4))


# Milestone 2: 365 Days Post-Transplantation (1-Year Long-Term Survival)
t_365 <- 365
P_matrix_1year <- as.matrix(expm(Q_matrix * t_365))
dimnames(P_matrix_1year) <- dimnames(Q_matrix)

print("Transition Probability Matrix P(t) at t = 365 Days (1 Year):")
print(round(P_matrix_1year, 4))

# ------------------------------------------------------------------------------
# 8. LIMITING BEHAVIOR & ABSORPTION PROBABILITIES
# ------------------------------------------------------------------------------
# As t approaches infinity, probabilities for all transient states limit to 0.
# We evaluate P(t) at an extreme horizon (t = 10000 days) to capture the exact
# absorption allocation between the competing absorbing states: Relapse and Death.
# ------------------------------------------------------------------------------

t_infinity <- 10000
P_matrix_infinity <- as.matrix(expm(Q_matrix * t_infinity))
dimnames(P_matrix_infinity) <- dimnames(Q_matrix)

print("Absorption / Limiting Probabilities Matrix (t -> Infinity):")
print(round(P_matrix_infinity, 5))

t_true_infinity <- 1000000
P_matrix_true_infinity <- as.matrix(expm(Q_matrix * t_true_infinity))
print(round(P_matrix_true_infinity, 5))

# ------------------------------------------------------------------------------
# 9. VISUALIZATION (ggplot2)
# ------------------------------------------------------------------------------
ctmc_plot <- ggplot(
  patient_projection,
  aes(x = t, y = Probability, color = Status)
) +
  geom_line(linewidth = 1.2, alpha = 0.85) +
  scale_color_manual(
    values = c(
      "Initial State (Tx)" = "#D4AF37",
      "Recovery Phase (Rec)" = "#2E7D32",
      "Adverse Event (AE)" = "#E65100",
      "Combined Phase (Rec+AE)" = "#303F9F",
      "Relapse (Rel)" = "#C62828",
      "Deceased (Death)" = "#212121"
    )
  ) +
  theme_minimal(base_size = 11, base_family = "sans") +
  labs(
    title = "Continuous-Time Markov Chain Dynamics in Bone Marrow Transplantation",
    subtitle = "Long-Term State Probability Trajectories Calibrated from EBMT Patient Data",
    x = "Days Post-Transplantation",
    y = "State Occupancy Probability P(t)",
    color = "Patient Status"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 13, margin = margin(b = 6)),
    plot.subtitle = element_text(
      size = 9.5,
      color = "gray30",
      margin = margin(b = 15)
    ),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 9.5),
    legend.text = element_text(size = 9),
    axis.title = element_text(face = "bold", size = 10.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray92")
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

print(ctmc_plot)

# ------------------------------------------------------------------------------
# 7. EXPORT
# ------------------------------------------------------------------------------
ggsave(
  filename = "figures/ctmc_trajectory.png",
  plot = ctmc_plot,
  width = 7,
  height = 5,
  dpi = 300
)
