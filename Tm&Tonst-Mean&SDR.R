# ---- Melt curves with unified Tonst (amp+slope + guardrails + pattern overrides) ----
# ---- PLUS replicate harmonization so nearly-identical wells share the same Tonst ----
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(purrr)
})

# >>> EDIT THIS if your file name/path differs
file <- "Melt Curve RFU Results_FRET.csv"

# -------------------------- Tunables ---------------------------------------
# Smoothing & prefiltering
loess_span      <- 0.20
use_runmed      <- TRUE
runmed_k        <- 5

# Amplitude (robust onset)
k_sigma_amp     <- 3.0
frac_fallback   <- 0.10
baseline_frac   <- 0.20
frac_floor_amp  <- 0.02
consec_pts      <- 2

# Slope-based criterion
k_sigma_slope   <- 3.0
slope_floor     <- 0.02
slope_win       <- 3
min_total_frac  <- 0.05

# Global guardrails for Tonst
early_cap_frac  <- 0.25
cap_delta_deg   <- 6

# --- Per-sample Tonst guardrail overrides (regex patterns for robustness) ---
caps_override_patterns <- list(
  "A0\\.1.*5" = list(early_cap_frac = 0.30, cap_delta_deg = 8),
  "A0\\.5.*9" = list(early_cap_frac = 0.30, cap_delta_deg = 8)
)

# Replicate harmonization (make near-identical replicates share Tonst)
replicate_locking       <- TRUE
rep_regex               <- "^(.*?)-\\d+$"  # capture base name before trailing -<rep>
sim_corr_threshold      <- 0.999
sim_rmse_frac_threshold <- 0.04            # 4% of dynamic range
use_pooled_baseline     <- TRUE

# Plot export sizes
plot_w_all <- 11.5; plot_h_all <- 6.7
plot_w_fac <- 12.5; plot_h_fac <- 8.5

# ---------------------------- Load & reshape -------------------------------
raw <- readr::read_csv(file, guess_max = 200000)
stopifnot(ncol(raw) >= 2)

# Auto-detect temperature column
temp_guess <- names(raw)[str_detect(tolower(names(raw)),
                                    "temp|°c|celsius|well\\s*temp|temperature")]
temp_col <- if (length(temp_guess)) temp_guess[1] else names(raw)[1]
message(sprintf("Using '%s' as temperature column.", temp_col))

# Coerce to numeric & keep useful columns
df_num <- raw %>%
  rename(Temperature_raw = !!temp_col) %>%
  mutate(Temperature = suppressWarnings(as.numeric(Temperature_raw))) %>%
  select(-Temperature_raw) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))) %>%
  select(Temperature, where(~ is.numeric(.x) && sum(!is.na(.x)) > 0))

if (!"Temperature" %in% names(df_num)) stop("Could not parse a numeric Temperature column.")
if (ncol(df_num) < 2) stop("No usable sample columns found (numeric with at least one non-NA value).")

# Long format and de-duplicate identical Temperatures per sample (average)
long <- df_num %>%
  pivot_longer(-Temperature, names_to = "Sample", values_to = "Signal") %>%
  filter(!is.na(Temperature), !is.na(Signal)) %>%
  group_by(Sample, Temperature) %>%
  summarise(Signal = mean(Signal, na.rm = TRUE), .groups = "drop") %>%
  arrange(Sample, Temperature)

# ------------------------------- Helpers -----------------------------------
lin_interp_x <- function(x1, y1, x2, y2, y_target) {
  if (!is.finite(x1) || !is.finite(x2) || !is.finite(y1) || !is.finite(y2) || (y2 - y1) == 0)
    return(NA_real_)
  x1 + (y_target - y1) * (x2 - x1) / (y2 - y1)
}

prefilter <- function(y, k = 5) {
  k <- max(3, as.integer(k) | 1L) # force odd >=3
  stats::runmed(y, k = k, endrule = "median")
}

robust_loess <- function(T, Y, span, do_runmed = TRUE, k_med = 5) {
  y_in <- if (do_runmed) prefilter(Y, k = k_med) else Y
  fit <- try(stats::loess(y_in ~ T, span = span, degree = 2,
                          family = "symmetric", surface = "direct"),
             silent = TRUE)
  Y_s <- if (inherits(fit, "try-error")) y_in else as.numeric(predict(fit, T))
  bad <- !is.finite(Y_s)
  if (any(bad)) Y_s[bad] <- y_in[bad]
  Y_s
}

moving_median <- function(x, k = 3) {
  k <- max(3, as.integer(k) | 1L)
  n <- length(x)
  out <- rep(NA_real_, n)
  h <- (k - 1) / 2
  for (i in seq_len(n)) {
    lo <- max(1, i - h); hi <- min(n, i + h)
    out[i] <- stats::median(x[lo:hi], na.rm = TRUE)
  }
  out
}

decide_mode <- function(T, Y, base_frac = 0.20) {
  o <- order(T); T <- T[o]; Y <- Y[o]
  n <- length(T)
  if (n < 6) return("unknown")
  cut_idx <- max(3, floor(base_frac * n))
  cut_idx <- min(cut_idx, n - 3)
  if (cut_idx < 3) return("unknown")
  base_med <- median(Y[1:cut_idx], na.rm = TRUE)
  up   <- max(Y, na.rm = TRUE) - base_med
  down <- base_med - min(Y, na.rm = TRUE)
  if (!is.finite(up) || !is.finite(down)) return("unknown")
  if (up >= down) "rising" else "minima"
}

# ---- Tm routines ----------------------------------------------------------
half_rise_tm <- function(T, Y) {
  o <- order(T); T <- T[o]; Y <- Y[o]
  if (length(T) < 2) return(list(Tm = NA_real_, Tm_kind = "half_rise",
                                 Y_marker = NA_real_, Y_base = NA_real_, Y_ext = NA_real_))
  j_max <- which.max(Y)
  T_pre <- T[seq_len(j_max)]
  Y_pre <- Y[seq_len(j_max)]
  if (length(T_pre) < 2) return(list(Tm = NA_real_, Tm_kind = "half_rise",
                                     Y_marker = NA_real_, Y_base = NA_real_, Y_ext = NA_real_))
  y_base <- min(Y_pre, na.rm = TRUE)
  y_max  <- max(Y_pre, na.rm = TRUE)
  if (!is.finite(y_base) || !is.finite(y_max) || y_max <= y_base)
    return(list(Tm = NA_real_, Tm_kind = "half_rise",
                Y_marker = NA_real_, Y_base = y_base, Y_ext = y_max))
  target <- y_base + 0.5 * (y_max - y_base)
  below_idx <- which(Y_pre <  target)
  above_idx <- which(Y_pre >= target)
  if (length(above_idx) == 0)
    return(list(Tm = NA_real_, Tm_kind = "half_rise",
                Y_marker = NA_real_, Y_base = y_base, Y_ext = y_max))
  j <- min(above_idx)
  if (length(below_idx) == 0 || max(below_idx) >= j)
    return(list(Tm = T_pre[j], Tm_kind = "half_rise",
                Y_marker = Y_pre[j], Y_base = y_base, Y_ext = y_max))
  i <- max(below_idx)
  T_half <- lin_interp_x(T_pre[i], Y_pre[i], T_pre[j], Y_pre[j], target)
  list(Tm = as.numeric(T_half), Tm_kind = "half_rise",
       Y_marker = as.numeric(target), Y_base = y_base, Y_ext = y_max)
}

tm_min <- function(T, Y) {
  o <- order(T); T <- T[o]; Y <- Y[o]
  if (!length(T)) return(list(Tm = NA_real_, Tm_kind = "minima", Y_marker = NA_real_))
  i <- which.min(Y)
  list(Tm = as.numeric(T[i]), Tm_kind = "minima", Y_marker = as.numeric(Y[i]))
}

# ---- Unified robust Tonst (amplitude + slope + guardrails + final fallback) ----
# dir = +1 (rising), dir = -1 (minima)
tonst_core <- function(T, Y, dir = +1,
                       span = 0.2,
                       # amplitude params
                       k_sigma_amp = 3, frac_fb = 0.10, base_frac = 0.20, frac_floor_amp = 0.02,
                       # slope params
                       k_sigma_slope = 3, slope_floor = 0.02, slope_win = 3,
                       # persistence & minimal effect size
                       consec = 2, min_total_frac = 0.05,
                       # guardrails
                       early_cap_frac = NA_real_,   # e.g., 0.25 caps Tonst at 25% of change
                       cap_delta_deg  = NA_real_,   # e.g., 6 caps Tonst ≤ T_ext - 6°C
                       # smoothing
                       do_runmed = TRUE, runmed_k = 5) {
  stopifnot(dir %in% c(+1, -1))
  o <- order(T); T <- T[o]; Y <- Y[o]
  n <- length(T)
  if (n < 7) return(list(Tonst = NA_real_, Y_on = NA_real_, method = "insufficient_points",
                         AmpThr = NA, SlopeThr = NA, BaseMed = NA, BaseMAD = NA))
  
  Y_s <- robust_loess(T, Y, span = span, do_runmed = do_runmed, k_med = runmed_k)
  
  idx_ext <- if (dir == +1) which.max(Y_s) else which.min(Y_s)
  if (idx_ext <= 4) return(list(Tonst = NA_real_, Y_on = NA_real_, method = "no_pretransition",
                                AmpThr = NA, SlopeThr = NA, BaseMed = NA, BaseMAD = NA))
  T_ext <- T[idx_ext]
  
  cut_idx <- max(3, floor(base_frac * n)); cut_idx <- min(cut_idx, idx_ext - 3)
  if (cut_idx < 3) return(list(Tonst = NA_real_, Y_on = NA_real_, method = "baseline_too_short",
                               AmpThr = NA, SlopeThr = NA, BaseMed = NA, BaseMAD = NA))
  base_med <- median(Y_s[1:cut_idx], na.rm = TRUE)
  base_mad <- stats::mad(Y_s[1:cut_idx], center = base_med, constant = 1.4826, na.rm = TRUE); if (!is.finite(base_mad)) base_mad <- 0
  
  y_ext <- Y_s[idx_ext]
  total_change <- if (dir == +1) (y_ext - base_med) else (base_med - y_ext)
  total_change <- max(0, total_change)
  
  if (!is.finite(total_change) || total_change < min_total_frac * max(abs(Y_s), na.rm = TRUE))
    return(list(Tonst = NA_real_, Y_on = NA_real_, method = "too_small_change",
                AmpThr = NA, SlopeThr = NA, BaseMed = base_med, BaseMAD = base_mad))
  
  amp_step <- max(k_sigma_amp * base_mad, frac_floor_amp * total_change)
  AmpThr <- if (dir == +1) base_med + amp_step else base_med - amp_step
  
  dT <- diff(T); dY <- diff(Y_s); slope <- dY / dT
  slope_s <- moving_median(slope, k = max(3, as.integer(slope_win) | 1L))
  s_base <- slope_s[1:max(1, cut_idx - 1)]
  s_med  <- median(s_base, na.rm = TRUE)
  s_mad  <- stats::mad(s_base, center = s_med, constant = 1.4826, na.rm = TRUE); if (!is.finite(s_mad)) s_mad <- 0
  max_abs_slope <- max(abs(slope_s), na.rm = TRUE)
  slope_step <- max(k_sigma_slope * s_mad, slope_floor * max_abs_slope)
  SlopeThr <- if (dir == +1) s_med + slope_step else s_med - slope_step
  
  y_pre <- Y_s[1:idx_ext]
  s_pre <- slope_s[1:max(1, idx_ext - 1)]
  
  cross_amp   <- if (dir == +1) (y_pre >= AmpThr)   else (y_pre <= AmpThr)
  cross_slope <- if (dir == +1) (s_pre >= SlopeThr) else (s_pre <= SlopeThr)
  
  sustained_idx <- function(v, need) {
    if (!any(v)) return(NA_integer_)
    r <- rle(v); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
    idx <- which(r$values & (r$lengths >= need))
    if (length(idx) == 0) return(NA_integer_)
    starts[idx[1]]
  }
  a_start <- sustained_idx(cross_amp, consec)
  s_start <- sustained_idx(cross_slope, consec)
  
  slope_on_points <- rep(FALSE, length(y_pre))
  if (length(s_pre)) {
    for (k in seq_along(s_pre)) if (isTRUE(cross_slope[k])) slope_on_points[k+1] <- TRUE
    if (consec > 1 && length(slope_on_points) >= consec) {
      roll_ok <- rep(FALSE, length(slope_on_points))
      for (i in seq_len(length(slope_on_points) - consec + 1)) {
        if (all(slope_on_points[i:(i+consec-1)])) roll_ok[i:(i+consec-1)] <- TRUE
      }
      slope_on_points <- roll_ok
    }
  }
  
  both <- cross_amp & slope_on_points
  Tonst_cand <- NA_real_; Y_on <- NA_real_; method <- NA_character_
  
  j <- which(both)[1]
  if (is.finite(j)) {
    if (j == 1) { Tonst_cand <- T[1]; Y_on <- y_pre[1]; method <- "amp+slope_edge" }
    else { i <- j - 1; Tonst_cand <- lin_interp_x(T[i], y_pre[i], T[i+1], y_pre[i+1], AmpThr); Y_on <- AmpThr; method <- "amp+slope" }
  } else if (is.finite(a_start)) {
    j <- a_start
    if (j == 1) { Tonst_cand <- T[1]; Y_on <- y_pre[1]; method <- "amp_only_edge" }
    else { i <- j - 1; Tonst_cand <- lin_interp_x(T[i], y_pre[i], T[i+1], y_pre[i+1], AmpThr); Y_on <- AmpThr; method <- "amp_only" }
  } else if (is.finite(s_start)) {
    j <- s_start + 1; Tonst_cand <- T[j]; Y_on <- y_pre[j]; method <- "slope_only"
  } else {
    frac_target <- if (dir == +1) base_med + frac_fallback * total_change else base_med - frac_fallback * total_change
    cross_frac <- if (dir == +1) (y_pre >= frac_target) else (y_pre <= frac_target)
    j <- which(cross_frac)[1]
    if (is.finite(j)) {
      if (j == 1) { Tonst_cand <- T[1]; Y_on <- y_pre[1]; method <- sprintf("frac_only(%.0f%%)-edge", 100*frac_fallback) }
      else { i <- j - 1; Tonst_cand <- lin_interp_x(T[i], y_pre[i], T[i+1], y_pre[i+1], frac_target); Y_on <- frac_target; method <- sprintf("frac_only(%.0f%%)", 100*frac_fallback) }
    } else {
      if (is.finite(early_cap_frac) && early_cap_frac > 0 && early_cap_frac < 1) {
        cap_target <- if (dir == +1) base_med + early_cap_frac * total_change
        else            base_med - early_cap_frac * total_change
        cross_cap <- if (dir == +1) (y_pre >= cap_target) else (y_pre <= cap_target)
        j2 <- which(cross_cap)[1]
        if (is.finite(j2)) {
          T_cap <- if (j2 == 1) T[1] else {
            i2 <- j2 - 1
            lin_interp_x(T[i2], y_pre[i2], T[i2+1], y_pre[i2+1], cap_target)
          }
          return(list(Tonst = as.numeric(T_cap), Y_on = as.numeric(cap_target),
                      method = sprintf("cap_fallback(%.0f%%)", 100*early_cap_frac),
                      AmpThr = AmpThr, SlopeThr = SlopeThr, BaseMed = base_med, BaseMAD = base_mad))
        }
      }
      return(list(Tonst = NA_real_, Y_on = NA_real_, method = "no_crossing",
                  AmpThr = AmpThr, SlopeThr = SlopeThr, BaseMed = base_med, BaseMAD = base_mad))
    }
  }
  
  Tonst_final <- Tonst_cand
  
  # Early-cap by fraction
  if (is.finite(early_cap_frac) && early_cap_frac > 0 && early_cap_frac < 1) {
    frac_target <- if (dir == +1) base_med + early_cap_frac * total_change
    else             base_med - early_cap_frac * total_change
    cross_frac2 <- if (dir == +1) (y_pre >= frac_target) else (y_pre <= frac_target)
    j2 <- which(cross_frac2)[1]
    if (is.finite(j2)) {
      T_frac <- if (j2 == 1) T[1] else {
        i2 <- j2 - 1
        lin_interp_x(T[i2], y_pre[i2], T[i2+1], y_pre[i2+1], frac_target)
      }
      if (is.finite(T_frac)) Tonst_final <- min(Tonst_final, T_frac, na.rm = TRUE)
    }
  }
  
  # Cap by minimum gap from extremum
  if (is.finite(cap_delta_deg) && cap_delta_deg > 0) {
    Tonst_final <- min(Tonst_final, T_ext - cap_delta_deg, na.rm = TRUE)
  }
  
  list(Tonst = as.numeric(Tonst_final), Y_on = as.numeric(Y_on), method = method,
       AmpThr = AmpThr, SlopeThr = SlopeThr, BaseMed = base_med, BaseMAD = base_mad)
}

tonst_rising <- function(T, Y, ...)  tonst_core(T, Y, dir = +1, ...)
tonst_minima <- function(T, Y, ...)  tonst_core(T, Y, dir = -1, ...)

# --- Pattern-based overrides helper ---
get_caps <- function(sample, default_early, default_delta, override_patterns) {
  e <- default_early; d <- default_delta
  for (pat in names(override_patterns)) {
    if (isTRUE(grepl(pat, sample))) {
      o <- override_patterns[[pat]]
      if (!is.null(o$early_cap_frac)) e <- o$early_cap_frac
      if (!is.null(o$cap_delta_deg))  d <- o$cap_delta_deg
      break
    }
  }
  list(early = e, delta = d)
}

# ----------------------- Compute per-sample outputs ------------------------
results <- long %>%
  group_by(Sample) %>%
  group_modify(~{
    df <- .x
    T <- df$Temperature
    Y <- df$Signal
    sm <- decide_mode(T, Y, base_frac = baseline_frac)
    
    # Per-sample caps via regex patterns
    smpl <- .y$Sample[[1]]
    caps <- get_caps(smpl, early_cap_frac, cap_delta_deg, caps_override_patterns)
    ecap <- caps$early
    dcap <- caps$delta
    
    if (identical(sm, "rising")) {
      tm <- half_rise_tm(T, Y)
      on <- tonst_rising(
        T, Y,
        span = loess_span,
        k_sigma_amp = k_sigma_amp, frac_fb = frac_fallback, base_frac = baseline_frac, frac_floor_amp = frac_floor_amp,
        k_sigma_slope = k_sigma_slope, slope_floor = slope_floor, slope_win = slope_win,
        consec = consec_pts, min_total_frac = min_total_frac,
        early_cap_frac = ecap, cap_delta_deg = dcap,
        do_runmed = use_runmed, runmed_k = runmed_k
      )
    } else if (identical(sm, "minima")) {
      tm <- tm_min(T, Y)
      on <- tonst_minima(
        T, Y,
        span = loess_span,
        k_sigma_amp = k_sigma_amp, frac_fb = frac_fallback, base_frac = baseline_frac, frac_floor_amp = frac_floor_amp,
        k_sigma_slope = k_sigma_slope, slope_floor = slope_floor, slope_win = slope_win,
        consec = consec_pts, min_total_frac = min_total_frac,
        early_cap_frac = ecap, cap_delta_deg = dcap,
        do_runmed = use_runmed, runmed_k = runmed_k
      )
    } else {
      tm <- list(Tm = NA_real_, Tm_kind = "unknown", Y_marker = NA_real_)
      on <- list(Tonst = NA_real_, Y_on = NA_real_, method = "unknown",
                 AmpThr = NA, SlopeThr = NA, BaseMed = NA, BaseMAD = NA)
    }
    
    tibble(
      Mode     = sm,
      Tm       = tm$Tm,
      Tm_kind  = tm$Tm_kind,
      Y_Tm     = tm$Y_marker,
      Tonst    = on$Tonst,
      Y_on     = on$Y_on,
      OnMethod = on$method,
      BaseMed  = on$BaseMed,
      BaseMAD  = on$BaseMAD,
      AmpThr   = on$AmpThr,
      SlopeThr = on$SlopeThr
    )
  }) %>%
  ungroup()

# ----------------- Replicate harmonization (force same Tonst if identical) -----------------
if (replicate_locking) {
  res2 <- results %>%
    mutate(RepGroup = ifelse(grepl(rep_regex, Sample),
                             sub(rep_regex, "\\1", Sample),
                             Sample))
  
  long_with_group <- long %>%
    mutate(RepGroup = ifelse(grepl(rep_regex, Sample),
                             sub(rep_regex, "\\1", Sample),
                             Sample))
  
  pool_tonst_for_group <- function(grp_name) {
    df_grp <- long_with_group %>% filter(RepGroup == grp_name)
    smpls  <- unique(df_grp$Sample)
    
    pooled <- df_grp %>%
      group_by(Temperature) %>%
      summarise(Signal = mean(Signal, na.rm = TRUE), .groups = "drop") %>%
      arrange(Temperature)
    
    Yp_s <- robust_loess(pooled$Temperature, pooled$Signal,
                         span = loess_span, do_runmed = use_runmed, k_med = runmed_k)
    
    mode_pooled <- decide_mode(pooled$Temperature, Yp_s, base_frac = baseline_frac)
    
    if (identical(mode_pooled, "rising")) {
      on_pooled <- tonst_rising(
        pooled$Temperature, Yp_s,
        span = loess_span,
        k_sigma_amp = k_sigma_amp, frac_fb = frac_fallback, base_frac = baseline_frac, frac_floor_amp = frac_floor_amp,
        k_sigma_slope = k_sigma_slope, slope_floor = slope_floor, slope_win = slope_win,
        consec = consec_pts, min_total_frac = min_total_frac,
        early_cap_frac = if (use_pooled_baseline) early_cap_frac else NA_real_,
        cap_delta_deg = cap_delta_deg,
        do_runmed = FALSE, runmed_k = runmed_k
      )
    } else if (identical(mode_pooled, "minima")) {
      on_pooled <- tonst_minima(
        pooled$Temperature, Yp_s,
        span = loess_span,
        k_sigma_amp = k_sigma_amp, frac_fb = frac_fallback, base_frac = baseline_frac, frac_floor_amp = frac_floor_amp,
        k_sigma_slope = k_sigma_slope, slope_floor = slope_floor, slope_win = slope_win,
        consec = consec_pts, min_total_frac = min_total_frac,
        early_cap_frac = if (use_pooled_baseline) early_cap_frac else NA_real_,
        cap_delta_deg = cap_delta_deg,
        do_runmed = FALSE, runmed_k = runmed_k
      )
    } else {
      on_pooled <- list(Tonst = NA_real_)
    }
    
    Tonst_pool <- on_pooled$Tonst
    
    similar_vec <- map_lgl(smpls, function(smpl) {
      df_s <- df_grp %>% filter(Sample == smpl) %>% arrange(Temperature)
      if (nrow(df_s) < 5) return(FALSE)
      Ys_s <- robust_loess(df_s$Temperature, df_s$Signal,
                           span = loess_span, do_runmed = use_runmed, k_med = runmed_k)
      joined <- inner_join(
        tibble(T = df_s$Temperature, Ys = Ys_s),
        tibble(T = pooled$Temperature, Yp = Yp_s),
        by = "T"
      )
      if (nrow(joined) < 5) return(FALSE)
      rng <- diff(range(joined$Yp, na.rm = TRUE)); rng <- if (is.finite(rng) && rng > 0) rng else 1
      rmse <- sqrt(mean((joined$Ys - joined$Yp)^2, na.rm = TRUE)) / rng
      cr   <- suppressWarnings(cor(joined$Ys, joined$Yp, use = "complete.obs"))
      isTRUE(cr >= sim_corr_threshold && rmse <= sim_rmse_frac_threshold)
    })
    
    if (length(similar_vec) >= 2 && all(similar_vec) && is.finite(Tonst_pool)) {
      tibble(RepGroup = grp_name, Tonst_pooled = Tonst_pool, LockApplied = TRUE)
    } else {
      tibble(RepGroup = grp_name, Tonst_pooled = NA_real_, LockApplied = FALSE)
    }
  }
  
  pool_tbl <- unique(long_with_group$RepGroup) %>%
    map_dfr(pool_tonst_for_group)
  
  results <- res2 %>%
    left_join(pool_tbl, by = "RepGroup") %>%
    mutate(Tonst = ifelse(LockApplied, Tonst_pooled, Tonst)) %>%
    select(-Tonst_pooled, -LockApplied, -RepGroup)
}

# ---- NEW: Replicate mean & SD (Tm and Tonst) -------------------------------
# Extract replicate base & number (assumes names like Base-1, Base-2, etc.)
extract_base <- function(x) ifelse(grepl(rep_regex, x), sub(rep_regex, "\\1", x), x)
extract_repnum <- function(x) { m <- regexpr("(\\d+)$", x, perl = TRUE); ifelse(m > 0, as.integer(regmatches(x, m)), NA_integer_) }

results_rep <- results %>%
  mutate(BaseName = extract_base(Sample),
         Rep      = extract_repnum(Sample))

# Wide table with Tm1, Tm2, Tonst1, Tonst2 (if present)
replicates_wide <- results_rep %>%
  filter(!is.na(Rep)) %>%
  select(BaseName, Rep, Tm, Tonst) %>%
  pivot_wider(names_from = Rep, values_from = c(Tm, Tonst), names_glue = "{.value}{Rep}") %>%
  arrange(BaseName)

# Mean & SD per BaseName (over all replicates present)
replicates_mean_sd <- results_rep %>%
  group_by(BaseName) %>%
  summarise(
    n_reps     = sum(!is.na(Rep)),
    Tm_mean    = if (n_reps > 0) mean(Tm, na.rm = TRUE) else NA_real_,
    Tm_sd      = if (n_reps > 1) sd(Tm, na.rm = TRUE)   else NA_real_,
    Tonst_mean = if (n_reps > 0) mean(Tonst, na.rm = TRUE) else NA_real_,
    Tonst_sd   = if (n_reps > 1) sd(Tonst, na.rm = TRUE)   else NA_real_,
    .groups = "drop"
  ) %>%
  arrange(BaseName)

# Save tables
readr::write_csv(results,              "Tm_Tonst_autodetect_withReplicateLock.csv")
readr::write_csv(replicates_wide,      "Replicates_wide.csv")
readr::write_csv(replicates_mean_sd,   "Replicates_mean_sd.csv")

# ------------------------------ Plotting -----------------------------------
tm_points <- results %>% select(Sample, Mode, Tm, Y_Tm)
on_points <- results %>% select(Sample, Mode, Tonst, Y_on)

p_all <- ggplot(long, aes(x = Temperature, y = Signal, group = Sample)) +
  geom_line(alpha = 0.55) +
  geom_point(data = tm_points, aes(x = Tm, y = Y_Tm, shape = Mode), size = 1.9) +
  geom_vline(data = tm_points, aes(xintercept = Tm, linetype = Mode)) +
  geom_point(data = on_points, aes(x = Tonst, y = Y_on, shape = Mode), size = 1.9) +
  geom_vline(data = on_points, aes(xintercept = Tonst, linetype = Mode), alpha = 0.7) +
  scale_shape_manual(values = c(rising = 16, minima = 17, unknown = 4)) +
  scale_linetype_manual(values = c(rising = "dashed", minima = "dotted", unknown = "longdash")) +
  labs(title = "Melt curves with unified Tonst (dual criteria + guardrails + replicate locking)",
       x = "Temperature (°C)", y = "Signal",
       caption = sprintf("LOESS=%.2f, runmed=%s(k=%d); Amp: k=%.1f floor=%.0f%%; Slope: k=%.1f floor=%.0f%%; consec=%d; early_cap=%.0f%%; Δcap=%g°C; lock=%s (corr≥%.3f & RMSE≤%.1f%% dyn.range)",
                         loess_span, use_runmed, runmed_k,
                         k_sigma_amp, 100*frac_floor_amp,
                         k_sigma_slope, 100*slope_floor, consec_pts,
                         100*early_cap_frac, cap_delta_deg,
                         replicate_locking, sim_corr_threshold, 100*sim_rmse_fr*_
                         