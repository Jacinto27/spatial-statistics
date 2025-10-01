#timeestimation

# -------- duration pretty-printer --------
.pretty_duration <- function(seconds) {
  if (seconds < 1) return(sprintf("%.3f s", seconds))
  mins <- floor(seconds / 60); secs <- seconds %% 60
  if (mins < 1) return(sprintf("%.2f s", secs))
  hrs <- floor(mins / 60); mins <- mins %% 60
  if (hrs < 1) return(sprintf("%dm %.0fs", mins, secs))
  sprintf("%dh %dm %.0fs", hrs, mins, secs)
}

# -------- simple per-"step" calibration on YOUR machine --------
# We define a "step" as a trivial scalar operation in a tight loop.
# This captures R loop overhead too (which is usually the dominating cost).
calibrate_step_time <- function(calibrate_iters = 2e6) {
  x <- 0L
  t <- system.time({
    for (i in 1:calibrate_iters) x <- x + 1L
  })["elapsed"]
  as.numeric(t) / calibrate_iters
}

# -------- 1) Generic estimator (any depth) --------
# loop_sizes: vector of sizes for each nested loop (e.g., c(n, m, k))
# steps_per_innermost_iter: how many "steps" happen inside the innermost body once
# t_step: optional seconds per step; if NULL we calibrate
# fixed_overhead: any constant setup/teardown cost you want to add (seconds)
estimate_runtime_generic <- function(loop_sizes,
                                     steps_per_innermost_iter,
                                     t_step = NULL,
                                     fixed_overhead = 0,
                                     calibrate_iters = 2e6) {
  if (is.null(t_step)) t_step <- calibrate_step_time(calibrate_iters)
  total_innermost_calls <- prod(as.numeric(loop_sizes))
  total_steps <- total_innermost_calls * steps_per_innermost_iter
  est_seconds <- total_steps * t_step + fixed_overhead
  list(
    total_innermost_calls = total_innermost_calls,
    steps_per_innermost_iter = steps_per_innermost_iter,
    total_steps = total_steps,
    t_step_seconds = t_step,
    est_seconds = est_seconds,
    pretty = .pretty_duration(est_seconds)
  )
}

# -------- 2) Two-loop convenience estimator --------
# Think:
# for (i in 1:n) {              # outer loop body cost = steps_outer
#   ... 4 "steps" ...
#   for (j in 1:m) {            # inner loop body cost = steps_inner
#     ... 85 "steps" ...
#   }
# }
estimate_runtime_2loops <- function(n, m,
                                    steps_outer,
                                    steps_inner,
                                    t_step = NULL,
                                    fixed_overhead = 0,
                                    calibrate_iters = 2e6) {
  if (is.null(t_step)) t_step <- calibrate_step_time(calibrate_iters)
  # Outer body runs n times; inner body runs n*m times
  total_steps <- n * steps_outer + (n * m) * steps_inner
  est_seconds <- total_steps * t_step + fixed_overhead
  list(
    n = n, m = m,
    steps_outer = steps_outer,
    steps_inner = steps_inner,
    total_steps = total_steps,
    t_step_seconds = t_step,
    est_seconds = est_seconds,
    pretty = .pretty_duration(est_seconds)
  )
}

# -------- optional: quick Big-O string helper (informational) --------
big_o_2loops <- function(n_sym = "n", m_sym = "m") sprintf("O(%s * %s)", n_sym, m_sym)




