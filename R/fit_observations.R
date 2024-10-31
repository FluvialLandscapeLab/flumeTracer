#' Fit a flume model to observed data
#'
#' Fits any model parameters to observed data, potentially allowing for fitting
#' a time offset between modeled and observed data.
#'
#' @param par A named vector of model parameters fit to the observed data
#'   (passed to \code{\link{optim}}; see that documentation for definition).  If
#'   one of the parameters is named "delta_t" an time offset error from t = 0
#'   will be included in the fit.
#' @param obs_times A vector of times where concentration values were observed.
#' @param obs_C_c A vetor of observed concentration values
#' @param m An environment describing a flume.  See \code{\link{createFlume}.
#' @param lower,upper A vector of lower and upper boundaries for par values.
#'   See "L-BFGS-B" method for \code{\link{optim}}, which is used for parameter
#'   estimates.
#'
#' @export
fitFlume <- function(par, obs_times, obs_C_c, m, lower, upper) {
  m$obs <- new.env()
  m$obs$times <- obs_times
  m$obs$C_c <- obs_C_c

  delta_t_idx <- which(names(par) == "delta_t")
  if(length(delta_t_idx) == 0) {
    result <-
      optim(par = par, simulationError, obj_fun = staticObjective, m = m,
            method = "L-BFGS-B", lower = lower, upper = upper)
  } else {
    delta_t_lower <- lower[delta_t_idx]
    delta_t_upper <- upper[delta_t_idx]
    par <- par[-delta_t_idx]
    lower <- lower[-delta_t_idx]
    upper <- upper[-delta_t_idx]
    result <- optim(par = par, simulationError, obj_fun = delta_tObjective, m = m,
                    delta_t_lower = delta_t_lower, delta_t_upper = delta_t_upper,
                    method = "L-BFGS-B", lower = lower, upper = upper)
  }
  result
}

# delta_t values are ignored with obj_fun is staticObjective.
simulationError <- function(par, obj_fun, m, delta_t_lower=NULL, delta_t_upper=NULL) {
  m_run <- as.environment(as.list(m))
  mapply(assign, x = names(par), value = par, MoreArgs = list(envir = m_run))
  simulateFlumeConcentration(m_run, debug = F)
  obj_fun(
    model_times = m_run$times,
    model_C_c = m_run$C_c,
    obs_times = m$obs$times,
    obs_C_c = m$obs$C_c,
    m  = m, delta_t_lower = delta_t_lower, delta_t_upper = delta_t_upper)
}

# m and delta_t parameters are ignored; must have same signature as delta_t_Objective
staticObjective <- function(delta_t = 0, model_times, model_C_c, obs_times, obs_C_c, m,
                            delta_t_lower = NULL, delta_t_upper = NULL) {
  adj_model_C_c <- approx(model_times, model_C_c, obs_times + delta_t)$y
  sqrt(mean((adj_model_C_c - obs_C_c)^2, na.rm = T))
}

delta_tObjective <- function(model_times, model_C_c, obs_times, obs_C_c, m,
                             delta_t_lower, delta_t_upper) {
  result <-
    optimize(
      f = staticObjective, interval = c(delta_t_lower, delta_t_upper),
      model_times = model_times, model_C_c = model_C_c,
      obs_times = obs_times, obs_C_c = obs_C_c)
  m$delta_t <- result$minimum
  result$objective
}

