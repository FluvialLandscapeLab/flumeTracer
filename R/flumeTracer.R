#' Model a flume tracer release.
#'
#' Simulates conservative tracer concentration in an annular flume following a
#' slug release.
#'
#' Uses and implicit solver to estimate flume channel water concentration at
#' each time step of a simulation.
#'
#' When the function if executed, the period of record for channel concentration
#' is built, one value at a time, via iteration.  The implicit solver takes
#' advantage of the fact that the mean concentration of water in the flume at
#' the end of the salt release is known:
#'
#' \code{C_final = (C_prerelease * V + C_slug * V_slug)/(V + V_slug) }
#'
#' where V is total water volume in flume.
#'
#' Alternatively, if the instantaneous flume concentration at the time of slug
#' injection can be estimated:
#'
#' \code{C_final = (C_prerelease * V_h + C_postrelease * V_c) / V}
#'
#' where V_h is the water volume in the hyporheic zone and V_c is the water
#' volume in the channel (surface water).
#'
#' These estimates can be confirmed with the observed concentrations at the end
#' of the experimental release.
#'
#' At all times during the release, the following equality should hold:
#'
#' \code{C_final = (C_h_t * V_h + C_c_t * V_c) / V}
#'
#' where "_t" designates a value at time (t) since the release. Importantly,
#' \code{C_c} at the end of the current time step affects the history of
#' interpolated \code{C_c} over the course of the time step and therefore
#' influences C_h at the end of the time step. Because C_c_t influences C_h_t,
#' there is only one value of C_c_t that will satisfy the equality.   Therefore,
#' at each time step of the model, the implicit solver finds the value of C_c
#' that minimizes the squared difference between each side of the equality. See
#' \code{\link{optimizationError}} for more details.
#' @param m An environment (often created with \code{\link{createFlume}}) in which the follow variables have been defined:
#'   \itemize{
#'
#'   \item{\code{timestep} the time period each iteration of the model
#'   represents}
#'
#'   \item{\code{duration} the total duration of the simulation}
#'
#'   \item{\code{tau_0} the minimum water age in the hyporheic zone}
#'
#'   \item{\code{tau_n} the maximum water age in the hyporheic zone}
#'
#'   \item{\code{shape} either "powerLaw" or "exponent" depending on the desired
#'   shape of the washout function.}
#'
#'   \item{\code{alpha} or \code{sigma} the exponent for shape "powerLaw"
#'   (\code{alpha} -- use a positive value. This value is negated within the
#'   model code...) or the decay rate (\code{sigma}) for shape "exponent".}
#'
#'   \item{\code{C_c} the initial channel (surface water) concentration
#'   immediately following the slug release}
#'
#'   \item{\code{C_h} the initial hyporheic water concentration (equal to
#'   concentration prior to the release)}
#'
#'   \item{\code{V} the total water volume in the flume}
#'
#'   \item{\code{V_h} the water volume within the hyporheic zone}
#'
#'   \item{\code{nSubdiv} the number of manual subdivisions used for
#'   integration.  Use 3 as an initial value.  Try increasing this value if
#'   "divergent integration" errors arise, or to improve the accuracy of the
#'   integration (at the cost of increased computational time).  A value above
#'   about 15 is not likely to improve things.}
#'
#'   }
#' @param debug Control debugging output. The model relies substantially on
#'   numerical integration.  If \code{debug} is set to TRUE, the model will
#'   record the integration quality information for every numerical integration
#'   to a variable called \code{debug} in the model environment passed to
#'   \code{m}.  This will slow the model execution somewhat and is extremely
#'   verbose.
#' @return TRUE if executed successfully.
#'
#'   As a side effect, the following values are added to the environment
#'   \code{m}:
#'
#'   \itemize{
#'
#'   \item{\code{nIterations} The number of iterations executed by the model.}
#'
#'   \item{\code{times} The time values at the end of each iteration.  This
#'   variable is useful as an x-axis for plotting time series results.}
#'
#'   \item{\code{nTimes} The length of any time series produced by the model
#'   (nIterations + 1, because initial conditions are include in all time series
#'   variables.}
#'
#'   \item{\code{V_c} The volume of surface water (difference between V and
#'   V_h)}
#'
#'   \item{\code{C_final} The expected final concentration in the flume;
#'   calculated as \code{(C_c*V_c + C_h*V_h) / V}}
#'
#'   \item{\code{debug}} Information about the quality of numerical integrations
#'   used in the model.  Only created if \code{debug} parameter is set to TRUE.
#'
#'   } Also, \code{C_c} is converted to a time series of length \code{nTimes}.
#' @examples
#' m <- createFlume(
#'   # TIME parameters
#'   timestep = 0.5,
#'   duration = 60, #hours
#'   # Water age parameters
#'   shape = "powerLaw",
#'   alpha = 1.4,
#'   tau_0 = 1/3600,
#'   tau_n = 60,
#'   # Initial concentrations and volumes in channel and hyporheic zone.
#'   C_c = 10,
#'   C_h = 5,
#'   V = 1,
#'   V_h = 0.5,
#'   # Integration parameter.  See help for pIntegrate Function
#'   # Set >1 if integration quality is problematic; slows the model.
#'   nSubdiv = 1
#' )
#'
#' #### RUN THE MODEL
#' simulateFlumeConcentration(m)
#'
#' # plot output
#' plot(m$times, m$C_c, xlab = "Time (hours)", ylab = "Concentration")
#'
#' # Assuming duration >= tau_n, compare simulated final conc to expected final conc
#' # to estimate of accuracy of model; smaller difference is better...
#' tail(m$C_c, 1) - m$C_final
#' @import hydrogeom
#' @importFrom purrr list_transpose
#' @importFrom stats approx integrate nlm
#' @export
simulateFlumeConcentration = function(m, debug = F) {

  if(m$tau_0 <= 0 & m$shape == "powerLaw") stop("tau_0 must be > 0.")
  if(m$tau_0 < 0) stop("tau_0 can't be < 0")

  # set up some values in the environment
  m$nIterations <- m$duration/m$timestep
  m$times <- seq(0, m$duration, length.out = m$nIterations+1)
  m$nTimes <- length(m$times)

  # Expected final concentration.
  m$V_c <- m$V - m$V_h
  m$C_final <- (m$C_c*m$V_c + m$C_h*m$V_h) / m$V
  m$C_c = rep(m$C_c, length(m$times))

  # get the functions associated with shape
  m$CCDF <- getFunction(paste0(m$shape, "CCDF"))
  m$IntCCDF <- getFunction(paste0(m$shape, "IntCCDF"))

  # copy "alpha" or "sigma" to shapeParam depending on "shape"
  CCDFFormals <- names(formals(m$CCDF))
  m$shapeParam <- get(CCDFFormals[length(CCDFFormals)], envir = m)

  #   index 1 when is t=0.  We already have all values for t=0, so we start
  # calculations with index 2.
  for(idx in 2:m$nTimes) {
    m$idx <- idx

    # calculate the integral for any t-tau that is less than zero. (this is
    # constant for the time step)
    if(m$times[idx] > m$tau_n) {
      preReleaseIntegral <- list(value = 0)
    } else {
      preReleaseBreaks <- logDistributedBreaks(min(m$times[idx], m$tau_n), m$tau_n, m$nSubdiv)
      preReleaseIntegral <- pIntegrate(C_hIntegrandStaticC, preReleaseBreaks, m = m, funName = "CCDF")
      if(debug) {
        m$debug = c(
          m$debug,
          c(list(type = "preRelease", ittr = idx), preReleaseIntegral)
        )
      }
    }

    # calculate integral for any t-tau that is between zero and last time step
    # (this is constant for the time step)
    if(idx == 2) {
      pastPostReleaseIntegral <- list(value = 0)
    } else {
      pastPostReleaseBreaks <- logDistributedBreaks(m$timestep, min(m$times[idx], m$tau_n), m$nSubdiv)
      pastPostReleaseIntegral <- pIntegrate(C_hIntegrandDynamicC, pastPostReleaseBreaks, t = m$times[m$idx], m = m, funName = "CCDF")
      if(debug) {
        m$debug <- c(
          m$debug,
          c(list(type = "pastPostRelease", ittr = idx), pastPostReleaseIntegral)
        )
      }
    }

    # calculate breaks for the "recent" interval (between last and current time
    # step, which needs to be optimized)
    recentPostReleaseBreaks <- logDistributedBreaks(m$tau_0, m$timestep, m$nSubdiv)

    # calculate the new concentration by descending on solution
    m$C_c[idx] <-
      nlm(
        optimizationError,
        m$C_c[idx-1],
        m = m,
        breaks = recentPostReleaseBreaks,
        staticIntegral = preReleaseIntegral$value,
        pastDynamicIntegral = pastPostReleaseIntegral$value,
        debug = debug
      )$estimate
  }

  # get rid of a few confusing things in the model environment.
  rm(list = c("CCDF", "IntCCDF", "shapeParam"), envir = m)

  TRUE
}

#' Integrands for determining the concentration of conservative tracer within
#' the hyporeheic zone
#'
#' If t represents simulation time, and t=0 is the time of the tracer release,
#' and tau is a water age of interest, these functions can be numerically
#' integrated to determine sum of past concentration at time t - tau, weighted
#' (multiplied) by the washout function (W(tau) -- the fraction of recharge
#' remaining in the hyporheic zone at water age tau).
#'
#' When the sum of the integrals of these functions are divided by the integral
#' of the washout function (W(tau)) from tau_0 to tau_n, the resulting value is
#' the mean concentration the hyporheic water.  See
#' \code{\link{optimizationError}} for details.
#'
#' @param tau A vector of water ages for which remaining-flow-weighted
#'   concentration is desired.
#' @param t Simulation time at which weighted concentration is desired
#' @param m An model environment.  See \code{\link{simulateFlumeConcentration}}
#'   for details.
#' @param funName Either "CCDF" or "PDF".  If "CCDF", the sum of the integrals
#'   of the C_hIntegrandDyanmic and C_hIntegrandStaticC will represent the mean
#'   concentration in the hyporheic zone.  If "PDF", the sum of the integrals
#'   will represent the upwelling concentration.
#' @return A vector of values representing past concentration weighted
#'   (multiplied) by the washout function:
#'
#'   \code{(W(tau) * C_c(t - tau))}.
#'
#'   \code{W(tau)} is the washout function (the CCDF of the age distribution of
#'   aquifer discharge). \code{W(tau)} assumes steady state hydrology and
#'   determines the fraction of recharge that is still in the aquifer at water
#'   age \code{tau}.
#'
#'   \code{C_hIntegrandStaticC} assumes \code{C_c} is a constant at the value of
#'   \code{C_c} prior to the release.  For any model iteration,
#'   \code{C_hIntegrandStaticC} should be integrated from the current model time
#'   to \code{tau_n}, but only when current model time is less than
#'   \code{tau_n}.
#'
#'   \code{C_hIntegrandDyanmic} looks backwards toward the time of release and
#'   estimates \code{C_c(t - tau)} using linear interpolation between known
#'   values of \code{C_c} at each time step.  For any model iteration,
#'   \code{C_hIntegrandDyanmic} should be integrated from \code{tau_0} to
#'   \code{min(current model time, tau_n)}.
#' @export
C_hIntegrandDynamicC = function(tau, t, m, funName){
  m[[funName]](tau, m$tau_0, m$tau_n, m$shapeParam) * approx(m$times[1:length(m$C_c)], m$C_c, t-tau)$y #* e^(m$k_h*tau)
}

#' @rdname  C_hIntegrandDynamicC
#' @export
C_hIntegrandStaticC = function(tau, m, funName) {
  # m$C_h[1] is the concentration prior to release
  m[[funName]](tau, m$tau_0, m$tau_n, m$shapeParam) * m$C_h[1] #* e^(m$k_h*tau)
}


#' Determine error associated with an estimate of channel concentration.
#'
#' Objective function for an implicit solution to determine the channel
#' concentration.  The measure of error returned is the difference between: 1)
#' the mean concentration of the channel and hyporheic zone, weighted by water
#' storage in each, at time t; and 2) the expected final concentration of water
#' in the flume.
#'
#' \code{optimazationError} sums \code{staticIntegral} with integration of
#' \code{\link{C_hIntegrandDynamicC}} from tau_0 to min(currentModelTime,
#' tau_n).  This way, the numerical integration doesn't span the discontinuity
#' in concentration associated with the slug addtion.  The two summed values are
#' then divided by the integral of W(tau) from tau_0 to tau_n, the washout
#' function to convert the sum to the mean conservative tracer concentraion in
#' the hyporehic zone.
#' @param C_c An estimate of channel concentration in the next model time step.
#'   Values are generated by nlm() to minimim the return value of
#'   \code{optimizationError}.
#' @param m The model environment, which must contain variables \code{tau_0,
#'   tau_n, alpha} or \code{sigma} (depending on "shape), \code{idx, times,} and \code{C_c}.  \code{idx} is used as a vector
#'   index, and points at values associated with the current model time.
#'   \code{idx} = (current model iteration + 1), because \code{t=0} is
#'   associated with \code{idx=1}.  \code{times} is a vector of simulation times
#'   corresponding to time steps and \code{C_c} is a vector of channel
#'   concentrations.
#' @param breaks See \code{\link{pIntegrate}}
#' @param staticIntegral The integral of \code{\link{C_hIntegrandStaticC}}
#'   (pre-release concentration weighted by W(tau), the washout function) from
#'   the current simulation time to tau_n, or zero if current simulation time is
#'   greater than tau_n.
#' @export
optimizationError <- function(C_c, m, breaks, staticIntegral, pastDynamicIntegral, debug = F) {
  # set the current channel concentration equal to the estimate.
  m$C_c[m$idx] <- C_c
  # determine what the current hyporheic concentration would be given the
  # estimate of C_c
  recentPostReleaseIntegral <- pIntegrate(C_hIntegrandDynamicC, breaks, t = m$times[m$idx], m = m, funName = "CCDF")
  C_h = (staticIntegral + pastDynamicIntegral + recentPostReleaseIntegral$value) /
    m$IntCCDF(m$tau_0, m$tau_n, m$tau_0, m$tau_n, m$shapeParam)
  # determine the error -- post release, the weighted mean of hyporheic conc and
  # channel conc should at all times be equal to final conc.
  result <- ((C_h*m$V_h + C_c*m$V_c)/m$V - m$C_final)^2
  if(debug) m$debug <- c(m$debug, c(list(type = "Optimize", ittr = m$idx), recentPostReleaseIntegral))
  result
}

# Generate a ln()-distributed sequence from tau_0 to tau_n
logDistributedBreaks = function(tau_0, tau_n, nBins){
  exp(seq(log(tau_0 + 1), log(tau_n + 1), length.out = nBins+1)) - 1
}

#' Numerical integration with manual subdivisons.
#'
#' A function that integrates by summing the numerical integration across
#' multiple consecutive ranges predetermined by the user.
#'
#' If the shape of an integrand is known ahead of time, choosing specific breaks
#' can ease or speed the integration.
#' @param funs A single function or a list of functions that will be used for
#'   each interval between \code{breaks}.  A list of multiple functions will be
#'   recycles for the intervals, with a warning if number of intervals between
#'   breaks is not an even multiple of length of the list of functions.  Summing
#'   the integral of multiple functions can be useful if, for instance, there is
#'   a discontinuity in the integrand, with one function describing the integrand
#'   to the left of the discontinuity, and the other describing the integrand to
#'   the right of the discontinuity.  In such a case, one of the breaks should
#'   be the x value of the discontinuity!
#' @param breaks A vector of values of representing the divisions among
#'   consecutive intervals to be integrated. The length of \code{breaks} is
#'   therefore the number of desired intervals + 1.  For instance, if
#'   \code{breaks} = c(1,2,5,10), \code{pIntegrate} will integration \code{funs}
#'   from 1 to 2, from 2 to 5, and from 5 to 10, and then sum the results. If
#'   there are only two values in the "breaks" list, this function is identical
#'   to integrate.
#' @param ... pIntegrate calls \code{\link{integrate}} for each subdivision.
#'   Values in '...' are passed to \code{\link{integrate}} unchanged, for each
#'   subdivision.
#' @export
pIntegrate = function(funs, breaks, ...) {
  if(is.function(funs)) funs = list(funs)
  if(length(breaks) == 2) {
    if(length(funs)>1) stop("number of functions is greater than number of intervals.")
    return(integrate(funs[[1]], breaks[1], breaks[2], ...))
    ## end of function if length(breaks == 2)
  }
  pieces =
    purrr::list_transpose(
      mapply(
        integrate,
        lower = breaks[-length(breaks)],
        upper = breaks[-1],
        f = funs,
        MoreArgs = list(
          ...
        ),
        SIMPLIFY = F
      )
    )[-5]
  structure(
    list(
      value = sum(unlist(pieces[["value"]])),
      abs.error = sum(unlist(pieces[["abs.error"]])),
      subdivisions = unlist(pieces[["subdivisions"]]),
      message = paste(unique(unlist(pieces[["message"]])), collapse = ", "),
      call = "pIntegrate" #deparse(match.call())
    ),
    class = "integrate"
  )
}

#' Create a flume model environment
#'
#' Creates a flume model environment for use with
#' \code{\link{simulateFlumeConcentration}}
#'
#' Pass named parameters to \code{...}.  Each parameter will be created as a
#' variable of the same name in an environment.  Pass one parameter for each
#' variable required by \code{\link{simulateFlumeConcentration}}. See
#' \code{\link{simulateFlumeConcentration}} documentation for description of
#' required parameters and see example under
#' \code{\link{simulateFlumeConcentration}} for example of using
#' \code{createFlume}
#'
#' @param ... Named parameter that will be used to create variable in an
#'   environment returned to the user.
#' @return An environment containing variables with names and values determined
#'   by the parameters passed to \code{...}.
#' @export
createFlume <- function(...) {
  m <- new.env()
  params = list(...)
  if(any(names(params) == "")) stop("All parameters must be named.")
  mapply(
    assign,
    x = names(params),
    value = params,
    MoreArgs = list(envir = m)
  )

  m
}

#' Post processing variables
#'
#' Add addtional hydrology descriptors by post-processing the model result.
#'
#' @param m Model environment that has been executed with
#'   \code{\link{simulateFlumeConcentration}}
#' @return Returns TRUE silently.  As a side effect, creates a \code{C_up} time series
#' in the \code{m} environment representing the concentration of upwelling water.
#' @export
postProcess <- function(m) {
  if (length(m$C_c) == 1)
    stop("You must execute the model before post processing.")
  m$PDF <- getFunction(paste0(m$shape, "PDF"))
  PDFFormals <- names(formals(m$PDF))
  m$shapeParam <- get(PDFFormals[length(PDFFormals)], envir = m)
  m$C_up <- sapply(m$times, function(t) {
    # print(t)
    # if(t >= m$tau_n-0.1) {
    #   print("yay")
    # }

    # When t >= tau_n there is no longer a pre-release influence
    if (t >= m$tau_n) {
      preReleaseIntegral <- list(value = 0)
    # Otherwise we calculate the pre-release influence.
    } else {
      # if t <= tau_0, then no solute has reached the HZ yet.  Therefore, the
      # upwelling concentration must be the pre-release concentration, which
      # will be returned C_hIntegrandStacticC when integrated from tau_0 to
      # tau_n
      if(t <= m$tau_0) {
        preReleaseBreaks <- logDistributedBreaks(min(m$tau_0, m$tau_n),
                                                 m$tau_n, m$nSubdiv)
      # any other time, we integrate from t to tau because tau_n affects only
      # the post release integral once t > tau_0
      } else {
        preReleaseBreaks <- logDistributedBreaks(min(t, m$tau_n),
                                                 m$tau_n, m$nSubdiv)
      }
      # preReleaseBreaks now contains the limits of integration, so just
      # integrate C_hIntegrandStaticC!
      preReleaseIntegral <- pIntegrate(C_hIntegrandStaticC,
                                       preReleaseBreaks, m = m, funName = "PDF")
    }

    # for post release, if t <= tau_0 there is no post-release influence
    # because, again, the salt hasn't reached the HZ
    if (t <= m$tau_0) {
      postReleaseIntegral <- list(value = 0)
    # otherwise, integrate C_hIntegrandDynamicC from tan_0 to min(t, taU_n)
    } else {
      postReleaseBreaks <- logDistributedBreaks(m$tau_0, min(t,
                                                             m$tau_n), m$nSubdiv)
      postReleaseIntegral <- pIntegrate(C_hIntegrandDynamicC,
                                        postReleaseBreaks, t = t, m = m, funName = "PDF")
    }
    preReleaseIntegral$value + postReleaseIntegral$value
  })
  rm(list = c("PDF", "shapeParam"), envir = m)
  invisible(TRUE)
}
