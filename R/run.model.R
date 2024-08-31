##' This will run the stochastic version of the model without reference to data
##'
##' TODO
##' @title Run dust version of model
##' @param parm 
##' @param times 
##' @param n.particles 
##' @param convert 
##' @return array
##' @author Pete Dodd
##' @import dust
##' @export
run.model <- function(parm,
                      times,
                      n.particles = 1,
                      convert = FALSE) {
  ## create dust model
  dmod <- BLASTtbmod:::stocm$new(
    pars = parm,
    time = min(times) + 1, # TODO
    n_particles = n.particles
  )
  ## model loop over time
  y_looped <- array(NA, dim = c(dmod$info()$len, n.particles, max(times)))
  for (t in 1:max(times)) {
    y_looped[, , t] <- dmod$run(t)
  }
  if (convert == TRUE) { # convert to df
    y_looped <- as.data.frame(t(y_looped[, 1, ])) # TODO cope with multiple parms?
    ## Make it readable
    if (ncol(y_looped) == length(BLASTtbmod::get_cols)) {
      colnames(y_looped) <- BLASTtbmod::get_cols
    }
  }
  y_looped
}
