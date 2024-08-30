##' This creates an object for running model as particle filtering etc.
##'
##' See mcstate docs
##' @title Create a particle filter
##' @param data must be mcstate particle filter data
##' @param comparefun mcstate compare function (log likelihood)
##' @param n_particles 
##' @param n_threads 
##' @param seed 
##' @param ... 
##' @return mcstate particle filter
##' @author Pete Dodd
##' @import mcstate
##' @export
create.particlefilter <- function(data,
                                  comparefun,
                                  n_particles,
                                  n_threads = 4,
                                  seed = 5L,
                                  ...){
  mcstate::particle_filter$new(
    data = data,
    model = BLASTtbmod:::stocm,
    compare = comparefun,
    n_particles = n_particles,
    n_threads = 4, # NOTE
    seed = 5L,
    ...
  )
}
