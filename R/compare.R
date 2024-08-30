##' Default compare function (log likelihood)
##'
##' See mcstate docs
##'
##' @title Comparison function with data
##' @param state 
##' @param observed 
##' @param pars 
##' @return real
##' @author Pete Dodd
##' @export
case_compare7 <- function(state, observed, pars = NULL) {
  ans <- rep(0, dim(state)[2]) # no particles long
  for (i in 1:7) {
    ## n particles long @ given timestep
    ## NOTE changed to calculate total note before aggregation and then calculate rates
    ## rate x pop = 7 total notes x 100k
    totnotes <- colSums(
      state[BLASTtbmod::n7[[i]], , drop = TRUE] *
      state[BLASTtbmod::N7[[i]], , drop = TRUE]
    )
    totpops <- colSums(state[BLASTtbmod::N7[[i]], , drop = TRUE])
    notes_modelled <- totnotes / totpops # back to rate
    notes_observed <- observed[[paste0("notifrate_", i)]]
    ans <- ans + dnorm(x = notes_modelled,
                       mean = notes_observed,
                       sd = 30,    #TODO - allow to vary?
                       log = TRUE) # sums densities across each patch
  }
  ans
}
## original sd was 50
