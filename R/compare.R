##' Default compare function (log likelihood)
##'
##' See mcstate docs
##'
##' @title Comparison function with data
##' @param state state from model
##' @param observed observed data
##' @param pars potential additional parameters
##' @return real
##' @author Pete Dodd
##' @export
case_compare7 <- function(state, observed, pars = NULL) {
  ans <- rep(0, dim(state)[2]) # no particles long
  for (i in 1:7) {
    ## n particles long @ given timestep
    ## ln7/bn7 index raw notes[] and N[] state summed over all age x HIV
    ## strata for zone i -- true per-100k rate is 1e5*sum(notes)/sum(N),
    ## NOT sum(notes*N) (which double-counts population and isn't a rate)
    totnotes <- colSums(state[BLASTtbmod::ln7[[i]], , drop = TRUE])
    totpops <- colSums(state[BLASTtbmod::bn7[[i]], , drop = TRUE])
    notes_modelled <- 1e5 * totnotes / totpops # per 100,000
    notes_observed <- observed[[paste0("notifrate_", i)]]
    ans <- ans + dnorm(x = notes_modelled,
                       mean = notes_observed,
                       sd = 30,    #TODO - allow to vary?
                       log = TRUE) # sums densities across each patch
  }
  ans
}
## original sd was 50
