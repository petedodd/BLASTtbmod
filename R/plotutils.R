## ## some utilities for plotting

##' Collect total populations as data table to be used by each plotter
##'
##' TODO
##' @title Postprocess model output
##' @param Y 
##' @param chain_step_num 
##' @param out_type 
##' @return data.table
##' @author Pete Dodd
##' @import data.table
##' @export
extract.pops <- function(Y, chain_step_num = 1, out_type = "N") {
  rownames(Y) <- BLASTtbmod::get_cols
  selnmz <- grep(pattern = paste0(out_type, "\\["), rownames(Y), value = TRUE)
  D <- data.table::as.data.table(Y[row.names(Y) %in% selnmz,
    chain_step_num, ,
    drop = TRUE
  ])
  D <- cbind(selnmz, D)
  D <- data.table::transpose(D, make.names = "selnmz")
  colnames(D) <- selnmz
  D[, step := 1:nrow(D)]
  D <- melt(D, id = "step")
  D[, c("patch", "acatno", "hivno") := data.table::tstrsplit(variable, split = ",")]
  D[, patch := as.integer(gsub(paste0(out_type, "\\["), "", patch))]
  D[, acatno := as.integer(acatno)]
  D[, hivno := as.integer(gsub("\\]", "", hivno))]
  D[, age := BLASTtbmod::agz[acatno]]
  D[, hiv := BLASTtbmod::hivz[hivno]]
  D[, patch := paste("Patch", patch)]
  D$age <- factor(D$age, levels = BLASTtbmod::agz)
  D$hiv <- factor(D$hiv, levels = BLASTtbmod::hivz)
  D[, c("acatno", "hivno") := NULL]
  D
}



##' Collect total populations as data table to be used by each plotter
##'
##' TODO check chains or particles?
##' @title Postprocess model outputs for multiple chains
##' @param Y 
##' @param n_chain_steps 
##' @param out_type 
##' @return data.table
##' @author Pete Dodd
##' @import data.table
##' @export
extract.pops.multi <- function( Y, n_chain_steps = 100, out_type = 'N' ){
  ## catch
  if( dim(Y)[2] < n_chain_steps ){
    n_chain_steps <- dim(Y)[2]
    print(paste0(
      "Warning: Setting N chain steps to max (",
      dim(Y)[2], ")"
    ))
  }
  samp <- sample(dim(Y)[2], n_chain_steps, replace = FALSE)
  ## Rearrange & sample desired output 
  rownames(Y) <- BLASTtbmod::get_cols
  D <- data.table::as.data.table( Y[, samp, ] )
  selnmz <- grep( pattern = paste0( out_type, '\\[' ), D$V1 )
  D <- D[ selnmz, ]
  colnames(D) <- c('subcomp','chain_step','t',out_type )
  ## columns for each index
  D[ ,c('patch','acatno','hivno'):=data.table::tstrsplit(subcomp,split=',')]
  D[ ,patch:=as.integer(gsub(paste0( out_type, '\\[' ),'',patch))]
  D[ ,acatno:=as.integer(acatno)]
  D[ ,hivno:=as.integer(gsub('\\]','',hivno))]
  D[ ,age:=BLASTtbmod::agz[acatno]]
  D[ ,hiv:=BLASTtbmod::hivz[hivno]]
  D[ ,patch:=paste( 'Patch', patch)]
  D$age <- factor(D$age,levels=BLASTtbmod::agz)
  D$hiv <- factor(D$hiv,levels=BLASTtbmod::hivz)
  D[, c('acatno', 'hivno'):=NULL]
  D
}



##' Default median N per age group per year (not by patch)
##'
##' TODO
##' @title Plot comparing total populations
##' @param Y 
##' @param n_chain_steps 
##' @param start_year 
##' @param years 
##' @param eps 
##' @param by_comp 
##' @return ggplot2
##' @author Pete Dodd
##' @import data.table
##' @import ggplot2
##' @export
plot_compare_demog <- function(Y,
                               n_chain_steps = 100,
                               start_year=2015,years=6,
                               eps=0.25, #quantile to use
                               by_comp="none"
                               ){
  if (!by_comp %in% c("age", "none")) stop("by_comp must be 'none' or 'age'!")
  dt <- 1 / 12 # months
  ## read MWI data
  real_dat <- BLASTtbmod::MWI$N[Year %in% start_year:(start_year + years)]
  real_dat[["age"]] <- real_dat$AgeGrp
  real_dat[["mid"]] <- real_dat$PopTotal * 1000
  real_dat[["lo"]] <- NA_real_
  real_dat[["hi"]] <- NA_real_
  tot_tmp_rl <- real_dat[Year == start_year, sum(mid)] # total population start year
  real_dat$mid_scaled <- real_dat$mid / tot_tmp_rl
  ## reshape simulation data
  D <- extract.pops.multi( Y, n_chain_steps )
  D[,'subcomp' := NULL ]
  ## sum over HIV
  D <- D[, list(tot_N = sum(N)), by = list(chain_step, t, age)]
  ## Scale to total pop at start by chainstep/particle at time step 1:
  tot_tmp <- D[t == 1, list(pop_start = sum(tot_N)), by = chain_step]
  D <- merge(D, tot_tmp, by = "chain_step")
  D$pop_scaled <- D$tot_N / D$pop_start
  ## aggregate if necessary & summarize
  if (by_comp == "age") {
        aggD <- D[, list(
          mid_scaled = median(pop_scaled),
          lo = quantile(pop_scaled, eps),
          hi = quantile(pop_scaled, 1 - eps)
        ),
        by = list(t, age)
        ]
  } else { #all ages
    real_dat <- real_dat[, list(mid_scaled = sum(mid_scaled)), by = Year]
    real_dat[, c("lo", "hi", "age") := list(NA_real_, NA_real_, "all")]
    ## NOTE chained:
    aggD <- D[, list(pop_scaled = sum(pop_scaled)), by = list(chain_step, t)
              ][, list(
      mid_scaled = median(pop_scaled),
      lo = quantile(pop_scaled, eps),
      hi = quantile(pop_scaled, 1 - eps)
    ),
    by = t
    ]
  }
  aggD[, "Year" := (t - 1) * dt + start_year]
  ## construct plot
  p <- ggplot2::ggplot( aggD, aes( Year, y = mid_scaled, ymin = lo, ymax = hi)) +
    ggplot2::geom_ribbon(alpha = 0.3) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~age) +
    ggplot2::geom_point(data = real_dat,col=2,shape = 16) +
    ggplot2::ylab("Population by year (scaled at t=1)") +
    ggplot2::xlab("Time") +
    ggplot2::ggtitle("Population, median + 50% PI vs real data (red)")
  p + ggplot2::expand_limits(y = c(0, NA))
}


## library(scales)
## library(plyr)
## library(dplyr)
## library(ggh4x)




## # Separate function to deal with aggregating rates
## # Defaults to providing notifrates over patches only
## # Currently will only do patch in line with current real data
## # Looks like real data has hiv/art/age available too, so can do this later?
## # Data md7 must be in env
## plot.compare.noterate.agrgt <- function( Y, n_chain_steps = 100,
##                                          agrgt_by = c( 'patch' )){
  
##   D1 <- extract.pops.multi( Y, n_chain_steps, out_type='N' )
##   D1[,'subcomp' := NULL ]
##   D2 <- extract.pops.multi( Y, n_chain_steps, out_type='notes' )
##   D2[,'subcomp' := NULL]
##   D <- merge( D1, D2, by = c( 'chain_step', 't', 'patch', 'age', 'hiv'))

##   agrgt_by <- c( agrgt_by, 'chain_step', 't' )
##   aggD <- D %>%
##     group_by(across(all_of(agrgt_by))) %>%
##     summarise( tot_N = sum(N),
##                tot_notes = sum(notes)) %>%
##     as.data.table()
##   aggD[, noterate := tot_notes/tot_N*1e5 ]
  
##   # Generate ribbon
##   eps = 0.25
##   aggD <- aggD[ ,.( mid=median( noterate ),
##                lo = quantile( noterate, eps ),
##                hi = quantile( noterate, 1-eps )),
##             by=.(t, patch )]
  
##   real_dat = md7
##   real_dat[[ 'patch' ]] <- paste( 'Patch', md7$comid )
  
##   p <- ggplot( aggD, aes( t, y = mid, ymin = lo, ymax = hi)) +
##     geom_ribbon(alpha = 0.3) +
##     geom_line() +
##     facet_wrap(~patch) +
##     geom_point(data = real_dat, col = 2, shape = 1) +
##     ylab("Notification rate per month") +
##     xlab("Month") +
##     ggtitle("Real data, median + 50% PI")
##   plot.new()
##   print(p)

## }





## # general plotter - takes state from mcstate
## # Either pmcmc_run$trajectories$state object or forecast$state object or filter.stocmmodr$history()
## # will plot any output (out_type) by patch, age and/or hiv 
## # (specified as character vector by_comp). Choose any chain_step.
## # Defaults to chain_step 1, total population, aggregates over compartments
## plot.pops.dynamic <- function( Y, chain_step_num=1, out_type='N',
##                        by_comp=NULL ){

##   D <- extract.pops( Y, chain_step_num, out_type )
##   grp_vars <- c( 'step', by_comp )
##   Dplot <- D %>% 
##     group_by( across( all_of( grp_vars ))) %>%
##     summarise( tot_pop = sum(value))
##   Dplot$default_col <- 1
  
##   if( 'age' %in% by_comp ){
##     by_comp <- by_comp[-which(by_comp=='age')]
##     col_age <- T
##   } else {
##     col_age <- F
##   }
##   out_lab <- out_type # default
##   if( out_type == 'N' ){
##     out_lab <- 'Total population'
##   } else {
##     if( out_type == 'incidence'){
##       out_lab <- 'Incidence'
##     } else {
##       if( out_type == 'notifrate'){
##         out_lab <- 'Notification rate'
##       }
##     }
##   }
  
##   p <- ggplot( Dplot, aes( x=step*dt+start_year, y = tot_pop )) +
##     xlab( 'Step' )+
##     ylab( out_lab )+
##     scale_y_continuous(label=comma)+
##     guides(color = guide_legend(title = "Age group"))+ 
##     theme_light()
##   if( col_age ){
##     p <- p + geom_line( aes( colour = age ))
##   }else{
##     p <- p + geom_line()
##   }
##   if( length( by_comp > 0)){
##     p <- p + facet_grid2( by_comp, 
##                           independent = TRUE,
##                           scales='free')
##   }
##   print(p)
## }


## ## demographic snapshot (bar)
## demo.snapshot <- function( Y, chain_step_num=1, out_type='N', timestep=NULL ){
  
##   D <- extract.pops( Y, chain_step_num, out_type )
##   if( length(timestep > 0)){
##     p <- ggplot( D[ step==timestep, ], aes(age,value))
##   } else {
##     p <- ggplot( D, aes( age, value))
##   }
##   out_lab <- out_type # default
##   if( out_type == 'N' ){
##     out_lab <- 'Total population'
##   } else {
##     if( out_type == 'incidence'){
##       out_lab <- 'Incidence'
##     } else {
##       if( out_type == 'notifrate'){
##         out_lab <- 'Notification rate'
##       }
##     }
##   }
  
##   p <- p + geom_bar(stat='identity')+
##     facet_wrap(~patch)+
##     scale_y_continuous(label=comma)+
##     theme_light()+
##     theme(axis.text.x = element_text(angle = 90, hjust = 1))+
##     xlab('Age')+
##     ylab(out_lab)
##   print(p)
## }


##' HIV/ART over time
##'
##' TODO
##' @title HIV by time
##' @param Y 
##' @param start_year 
##' @param chain_step_num 
##' @param out_type 
##' @param by_age 
##' @return ggplot
##' @author Pete Dodd
##' @import ggplot2
##' @import scales
##' @import data.table
##' @export
plot_HIV_dynamic <- function( Y, start_year, chain_step_num=1, out_type='N', by_age = FALSE){
  ## TODO total population
  D <- extract.pops( Y, chain_step_num, out_type )
  dt <- 1/12 #monthly
  ## reformat
  if( by_age == FALSE ){
    X <- data.table::dcast( D,patch+step ~ hiv, fun.aggregate = sum )
    vars <- c("HIVpc", "ARTpc", "patch", "step")
    id <- c("patch", "step")
  } else {
    if( by_age == TRUE ){
      X <- data.table::dcast(D, patch + step + age ~ hiv, fun.aggregate = sum)
      id <- c("patch", "step", "age")
      vars <- c("HIVpc", "ARTpc", "patch", "step", "age")
    } else { stop('by_age must be TRUE or FALSE') }
  }
  X[,tot:= `HIV-`+`HIV+/ART-`+`HIV+/ART+`]
  X[,HIVpc:=(`HIV+/ART-`+`HIV+/ART+`)/tot]
  X[,ARTpc:=`HIV+/ART+`/(`HIV+/ART-`+`HIV+/ART+` + 1e-10)]
  X <- data.table::melt(X[,..vars],id=id)
  out_lab <- out_type # default
  if( out_type == 'N' ){
    out_lab <- "Population"
  } else {
    if( out_type == 'incidence'){
      out_lab <- "Incidence"
    } else {
      if(out_type == "notifrate" ){
        out_lab <- "Notification rate"
      }
    }
  }
  ## plot
  wrap_vars <- id[-which(id=='step')]
  ggplot2::ggplot(X, aes(step*dt+start_year, value, col=variable))+
    ggplot2::geom_line()+
    ggplot2::facet_grid(wrap_vars)+
    ggplot2::scale_y_continuous(label=scales::percent)+
    ggplot2::theme_linedraw()+
    ggplot2::xlab('Year')+
    ggplot2::ylab( paste( out_lab, '(Proportion)'))+
    ggplot2::theme(legend.position = 'top',
          legend.title=element_blank())+
    ggplot2::scale_color_discrete(labels = c("HIV+/ART-", "HIV+/ART+"))
}



##' HIV snapshot by age at single time
##'
##' TODO
##' @title HIV snapshot
##' @param Y 
##' @param chain_step_num 
##' @param out_type 
##' @param timestep 
##' @return ggplot
##' @author Pete Dodd
##' @export
##' @import ggplot2
##' @import data.table
##' @import scales
plot_HIV_snapshot <- function( Y, chain_step_num=1, out_type='N', timestep=NULL ){
  ## TODO improve what this actually plots
  D <- extract.pops( Y, chain_step_num, out_type )
  ## reformat
  if( length( timestep > 0)){
    X <- data.table::dcast( D[ step==timestep, ], patch+age ~ hiv, fun.aggregate = sum)
  } else {
    X <- data.table::dcast(D, patch + age ~ hiv, fun.aggregate = sum)
  }
  X <- data.table::melt(X, id = c("patch", "age"))
  ## plot
  ggplot2::ggplot(X, aes( age, value, fill=variable))+
    ggplot2::geom_bar(stat='identity',position = position_fill())+
    ggplot2::facet_wrap(~patch)+
    ggplot2::scale_y_continuous(label=scales::percent)+
    ggplot2::theme_linedraw()+
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = 'top',legend.title=element_blank())+
    ggplot2::xlab("Age") + ggplot2::ylab("Population")
}


## ## TB incidence by time
## TB.dynamics <- function( Y, chain_step_num=1, out_type='incidence', separate=FALSE, by_age=F){
  
##   D <- extract.pops( Y, chain_step_num, out_type )
##   by_comp <- 'patch'

##   if( separate==FALSE ){
##     D[,hiv:=ifelse(hiv=='HIV-','HIV-','HIV+')]
##     D <- D[,.(value=sum(value)),by=.(step,patch,age,hiv)]
##   }
##   if( by_age==FALSE ){
##     D <- D[,.(value=sum(value)),by=.(step,patch,hiv)]
##   } else {
##     by_comp <- c( by_comp, 'age' )
##   }
  
##   ggplot(D, aes(step*dt+start_year, value, col=hiv))+
##     geom_line()+
##     facet_grid2( by_comp, 
##                  independent = TRUE,
##                  scales='free')+
##     theme_light()+
##     xlab('Year')+
##     ylab('TB incidence')
## }


# # To take either pmcmc_run$trajectories$state object or forecast$state object or filter.stocmmodr$history()
# # requires TBN to be in env
# #  To complete if patch-specific real N available
# plot.compare.dynamic <- function( Y, n_chain_steps = 100, out_type = 'N' ){
#   
#   D <- extract.pops.multi( Y, n_chain_steps, out_type )
#   D[ ,'subcomp':= NULL ]
#   
#   if( out_type == 'N'){
#     real_dat = TBN
#     real_dat[[ 't' ]] <- real_dat$month
#     real_dat[[ 'patch' ]] <- real_dat$zone
#   }
#   
#   agrgt_by <- c( 'chain_step', 't', 'patch' )
#   aggD <- D %>% 
#     group_by(across(all_of(agrgt_by))) %>%
#     summarise( tot_N = sum(N)) %>%
#     as.data.table()
#   
# }
