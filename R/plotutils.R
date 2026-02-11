## ## some utilities for plotting

##' Collect total populations as data table to be used by each plotter
##'
##' TODO
##' @title Postprocess model output
##' @param Y input data
##' @param chain_step_num chain step
##' @param out_type output type
##' @return data.table
##' @author Pete Dodd
##' @import data.table
##' @export
extract.pops <- function(Y, chain_step_num = 1, out_type = "N") {
  rownames(Y) <- BLASTtbmod::get_cols
  out_type <- paste0(out_type, "\\[")
  out_type <- paste0(out_type, collapse = "|") # handling vectors of out_type
  selnmz <- grep(pattern = out_type, rownames(Y), value = TRUE)
  D <- data.table::as.data.table(Y[row.names(Y) %in% selnmz,
    chain_step_num, ,
    drop = TRUE
    ])
  if (length(chain_step_num) == 1) {
    D <- cbind(selnmz, D)
    D <- data.table::transpose(D, make.names = "selnmz")
    colnames(D) <- selnmz
    D[, step := 1:nrow(D)]
    D <- data.table::melt(D, id = "step")
  } else {
    names(D)[1:3] <- c("variable","particle","step")
  }
  D[!grepl(",",variable),variable:=gsub("\\]",",,\\]",variable)]#to allow handing of patch totals
  D[, c("patch", "acatno", "hivno") := data.table::tstrsplit(variable, split = ",")]
  D[, patch := as.integer(gsub(out_type, "", patch))]
  D[, acatno := as.integer(acatno)]
  D[, hivno := as.integer(gsub("\\]", "", hivno))]
  D[, age := BLASTtbmod::agz[acatno]]
  D[, hiv := BLASTtbmod::hivz[hivno]]
  D[, patch := paste("Patch", patch)]
  D$age <- factor(D$age, levels = BLASTtbmod::agz)
  D$hiv <- factor(D$hiv, levels = BLASTtbmod::hivz)
  D[, varname := gsub("[[:punct:]]+|[[:digit:]]+", "", variable)]
  D[, c("acatno", "hivno") := NULL]
  D
}



##' Collect total populations as data table to be used by each plotter
##'
##' TODO check chains or particles?
##' @title Postprocess model outputs for multiple chains
##' @param Y input data
##' @param n_chain_steps chain step
##' @param out_type output type
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
##' @param Y input data
##' @param n_chain_steps chain step
##' @param start_year start year
##' @param years number of years
##' @param eps UI width
##' @param by_comp plot type
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


# Separate function to deal with aggregating rates
# Defaults to providing notifrates over patches only
# Currently will only do patch in line with current real data
# Looks like real data has hiv/art/age available too, so can do this later?
# Data md7 must be in env
#'
#' @title Plot average notification rates & compare with real data
#' @param Y input data
#' @param n_chain_steps chain steps
#' @param agrgt_by aggregate by
#' @param realdata logical (default TRUE) on whether to show real data points
#' @param start_year start for scale (default 2015)
#' @return ggplot2
#' @import data.table
#' @import ggplot2
#' @importFrom dplyr group_by summarise across all_of
#' @importFrom magrittr %>%
#' @export
plot_compare_noterate_agrgt <- function(Y,
                                        n_chain_steps = 100,
                                        agrgt_by = c("patch"),
                                        realdata = TRUE,
                                        start_year = 2015) {
  D1 <- extract.pops.multi(Y, n_chain_steps, out_type = "N")
  D1[, "subcomp" := NULL]
  D2 <- extract.pops.multi(Y, n_chain_steps, out_type = "notes")
  D2[, "subcomp" := NULL]
  D <- merge(D1, D2, by = c("chain_step", "t", "patch", "age", "hiv"))

  agrgt_by <- c(agrgt_by, "chain_step", "t")
  aggD <- D %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(agrgt_by))) %>%
    dplyr::summarise(
      tot_N = sum(N),
      tot_notes = sum(notes)
    ) %>%
    data.table::as.data.table()
  aggD[, noterate := tot_notes / tot_N * 1e5]

  # Generate ribbon
  eps = 0.25
  aggD <- aggD[, .(
    mid = median(noterate),
    lo = quantile(noterate, eps),
    hi = quantile(noterate, 1 - eps)
  ),
  by = .(t, patch)
  ]
  aggD <- aggD[t > 1]
  real_dat = BLASTtbmod::md7
  real_dat[["patch"]] <- paste("Patch", md7$comid)

  p <- ggplot2::ggplot(aggD, aes(t / 12 + start_year, y = mid, ymin = lo, ymax = hi)) +
    ggplot2::geom_ribbon(alpha = 0.3) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~patch) +
    ggplot2::ylab("Notification rate per month") +
    ggplot2::xlab("Month") +
    ggplot2::ggtitle("Real data, median + 50% PI") +
    ggplot2::theme_light() +
    ggplot2::expand_limits(y = c(0, NA))
  if (realdata) {
    p <- p + ggplot2::geom_point(data = real_dat, col = 2, shape = 1)
  }
  print(p)
}


# general plotter - takes state from mcstate
# Either pmcmc_run$trajectories$state object or forecast$state object or filter.stocmmodr$history()
# will plot any output (out_type) by patch, age and/or hiv
# (specified as character vector by_comp). Choose any chain_step.
# Defaults to chain_step 1, total population, aggregates over compartments
#' @title Plot selected populations/compartments over time (absolute numbers, no comparison w/real data)
#'
#' @param Y input data
#' @param chain_step_num chain step
#' @param out_type output type
#' @param by_comp comparison type
#' @importFrom dplyr group_by summarise across all_of
#' @import ggplot2
#' @importFrom ggh4x facet_grid2
#' @importFrom magrittr %>%
#' @return ggplot2
#' @export
plot_pops_dynamic <- function(Y,
                              chain_step_num = 1,
                              out_type = "N",
                              start_year = 2015,
                              by_comp = NULL) {
  D <- extract.pops(Y, chain_step_num, out_type)
  grp_vars <- c("step", by_comp)
  Dplot <- D %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grp_vars))) %>%
    dplyr::summarise(tot_pop = sum(value))
  Dplot$default_col <- 1
  dt <- 1 / 12

  if ("age" %in% by_comp) {
    by_comp <- by_comp[-which(by_comp == "age")]
    col_age <- T
  } else {
    col_age <- F
  }
  out_lab <- out_type # default
  if (out_type == "N") {
    out_lab <- "Total population"
  } else {
    if (out_type == "incidence") {
      out_lab <- "Incidence"
    } else {
      if (out_type == "notifrate") {
        out_lab <- "Notification rate"
      }
    }
  }
  p <- ggplot2::ggplot(Dplot, aes(x = step * dt + start_year, y = tot_pop)) +
    ggplot2::xlab("Step") +
    ggplot2::ylab(out_lab) +
    ggplot2::scale_y_continuous(label = comma) +
    ggplot2::guides(color = guide_legend(title = "Age group")) +
    ggplot2::theme_light()
  if (col_age) {
    p <- p + ggplot2::geom_line(aes(colour = age))
  } else {
    p <- p + ggplot2::geom_line()
  }
  if (length(by_comp > 0)) {
    p <- p + ggh4x::facet_grid2(by_comp,
      independent = TRUE,
      scales = "free"
    )
  }
  print(p)
}


##' HIV/ART over time
##'
##' Plotting of HIV/ART over time in various ways.
##' If `by_age = FALSE` then age 15-49 years subpopulation will be used.
##'
##' @title HIV by time
##' @param Y input data
##' @param start_year start year
##' @param chain_step_num chain step
##' @param out_type output type
##' @param by_age by age or not
##' @param by_patch split plot by patch
##' @param show_ART include/exclude ART coverage
##' @return ggplot
##' @author Pete Dodd
##' @import ggplot2
##' @import scales
##' @import data.table
##' @export
plot_HIV_dynamic <- function(Y,
                             start_year,
                             chain_step_num = 1,
                             out_type = "N",
                             by_age = FALSE,
                             by_patch = TRUE,
                             show_ART = TRUE) {
  D <- extract.pops(Y, chain_step_num, out_type)
  dt <- 1 / 12 # monthly
  ## reformat
  if (by_age == FALSE) {
    D <- D[age == "15-49"] # restrict to conventional HIV age group
    if (by_patch == TRUE) {
      X <- data.table::dcast(D, patch + step ~ hiv, fun.aggregate = sum)
      vars <- c("HIVpc", "ARTpc", "patch", "step")
      id <- c("patch", "step")
    } else {
      X <- data.table::dcast(D, step ~ hiv, fun.aggregate = sum)
      vars <- c("HIVpc", "ARTpc", "step")
      id <- c("step")
    }
  } else {
    if (by_age == TRUE) {
      X <- data.table::dcast(D, patch + step + age ~ hiv, fun.aggregate = sum)
      id <- c("patch", "step", "age")
      vars <- c("HIVpc", "ARTpc", "patch", "step", "age")
    } else {
      stop("by_age must be TRUE or FALSE")
    }
  }
  X[, tot := `HIV-` + `HIV+/ART-` + `HIV+/ART+`]
  X[, HIVpc := (`HIV+/ART-` + `HIV+/ART+`) / tot]
  X[, ARTpc := `HIV+/ART+` / (`HIV+/ART-` + `HIV+/ART+` + 1e-10)]
  if (show_ART == FALSE) {
    X[, ARTpc := NA_real_]
  }
  X <- data.table::melt(X[, ..vars], id = id)
  out_lab <- out_type # default
  if (out_type == "N") {
    out_lab <- "Population"
  } else {
    if (out_type == "incidence") {
      out_lab <- "Incidence"
    } else {
      if (out_type == "notifrate") {
        out_lab <- "Notification rate"
      }
    }
  }
  ## plot
  GP <- ggplot2::ggplot(X, aes(step * dt + start_year, value, col = variable)) +
    ggplot2::geom_line() +
    ggplot2::scale_y_continuous(label = scales::percent) +
    ggplot2::theme_linedraw() +
    ggplot2::xlab("Year") +
    ggplot2::ylab(paste(out_lab, "(proportion)")) +
    ggplot2::theme(
      legend.position = "top",
      legend.title = element_blank()
    ) +
    ggplot2::scale_color_discrete(labels = c("HIV+/ART-", "HIV+/ART+"))
  if (by_patch == TRUE) {
    wrap_vars <- id[-which(id == "step")]
    GP <- GP + ggplot2::facet_wrap(wrap_vars)
  }
  if (show_ART == FALSE) {
    GP <- GP + ggplot2::theme(legend.position = "none")
  }
  GP
}



## TB incidence by time
#' @title Plot sample TB incidence
#' @param Y input data
#' @param chain_step_num chain step number
#' @param out_type output type
#' @param separate separate plot?
#' @param by_age by age? (default FALSE)
#' @param by_hiv by HIV? (default FALSE)
#' @param wrap display with facet wrapping (default TRUE) or grids
#' @param start_year start for scale (default 2015)
#' @return ggplot2
#' @import data.table
#' @importFrom ggh4x facet_grid2
#' @importFrom ggh4x facet_wrap2
#' @export
plot_TB_dynamics <- function(Y, chain_step_num = 1,
                             out_type = "incidence", separate = FALSE,
                             by_age = FALSE, by_HIV = FALSE, wrap = TRUE,
                             start_year = 2015) {
  ## extract relevant data
  D <- extract.pops(Y, chain_step_num, out_type)

  ## ## version with multi
  ## D <- extract.pops.multi( Y, dim(Y)[2], out_type )
  ## names(D)[names(D)==out_type] <- "value"
  ## D[,varname:=out_type]
  ## D[,step:=t]
  ## eps <- 0.25
  ## D <- D[ ,.( mid=median( noterate ),
  ##            lo = quantile( noterate, eps ),
  ##            hi = quantile( noterate, 1-eps )),
  ##        by=.(t, patch )]
  ## TODO this will need adaptation of plotting code below (see noterate for eg)

  by_comp <- "patch"
  dt <- 1 / 12
  D <- D[step > 1]
  D[, yr := step * dt + start_year]

  if (separate == FALSE) {
    if (by_HIV) {
      D[, hiv := ifelse(hiv == "HIV-", "HIV-", "HIV+")]
      D <- D[, .(value = sum(value)), by = .(yr, patch, age, hiv)]
    } else {
      D <- D[, .(value = sum(value)), by = .(yr, patch, age)]
    }
  }
  if (by_age == FALSE) {
    if (by_HIV) {
      D <- D[, .(value = sum(value)), by = .(yr, patch, hiv)]
    } else {
      D <- D[, .(value = sum(value)), by = .(yr, patch)]
    }
  } else {
    by_comp <- c(by_comp, "age")
  }

  if (by_HIV) {
    plt <- ggplot2::ggplot(D, aes(yr, value, col = hiv))
  } else {
    plt <- ggplot2::ggplot(D, aes(yr, value))
  }
  plt <- plt +
    ggplot2::geom_line() +
    ggplot2::theme_light() +
    ggplot2::xlab("Year") +
    ggplot2::ylab(out_type) +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggplot2::expand_limits(y = c(0, NA))
  if (wrap) {
    plt <- plt + ggh4x::facet_wrap2(by_comp,
      scales = "free"
    )
  } else {
    plt <- plt + ggh4x::facet_grid2(by_comp,
      independent = TRUE,
      scales = "free"
    )
  }
  plt
  e
}



## Note - this may need to be made more flexible to align with
## different start-years for real-data comparison
## Plots HIV prevalence among TB notes -
## Choose % or per 100,000
## Choose to separate HIV into ART+/- or not
## Currently works for single particle/chain step.... add variation over chain??
#' @title Plot HIV prevalence among TB notifications#'
#' @param Y input data
#' @param chain_step_num chain step
#' @param separate separate HIV by ART status
#' @importFrom dplyr filter group_by summarise ungroup
#' @importFrom tidyr complete
#' @importFrom plyr mapvalues
#' @importFrom ggh4x facet_grid2
#' @return ggplot2
#' @export
HIV_prev_TB <- function(Y, n_chain_steps = 100, separate = FALSE) {
  D <- extract.pops.multi(Y, n_chain_steps, out_type = "notes")
  real_dat <- BLASTtbmod::TB_notes_HIV_patch

  if (separate == FALSE) {
    aggD$hiv <- plyr::mapvalues(aggD$hiv,
      from = c("HIV+/ART-", "HIV+/ART+"),
      to = c("HIV+", "HIV+")
    )
    real_dat$hiv <- plyr::mapvalues(real_dat$hiv,
      from = c("HIV+/ART-", "HIV+/ART+"),
      to = c("HIV+", "HIV+")
    )
    replen <- 2
  } else {
    replen <- 3
  }

  # Aggregate over appropriate compartments
  aggD <- D %>%
    dplyr::group_by(patch, hiv, chain_step) %>%
    dplyr::summarise(tot_notes = sum(notes)) %>%
    data.table::as.data.table()

  # Calculate prevalence
  D_notes <- aggD %>%
    dplyr::group_by(patch, chain_step) %>%
    dplyr::summarise(grandtot = sum(tot_notes, na.rm = T))
  D_out <- merge(aggD, D_notes, by = c("patch", "chain_step"), all.x = T)
  D_out$prev_pc <- D_out$tot_notes / D_out$grandtot * 100

  # Generate ribbon
  eps = 0.25
  D_out <- D_out[, .(
    midpc = median(prev_pc),
    lopc = quantile(prev_pc, eps),
    hipc = quantile(prev_pc, 1 - eps)
  ),
  by = .(hiv, patch)
  ]

  # Wrangle real data
  real_out <- real_dat %>%
    dplyr::filter(total > 0) %>%
    dplyr::group_by(patch, hiv) %>%
    dplyr::summarise(tot = sum(total)) %>%
    dplyr::ungroup() %>%
    tidyr::complete(patch, hiv, fill = list(sum_quantity = 0))
  real_notes <- real_out %>%
    dplyr::group_by(patch) %>%
    dplyr::summarise(grandtot = sum(tot, na.rm = T))
  real_out <- merge(real_out, real_notes, by = "patch", all.x = T)
  real_out$prev_pc_real <- real_out$tot / real_out$grandtot * 100

  D_out <- merge(D_out, real_out[, c("patch", "hiv", "prev_pc_real")],
    by = c("patch", "hiv"),
  )
  D_out$patch <- substr(D_out$patch, 7, 7)
  D_out$pad_hi <- sapply(1:nrow(D_out), function(x) {
    min(
      D_out$hipc[x] + 5,
      D_out$prev_pc_real + 5, 100
    )
  })
  D_out$pad_lo <- sapply(1:nrow(D_out), function(x) {
    max(
      D_out$lopc[x] - 5,
      D_out$prev_pc_real - 5, 0
    )
  })


  p <- ggplot2::ggplot(D_out, aes(x = patch, y = midpc)) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(aes(ymin = lopc, ymax = hipc), width = 0.4) +
    ggplot2::geom_blank(aes(x = patch, y = pad_lo)) +
    ggplot2::geom_blank(aes(x = patch, y = pad_hi)) +
    ggh4x::facet_grid2(~hiv,
      independent = T,
      scales = "free"
    ) +
    ggplot2::geom_point(aes(y = prev_pc_real), col = 2, shape = 16) +
    ggplot2::ylab("HIV status among TB notifications (%)") +
    ggplot2::xlab("Patch") +
    ggplot2::ggtitle("Real data (red), median + 50% PI simulated (black)")
  print(p)
}
