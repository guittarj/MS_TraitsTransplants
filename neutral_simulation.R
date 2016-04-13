# Script for simulation of neutral community dynamics

# -------------------------------------------------------
# User settings

# where should newly generated data go?
# wd.sim <- 'G:\\sim_bayes'
# setwd(wd.sim)
# dir.create(file.path(paste0(wd, wd.sim)))

# Define replacement (d) and immigration (m) rates
# Replacement: The number of 'deaths' and 'births' per simulated year
# d.vals <- seq(0, 100, by = 5)

# Or, pull turnover rates as calculated in TraitsTransplants.Rmd
# d.vals <- round(rates$turnover)

# Immigration: all new individuals will come from turf (m = 0) or site (m = 1).
# m.vals <- seq(0, 1, by = 0.025)

# Or, pull immigration rates as calculated from 'bayesian_immigration_estimate.R'
# m.vals <- signif(m.bayes$m, 3)

# Which rep numbers do you want to perform?
# rep0 <- 1
# rep1 <- 33

# -------------------------------------------------------
# Simulation

# initialize dataframes and ticker
nams <- c('turfID', 'year', 'm', 'd', 'rep')
len <- length(cover.meta$Year[cover.meta$Year != 2009])
sim <- matrix(nrow = len, ncol = length(colnames(cover)), dimnames = list(c(1:len), colnames(cover)))
sim.meta <- matrix(nrow = len, ncol = length(nams), dimnames = list(c(1:len), nams))
sim.empty <- sim
sim.meta.empty <- sim.meta
ticker <- 0

# Define simulation as a function
sim.fun <- function(m, repx, d) {

  # Initialize
  ticker <<- ticker + 1
  comm.row <- 1
  yrs <- c(2009, 2011:2013)

  for(turf in unique(cover.meta$turfID)) {

    # Locate turf position
    turf.pos <- with(cover.meta, c(1:nrow(cover))[turfID == turf & Year == 2009])
    
    # Initialize turf cover, with whole numbers
    simturf <- round(cover[turf.pos, ])

    # Loop over years
    for(yr in c(2009, 2011, 2012)) {
      controls <- with(cover.meta, 
                    colSums(cover[destSiteID == destSiteID[turf.pos] & 
                                  TTtreat  %in% c('TTC', 'TT1') & 
                                  Year       == yr, ]))
      # Perform simulation twice if year is 2009
      if(yr == 2009) { 
        iters <- 2 * d 
      } else {
        iters <- d
      }
      
      # Loop over replacement events ('iters')
      for(iter in 1:iters) {

          #select random species from target ('coming') and turf ('going')...
          if(runif(1) < m) {
            coming <- controls
          } else { 
            coming <- simturf
          }
          coming <- sample(names(coming), 1, prob = coming)
          going  <- sample(names(simturf), 1, prob = simturf)
          
          # Gain/lose cover units
          simturf[going] = simturf[going] - 1
          simturf[coming] = simturf[coming] + 1
      }

      # Perform additions/subtractions to match field cover
      yr1 <- yrs[match(yr, yrs) + 1]
      ideal.cover <- round(sum(cover[cover.meta$Year == yr1 & cover.meta$turfID == turf, ]))
      diff <- ideal.cover - sum(simturf)
      if (abs(diff) > 0) {
        for (i in 1:abs(diff)) {
          if (diff > 0) {
            if(runif(1) < m) {
              coming <- controls 
            } else {
              coming <- simturf
            }
            coming <- sample(names(coming), 1, prob = coming)
            simturf[coming] = simturf[coming] + 1
          } else {
            going <- sample(names(simturf), 1, prob = simturf)
            simturf[going] = simturf[going] - 1
      }}}

      sim[comm.row, ] <- simturf
      sim.meta[comm.row, ] <- c(turf, yr1, m, d, repx)
      comm.row <- comm.row + 1
    }
  }
  write.table(data.frame(sim.meta, sim), file = paste('m', m, 'd', d, 'rep', repx, '.csv', sep = '_'), 
                         sep       = ',', 
                         row.names = FALSE)
  sim <- sim.empty
  sim.meta <- sim.meta.empty

  # progress report
  print(ticker / (length(m.vals) * length(d.vals) * (rep1 - rep0)))
  flush.console()
}

# -----------------------------------------------------------

# Loop through parameter combinations, perform simulations
for (i in rep0:rep1) {
    for(d in d.vals) {
        for(m in m.vals) {
        sim.fun(d = d,
                m = m, 
                repx = i
)}}}

setwd(wd)