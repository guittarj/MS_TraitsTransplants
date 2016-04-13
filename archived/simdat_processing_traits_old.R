# Processing the simulation data produced by 'neutral_simulation.R' and 
# summarized by 'simdat_processing_veg.R'

# There are two ways of estimate parameter values. One way is by using the 
# 'bayesian_immigration_estimates.R' script. Another way is to choose best-fit parameters to 
# observed rates of vegetation change (currently performed in TraitsTransplants.Rmd -> "rates")

# Calculate trait distances in simulations that used desired parameters
if (pars.method == 'max') {
  wd.simdat <- 'E:\\sim_max'
  setwd(wd)
  simdat.veg <- fread("data\\simSummary_max_veg.csv")
  pars <- pars.max  
}

if (pars.method == 'bayes') {
  wd.simdat <- 'G:\\sim_bayes'
  setwd(wd)
  simdat.veg <- fread("data\\simSummary_bayes_veg.csv", drop = 1)
  pars <- pars.bayes
}

setwd(wd.simdat)
sim.list <- list.files()
iwant <- with(pars, paste(site, d, m))
dist.traits <- trait.list3

# Filter out simulations with unwanted parameters. First at parameter level
tmp <- plyr::ldply(strsplit(sim.list, '_'))

# sometimes necessary to fix units
tmp <- transmute(tmp, m = as.numeric(V2), d = as.integer(V4), rep = as.integer(V6))
sim.list <- sim.list[paste(tmp$d, tmp$m) %in% paste(pars$d, pars$m)]

# Filter out simulations with unwanted parameters at the turf*parameter resolution
simsites <- cover.meta$destSiteID[match(simdat.veg$turfID, cover.meta$turfID)]
simdat <- simdat.veg[paste(simsites, simdat.veg$m, simdat.veg$d) %in% paste(pars$site, pars$m, pars$d)]

# Calculate distances
for(i in sim.list) {
  for(trait in dist.traits) {
    setwd(wd.simdat)
    comm <- as.matrix(fread(i, sep = ',', header = TRUE, stringsAsFactors = FALSE))
    comm[, 'm'] <- signif(as.numeric(comm[, 'm']), 3)
    comm.sites <- cover.meta$destSiteID[match(comm[, 'turfID'], cover.meta$turfID)]
    comm <- comm[paste(comm.sites, comm[, 'd'], comm[, 'm']) %in% iwant, ]
    simdat <- rbind(simdat, process.comm(comm, trait))
  
    progress <- (match(i, sim.list) - 1) * length(dist.traits) + match(trait, dist.traits)
  
    # Summarize periodially and at the end...
    if(progress %% 50 == 0 | progress == length(sim.list) * length(dist.traits)) {
      simdat <- group_by(simdat, trait, turfID, year, m, d) %>%
                summarise(dist.tt1 = weighted.mean(dist.tt1, w = reps),
                          reps     = sum(reps))
      print(paste(progress, 'of', length(sim.list) * length(dist.traits)))
      flush.console()
    }
  }
}

# Add 2009 data for each simulation run

# finalize
setwd(wd)
if(pars.method == 'bayes') {
  write.csv(simdat, row.names = FALSE, file = 'data\\simSummary_bayes_traits.csv')
}

if(pars.method == 'max') {
  write.csv(simdat, row.names = FALSE, file = 'data\\simSummary_max_traits.csv')
}
