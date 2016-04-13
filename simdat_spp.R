# get simulation files based on bayes immigration estimates
wd.simdat <- 'G:\\sim_bayes'
setwd(wd.simdat)
sim.list <- list.files()

# make list of the wanted parameter combinations 
# (multiplying by 1000 and rounding simplifies everything).
wanted.pars <- paste(round(pars$m * 1000), pars$d)
simIDs <- as.data.frame(do.call('rbind', strsplit(sim.list, "_", fixed = TRUE)))
simIDs <- paste(round(as.numeric(as.character(simIDs$V2)) * 1000), simIDs$V4)

# filter the sim.list accordingly
sim.list <- sim.list[simIDs %in% wanted.pars]

# more filters. only warmer transplants
filt <- filter(cover.meta, Year > 2009)
filt$filt <- filt$TTtreat %in% c('TT2','TT4')
cov.filt <- cover.meta$Year > 2009 & cover.meta$TTtreat %in% c('TT2','TT4')

# prepare dat array and dat.meta
spp <- colnames(cover)
dat <- array(NA, c(sum(cov.filt), length(spp), length(sim.list)), dimnames = list(NULL, spp, NULL))
dat.meta <- matrix(NA, nrow = c(length(sim.list)), ncol = 3)

diff <- NULL
for (i in 1:length(sim.list)) {
  comm <- as.matrix(fread(sim.list[[i]], sep = ',', header = TRUE, stringsAsFactors = FALSE))
  dat.meta[i,] <- as.numeric(c(comm[1, c('m','d','rep')]))
  
  # ensure order
  if (!all(comm[, 'turfID'] == filt$turfID & comm[, 'year'] == filt$Year)) {
      stop("misaligned matrices")
  }
  
  comm <- comm[filt$filt, spp]
  class(comm) <- 'numeric'
  diff[i] <- sum(comm) - sum(cover[cov.filt, ])
  comm <- cover[cov.filt, ] - comm
  comm[cover[cov.filt, ] == 0] <- NA
  dat[,, i] <- comm
  if (i %% 50 == 0) print(signif(i/length(sim.list), 2))
}

# get metadata for each sim
tmp <- select(cover.meta[cov.filt, ], site = destSiteID, treat = TTtreat, year = Year)

# merge dat and metadata
wanted <- pars %>% transmute(id = paste(site, d, round(m * 1000)))
simdat.spp <- list()
for (i in 1:nrow(dat.meta)) {
  j <- as.data.table(cbind(tmp, m = dat.meta[i, 1], d = dat.meta[i, 2],
                    rep = dat.meta[i, 3], dat[, 1:length(spp), i]))
  j <- filter(j, paste(site, d, round(m * 1000)) %in% wanted$id)
  simdat.spp[[i]] <- j
}

# fold it into dataframe
simdat.spp <- rbindlist(simdat.spp)
simdat.spp <- gather_(simdat.spp, "sp", "abun", spp)

setwd(wd)
save(simdat.spp, file = 'data\\simdat_spp.rda')

