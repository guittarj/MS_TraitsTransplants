# Calculating bray-curtis distances from the simulation data produced by 'neutral_simulation.R'

# Overview: I divide the processing into 'chunks' because it takes too long to process
# the whole vegetation survey on one core. I get the list of all files to be processed ('sim.list')
# and then divide that list by the number of chunks I want. I then summarize the mean distances from
# those files and eventually write the chunk summary to the chunk directory. At the end, I merge the
# chunks into a single summary file, and then write that as an .rda (r project file). A .rda uses 
# factors when saving data, thereby saving lots of space, whereas a .csv saves every line as a new 
# string (thus, one survey was over 250mb.

wd.simdat <- 'E:\\sim_max'
setwd(wd.simdat)
sim.list <- list.files()

# divide file list by the number of chunks if desired, to speed up computation
chunk <- 1 # of
chunks <- 4 # total
chunk0 <- round(length(sim.list) / chunks * (chunk - 1))
chunk1 <- round(length(sim.list) / chunks * chunk)
sim.list <- sim.list[chunk0:chunk1]

simdat <- data.frame(trait = factor(), turfID = factor(), year = integer(), m = factor(), 
            d = factor(), dist.tt1 = numeric(), reps = integer())

for (i in sim.list) {
    comm <- as.matrix(fread(i, sep = ',', header = TRUE, stringsAsFactors = FALSE))
    simdat <- rbind(simdat, process.comm(comm, 'veg'))
    progress <- match(i, sim.list)
    
    # Summarize periodially and at the end...
    if(progress %% 50 == 0| progress == length(sim.list)) {
      simdat <- group_by(simdat, trait, turfID, year, m, d) %>%
                summarise(dist.tt1 = weighted.mean(dist.tt1, w = reps),
                          reps     = sum(reps))
      print(paste(progress, 'of', length(sim.list)))
      flush.console()
}}

# write chunks
chunks.dir <- paste0(wd, '\\data\\chunks')
dir.create(chunks.dir)
setwd(chunks.dir)
chunks.dir <- paste0(wd, '\\data\\chunks\\', as.character(Sys.Date()))
dir.create(chunks.dir)
setwd(chunks.dir)
write.csv(simdat, file = paste0('simSummary_veg_chunk', chunk, 'of', chunks, '.csv'))

# merge chunks
setwd(chunks.dir)
chunks <- list.files()
simdat <- plyr::ldply(chunks, fread, sep = ',', header = TRUE, drop = c('V1'))

# Convert to data.table to speed up aggregation
simdat.veg <- data.table(simdat)
simdat.veg <- group_by(simdat.veg, trait, turfID, year, m, d) %>%
              summarise(dist.tt1 = weighted.mean(dist.tt1, w = reps),
                        reps     = sum(reps))

# Write table as an .rda file because a .csv is too big 
setwd(wd)
save(simdat.veg, file = 'data\\simSummary_survey_veg.rda')
