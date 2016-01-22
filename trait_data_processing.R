# --------------------------------------------------------------------------------------------------

# A script that:
#   1. loads spreadsheets of raw trait data
#   2. connects to seedclim db and loads cover data
#   3. At one point it also loaded additional trait data (e.g. habitat, etc.), but not anymore

# Notes: 
# make sure the DB is running by starting it in: 
#   -  control Panel\System and Security\Administrative Tools\ODBC Data Sources (32-bit)
#   -  i tried the 64-bit driver version of ODBC and it doesn't work 
#   -  the DB is called 'seedclimDB' in this case
#   -  make sure you are running 32-bit R 

# --------------------------------------------------------------------------------------------------

# Set working directory if necessary
wd <- 'C:\\Users\\John\\Google Drive\\Documents\\michigan\\seedclim\\data'
setwd(wd)

# Load custom functions
source('C:\\Users\\John\\Google Drive\\Documents\\michigan\\seedclim\\MS_TraitsTransplants\\custom_functions.R')

# Load packages
loadpax(pkg = c('RODBC', 'dplyr', 'tidyr', 'reshape2'))

# --------------------------------------------------------------------------------------------------

# Create species list

splist <- read.csv('SpeciesList2013.csv', header = TRUE, stringsAsFactors = FALSE)
splist <- transmute(splist, speciescode      = species,
                            species          = speciesName,
                            family           = family,
                            graminoid        = ifelse(functionalGroup == 'graminoid', 1, 0),
                            perennial        = ifelse(lifeSpan == 'perennial', 1, 0))

# Load naming problem lists, correct names
source("C:\\Users\\John\\Google Drive\\Documents\\michigan\\seedclim\\MS_TraitsTransplants\\naming_problems.R")
splist$species <-  probfixes(splist$species, names(probs), probs)
  
# --------------------------------------------------------------------------------------------------

# Load trait data
my.sla         <- read.csv('RawTraitData_SLA.csv',       header = TRUE)
my.leaf.area   <- read.csv('RawTraitData_leafsize.csv',  header = TRUE)
my.leaf.chem   <- read.csv('raw_data_CN_2014Sept15.csv', header = TRUE)
leda.sla       <- read.csv('LEDA_SLA_RawData.csv',       header = TRUE)
leda.leaf.area <- read.csv('LEDA_leafSize_RawData.csv',  header = TRUE)
leda.seed.mass <- read.csv('LEDA_SeedMass_RawData.csv',  header = TRUE)
leda.ldmc      <- read.csv('LEDA_LDMC_RawData.csv',      header = TRUE)
norflor.height <- read.csv('norflorheight.csv',          header = TRUE)
clopla.root    <- read.csv('CLOPLA.csv',                 header = TRUE)

# Filter / adjust units / arrange trait data
sla <- rbind(my.sla %>%
               transmute(source  = 'my',
                         species = probfixes(Species, names(probs), probs),
                         val     = SLA..m2.kg.1.),
             leda.sla %>% 
               filter(general.method != 'laboratory/greenhouse/garden experiment') %>%
               transmute(source = 'leda',
                         species = probfixes(SBS.name, names(probs), probs),
                         val = Single.Value..LEDA.)) %>%
       group_by(source, species) %>%
       summarise(val = mean(val)) %>%
       group_by(species) %>%
       summarise(sla = mean(val, na.rm = T))

leaf.area <- rbind(my.leaf.area %>%
                     transmute(source  = 'my',
                               species = probfixes(Species, names(probs), probs),
                               val     = LeafArea * 100),
                   leda.leaf.area %>%
                     filter(general.method != 'laboratory/greenhouse/garden experiment') %>%
                     transmute(source  = 'leda',
                               species = probfixes(SBS.name, names(probs), probs),
                               val     = single.value..mm.2.) %>%
                     filter(val < 5000 & val != 0 & !is.na(val))) %>%
             group_by(source, species) %>%
             summarise(val = mean(val)) %>%
             group_by(species) %>%
             summarise(leaf.area = mean(val))

ldmc <- mutate(leda.ldmc, species = probfixes(species, names(probs), probs)) %>%
        group_by(species) %>%
        summarise(ldmc = weighted.mean(ldmc, w = n))

leaf.chem <- transmute(my.leaf.chem, species = splist$species[match(Species, splist$speciescode)],
                                     leafCN  = CN,
                                     d13C    = d13C.UCD,
                                     d15N    = d15N.UCD) %>%
             group_by(species) %>%
             summarise(leafCN = mean(leafCN, na.rm = TRUE),
                       d13C   = mean(d13C,   na.rm = TRUE),
                       d15N   = mean(d15N,   na.rm = TRUE))

seed.mass <- transmute(leda.seed.mass, species = probfixes(SBS.name, names(probs), probs),
                                       val     = single.value..mg.) %>%
             filter(val < 50 & val > 0) %>%
             group_by(species) %>%
             summarise(seed.mass = mean(val, na.rm = TRUE))

max.height <- transmute(norflor.height, species    = probfixes(speciesName, names(probs), probs),
                                        max.height = Max.height / 100) %>%
              filter(max.height < 2)

clopla <- mutate(clopla.root, species = probfixes(Species, names(probs), probs)) %>%
          filter(species %in% splist$species)
clopla <- mutate(clopla,
            cgo     = CGO_Type_1,
            cyc     = c('1', '1', '2+', '2+')[match(ShootCycl_1, 
                      c('1', '1, 2', '2', '>2'))],
            conper  = c('<2', '<2', '2+')[match(ConnectionPersist_1, 
                      c('1', '2', '>2'))],
            offs    = c(rep('0-1', 3), '2+', '2+')[match(OffspringPerParent_1, 
                      c('<1', '1', '1 to 1', '2 to 10', '>10'))],
            lat     = c('<0.01', rep('0.01+', 3))[match(LateralSpread_1, 
                      c('<0.01', '0.01-0.25', '>0.25', 'dispersable'))]) %>%
          melt(id.vars = 'species') %>%
          filter(!is.na(value) & !value == '') %>%
          group_by(species, variable) %>%
          summarise(value = names(which.max(table(value)))) %>%
          spread(variable, value, convert = TRUE)

clopla.buds <- select(clopla, species, BudBank_no_.10cm, BudBank_no_.0_to_10cm,
                 BudBank_no_0cm, BudBank_no_0_to_.10cm, BudBank_no_..10cm) %>%
               gather(bank, buds, BudBank_no_.10cm:BudBank_no_..10cm) %>%
               filter(buds != '') %>%
               mutate(buds = c(0:2)[match(buds, c('0', '1 to 10', '>10'))]) %>%
               group_by(species) %>%
               summarise(buds = sum(buds))

clopla <- transmute(clopla, species, cgo, cyc, conper, offs, lat, 
            buds = clopla.buds$buds[match(species, clopla.buds$species)]) %>%
          filter(species %in% splist$species)

# Do I want to convert clonal traits to a binary (buds remain numeric)
clopla <- mutate(clopla, 
       cyc = ifelse(cyc == '2+', 1, 0),
       conper = ifelse(conper == '2+', 1, 0),
       offs = ifelse(offs == '2+', 1, 0),
       lat = ifelse(lat == '0.01+', 1, 0))

# combine
trait.data <- list(splist, sla, leaf.area, ldmc, leaf.chem, seed.mass, max.height, clopla)
trait.data <- Reduce(function(...) merge(..., all = TRUE), trait.data)
trait.data <- filter(trait.data, !is.na(speciescode)) %>%
              mutate(sla        = log10(sla),
                     leaf.area  = log10(leaf.area),
                     ldmc       = log10(ldmc),
                     leafCN     = log10(leafCN),
                     seed.mass  = log10(seed.mass),
                     max.height = log10(max.height))
# --------------------------------------------------------------------------------------------------

# Connect to database
con <- odbcConnect('seedclimDB') 

# Cover data
coverQ <- "SELECT sites.siteID, blocks.blockID, turfs.TTtreat,turfs.turfID, dest_blocks.blockID AS destBlockID, (SELECT Count(subTurfEnvironment.bad) AS CountOfbad
FROM subTurfEnvironment where (subTurfEnvironment.year = new_TurfCommunity.year) AND (subTurfEnvironment.turfID = new_TurfCommunity.turfID)
 AND ( (subTurfEnvironment.bad)='')) AS notbad, sites.Temperature_level, sites.Precipitation_level, sites.SummerTemperature_gridded, sites.Annualprecipitation_gridded, new_TurfCommunity.Year, new_TurfCommunity.species, new_TurfCommunity.cover, turfEnvironment.recorder , dest_blocks.siteID as destSiteID
FROM (((blocks AS dest_blocks INNER JOIN plots AS dest_plots ON dest_blocks.blockID = dest_plots.blockID) INNER JOIN (((sites INNER JOIN blocks ON sites.siteID = blocks.siteID) INNER JOIN plots ON blocks.blockID = plots.blockID) 
INNER JOIN turfs ON plots.plotID = turfs.originPlotID) ON dest_plots.plotID = turfs.destinationPlotID) INNER JOIN new_TurfCommunity ON turfs.turfID = new_TurfCommunity.turfID) INNER JOIN turfEnvironment ON (turfEnvironment.year = new_TurfCommunity.Year) AND (turfs.turfID = turfEnvironment.turfID)
WHERE NOT turfs.TTtreat='' AND ((Not (new_TurfCommunity.Year)=2010));"

cover.thin <- sqlQuery(con, coverQ)
                                       
#correct for stomping
stompingQ <- "SELECT blocks.siteID, blocks.blockID, turfs.turfID, subTurfEnvironment.year, turfs.TTtreat, Count(subTurfEnvironment.bad) AS CountOfbad
FROM blocks INNER JOIN (plots INNER JOIN (turfs INNER JOIN subTurfEnvironment ON turfs.turfID = subTurfEnvironment.turfID) ON plots.plotID = turfs.destinationPlotID) ON blocks.blockID = plots.blockID
GROUP BY blocks.siteID, blocks.blockID, turfs.turfID, subTurfEnvironment.year, turfs.TTtreat, subTurfEnvironment.bad
HAVING (((subTurfEnvironment.bad)='x'));"

#delete turfs with too much stomping  
cover.thin <- cover.thin[cover.thin$notbad > 10, ]
                                 
#correct covers for stomping
cover.thin$cover <- cover.thin$cover * 25 / cover.thin$notbad
cover.thin$cover[cover.thin$cover > 80] <- 80 #stop doubtfully high values                                 

#correct for botanist effects
#PM
with(cover.thin, cover[recorder == 'PM'] <- cover[recorder == 'PM'] * 1.20)

#Siri
siriQ <- "SELECT turfs.turfID, new_TurfCommunity.Year, turfEnvironment.date, 
Sum(new_TurfCommunity.cover) AS SumOfcover, turfEnvironment.totalVascular, 
turfs.TTtreat, sites.Temperature_level FROM ((sites INNER JOIN blocks ON 
sites.siteID = blocks.siteID) INNER JOIN plots ON blocks.blockID = plots.blockID) 
INNER JOIN ((turfs INNER JOIN turfEnvironment ON turfs.turfID = 
turfEnvironment.turfID) INNER JOIN new_TurfCommunity ON (new_TurfCommunity.Year 
= turfEnvironment.year) AND (turfs.turfID = new_TurfCommunity.turfID)) 
ON plots.plotID = turfs.destinationPlotID WHERE (((turfEnvironment.recorder) = 
'Siri')) GROUP BY turfs.turfID, new_TurfCommunity.Year, turfEnvironment.date, 
turfEnvironment.totalVascular, turfs.TTtreat, sites.Temperature_level HAVING 
((Not (turfs.TTtreat)='')) ORDER BY new_TurfCommunity.Year, turfEnvironment.date, 
Sum(new_TurfCommunity.cover) DESC;"

siri <- sqlQuery(con, siriQ)
siriLOW<-siri[siri$SumOfcover / siri$totalV < 1.35, ]
siriLOW$turfID <- as.character(siriLOW$turfID)
siri.fix <- paste(as.character(cover.thin$turfID), cover.thin$Year) %in% 
  paste(siriLOW$turfID,siriLOW$Year)
cover.thin$cover[siri.fix] <- cover.thin$cover[siri.fix] * 1.3

# make fat table
cover <- xtabs(cover ~ paste(turfID, Year, sep = '_') + species, data = cover.thin)
cover <- as.data.frame(unclass(cover))
dim(cover)

#make meta data
cover.meta <- unique(cover.thin[, c('siteID', 'TTtreat', 'Year', 'blockID', 'turfID', 
  'Temperature_level', 'Precipitation_level', 'Annualprecipitation_gridded',
  'SummerTemperature_gridded', 'destBlockID', 'recorder', 'notbad', 'destSiteID')])
cover.meta <- cover.meta[order(paste(cover.meta$turfID, cover.meta$Year)), ]
cover.meta$TTtreat <- factor(as.character(cover.meta$TTtreat), 
  levels=c('TTC','TT1', 'TT2', 'TT3', 'TT4'))

#fix remaining data entry problems!?
cover['111 TT2 137_2011', 'Agr.cap'] <- 25
cover['32 TT3 109_2009',] <- cover['32 TT3 109_2009',] / 2
cover['32 TT3 109_2012',] <- cover['32 TT3 109_2012',] * 2 / 3
cover['33 TT2 58_2009',] <- cover['33 TT2 58_2009',] * 2 / 3
cover['34 TT1 32_2009',] <- cover['34 TT1 32_2009',] / 2
cover['40 TT2 62_2011',] <- cover['40 TT2 62_2011',] * 2 / 3 

#fix wonky row numbers
row.names(cover.meta) <- c(1:length(cover.meta[, 1]))

#verify
cover.meta$turf.year <- paste(cover.meta$turfID, cover.meta$Year, sep = '_')
all(cover.meta$turf.year == rownames(cover))

#ensure there are four years for each turf...
tmp <- as.data.frame(table(cover.meta$turfID))
tmp <- with(tmp, cover.meta[cover.meta$turfID %in% Var1[Freq == 4], ])
tmp <- with(tmp, tmp[blockID %in% destBlockID, ]) # remove turfs without destination/origin blocks 
cover <- cover[cover.meta$turfID %in% tmp$turfID, ] #apply the same filters to cover matrix
cover.meta <- tmp

# Load subplot data
subturf <- "SELECT newSubTurfCommunity.turfID, newSubTurfCommunity.subTurf, newSubTurfCommunity.Year, newSubTurfCommunity.species
FROM newSubTurfCommunity;"
subturf <- sqlQuery(con, subturf)
subturf <- subset(subturf, turfID%in%cover.meta$turfID & Year != 2010)
subturf$subturfID <- with(subturf, paste(turfID, subTurf, sep = '_'))
subturf$subturf.year <- with(subturf, paste(turfID, subTurf, Year, sep='_'))
tmp <- dcast(subturf, subturf.year ~ species, length, value.var = 'subturfID')
subturf.meta <- unique(mutate(subturf, species = NULL))
subturf.meta <- cbind(subturf.meta, 
                      cover.meta[match(paste(subturf.meta$turfID,subturf.meta$Year, sep = '_'), 
                                       cover.meta$turf.year), 
                                 !names(cover.meta) %in% names(subturf.meta)])
stopifnot(subturf.meta$subturfyear %in% tmp$subturf.year)
subturf.meta <- subturf.meta[match(tmp$subturf.year, subturf.meta$subturf.year), ]
row.names(tmp) <- tmp$subturf.year
tmp$subturf.year <- NULL
subturf <- tmp

# replace missing subturfs using random sampling of adjacent subturfs in that year
missing <- with(subturf.meta, expand.grid(unique(turfID), unique(subTurf), unique(Year)))
missing <- rename(missing, turfID = Var1, 
                           subturf = Var2, 
                           year = Var3) %>%
           mutate(subturf.year = paste(turfID, subturf, year, sep='_')) %>%
           filter(! subturf.year %in% subturf.meta$subturf.year)
for(i in 1:nrow(missing)) {
  tmp <- subturf[subturf.meta$turfID == missing$turfID[i] & subturf.meta$Year == missing$year[i], ]
  needed.cover <- round(mean(rowSums(tmp)))
  tmp <- colSums(tmp)
  tmp <- tmp[tmp > 0]
  subturf[missing$subturf.year[i], ] <- 0
  for(j in 1:needed.cover) {
    coming <- sample(names(tmp), 1, prob = tmp)
    subturf[missing$subturf.year[i], coming] <- subturf[missing$subturf.year[i], coming] + 1
  }
  tmp <- subset(subturf.meta, turfID == missing$turfID[i] & subTurf == missing$subturf[i])[1, ]
  subturf.meta <- rbind(subturf.meta,
                        mutate(tmp, Year         = missing$year[i], 
                                    subturf.year = missing$subturf.year[i], 
                                    turf.year    = paste(missing$turfID, missing$year,sep='_')[i]))
}

print(paste(nrow(missing), 'missing subturfs of' , nrow(subturf), 
  'total were created based on adjacent subturf compositions'))

subturf <- subturf[order(row.names(subturf)), ]
subturf.meta <- subturf.meta[order(subturf.meta$subturf.year), ]
stopifnot(identical(row.names(subturf), subturf.meta$subturf.year))

odbcCloseAll()

# First save unpruned versions
cover.all <- cover[, !names(cover) %in% c('NID.gram', 'NID.herb', 'NID.rosett', 'NID.seedling')]
trait.data.all <- trait.data[trait.data$speciescode %in% names(cover.all), ]

trait.data <- filter(trait.data, !speciescode %in% c('Caru.car', 'Pic.abi', 'Bet.pub', 'Jun.com')) %>% # no trees/shrubs
              filter(!speciescode %in% c('Car.sp', 'Dry.sp', 'Hie.sp', 'Sal.sp', 'Ver.sp', 'Tar.sp')) 
                      # no 'sp' ambiguities when other species in that genus are present

# Fill in some missing data with good guesses
trait.data <- within(trait.data, {
  seed.mass[speciescode == 'Alc.sp'] <- seed.mass[speciescode == 'Alc.alp']
  max.height[speciescode == 'Alc.sp'] <- max.height[speciescode == 'Alc.alp']
})


# Remove anything without trait data
cover <- cover.all[, names(cover.all) %in% trait.data$speciescode] # gets rid of UIDs
subturf <- subturf[, colnames(subturf) %in% trait.data$speciescode]
trait.data <- trait.data[trait.data$speciescode %in% names(cover), ]

# Ensure presence/absence data is represented by a '1' or '0'
subturf = apply(subturf, 2, function(x) ifelse(x > 0, 1, 0))

#verify again
all(with(cover.meta, paste(turfID, Year, sep = '_')) == rownames(cover))
all(with(subturf.meta, paste(turfID, subTurf, Year, sep = '_')) == row.names(subturf))

write.csv(cover, file = 'cover.csv')
write.csv(cover, file = 'cover_all.csv')
write.csv(cover.meta, row.names = FALSE, file = 'covermeta.csv')
write.csv(trait.data, row.names = FALSE, file = 'traitdataFull.csv')
write.csv(trait.data.all, row.names = FALSE, file = 'traitdataFull_all.csv')
write.csv(subturf, file = 'subturf.csv')  
write.csv(subturf.meta, row.names = FALSE, file = 'subturfmeta.csv')

# -----------------------------------------------------------------------------
# For some trait x clonal trait comparisons see 'some_trait_comparisons.R' in 
# the 'scripts' folder
# -----------------------------------------------------------------------------

