# --------------------------------------------------------------------------------------------------

# A script that:
#   1. loads spreadsheets of raw trait data
#   2. connects to seedclim db and loads cover data

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
leda.sla       <- read.csv('LEDA_SLA_RawData.csv',       header = TRUE)
leda.leaf.area <- read.csv('LEDA_leafSize_RawData.csv',  header = TRUE)
leda.seed.mass <- read.csv('LEDA_SeedMass_RawData.csv',  header = TRUE)
norflor.height <- read.csv('norflorheight.csv',          header = TRUE)
clopla.root    <- read.csv('CLOPLA.csv',                 header = TRUE)

# Filter / adjust units / arrange trait data
sla <- rbind(
  my.sla %>%
    transmute(
      source  = 'my',
      species = probfixes(Species, names(probs), probs),
      val     = SLA..m2.kg.1.),
  leda.sla %>%
    filter(general.method != 'laboratory/greenhouse/garden experiment') %>%
    transmute(
      source = 'leda',
      species = probfixes(SBS.name, names(probs), probs),
      val = Single.Value..LEDA.))

sla <- group_by(sla, source, species) %>% 
  summarise(val = mean(val))

tmp <- spread(sla, source, val) %>% filter(!is.na(my) & !is.na(leda))

print(paste("The Pearson Correlation between my SLA trait values and LEDA SLA trait values is", 
  signif(cor(tmp$my, tmp$leda), 3)))

sla <- sla %>% group_by(species) %>% summarise(sla = mean(val))

leaf.area <- rbind(
  transmute(my.leaf.area,
    source  = 'my',
    species = probfixes(Species, names(probs), probs),
    val     = LeafArea * 100),
  filter(leda.leaf.area, general.method != 'laboratory/greenhouse/garden experiment') %>%
  transmute(source  = 'leda',
    species = probfixes(SBS.name, names(probs), probs),
    val     = single.value..mm.2.) %>%
                     filter(val < 5000 & val != 0 & !is.na(val)))

leaf.area <- group_by(leaf.area, source, species) %>% summarise(val = mean(val))

tmp <- spread(leaf.area, source, val) %>% filter(!is.na(leda) & !is.na(my))
print(paste("The Pearson Correlation between my SLA leaf area trait values and 
  LEDA leaf area trait values is", signif(cor(tmp$my, tmp$leda), 3)))

leaf.area <- leaf.area %>% group_by(species) %>% summarise(leaf.area = mean(val))

seed.mass <- leda.seed.mass %>%
  transmute(
    species = probfixes(SBS.name, names(probs), probs),
    val     = single.value..mg.) %>%
  filter(val < 50 & val > 0) %>%
  group_by(species) %>%
  summarise(seed.mass = mean(val, na.rm = TRUE))

max.height <- norflor.height %>%
  transmute(species = probfixes(speciesName, names(probs), probs), max.height = Max.height / 100) %>%
  filter(max.height < 2)

clopla <- clopla.root %>%
  mutate(species = probfixes(Species, names(probs), probs)) %>%
  filter(species %in% splist$species)

clopla <- clopla %>%
  mutate(
    conper = c('<2', '<2', '2+')[match(ConnectionPersist_1, c('1', '2', '>2'))],
    lat = c('<0.01', rep('0.01+', 3))[match(LateralSpread_1, 
      c('<0.01', '0.01-0.25', '>0.25', 'dispersable'))],
    offs = c(rep('0-1', 3), '2+', '2+')[match(OffspringPerParent_1, 
      c('<1', '1', '1 to 1', '2 to 10', '>10'))]) %>%
  melt(id.vars = 'species') %>%
  filter(!is.na(value) & !value == '') %>%
  group_by(species, variable) %>%
  summarise(value = names(which.max(table(value)))) %>%
  spread(variable, value, convert = TRUE)

# turn bud number into a numeric rank
clopla.buds <- clopla %>%
  select(species, BudBank_no_.10cm, BudBank_no_.0_to_10cm, BudBank_no_0cm, BudBank_no_0_to_.10cm, 
    BudBank_no_..10cm) %>%
  gather(bank, buds, BudBank_no_.10cm:BudBank_no_..10cm) %>%
  filter(buds != '') %>%
  mutate(buds = c(0:2)[match(buds, c('0', '1 to 10', '>10'))]) %>%
  group_by(species) %>%
  summarise(buds = sum(buds))

# filter and convert clonal traits to binaries (high/low)
clopla <- clopla %>%
  filter(species %in% splist$species) %>%
  transmute(species,
    conper = ifelse(conper == '2+', 1, 0), 
    offs   = ifelse(offs == '2+', 1, 0), 
    lat    = ifelse(lat == '0.01+', 1, 0),
    buds   = clopla.buds$buds[match(species, clopla.buds$species)])

# combine
trait.data <- list(splist, sla, leaf.area, seed.mass, max.height, clopla)
trait.data <- Reduce(function(...) merge(..., all = TRUE), trait.data)
trait.data <- trait.data %>%
  filter(!is.na(speciescode)) %>%
  mutate(
    sla        = log10(sla),
    leaf.area  = log10(leaf.area),
    seed.mass  = log10(seed.mass),
    max.height = log10(max.height))

# --------------------------------------------------------------------------------------------------
# Connect to database

con <- odbcConnect('seedclimDB') 

# Cover data
coverQ <- "SELECT sites.siteID, blocks.blockID, turfs.TTtreat,turfs.turfID, dest_blocks.blockID AS 
destBlockID, (SELECT Count(subTurfEnvironment.bad) AS CountOfbad FROM subTurfEnvironment where 
(subTurfEnvironment.year = new_TurfCommunity.year) AND (subTurfEnvironment.turfID = 
new_TurfCommunity.turfID) AND ( (subTurfEnvironment.bad)='')) AS notbad, sites.Temperature_level, 
sites.Precipitation_level, sites.SummerTemperature_gridded, sites.Annualprecipitation_gridded, 
new_TurfCommunity.Year, new_TurfCommunity.species, new_TurfCommunity.cover, turfEnvironment.recorder, 
dest_blocks.siteID as destSiteID FROM (((blocks AS dest_blocks INNER JOIN plots AS dest_plots ON 
dest_blocks.blockID = dest_plots.blockID) INNER JOIN (((sites INNER JOIN blocks ON sites.siteID = 
blocks.siteID) INNER JOIN plots ON blocks.blockID = plots.blockID) INNER JOIN turfs ON plots.plotID 
= turfs.originPlotID) ON dest_plots.plotID = turfs.destinationPlotID) INNER JOIN new_TurfCommunity 
ON turfs.turfID = new_TurfCommunity.turfID) INNER JOIN turfEnvironment ON (turfEnvironment.year = 
new_TurfCommunity.Year) AND (turfs.turfID = turfEnvironment.turfID) WHERE NOT turfs.TTtreat='' AND 
((Not (new_TurfCommunity.Year)=2010));"

cover.thin <- sqlQuery(con, coverQ)
                                       
#correct for stomping
stompingQ <- "SELECT blocks.siteID, blocks.blockID, turfs.turfID, subTurfEnvironment.year, 
turfs.TTtreat, Count(subTurfEnvironment.bad) AS CountOfbad FROM blocks INNER JOIN (plots INNER JOIN 
(turfs INNER JOIN subTurfEnvironment ON turfs.turfID = subTurfEnvironment.turfID) ON plots.plotID = 
turfs.destinationPlotID) ON blocks.blockID = plots.blockID GROUP BY blocks.siteID, blocks.blockID, 
turfs.turfID, subTurfEnvironment.year, turfs.TTtreat, subTurfEnvironment.bad HAVING 
(((subTurfEnvironment.bad)='x'));"

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

odbcCloseAll()

# Disconnect from database
# --------------------------------------------------------------------------------------------------

# Save unpruned versions
cover.all <- cover[, !names(cover) %in% c('NID.gram', 'NID.herb', 'NID.rosett', 'NID.seedling')]
trait.data.all <- trait.data[trait.data$speciescode %in% names(cover.all), ]

# now filter trees/shrubs, 'sp' ambiguities when other species in that genus are present
trait.data <- trait.data %>%
  filter(!speciescode %in% c('Caru.car', 'Pic.abi', 'Bet.pub', 'Jun.com')) %>% 
  filter(!speciescode %in% c('Car.sp', 'Dry.sp', 'Hie.sp', 'Sal.sp', 'Ver.sp', 'Tar.sp'))

# Fill in missing data with good guesses (this might be negated later, as all 'Alc.alp' -> 'Alc.sp')
trait.data <- within(trait.data, {
  seed.mass[speciescode == 'Alc.sp'] <- seed.mass[speciescode == 'Alc.alp']
  max.height[speciescode == 'Alc.sp'] <- max.height[speciescode == 'Alc.alp']
})

# Now prune entries without trait data
cover <- cover.all[, names(cover.all) %in% trait.data$speciescode] # gets rid of UIDs
trait.data <- trait.data[trait.data$speciescode %in% names(cover), ]

#verify again
all(with(cover.meta, paste(turfID, Year, sep = '_')) == rownames(cover))

# save and finish
write.csv(cover, file = 'cover.csv')
write.csv(cover.all, file = 'cover_all.csv')
write.csv(cover.meta, row.names = FALSE, file = 'covermeta.csv')
write.csv(trait.data, row.names = FALSE, file = 'traitdata.csv')
