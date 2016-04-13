# A script that estimates immigration rates using JAGS, given the rates of turnover
# we observe in local controls, and the composition of controls at the site level

setwd(wd)

#load packages
loadpax(c("scales","R2jags", "lattice"))

# round cover to nearest unit
cover0 <- round((cover * 100) / rowSums(cover))

# initialize
mvals <- data.frame(site = character(), m = numeric(), rep = integer())
posts <- list()
rhat <- data.frame(deviance = numeric(), m = numeric())

# for each site...
for(site in unique(cover.meta$siteID)) {
  # a filter vector for TTCs by site
  site.filt <- cover.meta$TTtreat == 'TTC' & 
               cover.meta$destSiteID == site &
               cover.meta$Year %in% c(2011:2013)

  # a filter vector for TT1s by site
  turf.filt <- cover.meta$TTtreat == 'TT1' & 
               cover.meta$destSiteID == site &
               cover.meta$Year %in% c(2011:2013)
  
  # relative abundances in TT1s by year
  N <- as.data.frame(cbind(cover0[turf.filt, ], cover.meta[turf.filt, c('Year','turfID')]))
  N <- N %>%
    gather(sp, N0, -turfID, -Year) %>%
    group_by(turfID, sp) %>%
    mutate(N1 = ifelse(Year == 2013, NA, 
                ifelse(Year == 2012, N0[Year == 2013], N0[Year == 2012])))
  
  turfabun <- N %>%
    group_by(turfID, Year) %>%
    mutate(abun = N0 / sum(N0))
  
  # matrix of total cover for each turf * year (using TT1s)
  commN <- N %>%
    group_by(turfID, Year) %>%
    summarise(N = sum(N0))
  commN <- commN[match(N$turfID, commN$turfID), ]

  # Using mean relative abundances in TTCs for local flora source
  siteabun <- as.data.frame(cbind(cover0[site.filt, ], cover.meta[site.filt, c('turfID','Year')]))
  siteabun <- siteabun %>%
    gather(sp, abun, -turfID, -Year) %>%
    group_by(turfID, Year) %>%
    mutate(abun = abun / sum(abun)) %>%
    group_by(sp, Year) %>%
    summarize(abun = mean(abun))
  
  siteabun <- siteabun[match(N$sp, siteabun$sp), ]
  
  #filter out situations where the spp aren't in the site flora at all (bc errors)
  filt <- siteabun %>%
    group_by(sp) %>%
    mutate(filt = ifelse(sum(abun == 0) > 0, FALSE, TRUE))
  filt <- filt$filt & !is.na(N$N1)
  
  N1 <- N$N1[filt]
  commN <- commN$N[filt]
  turfabun <- turfabun$abun[filt]
  siteabun <- siteabun$abun[filt]
  N1P <- N1

  # organize into a list
  data <- list('N1','commN','turfabun','siteabun')

  # Initialize
  inits <- lapply(as.list(c(0.1, 0.9, 0.5)), function(x) list(m = x))
  reps <- length(inits)
  parameters <- c('m', 'N1P')
  
  sink(paste0(wd, "\\model.txt"))
  cat("
    model { 
      for (i in 1:length(N1)) {
        N1[i] ~ dpois(commN[i] * ((1 - m) * turfabun[i] + m * siteabun[i]))
      }
      m ~ dunif(0, 1) # nothing changes if m ~ dunif(0, 0.5)
      for (i in 1:length(N1)) {
        zN1P[i] ~ dpois(commN[i] * ((1 - m) * turfabun[i] + m * siteabun[i]))
        N1P[i] <- zN1P[i]
      }
    }
    ", fill=TRUE)
  sink()
  m1 <- jags(data, inits, parameters, "model.txt", n.thin = 50, 
    n.chains = reps, n.burnin = 1000, n.iter = 10000)
  rhat <- rbind(rhat, m1$BUGSoutput$summary[, 'Rhat'])
  m1.mcmc <- as.mcmc(m1)
  
  tmp <- sapply(m1.mcmc, function(x) mean(x[, 'm']))
  mvals <- rbind(mvals, data.frame(site, m = tmp, rep = 1:length(inits)))
  
  tmp <- as.matrix(m1.mcmc[[2]])[, -c(1,2)]
  tmp.ordr <- gsub(pattern = 'N1P\\[', '', colnames(tmp))
  tmp.ordr <- as.numeric(gsub(pattern = '\\]', '', tmp.ordr))
  tmp <- colMeans(tmp[, match(1:ncol(tmp), tmp.ordr)])
  tmp <- data.frame(N1 = N1, N1_predicted = tmp)
  posts[[site]] <- tmp
  
}

rsq <- do.call('rbind', posts)
rsq <- summary(lm(N1 ~ N1_predicted, rsq))$r.squared
rsq

m.bayes <- mvals %>% 
  group_by(site) %>% 
  summarise(m = mean(m))

write.csv(m.bayes, file = "data\\m.bayes.csv", row.names = FALSE)
