# A script that estimates immigration rates using JAGS, given the rates of turnover
# we observe in local controls, and the composition of controls at the site level

#load packages
loadpax(c("scales","R2jags", "lattice"))

# round cover to nearest unit
cover0 <- round((cover * 100) / rowSums(cover))

# initialize
dats <- data.frame(site = character(), m = numeric(), rep = integer())
rhat <- data.frame(deviance = numeric(), m = numeric())

# for each site...
for(site in unique(cover.meta$siteID)) {
  # a filter vector for TTCs by site
  filt.flora <- cover.meta$TTtreat == 'TTC' & 
                cover.meta$destSiteID == site &
                cover.meta$Year %in% c(2011:2013)

  # a filter vector for TT1s by site
  filt <- cover.meta$TTtreat == 'TT1' & 
          cover.meta$destSiteID == site &
          cover.meta$Year %in% c(2011:2013)
  
  ## list of site species (using TT1s)
  ## spp <- sort(colSums(cover0[filt, ]), decreasing = TRUE)
  ## spp <- names(spp[spp > 0])
  
  
  # for each species...
  # for(sp in spp) {
  # ticker, if desired
  #print(paste0(site, ': ', sp, ', sp.', match(sp, spp), ' of ', length(spp))) ; flush.console()
  
  # relative abundances in TT1s by year
  N <- as.data.frame(cbind(cover0[filt, ], cover.meta[filt, c('Year','turfID')])) %>%
    gather(sp, abun, -turfID, -Year)
  
  cover1 <- group_by(N, turfID, Year) %>%
    mutate(abun = abun / sum(abun))
  
  N <- spread(N, Year, abun)
  
  # matrix of total cover for each turf * year (using TT1s)
  totcov <- group_by(N, turfID) %>%
    summarise(`2011` = sum(`2011`), `2012` = sum(`2012`), `2013` = sum(`2013`))
  totcov <- totcov[match(N$turfID, totcov$turfID), ]

  cover1 <- spread(cover1, Year, abun) %>%
    select(-turfID, -sp)
  
  # Using mean relative abundances in TTCs for local flora source
  siteflora <- as.data.frame(cbind(cover0[filt.flora, ], cover.meta[filt.flora, c('turfID','Year')])) %>%
    gather(sp, abun, -turfID, -Year) %>%
    group_by(turfID, Year) %>%
    mutate(abun = abun / sum(abun)) %>%
    group_by(sp, Year) %>%
    summarize(abun = mean(abun)) %>%
    spread(Year, abun)
  
  siteflora <- siteflora[match(N$sp, siteflora$sp), -1]
  
  #filter out situations where the spp aren't in the site flora at all (bc errors)
  iwant <- apply(siteflora, 1, function(x) sum(x == 0) == 0)
  N <- as.matrix(N[iwant, -c(1, 2)])
  cover1 <- as.matrix(cover1[iwant, ])
  totcov <- as.matrix(totcov[iwant, -1])
  siteflora <- as.matrix(siteflora[iwant, ])
  
  # determine N
  rowz <- nrow(N)

  # organize into a list
  data <- list("cover1", "totcov", "siteflora", "N", "rowz")

  # Initialize
  inits <- lapply(as.list(c(0.1, 0.9, 0.5)), function(x) list(m = x))
  reps <- length(inits)
  #parameters <- c("m", "predicted_abun")
  parameters <- c('m')
  
  sink(paste0(wd, "\\model.txt"))
  cat("
    model { 
      for(i in 1:rowz) { #i is plot
        for(t in 1:2) {  #t is year
          N[i, t + 1] ~ dpois(lambda[i, t])
          lambda[i, t] <- totcov[i, t] * ((1 - m) * cover1[i, t] + m * siteflora[i, t])
          #predicted_abun[i, t] <- lambda[i, t] # predicted values of abundance for model fit
        }
     }  
      m ~ dunif(0,1)
    }
    ",fill=TRUE)
  sink()
  m1 <- jags(data, inits, parameters, "model.txt", n.thin = 50, n.chains = reps, n.burnin = 1000, n.iter = 10000)
  rhat <- rbind(rhat, m1$BUGSoutput$summary[, 'Rhat'])
  m1.mcmc <- as.mcmc(m1)
  
  mvals <- sapply(m1.mcmc, function(x) mean(x[, 'm']))
  dats <- rbind(dats, data.frame(site, mvals, rep = 1:length(inits)))

}
    
m.bayes <- dats %>% 
  group_by(site) %>% 
  summarise(m = mean(mvals))

# Visualize immigration rates by species as a running average
tmp <- group_by(m.bayes, site) %>%
       mutate(seq    = seq_along(m),
              cumave = cumsum(m) / seq)

figBayesMeans <- ggplot(tmp, aes(x = seq, y = cumave, color = site)) +
geom_line()
x11()
print(figBayesMeans)

pdf("figureBayesMeans.pdf", width = 10, height = 8)
print(figBayesMeans)
dev.off()

flush.console()

# write summary file
m.bayes <- group_by(m.bayes, site) %>%
            summarise(m = mean(m))

write.csv(m.bayes, file = "data\\m.bayes.csv", row.names = FALSE)