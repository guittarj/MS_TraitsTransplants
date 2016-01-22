loadpax(c("scales","R2jags", "lattice"))

# round cover to nearest unit
cover0 <- round((cover * 100) / rowSums(cover))

# initialize
dats <- list()
dats.meta <- data.frame(site = character(), sp = character(), rep = integer())
i <- 1
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
  
  # list of site species (using TT1s)
  spp <- sort(colSums(cover0[filt, ]), decreasing = TRUE)
  spp <- names(spp[spp > 0])
  
  # matrix of total cover for each turf * year (using TT1s)
  totcov <- rowSums(cover0[filt, ])
  totcov <- matrix(totcov, ncol = 3, byrow = TRUE)

  # for each species...
  for(sp in spp) {
    # ticker, if desired
    #print(paste0(site, ': ', sp, ', sp.', match(sp, spp), ' of ', length(spp))) ; flush.console()
    
    # relative abundances in TT1s by year
    N <- data.frame(sp = cover0[filt, sp] , yr = cover.meta$Year[filt], row.names = NULL)
    cover1 <- mutate(N, sp = sp / rowSums(cover0[filt, ]))
    N <- matrix(N$sp, ncol = 3, byrow = TRUE)
    cover1 <- matrix(cover1$sp, ncol = 3, byrow = TRUE)
    
    # Using mean relative abundances in TTCs for local flora source
    N.flora <- data.frame(sp = cover0[filt.flora, sp] , yr = cover.meta$Year[filt.flora], row.names = NULL)
    cover1.flora <- mutate(N.flora, sp = sp / rowSums(cover0[filt.flora, ]))
    N.flora <- matrix(N.flora$sp, ncol = 3, byrow = TRUE)
    cover1.flora <- matrix(cover1.flora$sp, ncol = 3, byrow = TRUE)
    regional <- colMeans(cover1.flora)

    # determine N
    rowz <- nrow(N)

    # organize into a list
    data <- list("cover1", "totcov", "regional", "N", "rowz")

    # Initialize
    inits <- seq(0.1, 0.9, 0.5)
    inits <- lapply(as.list(inits), function(x) list(m = x))
    reps <- length(inits)
    parameters <- c("m")

    if(!0 %in% colSums(N.flora)) {
      sink(paste0(wd, "\\model.txt"))
      cat("
        model { 
          for(i in 1:rowz) { #i is plot
            for(t in 1:2) {  #t is year         
              N[i, t + 1] ~ dpois(lambda[i, t])
              lambda[i, t] <- totcov[i, t] * ((1 - m) * cover1[i, t] + m * regional[t])
            }
         }  
          m ~ dunif(0,1)
        }
        ",fill=TRUE)
      sink()

      m1 <- jags(data, inits, parameters, "model.txt", n.thin = 50, n.chains = reps, n.burnin = 1000, n.iter = 10000)
      rhat <- rbind(rhat, m1$BUGSoutput$summary[, 'Rhat'])
      m1.mcmc <- as.mcmc(m1)

      for (j in 0:(reps - 1)) {
        dats[[i + j]] <- as.data.frame(m1.mcmc[[j + 1]])$m
        dats.meta <- rbind(dats.meta, data.frame(site = site, 
                                                 sp   = sp, 
                                                 rep  = j + 1))
      }
    }
    
    i <- i + reps

  }
}

lengs <- lapply(dats, length)
dats <- unlist(lapply(dats[lengs > 0], mean))
m.bayes <- mutate(dats.meta, m = dats)

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

# m.bayes <- group_by(m.bayes, site) %>%
#            summarise(m = mean(m))

# write.csv(m.bayes, file = "data\\m.bayes.csv", row.names = FALSE)