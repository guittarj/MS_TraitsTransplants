loadpax(c("scales","R2jags", "lattice"))
cover0 <- round((cover * 100) / rowSums(cover))
dats <- list()
dats.meta <- data.frame(site = character(), sp = character(), rep = integer())

i <- 1
rhat <- data.frame(deviance = numeric(), m = numeric())

for(site in unique(cover.meta$siteID)) {
  
  # So, throughout here, i use 'flora' to refer to the site flora. 
  filt.flora <- cover.meta$TTtreat == 'TTC' & 
                cover.meta$destSiteID == site &
                cover.meta$Year %in% c(2011:2013)
  filt <- cover.meta$TTtreat == 'TT1' & 
          cover.meta$destSiteID == site &
          cover.meta$Year %in% c(2011:2013)
  
  # list of site species
  spp <- sort(colSums(cover0[filt, ]), decreasing = TRUE)
  spp <- names(spp[spp > 0])
  
  # matrix of total cover for each turf * year
  totcov <- rowSums(cover0[filt, ])
  totcov <- matrix(totcov, ncol = 3, byrow = TRUE)

  for(sp in spp) {
    print(paste0(site, ': ', sp, ', sp.', match(sp, spp), ' of ', length(spp))) ; flush.console()
    
    # N is the number of individuals of sp. X in each plot*year (row)
    N <- data.frame(sp = cover0[filt, sp] , yr = cover.meta$Year[filt], row.names = NULL)
    
    # cover1 is the percent cover of sp X in each plot*year (row)
    cover1 <- mutate(N, sp = sp / rowSums(cover0[filt, ]))
    
    # turn them into 3 columns ~ year matrices
    N <- matrix(N$sp, ncol = 3, byrow = TRUE)
    cover1 <- matrix(cover1$sp, ncol = 3, byrow = TRUE)
    
    # Now do the same for the site floras
    N.flora <- data.frame(sp = cover0[filt.flora, sp] , yr = cover.meta$Year[filt.flora], row.names = NULL)
    cover1.flora <- mutate(N.flora, sp = sp / rowSums(cover0[filt.flora, ]))
    N.flora <- matrix(N.flora$sp, ncol = 3, byrow = TRUE)
    cover1.flora <- matrix(cover1.flora$sp, ncol = 3, byrow = TRUE)
    
    # site means
    regional <- colMeans(cover1.flora)
    
    # number of rows
    rowz <- nrow(N)

    data <- list("cover1", "totcov", "regional", "N", "rowz")

    # initial values for the MCMC
    inits <- lapply(as.list(c(0.1, 0.9, 0.5)), function(x) list(m = x))
    reps <- length(inits)
  
    parameters <- c("m")
    
    # note that I have to filter out situations where the species does not appear
    # in the site flora
    if(!0 %in% colSums(N.flora)) {
      sink(paste0(wd, "\\model.txt"))
      cat("
        model { 
          for(i in 1:rowz) { #i is plot
            for(t in 1:2) {  #t is year        

# Here, we are trying to predict N[, t+1] as a function of total.turf.cover[t] * 
# (prob of intra-turf replacement * sp X turf percent.cover[t] + 
#  prob of immigration * sp X * site percent.cover[t])
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