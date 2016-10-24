# My custom/collected functions

loadpax <- function(pkg){
  # (1) checks package installation, (2) installs them if not, then (3) loads them
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}

cwm <- function (cover.row, trait) {
  # Calculate Community Weighted Mean
  # row:   row of cover data
  # trait: trait of interest (must be a column name in trait.data)

  x <- trait.data[, trait]
  weighted.mean(x[!is.na(x)], w = cover.row[!is.na(x)])
}

dist.fun <- function(cover, trait) {
  # calculate distance matrix for 2+ communities
  # cover: a dataframe with communities (rows) and species abundances (cols)
  # trait: trait of interest, determines how distance is calculated

  if(trait == 'veg') {                                
    dmat <- as.matrix(vegdist(cover))                   # Bray-Curtis distance
  } else {                                             
    dmat <- apply(cover, 1, function(x) cwm(x, trait))  # Euclidean distance
    dmat <- as.matrix(dist(dmat))
  }
  return(dmat)
}

comp.fun <- function (dmat, dmat.meta, trait) {
  # determine mean distance between a turf and its local controls
  # dmat: distance matrix
  # dmat.meta: metadata for matrix
  # trait: trait used to calculate distance matrix
  
  comp <- transmute(dmat.meta, trait  = trait,
                               turfID = turfID, 
                               year   = Year)
  for (i in 1:nrow(dmat)) {
    controls <- with(dmat.meta, 
                  turf.year[destSiteID == destSiteID[i] & 
                            Year       == Year[i] &
                            TTtreat  %in% c('TT1', 'TTC') & 
                            turfID     != turfID[i]])
    comp$dist.tt1[i] <- mean(dmat[i, controls])
  }
  return(comp)
}

stat_sum_df <- function(fun, geom = 'crossbar', ...) {
  # For plotting
  stat_summary(fun.data = fun, geom = geom, width = 0.2, ...)
}

stat_sum_single <- function(fun, geom = 'point', ...) {
  # For plotting
  stat_summary(fun.y = fun, geom = geom, ...)
}

fmt <- function() {
  # For plotting
  function(x) as.character(round(x, 2))
}

process.comm = function(comm, trait) {
  # Takes a simulation file, and calculates distances to field controls

  comp.sims <- rep(NA, nrow(comm))

  for(i in 1:nrow(comm)) {

    # combine simturf and controls
    simturf <- paste(comm[i, 'turfID'], comm[i, 'year'], sep = '_')
    controls <- with(cover.meta, as.character(turf.year)[TTtreat  %in% c('TTC','TT1') &
                     destSiteID == destSiteID[turf.year == simturf] &
                     Year       == comm[i, 'year'] &
                     turf.year  != simturf])
    
    # Discard unnecessasry rows (won't affect distance measures)
    cover.sim <- cover[c(simturf, controls), ]
    cover.sim[simturf, ] <- as.numeric(comm[i, colnames(cover)])

    # Calculate distance
    comp.sim <- dist.fun(cover.sim, trait)
    comp.sims[i] <- mean(comp.sim[simturf, controls])
  }

  comp.sims <- data.frame(trait = trait, comm[, c('turfID', 'year', 'm', 'd')], dist.tt1 = comp.sims, reps = 1)
  return(comp.sims)
}

probfixes = function (vec, probs, fixes) {
  # A function that resolves any naming inconsitencies in a vector (vec)

  for(i in 1:length(probs)) {
    vec <- ifelse(as.character(vec) == as.character(probs[i]), as.character(fixes[i]), as.character(vec))
  }
  return(vec)
}

rao.fun <- function(traitvals, abuns) {
  abuns <- abuns[!is.na(traitvals)]
  traitvals <- traitvals[!is.na(traitvals)]
  p <- abuns / sum(abuns)
  d <- as.matrix(dist(traitvals))
  rao <- as.numeric(p %*% (d %*% p) / 2)
  return(rao)
}

grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  # A function that lets you put multiple plots together and have a single shared legend (from the first plot)
  # from hadley on the internet...
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)

  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)

}

# print list of loaded functions
print(data.frame(Custom_Functions = 
  c('loadpax: install+load multiple packages',
    'cwm: calculate community weighted means',
    'dist.fun: get distances among communities',
    'comp.fun: get distances to local control',
    'stat_sum_df: a plotting function', 
    'stat_sum_single: a plotting function', 
    'fmt: a plotting function',
    'process.comm: Calculates trait/veg distances to controls',
    'probfixes: corrects taxonomic inconsitencies',
    'grid_arrange_shared_legend: Multiple plots, one legend',
    'rao-fun: calculates Rao quadratic entropy index')))
