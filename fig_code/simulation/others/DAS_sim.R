set.seed(12345)
library(phyloseq)
library(doSNOW)
data("GlobalPatterns")

#The simulation code is downloaded from URL https://doi.org/10.1371/journal.pcbi.1003531.s001
#see the supplementary for detailed parameter Settings, and the main part of the code is the same as the original

# Define the number of true-positives per sample (10% in 100 taxon)
nTP = 10

# Minimum number of reads to consider an OTU 'observed' in a sample (about 5% of the sample size)
minobs = 1L

# The delimiter in the command parameter string
comdelim = "_"

# Define the different biological source templates to use sampletypes =
# c('Ocean', 'Soil')
sampletypes = levels(get_variable(GlobalPatterns, "SampleType"))[1]  #Feces

# Define the ceiling in the number of OTUs to consider in the template
nOTUs = 100L

# Define the number of samples in each class of a simulated experiment J =
J = 15#(n1=n2=25)

# The different values of effect size to apply foldeffect = c(1.5, 5, 20)
foldeffect = 5

## The number of reads per sample, ns ns = c(2000, 5E4)
ns = 10000

# Vector of the replicate numbers to repeat for each comb of simulation
# parameters (n, etc) reps=1:2
reps = 1:100

# The number of cores to use in this simulation Ncores = 7
Ncores = 50

# Define the simulation parameters combinations
simparams = apply(expand.grid(ns, sampletypes, reps, foldeffect, J), 1, paste0, 
                  collapse = comdelim)
# Define the labels to go with each element of the simulation parameter
# after splitting on the delimiter
simparamslabels = c("nreads", "SampleType", "Replicate", "EffectSize", "nsamples")

sampsums = sample_sums(GlobalPatterns)
keepsamples = sample_data(GlobalPatterns)$SampleType %in% sampletypes
template = prune_samples(keepsamples, GlobalPatterns)
# Make a list of source templates
templatelist = lapply(sampletypes, function(i, tempall, minobs, nOTUs) {
  cat(i, "\n")
  whtemp = (get_variable(tempall, "SampleType") %in% i)
  templatei = prune_samples(whtemp, tempall)
  samobs = apply(otu_table(templatei), 1, function(x, m) sum(x > m), m = minobs)
  otudf = data.frame(prev = samobs, sums = taxa_sums(templatei))
  otudf = otudf[order(-otudf$prev, -otudf$sums), ]
  # Trim all but the first nOTUs
  return(prune_taxa(rownames(otudf)[1:nOTUs], templatei))
}, template, minobs, nOTUs)
names(templatelist) <- sampletypes

microbesim = function(postfix = "sim", template, J, n = 10000) {
  # Generate `J` simulated microbiomes with `n` total reads each (all the
  # same, or n has length equal to the value of `J`), with subsamples drawn
  # from `template`.  `postfix` is a dummy idenitifer added to help
  # distinguish simulated samples in downstream code.
  require("phyloseq")
  # call the proporitions vector `pi`, similar to nomenclature from DMN
  pi = taxa_sums(template)
  # n must be a scalar (recycled as the number of reads for every simulation)
  # or it can be vector of length equal to J, the number of samples being
  # simulated.
  if (length(J) != 1) {
    stop("Length of J should be 1.")
  }
  if (length(n) != 1 & length(n) != J) {
    stop("n should be length 1, or length J.")
  }
  # Actually create the simulated abundance table
  simat = mapply(function(i, x, sample.size) {
    if (FALSE) 
    {
      print(i)
    }  # i is a dummy iterator
    phyloseq:::rarefaction_subsample(x, sample.size)
  }, i = 1:J, sample.size = n, MoreArgs = list(x = pi), SIMPLIFY = TRUE)
  simat = t(simat)
  # Add the OTU names to the OTU (column) indices
  colnames(simat) <- names(pi)
  # Add new simulated sample_names to the row (sample) indices
  rownames(simat) <- paste(i, "::", 1:nrow(simat), postfix, sep = "")
  # Put simulated abundances together with metadata as a phyloseq object
  OTU = otu_table(simat, taxa_are_rows = FALSE)
  # Define data.frame that will become sample_data
  SDF = data.frame(sample = sample_names(OTU), TableNumber = i, type = "simulated")
  SDF$postfix <- postfix
  rownames(SDF) <- sample_names(OTU)
  SD = sample_data(SDF)
  # Return a phyloseq object
  return(phyloseq(OTU, SD))
}

sumsim = function(n, sumtemplate, J) {
  # `n` - expected size target `sumtemplate` - the template vector of library
  # sizes observed in template `J` - The number of sample sizes to return
  scaledSums = round(n * (sumtemplate/median(sumtemplate)))
  return(sample(scaledSums, size = J, replace = TRUE))
}
cl <- makeCluster(Ncores, type = "SOCK") 
registerDoSNOW(cl)
simlist <- foreach(i = simparams, .packages = c("phyloseq")) %dopar% {
  # i = simparams[1]
  params = strsplit(i, comdelim)[[1]]
  names(params) <- simparamslabels
  # Initialize
  n = sim = sim1 = sim2 = n1 = n2 = NULL
  # cat(i, '\n')
  n = as.numeric(params["nreads"])
  sampletypei = params["SampleType"]
  # The number of samples to use for each class in this simulation
  Ji = as.integer(params["nsamples"])
  templatei = templatelist[[sampletypei]]
  # Rarely a simulation has a weird value and fails.  Catch these with `try`,
  # and repeat the simulation call if error (it will be a new seed)
  tryAgain = TRUE
  infiniteloopcounter = 1
  while (tryAgain & infiniteloopcounter < 5) {
    n1 = sumsim(n, sampsums, Ji)
    n2 = sumsim(n, sampsums, Ji)
    sim1 = microbesim(paste0(sampletypei, ";grp1"), templatei, Ji, n1)
    sim2 = microbesim(paste0(sampletypei, ";grp2"), templatei, Ji, n2)
    if (is.null(sim1) | is.null(sim2) | is.null(n1) | is.null(n2) | inherits(sim1, 
                                                                             "try-error") | inherits(sim2, "try-error")) {
      tryAgain = TRUE
      infiniteloopcounter = infiniteloopcounter + 1
    } else {
      tryAgain = FALSE
    }
  }
  if (infiniteloopcounter >= 5) {
    stop("Consistent error found during simulation. Need to investigate cause.")
  }
  # Merge the two simulated datasets together into one phyloseq object and add
  # back tree.
  sim = merge_phyloseq(sim1, sim2)
  sim = merge_phyloseq(sim, tax_table(GlobalPatterns), phy_tree(GlobalPatterns))
  return(sim)
}
names(simlist) <- simparams

simpletrim = function(physeq, minobs) {
  Ji = nsamples(physeq)
  # Force orientation to be sample-by-OTU
  if (taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  # `prevalence` is the fraction of total samples in which an OTU is observed
  # at least `minobs` times.
  prevalence = apply(as(otu_table(physeq), "matrix"), 2, function(x, minobs) {
    return(sum(x > minobs))
  }, minobs)/(Ji)
  # Will only keep OTUs that appear in more than X% of samples and have total
  # reads greater than half the number of samples.
  keepOTUs = prevalence > 0.05 & taxa_sums(physeq) > (0.5 * Ji)
  return(prune_taxa(keepOTUs, physeq))
}
# Simple prune
simlist0 <- lapply(simlist, simpletrim, minobs)


# Archive the simlist before adding an effect
simlistnoeffect = simlist0
# Apply the specified effect-size to randomly chosen OTUs
simlist0 <- foreach(i = simparams, .packages = c("phyloseq")) %dopar% {
  # physeq = simlist0[[3]]
  physeq = simlist0[[i]]
  params = strsplit(i, comdelim)[[1]]
  names(params) <- simparamslabels
  effectsize = as.numeric(params["EffectSize"])
  # Randomly sample from the available OTU names in physeq and assign them as
  # TP
  TPOTUs = sample(taxa_names(physeq), nTP, replace = FALSE)
  # Define the samples that will have effect applied
  effectsamples = grep(";grp1", sample_names(physeq), fixed = TRUE)
  # Apply effect (multiply abundances by effectsize scalar)
  otu_table(physeq)[effectsamples, TPOTUs] <- effectsize * otu_table(physeq)[effectsamples, 
                                                                             TPOTUs]
  # Rename these new 'true positive' OTUs, with 'TP'
  wh.TP = taxa_names(physeq) %in% TPOTUs
  newname = paste0(taxa_names(physeq)[wh.TP], "-TP")
  colnames(physeq@otu_table)[wh.TP] <- newname
  physeq@phy_tree$tip.label[wh.TP] <- newname
  rownames(physeq@tax_table)[wh.TP] <- newname
  return(physeq)
}
names(simlist0) <- names(simlist)

rngseeds = sample(20130921, size = length(simlist0), replace = TRUE)
rarelist = foreach(rngseedi = rngseeds, physeq = simlist0, .packages = "phyloseq") %dopar% 
  {
    # physeq = simlist0[[1]] rngseedi = rngseeds[[1]]
    physeq = rarefy_even_depth(physeq, rngseed = rngseedi)
    return(physeq)
  }
names(rarelist) <- names(simlist)
stopCluster(cl)
saveRDS(rarelist, file = "DATA_DAS_sig0.1.RDS")











set.seed(12345)
#The simulation code is downloaded from URL https://doi.org/10.1371/journal.pcbi.1003531.s001
#see the supplementary for detailed parameter Settings, and the main part of the code is the same as the original

# Define the number of true-positives per sample (20% in 100 taxon)
nTP = 20

# Minimum number of reads to consider an OTU 'observed' in a sample
minobs = 1L

# The delimiter in the command parameter string
comdelim = "_"

# Define the different biological source templates to use sampletypes =
# c('Ocean', 'Soil')
sampletypes = levels(get_variable(GlobalPatterns, "SampleType"))[1]  #Feces

# Define the ceiling in the number of OTUs to consider in the template
nOTUs = 100L

# Define the number of samples in each class of a simulated experiment J =
J = 15#(n1=n2=25)

# The different values of effect size to apply foldeffect = c(1.5, 5, 20)
foldeffect = 5

## The number of reads per sample, ns ns = c(2000, 5E4)
ns = 10000

# Vector of the replicate numbers to repeat for each comb of simulation
# parameters (n, etc) reps=1:2
reps = 1:100

# The number of cores to use in this simulation Ncores = 7
Ncores = 50

# Define the simulation parameters combinations
simparams = apply(expand.grid(ns, sampletypes, reps, foldeffect, J), 1, paste0, 
                  collapse = comdelim)
# Define the labels to go with each element of the simulation parameter
# after splitting on the delimiter
simparamslabels = c("nreads", "SampleType", "Replicate", "EffectSize", "nsamples")

sampsums = sample_sums(GlobalPatterns)
keepsamples = sample_data(GlobalPatterns)$SampleType %in% sampletypes
template = prune_samples(keepsamples, GlobalPatterns)
# Make a list of source templates
templatelist = lapply(sampletypes, function(i, tempall, minobs, nOTUs) {
  cat(i, "\n")
  whtemp = (get_variable(tempall, "SampleType") %in% i)
  templatei = prune_samples(whtemp, tempall)
  samobs = apply(otu_table(templatei), 1, function(x, m) sum(x > m), m = minobs)
  otudf = data.frame(prev = samobs, sums = taxa_sums(templatei))
  otudf = otudf[order(-otudf$prev, -otudf$sums), ]
  # Trim all but the first nOTUs
  return(prune_taxa(rownames(otudf)[1:nOTUs], templatei))
}, template, minobs, nOTUs)
names(templatelist) <- sampletypes

microbesim = function(postfix = "sim", template, J, n = 10000) {
  # Generate `J` simulated microbiomes with `n` total reads each (all the
  # same, or n has length equal to the value of `J`), with subsamples drawn
  # from `template`.  `postfix` is a dummy idenitifer added to help
  # distinguish simulated samples in downstream code.
  require("phyloseq")
  # call the proporitions vector `pi`, similar to nomenclature from DMN
  pi = taxa_sums(template)
  # n must be a scalar (recycled as the number of reads for every simulation)
  # or it can be vector of length equal to J, the number of samples being
  # simulated.
  if (length(J) != 1) {
    stop("Length of J should be 1.")
  }
  if (length(n) != 1 & length(n) != J) {
    stop("n should be length 1, or length J.")
  }
  # Actually create the simulated abundance table
  simat = mapply(function(i, x, sample.size) {
    if (FALSE) 
    {
      print(i)
    }  # i is a dummy iterator
    phyloseq:::rarefaction_subsample(x, sample.size)
  }, i = 1:J, sample.size = n, MoreArgs = list(x = pi), SIMPLIFY = TRUE)
  simat = t(simat)
  # Add the OTU names to the OTU (column) indices
  colnames(simat) <- names(pi)
  # Add new simulated sample_names to the row (sample) indices
  rownames(simat) <- paste(i, "::", 1:nrow(simat), postfix, sep = "")
  # Put simulated abundances together with metadata as a phyloseq object
  OTU = otu_table(simat, taxa_are_rows = FALSE)
  # Define data.frame that will become sample_data
  SDF = data.frame(sample = sample_names(OTU), TableNumber = i, type = "simulated")
  SDF$postfix <- postfix
  rownames(SDF) <- sample_names(OTU)
  SD = sample_data(SDF)
  # Return a phyloseq object
  return(phyloseq(OTU, SD))
}

sumsim = function(n, sumtemplate, J) {
  # `n` - expected size target `sumtemplate` - the template vector of library
  # sizes observed in template `J` - The number of sample sizes to return
  scaledSums = round(n * (sumtemplate/median(sumtemplate)))
  return(sample(scaledSums, size = J, replace = TRUE))
}
cl <- makeCluster(Ncores, type = "SOCK") 
registerDoSNOW(cl)
simlist <- foreach(i = simparams, .packages = c("phyloseq")) %dopar% {
  # i = simparams[1]
  params = strsplit(i, comdelim)[[1]]
  names(params) <- simparamslabels
  # Initialize
  n = sim = sim1 = sim2 = n1 = n2 = NULL
  # cat(i, '\n')
  n = as.numeric(params["nreads"])
  sampletypei = params["SampleType"]
  # The number of samples to use for each class in this simulation
  Ji = as.integer(params["nsamples"])
  templatei = templatelist[[sampletypei]]
  # Rarely a simulation has a weird value and fails.  Catch these with `try`,
  # and repeat the simulation call if error (it will be a new seed)
  tryAgain = TRUE
  infiniteloopcounter = 1
  while (tryAgain & infiniteloopcounter < 5) {
    n1 = sumsim(n, sampsums, Ji)
    n2 = sumsim(n, sampsums, Ji)
    sim1 = microbesim(paste0(sampletypei, ";grp1"), templatei, Ji, n1)
    sim2 = microbesim(paste0(sampletypei, ";grp2"), templatei, Ji, n2)
    if (is.null(sim1) | is.null(sim2) | is.null(n1) | is.null(n2) | inherits(sim1, 
                                                                             "try-error") | inherits(sim2, "try-error")) {
      tryAgain = TRUE
      infiniteloopcounter = infiniteloopcounter + 1
    } else {
      tryAgain = FALSE
    }
  }
  if (infiniteloopcounter >= 5) {
    stop("Consistent error found during simulation. Need to investigate cause.")
  }
  # Merge the two simulated datasets together into one phyloseq object and add
  # back tree.
  sim = merge_phyloseq(sim1, sim2)
  sim = merge_phyloseq(sim, tax_table(GlobalPatterns), phy_tree(GlobalPatterns))
  return(sim)
}
names(simlist) <- simparams

simpletrim = function(physeq, minobs) {
  Ji = nsamples(physeq)
  # Force orientation to be sample-by-OTU
  if (taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  # `prevalence` is the fraction of total samples in which an OTU is observed
  # at least `minobs` times.
  prevalence = apply(as(otu_table(physeq), "matrix"), 2, function(x, minobs) {
    return(sum(x > minobs))
  }, minobs)/(Ji)
  # Will only keep OTUs that appear in more than X% of samples and have total
  # reads greater than half the number of samples.
  keepOTUs = prevalence > 0.05 & taxa_sums(physeq) > (0.5 * Ji)
  return(prune_taxa(keepOTUs, physeq))
}
# Simple prune
simlist0 <- lapply(simlist, simpletrim, minobs)


# Archive the simlist before adding an effect
simlistnoeffect = simlist0
# Apply the specified effect-size to randomly chosen OTUs
simlist0 <- foreach(i = simparams, .packages = c("phyloseq")) %dopar% {
  # physeq = simlist0[[3]]
  physeq = simlist0[[i]]
  params = strsplit(i, comdelim)[[1]]
  names(params) <- simparamslabels
  effectsize = as.numeric(params["EffectSize"])
  # Randomly sample from the available OTU names in physeq and assign them as
  # TP
  TPOTUs = sample(taxa_names(physeq), nTP, replace = FALSE)
  # Define the samples that will have effect applied
  effectsamples = grep(";grp1", sample_names(physeq), fixed = TRUE)
  # Apply effect (multiply abundances by effectsize scalar)
  otu_table(physeq)[effectsamples, TPOTUs] <- effectsize * otu_table(physeq)[effectsamples, 
                                                                             TPOTUs]
  # Rename these new 'true positive' OTUs, with 'TP'
  wh.TP = taxa_names(physeq) %in% TPOTUs
  newname = paste0(taxa_names(physeq)[wh.TP], "-TP")
  colnames(physeq@otu_table)[wh.TP] <- newname
  physeq@phy_tree$tip.label[wh.TP] <- newname
  rownames(physeq@tax_table)[wh.TP] <- newname
  return(physeq)
}
names(simlist0) <- names(simlist)

rngseeds = sample(20130921, size = length(simlist0), replace = TRUE)
rarelist = foreach(rngseedi = rngseeds, physeq = simlist0, .packages = "phyloseq") %dopar% 
  {
    # physeq = simlist0[[1]] rngseedi = rngseeds[[1]]
    physeq = rarefy_even_depth(physeq, rngseed = rngseedi)
    return(physeq)
  }
names(rarelist) <- names(simlist)
stopCluster(cl)
saveRDS(rarelist, file = "DATA_DAS_sig0.2.RDS")

