library(gdata)
library(DOQTL)
library(GenomicRanges)
library(VariantAnnotation)
library(foreach)
library(doParallel)
library(RColorBrewer)
#library(abind)


## Help Documentation
describe = function(obj) {
  if ('help' %in% names(attributes(obj))) {
    writeLines(attr(obj, 'help'))
  }
}
attr(describe, 'help') = "
This function prints the contents of the 'help' attribute of any R object. 
It is meant to provide help documentation in the same vein as Docstrings in Python. 
"

geno_states = c("AA", "BB", "CC", "DD", "EE", "FF", "GG", "HH", "AB", "AC", 
"AD", "AE", "AF", "AG", "AH", "BC", "BD", "BE", "BF", "BG", "BH", 
"CD", "CE", "CF", "CG", "CH", "DE", "DF", "DG", "DH", "EF", "EG", 
"EH", "FG", "FH", "GH")


collapse_probs = function(probs_file, states=geno_states, gbuild='b38') {
  ## Load probabilities
  full_probs = read.csv(probs_file, as.is=T)
  
  ## Founder states
  founder_states = c('AA', 'BB', 'CC', 'DD', 'EE', 'FF', 'GG', 'HH')
  
  ## Get heterozygous states for each founder
  het_states = list()
  for (fs in founder_states) {
	letter = strsplit(fs, '')[[1]][1]
	pattern = paste(letter,'[^',letter,']|[^',letter,']',letter,sep='')
    het_states[[fs]] = states[grep(pattern, states)]
  }
  
  ## Initialize the collapsed dataframe with only the homozygous states (AA, BB, etc.)
  coll_probs = full_probs[, 1:11]
  
  ## Collapsed probabilities for each founder is the homozygous probability plus 0.5*probability of all the heterozygous states
  for (fs in founder_states) {
	coll_probs[, fs] = full_probs[, fs] + rowSums(full_probs[, het_states[[fs]]])/2
  }
  
  colnames(coll_probs) = c('marker', 'chromosome', paste0('position_', gbuild), 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H')
  rownames(coll_probs) = coll_probs$marker
  
  return(coll_probs)
}
attr(collapse_probs, 'help') = "
This function collapses 36-state probabilities to 8 probabilites (one for each founder strain).

Parameters:
probs_file: The filename of a .csv file containing the 36-state probabilities. Reading this file 
   with read.csv() should yield an n x 39 dataframe containing marker info (marker, chromosome, position) 
   for n markers, and the 36-state probabilities.
states: A character vector containing the 36 genotype states (default is the geno_states global variable
   as defined in this script).

Returns:
A n x 11 dataframe, where n is the number of markers and the columns represent the
marker, chromosome, position, and probabilities for each of the 8 founder strains. 
"


make_rix_probs = function(dam_file, sire_file) {
  ## Load probabilities
  dam = read.csv(dam_file, as.is=T)
  sire = read.csv(sire_file, as.is=T)
  
  ## Standardize column names, just in case
  colnames(dam)[1:2] = c('marker', 'chromosome')
  colnames(sire)[1:2] = c('marker', 'chromosome')
  
  ## Check that markers from parental lines match
  if (!identical(dam$marker, sire$marker)) {
    stop("Markers for dam and sire do not match!")
  }
  
  ## Check that genotype states from parental lines match
  dam_states = colnames(dam)[4:length(colnames(dam))]
  sire_states = colnames(sire)[4:length(colnames(sire))]
  if (!identical(dam_states, sire_states)) {
	stop("Genotype states for dam and sire do not match!")
  } else {
	states = dam_states
  }
  
  rownames(dam) = dam$marker
  rownames(sire) = sire$marker
  autosomes = sort(unique(dam$chromosome[!dam$chromosome %in% c('M', 'X', 'Y')]))
  
  ## Create 3D array
  probs = abind(dam[,states], sire[,states], along=3)
  dimnames(probs)[[3]] = c('dam', 'sire')
  names(dimnames(probs)) = c('markers', 'states', 'samples')
  
  ## Calculate marker probabilities for males
  ## Initialize probabilities as 0
  male_df = sire
  male_df[,states] = 0
  
  ## Autosomes
  ## Average the probabilities from dam and sire
  autosome_markers = male_df$marker[male_df$chromosome %in% autosomes]
  if (length(autosome_markers) > 0) {
    male_df[autosome_markers, states] = apply(probs[autosome_markers, states, ], c('markers', 'states'), function(x) {.Internal(mean(x))})
  }
  
  ## Y chromosome
  ## Copy the probabilities from sire
  y_markers = male_df$marker[male_df$chromosome == 'Y']
  if (length(y_markers) > 0) {
    male_df[y_markers, states] = probs[y_markers, states, 'sire']
  }
  
  ## X chromosome
  ## Copy the probabilities from dam
  x_markers = male_df$marker[male_df$chromosome == 'X']
  if (length(x_markers) > 0) {
    male_df[x_markers, states] = probs[x_markers, states, 'dam']
  }
  
  ## M chromosome
  ## Copy the probabilities from dam
  m_markers = male_df$marker[male_df$chromosome == 'M']
  if (length(m_markers) > 0) {
    male_df[m_markers, states] = probs[m_markers, states, 'dam']
  }

  ## Calculate marker probabilities for females
  ## Initialize probabilities as 0
  female_df = dam
  female_df[,states] = 0
  
  ## Autosomes and X chromosome
  ## Average the probabilities from dam and sire
  if (length(c(autosome_markers, x_markers)) > 0) {
    female_df[c(autosome_markers, x_markers), states] = apply(probs[c(autosome_markers, x_markers), states, ], c('markers', 'states'), function(x) {.Internal(mean(x))})
  }

  ## M chromosome
  ## Copy the probabilities from dam
  if (length(m_markers) > 0) {
    female_df[m_markers, states] = probs[m_markers, states, 'dam']
  }
  
  ## Y chromosome
  ## All probabilities = 0
  
  return(list(male=male_df, female=female_df))
}
attr(make_rix_probs, 'help') = "
This function calculates the genotype state probabilities for RIX lines using the probabilities 
from the two parental RI lines.

Parameters:
dam: An n x 3+s dataframe containing marker info (marker, chromosome, position) for n markers,
   and s genotype state probabilities.
sire: An n x 3+s dataframe containing marker info (marker, chromosome, position) for n markers,
   and s genotype state probabilities. 

Returns:
A list with two (male and female) n x 3+s dataframes containing the marker and probability 
information for a RIX (F1-hybrid of the given dam and sire).
"


make_rix_model_probs = function(samples, file_path='.', sex) {
  ## Founder states
  founder_states = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H')
  
  if (missing(sex)) {
    stop("Sex ('M' or 'F') must be specified!")
  }
  
  ## Get probabilities for the first sample
  ## Initialize the 3D array
  mating = unlist(strsplit(samples[1], '_'))[1]
  dam = unlist(strsplit(mating, 'x'))[1]
  sire = unlist(strsplit(mating, 'x'))[2] 
  dam_file = file.path(file_path, paste0(dam, '.csv'))
  sire_file = file.path(file_path, paste0(sire, '.csv'))
  
  if (sex == 'M') {
    rix_probs = make_rix_probs(dam_file, sire_file)$male[, founder_states, ]
  } else {
    rix_probs = make_rix_probs(dam_file, sire_file)$female[, founder_states, ]
  }
  model_prob_array = rix_probs
  
  mating_count = 1
  prev_mating = mating
  ## Iterate through the rest of the samples
  for (s in samples[2:length(samples)]) {
    mating = unlist(strsplit(s, '_'))[1]
    
	## For each mating
	if (mating != prev_mating) {
	  print(mating_count)
	  mating_count = mating_count + 1
      dam = unlist(strsplit(mating, 'x'))[1]
      sire = unlist(strsplit(mating, 'x'))[2]
      dam_file = file.path(file_path, paste0(dam, '.csv'))
      sire_file = file.path(file_path, paste0(sire, '.csv'))
	
	  ## Create RIX probabilities
    if (sex == 'M') {
      rix_probs = make_rix_probs(dam_file, sire_file)$male[, founder_states, ]
    } else {
      rix_probs = make_rix_probs(dam_file, sire_file)$female[, founder_states, ]
    }
	}
    
    model_prob_array = abind(model_prob_array, rix_probs, along=3)
	
	prev_mating = mating
  }
  
  ## Assign names to array
  names(dimnames(model_prob_array)) = c('markers', 'founders', 'samples')
  dimnames(model_prob_array)[['samples']] = samples
  
  ## Reshape the 3D array
  model_prob_array = aperm(model_prob_array, c('samples', 'founders', 'markers'))
  
  print(mating_count)
  return(model_prob_array)
}
attr(make_rix_model_probs, 'help') = "
This function calculates the genotype state probabilities for RIX lines using the probabilities 
from the two parental RI lines.

Parameters:
samples: A character vector containing sample names (matings will work as well).
sex: 'M' or 'F'
file_path: The path to the directory containing the RIX probability .csv files.

Returns:
A Nx8x77606 array, where N is the number of samples (or matings).
"


get_min_marker = function(doqtl, chr=NULL, window=0) {
  if (is.null(chr)) {
    lod = doqtl$lod$A
  } else {
	if (chr == 'X') {
      lod = doqtl$lod$X
	} else {
      lod = doqtl$lod$A[doqtl$lod$A[,2]==chr, ]
	}
  }
  min_idx = which.min(lod$p)
  lod[(min_idx-window):(min_idx+window),]
}
attr(get_min_marker, 'help') = "
This function returns the marker with the minimum p-value, either across the 
autosomes (if chr==NULL), or within a specific chromosome. Adjacent markers can 
also be returned.

Parameters:
doqtl: A QTL object returned by scanone()
chr: A chromosome (1-19, or 'X')
window: An integer specifying how many markers on either side of the minimum marker
   should be returned.

Returns:
A dataframe containing the minimum marker (and adjacent markers).   
"


get_lod_drop = function(doqtl, chr=NULL, drop=1.5, marker=NULL) {
  if (is.null(chr)) {
    stop("You must specify a chromosome.")
  }
  
  if (chr == 'X') {
    if (!is.null(marker)) {
      lod = doqtl$lod$X
      max_idx = which(lod$marker==marker)
      max_lod = lod$lod[max_idx]
      max_lod_pos = lod[max_idx, 3]
      print(max_lod)
	} else {
	  lod = doqtl$lod$X
	  max_idx = which.max(lod$lod)
      max_lod = lod$lod[max_idx]
      max_lod_pos = lod[max_idx, 3]
      print(max_lod)
	}
  } else {
    if (!is.null(marker)) {
      lod = doqtl$lod$A
      max_idx = which(lod$marker==marker)
	  chr = lod[max_idx, 2]
	  lod = lod[lod[,2]==chr,]
	  max_idx = which(lod$marker==marker)
      max_lod = lod$lod[max_idx]
      max_lod_pos = lod[max_idx, 3]
      print(max_lod)
    } else {
      lod = doqtl$lod$A[doqtl$lod$A[,2]==chr, ]  
      max_idx = which.max(lod$lod)
      max_lod = lod$lod[max_idx]
      max_lod_pos = lod[max_idx, 3]
      print(max_lod)
	}
  }
  print(max_lod-drop)
  drop_start = max(lod[lod$lod <= (max_lod-drop) & lod[,3] < max_lod_pos, 3])
  drop_end = min(lod[lod$lod <= (max_lod-drop) & lod[,3] > max_lod_pos, 3])
  return(c(drop_start, drop_end))
}
attr(get_lod_drop, 'help') = "
This function calculates the lod drop for a given marker, or the minimum 
marker in the given chromosome.

Parameters:
doqtl:
chr:
drop:
marker:

Returns:
The start and stop coordinates for the specified LOD drop.
"


get_bootstrap_interval = function(doqtl, chr, start=NULL, end=NULL, B=100, p=0.95, pheno, pheno.col, probs, K, addcovar, snps, parallel=F) {
  peaks = c()
  
  if (is.null(chr)) {
    stop("You must specify a chromosome.")
  }
  
  ## Get markers
  if (chr == 'X') {
    lod = doqtl$lod$X
  } else {
    lod = doqtl$lod$A
  }
  
  ## If interval on chr is specified
  if (!is.null(start) & !is.null(end)) {
	markers = lod[lod[,2]==chr & lod[,3] >= start & lod[,3] <= end, 1]
  } else {
    markers = lod[lod[,2]==chr,1]
  }
  probs = probs[, , markers]
  snps = snps[snps[,1] %in% markers,]
  
  ## Subset pheno and addcovar to include only samples in probs and K
  pheno = pheno[rownames(K), , drop=F]
  addcovar = addcovar[rownames(K), , drop=F]
  
  if (!parallel) {
    for (i in 1:B) {
	  ## Print progress
	  if (i==1) {
	    cat('\n Bootstrap sample:', i)
		flush.console()
	  } else {
	    cat('\r','Bootstrap sample:', i)
		flush.console()
	  }
	  
      ## Select sub-sample
      n = dim(model.probs)[1]
      sample_idx = sample(c(1:n), n, replace=T)
      pheno_sub = pheno[sample_idx, , drop=F]
      addcovar_sub = addcovar[sample_idx, , drop=F]
	  sample_names = gsub('\\..*', '', rownames(pheno_sub))
	  probs_sub = probs[sample_names, , ]
	  K_sub = K[sample_names, sample_names]
	  
	  ## Set unique sample names
	  dimnames(probs_sub)[[1]] = rownames(pheno_sub)
	  rownames(K_sub) = rownames(pheno_sub)
	  colnames(K_sub) = rownames(K_sub)
	
	  ## Perform QTL scan
      msg = capture.output({qtl_sub = scanone(pheno=pheno_sub, pheno.col=pheno.col, probs=probs_sub, K=K_sub, addcovar=addcovar_sub, snps=snps)})
	  peaks = c(peaks, get_min_marker(qtl_sub, chr=chr)[1,3])
	  gc(verbose=F)
    }
  } else {
    if (is.numeric(parallel)) {
	  num_cores = parallel
	} else {
	  num_cores = detectCores() - 1
	}
	registerDoParallel(num_cores)
	cat('\nNumber of parallel processes:', num_cores)
	peaks = foreach(i=1:B, .combine=c) %dopar% {
	  gc(verbose=F)
	  
	  ## Print progress
	  if (i==1) {
	    cat('\n Bootstrap sample:', i)
		flush.console()
	  } else if (i %% num_cores == 0) {
	    cat('\r','Bootstrap sample:', i)
		flush.console()
	  }
	  
      ## Select sub-sample
      n = dim(model.probs)[1]
      sample_idx = sample(c(1:n), n, replace=T)
      pheno_sub = pheno[sample_idx, , drop=F]
      addcovar_sub = addcovar[sample_idx, , drop=F]
	  sample_names = gsub('\\..*', '', rownames(pheno_sub))
	  #print(sample_names[!sample_names %in% dimnames(probs)[[1]]])
	  probs_sub = probs[sample_names, , ]
	  K_sub = K[sample_names, sample_names]
	  
	  ## Set unique sample names
	  dimnames(probs_sub)[[1]] = rownames(pheno_sub)
	  rownames(K_sub) = rownames(pheno_sub)
	  colnames(K_sub) = rownames(K_sub)
	
	  ## Perform QTL scan
      msg = capture.output({qtl_sub = scanone(pheno=pheno_sub, pheno.col=pheno.col, probs=probs_sub, K=K_sub, addcovar=addcovar_sub, snps=snps)})
	  #print(get_min_marker(qtl_sub, chr=chr))
	  get_min_marker(qtl_sub, chr=chr)[1,3]
	}
	gc(verbose=F)
	stopImplicitCluster()
  }
  
  ## Return interval on specified chromosome
  q = 1 - p
  interval = quantile(peaks, c(q/2, (1-(q/2))))
  cat(' ... Done. \n\n')
  print(interval)
  return(peaks)
}
attr(get_bootstrap_interval, 'help') = "
This function calculates the bootstrap confidence interval for a QTL. 

Parameters:
doqtl: results from scanone()
chr: the chromosome containing the QTL of interest
start: an optional start coordinate (in Mb) on the chr (default=NULL)
end: an optional end coordinate (in Mb) on the chr (default=NULL)
B: number of bootstrap samples
p: the confidence interval to calculate (default 0.95)
pheno: the following 6 parameters are those passed to scanone()
pheno.col:
probs:
K:
addcovar:
snps:

Returns:
A vector containing the location of the max peaks for each bootstrap sample 
(the specified interval is also printed to the screen).
"

scanone.perm_v2 = function(pheno, pheno.col, probs, K, addcovar, snps, value='lod', B=100, parallel=F, ...) {
  
  ## Check that samples and matings match across pheno, K, and probs
  if (!identical(rownames(pheno), rownames(K)) | !identical(rownames(pheno), dimnames(probs)[[1]])) {
    stop("Sample names do not match!")
  }
  
  ## Get sample IDs
  sample_ids = rownames(pheno)
  
  ## Get matings and check that they match IDs
  if ('Mating' %in% colnames(pheno)) {
    sample_matings = pheno$Mating
	if (!identical(sample_matings, sapply(strsplit(sample_ids, "_"), function(x) {x[1]}))) {
	  stop("Matings in pheno do not match sample IDs!")
	}
  } else {
    sample_matings = sapply(strsplit(sample_ids, "_"), function(x) {x[1]})
  }
  
  unique_matings = unique(sample_matings)
  matings_table = table(sample_matings)[unique_matings]
  
  max_vals = c()
  
  ## Assign matings as names for K and probability array so they can be subset easily during permutations
  dimnames(probs)[[1]] = sample_matings
  rownames(K) = sample_matings
  colnames(K) = sample_matings
  
  if (!parallel) {
    for (i in 1:B) {
	  ## Print progress
	  if (i==1) {
	    cat('\n Permutation:', i)
		flush.console()
	  } else {
	    cat('\r','Permutation:', i)
		flush.console()
	  }
	  
      ## Permute population and create new IDs
	  new_unique_matings = sample(unique_matings)
	  new_matings = rep(new_unique_matings, matings_table)
	  new_matings_table = table(new_matings)[new_unique_matings]
	  new_sample_ids = unname(unlist(sapply(names(new_matings_table), function(x) {paste(x, 1:new_matings_table[x], sep='_')})))
	  
	  ## Permute probability array and assign new IDs
	  new_probs = probs[new_matings, , ]
	  dimnames(new_probs)[[1]] = new_sample_ids
	  
	  ## Permute K matrix and assign new IDs
	  new_K = K[new_matings, new_matings]
	  rownames(new_K) = new_sample_ids
	  colnames(new_K) = new_sample_ids
	  
	  ## Assign new sample IDs to pheno dataframe
	  new_pheno = pheno
	  rownames(new_pheno) = new_sample_ids
	  
	  ## Assign new sample IDs to covariate dataframe
	  new_addcovar = addcovar
	  rownames(new_addcovar) = new_sample_ids
	  
      ## Perform QTL scan
      msg = capture.output({qtl_perm = scanone(pheno=new_pheno, pheno.col=pheno.col, probs=new_probs, K=new_K, addcovar=new_addcovar, snps=snps)})
	  
	  if (value == "lod") {
	    max_A = max(qtl_perm$lod$A[,'lod'])
		if ('X' %in% names(qtl_perm$lod)) {
 		  max_X = max(qtl_perm$lod$X[,'lod'])
		} else {
		  max_X = NA
		}
		perm_max_val = max(max_A, max_X, na.rm=T)
	  } else {
	    max_A = min(qtl_perm$lod$A[,'p'])
		if ('X' %in% names(qtl_perm$lod)) {
		  max_X = min(qtl_perm$lod$X[,'p'])
		} else {
		  max_X = NA
		}
		perm_max_val = min(max_A, max_X, na.rm=T)
	  }
	  max_vals = c(max_vals, perm_max_val)
	  gc(verbose=F)
    }
  } else {
    if (is.numeric(parallel)) {
	  num_cores = parallel
	} else {
	  num_cores = detectCores() - 1
	}
	registerDoParallel(num_cores)
	cat('\nNumber of parallel processes:', num_cores)
	max_vals = foreach(i=1:B, .combine=c, ...) %dopar% {
	  gc(verbose=F)
	  
	  ## Print progress
	  if (i==1) {
	    cat('\n Permutation:', i)
		flush.console()
	  } else if (i %% num_cores == 0) {
	    cat('\r','Permutation:', i)
		flush.console()
	  }
	  
	  ## Permute population and create new IDs
	  new_unique_matings = sample(unique_matings)
	  new_matings = rep(new_unique_matings, matings_table)
	  new_matings_table = table(new_matings)[new_unique_matings]
	  new_sample_ids = unname(unlist(sapply(names(new_matings_table), function(x) {paste(x, 1:new_matings_table[x], sep='_')})))
	  
	  ## Permute probability array and assign new IDs
	  new_probs = probs[new_matings, , ]
	  dimnames(new_probs)[[1]] = new_sample_ids
	  
	  ## Permute K matrix and assign new IDs
	  new_K = K[new_matings, new_matings]
	  rownames(new_K) = new_sample_ids
	  colnames(new_K) = new_sample_ids
	  
	  ## Assign new sample IDs to pheno dataframe
	  new_pheno = pheno
	  rownames(new_pheno) = new_sample_ids
	  
	  ## Assign new sample IDs to covariate dataframe
	  new_addcovar = addcovar
	  rownames(new_addcovar) = new_sample_ids
	  
      ## Perform QTL scan
      msg = capture.output({qtl_perm = scanone(pheno=new_pheno, pheno.col=pheno.col, probs=new_probs, K=new_K, addcovar=new_addcovar, snps=snps)})
	  
	  if (value == "lod") {
	    max_A = max(qtl_perm$lod$A[,'lod'])
		if ('X' %in% names(qtl_perm$lod)) {
 		  max_X = max(qtl_perm$lod$X[,'lod'])
		} else {
		  max_X = NA
		}
		perm_max_val = max(max_A, max_X, na.rm=T)
	  } else {
	    max_A = min(qtl_perm$lod$A[,'p'])
		if ('X' %in% names(qtl_perm$lod)) {
		  max_X = min(qtl_perm$lod$X[,'p'])
		} else {
		  max_X = NA
		}
		perm_max_val = min(max_A, max_X, na.rm=T)
	  }
	  perm_max_val
	}
	gc(verbose=F)
	stopImplicitCluster()
  }
  
  ## Return max values for all permutations
  cat(' ... Done. \n')
  return(max_vals)
}
attr(scanone.perm_v2, 'help') = "
This function permutes the genomes of the mapping population and returns the max LOD
(or min P-value) for each permutation.

Parameters:
pheno: the following 6 parameters are those passed to scanone()
pheno.col:
probs:
K:
addcovar:
snps:
value: the value to collect from each permutation ('lod' or 'p'; default='lod').
B: number of permutations to perform (default=100)
parallel: logical indicating whether parallel processes should be run (default=F), 
  or numeric indicating number of cores to use. If set to TRUE, detectCores() will be used
  to automatically start n-1 parallel processes.
...: optional parameters to pass to foreach

Returns:
A vector of the specified return values for all permutations
"
  

# Added marker_pos (GC)
prob.plot <- function(pheno, pheno.col, probs, marker, marker_position=NULL) {
  ## Sort the pheno dataframe
  pheno = pheno[order(pheno[, pheno.col]),]
  
  ## Remove NAs from pheno data
  keep = rownames(pheno[!is.na(pheno[,pheno.col]),])
  if (length(keep) < 3) {
	stop("Too many missing phenotype values!")
  }
  probmat = probs[keep, , marker]
  
  grays = gray(seq(from=0, to=0.94, by=0.01))
  oranges = brewer.pal(9, 'Oranges')[3:7]
  image(z=probmat, zlim=c(0, 1), col=rev(c(grays, oranges, '#FFFFFF')), axes=FALSE)
  #image(z=probmat, zlim=c(0, 1), col=rev(gray(seq(from=0, to=1, by=0.01))), axes=FALSE)
  states <- colnames(probmat)
  nk <- ncol(probmat)
  axis(2, at=seq(from=0, to=1, len=nk), labels=states, las=1, tck=0)
  
  ## Calculate z-scores
  y = pheno[keep, pheno.col]
  z <- (y - mean(y))/sd(y)
  
  axis(3, at=ecdf(z)(pretty(z)), labels=pretty(z))
  mtext("z-score scale", side=3, at=0.5, line=2)
  if (!is.null(marker_position)) {
    mtext(paste0("Marker: ", marker," at ", marker_position), side=1, at=0.5, line=3)
  } else {
    mtext(paste0("Marker: ", marker), side=1, at=0.5, line=3)
  }
  xpos <- c(0, mean(y<mean(y)), 1)
  xtxt <- signif(c(min(y), mean(y), max(y)), 4)
  xtxt <- ifelse(abs(xtxt) < 1e-8, 0, xtxt)
  axis(1, at=xpos, labels=xtxt)
  box()
}
attr(prob.plot, 'help') = "
This function creates a probability plot for the mapping population.

Parameters:
pheno: The phenotype dataframe
pheno.col: The column in the pheno dataframe to use as the phenotype
probs: The 3D probability array
marker: The ID of the marker to plot
marker_position: optional; text containing the marker position (e.g. 'Chr1: 54.1Mb')
"

coefplot_v2 = function(doqtl, chr = 1, start=NULL, end=NULL, stat.name = "LOD", remove.outliers=NULL, conf.int = FALSE, legend = TRUE, colors = "DO", sex, ...) {

  old.par = par(no.readonly = TRUE)

  cross = attr(doqtl, "cross")
  if(is.null(cross)) {
    if(colors[1] == "DO") {    
      colors = do.colors
    } else if(colors[1] == "HS") {
      colors = hs.colors
    } # else if(colors[1] == "HS")
  } else {
    if(cross == "DO") {    
      colors = do.colors
    } else if(cross == "HS") {
      colors = hs.colors
    } # else if(cross == "HS")
  } # else

  num.founders = nrow(colors)
  call = match.call()

  # Keep only the founder coefficients from the coef.matrix.
  lod = NULL
  coef = NULL
  if(chr == "X") {
    if(missing(sex)) {
      stop("Sex (either M or F) must be specified on X chromosome.")
    } # if(missing(sex))
    lod  = doqtl$lod$X
	if (is.null(start) | is.null(end)) {
	  lod = lod[lod[,2] == chr, ]
	} else {
      lod = lod[lod[,2] == chr & lod[,3] >= start & lod[,3] <= end,]
	}
    coef = doqtl$coef$X
    if(sex == "F") {
      columns = match(paste("F", colors[,1], sep = "."), colnames(coef))
    } else {
      columns = match(paste("M", colors[,1], sep = "."), colnames(coef))
    } # else
    columns = columns[!is.na(columns)]
    coef = coef[,c(1,columns)]
    colnames(coef)[1] = "A"
    colnames(coef) = sub("^[MF]\\.", "", colnames(coef))
    coef = coef[rownames(coef) %in% lod[,1],]
	
	## Remove outliers
  if (!is.null(remove.outliers) & is.numeric(remove.outliers)) {
	  coef[abs(coef) > remove.outliers | is.na(coef)] = 0
  }
  
    # Center the coefficient values.
    coef[,2:ncol(coef)] = coef[,2:ncol(coef)] + coef[,1]
    coef = coef - rowMeans(coef)
  } else {
    lod = doqtl$lod$A
	if (is.null(start) | is.null(end)) {
	  lod = lod[lod[,2] == chr, ]
	} else {
      lod = lod[lod[,2] == chr & lod[,3] >= start & lod[,3] <= end,]
	}
	## Remove markers with outlier effects
	#coef_subset = doqtl$coef$A[rownames(doqtl$coef$A) %in% lod[,1], 3:9]
	#markers_to_remove = apply(coef_subset, 1, function(x) {any(abs(x) > remove.outliers | is.na(x))})
	#lod = lod[!markers_to_remove,]
	
    intercept = doqtl$coef$A[,1]
    coef = doqtl$coef$A[,(ncol(doqtl$coef$A)-num.founders+1):ncol(doqtl$coef$A)]
    
	## Remove outliers
  if (!is.null(remove.outliers) & is.numeric(remove.outliers)) {
    coef[abs(coef) > remove.outliers | is.na(coef)] = 0
  }
	
	coef[,1] = intercept
    colnames(coef)[1] = "A"
    coef = coef[rownames(coef) %in% lod[,1],]
	
    # Center the coefficient values.
    coef[,2:ncol(coef)] = coef[,2:ncol(coef)] + coef[,1]
    coef = coef - rowMeans(coef)
  } # else 
  # Verify that the SNP IDs in the lod & coef matrices match.
  if(!all(lod[,1] == rownames(coef))) {
    stop(paste("The SNP IDs in column 1 of the qtl data frame must match",
         "the SNP IDs in the rownames of the coef matrix."))
  } # if(!all(lod[,1] == rownames(coef)))
  # Verify that the coefficient column names are in the colors.
  if(!all(colnames(coef) %in% colors[,1])) {
    stop(paste("The founder names in the colnames of the coefficient matrix",
         "must be in column 1 of the colors matrix."))
  } # if(!all(colnames(coef) %in% colors[,1]))
  # Convert the chromosome locations to Mb.
  if(max(lod[,3], na.rm = TRUE) > 200) {
    lod[,3] = lod[,3] * 1e-6
  } # if(max(lod[,3], na.rm = TRUE) > 200)
  # Set the layout to plot the coefficients on top and the p-values on the 
  # bottom.
  layout(mat = matrix(1:2, 2, 1), heights = c(0.66666667, 0.3333333))
  par(font = 2, font.lab = 2, font.axis = 2, las = 1, plt =
      c(0.12, 0.99, 0, 0.85), xaxs = "i", lwd = 2)
  # Plot the coefficients.
  plot(lod[,3], coef[,colors[1,1]], type = "l", col = colors[1,3], lwd = 2,
       ylim = c(min(coef, na.rm = TRUE), max(coef * 2, na.rm = TRUE)), xlab = 
       paste("Chr", chr, "(Mb)"), ylab = "Founder Effects", axes = FALSE, ...)
  #abline(v = 0:20 * 10, col = "grey80")
  for(i in 1:nrow(colors)) {
    points(lod[,3], coef[,colors[i,1]], type = "l", col = colors[i,3],
           lwd = 2)
  } # for(i)
  # Draw a legend for the founder names and colors.
  if(legend) {
    legend.side = "topleft"
    if(which.max(lod[,7]) < nrow(lod) / 2) {
      legend.side = "topright"
    } # if(which.max(apply(coef, 1, max)) < nrow(lod) / 2)
    legend(legend.side, colors[,2], col = colors[,3], lty = 1, lwd = 2,
           x.intersp = 0.75, y.intersp = 0.75, bg = "white", cex = 0.8)
  } # if(legend)
  # Add the axis.
  axis(2)
  # Plot a rectangle around the plot.
  par(xpd = NA)
  usr = par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], lwd = 2)
  par(xpd = FALSE)
  # Plot the mapping statistic.
  par(plt = c(0.12, 0.99, 0.35, 1))
  # Single chromosome plot.
  plot(lod[,3], lod[,7], type = "l", lwd = 2, xlab = "",
       ylab = stat.name, ...)
  #abline(v = 0:20 * 10, col = "grey80")
  points(lod[,3], lod[,7], type = "l", lwd = 2)
  # Shade the confidence interval.
  if(conf.int) {
    interval = bayesint_v2(doqtl, chr = chr)
    usr = par("usr")
    rect(interval[1,3], usr[3], interval[3,3], usr[4], col = rgb(0,0,1,0.1), 
         border = NA)
  } # if(!is.na(conf.int))
  mtext(paste("Chr", chr, "(Mb)"), 1, 2)
  usr = par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], lwd = 2)
  par(old.par)
}
attr(coefplot_v2, 'help') = "
A slight modification of the coefplot() DOQTL function that ignores very large effects from rare haplotypes.

Parameters:
doqtl: A doqtl object; results from scanone()
chr: the chromosome to plot
start: the start position (Mb) of the interval to be plotted
end: the end position (Mb) of the interval to be plotted
stat.name: default is 'LOD'
remove.outliers: numeric; sets coefficients above this value to zero (default is NULL--do nothing)
conf.int: draw a box showing the QTL confidence interval (calculated using bayesint_v2())
legend: draw a legend (default is TRUE)
colors: default is 'DO'
sex: required only for plotting the X chromosome
...: additional arguments passed to plot()
"


bayesint_v2 = function (qtl, chr, prob = 0.95, expandtomarkers = TRUE, ignore_genetic_dist=TRUE) 
{
    if (missing(qtl)) {
        stop("bayesint: The qtl cannot be null. Please supply a QTL object.")
    }
    if (missing(chr)) {
        stop(paste("bayesint: The chromosome cannot be null."))
    }
    else if (!chr %in% c(1:19, "X")) {
        stop(paste("bayesint: The chromosome must be 1 to 19 or X."))
    }
    if (prob < 0 | prob > 1) {
        stop(paste("bayesint: The probability must between 0 and 1."))
    }
    old.warn = options("warn")$warn
    options(warn = -1)
    if (is.numeric(chr)) {
        qtl = qtl$lod$A
    }
    else {
        qtl = qtl$lod$X
    }
    options(warn = old.warn)
    qtl[, 1] = as.character(qtl[, 1])
    qtl[, 2] = as.character(qtl[, 2])
    qtl[, 3] = as.numeric(qtl[, 3])
    qtl[, 7] = as.numeric(qtl[, 7])
    qtl = qtl[qtl[, 2] == chr, ]
    pos = qtl[, 3]
    if (any(is.na(pos))) {
        remove = which(is.na(pos))
        qtl = qtl[-remove, ]
        pos = pos[-remove]
    }
	
    breaks = approx(x = pos, y = 10^qtl[, 7], xout = seq(pos[1], 
        pos[length(pos)], length.out = 1e+05))
    widths = diff(breaks$x)
    heights = breaks$y[-1] + breaks$y[-length(breaks$y)]
    trapezoids = 0.5 * heights * widths
    trapezoids = trapezoids/sum(trapezoids)
    ord = order(breaks$y[-length(breaks$y)], decreasing = TRUE)
    wh = min(which(cumsum(trapezoids[ord]) >= prob))
    int = range(ord[1:wh])
	if (ignore_genetic_dist) {
		left.snp = c(NA, qtl[1, 2], breaks$x[int][1], NA, approx(qtl[, 3], qtl[, 
			5], breaks$x[int][1])$y, approx(qtl[, 3], qtl[, 6], breaks$x[int][1])$y, 
			approx(qtl[, 3], qtl[, 7], breaks$x[int][1])$y)
		max.snp = qtl[which.max(qtl[, 7]), ]
		right.snp = c(NA, qtl[1, 2], breaks$x[int][2], NA, approx(qtl[, 3], qtl[, 
			5], breaks$x[int][2])$y, approx(qtl[, 3], qtl[, 6], breaks$x[int][2])$y, 
			approx(qtl[, 3], qtl[, 7], breaks$x[int][2])$y)
			
			## TESTING/DEBUGGING
			#print(dim(qtl))
			#print(which.max(qtl[, 7]))
	} else {
		left.snp = c(NA, qtl[1, 2], breaks$x[int][1], approx(qtl[, 
			3], qtl[, 4], breaks$x[int][1])$y, approx(qtl[, 3], qtl[, 
			5], breaks$x[int][1])$y, approx(qtl[, 3], qtl[, 6], breaks$x[int][1])$y, 
			approx(qtl[, 3], qtl[, 7], breaks$x[int][1])$y)
		max.snp = qtl[which.max(qtl[, 7]), ]
		right.snp = c(NA, qtl[1, 2], breaks$x[int][2], approx(qtl[, 
			3], qtl[, 4], breaks$x[int][2])$y, approx(qtl[, 3], qtl[, 
			5], breaks$x[int][2])$y, approx(qtl[, 3], qtl[, 6], breaks$x[int][2])$y, 
			approx(qtl[, 3], qtl[, 7], breaks$x[int][2])$y)
	}
    if (expandtomarkers) {
        left.snp = qtl[max(which(breaks$x[int][1] >= qtl[, 3])), 
            ]
        max.snp = qtl[which.max(qtl[, 7]), ]
        right.snp = qtl[min(which(breaks$x[int][2] <= qtl[, 3])), 
            ]
    }
    retval = rbind(left.snp, max.snp, right.snp)
    retval[, 3] = round(as.numeric(retval[, 3]), digits = 6)
    retval[, 4] = round(as.numeric(retval[, 4]), digits = 6)
    retval[, 5] = round(as.numeric(retval[, 5]), digits = 6)
    retval[, 6] = round(as.numeric(retval[, 6]), digits = 6)
    retval$lod = as.numeric(retval[, 7])
    return(retval)
}
attr(bayesint_v2, 'help') = "
A slight modification of the bayesint() DOQTL function that does not require genetic distances for input.

Parameters:
qtl: A doqtl object; results from scanone()
chr: The column in the pheno dataframe to use as the phenotype
prob: A probability for the desired confidence interval
expandtomarkers: logical (default is TRUE)
ignore_genetic_dist: logical (default is TRUE)
"

## Below is the same code from DOQTL, except code for automatically writing results to a text file is commented out
scanone.perm = function (pheno, pheno.col = 1, probs, addcovar, intcovar, snps, model = c("additive", "full"), path = ".", nperm = 1000, return.val = c("lod", "p"))
{
  return.val = match.arg(return.val)
  if (!missing(intcovar)) {
    stop("Interactive covariates not yet implemented")
  }
  if (is.null(rownames(pheno))) {
    stop("rownames(pheno) is null. The sample IDs must be in rownames(pheno).")
  }
  if (is.character(pheno.col)) {
    pheno.col = match(pheno.col, colnames(pheno))
  }
  if (!file.exists(path)) {
    stop(paste("The path", path, "does not exist."))
  }
  probs = filter.geno.probs(probs)
  pheno = pheno[rownames(pheno) %in% dimnames(probs)[[1]], , drop = FALSE]
  probs = probs[dimnames(probs)[[1]] %in% rownames(pheno), , ]
  probs = probs[match(rownames(pheno), dimnames(probs)[[1]]), , ]
  probs = probs[, , dimnames(probs)[[3]] %in% snps[, 1]]
  if (any(dim(probs) == 0)) {
    stop(paste("There are no matching samples in pheno and probs. Please",
    "verify that the sample IDs in rownames(pheno) match the sample",
    "IDs in dimnames(probs)[[1]]."))
  }
  probs = filter.geno.probs(probs)
  snps = snps[snps[, 1] %in% dimnames(probs)[[3]], ]
  probs = probs[, , match(snps[, 1], dimnames(probs)[[3]])]
  if (!missing(addcovar)) {
    addcovar = as.matrix(addcovar)
    addcovar = addcovar[rowMeans(is.na(addcovar)) == 0, , drop = FALSE]
    pheno = pheno[rownames(pheno) %in% rownames(addcovar), , drop = FALSE]
    addcovar = addcovar[rownames(addcovar) %in% rownames(pheno), , drop = FALSE]
    addcovar = addcovar[match(rownames(pheno), rownames(addcovar)), , drop = FALSE]
    probs = probs[, , dimnames(probs)[[3]] %in% snps[, 1]]
    if (is.null(colnames(addcovar))) {
      colnames(addcovar) = paste("addcovar", 1:ncol(addcovar), sep = ".")
    }
  }
  if (!missing(intcovar)) {
    covar = as.matrix(intcovar)
    intcovar = intcovar[rowMeans(is.na(intcovar)) == 0, , drop = FALSE]
    pheno = pheno[rownames(pheno) %in% rownames(intcovar), , drop = FALSE]
    intcovar = intcovar[rownames(intcovar) %in% rownames(pheno), , drop = FALSE]
    intcovar = intcovar[match(rownames(pheno), rownames(intcovar)), , drop = FALSE]
    probs = probs[, , dimnames(probs)[[3]] %in% snps[, 1]]
    if (is.null(colnames(intcovar))) {
      colnames(intcovar) = paste("intcovar", 1:ncol(intcovar), sep = ".")
    }
  }
  for (i in pheno.col) {
    print(colnames(pheno)[i])
    p = pheno[, i]
    names(p) = rownames(pheno)
    keep = which(!is.na(p))
    if (!missing(addcovar)) {
      perms = permutations.qtl.LRS(pheno = p[keep], probs = probs[keep, , ], snps = snps, addcovar = addcovar[keep, , drop = FALSE], nperm = nperm, return.val = return.val)
    }
    else {
      perms = permutations.qtl.LRS(pheno = p[keep], probs = probs[keep, , ], snps = snps, nperm = nperm, return.val = return.val)
    }
    if (return.val == "lod") {
      perms = perms/(2 * log(10))
    }
    #write(perms, paste(path, "/", colnames(pheno)[i], ".perms.txt", sep = ""), sep = "\t")
  }
  return(perms)
}


abind = function (..., along = N, rev.along = NULL, new.names = NULL,
    force.array = TRUE, make.names = use.anon.names, use.anon.names = FALSE,
    use.first.dimnames = FALSE, hier.names = FALSE, use.dnns = FALSE)
{
    if (is.character(hier.names))
        hier.names <- match.arg(hier.names, c("before", "after",
            "none"))
    else hier.names <- if (hier.names)
        "before"
    else "no"
    arg.list <- list(...)
    if (is.list(arg.list[[1]]) && !is.data.frame(arg.list[[1]])) {
        if (length(arg.list) != 1) 
            stop("can only supply one list-valued argument for ...")
        if (make.names) 
            stop("cannot have make.names=TRUE with a list argument")
        arg.list <- arg.list[[1]]
        have.list.arg <- TRUE
    }
    else {
        N <- max(1, sapply(list(...), function(x) length(dim(x))))
        have.list.arg <- FALSE
    }
    if (any(discard <- sapply(arg.list, is.null))) 
        arg.list <- arg.list[!discard]
    if (length(arg.list) == 0) 
        return(NULL)
    N <- max(1, sapply(arg.list, function(x) length(dim(x))))
    if (!is.null(rev.along)) 
        along <- N + 1 - rev.along
    if (along < 1 || along > N || (along > floor(along) && along < 
        ceiling(along))) {
        N <- N + 1
        along <- max(1, min(N + 1, ceiling(along)))
    }
    if (length(along) > 1 || along < 1 || along > N + 1) 
        stop(paste("\"along\" must specify one dimension of the array,", 
            "or interpolate between two dimensions of the array", 
            sep = "\n"))
    if (!force.array && N == 2) {
        if (!have.list.arg) {
            if (along == 2) 
                return(cbind(...))
            if (along == 1) 
                return(rbind(...))
        }
        else {
            if (along == 2) 
                return(do.call("cbind", arg.list))
            if (along == 1) 
                return(do.call("rbind", arg.list))
        }
    }
    if (along > N || along < 0) 
        stop("along must be between 0 and ", N)
    pre <- seq(from = 1, len = along - 1)
    post <- seq(to = N - 1, len = N - along)
    perm <- c(seq(len = N)[-along], along)
    arg.names <- names(arg.list)
    if (is.null(arg.names)) 
        arg.names <- rep("", length(arg.list))
    if (is.character(new.names)) {
        arg.names[seq(along = new.names)[nchar(new.names) > 0]] <- new.names[nchar(new.names) > 
            0]
        new.names <- NULL
    }
    if (any(arg.names == "")) {
        if (make.names) {
            dot.args <- match.call(expand.dots = FALSE)$...
            if (is.call(dot.args) && identical(dot.args[[1]], 
                as.name("list"))) 
                dot.args <- dot.args[-1]
            arg.alt.names <- arg.names
            for (i in seq(along = arg.names)) {
                if (arg.alt.names[i] == "") {
                  if (object.size(dot.args[[i]]) < 1000) {
                    arg.alt.names[i] <- paste(deparse(dot.args[[i]], 
                      40), collapse = ";")
                  }
                  else {
                    arg.alt.names[i] <- paste("X", i, sep = "")
                  }
                  arg.names[i] <- arg.alt.names[i]
                }
            }
        }
        else {
            arg.alt.names <- arg.names
            arg.alt.names[arg.names == ""] <- paste("X", seq(along = arg.names), 
                sep = "")[arg.names == ""]
        }
    }
    else {
        arg.alt.names <- arg.names
    }
    use.along.names <- any(arg.names != "")
    names(arg.list) <- arg.names
    arg.dimnames <- matrix(vector("list", N * length(arg.names)), 
        nrow = N, ncol = length(arg.names))
    dimnames(arg.dimnames) <- list(NULL, arg.names)
    arg.dnns <- matrix(vector("list", N * length(arg.names)), 
        nrow = N, ncol = length(arg.names))
    dimnames(arg.dnns) <- list(NULL, arg.names)
    dimnames.new <- vector("list", N)
    arg.dim <- matrix(integer(1), nrow = N, ncol = length(arg.names))
    for (i in seq(len = length(arg.list))) {
        m <- arg.list[[i]]
        m.changed <- FALSE
        if (is.data.frame(m)) {
            m <- as.matrix(m)
            m.changed <- TRUE
        }
        else if (!is.array(m) && !is.null(m)) {
            if (!is.atomic(m)) 
                stop("arg '", arg.alt.names[i], "' is non-atomic")
            dn <- names(m)
            m <- as.array(m)
            if (length(dim(m)) == 1 && !is.null(dn)) 
                dimnames(m) <- list(dn)
            m.changed <- TRUE
        }
        new.dim <- dim(m)
        if (length(new.dim) == N) {
            if (!is.null(dimnames(m))) {
                arg.dimnames[, i] <- dimnames(m)
                if (use.dnns && !is.null(names(dimnames(m)))) 
                  arg.dnns[, i] <- as.list(names(dimnames(m)))
            }
            arg.dim[, i] <- new.dim
        }
        else if (length(new.dim) == N - 1) {
            if (!is.null(dimnames(m))) {
                arg.dimnames[-along, i] <- dimnames(m)
                if (use.dnns && !is.null(names(dimnames(m)))) 
                  arg.dnns[-along, i] <- as.list(names(dimnames(m)))
                dimnames(m) <- NULL
            }
            arg.dim[, i] <- c(new.dim[pre], 1, new.dim[post])
            if (any(perm != seq(along = perm))) {
                dim(m) <- c(new.dim[pre], 1, new.dim[post])
                m.changed <- TRUE
            }
        }
        else {
            stop("'", arg.alt.names[i], "' does not fit: should have `length(dim())'=", 
                N, " or ", N - 1)
        }
        if (any(perm != seq(along = perm))) 
            arg.list[[i]] <- aperm(m, perm)
        else if (m.changed) 
            arg.list[[i]] <- m
    }
    conform.dim <- arg.dim[, 1]
    for (i in seq(len = ncol(arg.dim))) {
        if (any((conform.dim != arg.dim[, i])[-along])) {
            stop("arg '", arg.alt.names[i], "' has dims=", paste(arg.dim[, 
                i], collapse = ", "), "; but need dims=", paste(replace(conform.dim, 
                along, "X"), collapse = ", "))
        }
    }
    if (N > 1) 
        for (dd in seq(len = N)[-along]) {
            for (i in (if (use.first.dimnames) 
                seq(along = arg.names)
            else rev(seq(along = arg.names)))) {
                if (length(arg.dimnames[[dd, i]]) > 0) {
                  dimnames.new[[dd]] <- arg.dimnames[[dd, i]]
                  if (use.dnns && !is.null(arg.dnns[[dd, i]])) 
                    names(dimnames.new)[dd] <- arg.dnns[[dd, 
                      i]]
                  break
                }
            }
        }
    for (i in seq(len = length(arg.names))) {
        if (arg.dim[along, i] > 0) {
            dnm.along <- arg.dimnames[[along, i]]
            if (length(dnm.along) == arg.dim[along, i]) {
                use.along.names <- TRUE
                if (hier.names == "before" && arg.names[i] != 
                  "") 
                  dnm.along <- paste(arg.names[i], dnm.along, 
                    sep = ".")
                else if (hier.names == "after" && arg.names[i] != 
                  "") 
                  dnm.along <- paste(dnm.along, arg.names[i], 
                    sep = ".")
            }
            else {
                if (arg.dim[along, i] == 1) 
                  dnm.along <- arg.names[i]
                else if (arg.names[i] == "") 
                  dnm.along <- rep("", arg.dim[along, i])
                else dnm.along <- paste(arg.names[i], seq(length = arg.dim[along, 
                  i]), sep = "")
            }
            dimnames.new[[along]] <- c(dimnames.new[[along]], 
                dnm.along)
        }
        if (use.dnns) {
            dnn <- unlist(arg.dnns[along, ])
            if (length(dnn)) {
                if (!use.first.dimnames) 
                  dnn <- rev(dnn)
                names(dimnames.new)[along] <- dnn[1]
            }
        }
    }
    if (!use.along.names) 
        dimnames.new[along] <- list(NULL)
    out <- array(unlist(arg.list, use.names = FALSE), dim = c(arg.dim[-along, 
        1], sum(arg.dim[along, ])), dimnames = dimnames.new[perm])
    if (any(order(perm) != seq(along = perm))) 
        out <- aperm(out, order(perm))
    if (!is.null(new.names) && is.list(new.names)) {
        for (dd in seq(len = N)) {
            if (!is.null(new.names[[dd]])) {
                if (length(new.names[[dd]]) == dim(out)[dd]) 
                  dimnames(out)[[dd]] <- new.names[[dd]]
                else if (length(new.names[[dd]])) 
                  warning(paste("Component ", dd, " of new.names ignored: has length ", 
                    length(new.names[[dd]]), ", should be ", 
                    dim(out)[dd], sep = ""))
            }
            if (use.dnns && !is.null(names(new.names)) && names(new.names)[dd] != 
                "") 
                names(dimnames(out))[dd] <- names(new.names)[dd]
        }
    }
    if (use.dnns && !is.null(names(dimnames(out))) && any(i <- is.na(names(dimnames(out))))) 
        names(dimnames(out))[i] <- ""
    out
}

