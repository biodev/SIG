
## Load R functions and libraries
source('./scripts/rix_qtl_mapping_functions.r')

## Read phenotype data
cleaned_data_dir = '/Users/mooneymi/Documents/SIG/WNV/Cleaned_Data_Releases/14-Jun-2016/'
pheno_dir = '/Users/mooneymi/Documents/SIG/Mapping/phenotypes'
mapping_dir = '/Users/mooneymi/Documents/SIG/Mapping'

## For this example we're working with qPCR data
pheno = read.xls(file.path(cleaned_data_dir, 'Gale_qPCR_byMouse_14-Jun-2016_final.xlsx'), header=T, as.is=T)
dim(pheno)

## The phenotype will be D12 WNV (viral transcript) in the brain
pheno = pheno[with(pheno, Virus=='WNV' & Timepoint==12 & Experiment=='WNV' & Tissue=='Brain'), ]
dim(pheno)

## Save the phenotype dataframe (optional)
#save(pheno, file=file.path(pheno_dir, 'wnv_dCt_brain_d12_wnv_17-Jun-2016.rda'))

## View the phenotype distribution accross the mapping population
library(lattice)
dotplot(reorder(pheno[,'UW_Line'], pheno[,'dCt'], mean, na.rm=T) ~ 
        pheno[,'dCt'] | pheno[,'Virus'], 
        panel = function(x,y,...) {panel.dotplot(x,y,...); panel.abline(v=0, col.line="red")}, 
        pch=19, ylab='UW Line', xlab="D12 WNV dCt")

## Sort pheno dataframe and set rownames
pheno = pheno[with(pheno, order(Mating, RIX_ID)),]
rownames(pheno) = pheno$ID

## Add sex column
pheno$sex = 'M'

## Create covariate dataframe (must include sex)
covar = data.frame(sex = as.numeric(pheno$sex == 'M'))
rownames(covar) = pheno$ID

## Get IDs and Matings for each sample
samples = pheno$ID
matings = unlist(lapply(strsplit(samples, '_'), function(x) {x[1]}))
matings = unique(matings)

## Read strain ID mapping file (with Oas1b status for WNV mapping analyses)
## Note: use the Mx1 file (linked at top) for Flu analyses
strain_map = read.delim(file.path('./data', 'Oas1b_status_recoded.txt'), header=T, as.is=T, sep='\t')
head(strain_map)

## Check that matings match the strain mapping file
setdiff(matings, strain_map$Mating)

## Update covar with Oas1b status (only for WNV analyses)
covar$Mating = pheno$Mating
covar$Oas1b = sapply(covar$Mating, function(x) {if (x %in% strain_map$Mating) strain_map$Oas1b[strain_map$Mating == x] else NA})
covar$Oas1b_High = sapply(covar$Mating, function(x) {if (x %in% strain_map$Mating) strain_map$Oas1b_High[strain_map$Mating == x] else NA})
covar$Oas1b_Mod = sapply(covar$Mating, function(x) {if (x %in% strain_map$Mating) strain_map$Oas1b_Mod[strain_map$Mating == x] else NA})
covar$Oas1b_Low = sapply(covar$Mating, function(x) {if (x %in% strain_map$Mating) strain_map$Oas1b_Low[strain_map$Mating == x] else NA})
pheno$Oas1b = covar$Oas1b
pheno$Oas1b_High = covar$Oas1b_High
pheno$Oas1b_Mod = covar$Oas1b_Mod
pheno$Oas1b_Low = covar$Oas1b_Low

## Load universal model probabilities (loads a model.probs object containing all RIX lines)
data_dir = './xdata'
load(file.path(data_dir, 'rix_universal_model_prob_males_27-Jun-2016.rda'))

## Create model.probs array specific to the mapping population
model.probs = model.probs[pheno$Mating, , ]
dimnames(model.probs)[[1]] = pheno$ID

## Check model.probs object
dim(model.probs)
names(dimnames(model.probs))
dim(model.probs)[1] == dim(pheno)[1]
model.probs[1,,1:5]

## Fix very small genotype probabilities
model.probs[model.probs < 0.005] = 1e-20

## Check model.probs object
model.probs[1,,1:5]

## Create kinship probability matrix
K = kinship.probs(model.probs)

## Garbage collection might help free memory
gc()

## Check kinship matrix
K[1:5, 1:5]

## First get marker positions
marker_pos = read.csv(file.path(data_dir, 'CC001-Uncb38V01.csv'), as.is=T)
marker_pos = marker_pos[,1:3]
marker_pos$position_cM = NA
head(marker_pos)

## Run QTL scan
qtl = scanone(pheno=pheno, pheno.col='dCt', probs=model.probs, K=K, 
              addcovar=covar[, c('sex'), drop=F], snps=marker_pos)

## Run permutations to calculate significance threshold
## Note: nperm=10 is for demonstration purposes only. Normally, at least 1000 permutations should be done.
perms = scanone.perm(pheno=pheno, pheno.col='dCt', probs=model.probs, 
                     addcovar=covar[, c('sex'), drop=F], snps=marker_pos, nperm = 10)

## Get the significance thresholds
sig.thr = get.sig.thr(perms[,,'dCt'], alpha=0.05)
sugg.thr = get.sig.thr(perms[,,'dCt'], alpha=0.1)

## Plot QTL results
plot(qtl, sig.thr = sig.thr, main = 'Day 12 Brain WNV dCt')

## Add lines for suggestive peaks (alpha = 0.1)
lines(x = c(0, 2462), y = rep(sugg.thr[1,"A"], 2), lwd = 2, lty=2, col = 'red')
lines(x = c(2462, 2633), y = rep(sugg.thr[1,"X"], 2), lwd = 2, lty=2, col = 'red')

## Identify the Bayes Credible Interval for the peak on chromosome 5
interval = bayesint_v2(qtl, chr = 5)
interval

## Create founder effects (coefficient) plot for peak interval
coefplot_v2(qtl, chr=5, start=interval[1,3], end=interval[3,3])

## Get marker at peak
get_min_marker(qtl, chr=5)

## Create probability plot
prob.plot(pheno=pheno, pheno.col='dCt', probs=model.probs, marker='UNC10043432', qtl=qtl)
