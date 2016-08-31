
## Load R functions and libraries
source('./scripts/rix_qtl_mapping_functions.r')

## First get all RIX matings
rix_list = read.xls('./data/BOBs_Corrected_16_06_22.xlsx', header=T, as.is=T)
rix_list = rix_list[,1:2]
colnames(rix_list) = c('Mating', 'CC_Mating')
head(rix_list)

## Set directory containing founder 36-state probability files, and get filenames
cc_dir = '/Users/mooneymi/Documents/SIG/Mapping/CC'
cc_prob_files = list.files(cc_dir, pattern="CC...-.*b38.*\\.csv")
length(cc_prob_files)
cc_prob_files[1:5]

## Get vectors of parental strains
rix_strains = unlist(strsplit(rix_list$Mating, 'x'))
cc_strains = unlist(strsplit(rix_list$CC_Mating, 'x'))

## Check strains (some cc_strains will be duplicated -- re-derived lines)
dup_idx = duplicated(rix_strains)
cc_strains = cc_strains[!dup_idx]
names(cc_strains) = rix_strains[!dup_idx]
length(cc_strains)
cc_strains[1:5]

## Check that all RIX strains mapped to a CC strains
sum(is.na(cc_strains))

## Set directory containing the 8-state probability files
## These files can be created using the collapse_probs() function (see next code block)
rix_dir = '/Users/mooneymi/Documents/SIG/Mapping/RIX'

## If you haven't already done so, the following code can be used to 
## create the 8-state probability files (you should only have to do this once)
## File names will be the RIX IDs (not CC IDs). 
## Note: you may have to make copies of some files for re-derived lines 
## (e.g. 8042 and 18042 are the same), if this wasn't handled above
# for (i in 1:length(cc_strains)) {
#     ## Get RIX ID
#     rix = names(cc_strains[i])
#     
#     ## Create and save the 8-state probabilities for each parental strain
#     cc_file = file.path(cc_dir, cc_prob_files[grep(paste0(cc_strains[i], '-.*\\.csv'), cc_prob_files)])
#     cc_8state = collapse_probs(cc_file)
#     
#     ## Save Y and M chromosomes separately
#     autosome_x_markers = cc_8state$marker[!cc_8state$chromosome %in% c('Y', 'M')]
#     y_m_markers = cc_8state$marker[cc_8state$chromosome %in% c('Y', 'M')]
#     file_name = file.path(rix_dir, paste0(rix, '.csv'))
#     write.table(cc_8state[autosome_x_markers, ], file=file_name, row.names=F, col.names=T, sep=',', quote=F)
#     file_name = file.path(rix_dir, paste0(rix, '_Y_M.csv'))
#     write.table(cc_8state[y_m_markers, ], file=file_name, row.names=F, col.names=T, sep=',', quote=F)
# }

## Create 3D probability array - Males
## Note: the number of matings included will be printed
model.probs = make_rix_model_probs(rix_list$Mating, rix_dir, 'M')
dim(model.probs)
names(dimnames(model.probs))

## Check model.probs object
model.probs[1,,1:5]

## Get vector of all markers
all_markers = dimnames(model.probs)[[3]]

## Check if any markers sum to 0
markers_zero_idx = which(apply(model.probs, 1, colSums) == 0)
markers_zero = rep(all_markers, length(samples))[markers_zero_idx]
markers_zero = markers_zero[!is.na(markers_zero)]
markers_zero = unique(markers_zero)
length(markers_zero)

## Remove markers with all zeros
model.probs = model.probs[, , setdiff(all_markers, markers_zero)]
dim(model.probs)

## Make all probabilities non-zero
model.probs[model.probs == 0] = 1e-20
model.probs[1,,1:5]

## Save genome cache - Males
mapping_dir = '/Users/mooneymi/Documents/SIG/Mapping'
save(model.probs, file=file.path(mapping_dir, 'rix_universal_model_prob_males_27-Jun-2016.rda'))

## Create 3D probability array - Females
model.probs = make_rix_model_probs(rix_list$Mating, rix_dir, 'F')
dim(model.probs)
names(dimnames(model.probs))

## Check model.probs object
model.probs[1,,1:5]

## Remove markers with all zeros
model.probs = model.probs[, , setdiff(all_markers, markers_zero)]
dim(model.probs)

## Make all probabilities non-zero
model.probs[model.probs == 0] = 1e-20
model.probs[1,,1:5]

## Save genome cache - Females
save(model.probs, file=file.path(mapping_dir, 'rix_universal_model_prob_females_27-Jun-2016.rda'))
