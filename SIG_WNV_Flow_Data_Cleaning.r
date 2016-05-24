
## Load functions for parsing the flow cytometry spreadsheets
## The gdata library is necessary for reading Excel spreadsheets; it will be loaded as well.
source('./scripts/flow_data_cleaning_functions.r')

## View help documentation on the functions listed above
describe(read_flow_exp_file)

## Load all expected flow variables (expected column names)
flow_cn = read.xls(xls="./data/WNV_Data_Dictionary.xlsx", sheet="Flow Data", as.is=T)
flow_cn = flow_cn[,1]

## Move to the directory holding the data
flow_dir = "/Users/mooneymi/Documents/MyDocuments/SystemsImmunogenetics/WNV/Lund_Flow_fixed_Apr_26"

## Get a list of data files to read (in this case all flow spreadsheets begin with the prefix 'Expt')
flow_files = list.files(flow_dir, pattern="Expt.*\\.xls")

print(flow_files)

## Iterate through all the files, parse each, and merge all data into a single dataframe
i = 1
for (file in flow_files) {
    print(file)
    flow_dat = read_flow_exp_file(file.path(flow_dir, file), flow_cn)
    
    ## Check if there are any unexpected columns
    new_columns = setdiff(colnames(flow_dat), flow_cn)
    if (length(new_columns) > 0) {
        flow_cn = c(flow_cn, new_columns)
    }
    if (i > 1) {
        ## Fill extra columns with NAs
        for (col in new_columns) {
            flow_all[,col] = NA
        }
        ## Merge data
        flow_all = rbind(flow_all[,flow_cn], flow_dat[,flow_cn])
    } else {
        flow_all = flow_dat
    }
    i = i + 1
}

## Check the dimensions of the dataframe
dim(flow_all)

## Check that all expected columns are present
setdiff(flow_cn, colnames(flow_all))

flow_all[1:10,1:11]

## Order columns, add Lab column and fix formatting
flow_all = flow_all[, flow_cn]
flow_all$Lab = "Lund"

flow_all$ID = gsub(" ", "", flow_all$ID)
flow_all$ID = gsub("X", "x", flow_all$ID)
flow_all$Mating = gsub(" ", "", flow_all$Mating)
flow_all$Mating = gsub("X", "x", flow_all$Mating) 
flow_all$UW_Line = as.numeric(flow_all$UW_Line)
flow_all$Timepoint = as.numeric(flow_all$Timepoint)

## For validation data mock animals do not have IDs, set to NA
flow_all$ID[is.na(flow_all$RIX_ID)] = NA

flow_all[1:10,1:11]

## Check for duplicate IDs
new_flow_ids = paste(flow_all$ID, flow_all$Tissue, sep='_')
sum(duplicated(new_flow_ids))

## Read the previously cleaned data
## Note: you will have to change the file path
cleaned_dir = '/Users/mooneymi/Documents/MyDocuments/SystemsImmunogenetics/WNV/Cleaned_Data_Releases/23-Mar-2016'
flow_prev = read.xls(file.path(cleaned_dir, 'Lund_Flow_21-Mar-2016_final.xlsx'), 
                     header=T, as.is=T, na.strings=c(""," ", "NA", "#DIV/0!"))

dim(flow_prev)

## Check for duplicate IDs
dup_ids1 = intersect(flow_prev$ID[flow_prev$Tissue=='brain'], flow_all$ID[flow_all$Tissue=='brain'])
dup_ids2 = intersect(flow_prev$ID[flow_prev$Tissue=='spleen'], flow_all$ID[flow_all$Tissue=='spleen'])

## Overwrite old data with new
idx1 = which(flow_prev$ID %in% dup_ids1 & flow_prev$Tissue=='brain')
idx2 = which(flow_prev$ID %in% dup_ids2 & flow_prev$Tissue=='spleen')
idx_dups = c(idx1, idx2)
print(length(idx_dups))

idx_to_keep = setdiff(1:nrow(flow_prev), idx_dups)
flow_prev = flow_prev[idx_to_keep,]
flow_all = rbind(flow_prev[, flow_cn], flow_all[, flow_cn])

## Check the dimensions of the dataframe
dim(flow_all)

## Change all data columns to numeric
for (i in 11:277) {
    flow_all[,i] = as.numeric(flow_all[,i])
}

## Calculate cell counts and ratios
flow_full = flow_all
flow_full = calc_treg_counts(flow_full)
flow_full = calc_tcell_counts(flow_full)
flow_full = calc_ics_counts(flow_full)
flow_full = calc_ics_percent_ratios(flow_full)
flow_full = calc_ics_count_ratios(flow_full)
flow_full = clean_inf_nan(flow_full)

dim(flow_full)

## Save the data file
write.table(flow_all, file=file.path(flow_dir, 'Lund_Flow_12-May-2016_final.txt'), 
            col.names=T, row.names=F, quote=T, sep='\t', na='')

## Save the full data file
write.table(flow_full, file=file.path(flow_dir, 'Lund_Flow_Full_12-May-2016_final.txt'), 
            col.names=T, row.names=F, quote=T, sep='\t', na='')
save(flow_full, file=file.path(flow_dir, 'lund_flow_full_12-May-2016_final.rda'))
