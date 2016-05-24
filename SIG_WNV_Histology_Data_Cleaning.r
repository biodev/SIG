
library(gdata)

## Set data directories
cleaned_dir = '~/Documents/SIG/WNV/Cleaned_Data_Releases/15-Jan-2016'
hist_dir = '~/Documents/SIG/WNV/Histology/fixed_files'

## Load column names
hist_cn = read.xls(xls="./data/WNV_Data_Dictionary.xlsx", sheet="Histology Data", as.is=T)
hist_cn = hist_cn[,1]

## Load cleaned weight data
gale_weight = read.xls(file.path(cleaned_dir, 'Gale_Weight_13-Jan-2016_final.xlsx'), header=T, 
                       as.is=T, na.strings=c(""," ", "NA", "na", "#DIV/0!"))

## Load previously cleaned histology data
hist_prev = read.xls(file.path(cleaned_dir, 'Gale_Histology_5-Jan-2016_final.xlsx'), header=T, 
                     as.is=T, na.strings=c(""," ", "NA", "na", "#DIV/0!"))

## Get files to be processed
hist_files = list.files(hist_dir, pattern=".*\\.xlsx")
print(hist_files)

## Read all files and combine into a single dataframe
print(hist_files[1])
hist_data_v2 = read.xls(file.path(hist_dir, hist_files[1]), header=T, 
                        as.is=T, na.strings=c(""," ", "NA", "na", "#DIV/0!"))
colnames(hist_data_v2) = hist_cn[c(9, c(13:35))]

for (file in hist_files[2:length(hist_files)]) {
  print(file)
  hist_data = read.xls(file.path(hist_dir, file), header=T, 
                       as.is=T, na.strings=c(""," ", "NA", "na", "#DIV/0!"))
  colnames(hist_data) = hist_cn[c(9, c(13:35))]
  hist_data_v2 = rbind(hist_data_v2, hist_data)
}

## Check dimensions of dataframe
dim(hist_data_v2)

## Create missing columns
for (cn in setdiff(hist_cn, colnames(hist_data_v2))) {
  print(cn)
  hist_data_v2[,cn] = NA
}

## Order columns
hist_data_v2 = hist_data_v2[, hist_cn]

## Remove any leading or trailing spaces from slide_label
hist_data_v2$slide_label = trim(hist_data_v2$slide_label)

## Update UW Line
hist_data_v2$UW_Line = sapply(hist_data_v2$slide_label, function(x){unlist(strsplit(x, " "))[1]})
hist_data_v2$UW_Line = as.numeric(hist_data_v2$UW_Line)

## Update UWID
hist_data_v2$UWID = sapply(hist_data_v2$slide_label, function(x){unlist(strsplit(x, " "))[2]})
hist_data_v2$UWID = gsub("m", "M", hist_data_v2$UWID)

## Update Lab
hist_data_v2$Lab = 'Gale'

## Update time points
hist_data_v2$Timepoint = sapply(hist_data_v2$UWID, function(x){unlist(strsplit(x, "\\."))[2]})

## Use time points to update virus
hist_data_v2$Virus = 'WNV'
hist_data_v2$Virus[grepl("M", hist_data_v2$Timepoint)] = 'Mock'
hist_data_v2$Virus[grepl("m", hist_data_v2$Timepoint)] = 'Mock'

## Update time points
hist_data_v2$Timepoint = gsub("M", "", hist_data_v2$Timepoint)
hist_data_v2$Timepoint = gsub("m", "", hist_data_v2$Timepoint)
hist_data_v2$Timepoint = as.numeric(hist_data_v2$Timepoint)

## Update Mating and RIX_ID
for (i in 1:dim(hist_data_v2)[1]) {
  line = hist_data_v2$UW_Line[i]
  if (!is.na(line)) {
    mating = gale_weight$Mating[with(gale_weight, UW_Line==hist_data_v2$UW_Line[i])]
    if (length(mating)>0) {
      hist_data_v2$Mating[i] = mating[1]
      rix_id = gale_weight$RIX_ID[with(gale_weight, UW_Line==hist_data_v2$UW_Line[i] & UWID==hist_data_v2$UWID[i])]
      if (length(rix_id)>0) {
        if (length(rix_id)>1) {
          ## If more than one ID matches the UWID, print the IDs
          print(paste0("line=", mating[1], "; rix_id=", rix_id))
        }
        hist_data_v2$RIX_ID[i] = rix_id[1]
      }
    }
  }
}

## Check Matings and IDs of the discrepancies
## Make manual corrections if necessary
hist_data_v2[hist_data_v2$Mating=='3260x1566' & !is.na(hist_data_v2$Mating), 1:12]

## Update ID
hist_data_v2$ID = paste(hist_data_v2$Mating, hist_data_v2$RIX_ID, sep="_")
hist_data_v2$ID[grepl("NA", hist_data_v2$ID)] = NA

## Update tissue
hist_data_v2$Tissue = 'Brain'
for (i in 1:length(hist_data_v2$slide_label)) {
  label = hist_data_v2$slide_label[i]
  if (grepl("Br", label)) {hist_data_v2$Tissue[i] = 'Brain'}
  if (grepl("LI", label)) {hist_data_v2$Tissue[i] = 'Large Intestine'}
  if (grepl("SI", label)) {hist_data_v2$Tissue[i] = 'Small Intestine'}
  if (grepl("Stom", label)) {hist_data_v2$Tissue[i] = 'Stomach'}
  if (grepl("[Cc]ec", label)) {hist_data_v2$Tissue[i] = 'Cecum'}
}

## Update GI_Lesions
hist_data_v2$GI_Lesions = NA
for (i in 1:length(hist_data_v2$Tissue)) {
  tissue = hist_data_v2$Tissue[i]
  if (tissue != 'Brain') {
    hist_data_v2$GI_Lesions[i] = hist_data_v2$Cortex_PV_Inf[i]
    hist_data_v2$Cortex_PV_Inf[i] = NA
  }
}

hist_data_v2$Cortex_PV_Inf = as.numeric(hist_data_v2$Cortex_PV_Inf)
hist_data_v2$GI_Lesions[hist_data_v2$GI_Lesions == 'NSL'] = 0

table(hist_data_v2$GI_Lesions)

## Fix STI sample
hist_data_v2[grepl("STI", hist_data_v2$slide_label),]

## Fix STI sample
hist_data_v2$UWID[397] = NA
hist_data_v2$Virus[397] = 'Mock'
hist_data_v2$Tissue[397] = 'Kidney'
hist_data_v2$Tissue[398] = 'Kidney'

## Check for duplicates
sum(duplicated(hist_data_v2$ID[!is.na(hist_data_v2$ID)]))

## Duplicated IDs
hist_data_v2$ID[duplicated(hist_data_v2$ID) & !is.na(hist_data_v2$ID)]

cortex_cols = c('Cortex_PV_Inf','Cortex_P_Inf','Cortex_hm','Cortex_nn')
hippocampus_cols = c('Hippocampus_PV_Inf','Hippocampus_P_Inf','Hippocampus_hm','Hippocampus_nn')
thalamus_cols = c('Thalamus_Midbrain_PV_Inf','Thalamus_Midbrain_P_Inf','Thalamus_Midbrain_hm','Thalamus_Midbrain_nn')
cerebellum_cols = c('Cerebellum_PV_Inf','Cerebellum_P_Inf','Cerebellum_hm','Cerebellum_nn')
subtotal_cols = c('Cortex_Subtotal','Hippocampus_Subtotal','Thalamus_Midbrain_Subtotal','Cerebellum_Subtotal','Meningitis')

## Fill in zeros if not NA
for (i in 1:dim(hist_data_v2)[1]) {
  if (!is.na(hist_data_v2[i, 'Cortex_Subtotal']) & hist_data_v2[i, 'Cortex_Subtotal'] == 0) {
    hist_data_v2[i,cortex_cols] = c(0,0,0,0)
  }
  if (!is.na(hist_data_v2[i, 'Hippocampus_Subtotal']) & hist_data_v2[i, 'Hippocampus_Subtotal'] == 0) {
    hist_data_v2[i,hippocampus_cols] = c(0,0,0,0)
  }
  if (!is.na(hist_data_v2[i, 'Thalamus_Midbrain_Subtotal']) & hist_data_v2[i, 'Thalamus_Midbrain_Subtotal'] == 0) {
    hist_data_v2[i,thalamus_cols] = c(0,0,0,0)
  }
  if (!is.na(hist_data_v2[i, 'Cerebellum_Subtotal']) & hist_data_v2[i, 'Cerebellum_Subtotal'] == 0) {
    hist_data_v2[i,cerebellum_cols] = c(0,0,0,0)
  }
  if (!is.na(hist_data_v2[i, 'Total']) & hist_data_v2[i, 'Total'] == 0) {
    hist_data_v2[i,subtotal_cols] = c(0,0,0,0,0)
  }
}

hist_data_v2$Cortex_Subtotal = apply(hist_data_v2[,cortex_cols], 1, function(x){sum(x)})
hist_data_v2$Hippocampus_Subtotal = apply(hist_data_v2[,hippocampus_cols], 1, function(x){sum(x)})
hist_data_v2$Thalamus_Midbrain_Subtotal = apply(hist_data_v2[,thalamus_cols], 1, function(x){sum(x)})
hist_data_v2$Cerebellum_Subtotal = apply(hist_data_v2[,cerebellum_cols], 1, function(x){sum(x)})
hist_data_v2$Total = apply(hist_data_v2[,subtotal_cols], 1, function(x){sum(x)})

head(hist_data_v2, 10)

## Check for mocks with scores > 0
hist_data_v2[!is.na(hist_data_v2$Total) & hist_data_v2$Total != 0 & hist_data_v2$Virus=='Mock',]

## Get IDs duplicated in new data
dup_ids = intersect(hist_data_v2$ID, hist_prev$ID)[!is.na(intersect(hist_data_v2$ID, hist_prev$ID))]
dup_ids

hist_data_v2[hist_data_v2$ID %in% dup_ids,]

hist_prev[hist_prev$ID %in% dup_ids,]

dim(hist_prev)

hist_prev = hist_prev[!hist_prev$ID %in% dup_ids,]

dim(hist_prev)

for (cn in setdiff(colnames(hist_data_v2), colnames(hist_prev))) {
  print(cn)
  hist_prev[,cn] = NA
}

gale_hist = rbind(hist_prev[,colnames(hist_data_v2)], hist_data_v2[,colnames(hist_data_v2)])

dim(gale_hist)

write.table(gale_hist, file=file.path(hist_dir, 'Gale_Histology_21-Mar-2016_final.txt'), 
            col.names=T, row.names=F, quote=F, sep='\t', na="")
