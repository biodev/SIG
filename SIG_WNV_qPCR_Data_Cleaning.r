
source('./scripts/qpcr_data_cleaning_functions.r')

## Set data directory
data_dir = "/Users/mooneymi/Documents/SIG/WNV/qPCR"

## Read in data (byLine)
qpcr_data = read.xls(file.path(data_dir, "16-May-2016/Gale_qPCR_byLine_5-16-16 %282%29.xlsx"), sheet=1)

#qpcr_data = read.xls(file.path(data_dir, "Gale_qPCR_byLine_5-16-16_MM_animalsRemoved.xlsx"), sheet=1)

head(qpcr_data)
dim(qpcr_data)

## Replace line 82 with fixed data (special case for May 16, 2016 data)
qpcr_data = qpcr_data[qpcr_data$UW_Line != 82, ]

line_82 = read.xls(file.path(data_dir, "18-May-2016/Gale_qPCR_byLine_5-18-16.xlsx"), sheet=1)
qpcr_data = rbind(qpcr_data, line_82)

dim(qpcr_data)

## Replace lines 54 and 58 with fixed data (special case for May 16, 2016 data)
qpcr_data = qpcr_data[!qpcr_data$UW_Line %in% c(54, 58), ]

lines_54_58 = read.xls(file.path(data_dir, "23-May-2016/Gale_qPCR_byLine_5-23-16.xlsx"), sheet=1)
qpcr_data = rbind(qpcr_data, lines_54_58)

dim(qpcr_data)

## Check for duplicates
sum(duplicated(qpcr_data))

## Check that each line has a single mating
line_matings = sapply(unique(qpcr_data$UW_Line), function(x){unique(qpcr_data$Mating[qpcr_data$UW_Line==x])})
names(line_matings) = unique(qpcr_data$UW_Line)
line_matings

## Add Data_Altered and Notes columns
qpcr_data$Data_Altered = NA
qpcr_data$Notes = NA

## Annotate Virus
qpcr_data$Virus = NA
qpcr_data[grep("M",qpcr_data[,"Timepoint"]),"Virus"] <- "Mock"
qpcr_data[-grep("M",qpcr_data[,"Timepoint"]),"Virus"] <- "WNV"

## Remove 'M' from time points and convert to numeric
qpcr_data[,"Timepoint"] <- as.numeric(as.character(gsub("M","",qpcr_data[,"Timepoint"])))

## Check experiment names
sum(names(summary(qpcr_data[,"Experiment"])) == c("IFIT1","IFITM1", "IFNb1", "IL12b", "WNV")) == 5
summary(qpcr_data[,"Experiment"])

## Add Group column: UW Line, Timepoint, Virus, Tissue, Experiment separated by "_"
qpcr_data$Group = paste(qpcr_data$UW_Line, qpcr_data$Timepoint, qpcr_data$Virus, qpcr_data$Tissue, qpcr_data$Experiment, sep="_")

## Add Lab column
qpcr_data$Lab = "Gale"

head(qpcr_data)

## Read in data (byMouse)
qpcr_data_mouse = read.xls(file.path(data_dir, "16-May-2016/Gale_qPCR_byMouse_5-16-16 %281%29.xlsx"), sheet=1)

#qpcr_data_mouse = read.xls(file.path(data_dir, "Gale_qPCR_byMouse_5-16-16_MM_animalsRemoved.xlsx"), sheet=1)

head(qpcr_data_mouse)
dim(qpcr_data_mouse)

## Replace line 82 with fixed data (special case for May 16, 2016 data)
qpcr_data_mouse = qpcr_data_mouse[qpcr_data_mouse$UW_Line != 82, ]

line_82_mouse = read.xls(file.path(data_dir, "18-May-2016/Gale_qPCR_byMouse_5-18-16.xlsx"), sheet=1)
qpcr_data_mouse = rbind(qpcr_data_mouse, line_82_mouse)

dim(qpcr_data_mouse)

## Replace lines 54 and 58 with fixed data (special case for May 16, 2016 data)
qpcr_data_mouse = qpcr_data_mouse[!qpcr_data_mouse$UW_Line %in% c(54, 58), ]

lines_54_58_mouse = read.xls(file.path(data_dir, "23-May-2016/Gale_qPCR_byMouse_5-23-16.xlsx"), sheet=1)
qpcr_data_mouse = rbind(qpcr_data_mouse, lines_54_58_mouse)

dim(qpcr_data_mouse)

## Check if any CT < 15
length(which(qpcr_data_mouse[,"Ct"] < 15))

## Check if any reference CT < 15
length(which(qpcr_data_mouse[,"ref.Ct"] < 15))

## Check if an reference CT == 40
length(which(qpcr_data_mouse[,"ref.Ct"] == 40))

## Check that each line has a single mating
mouse_line_matings = sapply(unique(qpcr_data_mouse$UW_Line), function(x){unique(qpcr_data_mouse$Mating[qpcr_data_mouse$UW_Line==x])})
names(mouse_line_matings) = unique(qpcr_data_mouse$UW_Line)
mouse_line_matings

## Remove 'M' from time points and convert to numeric
qpcr_data_mouse[,"Timepoint"] <- as.numeric(as.character(gsub("M","",qpcr_data_mouse[,"Timepoint"])))

## Change Condition column name to Virus
names(qpcr_data_mouse)[which(names(qpcr_data_mouse) == "Condition")] <- "Virus"

# Check experiment names
sum(names(summary(qpcr_data_mouse[,"Experiment"])) == c("IFIT1","IFITM1", "IFNb1", "IL12b", "WNV")) == 5
summary(qpcr_data_mouse[,"Experiment"])

## Add Group column: UW Line, Timepoint, Virus, Tissue, Experiment separated by "_"
qpcr_data_mouse$Group <- paste(qpcr_data_mouse$UW_Line, qpcr_data_mouse$Timepoint, qpcr_data_mouse$Virus, 
                               qpcr_data_mouse$Tissue, qpcr_data_mouse$Experiment, sep="_")

## Check that ByLine file contains data from only those animals in ByMouse file
length(setdiff(qpcr_data$Group, qpcr_data_mouse$Group))
length(setdiff(qpcr_data_mouse$Group, qpcr_data$Group))

## Add Data_Altered and Notes columns
qpcr_data_mouse$Data_Altered = NA
qpcr_data_mouse$Notes = NA

## Add Lab column
qpcr_data_mouse$Lab = "Gale"

head(qpcr_data_mouse)

## Calculate dCt mean from byMouse data
dct_mean = aggregate(formula=qpcr_data_mouse[,"dCt"]~qpcr_data_mouse[,"Group"], data=qpcr_data_mouse, FUN=mean)
names(dct_mean) <- c("Group","dCt.mean.V2")
dct_mean[order(dct_mean[,1]),] -> dct_mean_order

## Get dCt mean from byLine data
byline_dct_mean = unique(qpcr_data[,c("Group","dCt.mean")])
byline_dct_mean[order(byline_dct_mean[,1]),] -> byline_dct_mean_order

## Check that the dCt mean calculations are the same
check_dct_mean_v2 = sapply(1:dim(dct_mean_order)[1], 
                           function(x){isTRUE(all.equal(dct_mean_order[x,2], byline_dct_mean_order[x,2]))})

## Print discrepancies, if they exist
if(sum(check_dct_mean_v2) == nrow(dct_mean_order)){
  print("All dCt mean correct")
} else {
  print("Need to clean dCt mean")
  dct_mean_errs = cbind(byline_dct_mean_order[which(check_dct_mean_v2==F),], dct_mean_order[which(check_dct_mean_v2==F),])
  dct_mean_errs
}

if (dim(dct_mean_errs)[1] > 0) {
  for (i in 1:dim(dct_mean_errs)[1]) {
    qpcr_data$dCt.mean[qpcr_data$Group==dct_mean_errs[i,3]] = dct_mean_errs[i,4]
  }
}

## Calculate the N for each group from the byMouse data
data.frame(summary(as.factor(qpcr_data_mouse[,"Group"]),maxsum=8000)) -> bymouse_n
names(bymouse_n) <- c("N.V2")
bymouse_n[,2] <- row.names(bymouse_n)
names(bymouse_n)[2] <- "Group"
bymouse_n = bymouse_n[,c(2,1)]
bymouse_n[order(bymouse_n[,"Group"]),] -> bymouse_n_order

## Get the N for each group from the byLine data 
qpcr_data[,c("Group","N")] -> byline_n
byline_n[order(byline_n[,"Group"]),] -> byline_n_order

## Print discrepancies, if they exist
if(sum(bymouse_n_order[,2] == byline_n_order[,2]) == nrow(bymouse_n_order)){
  print("All N correct")
} else {
  print("Need to clean N")
  n_errs = cbind(byline_n_order[bymouse_n_order[,2] != byline_n_order[,2],], 
                 bymouse_n_order[bymouse_n_order[,2] != byline_n_order[,2],])
  n_errs
}

if (dim(n_errs)[1] > 0) {
  for (i in 1:dim(n_errs)[1]) {
    qpcr_data$N[qpcr_data$Group==n_errs[i,3]] = n_errs[i,4]
  }
}

## Calculate dCt SD from the byMouse data 
dct_sd = aggregate(formula=qpcr_data_mouse[,"dCt"]~qpcr_data_mouse[,"Group"], data=qpcr_data_mouse, FUN=sd)
names(dct_sd) <- c("Group","dCt.sd.V2")
dct_sd[order(dct_sd[,1]),] -> dct_sd_order

## Get dCt SD from the byLine data
byline_dct_sd = qpcr_data[,c("Group","dCt.sd")]
byline_dct_sd[order(byline_dct_sd[,1]),] -> byline_dct_sd_order

check_dct_sd = sapply(1:dim(dct_sd_order)[1],function(x){isTRUE(all.equal(dct_sd_order[x,2],byline_dct_sd_order[x,2]))})

## Print discrepancies, if they exist
if(sum(check_dct_sd) == nrow(dct_sd_order)) {
  print("All dCt SD correct")
} else {
  print("Need to clean dCt SD")
  dct_sd_errs = cbind(byline_dct_sd_order[which(check_dct_sd==F),], dct_sd_order[which(check_dct_sd==F),])
  dct_sd_errs
}

if (dim(dct_sd_errs)[1] > 0) {
  for (i in 1:dim(dct_sd_errs)[1]) {
    qpcr_data$dCt.sd[qpcr_data$Group==dct_sd_errs[i,3]] = dct_sd_errs[i,4]
  }
}

## Check baseline.dct
## Add new group column to annotate baseline
qpcr_data$Group_g <- paste(qpcr_data[,"UW_Line"],qpcr_data[,"Tissue"], qpcr_data[,"Experiment"],sep="_")

## Calculate baseline, baseline is 12 for this data
qpcr_data$baseline.dCt.V2 = NA
baseline = 12

for(i in unique(qpcr_data[,"Group_g"])){
  qpcr_data[which(qpcr_data[,"Group_g"] == i),"baseline.dCt.V2"] <- 
    qpcr_data[which(qpcr_data[,"Group_g"] == i & qpcr_data[,"Virus"] == "Mock" & qpcr_data[,"Timepoint"] == baseline),"dCt.mean"]
}

# Print discrepancies, if they exist
if(sum(qpcr_data[,"baseline.dCt"] == qpcr_data[,"baseline.dCt.V2"]) == nrow(qpcr_data)){
  print("All baseline dCt correct")
} else {
  print("Need to clean baseline dCt")
  baseline_dct_errs = qpcr_data[qpcr_data[,"baseline.dCt"] != qpcr_data[,"baseline.dCt.V2"], 
                                c("Group", "Group_g", "baseline.dCt", "baseline.dCt.V2")]
  baseline_dct_errs
}

if (dim(baseline_dct_errs)[1] > 0) {
  for (i in 1:dim(baseline_dct_errs)[1]) {
    qpcr_data$baseline.dCt[qpcr_data$Group==baseline_dct_errs[i,1]] = baseline_dct_errs[i,4]
  }
}

## Calculate ddCt mean 
qpcr_data$ddCt.mean.V2 <- as.numeric(as.character(qpcr_data[,"dCt.mean"])) - as.numeric(as.character(qpcr_data[,"baseline.dCt"]))

check_ddct_mean = sapply(1:dim(qpcr_data)[1],function(x)isTRUE(all.equal(qpcr_data[x,"ddCt.mean"], qpcr_data[x,"ddCt.mean.V2"], tolerance=5.5e-8)))

## Print discrepancies, if they exist
if(sum(check_ddct_mean) == dim(qpcr_data)[1]){
  print("All ddCt mean correct")
} else {
  print("Need to clean ddCt mean")
  ddct_errs = qpcr_data[which(check_ddct_mean==F), c("Group", "Group_g", "ddCt.mean", "ddCt.mean.V2")]
  ddct_errs
}

if (dim(ddct_errs)[1] > 0) {
  for (i in 1:dim(ddct_errs)[1]) {
    qpcr_data$ddCt.mean[qpcr_data$Group==ddct_errs[i,1]] = ddct_errs[i,4]
  }
}

# check fc mean correct
qpcr_data$fc.mean.V2 <- 2^-qpcr_data[,"ddCt.mean"]

check_fc_mean = sapply(1:dim(qpcr_data)[1],function(x)isTRUE(all.equal(qpcr_data[x,"fc.mean"], qpcr_data[x,"fc.mean.V2"])))

if(sum(check_fc_mean) == dim(qpcr_data)[1]){
  print("All FC mean correct")
} else {
  print("Need to clean FC mean")
  fc_errs = qpcr_data[which(check_fc_mean==F), c("Group", "Group_g", "fc.mean", "fc.mean.V2")]
  fc_errs
}

if (dim(fc_errs)[1] > 0) {
  for (i in 1:dim(fc_errs)[1]) {
    qpcr_data$fc.mean[qpcr_data$Group==fc_errs[i,1]] = fc_errs[i,4]
  }
}

## Remove extra (unused) columns
remove_cols = c("ddCt.sd", "fc.sd", "baseline.dCt.V2", "ddCt.mean.V2", "fc.mean.V2")
qpcr_data_v2 = qpcr_data[,-as.vector(unlist(sapply(remove_cols,function(x)which(x==names(qpcr_data)))))]

## Add baseline.dCt.sd column
qpcr_data_v2$baseline.dCt.sd = NA

# compute baseline sd, use baseline = 12 saved above
for (i in unique(qpcr_data_v2[,"Group_g"])){
  qpcr_data_v2[which(qpcr_data_v2[,"Group_g"] == i),"baseline.dCt.sd"] <- 
    qpcr_data_v2[which(qpcr_data_v2[,"Group_g"] == i & 
                       qpcr_data_v2[,"Virus"] == "Mock" & 
                       qpcr_data_v2[,"Timepoint"] == baseline),"dCt.sd"]
}

## Add ddCt.se
ddCt.se = sapply(1:dim(qpcr_data_v2)[1],function(x){
  sqrt((qpcr_data_v2[x,"dCt.sd"]^2/qpcr_data_v2[x,"N"]) + 
      (qpcr_data_v2[x,"baseline.dCt.sd"]^2/qpcr_data_v2[which(qpcr_data_v2[,"Group_g"] == qpcr_data_v2[x,"Group_g"] & 
                                                              qpcr_data_v2[,"Virus"] == "Mock" & 
                                                              qpcr_data_v2[,"Timepoint"] == baseline),"N"])
  )
})
qpcr_data_v2$ddCt.se <- ddCt.se

## Order by group
qpcr_data_v3_final = qpcr_data_v2[order(qpcr_data_v2[,"Group"]),]

## Double-check all calculations are correct
## dCt mean
isTRUE(all.equal(dct_mean_order[,2], qpcr_data_v3_final[,"dCt.mean"]))

## N
isTRUE(all.equal(bymouse_n_order[,2], qpcr_data_v3_final[,"N"]))

## dCt SD
isTRUE(all.equal(dct_sd_order[,2], qpcr_data_v3_final[,"dCt.sd"]))

## baseline dCt
isTRUE(all.equal(qpcr_data[order(qpcr_data[,"Group"]),"baseline.dCt.V2"],qpcr_data_v3_final[,"baseline.dCt"]))

## ddCt mean
isTRUE(all.equal(qpcr_data[order(qpcr_data[,"Group"]),"ddCt.mean.V2"],qpcr_data_v3_final[,"ddCt.mean"]))

## FC mean
isTRUE(all.equal(qpcr_data[order(qpcr_data[,"Group"]),"fc.mean.V2"],qpcr_data_v3_final[,"fc.mean"]))

## baseline SD, computed using corrected dCt SD at timepoint 12

## ddCt SE, computed using correct dCt SD and baseline SD

## Remove Group_g column
qpcr_data_final_format <- qpcr_data_v3_final[,!(names(qpcr_data_v3_final) %in% "Group_g")]

## Order columns according to data dictionary (byLine)
## Note: you may have to change the path to the data dictionary
data_dict <- read.xls(xls="./data/WNV_Data_Dictionary.xlsx", sheet="qPCR Data - By Line", as.is=T)
qpcr_data_final_format_order = qpcr_data_final_format[, data_dict[,1]]

## Order columns according to data dictionary (byMouse)
data_dict <- read.xls(xls="./data/WNV_Data_Dictionary.xlsx", sheet="qPCR Data - By Mouse", as.is=T)
qpcr_data_mouse_order = qpcr_data_mouse[, data_dict[,1]]

## Check both byLine and byMouse have the same UW lines
names(summary(as.factor(qpcr_data_mouse_order[,"UW_Line"]))) == 
names(summary(as.factor(qpcr_data_final_format_order[,"UW_Line"])))

dim(qpcr_data_final_format_order)
dim(qpcr_data_mouse_order)

## Read in previous version of byLine data
cleaned_data_dir = "~/Documents/SIG/WNV/Cleaned_Data_Releases/23-Mar-2016/"
prev_qpcr_byline = read.xls(xls=file.path(cleaned_data_dir, "Gale_qPCR_byLine_23-Mar-2016_final.xlsx"), sheet=1, as.is=T)
dim(prev_qpcr_byline)

## Check for duplicates (new data will overwrite old)
dup_groups = intersect(qpcr_data_final_format_order$Group, prev_qpcr_byline$Group)
length(dup_groups)

## Remove duplicated groups from previous data
prev_qpcr_byline = prev_qpcr_byline[!prev_qpcr_byline$Group %in% dup_groups, ]
dim(prev_qpcr_byline)

## Combine old and new data
qpcr_byline_updated = rbind(prev_qpcr_byline, qpcr_data_final_format_order)
dim(qpcr_byline_updated)

# Set blanks to NA
qpcr_byline_updated_cleaned = clean_na(qpcr_byline_updated)

# Remove duplicates
if(sum(duplicated(qpcr_byline_updated_cleaned)) != 0){
  qpcr_byline_updated_cleaned = qpcr_byline_updated_cleaned[!duplicated(qpcr_byline_updated_cleaned),]
}
dim(qpcr_byline_updated_cleaned)

## Read in previous version of byMouse data
prev_qpcr_bymouse = read.xls(xls=file.path(cleaned_data_dir, "Gale_qPCR_byMouse_23-Mar-2016_final.xlsx"), sheet=1)
dim(prev_qpcr_bymouse)

## Check for duplicates (new data will overwrite old)
prev_ids = paste(prev_qpcr_bymouse$ID, prev_qpcr_bymouse$Tissue, prev_qpcr_bymouse$Experiment, sep='_')
new_ids = paste(qpcr_data_mouse_order$ID, qpcr_data_mouse_order$Tissue, qpcr_data_mouse_order$Experiment, sep='_')
dup_ids = intersect(prev_ids, new_ids)
length(dup_ids)

## Remove duplicated IDs from previous data
prev_qpcr_bymouse = prev_qpcr_bymouse[which(!prev_ids %in% dup_ids), ]
dim(prev_qpcr_bymouse)

## Combine old and new data
qpcr_bymouse_updated = rbind(prev_qpcr_bymouse, qpcr_data_mouse_order)
dim(qpcr_bymouse_updated)

## Set blanks to NA
qpcr_bymouse_updated_cleaned = clean_na(qpcr_bymouse_updated)

## Remove duplicates
if(sum(duplicated(qpcr_bymouse_updated_cleaned)) != 0){
  qpcr_bymouse_updated_cleaned = qpcr_byline_updated_cleaned[!duplicated(qpcr_bymouse_updated_cleaned),]
}
dim(qpcr_bymouse_updated_cleaned)



## Save ByLine Data
write.table(qpcr_data_final_format_order, file=file.path(data_dir, "23-May-2016/Gale_qPCR_byLine_5-23-16_MM_updated.txt"), 
            col.names=T, row.names=F, sep='\t', quote=F, na="")
## Save ByMouse Data
write.table(qpcr_data_mouse_order, file=file.path(data_dir, "23-May-2016/Gale_qPCR_byMouse_5-23-16_MM_updated.txt"), 
            col.names=T, row.names=F, sep='\t', quote=F, na="")
