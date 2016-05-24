library(gdata)

describe = function(obj) {
  if ('help' %in% names(attributes(obj))) {
    writeLines(attr(obj, 'help'))
  }
}
attr(describe, 'help') = "
This function prints the contents of the 'help' attribute of any R object. 
It is meant to provide help documentation in the same vein as Docstrings in Python. 
"

read_flow_exp_file = function(f, cn_expected) {
  ### Parse the Treg panel - first sheet of the excel workbook
  ## Read the file
  dat1 = read.xls(f, sheet=1, verbose=F, header=F, as.is=T, na.strings=c("NA", "", " ", "#DIV/0!"))
  print(sheetNames(f)[1])
  
  ## Identify the first row containing information
  start_row = min(which(!is.na(dat1[,1])))
  ## Get the column names from the first row
  cn1 = dat1[start_row,]
  ## Subset the dataframe and list of column names to remove empty columns
  dat1 = dat1[,!is.na(cn1)]
  cn1 = cn1[!is.na(cn1)]
  
  ## DEBUGGING / TESTING
  #print(cn1)
  
  ## Parse nested column names
  cn_start = which(!is.na(dat1[1,]))
  lenTreg = cn_start[2]-cn_start[1]-1
  cn1[cn_start[1]:(cn_start[1]+lenTreg-1)] = paste(rep('CD4pos_Foxp3neg_', lenTreg), cn1[cn_start[1]:(cn_start[1]+lenTreg-1)], sep='')
  cn1[cn_start[2]:(cn_start[2]+lenTreg-1)] = paste(rep('Tregs_', lenTreg), cn1[cn_start[2]:(cn_start[2]+lenTreg-1)], sep='')
  
  ## Identify the time column
  start_col = which(cn1=='time')
  ## Standardize the column names
  cn1 = fix_column_names(cn1, 'treg', start_col)
  ## Identify the last row with data
  end_row = which(is.na(dat1[,1]))[2]
  ## Update the column names and subset the dataframe
  colnames(dat1) = cn1
  dat1 = dat1[(start_row+1):(end_row-1),]
  
  ## DEBUGGING / TESTING
  #print(end_row)
  #print(dim(dat1))
  
  ## Parse the T-cell panel - second and third sheets in the excel workbook
  ## Initialize counter to keep track of what sheet is being parsed
  i = 1
  ## Identify the T-cell sheets (since both D7 and D21 panels are not always included)
  CD8_sheets = grep("CD8", sheetNames(f))
  print(sheetNames(f)[CD8_sheets])
  ## Iterate through the sheets containing the T-cell panels
  for (sheet in CD8_sheets) {
    if (i == 1) {
      sheet_name = sheetNames(f)[sheet]
      ## Read the file, identify the first row containing information, and get column names
      dat2 = read.xls(f, sheet=sheet, verbose=F, header=F, as.is=T, na.strings=c("NA", "", " ", "#DIV/0!"))
	  
	  ## DEBUGGING / TESTING
      #print(dim(dat2))
	  
      start_row = min(which(dat2[,1]=='ID' | dat2[,1]=='UNC strain'))
      cn2 = dat2[start_row,]
      dat2 = dat2[,!is.na(cn2)]
      cn2 = cn2[!is.na(cn2)]
      
      ## DEBUGGING / TESTING
      #print(cn2)
      
      ## Parse nested column names
      cn_start = which(!is.na(dat2[1,]))
      lenCD8 = cn_start[2]-cn_start[1]-1
      lenCD4 = length(cn2)-cn_start[2]+1
      cn2[cn_start[1]:(cn_start[1]+lenCD8-1)] = paste(rep('CD8pos_', lenCD8), cn2[cn_start[1]:(cn_start[1]+lenCD8-1)], sep='')
      cn2[cn_start[2]:(cn_start[2]+lenCD4-1)] = paste(rep('CD4pos_', lenCD4), cn2[cn_start[2]:(cn_start[2]+lenCD4-1)], sep='')
      
      ## Identify the time column
      start_col = which(cn2=='time')
      ## Determine which T-cell panel is being parsed (D7 or D21) and standardize the column names
      if (length(grep("d7", sheet_name)) > 0 | length(grep("d12", sheet_name)) > 0) {
        cn2 = fix_column_names(cn2, 'tcell_d7', start_col)
      } else {
        cn2 = fix_column_names(cn2, 'tcell_d21', start_col)
      }
      ## Identify the last row with data
      end_row = which(is.na(dat2[,1]))[2]
      ## Update the column names and subset the dataframe
      colnames(dat2) = cn2
      dat2 = dat2[(start_row+1):(end_row-1),]
    } else {
      ## Repeat the above process for the second T-cell panel if it exists in the file
      sheet_name = sheetNames(f)[sheet]
      dat2b = read.xls(f, sheet=sheet, verbose=F, header=F, as.is=T, na.strings=c("NA", "", " ", "#DIV/0!"))
      start_row = min(which(dat2b[,1]=='ID' | dat2b[,1]=='UNC strain'))
      cn2 = dat2b[start_row,]
      dat2b = dat2b[,!is.na(cn2)]
      cn2 = cn2[!is.na(cn2)]
      
      ## DEBUGGING / TESTING
      #print(cn2)
      
      cn_start = which(!is.na(dat2b[1,]))
      lenCD8 = cn_start[2]-cn_start[1]-1
      lenCD4 = length(cn2)-cn_start[2]+1
      cn2[cn_start[1]:(cn_start[1]+lenCD8-1)] = paste(rep('CD8pos_', lenCD8), cn2[cn_start[1]:(cn_start[1]+lenCD8-1)], sep='')
      cn2[cn_start[2]:(cn_start[2]+lenCD4-1)] = paste(rep('CD4pos_', lenCD4), cn2[cn_start[2]:(cn_start[2]+lenCD4-1)], sep='')
      start_col = which(cn2=='time')
      if (length(grep("d7", sheet_name)) > 0 | length(grep("d12", sheet_name)) > 0) {
        cn2 = fix_column_names(cn2, 'tcell_d7', start_col)
      } else {
        cn2 = fix_column_names(cn2, 'tcell_d21', start_col)
      }
      end_row = which(is.na(dat2b[,1]))[2]
      colnames(dat2b) = cn2
      dat2b = dat2b[(start_row+1):(end_row-1),]
      
      ## Merge the two T-cell panels, since only some animals (timepoints) will be included in each
      dat2 = merge(dat2, dat2b, by=intersect(colnames(dat2), colnames(dat2b)), all=T)
    }
    i = i + 1
  }
  
  ## DEBUGGING / TESTING
  #print(end_row)
  #print(dim(dat2))
  
  ## Parse the ICS Panel - the third (or fourth) sheet in the excel workbook
  ## Read the file, identify the first row containing information, and get column names
  dat3 = read.xls(f, sheet=(max(CD8_sheets)+1), verbose=F, header=F, as.is=T, na.strings=c("NA", "", " ", "#DIV/0!"))
  print(sheetNames(f)[(max(CD8_sheets)+1)])
  start_row = min(which(dat3[,1]=='ID' | dat3[,1]=='UNC strain'))
  cn3 = dat3[start_row,]
  dat3 = dat3[,!is.na(cn3)]
  cn3 = cn3[!is.na(cn3)]
  
  ## Parse the nested column names
  cn_start1 = which(!is.na(dat3[1,]))
  cn_start2 = which(!is.na(dat3[2,]))
  lenICSCD8_1 = cn_start2[2]-cn_start2[1]-1
  lenICSCD4_1 = cn_start2[3]-cn_start2[2]-6
  lenICSCD8_2 = cn_start2[4]-cn_start2[3]-1
  lenICSCD4_2 = cn_start2[5]-cn_start2[4]-6
  lenICSCD8_3 = cn_start2[6]-cn_start2[5]-1
  lenICSCD4_3 = cn_start2[7]-cn_start2[6]-6
  lenICSCD8_4 = cn_start2[8]-cn_start2[7]-1
  lenICSCD4_4 = length(cn3)-cn_start2[8]+1
  cn3[cn_start2[1]:(cn_start2[1]+lenICSCD8_1-1)] = paste(rep('CD8pos_', lenICSCD8_1), cn3[cn_start2[1]:(cn_start2[1]+lenICSCD8_1-1)], sep='')
  cn3[cn_start2[2]:(cn_start2[2]+lenICSCD4_1-1)] = paste(rep('CD4pos_', lenICSCD4_1), cn3[cn_start2[2]:(cn_start2[2]+lenICSCD4_1-1)], sep='')
  cn3[cn_start2[3]:(cn_start2[3]+lenICSCD8_2-1)] = paste(rep('CD8_', lenICSCD8_2), cn3[cn_start2[3]:(cn_start2[3]+lenICSCD8_2-1)], sep='')
  cn3[cn_start2[4]:(cn_start2[4]+lenICSCD4_2-1)] = paste(rep('CD4_', lenICSCD4_2), cn3[cn_start2[4]:(cn_start2[4]+lenICSCD4_2-1)], sep='')
  cn3[cn_start2[5]:(cn_start2[5]+lenICSCD8_3-1)] = paste(rep('CD8_', lenICSCD8_3), cn3[cn_start2[5]:(cn_start2[5]+lenICSCD8_3-1)], sep='')
  cn3[cn_start2[6]:(cn_start2[6]+lenICSCD4_3-1)] = paste(rep('CD4_', lenICSCD4_3), cn3[cn_start2[6]:(cn_start2[6]+lenICSCD4_3-1)], sep='')
  cn3[cn_start2[7]:(cn_start2[7]+lenICSCD8_4-1)] = paste(rep('CD8_', lenICSCD8_4), cn3[cn_start2[7]:(cn_start2[7]+lenICSCD8_4-1)], sep='')
  cn3[cn_start2[8]:(cn_start2[8]+lenICSCD4_4-1)] = paste(rep('CD4_', lenICSCD4_4), cn3[cn_start2[8]:(cn_start2[8]+lenICSCD4_4-1)], sep='')
  
  lenDMSO = cn_start1[2]-cn_start1[1]
  lenNS4B = cn_start1[3]-cn_start1[2]
  lenHIWNV = cn_start1[4]-cn_start1[3]
  lenCD3CD28 = length(cn3)-cn_start1[4]+1
  cn3[cn_start1[1]:(cn_start1[1]+lenDMSO-1)] = paste(rep('DMSO_', lenDMSO), cn3[cn_start1[1]:(cn_start1[1]+lenDMSO-1)], sep='')
  cn3[cn_start1[2]:(cn_start1[2]+lenNS4B-1)] = paste(rep('NS4B_', lenNS4B), cn3[cn_start1[2]:(cn_start1[2]+lenNS4B-1)], sep='')
  cn3[cn_start1[3]:(cn_start1[3]+lenHIWNV-1)] = paste(rep('HIWNV_', lenHIWNV), cn3[cn_start1[3]:(cn_start1[3]+lenHIWNV-1)], sep='')
  cn3[cn_start1[4]:(cn_start1[4]+lenCD3CD28-1)] = paste(rep('CD3CD28_', lenCD3CD28), cn3[cn_start1[4]:(cn_start1[4]+lenCD3CD28-1)], sep='')
  
  ## Identify the first time column (DMSO ICS panel)
  start_col = which(cn3=='DMSO_time')
  ## Standardize the column names
  cn3 = fix_column_names(cn3, 'ics', start_col)
  ## Identify the last column with data
  end_row = which(is.na(dat3[,1]))[3]
  ## Update the column names and subset the dataframe
  colnames(dat3) = cn3
  dat3 = dat3[(start_row+1):(end_row-1),]
  
  ## DEBUGGING / TESTING
  #print(end_row)
  #print(dim(dat3))
  #print(colnames(dat1))
  #print(colnames(dat2))
  #print(colnames(dat3))
  
  ## Merge all the panels
  dat = merge(dat1, dat2 , by=c('UNC_strain', 'UW_strain', 'RIX_ID', 'Timepoint', 'Tissue', 'Total_Cell_Count'), suffixes=c('', ''), all=T)
  dat = merge(dat, dat3, by=c('UNC_strain', 'UW_strain', 'RIX_ID', 'Timepoint', 'Tissue', 'Total_Cell_Count'), suffixes=c('', ''), all=T)
  
  ## Fix final column names
  cn_final = colnames(dat)
  cn_final = gsub(" ", "_", cn_final)
  cn_final = gsub("UNC_strain", "Mating", cn_final)
  cn_final = gsub("UW_strain", "UW_Line", cn_final)
  colnames(dat) = cn_final
  
  ## Update the ID column
  ids = paste(dat$Mating, dat$RIX_ID, sep='_')
  dat$ID = ids
  
  ## Fix timepoints
  #dat$Timepoint[dat$Timepoint=='d7'] = '7'
  #dat$Timepoint[dat$Timepoint=='d12'] = '12'
  #dat$Timepoint[dat$Timepoint=='d21'] = '21'
  #dat$Timepoint[dat$Timepoint=='d28'] = '28'
  #dat$Timepoint[dat$Timepoint=='d12m'] = '12m'
  #dat$Timepoint[dat$Timepoint=='d28m'] = '28m'
  dat$Timepoint = gsub('d', '', dat$Timepoint)
  
  ## Create Virus column
  dat$Virus = 'WNV'
  #dat$Virus[dat$Timepoint=='12m'] = 'Mock'
  #dat$Virus[dat$Timepoint=='28m'] = 'Mock'
  dat$Virus[grepl('m', dat$Timepoint)] = 'Mock'
  
  #dat$Timepoint[dat$Timepoint=='12m'] = '12'
  #dat$Timepoint[dat$Timepoint=='28m'] = '28'
  dat$Timepoint = gsub('m', '', dat$Timepoint)
  
  ## Print unexpected column names, and add columns with NAs if necessary
  if (length(setdiff(colnames(dat), cn_expected)) > 0) {
    print("Unexpected columns: ")
    print(setdiff(colnames(dat), cn_expected))
  }
  for (col in cn_expected) {
    if (!(col %in% colnames(dat))) {
      dat[,col] = NA
    }
  }
  return(dat)
}
attr(read_flow_exp_file, 'help') = "
This function parses flow cytometry data from an Excel workbook.

Parameters:
f: The Excel file name.
cn_expected: A character vector containing the expected column names.

Returns:
A dataframe containing the processed flow cytometry data.
"

fix_column_names = function(cn, panel, start_col=8) {
  ## DEBUGGING / TESTING
  #print(start_col)
  #print(cn)
  
  ## Standardize the column names for all panels
  cn = trim(gsub("\\+([^[:space:]])", "+ \\1", cn))
  cn = gsub("\\+", "pos", cn)
  cn = gsub("IL-17", "IL_17", cn)
  cn = gsub("CTLA-4", "CTLA_4", cn)
  cn = gsub("PSGL-1", "PSGL_1", cn)
  cn = trim(gsub("-([^[:space:]])", "- \\1", cn))
  cn = gsub(" ", "_", cn)
  cn = gsub("-", "neg", cn)
  cn = gsub("cell", "Cell", cn)
  cn = gsub("count", "Count", cn)
  cn = gsub("live", "Live", cn)
  cn = gsub("singlets", "Singlets", cn)
  cn = gsub("lymphocytes", "Lymphocytes", cn)
  cn = gsub("total", "Total", cn)
  cn = gsub("time", "Time", cn)
  cn = gsub("Organ", "Tissue", cn)
  cn = gsub("Mouse_#", "RIX_ID", cn)
  cn = gsub("NS4B_CD4pos_Tbetpos", "NS4B_CD4_Tbetpos", cn)
  cn = gsub("HIWNV_CD4pos_Tbetpos", "HIWNV_CD4_Tbetpos", cn)
  cn = gsub("CD3CD28_CD4pos_Tbetpos", "CD3CD28_CD4_Tbetpos", cn)
  cn = gsub("%_CTLA_4pos", "CTLA_4pos", cn)
  cn = gsub("NS4b", "NS4B", cn)
  cn = gsub("TNFa", "TNFA", cn)
  cn = gsub("IFNg", "IFNG", cn)
  
  panel = paste(panel, '_', sep='')
  cn[start_col:length(cn)] = gsub("^", panel, cn[start_col:length(cn)])
  
  ## DEBUGGING / TESTING
  #print(cn)
  
  return(cn)
}
attr(fix_column_names, 'help') = "
This function standardizes the column names of each flow cytometry panel.

Parameters:
cn: A character vector containing the original column names.
panel: A string which will be appended to column names indicating the panel being processed (e.g. 'treg').
start_column: A number indicating in which column the flow data starts (preceding columns are sample identifiers) 

Returns:
A character vector containing the fixed column names.
"

calc_treg_counts = function(flow_df) {
	## Treg parent cell populations
	flow_df$treg_Live_count = flow_df$Total_Cell_Count*(flow_df$treg_Time/100)*(flow_df$treg_Singlets/100)*(flow_df$treg_Live/100)
	flow_df$treg_Lymphocytes_count = flow_df$treg_Live_count*(flow_df$treg_Lymphocytes/100)
	flow_df$treg_Total_CD4pos_count = flow_df$treg_Lymphocytes_count*(flow_df$treg_Total_CD4pos/100)
	flow_df$treg_CD4pos_Foxp3neg_count = flow_df$treg_Total_CD4pos_count*(flow_df$treg_CD4pos_Foxp3neg/100)
	flow_df$treg_T_regs_count = flow_df$treg_Total_CD4pos_count*(flow_df$treg_T_regs/100)
	
	## Treg sub-populations
	CD4posFoxp3neg_subpops = c('treg_CD4pos_Foxp3neg_CCR5pos','treg_CD4pos_Foxp3neg_CD25pos','treg_CD4pos_Foxp3neg_CD29pos','treg_CD4pos_Foxp3neg_CD44hi','treg_CD4pos_Foxp3neg_CD73pos','treg_CD4pos_Foxp3neg_CTLA_4pos','treg_CD4pos_Foxp3neg_CXCR3pos','treg_CD4pos_Foxp3neg_GITRpos','treg_CD4pos_Foxp3neg_ICOSpos')
	T_reg_subpops = c('treg_Tregs_CCR5pos','treg_Tregs_CD25pos','treg_Tregs_CD29pos','treg_Tregs_CD44hi','treg_Tregs_CD73pos','treg_Tregs_CTLA_4pos','treg_Tregs_CXCR3pos','treg_Tregs_GITRpos','treg_Tregs_ICOSpos')
	
	for (col in CD4posFoxp3neg_subpops) {
		new_col = paste(col, '_count', sep='')
		flow_df[,new_col] = flow_df$treg_CD4pos_Foxp3neg_count*(flow_df[,col]/100)
	}
	for (col in T_reg_subpops) {
		new_col = paste(col, '_count', sep='')
		flow_df[,new_col] = flow_df$treg_T_regs_count*(flow_df[,col]/100)
	}
	
	return(flow_df)
}
attr(calc_treg_counts, 'help') = "
This function calculates absolute cell counts for the Treg panel.

Parameters:
flow_df: The dataframe containing the flow cytometry data (total cell count, and percentages for all cell populations).

Returns:
A dataframe with columns for the cell counts added.
"

calc_tcell_counts = function(flow_df) {
	## T cell parent cell populations
	flow_df$tcell_d7_Live_count = flow_df$Total_Cell_Count*(flow_df$tcell_d7_Time/100)*(flow_df$tcell_d7_Singlets/100)*(flow_df$tcell_d7_Live/100)
	flow_df$tcell_d7_Lymphocytes_count = flow_df$tcell_d7_Live_count*(flow_df$tcell_d7_Lymphocytes/100)
	flow_df$tcell_d7_CD3_count = flow_df$tcell_d7_Lymphocytes_count*(flow_df$tcell_d7_CD3/100)
	flow_df$tcell_d7_CD8_count = flow_df$tcell_d7_CD3_count*(flow_df$tcell_d7_CD8/100)
	flow_df$tcell_d7_CD4_count = flow_df$tcell_d7_CD3_count*(flow_df$tcell_d7_CD4/100)
	
	flow_df$tcell_d21_Live_count = flow_df$Total_Cell_Count*(flow_df$tcell_d21_Time/100)*(flow_df$tcell_d21_Singlets/100)*(flow_df$tcell_d21_Live/100)
	flow_df$tcell_d21_Lymphocytes_count = flow_df$tcell_d21_Live_count*(flow_df$tcell_d21_Lymphocytes/100)
	flow_df$tcell_d21_CD3_count = flow_df$tcell_d21_Lymphocytes_count*(flow_df$tcell_d21_CD3/100)
	flow_df$tcell_d21_CD8_count = flow_df$tcell_d21_CD3_count*(flow_df$tcell_d21_CD8/100)
	flow_df$tcell_d21_CD4_count = flow_df$tcell_d21_CD3_count*(flow_df$tcell_d21_CD4/100)

	## T cell sup-populations
	d7_CD8pos_subpops = c("tcell_d7_CD8pos_CCR5pos", "tcell_d7_CD8pos_CCR7pos", "tcell_d7_CD8pos_CCR7neg_CD62Lneg", 
	"tcell_d7_CD8pos_CCR7neg_CD62Lpos", "tcell_d7_CD8pos_CCR7pos_CD62Lneg", 
	"tcell_d7_CD8pos_CCR7pos_CD62Lpos", "tcell_d7_CD8pos_CD44pos", 
	"tcell_d7_CD8pos_CD62Lneg", "tcell_d7_CD8pos_CXCR3pos", "tcell_d7_CD8pos_Ki67pos", 
	"tcell_d7_CD8pos_NS4Bpos", "tcell_d7_CD8pos_SLECs", "tcell_d7_CD8pos_MPECs", 
	"tcell_d7_CD8pos_CD44hi", "tcell_CD8pos_Trm")
	
	d7_CD4pos_subpops = c("tcell_d7_CD4pos_CCR5pos", "tcell_d7_CD4pos_CCR7pos", "tcell_d7_CD4pos_CCR7pos_CD62Lpos", 
	"tcell_d7_CD4pos_CCR7pos_CD62Lneg", "tcell_d7_CD4pos_CCR7neg_CD62Lpos", 
	"tcell_d7_CD4pos_CCR7neg_CD62Lneg", "tcell_d7_CD4pos_CD44pos", 
	"tcell_d7_CD4pos_CD62Lneg", "tcell_d7_CD4pos_CXCR3pos", "tcell_d7_CD4pos_Ki67pos", 
	"tcell_d7_CD4pos_CD44hi")
	
	d7_NS4Bpos_subpops = c("tcell_d7_CD8pos_NS4Bpos_MPECs", "tcell_d7_CD8pos_NS4Bpos_SLECs", "tcell_CD8pos_NS4Bpos_Trm")

	d21_CD8pos_subpops = c("tcell_d21_CD8pos_CCR5pos", "tcell_d21_CD8pos_CCR7pos", "tcell_d21_CD8pos_CCR7neg_CD62Lneg", 
	"tcell_d21_CD8pos_CCR7neg_CD62Lpos", "tcell_d21_CD8pos_CCR7pos_CD62Lneg", 
	"tcell_d21_CD8pos_CCR7pos_CD62Lpos", "tcell_d21_CD8pos_CD44pos", 
	"tcell_d21_CD8pos_CD62Lneg", "tcell_d21_CD8pos_CXCR3pos", "tcell_d21_CD8pos_Ki67pos", 
	"tcell_d21_CD8pos_NS4Bpos", "tcell_d21_CD8pos_CD69pos", "tcell_d21_CD8pos_CD103pos", 
	"tcell_d21_CD8pos_CD69neg_CD103pos", "tcell_d21_CD8pos_CD69pos_CD103pos", 
	"tcell_d21_CD8pos_CD69pos_CD103neg", "tcell_d21_CD8pos_CD69neg_CD103neg")
	
	d21_CD4pos_subpops = c("tcell_d21_CD4pos_CCR5pos", "tcell_d21_CD4pos_CCR7pos", "tcell_d21_CD4pos_CCR7pos_CD62Lpos", 
	"tcell_d21_CD4pos_CCR7pos_CD62Lneg", "tcell_d21_CD4pos_CCR7neg_CD62Lpos", 
	"tcell_d21_CD4pos_CCR7neg_CD62Lneg", "tcell_d21_CD4pos_CD44pos", 
	"tcell_d21_CD4pos_CD62Lneg", "tcell_d21_CD4pos_CXCR3pos", "tcell_d21_CD4pos_Ki67pos", 
	"tcell_d21_CD4pos_CD69pos", "tcell_d21_CD4pos_CD103pos", "tcell_d21_CD4pos_CD69neg_CD103pos", 
	"tcell_d21_CD4pos_CD69pos_CD103pos", "tcell_d21_CD4pos_CD69pos_CD103neg", 
	"tcell_d21_CD4pos_CD69neg_CD103neg")
	
	d21_NS4Bpos_subpops = c("tcell_d21_CD8pos_NS4Bpos_CD103pos")

	for (col in d7_CD8pos_subpops) {
		new_col = paste(col, '_count', sep='')
		flow_df[,new_col] = flow_df$tcell_d7_CD8_count*(flow_df[,col]/100)
	}
	for (col in d7_CD4pos_subpops) {
		new_col = paste(col, '_count', sep='')
		flow_df[,new_col] = flow_df$tcell_d7_CD4_count*(flow_df[,col]/100)
	}
	for (col in d7_NS4Bpos_subpops) {
		new_col = paste(col, '_count', sep='')
		flow_df[,new_col] = flow_df$tcell_d7_CD8pos_NS4Bpos_count*(flow_df[,col]/100)
	}
	
	for (col in d21_CD8pos_subpops) {
		new_col = paste(col, '_count', sep='')
		flow_df[,new_col] = flow_df$tcell_d21_CD8_count*(flow_df[,col]/100)
	}
	for (col in d21_CD4pos_subpops) {
		new_col = paste(col, '_count', sep='')
		flow_df[,new_col] = flow_df$tcell_d21_CD4_count*(flow_df[,col]/100)
	}
	for (col in d21_NS4Bpos_subpops) {
		new_col = paste(col, '_count', sep='')
		flow_df[,new_col] = flow_df$tcell_d21_CD8pos_NS4Bpos_count*(flow_df[,col]/100)
	}
	
	return(flow_df)
}
attr(calc_tcell_counts, 'help') = "
This function calculates absolute cell counts for the T-Cell panel.

Parameters:
flow_df: The dataframe containing the flow cytometry data (total cell count, and percentages for all cell populations).

Returns:
A dataframe with columns for the cell counts added.
"

calc_ics_counts = function(flow_df) {
	## ICS parent cell populations
	flow_df$ics_DMSO_Live_count = flow_df$Total_Cell_Count*(flow_df$ics_DMSO_Time/100)*(flow_df$ics_DMSO_Singlets/100)*(flow_df$ics_DMSO_Live/100)
	flow_df$ics_NS4B_Live_count = flow_df$Total_Cell_Count*(flow_df$ics_NS4B_Time/100)*(flow_df$ics_NS4B_Singlets/100)*(flow_df$ics_NS4B_Live/100)
	flow_df$ics_HIWNV_Live_count = flow_df$Total_Cell_Count*(flow_df$ics_HIWNV_Time/100)*(flow_df$ics_HIWNV_Singlets/100)*(flow_df$ics_HIWNV_Live/100)
	flow_df$ics_CD3CD28_Live_count = flow_df$Total_Cell_Count*(flow_df$ics_CD3CD28_Time/100)*(flow_df$ics_CD3CD28_Singlets/100)*(flow_df$ics_CD3CD28_Live/100)
	
	flow_df$ics_DMSO_Lymphocytes_count = flow_df$ics_DMSO_Live_count*(flow_df$ics_DMSO_Lymphocytes/100)
	flow_df$ics_NS4B_Lymphocytes_count = flow_df$ics_NS4B_Live_count*(flow_df$ics_NS4B_Lymphocytes/100)
	flow_df$ics_HIWNV_Lymphocytes_count = flow_df$ics_HIWNV_Live_count*(flow_df$ics_HIWNV_Lymphocytes/100)
	flow_df$ics_CD3CD28_Lymphocytes_count = flow_df$ics_CD3CD28_Live_count*(flow_df$ics_CD3CD28_Lymphocytes/100)
	
	flow_df$ics_DMSO_CD3_count = flow_df$ics_DMSO_Lymphocytes_count*(flow_df$ics_DMSO_CD3/100)
	flow_df$ics_NS4B_CD3_count = flow_df$ics_NS4B_Lymphocytes_count*(flow_df$ics_NS4B_CD3/100)
	flow_df$ics_HIWNV_CD3_count = flow_df$ics_HIWNV_Lymphocytes_count*(flow_df$ics_HIWNV_CD3/100)
	flow_df$ics_CD3CD28_CD3_count = flow_df$ics_CD3CD28_Lymphocytes_count*(flow_df$ics_CD3CD28_CD3/100)
	
	flow_df$ics_DMSO_CD8_count = flow_df$ics_DMSO_CD3_count*(flow_df$ics_DMSO_CD8/100)
	flow_df$ics_NS4B_CD8_count = flow_df$ics_NS4B_CD3_count*(flow_df$ics_NS4B_CD8/100)
	flow_df$ics_HIWNV_CD8_count = flow_df$ics_HIWNV_CD3_count*(flow_df$ics_HIWNV_CD8/100)
	flow_df$ics_CD3CD28_CD8_count = flow_df$ics_CD3CD28_CD3_count*(flow_df$ics_CD3CD28_CD8/100)
	
	flow_df$ics_DMSO_CD4_count = flow_df$ics_DMSO_CD3_count*(flow_df$ics_DMSO_CD4/100)
	flow_df$ics_NS4B_CD4_count = flow_df$ics_NS4B_CD3_count*(flow_df$ics_NS4B_CD4/100)
	flow_df$ics_HIWNV_CD4_count = flow_df$ics_HIWNV_CD3_count*(flow_df$ics_HIWNV_CD4/100)
	flow_df$ics_CD3CD28_CD4_count = flow_df$ics_CD3CD28_CD3_count*(flow_df$ics_CD3CD28_CD4/100)
	
	## ICS sub-populations
	ICS_DMSO_CD8pos_subpops = c("ics_DMSO_CD8pos_CD25pos", "ics_DMSO_CD8pos_CD44pos", "ics_DMSO_CD8pos_CD44neg_CD62Lpos", "ics_DMSO_CD8pos_CD44neg_CD62Lneg", 
	"ics_DMSO_CD8pos_CD44pos_CD62Lneg", "ics_DMSO_CD8pos_CD44pos_CD62Lpos", "ics_DMSO_CD8pos_CD62Lneg", "ics_DMSO_CD8pos_IFNGpos", 
	"ics_DMSO_CD8pos_IL_17pos", "ics_DMSO_CD8pos_Tbetpos", "ics_DMSO_CD8pos_TNFApos", "ics_DMSO_CD8pos_TNFAneg_IFNGpos", "ics_DMSO_CD8pos_TNFApos_IFNGneg", 
	"ics_DMSO_CD8pos_TNFApos_IFNGpos")
	
	ICS_DMSO_CD4pos_subpops = c("ics_DMSO_CD4pos_CD25pos", "ics_DMSO_CD4pos_CD44pos", "ics_DMSO_CD4pos_CD44neg_CD62Lpos", "ics_DMSO_CD4pos_CD44neg_CD62Lneg", 
	"ics_DMSO_CD4pos_CD44pos_CD62Lneg", "ics_DMSO_CD4pos_CD44pos_CD62Lpos", "ics_DMSO_CD4pos_CD62Lneg", "ics_DMSO_CD4pos_IFNGpos", 
	"ics_DMSO_CD4pos_IL_17pos", "ics_DMSO_CD4pos_Ly6Cpos", "ics_DMSO_CD4pos_Ly6Cneg_PSGL_1neg", "ics_DMSO_CD4pos_Ly6Cneg_PSGL_1pos", 
	"ics_DMSO_CD4pos_Ly6Cpos_PSGL_1neg", "ics_DMSO_CD4pos_Ly6Cpos_PSGL_1pos", "ics_DMSO_CD4pos_PSGL_1pos", "ics_DMSO_CD4pos_Tbetpos", 
	"ics_DMSO_CD4pos_TNFApos", "ics_DMSO_CD4pos_TNFAneg_IFNGpos", "ics_DMSO_CD4pos_TNFApos_IFNGneg", "ics_DMSO_CD4pos_TNFApos_IFNGpos")
	
	ICS_NS4B_CD8pos_subpops = c("ics_NS4B_CD8_CD25pos", "ics_NS4B_CD8_CD44pos", "ics_NS4B_CD8_CD44neg_CD62Lpos", "ics_NS4B_CD8_CD44neg_CD62Lneg", 
	"ics_NS4B_CD8_CD44pos_CD62Lneg", "ics_NS4B_CD8_CD44pos_CD62Lpos", "ics_NS4B_CD8_CD62Lneg", "ics_NS4B_CD8_IFNGpos", "ics_NS4B_CD8_IL_17pos", 
	"ics_NS4B_CD8_Tbetpos", "ics_NS4B_CD8_TNFApos", "ics_NS4B_CD8_TNFAneg_IFNGpos", "ics_NS4B_CD8_TNFApos_IFNGneg", "ics_NS4B_CD8_TNFApos_IFNGpos")
	
	ICS_NS4B_CD4pos_subpops = c("ics_NS4B_CD4_CD25pos", "ics_NS4B_CD4_CD44pos", "ics_NS4B_CD4_CD44neg_CD62Lpos", "ics_NS4B_CD4_CD44neg_CD62Lneg", 
	"ics_NS4B_CD4_CD44pos_CD62Lneg", "ics_NS4B_CD4_CD44pos_CD62Lpos", "ics_NS4B_CD4_CD62Lneg", "ics_NS4B_CD4_IFNGpos", "ics_NS4B_CD4_IL_17pos", 
	"ics_NS4B_CD4_Ly6Cpos", "ics_NS4B_CD4_Ly6Cneg_PSGL_1neg", "ics_NS4B_CD4_Ly6Cneg_PSGL_1pos", "ics_NS4B_CD4_Ly6Cpos_PSGL_1neg", 
	"ics_NS4B_CD4_Ly6Cpos_PSGL_1pos", "ics_NS4B_CD4_PSGL_1pos", "ics_NS4B_CD4_Tbetpos", "ics_NS4B_CD4_TNFApos", "ics_NS4B_CD4_TNFAneg_IFNGpos", 
	"ics_NS4B_CD4_TNFApos_IFNGneg", "ics_NS4B_CD4_TNFApos_IFNGpos")
	
	ICS_HIWNV_CD8pos_subpops = c("ics_HIWNV_CD8_CD25pos", "ics_HIWNV_CD8_CD44pos", "ics_HIWNV_CD8_CD44neg_CD62Lpos", "ics_HIWNV_CD8_CD44neg_CD62Lneg", 
	"ics_HIWNV_CD8_CD44pos_CD62Lneg", "ics_HIWNV_CD8_CD44pos_CD62Lpos", "ics_HIWNV_CD8_CD62Lneg", "ics_HIWNV_CD8_IFNGpos", "ics_HIWNV_CD8_IL_17pos", 
	"ics_HIWNV_CD8_Tbetpos", "ics_HIWNV_CD8_TNFApos", "ics_HIWNV_CD8_TNFAneg_IFNGpos", "ics_HIWNV_CD8_TNFApos_IFNGneg", "ics_HIWNV_CD8_TNFApos_IFNGpos")
	
	ICS_HIWNV_CD4pos_subpops = c("ics_HIWNV_CD4_CD25pos", "ics_HIWNV_CD4_CD44pos", "ics_HIWNV_CD4_CD44neg_CD62Lpos", "ics_HIWNV_CD4_CD44neg_CD62Lneg", 
	"ics_HIWNV_CD4_CD44pos_CD62Lneg", "ics_HIWNV_CD4_CD44pos_CD62Lpos", "ics_HIWNV_CD4_CD62Lneg", "ics_HIWNV_CD4_IFNGpos", "ics_HIWNV_CD4_IL_17pos", 
	"ics_HIWNV_CD4_Ly6Cpos", "ics_HIWNV_CD4_Ly6Cneg_PSGL_1neg", "ics_HIWNV_CD4_Ly6Cneg_PSGL_1pos", "ics_HIWNV_CD4_Ly6Cpos_PSGL_1neg", 
	"ics_HIWNV_CD4_Ly6Cpos_PSGL_1pos", "ics_HIWNV_CD4_PSGL_1pos", "ics_HIWNV_CD4_Tbetpos", "ics_HIWNV_CD4_TNFApos", "ics_HIWNV_CD4_TNFAneg_IFNGpos", 
	"ics_HIWNV_CD4_TNFApos_IFNGneg", "ics_HIWNV_CD4_TNFApos_IFNGpos")
	
	ICS_CD3CD28_CD8pos_subpops = c("ics_CD3CD28_CD8_CD25pos", "ics_CD3CD28_CD8_CD44pos", "ics_CD3CD28_CD8_CD44neg_CD62Lpos", 
	"ics_CD3CD28_CD8_CD44neg_CD62Lneg", "ics_CD3CD28_CD8_CD44pos_CD62Lneg", "ics_CD3CD28_CD8_CD44pos_CD62Lpos", "ics_CD3CD28_CD8_CD62Lneg", 
	"ics_CD3CD28_CD8_IFNGpos", "ics_CD3CD28_CD8_IL_17pos", "ics_CD3CD28_CD8_Tbetpos", "ics_CD3CD28_CD8_TNFApos", "ics_CD3CD28_CD8_TNFAneg_IFNGpos", 
	"ics_CD3CD28_CD8_TNFApos_IFNGneg", "ics_CD3CD28_CD8_TNFApos_IFNGpos")
	
	ICS_CD3CD28_CD4pos_subpops = c("ics_CD3CD28_CD4_CD25pos", "ics_CD3CD28_CD4_CD44pos", "ics_CD3CD28_CD4_CD44neg_CD62Lpos", 
	"ics_CD3CD28_CD4_CD44neg_CD62Lneg", "ics_CD3CD28_CD4_CD44pos_CD62Lneg", "ics_CD3CD28_CD4_CD44pos_CD62Lpos", "ics_CD3CD28_CD4_CD62Lneg", 
	"ics_CD3CD28_CD4_IFNGpos", "ics_CD3CD28_CD4_IL_17pos", "ics_CD3CD28_CD4_Ly6Cpos", "ics_CD3CD28_CD4_Ly6Cneg_PSGL_1neg", 
	"ics_CD3CD28_CD4_Ly6Cneg_PSGL_1pos", "ics_CD3CD28_CD4_Ly6Cpos_PSGL_1neg", "ics_CD3CD28_CD4_Ly6Cpos_PSGL_1pos", "ics_CD3CD28_CD4_PSGL_1pos", 
	"ics_CD3CD28_CD4_Tbetpos", "ics_CD3CD28_CD4_TNFApos", "ics_CD3CD28_CD4_TNFAneg_IFNGpos", "ics_CD3CD28_CD4_TNFApos_IFNGneg", 
	"ics_CD3CD28_CD4_TNFApos_IFNGpos")
	
	for (col in ICS_DMSO_CD8pos_subpops) {
		new_col = paste(col, '_count', sep='')
		flow_df[,new_col] = flow_df$ics_DMSO_CD8_count*(flow_df[,col]/100)
	}
	for (col in ICS_DMSO_CD4pos_subpops) {
		new_col = paste(col, '_count', sep='')
		flow_df[,new_col] = flow_df$ics_DMSO_CD4_count*(flow_df[,col]/100)
	}
	for (col in ICS_NS4B_CD8pos_subpops) {
		new_col = paste(col, '_count', sep='')
		flow_df[,new_col] = flow_df$ics_NS4B_CD8_count*(flow_df[,col]/100)
	}
	for (col in ICS_NS4B_CD4pos_subpops) {
		new_col = paste(col, '_count', sep='')
		flow_df[,new_col] = flow_df$ics_NS4B_CD4_count*(flow_df[,col]/100)
	}
	for (col in ICS_HIWNV_CD8pos_subpops) {
		new_col = paste(col, '_count', sep='')
		flow_df[,new_col] = flow_df$ics_HIWNV_CD8_count*(flow_df[,col]/100)
	}
	for (col in ICS_HIWNV_CD4pos_subpops) {
		new_col = paste(col, '_count', sep='')
		flow_df[,new_col] = flow_df$ics_HIWNV_CD4_count*(flow_df[,col]/100)
	}
	for (col in ICS_CD3CD28_CD8pos_subpops) {
		new_col = paste(col, '_count', sep='')
		flow_df[,new_col] = flow_df$ics_CD3CD28_CD8_count*(flow_df[,col]/100)
	}
	for (col in ICS_CD3CD28_CD4pos_subpops) {
		new_col = paste(col, '_count', sep='')
		flow_df[,new_col] = flow_df$ics_CD3CD28_CD4_count*(flow_df[,col]/100)
	}
	
	return(flow_df)
}
attr(calc_ics_counts, 'help') = "
This function calculates absolute cell counts for the ICS panels.

Parameters:
flow_df: The dataframe containing the flow cytometry data (total cell count, and percentages for all cell populations).

Returns:
A dataframe with columns for the cell counts added.
"

calc_ics_percent_ratios = function(flow_df) {
	ICS_NS4B_CD8pos_subpops = c("ics_NS4B_CD8_CD25pos", "ics_NS4B_CD8_CD44pos", "ics_NS4B_CD8_CD44neg_CD62Lpos", "ics_NS4B_CD8_CD44neg_CD62Lneg", 
	"ics_NS4B_CD8_CD44pos_CD62Lneg", "ics_NS4B_CD8_CD44pos_CD62Lpos", "ics_NS4B_CD8_CD62Lneg", "ics_NS4B_CD8_IFNGpos", "ics_NS4B_CD8_IL_17pos", 
	"ics_NS4B_CD8_Tbetpos", "ics_NS4B_CD8_TNFApos", "ics_NS4B_CD8_TNFAneg_IFNGpos", "ics_NS4B_CD8_TNFApos_IFNGneg", "ics_NS4B_CD8_TNFApos_IFNGpos")
	
	ICS_NS4B_CD4pos_subpops = c("ics_NS4B_CD4_CD25pos", "ics_NS4B_CD4_CD44pos", "ics_NS4B_CD4_CD44neg_CD62Lpos", "ics_NS4B_CD4_CD44neg_CD62Lneg", 
	"ics_NS4B_CD4_CD44pos_CD62Lneg", "ics_NS4B_CD4_CD44pos_CD62Lpos", "ics_NS4B_CD4_CD62Lneg", "ics_NS4B_CD4_IFNGpos", "ics_NS4B_CD4_IL_17pos", 
	"ics_NS4B_CD4_Ly6Cpos", "ics_NS4B_CD4_Ly6Cneg_PSGL_1neg", "ics_NS4B_CD4_Ly6Cneg_PSGL_1pos", "ics_NS4B_CD4_Ly6Cpos_PSGL_1neg", 
	"ics_NS4B_CD4_Ly6Cpos_PSGL_1pos", "ics_NS4B_CD4_PSGL_1pos", "ics_NS4B_CD4_Tbetpos", "ics_NS4B_CD4_TNFApos", "ics_NS4B_CD4_TNFAneg_IFNGpos", 
	"ics_NS4B_CD4_TNFApos_IFNGneg", "ics_NS4B_CD4_TNFApos_IFNGpos")
	
	ICS_HIWNV_CD8pos_subpops = c("ics_HIWNV_CD8_CD25pos", "ics_HIWNV_CD8_CD44pos", "ics_HIWNV_CD8_CD44neg_CD62Lpos", "ics_HIWNV_CD8_CD44neg_CD62Lneg", 
	"ics_HIWNV_CD8_CD44pos_CD62Lneg", "ics_HIWNV_CD8_CD44pos_CD62Lpos", "ics_HIWNV_CD8_CD62Lneg", "ics_HIWNV_CD8_IFNGpos", "ics_HIWNV_CD8_IL_17pos", 
	"ics_HIWNV_CD8_Tbetpos", "ics_HIWNV_CD8_TNFApos", "ics_HIWNV_CD8_TNFAneg_IFNGpos", "ics_HIWNV_CD8_TNFApos_IFNGneg", "ics_HIWNV_CD8_TNFApos_IFNGpos")
	
	ICS_HIWNV_CD4pos_subpops = c("ics_HIWNV_CD4_CD25pos", "ics_HIWNV_CD4_CD44pos", "ics_HIWNV_CD4_CD44neg_CD62Lpos", "ics_HIWNV_CD4_CD44neg_CD62Lneg", 
	"ics_HIWNV_CD4_CD44pos_CD62Lneg", "ics_HIWNV_CD4_CD44pos_CD62Lpos", "ics_HIWNV_CD4_CD62Lneg", "ics_HIWNV_CD4_IFNGpos", "ics_HIWNV_CD4_IL_17pos", 
	"ics_HIWNV_CD4_Ly6Cpos", "ics_HIWNV_CD4_Ly6Cneg_PSGL_1neg", "ics_HIWNV_CD4_Ly6Cneg_PSGL_1pos", "ics_HIWNV_CD4_Ly6Cpos_PSGL_1neg", 
	"ics_HIWNV_CD4_Ly6Cpos_PSGL_1pos", "ics_HIWNV_CD4_PSGL_1pos", "ics_HIWNV_CD4_Tbetpos", "ics_HIWNV_CD4_TNFApos", "ics_HIWNV_CD4_TNFAneg_IFNGpos", 
	"ics_HIWNV_CD4_TNFApos_IFNGneg", "ics_HIWNV_CD4_TNFApos_IFNGpos")
	
	ICS_CD3CD28_CD8pos_subpops = c("ics_CD3CD28_CD8_CD25pos", "ics_CD3CD28_CD8_CD44pos", "ics_CD3CD28_CD8_CD44neg_CD62Lpos", 
	"ics_CD3CD28_CD8_CD44neg_CD62Lneg", "ics_CD3CD28_CD8_CD44pos_CD62Lneg", "ics_CD3CD28_CD8_CD44pos_CD62Lpos", "ics_CD3CD28_CD8_CD62Lneg", 
	"ics_CD3CD28_CD8_IFNGpos", "ics_CD3CD28_CD8_IL_17pos", "ics_CD3CD28_CD8_Tbetpos", "ics_CD3CD28_CD8_TNFApos", "ics_CD3CD28_CD8_TNFAneg_IFNGpos", 
	"ics_CD3CD28_CD8_TNFApos_IFNGneg", "ics_CD3CD28_CD8_TNFApos_IFNGpos")
	
	ICS_CD3CD28_CD4pos_subpops = c("ics_CD3CD28_CD4_CD25pos", "ics_CD3CD28_CD4_CD44pos", "ics_CD3CD28_CD4_CD44neg_CD62Lpos", 
	"ics_CD3CD28_CD4_CD44neg_CD62Lneg", "ics_CD3CD28_CD4_CD44pos_CD62Lneg", "ics_CD3CD28_CD4_CD44pos_CD62Lpos", "ics_CD3CD28_CD4_CD62Lneg", 
	"ics_CD3CD28_CD4_IFNGpos", "ics_CD3CD28_CD4_IL_17pos", "ics_CD3CD28_CD4_Ly6Cpos", "ics_CD3CD28_CD4_Ly6Cneg_PSGL_1neg", 
	"ics_CD3CD28_CD4_Ly6Cneg_PSGL_1pos", "ics_CD3CD28_CD4_Ly6Cpos_PSGL_1neg", "ics_CD3CD28_CD4_Ly6Cpos_PSGL_1pos", "ics_CD3CD28_CD4_PSGL_1pos", 
	"ics_CD3CD28_CD4_Tbetpos", "ics_CD3CD28_CD4_TNFApos", "ics_CD3CD28_CD4_TNFAneg_IFNGpos", "ics_CD3CD28_CD4_TNFApos_IFNGneg", 
	"ics_CD3CD28_CD4_TNFApos_IFNGpos")
	
	## Calculate ICS Percentage Ratios
	flow_df$ics_NS4B_Live_ratio = flow_df$ics_NS4B_Live/flow_df$ics_DMSO_Live
	flow_df$ics_HIWNV_Live_ratio = flow_df$ics_HIWNV_Live/flow_df$ics_DMSO_Live
	flow_df$ics_CD3CD28_Live_ratio = flow_df$ics_CD3CD28_Live/flow_df$ics_DMSO_Live
	
	flow_df$ics_NS4B_Lymphocytes_ratio = flow_df$ics_NS4B_Lymphocytes/flow_df$ics_DMSO_Lymphocytes
	flow_df$ics_HIWNV_Lymphocytes_ratio = flow_df$ics_HIWNV_Lymphocytes/flow_df$ics_DMSO_Lymphocytes
	flow_df$ics_CD3CD28_Lymphocytes_ratio = flow_df$ics_CD3CD28_Lymphocytes/flow_df$ics_DMSO_Lymphocytes

	flow_df$ics_NS4B_CD3_ratio = flow_df$ics_NS4B_CD3/flow_df$ics_DMSO_CD3
	flow_df$ics_HIWNV_CD3_ratio = flow_df$ics_HIWNV_CD3/flow_df$ics_DMSO_CD3
	flow_df$ics_CD3CD28_CD3_ratio = flow_df$ics_CD3CD28_CD3/flow_df$ics_DMSO_CD3
	
	flow_df$ics_NS4B_CD8_ratio = flow_df$ics_NS4B_CD8/flow_df$ics_DMSO_CD8
	flow_df$ics_HIWNV_CD8_ratio = flow_df$ics_HIWNV_CD8/flow_df$ics_DMSO_CD8
	flow_df$ics_CD3CD28_CD8_ratio = flow_df$ics_CD3CD28_CD8/flow_df$ics_DMSO_CD8
	
	flow_df$ics_NS4B_CD4_ratio = flow_df$ics_NS4B_CD4/flow_df$ics_DMSO_CD4
	flow_df$ics_HIWNV_CD4_ratio = flow_df$ics_HIWNV_CD4/flow_df$ics_DMSO_CD4
	flow_df$ics_CD3CD28_CD4_ratio = flow_df$ics_CD3CD28_CD4/flow_df$ics_DMSO_CD4
	
	for (col in ICS_NS4B_CD8pos_subpops) {
		col_dmso = gsub("NS4B", "DMSO", col)
		col_dmso = gsub("CD8_", "CD8pos_", col_dmso)
		ratio_col = paste(col, '_ratio', sep='')
		flow_df[,ratio_col] = flow_df[,col]/flow_df[,col_dmso]
		flow_df[is.infinite(flow_df[,ratio_col]),ratio_col] = NA
		flow_df[is.nan(flow_df[,ratio_col]),ratio_col] = NA
	}
	for (col in ICS_NS4B_CD4pos_subpops) {
		col_dmso = gsub("NS4B", "DMSO", col)
		col_dmso = gsub("CD4_", "CD4pos_", col_dmso)
		ratio_col = paste(col, '_ratio', sep='')
		flow_df[,ratio_col] = flow_df[,col]/flow_df[,col_dmso]
		flow_df[is.infinite(flow_df[,ratio_col]),ratio_col] = NA
		flow_df[is.nan(flow_df[,ratio_col]),ratio_col] = NA
	}
	for (col in ICS_HIWNV_CD8pos_subpops) {
		col_dmso = gsub("HIWNV", "DMSO", col)
		col_dmso = gsub("CD8_", "CD8pos_", col_dmso)
		ratio_col = paste(col, '_ratio', sep='')
		flow_df[,ratio_col] = flow_df[,col]/flow_df[,col_dmso]
		flow_df[is.infinite(flow_df[,ratio_col]),ratio_col] = NA
		flow_df[is.nan(flow_df[,ratio_col]),ratio_col] = NA
	}
	for (col in ICS_HIWNV_CD4pos_subpops) {
		col_dmso = gsub("HIWNV", "DMSO", col)
		col_dmso = gsub("CD4_", "CD4pos_", col_dmso)
		ratio_col = paste(col, '_ratio', sep='')
		flow_df[,ratio_col] = flow_df[,col]/flow_df[,col_dmso]
		flow_df[is.infinite(flow_df[,ratio_col]),ratio_col] = NA
		flow_df[is.nan(flow_df[,ratio_col]),ratio_col] = NA
	}
	for (col in ICS_CD3CD28_CD8pos_subpops) {
		col_dmso = gsub("CD3CD28", "DMSO", col)
		col_dmso = gsub("CD8_", "CD8pos_", col_dmso)
		ratio_col = paste(col, '_ratio', sep='')
		flow_df[,ratio_col] = flow_df[,col]/flow_df[,col_dmso]
		flow_df[is.infinite(flow_df[,ratio_col]),ratio_col] = NA
		flow_df[is.nan(flow_df[,ratio_col]),ratio_col] = NA
	}
	for (col in ICS_CD3CD28_CD4pos_subpops) {
		col_dmso = gsub("CD3CD28", "DMSO", col)
		col_dmso = gsub("CD4_", "CD4pos_", col_dmso)
		ratio_col = paste(col, '_ratio', sep='')
		flow_df[,ratio_col] = flow_df[,col]/flow_df[,col_dmso]
		flow_df[is.infinite(flow_df[,ratio_col]),ratio_col] = NA
		flow_df[is.nan(flow_df[,ratio_col]),ratio_col] = NA
	}
	
	return(flow_df)
}
attr(calc_ics_percent_ratios, 'help') = "
This function calculates ratios of sub-population percentages for the ICS panels (DMSO as the reference).

Parameters:
flow_df: The dataframe containing the flow cytometry data (total cell count, and percentages for all cell populations).

Returns:
A dataframe with columns for the ratios added.
"

calc_ics_count_ratios = function(flow_df) {
	ICS_NS4B_CD8pos_subpops = c("ics_NS4B_CD8_CD25pos", "ics_NS4B_CD8_CD44pos", "ics_NS4B_CD8_CD44neg_CD62Lpos", "ics_NS4B_CD8_CD44neg_CD62Lneg", 
	"ics_NS4B_CD8_CD44pos_CD62Lneg", "ics_NS4B_CD8_CD44pos_CD62Lpos", "ics_NS4B_CD8_CD62Lneg", "ics_NS4B_CD8_IFNGpos", "ics_NS4B_CD8_IL_17pos", 
	"ics_NS4B_CD8_Tbetpos", "ics_NS4B_CD8_TNFApos", "ics_NS4B_CD8_TNFAneg_IFNGpos", "ics_NS4B_CD8_TNFApos_IFNGneg", "ics_NS4B_CD8_TNFApos_IFNGpos")
	
	ICS_NS4B_CD4pos_subpops = c("ics_NS4B_CD4_CD25pos", "ics_NS4B_CD4_CD44pos", "ics_NS4B_CD4_CD44neg_CD62Lpos", "ics_NS4B_CD4_CD44neg_CD62Lneg", 
	"ics_NS4B_CD4_CD44pos_CD62Lneg", "ics_NS4B_CD4_CD44pos_CD62Lpos", "ics_NS4B_CD4_CD62Lneg", "ics_NS4B_CD4_IFNGpos", "ics_NS4B_CD4_IL_17pos", 
	"ics_NS4B_CD4_Ly6Cpos", "ics_NS4B_CD4_Ly6Cneg_PSGL_1neg", "ics_NS4B_CD4_Ly6Cneg_PSGL_1pos", "ics_NS4B_CD4_Ly6Cpos_PSGL_1neg", 
	"ics_NS4B_CD4_Ly6Cpos_PSGL_1pos", "ics_NS4B_CD4_PSGL_1pos", "ics_NS4B_CD4_Tbetpos", "ics_NS4B_CD4_TNFApos", "ics_NS4B_CD4_TNFAneg_IFNGpos", 
	"ics_NS4B_CD4_TNFApos_IFNGneg", "ics_NS4B_CD4_TNFApos_IFNGpos")
	
	ICS_HIWNV_CD8pos_subpops = c("ics_HIWNV_CD8_CD25pos", "ics_HIWNV_CD8_CD44pos", "ics_HIWNV_CD8_CD44neg_CD62Lpos", "ics_HIWNV_CD8_CD44neg_CD62Lneg", 
	"ics_HIWNV_CD8_CD44pos_CD62Lneg", "ics_HIWNV_CD8_CD44pos_CD62Lpos", "ics_HIWNV_CD8_CD62Lneg", "ics_HIWNV_CD8_IFNGpos", "ics_HIWNV_CD8_IL_17pos", 
	"ics_HIWNV_CD8_Tbetpos", "ics_HIWNV_CD8_TNFApos", "ics_HIWNV_CD8_TNFAneg_IFNGpos", "ics_HIWNV_CD8_TNFApos_IFNGneg", "ics_HIWNV_CD8_TNFApos_IFNGpos")
	
	ICS_HIWNV_CD4pos_subpops = c("ics_HIWNV_CD4_CD25pos", "ics_HIWNV_CD4_CD44pos", "ics_HIWNV_CD4_CD44neg_CD62Lpos", "ics_HIWNV_CD4_CD44neg_CD62Lneg", 
	"ics_HIWNV_CD4_CD44pos_CD62Lneg", "ics_HIWNV_CD4_CD44pos_CD62Lpos", "ics_HIWNV_CD4_CD62Lneg", "ics_HIWNV_CD4_IFNGpos", "ics_HIWNV_CD4_IL_17pos", 
	"ics_HIWNV_CD4_Ly6Cpos", "ics_HIWNV_CD4_Ly6Cneg_PSGL_1neg", "ics_HIWNV_CD4_Ly6Cneg_PSGL_1pos", "ics_HIWNV_CD4_Ly6Cpos_PSGL_1neg", 
	"ics_HIWNV_CD4_Ly6Cpos_PSGL_1pos", "ics_HIWNV_CD4_PSGL_1pos", "ics_HIWNV_CD4_Tbetpos", "ics_HIWNV_CD4_TNFApos", "ics_HIWNV_CD4_TNFAneg_IFNGpos", 
	"ics_HIWNV_CD4_TNFApos_IFNGneg", "ics_HIWNV_CD4_TNFApos_IFNGpos")
	
	ICS_CD3CD28_CD8pos_subpops = c("ics_CD3CD28_CD8_CD25pos", "ics_CD3CD28_CD8_CD44pos", "ics_CD3CD28_CD8_CD44neg_CD62Lpos", 
	"ics_CD3CD28_CD8_CD44neg_CD62Lneg", "ics_CD3CD28_CD8_CD44pos_CD62Lneg", "ics_CD3CD28_CD8_CD44pos_CD62Lpos", "ics_CD3CD28_CD8_CD62Lneg", 
	"ics_CD3CD28_CD8_IFNGpos", "ics_CD3CD28_CD8_IL_17pos", "ics_CD3CD28_CD8_Tbetpos", "ics_CD3CD28_CD8_TNFApos", "ics_CD3CD28_CD8_TNFAneg_IFNGpos", 
	"ics_CD3CD28_CD8_TNFApos_IFNGneg", "ics_CD3CD28_CD8_TNFApos_IFNGpos")
	
	ICS_CD3CD28_CD4pos_subpops = c("ics_CD3CD28_CD4_CD25pos", "ics_CD3CD28_CD4_CD44pos", "ics_CD3CD28_CD4_CD44neg_CD62Lpos", 
	"ics_CD3CD28_CD4_CD44neg_CD62Lneg", "ics_CD3CD28_CD4_CD44pos_CD62Lneg", "ics_CD3CD28_CD4_CD44pos_CD62Lpos", "ics_CD3CD28_CD4_CD62Lneg", 
	"ics_CD3CD28_CD4_IFNGpos", "ics_CD3CD28_CD4_IL_17pos", "ics_CD3CD28_CD4_Ly6Cpos", "ics_CD3CD28_CD4_Ly6Cneg_PSGL_1neg", 
	"ics_CD3CD28_CD4_Ly6Cneg_PSGL_1pos", "ics_CD3CD28_CD4_Ly6Cpos_PSGL_1neg", "ics_CD3CD28_CD4_Ly6Cpos_PSGL_1pos", "ics_CD3CD28_CD4_PSGL_1pos", 
	"ics_CD3CD28_CD4_Tbetpos", "ics_CD3CD28_CD4_TNFApos", "ics_CD3CD28_CD4_TNFAneg_IFNGpos", "ics_CD3CD28_CD4_TNFApos_IFNGneg", 
	"ics_CD3CD28_CD4_TNFApos_IFNGpos")
	
	## Calculate ICS Count Ratios
	flow_df$ics_NS4B_Live_ratio_count = flow_df$ics_NS4B_Live_count/flow_df$ics_DMSO_Live_count
	flow_df$ics_HIWNV_Live_ratio_count = flow_df$ics_HIWNV_Live_count/flow_df$ics_DMSO_Live_count
	flow_df$ics_CD3CD28_Live_ratio_count = flow_df$ics_CD3CD28_Live_count/flow_df$ics_DMSO_Live_count
	
	flow_df$ics_NS4B_Lymphocytes_ratio_count = flow_df$ics_NS4B_Lymphocytes_count/flow_df$ics_DMSO_Lymphocytes_count
	flow_df$ics_HIWNV_Lymphocytes_ratio_count = flow_df$ics_HIWNV_Lymphocytes_count/flow_df$ics_DMSO_Lymphocytes_count
	flow_df$ics_CD3CD28_Lymphocytes_ratio_count = flow_df$ics_CD3CD28_Lymphocytes_count/flow_df$ics_DMSO_Lymphocytes_count
	
	flow_df$ics_NS4B_CD3_ratio_count = flow_df$ics_NS4B_CD3_count/flow_df$ics_DMSO_CD3_count
	flow_df$ics_HIWNV_CD3_ratio_count = flow_df$ics_HIWNV_CD3_count/flow_df$ics_DMSO_CD3_count
	flow_df$ics_CD3CD28_CD3_ratio_count = flow_df$ics_CD3CD28_CD3_count/flow_df$ics_DMSO_CD3_count
	
	flow_df$ics_NS4B_CD8_ratio_count = flow_df$ics_NS4B_CD8_count/flow_df$ics_DMSO_CD8_count
	flow_df$ics_HIWNV_CD8_ratio_count = flow_df$ics_HIWNV_CD8_count/flow_df$ics_DMSO_CD8_count
	flow_df$ics_CD3CD28_CD8_ratio_count = flow_df$ics_CD3CD28_CD8_count/flow_df$ics_DMSO_CD8_count
	
	flow_df$ics_NS4B_CD4_ratio_count = flow_df$ics_NS4B_CD4_count/flow_df$ics_DMSO_CD4_count
	flow_df$ics_HIWNV_CD4_ratio_count = flow_df$ics_HIWNV_CD4_count/flow_df$ics_DMSO_CD4_count
	flow_df$ics_CD3CD28_CD4_ratio_count = flow_df$ics_CD3CD28_CD4_count/flow_df$ics_DMSO_CD4_count
	
	for (col in ICS_NS4B_CD8pos_subpops) {
		count_col = paste(col, '_count', sep='')
		count_col_dmso = gsub("NS4B", "DMSO", count_col)
		count_col_dmso = gsub("CD8_", "CD8pos_", count_col_dmso)
		ratio_col = paste(col, '_ratio_count', sep='')
		flow_df[,ratio_col] = flow_df[,count_col]/flow_df[,count_col_dmso]
		flow_df[is.infinite(flow_df[,ratio_col]),ratio_col] = NA
		flow_df[is.nan(flow_df[,ratio_col]),ratio_col] = NA
	}
	for (col in ICS_NS4B_CD4pos_subpops) {
		count_col = paste(col, '_count', sep='')
		count_col_dmso = gsub("NS4B", "DMSO", count_col)
		count_col_dmso = gsub("CD4_", "CD4pos_", count_col_dmso)
		ratio_col = paste(col, '_ratio_count', sep='')
		flow_df[,ratio_col] = flow_df[,count_col]/flow_df[,count_col_dmso]
		flow_df[is.infinite(flow_df[,ratio_col]),ratio_col] = NA
		flow_df[is.nan(flow_df[,ratio_col]),ratio_col] = NA
	}
	for (col in ICS_HIWNV_CD8pos_subpops) {
		count_col = paste(col, '_count', sep='')
		count_col_dmso = gsub("HIWNV", "DMSO", count_col)
		count_col_dmso = gsub("CD8_", "CD8pos_", count_col_dmso)
		ratio_col = paste(col, '_ratio_count', sep='')
		flow_df[,ratio_col] = flow_df[,count_col]/flow_df[,count_col_dmso]
		flow_df[is.infinite(flow_df[,ratio_col]),ratio_col] = NA
		flow_df[is.nan(flow_df[,ratio_col]),ratio_col] = NA
	}
	for (col in ICS_HIWNV_CD4pos_subpops) {
		count_col = paste(col, '_count', sep='')
		count_col_dmso = gsub("HIWNV", "DMSO", count_col)
		count_col_dmso = gsub("CD4_", "CD4pos_", count_col_dmso)
		ratio_col = paste(col, '_ratio_count', sep='')
		flow_df[,ratio_col] = flow_df[,count_col]/flow_df[,count_col_dmso]
		flow_df[is.infinite(flow_df[,ratio_col]),ratio_col] = NA
		flow_df[is.nan(flow_df[,ratio_col]),ratio_col] = NA
	}
	for (col in ICS_CD3CD28_CD8pos_subpops) {
		count_col = paste(col, '_count', sep='')
		count_col_dmso = gsub("CD3CD28", "DMSO", count_col)
		count_col_dmso = gsub("CD8_", "CD8pos_", count_col_dmso)
		ratio_col = paste(col, '_ratio_count', sep='')
		flow_df[,ratio_col] = flow_df[,count_col]/flow_df[,count_col_dmso]
		flow_df[is.infinite(flow_df[,ratio_col]),ratio_col] = NA
		flow_df[is.nan(flow_df[,ratio_col]),ratio_col] = NA
	}
	for (col in ICS_CD3CD28_CD4pos_subpops) {
		count_col = paste(col, '_count', sep='')
		count_col_dmso = gsub("CD3CD28", "DMSO", count_col)
		count_col_dmso = gsub("CD4_", "CD4pos_", count_col_dmso)
		ratio_col = paste(col, '_ratio_count', sep='')
		flow_df[,ratio_col] = flow_df[,count_col]/flow_df[,count_col_dmso]
		flow_df[is.infinite(flow_df[,ratio_col]),ratio_col] = NA
		flow_df[is.nan(flow_df[,ratio_col]),ratio_col] = NA
	}
	
	return(flow_df)
}
attr(calc_ics_count_ratios, 'help') = "
This function calculates ratios of sub-population cell counts for the ICS panels (DMSO as the reference).

Parameters:
flow_df: The dataframe containing the flow cytometry data (total cell count, and percentages for all cell populations).

Returns:
A dataframe with columns for the ratios added.
"

clean_inf_nan = function(flow_df) {
	for (i in 1:dim(flow_df)[2]) {
		flow_df[is.infinite(flow_df[,i]),i] = NA
		flow_df[is.nan(flow_df[,i]),i] = NA
	}
	
	return(flow_df)
}
attr(calc_ics_percent_ratios, 'help') = "
This function cleans the results of the cell count and ratio calculations by removing an infinite or NaN values.

Parameters:
flow_df: The dataframe containing the full flow cytometry data.

Returns:
A cleaned dataframe.
"
