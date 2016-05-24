library(gdata)
library(RColorBrewer)

# clean NA function
clean_na = function(data_set){
  
  for (i in 1:dim(data_set)[2]){
    #print(i)    
    if(sum(na.omit(data_set[,i] == "")) > 0){
      data_set[which(data_set[,i] == ""),i] <- NA
    }
    
  }
  return(data_set)
}

# Error bar function
error.bar = function(x, y, upper, lower=upper, arrow.length=0.05, ...) {
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper)) {
    print('x')
    print(length(x))
    print('y')
    print(length(y))
    print('upper')
    print(length(upper))
    stop("vectors must be same length")
  }
  arrows(x, y+upper, x, y-lower, angle=90, code=3, length=arrow.length, ...)
}

# Format data for plotting
qpcr_barplot_data = function(qpcr_byline_updated_cleaned, baseline=12, line, tissue, line_colors=NULL) {
  tp = as.character(sort(as.numeric(unique(qpcr_byline_updated_cleaned$Timepoint[qpcr_byline_updated_cleaned$UW_Line==line & tolower(qpcr_byline_updated_cleaned$Tissue)==tolower(tissue) & qpcr_byline_updated_cleaned$Virus=='WNV']))))
  if (length(tp) == 0) {
    return(NULL)
  }
  mock_tp = as.character(sort(as.numeric(unique(qpcr_byline_updated_cleaned$Timepoint[qpcr_byline_updated_cleaned$UW_Line==line & tolower(qpcr_byline_updated_cleaned$Tissue)==tolower(tissue) & qpcr_byline_updated_cleaned$Virus=='Mock']))))
  mock_tp = mock_tp[!mock_tp %in% baseline]
  
  exp_names = sort(tolower(unique(qpcr_byline_updated_cleaned$Experiment[qpcr_byline_updated_cleaned$UW_Line==line & tolower(qpcr_byline_updated_cleaned$Tissue)==tolower(tissue) & qpcr_byline_updated_cleaned$Virus=='WNV'])))
  if ('wnv' %in% exp_names) {
    if (length(exp_names) > 1) {
      exp_names = exp_names[c(length(exp_names), c(1:(length(exp_names)-1)))]
    }
  }
  
  mock_line_mat = matrix(ncol=0, nrow=length(mock_tp))
  mock_se_mat = matrix(ncol=0, nrow=length(mock_tp))
  line_mat = matrix(ncol=0, nrow=length(tp))
  se_mat = matrix(ncol=0, nrow=length(tp))
  
  ## Mock data
  if (length(mock_tp) > 0) {
    for (exp in exp_names) {
      exp_dat = c()
      exp_se = c()
      for (t in mock_tp) {
        tp_dat = qpcr_byline_updated_cleaned$fc.mean[with(qpcr_byline_updated_cleaned, UW_Line==line & tolower(Tissue)==tolower(tissue) & tolower(Experiment)==exp & Timepoint==t & Virus=='Mock')]
        #print(tp_dat)
        se_dat = qpcr_byline_updated_cleaned$ddCt.se[with(qpcr_byline_updated_cleaned, UW_Line==line & tolower(Tissue)==tolower(tissue) & tolower(Experiment)==exp & Timepoint==t & Virus=='Mock')]
        #print(se_dat)
        se_dat = 2^(-se_dat)
        #print(tp_dat)
        if (length(tp_dat) == 0) {
          exp_dat = c(exp_dat, NA)
          exp_se = c(exp_se, NA)
        } else {
          exp_dat = c(exp_dat, tp_dat)
          exp_se = c(exp_se, se_dat)
        }
      }
      mock_line_mat = cbind(mock_line_mat, exp_dat)
      mock_se_mat = cbind(mock_se_mat, exp_se)
    }
  }
  #print(mock_line_mat)
  
  ## WNV data
  for (exp in exp_names) {
    exp_dat = c()
    exp_se = c()
    for (t in tp) {
      tp_dat = qpcr_byline_updated_cleaned$fc.mean[with(qpcr_byline_updated_cleaned, UW_Line==line & tolower(Tissue)==tolower(tissue) & tolower(Experiment)==exp & Timepoint==t & Virus=='WNV')]
      se_dat = qpcr_byline_updated_cleaned$ddCt.se[with(qpcr_byline_updated_cleaned, UW_Line==line & tolower(Tissue)==tolower(tissue) & tolower(Experiment)==exp & Timepoint==t & Virus=='WNV')]
      #print(se_dat)
      se_dat = 2^(-se_dat)
      #print(tp_dat)
      if (length(tp_dat) == 0) {
        exp_dat = c(exp_dat, NA)
        exp_se = c(exp_se, NA)
      } else {
        exp_dat = c(exp_dat, tp_dat)
        exp_se = c(exp_se, se_dat)
      }
    }
    line_mat = cbind(line_mat, exp_dat)
    se_mat = cbind(se_mat, exp_se)
  }
  
  if (dim(mock_line_mat)[1] > 0) {
    line_mat = rbind(mock_line_mat, line_mat)
    se_mat = rbind(mock_se_mat, se_mat)
  }

  mock_names = c()
  if (length(mock_tp) > 0) {
    mock_names = paste('D', mock_tp, 'M', sep='')
  } 
  wnv_names = paste('D', tp, sep='')
  rownames(line_mat) = c(mock_names, wnv_names)
  colnames(line_mat) = toupper(exp_names)
  rownames(se_mat) = c(mock_names, wnv_names)
  colnames(se_mat) = toupper(exp_names)
  
  ## Return barplot matrix
  return(list(bardat=line_mat, sedat=se_mat, uw_line=line, tissue=tissue))
  
}

# Bar plot function
qpcr_barplot = function(barplot_data) {
  title = paste('UW Line', barplot_data$uw_line, '-', barplot_data$tissue, sep=' ')
  tp_colors = brewer.pal(dim(barplot_data$bardat)[1],'Blues')
  #print(min(log2(barplot_data$bardat), na.rm=T))
  barx = barplot(log2(barplot_data$bardat), beside=T, legend=T, args.legend=list(x='top', ncol=dim(barplot_data$bardat)[1]), ylim=c(min(log2(barplot_data$bardat), na.rm=T)-2, max(log2(barplot_data$bardat), na.rm=T)+5), xpd=FALSE, ylab='Log2(Fold Change)', main=title, col=tp_colors)
  error.bar(barx, log2(barplot_data$bardat), log2(barplot_data$sedat))
}