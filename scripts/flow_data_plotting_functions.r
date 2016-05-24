library(gdata)
library(colorRamps)
library(RColorBrewer)
library(pheatmap)
options(warn=-1)

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
attr(error.bar, 'help') = "
This function adds error bars to plots.

Parameters:
x: A numeric vector of X values.
y: A numeric vector of Y values.
upper: A numeric vector containing the standard errors for the Y values.

Returns:
Plots arrows (error bars).
"

## Format the data for boxplots
flow_boxplot_data = function(flow_df, lines, tissue, flow_vars, tp, line_colors=NULL, mocks_only=FALSE, data_type=1) {
	## Initialize lists
	bdat = list()
	pdat = list()
	
	## Get Mock data
	n_mocks = 0
	mock_lines = c()
	for (l in lines) {
		label = paste(l, 'MOCK', sep='_')
		d = flow_df[flow_df$UW_Line==l & tolower(flow_df$Virus)=='mock' & flow_df$Tissue==tissue & flow_df$Timepoint %in% as.numeric(tp), flow_vars]
		if (!all(is.na(d))) {
			bdat[[label]] = d
			n_mocks = n_mocks + 1
			mock_lines = c(mock_lines, l)
		}
		tmp_list = list()
		for (t in tp) {
			d = flow_df[flow_df$UW_Line==l & tolower(flow_df$Virus)=='mock' & flow_df$Tissue==tissue & flow_df$Timepoint == as.numeric(t), flow_vars]
			if (!all(is.na(d))) {
				tmp_list[[t]] = d
			}
		}
		if (length(tmp_list) > 0) {
			pdat[[label]] = tmp_list
		}
	}
	
	## Get WNV data
	n_wnv = 0
	wnv_lines = c()
	if (!mocks_only) {
		for (l in lines) {
			label = paste(l, 'WNV', sep='_')
			d = flow_df[flow_df$UW_Line==l & tolower(flow_df$Virus)=='wnv' & flow_df$Tissue==tissue & flow_df$Timepoint %in% as.numeric(tp), flow_vars]
			if (!all(is.na(d))) {
				bdat[[label]] = d
				n_wnv = n_wnv + 1
				wnv_lines = c(wnv_lines, l)
			}
			tmp_list = list()
			for (t in tp) {
				d = flow_df[flow_df$UW_Line==l & tolower(flow_df$Virus)=='wnv' & flow_df$Tissue==tissue & flow_df$Timepoint == as.numeric(t), flow_vars]
				if (!all(is.na(d))) {
					tmp_list[[t]] = d
				}
			}
			if (length(tmp_list) > 0) {
				pdat[[label]] = tmp_list
			}
		}
	}
	
	## Create boxplot colors
	union_lines = union(mock_lines, wnv_lines)
	n_colors = length(union_lines)
	l_colors = colorRampPalette(brewer.pal(9,'Set1'))(n_colors)
	if (!mocks_only) {
		mock_colors = c()
		for (l in mock_lines) {
			mock_colors = c(mock_colors, l_colors[which(union_lines==l)])
		}
		wnv_colors = c()
		for (l in wnv_lines) {
			wnv_colors = c(wnv_colors, l_colors[which(union_lines==l)])
		}
		line_colors = c(mock_colors, wnv_colors)
	} else {
		line_colors = c(l_colors[1:n_mocks])
	}
	
	## DEBUGGING / TESTING
	#print(line_colors)
	
	## Determine data type (e.g. percentages vs. absolute counts)
	data_type = as.numeric(data_type)
	if (length(grep("^ics_.*", flow_vars)) > 0) {
		data_type = data_type + 2
	}
	
	## Create plot title with heritability annotation
	#icc_var1 = paste(tissue, '_mock_icc_gelman_hill', sep="")
	#icc_var2 = paste(tissue, '_icc_gelman_hill', sep="")
	#title = paste(flow_vars, " (ICC - Mock = ", round(flow_heritability[flow_vars, icc_var1], digits=3), "; ICC - All = ", round(flow_heritability[flow_vars, icc_var2], digits=3), ")", sep="")
	title = flow_vars
	
	## Return heatmap matrix
	return(list(bdat_list=bdat, pdat_list=pdat, colors=line_colors, num_mocks=n_mocks, time_points=tp, dtype=data_type, title=title))	
}
attr(flow_boxplot_data, 'help') = "
This function aggregates the data needed for creating boxplots of the flow cytometry data.

Parameters:
flow_df: The dataframe containing the flow data.
lines: A numeric vector containing the mouse lines that should be plotted.
tissue: A string indicating the tissue (e.g. 'brain' or 'spleen').
flow_vars: A string indicating the flow variable to be plotted.
tp: A character vector containing the time points to be included.
line_colors: A vector of colors (default=NULL; colors will be determined automatically).
mocks_only: Logical value indicating whether mocks only will be plotted.
data_type: A number indicating whether percentages (1) or absolute cell counts (2) will be plotted.

Returns:
A list containing the data to be plotted; input for the flow_boxplots() function.
"

## Boxplots
flow_boxplots = function(boxplot_data, ...) {
	tp_symbols = c("7"=21, "12"=22, "21"=23, "28"=24)
	if (boxplot_data$dtype == 1) {
		ylab="Percent"
	} else if (boxplot_data$dtype == 2) {
		ylab="Cell Count"
	} else if (boxplot_data$dtype == 3) {
		ylab="Percent Ratio"
	} else {
		ylab="Cell Count Ratio"
	}
	
	if (is.na(boxplot_data$y_min) | is.na(boxplot_data$y_max)) {
		if (boxplot_data$dtype == 1) {
			ylim = c(0,105)
		} else {
			ylim = NULL
		}
	} else {
		ylim = c(boxplot_data$y_min, boxplot_data$y_max)
	}
	
	if (boxplot_data$rm_outliers) {
		if (boxplot_data$dtype == 1) {
			boxplot(boxplot_data$bdat_list, col=boxplot_data$colors, xaxt = "n",  xlab = "", ylab=ylab, main=boxplot_data$title, outline=FALSE, ylim=ylim)
			if (boxplot_data$show_data) {
				for (t in boxplot_data$time_points) {
					d = lapply(boxplot_data$pdat_list, function(x){x[[t]]})
					stripchart(d, method='jitter', pch=tp_symbols[t], lwd=2, cex=1.25, vertical=T, add=T)
				}
			}
		} else {
			boxplot(boxplot_data$bdat_list, col=boxplot_data$colors, xaxt = "n",  xlab = "", ylab=ylab, main=boxplot_data$title, outline=FALSE, ylim=ylim)
			if (boxplot_data$show_data) {
				for (t in boxplot_data$time_points) {
					d = lapply(boxplot_data$pdat_list, function(x){x[[t]]})
					stripchart(d, method='jitter', pch=tp_symbols[t], lwd=2, cex=1.25, vertical=T, add=T)
				}
			}
		}
	} else {
		if (boxplot_data$dtype == 1) {
			boxplot(boxplot_data$bdat_list, col=boxplot_data$colors, xaxt = "n",  xlab = "", ylab=ylab, main=boxplot_data$title, pch=8, cex=1.5, ylim=ylim)
			if (boxplot_data$show_data) {
				for (t in boxplot_data$time_points) {
					d = lapply(boxplot_data$pdat_list, function(x){x[[t]]})
					stripchart(d, method='jitter', pch=tp_symbols[t], lwd=2, cex=1.25, vertical=T, add=T)
				}
			}
		} else {
			boxplot(boxplot_data$bdat_list, col=boxplot_data$colors, xaxt = "n",  xlab = "", ylab=ylab, main=boxplot_data$title, pch=8, cex=1.5, ylim=ylim)
			if (boxplot_data$show_data) {
				for (t in boxplot_data$time_points) {
					d = lapply(boxplot_data$pdat_list, function(x){x[[t]]})
					stripchart(d, method='jitter', pch=tp_symbols[t], lwd=2, cex=1.25, vertical=T, add=T)
				}
			}
		}
	}
	if (boxplot_data$show_data) {
		legend('topleft', legend=c('D7', 'D12', 'D21', 'D28'), pch=tp_symbols, pt.lwd=2, pt.cex=1.25, ncol=4)
	}
	labels = names(boxplot_data$bdat_list)
	axis(1, at=1:length(labels), labels=FALSE)
	y_interval = (par("usr")[4] - par("usr")[3])/675
	y_label_pos = par("usr")[3]-(80*y_interval)
	text(x=seq_along(labels), y=y_label_pos, srt=270, adj=1, labels=labels, xpd=TRUE, ...)
	abline(v=boxplot_data$num_mocks+0.5, lty=3)
	
	return(NULL)
}
attr(flow_boxplots, 'help') = "
This function plots the data returned by flow_boxplot_data().

Parameters:
boxplot_data: This should be a list combining the value returned by flow_boxplot_data() and additional options.

Returns:
NULL

Examples:
boxplot_data = flow_boxplot_data(flow_full, c(7,8,9), 'brain', 'treg_T_regs', c('7','12','21','28'))
opts = list(rm_outliers=F, show_data=F, y_min=0, y_max=100)
flow_boxplots(c(boxplot_data, opts))
"

flow_multiline_plot_data = function(flow_df, uw_lines, tissue, flow_vars, plot_type, FUN=NULL) {
	## Get time points
	## Only plots at time points where all variables exist
	if (plot_type == 1) {
		mock_lines = c()
		mock_tp = c()
		for (uw_line in uw_lines) {
			tp_existing = flow_df$Timepoint[flow_df$UW_Line==uw_line & flow_df$Virus == 'Mock' & flow_df$Tissue==tissue & !is.na(flow_df[,flow_vars])]
			if (length(tp_existing) > 0) {
				mock_tp = c(tp_existing, mock_tp)
				mock_lines = c(mock_lines, uw_line)
			}
		}
		mock_tp = sort(unique(mock_tp))
		wnv_lines = c()
		tp = c()
		for (uw_line in uw_lines) {
			tp_existing = flow_df$Timepoint[flow_df$UW_Line==uw_line & flow_df$Virus == 'WNV' & flow_df$Tissue==tissue & !is.na(flow_df[,flow_vars])]
			if (length(tp_existing) > 0) {
				tp = c(tp_existing, tp)
				wnv_lines = c(wnv_lines, uw_line)
			}
		}
		tp = sort(unique(tp))
		uw_lines = sort(as.numeric(unique(c(mock_lines, wnv_lines))))
	} else {
		new_vars = c()
		mock_tp = c()
		for (flow_var in flow_vars) {
			tp_existing = flow_df$Timepoint[flow_df$UW_Line==uw_lines & flow_df$Virus == 'Mock' & flow_df$Tissue==tissue & !is.na(flow_df[,flow_var])]
			if (length(tp_existing) > 0) {
				mock_tp = c(tp_existing, mock_tp)
				new_vars = c(new_vars, flow_var)
			}
		}
		mock_tp = sort(unique(mock_tp))
		tp = c()
		for (flow_var in flow_vars) {
			tp_existing = flow_df$Timepoint[flow_df$UW_Line==uw_lines & flow_df$Virus == 'WNV' & flow_df$Tissue==tissue & !is.na(flow_df[,flow_var])]
			if (length(tp_existing) > 0) {
				tp = c(tp_existing, tp)
				new_vars = c(new_vars, flow_var)
			}
		}
		tp = sort(unique(tp))
		flow_vars = unique(new_vars)
	}
	#print(mock_tp)
	#print(tp)
	#print(uw_lines)
	#print(flow_vars)
	
	if (length(mock_tp)+length(tp) == 0) {
		return(NULL)
	} 
	

	## Create data matrices
	if (plot_type == 1) {
		x = matrix(NA, nrow=(length(mock_tp)+length(tp)), ncol=length(uw_lines))
		y = matrix(NA, nrow=(length(mock_tp)+length(tp)), ncol=length(uw_lines))
		e = matrix(NA, nrow=(length(mock_tp)+length(tp)), ncol=length(uw_lines))
		l_mat = matrix(NA, nrow=(length(mock_tp)+length(tp)), ncol=length(uw_lines))
	} else {
		x = matrix(NA, nrow=(length(mock_tp)+length(tp)), ncol=length(flow_vars))
		y = matrix(NA, nrow=(length(mock_tp)+length(tp)), ncol=length(flow_vars))
		e = matrix(NA, nrow=(length(mock_tp)+length(tp)), ncol=length(flow_vars))
		l_mat = matrix(NA, nrow=(length(mock_tp)+length(tp)), ncol=length(flow_vars))
	}
	
	## Get mock data
	if (length(mock_tp) > 0) {
		if (plot_type == 1) {
		i = 1	
		end = length(mock_tp)
		for (uw_line in uw_lines) {
			if (uw_line %in% mock_lines) {
			## Get means for each group
			m = aggregate(flow_df[flow_df$UW_Line==uw_line & flow_df$Tissue==tissue & flow_df$Virus=='Mock', flow_vars], list(flow_df$Timepoint[flow_df$UW_Line==uw_line & flow_df$Tissue==tissue & flow_df$Virus=='Mock']), mean, na.rm=T)
			colnames(m) = c('group','y')
			for (t in setdiff(mock_tp, m$group)) {
				m = rbind(m, c(t, NA))
			}
			m = m[match(mock_tp,m[,1]),]
			m$x = c(1:end)
			#print(m)
			
			## Get SD for each group
			s = aggregate(flow_df[flow_df$UW_Line==uw_line & flow_df$Tissue==tissue & flow_df$Virus=='Mock', flow_vars], list(flow_df$Timepoint[flow_df$UW_Line==uw_line & flow_df$Tissue==tissue & flow_df$Virus=='Mock']), sd, na.rm=T)
			colnames(s) = c('group','y')
			for (t in setdiff(mock_tp, s$group)) {
				s = rbind(s, c(t, NA))
			}
			s = s[match(mock_tp,s[,1]),]
			
			## Get size of each group
			l = aggregate(flow_df[flow_df$UW_Line==uw_line & flow_df$Tissue==tissue & flow_df$Virus=='Mock', flow_vars], list(flow_df$Timepoint[flow_df$UW_Line==uw_line & flow_df$Tissue==tissue & flow_df$Virus=='Mock']), function(x){sum(!is.na(x))})
			colnames(l) = c('group','len')
			for (t in setdiff(mock_tp, l$group)) {
				l = rbind(l, c(t, NA))
			}
			l = l[match(mock_tp,l[,1]),]
			
			x[1:length(mock_tp),i] = m$x
			y[1:length(mock_tp),i] = m$y
			e[1:length(mock_tp),i] = s$y
			l_mat[1:length(mock_tp),i] = l$len
			} else {
				x[1:length(mock_tp),i] = NA
				y[1:length(mock_tp),i] = NA
				e[1:length(mock_tp),i] = NA
				l_mat[1:length(mock_tp),i] = NA
			}
			i = i + 1
		}
		} else {
		i = 1	
		end = length(mock_tp)
		for (flow_var in flow_vars) {
			#print(flow_var)
			## Get means for each group
			m = aggregate(flow_df[flow_df$UW_Line==uw_lines & flow_df$Tissue==tissue & flow_df$Virus=='Mock', flow_var], list(flow_df$Timepoint[flow_df$UW_Line==uw_lines & flow_df$Tissue==tissue & flow_df$Virus=='Mock']), mean, na.rm=T)
			#print(m)
			colnames(m) = c('group','y')
			for (t in setdiff(mock_tp, m$group)) {
				m = rbind(m, c(t, NA))
			}
			m = m[match(mock_tp,m[,1]),]
			m$x = c(1:end)
			
			## Get SD for each group
			s = aggregate(flow_df[flow_df$UW_Line==uw_lines & flow_df$Tissue==tissue & flow_df$Virus=='Mock', flow_var], list(flow_df$Timepoint[flow_df$UW_Line==uw_lines & flow_df$Tissue==tissue & flow_df$Virus=='Mock']), sd, na.rm=T)
			colnames(s) = c('group','y')
			for (t in setdiff(mock_tp, s$group)) {
				s = rbind(s, c(t, NA))
			}
			s = s[match(mock_tp,s[,1]),]
			
			## Get size of each group
			l = aggregate(flow_df[flow_df$UW_Line==uw_lines & flow_df$Tissue==tissue & flow_df$Virus=='Mock', flow_var], list(flow_df$Timepoint[flow_df$UW_Line==uw_lines & flow_df$Tissue==tissue & flow_df$Virus=='Mock']), function(x){sum(!is.na(x))})
			colnames(l) = c('group','len')
			for (t in setdiff(mock_tp, l$group)) {
				l = rbind(l, c(t, NA))
			}
			l = l[match(mock_tp,l[,1]),]
			
			x[1:length(mock_tp),i] = m$x
			y[1:length(mock_tp),i] = m$y
			e[1:length(mock_tp),i] = s$y
			l_mat[1:length(mock_tp),i] = l$len
			
			i = i + 1
		}
		}
	}
	
	#print(mock_tp)
	#print(tp)
	#print(x)
	#print(y)
	
	# Get WNV data
	if (length(tp) > 0) {
		if (plot_type == 1) {
		i = 1
		start = length(mock_tp)+1
		end = length(mock_tp)+length(tp)
		for (uw_line in uw_lines) {
			if (uw_line %in% wnv_lines) {
			## Get means for each group
			m = aggregate(flow_df[flow_df$UW_Line==uw_line & flow_df$Tissue==tissue & flow_df$Virus=='WNV', flow_vars], list(flow_df$Timepoint[flow_df$UW_Line==uw_line & flow_df$Tissue==tissue & flow_df$Virus=='WNV']), mean, na.rm=T)
			colnames(m) = c('group','y')
			for (t in setdiff(tp, m$group)) {
				m = rbind(m, c(t, NA))
			}
			m = m[match(tp,m[,1]),]
			m$x = c(start:end)
			
			## Get SD for each group
			s = aggregate(flow_df[flow_df$UW_Line==uw_line & flow_df$Tissue==tissue & flow_df$Virus=='WNV', flow_vars], list(flow_df$Timepoint[flow_df$UW_Line==uw_line & flow_df$Tissue==tissue & flow_df$Virus=='WNV']), sd, na.rm=T)
			colnames(s) = c('group','y')
			for (t in setdiff(tp, s$group)) {
				s = rbind(s, c(t, NA))
			}
			s = s[match(tp,s[,1]),]
			
			## Get size of each group
			l = aggregate(flow_df[flow_df$UW_Line==uw_line & flow_df$Tissue==tissue & flow_df$Virus=='WNV', flow_vars], list(flow_df$Timepoint[flow_df$UW_Line==uw_line & flow_df$Tissue==tissue & flow_df$Virus=='WNV']), function(x){sum(!is.na(x))})
			colnames(l) = c('group','len')
			for (t in setdiff(tp, l$group)) {
				l = rbind(l, c(t, NA))
			}
			l = l[match(tp,l[,1]),]
			
			x[start:end,i] = m$x
			y[start:end,i] = m$y
			e[start:end,i] = s$y
			l_mat[start:end,i] = l$len
			} else {
				x[start:end,i] = NA
				y[start:end,i] = NA
				e[start:end,i] = NA
				l_mat[start:end,i] = NA
			}
			
			i = i + 1
		}
		} else {
		i = 1
		start = length(mock_tp)+1
		end = length(mock_tp)+length(tp)
		for (flow_var in flow_vars) {
			## Get means for each group
			m = aggregate(flow_df[flow_df$UW_Line==uw_lines & flow_df$Tissue==tissue & flow_df$Virus=='WNV', flow_var], list(flow_df$Timepoint[flow_df$UW_Line==uw_lines & flow_df$Tissue==tissue & flow_df$Virus=='WNV']), mean, na.rm=T)
			colnames(m) = c('group','y')
			for (t in setdiff(tp, m$group)) {
				m = rbind(m, c(t, NA))
			}
			m = m[match(tp,m[,1]),]
			m$x = c(start:end)
			
			## Get SD for each group
			s = aggregate(flow_df[flow_df$UW_Line==uw_lines & flow_df$Tissue==tissue & flow_df$Virus=='WNV', flow_var], list(flow_df$Timepoint[flow_df$UW_Line==uw_lines & flow_df$Tissue==tissue & flow_df$Virus=='WNV']), sd, na.rm=T)
			colnames(s) = c('group','y')
			for (t in setdiff(tp, s$group)) {
				s = rbind(s, c(t, NA))
			}
			s = s[match(tp,s[,1]),]
			
			## Get size of each group
			l = aggregate(flow_df[flow_df$UW_Line==uw_lines & flow_df$Tissue==tissue & flow_df$Virus=='WNV', flow_var], list(flow_df$Timepoint[flow_df$UW_Line==uw_lines & flow_df$Tissue==tissue & flow_df$Virus=='WNV']), function(x){sum(!is.na(x))})
			colnames(l) = c('group','len')
			for (t in setdiff(tp, l$group)) {
				l = rbind(l, c(t, NA))
			}
			l = l[match(tp,l[,1]),]
			
			x[start:end,i] = m$x
			y[start:end,i] = m$y
			e[start:end,i] = s$y
			l_mat[start:end,i] = l$len
			
			i = i + 1
		}
		}
	}
	#print(x)
	#print(y)
	
	if (!is.null(FUN)) {
		y = FUN(y)
		e = FUN(e)
	}
	
	## return data list
	labels = c()
	if (length(mock_tp) > 0) {
		labels = paste(mock_tp, 'M', sep='')
	}
	labels = c(labels, as.character(tp))
	if (plot_type == 1) {
		return(list(x=x, y=y, e=e, l=l_mat, labels=labels, vars=uw_lines, num_vars=length(uw_lines), plot_type=plot_type, flow_var=flow_vars))
	} else {
		return(list(x=x, y=y, e=e, l=l_mat, labels=labels, vars=flow_vars, num_vars=length(flow_vars), plot_type=plot_type, uw_line=uw_lines))
	}
}
attr(flow_multiline_plot_data, 'help') = "
This function aggregates the data needed for creating time-series plots of the flow cytometry data. 

Parameters:
flow_df: The dataframe containing the flow data.
uw_lines: A numeric vector containing the mouse lines to be plotted.
tissue: A string indicating the tissue (e.g. 'brain' or 'spleen').
flow_vars: A character vector containing the variables to be plotted.
plot_type: A number indicating whether the plot will compare lines (1) or compare variables (2).
FUN: A function that can be applied to transform the data (default=NULL).

Returns:
A list containing the data to be plotted; input for the flow_multiline_plots() function.
"

flow_multiline_plots = function(lineplot_data) {
	## plot points and SE
	d_type = as.numeric(lineplot_data$data_type)
	if (length(grep("^ics_.*", lineplot_data$vars)) > 0) {
		d_type = d_type + 2
	}
	if (d_type == 1) {
		ylab="Percent"
	} else if (d_type == 2) {
		ylab="Cell Count"
	} else if (d_type == 3) {
		ylab="Percent Ratio"
	} else {
		ylab="Cell Count Ratio"
	}
	
	l_colors = colorRampPalette(brewer.pal(9,'Set1'))(lineplot_data$num_vars)
	
	if (is.na(lineplot_data$y_min) | is.na(lineplot_data$y_max)) {
		if (d_type == 1) {
			ylim = c(0,105)
		} else {
			ylim = NULL
		}
	} else {
		ylim = c(lineplot_data$y_min, lineplot_data$y_max)
	}
	
	if (lineplot_data$plot_type == 1) {
		title = lineplot_data$flow_var
		n_cols = lineplot_data$num_vars/2 + lineplot_data$num_vars%%2
	} else {
		title = paste("UW Line ", lineplot_data$uw_line, sep='')
		n_cols = 1
	}
	
	matplot(lineplot_data$x, lineplot_data$y, type="b", pch=19, lty=1, lwd=4, col=l_colors, ylim=ylim, xaxt="n", ylab=ylab, xlab="Time Point", main=title)
	#matplot(lineplot_data$x, lineplot_data$y, type="p", pch=19, col=l_colors, ylim=ylim, xaxt="n", ylab=ylab, xlab="Time Point", add=T)
	axis(1, at=lineplot_data$x[,1], labels=lineplot_data$labels)
	error.bar(lineplot_data$x, lineplot_data$y, lineplot_data$e/sqrt(lineplot_data$l), lwd=2, col=matrix(l_colors, nrow=length(lineplot_data$labels), ncol=lineplot_data$num_vars, byrow=T))
	legend("topleft", legend=lineplot_data$vars, col=l_colors, lty=1, lwd=4, ncol=n_cols)
	
	return(NULL)
}
attr(flow_multiline_plots, 'help') = "
This function plots the data returned by flow_multiline_plot_data().

Parameters:
lineplot_data: This should be a list combining the value returned by flow_multiline_plot_data() and additional options. 

Returns:
NULL

Examples:
lineplot_data = flow_multiline_plot_data(flow_full, c(7,8,9), 'brain', 'treg_T_regs', 1)
opts = list(data_type=1, y_min=0, y_max=60)
flow_multiline_plots(c(lineplot_data, opts))
"

flow_heatmap_data = function(flow_df, lines, line_labels=NULL, tissue, flow_vars, var_labels=NULL, tp, herit_df=NULL, line_colors=NULL, mocks_only=FALSE, collapse_mocks=FALSE, no_cluster=FALSE, cluster_all=FALSE, annotations=TRUE) {
	if (cluster_all) {
		if (mocks_only) {
			x = flow_df[flow_df$Tissue==tissue & tolower(flow_df$Virus)=='mock' & flow_df$Timepoint %in% tp,]
			if (collapse_mocks) {
				x$Timepoint[tolower(x$Virus)=='mock'] = 0
			}
		} else {
			x = flow_df[flow_df$Tissue==tissue,]
			if (collapse_mocks) {
				x$Timepoint[tolower(x$Virus)=='mock'] = 0
        x = x[x$Timepoint %in% c(0, tp), ]
			} else {
        x = x[x$Timepoint %in% tp,]
			}
		}
		y = aggregate(x[,flow_vars], list(x$UW_Line, x$Virus, x$Timepoint), mean, na.rm=T)
		colnames(y) = c('line', 'virus', 'tp', flow_vars)
		y = y[order(y$line, y$virus, y$tp),]
		rownames(y) = paste(y$line, y$virus, y$tp, sep='_')
		num_samples = dim(y)[1]
		var_indices = apply(y[,flow_vars], 2, function(x){sum(!is.na(x)) > num_samples*0.5})
		flow_vars_clust = flow_vars[var_indices]
		dist_matrix = dist(t(as.matrix(y[,flow_vars_clust])))
	} else {
		dist_matrix = "euclidean"
	}

	if (mocks_only) {
		x = flow_df[flow_df$UW_Line %in% lines & flow_df$Tissue==tissue & tolower(flow_df$Virus)=='mock' & flow_df$Timepoint %in% tp,]
		if (collapse_mocks) {
			x$Timepoint[tolower(x$Virus)=='mock'] = 0
		}
	} else {
		x = flow_df[flow_df$UW_Line %in% lines & flow_df$Tissue==tissue,]
		if (collapse_mocks) {
			x$Timepoint[tolower(x$Virus)=='mock'] = 0
      x = x[x$Timepoint %in% c(0, tp), ]
		} else {
      x = x[x$Timepoint %in% tp, ]
		}
	}
	y = aggregate(x[,flow_vars], list(x$UW_Line, x$Virus, x$Timepoint), mean, na.rm=T)
	colnames(y) = c('line', 'virus', 'tp', flow_vars)
	y$line = factor(y$line , levels = lines)
	y = y[order(y$line, y$virus, y$tp),]
	rownames(y) = paste(y$line, y$virus, y$tp, sep='_')
	
	blocks = c()
	for (l in lines) {blocks = c(blocks, max(which(y$line==l)))}

	## Column annotations	
	col_annot = y[,c('virus', 'line'), drop=F]
	col_annot$virus = as.factor(col_annot$virus)
	levels(col_annot$virus)[levels(col_annot$virus)=='Mock'] = 'M'
	levels(col_annot$virus)[levels(col_annot$virus)=='WNV'] = 'V'
	if (!is.null(line_labels)) {
		line_labels = factor(line_labels, levels=line_labels)
		col_annot$line = line_labels[col_annot$line]
	} else {
		col_annot$line = as.factor(col_annot$line)
	}
	#annot_colors = list(line=colorRampPalette(brewer.pal(9,'Set1'))(length(unique(y$line))), virus=c(M='gray', V='aquamarine'))
	annot_colors = list(virus=c(M='gray', V='white'))
	if (is.null(line_colors)) {
		#line_colors = sample(colorRampPalette(brewer.pal(9,'Set1'))(length(levels(col_annot$line))*3), length(levels(col_annot$line)))
		#line_colors = sample(primary.colors(length(levels(col_annot$line))*3), length(levels(col_annot$line)))
		#line_colors = c("#80FF00","#000000","#008080","#00FFFF","#800000","#FF00FF","#FF0000","#FF8000","#80FF80","#0080FF")
		line_colors = colorRampPalette(brewer.pal(9,'Set1'))(length(unique(y$line)))
		#line_colors = primary.colors(length(levels(col_annot$line)))
	}
		 
	names(line_colors) = levels(col_annot$line)
	annot_colors[['line']] = line_colors

	## Row and column labels
	num_samples = dim(y)[1]
	if (!no_cluster & !cluster_all) {
		var_indices = apply(y[,flow_vars], 2, function(x){sum(!is.na(x)) > num_samples*0.5})
		flow_vars = flow_vars[var_indices]
	} else if (cluster_all) {
		flow_vars = flow_vars_clust
	}
	#print(flow_vars)
	if (!is.null(labels)) {
		labels_row = var_labels
	} else {
		labels_row=gsub("_count","", flow_vars)
		labels_row=gsub("Lymphocytes","Lymph", labels_row)
	}
	labels_col = paste(y$virus, y$tp, sep='_')
	labels_col = gsub("WNV_", "  ", labels_col)
	if (collapse_mocks) {
		labels_col = gsub("Mock_.*", "  M", labels_col)
	} else {
		labels_col = gsub("Mock_", "  M", labels_col)
	}
	labels_col = gsub("$", "             ", labels_col)
	
	## Row annotations
	if (annotations & !is.null(herit_df)) {
		icc_var1 = paste(tissue, '_mock_icc_gelman_hill', sep="")
		icc_var2 = paste(tissue, '_icc_gelman_hill', sep="")
		row_annot = data.frame(ICC_mock=herit_df[flow_vars, icc_var1], ICC_all=herit_df[flow_vars, icc_var2])
		rownames(row_annot) = flow_vars
		annot_colors[['ICC_mock']] = colorRampPalette(brewer.pal(9,'YlGn'))(20)
		annot_colors[['ICC_all']] = colorRampPalette(brewer.pal(9,'YlGn'))(20)
	} else {
		row_annot = NULL
	}

	## Return heatmap matrix
	return(list(mat=t(as.matrix(y[,flow_vars])), cluster_rows=!no_cluster, dist_mat=dist_matrix, labels_rows=labels_row, labels_col=labels_col, col_annot=col_annot, row_annot=row_annot, blocks=blocks, annot_colors=annot_colors))
}
attr(flow_heatmap_data, 'help') = "
This function aggregates the data needed for creating heatmaps of the flow cytometry data. 

Parameters:
flow_df: The dataframe containing the flow data.
lines: A numeric vector containing the lines to be plotted.
line_labels: Alternate labels for the lines
tissue: A string indicating the tissue (e.g. 'brain' or 'spleen').
flow_vars: A character vector containing the variables to be plotted.
var_labels: Alternate labels for the flow variables
tp: A character vector containing the timepoints to be plotted.
herit_df: The dataframe containing the heritability estimates (default=NULL);
line_colors: A vector of color values (default=NULL; colors will be determined automatically).
mocks_only: A logical indicating whether mocks only should be plotted.
collapse_mocks: A logical indicating whether to combine the two mock timepoints.
no_cluster: A logical indicating whether variables should be clustered.
cluster_all: A logical indicating whether the variables should be clustered using all the data.
annotations: A logical indicating whether annotations should be added to the heatmap (default=TRUE).

Returns:
A list containing the data to be plotted; input for the flow_heatmap_plot() function.
"

weight_percents = c('D0_Percentage','D1_Percentage','D2_Percentage',
	'D3_Percentage','D4_Percentage','D5_Percentage','D6_Percentage',
	'D7_Percentage','D8_Percentage','D9_Percentage','D10_Percentage',
	'D11_Percentage','D12_Percentage','D13_Percentage','D14_Percentage',
	'D15_Percentage','D16_Percentage','D17_Percentage','D18_Percentage',
	'D19_Percentage','D20_Percentage','D21_Percentage','D22_Percentage',
	'D23_Percentage','D24_Percentage','D25_Percentage','D26_Percentage',
	'D27_Percentage','D28_Percentage')
	
cs_columns = c('D0','D1','D2','D3','D4','D5','D6','D7','D8','D9','D10',
	'D11','D12','D13','D14','D15','D16','D17','D18','D19','D20','D21',
	'D22','D23','D24','D25','D26','D27','D28')

flow_heatmap_plot = function(hm_data, weights_df=NULL, clinical_df=NULL, weight_cols=weight_percents, cs_cols=cs_columns, collapse_mocks=F, annotations=TRUE, ...) {
	## Add weights to heatmap
	if (annotations & !is.null(weights_df) & !is.null(clinical_df)) {
		if (collapse_mocks) {
			weights_df$Timepoint[tolower(weights_df$Virus)=='mock'] = 0
		}
		y_weights = aggregate(weights_df[, weight_cols], list(weights_df$UW_Line, weights_df$Virus, weights_df$Timepoint), mean, na.rm=T)
		y_weights$weight_change = NA
		for (i in 1:dim(y_weights)[1]) {
			w = y_weights[i,weight_cols]
			w = w[!is.nan(unlist(w)) & !is.na(unlist(w))]
			if (length(w) > 0) {
				y_weights$weight_change[i] = w[length(w)]
			}
		}
		colnames(y_weights) = c('line','virus','tp',weight_cols, 'weight_change')
		y_weights$tp[y_weights$tp=='d7'] = '7'
		y_weights$tp[y_weights$tp=='d12'] = '12'
		y_weights$tp[y_weights$tp=='d21'] = '21'
		y_weights$tp[y_weights$tp=='d28'] = '28'
		y_weights$tp[y_weights$tp=='d28m'] = '28'
		
		## DEBUGGING / TESTING
		#print(y_weights)
		
		annot_colnames = colnames(hm_data$col_annot)
		hm_data$col_annot$weight_loss = NA
		hm_data$col_annot = hm_data$col_annot[,c('weight_loss', annot_colnames)]
		for (i in 1:dim(hm_data$col_annot)[1]) {
			id_vector = unlist(strsplit(rownames(hm_data$col_annot)[i], "_"))
			w_change = unlist(y_weights$weight_change[y_weights$line==as.numeric(id_vector[1]) & tolower(y_weights$virus)==tolower(id_vector[2]) & y_weights$tp==id_vector[3]])
			if (!is.null(w_change) & length(w_change)>0) {
				hm_data$col_annot$weight_loss[i] = -w_change
			} else {
				hm_data$col_annot$weight_loss[i] = NA
			}
		}
		
		## DEBUGGING / TESTING
		#print(hm_data$col_annot)
		
		if (all(is.na(hm_data$col_annot$weight_loss))) {
			hm_data$col_annot = hm_data$col_annot[,!colnames(hm_data$col_annot) %in% c('weight_loss')]
		} else {
			hm_data$annot_colors[['weight_loss']] = colorRampPalette(brewer.pal(9,'RdYlBu'))(20)
		}
	
		## Add clinical score to heatmap
		if (collapse_mocks) {
			clinical_df$Timepoint[tolower(clinical_df$Virus)=='mock'] = 0
		}
		y_cs = aggregate(clinical_df[, cs_cols], list(clinical_df$UW_Line, clinical_df$Virus, clinical_df$Timepoint), max, na.rm=T)
		y_cs$clinical_score = NA
		for (i in 1:dim(y_cs)[1]) {
			cs = y_cs[i,cs_cols]
			cs = cs[!is.nan(unlist(cs)) & !is.na(unlist(cs)) & is.finite(unlist(cs))]
			if (length(cs) > 0) {
				y_cs$clinical_score[i] = max(cs, na.rm=T)
			}
		}
		colnames(y_cs) = c('line','virus','tp',cs_cols, 'clinical_score')
		y_cs$tp[y_cs$tp=='d7'] = '7'
		y_cs$tp[y_cs$tp=='d12'] = '12'
		y_cs$tp[y_cs$tp=='d21'] = '21'
		y_cs$tp[y_cs$tp=='d28'] = '28'
		y_cs$tp[y_cs$tp=='d28m'] = '28'
		
		## DEBUGGING / TESTING
		#print(y_cs)
		
		annot_colnames = colnames(hm_data$col_annot)
		hm_data$col_annot$clinical_score = NA
		hm_data$col_annot = hm_data$col_annot[,c('clinical_score', annot_colnames)]
		for (i in 1:dim(hm_data$col_annot)[1]) {
			id_vector = unlist(strsplit(rownames(hm_data$col_annot)[i], "_"))
			cs = unlist(y_cs$clinical_score[y_cs$line==as.numeric(id_vector[1]) & tolower(y_cs$virus)==tolower(id_vector[2]) & y_cs$tp==id_vector[3]])
			if (!is.null(cs) & length(cs)>0) {
				hm_data$col_annot$clinical_score[i] = cs
			} else {
				hm_data$col_annot$clinical_score[i] = NA
			}
		}
	
		if (all(is.na(hm_data$col_annot$clinical_score))) {
			hm_data$col_annot = hm_data$col_annot[,!colnames(hm_data$col_annot) %in% c('clinical_score')]
		} else {
			hm_data$col_annot$clinical_score = as.factor(hm_data$col_annot$clinical_score)
			cs_colors = colorRampPalette(brewer.pal(5,'Blues'))(length(levels(hm_data$col_annot$clinical_score)))
			names(cs_colors) = levels(hm_data$col_annot$clinical_score)
			hm_data$annot_colors[['clinical_score']] = cs_colors
		}
	} else {
		hm_data$col_annot = NULL
	}
	
	## Plot heatmap
	pheatmap(hm_data$mat, cluster_cols=F, cluster_rows=hm_data$cluster_rows, clustering_distance_rows=hm_data$dist_mat, scale='row', gaps_col=hm_data$blocks, labels_col=hm_data$labels_col, labels_row=hm_data$labels_row, color=colorRampPalette(brewer.pal(9,'RdYlBu'))(20), annotation_row=hm_data$row_annot, annotation_col=hm_data$col_annot, annotation_colors=hm_data$annot_colors, fontsize=11, fontsize_col=9, ...)

	return(NULL)
}
attr(flow_heatmap_plot, 'help') = "
This function plots the data returned by flow_heatmap_data().

Parameters:
hm_data: This should be a list returned by flow_heatmap_data().
weights_df: The dataframe containing the weight loss data (default=NULL).
clinical_df: The dataframe containing the clinical score data (default=NULL).
weight_cols: A character vector containing the columns names of the weight measurements (default=weight_percents). 
cs_cols: A character vector containing the columns names of the clinical scores (default=cs_columns). 
annotations: A logical indicating whether annotations should be added to the heatmap (default=TRUE).

Returns:
NULL

Examples:
flow_heatmap_plot(heatmap_data, weights, scores)

"
