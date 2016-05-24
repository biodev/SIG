
## Load functions for plotting the flow cytometry data
source('./scripts/flow_data_plotting_functions.r')

## Note: you may have to change the file paths
data_dir = '/Users/mooneymi/Documents/SIG/WNV/Cleaned_Data_Releases/23-Mar-2016/'

## Load data from an Excel spreadsheet (Warning: this can take a few minutes)
flow_full = read.xls(file.path(data_dir, 'Lund_Flow_Full_21-Mar-2016_final.xlsx'), 
                     header=T, as.is=T, na.strings=c(""," ", "NA", "#DIV/0!"))

describe(flow_boxplot_data)

## Aggregate the data for boxplots
boxplot_data = flow_boxplot_data(flow_full, c(7,8,9), 'brain', 'treg_T_regs', c('7','12','21','28'))

## Create a list of additional options for the boxplot
opts = list(rm_outliers=F, show_data=F, y_min=0, y_max=60)
## Create the boxplot (the 'cex' parameter controls the size of the x-axis text)
bp = flow_boxplots(c(boxplot_data, opts), cex=0.7)

describe(flow_multiline_plot_data)

lineplot_data = flow_multiline_plot_data(flow_full, c(30,8,36,38), 'brain', 'treg_T_regs_count', 1)

## Create a list of additional options for the lineplot
## data_type values: 1 = percentages, 2 = cell counts, 3 = percent ratio, 4 = count ratio
opts2 = list(data_type=2, y_min=NA, y_max=NA)
## Create a lineplot that compares a single variable across multiple lines
lp = flow_multiline_plots(c(lineplot_data, opts2))

lineplot_data2 = flow_multiline_plot_data(flow_full, c(9), 'brain', c('treg_T_regs', 'tcell_d7_CD8'), 2)

## Create a list of additional options for the lineplot
## data_type values: 1 = percentages, 2 = cell counts, 3 = percent ratio, 4 = count ratio
opts3 = list(data_type=1, y_min=0, y_max=50)
## Create a lineplot that compares multiple variables for a single line
lp2 = flow_multiline_plots(c(lineplot_data2, opts3))

## Load weight, clinical score, and heritability data from the latest data release
## Note: you may have to change the file paths
weights = read.xls(file.path(data_dir, 'Lund_Weight_22-Mar-2016_final.xlsx'),
                        header=T, as.is=T, na.strings=c(""," ", "NA", "#DIV/0!"))
scores = read.xls(file.path(data_dir, 'Lund_Scores_22-Mar-2016_final.xlsx'),
                       header=T, as.is=T, na.strings=c(""," ", "NA", "#DIV/0!"))
heritability = read.xls(file.path(data_dir, 'Lund_Flow_Heritability_21-Mar-2016_final.xlsx'),
                        header=T, as.is=T, na.strings=c(""," ", "NA", "#DIV/0!"))

## Set the rownames of the heritability dataframe
rownames(heritability) = heritability$variable

describe(flow_heatmap_data)

describe(flow_heatmap_plot)

## Heatmap with custom labels, mocks collapsed, and no heritability annotations
heatmap_data = flow_heatmap_data(flow_full, lines=c(11,12,14,30,8,36,38),
                line_labels=c('CC(017x004)F1','CC(011x042)F1','CC(032x017)F1','CC(032x013)F1','CC(005x001)F1','CC(061x026)F1','CC(016x038)F1'),
                tissue='brain', 
                flow_vars=c('treg_T_regs', 'tcell_d7_CD3', 'tcell_d7_CD4', 'tcell_d7_CD8'), 
                var_labels=c('Tregs', 'CD3+ Tcell', 'CD4+ Tcell', 'CD8+ Tcell'), 
                tp=c('7','12','21','28'), collapse_mocks=T)

## Create the heatmap
hm = flow_heatmap_plot(heatmap_data, weights, scores, collapse_mocks=T, annotations=T)

## The heatmap without any annotations
heatmap_data2 = flow_heatmap_data(flow_full, lines=c(7,8,9), tissue='brain', 
                                 flow_vars=c('treg_T_regs', 'tcell_d7_CD3', 'tcell_d7_CD4', 'tcell_d7_CD8'),
                                 tp=c('7','12','21','28'), annotations=F)
hm2 = flow_heatmap_plot(heatmap_data2, annotations=F)
