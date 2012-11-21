This is the program designed for make four panel plots


the example data are VIC60, ViC70, ERA, CRU


shellplot_cru.sh is the shell script to run and create the plot
Since the file name of different time periods are different,
the loop is divided into three groups: seasons, month, years


plot_cru.py is the script to tell python which variables you want


makeplot_cru.py is the script to set some setting of the plot,they are:
lorange, larange, colormp, labels, subtitle, title, pathout,
outnammed(outfile name), outformat(output format), projection_parameters

if the units are different among input files the number need 
to be adjusted in makeplot_cru.py script


layout.py is the basic settings for the plot. If you want to change
the position of subplot or some text, you can change it here.


cut_colmap.py is the code from internet to make the color map discrete

