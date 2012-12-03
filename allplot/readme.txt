
shellall.sh is the script to make all the plots

base on different input file names, there are 6 directories store the
6 groups of python and shell script:
./plot_surface_energy: for Latht, senht, swnet, Lwnet, swin, Lwin
./plot_swq_summer: for swq in summer time  
./plot_mosi_cru: for cru data's precipitation
./plot_tem_cru: for cru data's air temperature
./plot_surface_tem: for Tair Radt Swq
./plot_surface_moisture: for Precipitation Evap Runoff



Firstly it will run 4 background jobs and start the fifth on the screen

after the 5th is finished, the last jobs will run in the backgroud,
this is because the last jobs make plots to overwrite some of the fifth,
so please make sure the 5th(./plot_surface_tem/shellplot.sh)
and the 6th(./plot_swq_summer/shellplot_summer.sh) shellscript output 
directories(variable "figp") are the same.



so it will generate 5 logfile, they are: log_energy, log_mosi_cru,
log_tem_cru, log_surface_mos, log_swq_summer



in each of the 6 directories(plot_mosi_cru, plot_surface_moisture, plot_swq_summer,
plot_surface_energy, plot_surface_tem, plot_tem_cru), there is shell script
named shellplot*.sh. The input data and output directory of the figure are provided here.
"erap, figp, rasm2, rasm1, rasm3" are the variables for that information
 
please notice, the variables:'rasm1', 'rasm2', 'rasm3', and 'erap' stand for the data used
corresponding to the 1st, 2nd, 3rd and the last subplots in the figure. The
ERA and CRU data are only two dimensions but the rasm data are 3 dimensions
with one more time-dimension. please see the data used in these scripts as
example. variable 'figp' is the output directory of the figure

And please do NOT make the output directory of cru-plots the same with non-cru plots
because the output name of the figure file are the same so some figures will be
overwrited.





