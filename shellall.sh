#!/bin/bash



./plot_surface_energy/shellplot_energy.sh &>>  log_energy &

./plot_mosi_cru/shellplot_cru.sh          &>>  log_mosi_cru &

./plot_tem_cru/shellplot_cru.sh           &>>  log_tem_cru  &

./plot_surface_moisture/shellplot_mois.sh &>>  log_surface_mos &

./plot_surface_tem/shellplot.sh           

./plot_swq_summer/shellplot_summer.sh     &>>  log_swq_summer &

