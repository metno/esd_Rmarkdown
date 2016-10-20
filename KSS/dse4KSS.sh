#!/bin/bash

echo "Submit downscaling jobs for KSS as parallel tasks"
for PARAM in t2m mu fw
#for PARAM in t2m #TEST
do  
  for SEASON in djf mam jja son 
  #for SEASON in djf # TEST
  do
    for RCP in rcp45 rcp85 rcp26
    #for RCP in rcp45
    do
      echo dse4KSS.$PARAM.$SEASON.$RCP.job
      sed 's|param|'$PARAM'|;
           s|season|'$SEASON'|;
           s|rcp|'$RCP'|' dse4KSS.job > dse4KSS.$PARAM.$SEASON.$RCP.job
      qsub dse4KSS.$PARAM.$SEASON.$RCP.job
    done
  done
done 


