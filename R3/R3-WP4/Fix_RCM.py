###
# requires ncl and cdo
# and the first part python module load python at lus
# helenebe@met.no command list to average 3-hourly 4 km climate data
# dynamically downscaled using wrf over 9 years + to small files to use for esd
# will delete 3 hourly files, and keep variables Precip, t2, v10,u10,q2, tsk
# to daily means -- kg m2 to mm/day for precip

# makes monthly std of t2, v10,u10,q2, tsk
# makes monthly means of Precip, t2, v10,u10,q2, tsk
# makes monthly wet day mean and frequency

# finally hacks on dimensions to get lat lon from other wrf file and make file
# more cf-compatible, at least compatible for retrive.rcm in the
# R packege 'esd' 

import glob, os

folder=open('foldernameNORESM.txt', 'r').read().replace('\n',"") # filepath given in other file without cuotation marks

year=range(2000,2011,1)#  range(2000,2011,1) ##
header=open('http_server.txt', 'r').read().replace('\n',"")
for y in year:
  dirf=header+folder+str(y) +'/'  #+'/wrfxtrm*'
  os.system('wget -r --no-parent -A "wrfxtrm*" -nd -nc ' + dirf+ ' .')

# The remaining commands may be put in bash-script or included in
# python with os.system() or adding cdo and ncl
# but for transperancy I cp+paste and check files using ncdump and ncview
# as I go along now and then
# the commands, filenames and such may not be an order as old script is lost
# check carefully along the way!

'''
### merge all files with cat (works cause sorted alphanumerically) 
cdo -f nc4 cat wrfxtrm_d03* wrf_DS_NorESM_BC_2000_2010.nc

### Fix calendar for making daymean and monmean later  
cdo -f nc4 -z zip_2 selvar,RAINNCVMEAN,T2MEAN,U10MEAN,V10MEAN,SKINTEMPMEAN,Q2MEAN,T2MEAN wrf_DS_NorESM_BC_2000_2010.nc wrf_subset_DS_NorESM_BC_2000_2010.nc
# check if want >rm wrf_DS_NorESM_BC_2000_2010.nc
# calender is noleap 365_day
cdo settaxis,2000-12-01,00:00:00,3hour -setcalendar,365_day wrf_subset_DS_NorESM_BC_2000_2010.nc wrf_subset_NorESM_BC_2000_2010.nc

rm wrf_subset_DS_NorESM_BC_2000_2010.nc

cdo -f nc4 -z zip_2 daymean -selvar,T2MEAN,U10MEAN,V10MEAN,U10MEAN,SKINTEMPMEAN,Q2MEAN,T2MEAN wrf_subset_NorESM_BC_2000_2010.nc wrf_subsetday_TQUV_NorESM_BC_2000_2010.nc

cdo -f nc4 -z zip_2 monmean wrf_subsetday_TQUV_NorESM_BC_2000_2010.nc wrf_subsetmon_TQUV_NorESM_BC_2000_2010.nc

cdo -f nc4 -z zip_2 monstd wrf_subsetday_TQUV_NorESM_BC_2000_2010.nc wrf_subsetmon_std_TQUV_NorESM_BC_2000_2010.nc

### generating daily precip file in mm/day wrf_day_Pr_NorESM_BC_2000_2010.nc
cdo -f nc4 -z zip_2 mulc,86400 -daymean -selvar,RAINNCVMEAN wrf_subset_NorESM_BC_2000_2010.nc wrf_subsetday_P_NorESM_BC_2000_2010.nc
cdo chunit,'kg m-2 s-1','mm day^-1' wrf_subsetday_P_NorESM_BC_2000_2010.nc wrf_day_P_NorESM_BC_2000_2010.nc
rm wrf_subsetday_P_NorESM_BC_2000_2010.nc
ncatted -a description,RAINNCVMEAN,o,c,"Precipitation" wrf_day_P_NorESM_BC_2000_2010.nc wrf_day_Pr_NorESM_BC_2000_2010.nc
rm wrf_day_P_NorESM_BC_2000_2010.nc

### generate wet day frequency fw
cdo -f nc4 -z zip_2 divdpm -monsum -gec,1 wrf_day_Pr_NorESM_BC_2000_2010.nc wrf_mon_nwet_NorESM_BC_2000_2010.nc

cdo chunit,'mm day^-1','month^-1' -chname,"Precipitation",'fw' wrf_mon_nwet_NorESM_BC_2000_2010.nc wrf_mon_fwet_NorESM_BC_2000_2010.nc

ncatted -a description,fw,o,c,"Wet-day frequency (daily precip >= 1 mm day-1 )" wrf_mon_fwet_NorESM_BC_2000_2010.nc wrf_mon_fw_NorESM_BC_2000_2010.nc


#### generate wet day mean wrf_mon_mu_NorESM_BC_2000_2010.nc
cdo monmean -mul wrf_day_Pr_NorESM_BC_2000_2010.nc -gtc,1 wrf_day_Pr_NorESM_BC_2000_2010.nc wrf_mon_wetdaymean_NorESM_BC_2000_2010.nc

cdo chname,"Precipitation",'mu' wrf_mon_wetdaymean_NorESM_BC_2000_2010.nc tmp.nc

cdo monmean wrf_day_Pr_NorESM_BC_2000_2010.nc wrf_mon_Pr_NorESM_BC_2000_2010.nc
ncatted -a description,mu,o,c,"Wet-day mean precipitation" tmp.nc wrf_mon_mu_NorESM_BC_2000_2010.nc

#### check that these are here
wrf_mon_mu_NorESM_BC_2000_2010.nc #ok
wrf_mon_fw_NorESM_BC_2000_2010.nc
wrf_mon_Pr_NorESM_BC_2000_2010.nc
wrf_subsetmon_std_TQUV_NorESM_BC_2000_2010.nc


### std files change names and units
cdo chname,'T2MEAN','T2STD','Q2MEAN','Q2STD','SKINTEMPMEAN','SKINTEMPSTD','U10MEAN','U10STD','V10MEAN','V10STD' wrf_subsetmon_std_TQUV_NorESM_BC_2000_2010.nc wrf_mon_std_TQUV_NorESM_BC_2000_2010.nc

rm wrf_subsetmon_std_TQUV_NorESM_BC_2000_2010.nc

ncatted -O -a description,,a,c,'MONTHLY STANDARD DEVIATION' wrf_mon_std_TQUV_NorESM_BC_2000_2010.nc wrf_monstd_TQUV_NorESM_BC_2000_2010.nc

rm wrf_mon_std_TQUV_NorESM_BC_2000_2010.nc
#----------------------------------------


### finished steps merging
cdo merge wrf_monstd_TQUV_NorESM_BC_2000_2010.nc wrf_mon_mu_NorESM_BC_2000_2010.nc wrf_mon_fw_NorESM_BC_2000_2010.nc wrf_mon_Pr_NorESM_BC_2000_2010.nc wrf_subsetmon_TQUV_NorESM_BC_2000_2010.nc wrf_mon_sub_NorESM_BC_2000_2010.nc
#generate wrf_mon_sub_NorESM_BC_2000_2010.nc

cdo -f nc4 -z zip_2 merge wrf_day_Pr_NorESM_BC_2000_2010.nc wrf_subsetday_TQUV_NorESM_BC_2000_2010.nc wrf_day_sub_NorESM_BC_2000_2010.nc
#generate wrf_day_sub_NorESM_BC_2000_2010.nc

rm wrf_day_Pr_NorESM_BC_2000_2010.nc wrf_subsetday_TQUV_NorESM_BC_2000_2010.nc

rm wrf_monstd_TQUV_NorESM_BC_2000_2010.nc wrf_mon_mu_NorESM_BC_2000_2010.nc wrf_mon_fw_NorESM_BC_2000_2010.nc wrf_mon_Pr_NorESM_BC_2000_2010.nc wrf_subsetmon_TQUV_NorESM_BC_2000_2010.nc
#-----------------------------------------------

###  Get lat lon and fix coords in files 
ncrcat -d Time,0 -v XLAT,XLONG wrfpress_d03_2002-01-23_03\:00\:00 latlon.nc  
ncwa -a Time latlon.nc latlonn.nc
ncrename -h -d XLONG,longitude latlonn.nc 
ncrename -h -v XLONG,longitude latlonn.nc 
ncrename -h -v XLAT,latitude latlonn.nc 
ncrename -h -d west_east,x latlonn.nc 
ncrename -h -d south_north,y latlonn.nc 
ncatted -O -a cell_methods,,d,, latlonn.nc
ncatted -O -a coordinates,,d,, latlonn.nc
ncatted -O -a FieldType,,d,, latlonn.nc
ncatted -O -a stagger,,d,, latlonn.nc

#merging latlon with other files
cdo -f nc4 -z zip_2 merge wrf_day_sub_NorESM_BC_2000_2010.nc WRF_filetypes/latlonn.nc wrf_day_NorESM_BC_2000_2010.nc
cdo  -f nc4 merge wrf_mon_sub_NorESM_BC_2000_2010.nc WRF_filetypes/latlonn.nc wrf_mon_NorESM_BC_2000_2010.nc
cdo  -f nc4 -z zip_2 merge wrf_subset_NorESM_BC_2000_2010.nc WRF_filetypes/latlonn.nc wrf_sub3hr_NorESM_BC_2000_2010.nc

rm wrf_mon_sub_NorESM_BC_2000_2010.nc #is over
rm wrf_day_sub_NorESM_BC_2000_2010.nc #is over
rm wrf_subset_NorESM_BC_2000_2010.nc 
#------------------------------------------------------------
'''
