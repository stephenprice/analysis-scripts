#!/bin/bash

# This file uses nco to isolate and sum up the liquid water equivalent in the 
# snowpack over Greenland, which is extracted from ELM history files using an
# appropriate mask file. Once the filed is extracted, it is converted from mm
# of water equivalent to km of water equivalent thickness, multipled by the 
# area (in km), and then summed to get cubic km (or Gt) of total water. A few
# other time series variables are then calculated from that, e.g. the approximate
# sea level equivalent water stored (divde Gt by 360) and the change in equiv.
# sea level relative to the initial time step (as per above but subtracting off
# the initial value from the time series.

case='20260305.BGWCYCL2010.ne30pg2_r05_IcoswISC30E3r5_gis4to40.pm-cpu.testConfigNewSMBandIC'
mask='/global/cfs/cdirs/e3sm/sprice/temp/masks/glcmaskdata_0.5x0.5_GIS.20260310.nc'
year='2012'
fileout=h2osno${year}.nc

# source the necessary module to load nco, python
source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh

# use ncrcat to extract just snow liq water thickness (in mm) from 2005 elm h0 files
ncrcat -v H2OSNO,area ${case}.elm.h0.${year}*nc ${fileout}

# copy the mask file and update the dim names to be consistent with elm.h0 file
cp ${mask} ./gis-mask.nc

# replace dim names to be consistent w/ elm hist file dim names
ncrename -d lsmlat,lat gis-mask.nc 
ncrename -d lsmlon,lon gis-mask.nc

# append the mask field that isolates just GIS
ncks -A -v GLCMASK gis-mask.nc ${fileout}

# remove temp mask file
rm gis-mask.nc

# multiply the liqu water equiv field by the mask and convert from mm of liquid 
# to km of liquid and multiple by cell area (in km^2) to get km^3 = Gt of liq water
ncap2 -A -s 'maskedH2OSNOkm = H2OSNO/1000/1000*area*GLCMASK' ${fileout} 

# update attributes for new variable
ncatted -O -a long_name,maskedH2OSNOkm,o,c,"snow depth (liquid water) multipled by cell area and mask" ${fileout}
ncatted -O -a units,maskedH2OSNOkm,o,c,"cubic km or Gt" ${fileout} 

# sum previous over all cells to get time series of area integrated value 
ncap2 -A -s 'sumMaskedH2OSNOkm = maskedH2OSNOkm.total($lat,$lon);' ${fileout}

# update attributes for new variable
ncatted -O -a long_name,sumMaskedH2OSNOkm,o,c,"snow depth (liquid water) multipled by mask and cell area, summed over all cells" ${fileout}
ncatted -O -a units,sumMaskedH2OSNOkm,o,c,"cubic km or Gt" ${fileout} 

# convert timeseries in total km^3 = Gt of water to approx. global-mean sea leveli equiv.
ncap2 -A -s 'mmSeaLevelEquiv = sumMaskedH2OSNOkm/360' ${fileout}

# update attributes for new variable
ncatted -O -a long_name,mmSeaLevelEquiv,o,c,"snow depth (liquid water) multipled by mask and cell area, summed over all cells and converted to approx. global sea-level equiv." ${fileout}
ncatted -O -a units,mmSeaLevelEquiv,o,c,"mm of global sea-level equivalent" ${fileout} 
