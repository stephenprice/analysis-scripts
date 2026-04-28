# script to calc. integrated and area-normalized fluxes from E3SM hist. files
#
# Assumed is that we are looking at a single, daily time step from a coupled E3SM
# +MALI simulation (30 min lnd, atm timesteps assumed, for total of 48 per day)

# load unified to access NCO
source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_chrysalis.sh

# remove previous versions of output .nc files and final output text file
rm ./cplHist.nc ./cplHist_l2x.nc ./cplHist_x2g.nc ./maliHist.nc ./FluxSumsOut.txt

# input relevant names for e3sm history files 
# Note that for some reason, it looks like the cpl hist files are being written out one day later than
# the mali hist file outputs ... or at least that is the only way the sums approximately equal each other below. 
# Could this be due to the layout (e.g. when lnd, cpl, and mali are being called ... in sequence for simple layout?
export cplHistFilePath=./20250219.BGWCYCL20TR.ne30pg2_r05_IcoswISC30E3r5_gis4to40.chrysalis.gnu.test.cpl.hi.2005-09-03-00000.nc #interested in day 2-24 
export maliHistFilePath=./20250219.BGWCYCL20TR.ne30pg2_r05_IcoswISC30E3r5_gis4to40.chrysalis.gnu.test.mali.hist.2005-09-01_00000.nc
#export cplHistFilePath=./20250219.BGWCYCL20TR.ne30pg2_r05_IcoswISC30E3r5_gis4to40.chrysalis.gnu.test.cpl.hi.2005-07-30-00000.nc #interested in day 7-29 
#export maliHistFilePath=./20250219.BGWCYCL20TR.ne30pg2_r05_IcoswISC30E3r5_gis4to40.chrysalis.gnu.test.mali.hist.2005-07-01_00000.nc

# earth radius = 6.37e6 m
# earth surf area = 5.09904e14 m^2

# extract vars from the 2nd time field of mali hist (first contains init cond information we don't want included)
ncks -d Time,1 -v sfcMassBal,sfcMassBalApplied,areaCell $maliHistFilePath ./maliHist.nc 
#ncks -d Time,27 -v sfcMassBal,sfcMassBalApplied,areaCell $maliHistFilePath ./maliHist.nc 

# calculate budget terms from mali hist file
ncap2 -A -s 'qice = sfcMassBal * areaCell' ./maliHist.nc
ncap2 -A -s 'qiceSum = qice.total()' ./maliHist.nc
ncap2 -A -s 'qiceSumNorm = qiceSum / (5.09904e14)' ./maliHist.nc
ncap2 -A -s 'qiceApplied = sfcMassBalApplied * areaCell' ./maliHist.nc
ncap2 -A -s 'qiceSumApplied = qiceApplied.total()' ./maliHist.nc
ncap2 -A -s 'qiceSumNormApplied = qiceSumApplied / (5.09904e14)' ./maliHist.nc

# extract relevant fields from cpl hist file and rename
cp $cplHistFilePath ./cplHist.nc

# x2g_ fluxes
ncks -v domg_area,domg_aream,x2g_Flgl_qice ./cplHist.nc ./cplHist_x2g.nc

# l2x_ fluxes for elevation class 10
ncks -v doml_area,doml_aream,l2x_Flgl_qice10,x2l_Sg_ice_covered10 ./cplHist.nc ./cplHist_l2x.nc 

# Loop over and append remaining elev classes 0-9 
# tgt="${blah}"
# ncks -v "${tgt}"
for i in {0..9}; do
	ncks -A -v l2x_Flgl_qice0${i},x2l_Sg_ice_covered0${i} ./cplHist.nc ./cplHist_l2x.nc
done

# calculate x2g_ budget terms from cpl hist file
# using model cell areas
ncap2 -A -s 'qice = x2g_Flgl_qice * domg_area * (6.37e6)^2' ./cplHist_x2g.nc
ncap2 -A -s 'qiceSum = qice.total()' ./cplHist_x2g.nc
ncap2 -A -s 'qiceSumNorm = qiceSum / (5.09904e14)' ./cplHist_x2g.nc
# using mapping file cell areas
ncap2 -A -s 'qiceMap = x2g_Flgl_qice * domg_aream * (6.37e6)^2' ./cplHist_x2g.nc
ncap2 -A -s 'qiceMapSum = qiceMap.total()' ./cplHist_x2g.nc
ncap2 -A -s 'qiceMapSumNorm = qiceMapSum / (5.09904e14)' ./cplHist_x2g.nc

# calculate l2x_ budget terms from cpl hist file
# using model cell areas
ncap2 -A -s 'qice10 = l2x_Flgl_qice10*x2l_Sg_ice_covered10*doml_area*(6.37e6)^2' ./cplHist_l2x.nc  
ncap2 -A -s 'qice9 = l2x_Flgl_qice09*x2l_Sg_ice_covered09*doml_area*(6.37e6)^2' ./cplHist_l2x.nc
ncap2 -A -s 'qice8 = l2x_Flgl_qice08*x2l_Sg_ice_covered08*doml_area*(6.37e6)^2' ./cplHist_l2x.nc
ncap2 -A -s 'qice7 = l2x_Flgl_qice07*x2l_Sg_ice_covered07*doml_area*(6.37e6)^2' ./cplHist_l2x.nc
ncap2 -A -s 'qice6 = l2x_Flgl_qice06*x2l_Sg_ice_covered06*doml_area*(6.37e6)^2' ./cplHist_l2x.nc
ncap2 -A -s 'qice5 = l2x_Flgl_qice05*x2l_Sg_ice_covered05*doml_area*(6.37e6)^2' ./cplHist_l2x.nc
ncap2 -A -s 'qice4 = l2x_Flgl_qice04*x2l_Sg_ice_covered04*doml_area*(6.37e6)^2' ./cplHist_l2x.nc
ncap2 -A -s 'qice3 = l2x_Flgl_qice03*x2l_Sg_ice_covered03*doml_area*(6.37e6)^2' ./cplHist_l2x.nc
ncap2 -A -s 'qice2 = l2x_Flgl_qice02*x2l_Sg_ice_covered02*doml_area*(6.37e6)^2' ./cplHist_l2x.nc
ncap2 -A -s 'qice1 = l2x_Flgl_qice01*x2l_Sg_ice_covered01*doml_area*(6.37e6)^2' ./cplHist_l2x.nc
ncap2 -A -s 'qice0 = l2x_Flgl_qice00*x2l_Sg_ice_covered00*doml_area*(6.37e6)^2' ./cplHist_l2x.nc
ncap2 -A -s 'qiceSum0 = (qice0+qice1+qice2+qice3+qice4+qice5+qice6+qice7+qice8+qice9+qice10)' ./cplHist_l2x.nc 
ncap2 -A -s 'qiceSum = qiceSum0.total()' ./cplHist_l2x.nc
ncap2 -A -s 'qiceSumNorm = qiceSum / (5.09904e14)' ./cplHist_l2x.nc

# using mapping file cell areas
ncap2 -A -s 'qice10m = l2x_Flgl_qice10*x2l_Sg_ice_covered10*doml_aream*(6.37e6)^2' ./cplHist_l2x.nc
ncap2 -A -s 'qice9m = l2x_Flgl_qice09*x2l_Sg_ice_covered09*doml_aream*(6.37e6)^2' ./cplHist_l2x.nc
ncap2 -A -s 'qice8m = l2x_Flgl_qice08*x2l_Sg_ice_covered08*doml_aream*(6.37e6)^2' ./cplHist_l2x.nc
ncap2 -A -s 'qice7m = l2x_Flgl_qice07*x2l_Sg_ice_covered07*doml_aream*(6.37e6)^2' ./cplHist_l2x.nc
ncap2 -A -s 'qice6m = l2x_Flgl_qice06*x2l_Sg_ice_covered06*doml_aream*(6.37e6)^2' ./cplHist_l2x.nc
ncap2 -A -s 'qice5m = l2x_Flgl_qice05*x2l_Sg_ice_covered05*doml_aream*(6.37e6)^2' ./cplHist_l2x.nc
ncap2 -A -s 'qice4m = l2x_Flgl_qice04*x2l_Sg_ice_covered04*doml_aream*(6.37e6)^2' ./cplHist_l2x.nc
ncap2 -A -s 'qice3m = l2x_Flgl_qice03*x2l_Sg_ice_covered03*doml_aream*(6.37e6)^2' ./cplHist_l2x.nc
ncap2 -A -s 'qice2m = l2x_Flgl_qice02*x2l_Sg_ice_covered02*doml_aream*(6.37e6)^2' ./cplHist_l2x.nc
ncap2 -A -s 'qice1m = l2x_Flgl_qice01*x2l_Sg_ice_covered01*doml_aream*(6.37e6)^2' ./cplHist_l2x.nc
ncap2 -A -s 'qice0m = l2x_Flgl_qice00*x2l_Sg_ice_covered00*doml_aream*(6.37e6)^2' ./cplHist_l2x.nc
ncap2 -A -s 'qiceMapSum0 = (qice0m+qice1m+qice2m+qice3m+qice4m+qice5m+qice6m+qice7m+qice8m+qice9m+qice10m)' ./cplHist_l2x.nc 
ncap2 -A -s 'qiceMapSum = qiceMapSum0.total()' ./cplHist_l2x.nc
ncap2 -A -s 'qiceMapSumNorm = qiceMapSum / (5.09904e14)' ./cplHist_l2x.nc

maliOut=$(ncdump -v qiceSumNorm maliHist.nc | tail -2)
maliOutApplied=$(ncdump -v qiceSumNormApplied maliHist.nc | tail -2)
x2gOut=$(ncdump -v qiceSumNorm cplHist_x2g.nc | tail -2)
x2gOutM=$(ncdump -v qiceMapSumNorm cplHist_x2g.nc | tail -2)
l2xOut=$(ncdump -v qiceSumNorm cplHist_l2x.nc | tail -2)
l2xOutM=$(ncdump -v qiceMapSumNorm cplHist_l2x.nc | tail -2)

dayOfYear=$(ncks -v time $cplHistFilePath | tail -5 | grep "time =")
dateOfYear=$(ncks -d Time,1 -v xtime $maliHistFilePath | tail -5)

echo ' ' > FluxSumsOut.txt
echo 'Noramlized flux sums' >> FluxSumsOut.txt
echo ' ' >> FluxSumsOut.txt
echo 'for day of year:'
echo $dayOfYear 
echo ' ' >> FluxSumsOut.txt
echo 'and date:'
echo $dateOfYear 
echo ' ' >> FluxSumsOut.txt
echo 'SMB flux in kg / m^2 / sec from MALI: ' >> FluxSumsOut.txt 
echo $maliOut >> FluxSumsOut.txt
echo ' ' >> FluxSumsOut.txt
echo 'SMB applied flux in kg / m^2 / sec from MALI: ' >> FluxSumsOut.txt 
echo $maliOutApplied >> FluxSumsOut.txt
echo ' ' >> FluxSumsOut.txt
echo 'SMB flux in kg / m^2 / sec (x2g) from CPL(model): ' >> FluxSumsOut.txt 
echo $x2gOut >> FluxSumsOut.txt
echo ' ' >> FluxSumsOut.txt 
echo 'SMB flux in kg / m^2 / sec (x2g) from CPL(map): ' >> FluxSumsOut.txt 
echo $x2gOutM >> FluxSumsOut.txt  
echo ' ' >> FluxSumsOut.txt
echo 'SMB flux in kg / m^2 / sec (l2x) from CPL(model): ' >> FluxSumsOut.txt 
echo $l2xOut >> FluxSumsOut.txt  
echo ' ' >> FluxSumsOut.txt
echo 'SMB flux in kg / m^2 / sec (l2x) from CPL(map): ' >> FluxSumsOut.txt 
echo $l2xOutM   >> FluxSumsOut.txt
echo ' ' >> FluxSumsOut.txt

cat FluxSumsOut.txt

