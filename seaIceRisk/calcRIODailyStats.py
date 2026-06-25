# Note that this ia a .py script for the notebook script with the same name.
# It is useful for processing  the full dataset via a few years at a time, farmed
# out to multiple procs on a single node (e.g., by copying this scropt and changing 
# the dates it is being run for only). Use along with 'RIObatch.sh'

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import xarray
import os
import math

# list of ensembles to read in and calc stats on
ensembles = [ '0051', '0091', '0101', '0111', '0121',
              '0131', '0141', '0151', '0161', '0171',
              '0181', '0191', '0201', '0211', '0221',
              '0231', '0241', '0251', '0261', '0271',
              '0281', '0291', '0301', '0311', '0321' ]
#ensembles = [ '0051', '0151', '0251' ]

years = [ str(y) for y in range(2000, 2002) ]
#years = [ '2000' ]

months = [ '01', '02', '03', '04', '05', '06', '07',
           '08', '09', '10', '11', '12' ]
#months = [ '04' ]

classes = [ 'PC1', 'PC2', 'PC3', 'PC4',
            'PC5', 'PC6', 'PC7', 'IAsuper',
            'IA', 'IB', 'IC', 'NIS' ]
#classes = [ 'PC6' ]

path='/lustre/scratch5/sprice/data/'
#path='/Users/stephenprice/projects/CINS/Arctic-2024/ColdHarbor/testData2'
os.chdir(path)

# read in lat and lon information from mesh and find cell indices to restrict data to >=60 North
meshFileName = 'mpaso-IcoswISC30E3r5-restart.nc'
print('reading IN file: ', meshFileName)

dataIn = xarray.open_dataset(meshFileName)

latIn = dataIn.latCell.values
lonIn = dataIn.lonCell.values

rad2deg = 180.0 / np.pi

# find indices and generate new lat, lon data vectors based on these indices
ind = np.where( (latIn * rad2deg) >= 60.0)[0]

latOut = latIn[ind]
lonOut = lonIn[ind]

# save these for later use (or just re-do this step in other code?)

# set up some strings for manipulating paths and filenames
pathBase='/lustre/scratch5/sprice/'
pathWriteBase='/lustre/scratch5/sprice/data/'
#pathBase='/Users/stephenprice/projects/CINS/Arctic-2024/ColdHarbor/testData2/'
os.chdir(pathBase)

prefix='v3.LR.historical_'
dirs = '/hist/'
mpassiPrefix = '.mpassi.hist.am.timeSeriesStatsDaily.'
writePath = pathWriteBase+'v3.LR.historical_seaIceRIO_ensembleStatsNEW'

from polaris_rio import polar_rio
import time

t0 = time.perf_counter()


# calc. ensemble averages and stdev by looping over ensembles, for specific months, for specific years
# (resulting outputs will be written to montly netcdf files, as per the single ensemble datasets being
# read in)

for year in years:

    for month in months:

        first = True
        ensemble_index = 0

        for ensemble in ensembles:

            fileBase = pathBase+prefix+ensemble+dirs
            file = prefix+ensemble+mpassiPrefix+year+'-'+month+'-01.nc'

            # read in file and relevant vars
            print('reading IN file: ', fileBase+file)
            dataIn = xarray.open_dataset(fileBase+file)
            # dataIn = xarray.open_dataset(fileBase+file, drop_variables=[
            #     'timeDaily_avg_snowVolumeCell', 'timeDaily_avg_uVelocityGeo',
            #     'timeDaily_avg_vVelocityGeo' ] )

            concentration = dataIn.timeDaily_avg_iceAreaCell.values[:,ind]
            volume = dataIn.timeDaily_avg_iceVolumeCell.values[:,ind]
            #thickness = volume / concentration 
            thickness = np.where(concentration > 0.01, volume / concentration, 0.0)
            # concentration = dataIn["timeDaily_avg_iceAreaCell"].isel(nCells=ind).to_numpy()
            # thickness = dataIn["timeDaily_avg_iceVolumeCell"].isel(nCells=ind).to_numpy()
            # dataIn.close()

            ndays = np.size( thickness, axis=0 )
            ncells = np.size( thickness, axis=1 )
            nclasses = len( classes )
            nensembles = len( ensembles )

            if first:
                # dimension empty storage array 
                RIOArray = np.ones([ nensembles, nclasses, ndays, ncells] ) * 30
                first = False

            # change into relevant ensemble subdir and read in monthly .nc files
            fileBase = pathBase+prefix+ensemble+dirs
            file = prefix+ensemble+mpassiPrefix+year+'-'+month+'-01.nc'

            # loop over data and fill array w/ RIO values

            # loop over vessels, days, and cells to assign RIO value for combination of 
            # ice thickness and concentration for that vessel, at that day and location

            vessel_index = 0

            for vessel in classes:

                for day in np.arange(0,ndays):
                # for day in np.arange(0,1):

                    for c in np.arange(0,ncells):

                        rio, ice_type, level = polar_rio( vessel_ice_class=vessel,
                                                 ice_thickness_m=thickness[day,c],
                                                 concentration_tenths=concentration[day,c]*10
                                                )
                        # append data to pre-dimensioned array
                        RIOArray[ ensemble_index, vessel_index, day, c ] = rio

                vessel_index = vessel_index + 1

            ensemble_index = ensemble_index + 1

        # calc stats (collapse along ensemble dimension) and write to new .nc file
        RIOmedian = np.median( RIOArray, axis=0 )
        RIO5th = np.percentile( RIOArray, 5, axis=0 )
        RIO25th = np.percentile( RIOArray, 25, axis=0 )
        RIO75th = np.percentile( RIOArray, 75, axis=0 )
        RIO95th = np.percentile( RIOArray, 95, axis=0 )

        dataOutIndex = 0

        for vessel in classes:

            # write ensemble stats back out to (monthly) .nc file
            times = np.arange(1,ndays+1)
            cells = np.arange(1,ncells+1)
            coords = {'Time': times, 'nCells': cells}
            dims = ('Time', 'nCells')

            timeDaily_avg_RIO_ensembleMedian = xarray.DataArray(data=RIOmedian[dataOutIndex,:,:],
                coords=coords, dims=dims, name='timeDaily_avg_RIO_ensembleMedian')
            timeDaily_avg_RIO_ensemble5th = xarray.DataArray(data=RIO5th[dataOutIndex,:,:],
                coords=coords, dims=dims, name='timeDaily_avg_RIO_ensemble5th')
            timeDaily_avg_RIO_ensemble25th = xarray.DataArray(data=RIO25th[dataOutIndex,:,:],
                coords=coords, dims=dims, name='timeDaily_avg_RIO_ensemble25th')
            timeDaily_avg_RIO_ensemble75th = xarray.DataArray(data=RIO75th[dataOutIndex,:,:],
                coords=coords, dims=dims, name='timeDaily_avg_RIO_ensemble75th')
            timeDaily_avg_RIO_ensemble95th = xarray.DataArray(data=RIO95th[dataOutIndex,:,:],
                coords=coords, dims=dims, name='timeDaily_avg_RIO_ensemble95th')

            dataOut = xarray.Dataset({
                'timeDaily_avg_RIO_ensembleMedian':timeDaily_avg_RIO_ensembleMedian,
                'timeDaily_avg_RIO_ensemble5th':timeDaily_avg_RIO_ensemble5th,
                'timeDaily_avg_RIO_ensemble25th':timeDaily_avg_RIO_ensemble25th,
                'timeDaily_avg_RIO_ensemble75th':timeDaily_avg_RIO_ensemble75th,
                'timeDaily_avg_RIO_ensemble95th':timeDaily_avg_RIO_ensemble95th,
                })

            fileOut = writePath+'/'+prefix+vessel+'_EnsStats'+mpassiPrefix+year+'-'+month+'-01.nc'
            print('writing OUT file: ', fileOut)
            dataOut.to_netcdf(fileOut)

            dataOutIndex = dataOutIndex + 1

        ## write out all vessel types into single file (doesn't appear to save much time)
        # times = np.arange(1, ndays+1)
        # cells = np.arange(1, ncells+1)
        # 
        # dataOut = xarray.Dataset(
        #     data_vars=dict(
        #         RIO_ensembleMedian=(("vessel","Time","nCells"), RIOmedian),
        #         RIO_ensemble5th=(("vessel","Time","nCells"), RIO5th),
        #         RIO_ensemble95th=(("vessel","Time","nCells"), RIO95th),
        #     ),
        #     coords=dict(
        #         vessel=("vessel", np.array(classes, dtype="U")),
        #         Time=("Time", times),
        #         nCells=("nCells", cells),
        #     ),
        # )

        # fileOut = f"{writePath}/{prefix}EnsStats{mpassiPrefix}{year}-{month}-01.nc"
        # dataOut.to_netcdf(fileOut)

        ## clean up before next file is read in
        del dataIn, concentration, thickness, RIOArray, RIOmedian, RIO5th, RIO95th
        del timeDaily_avg_RIO_ensembleMedian, timeDaily_avg_RIO_ensemble5th, timeDaily_avg_RIO_ensemble95th
        del timeDaily_avg_RIO_ensemble25th, timeDaily_avg_RIO_ensemble75th
        del dataOut

elapsed = time.perf_counter() - t0
print(f"Loop time: {elapsed:.6f} s")
