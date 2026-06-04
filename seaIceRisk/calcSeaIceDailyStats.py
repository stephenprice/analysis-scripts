#!/usr/bin/env python3
"""
calcSeaIceDailyStats.py

Reads monthly MPAS-Sea Ice netCDF files for a set of ensemble members and
computes daily ensemble statistics (median, 5th and 95th percentile) for
ice area concentration and ice volume (effective thickness) across all members.
Results are written to monthly netCDF output files.

Usage:
    python calcSeaIceDailyStats.py
"""

import numpy as np
import xarray
import os


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# list of ensembles to read in and calc stats on
# ensembles = ['0051', '0101', '0111', '0121', '0131',
#              '0141', '0151', '0161', '0171', '0181',
#              '0191', '0201', '0211', '0221', '0231',
#              '0241', '0251', '0261', '0271', '0281',
#              '0301', '0311' ]
ensembles = ['0051', '0201']

# years = [ str(y) for y in range(2000, 2050) ]
# years = ['1980', '1981', '1982']
years = ['2000']
months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
# months = ['05']

# pathBase = '/lustre/scratch5/sprice/'
pathBase = '/Users/stephenprice/projects/CINS/ColdHarbor/'

prefix = 'v3.LR.historical_'
# dirs = '/archive/ice/hist/'
dirs = '/hist/'
mpassiPrefix = '.mpassi.hist.am.timeSeriesStatsDaily.'

# writePath = pathBase + 'data/v3.LR.historical_thkConc_ensembleStats'
writePath = pathBase + 'test-data-out/v3.LR.historical_thkConc_ensembleStats'

# ---------------------------------------------------------------------------
# Main processing loop
# ---------------------------------------------------------------------------

os.chdir(pathBase)
print('Output directory:', writePath)

# calc. ensemble stats by looping over years, months, and ensemble members
# (resulting outputs will be written to monthly netcdf files, as per the
# single ensemble datasets being read in)
for year in years:

    for month in months:

        index = 0

        for ensemble in ensembles:

            # build path and filename for this ensemble member / month
            fileBase = pathBase + prefix + ensemble + dirs
            file = prefix + ensemble + mpassiPrefix + year + '-' + month + '-01.nc'

            # read in file and relevant vars; use context manager to ensure
            # file handle is released after the block
            print('reading IN file: ', fileBase + file)
            with xarray.open_dataset(fileBase + file) as dataIn:

                concentration = dataIn.timeDaily_avg_iceAreaCell.values
                thickness = dataIn.timeDaily_avg_iceVolumeCell.values
                # can also read in other fields here (snow volume, u,v vel
                # components) if needed ...

                if index == 0:
                    # on first ensemble member: get dimensions, preserve the
                    # actual time coordinate from the input file, and allocate
                    # storage arrays
                    ndays = np.size(thickness, axis=0)
                    ncells = np.size(thickness, axis=1)
                    time_coord = dataIn.Time.values  # preserves actual datetime values
                    concentrationArray = np.zeros([len(ensembles), ndays, ncells])
                    thicknessArray = np.zeros([len(ensembles), ndays, ncells])

                concentrationArray[index, :, :] = concentration
                thicknessArray[index, :, :] = thickness
                index += 1

        # calc monthly stats (median, 5th, 95th percentile) across ensemble members
        concentrationMedian = np.median(concentrationArray, axis=0)
        concentration5th    = np.percentile(concentrationArray, 5,  axis=0)
        concentration95th   = np.percentile(concentrationArray, 95, axis=0)

        thicknessMedian = np.median(thicknessArray, axis=0)
        thickness5th    = np.percentile(thicknessArray, 5,  axis=0)
        thickness95th   = np.percentile(thicknessArray, 95, axis=0)

        # write monthly ensemble stats back out to (monthly) .nc file
        cells  = np.arange(1, ncells + 1)
        coords = {'Time': time_coord, 'nCells': cells}
        dims   = ('Time', 'nCells')

        timeDaily_avg_iceAreaCell_ensembleMedian = xarray.DataArray(
            data=concentrationMedian, coords=coords, dims=dims,
            name='timeDaily_avg_iceAreaCell_ensembleMedian')
        timeDaily_avg_iceAreaCell_ensemble5th = xarray.DataArray(
            data=concentration5th, coords=coords, dims=dims,
            name='timeDaily_avg_iceAreaCell_ensemble5th')
        timeDaily_avg_iceAreaCell_ensemble95th = xarray.DataArray(
            data=concentration95th, coords=coords, dims=dims,
            name='timeDaily_avg_iceAreaCell_ensemble95th')

        timeDaily_avg_iceVolumeCell_ensembleMedian = xarray.DataArray(
            data=thicknessMedian, coords=coords, dims=dims,
            name='timeDaily_avg_iceVolumeCell_ensembleMedian')
        timeDaily_avg_iceVolumeCell_ensemble5th = xarray.DataArray(
            data=thickness5th, coords=coords, dims=dims,
            name='timeDaily_avg_iceVolumeCell_ensemble5th')
        timeDaily_avg_iceVolumeCell_ensemble95th = xarray.DataArray(
            data=thickness95th, coords=coords, dims=dims,
            name='timeDaily_avg_iceVolumeCell_ensemble95th')

        dataOut = xarray.Dataset({
            'timeDaily_avg_iceAreaCell_ensembleMedian':   timeDaily_avg_iceAreaCell_ensembleMedian,
            'timeDaily_avg_iceAreaCell_ensemble5th':      timeDaily_avg_iceAreaCell_ensemble5th,
            'timeDaily_avg_iceAreaCell_ensemble95th':     timeDaily_avg_iceAreaCell_ensemble95th,
            'timeDaily_avg_iceVolumeCell_ensembleMedian': timeDaily_avg_iceVolumeCell_ensembleMedian,
            'timeDaily_avg_iceVolumeCell_ensemble5th':    timeDaily_avg_iceVolumeCell_ensemble5th,
            'timeDaily_avg_iceVolumeCell_ensemble95th':   timeDaily_avg_iceVolumeCell_ensemble95th,
        })

        fileOut = writePath + '/' + prefix + 'EnsStats' + mpassiPrefix + year + '-' + month + '-01.nc'
        print('writing OUT file: ', fileOut)
        dataOut.to_netcdf(fileOut)

        # clean up for next month
        del dataOut, fileOut
        del concentrationArray, thicknessArray, concentration, thickness, ndays, ncells
        del time_coord, cells, coords, dims
        del concentrationMedian, concentration5th, concentration95th
        del thicknessMedian, thickness5th, thickness95th
        del timeDaily_avg_iceAreaCell_ensembleMedian, \
            timeDaily_avg_iceAreaCell_ensemble5th, \
            timeDaily_avg_iceAreaCell_ensemble95th
        del timeDaily_avg_iceVolumeCell_ensembleMedian, \
            timeDaily_avg_iceVolumeCell_ensemble5th, \
            timeDaily_avg_iceVolumeCell_ensemble95th
