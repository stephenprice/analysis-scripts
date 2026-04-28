#
# Plots regional T,S profiles for ensemble members
# This breaks for more than one season or year/month conbination
#
from __future__ import absolute_import, division, print_function, \
    unicode_literals
import numpy as np
import xarray as xr
import matplotlib as mpl
#mpl.use('TkAgg') # avoid error associated w/ GUI support
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
import glob
import gsw


def haversine(lon1, lat1, lon2, lat2):
    # lon, lat should be in radians
    earthRadius = 6367.44 # km
    #earthRadius = 6371
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return 2 * earthRadius * np.arcsin(np.sqrt(a))

plotClimos = True
plotMonthly = False # not ready for prime time
if plotClimos==plotMonthly:
    raise ValueError('Variables plotClimos and plotMonthly cannot be identical')
plotPHCWOA = True # only works for monthly seasons for now (one season at a  time)
plotHighresMIP = True

#ensembleName = 'v3.LR.historical_'
ensembleName = 'HighresMIP'
#ensembleMemberNames = ['0051', '0101', '0151', '0201', '0251']
ensembleMemberNames = ['00']
colors = ['mediumblue', 'dodgerblue', 'deepskyblue', 'lightseagreen', 'teal'] # same length as ensembleMemberNames
#meshfile = '/lustre/scratch4/yellow/sprice/v3-lowRes-Historical/mesh/mpassi.IcoswISC30E3r5.rstFromG-chrysalis.20231121.nc'
meshfile = '/lustre/scratch4/yellow/sprice/v3-lowRes-Historical3/v3.LR.historical_0051/archive/rest/2025-01-01-00000/v3.LR.historical_0051.mpaso.rst.2025-01-01_00000.nc'

# Coordinates of point where to plot profiles (old list)
# Barents Sea:
#lonPoint = 45.5
#latPoint = 70.5
#pointTitle = 'Barents Sea South, 70N,37.5E'
#latPoint = 75
#pointTitle = 'Barents Sea Central, 75N,37.5E'
#latPoint = 80
#pointTitle = 'Barents Sea North, 80N,37.5E'
#latPoint = 75
#lonPoint = 27
#pointTitle = 'Barents Sea West, 75N,27E'
#lonPoint = 48
#pointTitle = 'Barents Sea East, 75N,48E'
#lonPoint = 35
#latPoint = 83
#pointTitle = 'Barents Sea Abyssal, 83N,35E'

# points to extract data from 
#coord_1 = [72.180254, 25.528254]
#coord_2 = [73.662038, 30.405605]
#coord_3 = [74.729154, 36.348002]
#coord_4 = [75.797133, 43.017310]
#coord_5 = [76.657394, 51.065280]
#coord_6 = [75.333924, 51.138969]
#coord_7 = [74.743833, 46.115986]
#coord_8 = [73.905472, 40.819224]
#coord_9 = [73.088944, 36.478335]
#coord_10 = [71.851702, 31.045571]
#coord_11 = [70.517602, 32.652502]
#coord_12 = [71.427471, 36.329257]
#coord_13 = [72.221160, 40.092688]
#coord_14 = [72.977877, 44.304213]
#coord_15 = [73.737446, 49.515908]
#coord_16 = [72.286483, 49.284013]
#coord_17 = [71.781488, 46.078221]
#coord_18 = [71.242994, 42.639499]
#coord_19 = [70.702002, 39.687635]
#coord_20= [69.891013, 36.447519]
#coord_21= [69.295241, 39.352997]
#coord_22= [69.954994, 42.287621]
#coord_23= [70.492933, 45.485927]
#coord_24= [71.002733, 48.590066]
#coord_25= [70.458097, 51.240771]
#coord_26= [69.813394, 47.698238]
#coord_27= [69.361516, 44.399446]
#coord_28= [75.000000, 25.700000]
#coord_29= [75.877452, 33.274699]
#coord_30 = [77.408830, 43.638144]
#coord_31 = [78.379693, 59.831029]

coord_1= [66.591467, 15.263025]
coord_2 = [62.457661, 13.023683]
coord_3 = [59.601244, 8.240783]


## relevant if plotClimos=True
climoyearStart = 1950
climoyearEnd = 1970
climoyearStart2 = 2000
climoyearEnd2 = 2014

# seasons options: '01'-'12', 'ANN', 'JFM', 'JAS', 'MAJ', 'OND'
# (depending on what was set in mpas-analysis)
seasons = ['02', '05', '08', '11']
#seasons = ['ANN']
#seasons = ['JFM', 'JAS']
modelClimodir1 = f'/lustre/scratch4/turquoise/sprice/v3-lowRes-ensemble-reduced2/{ensembleName}'
modelClimodir2 = f'mpas-analysis/Years{climoyearStart}-{climoyearEnd}'

# relevant if plotMonthly=True
years = [1950]
months = [9]
modeldir1 = f'/pscratch/sd/m/milena/e3sm_scratch/pm-cpu/{ensembleName}'
modeldir2 = f'archive/ocn/hist'

# relevant if plotPHCWOA=True
PHCfilename = '/lustre/scratch4/yellow/sprice/v3-lowRes-Historical3/analysis/ClimosObs/phc3.0_monthly_accessed08-08-2019.nc'
WOAfilename = '/lustre/scratch4/yellow/sprice/v3-lowRes-Historical3/analysis/ClimosObs/woa18_decav_04_TS_mon.nc'

# relevant if plotHighresMIP=True
HighresMIPdir = '/lustre/scratch4/yellow/sprice/v3-lowRes-Historical3/analysis/ClimosHighResMip/CESM1-CAM5-SE-HR/hist-1950/ncclimoFiles'
HighresMIP2dir = '/lustre/scratch4/yellow/sprice/v3-lowRes-Historical3/analysis/ClimosHighResMip/CESM1-CAM5-SE-HR/highres-future/ncclimoFiles'

# highresMIP bias corrections
biasCorrectT = -3.0
biasCorrectS = -0.4

# the two chunks below have been moved into the loop over points

#figdir = f'./TSprofiles/{ensembleName}'
#if not os.path.isdir(figdir):
#    os.makedirs(figdir)

#outdir0 = f'./TSprofiles_data'
#if not os.path.isdir(outdir0):
#    os.makedirs(outdir0)

figsize = (10, 15)
figdpi = 150
fontsize_smallLabels = 18
fontsize_labels = 20
fontsize_titles = 22
legend_properties = {'size':fontsize_smallLabels, 'weight':'bold'}

nEnsembles = len(ensembleMemberNames)

############### loop over points in list ##############

#for i in range(1, 32):
for i in range(1, 3):

    var_name = f"coord_{i}"

    pointTitle = 'Barents Sea, point ', str(i)

    figdir = f'./TSprofiles_{i}/{ensembleName}'
    if not os.path.isdir(figdir):
        os.makedirs(figdir)

    outdir0 = f'./TSprofiles_data_{i}'
    if not os.path.isdir(outdir0):
        os.makedirs(outdir0)

    # Read in relevant global mesh information
    if os.path.exists(meshfile):
        dsMesh = xr.open_dataset(meshfile)
    else:
        raise IOError(f'MPAS restart/mesh file {meshfile} not found')
    depth = dsMesh.refBottomDepth
    # Identify index of selected ocean cell, by computing the minimum
    # of the spherical distance between all points and lonPoint,latPoint
    coords = globals()[var_name]
    lonPoint = coords[1]
    latPoint = coords[0]

    nCells = dsMesh.dims['nCells']
    lonCell = dsMesh.lonCell
    latCell = dsMesh.latCell
    spherDist = haversine(lonCell, latCell, lonPoint*np.pi/180, latPoint*np.pi/180)
    indices = xr.DataArray(data=np.arange(nCells).astype(int), dims='nCells')
    iCell = indices.where(spherDist==np.min(spherDist), drop=True).values.astype(int)[0]

    lon_icell = lonCell.values[iCell]*180/np.pi
    lat_icell = latCell.values[iCell]*180/np.pi
    print(lonPoint, latPoint)
    print(lon_icell, lat_icell)
    pres = gsw.conversions.p_from_z(-depth, lat_icell)
    nLevels = dsMesh.dims['nVertLevels']
    maxLevelCell = dsMesh.maxLevelCell.isel(nCells=iCell)
    vertIndex = xr.DataArray.from_dict({'dims': ('nVertLevels',),
                                    'data': np.arange(nLevels)})
    vertMask = vertIndex < maxLevelCell
    
    if plotPHCWOA is True:
        # Read in PHC climo
        dsPHC = xr.open_dataset(PHCfilename, decode_times=False)
        # Identify index of selected ocean cell, by computing the minimum
        # of the spherical distance between all points and lonPoint,latPoint
        latPHC = dsPHC.lat.values
        lonPHC = dsPHC.lon.values
        [x, y] = np.meshgrid(lonPHC, latPHC)
        if lonPoint<0:
            spherDist = haversine(x*np.pi/180, y*np.pi/180, (lonPoint+360)*np.pi/180, latPoint*np.pi/180)
        else:
            spherDist = haversine(x*np.pi/180, y*np.pi/180, lonPoint*np.pi/180, latPoint*np.pi/180)
        x = x[np.where(spherDist==np.min(spherDist))][0]
        y = y[np.where(spherDist==np.min(spherDist))][0]
        dsPHC = dsPHC.sel(lat=y, lon=x)
        depthPHC = dsPHC.depth
        presPHC = gsw.conversions.p_from_z(-depthPHC, y)
    
        # Read in WOA climo
        dsWOA = xr.open_dataset(WOAfilename)
        # Identify index of selected ocean cell, by computing the minimum
        # of the spherical distance between all points and lonPoint,latPoint
        latWOA = dsWOA.lat.values
        lonWOA = dsWOA.lon.values
        [x, y] = np.meshgrid(lonWOA, latWOA)
        spherDist = haversine(x*np.pi/180, y*np.pi/180, lonPoint*np.pi/180, latPoint*np.pi/180)
        x = x[np.where(spherDist==np.min(spherDist))][0]
        y = y[np.where(spherDist==np.min(spherDist))][0]
        dsWOA = dsWOA.sel(lat=y, lon=x)
        depthWOA = dsWOA.depth
        presWOA = gsw.conversions.p_from_z(-depthWOA, y)
    
    if plotHighresMIP is True:
        # Read in data
        Tfiles = []
        Sfiles = []
        for im in range(1, 13):
            Tfiles.append(f'{HighresMIPdir}/thetao_Omon_CESM1-CAM5-SE-HR_hist-1950_r1i1p1f1_gn_{im:02d}_{climoyearStart:04d}{im:02d}_{climoyearEnd:04d}{im:02d}_climo.nc')
            Sfiles.append(f'{HighresMIPdir}/so_Omon_CESM1-CAM5-SE-HR_hist-1950_r1i1p1f1_gn_{im:02d}_{climoyearStart:04d}{im:02d}_{climoyearEnd:04d}{im:02d}_climo.nc')
        dsHighresMIPtemp = xr.open_mfdataset(Tfiles, combine='nested', concat_dim='time', decode_times=False)
        dsHighresMIPsalt = xr.open_mfdataset(Sfiles, combine='nested', concat_dim='time', decode_times=False)
        # Identify index of selected ocean cell, by computing the minimum
        # of the spherical distance between all points and lonPoint,latPoint
        lat = dsHighresMIPtemp.coords['lat'].values
        lon = dsHighresMIPtemp.coords['lon'].values
        if lonPoint<0:
            spherDist = haversine(lon*np.pi/180, lat*np.pi/180, (lonPoint+360)*np.pi/180, latPoint*np.pi/180)
        else:
            spherDist = haversine(lon*np.pi/180, lat*np.pi/180, lonPoint*np.pi/180, latPoint*np.pi/180)
        [nlat, nlon] = np.argwhere(spherDist==np.min(spherDist))[0]
        dsHighresMIPtemp = dsHighresMIPtemp.sel(nlat=nlat, nlon=nlon)
        dsHighresMIPsalt = dsHighresMIPsalt.sel(nlat=nlat, nlon=nlon)
        HighresMIPdepth = 1e-2 * dsHighresMIPtemp['lev']
        HighresMIPpres = gsw.conversions.p_from_z(-HighresMIPdepth, lat[nlat, nlon])
        #
        Tfiles = []
        Sfiles = []
        for im in range(1, 13):
            Tfiles.append(f'{HighresMIPdir}/thetao_Omon_CESM1-CAM5-SE-HR_hist-1950_r1i1p1f1_gn_{im:02d}_{climoyearStart2:04d}{im:02d}_{climoyearEnd2:04d}{im:02d}_climo.nc')
            Sfiles.append(f'{HighresMIPdir}/so_Omon_CESM1-CAM5-SE-HR_hist-1950_r1i1p1f1_gn_{im:02d}_{climoyearStart2:04d}{im:02d}_{climoyearEnd2:04d}{im:02d}_climo.nc')
        dsHighresMIPtemp2 = xr.open_mfdataset(Tfiles, combine='nested', concat_dim='time', decode_times=False)
        dsHighresMIPsalt2 = xr.open_mfdataset(Sfiles, combine='nested', concat_dim='time', decode_times=False)
        # Identify index of selected ocean cell, by computing the minimum
        # of the spherical distance between all points and lonPoint,latPoint
        lat = dsHighresMIPtemp2.coords['lat'].values
        lon = dsHighresMIPtemp2.coords['lon'].values
        if lonPoint<0:
            spherDist = haversine(lon*np.pi/180, lat*np.pi/180, (lonPoint+360)*np.pi/180, latPoint*np.pi/180)
        else:
            spherDist = haversine(lon*np.pi/180, lat*np.pi/180, lonPoint*np.pi/180, latPoint*np.pi/180)
        [nlat, nlon] = np.argwhere(spherDist==np.min(spherDist))[0]
        dsHighresMIPtemp2 = dsHighresMIPtemp2.sel(nlat=nlat, nlon=nlon)
        dsHighresMIPsalt2 = dsHighresMIPsalt2.sel(nlat=nlat, nlon=nlon)
        HighresMIPdepth2 = 1e-2 * dsHighresMIPtemp2['lev']
        HighresMIPpres2 = gsw.conversions.p_from_z(-HighresMIPdepth2, lat[nlat, nlon])
        #
        Tfiles = []
        Sfiles = []
        for im in range(1, 13):
            Tfiles.append(f'{HighresMIP2dir}/thetao_Omon_CESM1-CAM5-SE-HR_highres-future_r1i1p1f1_gn_{im:02d}_2031{im:02d}_2050{im:02d}_climo.nc')
            Sfiles.append(f'{HighresMIP2dir}/so_Omon_CESM1-CAM5-SE-HR_highres-future_r1i1p1f1_gn_{im:02d}_2031{im:02d}_2050{im:02d}_climo.nc')
        dsHighresMIPtemp3 = xr.open_mfdataset(Tfiles, combine='nested', concat_dim='time', decode_times=False)
        dsHighresMIPsalt3 = xr.open_mfdataset(Sfiles, combine='nested', concat_dim='time', decode_times=False)
        # Identify index of selected ocean cell, by computing the minimum
        # of the spherical distance between all points and lonPoint,latPoint
        lat = dsHighresMIPtemp3.coords['lat'].values
        lon = dsHighresMIPtemp3.coords['lon'].values
        if lonPoint<0:
            spherDist = haversine(lon*np.pi/180, lat*np.pi/180, (lonPoint+360)*np.pi/180, latPoint*np.pi/180)
        else:
            spherDist = haversine(lon*np.pi/180, lat*np.pi/180, lonPoint*np.pi/180, latPoint*np.pi/180)
        [nlat, nlon] = np.argwhere(spherDist==np.min(spherDist))[0]
        dsHighresMIPtemp3 = dsHighresMIPtemp3.sel(nlat=nlat, nlon=nlon)
        dsHighresMIPsalt3 = dsHighresMIPsalt3.sel(nlat=nlat, nlon=nlon)
        HighresMIPdepth3 = 1e-2 * dsHighresMIPtemp3['lev']
        HighresMIPpres3 = gsw.conversions.p_from_z(-HighresMIPdepth3, lat[nlat, nlon])

    if plotClimos is True:
        for season in seasons:
            # Initialize figure and axis objects
            fig_Tprofile = plt.figure(figsize=figsize, dpi=figdpi)
            ax_Tprofile = fig_Tprofile.add_subplot()
            for tick in ax_Tprofile.xaxis.get_ticklabels():
                tick.set_fontsize(fontsize_smallLabels)
                tick.set_weight('bold')
            for tick in ax_Tprofile.yaxis.get_ticklabels():
                tick.set_fontsize(fontsize_smallLabels)
                tick.set_weight('bold')
            ax_Tprofile.yaxis.get_offset_text().set_fontsize(fontsize_smallLabels)
            ax_Tprofile.yaxis.get_offset_text().set_weight('bold')
            #
            fig_Sprofile = plt.figure(figsize=figsize, dpi=figdpi)
            ax_Sprofile = fig_Sprofile.add_subplot()
            for tick in ax_Sprofile.xaxis.get_ticklabels():
                tick.set_fontsize(fontsize_smallLabels)
                tick.set_weight('bold')
            for tick in ax_Sprofile.yaxis.get_ticklabels():
                tick.set_fontsize(fontsize_smallLabels)
                tick.set_weight('bold')
            ax_Sprofile.yaxis.get_offset_text().set_fontsize(fontsize_smallLabels)
            ax_Sprofile.yaxis.get_offset_text().set_weight('bold')
            #
            fig_Cprofile = plt.figure(figsize=figsize, dpi=figdpi)
            ax_Cprofile = fig_Cprofile.add_subplot()
            for tick in ax_Cprofile.xaxis.get_ticklabels():
                tick.set_fontsize(fontsize_smallLabels)
                tick.set_weight('bold')
            for tick in ax_Cprofile.yaxis.get_ticklabels():
                tick.set_fontsize(fontsize_smallLabels)
                tick.set_weight('bold')
            ax_Cprofile.yaxis.get_offset_text().set_fontsize(fontsize_smallLabels)
            ax_Cprofile.yaxis.get_offset_text().set_weight('bold')

            Tfigtitle = f'Temperature ({pointTitle})\n{season} - years {climoyearStart:04d}-{climoyearEnd:04d}'
            Tfigfile = f'{figdir}/Tprofile_icell{iCell:d}_{ensembleName}_{season}_years{climoyearStart:04d}-{climoyearEnd:04d}.png'
            Sfigtitle = f'Salinity ({pointTitle})\n{season} - years {climoyearStart:04d}-{climoyearEnd:04d}'
            Sfigfile = f'{figdir}/Sprofile_icell{iCell:d}_{ensembleName}_{season}_years{climoyearStart:04d}-{climoyearEnd:04d}.png'
            Cfigtitle = f'Sound speed ({pointTitle})\n{season} - years {climoyearStart:04d}-{climoyearEnd:04d}'
            Cfigfile = f'{figdir}/Cprofile_icell{iCell:d}_{ensembleName}_{season}_years{climoyearStart:04d}-{climoyearEnd:04d}.png'

            ax_Tprofile.set_xlabel('Temperature ($^\circ$C)', fontsize=fontsize_labels, fontweight='bold')
            ax_Tprofile.set_ylabel('Depth (m)', fontsize=fontsize_labels, fontweight='bold')
            ax_Tprofile.set_title(Tfigtitle, fontsize=fontsize_titles, fontweight='bold')
            #ax_Tprofile.set_xlim(-1.85, 1.8)
            ax_Tprofile.set_ylim(-depth[maxLevelCell.values], 0)
            #ax_Tprofile.set_ylim(-800, 0)
            #
            ax_Sprofile.set_xlabel('Salinity (psu)', fontsize=fontsize_labels, fontweight='bold')
            ax_Sprofile.set_ylabel('Depth (m)', fontsize=fontsize_labels, fontweight='bold')
            ax_Sprofile.set_title(Sfigtitle, fontsize=fontsize_titles, fontweight='bold')
            #ax_Sprofile.set_xlim(27.8, 35)
            ax_Sprofile.set_ylim(-depth[maxLevelCell.values], 0)
            #ax_Sprofile.set_ylim(-800, 0)
            #
            ax_Cprofile.set_xlabel('C (m/s)', fontsize=fontsize_labels, fontweight='bold')
            ax_Cprofile.set_ylabel('Depth (m)', fontsize=fontsize_labels, fontweight='bold')
            ax_Cprofile.set_title(Cfigtitle, fontsize=fontsize_titles, fontweight='bold')
            #ax_Cprofile.set_xlim(1430., 1470.)
            ax_Cprofile.set_ylim(-depth[maxLevelCell.values], 0)
            #ax_Cprofile.set_ylim(-800, 0)

            for i in range(nEnsembles):
                ensembleMemberName = ensembleMemberNames[i]
                print(f'\nProcessing ensemble member {ensembleMemberName}, season {season}...')

#                modelfile = f'{modelClimodir1}{ensembleMemberName}/{modelClimodir2}/mpaso_{season}_{climoyearStart:04d}{season}_{climoyearEnd:04d}{season}_climo.nc'
#
#                dsIn = xr.open_dataset(modelfile).isel(Time=0, nCells=iCell)
#                dsIn = dsIn.where(vertMask)
    #            # Drop all variables but T and S, and mask bathymetry
#                allvars = dsIn.data_vars.keys()
#                dropvars = set(allvars) - set(['timeMonthly_avg_activeTracers_temperature',
#                                           'timeMonthly_avg_activeTracers_salinity'])
#                dsIn = dsIn.drop(dropvars)
#
#                Tprofile = dsIn.timeMonthly_avg_activeTracers_temperature.values
    #            Sprofile = dsIn.timeMonthly_avg_activeTracers_salinity.values
#                SA = gsw.conversions.SA_from_SP(Sprofile, pres, lon_icell, lat_icell)
#                CT = gsw.conversions.CT_from_pt(SA, Tprofile)
#                #sigma0profile = gsw.density.sigma0(SA, CT)
#                soundspeed = gsw.sound_speed(SA, CT, pres)
#
    #            ax_Tprofile.plot(Tprofile[::-1], -depth[::-1], '-', color=colors[i], linewidth=3, label=f'{ensembleMemberName}')
#                ax_Sprofile.plot(Sprofile[::-1], -depth[::-1], '-', color=colors[i], linewidth=3, label=f'{ensembleMemberName}')
#                ax_Cprofile.plot(soundspeed[::-1], -depth[::-1], '-', color=colors[i], linewidth=3, label=f'{ensembleMemberName}')

#                # Write to file
    #            outdir = f'{outdir0}/{ensembleName}/{ensembleMemberName}'
#                if not os.path.isdir(outdir):
#                    os.makedirs(outdir)
#                outfile = f'{outdir}/icell{iCell:d}_profiles_{ensembleName}{ensembleMemberName}_{season}_years{climoyearStart:04d}-{climoyearEnd:04d}.nc'
#                dsOut = xr.Dataset()
    #            dsOut['Tprofile'] = Tprofile
#                dsOut['Tprofile'].attrs['units'] = 'degC'
#                dsOut['Tprofile'].attrs['long_name'] = 'Potential temperature'
#                dsOut['Sprofile'] = Sprofile
#                dsOut['Sprofile'].attrs['units'] = 'psu'
    #            dsOut['Sprofile'].attrs['long_name'] = 'Salinity'
#                dsOut['CTprofile'] = CT
#                dsOut['CTprofile'].attrs['units'] = 'degC'
#                dsOut['CTprofile'].attrs['long_name'] = 'Conservative temperature'
#                dsOut['SAprofile'] = SA
    #            dsOut['SAprofile'].attrs['units'] = 'psu'
#                dsOut['SAprofile'].attrs['long_name'] = 'Absolute salinity'
#                dsOut['Cprofile'] = soundspeed
#                dsOut['Cprofile'].attrs['units'] = 'm/s'
#                dsOut['Cprofile'].attrs['long_name'] = 'Sound speed (computed with python gsw package)'
    #            dsOut['depth'] = depth
#                dsOut['depth'].attrs['units'] = 'm'
#                dsOut['depth'].attrs['long_name'] = 'depth levels'
#                dsOut['lon'] = lon_icell
#                dsOut['lon'].attrs['units'] = 'degrees_east'
    #            dsOut['lon'].attrs['long_name'] = 'point longitude'
#                dsOut['lat'] = lat_icell
#                dsOut['lat'].attrs['units'] = 'degrees_north'
#                dsOut['lat'].attrs['long_name'] = 'point latitude'
#                dsOut.to_netcdf(outfile)

            if plotPHCWOA is True:
                dsPHC_monthlyClimo = dsPHC.isel(time=int(season)-1)
                SA = gsw.conversions.SA_from_SP(dsPHC_monthlyClimo['salt'].values, presPHC, x, y)
                CT = gsw.conversions.CT_from_pt(SA, dsPHC_monthlyClimo['temp'].values)
                soundspeedPHC = gsw.sound_speed(SA, CT, presPHC)

                dsWOA_monthlyClimo = dsWOA.isel(month=int(season)-1)
                SA = gsw.conversions.SA_from_SP(dsWOA_monthlyClimo['s_an'].values, presWOA, x, y)
                CT = gsw.conversions.CT_from_pt(SA, dsWOA_monthlyClimo['t_an'].values)
                soundspeedWOA = gsw.sound_speed(SA, CT, presWOA)

                ax_Tprofile.plot(dsPHC_monthlyClimo['temp'][::-1], -depthPHC[::-1], '-', color='mediumvioletred',
                             linewidth=3, label='PHC climatology')
                ax_Sprofile.plot(dsPHC_monthlyClimo['salt'][::-1], -depthPHC[::-1], '-', color='mediumvioletred',
                             linewidth=3, label='PHC climatology')
                ax_Cprofile.plot(soundspeedPHC[::-1], -depthPHC[::-1], '-', color='mediumvioletred',
                             linewidth=3, label='PHC climatology')

                ax_Tprofile.plot(dsWOA_monthlyClimo['t_an'][::-1], -depthWOA[::-1], '-', color='salmon',
                             linewidth=3, label='WOA climatology')
                ax_Sprofile.plot(dsWOA_monthlyClimo['s_an'][::-1], -depthWOA[::-1], '-', color='salmon',
                             linewidth=3, label='WOA climatology')
                ax_Cprofile.plot(soundspeedWOA[::-1], -depthWOA[::-1], '-', color='salmon',
                             linewidth=3, label='WOA climatology')

            if plotHighresMIP is True:
                #HighresMIPtemp = dsHighresMIPtemp['thetao'].isel(time=int(season)-1)
                HighresMIPtemp = dsHighresMIPtemp['thetao'].isel(time=int(season)-1) + biasCorrectT
                HighresMIPsalt = dsHighresMIPsalt['so'].isel(time=int(season)-1) + biasCorrectS
                SA = gsw.conversions.SA_from_SP(HighresMIPsalt.values, HighresMIPpres, x, y)
                CT = gsw.conversions.CT_from_pt(SA, HighresMIPtemp.values)
                soundspeed = gsw.sound_speed(SA, CT, HighresMIPpres)

                ax_Tprofile.plot(HighresMIPtemp[::-1], -HighresMIPdepth[::-1], '-', color='gold',
                             linewidth=3, label='HighresMIP 1950-1970')
                ax_Sprofile.plot(HighresMIPsalt[::-1], -HighresMIPdepth[::-1], '-', color='gold',
                             linewidth=3, label='HighresMIP 1950-1970')
                ax_Cprofile.plot(soundspeed[::-1], -HighresMIPdepth[::-1], '-', color='gold',
                             linewidth=3, label='HighresMIP 1950-1970')

                # Write to file
                outdir = f'{outdir0}/HighresMIP/hist-1950'
                if not os.path.isdir(outdir):
                    os.makedirs(outdir)
                outfile = f'{outdir}/icell{iCell:d}_profiles_HighresMIP_hist-{season}_years{climoyearStart:04d}-{climoyearEnd:04d}.nc'
                dsOut = xr.Dataset()
                dsOut['Tprofile'] = HighresMIPtemp
                dsOut['Tprofile'].attrs['units'] = 'degC'
                dsOut['Tprofile'].attrs['long_name'] = 'Potential temperature'
                dsOut['Sprofile'] = HighresMIPsalt
                dsOut['Sprofile'].attrs['units'] = 'psu'
                dsOut['Sprofile'].attrs['long_name'] = 'Salinity'
                dsOut['CTprofile'] = CT
                dsOut['CTprofile'].attrs['units'] = 'degC'
                dsOut['CTprofile'].attrs['long_name'] = 'Conservative temperature'
                dsOut['SAprofile'] = SA
                dsOut['SAprofile'].attrs['units'] = 'psu'
                dsOut['SAprofile'].attrs['long_name'] = 'Absolute salinity'
                dsOut['Cprofile'] = soundspeed
                dsOut['Cprofile'].attrs['units'] = 'm/s'
                dsOut['Cprofile'].attrs['long_name'] = 'Sound speed (computed with python gsw package)'
                dsOut['depth'] = HighresMIPdepth
                dsOut['depth'].attrs['units'] = 'm'
                dsOut['depth'].attrs['long_name'] = 'depth levels'
                dsOut['lon'] = x
                dsOut['lon'].attrs['units'] = 'degrees_east'
                dsOut['lon'].attrs['long_name'] = 'point longitude'
                dsOut['lat'] = y
                dsOut['lat'].attrs['units'] = 'degrees_north'
                dsOut['lat'].attrs['long_name'] = 'point latitude'
                dsOut.to_netcdf(outfile)
                #
                #HighresMIPtemp2 = dsHighresMIPtemp2['thetao'].isel(time=int(season)-1)
                HighresMIPtemp2 = dsHighresMIPtemp2['thetao'].isel(time=int(season)-1) + biasCorrectT
                HighresMIPsalt2 = dsHighresMIPsalt2['so'].isel(time=int(season)-1) + biasCorrectS
                SA = gsw.conversions.SA_from_SP(HighresMIPsalt2.values, HighresMIPpres2, x, y)
                CT = gsw.conversions.CT_from_pt(SA, HighresMIPtemp2.values)
                soundspeed = gsw.sound_speed(SA, CT, HighresMIPpres2)
    
                ax_Tprofile.plot(HighresMIPtemp2[::-1], -HighresMIPdepth2[::-1], '-', color='darkgoldenrod',
                             linewidth=3, label='HighresMIP 2000-2014')
                ax_Sprofile.plot(HighresMIPsalt2[::-1], -HighresMIPdepth2[::-1], '-', color='darkgoldenrod',
                             linewidth=3, label='HighresMIP 2000-2014')
                ax_Cprofile.plot(soundspeed[::-1], -HighresMIPdepth2[::-1], '-', color='darkgoldenrod',
                             linewidth=3, label='HighresMIP 2000-2014')
                #
                outdir = f'{outdir0}/HighresMIP/hist-1950'
                if not os.path.isdir(outdir):
                    os.makedirs(outdir)
                outfile = f'{outdir}/icell{iCell:d}_profiles_HighresMIP_hist-{season}_years{climoyearStart2:04d}-{climoyearEnd2:04d}.nc'
                dsOut = xr.Dataset()
                dsOut['Tprofile'] = HighresMIPtemp2
                dsOut['Tprofile'].attrs['units'] = 'degC'
                dsOut['Tprofile'].attrs['long_name'] = 'Potential temperature'
                dsOut['Sprofile'] = HighresMIPsalt2
                dsOut['Sprofile'].attrs['units'] = 'psu'
                dsOut['Sprofile'].attrs['long_name'] = 'Salinity'
                dsOut['CTprofile'] = CT
                dsOut['CTprofile'].attrs['units'] = 'degC'
                dsOut['CTprofile'].attrs['long_name'] = 'Conservative temperature'
                dsOut['SAprofile'] = SA
                dsOut['SAprofile'].attrs['units'] = 'psu'
                dsOut['SAprofile'].attrs['long_name'] = 'Absolute salinity'
                dsOut['Cprofile'] = soundspeed
                dsOut['Cprofile'].attrs['units'] = 'm/s'
                dsOut['Cprofile'].attrs['long_name'] = 'Sound speed (computed with python gsw package)'
                dsOut['depth'] = HighresMIPdepth2
                dsOut['depth'].attrs['units'] = 'm'
                dsOut['depth'].attrs['long_name'] = 'depth levels'
                dsOut['lon'] = x
                dsOut['lon'].attrs['units'] = 'degrees_east'
                dsOut['lon'].attrs['long_name'] = 'point longitude'
                dsOut['lat'] = y
                dsOut['lat'].attrs['units'] = 'degrees_north'
                dsOut['lat'].attrs['long_name'] = 'point latitude'
                dsOut.to_netcdf(outfile)
                #
                #HighresMIPtemp3 = dsHighresMIPtemp3['thetao'].isel(time=int(season)-1)
                HighresMIPtemp3 = dsHighresMIPtemp3['thetao'].isel(time=int(season)-1) + biasCorrectT
                HighresMIPsalt3 = dsHighresMIPsalt3['so'].isel(time=int(season)-1) + biasCorrectS  
                SA = gsw.conversions.SA_from_SP(HighresMIPsalt3.values, HighresMIPpres3, x, y)
                CT = gsw.conversions.CT_from_pt(SA, HighresMIPtemp3.values)
                soundspeed = gsw.sound_speed(SA, CT, HighresMIPpres3)

                ax_Tprofile.plot(HighresMIPtemp3[::-1], -HighresMIPdepth3[::-1], '-', color='brown',
                             linewidth=3, label='HighresMIP 2031-2050')
                ax_Sprofile.plot(HighresMIPsalt3[::-1], -HighresMIPdepth3[::-1], '-', color='brown',
                             linewidth=3, label='HighresMIP 2031-2050')
                ax_Cprofile.plot(soundspeed[::-1], -HighresMIPdepth3[::-1], '-', color='brown',
                             linewidth=3, label='HighresMIP 2031-2050')
                #
                # Write to file
                outdir = f'{outdir0}/HighresMIP/highres-future'
                if not os.path.isdir(outdir):
                    os.makedirs(outdir)
                outfile = f'{outdir}/icell{iCell:d}_profiles_HighresMIP_highres-future_{season}_years2031-2050.nc'
                dsOut = xr.Dataset()
                dsOut['Tprofile'] = HighresMIPtemp3
                dsOut['Tprofile'].attrs['units'] = 'degC'
                dsOut['Tprofile'].attrs['long_name'] = 'Potential temperature'
                dsOut['Sprofile'] = HighresMIPsalt3
                dsOut['Sprofile'].attrs['units'] = 'psu'
                dsOut['Sprofile'].attrs['long_name'] = 'Salinity'
                dsOut['CTprofile'] = CT
                dsOut['CTprofile'].attrs['units'] = 'degC'
                dsOut['CTprofile'].attrs['long_name'] = 'Conservative temperature'
                dsOut['SAprofile'] = SA
                dsOut['SAprofile'].attrs['units'] = 'psu'
                dsOut['SAprofile'].attrs['long_name'] = 'Absolute salinity'
                dsOut['Cprofile'] = soundspeed
                dsOut['Cprofile'].attrs['units'] = 'm/s'
                dsOut['Cprofile'].attrs['long_name'] = 'Sound speed (computed with python gsw package)'
                dsOut['depth'] = HighresMIPdepth3
                dsOut['depth'].attrs['units'] = 'm'
                dsOut['depth'].attrs['long_name'] = 'depth levels'
                dsOut['lon'] = x
                dsOut['lon'].attrs['units'] = 'degrees_east'
                dsOut['lon'].attrs['long_name'] = 'point longitude'
                dsOut['lat'] = y
                dsOut['lat'].attrs['units'] = 'degrees_north'
                dsOut['lat'].attrs['long_name'] = 'point latitude'
                dsOut.to_netcdf(outfile)
    
            #ax_Tprofile.legend(prop=legend_properties)
            ax_Tprofile.legend(prop=legend_properties, loc='lower left', bbox_to_anchor=(1, 0.5))
            ax_Tprofile.grid(visible=True, which='both')
            fig_Tprofile.savefig(Tfigfile, bbox_inches='tight')
            plt.close(fig_Tprofile)

            #ax_Sprofile.legend(prop=legend_properties)
            ax_Sprofile.legend(prop=legend_properties, loc='lower left', bbox_to_anchor=(1, 0.5))
            ax_Sprofile.grid(visible=True, which='both')
            fig_Sprofile.savefig(Sfigfile, bbox_inches='tight')
            plt.close(fig_Sprofile)

            #ax_Cprofile.legend(prop=legend_properties)
            ax_Cprofile.legend(prop=legend_properties, loc='lower left', bbox_to_anchor=(1, 0.5))
            ax_Cprofile.grid(visible=True, which='both')
            fig_Cprofile.savefig(Cfigfile, bbox_inches='tight')
            plt.close(fig_Cprofile)
