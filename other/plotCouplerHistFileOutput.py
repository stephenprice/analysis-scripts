import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.interpolate import griddata

nc_file = './20240820.IGELM_MLI.ne30pg2_r05_IcoswISC30E3r5_gis20.chrysalis.gnu.DevelBudgets.cpl.hi.0001-01-02-00000.nc'
var_name = 'x2g_Flgl_qice'
lat_name = 'domg_lat'
lon_name = 'domg_lon'

# Open the NetCDF file
dataset = Dataset(nc_file)

# Read latitude, longitude, and variable data
#lats = dataset.variables[domg_lat][:]
#lons = dataset.variables[domg_lon][:]
#data = dataset.variables[x2g_Flgl_qice][:]

lats = dataset.variables[lat_name][:]
lons = dataset.variables[lon_name][:]
data = dataset.variables[var_name][:]

# Close the NetCDF file
dataset.close()

# reshape to vectors
shape = lats.shape
lats = lats.reshape(1,shape[2])
lons = lons.reshape(1,shape[2])
data = data.reshape(1,shape[2])

# Create a figure and axis with geographic projection
fig = plt.figure(figsize=(10, 7))
ax = plt.axes(projection=ccrs.PlateCarree())

# Add features to the map for better visualization
#ax.add_feature(cfeature.LAND)
#ax.add_feature(cfeature.OCEAN)
#ax.add_feature(cfeature.BORDERS, linestyle=':')
#ax.add_feature(cfeature.COASTLINE)

# Create a meshgrid for plotting
lon_grid, lat_grid = np.meshgrid(lons, lats)

# grid data onto mesh
dataI = np.griddata( coords, data, (lat_grid,lon_grid), method='cubic')

# Plot the data using pcolormesh
mesh = ax.pcolormesh(lon_grid, lat_grid, data, shading='auto', cmap='viridis')

# Add a colorbar
cbar = plt.colorbar(mesh, ax=ax, orientation='vertical')
cbar.set_label(var_name)

    # Add title
plt.title(f'Geographic Data: {var_name}')

# Show the plot
plt.show()

# Example usage

