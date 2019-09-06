import numpy as np
import sys
import xarray as xr
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')
from cartopy import config
import cartopy.crs as ccrs
import cartopy.feature as cfeature
# and a new matplotlib routine
import matplotlib.path as mpath
import ecco_v4_py as ecco
import pdb

def extractArctic(arr,depth):
    arr = arr[105300*depth:105300*(depth+1)]
    arr = arr.byteswap().reshape([105300,1])

    arr = xr.DataArray(arr[48600:56700,0].reshape([90,90]),coords=[np.arange(0,90,1),np.arange(0,90,1)],
                            dims=['j','i'])
    return arr

varnames = ["diffkr","kapredi","kapgm"]
for var in varnames:
    for depth in range(0,50,5):
        geoflx = np.fromfile('ecco/mixingdata/'+var+'.bin', dtype=np.float32)
        geoflx06 = extractArctic(geoflx,depth)
        geoflx = np.fromfile('ecco/mixingdata/'+var+'.data', dtype=np.float32)
        geoflx06 = extractArctic(geoflx,depth)+geoflx06

        ecco_grid = xr.open_dataset('ecco/TILEDATA/ARCTICGRID.nc')
        ecco_nvel = xr.open_dataset('ecco/TILEDATA/NVELARCTIC.nc',decode_times=False)
        tile_num=7
        lons = ecco_grid.XC
        lats = ecco_grid.YC

        tile_to_plot = geoflx06
        for i in range((ecco_grid.hFacC.shape[1])):
            for j in range((ecco_grid.hFacC.shape[2])):
                if not ecco_grid.hFacC[depth][i][j]:
                    tile_to_plot[i][j] = np.nan
        plt.figure(figsize=(8,6), dpi= 90)


        # Make a new projection, time of class "NorthPolarStereo"
        ax = plt.axes(projection=ccrs.NorthPolarStereo(true_scale_latitude=70))

        # here is here you tell Cartopy that the projection
        # of your 'x' and 'y' are geographic (lons and lats)
        # and that you want to transform those lats and lons
        #into 'x' and 'y' in the projection
        plt.pcolormesh(lons, lats, tile_to_plot,
                       transform=ccrs.PlateCarree());

        # plot land
        ax.add_feature(cfeature.LAND)
        ax.gridlines()
        ax.coastlines()
        #plt.colorbar()

        # Limit the map to -60 degrees latitude and below.
        ax.set_extent([-180, 180, 90, 60], ccrs.PlateCarree())

        # Compute a circle in axes coordinates, which we can use as a boundary
        # for the map. We can pan/zoom as much as we like - the boundary will be
        # permanently circular.
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        meterdepth = str(int(ecco_nvel.dep[depth]))
        ax.set_title(var + ": " + meterdepth)
        plt.savefig("refpics/surfaces/eccomixnoswap/"+var+"/"+meterdepth+".png")
        plt.close()
