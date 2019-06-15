import csv
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from numpy import linspace
from numpy import meshgrid
from netCDF4 import Dataset
import pickle



d = Dataset("data/ver1_netcdf_geo.nc")
z = d.variables["z"]
spacing = d.variables["spacing"][0]
startlon = d.variables["x_range"][0]
startlat = d.variables["y_range"][0]
lats=[]
lons = []
coords = [] 
depth = []
fig, ax1 = plt.subplots(1,1)
ret = d.variables["dimension"][:][0]
for i in range(0,d.variables["dimension"][:][1],100):
    for j in range(0,d.variables["dimension"][:][0],100):
        if not np.isnan(z[j+i*ret]):
            lats.append(90-(i*spacing))
            lons.append(j*spacing+startlon)
            coords.append((j*spacing+startlon,90-(i*spacing)))
            depth.append(int(z[j+i*ret]))

mapy = Basemap(projection='ortho', lat_0=90,lon_0=0)
mapy.drawmapboundary(fill_color='aqua')
mapy.fillcontinents(color='coral',lake_color='aqua')
mapy.drawcoastlines()
x,y = mapy(lons,lats)
plt.tricontourf(x,y,depth,levels=20)
plt.colorbar()
fig.suptitle("Bathemetry of the Arctic")
plt.show()

