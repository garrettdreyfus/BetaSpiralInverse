import scipy
import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt

linedata = sio.loadmat("data/aleutianline.mat")
lat = linedata["lat_BS"].T
lon = linedata["lon_BS"]
lat = list(lat[0])
lat.insert(0,56.084)
lat.insert(0,59.13)
lat = np.asarray(lat)
lon = list(lon[0])
lon.insert(0,-160.02)
lon.insert(0,-158.89)
lon = np.asarray(lon)


