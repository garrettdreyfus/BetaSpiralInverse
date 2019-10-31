import scipy
import scipy.io as sio

linedata = sio.loadmat("data/aleutianline.mat")
lat = linedata["lat_BS"].T
lon = linedata["lon_BS"]
