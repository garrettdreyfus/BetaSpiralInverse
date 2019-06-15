import numpy as np
import gsw 
from mpl_toolkits.mplot3d import Axes3D

class Profile:
    def __init__(self, eyed, data):
        ##id of profiles plus info
        self.eyed = eyed
        self.lat = data["lat"]
        self.lon = data["lon"]
        self.time = data["time"]
        #Temerature Salinity and Pressure
        self.temps = np.asarray(data["temp"])
        self.sals = np.asarray(data["sal"])
        self.pres = np.asarray(data["pres"])
        ##Interpolated Temperature, Salinity, and Pressure
        self.itemps = []
        self.isals = []
        self.ipres = []
        self.idensities = []
        self.neutraldepth = {}
        self.interpolate()

    #interpolate all quantities on a 1 dbar line
    def interpolate(self):
        self.ipres = range(int(min(self.pres)),int(max(self.pres)))
        self.itemps = np.interp(self.ipres,self.pres,self.temps)
        self.isals = np.interp(self.ipres,self.pres,self.sals)
        self.idensities = gsw.pot_rho_t_exact(self.isals, self.itemps, self.ipres, [0]*len(self.ipres))
            
    #
    def neutralDepthWrong(self,p2,depth,debug=False,searchrange=50):
        try:
            startindexself = depth-self.ipres[0]-searchrange
            startindexp2 = depth-p2.ipres[0]-searchrange
        except:
            print("This is off:   ", p2.eyed)
            return None
        if startindexp2 < 0 or startindexp2 + searchrange*2 > len(p2.ipres):
            return None
        Es = self.idensities[startindexself:startindexself + 2*searchrange]-p2.idensities[startindexp2:startindexp2 + 2*searchrange] 
        E = np.argmin(Es)
        self.neutraldepth[depth] = self.ipres[startindexself+E]
        if abs(Es[E])<0.01:
            return self.ipres[startindexself+E]
        else:
            return None
    def neutralDepth(self,p2,depth,debug=False,searchrange=50):
        try:
            startindexself = depth-self.ipres[0]-searchrange
            startindexp2 = depth-p2.ipres[0]-searchrange
        except:
            print("This is off:   ", p2.eyed)
            return None
        if startindexp2 < 0 or startindexp2 + searchrange*2 > len(p2.ipres):
            return None
        Es = p2.idensities[startindexp2:startindexp2 + 2*searchrange]-self.idensities[startindexself+searchrange]
        E = np.argmin(Es)
        p2.neutraldepth[depth] = p2.ipres[startindexp2+E]
        if abs(Es[E])<0.01:
            return p2.ipres[startindexp2+E]
        else:
            return None


