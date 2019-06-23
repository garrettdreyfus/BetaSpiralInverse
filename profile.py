import numpy as np
import gsw 
from mpl_toolkits.mplot3d import Axes3D
import datetime

class Profile:
    def __init__(self, eyed, data):
        ##id of profiles plus info
        self.eyed = eyed
        self.lat = data["lat"]
        self.lon = data["lon"]
        self.time = self.processDate(data["time"])
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

    def processDate(self,datestring):
        split = datestring.split("-")
        m = int(split[1])
        y = int(split[0])
        d = int(split[2].split("T")[0])
        h = int(split[2].split("T")[1].split(":")[0])
        return datetime.datetime(y,m,d,h)
        

    #interpolate all quantities on a 1 dbar line
    def interpolate(self):
        self.ipres = range(int(min(self.pres)),int(max(self.pres)))
        self.isals = np.interp(self.ipres,self.pres,self.sals)
        self.itemps = gsw.pt_from_t(self.isals,np.interp(self.ipres,self.pres,self.temps),self.ipres,0)
            
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
    def neutralDepthWronger(self,p2,depth,debug=False,searchrange=50):
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
    def calculateDensity(self,s,t,p,p_ref=0):
        return gsw.rho(s,t,p)

    def neutralDepth(self,p2,depth,debug=False,searchrange=50):
        try:
            startindexself = depth-self.ipres[0]-searchrange
            if abs(self.ipres[startindexself] - depth) != 50:
                print("What the")
            startindexp2 = depth-p2.ipres[0]-searchrange
            if abs(p2.ipres[startindexp2] - depth) != 50:
                print(p2.ipres[startindexp2],depth)
                print("What the")
        except:
            return None
        if startindexp2 < 0 or startindexp2 + searchrange*2 > len(p2.ipres):
            return None
        Es = []
        for index in range(2*searchrange):
            p2density=self.calculateDensity(p2.isals[index+startindexp2],p2.itemps[index+startindexp2],(p2.ipres[index+startindexp2]+depth)/2.0)
            selfdensity=self.calculateDensity(self.isals[startindexself+searchrange],self.itemps[startindexself+searchrange],(p2.ipres[index+startindexp2]+depth)/2.0)
            Es.append(p2density-selfdensity)
        zero_crossings = np.where(np.diff(np.sign(Es)))[0]
        zeroes = np.where(Es ==0)[0]
        if len(zero_crossings)==1 :
            p2.neutraldepth[depth] = p2.ipres[startindexp2+zero_crossings[0]]
            return p2.ipres[startindexp2+zero_crossings[0]]
        elif len(zeroes)==1:
            p2.neutraldepth[depth] = p2.ipres[startindexp2+zeroes[0]]
            return p2.ipres[startindexp2+zeroes[0]]

        else:
            #print("not 1 zero crossing",len(zero_crossings))
            return None
