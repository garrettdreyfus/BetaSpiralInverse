import numpy as np
import gsw 
from mpl_toolkits.mplot3d import Axes3D
import datetime
import matplotlib.pyplot as plt
from scipy import interpolate
import mygsw


class Profile:
    def __init__(self,eyed, data):
        ##id of profiles plus info
        self.eyed = eyed
        self.lat = data["lat"]
        self.f = gsw.f(self.lat)
        self.gamma = (9.8)/(f*1025.0)
        self.lon = data["lon"]
        self.time = self.processDate(data["time"])
        self.cruise = data["cruise"]+str(self.time.year)
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
        try:
            d = int(split[2].split("T")[0])
            #h = int(split[2].split("T")[1].split(":")[0])
            result = datetime.datetime(y,m,d)
        except:
            print(datestring,m,y)
            d = np.clip(int(split[2].split("T")[0]),1,27)
            result = datetime.datetime(y,m,d)
            print(datestring,result)
        return result
    
    def applyOffset(self,offset):
        self.sals =self.sals +offset
        self.ipres=[]
        self.isals=[]
        self.itemps=[]
        self.interpolate()
        

    #interpolate all quantities on a 1 dbar line
    def interpolate(self):
        self.ipres = range(int(min(self.pres)),int(max(self.pres)))
        tck = interpolate.splrep(self.pres,self.sals,s=20)
        self.isals = interpolate.splev(self.ipres,tck) 
        #self.isals = np.interp(self.ipres,self.pres,self.sals)
        self.temps = gsw.CT_from_t(self.sals,self.temps,self.pres)
        self.itemps = gsw.CT_from_t(self.isals,np.interp(self.ipres,self.pres,self.temps),self.ipres)
            
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

    def potentialVorticity(self,p,debug=False):
        index = np.where(np.asarray(self.ipres) == p)[0]
        if index!= None:
            index = index[0]
            n2 = np.mean(gsw.Nsquared(self.isals[index-5:index+5],self.itemps[index-5:index+5],
                    self.ipres[index-5:index+5], self.lat)[0])
            if n2 < 0 and debug:
                #plt.scatter(gsw.Nsquared(self.isals,self.itemps,self.ipres)[0],self.ipres[:-1],s=0.5)
                plt.plot(self.isals,self.ipres)
                plt.gca().invert_yaxis()
                plt.show()
            return (self.f/9.8)*(n2)
        else:
            return None

    def atPres(self,pres):
        i = np.where(np.asarray(self.ipres) == int(pres))[0][0]
        return self.itemps[i], self.isals[i]

    def densityAtPres(self,pres,ref=0):
        print(pres)
        i = np.where(np.asarray(self.ipres) == int(pres))[0][0]
        return gsw.rho(self.isals[i],self.itemps[i],ref)

    def sigma2(self,pres):
        i = np.where(np.asarray(self.ipres) == int(pres))[0][0]
        return gsw.sigma2(self.isals[i],self.itemps[i])
    
    def rhoZ(self,pres):
        return (densityAtPres(pres+1)-densityAtPres(pres-1))/3

    def neutralDepth(self,p2,depth,debug=False,searchrange=50,depthname=None):
        try:
            depth=int(depth)
            if depthname ==None:
                depthname=depth
            startindexself = depth-self.ipres[0]-searchrange
            if startindexself not in range(len(self.ipres)) or startindexself + 2*searchrange not in range(len(self.ipres)):
                return None
            if abs(self.ipres[startindexself] - depth) != searchrange:
                print("What the")
            startindexp2 = depth-p2.ipres[0]-searchrange
            if startindexp2 not in range(len(p2.ipres)) or startindexp2 + 2*searchrange not in range(len(self.ipres)):
                return None
            if abs(p2.ipres[startindexp2] - depth) != searchrange:
                print(p2.ipres[startindexp2],depth)
                print("What the")
        except Exception as e:
            print(str(e))
            print("strange exception")
            return None
        if startindexp2 < 0 or startindexp2 + searchrange*2 > len(p2.ipres):
            return None
        Es = []
        #print(depth,p2.ipres[startindexp2+searchrange],self.ipres[startindexself+searchrange])
        for index in range(2*searchrange):
            p2density=self.calculateDensity(p2.isals[index+startindexp2],p2.itemps[index+startindexp2],(p2.ipres[index+startindexp2]+depth)/2.0)
            selfdensity=self.calculateDensity(self.isals[startindexself+searchrange],self.itemps[startindexself+searchrange],(p2.ipres[index+startindexp2]+depth)/2.0)
            Es.append(p2density-selfdensity)
        zero_crossings = np.where(np.diff(np.sign(Es)))[0]
        zeroes = np.where(np.abs(Es) <0.1)[0]
        #plt.plot(Es,self.ipres[startindexself:startindexself+2*searchrange])
        #plt.show()
        if len(zero_crossings)==1 :
            p2.neutraldepth[depthname] = p2.ipres[startindexp2+zero_crossings[0]]
            #print("One zero crossing")
            return p2.ipres[startindexp2+zero_crossings[0]]
        elif len(zero_crossings) > 1 :
            p2.neutraldepth[depthname] = (p2.ipres[startindexp2+zero_crossings[0]] + p2.ipres[startindexp2+zero_crossings[-1]])/2.0
            #print("More than one crossing")
            return p2.neutraldepth[depthname]
        elif len(zeroes)==1:
            #print("One zero")
            p2.neutraldepth[depthname] = p2.ipres[startindexp2+zeroes[0]]
            return p2.ipres[startindexp2+zeroes[0]]

        else:
            #print("not a single zero crossing")
            #print(self.cruise,p2.cruise)
            return None
    def geoIsopycnal(self,ns,nsdensref):
        self.ipres = np.abs(np.asarray(self.ipres))
        dyn_height = gsw.geo_strf_dyn_height(self.isals,self.itemps,self.ipres,10.1235)
        nsuniques= np.unique(np.asarray(np.abs(ns)))
        sames = np.setdiff1d(range(len(ns)),nsuniques)
        ns[sames[1:]] = ns[sames[1:]]+1
        filt = np.where(np.isin(self.ipres,ns))
        #print("################")
        #print(ns)
        #print(self.ipres[filt])
        #print("################")
        nsdyn_height = dyn_height[filt]
        nstemps = self.itemps[filt]
        nssals = self.isals[filt]
        print(len(nstemps))
        ##robbing from gibbs
        #Things missing
        #   Iref cast variables
        #   enthalpy_SSO_0
        ###
        db2Pa = 1e4
        sa_iref_cast,ct_iref_cast,p_iref_cast = mygsw.interp_ref_cast(nsdensref,"s2")
        cp0 = 3991.86795711963
        print("##################")
        #print("py iref_cast: ",p_iref_cast,ct_iref_cast,sa_iref_cast)
        #print("py nssal and nstemps: ",nssals,nstemps,ns)
        part1 = 0.5 *db2Pa*(ns-p_iref_cast)*(gsw.specvol(nssals,nstemps,ns)-gsw.specvol(sa_iref_cast,ct_iref_cast,ns))
        part2=0
        part3 = (-0.225e-15)*(db2Pa*db2Pa)*(nstemps-ct_iref_cast)*(ns-p_iref_cast)*(ns-p_iref_cast)
        part4 = nsdyn_height - mygsw.enthalpy_SSO_0(ns)
        part5 = gsw.enthalpy(sa_iref_cast,ct_iref_cast,ns) -cp0*ct_iref_cast
        #print("specvol delta", (gsw.specvol(nssals,nstemps,ns)-gsw.specvol(sa_iref_cast,ct_iref_cast,ns)))
        #print("py part1: ",part1+part2)
        #print("py part2: ",part3)
        #print("dyn height", nsdyn_height)
        #print("SS0 enthalpy",  mygsw.enthalpy_SSO_0(ns))
        #print("py part3: ",part4+part5)
        #print("result:", part1+part2+part3+part4+part5)
        print("##################")
        return part1 + part2 + part3 + part4 + part5
