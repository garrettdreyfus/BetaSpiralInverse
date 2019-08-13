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
        self.gamma = (9.8)/(self.f*1025.0)
        self.lon = data["lon"]
        self.time = self.processDate(data["time"])
        self.cruise = data["cruise"]#+str(self.time.year)
        #Temerature Salinity and Pressure
        self.temps = np.asarray(data["temp"])
        self.sals = np.asarray(data["sal"])
        self.pres = np.asarray(data["pres"])
        s = np.argsort(self.pres)
        self.temps = self.temps[s]
        self.sals = self.sals[s]
        self.temps = gsw.CT_from_t(self.sals,self.temps,np.abs(self.pres))
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
        if len(self.pres)>4:
            tck = interpolate.splrep(self.pres,self.sals)
            self.isals = interpolate.splev(self.ipres,tck)
            self.itemps = np.interp(self.ipres,self.pres,self.temps)
            #self.irhos = gsw.rho(self.isals,self.itemps,self.ipres)
            #self.n2 = (9.8/1025.0)*np.gradient(self.irhos,-np.asarray(self.ipres))

            ###using gsw
            self.n2 = gsw.Nsquared(self.isals,self.itemps,self.ipres,self.lat)[0]
            #tck = interpolate.splrep(self.ipres[:-1],self.n2,s=0.01)
            #self.n2 = interpolate.splev(self.ipres,tck) 

    def calculateDensity(self,s,t,p,p_ref=0):
        return gsw.rho(s,t,p)


    def potentialVorticityAt(self,index,debug=False):
        index = np.where(np.asarray(self.ipres) == index)[0]
        if index:
            index = index[0]
            #if self.n2[index] < 0 and debug :
                #print(self.isals[index],self.itemps[index],self.ipres[index])
                
            return (self.f/9.8)*np.mean(self.n2[index])


    def potentialVorticityBetween(self,above,below,debug=False):
        above = np.where(np.asarray(self.ipres) == int(above))[0]
        below = np.where(np.asarray(self.ipres) == int(below))[0]
        if len(above)>0 and len(below) > 0:
            above = above[0]
            below = below[0]
        elif len(above)>0:
            above = above[0]
            below = min(self.ipres[-1],above+200)
            below = np.where(np.asarray(self.ipres) == int(below))[0][0]
        elif len(below)>0:
            above = max(self.ipres[0],above-200)
            above = np.where(np.asarray(self.ipres) == int(above))[0][0]
            below = below[0]
        else:
            return
        if np.isnan(self.n2[min(above,below):max(above,below)]).any():
            print("###################")
            print("######ipres#######")
            print(self.ipres[min(above,below):max(above,below)])
            print("###################")
            print("######itemps#######")
            print(self.itemps[min(above,below):max(above,below)])
            print("###################")
            print("######isals#######")
            print(self.isals[min(above,below):max(above,below)])
            print("###################")
            print("######n2#######")
            print(self.n2[min(above,below):max(above,below)])
        pv = (self.f/9.8)*np.mean(self.n2[min(above,below):max(above,below)])
        if pv<0 and False:
            p = np.where(self.n2>0)
            m = np.where(self.n2<0)
            print("\n########")
            print("\nless than 0",np.mean(np.abs(self.n2[m])))
            print("\nmore than 0",np.mean(np.abs(self.n2[p])))
            print("\n########")
        return pv

    def atPres(self,pres):
        i = np.where(np.asarray(self.ipres) == int(pres))[0][0]
        return self.itemps[i], self.isals[i]

    def betweenPres(self,above,below):
        above = np.where(np.asarray(self.ipres) == int(above))[0]
        below = np.where(np.asarray(self.ipres) == int(below))[0]
        if len(above)>0 and len(below) > 0:
            above = above[0]
            below = below[0]
        elif len(above)>0:
            above = above[0]
            below = len(self.ipres)-1
        elif len(below)>0:
            below = below[0]
            above = 0
        else:
            return
        return np.mean(self.itemps[min(above,below):max(above,below)]), np.mean(self.isals[min(above,below):max(above,below)])


    def densityAtPres(self,pres,ref=0):
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
        ns = ns[::-1]
        nsdensref = nsdensref[::-1]
        self.ipres = np.abs(np.asarray(self.ipres))

        dyn_height = gsw.geo_strf_dyn_height(self.isals,self.itemps,self.ipres,10.25)
        
        #print("###########")
        #print(ns)
        ns, idxs = np.unique(ns,return_index=True)
        nsdensref = nsdensref[idxs]
        #print("between")
        #print(ns)
        #print("###########")
        filt = np.where(np.isin(self.ipres,ns))
        #print("################")
        #print(ns)
        #print(self.ipres[filt])
        #print("################")
        nsdyn_height = dyn_height[filt]
        nstemps = self.itemps[filt]
        nssals = self.isals[filt]
        #print(len(nstemps))
        ###
        db2Pa = 1e4
        sa_iref_cast,ct_iref_cast,p_iref_cast = mygsw.interp_ref_cast(nsdensref,"s2")
        cp0 = 3991.86795711963
        #print("##################")
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
        #print("##################")
        return part1 + part2 + part3 + part4 + part5
