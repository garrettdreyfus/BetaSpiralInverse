import numpy as np
import operator
import gsw 
from mpl_toolkits.mplot3d import Axes3D
import datetime
import matplotlib.pyplot as plt
from scipy import interpolate
import mygsw
from scipy.interpolate import UnivariateSpline
#import gswmatlab.pyinterface as matgsw


class Profile:
    def __init__(self,eyed, data,ct=False,abssal=False):
        ##id of profiles plus info
        self.eyed = eyed
        self.lat = data["lat"]
        self.f = gsw.f(self.lat)
        self.gamma = (9.8)/(self.f*1025.0)
        self.lon = data["lon"]
        if "time" in data.keys():
            self.time = self.processDate(data["time"])
        if "cruise" in data.keys():
            self.cruise = data["cruise"]#+str(self.time.year)

        if "knownns" in data.keys():
            self.knownns = data["knownns"]
        else:
            self.knownns = {}


        #Temerature Salinity and Pressure
        self.temps = np.asarray(data["temp"])
        self.sals = np.asarray(data["sal"])
        self.pres = np.asarray(data["pres"])
        if not abssal:
            self.sals = gsw.SA_from_SP(self.sals,self.pres,self.lon,self.lat)
        if not ct:
            self.temps = gsw.CT_from_t(self.sals,self.temps,np.abs(self.pres))

        s = np.argsort(self.pres)
        self.temps = self.temps[s]
        self.sals = self.sals[s]
        self.pres = self.pres[s]
        ##Interpolated Temperature, Salinity, and Pressure
        self.itemps = []
        self.isals = []
        self.ipres = []
        theta = np.deg2rad(self.lat)
        r = (90-self.lon) *111*1000
        x = (r*np.cos(theta))
        y = (r*np.sin(theta))
        self.x = x
        self.y = y
        self.idensities = []
        self.neutraldepth = {}
        self.interpolate()
    #reformat date string
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
    ##apply a salinty offset 
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
            self.isals = np.interp(self.ipres,self.pres,self.sals)
            self.itemps = np.interp(self.ipres,self.pres,self.temps)
            self.ialpha = gsw.alpha(self.isals,self.itemps,self.ipres)
            self.ibeta = gsw.beta(self.isals,self.itemps,self.ipres)
            self.idalphadtheta = mygsw.cabbeling_CT_exact(self.isals,self.itemps,self.ipres)
            self.idalphadp = mygsw.thermobaric_CT_exact(self.isals,self.itemps,self.ipres)

            ###using gsw
            self.n2 = gsw.Nsquared(self.isals,self.itemps,self.ipres,self.lat)[0]

    def calculateDensity(self,s,t,p,p_ref=0):
        return gsw.rho(s,t,p)

    #finds PV at a depth
    def potentialVorticityAt(self,depth,debug=False):
        index = np.where(np.asarray(self.ipres) == depth)[0]
        if index:
            index = index[0]
                
            return (self.f/9.8)*np.mean(self.n2[index])

    def potentialVorticityAtHautala(self,depth,debug=False,halfdistance=35):
        index = np.where(np.asarray(self.ipres) == depth)[0]
        if index >36:
            index = index[0]
            densities = gsw.rho(self.isals,\
                    self.itemps,\
                    self.ipres[index])

            drhodz,notvalue = self.dz(depth,densities,35)
            pv = -(self.f/notvalue)*drhodz

            return pv

    def dthetads(self,depth,debug=False,halfdistance=35):
        index = np.where(np.asarray(self.ipres) == depth)[0]
        if index >36:
            index = index[0]
            dthetads,notvalue = self.dz(depth,self.itemps,35,xvalues=self.isals)

            return dthetads

    def dz(self,depth,values,halfdistance,xvalues=[]):
        if ~np.any(xvalues): xvalues=self.ipres

        index = np.where(np.asarray(self.ipres) == depth)[0]
        if index >halfdistance+1:
            index = index[0]
            values = values[index-halfdistance-1:index+halfdistance+1]
            xvalues = xvalues[index-halfdistance-1:index+halfdistance+1]

            L = sorted(zip(xvalues,values), key=operator.itemgetter(0))
            xvalues, values = zip(*L)

            if len(values)<5:
                print("what de")
            try:
                spl = UnivariateSpline(xvalues,values)
            except:
                print("THE SPLINE IS COOKED")
                print([xvalues,values])
                return [np.nan,np.nan]

            dvaluedz = (spl(xvalues[0])-spl(xvalues[-1]))/(xvalues[-1]-xvalues[0])
                
            return dvaluedz,spl(xvalues[halfdistance])
        return [np.nan,np.nan]

    def vertGrad(self,depth,quant,halfdistance=35):
        if quant == "s":
            values = self.isals
        if quant == "t":
            values = self.itemps
        if quant == "alpha":
            values = self.ialpha
        if quant == "beta":
            values = self.ibeta
        return self.dz(depth,values,halfdistance)[0]


    #finds average PV between two depths
    def potentialVorticityBetween(self,above,at,below,debug=False,method="interp"):
        above = np.where(np.asarray(self.ipres) == int(above))[0]
        below = np.where(np.asarray(self.ipres) == int(below))[0]
        at = np.where(np.asarray(self.ipres) == int(at))[0][0]
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
        #if method == "mean":
            #pv = (self.f/9.8)*np.mean(self.n2[min(above,below):max(above,below)])
        elif method == "interp":
            spl = UnivariateSpline(self.ipres[min(above,below):max(above,below)],self.n2[min(above,below):max(above,below)])
            pv = (self.f/9.8)*spl(self.ipres[at])
        if pv<0 and (debug ):
            p = np.where(self.n2>0)
            m = np.where(self.n2<0)
            spl = UnivariateSpline(self.ipres[min(above,below):max(above,below)],self.n2[min(above,below):max(above,below)])
            plt.plot(spl(self.ipres[min(above,below):max(above,below)]),self.ipres[min(above,below):max(above,below)])
            plt.gca().axhline(self.ipres[at])
            plt.show()
            print("\n########")
            print("\nless than 0",np.mean(np.abs(self.n2[m])))
            print("\nmore than 0",np.mean(np.abs(self.n2[p])))
            print("\n########")
        return pv
    
    ##returns the t s at a pressure
    def atPres(self,pres,full=False):
        i = np.where(np.asarray(self.ipres) == int(pres))[0][0]
        return self.itemps[i], self.isals[i], self.ialpha[i],self.ibeta[i],self.idalphadtheta[i],self.idalphadp[i]
    
    ##returns the index at pressure
    def presIndex(self,pres):
        i = np.argmin(np.abs(np.asarray(self.pres) - pres))#[0][0]
        return i

    ##returns the index at pressure
    def ipresIndex(self,pres):
        i = np.argmin(np.abs(np.asarray(self.ipres) - pres))#[0][0]
        return i

    ##returns t and s between two pressures
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

    #I think this makes sense if you read Trevor and Mcdougal
    #essentialy find depth between on one profile at which the 
    # at which potential density reference to the average pressure between
    ## the two points is equal
    def neutralDepthAlt(self,p2,depth,debug=False,searchrange=100,depthname=None):
        searchrange = int(min(abs(self.ipres[0]-depth),abs(self.ipres[-1]-depth),\
                            abs(p2.ipres[0]-depth),abs(p2.ipres[-1]-depth),searchrange))
        try:
            if searchrange <2:
                return None
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
        if depth>100:
            plt.plot(p2.ipres[startindexp2:startindexp2+len(Es)],Es)
            for j in zero_crossings:
                plt.gca().axvline(p2.ipres[startindexp2]+j)
            plt.show()
        smallest = np.argmin(np.abs(Es))
        #plt.plot(Es,self.ipres[startindexself:startindexself+2*searchrange])
        #plt.show()
        if len(zero_crossings)>=1 :
            if abs(p2.ipres[zero_crossings[0]] - p2.ipres[zero_crossings[-1]])>100:
                return None
            a  =np.asarray(startindexp2+zero_crossings)
            p2.neutraldepth[depthname] = np.mean(np.asarray(p2.ipres)[a])
            #print("More than one crossing")
            return p2.neutraldepth[depthname]
        elif abs(Es[smallest])<0.005:
            #print("One zero")
            #p2.neutraldepth[depthname] = p2.ipres[startindexp2+smallest]
            #return p2.neutraldepth[depthname]
            return None
        else:

            #print("not a single zero crossing")
            #print(self.cruise,p2.cruise)
            #p2.neutraldepth[depthname] = p2.ipres[startindexp2]
            #return p2.ipres[startindexp2]
            return None

    def neutralDepth(self,p2,depth,debug=False,searchrange=100,depthname=None):
        depth = int(depth)
        plowerbound = min(self.ipres[0],p2.ipres[0])
        pupperbound = min(self.ipres[-1],p2.ipres[-1])
        at = np.where(np.asarray(self.ipres) == depth)[0][0]
        prange = pupperbound - plowerbound
        p2offset = p2.ipres[0] - plowerbound
        selfoffset = self.ipres[0] - plowerbound
        if self.ipres[0] < p2.ipres[0]:
            p2offset = 0
            selfoffset = np.where(np.asarray(self.ipres) == p2.ipres[0])[0][0]
        elif p2.ipres[0] < self.ipres[0]:
            selfoffset = 0
            p2offset = np.where(np.asarray(p2.ipres) == self.ipres[0])[0][0]
        else:
            selfoffset = 0
            p2offset = 0

        if p2.ipres[p2offset] != self.ipres[selfoffset]:
            print(p2.ipres[p2offset],self.ipres[selfoffset])
    
        depths =  (np.asarray(p2.ipres[p2offset:p2offset+prange]) +self.ipres[at])/2.0
        
        p2densities = gsw.rho(p2.isals[p2offset:p2offset+len(depths)],\
                p2.itemps[p2offset:p2offset+len(depths)],\
                depths)


        selfdensities = gsw.rho([self.isals[at]]*(len(depths)),\
                [self.itemps[at]]*(len(depths)),\
                depths)

        minlen = min(len(p2densities),len(selfdensities))
        Es = p2densities[:minlen]-selfdensities[:minlen]
        if len(Es)<2:
            return None

        zero_crossings = np.where(np.diff(np.sign(Es)))[0]
        smallest = np.argmin(np.abs(Es))
        if len(zero_crossings)>=1 :
            if abs(p2.ipres[zero_crossings[0]] - p2.ipres[zero_crossings[-1]])>100:
                return None

            a  =np.asarray(p2offset+zero_crossings)
            p2.neutraldepth[depthname] = np.mean(np.asarray(p2.ipres)[a])
            #print("More than one crossing")
            return p2.neutraldepth[depthname]
        else:
            return None


    ## this thing is a little bit wild
    ## so in the matlab version of gsw they have a function that 
    ## gives the geostrophic stream function on neutral surfaces
    ## this does not exist in the python/c version of the gsw.
    ## so I created a port of it to python, annoying part is
    ## it requires some other functions also not in the python/c gsw
    ## but luckily they are in the python only gsw
    ## so I stole those functions and made a frankenstein that is like
    ## 100000 times faster than the matlab geostrophic function thing ;)
    def geoIsopycnal(self,ns,nsdensref):
        return mygsw.geo_strf_isopycnal(self.isals,self.itemps,self.ipres,1800,nsdensref,ns)
