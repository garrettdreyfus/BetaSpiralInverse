import numpy as np
import bathtools
import gsw
from progress.bar import Bar
import pickle
import inverttools as inv
import matplotlib.pyplot as plt

def bathVarTerm(lat,lon):
    d=bathtools.bathBox(lat,lon)
    dnot = 750
    return (np.var(d)/dnot)**(0.25)

def saveBathVarTermCache(surfaces,outfilename):
    bathVar = {}
    coords = []
    for k in surfaces.keys():
        coords += zip(surfaces[k]["lons"],surfaces[k]["lats"])
    coords = set(coords)
    for p in  Bar("var coord :").iter(coords):
        if not (np.isnan(p[0]) and np.isnan(p[1])):
            d=bathtools.bathBox(p[1],p[0])
            dnot = 750
            bathVar[p] = (np.var(d)/dnot)**(0.25)
    with open(outfilename, 'wb') as outfile:
        pickle.dump(bathVar, outfile)



def bathVarTermCache(lat,lon,filename):
    if not hasattr(bathVarTermCache,"p"):
        bathVarTermCache.p = pickle.load(open(filename, 'rb'))

        print("loaded bath variance file")
        
    return bathVarTermCache.p[(lon,lat)]



def Kv(lat,lon,pv,pres,cachename=None):
    if cachename:
       bVT = bathVarTermCache(lat,lon,cachename) 
    else:
        bVT = bathVarTerm(lat,lon)
    return bVT*np.exp(-(abs(bathtools.searchBath(lat,lon))-abs(pres))/200)

def kChecker(surfaces,k,found,debug=False):
    f = gsw.f(surfaces[k]["lats"][found])
    x = surfaces[k]["x"][found]
    y = surfaces[k]["y"][found]
    r = np.sqrt(x**2 + y**2)
    hx = surfaces[k]["data"]["hx"][found]
    hy = surfaces[k]["data"]["hy"][found]
    dsdx = surfaces[k]["data"]["dsdx"][found]
    dsdy = surfaces[k]["data"]["dsdy"][found]
    pres = surfaces[k]["data"]["pres"][found]
    alpha = surfaces[k]["data"]["alpha"][found] 
    betaTherm = surfaces[k]["data"]["beta"][found] 
    dsdz =  surfaces[k]["data"]["dsdz"][found] 
    d2sdx2 =  surfaces[k]["data"]["d2sdx2"][found] 
    d2sdy2 =  surfaces[k]["data"]["d2sdy2"][found] 
    dalphadtheta = surfaces[k]["data"]["dalphadtheta"][found] 
    dalphads = surfaces[k]["data"]["dalphads"][found] 
    dalphadp = surfaces[k]["data"]["dalphadp"][found] 
    dbetadp = surfaces[k]["data"]["dbetadp"][found] 
    dbetads = surfaces[k]["data"]["dbetads"][found] 
    dtdx = surfaces[k]["data"]["dtdx"][found] 
    dtdy = surfaces[k]["data"]["dtdy"][found] 
    dqnotdx = surfaces[k]["data"]["dqnotdx"][found] 
    dqnotdy = surfaces[k]["data"]["dqnotdy"][found] 
    dqdz = surfaces[k]["data"]["dqdz"][found] 
    d2qdz2 = surfaces[k]["data"]["d2qdz2"][found] 
    dpdx = surfaces[k]["data"]["dpdx"][found] 
    dpdy = surfaces[k]["data"]["dpdy"][found] 
    dqdx = surfaces[k]["data"]["dqdx"][found] 
    dqdy = surfaces[k]["data"]["dqdy"][found] 
    d2qdx2 = surfaces[k]["data"]["d2qdx2"][found] 
    d2qdy2 = surfaces[k]["data"]["d2qdy2"][found] 
    khpdz = surfaces[k]["data"]["khpdz"][found] 
    alphat = dalphadtheta+2*(alpha/betaTherm)*dalphads-(alpha**2/betaTherm**2)*dbetads
    alphap = dalphadp -(alpha/betaTherm)*dbetadp
    pv =  surfaces[k]["data"]["pv"][found] 
    doublets =  surfaces[k]["data"]["d2thetads2"][found] 
    CKVB =  surfaces[k]["data"]["CKVB"][found] 
    f = gsw.f(surfaces[k]["lats"][found])
    beta = inv.calcBeta(surfaces[k]["lats"][found])
    isitnan = [alpha,betaTherm,dsdz,hx,hy,dsdx,dsdy,pres,d2sdx2,d2sdy2,\
              dalphadtheta,dalphads,dalphadp,dbetadp,dbetads,dtdx,dtdy,\
              dqnotdx,dqnotdy,dpdx,dpdy,alphat,alphap,pv,doublets,CKVB,\
              beta,d2qdx2,d2qdy2,khpdz]
    kvbscale = 1.0/(5*(10**-5)) #(1/(10**-1.09))
#/630957)#5*(10**2)#5*(10**-5)
    kvoscale = 1.0/(5*(10**-6))#(1/0.001)#5*(10**1)#5*(10**-6)
    khscale = 1.0/500#1/(10**-2.4))
# 10**12#500

    if (np.isnan(isitnan).any()):
        if debug:
            print("pres is nan: ",np.isnan(pres))
            print("hx is nan: ",np.isnan(hx))
            print("hy is nan: ",np.isnan(hy))
            print("x is nan: ",np.isnan(x))
            print("y is nan: ",np.isnan(y))
            print("something here is nan")
        return None, None
    if not (np.isnan(isitnan).any()):
        eyed = surfaces[k]["ids"][found]
        above = np.where(np.asarray(surfaces[k-200]["ids"]) == eyed)[0]
        aboveabove = np.where(np.asarray(surfaces[k-400]["ids"]) == eyed)[0]
        belowbelow = np.where(np.asarray(surfaces[k+400]["ids"]) == eyed)[0]
        below= np.where(np.asarray(surfaces[k+200]["ids"]) == eyed)[0]
        depths = [k-400,k-200,k,k+200,k+400]
        matchingindex = [aboveabove,above,found,below,belowbelow]
        qs = []
        s = []
        t = []
        for j in range(len(depths)):
            qs.append(surfaces[depths[j]]["data"]["pv"][matchingindex[j]])
            s.append(surfaces[depths[j]]["data"]["s"][matchingindex[j]])
            t.append(surfaces[depths[j]]["data"]["t"][matchingindex[j]])
        kv0inspect = False
        if kv0inspect:
            plt.plot(qs,depths)
            plt.gca().set_title(str(str(k) + ": dqdz: "+str(dqdz)+"\ndq2dz2: "+str(d2qdz2)))
            plt.gca().invert_yaxis()
            plt.show()

            plt.plot(s,depths)
            plt.gca().set_title(str(str(k) + ": dsdz: "+str(dsdz)))
            plt.gca().invert_yaxis()
            plt.show()
            
            plt.plot(s,t)
            plt.gca().set_title(str(str(k) + ": doublets: "+str(doublets)))
            plt.show()

        pvkvb = (d2qdz2+2*(1/1000)*dqdz+(1/(1000**2))*pv)*CKVB
        pvkv0 = d2qdz2
        pvkh = (d2qdx2+d2qdy2)-2*(dqnotdx*dqdx+dqnotdy*dqdy)/pv -f*khpdz
        skvo = -alpha*f*(1/pv)*(dsdz**3)*doublets
        skvb = skvo*CKVB
        skhpart1 = (f/pv)*dsdz*(alphat*(dtdx**2 + dtdy**2)+alphap*(dtdx*dpdx+dtdy*dpdy))
        skhpart2 = (d2sdx2+d2sdy2)-2*(dqnotdx*dsdx + dqnotdy*dsdy)/pv
        skh = skhpart1 + skhpart2
        
        kv0breakdown=True
        if kv0breakdown:
            labels = ["d2qdz2","alpha","f","1/pv","dsdz**3","doublets"]
            values = [d2qdz2,alpha,f,1/pv,dsdz**3,doublets]
            plt.bar(labels,values)
            plt.yscale("log")
            plt.show()

        sixpartcompare = True
        if sixpartcompare:
            labels = ["pvkv0","pvkh","pvkvb","skv0","skh","skvb"]
            values = [pvkv0,pvkh,pvkvb,skvo,skh,skvb]
            plt.bar(labels,np.abs(values))
            plt.yscale("log")
            plt.show()
            
        return np.asarray([-pvkv0/kvoscale,-pvkvb/kvbscale,-pvkh/khscale]),np.asarray([-skvo/kvoscale,-skvb/kvbscale,-skh/khscale])


