import numpy as np
import bathtools
import gsw
from progress.bar import Bar
import pickle
import inverttools as inv
import matplotlib.pyplot as plt

def calcBeta(lat):
    omega =  (7.2921 * 10**-5)
    a = 6.357 * (10**6)
    return (2*omega*np.cos(lat))/a

## calcute the bathymetric variability term
def bathVarTerm(lat,lon):
    d=bathtools.bathBox(lat,lon)
    dnot = 750
    return (np.var(d)/dnot)**(0.25)

#one of the mixing terms is a mess and this calculates that
def calculateKHP(staggered,k,index):
    dalphadtheta = staggered[k]["data"]["dalphadtheta"][index]
    dalphadp = staggered[k]["data"]["dalphadp"][index]
    dalphads = staggered[k]["data"]["dalphads"][index]
    dbetads = staggered[k]["data"]["dbetads"][index]
    alpha = staggered[k]["data"]["alpha"][index]
    dbetadp = staggered[k]["data"]["dbetadp"][index]
    betaTherm = staggered[k]["data"]["beta"][index]
    alphat = dalphadtheta+2*(alpha/betaTherm)*dalphads-(alpha**2/betaTherm**2)*dbetads
    alphap = dalphadp -(alpha/betaTherm)*dbetadp
    magct = staggered[k]["data"]["dtdx"][index]**2 + staggered[k]["data"]["dtdy"][index]**2
    cdotp = staggered[k]["data"]["dtdx"][index]*staggered[k]["data"]["dpdx"][index]+staggered[k]["data"]["dtdy"][index]*staggered[k]["data"]["dpdy"][index]
    return alphat*magct+alphap*cdotp


## this stuff takes a long time to calculate so this
## throws everything into a pickle for later
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


## extracts bathVar data from pickle but only does so once
def bathVarTermCache(lat,lon,filename):
    if not hasattr(bathVarTermCache,"p"):
        bathVarTermCache.p = pickle.load(open(filename, 'rb'))

        print("loaded bath variance file")
        
    return bathVarTermCache.p[(lon,lat)]


## calulate the Kvb bathvar coefficient
def Kv(lat,lon,pv,pres,cachename=None):
    if cachename:
       bVT = bathVarTermCache(lat,lon,cachename) 
    else:
        bVT = bathVarTerm(lat,lon)
    return bVT*np.exp(-(abs(bathtools.searchBath(lat,lon))-abs(pres))/500)

#function for exploring k mixing term values
def kChecker(surfaces,k,found,scales,debug=False):
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
    beta = calcBeta(surfaces[k]["lats"][found])
    isitnan = [alpha,betaTherm,dsdz,hx,hy,dsdx,dsdy,pres,d2sdx2,d2sdy2,\
              dalphadtheta,dalphads,dalphadp,dbetadp,dbetads,dtdx,dtdy,\
              dqnotdx,dqnotdy,dpdx,dpdy,alphat,alphap,pv,doublets,CKVB,\
              beta,d2qdx2,d2qdy2,khpdz]
    kvoscale = scales["kvo"]
    kvbscale = scales["kvb"]
    khscale  = scales["kh"]

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
        
        kv0breakdown=False
        if kv0breakdown:
            labels = ["d2qdz2","alpha","f","1/pv","dsdz**3","doublets","skvo"]
            print(dsdz)
            values = [np.abs(d2qdz2),alpha,f,1/pv,np.abs(dsdz)**3,np.abs(doublets),np.abs(skvo)]
            plt.bar(labels,values)
            plt.yscale("log")
            plt.show()

        khbreakdown=False
        if khbreakdown:
            #(d2qdx2+d2qdy2)-2*(dqnotdx*dqdx+dqnotdy*dqdy)/pv -f*khpdz
            labels = ["d2qdx2+d2qdy2","2*(dqnotdx*dqdx+dqnotdy*dqdy)/pv","f*khpdz","f","khpdz"]
            print(dsdz)
            values = [np.abs(d2qdx2+d2qdy2),np.abs(2*(dqnotdx*dqdx+dqnotdy*dqdy)/pv),np.abs(f*khpdz),f,abs(khpdz)]
            plt.bar(labels,values)
            plt.yscale("log")
            plt.show()

        kvbbreakdown=True
        if kvbbreakdown:
            labels = ["d2qdz2","dqdz","pv"]
            print(dsdz)
            values = [np.abs(d2qdz2),2*(1/1000)*abs(dqdz),(1/(1000**2))*pv]
            plt.bar(labels,values)
            plt.yscale("log")
            plt.show()

        sixpartcompare = False
        if sixpartcompare:
            labels = ["pvkv0","pvkh","pvkvb","fakebeta","skv0","skh","skvb"]
            fakesal = (1/2*f)*(10**-5)*(dqnotdx-x*beta*pv/(f*r))
            fakebeta = (1/2*f)*(10**-5)*dsdx
            values = [pvkv0*kvoscale,pvkh*khscale,pvkvb*kvbscale,fakesal,skvo*kvoscale,skh*khscale,skvb*kvbscale]
            print(np.format_float_scientific(pvkv0/fakebeta, unique=False, precision=3),np.format_float_scientific(skvo/fakesal, unique=False, precision=3))
            plt.bar(labels,np.abs(values))
            plt.yscale("log")
            plt.show()
            
        return np.asarray([-pvkv0/kvoscale,-pvkvb/kvbscale,-pvkh/khscale]),np.asarray([-skvo/kvoscale,-skvb/kvbscale,-skh/khscale])

def fetchWithFallback(surfaces,k,q,found,fallback=None):
    r = surfaces[k]["data"][q][found]
    if np.isnan(r) and fallback:
        return surfaces[k]["data"][q][fallback]
    else:
        return r
## file that generates the mixing terms of the Fq and Fs
## given a surfaces object, a depth and an index
def kterms(surfaces,k,found,params,debug=False,fallback=None):
    scales = params["scalecoeffs"]
    f = gsw.f(surfaces[k]["lats"][found])
    x = surfaces[k]["x"][found]
    y = surfaces[k]["y"][found]
    r = np.sqrt(x**2 + y**2)
    hx = fetchWithFallback(surfaces,k,"hx",found,fallback)
    hy = fetchWithFallback(surfaces,k,"hy",found,fallback)
    toph = fetchWithFallback(surfaces,k,"toph",found,fallback)
    both = fetchWithFallback(surfaces,k,"both",found,fallback)
    dsdx = fetchWithFallback(surfaces,k,"dsdx",found,fallback)
    dsdy = fetchWithFallback(surfaces,k,"dsdy",found,fallback)
    pres = fetchWithFallback(surfaces,k,"pres",found,fallback)
    alpha = fetchWithFallback(surfaces,k,"alpha",found,fallback)
    betaTherm = fetchWithFallback(surfaces,k,"beta",found,fallback)
    dsdz =  fetchWithFallback(surfaces,k,"dsdz",found,fallback)
    d2sdx2 =  fetchWithFallback(surfaces,k,"d2sdx2",found,fallback)
    d2sdy2 =  fetchWithFallback(surfaces,k,"d2sdy2",found,fallback)
    dalphadtheta = fetchWithFallback(surfaces,k,"dalphadtheta",found,fallback)
    dalphads = fetchWithFallback(surfaces,k,"dalphads",found,fallback)
    dalphadp = fetchWithFallback(surfaces,k,"dalphadp",found,fallback)
    dbetadp = fetchWithFallback(surfaces,k,"dbetadp",found,fallback)
    dbetads = fetchWithFallback(surfaces,k,"dbetads",found,fallback)
    dtdx = fetchWithFallback(surfaces,k,"dtdx",found,fallback)
    dtdy = fetchWithFallback(surfaces,k,"dtdy",found,fallback)
    dqnotdx = fetchWithFallback(surfaces,k,"dqnotdx",found,fallback)
    dqnotdy = fetchWithFallback(surfaces,k,"dqnotdy",found,fallback)
    dqdz = fetchWithFallback(surfaces,k,"dqdz",found,fallback)
    d2qdz2 = fetchWithFallback(surfaces,k,"d2qdz2",found,fallback)
    dpdx = fetchWithFallback(surfaces,k,"dpdx",found,fallback)
    dpdy = fetchWithFallback(surfaces,k,"dpdy",found,fallback)
    dqdx = fetchWithFallback(surfaces,k,"dqdx",found,fallback)
    dqdy = fetchWithFallback(surfaces,k,"dqdy",found,fallback)
    d2qdx2 = fetchWithFallback(surfaces,k,"d2qdx2",found,fallback)
    d2qdy2 = fetchWithFallback(surfaces,k,"d2qdy2",found,fallback)
    khp = fetchWithFallback(surfaces,k,"khp",found,fallback)
    dkhpdz = fetchWithFallback(surfaces,k,"khpdz",found,fallback)
    alphat = dalphadtheta+2*(alpha/betaTherm)*dalphads-(alpha**2/betaTherm**2)*dbetads
    alphap = dalphadp -(alpha/betaTherm)*dbetadp
    pv = fetchWithFallback(surfaces,k,"pv",found,fallback)
    doublets = fetchWithFallback(surfaces,k,"d2thetads2",found,fallback)
    CKVB = fetchWithFallback(surfaces,k,"CKVB",found,fallback)
    uref = fetchWithFallback(surfaces,k,"uref",found,fallback)
    vref = fetchWithFallback(surfaces,k,"vref",found,fallback)
    f = gsw.f(surfaces[k]["lats"][found])
    beta = calcBeta(surfaces[k]["lats"][found])
    isitnan = [alpha,betaTherm,dsdz,hx,hy,dsdx,dsdy,pres,d2sdx2,d2sdy2,\
              dalphadtheta,dalphads,dalphadp,dbetadp,dbetads,dtdx,dtdy,\
              dqnotdx,dqnotdy,dpdx,dpdy,alphat,alphap,pv,doublets,CKVB,\
              beta,dkhpdz,d2qdx2,d2qdy2,khp,toph,both,uref,vref]
    kvoscale = scales["kvo"]
    kvbscale = scales["kvb"]
    khscale  = scales["kh"]
    if (np.isnan(isitnan).any()):
        if debug:
            print("pres is nan: ",np.isnan(pres))
            print("hx is nan: ",np.isnan(hx))
            print("hy is nan: ",np.isnan(hy))
            print("x is nan: ",np.isnan(x))
            print("y is nan: ",np.isnan(y))
            print("something here is nan")
        return {},{}
    if not (np.isnan(isitnan).any()):
        pvkvb = (d2qdz2+2*(1/200.0)*dqdz+(1/(200.0**2))*pv)*CKVB
        pvkvo = d2qdz2
        if params["modelmixing"]:
            pvkh = (d2qdx2+d2qdy2)-2*(dqnotdx*dqdx+dqnotdy*dqdy)/pv - f*dkhpdz
        else:
            pvkh = (d2qdx2+d2qdy2)-2*(dqnotdx*dqdx+dqnotdy*dqdy)/pv -((1.0/both-1.0/toph)*f*khp)
        skvo = -alpha*f*(1/pv)*(dsdz**3)*doublets
        skvb = skvo*CKVB
        skhpart1 = (f/pv)*dsdz*(alphat*(dtdx**2 + dtdy**2)+alphap*(dtdx*dpdx+dtdy*dpdy))
        skhpart2 = (d2sdx2+d2sdy2)-2*(dqnotdx*dsdx + dqnotdy*dsdy)/pv
        skh = skhpart1 + skhpart2
        kvs = {"kvo":pvkvo*kvoscale,"kvb":pvkvb*kvbscale,"kh":pvkh*khscale}
        ks = {"kvo":skvo*kvoscale,"kvb":skvb*kvbscale,"kh":skh*khscale}
        return kvs,ks

def calcFQ(surfaces,k,found,scales,kpvs,distances):
    dqdx = fetchWithFallback(surfaces,k,"dqdx",found)
    dqdy = fetchWithFallback(surfaces,k,"dqdy",found)
    dqnotdx = fetchWithFallback(surfaces,k,"dqnotdx",found)
    dqnotdy = fetchWithFallback(surfaces,k,"dqnotdy",found)
    pv = fetchWithFallback(surfaces,k,"pv",found)

    missingpiece = -(2*(dqnotdx*dqdx+dqnotdy*dqdy)/pv)
    print("Qkvterm: ",surfaces[k]["data"]["diffkr"][found]*kpvs["kvo"])
    print("Qkhterm: ",surfaces[k]["data"]["kapredi"][found]*(kpvs["kh"]-missingpiece)+\
            surfaces[k]["data"]["kapgm"][found]*(missingpiece))
    FQ = surfaces[k]["data"]["diffkr"][found]*kpvs["kvo"] +\
            surfaces[k]["data"]["kapredi"][found]*(kpvs["kh"]-missingpiece)+\
            surfaces[k]["data"]["kapgm"][found]*(missingpiece)
    return FQ


def calcFS(surfaces,k,found,scales,ks,distances):
    FS = surfaces[k]["data"]["diffkr"][found]*ks["kvo"] +\
            surfaces[k]["data"]["kapredi"][found]*ks["kh"] 
    print("Skvterm: ",surfaces[k]["data"]["diffkr"][found]*ks["kvo"])
    print("Skhterm: ",surfaces[k]["data"]["kapredi"][found]*ks["kh"] )
    if np.isinf(-ks["kvo"]):
        print("kvo inf")
        print(ks)
    if np.isinf(-ks["kh"]):
        print("kh inf")
    return FS


