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
def bathVarTerm(lat,lon,region):
    d=bathtools.bathBox(lat,lon)
    dnot = 750
    print("box",d)
    return (np.var(d)/dnot)**(0.25),np.mean(d)

#one of the mixing terms is a mess and this calculates that
def calculateKHP(staggered,k,index):
    dalphadtheta = staggered[k]["data"]["dalphadtheta"][index]
    dalphadp = staggered[k]["data"]["dalphadp"][index]
    dalphads = staggered[k]["data"]["dalphads"][index]
    alpha = staggered[k]["data"]["alpha"][index]
    alphat = staggered[k]["data"]["dalphadtheta"][index]
    alphap = staggered[k]["data"]["dalphadp"][index]
    magct = staggered[k]["data"]["dtdx"][index]**2 + staggered[k]["data"]["dtdy"][index]**2
    cdotp = staggered[k]["data"]["dtdx"][index]*staggered[k]["data"]["dpdx"][index]+staggered[k]["data"]["dtdy"][index]*staggered[k]["data"]["dpdy"][index]
    return alphat*magct+alphap*cdotp


## this stuff takes a long time to calculate so this
## throws everything into a pickle for later
def saveBathVarTermCache(surfaces,outfilename,region):
    bathVar = {}
    coords = []
    for k in surfaces.keys():
        coords += zip(surfaces[k]["lons"],surfaces[k]["lats"])
    coords = set(coords)
    for p in  Bar("var coord :").iter(coords):
        if not (np.isnan(p[0]) and np.isnan(p[1])):
            d=bathtools.bathBox(p[1],p[0])
            dnot = 750
            bathVar[p] = ((np.var(d)/dnot)**(0.25),np.mean(d))
    with open(outfilename, 'wb') as outfile:
        pickle.dump(bathVar, outfile)


## extracts bathVar data from pickle but only does so once
def bathVarTermCache(lat,lon,filename):
    if not hasattr(bathVarTermCache,"p") :
        bathVarTermCache.p = pickle.load(open(filename, 'rb'))

        print("loaded bath variance file")
    if not bathVarTermCache.p:
        bathVarTermCache.p = pickle.load(open(filename, 'rb'))

        print("loaded bath variance file")
        
    return bathVarTermCache.p[(lon,lat)]


## calulate the Kvb bathvar coefficient
def Kv(lat,lon,pv,pres,cachename=None):
    if cachename:
        bVT,mean = bathVarTermCache(lat,lon,cachename) 
    else:
        bVT,mean = bathVarTerm(lat,lon)
    return bVT*np.exp(-(abs(mean)-abs(pres))/1000),bVT,np.exp(-(abs(mean)-abs(pres))/5500)

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
    alphat = surfaces[k]["data"]["dalphadtheta"][found] 
    alphap = surfaces[k]["data"]["dalphadp"][found]
    pv =  surfaces[k]["data"]["pv"][found] 
    doublets =  surfaces[k]["data"]["d2thetads2"][found] 
    CKVB =  surfaces[k]["data"]["CKVB"][found] 
    f = gsw.f(surfaces[k]["lats"][found])
    beta = calcBeta(surfaces[k]["lats"][found])
    isitnan = [alpha,betaTherm,dsdz,hx,hy,dsdx,dsdy,pres,d2sdx2,d2sdy2,\
              dalphadtheta,dalphads,dalphadp,dtdx,dtdy,\
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
            values = [np.abs(d2qdz2),alpha,f,1/pv,np.abs(dsdz)**3,np.abs(doublets),np.abs(skvo)]
            plt.bar(labels,values)
            plt.yscale("log")
            plt.show()

        khbreakdown=False
        if khbreakdown:
            #(d2qdx2+d2qdy2)-2*(dqnotdx*dqdx+dqnotdy*dqdy)/pv -f*khpdz
            labels = ["d2qdx2+d2qdy2","2*(dqnotdx*dqdx+dqnotdy*dqdy)/pv","f*khpdz","f","khpdz"]
            values = [np.abs(d2qdx2+d2qdy2),np.abs(2*(dqnotdx*dqdx+dqnotdy*dqdy)/pv),np.abs(f*khpdz),f,abs(khpdz)]
            plt.bar(labels,values)
            plt.yscale("log")
            plt.show()

        kvbbreakdown=False
        if kvbbreakdown:
            labels = ["d2qdz2","dqdz","pv"]
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
            
        return np.asarray([-pvkv0*kvoscale,-pvkvb*kvbscale,-pvkh*khscale]),np.asarray([-skvo*kvoscale,-skvb*kvbscale,-skh*khscale])

def fetchWithFallback(surfaces,k,q,found,fallback=None):
    r = surfaces[k]["data"][q][found]
    if np.isnan(r) and fallback:
        return surfaces[k]["data"][q][fallback]
    else:
        return r
## file that generates the mixing terms of the Fq and Fs
## given a surfaces object, a depth and an index
def kterms(surfaces,k,found,params,fallback=None):
    debug = params["debug"]
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
    alphat = fetchWithFallback(surfaces,k,"dalphadtheta",found,fallback)
    alphap = fetchWithFallback(surfaces,k,"dalphadp",found,fallback)
    pv = fetchWithFallback(surfaces,k,"pv",found,fallback)
    doublets = fetchWithFallback(surfaces,k,"d2thetads2",found,fallback)
    CKVB = fetchWithFallback(surfaces,k,"CKVB",found,fallback)
    uref = fetchWithFallback(surfaces,k,"uref",found,fallback)
    vref = fetchWithFallback(surfaces,k,"vref",found,fallback)
    f = gsw.f(surfaces[k]["lats"][found])
    beta = calcBeta(surfaces[k]["lats"][found])
    isitnan = [alpha,betaTherm,dsdz,hx,hy,dsdx,dsdy,pres,d2sdx2,d2sdy2,\
              dalphadtheta,dalphads,dalphadp,dtdx,dtdy,\
              dqnotdx,dqnotdy,dpdx,dpdy,alphat,alphap,pv,doublets,CKVB,\
              beta,d2qdx2,d2qdy2,d2qdz2,khp,toph,both,uref,vref]
    isitnanstr = np.asarray(["alpha","betaTherm","dsdz","hx","hy","dsdx","dsdy","pres","d2sdx2","d2sdy2",\
              "dalphadtheta","dalphads","dalphadp","dtdx","dtdy",\
              "dqnotdx","dqnotdy","dpdx","dpdy","alphat","alphap","pv","doublets","CKVB",\
              "beta","d2qdx2","d2qdy2","d2qdz2","khp","toph","both","uref","vref"])
    kvoscale = scales["kvo"]
    kvbscale = scales["kvb"]
    khscale  = scales["kh"]
    if (np.isnan(isitnan).any()):
        if debug:
            print(isitnanstr[np.isnan(isitnan)])
            print("something here is nan")
        return {},{}
    if not (np.isnan(isitnan).any()):
        pvkvb = (d2qdz2+2*(-CKVB/5500.0)*dqdz+(CKVB/(5500.0**2))*pv)*CKVB
        #print("pvkvb: ",pvkvb," : ",d2qdz2,dqdz,pv,CKVB,)
        pvkvo = d2qdz2
        #print("pvkvo: ",pvkvo)
        if params["modelmixing"]:
            #print("NO")
            pvkh = (d2qdx2+d2qdy2)-2*(dqnotdx*dqdx+dqnotdy*dqdy)/pv - f*dkhpdz/(surfaces[k]["data"]["kapredi"][found])
        else:
            if pv == 0:
                print("wut pv")
            if toph==0:
                print("toph")
            if both==0:
                print("both")
            pvkh = (d2qdx2+d2qdy2)-2*(dqnotdx*dqdx+dqnotdy*dqdy)/pv -((1.0/both-1.0/toph)*f*khp)
            #print("pvkvh: ",pvkh)
        skvo = -alpha*f*(1/pv)*(dsdz**3)*doublets
        skvb = skvo*CKVB
        skhpart1 = (f/pv)*dsdz*(alphat*(dtdx**2 + dtdy**2)+alphap*(dtdx*dpdx+dtdy*dpdy))
        skhpart2 = (d2sdx2+d2sdy2)-2*(dqnotdx*dsdx + dqnotdy*dsdy)/pv
        skh = skhpart1 + skhpart2
        kvs = {"kvo":pvkvo*kvoscale,"kvb":pvkvb*kvbscale,"kh":pvkh*khscale}
        ks = {"kvo":skvo*kvoscale,"kvb":skvb*kvbscale,"kh":skh*khscale}
        return kvs,ks

def calcFQ(surfaces,k,found,scales,kpvs,distances,debug=False):
    dqdx = fetchWithFallback(surfaces,k,"dqdx",found)
    dqdy = fetchWithFallback(surfaces,k,"dqdy",found)
    dqdz = fetchWithFallback(surfaces,k,"dqdz",found)
    d2qdz2 = fetchWithFallback(surfaces,k,"d2qdz2",found)
    d2qdx2 = fetchWithFallback(surfaces,k,"d2qdx2",found)
    d2qdy2 = fetchWithFallback(surfaces,k,"d2qdy2",found)
    dqnotdx = fetchWithFallback(surfaces,k,"dqnotdx",found)
    dqnotdy = fetchWithFallback(surfaces,k,"dqnotdy",found)
    ddiffkrdz = fetchWithFallback(surfaces,k,"ddiffkrdz",found)
    d2diffkrdz2 = fetchWithFallback(surfaces,k,"d2diffkrdz2",found)
    khpdz = fetchWithFallback(surfaces,k,"khpdz",found)
    f = gsw.f(surfaces[k]["lats"][found])
    pv = fetchWithFallback(surfaces,k,"pv",found)

    missingpiece = -(2*(dqnotdx*dqdx+dqnotdy*dqdy)/pv)*scales["kh"]
    if debug:
        print("#"*5)
        print("diffkr: ",surfaces[k]["data"]["diffkr"][found])
        print("q: ", pv)
        print("Qkvterm part 1: ",surfaces[k]["data"]["diffkr"][found]*d2qdz2)
        print("Qkvterm part 2: ",2*ddiffkrdz*dqdz )
        print("Qkvterm part 3: ",d2diffkrdz2*pv )
        print("Qkhterm part 1 : ",f*khpdz)
        print("Qkhterm part 2 : ",surfaces[k]["data"]["kapgm"][found]*(missingpiece))
        print("kapgm: ",surfaces[k]["data"]["kapgm"][found])
        print("kapredi: ",surfaces[k]["data"]["kapredi"][found])
    FQ = surfaces[k]["data"]["diffkr"][found]*d2qdz2 +\
            2*ddiffkrdz*dqdz+ d2diffkrdz2*pv +  \
            f*khpdz+\
            surfaces[k]["data"]["kapgm"][found]*(missingpiece) +\
            surfaces[k]["data"]["kapredi"][found]*(d2qdx2+d2qdy2)

    if np.isnan(FQ):
        FQ=0
    surfaces[k]["data"]["FQ"][found] = FQ
    return FQ


def calcFS(surfaces,k,found,scales,ks,distances,debug=False):
    dsdy = fetchWithFallback(surfaces,k,"dqdy",found)
    dsdx = fetchWithFallback(surfaces,k,"dqdx",found)
    d2sdy2 = fetchWithFallback(surfaces,k,"d2sdy2",found)
    d2sdx2 = fetchWithFallback(surfaces,k,"d2sdx2",found)
    dqnotdx = fetchWithFallback(surfaces,k,"dqnotdx",found)
    dqnotdy = fetchWithFallback(surfaces,k,"dqnotdy",found)
    pv = fetchWithFallback(surfaces,k,"pv",found)
    missingpiece = -(2*(dqnotdx*dsdx+dqnotdy*dsdy)/pv)*scales["kh"]
    FS = surfaces[k]["data"]["diffkr"][found]*ks["kvo"] +\
            surfaces[k]["data"]["kapredi"][found]*(ks["kh"]-missingpiece) +\
            surfaces[k]["data"]["kapgm"][found]*missingpiece
    if debug or True:
        #print("diffkrterm: ",surfaces[k]["data"]["diffkr"][found]*ks["kvo"])
        #print("kapredi term: ",surfaces[k]["data"]["kapredi"][found]*(ks["kh"]-missingpiece))
        #print("kapgm term: ",surfaces[k]["data"]["kapgm"][found]*missingpiece)
        largest = np.max(np.abs([surfaces[k]["data"]["diffkr"][found]*ks["kvo"],surfaces[k]["data"]["kapredi"][found]*(ks["kh"]-missingpiece),surfaces[k]["data"]["kapgm"][found]*missingpiece]))
        ratio  = surfaces[k]["data"]["kapredi"][found]*(d2sdx2+d2sdy2)/largest
        if ratio > 100:
            print("how big",ratio, " where: ",k,surfaces[k]["lats"][found],surfaces[k]["lons"][found])
    surfaces[k]["data"]["FS"][found] = FS
    if np.isinf(-ks["kvo"]):
        print("kvo inf")
        #print(ks)
    if np.isinf(-ks["kh"]):
        print("kh inf")
    return FS

def calcFRho(surfaces):
    depths = sorted(surfaces.keys())
    for j in depths:
        surfaces[j]["data"]["e"] =  np.full(len(surfaces[j]["lons"]),np.nan)
    for j in range(1,len(depths)-1):
        k = depths[j]
        for index in range(len(surfaces[k]["lons"])):
            eyed = int(surfaces[k]["ids"][index])
            foundbelow = np.where(np.asarray(surfaces[depths[j+1]]["ids"])==eyed)
            found = index
            foundabove = np.where(np.asarray(surfaces[depths[j-1]]["ids"])==eyed)

            if eyed != -999 and len(foundbelow)!=0 and len(foundbelow[0]) != 0 and len(foundabove)!=0 and len(foundabove[0]) != 0:
                foundbelow = foundbelow[0][0]
                foundabove = foundabove[0][0]
                abovekv = surfaces[depths[j-1]]["data"]["kvo"][foundabove] + surfaces[depths[j-1]]["data"]["CKVB"][foundabove]*surfaces[depths[j-1]]["data"]["kvb"][foundabove]
                belowkv = surfaces[depths[j+1]]["data"]["kvo"][foundbelow] + surfaces[depths[j+1]]["data"]["CKVB"][foundbelow]*surfaces[depths[j+1]]["data"]["kvb"][foundbelow]
                firstterm = \
                        -(abovekv*surfaces[depths[j-1]]["data"]["drhodz"][foundabove] \
                        - belowkv*surfaces[depths[j+1]]["data"]["drhodz"][foundbelow])/\
                        (surfaces[depths[j-1]]["data"]["pres"][foundabove]-surfaces[depths[j+1]]["data"]["pres"][foundbelow])

                rhonot = 1025.0
                khp = surfaces[k]["data"]["khp"][index]
                kh = surfaces[k]["data"]["kh"][index]
                secondterm = rhonot*kh*khp
                e = (firstterm + secondterm)/surfaces[k]["data"]["drhodz"][index]
                surfaces[k]["data"]["e"][index]=e
    return surfaces


