import numpy as np
import bathtools
import gsw
from progress.bar import Bar
import pickle

def bathVarTerm(lat,lon):
    d=bathtools.bathBox(lat,lon)
    dnot = 750
    return (np.var(d)/dnot)**(0.25)

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
    f = gsw.f(lat)
    fnot = gsw.f(30)
    Nnot = 5.2*(10**-3)
    N = np.sqrt((9.8*pv)/f)
    part1 = (f*(np.cosh(N/f)**-1))/(fnot*(np.cosh(Nnot/fnot)**-1))
    if cachename:
       bVT = bathVarTermCache(lat,lon,cachename) 
    else:
        bVT = bathVarTerm(lat,lon)
    part2 = bVT*np.exp(-(abs(bathtools.searchBath(lat,lon))-abs(pres))/5500)
    return (part1,part1*part2)

