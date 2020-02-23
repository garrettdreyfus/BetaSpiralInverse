from regionlib import brasil
import os
import graph
import nstools
import pickle
import interptools
import inverttools
from functools import partial

profiles = brasil.extractArgoProfiles(os.path.abspath("data/brasilargonc"))
graph.plotProfiles(brasil,profiles,"profiles")
#profiles = profiles + brasil.extractWoceProfiles(os.path.abspath("data/brasilnc/"))
#profiles = brasil.extractWoceProfiles(os.path.abspath("data/brasilnc/"))

#with open('data/argoandwoce.pickle', 'wb') as outfile:
    #pickle.dump(profiles,outfile)
#with open('data/argoandwoce.pickle', 'rb') as outfile:
    #profiles = pickle.load(outfile)


#profilechoice = nstools.profileInBox(profiles,-31.16,-14.91,-24.48,-7.98,6000)
#profilechoice = profilechoice[0]
#preinterpsurfaces = nstools.runPeerSearch(profiles,range(100,6000,200),profilechoice,False,10**10)

#with open('data/annotatednepbprofilessingleref.pickle', 'wb') as outfile:
    #pickle.dump([preinterpsurfaces,profiles],outfile)
#with open('data/annotatednepbprofilessingleref.pickle', 'rb') as outfile:
    #preinterpsurfaces,profiles = pickle.load(outfile)

#surfaces = nstools.addDataToSurfaces(brasil,profiles,preinterpsurfaces)

#with open('data/annotatednepbprofilessingleref.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,profiles],outfile)
with open('data/annotatednepbprofilessingleref.pickle', 'rb') as outfile:
    preinterpsurfaces,profiles = pickle.load(outfile)

#surfaces,neighbors,distances = interptools.interpolateSurfaces(brasil,preinterpsurfaces,\
        #interpmethod="gam",smart=False,coord="latlon")

#with open('data/interpedbrasil.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,neighbors,distances], outfile)
with open('data/interpedbrasil.pickle', 'rb') as outfile:
    [surfaces,neighbors,distances] = pickle.load(outfile)
graph.graphSurfaces(brasil,surfaces,"pres",secondsurface=preinterpsurfaces,\
        contour=True,show=False,savepath = "refpics/surfaces/prescompbrasilcropped/")
#surfaces = nstools.addParametersToSurfaces(brasil,surfaces,\
        #neighbors,distances)
#nstools.inverseReady(surfaces)


#with open('data/interpedbrasil.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,neighbors,distances], outfile)
#with open('data/interpedbrasil.pickle', 'rb') as outfile:
    #[surfaces,neighbors,distances] = pickle.load(outfile)

#params = {"reflevel":1700,"upperbound":1000,"lowerbound":4000,\
        #"mixs":{"kvo":True,"kvb":True,"kh":True},"debug":False,\
        #"3point":True,"edgeguard":True}


#out= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)

#with open('data/inverseoutbrasil.pickle', 'wb') as outfile:
    #pickle.dump([out,neighbors,distances], outfile)
with open('data/inverseoutbrasil.pickle', 'rb') as outfile:
    [out,neighbors,distances] = pickle.load(outfile)

inv = nstools.streamFuncToUV(out["surfaces"],neighbors,distances)
#graph.saveAllQuants(brasil,inv,"refpics/surfaces/brasilandargocropped/")

#for lat in range(-30,-2,5):
    #graph.northSouthTransect(inv,"vabs",lat=lat,savepath="refpics/transects/")
#for lon in range(-35,-12,5):
    #graph.northSouthTransect(inv,"uabs",lon=lon,savepath="refpics/transects/")


#graph.graphSurfaces(brasil,inv,"psinew",stds=3,\
        #show=False, savepath="refpics/surfaces/inversesolutionnewcenter/",\
        #centerfunction = partial(nstools.domainMedianStdev,-30,-20,-25,-10))
#graph.graphVectorField(brasil,inv,"uabs","vabs","psi",\
        #metadata=out["metadata"],\
        #transform=False,show=False,
        #savepath="refpics/vectorfields/brasilmixcropped/4pointpv/")

