from regionlib import brasil
import os
import graph
import nstools
import pickle
import interptools
import inverttools
from functools import partial
import parametertools as ptools

#profiles = brasil.extractArgoProfiles(os.path.abspath("data/brasilargonc"))
#graph.plotProfiles(brasil,profiles,"profiles")
#profiles = profiles + brasil.extractWoceProfiles(os.path.abspath("data/brasilnc/"))
#profiles = brasil.extractWoceProfiles(os.path.abspath("data/brasilnc/"))

#with open('data/argoandwoce.pickle', 'wb') as outfile:
    #pickle.dump(profiles,outfile)
#with open('data/argoandwoce.pickle', 'rb') as outfile:
    #profiles = pickle.load(outfile)

#profilechoice = nstools.profileInBox(profiles,-40,-20,-31,-20,5000)
#profilechoice = profilechoice[0]
#print(profilechoice.lat,profilechoice.lon)
#preinterpsurfaces = nstools.runPeerSearch(profiles,range(100,6000,200),profilechoice,False,10**4)
#gammavals = [25.875,26.51,26.862,27.158,27.3605,27.526,27.6575,27.7825,27.8700, \
        #27.9275,27.965,27.99,28.015,28.03,28.0475,28.0625,28.08,28.108,28.136,28.164,28.2,28.33,28.36]
##print(range(100,6000,200)[len(gammavals)-1])
##preinterpsurfaces = nstools.runPeerSearch(profiles,range(100,4300,200),profilechoice,False,10**4,gammas=gammavals)
##preinterpsurfaces = nstools.runPeerSearch(profiles,range(100,4300,200),profilechoice,True,10**5)

#with open('data/annotatedbrasil.pickle', 'wb') as outfile:
    #pickle.dump([preinterpsurfaces,profiles],outfile)

#with open('data/annotatedbrasil.pickle', 'rb') as outfile:
    #preinterpsurfaces,profiles = pickle.load(outfile)
#surfaces = nstools.addDataToSurfaces(brasil,profiles,preinterpsurfaces)
#nstools.addGammaN(surfaces)

##graph.graphSurfaces(brasil,surfaces,"pres",show=False,\
        ##savepath="refpics/surfaces/ctdsouthernexension")

#with open('data/annotatednepbprofilessingleref.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,profiles],outfile)
#with open('data/annotatednepbprofilessingleref.pickle', 'rb') as outfile:
    #preinterpsurfaces,profiles = pickle.load(outfile)

#surfaces,neighbors,distances = interptools.interpolateSurfaces(brasil,preinterpsurfaces,\
        #interpmethod="gam",smart=False,coord="latlon")




#with open('data/interpedbrasil.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,neighbors,distances], outfile)
#with open('data/interpedbrasil.pickle', 'rb') as outfile:
    #[surfaces,neighbors,distances] = pickle.load(outfile)

#surfaces = nstools.addParametersToSurfaces(brasil,surfaces,\
        #neighbors,distances)
#nstools.inverseReady(surfaces)


#with open('data/interpedbrasil.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,neighbors,distances], outfile)
#with open('data/interpedbrasil.pickle', 'rb') as outfile:
    #[surfaces,neighbors,distances] = pickle.load(outfile)

##surfaces = nstools.addGammaN(surfaces)
##surfaces = nstools.neutralityError(surfaces)
##graph.graphSurfaces(brasil,surfaces,"nserror",stds=2,savepath="refpics/surfaces/neutralityerror/spread/",\
        ##boundfunc=graph.minmaxBound,show=False)
#params = {"reflevel":1700,"upperbound":1000,"lowerbound":4000,\
        #"mixs":{"kvo":True,"kvb":True,"kh":True},"debug":False,\
        #"3point":True,"edgeguard":True}


#out= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)

#with open('data/inverseoutbrasil.pickle', 'wb') as outfile:
    #pickle.dump([out,neighbors,distances], outfile)
with open('data/inverseoutbrasil.pickle', 'rb') as outfile:
    [out,neighbors,distances] = pickle.load(outfile)

inv = nstools.streamFuncToUV(out["surfaces"],neighbors,distances)
inv = ptools.calcFRho(inv)
#nstools.addGammaN(inv)
#inv = nstools.compareGammaN(inv)

#graph.saveAllQuants(brasil,inv,"refpics/surfaces/brasilandargocropped/")

#for lat in range(-30,-2,5):
    #graph.northSouthTransect(inv,"uabs",lat=lat,savepath="refpics/transects/noiseexpiriment/"+str(n)+"stdev"+str(i)+"/")
#for lon in range(-35,-12,5):
    #graph.northSouthTransect(inv,"vabs",lon=lon,savepath="refpics/transects/noiseexpiriment/"+str(n)+"stdev"+str(i)+"/")
#graph.graphSurfaces(brasil,inv,"gamma",stds=2,savepath="refpics/surfaces/gammaderivedsurface/",show=False)
#graph.graphVectorField(brasil,inv,"uabs","vabs","pres",\
        #transform=False)
graph.graphVectorField(brasil,inv,"uabs","vabs","pres",\
        metadata=out["metadata"],\
        transform=False,show=False,
        savepath="refpics/vectorfields/brasilwoceargo/",scale=0.1)

