from regionlib import brasil
import os
import graph
import nstools
import pickle
import regions
import interptools
import inverttools

#profiles = brasil.extractProfiles(os.path.abspath("data/brasilnc/"))


#profilechoice = nstools.profileInBox(profiles,-31.16,-14.91,-24.48,-7.98,4000)
#profilechoice = profilechoice[0]
#preinterpsurfaces = nstools.runPeerSearch(profiles,range(100,6000,200),profilechoice,False,10**10)

#with open('data/annotatednepbprofilessingleref.pickle', 'wb') as outfile:
    #pickle.dump([preinterpsurfaces,profiles],outfile)
#with open('data/annotatednepbprofilessingleref.pickle', 'rb') as outfile:
    #preinterpsurfaces,profiles = pickle.load(outfile)

#surfaces = nstools.addDataToSurfaces(profiles,preinterpsurfaces,region=regions.brasil)

#with open('data/annotatednepbprofilessingleref.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,profiles],outfile)
#with open('data/annotatednepbprofilessingleref.pickle', 'rb') as outfile:
    #preinterpsurfaces,profiles = pickle.load(outfile)

#surfaces,neighbors,distances = interptools.interpolateSurfaces(preinterpsurfaces,\
        #regions.brasil,interpmethod="gam",smart=False,coord="latlon")

#with open('data/interpedbrasil.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,neighbors,distances], outfile)
#with open('data/interpedbrasil.pickle', 'rb') as outfile:
    #[surfaces,neighbors,distances] = pickle.load(outfile)

#surfaces = nstools.addParametersToSurfaces(surfaces,\
        #neighbors,distances,regions.brasil,[])
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

#graph.northSouthTransect(inv,"vabs",lat=-31)
#graph.northSouthTransect(inv,"uabs",lon=-25)
#graph.graphSurfaces(regions.brasil,inv,"psinew",stds=0.6,\
        #show=False, savepath="refpics/surfaces/inversesolution/")
#graph.graphVectorField(regions.brasil,inv,"uabs","vabs","pv",\
        #metadata=out["metadata"],\
        #transform=False,show=False,
        #savepath="refpics/vectorfields/brasilmix/4pointpv/")

