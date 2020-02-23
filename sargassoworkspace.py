from regionlib import sargasso
import graph
import nstools
import pickle
import interptools
import inverttools

#profiles = sargasso.extractArgoProfiles("data/sargassonc")
#graph.plotProfiles(sargasso,profiles,"sargasso")

#profilechoice = profiles[nstools.deepestProfile(profiles)]
#preinterpsurfaces = nstools.runPeerSearch(profiles,range(100,3000,200),profilechoice,False)

#with open('data/annotatedsargassopbprofilessingleref.pickle', 'wb') as outfile:
    #pickle.dump([preinterpsurfaces,profiles],outfile)
#with open('data/annotatedsargassopbprofilessingleref.pickle', 'rb') as outfile:
    #preinterpsurfaces,profiles = pickle.load(outfile)

#surfaces = nstools.addDataToSurfaces(sargasso,profiles,preinterpsurfaces)

#with open('data/annotatedsargassoprofilessingleref.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,profiles],outfile)
with open('data/annotatedsargassoprofilessingleref.pickle', 'rb') as outfile:
    preinterpsurfaces,profiles = pickle.load(outfile)

surfaces,neighbors,distances = interptools.interpolateSurfaces(sargasso,preinterpsurfaces,\
        interpmethod="gam",smart=False,coord="latlon")


with open('data/interpedsargasso.pickle', 'wb') as outfile:
    pickle.dump([surfaces,neighbors,distances], outfile)
with open('data/interpedsargasso.pickle', 'rb') as outfile:
    [surfaces,neighbors,distances] = pickle.load(outfile)
    
graph.graphSurfaces(sargasso,surfaces,"pres",secondsurface=preinterpsurfaces,contour=True)

surfaces = nstools.addParametersToSurfaces(sargasso,surfaces,\
        neighbors,distances)
nstools.inverseReady(surfaces)


with open('data/interpedbrasil.pickle', 'wb') as outfile:
    pickle.dump([surfaces,neighbors,distances], outfile)
with open('data/interpedbrasil.pickle', 'rb') as outfile:
    [surfaces,neighbors,distances] = pickle.load(outfile)

params = {"reflevel":1700,"upperbound":1000,"lowerbound":2000,\
        "mixs":{"kvo":True,"kvb":True,"kh":True},"debug":False,\
        "3point":True,"edgeguard":True}


out= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)

with open('data/inverseoutbrasil.pickle', 'wb') as outfile:
    pickle.dump([out,neighbors,distances], outfile)
with open('data/inverseoutbrasil.pickle', 'rb') as outfile:
    [out,neighbors,distances] = pickle.load(outfile)

inv = nstools.streamFuncToUV(out["surfaces"],neighbors,distances)

#graph.graphSurfaces(sargasso,inv,"psinew",stds=3,\
        #show=False, savepath="refpics/surfaces/inversesolutionsargasso/",)
graph.graphVectorField(sargasso,inv,"uabs","vabs","psinew",\
        metadata=out["metadata"],\
        transform=False,show=False,scale=1,
        savepath="refpics/vectorfields/sargasso/")

