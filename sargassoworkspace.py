from regionlib import sargasso
import graph
import nstools

profiles = sargasso.extractArgoProfiles("data/sargassonc")
graph.plotProfiles(sargasso,profiles,"sargasso")

profilechoice = profiles[nstools.deepestProfile(profiles)]
preinterpsurfaces = nstools.runPeerSearch(profiles,range(100,3000,200),profilechoice,False,10**10)

with open('data/annotatedsargassopbprofilessingleref.pickle', 'wb') as outfile:
    pickle.dump([preinterpsurfaces,profiles],outfile)
with open('data/annotatedsargassoprofilessingleref.pickle', 'rb') as outfile:
    preinterpsurfaces,profiles = pickle.load(outfile)

surfaces = nstools.addDataToSurfaces(sargasso,profiles,preinterpsurfaces)

with open('data/annotatedsargassoprofilessingleref.pickle', 'wb') as outfile:
    pickle.dump([surfaces,profiles],outfile)
with open('data/annotatedsargassoprofilessingleref.pickle', 'rb') as outfile:
    preinterpsurfaces,profiles = pickle.load(outfile)

surfaces,neighbors,distances = interptools.interpolateSurfaces(sargasso,preinterpsurfaces,\
        interpmethod="gam",smart=False,coord="latlon")

with open('data/interpedbrasil.pickle', 'wb') as outfile:
    pickle.dump([surfaces,neighbors,distances], outfile)
with open('data/interpedbrasil.pickle', 'rb') as outfile:
    [surfaces,neighbors,distances] = pickle.load(outfile)

surfaces = nstools.addParametersToSurfaces(brasil,surfaces,\
        neighbors,distances)
nstools.inverseReady(surfaces)


with open('data/interpedbrasil.pickle', 'wb') as outfile:
    pickle.dump([surfaces,neighbors,distances], outfile)
with open('data/interpedbrasil.pickle', 'rb') as outfile:
    [surfaces,neighbors,distances] = pickle.load(outfile)

params = {"reflevel":1700,"upperbound":1000,"lowerbound":4000,\
        "mixs":{"kvo":True,"kvb":True,"kh":True},"debug":False,\
        "3point":True,"edgeguard":True}


out= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)

with open('data/inverseoutbrasil.pickle', 'wb') as outfile:
    pickle.dump([out,neighbors,distances], outfile)
with open('data/inverseoutbrasil.pickle', 'rb') as outfile:
    [out,neighbors,distances] = pickle.load(outfile)

inv = nstools.streamFuncToUV(out["surfaces"],neighbors,distances)
graph.saveAllQuants(brasil,inv,"refpics/surfaces/brasilandargocropped/")

for lat in range(-30,-2,5):
    graph.northSouthTransect(inv,"vabs",lat=lat,savepath="refpics/transects/")
for lon in range(-35,-12,5):
    graph.northSouthTransect(inv,"uabs",lon=lon,savepath="refpics/transects/")


graph.graphSurfaces(brasil,inv,"psinew",stds=3,\
        show=False, savepath="refpics/surfaces/inversesolutionnewcenter/",\
        centerfunction = partial(nstools.domainMedianStdev,-30,-20,-25,-10))
graph.graphVectorField(brasil,inv,"uabs","vabs","psi",\
        metadata=out["metadata"],\
        transform=False,show=False,
        savepath="refpics/vectorfields/brasilmixcropped/4pointpv/")

