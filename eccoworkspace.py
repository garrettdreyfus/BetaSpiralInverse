import nstools
import inverttools
import eccotools
import interptools
import bathtools
import saloffset
import graph
import parametertools as ptools
import pickle
import sensitivity
import random
    
with open('data/eccoprofiles.pickle', 'rb') as outfile:
    profiles = pickle.load(outfile)

#graph.plotProfiles(profiles,"yo",data="t",depth="200")

#deepestindex=nstools.deepestProfile(profiles)
#graph.plotProfiles(profiles,"ecco",specialprofile=profiles[deepestindex])

#profiles = nstools.filterCruises(profiles,offsets.keys())
#profiles = saloffset.applyOffsets(profiles,offsets)

#profilechoice = random.choice(nstools.profileInBox(profiles,-180,180,85,90))
#profilechoice = nstools.getProfileById(profiles,profiles[deepestindex].eyed)

#surfaces = nstools.runPeerSearch(profiles,200,4000,200,profilechoice,1000)

#with open('data/eccosurfaces.pickle', 'wb') as outfile:
  #pickle.dump(surfaces, outfile)

#with open('data/eccoannotatedprofiles.pickle', 'wb') as outfile:
    #pickle.dump(profiles, outfile)

#with open('data/eccoannotatedprofiles.pickle', 'rb') as outfile:
    #profiles=pickle.load(outfile)

with open('data/eccosurfaces.pickle', 'rb') as outfile:
    surfaces=pickle.load(outfile)
###graph.graphSurfaces(surfaces,"pres")
#####graph.tsNeutralExplore(profiles)

surfaces = nstools.addDataToSurfaces(profiles,surfaces,2)
with open('data/eccosurfwithd.pickle', 'wb') as outfile:
    pickle.dump(surfaces, outfile)

#with open('data/eccosurfwithd.pickle', 'rb') as outfile:
    #surfaces=pickle.load(outfile)

nstools.surfaceDiagnostic(surfaces)

surfaces,neighbors,distances = interptools.interpolateSurfaces(surfaces,gaminterpolate=False)

surfaces = eccotools.addModelEccoMix(surfaces)
surfaces = eccotools.addModelEccoUV(surfaces)

surfaces = nstools.addHeight(surfaces)

surfaces = nstools.addHorizontalGrad(surfaces,neighbors,distances)
surfaces = nstools.addBathAndMask(surfaces,neighbors)
surfaces = nstools.addVerticalGrad(surfaces)
ptools.saveBathVarTermCache(surfaces,"data/bathVarecco.pickle")
surfaces = nstools.addK(surfaces,"data/bathVarecco.pickle")
with open('data/ready4inverseecco.pickle', 'wb') as outfile:
    pickle.dump([surfaces,neighbors,distances], outfile)


#with open('data/ready4inverseecco.pickle', 'rb') as outfile:
    #[surfaces,neighbors,distances]=pickle.load(outfile)
#for q in surfaces[1000]["data"].keys():
    #graph.graphSurfaces(surfaces,q,savepath="refpics/eccoallquantsnewgrid/",show=False)


#sensitivity.mixSens("coupled",surfaces,neighbors,distances,savepath="refpics/sensitivity/bathmasked")
#params = {"mixs":{"kvo":False,"kvb":False,"kh":False},"modelmixing":False}
#params = {"mixs":{"kvo":False,"kvb":False,"kh":False},"debug":False,"modelmixing":True}
#sensitivity.conditionErrorRefLevel("coupled",surfaces,neighbors,distances,params=params)
##sensitivity.conditionError("coupled",staggeredsurfaces,neighbors,distances)

params = {"reflevel":1000,"upperbound":1000,"lowerbound":1400,"mixs":{"kvo":False,"kvb":False,"kh":False},\
        "debug":True,"modelmixing":True}

out = inverttools.invert("coupled",surfaces,neighbors,distances,params=params)
inv = out["surfaces"]
inv = nstools.streamFuncToUV(inv,neighbors,distances)
graph.graphVectorField(inv,"uabs","vabs","diffkr")
#graph.graphVectorField(inv,"uabs","vabs","kapredi")
#graph.graphVectorField(inv,"uabs","vabs","kapgm")
##surfaces = eccotools.addSSHToSurface(surfaces)
##graph.graphSurfaces(surfaces,"ssh")
##graph.graphVectorField(inv,"u","v","pv",savepath="refpics/vectorfields/geostroph/")
