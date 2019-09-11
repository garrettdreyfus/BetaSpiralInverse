import nstools
import numpy as np
import inverttools
import eccotools
import interptools
import bathtools
import saloffset
import graph
import matplotlib.pyplot as plt
import parametertools as ptools
import pickle
import sensitivity
import random
    
with open('data/ecconepbprofiles.pickle', 'rb') as outfile:
    profiles = pickle.load(outfile)


deepestindex=nstools.deepestProfile(profiles)
profilechoice = profiles[deepestindex]

#surfaces = nstools.runPeerSearch(profiles,200,4000,200,profilechoice,1000)

#with open('data/ecconepbsurfaces.pickle', 'wb') as outfile:
    #pickle.dump(surfaces,outfile)
#with open('data/ecconepbsurfaces.pickle', 'rb') as outfile:
    #surfaces=pickle.load(outfile)
#graph.graphSurfaces(surfaces,"pres",region="nepb",savepath="refpics/surfaces/ecconepbns/",show=False)
#########graph.tsNeutralExplore(profiles)

#surfaces = nstools.addDataToSurfaces(profiles,surfaces,2)

#with open('data/ecconepbsurfaceswithd.pickle', 'wb') as outfile:
    #pickle.dump(surfaces,outfile)

#with open('data/ecconepbsurfaceswithd.pickle', 'rb') as outfile:
    #surfaces=pickle.load(outfile)

####nstools.surfaceDiagnostic(surfaces)

#surfaces,neighbors,distances = interptools.interpolateSurfaces(surfaces,gaminterpolate=False,fixedgrid="nepb")

#surfaces = eccotools.addModelEccoMix(surfaces,"NEPB")
#surfaces = eccotools.addModelEccoUV(surfaces,"NEPB")

#with open('data/preparameter.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,neighbors,distances], outfile)

with open('data/preparameter.pickle', 'rb') as outfile:
    [surfaces,neighbors,distances ] = pickle.load(outfile)

surfaces = nstools.addParametersToSurfaces(surfaces,neighbors,distances)

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
        "debug":False,"modelmixing":False}

out = inverttools.invert("coupled",surfaces,neighbors,distances,params=params)
inv = out["surfaces"]
inv = nstools.streamFuncToUV(inv,neighbors,distances)

#surfs=[]
#eccovelsurfaces=[]
#velsurfaces=[]
#for k in surfaces.keys():
    #print("#"*5,k,"#"*5) 
    #eccovels=[]
    #vels = []
    #for l in range(len(surfaces[k]["data"]["knownu"])):
        #eccovels.append(np.sqrt(surfaces[k]["data"]["knownu"][l]**2+surfaces[k]["data"]["knownv"][l]**2))
    #for l in range(len(surfaces[k]["data"]["u"])):
        #vels.append(np.sqrt(surfaces[k]["data"]["u"][l]**2+surfaces[k]["data"]["v"][l]**2))
    #eccovelsurfaces.append(np.nanmean(eccovels))
    #velsurfaces.append(np.nanmean(vels))
    #surfs.append(k)
    #print("mean: ",np.nanmean(vels))
    #print("max: ",np.nanmax(vels))
    #print("min: ",np.nanmin(vels))
    #print("std: ",np.nanstd(vels))
#plt.plot(surfs,velsurfaces,label="GEOMEAN")
#plt.plot(surfs,eccovelsurfaces,label="ECCOMEAN")
#plt.legend()
#plt.show()
graph.graphVectorField(inv,"uabs","vabs","z")
#graph.saveAllQuants(inv,"refpics/surfaces/eccoallquants/")
#graph.graphVectorField(inv,"uref","vref",savepath="refpics/vectorfields/RefEcco/",show=False)
#graph.graphVectorField(inv,"u","v","z",show=False,savepath="refpics/vectorfields/nativegeostroph/")
#nstools.surfaceDiagnostic(surfaces)
#graph.graphSurfaces(surfaces,"psiref")
