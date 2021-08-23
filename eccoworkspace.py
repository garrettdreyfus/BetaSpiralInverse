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
    
with open('data/eccoprofiles.pickle', 'rb') as outfile:
    profiles = pickle.load(outfile)

##graph.plotProfiles(profiles,"yo",data="t",depth="200")

#deepestindex=nstools.deepestProfile(profiles)
##graph.plotProfiles(profiles,"ecco",specialprofile=profiles[deepestindex])

##profiles = nstools.filterCruises(profiles,offsets.keys())
##profiles = saloffset.applyOffsets(profiles,offsets)

##profilechoice = random.choice(nstools.profileInBox(profiles,-180,180,85,90))
#profilechoice = profiles[deepestindex]

#surfaces = nstools.runPeerSearch(profiles,200,4000,200,profilechoice,1000)

#with open('data/eccosurfaces.pickle', 'wb') as outfile:
  #pickle.dump(surfaces, outfile)

##with open('data/eccoannotatedprofiles.pickle', 'wb') as outfile:
    ##pickle.dump(profiles, outfile)
##with open('data/eccoannotatedprofiles.pickle', 'rb') as outfile:
    ##profiles=pickle.load(outfile)

with open('data/eccosurfaces.pickle', 'rb') as outfile:
    surfaces=pickle.load(outfile)
#######graph.graphSurfaces(surfaces,"pres")
#########graph.tsNeutralExplore(profiles)

preinterpsurfaces = nstools.addDataToSurfaces(profiles,surfaces,2)
##with open('data/eccosurfwithd.pickle', 'wb') as outfile:
    ##pickle.dump(surfaces, outfile)

##with open('data/eccosurfwithd.pickle', 'rb') as outfile:
    ##surfaces=pickle.load(outfile)

###nstools.surfaceDiagnostic(surfaces)
##with open('data/postinterp.pickle', 'wb') as outfile:
    ##pickle.dump([surfaces,neighbors,distances], outfile)
##with open('data/postinterp.pickle', 'rb') as outfile:
    ##[surfaces,neighbors,distances]=pickle.load(outfile)


nserrors = {}
for l in range(4,20):
    surfaces,neighbors,distances = interptools.interpolateSurfaces(brasil,preinterpsurfaces,\
                                                                   interpmethod="gam",smart=False,coord="latlon",splines=l)
    surfaces = nstools.addParametersToSurfaces(brasil,surfaces,\
        neighbors,distances)
    surfaces = nstools.neutralityError(surfaces)
    error = []
    for k in surfaces.keys():
        error += list(surfaces[k]["data"]["nserror"])
    nserrors[l] = np.asarray(error).flatten()
for k in nserrors.keys():
    sns.distplot(np.log10(nserrors[k]),kde_kws={"lw": 2, "label": str(k)})


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
