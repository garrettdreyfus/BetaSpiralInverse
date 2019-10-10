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
import nepbctdextract

#nepbctdextract.nepbCTDExtract("data/newnepbdata.mat","data/nepbctdprofiles.pickle")
#with open('data/nepbctdprofiles.pickle', 'rb') as outfile:
    #profiles = pickle.load(outfile)

######deepestindex=nstools.deepestProfile(profiles)
#######54.062°N, 157.477°W,

#profilechoice = nstools.profileInBox(profiles,-157.5,-157.4,54,54.1,4000)
#print(profilechoice[0].eyed)
#print(profilechoice[1].eyed)
##profilechoice = nstools.getProfileById(profiles,178)
#profilechoice = profilechoice[0]

#surfaces = nstools.runPeerSearch(profiles,100,4000,200,profilechoice,10**10)


#with open('data/annotatednepbprofilessingleref.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,profiles],outfile)

#with open('data/annotatednepbprofilessingleref.pickle', 'rb') as outfile:
    #surfaces,profiles = pickle.load(outfile)

############graph.tsNeutralExplore(profiles)
############graph.plotProfiles(profiles,"first ns search",specialprofile=profilechoice,region="nepb")
##print(len(surfaces[100]["ids"]))
##for k in surfaces[100]["data"]:
    ##print(len(surfaces[100]["data"][k]))
#surfaces = nstools.addDataToSurfaces(profiles,surfaces,2)

#with open('data/nepbsurfaceswithdata.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,profiles],outfile)
with open('data/nepbsurfaceswithdata.pickle', 'rb') as outfile:
    surfaces,profiles = pickle.load(outfile)
#graph.graphSurfaces(surfaces,"nsdiff",region="nepbmerc")

s = nepbctdextract.nepbCTDExtractPointSurfaces("data/newnepbdata.mat")
for k in s.keys():
    s[k]["data"]["psi"] = np.asarray(s[k]["data"]["psi"])*0.5
print(np.nanmax(s[100]["data"]["psi"])-np.nanmin(s[100]["data"]["psi"]))
print(np.nanmax(surfaces[100]["data"]["psi"])-np.nanmin(surfaces[100]["data"]["psi"]))
sminus = nstools.surfaceSubtract(nstools.normalizePSI(s),nstools.normalizePSI(surfaces),method="id")
#sminus = nstools.artificialPSIRef(sminus,reflevel=100)
graph.graphSurfaces(sminus,"psi",region="nepbmerc")
#graph.graphSurfaces(s,"psi",region="nepbmerc")


#nstools.artificialPSIRef(surfaces,reflevel=1700)
#graph.graphSurfaces(surfaces,"psiref",region="nepbmerc")
#graph.graphSurfaces(surfaces,"psiref",show=False,\
        #region="nepbmerc",savepath="refpics/surfaces/psiref1700/")
############################################
#graph.saveAllQuants(surfaces,"refpics/surfaces/nepballquantsnewref/",region="nepb")
#graph.graphSurfaces(surfaces,"pv",region="nepbmerc")
surfaces,neighbors,distances = interptools.interpolateSurfaces(surfaces,fixedgrid="nepb")
#graph.saveAllQuants(surfaces,"refpics/surfaces/nepballquantsnewrefinterp/",region="nepb")
#graph.graphSurfaces(surfaces,"pres",region="nepb")
#with open('data/readytoaddparamsnepb.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,neighbors,distances], outfile)
#with open('data/readytoaddparamsnepb.pickle', 'rb') as outfile:
    #surfaces,neighbors,distances = pickle.load(outfile)
surfaces = nstools.addParametersToSurfaces(surfaces,neighbors,distances)

##graph.saveAllQuants(surfaces,"refpics/surfaces/nepballquants/",region="nepb")
with open('data/ready4inversenepb.pickle', 'wb') as outfile:
    pickle.dump([surfaces,neighbors,distances], outfile)

with open('data/ready4inversenepb.pickle', 'rb') as outfile:
    surfaces,neighbors,distances = pickle.load(outfile)

############################################
#graph.saveAllQuants(surfaces,"refpics/surfaces/nepballquantsnewref/",region="nepb")
#graph.graphSurfaces(surfaces,"pres",region="nepb")
params = {"reflevel":1700,"upperbound":100,"lowerbound":3700,"mixs":{"kvo":False,"kvb":False,"kh":False}}
out= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)

##sensitivity.conditionErrorRefLevel("coupled",surfaces,neighbors,distances)
inv = nstools.streamFuncToUV(out["surfaces"],neighbors,distances)

###graph.graphSurfaces(inv,"e")
graph.graphVectorField(inv,"u","v","pv",metadata=out["metadata"],region="nepbmerc")



##graph.graphSurfaces(surfaces,"pres",region="nepb")
##graph.graphSurfaces(surfaces,"t",region="nepb")

