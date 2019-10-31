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
import pdb

#nepbctdextract.nepbCTDExtract("data/newnepbdata.mat","data/nepbctdprofiles.pickle")
with open('data/nepbctdprofiles.pickle', 'rb') as outfile:
    profiles = pickle.load(outfile)

##############54.062°N, 157.477°W,

#profilechoice = nstools.profileInBox(profiles,-157.5,-157.4,54,54.1,4000)
#profilechoice = profilechoice[0]

##surfaces = nstools.runPeerSearch(profiles,range(100,4000,200),profilechoice,10**10)
#ns = nepbctdextract.neutralSurfaces("data/newnepbdata.mat")
#print(ns)
#surfaces = nstools.runPeerSearch(profiles,ns,profilechoice,10**10)



#######surfaces = nstools.depthCopy(ref=s)

#####with open('data/annotatednepbprofilessingleref.pickle', 'wb') as outfile:
    ######pickle.dump([surfaces,profiles],outfile)
#with open('data/annotatednepbprofilessingleref.pickle', 'rb') as outfile:
    #preinterpsurfaces,profiles = pickle.load(outfile)

#surfaces = nstools.addDataToSurfaces(profiles,preinterpsurfaces,2)

#with open('data/nepbsurfaceswithdata.pickle', 'wb') as outfile:
  #pickle.dump([surfaces,profiles],outfile)
with open('data/nepbsurfaceswithdata.pickle', 'rb') as outfile:
  preinterpsurfaces,profiles = pickle.load(outfile)

#surfaces,neighbors,distances = interptools.interpolateSurfaces(surfaces,fixedgrid="nepb")

#with open('data/readytoaddparamsnepb.pickle', 'wb') as outfile:
  #pickle.dump([surfaces,neighbors,distances], outfile)
#with open('data/readytoaddparamsnepb.pickle', 'rb') as outfile:
  #surfaces,neighbors,distances = pickle.load(outfile)

#surfaces = nstools.addParametersToSurfaces(surfaces,neighbors,distances)

##############graph.saveAllQuants(surfaces,"refpics/surfaces/nepballquants/",region="nepb")
#with open('data/ready4inversenepb.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,neighbors,distances], outfile)

with open('data/ready4inversenepb.pickle', 'rb') as outfile:
    surfaces,neighbors,distances = pickle.load(outfile)

###print(nstools.surfaceDiagnostic(s))
###print(nstools.surfaceDiagnostic(surfaces))
preinterps = nepbctdextract.nepbCTDExtractPointSurfaces("data/newnepbdata.mat")

s = nepbctdextract.nepbCTDExtractInterpSurfaces("data/Run0.new.mat")
s = nstools.adddSdP(s)
surfaces = nstools.adddSdP(surfaces)
sminus = nstools.surfaceSubtract(preinterpsurfaces,preinterps,metric="%")
#graph.graphSurfaces(sminus,"s",region="nepbmerc")


#graph.saveAllQuants(surfaces,"refpics/surfaces/interpcompare/python/",region="nepb")
#graph.saveAllQuants(s,"refpics/surfaces/interpcompare/matlab/",region="nepb")

sminus = nstools.surfaceSubtract(surfaces,s,metric="%")
#graph.graphSurfaces(sminus,"dsdz",region="nepbmerc")
#graph.graphSurfaces(s,"t",region="nepbmerc")
#surfaces = nstools.artificialPSIRef(surfaces,reflevel=1700)
graph.graphSurfaces(sminus,"d2thetads2",region="nepbmerc")
graph.graphSurfacesOneContour(surfaces,preinterpsurfaces,"t",region="nepb",savepath="refpics/surfaces/gamcomp/python/",show=False)
graph.graphSurfacesOneContour(s,preinterps,"t",region="nepb",savepath="refpics/surfaces/gamcomp/matlab/",show=False)
def filterThisArray(surf):
    return np.asarray(surf["data"]["pres"])[np.logical_and.reduce((np.asarray(surf["lons"])<-100,np.asarray(surf["lons"])>-150,np.asarray(surf["lats"])>20,100>np.asarray(surf["lats"])))]

means = {"py":[],"interppy":[],"matlab":[],"interpmatlab":[],"depths":[]}
for k in surfaces.keys():
    if k in s.keys():
        means["py"].append(np.nanmean(filterThisArray(preinterpsurfaces[k]))-k)
        means["interppy"].append(np.nanmean(filterThisArray(surfaces[k]))-k)
        means["matlab"].append(np.nanmean(filterThisArray(preinterps[k]))-k)
        means["interpmatlab"].append(np.nanmean(filterThisArray(s[k]))-k)
        means["depths"].append(k)
plt.plot(means["depths"],means["py"],label="py")
plt.plot(means["depths"],means["interppy"],label="interppy")
plt.plot(means["depths"],means["interpmatlab"],label="interpmatlab")
plt.plot(means["depths"],means["matlab"],label="matlab")
plt.legend()
plt.xlabel("Neutral Surface")
plt.ylabel("Deviation from labelled neutral surface depth")
plt.show()

#s = nstools.artificialPSIRef(s,reflevel=1700)
#sminus = nstools.surfaceSubtract(surfaces,s,metric="%")
#graph.graphSurfaces(sminus,"dsdz",region="nepbmerc")
##print(nstools.surfaceDiagnostic(s))
#graph.graphSurfaces(sminus,"dsdz",region="nepbmerc")


#metricname = {"-":"difference","%":"percentdiff","/":"division"}
#for i in metricname.keys():
    #sminus = nstools.surfaceSubtract(surfaces,s,metric=i)
    #graph.saveAllQuants(sminus,"refpics/surfaces/nepballquantsmetrics/"+metricname[i]+"/",region="nepb")
#graph.saveAllQuants(sminus,"refpics/surfaces/interpercenterrors/",region="nepb")

#############################################
##graph.saveAllQuants(surfaces,"refpics/surfaces/nepballquantsnewref/",region="nepb")
##graph.graphSurfaces(surfaces,"pres",region="nepb")
#params = {"reflevel":1700,"upperbound":100,"lowerbound":3700,"mixs":{"kvo":True,"kvb":True,"kh":True}}
#out= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)

#with open('data/inverseout.pickle', 'wb') as outfile:
    #pickle.dump([out,neighbors,distances], outfile)
#with open('data/inverseout.pickle', 'rb') as outfile:
    #[out,neighbors,distances] = pickle.load(outfile)
#####sensitivity.conditionErrorRefLevel("coupled",surfaces,neighbors,distances)
##print(out.keys())
#inv = nstools.streamFuncToUV(out["surfaces"],neighbors,distances)

######graph.graphSurfaces(inv,"e")
#graph.graphVectorField(inv,"uabs","vabs","pv",metadata=out["metadata"],region="nepbmerc")
#graph.graphVectorField(inv,"uabs","vabs","pv",metadata=out["metadata"],\
        #region="nepbmerc",savepath="refpics/vectorfields/nolateralnepb/",show=False)



##graph.graphSurfaces(surfaces,"pres",region="nepb")
##graph.graphSurfaces(surfaces,"t",region="nepb")

