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
#with open('data/nepbsurfaceswithdata.pickle', 'rb') as outfile:
  #preinterpsurfaces,profiles = pickle.load(outfile)

#floats = nepbctdextract.floatExtract("data/Run0.new.mat")
#floats = nstools.artificialPSIRef(floats)
###graph.graphSurfaces(floats,"psi",region="nepbmerc")
#preinterps = nepbctdextract.extractPointSurfaces("data/newnepbdata.mat")
#preinterps = nstools.artificialPSIRef(preinterps)
#fulls = nstools.surfaceConcat(floats,preinterps)

#with open('data/fulls.pickle', 'wb') as outfile:
  #pickle.dump(fulls, outfile)
#with open('data/fulls.pickle', 'rb') as outfile:
  #fulls = pickle.load(outfile)

s = nepbctdextract.nepbCTDExtractInterpSurfaces("data/Run0.new.mat",calcDeriv=True)

s = nstools.removeAleutians(s)
#surfaces = nstools.depthCopy(ref= s,surfaces={})

#surfaces,neighbors,distances = interptools.interpolateSurfaces(surfaces,\
        #fixedgrid="hautala",gaminterpolate=False)

###nstools.surfaceDiagnostic(surfaces)

#with open('data/readytoaddparamsnepb.pickle', 'wb') as outfile:
  #pickle.dump([surfaces,neighbors,distances], outfile)
#with open('data/readytoaddparamsnepb.pickle', 'rb') as outfile:
  #surfaces,neighbors,distances = pickle.load(outfile)

#surfaces = nstools.addParametersToSurfaces(surfaces,\
        #neighbors,distances,["s","t"])

##nstools.inverseReady(surfaces)

#with open('data/ready4inversenepb.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,neighbors,distances], outfile)

with open('data/ready4inversenepb.pickle', 'rb') as outfile:
    surfaces,neighbors,distances = pickle.load(outfile)

#graph.graphSurfaces(surfaces,"d2sdx2",region="nepbmerc")
sminus = nstools.surfaceSubtract(surfaces,s,metric="/")
with open('data/sminus.pickle', 'wb') as outfile:
    pickle.dump(sminus, outfile)

with open('data/sminus.pickle', 'rb') as outfile:
    sminus = pickle.load(outfile)

#print(sminus[700]["data"].keys())
graph.graphSurfaces(sminus,"d2sdy2",region="nepbmerc")
#graph.graphSurfaces(sminus,"s",region="nepb")

#graph.graphVectorField(surfaces,"u","v","psi",region="nepb",transform=False)
#graph.graphVectorField(surfaces,"u","v","psi",region="nepbmerc")
#nstools.inverseReady(surfaces)

##s = nepbctdextract.nepbCTDExtractInterpSurfaces("data/Run0.new.mat")
##############################################
graph.saveAllQuants(sminus,"refpics/surfaces/sminushgrid/",
        region="nepb",select = [200,1000])
###graph.graphSurfaces(surfaces,"pres",region="nepb")

#params = {"reflevel":1700,"upperbound":1100,"lowerbound":3500,"mixs":{"kvo":True,"kvb":True,"kh":True}}

###nstools.surfaceDiagnostic(surfaces)
#out= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)

#with open('data/inverseout.pickle', 'wb') as outfile:
    #pickle.dump([out,neighbors,distances], outfile)
#with open('data/inverseout.pickle', 'rb') as outfile:
    #[out,neighbors,distances] = pickle.load(outfile)
#######sensitivity.conditionErrorRefLevel("coupled",surfaces,neighbors,distances)
####print(out.keys())
#inv = nstools.streamFuncToUV(out["surfaces"],neighbors,distances)
##print(out["metadata"])
########graph.graphSurfaces(inv,"e")
#graph.graphVectorField(inv,"uabs","vabs","psi",\
        #metadata=out["metadata"],region="nepbmerc",savepath="refpics/vectorfields/nomixnepb/pv/",show=False)
#graph.graphVectorField(inv,"uabs","vabs","psi",\
        #metadata=out["metadata"],region="nepbmerc",savepath="refpics/vectorfields/nomixnepb/psi/",show=False)
#graph.graphVectorField(inv,"uabs","vabs","pv",metadata=out["metadata"],\
        #region="nepbmerc",savepath="refpics/vectorfields/nolateralnepb/",show=False)



##graph.graphSurfaces(surfaces,"pres",region="nepb")
##graph.graphSurfaces(surfaces,"t",region="nepb")

