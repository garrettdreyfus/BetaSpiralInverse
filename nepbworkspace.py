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
#with open('data/nepbctdprofiles.pickle', 'rb') as outfile:
    #profiles = pickle.load(outfile)

###############54.062°N, 157.477°W,

#profilechoice = nstools.profileInBox(profiles,-157.5,-157.4,54,54.1,4000)
#profilechoice = profilechoice[0]

###surfaces = nstools.runPeerSearch(profiles,range(100,4000,200),profilechoice,10**10)
#ns = nepbctdextract.neutralSurfaces("data/newnepbdata.mat")
##print(ns)
#surfaces = nstools.runPeerSearch(profiles,ns,profilechoice,10**10)



########surfaces = nstools.depthCopy(ref=s)

#with open('data/annotatednepbprofilessingleref.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,profiles],outfile)
#with open('data/annotatednepbprofilessingleref.pickle', 'rb') as outfile:
    #preinterpsurfaces,profiles = pickle.load(outfile)

#surfaces = nstools.addDataToSurfaces(profiles,preinterpsurfaces,2)

#with open('data/nepbsurfaceswithdata.pickle', 'wb') as outfile:
  #pickle.dump([surfaces,profiles],outfile)
#with open('data/nepbsurfaceswithdata.pickle', 'rb') as outfile:
  #preinterpsurfaces,profiles = pickle.load(outfile)

###floats = nepbctdextract.floatExtract("data/Run0.new.mat")
###floats = nstools.artificialPSIRef(floats)
###preinterps = nepbctdextract.extractPointSurfaces("data/newnepbdata.mat")
###preinterps = nstools.artificialPSIRef(preinterps)
###fulls = nstools.surfaceConcat(floats,preinterps)

###with open('data/fulls.pickle', 'wb') as outfile:
  ###pickle.dump(fulls, outfile)
###with open('data/fulls.pickle', 'rb') as outfile:
  ###fulls = pickle.load(outfile)
####scompare = nepbctdextract.nepbCTDExtractInterpSurfaces("data/Run0.new.mat",calcDeriv=False)
####s = nepbctdextract.nepbCTDExtractInterpSurfaces("data/Run0.new.mat",calcDeriv=True)
###fulls = nstools.removeAleutians(fulls)
####sminus = nstools.surfaceSubtract(s,scompare,metric="/")

####s = nstools.removeAleutians(s)
####surfaces = nstools.depthCopy(ref= s,surfaces={})

#######graph.graphSurfaces(sminus,"dqdy",region="nepbmerc",select=[400,600])
#######graph.graphSurfaces(sminus,"dqdx",region="nepbmerc",select=[400,600])
#######graph.graphSurfaces(sminus,"dsdy",region="nepbmerc",select=[400,600])
#######graph.graphSurfaces(sminus,"dsdx",region="nepbmerc",select=[400,600])
#######graph.graphSurfaces(sminus,"dtdy",region="nepbmerc",select=[400,600])
#######graph.graphSurfaces(sminus,"dtdx",region="nepbmerc",select=[400,600])
#surfaces,neighbors,distances = interptools.interpolateSurfaces(preinterpsurfaces,\
        #fixedgrid="nepb",gaminterpolate=True)
##graph.graphSurfaces(surfaces,"pres",region="nepb")

#########nstools.surfaceDiagnostic(surfaces)

#with open('data/readytoaddparamsnepb.pickle', 'wb') as outfile:
  #pickle.dump([surfaces,neighbors,distances], outfile)
with open('data/readytoaddparamsnepb.pickle', 'rb') as outfile:
  surfaces,neighbors,distances = pickle.load(outfile)

surfaces = nstools.addParametersToSurfaces(surfaces,\
        neighbors,distances,[])

#graph.graphSurfaces(surfaces,"psi",region="nepb",select=[1500,2500,3600])
#nstools.inverseReady(surfaces)

with open('data/ready4inversenepb.pickle', 'wb') as outfile:
    pickle.dump([surfaces,neighbors,distances], outfile)

with open('data/ready4inversenepb.pickle', 'rb') as outfile:
    surfaces,neighbors,distances = pickle.load(outfile)
nstools.surfaceDiagnostic(surfaces)

#scompare = nepbctdextract.nepbCTDExtractInterpSurfaces("data/Run0.new.mat",calcDeriv=False)
graph.graphSurfaces(surfaces,"kvbdiagnostic",region="nepb",savepath="refpics/surfaces/bathvardepthmean/",show=False)
###scompare = nstools.domainChop(scompare)
###surfaces = nstools.domainChop(surfaces)
##sminus = nstools.surfaceSubtract(surfaces,scompare,metric="/")
##graph.graphSurfaces(sminus,"dqdy",region="nepb",select=[400,600])
##graph.graphSurfaces(sminus,"dqdx",region="nepb",select=[400,600])
##graph.graphSurfaces(sminus,"dqnotdy",region="nepb",select=[400,600])
##graph.graphSurfaces(sminus,"dqnotdx",region="nepb",select=[400,600])
##graph.graphSurfaces(sminus,"dsdx",region="nepb",select=[400,600])
##graph.graphSurfaces(sminus,"dsdy",region="nepb",select=[400,600])
###surfaces = artificialPSIRef(surfaces):
##graph.graphSurfaces(sminus,"dqnotdy",region="nepb",select=[400,600])
###graph.graphSurfaces(scompare,"bathvar",region="nepb",select=[400,600])
##surfaces = nstools.domainChop(surfaces)
####with open('data/sminus.pickle', 'wb') as outfile:
    ####pickle.dump(sminus, outfile)

###with open('data/sminus.pickle', 'rb') as outfile:
    ###sminus = pickle.load(outfile)

###graph.graphVectorField(surfaces,"u","v","psi",region="nepb",transform=False)
###graph.graphVectorField(surfaces,"u","v","psi",region="nepbmerc")
###nstools.inverseReady(surfaces)

####s = nepbctdextract.nepbCTDExtractInterpSurfaces("data/Run0.new.mat")
################################################
###graph.saveAllQuants(sminus,"refpics/surfaces/sminushgrid/",
        ###region="nepb",select = [200,1000])
#####graph.graphSurfaces(surfaces,"pres",region="nepb")

#params = {"reflevel":1700,"upperbound":1100,"lowerbound":3500,"mixs":{"kvo":False,"kvb":False,"kh":False},"3point":True,"edgeguard":True}

######nstools.surfaceDiagnostic(surfaces)
#out= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)

#with open('data/inverseout.pickle', 'wb') as outfile:
    #pickle.dump([out,neighbors,distances], outfile)
#with open('data/inverseout.pickle', 'rb') as outfile:
    #[out,neighbors,distances] = pickle.load(outfile)
###print(out["metadata"])
#inv = nstools.streamFuncToUV(out["surfaces"],neighbors,distances)

###graph.graphVectorField(inv,"uabs","vabs","pv",\
        ###metadata=out["metadata"],region="nepb",\
        ###transform=False,select=[1500,2500,3600])

#####graph.graphSurfaces(inv,"psinew",region="nepb",select=[1400,np.inf])
#####graph.graphSurfaces(inv,"psi",region="nepb",select=[1400,np.inf])
#####graph.velocityHeatMap(inv,"uabs",167)
####graph.graphVectorField(inv,"uabs","vabs","psi",\
        ####metadata=out["metadata"],region="nepb",\
        ####transform=False,select=[1400,np.inf])
#graph.graphVectorField(inv,"uabs","vabs","psinew",\
        #metadata=out["metadata"],region="nepb",\
        #savepath="refpics/vectorfields/allmixpyneutralsearch8splinexymeanfix/psinew/",show=False,transform=True)

