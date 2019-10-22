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

########deepestindex=nstools.deepestProfile(profiles)
#########54.062°N, 157.477°W,

##profilechoice = nstools.profileInBox(profiles,-157.5,-157.4,54,54.1,4000)
###print(profilechoice[0].eyed)
###print(profilechoice[1].eyed)
####profilechoice = nstools.getProfileById(profiles,178)
##profilechoice = profilechoice[0]

##surfaces = nstools.runPeerSearch(profiles,100,4000,200,profilechoice,10**10)

#s = nepbctdextract.nepbCTDExtractPointSurfaces("data/newnepbdata.mat")

##surfaces = nstools.depthCopy(ref=s)

##with open('data/depthcopy.pickle', 'wb') as outfile:
    ##pickle.dump([surfaces,profiles],outfile)

#with open('data/depthcopy.pickle', 'rb') as outfile:
    #bsurfaces,profiles = pickle.load(outfile)

##pdb.set_trace()
##############graph.tsNeutralExplore(profiles)
##############graph.plotProfiles(profiles,"first ns search",specialprofile=profilechoice,region="nepb")
####print(len(surfaces[100]["ids"]))
####for k in surfaces[100]["data"]:
    ####print(len(surfaces[100]["data"][k]))
#surfaces = nstools.addDataToSurfaces(profiles,bsurfaces,2)
#nstools.surfaceDiagnostic(bsurfaces)

#with open('data/nepbsurfaceswithdata.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,profiles],outfile)
#with open('data/nepbsurfaceswithdata.pickle', 'rb') as outfile:
    #surfaces,profiles = pickle.load(outfile)
#s = nepbctdextract.nepbCTDExtractPointSurfaces("data/newnepbdata.mat")
#sminus = nstools.surfaceSubtract(s,surfaces)
#graph.graphSurfaces(sminus,"psi",region="nepbmerc")
######sminus = nstools.artificialPSIRef(sminus,reflevel=100)
######graph.graphSurfaces(sminus,"psi",region="nepbmerc",savepath="refpics/psipercenterror/",show=False)
######graph.saveAllQuants(sminus,"refpics/surfaces/percenterrors/",region="nepb")
######graph.graphSurfaces(s,"psi",region="nepbmerc")


######nstools.artificialPSIRef(surfaces,reflevel=1700)
######graph.graphSurfaces(surfaces,"psiref",region="nepbmerc")
######graph.graphSurfaces(surfaces,"psiref",show=False,\
        ######region="nepbmerc",savepath="refpics/surfaces/psiref1700/")
#################################################
######graph.saveAllQuants(surfaces,"refpics/surfaces/nepballquantsnewref/",region="nepb")
######graph.graphSurfaces(surfaces,"pv",region="nepbmerc")
#surfaces,neighbors,distances = interptools.interpolateSurfaces(surfaces,fixedgrid="nepb")
####graph.graphSurfaces(surfaces,"t",region="nepbmerc")
#####graph.saveAllQuants(surfaces,"refpics/surfaces/nepballquantsnewrefinterp/",region="nepb")
#####graph.graphSurfaces(surfaces,"pres",region="nepb")
#with open('data/readytoaddparamsnepb.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,neighbors,distances], outfile)
#with open('data/readytoaddparamsnepb.pickle', 'rb') as outfile:
    #surfaces,neighbors,distances = pickle.load(outfile)
##print(surfaces.keys())
#surfaces = nstools.addParametersToSurfaces(surfaces,neighbors,distances)

######graph.saveAllQuants(surfaces,"refpics/surfaces/nepballquants/",region="nepb")
#with open('data/ready4inversenepb.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,neighbors,distances], outfile)

with open('data/ready4inversenepb.pickle', 'rb') as outfile:
    surfaces,neighbors,distances = pickle.load(outfile)

##print(nstools.surfaceDiagnostic(s))
##print(nstools.surfaceDiagnostic(surfaces))
#surfaces = nstools.artificialPSIRef(surfaces)
#s1 = nepbctdextract.nepbCTDExtractPointSurfaces("data/newnepbdata.mat")
#s2 = nepbctdextract.nepbCTDExtractInterpSurfaces("data/Run0.new.mat")
#sminus = nstools.surfaceSubtract(s1,s2,metric="%",offset=8.8)
#sminus = nstools.surfaceSubtract(surfaces,s,metric="%")
#graph.graphSurfaces(sminus,"psi",region="nepbmerc")

#graph.graphSurfaces(sminus,"dsdz",region="nepbmerc")
#graph.graphSurfaces(sminus,"pv",region="nepbmerc")
#metricname = {"-":"difference","%":"percentdiff","/":"division"}
#for i in metricname.keys():
    #sminus = nstools.surfaceSubtract(surfaces,s,metric=i)
    #graph.saveAllQuants(sminus,"refpics/surfaces/nepballquantsmetrics/"+metricname[i]+"/",region="nepb")
#graph.saveAllQuants(sminus,"refpics/surfaces/interpercenterrors/",region="nepb")

#############################################
##graph.saveAllQuants(surfaces,"refpics/surfaces/nepballquantsnewref/",region="nepb")
##graph.graphSurfaces(surfaces,"pres",region="nepb")
params = {"reflevel":1700,"upperbound":100,"lowerbound":3700,"mixs":{"kvo":True,"kvb":True,"kh":True}}
out= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)

with open('data/inverseout.pickle', 'wb') as outfile:
    pickle.dump([out,neighbors,distances], outfile)
with open('data/inverseout.pickle', 'rb') as outfile:
    [out,neighbors,distances] = pickle.load(outfile)
####sensitivity.conditionErrorRefLevel("coupled",surfaces,neighbors,distances)
#print(out.keys())
inv = nstools.streamFuncToUV(out["surfaces"],neighbors,distances)

#####graph.graphSurfaces(inv,"e")
graph.graphVectorField(inv,"uabs","vabs","pv",metadata=out["metadata"],region="nepbmerc")
#graph.graphVectorField(inv,"uabs","vabs","pv",metadata=out["metadata"],\
        #region="nepbmerc",savepath="refpics/vectorfields/nolateralnepb/",show=False)



##graph.graphSurfaces(surfaces,"pres",region="nepb")
##graph.graphSurfaces(surfaces,"t",region="nepb")

