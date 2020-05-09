import nstools
import numpy as np
import inverttools
#import eccotools
import interptools
import bathtools
import saloffset
import graph
import matplotlib.pyplot as plt
import parametertools as ptools
import pickle
import sensitivity
import random
import pdb
import scipy.io as sio
from regionlib import brasil, nepb
from progress.bar import Bar




#nepb.nepbCTDExtract("data/newnepbdata.mat","data/nepbctdprofiles.pickle")
# with open('data/ecconepbprofiles.pickle', 'rb') as outfile:
#     profiles = pickle.load(outfile)
# for p in profiles:
#   if p.lon>180: p.lon=p.lon-360
# #################54.062°N, 157.477°W,

# profilechoice = nstools.profileInBox(profiles,-157,-154,53,55,2000)
# profilechoice = profilechoice[0]
# #print(profilechoice)
# ##surfaces = nstools.runPeerSearch(profiles,range(100,4000,200),profilechoice,True,10**3)

# ##gammavals = [25.875,26.51,26.862,27.158,27.3605,27.526,27.6575,27.7825,27.8700, \
#         ##27.9275,27.965,27.99,28.015,28.03,28.0475,28.0625,28.08,28.108,28.136,28.164,28.2,28.33,28.36]
# ##surfaces = nstools.runPeerSearch(profiles,range(100,4000,200),profilechoice,True,gammas=gammavals)

# #####print(range(100,6000,200)[len(gammavals)-1])
# ##preinterpsurfaces = nstools.runPeerSearch(profiles,range(100,4300,200),profilechoice,False,10**4,gammas=gammavals)
# ##ns = nepbctdextract.neutralSurfaces("data/newnepbdata.mat")
# ##print(ns)
# surfaces = nstools.runPeerSearch(profiles,range(100,4300,200),profilechoice,False,10**10)

# with open('data/annotatednepbprofilessingleref.pickle', 'wb') as outfile:
#     pickle.dump([surfaces,profiles],outfile)
with open('data/annotatednepbprofilessingleref.pickle', 'rb') as outfile:
    preinterpsurfaces,profiles = pickle.load(outfile)

surfaces = nstools.addKnownPsiFromModelUV(nepb,preinterpsurfaces,profiles)

with open('data/annotatednepbprofilessingleref.pickle', 'wb') as outfile:
    pickle.dump([surfaces,profiles],outfile)

graph.graphSurfaces(nepb,surfaces,"knownpsi",show=False,savepath="../arcticcirc-pics/refpics/surfaces/modelpsi/",stds=1)

surfaces = nstools.addDataToSurfaces(nepb,profiles,preinterpsurfaces)

#with open('data/nepbsurfaceswithdata.pickle', 'wb') as outfile:
  #pickle.dump([surfaces,profiles],outfile)
#for k in range(5,25):
# with open('data/nepbsurfaceswithdata.pickle', 'rb') as outfile:
#   preinterpsurfaces,profiles = pickle.load(outfile)
# #graph.graphVectorField(nepb,preinterpsurfaces,"knownu","knownv","pres",show=False, \
#         #savepath="refpics/vectorfields/nepbecco/known/",refarrow=0.003,contour=False,scale=0.03,transform=False)


# surfaces,neighbors,distances = interptools.interpolateSurfaces(nepb,preinterpsurfaces,\
#         interpmethod="linear",smart=False,coord="latlon")

# graph.graphSurfaces(nepb,preinterpsurfaces,"psi",
#             contour=True,savepath="../arcticcirc-pics/surfaces/psisolution/",show=False,stds=1)

#graph.graphVectorField(nepb,surfaces,"knownu","knownv","pres",\
        #refarrow=0.003,\
        #transform=False,show=False,scale=0.05,\
        #savepath="../arcticcirc-pics/refpics/vectorfields/nepbecco/linearknown/")

#with open('data/readytoaddparamsnepb.pickle', 'wb') as outfile:
  #pickle.dump([surfaces,neighbors,distances], outfile)
#with open('data/readytoaddparamsnepb.pickle', 'rb') as outfile:
  #surfaces,neighbors,distances = pickle.load(outfile)

###graph.graphSurfaces(surfaces,"pres",region="nepb",select=range(1400,4000))
#surfaces = nstools.addParametersToSurfaces(nepb,surfaces,\
        #neighbors,distances)


###nstools.inverseReady(surfaces)

#with open('data/ready4inversenepb.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,neighbors,distances], outfile)

#with open('data/ready4inversenepb.pickle', 'rb') as outfile:
    #surfaces,neighbors,distances = pickle.load(outfile)

##graph.graphVectorField(nepb,surfaces,"knownu","knownv","z",\
        ##refarrow=0.003,\
        ##transform=False,show=False,scale=0.05,\
        ##savepath="refpics/vectorfields/nepbecco/knowninterp/")
##params = {"reflevel":1700,"upperbound":1100,"lowerbound":3500,\
        ##"mixs":{"kvo":True,"kvb":True,"kh":True},"debug":False,\
        ##"3point":True,"edgeguard":True}

#params = {"reflevel":1700,"upperbound":1100,"lowerbound":3300,\
        #"mixs":{"kvo":False,"kvb":False,"kh":False},"debug":False,\
        #"3point":True,"edgeguard":True,"modelmixing":True,\
        #"scalecoeffs":{"Ar":1,"kvo":1,"kvb":1,"kh":1}}

########nstools.surfaceDiagnostic(surfaces)
#out= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)
#with open('data/inverseout.pickle', 'wb') as outfile:
    #pickle.dump([out,neighbors,distances], outfile)
#with open('data/inverseout.pickle', 'rb') as outfile:
    #[out,neighbors,distances] = pickle.load(outfile)
###print(out["metadata"])
#inv = nstools.streamFuncToUV(out["surfaces"],neighbors,distances)
#inv = ptools.calcFRho(inv)
###print(out["metadata"])
##for lon in np.linspace(-190,-122,34):
    ##graph.northSouthTransect(inv,"uabs",lon=lon,show=False,coordrange=[25,55],\
            ##cbarrange=[-0.002,0.004],savepath="refpics/transects/eccocomp/modelmix/smart/")
##for lat in np.linspace(18,58,20):
    ##graph.northSouthTransect(inv,"uabs",lat=lat,show=False,coordrange=[25,55],\
            ##cbarrange=[-0.002,0.004],savepath="refpics/transects/eccocomp/modelmix/smart/")
###graph.northSouthTransect(inv,"uabs",lon=-145,show=False,coordrange=[25,55],\
        ###cbarrange=[-0.002,0.004],savepath="refpics/transects/eccocomp/modelmix/smart/")

#surfaces,neighbors,distances = interptools.interpolateSurfaces(nepb,inv,\
        #interpmethod="linear",smart=False,coord="latlon")

#graph.graphVectorField(nepb,inv,"uabs","vabs","pres",\
        #metadata=out["metadata"],refarrow=0.003,\
        #transform=False,show=False,scale=0.05,\
        #savepath="../arcticcirc-pics/refpics/vectorfields/nepbecco/inversemodelmixgaminterp/")

#graph.graphVectorField(inv,"uabs","vabs","psinew",\
        #metadata=out["metadata"],region="nepb",\
        #transform=True,show=False,
        #savepath="refpics/vectorfields/50varmesh/psinew/")

