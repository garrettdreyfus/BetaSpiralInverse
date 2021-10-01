from regionlib import brasil
import os, gsw
import graph
import nstools
import pickle
import interptools
import inverttools
from functools import partial
import parametertools as ptools
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pdb
from scipy.io import savemat
import sensitivity
from officialgamma import addOfficialGamma

#profilesWoce = brasil.extractWoceProfiles(os.path.abspath("data/brasilnc/"))

#profiles = profilesWoce
#with open('data/run0/argoandwoce.pickle', 'wb') as outfile:
    #pickle.dump(profiles,outfile)
with open('data/run0/argoandwoce.pickle', 'rb') as outfile:
   profiles = pickle.load(outfile)

seta = {0}
for k in profiles:
   seta.add(k.cruise)
print(seta)
## profilechoice = nstools.profileInBox(profiles,-40,-20,-31,-28,5000)
## #profilechoice = nstools.profileInBox(profiles,-45,-20,-40,-31,5000)
#profilechoice = nstools.profileInBox(profiles,-25,-20,-31,-28,4500)
#profilechoice = profilechoice[0]

#preinterpsurfaces = nstools.runPeerSearch(profiles,range(200,5000,200),profilechoice,False,10**10)
## preinterpsurfaces = nstools.runPeerSearch(profiles,range(200,5000,200),profilechoice,False,10**4)

#with open('data/run0/annotatedbrasilargowoce.pickle', 'wb') as outfile:
    #pickle.dump([preinterpsurfaces,profiles],outfile)

#with open('data/run0/annotatedbrasilargowoce.pickle', 'rb') as outfile:
     #preinterpsurfaces,profiles = pickle.load(outfile)

#preinterpsurfaces = nstools.addDataToSurfaces(brasil,profiles,preinterpsurfaces)

#with open('data/run0/brasilsurfaceswdata.pickle', 'wb') as infile:
    #pickle.dump([preinterpsurfaces,profiles],infile)
#with open('data/run0/brasilsurfaceswdata.pickle', 'rb') as outfile:
    #preinterpsurfaces,profiles = pickle.load(outfile)

#surfaces,neighbors,distances = \
     #interptools.interpolateSurfaces(brasil,preinterpsurfaces,\
                                     #interpmethod="gam",smart=False,coord="latlon")

#with open('data/run0/interpedbrasil.pickle', 'wb') as outfile:
     #pickle.dump([surfaces,neighbors,distances], outfile)
#with open('data/run0/interpedbrasil.pickle', 'rb') as outfile:
    #[surfaces,neighbors,distances] = pickle.load(outfile)

#surfaces, neighbors, distances = nstools.addParametersToSurfaces(brasil,surfaces,neighbors,distances,H_0=1000)

#with open('data/run0/withparams.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,neighbors,distances], outfile)
#with open('data/run0/withparams.pickle', 'rb') as outfile:
    #[surfaces,neighbors,distances] = pickle.load(outfile)

#params = {"reflevel":2000,"upperbound":1000,"lowerbound":4000,\
        #"mixs":{"kvo":True,"kvb":True,"kh":True},"debug":False,\
            #"3point":True,"edgeguard":True,"H_0":1000
          #}

#out= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)

#with open('data/run0/invertedbrasil.pickle', 'wb') as outfile:
    #pickle.dump([out,neighbors,distances], outfile)
with open('data/run0/invertedbrasil.pickle', 'rb') as outfile:
    [out,neighbors,distances] = pickle.load(outfile)

surfaces = addOfficialGamma(out["surfaces"])
graph.graphSurfaces(brasil,surfaces,"gamma",select=[2400])

