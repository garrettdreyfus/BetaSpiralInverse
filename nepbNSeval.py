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




#nepbctdextract.nepbCTDExtract("data/newnepbdata.mat","data/nepbctdprofiles.pickle")
with open('data/nepbctdprofiles.pickle', 'rb') as outfile:
    profiles = pickle.load(outfile)

################54.062°N, 157.477°W,

profilechoice = nstools.profileInBox(profiles,-157.5,-157.4,54,54.1,4000)
profilechoice = profilechoice[0]
#print(profilechoice)

gammavals = [25.875,26.51,26.862,27.158,27.3605,27.526,27.6575,27.7825,27.8700, \
27.9275,27.965,27.99,28.015,28.03,28.0475,28.0625,28.08,28.108,28.136,28.164,28.2,28.33,28.36]
#surfaces = nstools.runPeerSearch(profiles,range(100,4000,200),profilechoice,True,gammas=gammavals)
surfaces = nstools.runPeerSearch(profiles,range(100,4000,200),profilechoice,True,10**5)

surfaces = nstools.addDataToSurfaces(nepb,profiles,surfaces)


surfaces,neighbors,distances = interptools.interpolateSurfaces(nepb,surfaces,interpmethod="gam",smart=False,coord="latlon")


surfaces = nstools.addParametersToSurfaces(nepb,surfaces,\
        neighbors,distances)

surfaces = nstools.neutralityError(surfaces)

graph.graphSurfaces(nepb,surfaces,"nserror",show=False, savepath="../arcticcirc-pics/surfaces/nepbnserror/gamma/")


