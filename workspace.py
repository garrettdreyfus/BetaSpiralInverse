import nstools
import inverttools
import bathtools
import saloffset
import matplotlib.pyplot as plt
import pickle
from mpl_toolkits.basemap import Basemap
import numpy as np
import json
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pyproj
import graph
import copy
import parametertools as ptools
from pprint import pprint
import pdb
    
#profiles,deepestindex = nstools.extractProfilesBox(["data/1500mprofiles.json"],-180,180,65,90)
#profiles,deepestindex = nstools.removeNorwegianSea(profiles)

#fileObject = open("data/1500NoNorwegian.pickle",'rb')  
#offsets,badfiles,beepestindex = pickle.load(fileObject)

##profiles = nstools.filterCruises(profiles,offsets.keys())
##profiles = saloffset.applyOffsets(profiles,offsets)

#####profilechoice = random.choice(nstools.profileInBox(profiles,-180,180,85,90))
#####profilechoice = nstools.getProfileById(profiles,"286364")

#####surfaces = nstools.runPeerSearch(profiles,deepestindex,200,4000,200,profilechoice,1000)
#####fileObject = open(str(profilechoice.eyed)+"new.pickle",'wb')  
###### load the object from the file into var b
#####b = pickle.dump(surfaces,fileObject)  
#####fileObject.close()

#with open('data/286364new.pickle', 'rb') as outfile:
    #surfaces=pickle.load(outfile)

#surfaces = nstools.addDataToSurfaces(profiles,surfaces,2)
#surfaces = nstools.addStreamFunc(surfaces,profiles)

#surfaces =nstools.addXYToSurfaces(surfaces)
#interpolatedsurfaces,neighbors,lookups = nstools.interpolateSurfaces(surfaces)
#staggeredsurfaces = nstools.fillOutEmptyFields(interpolatedsurfaces)
#staggeredsurfaces = nstools.addHeight(staggeredsurfaces)

###with open('data/ready4horizontalgrad.pickle', 'wb') as outfile:
    ###pickle.dump([staggeredsurfaces,neighbors,lookups], outfile)

####with open('data/ready4horizontalgrad.pickle', 'rb') as outfile:
    ####[staggeredsurfaces,neighbors,lookups]=pickle.load(outfile)

###nstools.surfaceDiagnostic(staggeredsurfaces)
#staggeredsurfaces = nstools.addHorizontalGrad(staggeredsurfaces,neighbors,lookups)
##nstools.surfaceDiagnostic(staggeredsurfaces)
#staggeredsurfaces = nstools.addVerticalGrad(staggeredsurfaces)
##nstools.surfaceDiagnostic(staggeredsurfaces)
##ptools.saveBathVarTermCache(staggeredsurfaces,"data/bathVar.pickle")
#staggeredsurfaces = nstools.addK(staggeredsurfaces,"data/bathVar.pickle")

##graph.graphSurfaces(staggeredsurfaces,"t")
##graph.graphSurfaces(staggeredsurfaces,"CKVB",show=False,savepath = "refpics/CKs/")
##graph.graphSurfaces(staggeredsurfaces,"t")
##graph.graphSurfaces(staggeredsurfaces,"s")

with open('data/ready4inverse.pickle', 'rb') as outfile:
    [staggeredsurfaces,neighbors,lookups]=pickle.load(outfile)

#nstools.surfaceDiagnostic(staggeredsurfaces)
#with open('data/ready4inverse.pickle', 'wb') as outfile:
    #pickle.dump([staggeredsurfaces,neighbors,lookups], outfile)

coupleinvert,columndictionary,svds,A = inverttools.invert("coupled",staggeredsurfaces,neighbors,lookups)

coupleinvert = nstools.streamFuncToUV(coupleinvert,neighbors,lookups)
coupleinvert = bathtools.addBathToSurface(coupleinvert)


graph.graphVectorField(coupleinvert,"uabs","vabs","z")
graph.graphVectorField(coupleinvert,"u","v","z")
