import nstools
import inverttools
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
    
profiles,deepestindex = nstools.extractProfilesBox(["data/1500mprofiles.json"],-180,180,65,90)
profiles,deepestindex = nstools.removeNorwegianSea(profiles)

fileObject = open("data/1500NoNorwegian.pickle",'rb')  
offsets,badfiles,beepestindex = pickle.load(fileObject)

profiles = nstools.filterCruises(profiles,offsets.keys())
profiles = saloffset.applyOffsets(profiles,offsets)

#profilechoice = random.choice(nstools.profileInBox(profiles,-180,180,85,90))

#surfaces = nstools.runPeerSearch(profiles,deepestindex,200,4000,200,profilechoice,1000)
#fileObject = open(str(profilechoice.eyed)+".pickle",'wb')  
## load the object from the file into var b
#b = pickle.dump(surfaces,fileObject)  
#fileObject.close()

with open('data/287959.pickle', 'rb') as outfile:
    surfaces=pickle.load(outfile)

surfaces = nstools.addDataToSurfaces(profiles,surfaces,2)
surfaces = nstools.addStreamFunc(surfaces,profiles)

surfaces =nstools.addXYToSurfaces(surfaces)
interpolatedsurfaces,neighbors,lookups = nstools.interpolateSurfaces(surfaces)
#graph.graphSurfacesComparison(interpolatedsurfaces,surfaces,"psi",show=False,savepath="refpics/VersionTwo/")

#graph.graphSurfacesComparison(interpolatedsurfaces,surfaces,"psi")
#graph.graphAndSaveAll(interpolatedsurfaces,"refpics/VersionTwo/")

staggeredsurfaces = nstools.fillOutEmptyFields(interpolatedsurfaces)
staggeredsurfaces = nstools.addHeight(staggeredsurfaces)
staggeredsurfaces = nstools.addHorizontalGrad(staggeredsurfaces,neighbors,lookups)
staggeredsurfaces = nstools.addVerticalGrad(staggeredsurfaces)
staggeredsurfaces = nstools.addK(staggeredsurfaces,"data/bathVar.pickle")
nstools.surfaceDiagnostic(staggeredsurfaces)

#ptools.saveBathVarTermCache(staggeredsurfaces,"data/bathVar.pickle")

#graph.graphSurfaces(staggeredsurfaces,"CKVB",show=False,savepath = "refpics/CKs/")
#graph.graphSurfaces(staggeredsurfaces,"pv")
#graph.graphSurfaces(staggeredsurfaces,"t")
#graph.graphSurfaces(staggeredsurfaces,"s")

#with open('data/ready4inverse.pickle', 'rb') as outfile:
    #staggeredsurfaces=pickle.load(outfile)

#graph.graphSurfaces(staggeredsurfaces,"psi")
#simpleinvert = inverttools.invert("simple",staggeredsurfaces)
#simplesaltinvert = inverttools.invert("simplesalt",staggeredsurfaces)
#complexsaltinvert = inverttools.invert("complexsalt",staggeredsurfaces)
with open('data/ready4inverse.pickle', 'wb') as outfile:
    pickle.dump(staggeredsurfaces, outfile)

staggeredsurfaces = inverttools.invert("complex",staggeredsurfaces)

#graph.graphVectorField(simpleinvert,"uabs","vabs",savepath="refpics/VersionTwo/simple/",show=False)
#graph.graphVectorField(simplesaltinvert,"uabs","vabs",savepath="refpics/VersionTwo/simplesalt/",show=False)
#graph.graphVectorField(complexsaltinvert,"uabs","vabs",savepath="refpics/VersionTwo/complexsalt/",show=False)
#graph.graphVectorField(simplesaltinvert,"u","v")
##$$$$$$$$$$$$$$$44
#staggeredsurfaces = inverttools.invert("complexsalt",staggeredsurfaces,debug=False)
##$$$$$$$$$$$$$$$44

#with open('data/svdinverted.pickle', 'wb') as outfile:
    #pickle.dump(staggeredsurfaces, outfile)

#with open('data/svdinverted.pickle', 'rb') as outfile:
    #staggeredsurfaces=pickle.load(outfile)

#graph.graphSurfaces(staggeredsurfaces,"u",show=False,savepath = "refpics/ref1000/")
#graph.graphSurfaces(staggeredsurfaces,"v",show=False,savepath = "refpics/ref1000/")
#graph.graphSurfaces(staggeredsurfaces,"uabs",show=False,savepath = "refpics/ref1000/")
#graph.graphSurfaces(staggeredsurfaces,"vabs",show=False,savepath = "refpics/ref1000/")
#graph.graphSurfaces(staggeredsurfaces,"uprime",show=False,savepath = "refpics/ref1000/")
#graph.graphSurfaces(staggeredsurfaces,"vprime",show=False,savepath = "refpics/ref1000/")
#graph.graphSurfaces(staggeredsurfaces,"hx",show=False,savepath = "refpics/ref1000/")
#graph.graphSurfaces(staggeredsurfaces,"hy",show=False,savepath = "refpics/ref1000/")
#graph.graphSurfaces(staggeredsurfaces,"h",show=False,savepath = "refpics/ref1000/")
#graph.graphSurfaces(staggeredsurfaces,"CKVB")
#graph.graphSurfaces(staggeredsurfaces,"dsdx")
#graph.graphSurfaces(staggeredsurfaces,"d2sdx2")
#graph.graphSurfaces(staggeredsurfaces,"vabs")
#graph.graphVectorField(staggeredsurfaces,"u","v",show=False,savepath="refpics/NormRowsREF1000/simplesalt/")
#graph.graphVectorField(staggeredsurfaces,"uabs","vabs",show=False,savepath="refpics/NormRowsREF1000/simplesalt/")
#graph.graphVectorField(staggeredsurfaces,"uabs","vabs")
graph.graphVectorField(staggeredsurfaces,"uabs","vabs")
graph.graphVectorField(staggeredsurfaces,"u","v")
#graph.graphSurfaces(staggeredsurfaces,"v",show=False,savepath = "refpics/staggeredDerivatives/")
#graph.twentyRandomSpirals(staggeredsurfaces)

