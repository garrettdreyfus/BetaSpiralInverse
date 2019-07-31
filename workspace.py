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
    
#profiles,deepestindex = nstools.extractProfilesBox(["data/1500mprofiles.json"],-180,180,65,90)
#profiles,deepestindex = nstools.removeNorwegianSea(profiles)

#fileObject = open("data/1500NoNorwegian.pickle",'rb')  
#offsets,badfiles,beepestindex = pickle.load(fileObject)

#profiles = nstools.filterCruises(profiles,offsets.keys())
#profiles = saloffset.applyOffsets(profiles,offsets)

##profilechoice = random.choice(nstools.profileInBox(profiles,-180,180,85,90))
##profilechoice = nstools.getProfileById(profiles,"286364")

##surfaces = nstools.runPeerSearch(profiles,deepestindex,200,4000,200,profilechoice,1000)
##fileObject = open(str(profilechoice.eyed)+"new.pickle",'wb')  
### load the object from the file into var b
##b = pickle.dump(surfaces,fileObject)  
##fileObject.close()

#with open('data/286364new.pickle', 'rb') as outfile:
    #surfaces=pickle.load(outfile)

#surfaces = nstools.addDataToSurfaces(profiles,surfaces,2)
#surfaces = nstools.addStreamFunc(surfaces,profiles)

#surfaces =nstools.addXYToSurfaces(surfaces)
#interpolatedsurfaces,neighbors,lookups = nstools.interpolateSurfaces(surfaces)

#staggeredsurfaces = nstools.fillOutEmptyFields(interpolatedsurfaces)
#staggeredsurfaces = nstools.addHeight(staggeredsurfaces)

#with open('data/ready4horizontalgrad.pickle', 'wb') as outfile:
    #pickle.dump([staggeredsurfaces,neighbors,lookups], outfile)

#with open('data/ready4horizontalgrad.pickle', 'rb') as outfile:
    #[staggeredsurfaces,neighbors,lookups]=pickle.load(outfile)

#nstools.surfaceDiagnostic(staggeredsurfaces)
#staggeredsurfaces = nstools.addHorizontalGrad(staggeredsurfaces,neighbors,lookups)
#nstools.surfaceDiagnostic(staggeredsurfaces)
#staggeredsurfaces = nstools.addVerticalGrad(staggeredsurfaces)
#nstools.surfaceDiagnostic(staggeredsurfaces)
##ptools.saveBathVarTermCache(staggeredsurfaces,"data/bathVar.pickle")
#staggeredsurfaces = nstools.addK(staggeredsurfaces,"data/bathVar.pickle")


#graph.graphSurfaces(staggeredsurfaces,"CKVB",show=False,savepath = "refpics/CKs/")
#graph.graphSurfaces(staggeredsurfaces,"pv")
#graph.graphSurfaces(staggeredsurfaces,"t")
#graph.graphSurfaces(staggeredsurfaces,"s")

with open('data/ready4inverse.pickle', 'rb') as outfile:
    [staggeredsurfaces,neighbors,lookups]=pickle.load(outfile)

#with open('data/ready4inverse.pickle', 'wb') as outfile:
    #pickle.dump([staggeredsurfaces,neighbors,lookups], outfile)

#realinvert = inverttools.invert("coupled",staggeredsurfaces)
#simpleinverted = inverttools.invert("simple",staggeredsurfaces)
#graph.graphVectorField(simpleinverted,"uabs","vabs","s")
#saltinverted = inverttools.invert("simplesalt",staggeredsurfaces)
#graph.graphVectorField(saltinverted,"uabs","vabs","s")
#complexinverted = inverttools.invert("complex",staggeredsurfaces)
#graph.graphVectorField(complexinverted,"uabs","vabs","s")

#prime,coldict = inverttools.invert("coupled",staggeredsurfaces,neighbors,lookups)

#with open('data/coupledoutput.pickle', 'wb') as outfile:
    #pickle.dump([prime,coldict], outfile)

with open('data/coupledoutput.pickle', 'rb') as outfile:
    [prime,coldict]=pickle.load(outfile)

print(np.max(prime),np.min(prime))
for k in staggeredsurfaces.keys():
    staggeredsurfaces[k]["data"]["psinew"] =  np.full(len(staggeredsurfaces[k]["lons"]),np.nan)
    for i in range(len(staggeredsurfaces[k]["data"]["ids"])):
        eyed = staggeredsurfaces[k]["data"]["ids"][i] 
        if eyed in coldict.keys():
            print("modify")
            staggeredsurfaces[k]["data"]["psinew"][i] = staggeredsurfaces[k]["data"]["psi"][i] + prime[coldict[eyed]]

staggeredsurfaces = nstools.streamFuncToUV(staggeredsurfaces,neighbors,lookups)


graph.graphVectorField(staggeredsurfaces,"uabs","v","pv")

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

#graph.graphVectorField(simpleinverted,"u","v","s")
#graph.graphSurfaces(staggeredsurfaces,"v",show=False,savepath = "refpics/staggeredDerivatives/")
#graph.twentyRandomSpirals(staggeredsurfaces)

