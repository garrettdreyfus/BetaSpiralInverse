import nstools
import inverttools
import interptools
import bathtools
import saloffset
import graph
import parametertools as ptools
import pickle
import sensitivity
import random
    
#profiles,deepestindex = nstools.extractProfilesBox(["data/1500mprofiles.json"],-180,180,65,90)
#profiles,deepestindex = nstools.removeNorwegianSea(profiles)

#fileObject = open("data/1500NoNorwegian.pickle",'rb')  
#offsets,badfiles,beepestindex = pickle.load(fileObject)

#profiles = nstools.filterCruises(profiles,offsets.keys())
#profiles = saloffset.applyOffsets(profiles,offsets)

#profilechoice = random.choice(nstools.profileInBox(profiles,-180,180,85,90))
##profilechoice = nstools.getProfileById(profiles,"286364")

#surfaces = nstools.runPeerSearch(profiles,deepestindex,200,4000,200,profilechoice,1000)
#with open('data/annotatedprofiles.pickle', 'wb') as outfile:
    #pickle.dump(profiles, outfile)

#with open('data/annotatedprofiles.pickle', 'rb') as outfile:
    #profiles=pickle.load(outfile)


#graph.tsNeutralExplore(profiles)
######fileObject = open(str(profilechoice.eyed)+"new.pickle",'wb')  
####### load the object from the file into var b
######b = pickle.dump(surfaces,fileObject)  
######fileObject.close()

#with open('data/286364new.pickle', 'rb') as outfile:
    #surfaces=pickle.load(outfile)

#surfaces = nstools.addDataToSurfaces(profiles,surfaces,2)
#surfaces = nstools.addStreamFunc(surfaces,profiles)

#surfaces =interptools.addXYToSurfaces(surfaces)
#interpolatedsurfaces,neighbors,lookups = interptools.interpolateSurfaces(surfaces)
#staggeredsurfaces = nstools.fillOutEmptyFields(interpolatedsurfaces)
#staggeredsurfaces = nstools.addHeight(staggeredsurfaces)

#staggeredsurfaces = nstools.addHorizontalGrad(staggeredsurfaces,neighbors,lookups)
#staggeredsurfaces = nstools.addVerticalGrad(staggeredsurfaces)
#ptools.saveBathVarTermCache(staggeredsurfaces,"data/bathVar.pickle")
#staggeredsurfaces = nstools.addK(staggeredsurfaces,"data/bathVar.pickle")

with open('data/ready4inverse.pickle', 'rb') as outfile:
    [staggeredsurfaces,neighbors,lookups]=pickle.load(outfile)

###nstools.surfaceDiagnostic(staggeredsurfaces)
##with open('data/ready4inverse.pickle', 'wb') as outfile:
    ##pickle.dump([staggeredsurfaces,neighbors,lookups], outfile)
sensitivity.mConditionVsError("coupled",staggeredsurfaces,neighbors,lookups)
#coupleinvert,columndictionary,svds,A,errors= inverttools.invert("coupled",staggeredsurfaces,neighbors,lookups)

#coupleinvert = nstools.streamFuncToUV(coupleinvert,neighbors,lookups)
#coupleinvert = bathtools.addBathToSurface(coupleinvert)

#graph.graphVectorField(coupleinvert,"uabs","vabs","z")
#with open('data/couplednomix.pickle', 'wb') as outfile:
    #pickle.dump(coupleinvert, outfile)

#with open('data/couplednomix.pickle', 'rb') as outfile:
    #coupleinvert=pickle.load(outfile)

#graph.graphSurfaces(coupleinvert,"psi",idlabels=True)
#graph.barentsTransport(coupleinvert)
#graph.refinedTransport(coupleinvert,2633,1642)
#graph.refinedTransport(coupleinvert,375,2903)

#across basin
#graph.quantityLine(coupleinvert,(-143.04,70.46),(62.06,82.232),"pres",4000,-1)
##Fram strait
#graph.quantityLine(coupleinvert,(-11.46,81.469),(10.93,79.78),"t",2000)
#graph.quantityLine(coupleinvert,(-11.46,81.469),(10.93,79.78),"s",2000)
#graph.quantityLine(coupleinvert,(-11.46,81.469),(10.93,79.78),"pv",2000)

##Lomonosov Ridge
#graph.transportLine(coupleinvert,(143.03,79.89),(-47.91,84.04),2000)

##Fram strait
#graph.transportLine(coupleinvert,(-24.71,83.60),(27.496,80.06),2400,False)
#graph.transportLine(coupleinvert,(-24.71,83.60),(27.496,80.06),2400,True)
#graph.quantityLine(coupleinvert,(-24.71,83.60),(27.496,80.06),"h",2600)

##Barents
#graph.transportLine(coupleinvert,(-24.30,80.63),(94.79,81.26),1000)
#graph.quantityLine(coupleinvert,(-24.30,80.63),(94.79,81.26),"pres",1000)

##Bering
#graph.transportLine(coupleinvert,(-157.06,71.36),(170.98,69.99),1000)
    


#graph.framStraitTransport(coupleinvert)

#graph.graphVectorField(coupleinvert,"uabs","vabs","pv")
##graph.graphVectorField(coupleinvert,"uabs","vabs","z",show=False,savepath="refpics/nomixing/3pointbath/")
##graph.graphVectorField(coupleinvert,"uabs","vabs","pv",show=False,savepath="refpics/nomixing/3pointpv/")
#graph.graphVectorField(coupleinvert,"u","v","z")
