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

##profiles = nstools.filterCruises(profiles,offsets.keys())
##profiles = saloffset.applyOffsets(profiles,offsets)

###profilechoice = random.choice(nstools.profileInBox(profiles,-180,180,85,90))
#profilechoice = nstools.getProfileById(profiles,"286364")
#graph.plotProfiles(profiles,"UDASH DATA",nstools.getProfileById(profiles,"286364"))
#surfaces = nstools.runPeerSearch(profiles,deepestindex,200,4000,200,profilechoice,1000)

#with open('data/annotatedprofiles.pickle', 'wb') as outfile:
    #pickle.dump(profiles, outfile)

#with open('data/newsurfaces.pickle', 'wb') as outfile:
    #pickle.dump(surfaces, outfile)

with open('data/annotatedprofiles.pickle', 'rb') as outfile:
    profiles=pickle.load(outfile)

with open('data/newsurfaces.pickle', 'rb') as outfile:
    surfaces=pickle.load(outfile)

##graph.tsNeutralExplore(profiles)
#######fileObject = open(str(profilechoice.eyed)+"new.pickle",'wb')  
######## load the object from the file into var b
#######b = pickle.dump(surfaces,fileObject)  
#######fileObject.close()

#with open('data/286364new.pickle', 'rb') as outfile:
  #surfaces=pickle.load(outfile)

#surfaces = nstools.addDataToSurfaces(profiles,surfaces,2)
#surfaces = nstools.addStreamFunc(surfaces,profiles)

#surfaces =interptools.addXYToSurfaces(surfaces)
#surfaces,neighbors,distances = interptools.interpolateSurfaces(surfaces)
#surfaces = nstools.addBathAndMask(surfaces)
#surfaces = nstools.fillOutEmptyFields(surfaces)
#surfaces = nstools.addHeight(surfaces)

#surfaces = nstools.addHorizontalGrad(surfaces,neighbors,distances)
#surfaces = nstools.addVerticalGrad(surfaces)
#ptools.saveBathVarTermCache(surfaces,"data/bathVar.pickle")
#surfaces = nstools.addK(surfaces,"data/bathVar.pickle")

with open('data/ready4inverse.pickle', 'rb') as outfile:
    [surfaces,neighbors,distances]=pickle.load(outfile)

#nstools.surfaceDiagnostic(surfaces)
#with open('data/ready4inverse.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,neighbors,distances], outfile)
#for q in surfaces[1000]["data"].keys():
    #graph.graphSurfaces(surfaces,q,savepath="refpics/allquants1000/",show=False)
nstools.surfaceDiagnostic(surfaces)
#sensitivity.mixSens("coupled",staggeredsurfaces,neighbors,distances,savepath="refpics/fullmixexplore/")
#sensitivity.conditionErrorRefLevel("coupled",surfaces,neighbors,distances)
#sensitivity.conditionError("coupled",staggeredsurfaces,neighbors,distances)
#graph.graphSurfaces(surfaces,"khp")
params = {"reflevel":1000,"upperbound":1000,"lowerbound":2000,"mixs":[True,False,True],"debug":True}
inv,columndictionary,svds,A,errors,metadata= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)

coupleinvert = nstools.streamFuncToUV(inv,neighbors,distances)

#graph.graphSurfaces(inv,"e")
graph.graphVectorField(inv,"uabs","vabs","z",metadata=metadata)


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

