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
##offsets,badfiles,beepestindex = pickle.load(fileObject)

###profiles = nstools.filterCruises(profiles,offsets.keys())
###profiles = saloffset.applyOffsets(profiles,offsets)

####profilechoice = random.choice(nstools.profileInBox(profiles,-180,180,85,90))
#profilechoice = nstools.getProfileById(profiles,"286364")
##graph.plotProfiles(profiles,"UDASH DATA",nstools.getProfileById(profiles,"286364"))
#preinterpsurfaces = nstools.runPeerSearch(profiles,range(200,3001,200),profilechoice,1000)

#with open('data/preinterparctic.pickle', 'wb') as outfile:
    #pickle.dump([preinterpsurfaces,profiles],outfile)
with open('data/preinterparctic.pickle', 'rb') as outfile:
    preinterpsurfaces,profiles = pickle.load(outfile)

#graph.graphSurfaces(preinterpsurfaces,"pres",region="arctic")
preinterpsurfaces = nstools.addDataToSurfaces(profiles,preinterpsurfaces,2)


surfaces,neighbors,distances = interptools.interpolateSurfaces(preinterpsurfaces,\
        fixedgrid="arctic",interpmethod="gam")


nstools.surfaceDiagnostic(surfaces)
#graph.graphSurfaces(surfaces,"pres",region="arctic")


with open('data/postinterparctic.pickle', 'wb') as outfile:
    pickle.dump([surfaces,neighbors,distances], outfile)

with open('data/postinterparctic.pickle', 'rb') as outfile:
    surfaces,neighbors,distances = pickle.load(outfile)

surfaces = nstools.addParametersToSurfaces(surfaces,\
        neighbors,distances,"arctic",[])

with open('data/ready4inversearctic.pickle', 'wb') as outfile:
    pickle.dump([surfaces,neighbors,distances], outfile)

with open('data/ready4inversearctic.pickle', 'rb') as outfile:
    surfaces,neighbors,distances = pickle.load(outfile)



#with open('data/annotatedprofiles.pickle', 'wb') as outfile:
    #pickle.dump(profiles, outfile)

#with open('data/newsurfaces.pickle', 'wb') as outfile:
    #pickle.dump(surfaces, outfile)

#with open('data/annotatedprofiles.pickle', 'rb') as outfile:
    #profiles=pickle.load(outfile)

#with open('data/newsurfaces.pickle', 'rb') as outfile:
    #surfaces=pickle.load(outfile)

###graph.tsNeutralExplore(profiles)
########fileObject = open(str(profilechoice.eyed)+"new.pickle",'wb')  
######### load the object from the file into var b
########b = pickle.dump(surfaces,fileObject)  
########fileObject.close()

##with open('data/286364new.pickle', 'rb') as outfile:
  ##surfaces=pickle.load(outfile)
#with open('data/ready4inverse.pickle', 'rb') as outfile:
    #[surfaces,neighbors,distances]=pickle.load(outfile)

#nstools.surfaceDiagnostic(surfaces)
#with open('data/ready4inverse.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,neighbors,distances], outfile)
#for q in surfaces[1000]["data"].keys():
    #graph.graphSurfaces(surfaces,q,savepath="refpics/allquants1000/",show=False)

params = {"reflevel":1600,"upperbound":1000,"lowerbound":3500,\
        "mixs":{"kvo":True,"kvb":True,"kh":True},"debug":False,\
        "3point":True,"edgeguard":False}


out= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)

with open('data/inverseout.pickle', 'wb') as outfile:
    pickle.dump([out,neighbors,distances], outfile)
with open('data/inverseout.pickle', 'rb') as outfile:
    [out,neighbors,distances] = pickle.load(outfile)

###print(out["metadata"])
#graph.graphSurfaces(out["surfaces"],"pres",region="nepb")
inv = nstools.streamFuncToUV(out["surfaces"],neighbors,distances)


#graph.graphSurfaces(inv,"e")
graph.graphVectorField(inv,"uabs","vabs","z",metadata=out["metadata"])


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

