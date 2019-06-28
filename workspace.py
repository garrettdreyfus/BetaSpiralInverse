import nstools
import saloffset
import matplotlib.pyplot as plt
import pickle
import random

    

#offsets, profiles, deepestindex = saloffset.runSalinityOffsetTool(["data/3000mprofiles.json"],["ODEN_AGAVE"])
#print(offsets)
#profiles = nstools.filterCruises(profiles,offsets.keys())
#profiles = saloffset.applyOffsets(profiles,offsets)

#fileObject = open("corrected.pickle",'wb')  
## load the object from the file into var b
#b = pickle.dump([offsets,profiles,deepestindex],fileObject)  
#fileObject.close()


##USING CORRECTED PICKLED SURFACES TO MAP NEUTRAL SURFACES USING PEER SEARCH

#fileObject = open("corrected.pickle",'rb')  
#offsets,profiles,deepestindex = pickle.load(fileObject)
#fileObject.close()

#nstools.plotCruise(profiles,"IPY 2007")

##surfaces = nstools.peerSearch(profiles,deepestindex,1000)
#profilechoice = random.choice(nstools.profileInBox(profiles,-180,180,85,90))
#surfaces = {}
#for d in range(200,4000,200)[::-1]:
    #print(d)
    #surfaces.update(nstools.peerSearch(profiles.copy(),deepestindex,d,profilechoice,1000))

#nstools.graphSurfaces(surfaces)

#fileObject = open(str(profilechoice.eyed)+".pickle",'wb')  
## load the object from the file into var b
#b = pickle.dump(surfaces,fileObject)  
#fileObject.close()

##############################################
## 

fileObject = open("287335.pickle",'rb')  
surfaces = pickle.load(fileObject)
fileObject.close()
fileObject = open("corrected.pickle",'rb')  
offsets,profiles,deepestindex = pickle.load(fileObject)
fileObject.close()
tempSurfs = {}
for d in surfaces.keys():
    tempSurf = [[],[],[],[]]
    print(len(surfaces[d][0]))
    for l in range(len(surfaces[d][0])):
        p = nstools.getProfileById(profiles,surfaces[d][3][l])
        t,s = p.atPres(surfaces[d][2][l])
        pv = p.potentialVorticity(surfaces[d][2][l])
        tempSurf[0].append(surfaces[d][0][l])
        tempSurf[1].append(surfaces[d][1][l])
        tempSurf[2].append(surfaces[d][2][l])
        tempSurf[3].append(surfaces[d][3][l])
    if len(tempSurf[0])>5:
        tempSurfs[d] = tempSurf

nstools.graphSurfaces(tempSurfs)









#####################################################################
#profiles,deepestindex = nstools.extractProfiles(["data/3000mprofiles.json"])
#nstools.mapPV(profiles,100)
#fig, ax = plt.subplots(1,1)
#for p in profiles:
    #ax.plot(p.isals,p.itemps,linewidth=1)
#plt.show()
#print(nstools.cruiseCount(profiles))



#print(runSalinityOffsetTool(glob.glob("data/3000m2007profiles.json"),["Polarstern_ARK-XXIII_2"]))
#singleSalinityOffsetRun("data/2000mprofiles.json","LOUIS_S._ST._LAURENT_18SN940","HUDSON_HUDSON2")
#profiles,deepestindex = nstools.extractProfilesMonths("data/3000m2008profiles.json",range(13))
#nstools.plotCruise(nstools.cruiseSearch(profiles,"LOUIS_S._ST._LAURENT_18SN940",1994),"name")

         
#for name in nstools.transArcticSearch(profiles):
    #print(name)
        
    

#print("DONE WITH EXTRACTING PROFILES")
#surfaces = nstools.search(profiles,deepestindex)
#print("DONE FINDING SURFACES")
#with open('data/surfaces.json', 'w') as outfile:
    #json.dump(surfaces, outfile)
#json_file = open("data/surfaces.json") 
#surfaces = json.load(json_file)
#print("NOW GRAPHING")
#nstools.graphSurfaces(surfaces)
