import nstools
import inverttools
import eccotools
import interptools
import bathtools
import saloffset
import graph
import parametertools as ptools
import pickle
import sensitivity
import random
    
#with open('data/eccoprofiles.pickle', 'rb') as outfile:
    #profiles = pickle.load(outfile)

##graph.plotProfiles(profiles,"yo",data="t",depth="200")

#deepestindex=nstools.deepestProfile(profiles)
#graph.plotProfiles(profiles,"ecco",specialprofile=profiles[deepestindex])

##profiles = nstools.filterCruises(profiles,offsets.keys())
##profiles = saloffset.applyOffsets(profiles,offsets)

#profilechoice = random.choice(nstools.profileInBox(profiles,-180,180,85,90))
#profilechoice = nstools.getProfileById(profiles,profiles[deepestindex].eyed)

#surfaces = nstools.runPeerSearch(profiles,deepestindex,200,4000,200,profilechoice,1000)

#with open('data/eccosurfaces.pickle', 'wb') as outfile:
  #pickle.dump(surfaces, outfile)

#with open('data/eccoannotatedprofiles.pickle', 'wb') as outfile:
    #pickle.dump(profiles, outfile)

with open('data/eccoannotatedprofiles.pickle', 'rb') as outfile:
    profiles=pickle.load(outfile)

with open('data/eccosurfaces.pickle', 'rb') as outfile:
    surfaces=pickle.load(outfile)
#graph.graphSurfaces(surfaces,"pres")
###graph.tsNeutralExplore(profiles)

#surfaces = nstools.addDataToSurfaces(profiles,surfaces,2)
##with open('data/eccosurfwithd.pickle', 'wb') as outfile:
    ##pickle.dump(surfaces, outfile)

##with open('data/eccosurfwithd.pickle', 'rb') as outfile:
    ##surfaces=pickle.load(outfile)

#surfaces = nstools.addStreamFunc(surfaces,profiles)
#nstools.surfaceDiagnostic(surfaces)

#surfaces =interptools.addXYToSurfaces(surfaces)
#surfaces,neighbors,distances = interptools.interpolateSurfaces(surfaces,gaminterpolate=False)

#with open('data/prehist.pickle', 'wb') as outfile:
    #pickle.dump(profiles, outfile)

#graph.distanceHist(distances)

#with open('data/prebathmask.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,neighbors,distances], outfile)

with open('data/prebathmask.pickle', 'rb') as outfile:
    surfaces,neighbors,distances =pickle.load(outfile)

surfaces = nstools.addBathAndMask(surfaces,neighbors)
graph.graphSurfaces(surfaces,"psi")


#surfaces = eccotools.addModelEccoUV(surfaces)




surfaces = nstools.fillOutEmptyFields(surfaces)
surfaces = nstools.addHeight(surfaces)


#with open('data/eccopregradient.pickle', 'rb') as outfile:
    #surfaces,neighbors,distances =pickle.load(outfile)

surfaces = nstools.addHorizontalGrad(surfaces,neighbors,distances)
graph.graphSurfaces(surfaces,"psi")
surfaces = nstools.addVerticalGrad(surfaces)
ptools.saveBathVarTermCache(surfaces,"data/bathVarecco.pickle")
surfaces = nstools.addK(surfaces,"data/bathVarecco.pickle")

#with open('data/ready4inverseecco.pickle', 'rb') as outfile:
    #[surfaces,neighbors,distances]=pickle.load(outfile)


#with open('data/ready4inverseecco.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,neighbors,distances], outfile)

graph.graphSurfaces(surfaces,"psi")
graph.graphVectorField(surfaces,"u","v","psi",\
        savepath="refpics/vectorfields/eccogeo")

#for q in surfaces[1000]["data"].keys():
    #graph.graphSurfaces(surfaces,q,savepath="refpics/eccoallquantsbathmasked/",show=False)


#sensitivity.mixSens("coupled",surfaces,neighbors,distances,savepath="refpics/sensitivity/bathmasked")
#sensitivity.conditionErrorRefLevel("coupled",surfaces,neighbors,distances)
#sensitivity.conditionError("coupled",staggeredsurfaces,neighbors,distances)

params = {"reflevel":600,"upperbound":600,"lowerbound":2200,"mixs":[True,False,True],"debug":True}

inv,columndictionary,svds,A,errors,metadata= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)

inv = nstools.streamFuncToUV(inv,neighbors,distances)
#nstools.surfaceDiagnostic(surfaces)
#graph.graphVectorField(inv,"uabs","vabs","z")
##surfaces = eccotools.addSSHToSurface(surfaces)
##graph.graphSurfaces(surfaces,"ssh")
##graph.graphVectorField(inv,"u","v","pv",savepath="refpics/vectorfields/geostroph/")
