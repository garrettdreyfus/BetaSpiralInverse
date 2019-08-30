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

#deepestindex=nstools.deepestProfile(profiles)
#graph.plotProfiles(profiles,"ecco",specialprofile=profiles[deepestindex])

#profiles = nstools.filterCruises(profiles,offsets.keys())
#profiles = saloffset.applyOffsets(profiles,offsets)

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
#graph.tsNeutralExplore(profiles)
#######fileObject = open(str(profilechoice.eyed)+"new.pickle",'wb')  
######## load the object from the file into var b
#######b = pickle.dump(surfaces,fileObject)  
#######fileObject.close()

#with open('data/286364new.pickle', 'rb') as outfile:
  #surfaces=pickle.load(outfile)

#surfaces = nstools.addDataToSurfaces(profiles,surfaces,2)
#surfaces = nstools.addStreamFunc(surfaces,profiles)
#nstools.surfaceDiagnostic(surfaces)

#surfaces =interptools.addXYToSurfaces(surfaces)
#surfaces,neighbors,distances = interptools.interpolateSurfaces(surfaces)

#surfaces = eccotools.addModelEccoUV(surfaces)
#surfaces = nstools.addBathAndMask(surfaces)
###with open('data/betteruv.pickle', 'wb') as outfile:
    ###pickle.dump([surfaces,neighbors,distances], outfile)

#surfaces = nstools.fillOutEmptyFields(surfaces)
#surfaces = nstools.addHeight(surfaces)

#surfaces = nstools.addHorizontalGrad(surfaces,neighbors,distances)
#surfaces = nstools.addVerticalGrad(surfaces)
#ptools.saveBathVarTermCache(surfaces,"data/bathVarecco.pickle")
#surfaces = nstools.addK(surfaces,"data/bathVarecco.pickle")

with open('data/ready4inverseecco.pickle', 'rb') as outfile:
    [surfaces,neighbors,distances]=pickle.load(outfile)
surfaces = nstools.addBathAndMask(surfaces)

#with open('data/ready4inverseecco.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,neighbors,distances], outfile)


#for q in surfaces[1000]["data"].keys():
    #graph.graphSurfaces(surfaces,q,savepath="refpics/eccoallquantsbathmasked/",show=False)


#sensitivity.mixSens("coupled",surfaces,neighbors,distances,savepath="refpics/sensitivity/bathmasked")
#sensitivity.conditionErrorRefLevel("coupled",surfaces,neighbors,distances)
#sensitivity.conditionError("coupled",staggeredsurfaces,neighbors,distances)

params = {"reflevel":600,"upperbound":600,"lowerbound":2200,"mixs":[True,False,True],"debug":True}

inv,columndictionary,svds,A,errors,metadata= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)

inv = nstools.streamFuncToUV(inv,neighbors,distances)
#nstools.surfaceDiagnostic(surfaces)
graph.graphVectorField(inv,"uabs","vabs","z")
##surfaces = eccotools.addSSHToSurface(surfaces)
##graph.graphSurfaces(surfaces,"ssh")
##graph.graphVectorField(inv,"u","v","pv",savepath="refpics/vectorfields/geostroph/")
