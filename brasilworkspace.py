
from regionlib import brasil
import os, gsw
import graph
import nstools
import pickle
import interptools
import inverttools
from functools import partial
import parametertools as ptools
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pdb
from scipy.io import savemat
import sensitivity

#profiles = brasil.extractArgoProfiles(os.path.abspath("data/brasilargonc"))
#profiles = brasil.extractBodcProfiles(os.path.abspath("data/brasilmorectd/"))
#print("BODC: ",len(profiles))
# profilesDeep = brasil.extractDeepArgoProfiles(os.path.abspath("data/brasildeepargo.json"))
# print(len(profilesDeep))
# profilesWoce = brasil.extractWoceProfiles(os.path.abspath("data/brasilnc/"))
# for p in profilesWoce:
#     s=p.pres>3500
#     plt.plot(p.sals[s],p.temps[s],color="blue",label="Ship Based")
# for p in profilesDeep:
#     s=p.pres>3500
#     plt.plot(p.sals[s],p.temps[s],color="red",label="Deep Argo")
# plt.legend()
# plt.show()

# profiles = profilesWoce + profilesDeep
# print("WOCE and BODC: ",len(profiles))
# with open('data/argoandwoce.pickle', 'wb') as outfile:
#     pickle.dump(profiles,outfile)
# with open('data/argoandwoce.pickle', 'rb') as outfile:
#    profiles = pickle.load(outfile)

# profilechoice = nstools.profileInBox(profiles,-40,-20,-31,-28,5000)
# #profilechoice = nstools.profileInBox(profiles,-42,-37,-31,-28,4500)
# profilechoice = profilechoice[0]
# print("sals: ", list(gsw.SP_from_SA(profilechoice.sals,profilechoice.pres,profilechoice.lon,profilechoice.lat)))
# print("temps: ", list(gsw.t_from_CT(profilechoice.sals,profilechoice.temps,profilechoice.pres)))
# print("pres: ", list(profilechoice.pres))
# print(profilechoice.lon,profilechoice.lat)
# #graph.plotProfiles(brasil,profiles,"",specialprofile=profilechoice)
# #graph.layerChooser(profilechoice)
# # #print(profilechoice.lat,profilechoice.lon)
# HTLayers = [65,173,360,570,819,1172,1563,1887,2150,2802,3100,3400,3700,4000,4300,4600,4900]
# HTLE = [65]
# for i in range(1,len(HTLayers)-1):
#    HTLE.append(int(HTLayers[i] - (HTLayers[i]-HTLayers[i-1])/3.0))
#    HTLE.append(int(HTLayers[i] + (HTLayers[i+1]-HTLayers[i])/3.0))
# HTLE.append(HTLayers[-1])
# print(HTLE)
# preinterpsurfaces = nstools.runPeerSearch(profiles,HTLE,profilechoice,False,10**10)

# # # # #gammavals = [25.875,26.51,26.862,27.158,27.3605,27.526,27.6575,27.7825,27.8700, \
# # # #         #27.9275,27.965,27.99,28.015,28.03,28.0475,28.0625,28.08,28.108,28.136,28.164,28.2,28.33,28.36]
# # # # ##print(range(100,6000,200)[len(gammavals)-1])
# # # # ##preinterpsurfaces = nstools.runPeerSearch(profiles,range(100,4300,200),profilechoice,False,10**4,gammas=gammavals)
# # # # ##preinterpsurfaces = nstools.runPeerSearch(profiles,range(100,4300,200),profilechoice,True,10**5)

# with open('data/annotatedbrasilargowoce.pickle', 'wb') as outfile:
#     pickle.dump([preinterpsurfaces,profiles],outfile)

# with open('data/annotatedbrasilargowoce.pickle', 'rb') as outfile:
#      preinterpsurfaces,profiles = pickle.load(outfile)
# #graph.NSGAMCompareCruise(preinterpsurfaces,"A10")
# preinterpsurfaces = nstools.addDataToSurfaces(brasil,profiles,preinterpsurfaces)
# #nstools.addGammaN(surfaces)


# with open('data/brasilsurfaceswdata.pickle', 'wb') as infile:
#     pickle.dump([preinterpsurfaces,profiles],infile)
# with open('data/brasilsurfaceswdata.pickle', 'rb') as outfile:
#     preinterpsurfaces,profiles = pickle.load(outfile)
# # print("hi")
# # #graph.time_diagnostic(profiles,3100,-30,2.5)
# # #graph.time_diagnosuItic(profiles,3100,-25,2.5)
# # #graph.time_diagnostic(profiles,3100,-10,2.5)
# # # nserrors = {}
# surfaces,neighbors,distances = \
#     interptools.interpolateSurfaces(brasil,preinterpsurfaces,\
#     interpmethod="gam",smart=False,coord="latlon",splines=40)
# #nstools.neutralityError(surfaces)
# #graph.graphSurfaces(brasil,surfaces,"nserror",stds=2,show=False,savepath="../arcticcirc-pics/surfaces/normneutralerror/")
# #graph.nsHist(surfaces)


# with open('data/interpedbrasil.pickle', 'wb') as outfile:
#     pickle.dump([surfaces,neighbors,distances], outfile)
# with open('data/interpedbrasil.pickle', 'rb') as outfile:
#     [surfaces,neighbors,distances] = pickle.load(outfile)

# #sensitivity.decayScaleSensitivity(surfaces,neighbors,distances)
# surfaces = nstools.addParametersToSurfaces(brasil,surfaces,neighbors,distances)
# # graph.NSGAMCompare(preinterpsurfaces,surfaces,-30,-180,180)

# # print(surfaces.keys())
# params = {"reflevel":int(2062),"upperbound":1000,"lowerbound":4200,\
#         "mixs":{"kvo":True,"kvb":True,"kh":True},"debug":False,\
#           "3point":True,"edgeguard":True}
# # Conditions
# # All mixing: 201235
# # No mixing: 147
# # Kv0 only: 147
# # KvH and Kv0 only: 148
# # KvH and Kv0 only with out edgeguard (tm): 489

# out= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)
# print(out["metadata"])
# inv=out["surfaces"]
# inv = nstools.streamFuncToUV(inv,neighbors,distances)

# with open('data/invertedbrasil.pickle', 'wb') as outfile:
#     pickle.dump([inv,neighbors,distances], outfile)
with open('data/invertedbrasil.pickle', 'rb') as outfile:
    [inv,neighbors,distances] = pickle.load(outfile)
graph.AABWFinder(inv)
print(nstools.regionCurl(inv,3600,-39,-36,-38,-27))
transports = nstools.transportDiagnostics(inv)
print(transports)
# # u = [0,transports["northern"],0,transports["southern"]]
# # v = [transports["vema"],0,transports["hunter"],0]
# # x = [-1,0,1,0]
# # y = [0,1,0,-1]
# # transports = [transports["vema"],transports["northern"],transports["hunter"],transports["southern"]]
# # plt.plot(transports)
# # plt.show()
# #graph.pseudopolzin(inv,-20,-180,180,0,6000,quant="kv")
# #graph.HTtransports(inv,ht="data/hgt.pickle")
# #graph.northSouthTransect(inv,"kv",lat=-20,show=True)
# #nstools.inverseReady(inv)
# graph.meridionalHeatMap(inv,-30,-180,180,1000,6000,show=True,label="")
# result = nstools.transportDiagnostics(inv)
# print(result)
# #graph.latitudinalHeatMap(inv,-35.1,-40,-20,1000,6000,show=True,label="")
# # graph.meridionalSurfaces(inv,-30,-50,-25,1000,6000,show=True,label="")
# # graph.fourpanelVectorField(brasil,inv,"uabs","vabs",backgroundfield="s",\
# #                            select=[1054,1779,3201,4466],transform=False,scale=0.1)
# # graph.graphSurfaces(brasil,inv,"pres",stds=2,show=True,select=range(3200,5000),contour=True)
# # kvs = []
# # for l in inv.keys():
# #     kvs += list([np.nanmean(inv[l]["data"]["kh"])])
# # print("kh: ,",np.nanmean(np.asarray(kvs).flatten()))
# # kvs = []
# # for l in inv.keys():
# #     kvs += list([np.nanmean(inv[l]["data"]["kv"])])
# # print("kv: ,",np.nanmean(np.asarray(kvs).flatten()))
# # kvs = []
# # for l in inv.keys():
# #     kvs += list([np.nanmean(inv[l]["data"]["kvo"])])
# # print("kvo: ,",np.nanmean(np.asarray(kvs).flatten()))
# # #graph.graphSurfaces(brasil,inv,"bathvar",stds=1,show=True,select=range(200,5000))

# # with open('data/interpedbrasil.pickle', 'rb') as outfile:
# #     [surfaces,neighbors,distances] = pickle.load(outfile)
# # graph.graphVectorField(brasil,inv,"uabs","vabs","pres",# \
# #                         metadata=out["metadata"],\
# #                         transform=False,show=False,
# #                         savepath="../arcticcirc-pics/vectorfields/farthereast/",scale=0.1)
# ################### SAVE WITH DIFFERENT REF LEVELS
# with open('data/interpedbrasil.pickle', 'rb') as outfile:
#     [surfaces,neighbors,distances] = pickle.load(outfile)
# #graph.graphSurfaces(brasil,surfaces,"psi",stds=1,show=True,select=range(3000,5000))
# reflevels = surfaces.keys()
# levels = []
# hunters = []
# vemas = []
# conditions=[]
# curls = []
# errors = []
# for k in reflevels:
#     if int(k) >1000 and int(k) < 4000:
#         with open('data/interpedbrasil.pickle', 'rb') as outfile:
#             [surfaces,neighbors,distances] = pickle.load(outfile)
#         params = {"reflevel":int(k),"upperbound":1000,"lowerbound":4000,\
#                 "mixs":{"kvo":True,"kvb":True,"kh":True},"debug":False,\
#                 "3point":True,"edgeguard":True}
#         out= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)
#         inv = out["surfaces"]
#         inv = nstools.streamFuncToUV(inv,neighbors,distances)
#         graph.graphVectorField(brasil,inv,"uabs","vabs","pres",# \
#                                 metadata=out["metadata"],select=[3600],\
#                                 transform=False,show=False,scale=0.1,\
#                                 savepath="../articcirc-pics/vectorfields/curls/"+str(k))


#         levels.append(k)
#         result = nstools.transportDiagnostics(inv)
#         hunters.append(result["hunter"])
#         vemas.append(result["vema"])
#         curls.append(result["curl"])
#         conditions.append(out["metadata"]["condition"])
#         errors.append(out["metadata"]["error"])
# fig,((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2,3)
# ax1.plot(levels,conditions)
# ax2.plot(levels,errors)
# ax4.plot(levels,vemas)
# ax5.plot(levels,hunters)
# ax6.plot(levels,curls)
# plt.show()



# ################### PLOT DIFFERENT REF LEVELS
# with open('data/interpedbrasil.pickle', 'rb') as outfile:
#     [surfaces,neighbors,distances] = pickle.load(outfile)
# levels = surfaces.keys()
# for k in levels:
#     if int(k) >1000 and int(k) < 4000:
#         with open('data/reflevelchoice/{}.pickle'.format(k), 'rb') as outfile:
#             [out,neighbors,distances] = pickle.load(outfile)
#         inv = out["surfaces"]
#         inv = nstools.streamFuncToUV(inv,neighbors,distances)
#         graph.graphVectorField(brasil,inv,"uabs","vabs","pres",# \
#                                 metadata=out["metadata"],\
#                                 transform=False,show=False,
#                                 savepath="../articcirc-pics/vectorfields/reflevelchoice/{}".format(k),scale=0.1)




# inv = out["surfaces"]
# inv = nstools.streamFuncToUV(inv,neighbors,distances)
#for lat in range(-30,-2,5):
    #graph.northSouthTransect(inv,"uabs",lat=lat,savepath="refpics/transects/noiseexpiriment/"+str(n)+"stdev"+str(i)+"/")
#for lon in range(-35,-12,5):
    #graph.northSouthTransect(inv,"vabs",lon=lon,savepath="refpics/transects/noiseexpiriment/"+str(n)+"stdev"+str(i)+"/")
#graph.graphSurfaces(brasil,inv,"gamma",stds=2,savepath="refpics/surfaces/gammaderivedsurface/",show=False)
# graph.graphVectorField(brasil,out,"uabs","vabs","pres",\
#         transform=False)
# graph.graphVectorField(brasil,inv,"uabs","vabs","pres",# \
        # metadata=out["metadata"],\
        # transform=False,show=False,
        # savepath="refpics/vectorfields/brasilwoceargo/",scale=0.1)

