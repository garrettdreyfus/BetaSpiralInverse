from regionlib import brasil, laurent
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

# #profiles = brasil.extractArgoProfiles(os.path.abspath("data/brasilargonc"))
# #profiles = brasil.extractBodcProfiles(os.path.abspath("data/brasilmorectd/"))
# #print("BODC: ",len(profiles))
# profilesDeep = brasil.extractDeepArgoProfiles(os.path.abspath("data/brasildeepargo.json"))
# print(len(profilesDeep))
# profilesWoce = brasil.extractWoceProfiles(os.path.abspath("data/brasilnc/"))

# profiles = profilesWoce #+ profilesDeep
# with open('data/argoandwoce.pickle', 'wb') as outfile:
#     pickle.dump(profiles,outfile)
# with open('data/argoandwoce.pickle', 'rb') as outfile:
#    profiles = pickle.load(outfile)

# # # graph.plotProfiles(brasil,profiles,"p")
# # profilechoice = nstools.profileInBox(profiles,-40,-20,-31,-28,5000)
# # #profilechoice = nstools.profileInBox(profiles,-45,-20,-40,-31,5000)
# profilechoice = nstools.profileInBox(profiles,-25,-20,-31,-28,4500)
# profilechoice = profilechoice[0]
# # print(profilechoice.lon,profilechoice.maplat)
# # #profilechoice = nstools.smoothProfile(profilechoice)

# # # print("sals: ", list(gsw.SP_from_SA(profilechoice.sals,profilechoice.pres,profilechoice.lon,profilechoice.lat)))
# # # print("temps: ", list(gsw.t_from_CT(profilechoice.sals,profilechoice.temps,profilechoice.pres)))
# # # print("pres: ", list(profilechoice.pres))
# # # print(profilechoice.lon,profilechoice.lat)
# # # graph.plotProfiles(brasil,profiles,"",specialprofile=profilechoice)
# # # #graph.layerChooser(profilechoice)
# # # #print(profilechoice.lat,profilechoice.lon)
# # # HTLayers = [65,173,360,570,819,1172,1563,1887,2150,2802,3100,3400,3700,4000,4500]
# # # HTLE = [65]
# # # for i in range(1,len(HTLayers)-1):
# # #    HTLE.append(int(HTLayers[i] - (HTLayers[i]-HTLayers[i-1])/3.0))
# # #    HTLE.append(int(HTLayers[i] + (HTLayers[i+1]-HTLayers[i])/3.0))
# # # HTLE.append(HTLayers[-1])
# # # print(HTLE)
# # # preinterpsurfaces = nstools.runPeerSearch(profiles,HTLE,profilechoice,False,10**10)
# preinterpsurfaces = nstools.runPeerSearch(profiles,range(200,5000,200),profilechoice,False,10**10)
# # preinterpsurfaces = nstools.runPeerSearch(profiles,range(200,5000,200),profilechoice,False,10**4)

# # # gammavals = [25.875,26.51,26.862,27.158,27.3605,27.526,27.6575,27.7825,27.8700, \
# # #       27.9275,27.965,27.99,28.015,28.03,28.0475,28.0625,28.08,28.108,28.136,28.164,28.2,28.33,28.36]
# # # print(range(100,6000,200)[len(gammavals)-1])
# # # preinterpsurfaces = nstools.runPeerSearch(profiles,range(100,4300,200),profilechoice,False,10**4,gammas=gammavals)
# # #
# # # preinterpsurfaces = nstools.runPeerSearch(profiles,range(100,4300,200),profilechoice,True,10**5)

# with open('data/annotatedbrasilargowoce.pickle', 'wb') as outfile:
#     pickle.dump([preinterpsurfaces,profiles],outfile)

# with open('data/annotatedbrasilargowoce.pickle', 'rb') as outfile:
#      preinterpsurfaces,profiles = pickle.load(outfile)

# preinterpsurfaces = nstools.addDataToSurfaces(brasil,profiles,preinterpsurfaces)
# #nstools.addGammaN(surfaces)
# for k in preinterpsurfaces.keys():
#      total = len(preinterpsurfaces[k]["data"]["pv"])
#      negcount = np.nansum(np.asarray(preinterpsurfaces[k]["data"]["pv"])<0)
#      print(k,":",negcount/total)


# with open('data/brasilsurfaceswdata.pickle', 'wb') as infile:
#     pickle.dump([preinterpsurfaces,profiles],infile)
# with open('data/brasilsurfaceswdata.pickle', 'rb') as outfile:
#     preinterpsurfaces,profiles = pickle.load(outfile)

# #graph.graphSurfaces(brasil,preinterpsurfaces,"pv",select=range(2584,10000))
# # #graph.graphSurfaces(brasil,preinterpsurfaces,"pv",select=range(2000,10000),log=True)
# # for k in preinterpsurfaces.keys():
# #       total = len(preinterpsurfaces[k]["data"]["pv"])
# #       negcount = np.nansum(np.asarray(preinterpsurfaces[k]["data"]["pv"])<0)
# #       print(negcount/total)
# # print("*"*5)

# #  #print(preinterpsurfaces[65]["lats"])
# #  #graph.time_diagnostic(profiles,3100,-30,2.5)
# #  #graph.time_diagnosuItic(profiles,3100,-25,2.5)
# #  #graph.time_diagnostic(profiles,3100,-10,2.5)
# #  # nserrors = {}
# # #graph.graphSurfaces(brasil,preinterpsurfaces,"n^2")
# surfaces,neighbors,distances = \
#      interptools.interpolateSurfaces(laurent,preinterpsurfaces,\
#                                      interpmethod="gam",smart=False,coord="latlon",splines=30)

#graph.graphSurfaces(brasil,surfaces,"pv",contour=False,secondsurface=preinterpsurfaces,select=range(2584,10000))
# graph.nsHist(surfaces)
# with open('data/interpedlaurent.pickle', 'wb') as outfile:
#      pickle.dump([surfaces,neighbors,distances], outfile)
#with open('data/interpedlaurent.pickle', 'rb') as outfile:
    #[surfaces,neighbors,distances] = pickle.load(outfile)

## graph.graphSurfaces(brasil,surfaces,"pv")
## surfaces = nstools.neutralityError(surfaces)
## #sensitivity.weightingSensitivity(surfaces,neighbors,distances)
## sensitivity.decayScaleSensitivity(surfaces,neighbors,distances)
## # # print("a",surfaces["maplats"])
#surfaces = nstools.addParametersToSurfaces(laurent,surfaces,neighbors,distances,H_0=500)
## graph.NSGAMCompare(preinterpsurfaces,surfaces,-30.5,-180,180)

#with open('data/withparamslaurent.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,neighbors,distances], outfile)
#with open('data/withparamslaurent.pickle', 'rb') as outfile:
    #[surfaces,neighbors,distances] = pickle.load(outfile)

#params = {"reflevel":1800,"upperbound":1000,"lowerbound":4000,\
       #"mixs":{"kvo":True,"kvb":True,"kh":True},"debug":False,\
           #"3point":True,"edgeguard":True,"H_0":500
         #}
#out= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)
#print(out["metadata"])
#inv=out["surfaces"]
#with open('data/invertedlaurent.pickle', 'wb') as outfile:
   #pickle.dump([inv,neighbors,distances], outfile)
with open('data/invertedlaurent.pickle', 'rb') as outfile:
   [inv,neighbors,distances] = pickle.load(outfile)

for k in inv.keys():
     total = len(inv[k]["data"]["kv"])
     negcount = np.nansum(inv[k]["data"]["kv"]<0)
     m=inv[k]["data"]["kv"]<0
     print(k,":",negcount/total)
inv = nstools.addGammaN(inv)
graph.surfaceGammaRanges(inv)
inv = nstools.streamFuncToUV(inv,neighbors,distances)
# graph.graphVectorField(laurent,inv,"uabs","vabs","pres",\
#                         transform=False,show = False, scale=0.1,savepath="/home/gdreyfus/Projects/arcticcirc-pics/vectorfields/laurent/")
inv = nstools.twoCReference(inv)
#graph.graphSurfaces(brasil,surfaces,"bathvar",select=range(0,10000))
for k in inv.keys():
    plt.scatter(inv[k]["data"]["kv"],inv[k]["data"]["pres"])
plt.show()

k = 4400#list(inv.keys())[20]
print(k)
# plt.errorbar(range(len(inv[k]["data"]["kv"])),inv[k]["data"]["kv"],yerr=2*inv[k]["data"]["kverror"])
# plt.show()

#plt.scatter(inv[k]["data"]["uabs"],inv[k]["data"]["vabs"])
#plt.scatter(inv[k]["data"]["uerror"],inv[k]["data"]["verror"],c='red')
#plt.errorbar(inv[k]["data"]["uabs"],inv[k]["data"]["vabs"],yerr=inv[k]["data"]["verror"]*2,xerr=inv[k]["data"]["uerror"]*2,fmt='none')
###########################################
### error bubble plot##
###########################################
fig, [ax1,ax2,ax3] = plt.subplots(1,3)
d = inv[k]["data"]["vabs"]>0
c = ax1.scatter(inv[k]["lons"][d],inv[k]["maplats"][d],s=np.abs(inv[k]["data"]["vabs"][d])*50000,c=inv[k]["data"]["verror"][d]*2/np.abs(inv[k]["data"]["vabs"][d]),cmap="jet",vmin=0,vmax=2)
d = inv[k]["data"]["vabs"]<=0
c = ax1.scatter(inv[k]["lons"][d],inv[k]["maplats"][d],s=np.abs(inv[k]["data"]["vabs"][d])*50000,c=inv[k]["data"]["verror"][d]*2/np.abs(inv[k]["data"]["vabs"][d]),cmap="jet",vmin=0,vmax=2,marker="x")
plt.colorbar(c,ax=ax1)
c.set_clim(0,2)
ax1.set_title("V")
ax1.set_xlabel("lon")
ax1.set_ylabel("lat")
c = ax2.scatter(inv[k]["lons"],inv[k]["maplats"],s=np.abs(inv[k]["data"]["uabs"])*50000,c=inv[k]["data"]["uerror"]*2/np.abs(inv[k]["data"]["uabs"]),cmap="jet")
c.set_clim(0,2)
plt.colorbar(c,ax=ax2)
ax2.set_title("U")
ax2.set_xlabel("lon")
ax2.set_ylabel("lat")
d = inv[k]["data"]["kv"]>0
c = ax3.scatter(inv[k]["lons"][d],inv[k]["maplats"][d],s=np.abs(inv[k]["data"]["kv"][d])*500000,c=inv[k]["data"]["kverror"][d]*2/np.abs(inv[k]["data"]["kv"][d]),cmap="jet")
c.set_clim(0,2)
d = inv[k]["data"]["kv"]<0
c = ax3.scatter(inv[k]["lons"][d],inv[k]["maplats"][d],s=np.abs(inv[k]["data"]["kv"][d])*500000,c=inv[k]["data"]["kverror"][d]*2/np.abs(inv[k]["data"]["kv"][d]),cmap="jet",marker="x")
c.set_clim(0,2)
plt.colorbar(c,ax=ax3)
ax3.set_title("Kv")
ax3.set_xlabel("lon")
ax3.set_ylabel("lat")
plt.show()
###########################################
###########################################
###########################################

#graph.HTtransports(inv,ht="data/hgt.pickle")
# inv = nstools.externalReference(inv,"data/yomaha_1000.nc")
# print(nstools.transportDiagnostics(inv,["yomahau","yomahav"]))
# print(inv.keys())
# #graph.fourpanelVectorField(brasil,inv,"yomahau","yomahav",backgroundfield="s",\
#                            #select=[1054,1779,3200,4400],transform=False,scale=0.5)

# print(nstools.regionCurl(inv,3600,-39,-36,-38,-27))
# transports = nstools.transportDiagnostics(inv)
# print(transports)
# # u = [0,transports["northern"],0,transports["southern"]]
# # v = [transports["vema"],0,transports["hunter"],0]
# # x = [-1,0,1,0]
# # y = [0,1,0,-1]
# # transports = [transports["vema"],transports["northern"],transports["hunter"],transports["southern"]]
# # plt.plot(transports)
# # plt.show()
# #graph.pseudopolzin(inv,-20,-180,180,0,6000,quant="kv")
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
# with open('data/withparams.pickle', 'rb') as outfile:
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
#     if int(k) >1200 and int(k) < 4000:
#         with open('data/withparams.pickle', 'rb') as outfile:
#             [surfaces,neighbors,distances] = pickle.load(outfile)
#         params = {"reflevel":int(k),"upperbound":1000,"lowerbound":4000,\
#                 "mixs":{"kvo":True,"kvb":True,"kh":True},"debug":False,\
#                   "3point":True,"edgeguard":True,"H_0":1000}
#         out= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)

#         inv = out["surfaces"]
#         inv = nstools.streamFuncToUV(inv,neighbors,distances)
#         # graph.graphVectorField(brasil,inv,"uabs","vabs","pres",# \
#         #                         metadata=out["metadata"],select=[3600],\
#         #                         transform=False,show=False,scale=0.1,\
#         #                         savepath="../articcirc-pics/vectorfields/curls/"+str(k))


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



