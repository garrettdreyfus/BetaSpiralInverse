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

###############################
#### GAM VS OBSERVATION COMPARE
###############################
switch = {"fig1":False,"fig2":False,"fig3":False,"fig4":True,"fig5":False,}

if switch["fig2"]:
    with open('data/run0/annotatedbrasilargowoce.pickle', 'rb') as outfile:
        preinterpsurfaces,profiles = pickle.load(outfile)
    graph.NSGAMCompareCruise(preinterpsurfaces,"A10",brasil)
    plt.savefig("fig1.eps")


#with open('data/brasilsurfaceswdata.pickle', 'rb') as outfile:
     #preinterpsurfaces,profiles = pickle.load(outfile)
# with open('data/invertedbrasil.pickle', 'rb') as outfile:
#      [inv,neighbors,distances] = pickle.load(outfile)
# inv = nstools.twoCReference(inv)
# inv = nstools.streamFuncToUV(inv,neighbors,distances)
# inv = nstools.externalReference(inv,"data/yomaha_1000.nc")
# print(nstools.transportDiagnostics(inv))
#inv = nstools.addGammaN(inv)
#graph.surfaceGammaRanges(inv)
#graph.graphSurfaces(brasil,inv,"kvb",contour=True,secondsurface=preinterpsurfaces,select=range(4000,10000),show=False,savepath="../arcticcirc-pics/surfaces/finalgam/")
#graph.graphSurfaces(brasil,inv,"kvb",contour=False,select=range(4000,10000),log=True)
#graph.graphSurfaces(brasil,inv,"z",contour=False)


if switch["fig3"]:
    # ###############################
    # #### Bubble plot comparisons
    # ###############################
    with open('data/run0/invertedbrasil.pickle', 'rb') as outfile:
        [inv,neighbors,distances] = pickle.load(outfile)
    inv = inv["surfaces"]

    kvs = []
    for l in inv.keys():
        kvs += list(inv[l]["data"]["kvb"])
    print("kv: ,",np.sum(((np.asarray(kvs).flatten()<0)))/len(np.asarray(kvs).flatten()))

    inv = nstools.twoCReference(inv)
    inv = nstools.streamFuncToUV(inv,neighbors,distances)
    inv = nstools.externalReference(inv,"data/yomaha_1000.nc")
    graph.bubblePlot(inv,-35,-30)
    plt.savefig("fig2.eps")
    #therms = list(np.linspace(0,2,11))[::-1]
    #transports = []
    #transports_2 = []
    #for l in range(10):
        #transports.append(graph.transportTempClass(inv,-29.5,-30,-9,1000,6000,therms[l],therms[l+1],vfield="vabs")/1000000)
        #transports_2.append(graph.transportTempClass(inv,-29.5,-30,-9,1000,6000,therms[l],therms[l+1],vfield="v")/1000000)
    #print(transports)
    #therms = np.asarray(therms)
    #plt.barh((therms[1:] + therms[:-1]) / 2,transports,height=0.15,label="Run 0" )
    #plt.barh((therms[1:] + therms[:-1]) / 2,transports_2,height=0.15,label="2C reference" )
    #plt.gca().set_xlim(-0.2,0.6)
    #plt.gca().set_xlabel("Transport (Sv)")
    #plt.gca().set_ylabel("Potential Temperature Threshold")
    #plt.legend()
    #plt.show()

    # graph.graphVectorField(brasil,inv,"2CU","2CV","pres",\
    #                         transform=False, scale=0.1,select=[1200])
    # graph.graphSurfaces(brasil,inv,"2CU")


    plt.show()

if switch ["fig4"]:

    with open('data/invertedbrasil.pickle', 'rb') as outfile:
        [inv,neighbors,distances] = pickle.load(outfile)
    inv = nstools.twoCReference(inv)
    inv = nstools.streamFuncToUV(inv,neighbors,distances)
    inv = nstools.externalReference(inv,"data/yomaha_1000.nc")
    graph.fourpanelVectorField(brasil,inv,"uabs","vabs",backgroundfield="pres",\
                            select=[1000,1800,3400,4400],transform=False,scale=0.1)
    plt.savefig("fig3.eps")
    # print(inv.keys())
###############################
#### Figure 4 - Our AABW schematic
###############################
with open('data/invertedbrasil.pickle', 'rb') as outfile:
    [inv,neighbors,distances] = pickle.load(outfile)

inv = nstools.addOldUnits(inv)
inv = nstools.streamFuncToUV(inv,neighbors,distances)
graph.average_velocity_temp_class(brasil,inv,-30,-35)
# graph.graphVectorField(brasil,inv,"uabs","vabs","pv",\
#                         transform=False, scale=0.1,select=[3600])
# inv = nstools.addOldUnits(inv)
# graph.sigma4Plot(inv,-29.5,-180,-33)
# inv = nstools.twoCReference(inv)
# inv = nstools.streamFuncToUV(inv,neighbors,distances)
# inv = nstools.externalReference(inv,"data/yomaha_1000.nc")
# graph.transportRefIsotherm(inv,2,-29.5,-30,-9,1000,6000)
# print(inv.keys())
# graph.graphSurfaces(brasil,inv,"pv",contour=True,show=False, select=range(3600,3601),secondsurface=inv)
# fig = plt.gcf()
# left, bottom, width, height = [0.25, 0.6, 0.2, 0.2]
# ax2 = fig.add_axes([left, bottom, width, height])
# ax2.set_xlabel("Potential Temperature")
# ax2.set_ylabel("Practical Salinity")
# c = ax2.scatter(inv[3600]["data"]["pottemp"],inv[3600]["data"]["psu"],c=gsw.rho(inv[3600]["data"]["s"],inv[3600]["data"]["t"],2000))
# cbar = plt.colorbar(c,ax=ax2)
# cbar.ax.set_ylabel("$\sigma^2$")
# plt.show()


