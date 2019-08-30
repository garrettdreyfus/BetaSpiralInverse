import nstools
import inverttools
import interptools
import bathtools
import saloffset
import graph
import parametertools as ptools
import matplotlib.pyplot as plt
import numpy as np
import os

def distanceFromKnown(surfaces):
    s =  []
    for k in surfaces.keys():
        for j in range(len(surfaces[k]["data"]["knownu"])):
            if ~np.isnan(surfaces[k]["data"]["knownu"][j]) and ~np.isnan(surfaces[k]["data"]["uabs"][j]) :
                s.append(abs(surfaces[k]["data"]["knownu"][j]-surfaces[k]["data"]["uabs"][j]))
                s.append(abs(surfaces[k]["data"]["knownv"][j]-surfaces[k]["data"]["vabs"][j]))
    return np.sum(s)
def conditionError(inverse,surfaces,neighbors,distances,fname=False,disp=-1,params={},savepath=False,title=False,show=True):
    conditions = []
    levels = []
    errors = []
    fram = []
    print(params)
    params.setdefault("upperbound",1000)
    params.setdefault("reflevel",1000)
    if savepath:
        try:
            os.makedirs(savepath)
        except FileExistsError as e:
            print(e)
    for lowlevel in range(params["upperbound"]+1000,2600,200):
        params.update({"lowerbound":lowlevel})
        inv,columndict,svds,A,e,meta= inverttools.invert(inverse,surfaces,neighbors,distances,params=params)
        if inv:
            s = np.diag((1.0/svds[1]))
            condition = s[0]/s[-1]
            conditions.append(condition)
            errors.append(np.sum(np.abs(e[1])))
            levels.append(lowlevel)
            i,o = graph.transportLine(inv,(-24.71,83.60),(27.496,80.06),2400,False,show=False)
            #fram.append(abs(abs(i)-abs(o))) 
            inv = nstools.streamFuncToUV(inv,neighbors,distances)
            fram.append(abs(distanceFromKnown(inv))) 
            if lowlevel ==disp:
                coupleinvert = nstools.streamFuncToUV(inv,neighbors,distances)
                coupleinvert = bathtools.addBathToSurface(inv)
                #graph.transportLine(coupleinvert,(-24.71,83.60),(27.496,80.06),2400,True)
                graph.graphVectorField(inv,"uabs","vabs","z")
    fig, (ax1,ax3) = plt.subplots(2,1)
    ax1.scatter(levels,errors,label="error")
    ax1.set_ylabel("R")
    ax2 = ax1.twinx() 
    ax2.scatter(levels,conditions,color="red",label = "matrix condition")
    ax2.set_ylabel("Matrix Condition")
    ax1.set_xlabel("Lowest Neutral Surface")
    if not title:
        fig.suptitle(params["reflevel"])
    else:
        fig.suptitle(title)
    ax3.scatter(levels,fram)
    ax3.set_ylabel("Net transport through Fram Strait in Sv")
    ax3.set_xlabel("Lowest Neutral Surface")
    plt.legend()
    fig.set_size_inches(16.5,12)
    if savepath and not fname:
        plt.savefig(savepath+"ns"+str(params["reflevel"])+".png")
    if savepath and fname:
        plt.savefig(savepath+str(fname)+".png")
    if show:
        plt.show()
    plt.close()

def conditionErrorRefLevel(inverse,surfaces,neighbors,distances,disp=-1,savepath=False,params={}):
    for reflevel in range(600,1600,200):
        print(params)
        params.update({"mixs":[True,False,True],"reflevel":reflevel,"upperbound":reflevel})
        print(params)
        conditionError(inverse,surfaces,neighbors,distances,disp=disp,params=params,savepath=savepath)

def mixSens(inverse,surfaces,neighbors,distances,disp=-1,savepath=False,params={}):
    params = {"reflevel":600,"upperbound":600,"lowerbound":1600,"mixs":[True,True,True]}
    for i in np.arange(-5,-4,0.1):
        for j,name in enumerate(["Kvo","Kvb","Kvh"]):
            if  j== 1:
                newmix = mixs.copy()
                newmix[j] = newmix[j]*(10**i)
                params["mixcoeffs"] = newmix
                conditionError(inverse,surfaces,neighbors,distances,fname=name+str(i),disp=disp,params=params,savepath=savepath,title=name+": "+str(i),show=False)
