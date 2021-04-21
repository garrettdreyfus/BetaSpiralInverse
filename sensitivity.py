import nstools
import inverttools
import interptools
import bathtools
import saloffset
import graph
import parametertools as ptools
import matplotlib.pyplot as plt
import numpy as np
import os, copy
import pygam
import seaborn as sns
from regionlib import brasil

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
    for lowlevel in range(params["upperbound"]+200,2600,200):
        params.update({"lowerbound":lowlevel})
        print(params)
        out = inverttools.invert(inverse,surfaces,neighbors,distances,params=params)
        inv = out["surfaces"]
        columndict = out["coldict"]
        svds = out["svddecomp"]
        A = out["matrixsetup"]
        e = out["errors"]
        meta = out["metadata"]
        if inv:
            s = np.diag((1.0/svds[1]))
            condition = s[0]/s[-1]
            conditions.append(condition)
            errors.append(np.sum(np.abs(e[1])))
            levels.append(lowlevel)
            i,o = graph.transportLine(inv,(-24.71,83.60),(27.496,80.06),2400,False,show=False)
            fram.append(abs(abs(i)-abs(o))) 
            inv = nstools.streamFuncToUV(inv,neighbors,distances)
            #fram.append(abs(distanceFromKnown(inv))) 
            if lowlevel ==disp:
                coupleinvert = nstools.streamFuncToUV(inv,neighbors,distances)
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
        params.update({"mixs":{"kvo":True,"kvb":False,"kh":True},"reflevel":reflevel,"upperbound":reflevel})
        print(params)
        conditionError(inverse,surfaces,neighbors,distances,disp=disp,params=params,savepath=savepath)

def mixSens(inverse,surfaces,neighbors,distances,disp=-1,savepath=False,params={}):
    params = {"reflevel":600,"upperbound":600,"lowerbound":1600,"mixs":[True,True,True]}
    for i in np.arange(-5,-4,0.1):
        for j,name in enumerate(["Kvo","Kvb","Kh"]):
            if  j== 1:
                newmix = mixs.copy()
                newmix[j] = newmix[j]*(10**i)
                params["mixcoeffs"] = newmix
                conditionError(inverse,surfaces,neighbors,distances,fname=name+str(i),disp=disp,params=params,savepath=savepath,title=name+": "+str(i),show=False)

def gamSplineSens(preinterpsurfaces):
    sumerrors=[]
    for splines in range(4,16):
        error=[]
        for k in preinterpsurfaces.keys():
            if int(k)>2000:
                surface = preinterpsurfaces[k]
                X = np.zeros((len(surface["lons"]),2))
                X[:,0]=surface["lons"]
                X[:,1]=surface["lats"]
                #for d in Bar("Interpolating: ").iter(surface["data"].keys()):
                d = "pres"
                notnan = ~np.isnan(surface["data"][d])
                if np.count_nonzero(notnan)>10:
                    gam = pygam.GAM(pygam.te(0,1,n_splines=[splines,splines])).fit(X[notnan],np.asarray(surface["data"][d])[notnan])
                    #random_gam =  pygam.LinearGAM(pygam.s(0) + pygam.s(1) ).gridsearch(X, surface["data"][d])
                    error += list(np.log10(np.abs(surface["data"][d]-gam.predict(X))))
        sns.distplot(error,kde_kws={"fill":False,"label": str(splines)})
    #plt.plot(range(4,16),sumerrors)
    plt.legend()
    plt.show()


def cruiseSpline(preinterpsurfaces,cruise):
    sumerrors=[]
    splines=20
    error=[]
    for k in preinterpsurfaces.keys():
        surface = preinterpsurfaces[k]
        cruiseline = surface["data"]["cruise"] == cruise
        print(surface["data"]["cruise"])
        X = np.zeros((len(surface["lons"]),2))
        X[:,0]=surface["lons"]
        X[:,1]=surface["lats"]
        #for d in Bar("Interpolating: ").iter(surface["data"].keys()):
        d = "pres"
        notnan = ~np.isnan(surface["data"][d])
        if np.count_nonzero(notnan)>10:
            gam = pygam.GAM(pygam.te(0,1,n_splines=[splines,splines])).fit(X[notnan],np.asarray(surface["data"][d])[notnan])
            #random_gam =  pygam.LinearGAM(pygam.s(0) + pygam.s(1) ).gridsearch(X, surface["data"][d])
            error += list(np.log10(np.abs(surface["data"][d][cruiseline]-gam.predict(X[cruiseline]))))
    #sns.distplot(error,kde_kws={"fill":False,"label": str(splines)})
    #plt.plot(range(4,16),sumerrors)
    plt.legend()
    plt.show()

def decayScaleSensitivity(originalsurfaces,neighbors,distances):
    f = open("h0log.txt", "w")
    for H_0 in range(500,5500,750):
        surfaces = copy.deepcopy(originalsurfaces)
        surfaces = nstools.addParametersToSurfaces(brasil,surfaces,neighbors,distances,H_0=H_0)
        # graph.NSGAMCompare(preinterpsurfaces,surfaces,-30,-180,180)

        # print(surfaces.keys())
        params = {"reflevel":int(2062),"upperbound":1000,"lowerbound":4200,\
                "mixs":{"kvo":True,"kvb":True,"kh":True},"debug":False,\
                  "3point":True,"edgeguard":True,"H_0":H_0}
        # Conditions
        # All mixing: 201235
        # No mixing: 147
        # Kv0 only: 147
        # KvH and Kv0 only: 148
        # KvH and Kv0 only with out edgeguard (tm): 489

        out= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)
        inv = out["surfaces"]
        inv = nstools.streamFuncToUV(inv,neighbors,distances)
        print("#"*10+str(H_0),file=f)
        print(out["metadata"],file=f)
        result = nstools.transportDiagnostics(inv)
        print(result,file=f)
        kvs = []
        for l in inv.keys():
            kvs += list(inv[l]["data"]["kv"])
        kvs = np.asarray(kvs).flatten()
        print("kv: ,"+str(np.nanmean(kvs)),file=f)
        print("% negative: ,"+str(np.sum(kvs<0)/np.sum(~np.isnan(kvs))),file=f)
    f.close()
