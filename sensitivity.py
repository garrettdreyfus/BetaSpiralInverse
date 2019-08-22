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

def conditionError(inverse,surfaces,neighbors,lookups,disp=-1,reflevel=1000,savepath=False):
    conditions = []
    levels = []
    errors = []
    fram = []

    if savepath:
        try:
            os.makedirs(savepath)
        except FileExistsError as e:
            print(e)


    for lowlevel in range(reflevel+200,3600,200):
        inv,columndict,svds,A,e= inverttools.invert(inverse,surfaces,neighbors,lookups,lowlevel=lowlevel,highlevel=reflevel-400,reflevel=reflevel)
        if inv:
            s = np.diag((1.0/svds[1]))
            condition = s[0]/s[-1]
            conditions.append(condition)
            errors.append(np.sum(np.abs(e[1])))
            levels.append(lowlevel)
            i,o = graph.transportLine(inv,(-24.71,83.60),(27.496,80.06),2400,False,show=False)
            fram.append(abs(abs(i)-abs(o))) 
            if lowlevel ==disp:
                coupleinvert = nstools.streamFuncToUV(inv,neighbors,lookups)
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
    fig.suptitle(reflevel)
    ax3.scatter(levels,fram)
    ax3.set_ylabel("Net transport through Fram Strait in Sv")
    ax3.set_xlabel("Lowest Neutral Surface")
    plt.legend()
    fig.set_size_inches(16.5,12)
    if savepath:
        plt.savefig(savepath+"ns"+str(reflevel)+".png")
    plt.show()

def conditionErrorRefLevel(inverse,surfaces,neighbors,lookups,disp=-1,savepath=False):
    for reflevel in range(400,1600,200):
        conditionError(inverse,surfaces,neighbors,lookups,disp=disp,reflevel=reflevel,savepath=savepath)

