import nstools
import inverttools
import interptools
import bathtools
import saloffset
import graph
import parametertools as ptools
import matplotlib.pyplot as plt
import numpy as np

def mConditionVsError(inverse,surfaces,neighbors,lookups):
    conditions = []
    levels = []
    errors = []
    for lowlevel in range(1000,3600,200):
        inv,columndict,svds,A,e= inverttools.invert(inverse,surfaces,neighbors,lookups,lowlevel=lowlevel,highlevel=1000)
        s = np.diag((1.0/svds[1]))
        condition = s[0]/s[-1]
        conditions.append(condition)
        errors.append(np.sum(np.abs(e[1])))
        levels.append(lowlevel)
        #if lowlevel ==1800:
            #coupleinvert = nstools.streamFuncToUV(inv,neighbors,lookups)
            #coupleinvert = bathtools.addBathToSurface(inv)
            #graph.transportLine(coupleinvert,(-24.71,83.60),(27.496,80.06),2400,False)
            #graph.transportLine(coupleinvert,(-24.71,83.60),(27.496,80.06),2400,True)
            #graph.graphVectorField(inv,"uabs","vabs","z")
    fig, ax1 = plt.subplots()
    ax1.scatter(levels,errors,label="error")
    ax1.set_ylabel("R")
    ax2 = ax1.twinx() 
    ax2.scatter(levels,conditions,color="red",label = "matrix condition")
    ax2.set_ylabel("Matrix Condition")
    ax1.set_xlabel("Lowest Neutral Surface")
    plt.legend()
    plt.show()
