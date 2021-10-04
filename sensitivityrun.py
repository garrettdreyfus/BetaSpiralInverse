from regionlib import brasil
import os, gsw
import pickle
import graph, nstools, interptools, inverttools
from functools import partial
import parametertools as ptools
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pdb
from scipy.io import savemat
import sensitivity

switch = { "reflevel":False,
           "H_0":True,
           "no mix":False
    }
if switch["reflevel"]:
    ################################
    ### REFERENCE LEVEL TESTING 
    ################################
    with open('data/withparams.pickle', 'rb') as outfile:
        [surfaces,neighbors,distances] = pickle.load(outfile)
    reflevels = surfaces.keys()
    levels = []
    hunters = []
    vemas = []
    conditions=[]
    curls = []
    errors = []
    for k in reflevels:
        if int(k) >1200 and int(k) < 4000:
            with open('data/run0/withparams.pickle', 'rb') as outfile:
                [surfaces,neighbors,distances] = pickle.load(outfile)
            params = {"reflevel":int(k),"upperbound":1000,"lowerbound":4000,\
                    "mixs":{"kvo":True,"kvb":True,"kh":True},"debug":False,\
                    "3point":True,"edgeguard":True,"H_0":1000}
            out= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)

            inv = out["surfaces"]
            inv = nstools.streamFuncToUV(inv,neighbors,distances)
            with open('sens/reflevel/{}.pickle'.format(k), 'wb') as outfile:
                pickle.dump([out,neighbors,distances],outfile)


if switch["H_0"]:
    ######################
    ###### H_0 Testing
    ######################

    for H_0 in [200,500,1000,1500,2000]:
        with open('data/run0/interpedbrasil.pickle', 'rb') as outfile:
            [surfaces,neighbors,distances] = pickle.load(outfile)
        surfaces, neighbors, distances = nstools.addParametersToSurfaces(brasil,surfaces,neighbors,distances,H_0=H_0)
        params = {"reflevel":2000,"upperbound":1000,"lowerbound":4000,\
                "mixs":{"kvo":True,"kvb":True,"kh":True},"debug":False,\
                    "3point":True,"edgeguard":True,"H_0":H_0
                }

        out= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)
        with open('sens/H_0/{}.pickle'.format(H_0), 'wb') as outfile:
            pickle.dump([out,neighbors,distances],outfile)



if switch["nomix"]:
    ######################
    ###### No Mixing
    ######################

    with open('data/run0/withparams.pickle', 'rb') as outfile:
        [surfaces,neighbors,distances] = pickle.load(outfile)
    params = {"reflevel":int(k),"upperbound":1000,"lowerbound":4000,\
            "mixs":{"kvo":True,"kvb":True,"kh":True},"debug":False,\
                "3point":True,"edgeguard":True,"H_0":1000}
    out= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)

    with open('sens/nomix.pickle'.format(k), 'wb') as outfile:
        pickle.dump([out,neighbors,distances],outfile)

