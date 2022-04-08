from regionlib import brasil
import os, gsw
import pickle
import graph, nstools, interptools, inverttools
from functools import partial
import parametertools as ptools
import numpy as np
import matplotlib.pyplot as plt
import pdb
from scipy.io import savemat

switch = { "reflevel":False,
           "H_0":False,
           "no mix":False,
           "bounds":False,
           "column weighting":False,
           "gaussian noise": False,
           "reference station":True
}
if switch["reflevel"]:
    ################################
    ### REFERENCE LEVEL TESTING 
    ################################
    with open('data/withparams.pickle', 'rb') as outfile:
        [surfaces,neighbors,distances] = pickle.load(outfile)
    reflevels = surfaces.keys()
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



if switch["no mix"]:
    ######################
    ###### No Mixing
    ######################

    with open('data/run0/withparams.pickle', 'rb') as outfile:
        [surfaces,neighbors,distances] = pickle.load(outfile)

    params = {"reflevel":2000,"upperbound":1000,"lowerbound":4000,\
            "mixs":{"kvo":False,"kvb":False,"kh":False},"debug":False,\
                "3point":True,"edgeguard":True,"H_0":1000}

    out= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)

    with open('sens/nomix.pickle'.format(k), 'wb') as outfile:
        pickle.dump([out,neighbors,distances],outfile)


if switch["bounds"]:
    ######################
    ###### Change upper and lower bound of inversion
    ######################
    for delta in range(1,4):
        newupperbound = 1000+delta*200
        newlowerbound = 4000-delta*200

        with open('data/run0/withparams.pickle', 'rb') as outfile:
            [surfaces,neighbors,distances] = pickle.load(outfile)

        params = {"reflevel":2000,"upperbound":newupperbound,"lowerbound":newlowerbound,\
                "mixs":{"kvo":True,"kvb":True,"kh":True},"debug":False,\
                    "3point":True,"edgeguard":True,"H_0":1000}

        out= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)

        with open('sens/bounds{}-{}.pickle'.format(newlowerbound,newupperbound), 'wb') as outfile:
            pickle.dump([out,neighbors,distances],outfile)

if switch["column weighting"]:
    ######################
    ###### Change column weighting by a factor of magnitude
    ######################
    for delta in ((1,0,0),(0,1,0),(0,0,1),(-1,0,0),(0,-1,0),(0,0,-1)):
        with open('data/run0/withparams.pickle', 'rb') as outfile:
            [surfaces,neighbors,distances] = pickle.load(outfile)

        params = {"reflevel":2000,"upperbound":1000,"lowerbound":4000,\
                "mixs":{"kvo":True,"kvb":True,"kh":True},"debug":False,\
                    "3point":True,"edgeguard":True,"H_0":1000,
                  "scalecoeffs":{"Ar":0.05,"kvo":5*10**(-6+delta[0]),"kvb":5*10**(-4+delta[1]),"kh":5*10**(2+delta[2])}}

        out= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)

        with open('sens/columnweighting-{}-{}-{}.pickle'.format(delta[0],delta[1],delta[2]), 'wb') as outfile:
            pickle.dump([out,neighbors,distances],outfile)


if switch["gaussian noise"]:
    ######################
    ###### Add Gaussian noise to neutral surface determination
    ######################

    for noise in range(1,4):
        with open('data/run0/annotatedbrasilargowoce.pickle', 'rb') as outfile:
            preinterpsurfaces,profiles = pickle.load(outfile)

        preinterpsurfaces = nstools.addDataToSurfaces(brasil,profiles,preinterpsurfaces,noise=noise*15)

        surfaces,neighbors,distances = \
            interptools.interpolateSurfaces(brasil,preinterpsurfaces,\
                                            interpmethod="gam",smart=False,coord="latlon")

        surfaces, neighbors, distances = nstools.addParametersToSurfaces(brasil,surfaces,neighbors,distances,H_0=1000)

        params = {"reflevel":2000,"upperbound":1000,"lowerbound":4000,\
                "mixs":{"kvo":True,"kvb":True,"kh":True},"debug":False,\
                    "3point":True,"edgeguard":True,"H_0":1000
                }

        out= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)

        with open('sens/noise-{}.pickle'.format(noise*30), 'wb') as outfile:
            pickle.dump([out,neighbors,distances], outfile)

if switch["reference station"]:
    ######################
    ###### Choosing a different reference station
    ######################
    with open('data/run0/argoandwoce.pickle', 'rb') as outfile:
        profiles = pickle.load(outfile)


    ## profilechoice = nstools.profileInBox(profiles,-40,-20,-31,-28,5000)
    ## #profilechoice = nstools.profileInBox(profiles,-45,-20,-40,-31,5000)
    profilechoice = nstools.profileInBox(profiles,-25,-20,-25,-20,4500)
    profilechoice = profilechoice[0]

    preinterpsurfaces = nstools.runPeerSearch(profiles,range(200,5000,200),profilechoice,False,10**10)

    preinterpsurfaces = nstools.addDataToSurfaces(brasil,profiles,preinterpsurfaces,noise=0)

    surfaces,neighbors,distances = \
        interptools.interpolateSurfaces(brasil,preinterpsurfaces,\
                                        interpmethod="gam",smart=False,coord="latlon")

    surfaces, neighbors, distances = nstools.addParametersToSurfaces(brasil,surfaces,neighbors,distances,H_0=1000)

    params = {"reflevel":2000,"upperbound":1000,"lowerbound":4000,\
            "mixs":{"kvo":True,"kvb":True,"kh":True},"debug":False,\
                "3point":True,"edgeguard":True,"H_0":1000
            }

    out= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)

    with open('data/sens/differentref-{}-{}.pickle'.format(profilechoice.lon,profilechoice.lat), 'wb') as outfile:
        pickle.dump([out,neighbors,distances], outfile)



