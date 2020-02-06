from regionlib import brasil
import os
import graph
import nstools
import pickle
import regions
import interptools
import inverttools
import numpy as np


profiles = brasil.extractProfiles(os.path.abspath("data/brasilnc/"))
profiles = nstools.cruiseSearch(profiles,"A16S")
profiles = nstools.timeFilter(profiles,1.3*10**7,10**6)
ls = []
for p in profiles:
    ls.append(p.lat)
profilechoice = profiles[np.argmax(ls)]
preinterpsurfaces = nstools.runPeerSearch(profiles,range(100,6000,200),profilechoice,False,10**10)

with open('data/nsexperimentpicle.pickle', 'wb') as outfile:
    pickle.dump([preinterpsurfaces,profiles],outfile)
with open('data/nsexperimentpicle.pickle', 'rb') as outfile:
    preinterpsurfaces,profiles = pickle.load(outfile)
#nstools.stationDeltas(profiles)

#for c in nstools.cruiseCount(profiles):
    #graph.plotCruise(profiles,c,region="brasil")
#graph.plotCruise(profiles,"A16S",region="brasil")
a16s = nstools.cruiseSearch(profiles,"A16S")
graph.threedTransect(a16s,range(100,6000,200))


