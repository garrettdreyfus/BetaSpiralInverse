import pickle
import os
import glob
import nstools
#Utility to pair down our model output to what we care about
def surfaceCompress(surfaces):
    smallsurfaces = {}
    for k in surfaces.keys():
        smallsurfaces[k] = {}
        for j in surfaces[k].keys():
            if j != "data":
                smallsurfaces[k][j]=surfaces[k][j]
        smallsurfaces[k]["data"]={}
        for j in ["t","s","pres","pv","uabs","vabs","u","v","FQ","FS","kvb","kvo",'kvberror','kvoerror', 'kverror', 'usol', 'vsol', 'uerror', 'verror','h','z', \
'yomahau','yomahav','2CU','2CV','pottemp','psu']:
            smallsurfaces[k]["data"][j] = surfaces[k]["data"][j]
    return smallsurfaces

def sensitivityRunCompress(infilename,outfilename):
    with open(infilename, 'rb') as infile:
        [out, neighbors, distances] = pickle.load(infile)
    inv = out["surfaces"]
    inv = nstools.addOldUnits(inv)
    inv = nstools.twoCReference(inv)
    inv = nstools.streamFuncToUV(inv,neighbors,distances)
    inv = nstools.externalReference(inv,"data/yomaha_1000.nc")
    surfaces = surfaceCompress(inv)
    errors = out["errors"]
    metadata = out["metadata"]
    with open(outfilename, 'wb') as outfile:
        pickle.dump([surfaces,errors,metadata], outfile)

def compressFolder(infoldername,outfoldername):
    for f in glob.glob(infoldername+"*.pickle"):
        _, ext = os.path.splitext(f)
        name = os.path.basename(f)
        sensitivityRunCompress(f,outfoldername+name+"-small"+ext)

sensitivityRunCompress("sens/differentref--24-22.pickle","small-sens/differentref--24-22.pickle")
#compressFolder("sens/","small-sens/")
#compressFolder("sens/reflevel/","small-sens/reflevel/")
#compressFolder("sens/H_0/","small-sens/H_0/")
