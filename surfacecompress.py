import pickle
import os
#Utility to pair down our model output to what we care about
def surfaceCompress(surfaces):
    smallsurfaces = {}
    for k in surfaces.keys():
        smallsurfaces[k] = {}
        for j in surfaces[k].keys():
            if j != "data":
                smallsurfaces[k][j]=surfaces[k][j]
        smallsurfaces[k]["data"]={}
        for j in ["t","s","pres","pv","uabs","vabs","u","v","FQ","FS","kvb","kvo",'kvberror','kvoerror', 'kverror', 'usol', 'vsol', 'uerror', 'verror','h','z']:
            smallsurfaces[k]["data"][j] = surfaces[k]["data"][j]
    return smallsurfaces

def sensitivityRunCompress(infilename,outfilename):
    with open(infilename, 'rb') as infile:
        [out, neighbors, distances] = pickle.load(infile)
    surfaces = surfaceCompress(out["surfaces"])
    errors = out["errors"]
    metadata = out["metadata"]
    with open(outfilename, 'wb') as outfile:
        pickle.dump([surfaces,errors,metadata], outfile)

def compressFolder(foldername):
    for f in glob.glob(ncfolder+"*.pickle"):
        name, ext = os.path.splitext(f)
        sensitivitiyRunCompress(f,"small-sens/"+name+"-small"+ext)


