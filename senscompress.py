import pickle
# Extremely simple utility to pair down our model output to what we care about
def senscompress(surfaces,outpath):
    smallsurfaces = {}
    for k in surfaces.keys():
        smallsurfaces[k] = {}
        for j in surfaces[k].keys():
            if j != "data":
                smallsurfaces[k][j]=surfaces[k][j]
        smallsurfaces[k]["data"]={}
        for j in ["t","s","pres","pv","uabs","vabs","u","v","FQ","FS","kvb","kvo",'kvberror','kvoerror', 'kverror', 'usol', 'vsol', 'uerror', 'verror']:
            smallsurfaces[k]["data"][j] = surfaces[k]["data"][j]
    with open(outpath, 'wb') as outfile:
        pickle.dump(surfaces, outfile)

