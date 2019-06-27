import nstools

profiles,deepestindex = nstools.extractProfilesBox(["data/3000m2007profiles.json"],-100,100,80,90)

#nstools.plotCruise(profiles,"IPY 2007")

surfaces = nstools.search(profiles,deepestindex)
nstools.graphSurfaces(surfaces)
