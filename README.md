# Beta Spiral Inverse Toolbox

Hello and welcome to the beta spiral inverse toolbox. The purpose of this library is to provide a set of tools to allow people to conduct beta spiral inverses, in the fashion of "The abyssal and deep circulation of the Northeast Pacific Basin" (Hauta 2018)

## First steps - creating profiles!
This toolbox takes in data in the form of Profile objects. To construct a Profile object requires a set of potential temperatures, salinity readings, and pressures. Each profile should also be given a unique id, and can be provided with a cruise name.

Below is an example construction of a profile object
```
from betaspiral.profile import Profile
id = 1
data ={"lat":80.0,"lon":-120.0,pres":[10,100,1000], "temp":[4,3,2], "sal":[35.0,34.5,34.5]}
newprofile = Profile(id,data)
```

The hope is that if you are planning on conducting an inverse you have alot of data that you can turn into profile objects and throw into one big list.


## salinity offsets and quality control
## Finding Neutral Surfaces!
So you have a big bunch of profiles and now you want to find the neutral surfaces. LETS GET IT!

The traditional way to mark neutral surfaces on profiles is to find a depth of neutral density relative to a premarked station closest to you (MCDOUGAL PAPER). However in some places (like the arctic where this was made for) there is are not many readily available stations. Due to this, this library takes a slightly different approach to finding neutral surfaces.It begins with one profile that can be thought of as the reference profile (but really can be any one you choose). Arbitrary neutral surfaces are marked down the depth of the profile. At every depth, profiles are considered within a given range and then the neutral depth for that surface is labelled on them. The program then considers all profiles with a labelled neutral depth on this neutral surface, reference profiles from which it searches for more profiles. In doing so it crawls out and finds neutral surfaces across a basin.

```
from betaspiral import nstools

#let profiles be a large list of profile objects
#we will arbitrarily choose the deepest profile we have as the reference profile
deepestindex=nstools.deepestProfile(profiles)
profilechoice = profiles[deepestindex]
#kthis command searchs for surfaces from 200m to 4000m by steps of 200 meters, and only finds neutral depths between profiles which are within 100 km of each other.
surfaces = nstools.runPeerSearch(profiles,200,4000,200,profilechoice,100)
```
## Interpolating
We now have a surfaces dictionary and it seems like a good time to describe what that entails. The surfaces object is a dictionary with keys of depth, and values of surface dictionaries. Surface dictionaries have arrays of lat, lon, x, y, and id. X,y is a cartesian projection of lat,lon for the purpose of interpolation and id is what ties together points on different neutral surfaces. The surface dictionary also has a "data" key which stores a dictionary containing yet more arrays of all the measurements and values we have on the surface. This contains but is not limited to temperature,salinity, potential vorticity, height of neutral surface, gradients

Interpolating is quite easy. There are two options when interpolating, Gam interpolation, or linear interpolation. Loosely speaking if you are dealing with low nosie data like gridded data from a model then linear interpolation is a good choice, and if you are dealing with noisy observational data than Gam is for you.

```
surfaces,neighbors,distances = interptools.interpolateSurfaces(surfaces,gaminterpolate=True)

```
Neighbors is a dictionary which contains dictionaries for each neutral surface on which points are adjacent to which on the grid. Distances  is a dictionary which contains dictionaries for each neutral surface on the **real** distances between adjacent points on the grid.

## Adding parameters
Alot of things go into an inverse and they need to be added to your surfaces
```

surfaces = nstools.addHeight(surfaces)
surfaces = nstools.addHorizontalGrad(surfaces,neighbors,distances)
surfaces = nstools.addBathAndMask(surfaces,neighbors)
surfaces = nstools.addVerticalGrad(surfaces)
ptools.saveBathVarTermCache(surfaces,"data/bathVarecco.pickle")
surfaces = nstools.addK(surfaces,"data/bathVarecco.pickle")

```

## Performing the actual inverse!

## Sensitivity Stuff

## Graphing and Diagnosing

## Example inverse run
```

deepestindex=nstools.deepestProfile(profiles)
profilechoice = deepestindex

surfaces = nstools.runPeerSearch(profiles,200,4000,200,profilechoice,1000)
surfaces = nstools.addDataToSurfaces(profiles,surfaces,2)

nstools.surfaceDiagnostic(surfaces)

surfaces,neighbors,distances = interptools.interpolateSurfaces(surfaces,gaminterpolate=False)
surfaces = nstools.addParametersToSurfaces(surfaces,neighbors,distances)
params = {"reflevel":1000,"upperbound":1000,"lowerbound":1400,"mixs":{"kvo":False,"kvb":False,"kh":False},\
        "debug":True,"modelmixing":True}
out = inverttools.invert("coupled",surfaces,neighbors,distances,params=params)
inv = out["surfaces"]
inv = nstools.streamFuncToUV(inv,neighbors,distances)
graph.graphVectorField(inv,"uabs","vabs","z")

```

## Useful things that haven't been mentioned
