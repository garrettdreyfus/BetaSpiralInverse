# Beta Spiral Inverse Toolbox

Hello and welcome to the beta spiral inverse toolbox. The purpose of this library is to provide a set of tools to allow people to conduct beta spiral inverses, in the fashion of "The abyssal and deep circulation of the Northeast Pacific Basin" (Hauta 2018). If you have any questions or need any help please contact me at ( garretf at uw dot edu ). I would love to help anyone get an inverse up and running!

## Creating a Region
The first step to using this toolkit is to create a region object and add it to the region library. This step can be skipped if you are trying to invert a region which has already been inverted and plan to use the same data as whoever created that region file. If you do go through creating a region file it would be very nice of you to submit a pull request to this library so we can build an awesome library of regoins!


### Structure
For the purpose of this tutorial let us assume the region we want to add is the waters around Tahiti! First, we should create the file "tahiti.py" in the regionlib folder. This file should contain the following functions
- [ ] createMesh
- [ ] geoFilter
- [ ] some function which creates profile objects from your data
And a dictionary for graphing purposes:
- [ ] mapbounds

###createMesh
This function is what tells the toolbox the grid over which you would like to interpolate and then invert. It should take 4 parameters,the xvalues of the data, the y values of the data, the coordinate system, and a spacing scale. These four parameters are provided because they can be useful but none of them are required to actually make the function operate correctly. It should return a numpy meshgrid with XY coordinates. Below I will ilustrate the simple example of a 10 degree square around Tahiti.Tahiti is located at roughly 18 S, 149 W. Later steps in the toolbox will remove any points that intersect with land or the sea bed.
```
##this is a useful function to convert lat lon to x y
def singleXY(coord):
    theta = np.deg2rad(coord[0])
    r = ((90-coord[1]) *111*1000)
    x = (r*np.cos(theta))
    y = (r*np.sin(theta))
    return x,y

def createMesh(xvals,yvals,coord="xy",spacingscale=25):
	grd = np.meshgrid(np.linspace(-154,-144,20*spacingscale),np.linspace(-23,-13,20*spacingscale))
     for i in range(grd[0].shape[0]):
		for j in range(grd[0].shape[1]):
            x,y = singleXY((grd[0][i][j],grd[1][i][j]))
            grd[0][i][j] = x
            grd[1][i][j] = y
	return grd


```
Phew! Glad that's over! createMesh is the most confusing one! 


### geoFilter

geoFilter should take a latitude and longitude and return True if the point is in the domain and False if it is not. This function is used to limit the data we consider for interpolation, but also limit the points that we will eventually put into our inverse. It may seem confusing why we would need this function if we have already specified our grid in createMesh, but the fact of the matter is its easier to lay down a square and cut out what you dont want than to specify the perfect grid. If for example we wanted to lay a latlon grid over the Bering Sea we would find that part of it would peak out south over the Aleutian Islands. We could then crop these points out with geoFilter. However, for Tahiti no cropping is needed.

```
def geoFilter(lon,lat):
	return True


```

### The function that extracts profiles from your data.

This function is by far the loosest of the three. It is not used internally by the toolbox so it can have any number of parameters. The end goal is to produce a list of profiles. Well how do you create a profile you ask?  To construct a Profile object requires a set of potential temperatures, salinity readings, and pressures. Each profile should also be given a unique id, and can be provided with a cruise name, and time. The profile object assumes the data given is in practical salinity and in-situ temperature,and converts it to absolute salinity and conservative temperature. You can specify that the data has already been converted by setting ct and abssal to True instead of False as is done below.

Here is how you can construct a profile object
```
##if outside project directory
from betaspiral.profile import Profile
##if inside project directory
from profile import Profile

def extractProfiles():
	profiles = []
    id = 1
    data ={"lat":-20.0,"lon":-150.0,pres":[10,100,1000], "temp":[4,3,2], "sal":[35.0,34.5,34.5], "cruise":"CRUISE1", "time":1234122}
    profiles.append(Profile(id,data,ct=False,abssal=False))
    id = 2
    data ={"lat":-21.0,"lon":-151.0,pres":[10,110,1000], "temp":[5,3,2], "sal":[35.0,34.5,39.5], "cruise":"CRUISE1", "time":1234122}
    profiles.append(Profile(id,data,ct=False,abssal=False))
    return profiles
```

Any real example would hopefully include no hard-coded data and would rather pull from a large data set to generate a list of hundreds if not thousands of profiles.

### mapbounds
I use [Basemap](http://https://matplotlib.org/basemap/) to do all the graphing parts of this library. The dictionary contains 6 values, 4 describe the corners of the domain over which we'd like to map, and 2 which describe the center. Here is a simple example for the Brasil Basin. "llon" is "lower left longitude" and "urlat" is "upper right lat". lat_0 and lon_0 are the center of the domain.

```

mapbounds = {"lllon":-85,"urlon":16,"lllat":-50,"urlat":8,\
        "lat_0":-17,"lon_0":-33}

```

## Finding Neutral Surfaces!
So you have specified a region file, have a function to generate a big bunch of profiles and now you want to find the neutral surfaces. LETS DO IT!

Neutral Surfaces are found by comparing profiles against a reference profile. There are multiple ways to choose this reference profile, however, any reference profile you choose should be at least as deep as the deepest neutral surface which you hope to find. You could choose a profile that you trust to have accurate data, you could choose a profile within a given region, or even just a random profile. 

```
import nstools
## this won't actually run because our regionlib/tahiti.py does not exsist (yet!)
from regionlib import tahiti

profiles = tahiti.extractProfiles()

#we will arbitrarily choose the deepest profile we have as the reference profile
deepestindex=nstools.deepestProfile(profiles)
profilechoice = profiles[deepestindex]

#alternatively we could choose a profile between -150W - 148W, and 20S - 22S 
#that is at least 6000 dbar deep
profilechoice = nstools.profileInBox(profiles,-150,-148,-22,-20,6000)
profilechoice = profilechoice[0]

#this command searchs for surfaces from 200m to 4000m by steps of 200 meters.
surfaces = nstools.runPeerSearch(profiles,range(100,6000,200),profilechoice)
#Although we could also choose a completely arbitrary set of neutral surfaces
surfaces = nstools.runPeerSearch(profiles,[200,780,890,1021],profilechoice)

surfaces = nstools.addDataToSurfaces(tahiti,profiles,surfaces)


```
## Interpolating
We now have a surfaces dictionary and it seems like a good time to describe what that entails. The surfaces object is a dictionary with keys of depth, and values of surface dictionaries. Surface dictionaries have arrays of lat, lon, x, y, and id. X,y is a cartesian projection of lat,lon for the purpose of interpolation and id is what ties together points on different neutral surfaces. The surface dictionary also has a "data" key which stores a dictionary containing yet more arrays of all the measurements and values we have on the surface. This contains but is not limited to temperature,salinity, potential vorticity, height of neutral surface, gradients

Interpolating is quite easy. There are two options when interpolating, Gam interpolation, or linear interpolation. Loosely speaking if you are dealing with low nosie data like gridded data from a model then linear interpolation is a good choice, and if you are dealing with noisy observational data than the GAM is for you. (It is also worth noting that the choice of parameters to the GAM is not a trivial one. The library uses the library [pyGam](http://https://github.com/dswah/pyGAM) default number of splines because it did a fine job in the regions we considered)

```
import interptools
from regionlib import tahiti

surfaces,neighbors,distances = interptools.interpolateSurfaces(tahiti,surfaces,interpmethod="gam")

```
Neighbors is a set of dictionaries which describes which points on the grid are neighboring each other. Distances is a set of dictionaries which describes the distances between each neighboring point on the grid.

## Adding parameters
Alot of things go into an inverse and they need to be added to your surfaces
```
import nstools

##adds everything

surfaces = nstools.addParametersToSurfaces(tahiti,surfaces,neighbors,distances,[])

##nstools.inverseReady is a diagostic to show what fraction of your grid is Nan for every field and tell you if things are ready to be inverted. It is ok if NaN percentages are quite high. The process of masking out surfaces that intersect with land and the seabed can remove alot of points.

nstools.inverseReady(surfaces)
```
## Graphing and Diagnosing
We haven't gotten to invert yet but I figured I'd describe some of the functionality of the graphing part of this toolbox in case you want to see whats happened so far. All of the graphing functionality lies in graph.py. This file has a ton of functions that all make different kinds of graphs but is very cluttered. I will describe the three functions I think are most likely going to be useful to you the end user, graphSurfaces, saveAllQuants, and graphVectorField.

### graphSurfaces

This function creates graphs of a given quantity at all neutral surfaces in a surfaces object. There are alot of parameters. Lets go through them!



## Performing the actual inverse!
AMAZING!
## Sensitivity Stuff


## Example inverse run

## Useful things that haven't been mentioned
