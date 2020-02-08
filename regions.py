import numpy as np 
from interptools import singleXY
from netCDF4 import Dataset
import numpy as np
import scipy
from regionlib import brasil as brasillib
from regionlib import arctic as arcticlib
from regionlib import nepb as nepblib

arctic = {"createMesh":arcticlib.createMesh,"geofilter":arcticlib.geoFilter,\
        "mapbounds":arcticlib.mapbounds,"name":"arctic"}
nepb = {"createMesh":nepblib.createMesh,"geofilter":nepblib.geoFilter,\
        "mapbounds":nepblib.mapbounds,"name":"nepb"}
brasil = {"geofilter":brasillib.geoFilter,"name":"brasil",\
        "mapbounds":brasillib.mapbounds,"createMesh":brasillib.createMesh}





