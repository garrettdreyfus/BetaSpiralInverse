import numpy as np 
from interptools import singleXY

def arcticMesh(n,xvals,yvals,coord="xy"):
    xmin= -1793163
    xmax = 971927
    ymin = -1455096
    ymax = 1200385
    return np.meshgrid(np.linspace(xmin,xmax,n), np.linspace(ymin,ymax,n),indexing="xy")

arctic = {"createMesh":arcticMesh}

def hautalaGrid():
    grd = np.meshgrid(\
            np.concatenate((np.linspace(170,178,5),np.linspace(-180,-122,30))),

            np.linspace(18,58,21))
    for i in range(grd[0].shape[0]):
        for j in range(grd[0].shape[1]):
            x,y = singleXY((grd[0][i][j],grd[1][i][j]))
            grd[0][i][j] = x
            grd[1][i][j] = y
    return grd


def nepbMesh(n,xvals,yvals,coord="xy"):
    if coord=="latlon":
        return hautalaGrid()
    elif coord == "xy":
        x1,y1 = singleXY((-144,3))
        x2,y2 = singleXY((155,63))
        xmin = min(x1,x2)
        xmax = max(x1,x2)
        ymin = min(y1,y2)
        ymax = max(y1,y2)
        return np.meshgrid(np.linspace(xmin,xmax,n), np.linspace(ymin,ymax,n),indexing="xy")
    else:
        print("This grid type is not supported", sys.exc_info()[0])


nepb = {"createMesh":nepbMesh}
