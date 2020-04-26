from nstools import *
import pdb
import scipy
import copy
import parametertools as ptools

##calculate inverse from composite SVD matrixes
def SVDCalc(VT,D,U,n_elements=False):
    if n_elements:
        D = D[:, :n_elements]
        VT = VT[:n_elements, :]

    B = VT.T.dot(D.T).dot(U.T)
    return B

## generate psuedo inverse
## if full True return composite matrixes
def SVDdecomp(A,n_elements=2,full=True):
    U, s, VT = svd(A,full_matrices=True)
    # reciprocals of s
    d = 1.0 / s
    # create m x n D matrix
    D = np.zeros(A.shape)
    # populate D with n x n diagonal matrix
    D[:A.shape[1], :A.shape[1]] = np.diag(d)

    B = SVDCalc(VT,D,U,n_elements)
    if full:
        return B, [VT, D, U]
    else:
        return B

def generateMixingRows(columnindexs,mixing,k,n,i=0):
    Akvbrow = [0]*(columnindexs[i]+1)
    Akvbrow[columnindexs[i]]=-mixing["kvb"]/n
    Akhrow = [0]*(int(k/200))
    Akhrow[int(k/200-1)] = -mixing["kh"]/n
    Akvorow = [-mixing["kvo"]/n]
    return Akvbrow,Akhrow,Akvorow

def error(b,us,prime):
    b = np.asarray(b)
    us = np.asarray(us)
    prime = np.asarray(prime)
    errorafter = np.matmul(b,us+prime)
    errorbefore = np.matmul(b,us)
    return errorbefore, errorafter
 
#graph the error vector C before and after solution applied
def graphError(b,us,prime):
    b = np.asarray(b)
    us = np.asarray(us)
    prime = np.asarray(prime)
    errorbefore = np.matmul(b,us+prime/0.05)
    errorafter = np.matmul(b,us)
    plt.scatter(range(len(errorafter)),errorafter,label="after")
    plt.scatter(range(len(errorbefore)),errorbefore,label="before")
    plt.gca().legend()
    plt.show()
    delta = np.abs(errorbefore)-np.abs(errorafter)
    plt.scatter(range(len(delta[delta<0])),delta[delta<0],c="red")
    plt.scatter(range(len(delta[delta>=0])),delta[delta>=0],c="blue")
    plt.show()

##normalize a vector
def norm(v):
    return v/np.linalg.norm(v)



#get the column number of a parameter
#generates new number if necessary
def getColumnNumber(eyedict,eyed):
    #assign each eyed a column number 
    if "max" not in eyedict.keys():
        eyedict["max"]=-1
        eyedict[-999]=0
    if eyed not in eyedict.keys():
        eyedict[eyed]=eyedict["max"]+1 
        eyedict["max"] = eyedict["max"]+1
    return eyedict,eyedict[eyed]

#gets column numbers of point and neighbors
def neighborsColumnNumbers(surfaces,k,s,eyedict):
    ## return the column number of each neighbor
    columnnumbers = []
    for l in s:
        eyedict,col = getColumnNumber(eyedict,surfaces[k]["ids"][l])
        columnnumbers.append(col)
    return eyedict,columnnumbers

#subtract reference level from psi and velocities
def applyRefLevel(surfaces,reflevel=1000):
    for k in Bar("subtracting ref level: ").iter(surfaces.keys()):
        surfaces[k]["data"]["psiref"] = np.full(len(surfaces[k]["data"]["psi"]),np.nan)
        surfaces[k]["data"]["uref"] = np.full(len(surfaces[k]["data"]["psi"]),np.nan)
        surfaces[k]["data"]["vref"] = np.full(len(surfaces[k]["data"]["psi"]),np.nan)
        surfaces[k]["data"]["FQ"] = np.full(len(surfaces[k]["data"]["psi"]),np.nan)
        surfaces[k]["data"]["FS"] = np.full(len(surfaces[k]["data"]["psi"]),np.nan)
        for l in range(len(surfaces[k]["data"]["psi"])):
            found = np.where(np.asarray(surfaces[reflevel]["ids"])==surfaces[k]["ids"][l])
            if len(found) != 0 and len(found[0]) != 0:
                surfaces[k]["data"]["psiref"][l] = surfaces[k]["data"]["psi"][l] - surfaces[reflevel]["data"]["psi"][found[0][0]]
                surfaces[k]["data"]["uref"][l] = (surfaces[k]["data"]["u"][l]-surfaces[reflevel]["data"]["u"][found[0][0]])
                surfaces[k]["data"]["vref"][l] = (surfaces[k]["data"]["v"][l]-surfaces[reflevel]["data"]["v"][found[0][0]])
    return surfaces

def indexOfSurface(surfaces,params,k):
    cs=[]
    for i in surfaces.keys():
        if params["lowerbound"]>=i>=params["upperbound"]:
            cs.append(i)
    #print(cs)
    return sorted(cs).index(k)

## apply the solution to surfaces
def applyPrime(staggeredsurfaces,prime,coldict,params,widths,mixing=False):
    scales = params["scalecoeffs"]
    for k in Bar("adding solutions: ").iter(staggeredsurfaces.keys()):
        staggeredsurfaces[k]["data"]["psinew"] =  np.full(len(staggeredsurfaces[k]["lons"]),np.nan)
        staggeredsurfaces[k]["data"]["psisol"] =  np.full(len(staggeredsurfaces[k]["lons"]),np.nan)
        staggeredsurfaces[k]["data"]["kvb"] =  np.full(len(staggeredsurfaces[k]["lons"]),np.nan)
        staggeredsurfaces[k]["data"]["kvo"] =  np.full(len(staggeredsurfaces[k]["lons"]),np.nan)
        staggeredsurfaces[k]["data"]["kh"] =  np.full(len(staggeredsurfaces[k]["lons"]),np.nan)
        staggeredsurfaces[k]["data"]["usol"] =  np.full(len(staggeredsurfaces[k]["lons"]),np.nan)
        staggeredsurfaces[k]["data"]["vsol"] =  np.full(len(staggeredsurfaces[k]["lons"]),np.nan)
        #print(widths)
        #print(np.asarray(prime[widths[0]:widths[0]+widths[1]])*scales["kvb"])
        for i in range(len(staggeredsurfaces[k]["data"]["ids"])):
            eyed = staggeredsurfaces[k]["ids"][i] 
            if eyed in coldict.keys():
                staggeredsurfaces[k]["data"]["psinew"][i] = staggeredsurfaces[k]["data"]["psiref"][i] + prime[coldict[eyed]]*scales["Ar"]
                staggeredsurfaces[k]["data"]["psisol"][i] = prime[coldict[eyed]]*scales["Ar"]
                if params["mixs"]["kvb"] and params["mixs"]["kvo"] and params["mixs"]["kh"]:
                    staggeredsurfaces[k]["data"]["kvb"][i] = prime[widths[0]+coldict[eyed]]*scales["kvb"]
                    if params["lowerbound"]>=k>=params["upperbound"]:
                        #print(widths,len(prime))
                        #print(indexOfSurface(staggeredsurfaces,params,k))
                        staggeredsurfaces[k]["data"]["kh"][i] =\
                                prime[widths[0]+widths[1]+indexOfSurface(staggeredsurfaces,params,k)]*scales["kh"]
                    staggeredsurfaces[k]["data"]["kvo"][i] = prime[widths[0]+widths[1]+widths[2]]*scales["kvo"]
    return staggeredsurfaces

## return list of ids of points that are constrained by a given number of equations
def generateWhiteList(surfaces,neighbors,lowlevel,highlevel):
    idCount = {}
    for k in neighbors.keys():
        if lowlevel>k > highlevel:
            for s in neighbors[k]:
                for h in s:
                    l = surfaces[k]["ids"][h]
                    if l not in idCount.keys():
                        idCount[l] = 0
                    idCount[l] = idCount[l]+1
    whitelist = []
    for d in idCount.keys():
        if idCount[d] > 10:
            whitelist.append(d)
    return whitelist

## construct row of inverse that conserves PV, threepoint determines whether a 
## threepoint or four point grid setup will be used
def constructBetaRow(surfaces,k,distances,s,columnindexs,scales,threepoint=True,latlongrid=False):
    Arscale = scales["Ar"]
    Apsirow = [0]*(max(columnindexs)+1)

    beta = ptools.calcBeta(surfaces[k]["lats"][s[0]])
    pv =  surfaces[k]["data"]["pv"][s[0]] 
    dqnotdx = surfaces[k]["data"]["dqnotdx"][s[0]] 
    dqnotdy = surfaces[k]["data"]["dqnotdy"][s[0]] 
    u = surfaces[k]["data"]["uref"][s[0]] 
    v = surfaces[k]["data"]["vref"][s[0]] 
    f = gsw.f(surfaces[k]["lats"][s[0]])

    x = surfaces[k]["x"][s[0]] 
    y = surfaces[k]["y"][s[0]]
    r = np.sqrt(x**2 + y**2)
    #use mean gradient instead
    ## (-1/f)dAr/dy*dQnotdx
    if threepoint:
        values = [0]*3
        d01 = distances[(s[0],s[1])]
        d02 = distances[(s[0],s[2])]

        values[0] = (-1/(1*f*d02))*(dqnotdx-x*beta*pv/(f*r))\
                                   + (1/(1*f*d01))*(dqnotdy-y*beta*pv/(f*r))

        values[1] = (-1/(1*f*d01))*(dqnotdy-y*beta*pv/(f*r))
                                  

        values[2] = (1/(1*f*d02))*(dqnotdx-x*beta*pv/(f*r))\

        ## (-1/f)dAr/dx*dQnotdx
        Apsirow[columnindexs[0]] = -values[0]
        Apsirow[columnindexs[1]] = -values[1]
        Apsirow[columnindexs[2]] = -values[2]
        crow = (-u)*(dqnotdx-x*beta*pv/(f*r))+(-v)*(dqnotdy-y*beta*pv/(f*r))
 
    elif not threepoint:
        values = [0]*4
        d01 = distances[(s[0],s[1])]
        d02 = distances[(s[0],s[2])]
        d23 = distances[(s[2],s[3])]
        d13 = distances[(s[1],s[3])]

        values[0] = (-1/(2*f*d02))*(dqnotdx-x*beta*pv/(f*r))\
                                   + (1/(2*f*d01))*(dqnotdy-y*beta*pv/(f*r))

        values[1] = (-1/(2*f*d01))*(dqnotdy-y*beta*pv/(f*r))+(-1/(2*f*d13))*(dqnotdx-x*beta*pv/(f*r))\
                                  

        values[2] = (1/(2*f*d02))*(dqnotdx-x*beta*pv/(f*r))\
                                   + (1/(2*f*d23))*(dqnotdy-y*beta*pv/(f*r))
        ## (-1/f)dAr/dx*dQnotdx
        values[3] = (1/(2*f*d13))*(dqnotdx-x*beta*pv/(f*r))\
                                  + (-1/(2*f*d23))*(dqnotdy-y*beta*pv/(f*r))

        Apsirow[columnindexs[0]] = -values[0]
        Apsirow[columnindexs[1]] = -values[1]
        Apsirow[columnindexs[2]] = -values[2]
        Apsirow[columnindexs[3]] = -values[3]
        crow = (-u)*(dqnotdx-x*beta*pv/(f*r))+(-v)*(dqnotdy-y*beta*pv/(f*r))
    if np.count_nonzero(np.isnan(Apsirow)) or np.count_nonzero(np.isnan(values)) or np.count_nonzero(np.isnan(values)):
        print("boop")
    return np.asarray(Apsirow)*Arscale,np.asarray(values)*Arscale,np.asarray(crow)

## construct row of inverse that conserves salt
def constructSalRow(surfaces,k,distances,s,columnindexs,scales,threepoint=True):
    Arscale = scales["Ar"]
    Apsirow = [0]*(max(columnindexs)+1)
    dsdx = surfaces[k]["data"]["dsdx"][s[0]]
    dsdy = surfaces[k]["data"]["dsdy"][s[0]]
    f = gsw.f(surfaces[k]["lats"][s[0]])
    u = surfaces[k]["data"]["uref"][s[0]] 
    v = surfaces[k]["data"]["vref"][s[0]] 
    if threepoint: 
        values = [0]*3
        d01 = distances[(s[0],s[1])]
        d02 = distances[(s[0],s[2])]
        
        values[0] = (-1/(1*f*d02))*dsdx\
                                   + (1/(1*f*d01))*dsdy

        ## (-1/f)dAr/dx*dsdx
        values[1] =  (-1/(1*f*d01))*dsdy

        values[2] = (1/(1*f*d02))*dsdx\

        Apsirow[columnindexs[0]] = -values[0]
        Apsirow[columnindexs[1]] = -values[1]
        Apsirow[columnindexs[2]] = -values[2]
    else: 
        values = [0]*4
        d01 = distances[(s[0],s[1])]
        d02 = distances[(s[0],s[2])]
        d23 = distances[(s[2],s[3])]
        d13 = distances[(s[1],s[3])]
        
        values[0] = (-1/(2*f*d02))*dsdx\
                                   + (1/(2*f*d01))*dsdy

        ## (-1/f)dAr/dx*dsdx
        values[1] =(-1/(2*f*d13))*dsdx\
                                  + (-1/(2*f*d01))*dsdy

        values[2] = (1/(2*f*d02))*dsdx\
                                   + (1/(2*f*d23))*dsdy

        ## (1/f)dAr/dy*dsdy
        values[3] = (1/(2*f*d13))*dsdx\
                                  + (-1/(2*f*d23))*dsdy

        Apsirow[columnindexs[0]] = -values[0]
        Apsirow[columnindexs[1]] = -values[1]
        Apsirow[columnindexs[2]] = -values[2]
        Apsirow[columnindexs[3]] = -values[3]
    crow = -((u)*dsdx+(v)*dsdy)
    return np.asarray(Apsirow)*Arscale,np.asarray(values)*Arscale,np.asarray(crow)

def fillDefault(params):
    params.setdefault("debug",False)
    params.setdefault("lowerbound",2600)
    params.setdefault("upperbound",1000)
    params.setdefault("reflevel",1000)
    params.setdefault("mixs",{"kvh":True,"kvb":True,"kvo":True})
    params.setdefault("scalecoeffs",{"Ar":0.05,"kvo":5*10**-6,"kvb":5*10**-5,"kh":500})
    params.setdefault("3point",True)
    params.setdefault("edgeguard",True)
    params.setdefault("modelmixing",False)
    return params

def normOfRows(*argv):
    return np.linalg.norm(np.concatenate(argv))

#generate a list of profiles which are never the center
#at a neutral surface
def genOutliers(surfaces,neighbors,params):
    outliers={}
    for k in neighbors.keys():
        centerpoints = set()
        points = set()
        for s in neighbors[k]:
            kpv,ks = ptools.kterms(surfaces,k,s[0],params)
            if np.asarray(kpv.values()).any() and np.asarray(ks.values()).any():
                centerpoints.add(surfaces[k]["ids"][s[0]])
            for l in s[:3]:
                points.add(surfaces[k]["ids"][l])
        outliers[k] = points.difference(centerpoints)
    return outliers
    
def recenterRow(row,columnindexs,newcenter):
    centerval = row[columnindexs[0]]
    xval = row[columnindexs[1]]
    yval = row[columnindexs[2]]
    newrow = copy.copy(row)
    newrow[columnindexs[newcenter]] = centerval
    if newcenter == 1:
        newrow[columnindexs[0]] = xval
    if newcenter == 2:
        newrow[columnindexs[0]] = yval
    return row

def mixSwitch(mixing,switchs):
    scales = []
    for k in mixing.keys():
        if switchs[k]:
            scales.append(mixing[k])
    return scales
#full coupled inverse with mixing terms
def coupledInvert(surfaces,neighbors,distances,params={}):
    surfaces = copy.deepcopy(surfaces)
    Apsi,Akvb,Akh,Akvo,c=[],[],[],[],[]
    us = [0]*100000
    params = fillDefault(params)
    ##dictionary of ids to column numbers
    columndictionary = {}
    surfaces = applyRefLevel(surfaces,params["reflevel"])
    outliers = genOutliers(surfaces,neighbors,params)
    #surfaceDiagnostic(surfaces)
    totalcount = set()
    fqcrows= []
    fscrows= []
    fqs = []
    fss = []
    kvbpv=[]
    kvbs = []
    for k in Bar("adding to matrix: ").iter(sorted(neighbors.keys())[::-1]):
        for s in neighbors[k]:
            s=np.asarray(s)
            kpv,ks = ptools.kterms(surfaces,k,s[0],params)
            if len(kpv.values()) != 0 and len(ks)!=0 and \
                    (not (np.inf in np.abs(list(ks.values()))) or (np.inf in np.abs(list(kpv.values())))) \
                    and params["lowerbound"]>=k>=params["upperbound"]\
                    and surfaces[k]["ids"][s[0]]\
                    and surfaces[k]["ids"][s[0]] != -999 :
                kvbpv.append(kpv["kvb"])
                kvbs.append(ks["kvb"])
                ##find/store column indexs for each neighbor

                if params["3point"]:
                    columndictionary,columnindexs = neighborsColumnNumbers(surfaces,k,s[:3],columndictionary)
                else:
                    columndictionary,columnindexs = neighborsColumnNumbers(surfaces,k,s,columndictionary)

                ##this is a shorthad way of converting the NS to an index
                #######PVROW
                betarow,betavals, crow = constructBetaRow(surfaces,k,distances[k],s,columnindexs,params["scalecoeffs"],threepoint=params["3point"])
                n = normOfRows(np.asarray(betavals),mixSwitch(kpv,params["mixs"]))
                if params["modelmixing"]:
                    oldcrow= crow
                    fqcrows.append(crow)
                    crow = crow+ptools.calcFQ(surfaces,k,s[0],params["scalecoeffs"],kpv,distances)
                    fqs.append(ptools.calcFQ(surfaces,k,s[0],params["scalecoeffs"],kpv,distances))
                #print(betavals)
                #print(mixSwitch(kpv,params["mixs"]))
                #l = np.concatenate((betavals/n,mixSwitch(kpv,params["mixs"])/n))
                #plt.bar(range(len(l)),np.abs(l))
                ##plt.yscale("log")
                #plt.title("beta: "+str(np.sum(np.power(l,2))))
                #plt.show()
                ptools.kChecker(surfaces,k,s[0],params["scalecoeffs"])
                #Apsi.append(np.asarray(betarow)/n)

                #make rows that can fit it 
                Akvbrow, Akhrow, Akvorow = generateMixingRows(columnindexs,kpv,k,n)
                #print("pv: ",np.linalg.norm(betavals/n),np.linalg.norm(Akvbrow),np.linalg.norm(Akhrow),np.linalg.norm(Akvorow))
                #Akvb.append(Akvbrow)
                #Akh.append(Akhrow)
                #Akvo.append(Akvorow)
                plt.show()
                us[columnindexs[0]] = surfaces[k]["data"]["psiref"][s[0]]
                #c.append(crow/n)
                for i in [1,2]:
                    if surfaces[k]["ids"][s[i]] in outliers[k] and params["edgeguard"]:
                        newbetarow = recenterRow(betarow,columnindexs,i)
                        edgekpv,edgeks = ptools.kterms(surfaces,k,s[i],params,fallback=s[0])
                        if len(edgekpv.values())>0 and len(edgeks.values())>0 and not (np.inf in edgeks.values()) or (np.inf in edgekpv.values()):
                            n = normOfRows(np.asarray(betavals),mixSwitch(kpv,params["mixs"]))
                            newAkvbrow, newAkhrow, newAkvorow = generateMixingRows(columnindexs,edgekpv,k,n,i=i)
                            Apsi.append(np.asarray(newbetarow)/n)
                            Akvb.append(newAkvbrow)
                            Akh.append(newAkhrow)
                            Akvo.append(newAkvorow)
                            us[columnindexs[i]] = surfaces[k]["data"]["psiref"][s[0]]
                            c.append(crow/n)

                        

                #######SALROW
                ##make rows that can fit it 
                salrow, salvals, crow = constructSalRow(surfaces,k,distances[k],s,columnindexs,params["scalecoeffs"],threepoint=params["3point"])
                n = normOfRows(np.asarray(salvals),mixSwitch(ks,params["mixs"]))
                if params["modelmixing"]:
                    oldcrow = crow
                    fscrows.append(crow)
                    crow = crow+ptools.calcFS(surfaces,k,s[0],params["scalecoeffs"],ks,distances)
                    fss.append(ptools.calcFS(surfaces,k,s[0],params["scalecoeffs"],ks,distances))
                #l = np.concatenate((salvals/n,ks[params["mixs"]]/n))
                #plt.bar(range(len(l)),np.abs(l))
                ##plt.yscale("log")
                #plt.title("sal: "+str(np.sum(np.power(l,2))))
                #plt.show()
                Apsi.append(np.asarray(salrow)/n)

                ##im a rascal and this is a shorthad way of converting the NS to an index :P
                Akvbrow, Akhrow, Akvorow = generateMixingRows(columnindexs,ks,k,n)
                #print("sal: ",np.linalg.norm(salvals/n),np.linalg.norm(Akvbrow),np.linalg.norm(Akhrow),np.linalg.norm(Akvorow))
                Akvb.append(Akvbrow)
                Akh.append(Akhrow)
                Akvo.append(Akvorow)

                for i in [1,2]:
                    if surfaces[k]["ids"][s[i]] in outliers[k] and params["edgeguard"]:
                        newbetarow = recenterRow(salrow,columnindexs,i)
                        edgekpv,edgeks = ptools.kterms(surfaces,k,s[i],params,fallback=s[0])
                        if len(edgekpv.values())>0 and len(edgeks.values())>0 and not (np.inf in edgeks.values()) or (np.inf in edgekpv.values()):
                            n = normOfRows(np.asarray(salvals),mixSwitch(ks,params["mixs"]))
                            newAkvbrow, newAkhrow, newAkvorow = generateMixingRows(columnindexs,ks,k,n,i=i)
                            Apsi.append(np.asarray(newbetarow)/n)
                            Akvb.append(newAkvbrow)
                            Akh.append(newAkhrow)
                            Akvo.append(newAkvorow)
                            us[columnindexs[i]] = surfaces[k]["data"]["psiref"][s[0]]
                            c.append(crow/n)

                ######SAL Error row
                c.append(crow/n)
                if np.isnan(crow):
                    pdb.set_trace()
    #plt.plot(kvbpv, label="pv kvb")
    #print("kvbpv: ",np.nanmean(kvbpv))
    #plt.plot(kvbs, label="s kvb")
    #print("kvbs: ",np.nanmean(kvbs))
    #plt.legend()
    #plt.show()
    if len(Apsi)<1:
        return None, None, [None,None,None],None, None,None
    print("outliercount",len(totalcount))
    Apsi.insert(0,[1])
    Akvb.insert(0,[0])
    Akh.insert(0,[0])
    Akvo.insert(0,[0])
    c.insert(0,0)
    us.insert(0,0)
    #if params["debug"]:
        #plt.plot(np.abs(fqcrows),color="red",label="c")
        #plt.plot(np.abs(fqs),color="blue",label="FQ")
        #plt.legend()
        #plt.yscale("log")
        #plt.title("FQ vs c")
        #plt.show()
        #plt.title("FS vs c")
        #plt.plot(np.abs(fscrows),color="red",label="c")
        #plt.plot(np.abs(fss),color="blue",label="FS")
        #plt.legend()
        #plt.yscale("log")
        #plt.show()
    ##We now have all the terms we need, we just need to reshape and concatenate
    m = columndictionary["max"]+1
    us = us[:m]

    #for k in columndictionary.keys():
        #if columndictionary[k] == 6:
            #print(k)
    if params["mixs"]["kh"] and params["mixs"]["kvb"] and params["mixs"]["kvo"]:
        A,widths = combineAs([m,m,len(surfaces.keys()),1],[1,2],Apsi,Akvb,Akh,Akvo)
        print("ALL ON")
    elif params["mixs"]["kh"] and not params["mixs"]["kvb"] and params["mixs"]["kvo"]:
        print("ALL KVO, Kh")
        A,widths = combineAs([m,len(surfaces.keys()),1],[1],Apsi,Akh,Akvo)
    elif not params["mixs"]["kh"] and not params["mixs"]["kvb"] and params["mixs"]["kvo"]:
        print("ALL KVO")
        A,widths = combineAs([m,1],[],Apsi,Akvo)
    elif not params["mixs"]["kh"] and not params["mixs"]["kvb"] and not params["mixs"]["kvo"]:
        print("ALL OFF")
        A,widths = combineAs([m],[0],Apsi)
        

    if params["debug"]:
        rowCount(A)
        print("m: ",m)
        print("column count: ",np.count_nonzero(A, axis=0))
        print("#####A########")
        print(A.shape)
        print("mean Apsi: ", np.mean(np.abs(Apsi[Apsi!=0])))
        print("mean Akvb: ", np.mean(np.abs(Akvb[Akvb!=0])))
        print("mean Akh:  ", np.mean(np.abs(Akh[Akh!=0]  )))
        print("mean Akvo: ", np.mean(np.abs(Akvo[Akvo!=0])))

    c = np.matrix.transpose(np.asarray(c))
    us =  np.matrix.transpose(np.asarray(us))
    j,[VT, D, U] = SVDdecomp(A,n_elements=int(A.shape[1]))

    prime = np.matmul(j,c)

    if params["debug"]:
        graphError(A,np.concatenate((us,[0]*(A.shape[1]-len(us)))),prime)
        rdivc(D)

    errors = error(A,np.concatenate((us,[0]*(A.shape[1]-len(us)))),prime)
    #calculatediapycnalv(surfaces,prime,params,widths,columndictionary)
    metadata = copy.copy(params)
    metadata["condition"] = condition(D)
    metadata["error"] = np.sum(np.abs(errors[1]))
    surfaces = applyPrime(surfaces,prime,columndictionary,params,widths,mixing=True)

    return {"surfaces":surfaces,"coldict":columndictionary, "svddecomp":[VT,D,U],"matrixsetup":A,"errors":errors,"metadata":metadata}


#the number of values per row
def rowCount(A):
    number = []
    for i in range(A.shape[0]):
        number.append(np.count_nonzero(A[i]))
    print("in each row: ",np.unique(number))

def condition(D):
    d = D.diagonal()
    d = 1.0/d
    return d[0]/d[-1]

## essentially calculates singular values and matrix conditions
def rdivc(D):
    d = D.diagonal()
    d = 1.0/d
    plt.title("Singular Values")
    print("above 0.04: ",len(d[np.where(d>0.04)]))
    plt.scatter(range(len(d[:])),(d[:]))
    plt.show()
    plt.title("Largest Singular Values Divided by Smallest")
    plt.scatter(range(len(d[:])),(d[0]/d[:]))
    plt.show()

## export all necessary information about inverse as a mat file
def exportMat(surfaces,columndict,svds,A):
    outmat={"dqnotdx":[],"dqnotdy":[],"dsdx":[],\
        "dsdy":[],"pv":[],"VT":None,"D":None,"U":None,\
        "A":None,"A-1":None,"beta":[],"lats":[],"lons":[],"y":[],"x":[],\
        "psiref":[],"dpsidx":[],"dpsidy":[]}
    outmat["VT"] = svds[0]
    outmat["D"] = svds[1]
    outmat["U"] = svds[2]
    outmat["A"] =A
    outmat["A-1"] =SVDCalc(svds[0],svds[1],svds[2])
    for k in surfaces.keys():
        for l in outmat.keys():
            if type(outmat[l]) == list:
                outmat[l].append([])
    for eyed in columndict.keys():
        for k in sorted(surfaces.keys()):
            j = int(int(k-200)/200)
            print(j)
            if eyed in surfaces[k]["ids"]:
                found = np.where(np.asarray(surfaces[k]["ids"])==eyed)[0][0]
                outmat["dqnotdx"][j].append(surfaces[k]["data"]["dqnotdx"][found])
                outmat["dqnotdy"][j].append(surfaces[k]["data"]["dqnotdy"][found])
                outmat["dsdx"][j].append(surfaces[k]["data"]["dsdx"][found])
                outmat["dsdy"][j].append(surfaces[k]["data"]["dsdx"][found])
                outmat["pv"][j].append(surfaces[k]["data"]["pv"][found])
                outmat["beta"][j].append(surfaces[k]["data"]["beta"][found])
                outmat["lats"][j].append(surfaces[k]["lats"][found])
                outmat["lons"][j].append(surfaces[k]["lons"][found])
                outmat["x"][j].append(surfaces[k]["x"][found])
                outmat["y"][j].append(surfaces[k]["y"][found])
                outmat["psiref"][j].append(surfaces[k]["data"]["psiref"][found])
                outmat["dpsidx"][j].append(surfaces[k]["data"]["dpsidx"][found])
                outmat["dpsidy"][j].append(surfaces[k]["data"]["dpsidy"][found])
            else:
                outmat["dqnotdx"][j].append(np.nan)
                outmat["dqnotdy"][j].append(np.nan)
                outmat["dsdx"][j].append(np.nan)
                outmat["dsdy"][j].append(np.nan)
                outmat["pv"][j].append(np.nan)
                outmat["beta"][j].append(np.nan)
                outmat["lats"][j].append(np.nan)
                outmat["lons"][j].append(np.nan)
                outmat["x"][j].append(np.nan)
                outmat["y"][j].append(np.nan)
                outmat["psiref"][j].append(np.nan)
                outmat["dpsidx"][j].append(np.nan)
                outmat["dpsidy"][j].append(np.nan)
    scipy.io.savemat("outmats",outmat)




#a is an array to rectangularize
#l is the maximum length
def rectangularize(a,l):
    new = []
    lengths =[]
    for b in a:
        new.append(np.concatenate( (b,([0]*(l-len(b)))) ))
        lengths.append(len(b))
    return new

##add up a bunch of matrixes
def concat(argv):
    out = []
    for j in range(len(argv[0])):
        row = np.array([])
        for b in argv:
            row = np.concatenate((row,b[j]))
        out.append(row)
    return np.asarray(out)

def rectAndWidths(maxlengths,totrim,arrays):
    new = []
    widths = []
    for i in range(len(arrays)):
        x  = np.asarray(rectangularize(arrays[i],maxlengths[i]))
        print("====")
        print("before: ",x.shape)
        if np.all(x == 0, axis=0).any():
            print("removing 0 columns")
        if i in totrim:
            print(~np.all(x == 0, axis=0))
            x = x[:,~np.all(x == 0, axis=0)]
        print("after: ",x.shape)
        print("====")
        widths.append(x.shape[1])
        new.append(x)
    return new,widths

## square off all of the matrixes and combine them
def combineAs(maxlengths,totrim,*argv):
    new,widths = rectAndWidths(maxlengths,totrim,argv)
    return concat(new),widths

#graph solution of mixing terms
def graphKs(prime,m):
    kvb = prime[m:2*m]
    kh = prime[2*m:2*m+8]
    kvo =prime[-1]
    fig, (ax1,ax2) = plt.subplots(1,2)
    ax1.set_title("KVB")
    ax1.plot(range(len(kvb)),kvb+kvo)
    ax2.set_title("KH")
    ax2.plot(range(len(kh)),kh)
    plt.show()
    print(kvo)

#routing command to perform inverse 
def invert(kind,surfaces,neighbors=None,distances=None,params={},reflevel=1000,debug=False,lowlevel=2000,highlevel=1000):
    if kind == "simple":
        return simpleInvert(surfaces,reflevel,debug)
    if kind == "simplesalt":
        return simplesaltInvert(surfaces,reflevel,debug)
    if kind == "complexsalt":
        return complexSaltInvert(surfaces,reflevel,debug)
    if kind == "complex":
        return complexInvert(surfaces,reflevel,debug)
    if kind == "coupled":
        return coupledInvert(surfaces,neighbors,distances,params)
    else:
        print("Sorry I have no clue what inversion you are talking about")


