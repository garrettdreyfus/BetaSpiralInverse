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
    Akvbrow[columnindexs[i]]=mixing[1]/n
    Akhrow = [0]*(int(k/200))
    Akhrow[int(k/200-1)] = mixing[2]/n
    Akvorow = [mixing[0]/n]
    return Akvbrow,Akhrow,Akvorow


#simple pointwise inverse only conserves pv
def simpleInvert(surfaces,reflevel=1000,debug=False):
    outsurfaces = copy.deepcopy(surfaces)
    for k in surfaces.keys():
        outsurfaces[k]["data"]["u"].fill(np.nan)  
        outsurfaces[k]["data"]["v"].fill(np.nan)  
        outsurfaces[k]["data"]["uabs"].fill(np.nan)  
        outsurfaces[k]["data"]["vabs"].fill(np.nan)  
    for index in  Bar('Simple Invert: ').iter(range(len(surfaces[reflevel]["x"]))):
        eyed = int(surfaces[reflevel]["ids"][index])
        #print("id: ",eyed)
        ns = []
        us = [[],[]]
        b = []
        c = []
        prime = [[],[]]
        for k in surfaces.keys():
            found = np.where(np.asarray(surfaces[k]["ids"])==eyed)
            ns.append((k,found))
            if len(found)!=0 and len(found[0]) != 0:
                found = found[0][0]
                x = surfaces[k]["x"][found]
                y = surfaces[k]["y"][found]
                r = np.sqrt(x**2 + y**2)
                hx = surfaces[k]["data"]["hx"][found]
                hy = surfaces[k]["data"]["hy"][found]
                pres = surfaces[k]["data"]["pres"][found]
                f = gsw.f(surfaces[k]["lats"][found])
                beta = ptools.calcBeta(surfaces[k]["lats"][found])
                pv = surfaces[k]["data"]["pv"][found]
                dqnotdx = surfaces[k]["data"]["dqnotdx"][found]
                dqnotdy = surfaces[k]["data"]["dqnotdy"][found]

                if debug and (np.isnan(hx) or np.isnan(hy) or np.isnan(x) or np.isnan(y)):
                    print("pres is nan: ",np.isnan(pres))
                    print("hx is nan: ",np.isnan(hx))
                    print("hy is nan: ",np.isnan(hy))
                    print("x is nan: ",np.isnan(x))
                    print("y is nan: ",np.isnan(y))
                    print("something here is nan")
                if k>=1000 and not(np.isnan(hx) or np.isnan(hy) or np.isnan(x) or np.isnan(y)):
                    bvec = (hx+(beta*x)/(f*r),hy+(beta*y)/(f*r))
                    #bvec=(dqnotdx-x*beta*pv/(f*r),dqnotdy-y*beta*pv/(f*r))
                    b.append(norm(bvec))
                    u = (surfaces[k]["data"]["u"][found] - surfaces[reflevel]["data"]["u"][index])
                    v = (surfaces[k]["data"]["v"][found] - surfaces[reflevel]["data"]["v"][index])
                    us[0].append(u)
                    us[1].append(v)
                    c.append(np.dot(norm(bvec),(-u,-v)))

        if len(b)>0:
            b = np.asarray(b)
            c = np.matrix.transpose(np.asarray(c))
            us = np.asarray(us)
            #j = np.linalg.inv(np.matmul(np.matrix.transpose(b),b))
            #j = SVDdecomp(b,n_elements=4)
            j,[VT, D, U]  = SVDdecomp(b,n_elements=2,full=True)
            s = np.diag(D)
            print(b)
            print(1.0/s)
            print(s[-1]/s[0])

            prime = np.matmul(j,c)
            b = b.T
            errorbefore = []
            for i in range(len(b[0])):
                errorbefore.append(b[0][i]*(us[0][i]+prime[0]) + b[1][i]*(us[1][i]+prime[1]))
            errorafter = []
            for i in range(len(b[0])):
                errorafter.append(b[0][i]*(us[0][i]) + b[1][i]*(us[1][i]))

            for i in range(len(ns)):
                uref = surfaces[reflevel]["data"]["u"][index]
                vref = surfaces[reflevel]["data"]["v"][index]
                outsurfaces[ns[i][0]]["data"]["uprime"][ns[i][1]] = prime[0] 
                outsurfaces[ns[i][0]]["data"]["uabs"][ns[i][1]] = prime[0] + surfaces[ns[i][0]]["data"]["u"][ns[i][1]]-uref
                outsurfaces[ns[i][0]]["data"]["u"][ns[i][1]] = surfaces[ns[i][0]]["data"]["u"][ns[i][1]]-uref
                outsurfaces[ns[i][0]]["data"]["vprime"][ns[i][1]] = prime[1]
                outsurfaces[ns[i][0]]["data"]["vabs"][ns[i][1]] = prime[1] + surfaces[ns[i][0]]["data"]["v"][ns[i][1]]-vref
                outsurfaces[ns[i][0]]["data"]["v"][ns[i][1]] = surfaces[ns[i][0]]["data"]["v"][ns[i][1]]-vref

    return outsurfaces

#simpe pointwise inverse conserves pv and salt
def simplesaltInvert(surfaces,reflevel=1000,debug=False):
    outsurfaces = copy.deepcopy(surfaces)
    for k in surfaces.keys():
        outsurfaces[k]["data"]["u"].fill(np.nan)  
        outsurfaces[k]["data"]["v"].fill(np.nan)  
        outsurfaces[k]["data"]["uabs"].fill(np.nan)  
        outsurfaces[k]["data"]["vabs"].fill(np.nan)  
    for index in  Bar('Simple Salt Inv refevel'+str(reflevel)+': ').iter(range(len(surfaces[reflevel]["x"]))):
        eyed = int(surfaces[reflevel]["ids"][index])
        #print("id: ",eyed)
        ns = []
        us = []
        b = []
        c = []
        prime = [[],[]]
        for k in surfaces.keys():
            found = np.where(np.asarray(surfaces[k]["ids"])==eyed)
            if len(found)!=0 and len(found[0]) != 0:
                ns.append((k,found))
                found = found[0][0]
                x = surfaces[k]["x"][found]
                y = surfaces[k]["y"][found]
                r = np.sqrt(x**2 + y**2)
                hx = surfaces[k]["data"]["hx"][found]
                hy = surfaces[k]["data"]["hy"][found]
                dsdx = surfaces[k]["data"]["dsdx"][found]
                dsdy = surfaces[k]["data"]["dsdy"][found]
                dqnotdx = surfaces[k]["data"]["dqnotdx"][found]
                dqnotdy = surfaces[k]["data"]["dqnotdy"][found]
                pv = surfaces[k]["data"]["pv"][found]
                pres = surfaces[k]["data"]["pres"][found]
                f = gsw.f(surfaces[k]["lats"][found])
                beta = ptools.calcBeta(surfaces[k]["lats"][found])

                if debug and (np.isnan(hx) or np.isnan(hy) or np.isnan(x) or np.isnan(y)):
                    print("pres is nan: ",np.isnan(pres))
                    print("hx is nan: ",np.isnan(hx))
                    print("hy is nan: ",np.isnan(hy))
                    print("x is nan: ",np.isnan(x))
                    print("y is nan: ",np.isnan(y))
                    print("something here is nan")
                if k>= 1000 and not(np.isnan(hx) or np.isnan(hy) or np.isnan(x) or np.isnan(y)):
                    u = (surfaces[k]["data"]["u"][found] - surfaces[reflevel]["data"]["u"][index])
                    v = (surfaces[k]["data"]["v"][found] - surfaces[reflevel]["data"]["v"][index])
                    #surfaces[k]["data"]["uabs"][found]=0
                    pvvec=((hx+(beta*x)/(f*r)),(hy+(beta*y)/(f*r)))
                    #pvvec=(dqnotdx-x*beta*pv/(f*r),dqnotdy-y*beta*pv/(f*r))
                    #pvvec=pvvec/np.linalg.norm(pvvec)
                    b.append(norm(pvvec))
                    us.append((u,v))
                    c.append(np.dot((-u,-v),norm(pvvec)))

                    svec=(dsdx,dsdy)
                    #svec=svec/np.linalg.norm(svec)
                    b.append(norm(svec))
                    c.append(np.dot((-u,-v),norm(svec)))
                    us.append((u,v))
        if len(b)>0:
            b = np.asarray(b)
            c = np.matrix.transpose(np.asarray(c))
            us = np.asarray(us)
            #j = np.linalg.inv(np.matmul(np.matrix.transpose(b),b))
            #j = SVDdecomp(b,n_elements=4)
            j,[VT,D,U] = SVDdecomp(b,n_elements=3,full=True)
            prime = np.matmul(j,c)
            
            for i in range(len(ns)):
                uref = surfaces[reflevel]["data"]["u"][index]
                vref = surfaces[reflevel]["data"]["v"][index]
                outsurfaces[ns[i][0]]["data"]["uprime"][ns[i][1]] = prime[0] 
                outsurfaces[ns[i][0]]["data"]["uabs"][ns[i][1]] = prime[0] + surfaces[ns[i][0]]["data"]["u"][ns[i][1]]-uref
                outsurfaces[ns[i][0]]["data"]["u"][ns[i][1]] = surfaces[ns[i][0]]["data"]["u"][ns[i][1]]-uref
                outsurfaces[ns[i][0]]["data"]["vprime"][ns[i][1]] = prime[1]
                outsurfaces[ns[i][0]]["data"]["vabs"][ns[i][1]] = prime[1] + surfaces[ns[i][0]]["data"]["v"][ns[i][1]]-vref
                outsurfaces[ns[i][0]]["data"]["v"][ns[i][1]] = surfaces[ns[i][0]]["data"]["v"][ns[i][1]]-vref
    return outsurfaces

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

##performs a pointwise ivnerse that conserves pv and salt
## and also includes Fs mixing terms
def complexSaltInvert(surfaces,reflevel=1000,debug=False):
    outsurfaces = copy.deepcopy(surfaces)
    for k in surfaces.keys():
        outsurfaces[k]["data"]["u"].fill(np.nan)  
        outsurfaces[k]["data"]["v"].fill(np.nan)  
        outsurfaces[k]["data"]["uabs"].fill(np.nan)  
        outsurfaces[k]["data"]["vabs"].fill(np.nan)  
    for index in  Bar('Complex Salt Invert: ').iter(range(len(surfaces[reflevel]["x"]))):
        eyed = int(surfaces[reflevel]["ids"][index])
        #print("id: ",eyed)
        ns = []
        us = []
        b = []
        c = []
        prime = [[],[]]
        for k in sorted(list(surfaces.keys())):
            found = np.where(np.asarray(surfaces[k]["ids"])==eyed)
            ns.append((k,found))
            if len(found)!=0 and len(found[0]) != 0:
                found = found[0][0]
                x = surfaces[k]["x"][found]
                y = surfaces[k]["y"][found]
                r = np.sqrt(x**2 + y**2)
                hx = surfaces[k]["data"]["hx"][found]
                hy = surfaces[k]["data"]["hy"][found]
                dsdx = surfaces[k]["data"]["dsdx"][found]
                dsdy = surfaces[k]["data"]["dsdy"][found]
                pres = surfaces[k]["data"]["pres"][found]
                alpha = surfaces[k]["data"]["alpha"][found] 
                betaTherm = surfaces[k]["data"]["beta"][found] 
                dsdz =  surfaces[k]["data"]["dsdz"][found] 
                d2sdx2 =  surfaces[k]["data"]["d2sdx2"][found] 
                d2sdy2 =  surfaces[k]["data"]["d2sdy2"][found] 
                dalphadtheta = surfaces[k]["data"]["dalphadtheta"][found] 
                dalphads = surfaces[k]["data"]["dalphads"][found] 
                dalphadp = surfaces[k]["data"]["dalphadp"][found] 
                dbetadp = surfaces[k]["data"]["dbetadp"][found] 
                dbetads = surfaces[k]["data"]["dbetads"][found] 
                dtdx = surfaces[k]["data"]["dtdx"][found] 
                dtdy = surfaces[k]["data"]["dtdy"][found] 
                dqnotdx = surfaces[k]["data"]["dqnotdx"][found] 
                dqnotdy = surfaces[k]["data"]["dqnotdy"][found] 
                dpdx = surfaces[k]["data"]["dpdx"][found] 
                dpdy = surfaces[k]["data"]["dpdy"][found] 
                alphat = dalphadtheta+2*(alpha/betaTherm)*dalphads-(alpha**2/betaTherm**2)*dbetads
                alphap = dalphadp -(alpha/betaTherm)*dbetadp
                pv =  surfaces[k]["data"]["pv"][found] 
                doublets =  surfaces[k]["data"]["d2thetads2"][found] 
                CKVB =  surfaces[k]["data"]["CKVB"][found] 
                f = gsw.f(surfaces[k]["lats"][found])
                beta = ptools.calcBeta(surfaces[k]["lats"][found])
                isitnan = [alpha,betaTherm,dsdz,hx,hy,dsdx,dsdy,pres,d2sdx2,d2sdy2,\
                          dalphadtheta,dalphads,dalphadp,dbetadp,dbetads,dtdx,dtdy,\
                          dqnotdx,dqnotdy,dpdx,dpdy,alphat,alphap,pv,doublets,CKVB,\
                          beta]

                if debug and (np.isnan(isitnan).any()):
                    print("pres is nan: ",np.isnan(pres))
                    print("hx is nan: ",np.isnan(hx))
                    print("hy is nan: ",np.isnan(hy))
                    print("x is nan: ",np.isnan(x))
                    print("y is nan: ",np.isnan(y))
                    print("something here is nan")
                if debug and (np.isnan(isitnan[3:]).any()) and k!=200:
                    print(isitnan)
                if k>=1000 and not (np.isnan(isitnan).any()):
                    u = (surfaces[k]["data"]["u"][found] - surfaces[reflevel]["data"]["u"][index])
                    v = (surfaces[k]["data"]["v"][found] - surfaces[reflevel]["data"]["v"][index])
                    us.append((u,v,0,0,0))
                    us.append((u,v,0,0,0))
                    ######################################
                    ###Potential vorticity row
                    #left side
                    pvvec=((hx+(beta*x)/(f*r)),(hy+(beta*y)/(f*r)),0,0,0)
                    pvnorm=np.linalg.norm(pvvec)
                    b.append(pvvec/pvnorm)
                    #right side
                    c.append(np.dot((-u,-v),pvvec[0:2]/pvnorm))
                    ######################################
                    ###Salinity row
                    #left side
                    kvo = -alpha*f*(1/pv)*(dsdz**3)*doublets
                    kvb = kvo*CKVB

                    khpart1 = (f/pv)*dsdz*(alphat*(dtdx**2 + dtdy**2)+alphap*(dtdx*dpdx+dtdy*dpdy))
                    khpart2 = (d2sdx2+d2sdy2)-2*(dqnotdx*dsdx + dqnotdy*dsdy)/pv
                    kh = khpart1 + khpart2

                    svec=(dsdx,dsdy,-kvb,-kvo,-kh)
                    snorm = np.linalg.norm(svec)
                    b.append(svec/snorm)

                    #right side
                    c.append(np.dot((-u,-v),svec[0:2]/snorm))

        if len(b)>4:
            b = np.asarray(b)
            c = np.matrix.transpose(np.asarray(c))
            us = np.asarray(us)
            j = SVDdecomp(b,n_elements=3)
            prime = np.matmul(j,c)
            error = []
            #graphError(b,us,prime)
            #print(prime)
            for i in range(len(ns)):
                uref = surfaces[reflevel]["data"]["u"][index]
                vref = surfaces[reflevel]["data"]["v"][index]
                outsurfaces[ns[i][0]]["data"]["uprime"][ns[i][1]] = uref
                outsurfaces[ns[i][0]]["data"]["uabs"][ns[i][1]] = prime[0] + surfaces[ns[i][0]]["data"]["u"][ns[i][1]]-uref
                outsurfaces[ns[i][0]]["data"]["u"][ns[i][1]] = surfaces[ns[i][0]]["data"]["u"][ns[i][1]]-uref
                outsurfaces[ns[i][0]]["data"]["vprime"][ns[i][1]] = vref
                outsurfaces[ns[i][0]]["data"]["vabs"][ns[i][1]] = prime[1] + surfaces[ns[i][0]]["data"]["v"][ns[i][1]]-vref
                outsurfaces[ns[i][0]]["data"]["v"][ns[i][1]] = surfaces[ns[i][0]]["data"]["v"][ns[i][1]]-vref
    return outsurfaces

##performs a pointwise ivnerse that conserves pv and salt
## and also includes Fs and Fq mixing terms
def complexInvert(surfaces,reflevel=1000,debug=False):
    outsurfaces = copy.deepcopy(surfaces)
    for k in surfaces.keys():
        outsurfaces[k]["data"]["u"].fill(np.nan)  
        outsurfaces[k]["data"]["v"].fill(np.nan)  
        outsurfaces[k]["data"]["uabs"].fill(np.nan)  
        outsurfaces[k]["data"]["vabs"].fill(np.nan)  
    stats=[[],[],[]]
    for index in  Bar('Complex Salt Invert: ').iter(range(len(surfaces[reflevel]["x"]))):
        eyed = int(surfaces[reflevel]["ids"][index])
        #print("id: ",eyed)
        ns = []
        us = []
        b = []
        c = []
        prime = [[],[]]
        for k in sorted(list(surfaces.keys())):
            found = np.where(np.asarray(surfaces[k]["ids"])==eyed)
            ns.append((k,found))
            if len(found)!=0 and len(found[0]) != 0:
                found = found[0][0]
                kpv,ks = ptools.kterms(surfaces,k,found)
                if k>=1000 and kpv and ks:
                    x = surfaces[k]["x"][found]
                    y = surfaces[k]["y"][found]
                    r = np.sqrt(x**2 + y**2)
                    hx = surfaces[k]["data"]["hx"][found]
                    dsdx = surfaces[k]["data"]["dsdx"][found]
                    dsdy = surfaces[k]["data"]["dsdy"][found]
                    hy = surfaces[k]["data"]["hy"][found]
                    f = gsw.f(surfaces[k]["lats"][found])
                    beta = ptools.calcBeta(surfaces[k]["lats"][found])
                    u = (surfaces[k]["data"]["u"][found] - surfaces[reflevel]["data"]["u"][index])
                    v = (surfaces[k]["data"]["v"][found] - surfaces[reflevel]["data"]["v"][index])
                    us.append((u,v,0,0,0))
                    us.append((u,v,0,0,0))
                    ######################################
                    ###Potential vorticity row
                    #left side
                    pvvec=((hx+(beta*x)/(f*r)),(hy+(beta*y)/(f*r)),kpv[1],kpv[0],kpv[2])
                    pvnorm=np.linalg.norm(pvvec)
                    b.append(pvvec/pvnorm)
                    #right side
                    c.append(np.dot((-u,-v),pvvec[0:2]/pvnorm))
                    ######################################
                    ###Salinity row
                    #left side
                    svec=(dsdx,dsdy,ks[1],ks[0],ks[2])
                    snorm = np.linalg.norm(svec)
                    b.append(svec/snorm)

                    #right side
                    c.append(np.dot((-u,-v),svec[0:2]/snorm))

        if len(b)>4:
            b = np.asarray(b)
            c = np.matrix.transpose(np.asarray(c))
            us = np.asarray(us)
            j = SVDdecomp(b,n_elements=4)
            prime = np.matmul(j,c)
            stats[0].append(prime[2])
            stats[1].append(prime[3])
            stats[2].append(prime[4])
            error = []
            #graphError(b,us,prime)
            #print(prime)
            for i in range(len(ns)):
                uref = surfaces[reflevel]["data"]["u"][index]
                vref = surfaces[reflevel]["data"]["v"][index]
                outsurfaces[ns[i][0]]["data"]["uprime"][ns[i][1]] = uref
                outsurfaces[ns[i][0]]["data"]["uabs"][ns[i][1]] = prime[0] + surfaces[ns[i][0]]["data"]["u"][ns[i][1]]-uref
                outsurfaces[ns[i][0]]["data"]["u"][ns[i][1]] = surfaces[ns[i][0]]["data"]["u"][ns[i][1]]-uref
                outsurfaces[ns[i][0]]["data"]["vprime"][ns[i][1]] = vref
                outsurfaces[ns[i][0]]["data"]["vabs"][ns[i][1]] = prime[1] + surfaces[ns[i][0]]["data"]["v"][ns[i][1]]-vref
                outsurfaces[ns[i][0]]["data"]["v"][ns[i][1]] = surfaces[ns[i][0]]["data"]["v"][ns[i][1]]-vref
    print(np.mean(stats[0]),np.mean(stats[1]),np.mean(stats[2]))
    return outsurfaces

#get the column number of a parameter
#generates new number if necessary
def getColumnNumber(eyedict,eyed):
    #assign each eyed a column number 
    if "max" not in eyedict.keys():
        eyedict["max"]=-1
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
        for l in range(len(surfaces[k]["data"]["psi"])):
            found = np.where(np.asarray(surfaces[reflevel]["ids"])==surfaces[k]["ids"][l])
            if len(found) != 0 and len(found[0]) != 0:
                surfaces[k]["data"]["psiref"][l] = surfaces[k]["data"]["psi"][l] - surfaces[reflevel]["data"]["psi"][found[0][0]]
                surfaces[k]["data"]["uref"][l] = (surfaces[k]["data"]["u"][l]-surfaces[reflevel]["data"]["u"][found[0][0]])
                surfaces[k]["data"]["vref"][l] = (surfaces[k]["data"]["v"][l]-surfaces[reflevel]["data"]["v"][found[0][0]])
    return surfaces

## apply the solution to surfaces
def applyPrime(staggeredsurfaces,prime,coldict,params,widths,mixing=False):
    scales = params["scalecoeffs"]
    for k in Bar("adding solutions: ").iter(staggeredsurfaces.keys()):
        staggeredsurfaces[k]["data"]["psinew"] =  np.full(len(staggeredsurfaces[k]["lons"]),np.nan)
        staggeredsurfaces[k]["data"]["psisol"] =  np.full(len(staggeredsurfaces[k]["lons"]),np.nan)
        staggeredsurfaces[k]["data"]["kvb"] =  np.full(len(staggeredsurfaces[k]["lons"]),np.nan)
        staggeredsurfaces[k]["data"]["kvo"] =  np.full(len(staggeredsurfaces[k]["lons"]),np.nan)
        staggeredsurfaces[k]["data"]["usol"] =  np.full(len(staggeredsurfaces[k]["lons"]),np.nan)
        staggeredsurfaces[k]["data"]["vsol"] =  np.full(len(staggeredsurfaces[k]["lons"]),np.nan)
        for i in range(len(staggeredsurfaces[k]["data"]["ids"])):
            eyed = staggeredsurfaces[k]["data"]["ids"][i] 
            if eyed in coldict.keys():
                staggeredsurfaces[k]["data"]["psinew"][i] = staggeredsurfaces[k]["data"]["psiref"][i] + prime[coldict[eyed]]*scales[0]
                staggeredsurfaces[k]["data"]["psisol"][i] = prime[coldict[eyed]]*scales[0]
                #if params["mixs"] == [True,True,True]:
                    #staggeredsurfaces[k]["data"]["kvb"][i] = prime[widths[0]+coldict[eyed]]
                    #staggeredsurfaces[k]["data"]["kvh"][i] = prime[widths[0]+widths[1]+]
                    #staggeredsurfaces[k]["data"]["kvo"][i] = prime[widths[0]+widths[1](2*m)+int(k/200.0)]
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
def constructBetaRow(surfaces,k,distances,s,columnindexs,scales,threepoint=True):
    Arscale = scales[0]
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

        values[0] = (-1/(2*f*d02))*(dqnotdx-x*beta*pv/(f*r))\
                                   + (1/(2*f*d01))*(dqnotdy-y*beta*pv/(f*r))

        values[1] = (-1/(2*f*d01))*(dqnotdy-y*beta*pv/(f*r))
                                  

        values[2] = (1/(2*f*d02))*(dqnotdx-x*beta*pv/(f*r))\

        ## (-1/f)dAr/dx*dQnotdx
        Apsirow[columnindexs[0]] = values[0]
        Apsirow[columnindexs[1]] = values[1]
        Apsirow[columnindexs[2]] = values[2]
 
    else:
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

        Apsirow[columnindexs[0]] = values[0]
        Apsirow[columnindexs[1]] = values[1]
        Apsirow[columnindexs[2]] = values[2]
        Apsirow[columnindexs[3]] = values[3]
    crow = (-u)*(dqnotdx-x*beta*pv/(f*r))+(-v)*(dqnotdy-y*beta*pv/(f*r))
    return np.asarray(Apsirow)*Arscale,np.asarray(values)*Arscale,np.asarray(crow)

## construct row of inverse that conserves salt
def constructSalRow(surfaces,k,distances,s,columnindexs,scales,threepoint=True):
    Arscale = scales[0]
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
        
        values[0] = (-1/(2*f*d02))*dsdx\
                                   + (1/(2*f*d01))*dsdy

        ## (-1/f)dAr/dx*dsdx
        values[1] =  (-1/(2*f*d01))*dsdy

        values[2] = (1/(2*f*d02))*dsdx\

        Apsirow[columnindexs[0]] = values[0]
        Apsirow[columnindexs[1]] = values[1]
        Apsirow[columnindexs[2]] = values[2]
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

        Apsirow[columnindexs[0]] = values[0]
        Apsirow[columnindexs[1]] = values[1]
        Apsirow[columnindexs[2]] = values[2]
        Apsirow[columnindexs[3]] = values[3]
    crow = -((u)*dsdx+(v)*dsdy)
    return np.asarray(Apsirow)*Arscale,np.asarray(values)*Arscale,np.asarray(crow)

def fillDefault(params):
    params.setdefault("debug",False)
    params.setdefault("threepoint",True)
    params.setdefault("lowerbound",2600)
    params.setdefault("upperbound",1000)
    params.setdefault("reflevel",1000)
    params.setdefault("mixs",[True,False,True])
    params.setdefault("scalecoeffs",[0.05,5*10**-6,5*10**-5,500])
    params.setdefault("3point",True)
    return params

def normOfRows(*argv):
    return np.linalg.norm(np.concatenate(argv))

#generate a list of profiles which are never the center
#at a neutral surface
def genOutliers(surfaces,neighbors):
    outliers={}
    for k in neighbors.keys():
        centerpoints = set()
        points = set()
        for s in neighbors[k]:
            kpv,ks = ptools.kterms(surfaces,k,s[0],[1,1,1,1])
            if kpv.any() and ks.any():
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

#full coupled inverse with mixing terms
def coupledInvert(surfaces,neighbors,distances,params={},edgeguard=True):
    surfaces = copy.deepcopy(surfaces)
    Apsi,Akvb,Akh,Akvo,c=[],[],[],[],[]
    us = [0]*100000
    params = fillDefault(params)
    ##dictionary of ids to column numbers
    columndictionary = {}
    surfaces = applyRefLevel(surfaces,params["reflevel"])
    outliers = genOutliers(surfaces,neighbors)
    totalcount = set()
    for k in Bar("adding to matrix: ").iter(sorted(neighbors.keys())[::-1]):
        for s in neighbors[k]:
            s=np.asarray(s)
            kpv,ks = ptools.kterms(surfaces,k,s[0],params["scalecoeffs"])
            if kpv.any() and ks.any() and params["lowerbound"]>=k>=params["upperbound"] and surfaces[k]["ids"][s[0]]:
                ##find/store column indexs for each neighbor
                if params["3point"]:
                    columndictionary,columnindexs = neighborsColumnNumbers(surfaces,k,s[:3],columndictionary)
                else:
                    columndictionary,columnindexs = neighborsColumnNumbers(surfaces,k,s,columndictionary)

                ##this is a shorthad way of converting the NS to an index
                #######PVROW
                betarow,betavals, crow = constructBetaRow(surfaces,k,distances[k],s,columnindexs,params["scalecoeffs"],threepoint=params["3point"])
                n = normOfRows(np.asarray(betavals),kpv[params["mixs"]])
                #l = np.concatenate((betavals/n,kpv[params["mixs"]]/n))
                #plt.bar(range(len(l)),np.abs(l))
                ##plt.yscale("log")
                #plt.title("beta: "+str(np.sum(np.power(l,2))))
                #plt.show()
                #ptools.kChecker(surfaces,k,s[0],params["scalecoeffs"])
                Apsi.append(np.asarray(betarow)/n)

                ##make rows that can fit it 
                Akvbrow, Akhrow, Akvorow = generateMixingRows(columnindexs,kpv,k,n)
                Akvb.append(Akvbrow)
                Akh.append(Akhrow)
                Akvo.append(Akvorow)
                us[columnindexs[0]] = surfaces[k]["data"]["psiref"][s[0]]
                c.append(crow/n)

                for i in [1,2]:
                    if surfaces[k]["ids"][s[i]] in outliers[k]:
                        newbetarow = recenterRow(betarow,columnindexs,i)
                        kpv,ks = ptools.kterms(surfaces,k,s[i],params["scalecoeffs"],fallback=s[0])
                        if kpv.any() and ks.any():
                            n = normOfRows(np.asarray(betavals),kpv[params["mixs"]])
                            newAkvbrow, newAkhrow, newAkvorow = generateMixingRows(columnindexs,kpv,k,n,i=i)
                            Apsi.append(np.asarray(newbetarow)/n)
                            Akvb.append(newAkvbrow)
                            Akh.append(newAkhrow)
                            Akvo.append(newAkvorow)
                            us[columnindexs[i]] = surfaces[k]["data"]["psiref"][s[0]]
                            c.append(crow/n)

                        

                #######SALROW
                ##make rows that can fit it 
                salrow, salvals, crow = constructSalRow(surfaces,k,distances[k],s,columnindexs,params["scalecoeffs"],threepoint=params["3point"])
                n = normOfRows(np.asarray(salvals),ks[params["mixs"]])
                #l = np.concatenate((salvals/n,ks[params["mixs"]]/n))
                #plt.bar(range(len(l)),np.abs(l))
                ##plt.yscale("log")
                #plt.title("sal: "+str(np.sum(np.power(l,2))))
                #plt.show()
                Apsi.append(np.asarray(salrow)/n)

                ##im a rascal and this is a shorthad way of converting the NS to an index :P
                Akvbrow, Akhrow, Akvorow = generateMixingRows(columnindexs,ks,k,n)
                Akvb.append(Akvbrow)
                Akh.append(Akhrow)
                Akvo.append(Akvorow)

                for i in [1,2]:
                    if surfaces[k]["ids"][s[i]] in outliers[k]:
                        newbetarow = recenterRow(salrow,columnindexs,i)
                        kpv,ks = ptools.kterms(surfaces,k,s[i],params["scalecoeffs"],fallback=s[0])
                        n = normOfRows(np.asarray(salvals),ks[params["mixs"]])
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

    if len(Apsi)<1:
        return None, None, [None,None,None],None, None,None
    print("outliercount",len(totalcount))
    Apsi.insert(0,[1])
    Akvb.insert(0,[0])
    Akh.insert(0,[0])
    Akvo.insert(0,[0])
    c.insert(0,0)
    us.insert(0,0)

    ##We now have all the terms we need, we just need to reshape and concatenate
    m = columndictionary["max"]+1
    us = us[:m]

    if params["mixs"] == [True,True,True]:
        A,widths = combineAs([m,m,18,1],[1,2],Apsi,Akvb,Akh,Akvo)
        print("ALL ON")
    elif params["mixs"] == [True,False,True]:
        print("ALL KVO, Kh")
        A,widths = combineAs([m,18,1],[1],Apsi,Akh,Akvo)
    elif params["mixs"] == [True,False,False]:
        print("ALL KVO")
        A,widths = combineAs([m,1],[],Apsi,Akvo)
    elif params["mixs"] == [False,False,False]:
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
    if np.isnan(D).any():
        pdb.set_trace()
    print(c[np.isnan(c)])
    prime = np.matmul(j,c)

    if params["debug"]:
        graphError(A,np.concatenate((us,[0]*(A.shape[1]-len(us)))),prime)
        rdivc(D)

    errors = error(A,np.concatenate((us,[0]*(A.shape[1]-len(us)))),prime)
    #calculatediapycnalv(surfaces,prime,params,widths,columndictionary)
    metadata = copy.copy(params)
    metadata["condition"] = condition(D)
    metadata["error"] = np.sum(np.abs(errors[1]))
    surfacer = applyPrime(surfaces,prime,columndictionary,params,widths,mixing=True)
    return surfaces, columndictionary, [VT,D,U],A, errors,metadata

def calculatediapycnalv(surfaces,prime,params,widths,coldict):
    kvbs=np.asarray([])
    kvhs =np.asarray([])
    if params["mixs"] == [True,True,True]:
        kvbs = prime[widths[0]:widths[0]+widths[1]]
        kvhs = prime[widths[0]+widths[1]:widths[0]+widths[1]+widths[2]]
    if params["mixs"] == [True,False,True]:
        kvhs = prime[widths[0]:widths[0]+widths[1]]
    for k in surfaces.keys():
        surfaces[k]["data"]["e"] = np.full(len(surfaces[k]["data"]["psi"]),np.nan)
        if params["lowerbound"] >= k >= params["upperbound"]:
            print("bo")
            for point in range(len(surfaces[k]["data"]["e"])):
                kpvs,ks = ptools.kterms(surfaces,k,point,[1]*4)
                if kpvs.any():
                    rhonot = 1025.0
                    e=0
                    kvhindex = (k-params["upperbound"])/200
                    if kvhs.any() and kvhindex<len(kvhs) :
                        print("params scale: ",params["scalecoeffs"][3])
                        kvh = kvhs[int(kvhindex)]*params["scalecoeffs"][3]
                        e+=kvh*kpvs[2]
                    if kvbs.any():
                        kvb = kvbs[coldict[surfaces[k]["ids"][point]]]
                        e+=kvb*kpvs[1]
                    print(e)
                    surfaces[k]["data"]["e"][point]=e



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
        print(x.shape)
        if np.all(x == 0, axis=0).any():
            print("removing 0 columns")
        if i in totrim:
            x = x[:,~np.all(x == 0, axis=0)]
        print(x.shape)
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
    kvh = prime[2*m:2*m+8]
    kvo =prime[-1]
    fig, (ax1,ax2) = plt.subplots(1,2)
    ax1.set_title("KVB")
    ax1.plot(range(len(kvb)),kvb+kvo)
    ax2.set_title("KVH")
    ax2.plot(range(len(kvh)),kvh)
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


