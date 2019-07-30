from nstools import *
import copy

def SVDdecomp(A,n_elements=2):
    U, s, VT = svd(A,full_matrices=True)
    # reciprocals of s
    d = 1.0 / s
    # create m x n D matrix
    D = np.zeros(A.shape)
    # populate D with n x n diagonal matrix
    D[:A.shape[1], :A.shape[1]] = np.diag(d)
    D = D[:, :n_elements]
    VT = VT[:n_elements, :]
    # calculate pseudoinverse
    B = VT.T.dot(D.T).dot(U.T)
    return B

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
                beta = calcBeta(surfaces[k]["lats"][found])

                if debug and (np.isnan(hx) or np.isnan(hy) or np.isnan(x) or np.isnan(y)):
                    print("pres is nan: ",np.isnan(pres))
                    print("hx is nan: ",np.isnan(hx))
                    print("hy is nan: ",np.isnan(hy))
                    print("x is nan: ",np.isnan(x))
                    print("y is nan: ",np.isnan(y))
                    print("something here is nan")
                if k>=1000 and not(np.isnan(hx) or np.isnan(hy) or np.isnan(x) or np.isnan(y)):
                    bvec = (hx+(beta*x)/(f*r),hy+(beta*y)/(f*r))
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
            j = SVDdecomp(b,n_elements=2)
            prime = np.matmul(j,c)
            b = b.T
            errorbefore = []
            for i in range(len(b[0])):
                errorbefore.append(b[0][i]*(us[0][i]+prime[0]) + b[1][i]*(us[1][i]+prime[1]))
            errorafter = []
            for i in range(len(b[0])):
                errorafter.append(b[0][i]*(us[0][i]) + b[1][i]*(us[1][i]))
            #plt.plot(errorbefore,label="after")
            #plt.plot(errorafter,label="before")
            #plt.gca().legend()
            #plt.show()
            R = np.matmul((np.matmul(b.T,prime)-c),(np.matmul(b.T,prime)-c).T)
            #delta = np.sqrt(R/(b.shape[1]-2))
            #for i in b.shape[0]:
                #lon.alg.solv(np.matmul(b.T,b)[i][i]

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

def calcBeta(lat):
    omega =  (7.2921 * 10**-5)
    a = 6.357 * (10**6)
    return (2*omega*np.cos(lat))/a

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
                pres = surfaces[k]["data"]["pres"][found]
                f = gsw.f(surfaces[k]["lats"][found])
                beta = calcBeta(surfaces[k]["lats"][found])

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
            j = SVDdecomp(b,n_elements=3)
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

def graphError(b,us,prime):
    errorbefore = []
    errorafter = []
    for i in range(len(b)):
        errorafter.append(np.dot(b[i],(us[i]+prime)))
        errorbefore.append(np.dot(b[i],(us[i])))
    print("before: ",np.linalg.norm(errorbefore))
    print("after:",np.linalg.norm(errorafter))
    plt.plot(errorafter,label="after")
    plt.plot(errorbefore,label="before")
    plt.gca().legend()
    plt.show()



def norm(v):
    return v/np.linalg.norm(v)

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
                beta = calcBeta(surfaces[k]["lats"][found])
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


def complexInvert(surfaces,reflevel=1000,debug=False):
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
                dqdz = surfaces[k]["data"]["dqdz"][found] 
                d2qdz2 = surfaces[k]["data"]["d2qdz2"][found] 
                dpdx = surfaces[k]["data"]["dpdx"][found] 
                dpdy = surfaces[k]["data"]["dpdy"][found] 
                dqdx = surfaces[k]["data"]["dqdx"][found] 
                dqdy = surfaces[k]["data"]["dqdy"][found] 
                d2qdx2 = surfaces[k]["data"]["d2qdx2"][found] 
                d2qdy2 = surfaces[k]["data"]["d2qdy2"][found] 
                khpdz = surfaces[k]["data"]["khpdz"][found] 
                alphat = dalphadtheta+2*(alpha/betaTherm)*dalphads-(alpha**2/betaTherm**2)*dbetads
                alphap = dalphadp -(alpha/betaTherm)*dbetadp
                pv =  surfaces[k]["data"]["pv"][found] 
                doublets =  surfaces[k]["data"]["d2thetads2"][found] 
                CKVB =  surfaces[k]["data"]["CKVB"][found] 
                f = gsw.f(surfaces[k]["lats"][found])
                beta = calcBeta(surfaces[k]["lats"][found])
                isitnan = [alpha,betaTherm,dsdz,hx,hy,dsdx,dsdy,pres,d2sdx2,d2sdy2,\
                          dalphadtheta,dalphads,dalphadp,dbetadp,dbetads,dtdx,dtdy,\
                          dqnotdx,dqnotdy,dpdx,dpdy,alphat,alphap,pv,doublets,CKVB,\
                          beta,d2qdx2,d2qdy2,khpdz]

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
                    pvkvb = (d2qdz2+2*(1/1000)*dqdz+(1/(1000**2))*pv)*CKVB
                    pvkv0 = d2qdz2
                    pvkh = (d2qdx2+d2qdy2)-2*(dqnotdx*dqdx+dqnotdy*dqdy)/pv -f*khpdz
                    us.append((u,v,0,0,0))
                    us.append((u,v,0,0,0))
                    ######################################
                    ###Potential vorticity row
                    #left side
                    pvvec=((hx+(beta*x)/(f*r)),(hy+(beta*y)/(f*r)),-pvkvb,-pvkv0,-pvkh)
                    pvnorm=np.linalg.norm(pvvec)
                    b.append(pvvec/pvnorm)
                    #right side
                    c.append(np.dot((-u,-v),pvvec[0:2]/pvnorm))
                    ######################################
                    ###Salinity row
                    #left side
                    skvo = -alpha*f*(1/pv)*(dsdz**3)*doublets
                    skvb = skvo*CKVB

                    skhpart1 = (f/pv)*dsdz*(alphat*(dtdx**2 + dtdy**2)+alphap*(dtdx*dpdx+dtdy*dpdy))
                    skhpart2 = (d2sdx2+d2sdy2)-2*(dqnotdx*dsdx + dqnotdy*dsdy)/pv
                    skh = skhpart1 + skhpart2

                    svec=(dsdx,dsdy,-skvb,-skvo,-skh)
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


   
def invert(kind,surfaces,reflevel=1000,debug=False):
    if kind == "simple":
        return simpleInvert(surfaces,reflevel,debug)
    if kind == "simplesalt":
        return simplesaltInvert(surfaces,reflevel,debug)
    if kind == "complexsalt":
        return complexSaltInvert(surfaces,reflevel,debug)
    if kind == "complex":
        return complexInvert(surfaces,reflevel,debug)
    else:
        print("Sorry I have no clue what inversion you are talking about")


