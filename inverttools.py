from nstools import *
import scipy
import copy

def SVDCalc(VT,D,U,n_elements=False):
    if n_elements:
        D = D[:, :n_elements]
        VT = VT[:n_elements, :]

    B = VT.T.dot(D.T).dot(U.T)
    return B

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

def optimizeA(B,m):
    #justkvb
    bestscale = -100000
    bestsum = 0
    sings = []
    conditions = []
    scales = []
    for i in np.arange(-20,10,.5):
        print("optimize: ",i)
        scale = 10**i
        diags = [1]*m + [1]*m+[1]*18+[scale]
        scalematrix = np.diag(diags)
        A = np.matmul(B,scalematrix)
        A = A / np.linalg.norm(A, axis=1, keepdims=True)
        j,[VT, D, U] = SVDdecomp(A,n_elements=int(A.shape[1]))
        D = 1.0/np.diag(D)
        s = np.sum(D)
        print(s)
        conditions.append(s)
        scales.append(scale)
        if s >bestsum:
            bestscale = scale
            bestsum = s
            sings = D
    plt.plot(range(len(conditions)),conditions)
    plt.show()
    plt.plot(range(len(sings)),sings)
    plt.show()
    return bestscale

def scaleNormSVD(A,m,scale):
    diags = [1]*m + [scale]*m
    scalematrix = np.diag(diags)
    A = np.matmul(A,scalematrix)
    A = A / np.linalg.norm(A, axis=1, keepdims=True)
    j,[VT, D, U] = SVDdecomp(A,n_elements=int(A.shape[1]))
    D = 1.0/np.diag(D)
    s = D[0]/D[-1]
    return s,D

def gradAttack(A,m):
    prevsum = np.inf
    sings = []
    scale = 10**5
    sums = []
    while True:
        sums.append(prevsum)
        print(scale)
        plt.plot(range(len(sings)),sings)
        plt.show()
        forwardSum,forwardD = scaleNormSVD(A,m,scale*(1+10**-2))
        backwardSum,backwardD = scaleNormSVD(A,m,scale/(1+10**-2))
        print(prevsum,forwardSum,backwardSum)
        if forwardSum<backwardSum<prevsum:
            scale = scale*(1+10**-2)
            sings = forwardD
            prevsum=forwardSum
        elif backwardSum<forwardSum<prevsum: 
            scale = scale/(1+10**-2)
            sings = backwardD
            prevsum=backwardSum
        else:
            break
    plt.plot(sums)
    plt.show()
    return scale
    print("WINNER")

        
        


def kterms(surfaces,k,found,debug=False):
    f = gsw.f(surfaces[k]["lats"][found])
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
    kvbscale = 1.0/(5*(10**-5)) #(1/(10**-1.09))
#/630957)#5*(10**2)#5*(10**-5)
    kvoscale = 1.0/(5*(10**-6))#(1/0.001)#5*(10**1)#5*(10**-6)
    khscale = 1.0/500#1/(10**-2.4))
# 10**12#500

    if (np.isnan(isitnan).any()):
        if debug:
            print("pres is nan: ",np.isnan(pres))
            print("hx is nan: ",np.isnan(hx))
            print("hy is nan: ",np.isnan(hy))
            print("x is nan: ",np.isnan(x))
            print("y is nan: ",np.isnan(y))
            print("something here is nan")
        return None, None
    if not (np.isnan(isitnan).any()):
        pvkvb = (d2qdz2+2*(1/1000)*dqdz+(1/(1000**2))*pv)*CKVB
        pvkv0 = d2qdz2
        pvkh = (d2qdx2+d2qdy2)-2*(dqnotdx*dqdx+dqnotdy*dqdy)/pv -f*khpdz
        skvo = -alpha*f*(1/pv)*(dsdz**3)*doublets
        skvb = skvo*CKVB
        skhpart1 = (f/pv)*dsdz*(alphat*(dtdx**2 + dtdy**2)+alphap*(dtdx*dpdx+dtdy*dpdy))
        skhpart2 = (d2sdx2+d2sdy2)-2*(dqnotdx*dsdx + dqnotdy*dsdy)/pv
        skh = skhpart1 + skhpart2
        return (-pvkv0/kvoscale,-pvkvb/kvbscale,-pvkh/khscale),(-skvo/kvoscale,-skvb/kvbscale,-skh/khscale)


 

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
                dqnotdx = surfaces[k]["data"]["dqnotdx"][found]
                dqnotdy = surfaces[k]["data"]["dqnotdy"][found]
                pv = surfaces[k]["data"]["pv"][found]
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

def graphError(b,us,prime):
    b = np.asarray(b)
    us = np.asarray(us)
    prime = np.asarray(prime)
    errorbefore = np.matmul(b,us+prime)
    errorafter = np.matmul(b,us)
    plt.scatter(range(len(errorafter)),errorafter,label="after")
    plt.scatter(range(len(errorbefore)),errorbefore,label="before")
    plt.gca().legend()
    plt.show()
    delta = np.abs(errorbefore)-np.abs(errorafter)
    plt.scatter(range(len(delta[delta<0])),delta[delta<0],c="red")
    plt.scatter(range(len(delta[delta>=0])),delta[delta>=0],c="blue")
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
                kpv,ks = kterms(surfaces,k,found)
                if k>=1000 and kpv and ks:
                    x = surfaces[k]["x"][found]
                    y = surfaces[k]["y"][found]
                    r = np.sqrt(x**2 + y**2)
                    hx = surfaces[k]["data"]["hx"][found]
                    dsdx = surfaces[k]["data"]["dsdx"][found]
                    dsdy = surfaces[k]["data"]["dsdy"][found]
                    hy = surfaces[k]["data"]["hy"][found]
                    f = gsw.f(surfaces[k]["lats"][found])
                    beta = calcBeta(surfaces[k]["lats"][found])
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


def getColumnNumber(eyedict,eyed):
    #assign each eyed a column number 
    if "max" not in eyedict.keys():
        eyedict["max"]=-1
    if eyed not in eyedict.keys():
        eyedict[eyed]=eyedict["max"]+1 
        eyedict["max"] = eyedict["max"]+1
    return eyedict,eyedict[eyed]

def neighborsColumnNumbers(surfaces,k,s,eyedict):
    ## return the column number of each neighbor
    columnnumbers = []
    for l in s:
        eyedict,col = getColumnNumber(eyedict,surfaces[k]["ids"][l])
        columnnumbers.append(col)
    return eyedict,columnnumbers

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

def applyPrime(staggeredsurfaces,prime,coldict,mixing=False):
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
                staggeredsurfaces[k]["data"]["psinew"][i] = staggeredsurfaces[k]["data"]["psiref"][i] + prime[coldict[eyed]]
                staggeredsurfaces[k]["data"]["psisol"][i] = prime[coldict[eyed]]
                if mixing and False:
                    m = coldict["max"]
                    print(len(prime))
                    print(m)
                    staggeredsurfaces[k]["data"]["kvb"][i] = prime[m+coldict[eyed]]
                    staggeredsurfaces[k]["data"]["kvo"][i] = prime[(2*m)+int(k/200.0)]
                #staggeredsurfaces[k]["data"]["psinew"][i] =  prime[coldict[eyed]]
    return staggeredsurfaces


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

def constructBetaRow(surfaces,k,distances,s,columnindexs,Arscale=1/0.05):
    Apsirow = [0]*(max(columnindexs)+1)
    values = [0]*4

    beta = calcBeta(surfaces[k]["lats"][s[0]])
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

    d01 = distances[(s[0],s[1])]
    d02 = distances[(s[0],s[2])]
    d23 = distances[(s[2],s[3])]
    d13 = distances[(s[1],s[3])]

    values[0] = (-1/(2*f*d02))*(dqnotdx-x*beta*pv/(f*r))\
                               + (1/(2*f*d01))*(dqnotdy-y*beta*pv/(f*r))

    values[1] =(-1/(2*f*d13))*(dqnotdx-x*beta*pv/(f*r))\
                              + (-1/(2*f*d01))*(dqnotdy-y*beta*pv/(f*r))

    values[2] = (1/(2*f*d02))*(dqnotdx-x*beta*pv/(f*r))\
                               + (1/(2*f*d23))*(dqnotdy-y*beta*pv/(f*r))
    ## (-1/f)dAr/dx*dQnotdx
    values[3] = (1/(2*f*d13))*(dqnotdx-x*beta*pv/(f*r))\
                              + (-1/(2*f*d23))*(dqnotdy-y*beta*pv/(f*r))

    Apsirow[columnindexs[0]] = values[0]/Arscale
    Apsirow[columnindexs[1]] = values[1]/Arscale
    Apsirow[columnindexs[2]] = values[2]/Arscale
    Apsirow[columnindexs[3]] = values[3]/Arscale
    crow = (-u)*(dqnotdx-x*beta*pv/(f*r))+(-v)*(dqnotdy-y*beta*pv/(f*r))
    return Apsirow,values, crow/Arscale

def constructSalRow(surfaces,k,distances,s,columnindexs,Arscale=1/0.05):
    Apsirow = [0]*(max(columnindexs)+1)
    values = [0]*4
    dsdx = surfaces[k]["data"]["dsdx"][s[0]]
    dsdy = surfaces[k]["data"]["dsdy"][s[0]]
    f = gsw.f(surfaces[k]["lats"][s[0]])
    u = surfaces[k]["data"]["uref"][s[0]] 
    v = surfaces[k]["data"]["vref"][s[0]] 
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

    Apsirow[columnindexs[0]] = values[0]/Arscale
    Apsirow[columnindexs[1]] = values[1]/Arscale
    Apsirow[columnindexs[2]] = values[2]/Arscale
    Apsirow[columnindexs[3]] = values[3]/Arscale
    crow = -((u)*dsdx+(v)*dsdy)/Arscale

    return Apsirow,values,crow


def coupledInvertNoMixing(surfaces,reflevel,neighbors,distances,debug=False,single=False):
    outsurfaces = copy.deepcopy(surfaces)
    Apsi=[]
    Akvb=[]
    Akh=[]
    Akvo=[]
    c=[]
    us = [0]*100000
    ##dictionary of ids to column numbers
    columndictionary = {"max":-1}
    surfaces = applyRefLevel(surfaces)
    alldistances = distances 
    selects = False    
    for k in Bar("adding to matrix: ").iter(neighbors.keys()):
        distances = alldistances[k]
        #for s in Bar('coupled invert: ').iter(neighbors[k]):
        for s in neighbors[k]:
            if not selects:
                selects = surfaces[k]["ids"][s[0]]
            s=np.asarray(s)
            kpv,ks = kterms(surfaces,k,s[0])
            u = surfaces[k]["data"]["uref"][s[0]] 
            v = surfaces[k]["data"]["vref"][s[0]] 
            if (not np.isnan(u) and not np.isnan(v)) and 3600>k>=1000 and kpv and ks and (not single or selects == surfaces[k]["ids"][s[0]] ): #and surfaces[k]["ids"][s[0]] in whitelist :
                ##find/store column indexs for each neighbor
                columndictionary,columnindexs = neighborsColumnNumbers(surfaces,k,s,columndictionary)
                ##this is a shorthad way of converting the NS to an index
                #######PVROW
                betarow,betavals, crow = constructBetaRow(surfaces,k,distances,s,columnindexs,1)

                n = np.linalg.norm(betavals)
                Apsi.append(betarow/n)

                ##make rows that can fit it 
                us[columnindexs[0]] = surfaces[k]["data"]["psiref"][s[0]]
                c.append(crow/n)

                #######SALROW
                ##make rows that can fit it 
                salrow, salvals, crow = constructSalRow(surfaces,k,distances,s,columnindexs,1)
                n = np.linalg.norm(salvals)
                Apsi.append(salrow/n)
                ##im a rascal and this is a shorthad way of converting the NS to an index :P
                 ######SAL Error row
                c.append(crow/n)

    Apsi.insert(0,[1])
    c.insert(0,0)
    us.insert(0,0)
    ##We now have all the terms we need, we just need to reshape and concatenate
    m = columndictionary["max"]+1
    us = us[:m]
    #print(Apsi)

    A = combineAs([m,m,18,1],Apsi)#,Akvb,Akh,Akvo)
    rowCount(A)

        #print(np.matrix(A))
    print("#####A########")
    print(A.shape)
    print("#####Apsi########")
    print(len(Apsi))
    print("#####Akh########")
    print(len(Akh))
    print("#####Ako########")
    print(len(Akvo))
    print("#####c########")
    print(len(c))
    print("#############")

    c = np.matrix.transpose(np.asarray(c))
    us =  np.matrix.transpose(np.asarray(us))
    j,[VT, D, U] = SVDdecomp(A,n_elements=A.shape[1])

    prime = np.matmul(j,c)
    print("singular values: ", 1/(np.diag(D)))
    #def graphError(b,us,prime):
    graphError(A,us,prime)
    #parameterErrors(A,prime,c,D)
    rdivc(D)

    surfaces = applyPrime(surfaces,prime,columndictionary)
    return surfaces, columndictionary, [VT,D,U],A


def coupledInvert(surfaces,reflevel,neighbors,distances,debug=False,single=False,normalize=False):
    outsurfaces = copy.deepcopy(surfaces)
    Apsi=[]
    Akvb=[]
    Akh=[]
    Akvo=[]
    c=[]
    us = [0]*100000
    ##dictionary of ids to column numbers
    columndictionary = {"max":-1}
    surfaces = applyRefLevel(surfaces)
    alldistances = distances 
    selects = False    
    whitelist = generateWhiteList(surfaces,neighbors,3600,1000)
    for k in Bar("adding to matrix: ").iter(neighbors.keys()):
        distances = alldistances[k]
        #for s in Bar('coupled invert: ').iter(neighbors[k]):
        for s in neighbors[k]:
            if not selects:
                selects = surfaces[k]["ids"][s[0]]
            s=np.asarray(s)
            kpv,ks = kterms(surfaces,k,s[0])
            kpv = np.asarray(kpv)
            ks = np.asarray(ks)
            u = surfaces[k]["data"]["uref"][s[0]] 
            v = surfaces[k]["data"]["vref"][s[0]] 
            if (not np.isnan(u) and not np.isnan(v)) and kpv.any() and ks.any() and 3600>k>=1000 and (not single or selects == surfaces[k]["ids"][s[0]] ) and surfaces[k]["ids"][s[0]] in whitelist :
                ##find/store column indexs for each neighbor
                columndictionary,columnindexs = neighborsColumnNumbers(surfaces,k,s,columndictionary)
                ##this is a shorthad way of converting the NS to an index
                #######PVROW
                betarow,betavals, crow = constructBetaRow(surfaces,k,distances,s,columnindexs)
                betavals = np.asarray(betavals)
                n = np.linalg.norm(np.concatenate((np.asarray(betavals),kpv[[True,False,False]])))

                #n = np.linalg.norm(betavals)
                Apsi.append(np.asarray(betarow)/n)

                ##make rows that can fit it 
                Akvbrow = [0]*(columnindexs[0]+1)
                Akvbrow[columnindexs[0]]=kpv[1]/n
                #Akvbrow[columnindexs[0]]=0
                Akvb.append(Akvbrow)

                Akhrow = [0]*(int(k/200))
                Akhrow[int(k/200-1)] = kpv[2]/n
                Akh.append(Akhrow)


                Akvo.append([kpv[0]/n])
                #Akvo.append([kpv[0]/n])
                us[columnindexs[0]] = surfaces[k]["data"]["psiref"][s[0]]
                c.append(crow/n)

                #######SALROW
                ##make rows that can fit it 
                salrow, salvals, crow = constructSalRow(surfaces,k,distances,s,columnindexs)
                salvals = np.asarray(salvals)
                n = np.linalg.norm(np.concatenate((np.asarray(salvals),ks[[True,False,False]])))
                #n = np.linalg.norm(salvals)
                Apsi.append(np.asarray(salrow)/n)
                ##im a rascal and this is a shorthad way of converting the NS to an index :P
                Akvbrow = [0]*(columnindexs[0]+1)
                Akvbrow[columnindexs[0]]=ks[1]/n

                Akhrow = [0]*(int(k/200))
                Akhrow[int(k/200-1)] = ks[2]/n

                Akvb.append(Akvbrow)
                Akh.append(Akhrow)

                Akvo.append([ks[0]/n])

                ######SAL Error row
                c.append(crow/n)

    Apsi.insert(0,[1])
    Akvb.insert(0,[0])
    Akh.insert(0,[0])
    Akvo.insert(0,[0])
    c.insert(0,0)
    us.insert(0,0)
    ##We now have all the terms we need, we just need to reshape and concatenate
    m = columndictionary["max"]+1
    us = us[:m]

    #print(Apsi)

    #A = combineAs([m,m,18,1],Apsi,Akvb,Akh,Akvo)
    #A = combineAs([m,m,18,1],Apsi,Akvb,Akh,Akvo)
    A = combineAs([m,1],Apsi,Akvo)
    print("column count: ",np.count_nonzero(A, axis=0))
    rowCount(A)
    #print("Scale: ",optimizeA(A,m))

        #print(np.matrix(A))
    print("#####A########")
    print(A.shape)
    print("#####Apsi########")
    print(len(Apsi))
    print("#####Akh########")
    print(len(Akh))
    print("#####Ako########")
    print(len(Akvo))
    print("#####c########")
    print(len(c))
    print("#############")
    c = np.matrix.transpose(np.asarray(c))
    us =  np.matrix.transpose(np.asarray(us))
    j,[VT, D, U] = SVDdecomp(A,n_elements=int(A.shape[1]))
    prime = np.matmul(j,c)
    print(prime[-1])
    #graphError(A,np.concatenate((us,[0]*(len(us)+19))),prime)
    graphError(A,np.concatenate((us,[0]*(A.shape[1]-len(us)))),prime)
    graphKs(prime,m)
    print("singular values: ", 1/(np.diag(D)))

    rdivc(D)

    surfaces = applyPrime(surfaces,prime,columndictionary,mixing=True)
    return surfaces, columndictionary, [VT,D,U],A



def rowCount(A):
    number = []
    for i in range(A.shape[0]):
        number.append(np.count_nonzero(A[i]))
    print("in each row: ",np.unique(number))



def rdivc(D):
    d = D.diagonal()
    d = 1.0/d
    plt.title("Singular Values")
    plt.scatter(range(len(d[:])),(d[:]))
    plt.show()
    plt.title("Largest Singular Values Divided by Smallest")
    plt.scatter(range(len(d[:])),(d[0]/d[:]))
    plt.show()


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

def combineAs(maxlengths,*argv):
    new = []
    for i in range(len(argv)):
        x  = np.asarray(rectangularize(argv[i],maxlengths[i]))
        print(x.shape)
        if np.all(x == 0, axis=0).any():
            print("removing 0 columns")
        x = x[:,~np.all(x == 0, axis=0)]
        new.append(x)
    return concat(new)

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

     
def invert(kind,surfaces,neighbors=None,distances=None,reflevel=1000,debug=False):
    if kind == "simple":
        return simpleInvert(surfaces,reflevel,debug)
    if kind == "simplesalt":
        return simplesaltInvert(surfaces,reflevel,debug)
    if kind == "complexsalt":
        return complexSaltInvert(surfaces,reflevel,debug)
    if kind == "complex":
        return complexInvert(surfaces,reflevel,debug)
    if kind == "coupled":
        return coupledInvert(surfaces,reflevel,neighbors,distances,debug,single=False)
    if kind == "couplednomix":
        return coupledInvertNoMixing(surfaces,reflevel,neighbors,distances,debug,single=False)
    else:
        print("Sorry I have no clue what inversion you are talking about")


