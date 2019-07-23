from nstools import *

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
    omega =  (7.2921 * 10**-5)
    a = 6.357 * (10**6)
    for index in  Bar('Inverting: ').iter(range(len(surfaces[reflevel]["x"]))):
        eyed = int(surfaces[reflevel]["ids"][index])
        #print("id: ",eyed)
        ns = []
        us = [[],[]]
        b = [[],[]]
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
                    #surfaces[k]["data"]["uabs"][found]=0
                    b[0].append(hx+(beta*x)/(f*r))
                    b[1].append(hy+(beta*y)/(f*r))
                    u = (surfaces[k]["data"]["u"][found] - surfaces[reflevel]["data"]["u"][index])
                    v = (surfaces[k]["data"]["v"][found] - surfaces[reflevel]["data"]["v"][index])
                    us[0].append(u)
                    us[1].append(v)
                    c.append((-u)*hx + (-v)*hy - (beta/f)*(u*x+v*y)/(r))

                if debug:
                    print("k: ",k)
                    print("hx: ",hx)
                    print("hy: ",hy)
                    print("beta: ",beta)
                    print("x: ",x)
                    print("y: ",y)
                    print("f: ",f)
                    print("r: ",r)
                    print("c: ",u*hx + v*hy - (beta/f)*(-u*x-v*y)/(r))
        if len(b[0])>0:
            b = np.matrix.transpose(np.asarray(b))
            c = np.matrix.transpose(np.asarray(c))
            us = np.asarray(us)
            #j = np.linalg.inv(np.matmul(np.matrix.transpose(b),b))
            #j = SVDdecomp(b,n_elements=4)
            j = SVDdecomp(b,n_elements=1)
            prime = np.matmul(j,c)
            b = b.T
            error = []
            for i in range(len(b[0])):
                error.append(b[0][i]*(us[0][i]+prime[0]) + b[1][i]*(us[1][i]+prime[1]))
            #print("corrected: ",np.mean(error),",uncorrected: ",np.mean(c))
            #plt.plot(error,range(len(error)),label="with correction")
            #plt.plot(-c,range(len(error)),label= "unchanged")
            #plt.gca().invert_yaxis()
            #plt.gca().legend()
            #plt.show()

            if debug:
                print("######BBBBBBBBBBBB###############")
                print("b: ",b)
                print("c: ",c)
                print("j: ",j)
                print("j , b: ",j.shape,b.shape)
                print("prime: ",prime)
                print("########################")
            for i in range(len(ns)):
                surfaces[ns[i][0]]["data"]["uprime"][ns[i][1]] = prime[0] -surfaces[reflevel]["data"]["u"][index]
                surfaces[ns[i][0]]["data"]["uabs"][ns[i][1]] = prime[0] + surfaces[ns[i][0]]["data"]["u"][ns[i][1]]-surfaces[reflevel]["data"]["u"][index]
                surfaces[ns[i][0]]["data"]["u"][ns[i][1]] = surfaces[ns[i][0]]["data"]["u"][ns[i][1]]-surfaces[reflevel]["data"]["u"][index]
                surfaces[ns[i][0]]["data"]["vprime"][ns[i][1]] = prime[1]-surfaces[reflevel]["data"]["v"][index]
                surfaces[ns[i][0]]["data"]["vabs"][ns[i][1]] = prime[1] + surfaces[ns[i][0]]["data"]["v"][ns[i][1]]-surfaces[reflevel]["data"]["v"][index]
                surfaces[ns[i][0]]["data"]["v"][ns[i][1]] = surfaces[ns[i][0]]["data"]["v"][ns[i][1]]-surfaces[reflevel]["data"]["v"][index]

    return surfaces

def calcBeta(lat):
    omega =  (7.2921 * 10**-5)
    a = 6.357 * (10**6)
    return (2*omega*np.cos(lat))/a

def simplesaltInvert(surfaces,reflevel=1000,debug=False):
    for index in  Bar('Inverting: ').iter(range(len(surfaces[reflevel]["x"]))):
        eyed = int(surfaces[reflevel]["ids"][index])
        #print("id: ",eyed)
        ns = []
        us = [[],[]]
        b = [[],[]]
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
                if k>=1000 and not(np.isnan(hx) or np.isnan(hy) or np.isnan(x) or np.isnan(y)):
                    #surfaces[k]["data"]["uabs"][found]=0
                    b[0].append(hx+(beta*x)/(f*r)+dsdx)
                    b[1].append(hy+(beta*y)/(f*r)+dsdy)
                    u = (surfaces[k]["data"]["u"][found] - surfaces[reflevel]["data"]["u"][index])
                    v = (surfaces[k]["data"]["v"][found] - surfaces[reflevel]["data"]["v"][index])
                    us[0].append(u)
                    us[1].append(v)
                    c.append((-u)*hx + (-v)*hy - (beta/f)*(u*x+v*y)/(r) + dsdx*u + dsdy*v)

                if debug:
                    print("k: ",k)
                    print("hx: ",hx)
                    print("hy: ",hy)
                    print("beta: ",beta)
                    print("x: ",x)
                    print("y: ",y)
                    print("f: ",f)
                    print("r: ",r)
                    print("c: ",u*hx + v*hy - (beta/f)*(-u*x-v*y)/(r))
        if len(b[0])>0:
            b = np.matrix.transpose(np.asarray(b))
            c = np.matrix.transpose(np.asarray(c))
            us = np.asarray(us)
            #j = np.linalg.inv(np.matmul(np.matrix.transpose(b),b))
            #j = SVDdecomp(b,n_elements=4)
            j = SVDdecomp(b,n_elements=1)
            prime = np.matmul(j,c)
            b = b.T
            error = []
            for i in range(len(b[0])):
                error.append(b[0][i]*(us[0][i]+prime[0]) + b[1][i]*(us[1][i]+prime[1]))
            #print("corrected: ",np.mean(error),",uncorrected: ",np.mean(c))
            #plt.plot(error,range(len(error)),label="with correction")
            #plt.plot(-c,range(len(error)),label= "unchanged")
            #plt.gca().invert_yaxis()
            #plt.gca().legend()
            #plt.show()

            if debug:
                print("######BBBBBBBBBBBB###############")
                print("b: ",b)
                print("c: ",c)
                print("j: ",j)
                print("j , b: ",j.shape,b.shape)
                print("prime: ",prime)
                print("########################")
            for i in range(len(ns)):
                surfaces[ns[i][0]]["data"]["uprime"][ns[i][1]] = prime[0] -surfaces[reflevel]["data"]["u"][index]
                surfaces[ns[i][0]]["data"]["uabs"][ns[i][1]] = prime[0] + surfaces[ns[i][0]]["data"]["u"][ns[i][1]]-surfaces[reflevel]["data"]["u"][index]
                surfaces[ns[i][0]]["data"]["u"][ns[i][1]] = surfaces[ns[i][0]]["data"]["u"][ns[i][1]]-surfaces[reflevel]["data"]["u"][index]
                surfaces[ns[i][0]]["data"]["vprime"][ns[i][1]] = prime[1]-surfaces[reflevel]["data"]["v"][index]
                surfaces[ns[i][0]]["data"]["vabs"][ns[i][1]] = prime[1] + surfaces[ns[i][0]]["data"]["v"][ns[i][1]]-surfaces[reflevel]["data"]["v"][index]
                surfaces[ns[i][0]]["data"]["v"][ns[i][1]] = surfaces[ns[i][0]]["data"]["v"][ns[i][1]]-surfaces[reflevel]["data"]["v"][index]
    return surfaces

def saltInvert(surfaces,reflevel=1000,debug=False):
    for index in  Bar('Inverting: ').iter(range(len(surfaces[reflevel]["x"]))):
        eyed = int(surfaces[reflevel]["ids"][index])
        #print("id: ",eyed)
        ns = []
        us = [[],[]]
        b = [[],[]]
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
                if k>=1000 and not(np.isnan(hx) or np.isnan(hy) or np.isnan(x) or np.isnan(y)):
                    #surfaces[k]["data"]["uabs"][found]=0
                    b[0].append(hx+(beta*x)/(f*r)+dsdx)
                    b[1].append(hy+(beta*y)/(f*r)+dsdy)
                    u = (surfaces[k]["data"]["u"][found] - surfaces[reflevel]["data"]["u"][index])
                    v = (surfaces[k]["data"]["v"][found] - surfaces[reflevel]["data"]["v"][index])
                    us[0].append(u)
                    us[1].append(v)
                    c.append((-u)*hx + (-v)*hy - (beta/f)*(u*x+v*y)/(r) + dsdx*u + dsdy*v)

                if debug:
                    print("k: ",k)
                    print("hx: ",hx)
                    print("hy: ",hy)
                    print("beta: ",beta)
                    print("x: ",x)
                    print("y: ",y)
                    print("f: ",f)
                    print("r: ",r)
                    print("c: ",u*hx + v*hy - (beta/f)*(-u*x-v*y)/(r))
        if len(b[0])>0:
            b = np.matrix.transpose(np.asarray(b))
            c = np.matrix.transpose(np.asarray(c))
            us = np.asarray(us)
            #j = np.linalg.inv(np.matmul(np.matrix.transpose(b),b))
            #j = SVDdecomp(b,n_elements=4)
            j = SVDdecomp(b,n_elements=1)
            prime = np.matmul(j,c)
            b = b.T
            error = []
            for i in range(len(b[0])):
                error.append(b[0][i]*(us[0][i]+prime[0]) + b[1][i]*(us[1][i]+prime[1]))
            #print("corrected: ",np.mean(error),",uncorrected: ",np.mean(c))
            #plt.plot(error,range(len(error)),label="with correction")
            #plt.plot(-c,range(len(error)),label= "unchanged")
            #plt.gca().invert_yaxis()
            #plt.gca().legend()
            #plt.show()

            if debug:
                print("######BBBBBBBBBBBB###############")
                print("b: ",b)
                print("c: ",c)
                print("j: ",j)
                print("j , b: ",j.shape,b.shape)
                print("prime: ",prime)
                print("########################")
            for i in range(len(ns)):
                surfaces[ns[i][0]]["data"]["uprime"][ns[i][1]] = prime[0] -surfaces[reflevel]["data"]["u"][index]
                surfaces[ns[i][0]]["data"]["uabs"][ns[i][1]] = prime[0] + surfaces[ns[i][0]]["data"]["u"][ns[i][1]]-surfaces[reflevel]["data"]["u"][index]
                surfaces[ns[i][0]]["data"]["u"][ns[i][1]] = surfaces[ns[i][0]]["data"]["u"][ns[i][1]]-surfaces[reflevel]["data"]["u"][index]
                surfaces[ns[i][0]]["data"]["vprime"][ns[i][1]] = prime[1]-surfaces[reflevel]["data"]["v"][index]
                surfaces[ns[i][0]]["data"]["vabs"][ns[i][1]] = prime[1] + surfaces[ns[i][0]]["data"]["v"][ns[i][1]]-surfaces[reflevel]["data"]["v"][index]
                surfaces[ns[i][0]]["data"]["v"][ns[i][1]] = surfaces[ns[i][0]]["data"]["v"][ns[i][1]]-surfaces[reflevel]["data"]["v"][index]
    return surfaces



def invert(kind,surfaces,reflevel=1000,debug=False):
    if kind == "simple":
        return simpleInvert(surfaces,reflevel,debug)
    if kind == "simplesalt":
        return simplesaltInvert(surfaces,reflevel,debug)


