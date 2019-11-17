import scipy
import scipy.io as sio
import numpy as np
from progress.bar import Bar
from profile import Profile
import pickle
import gsw
import nstools
import graph


#generate a unique id
def idgenerator():
    if not hasattr(idgenerator,"id"):
        idgenerator.id=0
    idgenerator.id = idgenerator.id+1
    return idgenerator.id

def neutralSurfaces(fname):
    ctddata = sio.loadmat(fname)
    return np.asarray(ctddata["P_gref"])[:,0]
def nepbCTDExtract(fname,savepath):
    ctddata = sio.loadmat(fname)
    lats = np.asarray(ctddata["lat"])[0]
    lons = np.asarray(ctddata["lon"])[0]
    pres = np.asarray(ctddata["Pint"]).T[0]
    sals = np.asarray(ctddata["Sint"]).T
    thetas = np.asarray(ctddata["Tint"]).T
    CT = np.asarray(ctddata["CT"]).T
    SA = np.asarray(ctddata["Sint_abs"]).T
    nspres = np.asarray(ctddata["P_gamma"]).T
    PV = np.asarray(ctddata["PV"]).T
    ns = np.asarray(ctddata["P_gref"])
    profiles = []
    for p in Bar("profile").iter(range(int(len(lats)))):
        data = {}
        knownns = {}
        knownpv = {}
        if lats[p]>20 and lons[p]<0 or lons[p] > 170:
            data["lat"]=lats[p]
            data["lon"]=lons[p]
            data["temp"]=[]
            data["sal"]=[]
            data["pres"]=[]
            for j in range(len(pres)):
                if ~np.isnan(sals[p][j]) and ~np.isnan(thetas[p][j]) and ~np.isnan(pres[j])\
                        and abs(pres[j]) < 10000 and abs(thetas[p][j])<30 and abs(sals[p][j])<50:
                    insitutemp = thetas[p][j]
                    practicalsal = sals[p][j]
                    singlepres = abs(pres[j])
                    abssal = gsw.SA_from_SP(practicalsal,singlepres,lons[p],lats[p])
                    conservativetemp = gsw.CT_from_t(abssal,insitutemp,singlepres)
                    if abs(abssal-SA[p][j]) > 0.001:
                        print("saldiff = ",abssal-SA[p][j])
                    if abs(conservativetemp-CT[p][j]) > 0.001:
                        print("tempdiff = ",conservativetemp-CT[p][j])
                    data["temp"].append(insitutemp)
                    data["sal"].append(practicalsal)
                    data["pres"].append(singlepres)
            for a in range(len(ns)):
                knownns[ns[a][0]]=nspres[p][a]
            for a in range(len(ns)):
                knownpv[ns[a][0]]=PV[p][a]
            data["knownns"]=knownns
            data["knownpv"]=knownpv

            if len(data["pres"])>4 and max(data["pres"])>1500:
                eyed=int(p)
                prof=Profile(eyed,data)
                profiles.append(prof)
                
    with open(savepath, 'wb') as outfile:
        pickle.dump(profiles, outfile)

def extractPointSurfaces(fname,fields={}):
    if len(fields.keys())<1:
        fields = {"lat":"lat","lon":"lon","pres":"P_gamma",\
            "s":"NS_S","ct":"NS_CT","pv":"PV","psi":"accpo","dthetads":"dCTdS"}
    ctddata = sio.loadmat(fname)
    lats = np.asarray(ctddata[fields["lat"]])[0]
    lons = np.asarray(ctddata[fields["lon"]])[0]
    nspres = np.asarray(ctddata[fields["pres"]]).T
    PV = np.asarray(ctddata[fields["pv"]]).T
    NS_CT = np.asarray(ctddata[fields["ct"]]).T
    dCTdS = np.asarray(ctddata[fields["dthetads"]]).T
    NS_S = np.asarray(ctddata[fields["s"]]).T
    ACCPO = np.asarray(ctddata[fields["psi"]]).T
    ns = np.asarray(ctddata["P_gref"])
    profiles = []
    surfaces = {}

    for j in Bar("surface").iter(range(len(ns))):
        tempSurf = nstools.emptySurface()
        tempSurf["data"]["psi"] = []
        tempSurf["data"]["dthetads"] = []
        tempSurf["data"]["dalphadp"] = []
        tempSurf["data"]["dalphadtheta"] = []
        for p in range(len(lats)):
            if lats[p]>20 and lons[p]<0 or lons[p] > 170:
                tempSurf["lats"].append(lats[p])
                tempSurf["lons"].append(lons[p])
                tempSurf["ids"].append(p)
                tempSurf["data"]["pres"].append(nspres[p][j])
                tempSurf["data"]["t"].append(NS_CT[p][j])
                tempSurf["data"]["s"].append(NS_S[p][j])
                if ns[j][0] <1900:
                    tempSurf["data"]["psi"].append(ACCPO[p][j])
                else:
                    tempSurf["data"]["psi"].append(np.nan)

                tempSurf["data"]["pv"].append(PV[p][j])
                tempSurf["data"]["n^2"].append(np.nan)
                tempSurf["data"]["dalphadp"].append(np.nan)
                tempSurf["data"]["dalphadtheta"].append(np.nan)
                tempSurf["data"]["alpha"].append(gsw.alpha(NS_S[p][j],NS_CT[p][j],nspres[p][j]))
                tempSurf["data"]["beta"].append(gsw.beta(NS_S[p][j],NS_CT[p][j],nspres[p][j]))
                tempSurf["data"]["dthetads"].append(dCTdS[p][j])
        surfaces[ns[j][0]] = tempSurf
    return surfaces

def floatExtract(fname):
    fields = {"lat":"F_NS_lat","lon":"F_NS_lon","pres":"F_P_gamma",\
            "s":"F_NS_S","ct":"F_NS_CT","pv":"F_PV","psi":"F_accpo","dthetads":"F_dCTdS"}
    return extractPointSurfaces(fname,fields=fields)



def nepbCTDExtractInterpSurfaces(fname,calcDeriv = False):
    ctddata = sio.loadmat(fname)
    ##note: any x,y gradients wont be equivalent because of grid,
    ##			advisable to look at total gradient probably
    quantmap = {"CT_s":"t","S_s":"s","dQdz_s":"dqdz","dSdz_s":"dsdz",\
    "d2CTdS2_s":"d2thetads2","alpha_s":"alpha","beta_s":"beta",\
    "P_s":"pres","Q_s":"pv","d2Qdx2_s":"d2qdx2","d2Qdy2_s":"d2qdy2",\
    "d2Qdz2_s":"d2qdz2","d2Sdx2_s":"d2sdx2","d2Sdy2_s":"d2sdy2",\
    "aT":"dalphadtheta","aP":"dalphadp","A_s":"psi",\
    "lat_field":"khp","dCTdS_s":"dthetads","Sx":"dsdx",\
    "Sy":"dsdy","CTx":"dtdx","CTy":"dtdy","dCTdz_s":"dtdz"}

    latlist = range(20,60,2)
    lonlist = list(range(170,180,2))
    lonlist = lonlist+list(range(-180,-120,2))

    ns = np.asarray(ctddata["P_gref"])

    surfaces = {}

    for field in ctddata.keys():
        if field in list(quantmap.keys())+["F_Ssmooth","F_CTsmooth"]:
            if ctddata[field].shape == (20, 35, 30):
                ctddata[field] = np.transpose(ctddata[field],(2,0,1))
            if ctddata[field].shape == (41, 71, 30):
                ctddata[field] = np.transpose(ctddata[field],(2,0,1))

    for k in Bar("surface").iter(range(len(ns))):
        tempSurf = nstools.emptySurface()
        for j in quantmap.values():
            tempSurf["data"][j]=[]
        tempSurf["data"]["dsdx"]=[]
        tempSurf["data"]["dsdy"]=[]
        tempSurf["data"]["dtdx"]=[]
        tempSurf["data"]["dtdy"]=[]

        for j in range(len(latlist)):
            for l in range(len(lonlist)):
                tempSurf["lats"].append(latlist[j])
                tempSurf["lons"].append(lonlist[l])
                tempSurf["ids"].append(j*30+l)
                for field in ctddata.keys():
                    if field in quantmap.keys():
                        if calcDeriv and (field not in ["CTx","CTy","Sx","Sy"]):
                            tempSurf["data"][quantmap[field]].append(ctddata[field][k][j][l])
                        elif not calcDeriv:
                            tempSurf["data"][quantmap[field]].append(ctddata[field][k][j][l])

                    if calcDeriv and field == "F_Ssmooth":

                        if j == len(latlist)-1 or l == len(lonlist)-1: 
                            ydist = ((latlist[j] - latlist[j-1])/2.0)*111.0*1000.0
                            xdist = ydist *np.cos(np.deg2rad(latlist[j]+1))
                        else: 
                            ydist = ((latlist[j+1] - latlist[j])/2.0)*111.0*1000.0
                            xdist =ydist * np.cos(np.deg2rad(latlist[j]+1))

                        dx = ctddata[field][k][j*2][l*2+1] -  ctddata[field][k][j*2][l*2]
                        dx += ctddata[field][k][j*2+1][l*2+1] -  ctddata[field][k][j*2+1][l*2]
                        dy = ctddata[field][k][j*2+1][l*2] -  ctddata[field][k][j*2][l*2]
                        dy += ctddata[field][k][j*2+1][l*2+1] -  ctddata[field][k][j*2][l*2+1]

                        tempSurf["data"]["dsdx"].append(dx/(2*xdist))
                        tempSurf["data"]["dsdy"].append(dy/(ydist))


                    if calcDeriv and field == "F_CTsmooth":

                        if j == len(latlist)-1 or l == len(lonlist)-1: 
                            ydist = ((latlist[j] - latlist[j-1])/2.0)*111.0*1000.0
                            xdist = ydist *np.cos(np.deg2rad(latlist[j]+1))
                        else: 
                            ydist = ((latlist[j+1] - latlist[j])/2.0)*111.0*1000.0
                            xdist =ydist * np.cos(np.deg2rad(latlist[j]+1))

                        dx = ctddata[field][k][j*2][l*2+1] -  ctddata[field][k][j*2][l*2]
                        dx += ctddata[field][k][j*2+1][l*2+1] -  ctddata[field][k][j*2+1][l*2]
                        dy = ctddata[field][k][j*2+1][l*2] -  ctddata[field][k][j*2][l*2]
                        dy += ctddata[field][k][j*2+1][l*2+1] -  ctddata[field][k][j*2][l*2+1]

                        tempSurf["data"]["dtdx"].append(dx/(2*xdist))
                        tempSurf["data"]["dtdy"].append(dy/(ydist))



        surfaces[ns[k][0]] = tempSurf

    return surfaces


#nepbCTDExtract("data/newnepbdata.mat","data/nepbctdprofiles.pickle")
#nepbCTDExtract("data/newnepbdata.mat","data/nepbctdprofiles.pickle")

