# -*- coding: utf-8 -*-
"""
Created on Sat Oct 15 17:28:01 2016

@author: Sebasijan Mrak, Greg Starr
smrak@bu.edu
"""
import numpy as np
from datetime import datetime, timedelta
from scipy import interpolate
from pymap3d import ecef2geodetic,ecef2aer,aer2geodetic,geodetic2ecef
import xarray, os
import matplotlib.pyplot as plt
from pandas import Timestamp
import georinex as gr
from pyGnss import gnssUtils as uf
from pyGnss import scintillation
from glob import glob
from scipy.optimize import least_squares

#constnats for GPS
f1 = 1575420000
f2 = 1227600000
f5 = 1176450000
c0 = 299792458
freq = {'1': f1, '2': f2, '5': f5}

def getSatBias(fn, sv=None):
    import os
    if os.path.isdir(fn):
        fn = glob(fn + '*.**i')
        if len(fn) >= 1:
            fn = fn[0] 
        else: 
            raise('No satbias files found in the folder')
            return 0
    i = 0
    svbias = {}
    with open(fn, 'r') as f:
        while 1:
            line = f.readline()
            if 'DIFFERENTIAL CODE BIASES' in line:
                while i < 32:
                    l = f.readline().split()[:2]
                    svbias['G'+l[0]] = float(l[1])
                    i += 1
                f.close()
                break
    if sv is None:
        return svbias
    else:
        return svbias[sv]

# %% TEC
def getPRNSlantTEC(P1, P2, units='m'):
    """
    Sebsatijan Mrak
    Function returns slant TEC in TECU units. Input data are PRN information
    at two frequences, where f1 and f2 are difined as global variables. It assumes,
    that you use L1 and L2 frquencies. In case of different GNSS constellation or
    use of L5, correct the indexes. 
    Default config. assumes PRN distance in meters [m], otherwise, fulfill the 
    function parameter 'unit' to correct the units.
    Output units are by default in meters.
    """     
    if units == 'm':
        sTEC = ((1/40.3) * (( pow(f2, 2) * pow(f1, 2) ) / 
                (pow(f1, 2) - pow(f2, 2))) * (P2 - P1)) / pow(10,16)
    elif units == 'rad':
        sTEC = ((c0/(40.3*2*np.pi)) * (( pow(f2, 2) * pow(f1, 2) ) / 
                (pow(f1, 2) - pow(f2, 2))) * (P2/f2- P1/f1)) / pow(10,16)
            
    elif units == 'cycle':
        sTEC = ((c0/(40.3)) * (( pow(f2, 2) * pow(f1, 2) ) / 
                (pow(f1, 2) - pow(f2, 2))) * (P2/f2 - P1/f1)) / pow(10,16)        
        
    return sTEC
    
def getPSlantTEC(L1, L2, units = 'cycle'):
    """
    Sebsatijan Mrak
    Function returns slant TEC in TECU units. Input data are phase information
    at two frequences, where f1 and f2 are difined as global variables. It assumes,
    that you use L1 and L2 frquencies. In case of different GNSS constellation or
    use of L5, correct the indexes. 
    Default config. assumes phase information in cycles [cycle], otherwise, 
    fulfill the function parameter 'unit' to correct the units. 
    Output units are by default in meters.
    
    Use only if there is no cycle slips in raw phase file!
    """
    if units == 'cycle':
        sTEC = ((c0/40.3) * (( pow(f2, 2) * pow(f1, 2) ) / 
                (pow(f1, 2) - pow(f2, 2))) * (L1/f1 - L2/f2)) / pow(10,16)
    elif units == 'rad':
        sTEC = ((c0/(40.3*2*np.pi)) * (( pow(f2, 2) * pow(f1, 2) ) / 
                (pow(f1, 2) - pow(f2, 2))) * (L1/f1 - L2/f2)) / pow(10,16)
    elif units == 'm':
        sTEC = ((1/40.3) * (( pow(f2, 2) * pow(f1, 2) ) / 
                (pow(f1, 2) - pow(f2, 2))) * (L1 - L2)) / pow(10,16)
        
    return sTEC

def slantTEC(C1, C2, L1, L2, frequency = 2):
    global freq, c0
    F = (freq['1']**2 * freq[str(frequency)]**2) / (freq['1']**2 - freq[str(frequency)]**2)
    rangetec = F / 40.3 * (C2 - C1) / 1e16
    phasetec = F / 40.3 * c0 * (L1/freq['1'] - L2/freq[str(frequency)]) / 1e16
    N = np.nanmedian(rangetec - phasetec)
    
    return phasetec + N

def getPhaseCorrTEC(L1, L2, P1, P2, satbias=None, tec_err=False, channel=2,
                    intervals=None, fN = None, maxgap=3, maxjump=2):
    """
    Greg Starr
    Function returns a phase corrected TEC, following a paper by Coco el al:
    'Effect of GPS system biases on differential group delay measurements'.
    Imputs are raw data numpy arrays at two frequencies. Algorithm corrects
    the value of a phase TEC with a difference between mean values of a PRN and
    phase TEC, recpectively. This correction is performed on intervals between
    cycle slips. If the length of consequitive interval is 1, then NaN is inserted
    at this place.    
    """
    #
    global f1, f2, f5
    #Get intervals between nans and/or cycle slips    
    idx, ranges = getIntervals(L1, L2, P1, P2, maxgap=maxgap, maxjump=maxjump)
#    print ('Intervals between cycle slips ', ranges)
    ERR = np.nan * np.zeros(len(L1))
    TEC = np.nan * np.zeros(len(L1))
    for r in ranges:
        if (r[1] - r[0]) > 1:
            if fN is None: # GPS
                F1 = f1
                if channel == 2:
                    F2 = f2
                elif channel == 5:
                    F2 = f5
                range_tec = ((F1**2 * F2**2) / (F1**2 - F2**2)) * (P2[r[0] : r[1]] - 
                                          P1[r[0] : r[1]]) /40.3 / pow(10, 16)
                phase_tec = ((F1**2 * F2**2) / (F1**2 - F2**2)) * (c0/40.3) * \
                         (L1[r[0] : r[1]] / F1 - L2[r[0] : r[1]] / F2) / pow(10, 16)
            else: # GLONASS
                f1 = (1602 + fN*0.5625) * 1000000
                f2 = (1246 + fN*0.4375) * 1000000
                range_tec = ((f1**2 * f2**2) / (f1**2 - f2**2)) * (P2[r[0] : r[1]] - 
                                          P1[r[0] : r[1]]) /40.3 / pow(10, 16)
                phase_tec = ((f1**2 * f2**2) / (f1**2 - f2**2)) * (c0/40.3) * \
                         (L1[r[0] : r[1]] / f1 - L2[r[0] : r[1]] / f2) / pow(10, 16)
            tec_difference = np.array(sorted(phase_tec-range_tec))
            
            tec_difference = tec_difference[np.isfinite(tec_difference)]
            median_difference = tec_difference[int(len(tec_difference)/2)]
            difference_width = tec_difference[int(len(tec_difference)*.75)]-tec_difference[int(len(tec_difference)*.25)]
            median_error = difference_width/np.sqrt(len(tec_difference))
            tec = phase_tec - median_difference
            ERR[r[0]:r[1]] = median_error
            TEC[r[0]:r[1]] = tec
    if (tec_err):
        return TEC, ERR
    elif intervals:
        return TEC, ranges
    else:       
        return TEC 

def getPhaseCorrTECGLONASS(L1,L2,P1,P2, error=False):
    range_tec = ((f1**2 * f2**2) / (f1**2 - f2**2)) * (P2 - P1) / 40.3 / pow(10, 16)
    phase_tec = ((f1**2 * f2**2) / (f1**2 - f2**2)) * (c0/40.3) * (L1 / f1 - L2 / f2) / pow(10, 16)
    
    tec_difference = np.array(sorted(phase_tec-range_tec))
            
    tec_difference = tec_difference[np.isfinite(tec_difference)]
    median_difference = tec_difference[int(len(tec_difference)/2)]
    difference_width = tec_difference[int(len(tec_difference)*.75)]-tec_difference[int(len(tec_difference)*.25)]
    median_error = difference_width/np.sqrt(len(tec_difference))
    tec = phase_tec - median_difference
    
    if error:
        return tec, median_error
    else:
        return tec

def getVerticalTEC(tec, el, h, Fout=False):
    """
    Sebastijan Mrak
    Function takes the slant TEC numpy array, elevation as numpt array and an
    altitude 'h', to calculate the vertival TEC (vTEC). To map the sTEc to vTEC,
    we follow the thin shell approximation, described by Stefan Shaer in his
    desseration and is also illustrated in Brunini and Azpilicueta:
    'GPS slant total electron content accuracy using the single layer model 
    under different geomagnetic regions and ionospheric conditions'.
    Mapping function is defined as:
    
    F(el) = cos (arcsin (Re / (Re+h) * cos(el)))
    vTEC = sTEC * F(el)
    """
    Re = 6371.0
    h = h
    rc1 = (Re / (Re + h))
    vTEC =[]
    F = []
    for i in range(len(tec)):
        if np.isnan(tec[i]):
            vTEC.append(np.nan)
#            f = np.cos(np.arcsin(rc1*np.cos(np.radians(el[i]))))
            f = np.sqrt(1 - (np.cos(np.radians(el[i]))**2 * rc1**2))
            F.append(f)
        else:
#            f = np.cos(np.arcsin(rc1*np.cos(np.radians(el[i]))))
            f = np.sqrt(1 - (np.cos(np.radians(el[i]))**2 * rc1**2))
            vTEC.append(f * tec[i])
            F.append(f)
    
    if Fout:        
        return np.array(vTEC), np.array(F)
    else:
        return np.array(vTEC)
    
def singleFrequencyTEC(L1, C1, units='m', 
                       vertical:bool = False, 
                       el: np.ndarray = np.array([]), 
                       alt:int = 300):
    # L1 [cycle]: Convert C1 units to [cycles]
    if units == 'm':
        C1cycle = C1 * f1 / c0
    elif units == 'cycle':
        C1cycle = C1
    elif units == 'rad':
        C1cycle = C1 * f1 / c0 / (2*np.pi)
    else:
        raise ('Enter appropriate units "m", "rad", or "cycle"')
    dN = np.nanmean(L1 - C1cycle)
    L1corr = L1 - dN
    TECrad = C1cycle-L1corr
    if vertical:
        if el.shape[0] == L1.shape[0]:
            Fmap = getMappingFunction(el,alt)
            TECrad = TECrad * Fmap
    return TECrad

def retreiveDTECfromPhase(L, f=f1, units='cycle'):
    if units == 'cycle':
        dTEC = L * f**2 / 40.3 / 1E16
    elif units == 'rad':
        dTEC = L * c0 * f / 2 / np.pi / 40.3 /1e16
    else:
        raise ('Enter an appropriate name for units')
    return dTEC * -1
#%% Indexes
def getROTI(tec, length):
    """
    Sebastijan Mrak
    getROTI returns the rate of TEC Index calculated as the standard deviation 
    of the provided TEC on the moving window of the length 'length'. It returns 
    the ROTI as a numpy array data type.
    """
    roti = []    
    for i in range(len(tec)-length):
        roti.append(np.std(tec[i:i+length]))
    
    return np.array(roti)
    

def phaseScintillationIndex(data, N):
    """
    Sebastijan Mrak
    GNSS Phase scintillation index for the interval of the length 'N' samples
    """
    y = np.nan * np.zeros(data.shape[0]-N)
    for i in range(data.shape[0] - N):
        y[i] = np.std(data[i:i+N])
    return y
    
def AmplitudeScintillationIndex(data, N):
    """
    Sebastijan Mrak
    GNSS Amplitude scintillation index for the interval of the length 'N' samples
    """
    y = np.nan * np.zeros(data.shape[0])
    for i in range(data.shape[0] - N):
        y[i] = np.std(data[i:i+N] / np.mean(data[i:i+N]))
    return y

def gpsSatPositionSP3(fsp3, dt, sv=None, rx_position=None, coords='xyz'):
    assert sv is not None
    # Read in data
    D = gr.load(fsp3).sel(sv=sv)
    if isinstance(dt, list):
        dt = np.asarray(dt)
    dt = dt.astype('datetime64[s]')
    navtimes = D.time.values.astype('datetime64[s]')
    CSx = interpolate.CubicSpline(navtimes.astype(int), D.position.values[:,0])
    CSy = interpolate.CubicSpline(navtimes.astype(int), D.position.values[:,1])
    CSz = interpolate.CubicSpline(navtimes.astype(int), D.position.values[:,2])
    ecefxi = CSx(dt.astype(int))*1e3
    ecefyi = CSy(dt.astype(int))*1e3
    ecefzi = CSz(dt.astype(int))*1e3
    
    if coords == 'xyz':
        return np.array([ecefxi, ecefyi, ecefzi])
    else:
        AER = ecef2aer(x=ecefxi, y=ecefyi, z=ecefzi,
                       lon0=rx_position[1], lat0=rx_position[0], h0=rx_position[2])
        return np.array(AER)

def gpsSatPosition(fnav, dt, sv=None, rx_position=None, coords='xyz'):
    navdata = gr.load(fnav).sel(sv=sv)
    timesarray = np.asarray(dt,dtype='datetime64[ns]') #[datetime64 [ns]]
    # Manipulate with times, epochs and crap like this
    navtimes = navdata.time.values # [datetime64 [ns]]
    idnan = np.isfinite(navdata['Toe'].values)
    navtimes = navtimes[idnan]
    bestephind = []
    for t in timesarray:
        idt = abs(navtimes - t).argmin() #if t>navtimes[abs(navtimes - t).argmin()] else abs(navtimes - t).argmin()-1
        bestephind.append(idt)
#    bestephind = np.array([np.argmin(abs(navtimes-t)) for t in timesarray])
    gpstime = np.array([getGpsTime(t) for t in dt])
    t = gpstime - navdata['Toe'][idnan][bestephind].values # [datetime.datetime]
    # constants
    GM = 3986005.0E8 # universal gravational constant
    OeDOT = 7.2921151467E-5
    # Elements
    ecc = navdata['Eccentricity'][idnan][bestephind].values # Eccentricity
    Mk = navdata['M0'][idnan][bestephind].values + \
         t *(np.sqrt(GM / navdata['sqrtA'][idnan][bestephind].values**6) + 
             navdata['DeltaN'][idnan][bestephind].values)
    Ek = solveIter(Mk,ecc)
    Vk = np.arctan2(np.sqrt(1.0 - ecc**2) * np.sin(Ek), np.cos(Ek) - ecc)
    PhiK = Vk + navdata['omega'][idnan][bestephind].values
    # Perturbations
    delta_uk = navdata['Cuc'][idnan][bestephind].values * np.cos(2.0*PhiK) + \
              navdata['Cus'][idnan][bestephind].values * np.sin(2.0*PhiK)
    Uk = PhiK + delta_uk
    
    delta_rk = navdata['Crc'][idnan][bestephind].values * np.cos(2.0*PhiK) + \
               navdata['Crs'][idnan][bestephind].values * np.sin(2.0*PhiK)
    Rk = navdata['sqrtA'][idnan][bestephind].values**2 * (1.0 - ecc * np.cos(Ek)) + delta_rk
    
    delta_ik = navdata['Cic'][idnan][bestephind].values * np.cos(2.0*PhiK) + \
               navdata['Cis'][idnan][bestephind].values * np.sin(2.0*PhiK)
    Ik = navdata['Io'][idnan][bestephind].values + \
         navdata['IDOT'][idnan][bestephind].values * t + delta_ik
    
    #Compute the right ascension
    Omegak = navdata['Omega0'][idnan][bestephind].values + \
             (navdata['OmegaDot'][idnan][bestephind].values - OeDOT) * t - \
             (OeDOT * navdata['Toe'][idnan][bestephind].values)
             
    #X,Y coordinate corrections
    Xkprime = Rk * np.cos(Uk)
    Ykprime = Rk * np.sin(Uk)
    # ECEF XYZ
    X = Xkprime * np.cos(Omegak) - (Ykprime * np.sin(Omegak) * np.cos(Ik))
    Y = Xkprime * np.sin(Omegak) + (Ykprime * np.cos(Omegak) * np.cos(Ik))
    Z = Ykprime * np.sin(Ik)
    
    if coords == 'xyz':
        return np.array([X,Y,Z])
    elif coords == 'aer':
        assert rx_position is not None
        rec_lat, rec_lon, rec_alt = ecef2geodetic(rx_position[0], rx_position[1], rx_position[2])
        A,E,R = ecef2aer(X, Y, Z, rec_lat, rec_lon, rec_alt)
        return np.array([A,E,R])

#%% Sat positioning 
    
def aer2ipp(aer, rxp, H=350):
    H*=1e3
    aer_new = np.copy(aer)
    fm = np.sin(np.radians(aer[:,:,1]))
    aer_new[:,:,2] = (H / fm)
    rxp
    ipp = np.array(aer2geodetic(aer_new[:,:,0], aer_new[:,:,1], aer_new[:,:,2], 
                                rxp[0], rxp[1], rxp[2]))
    return ipp

def getIonosphericPiercingPoints(rx_xyz, sv, obstimes, ipp_alt, navfn,
                                 cs='wsg84', rx_xyz_coords='xyz', el0=0):
    """
    Sebastijan Mrak
    Function returns a list of Ionospheric Piersing Point (IPP) trajectory in WSG84
    coordinate system (CS). Function takes as parameter a receiver location in 
    ECEF CS, satellite number, times ob observation and desired altitude of a IPP
    trajectory. You also have to specify a full path to th navigation data file.
    It returns IPP location in either WSG84 or AER coordinate system.
    """
    
    ipp_alt = ipp_alt * 1E3
    if rx_xyz_coords == 'xyz':
        rec_lat, rec_lon, rec_alt = ecef2geodetic(rx_xyz[0], rx_xyz[1], rx_xyz[2])
    else:
        if not isinstance(rx_xyz, list): rx_xyz = list(rx_xyz)
        if len(rx_xyz) == 2: rx_xyz.append(0)
        assert len(rx_xyz) == 3
        
        rec_lat = rx_xyz[0]
        rec_lon = rx_xyz[1]
        rec_alt = rx_xyz[2]
        rx_xyz = geodetic2ecef(lat = rx_xyz[0], lon = rx_xyz[1], alt = rec_alt)
    
    if sv[0] == 'G':
        if navfn.endswith('n'):
            xyz = gpsSatPosition(navfn, obstimes, sv=sv, rx_position=rx_xyz, coords='xyz')
        elif navfn.endswith('sp3'):
            xyz = gpsSatPositionSP3(navfn, obstimes, sv=sv, rx_position=rx_xyz, coords='xyz')
        az,el,r = ecef2aer(xyz[0,:],xyz[1,:],xyz[2,:],rec_lat, rec_lon, rec_alt)
        aer_vector = np.array([az, el, r])
        
        r_new = []
        for i in range(len(el)):
            if el[i] > el0:
                fm = np.sin(np.radians(el[i]))
                r_new.append(ipp_alt / fm)
            else:
                r_new.append(np.nan)
        r_new = np.array(r_new)
        if cs == 'wsg84':
            idx = (el > el0)
            
            if el.size > 1:
                lla_vector = np.array(aer2geodetic(az[idx], el[idx], r_new[idx], rec_lat, rec_lon, rec_alt))
            else:
                if el > el0:
                    lla_vector = np.array(aer2geodetic(az, el, r_new, rec_lat, rec_lon, rec_alt))
                else:
                    lla_vector = np.array([np.nan, np.nan, np.nan])
        
    elif sv[0] == 'R':
        aer_vector = gloSatPosition(navfn=navfn, sv=sv, obstimes=obstimes, rx_position=[rec_lon, rec_lat, rec_alt], cs='aer')
        fm = np.sin(np.radians(aer_vector[1]))
        r_new = ipp_alt / fm
        lla_vector = np.array(aer2geodetic(aer_vector[0], aer_vector[1], r_new, rec_lat, rec_lon, rec_alt))
    else:
        print ('Type in valid sattype initial. "G" for GPS and "R" for GLONASS')
        
    if (cs == 'wsg84'):
        return lla_vector
    elif (cs == 'aer'):
        return aer_vector
    else:
        print ('Enter either "wsg84" or "aer" as coordinate system. "wsg84" is default one.')
# %%
def getMappingFunction(el, h):
    
    Re = 6371.0
    rc1 = (Re / (Re + h))
    arg = rc1 * np.cos(np.radians(el))
    F = np.cos(np.arcsin(arg))
#    F = np.sqrt(1 - (np.cos(np.radians(el))**2 * rc1**2))
    return np.array(F)

def gloSatPosition(navfn, sv, obstimes, rx_position=None, cs='xyz'):
    obstimes = obstimes.astype('datetime64[s]')
    D = gr.load(navfn).sel(sv = sv)
    navtimes = D.time.values.astype('datetime64[s]')
    x = D.X.values * 1e3
    y = D.Y.values * 1e3
    z = D.Z.values * 1e3
    # Rearange
#    deltaTin = np.diff(obstimes)[0]
    tmin = min(obstimes.min(), navtimes.min())
    tmax = max(obstimes.max(), navtimes.max())
    navtimes_interp = np.arange(tmin, tmax + 1, 
                                timedelta(seconds = 1)).astype('datetime64[s]') #deltaTin.item().total_seconds())).astype('datetime64[s]')
    x0 = np.linspace(0, 1, navtimes.shape[0])
    x1 = np.linspace(0, 1, navtimes_interp.shape[0])
    Xi = np.interp(x1, x0, x)
    Yi = np.interp(x1, x0, y)
    Zi = np.interp(x1, x0, z)
#    if isinstance(obstimes, (np.ndarray)):
    idt = np.isin(navtimes_interp, obstimes)#(navtimes_interp >= obstimes[0]) & (navtimes_interp <= obstimes[-1])
    
    if cs == 'xyz':
        xyz = np.array([Xi[idt], Yi[idt], Zi[idt]])
        return xyz
    elif cs == 'aer':
        assert rx_position is not None
        if isinstance(rx_position, list): rx_position = np.array(rx_position)
        assert rx_position.shape[0] == 3
        aer = ecef2aer(Xi[idt],Yi[idt],Zi[idt],
                       lon0 = rx_position[0], lat0 = rx_position[1],
                       h0 = rx_position[2])
        return aer
    else:
        return
    
# %%
#def getSatXYZ(nav, sv, times):
#    """
#    Greg Starr
#    getSatelliteXYZ returns the satellite XYZ as a tuple at the inputted times
#    inputs are rinex navigation data, satellite number, and list of times
#    Output: tuple of satellite position in ECEF coordinates (X,Y,Z)
#    Algorithm: Based on http://web.ics.purdue.edu/~ecalais/teaching/geodesy/EAS_591T_2003_lab_4.htm
#    also based on Bill Rideout's tec.py
#    """
#    allSvInfo = nav[nav['sv']==sv] 
#    
#    timesarray = np.asarray(times,dtype='datetime64[ms]')
#    navtimes = np.asarray(allSvInfo.index,dtype='datetime64[ms]')
#    bestephind = np.array([np.argmin(abs(navtimes-t)) for t in timesarray])
#    info = np.asarray(allSvInfo)[bestephind]
#    info = DataFrame(info,index=times,columns=allSvInfo.columns)
#    info['sv'] = sv
#    info['gpstime'] = np.array([getGpsTime(t) for t in times])
#    
#    # constants
#    GM = 3986005.0E8 # universal gravational constant
#    OeDOT = 7.2921151467E-5
#    
#    #Basic Parameters
#    t = info['gpstime']-info['TimeEph']
#    mu = info['M0']+t*(np.sqrt(GM/info['sqrtA']**6)+info['DeltaN'])
#    Ek = solveIter(mu,info['Eccentricity'])  
#    Vk = np.asarray(np.arctan2(np.sqrt(1.0-info['Eccentricity']**2)*np.sin(Ek),
#                               np.cos(Ek)-info['Eccentricity']),float)
#    PhiK = Vk + info['omega']
#    #Correct for orbital perturbations
#    omega = np.asarray(info['omega']+info['Cus']*np.sin(2.0*PhiK)
#             +info['Cuc']*np.cos(2.0*PhiK),float)
#    r = np.asarray((info['sqrtA']**2)*(1.0-info['Eccentricity']*np.cos(Ek))
#         +info['Crs']*np.sin(2.0*PhiK)+info['Crc']*np.cos(2.0*PhiK),float)
#    i = np.asarray(info['Io']+info['IDOT']*t+info['CIS']*np.sin(2.0*PhiK)
#         +info['Cic']*np.cos(2.0*PhiK),float)
#    
#    #Compute the right ascension
#    Omega = np.asarray(info['OMEGA']+(info['OMEGA DOT']-OeDOT)*t-
#        (OeDOT*info['TimeEph']),float)
#    #Convert satellite position from orbital frame to ECEF frame
#    cosOmega = np.cos(Omega)
#    sinOmega = np.sin(Omega)
#    cosomega = np.cos(omega)
#    sinomega = np.sin(omega)
#    cosi = np.cos(i)
#    sini = np.sin(i)
#    cosVk = np.cos(Vk)
#    sinVk = np.sin(Vk)
#    R11 = cosOmega*cosomega - sinOmega*sinomega*cosi
#    R12 = -1.0*cosOmega*sinomega - sinOmega*cosomega*cosi
#    #R13 = np.sin(Omega)*np.sin(i)
#    R21 = sinOmega*cosomega + cosOmega*sinomega*cosi
#    R22 = -1.0*sinOmega*sinomega + cosOmega*cosomega*cosi
#    #R23 = -1.0*np.cos(Omega)*np.sin(i)
#    R31 = sinomega*sini
#    R32 = cosomega*sini
#    #R33 = np.cos(i)
#          
#    xyz = np.zeros((len(times),3))
#    rv = np.column_stack((r*cosVk,r*sinVk,np.zeros(r.shape)))
#    
#    R = np.empty((rv.shape[0],3,3))
#    R[:,0,0] = R11
#    R[:,0,1] = R12
#    R[:,0,2] = 0
#    R[:,1,0] = R21
#    R[:,1,1] = R22
#    R[:,1,2] = 0
#    R[:,2,0] = R31
#    R[:,2,1] = R32
#    R[:,2,2] = 0
#    
#
#    for i in range(len(times)): #THIS IS THE SLOWEST PART NOW
#        xyz[i,:] = (R[i,:,:].dot(rv[i,:]))
#        
#    return xyz
# %% GeoRinex NAV file access and processing
#def getSatXYZ2(info, times):
#    """
#    Greg Starr
#    getSatelliteXYZ returns the satellite XYZ as a tuple at the inputted times
#    inputs are rinex navigation data, satellite number, and list of times
#    Output: tuple of satellite position in ECEF coordinates (X,Y,Z)
#    Algorithm: Based on http://web.ics.purdue.edu/~ecalais/teaching/geodesy/EAS_591T_2003_lab_4.htm
#    also based on Bill Rideout's tec.py
#    
#    Update:
#    nav == xarray from georinex = gr.load(brdc.**n)
#    """    
#    # Input argument: 'times' [datetime.datetime]
#    timesarray = np.asarray(times,dtype='datetime64[ns]') #[datetime64 [ns]]
#    # Manipulate with times, epochs and crap like this
#    navtimes = info.time.values # [datetime64 [ns]]
#    idnan = np.isfinite(info['Toe'].values)
#    bestephind = np.array([np.argmin(abs(navtimes[idnan]-t)) for t in timesarray])
#    gpstime = np.array([getGpsTime(t) for t in times])
#    t = gpstime - info['Toe'][idnan][bestephind].values # [datetime.datetime]
#    # constants
#    GM = 3986005.0E8 # universal gravational constant
#    OeDOT = 7.2921151467E-5
#    ##################3 REPLACEMENT ###########################################
#    # Basic parameters
#    ecc = info['Eccentricity'][idnan][bestephind].values # Eccentricity
#    mu = info['M0'][idnan][bestephind].values + \
#             t *(np.sqrt(GM / info['sqrtA'][idnan][bestephind].values**6) + 
#             info['DeltaN'][idnan][bestephind].values)
#    Ek = solveIter(mu,ecc)
#    Vk = np.asarray(np.arctan2(np.sqrt(1.0 - ecc**2) * np.sin(Ek),
#                               np.cos(Ek) - ecc), float)
#    PhiK = Vk + info['omega'][idnan][bestephind].values
#    
#    #Correct for orbital perturbations
#    omega = np.asarray(info['omega'][idnan][bestephind].values + 
#                       info['Cus'][idnan][bestephind].values * np.sin(2.0*PhiK) + 
#                       info['Cuc'][idnan][bestephind].values * np.cos(2.0*PhiK), float)
#    r = np.asarray((info['sqrtA'][idnan][bestephind].values**2) * 
#                   (1.0 - ecc * np.cos(Ek)) + info['Crs'][idnan][bestephind].values * 
#                   np.sin(2.0*PhiK) + info['Crc'][idnan][bestephind].values * 
#                   np.cos(2.0*PhiK), float)
#    i = np.asarray(info['Io'][idnan][bestephind].values + 
#                   info['IDOT'][idnan][bestephind].values * t + 
#                   info['Cis'][idnan][bestephind].values * np.sin(2.0*PhiK) + 
#                   info['Cic'][idnan][bestephind].values * np.cos(2.0*PhiK), float)
#    
#    #Compute the right ascension
#    Omega = np.asarray(info['Omega0'][idnan][bestephind].values + 
#                       (info['OmegaDot'][idnan][bestephind].values - OeDOT) * t -
#                       (OeDOT * info['Toe'][idnan][bestephind].values), float)
#    ###########################################################################
#    #Convert satellite position from orbital frame to ECEF frame
#    cosOmega = np.cos(Omega)
#    sinOmega = np.sin(Omega)
#    cosomega = np.cos(omega)
#    sinomega = np.sin(omega)
#    cosi = np.cos(i)
#    sini = np.sin(i)
#    cosVk = np.cos(Vk)
#    sinVk = np.sin(Vk)
#    R11 = cosOmega*cosomega - sinOmega*sinomega*cosi
#    R12 = -1.0*cosOmega*sinomega - sinOmega*cosomega*cosi
#    #R13 = np.sin(Omega)*np.sin(i)
#    R21 = sinOmega*cosomega + cosOmega*sinomega*cosi
#    R22 = -1.0*sinOmega*sinomega + cosOmega*cosomega*cosi
#    #R23 = -1.0*np.cos(Omega)*np.sin(i)
#    R31 = sinomega*sini
#    R32 = cosomega*sini
#    #R33 = np.cos(i)
#          
#    xyz = np.zeros((len(times),3))
#    rv = np.column_stack((r*cosVk,r*sinVk,np.zeros(r.shape)))
#    
#    R = np.empty((rv.shape[0],3,3))
#    R[:,0,0] = R11
#    R[:,0,1] = R12
#    R[:,0,2] = 0
#    R[:,1,0] = R21
#    R[:,1,1] = R22
#    R[:,1,2] = 0
#    R[:,2,0] = R31
#    R[:,2,1] = R32
#    R[:,2,2] = 0
#    
#
#    for i in range(len(times)): #THIS IS THE SLOWEST PART NOW
#        xyz[i,:] = (R[i,:,:].dot(rv[i,:]))
#        
#    return xyz
#%% Phase scintillation
#def phaseDetrend(y, order,polynom=False):
#    """
#    Sebastijan Mrak
#    Raw phase data detrending using N-th polinom approximation function.
#    Detrended output is input data subtracted with polinom approximation.
#    Output is of the same length as input data 'y'. 
#    """
#    x = np.arange(y.shape[0])
#    mask = np.isnan(y)
#    z = np.polyfit(x[~mask], y[~mask], order)
#    z = np.polyfit(x, y, order)
#    f = np.poly1d(z)
#    polyfit = f(x)
#    y_d = y-polyfit
#    
#    if polynom:
#        return y_d, polyfit
#    else:
#        return y_d


def phaseScintillation(data, fc=0.1, filt_order=6, polyfit_order=3, fs=1, skip=20):
    """
    Sebastijan Mrak
    Standard pocedure to obtain detrended carrier phase scintillation. Input 
    is a list of carrier phase measurements, than this list goes to a process of
    polynominal detrending and finaly detrended data is high-pass filterd (IIR).
    """
    L = np.nan*np.zeros(data.shape[0])
    idx = np.where(np.isfinite(data))[0]
    L1phi = data[idx]
    L1_d = uf.phaseDetrend(L1phi, polyfit_order)
    Y = uf.hpf(L1_d, fc=fc, order=filt_order, fs=fs)
    
    L[idx[skip:]] = Y[skip:]

    return L

# %% Utils
def getIntervals(L1, L2, P1, P2, maxgap=3,maxjump=2):
    """
    Greg Starr
    scans through the phase tec of a satellite and determines where "good"
    intervals begin and end
    inputs:
        L1, L2, P1, P2 - np arrays of the raw data used to calculate phase corrected
        TEC. 
        maxgap - maximum number of nans before starting new interval
        maxjump - maximum jump in phase TEC before starting new interval
    output:
        intervals - list of 2-tuples, beginning and end of each "good" interval
                    as a Pandas/numpy
    """
    f1 = 1575420000
    f2 = 1227600000
    
    r = np.array(range(len(P1)))
    idx = np.isfinite(L1) & np.isfinite(L2) & np.isfinite(P1) & np.isfinite(P2)
    r = r[idx]
    intervals=[]
    if len(r)==0:
        return idx, intervals
    phase_tec=2.85E9*(L1/f1-L2/f2)
    beginning=r[0]
    last=r[0]
    for i in r[1:]:
        if (i-last>maxgap) or (abs(phase_tec[i]-phase_tec[last])>maxjump):
            intervals.append((beginning,last))
            beginning=i
        last=i
        if i==r[-1]:
            intervals.append((beginning,last))
    return idx, intervals

def getIntervalsTEC(y, maxgap=3, maxjump=2):

    r = np.arange(y.size)
    idx = np.isfinite(y)
    r = r[idx]
    intervals=[]
    if len(r)==0:
        return idx, intervals

    beginning=r[0]
    last=r[0]
    for i in r[1:]:
        if (i-last > maxgap) or (abs(y[i] - y[last]) > maxjump):
            intervals.append((beginning, last))
            beginning=i
        last=i
        if i==r[-1]:
            intervals.append((beginning, last))
    return idx, intervals

def solveIter(mu,e):
    """__solvIter returns an iterative solution for Ek
    Mk = Ek - e sin(Ek)
    adapted to accept vectors instead of single values
    from Bill Rideout's tec.py
    """
    thisStart = np.asarray(mu-1.01*e)
    thisEnd = np.asarray(mu + 1.01*e)
    bestGuess = np.zeros(mu.shape)

    for i in range(5): 
        minErr = 10000*np.ones(mu.shape)
        for j in range(5):
            thisGuess = thisStart + j*(thisEnd-thisStart)/10.0
            thisErr = np.asarray(abs(mu - thisGuess + e*np.sin(thisGuess)))
            mask = thisErr<minErr
            minErr[mask] = thisErr[mask]
            bestGuess[mask] = thisGuess[mask]
        
        # reset for next loop
        thisRange = thisEnd - thisStart
        thisStart = bestGuess - thisRange/10.0
        thisEnd = bestGuess + thisRange/10.0
        
    return(bestGuess)


def getGpsTime(dt):
    """_getGpsTime returns gps time (seconds since midnight Sat/Sun) for a datetime
    """
    total = 0
    days = (dt.weekday()+ 1) % 7 # this makes Sunday = 0, Monday = 1, etc.
    total += days*3600*24
    total += dt.hour * 3600
    total += dt.minute * 60
    total += dt.second
    return(total)

def tecdPerLOS(stec, intervals, mask=None, eps=1, 
               tsps=30, polynom_list=None, zero_mean=False):
    tecd = np.nan * np.ones(stec.size)
    if mask is None:
        mask = np.zeros(stec.size, dtype=bool)
    for ir, r in enumerate(intervals):
        chunk = stec[r[0]+1 : r[1]-1]
        idf = np.isfinite(chunk)
        if np.sum(np.isfinite(chunk)) < (15 * (60/tsps)): 
            err_list = np.array([])
            continue
        if np.sum(np.isnan(chunk)) > 0:
            chunk = uf.cubicSplineFit(chunk, idf)
        
        res, err_list0, po  = uf.detrend(chunk, polynom_list=polynom_list, eps=eps, mask=mask[r[0]+1 : r[1]-1], polynomial_order=True)
        if ir == 0 or len(err_list) == 0:
            err_list = err_list0
        else:
            err_list = np.vstack((err_list, err_list0))
        res[~idf] = np.nan
        if zero_mean:
            if abs(np.nansum(res)) < 5:
                tecd[r[0]+1 : r[1]-1] = res
        else:
            tecd[r[0]+1 : r[1]-1] = res
    
    return tecd

# %%
def singleRx(obs, nav, sv='G23', args=['L1','S1'], tlim=None,rawplot=False,
             sp=True,s4=False,polyfit=False,indicator=False,
             alt=300,el_mask=30,skip=20,porder=8,tec_ch=2,
             forder=5,fc=0.1,fs=1):
    
    def _plot(x,y,title=''):
        plt.figure()
        plt.title(title)
        plt.plot(x,y,'b')
    
    D = xarray.open_dataset(obs, group='OBS')
    rx_xyz = D.position
    leap_seconds = uf.getLeapSeconds(nav)
    obstimes64 = D.time.values
    times = np.array([Timestamp(t).to_pydatetime() for t in obstimes64]) - \
                        timedelta(seconds = leap_seconds)
    # define tlim ie. skip
    if tlim is not None:
        s = ((times>=tlim[0]) & (times<=tlim[1]))
    else:
        s = np.full(times.shape[0], True, dtype=bool)
        s[:skip] = False
    times = times[s]
    # Get satellite position
    aer = gpsSatPosition(nav, times, sv=sv, rx_position=rx_xyz, coords='aer')
    # Elevation mask
    idel = aer[1] >= el_mask
    times = times[idel]
    # times to output dict: Y
    Y = {'times': times}
    Y['rx_xyz'] = rx_xyz
    Y['az'] = aer[0][idel]
    Y['el'] = aer[1][idel]
    for arg in args:
        if not arg[1:] == 'TEC':
            X = D.sel(sv=sv)[arg].values[s][idel]
            # To dict
            Y[arg] = X
            # Indicators?
            if indicator:
                try:
                    if arg[0] == 'L' and arg[-1] != '5':
                        argi = arg + 'lli'
                    elif arg[0] == 'C' or arg[0] == 'P' or arg == 'L5':
                        argi = arg + 'ssi'
                    else:
                        argi = ''
                    Xi = D.sel(sv=sv)[argi].values[s:][idel]
                    Y[argi] = Xi
                except Exception as e:
                    print (e)
            if rawplot:
                _plot(times,X,arg)
        if arg[0] == 'L' or arg[0] == 'C':
            if arg == 'C1':
                X = X * f1 / c0 # To cycles
            Xd = uf.phaseDetrend(X, order=porder)
                
            if polyfit:
                # To dict
                Y[arg+'polyfit'] = Xd
            if rawplot:
                if rawplot:
                    _plot(times,Xd,arg+'_detrend_poly')
            if sp:
                Xy = uf.hpf(Xd, order=forder,fc=fc,plot=False,fs=fs)
                # To dict
                Y['sp'] = Xy
                if rawplot:
                    _plot(times[100:],Xy[100:],arg+'_scint')
        if arg[0] == 'S':
            if s4:
                Xy = AmplitudeScintillationIndex(X,60)
                # To dict
                Y['s4'] = Xy
                if rawplot:
                    _plot(times,Xy,arg+'_S4')
                    
        if arg[1:] == 'TEC':
            if tec_ch == 2:
                C1 = D.sel(sv=sv)['C1'].values[s][idel]
                L1 = D.sel(sv=sv)['L1'].values[s][idel]
                C2 = D.sel(sv=sv)['P2'].values[s][idel]
                L2 = D.sel(sv=sv)['L2'].values[s][idel]
            elif tec_ch == 5:
                C1 = D.sel(sv=sv)['C1'].values[s][idel]
                L1 = D.sel(sv=sv)['L1'].values[s][idel]
                C2 = D.sel(sv=sv)['C5'].values[s][idel]
                L2 = D.sel(sv=sv)['L5'].values[s][idel]
            sTEC = getPhaseCorrTEC(L1, L2, C1, C2,channel=tec_ch)
            sTEC = sTEC - np.nanmin(sTEC)
            if arg == 'sTEC':
                # To dict
                Y[arg] = sTEC
                if rawplot:
                    _plot(times, sTEC,title=arg)
            elif arg == 'vTEC':
                vTEC = getVerticalTEC(sTEC,aer[1][idel],alt)
                # To dict
                Y[arg] = vTEC
                if rawplot:
                    _plot(times, vTEC, title=arg)
    # return dict with data/args
    return Y

def dataFromNC(fnc,fnav, sv, fsp3=None,
               el_mask=30,tlim=None, 
               satpos=False,
               ipp=False,
               ipp_alt=None):
    leap_seconds = uf.getLeapSeconds(fnav)
    D = gr.load(fnc, useindicators=True).sel(sv=sv)
    D['time'] = np.array([np.datetime64(ttt) for ttt in D.time.values])
    if tlim is not None:
        if len(tlim) == 2:
            t0 = tlim[0] + timedelta(seconds=leap_seconds)
            t1 = tlim[1] + timedelta(seconds=leap_seconds)
            D = D.where(np.logical_and(D.time >= np.datetime64(t0), 
                                       D.time <= np.datetime64(t1)), 
                                       drop=True)
    obstimes64 = D.time.values
    dt = np.array([Timestamp(t).to_pydatetime() for t in obstimes64]) - \
               timedelta(seconds = leap_seconds)
    if fsp3 is None:
        aer = gpsSatPosition(fnav, dt, sv=sv, rx_position=D.position, coords='aer')
    else:
        aer = gpsSatPositionSP3(fsp3, dt, sv=sv, rx_position=D.position_geodetic, coords='aer')
    idel = (aer[1,:] >= el_mask)
    aer[:,~idel] = np.nan
    D['time'] = dt

    if satpos:
        D['az'] = aer[0,:]
        D['el'] = aer[1,:]
        D['idel'] = idel
    if ipp:
        if ipp_alt is None:
            print ('Auto assigned altitude of the IPP: 250 km')
            ipp_alt = 250e3
        else:
            ipp_alt *= 1e3
        rec_lat, rec_lon, rec_alt = D.position_geodetic
        fm = np.sin(np.radians(aer[1]))
        r_new = ipp_alt / fm
        lla_vector = np.array(aer2geodetic(aer[0], aer[1], r_new, rec_lat, rec_lon, rec_alt))
        D['ipp_lon'] = lla_vector[1]
        D['ipp_lat'] = lla_vector[0]
        D['ipp_alt'] = ipp_alt
    return D

def processTEC(obs, sv, Ts = 30, frequency = 2, H=None, elevation=None, 
               rx_bias=None, vtec = False, F_out=False, sat_bias=None):
    if isinstance(obs, dict):
        tec = slantTEC(obs['C1'], obs['P2'], 
                        obs['L1'], obs['L2'], 
                        frequency = frequency)
    elif isinstance(obs, xarray.Dataset):
        tec = getPhaseCorrTEC(obs['L1'].values, obs['L2'].values, 
                        obs['C1'].values, obs['P2'].values, 
                        )
        
    if sat_bias is not None:
        try:
            b = 1 * getSatBias(sat_bias, sv=sv)
        except:
            b = 0
    else:
        b = 0
    if rx_bias is not None:
        b += rx_bias
    assert elevation is not None
    assert H is not None
    if vtec:
        tec += b
        F = getMappingFunction(elevation, h = H)
        tec = tec * F
    
    dtec = uf.getPlainResidual(tec, Ts=Ts)
    
    if F_out:
        return tec, dtec, F
    else:
        return tec, dtec
    
def getSTEC(fnc, fsp3 = None, el_mask=30, H=350, maxgap=1, maxjump=1.6):
    D = gr.load(fnc)
    stec = np.nan * np.zeros((D.time.values.size, D.sv.size))
    for isv, sv in enumerate(D.sv.values):
        
        if fsp3 is not None:
            dt = np.array([np.datetime64(ttt) for ttt in D.time.values])
            assert os.path.exists(fsp3)
            aer = getIonosphericPiercingPoints(D.position, sv, dt, 
                                                   ipp_alt=H, navfn=fsp3,
                                                   cs='aer', rx_xyz_coords='xyz')
            idel = (aer[1] >= el_mask)
        else:
            idel = np.ones(D.time.values.size, dtype=bool)
        stec[idel, isv] = getPhaseCorrTEC(L1=D.L1.values[idel,isv], L2=D.L2.values[idel,isv],
                                         P1=D.C1.values[idel,isv], P2=D.P2.values[idel,isv],
                                         maxgap=maxgap, maxjump=maxjump)
    return stec

def getVTEC(fnc, fsp3, dcb=None, jplg_file=None, el_mask=30, H=350,
            tskip=None, maxgap=1, maxjump=1.6, tec_shift=None,
            return_mapping_function=False, return_aer=False):
    if tec_shift is None: tec_shift = 0
    if dcb is None:
#        assert jplg_file is not None
#        assert os.path.exists(jplg_file)
        dcb = getDCB(fnc, fsp3, jplg_file=jplg_file, el_mask=30, H=350,
                 tskip=None, maxgap=1, maxjump=1.6)
    D = gr.load(fnc)
    assert dcb.size == D.sv.size
    dt = np.array([np.datetime64(ttt) for ttt in D.time.values])
    vtec = np.nan * np.zeros((dt.size, D.sv.size))
    if return_mapping_function:
        F = np.nan * np.zeros((dt.size, D.sv.size))
    if return_aer:
        AER = np.nan * np.zeros((dt.size, D.sv.size, 3))
    for isv, sv in enumerate(D.sv.values):
        aer = getIonosphericPiercingPoints(D.position, sv, dt, 
                                           ipp_alt=H, navfn=fsp3,
                                           cs='aer', rx_xyz_coords='xyz')
        idel = (aer[1] >= el_mask)
        if return_mapping_function:
            F[idel, isv] = getMappingFunction(aer[1][idel], h=H)
        vtec[idel, isv] = (getPhaseCorrTEC(L1=D.L1.values[idel,isv], L2=D.L2.values[idel,isv],
                                     P1=D.C1.values[idel,isv], P2=D.P2.values[idel,isv],
                                     maxgap=maxgap, maxjump=maxjump) - dcb[isv] + tec_shift) * getMappingFunction(aer[1][idel], h=H)
        if return_aer:
            AER[idel, isv, 0] = aer[0][idel]
            AER[idel, isv, 1] = aer[1][idel]
            AER[idel, isv, 2] = aer[2][idel]
    D.close()
    if return_mapping_function and (not return_aer):
        return vtec, F
    elif return_mapping_function and return_aer:
        return vtec, F, AER
    elif (not return_mapping_function) and return_aer:
        return vtec, AER
    else:
        return vtec

def getVTEC2(D, F, tskip=1, el_mask=30, maxgap=1, maxjump=1):
    
    def _fun(p, stec, F):
        vtec = (stec - p) * F
        ret = np.nansum(np.nanstd(vtec, axis=1)**2)
        return ret
    
    stec = np.nan * np.copy(F)
    x0 = np.ones(D.sv.size)
    for isv, sv in enumerate(D.sv.values):
        idel = (F[::tskip, isv] >= getMappingFunction(el_mask, 350))
        stec[idel, isv] = getPhaseCorrTEC(L1=D.L1.values[::tskip][idel,isv], L2=D.L2.values[::tskip][idel,isv],
                                     P1=D.C1.values[::tskip][idel,isv], P2=D.P2.values[::tskip][idel,isv],
                                     maxgap=maxgap, maxjump=maxjump)
        
    sb_lsq = least_squares(_fun, x0, args=(stec, F), xtol=1e-5, gtol=1e-5, 
                       diff_step=1e-4, loss='cauchy')
    dcb = sb_lsq.x
    D.close()
    
    vtec = (stec - dcb[isv]) * F
    
    return vtec
        
def getDCB(fnc, fsp3, jplg_file=None, el_mask=30, H=350, 
            tskip=None, maxgap=1, maxjump=1.6,
            return_mapping_function=False, return_aer=False):
    
    def _fun(p, stec, F):
        vtec = (stec - p) * F
        ret = np.nansum(np.nanstd(vtec, axis=1)**2)
        return ret
    D = gr.load(fnc)
    dt0 = np.array([np.datetime64(ttt) for ttt in D.time.values])
    if tskip is None:
        target = 60
        ts = np.int16(np.timedelta64(np.diff(dt0)[0], 's'))
        tskip = int(target/ts)
    dt = dt0[::tskip]
    stec = np.nan * np.zeros((dt.size, D.sv.size))
    F = np.nan * np.zeros((dt.size, D.sv.size))
    sb = np.ones(D.sv.size)
    if return_aer:
        AER = np.nan * np.zeros((dt0.size, D.sv.size, 3))
        
    if return_mapping_function:
        F0 = np.nan * np.zeros((dt0.size, D.sv.size))
    for isv, sv in enumerate(D.sv.values):
        aer = getIonosphericPiercingPoints(D.position, sv, dt, ipp_alt=H, 
                                           navfn=fsp3, cs='aer', 
                                           rx_xyz_coords='xyz')
        idel = (aer[1] >= el_mask)
        F[idel, isv] = getMappingFunction(aer[1][idel], h=H)
        if return_mapping_function or return_aer:
            aer0 = getIonosphericPiercingPoints(D.position, sv, dt0, ipp_alt=H, 
                                           navfn=fsp3, cs='aer', 
                                           rx_xyz_coords='xyz')
            idel0 = (aer0[1] >= el_mask)
        stec[idel, isv] = getPhaseCorrTEC(L1=D.L1.values[::tskip][idel,isv], L2=D.L2.values[::tskip][idel,isv],
                                     P1=D.C1.values[::tskip][idel,isv], P2=D.P2.values[::tskip][idel,isv],
                                     maxgap=maxgap, maxjump=maxjump)
#        penalty = np.nan * np.copy(stec)
#        for i in range(stec.shape[1]):
#            penalty[:,i] = scintillation.sigmaTEC(stec[:,i], 30)
#        penalty = np.nansum(penalty**2, axis=1)
#        idp = (penalty < 1)
        if jplg_file is not None:
            sb[isv] = getSatBias(jplg_file, sv)
        if return_mapping_function:
            F0[idel0, isv] = getMappingFunction(aer0[1][idel0], h=H)
        if return_aer:
            AER[idel0, isv, 0] = aer0[0][idel0]
            AER[idel0, isv, 1] = aer0[1][idel0]
            AER[idel0, isv, 2] = aer0[2][idel0]
    # LEAST sQUARES FIT
    x0 = np.empty(D.sv.size) if jplg_file is None else sb
#    sb_lsq = least_squares(_fun, x0, args=(stec, penalty, F), xtol=1e-5, gtol=1e-5, 
#                       diff_step=1e-4, loss='cauchy')
    sb_lsq = least_squares(_fun, x0, args=(stec, F), loss='soft_l1')
    D.close()
    
    if return_mapping_function and (not return_aer):
        return sb_lsq.x, F0
    elif return_mapping_function and return_aer:
        return sb_lsq.x, F0, AER
    elif (not return_mapping_function) and return_aer:
        return sb_lsq.x, AER
    else:
        return sb_lsq.x

def getCNR(D, fsp3=None, el_mask=30, H=350):
    if isinstance(D, str):
        D = gr.load(D)
    
    time = D.time.values
    try:
        CNO = D.S1.values

        if fsp3 is not None:
            assert os.path.exists(fsp3)
            dt = np.array([np.datetime64(ttt) for ttt in time]).astype('datetime64[s]').astype(datetime)
            for isv, sv in enumerate(D.sv.values):
                aer = getIonosphericPiercingPoints(D.position, sv, dt, 
                                                          ipp_alt=H, navfn=fsp3,
                                                          cs='aer', rx_xyz_coords='xyz')
                idel = (aer[1] >= el_mask)
                CNO[~idel, isv] = np.nan
    except:
        CNO = None
    return CNO

def getDTEC(fnc, fsp3, el_mask=30, maxjump=1.6, maxgap=1, eps=1, tsps=30):
    
    vtec, aer = getVTEC(fnc, fsp3, return_aer=True)
    dtec = np.nan * np.copy(vtec)
    
    el_mask0 = (el_mask - 20) if (el_mask - 20) > 5 else 5
    
    for i in range(vtec.shape[1]):
        idel0 = aer[:, i, 1] <= el_mask0
        vtec[idel0, i] = np.nan
        idx, intervals = getIntervalsTEC(vtec[:, i], maxgap=maxgap, maxjump=maxjump)
        dtec[:,i] = tecdPerLOS(vtec[:, i], intervals, eps=eps, tsps=tsps)
        idel = aer[:, i, 1] < el_mask
        dtec[idel, i] = np.nan
    
    return dtec

def getDTECra(fnc, fsp3, el_mask=30, maxjump=1.6, maxgap=1, N=40):
    
    vtec, aer = getVTEC(fnc, fsp3, return_aer=True)
    dtec = np.nan * np.copy(vtec)
    
    el_mask0 = (el_mask - 20) if (el_mask - 20) > 5 else 5
    
    for i in range(vtec.shape[1]):
        idel0 = aer[:, i, 1] <= el_mask0
        vtec[idel0, i] = np.nan
        dtec[:,i] = uf.detrend_running_mean(vtec[:,i], N=N)
        idel = aer[:, i, 1] < el_mask
        dtec[idel, i] = np.nan
    
    return dtec

def getDTEC2(vtec, maxgap=1, maxjump=1, eps=1, tsps=30):
    
    dtec = np.nan * np.copy(vtec)
    
    for i in range(vtec.shape[1]):
        idx, intervals = getIntervalsTEC(vtec[:, i], maxgap=maxgap, maxjump=maxjump)
        dtec[:,i] = tecdPerLOS(vtec[:, i], intervals, eps=eps, tsps=tsps)
    
    return dtec