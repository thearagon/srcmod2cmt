#!/usr/bin/env python
"""
    By Thea Ragon
    Geoazur
    ragon@geoazur.unice.fr
    
    Created      : Jul 22, 2016
    Last modified: Jul 22, 2016


    Translate finite source model in SRCMOD .fsp format to CMTSOLUTION format
"""


# Import Python Libraries
import numpy as np
import re

def fh(lon, lat, z, strike, distance):
    """
    
    Given a start point (lon lat), bearing (degrees), and distance (m),
    calculates the destination point (lon lat)

    """
    
    theta = strike
    delta = distance / 6371000.

    theta = theta * np.pi / 180.
    lat1 = lat * np.pi / 180.
    lon1 = lon * np.pi / 180.

    lat2 = np.arcsin( np.sin(lat1) * np.cos(delta) + \
     np.cos(lat1) * np.sin(delta) * np.cos(theta) )

    lon2 = lon1 + np.arctan2( np.sin(theta) * np.sin(delta) * np.cos(lat1), \
     np.cos(delta) - np.sin(lat1) * np.sin(lat2))

    return (lon2 * 180.0 / np.pi, lat2 * 180.0 / np.pi, z)


def fd(lon, lat, z, strike, dip, ddip):
    """
    
    Given a start point (lon lat z), strike, dip of the fault (degrees),
    and distance along dip (m), calculates the destination point (lon lat z)

    """
   
    theta = strike + 90
        
    z2 = z + ddip * np.sin( dip * np.pi / 180. )
    D = np.cos( dip * np.pi / 180. ) * ddip
   
    return fh(lon, lat, z2, theta, D)
    

def srcmod2cmtsol(file, cmtsol_1st_l, vel_rup = [], outfile = "", onefile = 1):
    
    '''
    Translate finite source model in SRCMOD .fsp format to CMTSOLUTION format
    
    SRCMOD database: http://equake-rc.info/SRCMOD/
    find CMTSOLUTION: http://www.globalcmt.org/CMTsearch.html
    
    
    IN ARGUMENT:
    You need to specify at least 2 arguments:
     --> file = SRCMOD eventTAG (string)
                [or] SRCMOD .fsp url (string)
                [or] path+filename for a local .fsp file (string)
          ex: s2011TOHOKU01WEIx
              [or] http://equake-rc.info/media/srcmod/_fsp_files/s2011TOHOKU01WEIx.fsp
              [or] ./s2011TOHOKU01WEIx.fsp
             
     --> cmtsol_1st_l = 1st line of cmt solution to write, string !!WITHOUT 1ST SPACE!!
          ex: 'PDEW2009  4  6  1 32 39.00  42.3300   13.3300   8.8 5.9 6.3 CENTRAL ITALY'
    
    OPTIONAL ARGUMENTS:
     --> vel_rup = rupture velocity in km/s (float). Default is None, if needed it is asked in the script.
     --> outfile = default is ./ + eventTAG +_cmtsol.dat (string)
     --> onefile = 1 (write to a single file -- default)
                   [or] 0 (write to individual files)

    USAGE:
    >>> import cmtsol
    >>> cmtsol.srcmod2cmtsol('s2009LAQUIL', '04CIRE', cmtsol_1st_l = 'PDEW2009  4  6  1 32 39.00  42.3300   13.3300   8.8 5.9 6.3 CENTRAL ITALY', event_name = '200904060132A')
            
    >>> cmtsol.srcmod2cmtsol('s2011TOHOKU', '01WEIx', cmtsol_1st_l = 'PDE 2011  3 11  5 46 23.00  38.3215  142.3693  24.4 8.8 8.8 P', event_name = '201103110546A')
    
    '''
    
    # ---------------------------------------------------------------------------
    # Get path to file
    if 'http' not in file and '/' not in file and '.' not in file:
        import urllib2
        data = urllib2.urlopen("http://equake-rc.info/media/srcmod/_fsp_files/{}.fsp".format(file)).read()
        print('\n \nData imported from http://equake-rc.info/media/srcmod/_fsp_files/{}.fsp'.format(file))
    elif 'http' in file:
        import urllib2
        data = urllib2.urlopen(file).read()
        print('\n \nData imported from {}'.format(file))
    else:
        data = open(file,'r').read()
        print('\n \nData imported from local file {}'.format(file))

    # ---------------------------------------------------------------------------
    # Extract data from file
    
    # Fault parameters
    param = ['LAT  = ', 'LON = ', 'DEP = ', 'LEN  = ', 'WID =  ', 'Mw = ', \
            'STRK = ', 'DIP = ', 'RAKE = ', 'Htop = ', 'HypX = ', \
            'HypZ = ', 'Nx  =  ', 'Nz  = ']
    
    res = []
    
    for pname in param:
        pattern = re.compile(pname+'([\-]{0,1}\d+[\.]{0,1}\d*)')
        
        for j in data.split("\t"):
            match = pattern.search(j)
            if match:
                 res.append(match.groups()[0])
                 break
    LAT = float(res[param.index('LAT  = ')])
    LON = float(res[param.index('LON = ')])
    DEP = float(res[param.index('DEP = ')])
    LEN = float(res[param.index('LEN  = ')])
    WID = float(res[param.index('WID =  ')])
    Mw = float(res[param.index('Mw = ')])
    STRK = np.float(res[param.index('STRK = ')])
    DIP = float(res[param.index('DIP = ')])
    RAKE = np.float(res[param.index('RAKE = ')])
    Htop = float(res[param.index('Htop = ')])
    HypX = float(res[param.index('HypX = ')]) # rupture nucleation point in fault-plane coordinates (starting at top-left corner)
    HypZ = float(res[param.index('HypZ = ')])
    Nx = float(res[param.index('Nx  =  ')]) # number of patches along strike
    Nz = float(res[param.index('Nz  = ')]) # number of patches along dip
       
    # Get seismic moment (scientific notation)
    pattern = re.compile('Mo = [\-]{0,1}\d+[\.]{0,1}\d*([e][+]\d+)')
    for j in data.split("\t"):
        match = pattern.search(j)
        if match:
            M = match.group()
    Mo = float(M.split()[2])
     
    # Get date
    pattern = re.compile('\d+/\d+/\d+')
    for j in data.split("\t"):
        match = pattern.search(j)
        if match:
            D = match.group()
    dd = D.split('/')[0]
    mm = D.split('/')[1]
    yyyy = D.split('/')[2]
    
    # Get event TAG
    pattern = re.compile('EventTAG: [\D][\d]{4}[\D]{6}[\d]{2}[\D]{4}')
    for j in data.split("\t"):
        match = pattern.search(j)
        if match:
            M = match.group()
    eventTAG = M.split()[1]
    
    ## Read SOURCE MODEL PARAMETERS of SRCMOD FSP file
    # Get coordinates of each patch top center, lon lat (deg), Z (km)
    tc_lat = []
    tc_lon = []
    tc_x = []
    tc_y = []
    tc_z = []
    
    # Get rake (deg), slip (m), rise (s) and trup (s)
    rake = []
    slip = []
    rise = []
    trup = []
    
    line = []
    line_0 = []
    
    # Read tab
    for i in data.split('\n'): # remove header
        if i.startswith('%    LAT       LON'):   # Get line with SOURCE MODEL PARAMETERS header --> line0      
            line_0.append(i)
        elif not i.startswith('%'):
            line.append(i)
       
    # Read lines
    line0 = line_0[0].split()   # first line with headers
    for j in range(len(line)-1): 
        splitline = line[j].split()
        
        if 'LAT' in line0:
            index = line0.index('LAT') - 1      # 1st charac of line0 is %
            tc_lat.append(np.float(splitline[index]))
        if 'LON' in line0:
            index = line0.index('LON') - 1      
            tc_lon.append(np.float(splitline[index]))
        if 'X==EW' in line0:
            index = line0.index('X==EW') - 1      
            tc_x.append(np.float(splitline[index]))
        if 'Y==NS' in line0:
            index = line0.index('Y==NS') - 1      
            tc_y.append(np.float(splitline[index]))
        if 'Z' in line0:
            index = line0.index('Z') - 1      
            tc_z.append(np.float(splitline[index]))
        if 'RAKE' in line0:
            index = line0.index('RAKE') - 1      
            rake.append(np.float(splitline[index]))
        if 'SLIP' in line0:
            index = line0.index('SLIP') - 1      
            slip.append(np.float(splitline[index]))
        if 'RISE' in line0:
            index = line0.index('RISE') - 1      
            rise.append(np.float(splitline[index]))
        if 'TRUP' in line0:
            index = line0.index('TRUP') - 1      
            trup.append(np.float(splitline[index]))
    
    
    
    # ---------------------------------------------------------------------------
    # Build subftaults point source coordinates from patches top center points
    
    lon = []
    lat= []
    z = []
    
    for i in range(len(tc_lat)):
        center = fd(tc_lon[i], tc_lat[i], tc_z[i], STRK, DIP, WID/Nz/2)
        lon.append(center[0])
        lat.append(center[1])
        z.append(center[2])
    
    
    # ---------------------------------------------------------------------------
    # Calculate moments for each subfault
    
    
    # Seismic moment
    cond = 0
    for i in range(len(lon)):
        if slip[i] < 0:
            slip[i] = 0
            cond = 1
    if cond == 1:
        print('\nSome slip values are < 0  -->  have been set to 0')
            
    
    S_fault = LEN*10**3 * WID*10**3     # LEN & WID in km
    mean_slip = np.mean(slip)       # slip in m
    
    mu = Mo / (S_fault * mean_slip)     # Mo in N.m
    
    
    S_sub = ((LEN*10**3)/Nx) * ((WID*10**3)/Nz)
    M0 = []
    for i in range(len(lon)):
        M0.append(mu * S_sub * slip[i])
        
        
    # Moment tensor
    d2r = np.pi/180.0
    S = STRK * d2r
    D = DIP * d2r
    
    sinS = np.sin(S)
    cosS = np.cos(S)
    sin2S = np.sin(2*S)
    cos2S = np.cos(2*S)
    sinD = np.sin(D)
    cosD = np.cos(D)
    sin2D = np.sin(2*D)
    cos2D = np.cos(2*D)
    
    M11 = np.empty(len(lon))
    M22 = np.empty(len(lon))
    M33 = np.empty(len(lon))
    M12 = np.empty(len(lon))
    M13 = np.empty(len(lon))
    M23 = np.empty(len(lon))
    

    for i in range(len(lon)):       # for each subfault
        
        if rake != []:
            R = rake[i] *d2r
        else:
            R = RAKE *d2r
            
        sinR = np.sin(R)
        cosR = np.cos(R)        
        
        # Moments
        M11[i] = (-sinS**2*sinR*sin2D-sin2S*cosR*sinD) *M0[i]*10**7
        M22[i] = (-cosS**2*sinR*sin2D+sin2S*cosR*sinD) *M0[i]*10**7
        M12[i] = (0.5*sin2S*sinR*sin2D+cos2S*cosR*sinD) *M0[i]*10**7
        M13[i] = (-cosS*cosR*cosD-sinS*sinR*cos2D) *M0[i]*10**7
        M23[i] = (cosS*sinR*cos2D-sinS*cosR*cosD) *M0[i]*10**7
        M33[i] = (sin2D*sinR) *M0[i]*10**7
    
    
    
    # ---------------------------------------------------------------------------
    # Calculate time shift and half duration
    
    half_dur = np.empty(len(lon))
    for i in range(len(lon)):
        half_dur[i] = 1.05e-8 * (M0[i]*10**7)**(1./3.)     # Mo in dyne.cm
    
    
    if trup != [] and [x for x in trup if x == 0] is not None: 
        time_shift = trup
    else:
        if vel_rup == []:# or rise_time == []:
            print('\n \nTime shif calculation  -->  no value in srcmod \n \nIn order to calculate it, please enter:')
        if vel_rup == []:
            vel_rup = input("---- Rupture velocity in km/s (commonly btw 2.5 and 3.5) ?  ")
        #if rise_time == []:
        #    rise_time = input("---- Rise time (s)?  \n        ex: btw 1 and 2s  :  ")
        
        time_shift = []
        
        for i in range(len(lon)):# for each patch
            # calculate distance between hypocentre and patch midle center
            d = np.tan(np.pi/2 - DIP *d2r) * (tc_z[i]-z[i]) # horizontal dist btw tc and mc
            H = np.sqrt( (tc_x[i] + d* np.cos(STRK*d2r))**2 + (tc_y[i] + d* np.sin(STRK*d2r))**2 ) # horiz dist btw mc and epicentre
            Phyp = HypZ * np.sin(DIP *d2r)      # Hypocenter depth
            D = np.sqrt( H**2 + (Phyp - z[i])**2 ) # dist btw hypocentre and mc
            
            # calculate time_shift
            T = D / vel_rup
            time_shift.append(T)
        
        m = min(time_shift)
        time_shift = time_shift - m   # Specfem need at least 1 time_shift set to zero
      
    
    
    # ---------------------------------------------------------------------------
    # Write CMTSOL file
    
    
    if outfile == '':    
        outfile = './'+eventTAG+'_cmtsol.dat'
        print('\nNo outfile name specified  -->  outfile is {}'.format(outfile))
    else: print('\n  outfile is {}'.format(outfile))
    

    if onefile == 1:
        
        f = open(outfile, 'w')
        
        for i in range(len(lon)):       # for each subfault
            
            f.write(" %s\n" % (cmtsol_1st_l) )
            f.write('event name: {:11s}_{}\n'.format(eventTAG, i+1))
            f.write('time shift: {:11.4f}\n'.format(time_shift[i]))  
            f.write('half duration: {:8.4f}\n'.format(half_dur[i]))
            f.write('latitude: {:13.4f}\n'.format(lat[i]))
            f.write('longitude: {:12.4f}\n'.format(lon[i]))
            f.write('depth: {:16.4f}\n'.format(z[i]))
            f.write('Mrr: {:18.6g}\n'.format(M33[i]))
            f.write('Mtt: {:18.6g}\n'.format(M11[i]))
            f.write('Mpp: {:18.6g}\n'.format(M22[i]))
            f.write('Mrt: {:18.6g}\n'.format(M13[i]))
            f.write('Mrp: {:18.6g}\n'.format(-M23[i]))
            f.write('Mtp: {:18.6g}\n'.format(-M12[i]))
            
        f.close()
        
        
    if onefile == 0:
    
        for i in range(len(lon)):       # for each subfault
            
            f = open(outfile+'{}'.format(i+1), 'w')
            
            f.write(" %s\n" % (cmtsol_1st_l) )
            f.write('event name: {:11s}_{}\n'.format(eventTAG, i+1))
            f.write('time shift: {:11.4f}\n'.format(time_shift[i]))  
            f.write('half duration: {:8.4f}\n'.format(half_dur[i]))
            f.write('latitude: {:13.4f}\n'.format(lat[i]))
            f.write('longitude: {:12.4f}\n'.format(lon[i]))
            f.write('depth: {:16.4f}\n'.format(z[i]))
            f.write('Mrr: {:18.6g}\n'.format(M33[i]))
            f.write('Mtt: {:18.6g}\n'.format(M11[i]))
            f.write('Mpp: {:18.6g}\n'.format(M22[i]))
            f.write('Mrt: {:18.6g}\n'.format(M13[i]))
            f.write('Mrp: {:18.6g}\n'.format(-M23[i]))
            f.write('Mtp: {:18.6g}\n'.format(-M12[i]))
            
            f.close()
    
    return
    
    
    
    
if __name__ == '__main__':
    
    '''
    Translate finite source model in SRCMOD .fsp format to CMTSOLUTION format
    
    SRCMOD database: http://equake-rc.info/SRCMOD/
    find CMTSOLUTION: http://www.globalcmt.org/CMTsearch.html
    
    
    IN ARGUMENT:
    You need to specify at least 2 arguments:
     --> file = SRCMOD eventTAG (string)
                [or] SRCMOD .fsp url (string)
                [or] path+filename for a local .fsp file (string)
          ex: s2011TOHOKU01WEIx
              [or] http://equake-rc.info/media/srcmod/_fsp_files/s2011TOHOKU01WEIx.fsp
              [or] ./s2011TOHOKU01WEIx.fsp
             
     --> cmtsol_1st_l = 1st line of cmt solution to write, string !!WITHOUT 1ST SPACE!!
          ex: 'PDEW2009  4  6  1 32 39.00  42.3300   13.3300   8.8 5.9 6.3 CENTRAL ITALY'
    
    OPTIONAL ARGUMENTS:
    You need to specify ALL optional arguments with their keyword.
     --> vel_rup = rupture velocity in km/s (float). Default is None, if needed it is asked in the script.
     --> outfile = default is ./ + eventTAG +_cmtsol.dat (string)
     --> onefile = 1 (write to a single file -- default)
                   [or] 0 (write to individual files)

    USAGE:
    > python cmtsol.py s2009LAQUIL04CIRE 'PDEW2009  4  6  1 32 39.00  42.3300   13.3300   8.8 5.9 6.3 CENTRAL ITALY'
           
    > python cmtsol.py 'http://equake-rc.info/media/srcmod/_fsp_files/s2011TOHOKU01WEIx.fsp' 'PDE 2011  3 11  5 46 23.00  38.3215  142.3693  24.4 8.8 8.8 P' vel_rup=1.2
    
    > python cmtsol.py './s2011TOHOKU01HAYE.fsp' 'PDE 2011  3 11  5 46 23.00  38.3215  142.3693  24.4 8.8 8.8 P' 'vel_rup =1.2' 'outfile = cmtsol.dat' onefile=1
    
    '''
    
    import sys
      
    # Define default values
    vel = []
    out = ""
    one = 1
    
    if len(sys.argv) == 3:
        srcmod2cmtsol(sys.argv[1], sys.argv[2], vel_rup = [], outfile = "", onefile = 1)
    elif len(sys.argv) >= 3:
        for i in [x for x in range(len(sys.argv)) if x>=3]:
            if 'vel_rup' in sys.argv[i]:
                vel = float(sys.argv[i].split('=')[1])
            if 'outfile' in sys.argv[i]:
                out = sys.argv[i].split('=')[1]
            if 'onefile' in sys.argv[i]:
                one = float(sys.argv[i].split('=')[1])
            if 'vel_rup' not in sys.argv[i] and 'outfile' not in sys.argv[i] and 'onefile' not in sys.argv[i]:
                print('\nYou need to specify ALL optional arguments with their keyword.  ex: vel_rup=2 ')
                exit()
        srcmod2cmtsol(sys.argv[1], sys.argv[2], vel_rup = vel, outfile = out, onefile = one )
        
    else: 
        print("\nYou need to specify at least 2 arguments: \n   --> event tag or url from SRCMOD or .fsp local file \n   --> 1st line of CMT solution \
                \n \n USAGE EXEMPLE: \n> python cmtsol.py s2009LAQUIL04CIRE 'PDEW2009  4  6  1 32 39.00  42.3300   13.3300   8.8 5.9 6.3 CENTRAL ITALY'\
                \n> python cmtsol.py 'http://equake-rc.info/media/srcmod/_fsp_files/s2011TOHOKU01WEIx.fsp' 'PDE 2011  3 11  5 46 23.00  38.3215  142.3693  24.4 8.8 8.8 P' vel_rup=1.2\
                \n> python cmtsol.py './s2011TOHOKU01HAYE.fsp' 'PDE 2011  3 11  5 46 23.00  38.3215  142.3693  24.4 8.8 8.8 P' 'vel_rup =1.2' 'outfile = cmtsol.dat' onefile=1")
        exit()

