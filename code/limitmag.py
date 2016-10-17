import numpy as np
import sys
import os

# Dictionaries of wavelengths(um), dlam (um), and sky background (mag) for each filter  
# Sky background in AB magnitudes (converted, original data from Cambridge Astronomy 
# Survey Unit and Cuby et al. 2000)
lamdict  = {'Y': 1.025, 'J': 1.22, 'H': 1.63,  'K': 2.193, 'JH': 1.445}
dlamdict = {'Y': 0.11,  'J': 0.22, 'H': 0.3,   'K': 0.39,  'JH': 0.67}
bgdmdict = {'Y': 18.18, 'J': 16.61, 'H': 14.99, 'K': 14.85, 'JH': 15.78}

dir_path = os.path.dirname(os.path.realpath(__file__))+'/'

# Calculated from filtereff(), lowresspeceff(), highresspeceff()
photoeff = {'Y': 0.447, 'J': 0.476, 'H': 0.484, 'K': 0.440}
lowreff  = {'Y': 0.130, 'J': 0.130, 'H': 0.133, 'K': 0.133}
highreff = {'Y': 0.134, 'J': 0.134, 'H': 0.065, 'K': 0.065}

def limit(detector, band, R, exptime,signos):

	"""
	***** VERIFIED PHOTOMETRY WITH RATIR AND GROND                  *****
	***** VERIFIED SPECTROSCOPY WITH CHART IN FILE BELOW            *****
	***** http://www.astro.wisc.edu/twiki/pub/Salt/RobertStobieSpectrograph/RSS_NIR_Sheinis_SPIE.pdf *****
	
	NAME:
	    limit
	PURPOSE:
	    Calculate limiting magnitudes given mode, exposure time, and SN
	INPUTS:
	    detector - detector to use [insb, h2rg]
	    band     - filter to calculate over [Y, J, H, K, JH]
	    R        - resolving power (0 for photometry, can handle anything above)
	    exptime  - exposure time in seconds
	    signos   - signal-to-noise ratio want to achieve
	OUTPUTS:
	    Limiting magnitude in AB
	EXAMPLE:
	    limit('h2rg', 'Y', 0, 400, 10)
	NOTES:
	    Only calculates the central wavelength of filter band for spectroscopy,
	    will eventually want to change this to match wavelength range
	    Need better efficiency measurements for spectroscopy modes
	"""
	
	if band.upper() not in ['Y','J', 'H', 'K', 'JH']:
		sys.exit('Please enter a valid photometric band (Y,J,H,K,JH)')

	if detector.lower() not in ['insb','h2rg']:
		sys.exit('Only configured for InSb and H2RG please enter a valid detector')

    # Actual DCT mirror data (nm and transmission perc.)
	skw, skp = np.loadtxt(dir_path+'skyback/skyback_900_2400_test.txt', unpack=True)
	
	# Upper aperture 1120mm, outer aperture 4280mm
	d_tel   = 428.0               # centimeters 
	f_num   = 2.47                # Actual telescope is 6.1, but demagnification on detector is less  
	aeff    = 0.93                # telescope effective collecting area
	obj_npix= 4.                  # number of pixel objects is smeared over  

	if detector.lower() == 'insb':
		pix     = 30.0            # pixel size in microns 
		ns_read = 30.             # detector read noise in e-
		drk     = 1.0             # electrons per second 
		eta     = filtereff_slitviewer() # overall instrument efficiency	
	
	if detector.lower() == 'h2rg':
		pix     = 18.0            # pixel size in microns 
		ns_read = 6.              # detector read noise in e-
		drk     = 0.005           # electrons per second 
		eta 	= photoeff[band] # overall instrument efficiency
	
	lam   = lamdict[band.upper()]	
	dlam  = dlamdict[band.upper()]
	bgd_m = bgdmdict[band.upper()]

	hi = lam + dlam/2.
	lo = lam - dlam/2.
	
	if R > 0:
		#lam = lam+dlam/2.0
				
		if detector.lower() == 'h2rg' and R == 30:
		    if (band.upper() == 'Y') or (band.upper() == 'J') :
			    eta = lowreff[band.upper()]
			    hi = (lamdict['J'] + dlamdict['J']/2.) 
			    lo = (lamdict['Y'] - dlamdict['Y']/2.)
			    lam = np.mean([lo,hi])
			    range = hi-lo
			    keep = np.where( (skw/1000. >= lo) & (skw/1000. <= hi) )

		    if (band.upper() == 'H') or (band.upper() == 'K') :
			    eta = lowreff[band.upper()]
			    hi = (lamdict['K'] + dlamdict['K']/2.) 
			    lo = (lamdict['H'] - dlamdict['H']/2.)
			    lam = np.mean([lo,hi])
			    range = hi-lo
			    keep = np.where( (skw/1000. >= lo) & (skw/1000. <= hi) )

		if detector.lower() == 'h2rg' and R == 4000:
		    if (band.upper() == 'Y') or (band.upper() == 'J') :
			    eta = highreff[band.upper()]
			    hi = (lamdict['J'] + dlamdict['J']/2.) 
			    lo = (lamdict['Y'] - dlamdict['Y']/2.)
			    lam = np.mean([lo,hi])
			    range = hi-lo
			    keep = np.where( (skw/1000. >= lo) & (skw/1000. <= hi) )
			    
		    if (band.upper() == 'H') or (band.upper() == 'K') :
			    eta = highreff[band.upper()]
			    hi = (lamdict['K'] + dlamdict['K']/2.) 
			    lo = (lamdict['H'] - dlamdict['H']/2.)
			    lam = np.mean([lo,hi])
			    range = hi-lo
			    keep = np.where( (skw/1000. >= lo) & (skw/1000. <= hi) )
	                     
		dlam = lam/R		
		skw = skw[keep]
		skp = skp[keep]
	else:
	    keep = np.where( (skw/1000. >= (lamdict[band.upper()] - dlamdict[band.upper()]/2.)) &\
	                     (skw/1000. <= (lamdict[band.upper()] + dlamdict[band.upper()]/2.)) )
	    skw = skw[keep]
	    skp = skp[keep]
	    
	# Make floats
	t_int   = 1.0*exptime         # integration time in seconds
	SN      = 1.0*signos          # limiting S/N ratio to do calculation for 

    # Calculate telescope area, effective focal length, pixel size in " and sr
	asrad   = 206264.8            # arcseconds per radian 
	a_tel   = np.pi/4*d_tel**2    # telescope area in cm^2
	EFL     = d_tel*f_num         # telescope effective focal length 
	pixas   = pix*1e-4/EFL*asrad  # pixel in arcseconds 
	p_omega = (pixas/asrad)**2    # pixel solid angle in steradians  
	
    # Just prints assumptions in easily readable format
	#print '  '
	#print  'ASSUMPTIONS: '
	#print  'T integration ', t_int 
	#print  'tel dia(cm) F/#   pixsize(um),  pixsize("),  scndry vgntn, instr effcy' 
	#print   d_tel, '     ', f_num,'',pix, '        ','%.3f'%(pixas), '      ',aeff, '        ',eta
	#print  'dark current ', drk
	#print  'read noise ', ns_read 
	#print  'objects spread over # pixels: ', obj_npix

	# Convert background flux from AB to photons/cm^2/s/um 
	# (extra factor of 1e-4 to convert m^2 to cm^2)
	#f_bgd_ph = i_conv(bgd_m, 'AB', 'ph', lam)*1e-4 
	
	# In ph/s/arcsec^2/nm/m^2
	#f=interpolate.interp1d(skw/1000., skp,kind='linear')
	
	if R==0:
	    x = np.arange(lo+0.001, hi, 0.001)
	else:
	    x = np.arange(lo+dlam, hi, dlam)
	    
	f = np.interp(x, skw/1000.,skp)

	# (extra factor of 1e-4 to convert m^2 to cm^2, and 1000 from nm to um)
	f_bgd_ph = f*1e-4*1000.
	
	# Sum number of total photon flux over whole bandpass
	if R == 0: f_bgd_ph = sum(f_bgd_ph*0.001)/dlam
	
	# Background signal per pixel over integration time 
	# (factors inefficiency-this includes QE, band, and area efficiency)
	# In electrons because eta includes QE
	# Turns into number of electrons per pixel
	arcsec = 1./pixas
	s_bgd_epix  = f_bgd_ph * dlam * a_tel * t_int * eta * aeff * 1.4*arcsec**2
	
	# NOISE: 
	ns_drk = np.sqrt(drk*t_int)                                 # dark noise (Poisson, e-)
	ns_det = np.sqrt(ns_drk**2+ns_read**2) * np.sqrt(obj_npix)  # total detector noise per total object area   
	ns_bgn = np.sqrt(s_bgd_epix * obj_npix)                     # photon noise from the background

	#print '  '
	#print 'NOISE:'
	#print 'S_bgnd: %.2f'%(np.median(s_bgd_epix))
	#print 'ns_bgn ns_det' 
	#print '%.2f  %.2f'%(np.median(ns_bgn), ns_det) 
	
	# noise background plus detector in quadrature
	ns_bgd_det   = np.sqrt(ns_bgn**2+ ns_det**2)
	
	# signal corresponding to the limiting magnitude at the specified S/N ratio assuming Poisson
	s_lmag_epix  =  0.5*(SN**2+ np.sqrt( SN**4 + 4*SN**2*ns_bgd_det**2))  
		
	# limiting object flux : photons/cm^2/s/um (dividing e- by QE, so in photons)
	# Make into flux by dividing e- by area of telescope, exposure time, dlam
	# (factoring in inefficiencies) then convert to AB mag
	f_lmag_pcmum = s_lmag_epix/(a_tel * dlam * t_int * eta * aeff)
	
	f_lmag_AB    = i_conv(f_lmag_pcmum*1e4, 'ph', 'AB', lam)  # 1e4 to convert to per m^2  
	
	if R != 0: 
	    f_lmag_AB = i_conv(sum(f_lmag_pcmum)*dlam/range*1e4, 'ph', 'AB', lam)
	
	# Print limiting magnitude of
	#print '    '
	#print 'limiting magnitude AB ', np.mean(f_lmag_AB), ' at signal to noise ratio of ', SN
	return np.mean(f_lmag_AB),eta


def exptimecalc(detector, band, R, objmag, signos):

	"""
	***** WORKS FOR PHOTOMETRY, NEED TO ADAPT FOR INDIVIDUAL DEL LAMBDA FOR SPECTROSCOPY *****
	NAME:
	    exptimecalc
	PURPOSE:
	    Calculate exposure time given mode, object magnitude (AB), and SN
	INPUTS:
	    detector - detector to use [insb, h2rg]
	    band     - filter to calculate over [Y, J, H, K, JH]
	    R        - resolving power (0 for photometry, can handle anything above)
	    objmag   - point source magnitude in AB
	    signos   - signal-to-noise ratio want to achieve
	OUTPUTS:
	    Exposure time in seconds
	EXAMPLE:
	    exptimecalc('h2rg', 'Y', 0, 23., 10)
	NOTES:
	    Only calculates the central wavelength of filter band for spectroscopy,
	    will eventually want to change this to match wavelength range
	    Need better efficiency measurements for spectroscopy modes
	"""
	
	if band.upper() not in ['Y','J', 'H', 'K', 'JH']:
		sys.exit('Please enter a valid photometric band (Y,J,H,K,JH)')

	if detector.lower() not in ['insb','h2rg']:
		sys.exit('Only configured for InSb and H2RG please enter a valid detector')

    # Actual DCT mirror data (nm and transmission perc.)
	skw, skp = np.loadtxt(dir_path+'skyback/skyback_900_2400_test.txt', unpack=True)
	
	# Upper aperture 1120mm, outer aperture 4280mm
	d_tel   = 428.0               # centimeters 
	f_num   = 2.47                # Actual telescope is 6.1, but demagnification on detector is less  
	aeff    = 0.93                # telescope effective collecting area
	obj_npix= 4.                  # number of pixel objects is smeared over  

	if detector.lower() == 'insb':
		pix     = 30.0            # pixel size in microns 
		ns_read = 30.             # detector read noise in e-
		drk     = 1.0             # electrons per second 
		eta     = filtereff_slitviewer() # overall instrument efficiency	
	
	if detector.lower() == 'h2rg':
		pix     = 18.0            # pixel size in microns 
		ns_read = 6.              # detector read noise in e-
		drk     = 0.005           # electrons per second 
		eta 	= photoeff[band] # overall instrument efficiency
	
	lam   = lamdict[band.upper()]	
	dlam  = dlamdict[band.upper()]
	bgd_m = bgdmdict[band.upper()]

	hi = lam + dlam/2.
	lo = lam - dlam/2.
	
	if R > 0:
		#lam = lam+dlam/2.0
				
		if detector.lower() == 'h2rg' and R == 30:
		    if (band.upper() == 'Y') or (band.upper() == 'J') :
			    eta = lowreff[band.upper()]
			    hi = (lamdict['J'] + dlamdict['J']/2.) 
			    lo = (lamdict['Y'] - dlamdict['Y']/2.)
			    lam = np.mean([lo,hi])
			    range = hi-lo
			    keep = np.where( (skw/1000. >= lo) & (skw/1000. <= hi) )

		    if (band.upper() == 'H') or (band.upper() == 'K') :
			    eta = lowreff[band.upper()]
			    hi = (lamdict['K'] + dlamdict['K']/2.) 
			    lo = (lamdict['H'] - dlamdict['H']/2.)
			    lam = np.mean([lo,hi])
			    range = hi-lo
			    keep = np.where( (skw/1000. >= lo) & (skw/1000. <= hi) )

		if detector.lower() == 'h2rg' and R == 4000:
		    if (band.upper() == 'Y') or (band.upper() == 'J') :
			    eta = highreff[band.upper()]
			    hi = (lamdict['J'] + dlamdict['J']/2.) 
			    lo = (lamdict['Y'] - dlamdict['Y']/2.)
			    lam = np.mean([lo,hi])
			    range = hi-lo
			    keep = np.where( (skw/1000. >= lo) & (skw/1000. <= hi) )
			    
		    if (band.upper() == 'H') or (band.upper() == 'K') :
			    eta = highresspeceff(band.upper())
			    hi = (lamdict['K'] + dlamdict['K']/2.) 
			    lo = (lamdict['H'] - dlamdict['H']/2.)
			    lam = np.mean([lo,hi])
			    range = hi-lo
			    keep = np.where( (skw/1000. >= lo) & (skw/1000. <= hi) )
	                     
		dlam = lam/R		
		skw = skw[keep]
		skp = skp[keep]
	else:
	    keep = np.where( (skw/1000. >= (lamdict[band.upper()] - dlamdict[band.upper()]/2.)) &\
	                     (skw/1000. <= (lamdict[band.upper()] + dlamdict[band.upper()]/2.)) )
	    skw = skw[keep]
	    skp = skp[keep]
			
	# Make floats
	SN      = 1.0*signos          # limiting S/N ratio to do calculation for 

    # Calculate telescope area, effective focal length, pixel size in " and sr
	asrad   = 206264.8            # arcseconds per radian 
	a_tel   = np.pi/4*d_tel**2    # telescope area in cm^2
	EFL     = d_tel*f_num         # telescope effective focal length 
	pixas   = pix*1e-4/EFL*asrad  # pixel in arcseconds 
	p_omega = (pix/asrad)**2      # pixel solid angle in steradians  

    # Just prints assumptions in easily readable format
	#print '  '
	#print  'ASSUMPTIONS: '
	#print  'Obj mag (AB) ', objmag 
	#print  'tel dia(cm) F/#   pixsize(um),  pixsize("),  scndry vgntn, instr effcy' 
	#print   d_tel, '     ', f_num,'',pix, '        ','%.3f'%(pixas), '      ',aeff, '        ',eta
	#print  'dark current ', drk
	#print  'read noise ', ns_read 
	#print  'objects spread over # pixels: ', obj_npix

	# In ph/s/arcsec^2/nm/m^2
	#f=interpolate.interp1d(skw/1000., skp,kind='linear')

	if R==0:
	    x = np.arange(lo+0.001, hi, 0.001)
	else:
	    x = np.arange(lo+dlam, hi, dlam)
	    
	f = np.interp(x, skw/1000.,skp)
	    
	# (extra factor of 1e-4 to convert m^2 to cm^2, and 1000 from nm to um)
	f_bgd_ph = f*1e-4*1000.
	
	# Sum number of total photon flux over whole bandpass
	if R == 0: f_bgd_ph = sum(f_bgd_ph*0.001)/dlam
	
	
    # Convert signal flux from AB to photons/cm^2/s/um
	f_sig_ph  = i_conv(objmag, 'AB', 'ph', lam)*1e-4  # 1e4 to convert to per cm^2 
	f_sig_e   = f_sig_ph  * (a_tel * dlam * eta * aeff)

	# Convert background flux from AB to photons/cm^2/s/um
	#f_bgd_ph  = i_conv(bgd_m, 'AB', 'ph', lam)*1e-4 
	arcsec = 1./pixas
	f_bgd_ph  *= 1.4*arcsec**2 
	f_bgd_e   = f_bgd_ph  * (a_tel * dlam * eta * aeff)
	
	quad_a = f_sig_e**2 - obj_npix*drk**2*SN**2
	quad_b = -SN**2*(f_sig_e + f_bgd_e*obj_npix)
	quad_c = -SN**2*(obj_npix*ns_read**2)
	
	exptime = (-quad_b + np.sqrt(quad_b**2 - 4*quad_a*quad_c))/(2.*quad_a)
	
	#print ' '
	#print 'exposure time (s)'
	#print '%.2f'%(np.median(exptime))
	return np.median(exptime),eta

def i_conv(i_in, i_name, o_name, lam):

	"""
	NAME:
	    i_conv
	PURPOSE:
	    Converts intensities between units
	INPUTS:
	    i_in   - input intensity
	    i_name - name of input intensity units
	    o_name - name of output intensity units
	    lam    - wavelength (um) to calculate output units
	OUTPUTS:
	    i_out  - output intensity 
	NOTES:
	    Allowed unit names for conversion ['AB', 'Jy', 'W', 'ph']
	    'W'  - W/m^2/um
	    'ph' - photons/m^2/s/um
	    'Jy' - Jansky (1Jy = 1e-26 W/m^2/Hz)
	    'AB' - AB magnitudes
	"""
	
	h=6.627e-34  #Planck constant J*s 
	c=2.997e8    # speed of light meters/s 
	hc=h*c       # h * c in J*m  

    # Convert any input to W/m^2/um
	if i_name == 'W':
		 ww = i_in                                  # Leave in same units
	if i_name == 'AB':
		 ww = 1/lam**2 * 1.116e-8 * 10**(-0.4*i_in) # convert from AB mag to W/m^2/um   
	if i_name == 'Jy':
		 ww = 2.997e-12 / lam**2 * i_in             # convert from Jy to W/m^2/um 
	if i_name == 'ph':
		 ww = 1/(5.03572e+18 * lam) * i_in          # convert from ph/m^2/s/um to W/m^2/um

    # Convert from W/m^2/um to output units
	if o_name == 'W':
		  i_out = ww                                # Leave in same units    
	if o_name == 'AB':
		  i_out = -2.5*np.log10(ww*lam**2/1.116e-8)    # convert from AB mag to W/m^2/um   
	if o_name == 'Jy':
		  i_out = lam**2/2.997e-12 * ww             # convert from W/m^2/um to Jy 
	if o_name == 'ph':
		  i_out = 5.03572e+18 * lam *ww             # convert from W/m^2/um to ph/m^2/s/um

	return i_out



################ UNNECESSARY FOR ONLINE CALCULATOR, BUT CAN REPLACE DICT ################
################ VALUES IF THERE ARE BETTER INSTRUMENT EFFICIENCY DATA   ################
################ USING THESE FUNCTIONS                                   ################

def filtereff(band, med=True):

    """
    NAME:
        filtereff
    PURPOSE:
        Finds effective efficiency of instrument from DCT primary and secondary mirror,
        fold mirror/dichroic, window, filter, dichroic, lenses, and detector
    INPUTS:
        band - photometric band (YJHK)
        med  - returns median of band, otherwise returns full array from 900-2400nm
    """
    
    # Actual RIMAS filter data (nm and transmission perc.)
    filtdict = {'Y': 'filterscurves/Y_filter_900_2400.txt',
                'J': 'filterscurves/J_filter_900_2400.txt',
                'H': 'filterscurves/H_filter_900_2400.txt',
                'K': 'filterscurves/K_filter_900_2400.txt'}
                
    # Lens transmission curves (NEED REPLACEMENT LENSES DATA) (nm and transmission perc.)
    lensdict = {'Y': 'Lenses/yj_T_data_interp.txt',
                'J': 'Lenses/yj_T_data_interp.txt',
                'H': 'Lenses/hk_T_data_interp.txt',
                'K': 'Lenses/hk_T_data_interp.txt'}
                
    # Load filter file for given band 
    filt = np.loadtxt(filtdict[band.upper()])
    
    #Initialize efficiency to 0 and keeping only the wavelengths from filter with over 50% transmission
    eta = np.copy(filt)
    ind = 0
    
    # Actual DCT mirror data (nm and transmission perc.)
    m = np.loadtxt('DCTmirrors/mirrors_900_2400.txt')
    
    # Dichroic from RIMAS theoretical model (nm, transmission perc., reflection perc.)
    d = np.loadtxt('dichroic/dich_900_2400.txt')
    
    # Kitt Peak atmospheric transmission data (nm and transmission perc.)
    t = np.loadtxt('atmosphere/atmos_900_2400.txt')
    
    l_fil = np.loadtxt(lensdict[band.upper()])
    l_col = np.loadtxt('Lenses/collimator_T_data_interp.txt')
    
    # Finds each wavelength in filter data and finds transmission of that wavelength
    # in the atmospheric transmission data and dichroic transmission.  If photometric
    # band is Y or J the light is reflected on dichroic so use reflection.
    #Multiply efficiencies to get the combined filter, dichroic, lenses, mirrors (DCT primary and secondary), and atmosphere efficiency
    for value in filt:
    
        dich = d[np.where(value[0] == d[:,0])][0]
        
        if band.upper() in ('Y', 'J'):
            dich[1] = dich[2] # For YJ use reflection, otherwise use transmission data
            
        atm = t[np.where(value[0] == t[:,0])][0]
        ref = m[np.where(value[0] == m[:,0])][0]
        
        len_fil = l_fil[np.where(value[0] == l_fil[:,0])][0]
        len_col = l_col[np.where(value[0] == l_col[:,0])][0]
        
        #Lenses
        len_fil_eta = (1e-2**5*len_fil[1]*len_fil[2]*len_fil[3]*len_fil[4]*len_fil[5])
        len_col_eta = (1e-2**5*len_col[1]*len_col[2]*len_col[3]*len_col[4]*len_col[5])
        
        feta = (dich[1]/100.)*(value[1]/100.)*(atm[1]/100.)*(ref[1]/100.)**2*len_fil_eta*len_col_eta
        
        eta[ind,1] = feta
        ind = ind + 1
        
    #feta (filter efficiency) includes atmospheric transmission, DCT primary 
    #(and secondary) reflectivity, filter transmission, dichroic transmission,
    #lens collimation transmission, optical arm lens transmission
    
    # Stuff that isn't measured precisely
    #Now include folding mirror at 95%, 2 windowed surfaces at 95%, and detector at 80%
    eta[:,1] = eta[:,1]*0.95*0.95**2*.8
    
    # Include mirror after dichroic beamsplitter for redder path assuming 95% reflectance
    if band.upper() in ('H', 'K'):
        eta[:,1] = eta[:,1]*0.95
        
    plt.plot(eta[:,0], eta[:,1])
    #plt.show()
    #plt.clf()
    
    # Either return the full array or just the median of the filter
    if med == True:
        keep = np.where(eta[:,1] != 0)
        return np.median(eta[keep,1])
    else:
        return eta	

#Finds efficiency of instrument given filter, dichroic, lenses, mirror (DCT primary and secondary), and atmospheric transmission curves
#ASSUMES efficiencies for window, folding mirror/dichroic, and detector
#for given band
#ASSUMES J filter
def filtereff_slitviewer(med=True):

	#Load filter file for given band
	filt = np.loadtxt('filterscurves/J_filter_900_2400.txt', skiprows=1)

    # Kitt Peak atmospheric transmission data (nm and transmission perc.)
	t = np.loadtxt('atmosphere/atmos_900_2400.txt')
	
    # Actual DCT mirror data (nm and transmission perc.)
	m = np.loadtxt('DCTmirrors/mirrors_900_2400.txt')
	
	# Lenses
	l_fil = np.loadtxt('Lenses/AR1.1-1.8_interp.txt')
	
	#Initialize efficiency to 0 and keeping only the wavelengths from filter with over 50% transmission
	eta = np.copy(l_fil)
	eta[:,1] = 0
	ind = 0	

	#NOTE: wavelength in transmission data given in um rest of files in nm and transmission given
	#in decimal whereas rest of files are given in percentage
	#Multiply efficiencies to get the combined filter, dichroic, lenses, mirrors (DCT primary and secondary), and atmosphere efficiency
	for value in filt:
		
		if value[0] < l_fil[0,0]: continue
		if value[0] > l_fil[-1,0]: continue
			
		atm = t[np.where(value[0] == t[:,0])][0]       # Atmospheric transmission
		ref = m[np.where(value[0] == m[:,0])][0]       # Reflection from primary (assume similar for secondary)
		
		len_fil = l_fil[np.where(value[0] == l_fil[:,0])][0]  # Lenses transmission
		
		#print len_fil
		
		#Lenses
		len_fil_eta = (len_fil[1]/100.)**8
		
		feta = (value[1]/100.)*(atm[1]/100.)*(ref[1]/100.)**2*len_fil_eta
		eta[ind,1] = feta
		ind = ind + 1
		
	#feta (final efficiency) includes atmospheric transmission, DCT primary 
	#(and secondary) reflectivity, filter transmission, dichroic transmission, 
	#lens collimation transmission, optical arm lens transmission
	
	#Now include folding mirror at 95%, 2 windowed surfaces at 95%, mirrored slits and 
	#mirror below collimator at 95% and detector at 70%
	eta[:,1] = eta[:,1]*0.95*0.95**2*0.95**2*.7
	
	plt.plot(eta[:,0], eta[:,1])
	#plt.show()
	#plt.clf()
	
	# Either return the full array or just the median of the filter
	if med == True:
	    keep = np.where(eta[:,1] != 0)
	    return np.median(eta[keep,1])
	else:
	    return eta

def highresspeceff(band, med=True):

    """
    NAME:
        highresspeceff
    PURPOSE:
        Finds effective efficiency of instrument from DCT primary and secondary mirror,
        fold mirror/dichroic, window, high-res spectroscopic elements, slit efficiency, 
        dichroic, lenses, and detector
    INPUTS:
        band - photometric band (YJHK)
        med  - returns median of band, otherwise returns full array from 900-2400nm
    """
    
    # Actual RIMAS filter data (nm and transmission perc.)
    filtdict = {'Y': 'highres_spec/yj_highres_trans.txt',
                'J': 'highres_spec/yj_highres_trans.txt',
                'H': 'highres_spec/hk_highres_trans.txt',
                'K': 'highres_spec/hk_highres_trans.txt'}
                
    # Lens transmission curves (NEED REPLACEMENT LENSES DATA) (nm and transmission perc.)
    lensdict = {'Y': 'Lenses/yj_T_data_interp.txt',
                'J': 'Lenses/yj_T_data_interp.txt',
                'H': 'Lenses/hk_T_data_interp.txt',
                'K': 'Lenses/hk_T_data_interp.txt'}
                
    # Load filter file for given band 
    filt = np.loadtxt(filtdict[band.upper()], skiprows=1)
    
    #Initialize efficiency to 0 and keeping only the wavelengths from filter with over 50% transmission
    eta = np.copy(filt)
    ind = 0
    
    # Actual DCT mirror data (nm and transmission perc.)
    m = np.loadtxt('DCTmirrors/mirrors_900_2400.txt')
    
    # Dichroic from RIMAS theoretical model (nm, transmission perc., reflection perc.)
    d = np.loadtxt('dichroic/dich_900_2400.txt')
    
    # Kitt Peak atmospheric transmission data (nm and transmission perc.)
    t = np.loadtxt('atmosphere/atmos_900_2400.txt')
    
    l_fil = np.loadtxt(lensdict[band.upper()])
    l_col = np.loadtxt('Lenses/collimator_T_data_interp.txt')
    
    # Finds each wavelength in filter data and finds transmission of that wavelength
    # in the atmospheric transmission data and dichroic transmission.  If photometric
    # band is Y or J the light is reflected on dichroic so use reflection.
    #Multiply efficiencies to get the combined filter, dichroic, lenses, mirrors (DCT primary and secondary), and atmosphere efficiency
    for value in filt:
    
        dich = d[np.where(value[0] == d[:,0])][0]
        
        if band.upper() in ('Y', 'J'):
            dich[1] = dich[2] # For YJ use reflection, otherwise use transmission data
            
        atm = t[np.where(value[0] == t[:,0])][0]
        ref = m[np.where(value[0] == m[:,0])][0]
        
        len_fil = l_fil[np.where(value[0] == l_fil[:,0])][0]
        len_col = l_col[np.where(value[0] == l_col[:,0])][0]
        
        #Lenses
        len_fil_eta = (1e-2**5*len_fil[1]*len_fil[2]*len_fil[3]*len_fil[4]*len_fil[5])
        len_col_eta = (1e-2**5*len_col[1]*len_col[2]*len_col[3]*len_col[4]*len_col[5])
        
        feta = (dich[1]/100.)*(value[1]/100.)*(atm[1]/100.)*(ref[1]/100.)**2*len_fil_eta*len_col_eta
        
        eta[ind,1] = feta
        ind = ind + 1
        
    #feta (filter efficiency) includes atmospheric transmission, DCT primary 
    #(and secondary) reflectivity, filter transmission, dichroic transmission,
    #lens collimation transmission, optical arm lens transmission
    
    # Stuff that isn't measured precisely
    #Now include folding mirror at 95%, 2 windowed surfaces at 95%, and detector at 80%
    eta[:,1] = eta[:,1]*0.95*0.95**2*.8
    
    # Include 50% slit efficiency
    eta[:,1] = eta[:,1]*0.5
    
    # Include mirror after dichroic beamsplitter for redder path assuming 95% reflectance
    if band.upper() in ('H', 'K'):
        eta[:,1] = eta[:,1]*0.95
    
    plt.plot(eta[:,0], eta[:,1])
    #plt.show()
    #plt.clf()
    
    # Either return the full array or just the median of the filter
    if med == True:
        keep = np.where(eta[:,1] != 0)
        return np.median(eta[keep,1])
    else:
        return eta

def lowresspeceff(band, med=True):

    """
    NAME:
        lowresspeceff
    PURPOSE:
        Finds effective efficiency of instrument from DCT primary and secondary mirror,
        fold mirror/dichroic, window, low-res spectroscopic elements, slit efficiency, 
        dichroic, lenses, and detector
    INPUTS:
        band - photometric band (YJHK)
        med  - returns median of band, otherwise returns full array from 900-2400nm
    """
    
    # Actual RIMAS filter data (nm and transmission perc.)
    filtdict = {'Y': 'lowres_spec/yj_lowres_trans.txt',
                'J': 'lowres_spec/yj_lowres_trans.txt',
                'H': 'lowres_spec/hk_lowres_trans.txt',
                'K': 'lowres_spec/hk_lowres_trans.txt'}
                
    # Lens transmission curves (NEED REPLACEMENT LENSES DATA) (nm and transmission perc.)
    lensdict = {'Y': 'Lenses/yj_T_data_interp.txt',
                'J': 'Lenses/yj_T_data_interp.txt',
                'H': 'Lenses/hk_T_data_interp.txt',
                'K': 'Lenses/hk_T_data_interp.txt'}
                
    # Load filter file for given band 
    filt = np.loadtxt(filtdict[band.upper()], skiprows=1)
    
    #Initialize efficiency to 0 and keeping only the wavelengths from filter with over 50% transmission
    eta = np.copy(filt)
    ind = 0
    
    # Actual DCT mirror data (nm and transmission perc.)
    m = np.loadtxt('DCTmirrors/mirrors_900_2400.txt')
    
    # Dichroic from RIMAS theoretical model (nm, transmission perc., reflection perc.)
    d = np.loadtxt('dichroic/dich_900_2400.txt')
    
    # Kitt Peak atmospheric transmission data (nm and transmission perc.)
    t = np.loadtxt('atmosphere/atmos_900_2400.txt')
    
    l_fil = np.loadtxt(lensdict[band.upper()])
    l_col = np.loadtxt('Lenses/collimator_T_data_interp.txt')
    
    # Finds each wavelength in filter data and finds transmission of that wavelength
    # in the atmospheric transmission data and dichroic transmission.  If photometric
    # band is Y or J the light is reflected on dichroic so use reflection.
    #Multiply efficiencies to get the combined filter, dichroic, lenses, mirrors (DCT primary and secondary), and atmosphere efficiency
    for value in filt:
    
        dich = d[np.where(value[0] == d[:,0])][0]
        
        if band.upper() in ('Y', 'J'):
            dich[1] = dich[2] # For YJ use reflection, otherwise use transmission data
            
        atm = t[np.where(value[0] == t[:,0])][0]
        ref = m[np.where(value[0] == m[:,0])][0]
        
        len_fil = l_fil[np.where(value[0] == l_fil[:,0])][0]
        len_col = l_col[np.where(value[0] == l_col[:,0])][0]
        
        #Lenses
        len_fil_eta = (1e-2**5*len_fil[1]*len_fil[2]*len_fil[3]*len_fil[4]*len_fil[5])
        len_col_eta = (1e-2**5*len_col[1]*len_col[2]*len_col[3]*len_col[4]*len_col[5])
        
        feta = (dich[1]/100.)*(value[1]/100.)*(atm[1]/100.)*(ref[1]/100.)**2*len_fil_eta*len_col_eta
        
        eta[ind,1] = feta
        ind = ind + 1
        
    #feta (filter efficiency) includes atmospheric transmission, DCT primary 
    #(and secondary) reflectivity, filter transmission, dichroic transmission,
    #lens collimation transmission, optical arm lens transmission
    
    # Stuff that isn't measured precisely
    #Now include folding mirror at 95%, 2 windowed surfaces at 95%, and detector at 80%
    eta[:,1] = eta[:,1]*0.95*0.95**2*.8
    
    # Include 50% slit efficiency
    eta[:,1] = eta[:,1]*0.5
    
    # Include mirror after dichroic beamsplitter for redder path assuming 95% reflectance
    if band.upper() in ('H', 'K'):
        eta[:,1] = eta[:,1]*0.95
    
    plt.plot(eta[:,0], eta[:,1])
    #plt.show()
    #plt.clf()
    
    # Either return the full array or just the median of the filter
    if med == True:
        keep = np.where(eta[:,1] != 0)
        return np.median(eta[keep,1])
    else:
        return eta

def quickplot(mode=''):

    if mode == 'photo':
        filtereff('Y')
        filtereff('J')
        filtereff('H')
        filtereff('K')
    if mode == 'lowres':
        lowresspeceff('Y')
        lowresspeceff('H')
    if mode == 'highres':
        highresspeceff('Y')
        highresspeceff('H')
    plt.ylabel('Efficiency')
    plt.xlabel('Wavelength (nm)')
    ax = plt.gca()
    plt.savefig('RIMAS_'+mode+'_efficiency.png', dpi=600,bbox_inches='tight', facecolor='None')


    #plt.show()
    plt.clf()
