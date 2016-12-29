#!/usr/bin/env python
from astropy.io import ascii
from astropy.table import Table
import numpy as np
from scipy.interpolate import interp1d

qedata = ascii.read('QE.txt')
qe = interp1d(qedata['wl'],qedata['QE-55C'],kind='cubic')

bands = ['u','g','r','i'] #,'z','zs','y']
for band in bands :
    data = ascii.read("%sband.csv"%(band))
    w = data['Wavelength (nm)']
    mask = (w>=qedata['wl'].min() )*(w<=qedata['wl'].max())
    data[band][mask] *= qe(w[mask])
    data[band][1-mask] *= 0.
    
    wmin,wmax = int(w.min())*10,int(w.max())*10
    ww = np.linspace(wmin,wmax,wmax-wmin)
    et = interp1d(w*10,data[band],kind='linear')
    
    
    ot = Table([ww,et(ww)],dtype=('f8','f8'))
    ascii.write(ot,'%sband_eff.dat'%(band),format='basic')
 
 # EOF
 
    