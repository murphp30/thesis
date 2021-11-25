#!/usr/bin/env python

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

#timerange = a.Time('2017/09/06 18:00', '2017/09/06 18:01')
#result = Fido.search(timerange, 
#        a.Instrument.aia, 
#        a.Physobs.intensity, 
#        a.Wavelength(4500*u.angstrom) | a.Wavelength(1600*u.angstrom) | a.Wavelength(304*u.angstrom) | a.Wavelength(211*u.angstrom))
#dlist=[]
#for r in result:
#    dl = r[0]
#    dlist.append(dl)
#
#downloads = Fido.fetch(*dlist, overwrite=True)
#
#smaps = sunpy.map.Map(downloads)
smaps = sunpy.map.Map([
    '/Users/murphp30/sunpy/data/aia_lev1_4500a_2017_09_06t18_00_05_68z_image_lev1.fits',
    '/Users/murphp30/sunpy/data/aia_lev1_1600a_2017_09_06t18_00_14_13z_image_lev1.fits',
    '/Users/murphp30/sunpy/data/aia_lev1_304a_2017_09_06t18_00_05_15z_image_lev1.fits',
    '/Users/murphp30/sunpy/data/aia_lev1_211a_2017_09_06t18_00_09_63z_image_lev1.fits'
    ])
fig, axs = plt.subplots(2,2, figsize=(7,7))
for smap, ax in zip(smaps, axs.flatten()):
    if smap.wavelength.value == 4500:
        vmin = 3000
        vmax = 12000 
    elif smap.wavelength.value == 304:
        vmin = 0 
        vmax = 1000
    else:
        vmin = 0 
        vmax = 10000
    smap.plot(axes=ax, vmin=vmin, vmax=vmax)
    ax.set_axis_off()
    ax.set_title("")
axs[0,0].annotate("a)", (-1000, 1000), color='white')
axs[0,1].annotate("b)", (-1000, 1000), color='white')
axs[1,0].annotate("c)", (-1000, 1000), color='white')
axs[1,1].annotate("d)", (-1000, 1000), color='white')
axs[0,0].annotate(r"{}$\AA$".format(int(smaps[0].wavelength.value)), (650,-1000), color='white')
axs[0,1].annotate(r"{}$\AA$".format(int(smaps[1].wavelength.value)), (650,-1000), color='white')
axs[1,0].annotate(r"{}$\AA$".format(int(smaps[2].wavelength.value)), (650,-1000), color='white')
axs[1,1].annotate(r"{}$\AA$".format(int(smaps[3].wavelength.value)), (650,-1000), color='white')
#plt.tight_layout()
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig("Images/pm_aia_quad.png")
plt.show()
