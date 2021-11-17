#!/usr/bin/env python

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

from astropy.constants import m_e, k_B

def f(v, T):
    """
    1D velocity distribution function normalised 
    to N(m/2pik_BT)**1/2
    """
    return np.exp(-((m_e*v**2)/(2*k_B*T)))

def f_shift(v,v0,N,T):
    """
    1D velocity distribution function normalised 
    to ratio of densities and shifted by v0
    """
    return N*np.exp(-((m_e*((v-v0)**2))/(2*k_B*T)))

T = 1e6*u.K #coronal temperature
mb_std = np.sqrt(k_B*T/m_e).decompose() #standard deviation of Maxwell Boltzman distribution
v = np.arange(-3*mb_std.value, 10*mb_std.value)*(u.m/u.s) #numpy doesn't like units sometimes
f_bg = f(v, T)
f_beam = f_shift(v, 4*mb_std, 0.4, T)
f_bot = f_bg + f_beam
pos_v = np.argwhere(v>=0)[0,0]
pos_fbot = f_bot[pos_v:]
dfdv = np.gradient(pos_fbot)
dfdv/=np.max(dfdv)
pos_dfdv0 = pos_v + np.where(dfdv > 0.01)[0][0]
pos_dfdv1 = pos_v + np.where(dfdv > 0.01)[0][-1]
fig, ax = plt.subplots(figsize=(10,6))
#plt.axis('off')
#plt.plot(v, f_bg)
ax.plot(v, f_bot)
ax.plot((1), (0), ls="", marker=">", ms=10, color="k",
        transform=ax.get_yaxis_transform(), clip_on=False)
ax.plot((0), (1), ls="", marker="^", ms=10, color="k",
        transform=ax.get_xaxis_transform(), clip_on=False)
#ax.vlines(v[pos_dfdv0].value, 0, f_bot[pos_dfdv0].value, color='r')
#ax.vlines(v[pos_dfdv1].value, 0, f_bot[pos_dfdv1].value, color='r')
fill_region = np.arange(v[pos_dfdv0].value,v[pos_dfdv1].value)
ax.fill_between(fill_region,f_bot[pos_dfdv0:pos_dfdv1].value, alpha=0.3)
ax.set_xlabel(r"$v_x$", fontsize=14)
ax.set_ylabel(r"$f(v_x)$", fontsize=14)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_xticks([])
ax.set_yticks([])
ax.spines['left'].set_position('zero')
#ax.axis['yzero'].set_axisline_style("-|>")
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_position('zero')
#ax.axis['xzero'].set_axisline_style("-|>")
ax.spines['top'].set_visible(False)

plt.tight_layout()
plt.savefig("Images/bump_on_tail.png")
plt.show()
