# This file is part of SEALEVEL - a tool to estimates future sea-level rise
# constrained by past obervations and long-term sea-level commitment
# Copyright (C) 2016 Matthias Mengel working at PIK Potsdam

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# LICENSE.txt for more details.

""" Code for Fig. S4 as in
    M. Mengel et al.
    Future sea-level rise constrained by observations and long-term commitment
    PNAS (2016)
    (C) Matthias Mengel working at Potsdam Institute for Climate Impact Research

    Note: You may not be able to fully plot this figure as it contains data that is not
          openly available. Request data from authors or comment out.
"""

import os
import glob
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
lib_path = os.path.abspath('../src')
sys.path.append(lib_path)
import get_calibration_data as gd
reload(gd)
import sealevel as sl
reload(sl)
import find_glacier_equi_coeffs as gicequi
reload(gicequi)
from scipy import optimize


plt.rcParams['xtick.major.pad'] = 6
plt.rcParams['font.size'] = 14
plt.rcParams['lines.markeredgewidth'] = 2
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['figure.figsize'] = 8, 10
plt.rcParams['figure.facecolor'] = "white"
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['pdf.fonttype'] = '42'

runname = __file__.split("/")[-1][0:-3]

plt.figure(8)
plt.clf()

cols = [cm.Accent(np.float(k) / 20) for k in np.arange(20)]

temp = np.linspace(0., 20., 100)

ax1 = plt.subplot(111)

i = 0
for dat in gicequi.gic_equi_marzeion12:
    lbl = "Marzeion et al. 2012" if i == 0 else ""
    ax1.plot(gicequi.gic_temp_marzeion12, dat + i / 10., lw=2, marker="|", color=cols[i], markersize=10,
             label=lbl)
    ax1.plot(temp, sl.gic_equi_functions[i](
        temp) + i / 10., lw=3, color=cols[i])
    i += 1

for dat in gicequi.gic_equi_radic:
    lbl = "Radic & Hock 2010" if i == 16 else ""
    ax1.plot(gicequi.gic_temp_radic, dat + i / 10., lw=2, marker="x", color=cols[i],
             markersize=8, label=lbl)
    ax1.plot(temp, sl.gic_equi_functions[i](
        temp) + i / 10., lw=3, color=cols[i])
    i += 1

l1 = ax1.legend(loc="lower right")

l1.draw_frame(0)
for ax in [ax1]:
    ax.grid()
    ax.set_xlim(0., 10.)
    ax.set_ylim(0., 2.3)

ax1.set_ylabel("Equilibrium sea-level contribution in m")
ax1.set_xlabel("Global temperature above preindustrial in K")

plt.draw()
plt.show()
plt.savefig("../figures/" + runname + ".pdf")
