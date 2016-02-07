
""" matthias.mengel@pik
"""

import os, glob, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
lib_path = os.path.abspath('../src')
sys.path.append(lib_path)
import get_data as gd; reload(gd)
from scipy import optimize


plt.rcParams['xtick.major.pad']  = 6
plt.rcParams['font.size']= 14
plt.rcParams['lines.markeredgewidth']=2
plt.rcParams['legend.fontsize']=14
plt.rcParams['figure.figsize'] = 8,10
plt.rcParams['figure.facecolor'] = "white"
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['pdf.fonttype'] = '42'

runname = __file__.split("/")[-1][0:-3]

plt.figure(8)
plt.clf()

cols = [cm.Accent(np.float(k)/20) for k in np.arange(20)]

temp = np.linspace(0.,20.,100)

ax1 = plt.subplot(111)

i=0
for dat in gd.gic_equi_marzeion12:
    lbl = "Marzeion et al. 2012" if i==0 else ""
    ax1.plot(gd.gic_temp_marzeion12,dat+i/10.,lw=2,marker="|",color=cols[i],markersize=10,
        label=lbl)
    ax1.plot(temp,gd.gic_equi_funcs[i](temp)+i/10.,lw=3,color=cols[i])
    i+=1

for dat in gd.gic_equi_radic:
    lbl = "Radic & Hock 2010" if i==16 else ""
    ax1.plot(gd.gic_temp_radic,dat+i/10.,lw=2,marker="x",color=cols[i],
        markersize=8, label=lbl)
    ax1.plot(temp,gd.gic_equi_funcs[i](temp)+i/10.,lw=3,color=cols[i])
    i+=1

l1 = ax1.legend(loc="lower right")
# for l in l1.get_lines():
#     l.set_color((1,1,1,1))
#     print l.get_markerfacecolor()
#     l.set_markeredgecolor((1,1,1,1))


l1.draw_frame(0)
for ax in [ax1]:
    ax.grid()
    ax.set_xlim(0.,10.)
    ax.set_ylim(0.,2.3)

ax1.set_ylabel("Equilibrium sea-level contribution in m")
ax1.set_xlabel("Global temperature above preindustrial in K")

plt.draw()
plt.show()
plt.savefig("../figures/"+runname+".pdf")
