# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 16:29:15 2014

@author: thomas
"""
fig, ax = plt.subplots()
ax.plot(xMax,VAvg)
ax.set_xlabel('x')
ax.set_xlim(0,np.max(xMax))
ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos:('%.01f')%(x*1e6)))
ax.set_ylabel('Potential [K]')
plt.savefig(__path__+'/Distance.png',dpi = 250)
plt.show()
