#!/usr/bin/env python

import math
import numpy as np
import mesh.patch as patch
import sys
import compressible.problems.vortex as vortex
import util.runparams as runparams
from util import msg

usage = """
      compare the output in file from the vortex problem to
      the analytic solution.

      usage: ./vortex_error.py file1 file2
"""

if not len(sys.argv) == 3:
    print usage
    sys.exit(2)


try: file1 = sys.argv[1]
except:
    print usage
    sys.exit(2)
try: file2 = sys.argv[2]
except:
    print usage
    sys.exit(2)

myg, myd = patch.read(file1)
myg2, myd2 = patch.read(file2)

U_analytic = myg2.scratch_array()


U_numerical = myg.scratch_array()

# compute total velocity from myd
px = myd.get_var("x-momentum")
py = myd.get_var("y-momentum")
density = myd.get_var("density")
u = px/density
v = py/density

px2 = myd2.get_var("x-momentum")
py2 = myd2.get_var("y-momentum")
u2 = px2/density
v2 = py2/density

U_numerical.d[:,:] = np.sqrt(u.d*u.d + v.d*v.d)


U_analytic.d[:,:] = np.sqrt(u2.d*u2.d + v2.d*v2.d)

# use myg to create the analytic total velocity
x_c = 0.5*(myg.xmin + myg.xmax)
y_c = 0.5*(myg.ymin + myg.ymax) 

for j in range(myg.jlo, myg.jhi+1):
    for i in range(myg.ilo, myg.ihi+1):
        
        r=np.sqrt((myg.x[i]-x_c)**2 + (myg.y[j]-y_c)**2)
        t_r = 1
        q_r=(0.4*np.pi)/t_r
        if r<0.2:
            u_phi = 5*r
        elif r < 0.4:
            u_phi = 2.0-5.0*r
        else:
            u_phi = 0
        u = -q_r*u_phi*((myg.y[j]-y_c)/r)
        v = q_r*u_phi*((myg.x[i]-x_c)/r)
        
        #U_analytic.d[i,j] = np.sqrt(u*u + v*v)





# compute the difference
e = U_numerical - U_analytic

# compute the norm
print "error norms (absolute, relative): ", e.norm(), e.norm()/U_analytic.norm()


# compare the error
#dens_numerical = myd.get_var("density")
#dens_analytic = analytic.get_var("density")

#print "mesh details"
#print myg

#e = dens_numerical - dens_analytic

#print "error norms (absolute, relative): ", e.norm(), e.norm()/dens_analytic.norm()
