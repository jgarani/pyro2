import math
import numpy as np
import mesh.patch as patch
import sys
import compressible.problems.vortex as vortex
import util.runparams as runparams
from util import msg

usage = """
      compare the outpue n file from the vortex problem to the analytic solution.

      usage: ./vortex_error.py file
"""

if not len(sys.argv) == 2:
    print usage
    sys.exit(2)

try: file 1 = sys.argv[1]
except:
    print usage
    sys.exit(2)

myg, myd = patch.read(file1)

U_analytic = myg.scratch_array()
U_numerical = myg.scratch_array()

#compute total velocity from myd
u = myd.get_var("x-velocity")
v = myd.get_var("y-velocity")
U_numerical.d[:,:] = np.sqrt(u.d*u.d + v.d*v.d)

#use myg to create the analytic total velocity
x_c = 0.5*(myg.xmin + myg.xmax)
y_c = 0.5*(myg.ymin + myg.ymax) #is this supposed to be ymax? or just max?

j = myg.jlo
while j<= myg.jhi:
    i = myg.ilo
    while i <= myg.ihi:

        r = np.sqrt((myg.x[i] - x_c)**2 + (myg.y[i] - y_c)**2
                    q_r = (0.4*np.pi)/t_r
                    if r < 0.2:
                        u_phi = 5*r
                    elif r < 0.4:
                        u_phi = 2.0 - 5.0*r
                    else:
                        u_phi = 0
                    
                    u = -q_r*u_phi*((myg.y[j] - y_c)/r)
                    v = q_r*u_phi*((myg.x[i] - x_c)/r)
                    
                    U_analytic.d[i,j] = np.sqrt(u.d*u.d + v.d*v.d)

#compute the difference
e = U_numerical.d[:,:] - U_analytic.d[i,j]

#compute the norm
print "error norms (absolute, relative): ", e.norm(), e.norm()/dens_analytic.norm()
