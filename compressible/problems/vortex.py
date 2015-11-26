from __future__ import print_function

import math
import numpy as np

import sys
import mesh.patch as patch
from util import msg

def init_data(my_data, rp):
    """ initialize the vortex problem """

    msg.bold("initializing the vortex problem...")
    
    
    
    
    # make sure that we are passed a valid patch object
    if not isinstance(my_data, patch.CellCenterData2d):
        print("ERROR: patch invalid in vortex.py")
        print(my_data.__class__)
        sys.exit()

    # get the density, momenta, and energy as separate variables
    dens = my_data.get_var("density")
    xmom = my_data.get_var("x-momentum")
    ymom = my_data.get_var("y-momentum")
    ener = my_data.get_var("energy")

    gamma = rp.get_param("eos.gamma")

    p0 = rp.get_param("vortex.p0")
    t_r = rp.get_param("vortex.t_r")
    mach = rp.get_param("vortex.mach")


    # initialize the components, remember, that ener here is
    # rho*eint + 0.5*rho*v**2, where eint is the specific
    # internal energy (erg/g)
    xmom.d[:,:] = 0.0
    ymom.d[:,:] = 0.0
    dens.d[:,:] = 0.0

    # set the density to be stratified in the y-direction
    myg = my_data.grid

    x_c = 0.5*(myg.xmin + myg.xmax)
    y_c = 0.5*(myg.ymin + myg.ymax)

    p = myg.scratch_array()

    dens.d[:,:] = 1.0
    
    nsub = 4
    
    j = myg.jlo
    while j <= myg.jhi:
        i = myg.ilo
        while i <= myg.ihi:

            uzone = 0.0
            vzone = 0.0
            pzone = 0.0

            for ii in range(nsub):
                for jj in range(nsub):
                    
                    xsub = my_data.grid.xl[i] + (my_data.grid.dx/nsub)*(ii + 0.5)
                    ysub = my_data.grid.yl[j] + (my_data.grid.dy/nsub)*(jj + 0.5)

            #initialize vortex problem
                    r=np.sqrt((xsub-x_c)**2 + (ysub-y_c)**2)
           
           
            # phi_hat = -np.sin(phi)*x_hat+cos(phi)*y_hat    
   
                    q_r=(0.4*np.pi)/t_r       
    
                    if r< 0.2:
                        u_phi = 5*r
                    elif r < 0.4:
                        u_phi = 2.0-5.0*r
                    else:
                        u_phi = 0

                    u = -q_r*u_phi*((myg.y[j]-y_c)/r)
                    v = q_r*u_phi*((myg.x[i]-x_c)/r)
                    
                    uzone += u
                    vzone += v


            
                    p_0 = (dens.d[i,j]/(gamma*(mach)**2)-(1.0/2.0)) #u_phimax is 1
  
                    if r<0.2:
                        p_r = p_0+(25.0/2)*(r**2)
                    elif r < 0.4:
                        p_r = p_0 + (25.0/2)*(r**2) + 4*(1.0-5.0*r-math.log(0.2)+math.log(r))
                    else:
                        p_r = p_0 - 2.0 + 4.0*math.log(2.0)
            
                    pzone += p_r

            u = uzone/(nsub*nsub)
            v = vzone/(nsub*nsub)

            xmom.d[i,j] = dens.d[i,j]*u
            ymom.d[i,j] = dens.d[i,j]*v


            p.d[i,j] = pzone/(nsub*nsub)

            i += 1
        j += 1


    # set the energy (P = cs2*dens)
    ener.d[:,:] = p.d[:,:]/(gamma - 1.0) + \
        0.5*(xmom.d[:,:]**2 + ymom.d[:,:]**2)/dens.d[:,:]


def finalize():
    """ print out any information to the user at the end of the run """
    pass
