        subroutine geniBCElas (x, iBC, BC, umesh )
c
c----------------------------------------------------------------------
c This routine generates iBC for the mesh-elastic solve.
c
c output: 
c  iBC   (nshg)        : Boundary Condition code
c
c         = 1 * iBC_1 + 2 * iBC_2 + 4 * iBC_3
c              density   temperature   pressure
c
c    if nsd = 3:
c
c        +  8 * iBC_4 +  16 * iBC_5 +  32 * iBC_6
c           x1-velocity   x2-velocity   x3-velocity
c
c        + 64 * iBC_7 + 128 * iBC_8 + 256 * iBC_9 + 512 * iBC_10
c          sclr1         sclr2        sclr3         sclr4
c
c        + 1024 * iBC_11  + 2048* iBC_12 + 4096* iBC_13 + 8192* iBC_14
c          perioidicity     spebc          axisym         deformwall
c
c        + 16384* iBC_15  + 32768*iBC_16 + 65536*iBC_17
c          x1-mesh_disp     x2-mesh_disp   x3-mesh_disp
c
c----------------------------------------------------------------------
c
        use pointer_data
        include "common.h"
c
c Arrays in the following 1 line are now dimensioned in readnblk
c        dimension iBCtmp(numpbc)
c
        dimension iBC(nshg), x(numnp, nsd),
     &            BC(nshg, nelas),
     &            umesh(numnp, nsd)
c
c.... set the iBC array for mesh-elastic displacement
c
        iBC = ibclr(iBC, 14)
        iBC = ibclr(iBC, 15)
        iBC = ibclr(iBC, 16)
        BC  = zero
c
        call tempgenibc(x,  iBC,  BC, umesh )
c          
c.... return
c
        return
        end
