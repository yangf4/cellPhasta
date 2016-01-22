        subroutine tempitrBC (y,ac, iBC, BC, iper, ilwork, x, umesh)
c
c----------------------------------------------------------------------
c
c This program satisfies the boundary conditions on the Y-variables.
c
c input:
c  y      (nshg,nflow)   : y variables 
c  iBC    (nshg)        : Boundary Condition Code
c  BC     (nshg,ndofBC) : boundary condition constraint parameters
c  ylimit (3,nflow)     : (1,:) limiting flag
c                         (2,:) lower bound
c                         (3,:) upper bound
c output:
c  y      (nshg,nflow)   : Adjusted V value(s) corresponding to a 
c                           constraint d.o.f.
c  x      (numnp,nsd)    : nodal coordinate. FOR ALE
c  umesh  (numnp,nsd)    : mesh velocity. FOR ALE 
c
c Farzin Shakib, Winter 1987.
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension y(nshg,nflow),             iBC(nshg),
     &            ac(nshg,nflow),            BC(nshg,ndofBC)

        dimension ilwork(nlwork),           iper(nshg)
     
        dimension x(numnp,nsd),             umesh(numnp,nsd)    !FOR ALE   
        integer   istp

        real*8 tmp(nshg), y1(nshg),q1(nshg)
        dimension rn1(nshg), rmagn1(nshg)
        real*8 limitcount(nflow)
        integer locmax(1),locmin(1)
c
c  limiting...ugly but sometimes the only way
c 
        limitcount=0
        do i=1,nflow
           if(ylimit(1,i).gt.0) then
              locmax=maxloc(y(:,i))
              locmin=minloc(y(:,i))
              ymax=maxval(y(:,i))
              ymin=minval(y(:,i))
              write(33,34)i,ymax,ymin,1.*locmax,1.*locmin
          call flush(33)
              do j=1,numnp      ! only limit linear modes for now
                 ypc=y(j,i)
                 y(j,i)=min(ylimit(3,i),max(ylimit(2,i),y(j,i)))
                 if(ypc.ne.y(j,i) )limitcount(i)=limitcount(i)+one
              enddo
           endif
        enddo
        if(maxval(limitcount).gt.0) then
           write(33,34)myrank,(limitcount(j)/numnp,j=1,nflow)
           call flush(33)
        endif
 34     format(i5,5(2x,e14.7))
c
c.... ------------------------>  Temperature  <------------------------
c.... temperature
c
        where (btest(iBC,1))
          y(:,5) =  BC(:,2)
        endwhere
c
c.... ------------------------->  Velocity  <--------------------------
c.... 3D
c
c.... x1-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 1)
            y(:,1) =  BC(:,3)  - BC(:,4) * y(:,2)
     &                         - BC(:,5) * y(:,3)
          endwhere
c
c.... x2-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 2)
            y(:,2) = BC(:,3)  - BC(:,4) * y(:,1)
     &                        - BC(:,5) * y(:,3)
          endwhere
c
c.... x1-velocity and x2-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 3)
            y(:,1) =  BC(:,3)  - BC(:,4) * y(:,3)
            y(:,2) =  BC(:,5)  - BC(:,6) * y(:,3)
          endwhere
c
c.... x3-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 4)
            y(:,3) = BC(:,3) - BC(:,4) * y(:,1)
     &                       - BC(:,5) * y(:,2)
          endwhere
c
c.... x1-velocity and x3-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 5)
            y(:,1) = BC(:,3) - BC(:,4) * y(:,2)
            y(:,3) = BC(:,5) - BC(:,6) * y(:,2)
          endwhere
c
c.... x2-velocity and x3-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 6)
            y(:,2) = BC(:,3)  - BC(:,4) * y(:,1)
            y(:,3) = BC(:,5)  - BC(:,6) * y(:,1)
          endwhere
c
c.... x1-velocity, x2-velocity and x3-velocity, 3D
c       
c.... ALE change takes place here
c          do i = 1,nshg 
c             if (ibits(iBC(i),3,3) .eq. 7) then
c               y(i,1:3) =umesh(i,1:3)
c             endif
c          enddo
c
c------> HARDCODED BC <---------
          where (ibits(iBC,3,3) .eq. 7)
            y(:,1) =  BC(:,3) !+ umesh(:,1)
            y(:,2) =  BC(:,4) !+ umesh(:,2)
            y(:,3) =  BC(:,5) !+ umesh(:,3)
          endwhere
c-------> END HARDCODED <----------
c
c       endif
c
c.... end of velocity
c
c.... -------------------------->  Density  <--------------------------
c
        if (any(btest(iBC,0))) then
c
c.... density
c
          where (btest(iBC,0))
            q1 = BC(:,1)  ! density
          elsewhere
            q1 = one
          endwhere
c
c
c
          npro = nshg
c
          ithm = 2  ! get pressure from rho and T 
c...when ithm=2 scalar is not used so tmp is in place
          call getthm (y1,        y(:,5),      tmp,
     &                 rk,         q1,         tmp,
     &                 tmp,        tmp,        tmp,
     &                 tmp,        tmp,        tmp,
     &                 tmp,        tmp)
c
          where (btest(iBC,0))
            y(:,1) = y1
          endwhere

c
        endif
c
c.... ------------------------->  Pressure  <--------------------------
c
        if (any(btest(iBC,2))) then
c
c.... pressure
c
          where (btest(iBC,2))
            y(:,4) = BC(:,1)  ! pressure here
          endwhere
c
        endif
c
c.... local periodic (and axisymmetric) boundary conditions (no communications)
c 
	do i = 1,nflow
           y(:,i) = y(iper(:),i)
           if(ires.ne.2) ac(:,i) = ac(iper(:),i)
	enddo
c
c.... communications
c 
        if (numpe > 1) then
           call commu (y, ilwork, nflow, 'out')
           if(ires.ne.2) call commu (ac, ilwork, nflow, 'out')
        endif
c
c       slave has masters value, for abc we need to rotate it
c
        if(iabc==1) then        !are there any axisym bc's
           call rotabc(y, iBC, 'out')
           if(ires.ne.2) call rotabc(ac, iBC, 'out')
        endif
     
c
c.... return
c
        return
        end


