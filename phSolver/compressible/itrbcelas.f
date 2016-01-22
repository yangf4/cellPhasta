        subroutine itrBCElas (disp, iBC, BC, iper, ilwork)
c
c----------------------------------------------------------------------
c
c This program satisfies the boundary conditions on the Y-variables.
c
c input:
c  disp   (nshg,nelas)  : mesh-elastic displacement 
c  iBC    (nshg)        : Boundary Condition Code
c  BC     (nshg,3)      : boundary condition constraint parameters
c                         pass only 3 componenets to this subroutine
c
c output:
c  disp   (nshg,nelas)  : Adjusted displacement value(s) corresponding
c                          to a constraint d.o.f.
c  
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension disp(nshg,nelas),          iBC(nshg),
     &            BC(nshg,3)
c
        dimension ilwork(nlwork),           iper(nshg)

c        real*8 tmp(nshg), y1(nshg),q1(nshg)
c        dimension rn1(nshg), rmagn1(nshg)
c        real*8 limitcount(nflow)
c        integer locmax(1),locmin(1)
c
c  limiting...ugly but sometimes the only way
c 
c        limitcount=0
c        do i=1,nflow
c           if(ylimit(1,i).gt.0) then
c              locmax=maxloc(y(:,i))
c              locmin=minloc(y(:,i))
c              ymax=maxval(y(:,i))
c              ymin=minval(y(:,i))
c              write(33,34)i,ymax,ymin,1.*locmax,1.*locmin
c          call flush(33)
c              do j=1,numnp      ! only limit linear modes for now
c                 ypc=y(j,i)
c                 y(j,i)=min(ylimit(3,i),max(ylimit(2,i),y(j,i)))
c                 if(ypc.ne.y(j,i) )limitcount(i)=limitcount(i)+one
c              enddo
c           endif
c        enddo
c        if(maxval(limitcount).gt.0) then
c           write(33,34)myrank,(limitcount(j)/numnp,j=1,nflow)
c           call flush(33)
c        endif
c 34     format(i5,5(2x,e14.7))
c
c.... ----------------------->  Displacement  <------------------------
c.... displacement
c
        where (btest(iBC,14))
          disp(:,1) =  BC(:,1)
        endwhere
c
        where (btest(iBC,15))
          disp(:,2) =  BC(:,2)
        endwhere
c
        where (btest(iBC,16))
          disp(:,3) =  BC(:,3)
        endwhere
c
c.... communications
c 
        if (numpe > 1) then
           call commu (disp, ilwork, nelas, 'out')
        endif
c
c       slave has masters value, for abc we need to rotate it
c
        if(iabc==1) then        !are there any axisym bc's
           call rotabc(disp, iBC, 'out')
        endif
c
c.... return
c
        return
        end

