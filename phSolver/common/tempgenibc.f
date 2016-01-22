        subroutine tempgenibc (x,  iBC,  BC, umesh )
c
c----------------------------------------------------------------------
c This routine temporally generates iBC and BC for the mesh-elastic 
c solve.
c 
c----------------------------------------------------------------------
c
        use pointer_data
        include "common.h"
c
        dimension iBC(nshg), x(numnp, nsd),
     &            BC(nshg, nelas), 
     &            umesh(numnp, nsd)
c
        integer   elascase, bit14, bit15, bit16
        real*8    rsq, rcr
c
c.... set parameters
c
      elascase = 2
      bit14 = 16384
      bit15 = 32768
      bit16 = 65536
c
c.... Testcase 1 droplet------------------------------------------
c
      if (elascase .eq. 1) then
c         tol =  0.001d0
c 
         do i = 1, numnp
c
c.... update iBC and BC of the big sphere
c
            rsq = x(i,1) * x(i,1) + x(i,2) * x(i,2) + x(i,3) * x(i,3) 
            if ( rsq .gt. 0.999  ) then
               iBC(i) = iBC(i) + bit14
               iBC(i) = iBC(i) + bit15
               iBC(i) = iBC(i) + bit16
               BC(i,1) = zero
               BC(i,2) = zero
               BC(i,3) = zero
            endif    
c
c.... update iBC and BC of the small sphere
c
c            if ( btest(iBC(i), 10) )   then
            if ( ibits(iBC(i), 3, 3) .eq. 7 )   then
               iBC(i) = iBC(i) + bit14
               iBC(i) = iBC(i) + bit15
               iBC(i) = iBC(i) + bit16
               BC(i,1) = - 0.01 * x(i,1)
               BC(i,2) = - 0.01 * x(i,2)
               BC(i,3) = - 0.01 * x(i,3)
c               write (*,'(f, f, f)') x(i,1), x(i,2), x(i,3) 
            endif 
         enddo
c
      endif
c
c.... End Testcase 1---------------------------------------------
c
c
c.... Testcase 2 droplet------------------------------------------
c
c.... tol =  0.001d0
c
      if (elascase .eq. 2) then
          rcr = x(290, 1) * x(290, 1)
     &        + x(290, 2) * x(290, 2)
     &        + x(290, 3) * x(290, 3)
c
c.... loop all nodes
c  
         do i = 1, numnp
c
c.... update iBC and BC of the big sphere
c
            rsq = x(i,1) * x(i,1) + x(i,2) * x(i,2) + x(i,3) * x(i,3) 
            if ( rsq .gt. 0.999  ) then
               iBC(i) = iBC(i) + bit14
               iBC(i) = iBC(i) + bit15
               iBC(i) = iBC(i) + bit16
               BC(i,1) = zero
               BC(i,2) = zero
               BC(i,3) = zero
            endif    
c
c.... update iBC and BC of the small sphere
c
            if ( rsq .gt. (rcr - 0.0001)
     &     .and. rsq .lt. (rcr + 0.0001) )   then
               iBC(i) = iBC(i) + bit14
               iBC(i) = iBC(i) + bit15
               iBC(i) = iBC(i) + bit16
               BC(i,1) = umesh(i,1) * Delt(1)
               BC(i,2) = umesh(i,2) * Delt(1)
               BC(i,3) = umesh(i,3) * Delt(1)
c               write (*,'(f, f, f)') x(i,1), x(i,2), x(i,3) 
            endif 
         enddo
c
      endif
c
c.... End Testcase 2---------------------------------------------
c

c.... return
c
      return
      end
