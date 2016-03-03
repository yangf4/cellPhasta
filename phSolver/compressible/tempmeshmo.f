         subroutine tempMeshMo (x, umesh, iBC, BC)          
c
c-----------------------------------------------------------------
c
c-----------------------------------------------------------------
c
      include "common.h"
c
      dimension x(numnp,nsd), disp(numnp,nsd),
     &          umesh(numnp,nsd)
      dimension iBC(nshg),    BC(nshg,4)
      integer   casenumbeir, offset, istp
      real*8    dyn_org   ! dynamic origin 
c
c.... Update ALE mesh coordinates x
c
      casenumber = 1
      disp       = 0.0d0
      dyn_org    = lstep * 0.04
c
c.... test case 1
c
      if (casenumber .eq. 1) then
        do i = 1,numnp
          if ( (ibits(iBC(i),3,3) .eq. 7) .and.
     &         (x(i,1) .lt. 23.9) .and. (x(i,1) .gt. -15.9) .and. 
     &         ((x(i,2)*x(i,2) + x(i,3)*x(i,3)) .lt. 63.9) ) then
            if( x(i,2) .ge. 2.0 ) then ! top
              if ( x(i,1) .le. (dyn_org-4.0) ) then ! top tail
                disp(i,1) = -0.005 * (x(i,1)+8.0) + 0.04
     &                    + (x(i,1)+8.0) * (cos(pi/600) - 1.0)
     &                    - (x(i,2)-4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)-4.0) + 0.01
     &                    + (x(i,2)-4.0) * (cos(pi/600) - 1.0)
     &                    + (x(i,1)+8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)

              else if ( x(i,1) .ge. (dyn_org+4.0) ) then ! top head
                disp(i,1) = -0.005 * (x(i,1)-8.0) + 0.04
     &                    + (x(i,1)-8.0) * (cos(pi/600) - 1.0)
     &                    - (x(i,2)-4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)-4.0) + 0.01
     &                    + (x(i,2)-4.0) * (cos(pi/600) - 1.0)
     &                    + (x(i,1)-8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)
 
              else ! top middle
                disp(i,1) = -0.005 * (x(i,1)) + 0.04
     &                    + (x(i,1)) * (cos(pi/600) - 1.0)
     &                    - (x(i,2)-4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)-4.0) + 0.01
     &                    + (x(i,2)-4.0) * (cos(pi/600) - 1.0)
     &                    + (x(i,1)) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)
 
              endif ! end top--switch head middle tail

            else if( x(i,2) .le. -2.0 ) then ! bottom
              if ( x(i,1) .le. (dyn_org-4.0) ) then ! bottom tail 
                disp(i,1) = -0.005 * (x(i,1)+8.0) + 0.04
     &                    - (x(i,1)+8.0) * (cos(pi/600) - 1.0)
     &                    + (x(i,2)+4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)+4.0) + 0.01
     &                    - (x(i,2)+4.0) * (cos(pi/600) - 1.0)
     &                    - (x(i,1)+8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)
 
              else if ( x(i,1) .ge. (dyn_org+4.0) ) then ! bottom head
                disp(i,1) = -0.005 * (x(i,1)-8.0) + 0.04
     &                    - (x(i,1)-8.0) * (cos(pi/600) - 1.0)
     &                    + (x(i,2)+4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)+4.0) + 0.01
     &                    - (x(i,2)+4.0) * (cos(pi/600) - 1.0)
     &                    - (x(i,1)-8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)
 
              else ! bottom middle
                disp(i,1) = -0.005 * (x(i,1)) + 0.04
     &                    - (x(i,1)) * (cos(pi/600) - 1.0)
     &                    + (x(i,2)+4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)+4.0) + 0.01
     &                    - (x(i,2)+4.0) * (cos(pi/600) - 1.0)
     &                    - (x(i,1)) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)
 
              endif ! end bottom--switch head middle tail

            else ! middle
              if ( x(i,1) .le. (dyn_org-4.0) ) then ! middle tail
                disp(i,1) = -0.005 * (x(i,1)+8.0) + 0.04
c     &                    + (x(i,1)+8.0) * (cos(pi/600) - 1.0)
c     &                    - (x(i,2)) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)) + 0.01
c     &                    + (x(i,2)) * (cos(pi/600) - 1.0)
c     &                    + (x(i,1)+8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)
 
              else if ( x(i,1) .ge. (dyn_org+4.0) ) then ! middle head
                disp(i,1) = -0.005 * (x(i,1)-8.0) + 0.04
c     &                    + (x(i,1)-8.0) * (cos(pi/600) - 1.0)
c     &                    - (x(i,2)) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)) + 0.01
c     &                    + (x(i,2)) * (cos(pi/600) - 1.0)
c     &                    + (x(i,1)-8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)
 
              else ! middle middle
                disp(i,1) = -0.005 * (x(i,1)) + 0.04
c     &                    + (x(i,1)) * (cos(pi/600) - 1.0)
c     &                    - (x(i,2)) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)) + 0.01
c     &                    + (x(i,2)) * (cos(pi/600) - 1.0)
c     &                    + (x(i,1)) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)
 
              endif ! end middle--switch head middle tail

            endif ! end if switch top middle bottom
            BC(i,1)   = disp(i,1) / Delt(1)
            BC(i,2)   = disp(i,2) / Delt(1)
            BC(i,3)   = disp(i,3) / Delt(1)
            umesh(i,1)= BC(i,1)
            umesh(i,2)= BC(i,2)
            umesh(i,3)= BC(i,3)
          endif ! end if inside cylinder
        enddo ! end loop numnp
      endif ! end if case 1
c
c.... end test case 1
c
      return
      end

