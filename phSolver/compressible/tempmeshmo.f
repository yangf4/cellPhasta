         subroutine tempMeshMo (x, umesh, iBC, BC)          
c
c-----------------------------------------------------------------
c
c-----------------------------------------------------------------
c
      include "common.h"
c
      real*8    x(numnp,nsd), disp(numnp,nsd),
     &          umesh(numnp,nsd)
      dimension iBC(nshg),    BC(nshg,4)
      integer   casenumbeir, offset, istp
      real*8    dyn_org   ! dynamic origin
      real*8    loc(444)
      real*8    norm 
c
c.... Update ALE mesh coordinates x
c
      if ( iMsCsNb .le. 0 ) then 
        write(*,*) "Change Mesh Input Source to GUI or Flow-driven"
      else
        casenumber = iMsCsNb
      endif
      disp       = zero
c
c.... test case 1
c.... original mov+rot+shrk
c
      if (casenumber .eq. 1) then
        dyn_org    = lstep * 0.04
        do i = 1,numnp
          if ( (ibits(iBC(i),14,3) .eq. 7) .and.
     &         (x(i,1) .lt. 23.9) .and. (x(i,1) .gt. -15.9) .and. 
     &         ((x(i,2)*x(i,2) + x(i,3)*x(i,3)) .lt. 63.9) ) then
            if( x(i,2) .ge. 2.0 ) then ! top
              if ( x(i,1) .le. (dyn_org-4.0) ) then ! top tail
                disp(i,1) = -0.005 * (x(i,1)+8.0) + 0.02
     &                    + (x(i,1)+8.0) * (cos(pi/600) - 1.0)
     &                    - (x(i,2)-4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)-4.0) + 0.01
     &                    + (x(i,2)-4.0) * (cos(pi/600) - 1.0)
     &                    + (x(i,1)+8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)

              else if ( x(i,1) .ge. (dyn_org+4.0) ) then ! top head
                disp(i,1) = -0.005 * (x(i,1)-8.0) + 0.02
     &                    + (x(i,1)-8.0) * (cos(pi/600) - 1.0)
     &                    - (x(i,2)-4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)-4.0) + 0.01
     &                    + (x(i,2)-4.0) * (cos(pi/600) - 1.0)
     &                    + (x(i,1)-8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)
 
              else ! top middle
                disp(i,1) = -0.005 * (x(i,1)) + 0.02
     &                    + (x(i,1)) * (cos(pi/600) - 1.0)
     &                    - (x(i,2)-4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)-4.0) + 0.01
     &                    + (x(i,2)-4.0) * (cos(pi/600) - 1.0)
     &                    + (x(i,1)) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)
 
              endif ! end top--switch head middle tail

            else if( x(i,2) .le. -2.0 ) then ! bottom
              if ( x(i,1) .le. (dyn_org-4.0) ) then ! bottom tail 
                disp(i,1) = -0.005 * (x(i,1)+8.0) + 0.02
     &                    - (x(i,1)+8.0) * (cos(pi/600) - 1.0)
     &                    + (x(i,2)+4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)+4.0) - 0.01
     &                    - (x(i,2)+4.0) * (cos(pi/600) - 1.0)
     &                    - (x(i,1)+8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)
 
              else if ( x(i,1) .ge. (dyn_org+4.0) ) then ! bottom head
                disp(i,1) = -0.005 * (x(i,1)-8.0) + 0.02
     &                    - (x(i,1)-8.0) * (cos(pi/600) - 1.0)
     &                    + (x(i,2)+4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)+4.0) - 0.01
     &                    - (x(i,2)+4.0) * (cos(pi/600) - 1.0)
     &                    - (x(i,1)-8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)
 
              else ! bottom middle
                disp(i,1) = -0.005 * (x(i,1)) + 0.02
     &                    - (x(i,1)) * (cos(pi/600) - 1.0)
     &                    + (x(i,2)+4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)+4.0) - 0.01
     &                    - (x(i,2)+4.0) * (cos(pi/600) - 1.0)
     &                    - (x(i,1)) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)
 
              endif ! end bottom--switch head middle tail

            else ! middle
              if ( x(i,1) .le. (dyn_org-4.0) ) then ! middle tail
                disp(i,1) = -0.005 * (x(i,1)+8.0) + 0.02
c     &                    + (x(i,1)+8.0) * (cos(pi/600) - 1.0)
c     &                    - (x(i,2)) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)) ! + 0.01
c     &                    + (x(i,2)) * (cos(pi/600) - 1.0)
c     &                    + (x(i,1)+8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)
 
              else if ( x(i,1) .ge. (dyn_org+4.0) ) then ! middle head
                disp(i,1) = -0.005 * (x(i,1)-8.0) + 0.02
c     &                    + (x(i,1)-8.0) * (cos(pi/600) - 1.0)
c     &                    - (x(i,2)) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)) ! + 0.01
c     &                    + (x(i,2)) * (cos(pi/600) - 1.0)
c     &                    + (x(i,1)-8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)
 
              else ! middle middle
                disp(i,1) = -0.005 * (x(i,1)) + 0.02
c     &                    + (x(i,1)) * (cos(pi/600) - 1.0)
c     &                    - (x(i,2)) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)) ! + 0.01
c     &                    + (x(i,2)) * (cos(pi/600) - 1.0)
c     &                    + (x(i,1)) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)
 
              endif ! end middle--switch head middle tail

            endif ! end if switch top middle bottom
            BC(i,1)   = disp(i,1) !/ Delt(1)
            BC(i,2)   = disp(i,2) !/ Delt(1)
            BC(i,3)   = disp(i,3) !/ Delt(1)
            umesh(i,1)= BC(i,1)
            umesh(i,2)= BC(i,2)
            umesh(i,3)= BC(i,3)
          endif ! end if inside cylinder
        enddo ! end loop numnp
      endif ! end if case 1
c
c.... end test case 1
c
c
c.... test case 2
c.... ellipical rotation
c
      if (casenumber .eq. 2) then
        do i = 1,numnp
          if ( (ibits(iBC(i),14,3) .eq. 7) .and.
     &         (abs(x(i,1)) .lt. 4.9) .and.
     &         (abs(x(i,2)) .lt. 4.9) .and. 
     &         (abs(x(i,3)) .lt. 4.9) ) then
            disp(i,1) = (x(i,1)) * (cos(pi/200) - 1.0)
     &                - (x(i,3)) *  sin(pi/200)
            disp(i,3) = (x(i,3)) * (cos(pi/200) - 1.0)
     &                + (x(i,1)) *  sin(pi/200)
            BC(i,1)   = disp(i,1) / Delt(1)
            BC(i,2)   = zero
            BC(i,3)   = disp(i,3) / Delt(1)
            umesh(i,1)= BC(i,1)
            umesh(i,2)= BC(i,2)
            umesh(i,3)= BC(i,3)
          endif ! end if inside box
        enddo ! end loop numnp
      endif ! end if case 1
c
c.... end test case 2
c
c
c.... test case 3
c
      if (casenumber .eq. 3) then
        dyn_org    = lstep * 0.04
        do i = 1,numnp
          if ( (ibits(iBC(i),14,3) .eq. 7) .and.
     &         (x(i,1) .lt. 23.9) .and. (x(i,1) .gt. -15.9) .and. 
     &         ((x(i,2)*x(i,2) + x(i,3)*x(i,3)) .lt. 63.9) ) then
            if( x(i,2) .ge. 2.0 ) then ! top
              if ( x(i,1) .le. (dyn_org-4.0) ) then ! top tail
                norm = SQRT((x(i,1)+8.0)*(x(i,1)+8.0)
     &                     +(x(i,2)-4.0)*(x(i,2)-4.0)+x(i,3)*x(i,3));
                disp(i,1) = -0.005 * (x(i,1)+8.0)/norm + 0.02
     &                    + (x(i,1)+8.0) * (cos(pi/600) - 1.0)
     &                    - (x(i,2)-4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)-4.0)/norm + 0.01
     &                    + (x(i,2)-4.0) * (cos(pi/600) - 1.0)
     &                    + (x(i,1)+8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)/norm

              else if ( x(i,1) .ge. (dyn_org+4.0) ) then ! top head
                norm = SQRT((x(i,1)-8.0)*(x(i,1)-8.0)
     &                     +(x(i,2)-4.0)*(x(i,2)-4.0)+x(i,3)*x(i,3));
                disp(i,1) = -0.005 * (x(i,1)-8.0)/norm + 0.02
     &                    + (x(i,1)-8.0) * (cos(pi/600) - 1.0)
     &                    - (x(i,2)-4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)-4.0)/norm + 0.01
     &                    + (x(i,2)-4.0) * (cos(pi/600) - 1.0)
     &                    + (x(i,1)-8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)/norm
 
              else ! top middle
                norm = SQRT((x(i,1))*(x(i,1))
     &                     +(x(i,2)-4.0)*(x(i,2)-4.0)+x(i,3)*x(i,3));
                disp(i,1) = -0.005 * (x(i,1))/norm + 0.02
     &                    + (x(i,1)) * (cos(pi/600) - 1.0)
     &                    - (x(i,2)-4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)-4.0)/norm + 0.01
     &                    + (x(i,2)-4.0) * (cos(pi/600) - 1.0)
     &                    + (x(i,1)) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)/norm
 
              endif ! end top--switch head middle tail

            else if( x(i,2) .le. -2.0 ) then ! bottom
              if ( x(i,1) .le. (dyn_org-4.0) ) then ! bottom tail 
                norm = SQRT((x(i,1)+8.0)*(x(i,1)+8.0)
     &                     +(x(i,2)+4.0)*(x(i,2)+4.0)+x(i,3)*x(i,3));
                disp(i,1) = -0.005 * (x(i,1)+8.0)/norm + 0.02
     &                    - (x(i,1)+8.0) * (cos(pi/600) - 1.0)
     &                    + (x(i,2)+4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)+4.0)/norm - 0.01
     &                    - (x(i,2)+4.0) * (cos(pi/600) - 1.0)
     &                    - (x(i,1)+8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)/norm
 
              else if ( x(i,1) .ge. (dyn_org+4.0) ) then ! bottom head
                norm = SQRT((x(i,1)-8.0)*(x(i,1)-8.0)
     &                     +(x(i,2)+4.0)*(x(i,2)+4.0)+x(i,3)*x(i,3));
                disp(i,1) = -0.005 * (x(i,1)-8.0)/norm + 0.02
     &                    - (x(i,1)-8.0) * (cos(pi/600) - 1.0)
     &                    + (x(i,2)+4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)+4.0)/norm - 0.01
     &                    - (x(i,2)+4.0) * (cos(pi/600) - 1.0)
     &                    - (x(i,1)-8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)/norm
 
              else ! bottom middle
                norm = SQRT((x(i,1))*(x(i,1))
     &                     +(x(i,2)+4.0)*(x(i,2)+4.0)+x(i,3)*x(i,3));
                disp(i,1) = -0.005 * (x(i,1))/norm + 0.02
     &                    - (x(i,1)) * (cos(pi/600) - 1.0)
     &                    + (x(i,2)+4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)+4.0)/norm - 0.01
     &                    - (x(i,2)+4.0) * (cos(pi/600) - 1.0)
     &                    - (x(i,1)) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)/norm
 
              endif ! end bottom--switch head middle tail

            else ! middle
              if ( x(i,1) .le. (dyn_org-4.0) ) then ! middle tail
                norm = SQRT((x(i,1)+8.0)*(x(i,1)+8.0)
     &                     +(x(i,2))*(x(i,2))+x(i,3)*x(i,3));
                disp(i,1) = -0.005 * (x(i,1)+8.0)/norm + 0.02
c     &                    + (x(i,1)+8.0) * (cos(pi/600) - 1.0)
c     &                    - (x(i,2)) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2))/norm ! + 0.01
c     &                    + (x(i,2)) * (cos(pi/600) - 1.0)
c     &                    + (x(i,1)+8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)/norm
 
              else if ( x(i,1) .ge. (dyn_org+4.0) ) then ! middle head
                norm = SQRT((x(i,1)-8.0)*(x(i,1)-8.0)
     &                     +(x(i,2))*(x(i,2))+x(i,3)*x(i,3));
                disp(i,1) = -0.005 * (x(i,1)-8.0)/norm + 0.02
c     &                    + (x(i,1)-8.0) * (cos(pi/600) - 1.0)
c     &                    - (x(i,2)) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2))/norm ! + 0.01
c     &                    + (x(i,2)) * (cos(pi/600) - 1.0)
c     &                    + (x(i,1)-8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)/norm
 
              else ! middle middle
                norm = SQRT((x(i,1))*(x(i,1))
     &                     +(x(i,2))*(x(i,2))+x(i,3)*x(i,3));
                disp(i,1) = -0.005 * (x(i,1))/norm + 0.02
c     &                    + (x(i,1)) * (cos(pi/600) - 1.0)
c     &                    - (x(i,2)) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2))/norm ! + 0.01
c     &                    + (x(i,2)) * (cos(pi/600) - 1.0)
c     &                    + (x(i,1)) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)/norm
 
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
c.... end test case 3
c
c
c.... test case 4
c
      if (casenumber .eq. 4) then
        dyn_org    = lstep * 0.04
        do i = 1,numnp
          if ( (ibits(iBC(i),14,3) .eq. 7) .and.
     &         (x(i,1) .lt. 23.9) .and. (x(i,1) .gt. -15.9) .and. 
     &         ((x(i,2)*x(i,2) + x(i,3)*x(i,3)) .lt. 63.9) ) then
            if( x(i,2) .ge. 2.0 ) then ! top
              if ( x(i,1) .le. (dyn_org-4.0) ) then ! top tail
                norm = SQRT((x(i,1)+8.0)*(x(i,1)+8.0)
     &                     +(x(i,2)-4.0)*(x(i,2)-4.0)+x(i,3)*x(i,3));
                disp(i,1) = -0.005 * (x(i,1)+8.0)/norm + 0.02
     &                    + (x(i,1)+8.0) * (cos(pi/600) - 1.0)
     &                    - (x(i,2)-4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)-4.0)/norm + 0.01
     &                    + (x(i,2)-4.0) * (cos(pi/600) - 1.0)
     &                    + (x(i,1)+8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)/norm

              else if ( x(i,1) .ge. (dyn_org+4.0) ) then ! top head
                norm = SQRT((x(i,1)-8.0)*(x(i,1)-8.0)
     &                     +(x(i,2)-4.0)*(x(i,2)-4.0)+x(i,3)*x(i,3));
                disp(i,1) = -0.005 * (x(i,1)-8.0)/norm + 0.02
     &                    + (x(i,1)-8.0) * (cos(pi/600) - 1.0)
     &                    - (x(i,2)-4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)-4.0)/norm + 0.01
     &                    + (x(i,2)-4.0) * (cos(pi/600) - 1.0)
     &                    + (x(i,1)-8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)/norm
 
              else ! top middle
                norm = SQRT((x(i,1))*(x(i,1))
     &                     +(x(i,2)-4.0)*(x(i,2)-4.0)+x(i,3)*x(i,3));
                disp(i,1) = -0.005 * (x(i,1))/norm + 0.02
     &                    + (x(i,1)) * (cos(pi/600) - 1.0)
     &                    - (x(i,2)-4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)-4.0)/norm + 0.01
     &                    + (x(i,2)-4.0) * (cos(pi/600) - 1.0)
     &                    + (x(i,1)) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)/norm
 
              endif ! end top--switch head middle tail

            else if( x(i,2) .le. -2.0 ) then ! bottom
              if ( x(i,1) .le. (dyn_org-4.0) ) then ! bottom tail 
                norm = SQRT((x(i,1)+8.0)*(x(i,1)+8.0)
     &                     +(x(i,2)+4.0)*(x(i,2)+4.0)+x(i,3)*x(i,3));
                disp(i,1) = -0.005 * (x(i,1)+8.0)/norm + 0.02
     &                    - (x(i,1)+8.0) * (cos(pi/600) - 1.0)
     &                    + (x(i,2)+4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)+4.0)/norm - 0.01
     &                    - (x(i,2)+4.0) * (cos(pi/600) - 1.0)
     &                    - (x(i,1)+8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)/norm
 
              else if ( x(i,1) .ge. (dyn_org+4.0) ) then ! bottom head
                norm = SQRT((x(i,1)-8.0)*(x(i,1)-8.0)
     &                     +(x(i,2)+4.0)*(x(i,2)+4.0)+x(i,3)*x(i,3));
                disp(i,1) = -0.005 * (x(i,1)-8.0)/norm + 0.02
     &                    - (x(i,1)-8.0) * (cos(pi/600) - 1.0)
     &                    + (x(i,2)+4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)+4.0)/norm - 0.01
     &                    - (x(i,2)+4.0) * (cos(pi/600) - 1.0)
     &                    - (x(i,1)-8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)/norm
 
              else ! bottom middle
                norm = SQRT((x(i,1))*(x(i,1))
     &                     +(x(i,2)+4.0)*(x(i,2)+4.0)+x(i,3)*x(i,3));
                disp(i,1) = -0.005 * (x(i,1))/norm + 0.02
     &                    - (x(i,1)) * (cos(pi/600) - 1.0)
     &                    + (x(i,2)+4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)+4.0)/norm - 0.01
     &                    - (x(i,2)+4.0) * (cos(pi/600) - 1.0)
     &                    - (x(i,1)) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)/norm
 
              endif ! end bottom--switch head middle tail

            else ! middle
              if ( x(i,1) .le. (dyn_org-4.0) ) then ! middle tail
                norm = SQRT((x(i,1)+8.0)*(x(i,1)+8.0)
     &                     +(x(i,2))*(x(i,2))+x(i,3)*x(i,3));
                disp(i,1) = -0.005 * (x(i,1)+8.0)/norm + 0.02
c     &                    + (x(i,1)+8.0) * (cos(pi/600) - 1.0)
c     &                    - (x(i,2)) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2))/norm ! + 0.01
c     &                    + (x(i,2)) * (cos(pi/600) - 1.0)
c     &                    + (x(i,1)+8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)/norm
 
              else if ( x(i,1) .ge. (dyn_org+4.0) ) then ! middle head
                norm = SQRT((x(i,1)-8.0)*(x(i,1)-8.0)
     &                     +(x(i,2))*(x(i,2))+x(i,3)*x(i,3));
                disp(i,1) = -0.005 * (x(i,1)-8.0)/norm + 0.02
c     &                    + (x(i,1)-8.0) * (cos(pi/600) - 1.0)
c     &                    - (x(i,2)) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2))/norm ! + 0.01
c     &                    + (x(i,2)) * (cos(pi/600) - 1.0)
c     &                    + (x(i,1)-8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)/norm
 
              else ! middle middle
                norm = SQRT((x(i,1))*(x(i,1))
     &                     +(x(i,2))*(x(i,2))+x(i,3)*x(i,3));
                disp(i,1) = -0.005 * (x(i,1))/norm + 0.02
c     &                    + (x(i,1)) * (cos(pi/600) - 1.0)
c     &                    - (x(i,2)) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2))/norm ! + 0.01
c     &                    + (x(i,2)) * (cos(pi/600) - 1.0)
c     &                    + (x(i,1)) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)/norm
 
              endif ! end middle--switch head middle tail

            endif ! end if switch top middle bottom
            disp(i,1) = disp(i,1) * Delt(1)
            disp(i,2) = disp(i,2) * Delt(1)
            disp(i,3) = disp(i,3) * Delt(1)
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
c.... end test case 4
c
c
c.... test case 5
c.... delay 200 steps then go x direction
c
      if (casenumber .eq. 5) then
        do i = 1,numnp
          if ( (ibits(iBC(i),14,3) .eq. 7) .and. 
     &        ((x(i,1)) .gt. -119e-6) .and. ((x(i,1)) .lt. 279e-6) .and.
     &        (abs(x(i,2)) .lt. 12e-6) .and.
     &        (abs(x(i,3)) .lt. 12e-6) ) then
            if ( lstep .gt. 200 ) then
              disp(i,1) = 1.0e-6
              BC(i,1)   = disp(i,1) / Delt(1)
              BC(i,2)   = zero
              BC(i,3)   = zero
              umesh(i,1)= BC(i,1)
              umesh(i,2)= BC(i,2)
              umesh(i,3)= BC(i,3)
            else 
              BC(i,1)   = zero
              BC(i,2)   = zero
              BC(i,3)   = zero
              umesh(i,1)= BC(i,1)
              umesh(i,2)= BC(i,2)
              umesh(i,3)= BC(i,3)
            endif ! end if larger than ramping time
          endif ! end if inside channel
        enddo ! end loop numnp
      endif ! end if case 1
c
c.... end test case 5
c
c.... test case 6
c.... ARO test case size mov+rot+shrk
c
      if (casenumber .eq. 6) then
        dyn_org    = lstep * 2e-5 
        do i = 1,numnp
          if ( (ibits(iBC(i),14,3) .eq. 7) .and.
     &         (x(i,1) .lt. 32.9e-3) .and. (x(i,1) .gt. 1e-4) .and. 
     &         ((x(i,2)*x(i,2) + x(i,3)*x(i,3)) .lt. 64e-6) ) then
            if( x(i,2) .ge. 1e-3 ) then ! top
              if ( x(i,1) .le. (dyn_org+9e-3) ) then ! top tail
                disp(i,1) = -0.005 * (x(i,1)-6.0e-3) + 2e-5
     &                    + (x(i,1)-6.0e-3) * (cos(pi/600) - 1.0)
     &                    - (x(i,2)-2.4e-3) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)-2.4e-3) + 1e-5
     &                    + (x(i,2)-2.4e-3) * (cos(pi/600) - 1.0)
     &                    + (x(i,1)-6.0e-3) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)

              else if ( x(i,1) .ge. (dyn_org+13e-3) ) then ! top head
                disp(i,1) = -0.005 * (x(i,1)-14.8e-3) + 2e-5
     &                    + (x(i,1)-14.8e-3) * (cos(pi/600) - 1.0)
     &                    - (x(i,2)-2.4e-3) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)-2.4e-3) + 1e-5
     &                    + (x(i,2)-2.4e-3) * (cos(pi/600) - 1.0)
     &                    + (x(i,1)-14.8e-3) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)
 
              else ! top middle
                disp(i,1) = -0.005 * (x(i,1)-9.4e-3) + 2e-5
     &                    + (x(i,1)-9.4e-3) * (cos(pi/600) - 1.0)
     &                    - (x(i,2)-2.4e-3) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)-2.4e-3) + 1e-5
     &                    + (x(i,2)-2.4e-3) * (cos(pi/600) - 1.0)
     &                    + (x(i,1)-9.4e-3) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)
 
              endif ! end top--switch head middle tail

            else if( x(i,2) .le. -1e-3 ) then ! bottom
              if ( x(i,1) .le. (dyn_org+9e-3) ) then ! bottom tail 
                disp(i,1) = -0.005 * (x(i,1)-6.0e-3) + 2e-5
     &                    - (x(i,1)-6.0e-3) * (cos(pi/600) - 1.0)
     &                    + (x(i,2)+2.4e-3) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)+2.4e-3) - 1e-5
     &                    - (x(i,2)+2.4e-3) * (cos(pi/600) - 1.0)
     &                    - (x(i,1)-6.0e-3) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)
 
              else if ( x(i,1) .ge. (dyn_org+13e-3) ) then ! bottom head
                disp(i,1) = -0.005 * (x(i,1)-14.8e-3) + 2e-5
     &                    - (x(i,1)-14.8e-3) * (cos(pi/600) - 1.0)
     &                    + (x(i,2)+2.4e-3) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)+2.4e-3) - 1e-5
     &                    - (x(i,2)+2.4e-3) * (cos(pi/600) - 1.0)
     &                    - (x(i,1)-14.8e-3) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)
 
              else ! bottom middle
                disp(i,1) = -0.005 * (x(i,1)-9.4e-3) + 2e-5
     &                    - (x(i,1)-9.4e-3) * (cos(pi/600) - 1.0)
     &                    + (x(i,2)+2.4e-3) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)+2.4e-3) - 1e-5
     &                    - (x(i,2)+2.4e-3) * (cos(pi/600) - 1.0)
     &                    - (x(i,1)-9.4e-3) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)
 
              endif ! end bottom--switch head middle tail

            else ! middle
              if ( x(i,1) .le. (dyn_org+9e-3) ) then ! middle tail
                disp(i,1) = -0.005 * (x(i,1)-6.0e-3) + 2e-5
c     &                    + (x(i,1)-6.0e-3) * (cos(pi/600) - 1.0)
c     &                    - (x(i,2)) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)) ! + 1e-5
c     &                    + (x(i,2)) * (cos(pi/600) - 1.0)
c     &                    + (x(i,1)-6.0e-3) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)
 
              else if ( x(i,1) .ge. (dyn_org+13e-3) ) then ! middle head
                disp(i,1) = -0.005 * (x(i,1)-14.8e-3) + 2e-5
c     &                    + (x(i,1)-14.8e-3) * (cos(pi/600) - 1.0)
c     &                    - (x(i,2)) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)) ! + 1e-5
c     &                    + (x(i,2)) * (cos(pi/600) - 1.0)
c     &                    + (x(i,1)-14.8e-3) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)
 
              else ! middle middle
                disp(i,1) = -0.005 * (x(i,1)-9.4e-3) + 2e-5
c     &                    + (x(i,1)-9.4e-3) * (cos(pi/600) - 1.0)
c     &                    - (x(i,2)) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)) ! + 1e-5
c     &                    + (x(i,2)) * (cos(pi/600) - 1.0)
c     &                    + (x(i,1)-9.4e-3) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)
 
              endif ! end middle--switch head middle tail

            endif ! end if switch top middle bottom
            BC(i,1)   = disp(i,1) !/ Delt(1)
            BC(i,2)   = disp(i,2) !/ Delt(1)
            BC(i,3)   = disp(i,3) !/ Delt(1)
            umesh(i,1)= BC(i,1)
            umesh(i,2)= BC(i,2)
            umesh(i,3)= BC(i,3)
          endif ! end if inside cylinder
        enddo ! end loop numnp
      endif ! end if case 6
c
c.... end test case 6
c
c
c.... test case 7
c.... delay 200 steps then go x direction
c.... mimic the real experiment. 
c
      if (casenumber .eq. 7) then
        do i = 1,numnp
          dyn_org = (lstep-300) * 2.0 * Delt(1)
          if ( (ibits(iBC(i),14,3) .eq. 7) .and. 
     &        ((x(i,1)) .gt. -119e-6) .and. ((x(i,1)) .lt. 279e-6) .and.
     &        (abs(x(i,2)) .lt. 12e-6) .and.
     &        (abs(x(i,3)) .lt. 12e-6) ) then
            if ( (lstep .ge. 300) .and. (lstep .lt. 420) ) then
              BC(i,1)   = 2.0
              BC(i,2)   = zero
              BC(i,3)   = zero
            else if ( lstep .ge. 420 ) then
              if ( (x(i,1) - dyn_org) .ge. 0.0 ) then            
                disp(i,1) = -0.05 * (x(i,1) - dyn_org) 
              else 
                disp(i,1) = 0.0
              endif
              disp(i,2) =  0.01274 *  x(i,2) 
              disp(i,3) =  0.01274 *  x(i,3)
              BC(i,1)   = disp(i,1) / Delt(1) + 2.0
              BC(i,2)   = disp(i,2) / Delt(1)
              BC(i,3)   = disp(i,3) / Delt(1)
            else 
              BC(i,1)   = zero
              BC(i,2)   = zero
              BC(i,3)   = zero
            endif ! end if larger than ramping time
            umesh(i,1)= BC(i,1)
            umesh(i,2)= BC(i,2)
            umesh(i,3)= BC(i,3)
          endif ! end if inside channel
        enddo ! end loop numnp
      endif ! end if case 1
c
c.... end test case 7
c
c
c.... test case 8
c.... delay 200 steps then go x direction
c.... mimic the real experiment. 
c
      if (casenumber .eq. 8) then
        loc(1:300)= 0.0
        loc(301)    = 1.1
        loc(302:375) = 1.8
        do i = 376,424
          loc(i) = loc(i-1) + 0.0142857
        enddo
        loc(425)  = 2.5
        do i = 426,444
          loc(i) = loc(i-1) - 0.1052631
        enddo
        dyn_org = 0.0
        do i = 300,lstep      
          dyn_org = dyn_org + loc(i-1) * Delt(1)
        enddo
        do i = 1,numnp
          if ( (ibits(iBC(i),14,3) .eq. 7) .and. 
     &        ((x(i,1)) .gt. -119e-6) .and. ((x(i,1)) .lt. 279e-6) .and.
     &        (abs(x(i,2)) .lt. 12e-6) .and.
     &        (abs(x(i,3)) .lt. 12e-6) ) then
            if ( (lstep .ge. 300) .and. (lstep .lt. 425) ) then
              BC(i,1)   = loc(lstep)
              BC(i,2)   = zero
              BC(i,3)   = zero
            else if ( lstep .ge. 425 ) then
              if ( x(i,1) .ge. dyn_org ) then            
                disp(i,1) = -0.05 * (x(i,1) - dyn_org) 
              else 
                disp(i,1) = 0.0
              endif
              disp(i,2) =  0.01274 *  x(i,2)
              if ( x(i,3) .ge. 0.0 ) then
                disp(i,3) = 0.02548 * x(i,3)
              else
                disp(i,3) = 0.0
              endif
              BC(i,1)   = disp(i,1) / Delt(1) + loc(lstep)
              BC(i,2)   = disp(i,2) / Delt(1)
              BC(i,3)   = disp(i,3) / Delt(1)
            else 
              BC(i,1)   = zero
              BC(i,2)   = zero
              BC(i,3)   = zero
            endif ! end if larger than ramping time
            umesh(i,1)= BC(i,1)
            umesh(i,2)= BC(i,2)
            umesh(i,3)= BC(i,3)
          endif ! end if inside channel
        enddo ! end loop numnp
      endif ! end if case 8
c
c.... end test case 8
c
      return
      end
c
