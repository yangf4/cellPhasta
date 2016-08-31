      subroutine shpMeasure (x, ien, meshq, tetVolm)
c---------------------------------------------------------
c This subroutine calculate the shape quality 
c based on the mean ratio method. 
c
c---------------------------------------------------------
c
        include "common.h"
c
        real*8    meshq(npro)     !   ,meshV(npro)
        dimension x(numnp, nsd),        ien(npro, nshl),
     &            xl(npro, nenl, nsd)
c
        dimension tetEdge(npro, 6, nsd), 
     &            tetVolm(npro),              tetSumEdge(npro),
     &            tetEdgeLength(npro, 6),     specArea(npro),
     &            surfArea(npro, 4) 
c
c---------------------------------------------------------
c          
       call localx(x,      xl,     ien,    nsd,    'gather  ')

       if ((melasModify .ge. 1) .and. (melasModify .le. 3)) then
       tetEdge(:,1,:) = xl(:,1,:) - xl(:,4,:)
       tetEdge(:,2,:) = xl(:,2,:) - xl(:,4,:)
       tetEdge(:,3,:) = xl(:,3,:) - xl(:,4,:)
       tetEdge(:,4,:) = xl(:,1,:) - xl(:,3,:)
       tetEdge(:,5,:) = xl(:,2,:) - xl(:,3,:)
       tetEdge(:,6,:) = xl(:,1,:) - xl(:,2,:)
c
       tetVolm(:) = 1.0/6.0 * ( 
     &       tetEdge(:,1,1) * ( tetEdge(:,2,2) * tetEdge(:,3,3)
     &                        - tetEdge(:,2,3) * tetEdge(:,3,2) )
     &     - tetEdge(:,1,2) * ( tetEdge(:,2,1) * tetEdge(:,3,3)
     &                        - tetEdge(:,2,3) * tetEdge(:,3,1) )
     &     + tetEdge(:,1,3) * ( tetEdge(:,2,1) * tetEdge(:,3,2)
     &                        - tetEdge(:,2,2) * tetEdge(:,3,1) ) )
c
c... if mesaure method = mean ratio
c
       if (melasModify .eq. 1) then
         tetSumEdge(:) = tetEdge(:,1,1) * tetEdge(:,1,1)
     &                 + tetEdge(:,1,2) * tetEdge(:,1,2)
     &                 + tetEdge(:,1,3) * tetEdge(:,1,3)
     &                 + tetEdge(:,2,1) * tetEdge(:,2,1)
     &                 + tetEdge(:,2,2) * tetEdge(:,2,2)
     &                 + tetEdge(:,2,3) * tetEdge(:,2,3)
     &                 + tetEdge(:,3,1) * tetEdge(:,3,1)
     &                 + tetEdge(:,3,2) * tetEdge(:,3,2)
     &                 + tetEdge(:,3,3) * tetEdge(:,3,3)
     &                 + tetEdge(:,4,1) * tetEdge(:,4,1)
     &                 + tetEdge(:,4,2) * tetEdge(:,4,2)
     &                 + tetEdge(:,4,3) * tetEdge(:,4,3)
     &                 + tetEdge(:,5,1) * tetEdge(:,5,1)
     &                 + tetEdge(:,5,2) * tetEdge(:,5,2)
     &                 + tetEdge(:,5,3) * tetEdge(:,5,3)
     &                 + tetEdge(:,6,1) * tetEdge(:,6,1)
     &                 + tetEdge(:,6,2) * tetEdge(:,6,2)
     &                 + tetEdge(:,6,3) * tetEdge(:,6,3)
c
         meshq(:) = 15552.0 * tetVolm(:) * abs(tetVolm(:))
     &            / ( tetSumEdge(:) * tetSumEdge(:) * tetSumEdge(:) ) 
c
         do i = 1, npro
            if (meshq(i) .ge. 1.0 ) write(*,*) 'Meshq larger than one.'
            if (meshq(i) .le. 0.0 ) then 
               write(*,*) 'Meshq smaller than zero.'
            endif
         enddo
c
       endif ! end melasModify = 1, mean ratio
c
c... if melas modify method = aspect ratio beta
c...  aspect ratio = CR/(3*IR) 
c
       if (melasModify .ge. 2) then
         tetEdgeLength(:,:) = sqrt(
     &                        tetEdge(:,:,1)*tetEdge(:,:,1)   
     &                      + tetEdge(:,:,2)*tetEdge(:,:,2) 
     &                      + tetEdge(:,:,3)*tetEdge(:,:,3) )
         specArea(:) = sqrt(( tetEdgeLength(:,1)*tetEdgeLength(:,5)  
     &                      + tetEdgeLength(:,2)*tetEdgeLength(:,4)
     &                      + tetEdgeLength(:,3)*tetEdgeLength(:,6) )
     &                    * (-tetEdgeLength(:,1)*tetEdgeLength(:,5)  
     &                      + tetEdgeLength(:,2)*tetEdgeLength(:,4) 
     &                      + tetEdgeLength(:,3)*tetEdgeLength(:,6) )
     &                    * ( tetEdgeLength(:,1)*tetEdgeLength(:,5)
     &                      - tetEdgeLength(:,2)*tetEdgeLength(:,4)
     &                      + tetEdgeLength(:,3)*tetEdgeLength(:,6) )
     &                    * ( tetEdgeLength(:,1)*tetEdgeLength(:,5)
     &                      + tetEdgeLength(:,2)*tetEdgeLength(:,4)
     &                      - tetEdgeLength(:,3)*tetEdgeLength(:,6) ))
         surfArea(:,1) = 0.25 * sqrt(
     &       (tetEdgeLength(:,4)+tetEdgeLength(:,5)+tetEdgeLength(:,6))
     &     * (tetEdgeLength(:,4)+tetEdgeLength(:,5)-tetEdgeLength(:,6))
     &     * (tetEdgeLength(:,4)-tetEdgeLength(:,5)+tetEdgeLength(:,6))
     &     *(-tetEdgeLength(:,4)+tetEdgeLength(:,5)+tetEdgeLength(:,6)))
         surfArea(:,2) = 0.25 * sqrt(
     &       (tetEdgeLength(:,1)+tetEdgeLength(:,2)+tetEdgeLength(:,6))
     &     * (tetEdgeLength(:,1)+tetEdgeLength(:,2)-tetEdgeLength(:,6))
     &     * (tetEdgeLength(:,1)-tetEdgeLength(:,2)+tetEdgeLength(:,6))
     &     *(-tetEdgeLength(:,1)+tetEdgeLength(:,2)+tetEdgeLength(:,6)))
         surfArea(:,3) = 0.25 * sqrt(
     &       (tetEdgeLength(:,2)+tetEdgeLength(:,5)+tetEdgeLength(:,3))
     &     * (tetEdgeLength(:,2)+tetEdgeLength(:,5)-tetEdgeLength(:,3))
     &     * (tetEdgeLength(:,2)-tetEdgeLength(:,5)+tetEdgeLength(:,3))
     &     *(-tetEdgeLength(:,2)+tetEdgeLength(:,5)+tetEdgeLength(:,3)))
         surfArea(:,4) = 0.25 * sqrt(
     &       (tetEdgeLength(:,3)+tetEdgeLength(:,4)+tetEdgeLength(:,1))
     &     * (tetEdgeLength(:,3)+tetEdgeLength(:,4)-tetEdgeLength(:,1))
     &     * (tetEdgeLength(:,3)-tetEdgeLength(:,4)+tetEdgeLength(:,1))
     &     *(-tetEdgeLength(:,3)+tetEdgeLength(:,4)+tetEdgeLength(:,1)))
c
         meshq(:) = specArea(:) * (surfArea(:,1)+surfArea(:,2)
     &     +surfArea(:,3)+surfArea(:,4)) /216.0 /tetVolm(:) /tetVolm(:)
c
       endif ! end melasModify = 2 or 3, aspect ratio = CR/(3*IR) 
       endif ! end if melasModify = 1 or 2 or 3
c
c.... return
c
       return
       end
