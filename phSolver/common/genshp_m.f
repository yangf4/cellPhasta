      module genshp_m
c
        contains

        subroutine genshpif
     &(            shpif0, shpif1,
     &             shgif0, shgif1
     &)
c
          include "common.h"
c
          real*8, dimension(MAXTOPIF,      maxsh, MAXQPT), intent(out) :: shpif0  ! shape function          for interface element0
          real*8, dimension(MAXTOPIF,      maxsh, MAXQPT), intent(out) :: shpif1  ! shape function          for interface element1
          real*8, dimension(MAXTOPIF, nsd, maxsh, MAXQPT), intent(out) :: shgif0  ! shape function gradient 
          real*8, dimension(MAXTOPIF, nsd, maxsh, MAXQPT), intent(out) :: shgif1  ! shape function gradient 
c
          integer :: i, iblk, id, lcsyst0, lcsyst1, nshl0, nshl1
c
c... loop over the interface blocks and 
c    generate shape functions and gradients
c
          do iblk = 1, nelblif
c
            lcsyst0 = lcblkif(3, iblk)
            lcsyst1 = lcblkif(4, iblk)
c
            id = iftpid(iblk)
c
            select case (id)
            case (1)
c... Tri
              nshl0   = lcblkif(13,iblk)
c
              do i = 1, nintif0(lcsyst0)
                call shpTet(ipord, Qptif0(id,1:3,i), shpif0(id,:,i), shgif0(id,:,:,i))
              enddo
c
              shgif0(id,:,1:nshl0,1:nintif0(lcsyst0)) = 
     &          shgif0(id,:,1:nshl0,1:nintif0(lcsyst0)) / two 
c
              nshl1   = lcblkif(14,iblk)
c
              do i = 1, nintif1(lcsyst1)
                call shpTet(ipord, Qptif1(id,1:3,i), shpif1(id,:,i), shgif1(id,:,:,i))
              enddo
c
              shgif1(id,:,1:nshl1,1:nintif1(lcsyst1)) = 
     &          shgif1(id,:,1:nshl1,1:nintif1(lcsyst1)) / two 
c
            case default
c
c.... nonexistent element
c
              call error ('genshpif  ', 'elem Cat', id)
c
            end select
c
          enddo
c
        end subroutine genshpif
c
      end module genshp_m
