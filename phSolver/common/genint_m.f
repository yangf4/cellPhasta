      module genint_m
c
        integer, dimension(4), parameter :: nint_tri (4) = (/1,3,6,12/)
        integer, dimension(5), parameter :: nint_quad(5) = (/1,4,9,16,25/)
c
        contains
c
        subroutine genint_if
c
          include "common.h"
c
          integer :: iblk, id
c
c... loop over interface blocks and generate 
c    quadrature rules for each...
c
          do iblk = 1, nelblif
c
            id = iftpid(iblk) 
c
            select case (id)
            case (1)
c
              nshapeif = (ipord+1)*(ipord+2)/2
c
              nintif0(id) = nint_tri(intg(3,1))
              nintif1(id) = nint_tri(intg(3,1))

              call symtri (nintif0(id), Qptif0(id,:,:), Qwtif0(id,:), nerr)
              call symtri (nintif1(id), Qptif1(id,:,:), Qwtif1(id,:), nerr)
c
c.... adjust quadrature weights to be consistent with the
c     design of tau. 
c
              Qwtif0 = two * Qwtif0
              Qwtif1 = two * Qwtif1
c
            case default
              call error ('genint_if ', 'elem Cat', id)
            end select
c
          enddo
c
        end subroutine genint_if
c
      end module genint_m
