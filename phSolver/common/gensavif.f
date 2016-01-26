      subroutine gensavif 
     &(
     &           ienif0, ienif1,
     &           ientmp, mattmp,
     &           mater
     &)
c
        use mattype_m
c
        include "common.h"
c
        integer, dimension(npro,nshl0), intent(out) :: ienif0
        integer, dimension(npro,nshl1), intent(out) :: ienif1
c
        integer :: i
c
        do i = 1, nshl0
          ienif0(:,i) = ienif0tmp(:,i)
        enddo
c
        do i = 1, nshl1
          ienif1(:,i) = ienif1tmp(:,i)
        enddo
c
      end subroutine gensavif
