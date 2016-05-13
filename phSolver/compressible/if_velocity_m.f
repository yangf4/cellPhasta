      module if_velocity_m
c
        use mpi_def_m
        use number_def_m
        use pointer_data
c
        implicit none
c
        real*8, pointer :: sum_vi_area(:,:)    ! interface velocity weighted by interfacea area
c
      contains
c
      subroutine init_sum_vi_area(nshg,nsd)
        integer, intent(in) :: nshg,nsd
        if (associated(sum_vi_area)) 
     &    deallocate (sum_vi_area)
        allocate (sum_vi_area(nshg,nsd+1))
        sum_vi_area = zero
      end subroutine init_sum_vi_area
c
      subroutine destruct_sum_vi_area
        if (associated(sum_vi_area))
     &    deallocate(sum_vi_area)
      end subroutine destruct_sum_vi_area
c
      subroutine set_if_velocity 
     & (
     &  BC, iBC, umesh, x, ilwork,
     &  lcblkif, nshg, ndofBC, nsd, nelblif, MAXBLK, nlwork
     & )
c
        include "mpif.h"
c
        real*8,  intent(inout) ::  BC(nshg,3)
        integer, intent(inout) :: iBC(nshg)
        real*8,  dimension(nshg,nsd), intent(inout)    :: umesh
        real*8,  dimension(nshg,nsd), intent(in)    :: x
        integer, intent(in)    :: lcblkif(14,MAXBLK+1), ilwork(nlwork)
        integer, intent(in) :: nshg, ndofBC, nsd, nelblif, MAXBLK, nlwork
c
        integer :: iblk, iel, npro,inode, inode0, inode1, n, ierr
        integer, pointer :: ienif0(:,:), ienif1(:,:)

      do inode = 1,nshg
c        write(*,100) 'BEFORE: ',myrank, inode, x(inode,:), sum_vi_area(inode,:),umesh(inode,:)
      enddo
 
        if (numpe > 1) then
          call commu (sum_vi_area(:,1:3), ilwork, nsd, 'in ')
          call commu (sum_vi_area(:,4), ilwork, 1, 'in ')
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
        endif
c
        do inode = 1,nshg
c
c ... NOT SURE IF THIS IS THE BEST IF :
c
          if (sum_vi_area(inode,nsd+1) > zero) then
            umesh(inode,:) = sum_vi_area(inode,:) / sum_vi_area(inode,nsd+1)
            BC(inode,:) = umesh(inode,:)
          endif
c      write(*,100) 'AFTER: ', myrank,inode, x(inode,:), sum_vi_area(inode,:),umesh(inode,:)
        enddo
c 
100   format(a,'[',i2,'] ',i6,3f7.3,x,7e24.16)
      end subroutine set_if_velocity
c
      end module if_velocity_m
