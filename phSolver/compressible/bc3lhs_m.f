      module bc3lhs_m
c
        use pointer_data
c
        implicit none
c
      contains
c
      subroutine bc3elas_if(BC, iBC, umesh, lcblkif, nshg, ndofBC, nsd, nelblif, MAXBLK)
c
        real*8,  intent(inout) ::  BC(nshg,3)
        integer, intent(inout) :: iBC(nshg)
        real*8,  intent(in)    :: umesh(nshg,nsd)
        integer, intent(in)    :: lcblkif(14,MAXBLK+1)
        integer, intent(in) :: nshg, ndofBC, nsd, nelblif, MAXBLK
c
        integer :: iblk, iel, npro, inode0, inode1, n
        integer, pointer :: ienif0(:,:), ienif1(:,:)
 
        integer, parameter :: bit14 = 16384,
     &                        bit15 = 32768,
     &                        bit16 = 65536

c
        do iblk = 1, nelblif
c
          npro = lcblkif(1,iblk+1) - lcblkif(1,iblk)
c
          ienif0 => mienif0(iblk)%p
          ienif1 => mienif1(iblk)%p
c
          do iel = 1, npro
c
c HARDCODED FOR TETs 
c
            do n = 1,3

              inode0 = ienif0(iel,n)
              inode1 = ienif1(iel,n)
c
              BC(inode0,:) = umesh(inode0,:)
              BC(inode1,:) = umesh(inode1,:)

      write(*,'(a,i4,6f12.4)') 'iel: ',iel,umesh(inode0,:),umesh(inode1,:)

              iBC(inode0) = ibclr(iBC(inode0), 14) 
              iBC(inode0) = ibclr(iBC(inode0), 15) 
              iBC(inode0) = ibclr(iBC(inode0), 16) 

              iBC(inode1) = ibclr(iBC(inode1), 14) 
              iBC(inode1) = ibclr(iBC(inode1), 15) 
              iBC(inode1) = ibclr(iBC(inode1), 16) 

              iBC(inode0) = iBC(inode0) + bit14
              iBC(inode0) = iBC(inode0) + bit15
              iBC(inode0) = iBC(inode0) + bit16

              iBC(inode1) = iBC(inode1) + bit14
              iBC(inode1) = iBC(inode1) + bit15
              iBC(inode1) = iBC(inode1) + bit16

            enddo
c
          enddo
c
        enddo
c
      end subroutine bc3elas_if
c
      end module bc3lhs_m
