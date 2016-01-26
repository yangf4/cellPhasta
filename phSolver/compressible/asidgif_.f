      subroutine asidgif_
     & (
     &   nshg_, nshl0_, nshl1_, nenl0_, nenl1_, lcsyst0_, lcsyst1_,
     &   nflow_, npro_, ndof_, nsd_, ipord_, numnp_, nqpt_,
     &   gbytes_, sbytes_, flops_,
     &   res,
     &   egmassif00, egmassif01, egmassif10, egmassif11,
     &   y,        x,       umesh,
     &   shpif0,   shpif1,  shgif0,  shgif1,
     &   qwtif0,   qwtif1,
     &   ienif0,   ienif1,
     &   materif0, materif1,
     &   time,
     &   sum_vi_area
     & )  
         use hierarchic_m
         use local_m
         use e3if_m
c
            implicit none
c
            integer, intent(in) :: nshg_, nshl0_, nshl1_, nenl0_, nenl1_,lcsyst0_,lcsyst1_
            integer, intent(in) :: nflow_, npro_, ndof_, nsd_, ipord_, numnp_, nqpt_
            integer, intent(inout) :: gbytes_, sbytes_, flops_
            real*8, dimension(nshg_,nflow_), intent(inout) :: res
            real*8, dimension(:,:,:), pointer, intent(out) :: egmassif00,egmassif01,egmassif10,egmassif11
            real*8, dimension(nshg_,ndof_),  intent(in)    :: y
            real*8, dimension(nshg_,nsd_),   intent(in)    :: x
            real*8, dimension(numnp_, nsd_), intent(inout) :: umesh
            real*8, dimension(nshl0_,nqpt_),intent(in)   :: shpif0
            real*8, dimension(nshl1_,nqpt_),intent(in)   :: shpif1
            real*8, dimension(nsd_,nshl0_,nqpt_), intent(in)  :: shgif0
            real*8, dimension(nsd_,nshl1_,nqpt_), intent(in)  :: shgif1
            real*8, dimension(nqpt_), intent(in) :: qwtif0, qwtif1
            integer, dimension(:,:), pointer, intent(in)   :: ienif0, ienif1
            integer, intent(in)   :: materif0, materif1
            real*8, intent(in) :: time
            real*8, intent(inout) :: sum_vi_area(:,:)
c
c#define debug_mode
#define droplet_case
c
      integer :: i, iel, inode
      real*8 :: r0,r1,n0(3),n1(3),t_ramp,c2
      real*8, parameter :: rmin = 0.01d0, rmax = 0.9d0,
     &                     vint = -1.0d0
c
            call setparam_e3if
     &      (
     &        nshg_,nshl0_,nshl1_,nenl0_,nenl1_,lcsyst0_,lcsyst1_,
     &        npro_,ndof_,nsd_,nflow_,ipord_,nqpt_,
     &        egmassif00,egmassif01,egmassif10,egmassif11,
     &        materif0, materif1
     &      )
c
          call malloc_e3if
c
c.... create the matrix of mode signs for the hierarchic basis 
c     functions. 
c
          if (ipord_ .gt. 1) then
           call getsgn(ienif0,sgn0,nshl0,nenl0)
           call getsgn(ienif1,sgn1,nshl1,nenl1)
        endif
c
c... localize
c
        call localy(y, ycl0, ienif0, ndof, 'gather  ', nshg, nshl0, npro, ipord, gbytes_)
        call localy(y, ycl1, ienif1, ndof, 'gather  ', nshg, nshl1, npro, ipord, gbytes_)
c
        call localx(x, xl0,  ienif0, nsd, 'gather  ', nshg, nenl0, npro, gbytes_)
        call localx(x, xl1,  ienif1, nsd, 'gather  ', nshg, nenl1, npro, gbytes_)
c
        call localx(umesh, umeshl0,  ienif0, nsd, 'gather  ', nshg, nenl0, npro, gbytes_)
        call localx(umesh, umeshl1,  ienif1, nsd, 'gather  ', nshg, nenl1, npro, gbytes_)
c
c
#ifdef debug_mode
c
c      write(*,*) 'materif0,1: ',materif0,materif1
c
      write(*,*) 'coordinates of nodes on gas side:'
      do iel = 1, 1
        do inode = 1, nshl0
         write(*,10) ienif0(iel,inode), x(ienif0(iel,inode),:)
       enddo
      enddo
c
      write(*,*) 'coordinates of nodes on liquid side:'
      do iel = 1, 1
        do inode = 1, nshl1
          write(*,10) ienif1(iel,inode), x(ienif1(iel,inode),:)
        enddo
      enddo
c
      write(*,'(a)') 'RHS residuals before e3if:'
      write(*,20) ienif0(1,1),res(ienif0(1,1),1:5)
      write(*,20) ienif0(1,4),res(ienif0(1,4),1:5)
      write(*,20) ienif1(1,1),res(ienif1(1,1),1:5)
      write(*,20) ienif1(1,4),res(ienif1(1,4),1:5)
c
      write(*,*) 'y variables at ienif0(1,:) :'
      write(*,20) ienif0(1,1), y(ienif0(1,1),:)
      write(*,20) ienif0(1,4), y(ienif0(1,4),:)
      write(*,*) 'y variables at ienif1(1,:) :'
      write(*,20) ienif1(1,1), y(ienif1(1,1),:)
      write(*,20) ienif1(1,4), y(ienif1(1,4),:)
c
c.... get the element residuals, LHS matrix, and preconditioner
c
#endif
 
       call e3if(shpif0,shpif1,shgif0,shgif1,qwtif0,qwtif1)
c
#ifdef debug_mode
c
      do iel = 1,npro
        do inode = 1,4
          if (ienif0(iel,inode) == 7) then
c            write(*,30) 'iel, inode, ien, rl:', iel, inode, ienif0(iel,inode), rl0(iel,inode,:)
          endif
          if (ienif1(iel,inode) == 465) then
c            write(*,30) 'iel, inode, ien, rl:', iel, inode, ienif1(iel,inode), rl1(iel,inode,:)
          endif
#ifdef droplet_case
          if (ienif0(iel,inode) == 98) then
            write(*,30) 'iel, inode, ien, rl:', iel, inode, ienif0(iel,inode), rl0(iel,inode,:)
          endif
          if (ienif1(iel,inode) == 3808) then
            write(*,30) 'iel, inode, ien, rl:', iel, inode, ienif1(iel,inode), rl1(iel,inode,:)
          endif
#endif
        enddo
      enddo
c
#endif
c
c.... assemble the local residual arrays
c
        call local (res, rl0, ienif0, nflow, 'scatter ', nshg,nshl0,npro,ipord,sbytes_,flops_)
        call local (res, rl1, ienif1, nflow, 'scatter ', nshg,nshl1,npro,ipord,sbytes_,flops_)
c
        call local (sum_vi_area, sum_vi_area_l0, ienif0, nsd+1, 'scatter ', nshg, nshl0,npro,ipord,sbytes_,flops_)
        call local (sum_vi_area, sum_vi_area_l1, ienif1, nsd+1, 'scatter ', nshg, nshl1,npro,ipord,sbytes_,flops_)
c
#ifdef debug_mode
c
      write(*,'(a)') 'RHS residuals after scatter e3if:'
      write(*,20) ienif0(1,1),res(ienif0(1,1),1:5)
      write(*,20) ienif0(1,4),res(ienif0(1,4),1:5)
      write(*,20) ienif1(1,1),res(ienif1(1,1),1:5)
      write(*,20) ienif1(1,4),res(ienif1(1,4),1:5)
c
#endif
c
10    format(i4,2x,3e24.16)
20    format(i4,2x,5e24.16)
30    format(a,2i3,i4,5e24.16)
40    format(i4,x,6e24.16)
50    format(a6,2x,5e24.16)
c
        call free_e3if
c
      end subroutine asidgif_
