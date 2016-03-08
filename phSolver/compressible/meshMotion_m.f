      module mesh_motion_m
c
        implicit none
c
      contains
c
      subroutine temp_mesh_motion(x,umesh,time,dt,istp,numnp,nsd,myrank)
c
        real*8, dimension(numnp,nsd), intent(out)   :: umesh
        real*8, dimension(numnp,nsd), intent(inout) :: x
        real*8, intent(in) :: time,dt
        integer, intent(in) :: istp, numnp, nsd, myrank
c
c HARDCODED mesh motion for the Stefan test cases
c
c...define options: no_ale, couette_case
c
c#define debug_mode
c#define no_ale
c#define stefan_case
c#define couette_case
c#define droplet_case
#define channel_case
c
        integer :: i
        integer, parameter :: nstp1 = 500
        real*8 :: c1,c2,c3,t1,d
        real*8, dimension(nsd) :: delta_x, xint(3)
c
#ifdef stefan_case
      integer, parameter :: idg = 7
      real*8, parameter :: vint = -2.5d0
     &,                    xmin = 0.0d0, xmax = 0.2d0
#endif
c
#ifdef droplet_case
        integer, parameter :: nd = 2
        integer, parameter :: idg = 65
c        integer, parameter :: nd = 3
c        integer, parameter :: idg = 290
        real*8 :: r,rint,n(3)
        real*8, parameter :: rin = 0.1d0, rout = 1.0d0
     &,                      rmin = 0.1d0*rin, rmax = 0.9d0*rout
     &,                      vint = -1.0d0
#endif
c
      umesh = 0.0d0
c
#ifdef  no_ale
      return
#endif
c
c      t1  = real(nstp1,8)*dt
      t1 = 0.d0
      if (time < t1) then
        c2 = time/t1
      else
        c2 = 1.d0
      endif
c
#ifdef stefan_case
c
      xint = x(idg,:)
c
      do i = 1,numnp
c
c      write(*,11) myrank
      cycle
c
        if (x(i,1) <= xint(1)) then
          c1 = (x(i,1) - xmin)/(xint(1) - xmin)
        else
          c1 = (x(i,1) - xmax)/(xint(1) - xmax)
        endif
c
        umesh(i,1) = c1*c2*vint
        x(i,1)     = x(i,1) + umesh(i,1)*dt
c        write(*,10) i,x(i,:),umesh(i,:)
      enddo
c
      return
c
#endif
c
#ifdef droplet_case
c
      xint = x(idg,:)
      rint = sqrt(sum(xint(1:nd)*xint(1:nd)))
c
      do i = 1,numnp
c
        r = sqrt(sum(x(i,1:nd)*x(i,1:nd)))
c      write(*,10) i,r,x(i,:)
c      cycle
        if (r < rmin .or. r > rmax) then
          umesh(i,:) = 0.d0
          cycle
        endif
c
        n = 0.d0
        n(1:nd) = x(i,1:nd) / r
        if (r <= rint) then
          c1 = (r - rmin) / (rint - rmin)
        else
          c1 = (r - rmax) / (rint - rmax)
        endif
c
        umesh(i,:) = c1*vint*n
c        umesh(i,:) = c1*c2*vint*n
        x(i,1:nd) = x(i,1:nd) + umesh(i,1:nd)*dt
c
      enddo
c
#endif
c
#ifdef channel_case
      do i = 1,numnp
        if (x(i,1) < 0.0d0) then
          d = abs(x(i,1)/0.05d0+1.d0)
          c1 = 0.025d0
        else
          d = abs(x(i,1)/0.05d0-1.d0)
          c1 = 0.05d0
        endif
        if (d <= 1.0d0) then
          umesh(i,1) = c1*c2*(1.0d0-d)
        endif
c      write(*,11) myrank,i,x(i,:),umesh(i,:)
        x(i,1:nsd) = x(i,1:nsd) + umesh(i,1:nsd)*dt
      enddo
#endif
10    format(i1,2x,3f7.2,2x,3f10.3)
11    format('[',i4,']',i6,x,6f10.3)
20    format(i6,2x,7e24.16)
c
      end subroutine temp_mesh_motion
c
      end module
