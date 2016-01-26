      module e3if_func_m
c
        use e3if_defs_m
c
        implicit none
c
      contains
c
      subroutine e3if_mtrx
c
        real*8, dimension(:), pointer :: rmu0, rlm0, rlm2mu0, con0
        real*8, dimension(:), pointer :: rmu1, rlm1, rlm2mu1, con1
c
c... set Ai and Kij matrices for elements 0 and 1
c
        call set_Ai (Ai0,u0(:,1),u0(:,2),u0(:,3),um0(:,1),um0(:,2),um0(:,3),rho0,T0,h0,cp0,alfaP0,betaT0)
        call set_Ai (Ai1,u1(:,1),u1(:,2),u1(:,3),um1(:,1),um1(:,2),um1(:,3),rho1,T1,h1,cp1,alfaP1,betaT1)
c
        allocate(rmu0(npro),rlm0(npro),rlm2mu0(npro),con0(npro))
        allocate(rmu1(npro),rlm1(npro),rlm2mu1(npro),con1(npro))
c
        call getdiff(rmu0, rlm0, rlm2mu0, con0, npro, mater0)
        call getdiff(rmu1, rlm1, rlm2mu1, con1, npro, mater1)
c
        call set_Kij(Kij0,rmu0,rlm0,rlm2mu0,con0,u0(:,1),u0(:,2),u0(:,3))
        call set_Kij(Kij1,rmu1,rlm1,rlm2mu1,con1,u1(:,1),u1(:,2),u1(:,3))
c
        deallocate(rmu0,rlm0,rlm2mu0,con0)
        deallocate(rmu1,rlm1,rlm2mu1,con1)
c
      end subroutine e3if_mtrx
c
      subroutine set_Ai (Ai,u1,u2,u3,um1,um2,um3,
     &                  rho,T,h,cp,alfaP,betaT)
c
        real*8, dimension(:,:,:,:), intent(out) :: Ai
        real*8, dimension(:), intent(in) :: u1,u2,u3,um1,um2,um3
        real*8, dimension(:), intent(in) :: rho,T,h,cp,alfaP,betaT
c
        real*8, dimension(npro) :: drdp,drdT,e2p,e3p,e4p,rk
c
        drdp = rho * betaT
        drdT = -rho * alfap
        rk   = pt50 * (u1*u1 + u2*u2 + u3*u3)
        e2p  = drdp * (h + rk)  - alfap * T + one
        e3p  = rho * ( h + rk)
        e4p  = drdT * (h + rk) + rho * cp
c
        Ai(:,1,1,1) = drdp * (u1 - um1)
        Ai(:,1,1,2) = rho
c       Ai(:,1,1,3) = zero
c       Ai(:,1,1,4) = zero
        Ai(:,1,1,5) = drdT * (u1 - um1)
c
        Ai(:,1,2,1) = drdp * (u1 - um1) * u1 + 1
        Ai(:,1,2,2) = rho  * (u1 - um1) + rho * u1
c       Ai(:,1,2,3) = zero
c       Ai(:,1,2,4) = zero
        Ai(:,1,2,5) = drdT * (u1 - um1) * u1
c
        Ai(:,1,3,1) = drdp * (u1 - um1) * u2 
        Ai(:,1,3,2) = rho  * u2
        Ai(:,1,3,3) = rho  * (u1 - um1)
c       Ai(:,1,3,4) = zero
        Ai(:,1,3,5) = drdT * (u1 - um1) * u2
c
        Ai(:,1,4,1) = drdp * (u1 - um1) * u3 
        Ai(:,1,4,2) = rho  * u3
c       Ai(:,1,4,3) = zero
        Ai(:,1,4,4) = rho  * (u1 - um1)
        Ai(:,1,4,5) = drdT * (u1 - um1) * u3
c
        Ai(:,1,5,1) = (u1 - um1) * e2p
        Ai(:,1,5,2) = e3p + rho * (u1 - um1) * u1
        Ai(:,1,5,3) = rho * (u1 - um1) * u2
        Ai(:,1,5,4) = rho * (u1 - um1) * u3
        Ai(:,1,5,5) = (u1 - um1) * e4p
c
        Ai(:,2,1,1) = drdp * (u2 - um2)
c       Ai(:,2,1,2) = zero
        Ai(:,2,1,3) = rho
c       Ai(:,2,1,4) = zero
        Ai(:,2,1,5) = drdT * (u2 - um2)
c
        Ai(:,2,2,1) = drdp * u1 * (u2 - um2)
        Ai(:,2,2,2) = rho  * (u2 - um2)
        Ai(:,2,2,3) = rho  * u1
c       Ai(:,2,2,4) = zero
        Ai(:,2,2,5) = drdT * u1 * (u2 - um2)
c
        Ai(:,2,3,1) = drdp * (u2 - um2) * u2 + 1
c       Ai(:,2,3,2) = zero
        Ai(:,2,3,3) = rho  * (u2 - um2) + rho * u2
c       Ai(:,2,3,4) = zero
        Ai(:,2,3,5) = drdT * (u2 - um2) * u2
c
        Ai(:,2,4,1) = drdp * (u2 - um2) * u3 
c       Ai(:,2,4,2) = zero
        Ai(:,2,4,3) = rho  * u3
        Ai(:,2,4,4) = rho  * (u2 - um2)
        Ai(:,2,4,5) = drdT * (u2 - um2) * u3
c
        Ai(:,2,5,1) = (u2 - um2) * e2p
        Ai(:,2,5,2) = rho * u1 * (u2 - um2)
        Ai(:,2,5,3) = e3p + rho * (u2 - um2) * u2
        Ai(:,2,5,4) = rho * (u2 - um2) * u3
        Ai(:,2,5,5) = (u2 - um2) * e4p
c
        Ai(:,3,1,1) = drdp * (u3 - um3)
c       Ai(:,3,1,2) = zero
c       Ai(:,3,1,3) = zero
        Ai(:,3,1,4) = rho
        Ai(:,3,1,5) = drdT * (u3 - um3)
c
        Ai(:,3,2,1) = drdp * u1 * (u3 - um3) 
        Ai(:,3,2,2) = rho  * (u3 - um3)
c       Ai(:,3,2,3) = zero
        Ai(:,3,2,4) = rho  * u1
        Ai(:,3,2,5) = drdT * u1 * (u3 - um3)
c
        Ai(:,3,3,1) = drdp * (u3 - um3) * u2 
c       Ai(:,3,3,2) = zero
        Ai(:,3,3,3) = rho  * (u3 - um3)
        Ai(:,3,3,4) = rho  * u2
        Ai(:,3,3,5) = drdT * (u3 - um3) * u2
c
        Ai(:,3,4,1) = drdp * (u3 - um3) * u3 + 1
c       Ai(:,3,4,2) = zero
c       Ai(:,3,4,3) = zero
        Ai(:,3,4,4) = rho  * u3 + rho * (u3 - um3)
        Ai(:,3,4,5) = drdT * (u3 - um3) * u3
c
        Ai(:,3,5,1) = (u3 - um3) * e2p
        Ai(:,3,5,2) = rho * u1 * (u3 - um3)
        Ai(:,3,5,3) = rho * u2 * (u3 - um3)
        Ai(:,3,5,4) = e3p + rho * (u3 - um3) * u3
        Ai(:,3,5,5) = (u3 - um3) * e4p
c
      end subroutine set_Ai
c
      subroutine set_Kij(Kij,rmu,rlm,rlm2mu,con,u1,u2,u3)
c
        real*8, dimension(:,:,:,:,:), pointer, intent(out) :: Kij
        real*8, dimension(:), pointer, intent(in) :: rmu,rlm,rlm2mu,con,u1,u2,u3
c
        Kij = zero
c
c.... K11
c
        Kij(:,1,1,2,2) = rlm2mu
        Kij(:,1,1,3,3) = rmu
        Kij(:,1,1,4,4) = rmu
        Kij(:,1,1,5,2) = rlm2mu * u1
        Kij(:,1,1,5,3) = rmu    * u2
        Kij(:,1,1,5,4) = rmu    * u3
        Kij(:,1,1,5,5) = con
c     
c.... K12
c     
        Kij(:,1,2,2,3) = rlm
        Kij(:,1,2,3,2) = rmu
        Kij(:,1,2,5,2) = rmu    * u2
        Kij(:,1,2,5,3) = rlm    * u1
c     
c.... K13
c     
        Kij(:,1,3,2,4) = rlm
        Kij(:,1,3,4,2) = rmu
        Kij(:,1,3,5,2) = rmu    * u3
        Kij(:,1,3,5,4) = rlm    * u1
c     
c.... K21
c     
        Kij(:,2,1,2,3) = rmu
        Kij(:,2,1,3,2) = rlm
        Kij(:,2,1,5,2) = rlm    * u2
        Kij(:,2,1,5,3) = rmu    * u1
c     
c.... K22
c     
        Kij(:,2,2,2,2) = rmu
        Kij(:,2,2,3,3) = rlm2mu
        Kij(:,2,2,4,4) = rmu
        Kij(:,2,2,5,2) = rmu    * u1
        Kij(:,2,2,5,3) = rlm2mu * u2
        Kij(:,2,2,5,4) = rmu    * u3
        Kij(:,2,2,5,5) = con
c     
c.... K23
c     
        Kij(:,2,3,3,4) = rlm
        Kij(:,2,3,4,3) = rmu
        Kij(:,2,3,5,3) = rmu    * u3
        Kij(:,2,3,5,4) = rlm    * u2
c     
c.... K31
c     
        Kij(:,3,1,2,4) = rmu
        Kij(:,3,1,4,2) = rlm
        Kij(:,3,1,5,2) = rlm    * u3
        Kij(:,3,1,5,4) = rmu    * u1
c     
c.... K32
c     
        Kij(:,3,2,3,4) = rmu
        Kij(:,3,2,4,3) = rlm
        Kij(:,3,2,5,3) = rlm    * u3
        Kij(:,3,2,5,4) = rmu    * u2
c     
c.... K33
c     
        kij(:,3,3,2,2) = rmu
        kij(:,3,3,3,3) = rmu
        kij(:,3,3,4,4) = rlm2mu
        kij(:,3,3,5,2) = rmu    * u1
        kij(:,3,3,5,3) = rmu    * u2
        kij(:,3,3,5,4) = rlm2mu * u3
        kij(:,3,3,5,5) = con
c
      end subroutine set_Kij
c
      end module e3if_func_m
