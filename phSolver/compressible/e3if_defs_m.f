      module e3if_defs_m
c
        use number_def_m
        use pointer_data
c
        implicit none
c
        integer :: nshg, nshl0, nshl1, nenl0, nenl1, lcsyst0, lcsyst1
        integer :: npro, ndof, nsd, nflow, ipord, nqpt
c
        integer, dimension(:,:),   pointer :: sgn0, sgn1
        integer :: mater0, mater1
c
        real*8,  dimension(:,:,:), pointer :: ycl0, ycl1
        real*8,  dimension(:,:,:), pointer :: acl0, acl1
        real*8,  dimension(:,:,:), pointer :: xl0, xl1      ! nodal coordinates
        real*8,  dimension(:,:,:), pointer :: umeshl0, umeshl1 ! mesh velocity
        real*8,  dimension(:,:),   pointer :: nv0, nv1      ! interface normal vectors
        real*8,  dimension(:,:),   pointer :: vi            ! interface velocity (at integration point)
        real*8,  dimension(:),     pointer :: area
        real*8,  dimension(:,:,:), pointer :: sum_vi_area_l0, sum_vi_area_l1
        real*8,  dimension(:,:,:), pointer :: cmtrx,ctc     ! kinematic continuity matrix C
        real*8,  dimension(:,:),   pointer :: shp0, shp1    ! element shape function at quadrature point
        real*8,  dimension(:,:,:), pointer :: shgl0, shgl1  ! element shape function gradient at a quadrature point
        real*8,  dimension(:,:,:), pointer :: shg0, shg1    ! shape function gradient at a quadrature point
        real*8,  dimension(:,:,:), pointer :: dxdxi0, dxdxi1  ! element deformation tensor
        real*8,  dimension(:),     pointer :: WdetJif0, WdetJif1
c
        real*8,  dimension(:,:,:), pointer :: rl0, rl1      ! residual over the element
        real*8,  dimension(:,:),   pointer :: ri0, ri1      ! residual at the integration point
c
        real*8, dimension(:,:,:,:), pointer :: Ai0, Ai1, AiNa0, AiNa1, KijNaj0, KijNaj1, KijNajC0, KijNajC1
        real*8, dimension(:,:,:,:,:), pointer :: Kij0, Kij1
        real*8, dimension(:,:,:), pointer :: egmass00, egmass01, egmass10, egmass11
c
c... evaluated parameters at the integration point
c
        real*8, pointer :: rho0(:), u0(:,:), pres0(:), T0(:), ei0(:), um0(:,:)  ! density, velocity, pressure, temperature, on elment 0
        real*8, pointer :: rho1(:), u1(:,:), pres1(:), T1(:), ei1(:), um1(:,:)  ! density, velocity, pressure, temperature, on elment 1
c
        real*8, dimension(:), pointer :: rk0, h0, cp0, alfaP0, betaT0
        real*8, dimension(:), pointer :: rk1, h1, cp1, alfaP1, betaT1
c
          real*8, parameter :: s = 1.0d0
     &,                        e = 1.0d-1
     &,                        h = 1.d-2
     &,                        mu = 1.d0
c
c
c... properties
c
        type prop_t
          integer :: mater
          real*8 :: rmu, rlm, rlm2mu, con
          real*8, dimension(15,15) :: stiff
        end type prop_t
c
        type qpt_t           ! type def for integration point
          integer :: nshl
          real*8, pointer :: shp(:), shg(:,:), shgl(:,:)  ! shape function and gradient at quadrature point
        end type qpt_t
c
        type var_t            ! type def for integration variables
          real*8 :: rho, u(3), p, T, ei, y(5), grad_y(3,5), grad_yl(3,5)
        end type var_t
c
        type element_t
          real*8 :: dxdxi(3,3)
          real*8 :: WdetJ
        end type element_t
c
        type(qpt_t), dimension(:), pointer :: qpt0, qpt1
        type(var_t),  dimension(:), pointer :: var0,  var1
        type(prop_t), dimension(:), pointer :: prop0, prop1
        type(element_t), dimension(:), pointer :: e0, e1
c
      contains
c
        subroutine setparam_e3if
     &  (
     &    nshg_,nshl0_,nshl1_,nenl0_,nenl1_,lcsyst0_,lcsyst1_,
     &    npro_,ndof_,nsd_,nflow_,ipord_,nqpt_,
     &    egmassif00,egmassif01,egmassif10,egmassif11,
     &    materif0, materif1
     &  )
c
          integer, intent(in) :: nshg_,nshl0_,nshl1_,nenl0_,nenl1_,lcsyst0_,lcsyst1_
          integer, intent(in) :: npro_,ndof_,nsd_,nflow_,ipord_,nqpt_
          real*8, dimension(:,:,:), pointer, intent(in) :: egmassif00,egmassif01,egmassif10,egmassif11
          integer, intent(in) :: materif0, materif1
c
          nshg  = nshg_
          nshl0 = nshl0_
          nshl1 = nshl1_
          nenl0 = nenl0_
          nenl1 = nenl1_
          lcsyst0 = lcsyst0_
          lcsyst1 = lcsyst1_
          npro  = npro_
          ndof  = ndof_
          nsd   = nsd_
          nflow = nflow_
          ipord = ipord_
          nqpt  = nqpt_
c
          mater0 = materif0
          mater1 = materif1
c
          egmass00 => egmassif00
          egmass01 => egmassif01
          egmass10 => egmassif10
          egmass11 => egmassif11
c
        end subroutine setparam_e3if
c
        subroutine     malloc_e3if
c
          integer :: i
c
          allocate(sgn0(npro,nshl0))
          allocate(sgn1(npro,nshl1))
c
          allocate(ycl0(npro,nshl0,ndof))
          allocate(ycl1(npro,nshl1,ndof))
          allocate(xl0(npro,nenl0,nsd))
          allocate(xl1(npro,nenl1,nsd))
          allocate(umeshl0(npro,nenl0,nsd))
          allocate(umeshl1(npro,nenl1,nsd))
          allocate(nv0(npro,nsd))
          allocate(nv1(npro,nsd))
          allocate(area(npro))
          allocate(vi(npro,nsd))
          allocate(sum_vi_area_l0(npro,nshl0,nsd+1))
          allocate(sum_vi_area_l1(npro,nshl1,nsd+1))
          allocate(cmtrx(npro,nflow,nflow))
          allocate(ctc(npro,nflow,nflow))
          allocate(acl0(npro,nshl0,ndof))
          allocate(acl1(npro,nshl1,ndof))
c
          allocate(shp0(npro,nshl0))
          allocate(shp1(npro,nshl1))
          allocate(shgl0(npro,nsd,nshl0))
          allocate(shgl1(npro,nsd,nshl1))
          allocate(shg0(npro,nshl0,nsd))
          allocate(shg1(npro,nshl1,nsd))
          allocate(dxdxi0(npro,nsd,nsd),dxdxi1(npro,nsd,nsd))
          allocate(WdetJif0(npro),WdetJif1(npro))
c
          allocate(rl0(npro,nshl0,nflow))
          allocate(rl1(npro,nshl1,nflow))
c
          allocate(ri0(npro,nflow*(nsd+1)))
          allocate(ri1(npro,nflow*(nsd+1)))
          allocate(rho0(npro),u0(npro,nsd),pres0(npro),T0(npro),ei0(npro),um0(npro,nsd))
          allocate(rho1(npro),u1(npro,nsd),pres1(npro),T1(npro),ei1(npro),um1(npro,nsd))
          allocate(rk0(npro),h0(npro),cp0(npro),alfaP0(npro),betaT0(npro))
          allocate(rk1(npro),h1(npro),cp1(npro),alfaP1(nprO),betaT1(npro))
c
          allocate(qpt0(npro),qpt1(npro))
          allocate(var0 (npro), var1 (npro))
          allocate(prop0(npro), prop1(npro))
          allocate(e0(npro),e1(npro))
c
          do i = 1,npro
            allocate(qpt0(i)%shp(nshl0))
            allocate(qpt1(i)%shp(nshl1))
            allocate(qpt0(i)%shg(nsd,nshl0))
            allocate(qpt1(i)%shg(nsd,nshl1))
            allocate(qpt0(i)%shgl(nsd,nshl0))
            allocate(qpt1(i)%shgl(nsd,nshl1))
          enddo
c
          allocate(Ai0(npro,nsd,nflow,nflow))
          allocate(Ai1(npro,nsd,nflow,nflow))
          allocate(AiNa0(npro,nsd,nflow,nflow))
          allocate(AiNa1(npro,nsd,nflow,nflow))
          allocate(Kij0(npro,nsd,nsd,nflow,nflow))
          allocate(Kij1(npro,nsd,nsd,nflow,nflow))
          allocate(KijNaj0(npro,nsd,nflow,nflow))
          allocate(KijNaj1(npro,nsd,nflow,nflow))
          allocate(KijNajC0(npro,nsd,nflow,nflow))
          allocate(KijNajC1(npro,nsd,nflow,nflow))
c
        end subroutine malloc_e3if
c
        subroutine     free_e3if
c
          integer :: i
c
          deallocate(sgn0,sgn1)
          deallocate(ycl0,ycl1)
          deallocate(xl0,xl1)
          deallocate(umeshl0,umeshl1)
          deallocate(nv0,nv1)
          deallocate(area)
          deallocate(vi)
          deallocate(sum_vi_area_l0,sum_vi_area_l1)
          deallocate(cmtrx)
          deallocate(ctc)
          deallocate(acl0,acl1)
          deallocate(shp0,shp1)
          deallocate(shgl0,shgl1)
          deallocate(shg0,shg1)
          deallocate(dxdxi0,dxdxi1)
          deallocate(WdetJif0,WdetJif1)
          deallocate(rl0,rl1)
          deallocate(ri0,ri1)
          deallocate(rho0,u0,pres0,T0,ei0,um0)
          deallocate(rho1,u1,pres1,T1,ei1,um1)
          deallocate(rk0,h0,cp0,alfaP0,betaT0)
          deallocate(rk1,h1,cp1,alfaP1,betaT1)
c
          do i = 1,npro
            deallocate(qpt0(i)%shp,qpt1(i)%shp)
            deallocate(qpt0(i)%shg,qpt1(i)%shg)
            deallocate(qpt0(i)%shgl,qpt1(i)%shgl)
          enddo
c
          deallocate(qpt0,qpt1)
          deallocate(var0, var1)
          deallocate(prop0, prop1)
          deallocate(e0,e1)
c
          deallocate(Ai0,Ai1)
          deallocate(AiNa0,AiNa1)
          deallocate(Kij0,Kij1)
          deallocate(KijNaj0,KijNaj1)
          deallocate(KijNajC0,KijNajC1)
c
          nullify(egmass00)
          nullify(egmass01)
          nullify(egmass10)
          nullify(egmass11)
c
        end subroutine free_e3if
c
      end module e3if_defs_m
