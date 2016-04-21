      module e3if_m
c
        use mpi_def_m
        use hierarchic_m
        use matdat_def_m
        use e3if_defs_m
        use e3if_func_m
        use e3if_diff_m
        use eqn_state_m
        use e3metric_m
        use e3if_lhs_m
        use e3if_vi_m
c
        implicit none
c
      contains

        subroutine e3if(shpif0,shpif1,shgif0,shgif1,qwtif0,qwtif1)
c
          real*8, dimension(nshl0,nqpt), intent(in) :: shpif0
          real*8, dimension(nshl1,nqpt), intent(in) :: shpif1
          real*8, dimension(nsd,nshl0,nqpt), intent(in) :: shgif0
          real*8, dimension(nsd,nshl1,nqpt), intent(in) :: shgif1
          real*8, dimension(nqpt), intent(in) :: qwtif0, qwtif1
c
          integer :: intp
      integer :: iel
c
c      write(*,*) 'In e3if...'
c
c
          rl0 = zero
          rl1 = zero
c
          sum_vi_area_l0 = zero
          sum_vi_area_l1 = zero
c
          do intp = 1, nqpt
c
            ri0 = zero
            ri1 = zero
c
            call calc_normal_vectors(nv0,area,WdetJif0,xl0,qwtif0,lcsyst0,intp)
            call calc_normal_vectors(nv1,area,WdetJif1,xl1,qwtif1,lcsyst1,intp)
c
c... do not include this quadrature point, if Det. .eq. 0
c
            if (qwtif0(intp) .eq. zero) cycle
            if (qwtif1(intp) .eq. zero) cycle
c
c... create a matrix of shape fucntions (and derivatives) for each
c    element at this quadrature point. These arrays will contain
c    the correct signs for the higher order hierarchic basis
c
            call  getshp(shp0, shgl0, shpif0, shgif0, sgn0, npro, nsd, nshl0, nqpt, nenl0, intp, ipord)
            call  getshp(shp1, shgl1, shpif1, shgif1, sgn1, npro, nsd, nshl1, nqpt, nenl1, intp, ipord)
c      call getshp_if(shp0,shp1,shgl0,shgl1,shpif0,shpif1,shgif0,shgif1,xl0,xl1,npro,nsd,nshl0,nqpt,nenl0,intp)
c
c... calculate the integration varibles
c
            call e3if_var
c
            call e3if_mtrx
c
c... calculate the contribution of the jump in the fluxes, across the interface
c
      prop0%mater = mater0
      prop1%mater = mater1
c
c... Element Metrics
c
            call e3metric(shg0, dxdxi0,shgl0,xl0)
            call e3metric(shg1, dxdxi1,shgl1,xl1)
c
            call e3var(var0, ycl0, shp0, shgl0, shg0, dxdxi0, nshl0) 
            call e3var(var1, ycl1, shp1, shgl1, shg1, dxdxi1, nshl1) 
c
            call calc_stiff(prop0, var0, mater0)
            call calc_stiff(prop1, var1, mater1)
c
            call e3if_flux
c
c... calculate the contribution of the interface velocity V_i
c
c            call e3if_vi
c
            call calc_cmtrx
            call kinematic_condition
c
c...LHS calculations...
c
            call set_lhs_matrices
c
c            call calc_egmass(egmass00,AiNa0,KijNaj0,KijNaj0,KijNajC0,shp0,shp0,nv0,nv0,WdetJif0,nshl0,nshl0)
c            call calc_egmass(egmass01,AiNa1,KijNaj0,KijNaj1,KijNajC0,shp0,shp1,nv0,nv1,WdetJif0,nshl0,nshl1)
c            call calc_egmass(egmass10,AiNa0,KijNaj1,KijNaj0,KijNajC1,shp1,shp0,nv1,nv0,WdetJif1,nshl1,nshl0)
c            call calc_egmass(egmass11,AiNa1,KijNaj1,KijNaj1,KijNajC1,shp1,shp1,nv1,nv1,WdetJif1,nshl1,nshl1)
c
            call e3if_wmlt(rl0, ri0, shp0, shg0, WdetJif0, nshl0)
            call e3if_wmlt(rl1, ri1, shp1, shg1, WdetJif1, nshl1)
c
c ... Interface velocity calculations...
c
            if     (mat_eos(mater0,1) == ieos_ideal_gas) then
              call calc_vi(pres0,nv0)
            elseif (mat_eos(mater1,1) == ieos_ideal_gas) then
              call calc_vi(pres1,nv1)
            else
              call error ('wrong mater: ', 'calc vi', 0)
            endif
c      do iel = 1,npro
c        write(*,'(a,i4,a,i4,3f12.4)') '[',myrank,'] iel: ',iel,vi(iel,:)
c      enddo
c
            call calc_vi_area_node(sum_vi_area_l0,shp0,nshl0)
            call calc_vi_area_node(sum_vi_area_l1,shp1,nshl1)
c
          enddo  ! end of integeration points loop 
c
        end subroutine e3if
c
        subroutine e3var(var,y,shp,shgl,shg,dxdxi,nshl)
c
          type(var_t), pointer, intent(out) :: var(:)
          real*8, pointer, intent(in) :: shp(:,:),shgl(:,:,:), shg(:,:,:)
          real*8, dimension(npro,nshl,nflow), intent(in) :: y
          real*8, dimension(:,:,:), intent(in) :: dxdxi
          integer, intent(in) :: nshl
c
          integer :: iel,ivar,isd,jsd,n
          real*8 :: grad_y(npro)
c
          do iel = 1,npro
            do ivar = 1,nflow
              var(iel)%y(ivar) = sum_qpt(nshl,y(iel,:,ivar),shp(iel,:))
              do isd = 1,nsd
c                var(iel)%grad_yl(isd,ivar) = sum_qpt(nshl,y(iel,:,ivar),shgl(iel,isd,:))
                var(iel)%grad_y(isd,ivar) = sum_qpt(nshl,y(iel,:,ivar),shg(iel,:,isd))
              enddo
            enddo
          enddo
      do n = 1,nshl
c         write(*,*)'y:     ', n,y(1,n,5)
c         write(*,*)'shg:   ', n,shg(1,n,:)
c         write(*,*)'grad_y:', n,var(1)%grad_y(:,5)
      enddo
c
      return
c
c...Multiply by the element deformation tensor to get the 
c...global Y-gradient:
c
          do ivar = 1, nflow
            do isd = 1, nsd
c
              grad_y = zero
c
              do jsd = 1, nsd
                grad_y = grad_y + dxdxi(:,jsd,isd)*var(:)%grad_yl(jsd,ivar)
              enddo
c
              var(:)%grad_y(isd,ivar) = grad_y(:)
c
            enddo
          enddo
c
        end subroutine e3var
c
        subroutine     e3if_var
c
          implicit none 
c
          integer :: iel
          real*8, dimension(npro) :: tmp
c
c... compute variables at the integration point
c
          call get_var(npro,nshl0,pres0,u0,T0,ycl0,shp0,var0)
          call get_var(npro,nshl1,pres1,u1,T1,ycl1,shp1,var1)
c
          call get_mesh_velocity(um0,umeshl0,shp0,npro,nshl0)
          call get_mesh_velocity(um1,umeshl1,shp1,npro,nshl1)
c
          call getthm(rho0, ei0, pres0, T0, npro, mater0
     &,               h0, tmp, cp0, alfaP0, betaT0, tmp, tmp)
          call getthm(rho1, ei1, pres1, T1, npro, mater1
     &,               h1, tmp, cp1, alfaP1, betaT1, tmp, tmp)
c
        end subroutine e3if_var
c
        subroutine get_var(npro,nshl,p,u,T,y,shp,var)
c
          integer, intent(in) :: npro,nshl
          real*8, dimension(npro), intent(out) :: p,T
          real*8, dimension(npro,nsd), intent(out) :: u
          real*8, dimension(npro,nshl,nflow), intent(in) :: y
          real*8, dimension(npro,nshl), intent(in) :: shp
          type(var_t), dimension(:), pointer, intent(inout) :: var
c
          integer :: iel
c
          do iel = 1,npro
            p(iel)   = sum_qpt(nshl,y(iel,:,1),shp(iel,:))
            u(iel,1) = sum_qpt(nshl,y(iel,:,2),shp(iel,:))
            u(iel,2) = sum_qpt(nshl,y(iel,:,3),shp(iel,:))
            u(iel,3) = sum_qpt(nshl,y(iel,:,4),shp(iel,:))
            T(iel)   = sum_qpt(nshl,y(iel,:,5),shp(iel,:))
          enddo
c
          var%p    = p
          var%u(1) = u(:,1)
          var%u(2) = u(:,2)
          var%u(3) = u(:,3)
          var%T    = T
c
        end subroutine get_var
c
        subroutine get_mesh_velocity(um,umeshl,shp,npro,nshl)
c
          integer, intent(in) :: npro, nshl
          real*8, dimension(npro,nsd), intent(out) :: um
          real*8, dimension(npro,nshl,nsd), intent(in)  :: umeshl
          real*8, dimension(npro,nshl), intent(in) :: shp
c
          integer :: n, isd
c
           um = zero
c
           do isd = 1,nsd
             do n = 1, nshl
                um(:,isd) = um(:,isd) + shp(:,n)*umeshl(:,n,isd)
             enddo
           enddo
c
        end subroutine get_mesh_velocity
c       
        real*8 function sum_qpt(nshl,y,shp)
c
c...      y(int)=SUM_{a=1}^nshl (N_a(int) Ya)
c
          integer :: n,nshl
          real*8  :: y(nshl),shp(nshl)
c
          sum_qpt = zero
c
          do n = 1,nshl
            sum_qpt = sum_qpt + shp(n)*y(n)
          enddo
c
        end function sum_qpt
c
        subroutine calc_normal_vectors(nv,area,WdetJ,xl,qwt,lcsyst,intp)
c
          real*8, dimension(:,:), pointer, intent(out) :: nv
          real*8, dimension(:),   pointer, intent(out) :: area, WdetJ
          real*8, dimension(:,:,:), pointer, intent(in) :: xl
          real*8, dimension(nqpt),           intent(in) :: qwt
          integer, intent(in) :: lcsyst,intp
c
          real*8 :: temp_len(npro)
          real*8, dimension(npro, nsd) :: v1, v2, temp_normal
          integer :: isd, iel
          character(len=8) :: err_msg
c
c      write(*,*) 'In calc_normal_vectors...'
c
c.... compute the normal to the boundary. This is achieved by taking
c     the cross product of two vectors in the plane of the 2-d 
c     boundary face.
c
          do isd = 1,nsd
            v1(:,isd) = xl(:,2,isd) - xl(:,1,isd)
            v2(:,isd) = xl(:,3,isd) - xl(:,1,isd)
          enddo
c
          select case (lcsyst)
          case (1)
c
            temp_normal(:,1) = + v1(:,2)*v2(:,3) - v1(:,3)*v2(:,2)
            temp_normal(:,2) = - v1(:,1)*v2(:,3) + v1(:,3)*v2(:,1)
            temp_normal(:,3) = + v1(:,1)*v2(:,2) - v1(:,2)*v2(:,1)
c
            area = pt50 * sqrt(temp_normal(:,1)*temp_normal(:,1)
     &                       + temp_normal(:,2)*temp_normal(:,2)
     &                       + temp_normal(:,3)*temp_normal(:,3))
c
          case default
c            write(err_msg,'(a,i1)') 'lcsyst ',lcsyst
            call error ('calc_normal_vectors ', err_msg, 0)
          end select
c
          do iel = 1,npro
            temp_len(iel) = sqrt(dot_product(temp_normal(iel,1:3),temp_normal(iel,1:3)))
            nv(iel,1:nsd) = temp_normal(iel,1:nsd) / temp_len(iel)
          enddo
c
c... also calculate WdetJ here...
c
          select case (lcsyst)
          case (1)
            WdetJ = qwt(intp) * temp_len * pt25
          case default
            write(err_msg,'(a,i1)') 'lcsyst ',lcsyst
            call error ('calc WdetJif ', err_msg, 0)
          end select
c
        end subroutine calc_normal_vectors
c
        subroutine e3if_flux
c
          integer :: iel,iflow,isd,jsd,jflow
          real*8, dimension(nsd,nflow) :: f0, f1, fconv0, fconv1, fdiff0, fdiff1
          real*8, dimension(nflow) :: f0n0, f0n1, f1n0, f1n1
c
          real*8 :: etot, diff_flux(nsd,nflow), dTdx0
c
          do iel = 1,npro
c
            call calc_conv_flux(fconv0,rho0(iel),u0(iel,:),um0(iel,:),pres0(iel),ei0(iel))
            call calc_conv_flux(fconv1,rho1(iel),u1(iel,:),um1(iel,:),pres1(iel),ei1(iel))
c
            call calc_diff_flux(fdiff0,var0(iel),prop0(iel))
            call calc_diff_flux(fdiff1,var1(iel),prop1(iel))
c
c... calculate flux in normal direction...
c
            do iflow = 1,nflow
c
              f0(:,iflow) = fconv0(:,iflow) - fdiff0(:,iflow)
              f1(:,iflow) = fconv1(:,iflow) - fdiff1(:,iflow)
c
              f0n0(iflow) = dot_product(f0(:,iflow),nv0(iel,:))
              f0n1(iflow) = dot_product(f0(:,iflow),nv1(iel,:))
              f1n0(iflow) = dot_product(f1(:,iflow),nv0(iel,:))
              f1n1(iflow) = dot_product(f1(:,iflow),nv1(iel,:))
c
            enddo
c
c      if (iel == 1) then
c        write(*,10) 'rho0, u0, p0, T0, ei0:', rho0(iel), u0(iel,1), pres0(iel), ei0(iel)
c        write(*,10) 'conv0: ',fconv0(1,:)
c        write(*,10) 'diff0: ',fdiff0(1,:)
c        write(*,10) 'grad_y u:', var0(iel)%grad_y(:,2)
c        write(*,10) 'grad_y T:', var0(iel)%grad_y(:,5)
c        write(*,10) 'f0: ',f0(1,:)
c       write(*,*) 'stiff0: ',prop0(1)%stiff(15,15)
c        write(*,10) 'rho1, u1, p1, ei1:', rho1(iel), u1(iel,1), pres1(iel), ei1(iel)
c        write(*,10) 'conv1: ',fconv1(1,:)
c        write(*,10) 'diff1: ',fdiff1(1,:)
c        write(*,10) 'grad_y u:', var1(iel)%grad_y(:,2)
c        write(*,10) 'grad_y T:', var1(iel)%grad_y(:,5)
c        write(*,10) 'f1: ',f1(1,:)
c       write(*,*) 'stiff1: ',prop1(1)%stiff(15,15)
c      endif
c
            ri0(iel,16:20) = ri0(iel,16:20) + 0.5 * ( f0n0(1:5) + f1n0(1:5) )
            ri1(iel,16:20) = ri1(iel,16:20) + 0.5 * ( f1n1(1:5) + f0n1(1:5) )
c
c...UPWIND????
c
c      ri0(iel,16:20) = ri0(iel,16:20) + f1n0(1:5)
c      ri1(iel,16:20) = ri1(iel,16:20) + f1n1(1:5)
c
          enddo
c
10    format(a,5e24.16)
20    format(a,1e24.16)
c
        end subroutine e3if_flux
c
        subroutine     e3if_vi
c
c... HARDCODED
c    interface progression velocity
c
          real*8, dimension(3), parameter :: vi = (/-1.0e-3,0.,0./)
          real*8 :: vini0,vini1
          real*8 :: etot0,etot1
c
          integer :: iel
c
          do iel = 1,npro
c
            vini0 = dot_product(vi,nv0(iel,:))
            vini1 = dot_product(vi,nv1(iel,:))
c
c...mass
            ri0(iel,1) = ri0(iel,1) + 0.5*(rho0(iel)-rho1(iel))*vini0
            ri1(iel,1) = ri1(iel,1) + 0.5*(rho1(iel)-rho0(iel))*vini1
c
c...momentum
            ri0(iel,2:4) = ri0(iel,2:4) + 0.5*(rho0(iel)*u0(iel,1:3)-rho1(iel)*u1(iel,1:3))*vini0
            ri1(iel,2:4) = ri1(iel,2:4) + 0.5*(rho1(iel)*u1(iel,1:3)-rho0(iel)*u0(iel,1:3))*vini1
c
c...energy
            etot0 = ei0(iel) + 0.5*dot_product(u0(iel,:),u0(iel,:))
            etot1 = ei1(iel) + 0.5*dot_product(u1(iel,:),u1(iel,:))
c
            ri0(iel,5) = ri0(iel,5) + 0.5*(rho0(iel)*etot0-rho1(iel)*etot1)*vini0
            ri1(iel,5) = ri1(iel,5) + 0.5*(rho1(iel)*etot1-rho0(iel)*etot0)*vini1
c
          enddo
c
        end subroutine e3if_vi
c
        subroutine kinematic_condition
c
          integer :: i,j,isd,jsd
          real*8 :: cy_jump(npro,nsd,nflow)
     &,             cy0(npro),cy1(npro)
          real*8, dimension(npro) :: kij_cy0, kij_cy1, cwcy0, cwcy1
c
c... calc jump in C*Y
c
          do i = 1,nflow
c
            cy0(:) = zero
            cy1(:) = zero
c
            do j = 1,nflow
              cy0 = cy0 + cmtrx(:,i,j)*var0(:)%y(j)
              cy1 = cy1 + cmtrx(:,i,j)*var1(:)%y(j)
            enddo
c
            do isd = 1,nsd
              cy_jump(:,isd,i) = cy0 * nv0(:,isd) + cy1 * nv1(:,isd)
            enddo
c
          enddo
c
c      write(*,*) 'cy_jump 1: ',cy_jump(1,1,:)
c      write(*,*) 'cy_jump 2: ',cy_jump(1,2,:)
c      write(*,*) 'cy_jump 3: ',cy_jump(1,3,:)
c
c... Multiply by diffusion matrix and add to RHS residual...
c
          do i = 1, nflow
c
            cwcy0 = zero
            cwcy1 = zero
c
            do isd = 1, nsd
c
c...another set of loops for matrix operation...
c
              kij_cy0 = zero
              kij_cy1 = zero
c
              do j = 1, nflow
                do jsd = 1, nsd
c
                  Kij_cy0 = Kij_cy0 + Kij0(:,isd,jsd,i,j)*cy_jump(:,jsd,j)
                  Kij_cy1 = Kij_cy1 + Kij1(:,isd,jsd,i,j)*cy_jump(:,jsd,j)
c
                  cwcy0 = cwcy0 + cmtrx(:,i,j)*nv0(:,jsd)*cy_jump(:,jsd,j)
                  cwcy1 = cwcy1 + cmtrx(:,i,j)*nv1(:,jsd)*cy_jump(:,jsd,j)
c
                enddo
              enddo
c
              ri0(:,nflow*(isd-1)+i) = ri0(:,nflow*(isd-1)+i) + pt50 * s * kij_cy0
              ri1(:,nflow*(isd-1)+i) = ri1(:,nflow*(isd-1)+i) + pt50 * s * kij_cy1
c
c      if (i == 5 .and. isd == 1) then
c        write(*,*) 'i, isd, kij_cy0: ',i,isd,kij_cy0(1)
c        write(*,*) 'i, isd, kij_cy1: ',i,isd,kij_cy1(1)
c      endif
c
            enddo
c
            ri0(:,nflow*3+i) = ri0(:,nflow*3+i) + e*mu/h * cwcy0
            ri1(:,nflow*3+i) = ri1(:,nflow*3+i) + e*mu/h * cwcy1
c
c      write(*,*) 'i, cwcy0: ',i,cwcy0(1)
c      write(*,*) 'i, cwcy1: ',i,cwcy1(1)
c
          enddo
c
c      write(*,20) 'kij_cy0, kij_cy1:', kij_cy0(1), kij_cy1(1)
20    format(a,2e24.16)
c
        end subroutine kinematic_condition
c
        subroutine calc_cmtrx
c
          integer :: p,q,r
c
          cmtrx = zero
c
          cmtrx(:,2,2) = one - nv0(:,1)*nv0(:,1)
          cmtrx(:,2,3) =     - nv0(:,1)*nv0(:,2)
          cmtrx(:,2,4) =     - nv0(:,1)*nv0(:,3)
c
          cmtrx(:,3,2) =     - nv0(:,2)*nv0(:,1)
          cmtrx(:,3,3) = one - nv0(:,2)*nv0(:,2)
          cmtrx(:,3,4) =     - nv0(:,2)*nv0(:,3)
c
          cmtrx(:,4,2) =     - nv0(:,3)*nv0(:,1)
          cmtrx(:,4,3) =     - nv0(:,3)*nv0(:,2)
          cmtrx(:,4,4) = one - nv0(:,3)*nv0(:,3)
c
          cmtrx(:,5,5) = one
c
c...NOTE: need to derive an expression for it, instead of this:
c
          do q = 1,nflow
            do p = 1,nflow
              ctc(:,p,q) = zero
              do r = 1,nflow
                ctc(:,p,q) = ctc(:,p,q) + cmtrx(:,r,p)*cmtrx(:,r,q)
              enddo
            enddo
          enddo
c
        end subroutine calc_cmtrx
c
        subroutine calc_conv_flux(flux,rho,u,um,p,ei)
c
          real*8, dimension(nsd,nflow), intent(out) :: flux
          real*8, intent(in) :: rho,p,ei,u(nsd),um(nsd)
          real*8 :: delta(3,3)
          data delta /1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
c
          integer :: i
          real*8  :: etot
c
          etot = ei + pt50 * dot_product(u,u)
c
          do i = 1,nsd
c
            flux(i,1) = rho*(u(i)-um(i))
            flux(i,2) = rho*(u(i)-um(i))*u(1) + p*delta(1,i)
            flux(i,3) = rho*(u(i)-um(i))*u(2) + p*delta(2,i)
            flux(i,4) = rho*(u(i)-um(i))*u(3) + p*delta(3,i)
            flux(i,5) = rho*(u(i)-um(i))*etot + p*u(i)
c
          enddo
c
        end subroutine calc_conv_flux
c
        subroutine calc_diff_flux(flux,var,prop)
c
          real*8, dimension(nsd,nflow), intent(out) :: flux
          type(var_t), intent(in) :: var
          type(prop_t), intent(in) :: prop
c
          integer :: isd,jsd,iflow,jflow
          real*8  :: this_sum
c
c... diffusive flux
c
            do iflow = 1,nflow
              do isd = 1,nsd
                this_sum = zero
                ! multiply matrix: K_ij * Y,j
                do jflow = 1,nflow
                  do jsd = 1,nsd
                    this_sum = this_sum +
     &                prop%stiff((isd-1)*nflow+iflow,(jsd-1)*nflow+jflow)*var%grad_y(jsd,jflow)
                  enddo
                enddo
                flux(isd,iflow) = this_sum
              enddo
            enddo
c
        end subroutine calc_diff_flux
c
        subroutine e3if_wmlt(rl,ri,shp,shg,WdetJ,nshl)
c
          integer, intent(in) :: nshl
          real*8, dimension(:,:,:), pointer, intent(out) :: rl
          real*8, dimension(:,:),   pointer, intent(in)  :: ri,shp
          real*8, dimension(:,:,:), pointer, intent(in)  :: shg
          real*8, dimension(:),     pointer, intent(in)  :: WdetJ
c
          integer :: n,iflow,isd
c
          do isd = 1,nsd
            do iflow = 1,nflow
              do n = 1, nshl
                rl(:,n,iflow) = rl(:,n,iflow) +
     &            WdetJ * shg(:,n,isd) * ri(:,nflow*(isd-1)+iflow)
              enddo
            enddo
          enddo
c
          do iflow = 1,nflow
            do n = 1,nshl
              rl(:,n,iflow) = rl(:,n,iflow) + shp(:,n)*WdetJ(:) * ri(:,3*nflow+iflow)
            enddo
          enddo
c
        end subroutine e3if_wmlt
c
      end module e3if_m
