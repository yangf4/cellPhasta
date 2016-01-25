      subroutine getdiff(rmu, rlm, rlm2mu, con, npro, mater)
c
        use number_def_m
        use matdat_def_m
c
        implicit none
c
        integer, intent(in) :: npro, mater
        real*8, dimension(npro), intent(out) :: rmu,rlm,rlm2mu,con
c
        integer :: ivisc, icon
c
c HARDCODED until the proper indexing system is implemented
c
        select case (mater)
        case (1)
          ivisc = 3
          icon = 4
        case (2)
          ivisc = 7
          icon = 8
        case default
          call error ('getdiff ', 'ERROR: index can not be set!',0)
        end select
c
        if (matflg(2,mater) == 0) then ! Shear Law: Constant Viscosity
          rmu = mat_prop(mater,ivisc,1)
        else
          call error ('getdiff ', 'ERROR: Constant Viscosity is supported ONLY!', 0)
        endif
c
        if (matflg(3,mater) == 0) then ! Bulk Viscosity Law: Constant Bulk Viscosity
          rlm = -pt66 * rmu
        else
          call error ('getdiff ', 'ERROR: Constant Bulk Viscosity is supported ONLY!', 0)
        endif
c
        rlm2mu = rlm + two * rmu
c
        con = mat_prop(mater,icon,1)
c
c        if (iLES .gt. 0 .or. iRANS.eq.0) then
c          call error ('getdiff ', 'ERROR: Turbulence Viscosity is NOT supported!', 0)
c        endif
c
      end subroutine getdiff
c
      subroutine getdiffsclr
      end subroutine getdiffsclr
