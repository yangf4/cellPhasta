        subroutine getthm (rho,    ei
     &,                    p,      T,     npro, mater
     &,                    h,      cv,    cp
     &,                    alphaP, betaT, gamb, c)
c
          use eqn_state_m
c
          implicit none
c
          integer, intent(in) :: npro, mater
          real*8, dimension(npro), intent(in) :: p,T
          real*8, dimension(npro), intent(out) :: rho,ei
          real*8, dimension(npro), optional, intent(out) :: h,cv,cp,alphaP,betaT,gamb,c
c 
          select case (mat_eos(mater,1))
          case (ieos_ideal_gas)
c
            call getthm_ideal_gas(rho,ei,p,T,npro,mater,
     &        h, cv, cp, alphaP, betaT, gamb, c)
          case (ieos_liquid_1)
            call getthm_liquid_1(rho,ei,p,T,npro,mater,
     &        h, cv, cp, alphaP, betaT, gamb, c)
          case default
            call error ('getthm  ', 'wrong material', mater)
          end select
c
        end subroutine getthm
