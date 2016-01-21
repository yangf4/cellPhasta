      module matdat_def_m
c
        implicit none
c
        integer, parameter :: MAXMAT  = 6
     &,                       MAXPROP = 10
     &,                       MAXTS   = 100
        real*8  :: datmat(MAXMAT,MAXPROP,MAXTS)
        integer :: nummat
     &,            matflg(MAXMAT,MAXTS)
     &,            mat_tag(MAXMAT,MAXTS)
     &,            mat_eos(MAXMAT,MAXTS)
        real*8  :: mat_prop(MAXMAT,MAXPROP,MAXTS)
        logical :: mexist
        real*8  :: datelas(1,2)
        common /matdat/ datmat,   matflg, mat_tag, mat_eos,
     &                  mat_prop, nummat, mexist,  datelas
        integer, parameter :: ieos_ideal_gas = 1
     &,                       ieos_liquid_1  = 2
c
        integer, parameter :: iprop_ideal_gas_mw    = 1  ! molecular weight
     &,                       iprop_ideal_gas_gamma = 2
     &,                       iprop_gas_visc  = 3
     &,                       iprop_liquid_1_rho_ref = 1
     &,                       iprop_liquid_1_p_ref   = 2
     &,                       iprop_liquid_1_T_ref   = 3
     &,                       iprop_liquid_1_cv      = 4
     &,                       iprop_liquid_1_alphaP  = 5
     &,                       iprop_liquid_1_betaT   = 6
c
      end module matdat_def_m
