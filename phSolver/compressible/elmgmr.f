        subroutine ElmGMRe (y,         ac,        x,         
     &                     shp,       shgl,      iBC,
     &                     BC,        shpb,      shglb,
     &                     res,       rmes,      BDiag,
     &                     iper,      ilwork,    EGmass, rerr)
c
c----------------------------------------------------------------------
c
c This routine computes the LHS mass matrix, the RHS residual 
c vector, and the preconditioning matrix, for use with the GMRES
c solver.
c
c Zdenek Johan, Winter 1991.      (Fortran 90)
c Chris Whiting, Winter 1998.     (Matrix EBE-GMRES)
c----------------------------------------------------------------------
c
        use pointer_data
        use timedataC
c
        include "common.h"
        include "mpif.h"
c
        dimension y(nshg,ndof),         ac(nshg,ndof),
     &            x(numnp,nsd),               
     &            iBC(nshg),
     &            BC(nshg,ndofBC),      
     &            res(nshg,nflow),
     &            rmes(nshg,nflow),      BDiag(nshg,nflow,nflow),
     &            iper(nshg),           EGmass(numel,nedof,nedof)
c
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
        dimension qres(nshg, idflx),     rmass(nshg)
c
        dimension ilwork(nlwork)

        real*8 Bdiagvec(nshg,nflow), rerr(nshg,10)

        real*8, allocatable :: tmpshp(:,:), tmpshgl(:,:,:)
        real*8, allocatable :: tmpshpb(:,:), tmpshglb(:,:,:)

	ttim(80) = ttim(80) - secs(0.0)
c
c.... set up the timer
c

        call timer ('Elm_Form')
c
c.... -------------------->   interior elements   <--------------------
c
c.... set up parameters
c
        ires   = 1
c
        if (idiff==1 .or. idiff==3 .or. isurf==1) then ! global reconstruction
                                                       ! of qdiff
c
c loop over element blocks for the global reconstruction
c of the diffusive flux vector, q, and lumped mass matrix, rmass
c
        qres = zero
        rmass = zero
        
        do iblk = 1, nelblk
c
c.... set up the parameters
c
          nenl   = lcblk(5,iblk)   ! no. of vertices per element
          iel    = lcblk(1,iblk)
          lelCat = lcblk(2,iblk)
          lcsyst = lcblk(3,iblk)
          iorder = lcblk(4,iblk)
          nenl   = lcblk(5,iblk)   ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          mattyp = lcblk(7,iblk)
          ndofl  = lcblk(8,iblk)
          nsymdl = lcblk(9,iblk)
          npro   = lcblk(1,iblk+1) - iel
          ngauss = nint(lcsyst)
c
c.... compute and assemble diffusive flux vector residual, qres,
c     and lumped mass matrix, rmass

          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))

          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)

          call AsIq (y,                x,                       
     &               tmpshp,              
     &               tmpshgl,
     &               mien(iblk)%p,     mxmudmi(iblk)%p,
     &               qres,                   
     &               rmass)

          deallocate ( tmpshp )
          deallocate ( tmpshgl ) 
       enddo
       
c
c.... take care of periodic boundary conditions
c

       call qpbc( rmass, qres, iBC,  iper, ilwork )       
c
      endif                     ! computation of global diffusive flux
c
c.... loop over element blocks to compute element residuals
c
c
c.... initialize the arrays
c
        res    = zero
        rmes   = zero ! to avoid trap_uninitialized
        if (lhs. eq. 1)   EGmass = zero
        if (iprec .ne. 0) BDiag = zero
        flxID = zero

c
c.... loop over the element-blocks
c
        do iblk = 1, nelblk
c
c.... set up the parameters
c
          iblkts = iblk          ! used in timeseries
          nenl   = lcblk(5,iblk) ! no. of vertices per element
          iel    = lcblk(1,iblk)
          lelCat = lcblk(2,iblk)
          lcsyst = lcblk(3,iblk)
          iorder = lcblk(4,iblk)
          nenl   = lcblk(5,iblk)   ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          mattyp = lcblk(7,iblk)
          ndofl  = lcblk(8,iblk)
          nsymdl = lcblk(9,iblk)
          npro   = lcblk(1,iblk+1) - iel 
          inum   = iel + npro - 1
          ngauss = nint(lcsyst)
c
c.... compute and assemble the residual and tangent matrix
c

          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))

          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)

          call AsIGMR (y,                   ac,
     &                 x,                   mxmudmi(iblk)%p,
     &                 tmpshp,
     &                 tmpshgl,             mien(iblk)%p,
     &                 mmat(iblk)%p,        res,
     &                 rmes,                BDiag,
     &                 qres,                EGmass(iel:inum,:,:),
     &                 rerr)
c
c.... satisfy the BC's on the implicit LHS
c     
          call bc3LHS (iBC,                  BC,  mien(iblk)%p, 
     &                 EGmass(iel:inum,:,:)  ) 

          deallocate ( tmpshp )
          deallocate ( tmpshgl )
c
c.... end of interior element loop
c
       enddo
c
c.... -------------------->   boundary elements   <--------------------
c
c.... loop over the boundary elements
c
        do iblk = 1, nelblb
c
c.... set up the parameters
c
          iel    = lcblkb(1,iblk)
          lelCat = lcblkb(2,iblk)
          lcsyst = lcblkb(3,iblk)
          iorder = lcblkb(4,iblk)
          nenl   = lcblkb(5,iblk)  ! no. of vertices per element
          nenbl  = lcblkb(6,iblk)  ! no. of vertices per bdry. face
          mattyp = lcblkb(7,iblk)
          ndofl  = lcblkb(8,iblk)
          nshl   = lcblkb(9,iblk)
          nshlb  = lcblkb(10,iblk)
          npro   = lcblkb(1,iblk+1) - iel 
          if(lcsyst.eq.3) lcsyst=nenbl
          ngaussb = nintb(lcsyst)
          
c
c.... compute and assemble the residuals corresponding to the 
c     boundary integral
c

          allocate (tmpshpb(nshl,MAXQPT))
          allocate (tmpshglb(nsd,nshl,MAXQPT))
          
          tmpshpb(1:nshl,:) = shpb(lcsyst,1:nshl,:)
          tmpshglb(:,1:nshl,:) = shglb(lcsyst,:,1:nshl,:)

          call AsBMFG (y,                       x,
     &                 tmpshpb,                 tmpshglb, 
     &                 mienb(iblk)%p,           mmatb(iblk)%p,
     &                 miBCB(iblk)%p,           mBCB(iblk)%p,
     &                 res,                     rmes)

          deallocate (tmpshpb)
          deallocate (tmpshglb)
c
c.... end of boundary element loop
c
        enddo
c
      ttim(80) = ttim(80) + secs(0.0)
c
c before the commu we need to rotate the residual vector for axisymmetric
c boundary conditions (so that off processor periodicity is a dof add instead
c of a dof combination).  Take care of all nodes now so periodicity, like
c commu is a simple dof add.
c
c      if(iabc==1)               !are there any axisym bc's
c     &     call rotabc(res(1,2), iBC, BC, nflow,  'in ')
      if(iabc==1) then               !are there any axisym bc's
          call rotabc(res(1,2), iBC,  'in ')
c          Bdiagvec(:,1)=BDiag(:,1,1)
c          Bdiagvec(:,2)=BDiag(:,2,2)
c          Bdiagvec(:,3)=BDiag(:,3,3)
c          Bdiagvec(:,4)=BDiag(:,4,4)
c          Bdiagvec(:,5)=BDiag(:,5,5)
c          call rotabc(Bdiagvec(1,2), iBC, BC, 2,  'in ')
c          BDiag(:,:,:)=zero
c          BDiag(:,1,1)=Bdiagvec(:,1)
c          BDiag(:,2,2)=Bdiagvec(:,2)
c          BDiag(:,3,3)=Bdiagvec(:,3)
c          BDiag(:,4,4)=Bdiagvec(:,4)
c          BDiag(:,5,5)=Bdiagvec(:,5)
       endif

c.... -------------------->   communications <-------------------------
c

      if (numpe > 1) then
        call commu (res  , ilwork, nflow  , 'in ')

        call MPI_BARRIER (MPI_COMM_WORLD,ierr)

        if(iprec .ne. 0) call commu (BDiag, ilwork, nflow*nflow, 'in ')
      endif

c
c.... ---------------------->   post processing  <----------------------
c
c.... satisfy the BCs on the residual
c
      call bc3Res (y,  iBC,  BC,  res,  iper, ilwork)
c
c.... satisfy the BCs on the block-diagonal preconditioner
c
      if (iprec .ne. 0) then
         call bc3BDg (y,  iBC,  BC,  BDiag, iper, ilwork)
      endif
c
c.... return
c
      call timer ('Back    ')
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccc       SPARSE   
c_______________________________________________________________

        subroutine ElmGMRs (y,         ac,        x,         
     &                     shp,       shgl,      iBC,
     &                     BC,        shpb,      shglb,
     &                     shpif0,    shpif1,    shgif0,   shgif1,
     &                     res,       rmes,      BDiag,
     &                     iper,      ilwork,    lhsK,  
     &                     col,       row,       rerr,     umesh)
c
c----------------------------------------------------------------------
c
c This routine computes the LHS mass matrix, the RHS residual 
c vector, and the preconditioning matrix, for use with the GMRES
c solver.
c
c Zdenek Johan, Winter 1991.      (Fortran 90)
c Chris Whiting, Winter 1998.     (Matrix EBE-GMRES)
c----------------------------------------------------------------------
c
        use pointer_data
        use timedataC
        use if_velocity_m
c
        include "common.h"
        include "mpif.h"
c#define DEBUG
c
        interface 
          subroutine asidgif_
     & (
     &   nshg_, nshl0_, nshl1_, nenl0_, nenl1_,lcsyst0_,lcsyst1_,
     &   nflow_, npro_, ndof_, nsd_, ipord_, numnp_, nqpt_,
     &   gbytes_, sbytes_, flops_,
     &   res,
     &   egmassif00,egmassif01,egmassif10,egmassif11,
     &   y,        x,       umesh,
     &   shpif0,   shpif1,  shgif0,  shgif1,
     &   qwtif0,   qwtif1,
     &   ienif0,   ienif1,
     &   mattypeif0, mattypeif1,
     &   time,
     &   sum_vi_area
     & )
         use hierarchic_m
         use local_m
         use e3if_m
            implicit none
            integer, intent(in) :: nshg_, nshl0_, nshl1_, nenl0_, nenl1_,lcsyst0_,lcsyst1_
            integer, intent(in) :: nflow_, npro_, ndof_, nsd_, ipord_, numnp_, nqpt_
            integer, intent(inout) :: gbytes_, sbytes_, flops_
            real*8, dimension(nshg_,nflow_), intent(inout) :: res
            real*8, dimension(:,:,:), pointer, intent(out) :: egmassif00,egmassif01,egmassif10,egmassif11
            real*8, dimension(nshg_,ndof_),  intent(in)    :: y
            real*8, dimension(nshg_,nsd_),   intent(in)    :: x
            real*8, dimension(nshl0_,nqpt_),intent(in)   :: shpif0
            real*8, dimension(nshl1_,nqpt_),intent(in)   :: shpif1
            real*8, dimension(nsd_,nshl0_,nqpt_), intent(in)  :: shgif0
            real*8, dimension(nsd_,nshl1_,nqpt_), intent(in)  :: shgif1
            real*8, dimension(nqpt_), intent(in) :: qwtif0, qwtif1
            real*8, dimension(numnp_, nsd_), intent(inout) :: umesh
            integer, dimension(:,:), pointer, intent(in)   :: ienif0, ienif1
            integer, intent(in)   :: mattypeif0, mattypeif1
            real*8, intent(in) :: time
            real*8, pointer, intent(inout) :: sum_vi_area(:,:)
          end subroutine asidgif_
          subroutine fillsparse_if
     &    ( lhsk,
     &      ienif0,ienif1,
     &      col,row,
     &      egmass,
     &      nflow,nshg,nnz,nnz_tot)
            implicit none
            real*8, intent(inout) :: lhsK(nflow*nflow,nnz_tot)
            integer, dimension(:,:), pointer, intent(in) :: ienif0,ienif1
            integer, intent(in) :: col(nshg+1), row(nnz*nshg)
            real*8, dimension(:,:,:), pointer, intent(in) :: egmass
            integer, intent(in) :: nflow,nshg,nnz,nnz_tot
          end subroutine fillsparse_if
        end interface
c
        integer col(nshg+1), row(nnz*nshg)
        real*8 lhsK(nflow*nflow,nnz_tot)
        
        dimension y(nshg,ndof),         ac(nshg,ndof),
     &            x(numnp,nsd),               
     &            iBC(nshg),
     &            BC(nshg,ndofBC),      
     &            res(nshg,nflow),
     &            rmes(nshg,nflow),      BDiag(nshg,nflow,nflow),
     &            iper(nshg)
c
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
        real*8, dimension(maxtopif,    maxsh,maxqpt) :: shpif0, shpif1
        real*8, dimension(maxtopif,nsd,maxsh,maxqpt) :: shgif0, shgif1
c
        dimension qres(nshg, idflx),     rmass(nshg)
c
        dimension ilwork(nlwork)
c  
        dimension umesh(numnp, nsd)
c
        real*8 Bdiagvec(nshg,nflow), rerr(nshg,10)

        real*8, allocatable :: tmpshp(:,:), tmpshgl(:,:,:)
        real*8, allocatable :: tmpshpb(:,:), tmpshglb(:,:,:)
        real*8, allocatable :: EGmass(:,:,:)
c
        real*8, dimension(:,:,:), pointer :: egmassif00,egmassif01,egmassif10,egmassif11
c
        ttim(80) = ttim(80) - secs(0.0)
c
c.... set up the timer
c

        call timer ('Elm_Form')
c
c.... -------------------->   interior elements   <--------------------
c
c.... set up parameters
c
        ires   = 1
c
        if (idiff==1 .or. idiff==3 .or. isurf==1) then ! global reconstruction
                                                       ! of qdiff
c
c loop over element blocks for the global reconstruction
c of the diffusive flux vector, q, and lumped mass matrix, rmass
c
        qres = zero
        rmass = zero
        
        do iblk = 1, nelblk
c
c.... set up the parameters
c
          nenl   = lcblk(5,iblk)   ! no. of vertices per element
          iel    = lcblk(1,iblk)
          lelCat = lcblk(2,iblk)
          lcsyst = lcblk(3,iblk)
          iorder = lcblk(4,iblk)
          nshl   = lcblk(10,iblk)
          mattyp = lcblk(7,iblk)
          ndofl  = lcblk(8,iblk)
          nsymdl = lcblk(9,iblk)
          npro   = lcblk(1,iblk+1) - iel
          ngauss = nint(lcsyst)
c
c.... compute and assemble diffusive flux vector residual, qres,
c     and lumped mass matrix, rmass

          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))

          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)

          call AsIq (y,                x,                       
     &               tmpshp,              
     &               tmpshgl,
     &               mien(iblk)%p,     mxmudmi(iblk)%p,
     &               qres,                   
     &               rmass)

          deallocate ( tmpshp )
          deallocate ( tmpshgl ) 
       enddo
       
c
c.... take care of periodic boundary conditions
c

       call qpbc( rmass, qres, iBC, iper, ilwork )       
c
      endif                     ! computation of global diffusive flux
c
c.... loop over element blocks to compute element residuals
c
c
c.... initialize the arrays
c
        res    = zero
        rmes   = zero ! to avoid trap_uninitialized
        if (lhs. eq. 1) lhsK = zero
        if (iprec .ne. 0) BDiag = zero
        flxID = zero
c
c.... loop over the element-blocks
c
        do iblk = 1, nelblk
c
c.... set up the parameters
c
          iblkts = iblk          ! used in timeseries
          nenl   = lcblk(5,iblk) ! no. of vertices per element
          iel    = lcblk(1,iblk)
          lelCat = lcblk(2,iblk)
          lcsyst = lcblk(3,iblk)
          iorder = lcblk(4,iblk)
          nenl   = lcblk(5,iblk)   ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          mattyp = lcblk(7,iblk)
          ndofl  = lcblk(8,iblk)
          nsymdl = lcblk(9,iblk)
          npro   = lcblk(1,iblk+1) - iel 
          inum   = iel + npro - 1
          ngauss = nint(lcsyst)
c
c.... compute and assemble the residual and tangent matrix
c

          if(lhs.eq.1) then
             allocate (EGmass(npro,nedof,nedof))
             EGmass = zero
          else
             allocate (EGmass(1,1,1))
          endif

          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))
          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)

          call AsIGMR (y,                   ac,
     &                 x,                   mxmudmi(iblk)%p,
     &                 tmpshp,
     &                 tmpshgl,             mien(iblk)%p,
     &                 mattyp,              res,
     &                 rmes,                BDiag,
     &                 qres,                EGmass,
     &                 rerr,                umesh )
          if(lhs.eq.1) then
c
c.... satisfy the BC's on the implicit LHS
c     
             call bc3LHS (iBC,                  BC,  mien(iblk)%p, 
     &                    EGmass  ) 

c
c.... Fill-up the global sparse LHS mass matrix
c
             call fillsparseC( mien(iblk)%p, EGmass,
     1                        lhsK, row, col)
          endif
c
          deallocate ( EGmass )
          deallocate ( tmpshp )
          deallocate ( tmpshgl )
c
c.... end of interior element loop
c
       enddo
!ifdef DEBUG !Nicholas Mati
!        call write_debug(myrank, 'res-afterAsIGMR'//char(0),
!     &                           'res'//char(0), res, 
!     &                           'd'//char(0), nshg, nflow, lstep)
!        call write_debug(myrank, 'y-afterAsIGMR'//char(0),
!     &                           'y'//char(0), y, 
!     &                           'd'//char(0), nshg, ndof, lstep)
!endif //DEBUG
c
!      write(*,998) '[',myrank,'] in elmgmr AFTER INTERIOR.'
!      do i = 1,nshg
!        if (x(i,1) < 0.1001 .and. x(i,1) > 0.0999 .and. 
!     &      x(i,2) < 0.0501 .and. x(i,2) > 0.0499 .and.
!     &      x(i,3) < 0.0001 .and. x(i,3) > -0.0001)
!     &   write(*,999) '[',myrank,'] :',i,x(i,:),res(i,:)
!      enddo
c
c.... -------------------->   boundary elements   <--------------------
c
c.... loop over the boundary elements
c
        do iblk = 1, nelblb
c
c.... set up the parameters
c
          iel    = lcblkb(1,iblk)
          lelCat = lcblkb(2,iblk)
          lcsyst = lcblkb(3,iblk)
          iorder = lcblkb(4,iblk)
          nenl   = lcblkb(5,iblk)  ! no. of vertices per element
          nenbl  = lcblkb(6,iblk)  ! no. of vertices per bdry. face
          mattyp = lcblkb(7,iblk)
          ndofl  = lcblkb(8,iblk)
          nshl   = lcblkb(9,iblk)
          nshlb  = lcblkb(10,iblk)
          npro   = lcblkb(1,iblk+1) - iel 
          if(lcsyst.eq.3) lcsyst=nenbl
          ngaussb = nintb(lcsyst)
          
c
c.... compute and assemble the residuals corresponding to the 
c     boundary integral
c

          allocate (tmpshpb(nshl,MAXQPT))
          allocate (tmpshglb(nsd,nshl,MAXQPT))
          if(lhs.eq.1 .and. iLHScond >= 1) then
             allocate (EGmass(npro,nshl,nshl))
             EGmass = zero
          else
             allocate (EGmass(1,1,1))
          endif
          
          tmpshpb(1:nshl,:) = shpb(lcsyst,1:nshl,:)
          tmpshglb(:,1:nshl,:) = shglb(lcsyst,:,1:nshl,:)

          call AsBMFG (y,                       x,
     &                 tmpshpb,                 tmpshglb, 
     &                 mienb(iblk)%p,           mattyp,
     &                 miBCB(iblk)%p,           mBCB(iblk)%p,
     &                 res,                     rmes, 
     &                 EGmass,                  umesh)
          if(lhs == 1 .and. iLHScond > 0) then
            call fillSparseC_BC(mienb(iblk)%p, EGmass, 
     &                   lhsk, row, col)
          endif

          deallocate (EGmass)
          deallocate (tmpshpb)
          deallocate (tmpshglb)
        enddo   !end of boundary element loop
    
!ifdef DEBUG !Nicholas Mati
!        call write_debug(myrank, 'res-afterAsBMFG'//char(0),
!     &                           'res'//char(0), res, 
!     &                           'd'//char(0), nshg, nflow, lstep)
!        call MPI_ABORT(MPI_COMM_WORLD) 
!endif //DEBUG
c
!      write(*,998) '[',myrank,'] in elmgmr AFTER BOUNDARY.'
!      do i = 1,nshg
!        if (x(i,1) < 0.1001 .and. x(i,1) > 0.0999 .and. 
!     &      x(i,2) < 0.0501 .and. x(i,2) > 0.0499 .and.
!     &      x(i,3) < 0.0001 .and. x(i,3) > -0.0001)
!     &   write(*,999) '[',myrank,'] :',i,x(i,:),res(i,:)
!      enddo
c
      ttim(80) = ttim(80) + secs(0.0)
c
c.... -------------------->   interface elements   <--------------------
c
c... loop over the interface element blocks
c
        sum_vi_area = zero
c
        if_blocks: do iblk = 1, nelblif
c
c... set up the parameters
c
          iblkts  = iblk                ! used in time series
          nenl0   = lcblkif(6, iblk)    ! number of vertices per element0
          nenl1   = lcblkif(7, iblk)    ! number of vertices per element1
          iel     = lcblkif(1, iblk)
          lelCat  = lcblkif(2, iblk)    ! ??? NOT USED?
          lcsyst0 = lcblkif(3, iblk)    ! element0 type
          lcsyst1 = lcblkif(4, iblk)    ! element1 type
          iorder  = lcblkif(5, iblk)    ! polynomial order
          nshl0   = lcblkif(13,iblk)
          nshl1   = lcblkif(14,iblk)
          mattyp0 = lcblkif(9, iblk)
          mattyp1 = lcblkif(10,iblk)
          ndof    = lcblkif(11,iblk)
          nsymdl  = lcblkif(12,iblk)    ! ???
          npro    = lcblkif(1,iblk+1) - iel
          inum    = iel + npro - 1
          ngaussif = nintif0(lcsyst0)   ! or nintif1(lcsyst1)? should be the same!
c
c... compute and assemble the residual and tangent matrix
c
          if (lhs .eq. 1) then
            allocate (egmassif00(npro,nflow*nshl0,nflow*nshl0))
            allocate (egmassif01(npro,nflow*nshl0,nflow*nshl1))
            allocate (egmassif10(npro,nflow*nshl1,nflow*nshl0))
            allocate (egmassif11(npro,nflow*nshl1,nflow*nshl1))
            egmassif00 = zero
            egmassif01 = zero
            egmassif10 = zero
            egmassif11 = zero
          else
            allocate (egmassif00(1,1,1))
            allocate (egmassif01(1,1,1))
            allocate (egmassif10(1,1,1))
            allocate (egmassif11(1,1,1))
          endif

      call asidgif_
     & (
     &   nshg, nshl0, nshl1, nenl0, nenl1, lcsyst0, lcsyst1,
     &   nflow, npro, ndof, nsd, iorder, numnp, ngaussif,
     &   gbytes, sbytes, flops,
     &         res,
     &         egmassif00,egmassif01,egmassif10,egmassif11,
     &         y, x, umesh,
     &         shpif0(lcsyst0,1:nshl0,:), 
     &         shpif1(lcsyst1,1:nshl1,:), 
     &         shgif0(lcsyst0,1:nsd,1:nshl0,:),
     &         shgif1(lcsyst1,1:nsd,1:nshl1,:),
     &         qwtif0(lcsyst0,:), qwtif1(lcsyst1,:),
     &         mienif0(iblk)%p, mienif1(iblk)%p,
     &         mattyp0, mattyp1,
     &         time,
     &         sum_vi_area
     & )
c
          if (lhs .eq. 1) then
c
c.... Fill-up the global sparse LHS mass matrix
c
c            call fillsparse_if( lhsk,mienif0(iblk)%p,mienif0(iblk)%p,col,row,egmassif00,nflow,nshg,nnz,nnz_tot)
c            call fillsparse_if( lhsk,mienif0(iblk)%p,mienif1(iblk)%p,col,row,egmassif01,nflow,nshg,nnz,nnz_tot)
c            call fillsparse_if( lhsk,mienif1(iblk)%p,mienif0(iblk)%p,col,row,egmassif10,nflow,nshg,nnz,nnz_tot)
c            call fillsparse_if( lhsk,mienif1(iblk)%p,mienif1(iblk)%p,col,row,egmassif11,nflow,nshg,nnz,nnz_tot)
c
          endif

          deallocate (egmassif00)
          deallocate (egmassif01)
          deallocate (egmassif10)
          deallocate (egmassif11)


c
        enddo if_blocks
c
c
c before the commu we need to rotate the residual vector for axisymmetric
c boundary conditions (so that off processor periodicity is a dof add instead
c of a dof combination).  Take care of all nodes now so periodicity, like
c commu is a simple dof add.
c
      if(iabc==1) then               !are there any axisym bc's
          call rotabc(res(1,2), iBC,  'in ')
c          Bdiagvec(:,1)=BDiag(:,1,1)
c          Bdiagvec(:,2)=BDiag(:,2,2)
c          Bdiagvec(:,3)=BDiag(:,3,3)
c          Bdiagvec(:,4)=BDiag(:,4,4)
c          Bdiagvec(:,5)=BDiag(:,5,5)
c          call rotabc(Bdiagvec(1,2), iBC,  'in ')
c          BDiag(:,:,:)=zero
c          BDiag(:,1,1)=Bdiagvec(:,1)
c          BDiag(:,2,2)=Bdiagvec(:,2)
c          BDiag(:,3,3)=Bdiagvec(:,3)
c          BDiag(:,4,4)=Bdiagvec(:,4)
c          BDiag(:,5,5)=Bdiagvec(:,5)
       endif

c.... -------------------->   communications <-------------------------
c
      if (numpe > 1) then
        call commu (res  , ilwork, nflow  , 'in ')

        call MPI_BARRIER (MPI_COMM_WORLD,ierr)

        if(iprec .ne. 0) call commu (BDiag, ilwork, nflow*nflow, 'in ')
      endif
c
c------> BEGIN DEBUG <---------
c
!      write(*,998) '[',myrank,'] in elmgmr AFTER commu.'
!      do i = 1,nshg
!        if (x(i,1) < 0.1001 .and. x(i,1) > 0.0999 .and. 
!     &      x(i,2) < 0.0501 .and. x(i,2) > 0.0499 .and.
!     &      x(i,3) < 0.0001 .and. x(i,3) > -0.0001)
!     &   write(*,999) '[',myrank,'] :',i,x(i,:),res(i,:)
!      enddo
c
!      do irank = 0,numpe-1
!        call MPI_Barrier (MPI_COMM_WORLD,ierr)
!        if (irank == myrank) then
!          numtask = ilwork(1)
!          itkbeg = 1
!          m = 0
!          write(*,990) myrank,numtask
!          do itask = 1, numtask
!            m = m + 1
!            iother = ilwork (itkbeg + 3)
!            numseg = ilwork (itkbeg + 4)
!            write(*,991) myrank,ilwork(itkbeg+1:itkbeg+5)
!      if (myrank == 0 .and. iother == 1 .or.
!     &    myrank == 1 .and. iother == 0) then
!        do is = 1,numseg
!          isgbeg = ilwork(itkbeg + 3 + 2*is)
!          lenseg = ilwork(itkbeg + 4 + 2*is)
!          isgend = isgbeg + lenseg - 1
!          do isg = isgbeg,isgend
!            write(*,801) myrank,is,isg,lenseg,x(isg,:)!,res(isg,:)
!          enddo
!        enddo
!      endif
!            itkbeg = itkbeg + 4 + 2*numseg
!          enddo
!        endif
!      enddo
c
801   format('[',i2,'] is,isgbeg,lenseg,isg,x:'3i4,x,3f8.3,x,5e24.16)
990   format('[',i2,'] numtask:',i3)
991   format('[',i2,'] itag, iacc, iother, numseg, isgbeg:',5i6)
998   format(a,i2,a)
999   format(a,i2,a,i6,3f7.3,5e24.16)
c
c--------> END DEBUG <--------------
c
c.... ---------------------->   post processing  <----------------------
c
c.... satisfy the BCs on the residual
c
      call bc3Res (y,  iBC,  BC,  res,  iper, ilwork)
c
c------> BEGIN DEBUG <-----------
c
!      write(*,998) '[',myrank,'] in elmgmr AFTER bc3res.'
!      do i = 1,nshg
!        if (x(i,1) < 0.1001 .and. x(i,1) > 0.0999 
!     & .and. x(i,2) < 0.0901 .and. x(i,2) > 0.0899
!     & .and. x(i,3) < 0.0101 .and. x(i,3) > 0.0099
!     &)
!     &   write(*,999) '[',myrank,'] :',i,x(i,:),res(i,:)
!      enddo
c------> END DEBUG <-----------
c
c
c.... satisfy the BCs on the block-diagonal preconditioner
c
c  This code fragment would extract the "on processor diagonal
c      blocks". lhsK alread has the BC's applied to it (using BC3lhs), 
c      though it was on an ebe basis. For now, we forgo this and still 
c      form BDiag before BC3lhs, leaving the need to still apply BC's
c      below.  Keep in mind that if we used the code fragment below we
c      would still need to make BDiag periodic since BC3lhs did not do
c      that part.
c
      if (iprec .ne. 0) then
c$$$         do iaa=1,nshg
c$$$            k = sparseloc( row(col(iaa)), col(iaa+1)-colm(iaa), iaa )
c$$$     &       + colm(iaa)-1
c$$$            do idof=1,nflow
c$$$               do jdof=1,nflow
c$$$                  idx=idof+(jdof-1)*nflow
c$$$                  BDiag(iaa,idof,jdof)=lhsK(idx,k)
c$$$               enddo
c$$$            enddo
c$$$         enddo
         call bc3BDg (y,  iBC,  BC,  BDiag, iper, ilwork)
      endif
c
c.... return
c
      call timer ('Back    ')
      return
      end
c
c

c
        subroutine ElmGMRSclr(y,      ac,
     &                        x,      elDw,         
     &                        shp,    shgl,   iBC,
     &                        BC,     shpb,   shglb,
     &                        rest,   rmest,  Diag,
     &                        iper,   ilwork, EGmasst)
c
c----------------------------------------------------------------------
c
c This routine computes the LHS mass matrix, the RHS residual 
c vector, and the preconditioning matrix, for use with the GMRES
c solver.
c
c Zdenek Johan, Winter 1991.      (Fortran 90)
c Chris Whiting, Winter 1998.     (Matrix EBE-GMRES)
c----------------------------------------------------------------------
c
        use pointer_data
c
        include "common.h"
        include "mpif.h"
c
        dimension y(nshg,ndof),         ac(nshg,ndof),
     &            x(numnp,nsd),              
     &            iBC(nshg),
     &            BC(nshg,ndofBC),      
     &            rest(nshg),           Diag(nshg),
     &            rmest(nshg),          BDiag(nshg,nflow,nflow),
     &            iper(nshg),           EGmasst(numel,nshape,nshape)
c
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
        dimension qrest(nshg),          rmasst(nshg)
c
        dimension ilwork(nlwork)
        dimension elDw(numel)
c
        real*8, allocatable :: tmpshp(:,:), tmpshgl(:,:,:)
        real*8, allocatable :: tmpshpb(:,:), tmpshglb(:,:,:)
        real*8, allocatable :: elDwl(:)
c
	ttim(80) = ttim(80) - tmr()
c
c.... set up the timer
c

        call timer ('Elm_Form')
c
c.... -------------------->   interior elements   <--------------------
c
c.... set up parameters
c
        intrul = intg  (1,itseq)
        intind = intpt (intrul)
c
        ires   = 1

c       if (idiff>=1) then  ! global reconstruction of qdiff
c
c loop over element blocks for the global reconstruction
c of the diffusive flux vector, q, and lumped mass matrix, rmass
c
c        qrest = zero
c        rmasst = zero
c        
c        do iblk = 1, nelblk
c
c.... set up the parameters
c
c          iel    = lcblk(1,iblk)
c          lelCat = lcblk(2,iblk)
c          lcsyst = lcblk(3,iblk)
c          iorder = lcblk(4,iblk)
c          nenl   = lcblk(5,iblk)   ! no. of vertices per element
c          mattyp = lcblk(7,iblk)
c          ndofl  = lcblk(8,iblk)
c          nsymdl = lcblk(9,iblk)
c          npro   = lcblk(1,iblk+1) - iel 
c
c          nintg  = numQpt (nsd,intrul,lcsyst)
c          
c
c.... compute and assemble diffusive flux vector residual, qres,
c     and lumped mass matrix, rmass
c
c          call AsIq (y,                x,                       
c     &               shp(1,intind,lelCat), 
c     &               shgl(1,intind,lelCat), 
c     &               mien(iblk)%p,    mxmudmi(iblk)%p,  
c     &               qres,           rmass )
c       
c       enddo
c       
c
c.... compute qi for node A, i.e., qres <-- qres/rmass
c
c       if (numpe > 1) then
c          call commu (qres  , ilwork, (ndof-1)*nsd  , 'in ')
c          call commu (rmass , ilwork,  1            , 'in ')
c       endif
c
c.... take care of periodic boundary conditions
c
c       call qpbc( rmass, qres, iBC, iper )       
c       
c       rmass = one/rmass
c       
c       do i=1, (nflow-1)*nsd
c          qres(:,i) = rmass*qres(:,i)
c       enddo
c
c       if(numpe > 1) then
c          call commu (qres, ilwork, (nflow-1)*nsd, 'out')    
c       endif
c
c      endif                     ! computation of global diffusive flux
c
c.... loop over element blocks to compute element residuals
c
c
c.... initialize the arrays
c
        rest    = zero
        rmest   = zero ! to avoid trap_uninitialized
        if (lhs .eq. 1)   EGmasst = zero
        if (iprec. ne. 0)   Diag  = zero 
c
c.... loop over the element-blocks
c
        do iblk = 1, nelblk
c
c
          nenl   = lcblk(5,iblk) ! no. of vertices per element
          iel    = lcblk(1,iblk)
          lelCat = lcblk(2,iblk)
          lcsyst = lcblk(3,iblk)
          iorder = lcblk(4,iblk)
          nenl   = lcblk(5,iblk)   ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          mattyp = lcblk(7,iblk)
          ndofl  = lcblk(8,iblk)
          nsymdl = lcblk(9,iblk)
          npro   = lcblk(1,iblk+1) - iel 
          inum   = iel + npro - 1
          ngauss = nint(lcsyst)
c
c.... compute and assemble the residual and tangent matrix
c
          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))

          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)

          allocate (elDwl(npro))
c
          call AsIGMRSclr(y,                   
     &                    ac,
     &                    x,               elDwl,                   
     &                    tmpshp,          tmpshgl,
     &                    mien(iblk)%p,
     &                    mmat(iblk)%p,    rest,
     &                    rmest,               
     &                    qrest,           EGmasst(iel:inum,:,:),
     &                    Diag )
c
 
c.... satisfy the BC's on the implicit LHS
c     
          call bc3LHSSclr (iBC, mien(iblk)%p, EGmasst(iel:inum,:,:) )
c
          elDw(iel:inum)=elDwl(1:npro)
          deallocate ( elDwl )
          deallocate ( tmpshp )
          deallocate ( tmpshgl )
c.... end of interior element loop
c
       enddo
c
c.... -------------------->   boundary elements   <--------------------
c
c.... set up parameters
c
        intrul = intg   (2,itseq)
        intind = intptb (intrul)
c
c.... loop over the boundary elements
c
        do iblk = 1, nelblb
c
c.... set up the parameters
c
          iel    = lcblkb(1,iblk)
          lelCat = lcblkb(2,iblk)
          lcsyst = lcblkb(3,iblk)
          iorder = lcblkb(4,iblk)
          nenl   = lcblkb(5,iblk)  ! no. of vertices per element
          nenbl  = lcblkb(6,iblk)  ! no. of vertices per bdry. face
          mattyp = lcblkb(7,iblk)
          ndofl  = lcblkb(8,iblk)
          nshl   = lcblkb(9,iblk)
          nshlb  = lcblkb(10,iblk)
          npro   = lcblkb(1,iblk+1) - iel 
          if(lcsyst.eq.3) lcsyst=nenbl
          ngaussb = nintb(lcsyst)          
c
c.... compute and assemble the residuals corresponding to the 
c     boundary integral
c

          allocate (tmpshpb(nshl,MAXQPT))
          allocate (tmpshglb(nsd,nshl,MAXQPT))
          
          tmpshpb(1:nshl,:) = shpb(lcsyst,1:nshl,:)
          tmpshglb(:,1:nshl,:) = shglb(lcsyst,:,1:nshl,:)
c
          call AsBMFGSclr (y,                  x,
     &                     tmpshpb,
     &                     tmpshglb, 
     &                     mienb(iblk)%p,      mmatb(iblk)%p,
     &                     miBCB(iblk)%p,      mBCB(iblk)%p,
     &                     rest,               rmest)
c
          deallocate ( tmpshpb )
          deallocate ( tmpshglb )

c.... end of boundary element loop
c
        enddo


      ttim(80) = ttim(80) + tmr()
c
c.... -------------------->   communications <-------------------------
c

      if (numpe > 1) then
        call commu (rest  , ilwork, 1  , 'in ')

        call MPI_BARRIER (MPI_COMM_WORLD,ierr)

        if(iprec .ne. 0) call commu (Diag, ilwork, 1, 'in ')
      endif

c
c.... ---------------------->   post processing  <----------------------
c
c.... satisfy the BCs on the residual
c
      call bc3ResSclr (y,  iBC,  BC,  rest,  iper, ilwork)
c
c.... satisfy the BCs on the preconditioner
c
      call bc3BDgSclr (iBC, Diag, iper, ilwork)
c
c.... return
c
      call timer ('Back    ')
      return
      end

