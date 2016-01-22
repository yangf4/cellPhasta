        subroutine ElmGMRElas (x,     disp,    shp,       shgl,  
     &                         iBC,   BC,      elasres,   elasBDiag,
     &                         iper,  ilwork,  elaslhsK,  
     &                         col,   row,     meshq)
c
c----------------------------------------------------------------------
c
c This routine computes the LHS mass matrix, the RHS residual 
c vector, and the preconditioning matrix, for mesh-elastic solve
c
c----------------------------------------------------------------------
c
         use pointer_data
c         use timedata
         use timedataC
c
        include "common.h"
        include "mpif.h"
c
        integer col(nshg+1), row(nnz*nshg)
        real*8  elaslhsK(nelas*nelas,nnz_tot),
     &          meshq(numel),
     &          meshV(numel)
c        
        dimension x(numnp,nsd),        disp(numnp,nsd),              
     &            iBC(nshg),
     &            BC(nshg,ndofBC),      
     &            elasres(nshg,nelas),
     &            elasBDiag(nshg,nelas,nelas),
     &            iper(nshg)
c
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT) 
c
        dimension ilwork(nlwork)
c
        real*8, allocatable :: tmpshp(:,:), tmpshgl(:,:,:)
        real*8, allocatable :: Estiff(:,:,:)
c
c.... -------------------->   interior elements   <--------------------
c
c.... loop over element blocks to compute element residuals
c
c
c.... initialize the arrays
c
        elasres = zero
        if (lhs. eq. 1)    elaslhsK  = zero
        if (iprec .ne. 0)  elasBDiag = zero
c
c.... loop over the element-blocks
c
        do iblk = 1, nelblk
c
c.... set up the parameters
c
          iblkts = iblk            ! used in timeseries
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
c     ndofelas = nshl * nelas
c
          allocate (Estiff(npro,ndofelas,ndofelas))
          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))
c
          Estiff = zero
          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)
c
c.... Shape measure. Calculate the shape quality
c
          call shpMeasure(x, mien(iblk)%p, meshq(iel:iel+npro-1),
     &                                     meshV(iel:iel+npro-1) ) 
c
          call AsIGMRElas (x,             disp,
     &                     tmpshp,        tmpshgl,    
     &                     mien(iblk)%p,  elasres,      
     &                     elasBDiag,     Estiff,
     &                     meshq(iel:iel+npro-1),
     &                     meshV(iel:iel+npro-1)   )
c
c.... satisfy the BCs on the implicit LHS
c     
          call bc3LHSElas (iBC, BC(:, ndof+2:ndof+4),
     &                     mien(iblk)%p, Estiff) 
c
c.... Fill-up the global sparse LHS mass matrix
c
          call fillsparseElas( mien(iblk)%p, Estiff,
     &                         elaslhsK, row, col)
c
          deallocate ( Estiff )
          deallocate ( tmpshp )
          deallocate ( tmpshgl )
c
c.... end of interior element loop
c
       enddo
c
      if(iabc==1) then               ! are there any axisym BCs
          call rotabc(elasres(1,2), iBC,  'in ')
      endif
c
c.... -------------------->   communications <-------------------------
c
      if (numpe > 1) then
          call commu (elasres  , ilwork, nelas  , 'in ')
c
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
c
          if (iprec .ne. 0) then
             call commu (elasBDiag,    ilwork,
     &                   nelas*nelas,  'in ')
          endif
      endif
c
c.... ---------------------->   post processing  <----------------------
c
c.... satisfy the BCs on the residual
c
      call bc3ResElas (iBC,     BC(:, ndof+2:ndof+4),
     &                 elasres, iper,    ilwork)
c
c.... satisfy the BCs on the block-diagonal preconditioner
c
      if (iprec .ne. 0) then
         call bc3BDgElas (iBC,       BC(:, ndof+2:ndof+4),
     &                    elasBDiag, iper,    ilwork)
      endif
c
c.... return
c
      call timer ('Back    ')
      return
      end

