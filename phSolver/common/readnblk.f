c  readnblk.f (pronounce "Reed and Block Dot Eff") contains:
c
c    module readarrays ("Red Arrays") -- contains the arrays that
c     are read in from binary files but not immediately blocked 
c     through pointers.
c
c    subroutine readnblk ("Reed and Block") -- allocates space for
c     and reads data to be contained in module readarrays.  Reads
c     all remaining data and blocks them with pointers.
c
      module readarrays
      
      real*8, allocatable :: point2x(:,:)
      real*8, allocatable :: qold(:,:)
      real*8, allocatable :: uold(:,:)
      real*8, allocatable :: acold(:,:)
      real*8, allocatable :: xn(:,:)
      real*8, allocatable :: xdotold(:,:)
      real*8, allocatable :: umesh(:,:)
      integer, allocatable :: iBCtmp(:)
      real*8, allocatable :: BCinp(:,:)

      integer, allocatable :: point2ilwork(:)
!      integer, allocatable :: fncorp(:)
      integer, allocatable :: twodncorp(:,:)
      integer, allocatable :: nBC(:)
      integer, allocatable :: point2iper(:)
      integer, target, allocatable :: point2ifath(:)
      integer, target, allocatable :: point2nsons(:)
      
      end module

      subroutine readnblk
      use iso_c_binding 
      use readarrays
      use fncorpmod
      use phio
      use phstr
      use syncio
      use posixio
      use streamio
      include "common.h"

      real*8, target, allocatable :: xread(:,:), qread(:,:), acread(:,:)
      real*8, target, allocatable :: uread(:,:), xnread(:,:)
      real*8, target, allocatable :: xdotoldread(:,:), umeshread(:,:)
      real*8, target, allocatable :: BCinpread(:,:)
      real*8 :: iotime
      integer, target, allocatable :: iperread(:), iBCtmpread(:)
      integer, target, allocatable :: ilworkread(:), nBCread(:)
      integer, target, allocatable :: fncorpread(:)
      integer fncorpsize
      character*10 cname2, cname2nd
      character*8 mach2
      character*30 fmt1
      character*255 fname1,fnamer,fnamelr
      character*255 warning
      integer igeomBAK, ibndc, irstin, ierr
      integer, target :: intfromfile(50) ! integers read from headers
      integer :: descriptor, descriptorG, GPID, color, nfields
      integer ::  numparts, nppf
      integer :: ierr_io, numprocs, itmp, itmp2
      integer :: ignored
      integer :: fileFmt
      character*255 fname2, temp2
      character*64 temp1
      type(c_ptr) :: handle
      character(len=1024) :: dataInt, dataDbl
      integer, target, allocatable :: itemp(:)
      dataInt = c_char_'integer'//c_null_char
      dataDbl = c_char_'double'//c_null_char
c
c.... determine the step number to start with
c
      open(unit=72,file='numstart.dat',status='old')
      read(72,*) irstart
      close(72)
      lstep=irstart ! in case restart files have no fields

      fname1='geombc.dat'
      fname1= trim(fname1)  // cname2(myrank+1)

      itmp=1
      if (irstart .gt. 0) itmp = int(log10(float(irstart)))+1
      write (fmt1,"('(''restart.'',i',i1,',1x)')") itmp
      write (fnamer,fmt1) irstart
      fnamer = trim(fnamer) // cname2(myrank+1)
c
c.... open input files
c.... input the geometry parameters
c
      numparts = numpe !This is the common settings. Beware if you try to compute several parts per process

      itwo=2
      ione=1
      ieleven=11
      itmp = int(log10(float(myrank+1)))+1

      iotime = TMRC()
      if( input_mode .eq. -1 ) then
        call streamio_setup_read(fhandle, geomRestartStream)
      else if( input_mode .eq. 0 ) then
        call posixio_setup(fhandle, c_char_'r')
      else if( input_mode .ge. 1 ) then
        call syncio_setup_read(nsynciofiles, fhandle)
      end if
      call phio_constructName(fhandle, 
     &        c_char_'geombc' // char(0), fname1)
      call phio_openfile(fname1, fhandle);

      call phio_readheader(fhandle,c_char_'number of nodes' // char(0),
     & c_loc(numnp),ione, dataInt, iotype)

      call phio_readheader(fhandle,c_char_'number of modes' // char(0),
     & c_loc(nshg),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     &  c_char_'number of interior elements' // char(0),
     &  c_loc(numel),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     &  c_char_'number of boundary elements' // char(0),
     &  c_loc(numelb),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     &  c_char_'number of interface elements' // char(0),
     &  c_loc(numelif),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     &  c_char_'maximum number of element nodes' // char(0),
     &  c_loc(nen),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     &  c_char_'number of interior tpblocks' // char(0),
     &  c_loc(nelblk),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     & c_char_'number of boundary tpblocks' // char(0),
     & c_loc(nelblb),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     & c_char_'number of interface tpblocks' // char(0),
     & c_loc(nelblif),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     & c_char_'number of nodes with Dirichlet BCs' // char(0),
     & c_loc(numpbc),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     & c_char_'number of shape functions' // char(0),
     & c_loc(ntopsh),ione, dataInt, iotype)

      write(fname2,"('material type interior')")
      call phio_readheader(fhandle, fname2 // char(0),
     & c_loc(intfromfile(1)), 1, dataInt, iotype)
      allocate(mattype_interior(intfromfile(1)))
      allocate(itemp(intfromfile(1)))
      call phio_readdatablock(fhandle,fname2 // char(0),
     & c_loc(itemp), intfromfile(1), dataInt, iotype)
      mattype_interior(:) = itemp(:)
      deallocate(itemp)
c
c.... calculate the maximum number of boundary element nodes
c     
      nenb = 0
      do i = 1, melCat
         if (nen .eq. nenCat(i,nsd)) nenb = max(nenCat(i,nsd-1), nenb)
      enddo
      if (myrank == master) then
         if (nenb .eq. 0) call error ('input   ','nen     ',nen)
      endif
c
c.... setup some useful constants
c
      I3nsd  = nsd / 3          ! nsd=3 integer flag
      E3nsd  = float(I3nsd)     ! nsd=3 real    flag
      if(matflg(1,1).lt.0) then
         nflow = nsd + 1
      else
         nflow = nsd + 2
      endif 
c
      ndof   = nsd + 2
      nsclr=impl(1)/100
      ndof=ndof+nsclr           ! number of sclr transport equations to solve
c      
      if (iALE .eq. 2) then     ! Mesh-elastic is ON
         nelas   = nsd              ! FOR mesh-elastic 
         ndofBC  = ndof + I3nsd     ! dimension of BC array
     &          + nelas + I3nsd    ! add nelas for mesh-elastic solve
         ndofBC2 = 3+2+4+7+8    ! (assuming 4 scalars to be ON) and 8 is for ec11 ec12 ec13 em1 ec21 ec22 ec23 em2 
      else
         nelas   = 0
         ndofBC  = ndof + I3nsd     ! dimension of BC array
         ndofBC2 = ndof + 7 
      endif
c
      ndiBCB = 2                ! dimension of iBCB array
      ndBCB  = ndof + 1         ! dimension of BCB array
      nsymdf = (ndof*(ndof + 1)) / 2 ! symm. d.o.f.'s
c
c.... ----------------------> Communication tasks <--------------------
c
      if(numpe > 1) then
         call phio_readheader(fhandle,
     &    c_char_'size of ilwork array' // char(0),
     &    c_loc(nlwork),ione, dataInt, iotype)

         call phio_readheader(fhandle,
     &    c_char_'ilwork' //char(0),
     &    c_loc(nlwork),ione, dataInt, iotype)

         allocate( point2ilwork(nlwork) )
         allocate( ilworkread(nlwork) )
         call phio_readdatablock(fhandle, c_char_'ilwork' // char(0),
     &      c_loc(ilworkread), nlwork, dataInt , iotype)

         point2ilwork = ilworkread
         call ctypes (point2ilwork)

       if((usingPETSc.eq.1).or.(svLSFlag.eq.1)) then
         fncorpsize = nshg
         allocate(fncorp(fncorpsize))
         call gen_ncorp(fncorp, ilworkread, nlwork, fncorpsize)
!
! the  following code finds the global range of the owned nodes
!
         maxowned=0
         minowned=maxval(fncorp)
         do i = 1,nshg      
            if(fncorp(i).gt.0) then  ! don't consider remote copies
               maxowned=max(maxowned,fncorp(i))
               minowned=min(minowned,fncorp(i))
            endif
         enddo
!
!  end of global range code
!
         call commuInt(fncorp, point2ilwork, 1, 'out')
         ncorpsize = fncorpsize 
       endif
       if(svLSFlag.eq.1) then
         allocate(ltg(ncorpsize))
         ltg(:)=fncorp(:)
       endif
      else
           allocate(ltg(nshg))
           do i =1,nshg
             ltg(i)=i
           enddo
           nlwork=1
           allocate( point2ilwork(1))
           nshg0 = nshg
      endif


      itwo=2

      call phio_readheader(fhandle,
     & c_char_'co-ordinates' // char(0),
     & c_loc(intfromfile),itwo, dataDbl, iotype)
      numnp=intfromfile(1)
      allocate( point2x(numnp,nsd) )
      allocate( xread(numnp,nsd) )
      ixsiz=numnp*nsd
      call phio_readdatablock(fhandle,
     & c_char_'co-ordinates' // char(0),
     & c_loc(xread),ixsiz, dataDbl, iotype)
      point2x = xread
c
c.... read in and block out the connectivity
c
      call genblk (IBKSIZ)
c
c.... read the boundary condition mapping array
c
      ione=1
      call phio_readheader(fhandle,
     & c_char_'bc mapping array' // char(0),
     & c_loc(nshg),ione, dataInt, iotype)

      allocate( nBC(nshg) )

      allocate( nBCread(nshg) )

      call phio_readdatablock(fhandle,
     & c_char_'bc mapping array' // char(0),
     & c_loc(nBCread), nshg, dataInt, iotype)

      nBC=nBCread
c
c.... read the temporary iBC array
c
      ione=1
      call phio_readheader(fhandle,
     & c_char_'bc codes array' // char(0),
     & c_loc(numpbc),ione, dataInt, iotype)

      if ( numpbc > 0 ) then
        allocate( iBCtmp(numpbc) )
        allocate( iBCtmpread(numpbc) )
      else
        allocate( iBCtmp(1) )
        allocate( iBCtmpread(1) )
      endif
      call phio_readdatablock(fhandle,
     & c_char_'bc codes array' // char(0),
     & c_loc(iBCtmpread), numpbc, dataInt, iotype)

      if ( numpbc > 0 ) then
         iBCtmp=iBCtmpread
      else  ! sometimes a partition has no BC's
         deallocate( iBCtmpread)
         iBCtmp=0
      endif
c
c.... read boundary condition data
c
      ione=1
      call phio_readheader(fhandle,
     & c_char_'boundary condition array' // char(0),
     & c_loc(intfromfile),ione, dataDbl, iotype)

      if ( numpbc > 0 ) then
c         if(intfromfile(1).ne.(ndof+7)*numpbc) then
c           warning='WARNING more data in BCinp than needed: keeping 1st'
c           write(*,*) warning, ndof+7, numpbc,intfromfile(1),(ndof+7)*numpbc
c         endif
         allocate( BCinp(numpbc,ndofBC2) )
         nsecondrank=intfromfile(1)/numpbc
         allocate( BCinpread(numpbc,nsecondrank) )
         iBCinpsiz=intfromfile(1)
      else
         allocate( BCinp(1,ndofBC2) )
         allocate( BCinpread(0,0) ) !dummy
         iBCinpsiz=intfromfile(1)
      endif

      call phio_readdatablock(fhandle,
     & c_char_'boundary condition array' // char(0),
     & c_loc(BCinpread), iBCinpsiz, dataDbl, iotype)

      if ( numpbc > 0 ) then
         BCinp(:,1:ndofBC2)=BCinpread(:,1:ndofBC2)
      else  ! sometimes a partition has no BC's
         deallocate(BCinpread)
         BCinp=0
      endif
c
c.... read periodic boundary conditions
c
      ione=1
      call phio_readheader(fhandle,
     & c_char_'periodic masters array' // char(0),
     & c_loc(nshg), ione, dataInt, iotype)
      allocate( point2iper(nshg) )
      allocate( iperread(nshg) )
      call phio_readdatablock(fhandle,
     & c_char_'periodic masters array' // char(0),
     & c_loc(iperread), nshg, dataInt, iotype)
      point2iper=iperread
c
c.... generate the boundary element blocks
c
      call genbkb (ibksiz)
c
      call genbkif (ibksiz)
c
c  Read in the nsons and ifath arrays if needed
c
c  There is a fundamental shift in the meaning of ifath based on whether
c  there exist homogenous directions in the flow.  
c
c  HOMOGENOUS DIRECTIONS EXIST:  Here nfath is the number of inhomogenous
c  points in the TOTAL mesh.  That is to say that each partition keeps a 
c  link to  ALL inhomogenous points.  This link is furthermore not to the
c  sms numbering but to the original structured grid numbering.  These 
c  inhomogenous points are thought of as fathers, with their sons being all
c  the points in the homogenous directions that have this father's 
c  inhomogeneity.  The array ifath takes as an arguement the sms numbering
c  and returns as a result the father.
c
c  In this case nsons is the number of sons that each father has and ifath
c  is an array which tells the 
c
c  NO HOMOGENOUS DIRECTIONS.  In this case the mesh would grow to rapidly
c  if we followed the above strategy since every partition would index its
c  points to the ENTIRE mesh.  Furthermore, there would never be a need
c  to average to a node off processor since there is no spatial averaging.
c  Therefore, to properly account for this case we must recognize it and
c  inerrupt certain actions (i.e. assembly of the average across partitions).
c  This case is easily identified by noting that maxval(nsons) =1 (i.e. no
c  father has any sons).  Reiterating to be clear, in this case ifath does
c  not point to a global numbering but instead just points to itself.
c
      nfath=1  ! some architectures choke on a zero or undeclared
                 ! dimension variable.  This sets it to a safe, small value.
      if(((iLES .lt. 20) .and. (iLES.gt.0))
     &                   .or. (itwmod.gt.0)  ) then ! don't forget same
                                                    ! conditional in proces.f
                                                    ! needed for alloc
         ione=1
         if(nohomog.gt.0) then
            call phio_readheader(fhandle,
     &       c_char_'number of father-nodes' // char(0),
     &       c_loc(nfath), ione, dataInt, iotype)

            call phio_readheader(fhandle,
     &       c_char_'number of son-nodes for each father' // char(0),
     &       c_loc(nfath), ione, dataInt, iotype)

            allocate (point2nsons(nfath))

            call phio_readdatablock(fhandle,
     &       c_char_'number of son-nodes for each father' // char(0),
     &       c_loc(point2nsons),nfath, dataInt, iotype)

            call phio_readheader(fhandle,
     &       c_char_'keyword ifath' // char(0),
     &       c_loc(nshg), ione, dataInt, iotype);

            allocate (point2ifath(nshg))

            call phio_readdatablock(fhandle,
     &       c_char_'keyword ifath' // char(0),
     &       c_loc(point2ifath), nshg, dataInt, iotype)
     
            nsonmax=maxval(point2nsons)
         else  ! this is the case where there is no homogeneity
               ! therefore ever node is a father (too itself).  sonfath
               ! (a routine in NSpre) will set this up but this gives
               ! you an option to avoid that.
            nfath=nshg
            allocate (point2nsons(nfath))
            point2nsons=1
            allocate (point2ifath(nshg))
            do i=1,nshg
               point2ifath(i)=i
            enddo
            nsonmax=1
         endif
      else
         allocate (point2nsons(1))
         allocate (point2ifath(1))
      endif

      call phio_closefile(fhandle);
      iotime = TMRC() - iotime
      if (myrank.eq.master) then
        write(*,*) 'time to read geombc (seconds)', iotime
      endif

c.... Read restart files
      iotime = TMRC()
      if( input_mode .eq. -1 ) then
        call streamio_setup_read(fhandle, geomRestartStream)
      else if( input_mode .eq. 0 ) then
        call posixio_setup(fhandle, c_char_'r')
      else if( input_mode .ge. 1 ) then
        call syncio_setup_read(nsynciofiles, fhandle)
      end if
      call phio_constructName(fhandle,
     &        c_char_'restart' // char(0), fnamer)
      call phstr_appendInt(fnamer, irstart)
      call phstr_appendStr(fnamer, c_char_'.'//c_null_char)
      call phio_openfile(fnamer, fhandle);

      ithree=3

      itmp = int(log10(float(myrank+1)))+1

      intfromfile=0
      call phio_readheader(fhandle,
     & c_char_'solution' // char(0), 
     & c_loc(intfromfile), ithree, dataInt, iotype)
c
c.... read the values of primitive variables into q
c
      allocate( qold(nshg,ndof) )
      if(intfromfile(1).ne.0) then
         nshg2=intfromfile(1)
         ndof2=intfromfile(2)
         lstep=intfromfile(3)
         if(ndof2.ne.ndof) then

         endif
        if (nshg2 .ne. nshg)
     &        call error ('restar  ', 'nshg   ', nshg)
         allocate( qread(nshg,ndof2) )
         iqsiz=nshg*ndof2
         call phio_readdatablock(fhandle,
     &    c_char_'solution' // char(0),
     &    c_loc(qread),iqsiz, dataDbl,iotype)
         qold(:,1:ndof)=qread(:,1:ndof)
         deallocate(qread)
      else
         if (myrank.eq.master) then
            if (matflg(1,1).eq.0) then ! compressible
               warning='Solution is set to zero (with p and T to one)'
            else
               warning='Solution is set to zero'
            endif
            write(*,*) warning
         endif
         qold=zero
         if (matflg(1,1).eq.0) then ! compressible
            qold(:,1)=one ! avoid zero pressure
            qold(:,nflow)=one ! avoid zero temperature
         endif
      endif

      intfromfile=0
      call phio_readheader(fhandle,
     & c_char_'time derivative of solution' // char(0),
     & c_loc(intfromfile), ithree, dataInt, iotype)
      allocate( acold(nshg,ndof) )
      if(intfromfile(1).ne.0) then
         nshg2=intfromfile(1)
         ndof2=intfromfile(2)
         lstep=intfromfile(3)

         if (nshg2 .ne. nshg)
     &        call error ('restar  ', 'nshg   ', nshg)
         allocate( acread(nshg,ndof2) )
         acread=zero
         iacsiz=nshg*ndof2
         call phio_readdatablock(fhandle,
     &    c_char_'time derivative of solution' // char(0),
     &    c_loc(acread), iacsiz, dataDbl,iotype)
         acold(:,1:ndof)=acread(:,1:ndof)
         deallocate(acread)
      else
         if (myrank.eq.master) then
            warning='Time derivative of solution is set to zero (SAFE)'
            write(*,*) warning
         endif
         acold=zero
      endif

c read in ALE stuff
c read in coordinate at n time step      
      intfromfile=0
      call phio_readheader(fhandle, 
     & c_char_'motion_coords' //char(0), 
     & c_loc(intfromfile), ithree, dataInt, iotype)
      allocate( xn(numnp,nsd) )
      if(intfromfile(1).ne.0) then 
         numnp2=intfromfile(1)
         nsd2=intfromfile(2)
         lstep=intfromfile(3)
         
         if (numnp2 .ne. numnp) 
     &        call error ('restar  ', 'numnp   ', numnp)
c     
         allocate( xnread(numnp,nsd2) )
         xnread=zero

         iacsiz=numnp*nsd2
         call phio_readdatablock(fhandle,
     &    c_char_'motion_coords' // char(0),
     &    c_loc(xnread), iacsiz, dataDbl,iotype)
         xn(:,1:nsd)=xnread(:,1:nsd)
         deallocate(xnread)
      else
         if (myrank.eq.master) then
           warning='Mesh coordinates set to original coordinates (SAFE)'
            write(*,*) warning
         endif
         xn=point2x
      endif
c read in xdotold
c      fname1='xdot?'
c      intfromfile=0
c      call phio_readheader(fhandle, 
c     & c_char_'xdot' //char(0), 
c     & c_loc(intfromfile), ithree, dataInt, iotype)
c      call readheader(irstin,fname1,intfromfile,
c     &     ithree,'integer', iotype)
c      allocate( xdotold(numnp,nsd) )
c      if(intfromfile(1).ne.0) then 
c         numnp2=intfromfile(1)
c         nsd2=intfromfile(2)
c         lstep=intfromfile(3)
c         
c         if (numnp2 .ne. numnp) 
c     &        call error ('restar  ', 'numnp   ', numnp)
c     
c         allocate( xdotoldread(numnp,nsd2) )
c         xdotoldread=zero
c
c         iacsiz=numnp*nsd2
c
c         call phio_readdatablock(fhandle,
c     &    c_char_'xdot' // char(0),
c     &    c_loc(xdotoldread), iacsiz, dataDbl,iotype)
c         call readdatablock(irstin,fname1,xdotoldread,iacsiz,
c     &                   'double',iotype)
c         xdotold(:,1:nsd)=xdotoldread(:,1:nsd)
c         deallocate(xdotoldread)
c      else
c         if (myrank.eq.master) then
c            warning='Time derivative of mesh disp is set to zero (SAFE)'
c            write(*,*) warning
c         endif
c         xdotold=zero
c      endif
c read in umesh
      fname1='umesh?'
      intfromfile=0
       call phio_readheader(fhandle, 
     & c_char_'mesh_vel' //char(0), 
     & c_loc(intfromfile), ithree, dataInt, iotype)
       allocate( umesh(numnp,nsd) )
       if(intfromfile(1).ne.0) then 
          numnp2=intfromfile(1)
          nsd2=intfromfile(2)
          lstep=intfromfile(3)
          
          if (numnp2 .ne. numnp) 
     &        call error ('restar  ', 'numnp   ', numnp)
      
          allocate( umeshread(numnp,nsd2) )
          umeshread=zero
 
          iacsiz=numnp*nsd2
          call phio_readdatablock(fhandle,
     &    c_char_'mesh_vel' // char(0),
     &    c_loc(umeshread), iacsiz, dataDbl,iotype)
          umesh(:,1:nsd)=umeshread(:,1:nsd)
          deallocate(umeshread)
       else
          if (myrank.eq.master) then
             warning='Mesh velocity is set to zero (SAFE)'
             write(*,*) warning
          endif
         umesh=zero
       endif
c
c.... end read ALE stuff
c
cc
cc.... read the header and check it against the run data
cc
      if (ideformwall.eq.1) then

          intfromfile=0
          call phio_readheader(fhandle,
     &     c_char_'displacement' // char(0),
     &     c_loc(intfromfile), ithree, dataInt, iotype)

         nshg2=intfromfile(1)
         ndisp=intfromfile(2)
         lstep=intfromfile(3)
         if(ndisp.ne.nsd) then
            warning='WARNING ndisp not equal nsd'
            write(*,*) warning , ndisp
         endif
         if (nshg2 .ne. nshg) 
     &        call error ('restar  ', 'nshg   ', nshg)
c
c.... read the values of primitive variables into uold
c
         allocate( uold(nshg,nsd) )
         allocate( uread(nshg,nsd) )
         
         iusiz=nshg*nsd

         call phio_readdatablock(fhandle,
     &    c_char_'displacement' // char(0),
     &    c_loc(uread), iusiz, dataDbl, iotype)

         uold(:,1:nsd)=uread(:,1:nsd)
       else
         allocate( uold(nshg,nsd) )
         uold(:,1:nsd) = zero
       endif
c
c.... close c-binary files
c
      call phio_closefile(fhandle)
      iotime = TMRC() - iotime
      if (myrank.eq.master) then
        write(*,*) 'time to read restart (seconds)', iotime
      endif

      deallocate(xread)
      if ( numpbc > 0 )  then
         deallocate(bcinpread)
         deallocate(ibctmpread)
      endif
      deallocate(iperread)
      if(numpe.gt.1)
     &     deallocate(ilworkread)
      deallocate(nbcread)

      return
 994  call error ('input   ','opening ', igeomBAK)
 995  call error ('input   ','opening ', igeomBAK)
 997  call error ('input   ','end file', igeomBAK)
 998  call error ('input   ','end file', igeomBAK)
      end
c
c No longer called but kept around in case....
c
      subroutine genpzero(iBC)

      use pointer_data
      include "common.h"
      integer iBC(nshg)
c
c....  check to see if any of the nodes have a dirichlet pressure
c
      pzero=1
      if (any(btest(iBC,2))) pzero=0  
      do iblk = 1, nelblb
         npro = lcblkb(1,iblk+1)-lcblkb(1,iblk)
         do i=1, npro
            iBCB1=miBCB(iblk)%p(i,1)
c     
c.... check to see if any of the nodes have a Neumann pressure 
c     but not periodic (note that 
c     
            if(btest(iBCB1,1)) pzero=0
         enddo
c     
c.... share results with other processors
c     
         pzl=pzero
         if (numpe .gt. 1)
     &        call MPI_ALLREDUCE (pzl, pzero, 1,
     &        MPI_DOUBLE_PRECISION,MPI_MIN, MPI_COMM_WORLD,ierr)
      enddo
      return
      end
