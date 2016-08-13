      subroutine genbkif (ibksz)
c
        use pointer_data
        use phio
        use iso_c_binding
        use mattype_m
c
        include "common.h"
        include "mpif.h" !Required to determine the max for itpblk
c
        integer, intent(in) :: ibksz
c
        integer :: iel, iblk, itpblk,
     &             n,neltp,ioprdl,nnface,
     &             iientpsiz0,iientpsiz1,
     &             nsymdl,ndofl
        integer, target :: intfromfile(50)
        integer, target, allocatable :: ientp(:,:)
        integer, dimension(ibksz) :: mater
        integer :: descriptor, descriptorG, GPID, color
        integer ::  numparts, writeLock
        integer :: ierr_io, numprocs
        integer, target :: itpblktot,ierr
        integer, parameter :: ione = 1, itwo = 2, iseven = 7, inine = 9
        character*255 :: fname1, fname2, tempchar
        character(len=30) :: dataInt
        dataInt = c_char_'integer'//c_null_char
c
c
        iel = 1              ! element index
        itpblk = nelblif     ! number of topological blocks
        nelblif = 0         ! set a counter for the interface element blocks (size of ibksz)
c
        ! Get the total number of different interface topologies in the whole domain. 
        ! Try to read from a field. If the field does not exist, scan the geombc file.
        itpblktot=1  ! hardwired to montopology for now
        call phio_readheader(fhandle,
     &   c_char_'total number of interface tpblocks' // char(0),
     &   c_loc(itpblktot), ione, dataInt, iotype) 
c
        if (itpblktot == -1) then 
          ! The field 'total number of different interface tpblocks' was not found in the geombc file.
          ! Scan all the geombc file for the 'connectivity interface' fields to get this information.
          iblk=0
          neltp=0
          do while(neltp .ne. -1) 

            ! intfromfile is reinitialized to -1 every time.
            ! If connectivity interior@xxx is not found, then 
            ! readheader will return intfromfile unchanged

            intfromfile(:)=-1
            iblk = iblk+1
            if(input_mode.ge.1) then
              write (fname2,"('connectivity interface',i1)") iblk
            else
              write (fname2,"('connectivity interface linear tetrahedron tetrahedron')") 
!              write (fname2,"('connectivity interface?')") 
            endif

            !write(*,*) 'rank, fname2',myrank, trim(adjustl(fname2))
            call phio_readheader(fhandle, fname2 // char(0),
     &       c_loc(intfromfile), iseven, dataInt, iotype)
            neltp = intfromfile(1) ! -1 if fname2 was not found, >=0 otherwise
          end do
          itpblktot = iblk-1   
        end if

        if (myrank == 0) then
          write(*,*) 'Number of interface topologies: ',itpblktot
        endif

        mattyp0 = 0
        mattyp1 = 1
        ndofl  = ndof
        nsymdl = nsymdf      ! ????
c
        do iblk = 1, itpblk
c
           writeLock=0;
           if(input_mode.ge.1) then
             write (fname2,"('connectivity interface',i1)") iblk
           else
             write (fname2,"('connectivity interface linear tetrahedron tetrahedron')") 
!              write (fname2,"('connectivity interior?')") 
           endif

           ! Synchronization for performance monitoring, as some parts do not include some topologies
c           call MPI_Barrier(MPI_COMM_WORLD,ierr) 
           call phio_readheader(fhandle, fname2 // char(0),
     &      c_loc(intfromfile), inine, dataInt, iotype)
c
           neltp  = intfromfile(1)       ! number of pair of elements in this block
           nenl0  = intfromfile(2)       ! number of nodes in element 0
           nenl1  = intfromfile(3)       ! number of nodes in element 1
           ipordl = intfromfile(4)       ! polynomial order
           nshl0  = intfromfile(5)
           nshl1  = intfromfile(6)
           nnface = intfromfile(7)       ! number of nodes on the interface
           lcsyst0= intfromfile(8)       ! element type 0
           lcsyst1= intfromfile(9)       ! element type 1

           if (neltp==0) then
              writeLock=1;
           endif
c
c... reads all the connectivity data in one array
c
           iientpsiz = neltp*(nshl0+nshl1)
           allocate (ientp(neltp,(nshl0+nshl1)))
           call phio_readdatablock(fhandle,fname2 // char(0),
     &      c_loc(ientp), iientpsiz, dataInt, iotype)
c
          if (nummat <= 1) then
            call error ('genbkif  ', 'Number of Materials', nummat)
          endif
c
           if(input_mode.ge.1) then
             write(fname2,"('material type interface',i1)") iblk
           else
             write(fname2,"('material type interface linear tetrahedron tetrahedron')")
           endif
c
c           call MPI_Barrier(MPI_COMM_WORLD,ierr) 
           call phio_readheader(fhandle, fname2 // char(0),
     &      c_loc(intfromfile), itwo, dataInt, iotype)
           allocate(mattypeif(intfromfile(1),intfromfile(2)))
           imattypesiz = intfromfile(1)*intfromfile(2)
           call phio_readdatablock(fhandle,fname2 // char(0),
     &      c_loc(mattypeif), imattypesiz, dataInt, iotype)
c
          allocate(ienif0tmp(neltp,nshl0))
          allocate(ienif1tmp(neltp,nshl1))
c
c... set material types:
c
          do imattype = 1,nummat
            ! check the material type from file against the input tags...
            if (mattypeif(1,1) == mat_tag(imattype,1)) then
              mattype0 = imattype
            endif
            if (mattypeif(1,2) == mat_tag(imattype,1)) then
              mattype1 = imattype
            endif
          enddo
c
c          if (mattype0 == mattype1) then
c            call error('genbkif ','Wrong interface mattype in geombc',0)
c          endif
c
          if(writeLock==0) then
c
c... make blocks of elements
c
          blocks_loop: do n = 1, neltp, ibksz
c
            ipro = 0
            npro = 0
            do
              npro = npro + 1
              if (mattypeif(iel+ipro,1) == mat_tag(mattype0,1)) then
                ienif0tmp(npro,1:nshl0) = ientp(iel+ipro,1:nshl0)
                ienif1tmp(npro,1:nshl1) = ientp(iel+ipro,nshl0+1:nshl0+nshl1)
              else if (mattypeif(iel+ipro,2) == mat_tag(mattype0,1)) then
                ienif0tmp(npro,1:nshl0) = ientp(iel+ipro,nshl0+1:nshl0+nshl1)
                ienif1tmp(npro,1:nshl1) = ientp(iel+ipro,1:nshl0)
              else
                call error ('genbkif  ', 'wrong material type', mattype(i))
              endif
c
              ipro = ipro + 1
              if (npro == ibksz .or. iel+ipro>neltp) exit
            enddo
c
            nelblif = nelblif + 1
c
            lcblkif(1,nelblif) = iel
c            lcblkif(2,nelblif) = iopen   ! ??? see genblk.f
            lcblkif(3,nelblif) = lcsyst0  ! local coordinate system
            lcblkif(4,nelblif) = lcsyst1
            lcblkif(5,nelblif) = ipordl
            lcblkif(6,nelblif) = nenl0
            lcblkif(7,nelblif) = nenl1
            lcblkif(8,nelblif) = nnface
            lcblkif(9,nelblif) = mattype0
            lcblkif(10,nelblif) = mattype1
            lcblkif(11,nelblif) = ndof
            lcblkif(12,nelblif) = nsymdl ! ????
            lcblkif(13,nelblif) = nshl0
            lcblkif(14,nelblif) = nshl1
c
c... fill the toplogical id
c    NOTE: the iftpid array can be replaced by using a function which will return
c          topological id 1-12, by taking the two element coordinate systems...
c          The current approach requires additional array of sizeof(int)*MAXBLK...
c
          if (lcsyst0 .eq. 1 .AND. lcsyst1 .eq. 1) then
            iftpid(nelblif) = 1
          else
            write(*,*) 'ERROR: THE INTERFACE TOPOLOGICAL ID HAS NOT COMPLETED!!!'
          endif          
c
c... allocate memory for stack arrays
c
            allocate(mienif0(nelblif)%p(npro,nshl0))
            allocate(mienif1(nelblif)%p(npro,nshl1))
c
c... save the interface elements block
c
            call gensavif(mienif0(nelblif)%p,mienif1(nelblif)%p)
c
            iel = iel + npro
c
          enddo blocks_loop
          endif
c
          deallocate(ientp)
          deallocate(mattypeif)
          deallocate(ienif0tmp,ienif1tmp)
c
        enddo
c
        lcblkif(1,nelblif+1) = iel
c
        return
c
      end subroutine genbkif
