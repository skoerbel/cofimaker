      module crysanpero

      implicit none

      contains

c	module constants
c	implicit none
c	double precision Bohr,ec,Pi
c	double precision cornerweight,centerweight,faceweight
c	double precision cornerdist,centerdist,facedist
c	double precision tol_cell
c	parameter (Bohr=0.529177D0,ec=1.60217648740E-19)
c	parameter (Pi=acos(-1.0D0))
c	parameter(cornerweight=0.125D0,centerweight=1.0D0,
c     &            faceweight=0.5D0)
c	parameter(cornerdist=0.6D0,centerdist=0.6D0,
c     &            facedist=0.6D0)
c	parameter(tol_cell=0.4D0)
c	end module constants

c      module functions
c      implicit none
c      contains

      logical function isincell(thisatom,centeratom,nuc)
      use defs
      implicit none
      
      double precision distvec(1:3)
      integer i,nele,iele
c      type patom
c        sequence
c        integer:: number
c        character:: name*2 
c        integer:: direc 
c        double precision:: where(3)
c        double precision:: charge(1:3,1:3)
c        double precision:: weight
c        double precision:: distance
c        character:: geometry*6
c      end type patom
      type(patom),intent(in) :: thisatom,centeratom
      integer, intent(in) :: nuc(1:3)
      do i=1,3
        distvec(i)=thisatom%where(i)-centeratom%where(i)
      end do
      isincell=.False.
      if (abs(distvec(1))*dble(nuc(1)).le.thisatom%distance.and.
     &      abs(distvec(2))*dble(nuc(2)).le.thisatom%distance.and. 
     &      abs(distvec(3))*dble(nuc(3)).le.thisatom%distance) 
     &      isincell=.True.
	end function isincell

	function nearestcubsite(site,geometry,
     &    cubcorners,cubcenters,cubfaces,ncubcorner,ncubcenter,
     &    ncubface)
	implicit none
	double precision nearestcubsite(1:3)
	integer, intent(in) :: ncubcorner,ncubcenter,ncubface
	double precision, intent(in) :: site(1:3),
     &    cubcorners(1:3,1:ncubcorner),cubcenters(1:3,1:ncubcenter),
     &    cubfaces(1:3,1:ncubface)
	character, intent(in) :: geometry*6
	double precision dist0(1:3),dist1(1:3)
	integer i,icart
	! initialize variables
	nearestcubsite=0.0D0
	dist0=2.0D0	
	! if atom is corneratom
	if (geometry.eq.'corner') then
	  do i=1,ncubcorner
	    do icart=1,3
	      dist1(icart)=site(icart)-cubcorners(icart,i)
	    end do
	    if (dist1(1)**2+dist1(2)**2+dist1(3)**2.lt.
     &        dist0(1)**2+dist0(2)**2+dist0(3)**2) then
              dist0=dist1
	      nearestcubsite(1:3)=cubcorners(1:3,i)
	    end if
	  end do		
	end if
	! if atom is centeratom
	if (geometry.eq.'center') then
	  do i=1,ncubcenter
	    do icart=1,3
	      dist1(icart)=site(icart)-cubcenters(icart,i)
	    end do
	    if (dist1(1)**2+dist1(2)**2+dist1(3)**2.lt.
     &        dist0(1)**2+dist0(2)**2+dist0(3)**2) then
              dist0=dist1
	      nearestcubsite(1:3)=cubcenters(1:3,i)
	    end if
	  end do		
	end if
	! if atom is faceatom
	if (geometry.eq.'face  ') then
	  do i=1,ncubface
	    do icart=1,3
	      dist1(icart)=site(icart)-cubfaces(icart,i)
	    end do
	    if (dist1(1)**2+dist1(2)**2+dist1(3)**2.lt.
     &        dist0(1)**2+dist0(2)**2+dist0(3)**2) then
              dist0=dist1
	      nearestcubsite(1:3)=cubfaces(1:3,i)
	    end if
	  end do		
	end if
	end function
c	end module functions

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine scrysanpero()

      use defs
!      use functions

      implicit none

	!--- structure related variables -----------------------------------------------------------------
	integer icorner,ncorner,ncenter,icenter,nface,iface
	integer nuc(1:3),ncubcorner,ncubcenter,
     &          ncubface,irun  !,direction
	integer natoms,nelements,natom
	integer iatom,iatom2,iele,icell,ncell,natomsXL
	integer iatomXL,iatomXL2,natomcell,iele2
	
	type(pelement), allocatable  :: pelements(:)
	
c	type patom
c   	  sequence
c   	  integer:: number
c	  character:: name*2 
c	  integer:: direc 
c   	  double precision:: where(3)
c   	  double precision:: charge(1:3,1:3)
c   	  double precision:: weight
c   	  double precision:: distance
c	  character:: geometry*6
c	end type patom
	type(patom), allocatable :: corners(:),centers(:),faces(:)
	type(patom), allocatable :: atoms(:),atomsref(:)
	type(patom), allocatable :: atomsXL(:),atomsrefXL(:)
	type(patom)  :: thisatom,thisrefatom,centeratom,cubatom

	double precision totalvol,avvol,alat(1:3),angles(1:3),
     &   avecs(1:3,1:3)
	double precision, allocatable :: cubcorners(:,:),
     &                    cubcenters(:,:),cubfaces(:,:)
	double precision x(1:3),distvec(1:3)
	double precision eleweight,elecharge(1:3,1:3),eledistance,
     &                   tol(1:3)
	
	character, allocatable :: centername(:)*2

	!--- polarization --------------------------------------------------------------------------------
	logical polarization
	double precision, allocatable :: dipol(:,:)
	double precision dipoltot(1:3),polfac

	!--- displacements -------------------------------------------------------------------------------
	logical print_disps
	double precision disp(1:3)

	! xsf file
	logical print_xsf

	!--- reference coordinates -----------------------------------------------------------------------
	integer natomsrefXL
	logical reference
	double precision xref(1:3)

	!--- octahedral rotation -------------------------------------------------------------------------
	integer alternate_rot(1:3)
	double precision rotangle(1:3)
	logical rotation

	!--- AFE displacements-- -------------------------------------------------------------------------
	integer alternate_afe(1:3)
	double precision afevec(1:3),afecenteramp,afefaceamp
	logical AFE
	character afemod*3

	!--- input/output --------------------------------------------------------------------------------
	integer inunit,outunit
	integer nlines,iread,ireadend
	logical found
	character name*2, line*80,nameat*2,elename*2

	!--- other variables -----------------------------------------------------------------------------
	integer n,icart,i
	integer n1,n2,n3
	double precision t_start,t_end
        logical isopen12,isopen33
	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call cpu_time(t_start)
      inquire(unit=12,opened=isopen12)
      if(isopen12) then
        write(12,fsubstart)trim(adjustl('crysanpero'))
        write(12,*) ' '
        write(12,*) '     output goes to OUTPUT.CRYSANPERO.'
      else
        print fsubstart, trim(adjustl('crysanpero'))
        print*, ' '
        print*, '     output goes to OUTPUT.CRYSANPERO.'
      end if

       ! the ouput of the program goes here.
       outunit=33 !13
       inquire(unit=33,opened=isopen33)
       if(isopen33) print*, '33 is already open.'
       open(outunit,file='OUTPUT.CRYSANPERO',status='replace')
       write(outunit,'("*****************************************")') 
       write(outunit,*) '* This is CrySAnPero, Version 0.18      *' 
       write(outunit,*) '* (2011-08-30).                         *' 
       write(outunit,*) '*                                       *' 
       write(outunit,*) '* New feature(s):                       *' 
       write(outunit,*) '* ------------------------------------- *' 
       write(outunit,*) '* None.                                 *'
       write(outunit,*) '* ------------------------------------- *' 
       write(outunit,*) '* Modification(s):                      *' 
       write(outunit,*) '* ------------------------------------- *' 
       write(outunit,*) '* coordinates from COORAT.IN are shifted*'
       write(outunit,*) '* between about 0 and 1.                *'
       write(outunit,*) '*****************************************' 
       write(outunit,*) ' ' 
       
       ! The input is read from file "INPUT.CRYSANPERO"
       inunit=10 !34  !14
       open(inunit, file='INPUT.CRYSANPERO',status='old')

       ! read 3 cartesian lattice constants and unit cell numbers from file INPUT.CRYSANPERO
       call findstring(inunit,'box',3,found)
       if (.not.found) then
         print*,'keyword box not found. Quitting.'
         close(outunit)
         stop
       end if
       read(inunit,*) alat(1:3),angles(1:3)
       read(inunit,*) nuc(1:3)
       do icart=1,3
         alat(icart)=alat(icart)*dble(nuc(icart))
         tol(icart)=0.2D0/dble(nuc(icart))
       end do
       write(outunit,*)'cell edges in Bohr:'
       write(outunit,'(3(F10.6))') alat(1:3)
       write(outunit,*)'cell angles in degree:'
       write(outunit,'(3(F10.6))') angles(1:3)
       avecs=0.0D0
       avecs(1,1)=alat(1)
       avecs(2,1)=alat(1)*cos(angles(3)*Pi/180.0D0)
       avecs(2,2)=alat(2)*sin(angles(3)*Pi/180.0D0)
       avecs(3,1)=alat(3)*cos(angles(2)*Pi/180.0D0)
       avecs(3,2)=alat(3)*(cos(angles(1)*Pi/180.0D0)
     &       -cos(angles(3)*Pi/180.0D0)*cos(angles(2)*Pi/180.0D0))
     &       /sin(angles(3)*Pi/180.0D0)
       avecs(3,3)=sqrt(alat(3)**2-avecs(3,1)**2-avecs(3,2)**2)
       write(outunit,*)'cell vectors:'
       write(outunit,'(3(F10.6))') avecs(1,1:3)
       write(outunit,'(3(F10.6))') avecs(2,1:3)
       write(outunit,'(3(F10.6))') avecs(3,1:3)
       write(outunit,*) ' ' 
       
       ! these atoms sit on the cube corners (A sites)
       call findstring(inunit,'corner',6,found)
       if (.not.found) then
         print*,'keyword "corner" not found. Quitting.'
         close(outunit)
         stop
       end if
       read(inunit,*) ncorner
       allocate(corners(ncorner))
       do icorner=1,ncorner
         read(inunit,*) corners(icorner)%name,
     &      corners(icorner)%charge(1,1:3),
     &      corners(icorner)%charge(2,1:3),
     &      corners(icorner)%charge(3,1:3) !,
!     &      corners(icorner)%weight, 
!     &      corners(icorner)%distance 
           corners(icorner)%weight=cornerweight 
           corners(icorner)%distance=cornerdist 
           corners(icorner)%geometry='corner' 
       end do
       write(outunit,*) 'A site atoms (unit cell corners):'
       write(outunit,'(A6,A90,A9,A9)') '  name','charge tensor',
     &  'weight','dist.'
       do icorner=1,ncorner
       write(outunit,'(A6,9(F10.3),2(F9.3))') corners(icorner)%name, 
     &          corners(icorner)%charge,corners(icorner)%weight,
     &          corners(icorner)%distance
       end do
       write(outunit,*) ' ' 
       
       ! these atoms sit in the cell centers (B sites)
       call findstring(inunit,'center',6,found)
       if (.not.found) then
         print*,'keyword "center" not found. Quitting.'
         close(outunit)
         stop
       end if
       read(inunit,*) ncenter
       allocate(centername(ncenter))
       allocate(centers(ncenter))
       do icenter=1,ncenter
         read(inunit,*) centers(icenter)%name,
     &      centers(icenter)%charge(1,1:3),
     &      centers(icenter)%charge(2,1:3),
     &      centers(icenter)%charge(3,1:3) !,
!     &      centers(icenter)%weight,
!     &      centers(icenter)%distance
           centers(icenter)%weight=centerweight
           centers(icenter)%distance=centerdist
           centers(icenter)%geometry='center'
       centername(icenter)=centers(icenter)%name
       end do
       write(outunit,*) 'B site atoms (unit cell centers):'
       write(outunit,'(A6,A90,A9,A9)') '  name','charge tensor',
     &  'weight','dist.'
       do icenter=1,ncenter
          write(outunit,'(A6,9(F10.3),2(F9.3))')centers(icenter)%name, 
     &          centers(icenter)%charge,
     &          centers(icenter)%weight,centers(icenter)%distance
          !write(outunit,*) centername(icenter)
       end do
       write(outunit,*) ' ' 

	! these atoms sit on the cube faces (C or oxygen sites)
	call findstring(inunit,'face',4,found)
	if (.not.found) then
	  print*,'keyword "face" not found. Quitting.'
	  close(outunit)
	  stop
	end if
	read(inunit,*) nface
	nface=nface*3
	allocate(faces(nface))
	do iface=1,nface
	  read(inunit,*) faces(iface)%name,faces(iface)%direc,
     &         faces(iface)%charge(1,1:3),
     &         faces(iface)%charge(2,1:3),
     &         faces(iface)%charge(3,1:3) !,
!     &         faces(iface)%weight,faces(iface)%distance	
	       faces(iface)%weight=faceweight
	       faces(iface)%distance=facedist
	       faces(iface)%geometry='face  '
	end do
	write(outunit,*) 'Oxygen site atoms (unit cell faces):'
	write(outunit,'(A6,A3,A87,A9,A9)') '  name','#',
     &  'charge tensor', 'weight','dist.'
     	do iface=1,nface
	  write(outunit,'(A3,I3,9(F10.3),2(F9.3),A8)') faces(iface)%name, 
     &          faces(iface)%direc,faces(iface)%charge,
     &          faces(iface)%weight,faces(iface)%distance
	end do
	write(outunit,*) ' ' 
	
	! read if octahedral rotation should be performed
	call findstring(inunit,'octahedral_rotation',19,found)
	if (found) then
	  rotation=.True.
	  read(inunit,*) rotangle(1:3),alternate_rot(1:3)
	  write(outunit,*)'Octahedral rotation will be performed.'
	  write(outunit,'(A13,3(F6.3))')'angle (°): ', rotangle(1:3) 
	  write(outunit,*) 'Alternating angle: ',alternate_rot(1:3)
	  write(outunit,*)' '
	else
	  rotation=.False.
	end if

	! read if an AFE structure should be constructed
	call findstring(inunit,'AFE',3,found)
	if (found) then
	  AFE=.True.
	  ! read AFE mode, direction, amplitude for centers, ampl. for faces, alternate direction?
	  read(inunit,*) afemod,afevec,afecenteramp,afefaceamp,
     &         alternate_afe(1:3)
	  ! if there is no tetragonal symm, alternate displacement direction in neighboring cells
	  if (afemod.ne.'001') then
	    if (alternate_afe(1).ne.1.or.alternate_afe(2).ne.1.or.
     &        alternate_afe(3).ne.1) then
	      alternate_afe=1
	      print*, 'Warning (AFE): alternate set to 1 1 1.'
	      write(outunit,*) 'Warning: alternate set to 1 1 1.'
	    end if
	  end if
	  write(outunit,*)'AFE displacements will be printed.'
	  write(outunit,*)'AFE mode:', ' "',afemod,'"' 
	  write(outunit,'(A7,3(F6.3))')'axis: ', afevec 
	  write(outunit,*)
     &      'amplitude (in single cell lattice constants):' 
	  write(outunit,*) '  center atoms: ',afecenteramp 
	  write(outunit,*) '    face atoms: ',afefaceamp 
	  write(outunit,*) 
     &        'Alternating direction of disps: ',alternate_afe(1:3)
	  write(outunit,*)' '
	else
	  AFE=.False.
	end if

	! Read if polarization should be calculated
	call findstring(inunit,'polarization',12,found)
	if (found) then
	  write(outunit,*) 
     &     'The electric polarization will be calculated.'
	  write(outunit,*) ' '
	  polarization=.True.
	else 
	  polarization=.False.
	end if

	! Read if reference coordinates should be read and used
	call findstring(inunit,'reference',9,found)
	if (found) then
	  write(outunit,*) 'Reference coordinates will be read 
     & from file "COORAT.REF".'
	  write(outunit,*) ' '
	  reference=.True.
	else 
	  reference=.False.
	end if

	! Read if displacements should be printed
	call findstring(inunit,'print_displacements',19,found)
	if (found) then
	  write(outunit,*) 'The atomic displacements will be printed to'
	  if (reference) write(outunit,*) 
     &      '"DISP.vscub.DAT" and "DISP.vsref.DAT".'
	  if (.not.reference) write(outunit,*) 
     &      '"DISP.vscub.DAT".'
	  write(outunit,*) ' '
	  print_disps=.True.
	else 
	  print_disps=.False.
	end if

	! Read if xsf file should be printed
	call findstring(inunit,'print_xsf',9,found)
	if (found) then
	  write(outunit,*) '(An) xsf file(s) will be printed to'
          write(outunit,*) '"GEOMETRY.IN.xsf".'
	  if (reference) write(outunit,*) 
     &      '"GEOMETRY.REF.xsf"'
	  if (rotation) write(outunit,*) 
     &      '"GEOMETRY.ROT.xsf"'
	  if (AFE) write(outunit,*) 
     &      '"GEOMETRY.AFE.xsf"'
	  write(outunit,*) ' '
	  print_xsf=.True.
	else 
	  print_xsf=.False.
	end if

	! Determine what would be the ideal cubic (paraelectric) positions
	ncubcorner=nuc(1)*nuc(2)*nuc(3)+(nuc(1)+1)*(nuc(2)+1)
     &             +(nuc(2)+1)*nuc(3)+(nuc(1))*nuc(3)
	ncubcenter=nuc(1)*nuc(2)*nuc(3)
	ncubface=3*nuc(1)*nuc(2)*nuc(3)+nuc(1)*nuc(3)+nuc(2)*nuc(3)
     &           +nuc(1)*nuc(2)  
	allocate(cubcorners(1:3,1:ncubcorner),
     &      cubcenters(1:3,1:ncubcenter),
     &      cubfaces(1:3,1:ncubface))
	! corners
	irun=1
	do n1=0,nuc(1)
	  do n2=0,nuc(2)
	    do n3=0,nuc(3)
	      cubcorners(1,irun)=dble(n1)/dble(nuc(1))
	      cubcorners(2,irun)=dble(n2)/dble(nuc(2))
	      cubcorners(3,irun)=dble(n3)/dble(nuc(3))
	      irun=irun+1
	    end do
	  end do
	end do
	! centers
	irun=1
	do n1=0,nuc(1)-1
	  do n2=0,nuc(2)-1
	    do n3=0,nuc(3)-1
	      cubcenters(1,irun)=(0.5D0+dble(n1))/dble(nuc(1))
	      cubcenters(2,irun)=(0.5D0+dble(n2))/dble(nuc(2))
	      cubcenters(3,irun)=(0.5D0+dble(n3))/dble(nuc(3))
	      irun=irun+1
	    end do
	  end do
	end do
	! faces
	! O1
	irun=1
	do n1=0,nuc(1)
	  do n2=0,nuc(2)-1
	    do n3=0,nuc(3)-1
	      cubfaces(1,irun)=dble(n1)/dble(nuc(1))
	      cubfaces(2,irun)=(0.5D0+dble(n2))/dble(nuc(2))
	      cubfaces(3,irun)=(0.5D0+dble(n3))/dble(nuc(3))
	      irun=irun+1
	    end do
	  end do
	end do
	! O2
	do n1=0,nuc(1)-1
	  do n2=0,nuc(2)
	    do n3=0,nuc(3)-1
	      cubfaces(1,irun)=(0.5D0+dble(n1))/dble(nuc(1))
	      cubfaces(2,irun)=dble(n2)/dble(nuc(2))
	      cubfaces(3,irun)=(0.5D0+dble(n3))/dble(nuc(3))
	      irun=irun+1
	    end do
	  end do
	end do
	! O3
	do n1=0,nuc(1)-1
	  do n2=0,nuc(2)-1
	    do n3=0,nuc(3)
	      cubfaces(1,irun)=(0.5D0+dble(n1))/dble(nuc(1))
	      cubfaces(2,irun)=(0.5D0+dble(n2))/dble(nuc(2))
	      cubfaces(3,irun)=dble(n3)/dble(nuc(3))
	      irun=irun+1
	    end do
	  end do
	end do
	write(outunit,*) 'For comparison: The corresponding 
     & cubic positions are' 
	write(outunit,*) 'corners:'
	do i=1,ncubcorner
	  write(outunit,'(3(F10.6))') cubcorners(1:3,i)
	end do 
	write(outunit,*) 'centers:'
	do i=1,ncubcenter
	  write(outunit,'(3(F10.6))') cubcenters(1:3,i)
	end do 
	write(outunit,*) 'faces:'
	do i=1,ncubface
	  write(outunit,'(3(F10.6))') cubfaces(1:3,i)
	end do 
	write(outunit,*) ' '
	  
	! determine number of atoms and elements in file COORAT.IN
	nlines=0
	natoms=0
	nelements=0
	!open(10,file='COORAT.IN',status='old')
	open(30,file='COORAT.IN',status='old')
!10	read(10,'(A80)',end=20) line
10	read(30,'(A80)',end=20) line
	nlines=nlines+1
	if (index(line,'natom').gt.0) nelements=nelements+1	
	if (index(line,'.').gt.0) natoms=natoms+1	
	goto 10

20      continue

	allocate(pelements(nelements))
	    
	allocate(atoms(1:natoms))
!	rewind(10)
	rewind(30)

	! get coordinates of atoms in the periodic box	
!	call getatoms(10,pelements,nelements,atoms,natoms,
	call getatoms(30,pelements,nelements,atoms,natoms,
     &                      natomsXL,tol,
     &                      corners,centers,faces,
     &                     ncorner,ncenter,nface,
     &                    cubcorners,cubcenters,cubfaces,
     &                    ncubcorner,ncubcenter,ncubface  )
	if (print_xsf) call write_xsf('GEOMETRY.IN.xsf',15,atoms,
     &                                natoms,avecs)
	do iatom=1,natoms
	  thisatom=atoms(iatom)
	  ! corners
	  do icorner=1,ncorner
	    if (thisatom%name==corners(icorner)%name) then  ! atoms is corner-atom?
              atoms(iatom)%direc=1 
              atoms(iatom)%charge=corners(icorner)%charge
              atoms(iatom)%weight=corners(icorner)%weight
              atoms(iatom)%distance=corners(icorner)%distance
              atoms(iatom)%geometry=corners(icorner)%geometry
	    end if
	  end do
	  ! centers
	  do icenter=1,ncenter
	    if (thisatom%name==centers(icenter)%name) then  ! atom is a center atom?
              atoms(iatom)%direc=1 
              atoms(iatom)%charge=centers(icenter)%charge
              atoms(iatom)%weight=centers(icenter)%weight
              atoms(iatom)%distance=centers(icenter)%distance
              atoms(iatom)%geometry=centers(icenter)%geometry
	    end if
	  end do
	  ! faces
	  do iface=1,nface
	    if (thisatom%name==faces(iface)%name) then   ! atom is a face atom?
	      thisatom%distance=faces(iface)%distance
	      do iatom2=1,natoms                         ! find central atom of cell it belongs to 
	        do icenter=1,ncenter
	          if (atoms(iatom2)%name.eq.
     &              centers(icenter)%name) then
		          centeratom=atoms(iatom2)
		    cubatom=thisatom
	 	    cubatom%where=nearestcubsite(thisatom%where,
     &                thisatom%geometry,
     &                cubcorners,cubcenters,cubfaces,ncubcorner,
     &                ncubcenter,ncubface)
	            if (isincell(cubatom,centeratom,nuc)) then
	              thisatom%direc=                    ! O1, O2 or O3?      
     &                  direction(thisatom,centeratom)
	            end if
	  	end if
	        end do
	      end do
	      if (thisatom%direc==faces(iface)%direc) then
                atoms(iatom)%direc=faces(iface)%direc
                atoms(iatom)%charge=faces(iface)%charge
                atoms(iatom)%weight=faces(iface)%weight
                atoms(iatom)%distance=faces(iface)%distance
                atoms(iatom)%geometry=faces(iface)%geometry
          end if
         end if
      end do
      end do

40      continue
    
      write(outunit,*) 'Information from COORAT.IN:'
      write(outunit,*) nlines, ' lines, ',natoms,' atoms,',
     &  nelements,' elements:'
      write(outunit,*) 'name     number of atoms'

      do iele=1,nelements
      write(outunit,5555) pelements(iele)%name,
     &          pelements(iele)%howmany
      end do
      write(outunit,*) ' '

      write(outunit,*) natomsXL,'atoms in extended UC:'
      write(outunit,'("(number, name, cart. direction in which the next   &
     &central atom sits (O only, else 0), 3 cart. coords, charge tensor,  &
     & weight, max. distance to central atom of the same cell)")')
!     close(10)
      close(30)

      ! get atoms of the extended unit cell
      allocate(atomsXL(natomsXL))
      call getatomsXL(atoms,natoms,atomsXL,natomsXL,pelements,
     &                        nelements,tol,
     &                      corners,centers,faces,
     &                     ncorner,ncenter,nface,
     &                    cubcorners,cubcenters,cubfaces,
     &                    ncubcorner,ncubcenter,ncubface  )
      do iatom=1,natomsXL
       thisatom=atomsXL(iatom)
       ! corners
       do icorner=1,ncorner
         if (thisatom%name==corners(icorner)%name) then  ! atoms is corner-atom?
              atomsXL(iatom)%direc=1 
              atomsXL(iatom)%charge=corners(icorner)%charge
              atomsXL(iatom)%weight=corners(icorner)%weight
              atomsXL(iatom)%distance=corners(icorner)%distance
              atomsXL(iatom)%geometry=corners(icorner)%geometry
         end if
       end do
       ! centers
       do icenter=1,ncenter
         if (thisatom%name==centers(icenter)%name) then  ! atom is a center atom?
              atomsXL(iatom)%direc=1 
              atomsXL(iatom)%charge=centers(icenter)%charge
              atomsXL(iatom)%weight=centers(icenter)%weight
              atomsXL(iatom)%distance=centers(icenter)%distance
              atomsXL(iatom)%geometry=centers(icenter)%geometry
         end if
       end do
       ! faces
       do iface=1,nface
         if (thisatom%name==faces(iface)%name) then   ! atom is a face atom?
           thisatom%distance=faces(iface)%distance
           do iatom2=1,natoms                         ! find central atom of cell it belongs to 
             do icenter=1,ncenter
               if (atoms(iatom2)%name.eq.
     &              centers(icenter)%name) then
              centeratom=atoms(iatom2)
        cubatom=thisatom
        cubatom%where=nearestcubsite(thisatom%where,
     &                thisatom%geometry,
     &                cubcorners,cubcenters,cubfaces,ncubcorner,
     &                ncubcenter,ncubface)
                if (isincell(cubatom,centeratom,nuc)) then
                  thisatom%direc=                    ! O1, O2 or O3?      
     &                  direction(thisatom,centeratom)
                end if
        end if
            end do
          end do
          if (thisatom%direc==faces(iface)%direc) then
                atomsXL(iatom)%direc=faces(iface)%direc
                atomsXL(iatom)%charge=faces(iface)%charge
                atomsXL(iatom)%weight=faces(iface)%weight
                atomsXL(iatom)%distance=faces(iface)%distance
                atomsXL(iatom)%geometry=faces(iface)%geometry
          end if
        end if
       end do
      end do
      ! write atom props to output file
      do iatom=1,natomsXL
            write(outunit,4444) atomsXL(iatom)  
      end do

      ! get reference coordinates
      call findstring(inunit,'reference',9,found)
      if (found) then
       allocate(atomsref(1:natoms))
!	  open(16,file='COORAT.REF',status='old')
      open(36,file='COORAT.REF',status='old')
!	  call getatoms(16,pelements,nelements,atomsref,natoms,
      call getatoms(36,pelements,nelements,atomsref,natoms,
     &                      natomsrefXL,tol,
     &                      corners,centers,faces,
     &                     ncorner,ncenter,nface,
     &                    cubcorners,cubcenters,cubfaces,
     &                    ncubcorner,ncubcenter,ncubface  )
      if (print_xsf) call write_xsf('GEOMETRY.REF.xsf',16,
     &               atomsref,natoms,avecs)
      do iatom=1,natoms
        thisatom=atomsref(iatom)
        ! corners
        do icorner=1,ncorner
          if (thisatom%name==corners(icorner)%name) then  ! atoms is corner-atom?
                atomsref(iatom)%direc=1 
                atomsref(iatom)%charge=corners(icorner)%charge
                atomsref(iatom)%weight=corners(icorner)%weight
                atomsref(iatom)%distance=corners(icorner)%distance
                atomsref(iatom)%geometry=corners(icorner)%geometry
          end if
        end do
        ! centers
        do icenter=1,ncenter
          if (thisatom%name==centers(icenter)%name) then  ! atom is a center atom?
                atomsref(iatom)%direc=1 
                atomsref(iatom)%charge=centers(icenter)%charge
                atomsref(iatom)%weight=centers(icenter)%weight
                atomsref(iatom)%distance=centers(icenter)%distance
                atomsref(iatom)%geometry=centers(icenter)%geometry
          end if
        end do
        ! faces
        do iface=1,nface
          if (thisatom%name==faces(iface)%name) then   ! atom is a face atom?
            thisatom%distance=faces(iface)%distance
            do iatom2=1,natoms                         ! find central atom of cell it belongs to 
              do icenter=1,ncenter
                if (atomsref(iatom2)%name.eq.
     &                centers(icenter)%name) then
                  centeratom=atomsref(iatom2)
                  cubatom=thisatom
                  cubatom%where=nearestcubsite(thisatom%where,
     &                  thisatom%geometry,
     &                  cubcorners,cubcenters,cubfaces,ncubcorner,
     &                  ncubcenter,ncubface)
                    if (isincell(cubatom,centeratom,nuc)) then
                      thisatom%direc=                    ! O1, O2 or O3?      
     &                    direction(thisatom,centeratom)
                  end if
                end if
              end do
            end do
            if (thisatom%direc==faces(iface)%direc) then
                  atomsref(iatom)%direc=faces(iface)%direc
                  atomsref(iatom)%charge=faces(iface)%charge
                  atomsref(iatom)%weight=faces(iface)%weight
                  atomsref(iatom)%distance=faces(iface)%distance
                  atomsref(iatom)%geometry=faces(iface)%geometry
            end if
          end if
        end do
      end do
      allocate(atomsrefXL(natomsrefXL))
      write(outunit,*)'There are',natomsrefXL,'atoms in the ext. 
     & ref. cell:'
      call getatomsXL(atomsref,natoms,atomsrefXL,natomsrefXL,
     &                    pelements,nelements,tol,
     &                      corners,centers,faces,
     &                     ncorner,ncenter,nface,
     &                    cubcorners,cubcenters,cubfaces,
     &                    ncubcorner,ncubcenter,ncubface  )
      do iatom=1,natomsXL
        thisatom=atomsrefXL(iatom)
        ! corners
        do icorner=1,ncorner
          if (thisatom%name==corners(icorner)%name) then  ! atom is corner-atom?
                atomsrefXL(iatom)%direc=1  
                atomsrefXL(iatom)%charge=corners(icorner)%charge
                atomsrefXL(iatom)%weight=corners(icorner)%weight
                atomsrefXL(iatom)%distance=corners(icorner)%distance
                atomsrefXL(iatom)%geometry=corners(icorner)%geometry
          end if
        end do
        ! centers
        do icenter=1,ncenter
          if (thisatom%name==centers(icenter)%name) then  ! atom is a center atom?
                atomsrefXL(iatom)%direc=1 
                atomsrefXL(iatom)%charge=centers(icenter)%charge
                atomsrefXL(iatom)%weight=centers(icenter)%weight
                atomsrefXL(iatom)%distance=centers(icenter)%distance
                atomsrefXL(iatom)%geometry=centers(icenter)%geometry
          end if
        end do
        ! faces
        do iface=1,nface
          if (thisatom%name==faces(iface)%name) then   ! atom is a face atom?
            thisatom%distance=faces(iface)%distance
            do iatom2=1,natoms                         ! find central atom of cell it belongs to 
              do icenter=1,ncenter
                if (atoms(iatom2)%name.eq.
     &                centers(icenter)%name) then
                 centeratom=atoms(iatom2)
                 cubatom=thisatom
                 cubatom%where=nearestcubsite(thisatom%where,
     &                  thisatom%geometry,
     &                  cubcorners,cubcenters,cubfaces,ncubcorner,
     &                  ncubcenter,ncubface)
                    if (isincell(cubatom,centeratom,nuc)) then
                    thisatom%direc=                    ! O1, O2 or O3?      
     &                    direction(thisatom,centeratom)
                  end if
               end if
              end do
            end do
            if (thisatom%direc==faces(iface)%direc) then
                  atomsrefXL(iatom)%direc=faces(iface)%direc
                  atomsrefXL(iatom)%charge=faces(iface)%charge
                  atomsrefXL(iatom)%weight=faces(iface)%weight
                  atomsrefXL(iatom)%distance=faces(iface)%distance
                  atomsrefXL(iatom)%geometry=faces(iface)%geometry
            end if
          end if
        end do
      end do
      do iatom=1,natomsrefXL
              write(outunit,4444) atomsrefXL(iatom)  
      end do
!	  open(18,file='COORDS.REF',status='replace')
      open(38,file='COORDS.REF',status='replace')
      do iatomXL=1,natomsrefXL
!	    write(18,'(A4,3(F12.6))') atomsrefXL(iatomXL)%name,
        write(38,'(A4,3(F12.6))') atomsrefXL(iatomXL)%name,
     &               atomsrefXL(iatomXL)%where(1:3)
      end do
!	  close(18)
      close(38)
!	  close(16)
      close(36)
      end if
    
      ! get atoms that belong to cell of a center atom 
      ncell=0
      do iele=1,nelements
       do icenter=1,ncenter
        if (pelements(iele)%name.eq.centername(icenter)) 
     &        ncell=ncell+pelements(iele)%howmany
       end do
      end do
        write(outunit,*) ncell,'cell(s) in the box.'
c	totalvol=alat(1)*alat(2)*alat(3) ! volume of the supercell in Bohr^3
      totalvol=avecs(1,1)*(avecs(2,2)*avecs(3,3)
     &          -avecs(2,3)*avecs(3,2))
     &          +avecs(1,2)*(avecs(2,3)*avecs(3,1)
     &          -avecs(2,1)*avecs(3,3))
     &          +avecs(1,3)*(avecs(2,1)*avecs(3,2)
     &          -avecs(2,2)*avecs(3,1))
      totalvol=abs(totalvol)
      avvol=totalvol/dble(ncell) ! averaged volume per cell 
      write(outunit,*)'The volume of the whole box is',totalvol,
     &             'Bohr^3.'      
      write(outunit,*)'The averaged volume per cell is ',avvol,
     &   'Bohr^3.'

      ! calc. local dipole moment
      if (polarization) then
          allocate(dipol(1:ncell,1:6))
       ! dipole moments from displacements w.r.t. central atom
       call dipole_center(natoms,natomsXL,atoms,atomsXL,avecs,         
     &   centername,ncenter,ncell,dipol,nuc,cubcorners,cubcenters,      
     &   cubfaces,ncubcorner,ncubcenter,ncubface)
 
      ! print dipole moment of each cell to file
!         open(12,file='DIPOL.CENTER.DAT',status='replace')
      open(32,file='DIPOL.CENTER.DAT',status='replace')
      do icell=1,ncell
       !write(12,1000) dipol(icell,1:6)
       write(32,1000) dipol(icell,1:6)
      end do
      !close(12)
      close(32)
      ! calculate total dipole moment and polarization
      dipoltot=0.0D0
      do icell=1,ncell
        dipoltot(1:3)=dipoltot+dipol(icell,4:6)
      end do
      write(outunit,*) 'using displ. from central atom:'
      write(outunit,6666) 'total dipole moment (e*Angs):',
     &                        dipoltot(1:3)
      polfac=ec*1.0E22/(totalvol*bohr**3)
      write(outunit,6666) 'total polarization (uC/cm^2):',
     &                        dipoltot(1:3)*polfac
      ! dipole moments from displacements w.r.t. cubic positions
      call dipole_cube(natoms,natomsXL,atoms,atomsXL,avecs,
     &                    centername,ncenter,ncell,dipol,nuc,
     &                    cubcorners,cubcenters,cubfaces,
     &                    ncubcorner,ncubcenter,ncubface  )

      ! print dipole moment of each cell to file
      !open(12,file='DIPOL.CUBE.DAT',status='replace')
      open(32,file='DIPOL.CUBE.DAT',status='replace')
      do icell=1,ncell
        !write(12,1000) dipol(icell,1:6)
        write(32,1000) dipol(icell,1:6)
      end do
      !close(12)
      close(32)
      ! calculate total dipole moment and polarization
      dipoltot=0.0D0
      do icell=1,ncell
        dipoltot(1:3)=dipoltot+dipol(icell,4:6)
      end do
      write(outunit,*) 'using displ. from cubic positions:'
      write(outunit,6666) 'total dipole moment (e*Angs):',
     &                        dipoltot(1:3)
      polfac=ec*1.0E22/(totalvol*bohr**3)
      write(outunit,6666) 'total polarization (uC/cm^2):',
     &                        dipoltot(1:3)*polfac
      end if  

    
      ! print displacements of all atoms with respect to cubic phase to file
      if (print_disps) then
      !open(15,file='DISP.vscub.DAT',status='replace')
      open(35,file='DISP.vscub.DAT',status='replace')
      do iatomXL=1,natomsXL
        thisatom=atomsXL(iatomXL)
        ! disp w.r.t. cubic phase
        disp=2.0D0
        x(1:3)=thisatom%where(1:3)
        ! if atom is corneratom, look for nearest cubic corner position
        if (thisatom%geometry.eq.'corner') then
          do i=1,ncubcorner
           if ((x(1)-cubcorners(1,i))**2+(x(2)
     &             -cubcorners(2,i))**2+(x(3)-cubcorners(3,i))**2
     &             .lt.disp(1)**2+disp(2)**2+disp(3)**2)
     &             disp(1:3)=x(1:3)-cubcorners(1:3,i)
          end do
        end if 
        ! if atom is centeratom, look for nearest cubic center position
        if (thisatom%geometry.eq.'center') then
          do i=1,ncubcenter
            if ((x(1)-cubcenters(1,i))**2+(x(2)
     &             -cubcenters(2,i))**2+(x(3)-cubcenters(3,i))**2
     &             .lt.disp(1)**2+disp(2)**2+disp(3)**2)
     &            disp(1:3)=x(1:3)-cubcenters(1:3,i)
          end do
        end if 
        ! if atom is faceatom, look for nearest cubic face position
        if (thisatom%geometry.eq.'face  ') then
          do i=1,ncubface
            if ((x(1)-cubfaces(1,i))**2+(x(2)
     &             -cubfaces(2,i))**2+(x(3)-cubfaces(3,i))**2
     &             .lt.disp(1)**2+disp(2)**2+disp(3)**2)
     &            disp(1:3)=x(1:3)-cubfaces(1:3,i)
          end do
        end if 
        !do i=1,ncubz
        !  if (abs(x(3)-cubicz(i)).lt.abs(disp(3))) 
     &      !      disp(3)=x(3)-cubicz(i)
        !end do 
        !write(15,2000) thisatom%name,x(1:3),disp(1:3)
        write(35,2000) thisatom%name,x(1:3),disp(1:3)
       end do
       !close(15)
       close(35)
      end if

      ! print displacements of all atoms with respect to reference coordinates to file
      if (reference.and.print_disps) then
       !open(17,file='DISP.vsref.DAT',status='replace')
       open(37,file='DISP.vsref.DAT',status='replace')
       do iatomXL=1,natomsXL
        thisatom=atomsXL(iatomXL)
        ! disp w.r.t. reference structure
        disp=1.0D0
        x(1:3)=thisatom%where(1:3)
        do iatomXL2=1,natomsrefXL
          thisrefatom=atomsrefXL(iatomXL2)
          xref(1:3)=thisrefatom%where(1:3)
          if (thisatom%name.eq.thisrefatom%name) then
            if ((x(1)-xref(1))**2+(x(2)-xref(2))**2
     &              +(x(3)-xref(3))**2.lt.disp(1)**2+disp(2)**2
     &              +disp(3)**2) disp=x-xref
          end if
        end do
        !write(17,2000) thisatom%name,x(1:3),disp(1:3)
        write(37,2000) thisatom%name,x(1:3),disp(1:3)
       end do
       !close(17)
       close(37)
      end if

      ! calculate positions of atoms in rotated octahedra
      if (rotation) then
      call rotate(atoms,natoms,pelements,nelements,
     &     rotangle,alternate_rot,faces,nface,centers,
     &     ncenter,nuc,outunit,cubcorners,cubcenters,cubfaces,
     &          ncubcorner,ncubcenter,ncubface,avecs,print_xsf)
      end if

      ! calculate positions of atoms for AFE structure
      if (AFE) then
       call AFE_DISP(atoms,natoms,pelements,nelements,afemod,
     &     afevec,afecenteramp,afefaceamp,alternate_afe,faces,
     &     nface,centers,ncenter,nuc,outunit,
     &                    cubcorners,cubcenters,cubfaces,
     &                    ncubcorner,ncubcenter,ncubface,
     &                    print_xsf,avecs  )
      end if

 1000 format(3(F12.6),3(F12.3))
 2000 format(A3,3(F12.3),3(F12.6))
 4444 format(I3, A3, I3, 14(F10.3),A8) 
 5555 format(A2,I22)
 6666 format(A30,3(F12.6))
      
      call cpu_time(t_end)
      write(outunit,*) ' '
      write(outunit,*) 'Computation time in s: ',t_end-t_start

      close(inunit)
      close(outunit)

      write(12,fsubend)
      return
      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine findstring(fileunit,string,stringlen,stringfound)
	implicit none
	integer, intent(in) :: fileunit,stringlen
	character (len=*), intent(in) :: string
	logical, intent(out) :: stringfound
	integer istring,icomment
	character line*120,comment
	stringfound=.False.
	comment='#'
	rewind(fileunit)
9000    read(fileunit,'(A120)',end=9001)line
	if (index(line,string).ge.1) then
	  stringfound=.True.  
	  istring=1
	  do while(line(istring:istring+stringlen-1).ne.string)
	    istring=istring+1
	  end do
	  if (index(line,comment).ge.1) then
	    icomment=1
	    do while(line(icomment:icomment).ne.comment)
	      icomment=icomment+1
	    end do
	    if (istring.gt.icomment) stringfound=.False.  
	  end if
  	  if (stringfound) goto 9001 ! formerly this line read: goto 9001
	end if
	goto 9000 	
9001 	continue
	end subroutine

      subroutine getatoms(fileunit,pelements,nelements,atoms,natoms,
     &                      natomsXL,tol,
     &                      corners,centers,faces,
     &                     ncorner,ncenter,nface,
     &                    cubcorners,cubcenters,cubfaces,
     &                    ncubcorner,ncubcenter,ncubface )
c	use functions
      use defs
      implicit none
      integer, intent(in) :: fileunit,ncorner,ncenter,nface
      integer, intent(in) :: ncubcorner,ncubcenter,ncubface
c      type patom
c        sequence
c        integer:: number
c        character:: name*2 
c        integer:: direc 
c        double precision:: where(3)
c        double precision:: charge(1:3,1:3)
c        double precision:: weight
c        double precision:: distance
c        character:: geometry*6
c      end type patom
c      type pelement
c        sequence
c        integer number
c        integer:: howmany
c        character:: name*2 
c        double precision:: charge(1:3,1:3)
c        double precision:: weight
c        double precision:: distance
c      end type pelement
      type(patom) :: atoms(1:natoms)
      type(patom),intent(in) :: corners(1:ncorner)
      type(patom),intent(in) :: centers(1:ncenter)
      type(patom),intent(in) :: faces(1:nface)
      type(pelement) :: pelements(1:nelements)
      double precision, intent(in) :: tol(1:3)
      double precision, intent(in) :: cubcorners(1:3,1:ncubcorner)
      double precision, intent(in) :: cubcenters(1:3,1:ncubcenter)
      double precision, intent(in) :: cubfaces(1:3,1:ncubface)
      integer,intent(in) :: natoms,nelements
      integer,intent(inout) :: natomsXL
      integer iatom,iele,iread,ireadend,natom,n1,n2,n3,i,n
      integer icorner,icenter,iface
      double precision x(1:3),dist(1:3),cubx(1:3)
      character line*80,nameat*2
      iatom=1
      iele=1
      natomsXL=0
 9010 read(fileunit,'(A80)',end=9011) line
        if (index(line,'.').le.0.and.index(line,'natom').le.0)
     &     goto 9010  
      ! find number of atoms of a species
      iread=1
      do while (line(iread:iread+4).ne.'natom')
        iread=iread+1
      end do
      do while (line(iread:iread).ne.'=')
        iread=iread+1
      end do
      iread=iread+1
      do while (line(iread:iread).eq.' ')
        iread=iread+1
      end do
      ireadend=iread+1
      do while (line(ireadend:ireadend).ne.' ')
        ireadend=ireadend+1
      end do
!        read(line(iread:ireadend),'(I)') natom
        read(line(iread:ireadend),*) natom
	! find name of a species
	iread=1
	do while (line(iread:iread+3).ne.'name')
	  iread=iread+1
	end do
	do while (line(iread:iread).ne.'=')
	  iread=iread+1
	end do
	iread=iread+1
	do while (line(iread:iread).eq.' ')
	  iread=iread+1
	end do
	ireadend=iread+1
	nameat=line(iread:iread+1)
	pelements(iele)%number=iele
	pelements(iele)%howmany=natom
	pelements(iele)%name=nameat
	iele=iele+1
	! read coordinates of this species
	do n=1,natom
	  read(fileunit,*) x(1:3)
          ! shift atoms to position between about 0 and 1
          do i=1,3
            do while (x(i).lt.-tol(i))
              x(i)=x(i)+1.0D0
            end do  
            do while (x(i).ge.1.0D0-tol(i))
              x(i)=x(i)-1.0D0
            end do  
          end do
	  atoms(iatom)%number=iatom
	  atoms(iatom)%name=nameat
	  atoms(iatom)%direc=0
	  atoms(iatom)%where(1:3)=x(1:3)
	  ! determine if atom is corner, center or face type
	  do icorner=1,ncorner
	    if (nameat.eq.corners(icorner)%name)
     &          atoms(iatom)%geometry='corner'
	  end do
	  do icenter=1,ncenter
	    if (nameat.eq.centers(icenter)%name)
     &          atoms(iatom)%geometry='center'
	  end do
	  do iface=1,nface
	    if (nameat.eq.faces(iface)%name)
     &        atoms(iatom)%geometry='face  '
	  end do
	  cubx=nearestcubsite(x,
     &      atoms(iatom)%geometry,
     &      cubcorners,cubcenters,cubfaces,ncubcorner,
     &      ncubcenter,ncubface)
	  ! are the atoms at the unit cell border? then add their images on the other side of the unit cell 
!	  do n1=0,1
!	    do n2=0,1
!	      do n3=0,1
!	        if (x(1)+dble(n1)<1.0D0+tol(1).and.x(2)+dble(n2)
!     &              <1.0D0+tol(2).and.x(3)+dble(n3)<1.0D0+tol(3)) 
!     &              natomsXL=natomsXL+1
!	      end do
!	    end do
!	  end do  	
	  do n1=-1,1
	    do n2=-1,1
	      do n3=-1,1
	        if (cubx(1)+dble(n1)>-tol(1).and.
     &              cubx(1)+dble(n1)<1.0D0+tol(1).and.
     &              cubx(2)+dble(n2)>-tol(2).and.
     &              cubx(2)+dble(n2)<1.0D0+tol(2).and.
     &              cubx(3)+dble(n3)>-tol(3).and.
     &              cubx(3)+dble(n3)<1.0D0+tol(3)) 
     &              natomsXL=natomsXL+1
	      end do
	    end do
	  end do  	
	  iatom=iatom+1
	end do
	goto 9010
9011	continue
	end subroutine

	subroutine getatomsXL(atoms,natoms,atomsXL,natomsXL,pelements,
     &                        nelements,tol,
     &                      corners,centers,faces,
     &                     ncorner,ncenter,nface,
     &                    cubcorners,cubcenters,cubfaces,
     &                    ncubcorner,ncubcenter,ncubface  )
c      use functions
      use defs
      implicit none
      integer, intent(in) :: natoms,natomsXL,nelements
      integer, intent(in) :: ncorner,ncenter,nface
      integer, intent(in) :: ncubcorner,ncubcenter,ncubface
      double precision, intent(in) :: tol(1:3)
      double precision, intent(in) :: cubcorners(1:3,1:ncubcorner)
      double precision, intent(in) :: cubcenters(1:3,1:ncubcenter)
      double precision, intent(in) :: cubfaces(1:3,1:ncubface)
c      type patom
c        sequence
c        integer:: number
c        character:: name*2 
c        integer:: direc 
c        double precision:: where(3)
c        double precision:: charge(1:3,1:3)
c        double precision:: weight
c        double precision:: distance
c        character:: geometry*6
c      end type patom
c	type pelement
c   	  sequence
c	  integer number
c	  integer:: howmany
c	  character:: name*2 
c   	  double precision:: charge(1:3,1:3)
c   	  double precision:: weight
c   	  double precision:: distance
c	end type pelement
      type(patom) atoms(1:natoms),atomsXL(1:natomsXL)
      type(patom),intent(in) :: corners(1:ncorner)
      type(patom),intent(in) :: centers(1:ncenter)
      type(patom),intent(in) :: faces(1:nface)
      type(pelement) pelements(1:nelements)
      integer iatom,iatomXL,n1,n2,n3,iele
      double precision x(1:3),cubx(1:3)
      iatomXL=0
      do iatom=1,natoms
      x=atoms(iatom)%where
      cubx=nearestcubsite(x,
     &      atoms(iatom)%geometry,
     &      cubcorners,cubcenters,cubfaces,ncubcorner,
     &      ncubcenter,ncubface)
!      do n1=0,1
!       do n2=0,1
!        do n3=0,1
!          if (x(1)+dble(n1)<1.0D0+tol(1).and.x(2)+dble(n2)
!      &              <1.0D0+tol(2).and.x(3)+dble(n3)<1.0D0+tol(3))then
!                     iatomXL=iatomXL+1
! 		    atomsXL(iatomXL)=atoms(iatom)
! 		    atomsXL(iatomXL)%number=iatomXL
! 		    atomsXL(iatomXL)%where(1)
!      &                 =x(1)+dble(n1)
! 		    atomsXL(iatomXL)%where(2)
!      &                 =x(2)+dble(n2)
! 		    atomsXL(iatomXL)%where(3)
!      &                 =x(3)+dble(n3)
! 		    do iele=1,nelements
! 		      if (atomsXL(iatomXL)%name(1:2).eq.
!      &                    pelements(iele)%name(1:2)) then
! 			atomsXL(iatomXL)%charge=pelements(iele)%charge
! 			atomsXL(iatomXL)%weight=pelements(iele)%weight
! 			atomsXL(iatomXL)%distance=pelements(iele)%distance
! 		      end if		
! 		    end do
! 		end if
! 	      end do
! 	    end do
! 	  end do  	
       do n1=-1,1
         do n2=-1,1
           do n3=-1,1
            if (cubx(1)+dble(n1)>-tol(1).and.
     &              cubx(1)+dble(n1)<1.0D0+tol(1).and.
     &              cubx(2)+dble(n2)>-tol(2).and.
     &              cubx(2)+dble(n2)<1.0D0+tol(2).and.
     &              cubx(3)+dble(n3)>-tol(3).and.
     &              cubx(3)+dble(n3)<1.0D0+tol(3)) then 
                    iatomXL=iatomXL+1
              atomsXL(iatomXL)=atoms(iatom)
              atomsXL(iatomXL)%number=iatomXL
              atomsXL(iatomXL)%where(1)
     &                 =x(1)+dble(n1)
              atomsXL(iatomXL)%where(2)
     &                 =x(2)+dble(n2)
              atomsXL(iatomXL)%where(3)
     &                 =x(3)+dble(n3)
              do iele=1,nelements
                if (atomsXL(iatomXL)%name(1:2).eq.
     &                    pelements(iele)%name(1:2)) then
                  atomsXL(iatomXL)%charge=pelements(iele)%charge
                  atomsXL(iatomXL)%weight=pelements(iele)%weight
                  atomsXL(iatomXL)%distance=pelements(iele)%distance
                  atomsXL(iatomXL)%geometry=atoms(iatom)%geometry
                end if
              end do
            end if
          end do
        end do
       end do  
      end do
      end subroutine getatomsXL
 
      subroutine dipole_center(natoms,natomsXL,atoms,atomsXL,avecs,       &
     &                    centername,ncenter,ncell,dipol,nuc,             &
     &                    cubcorners,cubcenters,cubfaces,                 &
     &                    ncubcorner,ncubcenter,ncubface  )               &
      use defs
c     use functions
      implicit none
      integer, intent(in) :: natoms,natomsXL,ncenter,ncell,nuc(1:3)
      integer, intent(in) :: ncubcorner,ncubcenter,ncubface
      double precision, intent(in) :: avecs(1:3,1:3)
      double precision, intent(in) :: cubcorners(1:3,1:ncubcorner)
      double precision, intent(in) :: cubcenters(1:3,1:ncubcenter)
      double precision, intent(in) :: cubfaces(1:3,1:ncubface)
      double precision,intent(inout) :: dipol(1:ncell,1:6)
      character (len=2), intent(in) :: centername(1:ncenter)
c	type patom
c   	  sequence
c   	  integer:: number
c	  character:: name*2 
c	  integer:: direc 
c   	  double precision:: where(3)
c   	  double precision:: charge(1:3,1:3)
c   	  double precision:: weight
c   	  double precision:: distance
c	  character:: geometry*6
c	end type patom
      type(patom),intent(in) :: atoms(1:natoms),atomsXL(1:natomsXL)
      type(patom) centeratom,thisatom,cubatom
      double precision dipolfac,distvec(1:3),locdipol(1:3)
      integer icell,iatom,iatom2,natomcell,icart,icenter,i,j
      dipolfac=bohr ! for dipole moments in e*Angs 
      icell=1
      do iatom2=1,natoms
        do icenter=1,ncenter
          if (atoms(iatom2)%name.eq.centername(icenter)) then
            centeratom=atoms(iatom2)
            locdipol=0.0D0
            natomcell=0
            !write(13,3333) 'atoms in the cell centered at',
            write(33,3333) 'atoms in the cell centered at',
     &         centeratom%where,':'
 3333 format(A31,3(F10.3),A3)
          do iatom=1,natomsXL
            thisatom=atomsXL(iatom)
            cubatom=thisatom
            cubatom%where=nearestcubsite(thisatom%where,
     &            thisatom%geometry,
     &            cubcorners,cubcenters,cubfaces,ncubcorner,
     &            ncubcenter,ncubface)
            if (isincell(cubatom,centeratom,nuc)) then
                  !write(13,2222)thisatom%number,thisatom%name,
                  write(33,2222)thisatom%number,thisatom%name,
     &                          thisatom%where
c	          distvec=thisatom%where-centeratom%where
          do i=1,3
            distvec(i)=0.0D0
            do j=1,3
              distvec(i)=distvec(i)+(thisatom%where(j)
     &                           -centeratom%where(j))*avecs(j,i)
            end do
          end do
c                  locdipol=locdipol+thisatom%weight*dipolfac*alat
c     &			 *matmul(thisatom%charge,distvec)
                  locdipol=locdipol+thisatom%weight*dipolfac
     &           *matmul(thisatom%charge,distvec)
      dipol(icell,1:3)=centeratom%where(1:3)
      dipol(icell,4:6)=locdipol(1:3)
      natomcell=natomcell+1
      end if
      end do
      icell=icell+1
      !write(13,*) 'this cell contains',natomcell,'atoms.'
      write(33,*) 'this cell contains',natomcell,'atoms.'
      end if
      end do
      end do
 2222 format(I3,A3,3(F10.3))
      end subroutine dipole_center
      
      subroutine dipole_cube(natoms,natomsXL,atoms,atomsXL,avecs,
     &                    centername,ncenter,ncell,dipol,nuc,
     &                    cubcorners,cubcenters,cubfaces,
     &                    ncubcorner,ncubcenter,ncubface  )
	use defs
c	use functions
	implicit none
	integer, intent(in) :: natoms,natomsXL,ncenter,ncell,nuc(1:3)
	integer, intent(in) :: ncubcorner,ncubcenter,ncubface
	double precision, intent(in) :: avecs(1:3,1:3)
	double precision, intent(in) :: cubcorners(1:3,1:ncubcorner)
	double precision, intent(in) :: cubcenters(1:3,1:ncubcenter)
	double precision, intent(in) :: cubfaces(1:3,1:ncubface)
	double precision,intent(inout) :: dipol(1:ncell,1:6)
	character (len=2), intent(in) :: centername(1:ncenter)
c	type patom
c   	  sequence
c   	  integer:: number
c	  character:: name*2 
c	  integer:: direc 
c   	  double precision:: where(3)
c   	  double precision:: charge(1:3,1:3)
c   	  double precision:: weight
c   	  double precision:: distance
c	  character:: geometry*6
c	end type patom
	type(patom),intent(in) :: atoms(1:natoms),atomsXL(1:natomsXL)
	type(patom) centeratom,thisatom,cubatom
	double precision dipolfac,distvec(1:3),locdipol(1:3)
	integer icell,iatom,iatom2,natomcell,icart,icenter,i,j
	dipolfac=bohr ! for dipole moments in e*Angs 
	icell=1
	do iatom2=1,natoms
	  do icenter=1,ncenter	
	    if (atoms(iatom2)%name.eq.centername(icenter)) then
	      centeratom=atoms(iatom2)
	      locdipol=0.0D0
	      natomcell=0	
	      !write(13,3333) 'atoms in the cell centered at',
	      write(33,3333) 'atoms in the cell centered at',
     &         centeratom%where,':'
 3333 format(A31,3(F10.3),A3)
	      do iatom=1,natomsXL
	        thisatom=atomsXL(iatom)
		cubatom=thisatom
	 	cubatom%where=nearestcubsite(thisatom%where,
     &            thisatom%geometry,
     &            cubcorners,cubcenters,cubfaces,ncubcorner,
     &            ncubcenter,ncubface)
	        if (isincell(cubatom,centeratom,nuc)) then
                  !write(13,2222)thisatom%number,thisatom%name,
                  write(33,2222)thisatom%number,thisatom%name,
     &                          thisatom%where
c                 distvec=thisatom%where-cubatom%where
	          do i=1,3
		    distvec(i)=0.0D0
		    do j=1,3
		      distvec(i)=distvec(i)+(thisatom%where(j)
     &                           -cubatom%where(j))*avecs(j,i)
		    end do
		  end do
c                  locdipol=locdipol+thisatom%weight*dipolfac*alat
c     &			 *matmul(thisatom%charge,distvec)
                  locdipol=locdipol+thisatom%weight*dipolfac
     &			 *matmul(thisatom%charge,distvec)
	 	  dipol(icell,1:3)=centeratom%where(1:3)
	 	  dipol(icell,4:6)=locdipol(1:3)
		  natomcell=natomcell+1
		end if
	      end do
	      icell=icell+1
	      !write(13,*) 'this cell contains',natomcell,'atoms.'
	      write(33,*) 'this cell contains',natomcell,'atoms.'
	    end if
	  end do
	end do
 2222 format(I3,A3,3(F10.3))
	end subroutine
	
	subroutine rotate(atoms,natoms,pelements,nelements,
     &     rotangle,alternate,faces,nface,centers,ncenter,nuc,
     &     outunit,
     &           cubcorners,cubcenters,cubfaces,
     &           ncubcorner,ncubcenter,ncubface,avecs,print_xsf)
	use defs
c	use functions
	implicit none
	integer,intent(in) :: natoms,nelements,nface,ncenter,
     &                        alternate(1:3)
	integer, intent(in) :: nuc(1:3),outunit
	integer, intent(in) :: ncubcorner,ncubcenter,ncubface
	integer iele,iatom,iatom2,iface,icenter,nx,ny,nz
	double precision, intent(in) :: rotangle(1:3),avecs(1:3,1:3)
	double precision, intent(in) :: cubcorners(1:3,1:ncubcorner)
	double precision, intent(in) :: cubcenters(1:3,1:ncubcenter)
	double precision, intent(in) :: cubfaces(1:3,1:ncubface)
	double precision distvec(1:3),matrix(1:3,1:3),signum(1:3)
	double precision matrix_1(1:3,1:3),matrix_2(1:3,1:3),
     &                   matrix_3(1:3,1:3)
c	type patom
c   	  sequence
c   	  integer:: number
c	  character:: name*2 
c	  integer:: direc 
c   	  double precision:: where(3)
c   	  double precision:: charge(1:3,1:3)
c   	  double precision:: weight
c   	  double precision:: distance
c	  character:: geometry*6
c	end type patom
	type(patom), intent(in) :: atoms(natoms),faces(nface),
     &               centers(ncenter)
	logical,intent(in) :: print_xsf
	type(patom) :: atomsrot(natoms),thisatom,centeratom,cubatom
c	type pelement
c   	  sequence
c	  integer number
c	  integer:: howmany
c	  character:: name*2 
c   	  double precision:: charge(1:3,1:3)
c   	  double precision:: weight
c   	  double precision:: distance
c	end type pelement
	type(pelement)  :: pelements(nelements)
	! which rotation mode?
	write(outunit,*) '*****************************************'
	write(outunit,*) 'This is subroutine "ROTATE".'
	write(outunit,*) 'Parameters of rotation:'
	write(outunit,'(A12,3(F8.3))') ' angle (°): ', rotangle(1:3)
	if (abs(avecs(2,1)).ge.1E-6.or.abs(avecs(3,1)).ge.1E-6.or.
     &      abs(avecs(3,2)).ge.1E-6) then
		write(outunit,*)'Warning: your cell edges are not'
		write(outunit,*)'orthogonal, in this case the the'
		write(outunit,*)'"rotated" coordinates will be only'
		write(outunit,*)'approximately correct.'
		write(outunit,*)' '
	end if
	write(outunit,*)' '
	! calculate rotated coordinates
	atomsrot=atoms
	! calculate transformation matrix (+ or - sign comes later)
	! separate matrices for rotation around x,y,z axis
	matrix_1=0.0D0
	matrix_1(2,3)=-1.0D0
	matrix_1(3,2)=1.0D0
	matrix_1=matrix_1*tan(rotangle(1)*Pi/180.0D0)
	matrix_2=0.0D0
	matrix_2(1,3)=1.0D0
	matrix_2(3,1)=-1.0D0
	matrix_2=matrix_2*tan(rotangle(2)*Pi/180.0D0)
	matrix_3=0.0D0
	matrix_3(1,2)=-1.0D0
	matrix_3(2,1)=1.0D0
	matrix_3=matrix_3*tan(rotangle(3)*Pi/180.0D0)
	write(outunit,*) 'transformation matrix:'
	write(outunit,'(9(F10.6))') matrix_1(1,1:3),matrix_1(2,1:3),
     &                             matrix_1(3,1:3)
	! loop over atoms
	do iatom=1,natoms
	  thisatom=atomsrot(iatom)
	  do iface=1,nface,3
	    if (thisatom%name==faces(iface)%name) then  ! does the atom belong to an octahedron?	
	      write(outunit,'(A34,A2,3(F10.3))')
     &            'octahedral atom (name, position): ',
     &              thisatom%name,thisatom%where
	      ! find a neighboring center atom
	      do iatom2=1,natoms
		do icenter=1,ncenter
		  if(atoms(iatom2)%name.eq.centers(icenter)%name)then
		      cubatom=thisatom
	 	      cubatom%where=nearestcubsite(thisatom%where,
     &                  thisatom%geometry,
     &                  cubcorners,cubcenters,cubfaces,ncubcorner,
     &                  ncubcenter,ncubface)
                    if (isincell(cubatom,atoms(iatom2),nuc)) then
		      centeratom=atoms(iatom2)
		      ! calculate sign of rotation
		      nx=nint(centeratom%where(1)*dble(nuc(1))+0.5D0)
		      ny=nint(centeratom%where(2)*dble(nuc(2))+0.5D0)
		      nz=nint(centeratom%where(3)*dble(nuc(3))+0.5D0)
		      signum(1)=(-1.0D0)**(alternate(1)*(nx-1)+ny+nz)
		      signum(2)=(-1.0D0)**(alternate(2)*(ny-1)+nx+nz)
		      signum(3)=(-1.0D0)**(alternate(3)*(nz-1)+nx+ny)
		      write(outunit,'(A25,A3,3(F10.6),3(I3),3(F8.3))') 
     &                'neighboring center atom: ',
     &                centeratom%name,centeratom%where,nx,ny,nz,
     &                signum(1:3)
		    end if  
		  end if  
		end do
	      end do	  
	      ! perform rotation
	      matrix=matrix_1*signum(1)+matrix_2*signum(2)
     &         +matrix_3*signum(3)
	      distvec=thisatom%where-centeratom%where
	      atomsrot(iatom)%where=thisatom%where
     &          +matmul(matrix,distvec)
	      write(outunit,*) ' '
	    end if
	  end do  
	end do
	! print resulting coordinates
	!open(18,file='COORAT.ROT',status='replace')
	open(38,file='COORAT.ROT',status='replace')
	do iele=1,nelements
	  !write(18,'(A7,A2,A14,I4)')'  name=',pelements(iele)%name,
	  write(38,'(A7,A2,A14,I4)')'  name=',pelements(iele)%name,
     &          '  natom=', pelements(iele)%howmany 
	  do iatom=1,natoms
	    if (atomsrot(iatom)%name==pelements(iele)%name) then
	      !write(18,'(3(F15.10))') atomsrot(iatom)%where
	      write(38,'(3(F15.10))') atomsrot(iatom)%where
	    end if	
	  end do
	end do
	!close(18)
	close(38)
	if (print_xsf) call write_xsf('GEOMETRY.ROT.xsf',16,
     &   atomsrot,natoms,avecs)
	end subroutine

	subroutine AFE_DISP(atoms,natoms,pelements,nelements,afemod,
     &     afevec,afecenteramp,afefaceamp,alternate,faces,
     &     nface,centers,ncenter,nuc,outunit,
     &                    cubcorners,cubcenters,cubfaces,
     &                    ncubcorner,ncubcenter,ncubface,
     &                    print_xsf,avecs)
	use defs
	implicit none
	integer,intent(in) :: natoms,nelements,nface,ncenter,
     &     alternate(1:3)
	integer, intent(in) :: nuc(1:3),outunit
	integer, intent(in) :: ncubcorner,ncubcenter,ncubface
	integer iele,iatom,iatom2,iface,icenter,nx,ny,nz
	double precision, intent(in) :: afecenteramp,afefaceamp,
     &                                  afevec(1:3),avecs(1:3,1:3)
	double precision, intent(in) :: cubcorners(1:3,1:ncubcorner)
	double precision, intent(in) :: cubcenters(1:3,1:ncubcenter)
	double precision, intent(in) :: cubfaces(1:3,1:ncubface)
	double precision distvec(1:3),matrix(1:3,1:3),signum
	character, intent(in) :: afemod*3
	type(patom), intent(in) :: atoms(natoms),faces(nface),
     &               centers(ncenter)
	logical, intent(in) :: print_xsf
	type(patom) :: atomsafe(natoms),thisatom,centeratom,cubatom
	type(pelement)  :: pelements(nelements)
	! which AFE mode?
	write(outunit,*) '*****************************************'
	write(outunit,*) 'This is subroutine "AFE".'
	write(outunit,*) 'Parameters of AFE disps:'
	write(outunit,'(A11,A3)') ' AFE mode: ', afemod
	write(outunit,*) 
     &  ' displacement amplitude (in single cell lattice constants): '
	write(outunit,*) ' center atoms: ',afecenteramp 
	write(outunit,*) ' face atoms:   ',afefaceamp 
	write(outunit,'(A21,3(F8.3))') ' direction of disps: ',afevec
	write(outunit,'(A31,3(I8))') 
     &       ' alternating direction (x,y,z): ',alternate(1:3)
	write(outunit,*)' '
	! calculate AFE coordinates
	atomsafe=atoms
	! loop over atoms
	do iatom=1,natoms
	  thisatom=atomsafe(iatom)
          !
	  ! begin only shift center atoms 
          !
	  if (thisatom%geometry.eq.'center') then
		      nx=nint(thisatom%where(1)*dble(nuc(1))+0.5D0)
		      ny=nint(thisatom%where(2)*dble(nuc(2))+0.5D0)
		      nz=nint(thisatom%where(3)*dble(nuc(3))+0.5D0)
		      signum=(-1.0D0)**(alternate(1)*(nx-1)
     &                  +alternate(2)*(ny-1)+alternate(3)*(nz-1))
		      write(outunit,*) 
     &               'center atom, x y z, nx,ny,nz,signum,alternate: ',
     &                thisatom%name,thisatom%where,nx,ny,nz,signum,
     &                alternate(1:3)
	      ! displace atoms
	      atomsafe(iatom)%where=thisatom%where
     &          +afevec*signum*afecenteramp   
	      !write(outunit,*) ' '
	  end if
          !
	  ! end only shift center atoms 
          !
	  !end do  
	end do
	! print resulting coordinates
	open(38,file='COORAT.AFE',status='replace')
	do iele=1,nelements
	  !write(18,'(A7,A2,A14,I4)')'  name=',pelements(iele)%name,
	  write(38,'(A7,A2,A14,I4)')'  name=',pelements(iele)%name,
     &          '  natom=', pelements(iele)%howmany 
	  do iatom=1,natoms
	    if (atomsafe(iatom)%name==pelements(iele)%name) then
	      !write(18,'(3(F15.10))') atomsafe(iatom)%where
	      write(38,'(3(F15.10))') atomsafe(iatom)%where
	    end if	
	  end do
	end do
	!close(18)
	close(38)
	if (print_xsf) call write_xsf('GEOMETRY.AFE.xsf',16,atomsafe,
     &                                natoms,avecs)
	end subroutine

c---------------------------------------------------------------------

	subroutine write_xsf(fname,fnamelen,atoms,natoms,avecs)
	
        use defs
        implicit none
	type(patom), intent(in) :: atoms(1:natoms)
	integer, intent(in) :: fnamelen,natoms
	character (len=*), intent(in) :: fname
	double precision, intent(in) :: avecs(1:3,1:3)
	double precision coords(1:3)
	integer i,j,icart
c	avecs=0.0D0
c	do i=1,3
c	  avecs(i,i)=alat(i)
c	end do
	!open(10,file=fname,status='replace')
	open(30,file=fname,status='replace')
	!write(10,*) ' CRYSTAL'
	!write(10,*) ' PRIMVEC'
	write(30,*) ' CRYSTAL'
	write(30,*) ' PRIMVEC'
	do i=1,3
	  !write(10,'(3(F15.10))') avecs(i,1:3)
	  write(30,'(3(F15.10))') avecs(i,1:3)
	end do
	!write(10,*) 'PRIMCOORD'
	!write(10,*) natoms, '  1'
	write(30,*) 'PRIMCOORD'
	write(30,*) natoms, '  1'
	do i=1,natoms
	  do icart=1,3
	    coords(icart)=0.0D0
	    do j=1,3
	      coords(icart)=coords(icart)+atoms(i)%where(j)*avecs(j,icart)
	    end do
	  end do
	  !write(10,'(A5,3(F15.10))') atoms(i)%name,coords
	  write(30,'(A5,3(F15.10))') atoms(i)%name,coords
	end do
	!close(10)
	close(30)
	end subroutine

      integer function direction(thisatom,centeratom)
      use defs
      implicit none
c	type patom
c   	  sequence
c   	  integer:: number
c	  character:: name*2 
c	  integer:: direc 
c   	  double precision:: where(3)
c   	  double precision:: charge(1:3,1:3)
c   	  double precision:: weight
c   	  double precision:: distance
c	  character:: geometry*6
c	end type patom
	type(patom):: thisatom,centeratom
	double precision x_this(1:3),x_center(1:3),distvec(1:3),largest
	x_this(1:3)=thisatom%where(1:3)
	x_center(1:3)=centeratom%where(1:3)
	distvec=x_center-x_this
	direction=1
	largest=abs(distvec(1))
	if (abs(distvec(2)).gt.largest) then
		largest=abs(distvec(2))
		direction=2
	end if
	if (abs(distvec(3)).gt.largest) then
		largest=abs(distvec(3))
		direction=3
	end if
	end function
		

      end module
