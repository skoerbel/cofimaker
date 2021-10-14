      module for_lammps

      contains

c---------------------------------------------------------------------
      subroutine lammps_xyz2sxyz()
      ! converts lammps xyz output file to Adham's special 
      ! xyz (sxyz) format
      use defs
      implicit none

	double precision a1,a2,b1,b2,c1,c2
	double precision x,y,z,mass
	integer i,j,natoms,ntypes,itype,imass,icomment
        integer iele,nele,numele
	integer, allocatable :: typenums(:)
	character finname*30,flinname*30,foutname*30,line*100
        character typename*2
	character,allocatable :: typenames(:)*2
        logical commented,newele
        type(element), allocatable :: eles(:)    
        type(atom), allocatable :: atoms(:)

        write(12,fsubstart)trim(adjustl('lammps_xyz2sxyz'))
        ! read input and output filenames, number of atom types, names of atom
        ! types from INPUT.COFIMA
        read(10,*) finname
        read(10,*) flinname
        read(10,*) foutname
	read(10,*) ntypes
	allocate(typenames(ntypes),typenums(ntypes))
	do i=1,ntypes
	  read(10,*) itype,typenames(itype)
	  do j=1,100
	    if (elements(j).eq.typenames(itype))
     &        typenums(itype)=j
	  end do
	end do
        ! read header of lammps xyz output file
	open(60,file=finname,status="old")
	open(61,file=foutname,status="replace")
	do i=1,3
	  read(60,'(A100)') line
	end do
	read(60,*) natoms
	read(60,'(A100)') line
	read(60,*) a1,a2
	read(60,*) b1,b2
	read(60,*) c1,c2
	write(61,*) natoms 
	write(61,*) "# Adham's special xyz file format"
	read(60,'(A100)') line
        ! read and write atom coordinates
        allocate(atoms(natoms))
	do i=1,natoms
	  !read(60,*) itype,x,y,z
          read(60,*) itype,atoms(i)%where !x,y,z
          atoms(i)%name='       '
          atoms(i)%name(1:2)=elements(typenums(itype))
          !write(61,'(A4,3(F12.6))') elements(typenums(itype)),x,y,z
          !write(61,'(A4,3(F12.6))') elements(typenums(itype)),
!     &       atoms(i)%where
	end do
      ! count elements and number of atoms of each element
      !atoms(1:natoms)%written=.false.
      nele=0
      do i=1,natoms
            newele=.true.
            do j=1,i-1
              if (atoms(i)%name(1:8).eq.atoms(j)%name(1:8))
     &         newele=.false.
            end do
            if (newele) nele=nele+1
      end do
      write(12,'(8x,"Found ",I4," elements.")')nele
      allocate(eles(1:nele))
      ! get names of elements
      iele=1
      do i=1,natoms
            newele=.true.
            do j=1,i-1
              if (atoms(i)%name(1:8).eq.atoms(j)%name(1:8))
     &         newele=.false.
            end do
            if (newele) then
                  eles(iele)%name=atoms(i)%name
                  iele=iele+1
            end if
      end do
      ! get numbers of elements
      do iele=1,nele
            numele=0
            do j=1,natoms
                  if(atoms(j)%name.eq.eles(iele)%name) numele=numele+1
            end do
            eles(iele)%howmany=numele
      end do
      ! write coordinates
      do iele=1,nele
            do i=1,natoms
      	          if (atoms(i)%name==eles(iele)%name) 
     &            write(61,'(A4,3(F12.6))') atoms(i)%name,
     &            atoms(i)%where
            end do
      	end do
        ! write end of sxyz file
	write(61,*) 'alat'
	write(61,*) '1.000'
	write(61,*) 'supercell'
	write(61,'(3(F10.6))') a2-a1,0.0,0.0
	write(61,'(3(F10.6))') 0.0,b2-b1,0.0
	write(61,'(3(F10.6))') 0.0,0.0,c2-c1
        ! scan lammps input file for masses 
        open(62,file=flinname,status='old')
        rewind(62)
10      read(62,'(A100)',end=15) line
        ! look for mass command 
        if(index(line,'mass').ge.1) then
            ! if the line is not commented, write mass to output file
            commented=.false.
            imass=0
            do while (line(imass:imass+3).ne.'mass')
              imass=imass+1
            end do
            icomment=0
            if (index(line,'#').ge.1) then
              do while (line(icomment:icomment).ne.'#')
                icomment=icomment+1
              end do
              if (icomment.lt.imass) commented=.true.
            end if  
            if (.not.commented) then
              read(line(imass+4:100),*) itype,mass
              write(61,'(A4,A6,F12.6)') 'mass',
     &               elements(typenums(itype)),mass
            end if  
        end if
        goto 10
15      continue
        close(60)
        close(62)

	write(61,*) 'reduced_coordinates'
	
	close(61)
      write(12,fsubend)
      end subroutine

c---------------------------------------------------------------------

      subroutine data2cfg()
      ! converts lammps structure input file (data) to cfg file
      use misc
      use defs
      implicit none

	double precision a1,a2,b1,b2,c1,c2,rdum
	double precision abscoords(1:3),mass,vecs(1:3,1:3)
	integer i,j,natoms,ntypes,itype,imass,icomment
        integer iele,nele,numele,idum
	integer, allocatable :: typenums(:)
	character finname*30,foutname*30,line*100
        character typename*2
	character,allocatable :: typenames(:)*2
        logical commented,newele
        type(element), allocatable :: eles(:)    
        type(atom), allocatable :: atoms(:)
        !character FMT1*1024,FMT2*1024

        write(12,fsubstart)trim(adjustl('data2cfg'))
        ! read input and output filenames, number of atom types, names of atom
        ! types from INPUT.COFIMA
        read(10,*,err=500) finname
        read(10,*,err=500) foutname
	read(10,*,err=500) ntypes
	allocate(typenames(ntypes),typenums(ntypes))
	do i=1,ntypes
	  read(10,*,err=500) itype,typenames(itype)
	  do j=1,100
	    if (elements(j).eq.typenames(itype))
     &        typenums(itype)=j
	  end do
	end do
        ! read header of lammps data file
	open(60,file=finname,status="old",err=600)
	open(61,file=foutname,status="replace",err=610)
	!do i=1,3
	!  read(60,'(A100)') line
	!end do
	read(60,'(A100)',err=600) line
	read(60,*,err=600) natoms
	read(60,'(A100)',err=600) line
	read(60,'(A100)',err=600) line
	read(60,'(A100)',err=600) line
	read(60,*,err=600) a1,a2
	read(60,*,err=600) b1,b2
	read(60,*,err=600) c1,c2
        vecs=0.0D0
        vecs(1,1)=a2-a1
        vecs(2,2)=b2-b1
        vecs(3,3)=c2-c1
!	write(61,'("Number of particles = ",I)') natoms 
        WRITE(FMT2,*) natoms
        WRITE(FMT1,*) '("Number of particles = ",I',
     &   len(trim(FMT2)),')'
        write(61,FMT1) natoms 
!	write(61,'("Number of particles = ",I)') natoms 
	write(61,'("A = 1 Angstrom (basic length-scale)")') 
	write(61,'("H0(1,1) = ",F10.5," A")') vecs(1,1)
	write(61,'("H0(1,2) = ",F10.5," A")') vecs(1,2)
	write(61,'("H0(1,3) = ",F10.5," A")') vecs(1,3)
        write(61,'("H0(2,1) = ",F10.5," A")') vecs(2,1)
	write(61,'("H0(2,2) = ",F10.5," A")') vecs(2,2)
	write(61,'("H0(2,3) = ",F10.5," A")') vecs(2,3)
	write(61,'("H0(3,1) = ",F10.5," A")') vecs(3,1)
	write(61,'("H0(3,2) = ",F10.5," A")') vecs(3,2)
	write(61,'("H0(3,3) = ",F10.5," A")') vecs(3,3)
	write(61,'(".NO_VELOCITY.")') 
	write(61,'("entry_count = 3")') 
        read(60,'(A100)',err=600) line
	read(60,'(A100)',err=600) line
	read(60,'(A100)',err=600) line
        ! read and write atom coordinates
        allocate(atoms(natoms))
	do i=1,natoms
          read(60,*,err=600) idum,itype,rdum,abscoords !absolute coords
          call abs2frac(abscoords,vecs,atoms(i)%where)
          atoms(i)%name='       '
          atoms(i)%name(1:2)=elements(typenums(itype))
          !write(61,'(A4,3(F12.6))') elements(typenums(itype)),x,y,z
          !write(61,'(A4,3(F12.6))') elements(typenums(itype)),
!     &       atoms(i)%where
	end do
      ! count elements and number of atoms of each element
      !atoms(1:natoms)%written=.false.
      nele=0
      do i=1,natoms
            newele=.true.
            do j=1,i-1
              if (atoms(i)%name(1:8).eq.atoms(j)%name(1:8))
     &         newele=.false.
            end do
            if (newele) nele=nele+1
      end do
      write(12,'(8x,"Found ",I4," elements.")')nele
      allocate(eles(1:nele))
      ! get names of elements
      iele=1
      do i=1,natoms
            newele=.true.
            do j=1,i-1
              if (atoms(i)%name(1:8).eq.atoms(j)%name(1:8))
     &         newele=.false.
            end do
            if (newele) then
                  eles(iele)%name=atoms(i)%name
                  iele=iele+1
            end if
      end do
      ! get numbers of elements
      do iele=1,nele
            numele=0
            do j=1,natoms
                  if(atoms(j)%name.eq.eles(iele)%name) numele=numele+1
            end do
            eles(iele)%howmany=numele
      end do
      ! write coordinates
      do iele=1,nele
            write(61,'("1.00")') 
            write(61,'(A8)') eles(iele)%name
            do i=1,natoms
      	          if (atoms(i)%name==eles(iele)%name) 
     &            write(61,'(3(F12.6))') atoms(i)%where
            end do
      	end do
!        ! write end of sxyz file
!	write(61,*) 'alat'
!	write(61,*) '1.000'
!	write(61,*) 'supercell'
!	write(61,'(3(F10.6))') a2-a1,0.0,0.0
!	write(61,'(3(F10.6))') 0.0,b2-b1,0.0
!	write(61,'(3(F10.6))') 0.0,0.0,c2-c1
!        ! scan lammps input file for masses 
!        open(62,file=flinname,status='old')
!        rewind(62)
!10      read(62,'(A100)',end=15) line
!        ! look for mass command 
!        if(index(line,'mass').ge.1) then
!            ! if the line is not commented, write mass to output file
!            commented=.false.
!            imass=0
!            do while (line(imass:imass+3).ne.'mass')
!              imass=imass+1
!            end do
!            icomment=0
!            if (index(line,'#').ge.1) then
!              do while (line(icomment:icomment).ne.'#')
!                icomment=icomment+1
!              end do
!              if (icomment.lt.imass) commented=.true.
!            end if  
!            if (.not.commented) then
!              read(line(imass+4:100),*) itype,mass
!              write(61,'(A4,A6,F12.6)') 'mass',
!     &               elements(typenums(itype)),mass
!            end if  
!        end if
!        goto 10
15      continue
        close(60)
!        close(62)

!	write(61,*) 'reduced_coordinates'
!	
	close(61)
      write(12,fsubend)
      return

500   print*,"Error: Something is wrong with INPUT.COFIMA!"
      write(12,'(8x,"Error: Something is wrong with INPUT.COFIMA!")')
      close(60)
      close(61)
      nerr=nerr+1
      return
600   print*,"Error: Something is wrong with the data file."
      write(12,'(8x,"Error: Something is wrong with the data file.")')
      close(60)
      close(61)
      nerr=nerr+1
      return
610   print*,"Error: Can't write the output file!"
      write(12,'(8x,"Error: Cannot write the output file!")')
      close(60)
      close(61)
      nerr=nerr+1
      return

      end subroutine


      end module
