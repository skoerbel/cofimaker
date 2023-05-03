      module readcoords

      implicit none

      contains

c---------------------------------------------------------------------
      
      subroutine read_coords(fincoords,finform,atoms,natoms,species,
     &                       nspecies,vecs)

      use defs
      implicit none
      
      ! variables
      character(len=*), intent(in) :: fincoords,finform
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms,nspecies
      character line*100,nameat*2
      double precision x(1:3),vecs(1:3,1:3)
      
      ! internal variables
      logical isopen12,fileexists
      integer natom,i,j,iatom,ispecies
      character errmssg*256
      
      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim('read_coords')
      else
            print fsubstart, trim(adjustl('read_coords'))
      end if

      ! check if the specified file exists
      inquire(file=fincoords, exist=fileexists)
      if (.not.fileexists) goto 1001

      ! according to specified format, choose routine to read
      ! coordinates
      select case (finform)
      case("ascii","ASCII")
      call read_ascii(fincoords,atoms,natoms,species,nspecies,
     &                    vecs)
      case("cfg","CFG")
      call read_cfg(fincoords,atoms,natoms,species,nspecies,
     &                    vecs)
      case("coorat","COORAT")
            call read_coorat(fincoords,atoms,natoms,species,nspecies,
     &                       vecs)
      case("data","DATA")
      call read_data(fincoords,atoms,natoms,species,nspecies,
     &                    vecs)
      case("gin","GIN")
      call read_gin(fincoords,atoms,natoms,species,nspecies,
     &                    vecs)
      case("lxyz","LXYZ")
      call read_lxyz(fincoords,atoms,natoms,species,nspecies,
     &                    vecs)
      case("out","OUT","mbpp","MBPP")
            call read_out(fincoords,atoms,natoms,species,nspecies,
     &                    vecs)
      case("xsf","XSF")
            call read_xsf(fincoords,atoms,natoms,species,nspecies,
     &                    vecs)
      case("xsf_forces","XSF_FORCES")
            call read_xsf_forces(fincoords,atoms,natoms,species,          &
     &                    nspecies,vecs)
      case("abinit","ABINIT")
            call read_abinit(fincoords,atoms,natoms,species,nspecies,
     &                    vecs)
      case("nwchem","NWCHEM")
            call read_nwchem(fincoords,atoms,natoms,species,nspecies,
     &                    vecs)
      case("xyz","XYZ")
            call read_xyz(fincoords,atoms,natoms,species,nspecies,
     &                    vecs)
      case("outcar","OUTCAR")
            call read_outcar(fincoords,atoms,natoms,species,nspecies,
     &                    vecs)
      case("poscar","POSCAR")
            call read_poscar(fincoords,atoms,natoms,species,nspecies,
     &                    vecs)
      case("poscar_bec","POSCAR_BEC")
            call read_poscar_bec(fincoords,atoms,natoms,species,          &
     &                           nspecies,vecs)
      case("poscar_forces","POSCAR_FORCES")
            call read_poscar_forces(fincoords,atoms,natoms,species,       &
     &                           nspecies,vecs)
      case("poscar_flags","POSCAR_FLAGS")
            call read_poscar_flags(fincoords,atoms,natoms,species,        &
     &                           nspecies,vecs)
      case("siesta","SIESTA")
            call read_siesta(fincoords,atoms,natoms,species,nspecies,     &
     &                    vecs)
      case("siesta_out","SIESTA_OUT")
            call read_siesta_out(fincoords,atoms,natoms,species,          &
     &     nspecies, vecs)
      case("cif","CIF")
            call read_cif(fincoords,atoms,natoms,species,nspecies,
     &                    vecs)
      case("vesta","VESTA")
            call read_vesta(fincoords,atoms,natoms,species,nspecies,
     &                    vecs)
      case("findsym","FINDSYM")
            call read_findsym(fincoords,atoms,natoms,species,nspecies,
     &                    vecs)
      case default
            goto 1000
!            if(isopen12) then
!                  write(12,ferrmssg) 
!     &             trim('Input coordinate file has unknown format')
!            else
!                  print ferrmssg, 
!     &             trim('Input coordinate file has unknown format')
!            end if
      end select
      !if(nerr.gt.0) goto 1002
      
      ! write some output ...
      if(isopen12) then
            write(12,'(8x,A)')'Species found in the file (name, howmany,
     & mass):' 
            do i=1,nspecies
!              write(12,'(8x,A2,1x,I,1x,F10.6)') species(i)%name,
!     &        species(i)%howmany,species(i)%mass         
              WRITE(FMT2,*) species(i)%howmany
              WRITE(FMT1,*) '(8x,A2,1x,I',len(trim(FMT2)),
     &          ',1x,F10.6)'
              if(isopen12) then
                write(12,FMT1)species(i)%name,species(i)%howmany,
     &            species(i)%mass
              else
                print FMT1,species(i)%name,species(i)%howmany,
     &            species(i)%mass
              end if
            end do
      else
            print '(8x,A)','Species found in the file (name, howmany, ma
     &ss):' 
            do i=1,nspecies
!              print '(8x,A2,1x,I,1x,F10.6)', species(i)%name,
!     &        species(i)%howmany,species(i)%mass         
              WRITE(FMT2,*) species(i)%howmany
              WRITE(FMT1,*) '(8x,A2,1x,I',len(trim(FMT2)),
     &          ',1x,F10.6)'
              print FMT1,species(i)%name,species(i)%howmany,
     &        species(i)%mass
            end do
      end if
                  
      if(isopen12) then
            write(12,fsubendext) 'read_coords'
      else
            print fsubendext, 'read_coords'
      end if
      return

      ! error handling
1000  if(isopen12) then
            write(12,ferrmssg) 
     &       trim('Input coordinate file has unknown format')
      else
            print ferrmssg, 
     &       trim('Input coordinate file has unknown format')
      end if
      nerr=nerr+1
      return
1001  nerr=nerr+1
      if(isopen12) then
            errmssg=''
            write(errmssg,'("file not found: ",A100)') fincoords
            write(12,ferrmssg) errmssg
            write(12,fsubendext) 'read_coords'
      else
            errmssg=''
            write(errmssg,'("file not found: ",A100)') fincoords
            print ferrmssg,errmssg
            print fsubendext, 'read_coords'
      end if
      return
1002  nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "The specified file could not be read" 
            write(12,fsubendext) 'read_coords'
      else
            print ferrmssg,"The specified file could not be read" 
            print fsubendext, 'read_coords'
      end if
      return
      
      end subroutine

c---------------------------------------------------------------------

      subroutine read_ascii(filename,atoms,natoms,species,nspecies,
     &                        vecs)
      ! read coordinates from ASE ascii file
      
      use defs
      use misc, only : abs2frac,getspecies

      implicit none
      
      ! variables
      character (LEN=*), intent(in) :: filename
      double precision vecs(1:3,1:3)
      integer natoms
      integer nspecies
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)

      ! local variables
      integer i,nlines,length
      logical isopen12
      character line*256

      INQUIRE (unit=12, opened=isopen12)
      if(talk) then
        if(isopen12) then
            write(12,fsubstart) trim(adjustl('read_ascii'))
        else
            print fsubstart, trim(adjustl('read_ascii'))
        end if
      end if
      !
      ! begin read ascii file
      open(51,file=filename,status='old',err=1000)
      ! determine number of meaningful lines
      nlines=0
      length=1
      do while ( length>0 )
        read(51,'(A256)',end=10) line
        length=len_trim(line)
        if (length>0) nlines=nlines+1
      end do
10    if (talk) print '(8x,"Found ",I0," lines")',nlines
      natoms=nlines-3
      allocate(atoms(natoms))
      rewind(51)
      read(51,'(A256)') line ! Title
      !
      ! begin read cell dimensions
      vecs=0.0d0
      read(51,*) vecs(1,1),vecs(2,1),vecs(2,2)
      read(51,*) vecs(3,1),vecs(3,2),vecs(3,3)
      ! end read cell dimensions
      !
      ! begin atom information
      atoms(:)%name=" "
      atoms(:)%occ=1.0d0
      atoms(:)%mass=0.0d0
      do i=1,natoms
        read(51,*) atoms(i)%abswhere,atoms(i)%name
        call abs2frac(atoms(i)%abswhere,vecs,atoms(i)%where)
        atoms(i)%core="core"
      end do
      close(51)
      call getspecies(atoms,species)
      nspecies=size(species)
      if(talk) then
        if(isopen12) then
            write(12,fsubendext),'read_ascii'
        else
            print fsubendext, 'read_ascii'
        end if
      end if
      return

1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Something went wrong during reading" 
      else
            print ferrmssg,"Something wrong during reading" 
      end if
      close(51)
      stop

      end subroutine read_ascii

c---------------------------------------------------------------------

      subroutine read_coorat(fincoords,atoms,natoms,species,nspecies,
     &                       vecs)

      use defs
      use misc, only : frac2abs,get_masses_of_atoms,                      &
     &                 get_masses_of_species
      implicit none
      
      ! variables
      character(len=*), intent(in) :: fincoords
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms,nspecies
      character line*100,nameat*2
      double precision x(1:3),vecs(1:3,1:3)
      
      ! internal variables
      logical isopen12
      integer natom,i,j,iatom,ispecies
      
      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim(adjustl('read_coorat'))
      else
            print fsubstart, trim(adjustl('read_coorat'))
      end if
     
      !! fake lattice vectors
      !vecs=0.0D0
      !vecs(1,1)=1.0D0
      !vecs(2,2)=1.0D0
      !vecs(3,3)=1.0D0
      ! read COORAT file
      open(21,file=fincoords,status='old',err=1000)
      ! get number of species and number of atoms
      natoms=0
      nspecies=0
10    read(21,'(A100)',end=15,err=1000) line
      line=adjustl(line)
      if (line(1:1).ne.'!') then
            if(index(line,'natom').le.0)natoms=natoms+1
            if(index(line,'natom').gt.0)nspecies=nspecies+1
      end if
      goto 10
15    continue
      if (isopen12) then
!        write(12,'(8x,"Number of elements and atoms found: ",I,I)'),
!     &        nspecies,natoms  
        WRITE(FMT2,*) nspecies
        WRITE(FMT3,*) natoms
        WRITE(FMT1,*) '(8x,"Number of elements and atoms found: ",I',
     &   len(trim(FMT2)),',I',len(trim(FMT3)),')'
        write(12,FMT1)  nspecies,natoms 
      else
!        print '(8x,"Number of elements and atoms found: ",I,I)',
!     &        nspecies,natoms  
        WRITE(FMT2,*) nspecies
        WRITE(FMT3,*) natoms
        WRITE(FMT1,*) '(8x,"Number of elements and atoms found: ",I',
     &   len(trim(FMT2)),',I',len(trim(FMT3)),')'
        print FMT1, nspecies,natoms 
      end if
      rewind(21)
      if (.not. allocated(species)) then
            allocate(species(nspecies))
            species%name=" "
      end if
      if (.not. allocated(atoms)) then
            allocate(atoms(natoms))
            atoms%name=" "
      end if
      ispecies=1
      iatom=1
20    read(21,'(A100)',end=25,err=1000) line
      !get number of atoms of each species
      line=adjustl(line)
      if (line(1:1).ne.'!') then
        if (index(line,'natom').gt.0) then
            i=index(line,'natom')+5
            do while (line(i:i).eq.' ') 
                  i=i+1
            end do
            do while (line(i:i).eq.'=') 
                  i=i+1
            end do
            do while (line(i:i).eq.' ') 
                  i=i+1
            end do
            j=i
            do while (line(j:j).ne.' ')
                  j=j+1
            end do
!            read(line(i:j),'(i)') natom
            read(line(i:j),*) natom
            species(ispecies)%howmany=natom
        end if
        !get name of each species
        if (index(line,'name').gt.0) then
            i=index(line,'name')+4
            do while (line(i:i).eq.' ') 
                  i=i+1
            end do
            do while (line(i:i).eq.'=') 
                  i=i+1
            end do
            do while (line(i:i).eq.' ') 
                  i=i+1
            end do
            nameat=line(i:i+1)
            species(ispecies)%name=nameat
            species(ispecies)%mass=0.0D0    ! fake mass
            ispecies=ispecies+1
        end if
        i=1
        do while (i.le.natom)
          read(21,'(A100)',end=23) line
          if(index(line,'natom').ge.1) goto 23
          line=adjustl(line)
          if(line(1:1).ne.'!') then
            read(line,*) x(1:3)
            atoms(iatom)%where=x
            call frac2abs(x,vecs,atoms(iatom)%abswhere)
            atoms(iatom)%name=nameat
            atoms(iatom)%core="core"   ! treat as core
            atoms(iatom)%charge=0.0D0  ! fake charge
            atoms(iatom)%occ=1.0D0  ! default occupation
            iatom=iatom+1
            i=i+1
          end if
        end do
      end if
      goto 20

23    continue
      if(isopen12) then
            write(12,ferrmssg)'coordinates section ended' 
      else
            print ferrmssg,'coordinates section ended unexpectedly'
      end if
      nerr=nerr+1
25    continue
      call get_masses_of_species(species)
      call get_masses_of_atoms(atoms)
      if(isopen12) then
            write(12,'(8x,"Species found in the file (name, howmany, mas
     &s):")') 
            do i=1,nspecies
!              write(12,'(8x,A2,1x,I,1x,F10.6)') species(i)%name,
!     &        species(i)%howmany,species(i)%mass         
              WRITE(FMT2,*) species(i)%howmany
              WRITE(FMT1,*) '(8x,A2,1x,I',len(trim(FMT2)),
     &          ',1x,F10.6)'
              write(12,FMT1)species(i)%name,species(i)%howmany,
     &        species(i)%mass
            end do
            !write(12,'(8x,A)')'Atoms found in the file:' 
            !do i=1,natoms
            !  write(12,*) atoms(i)
            !end do
      else
            print '(8x,"Species found in the file (name, howmany, mass):
     &")' 
            do i=1,nspecies
!              print '(8x,A2,1x,I,1x,F10.6)', species(i)%name,
!     &        species(i)%howmany,species(i)%mass         
              WRITE(FMT2,*) species(i)%howmany
              WRITE(FMT1,*) '(8x,A2,1x,I',len(trim(FMT2)),
     &          ',1x,F10.6)'
              print FMT1,species(i)%name,species(i)%howmany,
     &        species(i)%mass
            end do
            !print '(8x,A)','Atoms found in the file:' 
            !do i=1,natoms
            !  print *, atoms(i)
            !end do
      end if
      close(21)
      close(22)
      if(isopen12) then
            write(12,fsubendext) 'read_coorat'
      else
            print fsubendext, 'read_coorat'
      end if
      return
      
1000  nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Something wrong with COORAT file." 
            write(12,fsubendext) 'read_coorat'
      else
            print ferrmssg,"Something wrong with COORAT file." 
            print fsubendext, 'read_coorat'
      end if
      close(21)
      return

      end subroutine

c---------------------------------------------------------------------

      subroutine read_xsf(fincoords,atoms,natoms,species,nspecies,
     &           vecs)

      use defs
      use misc
      implicit none
      
      ! variables
      character(len=*), intent(in) :: fincoords
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms,nspecies
      character line*100,nameat*2
      double precision x(1:3),vecs(1:3,1:3)
      
      ! internal variables
      logical isopen12,newspecies
      integer natom,i,j,iatom,iatom2,ispecies
      integer iele
      character elechar*40,line2*256
      !character FMT1*1024,FMT2*1024,FMT3*1024 
      
      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim(adjustl('read_xsf'))
      else
            print fsubstart, trim(adjustl('read_xsf'))
      end if
      
      ! read xsf file
      open(21,file=fincoords,status='old',err=999)
      ! get lattice vectors
10    read(21,'(A100)',end=1000,err=1000) line
      line=adjustl(line)
      if (line(1:1).ne.'!'.and.line(1:1).ne.'#') then
            if(index(line,'PRIMVEC').gt.0) then
                  do i=1,3
                    read(21,*) vecs(i,1:3)
                  end do
                  goto 12
            end if
      end if
      goto 10

      ! get number of atoms
      natoms=0
      !nspecies=0
12    read(21,'(A100)',end=15,err=1000) line
      line=adjustl(line)
      if (line(1:1).ne.'!'.and.line(1:1).ne.'#') then
            if(index(line,'PRIMCOORD').gt.0) then
                  read(21,*,end=15,err=1000) natoms
                  if (.not. allocated(atoms))then
                        allocate(atoms(natoms))
                        atoms%name=' '
                  end if
                  do iatom=1,natoms
                    ! begin old colde
                    !read(21,*) ispecies,x(1:3)
                    ! end old code
                    ! begin new code
                    read(21,'(A256)') line2
                    read(line2,*) elechar
                    ispecies=-1
                    do iele=1,size(elements)
                      if(trim(elements(iele)).eq.
     &                   trim(elechar)) then
                         ispecies=iele
                      end if
                    end do
                    if (ispecies.lt.0) then
                          read(line2,*) ispecies,x(1:3)
                    else
                          read(line2,*) elechar,x(1:3)
                    end if
                    ! end new code  
                    ! store absolute coordinates
                    atoms(iatom)%abswhere=x
                    ! get fractional coordinates  
                    call abs2frac(x,vecs,atoms(iatom)%where)
                    atoms(iatom)%name(1:2)=elements(ispecies)
                    atoms(iatom)%core="core"
                    atoms(iatom)%charge=0.0D0  ! fake charge
                    atoms(iatom)%occ=1.0D0  ! default occupation
                    !newspecies=.true.
                    !do iatom2=1,iatom-1
                    !  if (atoms(iatom2)%name(1:2).eq.
     &              !      atoms(iatom)%name(1:2))newspecies=.false.
                    !end do
                    !if(newspecies) nspecies=nspecies+1
                  end do
                  goto 15
            end if
      end if
      goto 12
15    continue
      ! get number of species
!      nspecies=0
!      do iatom=1,natoms
!         newspecies=.true.
!         do iatom2=1,iatom-1
!           if (atoms(iatom2)%name(1:4).eq.
!     &         atoms(iatom)%name(1:4))newspecies=.false.
!         end do
!         if(newspecies) nspecies=nspecies+1
!      end do
      call getspecies(atoms,species)
      nspecies=size(species)
      call get_masses_of_species(species)
      call get_masses_of_atoms(atoms)
      if (isopen12) then
!        write(12,'(8x,"Number of elements and atoms found: ",I,I)'),
!     &        nspecies,natoms  
        WRITE(FMT2,*) nspecies
        WRITE(FMT3,*) natoms
        WRITE(FMT1,*) '(8x,"Number of elements and atoms found: ",I',
     &   len(trim(FMT2)),',I',len(trim(FMT3)),')'
        write(12,FMT1)  nspecies,natoms 
      else
!        print '(8x,"Number of elements and atoms found: ",I,I)',
!     &        nspecies,natoms  
        WRITE(FMT2,*) nspecies
        WRITE(FMT3,*) natoms
        WRITE(FMT1,*) '(8x,"Number of elements and atoms found: ",I',
     &   len(trim(FMT2)),',I',len(trim(FMT3)),')'
        print FMT1,  nspecies,natoms 
      end if
      close(21)
      !get name and number of atoms of each species
      if (.not. allocated(species)) then
            allocate(species(nspecies))
            species%name=" "
      end if
!      ispecies=0
!      do iatom=1,natoms
!        newspecies=.true.
!        do iatom2=1,iatom-1
!          if (atoms(iatom2)%name(1:4).eq.
!     &        atoms(iatom)%name(1:4))newspecies=.false.
!        end do
!        if(newspecies) then
!              ispecies=ispecies+1
!              species(ispecies)%name=atoms(iatom)%name
!              species(ispecies)%howmany=1
!        else
!              species(ispecies)%howmany=species(ispecies)%howmany+1
!        end if
!      end do
      
      ispecies=1
      do iatom=1,natoms
        newspecies=.true.
        do i=1,iatom-1
          if (atoms(i)%name(1:2).eq.atoms(iatom)%name(1:2))
     &        newspecies=.false.
        end do
        if(newspecies) then
            species(ispecies)%name=atoms(iatom)%name
            ispecies=ispecies+1
        end if
      end do  

      ! get numbers of atoms of each species
      do ispecies=1,nspecies
            species(ispecies)%howmany=0
            do i=1,natoms
               if(atoms(i)%name(1:2).eq.species(ispecies)%name(1:2)) 
     &            species(ispecies)%howmany=species(ispecies)%howmany+1
            end do
      end do

c      species(1)%name(1:2)=atoms(1)%name(1:2)
c      species(1)%howmany=1
c      do iatom=2,natoms
c        if(atoms(iatom)%name(1:2).eq.species(1)%name(1:2))
c     &     species(1)%howmany=species(1)%howmany+1
c      end do
c
c      do ispecies=2,nspecies
c        newspecies=.true.
c        do iatom=1,natom
c          do i=1,ispecies-1 
c            if (atoms(iatom)%name(1:2).eq.species(i)%name(1:2))
c     &          newspecies=.false.
c          end do
c        end do  
c     &     species(ispecies-1)%name(1:2)) iatom=iatom+1
c        end do
c        
c        
c        iatom=2,natoms
c          if(atoms(iatom)%name(1:2).eq.species(1)%name(1:2))
c     &       species(1)%howmany=species(1)%howmany+1
c        end do
c      end do
c
c
c        iatom=1
c        do while (
c      newspecies=.true
c        do iatom=1,natoms
c          if
c        do i=1,ispecies-1
c          if (species(i)%name(1:2).eq.species(ispecies)%name(1:2)
c     &       newspecies=.false.     

c23    continue
c      if(isopen12) then
c            write(12,ferrmssg)'coordinates section ended' 
c      else
c            print ferrmssg,'coordinates section ended unexpectedly'
c      end if
c      nerr=nerr+1
25    continue
      if(isopen12) then
            write(12,'(8x,"Species found in the file (name, howmany, mas
     &s):")') 
            do i=1,nspecies
!              write(12,'(8x,A2,1x,I,1x,F10.6)') species(i)%name,
!     &        species(i)%howmany,species(i)%mass         
               WRITE(FMT2,*) species(i)%howmany
               WRITE(FMT1,*) '(8x,A2,1x,I',len(trim(FMT2)),',1x,F10.6)'
               write(12,FMT1) species(i)%name,species(i)%howmany,
     &            species(i)%mass
            end do
            !write(12,'(8x,A)')'Atoms found in the file:' 
            !do i=1,natoms
            !  write(12,*) atoms(i)
            !end do
      else
            print '(8x,"Species found in the file (name, howmany, mass):
     &")' 
            do i=1,nspecies
!              print '(8x,A2,1x,I,1x,F10.6)', species(i)%name,
!     &        species(i)%howmany,species(i)%mass         
               WRITE(FMT2,*) species(i)%howmany
               WRITE(FMT1,*) '(8x,A2,1x,I',len(trim(FMT2)),',1x,F10.6)'
               print FMT1, species(i)%name,species(i)%howmany,
     &            species(i)%mass
            end do
            !print '(8x,A)','Atoms found in the file:' 
            !do i=1,natoms
            !  print *, atoms(i)
            !end do
      end if
      close(22)
 
      if(isopen12) then
            write(12,fsubendext) 'read_xsf'
      else
            print fsubendext, 'read_xsf'
      end if
      return
      
999   continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Cannot open xsf file" 
      else
            print ferrmssg,"Cannot open xsf file" 
      end if
      close(21)
      stop
1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Something wrong with xsf file" 
      else
            print ferrmssg,"Something wrong with xsf file" 
      end if
      close(21)
      stop
      end subroutine read_xsf

c---------------------------------------------------------------------

      subroutine read_xsf_forces(fincoords,atoms,natoms,species,          &
     &           nspecies,vecs)

      use defs
      use misc
      implicit none
      
      ! variables
      character(len=*), intent(in) :: fincoords
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms,nspecies
      character line*100,nameat*2
      double precision x(1:3),vecs(1:3,1:3),f(1:3)
      
      ! internal variables
      logical isopen12,newspecies
      integer natom,i,j,iatom,iatom2,ispecies
      integer iele
      character elechar*40,line2*256
      !character FMT1*1024,FMT2*1024,FMT3*1024 
      
      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim(adjustl('read_xsf_forces'))
      else
            print fsubstart, trim(adjustl('read_xsf_forces'))
      end if
      
      ! read xsf file
      open(21,file=fincoords,status='old',err=999)
      ! get lattice vectors
10    read(21,'(A100)',end=1000,err=1000) line
      line=adjustl(line)
      if (line(1:1).ne.'!'.and.line(1:1).ne.'#') then
            if(index(line,'PRIMVEC').gt.0) then
                  do i=1,3
                    read(21,*) vecs(i,1:3)
                  end do
                  goto 12
            end if
      end if
      goto 10

      ! get number of atoms
      natoms=0
      !nspecies=0
12    read(21,'(A100)',end=15,err=1000) line
      line=adjustl(line)
      if (line(1:1).ne.'!'.and.line(1:1).ne.'#') then
            if(index(line,'PRIMCOORD').gt.0) then
                  read(21,*,end=15,err=1000) natoms
                  if (.not. allocated(atoms))then
                        allocate(atoms(natoms))
                        atoms%name=' '
                  end if
                  do iatom=1,natoms
                    ! begin old colde
                    !read(21,*) ispecies,x(1:3)
                    ! end old code
                    ! begin new code
                    read(21,'(A256)') line2
                    read(line2,*) elechar
                    ispecies=-1
                    do iele=1,size(elements)
                      if(trim(elements(iele)).eq.
     &                   trim(elechar)) then
                         ispecies=iele
                      end if
                    end do
                    if (ispecies.lt.0) then
                          read(line2,*) ispecies,x(1:3),f(1:3)
                    else
                          read(line2,*) elechar,x(1:3),f(1:3)
                    end if
                    ! end new code  
                    ! store absolute coordinates
                    atoms(iatom)%abswhere=x
                    atoms(iatom)%force=f
                    ! get fractional coordinates  
                    call abs2frac(x,vecs,atoms(iatom)%where)
                    atoms(iatom)%name(1:2)=elements(ispecies)
                    atoms(iatom)%core="core"
                    atoms(iatom)%charge=0.0D0  ! fake charge
                    atoms(iatom)%occ=1.0D0  ! default occupation
                    !newspecies=.true.
                    !do iatom2=1,iatom-1
                    !  if (atoms(iatom2)%name(1:2).eq.
     &              !      atoms(iatom)%name(1:2))newspecies=.false.
                    !end do
                    !if(newspecies) nspecies=nspecies+1
                  end do
                  goto 15
            end if
      end if
      goto 12
15    continue
      ! get number of species
!      nspecies=0
!      do iatom=1,natoms
!         newspecies=.true.
!         do iatom2=1,iatom-1
!           if (atoms(iatom2)%name(1:4).eq.
!     &         atoms(iatom)%name(1:4))newspecies=.false.
!         end do
!         if(newspecies) nspecies=nspecies+1
!      end do
      call getspecies(atoms,species)
      nspecies=size(species)
      call get_masses_of_species(species)
      call get_masses_of_atoms(atoms)
      if (isopen12) then
!        write(12,'(8x,"Number of elements and atoms found: ",I,I)'),
!     &        nspecies,natoms  
        WRITE(FMT2,*) nspecies
        WRITE(FMT3,*) natoms
        WRITE(FMT1,*) '(8x,"Number of elements and atoms found: ",I',
     &   len(trim(FMT2)),',I',len(trim(FMT3)),')'
        write(12,FMT1)  nspecies,natoms 
      else
!        print '(8x,"Number of elements and atoms found: ",I,I)',
!     &        nspecies,natoms  
        WRITE(FMT2,*) nspecies
        WRITE(FMT3,*) natoms
        WRITE(FMT1,*) '(8x,"Number of elements and atoms found: ",I',
     &   len(trim(FMT2)),',I',len(trim(FMT3)),')'
        print FMT1,  nspecies,natoms 
      end if
      close(21)
      !get name and number of atoms of each species
      if (.not. allocated(species)) then
            allocate(species(nspecies))
            species%name=" "
      end if
!      ispecies=0
!      do iatom=1,natoms
!        newspecies=.true.
!        do iatom2=1,iatom-1
!          if (atoms(iatom2)%name(1:4).eq.
!     &        atoms(iatom)%name(1:4))newspecies=.false.
!        end do
!        if(newspecies) then
!              ispecies=ispecies+1
!              species(ispecies)%name=atoms(iatom)%name
!              species(ispecies)%howmany=1
!        else
!              species(ispecies)%howmany=species(ispecies)%howmany+1
!        end if
!      end do
      
      ispecies=1
      do iatom=1,natoms
        newspecies=.true.
        do i=1,iatom-1
          if (atoms(i)%name(1:2).eq.atoms(iatom)%name(1:2))
     &        newspecies=.false.
        end do
        if(newspecies) then
            species(ispecies)%name=atoms(iatom)%name
            ispecies=ispecies+1
        end if
      end do  

      ! get numbers of atoms of each species
      do ispecies=1,nspecies
            species(ispecies)%howmany=0
            do i=1,natoms
               if(atoms(i)%name(1:2).eq.species(ispecies)%name(1:2)) 
     &            species(ispecies)%howmany=species(ispecies)%howmany+1
            end do
      end do

c      species(1)%name(1:2)=atoms(1)%name(1:2)
c      species(1)%howmany=1
c      do iatom=2,natoms
c        if(atoms(iatom)%name(1:2).eq.species(1)%name(1:2))
c     &     species(1)%howmany=species(1)%howmany+1
c      end do
c
c      do ispecies=2,nspecies
c        newspecies=.true.
c        do iatom=1,natom
c          do i=1,ispecies-1 
c            if (atoms(iatom)%name(1:2).eq.species(i)%name(1:2))
c     &          newspecies=.false.
c          end do
c        end do  
c     &     species(ispecies-1)%name(1:2)) iatom=iatom+1
c        end do
c        
c        
c        iatom=2,natoms
c          if(atoms(iatom)%name(1:2).eq.species(1)%name(1:2))
c     &       species(1)%howmany=species(1)%howmany+1
c        end do
c      end do
c
c
c        iatom=1
c        do while (
c      newspecies=.true
c        do iatom=1,natoms
c          if
c        do i=1,ispecies-1
c          if (species(i)%name(1:2).eq.species(ispecies)%name(1:2)
c     &       newspecies=.false.     

c23    continue
c      if(isopen12) then
c            write(12,ferrmssg)'coordinates section ended' 
c      else
c            print ferrmssg,'coordinates section ended unexpectedly'
c      end if
c      nerr=nerr+1
25    continue
      if(isopen12) then
            write(12,'(8x,"Species found in the file (name, howmany, mas
     &s):")') 
            do i=1,nspecies
!              write(12,'(8x,A2,1x,I,1x,F10.6)') species(i)%name,
!     &        species(i)%howmany,species(i)%mass         
               WRITE(FMT2,*) species(i)%howmany
               WRITE(FMT1,*) '(8x,A2,1x,I',len(trim(FMT2)),',1x,F10.6)'
               write(12,FMT1) species(i)%name,species(i)%howmany,
     &            species(i)%mass
            end do
            !write(12,'(8x,A)')'Atoms found in the file:' 
            !do i=1,natoms
            !  write(12,*) atoms(i)
            !end do
      else
            print '(8x,"Species found in the file (name, howmany, mass):
     &")' 
            do i=1,nspecies
!              print '(8x,A2,1x,I,1x,F10.6)', species(i)%name,
!     &        species(i)%howmany,species(i)%mass         
               WRITE(FMT2,*) species(i)%howmany
               WRITE(FMT1,*) '(8x,A2,1x,I',len(trim(FMT2)),',1x,F10.6)'
               print FMT1, species(i)%name,species(i)%howmany,
     &            species(i)%mass
            end do
            !print '(8x,A)','Atoms found in the file:' 
            !do i=1,natoms
            !  print *, atoms(i)
            !end do
      end if
      close(22)
 
      if(isopen12) then
            write(12,fsubendext) 'read_xsf'
      else
            print fsubendext, 'read_xsf'
      end if
      return
      
999   continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Cannot open xsf file" 
      else
            print ferrmssg,"Cannot open xsf file" 
      end if
      close(21)
      stop
1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Something wrong with xsf file" 
      else
            print ferrmssg,"Something wrong with xsf file" 
      end if
      close(21)
      stop
      end subroutine read_xsf_forces

c---------------------------------------------------------------------

      subroutine read_cfg(infile,atoms,natoms,species,nspecies,
     &                    vecs)
      ! reads atomeye-readable cfg file 
      use defs
      use misc, only : frac2abs
      implicit none

      ! global variables
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      double precision vecs(1:3,1:3)
      integer natoms, nspecies
      
      ! local variables
      logical isopen12
      double precision rcut,propav,propprodlocav,binwidth
      double precision, allocatable ::
     & props(:,:),checkprops(:),propsloc(:)
      character(len=*) infile
      character ele*2,line*100,elestring*100,corrf*20
      integer iele,i,j,iprop,nprops,numele,ieleprop
      integer ncheckatoms,nproploc
      ! for correlation function:
      integer :: nr
      double precision rmax,corrint,maxint
      double precision, allocatable :: gcorr(:),r(:),rdistf(:)
      parameter(rmax=20.0D0)

      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim(adjustl('read_cfg'))
      else
            print fsubstart, trim(adjustl('read_cfg'))
      end if

      !read(10,*) infile
      !read(10,*) ele
      !read(10,*) iprop
      !read(10,*) corrf
      !read(10,*) rcut
      !binwidth=0.2D0
      !read(10,'(A100)',end=4,err=4) line
      !if (index(line,'binwidth').ge.1) then
      !      read(line(index(line,'binwidth')+9:100),*,end=4,err=4) 
!     &       binwidth
      !end if
!4     continue      
!      nr=int(rmax/binwidth)
!      allocate(gcorr(1:nr),r(1:nr),rdistf(nr))
      if(isopen12) then
        write(12,'(8x,"File to read is ",A,".")') trim(adjustl(infile))
      else
        print '(8x,"File to read is ",A,".")', trim(adjustl(infile))
      end if
!      write(12,'(8x,"Will check ",A2," for correlations")') 
!     & trim(adjustl(ele))
!      select case(corrf)
!      case('theta') 
!            write(12,'(8x,"using ",A," as corr. function")')
!     &          trim(adjustl(corrf)) 
!            write(12,'(8x,"with cutoff",F10.3," Angs.")') rcut
!            write(12,'(8x,"with binwidth",F10.3," Angs.")') binwidth
!      case default
!            write(12,'(8x,"Error: Unknown corr. funct.")')
!            print*,"Error: Unknown corr. funct."
!            nerr=nerr+1
!            return
!      end select      
      open(65,file=infile,status='old')
      ! read number of atoms from cfg file
      read(65,'(A100)') line
      line=adjustl(line)
      do while (index(line,'Number of particles').le.0
     &          .or.line(1:1).eq.'#')
            read(65,'(A100)') line
            line=adjustl(line)
      end do
      read(line(index(line,'=')+1:100),*) natoms
      if(isopen12) then
!        write(12,'(8x,"File ",A," contains ",I," atoms.")')
!     & trim(adjustl(infile)),natoms
        WRITE(FMT2,*) natoms
        WRITE(FMT1,*) '(8x,"File ",A',len(trim(infile)),', "contains ",
     &    I',len(trim(FMT2)),'," atoms.")'
        write(12,FMT1) trim(adjustl(infile)),natoms
      else
!        print '(8x,"File ",A," contains ",I," atoms.")',
!     & trim(adjustl(infile)),natoms
        WRITE(FMT2,*) natoms
        WRITE(FMT1,*) '(8x,"File ",A',len(trim(infile)),', "contains ",
     &    I',len(trim(FMT2)),'," atoms.")'
        print FMT1, trim(adjustl(infile)),natoms
      end if
      !write(12,'(8x,"Will check aux. prop. #",I3,".")'),iprop
      ! read cell vectors from  cfg file
      read(65,'(A100)') line
      line=adjustl(line)
      do while (index(line,'H0(1,1)').le.0.or.line(1:1).eq.'#')
            read(65,'(A100)') line
            line=adjustl(line)
      end do
      read(line(index(line,'=')+1:100),*) vecs(1,1)
      do j=2,3
            read(65,'(A100)') line
            line=adjustl(line)
            read(line(index(line,'=')+1:100),*) vecs(1,j)
      end do
      if(isopen12) then
        write(12,'(8x,"Cell vectors in Angs:")')
      else
        print '(8x,"Cell vectors in Angs:")'
      end if
      if(isopen12) then
        write(12,'(8x,3(F12.6))') vecs(1,1:3)
      else
        print '(8x,3(F12.6))', vecs(1,1:3)
      end if
      do i=2,3
            do j=1,3
                  read(65,'(A100)') line
                  line=adjustl(line)
                  read(line(index(line,'=')+1:100),*) vecs(i,j)
            end do
            if(isopen12) then
              write(12,'(8x,3(F12.6))') vecs(i,1:3)
            else
              print '(8x,3(F12.6))', vecs(i,1:3)
            end if
      end do
      ! find number of auxiliary properties
      do while(index(line,'entry_count').le.0.or.line(1:1).eq.'#')
            read(65,'(A100)') line
            line=adjustl(line)
      end do
      read(line(index(line,'=')+1:100),*) nprops
      nprops=nprops-3
      if(isopen12) then
        write(12,'(8x,"There are ",I3," auxil. properties.")')nprops
      else
        print '(8x,"There are ",I3," auxiliary properties.")',nprops
      end if
      ! find elements
      nspecies=0
      elestring(1:100)=' '
10    read(65,'(A100)',end=15) line
      line=adjustl(line)
      do i=1,100
            if (line(1:2).eq.elements(i)) then
                  elestring(2*nspecies+1:2*nspecies+2)=elements(i)
                  nspecies=nspecies+1
            end if
      end do
      goto 10
15    continue
      if (.not. allocated(species))then
            allocate(species(nspecies))
            species%name=" "
      end if
      if(isopen12) then
        write(12,'(8x,"Found ",I3," elements:")') nspecies
      else
        print '(8x,"Found ",I3," elements:")', nspecies
      end if
      do i=1,nspecies
            species(i)%name(1:2)=elestring(2*i-1:2*i)
      end do
      ! get numbers of atoms of each element
      rewind(65)
      do while (line(1:2).ne.species(1)%name(1:2))
            read(65,'(A100)')line
            line=adjustl(line)
      end do
      do i=2,nspecies
            numele=0
            do while (line(1:2).ne.species(i)%name(1:2))
                  read(65,'(A100)')line
                  line=adjustl(line)
                  numele=numele+1
            end do
            species(i-1)%howmany=numele-2
      end do
      numele=0
20    read(65,'(A100)',end=25) line    
      line=adjustl(line)
      numele=numele+1
      goto 20
25    continue
      species(nspecies)%howmany=numele
      do i=1,nspecies
            species(i)%name(1:2)=elestring(2*i-1:2*i)
            write(12,'(8x,A2,I10)')species(i)%name,species(i)%howmany
      end do
      ! read atom coordinates
      if (.not. allocated(atoms))then
            allocate(atoms(natoms))!,props(1:natoms,1:nprops))
            do i=1,natoms
              allocate(atoms(i)%properties(1:nprops))
            end do
            atoms%name=" "
      end if
      rewind(65)
      do while (line(1:2).ne.species(1)%name(1:2))
            read(65,'(A100)') line
            line=adjustl(line)
      end do
      i=1
      if(isopen12) then
        write(12,'(8x,"Coord. of first and last atom",
     &      " of each species:")')
      else
        print '(8x,"Coord. of first and last atom",
     &      " of each species:")'
      end if
      do iele=1,nspecies-1
            do j=1,species(iele)%howmany
                  read(65,*)atoms(i)%where(1:3),
     &               atoms(i)%properties(1:nprops)
                  atoms(i)%name=species(iele)%name
                  call frac2abs(atoms(i)%where,vecs,                      &
     &                          atoms(i)%abswhere)
                  if (j.eq.1.or.j.eq.species(iele)%howmany) then
                        if(isopen12) then
                          write(12,'(8x,A2,3(F15.6),F15.6)') 
     &                    atoms(i)%name(1:2),atoms(i)%where,
     &                    atoms(i)%properties(1:nprops)  
                        else
                          print '(8x,A2,3(F15.6),F15.6)', 
     &                    atoms(i)%name(1:2),atoms(i)%where,
     &                    atoms(i)%properties(1:nprops)  
                        end if
                  end if
                  i=i+1
            end do
            read(65,'(A100)') line
            read(65,'(A100)') line
      end do
      do j=1,species(nspecies)%howmany
            read(65,*)atoms(i)%where(1:3),atoms(i)%properties(1:nprops)
            atoms(i)%name=species(nspecies)%name
            call frac2abs(atoms(i)%where,vecs,atoms(i)%abswhere)
            if (j.eq.1.or.j.eq.species(nspecies)%howmany) then
                  if(isopen12) then
                    write(12,'(8x,A2,3(F15.6),F15.6)') 
     &              atoms(i)%name(1:2),atoms(i)%where,
     &            atoms(i)%properties(1:nprops)
                  else
                    print '(8x,A2,3(F15.6),F15.6)', 
     &              atoms(i)%name(1:2),atoms(i)%where,
     &              atoms(i)%properties(1:nprops)
                  end if
            end if
            i=i+1
      end do
      close(65)
      atoms%occ=1.0D0 !default occupation
      
      if(isopen12) then
            write(12,fsubendext),'read_cfg'
      else
            print fsubendext, 'read_cfg'
      end if
      return

      end subroutine

c---------------------------------------------------------------------
      
      subroutine read_lxyz(fincoords,atoms,natoms,species,nspecies,
     &           vecs)

      use defs
      use misc
      implicit none
      
      ! variables
      character(len=*), intent(in) :: fincoords
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms,nspecies
      character line*100,nameat*2
      double precision x(1:3),vecs(1:3,1:3),xy,xz,yz
      
      ! internal variables
      double precision xlo,xhi,ylo,yhi,zlo,zhi
      logical isopen12,newspecies
      integer natom,i,j,iatom,iatom2,ispecies,idum
      !character FMT1*1024,FMT2*1024,FMT3*1024
      
      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim(adjustl('read_lxyz'))
      else
            print fsubstart, trim(adjustl('read_lxyz'))
      end if
      
      ! read lammps xyz file
      open(21,file=fincoords,status='old',err=1000)
      ! get lattice vectors
10    read(21,'(A100)',end=1000,err=1000) line
      line=adjustl(line)
      if (line(1:1).ne.'!'.and.line(1:1).ne.'#') then
            if(index(line,'BOX BOUNDS').gt.0) then
                  xy=0.0D0
                  xz=0.0D0
                  yz=0.0D0
                  if(index(line,'xy').ge.1) then
                        read(21,*) xlo,xhi,xy
                  else
                        read(21,*) xlo,xhi
                  end if
                  if(index(line,'xz').ge.1) then
                        read(21,*) ylo,yhi,xz
                  else
                        read(21,*) ylo,yhi
                  end if
                  if(index(line,'yz').ge.1) then
                        read(21,*) zlo,zhi,yz
                  else
                        read(21,*) zlo,zhi
                  end if
                  vecs=0.0D0
                  vecs(1,1)=xhi-xlo-xy-xz
                  vecs(2,1)=xy
                  vecs(2,2)=yhi-ylo-yz
                  vecs(3,1)=xz
                  vecs(3,2)=yz
                  vecs(3,3)=zhi-zlo
                  rewind(21)
                  goto 12
            end if
      end if
      goto 10

      ! get number of atoms
      natoms=0
12    read(21,'(A100)',end=1000,err=1000) line
      line=adjustl(line)
      if (line(1:1).ne.'!'.and.line(1:1).ne.'#') then
            if(index(line,'ITEM: NUMBER OF ATOMS').gt.0) then
                  read(21,*,end=1000,err=1000) natoms
                  if (.not. allocated(atoms))then
                        allocate(atoms(natoms))
                        atoms%name=" "
                  end if
                  goto 13
            end if
      end if
      goto 12
      
      ! get atom coordinates
13    read(21,'(A100)',end=15,err=1000) line
      line=adjustl(line)
      if (line(1:1).ne.'!'.and.line(1:1).ne.'#') then
            if(index(line,'ITEM: ATOMS type x').gt.0) then
                  atoms%name=' '
                  do iatom=1,natoms
                    read(21,*) ispecies,x(1:3)
                    atoms(iatom)%abswhere=x
                    call abs2frac(x,vecs,atoms(iatom)%where)
                    atoms(iatom)%name(1:2)=elements(ispecies)
                    atoms(iatom)%core="core"
                    atoms(iatom)%charge=0.0D0  ! fake charge
                  end do
                  goto 15
            end if
            if(index(line,'ITEM: ATOMS id type x').gt.0) then
                  atoms%name=' '
                  do iatom=1,natoms
                    read(21,*) idum,ispecies,x(1:3)
                    atoms(iatom)%abswhere=x
                    call abs2frac(x,vecs,atoms(iatom)%where)
                    atoms(iatom)%name(1:2)=elements(ispecies)
                    atoms(iatom)%core="core"
                    atoms(iatom)%charge=0.0D0  ! fake charge
                  end do
                  goto 15
            end if
      end if
      goto 13

15    continue
      ! get number of species
      nspecies=0
      do iatom=1,natoms
         newspecies=.true.
         do iatom2=1,iatom-1
           if (atoms(iatom2)%name(1:4).eq.
     &         atoms(iatom)%name(1:4))newspecies=.false.
         end do
         if(newspecies) nspecies=nspecies+1
      end do
      if (isopen12) then
!        write(12,'(8x,"Number of elements and atoms found: ",I,I)'),
!     &        nspecies,natoms  
        WRITE(FMT2,*) nspecies
        WRITE(FMT3,*) natoms
        WRITE(FMT1,*) '(8x,"Number of elements and atoms found: ",I',
     &   len(trim(FMT2)),',I',len(trim(FMT3)),')'
        write(12,FMT1) nspecies,natoms 
      else
!        print '(8x,"Number of elements and atoms found: ",I,I)',
!     &        nspecies,natoms  
        WRITE(FMT2,*) nspecies
        WRITE(FMT3,*) natoms
        WRITE(FMT1,*) '(8x,"Number of elements and atoms found: ",I',
     &   len(trim(FMT2)),',I',len(trim(FMT3)),')'
        print FMT1,  nspecies,natoms 
      end if
      close(21)
      !get name and number of atoms of each species
      if (.not. allocated(species))then
            allocate(species(nspecies))
            species%name=" "
      end if
      
      ispecies=1
      do iatom=1,natoms
        newspecies=.true.
        do i=1,iatom-1
          if (atoms(i)%name(1:2).eq.atoms(iatom)%name(1:2))
     &        newspecies=.false.
        end do
        if(newspecies) then
            species(ispecies)%name=atoms(iatom)%name
            ispecies=ispecies+1
        end if
      end do  

      ! get numbers of atoms of each species
      do ispecies=1,nspecies
            species(ispecies)%howmany=0
            do i=1,natoms
               if(atoms(i)%name(1:2).eq.species(ispecies)%name(1:2)) 
     &            species(ispecies)%howmany=species(ispecies)%howmany+1
            end do
      end do

25    continue
      if(isopen12) then
            write(12,'(8x,"Species found in the file (name, howmany, mas
     &s):")') 
            do i=1,nspecies
!              write(12,'(8x,A2,1x,I,1x,F10.6)') species(i)%name,
!     &        species(i)%howmany,species(i)%mass         
              WRITE(FMT2,*) species(i)%howmany
              WRITE(FMT1,*) '(8x,A2,1x,I',len(trim(FMT2)),',1x,F10.6)'
              write (12,FMT1) species(i)%name,species(i)%howmany,
     &            species(i)%mass
            end do
      else
            print '(8x,"Species found in the file (name, howmany, mass):
     &")' 
            do i=1,nspecies
!              print '(8x,A2,1x,I,1x,F10.6)', species(i)%name,
!     &        species(i)%howmany,species(i)%mass         
              WRITE(FMT2,*) species(i)%howmany
              WRITE(FMT1,*) '(8x,A2,1x,I',len(trim(FMT2)),',1x,F10.6)'
              print FMT1, species(i)%name,species(i)%howmany,
     &            species(i)%mass
            end do
      end if
      close(22)
      atoms%occ=1.0D0  ! default occupation

      if(isopen12) then
            write(12,fsubendext) 'read_lxyz'
      else
            print fsubendext, 'read_lxyz'
      end if
      return
      
1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Something wrong with lxyz file." 
      else
            print ferrmssg,"Something wrong with lxyz file." 
      end if
      close(21)
      stop
      end subroutine

c---------------------------------------------------------------------
      
      subroutine read_gin(fincoords,atoms,natoms,species,nspecies,
     &           vecs)

      use defs
      use misc
      implicit none
      
      ! variables
      character(len=*), intent(in) :: fincoords
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms,nspecies
      character line*100,nameat*2
      double precision x(1:3),vecs(1:3,1:3)
      
      ! internal variables
      double precision alat(1:3),angles(1:3),counter,denominator
      logical isopen12,newspecies,cart
      integer natom,i,j,iatom,iatom2,ispecies,irun
      character coordchar1*20,coordchar2*20,coordchar3*20
      !character FMT1*1024,FMT2*1024,FMT3*1024
      
      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim(adjustl('read_gin'))
      else
            print fsubstart, trim(adjustl('read_gin'))
      end if
      
      ! read gulp input/restart file
      open(21,file=fincoords,status='old',err=1000)
      ! get lattice vectors
10    read(21,'(A100)',end=1000,err=1000) line
      line=adjustl(line)
      if (line(1:1).ne.'#') then
            if(line(1:4).eq.'cell') then
                  read(21,'(A100)',end=1000) line
                  line=adjustl(line)
                  do while (line(1:1).eq."#")
                    read(21,'(A100)',end=1000) line
                    line=adjustl(line)
                  end do
                  read(line,*) alat(1:3),angles(1:3)
                  vecs=0.0D0
!                  vecs(1,1)=alat(1)
!                  vecs(2,1)=alat(2)*cos(angles(3)*Pi/180.0D0)
!                  vecs(2,2)=alat(2)*sin(angles(3)*Pi/180.0D0)
!                  vecs(3,1)=alat(3)*cos(angles(2)*Pi/180.0D0)
!                  vecs(3,2)=alat(3)*(cos(angles(1)*Pi/180.0D0)
!     &            -cos(angles(3)*Pi/180.0D0)*cos(angles(2)*Pi/180.0D0))
!     &            /sin(angles(3)*Pi/180.0D0)
!                  vecs(3,3)=sqrt(alat(3)**2-vecs(3,1)**2-vecs(3,2)**2)
                  call cellpars2vecs((/alat,angles/),vecs)
                  rewind(21)
                  goto 11
            end if
            if(line(1:3).eq.'vec') then
                  read(21,'(A100)',end=1000) line
                  line=adjustl(line)
                  do while (line(1:1).eq."#")
                    read(21,'(A100)',end=1000) line
                    line=adjustl(line)
                  end do
                  read(line,*) vecs(1,1:3)
                  read(21,*) vecs(2,1:3)
                  read(21,*) vecs(3,1:3)
                  rewind(21)
                  goto 11
            end if
      end if
      goto 10

      ! get number of atoms
11    read(21,'(A100)',end=1000,err=1000) line
      line=adjustl(line)
      if(line(1:1).eq.'#') goto 11
      if (line(1:4).eq.'frac'.or.line(1:4).eq.'cart') then
         natoms=0
         read(21,'(A100)',end=1000,err=1000) line
         line=adjustl(line)
         do while (line(1:1).eq.'#'.or.any(elements.eq.line(1:2)))
           if(any(elements.eq.line(1:2)))
     &        natoms=natoms+1
           read(21,'(A100)',end=1000,err=1000) line
           line=adjustl(line)
         end do
         if (.not. allocated(atoms))then
               allocate(atoms(natoms))
               atoms%name=" "
         end if
         goto 12
      end if
      goto 11

      ! get atom coordinates
12    rewind(21)
      iatom=0
      cart=.false.
      atoms%name=' '
      atoms%occ=1.0D0  ! default occupation
      atoms%charge=0.0D0  ! default charge
13    read(21,'(A100)',end=15,err=1000) line
      !line=trim(line)
      line=trim(adjustl(line))
      if(line(1:1).eq.'#') goto 13
      if (line(1:4).eq.'frac'.or.line(1:4).eq.'cart') then
         if (line(1:4).eq.'cart') cart=.true.   
         read(21,'(A100)',end=1000,err=1000) line
         line=adjustl(line)
         do while (line(1:1).eq.'#'.or.any(elements.eq.line(1:2)))
           if(any(elements.eq.line(1:2))) then
              iatom=iatom+1
!              read(line,*,err=1000) atoms(iatom)%name,atoms(iatom)%core,
!     &        atoms(iatom)%where,atoms(iatom)%charge
              atoms(iatom)%BEC(1:3,1:3)=0.0D0  ! default BEC
              do i=1,len(line)
                if(line(i:i).eq.'/') line(i:i)='d'
              end do
              read(line,*,err=1000) atoms(iatom)%name,
     &            atoms(iatom)%core, coordchar1,coordchar2,               &
     &            coordchar3,atoms(iatom)%charge,atoms(iatom)%occ
              atoms(iatom)%BEC(1,1)=atoms(iatom)%charge
              atoms(iatom)%BEC(2,2)=atoms(iatom)%charge
              atoms(iatom)%BEC(3,3)=atoms(iatom)%charge
              ! if a '/' was found, calculate the number
              if (index(coordchar1,'d').gt.0) then
                    read(coordchar1(1:1),*)counter
                    read(coordchar1(3:3),*) denominator
                    atoms(iatom)%where(1)=counter/denominator
              else
                    read(coordchar1,*) atoms(iatom)%where(1)
              end if
              if (index(coordchar2,'d').gt.0) then
                    read(coordchar2(1:1),*)counter
                    read(coordchar2(3:3),*) denominator
                    atoms(iatom)%where(2)=counter/denominator
              else
                    read(coordchar2,*) atoms(iatom)%where(2)
              end if
              if (index(coordchar3,'d').gt.0) then
                    read(coordchar3(1:1),*)counter
                    read(coordchar3(3:3),*) denominator
                    atoms(iatom)%where(3)=counter/denominator
              else
                    read(coordchar3,*) atoms(iatom)%where(3)
              end if
              ! get absolute coordinates
              call frac2abs(atoms(iatom)%where,vecs,
     &          atoms(iatom)%abswhere)
              if(cart) then
                    x=atoms(iatom)%where
                    atoms(iatom)%abswhere=x  ! absolute coordinates
                    call abs2frac(x,vecs,atoms(iatom)%where)  ! fractional coordinates
              end if
           end if
           read(21,'(A100)',end=1000,err=1000) line
           line=adjustl(line)
         end do
         goto 15
      end if
      goto 13

15    continue
      ! get number of species
      nspecies=0
      do iatom=1,natoms
         newspecies=.true.
         do iatom2=1,iatom-1
           if (atoms(iatom2)%name(1:4).eq.
     &         atoms(iatom)%name(1:4))newspecies=.false.
         end do
         if(newspecies) nspecies=nspecies+1
      end do
      if (isopen12) then
!        write(12,'(8x,"Number of elements and atoms found: ",I,I)'),
!     &        nspecies,natoms  
        WRITE(FMT2,*) nspecies
        WRITE(FMT3,*) natoms
        WRITE(FMT1,*) '(8x,"Number of elements and atoms found: ",I',
     &   len(trim(FMT2)),',I',len(trim(FMT3)),')'
        write(12,FMT1)  nspecies,natoms 
      else
!        print '(8x,"Number of elements and atoms found: ",I,I)',
!     &        nspecies,natoms  
        WRITE(FMT2,*) nspecies
        WRITE(FMT3,*) natoms
        WRITE(FMT1,*) '(8x,"Number of elements and atoms found: ",I',
     &   len(trim(FMT2)),',I',len(trim(FMT3)),')'
        print FMT1,  nspecies,natoms 
      end if
      close(21)
      !get name and number of atoms of each species
      if (.not. allocated(species))then
            allocate(species(nspecies))
            species%name=" "
      end if
            
      ispecies=1
      do iatom=1,natoms
        newspecies=.true.
        do i=1,iatom-1
          if (atoms(i)%name(1:2).eq.atoms(iatom)%name(1:2))
     &        newspecies=.false.
        end do
        if(newspecies) then
            species(ispecies)%name=atoms(iatom)%name
            ispecies=ispecies+1
        end if
      end do  

      ! get numbers of atoms of each species
      do ispecies=1,nspecies
            species(ispecies)%howmany=0
            do i=1,natoms
               if(atoms(i)%name(1:2).eq.species(ispecies)%name(1:2)) 
     &            species(ispecies)%howmany=species(ispecies)%howmany+1
            end do
      end do

25    continue
      call get_masses_of_atoms(atoms)
      call get_masses_of_species(species)
      if(isopen12) then
            write(12,'(8x,A)')'Species found in the file (name, howmany,
     & mass):' 
            do i=1,nspecies
!              write(12,'(8x,A2,1x,I,1x,F10.6)') species(i)%name,
!     &        species(i)%howmany,species(i)%mass         
              WRITE(FMT2,*) species(i)%howmany
              WRITE(FMT1,*)'(8x,A2,1x,I',len(trim(FMT2)),',1x,F10.6)'
              write(12,FMT1) species(i)%name,species(i)%howmany,
     &            species(i)%mass
            end do
      else
            print '(8x,A)','Species found in the file (name, howmany, ma
     &ss):' 
            do i=1,nspecies
!              print '(8x,A2,1x,I,1x,F10.6)', species(i)%name,
!     &        species(i)%howmany,species(i)%mass         
              WRITE(FMT2,*) species(i)%howmany
              WRITE(FMT1,*)'(8x,A2,1x,I',len(trim(FMT2)),',1x,F10.6)'
              write(12,FMT1) species(i)%name,species(i)%howmany,
     &            species(i)%mass
            end do
      end if
      !close(22)
      if(isopen12) then
            write(12,fsubendext) 'read_gin'
      else
            print fsubendext, 'read_gin'
      end if
      return
      
1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Something wrong with gin file." 
      else
            print ferrmssg,"Something wrong with gin file." 
      end if
      close(21)
      stop

      end subroutine read_gin

c---------------------------------------------------------------------

      subroutine read_data(fincoords,atoms,natoms,species,nspecies,
     &           vecs)

      use defs
      use misc
      implicit none
      
      ! variables
      character(len=*), intent(in) :: fincoords
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms,nspecies
      character line*100,nameat*2
      double precision x(1:3),vecs(1:3,1:3)
      
      ! internal variables
      double precision xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz
      logical isopen12,newspecies
      integer natom,i,j,iatom,iatom2,ispecies,irun
      !!character FMT1*1024,FMT2*1024,FMT3*1024
      
      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim(adjustl('read_data'))
      else
            print fsubstart, trim(adjustl('read_data'))
      end if
      
      ! read lammps data file
      open(21,file=fincoords,status='old',err=1000)
      ! get lattice vectors
10    read(21,'(A100)',end=1000,err=1000) line
      line=adjustl(line)
      if (line(1:1).ne.'!'.and.line(1:1).ne.'#') then
      if(index(line,'atom types').gt.0) then
                  read(21,'(A100)',end=1000) line
                  read(21,*) xlo,xhi
                  read(21,*) ylo,yhi
                  read(21,*) zlo,zhi
                  xy=0.0D0
                  xz=0.0D0
                  yz=0.0D0
                  read(21,'(A100)',end=1000) line
                  if (index(line,'xy').gt.0) then
                        read(line(1:100),*,err=1000) xy,xz,yz
                  end if
                  vecs=0.0D0
                  vecs(1,1)=xhi-xlo
                  vecs(2,1)=xy
                  vecs(2,2)=yhi-ylo
                  vecs(3,1)=xz
                  vecs(3,2)=yz
                  vecs(3,3)=zhi-zlo
                  goto 12
            end if
      end if
      goto 10

      ! get number of atoms
12    rewind(21)
      read(21,'(A100)',end=1000,err=1000) line
      read(21,'(A100)',end=1000,err=1000) line
      read(21,*,end=1000,err=1000) natoms
      if (.not. allocated(atoms))then
            allocate(atoms(natoms))
            atoms%name=" "
      end if
      
      ! get atom coordinates
      rewind(21)
13    read(21,'(A100)',end=15,err=1000) line
      line=trim(line)
      if (line.eq.'Atoms'.or.line.eq.'atoms') then
          read(21,'(A100)',end=15,err=1000) line
                  atoms%name=' '
                  do iatom=1,natoms
                    read(21,*) irun,ispecies,atoms(iatom)%charge,x(1:3)
                    atoms(iatom)%abswhere=x
                    call abs2frac(x,vecs,atoms(iatom)%where)
                    atoms(iatom)%name(1:2)=elements(ispecies)
                    atoms(iatom)%core="core"
                    !atoms(iatom)%charge=0.0D0  ! fake charge
                  end do
                  goto 15
            !end if
      end if
      goto 13

15    continue
      ! get number of species
      nspecies=0
      do iatom=1,natoms
         newspecies=.true.
         do iatom2=1,iatom-1
           if (atoms(iatom2)%name(1:4).eq.
     &         atoms(iatom)%name(1:4))newspecies=.false.
         end do
         if(newspecies) nspecies=nspecies+1
      end do
      if (isopen12) then
!        write(12,'(8x,"Number of elements and atoms found: ",I,I)'),
!     &        nspecies,natoms  
        WRITE(FMT2,*) nspecies
        WRITE(FMT3,*) natoms
        WRITE(FMT1,*) '(8x,"Number of elements and atoms found: ",I',
     &   len(trim(FMT2)),',I',len(trim(FMT3)),')'
        write(12,FMT1)  nspecies,natoms 
      else
!        print '(8x,"Number of elements and atoms found: ",I,I)',
!     &        nspecies,natoms  
        WRITE(FMT2,*) nspecies
        WRITE(FMT3,*) natoms
        WRITE(FMT1,*) '(8x,"Number of elements and atoms found: ",I',
     &   len(trim(FMT2)),',I',len(trim(FMT3)),')'
        print FMT1,  nspecies,natoms 
      end if
      close(21)
      !get name and number of atoms of each species
      if (.not. allocated(species))then
            allocate(species(nspecies))
            species%name=" "
      end if
      
      ispecies=1
      do iatom=1,natoms
        newspecies=.true.
        do i=1,iatom-1
          if (atoms(i)%name(1:2).eq.atoms(iatom)%name(1:2))
     &        newspecies=.false.
        end do
        if(newspecies) then
            species(ispecies)%name=atoms(iatom)%name
            ispecies=ispecies+1
        end if
      end do  

      ! get numbers of atoms of each species
      do ispecies=1,nspecies
            species(ispecies)%howmany=0
            do i=1,natoms
               if(atoms(i)%name(1:2).eq.species(ispecies)%name(1:2)) 
     &            species(ispecies)%howmany=species(ispecies)%howmany+1
            end do
      end do

25    continue
      if(isopen12) then
            write(12,'(8x,A)')'Species found in the file (name, howmany,
     & mass):' 
            do i=1,nspecies
!              write(12,'(8x,A2,1x,I,1x,F10.6)') species(i)%name,
!     &        species(i)%howmany,species(i)%mass         
              WRITE(FMT2,*) species(i)%howmany
              WRITE(FMT1,*) '(8x,A2,1x,I',len(trim(FMT2)),',1x,F10.6)'
              write(12,FMT1) species(i)%name,species(i)%howmany,
     &            species(i)%mass
            end do
      else
            print '(8x,A)','Species found in the file (name, howmany, ma
     &ss):' 
            do i=1,nspecies
!              print '(8x,A2,1x,I,1x,F10.6)', species(i)%name,
!     &        species(i)%howmany,species(i)%mass         
              WRITE(FMT2,*) species(i)%howmany
              WRITE(FMT1,*) '(8x,A2,1x,I',len(trim(FMT2)),',1x,F10.6)'
              print FMT1, species(i)%name,species(i)%howmany,
     &            species(i)%mass
            end do
      end if
      !close(22)
      atoms%occ=1.0D0  ! default occupation
      if(isopen12) then
            write(12,fsubendext) 'read_data'
      else
            print fsubendext, 'read_data'
      end if
      return
      
1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Something wrong with data file." 
      else
            print ferrmssg,"Something wrong with data file." 
      end if
      close(21)
      stop
      end subroutine

c---------------------------------------------------------------------

      subroutine read_out(fincoords,atoms,natoms,species,nspecies,
     &           vecs)

      use defs
      use misc
      implicit none
      
      ! variables
      character(len=*), intent(in) :: fincoords
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms,nspecies
      character line*100,nameat*2
      double precision x(1:3),vecs(1:3,1:3)
      
      ! internal variables
      logical isopen12,newspecies
      integer natom,i,j,iatom,iatom2,ispecies
      integer iele
      character elechar*40,line2*256
      !character FMT1*1024,FMT2*1024,FMT3*1024
      
      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim(adjustl('read_out'))
      else
            print fsubstart, trim(adjustl('read_out'))
      end if
      
      ! read MBPP OUT file file
      open(21,file=fincoords,status='old',err=1000)
      ! get lattice vectors
10    read(21,'(A100)',end=1000,err=1000) line
      line=adjustl(line)
      if (line(1:1).ne.'!'.and.line(1:1).ne.'#') then
            if(index(line,'lattice vectors').gt.0) then
                  read(21,'(A100)',end=1000,err=1000) line
                  do i=1,3
                    read(21,'(A100)',end=1000,err=1000) line
                    read(line(5:100),*) vecs(i,1:3)
                  end do
                  vecs=vecs*bohr
                  print '(8x,"lattice vectors:",3(/8x,3(F20.10)))',vecs
                  rewind(21)
                  goto 11
            end if
      end if
      goto 10

      ! get number of atoms
11    continue
      natoms=0
      !nspecies=0
12    read(21,'(A100)',end=13,err=1000) line
      line=adjustl(line)
      if (line(1:1).ne.'!'.and.line(1:1).ne.'#') then
            if(index(line,'Atomic Positions').gt.0) then
                  !do while (index(line,'parameters').le.0)
                      !read(21,'(A100)',end=15,err=1000) line
                  read(21,'(A100)',err=1000) line
                  do while (index(line,'#').gt.0.or.len_trim(line)==0)
                      read(21,'(A100)',err=1000,end=13) line
                      !print*, line
                      if(index(line,'#').gt.0) natoms=natoms+1
                      !print*,"natoms=",natoms
                  end do
                  !if (.not. allocated(atoms))then
                  !      allocate(atoms(natoms))
                  !      atoms%name=' '
                  !end if
                  !print '(8x,I0," atoms")',natoms
                  goto 13
            end if
      end if
      goto 12

13    continue
      if (.not. allocated(atoms))then
            allocate(atoms(natoms))
            atoms%name=' '
      end if
      ! get atomic names and positions
      rewind(21)
14    read(21,'(A100)',end=15,err=1000) line
      line=adjustl(line)
      if (line(1:1).ne.'!'.and.line(1:1).ne.'#') then
            if(index(line,'Atomic Positions').gt.0) then
                  do iatom=1,natoms
                    read(21,'(A100)',end=15,err=1000) line
                    line=adjustl(line)
                    do while(len_trim(line).le.1) 
                        read(21,'(A100)',end=15,err=1000) line
                        line=adjustl(line)
                    end do
                    read(line(1:3),*) atoms(iatom)%name
                    read(line(9:len(line)),*) atoms(iatom)%where
                    call frac2abs(atoms(i)%where,vecs,atoms(i)%abswhere)
                    atoms(iatom)%core="core"
                    atoms(iatom)%charge=0.0D0
                  end do
                  goto 15
            end if
            goto 14
      end if

15    continue
      close(21)
      ! get number of species
      nspecies=0
      do iatom=1,natoms
         newspecies=.true.
         do iatom2=1,iatom-1
           if (atoms(iatom2)%name(1:4).eq.
     &         atoms(iatom)%name(1:4))newspecies=.false.
         end do
         if(newspecies) nspecies=nspecies+1
      end do
      if (isopen12) then
!        write(12,'(8x,"Number of elements and atoms found: ",I,I)'),
!     &        nspecies,natoms  
        WRITE(FMT2,*) nspecies
        WRITE(FMT3,*) natoms
        WRITE(FMT1,*) '(8x,"Number of elements and atoms found: ",I',
     &   len(trim(FMT2)),',I',len(trim(FMT3)),')'
        write(12,FMT1)  nspecies,natoms 
      else
!        print '(8x,"Number of elements and atoms found: ",I,I)',
!     &        nspecies,natoms  
        WRITE(FMT2,*) nspecies
        WRITE(FMT3,*) natoms
        WRITE(FMT1,*) '(8x,"Number of elements and atoms found: ",I',
     &   len(trim(FMT2)),',I',len(trim(FMT3)),')'
        print FMT1,  nspecies,natoms 
      end if
      !close(21)
      !get name of each species
      if (.not. allocated(species)) then
            allocate(species(nspecies))
            species%name=" "
      end if
      ispecies=1
      do iatom=1,natoms
        newspecies=.true.
        do i=1,iatom-1
          if (atoms(i)%name(1:2).eq.atoms(iatom)%name(1:2))
     &        newspecies=.false.
        end do
        if(newspecies) then
            species(ispecies)%name=atoms(iatom)%name
            ispecies=ispecies+1
        end if
      end do  

      ! get numbers of atoms of each species
      do ispecies=1,nspecies
            species(ispecies)%howmany=0
            do i=1,natoms
               if(atoms(i)%name(1:2).eq.species(ispecies)%name(1:2)) 
     &            species(ispecies)%howmany=species(ispecies)%howmany+1
            end do
      end do

25    continue
      call get_masses_of_atoms(atoms)
      call get_masses_of_species(species)
      if(isopen12) then
            write(12,'(8x,A)')'Species found in the file (name, howmany,
     & mass):' 
            do i=1,nspecies
!              write(12,'(8x,A2,1x,I,1x,F10.6)') species(i)%name,
!     &        species(i)%howmany,species(i)%mass         
              WRITE(FMT2,*) species(i)%howmany
              WRITE(FMT1,*) '(8x,A2,1x,I',len(trim(FMT2)),',1x,F10.6)'
              write(12,FMT1) species(i)%name,species(i)%howmany,
     &            species(i)%mass
            end do
      else
            print '(8x,A)','Species found in the file (name, howmany, ma
     &ss):' 
            do i=1,nspecies
!              print '(8x,A2,1x,I,1x,F10.6)', species(i)%name,
!     &        species(i)%howmany,species(i)%mass         
              WRITE(FMT2,*) species(i)%howmany
              WRITE(FMT1,*) '(8x,A2,1x,I',len(trim(FMT2)),',1x,F10.6)'
              print FMT1, species(i)%name,species(i)%howmany,
     &            species(i)%mass
            end do
      end if
      atoms%occ=1.0D0  ! default occupation
      if(isopen12) then
            write(12,fsubendext) 'read_out'
      else
            print fsubendext, 'read_out'
      end if
      return
      
1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Something wrong with OUT file." 
      else
            print ferrmssg,"Something wrong with OUT file." 
      end if
      close(21)
      return
      end subroutine

c---------------------------------------------------------------------

!      subroutine read_gin(nerr)
!      ! read gulp input/restart file  
!
!      use defs
!      implicit none
!      character infilename*20,outfilename*30,outfilenamec*30,
!     & outfilenames*30,line*100,name*10
!      character xchar*10,ychar*10,zchar*10,coreshel*10
!      double precision a,b,c,alpha,beta,gamm,vecs(1:3,1:3),x,y,z
!      double precision dumrvec(1:3),counter,denominator
!      integer natoms,nerr
!      integer i,i1
!      type(atom) , allocatable :: allatoms(:)
!      type(element), allocatable :: elesc(:),eless(:)
!      integer nelec,neles,numelec,numeles,j,ielec,ieles
!      logical newelec,neweles
!
!      write(12,fsubstart)trim(adjustl('res2coorat'))
!      read(10,'(A20)')infilename
!      read(10,'(A20)')outfilename
!      outfilename=adjustl(outfilename)
!      outfilenamec=outfilename
!      outfilenamec(len_trim(outfilename)+1:30)='.cores'
!      outfilenames=outfilename
!      outfilenames(len_trim(outfilename)+1:30)='.shels'
!      write(12,'(8x,"Inputfile: ",A)')infilename
!      write(12,'(8x,"Outputfile: ",A)') outfilename
!      write(12,*) '      Reading cell parameters...'
!      open(17,file=infilename,status='old')
!      open(18,file=outfilenamec,status='replace')
!      open(19,file=outfilenames,status='replace')
!5     read(17,'(A100)',end=6) line
!      line=adjustl(line)
!      if (index(line,'cell').ge.1.and.line(1:1).ne.'#') goto 7
!      goto 5
!      ! Stop if "cell" is not found
!6     print*,"Error: no cell parameters found."
!      write(12,'(8x,"Error: no cell parameters found.")')
!      nerr=nerr+1
!      close(17)
!      close(18)
!      close(19)
!      return
!7     continue
!      ! read cell parameters
!      read(17,*) a,b,c,alpha,beta,gamm
!      vecs=0.0D0
!      vecs(1,1)=a
!      vecs(2,1)=b*cos(gamm*pi/180.0D0)
!      vecs(2,2)=b*sin(gamm*Pi/180.0D0)
!      vecs(3,1)=c*cos(beta*Pi/180.0D0)
!      vecs(3,2)=c*(cos(alpha*pi/180.0D0)
!     &       -cos(gamm*pi/180.0D0)*cos(beta*pi/180.0D0))
!     &       /sin(gamm*pi/180.0D0)
!      vecs(3,3)=sqrt(c**2-vecs(3,1)**2-vecs(3,2)**2)
!      ! determine number of atoms
!      natoms=0
!8     read(17,'(A100)',end=9) line
!      line=adjustl(line)
!      if(index(line,'frac').ge.1.and.line(1:1).ne.'#') goto 10
!      goto 8
!      write(12,'(8x,"Reading atom info...")')
!      print*,"reading atom info ..."
!      ! Stop if "frac" is not found
!9     print*,"Error: no fractional coordinates found."
!      write(12,'(8x,"Error: no fractional coordiantes found.")')
!      nerr=nerr+1
!      close(17)
!      close(18)
!      close(19)
!      return
!10    read(17,'(A100)',end=15) line 
!      line=adjustl(line)
!      if(line(1:1).ne.'#') then
!        if (index(line,'core').ge.1) natoms=natoms+1
!        if (index(line,'shel').ge.1) natoms=natoms+1
!        if (index(line,'totalenergy').ge.1) goto 15
!        if(index(line,'core').le.0.and.index(line,'shel').le.0)goto 15
!      end if
!      goto 10
!15    continue
!      ! output
!      write(12,'(8x,"There are ",I5," atoms in the file.")')natoms
!      do i=1,3
!        write(18,'("! cell vecs:",3(F20.10))') vecs(i,1:3)
!        write(19,'("! cell vecs:",3(F20.10))') vecs(i,1:3)
!      end do
!      ! read coordinates
!      allocate(allatoms(natoms))
!      rewind(17)
!20    read(17,'(A100)') line  
!      line=adjustl(line)
!      if (index(line,'frac').ge.1.and.line(1:1).ne.'#')goto 30
!      goto 20
!30    continue 
!      do i=1,natoms
!        read(17,'(A100)') line
!        line=adjustl(line)
!        do while (line(1:1).eq.'#') 
!          read(17,'(A100)')line
!          line=adjustl(line)
!        end do 
!        do j=1,len(line)
!          if(line(j:j).eq.'/') line(j:j)='_'
!        end do
!        read(line,*) name,coreshel, xchar,ychar,zchar
!        if(index(xchar,'_').ge.1) then
!            read(xchar(1:index(xchar,'_')-1),*) counter
!            read(xchar(index(xchar,'_')+1:10),*) denominator
!            x=counter/denominator
!        else
!            read(xchar,'(F10.7)') x
!        end if
!        if(index(ychar,'_').ge.1) then
!            read(ychar(1:index(ychar,'_')-1),*) counter
!            read(ychar(index(ychar,'_')+1:10),*) denominator
!            y=counter/denominator
!        else
!            read(ychar,'(F10.7)') y
!        end if
!        if(index(zchar,'_').ge.1) then
!            read(zchar(1:index(zchar,'_')-1),*) counter
!            read(zchar(4:10),'(F7.0)') denominator
!            read(zchar(index(zchar,'_')+1:10),*) denominator
!            z=counter/denominator
!        else
!            read(zchar,'(F10.7)') z
!        end if
!        ! if close to 1, shift by -1
!        if (x.gt.1.0D0-tolfrac(1)) x=x-1.0d0
!        if (y.gt.1.0D0-tolfrac(2)) y=y-1.0d0
!        if (z.gt.1.0D0-tolfrac(3)) z=z-1.0d0
!        !if (index('0123456789 ',name(2:2)).ge.1) name(2:2)='_'
!        ! remove spaces from atomic names and write to file 
!        name(3:6)=name(7:10)
!        name(7:10)='    '
!        allatoms(i)%name=name(1:8)
!        allatoms(i)%core=trim(adjustl(coreshel))
!        allatoms(i)%where(1)=x
!        allatoms(i)%where(2)=y
!        allatoms(i)%where(3)=z
!      end do  
!!      ! sort coordinates
!!      ! sort x coordinates
!!      do i=1,natoms
!!        do i1=i+1,natoms
!!          if (allatoms(i1)%name.eq.allatoms(i)%name.and.
!!     &     allatoms(i1)%where(1).lt.allatoms(i)%where(1)-tolfrac(1))
!!     &    then
!!            dumrvec=allatoms(i)%where
!!            allatoms(i)%where=allatoms(i1)%where
!!            allatoms(i1)%where=dumrvec
!!          end if
!!        end do
!!      end do
!!      ! sort y coordinates
!!      do i=1,natoms
!!        do i1=i+1,natoms
!!          if (allatoms(i1)%name.eq.allatoms(i)%name.and.
!!     &     allatoms(i1)%where(2).lt.allatoms(i)%where(2)-tolfrac(2)
!!     &     .and.abs(allatoms(i1)%where(1)-allatoms(i)%where(1))
!!     &     .lt.tolfrac(1)) then
!!            dumrvec=allatoms(i)%where
!!            allatoms(i)%where=allatoms(i1)%where
!!            allatoms(i1)%where=dumrvec
!!          end if
!!        end do
!!      end do
!!      ! sort z coordinates
!!      do i=1,natoms
!!        do i1=i+1,natoms
!!          if (allatoms(i1)%name.eq.allatoms(i)%name.and.
!!     &     allatoms(i1)%where(3).lt.allatoms(i)%where(3)-tolfrac(3)
!!     &     .and.abs(allatoms(i1)%where(1)-allatoms(i)%where(1))
!!     &     .lt.tolfrac(1)
!!     &     .and.abs(allatoms(i1)%where(2)-allatoms(i)%where(2))
!!     &     .lt.tolfrac(2)) then
!!            dumrvec=allatoms(i)%where
!!            allatoms(i)%where=allatoms(i1)%where
!!            allatoms(i1)%where=dumrvec
!!          end if
!!        end do
!!      end do
!      !! transform to absolute coordinates 
!      !do i=1,natoms
!        !x=allatoms(i)%where(1)
!        !y=allatoms(i)%where(2)
!        !z=allatoms(i)%where(3)
!        !allatoms(i)%where(1:3)=x*vecs(1,1:3)+y*vecs(2,1:3)
!!     &  !        +z*vecs(3,1:3)
!        !write(18,'(A8,I12)') allatoms(i)%name,i
!        !write(18,'(F16.9,2(F20.9))') allatoms(i)%where(1:3)
!      !end do  
!      ! count elements and number of atoms of each element
!      allatoms%written=.false.
!      ! cores and shels
!      nelec=0
!      neles=0
!      do i=1,natoms
!            newelec=.true.
!            neweles=.true.
!            if (allatoms(i)%core(1:4).eq.'shel') newelec=.false.
!            if (allatoms(i)%core(1:4).eq.'core') neweles=.false.
!            do j=1,i-1
!              if (allatoms(i)%name(1:8).eq.allatoms(j)%name(1:8).and.
!     &           allatoms(i)%core(1:4).eq.'core'.and.
!     &           allatoms(j)%core(1:4).eq.'core') newelec=.false.
!              if (allatoms(i)%name(1:8).eq.allatoms(j)%name(1:8).and.
!     &           allatoms(i)%core(1:4).eq.'shel'.and.
!     &           allatoms(j)%core(1:4).eq.'shel') neweles=.false.
!            end do
!            if (newelec) nelec=nelec+1
!            if (neweles) neles=neles+1
!      end do
!      !write(12,'(8x,"Found ",I4," core elements.")')nelec
!      allocate(elesc(1:nelec))
!!      ! shels
!!      neles=0
!!      do i=1,natoms
!!            newele=.true.
!!            if (allatoms(i)%core(1:4).eq.'core') newele=.false.
!!            do j=1,i-1
!!              if (allatoms(i)%name(1:8).eq.allatoms(j)%name(1:8).and.
!!     &          allatoms(i)%core.eq.'shel'
!!     &          .and.allatoms(j)%core.eq.'shel')
!!     &         newele=.false.
!!            end do
!!            if (newele) neles=neles+1
!!      end do
!      !write(12,'(8x,"Found ",I4," shel elements.")')neles
!      allocate(eless(1:neles))
!      ! get names of elements
!      ! cores
!      ielec=1
!      ieles=1
!      do i=1,natoms
!            newelec=.true.
!            neweles=.true.
!            if (allatoms(i)%core(1:4).eq.'shel') newelec=.false.
!            if (allatoms(i)%core(1:4).eq.'core') neweles=.false.
!            do j=1,i-1
!              if (allatoms(i)%name(1:8).eq.allatoms(j)%name(1:8).and.
!     &           allatoms(i)%core(1:4).eq.'core'.and.
!     &           allatoms(j)%core(1:4).eq.'core') newelec=.false.
!              if (allatoms(i)%name(1:8).eq.allatoms(j)%name(1:8).and.
!     &           allatoms(i)%core(1:4).eq.'shel'.and.
!     &           allatoms(j)%core(1:4).eq.'shel') neweles=.false.
!!              if (allatoms(i)%name(1:8).eq.allatoms(j)%name(1:8).and.
!!     &            allatoms(i)%core.eq.'core')
!!     &         newele=.false.
!            end do
!            if (newelec) then
!                  elesc(ielec)%name=allatoms(i)%name
!                  ielec=ielec+1
!            end if
!            if (neweles) then
!                  eless(ieles)%name=allatoms(i)%name
!                  ieles=ieles+1
!            end if
!      end do
!!      ! shels
!!      iele=1
!!      do i=1,natoms
!!            newele=.true.
!!            do j=1,i-1
!!              if (allatoms(i)%name(1:8).eq.allatoms(j)%name(1:8).and.
!!     &            allatoms(i)%core.eq.'shel')
!!     &         newele=.false.
!!            end do
!!            if (newele) then
!!                  eless(iele)%name=allatoms(i)%name
!!                  iele=iele+1
!!            end if
!!      end do
!      ! get numbers of elements
!      ! cores
!      do ielec=1,nelec
!            numelec=0
!            do j=1,natoms
!               if(allatoms(j)%name.eq.elesc(ielec)%name) then
!                  if(allatoms(j)%core.eq.'core')numelec=numelec+1
!               end if
!            end do
!            elesc(ielec)%howmany=numelec
!      end do
!      ! shels
!      do ieles=1,neles
!            numeles=0
!            do j=1,natoms
!               if(allatoms(j)%name.eq.eless(ieles)%name) then
!                  if(allatoms(j)%core.eq.'shel')numeles=numeles+1
!               end if
!            end do
!            eless(ieles)%howmany=numeles
!      end do
!      write(12,'(8x,"Found ",I4," core elements.")')nelec
!      write(12,'(8x,"Elements, number of atoms of each element:")')
!      do i=1,nelec
!            write(12,'(8x,A8,I6)') elesc(i)%name,elesc(i)%howmany
!      end do
!      write(12,'(8x,"Found ",I4," shel elements.")')neles
!      write(12,'(8x,"Elements, number of atoms of each element:")')
!      do i=1,neles
!            write(12,'(8x,A8,I6)') eless(i)%name,eless(i)%howmany
!      end do
!      !write COORAT files
!      !cores
!      do ielec=1,nelec
!            write(18,'(4x,"natom= ",I6," name= ",A8)')
!     &       elesc(ielec)%howmany,elesc(ielec)%name
!            do i=1,natoms
!                  if(allatoms(i)%name.eq.elesc(ielec)%name) then
!                        if(allatoms(i)%core.eq.'core')   
!     &               write(18,'(3(F15.10))') allatoms(i)%where
!                  end if
!            end do
!      end do
!      ! shels
!      do ieles=1,neles
!            write(19,'(4x,"natom= ",I6," name= ",A8)')
!     &       eless(ieles)%howmany,eless(ieles)%name
!            do i=1,natoms
!                  if(allatoms(i)%name.eq.eless(ieles)%name) then
!                        if(allatoms(i)%core.eq.'shel')   
!     &               write(19,'(3(F15.10))') allatoms(i)%where
!                  end if
!            end do
!      end do
!           
!      close(17)
!      close(18)
!      close(19)
!      write(12,*) '     '
!
!      write(12,fsubend)
!      end subroutine

c---------------------------------------------------------------------

      subroutine read_abinit(fincoords,atoms,natoms,species,nspecies,
     &           vecs)

      use defs
      use misc
      implicit none
      
      ! variables
      character(len=*), intent(in) :: fincoords
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms,nspecies
      character line*100,nameat*2
      double precision x(1:3),vecs(1:3,1:3)
      
      ! internal variables
      logical isopen12,newspecies
      integer natom,i,j,iatom,iatom2,ispecies
      integer iele
      character elechar*40,line2*256
      !character FMT1*1024,FMT2*1024,FMT3*1024
      integer, allocatable :: typat(:)
      double precision, allocatable :: znucl(:)
      
      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim(adjustl('read_abinit'))
      else
            print fsubstart, trim(adjustl('read_abinit'))
      end if
      
      ! read abinit output file
      open(21,file=fincoords,status='old',err=1000)
      ! get lattice vectors
      rewind(21)
10    read(21,'(A100)',end=11,err=1001) line
      !line=adjustl(line)
      if (line(1:1).ne.'!'.and.line(1:1).ne.'#') then
        if(index(line,'eal space primitive translations (rprimd)').gt.0)  &
     &  then
              !read(21,'(A100)',end=1000,err=1000) line
              do i=1,3
                read(21,'(A100)',end=1001,err=1001) line
                read(line(1:len(line)),*) vecs(i,1:3)
              end do
              vecs=vecs*bohr
              ! rewind(21)
              ! goto 11
        end if
      end if
      goto 10

      ! get number of atoms
11    continue
      natoms=0
      rewind(21)
12    read(21,'(A100)',end=1002,err=1002) line
      !line=adjustl(line)
      if (line(1:1).ne.'!'.and.line(1:1).ne.'#') then
            if(index(line,'natom =').gt.0) then
                  !do while (index(line,'parameters').le.0)
                      !read(21,'(A100)',end=15,err=1000) line
                      !if(index(line,'#').gt.0) natoms=natoms+1
                  !end do
                  read(line(index(line,'natom =')+7:),*) natoms
                  if (.not. allocated(atoms))then
                        allocate(atoms(natoms))
                        atoms%name=' '
                        atoms%core='core'
                        atoms%charge=0.0D0
                  end if
                  goto 13
            end if
      end if
      goto 12

13    continue
      ! get atomic positions
      rewind(21)
14    read(21,'(A100)',end=1003,err=1003) line
      !line=adjustl(line)
      if (line(1:1).ne.'!'.and.line(1:1).ne.'#') then
            if(index(line,'reduced coordinates').gt.0) then
                  do iatom=1,natoms
                    read(21,*) atoms(iatom)%where
                    call frac2abs(atoms(i)%where,vecs,atoms(i)%abswhere)
                  end do
                  goto 15
            end if
            goto 14
      end if

15    continue
      ! get number of species
      rewind(21)
16    read(21,'(A100)',end=1004,err=1004) line
      !line=adjustl(line)
      if (line(1:1).ne.'!'.and.line(1:1).ne.'#') then
            if(index(line,'ntypat =').gt.0) then
                    read(line(index(line,'ntypat =')+8:),*) nspecies
                    if(.not.allocated(species))
     &                 allocate(species(nspecies))
                  !end do
                  goto 17
            end if
            goto 16
      end if
17    continue
      ! get types of species
      rewind(21)
18    read(21,'(A100)',end=1005,err=1005) line
      !line=adjustl(line)
      if (line(1:1).ne.'!'.and.line(1:1).ne.'#') then
            if(index(line,' typat ').gt.0) then
                  allocate(typat(1:natoms))
                    read(line(18:len(line)),*) typat(1:natoms)
                  goto 19
            end if
            goto 18
      end if
19    continue
      print*, "types read"       
      ! get nuclear charges of species
      rewind(21)
20    read(21,'(A100)',end=1006,err=1006) line
      !line=adjustl(line)
      if (line(1:1).ne.'!'.and.line(1:1).ne.'#') then
            if(index(line,'znucl ').gt.0) then
                  allocate(znucl(1:nspecies))
                    read(line(18:len(line)),*) znucl(1:nspecies)
                  goto 21
            end if
            goto 20
      end if
21    continue 
      !print*, "input read"
      !print*, "znucl=",znucl      
      !print*, "nspecies=",nspecies      
      close(21)
      ! get atom names
      do iatom=1,natoms
         do i=1,nspecies
           if(typat(iatom).eq.i) then
              species(i)%name(1:2)=elements(int(znucl(i)))
              !species(i)%mass=znucl(i)
              atoms(iatom)%name(1:2)=species(i)%name(1:2)
           end if 
         end do
      end do
      !print*, "atoms=",atoms(1:natoms)%name      
      call getspecies(atoms,species)           
      atoms%occ=1.0D0  ! default occupation
      call get_masses_of_species(species)
      call get_masses_of_atoms(atoms)
      if (isopen12) then
!        write(12,'(8x,"Number of elements and atoms found: ",I,I)'),
!     &        nspecies,natoms  
        WRITE(FMT2,*) nspecies
        WRITE(FMT3,*) natoms
        WRITE(FMT1,*) '(8x,"Number of elements and atoms found: ",I',
     &   len(trim(FMT2)),',I',len(trim(FMT3)),')'
        write(12,FMT1)  nspecies,natoms 
      else
!        print '(8x,"Number of elements and atoms found: ",I,I)',
!     &        nspecies,natoms  
        WRITE(FMT2,*) nspecies
        WRITE(FMT3,*) natoms
        WRITE(FMT1,*) '(8x,"Number of elements and atoms found: ",I',
     &   len(trim(FMT2)),',I',len(trim(FMT3)),')'
        print FMT1,  nspecies,natoms 
      end if

25    continue
      if(isopen12) then
            write(12,'(8x,A)')'Species found in the file (name, howmany,
     & mass):' 
            do i=1,nspecies
!              write(12,'(8x,A2,1x,I,1x,F10.6)') species(i)%name,
!     &        species(i)%howmany,species(i)%mass         
              WRITE(FMT2,*) species(i)%howmany
              WRITE(FMT1,*) '(8x,A2,1x,I',len(trim(FMT2)),',1x,F10.6)'
              write(12,FMT1) species(i)%name,species(i)%howmany,
     &            species(i)%mass
            end do
      else
            print '(8x,A)','Species found in the file (name, howmany, ma
     &ss):' 
            do i=1,nspecies
!              print '(8x,A2,1x,I,1x,F10.6)', species(i)%name,
!     &        species(i)%howmany,species(i)%mass         
              WRITE(FMT2,*) species(i)%howmany
              WRITE(FMT1,*) '(8x,A2,1x,I',len(trim(FMT2)),',1x,F10.6)'
              print FMT1, species(i)%name,species(i)%howmany,
     &            species(i)%mass
            end do
      end if
      if(isopen12) then
            write(12,fsubendext) 'read_abinit'
      else
            print fsubendext, 'read_abinit'
      end if
      return
      
1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not open abinit file." 
      else
            print ferrmssg,"Could not open abinit file." 
      end if
      close(21)
      return
1001  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read cell vectors." 
      else
            print ferrmssg,"Could not read cell vectors." 
      end if
      close(21)
      return
1002  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read atom number." 
      else
            print ferrmssg,"Could not read atom number." 
      end if
      close(21)
      return
1003  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read coordinates." 
      else
            print ferrmssg,"Could not read coordinates." 
      end if
      close(21)
      return
1004  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read species number." 
      else
            print ferrmssg,"Could not read species number." 
      end if
      close(21)
      return
1005  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read atom types." 
      else
            print ferrmssg,"Could not read atom types." 
      end if
      close(21)
      return
1006  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read nuclear charges." 
      else
            print ferrmssg,"Could not read nuclear charges." 
      end if
      close(21)
      return
      end subroutine

c---------------------------------------------------------------------

      subroutine read_nwchem(fincoords,atoms,natoms,species,nspecies,
     &           vecs)

      use defs
      use misc
      use writecoords
      implicit none
      
      ! variables
      character(len=*), intent(in) :: fincoords
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms,nspecies
      character line*100,nameat*2
      double precision x(1:3),vecs(1:3,1:3)
      
      ! internal variables
      logical isopen12,newspecies
      integer natom,i,j,iatom,iatom2,ispecies
      integer iele,idum,nstruc
      character elechar*40,line2*256
      !character FMT1*1024,FMT2*1024,FMT3*1024
      character foutcoords*1024
      double precision coords(1:3)
      
      talk=.true.
      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(talk) then
        if(isopen12) then
            write(12,fsubstart) trim(adjustl('read_nwchem'))
        else
            print fsubstart, trim(adjustl('read_nwchem'))
        end if
      end if
      
      ! read abinit output file
      open(21,file=fincoords,status='old',err=1000)
      rewind(21)
!      ! fake lattice vectors in case of a molecule
!      vecs=0.0D0
!      vecs(1,1)=1.0D0
!      vecs(2,2)=1.0D0
!      vecs(3,3)=1.0D0
      ! get lattice vectors
!10    read(21,'(A100)',end=1001,err=1001) line
!      !line=adjustl(line)
!      if (line(1:1).ne.'!'.and.line(1:1).ne.'#') then
!            if(index(line,'space primitive vectors').gt.0) then
!                  !read(21,'(A100)',end=1000,err=1000) line
!                  do i=1,3
!                    read(21,'(A100)',end=1001,err=1001) line
!                    read(line(7:len(line)),*) vecs(i,1:3)
!                  end do
!                  vecs=vecs*bohr
!                  rewind(21)
!                  goto 11
!            end if
!      end if
!      goto 10

      ! get number of atoms
11    continue
      natoms=0
      rewind(21)
12    read(21,'(A100)',end=1002,err=1002) line
      !line=adjustl(line)
      if (line(1:1).ne.'!'.and.line(1:1).ne.'#') then
            if(index(line,'Output coordinates in angstroms').gt.0) then
                  read (21,'(A100)',end=1002,err=1002) line
                  read (21,'(A100)',end=1002,err=1002) line
                  read (21,'(A100)',end=1002,err=1002) line
                  do while (index(line,'Atomic Mass').eq.0) 
                      read(21,'(A100)',end=1002,err=1002) line
                      natoms=natoms+1
                  end do
                  natoms=natoms-2
                  if (.not. allocated(atoms))then
                        allocate(atoms(natoms))
                        atoms%name=' '
                        atoms%core='core'
                        atoms%mass=0.0D0
                  end if
                  goto 13
            end if
      end if
      goto 12

13    continue
      ! get atomic positions and write structure file for each structure found in the file
      nstruc=0
      rewind(21)
14    read(21,'(A100)',end=15,err=1003) line
      !line=adjustl(line)
      if (line(1:1).ne.'!'.and.line(1:1).ne.'#') then
            if(index(line,'Output coordinates in angstroms').gt.0) then
                  read (21,'(A100)',end=1002,err=1002) line
                  read (21,'(A100)',end=1002,err=1002) line
                  read (21,'(A100)',end=1002,err=1002) line
                  do iatom=1,natoms
                    read(21,*) idum,atoms(iatom)%name(1:4),
     &                         atoms(iatom)%charge,coords(1:3)
                    call abs2frac(coords,vecs,
     &                   atoms(iatom)%where)
                    atoms(iatom)%abswhere(1:3)=coords(1:3)
                  end do
                  call getspecies(atoms,species)           
                  nspecies=size(species)
                  WRITE(FMT2,*) nstruc
                  WRITE(FMT1,*) 
     &            '("NWCHEMSTRUC.",I',len(trim(adjustl(FMT2))),
     &             ',".xsf")'
                  foutcoords=' '
                  write(foutcoords,FMT1) nstruc
                  call write_coords(foutcoords,"xsf",atoms,natoms,
     &                 species,nspecies,vecs)
                  nstruc=nstruc+1
                  !goto 15
            end if
            goto 14
      end if

15    continue
!      ! get number of species
!      rewind(21)
!16    read(21,'(A100)',end=1004,err=1004) line
!      !line=adjustl(line)
!      if (line(1:1).ne.'!'.and.line(1:1).ne.'#') then
!            if(index(line,'ntypat =').gt.0) then
!                    read(line(51:len(line)),*) nspecies
!                    if(.not.allocated(species))
!     &                 allocate(species(nspecies))
!                  !end do
!                  goto 17
!            end if
!            goto 16
!      end if
!17    continue
!      ! get types of species
!      rewind(21)
!18    read(21,'(A100)',end=1005,err=1005) line
!      !line=adjustl(line)
!      if (line(1:1).ne.'!'.and.line(1:1).ne.'#') then
!            if(index(line,' typat ').gt.0) then
!                  allocate(typat(1:natoms))
!                    read(line(18:len(line)),*) typat(1:natoms)
!                  goto 19
!            end if
!            goto 18
!      end if
!19    continue
!      print*, "types read"       
!      ! get nuclear charges of species
!      rewind(21)
!20    read(21,'(A100)',end=1006,err=1006) line
!      !line=adjustl(line)
!      if (line(1:1).ne.'!'.and.line(1:1).ne.'#') then
!            if(index(line,'znucl ').gt.0) then
!                  allocate(znucl(1:nspecies))
!                    read(line(18:len(line)),*) znucl(1:nspecies)
!                  goto 21
!            end if
!            goto 20
!      end if
!21    continue 
!      print*, "input read"
!      print*, "znucl=",znucl      
!      print*, "nspecies=",nspecies      
      close(21)
!      ! get atom names
!      do iatom=1,natoms
!         do i=1,nspecies
!           if(typat(iatom).eq.i) then
!              species(i)%name(1:2)=elements(int(znucl(i)))
!              species(i)%mass=znucl(i)
!              atoms(iatom)%name(1:2)=species(i)%name(1:2)
!           end if 
!         end do
!      end do
!      print*, "atoms=",atoms(1:natoms)%name      
      call getspecies(atoms,species)           
      nspecies=size(species)
      call get_masses_of_species(species)
      call get_masses_of_atoms(atoms)
      ! get number of species
      if (isopen12) then
!        write(12,'(8x,"Number of elements and atoms found: ",I,I)'),
!     &        nspecies,natoms  
        !print*, nspecies
        !print*, natoms
        WRITE(FMT2,*) nspecies
        WRITE(FMT3,*) natoms
        WRITE(FMT1,*) '(8x,"Number of elements and atoms found: ",I',
     &   len(trim(FMT2)),',I',len(trim(FMT3)),')'
        write(12,FMT1)  nspecies,natoms 
      else
!        print '(8x,"Number of elements and atoms found: ",I,I)',
!     &        nspecies,natoms  
        WRITE(FMT2,*) nspecies
        WRITE(FMT3,*) natoms
        WRITE(FMT1,*) '(8x,"Number of elements and atoms found: ",I',
     &   len(trim(FMT2)),',I',len(trim(FMT3)),')'
        print FMT1,  nspecies,natoms 
      end if

25    continue
      if(talk) then
        if(isopen12) then
            write(12,'(8x,A)')'Species found in the file (name, howmany,
     &   mass):' 
            do i=1,nspecies
!              write(12,'(8x,A2,1x,I,1x,F10.6)') species(i)%name,
!     &        species(i)%howmany,species(i)%mass         
              WRITE(FMT2,*) species(i)%howmany
              WRITE(FMT1,*) '(8x,A2,1x,I',len(trim(FMT2)),',1x,F10.6)'
              write(12,FMT1) species(i)%name,species(i)%howmany,
     &            species(i)%mass
            end do
        else
            print '(8x,A)','Species found in the file (name, howmany, ma
     &ss):' 
            do i=1,nspecies
!              print '(8x,A2,1x,I,1x,F10.6)', species(i)%name,
!     &        species(i)%howmany,species(i)%mass         
              WRITE(FMT2,*) species(i)%howmany
              WRITE(FMT1,*) '(8x,A2,1x,I',len(trim(FMT2)),',1x,F10.6)'
              print FMT1, species(i)%name,species(i)%howmany,
     &            species(i)%mass
            end do
        end if
      end if
      atoms%occ=1.0D0  ! default occupation
      if(talk) then
        if(isopen12) then
            write(12,fsubendext) 'read_nwchem'
        else
            print fsubendext, 'read_nwchem'
        end if
      end if
      return
      
1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not open nwchem file." 
      else
            print ferrmssg,"Could not open nwchem file." 
      end if
      close(21)
      return
1001  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read cell vectors." 
      else
            print ferrmssg,"Could not read cell vectors." 
      end if
      close(21)
      return
1002  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read atom number." 
      else
            print ferrmssg,"Could not read atom number." 
      end if
      close(21)
      return
1003  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read coordinates." 
      else
            print ferrmssg,"Could not read coordinates." 
      end if
      close(21)
      return
1004  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read species number." 
      else
            print ferrmssg,"Could not read species number." 
      end if
      close(21)
      return
1005  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read atom types." 
      else
            print ferrmssg,"Could not read atom types." 
      end if
      close(21)
      return
1006  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read nuclear charges." 
      else
            print ferrmssg,"Could not read nuclear charges." 
      end if
      close(21)
      return
      end subroutine

c---------------------------------------------------------------------

      subroutine read_xyz(fincoords,atoms,natoms,species,nspecies,
     &           vecs)

      use defs
      use misc
      implicit none
      
      ! variables
      character(len=*), intent(in) :: fincoords
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms,nspecies
      character line*100,nameat*2
      double precision x(1:3),vecs(1:3,1:3)
      
      ! internal variables
      logical isopen12,newspecies
      integer natom,i,j,iatom,iatom2,ispecies
      integer iele,idum
      character elechar*40,line2*256,line3*1024
      !character FMT1*1024,FMT2*1024,FMT3*1024
      double precision coords(1:3)
      character (len=1024), allocatable :: words(:)
      integer nprops
      
      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim(adjustl('read_xyz'))
      else
            print fsubstart, trim(adjustl('read_xyz'))
      end if
      
      open(21,file=fincoords,status='old',err=1000)
      rewind(21)

      ! get number of atoms
      read(21,*,end=1002,err=1002) natoms
      if (.not. allocated(atoms))then
          allocate(atoms(natoms))
      end if

      print'(8x,I0," atoms found")', natoms

      ! get atom names and positions
      read (21,'(A100)',end=1002,err=1002) line
      do iatom=1,natoms
          !
          ! begin check if there are properties to read (forces, ...)
          !
          read (21,'(A100)',end=1002,err=1002) line3
          line3=adjustl(trim(adjustl(line3)))
          call string2words(line3,words)
          nprops=size(words)-4
          if (.not.allocated(atoms(iatom)%properties)) then
            allocate(atoms(iatom)%properties(nprops))
          end if
          atoms(iatom)%name=" "
          read(line3,*) atoms(iatom)%name(1:4),coords(1:3),                  &
     &      atoms(iatom)%properties(1:nprops)
          atoms(iatom)%abswhere=coords
          atoms(iatom)%core="core"
          call abs2frac(coords,vecs,atoms(iatom)%where)
          !
          ! end check if there are properties to read (forces, ...)
          !
      end do
      close(21)
      call getspecies(atoms,species)           
      nspecies=size(species)
      atoms%occ=1.0D0  ! default occupation
      call get_masses_of_atoms(atoms)
      call get_masses_of_species(species)
      if (isopen12) then
        print*, nspecies
        print*, natoms
        WRITE(FMT2,*) nspecies
        WRITE(FMT3,*) natoms
        WRITE(FMT1,*) '(8x,"Number of elements and atoms found: ",I',
     &   len(trim(FMT2)),',I',len(trim(FMT3)),')'
        write(12,FMT1)  nspecies,natoms 
      else
!        print '(8x,"Number of elements and atoms found: ",I,I)',
!     &        nspecies,natoms  
        WRITE(FMT2,*) nspecies
        WRITE(FMT3,*) natoms
        WRITE(FMT1,*) '(8x,"Number of elements and atoms found: ",I',
     &   len(trim(FMT2)),',I',len(trim(FMT3)),')'
        print FMT1,  nspecies,natoms 
      end if

25    continue
      if(isopen12) then
            write(12,'(8x,A)')'Species found in the file (name, howmany,
     & mass):' 
            do i=1,nspecies
!              write(12,'(8x,A2,1x,I,1x,F10.6)') species(i)%name,
!     &        species(i)%howmany,species(i)%mass         
              WRITE(FMT2,*) species(i)%howmany
              WRITE(FMT1,*) '(8x,A2,1x,I',len(trim(FMT2)),',1x,F10.6)'
              write(12,FMT1) species(i)%name,species(i)%howmany,
     &            species(i)%mass
            end do
      else
            print '(8x,A)','Species found in the file (name, howmany, ma
     &ss):' 
            do i=1,nspecies
!              print '(8x,A2,1x,I,1x,F10.6)', species(i)%name,
!     &        species(i)%howmany,species(i)%mass         
              WRITE(FMT2,*) species(i)%howmany
              WRITE(FMT1,*) '(8x,A2,1x,I',len(trim(FMT2)),',1x,F10.6)'
              print FMT1, species(i)%name,species(i)%howmany,
     &            species(i)%mass
            end do
      end if
      if(isopen12) then
            write(12,fsubendext) 'read_xyz'
      else
            print fsubendext, 'read_xyz'
      end if
      return
      
1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not open xyz file." 
      else
            print ferrmssg,"Could not open xyz file." 
      end if
      close(21)
      return
1001  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read cell vectors." 
      else
            print ferrmssg,"Could not read cell vectors." 
      end if
      close(21)
      return
1002  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read atom number." 
      else
            print ferrmssg,"Could not read atom number." 
      end if
      close(21)
      return
1003  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read coordinates." 
      else
            print ferrmssg,"Could not read coordinates." 
      end if
      close(21)
      return
1004  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read species number." 
      else
            print ferrmssg,"Could not read species number." 
      end if
      close(21)
      return
1005  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read atom types." 
      else
            print ferrmssg,"Could not read atom types." 
      end if
      close(21)
      return
1006  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read nuclear charges." 
      else
            print ferrmssg,"Could not read nuclear charges." 
      end if
      close(21)
      return
      end subroutine

c---------------------------------------------------------------------

      subroutine read_outcar(fincoords,atoms,natoms,species,nspecies,
     &           vecs)

      use defs
      use misc, only : getspecies, get_masses_of_atoms,                   &
     &     get_masses_of_species,abs2frac,frac2abs
      implicit none
      
      ! variables
      character(len=*), intent(in) :: fincoords
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms,nspecies
      character line*256,nameat*2
      double precision x(1:3),vecs(1:3,1:3)
      
      ! internal variables
      logical isopen12,newspecies
      integer natom,i,j,iatom,jatom,ispecies
      integer iele,idum
      character elechar*40,line2*256
      !character FMT1*1024,FMT2*1024,FMT3*1024
      double precision coords(1:3)
      
      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim(adjustl('read_outcar'))
      else
            print fsubstart, trim(adjustl('read_outcar'))
      end if
      
      open(21,file=fincoords,status='old',err=1000)

      ! begin get number of atoms
      rewind(21)
      natoms=0
10    read(21,'(A256)',end=1002,err=1002) line
      if (index(line,'NIONS =').gt.0) then
        read(line(66:72),*) natoms
        if (.not. allocated(atoms))then
          allocate(atoms(natoms))
          atoms(:)%core='core'
        end if
        ! BEGIN DEBUG
        !if (talk) print '(8x,"Found ",I0," atoms.")',natoms
        ! END DEBUG
        goto 11 
      end if
      goto 10
      ! end get number of atoms
!
11    continue
      ! begin get number of species
      rewind(21)
      nspecies=0
12    read(21,'(A256)',end=13,err=1002) line
      if (index(line,'TITEL  =').gt.0) then
        nspecies=nspecies+1    
      end if
      goto 12
13    continue
      allocate(species(nspecies))
      species%name=' '
      ! BEGIN DEBUG
      !if (talk) print '(8x,"Found ",I0," species.")',nspecies
      ! END DEBUG
      ! end get number of species
!
      ! begin get species names
      rewind(21)  
      ispecies=0
14    read (21,'(A100)',end=15,err=1002) line
      ! BEGIN DEBUG
      !print*,line
      ! END DEBUG
      if (index(line,'VRHFIN =').gt.0) then
        ispecies=ispecies+1  
        species(ispecies)%name=line(index(line,'=')+1:index(line,':')-1)
        ! BEGIN DEBUG
        !print '(8x,"species name found.")'
        !print '(8x,A)', line(index(line,'=')+1:index(line,':')-1)
        ! END DEBUG
      end if
      !print '(8x,"species names read")'
      goto 14
      ! end get species names
!
      ! begin get atoms per species 
15    continue
      ! BEGIN DEBUG
      !if (talk) print '(8x,"Found these species:")'
      !do ispecies=1,nspecies
      !  if(talk) print '(8x,A2)',species(ispecies)%name(1:2)
      !end do
      ! END DEBUG
      rewind(21)  
16    read (21,'(A100)',end=1002,err=1002) line
      if (index(line,'ions per type =').gt.0) then
        read(line(19:len(line)),*) (species(ispecies)%howmany,            &
     &       ispecies=1,nspecies)
        goto 19
      end if
      !print '(8x,"atoms per species read")'
      goto 16  
      ! end get atoms per species 
!
      ! begin get atom positions 
19    continue
      rewind(21)  
20    read (21,'(A100)',end=21,err=1002) line
      ! vectors
      if (index(line,'direct lattice vectors').gt.0) then
        do i=1,3
          read (21,*,end=1002,err=1002) vecs(i,1:3)
        end do
      end if 
      !print '(8x,"lattice vecs read")'
      ! atom positions
      if (index(line,'POSITION').gt.0.and.index(line,'FORCE').eq.0) then
        !read (21,'(A100)',end=1002,err=1002) line
        read (21,'(A100)',err=1002) line
        jatom=0
        do ispecies=1,nspecies
          do iatom=1,species(ispecies)%howmany
            jatom=jatom+1
            atoms(jatom)%name=species(ispecies)%name
            read (21,*,err=1002,end=1002) atoms(jatom)%abswhere(1:3)
            call abs2frac(atoms(jatom)%abswhere,vecs,atoms(jatom)%where)
          end do
        end do
      end if
      !print '(8x,"absolute positions read")' ! DEBUG
      if (index(line,'POSITION').gt.0.and.index(line,'FORCE').gt.0) then
        !read (21,'(A100)',end=1002,err=1002) line
        read (21,'(A100)',err=1002) line
        jatom=0
        do ispecies=1,nspecies
          do iatom=1,species(ispecies)%howmany
            jatom=jatom+1
            atoms(jatom)%name=species(ispecies)%name
            read (21,*,err=1002,end=1002) atoms(jatom)%abswhere(1:3),     &
     &            atoms(jatom)%force
            call abs2frac(atoms(jatom)%abswhere,vecs,atoms(jatom)%where)
          end do
        end do
      end if
      !print '(8x,"forces read")' ! DEBUG
      if (index(line,'  position of ions in fractional coordinates (dire  &
     &ct lattice)').gt.0) then
        jatom=0
        do ispecies=1,nspecies
          do iatom=1,species(ispecies)%howmany
            jatom=jatom+1
            atoms(jatom)%name=species(ispecies)%name
            read (21,*,err=1002,end=1002) atoms(jatom)%where(1:3)
            call frac2abs(atoms(jatom)%where,vecs,atoms(jatom)%abswhere)
          end do
        end do
      end if
      !print '(8x,"fractional positions read")' ! DEBUG
      !
      ! begin read BEC (if polarization calculation)
      !
      if (index(line,'BORN EFFECTIVE CHARGES').gt.0) then 
        read(21,*) line
        do iatom=1,natoms
          read(21,*) line
          read(21,*) idum, atoms(iatom)%BEC(1,1:3)
          read(21,*) idum, atoms(iatom)%BEC(2,1:3)
          read(21,*) idum, atoms(iatom)%BEC(3,1:3)
        end do
        !print '(8x,"BORN EFFECTIVE CHARGES read")' ! DEBUG
      end if
      !
      ! end read BEC (if polarization calculation)
      !
      goto 20
!
!      do iatom=1,natoms
!          read(21,*) atoms(iatom)%name(1:4),coords(1:3)
!          atoms(iatom)%abswhere=coords
!          call abs2frac(coords,vecs,atoms(iatom)%where)
!      end do
21    continue
      close(21)
!      call getspecies(atoms,species)           
!     nspecies=size(species)
      call get_masses_of_atoms(atoms)
      call get_masses_of_species(species)
      !print '(8x,"masses attached")' ! DEBUG
      if (isopen12) then
        !print '(8x,I0," species")', nspecies
        !print '(8x,I0," atoms")', natoms
        WRITE(FMT2,*) nspecies
        WRITE(FMT3,*) natoms
        WRITE(FMT1,*) '(8x,"Number of elements and atoms found: ",I',
     &   len(trim(FMT2)),',I',len(trim(FMT3)),')'
        write(12,FMT1)  nspecies,natoms 
      else
!        print '(8x,"Number of elements and atoms found: ",I,I)',
!     &        nspecies,natoms  
        WRITE(FMT2,*) nspecies
        WRITE(FMT3,*) natoms
        WRITE(FMT1,*) '(8x,"Number of elements and atoms found: ",I',
     &   len(trim(FMT2)),',I',len(trim(FMT3)),')'
        print FMT1,  nspecies,natoms 
      end if

25    continue
      if(isopen12) then
            write(12,'(8x,A)')'Species found in the file (name, howmany,
     & mass):' 
            do i=1,nspecies
!              write(12,'(8x,A2,1x,I,1x,F10.6)') species(i)%name,
!     &        species(i)%howmany,species(i)%mass         
              WRITE(FMT2,*) species(i)%howmany
              WRITE(FMT1,*) '(8x,A2,1x,I',len(trim(FMT2)),',1x,F10.6)'
              write(12,FMT1) species(i)%name,species(i)%howmany,
     &            species(i)%mass
            end do
      else
            print '(8x,A)','Species found in the file (name, howmany, ma
     &ss):' 
            do i=1,nspecies
!              print '(8x,A2,1x,I,1x,F10.6)', species(i)%name,
!     &        species(i)%howmany,species(i)%mass         
              WRITE(FMT2,*) species(i)%howmany
              WRITE(FMT1,*) '(8x,A2,1x,I',len(trim(FMT2)),',1x,F10.6)'
              print FMT1, species(i)%name,species(i)%howmany,
     &            species(i)%mass
            end do
      end if
      atoms%occ=1.0D0  ! default occupation
      call getspecies(atoms,species)           
      nspecies=size(species)
      if(isopen12) then
            write(12,fsubendext) 'read_outcar'
      else
            print fsubendext, 'read_outcar'
      end if
      return
      
1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not open xyz file." 
      else
            print ferrmssg,"Could not open xyz file." 
      end if
      close(21)
      return
1001  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read cell vectors." 
      else
            print ferrmssg,"Could not read cell vectors." 
      end if
      close(21)
      return
1002  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read atom number." 
      else
            print ferrmssg,"Could not read atom number." 
      end if
      close(21)
      return
1003  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read coordinates." 
      else
            print ferrmssg,"Could not read coordinates." 
      end if
      close(21)
      return
1004  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read species number." 
      else
            print ferrmssg,"Could not read species number." 
      end if
      close(21)
      return
1005  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read atom types." 
      else
            print ferrmssg,"Could not read atom types." 
      end if
      close(21)
      return
1006  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read nuclear charges." 
      else
            print ferrmssg,"Could not read nuclear charges." 
      end if
      close(21)
      return
      end subroutine read_outcar

c---------------------------------------------------------------------

      subroutine read_poscar(fincoords,atoms,natoms,species,nspecies,
     &           vecs)

      use defs
      use misc
      implicit none
      
      ! variables
      character(len=*), intent(in) :: fincoords
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms,nspecies
      character line*256,nameat*2
      double precision x(1:3),vecs(1:3,1:3)
      
      ! internal variables
      logical isopen12,newspecies
      integer i,j,iatom,jatom,ispecies,i1,i2
      integer iele,idum
      character elechar*40,line2*256
      !character FMT1*1024,FMT2*1024,FMT3*1024
      character (len=256) formula
      double precision fac
      character*1024, allocatable :: words(:)
      
      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim(adjustl('read_poscar'))
      else
            print fsubstart, trim(adjustl('read_poscar'))
      end if
      
      open(21,file=fincoords,status='old',err=1000)
      ! begin read cell vectors
      rewind(21)
      read(21,'(A256)',end=1001,err=1001) line
      read(21,*,end=1001,err=1001) fac
      do i=1,3
        read(21,*) vecs(i,1:3)
      end do  
      vecs=vecs*fac
      ! end read cell vectors

      ! begin get number of species and atoms
      nspecies=0
      natoms=0
      read(21,'(A256)',end=1002,err=1002) line
      formula=''
      line=trim(adjustl(line))
      i=1
      j=1
      do while (i<=len(line)-1)
        if(line(i:i).ne.' '.or.line(i+1:i+1).ne.' ') then  
          formula(j:j)=line(i:i)
          j=j+1
        end if
        i=i+1
      end do
      formula=adjustl(formula)
      formula=trim(formula)
      print '(8x,"formula=",A)',trim(formula)
      call string2words(formula,words)
      !nspecies=countsubstring(trim(formula),' ')+1    
      nspecies=size(words)    
      print '(8x,"nspecies=",I0)',nspecies  
      allocate(species(nspecies))
      ! get name and occurrence of elements
      read(21,*) (species(ispecies)%howmany,ispecies=1,nspecies)
      read(formula,*) (species(ispecies)%name,ispecies=1,nspecies)
      !print*,(species(ispecies)%howmany,ispecies=1,nspecies)
      !print*,(species(ispecies)%name,ispecies=1,nspecies)
      natoms=sum(species(:)%howmany)
      !print*,"natoms=",natoms
      allocate(atoms(natoms))
      atoms(:)%core='core'
      iatom=1
      do ispecies=1,nspecies
        atoms(iatom:iatom+species(ispecies)%howmany-1)%name               &
     &       =species(ispecies)%name
        iatom=iatom+species(ispecies)%howmany
      end do
      ! end get number of species and atoms
!
      ! begin get atom positions 
      read (21,'(A100)',end=1003,err=1003) line
      if (index(line,'elective').gt.0) read(21,'(A100)') line
      if (index(line,'irect').gt.0) then
        do iatom=1,natoms 
          read(21,'(A100)') line
          read (line(1:len(line)),*,end=1003,err=1003)                    &
     &       atoms(iatom)%where
          call frac2abs(atoms(iatom)%where,vecs,atoms(iatom)%abswhere)
        end do ! iatom
      end if
      if (index(line,'artesian').gt.0) then
        do iatom=1,natoms 
          read(21,'(A100)') line
          read (line(1:len(line)),*,end=1003,err=1003)                    &
     &     atoms(iatom)%abswhere
        call abs2frac(atoms(iatom)%abswhere,vecs,atoms(iatom)%where)
        end do ! iatom
      end if
      !else
      !  goto 1003 
      !end if
!
      close(21)
      call get_masses_of_atoms(atoms)
      call get_masses_of_species(species)
      !
      ! begin assign atom labels
      !
      atoms(:)%label=" "
      j=1
      atoms(1)%label(1:len_trim(atoms(1)%name))=                          &
     &            trim(atoms(1)%name(1:2))
      write(atoms(1)%label(len_trim(atoms(1)%name)+1:),'(I0)'),j
      do iatom=2,natoms
        if(atoms(iatom)%name(1:2).ne.atoms(iatom-1)%name(1:2)) j=1
        atoms(iatom)%label(1:len_trim(atoms(iatom)%name))=                &
     &  atoms(iatom)%name(1:len_trim(atoms(iatom)%name))
        write(atoms(iatom)%label(len_trim(atoms(iatom)%name)+1:),'(I0)')  &
     &  ,j
        j=j+1
      end do
      atoms(:)%occ=1.0d0
      !
      ! end assign atom labels
      !
      if(isopen12) then
            write(12,fsubendext) 'read_poscar'
      else
            print fsubendext, 'read_poscar'
      end if
      return
      
1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not open file ",fincoords 
      else
            print ferrmssg,"Could not open file ",fincoords 
      end if
      close(21)
      return
1001  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read cell vectors." 
      else
            print ferrmssg,"Could not read cell vectors." 
      end if
      close(21)
      return
1002  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read species info." 
      else
            print ferrmssg,"Could not read species info." 
      end if
      close(21)
      return
1003  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read coordinates." 
      else
            print ferrmssg,"Could not read coordinates." 
      end if
      close(21)
      return
      end subroutine read_poscar

c---------------------------------------------------------------------

      subroutine read_poscar_bec(fincoords,atoms,natoms,species,          &
     &                           nspecies,vecs)

      use defs
      use misc
      implicit none
      
      ! variables
      character(len=*), intent(in) :: fincoords
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms,nspecies
      character line*256,nameat*2
      double precision x(1:3),vecs(1:3,1:3)
      
      ! internal variables
      logical isopen12,newspecies
      integer i,j,iatom,jatom,ispecies,i1,i2
      integer iele,idum
      character elechar*40,line2*256
      character (len=256) formula
      double precision fac
      character*1024, allocatable :: words(:)
      
      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim(adjustl('read_poscar_bec'))
      else
            print fsubstart, trim(adjustl('read_poscar_bec'))
      end if
      
      open(21,file=fincoords,status='old',err=1000)
      ! begin read cell vectors
      rewind(21)
      read(21,'(A256)',end=1001,err=1001) line
      read(21,*,end=1001,err=1001) fac
      do i=1,3
        read(21,*) vecs(i,1:3)
      end do  
      vecs=vecs*fac
      ! end read cell vectors

      ! begin get number of species and atoms
      nspecies=0
      natoms=0
      read(21,'(A256)',end=1002,err=1002) line
      formula=''
      line=trim(adjustl(line))
      i=1
      j=1
      do while (i<=len(line)-1)
        if(line(i:i).ne.' '.or.line(i+1:i+1).ne.' ') then  
          formula(j:j)=line(i:i)
          j=j+1
        end if
        i=i+1
      end do
      formula=adjustl(formula)
      formula=trim(formula)
      print '(8x,"formula=",A)',trim(formula)
      call string2words(formula,words)
      !nspecies=countsubstring(trim(formula),' ')+1    
      nspecies=size(words)    
      print '(8x,"nspecies=",I0)',nspecies  
      allocate(species(nspecies))
      ! get name and occurrence of elements
      read(21,*) (species(ispecies)%howmany,ispecies=1,nspecies)
      read(formula,*) (species(ispecies)%name,ispecies=1,nspecies)
      !print*,(species(ispecies)%howmany,ispecies=1,nspecies)
      !print*,(species(ispecies)%name,ispecies=1,nspecies)
      natoms=sum(species(:)%howmany)
      !print*,"natoms=",natoms
      allocate(atoms(natoms))
      atoms(:)%core='core'
      iatom=1
      do ispecies=1,nspecies
        atoms(iatom:iatom+species(ispecies)%howmany-1)%name               &
     &       =species(ispecies)%name
        iatom=iatom+species(ispecies)%howmany
      end do
      ! end get number of species and atoms
!
      ! begin get atom positions 
      read (21,'(A100)',end=1003,err=1003) line
      if (index(line,'elective').gt.0) read(21,*) line
      if (index(line,'irect').gt.0) then
        do iatom=1,natoms
          read (21,*,end=1003,err=1003) atoms(iatom)%where,               &
     &         (atoms(iatom)%BEC(i,1:3),i=1,3)
          call frac2abs(atoms(iatom)%where,vecs,atoms(iatom)%abswhere)
        end do
      else if (index(line,'artesian').gt.0) then
        do iatom=1,natoms
          read (21,*,end=1003,err=1003) atoms(iatom)%abswhere,            &
     &         (atoms(iatom)%BEC(i,1:3),i=1,3)
          call abs2frac(atoms(iatom)%abswhere,vecs,atoms(iatom)%where)
        end do
      else
        goto 1003 
      end if
!
      close(21)
      call get_masses_of_atoms(atoms)
      call get_masses_of_species(species)
      if(isopen12) then
            write(12,fsubendext) 'read_poscar_bec'
      else
            print fsubendext, 'read_poscar_bec'
      end if
      return
      
1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not open LOCPOT file." 
      else
            print ferrmssg,"Could not open LOCPOT file." 
      end if
      close(21)
      return
1001  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read cell vectors." 
      else
            print ferrmssg,"Could not read cell vectors." 
      end if
      close(21)
      return
1002  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read species info." 
      else
            print ferrmssg,"Could not read species info." 
      end if
      close(21)
      return
1003  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read coordinates." 
      else
            print ferrmssg,"Could not read coordinates." 
      end if
      close(21)
      return
      end subroutine read_poscar_bec

c---------------------------------------------------------------------

      subroutine read_poscar_forces(fincoords,atoms,natoms,species,       &
     &                           nspecies,vecs)

      use defs
      use misc
      implicit none
      
      ! variables
      character(len=*), intent(in) :: fincoords
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms,nspecies
      character line*256,nameat*2
      double precision x(1:3),vecs(1:3,1:3)
      
      ! internal variables
      logical isopen12,newspecies
      integer i,j,iatom,jatom,ispecies,i1,i2
      integer iele,idum
      character elechar*40,line2*256
      character (len=256) formula
      double precision fac
      character*1024, allocatable :: words(:)
      
      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim(adjustl('read_poscar_forces'))
      else
            print fsubstart, trim(adjustl('read_poscar_forces'))
      end if
      
      open(21,file=fincoords,status='old',err=1000)
      ! begin read cell vectors
      rewind(21)
      read(21,'(A256)',end=1001,err=1001) line
      read(21,*,end=1001,err=1001) fac
      do i=1,3
        read(21,*) vecs(i,1:3)
      end do  
      vecs=vecs*fac
      ! end read cell vectors

      ! begin get number of species and atoms
      nspecies=0
      natoms=0
      read(21,'(A256)',end=1002,err=1002) line
      formula=''
      line=trim(adjustl(line))
      i=1
      j=1
      do while (i<=len(line)-1)
        if(line(i:i).ne.' '.or.line(i+1:i+1).ne.' ') then  
          formula(j:j)=line(i:i)
          j=j+1
        end if
        i=i+1
      end do
      formula=adjustl(formula)
      formula=trim(formula)
      print '(8x,"formula=",A)',trim(formula)
      call string2words(formula,words)
      nspecies=size(words)    
      print '(8x,"nspecies=",I0)',nspecies  
      allocate(species(nspecies))
      ! get name and occurrence of elements
      read(21,*) (species(ispecies)%howmany,ispecies=1,nspecies)
      read(formula,*) (species(ispecies)%name,ispecies=1,nspecies)
      natoms=sum(species(:)%howmany)
      allocate(atoms(natoms))
      atoms(:)%core='core'
      iatom=1
      do ispecies=1,nspecies
        atoms(iatom:iatom+species(ispecies)%howmany-1)%name               &
     &       =species(ispecies)%name
        iatom=iatom+species(ispecies)%howmany
      end do
      ! end get number of species and atoms
!
      ! begin get atom positions 
      read (21,'(A100)',end=1003,err=1003) line
      if (index(line,'irect').gt.0) then
        do iatom=1,natoms
          read (21,*,end=1003,err=1003) atoms(iatom)%where,               &
     &         atoms(iatom)%force(1:3)
          call frac2abs(atoms(iatom)%where,vecs,atoms(iatom)%abswhere)
        end do
      else if (index(line,'artesian').gt.0) then
        do iatom=1,natoms
          read (21,*,end=1003,err=1003) atoms(iatom)%abswhere,            &
     &         atoms(iatom)%force(1:3)
          call abs2frac(atoms(iatom)%abswhere,vecs,atoms(iatom)%where)
        end do
      else
        goto 1003 
      end if
!
      close(21)
      call get_masses_of_atoms(atoms)
      call get_masses_of_species(species)
      if(isopen12) then
            write(12,fsubendext) 'read_poscar_forces'
      else
            print fsubendext, 'read_poscar_forces'
      end if
      return
      
1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not open LOCPOT file." 
      else
            print ferrmssg,"Could not open LOCPOT file." 
      end if
      close(21)
      return
1001  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read cell vectors." 
      else
            print ferrmssg,"Could not read cell vectors." 
      end if
      close(21)
      return
1002  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read species info." 
      else
            print ferrmssg,"Could not read species info." 
      end if
      close(21)
      return
1003  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read coordinates." 
      else
            print ferrmssg,"Could not read coordinates." 
      end if
      close(21)
      return
      end subroutine read_poscar_forces

c---------------------------------------------------------------------

      subroutine read_poscar_flags(fincoords,atoms,natoms,species,        &
     &                           nspecies,vecs)

      use defs
      use misc
      implicit none
      
      ! variables
      character(len=*), intent(in) :: fincoords
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms,nspecies
      character line*256,nameat*2
      double precision x(1:3),vecs(1:3,1:3)
      
      ! internal variables
      logical isopen12,newspecies
      integer i,j,iatom,jatom,ispecies,i1,i2
      integer iele,idum
      character elechar*40,line2*256
      character (len=256) formula
      double precision fac
      character*1024, allocatable :: words(:)
      
      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim(adjustl('read_poscar_flags'))
      else
            print fsubstart, trim(adjustl('read_poscar_flags'))
      end if
      
      open(21,file=fincoords,status='old',err=1000)
      ! begin read cell vectors
      rewind(21)
      read(21,'(A256)',end=1001,err=1001) line
      read(21,*,end=1001,err=1001) fac
      do i=1,3
        read(21,*) vecs(i,1:3)
      end do  
      vecs=vecs*fac
      ! end read cell vectors

      ! begin get number of species and atoms
      nspecies=0
      natoms=0
      read(21,'(A256)',end=1002,err=1002) line
      formula=''
      line=trim(adjustl(line))
      i=1
      j=1
      do while (i<=len(line)-1)
        if(line(i:i).ne.' '.or.line(i+1:i+1).ne.' ') then  
          formula(j:j)=line(i:i)
          j=j+1
        end if
        i=i+1
      end do
      formula=adjustl(formula)
      formula=trim(formula)
      print '(8x,"formula=",A)',trim(formula)
      call string2words(formula,words)
      nspecies=size(words)    
      print '(8x,"nspecies=",I0)',nspecies  
      allocate(species(nspecies))
      ! get name and occurrence of elements
      read(21,*) (species(ispecies)%howmany,ispecies=1,nspecies)
      read(formula,*) (species(ispecies)%name,ispecies=1,nspecies)
      natoms=sum(species(:)%howmany)
      allocate(atoms(natoms))
      atoms(:)%core='core'
      iatom=1
      do ispecies=1,nspecies
        atoms(iatom:iatom+species(ispecies)%howmany-1)%name               &
     &       =species(ispecies)%name
        iatom=iatom+species(ispecies)%howmany
      end do
      ! end get number of species and atoms
!
!      ! begin get atom positions 
!      read (21,'(A100)',end=1003,err=1003) line
!      if (index(line,'irect').gt.0) then
!        do iatom=1,natoms
!          read (21,*,end=1003,err=1003) atoms(iatom)%where,               &
!     &         atoms(iatom)%flags(1:3)
!          call frac2abs(atoms(iatom)%where,vecs,atoms(iatom)%abswhere)
!        end do
!      else
!        goto 1003 
!      end if
      ! begin get atom positions 
      read (21,'(A100)',end=1003,err=1003) line
      if (index(line,'elective').gt.0) read(21,'(A100)') line
      if (index(line,'irect').gt.0) then
        do iatom=1,natoms
          read(21,'(A100)') line
          read (line(1:len(line)),*,end=1003,err=1003)                    &
     &         atoms(iatom)%where,atoms(iatom)%flags(1:3)
          call frac2abs(atoms(iatom)%where,vecs,atoms(iatom)%abswhere)
        end do ! iatom
      end if
      if (index(line,'artesian').gt.0) then
        do iatom=1,natoms
          read(21,'(A100)') line
          read (line(1:len(line)),*,end=1003,err=1003)                    &
     &         atoms(iatom)%abswhere,atoms(iatom)%flags(1:3)
          call abs2frac(atoms(iatom)%abswhere,vecs,atoms(iatom)%where)
        end do ! iatom
      end if
!
      close(21)
      call get_masses_of_atoms(atoms)
      call get_masses_of_species(species)
      if(isopen12) then
            write(12,fsubendext) 'read_poscar_flags'
      else
            print fsubendext, 'read_poscar_flags'
      end if
      return
      
1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not open LOCPOT file." 
      else
            print ferrmssg,"Could not open LOCPOT file." 
      end if
      close(21)
      return
1001  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read cell vectors." 
      else
            print ferrmssg,"Could not read cell vectors." 
      end if
      close(21)
      return
1002  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read species info." 
      else
            print ferrmssg,"Could not read species info." 
      end if
      close(21)
      return
1003  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read coordinates." 
      else
            print ferrmssg,"Could not read coordinates." 
      end if
      close(21)
      return
      end subroutine read_poscar_flags

c---------------------------------------------------------------------

      subroutine read_siesta(fincoords,atoms,natoms,species,nspecies,
     &           vecs)

      use defs
      use misc
      implicit none
      
      ! variables
      character(len=*), intent(in) :: fincoords
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms,nspecies
      double precision vecs(1:3,1:3)
      
      ! internal variables
      logical isopen12
      integer i,j,ispecies,iatom
      integer iele
      
      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim(adjustl('read_siesta'))
      else
            print fsubstart, trim(adjustl('read_siesta'))
      end if
      
      open(21,file=fincoords,status='old',err=1000)
      ! begin read cell vectors
      rewind(21)
      do i=1,3
        read(21,*,err=1001,end=1001) vecs(i,1:3)
      end do  
      ! end read cell vectors

      ! begin get number of atoms
      read(21,*) natoms
      print '(8x,"natoms=",I0)',natoms  
      if (allocated(atoms)) deallocate(atoms)
      allocate(atoms(natoms))
      ! get name and position of atoms
      atoms(:)%name=''
      do iatom=1,natoms
        read(21,*, err=1003,end=1003) ispecies, iele, atoms(iatom)%where
        atoms(iatom)%name(1:2)=elements(iele)
        call frac2abs(atoms(iatom)%where,vecs,atoms(iatom)%abswhere)
      end do  
      close(21)
      atoms(:)%core='core'
      atoms(:)%occ=1.0d0
      call getspecies(atoms,species)
      nspecies=size(species)
      ! end get number of species and atoms

      call get_masses_of_atoms(atoms)
      call get_masses_of_species(species)
      !
      ! begin assign atom labels
      !
      atoms(:)%label=" "
      j=1
      atoms(1)%label(1:len_trim(atoms(1)%name))=                          &
     &            trim(atoms(1)%name(1:2))
      write(atoms(1)%label(len_trim(atoms(1)%name)+1:),'(I0)'),j
      do iatom=2,natoms
        if(atoms(iatom)%name(1:2).ne.atoms(iatom-1)%name(1:2)) j=1
        atoms(iatom)%label(1:len_trim(atoms(iatom)%name))=                &
     &  atoms(iatom)%name(1:len_trim(atoms(iatom)%name))
        write(atoms(iatom)%label(len_trim(atoms(iatom)%name)+1:),'(I0)')  &
     &  ,j
        j=j+1
      end do
      !
      ! end assign atom labels
      !
      if(isopen12) then
            write(12,fsubendext) 'read_siesta'
      else
            print fsubendext, 'read_siesta'
      end if
      return
      
1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not open file ",fincoords 
      else
            print ferrmssg,"Could not open file ",fincoords 
      end if
      close(21)
      return
1001  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read cell vectors." 
      else
            print ferrmssg,"Could not read cell vectors." 
      end if
      close(21)
      return
1003  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read coordinates." 
      else
            print ferrmssg,"Could not read coordinates." 
      end if
      close(21)
      return
      end subroutine read_siesta

c---------------------------------------------------------------------

      subroutine read_siesta_out(fincoords,atoms,natoms,species,          &
     &    nspecies,vecs)

      use defs
      use misc, only : getspecies, get_masses_of_atoms,                   &
     &     get_masses_of_species,frac2abs
      implicit none
      
      ! variables
      character(len=*), intent(in) :: fincoords
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms,nspecies
      character line*256,nameat*2
      double precision x(1:3),vecs(1:3,1:3)
      
      ! internal variables
      logical isopen12,newspecies
      integer natom,i,j,iatom,jatom,ispecies
      integer iele,idum
      character elechar*40,line2*256
      double precision coords(1:3)
      
      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim(adjustl('read_siesta_out'))
      else
            print fsubstart, trim(adjustl('read_siesta_out'))
      end if
      
      open(21,file=fincoords,status='old',err=1000)

      ! begin get number of atoms
      if(allocated(atoms)) deallocate(atoms)
      rewind(21)
      natoms=0
10    read(21,'(A256)',end=1002,err=1002) line
      if (index(line,'NumberOfAtoms').gt.0) then
        read(line(14:72),*) natoms
        if (.not. allocated(atoms))then
          allocate(atoms(natoms))
          atoms(:)%core='core'
          atoms(:)%occ=1.0d0
          atoms(:)%name=''
        end if
        if(talk) print '(8x,I0," atoms")',natoms
        goto 11 
      end if
      goto 10
      ! end get number of atoms
!
11    continue
      ! begin get number of species
      rewind(21)
      nspecies=0
12    read(21,'(A256)',end=13,err=1004) line
      if (index(line,'NumberOfSpecies').gt.0) then
        read(line(16:72),*) nspecies
        if(talk) print '(8x,I0," species")',nspecies
        goto 13
      end if
      goto 12

13    continue
      allocate(species(nspecies))
      ! end get number of species
!
      ! begin get species names
      rewind(21)  
      species(:)%name=''
14    read (21,'(A100)',end=19,err=1005) line
      if (index(line,'%block ChemicalSpeciesLabel').gt.0) then
        do ispecies=1,nspecies
          read(21,*,end=19,err=1005) idum,iele,                           &
     &      species(ispecies)%name(1:2)
            if (talk) print '(8x,A2)',species(ispecies)%name(1:2)
        end do
        goto 19
      end if
      goto 14
      ! end get species names

      ! begin get atom positions 
19    continue
      rewind(21)  
20    read (21,'(A100)',end=21,err=1001) line
      ! vectors
      if (index(line,'outcell: Unit cell vectors').gt.0) then
        !if (talk) print '(8x,"cell vectors:")'
        do i=1,3
          read (21,*,end=1002,err=1001) vecs(i,1:3)
          !if(talk) print '(8x,3(F12.6))',vecs(i,1:3)
        end do
      end if 
      ! atom positions
      if (index(line,'outcoor:').gt.0.and.index(line,'fractional').gt.0)  &
     &then
        !if(talk)print '(8x,"fractional positions:")' ! DEBUG
        do iatom=1,natoms
          read (21,*,err=1003,end=1003) atoms(iatom)%where(1:3),          &
     &          ispecies,idum,atoms(iatom)%name(1:2)
            !if (talk) print '(8x,3(F12.6))',atoms(iatom)%where(1:3)
        end do
      end if
      if (index(line,'siesta: Atomic forces').gt.0) then
        !if(talk)print '(8x,"forces:")' ! DEBUG
        do iatom=1,natoms
            read (21,*,err=1006,end=1006) idum,                           &
     &            atoms(iatom)%force(1:3)
            !if (talk) print '(8x,3(F12.6))',atoms(iatom)%force(1:3)
        end do
      end if
      !
      goto 20
!
21    continue
      close(21)
      call getspecies(atoms,species)           
      nspecies=size(species)
      call get_masses_of_atoms(atoms)
      call get_masses_of_species(species)
      print '(8x,"masses attached")' ! DEBUG
      do iatom=1,natoms
        call frac2abs(atoms(iatom)%where,vecs,atoms(iatom)%abswhere)
      end do
      if (isopen12) then
        WRITE(FMT2,*) nspecies
        WRITE(FMT3,*) natoms
        WRITE(FMT1,*) '(8x,"Number of elements and atoms found: ",I',
     &   len(trim(FMT2)),',I',len(trim(FMT3)),')'
        write(12,FMT1)  nspecies,natoms 
      else
        WRITE(FMT2,*) nspecies
        WRITE(FMT3,*) natoms
        WRITE(FMT1,*) '(8x,"Number of elements and atoms found: ",I',
     &   len(trim(FMT2)),',I',len(trim(FMT3)),')'
        print FMT1,  nspecies,natoms 
      end if

25    continue
      if(isopen12) then
            write(12,'(8x,A)')'Species found in the file (name, howmany,
     & mass):' 
            do i=1,nspecies
              WRITE(FMT2,*) species(i)%howmany
              WRITE(FMT1,*) '(8x,A2,1x,I',len(trim(FMT2)),',1x,F10.6)'
              write(12,FMT1) species(i)%name,species(i)%howmany,
     &            species(i)%mass
            end do
      else
            print '(8x,A)','Species found in the file (name, howmany, ma
     &ss):' 
            do i=1,nspecies
              WRITE(FMT2,*) species(i)%howmany
              WRITE(FMT1,*) '(8x,A2,1x,I',len(trim(FMT2)),',1x,F10.6)'
              print FMT1, species(i)%name,species(i)%howmany,
     &            species(i)%mass
            end do
      end if
      if(isopen12) then
            write(12,fsubendext) 'read_siesta_out'
      else
            print fsubendext, 'read_siesta_out'
      end if
      return
      
1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not open siesta output file." 
      else
            print ferrmssg,"Could not open siesta output file." 
      end if
      close(21)
      return
1001  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read cell vectors." 
      else
            print ferrmssg,"Could not read cell vectors." 
      end if
      close(21)
      return
1002  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read atom number." 
      else
            print ferrmssg,"Could not read atom number." 
      end if
      close(21)
      return
1003  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read coordinates." 
      else
            print ferrmssg,"Could not read coordinates." 
      end if
      close(21)
      return
1004  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read species number." 
      else
            print ferrmssg,"Could not read species number." 
      end if
      close(21)
      return
1005  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read atom types." 
      else
            print ferrmssg,"Could not read atom types." 
      end if
      close(21)
      return
1006  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read forces." 
      else
            print ferrmssg,"Could not read forces." 
      end if
      close(21)
      return
      end subroutine read_siesta_out

c---------------------------------------------------------------------

      subroutine read_cif(fincoords,atoms,natoms,species,nspecies,
     &           vecs)

      use defs
      use misc
      implicit none
      
      ! variables
      character(len=*), intent(in) :: fincoords
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms,nspecies
      character line*256,nameat*2
      double precision x(1:3),vecs(1:3,1:3)
      
      ! internal variables
      logical isopen12,newspecies
      integer i,j,iatom,jatom,ispecies,i1,i2
      integer iele,idum
      character elechar*40,line2*256
      !character FMT1*1024,FMT2*1024,FMT3*1024
      character (len=256) formula
      double precision cellpars(1:6)
      
      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim(adjustl('read_cif'))
      else
            print fsubstart, trim(adjustl('read_cif'))
      end if
      
      open(21,file=fincoords,status='old',err=1000)

11    continue
      ! begin get number of species and atoms
      rewind(21)
      nspecies=0
      natoms=0
12    read(21,'(A256)',end=13,err=1002) line
      if (index(line,'_chemical_formula_sum').gt.0) then
        formula=''
        formula=''
        j=index(line,'_chemical_formula_sum')                             &
     &        +len_trim('_chemical_formula_sum')
        formula=line(j:len_trim(line))
        do i=1,len_trim(formula)
          if (formula(i:i).eq."'") formula(i:i)=" "
        end do
        formula=trim(adjustl(formula))
        print '(8x,"chemical formula sum:",A)',formula
        !nspecies=countsubstring(line(42:len_trim(line)-1),' ')+1    
        !print '("line=",A)',line(42:len_trim(line)-1)
        !print '("formula=",A)',formula
        nspecies=countsubstring(trim(adjustl(formula)),' ')+1    
        allocate(species(nspecies))
        ! get name and occurrence of elements
        i1=1
        i2=1
        ispecies=1
        do while(i2<=len_trim(formula))
          do while (index('1234567890',formula(i2:i2)).le.0               &
     &              .and.formula(i2:i2).ne.' ')
            i2=i2+1 
          end do
          species(ispecies)%name=formula(i1:i2-1)
          i1=i2
          do while(index('1234567890',formula(i2:i2)).gt.0) 
            i2=i2+1  
          end do 
          if(formula(i1:i2-1).eq.'') then
            species(ispecies)%howmany=1  
          else 
            read(formula(i1:i2-1),*) species(ispecies)%howmany
          end if
          natoms=natoms+species(ispecies)%howmany
          ispecies=ispecies+1
          i1=i2+1
          i2=i2+1 
        end do  
        goto 13
      end if
      goto 12
13    continue
      print '(8x,I0,x,"species")',nspecies
      print '(8x,I0,x,"atoms")',natoms
      allocate(atoms(natoms))
      ! end get number of species and atoms
!
      ! begin get atom positions 
19    continue
      rewind(21)  
20    read (21,'(A100)',end=21,err=1002) line
      ! vectors
      if (index(line,'_cell_length_a').gt.0) then
         read (line(index(line,'_cell_length_a')+15:len(line)),*,         & 
     &         end=1002,err=1002) cellpars(1)
      end if 
      if (index(line,'_cell_length_b').gt.0) then
         read (line(index(line,'_cell_length_b')+15:len(line)),*,         & 
     &         end=1002,err=1002) cellpars(2)
      end if 
      if (index(line,'_cell_length_c').gt.0) then
         read (line(index(line,'_cell_length_c')+15:len(line)),*,         & 
     &         end=1002,err=1002) cellpars(3)
      end if 
      if (index(line,'_cell_angle_alpha').gt.0) then
         read (line(index(line,'_cell_length_alpha')+19:len(line)),*,     & 
     &         end=1002,err=1002) cellpars(4)
      end if 
      if (index(line,'_cell_angle_beta').gt.0) then
         read (line(index(line,'_cell_length_beta')+18:len(line)),*,      & 
     &         end=1002,err=1002) cellpars(5)
      end if 
      if (index(line,'_cell_angle_gamma').gt.0) then
         read (line(index(line,'_cell_length_gamma')+19:len(line)),*,     & 
     &         end=1002,err=1002) cellpars(6)
      end if 
      goto 20
21    continue
      print '(8x,"cell parameters: ",6(F12.6,1x))',cellpars
      call cellpars2vecs(cellpars,vecs)
      ! atom positions
      rewind(21)
22    read(21,'(A256)',err=1002,end=23) line 
      if (index(line,'_atom_site_occupancy').gt.0) then
        do iatom=1,natoms
          read (21,*,end=1002,err=1002) atoms(iatom)%name,elechar,        &         
     &      idum,atoms(iatom)%where,atoms(iatom)%occ
            call frac2abs(atoms(iatom)%where,vecs,atoms(iatom)%abswhere)
        end do
      end if
      goto 22
23    close(21)
      if(isopen12) then
            write(12,fsubendext) 'read_cif'
      else
            print fsubendext, 'read_cif'
      end if
      return
      
1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not open cif file." 
      else
            print ferrmssg,"Could not open cif file." 
      end if
      close(21)
      return
1001  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read cell vectors." 
      else
            print ferrmssg,"Could not read cell vectors." 
      end if
      close(21)
      return
1002  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read atom number." 
      else
            print ferrmssg,"Could not read atom number." 
      end if
      close(21)
      return
1003  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read coordinates." 
      else
            print ferrmssg,"Could not read coordinates." 
      end if
      close(21)
      return
1004  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read species number." 
      else
            print ferrmssg,"Could not read species number." 
      end if
      close(21)
      return
1005  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read atom types." 
      else
            print ferrmssg,"Could not read atom types." 
      end if
      close(21)
      return
1006  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read nuclear charges." 
      else
            print ferrmssg,"Could not read nuclear charges." 
      end if
      close(21)
      return
      end subroutine read_cif

c---------------------------------------------------------------------

      subroutine read_vesta(fincoords,atoms,natoms,species,nspecies,
     &           vecs)

      use defs
      use misc,only : cellpars2vecs,getspecies
      implicit none
      
      ! variables
      character(len=*), intent(in) :: fincoords
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms,nspecies
      character line*256,nameat*2
      double precision x(1:3),vecs(1:3,1:3)
      
      ! internal variables
      logical isopen12,newspecies
      integer i,j,iatom,jatom,ispecies,i1,i2
      integer iele,idum
      character elechar*40,line2*256,chardum*24
      !character FMT1*1024,FMT2*1024,FMT3*1024
      character (len=256) formula
      double precision cellpars(1:6)
      
      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim(adjustl('read_vesta'))
      else
            print fsubstart, trim(adjustl('read_vesta'))
      end if
      
      open(21,file=fincoords,status='old',err=1000)
      rewind(21)
      nspecies=0

      ! begin read cell parameters
12    read(21,'(A256)',end=1001,err=1001) line
      if (index(line,'CELLP').gt.0) then
        read(21,*,end=13,err=1001) cellpars(1:6)
        call cellpars2vecs(cellpars,vecs)
        goto 13
      end if
      goto 12
      ! end read cell parameters

13    continue      
      print '(8x,"cell vecs:",3(3(F10.6),1x))',vecs

      ! begin get # of atoms 
      natoms=0
14    read(21,'(A256)',end=1002,err=1002) line
      if (index(line,'STRUC').gt.0) then
        do while (index(line,'THERI').le.0)
          read(21,'(A256)',end=1002,err=1002) line
          natoms=natoms+1
        end do
        goto 15
      end if
      goto 14
      ! end read # of atoms
15    continue      
      natoms=(natoms-2)/2
      print '(8x,"# of atoms:",I0)',natoms
      allocate(atoms(natoms))
     
      ! begin get atom positions  
      rewind(21)
16    read(21,'(A256)',end=1003,err=1003) line
      if (index(line,'STRUC').gt.0) then
        do iatom=1,natoms
          atoms(iatom)%name=" "
          atoms(iatom)%label=" "
          read(21,*,end=1003,err=1003) idum,                              &
     &         atoms(iatom)%name,                                         &
     &      atoms(iatom)%label,atoms(iatom)%occ,atoms(iatom)%where
          read(21,'(A256)',end=1003,err=1003) line
          !print*, atoms(iatom)%label,atoms(iatom)%occ 
        end do
        goto 17
      end if
      goto 16
      ! end read positions of atoms

17    continue      
      call getspecies(atoms,species)
      nspecies=size(species)

23    close(21)
      if(isopen12) then
            write(12,fsubendext) 'read_vesta'
      else
            print fsubendext, 'read_vesta'
      end if
      return
      
1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not open vesta file." 
      else
            print ferrmssg,"Could not open vesta file." 
      end if
      close(21)
      return
1001  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read cell vectors." 
      else
            print ferrmssg,"Could not read cell vectors." 
      end if
      close(21)
      return
1002  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read atom number." 
      else
            print ferrmssg,"Could not read atom number." 
      end if
      close(21)
      return
1003  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read coordinates." 
      else
            print ferrmssg,"Could not read coordinates." 
      end if
      close(21)
      return
1004  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read species number." 
      else
            print ferrmssg,"Could not read species number." 
      end if
      close(21)
      return
1005  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read atom types." 
      else
            print ferrmssg,"Could not read atom types." 
      end if
      close(21)
      return
1006  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read nuclear charges." 
      else
            print ferrmssg,"Could not read nuclear charges." 
      end if
      close(21)
      return
      end subroutine read_vesta

c---------------------------------------------------------------------

      subroutine read_findsym(fincoords,atoms,natoms,species,nspecies,
     &           vecs)

      use defs
      use misc
      implicit none
      
      ! variables
      character(len=*), intent(in) :: fincoords
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms,nspecies
      double precision vecs(1:3,1:3)
      
      ! internal variables
      logical isopen12,read_wyck
      integer i,j,iatom,jatom,i_in_ele,idum
      character line*256
      character line2*256
      character (len=256) formula
      double precision cellpars(1:6)
      integer, allocatable :: elenums(:)
      character*12 wyckstring
      
      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim(adjustl('read_findsym'))
      else
            print fsubstart, trim(adjustl('read_findsym'))
      end if
      
      open(21,file=fincoords,status='old',err=1000)

11    continue
      ! begin get number of species and atoms
      rewind(21)
      nspecies=0
      natoms=0
12    read(21,'(A256)',end=13,err=1002) line
      !
      ! begin get lattice parameters
      !
      if (index(line,'Lattice parameters').gt.0) then
        read (21,*,end=1002,err=1002) cellpars(1:6)
        call cellpars2vecs(cellpars,vecs)
      end if 
      !
      ! end get lattice parameters
      !
      !
      ! begin read chemical formula (not needed)
      !
      if (index(line,'Brigham Young University').gt.0) then
        read(21,'(A256)') line2
        formula=''
        read(21,'(A256)') formula
        formula=trim(adjustl(formula))
        print '(8x,"chemical formula sum:",A)',formula
      end if
      !
      ! end read chemical formula (not needed)
      !
      if (index(line,'Number of atoms in unit cell').gt.0) then
        !
        ! begin read number of atoms
        !
        read(21,*) natoms
        print '(8x,I0," atoms")',natoms
        allocate(atoms(natoms))
        atoms(:)%name=' '
        atoms(:)%core='core'
        !
        ! end read number of atoms
        !
        !
        ! begin read atom types 
        !
        read(21,'(A256)') line2
        iatom=1
        allocate(elenums(20))
        do while (iatom<natoms)
          read(21,'(A256)') line2
          read(line2,'(20(I3))') elenums(1:20) 
          atoms(iatom:iatom+19)%name(1:2)=elements(elenums(1:20))
          iatom=iatom+20
        end do
        call getspecies(atoms,species)
        nspecies=size(species)
        if (isopen12) then
          WRITE(FMT2,*) nspecies
          WRITE(FMT3,*) natoms
          WRITE(FMT1,*) '(8x,"Number of elements and atoms found: ",I',
     &     len(trim(FMT2)),',I',len(trim(FMT3)),')'
          write(12,FMT1)  nspecies,natoms 
        else
          WRITE(FMT2,*) nspecies
          WRITE(FMT3,*) natoms
          WRITE(FMT1,*) '(8x,"Number of elements and atoms found: ",I',
     &     len(trim(FMT2)),',I',len(trim(FMT3)),')'
          print FMT1, nspecies,natoms 
        end if
        !
      end if
      !
      if (index(line,'Position of each atom').gt.0) then
        do iatom=1,natoms
          read(21,*) idum,atoms(iatom)%where
          call frac2abs(atoms(iatom)%where,vecs,atoms(iatom)%abswhere)
        end do
        goto 13
      end if
      goto 12
      !
13    continue
      print '(8x,I0,x,"species")',nspecies
      rewind(21)
      read_wyck=.false.
      !
14    read(21,'(A256)', end=15) line
      !
      ! begin read Wyckoff positions
      !
      if (index(line,'Atomic positions in terms of a,b,c:').gt.0) then
        read_wyck=.true.
        print '(8x,"Element, atom number, atom number within its species  &
     &, Wyckoff position:")'
      end if
      !
      if (read_wyck) then
        if (index(line,'Atomic positions in terms of a,b,c:').le.0) then
          !
          if (index(line,'Wyckoff').gt.0) then 
            read(line(18:),*) wyckstring 
          else
            read(line(1:len(line)),*) iatom 
            i_in_ele=num_in_ele(atoms,iatom)
            print '(8x, A2,1x,2(I10,1x),A12)', atoms(iatom)%name,iatom,   &
     &              i_in_ele, wyckstring
          end if
        end if 
        !
      end if
      !
      goto 14
      !
      ! end read Wyckoff positions
      !
15    continue
      close(21)
      call get_masses_of_atoms(atoms)
      call get_masses_of_species(species)
      !
      !
      if(isopen12) then
            write(12,fsubendext) 'read_findsym'
      else
            print fsubendext, 'read_findsym'
      end if
      return
      
1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not open findsym file." 
      else
            print ferrmssg,"Could not open findsym file." 
      end if
      close(21)
      return
1001  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read cell vectors." 
      else
            print ferrmssg,"Could not read cell vectors." 
      end if
      close(21)
      return
1002  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read atom number." 
      else
            print ferrmssg,"Could not read atom number." 
      end if
      close(21)
      return
1003  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read coordinates." 
      else
            print ferrmssg,"Could not read coordinates." 
      end if
      close(21)
      return
1005  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Could not read atom types." 
      else
            print ferrmssg,"Could not read atom types." 
      end if
      close(21)
      return
      !
      end subroutine read_findsym

c---------------------------------------------------------------------

      end module


