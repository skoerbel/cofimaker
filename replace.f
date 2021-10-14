      module replace
      ! contains routines that replace, remove, or add atoms

      implicit none

      contains 

c---------------------------------------------------------------------

      subroutine rmatoms(infile,informat,thespecies,nrm,
     &           outputfile)

      use defs
      use readcoords
      use misc
      use writecoords
      implicit none
      character(len=*) infile,informat,thespecies,outputfile
      integer nrm
      ! local
      type(atom),allocatable :: atoms(:),newatoms(:)
      type(element),allocatable :: species(:),newspecies(:)
      integer natoms,nspecies,nnewspecies,iatom,ispecies,nthere
      integer inewatom,ithespecies,n,n1
      integer, allocatable :: vacs(:)
      double precision vecs(1:3,1:3)
      logical isopen12,isdrawn
      !character FMT1*1024,FMT2*1024

      talk=.true.
      ! write hello
      inquire(unit=12,opened=isopen12)
      if(talk) then
        if(isopen12) then
            write(12,fsubstart)'rmatoms'
        else
            print fsubstart, 'rmatoms'
        end if
      end if
      thespecies=adjustl(thespecies)

      ! read the coordinate file
      call read_coords(infile,informat,atoms,natoms,species,
     &                       nspecies,vecs)
      
      ! get the total number of atoms of
      ! the species of which you want to remove atoms
      do ispecies=1,nspecies
        if(species(ispecies)%name(1:2).eq.thespecies(1:2)) 
     &     nthere=species(ispecies)%howmany
      end do
      if(talk) then
        if(isopen12) then
!            write(12,'(8x,"There are ",I,1x,A2," atoms altogether.")')
!     &            nthere,thespecies(1:2)       
            WRITE(FMT2,*) nthere
            WRITE(FMT1,*) '(8x,"There are ",I',
     &       len(trim(FMT2)),',1x,A2," atoms altogether.")'
            write(12,FMT1)  nthere,thespecies(1:2) 
!            write(12,'(8x,I," of them will be removed.")')nrm       
            WRITE(FMT2,*) nrm
            WRITE(FMT1,*) '(8x,I',len(trim(FMT2)),'" of them will be rem
     &oved.")'
            write(12,FMT1) nrm 
        else
!            print '(8x,"There are ",I,1x,A2," atoms altogether.")',
!     &            nthere,thespecies(1:2)       
            WRITE(FMT2,*) nthere
            WRITE(FMT1,*) '(8x,"There are ",I',
     &       len(trim(FMT2)),',1x,A2," atoms altogether.")'
            print FMT1, nthere,thespecies(1:2) 
!            print '(8x,I," of them will be removed.")',nrm       
            WRITE(FMT2,*) nrm
            WRITE(FMT1,*) '(8x,I',len(trim(FMT2)),'" of them will be rem
     &oved.")'
            print FMT1, nrm 
        end if
      end if

      ! 
      ! if positive number to be removed
      !
      if (nrm.ge.0) then
        !
        ! call subroutine that draws balls from an urn
        call drawballs(nthere,nrm,vacs)

        ! get new structure
        allocate(newatoms(natoms-nrm))
        inewatom=0
        ithespecies=0
        do iatom=1,natoms
          if (atoms(iatom)%name(1:2).ne.thespecies(1:2)) then
              inewatom=inewatom+1
              newatoms(inewatom)=atoms(iatom)
          end if
          if (atoms(iatom)%name(1:2).eq.thespecies(1:2)) then
             ithespecies=ithespecies+1
             if(.not.any(vacs.eq.ithespecies)) then
               inewatom=inewatom+1
               newatoms(inewatom)=atoms(iatom)
             end if
             !if (atoms(iatom)%name(1:2).eq.thespecies(1:2)) then
             !  ithespecies=ithespecies+1
             !  if(.not.any(vacs.eq.ithespecies)) then
             !    inewatom=inewatom+1
             !    newatoms(inewatom)=atoms(iatom)
             !  end if
             !end if
          end if
        end do
        ! 
        ! end if positive number to be removed
        !
      end if
      !
      ! if negative number to be removed
      !
      if (nrm.lt.0) then
        allocate(vacs(1))
        allocate(newatoms(natoms-1))
        inewatom=0
        ithespecies=0
        do iatom=1,natoms
          if (atoms(iatom)%name(1:2).ne.thespecies(1:2)) then
                inewatom=inewatom+1
                newatoms(inewatom)=atoms(iatom)
          end if
          if (atoms(iatom)%name(1:2).eq.thespecies(1:2)) then
                ithespecies=ithespecies+1
                if(-nrm.ne.ithespecies) then
                  inewatom=inewatom+1
                  newatoms(inewatom)=atoms(iatom)
                else
                  vacs(1)=iatom
                end if
          end if
        end do
      end if
      !
      ! end if negative number to be removed
      !

      ! get new numbers of species
      call getspecies(newatoms,newspecies)
      nnewspecies=size(newspecies)
      call write_coords(outputfile,informat,newatoms,size(newatoms),      &
     &     newspecies,nnewspecies,vecs)

      ! write a file with the atoms that were removed
      open(51,file="VACS.DAT",status="replace")
      write(51,'("# The ",A," atoms that were removed are:")')
     &  trim(thespecies)
      do n1=1,max(nrm,1)
        write(51,*) vacs(n1)
      end do
      close(51)

      ! write goodbye
      if(talk) then
        if(isopen12) then
            write(12,fsubendext)'rmatoms'
        else
            print fsubendext, 'rmatoms'
        end if
      end if
      end subroutine

c---------------------------------------------------------------------

      subroutine rpatoms(infile,informat,thespecies,nrp,
     &           theotherspecies,outputfile)

      use defs
      use readcoords
      use misc
      use writecoords
      implicit none
      character(len=*) infile,informat,thespecies,theotherspecies,
     &  outputfile
      integer nrp
      ! local
      type(atom),allocatable :: atoms(:)!,newatoms(:)
      type(element),allocatable :: species(:),newspecies(:)
      integer natoms,nspecies,nnewspecies,iatom,ispecies,nthere
      integer ithespecies,n,n1 !,inewatom
      integer, allocatable :: replaced(:)
      double precision vecs(1:3,1:3)
      logical isopen12,isdrawn
      !character FMT1*1024,FMT2*1024

      talk=.true.
      ! write hello
      inquire(unit=12,opened=isopen12)
      if(talk) then
        if(isopen12) then
            write(12,fsubstart)'rpatoms'
        else
            print fsubstart, 'rpatoms'
        end if
      end if
      thespecies=adjustl(thespecies)

      ! read the coordinate file
      call read_coords(infile,informat,atoms,natoms,species,
     &                       nspecies,vecs)
      
      ! get the total number of atoms of
      ! the species of which you want to remove atoms
      do ispecies=1,nspecies
        if(species(ispecies)%name(1:2).eq.thespecies(1:2)) 
     &     nthere=species(ispecies)%howmany
      end do
      if(talk) then
        if(isopen12) then
!            write(12,'(8x,"There are ",I,1x,A2," atoms altogether.")')
!     &            nthere,thespecies(1:2)       
            WRITE(FMT2,*) nthere
            WRITE(FMT1,*) '(8x,"There are ",I',len(trim(FMT2)),',A2," at
     &oms altogether.")'
            write(12,FMT1) nthere,thespecies(1:2) 
!            write(12,'(8x,I," of them will be replaced by ",A2,".")')
!     &            nrp,trim(adjustl(theotherspecies))       
            if(nrp.ge.0) then 
              WRITE(FMT2,*) nrp
              WRITE(FMT1,*) '(8x,I',len(trim(FMT2)),',"of them will be r
     &eplaced by ",A2,".")'
              write(12,FMT1) nrp,trim(adjustl(theotherspecies))  
            end if
            if(nrp.lt.0) then
              WRITE(FMT2,*) -nrp
              WRITE(FMT1,*) '(8x,A2, " number ",I',len(trim(FMT2)),
     &        '," will be replaced by ",A2,".")'
              write(12,FMT1) thespecies(1:2),-nrp,
     &        trim(adjustl(theotherspecies))  
            end if
        else
!            print '(8x,"There are ",I,1x,A2," atoms altogether.")',
!     &            nthere,thespecies(1:2)       
            WRITE(FMT2,*) nthere
            WRITE(FMT1,*) '(8x,"There are ",I',len(trim(FMT2)),',A2," at
     &oms altogether.")'
            print FMT1, nthere,thespecies(1:2) 
!            print '(8x,I," of them will be replaced by",A2,".")',nrp,
!     &           trim(adjustl(theotherspecies))       
            !WRITE(FMT2,*) nrp
            !WRITE(FMT1,*) '(8x,I',len(trim(FMT2)),',"of them will be rep
!     &laced by ",A2,".")'
            !print FMT1, nrp,trim(adjustl(theotherspecies)) 
            if(nrp.ge.0) then 
              WRITE(FMT2,*) nrp
              WRITE(FMT1,*) '(8x,I',len(trim(FMT2)),',"of them will be r
     &eplaced by ",A2,".")'
              print FMT1, nrp,trim(adjustl(theotherspecies))  
            end if
            if(nrp.lt.0) then
              WRITE(FMT2,*) -nrp
              WRITE(FMT1,*) '(8x,A2, " number ",I',len(trim(FMT2)),
     &        '," will be replaced by ",A2,".")'
              print FMT1, thespecies(1:2),-nrp,
     &        trim(adjustl(theotherspecies))  
            end if
        end if
      end if

!      ! if the species is new, increase the number of elements you have
!      ! in your system
!      nnewspecies=nspecies
!      if(.not.any(species%name(1:8).eq.theotherspecies(1:8))) then
!            nnewspecies=nspecies+1
!      end if
!      allocate(newspecies(nnewspecies))
!      do ispecies=1,nspecies
!            newspecies(ispecies)=species(ispecies)
!      end do
!      if (nnewspecies.gt.nspecies) then
!          newspecies(nnewspecies)%name=theotherspecies
!          newspecies(nnewspecies)%howmany=nrp
!          newspecies(nnewspecies)%mass=0.0D0
!          newspecies(nnewspecies)%charge=0.0D0
!      end if

      ! if random replacement:
      if(nrp.ge.0) then
        ! call subroutine that draws balls from an urn
        call drawballs(nthere,nrp,replaced)

        ! get new structure
        !allocate(newatoms(natoms-nrp))
        !inewatom=0
        ithespecies=0
        do iatom=1,natoms
          !if (atoms(iatom)%name(1:2).ne.thespecies(1:2)) then
          !      inewatom=inewatom+1
          !      newatoms(inewatom)=atoms(iatom)
          !end if
          if (atoms(iatom)%name(1:2).eq.thespecies(1:2)) then
                ithespecies=ithespecies+1
                !if(.not.any(replaced.eq.ithespecies)) then
                if(any(replaced.eq.ithespecies)) then
                  !inewatom=inewatom+1
                  !newatoms(inewatom)=atoms(iatom)
                  atoms(iatom)%name(1:8)
     &             =trim(adjustl(theotherspecies(1:8)))
                end if
          end if
        end do
      end if
      ! if a special atom i sto be replaced:
      if(nrp.lt.0) then
        ithespecies=0
        do iatom=1,natoms
          if (atoms(iatom)%name(1:2).eq.thespecies(1:2)) then
                ithespecies=ithespecies+1
                !if(any(replaced.eq.ithespecies)) then
                if (ithespecies.eq.-nrp) then
                  atoms(iatom)%name(1:2)
     &             =trim(adjustl(theotherspecies(1:2)))
                end if
          end if
        end do
      end if
      
      ! get new numbers of species
      call getspecies(atoms,newspecies)
      nnewspecies=size(newspecies)
      call write_coords(outputfile,informat,atoms,natoms,newspecies,
     &                       nnewspecies,vecs)

      ! write a file with the atoms that were replaced
      open(51,file="REPLACED.DAT",status="replace")
      write(51,'("# The ",A," atoms that were replaced are:")')
     &  trim(thespecies)
      do n1=1,nrp
        write(51,*) replaced(n1)
      end do
      close(51)

      ! write goodbye
      if(talk) then
        if(isopen12) then
            write(12,fsubendext)'rpatoms'
        else
            print fsubendext, 'rpatoms'
        end if
      end if
      end subroutine

c---------------------------------------------------------------------

      subroutine rpatoms2(atoms,thespecies,nrp,
     &           theotherspecies)
      ! replaces nrp atoms of species "thespecies" by "theotherspecies"
      ! in array "atoms". result: atoms      

      use defs
      use misc
      implicit none
      character(len=*) thespecies,theotherspecies
      integer nrp
      ! local
      type(atom),allocatable :: atoms(:)!,newatoms(:)
      type(element),allocatable :: species(:),newspecies(:)
      integer natoms,nspecies,nnewspecies,iatom,ispecies,nthere
      integer ithespecies,n,n1 !,inewatom
      integer, allocatable :: replaced(:)
      double precision vecs(1:3,1:3)
      logical isopen12,isdrawn
      !character FMT1*1024,FMT2*1024

      talk=.true.
      ! write hello
      inquire(unit=12,opened=isopen12)
      if(talk) then
        if(isopen12) then
            write(12,fsubstart)'rpatoms2'
        else
            print fsubstart, 'rpatoms2'
        end if
      end if
      thespecies=adjustl(thespecies)

      natoms=size(atoms)  
      call getspecies(atoms,species)      
      nspecies=size(species)

      ! get the total number of atoms of
      ! the species of which you want to remove atoms
      do ispecies=1,nspecies
        if(species(ispecies)%name(1:2).eq.thespecies(1:2)) 
     &     nthere=species(ispecies)%howmany
      end do
      if(talk) then
        if(isopen12) then
            WRITE(FMT2,*) nthere
            WRITE(FMT1,*) '(8x,"There are ",I',len(trim(FMT2)),',1x,A2,
     &    " atoms altogether.")'
            write(12,FMT1) nthere,thespecies(1:2) 
            if(nrp.ge.0) then 
              WRITE(FMT2,*) nrp
              WRITE(FMT1,*) '(8x,I',len(trim(FMT2)),',"of them will be r
     &eplaced by ",A2,".")'
              write(12,FMT1) nrp,trim(adjustl(theotherspecies))  
            end if
            if(nrp.lt.0) then
              WRITE(FMT2,*) -nrp
              WRITE(FMT1,*) '(8x,A2, " number ",I',len(trim(FMT2)),
     &        '," will be replaced by ",A2,".")'
              write(12,FMT1) thespecies(1:2),-nrp,
     &        trim(adjustl(theotherspecies))  
            end if
        else
            WRITE(FMT2,*) nthere
            WRITE(FMT1,*) '(8x,"There are ",I',len(trim(FMT2)),',1x,A2,
     &" atoms altogether.")'
            print FMT1, nthere,thespecies(1:2) 
            if(nrp.ge.0) then 
              WRITE(FMT2,*) nrp
              WRITE(FMT1,*) '(8x,I',len(trim(FMT2)),'," of them will be 
     &replaced by ",A2,".")'
              print FMT1, nrp,trim(adjustl(theotherspecies))  
            end if
            if(nrp.lt.0) then
              WRITE(FMT2,*) -nrp
              WRITE(FMT1,*) '(8x,A2, " number ",I',len(trim(FMT2)),
     &        '," will be replaced by ",A2,".")'
              print FMT1, thespecies(1:2),-nrp,
     &        trim(adjustl(theotherspecies))  
            end if
        end if
      end if

      ! if random replacement:
      if(nrp.ge.0) then
        ! call subroutine that draws balls from an urn
        call drawballs(nthere,nrp,replaced)

        ! get new structure
        !allocate(newatoms(natoms-nrp))
        !inewatom=0
        ithespecies=0
        do iatom=1,natoms
          !if (atoms(iatom)%name(1:2).ne.thespecies(1:2)) then
          !      inewatom=inewatom+1
          !      newatoms(inewatom)=atoms(iatom)
          !end if
          if (atoms(iatom)%name(1:2).eq.thespecies(1:2)) then
                ithespecies=ithespecies+1
                !if(.not.any(replaced.eq.ithespecies)) then
                if(any(replaced.eq.ithespecies)) then
                  !inewatom=inewatom+1
                  !newatoms(inewatom)=atoms(iatom)
                  atoms(iatom)%name(1:2)
     &             =theotherspecies(1:2)
                end if
          end if
        end do
      end if
      ! if a special atom i sto be replaced:
      if(nrp.lt.0) then
        ithespecies=0
        do iatom=1,natoms
          if (atoms(iatom)%name(1:2).eq.thespecies(1:2)) then
                ithespecies=ithespecies+1
                !if(any(replaced.eq.ithespecies)) then
                if (ithespecies.eq.-nrp) then
                  atoms(iatom)%name(1:2)
     &             =trim(adjustl(theotherspecies(1:2)))
                end if
          end if
        end do
      end if

      ! write goodbye
      if(talk) then
        if(isopen12) then
            write(12,fsubendext)'rpatoms2'
        else
            print fsubendext, 'rpatoms2'
        end if
      end if

      end subroutine

c---------------------------------------------------------------------

      subroutine rmeqv(infile,informat,symmformat,outfile,tol)

      use defs
      use readcoords
      use misc
      use writecoords
      implicit none
      character(len=*) infile,informat,symmformat,outfile
      ! local
      type(atom),allocatable :: atoms(:),irredatoms(:),redatoms(:)
      type(element),allocatable :: species(:),irredspecies(:),
     &                             redspecies(:) 
      type(symop), allocatable :: symops(:)
      integer natoms,nspecies,iatom,ispecies
      integer nsym,nred,nirred,ired,iirred
      integer n,k,i,idum
      double precision vecs(1:3,1:3),xtrans(1:3),distvec(1:3),tol
      character redoutfile*40,line*256
      logical isopen12

      talk=.true.
      ! write hello
      inquire(unit=12,opened=isopen12)
      if(talk) then
        if(isopen12) then
            write(12,fsubstart)'rmeqv'
            write(12,'(8x,"using as tolerance ",F10.6)'),tol
        else
            print fsubstart, 'rmeqv'
            print '(8x,"using as tolerance ",F10.6)',tol
        end if
      end if

      ! read the coordinate file
      call read_coords(infile,informat,atoms,natoms,species,
     &                       nspecies,vecs)
      
c     read number of symmetry operations to be read in,
c     and symmetry operations 
      open(51,file='SYMM.DAT',status='old')
      read(51,*) nsym
      allocate(symops(nsym))
      select case(symmformat)
      case('gulp','GULP')
c     (in gulp format, with 'symmetry_operator'
c     before each operation)
            do n=1,nsym
                  read(51,'(A40)')line
                  read(51,*) symops(n)%mat(1,1:3),symops(n)%trala(1)
                  read(51,*) symops(n)%mat(2,1:3),symops(n)%trala(2)
                  read(51,*) symops(n)%mat(3,1:3),symops(n)%trala(3)
            end do
      case('mbpp','MBPP')
c     (in MBPP format)
            do n=1,nsym
                  read(51,*) idum,symops(n)%mat(1,1:3),
     &              symops(n)%mat(2,1:3),symops(n)%mat(3,1:3),
     &              symops(n)%trala(1:3)
            end do
      case default
            if(isopen12) then
                  write(12,ferrmssg) "unknown format of symmetry operati
     &ons"
            else
                  print ferrmssg, "unknown format of symmetry operations
     &"
            end if
            close(51)
            nerr=nerr+1
            return
      end select
      close(51)

      ! find irreducible atoms
      ! at the beginning assume that every atom is irreducible:
      nred=0
      nirred=natoms
      atoms%reducible=.false.
c     for each atom check if it is reducible
      !xall(n,1:3)=x(1:3)
      !reducible(n)=.False.
      !nstart=1
c      if(n.gt.n_atom/2) nstart=n_atom/2+1
      do iatom=1,natoms  ! for all atoms
         do k=1,nsym  ! for all symmetry operations apply the symmetry operation to
             ! the atom coordinates 
             xtrans=matmul(symops(k)%mat,atoms(iatom)%where)
     &              +symops(k)%trala 
!            do i=1,3
!                  xtrans(i)=xall(iatom,1)*matrix(k,i,1)!+trala(k,i)
!     &                     +xall(iatom,2)*matrix(k,i,2)!+trala(k,i)
!     &                     +xall(iatom,3)*matrix(k,i,3)+trala(k,i)     
                  !do while (xtrans(i).gt.0.5)
                  !      xtrans(i)=xtrans(i)-1.0
                  !end do
                  !do while (xtrans(i).le.-0.5)
                  !      xtrans(i)=xtrans(i)+1.0
                  !end do
!            end do
             ! look if another atom is equivalent by symmetry
             do n=1,iatom-1 
!                  if (abs(xall(m,1)-xtrans(1)).lt.tol.and.
!     &                abs(xall(m,2)-xtrans(2)).lt.tol.and.
!     &                abs(xall(m,3)-xtrans(3)).lt.tol.and.
!     &                nameat(n)==nameat(m).and.cs(n)==cs(m)) 
!     &                  reducible(n)=.True.
                ! If the atoms belong to the same species,
                ! get the distance between the transposed original atom and
                ! the other atom
                if (atoms(iatom)%name.eq.atoms(n)%name.and.
     &           atoms(iatom)%core.eq.atoms(n)%core) then
                  distvec=xtrans-atoms(n)%where
                  ! use the distance modulo the lattice constants
                  do i=1,3
                     do while (distvec(i).gt.0.5D0)
                        distvec(i)=distvec(i)-1.0D0
                     end do
                     do while (distvec(i).le.-0.5D0)
                        distvec(i)=distvec(i)+1.0D0
                     end do
                  end do
                  ! if the distance is below a threshold, the atom
                  ! coordinates are equal
                  if (absvec(distvec).lt.tol) then
                      atoms(iatom)%reducible=.true.
                      !nred=nred+1
                  end if
                end if  
             end do
         end do
         if (atoms(iatom)%reducible) nred=nred+1
      end do
      if(talk) print '(8x,"number of red. atoms:",I5)',nred

      ! get reducible and irreducible atoms
      nirred=natoms-nred
      
      !print '(8x,"number of red./irred. atoms:",I,I)', nred,nirred
      !do i=1,natoms
      !  print*,atoms(i)%name,atoms(i)%reducible
      !end do


      allocate(redatoms(nred),irredatoms(nirred))
      ired=0
      iirred=0
      do iatom=1,natoms
        if (atoms(iatom)%reducible.eqv..true.) then
              ired=ired+1
              redatoms(ired)=atoms(iatom)
        else
              iirred=iirred+1
              irredatoms(iirred)=atoms(iatom)
        end if
      end do
      
      !do i=1,nred
      !  print*,redatoms(i)%name,redatoms(i)%reducible
      !end do
      !do i=1,nirred
      !  print*,irredatoms(i)%name,irredatoms(i)%reducible
      !end do

      ! get number of reducible and irreducible atoms of each species
!      ispecies=1
!      do iatom=1,nred
!        newspecies=.true.
!        do i=1,iatom-1
!          if (redatoms(i)%name(1:2).eq.redatoms(iatom)%name(1:2))
!     &        newspecies=.false.
!        end do
!        if(newspecies) then
!            species(ispecies)%name=redatoms(iatom)%name
!            ispecies=ispecies+1
!        end if
!      end do  

      ! get numbers of reducible atoms of each species
      allocate(redspecies(nspecies),irredspecies(nspecies))
      redspecies=species
      irredspecies=species
      do ispecies=1,nspecies
         redspecies(ispecies)%howmany=0
         do i=1,nred
            if(redatoms(i)%name(1:2).eq.redspecies(ispecies)%name(1:2)) 
     &         redspecies(ispecies)%howmany=
     &         redspecies(ispecies)%howmany+1
         end do
      end do
      ! get numbers of irreducible atoms of each species
      do ispecies=1,nspecies
         irredspecies(ispecies)%howmany=0
         do i=1,nirred
            if(irredatoms(i)%name(1:2).eq.
     &         irredspecies(ispecies)%name(1:2)) 
     &         irredspecies(ispecies)%howmany=
     &         irredspecies(ispecies)%howmany+1
         end do
      end do

!      ! get new structure
!      allocate(newatoms(natoms-nrm))
!      inewatom=0
!      ithespecies=0
!      do iatom=1,natoms
!        if (atoms(iatom)%name(1:2).ne.thespecies(1:2)) then
!              inewatom=inewatom+1
!              newatoms(inewatom)=atoms(iatom)
!        end if
!        if (atoms(iatom)%name(1:2).eq.thespecies(1:2)) then
!              ithespecies=ithespecies+1
!              if(.not.any(vacs.eq.ithespecies)) then
!                inewatom=inewatom+1
!                newatoms(inewatom)=atoms(iatom)
!              end if
!        end if
!      end do
      
      ! write a file with the irreducible atoms
      call write_coords(outfile,informat,irredatoms,nirred,
     &     irredspecies,nspecies,vecs)

      ! write a file with the reducible atoms
      write(redoutfile,*)"REDUCIBLE.",trim(adjustl(informat))
      redoutfile=adjustl(redoutfile)
      call write_coords(redoutfile,informat,redatoms,
     &     nred,redspecies,nspecies,vecs)

      ! write a file with the symmetry operations in other formats
      open(51,file="SYMM.OUT",status="replace")
      do k=1,nsym
        write(51,'("symmetry_operator")')
        do i=1,3
          write(51,'(4(F10.6))')
     &       symops(k)%mat(i,1:3),symops(k)%trala(i)
        end do
      end do   
      close(51)

      ! write goodbye
      if(talk) then
        if(isopen12) then
            write(12,fsubendext)'rmeqv'
        else
            print fsubendext, 'rmeqv'
        end if
      end if
      end subroutine rmeqv

c---------------------------------------------------------------------

      subroutine addeqv(infile,informat,symmformat,outfile,tol)

      use defs
      use readcoords
      use misc
      use writecoords
      implicit none
      character(len=*) infile,informat,symmformat,outfile
      ! local
      type(atom),allocatable :: atoms(:),irredatoms(:),redatoms(:)
      type(atom),allocatable :: addedatoms(:),dummyatoms(:),allatoms(:)
      type(element),allocatable :: species(:),irredspecies(:),
     &      redspecies(:),addedspecies(:),allspecies(:)
      type(symop), allocatable :: symops(:)
      integer natoms,nspecies,iatom,ispecies,naddedspecies,nirredspecies
      integer nredspecies,nallspecies
      integer nsym,nred,nirred,ired,iirred,iadded,nadded,nall
      integer n,k,i,idum
      double precision vecs(1:3,1:3),xtrans(1:3),distvec(1:3),tol
      character redoutfile*256,line*256
      logical isopen12,add

      talk=.true.
      ! write hello
      inquire(unit=12,opened=isopen12)
      if(talk) then
        if(isopen12) then
            write(12,fsubstart)'addeqv'
            write(12,'(8x,"using as tolerance ",F10.6)'),tol
        else
            print fsubstart, 'addeqv'
            print '(8x,"using as tolerance ",F10.6)',tol
        end if
      end if

      ! read the coordinate file
      call read_coords(infile,informat,atoms,natoms,species,
     &                       nspecies,vecs)
      
c     read number of symmetry operations to be read in,
c     and symmetry operations 
      open(51,file='SYMM.DAT',status='old')
      read(51,*) nsym
      allocate(symops(nsym))
      select case(symmformat)
      case('gulp','GULP')
c     (in gulp format, with 'symmetry_operator'
c     before each operation)
            do n=1,nsym
                  read(51,'(A40)')line
                  read(51,*) symops(n)%mat(1,1:3),symops(n)%trala(1)
                  read(51,*) symops(n)%mat(2,1:3),symops(n)%trala(2)
                  read(51,*) symops(n)%mat(3,1:3),symops(n)%trala(3)
            end do
      case('mbpp','MBPP')
c     (in MBPP format)
            do n=1,nsym
                  read(51,*) idum,symops(n)%mat(1,1:3),
     &              symops(n)%mat(2,1:3),symops(n)%mat(3,1:3),
     &              symops(n)%trala(1:3)
            end do
      case default
            if(isopen12) then
                  write(12,ferrmssg) "unknown format of symmetry operati
     &ons"
            else
                  print ferrmssg, "unknown format of symmetry operations
     &"
            end if
            close(51)
            nerr=nerr+1
            return
      end select
      close(51)

      ! find irreducible atoms
      ! at the beginning assume that every atom is irreducible:
      nred=0
      nirred=natoms
      nadded=0
      atoms%reducible=.false.
c     for each atom check if it is reducible
      do iatom=1,natoms  ! for all atoms
         do k=1,nsym  ! for all symmetry operations apply the symmetry operation to
             ! the atom coordinates 
             xtrans=matmul(symops(k)%mat,atoms(iatom)%where)
     &              +symops(k)%trala 
             do i=1,3
                do while (xtrans(i).ge.0.750D0) 
                  xtrans(i)=xtrans(i)-1.0d0
                end do
                do while (xtrans(i).le.-0.250D0) 
                  xtrans(i)=xtrans(i)+1.0d0
                end do
             end do !i
             ! look if another atom is equivalent by symmetry
             do n=1,iatom-1 
                ! If the atoms belong to the same species,
                ! get the distance between the transposed original atom and
                ! the other atom
                if (atoms(iatom)%name.eq.atoms(n)%name.and.
     &           atoms(iatom)%core.eq.atoms(n)%core) then
                  distvec=xtrans-atoms(n)%where
                  ! use the distance modulo the lattice constants
                  do i=1,3
                     do while (distvec(i).gt.0.50D0)
                        distvec(i)=distvec(i)-1.0D0
                     end do
                     do while (distvec(i).le.-0.50D0)
                        distvec(i)=distvec(i)+1.0D0
                     end do
                  end do
                  ! if the distance is below a threshold, the atom
                  ! coordinates are equal
                  if (absvec(distvec).lt.tol) then
                      atoms(iatom)%reducible=.true.
                  end if
                end if  
             end do
         end do
         if (atoms(iatom)%reducible) nred=nred+1
      end do
      if(talk) print'(8x,"number of red. atoms:",I5)',nred

      ! get reducible and irreducible atoms
      nirred=natoms-nred
      if(talk) print'(8x,"number of irred. atoms:",I5)',nirred
      
      allocate(redatoms(nred),irredatoms(nirred))
      ired=0
      iirred=0
      do iatom=1,natoms
        if (atoms(iatom)%reducible.eqv..true.) then
              ired=ired+1
              redatoms(ired)=atoms(iatom)
        else
              iirred=iirred+1
              irredatoms(iirred)=atoms(iatom)
        end if
      end do

      ! get numbers of reducible atoms of each species
      call getspecies(redatoms,redspecies)
      nredspecies=size(redspecies)

      ! get numbers of irreducible atoms of each species
      call getspecies(irredatoms,irredspecies)
      nirredspecies=size(irredspecies)

      !
      ! begin add atoms to irreducible ones according to symmetry 
      !
      allocate(addedatoms(nirred*nsym)) ! allocate the maximum possible
      ! number of new atoms to add. Excess space will be removed later. 
      iadded=0
      do iirred=1,nirred
        do k=1,nsym  ! for all symmetry operations apply the symmetry operation to the atom coordinates 
           xtrans=matmul(symops(k)%mat,irredatoms(iirred)%where)
     &              +symops(k)%trala 
           do i=1,3
              do while (xtrans(i).ge.0.750D0) 
                xtrans(i)=xtrans(i)-1.0d0
              end do
              do while (xtrans(i).le.-0.250D0) 
                xtrans(i)=xtrans(i)+1.0d0
              end do
           end do !i
           ! look if another atom with the same transformed atom coordinates exist
           ! already amog the irreducible atoms. If yes, do not add an atom 
           add=.true.
           do n=1,nirred
              ! If the atoms belong to the same species,
              ! get the distance between the transposed original atom and
              ! the other atom
              !if (n.ne.iirred.and.irredatoms(iirred)%name.eq.             &
              if (irredatoms(iirred)%name.eq.                             &
     &         irredatoms(n)%name.and.irredatoms(iirred)%core.eq.         &
     &         irredatoms(n)%core) then
                distvec=xtrans-irredatoms(n)%where
                ! use the distance modulo the lattice constants
                do i=1,3
                   do while (distvec(i).gt.0.5D0)
                      distvec(i)=distvec(i)-1.0D0
                   end do
                   do while (distvec(i).le.-0.5D0)
                      distvec(i)=distvec(i)+1.0D0
                   end do
                end do
                ! if the distance is above a threshold, add a new atom
                if (absvec(distvec).lt.tol) then
                  add=.false.
                end if ! (absvec(distvec).ge.tol)
              end if ! same atom name
           end do ! n
           ! look if another atom with the same transformed atom coordinates exist
           ! already among the added atoms. If yes, do not add an atom 
           do n=1,iadded
              ! If the atoms belong to the same species,
              ! get the distance between the transposed original atom and
              ! the other atom
              !if (n.ne.iirred.and.irredatoms(iirred)%name.eq.             &
              if (irredatoms(iirred)%name.eq.                             &
     &         addedatoms(n)%name.and.irredatoms(iirred)%core.eq.         &
     &         addedatoms(n)%core) then
                distvec=xtrans-addedatoms(n)%where
                ! use the distance modulo the lattice constants
                do i=1,3
                   do while (distvec(i).gt.0.5D0)
                      distvec(i)=distvec(i)-1.0D0
                   end do
                   do while (distvec(i).le.-0.5D0)
                      distvec(i)=distvec(i)+1.0D0
                   end do
                end do
                ! if the distance is above a threshold, add a new atom
                if (absvec(distvec).lt.tol) then
                  add=.false.
                end if ! (absvec(distvec).ge.tol)
              end if ! same atom name
           end do ! n
           if(add) then
             iadded=iadded+1
             addedatoms(iadded)=atoms(iirred)
             addedatoms(iadded)%where=xtrans
             call frac2abs(addedatoms(iadded)%where,vecs,                 &
     &              addedatoms(iadded)%abswhere)
           end if
        end do ! k
      end do ! iirred

      ! remove excess space in the added atom array. 
      nadded=iadded
      if (nadded.gt.0) then
        allocate(dummyatoms(nadded))
        dummyatoms(1:nadded)=addedatoms(1:nadded)
        deallocate(addedatoms)
        allocate(addedatoms(nadded))
        addedatoms=dummyatoms
        deallocate(dummyatoms) 
        call getspecies(addedatoms,addedspecies)
        naddedspecies=size(addedspecies)
      end if ! (iadded.gt.0)
      print '(8x,"number of new atoms:",I5)',nadded
      nall=nirred+nadded
      allocate(allatoms(nall))
      allatoms(1:nirred)=irredatoms(1:nirred)
      if(nadded.gt.0) allatoms(nirred+1:nall)=addedatoms(1:nadded)
      call getspecies(allatoms,allspecies)
      nallspecies=size(allspecies)
      !
      ! end add atoms to irreducible ones according to symmetry 
      !

      ! write a file with the added atoms
      if(nadded.gt.0) then
        redoutfile=""
        write(redoutfile,*)"ADDED.",trim(adjustl(informat))
        redoutfile=adjustl(redoutfile)
        call write_coords(redoutfile,informat,addedatoms,nadded,
     &     addedspecies,naddedspecies,vecs)
      end if
      
      ! write a file with all (the irreducible and the added) atoms
      call write_coords(outfile,informat,allatoms,nall,                   &
     &     allspecies,nallspecies,vecs)

      ! write a file with the reducible atoms
      redoutfile=""
      write(redoutfile,*)"REDUCIBLE.",trim(adjustl(informat))
      redoutfile=adjustl(redoutfile)
      open(51,file=redoutfile,status="replace")
      close(51)
      if (nred.gt.0)then
        call write_coords(redoutfile,informat,redatoms,                   &
     &     nred,redspecies,nspecies,vecs)
      end if
      
      ! write a file with the irreducible atoms
      write(redoutfile,*)"IRREDUCIBLE.",trim(adjustl(informat))
      redoutfile=adjustl(redoutfile)
      call write_coords(redoutfile,informat,irredatoms,
     &     nirred,irredspecies,nspecies,vecs)


      ! write a file with the symmetry operations in other formats
      open(51,file="SYMM.OUT",status="replace")
      do k=1,nsym
        write(51,'("symmetry_operator")')
        do i=1,3
          write(51,'(4(F10.6))')
     &       symops(k)%mat(i,1:3),symops(k)%trala(i)
        end do
      end do   
      close(51)

      ! write goodbye
      if(talk) then
        if(isopen12) then
            write(12,fsubendext)'addeqv'
        else
            print fsubendext, 'addeqv'
        end if
      end if
      end subroutine addeqv

c--------------------------------------------------------------------

      subroutine addatoms(atoms1,atoms2,atoms3)
      ! combines atoms1 and atoms2 to atoms3
      
      use defs
      !use misc
      implicit none

      ! variables
      type(atom) :: atoms1(:),atoms2(:)
      type(atom), allocatable :: atoms3(:)

      ! internal variables
      logical isopen12 !,newspecies
      integer natoms1,natoms2,natoms3
       
      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(talk) then
        if(isopen12) then
            write(12,fsubstart) trim('addatoms')
        else
            print fsubstart, trim(adjustl('addatoms'))
        end if
      end if

      natoms1=size(atoms1)
      natoms2=size(atoms2)
      natoms3=natoms1+natoms2
      if(allocated(atoms3)) deallocate(atoms3)
      allocate(atoms3(1:natoms3))
      atoms3(1:natoms1)=atoms1(1:natoms1)
      atoms3(natoms1+1:natoms3)=atoms2(1:natoms2)
                  
      if(talk) then
        if(isopen12) then
            write(12,fsubendext) 'addatoms'
        else
            print fsubendext, 'addatoms'
        end if
      end if
      return
      
      end subroutine 

c--------------------------------------------------------------------

      subroutine keepoccatoms(atoms1,atoms2)
      ! keeps only atoms with an occupation greater than a threshold
      ! (0.9)      

      use defs
      !use misc
      implicit none

      ! variables
      type(atom) :: atoms1(:) ! original atoms
      type(atom), allocatable :: atoms2(:)  ! final atoms

      ! internal variables
      logical isopen12 
      integer natoms1 ! original number of atoms
      integer nocc ! number of occupied sites 
      integer iatom,iatom2
       
      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim('keepoccatoms')
      else
            print fsubstart, trim(adjustl('keepoccatoms'))
      end if

      ! count occupied sites
      nocc=0
      natoms1=size(atoms1)
      do iatom=1,natoms1
        if(atoms1(iatom)%occ.gt.0.9d0) nocc=nocc+1
      end do
      ! if no atom is occupied, return an error
      if(nocc.eq.0) goto 1000

      ! store the occupied atoms from atoms1 in atoms2
      allocate(atoms2(1:nocc))
      iatom2=1
      do iatom=1,natoms1
        if(atoms1(iatom)%occ.gt.0.9d0) then
          atoms2(iatom2)=atoms1(iatom)
          iatom2=iatom2+1
        end if
      end do
      
      ! Normal exit      
      if(isopen12) then
            write(12,fsubendext) 'keepoccatoms'
      else
            print fsubendext, 'keepoccatoms'
      end if
      return
   
      ! Return with errors      
1000  nerr=nerr+1
      print ferrmssg,"No occupied site found"
      return

      end subroutine 

c---------------------------------------------------------------------

      subroutine diffatoms(atoms1,atoms2,vecs,atomssame,atomsdiff,tol)
      ! compares atoms1 and atoms2. Atoms that are (not) contained in both atoms1 and 
      ! atoms2 are written to atoms3 (atoms4)
      use misc
      use defs
      implicit none
      type(atom) atoms1(:),atoms2(:)
      type(atom), allocatable :: atomssame(:),atomsdiff(:)  
      double precision tol   ! tolerance for distance between atoms (if distance<tol, the positions are the same)  
      double precision vecs(3,3)
      ! internal variables:
      integer natoms1,natoms2,nsame,ndiff,iatom,iatom2,iatom3,iatom4  
      integer natoms
      double precision, allocatable :: coorddiffs(:,:) ! coordinate differences
      logical diff   

      if(talk) print fsubstart,"diffatoms"
      ! initialize
      natoms1=size(atoms1) 
      natoms2=size(atoms2)
      ! count atoms that are the same (different) in the two structures 
      nsame=0
      ndiff=0
      do iatom=1,natoms1
        do iatom2=1,natoms2
          if (atoms1(iatom)%name(1:2).eq.atoms2(iatom2)%name(1:2)
     &        .and.absvec(atoms1(iatom)%abswhere
     &        -atoms2(iatom2)%abswhere).lt.tol) nsame=nsame+1
        end do  
      end do
      ndiff=(natoms1-nsame)+(natoms2-nsame)
      print '(8x,I0," equal atoms")',nsame
      print '(8x,I0," different atoms")',ndiff
      allocate(atomssame(1:nsame),atomsdiff(1:ndiff))
      ! get the arrays of common and different atoms
      iatom3=0 ! counter for common atoms
      iatom4=0 ! counter for different atoms
      do iatom=1,natoms1
        diff=.true.
        do iatom2=1,natoms2
          if (atoms1(iatom)%name(1:2).eq.atoms2(iatom2)%name(1:2)
     &        .and.absvec(atoms1(iatom)%abswhere
     &        -atoms2(iatom2)%abswhere).lt.tol) then
            iatom3=iatom3+1
            atomssame(iatom3)=atoms1(iatom)
            diff=.false.
          end if 
        end do
        if(diff) then
          iatom4=iatom4+1
          atomsdiff(iatom4)=atoms1(iatom)
        end if
      end do  
      do iatom2=1,natoms2
        diff=.true.
        do iatom=1,natoms1
          if (atoms1(iatom)%name(1:2).eq.atoms2(iatom2)%name(1:2)
     &        .and.absvec(atoms1(iatom)%abswhere
     &        -atoms2(iatom2)%abswhere).lt.tol) then
            diff=.false.
          end if 
        end do
        if(diff) then
          iatom4=iatom4+1
          atomsdiff(iatom4)=atoms2(iatom2)
        end if
      end do  
      ! 
      ! begin get coordinate differences
      !
      natoms=max(natoms1,natoms2)
      allocate(coorddiffs(natoms,1:4))
      coorddiffs=100.0d0
      do iatom=1,natoms1
        do iatom2=1,natoms2
          if (atoms1(iatom)%name(1:2).eq.atoms2(iatom2)%name(1:2)
     &        .and.norm2(atoms1(iatom)%abswhere
     &        -atoms2(iatom2)%abswhere).lt.coorddiffs(iatom,4)) then
            coorddiffs(iatom,1:3)=atoms1(iatom)%abswhere                  &
     &                      -atoms2(iatom2)%abswhere
            coorddiffs(iatom,4)=norm2(atoms1(iatom)%abswhere              &
     &                      -atoms2(iatom2)%abswhere)
          end if
        end do  
      end do
      open(51,file="COORDDIFFS.OUT")
      open(52,file="COORDDIFFSFRAC.OUT")
      do iatom=1,natoms1
        write(51,'(4(F20.16))')  coorddiffs(iatom,1:4)
        ! fractional coordinates:
        call abs2frac(coorddiffs(iatom,1:3),vecs,coorddiffs(iatom,1:3))
        coorddiffs(iatom,4)=norm2(coorddiffs(iatom,1:3))
        write(52,'(4(F20.16))')  coorddiffs(iatom,1:4)
      end do
      close(51)
      close(52)
      deallocate(coorddiffs)
      !
      open(51,file="COORDDIFFS_ORIGINAL_ORDER.OUT")
      open(52,file="COORDDIFFSFRAC_ORIGINAL_ORDER.OUT")
      do iatom=1,min(natoms1,natoms2)
       write(51,'(4(F20.16))') atoms2(iatom)%abswhere                     &
     &   -atoms1(iatom)%abswhere,norm2(atoms2(iatom)%abswhere             &
     &   -atoms1(iatom)%abswhere)
       write(52,'(4(F20.16))') atoms2(iatom)%where-atoms1(iatom)%where,   &
     &  norm2(atoms2(iatom)%where-atoms1(iatom)%where)
      end do
      close(51)
      close(52)
      ! 
      ! end get coordinate differences
      !
      if(talk) print fsubendext,"diffatoms" 
      
      end subroutine diffatoms

c---------------------------------------------------------------------

      subroutine add_atoms(infile,informat,thespecies,nadd,
     &           coords,outputfile)

      use defs
      use readcoords
      use misc
      use writecoords
      implicit none
      character(len=*) infile,informat,thespecies,outputfile
      integer nadd
      ! local
      type(atom),allocatable :: atoms(:),newatoms(:)
      type(element),allocatable :: species(:),newspecies(:)
      integer natoms,nspecies,nnewspecies,iatom,ispecies,nthere
      integer inewatom,ithespecies,n,n1
      integer, allocatable :: added(:)
      double precision vecs(1:3,1:3)
      double precision coords(:,:)
      logical isopen12

      talk=.true.
      ! write hello
      inquire(unit=12,opened=isopen12)
      if(talk) then
        if(isopen12) then
            write(12,fsubstart)'addatoms'
        else
            print fsubstart, 'addatoms'
            print '(8x,"species to add: ", A2)', thespecies
            print '(8x,"number of atoms to add: ", I0)', nadd
            do n=1,nadd
              print '(8x,3(F10.6))',coords(n,1:3)
            end do
        end if
      end if
      thespecies=adjustl(thespecies)

      ! read the coordinate file
      call read_coords(infile,informat,atoms,natoms,species,
     &                       nspecies,vecs)
      
      ! get the total number of atoms of
      ! the species of which you want to add atoms
      nthere=0
      do ispecies=1,nspecies
        if(species(ispecies)%name(1:2).eq.thespecies(1:2)) 
     &     nthere=species(ispecies)%howmany
      end do
      if(talk) then
        if(isopen12) then
            WRITE(FMT2,*) nthere
            WRITE(FMT1,*) '(8x,"There are ",I',
     &       len(trim(FMT2)),',1x,A2," atoms already.")'
            write(12,FMT1)  nthere,thespecies(1:2) 
            WRITE(FMT2,*) nadd
            WRITE(FMT1,*) '(8x,I',len(trim(FMT2)),'" will be added.")'
            write(12,FMT1) nadd 
        else
            WRITE(FMT2,*) nthere
            WRITE(FMT1,*) '(8x,"There are ",I',
     &       len(trim(FMT2)),',1x,A2," atoms already.")'
            print FMT1, nthere,thespecies(1:2) 
            WRITE(FMT2,*) nadd
            WRITE(FMT1,*) '(8x,I',len(trim(FMT2)),'" will be added.")'
            print FMT1, nadd 
        end if
      end if

      ! 
      ! if positive number to be added
      !
      if (nadd.gt.0) then
        !
        ! get new structure
        allocate(newatoms(natoms+nadd))
        newatoms(1:natoms)=atoms(1:natoms)
        do n=1,nadd
          newatoms(natoms+n)%name(:)=" "
          newatoms(natoms+n)%core(:)="core"
          newatoms(natoms+n)%name(1:2)=thespecies(1:2)
          newatoms(natoms+n)%where=coords(n,1:3)
          call frac2abs(newatoms(natoms+n)%where,vecs,                    &
     &         newatoms(natoms+n)%abswhere)
        end do
        ! 
        ! end if positive number to be added
      else
        !
        allocate(newatoms(natoms))
        newatoms=atoms
        !
      end if
      !
      ! get new numbers of species
      call getspecies(newatoms,newspecies)
      nnewspecies=size(newspecies)
      call write_coords(outputfile,informat,newatoms,size(newatoms),
     &     newspecies,nnewspecies,vecs)
      !
      ! write goodbye
      if(talk) then
        if(isopen12) then
            write(12,fsubendext)'addatoms'
        else
            print fsubendext, 'addatoms'
        end if
      end if
      end subroutine add_atoms

c---------------------------------------------------------------------
     
      end module

