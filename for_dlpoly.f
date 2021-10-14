      module for_dlpoly
      implicit none

      contains

c---------------------------------------------------------------------

      subroutine sort_molecules()
      ! this subroutine reads coordinates in DL_POLY CONFIG file format
      ! (without the CONFIG file header, only the coordinates)
      ! and sorts them into molecules according to the order specified
      ! in "INPUT.COFIMA" after the keyword "sort molecules".

      use defs
      implicit none

      integer natoms,nmolecules,natmol

      !type atom
      !sequence
      !character name*8
      !logical written
      !double precision :: where(3)
      !end type atom

      type(atom), allocatable :: allatoms(:)
      type(atom), allocatable :: molecule(:)
      !---------------------------------------------------------------------
    
      character filename*20,outfile*20,line*100,filetype*20
      integer i,iallatoms, imol,iatom
      logical towrite
      integer atindex,imcon,levcfg
      double precision vecs(1:3,1:3)

      !---------------------------------------------------------------------
      !---------------------------------------------------------------------

      ! read coordinate file name, total number of atoms, number of molecules and number of
      ! atoms per molecule from 'INPUT.COFIMA' 
      ! and allocate 'allatoms', 'molecule'
      
      write(12,fsubstart)trim(adjustl('sort_molecules'))
      write(12,*)' '
      read(10,*) filename
      read(10,'(A20)') filetype
      read(10,*) outfile
      write(12,*) '     Input will be read from file: ',filename
      write(12,*) '     The input file is of type: ',filetype
      write(12,*) '     Output will be written to file: ',outfile
      write(12,*)'      Reading atom and molecule information...'
      open(19,file=outfile,status='replace')
      open(11,file=filename,status='old')
      select case(filetype)
      case('CONFIG')
            read(11,'(A100)') line
            write(19,'(A100)') line
            read(11,*) levcfg,imcon,natoms
            write(19,'(I10,I10,I16)') levcfg,imcon,natoms
            do i=1,3
              read(11,*) vecs(i,1:3)
              write(19,'(3(F20.10))') vecs(i,1:3)
            end do
      case('only coords')
      case default
            write(12,*) '       Unknown file format. Returning to main.'
            return
      end select
      read(10,*) natoms
      read(10,*) nmolecules
      read(10,*) natmol
      allocate(allatoms(natoms),molecule(natmol))
      write(12,*) '     ',natoms,' atoms, '
      write(12,*) '     ', nmolecules,' molecules, '
      write(12,*) '     ',natmol, ' atoms per molecule' 
      write(12,*) '      reading coords from ',filename
      write(12,*) ' '
      write(12,*) '      The molecule: '
      do i=1,natmol
        read(10,*) molecule(i)%name 
        write(12,*)'    ',molecule(i)%name 
      end do
      write(12,*) ' '
      
      !read coordinates of atoms 
      !write(12,*) '     Now come all atoms:'
      do i=1, natoms
        read(11,*) allatoms(i)%name
        read(11,*) allatoms(i)%where(1:3)
        !write(12,'(A8,3(F10.6))') allatoms(i)%name,
c     &   allatoms(i)%where(1:3)
        allatoms(i)%written=.false.
      end do
c      write(12,*) ' '
      close(11)

      ! divide 'allatoms' into molecules
      !write(12,*) 'Now come all molecules:'
      atindex=1
      do imol=1, nmolecules
        do iatom=1,natmol
          towrite=.false.
          do iallatoms=1,natoms
          if(allatoms(iallatoms)%name.eq.molecule(iatom)%name
     &       .and..not.towrite.and..not.allatoms(iallatoms)%written)then
!                write(12,'(A8,I10)') allatoms(iallatoms)%name,atindex
!                write(12,'(F16.9,2(F20.9))') 
                write(19,'(A8,I12)') allatoms(iallatoms)%name,atindex
                write(19,'(F16.9,2(F20.9))') 
     &           allatoms(iallatoms)%where(1:3)
            towrite=.true.
            allatoms(iallatoms)%written=.true.
            atindex=atindex+1
          end if
          end do
        end do
      end do

      write(12,fsubend)
      end subroutine

c---------------------------------------------------------------------

!      subroutine abs2frac()
!      ! this subroutine reads absolute coordinates and lattice vectors and
!      ! transforms the coordinates into fractional ones
!      
!      use misc 
!      use defs
!      implicit none
!
!      character infilename*20,
!     & infiletype*20,outfilename*20,line*100,name*8
!      integer levcfg,imcon,natms,i,j,iatom,atnum
!      double precision vecs(1:3,1:3),xabs(1:3),xfrac(1:3)
!      double precision matrix(1:3,1:3),inversematrix(1:3,1:3)
!
!      !---------------------------------------------------------------
!      write(12,fsubstart)trim(adjustl('abs2frac'))
!      read(10,'(A20)')infilename
!      read(10,'(A20)')infiletype
!      read(10,'(A20)')outfilename
!      write(12,*) '     Infilename: ',infilename
!      write(12,*) '     Infiletype: ',infiletype
!      write(12,*) '     Outfilename: ',outfilename
!      open(14,file=infilename,status='old')
!      open(15,file=outfilename,status='replace')
!      select case(infiletype)
!      case('CONFIG') 
!            read(14,'(A100)') line
!            write(12,'(A6,A100)') '      ',line
!            write(15,'(A100)') line
!            read(14,*) levcfg,imcon,natms
!            write(12,*)'      levcfg,imcon,natms: ',levcfg,imcon,natms 
!            write(15,'(I10,I10,I16)')levcfg,imcon,natms 
!!            write(12,*)'      lattice vecs:'
!            do i=1,3
!              read(14,*) vecs(i,1:3)
!!              write(12,'(3(F20.10))') vecs(i,1:3)
!              write(15,'(3(F20.10))') vecs(i,1:3)
!            end do
!            do i=1,3
!              do j=1,3
!                matrix(i,j)=vecs(j,i)
!              end do
!            end do
!            call invert_matrix(matrix,inversematrix)
!!            write(12,*)'      atomic positions: ' 
!            do iatom=1,natms
!              read(14,'(A8,I)') name,atnum
!!              write(12,'(A8,I)') name,atnum
!              write(15,'(A8,I)') name,atnum
!              read(14,*) xabs(1:3)
!!              write(12,'(3(F20.10))') xabs(1:3)
!              do i=1,3
!                xfrac(i)=0.0D0
!                do j=1,3
!                  xfrac(i)=xfrac(i)+inversematrix(i,j)*xabs(j)
!                end do
!              end do
!              write(15,'(3(F20.10))') xfrac(1:3)
!            end do
!            close(14)
!            close(15)
!      case default
!            print*,'        Error: Unknown filetype.'
!            write(12,*) '        Error: Unknown filetype.'
!            close(14)
!            close(15)
!            return
!      end select
!      
!      write(12,fsubend)
!      end subroutine
!
!c---------------------------------------------------------------------
!
!      subroutine frac2abs()
!      ! this subroutine reads fractional coordinates and lattice vectors and
!      ! transforms the coordinates into absolute ones
!      
!      use defs
!      implicit none
!
!      character infilename*20,
!     & infiletype*20,outfilename*20,line*100,name*8
!      integer levcfg,imcon,natms,i,j,iatom,atnum
!      double precision vecs(1:3,1:3),xabs(1:3),xfrac(1:3)
!      double precision matrix(1:3,1:3),inversematrix(1:3,1:3)
!
!      !---------------------------------------------------------------
!      write(12,fsubstart)trim(adjustl('frac2abs'))
!      write(12,*)'      '
!      read(10,'(A20)')infilename
!      read(10,'(A20)')infiletype
!      read(10,'(A20)')outfilename
!      write(12,*) '      Infilename: ',infilename
!      write(12,*) '      Infiletype: ',infiletype
!      write(12,*) '      Outfilename: ',outfilename
!      open(14,file=infilename,status='old')
!      open(15,file=outfilename,status='replace')
!      select case(infiletype)
!      case('CONFIG') 
!            read(14,'(A100)') line
!            write(12,'(A6,A100)') '      ',line
!            write(15,'(A100)') line
!            read(14,*) levcfg,imcon,natms
!            write(12,*)'      levcfg,imcon,natms: ',levcfg,imcon,natms 
!            write(15,'(I10,I10,I16)')levcfg,imcon,natms 
!!            write(12,*)'      lattice vecs:'
!            do i=1,3
!              read(14,*) vecs(i,1:3)
!!              write(12,'(3(F20.10))') vecs(i,1:3)
!              write(15,'(3(F20.10))') vecs(i,1:3)
!            end do
!            do i=1,3
!              do j=1,3
!                matrix(i,j)=vecs(j,i)
!              end do
!            end do
!!            write(12,*)'      atomic positions: ' 
!            do iatom=1,natms
!              read(14,'(A8,I)') name,atnum
!!              write(12,'(A8,I)') name,atnum
!              write(15,'(A8,I)') name,atnum
!              read(14,*) xfrac(1:3)
!!              write(12,'(3(F20.10))') xfrac(1:3)
!              do i=1,3
!                xabs(i)=0.0D0
!                do j=1,3
!                  xabs(i)=xabs(i)+matrix(i,j)*xfrac(j)
!                end do
!              end do
!              write(15,'(F16.9,2(F20.9))') xabs(1:3)
!            end do
!            close(14)
!            close(15)
!      case default
!            print*,'        Error: Unknown filetype.'
!            write(12,*) '       Error: Unknown filetype.'
!            close(14)
!            close(15)
!            return
!      end select
!      
!      write(12,fsubend)
!      end subroutine
!
c---------------------------------------------------------------------

      subroutine config2coorat()
      ! transforms DLPOLY CONFIG file into MBPP COORAT file
      
      use defs
      use misc
      implicit none
      character line*100,infile*20,outfile*20,nameat*8
      double precision vecs(1:3,1:3),xabs(1:3),xfrac(1:3)
      integer i,iinf,idum,natoms,nele,j,numele,iele
      logical newele
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: eles(:)

      write(12,fsubstart)trim(adjustl('config2coorat'))

      ! read input filename
      read(10,'(A100)',end=100) line
      line=adjustl(line)
      do while (line(1:1).eq.'#'.or.len_trim(line).le.1)
            read(10,'(A100)',end=100) line
            line=adjustl(line)
      end do
      read(line,*) infile
      ! read output filename
      read(10,'(A100)',end=100) line
      line=adjustl(line)
      do while (line(1:1).eq.'#'.or.len_trim(line).le.1)
            read(10,'(A100)',end=100) line
            line=adjustl(line)
      end do
      read(line,*) outfile
      ! read CONFIG file
      open(50,file=infile,status='old')
      open(51,file=outfile,status='replace')
      read(50,'(A100)',end=110) line      
      read(50,*,end=110) iinf,idum,natoms
      if(iinf.eq.0) write(12,'(8x,"I will read only coordinates.")')
      if(iinf.eq.1) write(12,'(8x,"I will read coordinates 
     &  & velocities .")')
      if(iinf.eq.2) write(12,'(8x,"I will read coordinates,
     &  velocities and forces.")')
      allocate(atoms(natoms))
      do i=1,3
            read(50,*,end=110) vecs(i,1:3)
      end do
      do i=1,natoms
            read(50,*,end=110)atoms(i)%name, idum
            read(50,*,end=110)xabs(1:3)
            do j=1,iinf
                  read(50,'(A100)',end=110) line
            end do
            call abs2frac(xabs,vecs,xfrac)
            atoms(i)%where(1:3)=xfrac(1:3)
      end do
      close(50)
      ! count elements and number of atoms of each element
      atoms%written=.false.
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
      write(12,'(8x,"Found ",I4," elements.")')nele
      write(12,'(8x,"Elements, number of atoms of each element:")')
      do i=1,nele
            write(12,'(8x,A8,I6)') eles(i)%name,eles(i)%howmany
      end do
      !write COORAT file
      do iele=1,nele
            write(51,'(4x,"natom= ",I6," name= ",A8)')
     &       eles(iele)%howmany,eles(iele)%name
            do i=1,natoms
                  if(atoms(i)%name.eq.eles(iele)%name) 
     &               write(51,'(3(F15.10))') atoms(i)%where
            end do
                  
      end do
      !atname=atoms(1)%name
      !do i=1,natoms
      !  if atoms(i)%written.eq.
      !  
      !end do
      close(51)

      write(12,fsubend)
      return

      ! errors:
100   write(12,'(8x,"Error: file INPUT.COFIMA ended unexpectedly")')
      nerr=nerr+1
      return
110   write(12,'(8x,"Error: CONFIG file ended unexpectedly")')
      nerr=nerr+1
      close(50)
      close(51)
      return

      end subroutine

c---------------------------------------------------------------------

      end module
