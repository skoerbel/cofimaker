      module writecoords

      implicit none

      contains

c---------------------------------------------------------------------

      subroutine write_coords(foutcoords,foutform,atoms,natoms,species,
     &           nspecies,vecs)
      
      use defs

      implicit none
      
      ! variables
      character (LEN=*), intent(in) :: foutcoords,foutform
      double precision vecs(1:3,1:3) !,props(natoms,nprops)
      integer i,natoms
      integer nspecies,nprops
      type(atom) :: atoms(natoms)
      type(element) :: species(nspecies)

      ! local variables
      logical isopen12

      INQUIRE (unit=12, opened=isopen12)
      if(talk) then
        if(isopen12) then
            write(12,fsubstart) trim(adjustl('write_coords'))
        else
            print fsubstart, trim(adjustl('write_coords'))
        end if
      end if
      write(FMT1,*) '(8x,"Writing file ",A',
     &   len(trim(adjustl(foutcoords))),')'
      if(talk) print FMT1,foutcoords

      select case(foutform)
      case("ascii","ASCII")
            call write_ascii(foutcoords,atoms,natoms,species,
     &                        nspecies,vecs)
      case("coorat","COORAT")
            call write_coorat(foutcoords,atoms,natoms,species,nspecies,
     &                        vecs)
      case("xsf","XSF")
            call write_xsf(foutcoords,atoms,natoms,species,nspecies,
     &                        vecs)
      case("xsf_forces","XSF_FORCES")
            call write_xsf_forces(foutcoords,atoms,natoms,species,        &
     &       nspecies,vecs)
      case("gulp","GULP","gin","GIN","res","RES")
            call write_gin(foutcoords,atoms,natoms,species,nspecies,
     &                        vecs)
      case("gulp_short","GULP_SHORT","gin_short","GIN_SHORT",             &
     &     "res_short","RES_SHORT")
            call write_gin_short(foutcoords,atoms,natoms,species,         &
     &      nspecies,vecs)
      case("data","DATA","lmpdata","LMPDATA","lmp","LMP")
            call write_lmpdata(foutcoords,atoms,natoms,species,nspecies,
     &                        vecs)
      case("cfg","CFG")
            call write_cfg(foutcoords,atoms,natoms,species,nspecies,
     &                        vecs)
      case("xyz","XYZ")
            call write_xyz(foutcoords,atoms,natoms,species,nspecies,
     &                        vecs)
      case("cell","CELL","fiesta","FIESTA")
            call write_cell(foutcoords,atoms,natoms,species,nspecies,
     &                        vecs)
      case("cell_old","CELL_OLD","fiesta_old","FIESTA_OLD")
            call write_cell_old(foutcoords,atoms,natoms,species,
     &                        nspecies,vecs)
      case("poscar","POSCAR","vasp","VASP")
            call write_poscar(foutcoords,atoms,natoms,species,
     &                        nspecies,vecs)
      case("poscar_bec","POSCAR_BEC","vasp_bec","VASP_BEC")
            call write_poscar_bec(foutcoords,atoms,natoms,species,
     &                        nspecies,vecs)
      case("poscar_forces","POSCAR_FORCES","vasp_forces","VASP_FORCES")
            call write_poscar_forces(foutcoords,atoms,natoms,species,
     &                        nspecies,vecs)
      case("poscar_flags","POSCAR_FLAGS","vasp_flags","VASP_FLAGS")
            call write_poscar_flags(foutcoords,atoms,natoms,species,
     &                        nspecies,vecs)
      case("poscar_spin","POSCAR_SPIN","vasp_spin","VASP_SPIN")
            call write_poscar_spin(foutcoords,atoms,natoms,species,
     &                        nspecies,vecs)
      case("cif","CIF")
            call write_cif(foutcoords,atoms,species,vecs)
      case("vesta","VESTA")
            call write_vesta(foutcoords,atoms,vecs)
      case default
            call error('Unknown output file format')
      end select
       
      if(talk) then
        if(isopen12) then
            write(12,fsubendext),'write_coords'
        else
            print fsubendext, 'write_coords'
        end if
      end if
      return


      end subroutine

c---------------------------------------------------------------------

      subroutine write_coorat(foutcoords,atoms,natoms,species,nspecies,
     &                        vecs)
      ! writes coordinates to MBPP COORAT file
      
      use defs

      implicit none
      
      ! variables
      character (LEN=*), intent(in) :: foutcoords
      double precision vecs(1:3,1:3)
      integer i,natoms
      integer nspecies
      type(atom) :: atoms(natoms)
      type(element) :: species(nspecies)

      ! local variables
      logical isopen12
      integer ispecies

      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim(adjustl('write_coorat'))
      else
            print fsubstart, trim(adjustl('write_coorat'))
      end if
      open(51,file=foutcoords,status='replace',err=1000)
      !write COORAT file
      do ispecies=1,nspecies
            write(51,'(4x,"natom= ",I6," name= ",A8)',err=1000)
     &       species(ispecies)%howmany,species(ispecies)%name
            do i=1,natoms
                  if(atoms(i)%name.eq.species(ispecies)%name) 
     &               write(51,'(3(F15.10))',err=1000) atoms(i)%where
            end do
                  
      end do
      
      ! write lattice vectors as comment
      do i=1,3
        write(51,'("! cell vecs (Angs) :",3(F20.10))',err=1000) 
     &   vecs(i,1:3)
      end do
      ! write lattice vectors in Bohrs as comment
      do i=1,3
        write(51,'("! cell vecs (Bohr):",3(F20.10))',err=1000) 
     &     vecs(i,1:3)/bohr
      end do
      
      close(51)
      if(isopen12) then
            write(12,fsubendext),'write_coorat'
      else
            print fsubendext, 'write_coorat'
      end if
      return

1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Something went wrong during writing" 
      else
            print ferrmssg,"Something wrong during writing" 
      end if
      close(51)
      stop

      end subroutine

c---------------------------------------------------------------------

      subroutine write_xsf(foutcoords,atoms,natoms,species,nspecies,
     &                        vecs)
      ! writes coordinates to MBPP COORAT file
      
      use defs
      use misc, only : frac2abs

      implicit none
      
      ! variables
      character (LEN=*), intent(in) :: foutcoords
      double precision vecs(1:3,1:3)
      integer i,natoms
      integer nspecies
      type(atom) :: atoms(natoms)
      type(element) :: species(nspecies)

      ! local variables
      logical isopen12
      integer j,ispecies,icart
      double precision coords(1:3)

      INQUIRE (unit=12, opened=isopen12)
      if(talk) then
        if(isopen12) then
            write(12,fsubstart) trim(adjustl('write_xsf'))
        else
            print fsubstart, trim(adjustl('write_xsf'))
        end if
      end if
      open(51,file=foutcoords,status='replace',err=1000)
      !write xsf file
      write(51,*) ' CRYSTAL'
      write(51,*) ' PRIMVEC'
      do i=1,3
        write(51,'(3(F15.10))') vecs(i,1:3)
      end do
      write(51,*) 'PRIMCOORD'
      write(51,*) natoms, '  1'
      do i=1,natoms
        call frac2abs(atoms(i)%where,vecs,coords)
        !write(51,'(A5,3(F15.10))') atoms(i)%name,coords
        do j=1,100
          if (atoms(i)%name(1:2).eq.elements(j))
     &      write(51,'(I5,3(F15.10))') j,coords
        end do
      end do
      close(51)
      if(talk) then
        if(isopen12) then
            write(12,fsubendext),'write_xsf'
        else
            print fsubendext, 'write_xsf'
        end if
      end if
      return

1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Something went wrong during writing" 
      else
            print ferrmssg,"Something wrong during writing" 
      end if
      close(51)
      stop

      end subroutine write_xsf

c---------------------------------------------------------------------

      subroutine write_xsf_forces(foutcoords,atoms,natoms,species,        &
     &     nspecies,vecs)
      ! writes coordinates to MBPP COORAT file
      
      use defs
      use misc, only : frac2abs

      implicit none
      
      ! variables
      character (LEN=*), intent(in) :: foutcoords
      double precision vecs(1:3,1:3)
      integer i,natoms
      integer nspecies
      type(atom) :: atoms(natoms)
      type(element) :: species(nspecies)

      ! local variables
      logical isopen12
      integer j,ispecies,icart
      double precision coords(1:3)

      INQUIRE (unit=12, opened=isopen12)
      if(talk) then
        if(isopen12) then
            write(12,fsubstart) trim(adjustl('write_xsf_forces'))
        else
            print fsubstart, trim(adjustl('write_xsf_forces'))
        end if
      end if
      open(51,file=foutcoords,status='replace',err=1000)
      !write xsf file
      write(51,*) ' CRYSTAL'
      write(51,*) ' PRIMVEC'
      do i=1,3
        write(51,'(3(F15.10))') vecs(i,1:3)
      end do
      write(51,*) 'PRIMCOORD'
      write(51,*) natoms, '  1'
      do i=1,natoms
        call frac2abs(atoms(i)%where,vecs,coords)
        !write(51,'(A5,3(F15.10))') atoms(i)%name,coords
        do j=1,100
          if (atoms(i)%name(1:2).eq.elements(j))
     &      write(51,'(I5,6(F15.10))') j,coords,atoms(i)%force
        end do
      end do
      close(51)
      if(talk) then
        if(isopen12) then
            write(12,fsubendext),'write_xsf_forces'
        else
            print fsubendext, 'write_xsf_forces'
        end if
      end if
      return

1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Something went wrong during writing" 
      else
            print ferrmssg,"Something wrong during writing" 
      end if
      close(51)
      stop

      end subroutine write_xsf_forces

c---------------------------------------------------------------------

      subroutine write_gin(foutcoords,atoms,natoms,species,nspecies,
     &                        vecs)
      ! this subroutine writes coordinates to a gin
      ! file (GULP input file). 
      
      use defs
      implicit none
     
      integer natoms,nspecies
      character(len=*) foutcoords
      double precision vecs(1:3,1:3)
      type(atom)  :: atoms(natoms)
      type(element) :: species(nspecies)

      ! local variables
      logical isopen12,done
      integer i,j
      double precision alat(1:3),angles(1:3)

      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim(adjustl('write_gin'))
      else
            print fsubstart, trim(adjustl('write_gin'))
      end if
      open(51,file=foutcoords,status='replace',err=100)
      write(51,'("#")')
      write(51,'("opti qok")')
      write(51,'("#")')
      write(51,'("#vectors")')
      do i=1,3
        write(51,'("#",3(F22.16))') vecs(i,1:3)
      end do
      write(51,'("#")') 
      write(51,'("cell")')
      do i=1,3
        alat(i)=sqrt(vecs(i,1)**2+vecs(i,2)**2+vecs(i,3)**2)
      end do
      angles(1)=acos((vecs(2,1)*vecs(3,1)+vecs(2,2)*vecs(3,2)
     &          +vecs(2,3)*vecs(3,3))/(alat(2)*alat(3)))*180.0D0/Pi
      angles(2)=acos((vecs(1,1)*vecs(3,1)+vecs(1,2)*vecs(3,2)
     &          +vecs(1,3)*vecs(3,3))/(alat(1)*alat(3)))*180.0D0/Pi
      angles(3)=acos((vecs(1,1)*vecs(2,1)+vecs(1,2)*vecs(2,2)
     &          +vecs(1,3)*vecs(2,3))/(alat(1)*alat(2)))*180.0D0/Pi
      write(51,'(6(F14.9)," 1 1 1 1 1 1")') alat(1:3),angles(1:3)
      write(51,'("frac")')
      !
      ! begin write coordinates
      !
      do i=1,natoms
        write(51,1000) atoms(i)%name(1:2),atoms(i)%core(1:4),             &
     &      atoms(i)%where(1:3),atoms(i)%charge
!     &      atoms(i)%where(1:3),atoms(i)%charge,'1.0 0.0 1 1 1'
      end do
      !
      ! end write coordinates
      !
      ! 
      ! begin write charges
      !
      write(51,'("species",I4)')nspecies*2
      do i=1,nspecies
        done=.false.
        do j=1,natoms
          if (atoms(j)%name(1:2).eq.species(i)%name(1:2).and..not.done)   &
     &    then
            if (atoms(j)%core(1:4).eq."core") write(51,'(A2,5x,A4,F12.6)  &
     &') atoms(j)%name(1:2),atoms(j)%core,atoms(j)%charge
              done=.true.
          end if
        end do 
      end do
      do i=1,nspecies
        done=.false.
        do j=1,natoms
          if (atoms(j)%name(1:2).eq.species(i)%name(1:2).and..not.done)   &
     &    then
            if (atoms(j)%core(1:4).eq."shel") write(51,'(A2,5x,A4,F12.6)  &
     &') atoms(j)%name(1:2),atoms(j)%core,atoms(j)%charge
              done=.true.
          end if
        end do 
      end do
      ! 
      ! end write charges
      !
      close(51)

! 1000 format(A2,2x,A4,2x,3(F20.16),2x,F10.6,1x,A13)
 1000 format(A2,4x,A4,1x,3(F19.16,1x),F10.7,1x,"1.00000 0.00000 1 1 1     &
     &  ")
      if(isopen12) then
            write(12,fsubendext),'write_gin'
      else
            print fsubendext,'write_gin'
      end if
      return
100   nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg),'Error during writing output file.'
      else
            print ferrmssg,'Error during writing output file.'
      end if

      end subroutine write_gin
      
c---------------------------------------------------------------------

      subroutine write_gin_short(foutcoords,atoms,natoms,species,         &
     &                 nspecies,vecs)
      ! this subroutine writes coordinates to a gin
      ! file (GULP input file). 
      
      use defs
      implicit none
     
      integer natoms,nspecies
      character(len=*) foutcoords
      double precision vecs(1:3,1:3)
      type(atom)  :: atoms(natoms)
      type(element) :: species(nspecies)

      ! local variables
      logical isopen12,done
      integer i,j
      double precision alat(1:3),angles(1:3)

      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim(adjustl('write_gin_short'))
      else
            print fsubstart, trim(adjustl('write_gin_short'))
      end if
      open(51,file=foutcoords,status='replace',err=100)
      write(51,'("#")')
      write(51,'("opti qok")')
      write(51,'("#")')
      write(51,'("#vectors")')
      do i=1,3
        write(51,'("#",3(F22.16))') vecs(i,1:3)
      end do
      write(51,'("#")') 
      write(51,'("cell")')
      do i=1,3
        alat(i)=sqrt(vecs(i,1)**2+vecs(i,2)**2+vecs(i,3)**2)
      end do
      angles(1)=acos((vecs(2,1)*vecs(3,1)+vecs(2,2)*vecs(3,2)
     &          +vecs(2,3)*vecs(3,3))/(alat(2)*alat(3)))*180.0D0/Pi
      angles(2)=acos((vecs(1,1)*vecs(3,1)+vecs(1,2)*vecs(3,2)
     &          +vecs(1,3)*vecs(3,3))/(alat(1)*alat(3)))*180.0D0/Pi
      angles(3)=acos((vecs(1,1)*vecs(2,1)+vecs(1,2)*vecs(2,2)
     &          +vecs(1,3)*vecs(2,3))/(alat(1)*alat(2)))*180.0D0/Pi
      write(51,'(6(F11.6)," 1 1 1 1 1 1")') alat(1:3),angles(1:3)
      write(51,'("frac")')
      !
      ! begin write coordinates
      !
      do i=1,natoms
        write(51,1000) atoms(i)%name(1:2),atoms(i)%core(1:4),             &
     &      atoms(i)%where(1:3),atoms(i)%charge
!     &      atoms(i)%where(1:3),atoms(i)%charge,'1.0 0.0 1 1 1'
      end do
      !
      ! end write coordinates
      !
      ! 
      ! begin write charges
      !
      write(51,'("species",I4)')nspecies*2
      do i=1,nspecies
        done=.false.
        do j=1,natoms
          if (atoms(j)%name(1:2).eq.species(i)%name(1:2).and..not.done)   &
     &    then
            if (atoms(j)%core(1:4).eq."core") write(51,'(A2,5x,A4,F12.6)  &
     &') atoms(j)%name(1:2),atoms(j)%core,atoms(j)%charge
              done=.true.
          end if
        end do 
      end do
      do i=1,nspecies
        done=.false.
        do j=1,natoms
          if (atoms(j)%name(1:2).eq.species(i)%name(1:2).and..not.done)   &
     &    then
            if (atoms(j)%core(1:4).eq."shel") write(51,'(A2,5x,A4,F12.6)  &
     &') atoms(j)%name(1:2),atoms(j)%core,atoms(j)%charge
              done=.true.
          end if
        end do 
      end do
      ! 
      ! end write charges
      !
      close(51)

! 1000 format(A2,2x,A4,2x,3(F20.16),2x,F10.6,1x,A13)
 1000 format(A2,4x,A4,1x,3(F9.6,1x),F10.7,1x,"1.00 0.00 1 1 1  ")
      if(isopen12) then
            write(12,fsubendext),'write_gin_short'
      else
            print fsubendext,'write_gin_short'
      end if
      return
100   nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg),'Error during writing output file.'
      else
            print ferrmssg,'Error during writing output file.'
      end if

      end subroutine write_gin_short
      
c---------------------------------------------------------------------

      subroutine write_lmpdata(foutcoords,atoms,natoms,species,nspecies,
     &                        vecs)
      ! writes coordinates to LAMMPS data file
      
      use defs

      implicit none
      
      ! variables
      character (LEN=*), intent(in) :: foutcoords
      double precision vecs(1:3,1:3)
      integer i,natoms
      integer nspecies
      type(atom) :: atoms(natoms)
      type(element) :: species(nspecies)

      ! local variables
      logical isopen12
      integer j,ispecies,icart
      double precision coords(1:3)

      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim(adjustl('write_lmpdata'))
      else
            print fsubstart, trim(adjustl('write_lmpdata'))
      end if
      if (abs(vecs(1,2)).gt.1E-4.or.abs(vecs(1,3)).gt.1E-4.or.
     &    abs(vecs(2,3)).gt.1E-4) then
        print (ferrmssg),"a must be along x, b in xy-plane"
        goto 1000
      end if
      open(51,file=foutcoords,status='replace',err=1000)
      !write LAMMPS data file
      write(51,*) ' coordinate file to be read by LAMMPS'
      write(51,*) ' '
      write(51,'(I8,1x,"atoms")') natoms
      write(51,*) ' '
      write(51,'(I8,1x,"atom types")') nspecies
      write(51,*) ' '
      write(51,'(2(F15.10)," xlo xhi")') 0.0,vecs(1,1)
      write(51,'(2(F15.10)," ylo yhi")') 0.0,vecs(2,2)
      write(51,'(2(F15.10)," zlo zhi")') 0.0,vecs(3,3)
      if (abs(vecs(2,1)).gt.1E-4.or.abs(vecs(3,1)).gt.1E-4.or.
     &    abs(vecs(3,2)).gt.1E-4) then
        write(51,'(3(F10.6)," xy xz yz")')vecs(2,1),vecs(3,1),vecs(3,2)
      end if
      write(51,*) ' '
      write(51,'("Atoms")') 
      write(51,*) ' '
      do i=1,natoms
        do icart=1,3
          coords(icart)=0.0D0
          do j=1,3
            coords(icart)=coords(icart)+atoms(i)%where(j)*vecs(j,icart)
          end do
        end do
        do j=1,nspecies
          if (atoms(i)%name(1:2).eq.species(j)%name(1:2)) ispecies=j
        end do
        write(51,'(I5,I3,F4.1,3(F15.10))') i,ispecies,0.0,coords
      end do
      close(51)
      if(isopen12) then
            write(12,fsubendext),'write_lmpdata'
      else
            print fsubendext, 'write_lmpdata'
      end if
      return

1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Something went wrong during writing" 
      else
            print ferrmssg,"Something wrong during writing" 
      end if
      close(51)
      stop

      end subroutine

c---------------------------------------------------------------------
      
      subroutine write_cfg(outfile,atoms,natoms,species,nspecies,
     &                        vecs)
      ! writes a cfg file readable by atomeye
      
      use defs
      implicit none
      
      integer natoms,nspecies,nprops
      character(len=*) outfile
      type(atom) atoms(natoms)
      type(element) species(nspecies)
      double precision vecs(1:3,1:3)
      !double precision :: props(natoms,nprops)

      ! internal variables
      logical isopen12
      !logical newspecies
      integer i,j
      !double precision alat
      !character line*100,chardum*10,propnames(1:4)*20
      !character FMT1*1024,FMT2*1024

      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim('write_cfg')
      else
            print fsubstart, trim(adjustl('write_cfg'))
      end if
      
      ! determine number of auxiliary properties
      nprops=0
      if(allocated(atoms(1)%properties))
     &   nprops=size(atoms(1)%properties)
      
      !open(55,file=infile,status='old')
      ! read number of atoms
      !read(55,*) natoms
      ! read atom coordinates and properties
      !allocate(atoms(natoms),props(1:natoms,1:4))
      !read(55,'(A100)') line
      !do i=1,natoms
      !      read(55,*) atoms(i)%name,atoms(i)%where,props(i,1:4)
      !end do
      ! count elements and number of atoms of each element
      !atoms(1:natoms)%written=.false.
      !nele=0
      !do i=1,natoms
      !      newele=.true.
      !      do j=1,i-1
      !        if (atoms(i)%name(1:8).eq.atoms(j)%name(1:8))
!     &!         newele=.false.
      !      end do
      !      if (newele) nele=nele+1
      !end do
      !write(12,'(8x,"Found ",I4," elements.")')nele
      !allocate(eles(1:nele))
      !! get names of elements
      !iele=1
      !do i=1,natoms
      !      newele=.true.
      !      do j=1,i-1
      !        if (atoms(i)%name(1:8).eq.atoms(j)%name(1:8))
!     &!         newele=.false.
      !      end do
      !      if (newele) then
      !            eles(iele)%name=atoms(i)%name
      !            iele=iele+1
      !      end if
      !end do
      !! get numbers of elements
      !do iele=1,nele
      !      numele=0
      !      do j=1,natoms
      !            if(atoms(j)%name.eq.eles(iele)%name) numele=numele+1
      !      end do
      !      eles(iele)%howmany=numele
      !end do
      ! read lattice constants, masses,property names from input file
      !read(55,'(A100)') line
      !read(55,*) alat
      !read(55,'(A100)') line
      !do i=1,3
      !      read(55,*) vecs(i,1:3)
      !end do
      !do iele=1,nele
      !      read(55,*) chardum,chardum,eles(iele)%mass
      !end do
      !read(55,'(A100)') line
      !do i=1,4
      !      read(55,*) chardum,idum,propnames(i)
      !end do
      !close(55)
      ! write cfg file
      open(56,file=outfile,status='replace')
!      write(56,'("Number of particles =",I)') natoms
      WRITE(FMT2,*) natoms
      WRITE(FMT1,*) '("Number of particles =",I',len(trim(FMT2)),')'
      write(56,FMT1) natoms 
      write(56,'("A = 1 Angstrom (basic length-scale)")')
      do i=1,3
          write(56,'("H0(",I1,",1) = ",F12.6," A",/,
     &     "H0(",I1,",2) = ",F12.6," A",/,
     &     "H0(",I1,",3) = ",F12.6," A")')
     &     i,vecs(i,1),i,vecs(i,2),i,vecs(i,3)
      end do
      write(56,'(".NO_VELOCITY.")')
!      write(56,'("entry_count = ",I)') nprops+3
      WRITE(FMT2,*) nprops+3
      WRITE(FMT1,*) '("entry_count = ",I',len(trim(FMT2)),')'
      write(56,FMT1) nprops+3 
      do i=1,nspecies
            write(56,'(F10.5)') species(i)%mass
            write(56,'(A8)') species(i)%name
            do j=1,natoms
                  if (atoms(j)%name.eq.species(i)%name)
     &           write(56,1000) atoms(j)%where,
     &             atoms(j)%properties(1:nprops) 
            end do
      end do
      close(56)
 1000 format(3(F15.7),10(F10.5))
      if(isopen12) then
            write(12,fsubendext) 'write_cfg'
      else
            print fsubendext, 'write_cfg'
      end if
      end subroutine
     
c---------------------------------------------------------------------
 
      subroutine write_xyz(outfile,atoms,natoms,species,nspecies,
     &                        vecs)
      ! writes a xyz file
      
      use defs
      use misc
      implicit none
      
      integer natoms,nspecies
      character(len=*) outfile
      type(atom) atoms(natoms)
      type(element) species(nspecies)
      double precision vecs(1:3,1:3)

      ! internal variables
      logical isopen12
      integer i
      !character FMT1*1024,FMT2*1024
      double precision coords(1:3)  

      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim('write_xyz')
      else
            print fsubstart, trim(adjustl('write_xyz'))
      end if
      
      ! write xyz file
      open(56,file=outfile,status='replace')
      WRITE(FMT2,*) natoms
      WRITE(FMT1,*) '(I',len(trim(FMT2)),')'
      write(56,FMT1) natoms 
      write(56,'("Atomic coordinates in Angstrom:")')
      do i=1,natoms
          !call frac2abs(atoms(i)%where,vecs,coords)
          !write(56,'(A4,3(F15.7))') atoms(i)%name(1:4),coords(1:3)
          write(56,'(A4,3(F15.7))') atoms(i)%name(1:4),                   &
     &                              atoms(i)%abswhere(1:3)
      end do
      close(56)
 1000 format(3(F15.7),10(F10.5))
      if(isopen12) then
            write(12,fsubendext) 'write_xyz'
      else
            print fsubendext, 'write_xyz'
      end if
      end subroutine write_xyz

c---------------------------------------------------------------------
 
      subroutine write_cell_old(outfile,atoms,natoms,species,nspecies,
     &                        vecs)
      ! writes a "cell.in" file for fiesta
      
      use defs
      use misc
      implicit none
      
      integer natoms,nspecies
      character(len=*) outfile
      type(atom) atoms(natoms)
      type(element) species(nspecies)
      double precision vecs(1:3,1:3)

      ! internal variables
      logical isopen12
      integer i,j,ispecies
      !character FMT1*1024,FMT2*1024
      double precision coords(1:3)  

      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim('write_cell_old')
      else
            print fsubstart, trim(adjustl('write_cell_old'))
      end if
      
      ! write xyz file
      open(56,file=outfile,status='replace')
      write(56,'("# number of atoms and n species")')
      WRITE(56,'(2x,I0,2x,I0)') natoms,nspecies
      write(56,'("# number of valence bands for xc and kxc")')
      WRITE(56,'(2x,"_nv_")')
      write(56,'("# number of energies points for sigma grid and grid sp
     &acing (in eV) ")')
      WRITE(56,'(2x,I0,2x,F3.1)') 6,0.5
      write(56,'("# number of COHSEX corrected valence and conduction ba
     &nds")')
      WRITE(56,'(2x,I0,2x,I0,2x,"H")') 0,0
      write(56,'("# number of COHSEX iterations, resolution of identity 
     &method, self consistence on wave function, mixing constant (unused
     & if no self consistence)")')
      WRITE(56,'(2x,I0,2x,"V",2x,I0,2x,F3.1)') 0,0,0.5
      write(56,'("# number of GW corrected valence and conduction bands"
     &)') 
      WRITE(56,'(2x,I0,2x,I0)') 20,20
      write(56,'("# number of GW iterations")') 
      WRITE(56,'(2x,I0)') 10
      write(56,'("# dumping for BSE and TDDFT")') 
      WRITE(56,'(2x,I0,2x,I0)') 1,0
      write(56,'("# number of valence and conduction states for bse used
     & only when calling bse)")') 
      WRITE(56,'(2x,I0,2x,I0)') 20,20
      write(56,'("# number of bse exitons and iterations  (used only whe
     &n calling bse)")') 
      WRITE(56,'(2x,I0,2x,I0)') 40,5
      write(56,'("# list of symbols in order")')
      do i=1,nspecies
        write(56,'(2x,A3)') species(i)%name
      end do              
      write(56,'("# scaling factor")')
      write(56,'(" 1.000")')
      write(56,'("# atoms x y z (times scales should give AA) and specie
     &s")')
      do i=1,natoms
        do j=1,nspecies
          if(species(j)%name(1:3).eq.atoms(i)%name(1:3)) ispecies=j
        end do
        write(56,'(2x,3(F14.8),1x,I0,1x,A3)') atoms(i)%abswhere,
     &       ispecies,atoms(i)%name
      end do              

!      WRITE(FMT1,*) '(I',len(trim(FMT2)),')'
!      write(56,FMT1) natoms 
!      write(56,'("Atomic coordinates in Angstrom:")')
!      do i=1,natoms
!          call frac2abs(atoms(i)%where,vecs,coords)
!          write(56,'(A4,3(F15.7))') atoms(i)%name(1:4),coords(1:3)
!      end do
!      close(56)
! 1000 format(3(F15.7),10(F10.5))
      if(isopen12) then
            write(12,fsubendext) 'write_cell_old'
      else
            print fsubendext, 'write_cell_old'
      end if
      end subroutine write_cell_old

c---------------------------------------------------------------------
 
      subroutine write_cell(outfile,atoms,natoms,species,nspecies,
     &                        vecs)
      ! writes a "cell.in" file for fiesta
      
      use defs
      use misc
      implicit none
      
      integer natoms,nspecies
      character(len=*) outfile
      type(atom) atoms(natoms)
      type(element) species(nspecies)
      double precision vecs(1:3,1:3)

      ! internal variables
      logical isopen12
      integer i,j,ispecies
      !character FMT1*1024,FMT2*1024
      double precision coords(1:3)  

      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim('write_cell')
      else
            print fsubstart, trim(adjustl('write_cell'))
      end if
      
      ! write xyz file
      open(56,file=outfile,status='replace')
      write(56,'("# number of atoms and n species")')
      WRITE(56,'(2x,I0,2x,I0)') natoms,nspecies
      write(56,'("# number of valence bands for xc and kxc")')
      WRITE(56,'(2x,"_nv_")')
      write(56,'("# number of energies points for sigma grid and grid sp
     &acing (in eV) ")')
      WRITE(56,'(2x,I0,2x,F3.1)') 6,0.5
      write(56,'("# 0. for Siesta routines calculating LDA Exc on grid o
     &r 1. for NWCHEM Vxcpsi.mat input file")')
      WRITE(56,'(2x,"1")')
      write(56,'("# number of COHSEX corrected valence and conduction ba
     &nds")')
      WRITE(56,'(2x,I0,2x,I0,2x,"H")') 0,0
      write(56,'("# number of COHSEX iterations, resolution of identity 
     &method, self consistence on wave function, mixing constant (unused
     & if no self consistence)")')
      WRITE(56,'(2x,I0,2x,"V",2x,I0,2x,F3.1)') 0,0,0.5
      write(56,'("# number of GW corrected valence and conduction bands"
     &)') 
      WRITE(56,'(2x,I0,2x,I0)') 20,20
      write(56,'("# number of GW iterations")') 
      WRITE(56,'(2x,I0)') 10
      write(56,'("# dumping for BSE and TDDFT")') 
      WRITE(56,'(2x,I0,2x,I0)') 1,0
      write(56,'("# number of valence and conduction states for bse used
     & only when calling bse)")') 
      WRITE(56,'(2x,I0,2x,I0)') 20,20
      write(56,'("# number of bse exitons and iterations  (used only whe
     &n calling bse)")') 
      WRITE(56,'(2x,I0,2x,I0)') 40,5
      write(56,'("# list of symbols in order")')
      do i=1,nspecies
        write(56,'(2x,A3)') species(i)%name
      end do              
      write(56,'("# scaling factor")')
      write(56,'(" 1.000")')
      write(56,'("# atoms x y z (times scales should give AA) and specie
     &s")')
      do i=1,natoms
        do j=1,nspecies
          if(species(j)%name(1:3).eq.atoms(i)%name(1:3)) ispecies=j
        end do
        write(56,'(2x,3(F14.8),1x,I0,1x,A3)') atoms(i)%abswhere,
     &       ispecies,atoms(i)%name
      end do              

!      WRITE(FMT1,*) '(I',len(trim(FMT2)),')'
!      write(56,FMT1) natoms 
!      write(56,'("Atomic coordinates in Angstrom:")')
!      do i=1,natoms
!          call frac2abs(atoms(i)%where,vecs,coords)
!          write(56,'(A4,3(F15.7))') atoms(i)%name(1:4),coords(1:3)
!      end do
!      close(56)
! 1000 format(3(F15.7),10(F10.5))
      if(isopen12) then
            write(12,fsubendext) 'write_cell'
      else
            print fsubendext, 'write_cell'
      end if
      end subroutine write_cell

c---------------------------------------------------------------------

      subroutine write_poscar(foutcoords,atoms,natoms,species,nspecies,
     &                        vecs)
      ! writes coordinates to VASP POSCAR file
      
      use defs

      implicit none
      
      ! variables
      character (LEN=*), intent(in) :: foutcoords
      double precision vecs(1:3,1:3)
      integer i,natoms
      integer nspecies
      type(atom) :: atoms(natoms)
      type(element) :: species(nspecies)
      character string*256

      ! local variables
      logical isopen12
      integer j,ispecies,icart
      double precision coords(1:3)

      INQUIRE (unit=12, opened=isopen12)
      if(talk) then
        if(isopen12) then
            write(12,fsubstart) trim(adjustl('write_poscar'))
        else
            print fsubstart, trim(adjustl('write_poscar'))
        end if
      end if
      ! begin write POSCAR file
      open(51,file=foutcoords,status='replace',err=1000)
      write(51,*) 'POSCAR generated with cofima.' ! Title
      ! begin lattice vecs
      write(51,'(F19.14)') 1.0d0
      do i=1,3
        write(51,'(1x,3(F22.16))') vecs(i,1:3)
      end do
      ! end lattice vecs
      string=''
      do i=1,nspecies
        write(string((i-1)*6+1:),'(2x,A3)')                               &
     &      adjustl(trim(adjustl(species(i)%name(1:2))))
      end do
      write(51,'(A128)') string
      string=''
      do i=1,nspecies
        write(string((i-1)*6+1:),'(I6)') species(i)%howmany
      end do
      write(51,'(A128)') string
      write(51,'("direct")') 
      do ispecies=1,nspecies  
        do i=1,natoms
          if (atoms(i)%name(1:2).eq.species(ispecies)%name(1:2))
     &      write(51,'(3(F20.16))') atoms(i)%where
        end do
      end do
      close(51)
      if(talk) then
        if(isopen12) then
            write(12,fsubendext),'write_poscar'
        else
            print fsubendext, 'write_poscar'
        end if
      end if
      return

1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Something went wrong during writing" 
      else
            print ferrmssg,"Something wrong during writing" 
      end if
      close(51)
      stop

      end subroutine write_poscar

c---------------------------------------------------------------------

      subroutine write_poscar_forces(foutcoords,atoms,natoms,species,     &
     &                            nspecies,vecs)
      ! writes coordinates to VASP POSCAR file containing forces
      
      use defs

      implicit none
      
      ! variables
      character (LEN=*), intent(in) :: foutcoords
      double precision vecs(1:3,1:3)
      integer i,natoms
      integer nspecies
      type(atom) :: atoms(natoms)
      type(element) :: species(nspecies)
      character string*128

      ! local variables
      logical isopen12
      integer j,ispecies,icart
      double precision coords(1:3)

      INQUIRE (unit=12, opened=isopen12)
      if(talk) then
        if(isopen12) then
            write(12,fsubstart) trim(adjustl('write_poscar_forces'))
        else
            print fsubstart, trim(adjustl('write_poscar_forces'))
        end if
      end if
      ! begin write POSCAR file
      open(51,file=foutcoords,status='replace',err=1000)
      write(51,*) 'POSCAR generated with cofima.' ! Title
      ! begin lattice vecs
      write(51,'(F19.14)') 1.0
      do i=1,3
        write(51,'(3(F22.16))') vecs(i,1:3)
      end do
      ! end lattice vecs
      string=''
      do i=1,nspecies
        write(string((i-1)*6+1:),'(A6)') trim(adjustl(species(i)%name))
      end do
      write(51,'(A128)') string
      string=''
      do i=1,nspecies
        write(string((i-1)*6+1:),'(I6)') species(i)%howmany
      end do
      write(51,'(A128)') string
      write(51,'("direct")')
      do ispecies=1,nspecies  
        do i=1,natoms
          if (atoms(i)%name(1:2).eq.species(ispecies)%name(1:2))
     &      write(51,'(3(F20.16),5x,3(F12.6))') atoms(i)%where,           &
     &      atoms(i)%force(1:3)
        end do
      end do
      close(51)
      if(talk) then
        if(isopen12) then
            write(12,fsubendext),'write_poscar_forces'
        else
            print fsubendext, 'write_poscar_forces'
        end if
      end if
      return

1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Something went wrong during writing" 
      else
            print ferrmssg,"Something wrong during writing" 
      end if
      close(51)
      stop

      end subroutine write_poscar_forces

c---------------------------------------------------------------------

      subroutine write_poscar_flags(foutcoords,atoms,natoms,species,      &
     &                            nspecies,vecs)
      ! writes coordinates to VASP POSCAR file containing flags
      
      use defs

      implicit none
      
      ! variables
      character (LEN=*), intent(in) :: foutcoords
      double precision vecs(1:3,1:3)
      integer i,natoms
      integer nspecies
      type(atom) :: atoms(natoms)
      type(element) :: species(nspecies)
      character string*128

      ! local variables
      logical isopen12
      integer j,ispecies,icart
      double precision coords(1:3)

      INQUIRE (unit=12, opened=isopen12)
      if(talk) then
        if(isopen12) then
            write(12,fsubstart) trim(adjustl('write_poscar_flags'))
        else
            print fsubstart, trim(adjustl('write_poscar_flags'))
        end if
      end if
      ! begin write POSCAR file
      open(51,file=foutcoords,status='replace',err=1000)
      write(51,*) 'POSCAR generated with cofima.' ! Title
      ! begin lattice vecs
      write(51,'(F19.14)') 1.0
      do i=1,3
        write(51,'(3(F22.16))') vecs(i,1:3)
      end do
      ! end lattice vecs
      string=''
      do i=1,nspecies
        write(string((i-1)*6+1:),'(A6)') trim(adjustl(species(i)%name))
      end do
      write(51,'(A128)') string
      string=''
      do i=1,nspecies
        write(string((i-1)*6+1:),'(I6)') species(i)%howmany
      end do
      write(51,'(A128)') string
      write(51,'("direct")')
      do ispecies=1,nspecies  
        do i=1,natoms
          if (atoms(i)%name(1:2).eq.species(ispecies)%name(1:2))
     &      write(51,'(3(F20.16),5x,3(A2,1x))') atoms(i)%where,           &
     &      atoms(i)%flags(1:3)
        end do
      end do
      close(51)
      if(talk) then
        if(isopen12) then
            write(12,fsubendext),'write_poscar_flags'
        else
            print fsubendext, 'write_poscar_flags'
        end if
      end if
      return

1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Something went wrong during writing" 
      else
            print ferrmssg,"Something wrong during writing" 
      end if
      close(51)
      stop

      end subroutine write_poscar_flags

c---------------------------------------------------------------------

      subroutine write_poscar_bec(foutcoords,atoms,natoms,species,        &
     &                            nspecies,vecs)
      ! writes coordinates to VASP POSCAR file containing BEC
      
      use defs

      implicit none
      
      ! variables
      character (LEN=*), intent(in) :: foutcoords
      double precision vecs(1:3,1:3)
      integer i,natoms
      integer nspecies
      type(atom) :: atoms(natoms)
      type(element) :: species(nspecies)
      character string*128

      ! local variables
      logical isopen12
      integer j,ispecies,icart
      double precision coords(1:3)

      INQUIRE (unit=12, opened=isopen12)
      if(talk) then
        if(isopen12) then
            write(12,fsubstart) trim(adjustl('write_poscar_bec'))
        else
            print fsubstart, trim(adjustl('write_poscar_bec'))
        end if
      end if
      ! begin write POSCAR file
      open(51,file=foutcoords,status='replace',err=1000)
      write(51,*) 'POSCAR generated with cofima.' ! Title
      ! begin lattice vecs
      write(51,'(F19.14)')  1.0
      do i=1,3
        write(51,'(3(F22.16))') vecs(i,1:3)
      end do
      ! end lattice vecs
      string=''
      do i=1,nspecies
        write(string((i-1)*6+1:),'(A6)'), trim(adjustl(species(i)%name))
      end do
      write(51,'(A128)') string
      string=''
      do i=1,nspecies
        write(string((i-1)*6+1:),'(I6)'), species(i)%howmany
      end do
      write(51,'(A128)') string
      write(51,'("direct")')
      do ispecies=1,nspecies  
        do i=1,natoms
          if (atoms(i)%name(1:2).eq.species(ispecies)%name(1:2))
     &      write(51,'(3(F20.16),5x,9(F12.6))') atoms(i)%where,           &
     &      (atoms(i)%BEC(j,1:3),j=1,3)
        end do
      end do
      close(51)
      if(talk) then
        if(isopen12) then
            write(12,fsubendext),'write_poscar_bec'
        else
            print fsubendext, 'write_poscar_bec'
        end if
      end if
      return

1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Something went wrong during writing" 
      else
            print ferrmssg,"Something wrong during writing" 
      end if
      close(51)
      stop

      end subroutine write_poscar_bec

c---------------------------------------------------------------------

      subroutine write_poscar_spin(foutcoords,atoms,natoms,species,       &
     &                            nspecies,vecs)
      ! writes coordinates to VASP POSCAR file containing SPIN
      
      use defs

      implicit none
      
      ! variables
      character (LEN=*), intent(in) :: foutcoords
      double precision vecs(1:3,1:3)
      integer i,natoms
      integer nspecies
      type(atom) :: atoms(natoms)
      type(element) :: species(nspecies)
      character string*128

      ! local variables
      logical isopen12
      integer j,ispecies,icart
      double precision coords(1:3)

      INQUIRE (unit=12, opened=isopen12)
      if(talk) then
        if(isopen12) then
            write(12,fsubstart) trim(adjustl('write_poscar_spin'))
        else
            print fsubstart, trim(adjustl('write_poscar_spin'))
        end if
      end if
      ! begin write POSCAR file
      open(51,file=foutcoords,status='replace',err=1000)
      write(51,*) 'POSCAR generated with cofima.' ! Title
      ! begin lattice vecs
      write(51,'(F19.14)') 1.0d0
      do i=1,3
        write(51,'(3(F22.16))') vecs(i,1:3)
      end do
      ! end lattice vecs
      string=''
      do i=1,nspecies
        write(string((i-1)*6+1:),'(A6)') trim(adjustl(species(i)%name))
      end do
      write(51,'(A128)') string
      string=''
      do i=1,nspecies
        write(string((i-1)*6+1:),'(I6)') species(i)%howmany
      end do
      write(51,'(A128)') string
      write(51,'("direct")')
      do ispecies=1,nspecies  
        do i=1,natoms
          if (atoms(i)%name(1:2).eq.species(ispecies)%name(1:2))
     &      write(51,'(3(F20.16),5x,1(F12.6))') atoms(i)%where,           &
     &      atoms(i)%spin
        end do
      end do
      close(51)
      if(talk) then
        if(isopen12) then
            write(12,fsubendext),'write_poscar_spin'
        else
            print fsubendext, 'write_poscar_spin'
        end if
      end if
      return

1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Something went wrong during writing" 
      else
            print ferrmssg,"Something wrong during writing" 
      end if
      close(51)
      stop

      end subroutine write_poscar_spin

c---------------------------------------------------------------------

      subroutine write_ascii(foutcoords,atoms,natoms,species,nspecies,
     &                        vecs)
      ! writes coordinates to ASE ascii file
      
      use defs
      use misc, only : vecs2cellpars, cellpars2vecs,frac2abs
      use transform, only : transformvectors

      implicit none
      
      ! variables
      character (LEN=*), intent(in) :: foutcoords
      double precision vecs(1:3,1:3)
      integer i,natoms
      integer nspecies
      type(atom) :: atoms(natoms)
      type(element) :: species(nspecies)

      ! local variables
      logical isopen12
      integer j,ispecies,icart
      double precision coords(1:3)
      double precision vecsrot(3,3),cellpars(1:6)
      type(atom), allocatable :: atomsrot(:)

      INQUIRE (unit=12, opened=isopen12)
      if(talk) then
        if(isopen12) then
            write(12,fsubstart) trim(adjustl('write_ascii'))
        else
            print fsubstart, trim(adjustl('write_ascii'))
        end if
      end if
      !
      ! begin write ascii file
      open(51,file=foutcoords,status='replace',err=1000)
      write(51,*) 'ASCII generated with cofima.' ! Title
      !
      ! begin cell dimensions
      call vecs2cellpars(vecs,cellpars)
      call cellpars2vecs(cellpars,vecsrot)
      call transformvectors(vecs,vecsrot,atoms,atomsrot,(/1,1,1/))
      write(51,'(x,3(F12.6,x))') vecsrot(1,1),vecsrot(2,1),vecsrot(2,2)
      write(51,'(x,3(F12.6,x))') vecsrot(3,1),vecsrot(3,2),vecsrot(3,3)
      ! end cell dimensions
      !
      ! begin atom information
      do i=1,natoms
         call frac2abs(atomsrot(i)%where,vecsrot,atomsrot(i)%abswhere)
         write(51,'(3(F15.10),2x,A3)') atomsrot(i)%abswhere,              &
     &         atomsrot(i)%name
      end do
      close(51)
      if(talk) then
        if(isopen12) then
            write(12,fsubendext),'write_ascii'
        else
            print fsubendext, 'write_ascii'
        end if
      end if
      return

1000  continue
      nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "Something went wrong during writing" 
      else
            print ferrmssg,"Something wrong during writing" 
      end if
      close(51)
      stop

      end subroutine write_ascii

c---------------------------------------------------------------------

      subroutine write_cif(foutcoords,atoms,species,vecs)
      ! this subroutine writes coordinates to a cif
      ! file (format like materials project). 
      
      use defs
      use misc
      implicit none
     
      character(len=*) foutcoords
      double precision vecs(1:3,1:3)
      type(atom)  :: atoms(:)
      type(element) :: species(:)

      ! local variables
      integer natoms,nspecies
      logical isopen12
      integer i,Z, ispecies
      !double precision alat(1:3),angles(1:3)
      double precision cellpars(1:6), vol
      character formula_sum*128,formula_struc*128

      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim(adjustl('write_cif'))
      else
            print fsubstart, trim(adjustl('write_cif'))
      end if

      natoms=size(atoms)
      nspecies=size(species)

      call vecs2cellpars(vecs, cellpars)
      open(51,file=foutcoords,status='replace',err=100)
      !
      ! begin get chemical_formula_sum
      !
      formula_sum="_chemical_formula_sum   '"
      do ispecies=1,nspecies
        if(ispecies==1) then
          write(formula_sum(len_trim(adjustl(formula_sum))+1:),'(2A)')    &
     &     adjustl(trim(adjustl(species(ispecies)%name(1:2))))
        else
          write(formula_sum(len_trim(adjustl(formula_sum))+2:),'(2A)')    &
     &     adjustl(trim(adjustl(species(ispecies)%name(1:2))))
        end if
        write(formula_sum(len_trim(formula_sum)+1:),'(I0)')               &
     &     species(ispecies)%howmany
      end do
      write(formula_sum(len_trim(formula_sum)+1:),'(1A)') "'"
      !
      ! end get chemical_formula_sum
      !
      !
      ! begin get chemical_formula_structural
      !
      Z=gcd(species(1:nspecies)%howmany) 
      formula_struc="_chemical_formula_structural    "
      do ispecies=1,nspecies
        if(ispecies==1) then
          write(formula_struc(len_trim(adjustl(formula_struc))+5:),       &
     &'(2A)') adjustl(trim(adjustl(species(ispecies)%name(1:2))))
        else
          write(formula_struc(len_trim(adjustl(formula_struc))+1:),       &
     &'(2A)') adjustl(trim(adjustl(species(ispecies)%name(1:2))))
        end if
        if (species(ispecies)%howmany.gt.Z) then
          write(formula_struc(len_trim(formula_struc)+1:),'(I0)')         &
     &     species(ispecies)%howmany/Z
        end if
      end do
      !
      ! end get chemical_formula_structural
      !
      call vecs2vol(vecs,vol)

      !
      ! begin write
      !
      write(51,'("# generated using cofima")')
      write(51,'("data_",128A)') trim(adjustl(formula_struc(32:))) 
      write(51,'(A38)'),"_symmetry_space_group_name_H-M   'P 1'"
      write(51,'("_cell_length_a",F13.8)') cellpars(1)
      write(51,'("_cell_length_b",F13.8)') cellpars(2)
      write(51,'("_cell_length_c",F13.8)') cellpars(3)
      write(51,'("_cell_angle_alpha",F13.8)') cellpars(4)
      write(51,'("_cell_angle_beta",F13.8)') cellpars(5)
      write(51,'("_cell_angle_gamma",F13.8)') cellpars(6)
      write(51,'("_symmetry_Int_Tables_number   1")')
      write(51,'(256A)') adjustl(trim(adjustl(formula_struc)))
      write(51,'(256A)') adjustl(trim(adjustl(formula_sum)))
      write(51,'("_cell_volume   ", F16.9)') vol
      write(51,'("_cell_formula_units_Z   ", I0)') Z
      write(51,'("loop_")')
      write(51,'(" _symmetry_equiv_pos_site_id")')
      write(51,'(" _symmetry_equiv_pos_as_xyz")')
      write(51,'(14A)')  "  1  'x, y, z'"
      write(51,'("loop_")')
      write(51,'(" _atom_site_type_symbol")')
      write(51,'(" _atom_site_label")')
      write(51,'(" _atom_site_symmetry_multiplicity")')
      write(51,'(" _atom_site_fract_x")')
      write(51,'(" _atom_site_fract_y")')
      write(51,'(" _atom_site_fract_z")')
      write(51,'(" _atom_site_occupancy")')
      do i=1,natoms
        write(51,1000) atoms(i)%name(1:2),                                &
     &    adjustl(trim(atoms(i)%name(1:2))),i, 1,                         &
     &    atoms(i)%where(1:3),1
      end do
      !
      ! end write
      !

      close(51)

 1000 format(2x,A2,2x,A2,I0,2x,I0, 3(F12.6),2x,I0)
      if(isopen12) then
            write(12,fsubendext),'write_cif'
      else
            print fsubendext,'write_cif'
      end if
      return
100   nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg),'Error during writing output file.'
      else
            print ferrmssg,'Error during writing output file.'
      end if
      end subroutine write_cif
      
c---------------------------------------------------------------------

      subroutine write_vesta(foutcoords,atoms,vecs)
      
      use defs
      use misc
      implicit none
     
      character(len=*) foutcoords
      double precision vecs(1:3,1:3)
      type(atom)  :: atoms(:)
      !type(element) :: species(:)

      ! local variables
      integer natoms
      logical isopen12
      integer i
      double precision cellpars(1:6)

      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) trim(adjustl('write_vesta'))
      else
            print fsubstart, trim(adjustl('write_vesta'))
      end if

      natoms=size(atoms)

      call vecs2cellpars(vecs, cellpars)
      open(51,file=foutcoords,status='replace',err=100)
      write(51,'("#VESTA_FORMAT_VERSION 3.3.0")')
      write(51,'("")')
      write(51,'("")')
      write(51,'("CRYSTAL")')
      write(51,'("")')
      write(51,'("TITLE")')
      write(51,'(" VESTA file generated with cofima.")')
      write(51,'("")')
      write(51,'("GROUP")')
      write(51,'("1 1 P 1")')
      write(51,'("SYMOP")')
      write(51,'(" 0.000000  0.000000  0.000000  1  0  0   0  1  0   0    &
     &0  1   1")')
      write(51,'(" -1.0 -1.0 -1.0  0 0 0  0 0 0  0 0 0")')
      write(51,'("TRANM 0")')
      write(51,'(" 0.000000  0.000000  0.000000  1  0  0   0  1  0   0    &
     &0  1")')
      write(51,'("LTRANSL")')
      write(51,'(" -1")')
      write(51,'(" 0.000000  0.000000  0.000000  0.000000  0.000000  0.0  &
     &00000")')
      write(51,'("LORIENT")')
      write(51,'(" -1   0   0   0   0")')
      write(51,'(" 1.000000  0.000000  0.000000  1.000000  0.000000  0.0  &
     &00000")')
      write(51,'(" 0.000000  0.000000  1.000000  0.000000  0.000000  1.0  &
     &00000")')
      write(51,'("LMATRIX")')
      write(51,'(" 1.000000  0.000000  0.000000  0.000000")')
      write(51,'(" 0.000000  1.000000  0.000000  0.000000")')
      write(51,'(" 0.000000  0.000000  1.000000  0.000000")')
      write(51,'(" 0.000000  0.000000  0.000000  1.000000")')
      write(51,'(" 0.000000  0.000000  0.000000")')
      write(51,'("CELLP")')
      write(51,'(5(F10.6,1x),F10.6)') cellpars(1:6)
      write(51,'("  0.000000   0.000000   0.000000   0.000000   0.000000  &
     &   0.000000")')
      write(51,'("STRUC")')
      do i=1,natoms
        !print*, atoms(i)%label,atoms(i)%occ 
        write(51,'(I3,A3,A11,F8.4,3(F11.6),4x,"1a",7x,"1")') i,           &
     &        adjustr(atoms(i)%name(1:2)),adjustr(atoms(i)%label(1:11)),  &
     &        atoms(i)%occ,atoms(i)%where
        write(51,'("                            0.000000   0.000000   0.  &
     &000000  0.00")')
      end do
      write(51,'("  0 0 0 0 0 0 0")')
      !
      ! end write
      !

      close(51)

 1000 format(2x,A2,2x,A2,I0,2x,I0, 3(F12.6),2x,I0)
      if(isopen12) then
            write(12,fsubendext),'write_vesta'
      else
            print fsubendext,'write_vesta'
      end if
      return
100   nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg),'Error during writing output file.'
      else
            print ferrmssg,'Error during writing output file.'
      end if
      end subroutine write_vesta
      
c---------------------------------------------------------------------

      end module
