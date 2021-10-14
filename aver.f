      module aver

      implicit none

      contains

c-----------------------------------------------------------

      subroutine av(infile,informat,firststep, step,laststep)

      use defs
      use readcoords

      implicit none

      integer firststep,step,laststep
      character(len=*) infile,informat

      ! local variables
      character fname*40,timestepchar*20
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer :: natoms,nspecies
      integer iatom,ispecies,istep,nsteps,timestep,ichar1
      double precision :: vecs(1:3,1:3),avvecs(1:3,1:3),vol,avvol
      !character FMT1*1024,FMT2*1024

      ! print greeting
      print fsubstart,"av"
      print '(8x,"Will calculate averaged lattice parameters")'
      print '(8x,"and write them to the file AV.DAT.")'

      ! determine number of files to read 
      nsteps=(laststep-firststep)/step + 1
!      print '(8x,"There should be ",I," files to read.")', nsteps
      WRITE(FMT2,*) nsteps
      WRITE(FMT1,*) '(8x,"There should be ",I',
     & len(trim(FMT2)),'," files to read".)'
      print FMT1, nsteps

      ! loop over timesteps
      avvecs=0.0D0
      avvol=0.0D0
      open(51,file="AV.DAT",status="replace",err=1000)
      write(51,'("# lattice vectors, averaged lattice vectors")')
      do istep=1,nsteps
            timestep=firststep+(istep-1)*step
            write(timestepchar,*) timestep
            timestepchar=adjustl(timestepchar)
            ichar1=len(trim(infile))
            fname(1:ichar1)=trim(infile)
            fname(ichar1+1:ichar1+len(trim(timestepchar)))
     &            =trim(timestepchar)
            ichar1=ichar1+1+len(trim(timestepchar))
            fname(ichar1:ichar1)='.'
            ichar1=ichar1+1
            fname(ichar1:ichar1+len(trim(informat)))=trim(informat)
            print '(8x,"Reading file ",A)',fname 
            ! read in  coordinates ...
            print '(" ")'
            print '(8x,"calling read_coords...")'
            call read_coords(fname,informat,atoms,natoms,species,
     &           nspecies,vecs)
            print '(8x,"Coords read in.")'
            avvecs=avvecs+vecs
            write(51,'(I20,18(F10.5))')timestep,vecs(1,1:3),
     &         vecs(2,1:3),vecs(3,1:3),avvecs(1,1:3)/dble(istep),
     &         avvecs(2,1:3)/dble(istep),avvecs(3,1:3)/dble(istep) 

      end do
      close(51)

      print fsubendext, 'av'
      return

1000  nerr=nerr+1
      print ferrmssg, "File 'AV.DAT' could not be written."
      close(51)
      return

      end subroutine

c---------------------------------------------------------------------

      subroutine traj(infile,informat,firststep,step,laststep,
     &      thespecies,thenumber,nndist,nnana,nnatom)

      use defs
      use readcoords
      use writecoords
      use misc

      implicit none

      integer firststep,step,laststep,thenumber
      character(len=*) infile,informat,thespecies,nnana,nnatom
      double precision nndist

      ! local variables
      character fname*40,timestepchar*20,outfile*40,numchar*8
      character writeformat*256
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer :: natoms,nspecies,atomnumber,i1,i2,i3
      integer,allocatable :: coordnums(:)
      integer iatom,ispecies,istep,nsteps,timestep,ichar1,inumber
      integer i,j,ithespecies
      integer startnum,endnum,currnum
      double precision :: vecs(1:3,1:3),vol,origin(1:3),place(1:3),
     &         abswhere(1:3),distvec(1:3),otherplace(1:3)  
      double precision, allocatable :: allvecs(:,:,:)
      logical isopen12
      ! variables used by "octet":
      character octchar*1,xsffilename*50,cfgfilename*50
      integer oct,nncoord
      integer, allocatable :: ntet(:),noct(:)
      logical, allocatable :: isoct(:,:)
      type(atom), allocatable :: octetatoms(:,:)
      !character FMT1*1024,FMT2*1024,FMT3*1024,FMT4*1024
      
      ! starting ...
      INQUIRE (unit=12, opened=isopen12)
      ! print greeting
      if(isopen12) then
        write(12,fsubstart)"traj"
      else
        print fsubstart,"traj"
      end if
      ! determine number of files to read
      if(step.gt.0)nsteps=(laststep-firststep)/step + 1
      if(step.eq.0)nsteps=1
      if(isopen12) then
!        write(12,'(8x,"There should be ",I," files to read.")') nsteps
        WRITE(FMT2,*) nsteps
        WRITE(FMT1,*) '(8x,"There should be ",I',
     &   len(trim(FMT2)),'," files to read".)'
        write(12,FMT1) nsteps
      else
!        print '(8x,"There should be ",I," files to read.")', nsteps
        WRITE(FMT2,*) nsteps
        WRITE(FMT1,*) '(8x,"There should be ",I',
     &   len(trim(FMT2)),'," files to read.")'
        print FMT1, nsteps
      end if
      ! if "octet" is to be performed, set numbers of octa/tetra sites
      ! to zero at the beginning
      allocate(ntet(1:nsteps),noct(1:nsteps))
      ntet=0
      noct=0
      ! loop over atom positions to analyzed
      ! if "thenumber" > 0, analyze only one atom position
      if(thenumber.gt.0) then
            startnum=thenumber
            endnum=thenumber
      end if
      ! get number of species from first coordinate file
      istep=1
      timestep=firststep+(istep-1)*step
      write(timestepchar,*) timestep
      timestepchar=adjustl(timestepchar)
      ichar1=len(trim(infile))
      fname(1:len(fname))=" "
      fname(1:ichar1)=trim(infile)
      fname(ichar1+1:ichar1+len(trim(timestepchar)))
     &      =trim(timestepchar)
      ichar1=ichar1+1+len(trim(timestepchar))
      fname(ichar1:ichar1)='.'
      ichar1=ichar1+1
      fname(ichar1:ichar1+len(trim(informat)))=trim(informat)
      !print '(8x,"Reading file ",A)',fname 
      ! read in  coordinates ...
      !print '(" ")'
      !print '(8x,"calling read_coords...")'
      !print*,"reading ",fname,"with format",informat
      call read_coords(fname,informat,atoms,natoms,species,
     &     nspecies,vecs)
      ! if "thenumber" = -1, analyze all atoms of "thespecies"
      if(thenumber.eq.-1) then
            startnum=1
            ! get the number of atoms of the interesting species
            do ispecies=1,nspecies
              if(trim(adjustl(species(ispecies)%name(1:2))).eq.
     &           trim(adjustl(thespecies(1:2)))) then
                endnum=species(ispecies)%howmany
                ithespecies=ispecies
              end if
            end do
      end if
      ! if octet is to be performed, we store for each atom if
      ! it is octa or tetra, for each timestep, in isoct
      allocate(allvecs(nsteps,1:3,1:3))
      if (nnana.eq."octet") then
            allocate(isoct(nsteps,startnum:endnum))
            isoct=.false.
            allocate(octetatoms(nsteps,startnum:endnum))
            octetatoms%name=thespecies
      end if
      ! loop over atoms
      !print*,"startnum,endnum:",startnum,endnum
      do currnum=startnum,endnum
        !print*,"currnum:",currnum
        ! get name of outputfile (TRAJ.speciesnumber.DAT)
        outfile=" "
        write(numchar(1:8),'(I8)') currnum
        numchar=trim(adjustl(numchar))
        write(outfile,'("TRAJ.",A,A,".DAT")')trim(adjustl(thespecies)),
     &     trim(numchar)
        if(isopen12) then
          write(12,'(8x,"Will get the trajectory of ",A,1x,A)')
     &       thespecies,numchar
          write(12,'(8x,"and write it to the file ",A,".")')
     &       trim(outfile)   
        else
          print '(8x,"Will get the trajectory of ",A,1x,A)',
     &       trim(thespecies),numchar
          print '(8x,"and write it to the file ",A,".")',
     &       trim(outfile)   
        end if
        ! loop over timesteps
        open(57,file=outfile,status="replace",err=1000)
        do istep=1,1
            ! find the right atom 
            inumber=0
            ! initially take some weird atom position. If the position later is
            ! still this one, the atom wasn't found.
            place(1)=111111.111D0
            place(2)=222222.222D0
            place(3)=333333.333D0
            timestep=firststep+(istep-1)*step
            write(timestepchar,*) timestep
            timestepchar=adjustl(timestepchar)
            ichar1=len(trim(infile))
            fname=" "
            fname(1:ichar1)=trim(infile)
            fname(ichar1+1:ichar1+len(trim(timestepchar)))
     &            =trim(timestepchar)
            ichar1=ichar1+1+len(trim(timestepchar))
            fname(ichar1:ichar1)='.'
            ichar1=ichar1+1
            fname(ichar1:ichar1+len(trim(informat)))=trim(informat)
            print*,"reading ",fname
            call read_coords(trim(fname),informat,atoms,natoms,species,
     &           nspecies,vecs)
            ! if octet is to be performed, store vectors in allvecs
            allvecs(istep,1:3,1:3)=vecs(1:3,1:3)
            do iatom=1,natoms
!              if(trim(adjustl(atoms(iatom)%name(1:2))).eq.
!     &           trim(adjustl(thespecies(1:2)))) then
              if(atoms(iatom)%name.eq.thespecies) then
                inumber=inumber+1
                if (inumber.eq.currnum) then
                      call frac2abs(atoms(iatom)%where,vecs,origin)
                      place=origin
                      atomnumber=iatom
                end if
              end if
            end do
            ! if the atom couldn't be found, quit with an error message.
            if (place(1).eq.111111.111D0.and.place(2).eq.222222.222D0
     &          .and.place(3).eq.333333.333D0) goto 2000
            ! get coordination number
            ! subtract one because the atom is counted as its own
            ! neighbor (distance=0)
            if (.not. allocated(coordnums))
     &         allocate(coordnums(1:nspecies+1))
            coordnums=0
            do ispecies=1,nspecies
              if(trim(adjustl(species(ispecies)%name(1:2))).eq.
     &           trim(adjustl(thespecies(1:2)))) then
                coordnums(ispecies)=-1
              end if
            end do
            coordnums(nspecies+1)=-1
            ! loop over atoms to find the neighbors
            do iatom=1,natoms
              do i1=-1,1
                do i2=-1,1
                  do i3=-1,1
                    call frac2abs(atoms(iatom)%where,vecs,otherplace)
                    otherplace=otherplace+dble(i1)*vecs(1,1:3)
     &                   +dble(i2)*vecs(2,1:3)+dble(i3)*vecs(3,1:3)
                    distvec=otherplace-place
                    if(absvec(distvec).le.nndist) then
                       coordnums(nspecies+1)=coordnums(nspecies+1)+1
                       do ispecies=1,nspecies
                         if (trim(adjustl(atoms(iatom)%name(1:2))).eq.
     &                     trim(adjustl(species(ispecies)%name(1:2))))
     &                     coordnums(ispecies)=coordnums(ispecies)+1
                       end do
                    end if
                  end do
                end do
              end do
            end do
            write(writeformat,'(A,I0,A)') "('#           timestep,      
     &    x,        y,        z,   (x-x0),   (y-y0),   (z-z0),  total, 
     &  ',",nspecies,"(A3,1x),' neighbors of ',A,1x,A,', oct/tet(1/0)')"
              write(57,writeformat)
     &          adjustr(species(1:nspecies)%name(1:2)),
     &          trim(adjustl(thespecies)),trim(adjustl(numchar))
              writeformat(1:len(writeformat))=' '
!              write(writeformat,'(A,I0,A)')"(I20,6(F10.5),I8,",nspecies,
!       &       "(I8))"
             write(writeformat,'(A,I0,A)')"(I20,6(F10.5),I8,1x,",
     &         nspecies,"(I4),22x,A)"
              ! if "nnana", analyze nearest neighbor sites
              if (nnana.eq."octet") then
                 do ispecies=1,nspecies
                    if(species(ispecies)%name(1:2).eq.nnatom(1:2))
     &                    nncoord=coordnums(ispecies)
                 end do
                 call octet(place,atoms,vecs,species,nnatom,nncoord,oct)
                 write(octchar,'(I1)') oct
                 if(oct.eq.1) noct(istep)=noct(istep)+1
                 if(oct.eq.0) ntet(istep)=ntet(istep)+1
                 ! store oct/tet in isoct
                 !octetatoms(istep,currnum)%where=place
                 call abs2frac(place,vecs,
     &            octetatoms(istep,currnum)%where)
                 allocate(octetatoms(istep,currnum)%properties(1:4))
                 octetatoms(istep,currnum)%properties(1:4)=0.0D0
                 if(oct.eq.1) then
                       isoct(istep,currnum)=.true.
                       octetatoms(istep,currnum)%properties(4)=1.0D0
                 end if
                 write(57,writeformat)timestep,place(1:3),
     &               place(1:3)-origin(1:3),
     &                    coordnums(nspecies+1),
     &                    coordnums(1:nspecies) ,octchar      
              else
                 write(57,writeformat)timestep,place(1:3),
     &                    place(1:3)-origin(1:3),
     &                    coordnums(nspecies+1),
     &                    coordnums(1:nspecies) ,' '      
              end if
        end do
        do istep=2,nsteps
              timestep=firststep+(istep-1)*step
              write(timestepchar,*) timestep
              timestepchar=adjustl(timestepchar)
              ichar1=len(trim(infile))
              fname(1:ichar1)=trim(infile)
              fname(ichar1+1:ichar1+len(trim(timestepchar)))
     &              =trim(timestepchar)
              ichar1=ichar1+1+len(trim(timestepchar))
              fname(ichar1:ichar1)='.'
              ichar1=ichar1+1
              fname(ichar1:ichar1+len(trim(informat)))=trim(informat)
              print*,"reading ",fname
              call read_coords(fname,informat,atoms,natoms,species,
     &             nspecies,vecs)
              ! if octet is to be performed, store vectors in allvecs
              allvecs(istep,1:3,1:3)=vecs(1:3,1:3)
              !print '(8x,"Coords read in.")'
              call frac2abs(atoms(atomnumber)%where,vecs,place)
              ! get coordination number
              ! subtract one because the atom is counted as its own
              ! neighbor (distance=0)
              coordnums=0
              do ispecies=1,nspecies
                if(trim(adjustl(species(ispecies)%name(1:2))).eq.
     &             trim(adjustl(thespecies(1:2)))) then
                  coordnums(ispecies)=-1
                end if
              end do
              coordnums(nspecies+1)=-1
              do iatom=1,natoms
                do i1=-1,1
                  do i2=-1,1
                    do i3=-1,1
                      call frac2abs(atoms(iatom)%where,vecs,otherplace)
                      otherplace=otherplace+dble(i1)*vecs(1,1:3)
     +                     +dble(i2)*vecs(2,1:3)+dble(i3)*vecs(3,1:3)
                      distvec=otherplace-place
                      if(absvec(distvec).le.nndist) then
                         coordnums(nspecies+1)=coordnums(nspecies+1)+1
                         do ispecies=1,nspecies
                           if (trim(adjustl(atoms(iatom)%name(1:2))).eq.
     &                       trim(adjustl(species(ispecies)%name(1:2))))
     &                       coordnums(ispecies)=coordnums(ispecies)+1
                         end do
                      end if
                    end do
                  end do
                end do
              end do
              ! if "nnana", analyze nearest neighbor sites
              if (nnana.eq."octet") then
                 do ispecies=1,nspecies
                  if(species(ispecies)%name(1:2).eq.nnatom(1:2))
     &                  nncoord=coordnums(ispecies)
                 end do
                 call octet(place,atoms,vecs,species,nnatom,nncoord,oct)
                 write(octchar,'(I1)') oct
                 if(oct.eq.1) noct(istep)=noct(istep)+1
                 if(oct.eq.0) ntet(istep)=ntet(istep)+1
                 ! store oct/tet in isoct
                 !octetatoms(istep,currnum)%where=place
                 call abs2frac(place,vecs,
     &            octetatoms(istep,currnum)%where)
                 allocate(octetatoms(istep,currnum)%properties(1:4))
                 octetatoms(istep,currnum)%properties(1:4)=0.0D0
                 if(oct.eq.1) then
                       isoct(istep,currnum)=.true.
                       octetatoms(istep,currnum)%properties(4)=1.0D0
                 end if
                 write(57,writeformat)timestep,place(1:3),
     &                  place(1:3)-origin(1:3),
     &                  coordnums(nspecies+1),
     &                  coordnums(1:nspecies) ,octchar      
              else
                 write(57,writeformat)timestep,place(1:3),
     &                  place(1:3)-origin(1:3),
     &                  coordnums(nspecies+1),
     &                  coordnums(1:nspecies),' '       
              end if
        end do
        close(57)   ! end loop over timesteps
      end do   ! end loop over atoms

      ! "nnana" and "all", print out
      if(nnana.eq."octet") then 
         ! print overall occupations 
         open(57,file="OCTET.DAT",status="replace")
         write(57,'("# timestep, number of octah. sites, number of tetra
     &hedral sites ")')
         do istep=1,nsteps
!            write(57,'(3(I))')(istep-1)*step+firststep,noct(istep),
!     &       ntet(istep)
            WRITE(FMT2,*) (istep-1)*step+firststep
            WRITE(FMT3,*) noct(istep) 
            WRITE(FMT4,*) ntet(istep)
            WRITE(FMT1,*) '(I',len(trim(FMT2)),',I',len(trim(FMT3)),
     &      ',I',len(trim(FMT4)),')'
            write(57, FMT1) (istep-1)*step+firststep,noct(istep),
     &       ntet(istep)
         end do
         close(57)
         ! print cfg files with tetra (0) and octa (1) as properties
         do istep=1,nsteps
            !do currnum=startnum,endnum
              !allocate(octetatoms(istep,currnum)%properties(1:4))
              !octetatoms(istep,currnum)%properties=0.0D0
              !if (isoct(istep,currnum)) then
              !     octetatoms(istep,currnum)%properties(4)=1.0D0
              !end if
            !end do
            cfgfilename(1:len(cfgfilename))=' '
            write(cfgfilename,'("OCTETATOMS.",I0,".cfg")')(istep-1)*step
     &            +firststep       
            !print*,cfgfilename
            !open(57,file=cfgfilename,status="replace")
            !write(57,*) ' CRYSTAL'
            !write(57,*) ' PRIMVEC'
            !do i=1,3
            !  write(57,'(3(F15.10))') allvecs(istep,i,1:3)
            !end do
            !write(57,*) 'PRIMCOORD'
            !write(57,*) noct(istep), '  1'
            !close(57)
            call write_cfg(cfgfilename,
     &           octetatoms(istep,startnum:endnum),
     &           size(octetatoms(istep,startnum:endnum)),
     &           species(ithespecies),1,allvecs(istep,1:3,1:3))
         end do
         ! print separate files with only tet and only oct atoms, resp.
         ! octa atoms:
         do istep=1,nsteps
            xsffilename(1:len(xsffilename))=' '
            write(xsffilename,'("OCTATOMS.",I0,".xsf")')(istep-1)*step
     &            +firststep       
            print*,xsffilename
            open(57,file=xsffilename,status="replace")
            write(57,*) ' CRYSTAL'
            write(57,*) ' PRIMVEC'
            do i=1,3
              write(57,'(3(F15.10))') allvecs(istep,i,1:3)
            end do
            write(57,*) 'PRIMCOORD'
            write(57,*) noct(istep), '  1'
            do currnum=startnum,endnum
              !write(51,'(A5,3(F15.10))') atoms(i)%name,coords
              if (isoct(istep,currnum)) then
                do j=1,100
                  if (octetatoms(istep,currnum)%name(1:2)
     &             .eq.elements(j)) then
                    call frac2abs(octetatoms(istep,currnum)%where,
     &                allvecs(istep,1:3,1:3),abswhere)
                    write(57,'(I5,3(F15.10))') j,abswhere
                  end if
                end do
              end if
            end do
            close(57)
         end do
         ! tetra atoms:
         do istep=1,nsteps
            xsffilename(1:len(xsffilename))=' '
            write(xsffilename,'("TETATOMS.",I0,".xsf")')(istep-1)*step
     &            +firststep       
            open(57,file=xsffilename,status="replace")
            write(57,*) ' CRYSTAL'
            write(57,*) ' PRIMVEC'
            do i=1,3
              write(57,'(3(F15.10))') allvecs(istep,i,1:3)
            end do
            write(57,*) 'PRIMCOORD'
            write(57,*) ntet(istep), '  1'
            do currnum=startnum,endnum
              !write(51,'(A5,3(F15.10))') atoms(i)%name,coords
              if (.not.isoct(istep,currnum)) then
                do j=1,100
                  if (octetatoms(istep,currnum)%name(1:2)
     &             .eq.elements(j)) then
                    call frac2abs(octetatoms(istep,currnum)%where,
     &                allvecs(istep,1:3,1:3),abswhere)
                   write(57,'(I5,3(F15.10))') j,abswhere
                 end if
                end do
              end if
            end do
            close(57)
         end do
       end if

      if(isopen12) then
        write(12,fsubendext) 'traj'
      else
        print fsubendext, 'traj'
      end if  
      return

1000  nerr=nerr+1
      if(isopen12) then
        write(12,ferrmssg) "File TRAJ.....DAT could not be written."
      else
        print ferrmssg, "File TRAJ.....DAT could not be written."
      end if
      close(51)
      return

2000  nerr=nerr+1
      if(isopen12) then
        write(12,ferrmssg) "The atom you specifiedcould not be found."
        write(12,'("Check species name and/or atom number.")')
      else
        print ferrmssg, "The atom you specifiedcould not be found."
        print '("Check species name and/or atom number.")'
      end if
      close(51)
      return

      end subroutine

c---------------------------------------------------------------------

      end module
