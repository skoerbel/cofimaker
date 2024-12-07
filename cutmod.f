      module cutmod

      contains

      subroutine cut(cutinfile,cutoutfile,cutinform,cutoutform,
     &           above,axisint,cuthere,origin,nndist,hdist,hydr)

      use defs
      use readcoords
      use writecoords
      use misc
      implicit none

      ! variables
      integer axisint
      double precision cuthere,origin(1:3),nndist,hdist
      character(len=*) cutinfile,cutoutfile,cutinform,cutoutform,
     &                 above
      logical, optional :: hydr

      ! internal variables
      logical isopen12,newspecies
      type(atom), allocatable :: atoms(:),keepatoms(:),rmatoms(:)
      type(atom),allocatable :: atomsH(:),keepatomsH(:)
      type(element), allocatable :: species(:),keepspecies(:),
     &      keepspeciesH(:)
      integer natoms,nspecies,nkeep,nspeckeep,
     &      nspeckeepH
      integer nremove,jrm
      double precision vecs(1:3,1:3),distvec(1:3),dist,gvecs(1:3,1:3)
      double precision coords1(3),coords2(3)
      integer i,j,k,l,ndangb !  i1,i2,i3,n_NN
      integer k1,k2,k3
      !character FMT1*1024,FMT2*1024,FMT3*1024
      character filename*256

      talk=.true.       
      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            if(talk)write(12,fsubstart) trim('cut')
      else
            if(talk)print fsubstart, trim(adjustl('cut'))
      end if

      ! read in  coordinates ...
      if (talk) then
            print '(" ")'
            print '(8x,"calling read_coords...")'
      end if
      call read_coords(cutinfile,cutinform,atoms,natoms,species,
     &           nspecies,vecs)
      if(talk)print '(8x,"Coords read in.")'

      ! remove atoms below/above threshold
      select case(above)
      case("above","ABOVE")
            nkeep=0
            do i=1,natoms
              if(atoms(i)%where(axisint).le.cuthere) then
                nkeep=nkeep+1
              end if
            end do
            nremove=natoms-nkeep
            allocate(keepatoms(nkeep),rmatoms(nremove))
            j=0
            jrm=0
            do i=1,natoms
              if(atoms(i)%where(axisint).le.cuthere) then
                j=j+1
                keepatoms(j)=atoms(i)
              end if
              if(atoms(i)%where(axisint).gt.cuthere) then
                jrm=jrm+1
                rmatoms(jrm)=atoms(i)
              end if
            end do
      case("below","BELOW")
            nkeep=0
            do i=1,natoms
              if(atoms(i)%where(axisint).ge.cuthere) then
                nkeep=nkeep+1
              end if
            end do
            nremove=natoms-nkeep
            allocate(keepatoms(nkeep),rmatoms(nremove))
            j=0
            jrm=0
            do i=1,natoms
              if(atoms(i)%where(axisint).gt.cuthere) then
                j=j+1
                keepatoms(j)=atoms(i)
              end if
              if(atoms(i)%where(axisint).le.cuthere) then
                jrm=jrm+1
                rmatoms(jrm)=atoms(i)
              end if
            end do
      case("sphere","SPHERE")
            nkeep=0
            select case(axisint)
            case(0)
              ! keep atoms inside sphere  
              do i=1,natoms
                !call frac2abs(atoms(i)%where,vecs,distvec)
                distvec=atoms(i)%abswhere 
                distvec=distvec-origin
                dist=absvec(distvec)
                if(dist.le.cuthere) then
                  nkeep=nkeep+1
                end if
              end do
              nremove=natoms-nkeep
              allocate(keepatoms(nkeep),rmatoms(nremove))
              j=0
              jrm=0
              do i=1,natoms
                !call frac2abs(atoms(i)%where,vecs,distvec)
                distvec=atoms(i)%abswhere 
                distvec=distvec-origin
                dist=absvec(distvec)
                if(dist.le.cuthere) then
                  j=j+1
                  keepatoms(j)=atoms(i)
                end if
                if(dist.gt.cuthere) then
                  jrm=jrm+1
                  rmatoms(jrm)=atoms(i)
                end if
              end do
            case(1) 
              ! keep atoms outside sphere  
              do i=1,natoms
                !call frac2abs(atoms(i)%where,vecs,distvec)
                distvec=atoms(i)%abswhere 
                distvec=distvec-origin
                dist=absvec(distvec)
                if(dist.gt.cuthere) then
                  nkeep=nkeep+1
                end if
              end do
              nremove=natoms-nkeep
              allocate(keepatoms(nkeep),rmatoms(nremove))
              j=0
              jrm=0
              do i=1,natoms
                !call frac2abs(atoms(i)%where,vecs,distvec)
                distvec=atoms(i)%abswhere 
                distvec=distvec-origin
                dist=absvec(distvec)
                if(dist.gt.cuthere) then
                  j=j+1
                  keepatoms(j)=atoms(i)
                end if
                if(dist.le.cuthere) then
                  jrm=jrm+1
                  rmatoms(jrm)=atoms(i)
                end if
              end do
            case default
              nerr=nerr+1
              print ferrmssg, 
     &           "Unknown parameter. Use 'above' or 'below'"
              if(isopen12) write(12,ferrmssg) 
     &           "Unknown parameter. Use 'above' or 'below'"     
              return
            end select
      case("plane","PLANE")
            nkeep=0
            if (talk) print'(8x,"normal vector:",3(F12.6))',origin
            select case(axisint)
            case(0)
              ! keep atoms below plane  
              ! "radius" is now distance from plane
              ! "origin" is now n1,n2,n3 (vector normal to plane)
              do i=1,natoms
                distvec=atoms(i)%abswhere
                distvec=distvec-cuthere*origin/norm2(origin)
                dist=dot_product(distvec,origin/norm2(origin))
                if(dist.le.0) then
                  nkeep=nkeep+1
                end if
              end do
              nremove=natoms-nkeep
              allocate(keepatoms(nkeep),rmatoms(nremove))
              j=0
              jrm=0
              do i=1,natoms
                distvec=atoms(i)%abswhere 
                distvec=distvec-cuthere*origin/norm2(origin)
                dist=dot_product(distvec,origin/norm2(origin))
                if(dist.le.0) then
                  j=j+1
                  keepatoms(j)=atoms(i)
                end if
                if(dist.gt.0) then
                  jrm=jrm+1
                  rmatoms(jrm)=atoms(i)
                end if
              end do
            case(1) 
              ! keep atoms above plane 
              do i=1,natoms
                distvec=atoms(i)%abswhere 
                distvec=distvec-cuthere*origin/norm2(origin)
                dist=dot_product(distvec,origin/norm2(origin))
                if(dist.gt.0) then
                  nkeep=nkeep+1
                end if
              end do
              nremove=natoms-nkeep
              allocate(keepatoms(nkeep),rmatoms(nremove))
              j=0
              jrm=0
              do i=1,natoms
                distvec=atoms(i)%abswhere 
                distvec=distvec-cuthere*origin/norm2(origin)
                dist=dot_product(distvec,origin/norm2(origin))
                if(dist.gt.0) then
                  j=j+1
                  keepatoms(j)=atoms(i)
                end if
                if(dist.le.0) then
                  jrm=jrm+1
                  rmatoms(jrm)=atoms(i)
                end if
              end do
            case default
              nerr=nerr+1
              print ferrmssg, 
     &           "Unknown parameter. cofima --help prints possible param  &
     &eters"
              if(isopen12) write(12,ferrmssg) 
     &           "Unknown parameter. cofima --help prints possible param  &
     &eters"
              return
            end select
      case("dirplane","DIRPLANE")
            nkeep=0
            origin=origin(1)*vecs(1,:)+origin(2)*vecs(2,:)                &
     &            +origin(3)*vecs(3,:)
            if (talk) print'(8x,"normal vector:",3(F12.6))',origin
            select case(axisint)
            case(0)
              ! keep atoms below plane  
              ! "radius" is now distance from plane
              ! "origin" is now n1,n2,n3 (vector normal to plane)
              do i=1,natoms
                distvec=atoms(i)%abswhere
                distvec=distvec-cuthere*origin/norm2(origin)
                dist=dot_product(distvec,origin/norm2(origin))
                if(dist.le.0) then
                  nkeep=nkeep+1
                end if
              end do
              nremove=natoms-nkeep
              allocate(keepatoms(nkeep),rmatoms(nremove))
              j=0
              jrm=0
              do i=1,natoms
                distvec=atoms(i)%abswhere 
                distvec=distvec-cuthere*origin/norm2(origin)
                dist=dot_product(distvec,origin/norm2(origin))
                if(dist.le.0) then
                  j=j+1
                  keepatoms(j)=atoms(i)
                end if
                if(dist.gt.0) then
                  jrm=jrm+1
                  rmatoms(jrm)=atoms(i)
                end if
              end do
            case(1) 
              ! keep atoms above plane 
              do i=1,natoms
                distvec=atoms(i)%abswhere 
                distvec=distvec-cuthere*origin/norm2(origin)
                dist=dot_product(distvec,origin/norm2(origin))
                if(dist.gt.0) then
                  nkeep=nkeep+1
                end if
              end do
              nremove=natoms-nkeep
              allocate(keepatoms(nkeep),rmatoms(nremove))
              j=0
              jrm=0
              do i=1,natoms
                distvec=atoms(i)%abswhere 
                distvec=distvec-cuthere*origin/norm2(origin)
                dist=dot_product(distvec,origin/norm2(origin))
                if(dist.gt.0) then
                  j=j+1
                  keepatoms(j)=atoms(i)
                end if
                if(dist.le.0) then
                  jrm=jrm+1
                  rmatoms(jrm)=atoms(i)
                end if
              end do
            case default
              nerr=nerr+1
              print ferrmssg, 
     &           "Unknown parameter. cofima --help prints possible param  &
     &eters"
              if(isopen12) write(12,ferrmssg) 
     &           "Unknown parameter. cofima --help prints possible param  &
     &eters"
              return
            end select
      case("hklplane","HKLPLANE")
            nkeep=0
            call recipr_latt_vecs(vecs,gvecs)  
            if (talk) print '(8x,"rec. latt. vecs:",/,8x,3(F12.6),/,8x,   &
     &        3(F12.6),/,8x,3(F12.6))',gvecs(1,:),gvecs(2,:),gvecs(3,:)
              if (talk) print '(8x,"h,k,l=",3(F12.6))',origin
            origin=origin(1)*gvecs(1,:)+origin(2)*gvecs(2,:)              &
     &            +origin(3)*gvecs(3,:)
            if (talk) print'(8x,"normal vector:",3(F12.6))',origin
            select case(axisint)
            case(0)
              ! keep atoms below plane  
              ! "radius" is now distance from plane
              ! "origin" is now h*b1,k*b2,l*b3 (b recipr latt vecs, vector normal to plane)
              do i=1,natoms
                distvec=atoms(i)%abswhere
                distvec=distvec-cuthere*origin/norm2(origin)
                dist=dot_product(distvec,origin/norm2(origin))
                if(dist.le.0) then
                  nkeep=nkeep+1
                end if
              end do
              nremove=natoms-nkeep
              allocate(keepatoms(nkeep),rmatoms(nremove))
              j=0
              jrm=0
              do i=1,natoms
                distvec=atoms(i)%abswhere 
                distvec=distvec-cuthere*origin/norm2(origin)
                dist=dot_product(distvec,origin/norm2(origin))
                if(dist.le.0) then
                  j=j+1
                  keepatoms(j)=atoms(i)
                end if
                if(dist.gt.0) then
                  jrm=jrm+1
                  rmatoms(jrm)=atoms(i)
                end if
              end do
            case(1) 
              ! keep atoms above plane 
              do i=1,natoms
                distvec=atoms(i)%abswhere 
                distvec=distvec-cuthere*origin/norm2(origin)
                dist=dot_product(distvec,origin/norm2(origin))
                if(dist.gt.0) then
                  nkeep=nkeep+1
                end if
              end do
              nremove=natoms-nkeep
              allocate(keepatoms(nkeep),rmatoms(nremove))
              j=0
              jrm=0
              do i=1,natoms
                distvec=atoms(i)%abswhere 
                distvec=distvec-cuthere*origin/norm2(origin)
                dist=dot_product(distvec,origin/norm2(origin))
                if(dist.gt.0) then
                  j=j+1
                  keepatoms(j)=atoms(i)
                end if
                if(dist.le.0) then
                  jrm=jrm+1
                  rmatoms(jrm)=atoms(i)
                end if
              end do
            case default
              nerr=nerr+1
              print ferrmssg, 
     &           "Unknown parameter. cofima --help prints possible param  &
     &eters"
              if(isopen12) write(12,ferrmssg) 
     &           "Unknown parameter. cofima --help prints possible param  &
     &eters"
              return
            end select
      case default
          nerr=nerr+1
          print ferrmssg, 
     &           "Unknown parameter. cofima --help prints possible param  &
     &eters"
          if(isopen12) write(12,ferrmssg) 
     &           "Unknown parameter. cofima --help prints possible param  &
     &eters"
          return
      end select

      ! determine new numbers of species asf
      ! get number of species
      nspeckeep=0
      do i=1,nkeep
         newspecies=.true.
         do j=1,i-1
           if (keepatoms(j)%name(1:2).eq.
     &         keepatoms(i)%name(1:2))newspecies=.false.
         end do
         if(newspecies) nspeckeep=nspeckeep+1
      end do
      if(talk) then
        if (isopen12) then
!        write(12,'(8x,"Number of elements and atoms kept: ",I,I)'),
!     &        nspeckeep,nkeep  
          WRITE(FMT1,*) '(8x,"Number of elements and atoms kept: ",I',
     &     len(trim(FMT2)),',I',len(trim(FMT3)),')'
          write(12,FMT1)  nspeckeep,nkeep 
        else
!        print '(8x,"Number of elements and atoms kept: ",I,I)',
!     &        nspeckeep,nkeep  
          WRITE(FMT1,*) '(8x,"Number of elements and atoms kept: ",I',
     &     len(trim(FMT2)),',I',len(trim(FMT3)),')'
          print FMT1, nspeckeep,nkeep 
        end if
      end if
      close(21)
      !get name and number of atoms of each species
      allocate(keepspecies(nspeckeep))
      k=1
      do j=1,nkeep
        newspecies=.true.
        do i=1,j-1
          if (keepatoms(i)%name(1:2).eq.keepatoms(j)%name(1:2))
     &        newspecies=.false.
        end do
        if(newspecies) then
            keepspecies(k)%name=keepatoms(j)%name
            k=k+1
        end if
      end do  

      ! get numbers of atoms of each species
      do k=1,nspeckeep
            keepspecies(k)%howmany=0
            do i=1,nkeep
               if(keepatoms(i)%name(1:2).eq.keepspecies(k)%name(1:2)) 
     &            keepspecies(k)%howmany=keepspecies(k)%howmany+1
            end do
      end do


      ! write
      call write_coords(cutoutfile,cutoutform,keepatoms,nkeep,
     &                  keepspecies,nspeckeep,vecs)
      
      if (hydr) then
       !!! hydrogenate dangling bonds
       allocate(atomsH(natoms))
       atomsH(1:nkeep)=keepatoms(1:nkeep)
       atomsH(nkeep+1:natoms)=rmatoms(1:nremove)
       !atomsH(nkeep+1:natoms)%name(1:2)="H "
       ! 1. determine number of H atoms (=number of dangling bonds)
       ndangb=0
       do i=1,nkeep
         call frac2abs(atomsH(i)%where,vecs,coords1)
         do j=nkeep+1,natoms
           call frac2abs(atomsH(j)%where,vecs,coords2)
           do k1=-1,1,1
           do k2=-1,1,1
           do k3=-1,1,1
             distvec=coords2-coords1+dble(k1)*vecs(1,:)                   &
     &            +dble(k2)*vecs(2,:)+dble(k3)*vecs(3,:)
             if (absvec(distvec).le.nndist) ndangb=ndangb+1
           end do ! k3
           end do ! k2
           end do ! k1
         end do
       end do 
       if(talk) 
     & print '(8x,"number of dangling bonds to hydrogenize:",I7)',ndangb
       ! place H at each dangling bond
       allocate(keepatomsH(nkeep+ndangb))
       keepatomsH(1:nkeep)=keepatoms(1:nkeep)
       keepatomsH(nkeep+1:nkeep+ndangb)%name(1:2)="H "
       k=nkeep+1
       do i=1,nkeep
         call frac2abs(atomsH(i)%where,vecs,coords1)
         do j=nkeep+1,natoms
           do k1=-1,1,1
           do k2=-1,1,1
           do k3=-1,1,1
             call frac2abs(atomsH(j)%where,vecs,coords2)
             distvec=coords2-coords1+dble(k1)*vecs(1,:)                   &
     &            +dble(k2)*vecs(2,:)+dble(k3)*vecs(3,:)
             if (absvec(distvec).le.nndist) then
               ! place a H at a fraction of the bond distance
               if(hdist.le.0.0D0) coords2=coords1+abs(hdist)*distvec
               ! place a H at a fixed distance (default: 1.25 Angs) along the former bond 
               if(hdist.gt.0.0D0) coords2=coords1+hdist*distvec
     &          /absvec(distvec)
               call abs2frac(coords2,vecs,keepatomsH(k)%where)
!               print*,"coords1,coords2,keepatomsH(k)%where:",
!      &         coords1,coords2,keepatomsH(k)%where
               k=k+1
             end if
           end do ! k3
           end do ! k2
           end do ! k1
         end do !j
       end do ! i
       
       ! write
       call getspecies(keepatomsH,keepspeciesH)
       nspeckeepH=size(keepspeciesH)
       filename(1:256)=" "
       write(filename(1:5),'(A5)') "HYDR."
       write(FMT1,*) '(A',len_trim(cutoutfile),')'
       write(filename(6:5+len_trim(cutoutfile)),FMT1) 
     &   trim(adjustl(cutoutfile))
       !print*,"outfile:",filename
       call write_coords(filename,cutoutform,keepatomsH,nkeep+ndangb,
     &             keepspeciesH,nspeckeepH,vecs)
      end if ! hydr
      
      
      if(talk) then       
        if(isopen12) then
            write(12,fsubendext) 'cut'
        else
            print fsubendext, 'cut'
        end if
      end if
      return
      
      end subroutine cut 

c--------------------------------------------------------------------

      subroutine mergefiles(filename1,format1,filename2,format2,
     &            filename3,format3)
      
      use defs
      use readcoords
      use writecoords
      use misc
      use replace
      implicit none

      ! variables
      character(len=*) filename1,filename2,filename3,format1,format2,
     &       format3

      ! internal variables
      logical isopen12,newspecies
      type(atom), allocatable :: atoms1(:),atoms2(:),atoms3(:)
      type(element), allocatable :: species1(:),species2(:),species3(:)
      integer natoms1,nspecies1,natoms2,nspecies2,
     &    natoms3,nspecies3
      double precision vecs(1:3,1:3)
      integer i,j,k,l
       
      talk=.true.
      ! write output to unit 12, if open, else print to terminal
      INQUIRE (unit=12, opened=isopen12)
      if(talk) then
        if(isopen12) then
            write(12,fsubstart) trim('merge')
        else
            print fsubstart, trim(adjustl('merge'))
        end if
      end if

      ! read in first coordinates ...
            print '(" ")'
            print '(8x,"calling read_coords...")'
            call read_coords(filename1,format1,atoms1,natoms1,species1,
     &           nspecies1,vecs)
            print '(8x,"First coords read in.")'
      ! read in second coordinates ...
            print '(" ")'
            print '(8x,"calling read_coords...")'
            call read_coords(filename2,format2,atoms2,natoms2,species2,
     &           nspecies2,vecs)
            print '(8x,"Second coords read in.")'
            print '(8x,"Merging ....")'

      ! old:
      !natoms3=natoms1+natoms2
      !allocate(atoms3(1:natoms3))
      !atoms3(1:natoms1)=atoms1(1:natoms1)
      !atoms3(natoms1+1:natoms3)=atoms2(1:natoms2)
      ! end old
      ! new: 
      call addatoms(atoms1,atoms2,atoms3)
      natoms3=size(atoms3)
      ! end new
      ! get the number of species of the merged structure 
      call getspecies(atoms3,species3)
      nspecies3=size(species3)

      ! write
      call write_coords(filename3,format3,atoms3,natoms3,species3,
     &                  nspecies3,vecs)
                  
      if(isopen12) then
            write(12,fsubendext) 'merge'
      else
            print fsubendext, 'merge'
      end if
      return
      
      end subroutine 

c---------------------------------------------------------------------

      end module
