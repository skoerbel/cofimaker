        module cluster
        ! contains subroutines that have to do with building clusters or cutting
        ! them out of the bulk 

        implicit none

        contains

c---------------------------------------------------------------------      

      subroutine cageZnO(nZnO)
      ! builds a cage-like cluster of nZno ZnO units.
      
      use defs
      implicit none
      
      integer nZnO
      ! local variables
      integer nZn,nO,i,j
      double precision dZnO,theta,rsphere,coordZn(1:16,1:3),
     &                 coordO(1:16,1:3),phi
      
      print(fsubstart), "cageZnO"
      select case(nZnO)
      case(16)
        nZn=16          ! number of Zn atoms
        nO=16           ! number of O atoms
        dZnO=2.0        ! distance between Zn and O
        
        ! 1. (uppermost) layer (6 atoms)
        theta=(1.0D0*Pi/6.0D0)
        rsphere=dZnO/sin(theta)
        phi=(4.0D0*Pi/6.0D0)
        j=0
        do i=1,3
          j=j+1
          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)))
          coordZn(j,2)=rsphere*sin(theta)*sin(phi*dble(i-1))
          coordZn(j,3)=rsphere*cos(theta)
          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)+0.5D0))
          coordO(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)+0.5D0))
          coordO(j,3)=rsphere*cos(theta)
        end do
        ! 2. layer (6 atoms)
        theta=(3.0D0*Pi/8.0D0)
        phi=(4.0D0*Pi/6.0D0)
        do i=1,3
          j=j+1
          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)))
          coordO(j,2)=rsphere*sin(theta)*sin(phi*dble(i-1))
          coordO(j,3)=rsphere*cos(theta)
          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)+0.5D0))
          coordZn(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)+0.5D0))
          coordZn(j,3)=rsphere*cos(theta)
        end do
        ! 3. layer (8 atoms)
        theta=(3.0D0*Pi/6.0D0)
        phi=(2.0D0*Pi/8.0D0)
        do i=1,2
          j=j+1
          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(2.0D0*dble(i-1)
     &               +0.5D0))
          coordZn(j,2)=rsphere*sin(theta)*sin(phi*(2.0D0*dble(i-1)
     &               +0.5D0))
          coordZn(j,3)=rsphere*cos(theta)
          coordO(j,1)=rsphere*sin(theta)*cos(phi*(2.0D0*dble(i)-0.5D0))
          coordO(j,2)=rsphere*sin(theta)*sin(phi*(2.0D0*dble(i)-0.5D0))
          coordO(j,3)=rsphere*cos(theta)
        end do 
        do i=3,4
          j=j+1
          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(2.0D0*dble(i)
     &                -0.5D0))
          coordZn(j,2)=rsphere*sin(theta)*sin(phi*(2.0D0*dble(i)
     &                -0.5D0))
          coordZn(j,3)=rsphere*cos(theta)
          coordO(j,1)=rsphere*sin(theta)*cos(phi*(2.0D0*dble(i-1)
     &                +0.5D0))
          coordO(j,2)=rsphere*sin(theta)*sin(phi*(2.0D0*dble(i-1)
     &                +0.5D0))
          coordO(j,3)=rsphere*cos(theta)
        end do
        ! 4. layer (6 atoms)
        theta=(5.0D0*Pi/8.0D0)
        phi=(4.0D0*Pi/6.0D0)
        do i=1,3
          j=j+1
          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)))
          coordO(j,2)=rsphere*sin(theta)*sin(phi*dble(i-1))
          coordO(j,3)=rsphere*cos(theta)
          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)+0.5D0))
          coordZn(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)+0.5D0))
          coordZn(j,3)=rsphere*cos(theta)
        end do
        ! 5th (lowest) layer (6 atoms)
        theta=(5.0D0*Pi/6.0D0)
        phi=(4.0D0*Pi/6.0D0)
        do i=1,3
          j=j+1
          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)))
          coordZn(j,2)=rsphere*sin(theta)*sin(phi*dble(i-1))
          coordZn(j,3)=rsphere*cos(theta)
          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)+0.5D0))
          coordO(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)+0.5D0))
          coordO(j,3)=rsphere*cos(theta)
        end do
      case(12)
        nZn=12          ! number of Zn atoms
        nO=12           ! number of O atoms
        dZnO=2.0        ! distance between Zn and O
        ! 1. (uppermost) layer (4 atoms)
        theta=(1.0D0*Pi/6.0D0)
        rsphere=dZnO/(sin(theta)*sqrt(2.0D0))
        phi=(8.0D0*Pi/8.0D0)
        j=0
        do i=1,2
          j=j+1
          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)))
          coordZn(j,2)=rsphere*sin(theta)*sin(phi*dble(i-1))
          coordZn(j,3)=rsphere*cos(theta)
          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)+0.5D0))
          coordO(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)+0.5D0))
          coordO(j,3)=rsphere*cos(theta)
        end do
        ! 2. layer (4 atoms)
        theta=(2.0D0*Pi/6.0D0)
        phi=(8.0D0*Pi/8.0D0)
        do i=1,2
          j=j+1
          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)))
          coordO(j,2)=rsphere*sin(theta)*sin(phi*dble(i-1))
          coordO(j,3)=rsphere*cos(theta)
          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)+0.5D0))
          coordZn(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)+0.5D0))
          coordZn(j,3)=rsphere*cos(theta)
        end do
        ! 3. (middle) layer (8 atoms)
        theta=(4.0D0*Pi/8.0D0)
        phi=(8.0D0*Pi/8.0D0)
        do i=1,2
          j=j+1
          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i)+0.625D0))
          coordO(j,2)=rsphere*sin(theta)*sin(phi*(dble(i)+0.625D0))
          coordO(j,3)=rsphere*cos(theta)
          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)+0.125D0))
          coordZn(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)+0.125D0))
          coordZn(j,3)=rsphere*cos(theta)
        end do
        do i=1,2
          j=j+1
          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i)+0.375D0))
          coordO(j,2)=rsphere*sin(theta)*sin(phi*(dble(i)+0.375D0))
          coordO(j,3)=rsphere*cos(theta)
          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)-0.125D0))
          coordZn(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)-0.125D0))
          coordZn(j,3)=rsphere*cos(theta)
        end do
        ! 4. layer (4 atoms)
        theta=(4.0D0*Pi/6.0D0)
        phi=(8.0D0*Pi/8.0D0)
        do i=1,2
          j=j+1
          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)))
          coordO(j,2)=rsphere*sin(theta)*sin(phi*dble(i-1))
          coordO(j,3)=rsphere*cos(theta)
          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)+0.5D0))
          coordZn(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)+0.5D0))
          coordZn(j,3)=rsphere*cos(theta)
        end do
        ! 5. (lowest) layer (4 atoms)
        theta=(5.0D0*Pi/6.0D0)
        phi=(8.0D0*Pi/8.0D0)
        do i=1,2
          j=j+1
          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)))
          coordZn(j,2)=rsphere*sin(theta)*sin(phi*dble(i-1))
          coordZn(j,3)=rsphere*cos(theta)
          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)+0.5D0))
          coordO(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)+0.5D0))
          coordO(j,3)=rsphere*cos(theta)
        end do
      case(8)
        nZn=8          ! number of Zn atoms
        nO=8           ! number of O atoms
        dZnO=2.0        ! distance between Zn and O
        ! 1. (uppermost) layer (4 atoms)
        theta=(1.75D0*Pi/6.0D0)
        rsphere=dZnO/sin(theta)
        phi=(8.0D0*Pi/8.0D0)
        j=0
        do i=1,2
          j=j+1
          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)))
          coordZn(j,2)=rsphere*sin(theta)*sin(phi*dble(i-1))
          coordZn(j,3)=rsphere*cos(theta)
          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)+0.5D0))
          coordO(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)+0.5D0))
          coordO(j,3)=rsphere*cos(theta)
        end do
c        ! 2. (middle) layer (8 atoms)
c        theta=(4.0D0*Pi/8.0D0)
c        phi=(4.0D0*Pi/8.0D0)
c        do i=1,4
c          j=j+1
c          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)))
c          coordO(j,2)=rsphere*sin(theta)*sin(phi*dble(i-1))
c          coordO(j,3)=rsphere*cos(theta)
c          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)+0.5D0))
c          coordZn(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)+0.5D0))
c          coordZn(j,3)=rsphere*cos(theta)
c        end do
        ! 2. (middle) layer (8 atoms)
        theta=(4.0D0*Pi/8.0D0)
        phi=(8.0D0*Pi/8.0D0)
        do i=1,2
          j=j+1
          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i)+0.625D0))
          coordZn(j,2)=rsphere*sin(theta)*sin(phi*(dble(i)+0.625D0))
          coordZn(j,3)=rsphere*cos(theta)
          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)+0.125D0))
          coordO(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)+0.125D0))
          coordO(j,3)=rsphere*cos(theta)
        end do
        do i=1,2
          j=j+1
          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i)+0.375D0))
          coordZn(j,2)=rsphere*sin(theta)*sin(phi*(dble(i)+0.375D0))
          coordZn(j,3)=rsphere*cos(theta)
          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)-0.125D0))
          coordO(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)-0.125D0))
          coordO(j,3)=rsphere*cos(theta)
        end do
        ! 3th (lowest) layer (4 atoms)
        theta=(4.25D0*Pi/6.0D0)
        phi=(8.0D0*Pi/8.0D0)
        do i=1,2
          j=j+1
          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)))
          coordZn(j,2)=rsphere*sin(theta)*sin(phi*dble(i-1))
          coordZn(j,3)=rsphere*cos(theta)
          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)+0.5D0))
          coordO(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)+0.5D0))
          coordO(j,3)=rsphere*cos(theta)
        end do
      case default 
        print(ferrmssg),"Can deal with 8, 12 or 16 ZnO units so far."
        nerr=nerr+1
        return
      end select
      
      ! write coordinate file
      open(51,file="CageZnO.xyz",status="replace")
      write(51,*) 2*nZnO
      write(51,*) "(ZnO) cage"
      do i=1,nZnO
        !write(51,'("Zn",1x,"core",3(F15.10))')coordZn(i,1:3)
        write(51,'("Zn",1x,3(F15.10))')coordZn(i,1:3)
      end do
      do i=1,nZnO
        !write(51,'("O ",1x,"core",3(F15.10))')coordO(i,1:3)
        !write(51,'("O ",1x,"shel",3(F15.10))')coordO(i,1:3)
        write(51,'("O ",1x,3(F15.10))')coordO(i,1:3)
      end do
      close(51)

      print(fsubend)
      return 

      end subroutine

c---------------------------------------------------------------------

        subroutine nicecluster(infile,informat,clspecies,radius,
     &    coordminpar,coordmeanpar,dipolepar)
        ! cuts a cluster out of bulk such that the average coordination number 
        ! is maximized for a given number of atoms in the cluster
        use defs
        use readcoords
        use writecoords
        use misc
        use cutmod

        implicit none
        character (len=*) infile,informat    ! name and format of input bulk structure file
        type(element), intent(in) :: clspecies(:)  ! a vector with the element names and atom numbers the cluster is supposed to contain. E.g. clspecies=(Cu 2,Se 4)
        double precision radius ! sphere radius,   
        double precision coordminpar,coordmeanpar,dipolepar ! thresholds for cluster acceptance 
        ! internal variables: 
        type(element), allocatable :: species(:),sphspecies(:)   ! the elements in the bulk structure, in the sphere
        type(atom), allocatable :: atoms(:),sphatoms(:),sphatoms2(:)  ! the atoms in the bulk, atoms in the sphere
        type(atom), allocatable :: clatoms(:) ! the atoms in the cluster
        integer :: natoms,nclatoms,nsphatoms  ! number of atoms in the bulk, cluster, in the sphere
        integer :: nspecies,nclspecies,nsphspecies  ! number of species in the bulk, cluster, in the sphere
        integer iatom,iatom2  ! atom index
        integer ispecies,jspecies,i1  ! species index
!        integer(kind=16) n,k,noverk ! number of different configurations, temporarily used integers
!        integer(kind=16), allocatable :: nconfigs(:) ! array with total number of configurations (nconfigs(0)), and numbers of configurations of each species (nconfigs(1:nclspecies))
        integer(8) n,k,noverk ! number of different configurations, temporarily used integers
        integer(8), allocatable :: nconfigs(:) ! array with total number of configurations (nconfigs(0)), and numbers of configurations of each species (nconfigs(1:nclspecies))
        double precision vecs(1:3,1:3)  ! lattice vectors of bulk
        double precision origin(1:3) ! origin of sphere  
        character filename*256
        integer , allocatable :: occ(:)  ! binary array, contains the combination of occupations
        integer, allocatable :: iconfig(:)   ! combination index
        integer j ! another index
        integer imove ! next species for which to change the configuration 
        integer coordmin, coordmax ! minimum, maximum coordination number in the cluster
        double precision coordmean,dipole(1:3),absdipole,ecoul ! average coordination number in the cluster, dipole moment, modulus of dipole moment, electrostatic energy of point charges
        logical talkinloop  ! if (talk), loops that call routines can lead to too much talking. Better call the routines in the loop in silent mode.

        if(talk)print fsubstart, "nicecluster"
        talkinloop=.false.

        ! read bulk structure file and check for
        ! compatibility with desired cluster properties
        call read_coords(infile,informat,atoms,natoms,species,
     &                       nspecies,vecs)
        ! write bulk structure to check
        filename(1:5)="BULK."
        filename(6:len_trim(infile)+5)=trim(infile)
        call write_coords(filename,informat,atoms,natoms,
     &     species,nspecies,vecs)
        ! check if enough species are present in the bulk structure
        nspecies=size(species)
        nclspecies=size(clspecies)
        if(nspecies.ne.nclspecies) goto 1000 
        ! if so, check if they are actually the correct ones
        do ispecies=1,nclspecies
          if(.not.any(species%name(1:2).eq.
     &      clspecies(ispecies)%name(1:2))) 
     &      goto 1000     
        end do   
        
        ! if the provided structure file is that of the bulk (radius>0), 
        ! cut a sphere that should be larger than the expected cluster radius,
        ! otherwise (radius<0 assume the structure file already contains
        ! the sphere
        if (radius.gt.0.0D0) then 
          filename=' '  ! output file for sphere cutting routine
          filename(1:7)="SPHERE."
          filename(8:len_trim(infile)+7)=trim(infile)
          origin(1:3)=0.5d0*(vecs(1,1:3)+vecs(2,1:3)+vecs(3,1:3)) ! place origin of sphere in center of bulk cell
          call cut(infile,filename,informat,informat,                     &
     &             "sphere",0,radius,origin,nndist_def,hdist_def)
          ! Read in the coordinate file of the sphere
          call read_coords(filename,informat,sphatoms,nsphatoms,
     &                     sphspecies,nsphspecies,vecs)
        else
          allocate(sphatoms(natoms),sphspecies(nspecies))
          sphatoms=atoms
          nsphatoms=size(sphatoms)
          sphspecies=species
          nsphspecies=size(sphspecies)
        end if !(radius.gt.0.D0)
        ! write sphere to check
        filename(1:17)="SPHERE.UNORDERED."
        filename(18:len_trim(infile)+17)=trim(infile)
        call write_coords(filename,informat,sphatoms,nsphatoms,
     &     sphspecies,nsphspecies,vecs)

        ! check if enough species are present in the sphere
        nsphspecies=size(sphspecies)
        if(nsphspecies.ne.nclspecies) goto 1001 
        ! if so, check if they are actually the correct ones
        do ispecies=1,nclspecies
          if(.not.any(sphspecies%name(1:2).eq.
     &      clspecies(ispecies)%name(1:2))) 
     &      goto 1001     
        end do   
        
        ! reorder species of sphere to get the same order as in cluster
        allocate(sphatoms2(nsphatoms))
        iatom2=1
        do ispecies=1,nclspecies
          do iatom=1,nsphatoms 
            if(sphatoms(iatom)%name(1:2).eq.
     &         clspecies(ispecies)%name(1:2)) then
              sphatoms2(iatom2)=sphatoms(iatom)
              iatom2=iatom2+1
            end if
          end do
        end do
        call getspecies(sphatoms2,sphspecies)
        sphatoms=sphatoms2
        deallocate(sphatoms2) 
        ! write ordered sphere to check
        filename(1:17)="SPHERE.REORDERED."
        filename(18:len_trim(infile)+17)=trim(infile)
        call write_coords(filename,informat,sphatoms,nsphatoms,
     &     sphspecies,nsphspecies,vecs)
        ! check if enough atoms of each species are present 
        do ispecies=1,nclspecies
          if(sphspecies(ispecies)%howmany.lt.
     &      clspecies(ispecies)%howmany) 
     &      goto 1001     
        end do

        ! assign to all species in the sphere an array of atoms of that
        ! species
        do ispecies=1,nsphspecies
          iatom2=1
          n=sphspecies(ispecies)%howmany
          allocate(sphspecies(ispecies)%atoms(1:n))
          do iatom=1,nsphatoms 
            if(sphatoms(iatom)%name(1:2).eq.
     &        sphspecies(ispecies)%name(1:2)) then
              sphspecies(ispecies)%atoms(iatom2)=sphatoms(iatom)
!              print*,ispecies,sphspecies(ispecies)%name(1:2),iatom2,
!     &           sphspecies(ispecies)%atoms(iatom2)%where
              iatom2=iatom2+1
            end if
          end do
        end do
!        do iatom=1,nsphatoms
!          print*,sphatoms(iatom)%name(1:2),iatom,sphatoms(iatom)%where
!        end do

        ! print the number of possible configurations
        allocate(nconfigs(0:nclspecies),iconfig(0:nclspecies))
        nconfigs(0)=1
        nclatoms=0
        do ispecies=1,nclspecies
          n=sphspecies(ispecies)%howmany
          k=clspecies(ispecies)%howmany
          nclatoms=nclatoms+k  ! get the number of atoms in the cluster
          noverk=binomial(n,k)
          nconfigs(ispecies)=noverk
          FMT1=' '
          FMT2=' '
          FMT3=' '
          FMT4=' '
          WRITE(FMT1,*) k
          FMT1=adjustl(FMT1)
          WRITE(FMT2,*) n
          FMT2=adjustl(FMT2)
          WRITE(FMT3,*) noverk
          FMT3=adjustl(FMT3)
          WRITE(FMT4,*) '(8x,"Number of sites and possible configuration
     &s of ",I',len_trim(FMT1),',1x,A2," over ",I',len_trim(FMT2),'," la
     &ttice sites: ",I',len_trim(FMT3),')'
          if(talk)print FMT4, k,clspecies(ispecies)%name(1:2),
     &      n,noverk 
          nconfigs(0)=nconfigs(0)*noverk
        end do
        ! Total number of configurations:
        WRITE(FMT1,*) nconfigs(0)
        write(FMT2,*) '(8x,"Total number of possible configurations: ",
     &    I',len_trim(adjustl(FMT1)),')'
        if(talk) print FMT2,nconfigs(0) 
 
        !!! test !!!
        ! now for each species distribute the atoms over the available sublattice
        ! sites of that species in the sphere. Do this for each species
        ! separately first as a test. For the real thing, every 
        ! configuration of each element has to be combined with every configuration of every other
        ! element.
        do ispecies=1,nclspecies ! (cluster species)
          ! go through all combinations
          n=sphspecies(ispecies)%howmany                               
          k=clspecies(ispecies)%howmany
          ! first combination: only the lowest numbers are occupied
          if(allocated(occ)) deallocate(occ)
          allocate(occ(1:n))
          occ=0
          occ(1:k)=1 
          ! now go successively upwards from one combination to the next
          iconfig(ispecies)=1  ! combination index, goes up to (n over k)=nconfigs(ispecies)
          do while (iconfig(ispecies).le.nconfigs(ispecies))
            ! if the occupations add up to k, this array corresponds to
            ! a combination, so print it and update the atomic
            ! occupations
            if(abs(sum(occ,1)-k).eq.0) then
              if(talk) print*,iconfig(ispecies),occ
              iconfig(ispecies)=iconfig(ispecies)+1
              sphspecies(ispecies)%atoms%occ=dble(occ)
!              !print occupations to check ... okay.
!              do iatom=1,n
!                print*,sphspecies(ispecies)%name(1:2),
!     &            sphspecies(ispecies)%atoms(iatom)%occ
!              end do
            end if
            ! add 1 to the binary number that is represented by the array
            ! occ
            if(n.gt.k) occ(1)=mod(occ(1)+1,2)
            j=1
            do while (occ(j)==0)
              occ(j+1)=mod(occ(j+1)+1,2)
              j=j+1
            end do
          end do ! iconfig(ispecies)
        end do ! ispecies
        !!! end test !!!
        
        ! now combine all combinations for all species 
        iconfig(0)=1   ! initial index of all configurations (runs up to the total number of configurations)
        iconfig(1:nclspecies)=1   ! initial index of all configurations of a species (runs up to the total number of configurations of that species)
        imove=nclspecies      ! initial species for which to get the next configuration
        do while (nconfigs(imove).eq.1)
!          if(talk) then
!            print*,"imove,iconfig(imove),nconfigs(imove):"
!            print*,imove,iconfig(imove),nconfigs(imove)
!          end if
          iconfig(imove)=1
          imove=imove-1
        end do
!        if(talk) then
!          print*,"imove,iconfig(imove),nconfigs(imove):"
!          print*,imove,iconfig(imove),nconfigs(imove)
!          print*,""
!        end if
        allocate(clatoms(nclatoms))  ! The cluster atoms
        
        ! get the initial occupations (first combination of each
        ! species)
        do ispecies=1,nclspecies ! 
          n=sphspecies(ispecies)%howmany                               
          k=clspecies(ispecies)%howmany
          ! first combination: only the lowest numbers are occupied
          if(allocated(occ)) deallocate(occ)
          allocate(occ(1:n))
          occ=0
          occ(1:k)=1
          ! pass the occupations to the atoms 
          sphspecies(ispecies)%atoms%occ=dble(occ)
        end do
     
        ! open file to which output is written
        open(52,file="CLUSTERS.DAT",status="replace")
        write(52,'("#     Config. number, min, max, avcoord.,   dipmom (
     &eV*Angs),   electrostatic energy (eV)")') 
        if(talk) then
          print  '("#     Config. number, min, max, avcoord.,   dipmom (
     &eV*Angs),   electrostatic energy (eV)")' 
        end if

        ! run over all combinations ! 
        do while (iconfig(0).le.nconfigs(0))  ! loop over all cluster configurations
          ! get the cluster geometry resulting from the current occupation combination
          iatom2=1
          do ispecies=1,nclspecies
            do iatom=1,size(sphspecies(ispecies)%atoms)
              if(sphspecies(ispecies)%atoms(iatom)%occ.gt.0.9D0) then
                clatoms(iatom2)=sphspecies(ispecies)%atoms(iatom)
                iatom2=iatom2+1
              end if
            end do
          end do
          ! determine the coordination numbers of the atoms
          call coord(clatoms,vecs,nndist_def,coordmin,coordmax,
     &      coordmean)
          ! calculate the dipole moment of the cluster
          call dipmom(clatoms,dipole)
          absdipole=absvec(dipole)
          ! calculate electrostatic energy of the cluster
          call ecoulomb(clatoms,ecoul)
           if(talk) then
             print '(I20,I5,I5,F10.6,2(F20.6))', iconfig(0),coordmin,
     &       coordmax,coordmean,absdipole,ecoul
           end if 
          ! write a structure file only if the rewuirements are met
          if(coordmin.ge.coordminpar.and.
     &     coordmean.ge.coordmeanpar.and. 
     &     absdipole.le.dipolepar) then
            write(52,'(I20,I5,I5,F10.6,2(F20.6))') iconfig(0),coordmin,
     &        coordmax,coordmean,absdipole,ecoul
            if(talk) print '(8x,"good cluster")'
            filename='CLUSTER.'
            FMT1=' '
            write(FMT1,*) iconfig(0)
            FMT1=trim(adjustl(FMT1))
            write(FMT2,*) '(I',len_trim(FMT1),')'
            write(filename(9:8+len_trim(FMT1)),FMT2)iconfig(0) !FMT1(1:len_trim(FMT1))
            write(filename(9+len_trim(FMT1):12+len_trim(FMT1)),'(A4)') 
     &         ".xsf"
            call write_coords(filename,"xsf",clatoms,nclatoms,clspecies,
     &           nclspecies,vecs)
          end if 
          ! determine next imove (next species for which to change the
          ! occupation combination)
          ! This is a preparatory step for the next cycle. For the
          ! last cylce it doesn't work but is alos not needed, so leave
          ! it out then:
          if(iconfig(0).lt.nconfigs(0)) then
            imove=nclspecies
            do while (iconfig(imove).eq.nconfigs(imove).or.
     &       nconfigs(imove).eq.1)
              iconfig(imove)=1
              ! set combination for that species to the initial one 
              !if(n.gt.k) then
                if(allocated(occ)) deallocate(occ)
                n=sphspecies(imove)%howmany                            
                k=clspecies(imove)%howmany
                allocate(occ(1:n))
                occ=0
                occ(1:k)=1
                sphspecies(imove)%atoms%occ=dble(occ) ! update the occupations of that species
              !end if
              imove=imove-1
            end do
            ! move the species "imove" to the next configuration
            ! add 1 to the binary number that is represented by the array
            ! occ. 
            n=sphspecies(imove)%howmany                               
            k=clspecies(imove)%howmany
            ! first copy the occupations of that species to an integer
            ! array "occ", just to simplify the writing
            if(allocated(occ)) deallocate(occ)
            allocate(occ(1:n))
            occ(1:n)=nint(sphspecies(imove)%atoms(1:n)%occ)
            if(n.gt.k) occ(1)=mod(occ(1)+1,2)
            j=1
            do while (occ(j)==0)
              occ(j+1)=mod(occ(j+1)+1,2)
              j=j+1
            end do
            ! if the occupations don't add up to k, repeat until they do
            do while(abs(sum(occ,1)-k).gt.0) 
              if(n.gt.k) occ(1)=mod(occ(1)+1,2)
              j=1
              do while (occ(j)==0)
                occ(j+1)=mod(occ(j+1)+1,2)
                j=j+1
              end do
            end do
            ! now the occupations should add up to k, and this array corresponds to
            ! a combination 
            sphspecies(imove)%atoms%occ=dble(occ) ! update the occupations of that species
            iconfig(imove)=iconfig(imove)+1 ! update the configuration number of that species and the total configuration number
          end if
          iconfig(0)=iconfig(0)+1
        end do   
        close(52) 

        ! normal exit
        if(talk) print fsubend
        return

        ! exit with error(s)
1000    nerr=nerr+1
        close(52) 
        print ferrmssg, "nicecluster: Your bulk cell contains too few or
     & the wrong atoms"
        return

1001    nerr=nerr+1
        close(52) 
        print ferrmssg, 'nicecluster: Your sphere is too small. Please i
     &ncrease "radius"'
        return

        end subroutine

c---------------------------------------------------------------------

        subroutine nicecluster2(infile,informat,clspecies,radius0,
     &    radius,origin,coordminpar,coordmeanpar,dipolepar,
     &    test)
        ! cuts a cluster out of bulk such that the average coordination number 
        ! is maximized for a given number of atoms in the cluster
        use defs
        use readcoords
        use writecoords
        use misc
        use replace
        use cutmod

        implicit none
        character (len=*) infile,informat    ! name and format of input bulk structure file
        type(element), intent(in) :: clspecies(:)  ! a vector with the element names and atom numbers the cluster is supposed to contain. E.g. clspecies=(Cu 2,Se 4)
        double precision radius,radius0 ! sphere radius, radius of fixed inner core of cluster 
        double precision origin(1:3) ! origin of sphere  
        double precision coordminpar,coordmeanpar,dipolepar ! thresholds for cluster acceptance 
        logical test ! talkativity, testmode
        ! internal variables: 
        type(element), allocatable :: species(:),sphspecies(:)   ! the elements in the bulk structure, in the sphere
        type(element), allocatable :: corespecies(:),sphbndspecies(:)   ! the elements in the core region, in the sphere boundary
        type(element), allocatable :: clbndspecies(:)   ! the elements in the cluster boundary
        type(atom), allocatable :: atoms(:),sphatoms(:),sphatoms2(:)  ! the atoms in the bulk, atoms in the sphere
        type(atom), allocatable :: clatoms(:) ! the atoms in the cluster
        type(atom), allocatable :: clbndatoms(:),sphbndatoms(:),
     &            coreatoms(:) ! the atoms in the boundary region of the cluster, in the boundary region of the sphere, in the fixed core region
        integer :: natoms,nclatoms,nsphatoms  ! number of atoms in the bulk, cluster, in the sphere
        integer :: ncoreatoms,nclbndatoms,nsphbndatoms  ! number of atoms in the fixed core region, in the boundary of teh cluster, of the sphere
        integer :: nspecies,nclspecies,nsphspecies  ! number of species in the bulk, cluster, in the sphere
        integer :: ncorespecies,nclbndspecies,nsphbndspecies  ! number of species in the bulk, cluster, in the sphere
        integer iatom,iatom2,iatom3  ! atom index
        integer ispecies,jspecies,i1  ! species index
        !integer(kind=16) n,k,noverk ! number of different configurations, temporarily used integers
        !integer(kind=16), allocatable :: nconfigs(:) ! array with total number of configurations (nconfigs(0)), and numbers of configurations of each species (nconfigs(1:nclspecies))
        integer(8) n,k,noverk ! number of different configurations, temporarily used integers
        integer(8), allocatable :: nconfigs(:) ! array with total number of configurations (nconfigs(0)), and numbers of configurations of each species (nconfigs(1:nclspecies))
        double precision vecs(1:3,1:3)  ! lattice vectors of bulk
        character filename*256
        integer , allocatable :: occ(:)  ! binary array, contains the combination of occupations
        integer, allocatable :: iconfig(:)   ! combination index
        integer j ! another index
        integer imove ! next species for which to change the configuration 
        integer coordmin, coordmax ! minimum, maximum coordination number in the cluster
        double precision coordmean,dipole(1:3),absdipole,ecoul ! average coordination number in the cluster, dipole moment, modulus of dipole moment, electrostatic energy of point charges
        logical talkinloop  ! if (talk), loops that call routines can lead to too much talking. Better call the routines in the loop in silent mode.
        logical iscore  ! true if atom belongs to the fixed inner core of the cluster

        if(talk)print fsubstart, "nicecluster2"
        talkinloop=.false.

        ! read bulk structure file and check for
        ! compatibility with desired cluster properties
        call read_coords(infile,informat,atoms,natoms,species,
     &                       nspecies,vecs)
        ! write bulk structure to check
        filename=' '
        filename(1:5)="BULK."
        filename(6:len_trim(infile)+5)=trim(infile)
        call write_coords(filename,informat,atoms,natoms,
     &     species,nspecies,vecs)
        ! check if enough species are present in the bulk structure
        nspecies=size(species)
        nclspecies=size(clspecies)
        if(nspecies.ne.nclspecies) goto 1000 
        ! if so, check if they are actually the correct ones
        do ispecies=1,nclspecies
          if(.not.any(species%name(1:2).eq.
     &      clspecies(ispecies)%name(1:2))) 
     &      goto 1000     
        end do   
        ! determine number of cluster atoms
        nclatoms=0
        do ispecies=1,nclspecies
          nclatoms=nclatoms+clspecies(ispecies)%howmany
        end do
        if(talk)print '(8x,"Your cluster contains ",I6," atoms.")',
     &      nclatoms     

        ! if the provided structure file is that of the bulk (radius>0), 
        ! cut a sphere that should be larger than the expected cluster radius,
        ! otherwise (radius<0) assume the structure file already contains
        ! the sphere
        if (radius.gt.0.0D0) then 
          filename=' '  ! output file for sphere cutting routine
          filename(1:7)="SPHERE."
          filename(8:len_trim(infile)+7)=trim(infile)
          !origin(1:3)=0.5d0*(vecs(1,1:3)+vecs(2,1:3)+vecs(3,1:3)) ! place origin of sphere in center of bulk cell
          call cut(infile,filename,informat,informat,
     &             "sphere",0,radius,origin,nndist_def,hdist_def)
          ! Read in the coordinate file of the sphere
          call read_coords(filename,informat,sphatoms,nsphatoms,
     &                     sphspecies,nsphspecies,vecs)
        else
          allocate(sphatoms(natoms),sphspecies(nspecies))
          sphatoms=atoms
          nsphatoms=size(sphatoms)
          sphspecies=species
          nsphspecies=size(sphspecies)
        end if !(radius.gt.0.D0)
        ! remove sphere atoms with a coordination lt coordmin
        ! copy sphereatoms to dummy varibale sphatoms2
        allocate(sphatoms2(nsphatoms))
        sphatoms2=sphatoms
        sphatoms2%occ=1.0D0
        ! get coordination numbers
        call coord(sphatoms2,vecs,nndist_def,coordmin,coordmax,
     &    coordmean) 
        ! count the number of sufficiently coordinated atoms
        nsphatoms=0
        do iatom=1,size(sphatoms2)
          if (sphatoms2(iatom)%coord.ge.coordminpar) 
     &      nsphatoms=nsphatoms+1
        end do
        ! store the sufficiently coordinated atoms in sphatoms
        deallocate(sphatoms)
        allocate(sphatoms(nsphatoms))
        iatom2=0
        do iatom=1,size(sphatoms2)
          if(sphatoms2(iatom)%coord.ge.coordminpar) then
            iatom2=iatom2+1
            sphatoms(iatom2)=sphatoms2(iatom)
          end if 
        end do
        if(talk) print '(8x,I5," low-coordinated atoms removed.")',
     &           size(sphatoms2)-nsphatoms
        deallocate(sphatoms2)
        ! write sphere to check
        filename(1:17)="SPHERE.UNORDERED."
        filename(18:len_trim(infile)+17)=trim(infile)
        call write_coords(filename,informat,sphatoms,nsphatoms,
     &     sphspecies,nsphspecies,vecs)

        ! check if enough species are present in the sphere
        nsphspecies=size(sphspecies)
        if(nsphspecies.ne.nclspecies) goto 1001 
        ! if so, check if they are actually the correct ones
        do ispecies=1,nclspecies
          if(.not.any(sphspecies%name(1:2).eq.
     &      clspecies(ispecies)%name(1:2))) 
     &      goto 1001     
        end do   
        
        ! reorder species of sphere to get the same order as in cluster
        allocate(sphatoms2(nsphatoms))
        iatom2=1
        do ispecies=1,nclspecies
          do iatom=1,nsphatoms 
            if(sphatoms(iatom)%name(1:2).eq.
     &         clspecies(ispecies)%name(1:2)) then
              sphatoms2(iatom2)=sphatoms(iatom)
              iatom2=iatom2+1
            end if
          end do
        end do
        call getspecies(sphatoms2,sphspecies)
        sphatoms=sphatoms2
        deallocate(sphatoms2) 
        ! write ordered sphere to check
        filename(1:17)="SPHERE.REORDERED."
        filename(18:len_trim(infile)+17)=trim(infile)
        call write_coords(filename,informat,sphatoms,nsphatoms,
     &     sphspecies,nsphspecies,vecs)
        ! check if enough atoms of each species are present 
        print*,"sphspecies,clspecies:",sphspecies%howmany,
     &            clspecies%howmany
        do ispecies=1,nclspecies
          if(sphspecies(ispecies)%howmany.lt.
     &      clspecies(ispecies)%howmany) 
     &      goto 1001     
        end do

        ! determine the atoms and species of the fixed core region of the sphere   
        filename=' '  ! output file for sphere cutting routine
        filename(1:8)="SPHERE0."
        filename(9:len_trim(infile)+8)=trim(infile)
        !origin(1:3)=0.5d0*(vecs(1,1:3)+vecs(2,1:3)+vecs(3,1:3)) ! place origin of sphere in center of bulk cell
        call cut(infile,filename,informat,informat,
     &           "sphere",0,radius0,origin,nndist_def,hdist_def)
        ! Read in the coordinate file of the inner sphere
        call read_coords(filename,informat,coreatoms,ncoreatoms,
     &                   corespecies,ncorespecies,vecs)
        nsphbndatoms=nsphatoms-ncoreatoms
        nclbndatoms=nclatoms-ncoreatoms
        deallocate(corespecies) ! Determine separately such that all cluster elements are contained, possibly with atom number 0
        
        ! determine the atoms and species in the boundary region
        ! and determine for each species how many atoms belong to the core
        ! and the boundary region
        allocate(clbndatoms(nclbndatoms),sphbndatoms(nsphbndatoms))
        iatom3=0 ! sphere boundary atom index
        do iatom=1,nsphatoms
          ! for each atom in the sphere check if it belongs to the core region
          iscore=.false.
          do iatom2=1,ncoreatoms
            if(absvec(sphatoms(iatom)%abswhere
     &        -coreatoms(iatom2)%abswhere).lt.0.01D0)
     &        iscore=.true.
          end do
          if(.not.iscore) then
            iatom3=iatom3+1
            sphbndatoms(iatom3)=sphatoms(iatom)
          end if
        end do
 
        ! get element numbers of sphere boundary
        nsphbndspecies=nsphspecies  
        allocate(sphbndspecies(nsphbndspecies))
        sphbndspecies=sphspecies
        sphbndspecies%howmany=0
        do ispecies=1,nsphbndspecies
          if(allocated(sphbndspecies(ispecies)%atoms)) 
     &      deallocate(sphbndspecies(ispecies)%atoms)
          do iatom=1,nsphbndatoms
            if(sphbndatoms(iatom)%name(1:2).eq.
     &         sphbndspecies(ispecies)%name(1:2))
     &         sphbndspecies(ispecies)%howmany
     &         =sphbndspecies(ispecies)%howmany+1
          end do        
        end do
        
        ! get element numbers of core
        ncorespecies=nsphspecies  
        allocate(corespecies(ncorespecies))
        corespecies=sphspecies
        do ispecies=1,nclspecies
          if(allocated(corespecies(ispecies)%atoms)) 
     &      deallocate(corespecies(ispecies)%atoms)
        end do
        corespecies%howmany=sphspecies%howmany
     &    -sphbndspecies%howmany
 
        ! check if inner core is small enough (does not contain too many
        ! atoms)
        do ispecies=1,nclspecies
          if(corespecies(ispecies)%howmany.gt.
     &      clspecies(ispecies)%howmany) goto 1002
        end do 
        
        ! get element numbers of cluster boundary
        nclbndspecies=nclspecies  
        allocate(clbndspecies(nclbndspecies))
        clbndspecies=clspecies
        do ispecies=1,nclspecies
          if(allocated(clbndspecies(ispecies)%atoms)) 
     &      deallocate(clbndspecies(ispecies)%atoms)
        end do
        clbndspecies%howmany=clspecies%howmany
     &    -corespecies%howmany

        ! assign to all species in the sphere an array of atoms of that
        ! species
        do ispecies=1,nsphspecies
          iatom2=1
          n=sphspecies(ispecies)%howmany
          allocate(sphspecies(ispecies)%atoms(1:n))
          do iatom=1,nsphatoms 
            if(sphatoms(iatom)%name(1:2).eq.
     &        sphspecies(ispecies)%name(1:2)) then
              sphspecies(ispecies)%atoms(iatom2)=sphatoms(iatom)
              iatom2=iatom2+1
            end if
          end do
        end do

        ! assign to all species in the sphere boundary an array of atoms of that
        ! species
        do ispecies=1,nsphbndspecies
          iatom2=1
          n=sphbndspecies(ispecies)%howmany
          allocate(sphbndspecies(ispecies)%atoms(1:n))
          do iatom=1,nsphbndatoms 
            if(sphbndatoms(iatom)%name(1:2).eq.
     &        sphbndspecies(ispecies)%name(1:2)) then
              sphbndspecies(ispecies)%atoms(iatom2)=sphbndatoms(iatom)
              iatom2=iatom2+1
            end if
          end do
        end do
        
        ! assign to all species in the core an array of atoms of that
        ! species
        do ispecies=1,ncorespecies
          iatom2=1
          n=corespecies(ispecies)%howmany
          allocate(corespecies(ispecies)%atoms(1:n))
          do iatom=1,ncoreatoms 
            if(coreatoms(iatom)%name(1:2).eq.
     &        corespecies(ispecies)%name(1:2)) then
              corespecies(ispecies)%atoms(iatom2)=coreatoms(iatom)
              iatom2=iatom2+1
            end if
          end do
        end do

        ! print the number of possible configurations
        allocate(nconfigs(0:nclspecies),iconfig(0:nclspecies))
        nconfigs(0)=1
        nclbndatoms=0
        do ispecies=1,nclspecies
          n=sphbndspecies(ispecies)%howmany
          k=clbndspecies(ispecies)%howmany
          if(n.gt.0.and.k.gt.0) then
            nclbndatoms=nclbndatoms+k  ! get the number of atoms in the core region of cluster
            noverk=binomial(n,k)
            nconfigs(ispecies)=noverk
            FMT1=' '
            FMT2=' '
            FMT3=' '
            FMT4=' '
            WRITE(FMT1,*) k
            FMT1=adjustl(FMT1)
            WRITE(FMT2,*) n
            FMT2=adjustl(FMT2)
            WRITE(FMT3,*) noverk
            FMT3=adjustl(FMT3)
            WRITE(FMT4,*) '(8x,"Number of sites and possible configurati
     &ons of   ",I',len_trim(FMT1),',1x,A2," over ",I',len_trim(FMT2),
     &        '," lattic  e sites: ",I',len_trim(FMT3),')'
            print FMT4, k,clspecies(ispecies)%name(1:2),
     &        n,noverk 
            nconfigs(0)=nconfigs(0)*noverk
          end if
        end do
        ! Total number of configurations:
        WRITE(FMT1,*) nconfigs(0)
        write(FMT2,*) '(8x,"Total number of possible configurations: ",
     &    I',len_trim(adjustl(FMT1)),')'
        print FMT2,nconfigs(0) 

        ! if test mode, stop here
        if(test) return   
     
        !!! test !!!
        ! now for each species distribute the atoms over the available sublattice
        ! sites of that species in the sphere. Do this for each species
        ! separately first as a test. For the real thing, every 
        ! configuration of each element has to be combined with every configuration of every other
        ! element.
        do ispecies=1,nclbndspecies ! (cluster species)
          if(clbndspecies(ispecies)%howmany.gt.0) then
            ! go through all combinations
            n=sphbndspecies(ispecies)%howmany                           
            k=clbndspecies(ispecies)%howmany
            ! first combination: only the lowest numbers are occupied
            if(allocated(occ)) deallocate(occ)
            allocate(occ(1:n))
            occ=0
            occ(1:k)=1 
            ! now go successively upwards from one combination to the next
            iconfig(ispecies)=1  ! combination index, goes up to (n over k)=nconfigs(ispecies)
            do while (iconfig(ispecies).le.nconfigs(ispecies))
              ! if the occupations add up to k, this array corresponds to
              ! a combination, so print it and update the atomic
              ! occupations
              if(abs(sum(occ,1)-k).eq.0) then
                if(talk) print*,iconfig(ispecies),occ
                iconfig(ispecies)=iconfig(ispecies)+1
                sphbndspecies(ispecies)%atoms%occ=dble(occ)
!                !print occupations to check ... okay.
!                do iatom=1,n
!                  print*,sphspecies(ispecies)%name(1:2),
!     &              sphspecies(ispecies)%atoms(iatom)%occ
!                end do
              end if
              ! add 1 to the binary number that is represented by the array
              ! occ
              if(n.gt.k) occ(1)=mod(occ(1)+1,2)
              j=1
              do while (occ(j)==0)
                occ(j+1)=mod(occ(j+1)+1,2)
                j=j+1
              end do
            end do ! iconfig(ispecies)
          end if ! (clbndspecies(ispecies)%howmany.gt.0)
        end do ! ispecies
        !!! end test !!!

!        ! if test mode, stop here (better stop earlier, the loop above can take already quite
!        ! long for moderate numbers of combinations
!        if(test) return   

        ! now combine all combinations for all species 
        iconfig(0)=1   ! initial index of all configurations (runs up to the total number of configurations)
        iconfig(1:nclbndspecies)=1   ! initial index of all configurations of a species (runs up to the total number of configurations of that species)
        imove=nclbndspecies      ! initial species for which to get the next configuration
        do while (nconfigs(imove).eq.1.or.
     &    clbndspecies(imove)%howmany.eq.0)
!          if(talk) then
!            print*,"imove,iconfig(imove),nconfigs(imove):"
!            print*,imove,iconfig(imove),nconfigs(imove)
!          end if
          iconfig(imove)=1
          imove=imove-1
        end do
!        if(talk) then
!          print*,"imove,iconfig(imove),nconfigs(imove):"
!          print*,imove,iconfig(imove),nconfigs(imove)
!          print*,""
!        end if
        allocate(clatoms(nclatoms))  ! The cluster atoms
        !allocate(clbndatoms(nclbndatoms))  ! The cluster boundary atoms
        
        ! get the initial occupations (first combination of each
        ! species)
        do ispecies=1,nclspecies ! 
          n=sphbndspecies(ispecies)%howmany   
          if(n.gt.0) sphbndspecies(ispecies)%atoms%occ=0.0D0           
          k=clbndspecies(ispecies)%howmany
          if(k.gt.0) then
            ! first combination: only the lowest numbers are occupied
            if(allocated(occ)) deallocate(occ)
            allocate(occ(1:n))
            occ=0
            occ(1:k)=1
            ! pass the occupations to the atoms 
            sphbndspecies(ispecies)%atoms%occ=dble(occ)
          end if !(k.gt.0)
        end do
     
        ! open file to which output is written
        open(52,file="CLUSTERS.DAT",status="replace")
        write(52,'("#     Config. number, min, max, avcoord.,   dipmom (
     &eV*Angs),   electrostatic energy (eV)")') 
        if(talk) then
          print  '("#     Config. number, min, max, avcoord.,   dipmom (
     &eV*Angs),   electrostatic energy (eV)")' 
        end if

        ! run over all combinations ! 
        do while (iconfig(0).le.nconfigs(0))  ! loop over all cluster configurations
          ! get the cluster geometry resulting from the current occupation combination
          iatom2=1
          do ispecies=1,nclbndspecies
            if(clbndspecies(ispecies)%howmany.gt.0) then
             do iatom=1,size(sphbndspecies(ispecies)%atoms)
               if(sphbndspecies(ispecies)%atoms(iatom)%occ.gt.0.9D0)then
                 clbndatoms(iatom2)=sphbndspecies(ispecies)%atoms(iatom)
                 iatom2=iatom2+1
              end if
             end do
            end if
          end do
          ! merge core and cluster boundary atoms
          deallocate(clatoms)
          call addatoms(coreatoms,clbndatoms,clatoms)
          ! determine the coordination numbers of the atoms
          call coord(clatoms,vecs,nndist_def,coordmin,coordmax,
     &      coordmean)
          ! calculate the dipole moment of the cluster
          call dipmom(clatoms,dipole)
          absdipole=absvec(dipole)
          ! calculate electrostatic energy of the cluster
          call ecoulomb(clatoms,ecoul)
           if(talk) then
             print '(I20,I5,I5,F10.6,2(F20.6))', iconfig(0),coordmin,
     &       coordmax,coordmean,absdipole,ecoul
           end if 
          ! write a structure file only if the rewuirements are met
          if(coordmin.ge.coordminpar.and.
     &     coordmean.ge.coordmeanpar.and. 
     &     absdipole.le.dipolepar) then
            write(52,'(I20,I5,I5,F10.6,2(F20.6))') iconfig(0),coordmin,
     &        coordmax,coordmean,absdipole,ecoul
            if(talk) print '(8x,"good cluster")'
            filename='CLUSTER.'
            FMT1=' '
            write(FMT1,*) iconfig(0)
            FMT1=trim(adjustl(FMT1))
            write(FMT2,*) '(I',len_trim(FMT1),')'
            write(filename(9:8+len_trim(FMT1)),FMT2)iconfig(0) !FMT1(1:len_trim(FMT1))
            write(filename(9+len_trim(FMT1):12+len_trim(FMT1)),'(A4)') 
     &         ".xsf"
            call write_coords(filename,"xsf",clatoms,nclatoms,clspecies,
     &           nclspecies,vecs)
          end if ! requirements met 
          ! determine next imove (next species for which to change the
          ! occupation combination)
          ! This is a preparatory step for the next cycle. For the
          ! last cylce it doesn't work but is alos not needed, so leave
          ! it out then:
          if(iconfig(0).lt.nconfigs(0)) then
            imove=nclbndspecies
            do while ((iconfig(imove).eq.nconfigs(imove)).or.
     &       (nconfigs(imove).eq.1).or.
     &       (clbndspecies(imove)%howmany.eq.0))
              iconfig(imove)=1
              ! set combination for that species to the initial one 
              if(allocated(occ)) deallocate(occ)
              n=sphbndspecies(imove)%howmany                         
              k=clbndspecies(imove)%howmany
              !if(n.gt.k.and.k.gt.0) then
              if(k.gt.0) then
                allocate(occ(1:n))
                occ=0
                occ(1:k)=1
                sphbndspecies(imove)%atoms%occ=dble(occ) ! update the occupations of that species
              end if !(k.gt.0)
              imove=imove-1
            end do
            ! move the species "imove" to the next configuration
            ! add 1 to the binary number that is represented by the array
            ! occ. 
            n=sphbndspecies(imove)%howmany                              
            k=clbndspecies(imove)%howmany
            ! first copy the occupations of that species to an integer
            ! array "occ", just to simplify the writing
            if(allocated(occ)) deallocate(occ)
            allocate(occ(1:n))
            occ(1:n)=nint(sphbndspecies(imove)%atoms(1:n)%occ)
            if(n.gt.k) occ(1)=mod(occ(1)+1,2)
            j=1
            do while (occ(j)==0)
              occ(j+1)=mod(occ(j+1)+1,2)
              j=j+1
            end do
            ! if the occupations don't add up to k, repeat until they do
            do while(abs(sum(occ,1)-k).gt.0) 
              if(n.gt.k) occ(1)=mod(occ(1)+1,2)
              j=1
              do while (occ(j)==0)
                occ(j+1)=mod(occ(j+1)+1,2)
                j=j+1
              end do
            end do
            ! now the occupations should add up to k, and this array corresponds to
            ! a combination 
            sphbndspecies(imove)%atoms%occ=dble(occ) ! update the occupations of that species
            iconfig(imove)=iconfig(imove)+1 ! update the configuration number of that species and the total configuration number
          end if
          iconfig(0)=iconfig(0)+1
        end do   
        close(52) 

        ! normal exit
        if(talk) print fsubend
        return

        ! exit with error(s)
1000    nerr=nerr+1
        close(52) 
        print ferrmssg, "nicecluster: Your bulk cell contains too few or
     & the wrong atoms"
        return

1001    nerr=nerr+1
        close(52) 
        print ferrmssg, 'nicecluster: Your sphere is too small. Please i
     &ncrease "radius"'
        return

1002    nerr=nerr+1
        close(52) 
        print ferrmssg, 'nicecluster: Your core region is too large. Ple
     &ase decrease "radius0"'
        return

        end subroutine

c---------------------------------------------------------------------

        subroutine Hcover(atoms,moreatoms,hydatoms)
        ! the routine covers the surface of a cluster ("atoms") with H. In order to find the 
        ! bonds to be saturated, it compares with "moreatoms", which
        ! should contain the cluster plus sufficiently
        ! many additional embedding atoms. The routine substitutes all additional
        ! atoms by H, then it removes the redundant H.
        use defs
        use replace
        use misc
        implicit none

        type(atom) atoms(:),moreatoms(:)  ! the cluster atoms, the cluster plus embedding atoms
        type(atom), allocatable :: hydatoms(:)  ! the atoms of the hydrogenated cluster
        ! internal variables ::
        type(atom), allocatable :: hatoms(:), dummyatoms1(:),             &
     &         dummyatoms2(:) ! the H atoms
        type(atom) hatom(1:1) ! a H atom
        type(element), allocatable :: embspecies(:)
        integer natoms,nmoreatoms,nhatoms,nhydatoms,iatom,ispecies,
     &         icoord
        double precision coordmean,distvec(1:3),hsite(1:3),nnsite(1:3)
        double precision hsiteold(1:3)
        integer coordmin,coordmax
        character*2 elename        
        logical firstH

        if(talk) print fsubstart,"Hcover"

        ! determine the embedding atoms (hatoms) in "moreatoms"
        natoms=size(atoms)
        nmoreatoms=size(moreatoms)
        call diffatoms(atoms,moreatoms,(/1.0d0,0.0d0,0.0d0,0.0d0,1.0d0,   &
     &     0.0d0,0.0d0,0.0d0,1.0d0/),dummyatoms1,hatoms,dummyatoms2,      &
     &     0.01D0)
        ! replace the embedding atoms by H
        nhatoms=size(hatoms)
        call getspecies(hatoms,embspecies)
        do ispecies=1,size(embspecies)
          call rpatoms2(hatoms,embspecies(ispecies)%name(1:2),
     &      embspecies(ispecies)%howmany,"H ")
        end do
        deallocate(embspecies)
        ! merge cluster and H atoms
        call addatoms(atoms,hatoms,dummyatoms1)
        ! remove all redundant H and all atoms that are only surrounded 
        ! by H 
        ! mark all redundat H with occ "0"
        ! determine coordination. Throw away H that is only H-coordinated 
        deallocate(hatoms)
        call coordfinite(dummyatoms1,nndist_def,coordmin,coordmax,
     &      coordmean) 
        dummyatoms1%occ=0.0d0
        do iatom=1,size(dummyatoms1) 
          do icoord=1,dummyatoms1(iatom)%coord
            elename(1:2)=dummyatoms1(iatom)%neighbors(icoord)
            if(elename.ne."H ") dummyatoms1(iatom)%occ=1.0d0
          end do
        end do
        ! remove atoms with occupation 0.0 (redundant H)
        call keepoccatoms(dummyatoms1,hydatoms)
        deallocate(dummyatoms1)
        nhydatoms=size(hydatoms) 

        ! shorten the H bonds  
        do iatom=1,nhydatoms  
          if(hydatoms(iatom)%name(1:2).eq."H ") then
            hsiteold=hydatoms(iatom)%abswhere
            firstH=.true.
            ! if atom=H: look for neighbor that is not H
            do icoord=1,hydatoms(iatom)%coord
              elename(1:2)=hydatoms(iatom)%neighbors(icoord)
              if(elename.ne."H ") then
                print*,"H-bond found"
                ! if neighbor.ne.H: move H closer to it
                hsite=hsiteold
                nnsite=hydatoms(iatom)%nnsites(icoord,1:3)
                distvec=hsite-nnsite
                distvec=distvec/absvec(distvec)
                distvec=distvec*hdist_def
                hsite=nnsite+distvec
                ! if it is the first NN, just move the H.
                if(firstH) then
                  hydatoms(iatom)%abswhere=hsite 
                  firstH=.false.
                ! for the next NN (if the H belongs to more than 1 NN), add an H
                else
                  Hatom=hydatoms(iatom)
                  Hatom(1)%abswhere=hsite 
                  call addatoms(hydatoms,Hatom,dummyatoms1)
                  deallocate(hydatoms)
                  allocate(hydatoms(1:size(dummyatoms1)))
                  hydatoms=dummyatoms1
                  deallocate(dummyatoms1)
                end if ! first H
              end if ! (elename.ne."H ")  
            end do ! icoord1
          end if ! (hydatoms(iatom)%name(1:2).eq."H ") 
        end do ! iatom

        if(talk) print fsubendext,"Hcover"
 
        end subroutine

c---------------------------------------------------------------------

        end module
