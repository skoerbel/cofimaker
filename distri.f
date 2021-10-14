      module distri

      implicit none

      contains


      subroutine dist(distinfile,distinform,
     &              axisint,distbin,distbroad)

      use defs
      use readcoords
      use writecoords
      use misc

      implicit none

      integer axisint
      double precision distbin,distbroad
      character*20 distinfile,distinform

      ! local variables
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer :: natoms,nspecies
      integer iatom,nbins,ibin,ispecies
      double precision :: vecs(1:3,1:3)
      double precision ximin,ximax
      double precision, allocatable :: distr(:,:),r(:)
      double precision, allocatable :: broaddistr(:,:)
      double precision, allocatable :: sumbroaddistr(:)
      ! lattice spacings
      integer nmax
      double precision maxlo,maxhi,maxhighest
      !character FMT1*1024,FMT2*1024

      print fsubstart,"dist"
      print '(8x,"Will calculate distribution along x",I1)',axisint

      ! read in  coordinates ...
            print '(" ")'
            print '(8x,"calling read_coords...")'
            call read_coords(distinfile,distinform,atoms,natoms,species,
     &           nspecies,vecs)
            print '(8x,"Coords read in.")'
      
!      ! get minimum and maximum along distribution direction
!      ximax=0.0D0
!      ximin=0.0D0
!      do iatom=1,natoms
!        if (atoms(iatom)%where(axisint).gt.ximax)
!     &     ximax=atoms(iatom)%where(axisint)
!        if (atoms(iatom)%where(axisint).lt.ximin)
!     &     ximin=atoms(iatom)%where(axisint)
!      end do
!      print '(8x,"smallest x",I1," coordinate:",F10.4)',axisint,ximin
!      print '(8x,"largest x",I1," coordinate:",F10.4)',axisint,ximax
      ximin=0.0D0
      ximax=1.0D0

      ! get number of bins
      !nbins=ceiling((ximax-ximin)/distbin)+2
      nbins=ceiling((ximax-ximin)/distbin)+1
!      print '(8x,"number of bins: ",I, "( ",F10.3,")")',nbins,
!     &     (ximax-ximin)/distbin
      WRITE(FMT2,*) nbins
      WRITE(FMT1,*) '(8x,"number of bins: ",I',len(trim(FMT2)),',F10.3)'
      print FMT1,nbins,(ximax-ximin)/distbin   

      allocate(distr(nspecies,nbins),r(nbins),
     &  broaddistr(nspecies,nbins),sumbroaddistr(nbins))
      distr=0.0D0
      broaddistr=0.0D0
      do ibin=1,nbins
        r(ibin)=ximin+dble(ibin-1)*distbin
        do ispecies=1,nspecies
          do iatom=1,natoms
            if(atoms(iatom)%where(axisint).ge.r(ibin)
     &       .and.atoms(iatom)%where(axisint).lt.r(ibin)+distbin
     &       .and.atoms(iatom)%name.eq.species(ispecies)%name)
     &        distr(ispecies,ibin)=distr(ispecies,ibin)+1.0D0
          end do  
        end do
      end do 
      ! broaden distribution, also include periodic images
      do ispecies=1,nspecies
        do ibin=1,nbins
          do iatom=1,natoms
            if(atoms(iatom)%name.eq.species(ispecies)%name)
     &        broaddistr(ispecies,ibin)=broaddistr(ispecies,ibin)
     &          +exp(-((atoms(iatom)%where(axisint)-r(ibin))
     &          /distbroad)**2)/distbroad          
     &          +exp(-((atoms(iatom)%where(axisint)-r(ibin)-1.0D0)
     &          /distbroad)**2)/distbroad
     &          +exp(-((atoms(iatom)%where(axisint)-r(ibin)+1.0D0)
     &          /distbroad)**2)/distbroad
          end do  
        end do
      end do
      broaddistr=broaddistr/sqrt(Pi)

      ! write file with distribution
      open(51,file="DIST.DAT",status='replace',err=1000)
      write(51,'("# distribution along x",I1," with binwidth",F6.3,
     &           " and broadening ",F10.3)')axisint,distbin,
     &      distbroad
      write(51,*) "# x,", species(1:nspecies)%name
      do ibin=1,nbins
        !write(51,*) r(ibin),distr(1:nspecies,ibin)
        write(51,'(F10.6,20(F20.3))') r(ibin),distr(1:nspecies,ibin)
      end do
      close(51)
      ! write file with broadened distr.
      open(51,file="DISTBRDND.DAT",status='replace',err=1000)
      write(51,'("# broadened distribution along x",I1," with binwidth",
     &  F6.3," and broadening ",F10.3)')axisint,distbin,
     &      distbroad
      do ibin=1,nbins
        write(51,'(F10.6,20(F20.3))') r(ibin),
     &     broaddistr(1:nspecies,ibin) !/sqrt(Pi)
      end do
      close(51)

      ! print integral over broadened distr.
      sumbroaddistr=0.0D0
      do ibin=1,nbins
        do ispecies=1,nspecies
          sumbroaddistr(ispecies)=sumbroaddistr(ispecies)
     &       +broaddistr(ispecies,ibin)
        end do
      end do
      sumbroaddistr=sumbroaddistr*distbin !/sqrt(Pi)
      print '(8x," ")'
      do ispecies=1,nspecies
        print '(8x,"integr. broadened dist.: ",A5,1x,F10.5)',
     &    species(ispecies)%name,sumbroaddistr(ispecies)
!        print '(8x,"should be equal to ",I)',species(ispecies)%howmany
      WRITE(FMT2,*) species(ispecies)%howmany
      WRITE(FMT1,*) '(8x,"should be equal to ",I',len(trim(FMT2)),')'
      print FMT1,species(ispecies)%howmany
      end do
      print '(8x,"If not, change broadening and/or binwidth factor.")'

      ! write file with lattice spacings
      open(51,file="LATTSPAC.DAT",status="replace")
      ! loop over species
      do ispecies=1,nspecies
        write(51,'("# Spacings d of ",A," layers (number of layer,",
     &   " fractional position of layer, absolute position of layer,",
     &   " fractional distance to previous layer, absolute distance to",
     &   " previous layer in Angs):")') 
     &     species(ispecies)%name
        ! find maxima of distribution ("layers")  
        nmax=0
        maxhi=-10.0D0
        maxlo=-10.0D0
        if ((broaddistr(ispecies,nbins).ge.broaddistr(ispecies,2)
     &      .and.broaddistr(ispecies,nbins).gt.
     &      broaddistr(ispecies,nbins-1)))then 
            nmax=nmax+1  
            maxlo=r(nbins)-1.0D0
            maxhi=maxlo
            write(51,'(I10,4(F15.6))')
     &        nmax,maxhi,maxhi*absvec(vecs(axisint,1:3)),
     &        (maxhi-maxlo),(maxhi-maxlo)*absvec(vecs(axisint,1:3))
            if(nmax.eq.1)maxhighest=maxhi+1.0D0
        end if
        do ibin=2,nbins-1
          if (broaddistr(ispecies,ibin).gt.broaddistr(ispecies,ibin-1)
     &        .and.broaddistr(ispecies,ibin).ge.
     &        broaddistr(ispecies,ibin+1))then
              nmax=nmax+1     
            maxhi=r(ibin)
!            if(nmax.gt.1)write(51,'(I10,4(F15.6))')
!     &         nmax,maxhi,maxhi*absvec(vecs(axisint,1:3)),maxhi-maxlo,
!     &         (maxhi-maxlo)*absvec(vecs(axisint,1:3))
!            if(nmax.le.1)
            if(nmax.le.1) maxlo=maxhi
            write(51,'(I10,4(F15.6))')
     &         nmax,maxhi,maxhi*absvec(vecs(axisint,1:3)),
     &         maxhi-maxlo,(maxhi-maxlo)*absvec(vecs(axisint,1:3))
            maxlo=maxhi
            if(nmax.eq.1)maxhighest=maxhi+1.0D0
        end if
        end do
        if (broaddistr(ispecies,nbins).ge.broaddistr(ispecies,2)
     &      .and.broaddistr(ispecies,nbins).gt.
     &      broaddistr(ispecies,nbins-1))then
            nmax=nmax+1     
            maxhi=r(ibin)
!            if(nmax.gt.1)write(51,'(I10,4(F15.6))') 
!     &        nmax,maxhi,maxhi*absvec(vecs(axisint,1:3)),maxhi-maxlo,
!     &        (maxhi-maxlo)*absvec(vecs(axisint,1:3))
!            if(nmax.le.1)
            write(51,'(I10,4(F15.6))') nmax,maxhi,
     &        maxhi*absvec(vecs(axisint,1:3)),maxhi-maxlo,
     &        (maxhi-maxlo)*absvec(vecs(axisint,1:3))
            maxlo=maxhi
            if(nmax.eq.1)maxhighest=maxhi+1.0D0
        end if
        if(maxhighest.gt.maxhi) then
            write(51,'(I10,4(F15.6))') nmax,maxhighest,
     &        maxhighest*absvec(vecs(axisint,1:3)),maxhighest-maxhi,
     &        (maxhighest-maxhi)*absvec(vecs(axisint,1:3))
        end if
      end do
      close(51)


!      ! write coordinates
!      call write_coords(file2,format2,atoms,natoms,
!     &                  species,nspecies,vecs,nerr)

      print fsubendext, 'dist'
      return

1000  nerr=nerr+1
      print ferrmssg, "File 'DIST.DAT' could not be written."
      close(51)
      return

      end subroutine

c---------------------------------------------------------------------

      subroutine rdf(infile,informat,rmax,binwidth,broad,vecs)
      use defs
      use readcoords
      !use writecoords
      use misc

      ! global vars
      implicit none
      character(len=*) infile,informat
      double precision rmax,binwidth,broad
      ! local vars
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer :: natoms,nspecies
      integer ir,jr,nr,iat,jat,nuc(1:3),ispecies,jspecies
      integer i,l,m,n
      double precision :: vecs(1:3,1:3),vol,vsphere,distvec(1:3),dist
      double precision ipr1,ipr2,ipr
      double precision angles(1:3),alats(1:3),rdfint,rdfbrdint,dens
      double precision, allocatable :: r(:),rdf0(:),rdfbrd(:) 
      character filename*40
      logical isopen12
      !character FMT1*1024,FMT2*1024,FMT3*1024,FMT4*1024

      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) "rdf"
      else
            print fsubstart,"rdf"
      end if

      if(isopen12) then
            write(12,'("rmax,binwidth,broad:",3(F10.5))')
     &       rmax,binwidth,broad
      else
            print '(8x,"rmax,binwidth,broad:",3(F10.5))',
     &       rmax,binwidth,broad
      end if
      
      ! read in the  structure
      call read_coords(infile,informat,atoms,natoms,species,
     &           nspecies,vecs)

      ! set up radial grid
      nr=int(rmax/binwidth)
      allocate(r(1:nr),rdf0(1:nr),rdfbrd(1:nr))
      do ir=1,nr
            r(ir)=dble(ir)*rmax/dble(nr)
      end do
      ! volume of the sphere with radius rmax (up to which the corr.
      ! fct. and rdf are calculated)
      vsphere=4.0D0*Pi*r(nr)**3/(3.0D0)

      ! extend the unit cell so that rmax lies inside the extended unit
      ! cell for each atom in the original unit cell.
      ! Take into account also tilted unit cells.
      do i=1,3
        alats(i)=sqrt(vecs(i,1)**2+vecs(i,2)**2+vecs(i,3)**2)
      end do
      angles(1)=acos((vecs(2,1)*vecs(3,1)+vecs(2,2)*vecs(3,2)
     &          +vecs(2,3)*vecs(3,3))/(alats(2)*alats(3)))*180.0D0/Pi
      angles(2)=acos((vecs(1,1)*vecs(3,1)+vecs(1,2)*vecs(3,2)
     &          +vecs(1,3)*vecs(3,3))/(alats(1)*alats(3)))*180.0D0/Pi
      angles(3)=acos((vecs(1,1)*vecs(2,1)+vecs(1,2)*vecs(2,2)
     &          +vecs(1,3)*vecs(2,3))/(alats(1)*alats(2)))*180.0D0/Pi
      do i=1,3
!        nuc(i)=1+ceiling(r(nr)*sin(angles(i)*Pi/180.0D0)/
!     &         (alats(i)*sin(angles(1)*Pi/180.0D0)
!     &        *sin(angles(2)*Pi/180.0D0)*sin(angles(3)*Pi/180.0D0)))
        nuc(i)=ceiling(r(nr)*sin(angles(i)*Pi/180.0D0)/
     &         (alats(i)*sin(angles(1)*Pi/180.0D0)
     &        *sin(angles(2)*Pi/180.0D0)*sin(angles(3)*Pi/180.0D0)))
      end do
      ! calculate the volume of the unit cell
      vol=vecs(1,1)*(vecs(2,2)*vecs(3,3)
     &          -vecs(2,3)*vecs(3,2))
     &          +vecs(1,2)*(vecs(2,3)*vecs(3,1)
     &          -vecs(2,1)*vecs(3,3))
     &          +vecs(1,3)*(vecs(2,1)*vecs(3,2)
     &          -vecs(2,2)*vecs(3,1))
      vol=abs(vol)
      if (isopen12) then
!            write(12,'(8x,"Extended supercell dimensions:",
!     &       I,1x,I,1x,I)') nuc
         WRITE(FMT2,*) nuc(1)
         WRITE(FMT3,*) nuc(2)
         WRITE(FMT4,*) nuc(3)
         WRITE(FMT1,*) '(8x,"Extended supercell dimensions: ",I',
     &     len(trim(FMT2)),',1x,I',len(trim(FMT3)),',1x,',
     &     len(trim(FMT4)),')'
         write(12,FMT1) nuc   
      else
!            print '(8x,"Extended supercell dimensions:",
!     &       I,1x,I,1x,I)', nuc
         WRITE(FMT2,*) nuc(1)
         WRITE(FMT3,*) nuc(2)
         WRITE(FMT4,*) nuc(3)
         WRITE(FMT1,*) '(8x,"Extended supercell dimensions: ",I',
     &     len(trim(FMT2)),
     &     ',1x,I',len(trim(FMT3)),',1x,I',len(trim(FMT4)),')'
         print FMT1, nuc   
      end if
      
      ! rdf for all atoms
      if (nspecies.gt.1) then
        rdf0=0.0D0
        rdfbrd=0.0D0
        do iat=1,natoms
          do jat=1,natoms
            if(jat.ne.iat) then
              do l=-nuc(1),nuc(1)
                do m=-nuc(2),nuc(2)
                  do n=-nuc(3),nuc(3)
                    ! calculate distance of first and second atom
                    ! distance in fractional coordinates
                    distvec(1)=atoms(jat)%where(1)+dble(l)
     &                         -atoms(iat)%where(1)
                    distvec(2)=atoms(jat)%where(2)+dble(m)
     &                         -atoms(iat)%where(2)
                    distvec(3)=atoms(jat)%where(3)+dble(n)
     &                         -atoms(iat)%where(3)
                    ! distance in absolute coordinates
                    distvec(1:3)=distvec(1)*vecs(1,1:3)
     &                      +distvec(2)*vecs(2,1:3)
     &                      +distvec(3)*vecs(3,1:3)       
                    dist=absvec(distvec)
                    if (dist.le.r(1).and.dist.gt.0.0D0) then
                        rdf0(1)=rdf0(1)+1.0D0/r(ir)**2   
                    end if
                    do ir=2,nr
                      if (dist.le.r(ir).and.dist.gt.r(ir-1)) then
                        rdf0(ir)=rdf0(ir)+1.0D0/r(ir)**2    
                      end if
                    end do
                  end do
                end do
              end do
            end if
          end do
        end do
        ! broaden rdf
        do ir=1,nr
          do jr=1,nr
            rdfbrd(ir)=rdfbrd(ir)+rdf0(jr)*r(jr)**2
     &        *exp(-((r(jr)-r(ir))/broad)**2)/((0.50D0*broad**3
     &         +r(jr)**2*broad)*sqrt(Pi)+r(jr)*broad**2) 
          end do
        end do
        ! normalize rdf
        rdf0(1:nr)=rdf0(1:nr)*vol/(dble(natoms**2)
     &       *4.0D0*Pi*(r(2)-r(1)))
        rdfbrd(1:nr)=rdfbrd(1:nr)*vol/(dble(natoms**2)
     &       *4.0D0*Pi) 
        
        ! calculate the inverse participation ratio (IPR) of the brd. rdf
        ipr1=0.0D0
        ipr2=0.0D0
        do ir=1,nr
          ipr1=ipr1+r(ir)**2*rdfbrd(ir)
          ipr2=ipr2+r(ir)**2*rdfbrd(ir)**2
        end do
        ipr=ipr1**2*(r(2)-r(1))*4.0D0*Pi/(ipr2*vsphere)
        ! write rdf and integrated for all atoms
        ! w/o, with broadening
        rdfint=0.0D0
        rdfbrdint=0.0D0
        open(51,file="RDF.ALL.DAT",status='replace',err=1000)
        write(51,'("# IPR of broadened RDF: ",F10.6)')ipr
        write(51,'("# r, rdf, int. rdf, broadened rdf, int. brdnd. rdf,
     &   normalized atom dens.")')
        do ir=1,nr-1
          !calculate integral of rdf
          rdfint=rdfint+rdf0(ir)*r(ir)**2*4.0D0*Pi*binwidth
          rdfbrdint=rdfbrdint+rdfbrd(ir)*r(ir)**2*4.0D0*Pi*binwidth
          dens=rdfbrdint/(4.0D0/3.0D0*Pi*r(ir)**3)
          write(51,2000) r(ir),rdf0(ir),rdfint,rdfbrd(ir),
     &         rdfbrdint,dens   
        end do
        close(51)
      end if

      ! rdf for each species
      do ispecies=1,nspecies
        rdf0=0.0D0
        rdfbrd=0.0D0
        do iat=1,natoms
          do jat=1,natoms
            if(jat.ne.iat.and.atoms(iat)%name(1:2)
     &        .eq.species(ispecies)%name(1:2).and.
     &        atoms(jat)%name(1:2).eq.species(ispecies)%name(1:2))then
              do l=-nuc(1),nuc(1)
                do m=-nuc(2),nuc(2)
                  do n=-nuc(3),nuc(3)
                    ! calculate distance of first and second atom
                    ! distance in fractional coordinates
                    distvec(1)=atoms(jat)%where(1)+dble(l)
     &                         -atoms(iat)%where(1)
                    distvec(2)=atoms(jat)%where(2)+dble(m)
     &                         -atoms(iat)%where(2)
                    distvec(3)=atoms(jat)%where(3)+dble(n)
     &                         -atoms(iat)%where(3)
                    ! distance in absolute coordinates
                    distvec(1:3)=distvec(1)*vecs(1,1:3)
     &                      +distvec(2)*vecs(2,1:3)
     &                      +distvec(3)*vecs(3,1:3)       
                    dist=absvec(distvec)
                    if (dist.le.r(1).and.dist.gt.0.0D0) then
                        rdf0(1)=rdf0(1)+1.0D0/r(ir)**2   
                    end if
                    do ir=2,nr
                      if (dist.le.r(ir).and.dist.gt.r(ir-1)) then
                        rdf0(ir)=rdf0(ir)+1.0D0/r(ir)**2    
                      end if
                    end do
                  end do
                end do
              end do
            end if
          end do
        end do
        ! broaden rdf
        do ir=1,nr
          do jr=1,nr
            rdfbrd(ir)=rdfbrd(ir)+rdf0(jr)*r(jr)**2
     &      *exp(-((r(jr)-r(ir))/broad)**2)
     &      /((0.50D0*broad**3+r(jr)**2*broad)*sqrt(Pi)+r(jr)*broad**2) 
          end do
        end do
        ! normalize rdf
        rdf0(1:nr)=rdf0(1:nr)*vol/
     &     (dble(species(ispecies)%howmany**2)*4.0D0*Pi*(r(2)-r(1)))
        rdfbrd(1:nr)=rdfbrd(1:nr)*vol/
     &     (dble(species(ispecies)%howmany**2)
     &       *4.0D0*Pi) 
        
        ! calculate the inverse participation ratio (IPR) of the brd. rdf
        ipr1=0.0D0
        ipr2=0.0D0
        do ir=1,nr
          ipr1=ipr1+r(ir)**2*rdfbrd(ir)
          ipr2=ipr2+r(ir)**2*rdfbrd(ir)**2
        end do
        ipr=ipr1**2*(r(2)-r(1))*4.0D0*Pi/(ipr2*vsphere)
        ! write rdf and integrated for all atoms
        ! w/o, with broadening
        rdfint=0.0D0
        rdfbrdint=0.0D0
        filename=" "
        write(filename,'("RDF.",A,".DAT")')trim(species(ispecies)%name)
        open(51,file=filename,status='replace',err=1000)
        write(51,'("# IPR of broadened RDF: ",F10.6)')ipr
        write(51,'("# r, rdf, int. rdf, broadened rdf, int. brdnd. rdf,
     &   normalized atom dens. for ",A,1x,A)') 
     &    trim(species(ispecies)%name), trim(species(ispecies)%name)
        do ir=1,nr-1
          !calculate integral of rdf
          rdfint=rdfint+rdf0(ir)*r(ir)**2*4.0D0*Pi*binwidth
          rdfbrdint=rdfbrdint+rdfbrd(ir)*r(ir)**2*4.0D0*Pi*binwidth
          dens=rdfbrdint/(4.0D0/3.0D0*Pi*r(ir)**3)
          write(51,2000) r(ir),rdf0(ir),rdfint,rdfbrd(ir),
     &       rdfbrdint,dens   
        end do
        close(51)
        if (nspecies==1) then
          rdfint=0.0D0
          rdfbrdint=0.0D0
          filename=" "
          write(filename,'("RDF.ALL.DAT")')
          open(51,file=filename,status='replace',err=1000)
          write(51,'("# IPR of broadened RDF: ",F10.6)')ipr
          write(51,'("# r, rdf, int. rdf, broadened rdf, int. brdnd. rdf
     &, normalized atom dens.")') 
          do ir=1,nr-1
            !calculate integral of rdf
            rdfint=rdfint+rdf0(ir)*r(ir)**2*4.0D0*Pi*binwidth
            rdfbrdint=rdfbrdint+rdfbrd(ir)*r(ir)**2*4.0D0*Pi*binwidth
            dens=rdfbrdint/(4.0D0/3.0D0*Pi*r(ir)**3)
            write(51,2000) r(ir),rdf0(ir),rdfint,rdfbrd(ir),
     &         rdfbrdint,dens   
          end do
          close(51)
        end if
      end do

      ! rdf for each species pair
      do ispecies=1,nspecies-1
       do jspecies=ispecies+1,nspecies
        rdf0=0.0D0
        rdfbrd=0.0D0
        do iat=1,natoms
          do jat=1,natoms
            if(jat.ne.iat.and.atoms(iat)%name(1:2)
     &        .eq.species(ispecies)%name(1:2).and.
     &        atoms(jat)%name(1:2).eq.species(jspecies)%name(1:2))then
              do l=-nuc(1),nuc(1)
                do m=-nuc(2),nuc(2)
                  do n=-nuc(3),nuc(3)
                    ! calculate distance of first and second atom
                    ! distance in fractional coordinates
                    distvec(1)=atoms(jat)%where(1)+dble(l)
     &                         -atoms(iat)%where(1)
                    distvec(2)=atoms(jat)%where(2)+dble(m)
     &                         -atoms(iat)%where(2)
                    distvec(3)=atoms(jat)%where(3)+dble(n)
     &                         -atoms(iat)%where(3)
                    ! distance in absolute coordinates
                    distvec(1:3)=distvec(1)*vecs(1,1:3)
     &                      +distvec(2)*vecs(2,1:3)
     &                      +distvec(3)*vecs(3,1:3)       
                    dist=absvec(distvec)
                    if (dist.le.r(1).and.dist.gt.0.0D0) then
                        rdf0(1)=rdf0(1)+1.0D0/r(ir)**2   
                    end if
                    do ir=2,nr
                      if (dist.le.r(ir).and.dist.gt.r(ir-1)) then
                        rdf0(ir)=rdf0(ir)+1.0D0/r(ir)**2    
                      end if
                    end do
                  end do
                end do
              end do
            end if
          end do
        end do
        ! broaden rdf
        do ir=1,nr
          do jr=1,nr
            rdfbrd(ir)=rdfbrd(ir)+rdf0(jr)*r(jr)**2
     &      *exp(-((r(jr)-r(ir))/broad)**2)
     &      /((0.50D0*broad**3+r(jr)**2*broad)*sqrt(Pi)+r(jr)*broad**2)          
          end do
        end do
        ! normalize rdf
        rdf0(1:nr)=rdf0(1:nr)*vol/
     &     (dble(species(ispecies)%howmany)
     &     *dble(species(jspecies)%howmany)*4.0D0*Pi*(r(2)-r(1)))
        rdfbrd(1:nr)=rdfbrd(1:nr)*vol/
     &     (dble(species(ispecies)%howmany)
     &     *dble(species(jspecies)%howmany)*4.0D0*Pi) 
        
        ! calculate the inverse participation ratio (IPR) of the brd. rdf
        ipr1=0.0D0
        ipr2=0.0D0
        do ir=1,nr
          ipr1=ipr1+r(ir)**2*rdfbrd(ir)
          ipr2=ipr2+r(ir)**2*rdfbrd(ir)**2
        end do
        ipr=ipr1**2*(r(2)-r(1))*4.0D0*Pi/(ipr2*vsphere)
        ! write rdf and integrated for all atoms
        ! w/o, with broadening
        rdfint=0.0D0
        rdfbrdint=0.0D0
        filename=" "
        write(filename,'("RDF.",A,A".DAT")')
     &     trim(species(ispecies)%name),trim(species(jspecies)%name)
        open(51,file=filename,status='replace',err=1000)
        write(51,'("# IPR of broadened RDF: ",F10.6)')ipr
        write(51,'("# r, rdf, int. rdf, broadened rdf, int. brdnd. rdf,
     &   normalized atom dens. for ",A,1x,A)') 
     &    trim(species(ispecies)%name), trim(species(jspecies)%name)
        do ir=1,nr-1
          !calculate integral of rdf
          rdfint=rdfint+rdf0(ir)*r(ir)**2*4.0D0*Pi*binwidth
          rdfbrdint=rdfbrdint+rdfbrd(ir)*r(ir)**2*4.0D0*Pi*binwidth
          dens=rdfbrdint/(4.0D0/3.0D0*Pi*r(ir)**3)
          write(51,2000) r(ir),rdf0(ir),rdfint,rdfbrd(ir),
     &       rdfbrdint,dens   
        end do
        close(51)
       end do
      end do

 2000 format(F15.8,4(F20.8),F15.8)

      if(isopen12) then
            write(12,fsubendext) "rdf"
      else
            print fsubendext, 'rdf'
      end if
      return

1000  nerr=nerr+1
      if(isopen12) then
            write(12,ferrmssg) "rdf"
      else
            print ferrmssg,"rdf"
      end if
      return
      end subroutine

c---------------------------------------------------------------------

      subroutine propcorr(infile,ele,iprop,corrf, 
     &                rcut,add,binwidth,rmax,brd)
      ! determines if atoms with the same property cluster in
      ! space
      use defs
      use readcoords
      use misc
      implicit none
      double precision, intent(in) :: brd,add
      double precision vecs(1:3,1:3)
      double precision rcut,propav,propprodlocav,binwidth
      double precision, allocatable ::
     & checkprops(:),propsloc(:)
      character(len=*) infile,corrf
      character ele*2,line*100,elestring*100
      type(atom), allocatable :: atoms(:),checkatoms(:)
      type(element),allocatable :: species(:)
      integer natoms, nspecies,iele,i,j,iprop,nprops,ieleprop
      integer ncheckatoms,nproploc
      ! for correlation function:
      integer :: nr
      double precision rmax,corrint,maxint
      double precision, allocatable :: gcorr(:),r(:),rdistf(:)
      !parameter(rmax=20.0D0)
      logical isopen12

      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart)trim(adjustl('propcorr'))
      else
            print fsubstart,trim(adjustl('propcorr'))
      end if
      
      nr=int(rmax/binwidth)
      allocate(gcorr(1:nr),r(1:nr),rdistf(nr))
      if (isopen12) then
        write(12,'(8x,"File to read is ",A,".")')trim(adjustl(infile))
        write(12,'(8x,"Will check ",A2," for correlations")') 
     &  trim(adjustl(ele))
      else
        print '(8x,"File to read is ",A,".")',trim(adjustl(infile))
        print '(8x,"Will check ",A2," for correlations")', 
     &  trim(adjustl(ele))
      end if
      select case(corrf)
      case('theta') 
            write(12,'(8x,"using ",A," as corr. function")')
     &          trim(adjustl(corrf)) 
            write(12,'(8x,"with cutoff",F10.3," Angs.")') rcut
            write(12,'(8x,"with binwidth",F10.3," Angs.")') binwidth
      case default
            write(12,'(8x,"Error: Unknown corr. funct.")')
            print*,"Error: Unknown corr. funct."
            nerr=nerr+1
            return
      end select      
      ! read number of atoms from cfg file
      call read_cfg(infile,atoms,natoms,species,nspecies,
     &                    vecs)
      nprops=size(atoms(1)%properties)

      do i=1,nspecies
          if(isopen12) then
            write(12,'(8x,A2,I10)')species(i)%name,species(i)%howmany
          end if
            if (species(i)%name(1:2).eq.ele(1:2))ieleprop=i
      end do
      i=1
      if(isopen12) then
            write(12,'(8x,"Coord. and props. of first and last atom",
     &      " of each species:")')
      end if
      do iele=1,nspecies-1
            do j=1,species(iele)%howmany
                  if (j.eq.1.or.j.eq.species(iele)%howmany) then
                    if(isopen12) then
                        write(12,'(8x,A2,3(F15.6),F15.6)') 
     &                  atoms(i)%name(1:2),atoms(i)%where,
     &                  atoms(i)%properties(1:nprops)           
                    end if
                  end if
                  i=i+1
            end do
      end do
      do j=1,species(nspecies)%howmany
            if (j.eq.1.or.j.eq.species(nspecies)%howmany) then
              if(isopen12) then
                  write(12,'(8x,A2,3(F15.6),F15.6)') 
     &            atoms(i)%name(1:2),atoms(i)%where,
     &            atoms(i)%properties(1:nprops)     
              end if
            end if
            i=i+1
      end do
      ! make list of atoms to check
      ncheckatoms=species(ieleprop)%howmany
      allocate(checkatoms(ncheckatoms),checkprops(ncheckatoms),
     &          propsloc(ncheckatoms))
      j=1
      do i=1,natoms
            if(atoms(i)%name(1:2).eq.ele(1:2)) then
                  checkatoms(j)=atoms(i)
                  checkprops(j)=atoms(i)%properties(iprop+1)+add
                  j=j+1
            end if
      end do
      propav=0.0D0
      do i=1,ncheckatoms
            propav=propav+checkprops(i)
      end do
      propav=propav/dble(ncheckatoms)
      ! calculate locally av. prop. and write to file
      open(66,file='PROPCORR.OUT',status='replace')
      write(66,'("# Atom no., prop., globally av.d prop,",
     &  "locally av. prop., # of atoms in local av., <p(0)p(r)>,",
     &  " r<rsph")')
!      j=1
      propprodlocav=0.0D0
      do i=1,ncheckatoms
            select case(corrf)
            case('theta') 
                propsloc(i)=
     &          propavloctheta(checkatoms(i)%where,ncheckatoms,
     &          checkatoms,checkprops,rcut,vecs,nproploc)
                propprodlocav=propprodlocav+propsloc(i)*checkprops(i)
            case default
                propsloc(i)=0.0D0
                nproploc=0
            end select      
      end do
      propprodlocav=propprodlocav/dble(ncheckatoms)
      do i=1,ncheckatoms
            write(66,1000) i,checkprops(i),propav,
     &          propsloc(i),nproploc,propprodlocav
      end do
! 1000 format(I5,3(F10.5),I3)
 1000 format(I5,3(F10.5),I3,F10.5)
      close(66)
      ! calculate correlation function
      do i=1,nr
            r(i)=dble(i)*rmax/dble(nr)
      end do
      call corr(nr,r,gcorr,ncheckatoms,checkatoms,checkprops,vecs,
     &          brd,corrint,maxint,rdistf)
      open(66,file='CORRF.DAT',status='replace')
      write(66,'("# distance, corr. funct., radial distr. f.")')
      do i=1,nr
            write(66,'(3(F15.8))') r(i),gcorr(i),rdistf(i)
      end do
      close(66)
      if(isopen12) then
        write(12,'(8x,"Average of prop.",I3,":",F15.7)') iprop,propav
        write(12,'(8x,"(<p(r)> av. over sim. box)")')
        write(12,'(8x,"<p(0)p(r)> av. over sphere:",F20.5)')
     &   propprodlocav
        write(12,'(8x,"Integrated corr. function corrint:",F20.5)')
     &   corrint
        write(12,'(8x,"Integrated volume maxint:",F20.5)')maxint
        write(12,'(8x,"corrint/maxint:",F10.5)')corrint/maxint
      else
        print'(8x,"Average of prop.",I3,":",F15.7)', iprop,propav
        print'(8x,"(<p(r)> av. over sim. box)")'
        print'(8x,"<p(0)p(r)> av. over sphere:",F20.5)',
     &   propprodlocav
        print'(8x,"Integrated corr. function corrint:",F20.5)',
     &   corrint
        print'(8x,"Integrated volume maxint:",F20.5)',maxint
        print'(8x,"corrint/maxint:",F10.5)',corrint/maxint
      end if
      if(isopen12) then
        write(12,fsubend)
      else
        print fsubend
      end if
      
      return

      end subroutine
                                                         
c---------------------------------------------------------------------

      subroutine corr(nr,r,gcorr,ncheckatoms,checkatoms,checkprops,
     &                vecs,brd,corrint,maxint,rdistf)
      ! calculates correlation function and rdf
      use defs
      implicit none

      double precision, intent(in) :: r(nr),vecs(1:3,1:3)
      double precision gcorr(nr),gcorrbrd(nr),rdistf(nr),rdfbrd(nr)
      double precision brd
      type (atom), intent(in) :: checkatoms(ncheckatoms)
      double precision, intent(in) :: checkprops(ncheckatoms)
      integer, intent(in) :: nr, ncheckatoms
      integer i,j,k,l,m,n,nsample(1:nr),nuc(1:3)
      double precision distvec(1:3),dist,gav,alats(1:3),corrint,
     & maxint,vol,corrintpart,corrintpartbrd,maxintpart,
     & corrintpartfrac,corrintpartfracbrd,maxintpartbrd,ipr1,ipr2,ipr,
     & vsphere 
      logical isopen12
      !character FMT1*1024,FMT2*1024,FMT3*1024,FMT4*1024

      inquire(unit=12,opened=isopen12)
      if(isopen12) then
            write(12,fsubstart)trim(adjustl('corr'))
      else
            print fsubstart,trim(adjustl('corr'))
      end if

      ! take into account also atoms in neighboring unit cells
      do i=1,3
        alats(i)=sqrt(vecs(i,1)**2+vecs(i,2)**2+vecs(i,3)**2)
        nuc(i)=ceiling(r(nr)/alats(i))
      end do
      if(isopen12) then
!            write(12,'(8x,"Extended supercell dimensions:",
!     &       I,1x,I,1x,I)') nuc
         WRITE(FMT2,*) nuc(1)
         WRITE(FMT3,*) nuc(2)
         WRITE(FMT4,*) nuc(3)
         WRITE(FMT1,*) '(8x,"Extended supercell dimensions: ",I',
     &     len(trim(FMT2)),
     &     ',1x,I',len(trim(FMT3)),',1x,I',len(trim(FMT4)),')'
         write(12,FMT1) nuc   
      else
            print*, "Extended supercell dimensions:", nuc
!            print (8x,"Extended supercell dimensions:",
!     &       I,1x,I,1x,I), nuc
      end if

      ! volume of the sphere with radius rmax (up to which the corr.
      ! fct. and rdf are calculated)
      vsphere=4.0D0*Pi*r(nr)**3/(3.0D0)
      !print*,'nuc=',nuc
      nsample=0
      gcorr=0.0D0
      gcorrbrd=0.0D0
      rdistf=0.0D0
      rdfbrd=0.0D0
      !brd=0.2D0
      gav=0.0D0
      do i=1,ncheckatoms
        gav=gav+checkprops(i)
        do j=1,ncheckatoms
         if(j.ne.i) then
          do l=-nuc(1),nuc(1)
            do m=-nuc(2),nuc(2)
              do n=-nuc(3),nuc(3)
                ! distance in fractional coordinates
                distvec(1)=checkatoms(j)%where(1)+dble(l)
     &                     -checkatoms(i)%where(1)
                distvec(2)=checkatoms(j)%where(2)+dble(m)
     &                     -checkatoms(i)%where(2)
                distvec(3)=checkatoms(j)%where(3)+dble(n)
     &                     -checkatoms(i)%where(3)
                ! distance in absolute coordinates
                distvec(1:3)=distvec(1)*vecs(1,1:3)
     &                  +distvec(2)*vecs(2,1:3)
     &                  +distvec(3)*vecs(3,1:3)       
                dist=sqrt(distvec(1)**2+distvec(2)**2
     &             +distvec(3)**2)
                ! corr. function
                if (dist.le.r(1)) then
                  gcorr(1)=gcorr(1)+checkprops(i)*checkprops(j)
                  rdistf(1)=rdistf(1)+1.0D0/r(1)**2
                  rdfbrd(1)=rdfbrd(1)+1.0D0/r(1)**2
                  nsample(1)=nsample(1)+1
                end if
                do k=2,nr
                  if (dist.gt.r(k-1).and.dist.le.r(k)) then
                    gcorr(k)=gcorr(k)+checkprops(i)*checkprops(j)
                    rdistf(k)=rdistf(k)+1.0D0/r(k)**2
                    nsample(k)=nsample(k)+1
                  end if
!                  rdfbrd(k)=rdfbrd(k)
!     &              +exp(-((dist-r(k))/brd)**2)
!     &              /(brd*sqrt(Pi)*r(k)**2)          
                end do
              end do  
            end do
          end do
         end if 
        end do 
      end do
      gav=gav/dble(ncheckatoms)
      ! subtract average from coorelation function and broaden rdf
      do k=1,nr
        if (nsample(k).gt.0) gcorr(k)=gcorr(k)/dble(nsample(k))-gav**2
        ! broaden rdf 
        do j=1,nr
!          rdfbrd(k)=rdfbrd(k)+rdistf(j)*r(j)**2
!     &       *exp(-((r(j)-r(k))/brd)**2)/(brd*sqrt(Pi)*r(k)**2)          
          rdfbrd(k)=rdfbrd(k)+rdistf(j)*r(j)**2
     &      *exp(-((r(j)-r(k))/brd)**2)
     &      /((0.50D0*brd**3+r(j)**2*brd)*sqrt(Pi)+r(j)*brd**2)          
        end do
      end do
      ! broaden corr. function 
      do k=1,nr
        do j=1,nr
!          gcorrbrd(k)=gcorrbrd(k)+gcorr(j)
!     &       *exp(-((r(j)-r(k))/brd)**2)/(brd*sqrt(Pi))          
          gcorrbrd(k)=gcorrbrd(k)+gcorr(j)
     &       *exp(-((r(j)-r(k))/brd)**2)*r(j)**2*(r(2)-r(1))
     &       /((0.50D0*brd**3+r(j)**2*brd)*sqrt(Pi)+r(j)*brd**2)          
        end do
      end do
      ! volume
      vol=vecs(1,1)*(vecs(2,2)*vecs(3,3)
     &          -vecs(2,3)*vecs(3,2))
     &          +vecs(1,2)*(vecs(2,3)*vecs(3,1)
     &          -vecs(2,1)*vecs(3,3))
     &          +vecs(1,3)*(vecs(2,1)*vecs(3,2)
     &          -vecs(2,2)*vecs(3,1))
      vol=abs(vol)
      rdistf=rdistf*vol/(dble(ncheckatoms**2)*4.0D0*Pi*(r(2)-r(1)))
      rdfbrd=rdfbrd*vol/(dble(ncheckatoms**2)*4.0D0*Pi)
      !rdistf=rdistf*vsphere/(dble(ncheckatoms**2)*4.0D0*Pi*(r(2)-r(1)))
      !rdfbrd=rdfbrd*vsphere/(dble(ncheckatoms**2)*4.0D0*Pi)
      !integrate correlation function
      corrint=0.0D0
      maxint=0.0D0
      do k=1,nr
        corrint=corrint+r(k)**2*gcorr(k)
        maxint=maxint+r(k)**2
      end do
      corrint=corrint*4.0D0*Pi*(r(2)-r(1))
      maxint=maxint*4.0D0*Pi*(r(2)-r(1))

      ! calculate the inverse participation ratio (IPR) of the rdf
      ipr1=0.0D0
      ipr2=0.0D0
      do k=1,nr
        !ipr1=ipr1+r(k)**2*1.0D0
        !ipr2=ipr2+r(k)**2*1.0D0
        ipr1=ipr1+r(k)**2*rdfbrd(k)
        ipr2=ipr2+r(k)**2*rdfbrd(k)**2
      end do
      ipr=ipr1**2*(r(2)-r(1))*4.0D0*Pi/(ipr2*vsphere)
      ! print broadened rdf
      open(65,file="RDF.DAT",status="replace")
      write(65,'("# IPR of broadened RDF: ",F10.6)')ipr
      write(65,'("# r, rdf, broadened rdf")')
      do k=1,nr-1
        write(65,'(3(F15.8))') r(k),rdistf(k),rdfbrd(k)
      end do
      close(65)
      ! print broadened corr
      corrintpart=0.0D0
      corrintpartbrd=0.0D0
      maxintpart=0.0D0
      maxintpartbrd=0.0D0
      open(65,file="CORR.DAT",status="replace")
      write(65,'("# r, corr, broadened corr, integrated corr., integrate
     &d broadened corr., relative integrals")')
      do k=1,nr-1
        corrintpart=corrintpart+r(k)**2*gcorr(k)*4.0D0*Pi*(r(2)-r(1))
        corrintpartbrd=corrintpartbrd+r(k)**2*gcorrbrd(k)
     &       *4.0D0*Pi*(r(2)-r(1))
        maxintpart=maxintpart+r(k)**2*rdistf(k)
     &       *4.0D0*Pi*(r(2)-r(1))
        maxintpartbrd=maxintpartbrd+r(k)**2*rdfbrd(k)
     &       *4.0D0*Pi*(r(2)-r(1))
        !print*,maxintpartbrd,corrintpartbrd
        corrintpartfrac=0.0D0
        corrintpartfracbrd=0.0D0
        if (maxintpart.gt.0.0D0)corrintpartfrac=corrintpart/maxintpart
        if (maxintpartbrd.gt.0.0D0)corrintpartfracbrd=corrintpartbrd
     &     /maxintpartbrd
        write(65,'(7(F15.8))') r(k),gcorr(k),gcorrbrd(k),corrintpart,
     &    corrintpartbrd,corrintpartfrac,corrintpartfracbrd 
      end do
      close(65)
      
      if(isopen12) then
            write(12,fsubendext)trim(adjustl('corr'))
      else
            print fsubendext,trim(adjustl('corr'))
      end if

      end subroutine
c---------------------------------------------------------------------

      end module
