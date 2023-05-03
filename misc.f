      module misc
      ! contains miscellaneous stuff 
      implicit none

      contains

c---------------------------------------------------------------------

      subroutine propcorr0(infile,ele,iprop,corrf, 
     &                rcut,binwidth)
      ! determines if atoms with the same property cluster in
      ! space
      use defs
      implicit none
      double precision vecs(1:3,1:3)
      double precision rcut,propav,propprodlocav,binwidth
      double precision, allocatable ::
     & props(:,:),checkprops(:),propsloc(:)
      character(len=*) infile,corrf
      !character infile*30,ele*2,line*100,elestring*100,corrf*20
      character ele*2,line*100,elestring*100
      type(atom), allocatable :: atoms(:),checkatoms(:)
      type(element),allocatable :: eles(:)
      integer natoms, nele,iele,i,j,iprop,nprops,numele,ieleprop
      integer ncheckatoms,nproploc
      ! for correlation function:
      integer :: nr
      double precision rmax,corrint,maxint
      double precision, allocatable :: gcorr(:),r(:),rdistf(:)
      parameter(rmax=20.0D0)
      logical isopen12
      !character FMT1*1024,FMT2*1024

      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
            write(12,fsubstart)trim(adjustl('propcorr0'))
      else
            print fsubstart,trim(adjustl('propcorr0'))
      end if
      
!      read(10,*) infile
!      read(10,*) ele
!      read(10,*) iprop
!      read(10,*) corrf
!      read(10,*) rcut
!      read(10,'(A100)',end=4,err=4) line
!      if (index(line,'binwidth').ge.1) then
!            read(line(index(line,'binwidth')+9:100),*,end=4,err=4) 
!     &       binwidth
!      end if
!4     continue      
      nr=int(rmax/binwidth)
      allocate(gcorr(1:nr),r(1:nr),rdistf(nr))
      write(12,'(8x,"File to read is ",A,".")') trim(adjustl(infile))
      write(12,'(8x,"Will check ",A2," for correlations")') 
     & trim(adjustl(ele))
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
!      write(12,'(8x,"File ",A," contains ",I," atoms.")')
!     & trim(infile),natoms
      WRITE(FMT2,*) natoms
      WRITE(FMT1,*) '(8x,"File ",A',len(trim(infile)),
     & '" contains ",I',len(trim(FMT2)),'" atoms.")' 
      WRITE(12,FMT1) infile,natoms
      write(12,'(8x,"Will check aux. prop. #",I3,".")'),iprop
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
      write(12,'(8x,"Cell vectors:")')
      write(12,'(8x,3(F12.6))') vecs(1,1:3)
      do i=2,3
            do j=1,3
                  read(65,'(A100)') line
                  line=adjustl(line)
                  read(line(index(line,'=')+1:100),*) vecs(i,j)
            end do
            write(12,'(8x,3(F12.6))') vecs(i,1:3)
      end do
      ! find number of auxiliary properties
      do while(index(line,'entry_count').le.0.or.line(1:1).eq.'#')
            read(65,'(A100)') line
            line=adjustl(line)
      end do
      read(line(index(line,'=')+1:100),*) nprops
      nprops=nprops-3
      write(12,'(8x,"There are ",I3," auxiliary properties.")')nprops
      ! find elements
      nele=0
      elestring(1:100)=' '
10    read(65,'(A100)',end=15) line
      line=adjustl(line)
      do i=1,100
            if (line(1:2).eq.elements(i)) then
                  elestring(2*nele+1:2*nele+2)=elements(i)
                  nele=nele+1
            end if
      end do
      goto 10
15    continue
      allocate(eles(nele))
      write(12,'(8x,"Found ",I3," elements:")') nele
      do i=1,nele
            eles(i)%name(1:2)=elestring(2*i-1:2*i)
      end do
      ! get numbers of atoms of each element
      rewind(65)
      do while (line(1:2).ne.eles(1)%name(1:2))
            read(65,'(A100)')line
            line=adjustl(line)
      end do
      do i=2,nele
            numele=0
            do while (line(1:2).ne.eles(i)%name(1:2))
                  read(65,'(A100)')line
                  line=adjustl(line)
                  numele=numele+1
            end do
            eles(i-1)%howmany=numele-2
      end do
      numele=0
20    read(65,'(A100)',end=25) line    
      line=adjustl(line)
      numele=numele+1
      goto 20
25    continue
      eles(nele)%howmany=numele
      do i=1,nele
            eles(i)%name(1:2)=elestring(2*i-1:2*i)
            write(12,'(8x,A2,I10)')eles(i)%name,eles(i)%howmany
            ! determine element number for which props are to be checked 
            if (eles(i)%name(1:2).eq.ele(1:2))ieleprop=i
      end do
      ! read atom coordinates
      allocate(atoms(natoms),props(1:natoms,1:nprops))
      rewind(65)
      do while (line(1:2).ne.eles(1)%name(1:2))
            read(65,'(A100)') line
            line=adjustl(line)
      end do
      i=1
      write(12,'(8x,"Coord. and props. of first and last atom",
     &      " of each species:")')
      do iele=1,nele-1
            do j=1,eles(iele)%howmany
                  read(65,*)atoms(i)%where(1:3),props(i,1:nprops)
                  atoms(i)%name=eles(iele)%name
                  if (j.eq.1.or.j.eq.eles(iele)%howmany)
     &                  write(12,'(8x,A2,3(F15.6),F15.6)') 
     &                  atoms(i)%name(1:2),atoms(i)%where,
     &                  props(i,1:nprops)             
                  i=i+1
            end do
            read(65,'(A100)') line
            read(65,'(A100)') line
      end do
      do j=1,eles(nele)%howmany
            read(65,*)atoms(i)%where(1:3),props(i,1:nprops)
            atoms(i)%name=eles(nele)%name
            if (j.eq.1.or.j.eq.eles(nele)%howmany)
     &            write(12,'(8x,A2,3(F15.6),F15.6)') 
     &            atoms(i)%name(1:2),atoms(i)%where,
     &            props(i,1:nprops)             
            i=i+1
      end do
      close(65)
      ! make list of atoms to check
      ncheckatoms=eles(ieleprop)%howmany
      allocate(checkatoms(ncheckatoms),checkprops(ncheckatoms),
     &          propsloc(ncheckatoms))
      j=1
      do i=1,natoms
            if(atoms(i)%name(1:2).eq.ele(1:2)) then
                  checkatoms(j)=atoms(i)
                  checkprops(j)=props(i,iprop+1)
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
!            write(66,1000) i,checkprops(i),propav,
!     &          propsloc(i),nproploc
!            j=j+1
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
      call corr0(nr,r,gcorr,ncheckatoms,checkatoms,checkprops,vecs,
     &          corrint,maxint,rdistf)
      open(66,file='CORRF.DAT',status='replace')
      write(66,'("# distance, corr. funct., radial distr. f.")')
      do i=1,nr
            write(66,'(3(F15.8))') r(i),gcorr(i),rdistf(i)
      end do
      close(66)
      write(12,'(8x,"Average of prop.",I3,":",F15.7)') iprop,propav
      write(12,'(8x,"(<p(r)> av. over sim. box)")')
      write(12,'(8x,"<p(0)p(r)> av. over sphere:",F20.5)')
     & propprodlocav
      write(12,'(8x,"Integrated corr. function corrint:",F20.5)')
     & corrint
      write(12,'(8x,"Integrated volume maxint:",F20.5)')maxint
      write(12,'(8x,"corrint/maxint:",F10.5)')corrint/maxint
      write(12,fsubend)
      end subroutine
                                                         
c---------------------------------------------------------------------

      subroutine prop_corr()
      ! determines if atoms with the same property cluster in
      ! space
      use defs
      implicit none
      double precision vecs(1:3,1:3)
      double precision rcut,propav,propprodlocav,binwidth
      double precision, allocatable ::
     & props(:,:),checkprops(:),propsloc(:)
      character infile*30,ele*2,line*100,elestring*100,corrf*20
      type(atom), allocatable :: atoms(:),checkatoms(:)
      type(element),allocatable :: eles(:)
      integer natoms, nele,iele,i,j,iprop,nprops,numele,ieleprop
      integer ncheckatoms,nproploc
      ! for correlation function:
      integer :: nr
      double precision rmax,corrint,maxint
      double precision, allocatable :: gcorr(:),r(:),rdistf(:)
      parameter(rmax=20.0D0)

      write(12,fsubstart)trim(adjustl('prop_corr'))
      
      read(10,*) infile
      read(10,*) ele
      read(10,*) iprop
      read(10,*) corrf
      read(10,*) rcut
      binwidth=0.2D0
      read(10,'(A100)',end=4,err=4) line
      if (index(line,'binwidth').ge.1) then
            read(line(index(line,'binwidth')+9:100),*,end=4,err=4) 
     &       binwidth
      end if
4     continue      
      nr=int(rmax/binwidth)
      allocate(gcorr(1:nr),r(1:nr),rdistf(nr))
      write(12,'(8x,"File to read is ",A,".")') trim(adjustl(infile))
      write(12,'(8x,"Will check ",A2," for correlations")') 
     & trim(adjustl(ele))
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
!      write(12,'(8x,"File ",A," contains ",I," atoms.")')
      WRITE(FMT2,*) natoms
      WRITE(FMT1,*) '(8x,"File ",A',len(trim(infile)),
     & '" contains ",I',len(trim(FMT2)),'" atoms.")' 
      WRITE(12,FMT1) infile,natoms
      write(12,'(8x,"Will check aux. prop. #",I3,".")'),iprop
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
      write(12,'(8x,"Cell vectors:")')
      write(12,'(8x,3(F12.6))') vecs(1,1:3)
      do i=2,3
            do j=1,3
                  read(65,'(A100)') line
                  line=adjustl(line)
                  read(line(index(line,'=')+1:100),*) vecs(i,j)
            end do
            write(12,'(8x,3(F12.6))') vecs(i,1:3)
      end do
      ! find number of auxiliary properties
      do while(index(line,'entry_count').le.0.or.line(1:1).eq.'#')
            read(65,'(A100)') line
            line=adjustl(line)
      end do
      read(line(index(line,'=')+1:100),*) nprops
      nprops=nprops-3
      write(12,'(8x,"There are ",I3," auxiliary properties.")')nprops
      ! find elements
      nele=0
      elestring(1:100)=' '
10    read(65,'(A100)',end=15) line
      line=adjustl(line)
      do i=1,100
            if (line(1:2).eq.elements(i)) then
                  elestring(2*nele+1:2*nele+2)=elements(i)
                  nele=nele+1
            end if
      end do
      goto 10
15    continue
      allocate(eles(nele))
      write(12,'(8x,"Found ",I3," elements:")') nele
      do i=1,nele
            eles(i)%name(1:2)=elestring(2*i-1:2*i)
      end do
      ! get numbers of atoms of each element
      rewind(65)
      do while (line(1:2).ne.eles(1)%name(1:2))
            read(65,'(A100)')line
            line=adjustl(line)
      end do
      do i=2,nele
            numele=0
            do while (line(1:2).ne.eles(i)%name(1:2))
                  read(65,'(A100)')line
                  line=adjustl(line)
                  numele=numele+1
            end do
            eles(i-1)%howmany=numele-2
      end do
      numele=0
20    read(65,'(A100)',end=25) line    
      line=adjustl(line)
      numele=numele+1
      goto 20
25    continue
      eles(nele)%howmany=numele
      do i=1,nele
            eles(i)%name(1:2)=elestring(2*i-1:2*i)
            write(12,'(8x,A2,I10)')eles(i)%name,eles(i)%howmany
            ! determine element number for which props are to be checked 
            if (eles(i)%name(1:2).eq.ele(1:2))ieleprop=i
      end do
      ! read atom coordinates
      allocate(atoms(natoms),props(1:natoms,1:nprops))
      rewind(65)
      do while (line(1:2).ne.eles(1)%name(1:2))
            read(65,'(A100)') line
            line=adjustl(line)
      end do
      i=1
      write(12,'(8x,"Coord. and props. of first and last atom",
     &      " of each species:")')
      do iele=1,nele-1
            do j=1,eles(iele)%howmany
                  read(65,*)atoms(i)%where(1:3),props(i,1:nprops)
                  atoms(i)%name=eles(iele)%name
                  if (j.eq.1.or.j.eq.eles(iele)%howmany)
     &                  write(12,'(8x,A2,3(F15.6),F15.6)') 
     &                  atoms(i)%name(1:2),atoms(i)%where,
     &                  props(i,1:nprops)             
                  i=i+1
            end do
            read(65,'(A100)') line
            read(65,'(A100)') line
      end do
      do j=1,eles(nele)%howmany
            read(65,*)atoms(i)%where(1:3),props(i,1:nprops)
            atoms(i)%name=eles(nele)%name
            if (j.eq.1.or.j.eq.eles(nele)%howmany)
     &            write(12,'(8x,A2,3(F15.6),F15.6)') 
     &            atoms(i)%name(1:2),atoms(i)%where,
     &            props(i,1:nprops)             
            i=i+1
      end do
      close(65)
      ! make list of atoms to check
      ncheckatoms=eles(ieleprop)%howmany
      allocate(checkatoms(ncheckatoms),checkprops(ncheckatoms),
     &          propsloc(ncheckatoms))
      j=1
      do i=1,natoms
            if(atoms(i)%name(1:2).eq.ele(1:2)) then
                  checkatoms(j)=atoms(i)
                  checkprops(j)=props(i,iprop+1)
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
!            write(66,1000) i,checkprops(i),propav,
!     &          propsloc(i),nproploc
!            j=j+1
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
      call corr0(nr,r,gcorr,ncheckatoms,checkatoms,checkprops,vecs,
     &          corrint,maxint,rdistf)
      open(66,file='CORRF.DAT',status='replace')
      write(66,'("# distance, corr. funct., radial distr. f.")')
      do i=1,nr
            write(66,'(3(F15.8))') r(i),gcorr(i),rdistf(i)
      end do
      close(66)
      write(12,'(8x,"Average of prop.",I3,":",F15.7)') iprop,propav
      write(12,'(8x,"(<p(r)> av. over sim. box)")')
      write(12,'(8x,"<p(0)p(r)> av. over sphere:",F20.5)')
     & propprodlocav
      write(12,'(8x,"Integrated corr. function corrint:",F20.5)')
     & corrint
      write(12,'(8x,"Integrated volume maxint:",F20.5)')maxint
      write(12,'(8x,"corrint/maxint:",F10.5)')corrint/maxint
      write(12,fsubend)
      end subroutine
                                                         
c---------------------------------------------------------------------

      subroutine corr0(nr,r,gcorr,ncheckatoms,checkatoms,checkprops,
     &                vecs,corrint,maxint,rdistf)
      ! calculates correlation function
      use defs
      implicit none

      double precision, intent(in) :: r(nr),vecs(1:3,1:3)
      double precision gcorr(nr),rdistf(nr)
      type (atom), intent(in) :: checkatoms(ncheckatoms)
      double precision, intent(in) :: checkprops(ncheckatoms)
      integer, intent(in) :: nr, ncheckatoms
      integer i,j,k,l,m,n,nsample(1:nr),nuc(1:3)
      double precision distvec(1:3),dist,gav,alats(1:3),corrint,
     & maxint,vol,vsphere

      do i=1,3
        alats(i)=sqrt(vecs(i,1)**2+vecs(i,2)**2+vecs(i,3)**2)
        nuc(i)=ceiling(r(nr)/alats(i))
      end do
      ! volume of the sphere with radius rmax (up to which the corr.
      ! fct. and rdf are calculated)
      vsphere=4.0D0*Pi*r(nr)**3/(3.0D0)
      !print*,'nuc=',nuc
      nsample=0
      gcorr=0.0D0
      rdistf=0.0D0
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
                  rdistf(1)=rdistf(k)+1.0D0/r(1)**2
                  nsample(1)=nsample(1)+1
                end if
                do k=2,nr
                  if (dist.gt.r(k-1).and.dist.le.r(k)) then
                    gcorr(k)=gcorr(k)+checkprops(i)*checkprops(j)
                    rdistf(k)=rdistf(k)+1.0D0/r(k)**2
                    nsample(k)=nsample(k)+1
                  end if
                end do
              end do  
            end do
          end do
         end if 
        end do 
      end do
      gav=gav/dble(ncheckatoms)
      !print*,gav
      do k=1,nr
        if (nsample(k).gt.0) gcorr(k)=gcorr(k)/dble(nsample(k))-gav**2
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
      !rdistf=rdistf*vsphere/(dble(ncheckatoms**2)*4.0D0*Pi*(r(2)-r(1)))
      !integrate correlation function
      corrint=0.0D0
      maxint=0.0D0
      do k=1,nr
        corrint=corrint+r(k)**2*gcorr(k)
        maxint=maxint+r(k)**2
      end do
      corrint=corrint*4.0D0*Pi*(r(2)-r(1))
      maxint=maxint*4.0D0*Pi*(r(2)-r(1))

      end subroutine
c---------------------------------------------------------------------

      subroutine abs2frac(abscoords,cell,fraccoords)
      ! transforms absolute to fractional coordinates

      implicit none

      double precision, intent(in) :: abscoords(1:3),cell(1:3,1:3)
      double precision fraccoords(1:3),matrix(1:3,1:3),
     &       inversematrix(1:3,1:3)
      integer i,j

      do i=1,3
        do j=1,3
          matrix(i,j)=cell(j,i)
        end do
      end do
      call invert_matrix(matrix,inversematrix)
      do i=1,3
        fraccoords(i)=0.0D0
        do j=1,3
          fraccoords(i)=fraccoords(i)+inversematrix(i,j)*abscoords(j)
        end do
      end do

      end subroutine

c---------------------------------------------------------------------

      subroutine frac2abs(fraccoords,vecs,abscoords)
      ! transforms fractional to absolute coordinates

      implicit none

      double precision, intent(in) :: fraccoords(1:3),vecs(1:3,1:3)
      double precision abscoords(1:3)
      integer i,j

      abscoords=0.0D0
      do i=1,3
        do j=1,3
          abscoords(i)=abscoords(i)+fraccoords(j)*vecs(j,i)
        end do
      end do

      end subroutine

c---------------------------------------------------------------------

      subroutine invert_matrix(inmatrix,outmatrix)
      ! this subroutine inverts a 3x3x3 matrix.
      implicit none
      double precision, intent(in) :: inmatrix(1:3,1:3)
      double precision, intent(out) :: outmatrix(1:3,1:3)
      double precision :: detin,identity(1:3,1:3)
      integer i,j,k
      !---------------------------------------------------------------
      detin=inmatrix(1,1)*inmatrix(2,2)*inmatrix(3,3)
     &     +inmatrix(1,2)*inmatrix(2,3)*inmatrix(3,1)
     &     +inmatrix(1,3)*inmatrix(2,1)*inmatrix(3,2)
     &     -inmatrix(3,1)*inmatrix(2,2)*inmatrix(1,3)
     &     -inmatrix(3,2)*inmatrix(2,3)*inmatrix(1,1)
     &     -inmatrix(3,3)*inmatrix(2,1)*inmatrix(1,2)
      outmatrix(1,1)=inmatrix(2,2)*inmatrix(3,3)
     &              -inmatrix(3,2)*inmatrix(2,3)
      outmatrix(1,2)=inmatrix(1,3)*inmatrix(3,2)
     &              -inmatrix(1,2)*inmatrix(3,3)
      outmatrix(1,3)=inmatrix(1,2)*inmatrix(2,3)
     &              -inmatrix(1,3)*inmatrix(2,2)
      outmatrix(2,1)=inmatrix(2,3)*inmatrix(3,1)
     &              -inmatrix(2,1)*inmatrix(3,3)
      outmatrix(2,2)=inmatrix(1,1)*inmatrix(3,3)
     &              -inmatrix(1,3)*inmatrix(3,1)
      outmatrix(2,3)=inmatrix(1,3)*inmatrix(2,1)
     &              -inmatrix(1,1)*inmatrix(2,3)
      outmatrix(3,1)=inmatrix(2,1)*inmatrix(3,2)
     &              -inmatrix(2,2)*inmatrix(3,1)
      !outmatrix(3,2)=inmatrix(2,1)*inmatrix(3,1) ! wrong!
      outmatrix(3,2)=inmatrix(1,2)*inmatrix(3,1)
     &              -inmatrix(1,1)*inmatrix(3,2)
      outmatrix(3,3)=inmatrix(1,1)*inmatrix(2,2)
     &              -inmatrix(1,2)*inmatrix(2,1)
      outmatrix=outmatrix/detin
      ! check if inmatrix*outmatrix=1
      !identity=0.0D0
      !do i=1,3
      !  do j=1,3
      !    do k=1,3
      !      identity(i,j)=identity(i,j)+inmatrix(i,k)*outmatrix(k,j)
      !    end do
      !  end do
      !end do
      !write(12,*) '     matrix*inverted_matrix:'
      !do i=1,3
      !  write(12,'      (3(F10.6))') identity(i,1:3)
      !end do
      !write(12,*)'       Subroutine finished successfully.'
      !write(12,*)'      *********************************************'
      end subroutine

c---------------------------------------------------------------------
      
      subroutine mirror(incoords,cell,nmirror,mirrorcoords,tolbound)
      ! finds out if atom sits on cell boundary and gives the number of 
      ! mirror atoms and their coordinates

      use defs
      implicit none
      
      double precision, intent(in) :: incoords(1:3),cell(1:3,1:3)
      double precision :: mirrorcoords(1:8,1:3),fraccoords(1:3)
      double precision :: mirrorfcoords(1:3),tolbound(1:3)
      integer nmirror,i,j,k,l,m
      
      mirrorfcoords=0.0D0
      call abs2frac(incoords,cell,fraccoords)
      !write(12,*)'      fract. coords:',fraccoords 
      nmirror=0
      do i=-1,1
        do j=-1,1
          do k=-1,1
            mirrorfcoords(1)=fraccoords(1)+dble(i)
            mirrorfcoords(2)=fraccoords(2)+dble(j)
            mirrorfcoords(3)=fraccoords(3)+dble(k)
            if (mirrorfcoords(1).ge.-tolbound(1)
     &        .and.mirrorfcoords(1).le.1.0D0+tolbound(1)
     &        .and.mirrorfcoords(2).ge.-tolbound(2)
     &        .and.mirrorfcoords(2).le.1.0D0+tolbound(2)
     &        .and.mirrorfcoords(3).ge.-tolbound(3)
     &        .and.mirrorfcoords(3).le.1.0D0+tolbound(3)) then
              nmirror=nmirror+1
              mirrorcoords(nmirror,1:3)=0.0D0
              do l=1,3
                do m=1,3
                  mirrorcoords(nmirror,l)=mirrorcoords(nmirror,l)
     &             +mirrorfcoords(m)*cell(m,l)
                end do
              end do
              end if
            end do
          end do
        end do
      !write(12,9000)nmirror,fraccoords
 9000 format('      nmirror=',I2,', xfrac=',3(F10.5))    
      !do i=1,nmirror
        !write(12,9001) mirrorcoords(i,1:3)
      !end do
 9001 format(3(F15.10))
      end subroutine

c---------------------------------------------------------------------

      subroutine writecoorat(atoms,natoms,vecs,outfile)
      ! writes coordinates to MBPP COORAT file
      
      use defs
      implicit none
      character line*100,nameat*8
      character (LEN=*), intent(in) :: outfile
      double precision vecs(1:3,1:3),xabs(1:3),xfrac(1:3)
      integer i,iinf,idum,natoms,nele,j,numele,iele
      logical newele
      type(atom) :: atoms(natoms)
      type(element), allocatable :: eles(:)

      open(51,file=outfile,status='replace')
      do i=1,natoms
            xabs(1:3)=atoms(i)%where(1:3)
            call abs2frac(xabs,vecs,xfrac)
            atoms(i)%where(1:3)=xfrac(1:3)
      end do
      close(50)
      ! count elements and number of atoms of each element
      atoms(1:natoms)%written=.false.
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
      !write(12,'(8x,"Found ",I4," elements.")')nele
      !write(12,'(8x,"Elements, number of atoms of each element:")')
      !do i=1,nele
      !      write(12,'(8x,A8,I6)') eles(i)%name,eles(i)%howmany
      !end do
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

      !write(12,fsubend)
!      return

      ! errors:
!100   write(12,'(8x,"Error: file INPUT.COFIMA ended unexpectedly")')
       !nerr=nerr+1
!      return
!110   write(12,'(8x,"Error: CONFIG file ended unexpectedly")')
!      close(50)
!      close(51)
       !nerr=nerr+1
!      return

      end subroutine

c---------------------------------------------------------------------

      subroutine sxyz2cfg()
      ! converts Adham's special xyz format to cfg
      use defs
      implicit none
      integer natoms,nele,i,j,iele,numele,idum
      character infile*25,outfile*25
      character line*100,chardum*10,propnames(1:4)*20
      logical newele
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: eles(:)
      double precision alat,vecs(1:3,1:3)
      double precision, allocatable :: props(:,:)
      !character FMT1*1024,FMT2*1024

      write(12,fsubstart)trim(adjustl('sxyz2cfg'))
      read(10,*) infile
      read(10,*) outfile
      open(55,file=infile,status='old')
      ! read number of atoms
      read(55,*) natoms
      ! read atom coordinates and properties
      allocate(atoms(natoms),props(1:natoms,1:4))
      read(55,'(A100)') line
      do i=1,natoms
            read(55,*) atoms(i)%name,atoms(i)%where,props(i,1:4)
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
      ! read lattice constants, masses,property names from input file
      read(55,'(A100)') line
      read(55,*) alat
      read(55,'(A100)') line
      do i=1,3
            read(55,*) vecs(i,1:3)
      end do
      do iele=1,nele
            read(55,*) chardum,chardum,eles(iele)%mass
      end do
      read(55,'(A100)') line
      do i=1,4
            read(55,*) chardum,idum,propnames(i)
      end do
      close(55)
      ! write cfg file
      open(56,file=outfile,status='replace')
!      write(56,'("Number of particles =",I)') natoms
      WRITE(FMT2,*) natoms
      WRITE(FMT1,*) '("Number of particles =",I',len(trim(FMT2)),')' 
      WRITE(56,FMT1) natoms
      write(56,'("A = 1 Angstrom (basic length-scale)")')
      do i=1,3
          write(56,'("H0(",I1,",1) = ",F12.6," A",/,
     &     "H0(",I1,",2) = ",F12.6," A",/,
     &     "H0(",I1,",3) = ",F12.6," A")')
     &     i,vecs(i,1),i,vecs(i,2),i,vecs(i,3)
      end do
      write(56,'(".NO_VELOCITY.")')
      write(56,'("entry_count = 7")')
      do iele=1,nele
            write(56,'(F10.5)') eles(iele)%mass
            write(56,'(8A)') eles(iele)%name
            do j=1,natoms
                  if (atoms(j)%name.eq.eles(iele)%name)
     &           write(56,1000) atoms(j)%where,props(j,1:4) 
            end do
      end do
      close(56)
 1000 format(3(F15.7),4(F12.5))
      write(12,fsubend)
      end subroutine

c---------------------------------------------------------------------

      subroutine sxyz2config()
      ! converts Adham's special xyz (sxyz) format to DLPOLY CONFIG
      ! format
      use defs
      implicit none

      integer natoms,nele,i,j,iele,numele,idum,imcon
      character infile*25,outfile*25
      character line*100,chardum*10,propnames(1:4)*20
      logical newele
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: eles(:)
      double precision alat,vecs(1:3,1:3),x,y,z
      double precision, allocatable :: props(:,:)

      write(12,fsubstart)trim(adjustl('sxyz2config'))
      read(10,*) infile
      read(10,*) outfile
      read(10,*)imcon
      open(55,file=infile,status='old')
      ! read number of atoms
      read(55,*) natoms
      ! read atom coordinates and properties
      allocate(atoms(natoms),props(1:natoms,1:4))
      read(55,'(A100)') line
      do i=1,natoms
            read(55,*) atoms(i)%name,atoms(i)%where,props(i,1:4)
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
      ! read lattice constants, masses,property names from input file
      read(55,'(A100)') line
      read(55,*) alat
      read(55,'(A100)') line
      do i=1,3
            read(55,*) vecs(i,1:3)
      end do
      do iele=1,nele
            read(55,*) chardum,chardum,eles(iele)%mass
      end do
      read(55,'(A100)') line
      do i=1,4
            read(55,*) chardum,idum,propnames(i)
      end do
      close(55)
      ! write cfg file
      open(56,file=outfile,status='replace')
!      write(56,'("Number of particles =",I)') natoms
!      write(56,'("A = 1 Angstrom (basic length-scale)")')
!      do i=1,3
!          write(56,'("H0(",I1,",1) = ",F12.6," A",/,
!     &     "H0(",I1,",2) = ",F12.6," A",/,
!     &     "H0(",I1,",3) = ",F12.6," A")')
!     &     i,vecs(i,1),i,vecs(i,2),i,vecs(i,3)
!      end do
!      write(56,'(".NO_VELOCITY.")')
!      write(56,'("entry_count = 7")')
!      do iele=1,nele
!            write(56,'(F10.5)') eles(iele)%mass
!            write(56,'(8A)') eles(iele)%name
!            do j=1,natoms
!                  if (atoms(j)%name.eq.eles(iele)%name)
!     &           write(56,1000) atoms(j)%where,props(j,1:4) 
!            end do
!      end do


      write(56,*)' '
      write(56,'(I10,I10,I16)') 0,imcon,natoms
      do i=1,3
        write(56,'(3(F20.10))') vecs(i,1:3)
      end do
      do i=1,natoms
        x=atoms(i)%where(1)
        y=atoms(i)%where(2)
        z=atoms(i)%where(3)
        atoms(i)%where(1:3)=x*vecs(1,1:3)+y*vecs(2,1:3)
     &          +z*vecs(3,1:3)
        write(56,'(A8,I12)') atoms(i)%name,i
        write(56,'(F16.9,2(F20.9))') atoms(i)%where(1:3)
      end do  
      close(56)

      write(12,fsubend)
      end subroutine

c---------------------------------------------------------------------

      subroutine cfg_del_species()
      ! removes one or more species from cfg file

      use defs
      implicit none
      character line*100,nextline*100,ladjline*100,
     &    nextladjline*100,infile*30,outfile*30
      integer i,n2rm,natoms,neles,iele,nread,nleft,j
      character, allocatable :: specnames(:)*2
      logical delete
      type(element),allocatable :: eles(:)
      !character FMT1*1024,FMT2*1024

      write(12,fsubstart)trim(adjustl('cfg_del_species'))
      read(10,*,err=500) infile
      read(10,*,err=500) outfile
      read(10,*,err=500) n2rm
      allocate(specnames(n2rm))
      do i=1,n2rm
            read(10,*,err=500) specnames(i)
      end do
      write(12,'(8x,"Will read file ",A,".")') trim(infile)
      write(12,'(8x,"Will write file ",A,".")') trim(outfile)
      write(12,'(8x,"Will delete the following species:")')
      do i=1,n2rm
            write(12,'(8x,A4)') specnames(i)
      end do

      ! read cfg file, how many elements are there?
      neles=0
      open(55,file=infile,status='old',err=600)
2     read(55,'(A100)',err=600,end=4) line
      ladjline=adjustl(line)
      do i=1,100
        if(ladjline(1:2).eq.elements(i)) neles=neles+1
      end do
      if(ladjline(1:1).ne.'#') then
          if (index(line,'Number of particles').ge.1) then
              read(line(index(line,'=')+1:100),*,err=600)natoms
!              write(12,'(8x,I," atoms in original cfg file.")')
!     &              natoms
              WRITE(FMT2,*) natoms
              WRITE(FMT1,*) '(I',len(trim(FMT2)),'," atoms in original c
     &fg file.")' 
              WRITE(12,FMT1) natoms
          end if
      end if
      goto 2

      ! get names and abundancies of elements
4     continue
      allocate(eles(neles))
      eles%howmany=0
!      write(12,'(8x,"Number of elements in cfg file: ",I)') neles
      WRITE(FMT2,*) neles
      WRITE(FMT1,*) '(8x, "Number of elements in cfg file: ",I',
     &    len(trim(FMT2)),')' 
      WRITE(12,FMT1) neles
      rewind(55)
      iele=0
6     read(55,'(A100)',end=8) line
      ladjline=adjustl(line)
      if(iele.gt.0) then
        if(ladjline(1:1).ne.'#') 
     &     eles(iele)%howmany=eles(iele)%howmany+1     
      end if
      do i=1,100
        if(ladjline(1:2).eq.elements(i)) then
              iele=iele+1
              eles(iele)%name='        '
              eles(iele)%name(1:2)=elements(i)
        end if
      end do
      goto 6

8     continue
      nleft=natoms
      do i=1,neles
        if(i.lt.neles)eles(i)%howmany=eles(i)%howmany-2     
!        write(12,'(8x,A2,I)') eles(i)%name(1:2),eles(i)%howmany
        WRITE(FMT2,*) eles(i)%howmany
        WRITE(FMT1,*) '(8x,A2,I',len(trim(FMT2)),')' 
        WRITE(12,FMT1) eles(i)%name(1:2),eles(i)%howmany
        do j=1,n2rm
          if(eles(i)%name(1:2).eq.specnames(j))
     &       nleft=nleft-eles(i)%howmany
        end do
      end do
!      write(12,'(8x,"# of atoms left: ",I)') nleft
      WRITE(FMT2,*) nleft
      WRITE(FMT1,*) '(8x,"# of atoms left: ",I',len(trim(FMT2)),')' 
      WRITE(12,FMT1) nleft

      ! read and write header of cfg file
      open(56,file=outfile,status='replace')
      rewind(55)
10    read(55,'(A100)',end=600,err=600)line
      ladjline=adjustl(line)
      if(index(line,'Number of particles').ge.1.and.
     &  ladjline(1:1).ne.'#')then 
!          write(56,'("Number of particles = ",I)') nleft 
          WRITE(FMT2,*) nleft
          WRITE(FMT1,*) '("Number of particles =",I',len(trim(FMT2)),')'
          WRITE(56,FMT1) nleft
      else
         write(56,'(A100)') line
      end if
      if(ladjline(1:1).ne.'#') then
          if(index(line,'entry_count').ge.1)goto 20
      end if
      goto 10

      ! read and write rest of cfg file 
20    continue
      ! loop through elements and determine if it is to be deleted
      do iele=1,neles
        delete=.false.
        do i=1,n2rm
          if(eles(iele)%name(1:2).eq.specnames(i)) delete=.true.
        end do  
        ! if this species is to be deleted, then go to the next element
        if(delete) then
              write(12,'(8x,"Found species to be deleted: ",A2)')
     &             eles(iele)%name(1:2)
              nread=0
              do while (nread.lt.eles(iele)%howmany+2)
                read(55,'(A100)',end=30) line
                ladjline=adjustl(line)
                if (ladjline(1:1).ne.'#') nread=nread+1
              end do  
        ! else read and write coordinates and properties     
        else
              write(12,'(8x,"This species is to be kept: ",A2)')
     &             eles(iele)%name(1:2)
              nread=0
              do while (nread.lt.eles(iele)%howmany+2)
                read(55,'(A100)',end=30) line
                ladjline=adjustl(line)
                if (ladjline(1:1).ne.'#') nread=nread+1
                write(56,'(A100)') line
              end do  
        end if
      end do

30    continue  
      close(55)
      close(56)
      write(12,fsubend)
      RETURN

      ! Error messages
500   print*,"Error when reading INPUT.COFIMA."
      write(12,'(8x,"Error when reading INPUT.COFIMA.")')
      nerr=nerr+1
      return
600   print*,"Error: something is wrong with the input cfg file."
      write(12,'(8x,"Error when reading the cfg file.")')
      nerr=nerr+1
      close(55)
      close(56)
      return

      end subroutine

c---------------------------------------------------------------------

      subroutine rm_double_atoms(atoms,tolat)
      !  
      use defs
      implicit none
      !
      type(atom), allocatable :: atoms(:)
      integer nallatoms
      double precision tolat(3)
      !
      type(atom), allocatable :: newatoms(:)
      integer nnewatoms
      integer n,p,iatom
      logical,allocatable :: delete(:)
      !
      nallatoms=size(atoms)
      allocate(delete(1:nallatoms)) 
      delete=.False.
      do n=1,nallatoms-1 
            do p=n+1,nallatoms 
              if(abs(atoms(p)%where(1)-atoms(n)%where(1))
     &          .le.tolat(1).and.abs(atoms(p)%where(2)
     &          -atoms(n)%where(2)).le.tolat(2).and.
     &               abs(atoms(p)%where(3)-atoms(n)%where(3))
     &               .le.tolat(3).and.
     &                  atoms(p)%core.eq.atoms(n)%core) then
                        delete(p)=.True.
              end if
            end do
      end do

      ! find final number of coordinates
      nnewatoms=0
      do n=1,nallatoms 
            if(delete(n).eqv..False.) nnewatoms=nnewatoms+1
      end do
      ! get new atoms
      allocate(newatoms(nnewatoms))
      iatom=0
      do n=1,nallatoms !natoms*nuc(1)*nuc(2)*nuc(3)
            if (delete(n).eqv..False.) then
                  iatom=iatom+1
                  newatoms(iatom)=atoms(n)
            end if
      end do
      !
      deallocate(atoms)
      allocate(atoms(nnewatoms))
      atoms=newatoms
      deallocate(newatoms)
      !
      !
      !
      end subroutine rm_double_atoms

c---------------------------------------------------------------------

      double precision function propavloctheta(place,ncheckatoms,
     &                  checkatoms,checkprops,rcut,vecs,nproploc)
      use defs
      implicit none
      double precision, intent(in) :: place(1:3),
     & checkprops(1:ncheckatoms),rcut,vecs(1:3,1:3)
      integer, intent(in) :: ncheckatoms
      type(atom), intent(in) :: checkatoms(1:ncheckatoms)
      integer i,j,k,l,m,nproploc
      double precision distvec(1:3),distance

      propavloctheta=0.D0
      nproploc=0
      do i=1,ncheckatoms
        do k=-1,1
          do l=-1,1
            do m=-1,1
              ! distance in fractional coordinates, take also atoms
              ! across the cell boundaries into account
              distvec(1)=(checkatoms(i)%where(1)+dble(k)-place(1))
              distvec(2)=(checkatoms(i)%where(2)+dble(l)-place(2))
              distvec(3)=(checkatoms(i)%where(3)+dble(m)-place(3))
              ! distance in absolute coordinates
              distvec(1:3)=distvec(1)*vecs(1,1:3)+distvec(2)*vecs(2,1:3)
     &                  +distvec(3)*vecs(3,1:3)       
              distance=sqrt(distvec(1)**2+distvec(2)**2+distvec(3)**2)
              if(distance.le.rcut.and.distance.gt.1E-3) then
                  propavloctheta=propavloctheta+checkprops(i)
                  nproploc=nproploc+1
              end if
            end do
          end do
        end do
      end do
      if(nproploc.gt.0) propavloctheta=propavloctheta/dble(nproploc)

      end function

c---------------------------------------------------------------------      
      
      subroutine refcoords()
      ! read a coordinate file and a reference coordinate file 
      ! and write a new coordinate file in which the original
      ! coordinates are replaced by the closest lying reference
      ! coordinates

      use defs
      implicit none
      integer nprops,nele,iele,iprop,i,j,natoms,ieleprop
      integer nrefatoms,nrefele,nrefprops,k
      character filetype*10,infile*30,reffile*30,outfile*30,line*100
      double precision vecs(1:3,1:3),refvecs(1:3,1:3),dumvecs(1:3,1:3)
      double precision, allocatable :: props(:,:),refprops(:,:)
      type(element),allocatable :: eles(:),refeles(:)
      type(atom), allocatable :: atoms(:),refatoms(:)
      !character FMT1*1024,FMT2*1024  

      write(12,fsubstart)trim(adjustl('refcoords'))
      ! read file format from input file
      read(10,*,end=100,err=100) line
      line=adjustl(line)
      do while (line(1:1).eq.'#')
            read(10,'(A100)',err=100,end=100) line
            line=adjustl(line)
      end do
      read(line,*,err=100) filetype
      write(12,'(8x,"Will read and write ",A," files.")')
     & trim(adjustl(filetype))
      ! read filenames from INPUT.COFIMA
      ! input file
      read(10,*,end=100,err=100) line
      line=adjustl(line)
      do while (line(1:1).eq.'#')
            read(10,'(A100)',err=100,end=100) line
            line=adjustl(line)
      end do
      read(line,*,err=100) infile
      write(12,'(8x,"Will read coordinates from file ",A,".")')
     &      trim(adjustl(infile)) 
      ! reference file
      read(10,*,end=100,err=100) line
      line=adjustl(line)
      do while (line(1:1).eq.'#')
            read(10,'(A100)',err=100,end=100) line
            line=adjustl(line)
      end do
      read(line,*,err=100) reffile
      write(12,'(8x,"Will read ref. coordinates from file ",A,".")')
     &      trim(adjustl(reffile)) 
      ! output file
      read(10,*,end=100,err=100) line
      line=adjustl(line)
      do while (line(1:1).eq.'#')
            read(10,'(A100)',err=100,end=100) line
            line=adjustl(line)
      end do
      read(line,*,err=100) outfile
      write(12,'(8x,"Will write coordinates to file ",A,".")')
     &      trim(adjustl(outfile))
      
      !
      select case(filetype)
            case('cfg')
            ! read input cfg file 
            open(65,file=infile,status='old')
            call readcfg(65,nele,eles,natoms,atoms,vecs,nprops,
     & props)
            close(65)
!            write(12,'(8x,"There are ",I," atoms in ",A,".")')
!     &            natoms,trim(infile)       
            WRITE(FMT2,*) natoms
            WRITE(FMT1,*) '(8x,"There are ",I',len(trim(FMT2)),
     &          '" atoms in ",A',len(trim(infile)),',".")'
            WRITE(12,FMT1) natoms,trim(infile)
            write(12,'(8x,"There are ",I3," auxiliary properties.")')
     &       nprops
            write(12,'(8x,"Found ",I3," elements:")') nele
            do i=1,nele
                  write(12,'(8x,A2,I10,F10.4)')eles(i)%name,
     &                  eles(i)%howmany,eles(i)%mass
            end do
            i=1
            write(12,'(8x,"Coord. and props. of first and last atom",
     &      " of each species:")')
            do iele=1,nele-1
                  do j=1,eles(iele)%howmany
                        if (j.eq.1.or.j.eq.eles(iele)%howmany)
     &                  write(12,'(8x,A2,3(F15.6),F15.6)') 
     &                  atoms(i)%name(1:2),atoms(i)%where,
     &                  props(i,1:nprops)             
                        i=i+1
                  end do
            end do
            do j=1,eles(nele)%howmany
                  if (j.eq.1.or.j.eq.eles(nele)%howmany)
     &            write(12,'(8x,A2,3(F15.6),F15.6)') 
     &            atoms(i)%name(1:2),atoms(i)%where,
     &            props(i,1:nprops)             
                  i=i+1
            end do
            ! read reference coordinates
            open(66,file=reffile,status='old', err=660)
            call readcfg(66,nrefele,refeles,nrefatoms,refatoms,
     &           refvecs,nrefprops,refprops)
            close(66)
!            write(12,'(8x,"There are",I," atoms in the file ",A,".")')
!     &            nrefatoms,trim(adjustl(reffile))
            WRITE(FMT2,*) nrefatoms
            WRITE(FMT1,*) '(8x,"There are ",I',len(trim(FMT2)),
     &          '" atoms in the file ",A',len(trim(infile)),',".")'
            WRITE(12,FMT1) nrefatoms,trim(adjustl(reffile))
            write(12,'(8x,"Ref. cell vectors:")')
            write(12,'(8x,3(F12.6))') refvecs(1,1:3)
            do i=2,3
                  write(12,'(8x,3(F12.6))') refvecs(i,1:3)
            end do
            write(12,'(8x,"Found ",I3," elements:")') nrefele
            do i=1,nrefele
                  write(12,'(8x,A2,I10)')refeles(i)%name,
     &                refeles(i)%howmany
            end do
            i=1
            write(12,'(8x,"Coord. and props. of first and last atom",
     &      " of each species:")')
            do iele=1,nrefele-1
                  do j=1,refeles(iele)%howmany
                        if (j.eq.1.or.j.eq.refeles(iele)%howmany)
     &                  write(12,'(8x,A2,3(F15.6),F15.6)') 
     &                  refatoms(i)%name(1:2),refatoms(i)%where
                        i=i+1
                  end do
            end do
            do j=1,refeles(nrefele)%howmany
                  if (j.eq.1.or.j.eq.refeles(nrefele)%howmany)
     &            write(12,'(8x,A2,3(F15.6),F15.6)') 
     &            refatoms(i)%name(1:2),refatoms(i)%where
                  i=i+1
            end do
            ! replace atom coordinates by reference coordinates
            ! for fractional coordinates, take cartesian unity vectors
            dumvecs=0.0D0
            do i=1,3
              dumvecs(i,i)=1.0D0
            end do  
            do i=1,natoms
              call nearestatom(atoms(i),refatoms,nrefatoms,dumvecs )
            end do
            ! write output cfg 
            open(67,file=outfile,status='replace')
            call writecfg(67,nele,eles,natoms,atoms,vecs,nprops,
     & props)
            close(67)
      case default
            ! error
            goto 110
      end select

      write(12,fsubend)
      return
      ! error messages
100   print*,"Error: Something is wrong with INPUT.COFIMA."
      write(12,'(8x,"Error: Something is wrong with INPUT.COFIMA.")')
      nerr=nerr+1
      return
110   print*,"Error: Unknown filetype."
      write(12,'(8x,"Error: Unknown filetype.")')
      nerr=nerr+1
      return
660   print*,"Error: something is wrong with the reference file."
      write(12,'(8x,"Error: something is wrong with the reference
     & file.")')
      nerr=nerr+1
      close(66)
      close(67)
      return
      end subroutine

c---------------------------------------------------------------------

      subroutine readcfg(funit,nele,eles,natoms,atoms,vecs,nprops,
     & props)
      ! reads a cfg file
      
      use defs
      implicit none
      integer natoms,nele,numele,iele,i,j,nprops
      integer funit,ieleprop
      character line*100,elestring*1000
      double precision vecs(1:3,1:3)
      double precision, allocatable :: props(:,:)
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: eles(:)

      ! read input cfg file and write header of output file
      ! read number of atoms from cfg file
      read(funit,'(A100)',end=700,err=700) line
      line=adjustl(line)
      do while (index(line,'Number of particles').le.0
     &          .or.line(1:1).eq.'#')
            read(funit,'(A100)',err=700,end=700) line
                  line=adjustl(line)
            end do
            read(line(index(line,'=')+1:100),*) natoms
            ! read cell vectors from  cfg file
            read(funit,'(A100)',end=700,err=700) line
            line=adjustl(line)
            do while (index(line,'H0(1,1)').le.0.or.line(1:1).eq.'#')
            read(funit,'(A100)',end=700,err=700) line
                  line=adjustl(line)
            end do
            read(line(index(line,'=')+1:100),*) vecs(1,1)
            do j=2,3
            read(funit,'(A100)',end=700,err=700) line
                  line=adjustl(line)
                  read(line(index(line,'=')+1:100),*) vecs(1,j)
            end do
            do i=2,3
                  do j=1,3
                  read(funit,'(A100)',err=700,end=700) line
                        line=adjustl(line)
                        read(line(index(line,'=')+1:100),*) vecs(i,j)
                  end do
            end do
            ! find number of auxiliary properties
            do while(index(line,'entry_count').le.0.or.line(1:1).eq.'#')
            read(funit,'(A100)',err=700,end=700) line
                  line=adjustl(line)
            end do
            read(line(index(line,'=')+1:100),*) nprops
            nprops=nprops-3
            ! find elements
            nele=0
            elestring(1:1000)=' '
10          read(funit,'(A100)',end=15,err=700) line
            line=adjustl(line)
            do i=1,100
                  if (line(1:2).eq.elements(i)) then
                        elestring(2*nele+1:2*nele+2)=elements(i)
                        nele=nele+1
                  end if
            end do
            goto 10
15          continue
            allocate(eles(nele))
            do i=1,nele
                  eles(i)%name(1:2)=elestring(2*i-1:2*i)
            end do
            ! get numbers of atoms of each element
            rewind(funit)
            do while (line(1:2).ne.eles(1)%name(1:2))
                  read(funit,'(A100)')line
                  line=adjustl(line)
            end do
            do i=2,nele
                  numele=0
                  do while (line(1:2).ne.eles(i)%name(1:2))
                  read(funit,'(A100)',err=700,end=700)line
                        line=adjustl(line)
                        numele=numele+1
                  end do
                  eles(i-1)%howmany=numele-2
            end do
            numele=0
20          read(funit,'(A100)',end=25,err=700) line    
            line=adjustl(line)
            numele=numele+1
            goto 20
25          continue
            eles(nele)%howmany=numele
            do i=1,nele
                  eles(i)%name(1:2)=elestring(2*i-1:2*i)
            end do
            ! read atom coordinates and masses
            allocate(atoms(natoms),props(1:natoms,1:nprops))
            rewind(funit)
            do while (line(1:2).ne.eles(1)%name(1:2))
                  read(funit,'(A100)',end=700,err=700) line
                  line=adjustl(line)
                  if(index(line,'entry_count').ge.1
     &               .and.line(1:1).ne.'#') then
                        read(funit,'(A100)',err=700,end=700) line
                        line=adjustl(line)
                        do while (line(1:1).eq.'#')
                              read(funit,'(A100)',err=700,end=700) 
     &                             line
                              line=adjustl(line)
                        end do
                        read(line,*)eles(1)%mass
                  end if
            end do
            i=1
            do iele=1,nele-1
                  do j=1,eles(iele)%howmany
                    read(funit,*)atoms(i)%where(1:3),props(i,1:nprops)
                        atoms(i)%name=eles(iele)%name
                        i=i+1
                  end do
                  read(funit,'(A100)',end=700,err=700) line
                  line=adjustl(line)
                  do while (line(1:1).eq.'#')
                        read(funit,'(A100)') line
                        line=adjustl(line)
                  end do
                  read(line,*) eles(iele+1)%mass
                  read(funit,'(A100)',err=700,end=700) line
                  do while (line(1:1).eq.'#')
                  read(funit,'(A100)',err=700,end=700) line
                        line=adjustl(line)
                  end do
            end do
            do j=1,eles(nele)%howmany
                  read(funit,*)atoms(i)%where(1:3),props(i,1:nprops)
                  atoms(i)%name=eles(nele)%name
                  i=i+1
            end do

      return
      
      ! Error messages
700   print*,"Error: something wrong with cfg file."
      close(funit)
      return
      end subroutine
      
c---------------------------------------------------------------------

      subroutine writecfg(funit,nele,eles,natoms,atoms,vecs,nprops,
     & props)
      ! writes a cfg file
      
      use defs
      implicit none
      integer natoms,nele,numele,iele,i,j,nprops
      integer funit,ieleprop
      character line*100,elestring*1000
      double precision vecs(1:3,1:3)
      double precision, allocatable :: props(:,:)
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: eles(:)

      ! write header of cfg file
!      write(funit,'("Number of particles = ",I)',err=700)
!     &       natoms 
      WRITE(FMT2,*) natoms
      WRITE(FMT1,*) '("Number of particles = ",I',len(trim(FMT2)),')'
      WRITE(funit,FMT1,err=700) natoms
      write(funit,'("A = 1 Angstrom (basic length-scale)")',err=700) 
      do i=1,3
        do j=1,3 
          write(funit,'("H0(",I1,",",I1,") = ",F12.6," A ")',err=700) 
     &       i,j,vecs(i,j)
        end do
      end do
      write(funit,'(".NO_VELOCITY.")',err=700)
      write(funit,'("entry_count = ",I4)',err=700)nprops+3
      ! write atom masses, names, and coordinates
      do iele=1,nele
            write(funit,'(F12.5)',err=700) eles(iele)%mass
            write(funit,'(A2)',err=700) eles(iele)%name(1:2)
            do j=1,natoms
              if(atoms(j)%name(1:2).eq.eles(iele)%name(1:2))then
!                write(funit,'(3(F15.8),F)',err=700)atoms(j)%where,
!     &                props(j,1:nprops)
                write(formatcfgcoords(11:13),'(I3)')nprops    
                write(funit,formatcfgcoords,err=700)atoms(j)%where,
     &                props(j,1:nprops)
              end if
            end do
      end do
      return
      
      ! Error messages
700   print*,"Error: Problem writing cfg file."
      nerr=nerr+1
      close(funit)
      return
      end subroutine

c---------------------------------------------------------------------

      subroutine nearestatom(thisatom,refatoms,nrefatoms,vecs)
      ! replaces atom coordinate by nearest reference coordinate

      use defs
      implicit none
      integer, intent(in) :: nrefatoms
      type(atom), intent(in) :: refatoms(1:nrefatoms)
      type(atom) :: thisatom
      integer i,i1,i2,i3
      double precision distvec(1:3),dist,distmin,refcoord(1:3)
      double precision, intent(in) :: vecs(1:3,1:3)

      distmin=1000.0D0
      refcoord=refatoms(1)%where
      do i=1,nrefatoms
        do i1=-1,1
          do i2=-1,1
            do i3=-1,1
              distvec(1:3)=thisatom%where(1:3)-refatoms(i)%where(1:3)
     &            +i1*vecs(1,1:3)+i2*vecs(2,1:3)+i3*vecs(3,1:3) 
              dist=sqrt(distvec(1)**2+distvec(2)**2+distvec(3)**2)
              if(dist.lt.distmin) then
                  distmin=dist
                  refcoord=refatoms(i)%where
              end if
            end do
          end do
        end do  
      end do
      thisatom%where=refcoord
      
      return
      end subroutine

c---------------------------------------------------------------------

      subroutine matchatoms(atoms1,atoms2,vecs1,vecs2)
      ! reorders atoms2 such that atom2 in atoms2 is closest to atom1 in
      ! atoms1 etc.

      use defs
      implicit none
      type(atom), intent(in) :: atoms1(:)
      type(atom), intent(inout) :: atoms2(:)
      double precision, intent(in) :: vecs1(1:3,1:3),vecs2(1:3,1:3)
      ! local variables
      type(atom) :: atoms0(1)
      integer i,j,j0
      integer natoms1,natoms2
      double precision dist,distmin,coords1(1:3),coords2(1:3)

      if (talk) then
        print fsubstart,'matchatoms'
      end if
      natoms1=size(atoms1)
      natoms2=size(atoms2)
      if (natoms2.lt.natoms1) call                                        &
     &   error_stop('atoms2 must have at least as many atoms as atoms1')
      if (natoms1.ne.natoms2) call warning('different numbers of atoms')
      do i=1,natoms1
        !
        ! begin find the atom j0 among atoms2 that matches the ith atom in
        ! atoms1
        !
        distmin=1000.0D0
        j0=-1
        do j=i,natoms2
          if (atoms2(j)%name(1:2).eq.atoms1(i)%name(1:2).and.             &
     &        atoms1(i)%core(1:4).eq.atoms2(j)%core(1:4)) then         
            coords1(1:3)=atoms1(i)%abswhere(1:3) 
            coords2(1:3)=atoms2(j)%abswhere(1:3) 
            dist=norm2(coords2-coords1) 
            if(dist.lt.distmin) then
              distmin=dist
              j0=j
            end if ! dist.lt.distmin 
          end if ! (same atom name and both are core or both are shell)
        end do ! j-loop over atoms2
        !
        if (j0.lt.1) call error_stop('atching atom could not be found')
        !
        ! end find the atom j0 among atoms2 that matches the ith atom in
        ! atoms1
        !
        !
        ! begin exchange atoms2(i) with atoms2(j0)
        ! 
        atoms0(1)=atoms2(i)
        atoms2(i)=atoms2(j0)
        atoms2(j0)=atoms0(1)
        !
        ! end exchange atoms2(i) with atoms2(j0)
        ! 
      end do ! i-loop over atoms1
      
      if (talk) then
        print fsubendext,'matchatoms'
      end if
      return

      end subroutine matchatoms

c---------------------------------------------------------------------

      subroutine mymatmul(inmat1,inmat2,outmat)
      
      implicit none
      
      double precision, intent(in) :: inmat1(1:3,1:3),inmat2(1:3,1:3)
      double precision  :: outmat(1:3,1:3)
      integer :: i,j,k  

      outmat=0.0D0

      do i=1,3
        do j=1,3
          do k=1,3
            outmat(i,j)=outmat(i,j)+inmat1(i,k)*inmat2(k,j)
          end do
        end do
      end do

      end subroutine

c---------------------------------------------------------------------
      function mat_times_vec(matrix,vec)
      ! multiplies a vector with a matrix.
     
      use defs
      implicit none

      ! variables
      double precision, DIMENSION(3) :: mat_times_vec
      double precision, intent(in) :: matrix(1:3,1:3)
      double precision, intent(inout) :: vec(1:3)
 
      mat_times_vec=0.0D0
      mat_times_vec(1)=vec(1)*matrix(1,1)
     &          +vec(2)*matrix(1,2)
     &          +vec(3)*matrix(1,3)
      mat_times_vec(2)=vec(1)*matrix(2,1)
     &          +vec(2)*matrix(2,2)
     &          +vec(3)*matrix(2,3)
      mat_times_vec(3)=vec(1)*matrix(3,1)
     &          +vec(2)*matrix(3,2)
     &          +vec(3)*matrix(3,3)
      
      end function mat_times_vec

c---------------------------------------------------
!      
!      subroutine mat_times_vec(matrix,vec,vec_out)
!      ! multiplies a vector with a matrix.
!     
!      use defs
!      implicit none
!
!      ! variables
!      double precision, DIMENSION(3) :: mat_times_vec
!      double precision, intent(in) :: matrix(1:3,1:3)
!      double precision, intent(in) :: vec(1:3)
!      double precision, intent(out) :: vec_out(1:3)
! 
!      vec_out=0.0D0
!      vec_out(1)=vec(1)*matrix(1,1)
!     &          +vec(2)*matrix(1,2)
!     &          +vec(3)*matrix(1,3)
!      vec_out(2)=vec(1)*matrix(2,1)
!     &          +vec(2)*matrix(2,2)
!     &          +vec(3)*matrix(2,3)
!      vec_out(3)=vec(1)*matrix(3,1)
!     &          +vec(2)*matrix(3,2)
!     &          +vec(3)*matrix(3,3)
!      
!      end subroutine mat_times_vec
!
c-------------------------------------------------------      

      subroutine transpon_matrix(inmat,outmat)
      
      implicit none
      
      double precision, intent(in) :: inmat(1:3,1:3)
      double precision  :: outmat(1:3,1:3)
      integer :: i,j  

      do i=1,3
        do j=1,3
          outmat(i,j)=inmat(j,i)
        end do
      end do

      end subroutine

c---------------------------------------------------------------------

      double precision function absvec(vec)
      implicit none
      double precision ,dimension(:),intent(in):: vec(:)
      integer i
      absvec=0.0D0
      do i=1,size(vec)
        absvec=absvec+vec(i)**2
      end do
      absvec=sqrt(absvec)
      end function

c---------------------------------------------------------------------

      double precision function matnorm(mat)
      implicit none
      double precision ,dimension(:,:),intent(in):: mat(:,:)
      integer i,j
      matnorm=0.0D0
      do i=1,size(mat,1)
        do j=1,size(mat,2)
          matnorm=matnorm+mat(i,j)**2
        end do     
      end do
      matnorm=sqrt(matnorm)

      end function

c---------------------------------------------------------------------    

        subroutine drawballs(nballs,ndraw,drawn)
	use defs
        implicit none
        ! local
	integer n,nlower,nballs,nrand,ndraw,ball,now(1:8),n1,idraw
	integer, allocatable :: drawn(:),drawnsorted(:)
	double precision tstart,tend,dummyfloat
        !double precision rand
        !intrinsic rand
	character line*24
	logical newball,isdrawn,isopen12
	integer, allocatable :: seed(:)

	inquire(unit=12, opened=isopen12) 
	if (isopen12) then
	  write(12,fsubstart) 'drawballs'
        else
	  print fsubstart,'drawballs'
	end if

        ! check if there are enough balls in the urn to draw from  
	if (ndraw.gt.nballs) then
          goto 1000
	end if

	allocate(drawn(1:ndraw),drawnsorted(1:ndraw))

	nrand=1
	n=1
	do while (n.le.ndraw) 
          !
          ! begin get random I8 number
          !
	  !call date_and_time(values=now)
          !dummyfloat=rand(now(8)*nrand)
          !write(line(1:14),'(e14.8)') dummyfloat
          !read(line(3:10),'(I8)')ball
	  
          ! get seed
          call init_random_seed()
          ! get random number
    	  call random_number(dummyfloat)  
          
          write(line(1:14),'(e14.8)') dummyfloat
          read(line(3:10),'(I8)')ball
          !
          ! end get random I8 number
          !
	  ball=modulo(ball,nballs)
	  if (ball.eq.0) ball=ball+nballs
	  !print*,'drawn ball:',ball
	  newball=.True.
	  do nlower=1,n-1
	    if (drawn(nlower).eq.ball) newball=.False.
	  end do
	  if (newball) then
	    drawn(n)=ball
	    !print*, 'new ball:', drawn(n)
	    n=n+1
	  end if
	  !nrand=nrand+1
	  nrand=nrand+37
	end do	
	
	!print*, 'The numbers of the balls you drew (in this order):'
	!do n=1,ndraw
	!  print*,drawn(n)
	!end do
       
        ! Sort
        idraw=0
        do n=1,nballs
            do n1=1,ndraw
              if (n.eq.drawn(n1)) then
                    idraw=idraw+1
                    drawnsorted(idraw)=n
              end if
            end do
        end do
        drawn=drawnsorted
        
!	!print*,'The numbers of the balls you drew in ascending order:'
!	do n=1,nballs
!	  isdrawn=.False.
!	  do n1=1,ndraw
!	    if (n.eq.drawn(n1)) isdrawn=.True.
!	  end do
!	  !if (isdrawn) print*,n
!	end do

        !print*, 'The numbers of the balls left in the urn are:'
	!do n=1,nballs
	!  isdrawn=.False.
	!  do n1=1,ndraw
	!    if (n.eq.drawn(n1)) isdrawn=.True.
	!  end do
	!  if (.not.isdrawn) print*,n
	!end do

        ! end normally
        if(isopen12) then
              write(12,fsubendext)'drawballs'
        else
              print fsubendext,'drawballs'
        end if
        return

        ! end with an error
1000    nerr=nerr+1
        if(isopen12) then
          write(12,ferrmssg)'You want to draw too many balls'
        else
          print ferrmssg,'You want to draw too many balls'
        end if
        return

        end subroutine

c---------------------------------------------------------------------

        subroutine init_random_seed()
        use iso_fortran_env, only: int64
        implicit none
        integer, allocatable :: seed(:)
        integer :: i, n, un, istat, dt(8), pid
        integer(int64) :: t
        
        call random_seed(size = n)
        allocate(seed(n))
        ! First try if the OS provides a random number generator
        open(newunit=un, file="/dev/urandom", access="stream",          &
     &    form="unformatted",action="read",status="old",iostat=istat)
        if (istat == 0) then
           read(un) seed
           close(un)
        else
           ! Fallback to date and time 
           call date_and_time(values=dt)
           seed(1)=dt(8)
           do i=2,n
             seed(i)=dt(mod(n,8))
           end do
        end if
        call random_seed(put=seed)
        end subroutine init_random_seed

c---------------------------------------------------------------------          

        subroutine octet(place,atoms,vecs,species,nnatom,nncoord,oct)
       !subroutine octet(place,atoms,natoms,vecs,species,nspecies,nnatom)

        use defs

        implicit none
        ! global variables
        double precision, intent(in):: place(1:3),vecs(1:3,1:3)
        type(atom), intent(in):: atoms(:)
        type(element), intent(in):: species(:)
        character(len=*), intent(in):: nnatom
        integer, intent(in):: nncoord
        integer oct
        ! local variables
        integer natoms,nspecies,i,j,numnn,inn
        double precision dist,dumvec(1:3),dumdble,tetcenter(1:3),
     &     tetdist,octcenter(1:3),octdist   
        logical isopen12
        type(atom), allocatable :: nnatoms(:)

        ! write greeting
        inquire(unit=12,opened=isopen12)
        if(isopen12) then
              write(12,fsubstart) "octet"
        else
              print fsubstart,"octet"
        end if

        ! determine number of atoms in cell and number of species
        natoms=size(atoms)
        nspecies=size(species)
        ! write an array with only the atoms of the species that
        ! interest us here (nnatom)
        ! Get number of nnatom atoms
        numnn=0
        do i=1,nspecies
            if(species(i)%name(1:2).eq.nnatom(1:2))
     &         numnn=species(i)%howmany
        end do
        allocate(nnatoms(numnn))
        ! get positions of nnatom species
        inn=0
        do i=1,natoms
            if(atoms(i)%name(1:2).eq.nnatom(1:2))then
               inn=inn+1
               nnatoms(inn)=atoms(i)
            end if
        end do
        ! get the distances
        do i=1,numnn
            ! shift by lattice vectors until closest position is found
            call frac2abs(nnatoms(i)%where,vecs,nnatoms(i)%abswhere)
            dist=absvec(nnatoms(i)%abswhere-place)
            do j=1,3
                do while
     &           (absvec(nnatoms(i)%abswhere-place+vecs(j,1:3)).lt.dist)
                    dist=absvec(nnatoms(i)%abswhere-place+vecs(j,1:3))
                    nnatoms(i)%abswhere=nnatoms(i)%abswhere+vecs(j,1:3)
                end do
                do while
     &           (absvec(nnatoms(i)%abswhere-place-vecs(j,1:3)).lt.dist)
                    dist=absvec(nnatoms(i)%abswhere-place-vecs(j,1:3))
                    nnatoms(i)%abswhere=nnatoms(i)%abswhere-vecs(j,1:3)
                end do
            end do
            nnatoms(i)%distance=dist
        end do
        ! sort according to distance
        do i=1,6       ! look only for the first 6 neighbors
        !do i=1,numnn-1
          do j=i+1,numnn
            if(nnatoms(j)%distance.lt.nnatoms(i)%distance) then
                  dumvec=nnatoms(i)%abswhere
                  dumdble=nnatoms(i)%distance
                  nnatoms(i)%abswhere=nnatoms(j)%abswhere
                  nnatoms(i)%distance=nnatoms(j)%distance
                  nnatoms(j)%abswhere=dumvec
                  nnatoms(j)%distance=dumdble
            end if
          end do
          !print*,i," dist=",nnatoms(i)%distance
        end do
        ! determine if the site is more octa than tetra from the
        ! coordination numbers and by comparing the
        ! distances from the tetra/octa centers
        if(nncoord.eq.6) then
              oct=1
        else
            if (nncoord.eq.4) then
                 oct=0
            else
                 oct=1
                 ! calculate distance from center of tetrahedron and octahedron
                 tetcenter=(nnatoms(1)%abswhere+nnatoms(2)%abswhere
     &              +nnatoms(3)%abswhere+nnatoms(4)%abswhere)/4.0D0
                 octcenter=(nnatoms(1)%abswhere+nnatoms(2)%abswhere
     &              +nnatoms(3)%abswhere+nnatoms(4)%abswhere
     &              +nnatoms(5)%abswhere+nnatoms(6)%abswhere)/6.0D0
                 tetdist=absvec(place-tetcenter)
                 octdist=absvec(place-octcenter)
                 if (tetdist.lt.octdist) oct=0
            end if
        end if
        !print '("distance from tet center: ",F)', tetdist
        !print '("distance from oct center: ",F)', octdist


        !deallocate(nnatoms)

        end subroutine

c---------------------------------------------------------------------

      subroutine getspecies(atoms,species)

      use defs
      implicit none 

      type(atom), intent(in) :: atoms(:)
      type(element), allocatable :: species(:)
      integer i,iatom,iatom2,natoms,nspecies,ispecies
      logical newspecies

      natoms=size(atoms)
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
      !get name mass of each species
      if (allocated(species)) deallocate(species)
      allocate(species(nspecies))
      species%name=" "
      ispecies=1
      do iatom=1,natoms
        newspecies=.true.
        do i=1,iatom-1
          if (atoms(i)%name(1:2).eq.atoms(iatom)%name(1:2))
     &        newspecies=.false.
        end do
        if(newspecies) then
            species(ispecies)%name=atoms(iatom)%name
            species(ispecies)%mass=atoms(iatom)%mass
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


      end subroutine

c---------------------------------------------------------------------

      subroutine get_masses_of_species(species)

      use defs
      implicit none 

      type(element) :: species(:)
      integer i,nspecies,ispecies

      nspecies=size(species)
      do ispecies=1,nspecies
         do i=1,100
           if (species(ispecies)%name(1:2).eq.elements(i)) then 
             species(ispecies)%mass=masses(i)
           end if
         end do ! i
      end do

      end subroutine get_masses_of_species

c---------------------------------------------------------------------

      subroutine get_masses_of_atoms(atoms)

      use defs
      implicit none 

      type(atom) :: atoms(:)
      integer i,natoms,iatom

      natoms=size(atoms)
      do iatom=1,natoms
         do i=1,100
           if (atoms(iatom)%name(1:2).eq.elements(i)) then 
             atoms(iatom)%mass=masses(i)
           end if
         end do ! i
      end do

      end subroutine get_masses_of_atoms

c---------------------------------------------------------------------

      subroutine getenergy(infile,informat,energy)

      use defs
      implicit none
      character(len=*), intent(in) :: infile,informat
      double precision energy
      character line*150
      logical isopen12

      select case(informat)
      case('OUT','out','MBPP','mbpp')
            open(54,file=infile,status='old')
10          read(54,'(A150)',end=11,err=100) line
            if(index(line,'total en').gt.0) then
                  read(line(34:49),*) energy
                  energy=energy*Ryd
            end if
            goto 10

11         continue      
      end select

      ! return normally
      close(54)
      return

      ! return with error
100   close(54)
      nerr=nerr+1
      inquire(unit=12,opened=isopen12)
      if(isopen12) then
            write(12,ferrmssg) "error during read (getenergy)"
      else
            print ferrmssg, "error during read (getenergy)"
      end if
      return

      end subroutine

c---------------------------------------------------------------------

      subroutine vecs2cellpars(vecs,cellpars)
      use defs
      implicit none
      double precision, intent(in) :: vecs(1:3,1:3)
      double precision cellpars(1:6)
      double precision alats(1:3),angles(1:3)
      integer i

      do i=1,3
        alats(i)=sqrt(vecs(i,1)**2+vecs(i,2)**2+vecs(i,3)**2)
      end do
      angles(1)=acos((vecs(2,1)*vecs(3,1)+vecs(2,2)*vecs(3,2)
     &          +vecs(2,3)*vecs(3,3))/(alats(2)*alats(3)))*180.0D0/Pi
      angles(2)=acos((vecs(1,1)*vecs(3,1)+vecs(1,2)*vecs(3,2)
     &          +vecs(1,3)*vecs(3,3))/(alats(1)*alats(3)))*180.0D0/Pi
      angles(3)=acos((vecs(1,1)*vecs(2,1)+vecs(1,2)*vecs(2,2)
     &          +vecs(1,3)*vecs(2,3))/(alats(1)*alats(2)))*180.0D0/Pi
      cellpars(1:3)=alats(1:3)
      cellpars(4:6)=angles(1:3)
      end subroutine

c---------------------------------------------------------------------

      subroutine cellpars2vecs(cellpars,vecs)
      use defs
      implicit none
      double precision, intent(in) :: cellpars(1:6)
      double precision vecs(1:3,1:3)
      double precision aa,bb,cc,alph,beta,gamm

      aa=cellpars(1)
      bb=cellpars(2)
      cc=cellpars(3)
      alph=cellpars(4)
      beta=cellpars(5)
      gamm=cellpars(6)

      vecs=0.0D0
      vecs(1,1)=aa
      vecs(2,1)=bb*cos(gamm*Pi/180.0D0)
      vecs(2,2)=bb*sin(gamm*Pi/180.0D0)
      vecs(3,1)=cc*cos(beta*Pi/180.0D0)
      vecs(3,2)=(cc*bb*cos(alph*Pi/180.0D0)
     &            -vecs(2,1)*vecs(3,1))/vecs(2,2)
      vecs(3,3)=sqrt(cc**2-vecs(3,1)**2-vecs(3,2)**2)

      end subroutine

c---------------------------------------------------------------------      
!
!      subroutine cageZnO(nZnO,nerr,nwarn,ncomm)
!      
!      use defs
!      implicit none
!      
!      integer nZnO
!      integer nerr,nwarn,ncomm
!      ! local variables
!      integer nZn,nO,i,j
!      double precision dZnO,theta,rsphere,coordZn(1:16,1:3),
!     &                 coordO(1:16,1:3),phi
!      
!      print(fsubstart), "cageZnO"
!      select case(nZnO)
!      case(16)
!        nZn=16          ! number of Zn atoms
!        nO=16           ! number of O atoms
!        dZnO=2.0        ! distance between Zn and O
!        
!        ! 1. (uppermost) layer (6 atoms)
!        theta=(1.0D0*Pi/6.0D0)
!        rsphere=dZnO/sin(theta)
!        phi=(4.0D0*Pi/6.0D0)
!        j=0
!        do i=1,3
!          j=j+1
!          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)))
!          coordZn(j,2)=rsphere*sin(theta)*sin(phi*dble(i-1))
!          coordZn(j,3)=rsphere*cos(theta)
!          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)+0.5D0))
!          coordO(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)+0.5D0))
!          coordO(j,3)=rsphere*cos(theta)
!        end do
!        ! 2. layer (6 atoms)
!        theta=(3.0D0*Pi/8.0D0)
!        phi=(4.0D0*Pi/6.0D0)
!        do i=1,3
!          j=j+1
!          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)))
!          coordO(j,2)=rsphere*sin(theta)*sin(phi*dble(i-1))
!          coordO(j,3)=rsphere*cos(theta)
!          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)+0.5D0))
!          coordZn(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)+0.5D0))
!          coordZn(j,3)=rsphere*cos(theta)
!        end do
!        ! 3. layer (8 atoms)
!        theta=(3.0D0*Pi/6.0D0)
!        phi=(2.0D0*Pi/8.0D0)
!        do i=1,2
!          j=j+1
!          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(2.0D0*dble(i-1)
!     &               +0.5D0))
!          coordZn(j,2)=rsphere*sin(theta)*sin(phi*(2.0D0*dble(i-1)
!     &               +0.5D0))
!          coordZn(j,3)=rsphere*cos(theta)
!          coordO(j,1)=rsphere*sin(theta)*cos(phi*(2.0D0*dble(i)-0.5D0))
!          coordO(j,2)=rsphere*sin(theta)*sin(phi*(2.0D0*dble(i)-0.5D0))
!          coordO(j,3)=rsphere*cos(theta)
!        end do 
!        do i=3,4
!          j=j+1
!          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(2.0D0*dble(i)
!     &                -0.5D0))
!          coordZn(j,2)=rsphere*sin(theta)*sin(phi*(2.0D0*dble(i)
!     &                -0.5D0))
!          coordZn(j,3)=rsphere*cos(theta)
!          coordO(j,1)=rsphere*sin(theta)*cos(phi*(2.0D0*dble(i-1)
!     &                +0.5D0))
!          coordO(j,2)=rsphere*sin(theta)*sin(phi*(2.0D0*dble(i-1)
!     &                +0.5D0))
!          coordO(j,3)=rsphere*cos(theta)
!        end do
!        ! 4. layer (6 atoms)
!        theta=(5.0D0*Pi/8.0D0)
!        phi=(4.0D0*Pi/6.0D0)
!        do i=1,3
!          j=j+1
!          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)))
!          coordO(j,2)=rsphere*sin(theta)*sin(phi*dble(i-1))
!          coordO(j,3)=rsphere*cos(theta)
!          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)+0.5D0))
!          coordZn(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)+0.5D0))
!          coordZn(j,3)=rsphere*cos(theta)
!        end do
!        ! 5th (lowest) layer (6 atoms)
!        theta=(5.0D0*Pi/6.0D0)
!        phi=(4.0D0*Pi/6.0D0)
!        do i=1,3
!          j=j+1
!          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)))
!          coordZn(j,2)=rsphere*sin(theta)*sin(phi*dble(i-1))
!          coordZn(j,3)=rsphere*cos(theta)
!          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)+0.5D0))
!          coordO(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)+0.5D0))
!          coordO(j,3)=rsphere*cos(theta)
!        end do
!      case(12)
!        nZn=12          ! number of Zn atoms
!        nO=12           ! number of O atoms
!        dZnO=2.0        ! distance between Zn and O
!        ! 1. (uppermost) layer (4 atoms)
!        theta=(1.0D0*Pi/6.0D0)
!        rsphere=dZnO/(sin(theta)*sqrt(2.0D0))
!        phi=(8.0D0*Pi/8.0D0)
!        j=0
!        do i=1,2
!          j=j+1
!          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)))
!          coordZn(j,2)=rsphere*sin(theta)*sin(phi*dble(i-1))
!          coordZn(j,3)=rsphere*cos(theta)
!          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)+0.5D0))
!          coordO(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)+0.5D0))
!          coordO(j,3)=rsphere*cos(theta)
!        end do
!        ! 2. layer (4 atoms)
!        theta=(2.0D0*Pi/6.0D0)
!        phi=(8.0D0*Pi/8.0D0)
!        do i=1,2
!          j=j+1
!          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)))
!          coordO(j,2)=rsphere*sin(theta)*sin(phi*dble(i-1))
!          coordO(j,3)=rsphere*cos(theta)
!          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)+0.5D0))
!          coordZn(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)+0.5D0))
!          coordZn(j,3)=rsphere*cos(theta)
!        end do
!        ! 3. (middle) layer (8 atoms)
!        theta=(4.0D0*Pi/8.0D0)
!        phi=(8.0D0*Pi/8.0D0)
!        do i=1,2
!          j=j+1
!          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i)+0.625D0))
!          coordO(j,2)=rsphere*sin(theta)*sin(phi*(dble(i)+0.625D0))
!          coordO(j,3)=rsphere*cos(theta)
!          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)+0.125D0))
!          coordZn(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)+0.125D0))
!          coordZn(j,3)=rsphere*cos(theta)
!        end do
!        do i=1,2
!          j=j+1
!          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i)+0.375D0))
!          coordO(j,2)=rsphere*sin(theta)*sin(phi*(dble(i)+0.375D0))
!          coordO(j,3)=rsphere*cos(theta)
!          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)-0.125D0))
!          coordZn(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)-0.125D0))
!          coordZn(j,3)=rsphere*cos(theta)
!        end do
!        ! 4. layer (4 atoms)
!        theta=(4.0D0*Pi/6.0D0)
!        phi=(8.0D0*Pi/8.0D0)
!        do i=1,2
!          j=j+1
!          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)))
!          coordO(j,2)=rsphere*sin(theta)*sin(phi*dble(i-1))
!          coordO(j,3)=rsphere*cos(theta)
!          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)+0.5D0))
!          coordZn(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)+0.5D0))
!          coordZn(j,3)=rsphere*cos(theta)
!        end do
!        ! 5. (lowest) layer (4 atoms)
!        theta=(5.0D0*Pi/6.0D0)
!        phi=(8.0D0*Pi/8.0D0)
!        do i=1,2
!          j=j+1
!          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)))
!          coordZn(j,2)=rsphere*sin(theta)*sin(phi*dble(i-1))
!          coordZn(j,3)=rsphere*cos(theta)
!          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)+0.5D0))
!          coordO(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)+0.5D0))
!          coordO(j,3)=rsphere*cos(theta)
!        end do
!      case(8)
!        nZn=8          ! number of Zn atoms
!        nO=8           ! number of O atoms
!        dZnO=2.0        ! distance between Zn and O
!        ! 1. (uppermost) layer (4 atoms)
!        theta=(1.75D0*Pi/6.0D0)
!        rsphere=dZnO/sin(theta)
!        phi=(8.0D0*Pi/8.0D0)
!        j=0
!        do i=1,2
!          j=j+1
!          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)))
!          coordZn(j,2)=rsphere*sin(theta)*sin(phi*dble(i-1))
!          coordZn(j,3)=rsphere*cos(theta)
!          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)+0.5D0))
!          coordO(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)+0.5D0))
!          coordO(j,3)=rsphere*cos(theta)
!        end do
!c        ! 2. (middle) layer (8 atoms)
!c        theta=(4.0D0*Pi/8.0D0)
!c        phi=(4.0D0*Pi/8.0D0)
!c        do i=1,4
!c          j=j+1
!c          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)))
!c          coordO(j,2)=rsphere*sin(theta)*sin(phi*dble(i-1))
!c          coordO(j,3)=rsphere*cos(theta)
!c          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)+0.5D0))
!c          coordZn(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)+0.5D0))
!c          coordZn(j,3)=rsphere*cos(theta)
!c        end do
!        ! 2. (middle) layer (8 atoms)
!        theta=(4.0D0*Pi/8.0D0)
!        phi=(8.0D0*Pi/8.0D0)
!        do i=1,2
!          j=j+1
!          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i)+0.625D0))
!          coordZn(j,2)=rsphere*sin(theta)*sin(phi*(dble(i)+0.625D0))
!          coordZn(j,3)=rsphere*cos(theta)
!          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)+0.125D0))
!          coordO(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)+0.125D0))
!          coordO(j,3)=rsphere*cos(theta)
!        end do
!        do i=1,2
!          j=j+1
!          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i)+0.375D0))
!          coordZn(j,2)=rsphere*sin(theta)*sin(phi*(dble(i)+0.375D0))
!          coordZn(j,3)=rsphere*cos(theta)
!          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)-0.125D0))
!          coordO(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)-0.125D0))
!          coordO(j,3)=rsphere*cos(theta)
!        end do
!        ! 3th (lowest) layer (4 atoms)
!        theta=(4.25D0*Pi/6.0D0)
!        phi=(8.0D0*Pi/8.0D0)
!        do i=1,2
!          j=j+1
!          coordZn(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)))
!          coordZn(j,2)=rsphere*sin(theta)*sin(phi*dble(i-1))
!          coordZn(j,3)=rsphere*cos(theta)
!          coordO(j,1)=rsphere*sin(theta)*cos(phi*(dble(i-1)+0.5D0))
!          coordO(j,2)=rsphere*sin(theta)*sin(phi*(dble(i-1)+0.5D0))
!          coordO(j,3)=rsphere*cos(theta)
!        end do
!      case default 
!        print(ferrmssg),"Can deal with 8, 12 or 16 ZnO units so far."
!        nerr=nerr+1
!        return
!      end select
!      
!      ! write coordinate file
!      open(51,file="CageZnO.xyz",status="replace")
!      write(51,*) 2*nZnO
!      write(51,*) "(ZnO) cage"
!      do i=1,nZnO
!        !write(51,'("Zn",1x,"core",3(F15.10))')coordZn(i,1:3)
!        write(51,'("Zn",1x,3(F15.10))')coordZn(i,1:3)
!      end do
!      do i=1,nZnO
!        !write(51,'("O ",1x,"core",3(F15.10))')coordO(i,1:3)
!        !write(51,'("O ",1x,"shel",3(F15.10))')coordO(i,1:3)
!        write(51,'("O ",1x,3(F15.10))')coordO(i,1:3)
!      end do
!      close(51)
!
!      print(fsubend)
!      return 
!
!      end subroutine

c---------------------------------------------------------------------

      subroutine sortdata(infile,ncolumns,thecolumn)
     
      use defs       
      implicit none
      character(*), intent(in) :: infile
      integer ncolumns,thecolumn
       
      ! local variables  
      character line*1024,lline*1024,filename*256
      character*1024 :: formatstring
      integer ndatalines,idataline,jdataline,icolumn 
      integer nnonredun
      double precision, allocatable :: dataset(:,:) , dumdata(:)
      double precision, allocatable :: nonredundata(:,:) 
      double precision avdata,mindata,maxdata,stddev
      double precision locmin,locmax
      character (len=1024), allocatable :: lines(:,:),dumlines(:)
      character (len=1024), allocatable :: nonredunlines(:,:),words(:) 
      logical redun

      print fsubstart,"sortdata"
      open(51,file=infile,status="old",err=100)
      rewind(51)
      ! count data lines
      ndatalines=0
10    read(51,*,end=20,err=100) line
      lline=adjustl(line)
      if(lline(1:1).eq."#".or.len_trim(lline).eq.0) goto 10   ! ignore commented and empty lines  
      ndatalines=ndatalines+1
      goto 10

20    continue
      allocate(dataset(ndatalines,ncolumns),dumdata(1:ncolumns))
      allocate(nonredundata(ndatalines,ncolumns))
      allocate(nonredunlines(ndatalines,ncolumns))
      allocate(lines(ndatalines,ncolumns))
      allocate(dumlines(ncolumns))
      rewind (51)

      ! read data
      idataline=1
30    read(51,'(A1024)',end=40,err=100) line(1:1024)
      lline=adjustl(line)
      if(lline(1:1).eq."#".or.len_trim(lline).eq.0) goto 30   ! ignore commented and empty lines  
      !read(line,*) dataset(idataline,1:ncolumns)
      call string2words(lline,words)
      if (.not.ncolumns.eq.size(words)) then
        close(51)
        call error_stop('sortdata: wrong number of columns in file?')
      end if
      !nwords=size(words)
      !print*,'words=',words
      lines(idataline,1:ncolumns)=' '
      lines(idataline,1:size(words))=words(1:size(words))
      !print*, 'lines(',idataline,')=',lines(idataline,1:ncolumns)
      read(words(thecolumn),*) dataset(idataline,thecolumn)
      idataline=idataline+1
      goto 30

40    continue
      close(51)
      ! begin write local extrema to file
      open(52,file="LOC_MIN_MAX",status="replace") 
      write(52,'("# local minima and maxima")')
      do idataline=1,1
          ! maximum 
          if(dataset(idataline+1,thecolumn).lt.
     &      dataset(idataline,thecolumn)) then
            !write(52,*) dataset(idataline,1:ncolumns),"  1 max"
            !write(52,*) lines(idataline,1:ncolumns),"  1 max"
            write(52,*) (trim(lines(idataline,icolumn)),' ',              &
     &               icolumn=1,ncolumns),'  1 max'
          end if
          ! minimum 
          if(dataset(idataline+1,thecolumn).gt.
     &      dataset(idataline,thecolumn)) then
            !write(52,*) dataset(idataline,1:ncolumns)," -1 min"
            !write(52,*) lines(idataline,1:ncolumns)," -1 min"
            write(52,*) (trim(lines(idataline,icolumn)),' ',              &
     &               icolumn=1,ncolumns),' -1 min'
          end if
          ! saddle 
          if(dataset(idataline+1,thecolumn).eq.
     &      dataset(idataline,thecolumn)) then
            !write(52,*) dataset(idataline,1:ncolumns),"  0 saddle"
            !write(52,*) lines(idataline,1:ncolumns),"  0 saddle"
            write(52,*) (trim(lines(idataline,icolumn)),' ',              &
     &               icolumn=1,ncolumns),'   0 saddle'
          end if
      end do
      do idataline=2,ndatalines-1
          ! maximum 
          if(dataset(idataline-1,thecolumn).lt.
     &      dataset(idataline,thecolumn).and.
     &      dataset(idataline+1,thecolumn).lt.
     &      dataset(idataline,thecolumn)) then
            !write(52,*) dataset(idataline,1:ncolumns),"  1 max"
            !write(52,*) lines(idataline,1:ncolumns),"  1 max"
            write(52,*) (trim(lines(idataline,icolumn)),' ',              &
     &               icolumn=1,ncolumns),'  1 max'
          end if
          ! minimum 
          if(dataset(idataline-1,thecolumn).gt.
     &      dataset(idataline,thecolumn).and.
     &      dataset(idataline+1,thecolumn).gt.
     &      dataset(idataline,thecolumn)) then
            !write(52,*) dataset(idataline,1:ncolumns)," -1 min"
            !write(52,*) lines(idataline,1:ncolumns)," -1 min"
            write(52,*) (trim(lines(idataline,icolumn)),' ',              &
     &               icolumn=1,ncolumns),' -1 min'
          end if
          ! saddle 
          if(dataset(idataline-1,thecolumn).eq.
     &      dataset(idataline,thecolumn).or.
     &      dataset(idataline+1,thecolumn).eq.
     &      dataset(idataline,thecolumn)) then
            !write(52,*) dataset(idataline,1:ncolumns),"  0 saddle"
            !write(52,*) lines(idataline,1:ncolumns),"  0 saddle"
            write(52,*) (trim(lines(idataline,icolumn)),' ',              &
     &               icolumn=1,ncolumns),'  0 saddle'
          end if
      end do
      do idataline=ndatalines,ndatalines
          ! maximum 
          if (dataset(idataline-1,thecolumn).lt.
     &      dataset(idataline,thecolumn)) then
            !write(52,*) dataset(idataline,1:ncolumns),"  1 max"
            !write(52,*) lines(idataline,1:ncolumns),"  1 max"
            write(52,*) (trim(lines(idataline,icolumn)),' ',              &
     &               icolumn=1,ncolumns),'  1 max'
          end if
          ! minimum 
          if (dataset(idataline-1,thecolumn).gt.
     &      dataset(idataline,thecolumn)) then
            !write(52,*) dataset(idataline,1:ncolumns)," -1 min"
            !write(52,*) lines(idataline,1:ncolumns)," -1 min"
            write(52,*) (trim(lines(idataline,icolumn)),' ',              &
     &               icolumn=1,ncolumns),' -1 min'
          end if
          ! saddle 
          if (dataset(idataline-1,thecolumn).eq.
     &      dataset(idataline,thecolumn)) then
            !write(52,*) dataset(idataline,1:ncolumns),"  0 saddle"
            !write(52,*) lines(idataline,1:ncolumns),"  0 saddle"
            write(52,*) (trim(lines(idataline,icolumn)),' ',              &
     &               icolumn=1,ncolumns),'  0 saddle'
          end if
      end do
      close(52)  
      ! end write local extrema to file
      !
      ! sort data
      do idataline=1,ndatalines
        do jdataline=1,ndatalines
          if(dataset(jdataline,thecolumn).gt.
     &      dataset(idataline,thecolumn)) then
            !dumdata(1:ncolumns)=dataset(idataline,1:ncolumns)
            dumdata(thecolumn)=dataset(idataline,thecolumn)
            dumlines(1:ncolumns)=lines(idataline,1:ncolumns)
            !dataset(idataline,1:ncolumns)=dataset(jdataline,1:ncolumns)
            dataset(idataline,thecolumn)=dataset(jdataline,thecolumn)
            lines(idataline,1:ncolumns)=lines(jdataline,1:ncolumns)
            !dataset(jdataline,1:ncolumns)=dumdata(1:ncolumns) 
            dataset(jdataline,thecolumn)=dumdata(thecolumn) 
            lines(jdataline,1:ncolumns)=dumlines(1:ncolumns) 
          end if
        end do
      end do
      ! begin average data
      ! and get minimum and maximum of data
      avdata=0.0D0
      mindata=1.0E12
      maxdata=-1.0E12 
      do idataline=1,ndatalines
        avdata=avdata+dataset(idataline,thecolumn)
        if(dataset(idataline,thecolumn).gt.maxdata) 
     &    maxdata=dataset(idataline,thecolumn)
        if(dataset(idataline,thecolumn).lt.mindata) 
     &    mindata=dataset(idataline,thecolumn)
      end do
      avdata=avdata/dble(ndatalines)
      ! end average data
      ! and get minimum and maximum of data
      !
      ! begin get standard deviation
      stddev=0.0d0
      do idataline=1,ndatalines
        stddev=stddev+(dataset(idataline,thecolumn)-avdata)**2
      end do
      stddev=sqrt(stddev/dble(ndatalines))
      ! end get standard deviation
      !
      ! begin get nonredundant data
      nnonredun=0
      do idataline=1,ndatalines
        redun=.false.
        do jdataline=idataline+1,ndatalines
          if(dataset(idataline,thecolumn).eq.dataset(jdataline,           &
     &       thecolumn)) redun=.true.
        end do
        if (.not.redun) then
          nnonredun=nnonredun+1
!          nonredundata(nnonredun,1:ncolumns)                              &
!     &      =dataset(idataline,1:ncolumns)
          nonredunlines(nnonredun,1:ncolumns)                              &
     &      =lines(idataline,1:ncolumns)
        end if
      end do
      ! end get nonredundant data
      !
      ! begin print sorted data
      filename=" "
      filename(1:len_trim(infile))=trim(infile)
      filename(len_trim(infile)+1:len_trim(infile)+7)=".sorted"
      open(51,file=filename,status="replace") 
      do idataline=1,ndatalines
        !write(51,*) dataset(idataline,1:ncolumns)
        !write(51,*) lines(idataline,1:ncolumns)
!        write(51,*) (trim(lines(idataline,icolumn)),' ',                  &
!     &               icolumn=1,ncolumns)
      formatstring=''
      write(formatstring,'("(",I0,"(A,x))")') ncolumns
      write(51,formatstring) (trim(lines(idataline,icolumn)),             &
     &               icolumn=1,ncolumns)
      end do
      write(51,'("# ndata, average, min, max, stddev : ")')
      !write(51,*) "#",ndatalines,avdata,mindata,maxdata
!      write(51,'("# ",I10,3(2x,e20.10))') ndatalines,avdata,mindata,        & 
!     &      maxdata 
      write(51,'("# ",I10,4(2x,F20.10))') ndatalines,avdata,mindata,      &
     &      maxdata, stddev 
      close(51)  
      ! end print sorted data
      !
      ! begin print nonredundant data
      filename=" "
      filename(1:len_trim(infile))=trim(infile)
      filename(len_trim(infile)+1:len_trim(infile)+20)                    &
     &        =".sorted.nonredundant"
      open(51,file=filename,status="replace") 
      do idataline=1,nnonredun
        !write(51,*) nonredundata(idataline,1:ncolumns)
        !write(51,*) nonredunlines(idataline,1:ncolumns)
!        write(51,*) (trim(nonredunlines(idataline,icolumn)),' ',          &
!     &               icolumn=1,ncolumns)
        write(51,formatstring) (trim(nonredunlines(idataline,icolumn)),   &
     &               icolumn=1,ncolumns)
      end do
      close(51)  
      ! end print nonredundant data
      !
      print fsubend
      deallocate(dataset)
      deallocate(lines)
      deallocate(nonredunlines)

      return   
 
100   nerr=nerr+1
      close(51)
      print ferrmssg, "could not read data file."
      return
      end subroutine sortdata

c---------------------------------------------------------------------

      subroutine slidwinav(infile,avlength,ncolumns,thecolumn)
     
      use defs       
      implicit none
      character(*), intent(in) :: infile
      double precision, intent(in) :: avlength 
      integer ncolumns,thecolumn,xcolumn
       
      ! local variables  
      character line*256,lline*256,filename*256,formatstring*128
      integer ndata,nav,nav_half,idata,jdata,jtrue 
      integer jstart,jend
      double precision xdatarange
      double precision, allocatable :: dataset(:,:),avdata(:)

      print fsubstart,"slidwinav"
      if(talk) then
        print '(8x,"infile: ",A)', adjustl(trim(infile))
        print '(8x,"avlength: ",F10.6)', avlength
        print '(8x,"ncolumns: ",I3)', ncolumns
        print '(8x,"icolumn: ",I3)', thecolumn
      end if 
      open(51,file=infile,status="old",err=100)
      rewind(51)
      ! count data lines
      ndata=0
10    read(51,*,end=20,err=100) line
      lline=adjustl(line)
      if(lline(1:1).eq."#".or.len_trim(lline).eq.0) goto 10   ! ignore commented and empty lines  
      ndata=ndata+1
      goto 10

20    continue
      allocate(dataset(ndata,ncolumns),avdata(1:ndata))
      rewind (51)

      ! read data
      idata=1
30    read(51,'(A256)',end=40,err=100) line(1:256)
      lline=adjustl(line)
      if(lline(1:1).eq."#".or.len_trim(lline).eq.0) goto 30   ! ignore commented and empty lines  
      read(line,*) dataset(idata,1:ncolumns)
      idata=idata+1
      goto 30

40    continue
      close(51)
      ! begin sliding window average
      print '(8x,"ndata: ",I8)', ndata
      avdata=0.0D0
      xcolumn=1 ! default: assume first column contains x
      xdatarange=dataset(ndata,xcolumn)-dataset(1,xcolumn)
      print '(8x,"x data range: ",F10.6)', xdatarange 
      nav=ceiling(dble(ndata)*xdatarange*avlength)
      print '(8x,"nav: ",I7)', nav 
      nav_half=ceiling(dble(nav/2.0))
      print '(8x,"nav/2: ",I7)', nav_half 
      do idata=1,ndata
        jstart=idata-nav_half 
        jend=jstart+nav-1  
        do jdata=jstart,jend
          jtrue=jdata
          do while (jtrue.gt.ndata) 
            jtrue=jtrue-ndata
          end do 
          do while (jtrue.lt.1) 
            jtrue=jtrue+ndata
          end do 
          avdata(idata)=avdata(idata)+dataset(jtrue,thecolumn)
        end do
      end do 
      avdata=avdata/dble(nav)                             
      dataset(:,thecolumn)=avdata(:)
      ! print data
      filename=" "
      filename(1:len_trim(infile))=trim(infile)
      filename(len_trim(infile)+1:len_trim(infile)+10)=".SLIDWINAV"
      open(51,file=filename,status="replace") 
      formatstring=''
      write(formatstring,'("(",I0,"(E24.15))")') ncolumns
      !print*,formatstring
      do idata=1,ndata
        write(51,formatstring) dataset(idata,1:ncolumns)
      end do
      close(51)  
      print fsubend
      return   
 
100   nerr=nerr+1
      close(51)
      print ferrmssg, "could not read data file."
      return
      end subroutine slidwinav

c---------------------------------------------------------------------

      subroutine spectrum(infile,informat,
     &            broad,binwidth)

      use defs
      implicit none
       
      character(len=*) infile,informat
      double precision broad,binwidth

      ! local variables
      character line*256,chardum*128
      integer nroots,iroot,ifg,jfg
      double precision integ
      double precision, allocatable :: exens(:),oscis(:),exenfs(:), 
     &      oscifs(:),oscifbrds(:),tmom(:,:),tmomf(:,:),tmomf_temp(:)
      double precision, allocatable :: exent(:),oscit(:),exenft(:), 
     &      oscift(:),oscifbrdt(:),tmom_brd(:,:),tmom_brd_temp(:)
      logical singlets,triplets 

      print fsubstart, "spectrum"

      select case(informat)
      ! NWCHEM !
      case("nwchem","NWCHEM")
        singlets=.false.
        triplets=.false.
        open(51,file=infile,status="old",err=101)
        ! determine number of transitions
10      read (51,'(A256)',err=101,end=101) line
        if(index(line,'Restricted').gt.0.and.index(line,'singlets')
     &     .gt.0) singlets=.true. 
        if(index(line,'Restricted').gt.0.and.index(line,'triplets')
     &     .gt.0) triplets=.true. 
        if(index(line,'Unrestricted').gt.0) triplets=.true. 
        if(index(line,'No. of roots :').gt.0) then
          read(line(29:34),*) nroots
          if (.not. allocated(exens))allocate(exens(nroots),
     &       exent(nroots),oscis(nroots),oscit(nroots),tmom(nroots,3))
          exens=0.0D0
          oscis=0.0D0
          exent=0.0D0
          oscit=0.0D0
          goto 11
        end if   
        goto 10   
11      continue  
        ! read singlet transitions, if available 
        if(singlets) then
          print '(8x,"Attempting to read singlet data.")'
          iroot=1
          tmom=0.0d0
          do while (iroot.le.nroots)
            read(51,'(A256)',err=101,end=101) line
            if(index(line,'Root ').gt.0.and.index(line,'singlet').gt.0) 
     &      then
                !read(line(43:53),*) exens(iroot)
                read(line(index(line,'a.u.')+4:),*) exens(iroot)
            end if
            if(index(line,' Transition Moments    X ').gt.0) 
     &      then
                read(line(index(line,'X')+1:),*) tmom(iroot,1),chardum,   &
     &          tmom(iroot,2),chardum,tmom(iroot,3)
            end if
            if(index(line,'Dipole Oscillator Strength').gt.0) then
              read(line(54:64),*) oscis(iroot)
              iroot=iroot+1 
            end if
          end do
        end if
        ! read triplet transitions, if available 
        if(triplets) then
          print '(8x,"Attempting to read triplet data.")'
          rewind(51)
          iroot=1
          do while (iroot.le.nroots)
            read(51,'(A256)',err=101,end=101) line
!            print*,line
            if(index(line,'Root ').gt.0.and.(index(line,'triplet').gt.0)  &
     &         ) then
                read(line(index(line,'eV')-12:),*,err=101) exent(iroot)
!            end if
!            read(51,'(A256)',err=101,end=101) line
!            read(51,'(A256)',err=101,end=101) line
!            read(51,'(A256)',err=101,end=101) line
!            if(index(line,'Oscillator Strength').gt.0) then
!              if (index(line,'Spin forbidden').le.0) then
!                 read(line(index(line,"trength")+10:),*,err=101)          &
!     &             oscit(iroot)
!              else  
                oscit(iroot)=100.0d0 ! if zero due to spin, fake finite oscillator strength
!              end if
              iroot=iroot+1 
            end if
          end do
        end if
        close(51)
        ! print spectrum as a collection of peaks
        ! singlets
        if(singlets) then
          open(51,file="SPECTRUM.SINGLETS.DAT",status="replace")
          write(51,'("# trans, Excit. energy (eV), Oscill. strength, tra  &
     &nsition moment")')
          do iroot=1,nroots
            write(51,'(I10,5(F15.10))')iroot,exens(iroot),oscis(iroot),   &
     &         tmom(iroot,1:3)
          end do
          close(51)
        end if 
        if(triplets) then
          open(51,file="SPECTRUM.TRIPLETS.DAT",status="replace")
          write(51,'("# trans, Excit. energy (eV), Oscill. strength")')
          do iroot=1,nroots
            write(51,'(I10,2(F15.10))')iroot,exent(iroot),oscit(iroot)
          end do
          close(51)
        end if 
        ! print output
        ! singlets 
        ! print spectrum as a collection of peaks, but now on a fine grid
        if (singlets) then
          call finegrid(exens,oscis,exens(1)-10.0D0*broad,
     &       exens(nroots)+10.0D0*broad,binwidth,
     &       exenfs,oscifs)
          call finegrid(exens,tmom(:,1),exens(1)-10.0D0*broad,
     &       exens(nroots)+10.0D0*broad,binwidth,
     &       exenfs,tmomf_temp)
          allocate(tmomf(size(tmomf_temp),3))
          tmomf(:,1)=tmomf_temp
          call finegrid(exens,tmom(:,2),exens(1)-10.0D0*broad,
     &       exens(nroots)+10.0D0*broad,binwidth,
     &       exenfs,tmomf_temp)
          tmomf(:,2)=tmomf_temp
          call finegrid(exens,tmom(:,3),exens(1)-10.0D0*broad,
     &       exens(nroots)+10.0D0*broad,binwidth,
     &       exenfs,tmomf_temp)
          tmomf(:,3)=tmomf_temp
          open(51,file="SPECTRUM.SINGLETS.FINE.DAT",status="replace")
          write(51,'("# Excit. energy (eV), Oscill. strength")')
          write(51,'("# Excit. energy (eV), Oscill. strength, transition  &
     & moment")')
          do ifg=1,size(exenfs)
            ! integral over spectrum up to this energy:
            integ=0.0d0
            do jfg=1,size(exenfs)
              if (exenfs(jfg).le.exenfs(ifg)) 
     &            integ=integ+oscifs(jfg)
            end do    
            !write(51,'(2(F15.10))') exenfs(ifg),oscifs(ifg)
            write(51,'(6(F15.10))') exenfs(ifg),oscifs(ifg),integ,        &
     &          tmomf(ifg,:)
          end do
          close(51)
        ! print distribution on fine grid with gaussian broadening
          call broaden(exenfs,oscifs,broad,.false.,0.0D0,oscifbrds,       &
     &         .false.)
          call broaden(exenfs,tmomf(:,1),broad,.false.,0.0D0,             &
     &        tmom_brd_temp,.false.)
          allocate(tmom_brd(size(tmom_brd_temp),3))
          tmom_brd(:,1)=tmom_brd_temp
          call broaden(exenfs,tmomf(:,2),broad,.false.,0.0D0,             &
     &        tmom_brd_temp,.false.)
          tmom_brd(:,2)=tmom_brd_temp
          call broaden(exenfs,tmomf(:,3),broad,.false.,0.0D0,             &
     &        tmom_brd_temp,.false.)
          tmom_brd(:,3)=tmom_brd_temp
          open(51,file="SPECTRUM.SINGLETS.BRD.DAT",status="replace")
          write(51,'("# Excit. energy (eV), Oscill. strength with Gaussi
     &an broadening by ",F10.7," transition moment")') broad
          do ifg=1,size(exenfs)
            ! integral over spectrum up to this energy:
            integ=0.0d0
            do jfg=1,size(exenfs)
              integ=integ+oscifbrds(jfg)*binwidth*0.5d0
     &             *( erf((exenfs(ifg)-exenfs(jfg))/broad)
     &               -erf(-exenfs(jfg)/broad))
            end do    
            write(51,'(6(F15.10))') exenfs(ifg),oscifbrds(ifg)*binwidth,
     &           integ,tmom_brd(ifg,:)
          end do
          close(51)
        end if
        ! triplets 
        if (triplets) then
          call finegrid(exent,oscit,exent(1)-10.0D0*broad,
     &       exent(nroots)+10.0D0*broad,binwidth,
     &       exenft,oscift)
          open(51,file="SPECTRUM.TRIPLETS.FINE.DAT",status="replace")
          write(51,'("# Excit. energy (eV), Oscill. strength")')
          do ifg=1,size(exenft)
            ! integral over spectrum up to this energy:
            integ=0.0d0
            do jfg=1,size(exenft)
              if (exenft(jfg).le.exenft(ifg)) 
     &            integ=integ+oscift(jfg)
            end do    
            write(51,'(3(F15.10))') exenft(ifg),oscift(ifg),integ
          end do
          close(51)
        ! print distribution on fine grid with gaussian broadening
          call broaden(exenft,oscift,broad,.false.,0.0D0,oscifbrdt,
     &      .false. )
          open(51,file="SPECTRUM.TRIPLETS.BRD.DAT",status="replace")
          write(51,'("# Excit. energy (eV), Oscill. strength with Gaussi
     &an broadnening by)",F10.7)')broad
          do ifg=1,size(exenft)
            ! integral over spectrum up to this energy:
            integ=0.0d0
            do jfg=1,size(exenft)
              integ=integ+oscifbrdt(jfg)*binwidth*0.5d0
     &             *( erf((exenft(ifg)-exenft(jfg))/broad)
     &               -erf(-exenft(jfg)/broad))
            end do    
            write(51,'(3(F15.10))') exenft(ifg),oscifbrdt(ifg)*binwidth,
     &        integ     
          end do
          close(51)
        end if
      case("abinit","ABINIT")
        singlets=.true.
        triplets=.true.
        open(51,file=infile,status="old",err=101)
        ! determine number of transitions
20      read (51,'(A256)',err=101,end=101) line
        if(index(line,'giving').gt.0.and.index(line,'excitations').gt.0)
     &    then
          read(line(8:13),*) nroots
          if (.not. allocated(exens))allocate(exens(nroots),
     &         exent(nroots),oscis(nroots),oscit(nroots))
          exens=0.0D0
          oscis=0.0D0
          exent=0.0D0
          oscit=0.0D0
          goto 21
        end if   
        goto 20   
21      continue
        ! read singlet transitions 
        do while (index(line,' Oscillator strengths ').le.0.)
          read(51,'(A256)',err=101,end=101) line  
        end do
        read(51,'(A256)',err=101,end=101) line  
        iroot=0
        do while (iroot.lt.nroots)
          read(51,*,err=101,end=101) iroot,exens(iroot),oscis(iroot)
        end do
        exens=exens*hartree
        ! read triplet transitions 
        do while (index(line,' Oscillator strengths ').le.0.)
          read(51,'(A256)',err=101,end=101) line  
        end do
        read(51,'(A256)',err=101,end=101) line  
        iroot=0
        do while (iroot.lt.nroots)
          read(51,*,err=101,end=101) iroot,exent(iroot),oscit(iroot)
        end do
        exent=exent*hartree
        close(51)
        ! print spectrum as a collection of peaks
        ! singlets
        if(singlets) then
          open(51,file="SPECTRUM.SINGLETS.DAT",status="replace")
          write(51,'("# trans, Excit. energy (eV), Oscill. strength")')
          do iroot=1,nroots
            write(51,'(I10,2(F15.10))')iroot,exens(iroot),oscis(iroot)
          end do
          close(51)
        end if 
        if(triplets) then
          open(51,file="SPECTRUM.TRIPLETS.DAT",status="replace")
          write(51,'("# trans, Excit. energy (eV), Oscill. strength")')
          do iroot=1,nroots
            write(51,'(I10,2(F15.10))')iroot,exent(iroot),oscit(iroot)
          end do
          close(51)
        end if 
        ! print output
        ! singlets 
        ! print spectrum as a collection of peaks, but now on a fine grid
        if (singlets) then
          call finegrid(exens,oscis,exens(1)-10.0D0*broad,
     &       exens(nroots)+10.0D0*broad,binwidth,
     &       exenfs,oscifs)
          open(51,file="SPECTRUM.SINGLETS.FINE.DAT",status="replace")
          write(51,'("# Excit. energy (eV), Oscill. strength")')
          do ifg=1,size(exenfs)
            ! integral over spectrum up to this energy:
            integ=0.0d0
            do jfg=1,size(exenfs)
              if (exenfs(jfg).le.exenfs(ifg)) 
     &            integ=integ+oscifs(jfg)
            end do    
            write(51,'(3(F15.10))') exenfs(ifg),oscifs(ifg),
     &         integ
          end do
          close(51)
        ! print distribution on fine grid with gaussian broadening
          call broaden(exenfs,oscifs,broad,.false.,0.0D0,oscifbrds,
     &     .false. )
          open(51,file="SPECTRUM.SINGLETS.BRD.DAT",status="replace")
          write(51,'("# Excit. energy (eV), Oscill. strength with Gaussi
     &an broadnening by)",F10.7)')broad
          do ifg=1,size(exenfs)
            ! integral over spectrum up to this energy:
            integ=0.0d0
            do jfg=1,size(exenfs)
              integ=integ+oscifbrds(jfg)*binwidth*0.5d0
     &             *( erf((exenfs(ifg)-exenfs(jfg))/broad)
     &               -erf(-exenfs(jfg)/broad))
            end do    
            write(51,'((F15.10))') exenfs(ifg),oscifbrds(ifg)*binwidth,
     &         integ
          end do
          close(51)
        end if
        ! triplets 
        if (triplets) then
          call finegrid(exent,oscit,exent(1)-10.0D0*broad,
     &       exent(nroots)+10.0D0*broad,binwidth,
     &       exenft,oscift)
          open(51,file="SPECTRUM.TRIPLETS.FINE.DAT",status="replace")
          write(51,'("# Excit. energy (eV), Oscill. strength")')
          do ifg=1,size(exenft)
            ! integral over spectrum up to this energy:
            integ=0.0d0
            do jfg=1,size(exenft)
              if (exenft(jfg).le.exenft(ifg)) 
     &            integ=integ+oscift(jfg)
            end do    
            write(51,'(2(F15.10))') exenft(ifg),oscift(ifg)
          end do
          close(51)
        ! print distribution on fine grid with gaussian broadening
          call broaden(exenft,oscift,broad,.false.,0.0D0,oscifbrdt,       &
     &         .false.)
          open(51,file="SPECTRUM.TRIPLETS.BRD.DAT",status="replace")
          write(51,'("# Excit. energy (eV), Oscill. strength with Gaussi
     &an broadnening by)",F10.7)')broad
          do ifg=1,size(exenft)
            ! integral over spectrum up to this energy:
            integ=0.0d0
            do jfg=1,size(exenft)
              integ=integ+oscifbrdt(jfg)*binwidth*0.5d0
     &             *( erf((exenft(ifg)-exenft(jfg))/broad)
     &               -erf(-exenft(jfg)/broad))
            end do    
            write(51,'(3(F15.10))') exenft(ifg),oscifbrdt(ifg)*binwidth,
     &         integ
          end do
          close(51)
        end if
      case("fiesta_bse","FIESTA_BSE")
        singlets=.true.
        triplets=.false.
        open(51,file=infile,status="old",err=101)
        ! determine number of transitions
30      read (51,'(A256)',err=101,end=101) line
        if(index(line,'Iterating on (n states) :').gt.0.)
     &    then
          read(line(26:33),*) nroots
          if (.not. allocated(exens))allocate(exens(nroots),
     &         oscis(nroots))
          exens=0.0D0
          oscis=0.0D0
          goto 31
        end if   
        goto 30   
31      continue
        ! read singlet transitions 
        do while (index(line,
     &     'FULL BSE eig.(eV), osc. strength and dipoles:').le.0.)
          read(51,'(A256)',err=101,end=101) line  
        end do
        iroot=0
        do while (iroot.lt.nroots)
          read(51,'(A256)',err=101,end=101) line  
          read(line(7:10),*,err=101,end=101) iroot
          read(line(12:21),*,err=101,end=101) exens(iroot)
          read(line(23:28),*,err=101,end=101) oscis(iroot)
        end do
        close(51)
        ! print spectrum as a collection of peaks
        ! singlets
        if(singlets) then
          open(51,file="SPECTRUM.SINGLETS.DAT",status="replace")
          write(51,'("# trans, Excit. energy (eV), Oscill. strength")')
          do iroot=1,nroots
            write(51,'(I10,2(F15.10))')iroot,exens(iroot),oscis(iroot)
          end do
          close(51)
        end if 
        if(triplets) then
          open(51,file="SPECTRUM.TRIPLETS.DAT",status="replace")
          write(51,'("# trans, Excit. energy (eV), Oscill. strength")')
          do iroot=1,nroots
            write(51,'(I10,2(F15.10))')iroot,exent(iroot),oscit(iroot)
          end do
          close(51)
        end if 
        ! print output
        ! singlets 
        ! print spectrum as a collection of peaks, but now on a fine grid
        if (singlets) then
          call finegrid(exens,oscis,exens(1)-10.0D0*broad,
     &       exens(nroots)+10.0D0*broad,binwidth,
     &       exenfs,oscifs)
          open(51,file="SPECTRUM.SINGLETS.FINE.DAT",status="replace")
          write(51,'("# Excit. energy (eV), Oscill. strength")')
          do ifg=1,size(exenfs)
            ! integral over spectrum up to this energy:
            integ=0.0d0
            do jfg=1,size(exenfs)
              if (exenfs(jfg).le.exenfs(ifg)) 
     &            integ=integ+oscifs(jfg)
            end do    
            write(51,'(3(F15.10))') exenfs(ifg),oscifs(ifg),integ
          end do
          close(51)
        ! print distribution on fine grid with gaussian broadening
          call broaden(exenfs,oscifs,broad,.false.,0.0D0,oscifbrds,       &
     & .false.)
          open(51,file="SPECTRUM.SINGLETS.BRD.DAT",status="replace")
          write(51,'("# Excit. energy (eV), Oscill. strength with Gaussi
     &an broadnening by)",F10.7)')broad
          do ifg=1,size(exenfs)
            ! integral over spectrum up to this energy:
            integ=0.0d0
            do jfg=1,size(exenfs)
              integ=integ+oscifbrds(jfg)*binwidth*0.5d0
     &             *( erf((exenfs(ifg)-exenfs(jfg))/broad)
     &               -erf(-exenfs(jfg)/broad))
            end do    
            !write(51,'(2(F15.10))') exenfs(ifg),oscifbrds(ifg)*binwidth
            write(51,'(3(F15.10))') exenfs(ifg),
     &       oscifbrds(ifg)*binwidth,integ 
          end do
          close(51)
        end if
        ! triplets 
        if (triplets) then
          call finegrid(exent,oscit,exent(1)-10.0D0*broad,
     &       exent(nroots)+10.0D0*broad,binwidth,
     &       exenft,oscift)
          open(51,file="SPECTRUM.TRIPLETS.FINE.DAT",status="replace")
          write(51,'("# Excit. energy (eV), Oscill. strength")')
          do ifg=1,size(exenft)
            ! integral over spectrum up to this energy:
            integ=0.0d0
            do jfg=1,size(exenft)
              if (exenft(jfg).le.exenft(ifg)) 
     &            integ=integ+oscift(jfg)
            end do    
            write(51,'(3(F15.10))') exenft(ifg),oscift(ifg),integ
          end do
          close(51)
        ! print distribution on fine grid with gaussian broadening
          call broaden(exenft,oscift,broad,.false.,0.0D0,oscifbrdt,       &
     &   .false.)
          open(51,file="SPECTRUM.TRIPLETS.BRD.DAT",status="replace")
          write(51,'("# Excit. energy (eV), Oscill. strength with Gaussi
     &an broadnening by)",F10.7)')broad
          do ifg=1,size(exenft)
            ! integral over spectrum up to this energy:
            integ=0.0d0
            do jfg=1,size(exenft)
              integ=integ+oscifbrdt(jfg)*binwidth*0.5d0
     &             *( erf((exenft(ifg)-exenft(jfg))/broad)
     &               -erf(-exenft(jfg)/broad))
            end do    
            write(51,'(3(F15.10))') exenft(ifg),oscifbrdt(ifg)*binwidth,
     &         integ
          end do
          close(51)
        end if
      ! UNKNOWN FORMAT !
      case default
        goto 100
      end select

      print fsubend
      return

! error messages
100   nerr=nerr+1
      print ferrmssg,"Unknown format. See help for known formats"
      return
101   nerr=nerr+1
      print ferrmssg,"A problem occurred while reading the input file."
      close(51)
      return

      end subroutine

c---------------------------------------------------------------------

      subroutine finegrid(cgrid,cval,fgridstart,fgridend,fgridstep,
     &     fgrid,fval)
      ! maps a 1d array defined on a 1d grid to a finer grid
     
      use defs
      implicit none

      double precision, intent(in) :: fgridstart,fgridend
      double precision, intent(in) :: cgrid(:), cval(:)
      double precision, allocatable :: fgrid(:),fval(:)
      double precision, intent(in) :: fgridstep

      ! local variables 
      integer dimcgrid,dimfgrid,icgrid,ifgrid

      print fsubstart,"finegrid"
      
      ! grid dimensions:
      dimcgrid =size(cgrid)
      dimfgrid=int((fgridend-fgridstart)/fgridstep)+2 
      !print*, "fgridstart, fgridend:",fgridstart, fgridend
      !print*, "fine grid step:",fgridstep
      print '(8x,"coarse grid dim.:",I10)',dimcgrid
      print '(8x,"fine grid dim.:",I10)',dimfgrid
      !if(dimfgrid.lt.dimcgrid) goto 100
      if (allocated(fgrid)) deallocate(fgrid)
      if (allocated(fval)) deallocate(fval)
      allocate(fgrid(dimfgrid),fval(dimfgrid))
      ! set up fine grid
      do ifgrid=1,dimfgrid
        fgrid(ifgrid)=fgridstart+fgridstep*dble(ifgrid-1)
      end do 
      !print*,"corase grid:"
      !do icgrid=1,dimcgrid
      !  print*,cgrid(icgrid)
      !end do
      !print*,"fine grid:"
      !do ifgrid=1,dimfgrid
      !  print*,fgrid(ifgrid)
      !end do
      ! map values of coarse grid to fine grid:
      fval=0.0D0   
      do icgrid=1,dimcgrid
        do ifgrid=1,dimfgrid-1
           if(fgrid(ifgrid).le.cgrid(icgrid).and.fgrid(ifgrid+1)
     &        .gt.cgrid(icgrid)) fval(ifgrid)=fval(ifgrid)+cval(icgrid)
        end do
      end do 
      !print*,"corase grid values:"
      !do icgrid=1,dimcgrid
      !  print*,cgrid(icgrid),cval(icgrid)
      !end do
      !print*,"fine grid values:"
      !do ifgrid=1,dimfgrid
      !  print*,fgrid(ifgrid),fval(ifgrid)
      !end do
      print fsubendext,'finegrid'
      return

! Error messages     
100   nerr=nerr+1
      print ferrmssg,"Fine grid dim. < coarse grid dim."
      return

      end subroutine finegrid

c---------------------------------------------------------------------------------------------------

      subroutine broaden(grid,val,brd,periodic,period,brdval,
     &  continuous )
      ! broadens a given distribution. If the distribution is periodic (periodic=.true.), 
      ! the periodic images at +- period are taken into account.
      use defs
      implicit none


      double precision, intent(in) :: grid(:),val(:),brd,period
      logical, intent(in) :: periodic,continuous
      ! local variables
      double precision, allocatable :: brdval(:)
      integer igrid1,igrid2,dimgrid
      double precision intorig,intbrd,sumorig
       
      print fsubstart,"broaden"
      print '(8x,"periodic distribution:",L8)',periodic 
      if(periodic) print '(8x,F15.10)',period
      print '(8x,"broadening by ",F15.10)',brd
      dimgrid=size(grid)
      if (size(val).ne.dimgrid) goto 100
      if (allocated(brdval)) deallocate(brdval)
      allocate(brdval(dimgrid))
      
      ! broaden distribution
      brdval=0.0D0
      do igrid1=1,dimgrid 
        do igrid2=1,dimgrid
          brdval(igrid1)=brdval(igrid1)
     &      +val(igrid2)*exp(-((grid(igrid2)-grid(igrid1))
     &      /brd)**2)/brd          
          if (periodic) then
            brdval(igrid1)=brdval(igrid1)
     &        +val(igrid2)*exp(-((grid(igrid2)-grid(igrid1)+period)
     &        /brd)**2)/brd          
            brdval(igrid1)=brdval(igrid1)
     &        +val(igrid2)*exp(-((grid(igrid2)-grid(igrid1)-period)
     &        /brd)**2)/brd          
          end if
        end do
      end do
      brdval=brdval/sqrt(Pi)
      ! print sum over original and integral over original and broadened distr.
      sumorig=0.0D0
      do igrid1=1,dimgrid
        sumorig=sumorig+val(igrid1)
      end do
      print '(8x,"Sum over original distr.:",F15.5)',sumorig
      intorig=0.0D0
      do igrid1=1,dimgrid-1
        intorig=intorig+val(igrid1)*(grid(igrid1+1)-grid(igrid1))
      end do
      intorig=intorig+val(dimgrid)*(grid(dimgrid)-grid(dimgrid-1))
      print '(8x,"Integ. over original distr.:",F15.5)',intorig
      intbrd=0.0D0
      do igrid1=1,dimgrid-1
        intbrd=intbrd+brdval(igrid1)*(grid(igrid1+1)-grid(igrid1))
      end do
      intbrd=intbrd+brdval(dimgrid)*(grid(dimgrid)-grid(dimgrid-1))
      print '(8x,"Integ. over broadened distr.:",F15.5)',intbrd
      if (continuous) brdval=brdval*intorig/intbrd
      print fsubendext,'broaden'
      return

! Error messages
100   nerr=nerr+1
      print ferrmssg,"Dim. of grid differs from that of field."
      return

      end subroutine broaden
      
c---------------------------------------------------------------------------------------------------

      subroutine broaden_lin(grid,val,brd1,brd2,periodic,period,brdval,
     &  continuous   )
      ! broadens a given distribution. If the distribution is periodic (periodic=.true.), 
      ! the periodic images at +- period are taken into account.
      use defs
      implicit none


      double precision, intent(in) :: grid(:),val(:),brd1,brd2,period
      logical, intent(in) :: periodic,continuous
      ! local variables
      double precision, allocatable :: brdval(:)
      integer igrid1,igrid2,dimgrid
      double precision intorig,intbrd,sumorig,brd
       
      print fsubstart,"broaden_lin"
      print '(8x,"periodic distribution:",L8)',periodic 
      if(periodic) print '(8x,F15.10)',period
      print '(8x,"broadening by ",F15.10," - ",F15.10)',brd1,brd2
      dimgrid=size(grid)
      if (size(val).ne.dimgrid) goto 100
      if (allocated(brdval)) deallocate(brdval)
      allocate(brdval(dimgrid))
      
      ! broaden distribution
      brdval=0.0D0
      do igrid1=1,dimgrid 
        do igrid2=1,dimgrid
          brd=brd1+(grid(igrid2)-grid(1))                                 &
     &       /(grid(dimgrid)-grid(1))*(brd2-brd1)
          if(brd.le.1D-4) then
            if(igrid2==igrid1)brdval(igrid1)=brdval(igrid1)+val(igrid2)
          else
            brdval(igrid1)=brdval(igrid1)
     &        +val(igrid2)*exp(-((grid(igrid2)-grid(igrid1))
     &        /brd)**2)/brd          
            if (periodic) then
              brdval(igrid1)=brdval(igrid1)
     &          +val(igrid2)*exp(-((grid(igrid2)-grid(igrid1)+period)
     &          /brd)**2)/brd          
              brdval(igrid1)=brdval(igrid1)
     &          +val(igrid2)*exp(-((grid(igrid2)-grid(igrid1)-period)
     &          /brd)**2)/brd          
            end if ! periodic
          end if ! brd<=0
        end do ! igrid2
      end do ! igrid1
      brdval=brdval/sqrt(Pi)
      ! print sum over original and integral over original and broadened distr.
      sumorig=0.0D0
      do igrid1=1,dimgrid
        sumorig=sumorig+val(igrid1)
      end do
      print '(8x,"Sum over original distr.:",F15.5)',sumorig
      intorig=0.0D0
      do igrid1=1,dimgrid-1
        intorig=intorig+val(igrid1)*(grid(igrid1+1)-grid(igrid1))
      end do
      intorig=intorig+val(dimgrid)*(grid(dimgrid)-grid(dimgrid-1))
      print '(8x,"Integ. over original distr.:",F15.5)',intorig
      intbrd=0.0D0
      do igrid1=1,dimgrid-1
        intbrd=intbrd+brdval(igrid1)*(grid(igrid1+1)-grid(igrid1))
      end do
      intbrd=intbrd+brdval(dimgrid)*(grid(dimgrid)-grid(dimgrid-1))
      print '(8x,"Integ. over broadened distr.:",F15.5)',intbrd
      if (continuous) brdval=brdval*intorig/intbrd
      print fsubendext,'broaden_lin'
      return

! Error messages
100   nerr=nerr+1
      print ferrmssg,"Dim. of grid differs from that of field."
      return

      end subroutine broaden_lin
      
c---------------------------------------------------------------------------------------------------

      subroutine broaden_quad(grid,val,brd1,brd2,periodic,period,brdval,
     & continuous)
      ! broadens a given distribution, broadening changes
      ! quadratically. If the distribution is periodic (periodic=.true.), 
      ! the periodic images at +- period are taken into account.
      use defs
      implicit none


      double precision, intent(in) :: grid(:),val(:),brd1,brd2,period
      logical, intent(in) :: periodic,continuous
      ! local variables
      double precision, allocatable :: brdval(:)
      integer igrid1,igrid2,dimgrid
      double precision intorig,intbrd,sumorig,brd
       
      print fsubstart,"broaden_quad"
      print '(8x,"periodic distribution:",L8)',periodic 
      if(periodic) print '(8x,F15.10)',period
      print '(8x,"broadening by ",F15.10," - ",F15.10)',brd1,brd2
      dimgrid=size(grid)
      if (size(val).ne.dimgrid) goto 100
      if (allocated(brdval)) deallocate(brdval)
      allocate(brdval(dimgrid))
      
      ! broaden distribution
      brdval=0.0D0
      do igrid1=1,dimgrid 
        do igrid2=1,dimgrid
          brd=sqrt(brd1)+(grid(igrid2)-grid(1))                           &
     &       /(grid(dimgrid)-grid(1))*(sqrt(brd2)-sqrt(brd1))
          brd=brd**2
          if(brd.le.1D-4) then
            if(igrid2==igrid1)brdval(igrid1)=brdval(igrid1)+val(igrid2)
          else
            brdval(igrid1)=brdval(igrid1)
     &        +val(igrid2)*exp(-((grid(igrid2)-grid(igrid1))
     &        /brd)**2)/brd          
            if (periodic) then
              brdval(igrid1)=brdval(igrid1)
     &          +val(igrid2)*exp(-((grid(igrid2)-grid(igrid1)+period)
     &          /brd)**2)/brd          
              brdval(igrid1)=brdval(igrid1)
     &          +val(igrid2)*exp(-((grid(igrid2)-grid(igrid1)-period)
     &          /brd)**2)/brd          
            end if ! periodic
          end if ! brd<=0
        end do ! igrid2
      end do ! igrid1
      brdval=brdval/sqrt(Pi)
      ! print sum over original and integral over original and broadened distr.
      sumorig=0.0D0
      do igrid1=1,dimgrid
        sumorig=sumorig+val(igrid1)
      end do
      print '(8x,"Sum over original distr.:",F15.5)',sumorig
      intorig=0.0D0
      do igrid1=1,dimgrid-1
        intorig=intorig+val(igrid1)*(grid(igrid1+1)-grid(igrid1))
      end do
      intorig=intorig+val(dimgrid)*(grid(dimgrid)-grid(dimgrid-1))
      print '(8x,"Integ. over original distr.:",F15.5)',intorig
      intbrd=0.0D0
      do igrid1=1,dimgrid-1
        intbrd=intbrd+brdval(igrid1)*(grid(igrid1+1)-grid(igrid1))
      end do
      intbrd=intbrd+brdval(dimgrid)*(grid(dimgrid)-grid(dimgrid-1))
      print '(8x,"Integ. over broadened distr.:",F15.5)',intbrd
      if (continuous) brdval=brdval*intorig/intbrd
      print fsubendext,'broaden_quad'
      return

! Error messages
100   nerr=nerr+1
      print ferrmssg,"Dim. of grid differs from that of field."
      return

      end subroutine broaden_quad
      
c---------------------------------------------------------------------

      subroutine broaden_dist(fname,columnx,columnv,periodic,period,      &
     &                        brd,lcontinuous)
      use defs, only : fsubstart,fsubendext,ferrmssg,nerr
      implicit none
      character(len=*), intent(in) :: fname
      integer, intent(in) :: columnx,columnv
      logical, intent(in) :: periodic
      logical, intent(in):: lcontinuous
      double precision, intent(in) :: period
      double precision, intent(in) :: brd
      ! internal variables :
      integer ndatalines,idataline,i
      character line*1024,lline*1024
      double precision, allocatable :: dataset(:,:)
      character*1024, allocatable :: lines(:,:),words(:)
      double precision, allocatable :: broadened_dist(:)
      double precision, allocatable :: equix(:),equiv(:)
      double precision xmin,xmax,binwidth
      character(len=1024) :: fname2
      !
      print fsubstart,'broaden_dist'
      !
      open(51,file=fname,status='old',err=100)
      rewind(51)
      ! count data lines
      ndatalines=0
10    read(51,'(A1024)',end=20,err=100) line
      lline=adjustl(line)
      if(lline(1:1).eq."#".or.len_trim(lline).eq.0) goto 10   ! ignore commented and empty lines  
      ndatalines=ndatalines+1
      goto 10

20    continue
      allocate(dataset(ndatalines,2))!,dumdata(1:ncolumns))
      !allocate(broadened_dist(ndatalines))
      !allocate(lines(ndatalines,2))
      rewind (51)

      ! read data
      idataline=1
30    read(51,'(A1024)',end=40,err=100) line(1:1024)
      lline=adjustl(line)
      if(lline(1:1).eq."#".or.len_trim(lline).eq.0) goto 30   ! ignore commented and empty lines  
      call string2words(lline,words)
      !lines(idataline,1)=words(columnx)
      !lines(idataline,2)=words(columnv)
      read(words(columnx),*) dataset(idataline,1)
      read(words(columnv),*) dataset(idataline,2)
      idataline=idataline+1
      goto 30

40    continue
      !
      close(51)
      ! 
      xmin=minval(dataset(:,1))
      xmax=maxval(dataset(:,1))
      !
      ! begin map to equidistant grid
      !
      binwidth=brd/5.0d0
      if (.not.lcontinuous) then
        call finegrid(dataset(:,1),dataset(:,2),xmin-10.0D0*brd,
     &      xmax+10.0D0*brd,binwidth,equix,equiv)
      ! 
      ! end map to equidistant grid
      !
      else
        allocate(equix(ndatalines),equiv(ndatalines))      
        equix=dataset(:,1)
        equiv=dataset(:,2)
      end if ! .not.lcontinuous
      !
!      call broaden(dataset(:,1),dataset(:,2),brd,periodic,period,         &
!     &     broadened_dist)
      call broaden(equix,equiv,brd,periodic,period,broadened_dist,        &
     & lcontinuous)
      !
      fname2=''
      fname2(1:len(fname))=fname(1:len(fname))
      fname2(len_trim(fname)+1:len_trim(fname)+10)='.BROADENED'
      print '(8x,"Writing broadened distibution to file ",A)',            &
     & trim(adjustl(fname2))
!      print '(8x,I0," points in broadened distribution")',                &
!     &          size(broadened_dist)
      open(51,file=fname2,status='replace')
!      do idataline=1,ndatalines
!        write(51,'(2(F20.10,5x))') dataset(idataline,1),                  &
!     &       broadened_dist(idataline)
!      end do
      do i=1,size(equix)
        write(51,'(2(F20.10,5x))') equix(i),broadened_dist(i)
      end do ! i
      !  
      close(51)
      !
      print fsubendext,'broaden_dist'
      !
      return
      !
      ! Error messages
100   nerr=nerr+1
      close(51)
      print ferrmssg,"Dim. of grid differs from that of field."
      return
      !
      end subroutine broaden_dist

!c---------------------------------------------------------------------
!
!      subroutine broaden_continuous_dist(fname,columnx,columnv,periodic,  &
!     &           period,brd)
!      use defs, only : fsubstart,fsubendext,ferrmssg,nerr
!      implicit none
!      character(len=*), intent(in) :: fname
!      integer, intent(in) :: columnx,columnv
!      logical, intent(in) :: periodic
!      double precision, intent(in) :: period
!      double precision, intent(in) :: brd
!      ! internal variables :
!      integer ndatalines,idataline,i
!      character line*1024,lline*1024
!      double precision, allocatable :: dataset(:,:)
!      character*1024, allocatable :: lines(:,:),words(:)
!      double precision, allocatable :: broadened_dist(:)
!      double precision, allocatable :: equix(:),equiv(:)
!      double precision xmin,xmax,binwidth
!      character(len=1024) :: fname2
!      !
!      print fsubstart,'broaden_continuous_dist'
!      !
!      open(51,file=fname,status='old',err=100)
!      rewind(51)
!      ! count data lines
!      ndatalines=0
!10    read(51,'(A1024)',end=20,err=100) line
!      lline=adjustl(line)
!      if(lline(1:1).eq."#".or.len_trim(lline).eq.0) goto 10   ! ignore commented and empty lines  
!      ndatalines=ndatalines+1
!      goto 10
!
!20    continue
!      allocate(dataset(ndatalines,2))!,dumdata(1:ncolumns))
!      !allocate(broadened_dist(ndatalines))
!      !allocate(lines(ndatalines,2))
!      rewind (51)
!
!      ! read data
!      idataline=1
!30    read(51,'(A1024)',end=40,err=100) line(1:1024)
!      lline=adjustl(line)
!      if(lline(1:1).eq."#".or.len_trim(lline).eq.0) goto 30   ! ignore commented and empty lines  
!      call string2words(lline,words)
!      !lines(idataline,1)=words(columnx)
!      !lines(idataline,2)=words(columnv)
!      read(words(columnx),*) dataset(idataline,1)
!      read(words(columnv),*) dataset(idataline,2)
!      idataline=idataline+1
!      goto 30
!
!40    continue
!      !
!      close(51)
!      ! 
!      ! begin map to equidistant grid
!      !
!      xmin=minval(dataset(:,1))
!      xmax=maxval(dataset(:,1))
!      binwidth=brd/5.0d0
!!      call finegrid(dataset(:,1),dataset(:,2),xmin-10.0D0*brd,
!!     &    xmax+10.0D0*brd,binwidth,equix,equiv)
!      ! 
!      ! end map to equidistant grid
!      !
!      !
!!      call broaden(dataset(:,1),dataset(:,2),brd,periodic,period,         &
!!     &     broadened_dist)
!      call broaden(dataset(:,1),dataset(:,2),brd,periodic,period,         &
!     &broadened_dist,.true.)
!      !
!      fname2=''
!      fname2(1:len(fname))=fname(1:len(fname))
!      fname2(len_trim(fname)+1:len_trim(fname)+10)='.BROADENED'
!      print '(8x,"Writing broadened distibution to file ",A)',            &
!     & trim(adjustl(fname2))
!!      print '(8x,I0," points in broadened distribution")',                &
!!     &          size(broadened_dist)
!      open(51,file=fname2,status='replace')
!!      do idataline=1,ndatalines
!!        write(51,'(2(F20.10,5x))') dataset(idataline,1),                  &
!!     &       broadened_dist(idataline)
!!      end do
!      do i=1,size(dataset,1)
!        write(51,'(2(F20.10,5x))') dataset(i,1),broadened_dist(i)
!      end do ! i
!      !  
!      close(51)
!      !
!      print fsubendext,'broaden_continuous_dist'
!      !
!      return
!      !
!      ! Error messages
!100   nerr=nerr+1
!      close(51)
!      print ferrmssg,"Dim. of grid differs from that of field."
!      return
!      !
!      end subroutine broaden_continuous_dist
!
c---------------------------------------------------------------------

      subroutine broaden_lin_dist(fname,columnx,columnv,periodic,period,  &
     &                        brd1,brd2,lcontinuous)
      use defs, only : fsubstart,fsubendext,ferrmssg,nerr
      implicit none
      character(len=*), intent(in) :: fname
      integer, intent(in) :: columnx,columnv
      logical, intent(in) :: periodic,lcontinuous
      double precision, intent(in) :: period
      double precision, intent(in) :: brd1,brd2
      ! internal variables :
      integer ndatalines,idataline,i
      character line*1024,lline*1024
      double precision, allocatable :: dataset(:,:)
      character*1024, allocatable :: lines(:,:),words(:)
      double precision, allocatable :: broadened_dist(:)
      double precision, allocatable :: equix(:),equiv(:)
      double precision xmin,xmax,binwidth,brd_av
      character(len=1024) :: fname2
      !
      print fsubstart,'broaden_lin_dist'
      !
      open(51,file=fname,status='old',err=100)
      rewind(51)
      ! count data lines
      ndatalines=0
10    read(51,'(A1024)',end=20,err=100) line
      lline=adjustl(line)
      if(lline(1:1).eq."#".or.len_trim(lline).eq.0) goto 10   ! ignore commented and empty lines  
      ndatalines=ndatalines+1
      goto 10

20    continue
      allocate(dataset(ndatalines,2))!,dumdata(1:ncolumns))
      !allocate(broadened_dist(ndatalines))
      !allocate(lines(ndatalines,2))
      rewind (51)

      ! read data
      idataline=1
30    read(51,'(A1024)',end=40,err=100) line(1:1024)
      lline=adjustl(line)
      if(lline(1:1).eq."#".or.len_trim(lline).eq.0) goto 30   ! ignore commented and empty lines  
      call string2words(lline,words)
      !lines(idataline,1)=words(columnx)
      !lines(idataline,2)=words(columnv)
      read(words(columnx),*) dataset(idataline,1)
      read(words(columnv),*) dataset(idataline,2)
      idataline=idataline+1
      goto 30

40    continue
      !
      close(51)
      !
      if (.not.lcontinuous) then
        !
        ! begin map to equidistant grid
        !
        xmin=minval(dataset(:,1))
        xmax=maxval(dataset(:,1))
        brd_av=(brd1+brd2)/2.0d0
        binwidth=brd_av/5.0d0
        call finegrid(dataset(:,1),dataset(:,2),xmin-10.0D0*brd1,
     &    xmax+10.0D0*brd2,binwidth,equix,equiv)
        ! 
        ! end map to equidistant grid
        !
      else
        allocate(equix(ndatalines),equiv(ndatalines))      
        equix=dataset(:,1)
        equiv=dataset(:,2)
      end if ! .not.lcontinous  
      !
!      call broaden(dataset(:,1),dataset(:,2),brd,periodic,period,         &
!     &     broadened_dist)
      call broaden_lin(equix,equiv,brd1,brd2,periodic,period,             &
     &broadened_dist,lcontinuous)
      !
      fname2=''
      fname2(1:len(fname))=fname(1:len(fname))
      fname2(len_trim(fname)+1:len_trim(fname)+10)='.BROADENED'
      print '(8x,"Writing broadened distibution to file ",A)',            &
     & trim(adjustl(fname2))
!      print '(8x,I0," points in broadened distribution")',                &
!     &          size(broadened_dist)
      open(51,file=fname2,status='replace')
!      do idataline=1,ndatalines
!        write(51,'(2(F20.10,5x))') dataset(idataline,1),                  &
!     &       broadened_dist(idataline)
!      end do
      do i=1,size(equix)
        write(51,'(2(F20.10,5x))') equix(i),broadened_dist(i)
      end do ! i
      !  
      close(51)
      !
      print fsubendext,'broaden_lin_dist'
      !
      return
      !
      ! Error messages
100   nerr=nerr+1
      close(51)
      print ferrmssg,"Dim. of grid differs from that of field."
      return
      !
      end subroutine broaden_lin_dist

c---------------------------------------------------------------------
!
!      subroutine broaden_lin_continuous_dist(fname,columnx,columnv,       &
!     &   periodic,period, brd1,brd2)
!      use defs, only : fsubstart,fsubendext,ferrmssg,nerr
!      implicit none
!      character(len=*), intent(in) :: fname
!      integer, intent(in) :: columnx,columnv
!      logical, intent(in) :: periodic
!      double precision, intent(in) :: period
!      double precision, intent(in) :: brd1,brd2
!      ! internal variables :
!      integer ndatalines,idataline,i
!      character line*1024,lline*1024
!      double precision, allocatable :: dataset(:,:)
!      character*1024, allocatable :: lines(:,:),words(:)
!      double precision, allocatable :: broadened_dist(:)
!      double precision, allocatable :: equix(:),equiv(:)
!      double precision xmin,xmax,binwidth,brd_av
!      character(len=1024) :: fname2
!      !
!      print fsubstart,'broaden_lin_continuous_dist'
!      !
!      open(51,file=fname,status='old',err=100)
!      rewind(51)
!      ! count data lines
!      ndatalines=0
!10    read(51,'(A1024)',end=20,err=100) line
!      lline=adjustl(line)
!      if(lline(1:1).eq."#".or.len_trim(lline).eq.0) goto 10   ! ignore commented and empty lines  
!      ndatalines=ndatalines+1
!      goto 10
!
!20    continue
!      allocate(dataset(ndatalines,2))!,dumdata(1:ncolumns))
!      !allocate(broadened_dist(ndatalines))
!      !allocate(lines(ndatalines,2))
!      rewind (51)
!
!      ! read data
!      idataline=1
!30    read(51,'(A1024)',end=40,err=100) line(1:1024)
!      lline=adjustl(line)
!      if(lline(1:1).eq."#".or.len_trim(lline).eq.0) goto 30   ! ignore commented and empty lines  
!      call string2words(lline,words)
!      !lines(idataline,1)=words(columnx)
!      !lines(idataline,2)=words(columnv)
!      read(words(columnx),*) dataset(idataline,1)
!      read(words(columnv),*) dataset(idataline,2)
!      idataline=idataline+1
!      goto 30
!
!40    continue
!      !
!      close(51)
!      ! 
!      ! begin map to equidistant grid
!      !
!!      xmin=minval(dataset(:,1))
!!      xmax=maxval(dataset(:,1))
!!      brd_av=(brd1+brd2)/2.0d0
!!      binwidth=brd_av/5.0d0
!!      call finegrid(dataset(:,1),dataset(:,2),xmin-10.0D0*brd1,
!!     &    xmax+10.0D0*brd2,binwidth,equix,equiv)
!!      ! 
!      ! end map to equidistant grid
!      !
!      !
!!      call broaden(dataset(:,1),dataset(:,2),brd,periodic,period,         &
!!     &     broadened_dist)
!      call broaden_lin(dataset(:,1),dataset(:,2),brd1,brd2,periodic,      &
!     &period,broadened_dist,.true.)
!      !
!      fname2=''
!      fname2(1:len(fname))=fname(1:len(fname))
!      fname2(len_trim(fname)+1:len_trim(fname)+10)='.BROADENED'
!      print '(8x,"Writing broadened distibution to file ",A)',            &
!     & trim(adjustl(fname2))
!!      print '(8x,I0," points in broadened distribution")',                &
!!     &          size(broadened_dist)
!      open(51,file=fname2,status='replace')
!!      do idataline=1,ndatalines
!!        write(51,'(2(F20.10,5x))') dataset(idataline,1),                  &
!!     &       broadened_dist(idataline)
!!      end do
!      do i=1,size(dataset,1)
!        write(51,'(2(F20.10,5x))') dataset(i,1),broadened_dist(i)
!      end do ! i
!      !  
!      close(51)
!      !
!      print fsubendext,'broaden_lin_continuous_dist'
!      !
!      return
!      !
!      ! Error messages
!100   nerr=nerr+1
!      close(51)
!      print ferrmssg,"Dim. of grid differs from that of field."
!      return
!      !
!      end subroutine broaden_lin_continuous_dist
!
c---------------------------------------------------------------------

      subroutine broaden_quad_dist(fname,columnx,columnv,                 &
     &   periodic,period, brd1,brd2,lcontinuous)
      use defs, only : fsubstart,fsubendext,ferrmssg,nerr
      implicit none
      character(len=*), intent(in) :: fname
      integer, intent(in) :: columnx,columnv
      logical, intent(in) :: periodic,lcontinuous
      double precision, intent(in) :: period
      double precision, intent(in) :: brd1,brd2
      ! internal variables :
      integer ndatalines,idataline,i
      character line*1024,lline*1024
      double precision, allocatable :: dataset(:,:)
      character*1024, allocatable :: lines(:,:),words(:)
      double precision, allocatable :: broadened_dist(:)
      double precision, allocatable :: equix(:),equiv(:)
      double precision xmin,xmax,binwidth,brd_av
      character(len=1024) :: fname2
      !
      print fsubstart,'broaden_quad_dist'
      !
      open(51,file=fname,status='old',err=100)
      rewind(51)
      ! count data lines
      ndatalines=0
10    read(51,'(A1024)',end=20,err=100) line
      lline=adjustl(line)
      if(lline(1:1).eq."#".or.len_trim(lline).eq.0) goto 10   ! ignore commented and empty lines  
      ndatalines=ndatalines+1
      goto 10

20    continue
      allocate(dataset(ndatalines,2))!,dumdata(1:ncolumns))
      !allocate(broadened_dist(ndatalines))
      !allocate(lines(ndatalines,2))
      rewind (51)

      ! read data
      idataline=1
30    read(51,'(A1024)',end=40,err=100) line(1:1024)
      lline=adjustl(line)
      if(lline(1:1).eq."#".or.len_trim(lline).eq.0) goto 30   ! ignore commented and empty lines  
      call string2words(lline,words)
      !lines(idataline,1)=words(columnx)
      !lines(idataline,2)=words(columnv)
      read(words(columnx),*) dataset(idataline,1)
      read(words(columnv),*) dataset(idataline,2)
      idataline=idataline+1
      goto 30

40    continue
      !
      close(51)
      ! 
      if (.not.lcontinuous) then
        !
        ! begin map to equidistant grid
        !
        xmin=minval(dataset(:,1))
        xmax=maxval(dataset(:,1))
        brd_av=(brd1+brd2)/2.0d0
        binwidth=brd_av/5.0d0
        call finegrid(dataset(:,1),dataset(:,2),xmin-10.0D0*brd1,
     &    xmax+10.0D0*brd2,binwidth,equix,equiv)
        ! 
        ! end map to equidistant grid
        !
      else
        allocate(equix(ndatalines),equiv(ndatalines))      
        equix=dataset(:,1)
        equiv=dataset(:,2)
      end if ! .not.lcontinous  
      !
      call broaden_quad(equix,equiv,brd1,brd2,periodic,period,            &
     &         broadened_dist,lcontinuous)
      !
      fname2=''
      fname2(1:len(fname))=fname(1:len(fname))
      fname2(len_trim(fname)+1:len_trim(fname)+10)='.BROADENED'
      print '(8x,"Writing broadened distibution to file ",A)',            &
     & trim(adjustl(fname2))
!      print '(8x,I0," points in broadened distribution")',                &
!     &          size(broadened_dist)
      open(51,file=fname2,status='replace')
!      do idataline=1,ndatalines
!        write(51,'(2(F20.10,5x))') dataset(idataline,1),                  &
!     &       broadened_dist(idataline)
!      end do
      do i=1,size(dataset,1)
        write(51,'(2(F20.10,5x))') equix(i),broadened_dist(i)
      end do ! i
      !  
      close(51)
      !
      print fsubendext,'broaden_quad_dist'
      !
      return
      !
      ! Error messages
100   nerr=nerr+1
      close(51)
      print ferrmssg,"Dim. of grid differs from that of field."
      return
      !
      end subroutine broaden_quad_dist

c---------------------------------------------------------------------
           
      subroutine bandgap(infile,informat)
      
      use defs
      implicit none
      
      character (len=*), intent(in) :: infile,informat
      ! internal variables 
      double precision gapa,gapb,homoa,homob,lumoa,lumob,occ,energy
      double precision efermi
      double precision fundgap,dirfundgap,dirfundhomo,dirfundlumo
      double precision directgapa,directgapb
      double precision homoaloc,homobloc,lumoaloc,lumobloc
      double precision homoloc,lumoloc
      double precision :: occtol=1.0E0 ! tolerance below which occupation is treated as zero (useful when there is smearing)
      double precision :: tol_en=5.0E-3 ! tolerance in eV for considering two energ eigenvalues the same 
      double precision, allocatable :: bandeminup(:),bandemaxup(:)
      double precision, allocatable :: bandemindown(:),bandemaxdown(:)
      double precision, allocatable :: kweights(:),levels_up(:,:)
      double precision, allocatable :: kpoints(:,:)
      double precision, allocatable :: levels_down(:,:)
      double precision, allocatable :: occups_up(:,:),occups_down(:,:)
      double precision magmom
      integer nkpoints,ikpoint,ikpoint_down
      integer nbands,nbands_printed,ibandup,ibanddown
      integer nelect
      double precision max_nele_k_up,min_nele_k_up
      double precision max_nele_k_down,min_nele_k_down
      double precision fdum
      integer idum,len_trim_line
      character line*256
      logical spinpol,spinorbit,read_k,qpgw

      print fsubstart, "bandgap"

      select case(informat)
      case("nwchem","NWCHEM")
        spinpol=.false.
        open(51,file=infile,status='old',err=101)
        !
        ! begin read number of printed vectors
        !
        nbands=0
6       read(51,'(A256)',end=7,err=102) line
        if(index(line,'DFT Final Alpha Molecular Orbital Analysis').gt.0  &
     &   .or.index(line,"DFT Final Alpha Molecular Orbital Analysis")
     &   .gt.0.
     &   .or.index(line,"ROHF Final Molecular Orbital Analysis")
     &   .gt.0) then
          nbands=0
        end if
        if(index(line,'Vector').gt.0.and.index(line,'Occ').gt.0) then
          nbands=nbands+1
        end if
        if(index(line,'center of mass').gt.0) then
          goto 7
        end if
        goto 6
        !
7       continue
        !
        print '(8x,I0," states")',nbands 
        rewind(51)
        if (nbands.gt.0) then
          allocate(bandeminup(nbands))
          allocate(bandemaxup(nbands))
        end if
        !
        ! end read number of printed vectors
        !
8       read(51,'(A256)',end=14,err=102) line
        if(index(line,"DFT Final Alpha Molecular Orbital Analysis")
     &     .gt.0) then
          spinpol=.true. 
          occtol=1.0d0/2.0d0
        end if
        if(index(line,"DFT Final Molecular Orbital Analysis").gt.0
     &   .or.index(line,"DFT Final Alpha Molecular Orbital Analysis")
     &   .gt.0.
     &   .or.index(line,"ROHF Final Molecular Orbital Analysis")
     &   .gt.0)
     &  then
          homoa=-1000.0D0
          lumoa=1000.0D0 
          homob=-1000.0D0
          lumob=1000.0D0 
          ibandup=1
          goto 10
        end if
        goto 8
10      read(51,'(A256)',end=14,err=102) line
        if(index(line,'Vector').gt.0.and.index(line,'Occ').gt.0) then
          read(line(19:30),*,err=102) occ
          read(line(35:47),*,err=102) energy
          if(occ.gt.occtol.and.energy.gt.homoa) then
             homoa=energy
          end if
          if(occ.le.occtol.and.energy.lt.lumoa) then
             lumoa=energy
          end if 
          bandeminup(ibandup)=energy*hartree
          bandemaxup(ibandup)=energy*hartree
          ibandup=ibandup+1
        end if 
        if(index(line," center of mass").gt.0) goto 8
        if(index(line," DFT Final Beta Molecular Orbital Analysis")
     &      .gt.0) then 
          if(.not.allocated(bandemindown))allocate(bandemindown(nbands))
          if(.not.allocated(bandemaxdown))allocate(bandemaxdown(nbands))
          ibanddown=1
          goto 12
        end if
        goto 10
12      read(51,'(A256)',end=14,err=102) line
        if(index(line,'Vector').gt.0.and.index(line,'Occ').gt.0) then
          read(line(19:30),*,err=102) occ
          read(line(35:47),*,err=102) energy
          if(occ.gt.occtol.and.energy.gt.homob) then
             homob=energy
          end if
          if(occ.le.occtol.and.energy.lt.lumob) then
             lumob=energy
          end if 
          bandemindown(ibanddown)=energy*hartree
          bandemaxdown(ibanddown)=energy*hartree
          ibanddown=ibanddown+1
        end if 
        if(index(line," center of mass").gt.0) goto 8
        goto 12
14      continue
        close(51)
        homoa=homoa*hartree
        lumoa=lumoa*hartree      
        gapa=lumoa-homoa
        homob=homob*hartree
        lumob=lumob*hartree      
        gapb=lumob-homob
        directgapa=gapa*Hartree
        directgapb=gapb*Hartree
      case ("vasp","VASP","outcar","OUTCAR")
        spinpol=.false.
        open(51,file=infile,status='old',err=101)
18      read(51,'(A256)',end=24,err=102) line
        if(index(line," irreducible k-points:").gt.0) then
          read(line(7:14),*) nkpoints
          print '(8x,I0,x,"kpoints")',nkpoints
          if (allocated(kweights)) deallocate(kweights)
          allocate(kweights(nkpoints))
          if (allocated(kpoints)) deallocate(kpoints)
          allocate(kpoints(nkpoints,1:3))
          read(51,'(A256)',end=24,err=102) line
          read(51,'(A256)',end=24,err=102) line
          read(51,'(A256)',end=24,err=102) line
          do ikpoint=1,nkpoints
            read(51,'(A256)',end=24,err=102) line
            read(line(32:44),*) kweights(ikpoint)
            read(line(1:31),*) kpoints(ikpoint,1:3)
          end do ! ikpoint
        end if ! (index(line," irreducible k-points:").gt.0)
        !
        ! begin read number of electrons
        !
        if(index(line," NELEC").gt.0.and.index(line,'total').gt.0)then
          read(line(index(line,'=')+1:index(line,'.')-1),*) nelect
          print '(8x,I0,x,"electrons")',nelect
        end if
        !  
        ! end read number of electrons
        !
        ! begin read magnetization
        !
        if(index(line,"number of electron").gt.0.and.                     &
     &       index(line,"magnetization").gt.0) then
            if (spinpol)  read(line(50:65),*) magmom
        end if
        !
        ! end read magnetization
        !
        if(index(line,"ISPIN  =      2")
     &     .gt.0) then
          spinpol=.true. 
          occtol=occtol/2.0d0
          print '(8x,"spinpol: ",L1)',spinpol
          if(.not.allocated(bandemindown))                                &
     &       allocate(bandemindown(nbands))
          if(.not.allocated(bandemaxdown))                                &
     &       allocate(bandemaxdown(nbands))
        end if !(index(line,"ISPIN  =      2")
        ! begin noncollinear calculation(OSC)
        if(index(line,"LNONCOLLINEAR =      T")                           &
     &     .gt.0) then
          occtol=occtol/2.0d0
          print '(8x,"spinpol: ",L1)',spinpol
          spinorbit=.true.
          print '(8x,"spinorbit: ",L1)',spinorbit
!          if(.not.allocated(bandemindown))                                &
!     &       allocate(bandemindown(nbands))
!          if(.not.allocated(bandemaxdown))                                &
!     &       allocate(bandemaxdown(nbands))
        end if !(index(line,"LNONCOLLINEAR =      T")
        ! end noncollinear calculation(SOC)
        if (index(line,"number of bands    NBANDS=").gt.0.and.            &
     &      index(line,"parameter").le.0) then
          read(line(index(line,"NBANDS")+7:len(line)),*) nbands
          print '(8x,I0,x,"bands")',nbands
          if(.not.allocated(bandeminup)) allocate(bandeminup(nbands))
          if(.not.allocated(bandemaxup)) allocate(bandemaxup(nbands))
          if (allocated(levels_up)) deallocate(levels_up)
          allocate(levels_up(nkpoints,nbands))
          if (allocated(levels_down)) deallocate(levels_down)
          allocate(levels_down(nkpoints,nbands))
          if (allocated(occups_up)) deallocate(occups_up)
          allocate(occups_up(nkpoints,nbands))
          if (allocated(occups_down)) deallocate(occups_down)
          allocate(occups_down(nkpoints,nbands))
        end if
        if(index(line,"E-fermi :").gt.0) then
          read (line(index(line,'E-fermi')+10:),*) efermi
          homoa=-1000.0D0
          lumoa=1000.0D0 
          homob=-1000.0D0
          lumob=1000.0D0 
          homoaloc=-1000.0D0
          lumoaloc=1000.0D0 
          homobloc=-1000.0D0
          lumobloc=1000.0D0 
          directgapa=2000.0d0
          directgapb=2000.0d0
          bandeminup=1.0d6
          bandemaxup=-1.0d6
          ibandup=1
          if (spinpol) then
            bandemindown=1.0d6
            bandemaxdown=-1.0d6
            ibanddown=1
            ikpoint_down=0
          end if
          ikpoint=0
          goto 20
        end if !(index(line,"E-fermi :").gt.0)
        goto 18
20      read(51,'(A256)',end=24,err=201) line
        if (index(line,'point').gt.0) then
          homoaloc=-1000.0D0
          lumoaloc=1000.0D0 
          ibandup=1
          ikpoint=ikpoint+1
        end if
        if(len_trim(line).gt.1.and.len_trim(line).lt.35.and.              &
     &     index(line,"-------------------------------").le.0.and.        &
     &     index(line,"spin component").le.0) then
          read(line(25:33),*,err=202) occ
          read(line(8:20),*,err=203) energy
          if(occ.gt.occtol.and.energy.gt.homoa) then
             homoa=energy
          end if
          if(occ.gt.occtol.and.energy.gt.homoaloc) then
             homoaloc=energy
          end if
          if(occ.le.occtol.and.energy.lt.lumoa) then
             lumoa=energy
          end if 
          if(occ.le.occtol.and.energy.lt.lumoaloc) then
             lumoaloc=energy
          end if 
          if (energy.gt.bandemaxup(ibandup)) bandemaxup(ibandup)=energy
          if (energy.lt.bandeminup(ibandup)) bandeminup(ibandup)=energy
          !print*,"ikpoint,level=",ikpoint,energy
          levels_up(ikpoint,ibandup)=energy
          occups_up(ikpoint,ibandup)=occ
          ibandup=ibandup+1
        end if 
        if (lumoaloc-homoaloc.lt.directgapa) then
          directgapa=lumoaloc-homoaloc
          dirfundhomo=homoaloc
          dirfundlumo=lumoaloc
        end if
        if(index(line,"-------------------------------").gt.0) goto 18
        if(index(line,"General timing and accounting").gt.0) goto 24
        if(index(line," spin component 2")
     &      .gt.0) then 
          ikpoint_down=0
          goto 22
        end if
        goto 20
22      read(51,'(A256)',end=24,err=221) line
        if (index(line,'point').gt.0) then
          homobloc=-1000.0D0
          lumobloc=1000.0D0 
          ibanddown=1
          ikpoint_down=ikpoint_down+1
        end if
        if(len_trim(line).gt.1.and.len_trim(line).lt.35.and.              &
     &     index(line,"-------------------------------").le.0) then
          read(line(25:33),*,err=222) occ
          read(line(8:20),*,err=223) energy
          if(occ.gt.occtol.and.energy.gt.homob) then
             homob=energy
          end if
          if(occ.gt.occtol.and.energy.gt.homobloc) then
             homobloc=energy
          end if
          if(occ.le.occtol.and.energy.lt.lumob) then
             lumob=energy
          end if 
          if(occ.le.occtol.and.energy.lt.lumobloc) then
             lumobloc=energy
          end if 
          if (energy.gt.bandemaxdown(ibanddown))                          &
     &        bandemaxdown(ibanddown)=energy
          if (energy.lt.bandemindown(ibanddown))                          &
     &        bandemindown(ibanddown)=energy
          levels_down(ikpoint_down,ibanddown)=energy
          occups_down(ikpoint_down,ibanddown)=occ
          ibanddown=ibanddown+1
        end if 
        if (lumobloc-homobloc.lt.directgapb)directgapb=lumobloc-homobloc
        if(index(line,"-------------------------------").gt.0) goto 18
        if(index(line,"General timing and accounting").gt.0) goto 24
        goto 22
24      continue
        close(51)
        gapa=lumoa-homoa
        gapb=lumob-homob
        !fundgap=min(lumoa,lumob)-max(homoa,homob)
        fundgap=gapa
        dirfundgap=directgapa
        !
        if (spinpol) then
         fundgap=min(lumoa,lumob)-max(homoa,homob)
         ! begin get fundamental direct gap (minimum over k-points
         ! (minimum over E(CB)-E(VB)) regradless of spin)
         !
         !dirfundgap=1000.0d0
         !dirfundhomo=1000.0d0
         !dirfundlumo=-1000.0d0
         do ikpoint=1,nkpoints
           !
           ! begin get local VBM and CBM for up and down spins
           !
           homoaloc=maxval(levels_up(ikpoint,:),occups_up(ikpoint,:)       &
     &         .gt.occtol)
           lumoaloc=minval(levels_up(ikpoint,:),occups_up(ikpoint,:)       &
     &         .le.occtol)
           homobloc=maxval(levels_down(ikpoint,:),occups_down(ikpoint,:)   &
     &         .gt.occtol)
           lumobloc=minval(levels_down(ikpoint,:),occups_down(ikpoint,:)   &
     &         .le.occtol)
           homoloc=max(homoaloc,homobloc)
           lumoloc=min(lumoaloc,lumobloc)
!           print*,"ikpoint, homoloc, lumoloc, locdirgap=",                 &
!     &       ikpoint,homoloc,lumoloc,lumoloc-homoloc
           if (lumoloc-homoloc.lt.dirfundgap) then                         &
              dirfundgap=lumoloc-homoloc
              dirfundhomo=homoloc
              dirfundlumo=lumoloc
              !print*, "test"
              !print*,dirfundgap
           end if
           !
           ! end get local VBM and CBM for up and down spins
           !
         end do ! ikpoint
        !
        end if ! spinpol
        !
        ! end get fundamental direct gap (minimum over k-points
        ! (minimum over E(CB)-E(VB)) regradless of spin)
        !
        if (mod(nelect,2).gt.0.and..not.spinpol.and..not.spinorbit)       &
     &      directgapa=gapa
        if (mod(nelect,2).gt.0.and.spinpol.and.abs(magmom).lt.1D-1) then
          directgapa=gapa
          directgapb=gapb
          fundgap=min(gapa,gapb)
          dirfundgap=min(gapa,gapb)
          dirfundhomo=max(homoa,homob)
          dirfundlumo=min(lumoa,lumob)
        end if
        if (any(abs(occups_up(:,:)-occtol).lt.0.5d0*occtol)) then
          gapa=0.0d0
          homoa=efermi
          lumoa=efermi
          directgapa=gapa
          fundgap=0.0d0
          dirfundgap=min(gapa,gapb)
          dirfundhomo=efermi !max(homoa,homob)
          dirfundlumo=efermi !min(lumoa,lumob)
        end if
        if (any(abs(occups_down(:,:)-occtol).lt.0.5d0*occtol)) then
          gapb=0.0d0
          homob=efermi
          lumob=efermi
          directgapb=gapb
          fundgap=0.0d0
          dirfundgap=min(gapa,gapb)
          dirfundhomo=efermi !max(homoa,homob)
          dirfundlumo=efermi !min(lumoa,lumob)
        end if
        max_nele_k_up=0.0d0
        min_nele_k_up=nelect
        if(spinpol) then
          max_nele_k_down=0.0d0
          min_nele_k_down=nelect
        end if
        do ikpoint=1,nkpoints
          if (sum(occups_up(ikpoint,:)).gt.max_nele_k_up)                 &
     &        max_nele_k_up=sum(occups_up(ikpoint,:))
          if (sum(occups_up(ikpoint,:)).lt.min_nele_k_up)                 &
     &        min_nele_k_up=sum(occups_up(ikpoint,:))
          if (spinpol) then
            if (sum(occups_down(ikpoint,:)).gt.max_nele_k_down)           &
     &        max_nele_k_down=sum(occups_down(ikpoint,:))
            if (sum(occups_down(ikpoint,:)).lt.min_nele_k_down)           &
     &        min_nele_k_down=sum(occups_down(ikpoint,:))
          end if
        end do
        ! if the system is a metal, any finite band gap is an artefact
        ! of sparse k-point sampling. Correct for this:
        do ikpoint=1,nkpoints
         if (.not.spinpol.and.(sum(occups_up(ikpoint,:)).lt.              &
     &     nelect-0.5d0*occtol.or.sum(occups_up(ikpoint,:)).gt.           &
     &     nelect+0.5d0*occtol)) then
           gapa=0.0d0
           homoa=efermi
           lumoa=efermi
           directgapa=gapa
           fundgap=0.0d0
           dirfundgap=min(gapa,gapb)
           dirfundhomo=efermi !max(homoa,homob)
           dirfundlumo=efermi !min(lumoa,lumob)
         end if
         if (spinpol.and.max_nele_k_up-min_nele_k_up.gt.                  &
     &       0.5d0*occtol) then
           gapa=0.0d0
           homoa=efermi
           lumoa=efermi
           directgapa=gapa
           fundgap=0.0d0
           dirfundgap=0.0d0
           dirfundhomo=efermi !max(homoa,homob)
           dirfundlumo=efermi !min(lumoa,lumob)
         end if
         if (spinpol.and.min_nele_k_down.lt.max_nele_k_down-0.5d0         &
     &       *occtol) then
           gapb=0.0d0
           homob=efermi
           lumob=efermi
           directgapb=gapb
           fundgap=0.0d0
           dirfundgap=0.0d0
           dirfundhomo=efermi !max(homoa,homob)
           dirfundlumo=efermi !min(lumoa,lumob)
         end if
        end do ! ikpoint 
      case ("vasp_gw","VASP_GW","outcar_gw","OUTCAR_GW")
        spinpol=.false.
        qpgw=.false.
        open(51,file=infile,status='old',err=101)
        read_k=.true.
28      read(51,'(A256)',end=34,err=102) line
        if(index(line," irreducible k-points:").gt.0) then
          if (read_k) read(line(7:14),*) nkpoints
          print '(8x,I0,x,"kpoints")',nkpoints
          if (allocated(kweights)) deallocate(kweights)
          allocate(kweights(nkpoints))
          if (allocated(kpoints)) deallocate(kpoints)
          allocate(kpoints(nkpoints,1:3))
          read(51,'(A256)',end=24,err=102) line
          read(51,'(A256)',end=24,err=102) line
          read(51,'(A256)',end=24,err=102) line
          do ikpoint=1,nkpoints
            read(51,'(A256)',end=34,err=102) line
            read(line(32:44),*) kweights(ikpoint)
            read(line(1:31),*) kpoints(ikpoint,1:3)
          end do ! ikpoint
          read_k=.false.
        end if ! (index(line," irreducible k-points:").gt.0)
        !
        ! begin read number of electrons
        !
        if(index(line," NELEC").gt.0) then
          read(line(index(line,'=')+1:index(line,'.')-1),*) nelect
          print '(8x,I0,x,"electrons")',nelect
        end if
        !  
        ! end read number of electrons
        !
        ! begin check if QPGW
        !
        if (index(line,'QPGW').gt.0) qpgw=.true.
        !
        ! end check if QPGW
        !
        ! begin read magnetization
        !
        if(index(line,"number of electron").gt.0.and.                     &
     &       index(line,"magnetization").gt.0) then
            if (spinpol)  read(line(50:65),*) magmom
        end if
        !
        ! end read magnetization
        !
        if(index(line,"ISPIN  =      2")
     &     .gt.0) then
          spinpol=.true. 
          occtol=occtol/2.0d0
          print '(8x,"spinpol: ",L1)',spinpol
          if(.not.allocated(bandemindown))                                &
     &       allocate(bandemindown(nbands))
          if(.not.allocated(bandemaxdown))                                &
     &       allocate(bandemaxdown(nbands))
        end if !(index(line,"ISPIN  =      2")
        ! begin noncollinear calculation(OSC)
        if(index(line,"LNONCOLLINEAR =      T")                           &
     &     .gt.0) then
          occtol=occtol/2.0d0
          print '(8x,"spinpol: ",L1)',spinpol
          spinorbit=.true.
          print '(8x,"spinorbit: ",L1)',spinorbit
!          if(.not.allocated(bandemindown))                                &
!     &       allocate(bandemindown(nbands))
!          if(.not.allocated(bandemaxdown))                                &
!     &       allocate(bandemaxdown(nbands))
        end if !(index(line,"LNONCOLLINEAR =      T")
        ! end noncollinear calculation(SOC)
        if (index(line,"number of bands    NBANDS=").gt.0.and.            &
     &      index(line,"parameter").le.0) then
          read(line(index(line,"NBANDS")+7:len(line)),*) nbands
          print '(8x,I0,x,"bands")',nbands
          if(.not.allocated(bandeminup)) allocate(bandeminup(nbands))
          if(.not.allocated(bandemaxup)) allocate(bandemaxup(nbands))
          if (allocated(levels_up)) deallocate(levels_up)
          allocate(levels_up(nkpoints,nbands))
          if (allocated(levels_down)) deallocate(levels_down)
          allocate(levels_down(nkpoints,nbands))
          if (allocated(occups_up)) deallocate(occups_up)
          allocate(occups_up(nkpoints,nbands))
          if (allocated(occups_down)) deallocate(occups_down)
          allocate(occups_down(nkpoints,nbands))
        end if
        if(index(line,"E-fermi :").gt.0) then
          read (line(index(line,'E-fermi')+10:),*) efermi
          print '(8x,"E-fermi=",F12.6)',efermi
          goto 30
        end if !(index(line,"E-fermi :").gt.0)
        goto 28
        !
30      read(51,'(A256)',end=34,err=201) line
        if (index(line,'spin component 1').gt.0.and..not.qpgw) then
          ! BEGIN DEBUG
          print '(8x,"going to read spin 1 eigenvalues")'
          ! END DEBUG
          bandeminup=1.0d6
          bandemaxup=-1.0d6
          read(51,'(A256)') line
          do ikpoint=1,nkpoints
            !ibandup=1
            read(51,'(A256)') line
            read(51,'(A256)') line
            read(51,'(A256)') line
            line='not empty'
            len_trim_line=len_trim(line)
            ! read EV for all bands for this kpoint:
            do while (len_trim_line.gt.1)
              read(51,'(A256)') line
              len_trim_line=len_trim(line)
            ! read iband, KS energy, QP-energies   sigma(KS)   V_xc(KS)     V^pw_x(r,r')   Z            occupation Imag(sigma)
            if(len_trim_line.gt.1) read(line,*) ibandup,fdum, energy,     &
     &                   fdum, fdum, fdum, fdum, occ
              if (energy.gt.bandemaxup(ibandup)) bandemaxup(ibandup)      &
     &                =energy
              if (energy.lt.bandeminup(ibandup)) bandeminup(ibandup)      &
     &                =energy
              levels_up(ikpoint,ibandup)=energy
              occups_up(ikpoint,ibandup)=occ
            end do ! while len_trim_line.gt.1
          end do ! ikpoint
          nbands_printed=ibandup ! Not all bands are printed
          ! BEGIN DEBUG
          print '(8x,"spin 1 eigenvalues read")'
          print '(8x,"1st and last EV:",1x,F12.6,1x,F12.6)',              &
     &     levels_up(1,1),levels_up(nkpoints,ibandup)
          ! END DEBUG
          if (.not.spinpol) call print_gap_info()
          if (.not.spinpol) goto 30
        end if  ! (index(line,'spin component 1').gt.0)
        if (index(line,'spin component 2').gt.0.and..not.qpgw) then
          ! BEGIN DEBUG
          print '(8x,"going to read spin 2 eigenvalues")'
          ! END DEBUG
          bandemindown=1.0d6
          bandemaxdown=-1.0d6
          read(51,'(A256)') line
          do ikpoint=1,nkpoints
            !ibanddown=1
            read(51,'(A256)') line
            read(51,'(A256)') line
            read(51,'(A256)') line
            line='not empty'
            len_trim_line=len_trim(line)
            ! read EV for all bands for this kpoint:
            do while (len_trim_line.gt.1)
              read(51,'(A256)') line
              len_trim_line=len_trim(line)
            ! read iband, KS energy, QP-energies   sigma(KS)   V_xc(KS)     V^pw_x(r,r')   Z            occupation Imag(sigma)
            if(len_trim_line.gt.1) read(line,*) ibanddown,fdum, energy,   &
     &         fdum, fdum, fdum, fdum, occ
              if (energy.gt.bandemaxdown(ibanddown))                      &
     &            bandemaxdown(ibanddown)=energy
              if (energy.lt.bandemindown(ibanddown))                      &
     &            bandemindown(ibanddown)=energy
              !print*,"ikpoint,level=",ikpoint,energy
              levels_down(ikpoint,ibanddown)=energy
              occups_down(ikpoint,ibanddown)=occ
            end do ! while len_trim_line.gt.1
          end do ! ikpoint
          ! BEGIN DEBUG
          print '(8x,"spin 2 eigenvalues read")'
          ! END DEBUG
          ! BEGIN DEBUG
          print '(8x,"going to print gap info")'
          ! END DEBUG
          call print_gap_info()
          ! BEGIN DEBUG
          print '(8x,"gap info printed")'
          ! END DEBUG
          goto 30
        end if  ! (index(line,'spin component 2').gt.0)
        !
        ! begin get eigenvalues of QPGW
        !  
        if (qpgw.and.index(line,'QP shifts').gt.0.and.index(line,'iter')  &
     &     .gt.0) then
          bandeminup=1.0d6
          bandemaxup=-1.0d6
          ! BEGIN DEBUG
          !print*,line
          !print '(8x,"going to read spin 1 QPGW eigenvalues")'
          ! END DEBUG
          read(51,'(A256)') line
          if (index(line,'for sc-GW calculations column').gt.0) goto 30
          ! BEGIN DEBUG
          !print*,line
          ! END DEBUG
          read(51,'(A256)') line
          ! BEGIN DEBUG
          !print*,line
          ! END DEBUG
          do ikpoint=1,nkpoints
            read(51,'(A256)') line
            read(51,'(A256)') line
            read(51,'(A256)') line
            line='not empty'
            len_trim_line=len_trim(line)
            ! read EV for all bands for this kpoint:
            do while (len_trim_line.gt.1)
              read(51,'(A256)') line
              len_trim_line=len_trim(line)
              ! BGEIN DEBUG
              !print*,line
              ! END DEBUG
              ! read iband, KS energy, QP-energies ,QP en diag, sigma, Z,  occupation 
              if(len_trim_line.gt.1) read(line,*) ibandup,fdum, energy,   &
     &                   fdum, fdum, fdum, occ
              if (energy.gt.bandemaxup(ibandup)) bandemaxup(ibandup)      &
     &                =energy
              if (energy.lt.bandeminup(ibandup)) bandeminup(ibandup)      &
     &                =energy
              levels_up(ikpoint,ibandup)=energy
              occups_up(ikpoint,ibandup)=occ
            end do ! while len_trim_line.gt.1
            read(51,'(A256)') line
            ! BEGIN DEBUG
            print '(8x,"kpoint, 1st and last QP EV:",1x,I0,1x,F12.6,1x,   &
     &        F12.6)', ikpoint,levels_up(ikpoint,1),                      &
     &        levels_up(ikpoint,ibandup)
            ! END DEBUG
          end do ! ikpoint
          nbands_printed=ibandup ! Not all bands are printed
          ! BEGIN DEBUG
          print '(8x,"spin 1 QPGW eigenvalues read")'
          print '(8x,"1st and last QP EV:",1x,F12.6,1x,F12.6)',           &
     &     levels_up(1,1),levels_up(nkpoints,ibandup)
          ! END DEBUG
          if (.not.spinpol) call print_gap_info()
          if (.not.spinpol) goto 30
          if (spinpol) then
            bandemindown=1.0d6
            bandemaxdown=-1.0d6
            ! BEGIN DEBUG
            print '(8x,"going to read spin 2 QPGW eigenvalues")'
            ! END DEBUG
            do ikpoint=1,nkpoints
              read(51,'(A256)') line
              read(51,'(A256)') line
              read(51,'(A256)') line
              line='not empty'
              len_trim_line=len_trim(line)
              ! read EV for all bands for this kpoint:
              do while (len_trim_line.gt.1)
                read(51,'(A256)') line
                len_trim_line=len_trim(line)
              ! read iband, KS energy, QP-energies   sigma(KS)   V_xc(KS)     V^pw_x(r,r')   Z            occupation Imag(sigma)
              if(len_trim_line.gt.1) read(line,*) ibanddown,fdum,         &
     &          energy, fdum, fdum, fdum, occ
                if (energy.gt.bandemaxdown(ibanddown))                    &
     &              bandemaxdown(ibanddown)=energy
                if (energy.lt.bandemindown(ibanddown))                    &
     &              bandemindown(ibanddown)=energy
                levels_down(ikpoint,ibanddown)=energy
                occups_down(ikpoint,ibanddown)=occ
              end do ! while len_trim_line.gt.1
              read(51,'(A256)') line
              ! BEGIN DEBUG
              print '(8x,"kpoint, 1st and last QP EV:",1x,I0,1x,F12.6,    &
     &          1x,F12.6)', ikpoint,levels_down(ikpoint,1),               &
     &          levels_down(ikpoint,ibanddown)
              ! BEGIN DEBUG
            end do ! ikpoint
            nbands_printed=ibanddown ! Not all bands are printed
            ! BEGIN DEBUG
            print '(8x,"spin 2 QPGW eigenvalues read")'
            print '(8x,"1st and last QP EV:",1x,F12.6,1x,F12.6)',         &
     &       levels_down(1,1),levels_down(nkpoints,ibanddown)
            ! END DEBUG
            call print_gap_info()
          end if !(spinpol)
        end if ! (qpgw.and.index(line,'QP shifts...).gt.0)  
        !
        ! end get eigenvalues of QPGW
        !  
        if(index(line,"-------------------------------").gt.0) goto 30
        if(index(line,"General timing and accounting").gt.0) goto 34
        goto 30
34      continue
        close(51)
        !
        !call print_gap_info()
        !
        ! begin get HOMO, LUMO, gap, direct gap
        !
        homoa=-1000.0d0
        lumoa=1000.0
        homob=-1000.0d0
        lumob=1000.0
        gapa=2000.0d0
        gapb=2000.0d0
        directgapa=2000.0d0
        directgapb=2000.0d0
        fundgap=2000.0d0
        dirfundgap=2000.0d0
        !
        do ikpoint=1,nkpoints
          !
          ! begin get local VBM and CBM for up and down spins
          !
          homoaloc=maxval(levels_up(ikpoint,1:nbands_printed),            &
     &            occups_up(ikpoint,1:nbands_printed) .gt.occtol)
          lumoaloc=minval(levels_up(ikpoint,1:nbands_printed),            &
     &            occups_up(ikpoint,1:nbands_printed) .le.occtol)
          if (homoaloc.gt.homoa) homoa=homoaloc
          if (lumoaloc.lt.lumoa) lumoa=lumoaloc
          if (lumoaloc-homoaloc.lt.directgapa) then  
             directgapa=lumoaloc-homoaloc
          end if
          if(directgapa.lt.dirfundgap) then
             dirfundgap=directgapa
             dirfundhomo=homoaloc
             dirfundlumo=lumoaloc
          end if
          if (spinpol) then
            homobloc=maxval(levels_down(ikpoint,1:nbands_printed),        &
     &              occups_down(ikpoint,1:nbands_printed).gt.occtol)
            lumobloc=minval(levels_down(ikpoint,1:nbands_printed),        &
     &              occups_down(ikpoint,1:nbands_printed).le.occtol)
            if (homobloc.gt.homob) homob=homobloc
            if (lumobloc.lt.lumob) lumob=lumobloc
            if (lumobloc-homobloc.lt.directgapb) then  
                directgapb=lumobloc-homobloc
            end if
            if(directgapb.lt.dirfundgap) then
               dirfundgap=directgapb
               dirfundhomo=homobloc
               dirfundlumo=lumobloc
            end if
          end if ! spinpol
          !
          ! end get local VBM and CBM for up and down spins
          !
        end do ! ikpoint
        gapa=lumoa-homoa
        fundgap=gapa
        if (spinpol) then
           gapb=lumob-homob
           fundgap=min(gapa,gapb)
        end if
        !
        ! end get HOMO, LUMO, gap, direct gap
        !
        ! begin set gap of metals to zero 
        !
        if (mod(nelect,2).gt.0.and..not.spinpol.and..not.spinorbit)       &
     &      directgapa=gapa
        if (mod(nelect,2).gt.0.and.spinpol.and.abs(magmom).lt.1D-1) then
          directgapa=gapa
          directgapb=gapb
          fundgap=min(gapa,gapb)
          dirfundgap=min(gapa,gapb)
          dirfundhomo=max(homoa,homob)
          dirfundlumo=min(lumoa,lumob)
        end if
        if (any(abs(occups_up(:,:)-occtol).lt.0.5d0*occtol)) then
          gapa=0.0d0
          homoa=efermi
          lumoa=efermi
          directgapa=gapa
          fundgap=0.0d0
          dirfundgap=min(gapa,gapb)
          dirfundhomo=efermi !max(homoa,homob)
          dirfundlumo=efermi !min(lumoa,lumob)
        end if
        if (any(abs(occups_down(:,:)-occtol).lt.0.5d0*occtol)) then
          gapb=0.0d0
          homob=efermi
          lumob=efermi
          directgapb=gapb
          fundgap=0.0d0
          dirfundgap=min(gapa,gapb)
          dirfundhomo=efermi !max(homoa,homob)
          dirfundlumo=efermi !min(lumoa,lumob)
        end if
        max_nele_k_up=0.0d0
        min_nele_k_up=nelect
        if(spinpol) then
          max_nele_k_down=0.0d0
          min_nele_k_down=nelect
        end if
        do ikpoint=1,nkpoints
          if (sum(occups_up(ikpoint,1:nbands_printed)).gt.max_nele_k_up)  &
     &        max_nele_k_up=sum(occups_up(ikpoint,1:nbands_printed))
          if (sum(occups_up(ikpoint,1:nbands_printed)).lt.min_nele_k_up)  &
     &        min_nele_k_up=sum(occups_up(ikpoint,1:nbands_printed))
          if (spinpol) then
            if (sum(occups_down(ikpoint,1:nbands_printed)).gt.            &
     &          max_nele_k_down) max_nele_k_down=sum(occups_down(         &
     &          ikpoint,1:nbands_printed))
            if (sum(occups_down(ikpoint,1:nbands_printed)).lt.            &
     &          min_nele_k_down) min_nele_k_down=sum(occups_down(         &
     &          ikpoint,1:nbands_printed))
          end if
        end do
        ! if the system is a metal, any finite band gap is an artefact
        ! of sparse k-point sampling. Correct for this:
        do ikpoint=1,nkpoints
         if (.not.spinpol.and.(sum(occups_up(ikpoint,1:nbands_printed))   &
     &       .lt.nelect-0.5d0*occtol.or.sum(occups_up(ikpoint,            &
     &       1:nbands_printed)).gt.nelect+0.5d0*occtol)) then
           gapa=0.0d0
           homoa=efermi
           lumoa=efermi
           directgapa=gapa
           fundgap=0.0d0
           dirfundgap=min(gapa,gapb)
           dirfundhomo=efermi !max(homoa,homob)
           dirfundlumo=efermi !min(lumoa,lumob)
         end if
         if (spinpol.and.max_nele_k_up-min_nele_k_up.gt.                  &
     &       0.5d0*occtol) then
           gapa=0.0d0
           homoa=efermi
           lumoa=efermi
           directgapa=gapa
           fundgap=0.0d0
           dirfundgap=0.0d0
           dirfundhomo=efermi !max(homoa,homob)
           dirfundlumo=efermi !min(lumoa,lumob)
         end if
         if (spinpol.and.min_nele_k_down.lt.max_nele_k_down-0.5d0         &
     &       *occtol) then
           gapb=0.0d0
           homob=efermi
           lumob=efermi
           directgapb=gapb
           fundgap=0.0d0
           dirfundgap=0.0d0
           dirfundhomo=efermi !max(homoa,homob)
           dirfundlumo=efermi !min(lumoa,lumob)
         end if
        end do ! ikpoint 
        !
        ! end set gap of metals to zero 
        !
      case default
        goto 100
      end select ! case(informat)
      ! 
      print '(8x,"occtol: ",F10.6)',occtol
      print '(8x,"tol_en: ",F10.6, " eV")',tol_en
      print '(8x,"")'
      print '(8x,"alpha homo in Hartree:",F10.6)',homoa/Hartree
      print '(8x,"alpha homo in eV:",F10.6)',homoa
      print '(8x,"")'
      !
      ! begin find kpoint(s) with VBM 
      !
      if (allocated(levels_up)) then
        do ibandup=1,nbands
         do ikpoint=1,nkpoints
          if (abs(levels_up(ikpoint,ibandup)-homoa).lt.tol_en             &
     &        .and.occups_up(ikpoint,ibandup).gt.0.5*occtol) then
            print '(8x,"alpha homo with eigenvalue ",F12.6,               &
     &            " at kpoint # ",I5," in band # ",I5,                    &
     &            " with fractional coordinates ", 3(F12.7))',            &
     &           levels_up(ikpoint,ibandup),ikpoint,ibandup,              &
     &           kpoints(ikpoint,1:3)
          end if           
         end do ! ikpoint
        end do ! ibandup
      end if ! allocated(levelsup)
      !
      ! end find kpoint(s) with VBM 
      !
      print '(8x,"")'
      print '(8x,"alpha lumo in Hartree:",F10.6)',lumoa/Hartree
      print '(8x,"alpha lumo in eV:",F10.6)',lumoa
      print '(8x,"")'
      !
      ! begin find kpoint(s) with CBM 
      !
      if (allocated(levels_up)) then
        do ibandup=1,nbands
         do ikpoint=1,nkpoints
          if (abs(levels_up(ikpoint,ibandup)-lumoa).lt.tol_en             &
     &        .and.occups_up(ikpoint,ibandup).lt.0.5*occtol) then
            print '(8x,"alpha lumo with eigenvalue ",F12.6,               &
     &            " at kpoint # ",I5," in band # ",I5,                    &
     &            " with fractional coordinates ", 3(F12.7))',            &
     &           levels_up(ikpoint,ibandup),ikpoint,ibandup,              &
     &           kpoints(ikpoint,1:3)
          end if           
         end do ! ikpoint
        end do ! ibandup
      end if ! allocated(levelsup)
      print '(8x,"")'
      !
      ! end find kpoint(s) with CBM 
      !
      print '(8x,"alpha bandgap in eV:",F10.6)',gapa
      print '(8x,"alpha direct bandgap in eV:",F10.6)',directgapa
      print '(8x,"")'
      !
      ! begin find kpoint(s) with direct band gap 
      !
      if (allocated(levels_up)) then
        do ibandup=1,nbands
         do ikpoint=1,nkpoints
         if (ibandup.lt.nbands) then
          if (abs(levels_up(ikpoint,ibandup+1)-levels_up(ikpoint,         &
     &     ibandup)-directgapa).lt.2.0d0*tol_en.and.directgapa.gt.tol_en  &
     &        .and.occups_up(ikpoint,ibandup+1).lt.0.5*occtol.and.        &
     &        occups_up(ikpoint,ibandup).gt.0.5*occtol) then
            print '(8x,"alpha direct gap of ",F12.6," eV at kpoint # ",   &
     &      I5," between bands # ",I5," and ",I5,                         &
     &      " with fractional coordinates ",3(F12.7))',                   &
     &      levels_up(ikpoint,ibandup+1)-levels_up(ikpoint,ibandup),      &
     &      ikpoint,ibandup,ibandup+1,kpoints(ikpoint,1:3)
          end if
         end if          
         end do ! ikpoint
        end do ! ibandup
      end if ! allocated(levelsup)
      print '(8x,"")'
      !
      ! end find kpoint(s) with direct band gap 
      !
      open(51,file='BANDEMINMAXUP.DAT',status='replace')
      do ibandup=1,nbands
        write(51,'(I10,5x,2(F20.8,x))') ibandup,bandeminup(ibandup),      &
     &       bandemaxup(ibandup)
      end do ! ibandup
      close(51)
      !
      if (allocated(levels_up)) then
        open(51,file='LEVELSUP.DAT',status='replace')
        write(51,'("# kpoint, band, eigenvalue, weight, occupation")')
        do ibandup=1,nbands
          do ikpoint=1,nkpoints
            write(51,'(2(I10,2x),5x,3(F20.8,x))') ikpoint,ibandup,        &
     &      levels_up(ikpoint,ibandup), kweights(ikpoint),                &
     &         occups_up(ikpoint,ibandup)
          end do ! ikpoint
        end do ! ibandup
        close(51)
      end if ! (allocated(levels_up))
      !  
      if(spinpol) then
        !
        print '(8x,"beta homo in Hartree:",F10.6)',homob/Hartree
        print '(8x,"beta homo in eV:",F10.6)',homob
        print '(8x,"")'
        !
        ! begin find kpoint(s) with VBM 
        !
        if (allocated(levels_down)) then
          do ibanddown=1,nbands
           do ikpoint=1,nkpoints
            if (abs(levels_down(ikpoint,ibanddown)-homob).lt.tol_en       &
     &        .and.occups_down(ikpoint,ibanddown).gt.0.5*occtol) then
              print '(8x,"beta homo with eigenvalue ",F12.6,              &
     &        " at kpoint # ",I5," in band # ",I5,                        &
     &        " with fractional coordinates ",3(F12.7))',                 &
     &        levels_down(ikpoint,ibanddown),ikpoint,ibanddown,           &
     &        kpoints(ikpoint,1:3)
            end if           
           end do ! ibanddown
          end do ! ikpoint
        end if ! (allocated(levels_down))
        print '(8x,"")'
        !
        ! end find kpoint(s) with VBM 
        !
        print '(8x,"beta lumo in Hartree:",F10.6)',lumob/Hartree
        print '(8x,"beta lumo in eV:",F10.6)',lumob
        print '(8x,"")'
        !
        ! begin find kpoint(s) with CBM 
        !
        if (allocated(levels_down)) then
          do ibanddown=1,nbands
           do ikpoint=1,nkpoints
            if (abs(levels_down(ikpoint,ibanddown)-lumob).lt.tol_en       &
     &        .and.occups_down(ikpoint,ibanddown).lt.0.5*occtol) then
              print '(8x,"beta lumo with eigenvalue ",F12.6,              &
     &              " at kpoint # ",I5," in band # ",I5,                  &
     &         " with fractional coordinates ",3(F12.7))',                &
     &       levels_down(ikpoint,ibanddown),ikpoint,ibanddown,            &
     &       kpoints(ikpoint,1:3)
            end if           
           end do ! ibanddown
          end do ! ikpoint
        end if ! (allocated(levels_down))
        !
        ! end find kpoint(s) with CBM 
        !
        print '(8x,"")'
        print '(8x,"beta bandgap in eV:",F10.6)',gapb
        print '(8x,"beta direct bandgap in eV:",F10.6)',directgapb
        print '(8x,"")'
        !
        ! begin find kpoint(s) with direct band gap 
        !
        if (allocated(levels_down)) then
          do ibanddown=1,nbands
           do ikpoint=1,nkpoints
           if (ibanddown.lt.nbands) then
            if(abs(levels_down(ikpoint,ibanddown+1)-levels_down(ikpoint,  &
     &       ibanddown)-directgapb).lt.2.0d0*tol_en.and.directgapb.gt.    &
     &       tol_en.and.occups_down(ikpoint,ibanddown+1).lt.0.5*occtol    &
     &       .and. occups_down(ikpoint,ibanddown).gt.0.5*occtol) then
              print '(8x,"beta direct gap of ",F12.6," eV at kpoint # ",  &
     &        I5," between bandss # ",I5," and ",I5,                      &
     &        " with fractional coordinates ",3(F12.7))',                 &
     &        levels_down(ikpoint,ibanddown+1)-levels_down(ikpoint,       &
     &        ibanddown),ikpoint,ibanddown,ibanddown+1,                   &
     &        kpoints(ikpoint,1:3)
            end if
           end if          
           end do ! ibanddown
          end do ! ikpoint
        end if ! (allocated(levels_down))
        print '(8x,"")'
        !
        ! end find kpoint(s) with direct band gap 
        !
        open(51,file='BANDEMINMAXDOWN.DAT',status='replace')
        do ibanddown=1,nbands
          write(51,'(I10,5x,2(F20.8,x))') ibanddown,                      &
     &      bandemindown(ibanddown),bandemaxdown(ibanddown)
        end do ! ibanddown
        close(51)
        !
        if (allocated(levels_down)) then
          open(51,file='LEVELSDOWN.DAT',status='replace')
          do ibanddown=1,nbands
            do ikpoint=1,nkpoints
              write(51,'(2(I10,2x),5x,3(F20.8,x))') ikpoint,ibanddown,    &
     &        levels_down(ikpoint,ibanddown), kweights(ikpoint),          &
     &           occups_down(ikpoint,ibanddown)
            end do ! ikpoint
          end do ! ibanddown
          close(51)
        end if ! (allocated(levels_down))
        !  
      end if ! (spinpol)
      !
      print '(8x,"fundamental bandgap in eV:",F10.6)',fundgap
      print '(8x,"fundamental direct bandgap in eV:",F10.6)',dirfundgap
      print '(8x,"direct fundamental HOMO in eV:",F10.6)',dirfundhomo
      print '(8x,"direct fundamental LUMO in eV:",F10.6)',dirfundlumo
      !
      ! begin write VASP EIGENVAL file
      !
      if (allocated(kpoints)) then
        open(51,file='COFIMA.EIGENVAL',status='replace')
        if (.not.spinpol) write(51,*) "header header header 1"
        if (spinpol) write(51,*) "header header header 2"
        write(51,*) "header"
        write(51,*) "header"
        write(51,*) "header"
        write(51,*) "unknown system"
        write(51,*) nelect, nkpoints,nbands 
        do ikpoint=1,nkpoints
          write(51,*) " "
          write(51,'(3(F16.12,1x))') kpoints(ikpoint,1:3)
          do ibandup=1,nbands
            if (spinpol) then
              write(51,'(1(I10,2x),5x,2(F20.8,x))') ibandup,              &
     &      levels_up(ikpoint,ibandup), levels_down(ikpoint,ibandup)    
            else
              write(51,'(1(I10,2x),5x,1(F20.8,x))') ibandup,              &
     &      levels_up(ikpoint,ibandup)    
            end if
          end do ! ibandup
        end do ! ikpoint
        close(51)
        !
        ! end write VASP EIGENVAL file
      end if ! allocated kpoints
      !
      !
      print fsubendext,"bandgap"
      return
      !
100   nerr=nerr+1
      print ferrmssg,"Unsupported file format"
      return
101   nerr=nerr+1
      print ferrmssg,"File could not be opened"
      close(51)
      return
102   nerr=nerr+1
      print ferrmssg,"Error during reading"
      close(51)
      return
201   nerr=nerr+1
      print ferrmssg,"201: Error during reading"
      close(51)
      return
202   nerr=nerr+1
      print ferrmssg,"202: Error during reading"
      close(51)
      return
203   nerr=nerr+1
      print ferrmssg,"203: Error during reading"
      close(51)
      return
221   nerr=nerr+1
      print ferrmssg,"221: Error during reading"
      close(51)
      return
222   nerr=nerr+1
      print ferrmssg,"222: Error during reading"
      close(51)
      return
223   nerr=nerr+1
      print ferrmssg,"223: Error during reading"
      close(51)
      return
      !
c---------------------------------------------------------------------
      !
      contains 
      !
!-----------------------------------------------------------------------

      subroutine print_gap_info()
        !
        ! begin get HOMO, LUMO, gap, direct gap
        !
        homoa=-1000.0d0
        lumoa=1000.0
        homob=-1000.0d0
        lumob=1000.0
        gapa=2000.0d0
        gapb=2000.0d0
        directgapa=2000.0d0
        directgapb=2000.0d0
        fundgap=2000.0d0
        dirfundgap=2000.0d0
        !
        ! BEGIN DEBUG
        !
        !print '(8x,"gaps initialized")'
        !
        ! END DEBUG
        !
        do ikpoint=1,nkpoints
          !
          ! begin get local VBM and CBM for up and down spins
          !
          homoaloc=maxval(levels_up(ikpoint,1:nbands_printed),            &
     &            occups_up(ikpoint,1:nbands_printed) .gt.occtol)
          lumoaloc=minval(levels_up(ikpoint,1:nbands_printed),            &
     &            occups_up(ikpoint,1:nbands_printed) .le.occtol)
          if (homoaloc.gt.homoa) homoa=homoaloc
          if (lumoaloc.lt.lumoa) lumoa=lumoaloc
          if (lumoaloc-homoaloc.lt.directgapa) then  
             directgapa=lumoaloc-homoaloc
          end if
          if(directgapa.lt.dirfundgap) then
             dirfundgap=directgapa
             dirfundhomo=homoaloc
             dirfundlumo=lumoaloc
          end if
          if (spinpol) then
            homobloc=maxval(levels_down(ikpoint,1:nbands_printed),        &
     &            occups_down(ikpoint,1:nbands_printed).gt.occtol)
            lumobloc=minval(levels_down(ikpoint,1:nbands_printed),        &
     &            occups_down(ikpoint,1:nbands_printed).le.occtol)
            if (homobloc.gt.homob) homob=homobloc
            if (lumobloc.lt.lumob) lumob=lumobloc
            if (lumobloc-homobloc.lt.directgapb) then  
                directgapb=lumobloc-homobloc
            end if
            if(directgapb.lt.dirfundgap) then
               dirfundgap=directgapb
               dirfundhomo=homobloc
               dirfundlumo=lumobloc
            end if
          end if ! spinpol
          !
          ! end get local VBM and CBM for up and down spins
          !
        end do ! ikpoint
        !
        ! BEGIN DEBUG
        !
        !print '(8x,"direct gaps obtained")'
        !
        ! END DEBUG
        !
        gapa=lumoa-homoa
        fundgap=gapa
        if (spinpol) then
           gapb=lumob-homob
           fundgap=min(gapa,gapb)
        end if
        !
        ! BEGIN DEBUG
        !
        !print '(8x,"fundamental gaps obtained")'
        !
        ! END DEBUG
        !
        !
        ! end get HOMO, LUMO, gap, direct gap
        !
        ! begin set gap of metals to zero 
        !
        if (mod(nelect,2).gt.0.and..not.spinpol.and..not.spinorbit)       &
     &      directgapa=gapa
        if (mod(nelect,2).gt.0.and.spinpol.and.abs(magmom).lt.1D-1) then
          directgapa=gapa
          directgapb=gapb
          fundgap=min(gapa,gapb)
          dirfundgap=min(gapa,gapb)
          dirfundhomo=max(homoa,homob)
          dirfundlumo=min(lumoa,lumob)
        end if
        if (any(abs(occups_up(:,:)-occtol).lt.0.5d0*occtol)) then
          gapa=0.0d0
          homoa=efermi
          lumoa=efermi
          directgapa=gapa
          fundgap=0.0d0
          dirfundgap=min(gapa,gapb)
          dirfundhomo=efermi !max(homoa,homob)
          dirfundlumo=efermi !min(lumoa,lumob)
        end if
        if (any(abs(occups_down(:,:)-occtol).lt.0.5d0*occtol)) then
          gapb=0.0d0
          homob=efermi
          lumob=efermi
          directgapb=gapb
          fundgap=0.0d0
          dirfundgap=min(gapa,gapb)
          dirfundhomo=efermi !max(homoa,homob)
          dirfundlumo=efermi !min(lumoa,lumob)
        end if
        max_nele_k_up=0.0d0
        min_nele_k_up=nelect
        if(spinpol) then
          max_nele_k_down=0.0d0
          min_nele_k_down=nelect
        end if
        do ikpoint=1,nkpoints
          if (sum(occups_up(ikpoint,1:nbands_printed)).gt.max_nele_k_up)  &
     &        max_nele_k_up=sum(occups_up(ikpoint,1:nbands_printed))
          if (sum(occups_up(ikpoint,1:nbands_printed)).lt.min_nele_k_up)  &
     &        min_nele_k_up=sum(occups_up(ikpoint,1:nbands_printed))
          if (spinpol) then
            if (sum(occups_down(ikpoint,1:nbands_printed)).gt.            &
     &          max_nele_k_down) max_nele_k_down=sum(occups_down(         &
     &          ikpoint,1:nbands_printed))
            if (sum(occups_down(ikpoint,1:nbands_printed)).lt.            &
     &          min_nele_k_down) min_nele_k_down=sum(occups_down(         &
     &          ikpoint,1:nbands_printed))
          end if
        end do
        ! if the system is a metal, any finite band gap is an artefact
        ! of sparse k-point sampling. Correct for this:
        do ikpoint=1,nkpoints
         if (.not.spinpol.and.(sum(occups_up(ikpoint,1:nbands_printed))   &
     &       .lt.nelect-0.5d0*occtol.or.sum(occups_up(ikpoint,            &
     &       1:nbands_printed)).gt.nelect+0.5d0*occtol)) then
           gapa=0.0d0
           homoa=efermi
           lumoa=efermi
           directgapa=gapa
           fundgap=0.0d0
           dirfundgap=min(gapa,gapb)
           dirfundhomo=efermi !max(homoa,homob)
           dirfundlumo=efermi !min(lumoa,lumob)
         end if
         if (spinpol.and.max_nele_k_up-min_nele_k_up.gt.                  &
     &       0.5d0*occtol) then
           gapa=0.0d0
           homoa=efermi
           lumoa=efermi
           directgapa=gapa
           fundgap=0.0d0
           dirfundgap=0.0d0
           dirfundhomo=efermi !max(homoa,homob)
           dirfundlumo=efermi !min(lumoa,lumob)
         end if
         if (spinpol.and.min_nele_k_down.lt.max_nele_k_down-0.5d0         &
     &       *occtol) then
           gapb=0.0d0
           homob=efermi
           lumob=efermi
           directgapb=gapb
           fundgap=0.0d0
           dirfundgap=0.0d0
           dirfundhomo=efermi !max(homoa,homob)
           dirfundlumo=efermi !min(lumoa,lumob)
         end if
        end do ! ikpoint 
        !
        ! end set gap of metals to zero 
        !
        !
        ! BEGIN DEBUG
        !
        !print '(8x,"gaps corrected for metals")'
        !
        ! END DEBUG
        !
      print '(8x,"occtol: ",F10.6)',occtol
      print '(8x,"tol_en: ",F10.6, " eV")',tol_en
      print '(8x,"")'
      print '(8x,"alpha homo in Hartree:",F10.6)',homoa/Hartree
      print '(8x,"alpha homo in eV:",F10.6)',homoa
      print '(8x,"")'
      !
      ! begin find kpoint(s) with VBM 
      !
      if (allocated(levels_up)) then
        do ibandup=1,nbands
         do ikpoint=1,nkpoints
          if (abs(levels_up(ikpoint,ibandup)-homoa).lt.tol_en             &
     &        .and.occups_up(ikpoint,ibandup).gt.0.5*occtol) then
            print '(8x,"alpha homo with eigenvalue ",F12.6,               &
     &            " at kpoint # ",I5," in band # ",I5,                    &
     &            " with fractional coordinates ", 3(F12.7))',            &
     &           levels_up(ikpoint,ibandup),ikpoint,ibandup,              &
     &           kpoints(ikpoint,1:3)
          end if           
         end do ! ikpoint
        end do ! ibandup
      end if ! allocated(levelsup)
      !
      ! end find kpoint(s) with VBM 
      !
      print '(8x,"")'
      print '(8x,"alpha lumo in Hartree:",F10.6)',lumoa/Hartree
      print '(8x,"alpha lumo in eV:",F10.6)',lumoa
      print '(8x,"")'
      !
      ! begin find kpoint(s) with CBM 
      !
      if (allocated(levels_up)) then
        do ibandup=1,nbands
         do ikpoint=1,nkpoints
          if (abs(levels_up(ikpoint,ibandup)-lumoa).lt.tol_en             &
     &        .and.occups_up(ikpoint,ibandup).lt.0.5*occtol) then
            print '(8x,"alpha lumo with eigenvalue ",F12.6,               &
     &            " at kpoint # ",I5," in band # ",I5,                    &
     &            " with fractional coordinates ", 3(F12.7))',            &
     &           levels_up(ikpoint,ibandup),ikpoint,ibandup,              &
     &           kpoints(ikpoint,1:3)
          end if           
         end do ! ikpoint
        end do ! ibandup
      end if ! allocated(levelsup)
      print '(8x,"")'
      !
      ! end find kpoint(s) with CBM 
      !
      print '(8x,"alpha bandgap in eV:",F10.6)',gapa
      print '(8x,"alpha direct bandgap in eV:",F10.6)',directgapa
      print '(8x,"")'
      !
      ! begin find kpoint(s) with direct band gap 
      !
      if (allocated(levels_up)) then
        do ibandup=1,nbands
         do ikpoint=1,nkpoints
         if (ibandup.lt.nbands) then
          if (abs(levels_up(ikpoint,ibandup+1)-levels_up(ikpoint,         &
     &     ibandup)-directgapa).lt.2.0d0*tol_en.and.directgapa.gt.tol_en  &
     &        .and.occups_up(ikpoint,ibandup+1).lt.0.5*occtol.and.        &
     &        occups_up(ikpoint,ibandup).gt.0.5*occtol) then
            print '(8x,"alpha direct gap of ",F12.6," eV at kpoint # ",   &
     &      I5," between bands # ",I5," and ",I5,                         &
     &      " with fractional coordinates ",3(F12.7))',                   &
     &      levels_up(ikpoint,ibandup+1)-levels_up(ikpoint,ibandup),      &
     &      ikpoint,ibandup,ibandup+1,kpoints(ikpoint,1:3)
          end if
         end if          
         end do ! ikpoint
        end do ! ibandup
      end if ! allocated(levelsup)
      print '(8x,"")'
      !
      ! end find kpoint(s) with direct band gap 
      !
      open(52,file='BANDEMINMAXUP.DAT',status='replace')
      do ibandup=1,nbands
        write(52,'(I10,5x,2(F20.8,x))') ibandup,bandeminup(ibandup),      &
     &       bandemaxup(ibandup)
      end do ! ibandup
      close(52)
      !
      if (allocated(levels_up)) then
        open(52,file='LEVELSUP.DAT',status='replace')
        write(52,'("# kpoint, band, eigenvalue, weight, occupation")')
        do ibandup=1,nbands
          do ikpoint=1,nkpoints
            write(52,'(2(I10,2x),5x,3(F20.8,x))') ikpoint,ibandup,        &
     &      levels_up(ikpoint,ibandup), kweights(ikpoint),                &
     &         occups_up(ikpoint,ibandup)
          end do ! ikpoint
        end do ! ibandup
        close(52)
      end if ! (allocated(levels_up))
      !  
      if(spinpol) then
        !
        print '(8x,"beta homo in Hartree:",F10.6)',homob/Hartree
        print '(8x,"beta homo in eV:",F10.6)',homob
        print '(8x,"")'
        !
        ! begin find kpoint(s) with VBM 
        !
        if (allocated(levels_down)) then
          do ibanddown=1,nbands
           do ikpoint=1,nkpoints
            if (abs(levels_down(ikpoint,ibanddown)-homob).lt.tol_en       &
     &        .and.occups_down(ikpoint,ibanddown).gt.0.5*occtol) then
              print '(8x,"beta homo with eigenvalue ",F12.6,              &
     &        " at kpoint # ",I5," in band # ",I5,                        &
     &        " with fractional coordinates ",3(F12.7))',                 &
     &        levels_down(ikpoint,ibanddown),ikpoint,ibanddown,           &
     &        kpoints(ikpoint,1:3)
            end if           
           end do ! ibanddown
          end do ! ikpoint
        end if ! (allocated(levels_down))
        print '(8x,"")'
        !
        ! end find kpoint(s) with VBM 
        !
        print '(8x,"beta lumo in Hartree:",F10.6)',lumob/Hartree
        print '(8x,"beta lumo in eV:",F10.6)',lumob
        print '(8x,"")'
        !
        ! begin find kpoint(s) with CBM 
        !
        if (allocated(levels_down)) then
          do ibanddown=1,nbands
           do ikpoint=1,nkpoints
            if (abs(levels_down(ikpoint,ibanddown)-lumob).lt.tol_en       &
     &        .and.occups_down(ikpoint,ibanddown).lt.0.5*occtol) then
              print '(8x,"beta lumo with eigenvalue ",F12.6,              &
     &              " at kpoint # ",I5," in band # ",I5,                  &
     &         " with fractional coordinates ",3(F12.7))',                &
     &       levels_down(ikpoint,ibanddown),ikpoint,ibanddown,            &
     &       kpoints(ikpoint,1:3)
            end if           
           end do ! ibanddown
          end do ! ikpoint
        end if ! (allocated(levels_down))
        !
        ! end find kpoint(s) with CBM 
        !
        print '(8x,"")'
        print '(8x,"beta bandgap in eV:",F10.6)',gapb
        print '(8x,"beta direct bandgap in eV:",F10.6)',directgapb
        print '(8x,"")'
        !
        ! begin find kpoint(s) with direct band gap 
        !
        if (allocated(levels_down)) then
          do ibanddown=1,nbands
           do ikpoint=1,nkpoints
           if (ibanddown.lt.nbands) then
            if(abs(levels_down(ikpoint,ibanddown+1)-levels_down(ikpoint,  &
     &       ibanddown)-directgapb).lt.2.0d0*tol_en.and.directgapb.gt.    &
     &       tol_en.and.occups_down(ikpoint,ibanddown+1).lt.0.5*occtol    &
     &       .and. occups_down(ikpoint,ibanddown).gt.0.5*occtol) then
              print '(8x,"beta direct gap of ",F12.6," eV at kpoint # ",  &
     &        I5," between bandss # ",I5," and ",I5,                      &
     &        " with fractional coordinates ",3(F12.7))',                 &
     &        levels_down(ikpoint,ibanddown+1)-levels_down(ikpoint,       &
     &        ibanddown),ikpoint,ibanddown,ibanddown+1,                   &
     &        kpoints(ikpoint,1:3)
            end if
           end if          
           end do ! ibanddown
          end do ! ikpoint
        end if ! (allocated(levels_down))
        print '(8x,"")'
        !
        ! end find kpoint(s) with direct band gap 
        !
        open(52,file='BANDEMINMAXDOWN.DAT',status='replace')
        do ibanddown=1,nbands
          write(52,'(I10,5x,2(F20.8,x))') ibanddown,                      &
     &      bandemindown(ibanddown),bandemaxdown(ibanddown)
        end do ! ibanddown
        close(52)
        !
        if (allocated(levels_down)) then
          open(52,file='LEVELSDOWN.DAT',status='replace')
          do ibanddown=1,nbands
            do ikpoint=1,nkpoints
              write(52,'(2(I10,2x),5x,3(F20.8,x))') ikpoint,ibanddown,    &
     &        levels_down(ikpoint,ibanddown), kweights(ikpoint),          &
     &           occups_down(ikpoint,ibanddown)
            end do ! ikpoint
          end do ! ibanddown
          close(52)
        end if ! (allocated(levels_down))
        !  
      end if ! (spinpol)
      !
      print '(8x,"fundamental bandgap in eV:",F10.6)',fundgap
      print '(8x,"fundamental direct bandgap in eV:",F10.6)',dirfundgap
      print '(8x,"direct fundamental HOMO in eV:",F10.6)',dirfundhomo
      print '(8x,"direct fundamental LUMO in eV:",F10.6)',dirfundlumo
      !
      ! begin write VASP EIGENVAL file
      !
      if (allocated(kpoints)) then
        open(52,file='COFIMA.EIGENVAL',status='replace')
        if (.not.spinpol) write(52,*) "header header header 1"
        if (spinpol) write(52,*) "header header header 2"
        write(52,*) "header"
        write(52,*) "header"
        write(52,*) "header"
        write(52,*) "unknown system"
        write(52,*) nelect, nkpoints,nbands 
        do ikpoint=1,nkpoints
          write(52,*) " "
          write(52,'(3(F16.12,1x))') kpoints(ikpoint,1:3)
          do ibandup=1,nbands
            if (spinpol) then
              write(52,'(1(I10,2x),5x,2(F20.8,x))') ibandup,              &
     &      levels_up(ikpoint,ibandup), levels_down(ikpoint,ibandup)    
            else
              write(52,'(1(I10,2x),5x,1(F20.8,x))') ibandup,              &
     &      levels_up(ikpoint,ibandup)    
            end if
          end do ! ibandup
        end do ! ikpoint
        close(52)
        !
        ! end write VASP EIGENVAL file
        !
      end if ! allocated kpoints
      !
      end subroutine print_gap_info
      !
c---------------------------------------------------------------------
      !
      end subroutine bandgap

c---------------------------------------------------------------------

      subroutine readxcme(infile,informat)

      use defs
      implicit none

      character (len=*) , intent(in) :: infile,informat
      logical hasxcme

      ! internal variables
      character line*256
      integer nbasis,i ,j ,iline,ispin
      double precision, allocatable :: xcme(:,:,:),occ(:,:),fock(:,:,:)
      double precision, allocatable :: evals(:,:),evecs(:,:,:),
     &     gvxc(:,:,:),hartreem(:,:)
      double precision esum

      print fsubstart,"readxcme"

      select case(informat)
      case("nwchem","NWCHEM")
        open(51,file=infile,status='old',err=100)  
        ! get number of basis functions (matrix dim=nbasis**2) 
8       read(51,'(A256)',err=101,end=102) line  
        !print*, line
        if(index(line,"AO basis - number of functions:").gt.0)then    
          read(line(42:47),*,err=102) nbasis
          print '(8x,"Number of basis functions:",I6)',nbasis
          allocate(xcme(2,nbasis,nbasis),occ(2,nbasis),
     &        fock(2,nbasis,nbasis),evals(2,nbasis),
     &        evecs(2,nbasis,nbasis),gvxc(2,nbasis,nbasis),
     &        hartreem(nbasis,nbasis))
          xcme=0.0D0
          fock=0.0D0
          gvxc=0.0D0
          hartreem=0.0D0
          evals=0.0D0
          evecs=0.0D0
          occ=0.0D0
          goto 10
        end if
        goto 8 
        ! read xc matrix elements
10      read(51,'(A256)',err=101,end=15) line
        if(index(line,"tmp matrix[").gt.0.or.
     &    index(line,"Vxcs[").gt.0) then
          if(index(line,"tmp matrix[").gt.0) ispin=1
          if(index(line,"Vxcs[").gt.0) ispin=2
          j=1
          do while (j.le.nbasis-5-mod(nbasis,6))
            do i=1,3
              read(51,'(A256)',err=101) line
            end do        
            do i=1,nbasis
              read(51,'(A256)',err=101) line
              read(line(8:78),*,err=103) xcme(ispin,i,j:j+5)
            end do
            j=j+6
          end do
          do i=1,3
            read(51,'(A256)',err=101) line
          end do        
          do i=1,nbasis
            read(51,'(A256)',err=101) line
            read(line(8:78),*,err=103) xcme(ispin,i,j:j+mod(nbasis,6)-1)
          end do
          j=j+mod(nbasis,6)
        end if    
        ! read Fock matrix elements
        if (index(line,"fock transf[").gt.0) then
          ispin=1
          j=1
          do while (j.le.nbasis-5-mod(nbasis,6))
            do i=1,3
              read(51,'(A256)',err=101) line
            end do        
            do i=1,nbasis
              read(51,'(A256)',err=101) line
              read(line(8:78),*,err=103) fock(ispin,i,j:j+5)
            end do
            j=j+6
          end do
          do i=1,3
            read(51,'(A256)',err=101) line
          end do        
          do i=1,nbasis
            read(51,'(A256)',err=101) line
            read(line(8:78),*,err=103) fock(ispin,i,j:j+mod(nbasis,6)-1)
          end do
          j=j+mod(nbasis,6)
        end if    
        ! read eigenvectors
        if(index(line,"alpha evecs[").gt.0.or.
     &    index(line,"beta evecs[").gt.0) then
          if(index(line,"alpha evecs[").gt.0) ispin=1
          if(index(line,"beta evecs[").gt.0) ispin=2
          j=1
          do while (j.le.nbasis-5-mod(nbasis,6))
            do i=1,3
              read(51,'(A256)',err=101) line
            end do        
            do i=1,nbasis
              read(51,'(A256)',err=101) line
              read(line(8:78),*,err=103) evecs(ispin,i,j:j+5)
            end do
            j=j+6
          end do
          do i=1,3
            read(51,'(A256)',err=101) line
          end do        
          do i=1,nbasis
            read(51,'(A256)',err=101) line
            read(line(8:78),*,err=103)evecs(ispin,i,j:j+mod(nbasis,6)-1)
          end do
          j=j+mod(nbasis,6)
        end if    
        ! read gvxc matrix elements
        if (index(line,"g vxc 1[").gt.0) then
          ispin=1
          j=1
          do while (j.le.nbasis-5-mod(nbasis,6))
            do i=1,3
              read(51,'(A256)',err=101) line
            end do        
            do i=1,nbasis
              read(51,'(A256)',err=101) line
              read(line(8:78),*,err=103) gvxc(ispin,i,j:j+5)
            end do
            j=j+6
          end do
          do i=1,3
            read(51,'(A256)',err=101) line
          end do        
          do i=1,nbasis
            read(51,'(A256)',err=101) line
            read(line(8:78),*,err=103) gvxc(ispin,i,j:j+mod(nbasis,6)-1)
          end do
          j=j+mod(nbasis,6)
        end if    
        if (index(line,"g vxc 2[").gt.0) then
          ispin=2
          j=1
          do while (j.le.nbasis-5-mod(nbasis,6))
            do i=1,3
              read(51,'(A256)',err=101) line
            end do        
            do i=1,nbasis
              read(51,'(A256)',err=101) line
              read(line(8:78),*,err=103) gvxc(ispin,i,j:j+5)
            end do
            j=j+6
          end do
          do i=1,3
            read(51,'(A256)',err=101) line
          end do        
          do i=1,nbasis
            read(51,'(A256)',err=101) line
            read(line(8:78),*,err=103) gvxc(ispin,i,j:j+mod(nbasis,6)-1)
          end do
          j=j+mod(nbasis,6)
        end if    
        ! read Hartree (?) matrix elements
        if (index(line,"global array: jk[1:7,1:7]").gt.0) then
          ispin=1
          j=1
          do while (j.le.nbasis-5-mod(nbasis,6))
            do i=1,3
              read(51,'(A256)',err=101) line
            end do        
            do i=1,nbasis
              read(51,'(A256)',err=101) line
              read(line(8:78),*,err=103) hartreem(i,j:j+5)
            end do
            j=j+6
          end do
          do i=1,3
            read(51,'(A256)',err=101) line
          end do        
          do i=1,nbasis
            read(51,'(A256)',err=101) line
            read(line(8:78),*,err=103) hartreem(i,j:j
     &       +mod(nbasis,6)-1)
          end do
          j=j+mod(nbasis,6)
        end if    
        ! read eigenvalues
        if (index(line,"alpha eigenvalues").gt.0.or. 
     &      index(line,"beta eigenvalues").gt.0) then
          if (index(line,"alpha eigenvalues").gt.0) ispin=1
          if (index(line,"beta eigenvalues").gt.0) ispin=2
          do iline=1,6
            read(51,'(A256)',err=101) line
          end do
          do while (index(line,'alpha foccs').le.0)
            read(line(1:256),*) i,evals(ispin,i)
            read(51,'(A256)',err=101) line
          end do
        end if    
        ! read occupations
        ! First spin channel. Since not always all vectors are printed,
        ! set occupations of the lower unprinted vectors to 1 
        ! (spin-polarized case) or 2 (unpolarized case)
        if(index(line,"DFT Final Molecular Orbital Analysis").gt.0
     &     .or.index(line,"DFT Final Alpha Molecular Orbital Analysis")
     &     .gt.0.
     &     .or.index(line,"DFT Final Beta Molecular Orbital Analysis")
     &     .gt.0) then
          if(index(line,"DFT Final Molecular Orbital Analysis").gt.0
     &      .or.index(line,"DFT Final Alpha Molecular Orbital Analysis")
     &      .gt.0 ) ispin=1
          if(index(line,"DFT Final Beta Molecular Orbital Analysis")
     &      .gt.0) ispin=2
          occ(ispin,1:nbasis)=0.0D0
          do iline=1,3
            read(51,'(A256)',err=101,end=15) line
          end do 
          read(line(8:14),*) i 
          read(line(19:30),*) occ(ispin,i)
          do j=1,i
            occ(ispin,j)=occ(ispin,i)
          end do
        end if 
        ! read the remaining printed occupations
        if(index(line,"Vector").gt.0.and.index(line,'Occ=').gt.0)then
          read(line(8:14),*) i 
          read(line(19:30),*) occ(ispin,i)
        end if
        goto 10
15      continue
        close(51)
      case default
        goto 900
      end select

      close(51)
      open(51,file="XCME.DAT",status="replace")
      ! print all XC matrix elements
      write(51,'("# ispin, i, j, <i|Vxc|j>")')
      do ispin=1,2 
        j=1
        do while (j.le.nbasis-5-mod(nbasis,6))
          do i=1,nbasis
            write(51,'(I5,I5,I5,"...",I5,6(F20.10))') 
     &          ispin,i,j,j+5,xcme(ispin,i,j:j+5)
          end do
          j=j+6
        end do
        do i=1,nbasis
          write(51,*) ispin,i,j,"...",j+mod(nbasis,6)-1,
     &    xcme(ispin,i,j:j+mod(nbasis,6)-1)
        end do
      end do
      ! print diagonal elements
      write(51,'("# ispin, i, occ(i), <i|Vxc|i>")')
      do ispin=1,2
        do i=1,nbasis
          write(51,'(2(I5),2(F20.10))') ispin,i,occ(ispin,i),
     &     xcme(ispin,i,i)
        end do
      end do 
      ! print Vxc (sum of diagonal ME's of occupied states)
      esum=0.0D0
      do i=1,nbasis
        do ispin=1,2
          esum=esum+xcme(ispin,i,i)*occ(ispin,i)
        end do
      end do
      write(51,'("# Vxc=",F20.10)') esum
      close(51)
      ! print GVXC matrix elements
      open(51,file="GVXC.DAT",status="replace")
      ! print all matrix elements
      write(51,'("# ispin, i, j, <i|Vxc|j>")')
      do ispin=1,2 
        j=1
        do while (j.le.nbasis-5-mod(nbasis,6))
          do i=1,nbasis
            write(51,'(I5,I5,I5,"...",I5,6(F20.10))') 
     &          ispin,i,j,j+5,gvxc(ispin,i,j:j+5)
          end do
          j=j+6
        end do
        do i=1,nbasis
          write(51,*) ispin,i,j,"...",j+mod(nbasis,6)-1,
     &    gvxc(ispin,i,j:j+mod(nbasis,6)-1)
        end do
      end do
      ! print diagonal elements
      write(51,'("# ispin, i, occ(i), <i|Vxc|i>")')
      do ispin=1,2
        do i=1,nbasis
          write(51,'(2(I5),2(F20.10))') ispin,i,occ(ispin,i),
     &     gvxc(ispin,i,i)
        end do
      end do 
      ! print Vxc (sum of diagonal ME's of occupied states)
      esum=0.0D0
      do i=1,nbasis
        do ispin=1,2
          esum=esum+gvxc(ispin,i,i)*occ(ispin,i)
        end do
      end do
      write(51,'("# Vxc=",F20.10)') esum
      close(51)
      ! print Hartree (?) matrix elements
      open(51,file="HARTREE.DAT",status="replace")
      ! print all matrix elements
      write(51,'("#  i, j, <i|V_Hartree|j>")')
      j=1
      do while (j.le.nbasis-5-mod(nbasis,6))
        do i=1,nbasis
          write(51,'(I5,I5,"...",I5,6(F20.10))') 
     &        i,j,j+5,Hartreem(i,j:j+5)
        end do
        j=j+6
      end do
      do i=1,nbasis
        write(51,*) i,j,"...",j+mod(nbasis,6)-1,
     &  Hartreem(i,j:j+mod(nbasis,6)-1)
      end do
      ! print diagonal elements
      write(51,'("# i, occ(i), <i|V_Hartree|i>")')
      do i=1,nbasis
        write(51,'(1(I5),2(F20.10))') i,occ(1,i)+occ(2,i),
     &   Hartreem(i,i)
      end do
      ! print VHartree (sum of diagonal ME's of occupied states)
      esum=0.0D0
      do i=1,nbasis
        do ispin=1,2
          esum=esum+Hartreem(i,i)*occ(ispin,i)
        end do
      end do
      write(51,'("# VHartree=",F20.10)') esum
      close(51)
      ! print Fock matrix
      open(51,file="FOCK.DAT",status="replace")
      ! print all matrix elements
      write(51,'("# ispin, i, j, <i|Fock|j>")')
      do ispin=1,2 
        j=1
        do while (j.le.nbasis-5-mod(nbasis,6))
          do i=1,nbasis
            write(51,'(I5,I5,I5,"...",I5,6(F20.10))') 
     &          ispin,i,j,j+5,fock(ispin,i,j:j+5)
          end do
          j=j+6
        end do
        do i=1,nbasis
          write(51,*) ispin,i,j,"...",j+mod(nbasis,6)-1,
     &    fock(ispin,i,j:j+mod(nbasis,6)-1)
        end do
      end do
      ! print diagonal elements
      write(51,'("# ispin, i, occ(i), <i|Fock|i>")')
      do ispin=1,2
        do i=1,nbasis
          write(51,'(2(I5),2(F20.10))') ispin,i,occ(ispin,i),
     &     fock(ispin,i,i)
        end do
      end do 
      ! print sum of Fock (sum of diagonal ME's of occupied states)
      esum=0.0D0
      do i=1,nbasis
        do ispin=1,2
          esum=esum+fock(ispin,i,i)*occ(ispin,i)
        end do
      end do
      write(51,'("# Focksum=",F20.10)') esum
      close(51)
      ! print eigenvector matrix elements
      open(51,file="EVECS.DAT",status="replace")
      ! print all matrix elements
      write(51,'("# ispin, i, j, <KS,j|GTO,i>")')
      do ispin=1,2 
        j=1
        do while (j.le.nbasis-5-mod(nbasis,6))
          do i=1,nbasis
            write(51,'(I5,I5,I5,"...",I5,6(F20.10))') 
     &          ispin,i,j,j+5,evecs(ispin,i,j:j+5)
          end do
          j=j+6
        end do
        do i=1,nbasis
          write(51,*) ispin,i,j,"...",j+mod(nbasis,6)-1,
     &    evecs(ispin,i,j:j+mod(nbasis,6)-1)
        end do
      end do
      ! print diagonal elements
      write(51,'("# ispin, i, occ(i), sum_j (M_ji*M_ji)")')
      do ispin=1,2
        do i=1,nbasis
          write(51,'(2(I5),2(F20.10))') ispin,i,occ(ispin,i),
     &     sum(evecs(ispin,1:nbasis,i)**2)
        end do
      end do 
      close(51)
      ! print eigenvalues
      open(51,file="EVALS.DAT",status="replace")
      ! print all eigenvalues
      write(51,'("# ispin, i, epsilon_i")')
      do ispin=1,2 
        do i=1,nbasis
          write(51,'(I5,I5,1(F20.10))') ispin,i,evals(ispin,i)
        end do
      end do
      ! print sum of eigenvalues of occupied states (band energy)
      esum=0.0D0
      do i=1,nbasis
        do ispin=1,2
          esum=esum+evals(ispin,i)*occ(ispin,i)
        end do
      end do
      write(51,'("# sum of eigenvalues of occupied states (band energy)=
     &",F20.10)') esum
      close(51)
      ! finish 
      print fsubend
      return

! Errors
100   nerr=nerr+1
      close(51)
      print ferrmssg,"Could not find/open inpt file"
      return
101   nerr=nerr+1
      close(51)
      print ferrmssg,"Could not read input file"
      return
102   nerr=nerr+1
      close(51)
      print ferrmssg,"Could not read number of functions"
      return
103   nerr=nerr+1
      close(51)
      print ferrmssg,"Could not read matrixelements"
      return
900   nerr=nerr+1
      close(51)
      print ferrmssg,"File seems not to contain xc me's"
      return
      end subroutine

c---------------------------------------------------------------------

      function binomial(n,k)
      ! calculates (n over k)=n!/(k!*(n-k)!) 
      implicit none
      !integer(kind=16) n,k,binomial
      !integer(kind=16) i,j,x,y
      integer(8) n,k,binomial
      integer(8) i,j,x,y
      
      ! n!/k!=(k+1)*...*n
      x=1
      do i=k+1,n
        x=x*i
      end do
      
      ! (n-k)!
      y=1
      do i=2,n-k
        y=y*i
      end do
    
      !print*,"n!/k!=",x,"(n-k)!=",y
      
      binomial=x/y  
  
      end function

c---------------------------------------------------------------------

      subroutine coord(atoms,vecs,nndist,coordmin,coordmax,coordmean)
      ! finds the coordination numbers of the atoms
      use defs      
      use combi
      implicit none

      type(atom) atoms(:)
      double precision coordmean,vecs(:,:),nndist  ! average coordination number, cell vectors, 
          ! nearest neighbor distance
      integer coordmin,coordmax ! minimum, maximum coordination number
      ! internal variables:
      double precision site1(1:3),site2(1:3),distvec(1:3)
      double precision site3(1:3),distvec2(1:3),angle
      integer natoms,iatom1,iatom2,i1,i2,i3
      integer coordnum,icoord,icoord2
      integer nbonds,nangles,iangle

      if(talk) print fsubstart,"coord"
      coordmin=10000 
      coordmax=0 
      coordmean=0.0D0 
      if(talk) print '(8x,"using cutoff of ",F10.6," Angs")',nndist
     
      ! get coordination numbers 
      natoms=size(atoms)
      do iatom1=1,natoms
        coordnum=0  ! set initial coordination number to 0
        call frac2abs(atoms(iatom1)%where,vecs,site1)
        do iatom2=1,natoms
          do i1=-1,1
            do i2=-1,1
              do i3=-1,1
                call frac2abs(atoms(iatom2)%where,vecs,site2)
                site2=site2+dble(i1)*vecs(1,1:3)
     &             +dble(i2)*vecs(2,1:3)+dble(i3)*vecs(3,1:3)
                distvec=site2-site1
                ! if the two atoms are not the same and their distance
                ! is below the threshold, add 1 to the coordination
                ! number of the first atom
                if(absvec(distvec).le.nndist.and.iatom1.ne.iatom2) 
     &            coordnum=coordnum+1
              end do  ! i3
            end do  ! i2
          end do  ! i1
        end do ! iatom2
        atoms(iatom1)%coord=coordnum
        if (coordnum.gt.coordmax) coordmax=coordnum
        if (coordnum.lt.coordmin) coordmin=coordnum
        coordmean=coordmean+dble(coordnum)
        !print '(8x,"iatom,coord=",I3,1x,I20)',iatom1,coordnum
      end do ! iatom1
      coordmean=coordmean/dble(natoms) 
      !print '(8x,"coordination number checked")'
      
      ! get names of neighbors 
      do iatom1=1,natoms
        if(allocated(atoms(iatom1)%neighbors)) 
     &     deallocate(atoms(iatom1)%neighbors) 
        allocate(atoms(iatom1)%neighbors(1:atoms(iatom1)%coord))
        if (allocated(atoms(iatom1)%neighborcoords))                      &
     &    deallocate(atoms(iatom1)%neighborcoords)
        allocate(atoms(iatom1)%neighborcoords(1:atoms(iatom1)%coord,      &
     &    1:3))
          if (allocated(atoms(iatom1)%ineighbors))                        &
     &     deallocate(atoms(iatom1)%ineighbors)
        allocate(atoms(iatom1)%ineighbors(1:atoms(iatom1)%coord))
        icoord=0
        call frac2abs(atoms(iatom1)%where,vecs,site1)
        do iatom2=1,natoms
          do i1=-1,1
            do i2=-1,1
              do i3=-1,1
                call frac2abs(atoms(iatom2)%where,vecs,site2)
                site2=site2+dble(i1)*vecs(1,1:3)
     &             +dble(i2)*vecs(2,1:3)+dble(i3)*vecs(3,1:3)
                distvec=site2-site1
                ! if the two atoms are not the same and their distance
                ! is below the threshold, they are neighbors
                if(absvec(distvec).le.nndist.and.iatom1.ne.iatom2) then
                   icoord=icoord+1
                  atoms(iatom1)%neighbors(icoord)
     &              =atoms(iatom2)%name(1:2)
                  atoms(iatom1)%neighborcoords(icoord,1:3)=site2
                  atoms(iatom1)%ineighbors(icoord)=iatom2
                end if
              end do  ! i3
            end do  ! i2
          end do  ! i1
        end do ! iatom2
      end do ! iatom1
      !print '(8x,"neighbor names checked")'

      !
      ! begin get bond lengths and angles
      !
      do iatom1=1,natoms
        !print*,'atom=',atoms(iatom1)%name
        nbonds=atoms(iatom1)%coord 
        !print '(8x,"nbonds=",I10)',nbonds
        !print*,'nbonds=',nbonds
        if (nbonds.gt.0) then
          if(allocated(atoms(iatom1)%bondlengths)) 
     &       deallocate(atoms(iatom1)%bondlengths)
          allocate(atoms(iatom1)%bondlengths(nbonds))
          !print*,facult(14)
          !print*,facult(20)
          !print*,n_over_k_large(14,2)
          !nangles=n_over_k_large(nbonds,2)
          nangles=0
          if (nbonds>0) nangles=nbonds*(nbonds-1)/2 
          if(allocated(atoms(iatom1)%bondangles)) 
     &       deallocate(atoms(iatom1)%bondangles)
          !print '(8x,"nangles checked for atom ",I20)',iatom1
          !print '(8x,"nangles=",I20)',nangles
          if(nangles.gt.0) allocate(atoms(iatom1)%bondangles(1:nangles))
          call frac2abs(atoms(iatom1)%where,vecs,site1)
          iangle=1
          ! bond length
          do icoord=1,nbonds
            site2=atoms(iatom1)%neighborcoords(icoord,1:3)
            distvec=site2-site1
            !print '(8x,"distvec=",3(1x,F10.6))',distvec
            atoms(iatom1)%bondlengths(icoord)=norm2(distvec)
            ! bond angle
            if (nangles.gt.0) then
              do icoord2=icoord+1,nbonds
                site3=atoms(iatom1)%neighborcoords(icoord2,1:3)
                distvec2=site3-site1
                !print '(8x,"distvec2=",3(1x,F10.6))',distvec2
                angle=dot_product(distvec,distvec2)                       &
     &               /(norm2(distvec)*norm2(distvec2))
                !print*, "cos(angle)=",angle
                ! correct for small numerical error that may result in
                ! cos(angle)>1:
                if (angle>1.0d0 .and. angle<1.0d0+1.0D-6) angle=1.0d0
                if (angle<-1.0d0 .and. angle>-1.0d0-1.0D-6) angle=-1.0d0
                angle=acos(angle)*180.0d0/Pi
                !print*, "iangle,nangles,angle=",iangle,nangles,angle
                atoms(iatom1)%bondangles(iangle)=angle
                iangle=iangle+1
              end do ! icoord2
            end if ! nangles > 0
          end do ! icoord
        end if ! nbonds > 0
        !print '(8x,"angles checked for atom",I20)',iatom1
      end do ! iatom1
      !
      ! end get bond lengths and angles
      !

      if(talk) print fsubendext,"coord"
      return
   
      end subroutine coord

c---------------------------------------------------------------------

      subroutine coordfinite(atoms,nndist,coordmin,coordmax,
     &    coordmean) 
      ! the same as coord, but for finite systems (clusters) 
      ! caution: needs absolute atom positions as input
      use defs      
      implicit none

      type(atom) atoms(:)
      double precision coordmean,nndist  ! average coordination number, cell vectors, 
          ! nearest neighbor distance
      integer coordmin,coordmax ! minimum, maximum coordination number
      ! internal variables:
      double precision site1(1:3),site2(1:3),distvec(1:3)
      integer natoms,iatom1,iatom2,i1,i2,i3
      integer coordnum,icoord

      if(talk) print fsubstart,"coordfinite"
      coordmin=10000 
      coordmax=0 
      coordmean=0.0D0 
     
      ! get coordination numbers 
      natoms=size(atoms)
      do iatom1=1,natoms
        coordnum=0  ! set initial coordination number to 0
        site1=atoms(iatom1)%abswhere
        do iatom2=1,natoms
           site2=atoms(iatom2)%abswhere 
           distvec=site2-site1
           ! if the two atoms are not the same and their distance
           ! is below the threshold, add 1 to the coordination
           ! number of the first atom
           if(absvec(distvec).le.nndist.and.iatom1.ne.iatom2) 
     &       coordnum=coordnum+1
        end do ! iatom2
        atoms(iatom1)%coord=coordnum
        if (coordnum.gt.coordmax) coordmax=coordnum
        if (coordnum.lt.coordmin) coordmin=coordnum
        coordmean=coordmean+dble(coordnum)
      end do ! iatom1
      coordmean=coordmean/dble(natoms) 
      
      ! get names of neighbors 
      do iatom1=1,natoms
        if(allocated(atoms(iatom1)%neighbors))
     &        deallocate(atoms(iatom1)%neighbors)
        if(allocated(atoms(iatom1)%nnsites))
     &        deallocate(atoms(iatom1)%nnsites)
        if(atoms(iatom1)%coord.gt.0) then
          allocate(atoms(iatom1)%neighbors(1:atoms(iatom1)%coord))
          allocate(atoms(iatom1)%nnsites(1:atoms(iatom1)%coord,1:3))
        end if
        icoord=0
        site1= atoms(iatom1)%abswhere
        do iatom2=1,natoms
           site2=atoms(iatom2)%abswhere
           distvec=site2-site1
           ! if the two atoms are not the same and their distance
           ! is below the threshold, they are neighbors
           if(absvec(distvec).le.nndist.and.iatom1.ne.iatom2)then 
            icoord=icoord+1
            atoms(iatom1)%neighbors(icoord)
     &        =atoms(iatom2)%name(1:2)
            atoms(iatom1)%nnsites(icoord,1:3)
     &        =atoms(iatom2)%abswhere(1:3)
           end if
        end do ! iatom2
      end do ! iatom1

      if(talk) print fsubendext,"coordfinite"
      return
   
      end subroutine 

c---------------------------------------------------------------------

      subroutine dipmom(atoms,dipole) 
      ! simple routine to calculate a dipole moment of a cluster 
      ! caution, this routine expects that the absolute coordinates
      ! are passed. it does not calculate them.       

      use defs
      implicit none

      type(atom) atoms(:)
      double precision dipole(1:3)
      ! internal variables
      integer natoms  ! number of atoms
      integer i  ! atom index

      if(talk) print fsubstart,"dipole" 
      
      natoms=size(atoms)
      dipole=0.0D0
      !the following simply adds over coordinates times charges.
      do i=1,natoms
        dipole=dipole+atoms(i)%charge*atoms(i)%abswhere
      end do 
      if(talk) then 
        print*,'        dipole moment (e*Angs):',dipole
      end if
      
      end subroutine

c---------------------------------------------------------------------

      subroutine ecoulomb(atoms,ecoul) 
      ! simple routine to calculate the electrostaticenergy of a cluster 
      ! caution, this routine expects that the absolute coordinates
      ! and charges are passed. it does not determine them.       

      use defs
      implicit none

      type(atom) atoms(:)
      double precision ecoul
      ! internal variables
      integer natoms  ! number of atoms
      integer i,j  ! atom index
      double precision distvec(1:3)

      if(talk) print fsubstart,"ecoulomb" 
      
      natoms=size(atoms)
      ecoul=0.0D0
      !the following simply adds over coordinates times charges.
      do i=1,natoms
        do j=i+1,natoms
          distvec=atoms(i)%abswhere-atoms(j)%abswhere
          !if(absvec(distvec).gt.0.0D0) then
            ecoul=ecoul+atoms(i)%charge*atoms(j)%charge
     &        /absvec(distvec)
          !end if 
        end do 
      end do 
      ecoul=ecoul*ec/(4.0D0*Pi*epsilon0*Angstrom)
      if(talk) then 
        print*,'        E_Coulomb (eV):',ecoul
      end if
      
      end subroutine

c---------------------------------------------------------------------

      subroutine compare(file1,file2)
      !
      use defs
      use integr
      implicit none
      ! calculates the overlap Ov between two fields on a grid, which are
      ! read from file1 and file2: Ov= < f1 | f2 > / ( |f1| * |f2| )
      !
      character(len=*) file1,file2
      ! internal variables
      double precision, allocatable :: rgrid1(:),rgrid2(:),field1(:),
     &   field2(:) 
      character(len=256) line
      integer ndata1,ndata2, ndata ! number of data in file 1 and 2, min(ndata1,ndata2)
      integer idata 
      double precision norm1,norm2 ! norms of the two fields, |f1|=sqrt(<f1|f1>)
      ! if the rgrids are slightly different, use the smaller grid as reference:
      double precision, allocatable :: rgrid1a(:),field1a(:),
     &   field2a(:),dumr(:) 
      double precision ov ! overlap
      character(len=256) file1a ! debug
      character(len=256) file2a ! debug
      character(len=256) file12 ! debug
      !
      if(talk) print fsubstart,"compare" 
      if (talk) print '(8x,"first file: ",A40)',file1
      if (talk) print '(8x,"second file: ",A40)',file2
      file12=" "
      file12=file1
      !
      ! read grid and field from first file
      open(51,file=file1,status="old",err=101)
      ndata1=0
 10   read(51,'(A256)',end=11,err=101) line
      line=adjustl(line)
      if(line(1:1).ne."#") ndata1=ndata1+1
      goto 10
 11   rewind(51)
      allocate(rgrid1(ndata1),field1(ndata1))
      idata=1
 12   read(51,'(A256)',end=13,err=101) line     
      line=adjustl(line)
      if(line(1:1).ne."#") then
        read(line(1:256),*) rgrid1(idata),field1(idata)
        idata=idata+1
      end if
      goto 12
 13   close(51)
      ! read grid and field from second file
      open(51,file=file2,status="old",err=102)
      ndata2=0
 20   read(51,'(A256)',end=21,err=102) line
      line=adjustl(line)
      if(line(1:1).ne."#") ndata2=ndata2+1
      goto 20
 21   rewind(51)
      allocate(rgrid2(ndata2),field2(ndata2))
      idata=1
 22   read(51,'(A256)',end=23,err=102) line     
      line=adjustl(line)
      if(line(1:1).ne."#") then
        read (line,*) rgrid2(idata),field2(idata)
        idata=idata+1
      end if
      goto 22
 23   close(51)
      !
      ! if rgrids are different, reshuffle grids and fields i.o.t. use nearest grid points
      ndata=min(ndata1,ndata2)
      allocate(rgrid1a(ndata),field1a(ndata),
     &         field2a(ndata))
      if (any(rgrid1.ne.rgrid2)) then
        nwarn=nwarn+1 
        print '(8x,"grids are not equal, using nearest points")'
        if (ndata==ndata1) then
          rgrid1a=rgrid1
          field1a=field1 
          do idata=1,ndata 
            field2a(idata)=field2(ithenearest(rgrid1a(idata),rgrid2))
          end do 
        else
          rgrid1a=rgrid2
          field1a=field2 
          do idata=1,ndata 
            field2a(idata)=field1(ithenearest(rgrid1a(idata),rgrid1))
          end do 
          allocate(dumr(ndata))
          dumr=field1a
          field1a=field2a
          field2a=dumr
          deallocate(dumr)
        end if
      else
        rgrid1a=rgrid1
        field1a=field1
        field2a=field2
      end if
      deallocate(rgrid1,rgrid2,field1,field2)
      !
      ! begin write fields to disk (debug)
      ! first field
      write(file1a,'(a,a)') adjustl(trim(file1)),"a" 
      open(51,file=trim(file1a),status="replace")
      do idata=1,ndata
        write(51,'(F10.6,F15.6)') rgrid1a(idata),field1a(idata)       
      end do 
      close(51)
      ! second field
      write(file2a,'(a,a)') adjustl(trim(file2)),"a" 
      open(51,file=trim(file2a),status="replace")
      do idata=1,ndata
        write(51,'(F10.6,F15.6)') rgrid1a(idata),field2a(idata)       
      end do 
      close(51)
      ! difference
      write(file12,'(a,a)') adjustl(trim(file12)),"diff" 
      open(51,file=trim(file12),status="replace")
      do idata=1,ndata
        write(51,'(F10.6,F15.6)') rgrid1a(idata),                         &
     &     field1a(idata)-field2a(idata)       
      end do 
      close(51)
      ! end write fields to disk (debug)
      !
      ! calculate norms
      !norm1=sqrt(int_1d_taylorn(rgrid1a,field1a**2,6) ) ! Taylor interpolation leads to numerical errors (norm>1)
      !norm2=sqrt(int_1d_taylorn(rgrid1a,field2a**2,6) ) ! Taylor interpolation leads to numerical errors (norm>1) 
      norm1=sqrt(int_1d_a(rgrid1a,field1a**2) ) 
      norm2=sqrt(int_1d_a(rgrid1a,field2a**2) ) 
      if(talk) print '(8x,"Norm of first field: ",F20.6)', norm1
      if(talk) print '(8x,"Norm of second field: ",F20.6)', norm2
      !
      ! calculate overlap
      !ov=int_1d_taylorn(rgrid1a,field1a*field2a,6)/(norm1*norm2)
      ov=int_1d_a(rgrid1a,field1a*field2a)/(norm1*norm2)
      deallocate(rgrid1a,field1a,field2a)
      ! relative overlap
      if(talk) print '(8x,"Overlap: ",F20.6)', ov
      ! absolute overlap wrt first field:
      if(talk) print '(8x,"Abs. overlap 1: ",F20.6)', ov*norm2/norm1
      ! absolute overlap wrt second field:
      if(talk) print '(8x,"Abs. overlap 2: ",F20.6)', ov*norm1/norm2
      !
      ! end normally:
      if(talk) print fsubendext, "compare"
      return
      
      ! Errors:
 101  nerr=nerr+1
      print ferrmssg," when opening/reading file1"
      close(51)
      return
 102  nerr=nerr+1
      print ferrmssg," when opening/reading file2"
      close(51)
      return
  
      end subroutine compare

c---------------------------------------------------------------------
c     
c      subroutine madelen(lattice,charge,epsil,alat,talk,nerr)
c      use defs
c      implicit none
c      ! calculates the madelung energy for some point charge lattices in
c      ! a homogeneous background charge  
c      character, intent(in) :: lattice*7 ! lattice type
c      double precision, intent(in) :: alat(1:3),charge,epsil ! lattice constant, point charge, permittivity
c      logical, intent(in) :: talk
c      integer :: nerr
c      ! internal variables
c      double precision madel,R,G,vm,dr,c_over_a,vol
c      double precision, parameter :: beta=1.0D0
c      integer i,j,k,l,m,n,imax,maxl
c     
c      if(talk) then 
c        print fsubstart,'madelen'
c        print '("Calculating Madelung energy for ",A7," lattice with poi
c     &nt charge ",F8.4,"e, epsilon=" , F12.6,", and lattice constant",
c     &   3(F10.6)," Bohr" )',lattice,charge,epsil,alat
c      end if
c      select case(lattice)
c      case('CUB','cub','sc','SC')
c        imax=500
c        maxl=500
c        madel=0.0D0
c        do i=-imax,imax
c        do j=-imax,imax
c        do k=-imax,imax
c              if (abs(i)+abs(j)+abs(k).gt.0) then
c                    R=sqrt(float(i**2+j**2+k**2))
c                    madel=madel+erfc(sqrt(beta)*R)/R
c              end if
c        end do
c        end do
c        end do
c        do l=-maxl,maxl
c        do m=-maxl,maxl
c        do n=-maxl,maxl
c              if (abs(l)+abs(m)+abs(n).gt.0) then
c                    G=2.0D0*pi*sqrt(float(l**2+m**2+n**2))
c                    madel=madel+4.0D0*pi*exp(-0.250D0*G**2/beta)/(G**2)
c              end if
c        end do
c        end do
c        end do
c        madel=madel-2.0D0*sqrt(beta/pi)-pi/beta
c        madel=-madel
c
c        print*,"Madelung constant=",madel
c        print '("Madelung energy in Ryd=",F20.10)',
c     &     -madel*charge**2/(epsil*alat(1))
c        print '("Madelung energy in eV=",F20.10)',
c     &     -madel*charge**2*Ryd/(epsil*alat(1))
c
c        !open(10,file="erfc.dat",status="replace")
c        !do i=-imax,imax
c        !            write(10,*)float(i),erf(float(i)),erfc(float(i))
c        !end do
c        !close(10)
c
c        !vm=0.0D0
c        !do i=0,10000
c        !      dr=0.010D0
c        !      R=float(i)*dr
c        !      vm=vm+R*erfc(R)*dr
c        !end do
c        !print*, "V_av=",vm
c      case('fcc','FCC')
c        imax=500
c        maxl=500
c        madel=0.0D0
c        do i=-imax,imax
c        do j=-imax,imax
c        do k=-imax,imax
c              if (abs(i)+abs(j)+abs(k).gt.0) then
c                    R=0.50D0*sqrt(float((i+j)**2+(j+k)**2+(i+k)**2))
c                    madel=madel+erfc(sqrt(beta)*R)/R
c              end if
c        end do
c        end do
c        end do
c        do l=-maxl,maxl
c        do m=-maxl,maxl
c        do n=-maxl,maxl
c              if (abs(l)+abs(m)+abs(n).gt.0) then
c                    G=2.0D0*pi*sqrt(float((l-m+n)**2+(-l+m+n)**2+(l+m-n)
c     &                **2))
c                    madel=madel+16.0D0*pi*exp(-0.250D0*G**2/beta)/(G**2)
c              end if
c        end do
c        end do
c        end do
c        madel=madel-2.0D0*sqrt(beta/pi)-4.0D0*pi/(beta)
c        madel=-madel
c
c        print*,"Madelung constant=",madel
c        print '("Madelung energy in Ryd=",F20.10)',
c     &     -madel*charge**2/(epsil*alat(1))
c        print '("Madelung energy in eV=",F20.10)',
c     &     -madel*charge**2*Ryd/(epsil*alat(1))
c
c      case('bcc','BCC')
c        imax=500
c        maxl=500
c        madel=0.0D0
c        do i=-imax,imax
c        do j=-imax,imax
c        do k=-imax,imax
c              if (abs(i)+abs(j)+abs(k).gt.0) then
c                    R=0.50D0*sqrt(float((i+j-k)**2+(i-j+k)**2+(-i+j+k)
c     &                **2))
c                    madel=madel+erfc(sqrt(beta)*R)/R
c              end if
c        end do
c        end do
c        end do
c        do l=-maxl,maxl
c        do m=-maxl,maxl
c        do n=-maxl,maxl
c              if (abs(l)+abs(m)+abs(n).gt.0) then
c                    G=2.0D0*pi*sqrt(float((l+m)**2+(l+n)**2+(m+n)
c     &                **2))
c                    madel=madel+8.0D0*pi*exp(-0.250D0*G**2/beta)/(G**2)
c              end if
c        end do
c        end do
c        end do
c        madel=madel-2.0D0*sqrt(beta/pi)-2.0D0*pi/(beta)
c        madel=-madel
c
c        print*,"Madelung constant=",madel
c        print '("Madelung energy in Ryd=",F20.10)',
c     &     -madel*charge**2/(epsil*alat(1))
c        print '("Madelung energy in eV=",F20.10)',
c     &     -madel*charge**2*Ryd/(epsil*alat(1))
c
c      case('tet','TET','tetra','TETRA')
c        c_over_a=alat(3)/(alat(1))
c        vol=c_over_a
c        imax=500
c        maxl=500
c        madel=0.0D0
c
c        do i=-imax,imax
c        do j=-imax,imax
c        do k=-imax,imax
c              if (abs(i)+abs(j)+abs(k).gt.0) then
c                    R=sqrt(float(i**2+j**2)+c_over_a**2*float(k**2))
c                    madel=madel+erfc(sqrt(beta)*R)/R
c              end if
c        end do
c        end do
c        end do
c        do l=-maxl,maxl
c        do m=-maxl,maxl
c        do n=-maxl,maxl
c              if (abs(l)+abs(m)+abs(n).gt.0) then
c                    G=2.0D0*pi*sqrt(float(l**2+m**2)+float(n**2)
c     &                /c_over_a**2)
c                    madel=madel+4.0D0*pi*exp(-0.250D0*G**2/beta)
c     &                    /(vol*G**2)
c              end if
c        end do
c        end do
c        end do
c        madel=madel-2.0D0*sqrt(beta/pi)-pi/(beta*vol)
c        madel=-madel
c        print*,"Madelung constant=",madel
c        print '("Madelung energy in Ryd=",F20.10)',
c     &     -madel*charge**2/(epsil*alat(1))
c        print '("Madelung energy in eV=",F20.10)',
c     &     -madel*charge**2*Ryd/(epsil*alat(1))
c      case default
c        goto 100
c      end select
c      if(talk) then 
c        print fsubendext,'madelen'
c      end if
c      return
c! errors
c100   nerr=nerr+1
c      print ferrmssg, "Sorry, lattice type not implemented"
c      return
c
c      end subroutine
c
cc---------------------------------------------------------------------
c     
c      subroutine dipoleen0(lattice,charge,dipole,epsil,alat,talk,nerr)
c      use defs
c      implicit none
c      ! calculates the dipole-dipole energy for some point charge lattices in
c      ! a homogeneous background charge  
c      character, intent(in) :: lattice*7 ! lattice type
c      double precision, intent(in) :: alat(1:3),charge,dipole(1:3),epsil! lattice constant, point charge, dipole vector (d/q), permittivity
c      logical, intent(in) :: talk
c      integer :: nerr
c      ! internal variables
c      double precision dipen,R,G,c_over_a,vol,rvec(1:3),dip,prod ! dipole-dipole energy, direct and recip. lattice vector length, ... , ... , lattice vector, |dipole|, scalar product dipole*rvec
c      double precision, parameter :: beta=1.0D0
c      integer i,j,k,imax
c     
c      if(talk) then 
c        print fsubstart,'dipoleen0'
c        print '("Calculating dipole-dipole energy for ",A7," lattice wit
c     &h point charge ",F8.4,"e, dipole vector (d/q)=",3(F12.6), " Bohr",
c     &   ",epsilon=" , F12.6,", and lattice constant",
c     &   3(F10.6)," Bohr" )',lattice,charge,dipole,epsil,alat
c      end if
c      select case(lattice)
c      case('CUB','cub','sc','SC')
c        imax=500
c        dipen=0.0D0
c        do i=-imax,imax
c        do j=-imax,imax
c        do k=-imax,imax
c              if (abs(i)+abs(j)+abs(k).gt.0) then
c                    R=sqrt(dble(i**2+j**2+k**2))
c                    dip=sqrt(dot_product(dipole,dipole))
c                    rvec(1)=dble(i)
c                    rvec(2)=dble(j)
c                    rvec(3)=dble(k)
c                    prod=3.0d0*dot_product(dipole,rvec)**2
c     &              /(dip**2 * R**2)
c                    dipen=dipen+(1.0d0-prod)/R**3
c              end if
c        end do
c        end do
c        end do
c        print*,"Dipole-dipole sum=",dipen
c        print '("Dipole-dipole energy in Ryd=",F20.10)',
c     &     dipen*dip**2*charge**2 / (epsil * alat(1)**3 )
c        print '("Dipole-dipole energy in eV=",F20.10)',
c     &     dipen*dip**2*charge**2 * Ryd / (epsil * alat(1)**3 )
c
c      case('fcc','FCC')
c        imax=500
c        dipen=0.0D0
c        do i=-imax,imax
c        do j=-imax,imax
c        do k=-imax,imax
c              if (abs(i)+abs(j)+abs(k).gt.0) then
c                    R=0.50D0*sqrt(float((i+j)**2+(j+k)**2+(i+k)**2))
c                    dip=sqrt(dot_product(dipole,dipole))
c                    rvec(1)=0.5d0*dble(j+k)
c                    rvec(2)=0.5d0*dble(i+k)
c                    rvec(3)=0.5d0*dble(i+j)
c                    prod=3.0d0*dot_product(dipole,rvec)**2
c     &              /(dip**2 * R**2)
c                    dipen=dipen+(1.0d0-prod)/R**3
c              end if
c        end do
c        end do
c        end do
c        print*,"Dipole-dipole sum=",dipen
c        print '("Dipole-dipole energy in Ryd=",F20.10)',
c     &     dipen*dip**2*charge**2 / (epsil * alat(1)**3 )
c        print '("Dipole-dipole energy in eV=",F20.10)',
c     &     dipen*dip**2*charge**2 * Ryd / (epsil * alat(1)**3 )
c
c      case('bcc','BCC')
c        imax=500
c        dipen=0.0D0
c        do i=-imax,imax
c        do j=-imax,imax
c        do k=-imax,imax
c              if (abs(i)+abs(j)+abs(k).gt.0) then
c                    R=0.50D0*sqrt(float((i+j-k)**2+(i-j+k)**2+(-i+j+k)
c     &                **2))
c                    dip=sqrt(dot_product(dipole,dipole))
c                    rvec(1)=0.5d0*dble(-i+j+k)
c                    rvec(2)=0.5d0*dble(i-j+k)
c                    rvec(3)=0.5d0*dble(i+j-k)
c                    prod=3.0d0*dot_product(dipole,rvec)**2
c     &              /(dip**2 * R**2)
c                    dipen=dipen+(1.0d0-prod)/R**3
c              end if
c        end do
c        end do
c        end do
c        print*,"Dipole-dipole sum=",dipen
c        print '("Dipole-dipole energy in Ryd=",F20.10)',
c     &     dipen*dip**2*charge**2 / (epsil * alat(1)**3 )
c        print '("Dipole-dipole energy in eV=",F20.10)',
c     &     dipen*dip**2*charge**2 * Ryd / (epsil * alat(1)**3 )
c
c      case('tet','TET','tetra','TETRA')
c        c_over_a=alat(3)/alat(1)
c        vol=c_over_a
c        imax=500
c        dipen=0.0D0
c        do i=-imax,imax
c        do j=-imax,imax
c        do k=-imax,imax
c              if (abs(i)+abs(j)+abs(k).gt.0) then
c                    R=sqrt(dble(i**2+j**2)+c_over_a**2*dble(k**2))
c                    dip=sqrt(dot_product(dipole,dipole))
c                    rvec(1)=dble(i)
c                    rvec(2)=dble(j)
c                    rvec(3)=c_over_a*dble(k)
c                    prod=3.0d0*dot_product(dipole,rvec)**2
c     &              /(dip**2 * R**2)
c                    dipen=dipen+(1.0d0-prod)/R**3
c              end if
c        end do
c        end do
c        end do
c        print*,"Dipole-dipole sum=",dipen
c        print '("Dipole-dipole energy in Ryd=",F20.10)',
c     &     dipen*dip**2*charge**2 / (epsil * alat(1)**3 )
c        print '("Dipole-dipole energy in eV=",F20.10)',
c     &     dipen*dip**2*charge**2 * Ryd / (epsil * alat(1)**3 )
c      case default
c        goto 100
c      end select
c      if(talk) then 
c        print fsubendext,'dipoleen0'
c      end if
c      return
c! errors
c100   nerr=nerr+1
c      print ferrmssg, "Sorry, lattice type not implemented"
c      return
c
c      end subroutine
c
c
cc---------------------------------------------------------------------
c     
c      subroutine dipoleen(lattice,charge,pvec,epsil,alat,talk,nerr)
c      use defs
c      implicit none
c      ! calculates the dipole-dipole energy for some point charge lattices in
c      ! a homogeneous background charge  
c      character, intent(in) :: lattice*7 ! lattice type
c      double precision, intent(in) :: alat(1:3),charge,pvec(1:3),
c     &                                epsil(1:3,1:3) ! lattice constant, point charge, dipole vector (d/q), permittivity
c      logical, intent(in) :: talk
c      integer :: nerr
c      ! internal variables
c      ! dipole-dipole energy, direct and recip. lattice vector length, ... , ... , lattice vector, |dipole|, scalar product dipole*rvec:
c      double precision dipen,rabs,G,c_over_a,vol,rvec(1:3),pabs
c      double precision epsil1(1:3,1:3) ! inverse epsilon tensor 
c      double precision p_epsil1_p ! dipole moment times inverse epsilon tensor times dipole moment 
c      double precision p_epsil1_r ! dipole moment times inverse epsilon tensor times lattice vector 
c      double precision p_r ! dipole moment times lattice vector 
c      double precision p_epsil1(1:3) ! dipole moment times inverse epsilon tensor  
c      !double precision, parameter :: beta=1.0D0 ! Ewald sum parameter
c      integer i,j,k,imax
c     
c      ! print the input parameters
c      if(talk) then 
c        print fsubstart,'dipoleen'
c        print '("Calculating dipole-dipole energy for ",A7," lattice wit
c     &h point charge ",F8.4,"e, dipole vector (d/q)=",3(F12.6), " Bohr",
c     &   //,",epsilon=" , 9(F12.6),", and lattice constant",
c     &   3(F10.6)," Bohr" )',lattice,charge,pvec,epsil,alat
c      end if
c      
c      ! invert epsilon
c      call invert_matrix(epsil,epsil1)
c      ! some abbreviations
c      p_epsil1=matmul(pvec,epsil1)
c      pabs=sqrt(dot_product(pvec,pvec))
c      p_epsil1_p=dot_product(p_epsil1,pvec)/pabs**2
c      
c      select case(lattice)
c      case('CUB','cub','sc','SC')
c        imax=500
c        dipen=0.0D0
c        do i=-imax,imax
c        do j=-imax,imax
c        do k=-imax,imax
c              if (abs(i)+abs(j)+abs(k).gt.0) then
c                    rvec(1)=dble(i)
c                    rvec(2)=dble(j)
c                    rvec(3)=dble(k)
c                    rabs=sqrt(dot_product(rvec,rvec))
c                    p_r=dot_product(pvec,rvec)/(pabs*rabs)
c                    p_epsil1_r=dot_product(p_epsil1,rvec)/(pabs * rabs)
c                    dipen=dipen+(p_epsil1_p-3.0d0*p_epsil1_r*p_r)
c     &                    /rabs**3
c              end if
c        end do
c        end do
c        end do
c        print*,"Dipole-dipole sum=",dipen
c        print '("Dipole-dipole energy in Ryd=",F20.10)',
c     &     dipen*pabs**2*charge**2 / alat(1)**3 
c        print '("Dipole-dipole energy in eV=",F20.10)',
c     &     dipen*pabs**2*charge**2 * Ryd /  alat(1)**3 
c
c      case('fcc','FCC')
c        imax=500
c        dipen=0.0D0
c        do i=-imax,imax
c        do j=-imax,imax
c        do k=-imax,imax
c              if (abs(i)+abs(j)+abs(k).gt.0) then
c                    rvec(1)=0.5d0*dble(j+k)
c                    rvec(2)=0.5d0*dble(i+k)
c                    rvec(3)=0.5d0*dble(i+j)
c                    rabs=sqrt(dot_product(rvec,rvec))
c                    p_r=dot_product(pvec,rvec)/(pabs*rabs)
c                    p_epsil1_r=dot_product(p_epsil1,rvec)/(pabs * rabs)
c                    dipen=dipen+(p_epsil1_p-3.0d0*p_epsil1_r*p_r)
c     &                    /rabs**3
c              end if
c        end do
c        end do
c        end do
c        print*,"Dipole-dipole sum=",dipen
c        print '("Dipole-dipole energy in Ryd=",F20.10)',
c     &     dipen*pabs**2*charge**2 / alat(1)**3
c        print '("Dipole-dipole energy in eV=",F20.10)',
c     &     dipen*pabs**2*charge**2 * Ryd / alat(1)**3 
c
c      case('bcc','BCC')
c        imax=500
c        dipen=0.0D0
c        do i=-imax,imax
c        do j=-imax,imax
c        do k=-imax,imax
c              if (abs(i)+abs(j)+abs(k).gt.0) then
c                    rvec(1)=0.5d0*dble(-i+j+k)
c                    rvec(2)=0.5d0*dble(i-j+k)
c                    rvec(3)=0.5d0*dble(i+j-k)
c                    rabs=sqrt(dot_product(rvec,rvec))
c                    p_r=dot_product(pvec,rvec)/(pabs*rabs)
c                    p_epsil1_r=dot_product(p_epsil1,rvec)/(pabs * rabs)
c                    dipen=dipen+(p_epsil1_p-3.0d0*p_epsil1_r*p_r)
c     &                    /rabs**3
c              end if
c        end do
c        end do
c        end do
c        print*,"Dipole-dipole sum=",dipen
c        print '("Dipole-dipole energy in Ryd=",F20.10)',
c     &     dipen*pabs**2*charge**2 / alat(1)**3 
c        print '("Dipole-dipole energy in eV=",F20.10)',
c     &     dipen*pabs**2*charge**2 * Ryd / alat(1)**3 
c
c      case('tet','TET','tetra','TETRA')
c        c_over_a=alat(3)/alat(1)
c        vol=c_over_a
c        imax=500
c        dipen=0.0D0
c        do i=-imax,imax
c        do j=-imax,imax
c        do k=-imax,imax
c              if (abs(i)+abs(j)+abs(k).gt.0) then
c                    rvec(1)=dble(i)
c                    rvec(2)=dble(j)
c                    rvec(3)=c_over_a*dble(k)
c                    rabs=sqrt(dot_product(rvec,rvec))
c                    p_r=dot_product(pvec,rvec)/(pabs*rabs)
c                    p_epsil1_r=dot_product(p_epsil1,rvec)/(pabs * rabs)
c                    dipen=dipen+(p_epsil1_p-3.0d0*p_epsil1_r*p_r)
c     &                    /rabs**3
c              end if
c        end do
c        end do
c        end do
c        print*,"Dipole-dipole sum=",dipen
c        print '("Dipole-dipole energy in Ryd=",F20.10)',
c     &     dipen*pabs**2*charge**2 / alat(1)**3 
c        print '("Dipole-dipole energy in eV=",F20.10)',
c     &     dipen*pabs**2*charge**2 * Ryd / alat(1)**3 
c      case default
c        goto 100
c      end select
c      if(talk) then 
c        print fsubendext,'dipoleen'
c      end if
c      return
c! errors
c100   nerr=nerr+1
c      print ferrmssg, "Sorry, lattice type not implemented"
c      return
c
c      end subroutine


c---------------------------------------------------------------------

      double precision function bessels(l,x)
      ! spherical Bessel function jl(x), l=0 ... 3
      implicit none
      integer, intent(in) :: l   
      double precision, intent(in) :: x
      double precision fl,gl

      bessels=0.0D0
      if(abs(x).gt.1.E-10) then
        select case(l)
        case(0)
                fl=1.0D0/x
                gl=0.0D0
        case(1)
                fl=1.0D0/x**2
                gl=1.0D0/x
        case(2)        
                fl=3.0D0/x**3-1.0D0/x
                gl=3.0D0/x**2
        case(3)
                fl=15.0D0/x**4-6.0D0/x**2
                gl=15.0D0/x**3-1.0D0/x
        end select
        bessels=fl*sin(x)-gl*cos(x)
      else  
        select case(l)
        case(0)
                bessels=1.0d0
        case default
                bessels=0.0d0
        end select
      end if
      end function

c---------------------------------------------------------------------

      double precision function legendre(l,z)
      ! Legendre polynome Pl(z), l=0 ... 3
      implicit none
      integer, intent(in) :: l   
      double precision, intent(in) :: z

      legendre=0.0D0
      select case(l)
      case(0)
              legendre=1.0D0
      case(1)
              legendre=z
      case(2)        
              legendre=0.50D0*(3.0D0*z**2-1.0D0)
      case(3)
              legendre=0.50D0*(5.0D0*z**3-3.0D0*z)
      end select

      end function

c---------------------------------------------------------------------
C
C        integer function facult(n)
C        implicit none
C        integer, intent(in) :: n
C        integer i
C
C        facult=1
C        do i=2,n
C          facult=facult*i
C        end do
C
C        end function
C
c---------------------------------------------------------------------
C
C        integer function dblfactorial(n)
C        implicit none
C        integer, intent(in) :: n
C        integer i
C
C        dblfactorial=1
C        do i=3,n,2
C          dblfactorial=dblfactorial*i
C        end do
C
C        end function
C
c---------------------------------------------------------------------

      double precision function scalar_product(vec1,vec2)
      
      use defs
      implicit none
      
      double precision, intent(in) :: vec1(:),vec2(:)
      integer :: i 

      scalar_product=0.0D0
      if(size(vec1).ne.size(vec2)) then
         print ferrmssg,"Scalar prod. needs 2 vecs with same dim"
         stop
      end if

      do i=1,size(vec1)
            scalar_product=scalar_product+vec1(i)*vec2(i)
      end do

      end function

c---------------------------------------------------------------------

      function cross_product(vec1,vec2)
      
      use defs
      implicit none
      
      double precision, DIMENSION(3) :: cross_product 
      !double precision cross_product(1:3)  
      double precision, intent(in) :: vec1(3),vec2(3)
      integer :: i 

      if(size(vec1).ne.3.or.size(vec2).ne.3) then
         print ferrmssg,"Cross prod. needs 2 vecs with dim 3"
         return
      end if

      cross_product(1)=vec1(2)*vec2(3)-vec1(3)*vec2(2)
      cross_product(2)=vec1(3)*vec2(1)-vec1(1)*vec2(3)
      cross_product(3)=vec1(1)*vec2(2)-vec1(2)*vec2(1)

      end function cross_product

c---------------------------------------------------------------------

      subroutine vecs2vol(vecs,vol)
      ! calculate cell volume from vectors
      implicit none
      double precision, intent(in) :: vecs(1:3,1:3)
      double precision vol
      vol=dot_product(vecs(1,1:3),cross_product(vecs(2,:),vecs(3,:)))
      vol=abs(vol)
      end subroutine vecs2vol

c---------------------------------------------------------------------

      subroutine recipr_latt_vecs(vecs,vecs2)
      ! calculates reciprocal lattice vectors (vecs2) from direct
      ! lattice vectors (vecs)
      use defs, only: pi
      implicit none
      double precision vol
      double precision vecs(1:3,1:3),vecs2(1:3,1:3)
      call vecs2vol(vecs,vol)
      ! calculate reciprocal lattice vectors as cross products
      vecs2(1,1)=vecs(2,2)*vecs(3,3)-vecs(2,3)*vecs(3,2)
      vecs2(1,2)=vecs(2,3)*vecs(3,1)-vecs(2,1)*vecs(3,3)
      vecs2(1,3)=vecs(2,1)*vecs(3,2)-vecs(2,2)*vecs(3,1)

      vecs2(2,1)=vecs(3,2)*vecs(1,3)-vecs(3,3)*vecs(1,2)
      vecs2(2,2)=vecs(3,3)*vecs(1,1)-vecs(3,1)*vecs(1,3)
      vecs2(2,3)=vecs(3,1)*vecs(1,2)-vecs(3,2)*vecs(1,1)

      vecs2(3,1)=vecs(1,2)*vecs(2,3)-vecs(1,3)*vecs(2,2)
      vecs2(3,2)=vecs(1,3)*vecs(2,1)-vecs(1,1)*vecs(2,3)
      vecs2(3,3)=vecs(1,1)*vecs(2,2)-vecs(1,2)*vecs(2,1)

      vecs2=vecs2*2.0d0*pi/vol

      end subroutine recipr_latt_vecs

c---------------------------------------------------------------------      

      integer function ithenearest(value,array)
      ! returns the index of the "array" element closest to "value" 
      implicit none
      double precision, intent(in) :: value, array(:)
      ! internal 
      integer i
      double precision diff
      !
      ithenearest=1
      diff=abs(array(1)-value)
      do i=2,size(array) 
        if (abs(array(i)-value).lt.diff) then
          diff=abs(array(i)-value)
          ithenearest=i
        end if
      end do 
      !      
      end function ithenearest 

c---------------------------------------------------------------------

      integer function inearestlower(value,array)
      ! returns the index of the "array" element closest below "value" 
      implicit none
      double precision, intent(in) :: value, array(:)
      ! internal 
      integer i
      double precision diff
      !
      inearestlower=1
      diff=abs(array(1)-value)
      do i=2,size(array) 
        if (abs(array(i)-value).lt.diff) then
          diff=abs(array(i)-value)
          inearestlower=i
        end if
      end do 
      if(array(inearestlower).gt.value) inearestlower=inearestlower-1
      !      
      end function inearestlower 

c---------------------------------------------------------------------

      integer function inearestupper(value,array)
      ! returns the index of the "array" element closest above "value" 
      implicit none
      double precision, intent(in) :: value, array(:)
      ! internal 
      integer i
      double precision diff
      !
      inearestupper=1
      diff=abs(array(1)-value)
      do i=2,size(array) 
        if (abs(array(i)-value).lt.diff) then
          diff=abs(array(i)-value)
          inearestupper=i
        end if
      end do 
      if(array(inearestupper).gt.value) inearestupper=inearestupper+1
      !      
      end function inearestupper

c---------------------------------------------------------------------

      integer function countsubstring(string,substring)
      implicit none
      character (len=*), intent(in) :: substring, string  
      character (len=len(substring)) :: linepart
      integer i  
      countsubstring=0
      if (len(string)==0) return
      if (len(substring)>len(string)) return  
      i=1
      do while (i<=len(string))
        if (index(string(i:i+len(substring)-1),substring).le.0) i=i+1
        if (index(string(i:i+len(substring)-1),substring).gt.0) then
          i=i+len(substring)
          countsubstring=countsubstring+1
        end if
      end do  
      end function countsubstring 

c---------------------------------------------------------------------

      subroutine string2words(string0,words)
      implicit none
      character (len=*), intent(in) :: string0
      character (len=1024), allocatable :: words(:)
      integer i,iword,nwords,istart,iend  
      character (len=len_trim(string0)) :: string1
      character (len=1024) formatstring
      string1=trim(adjustl(string0))
      !print*,'string="',string1,'"'
      nwords=0
      if (len(string1)==0) return
      i=1
      nwords=1
      if (len(string1).gt.1) then  
        i=2
        do while (i<=len(string1))
          if (string1(i:i).eq.' '.and.string1(i-1:i-1).ne.' ') then 
                  nwords=nwords+1
          end if
          i=i+1
        end do  
      end if ! len(string)>1
      !print*,'nwords=',nwords
      if(allocated(words)) deallocate(words)
      allocate(words(nwords))
      ! begin read words
      istart=1
      iend=1
      iword=1
      do while(istart.le.len(string1))
        !print*,'istart=',istart
        do while (iend.le.len(string1).and.string1(iend:iend).ne.' ')
          !print*,'iend,string(iend):',iend,string1(iend:iend)
          iend=iend+1
        end do
        !print*,'istart,iend:',istart,iend
        formatstring=''
        !write(formatstring,'(A3,I0,A2)') "'(A",iend-istart,")'"
        write(formatstring,'("(A",I0,")")') iend-istart
        !print*,formatstring
        !read(string1(istart:iend-1),*) words(iword)
        read(string1(istart:iend-1),formatstring) words(iword)
        !print*,'iword,words(iword):',iword,words(iword)
        iword=iword+1
        istart=iend
        do while(string1(istart:istart).eq.' '.and.istart.lt.             &
     &           len(string1))
          istart=istart+1
        end do
        iend=istart
      end do  
      !print*,words
      ! end read words
      !print*,'finished'
      end subroutine string2words 

c---------------------------------------------------------------------
    
      integer function num_in_ele(atoms,iatom)
      use defs
      implicit none
      type(atom), intent(in) :: atoms(:)
      integer, intent(in) :: iatom
      ! local variables:
      integer jatom

      num_in_ele=1
      do jatom=1,iatom-1
        if (atoms(iatom)%name(1:2)==atoms(jatom)%name(1:2).and.           &
     &      atoms(iatom)%core==atoms(jatom)%core)                         &
     &     num_in_ele=num_in_ele+1
      end do
    
      end function num_in_ele
      
c---------------------------------------------------------------------

      integer function gcd(intnums)
      ! greates common divisor
      use defs
      implicit none

      integer, intent(in) :: intnums(:)
      integer numint,iint,jint, maxnum
      logical is_common_divisor
      
      numint=size(intnums)
      
      maxnum=maxval(intnums)
      gcd=1
      
      iint=2

      do while (iint.le.maxnum) ! loop over possible common divisors iint
        is_common_divisor=.true.
        do jint=1, numint ! loop over numbers in vector
          if (mod(intnums(jint),iint).gt.0) is_common_divisor=.false.
        end do 
        if (is_common_divisor) gcd=iint
        iint=iint+1
      end do

      end function gcd

c---------------------------------------------------------------------
 
      subroutine vasp_kgrid(nk,shift) 
      use defs, only: talk, fsubstart,fsubendext
      implicit none

      integer nk(1:3)
      double precision shift(1:3)
      integer i1,i2,i3

      if (talk) print fsubstart,'vasp_kgrid'

      open(51,file='KPOINTS.GRID',status='replace')
      write(51,'("Grid generated by cofima")')
      write(51,'(I10)') nk(1)*nk(2)*nk(3)
      write(51,'("Reciprocal lattice")')
      do i1=1,nk(1)
      do i2=1,nk(2)
      do i3=1,nk(3)
        write(51,'(4F12.6)')                                              &
     &     dble(i1)/dble(nk(1))-0.5d0-dble(mod(nk(1),2))                  &
     &             /dble(2.0d0*nk(1)),                                    &
     &     dble(i2)/dble(nk(2))-0.5d0-dble(mod(nk(2),2))                  &
     &             /dble(2.0d0*nk(2)),                                    &
     &     dble(i3)/dble(nk(3))-0.5d0-dble(mod(nk(3),2))                  &
     &             /dble(2.0d0*nk(3)),  1.0d0
      end do ! i3
      end do ! i2
      end do ! i1
      close(51)
      !
      if (talk) print fsubendext,'vasp_kgrid'
      return

      end subroutine vasp_kgrid

c---------------------------------------------------------------------
 
      subroutine vasp_kpath(nk,startk,endk) 
      use defs, only: talk, fsubstart,fsubendext
      implicit none

      integer, intent(in) :: nk
      double precision startk(1:3),endk(1:3)
      integer i1

      if (talk) print fsubstart,'vasp_kpath'

      open(51,file='KPOINTS.PATH',status='replace')
      write(51,'("file generated by cofima")')
      write(51,'(I10)') nk+1
      write(51,'("Reciprocal lattice")')
      do i1=0,nk
        write(51,'(4F12.6)')                                              &
     &     startk(:)+(endk(:)-startk(:))*dble(i1)/dble(nk), 0.0d0
      end do ! i1
      close(51)
      !
      if (talk) print fsubendext,'vasp_kpath'
      return

      end subroutine vasp_kpath

c---------------------------------------------------------------------
 
      subroutine vasp_kaux(fname,dist) 
      use defs, only: pi, talk, fsubstart,fsubendext,error_stop
      implicit none

      double precision k(1:3),kaux(1:3),weight,dist(1:3)
      double precision kvecs(1:3,1:3),kabs(1:3),fdum
      character(len=*) fname
      integer nk,i1
      character line*128
      logical file_exists

      if (talk) print fsubstart,'vasp_kaux'

      ! if OUTCAR exists, read reciprocal lattice vectors
      kvecs=0.0d0
      do i1=1,3
        kvecs(i1,i1)=1.0d0
      end do
      kvecs=kvecs*2.0d0*pi
      INQUIRE(FILE="OUTCAR", EXIST=file_exists)
      if (file_exists) then
        open(51,file="OUTCAR",status='old',err=12)
        rewind(51)
10      read(51,'(A128)',end=12,err=12) line
        if (index(line,"reciprocal lattice vectors").gt.0) then
          read(51,*,err=12,end=12) fdum,fdum,fdum, kvecs(1,:)
          read(51,*,err=12,end=12) fdum,fdum,fdum, kvecs(2,:)
          read(51,*,err=12,end=12) fdum,fdum,fdum, kvecs(3,:)
        end if ! found rec. latt. vecs
        goto 10
12      continue
        close(51)        
      end if ! file_exists

      open(51,file=fname,status='old')
      open(52,file='KPOINTS.AUX',status='replace')
      read(51,'(A128)') line
      write(52,'("file generated by cofima")')
      read(51,*) nk 
      write(52,'(I10)') nk*7
      read(51,'(A128)') line
      if (index(line,'r').lt.1.and.index(line,'R').lt.1) goto 100
      write(52,'("Reciprocal lattice")')
      do i1=1,nk
        read(51,*) k(1:3), weight
        write(52,'(4F12.6)') k(1:3),weight                                &
        !write(52,'(4F12.6)') k(1:3)+(/dist(1),0.0d0,0.0d0/),0.0d0
        !write(52,'(4F12.6)') k(1:3)+(/-dist(1),0.0d0,0.0d0/),0.0d0
        !write(52,'(4F12.6)') k(1:3)+(/0.0d0,dist(2),0.0d0/),0.0d0
        !write(52,'(4F12.6)') k(1:3)+(/0.0d0,-dist(2),0.0d0/),0.0d0
        !write(52,'(4F12.6)') k(1:3)+(/0.0d0,0.0d0,dist(3)/),0.0d0
        !write(52,'(4F12.6)') k(1:3)+(/0.0d0,0.0d0,-dist(3)/),0.0d0
      end do ! i1
      rewind(51)
      read(51,'(A128)') line
      read(51,'(A128)') line
      read(51,'(A128)') line
      do i1=1,nk
        read(51,*) k(1:3), weight
        !write(52,'(4F12.6)') k(1:3),weight                                &
        call frac2abs(k,kvecs,kabs)
        !write(52,'(4F12.6)') k(1:3)+(/dist(1),0.0d0,0.0d0/),0.0d0
        !write(52,'(4F12.6)') k(1:3)+(/-dist(1),0.0d0,0.0d0/),0.0d0
        !write(52,'(4F12.6)') k(1:3)+(/0.0d0,dist(2),0.0d0/),0.0d0
        !write(52,'(4F12.6)') k(1:3)+(/0.0d0,-dist(2),0.0d0/),0.0d0
        !write(52,'(4F12.6)') k(1:3)+(/0.0d0,0.0d0,dist(3)/),0.0d0
        !write(52,'(4F12.6)') k(1:3)+(/0.0d0,0.0d0,-dist(3)/),0.0d0
        call abs2frac(kabs+(/dist(1),0.0d0,0.0d0/),kvecs,kaux)
        write(52,'(4F12.6)') kaux,0.0d0
        call abs2frac(kabs+(/-dist(1),0.0d0,0.0d0/),kvecs,kaux)
        write(52,'(4F12.6)') kaux,0.0d0
        call abs2frac(kabs+(/0.0d0,dist(2),0.0d0/),kvecs,kaux)
        write(52,'(4F12.6)') kaux,0.0d0
        call abs2frac(kabs+(/0.0d0,-dist(2),0.0d0/),kvecs,kaux)
        write(52,'(4F12.6)') kaux,0.0d0
        call abs2frac(kabs+(/0.0d0,0.0d0,dist(3)/),kvecs,kaux)
        write(52,'(4F12.6)') kaux,0.0d0
        call abs2frac(kabs+(/0.0d0,0.0d0,-dist(3)/),kvecs,kaux)
        write(52,'(4F12.6)') kaux,0.0d0
      end do ! i1
      close(51)
      close(52)
      !
      if (talk) print fsubendext,'vasp_kaux'
      return

100   continue
      close(51)
      close(52)
      call error_stop('KPOINT file has unexpected format.')      

      end subroutine vasp_kaux

c---------------------------------------------------------------------
 
      subroutine vasp_kmerge(fname1,fname2) 
      use defs, only: talk, fsubstart,fsubendext,error_stop
      implicit none

      double precision k(1:3),weight
      character(len=*) fname1,fname2
      integer nk1,nk2,i1
      character line*128

      if (talk) print fsubstart,'vasp_kmerge'

      open(51,file=fname1,status='old',err=100)
      open(52,file=fname2,status='old',err=100)
      open(53,file='KPOINTS.MERGED',status='replace',err=100)
      read(51,'(A128)', err=100) line
      read(52,'(A128)', err=100) line
      write(53,'("file generated by cofima")')
      read(51,*) nk1
      read(52,*) nk2 
      write(53,'(I10)') nk1+nk2
      read(51,'(A128)') line
      if (index(line,'r').lt.1.and.index(line,'R').lt.1) goto 100
      read(52,'(A128)') line
      if (index(line,'r').lt.1.and.index(line,'R').lt.1) goto 100
      write(53,'("Reciprocal lattice")')
      do i1=1,nk1
        read(51,*) k(1:3), weight
        write(53,'(4F12.6)') k(1:3),weight                                &
      end do ! i1
      close(51)
      do i1=1,nk2
        read(52,*) k(1:3), weight
        write(53,'(4F12.6)') k(1:3),weight                                &
      end do ! i1
      close(52)
      close(53)
      !
      if (talk) print fsubendext,'vasp_kmerge'
      return

100   continue
      close(51)
      close(52)
      close(53)
      call error_stop('reading/writing files went wrong.')      

      end subroutine vasp_kmerge

!--------------------------------------------------------------------------   

      double precision function theta_function(x)
      double precision x
      theta_function=0.0d0
      if (x==0.0d0) theta_function=0.5d0
      if (x>0.0d0) theta_function=1.0d0
      end function theta_function

!--------------------------------------------------------------------------   

      double precision function fermi_dist(energy,sigma)

      double precision energy,sigma
      fermi_dist=1.0d0/(1.0d0+exp(energy/sigma))
c      if (energy.lt.-10.0d0*sigma) fermi_dist=1.0d0
c      if (energy.gt.10.0d0*sigma) fermi_dist=0.0d0
      
      end function fermi_dist

!--------------------------------------------------------------------------   

      double precision function delta_function_fermi(energy,sigma)

      double precision energy,sigma
      delta_function_fermi=exp(-energy/sigma)/(exp(-energy/sigma)+1.0d0)  &
     &       **2/sigma
c      if (energy.lt.-10.0d0*sigma) delta_function_fermi=0.0d0
c      if (energy.gt.10.0d0*sigma) delta_function_fermi=0.0d0
      
      end function delta_function_fermi

c---------------------------------------------------------------------

      double precision function delta_function_gaussian(energy,sigma)

      use defs, only: pi
      double precision energy,sigma
      delta_function_gaussian=exp(-(energy/sigma)**2)/sigma/sqrt(pi)
c      if (energy.lt.-10.0d0*sigma) delta_function_gaussian=0.0d0
c      if (energy.gt.10.0d0*sigma) delta_function_gaussian=0.0d0

      end function delta_function_gaussian

c---------------------------------------------------------------------

      double precision function gaussian_dist(energy,sigma)

      double precision energy,sigma
      gaussian_dist=0.5d0-0.5d0*erf(energy/sigma)
c      if (energy.lt.-10.0d0*sigma) gaussian_dist=1.0d0
c      if (energy.gt.10.0d0*sigma) gaussian_dist=0.0d0

      end function gaussian_dist

c---------------------------------------------------------------------


      end module
