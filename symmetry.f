        module symmetry
        implicit none
        contains

c---------------------------------------------------------------------

      subroutine checksymgroup(symops,isymops,nsymm)
      use defs
      use misc  
      implicit none

      ! global variables
      type(symop) symops(nsymm),isymops(nsymm)
      integer nsymm
      ! local variables
      double precision iavec(1:3),cmat(1:3,1:3),cvec(1:3)
      integer isymm,jsymm,ksymm
      logical isgroup, haveidentity,haveinverse,havethisinverse,
     &        haveproduct,havethisproduct 
      double precision identity(1:3,1:3)

      ! initialize
      write(12,fsubstart) 'checksymgroup'
      write(12,'(8x,"Checking if the symmetry operations form a",
     &      " group...")')
      identity=0.0D0
      identity(1,1)=1.0D0
      identity(2,2)=1.0D0
      identity(3,3)=1.0D0
      isgroup=.true.
      haveidentity=.false.
      haveinverse=.true.
      haveproduct=.true.
      
      ! check if identity is included (neutral element)  
      do isymm=1,nsymm
        if (matnorm(symops(isymm)%mat-identity)
     &      /matnorm(symops(isymm)%mat).lt.1E-3.and.
     &      absvec(symops(isymm)%trala).lt.1E-3) 
     &     haveidentity=.true.
      end do 
      ! if identity is not there, the symmetry operations do not form a
      ! group.
      if (.not.haveidentity) isgroup=.false. 
      write(12,'(8x,"Identity included: ",L6)') haveidentity
      
      ! check for all symmetry operations A and B if A*B is included
      do isymm=1,nsymm
        do jsymm=1,nsymm
          havethisproduct=.false.
          cmat=matmul(symops(isymm)%mat,symops(jsymm)%mat)
          cvec=matmul(symops(isymm)%mat,symops(jsymm)%trala)
     &         +symops(isymm)%trala       
          do ksymm=1,nsymm
!            if     (symops(ksymm)%mat(1,1)==cmat(1,1)
!     &        .and.symops(ksymm)%mat(1,2)==cmat(1,2)
!     &        .and.symops(ksymm)%mat(1,3)==cmat(1,3)
!     &        .and.symops(ksymm)%mat(2,1)==cmat(2,1)
!     &        .and.symops(ksymm)%mat(2,2)==cmat(2,2)
!     &        .and.symops(ksymm)%mat(2,3)==cmat(2,3)
!     &        .and.symops(ksymm)%mat(3,1)==cmat(3,1)
!     &        .and.symops(ksymm)%mat(3,2)==cmat(3,2)
!     &        .and.symops(ksymm)%mat(3,3)==cmat(3,3)
!     &        .and. symops(ksymm)%trala(1)==cvec(1)
!     &        .and. symops(ksymm)%trala(2)==cvec(2)
!     &        .and. symops(ksymm)%trala(3)==cvec(3))
!     &        havethisproduct=.true.   
            if(matnorm(symops(ksymm)%mat-cmat)
     &         /matnorm(cmat).lt.1E-3
     &        .and.absvec(symops(ksymm)%trala-cvec)
     &        .lt.1E-3)
     &        havethisproduct=.true.
          end do  
          if (.not.havethisproduct) then
            haveproduct=.false.
            write(12,'(8x,"Product not included for SYMOPs ",I4,I4)') 
     &                 isymm,jsymm
          end if
        end do
      end do 
      if (.not.haveproduct) isgroup=.false.
      write(12,'(8x,"Products included: ",L6)') haveproduct

      ! check for all symmetry operations A if A^-1 is included (inverse
      ! element)
      do isymm=1,nsymm
        havethisinverse=.false.
        do jsymm=1,nsymm
          iavec=-matmul(isymops(isymm)%mat,symops(isymm)%trala)
!          if (symops(jsymm)%mat(1,1)==isymops(isymm)%mat(1,1)
!     &      .and.symops(jsymm)%mat(1,2)==isymops(isymm)%mat(1,2)
!     &      .and.symops(jsymm)%mat(1,3)==isymops(isymm)%mat(1,3)
!     &      .and.symops(jsymm)%mat(2,1)==isymops(isymm)%mat(2,1)
!     &      .and.symops(jsymm)%mat(2,2)==isymops(isymm)%mat(2,2)
!     &      .and.symops(jsymm)%mat(2,3)==isymops(isymm)%mat(2,3)
!     &      .and.symops(jsymm)%mat(3,1)==isymops(isymm)%mat(3,1)
!     &      .and.symops(jsymm)%mat(3,2)==isymops(isymm)%mat(3,2)
!     &      .and.symops(jsymm)%mat(3,3)==isymops(isymm)%mat(3,3)
!     &      .and. symops(jsymm)%trala(1)==iavec(1)
!     &      .and. symops(jsymm)%trala(2)==iavec(2)
!     &      .and. symops(jsymm)%trala(3)==iavec(3))
!     &      havethisinverse=.true.     
          if (matnorm(symops(jsymm)%mat-isymops(isymm)%mat)
     &        /matnorm(symops(jsymm)%mat).lt.1E-3     
     &        .and.absvec(symops(jsymm)%trala-iavec).lt.1E-3)
     &      havethisinverse=.true.     
        end do
        if (.not.havethisinverse) then
          haveinverse=.false.
          write(12,'(8x,"Inverse missing for SYMOP ",I4)') isymm
          write(12,'(8x,"iavec=",3(F7.2))') 
     &          -matmul(isymops(isymm)%mat,symops(isymm)%trala)
          write(12,'(8x,"matrix =",9(F7.2))') 
     &          symops(isymm)%mat(1,1:3),
     &          symops(isymm)%mat(2,1:3),
     &          symops(isymm)%mat(3,1:3)
          write(12,'(8x,"inverse=",9(F7.2))') 
     &          isymops(isymm)%mat(1,1:3),
     &          isymops(isymm)%mat(2,1:3),
     &          isymops(isymm)%mat(3,1:3)
        end if  
      end do 
      if (.not.haveinverse) isgroup=.false.
      write(12,'(8x,"Inverse element included: ",L6)') haveinverse

      ! if symmetry operations do not form a group, return an error
      ! message
      write(12,'(8x,"Symmetry operations form a group: ",L6)') isgroup
      if(.not.isgroup) goto 1000
      write(12,fsubend) 
      return

1000  write(12,ferrmssg) 'Symmetry operations do not form a group'    
      write(12,'(8x,"Make sure that your symmetry operations form a ")')
      write(12,'(8x,"group by adding or removing symmetry ")')
      write(12,'(8x,"operations. If not done yet, add the identity.")') 
      nerr=nerr+1  
      write(12,fsubend) 
      return

      end subroutine

c---------------------------------------------------------------------      
!
!      subroutine checksymmlat(symops,atoms,vecs,symlat)  
!      ! checks if the structure is compatible with the symmetry
!      use defs
!      use misc
!      implicit none
!      type(atom), intent(in) :: atoms(:)
!      type(symop), intent(in) :: symops(:)
!      double precision, intent(in) :: vecs(1:3,1:3)
!      !local
!      double precision xfrac(1:3),xabs(1:3),xtrans(1:3),
!     &     xtransfrac(1:3),distfrac(1:3)
!      integer natoms,nsymm,iatom,jatom,isymm,i
!      logical symlat,thissymlat  
!
!      write(12,fsubstart) 'checksymlat'
!      natoms=size(atoms)
!      nsymm=size(symops)
!      symlat=.true.
!
!      ! check if atoms get mapped onto same atoms by the 
!      ! symmetry operations
!      do iatom=1,natoms
!        do isymm=1,nsymm
!        call frac2abs(atoms(iatom)%where,vecs,xabs)
!          xtrans=matmul(symops(isymm)%mat,xabs)
!     &           +symops(isymm)%trala
!          call abs2frac(xtrans,vecs,xtransfrac)
!          thissymlat=.false.
!          do jatom=1,natoms
!            do i=1,3
!              distfrac(i)
!     &          =mod(abs(xtransfrac(i)-atoms(jatom)%where(i)),1.0D0)
!            end do
!            ! check if this atom is mapped onto another same atom: 
!            if(absvec(distfrac).lt.1E-3.and.atoms(iatom)%name
!     &         .eq.atoms(jatom)%name) thissymlat=.true.
!          end do
!        end do
!        ! if this atom is not mapped onto a same atom, the
!        ! symmetry operations are incompatible with the structure.
!        if(.not.thissymlat) symlat=.false.
!      end do
!
!      ! if symmetry operations incompatible with lattice, return an
!      ! error message
!      write(12,'(8x,"Symmetry operations compatible with structure:",
!     &      L4)') symlat
!      if (.not.symlat) goto 1000
!      write(12,fsubend) 
!      return
!
!1000  nerr= nerr+1     
!      write(12,ferrmssg)'Symmetries incompatible with structure'
!      write(12,fsubend) 
!      return
!
!      end subroutine
!
c---------------------------------------------------------------------      

      subroutine checksymmlat1(symops,atoms,vecs,symlat)  
      ! checks if the structure is compatible with the symmetry
      use defs
      use misc
      implicit none
      type(atom), intent(in) :: atoms(:)
      type(symop), intent(in) :: symops(:)
      double precision, intent(in) :: vecs(1:3,1:3)
      !local
      double precision xfrac(1:3),xabs(1:3),xtrans(1:3),
     &     xtransfrac(1:3),distfrac(1:3)
      integer natoms,nsymm,iatom,jatom,isymm,i
      logical symlat,thissymlat  

      if (talk) print fsubstart, 'checksymlat1'
      natoms=size(atoms)
      nsymm=size(symops)
      symlat=.true.

      ! check if atoms get mapped onto same atoms by the 
      ! symmetry operations
      do iatom=1,natoms
        xfrac=atoms(iatom)%where
        do isymm=1,nsymm
          xtransfrac=mod(matmul(symops(isymm)%mat,xfrac)                  &
     &           +symops(isymm)%trala,1.0d0)
          thissymlat=.false.
          do jatom=1,natoms
            do i=1,3
              distfrac(i)
     &          =mod(abs(xtransfrac(i)-atoms(jatom)%where(i)),1.0D0)
              do while (distfrac(i).le.-0.5d0)
                distfrac(i)=distfrac(i)+1.0d0  
              end do
              do while (distfrac(i).gt.0.5d0)
                distfrac(i)=distfrac(i)-1.0d0  
              end do
            end do ! i
            ! check if this atom is mapped onto another same atom: 
            if(absvec(distfrac).lt.1E-3.and.atoms(iatom)%name(1:2)        &
     &         .eq.atoms(jatom)%name(1:2)) thissymlat=.true.
!            ! begin debug
!            if(atoms(iatom)%name(1:2).eq.atoms(jatom)%name(1:2)) then
!              print '(8x,A3," atom ",I5," orig.: ",3(F8.3,x),"trans: ",   &
!     &         3(F8.3,x)," atom ",I5,": ",3(F8.3,x)," dist: ",F8.3)',     &
!     &         atoms(iatom)%name,iatom,atoms(iatom)%where,xtransfrac,     &
!     &         jatom,atoms(jatom)%where,absvec(distfrac)
!            end if
!            ! end debug
          end do ! jatom
!          print '(8x,"this atom is inv. under syms: ",L1)',thissymlat
        end do ! isymm
        ! if this atom is not mapped onto a same atom, the
        ! symmetry operations are incompatible with the structure.
        if(.not.thissymlat) symlat=.false.
      end do ! iatom

      if (talk) print fsubend
      return

      end subroutine checksymmlat1

c---------------------------------------------------------------------

      subroutine kreduce(symops,kmesh,kirred,kred)
        
      use defs
      use misc
      implicit none
      !double precision, dimension(:,:), intent(in) ::  kmesh
      type(kpoint), dimension(:), intent(inout) ::  kmesh
      type(symop), dimension (:), intent(in) :: symops(:)
      type(kpoint), allocatable :: kirred(:),kred(:)
      integer nkirred
      !local variables  
      type(kpoint), allocatable :: kpoints(:)
      integer ik,jk,lk,isymm,nsymm,nkpoints,nkred  

      write(12,fsubstart) 'kreduce'
      write(12,'(8x,"Looking for irreducible k-points/g-vectors...")')  
      nkpoints=size(kmesh)
      nsymm=size(symops)
      allocate (kpoints(nkpoints))
      ! check for each k-point if it is reducible 
      do ik=1,nkpoints
        !kpoints(ik)%kpt=kmesh(ik)%kpt(1:3)
        kpoints(ik)=kmesh(ik)
        kpoints(ik)%irred=.true.
        kmesh(ik)%irred=.true.
        kpoints(ik)%weight=1.0D0
        do jk=ik+1,nkpoints
          do isymm=1,nsymm
            ! k_j= Symop_i * k_i, then k_i reducible
            if(absvec(matmul(symops(isymm)%mat,kmesh(jk)%kpt(1:3))
     &        -kmesh(ik)%kpt(1:3))/absvec(kmesh(ik)%kpt(1:3))
     &        .lt.1E-3) then
               kpoints(ik)%irred=.false.    
               kmesh(ik)%irred=.false.
            end if   
          end do 
        end do  
      end do
      ! count reducible/irreducible k-points
      nkirred=0
      do ik=1,nkpoints
        if (kpoints(ik)%irred.eqv..true.) nkirred=nkirred+1
      end do
      nkred=nkpoints-nkirred
      allocate(kirred(nkirred))
      allocate(kred(nkred))
      jk=0
      lk=0
      do ik=1,nkpoints
        if (kpoints(ik)%irred) then
           jk=jk+1
           kirred(jk)=kpoints(ik)
        end if
        if (.not.kpoints(ik)%irred) then
           lk=lk+1
           kred(lk)=kpoints(ik)
        end if
      end do 
      ! count number of equivalent reduc. k-points for each irred. k-point  
      do ik=1,nkirred
        do jk=1,nkred
          do isymm=1,nsymm
            if(absvec(matmul(symops(isymm)%mat,kirred(ik)%kpt)
     &        -kred(jk)%kpt)/absvec(kred(jk)%kpt).lt.1E-3) 
     &        kirred(ik)%weight=kirred(ik)%weight+1.0D0      
          end do
        end do
      end do

      write(12,fsubend) 
      end subroutine

c---------------------------------------------------------------------

      subroutine gstar(girred,gred,symops)
      use defs
      use misc
      implicit none
      type(kpoint),intent(inout) ::girred(:),gred(:)
      type(symop),intent(in):: symops(:)
      !local
      integer nglatt,ngirred,ngred,ig,jg,isymm,nsymm,istar0
      double precision gtrans(1:3)               

      ! initializing
      write(12,fsubstart) 'gstar'
      write(12,'(8x,"Setting up stars of g-vectors...")')
      ngirred=size(girred)
      ngred=size(gred)  ! ngred=ngtotal-ngirred
      nsymm=size(symops)
      ! assign a star index to each irred. g-vector
      do ig=1,ngirred  
        girred(ig)%istar=ig
        girred(ig)%nstar=1
        girred(ig)%phase=0.0D0
      end do  
      do ig=1,ngred  
        gred(ig)%istar=0
        !gred(ig)%nstar=1
      end do  
      ! for each reducible g-vector find the irreducible one 
      ! and the symmetry operation 
      do ig=1,ngirred  ! for all irreducible gvecs
        do jg=1,ngred  ! for all reducible gvecs
          do isymm=1,nsymm  ! for all symmetry operations
            ! if the irreduc. gvec can be transformed into the reduc.
            ! one by the symmetry operation, the two belong to one star
            gtrans=matmul(symops(isymm)%mat,girred(ig)%kpt)
            if(absvec(gtrans-gred(jg)%kpt)/absvec(gtrans).lt.1E-3)
     &         then 
               gred(jg)%istar=ig
               gred(jg)%isymm=isymm
            end if   
          end do
        end do   
      end do
      !
      ! begin count members of each star
      do ig=1,ngirred
        do jg=1,ngred
          if(gred(jg)%istar==ig) then
               girred(ig)%nstar=girred(ig)%nstar+1    
          end if     
        end do ! jg
      end do ! ig
      ! end count members of each star
      ! 
      ! get phase for each red. g-vector
      do ig=1,ngred
        gred(ig)%phase=0.0D0
        istar0=gred(ig)%istar
        gred(ig)%phase=gred(ig)%phase+exp((0.0D0,1.0D0)
     &     *dot_product(matmul(symops(gred(ig)%isymm)%mat,
     &     girred(istar0)%kpt),symops(gred(ig)%isymm)%trala))
        gred(ig)%phase=gred(ig)%phase
     &     *dble(girred(istar0)%nstar)/dble(nsymm)
      end do
      write(12,fsubend)  
      return
     
      end subroutine  

!---------------------------------------------------------------------
      
      subroutine latsys(vecs0,atoms0,tol,lsys)
      use defs
      use misc, only : vecs2cellpars,getspecies,vecs2vol,rm_double_atoms
      use linalg, only : idet3,fdet3
      use transform, only : transformvectors
      use writecoords, only : write_coords,write_poscar
      implicit none
      double precision, intent(in):: vecs0(3,3),tol(2)
      double precision tolat(3)
      character(len=*), intent(inout) :: lsys
      type(atom), intent(in) :: atoms0(:)
      !local
      character(len=256) :: lsysp,lsysc,lsyst,lsys0
      double precision :: cellpars0(1:6),cellpars(1:6),vecs(3,3),volc
      double precision cellparsp(6),cellparst(6), cellparsc(1:6)
      double precision vecsc(3,3),vecst(3,3),vol0,volt,volp,vecsp(3,3)
      double precision :: bcc_angle=acos(-1.0d0/3.0d0)*180.0d0/Pi
      integer caxis,caxisp,caxisc,caxist,caxis0
      integer imat(3,3),imat0(3,3),i1,i2,i3,j1,j2,j3,k1,k2,k3
      integer imatt(3,3),imatc(3,3)
      integer :: imax=2
      integer ilat,ilat0,ilatp,ilatt,ilatc
      type(symop) :: symops(1)
      type(atom), allocatable :: atomsp(:),atomsc(:)
      type(atom), allocatable ::fewestatoms(:)
      type(element), allocatable :: species0(:),speciesp(:)
      type(element), allocatable :: speciesc(:)
      integer ispecies, nspecies,fewestspecies,iatom,natoms0
      logical lsym,lsym1,lsym2
      logical talk_old

      ! initializing
      print fsubstart, 'latsys'
      print '(8x)'
      print '(8x,75("#"))'
      print '(8x,"### Provided cell: ",56("#"))'
      print '(8x,75("#"))'
      print '(8x)'
      tolat=1.0E-3 ! tolerance in fractional coordinates to distinguish between atomic positions
      lsys='unknown'
      ilat0=0
      ilat=0
      call vecs2cellpars(vecs0,cellpars0)
      call vecs2vol(vecs0,vol0)
      volp=vol0
      vecsp=vecs0
      if(any(abs(cellpars0).le.1E-3)) then
            call error_stop('At least one cell parameter is about zero')
      end if
      print '(8x,"Provided cell vectors:",3(/,8x,3(F14.6)))',vecs0
      print '(8x,"Provided cell parameters:",/,8x,6(F14.6))',cellpars0
      print '(8x,"Provided volume: ",F10.4)',vol0
      print '(8x,"using tolerance for cell and angles: ",2(F10.6,1x))',   &
     &       tol(1:2)
      
      ! begin find symmetry of provided cell
      ilat0=0  
      lsys0='unknown'
      caxis0=0
      call icellsym(cellpars0,tol,ilat0,lsys0,caxis0)
      print '(8x,"Provided cell is ",A13," (Lattice system ",I1,")")',     &
     &       lsys0,ilat0
      ! end find symmetry of provided cell
      !
      ! begin find primitive cell
      print '(8x)'
      print '(8x,75("#"))'
      print '(8x,"### Primitive cell: ",55("#"))'
      print '(8x,75("#"))'
      print '(8x)'
      call getspecies(atoms0,species0)
      nspecies=size(species0)
      natoms0=size(atoms0)
      fewestspecies=1
      print '(8x,I3," species found")',nspecies
      do ispecies=1,nspecies
        if (species0(ispecies)%howmany.lt.fewestspecies) then
          fewestspecies=ispecies
        end if
      end do ! ispecies
      print '(8x,"species with fewest atoms: No",I3)', fewestspecies
      allocate(fewestatoms(species0(fewestspecies)%howmany))
      iatom=1
      print '(8x,"atoms of species with fewest atoms:")'
      do i1=1,natoms0
        if (species0(fewestspecies)%name==atoms0(i1)%name) then
          fewestatoms(iatom)=atoms0(i1)
          print '(8x,A4,x,3(F8.4))', fewestatoms(iatom)%name,             &
     &           fewestatoms(iatom)%where
          iatom=iatom+1
        end if
      end do
      symops(1)%mat=0.0d0
      symops(1)%mat(1,1)=1.0d0
      symops(1)%mat(2,2)=1.0d0
      symops(1)%mat(3,3)=1.0d0
      talk_old=talk
      talk=.false.
      do i1=-1,1,2
        do j1=2,size(fewestatoms)  
          symops(1)%trala=dble(i1)*(fewestatoms(j1)%where                 &
     &                              -fewestatoms(1)%where)
!          print '(8x,3(3(F8.4),/8x),3(F8.4))',symops(1)
          call checksymmlat1(symops,atoms0,vecs0,lsym)  
!          print '(8x,"invariant under translation by ",3(F8.4,x),L1)',    &
!     &       symops(1)%trala,lsym
          if(lsym) then
            do k1=1,3
              vecst=vecsp
              vecst(k1,1:3)=symops(1)%trala(1)*vecs0(1,:)                 &
     &                     +symops(1)%trala(2)*vecs0(2,:)                 &
     &                     +symops(1)%trala(3)*vecs0(3,:)              
              call vecs2vol(vecst,volt)
              call vecs2cellpars(vecst,cellparst)
              !print '(8x,"test cell:",/8x,3(3(F10.4,x),/8x))', vecst
              !print '(8x,"test volume: ",F10.4)',volt
              !print '(8x,"det (test cell): ",F10.4)',fdet3(vecst)
              if(fdet3(vecst).gt.tol(1)) then
                if(volt.lt.volp-tol(1).or.(volt.le.volp+tol(1).and.       &
     &             sum(cellparst(1:3)).lt.sum(cellparsp(1:3))+tol(1)))    &
     &          then
                   !print '(8x,"changing cell...")'
                   volp=volt
                   vecsp=vecst
                   cellparsp=cellparst
                   !print '(8x,"prim. volume: ",F10.4)',volp
!                   print '(8x,3(3(F10.4,x),/8x))',                          &
!     &                    vecsp(1,:),vecsp(2,:),vecsp(3,:)
                 end if
              end if ! det(vecst.gt.tol)
            end do ! k1
          end if ! lsym
        end do ! j1
      end do ! i1
      talk=talk_old
      print '(8x,"Primitive vectors:")'
      print '(8x,3(3(F10.4,x),/8x))',                                     &
     &                  vecsp(1,:),vecsp(2,:),vecsp(3,:)
      call vecs2cellpars(vecsp,cellparsp)
      print '(8x,"Primitive cell parameters:",/,8x,6(F14.6))',cellparsp
      print '(8x,"prim. volume: ",F10.4)',volp
      ! begin get atoms in primitive cell 
      call transformvectors(vecs0,vecsp,atoms0,atomsp,(/1,1,1/))
      !
      ! begin remove doubled atoms
      !
      call rm_double_atoms(atomsp,tolat)
      !
      ! end remove doubled atoms
      !
      print '(8x,"atoms in primitive cell:")'
      do iatom=1,size(atomsp)
        print '(8x,A4,x,3(F10.6,x))', atomsp(iatom)%name,                 &
     &         atomsp(iatom)%where
      end do
      ! end get atoms in primitive cell 
      if (volp.lt.vol0-tol(1)) then
         print fcomm,"You may want to switch to the primitive cell which  &
     & is smaller (see file PRIMCELL.vasp)"    
        ncomm=ncomm+1
        call getspecies(atomsp,speciesp)
        call write_coords("PRIMCELL.vasp","poscar",atomsp,size(atomsp),   &
     &                 speciesp,size(speciesp),vecsp)
      end if
      ! end find primitive cell
      !
      ! begin find symmetry of primitive cell
      ilatp=0  
      lsysp='unknown'
      call icellsym(cellparsp,tol,ilatp,lsysp,caxisp)
      print '(8x,"Primitive cell is ",A13," (Lattice system ",I1,")")',   &
     &       lsysp,ilatp
      ! end find symmetry of primitive cell
      !
      ! begin find conventional cell
      print '(8x)'
      print '(8x,75("#"))'
      print '(8x,"### Conventional cell: ",52("#"))'
      print '(8x,75("#"))'
      print '(8x)'
      ilatc=0
      lsysc="unknown"  
      caxisc=0
      imatc(1,1:3)=(/1,0,0/)
      imatc(2,1:3)=(/0,1,0/)
      imatc(3,1:3)=(/0,0,1/)
      vecsc=vecsp
      volc=volp
      print '(8x,"Looking for conventional cell...")'
      talk_old=talk
      talk=.false.
      do i1=-imax,imax
      do i2=-imax,imax
      do i3=-imax,imax
      do j1=-imax,imax
      do j2=-imax,imax
      do j3=-imax,imax
      do k1=-imax,imax
      do k2=-imax,imax
      do k3=-imax,imax
        imat0(1,1:3)=(/i1,i2,i3/)
        imat0(2,1:3)=(/j1,j2,j3/)
        imat0(3,1:3)=(/k1,k2,k3/)
        !print*,imat0
        !if (idet3(imat0).ne.0) then
        if (idet3(imat0).gt.0) then
          !vecs=matmul(imat0,vecsp)
          vecst=matmul(imat0,vecsp)
          !if (fdet3(vecs).gt.0.0d0.and.vecs(3,3).gt.0.0d0) then
          !if (fdet3(vecst).gt.0.0d0.and.vecst(3,3).gt.0.0d0) then
          !if (fdet3(vecst).gt.tol) then
          if (fdet3(vecst).gt.tol(1).and.vecst(1,1).gt.tol(1).and.        &
     &         vecst(2,2).gt.tol(1).and.vecst(3,3).gt.tol(1)) then
            call vecs2vol(vecst,volt)
            call vecs2cellpars(vecst,cellparst)
            call icellsym(cellparst,tol,ilatt,lsyst,caxist)
            if (ilatt>ilatc.or.(ilatt==ilatc.and.volt<volc-tol(1))) then
              ilatc=ilatt
              imatc=imat0
              vecsc=vecst
              caxisc=caxist
              lsysc=lsyst
              volc=volt
            end if  
          end if   
        end if
      end do !i1
      end do !i2
      end do !i3
      end do !j1
      end do !j2
      end do !j3
      end do !k1
      end do !k2
      end do !k3
      talk=talk_old
      ! begin get atoms in conventional cell 
      call transformvectors(vecsp,vecsc,atomsp,atomsc,(/1,1,1/))
      print '(8x,"atoms in conventional cell:")'
      do iatom=1,size(atomsc)
        print '(8x,A4,x,3(F10.6,x))', atomsc(iatom)%name,                 &
     &         atomsc(iatom)%where
      end do
      !
      ! begin remove doubled atoms
      !
      call rm_double_atoms(atomsc,tolat)
      !
      ! end remove doubled atoms
      !
      call getspecies(atomsc,speciesc)
      call write_coords("CONVCELL.vasp","poscar",atomsc,size(atomsc),   &
     &                 speciesc,size(speciesc),vecsc)
      ! end get atoms in conventional cell 
      !
      ! begin check if conventional cell is primitive, face-centered, ...
      !
      write(lsysc(len_trim(lsysc)+1:len_trim(lsysc)+2),'(" P")') 
      !
      ! monoclinic
      !
      if (ilatc==2) then
        !
        ! base-centered (C)
        !
        symops(1)%trala=(/0.5d0,0.5d0,0.0d0/)
        call checksymmlat1(symops,atomsc,vecsc,lsym)  
        if (lsym) then 
          write(lsysc(len_trim(lsysc)-1:len_trim(lsysc)+1),'(" C3")') 
        end if
        symops(1)%trala=(/0.5d0,0.0d0,0.5d0/)
        call checksymmlat1(symops,atomsc,vecsc,lsym)  
        if (lsym) then 
          write(lsysc(len_trim(lsysc)-1:len_trim(lsysc)+1),'(" C2")') 
        end if
        symops(1)%trala=(/0.0d0,0.5d0,0.5d0/)
        call checksymmlat1(symops,atomsc,vecsc,lsym)  
        if (lsym) then 
          write(lsysc(len_trim(lsysc)-1:len_trim(lsysc)+1),'(" C1")') 
        end if
      end if
      !
      ! orthorhombic
      !
      if (ilatc==3) then
        !
        ! base-centered (C)
        !
        symops(1)%trala=(/0.5d0,0.5d0,0.0d0/)
        call checksymmlat1(symops,atomsc,vecsc,lsym)  
        if (lsym) then 
          write(lsysc(len_trim(lsysc)-1:len_trim(lsysc)+1),'(" C3")') 
        end if
        symops(1)%trala=(/0.5d0,0.0d0,0.5d0/)
        call checksymmlat1(symops,atomsc,vecsc,lsym)  
        if (lsym) then 
          write(lsysc(len_trim(lsysc)-1:len_trim(lsysc)+1),'(" C2")') 
        end if
        symops(1)%trala=(/0.0d0,0.5d0,0.5d0/)
        call checksymmlat1(symops,atomsc,vecsc,lsym)  
        if (lsym) then 
          write(lsysc(len_trim(lsysc)-1:len_trim(lsysc)+1),'(" C1")') 
        end if
        !
        ! body-centered (I)
        !
        symops(1)%trala=(/0.5d0,0.5d0,0.5d0/)
        call checksymmlat1(symops,atomsc,vecsc,lsym)  
        if (lsym) then 
          write(lsysc(len_trim(lsysc)-1:len_trim(lsysc)),'(" I")') 
        end if
        !
        ! face-centered (F)
        !
        symops(1)%trala=(/0.5d0,0.5d0,0.0d0/)
        call checksymmlat1(symops,atomsc,vecsc,lsym)  
        if (lsym) then 
          symops(1)%trala=(/0.5d0,0.0d0,0.5d0/)
          call checksymmlat1(symops,atomsc,vecsc,lsym1)  
          if (lsym1) then 
            symops(1)%trala=(/0.0d0,0.5d0,0.5d0/)
            call checksymmlat1(symops,atomsc,vecsc,lsym2)  
            if (lsym2) then 
              write(lsysc(len_trim(lsysc)-1:len_trim(lsysc)),'(" F")')
            end if
          end if
        end if
        !
      end if ! ilatc==3
      !
      ! tetragonal
      !
      if (ilatc==4) then
        !
        ! body-centered (I)
        !
        symops(1)%trala=(/0.5d0,0.5d0,0.5d0/)
        call checksymmlat1(symops,atomsc,vecsc,lsym)  
        if (lsym) then 
          write(lsysc(len_trim(lsysc)-1:len_trim(lsysc)),'(" I")') 
        end if
        !
      end if ! ilatc==4
      !
      ! cubic
      !
      if (ilatc==7) then
        !
        ! body-centered (I)
        !
        symops(1)%trala=(/0.5d0,0.5d0,0.5d0/)
        call checksymmlat1(symops,atomsc,vecsc,lsym)  
        if (lsym) then 
          write(lsysc(len_trim(lsysc)-1:len_trim(lsysc)),'(" I")') 
        end if
        !
        ! face-centered (F)
        !
        symops(1)%trala=(/0.5d0,0.5d0,0.0d0/)
        call checksymmlat1(symops,atomsc,vecsc,lsym)  
        if (lsym) then 
          symops(1)%trala=(/0.5d0,0.0d0,0.5d0/)
          call checksymmlat1(symops,atomsc,vecsc,lsym1)  
          if (lsym1) then 
            symops(1)%trala=(/0.0d0,0.5d0,0.5d0/)
            call checksymmlat1(symops,atomsc,vecsc,lsym2)  
            if (lsym2) then 
              write(lsysc(len_trim(lsysc)-1:len_trim(lsysc)),'(" F")')
            end if
          end if
        end if
        !
      end if ! ilatc==7
      !
      ! end check if conventional cell is primitive, face-centered, ...
      !
      ! end find conventional cell
      !
      print '(8x,"Conventional lattice system is ",A13," (",I1,")")',     &
     &     lsysc,ilatc
      print '(8x,"Conventional cell vectors:",3(/,8x,3(F14.6)))',vecsc
      call vecs2cellpars(vecsc,cellparsc)
      call vecs2vol(vecsc,volc)
      print '(8x,"Conventional cell parameters:",/,8x,6(F14.6))',         &
     &       cellparsc
      print '(8x,"Conventional cell volume: ",F14.6)',volc
      print '(8x,"Transformation matrix:",3(/,8x,3(I4)))',imatc

      print fsubendext, 'latsys' 
      return
     
      end subroutine  

c---------------------------------------------------------------------
      
      subroutine icellsym(cellpars,tol,ilat,lsys,caxis)
      implicit none
      double precision, intent(in) :: cellpars(1:6),tol(2)
      character(len=*), intent(inout) :: lsys
      integer, intent(inout) :: ilat,caxis
      !
      ilat=0  
      if (lftriclinic(cellpars,tol)) then
        ilat=1      
      end if
      !
      if (lfmonoclinic(cellpars,tol,caxis)) then
        ilat=1      
      end if
      !
      if (lforthorhombic(cellpars,tol)) then
        ilat=3      
      end if
      !
      if (lftetragonal(cellpars,tol,caxis)) then
        ilat=4     
      end if
      !
      if (lfrhombohedral(cellpars,tol)) then
        ilat=5     
      end if
      !
      if (lfhexagonal(cellpars,tol,caxis)) then
        ilat=6    
      end if
      !
      if (lfcubic(cellpars,tol)) then
        ilat=7 
      end if   
      !
      !
      if (ilat==1) then
        write(lsys,'(A9)') 'triclinic'
      end if
      if (ilat==2) then
        write(lsys,'(A10,I1)') 'monoclinic',caxis
      end if
      if (ilat==3) then
        write(lsys,'(A12)') 'orthorhombic'
      end if
      if (ilat==4) then
        write(lsys,'(A10,I1)') 'tetragonal',caxis
      end if
      if (ilat==5) then
        write(lsys,'(A12)') 'rhombohedral'
      end if
      if (ilat==6) then
        write(lsys,'(A9,I1)') 'hexagonal',caxis
      end if
      if (ilat==7) then
        write(lsys,'(A5)') 'cubic'
      end if
      !
      end subroutine icellsym

c---------------------------------------------------------------------      

      logical function lftriclinic(cellpars,tol)
      implicit none
      double precision cellpars(6), tol(2)

      lftriclinic=.false.
      if ((abs(cellpars(4)-90.0d0).gt.tol(2).and.                         &
     &   abs(cellpars(5)-90.0d0).gt.tol(2)).or.                           &
     &    (abs(cellpars(4)-90.0d0).gt.tol(2).and.                         &
     &   abs(cellpars(6)-90.0d0).gt.tol(2)).or.                           &
     &    (abs(cellpars(5)-90.0d0).gt.tol(2).and.                         &
     &   abs(cellpars(6)-90.0d0).gt.tol(2))) then
         lftriclinic=.true.
      end if

      end function lftriclinic

c---------------------------------------------------------------------      

      logical function lfmonoclinic(cellpars,tol,caxis)
      implicit none
      double precision cellpars(6), tol(2)
      integer caxis

      lfmonoclinic=.false.
      caxis=0
      if (abs(cellpars(4)-90.0d0).le.tol(2).and.                         &
     &   abs(cellpars(5)-90.0d0).le.tol(2).and.                          &
     &   abs(cellpars(6)-90.0d0).gt.tol(2)) then
         lfmonoclinic=.true.
         caxis=3
      end if
      if (abs(cellpars(4)-90.0d0).le.tol(2).and.                         &
     &   abs(cellpars(5)-90.0d0).gt.tol(2).and.                          &
     &   abs(cellpars(6)-90.0d0).le.tol(2)) then                         &
         lfmonoclinic=.true.
         caxis=2
      end if
      if (abs(cellpars(4)-90.0d0).gt.tol(2).and.                         &
     &   abs(cellpars(5)-90.0d0).le.tol(2).and.                          &
     &   abs(cellpars(6)-90.0d0).le.tol(2)) then                        
         lfmonoclinic=.true.
         caxis=1
      end if
      !
      end function lfmonoclinic

c---------------------------------------------------------------------      

      logical function lforthorhombic(cellpars,tol)
      implicit none
      double precision cellpars(6), tol(2)

      lforthorhombic=.false.
      if (abs(cellpars(1)-cellpars(2)).gt.tol(1).and.                     &
     &   abs(cellpars(1)-cellpars(3)).gt.tol(1).and.                      &
     &   abs(cellpars(2)-cellpars(3)).gt.tol(1).and.                      &
     &   abs(cellpars(4)-90.0d0).le.tol(2).and.                           &
     &   abs(cellpars(5)-90.0d0).le.tol(2).and.                           &
     &   abs(cellpars(6)-90.0d0).le.tol(2)) then
              lforthorhombic=.true.
      end if
      !
      end function lforthorhombic

c---------------------------------------------------------------------      

      logical function lftetragonal(cellpars,tol,caxis)
      implicit none
      double precision cellpars(6), tol(2)
      integer caxis

      lftetragonal=.false.
      caxis=0
      if (abs(cellpars(4)-90.0d0).le.tol(2).and.                          &
     &   abs(cellpars(5)-90.0d0).le.tol(2).and.                           &
     &   abs(cellpars(6)-90.0d0).le.tol(2)) then
         if (abs(cellpars(1)-cellpars(2)).le.tol(1).and.                  &
     &      abs(cellpars(2)-cellpars(3)).gt.tol(1)) then    
            lftetragonal=.true.
            caxis=3         
         end if
         if (abs(cellpars(1)-cellpars(2)).gt.tol(1).and.                  &
     &      abs(cellpars(1)-cellpars(3)).le.tol(1)) then               
            lftetragonal=.true.
            caxis=2         
         end if
         if (abs(cellpars(1)-cellpars(2)).gt.tol(1).and.                  &
     &      abs(cellpars(2)-cellpars(3)).le.tol(1)) then                
            lftetragonal=.true.
            caxis=1         
         end if
      end if
      
      end function lftetragonal
      !
c---------------------------------------------------------------------      

      logical function lfrhombohedral(cellpars,tol)
      implicit none
      double precision cellpars(6), tol(2)

      lfrhombohedral=.false.
      if (abs(cellpars(1)-cellpars(2)).le.tol(1).and.                     &
     &   abs(cellpars(1)-cellpars(3)).le.tol(1).and.                      &
     &   abs(cellpars(2)-cellpars(3)).le.tol(1).and.                      &
     &   abs(cellpars(4)-90.0d0).gt.tol(2).and.                           &
     &   abs(cellpars(4)-cellpars(5)).le.tol(2).and.                      &
     &   abs(cellpars(5)-cellpars(6)).le.tol(2).and.                      &
     &   abs(cellpars(4)-cellpars(6)).le.tol(2)) then
         lfrhombohedral=.true.
      end if

      end function lfrhombohedral

!---------------------------------------------------------------------

      logical function lfhexagonal(cellpars,tol,caxis)
      implicit none
      double precision cellpars(6),tol(2)
      integer caxis
      
      lfhexagonal=.false.
      caxis=0
      if ((abs(cellpars(1)-cellpars(2)).le.tol(1).and.                    &
     &   abs(cellpars(4)-90.0d0).le.tol(2).and.                           &
     &   abs(cellpars(5)-90.0d0).le.tol(2).and.                           &
     &   (abs(cellpars(6)-120.0d0).le.tol(2).or.                          &
     &   abs(cellpars(6)-60.0d0).le.tol(2)))) then
         lfhexagonal=.true.
         caxis=3
      end if
      if  (abs(cellpars(1)-cellpars(3)).le.tol(1).and.                    &
     &   abs(cellpars(4)-90.0d0).le.tol(2).and.                           &
     &   abs(cellpars(6)-90.0d0).le.tol(2).and.                           &
     &   (abs(cellpars(5)-120.0d0).le.tol(2).or.                          &
     &   abs(cellpars(5)-60.0d0).le.tol(2))) then 
         lfhexagonal=.true.
         caxis=2
      end if
      if (abs(cellpars(2)-cellpars(3)).le.tol(1).and.                     &
     &   abs(cellpars(5)-90.0d0).le.tol(2).and.                           &
     &   abs(cellpars(6)-90.0d0).le.tol(2).and.                           &
     &   (abs(cellpars(4)-120.0d0).le.tol(2).or.                          &
     &   abs(cellpars(4)-60.0d0).le.tol(2))) then       
         lfhexagonal=.true.
         caxis=1
      end if
 
      end function lfhexagonal

!---------------------------------------------------------------------

      logical function lfcubic(cellpars,tol)
      implicit none
      double precision cellpars(6),tol(2)
      
      lfcubic=.false.
      ! begin sc primitive
      if (abs(cellpars(1)-cellpars(2)).le.tol(1).and.                     &
     &   abs(cellpars(1)-cellpars(3)).le.tol(1).and.                      &
     &   abs(cellpars(2)-cellpars(3)).le.tol(1).and.                      &
     &   abs(cellpars(4)-90.0d0).le.tol(2).and.                           &
     &   abs(cellpars(5)-90.0d0).le.tol(2).and.                           &
     &   abs(cellpars(6)-90.0d0).le.tol(2)) then
         lfcubic=.true.
      end if
      ! end sc primitive
!      ! begin fcc primitive
!      if (abs(cellpars(1)-cellpars(2)).le.tol.and.                        &
!     &   abs(cellpars(1)-cellpars(3)).le.tol.and.                         &
!     &   abs(cellpars(2)-cellpars(3)).le.tol.and.                         &
!     &   abs(cellpars(4)-cellpars(5)).le.tol.and.                         &
!     &   abs(cellpars(5)-cellpars(6)).le.tol.and.                         &
!     &   abs(cellpars(4)-60.0d0).le.tol) then
!         lfcubic=.true.
!      end if
!      ! end fcc primitive
!      ! begin bcc primitive
!      if (abs(cellpars(1)-cellpars(2)).le.tol.and.                        &
!     &   abs(cellpars(1)-cellpars(3)).le.tol.and.                         &
!     &   abs(cellpars(2)-cellpars(3)).le.tol.and.                         &
!     &   abs(cellpars(4)-cellpars(5)).le.tol.and.                         &
!     &   abs(cellpars(5)-cellpars(6)).le.tol.and.                         &
!     &   abs(cellpars(4)-bcc_angle).le.tol) then
!         lfcubic=.true.
!      end if
!      ! end bcc primitive

      end function lfcubic

!---------------------------------------------------------------------

      end module
