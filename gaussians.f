        module gaussians

        contains

c---------------------------------------------------------------------

        subroutine basis4fiesta(ele,zetamin,zetamax,nzeta,lmaxg,
     &                          rmax,npoints)

        use defs
        use integr
        use misc
        implicit none
        ! input(/output) variables
        character (len=*) ele
        double precision zetamin,zetamax
        integer nzeta ! number of functions per angular momentum
        integer lmaxg ! maximum angular momentum for which gaussians are to be made
        double precision :: rmax ! extension of radial grid 
        integer npoints ! number of radial gridpoints
        ! internal variables
        character filename*8 ! output filename
        integer norbnl  ! number of gaussian orbitals
        integer l       ! total angular momentum (l)
        integer iorb    ! running index for orbital
        integer izeta   ! index of zeta
        integer ng      ! number of gaussians in the contraction 
        integer ir      ! running index for radial grid point
        double precision zeta ! gaussian exponent
        double precision cg ! gaussian coefficient (normalization factor)
        double precision norm ! norm of gaussian basis function 
        double precision, allocatable :: rgrid(:) ! radial grid for normalization of basis function
        double precision, allocatable :: rfield(:) ! values of basis function on radial grid 

        if(talk) print fsubstart,'basis4fiesta'
        
        npoints=20000
        rmax=40.0d0
        allocate(rgrid(npoints),rfield(npoints))

        ! open file "el2.ion", where el is the elem,ent name, and write preamble
        filename(1:len(filename))=' '
        FMT2=' '
        FMT1=' '
        WRITE(FMT2,*) trim(adjustl(ele))
        WRITE(FMT1,*) '(A',len_trim(adjustl(FMT2)),',"2.ion")'
        write(filename,FMT1) trim(adjustl(ele)) 
        if(talk) print '(8x,"Writing to file ",A8)',filename
        open(51,file=filename,status='unknown')
        write(51,'("<preamble>")')
        write(51,'("zmin zmax nzeta lmax")')
        write(51,'(2(F12.6,1x),I0,1x,I0)') zetamin,zetamax,nzeta,lmaxg
        write(51,'("</preamble>")')
        ! write lmax, number of orbitals 
        norbnl=(lmaxg+1)*nzeta
        write(51,'(I0,1x,I0)') lmaxg,norbnl
        !set up radial grid
        do ir=1,size(rgrid)
          rgrid(ir)=rmax*dble(ir-1)/dble(size(rgrid))
        end do ! ir
        ng=1 ! so far, only 1 gaussian per basis function
        do l=0,lmaxg  ! start with lowest l
          do izeta=1,nzeta ! for each l, make nzeta basis functions 
            ! for each orbital, write 
            ! 0 1 1 # l zeta ng
            ! exponential decay constant, prefactor
            write(51,'(3(I0,1x),"# l zeta ng")') l,izeta,ng
            zeta=zetamin*(zetamax/zetamin)**(dble(izeta-1)
     &         /dble(nzeta-1))
            ! set up gaussian function
            rfield(:)=rgrid(:)**(2*l+2)*exp(-2.0d0*zeta*rgrid(:)**2)
!            norm=int_1d_taylorn(rgrid,rfield,4)
            norm=int_1d_a(rgrid,rfield)  ! same integration method as in Fiesta in order to obtain same norm
            ! analytical norm (Bronstein):
!            norm=sqrt(Pi)*dble(dblfactorial(2*l+1)) 
!     &           / ( 2**(l+2) * (2.0d0*zeta)**(dble(l)+1.50d0) )
            cg=1.0d0/sqrt(norm) ! prefactor or coefficient    
            write(51,'(2(1p,E25.16))')  zeta,cg
          end do ! izeta
        end do ! l
        close(51)

        if(talk) print fsubendext,'basis4fiesta'

        end subroutine

c---------------------------------------------------------------------

        subroutine expo_range(filename)

        use defs
        use misc
        implicit none
        ! input(/output) variables
        character (len=*) filename ! input filename
        ! internal variables
        double precision zetamin,zetamax
        double precision zeta ! gaussian exponent
        character line*1024

        if(talk) print fsubstart,'expo_range'
        
        ! open file "filename"
        open(51,file=filename,status='old',err=200)
        zetamin=1.0E20
        zetamax=1.0E-20
 10     read(51,'(A1024)',end=100) line
        if (index(line,"E+").gt.0.or.index(line,"E-").gt.0) then
          read(line,*) zeta
          if (zeta.gt.zetamax) zetamax=zeta
          if (zeta.lt.zetamin) zetamin=zeta
        end if
        goto 10

 100    close(51)
        if(talk) print '("zetamin=",(1p,E25.16),F20.10)',zetamin,zetamin
        if(talk) print '("zetamax=",(1p,E25.16),F20.10)',zetamax,zetamax
        if(talk) print '("range=",(1p,E25.16),F20.10)',zetamax-zetamin,
     &     zetamax-zetamin
        open(51,file="EXPO_RANGE.DAT",status="replace")
!        write(51,'("zetamin=",(1p,E25.16),F20.10)')zetamin,zetamin
!        write(51,'("zetamax=",(1p,E25.16),F20.10)')zetamax,zetamax
!        write(51,'("range=",(1p,E25.16),F20.10)')zetamax-zetamin,
!     &     zetamax-zetamin
        write(51,'("zetamin=",F24.12)')zetamin
        write(51,'("zetamax=",F24.12)')zetamax
        write(51,'("range=",F24.12)')zetamax-zetamin 
        close(51)
        if(talk) print fsubendext,'expo_range'
        return

! Error handling
 200    nerr=nerr+1
        close(51)
        print ferrmssg, "Could not open the input file"
        

        end subroutine

C---------------------------------------------------------------------

      subroutine overlapwf(file_wf1,file_wf2,file_smat)
      ! calculates the overlap between two sets of wavefunctions wf1 and
      ! wf2, whose coefficient wrt a Gaussian basis are in the files
      ! file_wf1/2, and where the overlap matrix is in file_smat.
      ! Overlaps are written to overlaps.dat
      use defs
      implicit none
      character(len=*), intent(in) :: file_wf1,file_wf2,file_smat 
      ! internal variables
      integer nbf,nbf2,nocc ! number of basis functions and occupied states
      double precision, allocatable ::wf1(:,:),wf2(:,:),smat(:,:),ov(:)
      double precision, allocatable :: en1(:),en2(:)
      double precision maov,ovhomo,ovlumo
      integer istate,ibf,jbf,n_nondegen
      logical, allocatable :: degen(:)
      ! 
      if(talk) print fsubstart,"overlapwf"
      !
      ! read first WF file
      open(51,file=file_wf1,status='old',err=101)
      rewind(51)
      read(51,*,err=101,end=101) nbf,nocc
      allocate(wf1(nbf,nbf),wf2(nbf,nbf),smat(nbf,nbf),ov(nbf),en1(nbf),
     &         en2(nbf),degen(nbf))
      do istate=1,nbf
        read(51,*,end=101,err=101) en1(istate)
        do ibf=1,nbf
          read(51,*,end=101,err=101) wf1(istate,ibf)
        end do ! i
      end do ! istate
      close(51)
      !
      ! read 2. WF file
      open(51,file=file_wf2,status='old',err=201)
      read(51,*,err=201,end=201) nbf2,nocc
      if (nbf2.ne.nbf) goto 202
      do istate=1,nbf
        read(51,*,end=201,err=201) en2(istate)
        do ibf=1,nbf
          read(51,*,end=201,err=201) wf2(istate,ibf)
        end do ! i
      end do ! istate
      close(51)
      ! 
      ! read smat file
      open(51,file=file_smat,status='old',err=301)
      read(51,*,err=301,end=301) nbf2,nocc
      if (nbf2.ne.nbf) goto 302
      do ibf=1,nbf
        do jbf=1,nbf
          read(51,*,end=301,err=301) smat(ibf,jbf)
        end do ! jbf
      end do ! ibf
      close(51)
      !
      ! find out which states are degenerate (these screw up overlaps)
      degen=.false.
      do ibf=1,nbf
        do jbf=1,nbf
          if ((abs(en1(ibf)-en1(jbf)).lt.1.0E-3
     &   .or.abs(en2(ibf)-en2(jbf)).lt.1.0E-3).and.ibf.ne.jbf) then
             degen(ibf)=.true.
          end if
        end do ! jbf
      end do !ibf
      !
      ! calculate overlap for each state and mean absolute overlap
      ov=0.0d0
      maov=0.0d0
      do istate=1,nbf
        do ibf=1,nbf
          do jbf=1,nbf
            ov(istate)=ov(istate)+wf1(istate,ibf)*wf2(istate,jbf)
     &                 *smat(ibf,jbf)
          end do ! jbf
        end do ! ibf
        maov=maov+abs(ov(istate))
      end do ! istate
      maov=maov/dble(nbf)
      ! print mean absolute overlap to standardout
      print'(8x,"mean abs. overlap: ",F10.6)',maov
      !
      ! calculate mean absolute overlap of valence states
      maov=0.0d0
      do istate=1,nocc
        maov=maov+abs(ov(istate))
      end do ! istate
      maov=maov/dble(nocc)
      ! print mean absolute valence overlap to standardout
      print'(8x,"mean abs. valence overlap: ",F10.6)',maov
      ! 
      ! calculate mean absolute overlap of conduction states
      maov=0.0d0
      do istate=nocc+1,nbf
        maov=maov+abs(ov(istate))
      end do ! istate
      maov=maov/dble(nbf-nocc)
      ! print mean absolute conduction overlap to standardout
      print'(8x,"mean abs. conduction overlap: ",F10.6)',maov
      ! 
      ! calculate mean absolute overlap of nondegenerate states
      maov=0.0d0
      n_nondegen=0
      do istate=1,nbf
        if(.not.degen(istate)) then
          n_nondegen=n_nondegen+1
          maov=maov+abs(ov(istate))
        end if
      end do ! istate
      if(n_nondegen.gt.0) maov=maov/dble(n_nondegen)
      ! print mean absolute nondegenerate overlap to standardout
      print'(8x,"mean abs. nondegenerate overlap: ",F10.6)',maov
      !
      ! calculate mean absolute overlap of nondegenerate valence states
      maov=0.0d0
      n_nondegen=0
      do istate=1,nocc
        if(.not.degen(istate)) then
          n_nondegen=n_nondegen+1
          maov=maov+abs(ov(istate))
        end if
      end do ! istate
      if(n_nondegen.gt.0) maov=maov/dble(n_nondegen)
      ! print mean absolute nondegenerate valence overlap to standardout
      print'(8x,"mean abs. nondegenerate valence overlap: ",F10.6)',maov
      !
      ! calculate mean absolute overlap of nondegenerate conduction states
      maov=0.0d0
      n_nondegen=0
      do istate=nocc+1,nbf
        if(.not.degen(istate)) then
          n_nondegen=n_nondegen+1
          maov=maov+abs(ov(istate))
        end if
      end do ! istate
      if(n_nondegen.gt.0) maov=maov/dble(n_nondegen)
      ! print mean absolute nondegenerate conduction overlap to standardout
      print'(8x,"mean abs. nondegenerate conduction overlap: ",F10.6)',
     &       maov
      !
      ! print overlaps of homo and lumo to standardout
      print'(8x,"abs. homo overlap: ",F10.6)',abs(ov(nocc))
      print'(8x,"abs. lumo overlap: ",F10.6)',abs(ov(nocc+1))
      ! 
      ! write overlap to file 
      open(51,file="OVERLAPS.DAT",status="replace")
      do istate=1,nbf
        write(51,'(I8,F16.12)') istate,ov(istate)
      end do
      close(51)
      ! 
      ! end normally
      deallocate(smat,wf1,wf2,ov)
      if(talk) print fsubendext,"overlapwf"
      return
      !
      ! errors 
 101  nerr=nerr+1
      print ferrmssg,"s.t. wrong with first WF file"
      close(51)
      return
 201  nerr=nerr+1
      print ferrmssg,"s.t. wrong with 2. WF file"
      close(51)
      return
 202  nerr=nerr+1
      print ferrmssg,"WF1 and WF2 differ in the number of basis fcts."
      close(51)
      return
 301  nerr=nerr+1
      print ferrmssg,"s.t. wrong with smat file"
      close(51)
      return
 302  nerr=nerr+1
      print ferrmssg,"Smat and WFs differ in the number of basis fcts."
      close(51)
      return
      !
      end subroutine overlapwf


c---------------------------------------------------------------------

        end module
