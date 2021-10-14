      module integr
      ! contains integration routines
      implicit none

      contains

c---------------------------------------------------------------------

        subroutine fourier_1d(file1,ftmode,qmax_chosen,nq_chosen)
        use defs
        implicit none
        ! reads data grid and field
        !
        character(len=*) file1,ftmode
        ! internal variables
        double precision, allocatable :: rgrid(:),rfield(:),qgrid(:)
        complex, allocatable :: qfield(:)
        character(len=256) line
        integer nr,nq ! number of data points
        integer ir,iq 
        !double precision h ! stepsize 
        character(len=256) file1a
        double precision, optional :: qmax_chosen
        integer, optional :: nq_chosen
        !
        if(talk) print fsubstart,"fourier_1d" 
        if (talk) print '(8x,"file: ",A40)',file1
        if (talk) print '(8x,"FT mode: ",A40)',ftmode
        !
        ! begin read grid and field from file
        open(51,file=file1,status="old",err=101)
        nr=0
        rewind(51)
 10     read(51,'(A256)',end=11,err=101) line
        line=adjustl(line)
        if(len_trim(line).gt.0.and.line(1:1).ne."#") nr=nr+1
        goto 10
 11     rewind(51)
        allocate(rgrid(nr),rfield(nr))
        ir=1
 12     read(51,'(A256)',end=13,err=101) line     
        line=adjustl(line)
        if(len_trim(line).gt.0.and.line(1:1).ne."#") then
          read(line(1:256),*) rgrid(ir),rfield(ir)
          ir=ir+1
        end if
        goto 12
 13     close(51)
        ! end read grid and field from file
        !
         call fourier_1d_simple(rgrid,rfield,qgrid,qfield,qmax_chosen,    &
     &        nq_chosen)
        !
        ! begin write FT to disk 
        nq=size(qfield)
        write(file1a,'(a,a)') adjustl(trim(file1)),".FT" 
        open(51,file=trim(file1a),status="replace")
        do iq=1,nq
          write(51,'(F20.6,2F15.6)') qgrid(iq),qfield(iq)       
        end do 
        close(51)
        ! end write deriv. to disk 
        !
        ! end normally:
        if(talk) print fsubendext, "fourtier_1d"
        return
        !
        ! Errors:
 101    nerr=nerr+1
        print ferrmssg," when opening/reading file"
        close(51)
        return

        end subroutine fourier_1d

c---------------------------------------------------------------------

        subroutine fourier_1d_reverse(file1,ftmode)
        use defs
        implicit none
        ! reads data grid and field in Fourier space
        !
        character(len=*) file1,ftmode
        ! internal variables
        double precision, allocatable :: qgrid(:),rgrid(:)
        double precision z1,z2
        complex, allocatable :: qfield(:),rfield(:)
        character(len=256) line
        integer nr,nq ! number of data points
        integer ir,iq 
        character(len=256) file1a
        !
        if(talk) print fsubstart,"fourier_1d_reverse" 
        if (talk) print '(8x,"file: ",A40)',file1
        if (talk) print '(8x,"FT mode: ",A40)',ftmode
        !
        ! begin read q grid and field from file
        open(51,file=file1,status="old",err=101)
        nq=0
        rewind(51)
 10     read(51,'(A256)',end=11,err=101) line
        line=adjustl(line)
        if(len_trim(line).gt.0.and.line(1:1).ne."#") nq=nq+1
        goto 10
 11     rewind(51)
        allocate(qgrid(nq),qfield(nq))
        iq=1
 12     read(51,'(A256)',end=13,err=101) line     
        line=adjustl(line)
        if(len_trim(line).gt.0.and.line(1:1).ne."#") then
          read(line(1:256),*) qgrid(iq),z1,z2
          qfield(iq)=z1+icomp*z2
          iq=iq+1
        end if
        goto 12
 13     close(51)
        ! end read q grid and field from file
        !
        call fourier_1d_reverse_simple(qgrid,qfield,rgrid,rfield)
        !
        ! begin write FTR to disk 
        nr=size(rfield)
        write(file1a,'(a,a)') adjustl(trim(file1)),".FTR" 
        open(51,file=trim(file1a),status="replace")
        do ir=1,nr
          write(51,'(F20.6,2F15.6)') rgrid(ir),rfield(ir)       
        end do 
        close(51)
        ! end write FTR to disk 
        !
        ! end normally:
        if(talk) print fsubendext, "fourtier_1d_reverse"
        return
        !
        ! Errors:
 101    nerr=nerr+1
        print ferrmssg," when opening/reading file"
        close(51)
        return

        end subroutine fourier_1d_reverse

c---------------------------------------------------------------------

        double precision function int_1d(rgrid,rfield)
        ! performs a simple 1 d integration of a 1d function "rfield" on a 1d
        ! grid "rgrid"  
        use defs
        !use misc
        implicit none
        double precision, intent(in) :: rgrid(:),rfield(:)
        double precision dr
        integer nr,ir

        int_1d=0.0D0
        nr=size(rgrid)
        do ir=1,nr-1
          dr=rgrid(ir+1)-rgrid(ir)
          int_1d=int_1d+rfield(ir)*dr
        end do
        return
        end function

c---------------------------------------------------------------------

        double precision function int_1d_a(rgrid,rfield)
        ! performs a simple 1 d integration of a 1d function "rfield" on a 1d
        ! grid "rgrid"  
        use defs
        !use misc
        implicit none
        double precision, intent(in) :: rgrid(:),rfield(:)
        double precision dr
        integer nr,ir

        int_1d_a=0.0D0
        nr=size(rgrid)
        do ir=1,nr-1
          dr=rgrid(ir+1)-rgrid(ir)
          int_1d_a=int_1d_a+rfield(ir)*dr
        end do
        do ir=nr,nr
          dr=rgrid(ir)-rgrid(ir-1)
          int_1d_a=int_1d_a+rfield(ir)*dr
        end do
        return
        end function

c---------------------------------------------------------------------

        double precision function int_1d_taylorn0(rgrid,rfield,ntaylor)
        ! performs a 1 d integration of a 1d function "rfield" on a 1d
        ! grid "rgrid" taylor-expanding the rfield to ntaylorth order 
        use defs
        !use misc
        use combi
        implicit none
        double precision, intent(in) :: rgrid(:),rfield(:)
        integer, intent(in) :: ntaylor
        double precision dr
        double precision, allocatable :: fprime(:,:)
        integer nr,ir,itaylor

        int_1d_taylorn0=0.0D0
        nr=size(rgrid)
        ! calculate derivatives
        allocate(fprime(0:ntaylor,nr))
        fprime=0.0D0
        fprime(0,1:nr)=rfield(1:nr)   ! 0th order = rfield
        do itaylor=1,ntaylor
          ! nth derivative
          do ir=2,nr-1
            fprime(itaylor,ir)=(fprime(itaylor-1,ir+1)
     &        -fprime(itaylor-1,ir-1))
     &        /(rgrid(ir+1)-rgrid(ir-1))
          end do
          ! set derivative constant at the integral borders
          fprime(itaylor,1)=fprime(itaylor,2)
          fprime(itaylor,nr)=fprime(itaylor,nr-1)
        end do  ! itaylor
        ! start integration
        do ir=1,nr-1
          dr=rgrid(ir+1)-rgrid(ir)
          do itaylor=0,ntaylor
            int_1d_taylorn0=int_1d_taylorn0 
     &      +fprime(itaylor,ir)*dr**(itaylor+1)/dble(facult(itaylor+1))! nth order Taylor
          end do
        end do
        return
        end function ! int_1d_taylorn0

c---------------------------------------------------------------------

        double precision function int_1d_taylorn(rgrid,rfield,ntaylor)
        ! performs a 1 d integration of a 1d function "rfield" on a 1d
        ! grid "rgrid" taylor-expanding the rfield to ntaylorth order 
        use defs
        !use misc
        use combi
        implicit none
        double precision, intent(in) :: rgrid(:),rfield(:)
        integer, intent(in) :: ntaylor
        !double precision dr
        double precision, allocatable :: drgrid(:)
        double precision, allocatable :: fprime(:,:)
        integer nr,ir,itaylor

        int_1d_taylorn=0.0D0
        nr=size(rgrid)
        ! calculate derivatives
        allocate(fprime(0:ntaylor,nr))
        allocate(drgrid(1:nr))
        fprime=0.0D0
        drgrid=0.0D0
        int_1d_taylorn=0.0d0
        fprime(0,1:nr)=rfield(1:nr)   ! 0th order = rfield
        do itaylor=1,ntaylor
          ! nth derivative
!          do ir=2,nr-1
!            fprime(itaylor,ir)=(fprime(itaylor-1,ir+1)
!     &        -fprime(itaylor-1,ir-1))
!     &        /(rgrid(ir+1)-rgrid(ir-1))
          fprime(itaylor,2:nr-1)=(fprime(itaylor-1,3:nr)
     &      -fprime(itaylor-1,1:nr-2))
     &      /(rgrid(3:nr)-rgrid(1:nr-2))
!          end do
          ! set derivative constant at the integral borders
          fprime(itaylor,1)=fprime(itaylor,2)
          fprime(itaylor,nr)=fprime(itaylor,nr-1)
        end do  ! itaylor
        ! start integration
!        do ir=1,nr-1
!          dr=rgrid(ir+1)-rgrid(ir)
          drgrid(1:nr-1)=rgrid(2:nr)-rgrid(1:nr-1)
          do itaylor=0,ntaylor
!            int_1d_taylorn=int_1d_taylorn 
!     &      +fprime(itaylor,ir)*dr**(itaylor+1)/dble(facult(itaylor+1))! nth order Taylor
            int_1d_taylorn=int_1d_taylorn 
     &      +dot_product(fprime(itaylor,1:nr-1),
     &                   drgrid(1:nr-1)**(itaylor+1))
     &      /dble(facult(itaylor+1))! nth order Taylor
          end do
!        end do
        return
        end function

c---------------------------------------------------------------------

        double precision function int_radial(rgrid,rfield)
        ! performs a simple 1 d integration of a radial-dependent function "rfield" on a radial
        ! grid "rgrid". Differs from int_1d by the factor r^2 in the
        ! integral 
        use defs
        !use misc
        implicit none
        double precision, intent(in) :: rgrid(:),rfield(:)
        double precision dr
        integer nr,ir

        int_radial=0.0D0
        nr=size(rgrid)
        ! integrate from 0 to first radial gridpoint assuming that
        ! integrand rfield is constant here
        ir=1
        int_radial=int_radial+rfield(ir)*rgrid(ir)**3/3.0D0
        do ir=1,nr-1
          dr=rgrid(ir+1)-rgrid(ir)
          int_radial=int_radial+rfield(ir)*rgrid(ir)**2*dr
        end do
        return
        end function

c---------------------------------------------------------------------

        double precision function int_radial_taylorn(rgrid,rfield,
     &    ntaylor)
        ! performs a 1 d integration of a 1d function "rfield" on a 1d
        ! grid "rgrid" taylor-expanding the rfield to ntaylorth order 
        use defs
        !use misc
        use combi
        implicit none
        double precision, intent(in) :: rgrid(:),rfield(:)
        integer, intent(in) :: ntaylor
        double precision dr
        double precision, allocatable :: fprime(:,:)
        integer nr,ir,itaylor

        int_radial_taylorn=0.0D0
        nr=size(rgrid)
        ! calculate derivatives
        allocate(fprime(0:ntaylor,nr))
        fprime=0.0D0
        fprime(0,1:nr)=rfield(1:nr)   ! 0th order = rfield
        do itaylor=1,ntaylor
          ! nth derivative
          do ir=2,nr-1
            fprime(itaylor,ir)=(fprime(itaylor-1,ir+1)
     &        -fprime(itaylor-1,ir-1))
     &        /(rgrid(ir+1)-rgrid(ir-1))
          end do
          ! set derivative constant at the integral borders
          fprime(itaylor,1)=fprime(itaylor,2)
          fprime(itaylor,nr)=fprime(itaylor,nr-1)
        end do  ! itaylor
        ! start integration
        ! integrate from 0 to first radial gridpoint assuming that
        ! integrand rfield is constant here
        ir=1
        int_radial_taylorn=int_radial_taylorn
     &    +rfield(ir)*rgrid(ir)**3/3.0D0
        do ir=1,nr-1
          dr=rgrid(ir+1)-rgrid(ir)
          do itaylor=0,ntaylor
            int_radial_taylorn=int_radial_taylorn 
     &      +fprime(itaylor,ir)*dr**(itaylor+1)/dble(facult(itaylor+1))
     &       *(rgrid(ir+1)**2-2.0D0/dble(itaylor+2)*dr*rgrid(ir+1)
     &         +2.0D0/dble((itaylor+2)*(itaylor+3))*dr**2)! nth order Taylor
          end do
        end do
        return
        end function int_radial_taylorn

c---------------------------------------------------------------------

        subroutine fourier_1d_simple(rgrid,rfield,qgrid,qfield,           &
     &       qmax_chosen,nq_chosen)
        ! performs a simple Fourier trafo of a 1d function "rfield" on a 1d
        ! grid "rgrid" to "qfield" which lives on "qgrid" 
        use defs
        implicit none
        double precision, intent(in) :: rgrid(:),rfield(:)
        double precision, allocatable :: qgrid(:)
        complex, allocatable :: qfield(:)
        double precision dr,dq,rrange,qmax
        double precision, optional :: qmax_chosen
        integer nr,ir,nq,iq
        integer, optional :: nq_chosen

        nr=size(rgrid)
        nq=nr
        if (present(nq_chosen)) nq=nq_chosen
        allocate(qgrid(nq))
        allocate(qfield(nq))
        qfield=0.0D0
        rrange=maxval(rgrid)-minval(rgrid)
        qmax=2.0d0*6.0d0*2.0d0*pi/rrange
        if (present(qmax_chosen)) qmax=qmax_chosen
        dq=qmax/dble(nq)
        do iq=1,nq
          qgrid(iq)=dq*dble(iq-1)
        end do
        qgrid=qgrid-0.5*maxval(qgrid)
        do iq=1,nq
          do ir=1,nr-1
            dr=rgrid(ir+1)-rgrid(ir)
            qfield(iq)=qfield(iq)+rfield(ir)                              &
     &            *exp(-icomp*qgrid(iq)*rgrid(ir))*dr
          end do ! ir
          do ir=nr,nr
            dr=rgrid(ir)-rgrid(ir-1)
            qfield(iq)=qfield(iq)+rfield(ir)                              &
     &            *exp(-icomp*qgrid(iq)*rgrid(ir))*dr
          end do ! ir
          !print*,iq,qgrid(iq)
        end do ! iq
        qfield=qfield/sqrt(2.0d0*pi)
        return
        end subroutine fourier_1d_simple

c---------------------------------------------------------------------

        subroutine fourier_1d_reverse_simple(qgrid,qfield,rgrid,rfield)
          ! performs a simple Fourier trafo of a 1d function "qfield" on a 1d
        ! grid "qgrid" to "rfield" which lives on "rgrid" 
        use defs
        implicit none
        double precision, intent(in) :: qgrid(:)
        double precision, allocatable :: rgrid(:)
        complex, allocatable :: qfield(:),rfield(:)
        double precision dr,dq,qrange
        integer nr,ir,nq,iq

        nq=size(qgrid)
        nr=nq
        allocate(rgrid(nr))
        allocate(rfield(nr))
        rfield=0.0D0
        !qrange=maxval(qgrid)-minval(qgrid)
        dr=1.0d0/dble(nq) !qrange*(2.0d0*pi)
        do ir=1,nr
          rgrid(ir)=dr*dble(ir-1)
        end do
        do ir=1,nr
          do iq=1,nq-1
            dq=qgrid(iq+1)-qgrid(iq)
            rfield(ir)=rfield(ir)+qfield(iq)                              &
     &            *exp(icomp*qgrid(iq)*rgrid(ir))*dq
          end do ! iq
          do iq=nq,nq
            dq=qgrid(iq)-qgrid(iq-1)
            rfield(ir)=rfield(ir)+qfield(iq)                              &
     &            *exp(icomp*qgrid(iq)*rgrid(ir))*dq
          end do ! iq
          !print*,ir,rgrid(ir)
        end do ! ir
        rfield=rfield/sqrt(2.0d0*pi)
        return
        end subroutine fourier_1d_reverse_simple

c---------------------------------------------------------------------
c---------------------------------------------------------------------
        end module
