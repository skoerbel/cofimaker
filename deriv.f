      module deriv

      contains

      subroutine derivs(file1,dmode)
      use defs
      implicit none
      ! calculates the first derivative of the second column wrt the
      ! first one from finite differences 
      !
      character(len=*) file1,dmode
      ! internal variables
      double precision, allocatable :: rgrid(:),field(:)
      double precision, allocatable :: field_deriv(:) ! derivative
      character(len=256) line
      integer ndata ! number of data in file1 
      integer idata 
      double precision h ! stepsize 
      character(len=256) file1a
      !
      if(talk) print fsubstart,"derivs" 
      if (talk) print '(8x,"file: ",A40)',file1
      if (talk) print '(8x,"derivative mode: ",A40)',dmode
      !
      ! begin read grid and field from file
      open(51,file=file1,status="old",err=101)
      ndata=0
 10   read(51,'(A256)',end=11,err=101) line
      line=adjustl(line)
      if(len_trim(line).gt.0.and.line(1:1).ne."#") ndata=ndata+1
      goto 10
 11   rewind(51)
      allocate(rgrid(ndata),field(ndata))
      idata=1
 12   read(51,'(A256)',end=13,err=101) line     
      line=adjustl(line)
      if(len_trim(line).gt.0.and.line(1:1).ne."#") then
        read(line(1:256),*) rgrid(idata),field(idata)
        idata=idata+1
      end if
      goto 12
 13   close(51)
      ! end read grid and field from file
      !
      ! calculate derivative
      allocate(field_deriv(ndata))
      field_deriv=0.0D0
      select case(dmode)
      case('cd') ! first derivative, central differences 
        if (ndata.le.1) goto 300
        ! central difference derivative except for the borders:
        do idata=2,ndata-1
          field_deriv(idata)=(field(idata+1)
     &      -field(idata-1))
     &      /(rgrid(idata+1)-rgrid(idata-1))
        end do
        ! one-sided derivative at the borders:
        do idata=1,1
          field_deriv(idata)=(field(idata+1)
     &      -field(idata))
     &      /(rgrid(idata+1)-rgrid(idata))
        end do
        do idata=ndata,ndata
          field_deriv(idata)=(field(idata)
     &      -field(idata-1))
     &      /(rgrid(idata)-rgrid(idata-1))
        end do
      case('cd2') ! second derivative, central diff.
        if (ndata.le.1) goto 300
        ! central second order difference derivative except for the borders:
        do idata=2,ndata-1
          h=0.5d0*(rgrid(idata+1)-rgrid(idata-1))
          field_deriv(idata)=(field(idata+1)      
     &      -2.0d0*field(idata)+field(idata-1))/h**2 
        end do
        ! one-sided second order derivative at the borders:
        do idata=1,1
          h=0.5d0*(rgrid(idata+2)-rgrid(idata))
          field_deriv(idata)=(field(idata+2)      
     &      -2.0d0*field(idata+1)+field(idata))/h**2           
        end do
        do idata=ndata,ndata
          h=0.5d0*(rgrid(idata)-rgrid(idata-2))
          field_deriv(idata)=(field(idata-2)      
     &      -2.0d0*field(idata-1)+field(idata))/h**2           
        end do
      case('fps') ! 5-point stencil for 1. deriv.
        if (ndata.le.2) goto 300
        ! central difference first derivative except for the borders:
        do idata=3,ndata-2
          h=0.5d0*(rgrid(idata+1)-rgrid(idata-1))
          field_deriv(idata)=(-field(idata+2)+8.0d0*field(idata+1)      
     &      -8.0d0*field(idata-1)+field(idata-2))/(12.0d0*h) 
        end do
        ! one-sided second order derivative at the borders:
        do idata=1,2
          h=0.5d0*(rgrid(idata+2)-rgrid(idata))
          field_deriv(idata)=(field(idata+1)
     &      -field(idata))
     &      /(rgrid(idata+1)-rgrid(idata))
        end do
        do idata=ndata-1,ndata
          h=0.5d0*(rgrid(idata)-rgrid(idata-2))
          field_deriv(idata)=(field(idata)
     &      -field(idata-1))
     &      /(rgrid(idata)-rgrid(idata-1))
        end do
      case default
        goto 200 
      end select  
      !
      ! begin write deriv. to disk 
      write(file1a,'(a,a)') adjustl(trim(file1)),".deriv" 
      open(51,file=trim(file1a),status="replace")
      do idata=1,ndata
        write(51,'(F10.6,F15.6)') rgrid(idata),field_deriv(idata)       
      end do 
      close(51)
      ! end write deriv. to disk 
      !
      ! end normally:
      if(talk) print fsubendext, "derivs"
      return
      !
      ! Errors:
 101  nerr=nerr+1
      print ferrmssg," when opening/reading file1"
      close(51)
      return
 102  nerr=nerr+1
      print ferrmssg," when opening/reading file2"
      close(51)
      return
 200  nerr=nerr+1
      print ferrmssg," unknown deriv_mode. Consult help (cofima --help)"
      close(51)
      return
 300  nerr=nerr+1
      print ferrmssg," data array too short for this deriv_mode."
      return
      !
      end subroutine derivs

!----------------------------------------------------------------------------

      subroutine integrate(file1,imode)
      use defs
      implicit none
      ! calculates the integral of the second column wrt the
      ! first one  
      !
      character(len=*) file1,imode
      ! internal variables
      double precision, allocatable :: rgrid(:),field(:)
      double precision, allocatable :: field_integ(:) ! integrated function
      character(len=256) line
      integer ndata ! number of data in file1 
      integer idata 
      double precision h ! stepsize 
      character(len=256) file1a
      !
      if(talk) print fsubstart,"integrate" 
      if (talk) print '(8x,"file: ",A40)',file1
      if (talk) print '(8x,"integration mode: ",A40)',imode
      !
      ! begin read grid and field from file
      open(51,file=file1,status="old",err=101)
      ndata=0
 10   read(51,'(A256)',end=11,err=101) line
      line=adjustl(line)
      if(line(1:1).ne."#") ndata=ndata+1
      goto 10
 11   rewind(51)
      allocate(rgrid(ndata),field(ndata))
      idata=1
 12   read(51,'(A256)',end=13,err=101) line     
      line=adjustl(line)
      if(line(1:1).ne."#") then
        read(line(1:256),*) rgrid(idata),field(idata)
        idata=idata+1
      end if
      goto 12
 13   close(51)
      ! end read grid and field from file
      !
      ! calculate integrated function
      allocate(field_integ(ndata))
      field_integ=0.0D0
      select case(imode)
      case('riemann') ! Riemann sum
        if (ndata.lt.1) goto 300
        do idata=2,ndata
          field_integ(idata)=field_integ(idata-1)+(rgrid(idata)           &
     &      -rgrid(idata-1))*field(idata)
        end do
      case default
        goto 200 
      end select  
      !
      ! begin write integ. to disk 
      write(file1a,'(a,a)') adjustl(trim(file1)),".integ" 
      open(51,file=trim(file1a),status="replace")
      do idata=1,ndata
        write(51,'(F10.6,F15.6)') rgrid(idata),field_integ(idata)       
      end do 
      close(51)
      ! end write integ. to disk 
      !
      ! end normally:
      if(talk) print fsubendext, "integrate"
      return
      !
      ! Errors:
 101  nerr=nerr+1
      print ferrmssg," when opening/reading file1"
      close(51)
      return
 102  nerr=nerr+1
      print ferrmssg," when opening/reading file2"
      close(51)
      return
 200  nerr=nerr+1
      print ferrmssg," unknown integ_mode. Consult help (cofima --help)"
      close(51)
      return
 300  nerr=nerr+1
      print ferrmssg," data array too short for this integ_mode."
      return
      !
      end subroutine integrate

      end module 
