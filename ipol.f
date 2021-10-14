        module ipol
        implicit none

        contains 

        SUBROUTINE spline(x,y,n,yp1,ypn,y2)
        INTEGER n,NMAX
        DOUBLE PRECISION yp1,ypn,x(n),y(n),y2(n)
        PARAMETER (NMAX=500)
        ! Given arrays x(1:n) and y(1:n) containing a tabulated
        ! function, i.e., yi = f (xi ), with
        ! x1 < x2 < . . . < xN , and given values yp1 and ypn for
        ! the first derivative of the inter-
        ! polating function at points 1 and n, respectively, this
        ! routine returns an array y2(1:n) of
        ! length n which contains the second derivatives of the
        ! interpolating function at the tabulated
        ! points xi . If yp1 and/or ypn are equal to 1 × 1030 or
        ! larger, the routine is signaled to set
        ! the corresponding boundary condition for a natural
        ! spline, with zero second derivative on
        ! that boundary.
        ! Parameter: NMAX is the largest anticipated value of n.
        INTEGER i,k
        DOUBLE PRECISION p,qn,sig,un,u(NMAX)
        if (yp1.gt..99d30) then
                !The lower boundary condition is set either to be
                !“natural”
                y2(1)=0.
                u(1)=0.
        else
                !or else to have a specified first derivative.
                y2(1)=-0.5
                u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
        endif
        !do 11 i=2,n-1
        do i=2,n-1
                ! This is the decomposition loop of the tridiagonal
                ! algorithm. y2 and u are used for temporary
                ! storage of the decomposed factors.
                sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
                p=sig*y2(i-1)+2.
                y2(i)=(sig-1.)/p
                u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     &               /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
        !enddo 11 
        enddo 
        if (ypn.gt..99d30) then
                !The upper boundary condition is set either to be
                !“natural”
                qn=0.
                un=0.
        else
                !or else to have a specified first derivative.
                qn=0.5
                un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
        endif
        y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
        !do 12 k=n-1,1,-1
        do k=n-1,1,-1
                !This is the backsubstitution loop of the tridiago-
                !nal algorithm.
                y2(k)=y2(k)*y2(k+1)+u(k)
        !enddo 12
        enddo 
        return
        END SUBROUTINE

        SUBROUTINE splie2(x1a,x2a,ya,m,n,y2a)
        INTEGER m,n,NN
        DOUBLE PRECISION x1a(m),x2a(n),y2a(m,n),ya(m,n)
        PARAMETER (NN=300)
        !Maximum expected value of n and m.
        !USES spline
        !Given an m by n tabulated function ya(1:m,1:n), and tabulated
        !independent variables
        !x2a(1:n), this routine constructs one-dimensional natural cubic
        !splines of the rows of ya
        !and returns the second-derivatives in the array y2a(1:m,1:n).
        !(The array x1a is included
        !in the argument list merely for consistency with routine
        !splin2.)
        INTEGER j,k
        DOUBLE PRECISION y2tmp(NN),ytmp(NN)
        !do 13 j=1,m
        do j=1,m
                !do 11 k=1,n
                !do k=1,n
                !        ytmp(k)=ya(j,k)
                !enddo !11
                ytmp(1:n)=ya(j,1:n)
                !call spline(x2a,ytmp,n,1.e30,1.e30,y2tmp)
                call spline(x2a,ytmp,n,1.d30,1.d30,y2tmp)
                !Values 1×1030 signal a natural spline.
                !do 12 k=1,n
                !do k=1,n
                !        y2a(j,k)=y2tmp(k)
                !enddo !12
                y2a(j,1:n)=y2tmp(1:n)
        enddo !13
        return
        END SUBROUTINE

        SUBROUTINE splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
        INTEGER m,n,NN
        DOUBLE PRECISION x1,x2,y,x1a(m),x2a(n),y2a(m,n),ya(m,n)
        PARAMETER (NN=300)
        !Maximum expected value of n and m.
        !  USES spline,splint
        !  Given x1a, x2a, ya, m, n as described in splie2 and y2a as
        !produced by that routine;
        ! and given a desired interpolating point x1,x2; this routine
        !returns an interpolated function
        !value y by bicubic spline interpolation.
        INTEGER j,k
        DOUBLE PRECISION y2tmp(NN),ytmp(NN),yytmp(NN) 
        !do 12 j=1,m
        do j=1,m
        !Perform m evaluations of the row splines
                !do 11 k=1,n
                !do k=1,n
                !        !constructed by splie2, using the one-
                !        !dimensional spline evaluator splint.
                !        ytmp(k)=ya(j,k)
                !        y2tmp(k)=y2a(j,k)
                !enddo !11
                ytmp(1:n)=ya(j,1:n)
                y2tmp(1:n)=y2a(j,1:n)
                call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
        enddo !12
        !Construct the one-dimensional column spline
        !and evaluate it.
        !call spline(x1a,yytmp,m,1.e30,1.e30,y2tmp)
        call spline(x1a,yytmp,m,1.d30,1.d30,y2tmp)
        call splint(x1a,yytmp,y2tmp,m,x1,y)
        return
        END SUBROUTINE

!---------------------------------------------------------------------        

        SUBROUTINE splint(xa,ya,y2a,n,x,y)
        use defs
        implicit none        
        INTEGER n
        DOUBLE PRECISION x,y,xa(n),y2a(n),ya(n)
        ! Given the arrays xa(1:n) and ya(1:n) of length n, which
        ! tabulate a function (with the
        ! xai ’s in order), and given the array y2a(1:n), which is
        ! the output from spline above,
        ! and given a value of x, this routine returns a
        ! cubic-spline interpolated value y.
        INTEGER k,khi,klo
        DOUBLE PRECISION a,b,h
        klo=1
        ! We will find the right place in the table by means of
        ! bisection.
        khi=n
        ! This is optimal if sequential calls to this routine are
        ! at random
 1      if (khi-klo.gt.1) then
                ! values of x. If sequential calls are in order,
                ! and closely
                ! spaced, one would do better to store previous
                ! values of
                ! klo and khi and test if they remain
                ! appropriate on the
                ! next call.
                k=(khi+klo)/2
                if(xa(k).gt.x)then
                        khi=k
                else
                        klo=k
                endif
                        goto 1
        endif
        ! klo and khi now bracket the input value of x.
        h=xa(khi)-xa(klo)
        if (h.eq.0.) then
          print ferrmssg,"bad xa input in splint" 
          return 
        end if   
        ! The xa’s
        ! must be distinct.
        a=(xa(khi)-x)/h
        ! Cubic spline polynomial is now evaluated.
        b=(x-xa(klo))/h
        y=a*ya(klo)+b*ya(khi)+
     &   ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
        return
        END SUBROUTINE
        
!---------------------------------------------------------------------

        subroutine qspline_setup(x,y,yp1,y2)
        use linalg
        implicit none
        ! Given arrays x(1:n) and y(1:n) containing a tabulated
        ! function, i.e., yi = f (xi ), with
        ! x1 < x2 < . . . < xN , this
        ! routine returns an array d(1:n) of
        ! length n which contains the first derivatives (*2) of the
        ! interpolating function at the tabulated
        ! points xi . If yp1 is equal to 1 × 10^30 or
        ! larger, the routine is signaled to set
        ! the corresponding boundary condition for a natural
        ! spline, with zero derivative on
        ! that boundary.
        INTEGER n
        DOUBLE PRECISION x(:),y(:),yp1,y2(:)
        ! internal variables
        double precision, allocatable :: mat(:,:),mati(:,:),d(:) 
        INTEGER i,k
        !
        ! begin test consistence of input parameters
        n=size(x)
        if(size(y).ne.n) stop
        if(size(y2).ne.n) stop
        ! end test consistence of input parameters
        !
        ! begin set up and invert bidiagonal matrix
        allocate(mat(n,n),mati(n,n),d(n))
        mat=0.0d0
        mat(1,1)=1.0d0
        do i=2,n
          mat(i,i)=1.0d0
          mat(i,i-1)=1.0d0
        end do
        call inverse(mat,mati,n)
        print '(8x,"Inverse matrix set up.")'
        ! end set up and invert bidiagonal matrix
        !
        ! begin get d
        if (yp1.gt..99d30) then
                !The lower boundary condition is set either to be
                !“natural”
                d(1)=0.
        else
                !or else to have a specified first derivative.
                d(1)=2.0d0*(y(2)-y(1))/(x(2)-x(1))
        endif
        do i=2,n
                d(i)=2.0d0*(y(i)-y(i-1))/(x(i)-x(i-1))
        enddo 
        ! end get d
        !
        ! begin get z (y2)
        y2=matmul(mati,d)
        ! end get z (y2)
        !
        deallocate(mat,mati,d)
        return
        END SUBROUTINE qspline_setup

!---------------------------------------------------------------------

        SUBROUTINE qspline(xa,ya,y2a,x,y)
        use defs
        implicit none
        ! Given the arrays xa(1:n) and ya(1:n) of length n, which
        ! tabulate a function (with the
        ! xai ’s in order), and given the array y2a(1:n), which is
        ! the output from qspline_setup above,
        ! and given a value of x, this routine returns a
        ! quadratic-spline interpolated value y.
        DOUBLE PRECISION x,y,xa(:),y2a(:),ya(:)
        INTEGER k,khi,klo
        DOUBLE PRECISION a,b,h
        ! internal variables
        INTEGER n
        !
        n=size(xa)
        klo=1
        ! We will find the right place in the table by means of
        ! bisection.
        khi=n
        ! This is optimal if sequential calls to this routine are
        ! at random
 1      if (khi-klo.gt.1) then
                ! values of x. If sequential calls are in order,
                ! and closely
                ! spaced, one would do better to store previous
                ! values of
                ! klo and khi and test if they remain
                ! appropriate on the
                ! next call.
                k=(khi+klo)/2
                if(xa(k).gt.x)then
                        khi=k
                else
                        klo=k
                endif
                        goto 1
        endif
        ! klo and khi now bracket the input value of x.
        h=xa(khi)-xa(klo)
        if (h.eq.0.) then
            call error("bad xa input in qspline") 
        end if
        ! The xa’s
        ! must be distinct.
        ! Quadratic spline polynomial is now evaluated.
        y=y2a(klo)*(x-xa(klo))+(x-xa(klo))**2*(y2a(khi)-y2a(klo))         &
     &         /(2.0d0*h) +ya(klo)         
        return
        END SUBROUTINE qspline

!---------------------------------------------------------------------

        subroutine akima_setup(x,y,y2)
        use linalg
        implicit none
        ! from https://www.ads.tuwien.ac.at/docs/lva/mmgdv/k1___011.htm
        ! Given arrays x(1:n) and y(1:n) containing a tabulated
        ! function, i.e., yi = f (xi ), with
        ! x1 < x2 < . . . < xN , this
        ! routine returns an array y2(1:n) of
        ! length n which contains the estimated first derivatives of the
        ! function at the tabulated
        ! points xi . 
        DOUBLE PRECISION x(:),y(:),y2(:)
        ! internal variables
        INTEGER n
        INTEGER i,k
        double precision, allocatable :: q(:)
        double precision dq12,dq01
        !
        ! begin test consistence of input parameters
        n=size(x)
        if(size(y).ne.n) stop
        if(size(y2).ne.n) stop
        ! end test consistence of input parameters
        !
        ! begin get q (slopes)
        allocate(q(n))
        q(1)=(y(2)-y(1))/(x(2)-x(1))
        do i=2,n
          q(i)=(y(i)-y(i-1))/(x(i)-x(i-1))
        enddo 
        ! end get q (slopes)
        !
        ! begin get y2 (estimated slope)
        y2(1)=q(1)
        do i=2,n-2
          dq12=abs(q(i+2)-q(i+1))
          dq01=abs(q(i)-q(i-1))
          if(dq12.gt.1.0d-6.or.dq01.gt.1.0d-6) then
             y2(i)=(q(i)*dq12+q(i+1)*dq01)/(dq12+dq01)
          else
             y2(i)=0.5d0*(q(i)+q(i+1))   
          end if
        enddo 
        y2(n-1)=0.5d0*(q(n)+q(n-1))   
        y2(n)=q(n)
        deallocate(q)
        ! end get y2 (estimated slope)
        !
        return
        END SUBROUTINE akima_setup

!---------------------------------------------------------------------

        SUBROUTINE akima(xa,ya,y2a,x,y)
        use defs
        implicit none
        ! Given the arrays xa(1:n) and ya(1:n) of length n, which
        ! tabulate a function (with the
        ! xai ’s in order), and given the array y2a(1:n), which is
        ! the output from akima_setup above,
        ! and given a value of x, this routine returns a
        ! cubic Akima-spline interpolated value y.
        DOUBLE PRECISION x,y,xa(:),y2a(:),ya(:)
        INTEGER k,khi,klo
        ! internal variables
        DOUBLE PRECISION h,t,g0,g1,h0,h1,pdotl,pdotr
        INTEGER n
        !
        n=size(xa)
        klo=1
        ! We will find the right place in the table by means of
        ! bisection.
        khi=n
        ! This is optimal if sequential calls to this routine are
        ! at random
 1      if (khi-klo.gt.1) then
                ! values of x. If sequential calls are in order,
                ! and closely
                ! spaced, one would do better to store previous
                ! values of
                ! klo and khi and test if they remain
                ! appropriate on the
                ! next call.
                k=(khi+klo)/2
                if(xa(k).gt.x)then
                        khi=k
                else
                        klo=k
                endif
                        goto 1
        endif
        ! klo and khi now bracket the input value of x.
        h=xa(khi)-xa(klo)
        t=(x-xa(klo))/(xa(khi)-xa(klo))
        pdotl=h*y2a(klo)
        pdotr=h*y2a(khi)
        g0=0.5d0-1.5d0*(t-0.5d0)+2.0d0*(t-0.5d0)**3
        g1=1.0d0-g0 
        h0=t*(t-1.0d0)**2 
        h1=t**2*(t-1.0d0)
        if (h.eq.0.) then
                call error("bad xa input in akima")
        end if        
        ! The xa’s
        ! must be distinct.
        ! Cubic Akima spline polynomial is now evaluated.
        y=ya(klo)*g0+ya(khi)*g1+pdotl*h0+pdotr*h1
        return
        END SUBROUTINE akima

!---------------------------------------------------------------------
        

!---------------------------------------------------------------------
        
        SUBROUTINE bcucof(y,y1,y2,y12,d1,d2,c)
        DOUBLE PRECISION d1,d2,c(4,4),y(4),y1(4),y12(4),y2(4)
        !Given arrays y,y1,y2, and y12, each of length 4,
        !containing the function, gradients, and
        !cross derivative at the four grid points of a
        !rectangular grid cell (numbered counterclockwise
        !from the lower left), and given d1 and d2, the length of
        !the grid cell in the 1- and 2-
        !directions, this routine returns the table c(1:4,1:4)
        !that is used by routine bcuint for
        !bicubic interpolation.
        INTEGER i,j,k,l
        DOUBLE PRECISION d1d2,xx,cl(16),wt(16,16),x(16)
        SAVE wt
        DATA wt/1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,8*0,3,0,-9,6,-2,0,6,-4
     *       ,10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,2*0,6,-4
     *       ,4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,2
     *       ,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2
     *       ,0,1,-2,1,5*0,-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2
     *       ,10*0,-3,3,2*0,2,-2,2*0,-1,1,6*0,3,-3,2*0,-2,2
     *       ,5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,-1,0,1,-2,1
     *       ,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/
        d1d2=d1*d2
        !do 11 i=1,4
        do i=1,4
                ! Pack a temporary vector x.
                x(i)=y(i)
                x(i+4)=y1(i)*d1
                x(i+8)=y2(i)*d2
                x(i+12)=y12(i)*d1d2
        enddo !11
        !do 13 i=1,16 ! Matrix multiply by the stored table.
        do i=1,16
                xx=0.
                !do 12 k=1,16
                do k=1,16
                        xx=xx+wt(i,k)*x(k)
                enddo !12
                cl(i)=xx
        enddo !13
        l=0
        !do 15 i=1,4 ! Unpack the result into the output table.
        do i=1,4
                !do 14 j=1,4
                do j=1,4
                        l=l+1
                        c(i,j)=cl(l)
                enddo !14
        enddo !15
        return
        END SUBROUTINE

        
        SUBROUTINE bcuint(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,ansy,
     &           ansy1,ansy2)
        use defs
        implicit none        
        DOUBLE PRECISION ansy,ansy1,ansy2,x1,x1l,x1u,x2,x2l,x2u,y(4),
     &                   y1(4),y12(4),y2(4)
        !USES bcucof
        !Bicubic interpolation within a grid square. Input
        !quantities are y,y1,y2,y12 (as described
        !in bcucof); x1l and x1u, the lower and upper coordinates
        !of the grid square in the 1-
        !direction; x2l and x2u likewise for the 2-direction; and
        !x1,x2, the coordinates of the
        !desired point for the interpolation. The interpolated
        !function value is returned as ansy,
        !and the interpolated gradient values as ansy1 and ansy2.
        !This routine calls bcucof.
        INTEGER i
        DOUBLE PRECISION t,u,c(4,4)
        call bcucof(y,y1,y2,y12,x1u-x1l,x2u-x2l,c)
        !Get the c’s.
        if(x1u.eq.x1l.or.x2u.eq.x2l) then
          print ferrmssg, "bad input in bcuint"
          return
        end if  
        t=(x1-x1l)/(x1u-x1l)
        !Equation (3.6.4).
        u=(x2-x2l)/(x2u-x2l)
        ansy=0.
        ansy2=0.
        ansy1=0.
        !do 11 i=4,1,-1
        do i=4,1,-1
                !Equation (3.6.6).
                ansy=t*ansy+((c(i,4)*u+c(i,3))*u+c(i,2))*u+c(i,1)
                ansy2=t*ansy2+(3.*c(i,4)*u+2.*c(i,3))*u+c(i,2)
                ansy1=u*ansy1+(3.*c(4,i)*t+2.*c(3,i))*t+c(2,i)
        enddo !11
        ansy1=ansy1/(x1u-x1l)
        ansy2=ansy2/(x2u-x2l)
        return
        END SUBROUTINE

c---------------------------------------------------------------------
        !
        subroutine polint1d(x,y,x0,y0)
        implicit none
        double precision, intent(in) :: x(:),y(:,:),x0
        double precision, intent(inout) :: y0(:)
        ! internal:
        integer nx,ix,jx,dimy
        double precision prod
        !
        ! begin testing input values:
        nx=size(x)
        if (size(y,1).ne.nx) then
          goto 100
        end if
        dimy=size(y,2)
        if (size(y0,1).ne.dimy) then
          goto 100
        end if
        ! end testing input values
        !
        ! begin calculate interpolated function value y0 at x0
        y0=0.0d0
        do ix=1,nx
          prod=1.0d0
          do jx=1,nx
            if (jx.ne.ix) then
              prod=prod*(x0-x(jx))/(x(ix)-x(jx))
            end if
          end do ! jx
          y0(:)=y0(:)+y(ix,:)*prod
        end do ! ix
        ! end calculate interpolated function value y0 at x0
        !
        ! end normally:
        return
        !
        ! errors
 100    print*, "Error in polint1d: wrong dimension of function array"
        stop
        !
        end subroutine polint1d
        !
c---------------------------------------------------------------------
        !
        !
        subroutine polint1dprime(x,y,x0,y0prime)
        implicit none
        double precision, intent(in) :: x(:),y(:,:),x0
        double precision, intent(inout) :: y0prime(:)
        ! internal:
        integer nx,ix,jx,kx,dimy
        double precision prod,summ
        !
        ! begin testing input values:
        nx=size(x)
        if (size(y,1).ne.nx) then
          goto 100
        end if
        dimy=size(y,2)
        if (size(y0prime,1).ne.dimy) then
          goto 100
        end if
        ! end testing input values
        !
        ! begin calculate interpolated function value y0 at x0
        y0prime=0.0d0
        do ix=1,nx
          summ=0.0d0
          do jx=1,nx
            if (jx.ne.ix) then
              prod=1.0d0
              do kx=1,nx
                if (kx.ne.ix.and.kx.ne.jx) then
                  prod=prod*(x0-x(kx))/(x(ix)-x(kx))
                end if
              end do ! kx
              summ=summ+prod/(x(ix)-x(jx))
            end if ! jx/=ix
          end do ! jx
          y0prime(:)=y0prime(:)+y(ix,:)*summ
        end do ! ix
        ! end calculate interpolated function derivative y0prime at x0
        !
        ! end normally:
        return
        !
        ! errors
 100    print*, "Error in polint1d: wrong dimension of function array"
        stop
        !
        end subroutine polint1dprime
        !
c---------------------------------------------------------------------
        !
        subroutine linint1d(grid0,field0,x0,y0)
        ! linear interpolation on 1d grid
        use defs
        implicit none
        double precision, intent(in) :: grid0(:),field0(:,:),x0
        double precision, intent(inout) :: y0(:)
        ! internal:
        double precision, allocatable :: grid(:),field(:,:) ! like grid0 and field0, but sorted
        double precision, allocatable :: b(:) 
        double precision a,h
        integer nx,ix,jx,dimy,klo,khi,k
        double precision rtemp,ftemp
        !
        ! begin testing input values:
        nx=size(grid0)
        if (size(field0,1).ne.nx) then
          goto 100
        end if
        dimy=size(field0,2)
        if (size(y0,1).ne.dimy) then
          goto 100
        end if
        ! end testing input values
        !
        ! begin sort grid points 
        allocate(grid(nx),field(nx,dimy))
        allocate(b(dimy))
        grid=grid0
        field=field0
        do ix=1,nx       
          do jx=ix+1,nx
            if(grid(jx).lt.grid(ix)) then
                  rtemp=grid(ix)
                  ftemp=field(ix,1)
                  grid(ix)=grid(jx)
                  field(ix,1)=field(jx,1)
                  grid(jx)=rtemp
                  field(jx,1)=ftemp
            end if
          end do
        end do
        ! end sort grid points
        !
        ! begin calculate interpolated function value y0 at x0
        y0=0.0d0
        ! find the nearest two grid points by bysection
        klo=1
        ! We will find the right place in the table by means of
        ! bisection.
        khi=nx
        ! This is optimal if sequential calls to this routine are
        ! at random
 1      if (khi-klo.gt.1) then
                ! values of x. If sequential calls are in order,
                ! and closely
                ! spaced, one would do better to store previous
                ! values of
                ! klo and khi and test if they remain
                ! appropriate on the
                ! next call.
                k=(khi+klo)/2
                if(grid(k).gt.x0)then
                        khi=k
                else
                        klo=k
                endif
                        goto 1
        endif
        ! klo and khi now bracket the input value of x.
        h=grid(khi)-grid(klo)
        ! The grid points
        ! must be distinct.
        if (h.eq.0.) then
          print ferrmssg, "bad grid input in linint1d" 
          return
        end if  
        !
        ! Linear interpolation is now evaluated.
        a=(x0-grid(klo))/h
        b(:)=field(khi,:)-field(klo,:)
        y0(:)=field(klo,:)+a*b(:) 
        ! end calculate interpolated function value y0 at x0
        !
        ! end normally:
        return
        !
        ! errors
 100    print*, "Error in linint1d: wrong dimension of function array"
        stop
        !
        end subroutine linint1d
        !
!---------------------------------------------------------------------
        !
        subroutine gaussint(x,y,gc,lambda,x0,y0)
        implicit none
        double precision, intent(in) :: x(:),y(:,:),x0,gc(:),lambda ! grid, field, position, gaussian coefficients, gaussian decay constant 
        double precision, intent(inout) :: y0(:) ! interpolated function value
        ! internal:
        integer nx,ix,jx,dimy
        !
        ! begin testing input values:
        nx=size(x)
        if (size(y,1).ne.nx) then
          goto 100
        end if
        dimy=size(y,2)
        if (size(y0,1).ne.dimy) then
          goto 100
        end if
        if (size(gc).ne.nx) then
          goto 100
        end if
        ! end testing input values
        !
        ! begin calculate interpolated function value y0 at x0
        y0=0.0d0
        do ix=1,nx
          y0(:)=y0(:)+y(ix,:)*gc(ix)*exp(-(x0-x(ix))**2/lambda**2)
        end do ! ix
        ! end calculate interpolated function value y0 at x0
        !
        ! end normally:
        return
        !
        ! errors
 100    print*, "Error in gaussint: wrong dimension of function array"
        stop
        !
        end subroutine gaussint
        !
!---------------------------------------------------------------------
        !
        subroutine fourierinteven(gc,gc2,length,ndata,x,x0,y0)
        use defs
        implicit none
        double precision, intent(in) :: x,x0,gc(:),gc2(:) ! position, first grid point, cosine coefficients, sine coefficients, position
        double precision, intent(in) :: length  ! periodicity length
        integer, intent(in) :: ndata ! number of original data points
        double precision, intent(inout) :: y0(:) ! interpolated function value
        ! internal:
        integer nc,ic,jx,dimy
        !
        nc=size(gc)
        ! begin calculate interpolated function value y0 at x0
        y0=0.0d0
        ic=1 
        y0(:)=y0(:)+gc(ic) 
        do ic=2,nc-1
          y0(:)=y0(:)+2.0d0*gc(ic)*cos(2.0d0*Pi*dble(ic-1)              &
     &         *(x-x0)/length)
          y0(:)=y0(:)+2.0d0*gc2(ic)*sin(2.0d0*Pi*dble(ic-1)             &
     &         *(x-x0)/length)
        end do ! ic
        ic=nc
        y0(:)=y0(:)+gc(ic)*cos(Pi*dble(ndata)*(x-x0)/length)
        ! end calculate interpolated function value y0 at x0
        !
        ! end normally:
        return
        !
        end subroutine fourierinteven
        !
!---------------------------------------------------------------------
        !
        subroutine fourierintodd(gc,gc2,length,ndata,x,x0,y0)
        use defs
        implicit none
        double precision, intent(in) :: x,x0,gc(:),gc2(:) ! position, first grid point, cosine coefficients, sine coefficients, position
        double precision, intent(in) :: length  ! periodicity length
        integer, intent(in) :: ndata ! number of original data points
        double precision, intent(inout) :: y0(:) ! interpolated function value
        ! internal:
        integer nc,ic,jx,dimy
        !
        nc=size(gc)
        ! begin calculate interpolated function value y0 at x0
        y0=0.0d0
        ic=1 
        y0(:)=y0(:)+gc(ic) 
        do ic=2,nc
          y0(:)=y0(:)+2.0d0*gc(ic)*cos(2.0d0*Pi*dble(ic-1)              &
     &         *(x-x0)/length)
          y0(:)=y0(:)+2.0d0*gc2(ic)*sin(2.0d0*Pi*dble(ic-1)             &
     &         *(x-x0)/length)
        end do ! ic
        ! end calculate interpolated function value y0 at x0
        !
        ! end normally:
        return
        !
        end subroutine fourierintodd
        !
c---------------------------------------------------------------------
      !
      subroutine interpol(file1,ipolmeth,ipoln)
      use defs
      use linalg
      implicit none
      ! interpolates the second column wrt the
      ! first one using ipolmeth 
      !
      character(len=*) file1,ipolmeth
      integer ipoln
      ! internal variables
      double precision, allocatable :: rgrid(:),field(:,:)
      double precision, allocatable :: fieldi(:,:) ! interpolated field
      character(len=256) line
      integer ndata,ndatai ! number of data in file1, interpolated grid 
      integer idata,jdata,kdata
      double precision, allocatable :: rgridi(:) ! fine grid for interpolation
      double precision, allocatable :: y2(:) ! 2. derivatives, needed for cspline
      character(len=256) file1i ! file with interpolated values
      double precision, allocatable :: gmat(:,:),gmati(:,:),gcoeffs(:) ! for gaussian interpolation 
      double precision, allocatable :: gmat2(:,:),gcoeffs2(:) ! for Fourier interpolation 
      double precision lambda ! for gaussian interpolation
      double precision plength ! periodicity length for Fourier interpolation
      double precision rtemp,ftemp
      !
      if(talk) print fsubstart,"interpol" 
      if (talk) print '(8x,"file: ",A40)',file1
      if (talk) print '(8x,"ipolmeth: ",A40)',ipolmeth
      !
      ! begin read grid and field from file
      open(51,file=file1,status="old",err=101)
      ndata=0
 10   read(51,'(A256)',end=11,err=101) line
      line=adjustl(line)
      if(line(1:1).ne."#") ndata=ndata+1
      goto 10
 11   rewind(51)
      allocate(rgrid(ndata),field(ndata,1))
      idata=1
 12   read(51,'(A256)',end=13,err=101) line     
      line=adjustl(line)
      if(line(1:1).ne."#") then
        read(line(1:256),*) rgrid(idata),field(idata,1)
        idata=idata+1
      end if
      goto 12
 13   close(51)
      ! end read grid and field from file
      !
      ! begin sort grid points 
      do idata=1,ndata       
        do jdata=idata+1,ndata
          if(rgrid(jdata).lt.rgrid(idata)) then
                rtemp=rgrid(idata)
                ftemp=field(idata,1)
                rgrid(idata)=rgrid(jdata)
                field(idata,1)=field(jdata,1)
                rgrid(jdata)=rtemp
                field(jdata,1)=ftemp
          end if
        end do
      end do
      ! end sort grid points
      !
      ! begin set up fine interpolation grid
      !ndatai=ndata*ipoln
      ndatai=(ndata-1)*ipoln+1
      allocate(rgridi(ndatai),fieldi(ndatai,1))
      kdata=0
      do idata=1,ndata-1
        do jdata=1,ipoln
          kdata=kdata+1
          rgridi(kdata)=rgrid(idata)+(rgrid(idata+1)-rgrid(idata))
     &                 /dble(ipoln) * dble(jdata-1)
        end do
      end do
      kdata=kdata+1
      rgridi(kdata)=rgrid(ndata)
      ! end set up fine interpolation grid
      ! 
      ! begin interpolate
      fieldi=0.0D0
      select case(ipolmeth)
      case('lin')
        do kdata=1,ndatai 
          call linint1d(rgrid,field(:,:),rgridi(kdata),fieldi(kdata,:))
        end do
      case('polint')
        do kdata=1,ndatai 
          call polint1d(rgrid,field(:,:),rgridi(kdata),fieldi(kdata,:))
        end do
      case('cspline')
        ! set up cubic spline 
        allocate(y2(ndata))
        call spline(rgrid(:),field(:,1),ndata,99.d30,99.d30,y2(:))
        if(talk) print '(8x,"spline set up.")' 
        do kdata=1,ndatai 
          call splint(rgrid(:),field(:,1),y2(:),ndata,rgridi(kdata),
     &                fieldi(kdata,1))
        end do
        deallocate(y2)
      case('qspline')
        ! set up quadratic spline 
        allocate(y2(ndata))
        !call qspline_setup(rgrid(:),field(:,1),99.d30,y2(:)) ! natural spline
        call qspline_setup(rgrid(:),field(:,1),1.0d0,y2(:)) ! unnatural spline
        if(talk) print '(8x,"quadratic spline set up.")' 
        do kdata=1,ndatai 
          call qspline(rgrid(:),field(:,1),y2(:),rgridi(kdata),           &
     &                fieldi(kdata,1))
        end do
        deallocate(y2)
      case('akima')
        ! set up cubic Akima spline 
        allocate(y2(ndata))
        call akima_setup(rgrid(:),field(:,1),y2(:))
        if(talk) print '(8x,"quadratic spline set up.")' 
        do kdata=1,ndatai 
          call akima(rgrid(:),field(:,1),y2(:),rgridi(kdata),           &
     &                fieldi(kdata,1))
        end do
        deallocate(y2)
      case('gauss') ! interpolate with sum of Gaussians:
        ! p(x)=sum_j a_j f_j exp(-(x-x_j)**2/lambda**2), j=1,ndata
        ! 
        ! begin set up matrix M (gmat)
        ! M_ij=f_j exp(-(x-x_j)**2/lambda**2), such that sum_j M_ij * a_j = f_i
        lambda=(rgrid(ndata)-rgrid(1))/dble(ndata) ! Gaussian decay parameter
        !lambda=(rgrid(ndata)-rgrid(1))*1.1d0/dble(ndata) ! Gaussian decay parameter, larger (smoother)
        !lambda=(rgrid(ndata)-rgrid(1))*0.9d0/dble(ndata) ! Gaussian decay parameter, smaller (less smooth)
        allocate(gmat(ndata,ndata),gmati(ndata,ndata),gcoeffs(ndata))
        gmat=0.0d0
        do idata=1,ndata
          do jdata=1,ndata
           gmat(idata,jdata)=field(jdata,1)
     &                       *exp(-(rgrid(idata)-rgrid(jdata))**2        
     &                 /lambda**2)
          end do
        end do
        ! end set up matrix M (gmat)
        ! 
        ! begin invert gmat and get coefficients  
        call inverse(gmat,gmati,ndata)
        gcoeffs=matmul(gmati,field(:,1))
        deallocate(gmat,gmati)
        ! end invert gmat and get coefficients 
        ! 
        ! begin interpolate
        do kdata=1,ndatai 
          call gaussint(rgrid,field(:,:),gcoeffs,lambda,rgridi(kdata),
     &                  fieldi(kdata,:))
        end do
        ! end interpolate
        !
        deallocate(gcoeffs)
        !
      case('rdfi') ! real discrete Fourier interpolation 
        ! N even:
        ! p(x)=a0 + 1/N 2*sum_k (a_k cos(2 pi k x/L ) + b_k sin(2 pi k x / L) ) + A_N/2 cos(N pi x/L), j=0,N/2-1
        if(mod(ndata,2).eq.0) then
          ! 
          ! begin set up matrices M (gmat, gmat2) such that
          ! M_ij=f_j cos(2 pi (i-1)(j-1) / L) or M_ij=f_j sin(2 pi (i-1)(j-1) / L), 
          ! with a = M f or b = M f
          allocate(gmat(ndata/2+1,ndata),gmat2(ndata/2+1,ndata),              &
     &             gcoeffs(ndata/2+1),gcoeffs2(ndata/2+1))
          gmat=0.0d0
          gmat2=0.0d0
          do idata=1,ndata/2+1
            do jdata=1,ndata
             gmat(idata,jdata)=1.0d0/dble(ndata) *                        &
     &          cos(2.0d0*Pi*dble((idata-1)*(jdata-1))/dble(ndata))
             gmat2(idata,jdata)=1.0d0/dble(ndata) *                       &
     &          sin(2.0d0*Pi*dble((idata-1)*(jdata-1))/dble(ndata))
            end do
          end do
          ! end set up matrices M (gmat,gmat2)
          ! 
          ! begin get coefficients  
          gcoeffs=matmul(gmat,field(:,1))
          gcoeffs2=matmul(gmat2,field(:,1))
          deallocate(gmat,gmat2)
          ! end get coefficients 
          ! 
          ! begin interpolate
          plength=(rgrid(ndata)-rgrid(1))*dble(ndata)/dble(ndata-1)
          do kdata=1,ndatai 
            call fourierinteven(gcoeffs,gcoeffs2,plength,ndata,           &
     &                      rgridi(kdata),rgrid(1),fieldi(kdata,:))
          end do
          ! end interpolate
          !
          deallocate(gcoeffs)
        end if ! N even
        !
        ! N odd:
        ! p(x)=p0 + 1/N 2*sum_k (a_k cos(2 pi k x/L ) + b_k sin(2 pi k x / L) ) , j=0,(N-1)/2
        if(mod(ndata,2).ne.0) then
          ! 
          ! begin set up matrices M (gmat, gmat2) such that
          ! M_ij=f_j cos(2 pi (i-1)(j-1) / L) or M_ij=f_j sin(2 pi (i-1)(j-1) / L), 
          ! with a = M f or b = M f
          allocate(gmat((ndata+1)/2,ndata),gmat2((ndata+1)/2,ndata),      &
     &             gcoeffs((ndata+1)/2),gcoeffs2((ndata+1)/2))
          gmat=0.0d0
          gmat2=0.0d0
          do idata=1,(ndata+1)/2
            do jdata=1,ndata
             gmat(idata,jdata)=1.0d0/dble(ndata) *                        &
     &          cos(2.0d0*Pi*dble((idata-1)*(jdata-1))/dble(ndata))
             gmat2(idata,jdata)=1.0d0/dble(ndata) *                       &
     &          sin(2.0d0*Pi*dble((idata-1)*(jdata-1))/dble(ndata))
            end do
          end do
          ! end set up matrices M (gmat,gmat2)
          ! 
          ! begin get coefficients  
          gcoeffs=matmul(gmat,field(:,1))
          gcoeffs2=matmul(gmat2,field(:,1))
          deallocate(gmat,gmat2)
          ! end get coefficients 
          ! 
          ! begin interpolate
          plength=(rgrid(ndata)-rgrid(1))*dble(ndata)/dble(ndata-1)
          do kdata=1,ndatai 
            call fourierintodd(gcoeffs,gcoeffs2,plength,ndata,           &
     &                      rgridi(kdata),rgrid(1),fieldi(kdata,:))
          end do
          ! end interpolate
          !
          deallocate(gcoeffs)
        end if ! N odd
        ! 
      case default
        nerr=nerr+1
        print ferrmssg,
     &     "unknown ipolmeth. Please consult help (cofima --help)" 
        return
      end select 
      !
      ! begin write ipol to disk 
      write(file1i,'(a,a)') adjustl(trim(file1)),".ipol" 
      open(51,file=trim(file1i),status="replace")
      do kdata=1,ndatai
        write(51,'(F10.6,F15.6)') rgridi(kdata),fieldi(kdata,1)       
      end do 
      close(51)
      ! end write ipol to disk 
      !
      ! end normally:
      if(talk) print fsubendext, "interpol"
      return
      !
      ! Errors:
 101  nerr=nerr+1
      print ferrmssg," when opening/reading file1"
      close(51)
      return
      !
      end subroutine
      !
!********************************************************************* 
      !
      subroutine sort(rgrid0,field0)
      use defs
      implicit none
      ! calculates the first derivative of the second column wrt the
      ! first one from finite differences 
      !
      double precision :: rgrid0(:),field0(:)
      ! internal variables
      double precision, allocatable :: rgrid(:),field(:)
      integer ndata ! dimension of input data array 
      integer idata,jdata 
      double precision rtemp,ftemp 
      !
      if(talk) print fsubstart,"sort" 
      ndata=size(rgrid0)
      if(size(field0).ne.ndata) goto 100
      allocate(rgrid(ndata),field(ndata))
      rgrid=rgrid0
      field=field0
      !
      ! begin sort grid points
      do idata=1,ndata       
        do jdata=idata+1,ndata
          if(rgrid(jdata).lt.rgrid(idata)) then
                rtemp=rgrid(idata)
                ftemp=field(idata)
                rgrid(idata)=rgrid(jdata)
                field(idata)=field(jdata)
                rgrid(jdata)=rtemp
                field(jdata)=ftemp
          end if
        end do
      end do
      rgrid0=rgrid
      field0=field
      ! end sort grid points
      !
      ! end normally:
      deallocate(rgrid,field)
      if(talk) print fsubendext, "sort"
      return
      !
      ! Errors:
 100  nerr=nerr+1
      print ferrmssg," inconsistent array dimensions"
      return
      !
      end subroutine sort
      !

        !
        end module
