        module linalg

        type point
          sequence
          integer onhull ! -1 : inside, 0 : on hull, 1 : outside
          logical dealtwith
          double precision, allocatable :: coord(:)
          logical hasfacet
          integer pointnumber
          character*1024 label
        end type point
        type facet
          sequence
          integer facetnumber
          logical visible
          integer, allocatable :: facetpoints(:)
          integer noutside
          integer, allocatable :: outsidepoints(:)
          double precision, allocatable :: thisfacet(:,:)
          character*1024, allocatable :: labels(:)
        end type facet
        type ridge
          sequence
          integer belongs_to_facets(2)
          integer, allocatable  :: ridgepoints(:)
          logical on_horizon
          double precision, allocatable :: thisridge(:,:)
        end type ridge

        contains 

!        subroutine invert_nxn_matrix(
!        !====================================================================
!        !  Computing Inverse matrix
!        !  Method: Based on the Doolittle LU method
!        !====================================================================
!        implicit none
!        integer, parameter :: n=3
!        double precision a(n,n), c(n,n)
!        integer i,j
!        ! matrix A
!          data (a(1,i), i=1,3) /  3.0,  2.0,  4.0 /
!          data (a(2,i), i=1,3) /  2.0, -3.0,  1.0 /
!          data (a(3,i), i=1,3) /  1.0,  1.0,  2.0 /
!        
!        ! print a header and the original matrix
!          write (*,200)
!          do i=1,n
!             write (*,201) (a(i,j),j=1,n)
!          end do
!        
!          call inverse(a,c,n)
!        
!        ! print the inverse matrix C = A^{-1} 
!          write (*,202)
!          do i = 1,n
!             write (*,201)  (c(i,j),j=1,n)
!          end do
! 200 format (' Computing Inverse matrix ',/,/, &
!                    ' Matrix A')
! 201 format (6f12.6)
! 202 format (/,' Inverse matrix A^{-1}')
!          end subroutine invert_nxn_matrix

          subroutine inverse(a,c,n)
        !============================================================
        ! Inverse matrix
        ! Method: Based on Doolittle LU factorization for Ax=b
        ! Alex G. December 2009
        !-----------------------------------------------------------
        ! input ...
        ! a(n,n) - array of coefficients for matrix A
        ! n      - dimension
        ! output ...
        ! c(n,n) - inverse matrix of A
        ! comments ...
        ! the original matrix a(n,n) will be destroyed 
        ! during the calculation
        !===========================================================
        implicit none 
        integer n
        double precision a(n,n), c(n,n)
        double precision L(n,n), U(n,n), b(n), d(n), x(n)
        double precision coeff
        integer i, j, k
        
        ! step 0: initialization for matrices L and U and b
        ! Fortran 90/95 aloows such operations on matrices
        L=0.0
        U=0.0
        b=0.0
        
        ! step 1: forward elimination
        do k=1, n-1
           do i=k+1,n
              coeff=a(i,k)/a(k,k)
              L(i,k) = coeff
              do j=k+1,n
                 a(i,j) = a(i,j)-coeff*a(k,j)
              end do
           end do
        end do
        
        ! Step 2: prepare L and U matrices 
        ! L matrix is a matrix of the elimination coefficient
        ! + the diagonal elements are 1.0
        do i=1,n
          L(i,i) = 1.0
        end do
        ! U matrix is the upper triangular part of A
        do j=1,n
          do i=1,j
            U(i,j) = a(i,j)
          end do
        end do
        
        ! Step 3: compute columns of the inverse matrix C
        do k=1,n
          b(k)=1.0
          d(1) = b(1)
        ! Step 3a: Solve Ld=b using the forward substitution
          do i=2,n
            d(i)=b(i)
            do j=1,i-1
              d(i) = d(i) - L(i,j)*d(j)
            end do
          end do
        ! Step 3b: Solve Ux=d using the back substitution
          x(n)=d(n)/U(n,n)
          do i = n-1,1,-1
            x(i) = d(i)
            do j=n,i+1,-1
              x(i)=x(i)-U(i,j)*x(j)
            end do
            x(i) = x(i)/u(i,i)
          end do
        ! Step 3c: fill the solutions x(n) into column k of C
          do i=1,n
            c(i,k) = x(i)
          end do
          b(k)=0.0
        end do
        end subroutine inverse

c---------------------------------------------------------------------

        subroutine rotate_matrix(n,A,S,B)
        ! 
        ! rotate matrix A to new coordinate system: B = S⁺ A S
        !
        use defs
        implicit none 
        integer n,j
        double precision A(n,n), B(n,n), S(n,n), S_trans(n,n)
        ! 
        if (talk) then
          print fsubstart, 'rotate_matrix'
          ! print matrix
          print '(8x," ")'
          print '(8x,"matrix dimension: ",I0)',n
          print '(8x," ")'
          print '(8x,"Matrix:")'
          do j=1,n
            print '(8x,3(F10.6))',A(j,1:3)
          end do  
          print '(8x," ")'
          ! print rotation matrix
          print '(8x," ")'
          print '(8x,"Rotation matrix:")'
          do j=1,n
            print '(8x,3(F10.6))',S(j,1:3)
          end do  
          print '(8x," ")'
        end if ! talk
        !
        S_trans=transpose(S)
        !
        if (talk) then
          print '(8x," ")'
          print '(8x,"Rotation matrix transposed :")'
          do j=1,n
            print '(8x,3(F10.6))',S_trans(j,1:3)
          end do  
          print '(8x," ")'
        end if ! talk

        B=matmul(A,S)
        B=matmul(S_trans,B)
        !
        if (talk) then
          print '(8x," ")'
          print '(8x,"Rotated matrix:")'
          do j=1,n
            print '(8x,3(F10.6))',B(j,1:3)
          end do  
          print '(8x," ")'
        end if ! talk
        !
        if (talk) then      
          print fsubendext, 'rotate_matrix'
        end if  
        !  
        end subroutine rotate_matrix

c---------------------------------------------------------------------
      
        subroutine quickhull(filename,d)
        ! determines the convex hull of a set of points in 3d using the 
        ! algorithm of Barber, Barber, Dobkin, Huhdanpaa, ACM
        ! Transactions on Mathematical software volume 22 no 4, pages
        ! 469-483, 1996
        use defs
        use combi, only : facult
        implicit none
        character(len=*), intent (in) :: filename
        integer, intent(in) :: d ! dimension
        ! local variables
        integer npoints,nhullpoints,nfacets
        integer ipoint,jpoint,ifarthestpoint
        integer i,ifacet,jfacet,nvisible
        integer i1,i2,i3
        character line*256,pointformat*256
        double precision distance,distance2,distance3
        double precision, allocatable :: vec(:),vecpar(:),vecperp(:)
        double precision, allocatable :: vecperp2(:),vec12(:),vec13(:)
        integer, allocatable :: newfacetpoints(:)
        type(point), allocatable :: points(:)
        type(point), allocatable :: hullpoints(:)
        type(facet), allocatable :: facets(:)
        type(facet), allocatable :: facets_temp(:)
        type(ridge), allocatable :: ridges(:)
        integer nridges,iridge
        double precision, allocatable :: midpoint(:)
        integer, allocatable :: outsidepoints_temp(:)
        type(facet), allocatable :: V(:),V_temp(:)
        integer ntotoutside
        integer nrp,nfp
        integer idum
        double precision, allocatable :: crosspoint(:)
        !double precision, parameter :: tol=1.0d-8
        double precision, parameter :: tol=1.0d-12
        !
        if (talk) print fsubstart,"quickhull"
        if (talk) print'(8x,"convex hull in ",I0," dimensions")',d
        !
        ! begin check if dimension is implemented:
        !if (d.le.1) call error_stop('dim too small')
        !if (d.gt.3) call error_stop('dim not implemented')
        ! end check if dimension is implemented:
        !
        !
        ! begin read points from file
        !
        call read_points()
        !
        ! end read points from file
        !
        !
        ! if there are only d+1 points or less, we are done:
        if (npoints.le.d+1) then
          nhullpoints=npoints
          allocate(hullpoints(nhullpoints))
          do ipoint=1,nhullpoints
            allocate(hullpoints(ipoint)%coord(1:d))
            hullpoints(ipoint)%coord(:)=points(ipoint)%coord(:)
            hullpoints(ipoint)%pointnumber=ipoint
            hullpoints(ipoint)%label=points(ipoint)%label
          end do
          !call get_initial_facets()
          !call get_initial_facets_2()
          call get_facets(hullpoints,facets)
          nfacets=size(facets)
          goto 50
        end if 
        !
        !
        ! begin select initial convex (d+1)hedron "simplex"
        !
        !call get_simplex()
        call get_simplex_2(points,hullpoints)
        ! 
        ! begin write simplex
        nhullpoints=size(hullpoints)
        open(51,file='SIMPLEX_POINTS',status='replace')
        pointformat=' '
        write(pointformat,'("(",I0,"(F10.6),1x,A)")') d
        do ipoint=1,nhullpoints
          write(51,pointformat) hullpoints(ipoint)%coord(1:d),            &
     &          adjustl(trim(hullpoints(ipoint)%label))
        end do ! ipoint
        close(51)
        !
        ! end write simplex
        !
        ! end select initial convex (d+1)hedron "simplex"
        !
        !
        ! begin get mid point
        !
        allocate(midpoint(d))
        midpoint=0.0d0
        do ipoint=1,nhullpoints
          midpoint=midpoint+points(hullpoints(ipoint)%pointnumber)%coord
        end do ! ipoint
        midpoint=midpoint/dble(nhullpoints)
        !
        ! end get mid point
        !
        !
        !
        ! begin get initial facets
        !
        !call get_initial_facets()
        !call get_initial_facets_2()
        call get_facets(hullpoints,facets)
        nfacets=size(facets)
        !
        ! begin write simplex facets
        !
        open(51,file='SIMPLEX_FACETS',status='replace')
        pointformat=' '
!        write(pointformat,'("(",I0,"(",I0,"(F10.6),3x),",I0,              &
!     &          "(I0,x)",")")') size(facets(1)%facetpoints),d,            &
!     &                          size(facets(1)%facetpoints) 
        write(pointformat,'("(",I0,"(",I0,"(F10.6),3x),",I0,              &
     &        "(I0,x),",I0,"(A,x)",")")') size(facets(1)%facetpoints),d,  &
     &          size(facets(1)%facetpoints),size(facets(1)%facetpoints) 
        do ifacet=1,nfacets
          write(51,pointformat) (points(facets(ifacet)                    &
     &    %facetpoints(i))%coord(1:d),i=1,size(facets(1)%facetpoints)),   &
     &    facets(ifacet)%facetpoints(1:size(facets(1)%facetpoints)),      &
     &    (adjustl(trim(points(facets(ifacet)%facetpoints(i))%label)),    &
     &     i=1,size(facets(1)%facetpoints))
        end do ! ifacet
        close(51)
        !
        ! end write simplex facets
        !
        ! end get initial facets
        !
        !
        !
        ! begin get points outside facets
        !
        call get_outside_points(facets,points,midpoint)
        ntotoutside=sum(facets(:)%noutside)
        !
        ! begin DEBUG
        !
        open(52,file="DEBUG.OUT",status="replace")
        write(52,'("number and index of points outside each current face  &
     &t:")')
        do ifacet=1,nfacets
          do ipoint=1,facets(ifacet)%noutside
            write(52,*) facets(ifacet)%outsidepoints(ipoint)
          end do ! ipoint
        end do ! ifacet
        !
        ! end DEBUG
        !
        !
        ! end get points outside facets
        !
        !
        !
        !
        ! begin rebuild facets and hull
        !
        do while (ntotoutside>0)
          write(52,'(" ")') ! DEBUG
          write(52,'(8x,"ntotoutside=",I0)') ntotoutside ! DEBUG
          !
          !do ifacet=1,nfacets
          nfacets=size(facets)
          ifacet=1 
          do while (ifacet.le.nfacets)
            write(52,'(" ")') ! DEBUG
            write(52,'(8x,"ifacet=",I0)') ifacet ! DEBUG
            write(52,'(8x,"facetpoints:")') ! DEBUG
            write(52,*) facets(ifacet)%facetpoints(:) ! DEBUG
            write(52,'(8x,"noutside(ifacet)=",I0)')                       &
     &          facets(ifacet)%noutside ! DEBUG
            !
            if (facets(ifacet)%noutside>0) then
              !
              ! begin find farthest point p in outside set of F
              !
              write(52,'(8x,"outsidepoints(ifacet)=",10(I4))')            &
     &              facets(ifacet)%outsidepoints ! DEBUG
              ifarthestpoint=facets(ifacet)%outsidepoints(1)
!              distance2=distance_fp(facets(ifacet),                       &
!     &        points(facets(ifacet)%outsidepoints(1)))
              call distance_from_surface(                                 &
     &             points(facets(ifacet)%outsidepoints(1))%coord,         &
     &             facets(ifacet)%thisfacet,crosspoint,distance2)
              !
              do ipoint=2,facets(ifacet)%noutside
                ! 
!                if (distance_fp(facets(ifacet),                           &
!     &            points(facets(ifacet)%outsidepoints(ipoint))).gt.       &
!     &            distance2) then 
!                  distance2=distance_fp(facets(ifacet),                   &
!     &            points(facets(ifacet)%outsidepoints(ipoint)))    
!                    ifarthestpoint=facets(ifacet)%outsidepoints(ipoint)
!                end if
                call distance_from_surface(                               &
     &             points(facets(ifacet)%outsidepoints(ipoint))%coord,    &
     &             facets(ifacet)%thisfacet,crosspoint,distance3)
                if (distance3>distance2) then
                  distance2=distance3
                  ifarthestpoint=facets(ifacet)%outsidepoints(ipoint)
                end if
                write(52,'(8x,"ifarthestpoint(ifacet)=",I0)')             &
     &          ifarthestpoint ! DEBUG
              !
              end do ! ipoint
              !print*, "adding point ",ifarthestpoint
              !
              ! end find farthest point p in outside set of F
              !
              !
              !
              ! begin find facets V visible by p
              !
              ! initialize V:
              if(allocated(V)) deallocate(V)
              allocate(V(1))
              write(52,'(8x,"setting visibility for facet ifacet=",I0)')  &
     &             ifacet ! DEBUG
              V(1)=facets(ifacet)
              if(allocated(V_temp)) deallocate(V_temp)
              allocate(V_temp(1))
              V_temp(1)=facets(ifacet)
              nvisible=1
              facets(:)%visible=.false.
              facets(ifacet)%visible=.true.
              ! check other facets: (TODO: only unvisited neighbors)
              nfp=size(facets(ifacet)%facetpoints)
              !do jfacet=1,nfacets ! DEBUG
              do jfacet=1,size(facets) ! DEBUG
                if (jfacet.ne.ifacet) then
!                  print*,"ifacet,jfacet=",ifacet,jfacet
!                  do i17=1,size(facets(jfacet)%facetpoints)
!                  print*,"jfacet,ifp,facetpoints=",                       &
!     &                jfacet,i17,facets(jfacet)%thisfacet(i17,:)
!                  end do
!                  if (is_visible(facets(jfacet)%thisfacet,                &
!     &                points(ifarthestpoint)%coord,midpoint)) then
!                  if (is_visible_3(facets(jfacet)%thisfacet,              &
!     &                points(ifarthestpoint)%coord,midpoint)) then
                  if (is_visible_4(facets(jfacet)%thisfacet,              &
     &               points(ifarthestpoint)%coord,midpoint,-tol))then
                    write(52,'(8x,"point is also visible for facet",I0)   &
     &')               jfacet ! DEBUG
                    write(52,'(8x,"facet ",I0," has the points ",10(I4))  &
     &')               jfacet,facets(jfacet)%facetpoints ! DEBUG
                    nvisible=nvisible+1
                    facets(jfacet)%visible=.true.
                    deallocate(V_temp)
                    allocate(V_temp(nvisible))
                    V_temp(1:nvisible-1)=V(:)
                    V_temp(nvisible)=facets(jfacet)
                    deallocate(V)
                    allocate(V(nvisible))
                    V=V_temp
                  end if ! is_outside
                  !print*, jfacet,facets(jfacet)%visible
                end if ! jfacet.ne.ifacet
              end do ! jfacet
              !print*,"V (visible by point)::"
              !do i=1,nvisible
              !  print*, V(i)%facetpoints 
              !end do
              !
              ! end find facets V visible by p
              !
              !
              !
              ! begin find horizon of V
              !
              !call get_ridges(V,ridges)
              call get_all_ridges(V,ridges)
              !call get_ridges(facets_new,ridges)
              ridges(:)%on_horizon=.false.
              nridges=size(ridges)
              !print*,"nridges=",nridges
!              do iridge=1,nridges
!                print*,"iridge,ridgepoint:",iridge,                       &
!     &             ridges(iridge)%ridgepoints
!                print*,ridges(iridge)%ridgepoints," belongs_to_facet ",   &
!     &          facets(ridges(iridge)%belongs_to_facets(1))%facetpoints,  &
!     &          facets(ridges(iridge)%belongs_to_facets(2))%facetpoints
!                print*,ridges(iridge)%ridgepoints," belongs_to_facet ",   &
!     &                (ridges(iridge)%belongs_to_facets(1)) ,             &
!     &                (ridges(iridge)%belongs_to_facets(2)) 
!              end do
              write(52,'("horizon (ridges):")') ! DEBUG
              do iridge=1,nridges
!                if (facets(ridges(iridge)%belongs_to_facets(1))%visible   &
!     &             .neqv.facets(ridges(iridge)%belongs_to_facets(2))      &
!     &              %visible) then
                if (any(ridges(iridge)%belongs_to_facets(:)==0)) then 
                     ridges(iridge)%on_horizon=.true.
                     !print*,ridges(iridge)%ridgepoints
                end if
              end do ! iridge
              !
              ! end find horizon of V
              !
              !
              !
              ! begin add new facets from horizon to p
              !
              if (.not.allocated(newfacetpoints))                         &
     &            allocate(newfacetpoints(d))
              do iridge=1,nridges
                if (ridges(iridge)%on_horizon) then
                  newfacetpoints=(/ridges(iridge)%ridgepoints,            &
     &              ifarthestpoint/)
                  call add_to_facets(facets,newfacetpoints,points,        &
     &              facets_temp)
                  deallocate(facets)
                  allocate(facets(size(facets_temp)))
                  facets=facets_temp
                  !
                  ! BEGIN NEW
                  !
                  ! begin find points outside new facets
                  !
                   call get_outside_points(facets(size(facets):           &
     &               size(facets)),points,midpoint)
                  !    
                  ! end find points outside new facets 
                  !
                  ! END NEW
                  !
                  deallocate(facets_temp)
                end if ! on_horizon
              end do ! iridge
              points(ifarthestpoint)%dealtwith=.true.
              write(52,'("old facets plus new facets:")') ! DEBUG
              do idum=1,size(facets) ! DEBUG
                write(52,*) facets(idum)%facetpoints(:) ! DEBUG
              end do  
              !
              ! end add new facets from horizon to p
              !
              !
              ! BEGIN OLD
              !
!              ! begin find points outside new facets
!              !
!              call get_outside_points(facets(nfacets+1:                   &
!     &               size(facets)),points,midpoint)
!              !    
!              ! end find points outside new facets 
              !
              ! END OLD
              !
              !
              ! begin remove inner facets
              !
              do i=1,nvisible
                call remove_from_facets(facets,V(i),facets_temp)
                deallocate(facets)
                allocate(facets(size(facets_temp)))
                facets=facets_temp
                deallocate(facets_temp)
              end do ! i
              !
              write(52,'(8x,"updated facets (inner facets removed):")') ! DEBUG
              do idum=1,size(facets) ! DEBUG
                write(52,*) facets(idum)%facetpoints(:) ! DEBUG
              end do ! DEBUG  
!              !
              ! end remove inner facets
              !
              !
            end if ! while facets(ifacet)%noutside>0
            !
            ifacet=ifacet+1
            !
          end do ! ifacet
          !
          ! BEGIN OLD
          !
!          call get_outside_points(facets,points,midpoint) ! DEBUG 
          ntotoutside=sum(facets(:)%noutside)
          nfacets=size(facets)
          !
          ! END OLD
          !
        end do ! while ntotoutside.gt.0
        close(52) ! DEBUG file
        !
        ! end rebuild facets and hull
        !  
        !  
        !  
 50     continue
        ! 
        ! begin write convex hull points
        !  
        if(allocated(hullpoints)) deallocate(hullpoints)
        call get_points(facets,hullpoints)
        nhullpoints=size(hullpoints)
        open(51,file='HULL_POINTS',status='replace')
        pointformat=' '
        !write(pointformat,'("(",I0,"(F10.6))")') d
        write(pointformat,'("(",I0,"(F20.12),2x,I0,1x,A)")') d
        !print*,pointformat
        do ipoint=1,nhullpoints
          !print*,hullpoints(ipoint)%pointnumber
          hullpoints(ipoint)%label=                                       &
     &       points(hullpoints(ipoint)%pointnumber)%label
          write(51,pointformat) hullpoints(ipoint)%coord(1:d),            &
     &           hullpoints(ipoint)%pointnumber,                          &
     &           adjustl(trim(hullpoints(ipoint)%label))
        end do ! ipoint
        close(51)
        !
        ! end write convex hull points
        !
        !
        !
        ! begin write convex hull facets
        !
        open(51,file='HULL_FACETS',status='replace')
        pointformat=' '
!        write(pointformat,'("(",I0,"(",I0,"(F20.12),3x),",I0,             &
!     &          "(I0,x)",")")') size(facets(1)%facetpoints),d,            &
!     &                          size(facets(1)%facetpoints) 
        write(pointformat,'("(",I0,"(",I0,"(F20.12),3x),",I0,             &
     &        "(I0,x),",I0,"(A,x)",")")') size(facets(1)%facetpoints),d,  &
     &          size(facets(1)%facetpoints),size(facets(1)%facetpoints) 
        !print*,pointformat
        do ifacet=1,nfacets
          write(51,pointformat) (points(facets(ifacet)                    &
     &    %facetpoints(i))%coord(1:d),i=1,size(facets(1)%facetpoints)),   &
     &       facets(ifacet)%facetpoints(1:size(facets(1)%facetpoints)),   &
     &    (adjustl(trim(points(facets(ifacet)%facetpoints(i))%label)),    &
     &     i=1,size(facets(1)%facetpoints))
        end do ! ifacet
        close(51)
        !
        ! end write convex hull facets
        !
        !
        !
        ! begin write convex hull ridges
        !
        if (allocated(ridges)) deallocate(ridges)
        !call get_ridges(facets,ridges)
        call get_all_ridges(facets,ridges)
        !print*, "ridges obtained for writing"
        nridges=size(ridges)
        open(51,file='HULL_RIDGES',status='replace')
        pointformat=' '
        !write(pointformat,'("(",I0,"(",I0,"(F10.6),3x))")') d-1,d
        nrp=size(ridges(1)%ridgepoints)
        write(pointformat,'("(",I0,"(",I0,"(F20.12),3x))")') nrp,d
        do iridge=1,nridges
          write(51,pointformat)                                           &
     &          (points(ridges(iridge)%ridgepoints(i))                    &
     &         %coord(1:d),i=1,nrp)
        end do ! iridge
        close(51)
        !
        ! end write convex hull ridges
        !
        !
        !
        if (talk) print fsubendext,"quickhull"
        return
        ! 
        !
c---------------------------------------------------------------------

        contains

c---------------------------------------------------------------------

        subroutine read_points()
        implicit none
        type(point), allocatable :: points_temp(:)
        !!double precision, parameter :: tol=1.0E-12
        !double precision, parameter :: tol=1.0D-8
        !double precision, parameter :: tol=1.0D-12
        !double precision, parameter :: tol=5.0D-3 ! absolute tolerance to distinguish two points
        !double precision, parameter :: tol=1.0D-2 ! absolute tolerance to distinguish two points
        double precision, parameter :: tol=1.0D-6 ! absolute tolerance to distinguish two points
        double precision :: maxabspoint ! maximum distance between origin and point
        double precision :: tolrel ! relative tolerance to distinguish two points, tolrel=maxabspoint*tol
        integer ipoint,jpoint,kpoint
        logical newpoint
        double precision dist
        !
        !
        open(51,file=filename,status='old',err=101) 
        npoints=0
 10     read(51,'(A256)', err=102,end=12) line  
        line=trim(adjustl(line))
        if(len_trim(line).gt.0.and.line(1:1).ne.'#'                       &
     &     .and.line(1:1).ne.'!'.and.line(1:1).ne.'%') npoints=npoints+1
        goto 10
 12     continue
        !
!        if (npoints.lt.d+1) then
!          close(51)
!          call error_stop('too few points')
!        end if
        !
        allocate(points(npoints))
        do ipoint=1,npoints
          allocate(points(ipoint)%coord(1:d))
        end do
        ipoint=0
        maxabspoint=0.0d0
        rewind(51)
 14     read(51,'(A256)', err=102,end=16) line  
        line=trim(adjustl(line))
        if(len_trim(line).gt.0.and.line(1:1).ne.'#'                       &
     &     .and.line(1:1).ne.'!'.and.line(1:1).ne.'%') then
          ipoint=ipoint+1
          read(line,*) points(ipoint)%coord(1:d),points(ipoint)%label
!          print '(8x,A)',adjustl(trim(points(ipoint)%label))
          points(ipoint)%dealtwith=.false.
          points(ipoint)%onhull=1
          points(ipoint)%hasfacet=.false.
          if (norm2(points(ipoint)%coord).gt.maxabspoint)                 &
     &      maxabspoint=norm2(points(ipoint)%coord)
        end if
        goto 14
 16     continue
        close(51)
        !
        !
        tolrel=tol*maxabspoint 
        !print '(8x,"tol=",F20.16)',tol
        !print '(8x,"maxabspoint=",F20.16)',maxabspoint
        print '(8x,"tolrel=",F20.16)',tolrel
        !
        ! begin remove double points
        !
        !print*,"initial npoints:",npoints
        allocate(points_temp(npoints))
        open(51,file="POINTS.DISCARDED",status="replace")
        ipoint=1
        do kpoint=1,npoints
          newpoint=.true.
          do jpoint=1,ipoint-1
            dist=norm2(points(kpoint)%coord-points(jpoint)%coord)
            if (dist.lt.tolrel)  newpoint=.false.
          end do ! jpoint
          if (newpoint) then
            points_temp(ipoint)=points(kpoint)
            ipoint=ipoint+1
          else
            write(51,*) points(kpoint)%coord
          end if
        end do ! kpoint
        npoints=ipoint-1
        !print*,"final npoints:",npoints
        deallocate(points)
        allocate(points(npoints))
        points(1:npoints)=points_temp(1:npoints)
        deallocate(points_temp)
        close(51)
        !
        ! end remove double points
        !
        !
        !
        return

 101    close(51)
        call error_stop('file not found')
 102    close(51)
        call error_stop('file empty?')

      end subroutine read_points

c---------------------------------------------------------------------
!
!        subroutine get_simplex()        
!        implicit none
!        double precision :: norm2
!        !
!        ! begin select initial convex (d+1)hedron "simplex"
!        !
!        nhullpoints=d+1
!        allocate(hullpoints(nhullpoints))
!        do ipoint=1,nhullpoints
!          allocate(hullpoints(ipoint)%coord(1:d))
!        end do
!        distance=0.0d0
!        allocate(vec(1:d))
!        !
!        !
!        !
!        ! begin find first 2 points (the ones farthest apart)
!        !
!        do ipoint=1,npoints
!          do jpoint=1,npoints
!            vec=points(ipoint)%coord(:)-points(jpoint)%coord(:)
!            if (norm2(vec).gt.distance) then
!              distance=norm2(vec)
!              hullpoints(1)%coord(1:d)=points(ipoint)%coord(1:d)
!              hullpoints(1)%pointnumber=ipoint
!              hullpoints(2)%coord(1:d)=points(jpoint)%coord(1:d)
!              hullpoints(2)%pointnumber=jpoint
!            end if
!          end do ! jpoint
!        end do ! ipoint
!        !
!        ! end find first 2 points (the ones farthest apart)
!        !
!        !
!        !
!        ! begin find third point (largest distance from connection
!        ! between first two points)
!        !
!        distance2=0.0d0
!        allocate(vec12(1:d),vecpar(1:d),vecperp(1:d))
!        vec12=hullpoints(2)%coord(:)-hullpoints(1)%coord(:)
!        !
!        do ipoint=1,npoints
!          !
!          vec=points(ipoint)%coord(:)-hullpoints(1)%coord(:)
!          vecpar=vec12*dot_product(vec,vec12)/norm2(vec12)**2        
!          vecperp=vec-vecpar
!          !
!          if (norm2(vecperp).gt.distance2) then
!            distance2=norm2(vecperp)
!            hullpoints(3)%coord(:)=points(ipoint)%coord(:)
!            hullpoints(3)%pointnumber=ipoint
!          end if ! norm2(vecperp).gt.distance2
!          !
!        end do ! ipoint
!        !
!        !
!        ! end find third point (largest distance from connection
!        ! between first two points)
!        !
!        !
!        ! begin if d>=3, get fourth point (largest distance from base
!        ! triangle)
!        !
!        if (d>=3) then
!          !
!          distance3=0.0d0
!          allocate(vec13(1:d),vecperp2(1:d))
!          vec13=hullpoints(3)%coord(:)-hullpoints(1)%coord(:)
!          !
!          do ipoint=1,npoints
!            !
!            vec=points(ipoint)%coord(:)-hullpoints(1)%coord(:)
!            vecpar=vec12*dot_product(vec,vec12)/norm2(vec12)**2   
!            vecperp=vec13-dot_product(vec13,vec12)*vec12                  &
!     &              /norm2(vec12)**2
!              vecperp2=vec-vecpar-vecperp*dot_product(vec,vecperp)        &
!     &               /norm2(vecperp)**2
!            !
!            if (norm2(vecperp2).gt.distance3) then
!              !
!              distance3=norm2(vecperp2)
!              hullpoints(4)%coord(:)=points(ipoint)%coord(:)
!              hullpoints(4)%pointnumber=ipoint
!              !
!            end if ! norm2(distvecperp2).gt.distance3
!            !
!          end do ! ipoint
!          !
!        end if ! d>=3
!        !
!        ! end if d>=3, get fourth point (largest distance from base
!        ! triangle)
!        !
!        !
!        !
!        ! begin set hullpoints to "dealtwith"
!        !
!        do ipoint=1,nhullpoints
!          !
!          points(hullpoints(ipoint)%pointnumber)%dealtwith=.true.
!          points(hullpoints(ipoint)%pointnumber)%onhull=0
!          ! 
!        end do
!        !
!        ! end set hullpoints to "dealtwith"
!
!        end subroutine get_simplex
!        
c---------------------------------------------------------------------
!
!        subroutine get_initial_facets()
!        implicit none
!
!        nfacets=d+1
!        do i=0,d-2
!          nfacets=nfacets*(d-i)
!        end do
!        nfacets=nfacets/facult(d)
!        allocate(facets(nfacets))
!        do ifacet=1,nfacets
!          allocate(facets(ifacet)%facetpoints(1:d))
!          facets(ifacet)%visible=.false.
!          allocate(facets(ifacet)%thisfacet(d,d))
!        end do
!        !
!        select case(d)
!        case(2)
!          ifacet=1
!          do i1=1,nhullpoints
!            do i2=i1+1,nhullpoints
!              facets(ifacet)%facetpoints(1:2)=                            &
!     &                   (/hullpoints(i1)%pointnumber,                    &
!     &                     hullpoints(i2)%pointnumber/)
!              ifacet=ifacet+1
!            end do
!          end do
!        case(3)
!          ifacet=1
!          do i1=1,nhullpoints
!            do i2=i1+1,nhullpoints
!              do i3=i2+1,nhullpoints
!                facets(ifacet)%facetpoints(1:3)=                          &
!     &                 (/hullpoints(i1)%pointnumber,                      &
!     &                   hullpoints(i2)%pointnumber,                      &
!     &                   hullpoints(i3)%pointnumber/)
!                ifacet=ifacet+1
!              end do
!            end do
!          end do
!        case default
!        end select
!
!        do ifacet=1,nfacets
!          facets(ifacet)%facetnumber=ifacet
!          do ipoint=1,d
!            facets(ifacet)%thisfacet(ipoint,1:d)=                         &
!     &           points(facets(ifacet)%facetpoints(ipoint))%coord(1:d)
!          end do ! ipoint
!          facets(ifacet)%noutside=0
!        end do ! ifacet
!
!        end subroutine get_initial_facets

c---------------------------------------------------------------------
!
!        subroutine get_initial_facets_2()
!        use combi, only : n_over_k
!        implicit none
!        !
!        ! local :
!        integer nhp,nfp
!        integer, allocatable :: indices(:)
!        integer i,i_ind,j_ind,k_ind,ifacet
!        logical newfacet
!
!        nhp=size(hullpoints)
!        nfp=min(d,nhp)
!        nfacets=n_over_k(nhp,nfp)
!        !
!        !
!        allocate(indices(nfp))
!        do i_ind=1,size(indices)
!            !indices(i_ind)=i_ind
!            indices(i_ind)=1
!        end do ! i_ind
!        i_ind=size(indices)
!
!        allocate(facets(nfacets))
!        do ifacet=1,nfacets
!          allocate(facets(ifacet)%facetpoints(1:nfp))
!          facets(ifacet)%visible=.false.
!          allocate(facets(ifacet)%thisfacet(nfp,d))
!        end do
!        !
!        ifacet=1
!        do while (ifacet<=nfacets)
!            newfacet=.true.
!            do j_ind=1,size(indices)
!              do k_ind=j_ind+1,size(indices)
!                if (indices(j_ind)>=indices(k_ind)) newfacet=.false.
!              end do ! k_ind
!            end do ! j_ind
!            if (newfacet) then
!              facets(ifacet)%facetpoints(1:nfp)=                           &
!     &          hullpoints(indices(1:nfp))%pointnumber
!              ifacet=ifacet+1
!            end if
!            do while (all(indices(i_ind:nfp)==nhp.and.i_ind>1))
!                i_ind=i_ind-1
!            end do
!            do j_ind=i_ind+1,size(indices)
!                indices(j_ind)=1
!            end do
!            indices(i_ind)=indices(i_ind)+1
!            i_ind=size(indices)
!        end do ! ifacet <= nfacets
!        !
!        !
!        do ifacet=1,nfacets
!          facets(ifacet)%facetnumber=ifacet
!          do ipoint=1,nfp
!            facets(ifacet)%thisfacet(ipoint,1:d)=                         &
!     &           points(facets(ifacet)%facetpoints(ipoint))%coord(1:d)
!          end do ! ipoint
!          facets(ifacet)%noutside=0
!        end do ! ifacet
!        !
!        end subroutine get_initial_facets_2
!
c---------------------------------------------------------------------

        end subroutine quickhull

c---------------------------------------------------------------------
      
        subroutine quickhull2(filename,d)
        ! determines the convex hull of a set of points in 3d using the 
        ! algorithm of Barber, Barber, Dobkin, Huhdanpaa, ACM
        ! Transactions on Mathematical software volume 22 no 4, pages
        ! 469-483, 1996
        ! Like quickhull, but hopefully optimized
        !
        use defs
        use combi, only : facult
        implicit none
        character(len=*), intent (in) :: filename
        integer, intent(in) :: d ! dimension
        ! local variables
        integer npoints,nhullpoints,nfacets
        integer ipoint,jpoint,ifarthestpoint
        integer i,ifacet,jfacet,nvisible
        integer i1,i2,i3
        character line*256,pointformat*256
        double precision distance,distance2,distance3
        double precision, allocatable :: vec(:),vecpar(:),vecperp(:)
        double precision, allocatable :: vecperp2(:),vec12(:),vec13(:)
        integer, allocatable :: newfacetpoints(:)
        type(point), allocatable :: points(:)
        type(point), allocatable :: hullpoints(:)
        type(facet), allocatable :: facets(:)
        type(facet), allocatable :: facets_temp(:)
        type(ridge), allocatable :: ridges(:)
        integer nridges,iridge
        double precision, allocatable :: midpoint(:)
        integer, allocatable :: outsidepoints_temp(:)
        type(facet), allocatable :: V(:),V_temp(:)
        integer ntotoutside
        integer nrp,nfp
        integer idum
        double precision, allocatable :: crosspoint(:)
        !double precision, parameter :: tol=1.0d-8
        double precision, parameter :: tol=1.0d-12
        !
        if (talk) print fsubstart,"quickhull2"
        if (talk) print'(8x,"convex hull in ",I0," dimensions")',d
        !
        !
        ! begin read points from file
        !
        call read_points()
        !
        ! end read points from file
        !
        !
        ! if there are only d+1 points or less, we are done:
        if (npoints.le.d+1) then
          nhullpoints=npoints
          allocate(hullpoints(nhullpoints))
          do ipoint=1,nhullpoints
            allocate(hullpoints(ipoint)%coord(1:d))
            hullpoints(ipoint)%coord(:)=points(ipoint)%coord(:)
            hullpoints(ipoint)%pointnumber=ipoint
            hullpoints(ipoint)%label=points(ipoint)%label
          end do
          call get_facets(hullpoints,facets)
          nfacets=size(facets)
          goto 50
        end if 
        !
        !
        ! begin select initial convex (d+1)hedron "simplex"
        !
        !call get_simplex()
        call get_simplex_2(points,hullpoints)
        ! 
        ! begin write simplex
        nhullpoints=size(hullpoints)
        open(51,file='SIMPLEX_POINTS',status='replace')
        pointformat=' '
        write(pointformat,'("(",I0,"(F10.6),1x,A)")') d
        do ipoint=1,nhullpoints
          write(51,pointformat) hullpoints(ipoint)%coord(1:d),            &
     &          adjustl(trim(hullpoints(ipoint)%label))
        end do ! ipoint
        close(51)
        !
        ! end write simplex
        !
        ! end select initial convex (d+1)hedron "simplex"
        !
        !
        ! begin get mid point
        !
        allocate(midpoint(d))
        midpoint=0.0d0
        do ipoint=1,nhullpoints
          midpoint=midpoint+points(hullpoints(ipoint)%pointnumber)%coord
        end do ! ipoint
        midpoint=midpoint/dble(nhullpoints)
        !
        ! end get mid point
        !
        !
        !
        ! begin get initial facets
        !
        !call get_initial_facets()
        !call get_initial_facets_2()
        call get_facets(hullpoints,facets)
        nfacets=size(facets)
        !
        ! begin write simplex facets
        !
        open(51,file='SIMPLEX_FACETS',status='replace')
        pointformat=' '
        write(pointformat,'("(",I0,"(",I0,"(F10.6),3x),",I0,              &
     &        "(I0,x),",I0,"(A,x)",")")') size(facets(1)%facetpoints),d,  &
     &          size(facets(1)%facetpoints),size(facets(1)%facetpoints) 
        do ifacet=1,nfacets
          write(51,pointformat) (points(facets(ifacet)                    &
     &    %facetpoints(i))%coord(1:d),i=1,size(facets(1)%facetpoints)),   &
     &    facets(ifacet)%facetpoints(1:size(facets(1)%facetpoints)),      &
     &    (adjustl(trim(points(facets(ifacet)%facetpoints(i))%label)),    &
     &     i=1,size(facets(1)%facetpoints))
        end do ! ifacet
        close(51)
        !
        ! end write simplex facets
        !
        ! end get initial facets
        !
        !
        !
        ! begin get points outside facets
        !
        call get_outside_points(facets,points,midpoint)
        ntotoutside=sum(facets(:)%noutside)
        !
        ! end get points outside facets
        !
        !
        !
        !
        ! begin rebuild facets and hull
        !
        do while (ntotoutside>0)
          !
          !do ifacet=1,nfacets
          nfacets=size(facets)
          ifacet=1 
          do while (ifacet.le.nfacets)
            !
            if (facets(ifacet)%noutside>0) then
              !
              ! begin find farthest point p in outside set of F
              !
              ifarthestpoint=facets(ifacet)%outsidepoints(1)
!              distance2=distance_fp(facets(ifacet),                       &
!     &        points(facets(ifacet)%outsidepoints(1)))
              call distance_from_surface(                                 &
     &             points(facets(ifacet)%outsidepoints(1))%coord,         &
     &             facets(ifacet)%thisfacet,crosspoint,distance2)
              !
              do ipoint=2,facets(ifacet)%noutside
                ! 
!                if (distance_fp(facets(ifacet),                           &
!     &            points(facets(ifacet)%outsidepoints(ipoint))).gt.       &
!     &            distance2) then 
!                  distance2=distance_fp(facets(ifacet),                   &
!     &            points(facets(ifacet)%outsidepoints(ipoint)))    
!                    ifarthestpoint=facets(ifacet)%outsidepoints(ipoint)
!                end if
                call distance_from_surface(                               &
     &             points(facets(ifacet)%outsidepoints(ipoint))%coord,    &
     &             facets(ifacet)%thisfacet,crosspoint,distance3)
                if (distance3>distance2) then
                  distance2=distance3
                  ifarthestpoint=facets(ifacet)%outsidepoints(ipoint)
                end if
              !
              end do ! ipoint
              !print*, "adding point ",ifarthestpoint
              !
              ! end find farthest point p in outside set of F
              !
              !
              !
              ! begin find facets V visible by p
              !
              ! initialize V:
              if(allocated(V)) deallocate(V)
              allocate(V(1))
              V(1)=facets(ifacet)
              if(allocated(V_temp)) deallocate(V_temp)
              allocate(V_temp(1))
              V_temp(1)=facets(ifacet)
              nvisible=1
              facets(:)%visible=.false.
              facets(ifacet)%visible=.true.
              ! check other facets: (TODO: only unvisited neighbors)
              nfp=size(facets(ifacet)%facetpoints)
              !do jfacet=1,nfacets ! OLD
              do jfacet=1,size(facets) ! NEW
                if (jfacet.ne.ifacet) then
!                  print*,"ifacet,jfacet=",ifacet,jfacet
!                  do i17=1,size(facets(jfacet)%facetpoints)
!                  print*,"jfacet,ifp,facetpoints=",                       &
!     &                jfacet,i17,facets(jfacet)%thisfacet(i17,:)
!                  end do
!                  if (is_visible(facets(jfacet)%thisfacet,                &
!     &                points(ifarthestpoint)%coord,midpoint)) then
!                  if (is_visible_3(facets(jfacet)%thisfacet,              &
!     &                points(ifarthestpoint)%coord,midpoint)) then
                  if (is_visible_4(facets(jfacet)%thisfacet,              &
     &               points(ifarthestpoint)%coord,midpoint,-tol))then
                    nvisible=nvisible+1
                    facets(jfacet)%visible=.true.
                    deallocate(V_temp)
                    allocate(V_temp(nvisible))
                    V_temp(1:nvisible-1)=V(:)
                    V_temp(nvisible)=facets(jfacet)
                    deallocate(V)
                    allocate(V(nvisible))
                    V=V_temp
                  end if ! is_outside
                  !print*, jfacet,facets(jfacet)%visible
                end if ! jfacet.ne.ifacet
              end do ! jfacet
              !print*,"V:"
              !do i=1,nvisible
              !  print*, V(i)%facetpoints 
              !end do
              !
              ! end find facets V visible by p
              !
              !
              !
              ! begin find horizon of V
              !
              !call get_ridges(V,ridges)
              call get_all_ridges(V,ridges)
              !call get_ridges(facets_new,ridges)
              ridges(:)%on_horizon=.false.
              nridges=size(ridges)
              !print*,"nridges=",nridges
!              do iridge=1,nridges
!              end do
              do iridge=1,nridges
!                if (facets(ridges(iridge)%belongs_to_facets(1))%visible   &
!     &             .neqv.facets(ridges(iridge)%belongs_to_facets(2))      &
!     &              %visible) then
                if (any(ridges(iridge)%belongs_to_facets(:)==0)) then 
                     ridges(iridge)%on_horizon=.true.
                     !print*,ridges(iridge)%ridgepoints
                end if
              end do ! iridge
              !
              ! end find horizon of V
              !
              !
              !
              ! begin add new facets from horizon to p
              !
              if (.not.allocated(newfacetpoints))                         &
     &            allocate(newfacetpoints(d))
              do iridge=1,nridges
                if (ridges(iridge)%on_horizon) then
                  newfacetpoints=(/ridges(iridge)%ridgepoints,            &
     &              ifarthestpoint/)
                  call add_to_facets(facets,newfacetpoints,points,        &
     &              facets_temp)
                  deallocate(facets)
                  allocate(facets(size(facets_temp)))
                  facets=facets_temp
                  !
                  ! BEGIN NEW
                  !
                  ! begin find points outside new facets
                  !
                   call get_outside_points(facets(size(facets):           &
     &               size(facets)),points,midpoint)
                  !    
                  ! end find points outside new facets 
                  !
                  ! END NEW
                  !
                  deallocate(facets_temp)
                end if ! on_horizon
              end do ! iridge
              points(ifarthestpoint)%dealtwith=.true.
              !
              ! end add new facets from horizon to p
              !
              !
              !
!              ! begin find points outside new facets
!              !
!              call get_outside_points(facets(nfacets+1:                   &
!     &               size(facets)),points,midpoint)
!              !    
!              ! end find points outside new facets 
              !
              !
              !
              ! begin remove inner facets
              !
              do i=1,nvisible
                call remove_from_facets(facets,V(i),facets_temp)
                deallocate(facets)
                allocate(facets(size(facets_temp)))
                facets=facets_temp
                deallocate(facets_temp)
              end do ! i
              !
              ! end remove inner facets
              !
              !
            end if ! while facets(ifacet)%noutside>0
            !
            ifacet=ifacet+1
            !
          end do ! ifacet
          !
          ntotoutside=sum(facets(:)%noutside)
          nfacets=size(facets)
          !
        end do ! while ntotoutside.gt.0
        !
        ! end rebuild facets and hull
        !  
        !  
        !  
 50     continue
        ! 
        ! begin write convex hull points
        !  
        if(allocated(hullpoints)) deallocate(hullpoints)
        call get_points(facets,hullpoints)
        nhullpoints=size(hullpoints)
        open(51,file='HULL_POINTS',status='replace')
        pointformat=' '
        !write(pointformat,'("(",I0,"(F10.6))")') d
        write(pointformat,'("(",I0,"(F20.12),2x,I0,1x,A)")') d
        !print*,pointformat
        do ipoint=1,nhullpoints
          !print*,hullpoints(ipoint)%pointnumber
          hullpoints(ipoint)%label=                                       &
     &       points(hullpoints(ipoint)%pointnumber)%label
          write(51,pointformat) hullpoints(ipoint)%coord(1:d),            &
     &           hullpoints(ipoint)%pointnumber,                          &
     &           adjustl(trim(hullpoints(ipoint)%label))
        end do ! ipoint
        close(51)
        !
        ! end write convex hull points
        !
        !
        !
        ! begin write convex hull facets
        !
        open(51,file='HULL_FACETS',status='replace')
        pointformat=' '
!        write(pointformat,'("(",I0,"(",I0,"(F20.12),3x),",I0,             &
!     &          "(I0,x)",")")') size(facets(1)%facetpoints),d,            &
!     &                          size(facets(1)%facetpoints) 
!        !print*,pointformat
!        do ifacet=1,nfacets
!          write(51,pointformat) (points(facets(ifacet)                    &
!     &    %facetpoints(i))%coord(1:d),i=1,size(facets(1)%facetpoints)),   &
!     &         facets(ifacet)%facetpoints(1:size(facets(1)%facetpoints))
!        end do ! ifacet
        write(pointformat,'("(",I0,"(",I0,"(F20.12),3x),",I0,             &
     &        "(I0,x),",I0,"(A,x)",")")') size(facets(1)%facetpoints),d,  &
     &          size(facets(1)%facetpoints),size(facets(1)%facetpoints) 
        !print*,pointformat
        do ifacet=1,nfacets
          write(51,pointformat) (points(facets(ifacet)                    &
     &    %facetpoints(i))%coord(1:d),i=1,size(facets(1)%facetpoints)),   &
     &       facets(ifacet)%facetpoints(1:size(facets(1)%facetpoints)),   &
     &    (adjustl(trim(points(facets(ifacet)%facetpoints(i))%label)),    &
     &     i=1,size(facets(1)%facetpoints))
        end do ! ifacet
        close(51)
        !
        ! end write convex hull facets
        !
        !
        !
        ! begin write convex hull ridges
        !
        if (allocated(ridges)) deallocate(ridges)
        !call get_ridges(facets,ridges)
        call get_all_ridges(facets,ridges)
        !print*, "ridges obtained for writing"
        nridges=size(ridges)
        open(51,file='HULL_RIDGES',status='replace')
        pointformat=' '
        !write(pointformat,'("(",I0,"(",I0,"(F10.6),3x))")') d-1,d
        nrp=size(ridges(1)%ridgepoints)
        write(pointformat,'("(",I0,"(",I0,"(F20.12),3x))")') nrp,d
        do iridge=1,nridges
          write(51,pointformat)                                           &
     &          (points(ridges(iridge)%ridgepoints(i))                    &
     &         %coord(1:d),i=1,nrp)
        end do ! iridge
        close(51)
        !
        ! end write convex hull ridges
        !
        !
        !
        if (talk) print fsubendext,"quickhull2"
        return
        ! 
        !
c---------------------------------------------------------------------

        contains

c---------------------------------------------------------------------

        subroutine read_points()
        implicit none
        type(point), allocatable :: points_temp(:)
        !!double precision, parameter :: tol=1.0E-12
        !double precision, parameter :: tol=1.0D-8
        !double precision, parameter :: tol=1.0D-12
        !double precision, parameter :: tol=5.0D-3 ! absolute tolerance to distinguish two points
        !double precision, parameter :: tol=4.0D-2 ! absolute tolerance to distinguish two points
        double precision, parameter :: tol=5.0D-2 ! absolute tolerance to distinguish two points
        !double precision, parameter :: tol=9.0D-2 ! absolute tolerance to distinguish two points
        double precision :: maxabspoint ! maximum distance between origin and point
        double precision :: tolrel ! relative tolerance to distinguish two points, tolrel=maxabspoint*tol
        integer ipoint,jpoint,kpoint
        logical newpoint
        double precision dist
        !
        !
        open(51,file=filename,status='old',err=101) 
        npoints=0
 10     read(51,'(A256)', err=102,end=12) line  
        line=trim(adjustl(line))
        if(len_trim(line).gt.0.and.line(1:1).ne.'#'                       &
     &     .and.line(1:1).ne.'!'.and.line(1:1).ne.'%') npoints=npoints+1
        goto 10
 12     continue
        !
!        if (npoints.lt.d+1) then
!          close(51)
!          call error_stop('too few points')
!        end if
        !
        allocate(points(npoints))
        do ipoint=1,npoints
          allocate(points(ipoint)%coord(1:d))
        end do
        ipoint=0
        maxabspoint=0.0d0
        rewind(51)
 14     read(51,'(A256)', err=102,end=16) line  
        line=trim(adjustl(line))
        if(len_trim(line).gt.0.and.line(1:1).ne.'#'                       &
     &     .and.line(1:1).ne.'!'.and.line(1:1).ne.'%') then
          ipoint=ipoint+1
          !read(line,*) points(ipoint)%coord(1:d)
          read(line,*) points(ipoint)%coord(1:d),points(ipoint)%label
          points(ipoint)%dealtwith=.false.
          points(ipoint)%onhull=1
          points(ipoint)%hasfacet=.false.
          if (norm2(points(ipoint)%coord).gt.maxabspoint)                 &
     &      maxabspoint=norm2(points(ipoint)%coord)
        end if
        goto 14
 16     continue
        close(51)
        !
        !
        tolrel=tol*maxabspoint 
        !print '(8x,"tol=",F20.16)',tol
        !print '(8x,"maxabspoint=",F20.16)',maxabspoint
        print '(8x,"tolrel=",F20.16)',tolrel
        !
        ! begin remove double points
        !
        !print*,"initial npoints:",npoints
        allocate(points_temp(npoints))
        open(51,file="POINTS.DISCARDED",status="replace")
        ipoint=1
        do kpoint=1,npoints
          newpoint=.true.
          do jpoint=1,ipoint-1
            dist=norm2(points(kpoint)%coord-points(jpoint)%coord)
            if (dist.lt.tolrel)  newpoint=.false.
          end do ! jpoint
          if (newpoint) then
            points_temp(ipoint)=points(kpoint)
            ipoint=ipoint+1
          else
            write(51,*) points(kpoint)%coord
          end if
        end do ! kpoint
        npoints=ipoint-1
        !print*,"final npoints:",npoints
        deallocate(points)
        allocate(points(npoints))
        points(1:npoints)=points_temp(1:npoints)
        deallocate(points_temp)
        close(51)
        !
        ! end remove double points
        !
        !
        !
        return

 101    close(51)
        call error_stop('file not found')
 102    close(51)
        call error_stop('file empty?')

      end subroutine read_points

c---------------------------------------------------------------------

        end subroutine quickhull2

c---------------------------------------------------------------------

        subroutine get_simplex_2(points,simplex)        
        !
        ! sets up a convex polyhedron (simplex) of d+1 points from of a set of points (not
        ! identical with the convex hull of the points, unless there are
        ! only <= d+1 points, d=dimension)
        !
        use defs, only: fsubstart,fsubendext,ncomm,fcomm
        implicit none
        !
        type(point) :: points(:)
        type(point), allocatable :: simplex(:)
        ! local:
        type(point), allocatable :: simplex_temp(:)
        type(facet), allocatable :: facets(:)
        integer nsp,np,nsp_temp,d
        integer nfacets,ifacet,jfacet,nfp
        integer ipoint,jpoint
        double precision distance,distance2
        double precision, allocatable :: midpoint(:),crosspoint(:)
        double precision, allocatable :: vec(:)
        double precision :: norm2
        logical isoutside
        !double precision :: tol=1D-8
        double precision :: tol=1D-12
        !
        ! begin initialize
        !
        print fsubstart,'get_simplex_2'
        !
        np=size(points,1) ! number of points
        d=size(points(1)%coord)
        print '(8x,I0," points")',np
        nsp=min(d+1,np) ! number of simplex points
        print '(8x,I0," simplex points")',nsp
        if (nsp==0) then
          ncomm=ncomm+1
          print fcomm,'You supplied zero points'
          return
        end if
        if (nsp==np) then
          allocate(simplex(nsp))
          simplex=points
          return
        end if
        !print '(8x,"going to get simplex points")'
        !
        ! end initialize
        !
        allocate(simplex(2))
        allocate(simplex(1)%coord(d))
        allocate(simplex(2)%coord(d))
        nsp_temp=0
        distance=0.0d0
        allocate(vec(d))
        !
        ! begin find first 2 points (the ones farthest apart)
        !print '(8x,"going to get first 2 simplex points")'
        !
        do ipoint=1,np
          do jpoint=ipoint+1,np
            vec=points(ipoint)%coord(:)-points(jpoint)%coord(:)
            if (norm2(vec).gt.distance) then
              distance=norm2(vec)
              simplex(1)%coord(1:d)=points(ipoint)%coord(1:d)
              simplex(1)%label=points(ipoint)%label
              simplex(1)%pointnumber=ipoint
              simplex(2)%coord(1:d)=points(jpoint)%coord(1:d)
              simplex(2)%label=points(jpoint)%label
              simplex(2)%pointnumber=jpoint
              !simplex(1)=points(ipoint)
              !simplex(1)%pointnumber=ipoint
              !simplex(2)=points(jpoint)
              !simplex(2)%pointnumber=jpoint
            end if
          end do ! jpoint
        end do ! ipoint
        !print '(8x,"simplex initialized with 2 points")'
        !
        ! end find first 2 points (the ones farthest apart)
        !
        nsp_temp=size(simplex)
        allocate(midpoint(d))
        !
        ! begin successively add most distance point until simplex has
        ! enough points
        !
        do while (nsp_temp<nsp)
          !print*," "
          !print '(8x,"simplex points:")'
          !do jpoint=1,nsp_temp
            !print '(8x,10(F10.5))',simplex(jpoint)%coord
          !end do
          !print '(8x,"nsp_temp=",I0)',nsp_temp
          !
          allocate(simplex_temp(nsp_temp+1))
          simplex_temp(1:nsp_temp)=simplex
          if(allocated(facets)) deallocate(facets)
          call get_facets(simplex,facets)
          print '(8x,"simplex facets obtained")'
          nfacets=size(facets)
          print '(8x,"nfacets=",I0)',nfacets
          do jfacet=1,nfacets
            !print '(8x,I10)',facets(jfacet)%facetpoints
            nfp=size(facets(jfacet)%facetpoints)
            do jpoint=1,nfp
              !print '(8x,10(F10.5))',facets(jfacet)%thisfacet(jpoint,:)
            end do
          end do
          !
          ! begin get COM of simplex
          !
          midpoint=0.0d0
          do jpoint=1,nsp_temp
            midpoint=midpoint+simplex(jpoint)%coord
          end do
          midpoint=midpoint/dble(nsp_temp)  
          !print '(8x,"midpoint obtained")'
          !print '(8x,"midpoint=",10(F10.5))',midpoint
          !
          ! end get COM of simplex
          !
          !
          distance=0.0d0
          do ipoint=1,np
            !print '(8x,"ipoint=",I0)',ipoint
            !print '(8x,10(F10.5))',points(ipoint)%coord
            !
            isoutside=.true.
            !
            if (.not.any(simplex%pointnumber==ipoint)) then
              !print '(8x,"point is new")'
              do ifacet=1,nfacets
                !print '(8x,"ifacets=",I0)',ifacet
                nfp=size(facets(ifacet)%facetpoints)
                !print '(8x,"nfp=",I0)',nfp
                !
                ! begin get distance of point from facet
                !
                call distance_from_surface(points(ipoint)%coord,          &
     &           facets(ifacet)%thisfacet,crosspoint,distance2)
                !print '(8x,"distance=",F10.5)',distance2
                !
                ! end get distance of point from facet
                !
                ! begin check if point sees facet
                !
                if (.not.is_visible_4(facets(ifacet)%thisfacet,           &
     &            points(ipoint)%coord,midpoint,tol)) isoutside=.false.
                !print '(8x,"point sees facet:",L1)',isoutside
                !
                ! end check if point sees facet
                !
              end do ! ifacet  
              ! 
              ! begin if point is outside simplex and in larger
              ! distance to it, take it
              !
              if (isoutside.and.distance2>=distance) then
                !print '(8x,"taking point ")'
                !
                !simplex_temp(nsp_temp+1)=points(ipoint)
                if(allocated(simplex_temp(nsp_temp+1)%coord))             &
     &             deallocate(simplex_temp(nsp_temp+1)%coord)           
                allocate(simplex_temp(nsp_temp+1)%coord(d))           
                simplex_temp(nsp_temp+1)%coord=points(ipoint)%coord
                simplex_temp(nsp_temp+1)%label=points(ipoint)%label
                simplex_temp(nsp_temp+1)%pointnumber=ipoint
                distance=distance2
                !
              end if ! isoutside and in larger distance
              !
            end if ! point not yet part of simplex
            !
          end do ! ipoint
          !
          ! begin update simplex
          !
          deallocate(simplex)
          allocate(simplex(size(simplex_temp)))
          simplex=simplex_temp
          nsp_temp=size(simplex)
          deallocate(simplex_temp)
          !
          ! end update simplex
          !
        end do ! nsp_temp<nsp
        !
        ! begin successively add most distance point until simplex has
        ! enough points
        !
        ! begin set simplex points to "dealtwith"
        !
        do ipoint=1,nsp
          !
          points(simplex(ipoint)%pointnumber)%dealtwith=.true.
          points(simplex(ipoint)%pointnumber)%onhull=0
          ! 
        end do
        !
        ! end set simplex points to "dealtwith"

        print*," "
        print '(8x,"simplex points:")'
        do jpoint=1,nsp_temp
          print '(8x,10(F10.5))',simplex(jpoint)%coord
        end do
        !
        print fsubendext,'get_simplex_2'
        !
        end subroutine get_simplex_2
        
c---------------------------------------------------------------------

        subroutine get_facets(points,facets)
        use combi, only : n_over_k
        implicit none
        type(point), intent(in) :: points(:)
        type(facet), allocatable :: facets(:)
        !
        ! local :
        integer nhp,nfp,nfacets
        integer, allocatable :: indices(:)
        integer i,i_ind,j_ind,k_ind,ifacet,d,ipoint
        logical newfacet

        d=size(points(1)%coord)
        nhp=size(points)
        nfp=min(d,nhp)
        nfacets=n_over_k(nhp,nfp)
        !
        !
        allocate(indices(nfp))
        do i_ind=1,size(indices)
            indices(i_ind)=1
        end do ! i_ind
        i_ind=size(indices)

        allocate(facets(nfacets))
!        do ifacet=1,nfacets
!          allocate(facets(ifacet)%facetpoints(1:nfp))
!          facets(ifacet)%visible=.false.
!          allocate(facets(ifacet)%thisfacet(nfp,d))
!        end do
        !
        ifacet=1
        do while (ifacet<=nfacets)
            !
            ! begin consider only tuples of increasing numbers, e.g. (1 2 8 in 3D), (2 3 9 26 in 4D), ...
            !
            newfacet=.true.
            do j_ind=1,size(indices)
              do k_ind=j_ind+1,size(indices)
                if (indices(j_ind)>=indices(k_ind)) newfacet=.false.
              end do ! k_ind
            end do ! j_ind
            !
            ! begin consider only tuples of increasing numbers, e.g. (1 2 8 in 3D), (2 3 9 26 in 4D), ...
            !
            if (newfacet) then
              !
              ! begin create facet
              !
              allocate(facets(ifacet)%facetpoints(1:nfp))
              allocate(facets(ifacet)%thisfacet(nfp,d))
              facets(ifacet)%facetpoints(1:nfp)=                          &
     &          points(indices(1:nfp))%pointnumber
              facets(ifacet)%facetnumber=ifacet
              do j_ind=1,nfp
                facets(ifacet)%thisfacet(j_ind,1:d)=                      &
     &            points(indices(j_ind))%coord(1:d)
              end do ! ipoint
              facets(ifacet)%visible=.false.
              facets(ifacet)%noutside=0
              !
              ! end create facet
              !
              ifacet=ifacet+1
              !
            end if ! newfacet
            !
            do while (all(indices(i_ind:nfp)==nhp.and.i_ind>1))
                i_ind=i_ind-1
            end do
            !
            do j_ind=i_ind+1,size(indices)
                indices(j_ind)=1
            end do
            !
            indices(i_ind)=indices(i_ind)+1
            i_ind=size(indices)
            !
        end do ! ifacet <= nfacets
        !
        !
!        do ifacet=1,nfacets
!          facets(ifacet)%facetnumber=ifacet
!          do ipoint=1,nfp
!            facets(ifacet)%thisfacet(ipoint,1:d)=                         &
!     &           points(facets(ifacet)%facetpoints(ipoint))%coord(1:d)
!          end do ! ipoint
!          facets(ifacet)%noutside=0
!        end do ! ifacet
        !
        end subroutine get_facets

c---------------------------------------------------------------------

        subroutine distance_from_surface(point,surfpoints,                &
     &                       crosspoint,dist)
        ! calculates the distance of a point from a surface. Assumes
        ! that the surface is plane, so that there exists a direction 
        ! perpendicular to the surface. The surface can have lower
        ! dimension (for example, a straight line in 3D is possible) 
        use defs, only : fsubstart, fsubendext, error_stop
        implicit none
        double precision, intent(in) :: point(:)
        double precision, intent(in) :: surfpoints(:,:)
        double precision, allocatable :: crosspoint(:)
        double precision :: dist
        !
        ! local variables :
        integer :: d,nsp
        integer i,j,k,ds,nvpar
        double precision, allocatable :: vperp(:),vunitpar(:,:)
        !double precision, parameter :: tol=1.0D-8
        double precision, parameter :: tol=1.0D-12
        !
        !print fsubstart,'distance_from_surface'
        !
        ! begin get dimension and surface dimension
        ds=size(surfpoints,1)
        d=size(point)
        if (d.ne.size(surfpoints,2))   then
            call error_stop('distance_from_surface: inconsistent dimensi  &
     &ons')
        end if
        allocate(vperp(d))
        if (allocated(crosspoint)) deallocate(crosspoint)
        allocate(crosspoint(d))
        !nvpar=ds*(ds-1)/2
        !nvpar=min(d-1,ds*(ds-1)/2)
        nvpar=min(d-1,ds-1)
        !print '(8x,I0," dimensions")',d
        !print '(8x,I0," surface points")',ds
        !print '(8x,I0," vectors")',nvpar
        allocate(vunitpar(nvpar,d))
        ! end get dimension and surface dimension
        !
        ! begin set up unit vectors that span surface
        !
        k=1
        !
        do i=1, ds
            do j=i+1,ds
                if (k<=nvpar) then
                  vunitpar(k,:)=surfpoints(i,:)-surfpoints(j,:)
                  if (norm2(vunitpar(k,:)).gt.tol) then
                    vunitpar(k,:)=vunitpar(k,:)/norm2(vunitpar(k,:))
                  end if
                end if ! k<=nvpar
                k=k+1
            end do ! j
        end do ! i
        !
        !print '(8x,"vectors set up")'
        ! begin orthogonalize vectors
        !
        do i=2,nvpar
            do j=1,i-1
                vunitpar(i,:)=vunitpar(i,:)                               &
     &           -dot_product(vunitpar(i,:),vunitpar(j,:))*vunitpar(j,:)
                vunitpar(i,:)=vunitpar(i,:)/norm2(vunitpar(i,:))
            end do ! j
        end do ! i
        !print '(8x,"vectors orthogonalized")'
        ! end orthogonalize vectors
        !
        ! end set up unit vectors that span surface
        !
        ! begin subtract components in plane
        vperp=surfpoints(1,:)-point(:)
        !
        do k=1,nvpar
          vperp=vperp-dot_product(vperp,vunitpar(k,:))*vunitpar(k,:)
        end do ! k
        !
        dist=norm2(vperp)
        crosspoint=point+vperp
        !
        !print fsubendext,'distance_from_surface'
        !
        end subroutine distance_from_surface

c---------------------------------------------------------------------

        subroutine get_outside_points(facets,points,midpoint)
        implicit none
        type(facet), intent(inout) :: facets(:)
        type(point), intent(inout) :: points(:)
        double precision, intent(in) :: midpoint(:)
        ! internal
        integer ipoint,ifacet,nfacets,npoints,d,nfp
        integer, allocatable :: outsidepoints_temp(:)

        npoints=size(points)
        nfacets=size(facets)
        d=size(midpoint)

        points(:)%hasfacet=.false.
        do ifacet=1,nfacets
          facets(ifacet)%noutside=0
          nfp=size(facets(ifacet)%facetpoints)
          do ipoint=1,npoints
            if (.not.points(ipoint)%dealtwith) then
              if (.not.points(ipoint)%hasfacet) then
!                if (is_visible(facets(ifacet)%thisfacet,                  &
!     &              points(ipoint)%coord,midpoint)) then
!                if (is_visible_3(facets(ifacet)%thisfacet,                &
!     &              points(ipoint)%coord,midpoint)) then
                if (is_visible_4(facets(ifacet)%thisfacet,                &
     &              points(ipoint)%coord,midpoint,-1.0D-12)) then
                  points(ipoint)%hasfacet=.true.
                  ! add point to outsidepoints(facet):
                  if (facets(ifacet)%noutside.gt.0) then
                    if(allocated(outsidepoints_temp))                     &
     &                 deallocate(outsidepoints_temp)
                    allocate(outsidepoints_temp(facets(ifacet)%           &
     &                noutside))
                    outsidepoints_temp=facets(ifacet)%outsidepoints
                    deallocate(facets(ifacet)%outsidepoints)
                  end if ! noutside > 0
                  facets(ifacet)%noutside=facets(ifacet)%noutside+1
                  if (allocated(facets(ifacet)%outsidepoints))            &
     &               deallocate(facets(ifacet)%outsidepoints)
                  allocate(facets(ifacet)%outsidepoints                   &
     &                  (facets(ifacet)%noutside))
                  if (facets(ifacet)%noutside.gt.1) then
                    facets(ifacet)%outsidepoints=(/outsidepoints_temp,    &
     &                        ipoint /)              
                    deallocate(outsidepoints_temp)
                  else
                    facets(ifacet)%outsidepoints(1)=ipoint
                  end if ! noutside > 1
                end if ! is visible
                !print*,ifacet,ipoint,points(ipoint)%hasfacet
              end if ! not has_facet
            end if ! point not dealt with
          end do ! ipoint
        end do ! ifacet

        end subroutine get_outside_points

c---------------------------------------------------------------------

c---------------------------------------------------------------------
!
!        logical function is_visible(facet,point,refpoint) 
!        ! point sees facet?  
!        use defs, only : error_stop
!        implicit none
!        double precision :: norm2
!        double precision facet(:,:),point(:),refpoint(:)
!        double precision, allocatable :: vecpf(:),vecrf(:)
!        double precision, allocatable :: vec12(:),vec13(:)
!        double precision, allocatable :: vecpfpar(:),vecpfperp(:)
!        double precision, allocatable :: vecrfpar(:),vecrfperp(:)
!        double precision, allocatable :: vecperp(:)
!        integer d
!        is_visible=.false.
!        d=size(point)
!        select case(d)
!        case(2)
!          allocate(vec12(1:d))
!          allocate(vecpf(1:d),vecpfpar(1:d),vecpfperp(1:d))
!          allocate(vecrf(1:d),vecrfpar(1:d),vecrfperp(1:d))
!          vecpf(:)=-facet(1,:)+point(:)
!          vecrf(:)=-facet(1,:)+refpoint(:)
!          vec12=facet(2,:)-facet(1,:)
!          vecpfpar=vec12*dot_product(vecpf,vec12)/norm2(vec12)**2      
!          vecrfpar=vec12*dot_product(vecrf,vec12)/norm2(vec12)**2      
!          vecpfperp=vecpf-vecpfpar
!          vecrfperp=vecrf-vecrfpar
!        case(3)
!          allocate(vec12(1:d))
!          allocate(vecpf(1:d),vecpfpar(1:d),vecpfperp(1:d))
!          allocate(vecrf(1:d),vecrfpar(1:d),vecrfperp(1:d))
!          allocate(vec13(1:d),vecperp(1:d))
!          vecpf(:)=-facet(1,:)+point(:)
!          vecrf(:)=-facet(1,:)+refpoint(:)
!          vec12=facet(2,:)-facet(1,:)
!          vec13=facet(3,:)-facet(1,:)
!          vecpfpar=vec12*dot_product(vecpf,vec12)/norm2(vec12)**2      
!          vecrfpar=vec12*dot_product(vecrf,vec12)/norm2(vec12)**2      
!          vecperp=vec13-dot_product(vec13,vec12)*vec12                    &
!     &            /norm2(vec12)**2
!          vecpfperp=vecpf-vecpfpar-vecperp*dot_product(vecpf,vecperp)     &
!     &             /norm2(vecperp)**2
!          vecrfperp=vecrf-vecrfpar-vecperp*dot_product(vecrf,vecperp)     &
!     &             /norm2(vecperp)**2
!        case default
!          call error_stop('dimension not implemented')
!        end select
!        if (dot_product(vecpfperp,vecrfperp).lt.0.0d0) is_visible=.true.
!        end function is_visible
!
c---------------------------------------------------------------------
!
!        logical function is_visible_2(ridge,point,refpoint) 
!        ! point sees ridge?  
!        use defs, only : error_stop
!        implicit none
!        double precision :: norm2
!        double precision ridge(:,:),point(:),refpoint(:)
!        double precision, allocatable :: vecp1(:),vecr1(:)
!        double precision, allocatable :: vec12(:),vecperp1(:)
!        double precision, allocatable :: vecperp2(:)
!        double precision, allocatable :: vecpar1(:),vecpar2(:)
!        integer d
!        double precision, parameter :: tol=1.0E-10
!        is_visible_2=.false.
!        d=size(point)
!        select case(d)
!        case(2)
!          allocate(vecr1(1:d))
!          allocate(vecp1(1:d))
!          vecp1(:)=-point(:)+ridge(1,:)
!          vecr1(:)=-refpoint(:)+ridge(1,:)
!          !print*,"midpoint:",refpoint
!          !print*,"ridge:",ridge(1,:)
!          !print*,"dot_product:",dot_product(vecp1,vecr1)
!          if (dot_product(vecp1,vecr1).lt.-tol) is_visible_2=.true.
!        case(3)
!          allocate(vecr1(1:d))
!          allocate(vecp1(1:d))
!          allocate(vec12(1:d))
!          allocate(vecperp1(1:d))     
!          allocate(vecperp2(1:d))     
!          allocate(vecpar1(1:d))     
!          allocate(vecpar2(1:d))     
!          vec12=-ridge(1,:)+ridge(2,:)
!          vecr1=-refpoint+ridge(1,:)
!          vecp1=-point+ridge(1,:)
!          vecpar1=vec12*dot_product(vecr1,vec12)/norm2(vec12)**2
!          vecperp1=vecr1-vecpar1
!          vecpar2=vec12*dot_product(vecp1,vec12)/norm2(vec12)**2
!          vecperp2=vecp1-vecpar2
!          if (dot_product(vecperp1,vecperp2).lt.-tol) then
!            is_visible_2=.true.
!          end if
!        case default
!          call error_stop('dimension not implemented')
!        end select
!        !print*,"is_visible:",is_visible_2
!        end function is_visible_2
!
c---------------------------------------------------------------------
!
!        logical function is_visible_3(surfpoints,point,refpoint) 
!        ! point sees ridge?  
!        use defs, only : error_stop,fsubstart,fsubendext
!        implicit none
!        !integer, intent(in) :: nsp,d
!        double precision :: norm2
!        double precision :: surfpoints(:,:),point(:),                     &
!     &         refpoint(:)
!!        double precision, intent(in) :: surfpoints(nsp,d),point(d),       &
!!     &         refpoint(d)
!        double precision, allocatable :: crosspoint_P(:),crosspoint_M(:)
!        double precision, allocatable :: nvec_P(:), nvec_M(:)
!        double precision dist
!        integer dsp,ds,d!,nsp
!        double precision, parameter :: tol=1.0E-8
!        !
!        ! begin initialize
!        !
!        !print fsubstart,'is_visible_3'
!        !
!        is_visible_3=.false.
!        d=size(point)
!        ds=size(surfpoints,1)
!        dsp=size(surfpoints,2)
!        !print '(8x,I0, "dimensions")',d
!        !print '(8x,I0, " surface points")', size(surfpoints,1)
!        !print '(8x,I0, " surface dimension")', size(surfpoints,2)
!        !print*,"surface points=",surfpoints
!        if (.not.d==dsp) then
!        !if (.not.d==size(surfpoints,2)) then
!          !print '(8x,I0, " dimensions")',d
!          !print '(8x,I0, " surface points")', size(surfpoints,1)
!          !print '(8x,I0, " surface dimension")', size(surfpoints,2)
!          !print*,"surface points=",surfpoints
!          call error_stop('inconsistent dimensions')
!        end if
!        allocate(nvec_P(d),nvec_M(d))
!        !allocate(crosspoint_P(d),crosspoint_M(d))
!        !
!        ! end initialize
!        !
!        call distance_from_surface(point,surfpoints,                      &
!     &        crosspoint_P,dist)
!        nvec_P=crosspoint_P-point
!        call distance_from_surface(refpoint,surfpoints,                   &
!     &        crosspoint_M,dist)
!        nvec_M=crosspoint_M-refpoint
!        !print '(8x,"nvec_P=",10(F10.5))',nvec_P
!        !print '(8x,"crosspoint_P=",10(F10.5))',crosspoint_P
!        !print '(8x,"nvec_M=",10(F10.5))',nvec_M
!        !print '(8x,"crosspoint_M=",10(F10.5))',crosspoint_M
!        !print '(8x,"nvec_P*nvec_M=",F10.5)',dot_product(nvec_M,nvec_P)
!        !
!        if (dot_product(nvec_P,nvec_M).lt.-tol) is_visible_3=.true.
!        !
!        !print fsubendext,'is_visible_3'
!        !
!        end function is_visible_3
!
c---------------------------------------------------------------------

        logical function is_visible_4(surfpoints,point,refpoint,tol) 
        ! point sees ridge?  
        use defs, only : error_stop,fsubstart,fsubendext
        implicit none
        !integer, intent(in) :: nsp,d
        double precision :: norm2
        double precision :: surfpoints(:,:),point(:),                     &
     &         refpoint(:)
!        double precision, intent(in) :: surfpoints(nsp,d),point(d),       &
!     &         refpoint(d)
        double precision, allocatable :: crosspoint_P(:),crosspoint_M(:)
        double precision, allocatable :: nvec_P(:), nvec_M(:)
        double precision dist
        integer dsp,ds,d!,nsp
        double precision :: tol
        !
        ! begin initialize
        !
        !print fsubstart,'is_visible_4'
        !
        is_visible_4=.false.
        d=size(point)
        ds=size(surfpoints,1)
        dsp=size(surfpoints,2)
        !print '(8x,I0, "dimensions")',d
        !print '(8x,I0, " surface points")', size(surfpoints,1)
        !print '(8x,I0, " surface dimension")', size(surfpoints,2)
        !print*,"surface points=",surfpoints
        if (.not.d==dsp) then
        !if (.not.d==size(surfpoints,2)) then
          !print '(8x,I0, " dimensions")',d
          !print '(8x,I0, " surface points")', size(surfpoints,1)
          !print '(8x,I0, " surface dimension")', size(surfpoints,2)
          !print*,"surface points=",surfpoints
          print '(8x,I0,", ",I0)', d,dsp
          call error_stop('is_visible_4: inconsistent dimensions')
        end if
        allocate(nvec_P(d),nvec_M(d))
        !allocate(crosspoint_P(d),crosspoint_M(d))
        !
        ! end initialize
        !
        call distance_from_surface(point,surfpoints,                      &
     &        crosspoint_P,dist)
        nvec_P=crosspoint_P-point
        call distance_from_surface(refpoint,surfpoints,                   &
     &        crosspoint_M,dist)
        nvec_M=crosspoint_M-refpoint
        !print '(8x,"nvec_P=",10(F10.5))',nvec_P
        !print '(8x,"crosspoint_P=",10(F10.5))',crosspoint_P
        !print '(8x,"nvec_M=",10(F10.5))',nvec_M
        !print '(8x,"crosspoint_M=",10(F10.5))',crosspoint_M
        !print '(8x,"nvec_P*nvec_M=",F10.5)',dot_product(nvec_M,nvec_P)
        !
        if (dot_product(nvec_P,nvec_M).lt.tol) is_visible_4=.true.
        !
        !print fsubendext,'is_visible_4'
        !
        end function is_visible_4

c---------------------------------------------------------------------
!
!        double precision function distance_fp(thefacet,thepoint)
!          use defs, only : error_stop
!        implicit none
!        double precision  :: norm2
!        type(point), intent(in) :: thepoint
!        type(facet), intent(in) :: thefacet
!        integer d
!        double precision, allocatable :: vec12(:),vecpar(:),vecperp(:)
!        double precision, allocatable :: vec(:),vec13(:),vecperp2(:)
!        d=size(thepoint%coord)
!        select case(d)
!        case(2)
!          distance_fp=0.0d0
!          allocate(vec12(1:d),vecpar(1:d),vecperp(1:d))
!          vec12=thefacet%thisfacet(2,1:d)-thefacet%thisfacet(1,1:d)
!          vec=thepoint%coord(:)-thefacet%thisfacet(1,:)
!          vecpar=vec12*dot_product(vec,vec12)/norm2(vec12)**2        
!          vecperp=vec-vecpar
!          distance_fp=norm2(vecperp)
!        case(3)
!          allocate(vec13(1:d),vecperp2(1:d))
!          vec12=thefacet%thisfacet(2,:)-thefacet%thisfacet(1,:)
!          vec13=thefacet%thisfacet(3,:)-thefacet%thisfacet(1,:)
!          vec=thepoint%coord(:)-thefacet%thisfacet(1,:)
!          vecpar=vec12*dot_product(vec,vec12)/norm2(vec12)**2   
!          vecperp=vec13-dot_product(vec13,vec12)*vec12                    &
!     &            /norm2(vec12)**2
!          vecperp2=vec-vecpar-vecperp*dot_product(vec,vecperp)            &
!     &             /norm2(vecperp)**2
!          distance_fp=norm2(vecperp2)
!        case default
!          call error_stop('dimension not implemented')
!        end select
!        end function distance_fp
!
c---------------------------------------------------------------------


        subroutine add_to_facets(facets,facetpoints,points,newfacets)
        use defs, only : error_stop, fsubendext,fsubstart
        implicit none
        type(point), intent(in) :: points(:)
        type(facet), intent(in) :: facets(:)
        integer, intent(in) ::  facetpoints(:)
        ! local variables:
        type(facet), allocatable :: newfacets(:)
        integer nfacets,d,ipoint,nfacetpoints,ifacet
        logical newfacet
        !
        !print fsubstart,'add_to_facets'
        d=size(points(1)%coord)
        !print*,"d=",d
        !if (d>3) call error_stop('dim not implemented')
        nfacets=size(facets)
        nfacetpoints=size(facetpoints)
        !print*,"nfacets=",nfacets
        !print*,"newfacetpoints=",nfacetpoints
        !print*,"new facet points=",facetpoints
        newfacet=.true.
        do ifacet=1,nfacets
          if (all(facets(ifacet)%facetpoints==facetpoints)) then
            newfacet=.false.
          end if
        end do ! ifacet

        if (newfacet) then
          allocate(newfacets(nfacets+1))
          newfacets(1:nfacets)=facets(1:nfacets)
          allocate(newfacets(nfacets+1)%facetpoints(nfacetpoints))
          newfacets(nfacets+1)%facetpoints=facetpoints
          allocate(newfacets(nfacets+1)%thisfacet(nfacetpoints,1:d))
          !print*,"facets copied"
          do ipoint=1,nfacetpoints
            newfacets(nfacets+1)%thisfacet(ipoint,:)                      &
     &        =points(facetpoints(ipoint))%coord(:)
          end do ! ipoint
          newfacets(nfacets+1)%visible=.false.
          newfacets(nfacets+1)%facetnumber=nfacets+1
          newfacets(nfacets+1)%noutside=0
        else
          allocate(newfacets(nfacets))
          newfacets(1:nfacets)=facets(1:nfacets)
        end if
        !print fsubendext,'add_to_facets'
        end subroutine add_to_facets

c---------------------------------------------------------------------        

        subroutine remove_from_facets(facets,afacet,newfacets)
        use defs, only : error_stop, fsubendext,fsubstart
        implicit none
        type(facet), intent(in) :: facets(:),afacet
        type(facet), allocatable :: newfacets(:)
        ! local variables:
        integer nfacets,d,ipoint,nfacetpoints
        integer jfacet,ifacet
        !
        !print fsubstart,'remove_from_facets'
        d=size(afacet%thisfacet,2)
        !print*,"d=",d
        !if (d>3) call error_stop('dim not implemented')
        nfacets=size(facets)
        nfacetpoints=size(afacet%thisfacet,1)
        !print*,"nfacets=",nfacets
        !print*,"n_facetpoints_to_remove=",nfacetpoints
        !print*,"facetpoints",afacet%facetpoints
        allocate(newfacets(nfacets-1))
        ifacet=1
        do jfacet=1,nfacets
          !print*,"facet",facets(jfacet)%facetpoints
          if (.not.all(facets(jfacet)%facetpoints==                       &
     &        afacet%facetpoints)) then
            newfacets(ifacet)=facets(jfacet)
            ifacet=ifacet+1
            !print*,"keeping facet."
          else
            !print*,"removing facet."
          end if
        end do ! jfacet
        !print fsubendext,'remove_from_facets'
        end subroutine remove_from_facets

c---------------------------------------------------------------------        

        subroutine get_points(facets,points)
        use defs, only : error_stop, fsubstart
        implicit none
        type(facet), intent(in) :: facets(:)
        type(point), allocatable, intent(out) :: points(:)
        ! local variables:
        type(point), allocatable :: points_temp(:)
        integer npoints,nfacets,d,ipoint,ifacet,jpoint,kpoint
        logical newpoint

        !print fsubstart,'get_points'
        nfacets=size(facets)
        d=size(facets(1)%thisfacet(1,:))
        npoints=nfacets*d
        if (allocated(points)) deallocate(points)
        allocate(points(1:npoints))
        do ipoint=1,npoints
          allocate(points(ipoint)%coord(d))
        end do
        !print '("initialized...")'
        ipoint=1
        !select case(d)
        !case(2,3)
          do ifacet=1,nfacets
            do jpoint=1,size(facets(ifacet)%facetpoints)
              newpoint=.true.
              do kpoint=1,ipoint-1
                if (all(points(kpoint)%coord(:)==                         &
     &                facets(ifacet)%thisfacet(jpoint,:)) ) then 
                  newpoint=.false.
                end if  
              end do ! kpoint
              if (newpoint) then
                points(ipoint)%coord(:)                                   &
     &              =facets(ifacet)%thisfacet(jpoint,:)
                points(ipoint)%pointnumber                                &
     &              =facets(ifacet)%facetpoints(jpoint)
                ipoint=ipoint+1
              end if ! new point
            end do ! jpoint
          end do ! ifacet
        !case default
        !  call error_stop('dim not implemented')
        !end select
        npoints=ipoint-1
        allocate(points_temp(npoints))
        points_temp(1:npoints)=points(1:npoints)
        deallocate(points)
        allocate(points(npoints))
        points=points_temp
        deallocate(points_temp)
        end subroutine get_points

c---------------------------------------------------------------------        
!
!        subroutine get_ridges(facets,ridges)
!        use defs, only : error_stop, fsubstart
!        implicit none
!        type(facet), intent(in) :: facets(:)
!        type(ridge), allocatable, intent(out) :: ridges(:)
!        ! local variables:
!        type(ridge), allocatable :: ridges_temp(:)
!        integer nridges,nfacets,d,iridge,ifacet,ipoint,jpoint
!        integer jridge
!        logical newridge
!
!        !print fsubstart,'get_ridges'
!        nfacets=size(facets)
!        d=size(facets(1)%thisfacet(1,:))
!        nridges=nfacets*d
!        if (allocated(ridges)) deallocate(ridges)
!        allocate(ridges(1:nridges))
!        do iridge=1,nridges
!          allocate(ridges(iridge)%ridgepoints(d-1))
!        end do
!        !print '("initialized...")'
!        iridge=1
!        select case(d)
!        case(2)
!          do ifacet=1,nfacets
!            do ipoint=1,size(facets(ifacet)%facetpoints)
!              newridge=.true.
!              do jridge=1,iridge-1
!!                print*,"iridge,ifacet,ipoint,jridge:",iridge,ifacet,      &
!!     &            ipoint,jridge
!                if (ridges(jridge)%ridgepoints(1)==                       &
!     &                facets(ifacet)%facetpoints(ipoint) ) then 
!                  newridge=.false.
!                  ridges(jridge)%belongs_to_facets(2)=ifacet
!                end if  
!              end do ! jridge
!              !print*,"newridge:",newridge
!              if (newridge) then
!                ridges(iridge)%ridgepoints(1:d-1)                         &
!     &              =facets(ifacet)%facetpoints(ipoint)
!                !print*,"new ridge added."
!                ridges(iridge)%belongs_to_facets(1)=ifacet
!                !print*,"facet assigned."
!                ! begin get ridge point coordinates
!!                allocate(ridges(iridge)%thisridge(d-1,d))
!!                ridges(iridge)%thisridge(1,:)=                          &
!!     &          facets(ifacet)%thisfacet(ipoint,:)                
!                ! end get ridge point coordinates
!                iridge=iridge+1
!                !print*,"counter raised."
!              end if ! new ridge
!            end do ! ipoint
!          end do ! ifacet
!        case(3)
!          do ifacet=1,nfacets
!            do ipoint=1,size(facets(ifacet)%facetpoints)
!              do jpoint=ipoint+1,size(facets(ifacet)%facetpoints)
!                newridge=.true.
!                do jridge=1,iridge-1
!                  if (all(ridges(jridge)%ridgepoints(1:d-1)==           &
!     &                (/facets(ifacet)%facetpoints(ipoint),             &
!     &                  facets(ifacet)%facetpoints(jpoint) /))) then 
!                    newridge=.false.
!                    ridges(jridge)%belongs_to_facets(2)=ifacet
!                  end if  
!                end do ! jridge
!                if (newridge) then
!                  ridges(iridge)%ridgepoints(1:d-1)                     &
!     &              =(/facets(ifacet)%facetpoints(ipoint),              &
!     &                 facets(ifacet)%facetpoints(jpoint)/)
!                  ridges(iridge)%belongs_to_facets(1)=ifacet
!                  ! begin get ridge point coordinates
!!                  allocate(ridges(iridge)%thisridge(d-1,d))
!!                  ridges(iridge)%thisridge(1,:)=                          &
!!     &            facets(ifacet)%thisfacet(ipoint,:)                
!!                  ridges(iridge)%thisridge(2,:)=                          &
!!     &            facets(ifacet)%thisfacet(jpoint,:)
!                  ! end get ridge point coordinates
!                  iridge=iridge+1
!                end if ! new ridge
!              end do ! jpoint
!            end do ! ipoint
!          end do ! ifacet
!        case default
!          call error_stop('dim not implemented')
!        end select
!        nridges=iridge-1
!        allocate(ridges_temp(nridges))
!        ridges_temp(1:nridges)=ridges(1:nridges)
!        deallocate(ridges)
!        allocate(ridges(nridges))
!        ridges=ridges_temp
!        deallocate(ridges_temp)
!        end subroutine get_ridges
!
c---------------------------------------------------------------------        

        subroutine get_facet_ridges(thefacet,ridges)
        ! 
        ! gets the ridges (edges) of a facet. 
        ! Assumes that the facet has maximum D points (e.g., maximum 3
        ! points in 3D).
        !
        use defs, only : error_stop, fsubstart,fsubendext
        implicit none
        type(facet), intent(in) :: thefacet(:)
        type(ridge), allocatable :: ridges(:)
        ! local variables:
        integer, allocatable :: indices(:)
        integer i_ind,j_ind,k_ind
        integer nrp,nridges,nfacets,nfp,d,iridge,jridge,ipoint
        logical newridge

        !print fsubstart,'get_facet_ridges'
        !
        ! begin initialize
        !
        nfacets=size(thefacet)
        !print '(8x,I0," facets")',nfacets
        if (nfacets.ne.1)call error_stop('Please supply only one facet')
        d=size(thefacet(1)%thisfacet(1,:))
        !print '(8x,I0," dimensions")',d
        nfp=size(thefacet(1)%facetpoints) ! NUMBER OF POINTS OF THE FACET
        !print '(8x,I0," facetpoints")',nfp
        nrp=max(nfp-1,1) ! number of points per ridge
        !print '(8x,I0," pointsper ridge")',nrp
        nridges=nfp ! number of ridges
        !print '(8x,I0," ridges")',nridges
        if (allocated(ridges)) deallocate(ridges)
        allocate(ridges(1:nridges))
        !print '(8x,"ridges allocated")'
        do iridge=1,nridges
          if (allocated(ridges(iridge)%ridgepoints)) then
            deallocate(ridges(iridge)%ridgepoints)
            !print '(8x,"ridgepoints(iridge) deallocated")'
          end if
          allocate(ridges(iridge)%ridgepoints(nrp))
          allocate(ridges(iridge)%thisridge(nrp,d))
          ridges(iridge)%belongs_to_facets=0
          !print '(8x,"ridgepoints(iridge) allocated")'
        end do
        !print '(8x,"ridges%ridgepoints allocated")'
        !
        !
        allocate(indices(nrp))
        do i_ind=1,size(indices)
            indices(i_ind)=1
        end do ! i_ind
        i_ind=size(indices)
        !print '("indices initialized")'

        !print '("initialized...")'
        !
        ! end initialize
        !
        iridge=1
        do while (iridge<=nridges)
            !print '(8x,"iridge= ",I0)',iridge
            newridge=.true.
            do j_ind=1,size(indices)
              do k_ind=j_ind+1,size(indices)
                if (indices(j_ind)>=indices(k_ind)) newridge=.false.
              end do ! k_ind
            end do ! j_ind
            !print '(8x,"indices= ",10(I0,x))',indices
            !print '(8x,"newridge:",L1)',newridge
            if (newridge) then
              !print '(8x,"indices= ",10(I0,x))',indices
              !
              ! begin create ridge
              !
              do ipoint=1,nrp
                ridges(iridge)%ridgepoints(ipoint)=                       &
     &            thefacet(1)%facetpoints(indices(ipoint))
              end do ! ipoint
              !print '(8x,"ridgepoints assigned")'
              ridges(iridge)%belongs_to_facets(1)                         &
     &           =thefacet(1)%facetnumber
              !print '(8x,"facet assigned")'
              ridges(iridge)%on_horizon=.false.
              !print '(8x,"on_horizon initialized")'
              do j_ind=1,nrp
                ridges(iridge)%thisridge(j_ind,1:d)=                      &
     &            thefacet(1)%thisfacet(indices(j_ind),1:d)
              end do ! j_ind
              !print '(8x,"thisridge assigned")'
              !
              ! end create ridge
              !
              iridge=iridge+1
              !
            end if ! newridge
            !
            do while (all(indices(i_ind:nrp)==nfp.and.i_ind>1))
                i_ind=i_ind-1
            end do
            !
            do j_ind=i_ind+1,size(indices)
                indices(j_ind)=1
            end do
            !
            indices(i_ind)=indices(i_ind)+1
            i_ind=size(indices)
            !
        end do ! iridge <= nridges
        !
        !print fsubendext,'get_facet_ridges'
        !
        !
        end subroutine get_facet_ridges

c---------------------------------------------------------------------        

        subroutine get_all_ridges(facets,ridges)
        ! 
        ! gets the ridges (edges) of a set of facets. 
        ! Assumes that the facets all have the same number of points, 
        ! and that the facets have maximum D points (e.g., maximum 3
        ! points in 3D).
        !
        use defs, only : error_stop, fsubstart
        implicit none
        type(facet), intent(in) :: facets(:)
        type(ridge), allocatable, intent(out) :: ridges(:)
        ! local variables:
        type(ridge), allocatable :: facetridges(:),ridges_temp(:)
        integer nridges,nfacets,d,iridge,jridge,kridge,ifacet
        integer nfp,nrp
        logical newridge

        !print fsubstart,'get_all_ridges'
        !
        ! begin initialize
        !
        nfacets=size(facets)
        d=size(facets(1)%thisfacet(1,:))
        nfp=size(facets(1)%facetpoints) ! NUMBER OF POINTS OF EACH FACET
        nrp=max(nfp-1,1) ! number of points per ridge
        nridges=nfacets*nfp ! maximum possible number of ridges
        if (allocated(ridges)) deallocate(ridges)
        allocate(ridges(1:nridges))
        do iridge=1,nridges
            allocate(ridges(iridge)%ridgepoints(nrp))
        end do
        !print '("initialized...")'
        !
        ! end initialize
        !
        iridge=1
        do ifacet=1,nfacets
          call get_facet_ridges(facets(ifacet:ifacet),facetridges)
          do kridge=1,size(facetridges)
            newridge=.true.
            do jridge=1,iridge-1
                if (all(ridges(jridge)%ridgepoints                        &
     &                 ==facetridges(kridge)%ridgepoints)) then   
                  newridge=.false.
                  ridges(jridge)%belongs_to_facets(2)=ifacet
                end if  
            end do ! jridge
            !print*,"newridge:",newridge
            if (newridge) then
              ridges(iridge)=facetridges(kridge)
              !print*,"new ridge added."
              ridges(iridge)%belongs_to_facets(1)=ifacet
              iridge=iridge+1
              !print*,"counter raised."
            end if ! new ridge
          end do ! kridge
        end do ! ifacet
        !
        nridges=iridge-1
        allocate(ridges_temp(nridges))
        ridges_temp(1:nridges)=ridges(1:nridges)
        deallocate(ridges)
        allocate(ridges(nridges))
        ridges=ridges_temp
        deallocate(ridges_temp)
        !
        end subroutine get_all_ridges

c---------------------------------------------------------------------        
!
!        subroutine get_ridges_2(facets,ridges)
!        use defs, only : error_stop, fsubstart
!        implicit none
!        type(facet), intent(in) :: facets(:)
!        type(ridge), allocatable, intent(out) :: ridges(:)
!        ! local variables:
!        type(ridge), allocatable :: ridges_temp(:)
!        integer nridges,nfacets,d,iridge,ifacet,ipoint,jpoint
!        integer jridge
!        logical newridge
!
!        !print fsubstart,'get_ridges_2'
!        nfacets=size(facets)
!        d=size(facets(1)%thisfacet(1,:))
!        nridges=nfacets*d
!        if (allocated(ridges)) deallocate(ridges)
!        allocate(ridges(1:nridges))
!        do iridge=1,nridges
!          allocate(ridges(iridge)%ridgepoints(d-1))
!        end do
!        !print '("initialized...")'
!        iridge=1
!        select case(d)
!        case(2)
!          do ifacet=1,nfacets
!            do ipoint=1,size(facets(ifacet)%facetpoints)
!              newridge=.true.
!              do jridge=1,iridge-1
!!                print*,"iridge,ifacet,ipoint,jridge:",iridge,ifacet,      &
!!     &            ipoint,jridge
!                if (ridges(jridge)%ridgepoints(1)==                       &
!     &                facets(ifacet)%facetpoints(ipoint) ) then 
!                  newridge=.false.
!                  ridges(jridge)%belongs_to_facets(2)=ifacet
!                end if  
!              end do ! jridge
!              !print*,"newridge:",newridge
!              if (newridge) then
!                ridges(iridge)%ridgepoints(1:d-1)                         &
!     &              =facets(ifacet)%facetpoints(ipoint)
!                !print*,"new ridge added."
!                ridges(iridge)%belongs_to_facets(1)=ifacet
!                !print*,"facet assigned."
!                ! begin get ridge point coordinates
!                allocate(ridges(iridge)%thisridge(d-1,d))
!                ridges(iridge)%thisridge(1,:)=                          &
!     &          facets(ifacet)%thisfacet(ipoint,:)                
!                ! end get ridge point coordinates
!                iridge=iridge+1
!                !print*,"counter raised."
!              end if ! new ridge
!            end do ! ipoint
!          end do ! ifacet
!        case(3)
!          do ifacet=1,nfacets
!            do ipoint=1,size(facets(ifacet)%facetpoints)
!              do jpoint=ipoint+1,size(facets(ifacet)%facetpoints)
!                newridge=.true.
!                do jridge=1,iridge-1
!                  if (all(ridges(jridge)%ridgepoints(1:d-1)==           &
!     &                (/facets(ifacet)%facetpoints(ipoint),             &
!     &                  facets(ifacet)%facetpoints(jpoint) /))) then 
!                    newridge=.false.
!                    ridges(jridge)%belongs_to_facets(2)=ifacet
!                  end if  
!                end do ! jridge
!                if (newridge) then
!                  ridges(iridge)%ridgepoints(1:d-1)                     &
!     &              =(/facets(ifacet)%facetpoints(ipoint),              &
!     &                 facets(ifacet)%facetpoints(jpoint)/)
!                  ridges(iridge)%belongs_to_facets(1)=ifacet
!                  ! begin get ridge point coordinates
!                  allocate(ridges(iridge)%thisridge(d-1,d))
!                  ridges(iridge)%thisridge(1,:)=                          &
!     &            facets(ifacet)%thisfacet(ipoint,:)                
!                  ridges(iridge)%thisridge(2,:)=                          &
!     &            facets(ifacet)%thisfacet(jpoint,:)
!                  ! end get ridge point coordinates
!                  iridge=iridge+1
!                end if ! new ridge
!              end do ! jpoint
!            end do ! ipoint
!          end do ! ifacet
!        case default
!          call error_stop('dim not implemented')
!        end select
!        nridges=iridge-1
!        allocate(ridges_temp(nridges))
!        ridges_temp(1:nridges)=ridges(1:nridges)
!        deallocate(ridges)
!        allocate(ridges(nridges))
!        ridges=ridges_temp
!        deallocate(ridges_temp)
!        end subroutine get_ridges_2
!
c---------------------------------------------------------------------

        subroutine get_dfh(d,point,dir)
        !
        use defs, only : fsubstart,fsubendext, infty
        implicit none
        double precision :: norm2
        !
        integer, intent(in) :: d
        double precision, intent(in) :: point(:), dir(:)
        double precision dfh
        !
        ! local : 
        integer ifacet, nfacets,ipoint,iridge
        type(facet), allocatable :: facets(:)
        type(facet) closestfacet
        double precision ,allocatable :: nvec(:),crosspoint(:)
        double precision, allocatable :: midpoint(:)
        double precision :: dfh_old,crosspoint_dot_dir               
        double precision :: crosspoint_dot_dir_old,l
        type(ridge), allocatable :: ridges(:)
        double precision, allocatable :: theridge(:,:)
        logical is_inside
        character pointformat*256
        integer nrp
        !
        !
        print fsubstart, 'get_dfh'
        !
        pointformat=' '
        write(pointformat,'("(8x,",A,",",I0,"(F12.6,x))")'),              &
     &     "'point: '",d
        print pointformat, point
        !
        ! begin initialize
        !
        dfh=infty
        !allocate(nvec(d),crosspoint(d),midpoint(d),theridge(d,d-1)) ! OLD
        allocate(nvec(d),crosspoint(d),midpoint(d),theridge(d-1,d)) ! NEW
        crosspoint_dot_dir_old=infty
        dfh_old=infty
        !
        ! end initialize
        !
        !
        ! begin read facets from file HULL_FACETS
        !
        call read_facets(d,facets)
        nfacets=size(facets)
        !
        ! end read facets from file HULL_FACETS
        !
        !
        !print*,"point:",point
        !print*,"dir:",dir
        !print*," "
        do ifacet=1,nfacets
          !print*,"---------------------"
          !print*,"facetpoints:"
          !do ipoint=1,d
          !  print*,facets(ifacet)%thisfacet(ipoint,:)
          !end do ! ipoint
          !
          ! begin get nromal vector for facet
          !
          !call get_normal_vec(facets(ifacet)%thisfacet,nvec)
          call get_normal_vec_2(facets(ifacet)%thisfacet,nvec)
          !print*,"normal vector: unnormed",nvec
          !print*,"normal vector: normed",nvec/norm2(nvec)
          !
          ! end get nromal vector for facet
          !
          !
          ! begin find crossing point with facet along dir
          !
          !if (abs(dot_product(dir,nvec)).gt.0.0d0.or.                     &
          if (abs(dot_product(dir,nvec)).gt.1.0d-18.or.                   &
     &        norm2(dir).eq.0.0d0) then
            !
            call get_crossing_point(facets(ifacet)%thisfacet,nvec,point,  &
     &                              dir,crosspoint)
            !print*,"crosspoint:",crosspoint
            crosspoint_dot_dir=dot_product(crosspoint,dir)
            !
            ! end find crossing point with facet along dir
            !
            !
            !
            ! begin calculate midpoint of facet
            !
            midpoint=0.0d0
            do ipoint=1,d
              midpoint=midpoint+facets(ifacet)%thisfacet(ipoint,:)
            end do ! ipoint
            midpoint=midpoint/dble(d)
            !print*, "midpoint obtained"
            !call get_ridges_2(facets(ifacet:ifacet),ridges)
            call get_facet_ridges(facets(ifacet:ifacet),ridges)
            !print*, "ridges obtained"
            !
            ! end calculate midpoint of facet
            !
            !
            !
            ! begin check if crosspoin is inside facet
            !
            is_inside=.true.
            do iridge=1,d
              theridge=ridges(iridge)%thisridge
              !print*,"the ridge:",theridge
              !if (is_visible_2(theridge,crosspoint,midpoint)) then
              !nrp=size(theridge,1)
              !if (is_visible_3(theridge,crosspoint,midpoint)) then
              !if (is_visible_4(theridge,crosspoint,midpoint,-1.0D-8))     &
              if (is_visible_4(theridge,crosspoint,midpoint,-1.0D-12))     &
     &         then
                is_inside=.false.
              end if
            end do ! iridge
            !print*,"crosspoint inside:",is_inside
            !
            ! end check if crosspoin is inside facet
            !
            !
            if (is_inside) then
              !
              !l=norm2(crosspoint-point)
              !
              if (norm2(dir).gt.0.0d0) then
                l=dot_product(crosspoint-point,dir)/norm2(dir)
                if (crosspoint_dot_dir.lt.crosspoint_dot_dir_old) then
                  crosspoint_dot_dir_old=crosspoint_dot_dir
                  dfh=-l
                  dfh_old=dfh
                  closestfacet=facets(ifacet)
                end if 
              else
                l=norm2(crosspoint-point)
                if (l.lt.dfh_old) then
                  dfh_old=l
                  dfh=l
                  closestfacet=facets(ifacet)
                end if ! l < l_old
              end if ! | dir | > 0
              !print*,"l=",l
              !
            end if ! crosspoint inside
            !
          end if ! |(dir,nvec)| > 0 or | dir | == 0         
          !print*,"new dfh:",dfh
          !
        end do ! ifacet
        !
        print '(8x,"dfh=",F24.6)', dfh
        pointformat=' '
        write(pointformat,'("(8x,",A,",",I0,"(I0,x))")')                  &
     &     "'closest facet point numbers: '",d
        !print*, pointformat
        if(dfh.lt.infty) print pointformat,closestfacet%facetpoints(1:d)
        write(pointformat,'("(8x,",A,",",I0,"(A,x))")')                   &
     &     "'closest facet point labels: '",d
        if(dfh.lt.infty) print pointformat,                               &
     &     (adjustl(trim(closestfacet%labels(ipoint))),ipoint=1,d)
        !
        !
        print fsubendext, 'get_dfh'
        return
        !
        end subroutine get_dfh

c---------------------------------------------------------------------

        subroutine read_facets(d,facets)   
        use defs, only : error_stop, fsubstart,fsubendext
        implicit none
        !
        type(facet), allocatable, intent(inout) :: facets(:)     
        integer, intent(in) :: d
        ! local: 
        integer nfacets, ifacet, ipoint
        character line*1024
        !
        print fsubstart,'read_facets'
        !
        open(51,file='HULL_FACETS', status='old', err=101)
        nfacets=0
 10     read(51,'(A1024)', end=12) line
        line=trim(adjustl(line))
        if (line(1:1).ne.'#'.and.line(1:1).ne.'!'.and.line(1:1).ne.'%')   &
     &  then
          nfacets=nfacets+1
        end if   
        goto 10
        !
 12     continue
        print '(8x,I0," facets")', nfacets
        if (nfacets.le.0) goto 102
        rewind(51)
        allocate(facets(nfacets))
        facets(:)%visible=.false.
        facets(:)%noutside=0
        do ifacet=1, nfacets
          allocate(facets(ifacet)%facetpoints(d))
          allocate(facets(ifacet)%labels(d))
          allocate(facets(ifacet)%thisfacet(d,d))
        end do
        ifacet=0
        !print '(8x,"facets allocated")'
        !
 14     read(51,'(A1024)', end=16,err=101) line
        line=trim(adjustl(line))
        if (line(1:1).ne.'#'.and.line(1:1).ne.'!'.and.line(1:1).ne.'%')   &
     &  then
          ifacet=ifacet+1
          read(line,*,err=101,end=101)                                    &
     &                 (facets(ifacet)%thisfacet(ipoint,1:d),             &
     &                  ipoint=1,d),facets(ifacet)%facetpoints(1:d),      &
     &                  facets(ifacet)%labels(1:d)
          facets(ifacet)%facetnumber=ifacet
        end if   
        goto 14
        !
        !
 16     close(51)       
        print fsubendext,'read_facets'
        !
        return
 
 101    close(51)
        call error_stop('Error when reading HULL_FACETS')             
 102    close(51)
        call error_stop('No facets found in HULL_FACETS')             

        end subroutine read_facets

c---------------------------------------------------------------------

!        subroutine get_normal_vec(facet,nvec)
!        use defs, only : error_stop
!        use misc, only : cross_product
!        implicit none
!        double precision, intent(in) :: facet(:,:)
!        double precision, intent(inout) :: nvec(:)
!        ! local
!        integer d,ipoint
!        double precision, allocatable :: vec12(:),vec13(:)
!        !
!        d=size(nvec)
!        print*,"facet=",facet
!        !
!        select case(d)
!        case(2)
!          allocate(vec12(d))
!          vec12=facet(1,:)-facet(2,:)
!          nvec=(/vec12(2),-vec12(1)/) 
!          ! begin check if really normal
!!          print '(8x,"nvec*facet:",F15.8)',                               &
!!     &       dot_product(nvec,vec12)
!          ! end check if really normal
!        case(3)
!          allocate(vec12(d),vec13(d))
!          vec12=facet(1,:)-facet(2,:)
!          vec13=facet(1,:)-facet(3,:)
!          nvec=cross_product(vec12,vec13)
!          !print*,"vunitpar:",vec12/norm2(vec12)
!          !print*,"vunitpar:",vec13/norm2(vec13)
!!          print '(8x,"nvec*facet:",F15.8)',                               &
!!     &       dot_product(nvec,vec12)
!!          print '(8x,"nvec*facet:",F15.8)',                               &
!!     &       dot_product(nvec,vec13)
!        case default
!          call error_stop('dim not implemented')
!        end select
!        !
!
!        end subroutine get_normal_vec
!
c---------------------------------------------------------------------

        subroutine get_normal_vec_2(facet,nvec)
        use defs, only : error_stop
        use misc, only : init_random_seed 
        implicit none
        double precision, intent(in) :: facet(:,:)
        double precision, intent(inout) :: nvec(:)
        ! local
        integer d,ds,ipoint,nvpar,i,j,k
        double precision, allocatable :: vunitpar(:,:)
        !double precision, parameter :: tol=1D-8
        double precision, parameter :: tol=1D-12
        !
        ! begin get dimension and number of facet point
        d=size(facet,2)
        ds=size(facet,1)
        !allocate(nvec(d))
        !nvpar=ds*(ds-1)/2-1
        nvpar=min(d-1,ds*(ds-1)/2) 
!        print*,"d=",d
!        print*,"ds=",ds
!        print*,"nvpar=",nvpar
!        print*,"facet=",facet
        allocate(vunitpar(1:nvpar,1:d))
        ! end get dimension and number of facet points
        !
        ! begin set up unit vectors that span surface
        !
        k=1
        !
        do i=1,ds
            do j=i+1,ds
                if (k<=nvpar) then
                  vunitpar(k,:)=facet(i,:)-facet(j,:)
                  vunitpar(k,:)=vunitpar(k,:)/norm2(vunitpar(k,:))
!                  print*,"k,vunitpar(k):",k,vunitpar(k,:)                 &
!     &              /norm2(vunitpar(k,:))
                end if ! ky=nvpar
                k=k+1
            end do ! j
        end do ! i
        !
        ! begin orthogonalize vectors
        !
        do i=2,nvpar
            do j=1,i-1
                vunitpar(i,:)=vunitpar(i,:)                               &
     &           -dot_product(vunitpar(i,:),vunitpar(j,:))*vunitpar(j,:)
                vunitpar(i,:)=vunitpar(i,:)/norm2(vunitpar(i,:))
            end do ! j
        end do ! i
        ! end orthogonalize vectors
        !
        ! end set up unit vectors that span surface
        !
        nvec=0.0d0
        !
        do while (norm2(nvec).lt.tol)
            !
            ! begin set up trial vector
            !
            do i=1,size(nvec)
              ! get seed
              call init_random_seed()
              ! get random number
              call random_number(nvec(i))  
            end do
            !
            ! end set up trial vector

            ! begin subtract components in plane
            !
            do k=1,nvpar
              nvec(:)=nvec(:)-dot_product(nvec(:),vunitpar(k,:))          &
     &             *vunitpar(k,:)
            end do ! k
            !
        end do
        !
        if (norm2(nvec).ge.tol) nvec=nvec/norm2(nvec)
        !
        ! begin check if really normal
!        do k=1,nvpar
!          print '(8x,"nvec*facet:",F15.8)',                               &
!     &       dot_product(nvec,vunitpar(k,:))
!        end do ! k
        ! end check if really normal
        !
        end subroutine get_normal_vec_2

c---------------------------------------------------------------------

        subroutine get_crossing_point(facet,nvec,point,dir,crosspoint)
        implicit none
        double precision :: norm2
        !
        double precision, intent(in) :: facet(:,:), point(:), dir(:)
        double precision, intent(in) :: nvec(:)
        double precision, allocatable, intent(inout) :: crosspoint(:)
        ! local:
        integer d
        double precision l, cos_theta
        double precision, allocatable :: X(:) ! vertical crossing point
        !double precision, parameter :: tol=1D-8
        double precision, parameter :: tol=1D-12
        !
        d=size(point)
        allocate(X(d))
        !
        X=point+dot_product(facet(1,:)-point,nvec)*nvec/norm2(nvec)**2
        !
        cos_theta=1.0d0
        l=norm2(X-point)
        crosspoint=X
        if (norm2(dir).gt.tol.and.norm2(X-point).gt.tol) then
          cos_theta=dot_product(X-point,dir)/(norm2(X-point)*norm2(dir))
          l=norm2(X-point)/(norm2(dir)*cos_theta)
          crosspoint=point+l*dir
        end if

        end subroutine get_crossing_point

c---------------------------------------------------------------------

        double precision function fdet3(matrix) 
        implicit none
        double precision matrix(3,3)
        fdet3=0.0d0
        fdet3=fdet3+matrix(1,1)*matrix(2,2)*matrix(3,3)                   &
     &             +matrix(1,2)*matrix(2,3)*matrix(3,1)                   &
     &             +matrix(1,3)*matrix(2,1)*matrix(3,2)                   &
     &             -matrix(3,1)*matrix(2,2)*matrix(1,3)                   &
     &             -matrix(3,2)*matrix(2,3)*matrix(1,1)                   &
     &             -matrix(3,3)*matrix(2,1)*matrix(1,2)                 
        end function fdet3

c---------------------------------------------------------------------

        integer function idet3(matrix) 
        implicit none
        integer matrix(3,3)
        idet3=0
        idet3=idet3+matrix(1,1)*matrix(2,2)*matrix(3,3)                   &
     &             +matrix(1,2)*matrix(2,3)*matrix(3,1)                   &
     &             +matrix(1,3)*matrix(2,1)*matrix(3,2)                   &
     &             -matrix(3,1)*matrix(2,2)*matrix(1,3)                   &
     &             -matrix(3,2)*matrix(2,3)*matrix(1,1)                   &
     &             -matrix(3,3)*matrix(2,1)*matrix(1,2)                 
        end function idet3
  
c---------------------------------------------------------------------

        subroutine unitmat(n,outmat)
        ! create nxn unit matrix
        implicit none
        integer n,i
        double precision, allocatable :: outmat(:,:)
        allocate(outmat(n,n))
        outmat=0.0d0
        do i=1,n
          outmat(i,i)=1.0d0
        end do
        end subroutine unitmat
        
c---------------------------------------------------------------------

        subroutine get_mateigs(A,ALPHAR,printvecs)
        use defs, only : error_stop
        implicit none
        integer i
        logical printvecs
        CHARACTER          JOBVL, JOBVR
        INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N,M
*       ..
*       .. Array Arguments ..
        DOUBLE PRECISION ::  A(:, : )
        DOUBLE PRECISION, allocatable ::  ALPHAI( : ),                    &
     $                   ALPHAR( : ),
     $                   B( :, : ), BETA( : ), VL( :, : ),
     $                   VR( :, : ), WORK( : )
*      
        !
        N=size(A,1)
        M=size(A,2)
        if (.not.N==M) then
          call error_stop('Square matrices only.')
        end if
        LDA=max(N,1)
        LDB=max(N,1)
        LDVL=N
        LDVR=N
        LWORK=max(1,10*N)
        allocate(ALPHAR(N),ALPHAI(N),BETA(N),VL(LDVL,N),VR(LDVR,N))
        allocate(WORK(LWORK))
        JOBVL='N'
        !JOBVR='N'
        JOBVR='V'
        call unitmat(N,B)
        call DGGEV(JOBVL,JOBVR,N,A,LDA,B,LDB,ALPHAR,ALPHAI,               &
     &              BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
        if (printvecs) then
          print'(8x," ")'
          do i=1,N
              print *,VR(:,i)
          end do
          print'(8x," ")'
        end if
        !
        end subroutine get_mateigs
        !
c---------------------------------------------------------------------

        subroutine get_mateigs_full(A,ALPHAR,VR)
        use defs, only : error_stop
        implicit none
        integer i
        CHARACTER          JOBVL, JOBVR
        INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N,M
*       ..
*       .. Array Arguments ..
        DOUBLE PRECISION ::  A(:, : )
        DOUBLE PRECISION, allocatable ::  ALPHAI( : ),                    &
     $                   ALPHAR( : ),
     $                   B( :, : ), BETA( : ), VL( :, : ),
     $                   VR( :, : ), WORK( : )
*      
        !
        N=size(A,1)
        M=size(A,2)
        if (.not.N==M) then
          call error_stop('Square matrices only.')
        end if
        LDA=max(N,1)
        LDB=max(N,1)
        LDVL=N
        LDVR=N
        LWORK=max(1,10*N)
        allocate(ALPHAR(N),ALPHAI(N),BETA(N),VL(LDVL,N),VR(LDVR,N))
        allocate(WORK(LWORK))
        JOBVL='N'
        !JOBVR='N'
        JOBVR='V'
        call unitmat(N,B)
        call DGGEV(JOBVL,JOBVR,N,A,LDA,B,LDB,ALPHAR,ALPHAI,               &
     &              BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
        !
        ! begin normalize EV
        !
        do i=1,N
          VR(:,i)=VR(:,i)/norm2(VR(:,i))
        end do
        !
        ! end normalize EV
        !
        !print'(8x," ")'
        !do i=1,N
        !    print *,VR(:,i)
        !end do
        !print'(8x," ")'
        !
        end subroutine get_mateigs_full
        !
c---------------------------------------------------------------------

        subroutine get_mateigs_full_symm(A,ALPHAR,VR)
        use defs, only : error_stop
        implicit none
        integer i
        CHARACTER          JOBVL, JOBVR
        INTEGER            INFO, LDA, LDVR, LIWORK, LWORK, N,M
        INTEGER LWMAX
        PARAMETER(LWMAX=1000)  
*       ..
*       .. Array Arguments ..
        DOUBLE PRECISION ::  A(:, : )
        DOUBLE PRECISION, allocatable ::  ALPHAR( : ),  VR( :, : )
        double PRECISION WORK( LWMAX)
        INtegEr :: IWORK(1000)
*      
        !
        N=size(A,1)
        M=size(A,2)
        if (.not.N==M) then
          call error_stop('Square matrices only.')
        end if
        LDA=max(N,1)
        LDVR=N
        !LWORK=1+6*N+2*N**2
        LWORK=-1
        !LIWORK=3+5*N
        LIWORK=-1
        allocate(ALPHAR(N),VR(LDVR,N))
        !allocate(WORK(1000))
        call DSYEVD('V','U',N,A,LDA,ALPHAR, WORK, LWORK, IWORK,           &
     &              LIWORK, INFO)

        LWORK =  INT( WORK( 1 )  )
        LIWORK = IWORK( 1 ) 
        call DSYEVD('V','U',N,A,LDA,ALPHAR, WORK, LWORK, IWORK,           &
     &              LIWORK, INFO)
        !
        ! begin normalize EV
        !
        do i=1,N
          VR(:,i)=A(:,i)/norm2(A(:,i))
        end do
        !
        ! end normalize EV
        !
        end subroutine get_mateigs_full_symm
        !
c---------------------------------------------------------------------

        subroutine get_mateigs_full_cmplx(A,ALPHAR,ALPHAI,VR)
        use defs, only : error_stop
        implicit none
        integer i
        CHARACTER          JOBVL, JOBVR
        INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N,M
*       ..
*       .. Array Arguments ..
        DOUBLE PRECISION ::  A(:, : )
        DOUBLE PRECISION, allocatable ::  ALPHAI( : ),                    &
     $                   ALPHAR( : ),
     $                   B( :, : ), BETA( : ), VL( :, : ),
     $                   VR( :, : ), WORK( : )
*      
        !
        N=size(A,1)
        M=size(A,2)
        if (.not.N==M) then
          call error_stop('Square matrices only.')
        end if
        LDA=max(N,1)
        LDB=max(N,1)
        LDVL=N
        LDVR=N
        LWORK=max(1,10*N)
        allocate(ALPHAR(N),ALPHAI(N),BETA(N),VL(LDVL,N),VR(LDVR,N))
        allocate(WORK(LWORK))
        JOBVL='N'
        !JOBVR='N'
        JOBVR='V'
        call unitmat(N,B)
        call DGGEV(JOBVL,JOBVR,N,A,LDA,B,LDB,ALPHAR,ALPHAI,               &
     &              BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
        !
        ! begin normalize EV
        !
        do i=1,N
          VR(:,i)=VR(:,i)/norm2(VR(:,i))
        end do
        !
        ! end normalize EV
        !
        !print'(8x," ")'
        !do i=1,N
        !    print *,VR(:,i)
        !end do
        !print'(8x," ")'
        !
        end subroutine get_mateigs_full_cmplx
        !
!-------------------------------------------------------------
        !
        subroutine sort_mateigs(mateigs,matev)
        use defs
        implicit none
        double precision :: mateigs(:), matev(:,:)
        double precision,allocatable :: mateigdum, matevdum(:)
        double precision, allocatable:: unitvecs(:,:)
        double precision scalarprod0,scalarprod1
        integer matdim,i,j
        !
        !print fsubstart,'sort_mateigs'  
        matdim=size(matev,1)
        !print '(8x,"matrix dimension: ",I0)',matdim
        allocate(unitvecs(matdim,matdim))
        unitvecs=0.0d0
        do i=1,matdim
          unitvecs(i,i)=1.0d0
        end do
        !
        do i=1,matdim
          scalarprod0=abs(dot_product(matev(:,i),unitvecs(i,:)))
          do j=i+1,matdim
            scalarprod1=abs(dot_product(matev(:,j),unitvecs(i,:)))
            if (scalarprod1.gt.scalarprod0) then
              mateigdum=mateigs(i)
              mateigs(i)=mateigs(j)
              mateigs(j)=mateigdum
              matevdum=matev(:,i)
              matev(:,i)=matev(:,j)
              matev(:,j)=matevdum
              scalarprod0=abs(dot_product(matev(:,i),unitvecs(i,:)))
            end if
          end do ! j
          if (dot_product(matev(:,i),unitvecs(i,:)).lt.0.0d0) then
            matev(:,i)=-matev(:,i)
          end  if
        end do ! i
        !
        deallocate(unitvecs)
        !print fsubendext,'sort_mateigs'  
        !
        end subroutine sort_mateigs
!-------------------------------------------------------------
        !
        subroutine angle_axis_2_rotmatrix(angle,axis,matrix)
          !
          use defs
          implicit none
          double precision matrix(3,3)
          double precision, intent(in) :: angle
          double precision axis(1:3)
          ! local
          double precision anglerad
          integer i

        anglerad=angle*Pi/180.0d0
        ! biuld the rotation matrix
        matrix=0.0D0
        axis=axis/norm2(axis)
        matrix(1,1)=cos(anglerad)+axis(1)**2*(1.0d0-cos(anglerad))
        matrix(1,2)=axis(1)*axis(2)*(1.0d0-cos(anglerad))                 &
     &            -axis(3)*sin(anglerad)
        matrix(1,3)=axis(1)*axis(3)*(1.0d0-cos(anglerad))                 &
     &            +axis(2)*sin(anglerad)
        matrix(2,1)=axis(2)*axis(1)*(1.0d0-cos(anglerad))                 &
     &            +axis(3)*sin(anglerad)
        matrix(2,2)=cos(anglerad)+axis(2)**2*(1.0d0-cos(anglerad)) 
        matrix(2,3)=axis(2)*axis(3)*(1.0d0-cos(anglerad))                 &
     &            -axis(1)*sin(anglerad)
        matrix(3,1)=axis(3)*axis(1)*(1.0d0-cos(anglerad))                 &
     &            -axis(2)*sin(anglerad)
        matrix(3,2)=axis(3)*axis(2)*(1.0d0-cos(anglerad))                 &
     &            +axis(1)*sin(anglerad)
        matrix(3,3)=cos(anglerad)+axis(3)**2*(1.0d0-cos(anglerad)) 

        end subroutine angle_axis_2_rotmatrix

!-------------------------------------------------------------
        !
        subroutine angle_axis_2_mirrmatrix(angle,axis,matrix)
          !
          use defs
          implicit none
          double precision matrix(3,3)
          double precision, intent(in) :: angle
          double precision axis(1:3)
          ! local
          double precision anglerad
          integer i

        anglerad=angle*Pi/180.0d0
        ! biuld the mirror matrix
        matrix=0.0D0
        axis=axis/norm2(axis)
        matrix(1,1)=1.0d0-2.0d0*axis(1)**2
        matrix(1,2)=-2.0*axis(1)*axis(2)
        matrix(1,3)=-2.0d0*axis(1)*axis(3)
        matrix(2,1)=-2.0d0*axis(1)*axis(2)
        matrix(2,2)=1.0d0-2.0d0*axis(2)**2 
        matrix(2,3)=-2.0d0*axis(2)*axis(3)
        matrix(3,1)=-2.0d0*axis(1)*axis(3)
        matrix(3,2)=-2.0d0*axis(2)*axis(3)
        matrix(3,3)=1.0d0-2.0d0*axis(3)**2 

        end subroutine angle_axis_2_mirrmatrix

!--------------------------------------------------------------

        subroutine polyhedral_volume(ndim)
        use defs  
        use misc, only : cross_product
        implicit none
        integer, intent(in) :: ndim
        !character(len=*), intent (in) :: filename
        logical isopen,useunit
        integer nunit,nfacets,ifacet,npoints,ipoint
        character line*1024
        double precision, allocatable :: facets(:,:,:),origin(:)
        double precision, allocatable :: facetvec1(:),facetvec2(:)
        double precision, allocatable :: facetnormal(:)
        double precision vol,signum, facetarea
  
        if (talk) print fsubstart,'polyhedral_volume'  
        ! begin find free file unit
        nunit=51
        useunit=.false.
        do while ((.not.useunit).and.nunit.lt.100)
          INQUIRE (unit=nunit, opened=isopen) 
          if (.not.isopen) useunit=.true.
          nunit=nunit+1
        end do
        if (.not.useunit) call error_stop("no free file unit found")
        ! end find free file unit
        !
        open(nunit,file="HULL_FACETS",status="old")
        rewind(nunit)
        ! begin count facet points
        nfacets=0
10      read(nunit,'(A1024)',end=11) line
        line=adjustl(trim(line))
        if (line(1:1).eq."#") goto 10
        if (line(1:1).eq."!") goto 10
        nfacets=nfacets+1
        goto 10
11      rewind(nunit)
        if (nfacets.lt.ndim+1) call error_stop("too few facets found")  
        if (talk) print '(8x,I0," facets")',nfacets
        ! end count facet points
        !
        ! begin read facet points
        npoints=ndim
        allocate(facets(nfacets,npoints,ndim))
        allocate(origin(ndim))
        ifacet=1
        do while (ifacet.le.nfacets)
          read(nunit,'(A1024)',err=100) line
          line=adjustl(trim(line))
          do while(line(1:1).eq."#".or.line(1:1).eq."!") 
           read(nunit,'(A1024)',end=100,err=100) line
          end do
          read(line,*) (facets(ifacet,ipoint,1:ndim),ipoint=1,npoints)
          ifacet=ifacet+1
        end do
        close(nunit)
        if (talk) print '(8x,"facets read.")' 
        ! end read facet points
        !
        ! begin calculate origin (center of polyhedron)
        !
        origin=0.0d0
        do ifacet=1,nfacets
          do ipoint=1,npoints
            origin=origin + facets(ifacet,ipoint,:)
          end do ! ipoint
        end do ! ifacet
        origin=origin/dble(npoints*nfacets)
        if (talk) print '(8x,"polyhedral center is at ",120(F12.6))',     &
     &    origin
        !
        ! end calculate origin (center of polyhedron)
        !
        ! begin calculate volume (V=1/3 sum_F (df * P) A_F, F facet, df
        ! facet normal unit vector pointing outside, A_F facet area)
        ! https://en.wikipedia.org/wiki/Polyhedron#Volume
        !
        vol=0.0d0
        allocate(facetvec1(ndim),facetvec2(ndim),facetnormal(ndim))
        do ifacet=1,nfacets
          ! To do: this needs to be generalized to any dimension. 
          if (ndim==3) then
            facetvec1=facets(ifacet,2,:)-facets(ifacet,1,:)
            facetvec2=facets(ifacet,3,:)-facets(ifacet,1,:)
            facetarea=norm2(cross_product(facetvec1,facetvec2))/2.0d0
            facetnormal=cross_product(facetvec1,facetvec2)
            facetnormal=facetnormal/norm2(facetnormal)
            signum=dot_product(facetnormal,facets(ifacet,1,:)-origin)
            if (signum.lt.0.0d0) facetnormal=-facetnormal
          else 
            call error_stop('Only implemented for 3D')
          end if
          vol=vol+facetarea*dot_product(facets(ifacet,1,:)                &
     &      -origin,facetnormal)
        end do ! ifacet
        vol=vol/3.0d0
        print '(8x,"volume=",F24.6)',vol
        !
        ! end calculate volume
        !
        ! end normally
        !
        if (talk) print fsubendext,'polyhedral_volume'  
        return

100     nerr=nerr+1
        call error_stop("could not read facet points")        

        end subroutine polyhedral_volume        

!--------------------------------------------------------------

        subroutine matmul_fortran(n,A,B,C)
        ! 
        ! multiply matrix A with matrix B using fortran's matmul
        !
        use defs
        use misc, only : mymatmul
        implicit none 
        integer n,j
        double precision A(n,n), B(n,n), C(n,n)
        double precision D(n,n) ! DEBUG
        ! 
        if (talk) then
          print fsubstart, 'matmul_fortran'
          ! print matrix A 
          print '(8x," ")'
          print '(8x,"matrix dimension: ",I0)',n
          print '(8x," ")'
          print '(8x,"Matrix A:")'
          do j=1,n
            print '(8x,3(F10.6))',A(j,1:3)
          end do  
          print '(8x," ")'
          ! print matrix B
          print '(8x," ")'
          print '(8x,"matrix B:")'
          do j=1,n
            print '(8x,3(F10.6))',B(j,1:3)
          end do  
          print '(8x," ")'
        end if ! talk
        !
        C=matmul(A,B)
        call mymatmul(A,B,D)
        !
        if (talk) then
          print '(8x," ")'
          print '(8x,"fortran matmul: matrix C (C=matmul(AB)):")'
          do j=1,n
            print '(8x,3(F10.6))',C(j,1:3)
          end do  
          print '(8x," ")'
        end if ! talk
!        !
!        ! BEGIN DEBUG
!        !
!        if (talk) then
!          print '(8x," ")'
!          print '(8x,"homemade matmul: matrix D (D=AB):")'
!          do j=1,n
!            print '(8x,3(F10.6))',D(j,1:3)
!          end do  
!          print '(8x," ")'
!        end if ! talk
!        !
!        ! END DEBUG
!        !
        if (talk) then      
          print fsubendext, 'matmul_fortran'
        end if  
        !  
        end subroutine matmul_fortran

!-------------------------------------------------------------

        end module 
