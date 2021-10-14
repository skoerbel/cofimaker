        module madelung
        implicit none

        contains

c---------------------------------------------------------------------
!     
!      subroutine madelen(lattice,charge,epsil,alat,talk,nerr)
!      use defs
!      implicit none
!      ! calculates the madelung energy for some point charge lattices in
!      ! a homogeneous background charge  
!      character, intent(in) :: lattice*7 ! lattice type
!      double precision, intent(in) :: alat(1:3),charge,epsil ! lattice constant, point charge, permittivity
!      logical, intent(in) :: talk
!      integer :: nerr
!      ! internal variables
!      ! madelung constant, lattice vector modulus, c over a (tetra),
!      ! volume factor (tetra), wigner-seitz radius:
!      double precision madel,R,G,c_over_a,vol,r_ws
!      !double precision vm,dr
!      double precision, parameter :: beta=1.0D0
!      integer i,j,k,l,m,n,imax,maxl
!     
!      if(talk) then 
!        print fsubstart,'madelen'
!        print '("Calculating Madelung energy for ",A7," lattice with poi
!     &nt charge ",F8.4,"e, epsilon=" , F12.6,", and lattice constant",
!     &   3(F10.6)," Bohr" )',lattice,charge,epsil,alat
!      end if
!      select case(lattice)
!      case('CUB','cub','sc','SC')
!        vol=1.0d0
!        r_ws=(3.0d0/(4.0d0*pi)*vol)**(1.0d0/3.0d0) 
!        imax=500
!        maxl=500
!        madel=0.0D0
!        do i=-imax,imax
!        do j=-imax,imax
!        do k=-imax,imax
!              if (abs(i)+abs(j)+abs(k).gt.0) then
!                    R=sqrt(float(i**2+j**2+k**2))
!                    madel=madel+erfc(sqrt(beta)*R)/R
!              end if
!        end do
!        end do
!        end do
!        do l=-maxl,maxl
!        do m=-maxl,maxl
!        do n=-maxl,maxl
!              if (abs(l)+abs(m)+abs(n).gt.0) then
!                    G=2.0D0*pi*sqrt(float(l**2+m**2+n**2))
!                    madel=madel+4.0D0*pi*exp(-0.250D0*G**2/beta)/(G**2)
!              end if
!        end do
!        end do
!        end do
!        madel=madel-2.0D0*sqrt(beta/pi)-pi/beta
!        madel=-madel
!
!!        print*,"Madelung constant=",madel
!!        print '("Madelung energy in Ryd=",F20.10)',
!!     &     -madel*charge**2/(epsil*alat(1))
!!        print '("Madelung energy in eV=",F20.10)',
!!     &     -madel*charge**2*Ryd/(epsil*alat(1))
!
!        !open(10,file="erfc.dat",status="replace")
!        !do i=-imax,imax
!        !            write(10,*)float(i),erf(float(i)),erfc(float(i))
!        !end do
!        !close(10)
!
!        !vm=0.0D0
!        !do i=0,10000
!        !      dr=0.010D0
!        !      R=float(i)*dr
!        !      vm=vm+R*erfc(R)*dr
!        !end do
!        !print*, "V_av=",vm
!      case('fcc','FCC')
!        vol=0.25d0
!        r_ws=(3.0d0/(4.0d0*pi)*vol)**(1.0d0/3.0d0) 
!        imax=500
!        maxl=500
!        madel=0.0D0
!        do i=-imax,imax
!        do j=-imax,imax
!        do k=-imax,imax
!              if (abs(i)+abs(j)+abs(k).gt.0) then
!                    R=0.50D0*sqrt(float((i+j)**2+(j+k)**2+(i+k)**2))
!                    madel=madel+erfc(sqrt(beta)*R)/R
!              end if
!        end do
!        end do
!        end do
!        do l=-maxl,maxl
!        do m=-maxl,maxl
!        do n=-maxl,maxl
!              if (abs(l)+abs(m)+abs(n).gt.0) then
!                    G=2.0D0*pi*sqrt(float((l-m+n)**2+(-l+m+n)**2+(l+m-n)
!     &                **2))
!                    madel=madel+16.0D0*pi*exp(-0.250D0*G**2/beta)/(G**2)
!              end if
!        end do
!        end do
!        end do
!        madel=madel-2.0D0*sqrt(beta/pi)-4.0D0*pi/(beta)
!        madel=-madel
!
!!        print*,"Madelung constant=",madel
!!        print '("Madelung energy in Ryd=",F20.10)',
!!     &     -madel*charge**2/(epsil*alat(1))
!!        print '("Madelung energy in eV=",F20.10)',
!!     &     -madel*charge**2*Ryd/(epsil*alat(1))
!
!      case('bcc','BCC')
!        vol=0.5d0
!        r_ws=(3.0d0/(4.0d0*pi)*vol)**(1.0d0/3.0d0) 
!        imax=500
!        maxl=500
!        madel=0.0D0
!        do i=-imax,imax
!        do j=-imax,imax
!        do k=-imax,imax
!              if (abs(i)+abs(j)+abs(k).gt.0) then
!                    R=0.50D0*sqrt(float((i+j-k)**2+(i-j+k)**2+(-i+j+k)
!     &                **2))
!                    madel=madel+erfc(sqrt(beta)*R)/R
!              end if
!        end do
!        end do
!        end do
!        do l=-maxl,maxl
!        do m=-maxl,maxl
!        do n=-maxl,maxl
!              if (abs(l)+abs(m)+abs(n).gt.0) then
!                    G=2.0D0*pi*sqrt(float((l+m)**2+(l+n)**2+(m+n)
!     &                **2))
!                    madel=madel+8.0D0*pi*exp(-0.250D0*G**2/beta)/(G**2)
!              end if
!        end do
!        end do
!        end do
!        madel=madel-2.0D0*sqrt(beta/pi)-2.0D0*pi/(beta)
!        madel=-madel
!
!!        print*,"Madelung constant=",madel
!!        print '("Madelung energy in Ryd=",F20.10)',
!!     &     -madel*charge**2/(epsil*alat(1))
!!        print '("Madelung energy in eV=",F20.10)',
!!     &     -madel*charge**2*Ryd/(epsil*alat(1))
!
!      case('tet','TET','tetra','TETRA')
!        c_over_a=alat(3)/(alat(1))
!        vol=c_over_a
!        r_ws=(3.0d0/(4.0d0*pi)*vol)**(1.0d0/3.0d0) 
!        imax=500
!        maxl=500
!        madel=0.0D0
!
!        do i=-imax,imax
!        do j=-imax,imax
!        do k=-imax,imax
!              if (abs(i)+abs(j)+abs(k).gt.0) then
!                    R=sqrt(float(i**2+j**2)+c_over_a**2*float(k**2))
!                    madel=madel+erfc(sqrt(beta)*R)/R
!              end if
!        end do
!        end do
!        end do
!        do l=-maxl,maxl
!        do m=-maxl,maxl
!        do n=-maxl,maxl
!              if (abs(l)+abs(m)+abs(n).gt.0) then
!                    G=2.0D0*pi*sqrt(float(l**2+m**2)+float(n**2)
!     &                /c_over_a**2)
!                    madel=madel+4.0D0*pi*exp(-0.250D0*G**2/beta)
!     &                    /(vol*G**2)
!              end if
!        end do
!        end do
!        end do
!        madel=madel-2.0D0*sqrt(beta/pi)-pi/(beta*vol)
!        madel=-madel
!!        print*,"Madelung constant=",madel
!!        print '("Madelung energy in Ryd=",F20.10)',
!!     &     -madel*charge**2/(epsil*alat(1))
!!        print '("Madelung energy in eV=",F20.10)',
!!     &     -madel*charge**2*Ryd/(epsil*alat(1))
!      case default
!        goto 100
!      end select
!
!      print*,"Madelung constant (a)=",madel
!      print*,"Madelung constant (r_ws)=",madel*r_ws
!      print '("Madelung energy in Ryd=",F20.10)',
!     &   -madel*charge**2/(epsil*alat(1))
!      print '("Madelung energy in eV=",F20.10)',
!     &   -madel*charge**2*Ryd/(epsil*alat(1))
!
!      if(talk) then 
!        print fsubendext,'madelen'
!      end if
!      return
!! errors
!100   nerr=nerr+1
!      print ferrmssg, "Sorry, lattice type not implemented"
!      return
!
!      end subroutine

!c---------------------------------------------------------------------
!     
!      subroutine madelen1(cell,charge,epsil,talk,nerr)
!      use defs
!      use misc
!      implicit none
!      ! calculates the madelung energy for some point charge lattices in
!      ! a homogeneous background charge  
!      double precision, intent(in) :: cell(1:3,1:3),charge,
!     &                                epsil(1:3,1:3) ! lattice vectors, point charge, epsilon tensor
!      logical, intent(in) :: talk
!      integer :: nerr
!      ! internal variables
!      ! madelung constant, lattice vector modulus, c over a (tetra),
!      ! volume factor (tetra), wigner-seitz radius, inverse epsilon, det
!      ! (epsilon):
!      double precision madel,rvec(1:3),rabs,gvecs(1:3,1:3),gvec(1:3),
!     &        gabs,vol,r_ws,epsil1(1:3,1:3),det_eps
!      double precision, parameter :: beta=1.0D0 ! sqrt(Ewald parameter)
!      integer i,j,k,l,m,n,imax,maxl
!     
!      if(talk) then 
!        print fsubstart,'madelen1'
!        print '("Calculating Madelung energy for a lattice of point char
!     &ges ",F8.4,"e,"/"epsilon="/,3(3(F12.6)/),
!     & ", and lattice constant"/,
!     &   3(3(F10.6)/)," Bohr" )',charge,epsil,cell(1:3,1:3)
!      end if
!
!      ! calculate volume, reciprocal lattice vecs 
!      vol=abs(dot_product(cell(1,:),cross_product(cell(2,:),cell(3,:))))
!      print '("volume=",F12.6)',vol
!      call invert_matrix(epsil,epsil1)
!      print '("epsilon^-1="/,3(3(F12.6)/))',epsil1
!      det_eps=abs(dot_product(epsil(:,1),
!     &     cross_product(epsil(:,2),epsil(:,3))))
!      print '("det(epsilon)="/,F25.6)', det_eps
!      r_ws=(3.0d0/(4.0d0*pi)*vol)**(1.0d0/3.0d0) 
!      print '("r_ws=",F12.6)',r_ws
!      gvecs(1,1:3)=2.0d0*pi*cross_product(cell(2,:),cell(3,:))/vol  
!      gvecs(2,1:3)=2.0d0*pi*cross_product(cell(3,:),cell(1,:))/vol  
!      gvecs(3,1:3)=2.0d0*pi*cross_product(cell(1,:),cell(2,:))/vol  
!     
!      ! limits of direct and reciprocal sum 
!      imax=50
!      maxl=50
!
!      madel=0.0D0
!      ! sum over real lattice vectors
!      do i=-imax,imax
!        do j=-imax,imax
!          do k=-imax,imax
!            if (abs(i)+abs(j)+abs(k).gt.0) then
!              rvec=dble(i)*cell(1,:)+dble(j)*cell(2,:)+dble(k)*cell(3,:)
!              rabs=sqrt(dot_product(rvec,rvec))
!              madel=madel+erfc(sqrt(beta)*rabs)/rabs
!            end if
!          end do
!        end do
!      end do
!      ! sum over reciprocal lattice vectors
!      do l=-maxl,maxl
!        do m=-maxl,maxl
!          do n=-maxl,maxl
!            if (abs(l)+abs(m)+abs(n).gt.0) then
!              gvec=dble(l)*gvecs(1,:)+dble(m)*gvecs(2,:)
!     &            +dble(n)*gvecs(3,:)
!              gabs=sqrt(dot_product(gvec,gvec))
!              madel=madel+4.0D0*pi/vol
!     &             *exp(-0.250D0*gabs**2/beta)/(gabs**2)
!            end if
!          end do
!        end do
!      end do
!      madel=madel-2.0D0*sqrt(beta/pi)-pi/(beta*vol)
!      madel=-madel
!
!      !print*,"Madelung constant (a)=",madel*vol**(1.0d0/3.0d0)
!      print*,"Madelung constant (r_ws)=",madel*r_ws
!      print '("Madelung energy in Ryd=",F20.10)',
!     &   -madel*charge**2/(2.0d0*epsil(1,1))
!      print '("Madelung energy in eV=",F20.10)',
!     &   -madel*charge**2*Ryd/(2.0d0*epsil(1,1))
!
!      if(talk) then 
!        print fsubendext,'madelen1'
!      end if
!      return
!! errors
!100   nerr=nerr+1
!      return
!
!      end subroutine

c---------------------------------------------------------------------
     
      subroutine madelen1(cell,charge,epsil)
      use defs
      use misc
      implicit none
      ! calculates the madelung energy for some point charge lattices in
      ! a homogeneous background charge  
      double precision, intent(in) :: cell(1:3,1:3),charge,
     &                                epsil(1:3,1:3) ! lattice vectors, point charge, epsilon tensor
      ! internal variables
      ! madelung constant, lattice vector modulus, 
      ! volume, wigner-seitz radius, inverse epsilon, det
      ! (epsilon):
      double precision madel,rvec(1:3),rabs,gvecs(1:3,1:3),gvec(1:3),
     &        gabs,vol,r_ws,epsil1(1:3,1:3),det_eps
      double precision r_epsil1(1:3),sqrt_r_epsil1_r,epsil_g(1:3),
     &        g_epsil_g
      double precision, parameter :: beta=1.0D0 ! sqrt(Ewald parameter)
      integer i,j,k,l,m,n,imax,maxl
     
      if(talk) then 
        print fsubstart,'madelen1'
        print '("Calculating Madelung energy for a lattice of point char
     &ges ",F8.4,"e,"/"epsilon="/,3(3(F12.6)/),
     & ", and lattice constant"/,
     &   3(3(F10.6)/)," Bohr" )',charge,epsil,cell(1:3,1:3)
      end if

      ! calculate volume, reciprocal lattice vecs 
      vol=abs(dot_product(cell(1,:),cross_product(cell(2,:),cell(3,:))))
      print '("volume=",F12.6)',vol
      call invert_matrix(epsil,epsil1)
      print '("epsilon^-1="/,3(3(F12.6)/))',epsil1
      det_eps=abs(dot_product(epsil(:,1),
     &     cross_product(epsil(:,2),epsil(:,3))))
      print '("det(epsilon)="/,F25.6)', det_eps
      r_ws=(3.0d0/(4.0d0*pi)*vol)**(1.0d0/3.0d0) 
      print '("r_ws=",F12.6)',r_ws
      gvecs(1,1:3)=2.0d0*pi*cross_product(cell(2,:),cell(3,:))/vol  
      gvecs(2,1:3)=2.0d0*pi*cross_product(cell(3,:),cell(1,:))/vol  
      gvecs(3,1:3)=2.0d0*pi*cross_product(cell(1,:),cell(2,:))/vol  
     
      ! limits of direct and reciprocal sum 
      imax=50
      maxl=50

      madel=0.0D0
      ! sum over real lattice vectors
      do i=-imax,imax
        do j=-imax,imax
          do k=-imax,imax
            if (abs(i)+abs(j)+abs(k).gt.0) then
              rvec=dble(i)*cell(1,:)+dble(j)*cell(2,:)+dble(k)*cell(3,:)
              r_epsil1=matmul(epsil1,rvec)
              sqrt_r_epsil1_r=sqrt(dot_product(r_epsil1,rvec))
              madel=madel+erfc(sqrt(beta)*sqrt_r_epsil1_r)/
     &             (sqrt(det_eps)*sqrt_r_epsil1_r)
            end if
          end do
        end do
      end do
      ! sum over reciprocal lattice vectors
      do l=-maxl,maxl
        do m=-maxl,maxl
          do n=-maxl,maxl
            if (abs(l)+abs(m)+abs(n).gt.0) then
              gvec=dble(l)*gvecs(1,:)+dble(m)*gvecs(2,:)
     &            +dble(n)*gvecs(3,:)
              epsil_g=matmul(epsil,gvec)
              g_epsil_g=dot_product(gvec,epsil_g)
              madel=madel+4.0D0*pi/vol
     &             *exp(-0.250D0*g_epsil_g/beta)/(g_epsil_g)
            end if
          end do
        end do
      end do
      madel=madel-2.0D0*sqrt(beta/(pi*det_eps))-pi/(beta*vol)

      print*,"Madelung constant (r_ws)=",-madel*r_ws
      print '("Madelung energy in Ryd=",F20.10)',
     &   madel*charge**2/(1.0d0)
      print '("Madelung energy in eV=",F20.10)',
     &   madel*charge**2*Ryd/(1.0d0)

      if(talk) then 
        print fsubendext,'madelen1'
      end if
      return
! errors
100   nerr=nerr+1
      return

      end subroutine madelen1

c---------------------------------------------------------------------
     
      subroutine dipen(cell,pvec,epsil)
      use defs
      use misc
      implicit none
      ! calculates the dipole-dipole energy for a lattice of point
      ! dipoles with the formula of Kantorovic, PRB 60, 15477 
      double precision, intent(in) :: cell(1:3,1:3),pvec(1:3), epsil ! lattice vectors, dipole moment, epsilon (scalar)
      ! internal variables
      ! madelung constant, lattice vector modulus, 
      ! volume, wigner-seitz radius, matrix "psi" from Kantotovic's
      ! paper :
      double precision energy,rvec(1:3),gvecs(1:3,1:3),gvec(1:3),
     &        gabs,vol,r_ws,psimat(1:3,1:3)
      double precision p_psi_p
      double precision, parameter :: beta=1.0D0 ! Ewald parameter
      integer i,j,k,l,m,n,imax,maxl,imat,jmat
     
      if(talk) then 
        print fsubstart,'dipen'
        print '("Calculating dipole-dipole energy of point dipoles with 
     &dipole moment ",3(F12.6)," e Bohr,"/"epsilon=",F12.6,
     & ", and cell vectors in Bohr"/,
     &   3(3(F10.6)/))',pvec(1:3),epsil,cell(1,1:3),cell(2,:),cell(3,:)
      end if
      
      ! limits of direct and reciprocal sum 
      imax=50
      maxl=50

      ! calculate volume, reciprocal lattice vecs , r_ws
      vol=abs(dot_product(cell(1,:),cross_product(cell(2,:),cell(3,:))))
      if (talk) print '("volume=",F12.6)',vol
      r_ws=(3.0d0/(4.0d0*pi)*vol)**(1.0d0/3.0d0) 
      if (talk) print '("r_ws=",F12.6)',r_ws
      gvecs(1,1:3)=2.0d0*pi*cross_product(cell(2,:),cell(3,:))/vol  
      gvecs(2,1:3)=2.0d0*pi*cross_product(cell(3,:),cell(1,:))/vol  
      gvecs(3,1:3)=2.0d0*pi*cross_product(cell(1,:),cell(2,:))/vol  
     
      psimat=0.0D0
      ! sum over real lattice vectors
      do i=-imax,imax
        do j=-imax,imax
          do k=-imax,imax
            if (abs(i)+abs(j)+abs(k).gt.0) then
              rvec=dble(i)*cell(1,:)+dble(j)*cell(2,:)+dble(k)*cell(3,:)
              psimat=psimat-beta**3*Hmat(beta*rvec)
            end if
          end do
        end do
      end do
      ! sum over reciprocal lattice vectors
      do l=-maxl,maxl
        do m=-maxl,maxl
          do n=-maxl,maxl
            if (abs(l)+abs(m)+abs(n).gt.0) then
              gvec=dble(l)*gvecs(1,:)+dble(m)*gvecs(2,:)
     &            +dble(n)*gvecs(3,:)
              gabs=sqrt(dot_product(gvec,gvec))
              do imat=1,3
                do jmat=1,3
                  psimat(imat,jmat)=psimat(imat,jmat)
     &                 + (4.0d0*pi/vol)*gvec(imat)*gvec(jmat)/gabs**2
     &                  *exp(-gabs**2/(4.0d0*beta**2))
                end do
              end do
            end if
          end do
        end do
      end do
      
      p_psi_p=dot_product(pvec,matmul(psimat,pvec))
!      energy= - 2.0d0 * beta**3 / ( 3.0d0 * sqrt(pi) ) 
!     &        * dot_product(pvec,pvec)  ! Makov-Payne term
!     &        + 0.5d0 * p_psi_p 
      ! factor two is unclear, maybe comes from Bohr/Hartree factor?
      energy= - 4.0d0 * beta**3 / ( 3.0d0 * sqrt(pi) ) 
     &        * dot_product(pvec,pvec)  ! Makov-Payne term
     &        + 1.0d0 * p_psi_p 
      energy=energy/epsil

      print '("Dipole-dipole energy in Ryd=",F20.10)',
     &   energy
      print '("Dipole-dipole energy in eV=",F20.10)',
     &   energy*Ryd
      print '("Makov-Payne dipole-dipole energy in Ryd=",F20.10)',
     &  - 2.0d0 * pi / ( 3.0d0 * vol ) 
     &        * dot_product(pvec,pvec) / epsil  ! Makov-Payne term
      print '("Makov-Payne dipole-dipole energy in eV=",F20.10)',
     &  - 2.0d0 * pi / ( 3.0d0 * vol ) 
     &        * dot_product(pvec,pvec) / epsil * Ryd  ! Makov-Payne term

      if(talk) then 
        print fsubendext,'dipen'
      end if
      return
! errors
100   nerr=nerr+1
      return

      contains 

      function Hmat(y)
      implicit none
      double precision, dimension (1:3,1:3) :: Hmat
      double precision y(1:3)
      double precision y_y
      integer i,j
      y_y=dot_product(y,y)
      Hmat=0.0d0
      do i=1,3
        Hmat(i,i)=-1.0d0*H(y)
        do j=1,3
          Hmat(i,j)=Hmat(i,j)+y(i)*y(j)/y_y
     &              *(3.0d0*H(y)+4.0d0/sqrt(pi)*exp(-y_y))
        end do
      end do
      end function

      double precision function H(y)
      implicit none
      double precision y(1:3)
      double precision y_y,yabs
      y_y=dot_product(y,y)
      yabs=sqrt(y_y)
      H=2.0d0/sqrt(pi)/y_y*exp(-y_y)+erfc(yabs)/yabs**3
      end function 

      end subroutine

c---------------------------------------------------------------------
     
      subroutine madeldipen2(cell,charge,pvec,epsil)
      use defs
      use misc
      implicit none
      ! calculates the madelung energy for point charge lattices whose
      ! basis is a dipole   
      double precision, intent(in) :: cell(1:3,1:3),charge,pvec(1:3),
     &                                epsil(1:3,1:3) ! lattice vectors, point charge, dipole moment/charge (the distance vector), epsilon tensor
      ! internal variables
      ! madelung constant, lattice vector modulus, 
      ! volume, wigner-seitz radius, inverse epsilon, det
      ! (epsilon):
      double precision madel,rvec(1:3),rabs,gvecs(1:3,1:3),gvec(1:3),
     &        gabs,vol,r_ws,epsil1(1:3,1:3),det_eps,pabs
      double precision r_epsil1(1:3),sqrt_r_epsil1_r,epsil_g(1:3),
     &        g_epsil_g,g_p
      double precision, parameter :: beta=1.0D0 ! sqrt(Ewald parameter)
      integer i,j,k,l,m,n,imax,maxl
     
      if(talk) then 
        print fsubstart,'madeldipen2'
        print '("Calculating Madelung energy for a lattice of point char
     &ges ",F8.4,"e,"/"epsilon="/,3(3(F12.6)/),
     & ", and lattice constant"/,
     &   3(3(F10.6)/)," Bohr" )',charge,epsil,cell(1:3,1:3)
      end if

      ! limits of direct and reciprocal sum 
      imax=50
      maxl=50

      ! calculate volume, reciprocal lattice vecs 
      vol=abs(dot_product(cell(1,:),cross_product(cell(2,:),cell(3,:))))
      print '("volume=",F12.6)',vol
      call invert_matrix(epsil,epsil1)
      print '("epsilon^-1="/,3(3(F12.6)/))',epsil1
      det_eps=abs(dot_product(epsil(:,1),
     &     cross_product(epsil(:,2),epsil(:,3))))
      print '("det(epsilon)="/,F25.6)', det_eps
      r_ws=(3.0d0/(4.0d0*pi)*vol)**(1.0d0/3.0d0) 
      print '("r_ws=",F12.6)',r_ws
      gvecs(1,1:3)=2.0d0*pi*cross_product(cell(2,:),cell(3,:))/vol  
      gvecs(2,1:3)=2.0d0*pi*cross_product(cell(3,:),cell(1,:))/vol  
      gvecs(3,1:3)=2.0d0*pi*cross_product(cell(1,:),cell(2,:))/vol  
      pabs=sqrt(dot_product(pvec,pvec))
     
      madel=0.0D0
      ! sum over real lattice vectors
      do i=-imax,imax
        do j=-imax,imax
          do k=-imax,imax
            if (abs(i)+abs(j)+abs(k).gt.0) then
              ! positive charge at rvec :
              rvec=dble(i)*cell(1,:)+dble(j)*cell(2,:)+dble(k)*cell(3,:)
              r_epsil1=matmul(epsil1,rvec)
              sqrt_r_epsil1_r=sqrt(dot_product(r_epsil1,rvec))
              madel=madel+erfc(sqrt(beta)*sqrt_r_epsil1_r)/
     &             (sqrt(det_eps)*sqrt_r_epsil1_r)
              ! now the negative one at rvec + pvec  
              rvec=rvec+pvec 
              r_epsil1=matmul(epsil1,rvec)
              sqrt_r_epsil1_r=sqrt(dot_product(r_epsil1,rvec))
              madel=madel-erfc(sqrt(beta)*sqrt_r_epsil1_r)/
     &             (sqrt(det_eps)*sqrt_r_epsil1_r)
            end if
          end do
        end do
      end do
      ! self-interaction of dipole (R=0 term in direct lattice sum):
      rvec=pvec 
      r_epsil1=matmul(epsil1,rvec)
      sqrt_r_epsil1_r=sqrt(dot_product(r_epsil1,rvec))
      madel=madel-erfc(sqrt(beta)*sqrt_r_epsil1_r)/
     &     (sqrt(det_eps)*sqrt_r_epsil1_r)

      ! sum over reciprocal lattice vectors
      do l=-maxl,maxl
        do m=-maxl,maxl
          do n=-maxl,maxl
            if (abs(l)+abs(m)+abs(n).gt.0) then
              ! positive charge at rvec:      
              gvec=dble(l)*gvecs(1,:)+dble(m)*gvecs(2,:)
     &            +dble(n)*gvecs(3,:)
              epsil_g=matmul(epsil,gvec)
              g_epsil_g=dot_product(gvec,epsil_g)
              madel=madel+4.0D0*pi/vol
     &             *exp(-0.250D0*g_epsil_g/beta)/(g_epsil_g)
              ! negative charge at rvec + pvec :
              g_p=dot_product(gvec,pvec)    
              !madel=madel-4.0D0*pi*exp(-0.250D0*gabs**2/beta)
!     &                    *cos(G_r)/(vol*gabs**2)
              madel=madel-(4.0D0*pi/vol)
     &             *exp(-0.250D0*g_epsil_g/beta)
     &             * cos (g_p) /(g_epsil_g)
            end if
          end do
        end do
      end do
      madel=madel-2.0D0*sqrt(beta/(pi*det_eps))
      
      print*,"Madelung constant (r_nn)=",-madel*pabs
      print '("Madelung energy in Ryd=",F20.10)',
     &   2.0d0*madel*charge**2/(1.0d0)
      print '("Madelung energy in eV=",F20.10)',
     &   2.0d0*madel*charge**2*Ryd/(1.0d0)
      
      ! for the dipole-dipole energy, subtract the dipole self energy
      rvec=pvec 
      r_epsil1=matmul(epsil1,rvec)
      sqrt_r_epsil1_r=sqrt(dot_product(r_epsil1,rvec))
      madel=madel + 1.0d0/(sqrt(det_eps)*sqrt_r_epsil1_r)
      
      print '("self-energy of one dipole in Ryd=",F20.10)',
     &   -2.0d0*charge**2/(sqrt(det_eps)*sqrt_r_epsil1_r)
      print '("self-energy of one dipole in eV=",F20.10)',
     &   -2.0d0*charge**2/(sqrt(det_eps)*sqrt_r_epsil1_r) * Ryd
      
      print '("dipole-dipole energy in Ryd=",F20.10)',
     &   (2.0d0*madel*charge**2)/(1.0d0)
      print '("dipole-dipole energy in eV=",F20.10)',
     &   (2.0d0*madel*charge**2)*Ryd/(1.0d0)

      if(talk) then 
        print fsubendext,'madeldipen2'
      end if
      return
! errors
100   nerr=nerr+1
      return

      end subroutine

c---------------------------------------------------------------------
     
      subroutine madelen(cell,atoms,epsil,madel)
      use defs
      use misc
      implicit none
      ! calculates the madelung energy for any point charge lattice 
      ! in a medium with any dielectric tensor using Eq. F.5 in Richard
      ! Martin, Electronic structure  
      double precision, intent(inout) :: cell(1:3,1:3),
     &                                epsil(1:3,1:3) ! lattice vectors, epsilon tensor
      type(atom), intent(inout) :: atoms(:) 
      ! internal variables
      ! madelung constant, lattice vector modulus, 
      ! volume, wigner-seitz radius, inverse epsilon, det
      ! (epsilon):
      double precision madel,rvec(1:3),rabs,gvecs(1:3,1:3),gvec(1:3),
     &        gabs,vol,r_ws,epsil1(1:3,1:3),det_eps,g_r,madeld,madelr
      double precision r_epsil1(1:3),sqrt_r_epsil1_r,epsil_g(1:3),
     &        g_epsil_g,sum_q,sum_q2
      double precision :: beta=10.0D0 ! sqrt(Ewald parameter)
      integer i,j,k,l,m,n,imax,maxl
      integer iatom,jatom,natoms

      if(talk) then 
        print fsubstart,'madelen'
      end if
      !
      ! number of atoms
      natoms=size(atoms)
      !
      ! get absolute positions
      do iatom=1,natoms
        atoms(iatom)%abswhere=atoms(iatom)%where(1)*cell(1,1:3)
     &               +atoms(iatom)%where(2)*cell(2,1:3)
     &               +atoms(iatom)%where(3)*cell(3,1:3)
      end do
      !
      ! convert lengths to Bohr
      cell=cell/bohr
      do iatom=1,natoms
        atoms(iatom)%abswhere(:)=atoms(iatom)%abswhere(:)/bohr
      end do
      !
      if(talk) then 
        print '(8x,"Calculating Madelung energy for a lattice of point c
     &harges with epsilon="/,3(8x,3(F12.6)/),8x,
     & " and lattice constant"/,
     &   3(8x,3(F10.6)/),8x," Bohr" )',epsil,cell(1:3,1:3)
        do iatom=1,natoms
          print '(8x,A2,1x,3(F8.3,1x),F7.3)',atoms(iatom)%name(1:2),
     &           atoms(iatom)%abswhere,atoms(iatom)%charge
        end do
      end if
      !
      ! calculate sum of charges and squared charges:
      sum_q=0.0d0
      sum_q2=0.0d0
      do iatom=1,natoms
        sum_q=sum_q+atoms(iatom)%charge
        sum_q2=sum_q2+atoms(iatom)%charge**2
      end do
      if(talk)print '(8x,"sum of charges=",F12.6)',sum_q
      if(talk)print '(8x,"sum of squared charges=",F12.6)',sum_q2

      ! limits of direct and reciprocal sum 
      imax=40
      maxl=40

      ! calculate volume, reciprocal lattice vecs 
      vol=abs(dot_product(cell(1,:),cross_product(cell(2,:),cell(3,:))))
      if(talk) print '(8x,"volume=",F12.6)',vol
      call invert_matrix(epsil,epsil1)
      if(talk) print '(8x,"epsilon^-1="/,3(8x,3(F12.6)/))',epsil1
      det_eps=abs(dot_product(epsil(:,1),
     &     cross_product(epsil(:,2),epsil(:,3))))
      if(talk) print '(8x,"det(epsilon)="/,8x,F25.6)', det_eps
      r_ws=(3.0d0/(4.0d0*pi)*vol)**(1.0d0/3.0d0) 
      if(talk) print '(8x,"r_ws=",F12.6)',r_ws
      ! reset beta to beta/epsilon in order to make the Gaussian decay
      ! independent of epsilon:
      beta=beta/(det_eps**(1.0d0/3.0d0))

      gvecs(1,1:3)=2.0d0*pi*cross_product(cell(2,:),cell(3,:))/vol  
      gvecs(2,1:3)=2.0d0*pi*cross_product(cell(3,:),cell(1,:))/vol  
      gvecs(3,1:3)=2.0d0*pi*cross_product(cell(1,:),cell(2,:))/vol  
     
      ! sum over real lattice vectors
      madeld=0.0D0
      ! sum over atoms:
      do iatom=1,natoms
       do jatom=1,natoms
        do i=-imax,imax
          do j=-imax,imax
            do k=-imax,imax
                ! ion at lattice vec + distance from origin of cell  
                rvec=dble(i)*cell(1,:)+dble(j)*cell(2,:)
     &              +dble(k)*cell(3,:)
                rvec=rvec+atoms(jatom)%abswhere-atoms(iatom)%abswhere
                r_epsil1=matmul(epsil1,rvec)
                sqrt_r_epsil1_r=sqrt(dot_product(r_epsil1,rvec))
                rabs=sqrt(rvec(1)**2+rvec(2)**2+rvec(3)**2)
                ! only sum over different atoms, or the same atom in
                ! different cells:
                if(rabs.gt.1.0E-5) then
                  madeld=madeld+atoms(iatom)%charge*atoms(jatom)%charge
     &                 *erfc(sqrt(beta)*sqrt_r_epsil1_r)/
     &                 (sqrt(det_eps)*sqrt_r_epsil1_r)
                end if   
!              end if
            end do ! k
          end do ! j
        end do ! i
       end do ! jatom
      end do  ! iatom
      ! correct for double-counting in the sum over atoms:
      madeld=madeld/2.0d0
      madeld=madeld-sum_q2*sqrt(beta/(pi*det_eps))
      if(talk) print '(8x,"Direct part: ",F10.4," eV")',                  &
     &                 madeld*Ryd*2.0d0

      ! sum over reciprocal lattice vectors
      madelr=0.0d0
      do iatom=1,natoms
       do jatom=1,natoms
        do l=-maxl,maxl
          do m=-maxl,maxl
            do n=-maxl,maxl
              if (abs(l)+abs(m)+abs(n).gt.0) then
!                ! positive charge at rvec:      
                gvec=dble(l)*gvecs(1,:)+dble(m)*gvecs(2,:)
     &              +dble(n)*gvecs(3,:)
                epsil_g=matmul(epsil,gvec)
                g_epsil_g=dot_product(gvec,epsil_g)
                ! charge at rvec :
                rvec=atoms(jatom)%abswhere-atoms(iatom)%abswhere
                g_r=dot_product(gvec,rvec)    
                madelr=madelr+atoms(iatom)%charge*atoms(jatom)%charge
     &               *(4.0D0*pi/vol)
     &               *exp(-0.250D0*g_epsil_g/beta)
     &               * cos (g_r) /(g_epsil_g)
              end if
            end do !n
          end do ! m
        end do ! l
       end do ! jatom   
      end do ! iatom   
      !
      ! correct for double-counting in the sum over atoms:
      madelr=madelr/2.0d0
      madelr=madelr-pi*sum_q**2/(2.0d0*beta*vol)
      if(talk) print '(8x,"Reciprocal part: ",F10.4," eV")',            &
     &               madelr*Ryd*2.0d0
      !  
      madel=madeld+madelr
      madel=madel*2.0d0
      if(talk) print '(8x,"Madelung energy in Ryd=",F20.10)',           &
     &   madel
      if(talk) print '(8x,"Madelung energy in eV=",F20.10)',            &
     &   madel*Ryd
      if(talk) print '(8x,"Madelung constant (r_ws)=",F20.10)',         &
     &                -madel*r_ws/2.0d0
      !
      ! convert lengths back to Angs
      cell=cell*bohr
      do iatom=1,natoms
        atoms(iatom)%abswhere(:)=atoms(iatom)%abswhere(:)*bohr
      end do
      !
      if(talk) then 
        print fsubendext,'madelen'
      end if
      return
! errors
100   nerr=nerr+1
      return

      end subroutine madelen

c---------------------------------------------------------------------
     
      subroutine madeldipen(lattice,charge,pvec,epsil,alat)
      use defs
      implicit none
      ! calculates the madelung energy for dipole lattices in
      ! an isotropic dielectric   
      character, intent(in) :: lattice*7 ! lattice type
      double precision, intent(in) :: alat(1:3),charge,epsil,pvec(1:3) ! lattice constant in Bohr, point charge, permittivity, dipole vector (p/q) in Bohr
      ! internal variables
      double precision madel,rvec(1:3),pabs,rabs,gvec(1:3),gabs,
     &                 c_over_a,vol,G_r
      ! double precision vm,dr
      double precision, parameter :: beta=1.0d0 !beta=1.0D0
      integer i,j,k,l,m,n,imax,maxl
     
      if(talk) then 
        print fsubstart,'madeldipen'
        print '("Calculating Madelung energy for ",A7," lattice with poi
     &nt dipole charges +-",F8.4,"e, epsilon=" ,F12.6/,
     &", and lattice constant",
     &   3(F10.6)," Bohr" )',lattice,charge,epsil,alat
      end if

      pabs=sqrt(dot_product(pvec,pvec))

      select case(lattice)
      case('CUB','cub','sc','SC')
        vol=1.0d0 ! in units of alat
        imax=50
        maxl=50
        madel=0.0D0
        do i=-imax,imax
        do j=-imax,imax
        do k=-imax,imax
              if (abs(i)+abs(j)+abs(k).gt.0) then
                    ! first the positive charge  
                    rvec(1)=dble(i)
                    rvec(2)=dble(j)
                    rvec(3)=dble(k)
                    rabs=sqrt(dot_product(rvec,rvec))
                    madel=madel+erfc(sqrt(beta)*rabs)/rabs
                    ! now the negative one at rvec + pvec  
                    rvec=rvec+pvec/alat(1) 
                    rabs=sqrt(dot_product(rvec,rvec))
                    madel=madel-erfc(sqrt(beta)*rabs)/rabs
              end if
        end do
        end do
        end do
        do l=-maxl,maxl
        do m=-maxl,maxl
        do n=-maxl,maxl
              if (abs(l)+abs(m)+abs(n).gt.0) then
                    gvec(1)=2.0d0*pi*dble(l)
                    gvec(2)=2.0d0*pi*dble(m)
                    gvec(3)=2.0d0*pi*dble(n)
                    gabs=sqrt(dot_product(gvec,gvec))
                    ! positive charge
                    madel=madel+4.0D0*pi*exp(-0.250D0*gabs**2/beta)
     &                    /(vol*gabs**2)
                    ! negative charge 
                    G_r=dot_product(gvec,pvec)/alat(1)    
                    madel=madel-4.0D0*pi*exp(-0.250D0*gabs**2/beta)
     &                    *cos(G_r)/(vol*gabs**2)
              end if
        end do
        end do
        end do
        ! positive charge 
        madel=madel-2.0D0*sqrt(beta/pi) !-pi/beta
        ! negative charge 
        !madel=madel+2.0D0*sqrt(beta/pi)+pi/beta
        ! the contribution of the negative charge at pvec (R=0):
        madel=madel-1.0d0/(1.0d0*pabs/alat(1))
        madel=-madel

        print*,"Madelung constant=",madel
        print '("Madelung energy in Ryd=",F20.10)',
     &     -madel*charge**2/(2.0d0*epsil*alat(1))
        print '("Madelung energy in eV=",F20.10)',
     &     -madel*charge**2*Ryd/(2.0d0*epsil*alat(1))
        
        print '("Energy per dipole in Ryd=",F20.10)',
     &     (-1.0d0*madel+1.0d0/(pabs/alat(1)))*charge**2/(epsil*alat(1))
        print '("Madelung energy in eV=",F20.10)',
     &     (-1.0d0*madel+1.0d0/(pabs/alat(1)))*charge**2*Ryd
     &     /(epsil*alat(1))

      case('fcc','FCC')
        !vol=0.25d0
        vol=1.0d0
        imax=50
        maxl=50
        madel=0.0D0
        do i=-imax,imax
        do j=-imax,imax
        do k=-imax,imax
              if (abs(i)+abs(j)+abs(k).gt.0) then
                    ! first the positive charge  
                    rvec(1)=0.5d0*dble(j+k)
                    rvec(2)=0.5d0*dble(i+k)
                    rvec(3)=0.5d0*dble(i+j)
                    rabs=sqrt(dot_product(rvec,rvec))
                    madel=madel+erfc(sqrt(beta)*rabs)/rabs
                    ! now the negative one at rvec + pvec  
                    rvec=rvec+pvec/alat(1) 
                    rabs=sqrt(dot_product(rvec,rvec))
                    madel=madel-erfc(sqrt(beta)*rabs)/rabs
              end if
        end do
        end do
        end do
        do l=-maxl,maxl
        do m=-maxl,maxl
        do n=-maxl,maxl
              if (abs(l)+abs(m)+abs(n).gt.0) then
                    gvec(1)=1.0d0*pi*dble(-l+m+n)
                    gvec(2)=1.0d0*pi*dble(l-m+n)
                    gvec(3)=1.0d0*pi*dble(l+m-n)
                    gabs=sqrt(dot_product(gvec,gvec))
                    ! positive charge
                    madel=madel+4.0D0*pi*exp(-0.250D0*gabs**2/beta)
     &                    /(vol*gabs**2)
                    ! negative charge
                    G_r=dot_product(gvec,pvec)/alat(1)    
                    madel=madel-4.0D0*pi*exp(-0.250D0*gabs**2/beta)
     &                    *cos(G_r)/(vol*gabs**2)
              end if
        end do
        end do
        end do
        madel=madel-2.0D0*sqrt(beta/pi) !-4.0D0*pi/(beta)
        !madel=madel-1.0d0/(1.0d0*pabs/alat(1))
        madel=madel-erfc(sqrt(beta)*pabs/alat(1))/(pabs/alat(1))
        madel=-madel

        print*,"Madelung constant=",madel
        print '("Madelung energy in Ryd=",F20.10)',
     &     -madel*charge**2/(2.0d0*epsil*alat(1))
        print '("Madelung energy in eV=",F20.10)',
     &     -madel*charge**2*Ryd/(2.0d0*epsil*alat(1))
        
        print '("Energy per dipole in Ryd=",F20.10)',
     &     (-1.0d0*madel+1.0d0/(pabs/alat(1)))*charge**2/(epsil*alat(1))
        print '("Energy per dipole in eV=",F20.10)',
     &     (-1.0d0*madel+1.0d0/(pabs/alat(1)))*charge**2*Ryd
     &     /(epsil*alat(1))

      case('bcc','BCC')
        vol=0.5d0
        imax=50
        maxl=50
        madel=0.0D0
        do i=-imax,imax
        do j=-imax,imax
        do k=-imax,imax
              if (abs(i)+abs(j)+abs(k).gt.0) then
                    ! first the positive charge  
                    rvec(1)=0.5d0*dble(-i+j+k)
                    rvec(2)=0.5d0*dble(i-j+k)
                    rvec(3)=0.5d0*dble(i+j-k)
                    rabs=sqrt(dot_product(rvec,rvec))
                    madel=madel+erfc(sqrt(beta)*rabs)/rabs
                    ! now the negative one at rvec + pvec  
                    rvec=rvec+pvec/alat(1) 
                    rabs=sqrt(dot_product(rvec,rvec))
                    madel=madel-erfc(sqrt(beta)*rabs)/rabs
              end if
        end do
        end do
        end do
        do l=-maxl,maxl
        do m=-maxl,maxl
        do n=-maxl,maxl
              if (abs(l)+abs(m)+abs(n).gt.0) then
                    gvec(1)=1.0d0*pi*dble(m+n)
                    gvec(2)=1.0d0*pi*dble(l+n)
                    gvec(3)=1.0d0*pi*dble(l+m)
                    gabs=sqrt(dot_product(gvec,gvec))
                    ! positive charge
                    madel=madel+4.0D0*pi*exp(-0.250D0*gabs**2/beta)
     &                    /(vol*gabs**2)
                    ! negative charge
                    G_r=dot_product(gvec,pvec)/alat(1)    
                    madel=madel-4.0D0*pi*exp(-0.250D0*gabs**2/beta)
     &                    *cos(G_r)/(vol*gabs**2)
              end if
        end do
        end do
        end do
        madel=madel-2.0D0*sqrt(beta/pi) !-2.0D0*pi/(beta)
        madel=madel-1.0d0/(1.0d0*pabs/alat(1))
        madel=-madel

        print*,"Madelung constant=",madel
        print '("Madelung energy in Ryd=",F20.10)',
     &     -madel*charge**2/(2.0d0*epsil*alat(1))
        print '("Madelung energy in eV=",F20.10)',
     &     -madel*charge**2*Ryd/(2.0d0*epsil*alat(1))
        
        print '("Energy per dipole in Ryd=",F20.10)',
     &     (-1.0d0*madel+1.0d0/(pabs/alat(1)))*charge**2/(epsil*alat(1))
        print '("Madelung energy in eV=",F20.10)',
     &     (-1.0d0*madel+1.0d0/(pabs/alat(1)))*charge**2*Ryd
     &     /(epsil*alat(1))

      case('tet','TET','tetra','TETRA')
        c_over_a=alat(3)/(alat(1))
        vol=c_over_a
        imax=50
        maxl=50
        madel=0.0D0

        do i=-imax,imax
        do j=-imax,imax
        do k=-imax,imax
              if (abs(i)+abs(j)+abs(k).gt.0) then
                    ! first the positive charge  
                    rvec(1)=dble(i)
                    rvec(2)=dble(j)
                    rvec(3)=dble(k)*c_over_a
                    rabs=sqrt(dot_product(rvec,rvec))
                    madel=madel+erfc(sqrt(beta)*rabs)/rabs
                    ! now the negative one at rvec + pvec  
                    rvec=rvec+pvec/alat(1) 
                    rabs=sqrt(dot_product(rvec,rvec))
                    madel=madel-erfc(sqrt(beta)*rabs)/rabs
              end if
        end do
        end do
        end do
        do l=-maxl,maxl
        do m=-maxl,maxl
        do n=-maxl,maxl
              if (abs(l)+abs(m)+abs(n).gt.0) then
                    gvec(1)=2.0d0*pi*dble(l)
                    gvec(2)=2.0d0*pi*dble(m)
                    gvec(3)=2.0d0*pi*dble(n)/c_over_a
                    gabs=sqrt(dot_product(gvec,gvec))
                    ! positive charge
                    madel=madel+4.0D0*pi*exp(-0.250D0*gabs**2/beta)
     &                   /(vol*gabs*2)
                    ! negative charge
                    G_r=dot_product(gvec,pvec)/alat(1)    
                    madel=madel-4.0D0*pi*exp(-0.250D0*gabs**2/beta)
     &                    *cos(G_r)/(vol*gabs**2)
              end if
        end do
        end do
        end do
        madel=madel-2.0D0*sqrt(beta/pi) !-pi/(beta*vol)
        madel=madel-1.0d0/(1.0d0*pabs/alat(1)) ! contribution of R=0,tau=pvec
        !madel=madel-erfc(sqrt(beta)*pabs/alat(1))/(pabs/alat(1)) ! contribution of R=0,tau=pvec
        madel=-madel
        print*,"Madelung constant=",madel
        print '("Madelung energy of the lattice in Ryd=",F20.10)',
     &     -madel*charge**2/(2.0d0*epsil*alat(1))
        print '("Madelung energy of the lattice in eV=",F20.10)',
     &     -madel*charge**2*Ryd/(2.0d0*epsil*alat(1))
        print '("Energy of one dipole in eV=",F20.10)',
     &     -charge**2*Ryd/(2.0d0*pabs*epsil)
        print '("Dipole-dipole energy per dipole in Ryd=",F20.10)',
     &     (-1.0d0*madel+1.0d0/(pabs/alat(1)))*charge**2/(epsil*alat(1))
        print '("Dipole-dipole energy per dipole in eV=",F20.10)',
     &     (-1.0d0*madel+1.0d0/(pabs/alat(1)))*charge**2*Ryd
     &     /(epsil*alat(1))

      case default
        goto 100
      end select
      if(talk) then 
        print fsubendext,'madeldipolen'
      end if
      return
! errors
100   nerr=nerr+1
      print ferrmssg, "Sorry, lattice type not implemented"
      return

      end subroutine

c---------------------------------------------------------------------
!     
!      subroutine dipoleen0(lattice,charge,dipole,epsil,alat,talk,nerr)
!      use defs
!      implicit none
!      ! calculates the dipole-dipole energy for some point charge lattices in
!      ! a homogeneous background charge  
!      ! *** WRONG *** only direct sum, depends on shape of box. Do not
!      ! use. 
!      character, intent(in) :: lattice*7 ! lattice type
!      double precision, intent(in) :: alat(1:3),charge,dipole(1:3),epsil! lattice constant, point charge, dipole vector (d/q), permittivity
!      logical, intent(in) :: talk
!      integer :: nerr
!      ! internal variables
!      double precision dipen,R,G,c_over_a,vol,rvec(1:3),dip,prod ! dipole-dipole energy, direct and recip. lattice vector length, ... , ... , lattice vector, |dipole|, scalar product dipole*rvec
!      double precision, parameter :: beta=1.0D0
!      integer i,j,k,imax
!     
!      if(talk) then 
!        print fsubstart,'dipoleen0'
!        print '("Calculating dipole-dipole energy for ",A7," lattice wit
!     &h point charge ",F8.4,"e, dipole vector (d/q)=",3(F12.6), " Bohr",
!     &   ",epsilon=" , F12.6,", and lattice constant",
!     &   3(F10.6)," Bohr" )',lattice,charge,dipole,epsil,alat
!      end if
!      select case(lattice)
!      case('CUB','cub','sc','SC')
!        imax=500
!        dipen=0.0D0
!        do i=-imax,imax
!        do j=-imax,imax
!        do k=-imax,imax
!              if (abs(i)+abs(j)+abs(k).gt.0) then
!                    R=sqrt(dble(i**2+j**2+k**2))
!                    dip=sqrt(dot_product(dipole,dipole))
!                    rvec(1)=dble(i)
!                    rvec(2)=dble(j)
!                    rvec(3)=dble(k)
!                    prod=3.0d0*dot_product(dipole,rvec)**2
!     &              /(dip**2 * R**2)
!                    dipen=dipen+(1.0d0-prod)/R**3
!              end if
!        end do
!        end do
!        end do
!        print*,"Dipole-dipole sum=",dipen
!        print '("Dipole-dipole energy in Ryd=",F20.10)',
!     &     dipen*dip**2*charge**2 / (epsil * alat(1)**3 )
!        print '("Dipole-dipole energy in eV=",F20.10)',
!     &     dipen*dip**2*charge**2 * Ryd / (epsil * alat(1)**3 )
!
!      case('fcc','FCC')
!        imax=500
!        dipen=0.0D0
!        do i=-imax,imax
!        do j=-imax,imax
!        do k=-imax,imax
!              if (abs(i)+abs(j)+abs(k).gt.0) then
!                    R=0.50D0*sqrt(float((i+j)**2+(j+k)**2+(i+k)**2))
!                    dip=sqrt(dot_product(dipole,dipole))
!                    rvec(1)=0.5d0*dble(j+k)
!                    rvec(2)=0.5d0*dble(i+k)
!                    rvec(3)=0.5d0*dble(i+j)
!                    prod=3.0d0*dot_product(dipole,rvec)**2
!     &              /(dip**2 * R**2)
!                    dipen=dipen+(1.0d0-prod)/R**3
!              end if
!        end do
!        end do
!        end do
!        print*,"Dipole-dipole sum=",dipen
!        print '("Dipole-dipole energy in Ryd=",F20.10)',
!     &     dipen*dip**2*charge**2 / (epsil * alat(1)**3 )
!        print '("Dipole-dipole energy in eV=",F20.10)',
!     &     dipen*dip**2*charge**2 * Ryd / (epsil * alat(1)**3 )
!
!      case('bcc','BCC')
!        imax=500
!        dipen=0.0D0
!        do i=-imax,imax
!        do j=-imax,imax
!        do k=-imax,imax
!              if (abs(i)+abs(j)+abs(k).gt.0) then
!                    R=0.50D0*sqrt(float((i+j-k)**2+(i-j+k)**2+(-i+j+k)
!     &                **2))
!                    dip=sqrt(dot_product(dipole,dipole))
!                    rvec(1)=0.5d0*dble(-i+j+k)
!                    rvec(2)=0.5d0*dble(i-j+k)
!                    rvec(3)=0.5d0*dble(i+j-k)
!                    prod=3.0d0*dot_product(dipole,rvec)**2
!     &              /(dip**2 * R**2)
!                    dipen=dipen+(1.0d0-prod)/R**3
!              end if
!        end do
!        end do
!        end do
!        print*,"Dipole-dipole sum=",dipen
!        print '("Dipole-dipole energy in Ryd=",F20.10)',
!     &     dipen*dip**2*charge**2 / (epsil * alat(1)**3 )
!        print '("Dipole-dipole energy in eV=",F20.10)',
!     &     dipen*dip**2*charge**2 * Ryd / (epsil * alat(1)**3 )
!
!      case('tet','TET','tetra','TETRA')
!        c_over_a=alat(3)/alat(1)
!        vol=c_over_a
!        imax=500
!        dipen=0.0D0
!        do i=-imax,imax
!        do j=-imax,imax
!        do k=-imax,imax
!              if (abs(i)+abs(j)+abs(k).gt.0) then
!                    R=sqrt(dble(i**2+j**2)+c_over_a**2*dble(k**2))
!                    dip=sqrt(dot_product(dipole,dipole))
!                    rvec(1)=dble(i)
!                    rvec(2)=dble(j)
!                    rvec(3)=c_over_a*dble(k)
!                    prod=3.0d0*dot_product(dipole,rvec)**2
!     &              /(dip**2 * R**2)
!                    dipen=dipen+(1.0d0-prod)/R**3
!              end if
!        end do
!        end do
!        end do
!        print*,"Dipole-dipole sum=",dipen
!        print '("Dipole-dipole energy in Ryd=",F20.10)',
!     &     dipen*dip**2*charge**2 / (epsil * alat(1)**3 )
!        print '("Dipole-dipole energy in eV=",F20.10)',
!     &     dipen*dip**2*charge**2 * Ryd / (epsil * alat(1)**3 )
!      case default
!        goto 100
!      end select
!      if(talk) then 
!        print fsubendext,'dipoleen0'
!      end if
!      return
!! errors
!100   nerr=nerr+1
!      print ferrmssg, "Sorry, lattice type not implemented"
!      return
!
!      end subroutine
!

c---------------------------------------------------------------------
!     
!      subroutine dipoleen(lattice,charge,pvec,epsil,alat,talk,nerr)
!      use defs
!      use misc ! invert_matrix
!      implicit none
!      ! calculates the dipole-dipole energy for some point charge lattices in
!      ! a homogeneous background charge  
!      ! *** WRONG *** only direct sum, depends on shape of summation
!      ! box. Do not use!
!      character, intent(in) :: lattice*7 ! lattice type
!      double precision, intent(in) :: alat(1:3),charge,pvec(1:3),
!     &                                epsil(1:3,1:3) ! lattice constant, point charge, dipole vector (d/q), permittivity
!      logical, intent(in) :: talk
!      integer :: nerr
!      ! internal variables
!      ! dipole-dipole energy, direct and recip. lattice vector length, ... , ... , lattice vector, |dipole|, scalar product dipole*rvec:
!      double precision dipen,rabs,G,c_over_a,vol,rvec(1:3),pabs
!      double precision epsil1(1:3,1:3) ! inverse epsilon tensor 
!      double precision p_epsil1_p ! dipole moment times inverse epsilon tensor times dipole moment 
!      double precision p_epsil1_r ! dipole moment times inverse epsilon tensor times lattice vector 
!      double precision p_r ! dipole moment times lattice vector 
!      double precision p_epsil1(1:3) ! dipole moment times inverse epsilon tensor  
!      !double precision, parameter :: beta=1.0D0 ! Ewald sum parameter
!      integer i,j,k,imax
!     
!      ! print the input parameters
!      if(talk) then 
!        print fsubstart,'dipoleen'
!        print '("Calculating dipole-dipole energy for ",A7," lattice wit
!     &h point charge ",F8.4,"e, dipole vector (d/q)=",3(F12.6), " Bohr",
!     &   //,",epsilon=" , 9(F12.6),", and lattice constant",
!     &   3(F10.6)," Bohr" )',lattice,charge,pvec,epsil,alat
!      end if
!      
!      ! invert epsilon
!      call invert_matrix(epsil,epsil1)
!      ! some abbreviations
!      p_epsil1=matmul(pvec,epsil1)
!      pabs=sqrt(dot_product(pvec,pvec))
!      p_epsil1_p=dot_product(p_epsil1,pvec)/pabs**2
!      
!      select case(lattice)
!      case('CUB','cub','sc','SC')
!        imax=500
!        dipen=0.0D0
!        do i=-imax,imax
!        do j=-imax,imax
!        do k=-imax,imax
!              if (abs(i)+abs(j)+abs(k).gt.0) then
!                    rvec(1)=dble(i)
!                    rvec(2)=dble(j)
!                    rvec(3)=dble(k)
!                    rabs=sqrt(dot_product(rvec,rvec))
!                    p_r=dot_product(pvec,rvec)/(pabs*rabs)
!                    p_epsil1_r=dot_product(p_epsil1,rvec)/(pabs * rabs)
!                    dipen=dipen+(p_epsil1_p-3.0d0*p_epsil1_r*p_r)
!     &                    /rabs**3
!              end if
!        end do
!        end do
!        end do
!        print*,"Dipole-dipole sum=",dipen
!        print '("Dipole-dipole energy in Ryd=",F20.10)',
!     &     dipen*pabs**2*charge**2 / alat(1)**3 
!        print '("Dipole-dipole energy in eV=",F20.10)',
!     &     dipen*pabs**2*charge**2 * Ryd /  alat(1)**3 
!
!      case('fcc','FCC')
!        imax=500
!        dipen=0.0D0
!        do i=-imax,imax
!        do j=-imax,imax
!        do k=-imax,imax
!              if (abs(i)+abs(j)+abs(k).gt.0) then
!                    rvec(1)=0.5d0*dble(j+k)
!                    rvec(2)=0.5d0*dble(i+k)
!                    rvec(3)=0.5d0*dble(i+j)
!                    rabs=sqrt(dot_product(rvec,rvec))
!                    ! test sum over sphere
!                    !if(rabs.lt.0.50d0*dble(imax)) then
!                     p_r=dot_product(pvec,rvec)/(pabs*rabs)
!                     p_epsil1_r=dot_product(p_epsil1,rvec)/(pabs * rabs)
!                     dipen=dipen+(p_epsil1_p-3.0d0*p_epsil1_r*p_r)
!     &                    /rabs**3
!                    !end if    
!              end if
!        end do
!        end do
!        end do
!        print*,"Dipole-dipole sum=",dipen
!        print '("Dipole-dipole energy in Ryd=",F20.10)',
!     &     dipen*pabs**2*charge**2 / alat(1)**3
!        print '("Dipole-dipole energy in eV=",F20.10)',
!     &     dipen*pabs**2*charge**2 * Ryd / alat(1)**3 
!
!      case('bcc','BCC')
!        imax=500
!        dipen=0.0D0
!        do i=-imax,imax
!        do j=-imax,imax
!        do k=-imax,imax
!              if (abs(i)+abs(j)+abs(k).gt.0) then
!                    rvec(1)=0.5d0*dble(-i+j+k)
!                    rvec(2)=0.5d0*dble(i-j+k)
!                    rvec(3)=0.5d0*dble(i+j-k)
!                    rabs=sqrt(dot_product(rvec,rvec))
!                    p_r=dot_product(pvec,rvec)/(pabs*rabs)
!                    p_epsil1_r=dot_product(p_epsil1,rvec)/(pabs * rabs)
!                    dipen=dipen+(p_epsil1_p-3.0d0*p_epsil1_r*p_r)
!     &                    /rabs**3
!              end if
!        end do
!        end do
!        end do
!        print*,"Dipole-dipole sum=",dipen
!        print '("Dipole-dipole energy in Ryd=",F20.10)',
!     &     dipen*pabs**2*charge**2 / alat(1)**3 
!        print '("Dipole-dipole energy in eV=",F20.10)',
!     &     dipen*pabs**2*charge**2 * Ryd / alat(1)**3 
!
!      case('tet','TET','tetra','TETRA')
!        c_over_a=alat(3)/alat(1)
!        vol=c_over_a
!        imax=500
!        dipen=0.0D0
!        do i=-imax,imax
!        do j=-imax,imax
!        do k=-imax,imax
!              if (abs(i)+abs(j)+abs(k).gt.0) then
!                    rvec(1)=dble(i)
!                    rvec(2)=dble(j)
!                    rvec(3)=c_over_a*dble(k)
!                    rabs=sqrt(dot_product(rvec,rvec))
!                    p_r=dot_product(pvec,rvec)/(pabs*rabs)
!                    p_epsil1_r=dot_product(p_epsil1,rvec)/(pabs * rabs)
!                    dipen=dipen+(p_epsil1_p-3.0d0*p_epsil1_r*p_r)
!     &                    /rabs**3
!              end if
!        end do
!        end do
!        end do
!        print*,"Dipole-dipole sum=",dipen
!        print '("Dipole-dipole energy in Ryd=",F20.10)',
!     &     dipen*pabs**2*charge**2 / alat(1)**3 
!        print '("Dipole-dipole energy in eV=",F20.10)',
!     &     dipen*pabs**2*charge**2 * Ryd / alat(1)**3 
!      case default
!        goto 100
!      end select
!      if(talk) then 
!        print fsubendext,'dipoleen'
!      end if
!      return
!! errors
!100   nerr=nerr+1
!      print ferrmssg, "Sorry, lattice type not implemented"
!      return
!
!      end subroutine

c---------------------------------------------------------------------
        
        end module 
