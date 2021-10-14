        module energy
        ! contains routines that calculate energy contributions


        contains

c---------------------------------------------------------------------

        double precision function e_hartree(rhog,glatt,vol)
        ! calculates the Hartree energy (in Hartree)
        use defs
        use misc
        implicit none
        complex*16 :: rhog(:)  ! electronic density in Fourier-space
        type(kpoint) :: glatt(:) ! reciprocal lattice vectors
        double precision vol, ehartree
        ! internal variables:
        integer ig,ng
        double precision gabs ! modulus of G vector
        
        ng=size(glatt,1)
        e_hartree=0.0D0
        do ig=1,ng
          gabs=absvec(glatt(ig)%kpt(1:3))
          if (gabs.gt.0.0D0) e_hartree=e_hartree+abs(rhog(ig))**2
     &    /gabs**2
        end do
        e_hartree=e_hartree*2.0D0*Pi/vol

        end function

c---------------------------------------------------------------------

        double precision function e_electron_ion(rhog,vpslocg,glatt,
     &     vol)
        ! calculates the local electron-ion energy
        use defs
        use misc
        implicit none
        complex*16 :: rhog(:)  ! electronic density in Fourier-space
        complex*16 :: vpslocg(:)  ! local pseudopotential in Fourier-space
        type(kpoint) :: glatt(:) ! reciprocal lattice vectors
        double precision vol ! volume
        ! internal variables:
        integer ig,ng
        double precision gabs ! modulus of G vector
        
        ng=size(glatt,1)
        e_electron_ion=0.0D0
        do ig=1,ng
          gabs=absvec(glatt(ig)%kpt(1:3))
          if (gabs.gt.0.0D0) e_electron_ion=e_electron_ion                &
     &        +conjg(rhog(ig))*vpslocg(ig)
        end do
        e_electron_ion=e_electron_ion/vol

        end function

c---------------------------------------------------------------------        

        end module
