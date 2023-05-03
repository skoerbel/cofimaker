      module defs
      ! definitions, constants and type declarations
     
      implicit none
      
      ! variable types  
      type atom
        sequence
        character name*8     ! element name, a bit longer to allow for e.g. additional numbering, though most routines use only the first 2 characters
        character core*4      ! used for GULP (shell model)
        double precision :: where(3)  ! fractional coordinates
        double precision :: abswhere(3) ! absolute coordinates in Angstrom
        double precision charge  ! used for GULP: ionic charges in classical potentials
        double precision mass  ! used in cfg files
        double precision distance  
        double precision :: occ  ! occupation (in GULP fractional occupations are possible) 
        double precision, allocatable :: properties(:) ! in cfg files atomic properties are included (forces,...)
        logical written
        logical reducible  ! equivalent to another atom by symmetry
        integer coord ! number of nearest neighbors
        double precision, allocatable :: neighborcoords(:,:)
        double precision, allocatable :: bondlengths(:)
        double precision, allocatable :: bondangles(:)
        character*2, allocatable :: neighbors(:)     ! element name, a bit longer to allow for e.g. additional numbering, though most routines use only the first 2 characters
        integer, allocatable :: ineighbors(:)     ! atom numbers of neighbors
        double precision, allocatable :: nnsites(:,:)     ! element name, a bit longer to allow for e.g. additional numbering, though most routines use only the first 2 characters
        character label*256
        double precision BEC(1:3,1:3)
        double precision spin
        double precision force(1:3)
        character*2 flags(1:3) 
      end type atom
      type element
            sequence
            character name*8
            integer howmany  ! number of atoms of this species in the structure file
            double precision mass  ! atomic mass, usuaaly not used, is contained in cfg files
            double precision charge ! ionic charge in classical potentials
            type(atom), allocatable :: atoms(:)
      end type
      type symop ! symmetry operations
        sequence
        double precision mat(1:3,1:3)  ! e.g. rotation matrix
        double precision trala(1:3)  ! translation vector
      end type
      type kpoint        
        sequence
        double precision kpt(1:3)
        double precision kptfrac(1:3)
        double precision weight
        double precision phase
        integer istar
        integer nstar
        logical irred
        integer isymm ! kred=symops(isymm)*kirred
        integer n_neighbors
        integer, allocatable :: i_neighbors(:)
        double precision, allocatable :: e_up(:),e_down(:)
        double precision, allocatable :: neighbor_kptsfrac(:,:) 
        double precision, allocatable :: neighbor_kpts(:,:) 
        double precision :: nn_dist_min(1:3,1:3)  ! smallest distance to neighbor kpoints for each direction
        double precision, allocatable :: dE_dk(:,:,:) ! energy gradient for band 1:n in kspace (3 coordinates) for 2 spin channels
      end type  
      type pspot
        sequence
        character name*8
        double precision, allocatable :: r(:) ! radial grid
        double precision, allocatable :: wf(:,:) ! pseudo wavefunction (l,r)
        double precision, allocatable :: vps(:,:) ! pseudo potential (l,r)
        double precision, allocatable :: pvdens(:) ! atomic pseudo valence electron density (r)
        double precision, allocatable :: pcdens(:) ! partial core electron density (r)
        double precision, allocatable :: occ(:) ! atomic pseudo orbital occupation
        double precision zion ! ionic charge (=number of valence electrons)
        integer lmaxps ! lmax for this species, should be >= lmax (crystal)
        integer lloc ! l of local pseudopotential
        integer ixc ! xc functional index (like in abinit)
        double precision znuc ! nuclear charge, identifies species
      end type pspot

      ! constants
      character elements(100)*2
      data elements / 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ',
     &                'F ', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ',
     &                'Cl', 'Ar', 'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr',
     &                'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge',
     &                'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr',
     &                'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
     &                'In', 'Sn', 'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba',
     &                'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
     &                'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',
     &                'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
     &                'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra',
     &                'Ac', 'Th', 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm',
     &                'Bk', 'Cf', 'Es', 'Fm' /
      double precision masses(100)
      data masses/1.01d0,4.00d0,6.94d0,9.01d0,10.81d0,12.01d0,14.01d0,    &
     &          16.00d0,                                                  &
     &          19.00d0,20.18d0,22.99d0,24.31d0,26.98d0,28.09d0,30.97d0,  &
     &          32.07d0,35.45d0,39.95d0,39.10d0,40.08d0,44.96d0,47.88d0,  &
     &          50.94d0,52.00d0,54.94d0,55.85d0,58.93d0,58.69d0,63.55d0,  &
     &          65.39d0,69.72d0,72.61d0,74.92d0,78.96d0,79.90d0,83.80d0,  &
     &          85.47d0,87.62d0,88.91d0,91.22d0,92.91d0,95.94d0,98.00d0,  &
     &          101.07d0,102.91d0,106.42d0,107.87d0,112.41d0,114.82d0,    &
     &          118.71d0,121.75d0,127.60d0,126.91d0,131.29d0,132.91d0,    &
     &          137.33d0,138.91d0,140.12d0,140.91d0,144.24d0,145.00d0,    &
     &          150.36d0,151.97d0,157.25d0,158.93d0,162.50d0,164.93d0,    &
     &          167.26d0,168.93d0,173.04d0,174.97d0,178.49d0,180.95d0,    &
     &          183.85d0,186.21d0,190.20d0,192.22d0,195.08d0,196.97d0,    &
     &          200.59d0,204.38d0,207.20d0,208.98d0,209.00d0,210.00d0,    &
     &          222.00d0,223.00d0,226.03d0,227.03d0,232.04d0,231.04d0,    &
     &          238.03d0,237.05d0,244.00d0,243.00d0,247.00d0,247.00d0,    &
     &          251.00d0,252.00d0,257.00d0/
      INTEGER, DIMENSION (1:19) :: n_list,l_list
      DOUBLE PRECISION, DIMENSION (1:19) :: occ_list
      data n_list /1,      2,     2,     3,     3,     4,     3, 
     &             4,      5,     4,     5,     6,     4,     5, 
     &             6,      7,     5,     6,     7/
      data l_list /0,      0,     1,     0,     1,     0,     2,
     &             1,      0,     2,     1,     0,     3,     2,  
     &             1,      0,     3,     2,     1/
      data occ_list /2.0D0,  2.0D0, 6.0D0, 2.0D0, 6.0D0, 2.0D0,10.0D0,
     &             6.0D0,  2.0D0,10.0D0, 6.0D0, 2.0D0,14.0D0,10.0D0,
     &             6.0D0,  2.0D0,14.0D0,10.0D0, 6.0D0/
      !
      ! definitions
      integer nerr,nwarn,ncomm
      logical talk ! if (.not.talk), the program talks less 
      double precision tolfrac(1:3)
      parameter(tolfrac=0.05D0)
      integer, parameter :: lmax=2     ! lmax for crystal states
      double precision, parameter :: hdist_def=1.250D0,nndist_def=3.0D0 ! defaults for hydrogen bond length and upper limit for NN bond length 
      double precision :: vecs_def(1:3,1:3) ! default (dummy) lattice vectors, if none are contained in the file 
      data vecs_def / 1.0D0, 0.0D0, 0.0D0,
     &                0.0D0, 1.0D0, 0.0D0, 
     &                0.0D0, 0.0D0, 1.0D0  /
      !
      ! constants
      complex icomp
      parameter(icomp=(0.0d0,1.0d0))
      double precision pi,bohr,polfac_Angs_e,ec,Ryd,hartree,kcalpermole
      parameter(pi=acos(-1.0D0),bohr=0.52917721092D0,
     & Ryd=13.60569253D0,hartree=2.0D0*Ryd,kcalpermole=0.043d0)
      parameter(ec=1.60217648740E-19)
      double precision, parameter :: em=9.1093829140E-31
      double precision, parameter :: hbar=1.054571726E-34     
      parameter(polfac_Angs_e=1.60217656535E3)
      double precision, parameter :: Angstrom=1.0E-10,  
     &     epsilon0=8.8541878170E-12   ! vacuum permittivity
      double precision, parameter:: c_light=299792458.0d0

      ! special types and defs. for subroutine crysanpero
      double precision cornerweight,centerweight,faceweight
      double precision cornerdist,centerdist,facedist
      double precision tol_cell
      parameter(cornerweight=0.125D0,centerweight=1.0D0,
     &            faceweight=0.5D0)
      parameter(cornerdist=0.7D0,centerdist=0.6D0,
     &            facedist=0.7D0)
      parameter(tol_cell=0.4D0)
      type patom
        sequence
        integer:: number
        character:: name*2 
        integer:: direc 
        double precision:: where(3)
        double precision:: charge(1:3,1:3)
        double precision:: weight
        double precision:: distance
        character:: geometry*6
      end type patom
      type pelement
        sequence
        integer number
        integer:: howmany
        character:: name*2 
        double precision:: charge(1:3,1:3)
        double precision:: weight
        double precision:: distance
      end type pelement
      
      ! formatting 
      character fsubend*100,fsubendext*100,fsubstart*100,ferrmssg*100,
     &   fwarn*100,fcomm*100,fnumexp*256
      parameter(fsubstart="(8x,'This is subroutine ',A,'.')")
      parameter(fsubend="(8x,'Subroutine finished.')")
      parameter(fsubendext="(8x,'Subroutine ',A,
     &              ' finished.')")
      parameter(ferrmssg="(//8x,'Error: ',A,'.',//8x)")
      parameter(fwarn="(//8x,'Warning: ',A,'.',//8x)")
      parameter(fcomm="(//8x,'Comment: ',A,'.',//8x)")
      parameter(fnumexp="(/8x,'Your result:',1x,1p,e17.10)")
      character formatcfgcoords*21
      data formatcfgcoords /"(3(F15.8),   (F10.4))"/
      character FMT1*1024,FMT2*1024,FMT3*1024,FMT4*1024,FMT5*1024
      character formatvar*100
      character ftime*256,fcurrtime*256
      parameter(fcurrtime="(//2x,'Current time:',F20.6,' s')")
      parameter(ftime="(//2x,'Time needed for ',A50,':',F20.6,' s')")
      character*1024, parameter :: fsection="(/2(100('#')/),
     & '###',1x,a92,1x,'###'/2(100('#')/))"
      !
      ! misc.:
      double precision, parameter :: infty=1.0D12


c---------------------------------------------------------------------

      contains
 
!---------------------------------------------------------------------

      subroutine error(errmssg)
      implicit none
      character(len=*), intent(in) :: errmssg
      nerr=nerr+1   
      print ferrmssg,errmssg
              
      end subroutine error

!---------------------------------------------------------------------

      subroutine error_stop(errmssg)
      implicit none
      character(len=*), intent(in) :: errmssg
      nerr=nerr+1   
      print ferrmssg,errmssg
      stop  
              
      end subroutine error_stop

!---------------------------------------------------------------------

      subroutine warning(mssg)
      implicit none
      character(len=*), intent(in) :: mssg
      nwarn=nwarn+1   
      print fwarn,mssg
              
      end subroutine warning

!---------------------------------------------------------------------

      end module
     
