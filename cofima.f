      program cofima 
      ! this program manipulates coordinate files. Usage: cofima [options] [file(s)]
      ! for a list of options, type "cofima --help".
      
      use defs
      use for_lammps
      use for_mbpp
      use for_gulp
      use for_dlpoly
      use misc
      use postproc
      use crysanpero
      use transform
      use mult
      use convert
      use cutmod
      use distri
      use aver
      use md
      use replace
      use cluster
      use readcoords
      use writecoords
      use gaussians
      use madelung
      use deriv
      use integr
      use ipol
      use help
      use symmetry, only : latsys
      use linalg, only : quickhull, quickhull2, get_dfh, fdet3,           &
     &                   get_mateigs, rotate_matrix, polyhedral_volume,   &
     &                   matmul_fortran
     
 
      implicit none

      ! variables for the command line options 
      character(len=*), parameter :: version = '3.4.0'
      character(len=32) :: arg
      character(len=8) :: date
      character(len=10) :: time
      character(len=40) :: inputfile
      character(len=40) :: outputfile
      character(len=40) :: fincoords
      character(len=40) :: foutcoords
      !character(len=40) :: finform          ! format of input coordinate file
      character(len=40) :: foutform         ! format of output coordinate file
      integer :: i   ! command line argument index
      
      ! used by "convert":
      character(len=256) :: convfile1          ! input filename
      character(len=256) :: convfile2         ! output filename
      character(len=256) :: convformat1       ! format of input file
      character(len=256) :: convformat2          ! format of outputfile

      ! used by "transform":
      integer :: nsc(1:3)
      character(len=10) :: trans,transoption
      character(len=20) :: matchar,vecchar,intchar
      double precision matrix(1:3,1:3),shift(1:3),tol(1:3),tol1,tol2(2)
      double precision inmat(1:3,1:3),outmat(1:3,1:3),outimat(1:3,1:3)
      double precision invecs(1:3,1:3),outvecs(1:3,1:3)


      ! used by ansg2borh, bohr2angs:
      double precision  numinbohr,numinangs,numin,numout

      ! used by "av":
      integer firststep, step, laststep

      ! used by coord:
      double precision coordmean
      integer coordmin,coordmax

      ! used by "cut":
      integer axisint
      double precision cuthere,hdist
      character*20 cutinfile,cutoutfile,cutinform,cutoutform,abovechar,
     &             axischar,cutchar

      ! used by nicecluster:
      double precision coordminpar,coordmeanpar,dipolepar  ! thresholds for acceptance of a cluster (see below)
  
      ! used by "dist" and "rdf":
      !integer axisint
      double precision binwidth,broad,rmax,rcut
      character*20 distinfile,distinform,
     &             distchar

      ! used by dipole
      double precision dipole(1:3) 

      !used by divide:
      double precision number1,number2
      
      ! used by ecoulomb
      double precision ecoul 

      ! used by help
      character*20 morehelpon

      ! used by lmpthermo
      character lmplog*40

      ! used by mbppnebana and nebchain:
      integer nimg

      ! used by "mdana":
      !integer :: j,k,nsc(1:3)
      character(len=10) :: bndrychar
      character(len=40) :: trjfile,trjatomfile,trjformat
      double precision trjbndry(1:3)
      
      ! used by traj:
      integer thenumber
      character thespecies*8,nnana*20,nnatom*20
      double precision nndist

      ! used by replace
      character theotherspecies*8 

      ! used by rmeqv:
      character symmformat*40

      ! other variables
      character line*100
      double precision t_start,t_end,add,angle,radius,radius0,afloat
      integer intnum,intnum2,intnum3,intnum4,intnum5
      integer :: j,k,iatom,step2
      integer, allocatable :: intvec(:),intvec2(:),intvec3(:),intvec4(:)
      logical cmdlonly,test ! command-line-only mode, test-mode
      
      character infile*1024,infile2*1024,informat*1024,outfile*1024,
     &           outformat*1024,charac*1024,charac2*1024,charac3*1024
      character format1*1024,format2*40,format3*40,
     &          filename1*1024,filename2*80,filename3*80 
      character thetype*20, lattice*7
      !character FMT1*1024,FMT2*1024
      logical isopen12  ! checks if an output file is open
      type(atom), allocatable :: atoms(:),atoms2(:),atoms3(:),atoms4(:)
      type(element), allocatable :: species(:),species2(:),species3(:),
     &     species4(:) 
      integer natoms,natoms2,natoms3,nspecies,nspecies2
      double precision vecs(1:3,1:3),vecs2(1:3,1:3),origin(1:3),          &
     &     vector(1:3),epsil,energy
      double precision charge
      double precision, allocatable :: allovec1(:), allovec2(:)
      double precision, allocatable :: allomat1(:,:),allomat2(:,:)
      double precision, allocatable :: allomat3(:,:)
      logical lnum,lnum1

c---------------------------------------------------------------------

      ! Initializing...
      call cpu_time(t_start)
      nerr=0                            ! error counter
      nwarn=0                            ! error counter
      ncomm=0                            ! error counter
      inputfile='INPUT.COFIMA'          ! default input file
      outputfile='OUTPUT.COFIMA'        ! DEFAULT OUTPUT FILE
      cmdlonly=.false.                  ! default: command-line only mode off.
      talk=.true.                    ! default: the program talks a lot.  
      
!      print*,'***************************************************'
!      print '(" * This is cofima, v ",A,
!     & ".                        *")', version
!      print*,'*                                                 *'
!      print*,'* New capabilities/changes:                       *'
!      print*,'* -compatible with gfortran since v 3.4.0         *'
!      print*,'* -new option "sphere" for the cut command        *'
!      print*,'* -new format for coordinate files: "abinit"      *'
!      print*,'*  (abinit output file)                           *'
!      print*,'***************************************************'

      print*,'***************************************************'
!      print '(" * This is cofima, v ",A,
!     & ".                        *")', version
!      print*,'*                                                 *'
!      print*,'* New capabilities/changes: see file CHANGES      *'
      print '(" * This is cofima, the beta version                *")'
      print*,'***************************************************'
      
      ! read command line options
      i=1
      do while (i.le.command_argument_count())
         call get_command_argument(i, arg)
    
         select case (arg)
         case ('--add')
                 ! adds a number to another number
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) number1
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) number2
!                 print '(8x,"Your result: ",F)',number1+number2
                 !print '(8x,"Your result: ",1p,e12.3)',number1+number2
                 !print '(8x,"Your result: ",1p,e12.6)',number1+number2
                 print fnumexp,number1+number2
         case("--addatoms")
               ! add atoms of one species
               i=i+1
               call get_command_argument(i,infile)
               i=i+1
               call get_command_argument(i,informat)
               i=i+1
               call get_command_argument(i,thespecies)
               i=i+1
               call get_command_argument(i,charac)
               read(charac,*) thenumber
               i=i+1
               call get_command_argument(i,outfile)
               allocate(allomat1(1:thenumber,1:3))
               do j=1,thenumber
                 do k=1,3
                   i=i+1
                   call get_command_argument(i,charac)
                   read(charac,*) allomat1(j,k)
                 end do ! k=cartesian coordinate 
               end do ! j=atom number
               call add_atoms(infile,informat,thespecies,thenumber,
     &              allomat1,outfile)
         case("--addeqv")
               ! add atoms to a set of irreducible atoms by applying
               ! symmetries read from a file 
                 i=i+1
                 call get_command_argument(i,infile)
                 i=i+1
                 call get_command_argument(i,informat)
                 i=i+1
                 call get_command_argument(i,symmformat)
                 i=i+1
                 call get_command_argument(i,outfile)
                 tol1=1.0D-5
                 do while (i+1.le.command_argument_count()) 
                   i=i+1
                   call get_command_argument(i,charac)
                   read(charac,*) tol1
                 end do
                 call addeqv(infile,informat,symmformat,outfile,tol1)
         case ("--angs2bohr")
                 ! converts between Angstrom and Bohr
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) numinangs
                 numinbohr=numinangs/bohr
!                 print '(8x,"Your number in Bohr: ",F)',numinbohr
                 print '(8x,"Your number in Bohr: ",1p,e12.6)',numinbohr
         case("--av")
                ! average lattice parameters
                i=i+1
                call get_command_argument(i,infile)
                i=i+1
                call get_command_argument(i,informat)
                i=i+1
                call get_command_argument(i,charac)
                read(charac(1:20),*) firststep
                i=i+1
                call get_command_argument(i,charac)
                read(charac(1:20),*) step
                i=i+1
                call get_command_argument(i,charac)
                read(charac(1:20),*) laststep
                call av(infile,informat,
     &               firststep,step,laststep)
         case ("--bandgap")
                 ! reads the bandgap from a DFT output file
                 i=i+1
                 call get_command_argument(i,infile)
                 i=i+1
                 call get_command_argument(i,informat)
                 call bandgap(infile,informat)
         case ("--basis4fiesta")
                 ! creates an auxiliary basis (el2.ion) to be read by
                 ! fiesta, where el is an element. Input: element name,
                 ! zetamin, zetamax, nzeta, lmax. E.g.: 
                 ! cofima --basis4fiesta He 0.8 4.0 8 2 
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) thespecies! element name
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) number1   ! zetamin
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) number2   ! zetamax
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) intnum    ! nzeta
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) intnum2   ! lmax
                 rmax=20.0d0 ! default
                 intnum3=10000 ! npoints, default
                 do while (i+1.le.command_argument_count()) 
                   i=i+1
                   call get_command_argument(i,charac)
                   select case(charac)
                   case('rmax') 
                     i=i+1
                     call get_command_argument(i,charac2)
                     read(charac2,*) rmax   ! rmax
                   case('npoints')
                     i=i+1
                     call get_command_argument(i,charac2)
                     read(charac2,*) intnum3   ! npoints
                   case default
                   end select
                 end do
                 call basis4fiesta(thespecies,number1,number2,intnum,
     &             intnum2,rmax,intnum3)
         case ("--bohr2angs")
                 ! converts between Bohr and Angstrom
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) numinbohr
                 numinangs=numinbohr*bohr
!                 print '(8x,"Your number in Angstrom: ",F)',numinangs
                 print '(8x,"Your number in Angsstrom: ",1p,e12.3)',
     &                numinangs
         case ("--broaden")
                 ! applies a gaussian broadening to a distribution read
                 ! from a file
                 ! 
                 ! filename:
                 i=i+1
                 call get_command_argument(i,infile)
                 !
                 ! column with x values:
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) intnum 
                 !                  
                 ! column with function values:
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) intnum2 
                 !
                 ! periodic?:
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) lnum 
                 !
                 ! periodicity length (any number will do if not
                 ! periodic):
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) afloat 
                 !
                 ! broadening to be applied:
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) broad
                 ! 
                 call broaden_dist(infile,intnum,intnum2,lnum,            &
     &                             afloat,broad)
                 !
         case ("--cageZnO")
                 ! creates a cage-like cluster of n ZnO units. n is read from the command line
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) intnum
                 call cageZnO(intnum)
         case ('--cmdlonly')
               cmdlonly=.true.
         case ('--com','--COM')
               ! calculate the Center Of the Molecule 
            i=i+1
            call get_command_argument(i,infile)
            i=i+1
            call get_command_argument(i,informat)
            call get_com(infile,informat)
         case ('--comass','--COMASS')
               ! calculate the Center Of Mass 
            i=i+1
            call get_command_argument(i,infile)
            i=i+1
            call get_command_argument(i,informat)
            call get_comass(infile,informat)
         case ('--compare','--overlaps')
               ! compare two 1d fields by calculating their overlap
               ! Ov=<f1|f2> / (|f1| * |f2| )
            i=i+1
            call get_command_argument(i,filename1)
            i=i+1
            call get_command_argument(i,filename2)
            call compare(filename1,filename2)
         case ('--convert')
            i=i+1
            call get_command_argument(i,convfile1)
            i=i+1
            call get_command_argument(i, convfile2)
            i=i+1
            call get_command_argument(i,convformat1)
            i=i+1
            call get_command_argument(i, convformat2)
            if(talk) print '(8x,"File ",A," will be converted to file ",  &
     &          A,".")',trim(convfile1),trim(convfile2)
            ! fake lattice vectors in case they are not included in the
            ! input coordinate file
            vecs=0.0D0
            vecs(1,1)=1.0D0
            vecs(2,2)=1.0D0
            vecs(3,3)=1.0D0
            ! read lattice vectors, if needed
            select case(convformat1)
            case("coorat","COORAT",'nwchem','NWCHEM','xyz','XYZ')
              if (i+9.le.command_argument_count()) then
                do j=1,3
                  i=i+1
                  call get_command_argument(i,vecchar)
                  read(vecchar(1:20),*)vecs(j,1)
                  i=i+1
                  call get_command_argument(i,vecchar)
                  read(vecchar(1:20),*)vecs(j,2)
                  i=i+1
                  call get_command_argument(i,vecchar)
                  read(vecchar(1:20),*)vecs(j,3)
                end do
              else  
                ncomm=ncomm+1  
                print fcomm, "Missing lattice vectors. The final coordin  &
     &ates may be wrong"
              end if  
            case default
            end select
            call conv(convfile1,convfile2,convformat1,
     &           convformat2,vecs)
         case ('--coord')
            ! determine the coordination numbers of the atoms in a file
            ! and print min, max, av 
            i=i+1
            call get_command_argument(i,infile)
            i=i+1
            call get_command_argument(i, informat)
            vecs=vecs_def
            ! read lattice vectors, if needed
            select case(informat)
            case("coorat","COORAT",'nwchem','NWCHEM','xyz','XYZ')
              if (i+9.le.command_argument_count()) then
                do j=1,3
                  i=i+1
                  call get_command_argument(i,vecchar)
                  read(vecchar(1:20),*)vecs(j,1)
                  i=i+1
                  call get_command_argument(i,vecchar)
                  read(vecchar(1:20),*)vecs(j,2)
                  i=i+1
                  call get_command_argument(i,vecchar)
                  read(vecchar(1:20),*)vecs(j,3)
                end do
              else  
                call error_stop('Please provide lattice vectors after fo  &
     &rmat, or use a format that contains lattice vectors')
              end if  
            case default
            end select
            nndist=nndist_def
            if (i+1.lt.command_argument_count()) then
              i=i+1
              call get_command_argument(i,charac)
              select case(charac)
              case('-cutoff')
                i=i+1
                call get_command_argument(i,charac)
                read(charac,*) nndist
              case default
                call error_stop('unknown option')
              end select
              !
            end if ! i+1.lt.command_argument_count()
            call read_coords(infile,informat,atoms,natoms,species,
     &                       nspecies,vecs)
            call coord(atoms,vecs,nndist,coordmin,coordmax,
     &        coordmean)
            print '(8x,"min, max, av coord. num.: ",2(I6),F10.6)',
     &              coordmin,coordmax,coordmean 
            open(51,file="COORDINATION",status='replace')
            do iatom=1,size(atoms)
              write(51,*) atoms(iatom)%name,atoms(iatom)%coord
            end do ! iatom
            close(51)
            open(51,file="COORDINATION_DETAILED",status='replace')
            do iatom=1,size(atoms)
             if (1.lt.atoms(iatom)%coord) then
              write(51,*) atoms(iatom)%name,atoms(iatom)%coord,           &
     &         (atoms(iatom)%neighbors(j),' ',j=1,atoms(iatom)%coord),    &
     &         (atoms(iatom)%bondlengths(j),' ',j=1,atoms(iatom)%coord),  &
     &          size(atoms(iatom)%bondangles),                            &
     &         (atoms(iatom)%bondangles(j),' ',j=1,                       &
     &          size(atoms(iatom)%bondangles))
             else
              if (1.eq.atoms(iatom)%coord) then 
               write(51,*) atoms(iatom)%name,atoms(iatom)%coord,          &
     &         (atoms(iatom)%neighbors(j),' ',j=1,atoms(iatom)%coord),    &
     &         (atoms(iatom)%bondlengths(j),' ',j=1,atoms(iatom)%coord)
              else
               write(51,*) atoms(iatom)%name,atoms(iatom)%coord         
              end if
             end if
            end do ! iatom
            close(51)
            ! begin print coordinates of neighbors 
            open(51,file="COORDINATION_NEIGHBOR_COORDS",                  &
     &           status='replace')
            do iatom=1,natoms
             if (1.lt.atoms(iatom)%coord) then
               write(51,*) atoms(iatom)%name,atoms(iatom)%where,          &
     &         atoms(iatom)%abswhere,atoms(iatom)%coord,                  &
     &         (atoms(iatom)%neighbors(j),' ',j=1,atoms(iatom)%coord),    &
     &         (atoms(iatom)%neighborcoords(j,1:3),' ',                   &
     &        j=1,atoms(iatom)%coord)        
             else
              if (1.eq.atoms(iatom)%coord) then 
                write(51,*) atoms(iatom)%name,atoms(iatom)%where,         &
     &         atoms(iatom)%abswhere,atoms(iatom)%coord,                  &
     &         (atoms(iatom)%neighbors(j),' ',j=1,atoms(iatom)%coord),    &
     &         (atoms(iatom)%neighborcoords(j,1:3),' ',                   &
     &         j=1,atoms(iatom)%coord)
              else
               write(51,*) atoms(iatom)%name,atoms(iatom)%coord         
              end if
             end if
            end do ! iatom
            close(51)
            ! end print coordinates of neighbors 
         case('--corr0')
                 i=i+1
                 call get_command_argument(i,infile)
                 i=i+1
                 call get_command_argument(i,thespecies)
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) thenumber
                 i=i+1
                 call get_command_argument(i,thetype)
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) rmax
                 binwidth=0.2D0
                 if (i+1.le.command_argument_count()) then
                   i=i+1
                   call get_command_argument(i,charac)
                   read(charac,*) binwidth
                 end if
                 call propcorr0(infile,thespecies,thenumber,thetype,
     &                rmax,binwidth)
         case('--corr')
                 i=i+1
                 call get_command_argument(i,infile)
                 i=i+1
                 call get_command_argument(i,thespecies)
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) thenumber
                 i=i+1
                 call get_command_argument(i,thetype)
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) rcut
                 ! read an additive constant (optional) to be added to
                 ! the property that interests you
                 add=0.0D0
                 if (i+1.le.command_argument_count()) then
                       i=i+1
                       call get_command_argument(i,charac)
                       if (charac.eq."add") then
                             i=i+1
                             call get_command_argument(i,charac)
                             read(charac,*) add
                       else
                             i=i-1
                       end if
                 end if
                 ! read binwidth, maximum distance for which the
                 ! correlation function is calculated, and the gaussian
                 ! broadening (optional)
                 binwidth=0.2D0
                 rmax=20.0D0
                 broad=0.2D0 
                 if (i+1.le.command_argument_count()) then
                   i=i+1
                   call get_command_argument(i,charac)
                   read(charac,*) binwidth
                 end if
                 if (i+1.le.command_argument_count()) then
                   i=i+1
                   call get_command_argument(i,charac)
                   read(charac,*) rmax
                 end if
                 if (i+1.le.command_argument_count()) then
                   i=i+1
                   call get_command_argument(i,charac)
                   read(charac,*) broad
                 end if
                 call propcorr(infile,thespecies,thenumber,thetype,
     &                rcut,add,binwidth,rmax,broad)
         case('--create_neb_chain')        
                 ! read two coordinate files with starting and end
                 ! coordinates and create linearly interpolated chain of
                 ! coordinates
                 i=i+1
                 call get_command_argument(i,filename1)
                 i=i+1
                 call get_command_argument(i,filename2)
                 i=i+1
                 call get_command_argument(i,format1)
                 nimg=10 ! default
                 if (i+1.le.command_argument_count()) then
                     i=i+1
                     call get_command_argument(i,charac)
                     read(charac,*) nimg
                 end if
                 if(talk) print '(8x,"NEB chain with ",I0, "images will   &
     &be created between files ",A," and ",A)',nimg,trim(filename1),      &
     &           trim(filename2)
                 ! fake lattice vectors in case they are not included in the
                 ! input coordinate file
                 vecs=0.0D0
                 vecs(1,1)=1.0D0
                 vecs(2,2)=1.0D0
                 vecs(3,3)=1.0D0
                 vecs2=0.0D0
                 vecs2(1,1)=1.0D0
                 vecs2(2,2)=1.0D0
                 vecs2(3,3)=1.0D0
                 ! read lattice vectors, if needed
                 select case(format1)
                 case("coorat","COORAT",'nwchem','NWCHEM','xyz','XYZ')
                   if (i+18.le.command_argument_count()) then
                     do j=1,3
                       i=i+1
                       call get_command_argument(i,vecchar)
                       read(vecchar(1:20),*)vecs(j,1)
                       i=i+1
                       call get_command_argument(i,vecchar)
                       read(vecchar(1:20),*)vecs(j,2)
                       i=i+1
                       call get_command_argument(i,vecchar)
                       read(vecchar(1:20),*)vecs(j,3)
                     end do
                     do j=1,3
                       i=i+1
                       call get_command_argument(i,vecchar)
                       read(vecchar(1:20),*)vecs2(j,1)
                       i=i+1
                       call get_command_argument(i,vecchar)
                       read(vecchar(1:20),*)vecs2(j,2)
                       i=i+1
                       call get_command_argument(i,vecchar)
                       read(vecchar(1:20),*)vecs2(j,3)
                     end do
                   else  
                     ncomm=ncomm+1  
                     print fcomm, "Missing lattice vectors. The final co  &
     &ordinates may be wrong"
                   end if  
                 case default
                 end select
                 call create_neb_chain(filename1,filename2,format1,vecs,  &
     &                            vecs2,nimg)
         case ('--crysanpero')
            call scrysanpero()
         case("--cut")
               ! remove atoms above/below x/y/z=number
               ! or inside/outside a sphere
                 i=i+1
                 call get_command_argument(i,infile)
                 i=i+1
                 call get_command_argument(i,outfile)
                 i=i+1
                 call get_command_argument(i,cutinform)
                 i=i+1
                 call get_command_argument(i,cutoutform)
                 i=i+1
                 call get_command_argument(i,abovechar)
                 i=i+1
                 call get_command_argument(i,axischar)
                 i=i+1
                 read(axischar(1:20),*) axisint
                 call get_command_argument(i,cutchar)
                 read(cutchar(1:20),*) cuthere
                 if(abovechar.eq."sphere".or.abovechar.eq.'SPHERE'.or.    &
     &              abovechar.eq."plane".or.abovechar.eq."PLANE") then
                   i=i+1
                   call get_command_argument(i,cutchar)
                   read(cutchar(1:20),*) shift(1)
                   i=i+1
                   call get_command_argument(i,cutchar)
                   read(cutchar(1:20),*) shift(2)
                   i=i+1
                   call get_command_argument(i,cutchar)
                   read(cutchar(1:20),*) shift(3)
                 end if
                 ! nearest neighbor distance (used for hydrogenation)
                 !nndist=3.0d0 ! default value in Angstrom
                 nndist=nndist_def ! default value in Angstrom
                 if (i+1.le.command_argument_count()) then
                   i=i+1
                   call get_command_argument(i,charac)
                   read(charac,*) nndist
                 end if
                 ! distance along broken in which to place H
                 !hdist=1.25d0  ! default: 1.25 Angstrom
                 hdist=hdist_def  ! default: 1.25 Angstrom
                 if (i+1.le.command_argument_count()) then
                   i=i+1
                   call get_command_argument(i,charac)
                   read(charac,*) hdist
                 end if
                 if (i+1.le.command_argument_count()) then
                   i=i+1
                   call get_command_argument(i,charac)
                   select case(charac)
                   case('-hydr')
                     lnum=.true.
                   case default
                   end select
                 end if
                 if(talk) 
     &             print '(4x,"Using as H-bond length:",F10.6)', hdist
                 call cut(infile,outfile,cutinform,cutoutform,
     &                abovechar,axisint,cuthere,shift,nndist,hdist,lnum)
         case ('--dfh','--DFH')
            ! calculate distance from hull
            !
            ! read dimension:
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) k
            !
            ! read point:
            if (allocated(allovec1)) deallocate(allovec1)
            allocate(allovec1(k))
            do j=1,k
              i=i+1
              call get_command_argument(i,charac)
              read(charac,*) allovec1(j)
            end do  
            !
            ! read direction:
            if (allocated(allovec2)) deallocate(allovec2)
            allocate(allovec2(k))
            do j=1,k
              i=i+1
              call get_command_argument(i,charac)
              read(charac,*) allovec2(j)
            end do  
            !
            call get_dfh(k,allovec1,allovec2)
            !
            do while (i+1.le.command_argument_count())
              i=i+1
              ncomm=ncomm+1
              print fcomm,'addit. comm. line options found and ignored'
            end do
            !
         case ('--diffatoms')
            ! determine which atoms are common and which differ between
            ! two files
            i=i+1
            call get_command_argument(i,infile)
            i=i+1
            call get_command_argument(i,infile2)
            i=i+1
            call get_command_argument(i, informat)
            vecs=vecs_def
            ! read lattice vectors, if needed
            select case(informat)
            case("coorat","COORAT",'nwchem','NWCHEM','xyz','XYZ')
              if (i+9.le.command_argument_count()) then
                do j=1,3
                  i=i+1
                  call get_command_argument(i,vecchar)
                  read(vecchar(1:20),*)vecs(j,1)
                  i=i+1
                  call get_command_argument(i,vecchar)
                  read(vecchar(1:20),*)vecs(j,2)
                  i=i+1
                  call get_command_argument(i,vecchar)
                  read(vecchar(1:20),*)vecs(j,3)
                end do
              else  
                nwarn=nwarn+1  
                print fwarn, "Missing lattice vectors. The final 
     &                coordinates may be wrong."
              end if  
            case default
            end select
            ! tolerance for "equal" positions
            tol1=0.01d0 
            if (i+2.le.command_argument_count()) then
              i=i+1
              call get_command_argument(i,charac)
              select case (charac)
              case ("-tol")
                i=i+1
                call get_command_argument(i,charac)
                read(charac,*) tol1 
              case default
              end select
            end if
            call read_coords(infile,informat,atoms,natoms,species,
     &                       nspecies,vecs)
            call read_coords(infile2,informat,atoms2,natoms2,species2,
     &                       nspecies2,vecs)
            call diffatoms(atoms,atoms2,vecs,atoms3,atoms4,tol1)
            call getspecies(atoms3,species3)
            call getspecies(atoms4,species4)
            filename1=' '
            filename1(1:12)="COMMONATOMS."
            filename1(13:12+len_trim(informat))= trim(informat)
            call write_coords(filename1,informat,atoms3,size(atoms3),
     &        species3,size(species3),vecs)
            filename1=' '
            filename1(1:10)= "DIFFATOMS."
            filename1(11:10+len_trim(informat))= trim(informat)
            call write_coords(filename1,informat,atoms4,size(atoms4),
     &        species4,size(species4),vecs)
         case ('--deriv')
               ! calculates the first derivative of the second column wrt the
               ! first one from finite differences 
            i=i+1
            call get_command_argument(i,filename1)
            i=i+1
            call get_command_argument(i,charac) ! derivative mode (central differences,...)
            call derivs(filename1,charac)
         case ('--detmat3')
            do j=1,3
              do k=1,3
                i=i+1
                call get_command_argument(i,charac)
                read(charac,*) inmat(j,k)
              end do
            end do
            print '(8x,"matrix: ",/3(8x,3(E32.16)/))',                    &
     &           inmat(1,1:3),inmat(2,1:3),inmat(3,1:3)
            print '(8x,"determinant of matrix: ",E32.16)',fdet3(inmat)
         case("--dipole")
            ! calculate the dipole moment of a cluster of point charges.
            ! Caution: Works only for finite systems. 
            i=i+1
            call get_command_argument(i,infile)
            i=i+1
            call get_command_argument(i,informat)
            call read_coords(infile,informat,atoms,natoms,species,
     &                       nspecies,vecs)
            call dipmom(atoms,dipole) 
            print '(2x,"dipole moment in eV*Angs: ",3(F15.6),
     &              ", modulus:",F15.6)', dipole,absvec(dipole) 
!         case("--dipoleen0")
!               ! calculate the dipole-dipole energy for some point dipole
!               ! lattices, use scalar epsilon 
!                 i=i+1
!                 call get_command_argument(i,lattice)
!                 i=i+1
!                 call get_command_argument(i,charac)
!                 read(charac,*) charge
!                 do j=1,3
!                   i=i+1
!                   call get_command_argument(i,charac)
!                   read(charac,*) dipole(j)
!                 end do
!                 i=i+1
!                 call get_command_argument(i,charac)
!                 read(charac,*) epsil
!                 i=i+1
!                 call get_command_argument(i,charac)
!                 read(charac,*) vector(1)
!                 vector(2)=vector(1)
!                 vector(3)=vector(1)
!                 select case (lattice)
!                 case('tet','TET','tetra','TETRA')
!                   i=i+1
!                   call get_command_argument(i,charac)
!                   read(charac,*) vector(3)
!                 end select
!                 call dipoleen0(lattice,charge,dipole,epsil,vector,talk,
!     &                nerr)
!         case("--dipoleen")
!               ! calculate the dipole-dipole energy for some point dipole
!               ! lattices, use epsilon tensor 
!                 i=i+1
!                 call get_command_argument(i,lattice)
!                 i=i+1
!                 call get_command_argument(i,charac)
!                 read(charac,*) charge
!                 do j=1,3
!                   i=i+1
!                   call get_command_argument(i,charac)
!                   read(charac,*) dipole(j)
!                 end do
!                 ! epsilon tensor
!                 do j=1,3
!                   do k=1,3
!                     i=i+1
!                     call get_command_argument(i,charac)
!                     read(charac,*) matrix(j,k)
!                   end do
!                 end do
!                 ! lattice constant(s)
!                 i=i+1
!                 call get_command_argument(i,charac)
!                 read(charac,*) vector(1)
!                 vector(2)=vector(1)
!                 vector(3)=vector(1)
!                 select case (lattice)
!                 case('tet','TET','tetra','TETRA')
!                   i=i+1
!                   call get_command_argument(i,charac)
!                   read(charac,*) vector(3)
!                 end select
!                 call dipoleen(lattice,charge,dipole,matrix,vector,talk,
!     &                nerr)
         case("--dipen")
               ! calculate the dipole-dipole energy for a point dipole
               ! lattice, handles only scalar epsilon 
                 ! lattice vectors
                 do j=1,3
                   do k=1,3
                     i=i+1
                     call get_command_argument(i,charac)
                     read(charac,*) matrix(j,k)
                   end do
                 end do
                 ! dipole moment
                 do j=1,3
                   i=i+1
                   call get_command_argument(i,charac)
                   read(charac,*) dipole(j)
                 end do
                 ! epsilon 
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) epsil
                 call dipen(matrix,dipole,epsil)
         case("--dist")
               ! calculate distribution along x/y/z (optional: binwidth,
               ! broadening)
                 i=i+1
                 call get_command_argument(i,distinfile)
                 i=i+1
                 call get_command_argument(i,distinform)
                 i=i+1
                 call get_command_argument(i,axischar)
                 read(axischar(1:20),*) axisint
                 binwidth=0.01D0 ! fractional
                 broad=0.1D0 ! fractional
                 do while(i+1.lt.command_argument_count())
                   i=i+1
                   call get_command_argument(i,distchar)
                   select case(distchar)
                   case("-bin")
                     i=i+1
                     call get_command_argument(i,distchar)
                     read(distchar(1:20),*) binwidth
                   case("-broad")
                     i=i+1
                     call get_command_argument(i,distchar)
                     read(distchar(1:20),*) broad
                   case default
                    ncomm=ncomm+1  
                    if(talk) 
     &                print fcomm, "Default parameters for distribution 
     &will be used."
                  end select      
                end do  
                call dist(distinfile,distinform,
     &               axisint,binwidth,broad)
         case ('--divide')
                 ! divides a number by another number
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) number1
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) number2
!                 print '(8x,"Your result: ",F)',number1/number2
                 !print '(8x,"Your result: ",1p,e12.3)',
!     &                number1/number2
                  print fnumexp,number1/number2
         case("--ecoulomb")
            ! calculate the electrostatic energy of a cluster of point charges.
            ! Caution: Works only for finite systems. 
            i=i+1
            call get_command_argument(i,infile)
            i=i+1
            call get_command_argument(i,informat)
            call read_coords(infile,informat,atoms,natoms,species,
     &                       nspecies,vecs)
            call ecoulomb(atoms,ecoul) 
            print '(2x,"Electrostatic energy in eV: ",1(F15.6))',ecoul 
         case ('--ekin','--EKIN')
               ! calculates the kinetic energy from the density (works
               ! only for real WF) in cube files 
            i=i+1
            call get_command_argument(i,infile)
            i=i+1
            call get_command_argument(i,informat)
            ! read optional arguments
            step=1
            !rcut=0.650d0
            rcut=1.10d0
            do while (i+1.lt.command_argument_count()) 
              i=i+1
              call get_command_argument(i,charac)
              select case(charac)
                case("-step")
                 i=i+1
                 call get_command_argument(i,charac2)
                 read(charac2,*) step
                case default 
              end select
            end do
            call get_ekin(infile,informat,step)
         case ("--ev2hart","evolt2hartree","eV2Hart",
     &         "eVolt2Hartree")
                 ! converts between eV and Hartree
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) numin
                 numout=numin/hartree
                 !print '(8x,"Your number in Hartree: ",1p,e12.3)',numout
                 print fnumexp,numout
         case ("--ev2ryd","evolt2rydberg","eV2Ryd",
     &         "eVolt2Rydberg")
                 ! converts between eV and Rydberg
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) numin
                 numout=numin/Ryd
                 !print '(8x,"Your number in Rydberg: ",1p,e12.3)',numout
                 print fnumexp,numout
         case("--expand")
                 ! uniformly expands or shrinks a unit cell
                 i=i+1
                 call get_command_argument(i,infile)
                 i=i+1
                 call get_command_argument(i,informat)
                 i=i+1
                 call get_command_argument(i,outfile)
                 i=i+1
                 call get_command_argument(i,outformat)
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) number1 ! expansion factor
                 if(talk)  then
                   print '(2x,"Expanding (shrinking) cell by a factor of  &
     &                     ",F10.6)',number1
                 end if
                call read_coords(infile,informat,atoms,natoms,species,
     &                 nspecies,vecs)
                vecs=vecs*number1
                do i=1,natoms
                  call frac2abs(atoms(i)%where,vecs,atoms(i)%abswhere)
                end do
                call write_coords(outfile,outformat,atoms,natoms,         &
     &                 species,nspecies,vecs)
         case ("--expo_range")
                 ! determines the range of gaussian exponents zeta (range=zetamax-zetamin) 
                 ! from a file which should have the fiesta basis format, e.g. NWCHEM_Ga.ion , He.ion, Cu2.ion etc.  
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) filename1 ! filename
                 call expo_range(filename1)
         case ('--fourier_1d')
               ! calculates the Fourier Transform of the second column in a data file
               i=i+1
               call get_command_argument(i,filename1)
               i=i+1
               call get_command_argument(i,charac)
               lnum=.false.
               lnum1=.false.
               do while(i+1.lt.command_argument_count())
                 i=i+1
                 call get_command_argument(i,distchar)
                 select case(distchar)
                 case("-nq")
                   i=i+1
                   call get_command_argument(i,distchar)
                   read(distchar(1:20),*) intnum
                   lnum=.true.
                 case("-qmax")
                   i=i+1
                   call get_command_argument(i,distchar)
                   read(distchar(1:20),*) number1
                   lnum1=.true.
                 case default
                 ncomm=ncomm+1  
                 if(talk) print fcomm,"Default parameters will be used."
                end select      
               end do  
               if (.not.lnum.and..not.lnum1) then
                call fourier_1d(filename1,charac)
               end if
               if (lnum.and..not.lnum1) then
                call fourier_1d(filename1,charac,nq_chosen=intnum)
               end if
               if (.not.lnum.and.lnum1) then
                call fourier_1d(filename1,charac,qmax_chosen=number1)
               end if
               if (lnum.and.lnum1) then
                call fourier_1d(filename1,charac,qmax_chosen=number1,     &
     &                          nq_chosen=intnum)
               end if
         case ('--fourier_1d_reverse')
               ! calculates the reverse (q -> r) Fourier Transform of the second column in a data file
               i=i+1
               call get_command_argument(i,filename1)
               i=i+1
               call get_command_argument(i,charac)
               call fourier_1d_reverse(filename1,charac)
         case ("--gcd")
                 ! determines greatest common divisor (GCD) of a vector of integers  
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) intnum ! filename
                 allocate(intvec(intnum))
                 do j=1, intnum
                   i=i+1
                   call get_command_argument(i,charac)
                   read(charac,*) intvec(j) ! filename
                 end do 
                 print '(8x,"GDC=",I0)', gcd(intvec)
         case ('--nicecluster','--goodcluster')
            i=i+1
            call get_command_argument(i,infile)
            i=i+1
            call get_command_argument(i,informat)
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) radius  ! radius of sphere that contains cluster
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) coordminpar  ! minimum accepted coordination of atoms in cluster (e.g. 2)
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) coordmeanpar ! minimum accepted average coordination number in cluster
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) dipolepar  ! maximum accepted cluster dipole moment 
            ! determine number of species
            j=1
            if(talk) then
              print*,"        infile:",infile
              print*,"        informat:",informat
            end if
            !print*,"i:",i
            !print*, command_argument_count()," arguments"
            do while (i+j.lt.command_argument_count()) 
              call get_command_argument(i+j,charac)
              !print*, "argument",i+j,charac
              j=j+1
            end do
            ! There should be two informations per species. If not,
            ! something is wrong.
            if(mod(j,2).gt.0) then
              print(ferrmssg), "main: bad options. Consult help"
              nerr=nerr+1
              goto 4000 
            end if
            k=j/2 ! number of species
            if(talk) 
     &        print '(2x,"Your cluster will contain ",I6," species:")',k
            ! for each species read name and number of atoms to be
            ! contained in the cluster
            allocate(species(k))
            j=0
            !do while (i.lt.command_argument_count()) 
            do while (j.lt.k) 
              j=j+1 ! element index
              i=i+1 ! command line argument index
              call get_command_argument(i,charac)
              !print*,charac
              read(charac,'(A2)') species(j)%name(1:2)
              i=i+1
              call get_command_argument(i,charac)
              !print*,charac
              read(charac,*) species(j)%howmany
              if(talk) 
     &          print '(2x,A2," (",I6,")")', species(j)%name(1:2),
     &          species(j)%howmany
            end do
            call nicecluster(infile,informat,species,radius,
     &        coordminpar,coordmeanpar,dipolepar)
         case ('--nicecluster2','--goodcluster2')
            ! the same as nicecluster, with some additional things
            ! defaults
            test=.false.
            origin=0.0
            coordminpar=2
            coordmeanpar=2.0
            do while (i+1.lt.command_argument_count()) 
              i=i+1
              call get_command_argument(i,charac)
              select case(charac)
              case('-test')
                test=.true.
              case('-origin')
                do j=1,3
                  i=i+1
                  call get_command_argument(i,charac)
                  read(charac,*) origin(j)  ! origin of sphere that contains cluster and of inner core region
                end do
              case('-file')
                i=i+1
                call get_command_argument(i,infile)
              case('-format')
                i=i+1
                call get_command_argument(i,informat)
              case('-radius0')
                i=i+1
                call get_command_argument(i,charac)
                read(charac,*) radius0  ! radius of core of sphere that contains cluster
              case('-radius')
                i=i+1
                call get_command_argument(i,charac)
                read(charac,*) radius  ! radius of sphere that contains cluster
              case('-coordmin')
                i=i+1
                call get_command_argument(i,charac)
                read(charac,*) coordminpar  ! minimum accepted coordination of atoms in cluster (e.g. 2)
              case('-coordmean')
                i=i+1
                call get_command_argument(i,charac)
                read(charac,*) coordmeanpar ! minimum accepted average coordination number in cluster
              case('-dipmax')
                i=i+1
                call get_command_argument(i,charac)
                read(charac,*) dipolepar  ! maximum accepted cluster dipole moment 
              case('-atoms')
                ! determine number of species
                j=1
                do while (i+j.lt.command_argument_count()) 
                  call get_command_argument(i+j,charac)
                  j=j+1
                end do
                ! There should be two informations per species. If not,
                ! something is wrong.
                if(mod(j,2).gt.0) then
                  print(ferrmssg), "main: bad options. Consult help"
                  nerr=nerr+1
                  goto 4000 
                end if
                k=j/2 ! number of species
                if(talk) 
     &            print '(2x,"Your cluster will contain ",I6," species:"
     &)',k
                ! for each species read name and number of atoms to be
                ! contained in the cluster
                allocate(species(k))
                j=0
                !do while (i.lt.command_argument_count()) 
                do while (j.lt.k) 
                  j=j+1 ! element index
                  i=i+1 ! command line argument index
                  call get_command_argument(i,charac)
                  !print*,charac
                  read(charac,'(A2)') species(j)%name(1:2)
                  i=i+1
                  call get_command_argument(i,charac)
                  !print*,charac
                  read(charac,*) species(j)%howmany
                  if(talk) 
     &              print '(2x,A2," (",I6,")")', species(j)%name(1:2),
     &              species(j)%howmany
                end do ! while(j.lt.k)
                case default
                  print ferrmssg,"Unknown cmd line option"
                  nerr=nerr+1
                  stop
                end select
              end do
            if(talk) then
              print*,"        infile:",infile
              print*,"        informat:",informat
            end if
            call nicecluster2(infile,informat,species,radius0,radius,
     &        origin,coordminpar,coordmeanpar,dipolepar,
     &        test)
         case("--gulpnebana")
                 ! reads output from gulp 3.4 NEB calculation and writes GEOMETRY
                 ! files for all images readable by xcrysden/VESTA. 
                 ! The routine assumes that coordinates are fractional.
                 ! writes also plottable file NEBENERGIES.DAT with energies of images.
                 i=i+1
                 call get_command_argument(i,infile)
                 call gulpnebana(infile)
         case ('-h', '--help')
            if (i+1.le.command_argument_count()) then
              i=i+1
              call get_command_argument(i, morehelpon)
            else 
              morehelpon='nomorehelp'
            end if
            call print_help(morehelpon)
            stop
         case ("--hart2ev","hartree2evolt","Hart2eV",
     &         "Hartree2eVolt")
                 ! converts between Hartree and eV
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) numin
                 numout=numin*hartree
                 !print '(8x,"Your number in eV: ",1p,e12.3)',numout
                 print fnumexp,numout
         case ("--hart2ryd","hartree2rydberg","Hart2Rydt",
     &         "Hartree2Rydberg")
                 ! converts between Hartree and Rydberg
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) numin
                 numout=numin*hartree/Ryd
                 !print '(8x,"Your number in Rydberg: ",1p,e12.3)',numout
                 print fnumexp,numout
         case("--hcover","--Hcover")
               ! cover a cluster with H. Please provide the cluster file
               ! and a structure file that contains the cluster and
               ! enough embedding atoms (those will be replaced by H)
                 i=i+1
                 call get_command_argument(i,infile)
                 i=i+1
                 call get_command_argument(i,infile2)
                 i=i+1
                 call get_command_argument(i,informat)
                 call read_coords(infile,informat,atoms,natoms,species,
     &                       nspecies,vecs)
                 call read_coords(infile2,informat,atoms2,natoms2,
     &                 species2,nspecies2,vecs)
                 call Hcover(atoms,atoms2,atoms3)
                 ! get correct fractional coordinates
                 do iatom=1,size(atoms3)
                   call abs2frac(atoms3(iatom)%abswhere,vecs,
     &               atoms3(iatom)%where)
                 end do
                 call getspecies(atoms3,species3)
                 filename1=' '
                 filename1(1:5)="HYDR."
                 filename1(6:5+len_trim(infile))=trim(infile)
                 call write_coords(filename1,informat,atoms3,
     &             size(atoms3),species3,size(species3),vecs)
         case ('--input')
            i=i+1
            call get_command_argument(i, inputfile)
            if(talk) 
     &        print '(a,x,a,a)', 'Input will be read from file',
     &            trim(inputfile),'.'       
         case ('--integrate')
               ! calculates the integral of the second column wrt the
               ! first one 
            i=i+1
            call get_command_argument(i,filename1)
            ! TODO (?) : implement different integration schemes,
            ! Runge-Kutta? Leap-frog?
            !i=i+1
            !call get_command_argument(i,charac)
            charac="riemann"
            call integrate(filename1,charac)
         case ('--invmat3')
            do j=1,3
              do k=1,3
                i=i+1
                call get_command_argument(i,charac)
                read(charac,*) inmat(j,k)
              end do
            end do
            print '(8x,"matrix: ",/3(8x,3(E32.16)/))',                    &
     &           inmat(1,1:3),inmat(2,1:3),inmat(3,1:3)
            call invert_matrix(inmat,outmat)
            print '(8x,"inverse matrix: ",/3(8x,3(E32.16)/))',            &
     &           outmat(1,1:3),outmat(2,1:3),outmat(3,1:3)
            print '(8x,"inverse matrix (float): ",/3(8x,3(F32.16)/))',    &
     &           outmat(1,1:3),outmat(2,1:3),outmat(3,1:3)
         case ('--ipol')
            ! interpolates the second column wrt the
            ! first one  
            i=i+1
            call get_command_argument(i,filename1) ! filename
            i=i+1
            call get_command_argument(i,charac) ! ipol finegrid factor (finegrid=grid/step)
            read(charac,*) step
            i=i+1
            call get_command_argument(i,charac) ! ipolmeth
            call interpol(filename1,charac,step)
         case ('--isperovskite')
            ! checks whether a structure is a (distorted) perovskite structure   
            i=i+1
            call get_command_argument(i,infile) ! filename
            i=i+1
            call get_command_argument(i,informat) ! filename
            nndist=nndist_def
            if (i+1.lt.command_argument_count()) then
              i=i+1
              call get_command_argument(i,charac)
              select case(charac)
              case('-cutoff')
                i=i+1
                call get_command_argument(i,charac)
                read(charac,*) nndist
              case default
                call error_stop('unknown option')
              end select
              !
            end if ! i+1.lt.command_argument_count()
            call check_if_perovskite(infile,informat,nndist,lnum)
            print '(8x,"structure is a perovskite: ",L1)',lnum
         case ('--isperovskite_loose')
            ! checks whether a structure is a (distorted) perovskite structure.
            ! The same as is_perovskite, but with looser criteria (more
            ! distortion allowed)  
            i=i+1
            call get_command_argument(i,infile) ! filename
            i=i+1
            call get_command_argument(i,informat) ! filename
            nndist=nndist_def
            if (i+1.lt.command_argument_count()) then
              i=i+1
              call get_command_argument(i,charac)
              select case(charac)
              case('-cutoff')
                i=i+1
                call get_command_argument(i,charac)
                read(charac,*) nndist
              case default
                call error_stop('unknown option')
              end select
              !
            end if ! i+1.lt.command_argument_count()
            call check_if_perovskite_loose(infile,informat,nndist,lnum)
            print '(8x,"structure is a perovskite: ",L1)',lnum
         case ('--isperovskite_loose_analyze')
            ! The same as is_perovskite_loose, but with additional
            ! analysis  
            i=i+1
            call get_command_argument(i,infile) ! filename
            i=i+1
            call get_command_argument(i,informat) ! filename
            nndist=nndist_def
            matrix=0.0d0 ! default: axes == cell vectors (axes==0)
            do while (i+1.lt.command_argument_count()) 
              i=i+1
              call get_command_argument(i,charac)
              select case(charac)
              case('-cutoff')
                i=i+1
                call get_command_argument(i,charac)
                read(charac,*) nndist
              case('-axes')
                ! read axes wrt which to calculate octahedral tilts.
                ! Default: Supercell vectors
                if (i+9.le.command_argument_count()) then
                  do j=1,3
                    do k=1,3
                      i=i+1
                      call get_command_argument(i,charac)
                      read(charac,*) matrix(j,k)
                    end do
                  end do
                end if
              case('-autoaxes')
                ! automatically determine axes from A or B sublattice 
                if (i+1.le.command_argument_count()) then
                    i=i+1
                    call get_command_argument(i,charac)
                    if (charac.eq."A") matrix=-1.0d0
                    if (charac.eq."B") matrix=-2.0d0
                end if
              case default
                call error_stop('unknown option')
              end select
              !
            end do ! i+1.lt.command_argument_count()
            call check_if_perovskite_loose_analyze(infile,informat,       &
     &           matrix,nndist,lnum)
            print '(8x,"structure is a perovskite: ",L1)',lnum
         case('--kcalpermole')
           ! convert kcalpermole to eV per molecule
           i=i+1
           call get_command_argument(i,charac) ! number in kcal/mole
           read(charac,*) afloat
           print '(2x,F12.6,x,"eV per molecule")',afloat*kcalpermole
         case ("--latsys")
                 ! determines the crystal system
                 i=i+1
                 call get_command_argument(i,infile)
                 i=i+1
                 call get_command_argument(i,informat)
                 tol2(1)=1E-3
                 if (i+1.le.command_argument_count()) then
                   i=i+1
                   call get_command_argument(i,charac)
                   read(charac,*) tol2(1)
                 end if
                 tol2(2)=tol2(1)
                 if (i+1.le.command_argument_count()) then
                   i=i+1
                   call get_command_argument(i,charac)
                   read(charac,*) tol2(2)
                 end if
                 call read_coords(infile,informat,atoms,natoms,species,
     &                       nspecies,vecs)
                 call latsys(vecs,atoms,tol2,charac)
         case("--lmpthermo")
               ! write MD data from LAMMPS logfile to "MD.DAT"
                 i=i+1
                 call get_command_argument(i,lmplog)
                 call lmpthermo(lmplog)
         case ('--locpotav')
            ! average VASP LOCPOT file
            i=i+1
            call get_command_argument(i,infile)
            ! default direction: 3 (z)
            axisint=3
            if (i+1.le.command_argument_count()) then
              i=i+1
              call get_command_argument(i,charac)
              read(charac,*) numin
              if(numin.ge.0.99d0) then
                read(charac(1:len(charac)),*) axisint
              end if
            end if 
            ! default length for macroscopic average: 0.5 unit cell
            ! length
            rcut=0.5d0
            if (i+1.le.command_argument_count()) then
              i=i+1
              call get_command_argument(i,charac)
              read(charac,*) rcut
            end if 
            if(talk) print '(a,x,a,x,a,i0,x,a,F6.3,a)', 'File ',
     &        trim(infile),'will be averaged along x', axisint,           &
     &        ', stepsize ',rcut,' unit cell length.'
            call locpotav(infile,'poscar',axisint,rcut)
         case("--madeldipen")
               ! calculate the madelung energy for some point charge
               ! lattices in a homogeneous background charge
                 ! lattice type
                 i=i+1
                 call get_command_argument(i,lattice)
                 ! charge
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) charge
                 ! dipole vector
                 do j=1,3
                   i=i+1
                   call get_command_argument(i,charac)
                   read(charac,*) dipole(j)
                 end do
                 ! epsilon (scalar)
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) epsil
                 ! lattice constant a (if tet, also c)
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) vector(1)
                 vector(2)=vector(1)
                 vector(3)=vector(1)
                 select case (lattice)
                 case('tet','TET','tetra','TETRA')
                   i=i+1
                   call get_command_argument(i,charac)
                   read(charac,*) vector(3)
                 end select
                 call madeldipen(lattice,charge,dipole,epsil,vector)
         case("--madeldipen2")
               ! calculate the madelung energy for a point charge
               ! lattice whose basis is a dipole 
                 ! cell vectors
                 do j=1,3
                   do k=1,3
                     i=i+1
                     call get_command_argument(i,charac)
                     read(charac,*) vecs(j,k)
                   end do
                 end do
                 ! charge
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) charge
                 ! dipole vector
                 do j=1,3
                   i=i+1
                   call get_command_argument(i,charac)
                   read(charac,*) dipole(j)
                 end do
                 ! epsilon (tensor)
                 do j=1,3
                   do k=1,3
                     i=i+1
                     call get_command_argument(i,charac)
                     read(charac,*) matrix(j,k)
                   end do
                 end do
                 call madeldipen2(vecs,charge,dipole,matrix)
!         case("--madelen")
!               ! calculate the madelung energy for some point charge
!               ! lattices in a homogeneous background charge
!                 i=i+1
!                 call get_command_argument(i,lattice)
!                 i=i+1
!                 call get_command_argument(i,charac)
!                 read(charac,*) charge
!                 i=i+1
!                 call get_command_argument(i,charac)
!                 read(charac,*) epsil
!                 i=i+1
!                 call get_command_argument(i,charac)
!                 read(charac,*) vector(1)
!                 vector(2)=vector(1)
!                 vector(3)=vector(1)
!                 select case (lattice)
!                 case('tet','TET','tetra','TETRA')
!                   i=i+1
!                   call get_command_argument(i,charac)
!                   read(charac,*) vector(3)
!                 end select
!                 call madelen(lattice,charge,epsil,vector,talk,nerr)
!         case("--madelen")
!               ! calculate the madelung energy for some point charge
!               ! lattices in a homogeneous background charge
!                 ! lattice vectors
!                 do j=1,3
!                   do k=1,3
!                     i=i+1
!                     call get_command_argument(i,charac)
!                     read(charac,*) vecs(j,k)
!                   end do
!                 end do
!                 ! charge
!                 i=i+1
!                 call get_command_argument(i,charac)
!                 read(charac,*) charge
!                 ! epsilon tensor
!                 do j=1,3
!                   do k=1,3
!                     i=i+1
!                     call get_command_argument(i,charac)
!                     read(charac,*) matrix(j,k)
!                   end do
!                 end do
!                 call madelen(vecs,charge,matrix,talk,nerr)
         case("--madelen1")
               ! calculate the madelung energy for some point charge
               ! lattices in a homogeneous background charge
                 ! lattice vectors
                 do j=1,3
                   do k=1,3
                     i=i+1
                     call get_command_argument(i,charac)
                     read(charac,*) vecs(j,k)
                   end do
                 end do
                 ! charge
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) charge
                 ! epsilon tensor
                 do j=1,3
                   do k=1,3
                     i=i+1
                     call get_command_argument(i,charac)
                     read(charac,*) matrix(j,k)
                   end do
                 end do
                 call madelen1(vecs,charge,matrix)
         case("--madelen")
               ! calculate the madelung energy for any point charge
               ! lattice in a medium with any dielectic tensor
                 i=i+1
                 call get_command_argument(i,infile)
                 i=i+1
                 call get_command_argument(i,informat)
                 ! epsilon tensor
                 do j=1,3
                   do k=1,3
                     i=i+1
                     call get_command_argument(i,charac)
                     read(charac,*) matrix(j,k)
                   end do
                 end do
                 call read_coords(infile,informat,atoms,natoms,species,
     &                       nspecies,vecs)
                 call madelen(vecs,atoms,matrix,energy)
                 if(.not.talk) then 
                   print '(8x,"Madelung energy in Ryd=",F20.10)',energy
                   print '(8x,"Madelung energy in eV=",F20.10)',          &
     &                    energy*Ryd
                 end if          
         case ('--matchatoms')
            ! input: two geometris atoms1 and atoms2. reorders 
            ! atoms2 such that their positions match best those of
            ! atoms1
            i=i+1
            call get_command_argument(i,infile)
            i=i+1
            call get_command_argument(i,format1)
            i=i+1
            call get_command_argument(i,infile2)
            i=i+1
            call get_command_argument(i,format2)
            i=i+1
            call get_command_argument(i,outfile)
            i=i+1
            call get_command_argument(i,format3)
            call read_coords(infile,format1,atoms,size(atoms),            &
     &           species,size(species),vecs)
            call read_coords(infile2,format2,atoms2,size(atoms2),         &
     &           species2,size(species2),vecs2)
            call matchatoms(atoms,atoms2,vecs,vecs2)
            call write_coords(outfile,format3,atoms2,size(atoms2),        &
     &                        species2,size(species2),vecs2)
         case ('--mateigs')
            ! diagonalize real nxn matrix
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) intnum ! dimension of matrix
            print '(8x,"matrix dimension: ",I0," x ",I0)',intnum,intnum 
            if (allocated(allomat1)) deallocate(allomat1)
            allocate(allomat1(intnum,intnum))
            print '(8x,"Matrix:")'
            format1=''
            write(format1,*) "(8x,",intnum,"(F20.10,1x))"
            do k=1,intnum
              do j=1,intnum
                i=i+1
                call get_command_argument(i,charac)
                read(charac,*) allomat1(k,j) ! matrix
              end do
              print format1, allomat1(k,1:intnum)
            end do
            if (allocated(allovec1)) deallocate(allovec1)
            call get_mateigs(allomat1,allovec1,.true.)
            print '(8x,"Eigenvalues:")'
            do k=1,size(allovec1)
              print '(8x,E20.10)', allovec1(k)
            end do    
            print '(8x,"Eigenvalues (float):")'
            do k=1,size(allovec1)
              print '(8x,F20.10)', allovec1(k)
            end do    
         case("--matrix")
               ! multiply all coordinates with a matrix
               i=i+1
               call get_command_argument(i,infile)
               i=i+1
               call get_command_argument(i,outfile)
               i=i+1
               call get_command_argument(i,informat)
               do j=1,3
                 do k=1,3
                   i=i+1
                   call get_command_argument(i,matchar)
                   read(matchar(1:20),*) matrix(j,k)
                 end do
               end do
               call read_coords(infile,informat,atoms,natoms,species,
     &                     nspecies,vecs)
               call matrixmult(atoms,vecs,matrix)
               call write_coords(outfile,informat,atoms,natoms,species,
     &                     nspecies,vecs)
         case("--matrix_mult")
                ! 
                ! multiplies (n,n) matrix A with (n,n) matrix B: C=AB 
                ! 
                ! read dimension    
                i=i+1
                call get_command_argument(i,charac)
                read(charac,*) intnum
                if (allocated(allomat1)) deallocate(allomat1)
                if (allocated(allomat2)) deallocate(allomat2)
                allocate(allomat1(intnum,intnum))
                allocate(allomat2(intnum,intnum))
                allocate(allomat3(intnum,intnum))
                !
                ! begin read A
                !
                do j=1,intnum
                  do k=1,intnum
                    i=i+1
                    call get_command_argument(i,charac)
                    read(charac,*) allomat1(j,k)
                  end do ! k
                end do ! j
                !
                ! end read A
                !
                !
                ! begin read B
                !
                do j=1,intnum
                  do k=1,intnum
                    i=i+1
                    call get_command_argument(i,charac)
                    read(charac,*) allomat2(j,k)
                  end do ! k
                end do ! j
                !
                ! end read B
                !
                call matmul_fortran(intnum,allomat1,allomat2,allomat3)
                !
                print '(8x,"AB=")'
                do j=1,intnum
                  print '(8x,3(F10.6))',allomat3(j,1:3)
                end do  
                !
         case("--mbppnebana")
                 ! reads output from MBPP NEB calculation and writes GEOMETRY
                 ! files for all images readable by xcrysden/VESTA. 
                 ! The routine assumes that coordinates are fractional.
                 ! writes also plottable file NEBENERGIES.DAT with energies of images.
                 i=i+1
                 call get_command_argument(i,infile)
                 infile=trim(adjustl(infile))
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) nimg
                 !call mbppnebana(infile,nimg,ncomm,nwarn,nerr)
         case ('--mdana')
            i=i+1
            call get_command_argument(i,trjfile)
            i=i+1
            call get_command_argument(i,trjformat)
            if(trjformat.eq.'trg') then
              i=i+1
              call get_command_argument(i,trjatomfile)
            else
                  trjatomfile='notneeded'    
            end if
            trjbndry=tolfrac
            ! read boundary parameters, if specified
            if (i+1.le.command_argument_count()) then
                  i=i+1
                  call get_command_argument(i,bndrychar)
                  select case (bndrychar)
                  case("-bndry")
                      if (i+3.le.command_argument_count()) then
                        do j=1,3
                          i=i+1
                          call get_command_argument(i,vecchar)
                          read(vecchar(1:20),*) trjbndry(j)
                        end do
                      end if
                  case default  
                      if(talk) then
                        print '(A)', "No boundary parameter found."
                        print '(A)', "Default param. will be used."
                      end if
                      i=i-1
                  end select
            end if
            if(talk) then 
              print'(A)', "Will use as boundary parameters:"
              print'(3(F10.6))', trjbndry
              print '(a,x,a,x,a)', 'File',trim(trjfile),
     &            'will be analyzed.'
            end if
            call mdana(trjfile,trjformat,trjatomfile,trjbndry)
         case("--merge")
            i=i+1
            call get_command_argument(i,filename1)
            i=i+1
            call get_command_argument(i,format1)
            i=i+1
            call get_command_argument(i,filename2)
            i=i+1
            call get_command_argument(i,format2)
            i=i+1
            call get_command_argument(i,filename3)
            i=i+1
            call get_command_argument(i,format3)
            if(talk) 
     &        print '(a,x,a,x,a)', 'File',trim(filename1),
     &            ' will be merged with file',trim(filename2),'.'
            call mergefiles(filename1,format1,filename2,format2,
     &            filename3,format3)
         case("--mult")
            !if (i+3.le.command_argument_count()) then
            do j=1,3
               i=i+1
               call get_command_argument(i,intchar)
               read(intchar(1:20),*) nsc(j)
               !matrix=0.0D0
               !matrix(1,1)=1.0D0
               !matrix(2,2)=1.0D0
               !matrix(3,3)=1.0D0
            end do
            i=i+1
            call get_command_argument(i, fincoords)
            i=i+1
            call get_command_argument(i, foutcoords)
            i=i+1
            call get_command_argument(i, informat)
            i=i+1
            call get_command_argument(i, foutform)
            if(talk) then
              print '("Will read coordinates from file ",A)',
     &             trim(fincoords)
              print '("and write output to file ",A,".")',
     &             trim(foutcoords)
              print '("Format of input file: ",A)',
     &             trim(informat)
              print '("Format of output file: ",A)',
     &             trim(foutform)
            end if
            ! read supercell dimension, if specified
            !if (i+1.le.command_argument_count()) then
            !      i=i+1
            !      call get_command_argument(i,transoption)
            !      select case (transoption)
            !      case("-mult")
            !          if (i+3.le.command_argument_count()) then
            !            do j=1,3
            !              i=i+1
            !              call get_command_argument(i,intchar)
            !              read(intchar(1:20),*) nsc(j)
            !            end do
            !          end if
            !      case default  
            !          print '(A)', "No trans option found."
            !          i=i-1
            !      end select
            !end if
            if(talk) then
              print'(A)', "Final supercell dimension:"
              print'(3(I4))', nsc
            end if
            call multi(fincoords,foutcoords,
     &           informat,foutform,nsc)
         case ('--multiply')
                 ! multiplies a number by another number
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) number1
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) number2
!                 print '(8x,"Your result: ",F)',number1*number2
!                 print '(8x,"Your result: ",1p,e12.3)',
!     &                number1*number2
                 print fnumexp,number1*number2 
         case("--nebana")
                 ! reads output from a NEB calculation and writes GEOMETRY
                 ! files for all images readable by xcrysden/VESTA. 
                 ! The routine assumes that coordinates are fractional.
                 ! writes also plottable file NEBENERGIES.DAT with energies of images.
                 i=i+1
                 call get_command_argument(i,infile)
                 i=i+1
                 call get_command_argument(i,informat)
                 call nebana(infile,informat)
         case("--nebchain")
                 ! reads a chain of structure files and sets up an input file for a NEB calculation.
                 i=i+1
                 call get_command_argument(i,infile)
                 i=i+1
                 call get_command_argument(i,informat)
                 i=i+1
                 call get_command_argument(i,outformat)
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) nimg
                 call nebchain(infile,informat,outformat,
     &            nimg)
         case ('--output')
            i=i+1
            call get_command_argument(i, outputfile)
            if(talk) 
     &       print '(a,x,a,a)', 'Output will be written to file',
     &            trim(outputfile),'.'    
         case ('--overlap','--overlapwf')
            i=i+1
            call get_command_argument(i,filename1)
            i=i+1
            call get_command_argument(i,filename2)
            i=i+1
            call get_command_argument(i,filename3)
            call overlapwf(filename1,filename2,filename3)
         case ('--phonopy_print_bands')
            ! reads phonon band structure from specified file (e.g.
            ! band.yaml) and prints phonon eigenvalues and if
            ! present eigenvectors to specified outputfile
            i=i+1
            call get_command_argument(i,filename1)
            i=i+1
            call get_command_argument(i,filename2)
            call phonopy_print_bands(filename1,filename2)
         case ('--phonopy_print_eigenvecs')
            ! reads phonon eigenvectors from specified file (e.g.
            ! band.yaml) and prints xsf file with phonon eigenvectors 
            ! to specified outputfile
            i=i+1
            call get_command_argument(i,filename1)
            i=i+1
            call get_command_argument(i,filename2)
            intnum=1
            intnum2=1
            do while (i+2.le.command_argument_count()) 
              i=i+1
              call get_command_argument(i,charac)
              select case(charac)
              case('-qpoint')
                i=i+1
                call get_command_argument(i,charac)
                read(charac,*) intnum
              case('-band')
                i=i+1
                call get_command_argument(i,charac)
                read(charac,*) intnum2
              case default
              end select
            end do
            call phonopy_print_eigenvecs(filename1,filename2,intnum,      &
     &           intnum2)
         case ('--power')
                 ! takes a number to the power of another number
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) number1
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) number2
!                 print '(8x,"Your result: ",F)',number1**number2
!                 print '(8x,"Your result: ",1p,e20.11)',
!     &                number1**number2
                 print fnumexp,number1**number2
         case ('--pol')
            i=i+1
            call get_command_argument(i,infile)
            i=i+1
            call get_command_argument(i,informat)
            ! calculate electric polarization from positions and charges 
            ! given in the input file.
            call get_pol(infile,informat)
         case ('--pol_vol','--polyhedral_volume')
           ! calculates polyhedral volume. Polyhedron facets 
           ! are read from file "HUlL_FACETS" (see quickhull).
           i=i+1
           call get_command_argument(i,charac)
           read(charac,*) intnum ! dimension 
           if (.not.intnum==3) call error_stop('Only implemented for 3D'  &
     &)
           call polyhedral_volume(intnum)
         case('--quickhull')
           ! determine the convex hull of a set of points in 3d using the 
           ! algorithm of Barber, Barber, Dobkin, Huhdanpaa, ACM
           ! Transactions on Mathematical software volume 22 no 4, pages
           ! 469-483, 1996
            i=i+1
            call get_command_argument(i,filename1)
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) intnum ! dimension
            call quickhull(filename1,intnum)
         case('--quickhull2')
           ! determine the convex hull of a set of points in 3d using the 
           ! algorithm of Barber, Barber, Dobkin, Huhdanpaa, ACM
           ! Transactions on Mathematical software volume 22 no 4, pages
           ! 469-483, 1996
            i=i+1
            call get_command_argument(i,filename1)
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) intnum ! dimension
            call quickhull2(filename1,intnum)
         case("--randisp")
            ! apply random displacements to the atom positions   
            i=i+1
            call get_command_argument(i,infile)
            i=i+1
            call get_command_argument(i,informat)
            i=i+1
            call get_command_argument(i,outfile)
            i=i+1
            call get_command_argument(i,outformat)
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) shift(1) ! (frac) maximum fractional displacement (scaling factor)
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) shift(2) ! (frac) maximum fractional displacement (scaling factor)
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) shift(3) ! (frac) maximum fractional displacement (scaling factor)
            call read_coords(infile,informat,atoms,natoms,species,
     &                     nspecies,vecs)
            !call randisp(infile,informat,outfile,outformat,shift)
            call randisp(atoms,vecs,shift)
            call write_coords(outfile,outformat,atoms,natoms,species,
     &                     nspecies,vecs)
         case("--rdf")
            ! calculate radial distribution functions
            i=i+1
            call get_command_argument(i,infile)
            i=i+1
            call get_command_argument(i,informat)
            ! read lattice vectors, if needed
            select case(informat)
            case("coorat","COORAT",'nwchem','NWCHEM','xyz','XYZ')
              if (i+9.le.command_argument_count()) then
                do j=1,3
                  i=i+1
                  call get_command_argument(i,vecchar)
                  read(vecchar(1:20),*)vecs(j,1)
                  i=i+1
                  call get_command_argument(i,vecchar)
                  read(vecchar(1:20),*)vecs(j,2)
                  i=i+1
                  call get_command_argument(i,vecchar)
                  read(vecchar(1:20),*)vecs(j,3)
                end do
             else  
               call error_stop('Please provide lattice vectors after for  &
     &mat, or use a format that contains lattice vectors')
             end if  
             case default
             end select
             i=i+1
             call get_command_argument(i,charac)
             read(charac,*)rmax
             i=i+1
             call get_command_argument(i,charac)
             read(charac,*) binwidth
             i=i+1
             call get_command_argument(i,charac)
             read(charac,*) broad
             call rdf(infile,informat,rmax,binwidth,broad,vecs)
         case("--read_symops_vasp")
               ! read Symmetry operations from vasp OUTCAR 
               call read_symops_vasp()
         case("--readxcme")
               ! read xc matrix elements from nwchem output file
               i=i+1
               call get_command_argument(i,infile)
               i=i+1
               call get_command_argument(i,informat)
               i=i+1
               call get_command_argument(i,charac)
               call readxcme(infile,informat)
         case ('--relax_atoms')
               ! move a bunch of atoms along forces 
            i=i+1
            call get_command_argument(i,infile)
            i=i+1
            call get_command_argument(i,informat)
            i=i+1
            call get_command_argument(i,outfile)
            i=i+1
            call get_command_argument(i,outformat)
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) number1
            call read_coords(infile,informat,atoms,natoms,species,
     &                     nspecies,vecs)
            call relax_atoms(atoms,vecs,number1)
            call write_coords(outfile,outformat,atoms,natoms,species,
     &                     nspecies,vecs)
         case ('--relax_rigid_molecule')
               ! move a bunch of atoms rigidly (translation and rotation
               ! of a molecule)
            i=i+1
            call get_command_argument(i,infile)
            i=i+1
            call get_command_argument(i,informat)
            i=i+1
            call get_command_argument(i,outfile)
            i=i+1
            call get_command_argument(i,outformat)
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) number1 ! this is the step size for relaxation
            call read_coords(infile,informat,atoms,natoms,species,
     &                     nspecies,vecs)
            call relax_rigid_molecule(atoms,vecs,number1)
            call write_coords(outfile,outformat,atoms,natoms,species,
     &                     nspecies,vecs)
         case ('--relax_rigid_molecule_trala')
               ! relax a bunch of atoms rigidly along forces (rigid translation only)
            i=i+1
            call get_command_argument(i,infile)
            i=i+1
            call get_command_argument(i,informat)
            i=i+1
            call get_command_argument(i,outfile)
            i=i+1
            call get_command_argument(i,outformat)
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) number1 ! this is the step size for relaxation
            call read_coords(infile,informat,atoms,natoms,species,
     &                     nspecies,vecs)
            call relax_rigid_molecule_trala(atoms,vecs,number1)
            call write_coords(outfile,outformat,atoms,natoms,species,
     &                     nspecies,vecs)
         case ('--relax_rigid_molecule_rot')
               ! relax a bunch of atoms rigidly along forces (rigid
               ! rotation only)
            i=i+1
            call get_command_argument(i,infile)
            i=i+1
            call get_command_argument(i,informat)
            i=i+1
            call get_command_argument(i,outfile)
            i=i+1
            call get_command_argument(i,outformat)
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) number1 ! this is the step size for relaxation
            call read_coords(infile,informat,atoms,natoms,species,
     &                     nspecies,vecs)
            call relax_rigid_molecule_rot(atoms,vecs,number1)
            call write_coords(outfile,outformat,atoms,natoms,species,
     &                     nspecies,vecs)
         case ('--relax_rigid_molecule_rot2')
               ! relax a bunch of atoms rigidly along forces (rigid
               ! rotation only). Same as above, but different averageing
            i=i+1
            call get_command_argument(i,infile)
            i=i+1
            call get_command_argument(i,informat)
            i=i+1
            call get_command_argument(i,outfile)
            i=i+1
            call get_command_argument(i,outformat)
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) number1 ! this is the step size for relaxation
            call read_coords(infile,informat,atoms,natoms,species,
     &                     nspecies,vecs)
            call relax_rigid_molecule_rot2(atoms,vecs,number1)
            call write_coords(outfile,outformat,atoms,natoms,species,
     &                     nspecies,vecs)
         case ('--relax_rigid_molecule_breath')
               ! relax a bunch of atoms along forces (breathing mode only)
            i=i+1
            call get_command_argument(i,infile)
            i=i+1
            call get_command_argument(i,informat)
            i=i+1
            call get_command_argument(i,outfile)
            i=i+1
            call get_command_argument(i,outformat)
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) number1 ! this is the step size for relaxation
            call read_coords(infile,informat,atoms,natoms,species,
     &                     nspecies,vecs)
            call relax_rigid_molecule_breath(atoms,vecs,number1)
            call write_coords(outfile,outformat,atoms,natoms,species,
     &                     nspecies,vecs)
         case("--replace")
               ! replace atoms of one species by atoms of a different
               ! species
               i=i+1
               call get_command_argument(i,infile)
               i=i+1
               call get_command_argument(i,informat)
               i=i+1
               call get_command_argument(i,thespecies)
               i=i+1
               call get_command_argument(i,charac)
               read(charac,*) thenumber
               i=i+1
               call get_command_argument(i,theotherspecies)
               theotherspecies=adjustl(theotherspecies)
               i=i+1
               call get_command_argument(i,outfile)
               call rpatoms(infile,informat,thespecies,thenumber,
     &              theotherspecies,outfile)
         case("--rmeqv")
               ! remove all atoms but one in each group of atoms that are equivalent by symmetry
                 i=i+1
                 call get_command_argument(i,infile)
                 i=i+1
                 call get_command_argument(i,informat)
                 i=i+1
                 call get_command_argument(i,symmformat)
                 i=i+1
                 call get_command_argument(i,outfile)
                 tol1=1.0D-5
                 do while (i+1.le.command_argument_count()) 
                   i=i+1
                   call get_command_argument(i,charac)
                   read(charac,*) tol1
                 end do
                 call rmeqv(infile,informat,symmformat,outfile,tol1)
         case("--rotate")
               ! rotate all coordinates by a chosen angle around a chosen axis
                i=i+1
                call get_command_argument(i,infile)
                i=i+1
                call get_command_argument(i,outfile)
                i=i+1
                call get_command_argument(i,informat)
                i=i+1
                call get_command_argument(i,charac)
                read(charac,*) axisint
                if (axisint==0) then
                  i=i+1
                  call get_command_argument(i,charac)
                  read(charac,*) vector(1)
                  i=i+1
                  call get_command_argument(i,charac)
                  read(charac,*) vector(2)
                  i=i+1
                  call get_command_argument(i,charac)
                  read(charac,*) vector(3)
                end if
                i=i+1
                call get_command_argument(i,charac)
                read(charac,*) angle
                call read_coords(infile,informat,atoms,natoms,species,
     &                      nspecies,vecs)
                if(axisint==1.or.axisint==2.or.axisint==3) then         
                  call rotatecoords(atoms,vecs,axisint,angle)
                end if
                if(axisint==0) then         
                  call rotatecoords(atoms,vecs,axisint,angle,             &
     &                 vector=vector)
                end if
                call write_coords(outfile,informat,atoms,natoms,species,
     &                       nspecies,vecs)
         case("--rotate_matrix")
                ! 
                ! rotates matrix A to new coordinate system: B = S⁺ A S 
                ! 
                ! read dimension    
                i=i+1
                call get_command_argument(i,charac)
                read(charac,*) intnum
                if (allocated(allomat1)) deallocate(allomat1)
                if (allocated(allomat2)) deallocate(allomat2)
                allocate(allomat1(intnum,intnum))
                allocate(allomat2(intnum,intnum))
                allocate(allomat3(intnum,intnum))
                !
                ! begin read A
                !
                do j=1,intnum
                  do k=1,intnum
                    i=i+1
                    call get_command_argument(i,charac)
                    read(charac,*) allomat1(j,k)
                  end do ! k
                end do ! j
                !
                ! end read A
                !
                !
                ! begin read S
                !
                do j=1,intnum
                  do k=1,intnum
                    i=i+1
                    call get_command_argument(i,charac)
                    read(charac,*) allomat2(j,k)
                  end do ! k
                end do ! j
                !
                ! end read S
                !
                call rotate_matrix(intnum,allomat1,allomat2,allomat3)
                !
                print '(8x,"Rotated matrix:")'
                do j=1,intnum
                  print '(8x,3(F10.6))',allomat3(j,1:3)
                end do  
                !
         case ("--ryd2ev","rydberg2evolt","Ryd2eV",
     &         "Rydberg2eVolt")
                 ! converts between Rydberg and eV
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) numin
                 numout=numin*Ryd
                 !print '(8x,"Your number in eV: ",1p,e12.3)',numout
                 print fnumexp,numout
         case ("--ryd2hart","rydberg2hartree","Ryd2Hart",
     &         "Rydberg2Hartree")
                 ! converts between Rydberg and Hartree
                 i=i+1
                 call get_command_argument(i,charac)
                 read(charac,*) numin
                 numout=numin*Ryd/hartree
                 !print '(8x,"Your number in Hartree: ",1p,e12.3)',numout
                 print fnumexp,numout
         case ('--silent')  
            talk=.false.
         case ('--slidwinav')
            ! sliding window average of data in file
            i=i+1
            call get_command_argument(i,infile)
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) rcut ! length to average 
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) intnum ! number of columns in the file
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) intnum2 ! the number of the column to average
            call slidwinav(infile,rcut,intnum,intnum2)
         case ('--sortdata')
            ! sorts numbers in a data file
            i=i+1
            call get_command_argument(i,infile)
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) intnum ! number of columns in the file
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) intnum2 ! the number of the column to order 
            call sortdata(infile,intnum,intnum2)
         case("--spectrum")
                 ! extracts excitation energies from a nwchem output file.
                 i=i+1
                 call get_command_argument(i,infile)
                 i=i+1
                 call get_command_argument(i,informat)
                 ! default values
                 broad=0.05D0
                 binwidth=0.01D0 
                 do while (i+2.le.command_argument_count()) 
                   i=i+1
                   call get_command_argument(i,charac)
                   select case(charac)
                   case("-brd")
                     i=i+1
                     call get_command_argument(i,charac)
                     read(charac,*) broad
                   case("-bin")
                     i=i+1
                     call get_command_argument(i,charac)
                     read(charac,*) binwidth
                   case default
                     print '("Unknown option. Will be ignored.")'
                   end select 
                 end do
                 if(talk)  then
                   print '(2x,"Calling spectrum with the following param
     &eters:")'
                   print '(2x,"broadening=",F10.7)',broad
                   print '(2x,"binwidth=",F10.7)',binwidth
                 end if
                 call spectrum(infile,informat,
     &            broad,binwidth)
         case ('-t', '--time')
            call date_and_time(DATE=date, TIME=time)
            write (*, '(a,x,a,"-",a,"-",a,x,a,":",a,a)')
     &          'The current date and time is',date(1:4), 
     &          date(5:6),date(7:8),time(1:2),time(3:4),'.'
            stop
         case("--traj")
                ! get trajectory of single atom
                i=i+1
                call get_command_argument(i,infile)
                i=i+1
                call get_command_argument(i,informat)
                i=i+1
                call get_command_argument(i,charac)
                read(charac(1:20),*) firststep
                i=i+1
                call get_command_argument(i,charac)
                read(charac(1:20),*) step
                i=i+1
                call get_command_argument(i,charac)
                read(charac(1:20),*) laststep
                i=i+1
                call get_command_argument(i,charac)
                read(charac(1:20),*) thespecies
                thespecies=trim(thespecies)
                i=i+1
                call get_command_argument(i,charac)
                ! if a single atom position is to be analyzed, read the
                ! atom number ("thenumber), if all atoms of this species are to
                ! analyzed, "thenumber" is -1.
                if(index('0123456789',charac(1:1)).gt.0) then
                   read(charac(1:20),*) thenumber
                else
                   if(trim(adjustl(charac)).eq.'all') thenumber=-1
                end if
                i=i+1
                call get_command_argument(i,charac)
                read(charac(1:20),*) nndist
                ! determine if site type(octa or tetra) should be
                ! determined (default: no)
                nnana="none"
                nnatom="25"
                if (i+3.le.command_argument_count()) then
                  i=i+1
                  call get_command_argument(i,charac)
                  charac=trim(adjustl(charac))
                  if(charac.eq."-nnana") then
                        i=i+1
                        call get_command_argument(i,nnana)
                        nnana=trim(adjustl(nnana))
                        ! read species that forms the octa-/tetrahedra, e.g.
                        ! O
                        i=i+1
                        call get_command_argument(i,nnatom)
                        nnatom=trim(adjustl(nnatom))
                  else
                        i=i-1
                  end if
                end if
                call traj(infile,informat,firststep,step,laststep,
     &               thespecies,thenumber,nndist,nnana,nnatom)
         case ('--transform')
            i=i+1
            call get_command_argument(i, trans)
            if(talk)
     &       print '("Will apply a ",a," to the coordinates.")',
     &             trim(trans)
            select case (trans)
            case("vectors")
                  nsc=1
                  do j=1,3
                    do k=1,3
                      i=i+1
                      call get_command_argument(i,vecchar)
                      read(vecchar(1:20),*) invecs(j,k)
                    end do
                  end do
                  do j=1,3
                    do k=1,3
                      i=i+1
                      call get_command_argument(i,vecchar)
                      read(vecchar(1:20),*) outvecs(j,k)
                    end do
                  end do
                  i=i+1
                  call get_command_argument(i, fincoords)
                  i=i+1
                  call get_command_argument(i, foutcoords)
                  i=i+1
                  call get_command_argument(i, informat)
                  if(talk) then
                    print '("Will read coordinates from file ",A)',
     &                     trim(fincoords)
                    print '("and write output to file ",A,".")',
     &                     trim(foutcoords)
                    print '("Format of input file: ",A)',
     &                   trim(informat)
                  end if
                  ! read supercell dimension, if specified
                  if (i+1.le.command_argument_count()) then
                        i=i+1
                        call get_command_argument(i,transoption)
                        select case (transoption)
                        case("-mult")
                            if (i+3.le.command_argument_count()) then
                              do j=1,3
                                i=i+1
                                call get_command_argument(i,intchar)
                                read(intchar(1:20),*) nsc(j)
                              end do
                            end if
                        case default  
                            if(talk) 
     &                        print '(A)', "No trans option found."
                            i=i-1
                        end select
                  end if
                  if(talk) then
                    print'(A)', "Final supercell dimension:"
                    print'(3(I4))', nsc
                  end if
                  call read_coords(fincoords,informat,atoms,natoms,       &
     &                             species,nspecies,vecs)
                  call transformvectors(invecs,outvecs,atoms,atoms2,nsc)
                  call getspecies(atoms2,species)
                  call write_coords(foutcoords,informat,atoms2,           &
     &                    size(atoms2),species,size(species),outvecs)
            case("matrix","mult","sc2hex","cf2hex","hex2sc")
                  nsc=1
                  select case (trans)
                  case("matrix")
                      do j=1,3
                        do k=1,3
                          i=i+1
                          call get_command_argument(i,matchar)
                          read(matchar(1:20),*) matrix(j,k)
                        end do
                      end do
                  case("vectors")
                      do j=1,3
                        do k=1,3
                          i=i+1
                          call get_command_argument(i,vecchar)
                          read(vecchar(1:20),*) invecs(j,k)
                        end do
                      end do
                      do j=1,3
                        do k=1,3
                          i=i+1
                          call get_command_argument(i,vecchar)
                          read(vecchar(1:20),*) outvecs(j,k)
                        end do
                      end do
                      !do n=1,3
                      !  read(10,*,end=150) invecs(n,1:3)
                      !end do
                      !do n=1,3
                      !  read(10,*,end=150) outvecs(n,1:3)
                      !end do
                      do j=1,3
                        inmat(1:3,j)=invecs(j,1:3)
                        outmat(1:3,j)=outvecs(j,1:3)
                      end do
                      call invert_matrix(outmat,outimat)
                      matrix=matmul(outimat,inmat)
                  case("mult")
                      !if (i+3.le.command_argument_count()) then
                        do j=1,3
                          i=i+1
                          call get_command_argument(i,intchar)
                          read(intchar(1:20),*) nsc(j)
                          matrix=0.0D0
                          matrix(1,1)=1.0D0
                          matrix(2,2)=1.0D0
                          matrix(3,3)=1.0D0
                        end do
                      !else
                      !    print*, "Error:missing supercell dimensions."
                      !    stop
                      !end if
                  case("sc2hex")
                      !matrix(1,1)=2.0D0
                      !matrix(1,2)=2.0D0
                      !matrix(1,3)=-4.0D0
                      !matrix(2,1)=-4.0D0
                      !matrix(2,2)=2.0D0
                      !matrix(2,3)=2.0D0
                      !matrix(3,1)=1.0D0
                      !matrix(3,2)=1.0D0
                      !matrix(3,3)=1.0D0
                      !matrix=matrix/3.0D0
                  !    matrix(1,1)=1.4142140D0
                  !    matrix(1,2)=-0.8164970D0
                  !    matrix(1,3)=0.000000D0
                  !    matrix(2,1)=0.000000D0
                  !    matrix(2,2)=1.632993D0
                  !    matrix(2,3)=0.000000D0
                  !    matrix(3,1)=0.000000D0
                  !    matrix(3,2)=0.000000D0
                  !    matrix(3,3)=0.577350D0
                      ! sc : 
                      matrix(1,1)=-2.0D0/3.0D0
                      matrix(1,2)=1.0D0/3.0D0
                      matrix(1,3)=1.0D0/3.0D0
                      matrix(2,1)=1.0D0/3.0D0
                      matrix(2,2)=1.0D0/3.0D0
                      matrix(2,3)=-2.0D0/3.0D0
                      matrix(3,1)=1.0D0/3.0D0
                      matrix(3,2)=1.0D0/3.0D0
                      matrix(3,3)=1.0D0/3.0D0
                      !matrix=matrix*2.0D0 !/3.0D0
                  case("cf2hex")
                      ! fcc in conventional sc unit cell (cf): 
                      matrix(1,1)=-4.0D0/3.0D0
                      matrix(1,2)=2.0D0/3.0D0
                      matrix(1,3)=2.0D0/3.0D0
                      matrix(2,1)=2.0D0/3.0D0
                      matrix(2,2)=2.0D0/3.0D0
                      matrix(2,3)=-4.0D0/3.0D0
                      matrix(3,1)=1.0D0/3.0D0
                      matrix(3,2)=1.0D0/3.0D0
                      matrix(3,3)=1.0D0/3.0D0
                  case("hex2sc")
                      ! hex to simple cubic
                      matrix(1,:)=(/-1.0D0,-0.00D0,1.0D0/)
                      matrix(2,:)=(/1.0D0, 1.0D0, 1.00D0 /)
                      matrix(3,:)=(/0.0D0,  -1.000D0, 1.0D0/)
                  end select
                  if(talk) then
                  print'("Will multiply coordinates with this matrix:")'
                    do j=1,3
                      print '(3(F10.6))',matrix(j,1),matrix(j,2),
     &                   matrix(j,3)
                    end do
                  end if
                  i=i+1
                  call get_command_argument(i, fincoords)
                  i=i+1
                  call get_command_argument(i, foutcoords)
                  i=i+1
                  call get_command_argument(i, informat)
                  if(talk) then
                    print '("Will read coordinates from file ",A)',
     &                     trim(fincoords)
                    print '("and write output to file ",A,".")',
     &                     trim(foutcoords)
                    print '("Format of input file: ",A)',
     &                   trim(informat)
                  end if
                  ! read supercell dimension, if specified
                  if (i+1.le.command_argument_count()) then
                        i=i+1
                        call get_command_argument(i,transoption)
                        select case (transoption)
                        case("-mult")
                            if (i+3.le.command_argument_count()) then
                              do j=1,3
                                i=i+1
                                call get_command_argument(i,intchar)
                                read(intchar(1:20),*) nsc(j)
                              end do
                            end if
                        case default  
                            if(talk) 
     &                        print '(A)', "No trans option found."
                            i=i-1
                        end select
                  end if
                  if(talk) then
                    print'(A)', "Final supercell dimension:"
                    print'(3(I4))', nsc
                  end if
                  call read_coords(fincoords,informat,atoms,natoms,       &
     &                             species,nspecies,vecs)
                  call transformmatrix(matrix,atoms,vecs,atoms2,          &
     &                                 outvecs,nsc)
                  call getspecies(atoms2,species)
                  call write_coords(foutcoords,informat,atoms2,           &
     &                    size(atoms2),species,nspecies,outvecs)
            case("shift")
                  do j=1,3
                    i=i+1
                    call get_command_argument(i,vecchar)
                    read(vecchar(1:20),*) shift(j)
                  end do  
                  if(talk) then
                    print'(A)', "Will shift coordinates by this vector (
     &x -> x+v):"
                    print '(3(F10.6))',shift(1:3)
                  end if
                  i=i+1
                  call get_command_argument(i, fincoords)
                  i=i+1
                  call get_command_argument(i, foutcoords)
                  i=i+1
                  call get_command_argument(i, informat)
                  if(talk) then
                    print '("Will read coordinates from file ",A)',
     &                   trim(fincoords)
                    print '("and write output to file ",A,".")',
     &                   trim(foutcoords)
                    print '("Format of input file: ",A)',
     &                   trim(informat)
                  end if
                  call read_coords(fincoords,informat,atoms,natoms,       &
     &                             species,nspecies,vecs)     
                  call transformshift(shift,atoms,vecs,atoms2)
                  call write_coords(foutcoords,informat,atoms2,           &
     &                    size(atoms2),species,nspecies,vecs)
            case("shiftabs")
                  do j=1,3
                    i=i+1
                    call get_command_argument(i,vecchar)
                    read(vecchar(1:20),*) shift(j)
                  end do  
                  if(talk) then
                    print'(A)', "Will shift coordinates by this vector (
     &x -> x+v):"
                    print '(3(F10.6))',shift(1:3)
                  end if
                  i=i+1
                  call get_command_argument(i, fincoords)
                  i=i+1
                  call get_command_argument(i, foutcoords)
                  i=i+1
                  call get_command_argument(i, informat)
                  if(talk) then
                    print '("Will read coordinates from file ",A)',
     &                   trim(fincoords)
                    print '("and write output to file ",A,".")',
     &                   trim(foutcoords)
                    print '("Format of input file: ",A)',
     &                   trim(informat)
                  end if
                  call read_coords(fincoords,informat,atoms,natoms,       &
     &                             species,nspecies,vecs)     
                  call transformshiftabs(shift,atoms,vecs,atoms2)
                  call write_coords(foutcoords,informat,atoms2,           &
     &                    size(atoms2),species,nspecies,vecs)
            case("sort")
                  if(talk) print'(A)', "Will sort coordinates."
                  i=i+1
                  call get_command_argument(i, fincoords)
                  i=i+1
                  call get_command_argument(i, foutcoords)
                  i=i+1
                  call get_command_argument(i, informat)
                  if(talk) then
                    print '("Will read coordinates from file ",A)',
     &                   trim(fincoords)
                    print '("and write output to file ",A,".")',
     &                   trim(foutcoords)
                   print '("Format of input file: ",A)',
     &                   trim(informat)
                  end if
                  tol=tolfrac
                  ! read tolerance factors, if specified
                  if (i+1.le.command_argument_count()) then
                        i=i+1
                        call get_command_argument(i,transoption)
                        select case (transoption)
                        case("-tol")
                            if (i+3.le.command_argument_count()) then
                              do j=1,3
                                i=i+1
                                call get_command_argument(i,vecchar)
                                read(vecchar(1:20),*) tol(j)
                              end do
                            end if
                        case default  
                        print '(A)', 
     &                    "No trans option found."
                            i=i-1
                        end select
                  end if
                  if(talk) then
                    print'(A)', "Will use as tolerance:"
                    print'(3(F10.6))', tol
                  end if  
                  call read_coords(fincoords,informat,atoms,natoms,       &
     &                             species,nspecies,vecs)     
                  call transformsort(atoms,vecs,tol)
                  call write_coords(foutcoords,informat,atoms,            &
     &                    size(atoms),species,nspecies,vecs)
            case("mapback")
                  print'(A)', "Will map coordinates back to region betwe
     &en 0 and 1. Use only for fractional coordinates."
                  i=i+1
                  call get_command_argument(i, fincoords)
                  i=i+1
                  call get_command_argument(i, foutcoords)
                  i=i+1
                  call get_command_argument(i, informat)
                  if(talk) then
                    print '("Will read coordinates from file ",A)',
     &                   trim(fincoords)
                    print '("and write output to file ",A,".")',
     &                   trim(foutcoords)
                    print '("Format of input file: ",A)',
     &                   trim(informat)
                  end if
                  tol=tolfrac
                  ! read tolerance factors, if specified
                  if (i+1.le.command_argument_count()) then
                        i=i+1
                        call get_command_argument(i,transoption)
                        select case (transoption)
                        case("-tol")
                            if (i+3.le.command_argument_count()) then
                              do j=1,3
                                i=i+1
                                call get_command_argument(i,vecchar)
                                read(vecchar(1:20),*) tol(j)
                              end do
                            end if
                        case default  
                            print '(A)', "No trans option found."
                            i=i-1
                        end select
                  end if
                  if(talk) then
                    print'(A)', "Will use as tolerance:"
                    print'(3(F10.6))', tol
                  end if
                  call read_coords(fincoords,informat,atoms,natoms,       &
     &                             species,nspecies,vecs)     
                  call transformmapback(atoms,vecs,tol)
                  call write_coords(foutcoords,informat,atoms,            &
     &                    size(atoms),species,nspecies,vecs)
            case("mirror")
                  print'(A)', "Will mirror coordinates at plane orthogon  &
     &al to specified axis. Use only for fractional coordinates."
                  i=i+1
                  call get_command_argument(i, fincoords)
                  i=i+1
                  call get_command_argument(i, foutcoords)
                  i=i+1
                  call get_command_argument(i, informat)
                  if(talk) then
                    print '("Will read coordinates from file ",A)',
     &                   trim(fincoords)
                    print '("and write output to file ",A,".")',
     &                   trim(foutcoords)
                    print '("Format of input file: ",A)',
     &                   trim(informat)
                  end if
                  ! read axis number 
                  i=i+1
                  call get_command_argument(i,charac)
                  read(charac,*) axisint
                  select case (axisint)
                  case(1)
                  case(2)
                  case(3)
                  case default
                    call error_stop('Invalid axis')
                  end select  
                  call read_coords(fincoords,informat,atoms,natoms,       &
     &                             species,nspecies,vecs)     
                  call transformmirror(atoms,vecs,axisint)
                  call write_coords(foutcoords,informat,atoms,            &
     &                    size(atoms),species,nspecies,vecs)
            end select
         case ('-v', '--version')
!            print '(2a)', 'cofima version ', version
            print '(2a)', 'cofima beta version '
            stop
         case("--vacs")
               ! remove atoms of one species, i.e. create vacancies
               i=i+1
               call get_command_argument(i,infile)
               i=i+1
               call get_command_argument(i,informat)
               i=i+1
               call get_command_argument(i,thespecies)
               i=i+1
               call get_command_argument(i,charac)
               read(charac,*) thenumber
               i=i+1
               call get_command_argument(i,outfile)
               call rmatoms(infile,informat,thespecies,thenumber,
     &              outfile)
         case("--vasp_bs")
               ! get vasp band structure from EIGENVAL file 
               lattice='cubic' ! default
               thenumber=0 ! default number of kpoints to skip
               do while (i+1.lt.command_argument_count()) 
                i=i+1
                call get_command_argument(i,charac2)
                select case(charac2)
                case("-lattice") ! specify lattice type (hex, ...) to get special kpoints
                  i=i+1
                  call get_command_argument(i,charac)
                  read(charac,*) lattice  
                case("-nskip") ! skip the first nskip kpoints in KPOINTS
                  i=i+1
                  call get_command_argument(i,charac)
                  read(charac,*) thenumber  ! origin of sphere that contains cluster and of inner core region
                end select
               end do
               call get_vasp_bandstructure(lattice,thenumber)
         case("--vasp_bs_pr")
               ! get projected vasp band structure from PROCAR file 
               call get_vasp_projected_bandstructure()
         case("--vasp_CHG_cut_sphere")
             ! print charge density inside and outside a sphere with
             ! requested origin and radius
             !
             ! defaults:
             rcut=1.50d0 ! Angstrom
             origin=(/0.5d0,0.5d0,0.5d0/)
             infile="CHGCAR"
             do while (i+1.lt.command_argument_count()) 
              i=i+1
              call get_command_argument(i,charac2)
              select case(charac2)
                case("-origin") ! direction for averaging. 1,2,3 correspond to average perpendicular to a1, a2, a3 lattice vector
                do j=1,3
                  i=i+1
                  call get_command_argument(i,charac)
                  read(charac,*) origin(j)  ! origin of sphere that contains cluster and of inner core region
                end do
                case("-radius") ! position of cut along lattice vector cutdir in fractional coord..
                 i=i+1
                 call get_command_argument(i,charac2)
                 read(charac2,*) radius
               case("-file") ! filename (optional, default: CHGCAR)
                 i=i+1
                 call get_command_argument(i,infile)
                case default 
              end select
            end do
            call vasp_CHG_cut_sphere(origin,radius,infile)
            !
         case("--vasp_CHG_overlap")
             ! read 2 CHGCAR files (VASP output) and calculate their
             ! overlap
             i=i+1
             call get_command_argument(i,infile)
             i=i+1
             call get_command_argument(i,infile2)
            call vasp_CHG_overlap(infile,infile2)
         case("--vasp_dos_pr")
               ! get projected vasp DOS from DOSCAR and OUTCAR file 
               intnum=10   ! default number of layers if not specified
               intnum2=1 ! default direction if not specified
               t_start=0.0d0 ! default starting point if not specified
               tol=1.0E-4 ! default DOS tolerance used to detect band edges
               do while (i+1.lt.command_argument_count()) 
                 i=i+1
                 call get_command_argument(i,charac)
                 select case(charac)
                   case("-layers")
                    i=i+1
                    call get_command_argument(i,charac2)
                    read(charac2,*) intnum
                   case("-direction")
                    i=i+1
                    call get_command_argument(i,charac2)
                    read(charac2,*) intnum2
                   case("-origin")
                    i=i+1
                    call get_command_argument(i,charac2)
                    read(charac2,*) t_start
                  case("-dos_tol")
                    i=i+1
                    call get_command_argument(i,charac2)
                    read(charac2,*) tol1
                   case default 
                 end select
               end do
               call get_vasp_projected_dos(intnum,intnum2,t_start,tol1)
         case("--vasp_eps2_from_WAVEDER")
             ! read matrix elements and other properties from WAVEDER,
             ! OUTCAR, IBZKPT, and EIGENVAL and calculate imag(epsilon) 
             ! in the IPA
             nwarn=nwarn+1
             print fwarn,"only use without symmetries (ISYM=0 or spacegr  &
     &oup #1 / P1 / C_1)"          
             !
             broad=0.025d0 ! electronic smearing
             intnum=1001 ! number of frequencies
             number1=0.0d0 ! minimum frequency considered
             number2=10.0d0 ! maximum frequency considered
             charac3='fermi' ! electronic smearing type
             ! default: all valence bands (size(vbands)=1,vbands(1)=-1)
             allocate(intvec(1:1)) 
             intvec=-1
             ! default: all conduction bands (size(cbands)=1,cbands(1)=-1)
             allocate(intvec2(1:1)) 
             intvec2=-1
             ! default: all spins (size(spins)=1,spins(1)=-1)
             allocate(intvec3(1:1)) 
             intvec3=-1
             allocate(intvec4(1:1)) 
             ! default: all kpoints (size(kpoints)=1,kpoints(1)=-1)
             intvec4=-1
             do while (i+1.lt.command_argument_count()) 
               i=i+1
               call get_command_argument(i,charac)
               select case(charac)
                 case("-broad")
                  i=i+1
                  call get_command_argument(i,charac2)
                  read(charac2,*) broad
                 case("-nomega")
                  i=i+1
                  call get_command_argument(i,charac2)
                  read(charac2,*) intnum
                 case("-omegamin")
                  i=i+1
                  call get_command_argument(i,charac2)
                  read(charac2,*) number1
                 case("-omegamax")
                  i=i+1
                  call get_command_argument(i,charac2)
                  read(charac2,*) number2
                 case("-vbands")
                  i=i+1
                  call get_command_argument(i,charac2)
                  read(charac2,*) intnum4 ! number of valence bands
                  deallocate(intvec)
                  allocate(intvec(1:intnum4))
                  do j=1,intnum4
                    i=i+1
                    call get_command_argument(i,charac2)
                    read(charac2,*) intvec(j) ! valence band number
                  end do ! j
                 case("-vbandrange")
                  i=i+1
                  call get_command_argument(i,charac2)
                  read(charac2,*) intnum4 ! first valence band to consider
                  i=i+1
                  call get_command_argument(i,charac2)
                  read(charac2,*) intnum5 ! last valence band to consider
                  deallocate(intvec)
                  allocate(intvec(1:intnum5-intnum4+1))
                  do j=1,size(intvec)
                    intvec(j)=intnum4+j-1 ! valence band number
                  end do ! j
                 case("-cbands")
                  i=i+1
                  call get_command_argument(i,charac2)
                  read(charac2,*) intnum4 ! number of conduction bands
                  deallocate(intvec2)
                  allocate(intvec2(1:intnum4))
                  do j=1,intnum4
                    i=i+1
                    call get_command_argument(i,charac2)
                    read(charac2,*) intvec2(j) ! conduction band number
                  end do ! j
                 case("-cbandrange")
                  i=i+1
                  call get_command_argument(i,charac2)
                  read(charac2,*) intnum4 ! first conduction band to consider
                  i=i+1
                  call get_command_argument(i,charac2)
                  read(charac2,*) intnum5 ! last conduction band to consider
                  deallocate(intvec2)
                  allocate(intvec2(1:intnum5-intnum4+1))
                  do j=1,size(intvec2)
                    intvec2(j)=intnum4+j-1 ! valence band number
                  end do ! j
                 case("-spins")
                  i=i+1
                  call get_command_argument(i,charac2)
                  read(charac2,*) intnum4 ! number of spins (1 for nonspinpolarized, 2 for spin-polarized)
                  deallocate(intvec3)
                  allocate(intvec3(1:intnum4))
                  do j=1,intnum4
                    i=i+1
                    call get_command_argument(i,charac2)
                    read(charac2,*) intvec3(j) ! spin number
                  end do ! j
                 case("-kpoints")
                  i=i+1
                  call get_command_argument(i,charac2)
                  read(charac2,*) intnum4 ! number of kpoints to be considered
                  deallocate(intvec4)
                  allocate(intvec4(1:intnum4))
                  do j=1,intnum4
                    i=i+1
                    call get_command_argument(i,charac2)
                    read(charac2,*) intvec4(j) ! kpoint number
                  end do ! j
                 case("-smearing")
                  i=i+1
                  call get_command_argument(i,charac3)
               end select
             end do
             call vasp_eps2_from_WAVEDER(intnum,number1,number2,broad,    &
     &            intvec,intvec2,intvec3,intvec4,charac3) ! nomega,omegamin,omegamax, broadening, 
                  !considered valence bands,considered conduction bands,
                  !considered spins, considered kpoints,smearing
         case("--vasp_get_eps_vs_omega")
             ! get frequency-dependent dielectric
             ! matrix from vasp output file. Basically a grep command.
             ! Optionally rotate epsilon to another coordinate system.
             ! In that case supply a rotation matrix
             if (i+1.lt.command_argument_count()) then
               i=i+1
               call get_command_argument(i,charac)
               select case(charac)
                 case('-rotate')
                 if (i+9.le.command_argument_count()) then
                   if (allocated(allomat1)) deallocate(allomat1)
                   allocate(allomat1(1:3,1:3))
                   do j=1,3
                     do k=1,3
                       i=i+1
                       call get_command_argument(i,charac)
                       read(charac,*)  allomat1(j,k)
                     end do
                   end do
                   call vasp_get_eps_vs_omega(allomat1)
                 end if 
               end select
             else
               call vasp_get_eps_vs_omega()
             end if
         case("--vasp_get_eps_vs_omega_xml")
             ! get frequency-dependent dielectric
             ! matrix from a TDDFT calculation from vasp_run.xml. Basically a grep command.
             ! Optionally rotate epsilon to another coordinate system.
             ! In that case supply a rotation matrix
             if (i+1.lt.command_argument_count()) then
               i=i+1
               call get_command_argument(i,charac)
               select case(charac)
                 case('-rotate')
                 if (i+9.le.command_argument_count()) then
                   if (allocated(allomat1)) deallocate(allomat1)
                   allocate(allomat1(1:3,1:3))
                   do j=1,3
                     do k=1,3
                       i=i+1
                       call get_command_argument(i,charac)
                       read(charac,*)  allomat1(j,k)
                     end do
                   end do
                   call vasp_get_eps_vs_omega_xml(allomat1)
                 end if 
               end select
             else
               call vasp_get_eps_vs_omega_xml()
             end if
         case("--vasp_get_pol")
             ! get ferroelectric polarization from a Berry phase
             ! calculation in vasp (LCALCPOL=.TRUE.) 
             call vasp_get_pol()
         case("--vasp_plot_CHG")
             ! get CHG from VASP CHGCAR file and write to gnuplottable
             ! file CHG4PLOTTING
             ! defaults:
             intnum=1
             intnum2=1
             rcut=0.0
             infile="CHGCAR"
             do while (i+1.lt.command_argument_count()) 
              i=i+1
              call get_command_argument(i,charac2)
              select case(charac2)
                case("-avdir") ! direction for averaging. 1,2,3 correspond to average perpendicular to a1, a2, a3 lattice vector
                 i=i+1
                 call get_command_argument(i,charac2)
                 read(charac2,*) intnum
                case("-cutdir") ! direction for cut. 1,2,3 correspond to cut perpendicular to a1, a2, a3 lattice vector
                 i=i+1
                 call get_command_argument(i,charac2)
                 read(charac2,*) intnum2
                case("-cutpos") ! position of cut along lattice vector cutdir in fractional coord..
                 i=i+1
                 call get_command_argument(i,charac2)
                 read(charac2,*) rcut
               case("-file") ! filename (optional, default: CHGCAR)
                 i=i+1
                 call get_command_argument(i,infile)
                case default 
              end select
            end do
            call vasp_plot_CHG(intnum,intnum2,rcut,infile)
         case('--vasp_dE_dk')   
            ! prints a file VASP_E_vs_K_n with the energy eigenvalues vs
            ! absolute kpoint coordinates. Needs as input parameter the
            ! band index n, and the files EIGENVAL, OUTCAR need to be
            ! present 
            ! Another file is written with the energy on a k plane.
            ! Defaults: plane at k1=0. 
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) intnum ! band index
            ! defaults:
            intnum2=1 ! direction for cutting (1: cut at constant fractional k x1, etc.)
            rcut=0.0d0 ! cut at x1=rcut, or x2=rcut, ...
            do while (i+1.lt.command_argument_count()) 
             i=i+1
             call get_command_argument(i,charac)
             select case(charac)
               case("-cutpos") ! direction for averaging. 1,2,3 correspond to average perpendicular to a1, a2, a3 lattice vector
                i=i+1
                call get_command_argument(i,charac2)
                read(charac2,*) rcut
               case("-cutdir") ! direction for averaging. 1,2,3 correspond to average perpendicular to a1, a2, a3 lattice vector
                i=i+1
                call get_command_argument(i,charac2)
                read(charac2,*) intnum2
              case default
              end select
            end do ! while 
            call vasp_dE_dk(intnum,intnum2,rcut)
         case('--vasp_kgrid')   
            ! prints a file KPOINTS
            ! with a Gamma-centered Monkhorst-Pack k grid
            allocate(intvec(3))
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) intvec(1) ! nk1
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) intvec(2) ! nk1
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) intvec(3) ! nk1
            ! defaults:
            shift=0.0d0 ! cut at x1=rcut, or x2=rcut, ...
            do while (i+1.lt.command_argument_count()) 
             i=i+1
             call get_command_argument(i,charac)
             select case(charac)
               case("-shift") ! direction for averaging. 1,2,3 correspond to average perpendicular to a1, a2, a3 lattice vector
                i=i+1
                call get_command_argument(i,charac2)
                read(charac2,*) shift(1)
                i=i+1
                call get_command_argument(i,charac2)
                read(charac2,*) shift(2)
                i=i+1
                call get_command_argument(i,charac2)
                read(charac2,*) shift(3)
              case default
              end select
            end do ! while 
            print '(8x,"creating ",I0,"x",I0,"x",I0," k-point grid")',    &
     &        intvec(1:3)         
            call vasp_kgrid(intvec,shift)
            deallocate(intvec)
         case('--vasp_kpath')   
            ! prints a file KPOINTS.PATH with a k-path
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) intnum ! number of points in the path
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) origin(1) ! starting k-point 1. coord. (fractional)
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) origin(2) ! starting k-point 2. coord. (fractional)
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) origin(3) ! starting k-point 3. coord. (fractional)
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) vector(1) ! end k-point 1. coord. (fractional)
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) vector(2) ! end k-point 2. coord. (fractional)
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) vector(3) ! end k-point 3. coord. (fractional)
            call vasp_kpath(intnum,origin,vector)
         case('--vasp_kaux')   
            ! prints a file KPOINTS.AUX with 6 auxiliary k-points for
            ! each kpoint in provided KPOINT file.
            i=i+1
            call get_command_argument(i,infile)
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) tol(1) ! end k-point 1. coord. (fractional)
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) tol(2) ! end k-point 2. coord. (fractional)
            i=i+1
            call get_command_argument(i,charac)
            read(charac,*) tol(3) ! end k-point 3. coord. (fractional)
            call vasp_kaux(infile,tol)
         case('--vasp_kmerge')   
            ! reads two KPOINT files and merges them to KPOINT.MERGED
            i=i+1
            call get_command_argument(i,infile)
            i=i+1
            call get_command_argument(i,infile2)
            i=i+1
            call vasp_kmerge(infile,infile2)
         case("--vasp_read_WAVEDER")
            ! reads binary WAVEDER file and prints to formatted file
            ! WAVEDER.dat
            call vasp_read_WAVEDER()
         case("--vasp_write_BORN")
            ! reads dielectric tensor and Born effective charges from OUTCAR
            ! file and writes them to file BORN (used by phonopy).
            !
            ! begin optional : only BEC of irreducible atoms
            !
            natoms2=-1
            do while (i+1.lt.command_argument_count()) 
             i=i+1
             call get_command_argument(i,charac)
             select case(charac)
               case("-natoms") ! direction for averaging. 1,2,3 correspond to average perpendicular to a1, a2, a3 lattice vector
                i=i+1
                call get_command_argument(i,charac2)
                read(charac2,*) natoms2
                allocate(intvec(natoms2))
               case("-atoms") ! list indices of atoms for which BEC should be written to BORN
                do j=1,natoms2
                  i=i+1
                  call get_command_argument(i,charac2)
                  read(charac2,*) intvec(j)
                end do ! j
              case default
              end select
            end do ! while 
            !
            ! end optional : only BEC of irreducible atoms
            !
            call vasp_write_BORN(natoms2,intvec)
         case("--WAVECAR_ana")
            ! reads binary WAVECAR file header 
            call WAVECAR_ana()
         case ('--xcoulomb','--XCOULOMB')
               ! calculates the classical Coulomb energy of an exciton
               ! from e and h wave functions in cube files 
            i=i+1
            call get_command_argument(i,infile)
            i=i+1
            call get_command_argument(i,infile2)
            i=i+1
            call get_command_argument(i,informat)
            ! read optional arguments
            step=4
            step2=2
            !rcut=0.650d0
            rcut=1.10d0
            do while (i+1.lt.command_argument_count()) 
              i=i+1
              call get_command_argument(i,charac)
              select case(charac)
                case("-step")
                 i=i+1
                 call get_command_argument(i,charac2)
                 read(charac2,*) step
                case("-finestep")
                 i=i+1
                 call get_command_argument(i,charac2)
                 read(charac2,*) step2
                case("-rmin")
                 i=i+1
                 call get_command_argument(i,charac2)
                 read(charac2,*) rcut
                case default 
              end select
            end do
            call get_xcoulomb(infile,infile2,informat,step,step2,         &
     &                        rcut)
         case default
            print '(a,a,/)', 'Unrecognized command-line option: ', arg
            call print_help(morehelpon)
            stop
         end select
         i=i+1
      end do

      !if "command-line only" mode, stop here.
      if (cmdlonly) goto 4000

      ! read commands from input file
      open(10,file=inputfile,status='old', err=3000)
      open(12,file=outputfile,status='replace')
      if(talk) then
      write(12,*)'***************************************************'
      write(12,'(" * This is cofima, v ",A,
     & ".                          *")') version
      write(12,*)'*                                                 *'
      write(12,*)'* New capabilities/changes:                       *'
      write(12,*)'* -new option "sphere" for the cut command        *'
      write(12,*)'* -new formats supported: xyz (r/w), abinit (ro), *'
      write(12,*)'*  nwchem (ro)                                    *'
      write(12,*)'***************************************************'
      end if

10    read(10,'(A100)',end=1000) line  
      if (len(line).ge.1) line=adjustl(line)
      if (line(1:1).eq.'#') goto 10

      ! check if subroutine sort_molecules should be performed
      if (index(line,'sort molecules').ge.1) then
            write(12,*)' '
            write(12,*)' Calling subroutine sort_molecules '
            write(12,*)' '
            call sort_molecules()
      end if
      
!      ! check if subroutine abs2frac should be performed
!      if (index(line,'absolute to fractional').ge.1) then
!            write(12,*)' '
!            write(12,*)' Calling subroutine abs2frac '
!            write(12,*)' '
!            call abs2frac()
!      end if
      
!      ! check if subroutine frac2abs should be performed
!      if (index(line,'fractional to absolute').ge.1) then
!            write(12,*)' '
!            write(12,*)' Calling subroutine frac2abs '
!            write(12,*)' '
!            call frac2abs()
!      end if
      
      ! check if subroutine res2config should be performed
      if(index(line,'res to config').ge.1) then
            write(12,*)' '
            write(12,*)' Calling subroutine res2config '
            write(12,*)' '
            call res2config()
      end if
      
      ! check if subroutine coorat2xsf should be performed
      if(index(line,'coorat to xsf').ge.1) then
            write(12,*)' '
            write(12,*)' Calling subroutine coorat2xsf '
            write(12,*)' '
            call coorat2xsf()
      end if
      
      ! check if subroutine trg2data should be performed
      if(index(line,'trg to data').ge.1) then
            write(12,*)' '
            write(12,*)' Calling subroutine trg2data '
            write(12,*)' '
            call trg2data()
      end if
      
      ! check if subroutine gulp_md_data should be performed
      line=adjustl(line)
      if(index(line,'gulp').ge.1.and.index(line,'md').ge.1.and.
     &   index(line,'data').ge.1) then
            write(12,*)' '
            write(12,*)' Calling subroutine gulp_md_data '
            write(12,*)' '
            call gulp_md_data()
      end if
      
      ! check if subroutine crysanpero should be performed
      line=adjustl(line)
      if(index(line,'crysanpero').ge.1) then
            write(12,*)' '
            write(12,*)' Calling subroutine scrysanpero '
            write(12,*)' '
            call scrysanpero()
            write(12,'(//1x,"reading input file...")')
      end if
      
      ! check if subroutine MBPP_transform_coords should be performed
      line=adjustl(line)
      if(index(line,'mbpp').ge.1.and.index(line,'transform').ge.1
     &  .and.index(line,'coo').ge.1) then
            write(12,*)' '
            write(12,*)' Calling subroutine MBPP_transform_coords '
            write(12,*)' '
            call MBPP_transform_coords()
      end if
      
      ! check if subroutine config2coorat should be performed
      line=adjustl(line)
      if(index(line,'config to coorat').ge.1) then
            write(12,*)' '
            write(12,*)' Calling subroutine config2coorat. '
            write(12,*)' '
            call config2coorat()
      end if
      
      ! check if subroutine config2coorat should be performed
      line=adjustl(line)
      if(index(line,'trg to coorat').ge.1) then
            write(12,*)' '
            write(12,*)' Calling subroutine trg2coorat. '
            write(12,*)' '
            call trg2coorat()
      end if
      
      ! check if subroutine sxyz2cfg should be performed
      line=adjustl(line)
      if(index(line,'sxyz to cfg').ge.1) then
            write(12,*)' '
            write(12,*)' Calling subroutine sxyz2cfg. '
            write(12,*)' '
            call sxyz2cfg()
      end if
      
      ! check if subroutine lammps_xyz2sxyz should be performed
      line=adjustl(line)
      if(index(line,'lammps xyz to sxyz').ge.1) then
            write(12,*)' '
            write(12,*)' Calling subroutine lammps_xyz2sxyz. '
            write(12,*)' '
            call lammps_xyz2sxyz()
      end if
      
      ! check if subroutine prop_corr should be performed
      line=adjustl(line)
      if(index(line,'prop').ge.1.and.index(line,'corr').ge.1) then
            write(12,*)' '
            write(12,*)' Calling subroutine prop_corr. '
            write(12,*)' '
            call prop_corr()
      end if
      
      ! check if subroutine sxyz2config should be performed
      line=adjustl(line)
      if(index(line,'sxyz to config').ge.1) then
            write(12,*)' '
            write(12,*)' Calling subroutine sxyz2config. '
            write(12,*)' '
            call sxyz2config()
      end if
      
      ! check if subroutine res2coorat should be performed
      line=adjustl(line)
      if(index(line,'res to coorat').ge.1) then
            write(12,*)' '
            write(12,*)' Calling subroutine res2coorat. '
            write(12,*)' '
            call res2coorat()
      end if
      
      ! check if subroutine coorat2gin should be performed
      line=adjustl(line)
      if(index(line,'coorat to gin').ge.1) then
            write(12,*)' '
            write(12,*)' Calling subroutine coorat2gin. '
            write(12,*)' '
            call coorat2gin()
      end if
      
      ! check if subroutine cfg_del_species should be performed
      line=adjustl(line)
      if(index(line,'cfg del').ge.1.and.
     & index(line,'spec').ge.1) then
            write(12,*)' '
            write(12,*)' Calling subroutine cfg_del_species. '
            write(12,*)' '
            call cfg_del_species()
      end if
      
      ! check if subroutine data2cfg should be performed
      line=adjustl(line)
      if(index(line,'lammps data to cfg').ge.1) then
            write(12,*)' '
            write(12,*)' Calling subroutine data2cfg. '
            write(12,*)' '
            call data2cfg()
      end if
      
      ! check if subroutine refcoords should be performed
      line=adjustl(line)
      if(index(line,'refcoords').ge.1) then
            write(12,*)' '
            write(12,*)' Calling subroutine refcoords. '
            write(12,*)' '
            call refcoords()
      end if
      
      goto 10
      
1000  continue      
      inquire(unit=12,opened=isopen12)
      close(10)
      call cpu_time(t_end)
      !write(12,*) ' '
      !write(12,'(2x,"Errors:",I)') nerr
      !write(12,1111) t_end-t_start
!      write(12,'("  Comments             :",I)')ncomm
      if(isopen12) then
        WRITE(FMT2,*) ncomm
        WRITE(FMT1,*) '(8x,"     Comments             :",I',
     &   len(trim(FMT2)),')'
        write(12,FMT1)  ncomm 
!      write(12,'("  Warnings             :",I)')nwarn
        WRITE(FMT2,*) nwarn
        WRITE(FMT1,*) '(8x,"     Warnings             :",I',
     &   len(trim(FMT2)),')'
        write(12,FMT1)  nwarn 
!      write(12,'("  Errors               :",I)')nerr
        WRITE(FMT2,*) nerr
        WRITE(FMT1,*) '(8x,"     Errors               :",I',
     &   len(trim(FMT2)),')'
        write(12,FMT1)  nerr 
        write(12,'("  Computation time in s:",F15.6)')t_end-t_start  
! 1111 format('  Computation time in s:',F15.6) 
        close(12)
      else
        WRITE(FMT2,*) ncomm
        WRITE(FMT1,*) '(8x,"     Comments             :",I',
     &   len(trim(FMT2)),')'
        print FMT1,  ncomm 
!      write(12,'("  Warnings             :",I)')nwarn
        WRITE(FMT2,*) nwarn
        WRITE(FMT1,*) '(8x,"     Warnings             :",I',
     &   len(trim(FMT2)),')'
        print FMT1, nwarn 
!      write(12,'("  Errors               :",I)')nerr
        WRITE(FMT2,*) nerr
        WRITE(FMT1,*) '(8x,"     Errors               :",I',
     &   len(trim(FMT2)),')'
        print FMT1, nerr 
        print '("  Computation time in s:",F15.6)',t_end-t_start  
      end if 
      stop

3000  continue    
      !nerr=nerr+1
      !print ferrmssg,trim(adjustl('Input file not found. Quitting'))
!      ncomm=ncomm+1
!      print fcomm,trim(adjustl(
!     & 'I will only execute the commands passed on the command line.'))
4000  call cpu_time(t_end)
      inquire(unit=12,opened=isopen12)
      if(talk) then 
        WRITE(FMT2,*) ncomm
        WRITE(FMT1,*) '(8x,"     Comments             :",I',
     &   len(trim(FMT2)),')'
        print FMT1, ncomm 
!        write(12,'("  Warnings             :",I)')nwarn
        WRITE(FMT2,*) nwarn
        WRITE(FMT1,*) '(8x,"     Warnings             :",I',
     &   len(trim(FMT2)),')'
        print FMT1, nwarn 
!        write(12,'("  Errors               :",I)')nerr
        WRITE(FMT2,*) nerr
        WRITE(FMT1,*) '(8x,"     Errors               :",I',
     &   len(trim(FMT2)),')'
        print FMT1, nerr 
        print '("  Computation time in s:",F15.6)',t_end-t_start  
!        print '("  Comments             :",I)', ncomm
!        print '("  Warnings             :",I)', nwarn
!        print '("  Errors               :",I)', nerr
!        print '("  Computation time in s:",F15.6)',t_end-t_start  
      end if
      if(isopen12) close(12)
      !print 1111, t_end-t_start
      stop
    

      end program

