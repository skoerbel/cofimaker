      module postproc
      implicit none

      contains
c---------------------------------------------------------------------

      subroutine locpotav(fincoords,finformat,direction,avlength)
      ! averages VASP LOCPOT file in a chosen direction (1,2,3) 
      use defs
      use misc 
      use readcoords
      implicit none
      character(len=*), intent(in) :: fincoords,finformat
      integer, intent(in) :: direction
      double precision, intent(in) :: avlength
      ! local variables:
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms, nspecies,ngrid(1:3)
      integer i,j,jstart,jend,navmacro,navmacro_half,jtrue
      double precision vecs(1:3,1:3)!,vol
      double precision, allocatable :: potarr(:,:,:),potav(:)
      double precision, allocatable :: potavmacro(:)
      character(len=256) :: line
!      
      if(talk) print fsubstart,"locpotav" 
      if (talk) print '(8x,"averageing file: ",A40)',fincoords
      if (direction.ne.1.and.direction.ne.2.and.direction.ne.3) then
        goto 100
      end if
!
      ! read coordinate and cell information     
      call read_coords(fincoords,finformat,atoms,natoms,species,
     &                       nspecies,vecs)
      open(51,file=fincoords,status='old',err=101)
      ! skip "header"
      do i=1,natoms+9
        read(51,'(A256)',err=101,end=101) line
      end do 
      ! read grid dimensions
      read(51,*) ngrid(1:3)
      print*, "grid dimensions:",ngrid(1:3)  
      allocate(potarr(ngrid(1),ngrid(2),ngrid(3)))
      read(51,*) potarr(1:ngrid(1),1:ngrid(2),1:ngrid(3))
      close(51)
      !
      ! begin average
      allocate(potav(ngrid(direction))) 
      do i=1,ngrid(direction)
        potav(i)=0.0d0
        if(direction==1) then
          do j=1,ngrid(2)
            potav(i)=potav(i)+sum(potarr(i,j,:))
          end do
        end if
        if(direction==2) then
          do j=1,ngrid(1)
            potav(i)=potav(i)+sum(potarr(j,i,:))
          end do
        end if
        if(direction==3) then
          do j=1,ngrid(1)
            potav(i)=potav(i)+sum(potarr(j,:,i))
          end do
        end if
      end do
      ! end average
      ! begin write average
      !call vecs2vol(vecs,vol)
      potav=potav*dble(ngrid(direction))                                  &
     &     /(dble(ngrid(1)*ngrid(2)*ngrid(3)))
      open(51,file="LOCPOTAV",status='replace')
      do i=1,ngrid(direction)
        write(51,*) dble(i)/dble(ngrid(direction)),potav(i)
      end do
      close(51)
      ! end write average
      deallocate(potarr)
      !
      ! begin calculate and write macroscopic (sliding window) average
      allocate(potavmacro(ngrid(direction))) 
      navmacro=ceiling(ngrid(direction)*avlength)
      navmacro_half=int(navmacro/2)
      do i=1,ngrid(direction)
        !print*, "i=",i
        potavmacro(i)=0.0d0
        jstart=i-navmacro_half 
        jend=jstart+navmacro-1  
        !print*, "jstart,jend=",jstart,jend
        do j=jstart,jend
          jtrue=j
          do while (jtrue.gt.ngrid(direction)) 
            jtrue=jtrue-ngrid(direction)
          end do 
          do while (jtrue.lt.1) 
            jtrue=jtrue+ngrid(direction)
          end do 
          potavmacro(i)=potavmacro(i)+potav(jtrue)
        end do
      end do 
      potavmacro=potavmacro/dble(navmacro)                             
      open(51,file="LOCPOTAVMACRO",status='replace')
      do i=1,ngrid(direction)
        write(51,*) dble(i)/dble(ngrid(direction)),potavmacro(i)
      end do
      close(51)
      ! end calculate and write macroscopic average
      deallocate(potav)
      deallocate(potavmacro)
      !
      ! end normally:
      if(talk) print fsubendext, "locpotav"
      return
      !
      ! Errors:
 100  nerr=nerr+1
      print ferrmssg," direction has to be 1,2, or 3 (cofima --help)"
      close(51)
      return
 101  nerr=nerr+1
      print ferrmssg," when opening/reading LOCPOT file"
      close(51)
      return
      end subroutine locpotav  
      
c---------------------------------------------------------------------

      subroutine get_com(file1,format1)
      use defs
      use misc
      use readcoords
      implicit none
      ! calculates the Center Of Molecule
      !
      ! variables
      character(len=*), intent(in) :: file1,format1
      ! internal variables
      type(atom),allocatable :: atoms(:)
      type(element),allocatable :: species(:)
      double precision vecs(1:3,1:3),com(1:3)
      integer natoms,nspecies
      integer iatom 
      !
      if(talk) print fsubstart,"get_com" 
      if (talk) print '(8x,"file: ",A40)',file1
      if (talk) print '(8x,"format: ",A40)',format1
      !
      ! begin read coordinates from file
      call read_coords(file1,format1,atoms,natoms,species,
     &                       nspecies,vecs)
      !print*,vecs
      ! end read coordinates from file
      !
      ! begin calculate com
      com=0.0d0
      do iatom=1,natoms
        !print*,atoms(iatom)%abswhere
        !call frac2abs(atoms(iatom)%where,vecs,atoms(iatom)%abswhere)
        com=com+atoms(iatom)%abswhere
      end do 
      com=com/dble(natoms)
      print '(8x,"COM in Angs (probably):",3(F10.6))',com(1:3)
      ! end calculate com
      !
      ! end normally:
      deallocate(atoms,species)
      if(talk) print fsubendext, "get_com"
      return
      
      end subroutine get_com

c---------------------------------------------------------------------

      subroutine get_comass(file1,format1)
      use defs
      use misc
      use readcoords
      implicit none
      ! calculates the Center Of Molecule
      !
      ! variables
      character(len=*), intent(in) :: file1,format1
      ! internal variables
      type(atom),allocatable :: atoms(:)
      type(element),allocatable :: species(:)
      double precision vecs(1:3,1:3),comass(1:3)
      integer natoms,nspecies
      integer iatom 
      !
      if(talk) print fsubstart,"get_comass" 
      if (talk) print '(8x,"file: ",A40)',file1
      if (talk) print '(8x,"format: ",A40)',format1
      !
      ! begin read coordinates from file
      call read_coords(file1,format1,atoms,natoms,species,
     &                       nspecies,vecs)
      !print*,vecs
      ! end read coordinates from file
      !
      ! begin calculate comass
      comass=0.0d0
      do iatom=1,natoms
        !print*,atoms(iatom)%abswhere
        !call frac2abs(atoms(iatom)%where,vecs,atoms(iatom)%abswhere)
        comass=comass+atoms(iatom)%abswhere*atoms(iatom)%mass
      end do 
      comass=comass/sum(atoms(:)%mass)
      print '(8x,"Center of mass in Angs:",3(F10.6))',comass(1:3)
      ! end calculate comass
      !
      ! end normally:
      deallocate(atoms,species)
      if(talk) print fsubendext, "get_comass"
      return
      
      end subroutine get_comass

c---------------------------------------------------------------------

      subroutine relax_atoms(atoms,vecs,step)
      use defs  
      use misc, only : abs2frac
      implicit none
      double precision, intent(in) :: step,vecs(1:3,1:3) 
      type(atom), intent(inout) :: atoms(:)
      integer i,natoms

      if(talk) print fsubstart,"relax_atoms" 
      if(talk) print '(8x,"forces are multiplied by ",F10.6)',step 
      natoms=size(atoms)
      if(talk) print '(8x,I6," atoms")',natoms

      do i=1,natoms
        !print*,atoms(i)%force(1:3)
        atoms(i)%abswhere=atoms(i)%abswhere+atoms(i)%force*step
        call abs2frac(atoms(i)%abswhere,vecs,atoms(i)%where)
      end do

      if(talk) print fsubendext, "relax_atoms"
      return

      end subroutine relax_atoms

c---------------------------------------------------------------------

      subroutine relax_rigid_molecule(atoms,vecs,step)
      use defs  
      use misc, only : abs2frac,cross_product
      use transform, only : rotatecoords
      implicit none
      double precision, intent(in) :: step,vecs(1:3,1:3) 
      double precision av_radial_force, delta_r(1:3), radial_force
      double precision sum_of_forces(1:3),sum_of_torques(1:3)
      double precision angle,vector(1:3),com(1:3)
      type(atom), intent(inout) :: atoms(:)
      integer i,natoms

      if(talk) print fsubstart,"relax_rigid_molecule" 

      natoms=size(atoms)
      sum_of_forces=0.0d0

      do i=1,natoms
        sum_of_forces=sum_of_forces+atoms(i)%force
      end do ! i

      ! rigid tranlation
      print '(8x,"translating molecule by ",3(F12.6,x)," Angs ")',        &
     &  sum_of_forces*step/dble(natoms) 
      do i=1,natoms
        atoms(i)%abswhere=atoms(i)%abswhere+sum_of_forces*step            &
     &    /dble(natoms)
        call abs2frac(atoms(i)%abswhere,vecs,atoms(i)%where)
      end do
      !  
      ! move center of molecule (com) to 0
      com=0.0d0
      do i=1,natoms
        com=com+atoms(i)%abswhere
      end do 
      com=com/dble(natoms)
      do i=1,natoms
        atoms(i)%abswhere=atoms(i)%abswhere-com
        call abs2frac(atoms(i)%abswhere,vecs,atoms(i)%where)
      end do 
      !
      ! calculate total torque
      sum_of_torques=0.0d0
      do i=1,natoms
        sum_of_torques=sum_of_torques+cross_product(atoms(i)%abswhere,    &
     &    atoms(i)%force)
      end do 
      !
      ! rigid rotation
      !angle=norm2(sum_of_torques)*step/dble(natoms)
      angle=norm2(sum_of_torques)*step
      if (angle.gt.0.0d0) then
        vector=sum_of_torques/norm2(sum_of_torques) ! not really needed, just for neat output 
      else
        vector=0.0d0
        vector(1)=1.0d0
      end if
      print '(8x,"rotating molecule by ",F9.3," degrees around ",         &
     &   3(F9.3,x))', angle, vector
      call rotatecoords(atoms,vecs,0,angle,vector)
      !
      ! move com back
      do i=1,natoms
        atoms(i)%abswhere=atoms(i)%abswhere+com
        call abs2frac(atoms(i)%abswhere,vecs,atoms(i)%where)
      end do 

      ! expand or shrink molecule
      av_radial_force=0.0d0

      ! calculate average radial force
      do i=1,natoms
        if (norm2(atoms(i)%abswhere-com).gt.0.0d0) then
          radial_force=dot_product(atoms(i)%force,atoms(i)%abswhere-com)  &
     &               /(norm2(atoms(i)%abswhere-com))**2
        else 
          radial_force=0.0d0
        end if
        av_radial_force=av_radial_force + radial_force
      end do ! i
      av_radial_force=av_radial_force/dble(natoms)

      ! rigid breathing
      print '(8x," expanding molecule by ",3(F12.6,x)," Angs ")',         &
     &  av_radial_force*step
      do i=1,natoms
        delta_r=(atoms(i)%abswhere-com)*av_radial_force*step           
        atoms(i)%abswhere=atoms(i)%abswhere+delta_r
        call abs2frac(atoms(i)%abswhere,vecs,atoms(i)%where)
      end do

      if(talk) print fsubendext, "relax_rigid_molecule"
      return

      end subroutine relax_rigid_molecule

c---------------------------------------------------------------------

      subroutine relax_rigid_molecule_rot(atoms,vecs,step)
      use defs  
      use misc, only : abs2frac,cross_product
      use transform, only : rotatecoords
      implicit none
      double precision, intent(in) :: step,vecs(1:3,1:3) 
      double precision sum_of_torques(1:3)
      double precision angle,vector(1:3),com(1:3)
      type(atom), intent(inout) :: atoms(:)
      integer i,natoms

      if(talk) print fsubstart,"relax_rigid_molecule_rot" 

      natoms=size(atoms)
      !  
      ! move center of molecule (com) to 0
      com=0.0d0
      do i=1,natoms
        com=com+atoms(i)%abswhere
      end do 
      com=com/dble(natoms)
      do i=1,natoms
        atoms(i)%abswhere=atoms(i)%abswhere-com
        call abs2frac(atoms(i)%abswhere,vecs,atoms(i)%where)
      end do 
      !
      ! calculate total torque
      sum_of_torques=0.0d0
      do i=1,natoms
        sum_of_torques=sum_of_torques+cross_product(atoms(i)%abswhere,    &
     &    atoms(i)%force)
      end do 
      !
      ! rigid rotation
      !angle=norm2(sum_of_torques)*step/dble(natoms)
      angle=norm2(sum_of_torques)*step
      if (angle.gt.0.0d0) then
        vector=sum_of_torques/norm2(sum_of_torques) ! not really needed, just for neat output 
      else
        vector=0.0d0
        vector(1)=1.0d0
      end if
      print '(8x,"rotating molecule by ",F9.3," degrees around ",         &
     &   3(F9.3,x))', angle, vector
      call rotatecoords(atoms,vecs,0,angle,vector)
      !
      ! move com back
      do i=1,natoms
        atoms(i)%abswhere=atoms(i)%abswhere+com
        call abs2frac(atoms(i)%abswhere,vecs,atoms(i)%where)
      end do 

      if(talk) print fsubendext, "relax_rigid_molecule_rot"
      return

      end subroutine relax_rigid_molecule_rot

c---------------------------------------------------------------------

      subroutine relax_rigid_molecule_rot2(atoms,vecs,step)
      ! same as relax_rigid_molecule_rot, but with a different way to
      ! calculate the total rotation angle and axis
      use defs  
      use misc, only : abs2frac,cross_product
      use transform, only : rotatecoords
      implicit none
      double precision, intent(in) :: step,vecs(1:3,1:3) 
      double precision sum_of_torques(1:3)
      double precision angle,vector(1:3),com(1:3)
      type(atom), intent(inout) :: atoms(:)
      integer i,natoms

      if(talk) print fsubstart,"relax_rigid_molecule_rot2" 

      natoms=size(atoms)
      !  
      ! move center of molecule (com) to 0
      com=0.0d0
      do i=1,natoms
        com=com+atoms(i)%abswhere
      end do 
      com=com/dble(natoms)
      do i=1,natoms
        atoms(i)%abswhere=atoms(i)%abswhere-com
        call abs2frac(atoms(i)%abswhere,vecs,atoms(i)%where)
      end do 
      !
      ! calculate average rotation angle and axis
      sum_of_torques=0.0d0
      do i=1,natoms
        sum_of_torques=sum_of_torques+cross_product(atoms(i)%abswhere,    &
     &    atoms(i)%force)/norm2(atoms(i)%abswhere)**2
      end do 
      sum_of_torques=sum_of_torques/dble(natoms)
      !
      ! rigid rotation
      angle=norm2(sum_of_torques)*step*360.0d0/(2.0d0*Pi)
      if (angle.gt.0.0d0) then
        vector=sum_of_torques/norm2(sum_of_torques) ! not really needed, just for neat output 
      else
        vector=0.0d0
        vector(1)=1.0d0
      end if
      print '(8x,"rotating molecule by ",F9.3," degrees around ",         &
     &   3(F9.3,x))', angle, vector
      call rotatecoords(atoms,vecs,0,angle,vector)
      !
      ! move com back
      do i=1,natoms
        atoms(i)%abswhere=atoms(i)%abswhere+com
        call abs2frac(atoms(i)%abswhere,vecs,atoms(i)%where)
      end do 

      if(talk) print fsubendext, "relax_rigid_molecule_rot2"
      return

      end subroutine relax_rigid_molecule_rot2

c---------------------------------------------------------------------

      subroutine relax_rigid_molecule_trala(atoms,vecs,step)
      use defs  
      use misc, only : abs2frac
      implicit none
      double precision, intent(in) :: step,vecs(1:3,1:3) 
      double precision sum_of_forces(1:3)
      double precision angle,vector(1:3)
      type(atom), intent(inout) :: atoms(:)
      integer i,natoms

      if(talk) print fsubstart,"relax_rigid_molecule_trala" 

      natoms=size(atoms)
      sum_of_forces=0.0d0

      do i=1,natoms
        sum_of_forces=sum_of_forces+atoms(i)%force
      end do ! i

      ! rigid tranlation
      print '(8x,"translating molecule by ",3(F12.6,x)," Angs ")',        &
     &  sum_of_forces*step/dble(natoms) 
      do i=1,natoms
        atoms(i)%abswhere=atoms(i)%abswhere+sum_of_forces*step            &
     &    /dble(natoms)
        call abs2frac(atoms(i)%abswhere,vecs,atoms(i)%where)
      end do
      !  
      if(talk) print fsubendext, "relax_rigid_molecule_trala"
      return

      end subroutine relax_rigid_molecule_trala

c---------------------------------------------------------------------

      subroutine relax_rigid_molecule_breath(atoms,vecs,step)
      use defs  
      use misc, only : abs2frac
      implicit none
      double precision, intent(in) :: step,vecs(1:3,1:3) 
      double precision av_radial_force,radial_force,com(1:3)
      double precision delta_r(1:3),tol
      type(atom), intent(inout) :: atoms(:)
      integer i,natoms

      if(talk) print fsubstart,"relax_rigid_molecule_breath" 

      natoms=size(atoms)
      av_radial_force=0.0d0
      tol=1.0E-3 ! in Angs. Only atoms with a minimum distance from the center of molecule are relaxed (numerical reasons)

      ! get center of molecule (com) 
      com=0.0d0
      do i=1,natoms
        com=com+atoms(i)%abswhere
      end do 
      com=com/dble(natoms)
      print '(8x,"center of molecule at ",3(F12.6))',com

      ! calculate average radial force
      do i=1,natoms
        if (norm2(atoms(i)%abswhere-com).gt.tol) then
          radial_force=dot_product(atoms(i)%force,atoms(i)%abswhere-com)  &
     &               /(norm2(atoms(i)%abswhere-com))**2
          print '(8x,"radial force acts on atom ",I0,x,A2)',              &
     &          i, atoms(i)%name
        else 
          radial_force=0.0d0
        end if
        av_radial_force=av_radial_force + radial_force
      end do ! i
      av_radial_force=av_radial_force/dble(natoms)

      ! rigid breathing
      print '(8x," expanding molecule by ",3(F12.6,x)," Angs ")',         &
     &  av_radial_force*step
      do i=1,natoms
        delta_r=(atoms(i)%abswhere-com)*av_radial_force*step           
        atoms(i)%abswhere=atoms(i)%abswhere+delta_r
        call abs2frac(atoms(i)%abswhere,vecs,atoms(i)%where)
      end do
      !  
      if(talk) print fsubendext, "relax_rigid_molecule_breath"
      return

      end subroutine relax_rigid_molecule_breath

c---------------------------------------------------------------------

      subroutine get_xcoulomb(file1,file2,format1,step,finestep,rmin)
      use defs
      use misc
      use readcoords
      implicit none
      ! calculates the classical Coulomb energy of an exciton
      !
      ! variables
      character(len=*), intent(in) :: file1,file2,format1
      ! internal variables
      type(atom),allocatable :: atoms(:)
      !type(element),allocatable :: species(:)
      !double precision vecs(1:3,1:3),
      double precision xcoul,enorm,hnorm,rvec1(1:3),rvec2(1:3)
      double precision overlap,emom1,emom2,hmom1,hmom2
      double precision ecom(1:3),hcom(1:3),e_delta_r(1:3),h_delta_r(1:3)
      double precision dx,dy,dz,d3r,origin(1:3),rmin,reh
      double precision d3rfine,r0
      double precision, allocatable :: grid(:,:,:,:),wf1(:,:,:)
      double precision, allocatable :: wf2(:,:,:)
      integer i1,i2,i3,j1,j2,j3,nx,ny,nz,step
      integer k1,k2,k3,l1,l2,l3,finestep
      !integer natoms,nspecies
      !integer iatom 
      !
      if(talk) print fsubstart,"get_xcoulomb" 
      if (talk) print '(8x,"files: ",2(A40))',file1,file2
      if (talk) print '(8x,"format: ",A40)',format1
      !
      ! begin definitions 
      print '(8x,"step: ",I0)', step 
      print '(8x,"fine step: ",I0)', finestep 
      print '(8x,"rmin: ",F10.6," * coarse grid spacing")', rmin 
      ! end definitions 
      !
      ! begin read cube file 
      call read_cube(file1,atoms,grid,wf1)
      call read_cube(file2,atoms,grid,wf2)
      deallocate(atoms)
      ! end read cube file
      !
      ! begin initialize
      xcoul=0.0d0
      enorm=0.0d0
      ecom=0.0d0 ! center of mass for e density
      emom1=0.0d0 ! first radial moment for e density
      emom2=0.0d0 ! 2nd radial moment for e density
      hnorm=0.0d0
      hcom=0.0d0 ! center of mass for h density
      hmom1=0.0d0 ! first radial moment for h density
      hmom2=0.0d0 ! 2nd radial moment for h density
      overlap=0.0d0
      nx=size(grid,1)
      ny=size(grid,2)
      nz=size(grid,3)
      dx=(grid(2,1,1,1)-grid(1,1,1,1))*dble(step)  
      dy=(grid(1,2,1,2)-grid(1,1,1,2))*dble(step)  
      dz=(grid(1,1,2,3)-grid(1,1,1,3))*dble(step)  
      origin(1:3)=grid(1,1,1,1:3)
      rmin=rmin*max(max(dx,dy),dz)
      d3r=dx*dy*dz
      d3rfine=d3r/dble((2*finestep)**3)
      print '(8x,"grid dimension: ",3(I0,1x))',nx,ny,nz
      ! end initialize
      !
      ! begin calculate density norms and center of mass 
      do i1=1,nx,step
!       print '(8x,"current progress: ",1(F6.1)," %")',                    &
!     &        dble(i1-1)/dble(nx)*100.0d0
       do i2=1,ny,step
        do i3=1,nz,step 
         rvec1(1:3) = grid(i1,i2,i3,1:3) 
         enorm=enorm+wf1(i1,i2,i3)*d3r
         ecom=ecom+wf1(i1,i2,i3)*rvec1(1:3)*d3r
         hnorm=hnorm+wf2(i1,i2,i3)*d3r 
         hcom=hcom+wf2(i1,i2,i3)*rvec1(1:3)*d3r
        end do ! i3
       end do ! i2
      end do ! i1
      ! end calculate density norms and center of mass 
      !
      ! begin calculate density moments 
      do i1=1,nx,step
!       print '(8x,"current progress: ",1(F6.1)," %")',                    &
!     &        dble(i1-1)/dble(nx)*100.0d0
       do i2=1,ny,step
        do i3=1,nz,step 
         rvec1(1:3) = grid(i1,i2,i3,1:3) 
         e_delta_r(1:3)=rvec1(1:3)-ecom(1:3)
         h_delta_r(1:3)=rvec1(1:3)-hcom(1:3)
         emom1=emom1+wf1(i1,i2,i3)*absvec(e_delta_r)*d3r
         emom2=emom2+wf1(i1,i2,i3)*dot_product(e_delta_r,e_delta_r)*d3r
         hmom1=hmom1+wf2(i1,i2,i3)*absvec(h_delta_r)*d3r
         hmom2=hmom2+wf2(i1,i2,i3)*dot_product(h_delta_r,h_delta_r)*d3r
        end do ! i3
       end do ! i2
      end do ! i1
      !
      print '(8x,"e WF norm:",1(F10.6))',enorm
      print '(8x,"e WF com:",3(F10.6))',ecom
      print '(8x,"e WF mom1:",1(F10.6))',emom1
      print '(8x,"e WF mom2:",3(F10.6))',emom2
      print '(8x,"h WF norm:",1(F10.6))',hnorm
      print '(8x,"h WF com:",3(F10.6))',hcom
      print '(8x,"h WF mom1:",1(F10.6))',hmom1
      print '(8x,"h WF mom2:",3(F10.6))',hmom2
      ! end calculate density moments 
      !
      !
      if (finestep<0) then
        if(allocated(atoms)) deallocate(atoms)
        if(allocated(grid)) deallocate(grid)
        if(allocated(wf1)) deallocate(wf1)
        if(allocated(wf2)) deallocate(wf2)
        if(talk) print fsubendext, "get_xcoulomb"
        return
      end if 
      !
      ! begin calculate Coulomb energy and overlaps 
!      enorm=0.0d0
!      hnorm=0.0d0 
      do i1=1,nx,step
       print '(8x,"current progress: ",1(F6.1)," %")',                    &
     &        dble(i1-1)/dble(nx)*100.0d0
       do i2=1,ny,step
        !print '(8x,"current y grid index: ",1(I0))',i2
        do i3=1,nz,step
         rvec1(1:3) = grid(i1,i2,i3,1:3) 
!         enorm=enorm+wf1(i1,i2,i3)*d3r
!         hnorm=hnorm+wf2(i1,i2,i3)*d3r 
         overlap=overlap+sqrt(abs(wf1(i1,i2,i3)*wf2(i1,i2,i3)))*d3r
         do j1=1,nx,step
          do j2=1,ny,step
           do j3=1,nz,step
!            if (.not.(i1==j1.and.i2==j2.and.i3==j3)) then
              rvec2(1:3) = grid(j1,j2,j3,1:3) 
              reh=dot_product(rvec1-rvec2,rvec1-rvec2)**0.5d0
              if (reh.gt.rmin) then
                xcoul=xcoul-wf1(i1,i2,i3)*wf2(j1,j2,j3)*d3r**2/reh      
              else
               if (.not.(i1==j1.and.i2==j2.and.i3==j3)) then
                ! treat region with r_eh<=rmin with fine grid 
                !print '(8x,"Entering fine grid")'
                do k1=-finestep,finestep-1                        
                 do k2=-finestep,finestep-1                        
                  do k3=-finestep,finestep-1                        
                   do l1=-finestep,finestep-1                        
                    do l2=-finestep,finestep-1                        
                     do l3=-finestep,finestep-1                        
                      rvec1(1:3)=grid(i1,i2,i3,1:3)                       &
     &                       +(/dble(k1)*dx,dble(k2)*dy,dble(k3)*dz/)     &
     &                       /dble(2*finestep) 
                      rvec2(1:3)=grid(j1,j2,j3,1:3)                       &
     &                       +(/dble(l1)*dx,dble(l2)*dy,dble(l3)*dz/)     &
     &                       /dble(2*finestep) 
                      reh=dot_product(rvec1-rvec2,rvec1-rvec2)**0.5d0
                      if(reh.gt.0.0d0) then
                        xcoul=xcoul-wf1(i1,i2,i3)*wf2(j1,j2,j3)           &
     &                     *d3rfine**2/reh
                      end if
                     end do
                    end do
                   end do
                  end do
                 end do
                end do
               else
                ! use spherical average near r=0:
                r0=(3.0d0*d3r/(4.0d0*Pi))**(1.0d0/3.0d0)
                xcoul=xcoul-wf1(i1,i2,i3)*wf2(j1,j2,j3)*2.0d0*Pi*r0**2
               end if ! (.not.(i1==j1.and.i2==j2.and.i3==j3))
!              end if
            end if
           end do
          end do
         end do
        end do
       end do
      end do
      print '(8x,"E_Coulomb(X) in Hartree:",1(F10.6))',xcoul
      print '(8x,"E_Coulomb(X) in eV:",1(F10.6))',xcoul*hartree
      print '(8x,"e WF norm:",1(F10.6))',enorm
      print '(8x,"h WF norm:",1(F10.6))',hnorm
      print '(8x,"e-h density overlap:",1(F10.6))',overlap
      ! begin calculate Coulomb energy and overlaps 
      !
      ! end normally:
      deallocate(grid,wf1,wf2)
      if(talk) print fsubendext, "get_xcoulomb"
      return
      
      end subroutine get_xcoulomb

c---------------------------------------------------------------------

      subroutine get_ekin(file1,format1,step)
      use defs
      use misc
      use readcoords
      implicit none
      ! calculates the kinetic energy of a WF
      !
      ! variables
      character(len=*), intent(in) :: file1,format1
      ! internal variables
      type(atom),allocatable :: atoms(:)
      !type(element),allocatable :: species(:)
      !double precision vecs(1:3,1:3),
      double precision ekin,enorm
      double precision dx,dy,dz,d3r,origin(1:3)
      double precision r0
      double precision nabla_psi(1:3)
      double precision, allocatable :: grid(:,:,:,:)
      double precision, allocatable :: wf1(:,:,:)
      integer i1,i2,i3,nx,ny,nz,step
      !
      if(talk) print fsubstart,"get_ekin" 
      if (talk) print '(8x,"file: ",1(A40))',file1
      if (talk) print '(8x,"format: ",A40)',format1
      !
      ! begin definitions 
      print '(8x,"step: ",I0)', step 
      ! end definitions 
      !
      ! begin read cube file 
      call read_cube(file1,atoms,grid,wf1)
      deallocate(atoms)
      ! end read cube file
      !
      ! begin calculate ekin
      ekin=0.0d0
      enorm=0.0d0
      nx=size(grid,1)
      ny=size(grid,2)
      nz=size(grid,3)
      dx=(grid(2,1,1,1)-grid(1,1,1,1))*dble(step)  
      dy=(grid(1,2,1,2)-grid(1,1,1,2))*dble(step)  
      dz=(grid(1,1,2,3)-grid(1,1,1,3))*dble(step)  
      origin(1:3)=grid(1,1,1,1:3)
      d3r=dx*dy*dz
      print '(8x,"grid dimension: ",3(I0,1x))',nx,ny,nz
      do i1=2,nx-1,step
       print '(8x,"current progress: ",1(F6.1)," %")',                    &
     &        dble(i1-1)/dble(nx)*100.0d0
       do i2=2,ny-1,step
        do i3=2,nz-1,step
         enorm=enorm+wf1(i1,i2,i3)*d3r
         nabla_psi(1:3)=(/(sqrt(abs(wf1(i1+1,i2,i3)))                     &
     &                   -sqrt(abs(wf1(i1-1,i2,i3))))/dx,                 &
     &                    (sqrt(abs(wf1(i1,i2+1,i3)))                     &
     &                   -sqrt(abs(wf1(i1,i2-1,i3))))/dy,                 &
     &                    (sqrt(abs(wf1(i1,i2,i3+1)))                     &
     &                   -sqrt(abs(wf1(i1,i2,i3-1))))/dz/)/2.0d0
         ekin=ekin+0.50d0*dot_product(nabla_psi,nabla_psi)*d3r  
        end do
       end do
      end do
      print '(8x,"E_kin in Hartree:",1(F10.6))',ekin
      print '(8x,"E_kin in eV:",1(F10.6))',ekin*hartree
      print '(8x,"e WF norm:",1(F10.6))',enorm
      ! end calculate xcoul
      !
      ! end normally:
      deallocate(grid,wf1)
      if(talk) print fsubendext, "get_ekin"
      return
      
      end subroutine get_ekin

c---------------------------------------------------------------------

      subroutine read_cube(file1,atoms,grid,wf)
      use defs
      implicit none
      ! reads WF or density from cube file
      ! variables
      character(len=*), intent(in) :: file1
      double precision, allocatable :: wf(:,:,:),grid(:,:,:,:)
      type(atom),allocatable :: atoms(:)
      !type(element),allocatable :: species(:)
      ! internal variables
      character line*1024
      integer natoms,nx(1:3)
      integer iatom,atnum,i1,i2,i3
      double precision xmin(1:3),amin,rdum
      !double precision xgrid(1:3)
     
      ! ************************************************ 
      ! read cube file header
      !
      open(51,file=file1,err=100)
      read(51,'(A256)',err=100,end=100) line
      read(51,'(A256)',err=100,end=100) line
      !write (10+it,'(A)') "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
      read(51,'(I6,F10.5,F10.5,F10.5)') natoms,xmin(1),xmin(2),xmin(3)
      read(51,'(I6,F10.5,F10.5,F10.5)') nx(1),amin,rdum,rdum
      read(51,'(I6,F10.5,F10.5,F10.5)') nx(2),rdum,amin,rdum
      read(51,'(I6,F10.5,F10.5,F10.5)') nx(3),rdum,rdum,amin
      ! read atoms positions
      if(allocated(atoms)) deallocate(atoms)
      allocate(atoms(natoms)) 
      do iatom=1,natoms
        read(51,'(I6,F10.5,F10.5,F10.5,F10.5)') atnum,rdum,               & 
     &       atoms(iatom)%abswhere(1:3)
        atoms(iatom)%abswhere=atoms(iatom)%abswhere*bohr
        atoms(iatom)%name=elements(atnum)
      end do
      !
      print '(8x,"natoms: ",I0)', natoms
      !
      ! End read cube file header
      ! ************************************************ 
      
            
      ! ************************************************ 
      ! loop on the grid
      if(allocated(grid)) deallocate(grid) 
      allocate(grid(nx(1),nx(2),nx(3),1:3))
      allocate(wf(nx(1),nx(2),nx(3)))
      !
      do i1=1,nx(1)
        do i2=1,nx(2)
          do i3=1,nx(3)
            
            ! get grid point position
            grid(i1,i2,i3,1) = xmin(1) + dble(i1-1)*amin
            grid(i1,i2,i3,2) = xmin(2) + dble(i2-1)*amin
            grid(i1,i2,i3,3) = xmin(3) + dble(i3-1)*amin
    
            ! get wf or density
            read(51,'(E14.6)') wf(i1,i2,i3)
          
          end do
        end do
      end do
      
      close(51)
      return
   
100   nerr=nerr+1
      close(51)
      print ferrmssg,'Coul not open/read cube file.'
      return

      end subroutine read_cube

c---------------------------------------------------------------------

      subroutine get_vasp_bandstructure(lattice,nskip)
      use defs
      use misc, only: frac2abs
      implicit none 
      character*1024 :: eigenval,outcar,line,xtics
      character*1024 :: formatstring
      integer nskip ! number of kpoints to skip in plot at the beginning of the kpoint list
      character lattice ! lattice type, e.g. hex,cubic..
      integer nbands,iband,nkpoints,ikpoint,nele,iline,idum
      integer iwrite,iread,i1
      logical spinpol,has_G,has_X,has_R,has_M
      double precision, allocatable :: eigenvalues_up(:,:)
      double precision, allocatable :: eigenvalues_down(:,:)
      double precision, allocatable :: k(:,:),kabs(:,:),svec(:)
      double precision s,tol,efermi,gvecs(3,3),fdum
      !
      print fsubstart,'get_vasp_bandstructure'  
      !
      ! begin initialize
      !
      eigenval="EIGENVAL"
      outcar="OUTCAR"
      s=0.0d0
      spinpol=.false.
      has_G=.false.
      has_X=.false.
      has_R=.false.
      has_M=.false.
      tol=1D-4
      !
      ! end initialize
      !
      !
      ! begin read Fermi energy and reciprocal lattice vectors
      !
      open(51,file=outcar,status='old')
10    read(51,'(A256)',err=100,end=20) line
      if(index(line,'E-fermi').gt.0) then
        iread=index(line,'E-fermi')+9
        read(line(iread:),*) efermi
      end if
      if (index(line,'reciprocal lattice vectors').gt.0) then
        do i1=1,3
          read (51,*,end=20,err=100) fdum,fdum,fdum,gvecs(i1,1:3)
        end do
        gvecs=gvecs*2.0d0*pi
      end if 
      goto 10
      !
20    continue
      close(51)
      print '(8x,"reciprocal lattice vecs:",3(/,8x,3(F12.6)))',           &
     &     (gvecs(i1,1:3),i1=1,3)
      !
      ! end read Fermi energy and reciprocal lattice vectors
      !

      open(51,file=eigenval, status='old')
      !
      ! begin read nbands, nkpoints
      !    
      read(51,*,end=100,err=100) idum,idum,idum, idum  
      if (idum==2) spinpol=.true.   
      do iline=1,4
        read(51,'(A1024)',end=100,err=100) line    
      end do
      read(51,*,end=100,err=100) nele,nkpoints,nbands
      !print '(8x,3(1x,I0))',nele,nkpoints,nbands
      allocate(k(1:nkpoints,1:3),eigenvalues_up(nkpoints,nbands))
      allocate(kabs(1:nkpoints,1:3))
      allocate(svec(1:nkpoints))
      svec=0.0d0
      if(spinpol) allocate(eigenvalues_down(nkpoints,nbands))
      !
      ! end read nbands, nkpoints
      !
      !
      ! begin read kpoints and eigenvalues
      !    
      do ikpoint=1,nkpoints
        read(51,'(A1024)',end=100,err=100) line    
        read(51,*,end=100,err=100) k(ikpoint,1:3)
        call frac2abs(k(ikpoint,:),gvecs,kabs(ikpoint,:))
        do iband=1,nbands
          if (.not.spinpol) then
            read(51,*,err=100,end=100) idum,                              &
     &         eigenvalues_up(ikpoint,iband)
          else
            read(51,*,err=100,end=100) idum,                              &
     &         eigenvalues_up(ikpoint,iband),                             &
     &         eigenvalues_down(ikpoint,iband)
          end if
        end do
      end do
      close(51)
      !
      ! end read kpoints and eigenvalues
      !
      !
      ! begin print eigenvalues
      !    
      open(51,file='BS_UP.DAT',status='replace')
      formatstring=''
      write(formatstring,'("(",I0,"(F20.12))")') nbands+4
      do ikpoint=1+nskip,nkpoints
        !if (ikpoint.gt.1) s=s+norm2(k(ikpoint,1:3)-k(ikpoint-1,1:3))
        if (ikpoint.gt.1) s=s+norm2(kabs(ikpoint,:)-kabs(ikpoint-1,:))
        if (ikpoint.gt.1) svec(ikpoint)=s
!        write(51,*) svec(ikpoint),k(ikpoint,1:3),                         &
!     & eigenvalues_up(ikpoint,1:nbands)-efermi
        write(51,formatstring) svec(ikpoint),k(ikpoint,1:3),              &
     & eigenvalues_up(ikpoint,1:nbands)-efermi
      end do
      close(51)
      if (spinpol) then
        open(51,file='BS_DOWN.DAT',status='replace')
        do ikpoint=1+nskip,nkpoints
          !if (ikpoint.gt.1) s=s+norm2(k(ikpoint,1:3)-k(ikpoint-1,1:3))
          !if (ikpoint.gt.1) s=s+norm2(kabs(ikpoint,:)-kabs(ikpoint-1,:))
!          write(51,*) svec(ikpoint),k(ikpoint,1:3),                       &
!     &   eigenvalues_down(ikpoint,1:nbands)-efermi
          write(51,formatstring) svec(ikpoint),k(ikpoint,1:3),            &
     &   eigenvalues_down(ikpoint,1:nbands)-efermi
        end do
        close(51)
      end if
      !
      ! end print eigenvalues
      !    
      ! 
      ! begin print gnuplot file
      !
      open(51,file='BS.gpi',status='replace')
      write(51,'(256A)') 'set term epslatex standalone color "Times-Roma  &
     &n" font 8' 
      write(51,'(256A)') 'set out "BS.tex"' 
      write(51,'(256A)') '#set term post enhanced color "Times-Roman"' 
      write(51,'(256A)') '#set out BS.eps' 
      write(51,'(256A)') "set size 0.75" 
      write(51,'(256A)') "set ylabel '$E$ (eV)' "
      write(51,'(256A)') "ymin=-10" 
      write(51,'(256A)') "ymax=10" 
      write(51,'(256A)') "set yrange [ymin:ymax]" 
      ! 
      ! begin get special points
      !
      xtics='('
      do ikpoint=1+nskip,nkpoints
        !
        ! Gamma
        !
        if (norm2(k(ikpoint,1:3)-(/0.0,0.0,0.0/)).lt.tol)  then         
          has_G=.true.
          iwrite=len_trim(xtics)+1
          write(xtics(iwrite:),                                           &
!     &     '(A15,1x,F10.6,A1)') '"{/Symbol G}" ',svec(ikpoint),','  
     &     '(A15,1x,F10.6,A1)') "'$\Gamma$ '",svec(ikpoint),','  
          write(51,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
     &      ",ymin to ",F10.6,",ymax")')svec(ikpoint),svec(ikpoint) 
        end if
        !
        ! X
        !
        if (norm2(k(ikpoint,1:3)-(/0.0,0.0,0.5/)).lt.tol)  then         
          has_X=.true.
          iwrite=len_trim(xtics)+1
          write(xtics(iwrite:),                                           &
     &     '(A8,1x,F10.6,A1)') '"X" ',svec(ikpoint),','  
          write(51,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
     &      ",ymin to ",F10.6,",ymax")')svec(ikpoint),svec(ikpoint) 
        end if
        if (norm2(k(ikpoint,1:3)-(/0.0,0.5,0.0/)).lt.tol)  then         
          has_X=.true.
          iwrite=len_trim(xtics)+1
          write(xtics(iwrite:),                                           &
     &     '(A8,1x,F10.6,A1)') '"X" ',svec(ikpoint),','  
          write(51,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
     &      ",ymin to ",F10.6,",ymax")')svec(ikpoint),svec(ikpoint) 
        end if
        if (norm2(k(ikpoint,1:3)-(/0.5,0.0,0.0/)).lt.tol)  then         
          has_X=.true.
          iwrite=len_trim(xtics)+1
          write(xtics(iwrite:),                                           &
     &     '(A8,1x,F10.6,A1)') '"X" ',svec(ikpoint),','  
          write(51,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
     &      ",ymin to ",F10.6,",ymax")')svec(ikpoint),svec(ikpoint) 
        end if
        !
        ! M
        !
        if (norm2(k(ikpoint,1:3)-(/0.0,0.5,0.5/)).lt.tol)  then         
          has_M=.true.
          iwrite=len_trim(xtics)+1
          write(xtics(iwrite:),                                           &
     &     '(A8,1x,F10.6,A1)') '"M" ',svec(ikpoint),','  
          write(51,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
     &      ",ymin to ",F10.6,",ymax")')svec(ikpoint),svec(ikpoint) 
        end if
        if (norm2(k(ikpoint,1:3)-(/0.5,0.0,0.5/)).lt.tol)  then         
          has_M=.true.
          iwrite=len_trim(xtics)+1
          write(xtics(iwrite:),                                           &
     &     '(A8,1x,F10.6,A1)') '"M" ',svec(ikpoint),','  
          write(51,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
     &      ",ymin to ",F10.6,",ymax")')svec(ikpoint),svec(ikpoint) 
        end if
        if (norm2(k(ikpoint,1:3)-(/0.5,0.5,0.0/)).lt.tol)  then         
          has_M=.true.
          iwrite=len_trim(xtics)+1
          write(xtics(iwrite:),                                           &
     &     '(A8,1x,F10.6,A1)') '"M" ',svec(ikpoint),','  
          write(51,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
     &      ",ymin to ",F10.6,",ymax")')svec(ikpoint),svec(ikpoint) 
        end if
        !
        ! R
        !
        if (norm2(k(ikpoint,1:3)-(/0.5,0.5,0.5/)).lt.tol)  then         
          has_R=.true.
          iwrite=len_trim(xtics)+1
          write(xtics(iwrite:),                                           &
     &     '(A8,1x,F10.6,A1)') '"R" ',svec(ikpoint),','  
          write(51,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
     &      ",ymin to ",F10.6,",ymax")')svec(ikpoint),svec(ikpoint) 
        end if
      end do
      iwrite=len_trim(xtics)
      write(xtics(iwrite:),'(A1)') ')'  
      !print*, xtics
      ! 
      ! end get special points
      !
      !
      write(51,'(256A)') "set xtics ",xtics 
      write(51,'(256A)') "set samples 1000" 
      write(51,'(256A)') "N=`awk 'NR==1 {print NF}' BS_UP.DAT`"
      write(51,'(256A)')"plot for [i=5:N] 'BS_UP.DAT' u 1:i smooth cspli  &
     &nes w l lt 1 lc 1 t ''"
      write(51,'(256A)') "set out" 
      write(51,'(256A)') "system('pdflatex BS.tex')" 
      close(51)
      !
      ! end print gnuplot file
      !
      print fsubendext,'get_vasp_bandstructure'  
      ! 
      return
      !
100   continue    
      close(51)
      nerr=nerr+1
      print ferrmssg,'Problem with EIGENVAL file'
      !
      end subroutine get_vasp_bandstructure

c---------------------------------------------------------------------

      subroutine get_vasp_projected_bandstructure()
      use defs
      use misc, only: frac2abs
      use readcoords
      implicit none 
      character*1024 :: eigenval,procar,outcar,line,xtics
      character*1024 :: formatstring
      character filename*256,filename2*256
      integer nbands,iband,nkpoints,ikpoint,nele,iline,idum
      integer iwrite,iread,nspecies,ispecies,natoms,iatom,i1
      logical spinpol,has_G,has_X,has_R,has_M,noncollinear
      double precision, allocatable :: eigenvalues_up(:,:)
      double precision, allocatable :: eigenvalues_down(:,:)
      double precision, allocatable :: proj_up(:,:,:,:)
      double precision, allocatable :: proj_ele_up(:,:,:,:)
      double precision, allocatable :: proj_down(:,:,:,:)
      double precision, allocatable :: proj_ele_down(:,:,:,:)
      double precision, allocatable :: k(:,:),kabs(:,:),svec(:)
      double precision s,tol,efermi,vecs(1:3,1:3),gvecs(1:3,1:3),fdum
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      !
      print fsubstart,'get_vasp_projected_bandstructure'  
      !
      ! begin initialize
      !
      eigenval="EIGENVAL"
      procar="PROCAR"
      outcar="OUTCAR"
      s=0.0d0
      spinpol=.false.
      noncollinear=.false.
      has_G=.false.
      has_X=.false.
      has_R=.false.
      has_M=.false.
      tol=1D-4
      !
      ! end initialize
      !
      !
      ! begin read atom info
      !
      call read_coords(outcar,"outcar",atoms,natoms,species,              &
     &                 nspecies,vecs)
      !
      ! end read atom info
      !
      ! begin read Fermi energy and reciprocal lattice vectors
      !
      open(51,file=outcar,status='old')
10    read(51,'(A256)',err=100,end=20) line
      if(index(line,'E-fermi').gt.0) then
        iread=index(line,'E-fermi')+9
        read(line(iread:),*) efermi
        !print '(8x,"E-fermi = ",F12.6, " eV")',efermi ! DEBUG
      end if
      if (index(line,'reciprocal lattice vectors').gt.0) then
        do i1=1,3
          read (51,*,end=20,err=100) fdum,fdum,fdum,gvecs(i1,1:3)
        end do
        gvecs=gvecs*2.0d0*pi
      end if 
      if(index(line,' LNONCOLLINEAR ').gt.0) then
        iread=index(line,'=')+1
        read(line(iread:),*) noncollinear
      end if
      goto 10
      !
20    continue
      close(51)
      print '(8x,"reciprocal lattice vecs:",3(/,8x,3(F12.6)))',           &
     &     (gvecs(i1,1:3),i1=1,3)
      !
      ! end read Fermi energy and reciprocal lattice vectors
      !
      !
      open(51,file=eigenval, status='old')
      !
      ! begin read nbands, nkpoints
      !    
      read(51,*,end=102,err=102) idum,idum,idum, idum  
      if (idum==2) spinpol=.true.   
!      print '(8x,"spinpol=",L1)',spinpol ! DEBUG
      do iline=1,4
        read(51,'(A1024)',end=102,err=102) line    
      end do
      read(51,*,end=102,err=102) nele,nkpoints,nbands
!      print '(8x,"nele,nkpoints,nbands: ",3(1x,I0))',nele,nkpoints,       &
!     &      nbands ! DEBUG
      allocate(k(1:nkpoints,1:3),eigenvalues_up(nkpoints,nbands))
      allocate(kabs(1:nkpoints,1:3))
      allocate(proj_up(nkpoints,nbands,natoms,1:10))
      allocate(proj_ele_up(nkpoints,nbands,nspecies,1:10))
      proj_ele_up(1:nkpoints,1:nbands,1:nspecies,1:10)=0.0d0
      allocate(svec(1:nkpoints))
      svec=0.0d0
      if (spinpol) then
        allocate(eigenvalues_down(nkpoints,nbands))
        allocate(proj_down(nkpoints,nbands,natoms,1:10))
        allocate(proj_ele_down(nkpoints,nbands,nspecies,1:10))
        proj_ele_down(1:nkpoints,1:nbands,1:nspecies,1:10)=0.0d0
      end if
      !
      ! end read nbands, nkpoints
      !
      !
      ! begin read kpoints and eigenvalues
      !    
      do ikpoint=1,nkpoints
        read(51,'(A1024)',end=102,err=102) line    
        read(51,*,end=102,err=102) k(ikpoint,1:3)
        call frac2abs(k(ikpoint,:),gvecs,kabs(ikpoint,:))
!        print '(8x,"reading EV for kpoint # ",I0)',ikpoint ! DEBUG
        do iband=1,nbands
!          print '(8x,"reading EV for band # ",I0)',iband ! DEBUG
          if (.not.spinpol) then
            read(51,*,err=102,end=102) idum,                              &
     &         eigenvalues_up(ikpoint,iband)
!               print '(8x,"EV =",F12.6," eV")',eigenvalues_up(ikpoint,    &
!     &iband) ! DEBUG
          else
            read(51,*,err=102,end=102) idum,                              &
     &         eigenvalues_up(ikpoint,iband),                             &
     &         eigenvalues_down(ikpoint,iband)
               print '(8x,"EV =",2(F12.6,x)," eV")',eigenvalues_up(       &
     &ikpoint,iband),eigenvalues_down(ikpoint,iband) ! DEBUG
          end if
        end do
      end do
      close(51)
      !
      ! end read kpoints and eigenvalues
      !
      !
      ! begin print eigenvalues
      !    
      open(51,file='BS_UP.DAT',status='replace')
      formatstring=''
      write(formatstring,'("(",I0,"(F20.12))")') nbands+4
      !print '(8x,"formatstring=",A50)',formatstring
      do ikpoint=1,nkpoints
        !if (ikpoint.gt.1) s=s+norm2(k(ikpoint,1:3)-k(ikpoint-1,1:3))
        if (ikpoint.gt.1) s=s+norm2(kabs(ikpoint,:)-kabs(ikpoint-1,:))
        if (ikpoint.gt.1) svec(ikpoint)=s
        write(51,formatstring) svec(ikpoint),k(ikpoint,1:3),              &
     & eigenvalues_up(ikpoint,1:nbands)-efermi
      end do
      close(51)
      if (spinpol) then
        open(51,file='BS_DOWN.DAT',status='replace')
        do ikpoint=1,nkpoints
          write(51,formatstring) svec(ikpoint),k(ikpoint,1:3),            &
     &   eigenvalues_down(ikpoint,1:nbands)-efermi
        end do
        close(51)
      end if
      !
      ! end print eigenvalues
      !    
      ! 
      ! begin print gnuplot file
      !
      open(51,file='BS.gpi',status='replace')
      open(52,file='BS_PROJ.gpi',status='replace')
      write(51,'(256A)') 'set term x11 "Times-Roman"' 
      write(52,'(256A)') 'set term x11 "Times-Roman"' 
      write(51,'(256A)') '#set term post enhanced color "Times-Roman"' 
      write(52,'(256A)') 'set term post enhanced color "Times-Roman"' 
      write(51,'(256A)') '#set out "BS.eps"' 
      !write(52,'(256A)') '#set out "BS_PROJ.eps"' 
      write(51,'(256A)') "set size 0.25, 0.5" 
      write(52,'(256A)') "set size 0.25, 0.5" 
      write(51,'(256A)') "set ylabel 'E-E_F (eV)' offset 1"
      write(52,'(256A)') "set ylabel 'E-E_F (eV)' offset 1"
      write(51,'(256A)') "ymin=-10" 
      write(52,'(256A)') "ymin=-10" 
      write(51,'(256A)') "ymax=10" 
      write(52,'(256A)') "ymax=10" 
      write(51,'(256A)') "set yrange [ymin:ymax]" 
      write(52,'(256A)') "set yrange [ymin:ymax]" 
      write(52,'(256A)') "set palette negative" 
      ! 
      ! begin get special points
      !
      xtics='('
      do ikpoint=1,nkpoints
        if (norm2(k(ikpoint,1:3)-(/0.0,0.0,0.0/)).lt.tol)  then         
          has_G=.true.
          iwrite=len_trim(xtics)+1
          write(xtics(iwrite:),                                           &
     &     '(A15,1x,F10.6,A1)') '"{/Symbol G}" ',svec(ikpoint),','  
          write(51,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
     &      ",ymin to ",F10.6,",ymax")')svec(ikpoint),svec(ikpoint) 
          write(52,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
     &      ",ymin to ",F10.6,",ymax")')svec(ikpoint),svec(ikpoint) 
        end if
        if (norm2(k(ikpoint,1:3)-(/0.0,0.0,0.5/)).lt.tol)  then         
          has_X=.true.
          iwrite=len_trim(xtics)+1
          write(xtics(iwrite:),                                           &
     &     '(A8,1x,F10.6,A1)') '"X" ',svec(ikpoint),','  
          write(51,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
     &      ",ymin to ",F10.6,",ymax")')svec(ikpoint),svec(ikpoint) 
          write(52,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
     &      ",ymin to ",F10.6,",ymax")')svec(ikpoint),svec(ikpoint) 
        end if
        if (norm2(k(ikpoint,1:3)-(/0.0,0.5,0.0/)).lt.tol)  then         
          has_X=.true.
          iwrite=len_trim(xtics)+1
          write(xtics(iwrite:),                                           &
     &     '(A8,1x,F10.6,A1)') '"X" ',svec(ikpoint),','  
          write(51,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
     &      ",ymin to ",F10.6,",ymax")')svec(ikpoint),svec(ikpoint) 
          write(52,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
     &      ",ymin to ",F10.6,",ymax")')svec(ikpoint),svec(ikpoint) 
        end if
        if (norm2(k(ikpoint,1:3)-(/0.5,0.0,0.0/)).lt.tol)  then         
          has_X=.true.
          iwrite=len_trim(xtics)+1
          write(xtics(iwrite:),                                           &
     &     '(A8,1x,F10.6,A1)') '"X" ',svec(ikpoint),','  
          write(51,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
     &      ",ymin to ",F10.6,",ymax")')svec(ikpoint),svec(ikpoint) 
          write(52,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
     &      ",ymin to ",F10.6,",ymax")')svec(ikpoint),svec(ikpoint) 
        end if
        if (norm2(k(ikpoint,1:3)-(/0.0,0.5,0.5/)).lt.tol)  then         
          has_M=.true.
          iwrite=len_trim(xtics)+1
          write(xtics(iwrite:),                                           &
     &     '(A8,1x,F10.6,A1)') '"M" ',svec(ikpoint),','  
          write(51,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
     &      ",ymin to ",F10.6,",ymax")')svec(ikpoint),svec(ikpoint) 
          write(52,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
     &      ",ymin to ",F10.6,",ymax")')svec(ikpoint),svec(ikpoint) 
        end if
        if (norm2(k(ikpoint,1:3)-(/0.5,0.0,0.5/)).lt.tol)  then         
          has_M=.true.
          iwrite=len_trim(xtics)+1
          write(xtics(iwrite:),                                           &
     &     '(A8,1x,F10.6,A1)') '"M" ',svec(ikpoint),','  
          write(51,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
     &      ",ymin to ",F10.6,",ymax")')svec(ikpoint),svec(ikpoint) 
          write(52,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
     &      ",ymin to ",F10.6,",ymax")')svec(ikpoint),svec(ikpoint) 
        end if
        if (norm2(k(ikpoint,1:3)-(/0.5,0.5,0.0/)).lt.tol)  then         
          has_M=.true.
          iwrite=len_trim(xtics)+1
          write(xtics(iwrite:),                                           &
     &     '(A8,1x,F10.6,A1)') '"M" ',svec(ikpoint),','  
          write(51,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
     &      ",ymin to ",F10.6,",ymax")')svec(ikpoint),svec(ikpoint) 
          write(52,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
     &      ",ymin to ",F10.6,",ymax")')svec(ikpoint),svec(ikpoint) 
        end if
        if (norm2(k(ikpoint,1:3)-(/0.5,0.5,0.5/)).lt.tol)  then         
          has_R=.true.
          iwrite=len_trim(xtics)+1
          write(xtics(iwrite:),                                           &
     &     '(A8,1x,F10.6,A1)') '"R" ',svec(ikpoint),','  
          write(51,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
     &      ",ymin to ",F10.6,",ymax")')svec(ikpoint),svec(ikpoint) 
          write(52,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
     &      ",ymin to ",F10.6,",ymax")')svec(ikpoint),svec(ikpoint) 
        end if
      end do
      iwrite=len_trim(xtics)
      write(xtics(iwrite:),'(A1)') ')'  
      !print*, xtics
      ! 
      ! end get special points
      !
      !
      write(51,'(256A)') "set xtics ",xtics 
      write(52,'(256A)') "set xtics ",xtics 
      write(51,'(256A)') "N=`awk 'NR==1 {print NF}' BS_UP.DAT`"
      write(51,'(256A)')"plot for [i=5:N] 'BS_UP.DAT' u 1:i smooth cspli  &
     &nes w l lt 1 lc 1 lw 2 t ''"
      close(51)
      do ispecies=1,nspecies
        filename="BS_PROJ_UP_"
        filename2="BS_PROJ_UP_"
        write(filename(12:),'(2A)')                                      &
     &    trim(adjustl(species(ispecies)%name ) )
        write(filename2(12:),'(2A)')                                     &
     &    trim(adjustl(species(ispecies)%name ) )
        iwrite=len_trim(filename)+1
        write(filename(iwrite:),'(".DAT")')
        write(filename2(iwrite:),'(".EPS",1A)')"'"
        write(52,'(256A)') "set out '",filename2
        write(52,'(256A,1A)')"N=`awk 'NR==1 {print NF}' ",filename ,"`"
        write(52,'(256A)')"plot for [i=5:N-1:2] '",                      &
     &trim(adjustl(filename)),"' u 1:i:i+1 with lines palette lw 2 t ''"
        write(52,'(256A,1A)') "system('epstopdf ",                        &
     &trim(adjustl(filename2)),")" 
      end do
      close(52)
      !
      ! end print gnuplot file
      !
      !
      ! begin read projections
      !
      open(51,file=procar, status='old')
      !    
      read(51,'(A1024)',end=101,err=101) line    
      !print*,line
      !
      do ikpoint=1,nkpoints
        do iline=1,3
          read(51,'(A1024)',end=101,err=101) line    
          !print*,line
        end do
        do iband=1,nbands
          do iline=1,4
            read(51,'(A1024)',end=101,err=101) line    
            !print*,line
          end do
          do iatom=1,natoms
            read(51,*,end=101,err=101) idum,                              &
     &         proj_up(ikpoint,iband,iatom,1:10)
            !print*,idum,proj_up(ikpoint,iband,iatom,1:10)
          end do ! iatom
          read(51,'(A1024)',end=101,err=101) line    
          ! ifnoncollinear, then the magnetization DOS (3 directions) is
          ! also contained in PROCAR - here we skip this information.
          if (noncollinear) then
            do iatom=1,3*natoms+3
              read(51,*,end=101,err=101) line                             &
            end do ! iatom
          end if ! noncollinear
!          if (.not.spinpol) then
!            read(51,*,err=100,end=100) idum,                              &
!     &         eigenvalues_up(ikpoint,iband)
!          else
!            read(51,*,err=100,end=100) idum,                              &
!     &         eigenvalues_up(ikpoint,iband),                             &
!     &         eigenvalues_down(ikpoint,iband)
!          end if
        end do ! iband
      end do ! ikpoint
      close(51)
      !
      ! end read projections
      !
      !
      ! begin print total atom-projected band structures 
      !    
      formatstring=''
      write(formatstring,'("(",I0,"(F20.12))")') 2*nbands+4
      do iatom=1,natoms
        write(filename,'("BS_PROJ_UP_",I0,".DAT")')iatom
        open(51,file=filename,status='replace')
        do ikpoint=1,nkpoints
          write(51,formatstring) svec(ikpoint),k(ikpoint,1:3),            &
     &    (eigenvalues_up(ikpoint,iband)-efermi,                          &
     &     proj_up(ikpoint,iband,iatom,10),iband=1,nbands)
        end do
        close(51)
      end do
!      if (spinpol) then
!        open(51,file='BS_PROJ_DOWN.DAT',status='replace')
!        do ikpoint=1,nkpoints
!          write(51,*) svec(ikpoint),k(ikpoint,1:3),                       &
!     &   eigenvalues_down(ikpoint,1:nbands)-efermi
!        end do
!        close(51)
!      end if
      !
      ! end print total atom-projected band structures 
      !
      ! begin print atom-projected s band structures 
      !    
      formatstring=''
      write(formatstring,'("(",I0,"(F20.12))")') 2*nbands+4
      do iatom=1,natoms
        write(filename,'("BS_PROJ_UP_s_",I0,".DAT")')iatom
        open(51,file=filename,status='replace')
        do ikpoint=1,nkpoints
          write(51,formatstring) svec(ikpoint),k(ikpoint,1:3),            &
     &    (eigenvalues_up(ikpoint,iband)-efermi,                          &
     &     proj_up(ikpoint,iband,iatom,1),iband=1,nbands)
        end do
        close(51)
      end do
      !
      ! end print atom-projected s band structures 
      !    
      ! begin print atom-projected px band structures 
      !    
      formatstring=''
      write(formatstring,'("(",I0,"(F20.12))")') 2*nbands+4
      do iatom=1,natoms
        write(filename,'("BS_PROJ_UP_px_",I0,".DAT")')iatom
        open(51,file=filename,status='replace')
        do ikpoint=1,nkpoints
          write(51,formatstring) svec(ikpoint),k(ikpoint,1:3),            &
     &    (eigenvalues_up(ikpoint,iband)-efermi,                          &
     &     proj_up(ikpoint,iband,iatom,2),iband=1,nbands)
        end do
        close(51)
      end do
      !
      ! end print atom-projected px band structures 
      !
      ! begin print atom-projected py band structures 
      !    
      formatstring=''
      write(formatstring,'("(",I0,"(F20.12))")') 2*nbands+4
      do iatom=1,natoms
        write(filename,'("BS_PROJ_UP_py_",I0,".DAT")')iatom
        open(51,file=filename,status='replace')
        do ikpoint=1,nkpoints
          write(51,formatstring) svec(ikpoint),k(ikpoint,1:3),            &
     &    (eigenvalues_up(ikpoint,iband)-efermi,                          &
     &     proj_up(ikpoint,iband,iatom,3),iband=1,nbands)
        end do
        close(51)
      end do
      !
      ! end print atom-projected py band structures 
      !
      ! begin print atom-projected pz band structures 
      !    
      formatstring=''
      write(formatstring,'("(",I0,"(F20.12))")') 2*nbands+4
      do iatom=1,natoms
        write(filename,'("BS_PROJ_UP_pz_",I0,".DAT")')iatom
        open(51,file=filename,status='replace')
        do ikpoint=1,nkpoints
          write(51,formatstring) svec(ikpoint),k(ikpoint,1:3),            &
     &    (eigenvalues_up(ikpoint,iband)-efermi,                          &
     &     proj_up(ikpoint,iband,iatom,4),iband=1,nbands)
        end do
        close(51)
      end do
      !
      ! end print atom-projected pz band structures 
      !
      ! begin print atom-projected d_xy band structures 
      !    
      formatstring=''
      write(formatstring,'("(",I0,"(F20.12))")') 2*nbands+4
      do iatom=1,natoms
        write(filename,'("BS_PROJ_UP_d_xy_",I0,".DAT")')iatom
        open(51,file=filename,status='replace')
        do ikpoint=1,nkpoints
          write(51,formatstring) svec(ikpoint),k(ikpoint,1:3),            &
     &    (eigenvalues_up(ikpoint,iband)-efermi,                          &
     &     proj_up(ikpoint,iband,iatom,5),iband=1,nbands)
        end do
        close(51)
      end do
      !
      ! end print atom-projected d_xy band structures 
      !
      ! begin print atom-projected d_yz band structures 
      !    
      formatstring=''
      write(formatstring,'("(",I0,"(F20.12))")') 2*nbands+4
      do iatom=1,natoms
        write(filename,'("BS_PROJ_UP_d_yz_",I0,".DAT")')iatom
        open(51,file=filename,status='replace')
        do ikpoint=1,nkpoints
          write(51,formatstring) svec(ikpoint),k(ikpoint,1:3),            &
     &    (eigenvalues_up(ikpoint,iband)-efermi,                          &
     &     proj_up(ikpoint,iband,iatom,6),iband=1,nbands)
        end do
        close(51)
      end do
      !
      ! end print atom-projected d_yz band structures 
      !
      ! begin print atom-projected d_z2 band structures 
      !    
      formatstring=''
      write(formatstring,'("(",I0,"(F20.12))")') 2*nbands+4
      do iatom=1,natoms
        write(filename,'("BS_PROJ_UP_d_z2_",I0,".DAT")')iatom
        open(51,file=filename,status='replace')
        do ikpoint=1,nkpoints
          write(51,formatstring) svec(ikpoint),k(ikpoint,1:3),            &
     &    (eigenvalues_up(ikpoint,iband)-efermi,                          &
     &     proj_up(ikpoint,iband,iatom,7),iband=1,nbands)
        end do
        close(51)
      end do
      !
      ! end print atom-projected d_z2 band structures 
      !
      ! begin print atom-projected d_xz band structures 
      !    
      formatstring=''
      write(formatstring,'("(",I0,"(F20.12))")') 2*nbands+4
      do iatom=1,natoms
        write(filename,'("BS_PROJ_UP_d_xz_",I0,".DAT")')iatom
        open(51,file=filename,status='replace')
        do ikpoint=1,nkpoints
          write(51,formatstring) svec(ikpoint),k(ikpoint,1:3),            &
     &    (eigenvalues_up(ikpoint,iband)-efermi,                          &
     &     proj_up(ikpoint,iband,iatom,8),iband=1,nbands)
        end do
        close(51)
      end do
      !
      ! end print atom-projected d_xz band structures 
      !
      ! begin print atom-projected d_z2 band structures 
      !    
      formatstring=''
      write(formatstring,'("(",I0,"(F20.12))")') 2*nbands+4
      do iatom=1,natoms
        write(filename,'("BS_PROJ_UP_d_z2_",I0,".DAT")')iatom
        open(51,file=filename,status='replace')
        do ikpoint=1,nkpoints
          write(51,formatstring) svec(ikpoint),k(ikpoint,1:3),            &
     &    (eigenvalues_up(ikpoint,iband)-efermi,                          &
     &     proj_up(ikpoint,iband,iatom,9),iband=1,nbands)
        end do
        close(51)
      end do
      !
      ! end print atom-projected d_z2 band structures 
      !    
      !    
      !
      ! begin print element-projected band structures 
      !
      ! begin total BS projected to species
      !
      do ispecies=1,nspecies
        do ikpoint=1,nkpoints
          do iband=1,nbands
            do iatom=1,natoms
              if(atoms(iatom)%name.eq.species(ispecies)%name)            &
     &        then
                !print*,"species name=",species(ispecies)%name
                proj_ele_up(ikpoint,iband,ispecies,1:10)                  &
     &            =proj_ele_up(ikpoint,iband,ispecies,1:10)               &
     &            +proj_up(ikpoint,iband,iatom,1:10)
                !print*,proj_ele_up(ikpoint,iband,ispecies,1:10)
              end if
            end do ! iatom
          end do ! iband
        end do ! ikpoint
        filename="BS_PROJ_UP_"
        write(filename(12:),'(2A)')                                       &
     &    trim(adjustl(species(ispecies)%name ) )
        iwrite=len_trim(filename)+1
        write(filename(iwrite:),'(".DAT")')
        open(51,file=filename,status='replace')
        do ikpoint=1,nkpoints
          write(51,formatstring) svec(ikpoint),k(ikpoint,1:3),            &
     &    (eigenvalues_up(ikpoint,iband)-efermi,                          &
     &     proj_ele_up(ikpoint,iband,ispecies,10),iband=1,nbands)
        end do
        close(51)
      end do ! ispecies
      !
      ! end total BS projected to atom
      !
      ! begin s BS projected to species
      !
      do ispecies=1,nspecies
        filename="BS_PROJ_UP_s"
        write(filename(12:),'(2A)')                                       &
     &    trim(adjustl(species(ispecies)%name ) )
        iwrite=len_trim(filename)+1
        write(filename(iwrite:),'(".DAT")')
        open(51,file=filename,status='replace')
        do ikpoint=1,nkpoints
          write(51,formatstring) svec(ikpoint),k(ikpoint,1:3),            &
     &    (eigenvalues_up(ikpoint,iband)-efermi,                          &
     &     proj_ele_up(ikpoint,iband,ispecies,1),iband=1,nbands)
        end do
        close(51)
      end do ! ispecies
      !
      ! end s BS projected to species
      !
      ! begin px BS projected to species
      !
      do ispecies=1,nspecies
        filename="BS_PROJ_UP_px"
        write(filename(12:),'(2A)')                                       &
     &    trim(adjustl(species(ispecies)%name ) )
        iwrite=len_trim(filename)+1
        write(filename(iwrite:),'(".DAT")')
        open(51,file=filename,status='replace')
        do ikpoint=1,nkpoints
          write(51,formatstring) svec(ikpoint),k(ikpoint,1:3),            &
     &    (eigenvalues_up(ikpoint,iband)-efermi,                          &
     &     proj_ele_up(ikpoint,iband,ispecies,2),iband=1,nbands)
        end do
        close(51)
      end do ! ispecies
      !
      ! end px BS projected to species
      !
      ! begin py BS projected to species
      !
      do ispecies=1,nspecies
        filename="BS_PROJ_UP_py"
        write(filename(12:),'(2A)')                                       &
     &    trim(adjustl(species(ispecies)%name ) )
        iwrite=len_trim(filename)+1
        write(filename(iwrite:),'(".DAT")')
        open(51,file=filename,status='replace')
        do ikpoint=1,nkpoints
          write(51,formatstring) svec(ikpoint),k(ikpoint,1:3),            &
     &    (eigenvalues_up(ikpoint,iband)-efermi,                          &
     &     proj_ele_up(ikpoint,iband,ispecies,3),iband=1,nbands)
        end do
        close(51)
      end do ! ispecies
      !
      ! end py BS projected to species
      !
      ! begin pz BS projected to species
      !
      do ispecies=1,nspecies
        filename="BS_PROJ_UP_pz"
        write(filename(12:),'(2A)')                                       &
     &    trim(adjustl(species(ispecies)%name ) )
        iwrite=len_trim(filename)+1
        write(filename(iwrite:),'(".DAT")')
        open(51,file=filename,status='replace')
        do ikpoint=1,nkpoints
          write(51,formatstring) svec(ikpoint),k(ikpoint,1:3),            &
     &    (eigenvalues_up(ikpoint,iband)-efermi,                          &
     &     proj_ele_up(ikpoint,iband,ispecies,4),iband=1,nbands)
        end do
        close(51)
      end do ! ispecies
      !
      ! end pz BS projected to species
      !
      ! begin d_xy BS projected to species
      !
      do ispecies=1,nspecies
        filename="BS_PROJ_UP_d_xy"
        write(filename(12:),'(2A)')                                       &
     &    trim(adjustl(species(ispecies)%name ) )
        iwrite=len_trim(filename)+1
        write(filename(iwrite:),'(".DAT")')
        open(51,file=filename,status='replace')
        do ikpoint=1,nkpoints
          write(51,formatstring) svec(ikpoint),k(ikpoint,1:3),            &
     &    (eigenvalues_up(ikpoint,iband)-efermi,                          &
     &     proj_ele_up(ikpoint,iband,ispecies,5),iband=1,nbands)
        end do
        close(51)
      end do ! ispecies
      !
      ! end d_xy BS projected to species
      !
      ! begin d_yz BS projected to species
      !
      do ispecies=1,nspecies
        filename="BS_PROJ_UP_d_yz"
        write(filename(12:),'(2A)')                                       &
     &    trim(adjustl(species(ispecies)%name ) )
        iwrite=len_trim(filename)+1
        write(filename(iwrite:),'(".DAT")')
        open(51,file=filename,status='replace')
        do ikpoint=1,nkpoints
          write(51,formatstring) svec(ikpoint),k(ikpoint,1:3),            &
     &    (eigenvalues_up(ikpoint,iband)-efermi,                          &
     &     proj_ele_up(ikpoint,iband,ispecies,6),iband=1,nbands)
        end do
        close(51)
      end do ! ispecies
      !
      ! end d_yz BS projected to species
      !
      ! begin d_z2 BS projected to species
      !
      do ispecies=1,nspecies
        filename="BS_PROJ_UP_d_z2"
        write(filename(12:),'(2A)')                                       &
     &    trim(adjustl(species(ispecies)%name ) )
        iwrite=len_trim(filename)+1
        write(filename(iwrite:),'(".DAT")')
        open(51,file=filename,status='replace')
        do ikpoint=1,nkpoints
          write(51,formatstring) svec(ikpoint),k(ikpoint,1:3),            &
     &    (eigenvalues_up(ikpoint,iband)-efermi,                          &
     &     proj_ele_up(ikpoint,iband,ispecies,7),iband=1,nbands)
        end do
        close(51)
      end do ! ispecies
      !
      ! end d_z2 BS projected to species
      !
      ! begin d_xz BS projected to species
      !
      do ispecies=1,nspecies
        filename="BS_PROJ_UP_d_xz"
        write(filename(12:),'(2A)')                                       &
     &    trim(adjustl(species(ispecies)%name ) )
        iwrite=len_trim(filename)+1
        write(filename(iwrite:),'(".DAT")')
        open(51,file=filename,status='replace')
        do ikpoint=1,nkpoints
          write(51,formatstring) svec(ikpoint),k(ikpoint,1:3),            &
     &    (eigenvalues_up(ikpoint,iband)-efermi,                          &
     &     proj_ele_up(ikpoint,iband,ispecies,8),iband=1,nbands)
        end do
        close(51)
      end do ! ispecies
      !
      ! end d_xz BS projected to species
      !
      ! begin d_x2 BS projected to species
      !
      do ispecies=1,nspecies
        filename="BS_PROJ_UP_d_x2"
        write(filename(12:),'(2A)')                                       &
     &    trim(adjustl(species(ispecies)%name ) )
        iwrite=len_trim(filename)+1
        write(filename(iwrite:),'(".DAT")')
        open(51,file=filename,status='replace')
        do ikpoint=1,nkpoints
          write(51,formatstring) svec(ikpoint),k(ikpoint,1:3),            &
     &    (eigenvalues_up(ikpoint,iband)-efermi,                          &
     &     proj_ele_up(ikpoint,iband,ispecies,9),iband=1,nbands)
        end do
        close(51)
      end do ! ispecies
      !
      ! end d_x2 BS projected to species
      !
      ! end print element-projected band structures 
      !
      print fsubendext,'get_vasp_projected_bandstructure'  
      ! 
      return
      !
100   continue    
      close(51)
      nerr=nerr+1
      print ferrmssg,'Problem with OUTCAR file'
      !
101   continue    
      close(51)
      nerr=nerr+1
      print ferrmssg,'Problem with PROCAR file'
      !
102   continue    
      close(51)
      nerr=nerr+1
      print ferrmssg,'Problem with EIGENVAL file'
      !
      end subroutine get_vasp_projected_bandstructure

c---------------------------------------------------------------------

      subroutine get_vasp_projected_dos(nlayers,direction,origin,         &
     &             dos_tol)
      use defs
      use readcoords
      use misc, only : theta_function
      implicit none 
      character*1024 :: eigenval,procar,outcar,doscar,line,xtics
      character filename*256,filename2*256, string*256
      integer nlayers,direction
      double precision origin,xstart,xend
      integer nbands,iband,nkpoints,ikpoint,nele,iline,idum
      integer iwrite,iread,nspecies,ispecies,natoms,iatom
      integer numener,iener,jener,ilayer,kener
      integer ispin,iorb
      logical spinpol,has_G,has_X,has_R,has_M,noncollinear
      double precision, allocatable :: eigenvalues_up(:,:)
      double precision, allocatable :: eigenvalues_down(:,:)
      double precision, allocatable :: energies(:),delta_E
      double precision, allocatable :: dos_tot(:,:)
      double precision, allocatable :: dos_lp(:,:)
      double precision, allocatable :: dos_lp_proj(:,:,:,:)
      double precision, allocatable :: dos_proj(:,:,:,:)
!      double precision, allocatable :: proj_up(:,:,:,:)
!      double precision, allocatable :: proj_ele_up(:,:,:,:)
!      double precision, allocatable :: proj_down(:,:,:,:)
!      double precision, allocatable :: proj_ele_down(:,:,:,:)
      double precision, allocatable :: k(:,:),svec(:)
      double precision s,tol,efermi,vecs(1:3,1:3),fdum,jdos_tot(1:2)
      double precision tot_dos_lp_proj_1,tot_dos_lp_proj_2,dos_tol
      logical have_vbm,have_cbm
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      !
      print fsubstart,'get_vasp_projected_dos'  
      print '(8x,"Using as DOS tolerance ",F12.6)',dos_tol
      !
      ! begin initialize
      !
      eigenval="EIGENVAL"
      procar="PROCAR"
      outcar="OUTCAR"
      doscar="DOSCAR"
      s=0.0d0
      spinpol=.false.
      noncollinear=.false.
      ispin=1
      has_G=.false.
      has_X=.false.
      has_R=.false.
      has_M=.false.
      tol=1D-4
      !
      ! end initialize
      !
      !
      ! begin read atom info
      !
      call read_coords(outcar,"outcar",atoms,natoms,species,              &
     &                 nspecies,vecs)
      if (talk) print '(8x,"atom info read")'
      !
      ! end read atom info
      !
      ! begin read Fermi energy and spin polarization
      !
      open(51,file=outcar,status='old')
10    read(51,'(A256)',err=100,end=20) line
      if(index(line,'E-fermi').gt.0) then
        iread=index(line,'E-fermi')+9
        read(line(iread:),*) efermi
      end if
      if(index(line,'   ISPIN  =    ').gt.0) then
        iread=index(line,'=')+1
        read(line(iread:),*) ispin
        if (ispin==2) spinpol=.true.
      end if
      if(index(line,' LNONCOLLINEAR ').gt.0) then
        iread=index(line,'=')+1
        read(line(iread:),*) noncollinear
      end if
      goto 10
      !
20    continue
      close(51)
      !
      ! end read Fermi energy and spin polarization
      !
      ! begin read dos
      !
      open(51,file=doscar,status="old", err=101)
      rewind(51)
      do iline=1,5
        read(51,*) line
      end do
      read(51,*) fdum,fdum,numener
      print '(8x,"Using ",I0," layers in direction ",I0,", starting at "  &
     &        ,F10.6)',nlayers, direction, origin
      print '(8x,I0," energies in DOS")' , numener
      print '(8x,I0," spin channels")' , ispin
      allocate(dos_tot(numener,1:2))
      ! allocate projected dos (atoms, energies, total s px py pz d1 d2 d3 d4 d5 , spin)
      allocate(dos_proj(1:natoms,1:numener,0:9,1:2))
      allocate(energies(numener))
      ! red total dos:
      do iener=1,numener
        read(51,*) energies(iener),dos_tot(iener,1:ispin)
        energies(iener)=energies(iener)-efermi
      end do
      ! read atom-resolved DOSes:
      do iatom=1,natoms
       read(51,*) line
        do iener=1,numener
          !read(51,*) fdum, dos_proj(iatom, iener, 1:9,1:ispin)
          if (.not.noncollinear) then 
            read(51,*) fdum, (dos_proj(iatom, iener,iorb,1:ispin),        &
     &                      iorb=1,9)
          else
            read(51,*) fdum, (dos_proj(iatom, iener,iorb,1:ispin), fdum,  &
     &           fdum,fdum,iorb=1,9)
          end if
          dos_proj(iatom,iener,0,1)=sum(dos_proj(iatom,iener,1:9,1),1)
          if (ispin ==2 ) then
            dos_proj(iatom,iener,0,2)=sum(dos_proj(iatom,iener,1:9,2),1)
          end if
        end do
      end do
      close(51)
      !
      ! end read dos
      !
      ! begin print DOS
      !
      ! total DOS :
      filename="DOS_TOT.DAT"
      open(51,file=filename, status="replace")
      do iener=1,numener
        write(51,'(12(F12.6))') energies(iener),                          &
     &      dos_tot(iener,1:ispin)
      end do
      close(51)
      !
      ! joint DOS :
      filename="JDOS_TOT.DAT"
      open(51,file=filename, status="replace")
      delta_E=(energies(numener)-energies(1))/(dble(numener-1))
      do iener=1,numener
        jdos_tot=0.0d0
        do jener=1,numener-iener+1
          kener=jener+iener-1
          jdos_tot(1:ispin)=jdos_tot(1:ispin)+dos_tot(jener,1:ispin)      &
     &      *theta_function(-energies(jener))                             &
     &      *dos_tot(kener,1:ispin)                                       &
     &      *theta_function(energies(kener))*delta_E                
        end do! jener
        write(51,'(3(F20.6))') energies(iener)-energies(1),               &
     &    jdos_tot(1:ispin)
      end do ! iener
      close(51)
      !
      !
      ! atom-projected DOS :
      do iatom=1,natoms
        filename="DOS_PROJ_UP_"
        iwrite=len_trim(filename)+1
        write(filename(iwrite:),'(I0)') iatom           
        iwrite=len_trim(filename)+1
        write(filename(iwrite:),'(A1,2A)')  '_',                          &
     &    trim(adjustl(atoms(iatom)%name ) )
        iwrite=len_trim(filename)+1
        write(filename(iwrite:),'(".DAT")')
        open(51,file=filename, status="replace")
        write(51,'("# E, tot, s, py, pz, px, dxy, dyz, dz2, dxz, dx2")')
        do iener=1,numener
          write(51,'(11(F12.6))') energies(iener),                        &
     &        dos_proj(iatom,iener,0:9,1)
        end do
        close(51)
        if (ispin==2) then
          filename="DOS_PROJ_DOWN_"
          iwrite=len_trim(filename)+1
          write(filename(iwrite:),'(I0)')  iatom           
          iwrite=len_trim(filename)+1
          write(filename(iwrite:),'(A1,2A)')  '_',                        &
     &      trim(adjustl(atoms(iatom)%name ) )
          iwrite=len_trim(filename)+1
          write(filename(iwrite:),'(".DAT")')
          open(51,file=filename, status="replace")
          write(51,'("# E, tot, s, py, pz, px, dxy, dyz, dz2, dxz, dx2")  &
     &')
          do iener=1,numener
            write(51,'(11(F12.6))') energies(iener),                        &
     &          dos_proj(iatom,iener,0:9,2)
          end do
          close(51)
        end if
      end do
      !
      !
      ! layer-projected DOS :
      allocate(dos_lp(numener,2))
      allocate(dos_lp_proj(natoms,numener,0:9,2))
      do ilayer=1,nlayers
        dos_lp=0.0d0
        dos_lp_proj=0.0d0
        xstart=origin+dble(ilayer-1)/dble(nlayers)
        xend=origin+dble(ilayer)/dble(nlayers)
        filename="DOS_PROJ_UP_LAYER_"
        filename2="DOS_PROJ_UP_LAYER_"
        iwrite=len_trim(filename)+1
        write(filename(iwrite:),'(I0)') ilayer           
        write(filename2(iwrite:),'(I0)') ilayer           
        iwrite=len_trim(filename)+1
        write(filename(iwrite:),'(".DAT")')
        write(filename2(iwrite:),'("_atoms.DAT")')
        open(51,file=filename, status="replace")
        open(52,file=filename2, status="replace")
        write(51,'("# E, tot, s, py, pz, px, dxy, dyz, dz2, dxz, dx2")')
        do iatom=1,natoms
          do idum=-1,1
          if (atoms(iatom)%where(direction)+float(idum).ge.xstart.and.    &
     &        atoms(iatom)%where(direction)+float(idum).lt.xend) then
            dos_lp(:,1:ispin)=dos_lp(:,1:ispin)+dos_proj(iatom,:,         &
     &           0,1:ispin)
            dos_lp_proj(iatom,:,0:9,1:ispin)=dos_lp_proj(iatom,:,0:9,     &
     &           1:ispin) +dos_proj(iatom,:,0:9,1:ispin)
            write(52,'(I5,1x,A2)') iatom,atoms(iatom)%name(1:2)   
          end if
          end do ! idum
        end do ! iatom
        close(52)
        have_cbm=.false.
        do iener=1,numener
!          write(51,'(2(F12.6))') energies(iener),                        &
!     &        dos_lp(iener,1)
          write(51,'(12(F12.6))') energies(iener),                        &
     &        sum(dos_lp_proj(1:natoms,iener,0,1)),                       &
     &        sum(dos_lp_proj(1:natoms,iener,1,1)),                       &
     &        sum(dos_lp_proj(1:natoms,iener,2,1)),                       &
     &        sum(dos_lp_proj(1:natoms,iener,3,1)),                       &
     &        sum(dos_lp_proj(1:natoms,iener,4,1)),                       &
     &        sum(dos_lp_proj(1:natoms,iener,5,1)),                       &
     &        sum(dos_lp_proj(1:natoms,iener,6,1)),                       &
     &        sum(dos_lp_proj(1:natoms,iener,7,1)),                       &
     &        sum(dos_lp_proj(1:natoms,iener,8,1)),                       &
     &        sum(dos_lp_proj(1:natoms,iener,9,1))
            !
            ! begin check for CBM
            !
            if (iener.lt.numener) then
                tot_dos_lp_proj_1=sum(dos_lp_proj(1:natoms,iener,0,1))    &
                tot_dos_lp_proj_2=sum(dos_lp_proj(1:natoms,iener+1,0,1))  &
               !
               ! begin check for CBM
               !
               if       (tot_dos_lp_proj_1.le.dos_tol.and.                &
     &                    tot_dos_lp_proj_2.gt.dos_tol)   then    
                 if (.not.have_cbm) then
                   if (energies(iener).gt.0.01d0) then ! if above Fermi level
                      print '(8x,"Layer ",I0,", CBM up", F12.6)',         &
     &                       ilayer,energies(iener+1)
                      have_cbm=.true.
                   end if ! energy above EF
                 end if ! not have cbm
               end if ! DOS crosses threshold
               !
            end if ! iener.lt.numener
            !
            ! end check for CBM
            !
        end do ! iener
        close(51)
        !
        ! begin check for VBM
        !
        have_vbm=.false.
        do iener=numener-1,1,-1
            tot_dos_lp_proj_1=sum(dos_lp_proj(1:natoms,iener,0,1))        &
            tot_dos_lp_proj_2=sum(dos_lp_proj(1:natoms,iener+1,0,1))      &
           !
           ! begin check for VBM
           !
           if       (tot_dos_lp_proj_1.gt.dos_tol.and.                    &
     &                tot_dos_lp_proj_2.le.dos_tol)   then    
             if (.not.have_vbm) then
               if (energies(iener).lt.0.01d0) then ! if below Fermi level
                  print '(8x,"Layer ",I0,", VBM up", F12.6)',             &
     &                   ilayer,energies(iener+1)
                  have_vbm=.true.
               end if
             end if
           end if
           !
        end do ! iener
        !
        ! end check for VBM
        !
        if (ispin==2) then
          filename="DOS_PROJ_DOWN_LAYER_"
          iwrite=len_trim(filename)+1
          write(filename(iwrite:),'(I0)') ilayer           
          iwrite=len_trim(filename)+1
          write(filename(iwrite:),'(".DAT")')
          open(51,file=filename, status="replace")
          !do iatom=1,natoms
          !  if (atoms(iatom)%where(direction).ge.xstart.and.                &
     &    !      atoms(iatom)%where(direction).lt.xend) then
          !    dos_lp(:,1:ispin)=dos_lp(:,1:ispin)+dos_proj(iatom,:,         &
     &    !         0,1:ispin)
          !    dos_lp_proj(iatom,:,0:9,1:ispin)=dos_lp_proj(iatom,:,0:9,     &
     &    !         1:ispin) +dos_proj(iatom,:,0:9,1:ispin)
          !  end if
          !end do ! iatom
          have_cbm=.false.
          write(51,'("# E, tot, s, py, pz, px, dxy, dyz, dz2, dxz, dx2")  &
     &')
          do iener=1,numener
!            write(51,'(2(F12.6))') energies(iener),                        &
!     &          dos_lp(iener,1)
            write(51,'(12(F12.6))') energies(iener),                        &
     &          sum(dos_lp_proj(1:natoms,iener,0,2)),                       &
     &          sum(dos_lp_proj(1:natoms,iener,1,2)),                       &
     &          sum(dos_lp_proj(1:natoms,iener,2,2)),                       &
     &          sum(dos_lp_proj(1:natoms,iener,3,2)),                       &
     &          sum(dos_lp_proj(1:natoms,iener,4,2)),                       &
     &          sum(dos_lp_proj(1:natoms,iener,5,2)),                       &
     &          sum(dos_lp_proj(1:natoms,iener,6,2)),                       &
     &          sum(dos_lp_proj(1:natoms,iener,7,2)),                       &
     &          sum(dos_lp_proj(1:natoms,iener,8,2)),                       &
     &          sum(dos_lp_proj(1:natoms,iener,9,2))
            !
            ! begin check for CBM
            !
            if (iener.lt.numener) then
                tot_dos_lp_proj_1=sum(dos_lp_proj(1:natoms,iener,0,2))    &
                tot_dos_lp_proj_2=sum(dos_lp_proj(1:natoms,iener+1,0,2))  &
               !
               ! begin check for CBM
               !
               if       (tot_dos_lp_proj_1.le.dos_tol.and.                &
     &                    tot_dos_lp_proj_2.gt.dos_tol)   then    
                 if (.not.have_cbm) then
                   if (energies(iener).gt.0.01d0) then ! if above Fermi level
                      print '(8x,"Layer ",I0,", CBM down", F12.6)',       &
     &                       ilayer,energies(iener+1)
                      have_cbm=.true.
                   end if
                 end if
               end if
               !
               ! end check for CBM
               !
            end if ! iener.lt.numener
            !
            ! end check for CBM
            !
          end do ! iener
          close(51)
          !
          ! begin check for VBM
          !
          have_vbm=.false.
          do iener=numener-1,1,-1
                tot_dos_lp_proj_1=sum(dos_lp_proj(1:natoms,iener,0,2))    &
                tot_dos_lp_proj_2=sum(dos_lp_proj(1:natoms,iener+1,0,2))  &
                if       (tot_dos_lp_proj_2.le.dos_tol.and.               &
     &                    tot_dos_lp_proj_1.gt.dos_tol)  then     
                  if (energies(iener).lt.0.01d0) then ! if below Fermi level
                    if (.not.have_vbm) then
                      print '(8x,"Layer ",I0,", VBM down", F12.6)',       &
     &                       ilayer,energies(iener)
                      have_vbm=.true.
                    end if ! .not.have_vbm
                  end if ! if close to Fermi level
                end if       
            !
            ! end check for VBM, CBM
            !
          end do ! iener
        end if ! ispin==2
      end do ! ilayer
      !
      ! end print DOS
      !
      ! begin print gnuplot file
      !
      ! total DOS:
      open(52,file='plot_dos.gpi',status='replace')
      write(52,'(256A)') '#set term x11 "Times-Roman"' 
      write(52,'(256A)') 'set term post enhanced color "Times-Roman"' 
      write(52,'(256A)') 'set out "DOS_PROJ.eps"' 
      write(52,'(256A)') "set size 0.25, 0.3" 
      write(52,'(256A)') "set xlabel 'E-E_F (eV)' offset 1"
      write(52,'(256A)') "xmin=-5" 
      write(52,'(256A)') "xmax=5" 
      write(52,'(256A)') "set xrange [xmin:xmax]" 
      write(52,'(256A)') "shift=40" 
      write(52,'(256A)') "fac=10" 
      write(52,'(256A)') "" 
      do ilayer=1,1
        ! get filename for that layer:
        filename="DOS_PROJ_UP_LAYER_"
        iwrite=len_trim(adjustl(filename))+1
        write(filename(iwrite:),'(I0)') ilayer           
        iwrite=len_trim(filename)+1
        write(filename(iwrite:),'(".DAT")')
        string='plot '
        iwrite=len_trim(string)+1
        write(string(iwrite:),*)'"',trim(adjustl(filename)),'"'
        !print*, string
        iwrite=len_trim(string)+1
        write(string(iwrite:),'(" u 1:($2+",I0,"*shift) w l lt 1 lc 1 t   &
     & ",A5)') ilayer,' "",\'           
        !print*, string
        write(52,*) adjustl(trim(adjustl(string)))
        if (ispin==2) then
          ! get filename for that layer:
          filename="DOS_PROJ_DOWN_LAYER_"
          iwrite=len_trim(adjustl(filename))+1
          write(filename(iwrite:),'(I0)') ilayer           
          iwrite=len_trim(filename)+1
          write(filename(iwrite:),'(".DAT")')
          string="     "
          iwrite=len_trim(string)+1
          write(string(iwrite:),*) '"',trim(adjustl(filename)),'"'
          !print*, string
          iwrite=len_trim(string)+1
          write(string(iwrite:),'(" u 1:(-$2+",I0,"*shift) w l lt 1 lc 1  &
     & t ",A5)') ilayer,' "",\'           
          !print*, string
          write(52,*) adjustl(trim(adjustl(string)))
        end if
      end do ! ilayer
!      do ilayer=2,nlayers-1
      do ilayer=2,nlayers
        ! get filename for that layer:
        filename="DOS_PROJ_UP_LAYER_"
        iwrite=len_trim(filename)+1
        write(filename(iwrite:),'(I0)') ilayer           
        iwrite=len_trim(filename)+1
        write(filename(iwrite:),'(".DAT")')
        string="     "
        iwrite=len_trim(string)+1
        write(string(iwrite:),*)'"',adjustl(trim(filename)),'"'
        !print*, string
        iwrite=len_trim(string)+1
        write(string(iwrite:),'(" u 1:($2+",I0,"*shift) w l lt 1 lc 1     &
     & t ",A5)') ilayer,' "",\'           
        !print*, string
        write(52,*) adjustl(trim(adjustl(string)))
        if (ispin==2) then
          ! get filename for that layer:
          filename="DOS_PROJ_DOWN_LAYER_"
          iwrite=len_trim(filename)+1
          write(filename(iwrite:),'(I0)') ilayer           
          iwrite=len_trim(filename)+1
          write(filename(iwrite:),'(".DAT")')
          string="     "
          iwrite=len_trim(string)+1
          write(string(iwrite:),*)'"',adjustl(trim(filename)),'"'
          !print*, string
          iwrite=len_trim(string)+1
          write(string(iwrite:),'(" u 1:(-$2+",I0,"*shift) w l lt 1 lc 1  &
     & t ",A5)') ilayer,' "",\'           
          !print*, string
          write(52,*) adjustl(trim(adjustl(string)))
        end if
      end do ! ilayer
!      do ilayer=nlayers,nlayers
!        ! get filename for that layer:
!        filename="DOS_PROJ_UP_LAYER_"
!        iwrite=len_trim(filename)+1
!        write(filename(iwrite:),'(I0)'), ilayer           
!        iwrite=len_trim(filename)+1
!        write(filename(iwrite:),'(".DAT")')
!        string="     "
!        iwrite=len_trim(string)+1
!        write(string(iwrite:),*),'"',adjustl(trim(filename)),'"'
!        !print*, string
!        iwrite=len_trim(string)+1
!        if (ispin==1) then
!          write(string(iwrite:),'(" u 1:($2+",I0,"*shift) w l t ",A5)')   &
!     &        ilayer,' ""'           
!        else
!          write(string(iwrite:),'(" u 1:($2+",I0,"*shift) w l t ",A5)')   &
!     &        ilayer,' "",\'
!        end if
!        !print*, string
!        write(52,*) adjustl(trim(adjustl(string)))
!        if (ispin==2) then
!          ! get filename for that layer:
!          filename="DOS_PROJ_DOWN_LAYER_"
!          iwrite=len_trim(filename)+1
!          write(filename(iwrite:),'(I0)'), ilayer           
!          iwrite=len_trim(filename)+1
!          write(filename(iwrite:),'(".DAT")')
!          string="     "
!          iwrite=len_trim(string)+1
!          write(string(iwrite:),*),'"',adjustl(trim(filename)),'"'
!          !print*, string
!          iwrite=len_trim(string)+1
!          write(string(iwrite:),'(" u 1:(-$2+",I0,"*shift) w l t ",A5)')  &
!     &          ilayer,' ""'           
!          !print*, string
!          write(52,*) adjustl(trim(adjustl(string)))
!        end if
!      end do ! ilayer
        ! get filename for that layer:
        filename="DOS_TOT.DAT"
        string="     "
        iwrite=len_trim(string)+1
        write(string(iwrite:),*)'"',adjustl(trim(filename)),'"'
        !print*, string
        iwrite=len_trim(string)+1
        if (ispin==1) then
          write(string(iwrite:),'(" u 1:($2*fac+",I0,"*shift) w l lt 1    &
     &lc 0 t ",A5)') ilayer,' ""'           
        else
          write(string(iwrite:),'(" u 1:($2*fac+",I0,"*shift) w l lt 1    &
     &lc 0 t ",A5)') ilayer,' "",\'
        end if
        !print*, string
        write(52,*) adjustl(trim(adjustl(string)))
        if (ispin==2) then
          filename="DOS_TOT.DAT"
          string="     "
          iwrite=len_trim(string)+1
          write(string(iwrite:),*)'"',adjustl(trim(filename)),'"'
          !print*, string
          iwrite=len_trim(string)+1
          write(string(iwrite:),'(" u 1:(-$3*fac+",I0,"*shift) w l lt 1   &
     & lc 0 t ",A5)') ilayer,' ""'           
          !print*, string
          write(52,*) adjustl(trim(adjustl(string)))
        end if
      write(52,'(256A)') ' ' 
      write(52,'(256A)') 'system("epstopdf DOS_PROJ.eps")' 
      close(52)
      ! d DOS:
      open(52,file='plot_dos_d.gpi',status='replace')
      write(52,'(256A)') '#set term x11 "Times-Roman"' 
      write(52,'(256A)') 'set term post enhanced color "Times-Roman"' 
      write(52,'(256A)') 'set out "DOS_PROJ_d.eps"' 
      write(52,'(256A)') "set size 0.25, 0.3" 
      write(52,'(256A)') "set xlabel 'E-E_F (eV)' offset 1"
      write(52,'(256A)') "xmin=-3" 
      write(52,'(256A)') "xmax=3" 
      write(52,'(256A)') "set xrange [xmin:xmax]" 
      write(52,'(256A)') "shift=0.2" 
      write(52,'(256A)') "" 
      do ilayer=1,1
        ! get filename for that layer:
        filename="DOS_PROJ_UP_LAYER_"
        iwrite=len_trim(adjustl(filename))+1
        write(filename(iwrite:),'(I0)') ilayer           
        iwrite=len_trim(filename)+1
        write(filename(iwrite:),'(".DAT")')
        string='plot '
        iwrite=len_trim(string)+1
        write(string(iwrite:),*)'"',trim(adjustl(filename)),'"'
        !print*, string
        iwrite=len_trim(string)+1
        write(string(iwrite:),                                            &
     &          '(" u 1:($7+$8+$9+$10+$11+",I0,"*shift) w l t ",A5)')     &
     &        ilayer,' "",\'           
        !print*, string
        write(52,*) adjustl(trim(adjustl(string)))
        if (ispin==2) then
          ! get filename for that layer:
          filename="DOS_PROJ_DOWN_LAYER_"
          iwrite=len_trim(adjustl(filename))+1
          write(filename(iwrite:),'(I0)') ilayer           
          iwrite=len_trim(filename)+1
          write(filename(iwrite:),'(".DAT")')
          string="     "
          iwrite=len_trim(string)+1
          write(string(iwrite:),*)'"',trim(adjustl(filename)),'"'
          !print*, string
          iwrite=len_trim(string)+1
          write(string(iwrite:),                                          &
     &            '(" u 1:(-$7-$8-$9-$10-$11+",I0,"*shift) w l t ",A5)')  &
     &          ilayer,' "",\'           
          !print*, string
          write(52,*) adjustl(trim(adjustl(string)))
        end if
      end do ! ilayer
      do ilayer=2,nlayers-1
        ! get filename for that layer:
        filename="DOS_PROJ_UP_LAYER_"
        iwrite=len_trim(filename)+1
        write(filename(iwrite:),'(I0)') ilayer           
        iwrite=len_trim(filename)+1
        write(filename(iwrite:),'(".DAT")')
        string="     "
        iwrite=len_trim(string)+1
        write(string(iwrite:),*)'"',adjustl(trim(filename)),'"'
        !print*, string
        iwrite=len_trim(string)+1
        write(string(iwrite:),                                            &
     &                  '(" u 1:($7+$8+$9+",I0,"*shift) w l t ",A5)')     &
     &        ilayer,' "",\'           
        !print*, string
        write(52,*) adjustl(trim(adjustl(string)))
        if (ispin==2) then
          ! get filename for that layer:
          filename="DOS_PROJ_DOWN_LAYER_"
          iwrite=len_trim(filename)+1
          write(filename(iwrite:),'(I0)') ilayer           
          iwrite=len_trim(filename)+1
          write(filename(iwrite:),'(".DAT")')
          string="     "
          iwrite=len_trim(string)+1
          write(string(iwrite:),*)'"',adjustl(trim(filename)),'"'
          !print*, string
          iwrite=len_trim(string)+1
          write(string(iwrite:),                                          &
     &            '(" u 1:(-$7-$8-$9-$10-$11+",I0,"*shift) w l t ",A5)')  &
     &          ilayer,' "",\'           
          !print*, string
          write(52,*) adjustl(trim(adjustl(string)))
        end if
      end do ! ilayer
      do ilayer=nlayers,nlayers
        ! get filename for that layer:
        filename="DOS_PROJ_UP_LAYER_"
        iwrite=len_trim(filename)+1
        write(filename(iwrite:),'(I0)') ilayer           
        iwrite=len_trim(filename)+1
        write(filename(iwrite:),'(".DAT")')
        string="     "
        iwrite=len_trim(string)+1
        write(string(iwrite:),*)'"',adjustl(trim(filename)),'"'
        !print*, string
        iwrite=len_trim(string)+1
        if (ispin==1) then
          write(string(iwrite:),                                          &
     &            '(" u 1:($7+$8+$9+$10+$11+",I0,"*shift) w l t ",A5)')   &
     &          ilayer,' ""'           
        else
          write(string(iwrite:),                                          &
     &            '(" u 1:($7+$8+$9+$10+$11+",I0,"*shift) w l t ",A5)')   &
     &          ilayer,' "",\'           
        end if
        !print*, string
        write(52,*) adjustl(trim(adjustl(string)))
        if (ispin==2) then
          ! get filename for that layer:
          filename="DOS_PROJ_DOWN_LAYER_"
          iwrite=len_trim(filename)+1
          write(filename(iwrite:),'(I0)') ilayer           
          iwrite=len_trim(filename)+1
          write(filename(iwrite:),'(".DAT")')
          string="     "
          iwrite=len_trim(string)+1
          write(string(iwrite:),*)'"',adjustl(trim(filename)),'"'
          !print*, string
          iwrite=len_trim(string)+1
          write(string(iwrite:),                                          &
     &            '(" u 1:(-$7-$8-$9-$10-$11+",I0,"*shift) w l t ",A5)')  &
     &          ilayer,' ""'           
          !print*, string
          write(52,*) adjustl(trim(adjustl(string)))
        end if
      end do ! ilayer
      write(52,'(256A)') ' ' 
      write(52,'(256A)') 'system("epstopdf DOS_PROJ_d.eps")' 
      close(52)
      !
      !
      ! end print gnuplot file
      !
      !
!      open(51,file=eigenval, status='old')
!      !
!      ! begin read nbands, nkpoints
!      !    
!      read(51,*,end=100,err=100) idum,idum,idum, idum  
!      if (idum==2) spinpol=.true.   
!      do iline=1,4
!        read(51,'(A1024)',end=100,err=100) line    
!      end do
!      read(51,*,end=100,err=100) nele,nkpoints,nbands
!      !print '(8x,3(1x,I0))',nele,nkpoints,nbands
!      allocate(k(1:nkpoints,1:3),eigenvalues_up(nkpoints,nbands))
!      allocate(proj_up(nkpoints,nbands,natoms,1:10))
!      allocate(proj_ele_up(nkpoints,nbands,nspecies,1:10))
!      proj_ele_up(1:nkpoints,1:nbands,1:nspecies,1:10)=0.0d0
!      allocate(svec(1:nkpoints))
!      svec=0.0d0
!      if (spinpol) then
!        allocate(eigenvalues_down(nkpoints,nbands))
!        allocate(proj_down(nkpoints,nbands,natoms,1:10))
!        allocate(proj_ele_down(nkpoints,nbands,nspecies,1:10))
!        proj_ele_down(1:nkpoints,1:nbands,1:nspecies,1:10)=0.0d0
!      end if
!      !
!      ! end read nbands, nkpoints
!      !
!      !
!      ! begin read kpoints and eigenvalues
!      !    
!      do ikpoint=1,nkpoints
!        read(51,'(A1024)',end=100,err=100) line    
!        read(51,*,end=100,err=100) k(ikpoint,1:3)
!        do iband=1,nbands
!          if (.not.spinpol) then
!            read(51,*,err=100,end=100) idum,                              &
!     &         eigenvalues_up(ikpoint,iband)
!          else
!            read(51,*,err=100,end=100) idum,                              &
!     &         eigenvalues_up(ikpoint,iband),                             &
!     &         eigenvalues_down(ikpoint,iband)
!          end if
!        end do
!      end do
!      close(51)
!      !
!      ! end read kpoints and eigenvalues
!      !
!      !
!      ! begin print eigenvalues
!      !    
!      open(51,file='BS_UP.DAT',status='replace')
!      do ikpoint=1,nkpoints
!        if (ikpoint.gt.1) s=s+norm2(k(ikpoint,1:3)-k(ikpoint-1,1:3))
!        if (ikpoint.gt.1) svec(ikpoint)=s
!        write(51,*) svec(ikpoint),k(ikpoint,1:3),                         &
!     & eigenvalues_up(ikpoint,1:nbands)-efermi
!      end do
!      close(51)
!      if (spinpol) then
!        open(51,file='BS_DOWN.DAT',status='replace')
!        do ikpoint=1,nkpoints
!          write(51,*) svec(ikpoint),k(ikpoint,1:3),                       &
!     &   eigenvalues_down(ikpoint,1:nbands)-efermi
!        end do
!        close(51)
!      end if
!      !
!      ! end print eigenvalues
!      !    
!      ! 
      ! begin print gnuplot file
      !
!      open(51,file='DOS.gpi',status='replace')
!      write(51,'(256A)') 'set term x11 "Times-Roman"' 
!      write(51,'(256A)') '#set term post enhanced color "Times-Roman"' 
!      write(51,'(256A)') '#set out "DOS.eps"' 
      !write(52,'(256A)') '#set out "DOS_PROJ.eps"' 
!      write(51,'(256A)') "set size 0.25, 0.5" 
!      write(51,'(256A)') "set xlabel 'E-E_F (eV)' offset 1"
!      write(51,'(256A)') "xmin=-10" 
!      write(51,'(256A)') "xmax=10" 
!      write(51,'(256A)') "set xrange [ymin:ymax]" 
!      write(52,'(256A)') "set palette negative" 
      ! 
!      ! begin get special points
!      !
!      xtics='('
!      do ikpoint=1,nkpoints
!        if (norm2(k(ikpoint,1:3)-(/0.0,0.0,0.0/)).lt.tol)  then         
!          has_G=.true.
!          iwrite=len_trim(xtics)+1
!          write(xtics(iwrite:),                                           &
!     &     '(A15,1x,F10.6,A1)'), '"{/Symbol G}" ',svec(ikpoint),','  
!          write(51,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
!     &      ",ymin to ",F10.6,",ymax")'),svec(ikpoint),svec(ikpoint) 
!          write(52,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
!     &      ",ymin to ",F10.6,",ymax")'),svec(ikpoint),svec(ikpoint) 
!        end if
!        if (norm2(k(ikpoint,1:3)-(/0.0,0.0,0.5/)).lt.tol)  then         
!          has_X=.true.
!          iwrite=len_trim(xtics)+1
!          write(xtics(iwrite:),                                           &
!     &     '(A8,1x,F10.6,A1)'), '"X" ',svec(ikpoint),','  
!          write(51,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
!     &      ",ymin to ",F10.6,",ymax")'),svec(ikpoint),svec(ikpoint) 
!          write(52,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
!     &      ",ymin to ",F10.6,",ymax")'),svec(ikpoint),svec(ikpoint) 
!        end if
!        if (norm2(k(ikpoint,1:3)-(/0.0,0.5,0.0/)).lt.tol)  then         
!          has_X=.true.
!          iwrite=len_trim(xtics)+1
!          write(xtics(iwrite:),                                           &
!     &     '(A8,1x,F10.6,A1)'), '"X" ',svec(ikpoint),','  
!          write(51,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
!     &      ",ymin to ",F10.6,",ymax")'),svec(ikpoint),svec(ikpoint) 
!          write(52,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
!     &      ",ymin to ",F10.6,",ymax")'),svec(ikpoint),svec(ikpoint) 
!        end if
!        if (norm2(k(ikpoint,1:3)-(/0.5,0.0,0.0/)).lt.tol)  then         
!          has_X=.true.
!          iwrite=len_trim(xtics)+1
!          write(xtics(iwrite:),                                           &
!     &     '(A8,1x,F10.6,A1)'), '"X" ',svec(ikpoint),','  
!          write(51,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
!     &      ",ymin to ",F10.6,",ymax")'),svec(ikpoint),svec(ikpoint) 
!          write(5,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
!     &      ",ymin to ",F10.6,",ymax")'),svec(ikpoint),svec(ikpoint) 
!        end if
!        if (norm2(k(ikpoint,1:3)-(/0.0,0.5,0.5/)).lt.tol)  then         
!          has_M=.true.
!          iwrite=len_trim(xtics)+1
!          write(xtics(iwrite:),                                           &
!     &     '(A8,1x,F10.6,A1)'), '"M" ',svec(ikpoint),','  
!          write(51,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
!     &      ",ymin to ",F10.6,",ymax")'),svec(ikpoint),svec(ikpoint) 
!          write(52,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
!     &      ",ymin to ",F10.6,",ymax")'),svec(ikpoint),svec(ikpoint) 
!        end if
!        if (norm2(k(ikpoint,1:3)-(/0.5,0.0,0.5/)).lt.tol)  then         
!          has_M=.true.
!          iwrite=len_trim(xtics)+1
!          write(xtics(iwrite:),                                           &
!     &     '(A8,1x,F10.6,A1)'), '"M" ',svec(ikpoint),','  
!          write(51,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
!     &      ",ymin to ",F10.6,",ymax")'),svec(ikpoint),svec(ikpoint) 
!          write(52,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
!     &      ",ymin to ",F10.6,",ymax")'),svec(ikpoint),svec(ikpoint) 
!        end if
!        if (norm2(k(ikpoint,1:3)-(/0.5,0.5,0.0/)).lt.tol)  then         
!          has_M=.true.
!          iwrite=len_trim(xtics)+1
!          write(xtics(iwrite:),                                           &
!     &     '(A8,1x,F10.6,A1)'), '"M" ',svec(ikpoint),','  
!          write(51,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
!     &      ",ymin to ",F10.6,",ymax")'),svec(ikpoint),svec(ikpoint) 
!          write(52,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
!     &      ",ymin to ",F10.6,",ymax")'),svec(ikpoint),svec(ikpoint) 
!        end if
!        if (norm2(k(ikpoint,1:3)-(/0.5,0.5,0.5/)).lt.tol)  then         
!          has_R=.true.
!          iwrite=len_trim(xtics)+1
!          write(xtics(iwrite:),                                           &
!     &     '(A8,1x,F10.6,A1)'), '"R" ',svec(ikpoint),','  
!          write(51,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
!     &      ",ymin to ",F10.6,",ymax")'),svec(ikpoint),svec(ikpoint) 
!          write(52,'("set arrow nohead lt 1 lc 0 from ",F10.6,            &
!     &      ",ymin to ",F10.6,",ymax")'),svec(ikpoint),svec(ikpoint) 
!        end if
!      end do
!      iwrite=len_trim(xtics)
!      write(xtics(iwrite:),'(A1)'), ')'  
!      !print*, xtics
!      ! 
!      ! end get special points
!      !
!      !
!      write(51,'(256A)') "set xtics ",xtics 
!      write(52,'(256A)') "set xtics ",xtics 
!      write(51,'(256A)') "N=`awk 'NR==1 {print NF}' BS_UP.DAT`"
!      write(51,'(256A)'),"plot for [i=5:N] 'BS_UP.DAT' u 1:i smooth cspl  &
!     &ines w l lt 1 lc 1 lw 2 t ''"
!      close(51)
!      do ispecies=1,nspecies
!        filename="BS_PROJ_UP_"
!        filename2="BS_PROJ_UP_"
!        write(filename(12:),'(2A)'),                                      &
!     &    trim(adjustl(species(ispecies)%name ) )
!        write(filename2(12:),'(2A)'),                                     &
!     &    trim(adjustl(species(ispecies)%name ) )
!        iwrite=len_trim(filename)+1
!        write(filename(iwrite:),'(".DAT")')
!        write(filename2(iwrite:),'(".EPS",1A)'),"'"
!        write(52,'(256A)') "set out '",filename2
!        write(52,'(256A,1A)')"N=`awk 'NR==1 {print NF}' ",filename ,"`"
!        write(52,'(256A)'),"plot for [i=5:N-1:2] '",                      &
!     &trim(adjustl(filename)),"' u 1:i:i+1 with lines palette lw 2 t ''"
!        write(52,'(256A,1A)') "system('epstopdf ",                        &
!     &trim(adjustl(filename2)),")" 
!      end do
!      close(52)
!      !
!      ! end print gnuplot file
!      !
!      !
!      ! begin read projections
!      !
!      open(51,file=procar, status='old')
!      !    
!      read(51,'(A1024)',end=100,err=100) line    
!      !print*,line
!      !
!      do ikpoint=1,nkpoints
!        do iline=1,3
!          read(51,'(A1024)',end=100,err=100) line    
!          !print*,line
!        end do
!        do iband=1,nbands
!          do iline=1,4
!            read(51,'(A1024)',end=100,err=100) line    
!            !print*,line
!          end do
!          do iatom=1,natoms
!            read(51,*,end=100,err=100) idum,                              &
!     &         proj_up(ikpoint,iband,iatom,1:10)
!            !print*,idum,proj_up(ikpoint,iband,iatom,1:10)
!          end do
!          read(51,'(A1024)',end=100,err=100) line    
!!          if (.not.spinpol) then
!!            read(51,*,err=100,end=100) idum,                              &
!!     &         eigenvalues_up(ikpoint,iband)
!!          else
!!            read(51,*,err=100,end=100) idum,                              &
!!     &         eigenvalues_up(ikpoint,iband),                             &
!!     &         eigenvalues_down(ikpoint,iband)
!!          end if
!        end do
!      end do
!      close(51)
!      !
!      ! end read projections
!      !
!      !
!      ! begin print atom-projected DOS 
!      !    
!      do iatom=1,natoms
!        write(filename,'("BS_PROJ_UP_",I0,".DAT")'),iatom
!        open(51,file=filename,status='replace')
!        do ikpoint=1,nkpoints
!          write(51,*) svec(ikpoint),k(ikpoint,1:3),                       &
!     &    (eigenvalues_up(ikpoint,iband)-efermi,                          &
!     &     proj_up(ikpoint,iband,iatom,10),iband=1,nbands)
!        end do
!        close(51)
!      end do
!!      if (spinpol) then
!!        open(51,file='BS_PROJ_DOWN.DAT',status='replace')
!!        do ikpoint=1,nkpoints
!!          write(51,*) svec(ikpoint),k(ikpoint,1:3),                       &
!!     &   eigenvalues_down(ikpoint,1:nbands)-efermi
!!        end do
!!        close(51)
!!      end if
!      !
!      !
!      ! end print atom-projected DOS
!      !    
!      !
!      ! begin print element-projected DOS 
!      !
!      do ispecies=1,nspecies
!        do ikpoint=1,nkpoints
!          do iband=1,nbands
!            do iatom=1,natoms
!              if(atoms(iatom)%name.eq.species(ispecies)%name)            &
!     &        then
!                !print*,"species name=",species(ispecies)%name
!                proj_ele_up(ikpoint,iband,ispecies,1:10)                  &
!     &            =proj_ele_up(ikpoint,iband,ispecies,1:10)               &
!     &            +proj_up(ikpoint,iband,iatom,1:10)
!                !print*,proj_ele_up(ikpoint,iband,ispecies,1:10)
!              end if
!            end do ! iatom
!          end do ! iband
!        end do ! ikpoint
!        filename="BS_PROJ_UP_"
!        write(filename(12:),'(2A)'),                                      &
!     &    trim(adjustl(species(ispecies)%name ) )
!        iwrite=len_trim(filename)+1
!        write(filename(iwrite:),'(".DAT")')
!        open(51,file=filename,status='replace')
!        do ikpoint=1,nkpoints
!          write(51,*) svec(ikpoint),k(ikpoint,1:3),                       &
!     &    (eigenvalues_up(ikpoint,iband)-efermi,                          &
!     &     proj_ele_up(ikpoint,iband,ispecies,10),iband=1,nbands)
!        end do
!        close(51)
!      end do ! ispecies
!      !
!      !
!      ! end print element-projected DOS 
      !
      print fsubendext,'get_vasp_projected_dos'  
      ! 
      return
      !
100   continue    
      close(51)
      nerr=nerr+1
      print ferrmssg,'Problem with PROCAR file'
      !
101   continue    
      close(51)
      nerr=nerr+1
      print ferrmssg,'DOSCAR not found'
      return
      !
      end subroutine get_vasp_projected_dos

c---------------------------------------------------------------------

      subroutine vasp_get_eps_vs_omega(rotmat)
      use defs
      use linalg
      implicit none 
      integer i,j
      character*1024 :: line
      double precision :: eps(1:3,1:3),epsi(1:3,1:3),energy
      double precision, allocatable :: epseigs(:),epsev(:,:)
      double precision, allocatable :: epseigsim(:)
      double precision, optional :: rotmat(1:3,1:3)
      double precision epsrot(3,3)
      !
      print fsubstart,'get_vasp_eps_vs_omega'  
      !
      ! begin initialize
      !
      open(51,file='OUTCAR',status='old', err=100)  
      rewind(51)
      !
      ! begin read eps2
      !
      open(52,file='EPS2.DAT',status='replace', err=102)  
10    read(51,'(A1024)',end=101,err=101) line
      if (index(line,'frequency dependent IMAGINARY DIELECTRIC FUNCTION   &
     &(independent particle, no local field effects)').gt.0) then
        write(52,'("# ",A128)') adjustl(trim(line)) 
        read(51,'(A1024)',err=102,end=102) line
        write(52,'("# ",A128," eigenvalues & eigenvectors")')             &
     &            adjustl(trim(line)) 
        read(51,'(A1024)',err=102,end=102) line
        write(52,'("# ",A128)') adjustl(trim(line)) 
        read(51,'(A1024)',err=102,end=102) line
        read(line(1:),*,err=102) energy, eps(1,1),eps(2,2),               &
     &       eps(3,3),eps(1,2),eps(2,3),eps(3,1)
        eps(2,1)=eps(1,2)
        eps(1,3)=eps(3,1)
        do i=1,3
            do j=1,3
                if (abs(eps(i,j)).lt.1E-10) eps(i,j)=0.0d0
            end do
        end do
        epsi=eps ! epsi is overwritten during diagonalization
        epsrot=eps
        if (allocated(epseigs)) deallocate(epseigs)
        if (allocated(epsev)) deallocate(epsev)
        !call get_mateigs(epsi,epseigs,.false.)
!        call get_mateigs(epsi,epseigs,.true.)
!        write(52,'(A128,3x,3(F12.6))') adjustl(trim(line)),               &
!     &                                 epseigs(1:3) 
        !call get_mateigs_full(epsi,epseigs,epsev)
        call get_mateigs_full_symm(epsi,epseigs,epsev)
        call sort_mateigs(epseigs,epsev)
        if (present(rotmat)) then
          call rotate_matrix(3,eps,rotmat,epsrot)
        end if
!        print '(8x,3(F7.3))',dot_product(epsev(:,1),epsev(:,2)),          &
!     &                       dot_product(epsev(:,2),epsev(:,3)),          &
!     &                       dot_product(epsev(:,3),epsev(:,1))
!        write(52,'(A96,1x,3(F12.6),1x,9(F7.3))')                          &
!     &            adjustl(trim(line)),epseigs(1:3),(epsev(1:3,i),i=1,3) 
        write(52,'(A96,1x,3(F12.6),1x,9(F7.3),5x,9(F12.6))')              &
     &      adjustl(trim(line)),epseigs(1:3),(epsev(1:3,i),i=1,3),        &
     &      (epsrot(i,1:3),i=1,3)           
        do while (index(line,'.').gt.0)
            read(51,'(A1024)',err=101,end=101) line
            if (index(line,'.').gt.0) then
              read(line(1:),*,err=102) energy, eps(1,1),eps(2,2),         &
     &             eps(3,3),eps(1,2),eps(2,3),eps(3,1)
              eps(2,1)=eps(1,2)
              eps(1,3)=eps(3,1)
              do i=1,3
                  do j=1,3
                      if (abs(eps(i,j)).lt.1E-10) eps(i,j)=0.0d0
                  end do
              end do
              epsi=eps ! epsi is overwritten during diagonalization
              epsrot=eps
              if (allocated(epseigs)) deallocate(epseigs)
              if (allocated(epseigsim)) deallocate(epseigsim)
              if (allocated(epsev)) deallocate(epsev)
              !call get_mateigs(epsi,epseigs,.false.)
!              call get_mateigs(epsi,epseigs,.true.)
!              write(52,'(A128,3x,3(F12.6))') adjustl(trim(line)),         &
!     &                                       epseigs(1:3) 
              !call get_mateigs_full(epsi,epseigs,epsev)
              !call get_mateigs_full_cmplx(epsi,epseigs,epseigsim,epsev)
              call get_mateigs_full_symm(epsi,epseigs,epsev)
              call sort_mateigs(epseigs,epsev)
!              print '(8x,6(F7.3))',dot_product(epsev(:,1),epsev(:,2)),    &
!     &                       dot_product(epsev(:,2),epsev(:,3)),          &
!     &          dot_product(epsev(:,3),epsev(:,1)),epseigsim(:)
              if (present(rotmat)) then
                call rotate_matrix(3,eps,rotmat,epsrot)
              end if
!              write(52,'(A96,1x,3(F12.6),1x,9(F7.3))')                    &
!     &            adjustl(trim(line)),epseigs(1:3),(epsev(1:3,i),i=1,3) 
              write(52,'(A96,1x,3(F12.6),1x,9(F7.3),5x,9(F12.6))')        &
     &            adjustl(trim(line)),epseigs(1:3),(epsev(1:3,i),i=1,3),  &
     &            (epsrot(i,1:3),i=1,3)           
            end if
        end do
        goto 11
      end if
      goto 10

11    continue      
      close(52)
      !
      ! end read eps2
      !
      !
      ! begin read eps1
      !
      open(52,file='EPS1.DAT',status='replace', err=103)  
12    read(51,'(A1024)',end=101,err=101) line
      if (index(line,'frequency dependent      REAL DIELECTRIC FUNCTION   &
     &(independent particle, no local field effects)').gt.0) then
        write(52,'("# ",A128)') adjustl(trim(line)) 
        read(51,'(A1024)',err=101,end=101) line
        write(52,'("# ",A128," eigenvalues & eigenvectors")')             &
     &      adjustl(trim(line)) 
        read(51,'(A1024)',err=101,end=101) line
        write(52,'("# ",A128)') adjustl(trim(line)) 
        read(51,'(A1024)',err=101,end=101) line
        read(line(1:),*,err=101) energy, eps(1,1),eps(2,2),               &
     &       eps(3,3),eps(1,2),eps(2,3),eps(3,1)
        eps(2,1)=eps(1,2)
        eps(1,3)=eps(3,1)
        do i=1,3
            do j=1,3
                if (abs(eps(i,j)).lt.1E-10) eps(i,j)=0.0d0
            end do
        end do
        epsi=eps ! epsi is overwritten during diagonalization
        epsrot=eps
        if (allocated(epseigs)) deallocate(epseigs)
        if (allocated(epsev)) deallocate(epsev)
        !call get_mateigs(epsi,epseigs,.false.)
!        call get_mateigs(epsi,epseigs,.true.)
!        write(52,'(A128,3x,3(F12.6))') adjustl(trim(line)),               &
        !call get_mateigs_full(eps,epseigs,epsev)
        call get_mateigs_full_symm(epsi,epseigs,epsev)
        call sort_mateigs(epseigs,epsev)
        if (present(rotmat)) then
          call rotate_matrix(3,eps,rotmat,epsrot)
        end if
!              print '(8x,3(F7.3))',dot_product(epsev(:,1),epsev(:,2)),    &
!     &                       dot_product(epsev(:,2),epsev(:,3)),          &
!     &                       dot_product(epsev(:,3),epsev(:,1))
!        write(52,'(A96,1x,3(F12.6),1x,9(F7.3))')                          &
!     &     adjustl(trim(line)),  epseigs(1:3),                            &
!     &       (epsev(1:3,i),i=1,3) 
         write(52,'(A96,1x,3(F12.6),1x,9(F7.3),5x,9(F12.6))')             &
     &       adjustl(trim(line)),epseigs(1:3),(epsev(1:3,i),i=1,3),       &
     &       (epsrot(i,1:3),i=1,3)           
        do while (index(line,'.').gt.0)
            read(51,'(A1024)',err=101,end=101) line
            if (index(line,'.').gt.0) then
              read(line(1:),*,err=101) energy, eps(1,1),eps(2,2),         &
     &             eps(3,3),eps(1,2),eps(2,3),eps(3,1)
              eps(2,1)=eps(1,2)
              eps(1,3)=eps(3,1)
              do i=1,3
                  do j=1,3
                      if (abs(eps(i,j)).lt.1E-10) eps(i,j)=0.0d0
                  end do
              end do
              epsi=eps ! epsi is overwritten during diagonalization
              epsrot=eps
              if (allocated(epseigs)) deallocate(epseigs)
              if (allocated(epsev)) deallocate(epsev)
!              call get_mateigs(epsi,epseigs,.false.)
!              write(52,'(A128,3x,3(F12.6))') adjustl(trim(line)),         &
!     &                                       epseigs(1:3) 
              !call get_mateigs_full(epsi,epseigs,epsev)
              call get_mateigs_full_symm(epsi,epseigs,epsev)
              call sort_mateigs(epseigs,epsev)
              if (present(rotmat)) then
                call rotate_matrix(3,eps,rotmat,epsrot)
              end if
!              print '(8x,3(F7.3))',dot_product(epsev(:,1),epsev(:,2)),    &
!     &                       dot_product(epsev(:,2),epsev(:,3)),          &
!     &                       dot_product(epsev(:,3),epsev(:,1))
!              write(52,'(A96,1x,3(F12.6),1x,9(F7.3))')                    &
!     &            adjustl(trim(line)),epseigs(1:3),(epsev(1:3,i),i=1,3) 
              write(52,'(A96,1x,3(F12.6),1x,9(F7.3),5x,9(F12.6))')        &
     &            adjustl(trim(line)),epseigs(1:3),(epsev(1:3,i),i=1,3),  &
     &            (epsrot(i,1:3),i=1,3)           
            end if
        end do
        goto 13
      end if
      goto 12

13    continue      
      close(52)
      !
      ! end read eps1
      !
      close(51)
      print fsubendext,'vasp_get_eps_vs_omega'  
      return
      !
100   continue    
      close(51)
      nerr=nerr+1
      print ferrmssg,'OUTCAR not found'
      return
101   continue    
      close(51)
      nerr=nerr+1
      print ferrmssg,'EPS(omega) not found'
      return
102   continue    
      close(51)
      nerr=nerr+1
      print ferrmssg,'EPS2.DAT could not be written'
      return
103   continue    
      close(51)
      nerr=nerr+1
      print ferrmssg,'EPS1.DAT could not be written'
      return
      !
      end subroutine vasp_get_eps_vs_omega

!--------------------------------------------------------------------------

      subroutine vasp_get_eps_vs_omega_xml(rotmat)
      use defs
      use linalg
      implicit none 
      integer i,j,iline,nomega,istart
      character*1024 :: line
      double precision :: eps(1:3,1:3),epsi(1:3,1:3),energy
      double precision, allocatable :: epseigs(:),epsev(:,:)
      double precision, allocatable :: epseigsim(:)
      double precision, optional :: rotmat(1:3,1:3)
      double precision epsrot(3,3)
      !
      print fsubstart,'get_vasp_eps_vs_omega_xml'  
      !
      ! begin initialize
      !
      open(51,file='vasprun.xml',status='old', err=100)  
      rewind(51)
      !
      ! begin read eps2
      !
      open(52,file='EPS2.DAT',status='replace', err=102)  
10    read(51,'(A1024)',end=101,err=101) line
      if (index(line,'NEDOS').gt.0) then
        istart=index(line,'NEDOS')+7
        read(line(istart:istart+5),*) nomega
        print '(8x,I0," frequencies")',nomega
      end if
      if (index(line,'<dielectricfunction>').gt.0) then
        read(51,'(A1024)',err=102,end=102) line
        if (index(line,'imag').le.0) goto 102 ! error handling
        write(52,'("# Dielectric function, imaginary part")') 
        write(52,'("# E (eV), XX, YY,ZZ,XY,YZ,ZX,  eigenvalues & eigenve  &
     &ctors")')  
        do iline=1,10
          read(51,'(A1024)',err=102,end=102) line ! header
        end do
        do iline=1,nomega
          read(51,'(A1024)',err=102,end=102) line 
          read(line(9:85),*,err=102) energy, eps(1,1),eps(2,2),           &
     &       eps(3,3),eps(1,2),eps(2,3),eps(3,1)
          eps(2,1)=eps(1,2)
          eps(1,3)=eps(3,1)
          do i=1,3
              do j=1,3
                  if (abs(eps(i,j)).lt.1E-10) eps(i,j)=0.0d0
              end do
          end do
          epsi=eps ! epsi is overwritten during diagonalization
          epsrot=eps
          if (allocated(epseigs)) deallocate(epseigs)
          if (allocated(epseigsim)) deallocate(epseigsim)
          if (allocated(epsev)) deallocate(epsev)
          call get_mateigs_full_symm(epsi,epseigs,epsev)
          call sort_mateigs(epseigs,epsev)
          if (present(rotmat)) then
            call rotate_matrix(3,eps,rotmat,epsrot)
          end if
          write(52,'(A96,1x,3(F12.6),1x,9(F7.3),5x,9(F12.6))')            &
     &      adjustl(trim(line(9:85))),epseigs(1:3),(epsev(1:3,i),i=1,3),  &
     &      (epsrot(i,1:3),i=1,3)           
        end do
        close(52)
        do iline=1,14
          read(51,'(A1024)',err=102,end=102) line ! header
        end do
        open(52,file='EPS1.DAT',status='replace', err=102)  
        write(52,'("# Dielectric function, real part")') 
        write(52,'("# E (eV), XX, YY,ZZ,XY,YZ,ZX,  eigenvalues & eigenve  &
     &ctors")')  
        do iline=1,nomega
          read(51,'(A1024)',err=102,end=102) line 
          read(line(9:85),*,err=102) energy, eps(1,1),eps(2,2),           &
     &       eps(3,3),eps(1,2),eps(2,3),eps(3,1)
          eps(2,1)=eps(1,2)
          eps(1,3)=eps(3,1)
          do i=1,3
              do j=1,3
                  if (abs(eps(i,j)).lt.1E-10) eps(i,j)=0.0d0
              end do
          end do
          epsi=eps ! epsi is overwritten during diagonalization
          epsrot=eps
          if (allocated(epseigs)) deallocate(epseigs)
          if (allocated(epseigsim)) deallocate(epseigsim)
          if (allocated(epsev)) deallocate(epsev)
          call get_mateigs_full_symm(epsi,epseigs,epsev)
          call sort_mateigs(epseigs,epsev)
          if (present(rotmat)) then
            call rotate_matrix(3,eps,rotmat,epsrot)
          end if
          write(52,'(A96,1x,3(F12.6),1x,9(F7.3),5x,9(F12.6))')            &
     &      adjustl(trim(line(9:85))),epseigs(1:3),(epsev(1:3,i),i=1,3),  &
     &      (epsrot(i,1:3),i=1,3)           
        end do
        goto 11
      end if
      goto 10

11    continue      
      close(52)
      !
      ! end read eps2 and eps1
      !
      close(51)
      print fsubendext,'vasp_get_eps_vs_omega_xml'  
      return
      !
100   continue    
      close(51)
      nerr=nerr+1
      print ferrmssg,'vasprun.xml not found'
      return
101   continue    
      close(51)
      nerr=nerr+1
      print ferrmssg,'eps(omega) not found in file'
      return
102   continue    
      close(51)
      nerr=nerr+1
      print ferrmssg,'EPS2.DAT could not be written'
      return
103   continue    
      close(51)
      nerr=nerr+1
      print ferrmssg,'EPS1.DAT could not be written'
      return
      !
      end subroutine vasp_get_eps_vs_omega_xml

c---------------------------------------------------------------------

      subroutine read_symops_vasp()
      use defs
      use linalg, only : angle_axis_2_rotmatrix,angle_axis_2_mirrmatrix
      implicit none 
      type(symop), allocatable :: symops(:)
      ! local
      character*1024 :: line
      integer nsymop,isymop,idum,i
      double precision matrix(1:3,1:3),axis(1:3),angle,det_mat
      !
      print fsubstart,'read_symops_vasp'  
      !
      ! begin initialize
      !
      open(51,file='OUTCAR',status='old', err=100)  
      rewind(51)
      !
      ! begin read eps2
      !
      open(52,file='SYMOPS_VASP.DAT',status='replace', err=102)  
      !
      ! Begin Find number of symmetry operations
      !
10    read(51,'(A1024)',end=101,err=101) line
      if (index(line,'Subroutine GETGRP returns: Found').gt.0) then
        read(line(index(line,"Found")+5:),*) nsymop
        allocate(symops(nsymop))
        goto 11
      end if
      goto 10
      !
      ! End Find number of symmetry operations
      !
11    continue
      !
      ! Begin read symmetry operations
      !
12    read(51,'(A1024)',end=101,err=101) line
      if (index(line,'Space group operators:').gt.0) then
        isymop=isymop+1
        !
        ! begin read/write header (skip later)
        !
        read(51,'(A1024)',end=101,err=101) line
        write(52,'("# ",A128)') adjustl(trim(line)) 
        !
        ! end read/write header (skip later)
        !
        ! begin read symops
        !
        do isymop=1,nsymop
            read(51,*,err=101,end=101) idum,det_mat,angle,axis(1:3),      &
     &              symops(isymop)%trala(1:3)
             ! print*, norm2(axis)
            !  
            ! begin if det(A)==1, it is a rotation
            !  
            !if(det_mat.gt.0.5d0) then
              call angle_axis_2_rotmatrix(angle,axis,matrix)
            !end if
            !  
            ! end if det(A)==1, it is a rotation
            !  
            ! begin if det(A)==-1, it is a mirror operation
            !
            !if (det_mat.lt.-0.5d0) then
            !  call angle_axis_2_mirrmatrix(angle,axis,matrix)
            !end if
            !  
            ! end if det(A)==-1, it is a mirror operation
            !
            symops(isymop)%mat=matrix
            write(52,'(I0,8(F15.8))',err=102) idum,det_mat,              &
     &              angle,axis(1:3), symops(isymop)%trala(1:3)
!            write(52,'(I0,12(F15.8))',err=102) idum,symops(isymop)%mat,   &
!     &              symops(isymop)%trala(1:3)
        end do
        !
        ! end read symops
        !
        goto 13
      end if
      goto 12
      !
13    continue
      !
      close(51)
      !      
      ! begin transform symops from vasp format to matrix format
      !
      !      
      ! end transform symops from vasp format to matrix format
      !
      close(52)
      open(52,file='SYMOPS.DAT',status='replace', err=102)  
      do isymop=1,nsymop
            write(52,1000,err=102) (symops(isymop)%mat(i,:),              &
     &              symops(isymop)%trala(i),i=1,3)
      end do
      !
      close(52)
      !
      print fsubendext,'read_symops_vasp'  
      return
      !
100   continue    
      close(51)
      nerr=nerr+1
      print ferrmssg,'OUTCAR not found'
      return
101   continue    
      close(51)
      nerr=nerr+1
      print ferrmssg,'symops not found'
      return
102   continue    
      close(51)
      nerr=nerr+1
      print ferrmssg,'SYMOPS_VASP.DAT could not be written'
      return
      !
 1000 format(3(3(F15.8),5x,F15.8,/))
      !
      end subroutine read_symops_vasp

c---------------------------------------------------------------------

      subroutine vasp_plot_CHG(avdir,cutdir,cutpos,filename)
      ! needs as input VASP CHGCAR file. Writes spin-resolved CHGCAR files
      ! and gnuplottable 
      ! charge density file CHG.4PLOTTING
      use defs
      use readcoords, only : read_poscar
      use writecoords
      use misc, only : vecs2vol, ithenearest,cross_product
      !use linalg
      implicit none 
      integer i,j,i0
      integer, intent(in) :: avdir,cutdir
      double precision :: cutpos
      character*1024 :: line,filename
      character*24 FORMA
      double precision :: vecs(1:3,1:3),vol,area
      double precision pos(3),posfrac(1:3)
      double precision, allocatable :: CD(:,:,:),abspos(:,:,:,:),         &
     &                          fracpos(:,:,:,:),CDav(:),MCD(:,:,:)
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms, nspecies,ngxf,ngyf,ngzf ! grid point numbers
      integer ix,iy,iz
      !
      print fsubstart,'vasp_plot_CHG'  
      !
      ! read coordinates, cell vectors, ... from header of CHGCAR.
      ! get cell volume to scale CD
      print '(8x,"Reading file ",A)',adjustl(trim(adjustl(filename)))
      call read_poscar(filename,atoms,natoms,species,nspecies,
     &           vecs)
      call vecs2vol(vecs,vol)
      print '(8x,"Cell volume (Angs³): ",F12.6)',vol
      !
      !
      ! begin read charge density and positions
      !
      open(51,file=filename,status='old', err=100)  
      rewind(51)
      read(51,'(A1024)', end=101,err=101) line
      do while (line.ne.'')
        read(51,'(A1024)', end=101,err=101) line
      end do
      read(51,*) ngxf,ngyf,ngzf
      print '(8x,"NGXF, NGYF, NGZF: ",3(I0,", "))', ngxf,ngyf,ngzf
      allocate(CD(ngxf,ngyf,ngzf),MCD(ngxf,ngyf,ngzf))
      allocate(abspos(ngxf,ngyf,ngzf,1:3))
      allocate(fracpos(ngxf,ngyf,ngzf,1:3))
      read(51,*) (((CD(iX,iY,iZ),iX=1,NGXF),iY=1,NGYF),iZ=1,NGZF) ! Charge density
      print '(8x,"Charge density read.")'
      !
      read(51,*,err=104,end=104) ngxf,ngyf,ngzf
      print '(8x,"NGXF, NGYF, NGZF: ",3(I0,", "))', ngxf,ngyf,ngzf
      read(51,*,end=104,err=104)(((MCD(iX,iY,iZ),iX=1,NGXF),iY=1,         &
     &      NGYF),iZ=1,NGZF) ! Magnetization charge density
      print '(8x,"Magnetization charge density read.")'
      close(51)
      !
      ! begin write magnetization and spin-polarized CHGCAR
      !
      ! MAG
      !
      call write_coords('CHG.MAG','poscar',atoms,natoms,species,          &
     &   nspecies,vecs)
      open(52,file='CHG.MAG',status='old',position='append',err=103) 
      write(52,'(" ")')
      write(52,'(3I5)') ngxf,ngyf,ngzf
      FORMA='(10(1X,G11.5))'
      write(52,FORMA) (((MCD(iX,iY,iZ),iX=1,NGXF),iY=1,NGYF),iZ=1,NGZF) 
      close(52)
      !
      ! UP
      !
      call write_coords('CHG.UP','poscar',atoms,natoms,species,           &
     &   nspecies,vecs)
      open(52,file='CHG.UP',status='old',position='append',err=103) 
      write(52,'(" ")')
      write(52,'(3I5)') ngxf,ngyf,ngzf
      FORMA='(10(1X,G11.5))'
      write(52,FORMA) ((((CD(iX,iY,iZ)+MCD(iX,iY,iZ))/2.0d0,iX=1,NGXF),    &
     &     iY=1,NGYF),iZ=1,NGZF) 
      close(52)
      !
      ! DOWN
      !
      call write_coords('CHG.DOWN','poscar',atoms,natoms,species,           &
     &   nspecies,vecs)
        open(52,file='CHG.DOWN',status='old',position='append',err=103) 
      write(52,'(" ")')
      write(52,'(3I5)') ngxf,ngyf,ngzf
      FORMA='(10(1X,G11.5))'
      write(52,FORMA) ((((CD(iX,iY,iZ)-MCD(iX,iY,iZ))/2.0d0,iX=1,NGXF),    &
     &     iY=1,NGYF),iZ=1,NGZF) 
      close(52)
      goto 105
      !
      ! end write magnetization and spin-polarized CHGCAR
      !
104   continue    
      close(51)
      close(52)
      nwarn=nwarn+1
      print fwarn,'cannot read mag. density. Ok for PARCHG or ISPIN=1'
      deallocate(MCD)
      !
105   continue 
      !     
      do ix=1,ngxf
       do iy=1,ngyf
        do iz=1,ngzf
         abspos(ix,iy,iz,1)=dble(mod(ix,ngxf))/dble(ngxf)*vecs(1,1)       &
     &                     +dble(mod(iy,ngyf))/dble(ngyf)*vecs(2,1)       &
     6                     +dble(mod(iz,ngzf))/dble(ngzf)*vecs(3,1)   
         abspos(ix,iy,iz,2)=dble(mod(ix,ngxf))/dble(ngxf)*vecs(1,2)       &
     &                     +dble(mod(iy,ngyf))/dble(ngyf)*vecs(2,2)       &
     &                     +dble(mod(iz,ngzf))/dble(ngzf)*vecs(3,2)
         abspos(ix,iy,iz,3)=dble(mod(ix,ngxf))/dble(ngxf)*vecs(1,3)       &
     &                     +dble(mod(iy,ngyf))/dble(ngyf)*vecs(2,3)       &
     &                     +dble(mod(iz,ngzf))/dble(ngzf)*vecs(3,3)
         fracpos(ix,iy,iz,1:3)=(/dble(mod(ix,ngxf))/dble(ngxf),           &
     &                           dble(mod(iy,ngyf))/dble(ngyf),           &
     &                           dble(mod(iz,ngzf))/dble(ngzf)/)
        end do ! igzf
       end do ! igyf
      end do ! igxf
      !
      ! end read charge density and positions
      !
!      !
!      ! begin write whole CHGCAR file in gnuplottable format. Huge
!      file, takes ages to write...
!      !
!      open(52,file='CHG4PLOTTING',status='replace', err=102)  
!      write(52,'("# IGXF/ngxf, IGYF/ngyf, IGZF/NGZF (fractional), the sa  &
!     &me absolute, CD(IGXF,IGYF,IGZF)")') 
!      do ix=1,ngxf
!       do iy=1,ngyf
!        do iz=1,ngzf
!         write(52,'(3(F9.6,1x),4(F12.6))') fracpos(ix,iy,iz,1:3),         &
!     &            abspos(ix,iy,iz,1:3),CD(ix,iy,iz)/vol 
!        end do ! igzf
!        if(ix.ne.ngxf.and.iy.ne.ngyf) write(52,'("")')
!       end do ! igyf
!       if(.not.ix==ngxf) write(52,'("")')
!      end do ! igxf
!      close(52)
      !
      ! end write whole CHGCAR file in gnuplottable format
      !
      ! begin average CHG perpendicular to avdir
      !
      select case(avdir) 
        case(1)
          print '(8x,"Averaging over cell vectors 2 and 3")' 
          area=norm2(cross_product(vecs(2,1:3),vecs(3,1:3)))
          allocate(CDav(ngxf))
          CDav=0.0d0
          do ix=1,ngxf
           do iy=1,ngyf
            do iz=1,ngzf
             CDav(ix)=CDav(ix)+CD(ix,iy,iz) ! could be done with sum(..) 
            end do ! iz
           end do ! iy
          end do ! ix
          CDav=CDav/(dble(ngyf*ngzf))
          open(52,file='CHGAV',status='replace', err=202)  
          write(52,'("# fractional x, absolute x y z, CDav(x), area, vol  &
     &/area")') 
          do ix=ngxf,ngxf ! first write last component (corresponds to x=0)
             write(52,'(1(F9.6,1x),6(F15.9))') fracpos(ix,1,1,1),         &
     &                abspos(ix,1,1,1:3),CDav(ix)/vol,area,vol/area
          end do ! igxf
          do ix=1,ngxf-1
             write(52,'(1(F9.6,1x),6(F15.9))') fracpos(ix,1,1,1),         &
     &                abspos(ix,1,1,1:3),CDav(ix)/vol,area,vol/area
          end do ! igxf
          close(52)
        case(2)
          print '(8x,"Averaging over cell vectors 1 and 3")' 
          area=norm2(cross_product(vecs(1,1:3),vecs(3,1:3)))
          allocate(CDav(ngyf))
          CDav=0.0d0
          do iy=1,ngyf
           do ix=1,ngxf
            do iz=1,ngzf
             CDav(iy)=CDav(iy)+CD(ix,iy,iz) ! could be done with sum(..) 
            end do ! iz
           end do ! ix
          end do ! iy
          CDav=CDav/(dble(ngxf*ngzf))
          open(52,file='CHGAV',status='replace', err=202)  
          write(52,'("# fractional y, absolute x y z, CDav(y), area, vol  &
     &/area")') 
          do iy=ngyf,ngyf ! first write last component (corresponds to y=0)
             write(52,'(1(F9.6,1x),6(F15.9))') fracpos(1,iy,1,2),         &
     &                abspos(1,iy,1,1:3),CDav(iy)/vol,area,vol/area
          end do ! igyf
          do iy=1,ngyf-1
             write(52,'(1(F9.6,1x),6(F15.9))') fracpos(1,iy,1,2),         &
     &                abspos(1,iy,1,1:3),CDav(iy)/vol,area,vol/area
          end do ! igyf
          close(52)
        case(3)
          print '(8x,"Averaging over cell vectors 1 and 2")' 
          area=norm2(cross_product(vecs(1,1:3),vecs(2,1:3)))
          allocate(CDav(ngzf))
          CDav=0.0d0
          do iz=1,ngzf
           do ix=1,ngxf
            do iy=1,ngyf
             CDav(iz)=CDav(iz)+CD(ix,iy,iz) ! could be done with sum(..) 
            end do ! iy
           end do ! ix
          end do ! iz
          CDav=CDav/(dble(ngxf*ngyf))
          open(52,file='CHGAV',status='replace', err=202)  
          write(52,'("# fractional z, absolute x y z, CDav(z), area, vol  &
     &/area")') 
          do iz=ngzf,ngzf ! first write last component (corresponds to z=0)
             write(52,'(1(F9.6,1x),6(F15.9))') fracpos(1,1,iz,3),         &
     &                abspos(1,1,iz,1:3),CDav(iz)/vol,area,vol/area
          end do ! igzf
          do iz=1,ngzf-1
             write(52,'(1(F9.6,1x),6(F15.9))') fracpos(1,1,iz,3),         &
     &                abspos(1,1,iz,1:3),CDav(iz)/vol,area,vol/area
          end do ! igzf
          close(52)
        case default
          goto 200
      end select
      print '(8x,"Charge average = ",F20.6," e")',                        &
     &        sum(CDav)/dble(size(CDav))
      !
      ! end average CHGCAR perpendicular to avdir
      !
      ! begin get CHGCAR on cut plane
      !
      do while (cutpos.lt.0.0d0) 
        cutpos=cutpos+1.0
      end do
      do while (cutpos.gt.1.0d0) 
        cutpos=cutpos-1.0
      end do
      select case(cutdir)
        case(1)
          print '(8x,"Getting CHG on plane spanned by vectors 2 and 3")'
          print '(8x,"with x1(frac)=",F9.6)',cutpos
          i0=ithenearest(cutpos,fracpos(:,1,1,1))
          print '(8x,"with x1(frac) really=",F9.6)',fracpos(i0,1,1,1)
          ! write
          open(52,file='CHGPLANE',status='replace', err=102)  
          write(52,'("# fractional x y z, absolute x y z, CD(y,z)")') 
          do ix=i0,i0
           do iy=1,ngyf
            do iz=1,ngzf
             write(52,'(3(F9.6,1x),4(F12.6))') fracpos(ix,iy,iz,1:3),     &
     &                abspos(ix,iy,iz,1:3),CD(ix,iy,iz)/vol 
            end do ! igzf
            if(iy.ne.ngyf) write(52,'("")')
           end do ! igyf
          end do ! igxf
          close(52)
        case(2)
          print '(8x,"Getting CHG on plane spanned by vectors 1 and 3")'
          print '(8x,"with x2(frac)=",F9.6)',cutpos
          i0=ithenearest(cutpos,fracpos(1,:,1,2))
          print '(8x,"with x2(frac) really=",F9.6)',fracpos(1,i0,1,2)
          ! write
          open(52,file='CHGPLANE',status='replace', err=102)  
          write(52,'("# fractional x y z, absolute x y z, CD(x,z)")') 
          do iy=i0,i0
           do ix=1,ngxf
            do iz=1,ngzf
             write(52,'(3(F9.6,1x),4(F12.6))') fracpos(ix,iy,iz,1:3),     &
     &                abspos(ix,iy,iz,1:3),CD(ix,iy,iz)/vol 
            end do ! igzf
            if(ix.ne.ngxf) write(52,'("")')
           end do ! igxf
          end do ! igyf
          close(52)
        case(3)
          print '(8x,"Getting CHG on plane spanned by vectors 1 and 2")'
          print '(8x,"with x3(frac)=",F9.6)',cutpos
          i0=ithenearest(cutpos,fracpos(1,1,:,3))
          print '(8x,"with x3(frac) really=",F9.6)',fracpos(1,1,i0,3)
          ! write
          open(52,file='CHGPLANE',status='replace', err=102)  
          write(52,'("# fractional x y z, absolute x y z, CD(x,y)")') 
          do iz=i0,i0
           do ix=1,ngxf
            do iy=1,ngyf
             write(52,'(3(F9.6,1x),4(F12.6))') fracpos(ix,iy,iz,1:3),         &
     &                abspos(ix,iy,iz,1:3),CD(ix,iy,iz)/vol 
            end do ! igyf
            if(ix.ne.ngxf) write(52,'("")')
           end do ! igxf
          end do ! igzf
          close(52)
        case default
          goto 300
      end select
      !
      ! end get CHG on cut plane
      !
      !
      deallocate(CD)
      deallocate(CDav)
      if(allocated(MCD)) deallocate(MCD)
      print fsubendext,'vasp_plot_CHG'  
      return
      !
100   continue    
      close(51)
      if(allocated(CD)) deallocate(CD)
      deallocate(MCD)
      nerr=nerr+1
      print ferrmssg,'CHGCAR not found'
      return
      !
101   continue    
      close(51)
      deallocate(CD)
      deallocate(MCD)
      nerr=nerr+1
      print ferrmssg,'something wrong with CHGCAR...'
      return
      !
102   continue    
      close(51)
      close(52)
      deallocate(CD)
      deallocate(MCD)
      nerr=nerr+1
      print ferrmssg,'cannot write CHG.4PLOTTING'
      return
      !
103   continue    
      close(51)
      close(52)
      deallocate(CD)
      deallocate(MCD)
      nerr=nerr+1
      print ferrmssg,'cannot write CHG.MAG'
      return
      !
!104   continue    
!      close(51)
!      close(52)
!      nerr=nerr+1
!      print ferrmssg,'cannot read mag. density'
!      deallocate(CD)
!      deallocate(MCD)
!      return
!      !
200   continue
      nerr=nerr+1
      deallocate(CD)
      deallocate(MCD)
      print ferrmssg,'please choose avdir between 1 and 3'
      return
202   continue    
      close(51)
      close(52)
      deallocate(CD)
      deallocate(MCD)
      deallocate(CDav)
      nerr=nerr+1
      print ferrmssg,'cannot write CHGAV'
      return
300   continue
      nerr=nerr+1
      deallocate(CD)
      deallocate(MCD)
      print ferrmssg,'please choose cutdir between 1 and 3'
      return
      !
       
      end subroutine vasp_plot_CHG

c---------------------------------------------------------------------

      subroutine vasp_CHG_cut_sphere(origin,radius,filename)
      !  
      ! needs as input VASP CHGCAR or PARCHG file. 
      ! Prints the charge density inside and outside a sphere with requested
      ! origin and radius
      ! 
      use defs
      use readcoords, only : read_poscar
      use writecoords
      use misc, only : vecs2vol, ithenearest,cross_product
      !
      implicit none 
      integer i,j,i0
      double precision :: origin(3),radius
      character*1024 :: line,filename
      character*24 FORMA
      double precision :: vecs(1:3,1:3) !,vol
      double precision pos(3),posfrac(1:3)
      double precision, allocatable :: CD(:,:,:),abspos(:,:,:,:),         &
     &                          fracpos(:,:,:,:),MCD(:,:,:)
      double precision, allocatable ::CD_inside(:,:,:),CD_outside(:,:,:)
      double precision, allocatable :: MCD_inside(:,:,:),                 &
     &MCD_outside(:,:,:)
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms, nspecies,ngxf,ngyf,ngzf ! grid point numbers
      integer ix,iy,iz,i1,i2,i3
      !
      print fsubstart,'vasp_CHG_cut_sphere'  
      !
      ! read coordinates, cell vectors, ... from header of CHGCAR.
      ! get cell volume to scale CD
      print '(8x,"Reading file ",A)',adjustl(trim(adjustl(filename)))
      call read_poscar(filename,atoms,natoms,species,nspecies,
     &           vecs)
      print '(8x,"atomic positions and lattice vectors read in.")'

      print '(8x,"Cutting out a sphere at fractional position ",          &
     &     3(F12.6)))', origin
      origin=matmul(transpose(vecs),origin)
      print '(8x,"Cutting out a sphere at absolute position ",            &
     &     3(F12.6)))', origin
      print '(8x,"Cutting out a sphere with radius ",F12.6," Angs")',     &
     &     radius
      !
      !
      ! begin read charge density and positions
      !
      open(51,file=filename,status='old', err=101)  
      rewind(51)
      read(51,'(A1024)', end=101,err=101) line
      do while (line.ne.'')
        read(51,'(A1024)', end=101,err=101) line
      end do
      read(51,*, err=101) ngxf,ngyf,ngzf
      print '(8x,"NGXF, NGYF, NGZF: ",3(I0,", "))', ngxf,ngyf,ngzf
      allocate(CD(ngxf,ngyf,ngzf),MCD(ngxf,ngyf,ngzf))
      allocate(CD_inside(ngxf,ngyf,ngzf),MCD_inside(ngxf,ngyf,ngzf))
      allocate(CD_outside(ngxf,ngyf,ngzf),MCD_outside(ngxf,ngyf,ngzf))
      allocate(abspos(ngxf,ngyf,ngzf,1:3))
      allocate(fracpos(ngxf,ngyf,ngzf,1:3))
      print '(8x,"arrays allocated.")'
      read(51,*,err=101) (((CD(iX,iY,iZ),iX=1,NGXF),iY=1,NGYF),iZ=1,      &
     &     NGZF) ! Charge density
      print '(8x,"Charge density read.")'
      !
      read(51,*,err=14,end=14) ngxf,ngyf,ngzf
      print '(8x,"NGXF, NGYF, NGZF: ",3(I0,", "))', ngxf,ngyf,ngzf
      read(51,*,end=14,err=14)(((MCD(iX,iY,iZ),iX=1,NGXF),iY=1,           &
     &      NGYF),iZ=1,NGZF) ! Magnetization charge density
      print '(8x,"Magnetization charge density read.")'
      close(51)
      goto 15
      !
14    continue    
      close(51)
      close(52)
      nwarn=nwarn+1
      print fwarn,'cannot read mag. density. Ok for PARCHG or ISPIN=1'
      deallocate(MCD)
      deallocate(MCD_inside)
      deallocate(MCD_outside)
      !
15    continue 
      !     
      print '(8x,"setting up spatial grid...")'
      do ix=1,ngxf
       do iy=1,ngyf
        do iz=1,ngzf
         abspos(ix,iy,iz,1)=dble(mod(ix,ngxf))/dble(ngxf)*vecs(1,1)       &
     &                     +dble(mod(iy,ngyf))/dble(ngyf)*vecs(2,1)       &
     6                     +dble(mod(iz,ngzf))/dble(ngzf)*vecs(3,1)   
         abspos(ix,iy,iz,2)=dble(mod(ix,ngxf))/dble(ngxf)*vecs(1,2)       &
     &                     +dble(mod(iy,ngyf))/dble(ngyf)*vecs(2,2)       &
     &                     +dble(mod(iz,ngzf))/dble(ngzf)*vecs(3,2)
         abspos(ix,iy,iz,3)=dble(mod(ix,ngxf))/dble(ngxf)*vecs(1,3)       &
     &                     +dble(mod(iy,ngyf))/dble(ngyf)*vecs(2,3)       &
     &                     +dble(mod(iz,ngzf))/dble(ngzf)*vecs(3,3)
         fracpos(ix,iy,iz,1:3)=(/dble(mod(ix,ngxf))/dble(ngxf),           &
     &                           dble(mod(iy,ngyf))/dble(ngyf),           &
     &                           dble(mod(iz,ngzf))/dble(ngzf)/)
        end do ! igzf
       end do ! igyf
      end do ! igxf
      print '(8x,"...spatial grid has been set up.")'
      !
      ! end read charge density and positions
      !
      print '(8x,"partitioning charge density...")'
      !
      CD_inside=0.0d0
      if (allocated(MCD_inside)) MCD_inside=0.0d0
      CD_outside=CD
      if (allocated(MCD_outside)) MCD_outside=MCD
      do iX=1,ngxf
       do iY=1,ngyf
        do iZ=1,ngzf
         do i1=-1,1
         do i2=-1,1
         do i3=-1,1
         if (norm2(abspos(iX,iY,iZ,:)+dble(i1)*vecs(1,:)+dble(i2)         &
     &     *vecs(2,:)+dble(i3)*vecs(3,:)-origin(:)).le.radius) then
          CD_inside(iX,iY,iZ)=CD(iX,iY,iZ)
          if (allocated(MCD)) MCD_inside(iX,iY,iZ)=MCD(iX,iY,iZ)
          CD_outside(iX,iY,iZ)=0.0d0
          if (allocated(MCD)) MCD_outside(iX,iY,iZ)=0.0d0
         !else
          !CD_inside(iX,iY,iZ)=0.0d0
          !if (allocated(MCD)) MCD_inside(iX,iY,iZ)=0.0d0
          !CD_outside(iX,iY,iZ)=CD(iX,iY,iZ)
          !if (allocated(MCD)) MCD_outside(iX,iY,iZ)=MCD(iX,iY,iZ)
         end if
         end do
         end do
         end do
        end do ! igzf
       end do ! igyf
      end do ! igxf
      print '(8x,"...done.")'
      !
      ! begin write CD inside sphere
      !
      print '(8x,"writing charge density inside sphere...")'
      call write_coords('CHG_INSIDE_SPHERE','poscar',atoms,natoms,        &
     &                  species,nspecies,vecs)
      open(52,file='CHG_INSIDE_SPHERE',status='old',position='append',    &
     &     err=102) 
      write(52,'(" ")')
      write(52,'(3I5)') ngxf,ngyf,ngzf
      FORMA='(10(1X,G11.5))'
      write(52,FORMA) (((CD_inside(iX,iY,iZ),iX=1,NGXF),                  &
     &     iY=1,NGYF),iZ=1,NGZF) 
      if (allocated(MCD_inside)) then
        write(52,'(3I5)') ngxf,ngyf,ngzf
        FORMA='(10(1X,G11.5))'
        write(52,FORMA) (((MCD_inside(iX,iY,iZ),iX=1,NGXF),               &
     &       iY=1,NGYF),iZ=1,NGZF) 
      end if
      close(52)
      print '(8x,"...done.")'
      !
      print '(8x,"and outside sphere...")'
      call write_coords('CHG_OUTSIDE_SPHERE','poscar',atoms,natoms,       &
     &                  species,nspecies,vecs)
      open(52,file='CHG_OUTSIDE_SPHERE',status='old',position='append',   &
     &     err=103) 
      write(52,'(" ")')
      write(52,'(3I5)') ngxf,ngyf,ngzf
      FORMA='(10(1X,G11.5))'
      write(52,FORMA) (((CD_outside(iX,iY,iZ),iX=1,NGXF),                 &
     &     iY=1,NGYF),iZ=1,NGZF) 
      if (allocated(MCD_outside)) then
        write(52,'(3I5)') ngxf,ngyf,ngzf
        FORMA='(10(1X,G11.5))'
        write(52,FORMA) (((MCD_outside(iX,iY,iZ),iX=1,NGXF),                &
     &       iY=1,NGYF),iZ=1,NGZF) 
      end if
      close(52)
      print '(8x,"...done.")'
      !
      !
      !
      deallocate(CD)
      deallocate(CD_inside)
      deallocate(CD_outside)
      if(allocated(MCD)) deallocate(MCD)
      if(allocated(MCD_inside)) deallocate(MCD_inside)
      if(allocated(MCD_outside)) deallocate(MCD_outside)
      print fsubendext,'vasp_CHG_cut_sphere'  
      return
      !
      !
101   continue    
      close(51)
      if(allocated(CD)) deallocate(CD)
      if(allocated(CD_inside)) deallocate(CD_inside)
      if(allocated(CD_outside)) deallocate(CD_outside)
      if(allocated(MCD)) deallocate(MCD)
      if(allocated(MCD_inside)) deallocate(MCD_inside)
      if(allocated(MCD_outside)) deallocate(MCD_outside)
      nerr=nerr+1
      print ferrmssg,'something wrong with Charge density file...'
      return
      !
102   continue    
      close(51)
      close(52)
      deallocate(CD)
      deallocate(CD_inside)
      deallocate(CD_outside)
      deallocate(MCD)
      deallocate(MCD_inside)
      deallocate(MCD_outside)
      nerr=nerr+1
      print ferrmssg,'cannot write CHG_INSIDE_SPHERE'
      return
      !
103   continue    
      close(51)
      close(52)
      deallocate(CD)
      deallocate(CD_inside)
      deallocate(CD_outside)
      deallocate(MCD)
      deallocate(MCD_inside)
      deallocate(MCD_outside)
      nerr=nerr+1
      print ferrmssg,'cannot write CHG_OUTSIDE_SPHERE'
      return
      !
      end subroutine vasp_CHG_cut_sphere

c---------------------------------------------------------------------

      subroutine vasp_CHG_overlap(filename1,filename2)
      ! needs as input two VASP CHGCAR files. 
      ! Writes overlap = int (CHG1 * CHG2) d³r
      use defs
      use readcoords, only : read_poscar
      use writecoords, only : write_poscar
      use misc, only : vecs2vol, ithenearest,cross_product
      !use linalg
      implicit none 
      integer i,j,i0
      !integer, intent(in) :: avdir,cutdir
      !double precision :: cutpos
      character*1024 :: line,filename1,filename2
      double precision :: vecs1(1:3,1:3),vecs2(1:3,1:3),vol1,vol2
      double precision pos1(3),posfrac1(1:3)
      double precision pos2(3),posfrac2(1:3)
      double precision, allocatable :: CD1(:,:,:),CD2(:,:,:),             &
     &     CDdiff(:,:,:)  ,                                                &
     &     abspos1(:,:,:,:),fracpos1(:,:,:,:),abspos2(:,:,:,:),           &
     &     fracpos2(:,:,:,:)
      type(atom), allocatable :: atoms1(:),atoms2(:)
      type(element), allocatable :: species1(:),species2(:)
      integer natoms1,nspecies1,ngxf1,ngyf1,ngzf1 ! # of atoms, # of species, grid point numbers in first file
      integer natoms2,nspecies2,ngxf2,ngyf2,ngzf2 ! # of atoms, # of species, grid point numbers in 2. file
      integer ix,iy,iz,ix2,iy2,iz2
      double precision overlap, CD1norm,CD2norm
      !
      print fsubstart,'vasp_CHG_overlap'  
      !
      ! read coordinates, cell vectors, ... from header of CHGCAR1.
      ! get cell volume1 to scale CD1
      print '(8x,"Reading file ",A)',adjustl(trim(adjustl(filename1)))
      call read_poscar(filename1,atoms1,natoms1,species1,nspecies1,
     &           vecs1)
      call vecs2vol(vecs1,vol1)
      print '(8x,"1. Cell volume (Angs³): ",F12.6)',vol1
      !
      ! read coordinates, cell vectors, ... from header of CHGCAR2.
      ! get cell volume2 to scale CD2
      print '(8x,"Reading file ",A)',adjustl(trim(adjustl(filename2)))
      call read_poscar(filename2,atoms2,natoms2,species2,nspecies2,
     &           vecs2)
      call vecs2vol(vecs2,vol2)
      print '(8x,"2. Cell volume (Angs³): ",F12.6)',vol2
      !
      !
      ! begin read charge density and positions from 1. file
      !
      print '(8x," ")'
      open(51,file=filename1,status='old', err=100)  
      rewind(51)
      read(51,'(A1024)', end=101,err=101) line
      do while (line.ne.'')
        read(51,'(A1024)', end=101,err=101) line
        !print *,line
      end do
      read(51,*) ngxf1,ngyf1,ngzf1
      print '(8x,"NGXF1, NGYF1, NGZF1: ",3(I0,", "))', ngxf1,ngyf1,ngzf1
      allocate(CD1(ngxf1,ngyf1,ngzf1))
      allocate(abspos1(ngxf1,ngyf1,ngzf1,1:3))
      allocate(fracpos1(ngxf1,ngyf1,ngzf1,1:3))
      read(51,*) (((CD1(iX,iY,iZ),iX=1,NGXF1),iY=1,NGYF1),iZ=1,NGZF1)
      close(51)
      CD1norm=0.0d0
      do ix=1,ngxf1
       do iy=1,ngyf1
        do iz=1,ngzf1
         abspos1(ix,iy,iz,1)=dble(mod(ix,ngxf1))/dble(ngxf1)*vecs1(1,1)   &
     &                     +dble(mod(iy,ngyf1))/dble(ngyf1)*vecs1(2,1)    &
     6                     +dble(mod(iz,ngzf1))/dble(ngzf1)*vecs1(3,1)  
         abspos1(ix,iy,iz,2)=dble(mod(ix,ngxf1))/dble(ngxf1)*vecs1(1,2)   &
     &                     +dble(mod(iy,ngyf1))/dble(ngyf1)*vecs1(2,2)    &
     &                     +dble(mod(iz,ngzf1))/dble(ngzf1)*vecs1(3,2)
         abspos1(ix,iy,iz,3)=dble(mod(ix,ngxf1))/dble(ngxf1)*vecs1(1,3)   &
     &                     +dble(mod(iy,ngyf1))/dble(ngyf1)*vecs1(2,3)    &
     &                     +dble(mod(iz,ngzf1))/dble(ngzf1)*vecs1(3,3)
         fracpos1(ix,iy,iz,1:3)=(/dble(mod(ix,ngxf1))/dble(ngxf1),        &
     &                           dble(mod(iy,ngyf1))/dble(ngyf1),         &
     &                           dble(mod(iz,ngzf1))/dble(ngzf1)/)
         CD1norm=CD1norm+CD1(ix,iy,iz)
        end do ! igzf1
       end do ! igyf1
      end do ! igxf1
      CD1norm=CD1norm/dble(ngxf1*ngyf1*ngzf1)
      CD1=CD1/CD1norm
      print '(8x,"CD1 norm=",F12.6)', CD1norm
      !
      ! end read charge density and positions from first file
      !
      !
      ! begin read charge density and positions from 2. file
      !
      print '(8x," ")'
      open(51,file=filename2,status='old', err=100)  
      rewind(51)
      read(51,'(A1024)', end=101,err=101) line
      do while (line.ne.'')
        read(51,'(A1024)', end=101,err=101) line
        !print *,line
      end do
      read(51,*) ngxf2,ngyf2,ngzf2
      print '(8x,"NGXF2, NGYF2, NGZF2: ",3(I0,", "))', ngxf2,ngyf2,ngzf2
      allocate(CD2(ngxf2,ngyf2,ngzf2))
      allocate(abspos2(ngxf2,ngyf2,ngzf2,1:3))
      allocate(fracpos2(ngxf2,ngyf2,ngzf2,1:3))
      read(51,*) (((CD2(iX,iY,iZ),iX=1,NGXF2),iY=1,NGYF2),iZ=1,NGZF2)
      close(51)
      CD2norm=0.0d0
      do ix=1,ngxf2
       do iy=1,ngyf2
        do iz=1,ngzf2
         abspos2(ix,iy,iz,1)=dble(mod(ix,ngxf2))/dble(ngxf2)*vecs2(1,1)   &
     &                     +dble(mod(iy,ngyf2))/dble(ngyf2)*vecs2(2,1)    &
     6                     +dble(mod(iz,ngzf2))/dble(ngzf2)*vecs2(3,1)  
         abspos2(ix,iy,iz,2)=dble(mod(ix,ngxf2))/dble(ngxf2)*vecs2(1,2)   &
     &                     +dble(mod(iy,ngyf2))/dble(ngyf2)*vecs2(2,2)    &
     &                     +dble(mod(iz,ngzf2))/dble(ngzf2)*vecs2(3,2)
         abspos2(ix,iy,iz,3)=dble(mod(ix,ngxf2))/dble(ngxf2)*vecs2(1,3)   &
     &                     +dble(mod(iy,ngyf2))/dble(ngyf2)*vecs2(2,3)    &
     &                     +dble(mod(iz,ngzf2))/dble(ngzf2)*vecs2(3,3)
         fracpos2(ix,iy,iz,1:3)=(/dble(mod(ix,ngxf2))/dble(ngxf2),        &
     &                           dble(mod(iy,ngyf2))/dble(ngyf2),         &
     &                           dble(mod(iz,ngzf2))/dble(ngzf2)/)
         CD2norm=CD2norm+CD2(ix,iy,iz)
        end do ! igzf2
       end do ! igyf2
      end do ! igxf2
      CD2norm=CD2norm/dble(ngxf2*ngyf2*ngzf2)
      CD2=CD2/CD2norm
      print '(8x,"CD2 norm=",F12.6)', CD2norm
      print '(8x," ")'
      !
      ! end read charge density and positions from 2. file
      !
      ! begin integrate sqrt(CD1 * CD2 )
      !
      overlap=0.0d0
      !
      do ix=1,ngxf1
       do iy=1,ngyf1
        do iz=1,ngzf1
          ix2=ithenearest(fracpos1(ix,iy,iz,1),fracpos2(:,1,1,1))
          iy2=ithenearest(fracpos1(ix,iy,iz,2),fracpos2(1,:,1,2))
          iz2=ithenearest(fracpos1(ix,iy,iz,3),fracpos2(1,1,:,3))
          overlap=overlap+sqrt(max(CD1(ix,iy,iz),0.0d0)*                  &
     &                         max(CD2(ix2,iy2,iz2),0.0d0))
        end do ! iz 
       end do ! iy
      end do ! ix  
      overlap=overlap/dble(ngxf1*ngyf1*ngzf1)
      print '(8x,"overlap=",F12.6)', overlap
      !
      ! end integrate sqrt(CD1 * CD2 )
      !
      ! begin calculate difference CD1-CD2
      !
      allocate(CDdiff(ngxf1,ngyf1,ngzf1))
      do ix=1,ngxf1
       do iy=1,ngyf1
        do iz=1,ngzf1
          ix2=ithenearest(fracpos1(ix,iy,iz,1),fracpos2(:,1,1,1))
          iy2=ithenearest(fracpos1(ix,iy,iz,2),fracpos2(1,:,1,2))
          iz2=ithenearest(fracpos1(ix,iy,iz,3),fracpos2(1,1,:,3))
          !CDdiff(ix,iy,iz)=abs(CD1(ix,iy,iz)-CD2(ix2,iy2,iz2))
          CDdiff(ix,iy,iz)=CD1(ix,iy,iz)-CD2(ix2,iy2,iz2)
        end do ! iz 
       end do ! iy
      end do ! ix  
      !
      call write_poscar("CD_diff.vasp",atoms1,natoms1,species1,           &
     &                    nspecies1,vecs1)
      open(51,file="CD_diff.vasp",status='old', position='append',        &
     &      err=400)  
      write(51,'("")')
      write(51,'(3(I5))') ngxf1,ngyf1,ngzf1
      write(51,'(5E19.11)')(((CDdiff(iX,iY,iZ),iX=1,NGXF1),iY=1,NGYF1),   &
     &   iZ=1,NGZF1)
      close(51)
      !
      call write_poscar("CD_abs_diff.vasp",atoms1,natoms1,species1,       &
     &                    nspecies1,vecs1)
      open(51,file="CD_abs_diff.vasp",status='old', position='append',    &
     &      err=400)  
      write(51,'("")')
      write(51,'(3(I5))') ngxf1,ngyf1,ngzf1
      !write(51,'(5E18.11)')(((CDdiff(iX,iY,iZ),iX=1,NGXF1),iY=1,NGYF1),   &
!     &   iZ=1,NGZF1)
      write(51,'(5E18.11)')(((abs(CDdiff(iX,iY,iZ)),iX=1,NGXF1),iY=1,     &
     &   NGYF1),iZ=1,NGZF1)
      close(51)
      !
      ! end calculate difference CD2-CD1
      !
      deallocate(CD1,CD2,CDdiff)
      !
      print fsubendext,'vasp_CHG_overlap'  
      return
      !
100   continue    
      close(51)
      deallocate(CD1,CD2)
      nerr=nerr+1
      print ferrmssg,'CHGCAR not found'
      return
      !
101   continue    
      close(51)
      deallocate(CD1,CD2)
      nerr=nerr+1
      print ferrmssg,'something wrong with CHGCAR...'
      return
      !
102   continue    
      close(51)
      close(52)
      deallocate(CD1,CD2)
      nerr=nerr+1
      print ferrmssg,'cannot write CHG.4PLOTTING'
      return
      !
200   continue
      nerr=nerr+1
      deallocate(CD1,CD2)
      print ferrmssg,'please choose avdir between 1 and 3'
      return
202   continue    
      close(51)
      close(52)
      deallocate(CD1,CD2)
      nerr=nerr+1
      print ferrmssg,'cannot write CHGAV'
      return
300   continue
      nerr=nerr+1
      deallocate(CD1,CD2)
      print ferrmssg,'please choose cutdir between 1 and 3'
      return
400   continue
      nerr=nerr+1
      deallocate(CD1,CD2,CDdiff)
      print ferrmssg,'Problem with writing CHG difference'
      return
      !
       
      end subroutine vasp_CHG_overlap

c---------------------------------------------------------------------

      subroutine vasp_dE_dk(n,cutdir,cutpos)
      use defs
      use misc, only : ithenearest,frac2abs
      implicit none    
      !  
      integer n ! band index
      integer cutdir ! cutting direction
      double precision cutpos ! cutting position along cutdir for 2D cut
      double precision cutposreal ! real cutting position
      double precision fdum
      character*1024 :: eigenval,outcar,line,filename
      integer nbands,iband,nkpoints,nele,iline,idum,i0
      integer ikpoint,jkpoint,ikpointlast
      integer icart,innkpt,i1,i2,i3
      integer,  allocatable :: innk(:,:)
      double precision kdistvec(1:3), kdists_old(1:3),kdists(1:3)
      integer iwrite,iread,i,inn_temp(100)
      logical spinpol,new
      type(kpoint), allocatable :: k(:)
      type(kpoint) kstore
      double precision tol,tol3(1:3),tolfac_k ! tolfac_k: tolerance factor for NN search. k-points within tolfac_k*shortest_NN_dist are nearest neighbors
      double precision ktemp(1:3),nn_temp(100,1:3)
      ! ddk:
      integer nddk
      double precision e_up, e_down, E_nnk_up,E_nnk_down,dE(1:2)
      double precision gvecs(3,3),kcoords(3),knn_coords(3),dk,s
      !
      print fsubstart,'vasp_dE_dk'
      !
      ! begin initialize
      !
      eigenval="EIGENVAL"
      outcar="OUTCAR"
      spinpol=.false.
      tol=1D-5
      tol3=tol
      tolfac_k=1.1d0
      open(54,file="DEBUG.out",status="replace")
      !
      ! end initialize
      !
      !
      ! begin read Fermi energy & reciprocal lattice vecs
      !
      open(51,file=outcar,status='old', err=100)
10    read(51,'(A256)',err=100,end=20) line
      if (index(line,'reciprocal lattice vectors').gt.0) then
        do i1=1,3
          read (51,*,end=20,err=100) fdum,fdum,fdum,gvecs(i1,1:3)
        end do
        gvecs=gvecs*2.0d0*pi
      end if 
      goto 10
      !
20    continue
      close(51)
      print '(8x,"reciprocal lattice vecs:",3(/,8x,3(F12.6)))',           &
     &     (gvecs(i1,1:3),i1=1,3)
      !
      ! end read Fermi energy & reciprocal lattice vecs
      !
      open(51,file=eigenval, status='old', err=101)
      !
      ! begin read nbands, nkpoints
      !    
      read(51,*,end=101,err=101) idum,idum,idum, idum  
      if (idum==2) spinpol=.true.   
      do iline=1,4
        read(51,'(A1024)',end=100,err=100) line    
      end do
      read(51,*,end=101,err=101) nele,nkpoints,nbands
      print '(8x,I0," electrons, ",I0," kpoints, ",I0," bands")',nele,    &
     &  nkpoints,nbands
      allocate(k(1:nkpoints))
      !
      ! end read nbands, nkpoints
      !
      !
      ! begin read fractional kpoints and eigenvalues
      !    
      do ikpoint=1,nkpoints
        if (.not.allocated(k(ikpoint)%e_up))                              &
     &       allocate(k(ikpoint)%e_up(nbands))
        if (.not.allocated(k(ikpoint)%e_down))                            &
     &       allocate(k(ikpoint)%e_down(nbands))
        read(51,'(A1024)',end=101,err=101) line    
        read(51,*,end=101,err=101) k(ikpoint)%kptfrac(1:3),               &
     &                             k(ikpoint)%weight
        do iband=1,nbands
          if (.not.spinpol) then
            read(51,*,err=101,end=101) idum,                              &
     &         k(ikpoint)%e_up(iband)
            k(ikpoint)%e_down(iband)=k(ikpoint)%e_up(iband)
          else
            read(51,*,err=101,end=101) idum,                              &
     &         k(ikpoint)%e_up(iband),                                    &
     &         k(ikpoint)%e_down(iband)
          end if
        end do
      end do
      close(51)
      !
      ! end read fractional kpoints and eigenvalues
      !
      !
      ! begin calculate cartesian kpoint coordinates
      !
      do ikpoint=1,nkpoints
        call frac2abs(k(ikpoint)%kptfrac,gvecs,k(ikpoint)%kpt)
      end do ! ikpoint
      !
      ! end calculate cartesian kpoint coordinates
      !
      ! begin find nearest neighbor indices in cartesian k1,k2,k3 direction
      !
      allocate(innk(nkpoints,1:3))
      innk=0
      do ikpoint=1,nkpoints
        do icart=1,3 ! reciprocal lattice direction
          kdists_old=1.0E6
          do jkpoint=1,nkpoints
            if (jkpoint.ne.ikpoint) then
              do i1=-1,1 ! periodic images
              do i2=-1,1
              do i3=-1,1
                    kdistvec=k(jkpoint)%kpt                               &
     &                      +dble(i1)*gvecs(1,:)+dble(i2)*gvecs(2,:)      &
     &                      +dble(i3)*gvecs(3,:)-k(ikpoint)%kpt
                    kdists=abs(kdistvec)
                    if(abs(kdistvec(icart)).le.kdists_old(icart)          &
     &                 +tol3(icart).and.abs(kdistvec(icart))              &
     &                   .gt.tol3(icart).and.                             &
     &                   all(kdists.le.kdists_old+tol3)) then
                        kdists_old=kdists
                        innk(ikpoint,icart)=jkpoint
                        k(ikpoint)%nn_dist_min(icart,1:3)=kdists
                    end if ! 
              end do ! i3
              end do ! i2
              end do ! i1
            end if ! ikpoint.ne.jkpoint
          end do ! jkpoint
        end do ! icart
        do icart=1,3
          write(54,'("kpt ",I0,": fcoords=",3(F10.6),", acoords=",        &
     &3(F10.6),", NN(x",I0,"): kpt ",I0,", fcoords=",                     &
     &3(F10.6),", acoords=",3(F10.6),", abs dk=",3(F10.6))'               &
     &) ikpoint,k(ikpoint)%kptfrac,k(ikpoint)%kpt,icart,                  &
     &innk(ikpoint,icart),k(innk(ikpoint,icart))%kptfrac,k(innk(ikpoint,  &
     &icart))%kpt,k(innk(ikpoint,icart))%kpt-k(ikpoint)%kpt ! DEBUG
        end do ! icart
      end do ! ikpoint
      !
      ! end find nearest neighbor indices in cartesian k1,k2,k3 direction
      !
      ! begin find nearest neighbors in cartesian k1,k2,k3 direction
      !
      do ikpoint=1,nkpoints
        k(ikpoint)%n_neighbors=0
        nn_temp=-10000.0d0
        do icart=1,3 ! cartesian reciprocal lattice direction
          !kdists_old=abs(k(innk(ikpoint,icart))%kpt-k(ikpoint)%kpt)
          kdists_old=k(ikpoint)%nn_dist_min(icart,:)
          do jkpoint=1,nkpoints
            if (jkpoint.ne.ikpoint) then
              do i1=-1,1 ! periodic images
              do i2=-1,1
              do i3=-1,1
                    ktemp=k(jkpoint)%kpt                                  &
     &                      +dble(i1)*gvecs(1,:)+dble(i2)*gvecs(2,:)      &
     &                      +dble(i3)*gvecs(3,:)
                    kdistvec=ktemp-k(ikpoint)%kpt
                    kdists=abs(kdistvec)
                    if (all(kdists.le.kdists_old*tolfac_k)                &
     &                  .and.kdists(icart).gt.tol3(icart)) then
                      new=.true.
                      do i=1,k(ikpoint)%n_neighbors
                        if (all(nn_temp(i,1:3)==ktemp(1:3))) new=.false.
                      end do  
                      if (new) then
                        k(ikpoint)%n_neighbors=k(ikpoint)%n_neighbors+1
                        nn_temp(k(ikpoint)%n_neighbors,1:3)=ktemp  
                        inn_temp(k(ikpoint)%n_neighbors)=jkpoint
                      end if  
                    end if ! new
              end do ! i3
              end do ! i2
              end do ! i3
            end if ! ikpoint.ne.jkpoint
          end do ! jkpoint
        end do ! icart
        write(54,'(8x,"kpt no. ",I0,": # of nearest neighbors=",          &
     &          I0)') ikpoint,k(ikpoint)%n_neighbors ! DEBUG
        allocate(k(ikpoint)%neighbor_kpts(k(ikpoint)%n_neighbors,1:3))
        allocate(k(ikpoint)%i_neighbors(k(ikpoint)%n_neighbors))
        k(ikpoint)%i_neighbors=inn_temp
        k(ikpoint)%neighbor_kpts=nn_temp(1:k(ikpoint)%n_neighbors,1:3)
      end do  ! ikpoint
      !
      ! end find number of nearest neighbors in k1,k2,k3 direction
      !
!      close(51)
!      close(52)
!      close(53)
!      close(54)
!      return
      !
      ! begin store nearest neighbors in k1,k2,k3 direction
      !
      do ikpoint=1,nkpoints
!        innkpt=0
!        allocate(k(ikpoint)%neighbor_kpts(k(ikpoint)%n_neighbors,1:3))
!        allocate(k(ikpoint)%i_neighbors(k(ikpoint)%n_neighbors))
!        do icart=1,3 ! reciprocal lattice direction
!          kdists_old=abs(k(innk(ikpoint,icart))%kpt-k(ikpoint)%kpt)
!          do i1=-1,1 ! periodic images
!          do i2=-1,1
!          do i3=-1,1
!            kdists=abs(k(innk(ikpoint,icart))%kpt                         &
!     &                      +dble(i1)*gvecs(1,:)+dble(i2)*gvecs(2,:)      &
!     &                      +dble(i3)*gvecs(3,:)-k(ikpoint)%kpt)
!            if (all(kdists.lt.kdists_old+tol3)) kdists_old=kdists  
!          end do ! i3
!          end do ! i2
!          end do ! i3
!          do jkpoint=1,nkpoints
!            if (jkpoint.ne.ikpoint) then
!              do i1=-1,1 ! periodic images
!              do i2=-1,1
!              do i3=-1,1
!                  kdistvec=k(jkpoint)%kpt                                 &
!     &                      +dble(i1)*gvecs(1,:)+dble(i2)*gvecs(2,:)      &
!     &                      +dble(i3)*gvecs(3,:)-k(ikpoint)%kpt
!                  kdists=abs(kdistvec)
!                  if (all(kdists.le.kdists_old*tolfac_k)) then
!                    innkpt=innkpt+1  
!                    k(ikpoint)%neighbor_kpts(innkpt,1:3)                  &
!     &                =k(jkpoint)%kpt(1:3)+dble(i1)*gvecs(1,:)+dble(i2)   &
!     &                  +dble(i3)*gvecs(3,:)       
!                    k(ikpoint)%i_neighbors(innkpt)=jkpoint
!                  end if ! 
!              end do ! i3 
!              end do ! i2
!              end do ! i1 
!            end if ! ikpoint.ne.jkpoint
!          end do ! jkpoint
!        end do ! icart
        write(54,'("kpt ",I0,": fcoords=",3(F10.6),", acoords=",          &
     &3(F10.6),", ",I0," NN:"," acoords= ",20(3(F10.6),1x))'              &
     &)    ikpoint,k(ikpoint)%kptfrac,k(ikpoint)%kpt,                     &
     &     k(ikpoint)%n_neighbors,(k(ikpoint)%neighbor_kpts(innkpt,:),    &
     &           innkpt=1,k(ikpoint)%n_neighbors) ! DEBUG
      end do  ! ikpoint
      deallocate(innk)
      !
      ! end store nearest neighbors in k1,k2,k3 direction
      !
      ! begin calculate k gradient of E_n(k)
      !
      do ikpoint=1,nkpoints !
        allocate(k(ikpoint)%dE_dk(n:n,1:3,1:2)) ! band n only, 3 k directions, 2 spin channels
        k(ikpoint)%dE_dk=0.0d0
        do i1=1,3 ! cartesian directions
          nddk=0 ! count k-neighbors in this direction
          do innkpt=1,k(ikpoint)%n_neighbors
            jkpoint=k(ikpoint)%i_neighbors(innkpt)
            !call frac2abs(k(ikpoint)%kptfrac,gvecs,kcoords)               ! cartesian k-neighbor coords
            !call frac2abs(k(ikpoint)%neighbor_kptsfrac(innkpt,:),         &
!     &           gvecs,knn_coords)  ! cartesian k-neighbor coords
            kcoords=k(ikpoint)%kpt
            knn_coords=k(ikpoint)%neighbor_kpts(innkpt,:)
            dk=knn_coords(i1)-kcoords(i1)                                 ! k change in cartesian direction i1
            e_up=k(ikpoint)%e_up(n) ! energy eigenvalues at kpoint
            e_down=k(ikpoint)%e_down(n)
            E_nnk_up=k(jkpoint)%e_up(n) ! energy eigenvalues at neighbor kpoint
            E_nnk_down=k(jkpoint)%e_down(n)
            dE=(/E_nnk_up,E_nnk_down/)-(/e_up,e_down/) ! energy change 
            if (abs(dk).gt.tol3(i1)) then                                 ! consider only neighbors with substantial distance in i1 direction
              k(ikpoint)%dE_dk(n,i1,1:2)=k(ikpoint)%dE_dk(n,i1,1:2)       &
     &              +dE/dk   
              write(54,'("kpt ",I0,": acoords=",3(F10.6),", NN(x",        &
     &I0,")= kpt ",I0,": acoords=",3(F10.6),                              &
     &", dk(abs)=",3(F10.6))') ikpoint,kcoords,i1,jkpoint,knn_coords,     &
     &        knn_coords-kcoords           
              nddk=nddk+1
            end if  
          end do ! innkpt
          if(nddk.gt.0) k(ikpoint)%dE_dk(n,i1,1:2)                        &
     &       =k(ikpoint)%dE_dk(n,i1,1:2)/dble(nddk)
          write(54,'(8x,"kpt ",I5,", dir ",I0,": nddk=",I0)')             &
     &         ikpoint,i1,nddk
        end do ! i1  
      end do ! ikpoint
      !
      ! end calculate k gradient of E_n(k)
      !
      ! begin print file with kpoints, energies, and energy gradients for band n
      !
      write(filename,'("VASP_E_and_dE_dK_band_",I0)') n
      open(51,file=filename,status='replace')
      write(51,'("# ikpoint, fract. kpoint, abs. kpoint, E for spin up/d  &
     &own, dE/dK for spin up/down for band ",I0,", k-weight, distance tr  &
     &avelled along k-path")') n
      s=0.0d0
      do ikpoint=1,nkpoints
        if (ikpoint.gt.1) s=s+norm2(k(ikpoint)%kpt(1:3)                   &
     &     -k(ikpoint-1)%kpt(1:3))
        write(51,'(I5,1x,3(F12.6),4(1x,3(F12.6)),1x,F12.6)')              &
     &  ikpoint,k(ikpoint)%kptfrac,k(ikpoint)%kpt,                        &
     &  k(ikpoint)%e_up(n),k(ikpoint)%e_down(n),                          &
     &  k(ikpoint)%dE_dk(n,1:3,1),k(ikpoint)%dE_dk(n,1:3,2),              &
     &  k(ikpoint)%weight,s
      end do
      close(51)
      !
      ! end print file with kpoints, energies, and energy gradients for band n
      !
      ! begin sort wrt. k
      !
      ! sort with respect to k1 coordinate
      do ikpoint=1,nkpoints
        do jkpoint=ikpoint+1,nkpoints
          if (k(jkpoint)%kptfrac(1).lt.k(ikpoint)%kptfrac(1)
     &        -tol3(1)) then
                kstore=k(ikpoint)
                k(ikpoint)=k(jkpoint)
                k(jkpoint)=kstore
          end if
        end do
      end do
      ! sort with respect to k2 coordinate
      do ikpoint=1,nkpoints
        do jkpoint=ikpoint+1,nkpoints
          if (abs(k(jkpoint)%kptfrac(1)-k(ikpoint)%kptfrac(1))
     &          .le.tol3(1).and.k(jkpoint)%kptfrac(2).lt.
     &          k(ikpoint)%kptfrac(2)-tol3(2)) then
                  kstore=k(ikpoint)
                  k(ikpoint)=k(jkpoint)
                  k(jkpoint)=kstore
          end if
        end do
      end do
      ! sort with respect to k3 coordinate
      do ikpoint=1,nkpoints
        do jkpoint=ikpoint+1,nkpoints
          if (abs(k(jkpoint)%kptfrac(1)-k(ikpoint)%kptfrac(1))
     &        .le.tol3(1).and.abs(k(jkpoint)%kptfrac(2)
     &        -k(ikpoint)%kptfrac(2)).le.tol3(2).and.
     &        k(jkpoint)%kptfrac(3).lt.k(ikpoint)%kptfrac(3)
     &        -tol3(3)) then
                kstore=k(ikpoint)
                k(ikpoint)=k(jkpoint)
                k(jkpoint)=kstore
          end if
        end do
      end do
      !
      ! end sort wrt. k
      !
      ! begin print E_nk and ddk_E_nk on cut plane
      !
      do while (cutpos.lt.0.0d0) 
        cutpos=cutpos+1.0
      end do
      do while (cutpos.gt.1.0d0) 
        cutpos=cutpos-1.0
      end do
      if (cutdir==1) then
        print '(8x,"Printing E_n and dE_n/dK on plane spanned by G vecto  &
     &rs 2 and 3")'
      else if (cutdir==2) then
        print '(8x,"Printing E_n and dE_n/dK on plane spanned by G vecto  &
     &rs 1 and 3")'
      else if (cutdir==3) then
        print '(8x,"Printing E_n and dE_n/dK on plane spanned by G vecto  &
     &rs 1 and 2")'
      end if
      print '(8x,"with requested k1(frac)=",F9.6)',cutpos
      i0=ithenearest(cutpos,k(:)%kptfrac(1))
      cutposreal=k(i0)%kptfrac(cutdir)
      print '(8x,"with k1(frac) really=",F9.6)',                          &
     &       cutposreal
      ! write
      write(filename,'("VASP_E_and_dE_dK_band_",I0,"_k",I0,"_",           &
     &F6.4)') n,cutdir,cutpos
      open(52,file=filename,status='replace', err=102)  
      write(52,'("# fractional k, absolute k, E_n(k) for spin up/down, d  &
     &E(n)/dk for spin up, spin down, kpoint weight")') 
      ikpointlast=1
      do ikpoint=1,nkpoints
        if (abs(k(ikpoint)%kptfrac(cutdir)-cutposreal).lt.tol) then ! if kpoint is more or less on the plane
          select case (cutdir)
          case(1)
          if(k(ikpointlast)%kptfrac(2).gt.k(ikpoint)%kptfrac(2)+tol3(2)   &
     &   .or.k(ikpointlast)%kptfrac(3).gt.k(ikpoint)%kptfrac(3)+tol3(3))  & ! are we at the end of a data block ? 
     &    then
            write(52,'("")') ! for gnuplot, need to insert blank lines after each data block
          end if ! are we at the end of a data block ?
          case(2)
          if(k(ikpointlast)%kptfrac(1).gt.k(ikpoint)%kptfrac(1)+tol3(1)   &
     &   .or.k(ikpointlast)%kptfrac(3).gt.k(ikpoint)%kptfrac(3)+tol3(3))  & ! are we at the end of a data block ? 
     &    then
            write(52,'("")') ! for gnuplot, need to insert blank lines after each data block
          end if ! are we at the end of a data block ?
          case(3)
          if(k(ikpointlast)%kptfrac(1).gt.k(ikpoint)%kptfrac(1)+tol3(1)   &
     &   .or.k(ikpointlast)%kptfrac(2).gt.k(ikpoint)%kptfrac(2)+tol3(2))  & ! are we at the end of a data block ? 
     &    then
            write(52,'("")') ! for gnuplot, need to insert blank lines after each data block
          end if ! are we at the end of a data block ?
          end select
          write(52,'(2(1x,3F12.6),2(1x,F12.6),2(1x,3F12.6),1x,F12.6)')    &
     &         k(ikpoint)%kptfrac(1:3),k(ikpoint)%kpt(1:3),               &
     &         k(ikpoint)%e_up(n),k(ikpoint)%e_down(n),                   &
!     &         ddk_E_n(ikpoint,1:3,1),ddk_E_n(ikpoint,1:3,2),             &
     &         k(ikpoint)%dE_dk(n,1:3,1),k(ikpoint)%dE_dk(n,1:3,2),       &
     &         k(ikpoint)%weight       
          ikpointlast=ikpoint
          !print '(8x,"ikpointlast=",I0)',ikpointlast ! DEBUG
        end if ! kpoint is more or less on the plane
      end do ! ikpoint
      close(52)
      !
      ! end print E_nk and ddk_E_nk on cut plane
      !
      ! end normally:
      print fsubendext,'vasp_dE_dk'
      close(54)
      return

 100  close(51)
      n=n+1
      close(54)
      call error_stop('problem with OUTCAR file')
 101  close(51)
      n=n+1
      close(54)
      call error_stop('problem with EIGENVAL file')
 102  close(52)
      n=n+1
      close(54)
      call error_stop('problem with writing')
      !
      end subroutine vasp_dE_dk

c---------------------------------------------------------------------
    
      subroutine check_if_perovskite(infile,informat,nndist,              &
     &               is_perovskite)
      use defs
      use misc
      use readcoords
      implicit none
      character(len=*), intent(in) :: infile,informat
      logical is_perovskite
      type(atom),allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms,nspecies,iatom,j,ibond,nbonds,nangles,nbonds2
      integer nangles_180,nangles_90,iangle,iatom2
      double precision vecs(3,3),coordmean
      integer coordmin, coordmax 
      integer n_A, n_B, n_X
      double precision nndist,a_pseudocubic_half
      double precision tol_octahedral_angle,tol_lattice_distortion
      double precision tol_lattice_angle,tol_octahedral_tilt
      double precision B_B_dist_min,B_B_dist_max,angle
      double precision distvec(1:3),distvec2(1:3)
      character X_element*2
      logical good_octahedral_angles,almost_cubic_lattice,has_B,has_X
      logical good_tilt_angles,stoichio
      !
      ! begin initialize
      !
      print fsubstart,'check_if_perovskite'
      is_perovskite=.false.  
      !nndist=nndist_def
      !nndist=2.8d0
      tol_octahedral_angle=20.0
      tol_octahedral_tilt=30.0
      tol_lattice_distortion=0.20d0 
      tol_lattice_angle=10.0d0 
      print '(8x,"using ",F10.6," as cutoff")',nndist
      print '(8x,"using ",F10.6," as octahedral angle tolerance")',       &
     &        tol_octahedral_angle
      print '(8x,"using ",F10.6," as octahedral tilt tolerance")',        &
     &        tol_octahedral_tilt
      print '(8x,"using ",F10.6," as lattice distortion tolerance")',     &
     &        tol_lattice_distortion
      print '(8x,"using ",F10.6," as lattice angle tolerance")',          &
     &        tol_lattice_angle
      !
      ! end initialize
      !
      ! begin find octahedral atom
      !
      call read_coords(infile,informat,atoms,natoms,species,
     &                       nspecies,vecs)
      call coord(atoms,vecs,nndist,coordmin,coordmax,                     &
     &        coordmean)
      !
      a_pseudocubic_half=0.0d0
      do iatom=1,natoms
        atoms(iatom)%label="A" 
        !if (atoms(iatom)%coord==6) then ! OLD
        if (atoms(iatom)%coord==6) then                                 
          if(all(abs(atoms(iatom)%bondangles-90.0d0)                      &
     &    .lt.tol_octahedral_angle.or.abs(atoms(iatom)%bondangles         &
     &    -180.0d0).lt.tol_octahedral_angle)) then ! NEW
            atoms(iatom)%label="B" 
            !print*,atoms(iatom)%name,atoms(iatom)%label
            X_element=atoms(iatom)%neighbors(1)
            do ibond=1,atoms(iatom)%coord
              if (atoms(iatom)%bondlengths(ibond).gt.a_pseudocubic_half)  &
     &          a_pseudocubic_half=atoms(iatom)%bondlengths(ibond)
            end do
          end if
        end if
      end do ! iatom 
      ! use as cutoff for nearest neighbor the largest octahedral
      ! distance plus 10%
      a_pseudocubic_half=a_pseudocubic_half*1.1d0
      !
      !
      n_A=0
      n_B=0
      n_X=0
      do iatom=1,natoms
        if (atoms(iatom)%name(1:2)==X_element(1:2)) then
          atoms(iatom)%label="X" 
        end if
        print '(8x,"Atom ",A2," has type ",A2)', atoms(iatom)%name,       &
     &      atoms(iatom)%label
        if (atoms(iatom)%label=="A") n_A=n_A+1
        if (atoms(iatom)%label=="B") n_B=n_B+1
        if (atoms(iatom)%label=="X") n_X=n_X+1
      end do ! iatom 
      !
      ! end find octahedral atoms
      !
      ! begin get coordination again, this time with optimized cutoff
      !
      nndist=a_pseudocubic_half
      call coord(atoms,vecs,nndist,coordmin,coordmax,                     &
     &        coordmean)
      !
      ! end get coordination again, this time with optimized cutoff
      !
      ! begin check if octahedra are close to ideal (angles near 90 and
      ! 180 degree)
      !
      good_octahedral_angles=.true.
      good_tilt_angles=.true.
      has_B=.false.
      has_X=.false.
      do iatom=1,natoms
        if (atoms(iatom)%label=="B") then
          has_B=.true.
          nangles_90=0
          nangles_180=0
          nangles=size(atoms(iatom)%bondangles)
          do iangle=1,nangles
            if (abs(atoms(iatom)%bondangles(iangle)-180.0d0).lt.          &
     &          tol_octahedral_angle) nangles_180=nangles_180+1
            if (abs(atoms(iatom)%bondangles(iangle)- 90.0d0).lt.          &
     &          tol_octahedral_angle) nangles_90=nangles_90+1
          end do ! iangle
          if (nangles_90.ne.12.or.nangles_180.ne.3)                       &
     &        good_octahedral_angles=.false.
        end if
        if (atoms(iatom)%label=="X") then
          has_X=.true.
          nangles_180=0
          nangles=size(atoms(iatom)%bondangles)
          do iangle=1,nangles
            if (abs(atoms(iatom)%bondangles(iangle)-180.0d0).lt.          &
     &          tol_octahedral_tilt) nangles_180=nangles_180+1
!            print*,'atom,nangles,angle=',atoms(iatom)%name,nangles,       &
!     &          atoms(iatom)%bondangles(iangle)
!            print*,nangles,nangles_180
          end do ! iangle
          !if ((nangles.ne.1).or.(nangles_180.ne.1))                       &
          if ((nangles.lt.1).or.(nangles_180.lt.1))                       &
     &        good_tilt_angles=.false.
        end if
      end do ! iatom
      print '(8x,"good_octahedral_angles:",L1)',good_octahedral_angles
      print '(8x,"good_tilt_angles:",L1)',good_tilt_angles
      !
      ! end check if octahedra are close to ideal (angles near 90 and
      ! 180 degree)
      !
      !
      ! begin check if lattice parameters have correct angles and size
      ! ratios
      !
      almost_cubic_lattice=.true.
      B_B_dist_max=0.0d0
      B_B_dist_min=1000000.0d0
      do iatom=1,natoms
        if (atoms(iatom)%label=="X") then
          nbonds=atoms(iatom)%coord
          !if (nbonds.ne.2) then
          !  almost_cubic_lattice=.false.
          !  print '(8x,"found X atom with number of bonds/=2")'
          !else
            distvec=atoms(iatom)%neighborcoords(2,1:3)                    &
     &        -atoms(iatom)%neighborcoords(1,1:3)
            if (norm2(distvec).gt.B_B_dist_max)                           &
     &          B_B_dist_max=norm2(distvec)
            if (norm2(distvec).lt.B_B_dist_min)                           &
     &          B_B_dist_min=norm2(distvec)
          !end if
        end if ! B atom
      end do ! iatom
      print '(8x,"max,min B-B dist: ",F10.6,", ",F10.6)' ,B_B_dist_max,   &
     &        B_B_dist_min
      if(2.0d0*(B_B_dist_max-B_B_dist_min)/(B_B_dist_max+B_B_dist_min)    &
     &    .gt.tol_lattice_distortion) almost_cubic_lattice=.false.
      !
      ! angles:
      do iatom=1,natoms
        if (atoms(iatom)%label=="X") then
          nbonds=atoms(iatom)%coord
          do iatom2=iatom+1,natoms
            if (atoms(iatom2)%label=="X") then
              nbonds2=atoms(iatom2)%coord
              if (nbonds.eq.2.and.nbonds2.eq.2) then
                distvec=atoms(iatom)%neighborcoords(2,1:3)                &
     &           -atoms(iatom)%neighborcoords(1,1:3)
                distvec2=atoms(iatom2)%neighborcoords(2,1:3)              &
     &           -atoms(iatom2)%neighborcoords(1,1:3)
                angle=dot_product(distvec,distvec2)/                      &
     &               (norm2(distvec)*norm2(distvec2))
                if (angle.gt.1.0d0.and.angle.lt.1.0d0+1.0E-6)             &
     &            angle=1.0d0
                if (angle.lt.-1.0d0.and.angle.gt.-1.0d0-1.0E-6)           &
     &            angle=-1.0d0
                angle=acos(angle)*180.0d0/Pi
                if (abs(angle-90.0d0).gt.tol_lattice_angle.and.           &
     &              abs(angle-0.0d0).gt.tol_lattice_angle .and.           &
     &              abs(angle-180.0d0).gt.tol_lattice_angle) then      
                    almost_cubic_lattice=.false.
!                  print '(8x,"Found lattice angle of ",F10.6,             &
!     &                    " degree")',angle
                end if
              end if ! 2 bonds each
            end if ! X atom
          end do ! iatom2
        end if ! X atom
      end do ! iatom
      print '(8x,"almost_cubic_lattice:",L1)', almost_cubic_lattice 
      !
      ! end check if lattice parameters have correct angles and size
      ! ratios
      !
      ! begin check if stoichiometry is ~ ABX3
      !
      stoichio=.true.
      if ( n_B==0 .or. n_X==0) then
        stoichio=.false.
      else
        if (dble(n_X) .lt. 2.5d0*dble(n_B)) stoichio=.false.
        if (dble(n_X) .gt. 3.5d0*dble(n_B)) stoichio=.false.
      end if
      !
      ! end check if stoichiometry is A_yBX_2.5 <= stoichio <= A_yBX_3.5 
      !
      ! check if all requirements for perovskites are fulfilled:
      if (good_octahedral_angles .and. almost_cubic_lattice.and.has_B     &
     &    .and. has_X.and.good_tilt_angles.and. stoichio)                 &
     &     is_perovskite=.true.
      !
      ! begin print
      ! 
      print '(8x,"min, max, av coord. num.: ",2(I6),F10.6)',
     &        coordmin,coordmax,coordmean 
      print '(8x,"N_A=",I0)', n_A
      print '(8x,"N_B=",I0)', n_B
      print '(8x,"N_X=",I0)', n_X
      print '(8x,"~ stoichiometric: ",L1)', stoichio
      
      open(51,file="COORDINATION",status='replace')
      do iatom=1,natoms
        write(51,*) atoms(iatom)%name,atoms(iatom)%coord
      end do ! iatom
      close(51)
      open(51,file="COORDINATION_DETAILED",status='replace')
      do iatom=1,natoms
       if (1.lt.atoms(iatom)%coord) then
        write(51,*) atoms(iatom)%name,atoms(iatom)%coord,                 &
     &   (atoms(iatom)%neighbors(j),' ',j=1,atoms(iatom)%coord),          &
     &   (atoms(iatom)%bondlengths(j),' ',j=1,atoms(iatom)%coord),        &
     &    size(atoms(iatom)%bondangles),                                  &
     &   (atoms(iatom)%bondangles(j),' ',j=1,                             &
     &    size(atoms(iatom)%bondangles))
       else
        if (1.eq.atoms(iatom)%coord) then 
         write(51,*) atoms(iatom)%name,atoms(iatom)%coord,                &
     &   (atoms(iatom)%neighbors(j),' ',j=1,atoms(iatom)%coord),          &
     &   (atoms(iatom)%bondlengths(j),' ',j=1,atoms(iatom)%coord)
        else
         write(51,*) atoms(iatom)%name,atoms(iatom)%coord            
        end if
       end if
      end do ! iatom
      close(51)
      ! begin print coordinates of neighbors 
      open(51,file="COORDINATION_NEIGHBOR_COORDS",                        &
     &     status='replace')
      do iatom=1,natoms
       if (1.lt.atoms(iatom)%coord) then
         write(51,*) atoms(iatom)%name,atoms(iatom)%where,                &
     &   atoms(iatom)%abswhere,atoms(iatom)%coord,                        &
     &   (atoms(iatom)%neighbors(j),' ',j=1,atoms(iatom)%coord),          &
     &   (atoms(iatom)%neighborcoords(j,1:3),' ',                         &
     &  j=1,atoms(iatom)%coord)        
       else
        if (1.eq.atoms(iatom)%coord) then 
          write(51,*) atoms(iatom)%name,atoms(iatom)%where,               &
     &   atoms(iatom)%abswhere,atoms(iatom)%coord,                        &
     &   (atoms(iatom)%neighbors(j),' ',j=1,atoms(iatom)%coord),          &
     &   (atoms(iatom)%neighborcoords(j,1:3),' ',                         &
     &   j=1,atoms(iatom)%coord)
        else
         write(51,*) atoms(iatom)%name,atoms(iatom)%coord         
        end if
       end if
      end do ! iatom
      close(51)
      ! end print coordinates of neighbors 
      !
      ! end print
      ! 
      !
      print fsubendext,'check_if_perovskite'
      return
      !
      end subroutine check_if_perovskite

c---------------------------------------------------------------------
    
      subroutine check_if_perovskite_loose(infile,informat,nndist,        &
     &               is_perovskite)
      use defs
      use misc
      use readcoords
      implicit none
      character(len=*), intent(in) :: infile,informat
      logical is_perovskite
      type(atom),allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms,nspecies,iatom,j,ibond,nbonds,nangles,nbonds2
      integer nangles_180,nangles_90,iangle,iatom2
      integer n_A, n_B, n_X
      double precision vecs(3,3),coordmean
      integer coordmin, coordmax 
      double precision nndist,a_pseudocubic_half
      double precision tol_octahedral_angle,tol_lattice_distortion
      double precision tol_lattice_angle,tol_octahedral_tilt
      double precision B_B_dist_min,B_B_dist_max,angle
      double precision distvec(1:3),distvec2(1:3)
      character X_element*2
      logical good_octahedral_angles,almost_cubic_lattice,has_B,has_X
      logical good_tilt_angles, stoichio
      !
      ! begin initialize
      !
      print fsubstart,'check_if_perovskite_loose'
      is_perovskite=.false.  
      tol_octahedral_angle=30.0d0
      tol_octahedral_tilt=30.0d0
      tol_lattice_distortion=1.0d0/3.0d0 
      tol_lattice_angle=15.0d0 
      print '(8x,"using ",F10.6," as cutoff")',nndist
      print '(8x,"using ",F10.6," as octahedral angle tolerance")',       &
     &        tol_octahedral_angle
      print '(8x,"using ",F10.6," as octahedral tilt tolerance")',        &
     &        tol_octahedral_tilt
      print '(8x,"using ",F10.6," as lattice distortion tolerance")',     &
     &        tol_lattice_distortion
      print '(8x,"using ",F10.6," as lattice angle tolerance")',          &
     &        tol_lattice_angle
      !
      ! end initialize
      !
      ! begin find octahedral atom
      !
      call read_coords(infile,informat,atoms,natoms,species,
     &                       nspecies,vecs)
      call coord(atoms,vecs,nndist,coordmin,coordmax,                     &
     &        coordmean)
      !
      a_pseudocubic_half=0.0d0
      do iatom=1,natoms
        atoms(iatom)%label="A" 
        !if (atoms(iatom)%coord==6) then ! OLD
        if (atoms(iatom)%coord==6) then                                 
          if(all(abs(atoms(iatom)%bondangles-90.0d0)                      &
     &    .lt.tol_octahedral_angle.or.abs(atoms(iatom)%bondangles         &
     &    -180.0d0).lt.tol_octahedral_angle)) then ! NEW
            atoms(iatom)%label="B" 
            !print*,atoms(iatom)%name,atoms(iatom)%label
            X_element=atoms(iatom)%neighbors(1)
            do ibond=1,atoms(iatom)%coord
              if (atoms(iatom)%bondlengths(ibond).gt.a_pseudocubic_half)  &
     &          a_pseudocubic_half=atoms(iatom)%bondlengths(ibond)
            end do
          end if
        end if
      end do ! iatom 
      ! use as cutoff for nearest neighbor the largest octahedral
      ! distance plus 10%
      a_pseudocubic_half=a_pseudocubic_half*1.1d0
      !
!      !
!      do iatom=1,natoms
!        if (atoms(iatom)%name(1:2)==X_element(1:2)) then
!          atoms(iatom)%label="X" 
!        end if
!        print '(8x,"Atom ",A2," has type ",A2)', atoms(iatom)%name,       &
!     &      atoms(iatom)%label
!      end do ! iatom 
      !
      n_A=0
      n_B=0
      n_X=0
      do iatom=1,natoms
        if (atoms(iatom)%name(1:2)==X_element(1:2)) then
          atoms(iatom)%label="X" 
        end if
        print '(8x,"Atom ",A2," has type ",A2)', atoms(iatom)%name,       &
     &      atoms(iatom)%label
        if (atoms(iatom)%label=="A") n_A=n_A+1
        if (atoms(iatom)%label=="B") n_B=n_B+1
        if (atoms(iatom)%label=="X") n_X=n_X+1
      end do ! iatom 
      !
      ! end find octahedral atoms
      !
      ! end find octahedral atoms
      !
      ! begin get coordination again, this time with optimized cutoff
      !
      !nndist=a_pseudocubic_half
      nndist=nndist
      call coord(atoms,vecs,nndist,coordmin,coordmax,                     &
     &        coordmean)
      !
      ! end get coordination again, this time with optimized cutoff
      !
      ! begin check if octahedra are close to ideal (angles near 90 and
      ! 180 degree)
      !
      good_octahedral_angles=.true.
      good_tilt_angles=.true.
      has_B=.false.
      has_X=.false.
      do iatom=1,natoms
        if (atoms(iatom)%label=="B") then
          has_B=.true.
          nangles_90=0
          nangles_180=0
          nangles=size(atoms(iatom)%bondangles)
          do iangle=1,nangles
            if (abs(atoms(iatom)%bondangles(iangle)-180.0d0).lt.          &
     &          tol_octahedral_angle) nangles_180=nangles_180+1
            if (abs(atoms(iatom)%bondangles(iangle)- 90.0d0).lt.          &
     &          tol_octahedral_angle) nangles_90=nangles_90+1
          end do ! iangle
          if (nangles_90.ne.12.or.nangles_180.ne.3)                       &
     &        good_octahedral_angles=.false.
        end if
        if (atoms(iatom)%label=="X") then
          has_X=.true.
          nangles_180=0
          nangles=size(atoms(iatom)%bondangles)
          do iangle=1,nangles
            if (abs(atoms(iatom)%bondangles(iangle)-180.0d0).lt.          &
     &          tol_octahedral_tilt) nangles_180=nangles_180+1
!            print*,'atom,nangles,angle=',atoms(iatom)%name,nangles,       &
!     &          atoms(iatom)%bondangles(iangle)
!            print*,nangles,nangles_180
          end do ! iangle
          !if ((nangles.ne.1).or.(nangles_180.ne.1))                       &
          if ((nangles.lt.1).or.(nangles_180.lt.1))                       &
     &        good_tilt_angles=.false.
        end if
      end do ! iatom
      print '(8x,"good_octahedral_angles:",L1)',good_octahedral_angles
      print '(8x,"good_tilt_angles:",L1)',good_tilt_angles
      !
      ! end check if octahedra are close to ideal (angles near 90 and
      ! 180 degree)
      !
      !
      ! begin check if lattice parameters have correct angles and size
      ! ratios
      !
      almost_cubic_lattice=.true.
      B_B_dist_max=0.0d0
      B_B_dist_min=1000000.0d0
      do iatom=1,natoms
        if (atoms(iatom)%label=="X") then
          nbonds=atoms(iatom)%coord
          !if (nbonds.ne.2) then
          !  almost_cubic_lattice=.false.
          !  print '(8x,"found X atom with number of bonds/=2")'
          !else
            distvec=atoms(iatom)%neighborcoords(2,1:3)                    &
     &        -atoms(iatom)%neighborcoords(1,1:3)
            if (norm2(distvec).gt.B_B_dist_max)                           &
     &          B_B_dist_max=norm2(distvec)
            if (norm2(distvec).lt.B_B_dist_min)                           &
     &          B_B_dist_min=norm2(distvec)
          !end if
        end if ! B atom
      end do ! iatom
      print '(8x,"max,min B-B dist: ",F10.6,", ",F10.6)' ,B_B_dist_max,   &
     &        B_B_dist_min
      if(2.0d0*(B_B_dist_max-B_B_dist_min)/(B_B_dist_max+B_B_dist_min)    &
     &    .gt.tol_lattice_distortion) almost_cubic_lattice=.false.
      !
      ! angles:
      do iatom=1,natoms
        if (atoms(iatom)%label=="X") then
          nbonds=atoms(iatom)%coord
          do iatom2=iatom+1,natoms
            if (atoms(iatom2)%label=="X") then
              nbonds2=atoms(iatom2)%coord
              if (nbonds.eq.2.and.nbonds2.eq.2) then
                distvec=atoms(iatom)%neighborcoords(2,1:3)                &
     &           -atoms(iatom)%neighborcoords(1,1:3)
                distvec2=atoms(iatom2)%neighborcoords(2,1:3)              &
     &           -atoms(iatom2)%neighborcoords(1,1:3)
                angle=dot_product(distvec,distvec2)/                      &
     &               (norm2(distvec)*norm2(distvec2))
                if (angle.gt.1.0d0.and.angle.lt.1.0d0+1.0E-6)             &
     &            angle=1.0d0
                if (angle.lt.-1.0d0.and.angle.gt.-1.0d0-1.0E-6)           &
     &            angle=-1.0d0
                angle=acos(angle)*180.0d0/Pi
                if (abs(angle-90.0d0).gt.tol_lattice_angle.and.           &
     &              abs(angle-0.0d0).gt.tol_lattice_angle .and.           &
     &              abs(angle-180.0d0).gt.tol_lattice_angle) then      
                    almost_cubic_lattice=.false.
!                  print '(8x,"Found lattice angle of ",F10.6,             &
!     &                    " degree")',angle
                end if
              end if ! 2 bonds each
            end if ! X atom
          end do ! iatom2
        end if ! X atom
        if (atoms(iatom)%label=="A") n_A=n_A+1
        if (atoms(iatom)%label=="B") n_B=n_B+1
        if (atoms(iatom)%label=="X") n_X=n_X+1
      end do ! iatom
      print '(8x,"almost_cubic_lattice:",L1)', almost_cubic_lattice 
      !
      ! end check if lattice parameters have correct angles and size
      ! ratios
      !
      ! begin check if stoichiometry is ~ ABX3
      !
      stoichio=.true.
      if ( n_B==0 .or. n_X==0) then
        stoichio=.false.
      else
        if (dble(n_X) .lt. 2.5d0*dble(n_B)) stoichio=.false.
        if (dble(n_X) .gt. 3.5d0*dble(n_B)) stoichio=.false.
      end if
      !
      ! end check if stoichiometry is A_yBX_2.5 <= stoichio <= A_yBX_3.5 
      !
      ! check if all requirements for perovskites are fulfilled:
      if (good_octahedral_angles .and. almost_cubic_lattice.and.has_B     &
     &    .and. has_X.and.good_tilt_angles .and. stoichio)                &
     &    is_perovskite=.true.
      !
      ! begin print
      ! 
      print '(8x,"min, max, av coord. num.: ",2(I6),F10.6)',
     &        coordmin,coordmax,coordmean 
      print '(8x,"N_A=",I0)', n_A
      print '(8x,"N_B=",I0)', n_B
      print '(8x,"N_X=",I0)', n_X
      print '(8x,"~ stoichiometric: ",L1)', stoichio
      
      open(51,file="COORDINATION",status='replace')
      do iatom=1,natoms
        write(51,*) atoms(iatom)%name,atoms(iatom)%coord
      end do ! iatom
      close(51)
      open(51,file="COORDINATION_DETAILED",status='replace')
      do iatom=1,natoms
       if (1.lt.atoms(iatom)%coord) then
        write(51,*) atoms(iatom)%name,atoms(iatom)%coord,                 &
     &   (atoms(iatom)%neighbors(j),' ',j=1,atoms(iatom)%coord),          &
     &   (atoms(iatom)%bondlengths(j),' ',j=1,atoms(iatom)%coord),        &
     &    size(atoms(iatom)%bondangles),                                  &
     &   (atoms(iatom)%bondangles(j),' ',j=1,                             &
     &    size(atoms(iatom)%bondangles))
       else
        if (1.eq.atoms(iatom)%coord) then 
         write(51,*) atoms(iatom)%name,atoms(iatom)%coord,                &
     &   (atoms(iatom)%neighbors(j),' ',j=1,atoms(iatom)%coord),          &
     &   (atoms(iatom)%bondlengths(j),' ',j=1,atoms(iatom)%coord)
        else
         write(51,*) atoms(iatom)%name,atoms(iatom)%coord            
        end if
       end if
      end do ! iatom
      close(51)
      ! begin print coordinates of neighbors 
      open(51,file="COORDINATION_NEIGHBOR_COORDS",                        &
     &     status='replace')
      do iatom=1,natoms
       if (1.lt.atoms(iatom)%coord) then
         write(51,*) atoms(iatom)%name,atoms(iatom)%where,                &
     &   atoms(iatom)%abswhere,atoms(iatom)%coord,                        &
     &   (atoms(iatom)%neighbors(j),' ',j=1,atoms(iatom)%coord),          &
     &   (atoms(iatom)%neighborcoords(j,1:3),' ',                         &
     &  j=1,atoms(iatom)%coord)        
       else
        if (1.eq.atoms(iatom)%coord) then 
          write(51,*) atoms(iatom)%name,atoms(iatom)%where,               &
     &   atoms(iatom)%abswhere,atoms(iatom)%coord,                        &
     &   (atoms(iatom)%neighbors(j),' ',j=1,atoms(iatom)%coord),          &
     &   (atoms(iatom)%neighborcoords(j,1:3),' ',                         &
     &   j=1,atoms(iatom)%coord)
        else
         write(51,*) atoms(iatom)%name,atoms(iatom)%coord         
        end if
       end if
      end do ! iatom
      close(51)
      ! end print coordinates of neighbors 
      !
      ! end print
      ! 
      !
      print fsubendext,'check_if_perovskite_loose'
      return
      !
      end subroutine check_if_perovskite_loose

c---------------------------------------------------------------------
    
      subroutine check_if_perovskite_loose_analyze(infile,informat,axes,  &
     &                      nndist,is_perovskite)
      use defs
      use misc
      use readcoords
      use writecoords

      implicit none
      character(len=*), intent(in) :: infile,informat
      logical is_perovskite
      type(atom),allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms,nspecies,iatom,j,ibond,nbonds,nangles,nbonds2
      integer nangles_180,nangles_90,iangle,iatom2,iatom3
      integer n_A, n_B, n_X
      integer coordmin, coordmax 
      integer k,l,m,iatom4,n_plus_x,n_minus_x,n_plus_y,n_minus_y
      integer n_plus_z,n_minus_z,n_A_neighbors,i_A_neighbors
      integer n_B_neighbors,i_B_neighbors
      double precision vecs(3,3),coordmean,axes(3,3)
      double precision nndist,a_pseudocubic_half
      double precision tol_octahedral_angle,tol_lattice_distortion
      double precision tol_lattice_angle,tol_octahedral_tilt
      double precision B_B_dist_min,B_B_dist_max,angle
      double precision distvec(1:3),distvec2(1:3)
      double precision locpol(1:3),dist_AB,coords(1:3),vol,totalpol(1:3)
      double precision dist_BB,dist_AA
      double precision axis(1:3),angle_a,angle_b,angle_c,dist_AX
      double precision octa_angles(1:3),proj(1:3),octa_tilt_vector(1:3)
      double precision octa_axes(1:3,1:3),prod1,prod2,prod3
      double precision cross_prod(1:3),dot_prod
      double precision avdist_plus_x,avdist_minus_x,avdist_plus_y
      double precision avdist_minus_y,avdist_plus_z,avdist_minus_z
      double precision, allocatable:: A_neighbors(:,:),B_neighbors(:,:)
      character X_element*2,string*256
      logical good_octahedral_angles,almost_cubic_lattice,has_B,has_X
      logical good_tilt_angles, stoichio
      !
      ! begin initialize
      !
      print fsubstart,'check_if_perovskite_loose_analyze'
      is_perovskite=.false.  
      tol_octahedral_angle=30.0d0
      tol_octahedral_tilt=30.0d0
      !tol_octahedral_tilt=45.0d0
      !tol_lattice_distortion=1.0d0/3.0d0 
      tol_lattice_distortion=1.1d0/3.0d0 
      !tol_lattice_angle=15.0d0 
      tol_lattice_angle=20.0d0 
      print '(8x,"using ",F10.6," as cutoff")',nndist
      print '(8x,"using ",F10.6," as octahedral angle tolerance")',       &
     &        tol_octahedral_angle
      print '(8x,"using ",F10.6," as octahedral tilt tolerance")',        &
     &        tol_octahedral_tilt
      print '(8x,"using ",F10.6," as lattice distortion tolerance")',     &
     &        tol_lattice_distortion
      print '(8x,"using ",F10.6," as lattice angle tolerance")',          &
     &        tol_lattice_angle
      !
      ! end initialize
      !
      ! begin find octahedral atom
      !
      call read_coords(infile,informat,atoms,natoms,species,
     &                       nspecies,vecs)
      call coord(atoms,vecs,nndist,coordmin,coordmax,                     &
     &        coordmean)
      !
      a_pseudocubic_half=0.0d0
      do iatom=1,natoms
        atoms(iatom)%label="A" 
        !if (atoms(iatom)%coord==6) then ! OLD
        if (atoms(iatom)%coord==6) then                                 
          if(all(abs(atoms(iatom)%bondangles-90.0d0)                      &
     &    .lt.tol_octahedral_angle.or.abs(atoms(iatom)%bondangles         &
     &    -180.0d0).lt.tol_octahedral_angle)) then ! NEW
            atoms(iatom)%label="B" 
            !print*,atoms(iatom)%name,atoms(iatom)%label
            X_element=atoms(iatom)%neighbors(1)
            do ibond=1,atoms(iatom)%coord
              if (atoms(iatom)%bondlengths(ibond).gt.a_pseudocubic_half)  &
     &          a_pseudocubic_half=atoms(iatom)%bondlengths(ibond)
            end do
          end if
        end if
      end do ! iatom 
      ! use as cutoff for nearest neighbor the largest octahedral
      ! distance plus 10%
      a_pseudocubic_half=a_pseudocubic_half*1.1d0
      !
      !
      n_A=0
      n_B=0
      n_X=0
      do iatom=1,natoms
        if (atoms(iatom)%name(1:2)==X_element(1:2)) then
          atoms(iatom)%label="X" 
        end if
        print '(8x,"Atom ",A2," has type ",A2)', atoms(iatom)%name,       &
     &      atoms(iatom)%label
        if (atoms(iatom)%label=="A") n_A=n_A+1
        if (atoms(iatom)%label=="B") n_B=n_B+1
        if (atoms(iatom)%label=="X") n_X=n_X+1
      end do ! iatom 
      !
      ! end find octahedral atoms
      !
      ! end find octahedral atoms
      !
      ! begin get coordination again, this time with optimized cutoff
      !
      !nndist=a_pseudocubic_half
      nndist=nndist
      call coord(atoms,vecs,nndist,coordmin,coordmax,                     &
     &        coordmean)
      !
      ! end get coordination again, this time with optimized cutoff
      !
      ! begin check if octahedra are close to ideal (angles near 90 and
      ! 180 degree)
      !
      good_octahedral_angles=.true.
      good_tilt_angles=.true.
      has_B=.false.
      has_X=.false.
      do iatom=1,natoms
        if (atoms(iatom)%label=="B") then
          has_B=.true.
          nangles_90=0
          nangles_180=0
          nangles=size(atoms(iatom)%bondangles)
          do iangle=1,nangles
            if (abs(atoms(iatom)%bondangles(iangle)-180.0d0).lt.          &
     &          tol_octahedral_angle) nangles_180=nangles_180+1
            if (abs(atoms(iatom)%bondangles(iangle)- 90.0d0).lt.          &
     &          tol_octahedral_angle) nangles_90=nangles_90+1
          end do ! iangle
          if (nangles_90.ne.12.or.nangles_180.ne.3)                       &
     &        good_octahedral_angles=.false.
        end if
        if (atoms(iatom)%label=="X") then
          has_X=.true.
          nangles_180=0
          nangles=size(atoms(iatom)%bondangles)
          do iangle=1,nangles
            if (abs(atoms(iatom)%bondangles(iangle)-180.0d0).lt.          &
     &          tol_octahedral_tilt) nangles_180=nangles_180+1
          end do ! iangle
          if ((nangles.lt.1).or.(nangles_180.lt.1))                       &
     &        good_tilt_angles=.false.
        end if
      end do ! iatom
      print '(8x,"good_octahedral_angles:",L1)',good_octahedral_angles
      print '(8x,"good_tilt_angles:",L1)',good_tilt_angles
      !
      ! end check if octahedra are close to ideal (angles near 90 and
      ! 180 degree)
      !
      !
      ! begin check if lattice parameters have correct angles and size
      ! ratios
      !
      almost_cubic_lattice=.true.
      B_B_dist_max=0.0d0
      B_B_dist_min=1000000.0d0
      do iatom=1,natoms
        if (atoms(iatom)%label=="X") then
          nbonds=atoms(iatom)%coord
            distvec=atoms(iatom)%neighborcoords(2,1:3)                    &
     &        -atoms(iatom)%neighborcoords(1,1:3)
            if (norm2(distvec).gt.B_B_dist_max)                           &
     &          B_B_dist_max=norm2(distvec)
            if (norm2(distvec).lt.B_B_dist_min)                           &
     &          B_B_dist_min=norm2(distvec)
        end if ! B atom
      end do ! iatom
      print '(8x,"max,min B-B dist: ",F10.6,", ",F10.6)' ,B_B_dist_max,   &
     &        B_B_dist_min
      if(2.0d0*(B_B_dist_max-B_B_dist_min)/(B_B_dist_max+B_B_dist_min)    &
     &    .gt.tol_lattice_distortion) almost_cubic_lattice=.false.
      !
      ! angles:
      do iatom=1,natoms
        if (atoms(iatom)%label=="X") then
          nbonds=atoms(iatom)%coord
          do iatom2=iatom+1,natoms
            if (atoms(iatom2)%label=="X") then
              nbonds2=atoms(iatom2)%coord
              if (nbonds.eq.2.and.nbonds2.eq.2) then
                distvec=atoms(iatom)%neighborcoords(2,1:3)                &
     &           -atoms(iatom)%neighborcoords(1,1:3)
                distvec2=atoms(iatom2)%neighborcoords(2,1:3)              &
     &           -atoms(iatom2)%neighborcoords(1,1:3)
                angle=dot_product(distvec,distvec2)/                      &
     &               (norm2(distvec)*norm2(distvec2))
                if (angle.gt.1.0d0.and.angle.lt.1.0d0+1.0E-6)             &
     &            angle=1.0d0
                if (angle.lt.-1.0d0.and.angle.gt.-1.0d0-1.0E-6)           &
     &            angle=-1.0d0
                angle=acos(angle)*180.0d0/Pi
                if (abs(angle-90.0d0).gt.tol_lattice_angle.and.           &
     &              abs(angle-0.0d0).gt.tol_lattice_angle .and.           &
     &              abs(angle-180.0d0).gt.tol_lattice_angle) then      
                    almost_cubic_lattice=.false.
                end if
              end if ! 2 bonds each
            end if ! X atom
          end do ! iatom2
        end if ! X atom
      end do ! iatom
      print '(8x,"almost_cubic_lattice:",L1)', almost_cubic_lattice 
      !
      ! end check if lattice parameters have correct angles and size
      ! ratios
      !
      ! begin check if stoichiometry is ~ ABX3
      !
      stoichio=.true.
      if ( n_B==0 .or. n_X==0) then
        stoichio=.false.
      else
        if (dble(n_X) .lt. 2.5d0*dble(n_B)) stoichio=.false.
        if (dble(n_X) .gt. 3.5d0*dble(n_B)) stoichio=.false.
      end if
      !
      ! end check if stoichiometry is A_yBX_2.5 <= stoichio <= A_yBX_3.5 
      !
      ! check if all requirements for perovskites are fulfilled:
      if (good_octahedral_angles .and. almost_cubic_lattice.and.has_B     &
     &    .and. has_X.and.good_tilt_angles .and. stoichio)                &
     &    is_perovskite=.true.
      !
      ! begin print
      ! 
      print '(8x,"min, max, av coord. num.: ",2(I6),F10.6)',
     &        coordmin,coordmax,coordmean 
      print '(8x,"N_A=",I0)', n_A
      print '(8x,"N_B=",I0)', n_B
      print '(8x,"N_X=",I0)', n_X
      print '(8x,"~ good octahedral angles: ",L1)',                       &
     &         good_octahedral_angles
      print '(8x,"~ almost cubic lattice: ",L1)',                         &
     &         almost_cubic_lattice
      print '(8x,"~ has B: ",L1)', has_B
      print '(8x,"~ has X: ",L1)', has_X
      print '(8x,"~ stoichiometric: ",L1)', stoichio
      print '(8x,"~ good tilt angles: ",L1)', good_tilt_angles
      
      open(51,file="COORDINATION",status='replace')
      do iatom=1,natoms
        write(51,*) atoms(iatom)%name,atoms(iatom)%coord
      end do ! iatom
      close(51)
      open(51,file="COORDINATION_DETAILED",status='replace')
      do iatom=1,natoms
       if (1.lt.atoms(iatom)%coord) then
        write(51,*) atoms(iatom)%name,atoms(iatom)%coord,                 &
     &   (atoms(iatom)%neighbors(j),' ',j=1,atoms(iatom)%coord),          &
     &   (atoms(iatom)%bondlengths(j),' ',j=1,atoms(iatom)%coord),        &
     &    size(atoms(iatom)%bondangles),                                  &
     &   (atoms(iatom)%bondangles(j),' ',j=1,                             &
     &    size(atoms(iatom)%bondangles))
       else
        if (1.eq.atoms(iatom)%coord) then 
         write(51,*) atoms(iatom)%name,atoms(iatom)%coord,                &
     &   (atoms(iatom)%neighbors(j),' ',j=1,atoms(iatom)%coord),          &
     &   (atoms(iatom)%bondlengths(j),' ',j=1,atoms(iatom)%coord)
        else
         write(51,*) atoms(iatom)%name,atoms(iatom)%coord            
        end if
       end if
      end do ! iatom
      close(51)
      ! begin print coordinates of neighbors 
      open(51,file="COORDINATION_NEIGHBOR_COORDS",                        &
     &     status='replace')
      do iatom=1,natoms
       if (1.lt.atoms(iatom)%coord) then
         write(51,*) atoms(iatom)%name,atoms(iatom)%where,                &
     &   atoms(iatom)%abswhere,atoms(iatom)%coord,                        &
     &   (atoms(iatom)%neighbors(j),' ',j=1,atoms(iatom)%coord),          &
     &   (atoms(iatom)%neighborcoords(j,1:3),' ',                         &
     &  j=1,atoms(iatom)%coord)        
       else
        if (1.eq.atoms(iatom)%coord) then 
          write(51,*) atoms(iatom)%name,atoms(iatom)%where,               &
     &   atoms(iatom)%abswhere,atoms(iatom)%coord,                        &
     &   (atoms(iatom)%neighbors(j),' ',j=1,atoms(iatom)%coord),          &
     &   (atoms(iatom)%neighborcoords(j,1:3),' ',                         &
     &   j=1,atoms(iatom)%coord)
        else
         write(51,*) atoms(iatom)%name,atoms(iatom)%coord         
        end if
       end if
      end do ! iatom
      close(51)
      ! end print coordinates of neighbors 
      !
      ! begin set G-type spin on B atoms
      !
      atoms(:)%spin=0.0d0
      !
      ! begin set first B atom spin up
      !
      iatom=1
10    if (atoms(iatom)%label=="B") then
        atoms(iatom)%spin=1.0d0
        goto 20
      end if
      iatom=iatom+1
      goto 10
20    continue
      !
      ! end set first B atom spin up
      !
      do while (any(atoms(:)%spin==0.0d0))
        do iatom=1,natoms
          if (.not.atoms(iatom)%label=="B") atoms(iatom)%spin=2.0d0
          if (atoms(iatom)%label=="B") then
            !
            ! begin collect coordinates of B neighbors of B
            !
            ! begin find B atoms
            n_B_neighbors=0
            do iatom2=1,natoms
              if (atoms(iatom2)%label=="B".and.iatom2.ne.iatom) then
                do k=-1,1
                do l=-1,1
                do m=-1,1
                coords=atoms(iatom2)%abswhere+dble(k)*vecs(1,:)           &
     &            +dble(l)*vecs(2,:)+dble(m)*vecs(3,:)                 
                distvec=coords-atoms(iatom)%abswhere
                dist_BB=norm2(distvec)
                !if (dist_BB.le.B_B_dist_max*1.01d0) then
                if (dist_BB.le.nndist*2.0d0) then
                  n_B_neighbors=n_B_neighbors+1 ! count neighbors of B
                  if (atoms(iatom2)%spin.ne.0.0d0) then
                    atoms(iatom)%spin=-atoms(iatom2)%spin
                  end if
                end if ! dist_BB.le.nndist*...
                end do ! m
                end do ! l
                end do ! k
              end if !atoms(iatom2)%label=="B".and.iatom2.ne.iatom
            end do ! iatom2
            if (n_B_neighbors.ne.6) then
              string=""
              write(string,'("B atom has wrong number of neighbors (",    &
     &              I0,")")') n_B_neighbors
              call warning(adjustl(trim(string)))
            end if
            !
            ! end collect coordinates of B neighbors of B
            !
          end if ! label=B
        end do ! iatoms
      end do ! while any spin not set  
      do iatom=1,natoms
        if(atoms(iatom)%label.ne."B") atoms(iatom)%spin=0.0d0
      end do
      !
      ! write
      call write_coords("COORDS_AND_SPINS","poscar_spin",atoms,natoms,    &
     &                  species,nspecies,vecs)
      !
      ! end set G-type spin on B atoms
      !
      if (all(axes==0.0d0)) axes=vecs ! if no axes specified, calculate octahedral tilts wrt. cell vectors
      !
      ! begin determine axis automatically from A lattice
      !
      if (all(axes==-1.0d0)) then 
        ! 
        ! begin find first A atom
        !
        iatom=1
        do while (atoms(iatom)%label.ne."A") 
          iatom=iatom+1
        end do
        ! 
        ! end find first A atom
        !
        ! begin search for closest 6 A neighbors
        !
        print '(8x,"Looking for nearest A-neighbors of A...")'
        if(allocated(A_neighbors)) deallocate(A_neighbors)
        allocate(A_neighbors(6,1:3))
        i_A_neighbors=1
        do iatom2=1,natoms
          if (atoms(iatom2)%label=="A".and.iatom2.ne.iatom) then
            do k=-1,1
            do l=-1,1
            do m=-1,1
            coords=atoms(iatom2)%abswhere+dble(k)*vecs(1,:)             &
     &        +dble(l)*vecs(2,:)+dble(m)*vecs(3,:)                   
            distvec=coords-atoms(iatom)%abswhere
            dist_AA=norm2(distvec)
            !print*, dist_AA
            if (dist_AA.le.nndist*2.0d0) then
              if(i_A_neighbors.le.6)                                      &
     &             A_neighbors(i_A_neighbors,1:3)=coords
              i_A_neighbors=i_A_neighbors+1
            end if ! dist_AA.le.nndist*...
            end do ! m
            end do ! l
            end do ! k
          end if !atoms(iatom2)%label=="A".and.iatom2.ne.iatom
        end do ! iatom2
        print '(8x,"...done.")'
        print '(8x,"Found ",I0," A-neighbors within ", F6.3, " Angs.")',  &
     &          i_A_neighbors-1,nndist*2.0d0
        !
        ! end search for closest 6 A neighbors
        !
        ! begin determine axes
        !
        axes(1,1:3)=A_neighbors(1,1:3)-atoms(iatom)%abswhere(1:3) ! first axis
        k=2
        dot_prod=1.0d0
        do while (dot_prod.gt.0.1d0.or.dot_prod.lt.-0.1d0)
          k=k+1
          distvec=A_neighbors(k,1:3)-atoms(iatom)%abswhere(1:3)
          dot_prod=dot_product(distvec,axes(1,1:3))/(norm2(distvec)       &
     &             *norm2(axes(k,1:3)))
        end do
        axes(2,1:3)=distvec ! 2nd axis
        axes(3,1:3)=cross_product(axes(1,:),axes(2,:)) ! 3rd axis
        !
        ! end determine axes
        !
      end if ! (all(axes==-1.0D0))
      !
      ! end determine axis automatically from A lattice
      !
      !
      ! begin determine axis automatically from B lattice
      !
      if (all(axes==-2.0d0)) then 
        ! 
        ! begin find first B atom
        !
        iatom=1
        do while (atoms(iatom)%label.ne."B") 
          iatom=iatom+1
        end do
        ! 
        ! end find first B atom
        !
        ! begin search for closest 6 B neighbors
        !
        print '(8x,"Looking for nearest B-neighbors of B...")'
        if(allocated(B_neighbors)) deallocate(B_neighbors)
        allocate(B_neighbors(6,1:3))
        i_B_neighbors=1
        do iatom2=1,natoms
          if (atoms(iatom2)%label=="B".and.iatom2.ne.iatom) then
            do k=-1,1
            do l=-1,1
            do m=-1,1
            coords=atoms(iatom2)%abswhere+dble(k)*vecs(1,:)             &
     &        +dble(l)*vecs(2,:)+dble(m)*vecs(3,:)                   
            distvec=coords-atoms(iatom)%abswhere
            dist_BB=norm2(distvec)
            !print*, dist_BB
            if (dist_BB.le.nndist*2.0d0) then
              if(i_B_neighbors.le.6)                                      &
     &             B_neighbors(i_B_neighbors,1:3)=coords
              i_B_neighbors=i_B_neighbors+1
            end if ! dist_BB.le.nndist*...
            end do ! m
            end do ! l
            end do ! k
          end if !atoms(iatom2)%label=="B".and.iatom2.ne.iatom
        end do ! iatom2
        print '(8x,"...done.")'
        print '(8x,"Found ",I0," B-neighbors within ", F6.3, " Angs.")',  &
     &          i_B_neighbors-1,nndist*2.0d0
        !
        ! end search for closest 6 B neighbors
        !
        ! begin determine axes
        !
        axes(1,1:3)=B_neighbors(1,1:3)-atoms(iatom)%abswhere(1:3) ! first axis
        k=2
        dot_prod=1.0d0
        do while (dot_prod.gt.0.1d0.or.dot_prod.lt.-0.1d0)
          k=k+1
          distvec=B_neighbors(k,1:3)-atoms(iatom)%abswhere(1:3)
          dot_prod=dot_product(distvec,axes(1,1:3))/(norm2(distvec)       &
     &             *norm2(axes(k,1:3)))
        end do
        axes(2,1:3)=distvec ! 2nd axis
        axes(3,1:3)=cross_product(axes(1,:),axes(2,:)) ! 3rd axis
        !
        ! end determine axes
        !
      end if ! (all(axes==-1.0D0))
      !
      ! end determine axis automatically from A lattice
      !
      axes(1,:)=axes(1,:)/norm2(axes(1,:))
      axes(2,:)=axes(2,:)/norm2(axes(2,:))
      axes(3,:)=axes(3,:)/norm2(axes(3,:))
      print '(8x,"using as reference axes",/,3(8x,3(F12.6)))',            &
     &       (axes(k,:),k=1,3)
      !
      ! begin print polarization, angles for each unit cell
      !
      open(51,file="PEROVSKITE_DETAILS",status='replace')
      open(52,file="POLARIZATION_DETAILS",status='replace')
      open(53,file="OCTAHEDRA_DETAILS",status='replace')
      write(53,'("# B site:  x, y z (Angs), 1. octahedral axis, tilt ang  &
     &les around 1., 2., 3. axis, spin, octahedral tilt vector")')
      write(53,'("# B site:  x, y z (Angs), 2. octahedral axis, tilt ang  &
     &les around 1., 2., 3. axis, spin, octahedral tilt vector")')
      write(53,'("# B site:  x, y z (Angs), 3. octahedral axis, tilt ang  &
     &les around 1., 2., 3. axis, spin, octahedral tilt vector, B pos. (  &
     &frac, neighbor coords (abs))")')
      call vecs2vol(vecs,vol)
      totalpol=0.0d0
      do iatom=1,natoms
        locpol=0.0d0
        if (atoms(iatom)%label=="B") then
          write(51,'(A2,1x,6(F12.6))') atoms(iatom)%name(1:2),            &
     &     atoms(iatom)%where, atoms(iatom)%abswhere
          write(51,*) atoms(iatom)%ineighbors(:)
          ! begin add polarization contribution of B atoms 
          locpol=locpol+atoms(iatom)%abswhere*atoms(iatom)%charge
          ! end add polarization contribution of B atoms 
          ibond=1
          octa_axes=0.0d0
          do iatom2=1,atoms(iatom)%coord
            iatom3=atoms(iatom)%ineighbors(iatom2)
            ! begin add polarization contribution of X atoms 
            locpol=locpol+atoms(iatom)%neighborcoords(iatom2,:)           &
     &        *atoms(iatom3)%charge/2.0d0
            ! end add polarization contribution of X atoms 
            write(51,'(A2,1x,12(F12.6))') atoms(iatom3)%name(1:2),        &
     &       atoms(iatom3)%where, atoms(iatom3)%abswhere,                 &
     &       atoms(iatom)%neighborcoords(iatom2,1:3),                     &
     &       atoms(iatom)%neighborcoords(iatom2,1:3)                      &
     &       -atoms(iatom)%abswhere
            ! begin analyze octahedron
            do iatom4=iatom2+1,atoms(iatom)%coord
              ! begin look for 180 degree angle (octahedral axis)
              angle=atoms(iatom)%bondangles(ibond)
              write(51,'(" O2: ",6(F12.6))')                              &
     &          atoms(iatom)%neighborcoords(iatom4,1:3),angle
              if (abs(angle-180.0d0).le.tol_lattice_angle) then
                axis=atoms(iatom)%neighborcoords(iatom4,1:3)              &
     &              -atoms(iatom)%neighborcoords(iatom2,1:3)           
                axis=axis/norm2(axis)
                !
                ! begin sort axes
                !
                prod1=abs(dot_product(axis,axes(1,:)))
                prod2=abs(dot_product(axis,axes(2,:)))
                prod3=abs(dot_product(axis,axes(3,:)))
                if (prod1.gt.prod2.and.prod1.gt.prod3) then            
                        octa_axes(1,:)=axis(:)
                        if (dot_product(octa_axes(1,:),axes(1,:)).lt.     &
     &                      0.0d0) octa_axes(1,:)=-octa_axes(1,:)
                end if
                if (prod2.gt.prod1.and.prod2.gt.prod3) then 
                        octa_axes(2,:)=axis(:)
                        if (dot_product(octa_axes(2,:),axes(2,:)).lt.     &
     &                      0.0d0) octa_axes(2,:)=-octa_axes(2,:)
                end if
                if (prod3.gt.prod1.and.prod3.gt.prod2) then 
                        octa_axes(3,:)=axis(:)
                        if (dot_product(octa_axes(3,:),axes(3,:)).lt.     &
     &                      0.0d0) octa_axes(3,:)=-octa_axes(3,:)
                end if
                !
                ! end sort axes
                !
              end if  
              ! end look for 180 degree angle (octahedral axis)
              ibond=ibond+1
            end do ! iatom4
            ! end analyze octahedron
          end do ! iatom2
          !
          ! begin print axes and rotations
          !
          octa_angles=0.0d0
          do k=1,3 ! axis index
!            print*, "k=",k 
            !  
            ! first axis  
            !  
            l=mod(k,3)+1 
!            print*, "l1=",l 
            proj(:)=octa_axes(l,:)-dot_product(octa_axes(l,:),axes(k,:))  &
     &        * axes(k,:)
            if (norm2(proj).gt.0.0d0) proj=proj/(norm2(proj))
            cross_prod=cross_product(axes(l,:),proj(:))
!            print*, "axis, octa_axis,cross_prod=",axes(l,:),              &
!     &          proj(:), cross_prod
            dot_prod=dot_product(cross_prod(:),axes(k,:))
            octa_angles(k)=octa_angles(k)+asin(norm2(cross_prod))         &
     &                     *sign(1.0d0,dot_prod)*180.0d0/Pi
            !  
            ! second axis  
            !  
            l=mod(l,3)+1
!            print*, "l2=",l
            proj=octa_axes(l,:)-dot_product(octa_axes(l,:),axes(k,:))     &
     &        * axes(k,:)
            if (norm2(proj).gt.0.0d0) proj=proj/(norm2(proj))
            cross_prod=cross_product(axes(l,:),proj(:))
            dot_prod=dot_product(cross_prod(:),axes(k,:))
            octa_angles(k)=octa_angles(k)+asin(norm2(cross_prod))         &
     &                     *sign(1.0d0,dot_prod)*180.0d0/Pi
!            print*, "axis, octa_axis,cross_prod=",axes(l,:),              &
!     &          octa_axes(l,:), cross_prod
!              print*, "dot_prod,sign=",dot_prod, sign(1.0d0,dot_prod)     &
            !  
            ! second axis  
            !  
            octa_angles(k)=octa_angles(k)/2.0d0
            !print*," "
            octa_tilt_vector=octa_angles(1)*axes(1,:)+octa_angles(2)*     &
     &               axes(2,:)+octa_angles(3)*axes(3,:)                
          end do ! k (axis index
          do k=1,3
            write(51,'("oct. axis: ",3(F12.6),8x,3(F12.6))')              &
     &        octa_axes(k,:),octa_angles(1:3)
              write(53,'(A2," site:",3(F12.6)," oct. axis:",3(F12.6),     &
     &              1x,"angles:",3(F12.6),F8.3,5x,3(F8.3),126(F12.6))')   &
     &               atoms(iatom)%name(1:2),                              &
     &               atoms(iatom)%abswhere,                               &
     &               octa_axes(k,:) ,octa_angles(1:3),atoms(iatom)%spin,  &
     &               octa_tilt_vector(1:3),atoms(iatom)%where,            &
     &               (atoms(iatom)%neighborcoords(iatom2,1:3),iatom2=1,   &
     &               atoms(iatom)%coord)             
          end do ! k (axis index
          !
          ! end print axes and rotations
          !
          ! begin find A atoms
          do iatom2=1,natoms
            if (atoms(iatom2)%label=="A") then
              do k=-1,1
              do l=-1,1
              do m=-1,1
              coords=atoms(iatom2)%abswhere+dble(k)*vecs(1,:)             &
     &          +dble(l)*vecs(2,:)+dble(m)*vecs(3,:)                   
              dist_AB=norm2(coords-atoms(iatom)%abswhere)
              if (dist_AB.le.nndist*sqrt(3.0d0)*1.25d0) then ! totally arbitrary threshold that seems to work for BFO in PBE
                write(51,'(A2,1x,6(F12.6))') atoms(iatom2)%name(1:2),     &
     &          coords,coords-atoms(iatom)%abswhere
                ! begin add polarization contribution of A atoms 
                locpol=locpol+coords*atoms(iatom2)%charge/8.0d0
                ! end add polarization contribution of A atoms 
              end if ! dist_AB.le.nndist*sqrt(3.0d0)
              end do ! m
              end do ! l
              end do ! k
            end if !atoms(iatom2)%label=="A"
          end do ! iatom2
          ! end find A atoms
          locpol=locpol*ec*1.0E22/(vol/dble(n_B))
          totalpol=totalpol+locpol/dble(n_B)
          write(51,'(3(F12.6))')  locpol
          write(52,'(A2,1x,9(F12.6))') atoms(iatom)%name(1:2),            &
     &     atoms(iatom)%where, atoms(iatom)%abswhere, locpol
        end if !atoms(iatom)%label=="B" 
      end do ! iatom  
      write(51,'("total polarization in uC/cm²: ",3(F12.6))') totalpol
      print '(8x,"total polarization in uC/cm²: ",3(F12.6))', totalpol
      close(51)
      close(52)
      close(53)
      !
      ! end print polarization, angles for each unit cell
      ! 
      ! begin analyze A-site environment
      !
      open(51,file="ENVIRONMENT_OF_A",status="replace")
      open(52,file="ENVIRONMENT_OF_A_x",status="replace")
      open(53,file="ENVIRONMENT_OF_A_y",status="replace")
      open(54,file="ENVIRONMENT_OF_A_z",status="replace")
      open(55,file="ENVIRONMENT_OF_A_dist",status="replace")
      open(56,file="ENVIRONMENT_OF_A_x_hist",status="replace")
      open(57,file="ENVIRONMENT_OF_A_y_hist",status="replace")
      open(58,file="ENVIRONMENT_OF_A_z_hist",status="replace")
      open(59,file="ENVIRONMENT_OF_A_dist_hist",status="replace")
      do iatom=1,natoms
        if (atoms(iatom)%label=="A") then
          avdist_plus_x=0.0d0
          avdist_minus_x=0.0d0
          avdist_plus_y=0.0d0
          avdist_minus_y=0.0d0
          avdist_plus_z=0.0d0
          avdist_minus_z=0.0d0
          n_plus_x=0
          n_minus_x=0
          n_plus_y=0
          n_minus_y=0
          n_plus_z=0
          n_minus_z=0
          ! begin find X atoms
          if(allocated(A_neighbors)) deallocate(A_neighbors)
          n_A_neighbors=0
          do iatom2=1,natoms
            if (atoms(iatom2)%label=="X") then
              do k=-1,1
              do l=-1,1
              do m=-1,1
              coords=atoms(iatom2)%abswhere+dble(k)*vecs(1,:)             &
     &          +dble(l)*vecs(2,:)+dble(m)*vecs(3,:)                   
              distvec=coords-atoms(iatom)%abswhere
              dist_AX=norm2(distvec)
              !print*, dist_AX
              if (dist_AX.le.nndist*sqrt(2.0d0)) then
                n_A_neighbors=n_A_neighbors+1 ! count neighbors of A
              end if ! dist_AX.le.nndist*...
              end do ! m
              end do ! l
              end do ! k
            end if !atoms(iatom2)%label=="X"
          end do ! iatom2
          !
          ! begin collect coordinates of neighbors of A
          !
          allocate(A_neighbors(n_A_neighbors,1:3))
          i_A_neighbors=1
          do iatom2=1,natoms
            if (atoms(iatom2)%label=="X") then
              do k=-1,1
              do l=-1,1
              do m=-1,1
              coords=atoms(iatom2)%abswhere+dble(k)*vecs(1,:)             &
     &          +dble(l)*vecs(2,:)+dble(m)*vecs(3,:)                   
              distvec=coords-atoms(iatom)%abswhere
              dist_AX=norm2(distvec)
              if (dist_AX.le.nndist*sqrt(2.0d0)) then
                A_neighbors(i_A_neighbors,1:3)=coords
                i_A_neighbors=i_A_neighbors+1
              end if ! dist_AX.le.nndist*...
              end do ! m
              end do ! l
              end do ! k
            end if !atoms(iatom2)%label=="X"
          end do ! iatom2
          !
          ! end collect coordinates of neighbors of A
          !
          ! end find X atoms
          !
!          write(51,'(6(F12.6))') avdist_plus_x/dble(n_plus_x),            &
!     &    avdist_minus_x/dble(n_minus_x),avdist_plus_y/dble(n_plus_y),    &
!     &    avdist_minus_y/dble(n_minus_y),avdist_plus_z/dble(n_plus_z),    &
!     &    avdist_minus_z/dble(n_minus_z)
          write(51,'(A2,1x,6(F12.6),1x,6(F12.6),48(F12.6))')              &
     &    atoms(iatom)%name(1:2),atoms(iatom)%where,                      &
     &    atoms(iatom)%abswhere,                                          &
     &    (A_neighbors(i_A_neighbors,1:3)-atoms(iatom)%abswhere,          &
     &     i_A_neighbors=1,n_A_neighbors),                                &
     &    (norm2(A_neighbors(i_A_neighbors,1:3)-atoms(iatom)%abswhere),   &
     &     i_A_neighbors=1,n_A_neighbors)
          write(52,'(A2,1x,6(F12.6),1x,6(F12.6),48(F12.6))')              &
     &    atoms(iatom)%name(1:2),atoms(iatom)%where,                      &
     &    atoms(iatom)%abswhere,                                          &
     &    (A_neighbors(i_A_neighbors,1)-atoms(iatom)%abswhere(1),         &
     &     i_A_neighbors=1,n_A_neighbors)                               
          write(53,'(A2,1x,6(F12.6),1x,6(F12.6),48(F12.6))')              &
     &    atoms(iatom)%name(1:2),atoms(iatom)%where,                      &
     &    atoms(iatom)%abswhere,                                          &
     &    (A_neighbors(i_A_neighbors,2)-atoms(iatom)%abswhere(2),         &
     &     i_A_neighbors=1,n_A_neighbors)                               
          write(54,'(A2,1x,6(F12.6),1x,6(F12.6),48(F12.6))')              &
     &    atoms(iatom)%name(1:2),atoms(iatom)%where,                      &
     &    atoms(iatom)%abswhere,                                          &
     &    (A_neighbors(i_A_neighbors,3)-atoms(iatom)%abswhere(3),         &
     &     i_A_neighbors=1,n_A_neighbors)                               
          write(55,'(A2,1x,6(F12.6),1x,6(F12.6),48(F12.6))')              &
     &    atoms(iatom)%name(1:2),atoms(iatom)%where,                      &
     &    atoms(iatom)%abswhere,                                          &
     &    (norm2(A_neighbors(i_A_neighbors,1:3)-atoms(iatom)%abswhere),   &
     &     i_A_neighbors=1,n_A_neighbors)
          do i_A_neighbors=1,n_A_neighbors
            write(56,'(A2,1x,6(F12.6),1x,3(F12.6))')                      &
     &      atoms(iatom)%name(1:2),atoms(iatom)%where,                    &
     &      atoms(iatom)%abswhere,                                        &
     &      A_neighbors(i_A_neighbors,1)-atoms(iatom)%abswhere(1)    
            write(57,'(A2,1x,6(F12.6),1x,3(F12.6))')                      &
     &      atoms(iatom)%name(1:2),atoms(iatom)%where,                    &
     &      atoms(iatom)%abswhere,                                        &
     &      A_neighbors(i_A_neighbors,2)-atoms(iatom)%abswhere(2)    
            write(58,'(A2,1x,6(F12.6),1x,3(F12.6))')                      &
     &      atoms(iatom)%name(1:2),atoms(iatom)%where,                    &
     &      atoms(iatom)%abswhere,                                        &
     &      A_neighbors(i_A_neighbors,3)-atoms(iatom)%abswhere(3)    
            write(59,'(A2,1x,6(F12.6),1x,3(F12.6))')                      &
     &      atoms(iatom)%name(1:2),atoms(iatom)%where,                    &
     &      atoms(iatom)%abswhere,                                        &
     &      norm2(A_neighbors(i_A_neighbors,:)-atoms(iatom)%abswhere(:))
          end do  
          write(56,'(" ")')
          write(57,'(" ")')
          write(58,'(" ")')
          write(59,'(" ")')
          !
        end if ! label==A
      end do ! iatom
      close(51)
      close(52)
      close(53)
      close(54)
      close(55)
      close(56)
      close(57)
      close(58)
      close(59)
      !
      ! end analyze A-site environment
      !
      ! end print
      !
      print fsubendext,'check_if_perovskite_loose_analyze'
      return
      !
      end subroutine check_if_perovskite_loose_analyze

c---------------------------------------------------------------------
    
      subroutine get_pol(infile, informat)
      use defs
      use readcoords
      use misc, only: vecs2vol

      implicit none
      character(len=*) infile, informat
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms, nspecies
      double precision vecs(1:3,1:3)
      integer iatom,i
      character string*128
      double precision pol(1:3),vol

      print fsubstart, 'get_pol'
      call read_coords(infile,informat,atoms,natoms,species,
     &                       nspecies,vecs)
      
      ! test: write BEC to file
      open(51,file="POSCAR_BEC.vasp",status="replace")
      write(51,'("POSCAR containing BEC generated with cofima")')
      write(51,*) ' 1.0'
      write(51,'(3(F20.10))') vecs(1,1:3)
      write(51,'(3(F20.10))') vecs(2,1:3)
      write(51,'(3(F20.10))') vecs(3,1:3)
      string=''
      do i=1,nspecies
        write(string((i-1)*6+1:),'(A6)') trim(adjustl(species(i)%name))
      end do
      write(51,'(A128)') string
      string=''
      do i=1,nspecies
        write(string((i-1)*6+1:),'(I6)') species(i)%howmany
      end do
      write(51,'(A128)') string
      write(51,*) 'direct'
      do iatom=1,natoms
!        write(51, '(3(F20.10),5x,3(F20.10),5x,9(F12.6))')                 &
!     &         atoms(iatom)%where,atoms(iatom)%abswhere,
!     &        (atoms(iatom)%BEC(i,1:3),i=1,3)
        write(51, '(3(F20.10),5x,9(F12.6))')                              &
     &         atoms(iatom)%where,(atoms(iatom)%BEC(i,1:3),i=1,3)
      end do
      close(51)
      !
      ! begin calculate polarization as P=sum_i r_i Z_i, Z_i=BEC tensor
      !
      pol=0.0d0
      do iatom=1,natoms
        pol=pol+matmul(atoms(iatom)%BEC,atoms(iatom)%abswhere)
      end do
      call vecs2vol(vecs,vol)
      pol=pol*ec*1.0E22/(vol)
      !
      ! end calculate polarization as P=sum_i r_i Z_i, Z_i=BEC tensor
      !
      print '(8x,"electric polarization = ",3(F12.6))', pol
      !
      print fsubendext, 'get_pol'

      end subroutine get_pol

c---------------------------------------------------------------------

      subroutine vasp_read_WAVEDER()
      ! 
      ! reads binary vasp output file WAVEDER and prints content to a formatted
      ! file WAVEDER.dat
      !  
      use defs, only : fsubstart,fsubendext,error_stop
      !
      implicit none
      !
      logical file_found
      integer IU
      INTEGER, PARAMETER :: q =SELECTED_REAL_KIND(10)
      INTEGER, PARAMETER :: qs =SELECTED_REAL_KIND(5)
      !
      TYPE wavedes ! the wavedes type in vasp contains a lot more variables. Here we keep only the ones needed to read WAVEDER
         INTEGER NB_TOT                ! total number bands
         INTEGER NKPTS                 ! number of k-points in the irreducable wedge of the BZ (IBZ)
         INTEGER ISPIN                 ! number of spins
      END TYPE wavedes
      !
      TYPE (wavedes)     WDES
      REAL(q) :: NODES_IN_DIELECTRIC_FUNCTION
      INTEGER :: NBANDS_CDER        ! number of bands n considered 
      REAL(q) :: WPLASMON(3,3)
      !REAL(qs), ALLOCATABLE :: CDER_BETWEEN_STATES(:,:,:,:,:)
      COMPLEX(qs), ALLOCATABLE :: CDER_BETWEEN_STATES(:,:,:,:,:)
      INTEGER ISP,NB1,NB2,NK
      !
      print fsubstart,'vasp_read_WAVEDER'  
      !
      INQUIRE(file="WAVEDER", exist=file_found)
      if (.not.file_found) call error_stop('WAVEDER not found')
      ! 
      !
      ! begin read WAVEDER
      !
      IU=51
      OPEN(UNIT=IU,FILE='WAVEDER',FORM='UNFORMATTED',STATUS='UNKNOWN')
      READ(IU) WDES%NB_TOT, NBANDS_CDER, WDES%NKPTS, WDES%ISPIN
      ALLOCATE(CDER_BETWEEN_STATES(WDES%NB_TOT, NBANDS_CDER, WDES%NKPTS,  &
     &         WDES%ISPIN, 3))
      READ(IU) NODES_IN_DIELECTRIC_FUNCTION
      READ(IU) WPLASMON
      READ(IU) CDER_BETWEEN_STATES
      CLOSE(IU)
      !
      ! end read WAVEDER
      !
      ! begin write WAVEDER.dat
      !
      ! open file
      OPEN(IU,FILE='WAVEDER.dat',STATUS='REPLACE')
      ! write header information
      WRITE(IU,'(3I6," # ISPIN NKPTS CONSIDERED_BANDS")') WDES%ISPIN,     &
     &      WDES%NKPTS,NBANDS_CDER
      ! write CDER_BETWEEN_STATES
!      WRITE(IU,'("# ISP KPT IB1 IB2 WAVEDER(x) WAVEDER(y) WAVEDER(z)")')
      DO ISP=1,WDES%ISPIN
      DO NK=1,WDES%NKPTS
         DO NB1=1,NBANDS_CDER
         DO NB2=1,NBANDS_CDER
!           WRITE(IU,'(4I6,3(2F20.10))')                                   &
!     &       ISP,NK,NB1,NB2,CDER_BETWEEN_STATES(NB1,NB2,NK,ISP,1:3)
           WRITE(IU,'(2I6,3(2F20.10))')                                   &
     &       NB1,NB2,CDER_BETWEEN_STATES(NB1,NB2,NK,ISP,1:3)
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      ! close
      CLOSE(IU)
      !
      ! end write WAVEDER.dat
      !
      !
      print fsubendext,'vasp_read_WAVEDER'  
      !
      end subroutine vasp_read_WAVEDER

c---------------------------------------------------------------------

      subroutine vasp_eps2_from_WAVEDER(nomega,omega_min,omega_max,       &
     &           sigma,vbands,cbands,spins,kpoints, smearing)
      ! 
      ! reads binary vasp output file WAVEDER and prints content to a formatted
      ! file WAVEDER.dat
      !  
      use defs, only : fsubstart,fsubendext,error_stop,pi,ec,bohr,Ryd
      use misc, only : fermi_dist,delta_function_fermi
      use misc, only : gaussian_dist,delta_function_gaussian
      !
      implicit none
      !
      logical file_found
      integer fileunit
      integer, parameter :: q=selected_real_kind(10)
      integer,parameter :: qs=selected_real_kind(5)
      !
      integer nomega
      integer, allocatable :: vbands(:),cbands(:),spins(:),kpoints(:)
      double precision omega_min,omega_max,efermi,sigma,occ_v,occ_c,vol
      double precision,allocatable :: eigenvalues(:,:,:) ! energy eigenvalues for each kpoint, band, spin
      double precision, allocatable :: kweights(:)
      real :: nodes_in_epsilon
      integer :: nkpoints,nbands_tot,nspins,ncbands_eps     ! number of bands considered for epsilon
      real omega_plasmon(1:3,1:3)
      complex(qs), allocatable :: matrix_elements(:,:,:,:,:)
      double precision, allocatable :: omega(:)
      integer ispin,ik,ib1,ib2,iomega,iomega2,idir,jdir
      double precision, allocatable :: epsilon2(:,:,:) ! epsilon(direction1,direction2,omega)
      character*1024 smearing
      !
      print fsubstart,'vasp_eps2_from_WAVEDER'  
      !
      print '(8x,"Using ",I0," frequencies from ",F12.6," to ",F12.6,     &
     &        " eV")',nomega,omega_min,omega_max
      print '(/8x,"Using ",A10," smearing by ",F9.6," eV")', adjustl(     &
     &   trim(adjustl(smearing))),sigma
      !
      allocate(omega(1:nomega))
      do iomega=1,nomega
        omega(iomega)=omega_min+dble(iomega-1)*(omega_max-omega_min)      &
     &    /dble(nomega-1)
      end do
      !
      ! begin get volume
      !
      call vasp_read_volume(vol)
      !
      ! end get volume
      !
      ! begin read Fermi energy
      !
      call vasp_read_fermi_energy(efermi)
      efermi=efermi+0.02d0 ! move Fermi energy a bit above the VBM
      !
      ! end read Fermi energy
      !
      ! begin read WAVEDER
      !
      inquire(file="WAVEDER", exist=file_found)
      if (.not.file_found) call error_stop('WAVEDER not found')
      !
      fileunit=51
      open(unit=fileunit,file='WAVEDER',form='unformatted',               &
     &   status='old')
      read(fileunit) nbands_tot, ncbands_eps, nkpoints, nspins
      print '(8x,I0," bands in total, ",I0," cond. bands for epsilon, ",  &
     &        I0," kpoints, ",I0," spins")',nbands_tot,ncbands_eps,       &
     &        nkpoints,nspins
      if (any(vbands==-1)) then
        deallocate(vbands)
        allocate(vbands(nbands_tot))
        do ib1=1,nbands_tot
          vbands(ib1)=ib1
        end do
      end if
      print '(8x,"Using ",I0," valence bands")',size(vbands)
      if (any(cbands==-1)) then
        deallocate(cbands)
        allocate(cbands(nbands_tot))
        do ib1=1,nbands_tot
          cbands(ib1)=ib1
        end do
      end if
      print '(8x,"Using ",I0," conduction bands")',size(cbands)
      if (any(spins==-1)) then
        deallocate(spins)
        allocate(spins(nspins))
        do ispin=1,nspins
          spins(ispin)=ispin
        end do
      end if
      if (any(kpoints==-1)) then
        deallocate(kpoints)
        allocate(kpoints(nkpoints))
        do ik=1,nkpoints
          kpoints(ik)=ik
        end do
      end if
      print '(8x,"Using ",I0," spins")',size(spins)
      allocate(matrix_elements(nbands_tot, ncbands_eps, nkpoints,         &
     &         nspins, 3))
      matrix_elements=(0.0d0,0.0d0)
      read(fileunit) nodes_in_epsilon
      read(fileunit) omega_plasmon
      read(fileunit) matrix_elements(:,:,1:nkpoints,:,:)
      close(fileunit)
      !
      ! end read WAVEDER
      !
      ! begin read eigenvalues and weights
      !
      call vasp_read_eigenvalues_and_weights(eigenvalues,kweights)
      !
      ! end read eigenvalues and weights
      !
      ! begin write WAVEDER.dat
      !
      ! open file
      open(fileunit,FILE='WAVEDER_EXTENDED.dat',STATUS='REPLACE')
      ! write header information
      write(fileunit,'(3I6," # nspins nkpoints nbands_tot ncbands_eps")'  &
     &     ) nspins, nkpoints,nbands_tot,ncbands_eps
      ! write matrix elements
      write(fileunit,'("# ispin ik ib1 ib2 EV(ib1) EV(ib2) complex_matri  &
     &x_ele(x) complex_matrix_ele(y) complex_matrix_ele(z)")')
      do ispin=1,nspins
       do ik=1,nkpoints
        do ib1=1,nbands_tot
         do ib2=1,ncbands_eps
           write(fileunit,'(4I6,2(1x,F15.6),1x,3(2F20.10))')              &
     &        ispin,ik,ib1,ib2,eigenvalues(ik,ib1,ispin),                 &
     &        eigenvalues(ik,ib2,ispin),                                  &
     &        matrix_elements(ib1,ib2,ik,ispin,1:3)
         end do ! ib2
        end do ! ib1
       end do ! ik
      end do ! ispin
      ! close
      close(fileunit)
      !
      ! end write WAVEDER.dat
      !
      !
      ! begin calculate epsilon2
      !
      allocate(epsilon2(1:3,1:3,nomega))
      epsilon2=0.0d0
      !
      do idir=1,3
      do jdir=1,3
       do iomega=1,nomega
         do ispin=1,nspins
          if (.not.any(spins==ispin)) cycle
          do ik=1,nkpoints 
            if (.not.any(kpoints==ik)) cycle
           do ib1=1,nbands_tot ! valence bands
            if (.not.(any(vbands==ib1).or.any(cbands==ib1))) cycle
            do ib2=ib1+1,nbands_tot ! conduction bands
             if (.not.(any(vbands==ib2).or.any(cbands==ib2))) cycle
             select case(smearing)
             case('fermi')
             occ_v=fermi_dist(eigenvalues(ik,ib1,ispin)-efermi,sigma)
             occ_c=fermi_dist(eigenvalues(ik,ib2,ispin)-efermi,sigma)
             epsilon2(idir,jdir,iomega)=epsilon2(idir,jdir,iomega)         &
     &         +conjg(matrix_elements(ib2,ib1,ik,ispin,idir))              &
     &         *matrix_elements(ib2,ib1,ik,ispin,jdir)                     &
     &         *delta_function_fermi(eigenvalues(ik,ib2,ispin)             &
     &         -eigenvalues(ik,ib1,ispin)-omega(iomega),sigma)             &
     &         *kweights(ik)*(occ_v-occ_c)
               ! contributions from negative frequencies:
             epsilon2(idir,jdir,iomega)=epsilon2(idir,jdir,iomega)         &
     &         +conjg(matrix_elements(ib2,ib1,ik,ispin,idir))              &
     &         *matrix_elements(ib2,ib1,ik,ispin,jdir)                     &
     &         *delta_function_fermi(-eigenvalues(ik,ib2,ispin)            &
     &         +eigenvalues(ik,ib1,ispin)-omega(iomega),sigma)             &
     &         *kweights(ik)*(occ_v-occ_c)
             case('gaussian')
             occ_v=gaussian_dist(eigenvalues(ik,ib1,ispin)-efermi,sigma)
             occ_c=gaussian_dist(eigenvalues(ik,ib2,ispin)-efermi,sigma)
             epsilon2(idir,jdir,iomega)=epsilon2(idir,jdir,iomega)        &
     &         +conjg(matrix_elements(ib2,ib1,ik,ispin,idir))             &
     &         *matrix_elements(ib2,ib1,ik,ispin,jdir)                    &
     &         *delta_function_gaussian(eigenvalues(ik,ib2,ispin)         &
     &         -eigenvalues(ik,ib1,ispin)-omega(iomega),sigma)            &
     &         *kweights(ik)*(occ_v-occ_c)
               ! contributions from negative frequencies:
             epsilon2(idir,jdir,iomega)=epsilon2(idir,jdir,iomega)        &
     &         +conjg(matrix_elements(ib2,ib1,ik,ispin,idir))             &
     &         *matrix_elements(ib2,ib1,ik,ispin,jdir)                    &
     &         *delta_function_gaussian(-eigenvalues(ik,ib2,ispin)        &
     &         +eigenvalues(ik,ib1,ispin)-omega(iomega),sigma)            &
     &         *kweights(ik)*(occ_v-occ_c)
           case default
             goto 100
           end select
            end do ! ib2
           end do ! ib1
          end do ! ik
         end do ! ispin
         if (mod(iomega*10,nomega).lt.10) print'(8x,"i,j,omega=",x,I0,x,  &
     &     I0,x,F12.6)',idir,jdir,omega(iomega)
       end do !iomega
      end do ! jdir
      end do ! idir
      epsilon2=epsilon2*4.0d0*pi**2/vol*bohr*Ryd*2.0d0
      if (nspins==1) epsilon2=2.0d0*epsilon2
      !
      ! end calculate epsilon
      
      ! begin write epsilon
      !
      ! open file
      open(fileunit,FILE='EPS2_FROM_WAVEDER.DAT',STATUS='REPLACE')
      ! write header information
      write(fileunit,'(" # nspins nkpoints n_vbands n_cbands",4I6)')      &
     &       size(spins),nkpoints,size(vbands),size(cbands)
      ! write epsilon
      write(fileunit,'("# omega eps2(1,1) eps2(2,2), eps2(3,3), eps2(1,2  &
     &) eps2(2,3) eps2(1,3)")')
      do iomega=1,nomega
!           write(fileunit,'(F15.6,1x,6(e13.5))')                          &
!           write(fileunit,'(F15.6,1x,6(e16.6))')                          &
           write(fileunit,'(F15.6,1x,6(E15.6E3))')                        &
     &        omega(iomega),epsilon2(1,1,iomega),epsilon2(2,2,iomega),    &
     &        epsilon2(3,3,iomega),epsilon2(1,2,iomega),epsilon2(2,3,     &
     &        iomega),epsilon2(3,1,iomega)
      end do ! iomega
      ! close
      close(fileunit)
      !
      ! end write epsilon2
      !
      !
      print fsubendext,'vasp_eps2_from_WAVEDER'  
      return
      !
100   close(fileunit)
      call error_stop('unknown electronic smearing')      
      !
      end subroutine vasp_eps2_from_WAVEDER

c---------------------------------------------------------------------

      subroutine WAVECAR_ana()

      use defs, only : fsubstart,fsubendext,error_stop

      IMPLICIT NONE

      INTEGER, PARAMETER :: q =SELECTED_REAL_KIND(10)
      LOGICAL  :: WAVECAR_exists
      REAL(q)  :: RDUM, RISPIN
      integer  :: ICMPLX,IRECLW
      parameter(ICMPLX=16)

!      FROM INWAV:
!      USE prec
!      USE wave
!      USE lattice
!      IMPLICIT REAL(q) (A-H,O-Z)
      REAL(q) RKPTSF,RBANDF,ENMAXI
      integer I,J,NKPTSF,NBANDF
      TYPE latt
         REAL(q) :: SCALE
         REAL(q) :: A(3,3),B(3,3)
         REAL(q) :: ANORM(3),BNORM(3)
         REAL(q) :: OMEGA
!tb start
         REAL(q) AVEL(3,3)             ! lattice velocities
         INTEGER INITlatv              ! lattice velocities initialized                  !
!tb end
      END TYPE

      !TYPE (wavedes)  WDES
      TYPE (latt)     LATT_INI !,LATT_CUR
      !LOGICAL     LDIFF

      print fsubstart,'WAVECAR_ana'  

      INQUIRE(FILE='WAVECAR',EXIST=WAVECAR_exists)
      if (.not.WAVECAR_exists) call error_stop('WAVECAR file not found')

! first reopen with assumed (wrong) record length ICMPLX
      OPEN(UNIT=51,FILE='WAVECAR',ACCESS='DIRECT',                        &
     &         FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=ICMPLX)
! the first record contains the record length, get it ...
      RDUM=0._q
      READ(51,REC=1,ERR=17421) RDUM,RISPIN 
      ! print*,RDUM,RISPIN ! DEBUG
      IRECLW=NINT(RDUM)
! in the worst case IRECLW could contain completely useless data and useless is
! all <=0 or all >10000000 (since on output we use RECL=IO%ICMPLX*MAX(NRPLWV,6)
! or RECL=(NB_TOT+2)*ICMPLX more than ten millions sounds not very realistic)
      IF (IRECLW<=0) then
        close(51)
        call error_stop('RECL from WAVECAR makes no sense')
      end if
      GOTO 17422

17421 CONTINUE
      close(51)
      !print*,RDUM,RISPIN
      call error_stop('error reading RECL from WAVECAR')
17422 CONTINUE
      CLOSE(51)
      print'(8x,"Record length (RECL):",I30)',IRECLW
      print'(8x,"Number of spins (ISPIN):",I30)',NINT(RISPIN)
! reopen with correct record length (clumsy all that, I know ...)
      OPEN(UNIT=51,FILE='WAVECAR',ACCESS='DIRECT',                        &
     & FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=IRECLW)
!
!      The following is from VASP 
!      SUBROUTINE INWAV_HEAD(WDES, LATT_INI, LATT_CUR, ENMAXI, ISTART, IU0)

      READ(51,REC=2,ERR=100) RKPTSF,RBANDF,ENMAXI,                        &
     &                       ((LATT_INI%A(I,J),I=1,3),J=1,3)

      NKPTSF=NINT(RKPTSF)
      print '(8x,"Number of k-points (NKPTS):",I5)',NKPTSF
      NBANDF=NINT(RBANDF)
      print '(8x,"Number of bands (NBANDF):",I5)',NBANDF
      print '(8x,"Energy cutoff in eV (ENMAXI):",F20.6)',ENMAXI
      print '(8x,"Lattice vectors ine Angs (LATT_INI%A):",                &
     &     /,3(8x,3(F20.6),/))',LATT_INI%A


  100 CONTINUE
      close(51)

      print fsubendext,'WAVECAR_ana'

      end subroutine WAVECAR_ana

c---------------------------------------------------------------------

      subroutine vasp_write_BORN(natoms_BEC,the_atoms_BEC)
      ! 
      ! reads dielectric tensor and Born effective charges from OUTCAR
      ! file and writes them to file BORN (used by phonopy).
      !  
      use defs, only : fsubstart,fsubendext,error_stop, atom, element
      use readcoords
      !
      implicit none
      !
      logical file_found
      integer IU
      !
      double precision  epsilon_tensor(3,3)
      double precision vecs(1:3,1:3)
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer natoms,nspecies, iatom, i
      integer natoms_BEC
      integer, intent(in) :: the_atoms_BEC(:)
      logical found_epsilon, found_BEC
      character line*256
      !
      print fsubstart,'vasp_write_BORN'  
      !
      INQUIRE(file="OUTCAR", exist=file_found)
      if (.not.file_found) call error_stop('OUTCAR not found')
      !
      ! begin read BEC
      !
      call read_coords("OUTCAR","outcar",atoms,natoms,species,nspecies,   &
     &  vecs)
      !
      ! end read BEC
      ! 
      ! begin read epsilon from OUTCAR
      !
      IU=51
      OPEN(UNIT=IU,FILE='OUTCAR',FORM='FORMATTED',STATUS='OLD')
      rewind(IU)
      
      found_epsilon=.false.
      found_BEC=.false.
10    read (IU,'(A100)',end=100,err=101) line
      if (index(line,' MACROSCOPIC STATIC DIELECTRIC TENSOR (including l  &
     &ocal field effects in DFT)').gt.0) then
        found_epsilon=.true.
        read (IU,'(A100)',end=101,err=101) line
        do i=1,3
          read (IU,*,end=101,err=101) epsilon_tensor(i,1:3)
        end do ! i
        goto 10
      end if
      if (index(line,' BORN EFFECTIVE CHARGES (in e, cummulative output)  &
     &').gt.0) then
        found_BEC=.true.
      end if
      goto 10  
100   continue
      close(IU)
      if (.not.found_epsilon) goto 101
      if (.not.found_BEC) goto 101
      print fsubendext,'vasp_write_BORN'  
      !
      ! begin write BORN
      !
      open(IU, file="BORN", status="replace")
      !write(IU,'("default value")')
      write(IU,'("14.399652")')
      write(IU,'(9(F15.6))') epsilon_tensor(1,1:3),epsilon_tensor(2,:),   &
     &   epsilon_tensor(3,:)
      if (natoms_BEC.le.0) then
        do iatom=1,natoms
          write(IU,'(9(F15.6))') atoms(iatom)%BEC(1,:),                   &
     &    atoms(iatom)%BEC(2,:),atoms(iatom)%BEC(3,:)
        end do ! iatom
      else
        do iatom=1,natoms_BEC
          write(IU,'(9(F15.6))') atoms(the_atoms_BEC(iatom))%BEC(1,:),    &
     &    atoms(the_atoms_BEC(iatom))%BEC(2,:),                           &
     &    atoms(the_atoms_BEC(iatom))%BEC(3,:)
        end do ! iatom
      end if
      close(IU)
      !
      ! end write BORN
      !
      return

101   continue
      close(IU)
      call error_stop("could not read epsilon or BEC from OUTCAR") 
      !
      end subroutine vasp_write_BORN

!----------------------------------------------------      

      subroutine phonopy_print_bands(infile,outfile)
      !  
      use defs, only: fsubstart,fsubendext, error_stop
      implicit none
      !
      character(len=*) :: infile, outfile
      integer nqpoint,npath,natom
      integer, allocatable :: segment_nqpoint(:)
      character*64, allocatable :: labels(:,:)
      double precision, allocatable :: qpos(:,:),freq(:,:),distance(:)
      complex, allocatable :: eigenvec(:,:,:,:)
      !
      character*128 line,chardum
      integer begin_to_read, ipath,iband,iatom,iqpoint,segment_iqpoint,i
      logical has_eigenvec
      double precision re_ev,im_ev
      !
      print fsubstart,'phonopy_print_bands'  
      !
      ! read input file (band.yaml):
      !
      print '(8x,"Going to read phonon band structure from file ",A90)',  &
     &        adjustl(trim(adjustl(infile)))
      open(51,file=infile,status="old")
      rewind(51)      
      !
      ! check if eigenvectors are there:
      has_eigenvec=.false.
8     read(51,'(A128 )', end=9, err=100) line 
      if (index(line,'eigenvector').gt.0) then
        has_eigenvec=.true.
        goto 9
      end if
      goto 8

9     rewind(51)      
      print '(8x,"has eigenvectors: ",L1)',has_eigenvec 
      !
10    read(51,'(A128 )', end=20, err=100) line 
      ! read number of qpoints:
      if (index(line,'nqpoint:').gt.0 .and.index(line,'segment').le.0)    &
     &    then
        begin_to_read=index(line,'nqpoint:')+8
        read(line(begin_to_read:),*) nqpoint
        print '(8x,I10," qpoints")', nqpoint
        allocate(qpos(nqpoint,1:3),distance(nqpoint))
      end if    
      ! read number of paths:
      if (index(line,'npath:').gt.0) then
        begin_to_read=index(line,'npath:')+6
        read(line(begin_to_read:),*) npath
        print '(8x,I10," paths")', npath
        allocate(segment_nqpoint(npath),labels(npath,2))
        ! points per path
        !print*,nqpoint,npath
        print '(8x,"On average ",I0," qpoints/path")',int(nqpoint/npath)
        segment_nqpoint(:)=int(nqpoint/npath)
        print '(8x,"Initial path labels:")'
        do ipath=1,npath
          chardum=''
          write(chardum,'(I0,"i")') ipath
          labels(ipath,1)=chardum
          write(chardum,'(I0,"f")') ipath
          labels(ipath,2)=chardum
          print '(8x,A24,", ",A24)', labels(ipath,1), labels(ipath,2)
        end do
      end if    
      ! if present, read number of qpoints per path:
      if (index(line,'segment_nqpoint:').gt.0) then
        print '(8x,"segment_nqpoints:")'
        do ipath=1,npath
          read(51,*) chardum, segment_nqpoint(ipath)
          print '(8x,I10)', segment_nqpoint(ipath)
        end do
      end if    
      ! if present, read labels of paths:
      if (index(line,'labels:').gt.0) then
        print '(8x,"path labels:")'
        do ipath=1,npath
          read(51,*) chardum, chardum, labels(ipath,1), labels(ipath,2)
          print '(8x,A24,", ",A24)', labels(ipath,1), labels(ipath,2)
        end do
      end if    
      ! read number of atoms:
      if (index(line,'natom:').gt.0) then
        begin_to_read=index(line,'natom:')+6
        read(line(begin_to_read:),*) natom
        print '(8x,I10," atoms")',natom
        allocate(freq(nqpoint,3*natom),eigenvec(nqpoint,3*natom,natom,    &
     &           1:3))
      end if    
      !do ipath=1,npath
      !  print '(8x,I10)', segment_nqpoint(ipath)
      !end do
      !print '(8x,"path labels:")'
      !do ipath=1,npath
      !  print '(8x,A24,", ",A24)', labels(ipath,1), labels(ipath,2)
      !end do
      ! read phonon frequencies:
      if (index(line,'phonon:').gt.0) then
        do iqpoint=1,nqpoint
          !print '(8x,"iqpoint=",I0)',iqpoint
          read(51,*) chardum,chardum,chardum,qpos(iqpoint,1:3)
          !print '(8x,"qpos=",3(F12.6))',qpos(iqpoint,1:3)
          read(51,*) chardum, distance(iqpoint)
          !print '(8x,"distance=",F12.6)',distance(iqpoint)
          read(51,*) line ! "band"
          do iband=1,3*natom
            !print '(8x,"iband=",I0)',iband
            read(51,*) line ! band number
            read(51,*) chardum, freq(iqpoint,iband)
            !print '(8x,"frequency=",F12.6)',freq(iqpoint,iband)
            if(has_eigenvec) then
              read(51,*) line ! "eigenvector"
              !print*,line
              do iatom=1,natom
                !print '(8x,"iatom=",I0)',iatom
                read(51,*) line ! "atom #"
                !print*,line
                read(51,*) chardum,chardum,re_ev,im_ev
                !print*,re_ev,im_ev
                eigenvec(iqpoint,iband,iatom,1)=complex(re_ev,im_ev)
                read(51,*) chardum,chardum,re_ev,im_ev
                !print*,re_ev,im_ev
                eigenvec(iqpoint,iband,iatom,2)=complex(re_ev,im_ev)
                read(51,*) chardum,chardum,re_ev,im_ev
                !print*,re_ev,im_ev
                eigenvec(iqpoint,iband,iatom,3)=complex(re_ev,im_ev)
                !print '(6(F12.6))',eigenvec(iqpoint,iband,iatom,:)
              end do ! iatom
            end if ! has_eigenvec
          end do ! iband
          read(51,'(A64)') line
          !print'(A64)', line
        end do ! iqpoint
      end if    
      goto 10
      !
20    continue      
      close(51)
      !
      ! write output file:
      open(52,file=outfile,status='replace')
      if (has_eigenvec) then
      write(52,'("# distance  frequency (bands are separated by blank li  &
     &nes)  atomic weight in EV atom1...natom  full eigenvecs (atom 1 Re  &
     &(v_x) Im(v_x) Re(v_y) Im(v_y) Re(v_z) Im(v_z)  atom 2 Re(v_x)...)   &
     &")')
      else
      write(52,'("# distance  frequency (bands are separated by blank li  &
     &nes)")')
      end if
      write(52,'("# End points of segments: ")')
      write(52,'("#",F13.8)',advance="no") distance(1)
      iqpoint=0
      do ipath=1,npath
        iqpoint=iqpoint+segment_nqpoint(ipath)
        write(52,'(F11.8)',advance="no") distance(iqpoint)
      end do
      write(52,'(" ")',advance="yes")
      do iband=1,3*natom
        iqpoint=1
        do ipath=1,npath
          do segment_iqpoint=1,segment_nqpoint(ipath)
            if (has_eigenvec) then
              write(52,'(F9.6,x,F12.6)',advance="no") distance(iqpoint),  &
     &          freq(iqpoint,iband)
              do iatom=1,natom
                write(52,'(6(F12.6))',advance="no")                       &
     &           norm2(cabs(eigenvec(iqpoint,iband,iatom,:)))
              end do
              do iatom=1,natom
                write(52,'(6(F12.6))',advance="no")                       &
     &           eigenvec(iqpoint,iband,iatom,:)
              end do
              write(52,'(" ")',advance="yes")
            else
              write(52,'(F9.6,x,F12.6)') distance(iqpoint),               &
     &          freq(iqpoint,iband)
            end if
            iqpoint=iqpoint+1
          end do ! segment_iqpoint
          write(52,'("")')
        end do ! ipath
        write(52,'("")')
      end do ! iband  
      write(52,'("")')
      close(52)
      !
      !
      ! end normally :
      !
      print fsubendext,'phonopy_print_bands'  
      !
      return
      !
100   continue
      close(51)
      call error_stop('problem with reading input file')     
      !
      end subroutine phonopy_print_bands

!------------------------------------------------------

      subroutine phonopy_print_eigenvecs(infile,outfile,qpoint,band)
      !  
      use defs, only: fsubstart,fsubendext, error_stop,atom,element
      use misc, only: getspecies,frac2abs
      use writecoords, only: write_coords,write_xsf_forces
      implicit none
      !
      character(len=*) :: infile, outfile
      integer nqpoint,npath,natom,nspecies,qpoint,band
      integer, allocatable :: segment_nqpoint(:)
      character*64, allocatable :: labels(:,:)
      double precision, allocatable :: qpos(:,:),freq(:,:),distance(:)
      complex, allocatable :: eigenvec(:,:,:,:)
      double precision vecs(3,3)
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      !
      character*128 line,chardum,outfile_re,outfile_im
      integer begin_to_read, ipath,iband,iatom,iqpoint,segment_iqpoint,i
      logical has_eigenvec
      double precision re_ev,im_ev
      !
      print fsubstart,'phonopy_print_eigenvecs'  
      !
      ! read input file (band.yaml):
      !
      print '(8x,"reading phonon eigenvecs from file ",A64)',             &
     &     adjustl(trim(infile))
      open(51,file=infile,status="old")
      rewind(51)
      !
      ! check if eigenvectors are there:
      has_eigenvec=.false.
8     read(51,'(A)', end=100, err=100) line 
      if (index(line,'eigenvector').gt.0) then
        has_eigenvec=.true.
        print '(8x,"Eigenvecs found")'
        goto 9
      end if
      goto 8

9     rewind(51)      
      !
10    read(51,'(A128 )', end=20, err=100) line 
      ! read number of qpoints:
      if (index(line,'nqpoint:').gt.0 .and.index(line,'segment').le.0)    &
     &    then
        begin_to_read=index(line,'nqpoint:')+8
        read(line(begin_to_read:),*) nqpoint
        print '(8x,I10," qpoints")', nqpoint
        allocate(qpos(nqpoint,1:3),distance(nqpoint))
      end if    
      ! read number of paths:
      if (index(line,'npath:').gt.0) then
        begin_to_read=index(line,'npath:')+6
        read(line(begin_to_read:),*) npath
        print '(8x,I10," paths")', npath
        allocate(segment_nqpoint(npath),labels(npath,2))
      end if    
      ! read number of qpoints per path:
      if (index(line,'segment_nqpoint:').gt.0) then
        print '(8x,"segment_nqpoints:")'
        do ipath=1,npath
          read(51,*) chardum, segment_nqpoint(ipath)
          print '(8x,I10)', segment_nqpoint(ipath)
        end do
      end if    
      ! read labels of paths:
      if (index(line,'labels:').gt.0) then
        print '(8x,"path labels:")'
        do ipath=1,npath
          read(51,*) chardum, chardum, labels(ipath,1), labels(ipath,2)
          print '(8x,A24,", ",A24)', labels(ipath,1), labels(ipath,2)
        end do
      end if    
      ! read number of atoms:
      if (index(line,'natom:').gt.0) then
        begin_to_read=index(line,'natom:')+6
        read(line(begin_to_read:),*) natom
        print '(8x,I10," atoms")',natom
        allocate(freq(nqpoint,3*natom),eigenvec(nqpoint,3*natom,natom,    &
     &           1:3))
        if(allocated(atoms)) deallocate(atoms)
        allocate(atoms(natom))
      end if    
      ! read lattice vecs:
      if (index(line,'lattice:').gt.0) then
        print '(8x,"lattice vecs:")'
        do i=1,3
          read(51,*) chardum,chardum,vecs(i,1:3)
          print '(8x,3(F12.6))',vecs(i,1:3)
        end do
      end if    
      ! read atom positions:
      if (index(line,'points:').gt.0) then
        print '(8x,"atom positions:")'
        do iatom=1,natom
          read(51,*) chardum,chardum,atoms(iatom)%name
          read(51,*) chardum,chardum,atoms(iatom)%where(1:3)
          read(51,*) chardum,atoms(iatom)%mass
          print '(8x,A2,x,3(F12.6))',atoms(iatom)%name(1:2),              &
     &       atoms(iatom)%where(1:3)
        end do
      end if    
      ! read phonon frequencies:
      if (index(line,'phonon:').gt.0) then
        do iqpoint=1,nqpoint
          read(51,*) chardum,chardum,chardum,qpos(iqpoint,1:3)
          read(51,*) chardum, distance(iqpoint)
          read(51,*) line ! "band"
          do iband=1,3*natom
            read(51,*) line ! band number
            read(51,*) chardum, freq(iqpoint,iband)
            if(has_eigenvec) then
              read(51,*) line ! "eigenvector"
              do iatom=1,natom
                read(51,*) line ! "atom #"
                read(51,*) chardum,chardum,re_ev,im_ev
                eigenvec(iqpoint,iband,iatom,1)=complex(re_ev,im_ev)
                read(51,*) chardum,chardum,re_ev,im_ev
                eigenvec(iqpoint,iband,iatom,2)=complex(re_ev,im_ev)
                read(51,*) chardum,chardum,re_ev,im_ev
                eigenvec(iqpoint,iband,iatom,3)=complex(re_ev,im_ev)
              end do ! iatom
            end if ! has_eigenvec
          end do ! iband
          read(51,'(A64)') line
        end do ! iqpoint
      end if    
      goto 10
      !
20    continue      
      close(51)
      !
      print '(/,8x,"printing eigenvector for qpoint ",I0,", band ",I0,/)  &
     &',qpoint,band
      outfile_re=trim(outfile)
      outfile_im=trim(outfile)
      write(outfile_re(len_trim(outfile)+1:),'("_RE.xsf")')
      write(outfile_im(len_trim(outfile)+1:),'("_IM.xsf")')
      do iatom=1,natom
        call frac2abs(atoms(iatom)%where,vecs,atoms(iatom)%abswhere)
        atoms(iatom)%force(:)=real(eigenvec(qpoint,band,iatom,:))
      end do ! iatom
      call getspecies(atoms,species)
      call write_coords(outfile_re,'xsf_forces',atoms,natom,species,      &
     &   nspecies,vecs)
      do iatom=1,natom
        atoms(iatom)%force(:)=aimag(eigenvec(qpoint,band,iatom,:))
      end do ! iatom
      call write_coords(outfile_im,'xsf_forces',atoms,natom,species,      &
     &   nspecies,vecs)
!      ! write output file:
!      open(52,file=outfile,status='replace')
!      if (has_eigenvec) then
!      write(52,'("# distance  frequency (bands are separated by blank li  &
!     &nes)  atomic weight in EV atom1...natom  full eigenvecs (atom 1 Re  &
!     &(v_x) Im(v_x) Re(v_y) Im(v_y) Re(v_z) Im(v_z)  atom 2 Re(v_x)...)   &
!     &")')
!      else
!      write(52,'("# distance  frequency (bands are separated by blank li  &
!     &nes)")')
!      end if
!      write(52,'("# End points of segments: ")')
!      write(52,'("#",F13.8)',advance="no") distance(1)
!      iqpoint=0
!      do ipath=1,npath
!        iqpoint=iqpoint+segment_nqpoint(ipath)
!        write(52,'(F11.8)',advance="no") distance(iqpoint)
!      end do
!      write(52,'(" ")',advance="yes")
!      do iband=1,3*natom
!        iqpoint=1
!        do ipath=1,npath
!          do segment_iqpoint=1,segment_nqpoint(ipath)
!            if (has_eigenvec) then
!              write(52,'(F9.6,x,F12.6)',advance="no") distance(iqpoint),  &
!     &          freq(iqpoint,iband)
!              do iatom=1,natom
!                write(52,'(6(F12.6))',advance="no")                       &
!     &           norm2(cabs(eigenvec(iqpoint,iband,iatom,:)))
!              end do
!              do iatom=1,natom
!                write(52,'(6(F12.6))',advance="no")                       &
!     &           eigenvec(iqpoint,iband,iatom,:)
!              end do
!              write(52,'(" ")',advance="yes")
!            else
!              write(52,'(F9.6,x,F12.6)') distance(iqpoint),               &
!     &          freq(iqpoint,iband)
!            end if
!            iqpoint=iqpoint+1
!          end do ! segment_iqpoint
!          write(52,'("")')
!        end do ! ipath
!        write(52,'("")')
!      end do ! iband  
!      write(52,'("")')
!      close(52)
      !
      !
      ! end normally :
      !
      print fsubendext,'phonopy_print_eigenvecs'  
      !
      return
      !
100   continue
      close(51)
      call error_stop('no eigenvectors found')     
      !
      end subroutine phonopy_print_eigenvecs

c---------------------------------------------------------------------

      subroutine vasp_get_pol()
      use defs
      !use linalg
      use misc, only : vecs2vol
      use readcoords, only : read_outcar
      implicit none 
      integer i,j
      character*1024 :: line,chardum
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      double precision :: pol_tot(1:3),pol_ele(1:3),pol_ion(1:3),         &
     &                    q_of_pol(1:3,1:3)
      double precision dip_ele(1:3),dip_ion(1:3),dip_tot(1:3)
      double precision vecs(1:3,1:3),vol,polfac1,polfac2
      integer nspecies, natoms
      !double precision, allocatable :: epseigs(:),epsev(:,:)
      !double precision, allocatable :: epseigsim(:)
      !double precision, optional :: rotmat(1:3,1:3)
      !double precision epsrot(3,3)
      !
      print fsubstart,'get_vasp_pol'  
      !
      ! begin read vecs, volume 
      !
      call read_outcar('OUTCAR',atoms,natoms,species,nspecies,
     &           vecs)
      call vecs2vol(vecs,vol)
      print '(8x,"cell volume = ", F15.6, " Angs³")', vol
      polfac1=16.02176620000000D0! ec/Angs² --> C/m²
      polfac2=1602.176620000000D0! ec/Angs² --> uC/cm²
      !
      ! end read vecs, volume 
      !
      !
      open(51,file='OUTCAR',status='old', err=100)  
      rewind(51)
      !
10    read(51,'(A1024)',end=12,err=101) line
      if (index(line,'Ionic dipole moment:').gt.0) then
        read(line,*) chardum,chardum,chardum,chardum,dip_ion(1:3)
        dip_ion=-dip_ion
        if (talk) print '(8x,"ionic dipole moment=",1x,3(F12.6),1x,       &
     &                  "|e| Angs")',dip_ion
        pol_ion=dip_ion/vol
        print '(8x,"ionic polarization= ",3(F15.6)," |e| / Angs²")',      &
     &  pol_ion
        print '(8x,"ionic polarization=",1x,3(F15.6),1x," C/m² ")',       &
     &  pol_ion*polfac1     
        print '(8x,"ionic polarization=",1x,3(F15.6),1x," uC/cm² ")',     &
     &  pol_ion*polfac2     
      end if
      if (index(line,'Total electronic dipole moment:').gt.0) then
        read(line,*)chardum,chardum,chardum,chardum,chardum,dip_ele(1:3)
        dip_ele=-dip_ele
        print '(8x,"electronic dipole moment=",1x,3(F12.6),1x,            &
     &                  "|e| Angs")',dip_ele
        pol_ele=dip_ele/vol
        print '(8x,"electronic polarization= ",3(F15.6),                  &
     &   " |e| / Angs²")',  pol_ele
        print '(8x,"electronic polarization= ",3(F15.6)," C/m² ")',       &
     &  pol_ele*polfac1   
        print '(8x,"electronic polarization= ",3(F15.6)," uC/cm² ")',     &
     &  pol_ele*polfac2     
      end if
      goto 10

12    continue
      close(51) 
      dip_tot=dip_ion+dip_ele
        print '(8x,"Total dipole moment=",1x,3(F12.6),1x,"|e| Angs")',    &
     &dip_tot
      pol_tot=pol_ion+pol_ele
      print '(8x,"total polarization= ",3(F15.6)," |e| / Angs²")',        &
     &pol_tot
      print '(8x,"total polarization= ",3(F15.6)," C/m² ")',              &
     &pol_tot*polfac1   
      print '(8x,"total polarization= ",3(F15.6)," uC/cm² ")',            &
     &pol_tot*polfac2     
      q_of_pol(1,:)=ec*vecs(1,:)/vol*1.0E20
      q_of_pol(2,:)=ec*vecs(2,:)/vol*1.0E20
      q_of_pol(3,:)=ec*vecs(3,:)/vol*1.0E20
      print '(8x,"Quantum of polarization 1 = ",3(F15.6)," C/m² ")',      &
     &q_of_pol(1,:)
      print '(8x,"Quantum of polarization 2 = ",3(F15.6)," C/m² ")',      &
     &q_of_pol(2,:)
      print '(8x,"Quantum of polarization 3 = ",3(F15.6)," C/m² ")',      &
     &q_of_pol(3,:)
      print '(8x,"Quantum of polarization 1 = ",3(F15.6)," uC/cm² ")',    &
     &q_of_pol(1,:)*100.0D0
      print '(8x,"Quantum of polarization 2 = ",3(F15.6)," uC/cm² ")',    &
     &q_of_pol(2,:)*100.0D0
      print '(8x,"Quantum of polarization 3 = ",3(F15.6)," uC/cm² ")',    &
     &q_of_pol(3,:)*100.0D0

      !
      print fsubendext,'vasp_get_pol'  
      return
      !
100   continue    
      close(51)
      nerr=nerr+1
      print ferrmssg,'OUTCAR not found'
      return
101   continue    
      close(51)
      nerr=nerr+1
      print ferrmssg,'EPS(omega) not found'
      return
      !
      end subroutine vasp_get_pol

!------------------------------------------------------------------

      subroutine vasp_read_fermi_energy(efermi)
      !
      use defs, only : error_stop,fsubstart,fsubendext
      implicit none

      integer fileunit
      logical file_is_open
      double precision efermi
      character(len=256) line
      integer iread
      !
      print fsubstart,"vasp_read_fermi_energy"
      !
      ! begin read Fermi energy
      !
      efermi=1000000.0d0
      fileunit=51
      INQUIRE (unit=fileunit, opened=file_is_open)
      do while (file_is_open)
        fileunit=fileunit+1
        INQUIRE (unit=fileunit, opened=file_is_open)
      end do
      print '(8x,"file unit=",I0)',fileunit
      open(fileunit,file="OUTCAR",status='old', err=21)
      rewind(fileunit)
10    read(fileunit,'(A256)',err=21,end=20) line
      if(index(line,'E-fermi').gt.0) then
        iread=index(line,'E-fermi')+9
        read(line(iread:),*) efermi
      end if
      goto 10
      !
20    continue
      close(fileunit)
      print '(8x,"E_Fermi=",F12.6," eV")',efermi
      print fsubendext,"vasp_read_fermi_energy"
      return

21    continue
      close(fileunit)
      call error_stop("problem with OUTCAR") 
      
      end subroutine vasp_read_fermi_energy 

!-------------------------------------------------------------------      

      subroutine vasp_read_eigenvalues_and_weights(eigenvalues,kweights)

      use defs, only : error_stop, fsubstart,fsubendext
      implicit none

      integer fileunit
      logical file_is_open
      integer idum,iline,ikpoint,iband
      integer nele,nkpoints,nbands,nspins
      character(len=1024) line
      double precision, allocatable :: eigenvalues(:,:,:),kweights(:)
      real rdum

      print fsubstart,"vasp_read_eigenvalues"
      fileunit=51
      INQUIRE (unit=fileunit, opened=file_is_open)
      do while (file_is_open)
        fileunit=fileunit+1
        INQUIRE (unit=fileunit, opened=file_is_open)
      end do
      print '(8x,"file unit=",I0)',fileunit
      open(fileunit,file="EIGENVAL", status='old')
      !
      ! begin read nbands, nkpoints
      !    
      read(fileunit,*,end=21,err=21) idum,idum,idum, nspins
      do iline=1,4
        read(fileunit,'(A1024)',end=21,err=21) line    
      end do
      read(fileunit,*,end=21,err=21) nele,nkpoints,nbands
      print '(8x,"nele, nkpts, nbands=",3(1x,I0))',nele,nkpoints,nbands
      allocate(kweights(1:nkpoints),eigenvalues(nkpoints,nbands,nspins))
      !
      ! end read nbands, nkpoints
      !
      !
      ! begin read kpoints and eigenvalues
      !    
      do ikpoint=1,nkpoints
        read(fileunit,'(A1024)',end=21,err=21) line    
        read(fileunit,*,end=21,err=21) rdum,rdum,rdum,kweights(ikpoint)
        !print '(8x,"kweight=",F12.6)',kweights(ikpoint)
        do iband=1,nbands
            read(fileunit,*,err=21,end=21) idum,                          &
     &         eigenvalues(ikpoint,iband,1:nspins)
!            print '(8x,"iband, eigenvalue=",I0,x,F12.6)', iband,          &
!     &          eigenvalues(ikpoint,iband,1:nspins)
        end do
      end do
      close(fileunit)
      !
!      open(fileunit, file="OUTCAR",status="old")
!      rewind(fileunit)
!      line=" "
!      do while (index(line,'Following reciprocal coordinates').le.0)
!        read(fileunit,'(A256)') line
!      end do
!      read (fileunit,'(A256)') line
!      !print*,"line before=",line ! DEBUG
!      do ikpoint=1,nkpoints
!        read(fileunit,*) rdum,rdum,rdum,kweights(ikpoint)
!        print*,kweights(ikpoint)
!      end do ! ikpoint
!      close(fileunit)
      !
      ! end read kpoints and eigenvalues
      !
      !
20    continue      
      close(fileunit)
      print fsubendext,"vasp_read_eigenvalues"
      return
      !
21    continue
      close(fileunit)
      call error_stop("problem with EIGENVAL")

      end subroutine vasp_read_eigenvalues_and_weights

!------------------------------------------------

      subroutine vasp_read_volume(vol)
      use defs, only : error_stop,fsubstart,fsubendext
      use misc, only : vecs2vol
      implicit none

      integer fileunit
      logical file_is_open
      double precision vol,vecs(1:3,1:3)
      character(len=256) line
      integer iread
      !
      print fsubendext,"vasp_read_volume"
      !
      ! begin read volume
      !
      vol=-1.0d0
      fileunit=51
      INQUIRE (unit=fileunit, opened=file_is_open)
      do while (file_is_open)
        fileunit=fileunit+1
        INQUIRE (unit=fileunit, opened=file_is_open)
      end do
      print '(8x,"file unit=",I0)',fileunit
      open(fileunit,file="OUTCAR",status='old', err=21)
      rewind(fileunit)
10    read(fileunit,'(A256)',err=21,end=20) line
      if(index(line,'direct lattice vectors').gt.0) then
        read(fileunit,*) vecs(1,1:3)
        read(fileunit,*) vecs(2,1:3)
        read(fileunit,*) vecs(3,1:3)
        call vecs2vol(vecs,vol)
      end if
      goto 10
      !
20    continue
      close(fileunit)
      print '(8x,"Volume=",F15.6," Angstrom³")',vol
      print fsubendext,"vasp_read_volume"
      return

21    continue
      close(fileunit)
      call error_stop("problem with OUTCAR") 
      
      end subroutine vasp_read_volume 

!------------------------------------------------

      end module 
