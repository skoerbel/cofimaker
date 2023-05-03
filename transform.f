      module transform

      implicit none

      contains

c---------------------------------------------------------------------

      subroutine transformvectors(vecs1,vecs2,atoms1,newatoms,nsc)
     
      use defs
      !use readcoords
      !use writecoords
      use misc
      implicit none

      ! variables
      type(atom), intent(in) :: atoms1(:)
      type(atom), allocatable :: newatoms(:)
      integer :: natoms,nnewatoms,nsc(1:3)
      double precision, intent(in) :: vecs1(3,3),vecs2(3,3)
 
      ! local variables
      double precision :: matrix(1:3,1:3),imatrix(1:3,1:3)
      double precision :: inmat(1:3,1:3),outmat(3,3),outimat(1:3,1:3)
      integer iatom,i,j,k,n,p,q,r,s,nuc(1:3),iallatom,nallatoms
      double precision xold(1:3),xnew(1:3)
      double precision newvecs(1:3,1:3)
      double precision invinmat(1:3,1:3),inewvecs(1:3,1:3)
      double precision :: tmatrix(1:3,1:3)
      double precision alat(1:3),angles(1:3)
      double precision, parameter :: tol=1.0D-6
      logical,allocatable :: delete(:)
      type(atom), allocatable :: allatoms(:)
      double precision tolat(3)

      talk=.true.
      if(talk)print fsubstart,"transformvectors"
      print '(8x,"old vecs: ")'
      print '(8x,3(F10.6,1x))', vecs1(1,1:3)
      print '(8x,3(F10.6,1x))', vecs1(2,1:3)
      print '(8x,3(F10.6,1x))', vecs1(3,1:3)
      print '(8x,"new vecs: ")'
      print '(8x,3(F10.6,1x))', vecs2(1,1:3)
      print '(8x,3(F10.6,1x))', vecs2(2,1:3)
      print '(8x,3(F10.6,1x))', vecs2(3,1:3)
      
      natoms=size(atoms1)
      nuc(1:3)=1
      tolat=1.0E-3 ! tolerance in fractional coordinates to distinguish between atomic positions
      !print '(8x,"old cell: ",I0, " atoms")',natoms
      do j=1,3
        inmat(1:3,j)=vecs1(j,1:3)
        outmat(1:3,j)=vecs2(j,1:3)
        !inmat(1:3,j)=vecs1(1:3,j)
        !outmat(1:3,j)=vecs2(1:3,j)
      end do
      call invert_matrix(outmat,outimat)
      matrix=matmul(outimat,inmat)
      ! take into account desired supercell dimension
      do n=1,3
        matrix(n,1:3)=matrix(n,1:3)/dble(nsc(n))
      end do  
      nallatoms=natoms*nuc(1)*nuc(2)*nuc(3)*8
      allocate(allatoms(1:nallatoms)) ! natoms*nuc(1)*nuc(2)*nuc(3)))
      allocate(delete(1:nallatoms)) !natoms*nuc(1)*nuc(2)*nuc(3)))
      
      ! multiply with matrix
      !print '(8x,"transformation matrix: ")'
      !print '(8x,3(F10.6,1x))', matrix(1,1:3)
      !print '(8x,3(F10.6,1x))', matrix(2,1:3)
      !print '(8x,3(F10.6,1x))', matrix(3,1:3)
      iallatom=0
      !print '(8x,"nuc=",3(I1,x))',nuc
      do iatom=1,natoms
            !print '(8x,"iatom=",I0)', iatom
            xold=atoms1(iatom)%where
            !print '(8x,"atom(iatom)=",A4)', atoms1(iatom)%name
            xnew=0.0D0
            do q=-nuc(1),nuc(1)-1  !0
              do r=-nuc(2),nuc(2)-1 !0
                do s=-nuc(3),nuc(3)-1 !0
                  !print '(8x,"r,s,t=",3(I0,1x))', q,r,s
                  iallatom=iallatom+1
                  xnew(1)=(xold(1)+float(q))*matrix(1,1)
     &                      +(xold(2)+float(r))*matrix(1,2)
     &                      +(xold(3)+float(s))*matrix(1,3)
                  xnew(2)=(xold(1)+float(q))*matrix(2,1)
     &                      +(xold(2)+float(r))*matrix(2,2)
     &                      +(xold(3)+float(s))*matrix(2,3)
                  xnew(3)=(xold(1)+float(q))*matrix(3,1)
     &                      +(xold(2)+float(r))*matrix(3,2)
     &                      +(xold(3)+float(s))*matrix(3,3)
                  ! shift back to between 0 and 1
                  !print '(8x,"going to shift back...")'
                  do j=1,3
                    do while(xnew(j).ge.1.0D0-tolfrac(j))
                      xnew(j)=xnew(j)-1.0D0
                    end do
                    do while(xnew(j).lt.-tolfrac(j))
                      xnew(j)=xnew(j)+1.0D0
                    end do
                  end do
!                 print '(8x,"going to copy atom ",A4, " to ",I1)', 
!     &               atoms1(iatom)%name,iallatom
                  allatoms(iallatom)=atoms1(iatom)
                  !print '(8x,"going to update atom position...")'
                  allatoms(iallatom)%where=xnew
                  !print '(8x,"end of (rst) tuple")'
                end do
              end do
            end do
      end do
      !print '(8x,"size(allatoms): ",I0)', size(allatoms)
      
      ! delete doubled coordinates
      delete=.False.
      do n=1,nallatoms-1 !natoms*nuc(1)*nuc(2)*nuc(3)-1
            do p=n+1,nallatoms !natoms*nuc(1)*nuc(2)*nuc(3)
              if(abs(allatoms(p)%where(1)-allatoms(n)%where(1))
     &          .le.tol.and.abs(allatoms(p)%where(2)
     &          -allatoms(n)%where(2)).le.tol.and.
     &               abs(allatoms(p)%where(3)-allatoms(n)%where(3))
     &               .le.tol.and.
     &                  allatoms(p)%core.eq.allatoms(n)%core) then
                        delete(p)=.True.
              end if
            end do
      end do

      ! find final number of coordinates
      nnewatoms=0
      do n=1,nallatoms !natoms*nuc(1)*nuc(2)*nuc(3)
            if(delete(n).eqv..False.) nnewatoms=nnewatoms+1
      end do
      ! get new atoms
      allocate(newatoms(nnewatoms))
      !print '(8x,"found ",I0, " new atoms")',nnewatoms
      iatom=0
      do n=1,nallatoms !natoms*nuc(1)*nuc(2)*nuc(3)
            if (delete(n).eqv..False.) then
                  iatom=iatom+1
                  newatoms(iatom)=allatoms(n)
            end if
      end do
      !
      call rm_double_atoms(newatoms,tolat)
      !
      !
      if(talk)print fsubendext, 'transformvectors'
      !
      end subroutine transformvectors

c---------------------------------------------------------------------
!
!      subroutine transformmatrix(matrix,fincoords,foutcoords,finform,
!     &                           nsc)
!     
!      use defs
!      use readcoords
!      use writecoords
!      use misc
!      implicit none
!
!      ! variables
!      double precision :: matrix(1:3,1:3),imatrix(1:3,1:3)
!      character(len=*) :: fincoords,foutcoords,finform
!      type(atom), allocatable :: atoms(:),newatoms(:)
!      type(element), allocatable :: species(:)
!      integer :: natoms,nspecies,nnewatoms,nsc(1:3)
! 
!      ! local variables
!      integer iatom,i,j,k,n,p,q,r,s,nuc(1:3),iallatom,nallatoms
!      integer ispecies
!      double precision xold(1:3),xnew(1:3),vecs(1:3,1:3)
!      double precision newvecs(1:3,1:3),inmat(1:3,1:3)
!      double precision invinmat(1:3,1:3),inewvecs(1:3,1:3)
!      double precision :: tmatrix(1:3,1:3)
!      double precision alat(1:3),angles(1:3)
!      logical,allocatable :: delete(:)
!      type(atom), allocatable :: allatoms(:)
!
!      talk=.true.
!      if(talk)print fsubstart,"transformmatrix"
!      
!      ! read in  coordinates ...
!      if(talk) then
!            print '(" ")'
!            print '(8x,"calling read_coords...")'
!      end if
!      call read_coords(fincoords,finform,atoms,natoms,species,
!     &           nspecies,vecs)
!      if(talk) print '(8x,"Coords read in.")'
!      
!      !! take also atoms of neighboring old unit cells into account, if
!      !! the new cell lies partly outside the old one
!      !nuc(1)=3
!      !nuc(2)=3
!      !nuc(3)=3
!      nuc(1)=1
!      nuc(2)=1
!      nuc(3)=1
!      
!      ! take into account desired supercell dimension
!      do n=1,3
!        matrix(n,1:3)=matrix(n,1:3)/dble(nsc(n))
!        !vecs(n,1:3)=vecs(n,1:3)*dble(nsc(n))
!      end do  
!      nallatoms=natoms*nuc(1)*nuc(2)*nuc(3)*8
!      allocate(allatoms(1:nallatoms)) ! natoms*nuc(1)*nuc(2)*nuc(3)))
!      !allocate(x_new(1:natom*nuc(1)*nuc(2)*nuc(3),1:3))
!      allocate(delete(1:nallatoms)) !natoms*nuc(1)*nuc(2)*nuc(3)))
!      
!      ! multiply with matrix
!      iallatom=0
!      do iatom=1,natoms
!            xold=atoms(iatom)%where
!            xnew=0.0D0
!            do q=-nuc(1),nuc(1)-1  !0
!              do r=-nuc(2),nuc(2)-1 !0
!                do s=-nuc(3),nuc(3)-1 !0
!                  iallatom=iallatom+1
!                  xnew(1)=(xold(1)+float(q))*matrix(1,1)
!     &                      +(xold(2)+float(r))*matrix(1,2)
!     &                      +(xold(3)+float(s))*matrix(1,3)
!                  xnew(2)=(xold(1)+float(q))*matrix(2,1)
!     &                      +(xold(2)+float(r))*matrix(2,2)
!     &                      +(xold(3)+float(s))*matrix(2,3)
!                  xnew(3)=(xold(1)+float(q))*matrix(3,1)
!     &                      +(xold(2)+float(r))*matrix(3,2)
!     &                      +(xold(3)+float(s))*matrix(3,3)
!                  ! shift back to between 0 and 1
!                  do j=1,3
!                    do while(xnew(j).ge.1.0D0-tolfrac(j))
!                      xnew(j)=xnew(j)-1.0D0
!                    end do
!                    do while(xnew(j).lt.-tolfrac(j))
!                      xnew(j)=xnew(j)+1.0D0
!                    end do
!                  end do
!                  !do j=1,3
!                  !  do k=1,3
!                  !    xnew(j)=xnew(j)+matrix(j,k)*xold(k)
!                  !  end do
!                  !end do
!                  allatoms(iallatom)=atoms(iatom)
!                  allatoms(iallatom)%where=xnew
!                end do
!              end do
!            end do
!      end do
!      
!      ! delete doubled coordinates
!      delete=.False.
!      do n=1,nallatoms-1 !natoms*nuc(1)*nuc(2)*nuc(3)-1
!            do p=n+1,nallatoms !natoms*nuc(1)*nuc(2)*nuc(3)
!              if(abs(allatoms(p)%where(1)-allatoms(n)%where(1))
!     &          .le.1E-6.and.abs(allatoms(p)%where(2)
!     &          -allatoms(n)%where(2)).le.1E-6.and.
!     &               abs(allatoms(p)%where(3)-allatoms(n)%where(3))
!     &               .le.1E-6.and.
!     &                  allatoms(p)%core.eq.allatoms(n)%core) then
!
!                        delete(p)=.True.
!              end if
!            end do
!      end do
!!      ! delete doubled coordinates
!!      delete=.False.
!!      do n=1,nallatoms-1 !natoms*nuc(1)*nuc(2)*nuc(3)-1
!!            do p=n+1,nallatoms !natoms*nuc(1)*nuc(2)*nuc(3)
!!              if(abs(allatoms(p)%where(1)-allatoms(n)%where(1))
!!     &          .le.tolfrac(1).and.abs(allatoms(p)%where(2)
!!     &          -allatoms(n)%where(2)).le.tolfrac(2).and.
!!     &               abs(allatoms(p)%where(3)-allatoms(n)%where(3))
!!     &               .le.tolfrac(3)) then
!!                        delete(p)=.True.
!!              end if
!!            end do
!!      end do
!
!
!      ! find final number of coordinates
!      nnewatoms=0
!      do n=1,nallatoms !natoms*nuc(1)*nuc(2)*nuc(3)
!            if(delete(n).eqv..False.) nnewatoms=nnewatoms+1
!      end do
!      ! get new atoms
!      allocate(newatoms(nnewatoms))
!      iatom=0
!      do n=1,nallatoms !natoms*nuc(1)*nuc(2)*nuc(3)
!            if (delete(n).eqv..False.) then
!                  iatom=iatom+1
!                  newatoms(iatom)=allatoms(n)
!            end if
!      end do
!      
!      call getspecies(newatoms,species)
!      nspecies=size(species)
!
!      ! transform vectors
!      call transpon_matrix(vecs,inmat)
!!      print*, "Vi:"
!!      do j=1,3
!!        print*, inmat(j,1:3)
!!      end do
!      call invert_matrix(inmat,invinmat)
!!      print*, "(Vi)^-1:"
!!      do j=1,3
!!        print*, invinmat(j,1:3)
!!      end do
!      !call transpon_matrix(matrix,tmatrix)
!      call mymatmul(matrix,invinmat,inewvecs)
!!      print*, "M*(Vi)^-1:"
!!      do j=1,3
!!        print*, inewvecs(j,1:3)
!!      end do
!      call invert_matrix(inewvecs,newvecs)
!!      print*, "(M*(Vi)^-1)^-1:"
!!      do j=1,3
!!        print*, newvecs(j,1:3)
!!      end do
!      call transpon_matrix(newvecs,inewvecs)
!!      print*, "((M*(Vi)^-1)^-1)^t:"
!!      do j=1,3
!!        print*, inewvecs(j,1:3)
!!      end do
!      vecs=inewvecs
!      ! rotate vectors so that they are conform with lammps
!      do i=1,3
!        alat(i)=sqrt(vecs(i,1)**2+vecs(i,2)**2+vecs(i,3)**2)
!      end do
!      angles(1)=acos((vecs(2,1)*vecs(3,1)+vecs(2,2)*vecs(3,2)
!     &          +vecs(2,3)*vecs(3,3))/(alat(2)*alat(3)))*180.0D0/Pi
!      angles(2)=acos((vecs(1,1)*vecs(3,1)+vecs(1,2)*vecs(3,2)
!     &          +vecs(1,3)*vecs(3,3))/(alat(1)*alat(3)))*180.0D0/Pi
!      angles(3)=acos((vecs(1,1)*vecs(2,1)+vecs(1,2)*vecs(2,2)
!     &          +vecs(1,3)*vecs(2,3))/(alat(1)*alat(2)))*180.0D0/Pi
!      vecs=0.0D0
!      vecs(1,1)=alat(1)
!      vecs(2,1)=alat(1)*cos(angles(3)*Pi/180.0D0)
!      vecs(2,2)=alat(2)*sin(angles(3)*Pi/180.0D0)
!      vecs(3,1)=alat(3)*cos(angles(2)*Pi/180.0D0)
!      vecs(3,2)=alat(3)*(cos(angles(1)*Pi/180.0D0)
!     &       -cos(angles(3)*Pi/180.0D0)*cos(angles(2)*Pi/180.0D0))
!     &       /sin(angles(3)*Pi/180.0D0)
!      vecs(3,3)=sqrt(alat(3)**2-vecs(3,1)**2-vecs(3,2)**2)
!
!      ! write
!      call write_coords(foutcoords,finform,newatoms,nnewatoms,
!     &                  species,nspecies,vecs)
!
!      if(talk)print fsubendext, 'transformmatrix'
!
!      end subroutine
!
c---------------------------------------------------------------------

      subroutine transformmatrix(matrix,atoms,vecs,newatoms,              &
     &                           newvecs,nsc)
     
      use defs
      use misc
      implicit none

      ! variables
      double precision :: matrix(1:3,1:3)
      type(atom) :: atoms(:)
      type(atom), allocatable :: newatoms(:)
      type(element), allocatable :: species(:)
      integer :: natoms,nspecies,nnewatoms,nsc(1:3)
      double precision, intent(in) :: vecs(3,3)
      double precision :: newvecs(3,3)
 
      ! local variables
      double precision :: imatrix(1:3,1:3)
      integer iatom,i,j,k,n,p,q,r,s,nuc(1:3),iallatom,nallatoms
      integer ispecies
      double precision xold(1:3),xnew(1:3)
      double precision inmat(1:3,1:3)
      double precision invinmat(1:3,1:3),inewvecs(1:3,1:3)
      double precision :: tmatrix(1:3,1:3)
      double precision :: cellpars(1:6)
      logical,allocatable :: delete(:)
      type(atom), allocatable :: allatoms(:)

      talk=.true.
      if(talk)print fsubstart,"transformmatrix"
      !
      natoms=size(atoms)
      
      !! take also atoms of neighboring old unit cells into account, if
      !! the new cell lies partly outside the old one
      !nuc(1)=3
      !nuc(2)=3
      !nuc(3)=3
      nuc(1)=1
      nuc(2)=1
      nuc(3)=1
      
      ! take into account desired supercell dimension
      do n=1,3
        matrix(n,1:3)=matrix(n,1:3)/dble(nsc(n))
        !vecs(n,1:3)=vecs(n,1:3)*dble(nsc(n))
      end do  
      nallatoms=natoms*nuc(1)*nuc(2)*nuc(3)*8
      allocate(allatoms(1:nallatoms)) ! natoms*nuc(1)*nuc(2)*nuc(3)))
      !allocate(x_new(1:natom*nuc(1)*nuc(2)*nuc(3),1:3))
      allocate(delete(1:nallatoms)) !natoms*nuc(1)*nuc(2)*nuc(3)))
      
      ! multiply with matrix
      iallatom=0
      do iatom=1,natoms
            xold=atoms(iatom)%where
            xnew=0.0D0
            do q=-nuc(1),nuc(1)-1  !0
              do r=-nuc(2),nuc(2)-1 !0
                do s=-nuc(3),nuc(3)-1 !0
                  iallatom=iallatom+1
                  xnew(1)=(xold(1)+float(q))*matrix(1,1)
     &                      +(xold(2)+float(r))*matrix(1,2)
     &                      +(xold(3)+float(s))*matrix(1,3)
                  xnew(2)=(xold(1)+float(q))*matrix(2,1)
     &                      +(xold(2)+float(r))*matrix(2,2)
     &                      +(xold(3)+float(s))*matrix(2,3)
                  xnew(3)=(xold(1)+float(q))*matrix(3,1)
     &                      +(xold(2)+float(r))*matrix(3,2)
     &                      +(xold(3)+float(s))*matrix(3,3)
                  ! shift back to between 0 and 1
                  do j=1,3
                    do while(xnew(j).ge.1.0D0-tolfrac(j))
                      xnew(j)=xnew(j)-1.0D0
                    end do
                    do while(xnew(j).lt.-tolfrac(j))
                      xnew(j)=xnew(j)+1.0D0
                    end do
                  end do
                  !do j=1,3
                  !  do k=1,3
                  !    xnew(j)=xnew(j)+matrix(j,k)*xold(k)
                  !  end do
                  !end do
                  allatoms(iallatom)=atoms(iatom)
                  allatoms(iallatom)%where=xnew
                end do
              end do
            end do
      end do
      
      ! delete doubled coordinates
      delete=.False.
      do n=1,nallatoms-1 !natoms*nuc(1)*nuc(2)*nuc(3)-1
            do p=n+1,nallatoms !natoms*nuc(1)*nuc(2)*nuc(3)
              if(abs(allatoms(p)%where(1)-allatoms(n)%where(1))
     &          .le.1E-6.and.abs(allatoms(p)%where(2)
     &          -allatoms(n)%where(2)).le.1E-6.and.
     &               abs(allatoms(p)%where(3)-allatoms(n)%where(3))
     &               .le.1E-6.and.
     &                  allatoms(p)%core.eq.allatoms(n)%core) then

                        delete(p)=.True.
              end if
            end do
      end do
!      ! delete doubled coordinates
!      delete=.False.
!      do n=1,nallatoms-1 !natoms*nuc(1)*nuc(2)*nuc(3)-1
!            do p=n+1,nallatoms !natoms*nuc(1)*nuc(2)*nuc(3)
!              if(abs(allatoms(p)%where(1)-allatoms(n)%where(1))
!     &          .le.tolfrac(1).and.abs(allatoms(p)%where(2)
!     &          -allatoms(n)%where(2)).le.tolfrac(2).and.
!     &               abs(allatoms(p)%where(3)-allatoms(n)%where(3))
!     &               .le.tolfrac(3)) then
!                        delete(p)=.True.
!              end if
!            end do
!      end do


      ! find final number of coordinates
      nnewatoms=0
      do n=1,nallatoms !natoms*nuc(1)*nuc(2)*nuc(3)
            if(delete(n).eqv..False.) nnewatoms=nnewatoms+1
      end do
      ! get new atoms
      allocate(newatoms(nnewatoms))
      iatom=0
      do n=1,nallatoms !natoms*nuc(1)*nuc(2)*nuc(3)
            if (delete(n).eqv..False.) then
                  iatom=iatom+1
                  newatoms(iatom)=allatoms(n)
            end if
      end do
      !
      !
      !
      ! transform vectors
      call transpon_matrix(vecs,inmat)
      !print*, "Vi:"
      !do j=1,3
      !  print*, inmat(j,1:3)
      !end do
      call invert_matrix(inmat,invinmat)
      !print*, "(Vi)^-1:"
      !do j=1,3
      !  print*, invinmat(j,1:3)
      !end do
      !call transpon_matrix(matrix,tmatrix)
      call mymatmul(matrix,invinmat,inewvecs)
      !print*, "M*(Vi)^-1:"
      !do j=1,3
      !  print*, inewvecs(j,1:3)
      !end do
      call invert_matrix(inewvecs,newvecs)
      print*, "(M*(Vi)^-1)^-1:"
      !do j=1,3
      !  print*, newvecs(j,1:3)
      !end do
      call transpon_matrix(newvecs,inewvecs)
      !print*, "((M*(Vi)^-1)^-1)^t:"
      !do j=1,3
      !  print*, inewvecs(j,1:3)
      !end do
      newvecs=inewvecs
      ! rotate vectors so that they are conform with lammps
      call vecs2cellpars(newvecs,cellpars)
      call cellpars2vecs(cellpars, newvecs)
      !print*,"vecs:",newvecs
      !
      ! get correct numbers of atoms of each species of final structure:  
      call getspecies(newatoms,species)
      !
      if(talk)print fsubendext, 'transformmatrix'
      !
      end subroutine transformmatrix

c---------------------------------------------------------------------
!
!      subroutine transformshift(shift,fincoords,foutcoords,finform)
!     
!      use defs
!      use readcoords
!      use writecoords
!      implicit none
!
!      ! variables
!      double precision :: shift(1:3)
!      character(len=*) :: fincoords,foutcoords,finform
!      type(atom), allocatable :: atoms(:)
!      type(element), allocatable :: species(:)
!      integer :: natoms,nspecies
!
!      ! local variables
!      integer iatom,j
!      double precision xold(1:3),xnew(1:3), vecs(1:3,1:3)
!
!      talk=.true.
!      if(talk)print '(" ")'
!      if(talk)print fsubstart,"transformshift"
!      
!      ! read in  coordinates ...
!      if(talk)      print '(8x,"calling read_coords...")'
!            call read_coords(fincoords,finform,atoms,natoms,species,
!     &           nspecies,vecs)
!      
!      ! shift atoms
!      do iatom=1,natoms
!            xold=atoms(iatom)%where
!            do j=1,3
!              xnew(j)=xold(j)+shift(j)
!            end do
!            atoms(iatom)%where=xnew
!      end do
!      
!      ! write
!      if(talk)print '(8x,"writing coordinates ...")'
!      call write_coords(foutcoords,finform,atoms,natoms,species,
!     &                  nspecies,vecs)
!
!      if(talk)print fsubendext, 'transformshift'
!
!      end subroutine
!
c---------------------------------------------------------------------

      subroutine transformshift(shift,atoms,vecs,newatoms)
      ! shifts by relative vector
     
      use defs
      implicit none

      ! variables
      type(atom), intent(in) :: atoms(:)
      type(atom), allocatable :: newatoms(:)
      double precision, intent(in) :: shift(1:3),vecs(3,3)

      ! local variables
      integer iatom
      integer :: natoms

      talk=.true.
      if(talk)print '(" ")'
      if(talk)print fsubstart,"transformshift"
      
      natoms=size(atoms)
      allocate(newatoms(natoms))
      newatoms=atoms
      ! shift atoms
      do iatom=1,natoms
            newatoms(iatom)%where=atoms(iatom)%where+shift
            newatoms(iatom)%abswhere=atoms(iatom)%abswhere                &
     &                 +shift(1)*vecs(1,1:3)+shift(2)*vecs(2,1:3)         &
     &                 +shift(3)*vecs(3,1:3)
      end do
      ! 
      if(talk)print fsubendext, 'transformshift'
      !
      end subroutine transformshift

c---------------------------------------------------------------------

      !subroutine randisp(infile,informat,outfile,outformat,maxdisp)
      subroutine randisp(atoms,vecs,maxdisp)
      ! apply random displacements
     
      use defs
      !use readcoords
      !use writecoords
      use misc, only : abs2frac,frac2abs,init_random_seed
      implicit none

      ! variables
      double precision, intent(in) :: maxdisp(1:3)
      !character(len=*), intent(in) :: informat,infile,outformat,outfile
      type(atom) :: atoms(:)
      double precision :: vecs(1:3,1:3)

      ! local variables
      type(element), allocatable :: species(:)
      integer iatom
      integer :: natoms,nspecies
      double precision shift(1:3)

      talk=.true.
      if(talk)print '(" ")'
      if(talk)print fsubstart,"randisp"
      
      !call read_coords(infile,informat,atoms,natoms,species,
      !&                       nspecies,vecs)
      natoms=size(atoms)
      ! shift atoms
      if (all(maxdisp.ge.0.0d0)) then
        print '(8x,"Shifting by fractional random number up to ",3(F12.8  &
     &))',maxdisp(1:3)
        do iatom=1,natoms
            ! get seed
            call init_random_seed()
            ! get random number
    	    call random_number(shift(1))
            ! get seed
            call init_random_seed()
            ! get random number
    	    call random_number(shift(2))
            ! get seed
            call init_random_seed()
            ! get random number
    	    call random_number(shift(3))
            ! random numbers are between 0 and 1. Shift by -0.5 to get
            ! random numbers between -0.5 and 0.5. Multiply with 2 to
            ! get random numbers between -1 and 1. Finally scale with
            ! maxdisp to get random nnumbers between -maxdisp and
            ! maxdisp.
            shift(1)=(shift(1)-0.5d0)*2.0d0*maxdisp(1)  
            shift(2)=(shift(2)-0.5d0)*2.0d0*maxdisp(2)  
            shift(3)=(shift(3)-0.5d0)*2.0d0*maxdisp(3)  
            atoms(iatom)%where=atoms(iatom)%where+shift
            call frac2abs(atoms(iatom)%where,vecs,atoms(iatom)%abswhere)
        end do
      else
        print '(8x,"Shifting by absolute random number up to ",3(F12.8))  &
     &',  abs(maxdisp(1)),abs(maxdisp(2)),abs(maxdisp(3))
        do iatom=1,natoms
            ! get seed
            call init_random_seed()
            ! get random number
    	    call random_number(shift(1))
            ! get seed
            call init_random_seed()
            ! get random number
    	    call random_number(shift(2))
            ! get seed
            call init_random_seed()
            ! get random number
    	    call random_number(shift(3))
            ! random numbers are between 0 and 1. Shift by -0.5 to get
            ! random numbers between -0.5 and 0.5. Multiply with 2 to
            ! get random numbers between -1 and 1. Finally scale with
            ! maxdisp to get random nnumbers between -maxdisp and
            ! maxdisp.
            shift(1)=(shift(1)-0.5d0)*2.0d0*abs(maxdisp(1))  
            shift(2)=(shift(2)-0.5d0)*2.0d0*abs(maxdisp(2))  
            shift(3)=(shift(3)-0.5d0)*2.0d0*abs(maxdisp(3))  
            atoms(iatom)%abswhere=atoms(iatom)%abswhere+shift             &
            call abs2frac(atoms(iatom)%abswhere,vecs,atoms(iatom)%where)
        end do
      end if
!      call write_coords(outfile,outformat,atoms,natoms,species,
!     &                       nspecies,vecs)
      ! 
      if(talk)print fsubendext, 'randisp'
      !
      end subroutine randisp

c---------------------------------------------------------------------

      subroutine transformshiftabs(shift,atoms,vecs,newatoms)
      ! shifts by absolute vector
     
      use defs
      use misc, only : abs2frac
      implicit none

      ! variables
      type(atom), intent(in) :: atoms(:)
      type(atom), allocatable :: newatoms(:)
      double precision, intent(in) :: shift(1:3),vecs(3,3)

      ! local variables
      integer iatom
      integer :: natoms

      talk=.true.
      if(talk)print '(" ")'
      if(talk)print fsubstart,"transformshiftabs"
      
      natoms=size(atoms)
      allocate(newatoms(natoms))
      newatoms=atoms
      ! shift atoms
      do iatom=1,natoms
        newatoms(iatom)%abswhere=atoms(iatom)%abswhere+shift
        call abs2frac(newatoms(iatom)%abswhere,vecs,                      &
     &        newatoms(iatom)%where)
      end do
      ! 
      if(talk)print fsubendext, 'transformshiftabs'
      !
      end subroutine transformshiftabs

c---------------------------------------------------------------------
!
!      subroutine transformsort(fincoords,foutcoords,finform,
!     &                         tol)
!     
!      use defs
!      use readcoords
!      use writecoords
!      implicit none
!
!      ! variables
!      character(len=*) :: fincoords,foutcoords,finform
!      type(atom), allocatable :: atoms(:)
!      type(element), allocatable :: species(:)
!      integer :: natoms,nspecies
!      double precision tol(1:3)
!
!      ! local variables
!      integer iatom,jatom,i,j,k,l
!      type(atom) atomstore
!      double precision vecs(1:3,1:3)
!
!      talk=.true.
!      if(talk)print '(" ")'
!      if(talk)print fsubstart,"transformsort"
!      
!      ! read in  coordinates ...
!      if(talk)      print '(8x,"calling read_coords...")'
!      call read_coords(fincoords,finform,atoms,natoms,species,
!     &           nspecies,vecs)
!      
!      if(talk)print '(8x,"sorting coordinates...")'
!      ! sort atoms with respect to x1 coordinate
!      do iatom=1,natoms
!        do jatom=iatom+1,natoms
!          if (atoms(iatom)%name.eq.atoms(jatom)%name.and.
!     &        atoms(iatom)%core.eq.atoms(jatom)%core) then
!                if (atoms(jatom)%where(1).lt.atoms(iatom)%where(1)
!     &              -tol(1)) then
!                      atomstore=atoms(iatom)
!                      atoms(iatom)=atoms(jatom)
!                      atoms(jatom)=atomstore
!                end if
!          end if
!        end do
!      end do
!      ! sort atoms with respect to x2 coordinate
!      do iatom=1,natoms
!        do jatom=iatom+1,natoms
!          if (atoms(iatom)%name.eq.atoms(jatom)%name.and.
!     &        atoms(iatom)%core.eq.atoms(jatom)%core) then
!                if (abs(atoms(jatom)%where(1)-atoms(iatom)%where(1))
!     &              .le.tol(1).and.atoms(jatom)%where(2).lt.
!     &              atoms(iatom)%where(2)-tol(2)) then
!                      atomstore=atoms(iatom)
!                      atoms(iatom)=atoms(jatom)
!                      atoms(jatom)=atomstore
!                end if
!          end if
!        end do
!      end do
!      ! sort atoms with respect to x3 coordinate
!      do iatom=1,natoms
!        do jatom=iatom+1,natoms
!          if (atoms(iatom)%name.eq.atoms(jatom)%name.and.
!     &        atoms(iatom)%core.eq.atoms(jatom)%core) then
!                if (abs(atoms(jatom)%where(1)-atoms(iatom)%where(1))
!     &              .le.tol(1).and.abs(atoms(jatom)%where(2)
!     &              -atoms(iatom)%where(2)).le.tol(2).and.
!     &              atoms(jatom)%where(3).lt.atoms(iatom)%where(3)
!     &              -tol(3)) then
!                      atomstore=atoms(iatom)
!                      atoms(iatom)=atoms(jatom)
!                      atoms(jatom)=atomstore
!                end if
!          end if
!        end do
!      end do
!      
!      ! write
!      if(talk)print '(8x,"writing coordinates ...")'
!      call write_coords(foutcoords,finform,atoms,natoms,species,
!     &                  nspecies,vecs)
!
!      if(talk)print fsubendext, 'transformsort'
!
!      end subroutine
!
c---------------------------------------------------------------------

      subroutine transformsort(atoms,vecs,tol)
     
      use defs
      implicit none

      ! variables
      type(atom), intent(inout) :: atoms(:)
      double precision, intent(in) :: tol(1:3)
      double precision, intent(in) :: vecs(1:3,1:3)

      ! local variables
      integer :: natoms
      integer iatom,jatom,i,j,k,l
      type(atom) atomstore

      if(talk)print fsubstart,"transformsort"
      
      if(talk)print '(8x,"sorting coordinates...")'
      !
      natoms=size(atoms)
      ! sort atoms with respect to x1 coordinate
      do iatom=1,natoms
        do jatom=iatom+1,natoms
          if (atoms(iatom)%name.eq.atoms(jatom)%name.and.
     &        atoms(iatom)%core.eq.atoms(jatom)%core) then
                if (atoms(jatom)%where(1).lt.atoms(iatom)%where(1)
     &              -tol(1)) then
                      atomstore=atoms(iatom)
                      atoms(iatom)=atoms(jatom)
                      atoms(jatom)=atomstore
                end if
          end if
        end do
      end do
      ! sort atoms with respect to x2 coordinate
      do iatom=1,natoms
        do jatom=iatom+1,natoms
          if (atoms(iatom)%name.eq.atoms(jatom)%name.and.
     &      atoms(iatom)%core.eq.atoms(jatom)%core) then
              if (abs(atoms(jatom)%where(1)-atoms(iatom)%where(1))
     &              .le.tol(1).and.atoms(jatom)%where(2).lt.
     &              atoms(iatom)%where(2)-tol(2)) then
                      atomstore=atoms(iatom)
                      atoms(iatom)=atoms(jatom)
                      atoms(jatom)=atomstore
              end if
          end if
        end do
      end do
      ! sort atoms with respect to x3 coordinate
      do iatom=1,natoms
        do jatom=iatom+1,natoms
          if (atoms(iatom)%name.eq.atoms(jatom)%name.and.
     &        atoms(iatom)%core.eq.atoms(jatom)%core) then
                if (abs(atoms(jatom)%where(1)-atoms(iatom)%where(1))
     &              .le.tol(1).and.abs(atoms(jatom)%where(2)
     &              -atoms(iatom)%where(2)).le.tol(2).and.
     &              atoms(jatom)%where(3).lt.atoms(iatom)%where(3)
     &              -tol(3)) then
                      atomstore=atoms(iatom)
                      atoms(iatom)=atoms(jatom)
                      atoms(jatom)=atomstore
                end if
          end if
        end do
      end do
      ! 
      if(talk)print fsubendext, 'transformsort'
      !
      end subroutine transformsort

c---------------------------------------------------------------------
!
!      subroutine transformmapback(fincoords,foutcoords,finform,
!     &                           tol)
!     
!      use defs
!      use readcoords
!      use writecoords
!      implicit none
!
!      ! variables
!      double precision :: tol(1:3)
!      character(len=*) :: fincoords,foutcoords,finform
!      type(atom), allocatable :: atoms(:)
!      type(element), allocatable :: species(:)
!      integer :: natoms,nspecies
!
!      ! local variables
!      integer iatom,j
!      double precision x(1:3), vecs(1:3,1:3)
!
!      talk=.true.
!      if(talk)print '(" ")'
!      if(talk)print fsubstart,"transformmapback"
!      
!      ! read in  coordinates ...
!      if(talk) print '(8x,"calling read_coords...")'
!            call read_coords(fincoords,finform,atoms,natoms,species,
!     &           nspecies,vecs)
!      
!      ! shift atoms
!      if(talk)print '(8x,"mapping back atoms...")'
!      do iatom=1,natoms
!            x=atoms(iatom)%where
!            do j=1,3
!              do while (x(j).ge.1.0D0-tol(j))
!                x(j)=x(j)-1.0D0
!              end do
!              do while (x(j).lt.0.0D0-tol(j))
!                x(j)=x(j)+1.0D0
!              end do
!            end do
!            atoms(iatom)%where=x
!      end do
!      
!      ! write
!      if(talk)print '(8x,"writing coordinates ...")'
!      call write_coords(foutcoords,finform,atoms,natoms,species,
!     &                  nspecies,vecs)
!
!      if(talk)print fsubendext, 'transformmapback'
!
!      end subroutine
!
c---------------------------------------------------------------------

      subroutine transformmapback(atoms,vecs,tol)
     
      use defs
      implicit none

      ! variables
      double precision :: tol(1:3)
      type(atom), intent(inout) :: atoms(:)
      double precision, intent(in) :: vecs(3,3)

      ! local variables
      integer natoms,iatom,j

      if(talk)print fsubstart,"transformmapback"
      
      ! shift atoms
      natoms=size(atoms)
      if(talk)print '(8x,"mapping back atoms...")'
      do iatom=1,natoms
        do j=1,3
          do while (atoms(iatom)%where(j).ge.1.0D0-tol(j))
            atoms(iatom)%where(j)=atoms(iatom)%where(j)-1.0d0
            atoms(iatom)%abswhere=atoms(iatom)%abswhere-vecs(j,1:3)
          end do
          do while (atoms(iatom)%where(j).lt.0.0D0-tol(j))
            atoms(iatom)%where(j)=atoms(iatom)%where(j)+1.0d0
            atoms(iatom)%abswhere=atoms(iatom)%abswhere+vecs(j,1:3)
          end do 
        end do ! j
      end do ! iatom
      ! 
      if(talk)print fsubendext, 'transformmapback'
      !
      end subroutine transformmapback

c---------------------------------------------------------------------

      subroutine transformmirror(atoms,vecs,axisint)
     
      use defs
      use misc, only : frac2abs
      implicit none

      ! variables
      integer :: axisint
      type(atom), intent(inout) :: atoms(:)
      double precision, intent(in) :: vecs(3,3)

      ! local variables
      integer natoms,iatom

      if(talk)print fsubstart,"transformmirror"
      
      ! shift atoms
      natoms=size(atoms)
      if(talk)print '(8x,"mirroring atoms...")'
      atoms%where(axisint)=-atoms%where(axisint)
      do iatom=1,natoms
        call frac2abs(atoms(iatom)%where,vecs,atoms(iatom)%abswhere)
      end do ! iatom
      ! 
      if(talk)print fsubendext, 'transformmirror'
      !
      end subroutine transformmirror

c---------------------------------------------------------------------
!
!      subroutine matrixmult(fincoords,foutcoords,finform,matrix)
!     
!      use defs
!      use readcoords
!      use writecoords
!      use misc
!      implicit none
!
!      ! variables
!      double precision :: matrix(1:3,1:3)!,imatrix(1:3,1:3)
!      character(len=*) :: fincoords,foutcoords,finform
! 
!      ! local variables
!      type(atom), allocatable :: atoms(:),newatoms(:)
!      type(element), allocatable :: species(:)
!      integer :: natoms,nspecies
!      integer iatom,i,j,k,n
!      integer ispecies
!      double precision xold(1:3),xnew(1:3),vecs(1:3,1:3)
!      !double precision newvecs(1:3,1:3),inmat(1:3,1:3)
!      !double precision invinmat(1:3,1:3),inewvecs(1:3,1:3)
!      !double precision :: tmatrix(1:3,1:3)
!      !double precision alat(1:3),angles(1:3)
!
!      if(talk)print fsubstart,"matrixmult"
!      
!      ! read in  coordinates ...
!      if(talk)print '(" ")'
!      if(talk)print '(8x,"calling read_coords...")'
!      call read_coords(fincoords,finform,atoms,natoms,species,
!     &     nspecies,vecs)
!      if(talk)print '(8x,"Coords read in.")'
!      
!      ! multiply with matrix
!      do iatom=1,natoms
!            ! obtains cartesian coordinates from fractional ones
!            call frac2abs(atoms(iatom)%where,vecs,xold)
!            ! now apply the matrix to the cartesian coordinates
!            xnew=0.0D0
!            xnew(1)=xold(1)*matrix(1,1)
!     &                +xold(2)*matrix(1,2)
!     &                +xold(3)*matrix(1,3)
!            xnew(2)=xold(1)*matrix(2,1)
!     &                +xold(2)*matrix(2,2)
!     &                +xold(3)*matrix(2,3)
!            xnew(3)=xold(1)*matrix(3,1)
!     &                +xold(2)*matrix(3,2)
!     &                +xold(3)*matrix(3,3)
!            call abs2frac(xnew,vecs,atoms(iatom)%where)
!      end do
!      
!
!!      ! transform vectors
!!      call transpon_matrix(vecs,inmat)
!!!      print*, "Vi:"
!!!      do j=1,3
!!!        print*, inmat(j,1:3)
!!!      end do
!!      call invert_matrix(inmat,invinmat)
!!!      print*, "(Vi)^-1:"
!!!      do j=1,3
!!!        print*, invinmat(j,1:3)
!!!      end do
!!      !call transpon_matrix(matrix,tmatrix)
!!      call mymatmul(matrix,invinmat,inewvecs)
!!!      print*, "M*(Vi)^-1:"
!!!      do j=1,3
!!!        print*, inewvecs(j,1:3)
!!!      end do
!!      call invert_matrix(inewvecs,newvecs)
!!!      print*, "(M*(Vi)^-1)^-1:"
!!!      do j=1,3
!!!        print*, newvecs(j,1:3)
!!!      end do
!!      call transpon_matrix(newvecs,inewvecs)
!!!      print*, "((M*(Vi)^-1)^-1)^t:"
!!!      do j=1,3
!!!        print*, inewvecs(j,1:3)
!!!      end do
!!      vecs=inewvecs
!!      ! rotate vectors so that they are conform with lammps
!!      do i=1,3
!!        alat(i)=sqrt(vecs(i,1)**2+vecs(i,2)**2+vecs(i,3)**2)
!!      end do
!!      angles(1)=acos((vecs(2,1)*vecs(3,1)+vecs(2,2)*vecs(3,2)
!!     &          +vecs(2,3)*vecs(3,3))/(alat(2)*alat(3)))*180.0D0/Pi
!!      angles(2)=acos((vecs(1,1)*vecs(3,1)+vecs(1,2)*vecs(3,2)
!!     &          +vecs(1,3)*vecs(3,3))/(alat(1)*alat(3)))*180.0D0/Pi
!!      angles(3)=acos((vecs(1,1)*vecs(2,1)+vecs(1,2)*vecs(2,2)
!!     &          +vecs(1,3)*vecs(2,3))/(alat(1)*alat(2)))*180.0D0/Pi
!!      vecs=0.0D0
!!      vecs(1,1)=alat(1)
!!      vecs(2,1)=alat(1)*cos(angles(3)*Pi/180.0D0)
!!      vecs(2,2)=alat(2)*sin(angles(3)*Pi/180.0D0)
!!      vecs(3,1)=alat(3)*cos(angles(2)*Pi/180.0D0)
!!      vecs(3,2)=alat(3)*(cos(angles(1)*Pi/180.0D0)
!!     &       -cos(angles(3)*Pi/180.0D0)*cos(angles(2)*Pi/180.0D0))
!!     &       /sin(angles(3)*Pi/180.0D0)
!!      vecs(3,3)=sqrt(alat(3)**2-vecs(3,1)**2-vecs(3,2)**2)
!
!      ! write
!      call write_coords(foutcoords,finform,atoms,natoms,
!     &                  species,nspecies,vecs)
!
!      if(talk)print fsubendext, 'matrixmult'
!
!      end subroutine
!
c---------------------------------------------------------------------

      subroutine matrixmult(atoms,vecs,matrix)
      ! multiplies absolute coordinates with a given matrix.
     
      use defs
      use misc
      implicit none

      ! variables
      double precision, intent(in) :: matrix(1:3,1:3)
      type(atom), intent(inout) :: atoms(:)
      double precision, intent(inout) :: vecs(3,3)
 
      ! local variables
      !type(element), allocatable :: species(:)
      integer :: natoms!,nspecies
      integer iatom,i,j,k,n
      !integer ispecies
      double precision xold(1:3),xnew(1:3)

      if(talk)print fsubstart,"matrixmult"
      
      if(talk)print '(8x,"Coords read in.")'
      
      ! multiply with matrix
      do iatom=1,natoms
            ! obtains cartesian coordinates from fractional ones
            !call frac2abs(atoms(iatom)%where,vecs,xold)
            xold=atoms(iatom)%abswhere
            ! now apply the matrix to the cartesian coordinates
            xnew=0.0D0
            xnew(1)=xold(1)*matrix(1,1)
     &                +xold(2)*matrix(1,2)
     &                +xold(3)*matrix(1,3)
            xnew(2)=xold(1)*matrix(2,1)
     &                +xold(2)*matrix(2,2)
     &                +xold(3)*matrix(2,3)
            xnew(3)=xold(1)*matrix(3,1)
     &                +xold(2)*matrix(3,2)
     &                +xold(3)*matrix(3,3)
            atoms(iatom)%abswhere=xnew
            call abs2frac(xnew,vecs,atoms(iatom)%where)
      end do
      
      if(talk)print fsubendext, 'matrixmult'

      end subroutine matrixmult

c---------------------------------------------------------------------
!
!      subroutine rotatecoords(fincoords,foutcoords,finform,axisint,
!     &                  angle)
!     
!      use defs
!      use readcoords
!      use writecoords
!      use misc
!      implicit none
!
!      ! variables
!      character(len=*) :: fincoords,foutcoords,finform
!      integer axisint
!      double precision angle
!
!      ! local variables
!      double precision :: matrix(1:3,1:3)
!      type(atom), allocatable :: atoms(:),newatoms(:)
!      type(element), allocatable :: species(:)
!      integer :: natoms,nspecies
!      integer iatom,i,j,k,n
!      integer ispecies
!      double precision xold(1:3),xnew(1:3),vecs(1:3,1:3)
!
!      talk=.true.
!      if(talk)print fsubstart,"rotatecoords"
!      
!      ! read in  coordinates ...
!      if(talk)print '(" ")'
!      if(talk)print '(8x,"calling read_coords...")'
!      call read_coords(fincoords,finform,atoms,natoms,species,
!     &     nspecies,vecs)
!      if(talk)print '(8x,"Coords read in.")'
!      
!      ! biuld the rotation matrix
!      matrix=0.0D0
!      do i=1,3
!        matrix(i,i)=1.0D0
!      end do
!      select case(axisint)
!      case(1)
!        matrix(2,2)=cos(angle*Pi/180.0D0)
!        matrix(2,3)=-sin(angle*Pi/180.0D0)
!        matrix(3,2)=sin(angle*Pi/180.0D0)
!        matrix(3,3)=cos(angle*Pi/180.0D0)
!      case(2)
!        matrix(1,1)=cos(angle*Pi/180.0D0)
!        matrix(1,3)=sin(angle*Pi/180.0D0)
!        matrix(3,1)=-sin(angle*Pi/180.0D0)
!        matrix(3,3)=cos(angle*Pi/180.0D0)
!      case(3)
!        matrix(1,1)=cos(angle*Pi/180.0D0)
!        matrix(1,2)=-sin(angle*Pi/180.0D0)
!        matrix(2,1)=sin(angle*Pi/180.0D0)
!        matrix(2,2)=cos(angle*Pi/180.0D0)
!      case default
!        goto 1000
!      end select
!      ! multiply coordinates with matrix
!      do iatom=1,natoms
!            ! obtains cartesian coordinates from fractional ones
!            call frac2abs(atoms(iatom)%where,vecs,xold)
!            ! now apply the matrix to the cartesian coordinates
!            xnew=0.0D0
!            xnew(1)=xold(1)*matrix(1,1)
!     &                +xold(2)*matrix(1,2)
!     &                +xold(3)*matrix(1,3)
!            xnew(2)=xold(1)*matrix(2,1)
!     &                +xold(2)*matrix(2,2)
!     &                +xold(3)*matrix(2,3)
!            xnew(3)=xold(1)*matrix(3,1)
!     &                +xold(2)*matrix(3,2)
!     &                +xold(3)*matrix(3,3)
!            call abs2frac(xnew,vecs,atoms(iatom)%where)
!      end do
!      
!      ! write
!      call write_coords(foutcoords,finform,atoms,natoms,
!     &                  species,nspecies,vecs)
!
!      if(talk)print fsubendext, 'rotatecoords'
!      return
!
!1000  nerr=nerr+1
!      print ferrmssg, "Please specify a cartesian axis (1,2, or 3)"
!      return
!
!      end subroutine
!
c---------------------------------------------------------------------

      subroutine rotatecoords(atoms,vecs,axisint,angle,vector)
     
      use defs
      use misc
      implicit none

      ! variables
      integer, intent(in) :: axisint
      double precision, intent(in) :: angle
      type(atom), intent(inout) :: atoms(:)
      type(atom), allocatable :: newatoms(:)
      double precision, intent(in) :: vecs(3,3)
      double precision, optional :: vector(1:3)

      ! local variables
      double precision :: matrix(1:3,1:3)
      integer :: natoms
      integer iatom,i,j,k,n
      double precision xold(1:3),xnew(1:3),anglerad

      if(talk)print fsubstart,"rotatecoords"

      natoms=size(atoms)
      
      anglerad=angle*Pi/180.0d0
      ! biuld the rotation matrix
      matrix=0.0D0
      do i=1,3
        matrix(i,i)=1.0D0
      end do
      select case(axisint)
      case(0)
        if (.not.present(vector)) then
          call error_stop('no rotation axis specified')
        end if
        vector=vector/norm2(vector)
        matrix(1,1)=cos(anglerad)+vector(1)**2*(1.0d0-cos(anglerad))
        matrix(1,2)=vector(1)*vector(2)*(1.0d0-cos(anglerad))             &
     &            -vector(3)*sin(anglerad)
        matrix(1,3)=vector(1)*vector(3)*(1.0d0-cos(anglerad))             &
     &            +vector(2)*sin(anglerad)
        matrix(2,1)=vector(2)*vector(1)*(1.0d0-cos(anglerad))             &
     &            +vector(3)*sin(anglerad)
        matrix(2,2)=cos(anglerad)+vector(2)**2*(1.0d0-cos(anglerad)) 
        matrix(2,3)=vector(2)*vector(3)*(1.0d0-cos(anglerad))             &
     &            -vector(1)*sin(anglerad)
        matrix(3,1)=vector(3)*vector(1)*(1.0d0-cos(anglerad))             &
     &            -vector(2)*sin(anglerad)
        matrix(3,2)=vector(3)*vector(2)*(1.0d0-cos(anglerad))             &
     &            +vector(1)*sin(anglerad)
        matrix(3,3)=cos(anglerad)+vector(3)**2*(1.0d0-cos(anglerad)) 
      case(1)
        matrix(2,2)=cos(angle*Pi/180.0D0)
        matrix(2,3)=-sin(angle*Pi/180.0D0)
        matrix(3,2)=sin(angle*Pi/180.0D0)
        matrix(3,3)=cos(angle*Pi/180.0D0)
      case(2)
        matrix(1,1)=cos(angle*Pi/180.0D0)
        matrix(1,3)=sin(angle*Pi/180.0D0)
        matrix(3,1)=-sin(angle*Pi/180.0D0)
        matrix(3,3)=cos(angle*Pi/180.0D0)
      case(3)
        matrix(1,1)=cos(angle*Pi/180.0D0)
        matrix(1,2)=-sin(angle*Pi/180.0D0)
        matrix(2,1)=sin(angle*Pi/180.0D0)
        matrix(2,2)=cos(angle*Pi/180.0D0)
      case default
        goto 1000
      end select
      ! print rotation matrix 
      print '(8x,"rotation matrix: ",9(F15.9,x))',                        &
     &    matrix(1,:),matrix(2,:),matrix(3,:)
      ! multiply coordinates with matrix
      do iatom=1,natoms
            ! obtains cartesian coordinates from fractional ones
            !call frac2abs(atoms(iatom)%where,vecs,xold)
            xold=atoms(iatom)%abswhere
            ! now apply the matrix to the cartesian coordinates
            xnew=0.0D0
            xnew(1)=xold(1)*matrix(1,1)
     &                +xold(2)*matrix(1,2)
     &                +xold(3)*matrix(1,3)
            xnew(2)=xold(1)*matrix(2,1)
     &                +xold(2)*matrix(2,2)
     &                +xold(3)*matrix(2,3)
            xnew(3)=xold(1)*matrix(3,1)
     &                +xold(2)*matrix(3,2)
     &                +xold(3)*matrix(3,3)
            atoms(iatom)%abswhere=xnew
            call abs2frac(xnew,vecs,atoms(iatom)%where)
      end do
      
      if(talk)print fsubendext, 'rotatecoords'
      return

1000  nerr=nerr+1
      print ferrmssg, "Please specify a cartesian axis (1,2, or 3)"
      return

      end subroutine rotatecoords

c---------------------------------------------------------------------

      end module
