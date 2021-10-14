      module mult

      implicit none

      contains

      subroutine multi(fincoords,foutcoords,finform,foutform,nsc)
      use defs
      use readcoords
      use writecoords
      use misc
      implicit none
      
      ! variables
      character(len=*) :: fincoords,foutcoords,finform,foutform
      type(atom), allocatable :: atoms(:)!,newatoms(:)
      type(element), allocatable :: species(:)
      integer :: natoms,nspecies,nsc(1:3) !,nnewatoms

      ! local variables
      integer iatom,i,j,k,n,p,q,r,s,nuc(1:3),iallatom,nallatoms
      integer ispecies
      double precision xold(1:3),xnew(1:3),vecs(1:3,1:3)
      !double precision newvecs(1:3,1:3)
      !double precision inewvecs(1:3,1:3)
      double precision alat(1:3),angles(1:3),cellpars(1:6)
      !logical,allocatable :: delete(:)
      type(atom), allocatable :: allatoms(:)

      talk=.true.
      if(talk)print fsubstart,"multi"
      
      ! read in  coordinates ...
      if(talk)      print '(" ")'
      if(talk)      print '(8x,"calling read_coords...")'
      call read_coords(fincoords,finform,atoms,natoms,species,
     &           nspecies,vecs)
      if(talk)      print '(8x,"Coords read in.")'
      
!      ! take also atoms of neighboring old unit cells into account, if
!      ! the new cell lies partly outside the old one
!      nuc(1)=3
!      nuc(2)=3
!      nuc(3)=3
      
      ! take into account desired supercell dimension
!      do n=1,3
!        matrix(n,1:3)=matrix(n,1:3)/dble(nsc(n))
!        !vecs(n,1:3)=vecs(n,1:3)*dble(nsc(n))
!      end do  
      nallatoms=natoms*nsc(1)*nsc(2)*nsc(3)
      allocate(allatoms(1:nallatoms)) ! natoms*nuc(1)*nuc(2)*nuc(3)))
      !allocate(x_new(1:natom*nuc(1)*nuc(2)*nuc(3),1:3))
!      allocate(delete(1:nallatoms)) !natoms*nuc(1)*nuc(2)*nuc(3)))
      
      ! multiply with matrix
      iallatom=0
      do iatom=1,natoms
            xold=atoms(iatom)%where
            xnew=0.0D0
            do q=0,nsc(1)-1  !0
              do r=0,nsc(2)-1 !0
                do s=0,nsc(3)-1 !0
                  iallatom=iallatom+1
                  xnew(1)=(xold(1)+float(q))/dble(nsc(1))
                  xnew(2)=(xold(2)+float(r))/dble(nsc(2))
                  xnew(3)=(xold(3)+float(s))/dble(nsc(3))
                  !xnew(1)=(xold(1)+float(q))/dble(nsc(1))
     &            !          +(xold(2)+float(r))/dble(nsc(2))
     &            !          +(xold(3)+float(s))/dble(nsc(3))
                  !xnew(2)=(xold(1)+float(q))/dble(nsc(1))
     &            !          +(xold(2)+float(r))/dble(nsc(2))
     &            !          +(xold(3)+float(s))/dble(nsc(3))
                  !xnew(3)=(xold(1)+float(q))/dble(nsc(1))
     &            !          +(xold(2)+float(r))/dble(nsc(2))
     &            !          +(xold(3)+float(s))/dble(nsc(3))
                  !! shift back to between 0 and 1
                  !do j=1,3
                  !  do while(xnew(j).ge.1.0D0-tolfrac(j))
                  !    xnew(j)=xnew(j)-1.0D0
                  !  end do
                  !  do while(xnew(j).lt.-tolfrac(j))
                  !    xnew(j)=xnew(j)+1.0D0
                  !  end do
                  !end do
                  !!do j=1,3
                  !!  do k=1,3
                  !!    xnew(j)=xnew(j)+matrix(j,k)*xold(k)
                  !!  end do
                  !!end do
                  allatoms(iallatom)=atoms(iatom)
                  allatoms(iallatom)%where=xnew
                end do
              end do
            end do
      end do
      
!      ! delete doubled coordinates
!      delete=.False.
!      do n=1,nallatoms-1 !natoms*nuc(1)*nuc(2)*nuc(3)-1
!            do p=n+1,nallatoms !natoms*nuc(1)*nuc(2)*nuc(3)
!              if(abs(allatoms(p)%where(1)-allatoms(n)%where(1))
!     &          .le.1E-6.and.abs(allatoms(p)%where(2)
!     &          -allatoms(n)%where(2)).le.1E-6.and.
!     &               abs(allatoms(p)%where(3)-allatoms(n)%where(3))
!     &               .le.1E-6) then
!                        delete(p)=.True.
!              end if
!            end do
!      end do
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
!      ! get new numbers of atoms of each species
!      do ispecies=1,nspecies
!        species(ispecies)%howmany=0
!        do iatom=1,nnewatoms
!          if (newatoms(iatom)%name(1:2).eq.species(ispecies)%name(1:2))
!     &      species(ispecies)%howmany=species(ispecies)%howmany+1
!        end do
!      end do  

      ! transform vectors
      do i=1,3
        vecs(i,1:3)=vecs(i,1:3)*dble(nsc(i))
      end do
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
      ! rotate vectors so that they are conform with lammps
      call vecs2cellpars(vecs,cellpars)
      call cellpars2vecs(cellpars,vecs) 

      ! write
      call getspecies(allatoms,species)
      nspecies=size(species)
      call write_coords(foutcoords,foutform,allatoms,nallatoms,
     &                  species,nspecies,vecs)

      if(talk) print fsubendext, 'multi'

      end subroutine


      end module
