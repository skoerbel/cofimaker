      module for_mbpp
      implicit none

      contains

c---------------------------------------------------------------------

      subroutine coorat2xsf()
      ! this subroutine transforms a COORAT file (MBPP input) to an xsf
      ! file (readable by VESTA). In INPUT.COFIMA the lattice vectors
      ! must be given.
      
      use defs
      implicit none
     
      integer i,j,natoms,natom
      character infilename*20,outfilename*20,units*4,line*100,nameat*2
      double precision xfrac(1:3),xabs(1:3),vecs(1:3,1:3)

      type(atom) , allocatable :: allatoms(:)

      write(12,fsubstart)trim(adjustl('coorat2xsf'))
      read(10,'(A20)') infilename
      read(10,'(A20)') outfilename
      read(10,*) units
      read(10,*) vecs(1,1:3)
      read(10,*) vecs(2,1:3)
      read(10,*) vecs(3,1:3)
      select case(units)
      case('bohr')
            vecs=vecs*bohr
      case('angs')
      case default
            print*,'        Error: Unknown units for lattice '
            print*,'        vectors. Options are: bohr, angs.'
            write(12,*) '        Error: Unknown units for lattice '
            write(12,*) '        vectors. Options are: bohr, angs.'
            nerr=nerr+1
      end select
      write(12,'(8x,"I will read from ",A)')infilename
      write(12,'(8x,"and write to ",A,".")')trim(adjustl(outfilename))
      open(21,file=infilename,status='old')
      open(22,file=outfilename,status='replace')
      write(22,*) 'CRYSTAL'
      write(22,*) 'PRIMVEC'
      do i=1,3
        write(22,'(3(F20.10))') vecs(i,1:3)
      end do
      write(22,*) 'PRIMCOORD'
      ! count atoms
      natoms=0
10    read(21,'(A100)',end=15) line
      line=adjustl(line)
      if (line(1:1).ne.'!') then
            if(index(line,'natom').le.0)natoms=natoms+1
      end if
      goto 10
15    continue
      rewind(21)
      write(22,*) natoms,1
20    read(21,'(A100)',end=25) line
      !get number of atoms of each species
      line=adjustl(line)
      if (line(1:1).ne.'!') then
        if (index(line,'natom').gt.0) then
            i=index(line,'natom')+5
            do while (line(i:i).eq.'=') 
                  i=i+1
            end do
            do while (line(i:i).eq.' ') 
                  i=i+1
            end do
            j=i
            do while (line(j:j).ne.' ')
                  j=j+1
            end do
!            read(line(i:j),'(i)') natom
            read(line(i:j),*) natom
        end if
        !get name of each species
        if (index(line,'name').gt.0) then
            i=index(line,'name')+4
            do while (line(i:i).eq.'=') 
                  i=i+1
            end do
            do while (line(i:i).eq.' ') 
                  i=i+1
            end do
            nameat=line(i:i+1)
        end if
        i=1
        do while (i.le.natom)
          read(21,'(A100)',end=23) line
          if(index(line,'natom').ge.1) goto 23
          line=adjustl(line)
          if(line(1:1).ne.'!') then
            read(line,*) xfrac(1:3)
            i=i+1
            do j=1,3
              xabs(j)=xfrac(1)*vecs(1,j)+xfrac(2)*vecs(2,j)
     &               +xfrac(3)*vecs(3,j)
            end do  
            do j=1,100
              if (nameat.eq.elements(j)) 
     &          write(22,'(I5,3(F20.10))') j,xabs(1:3)
            end do
          end if
        end do
      end if
      goto 20

23    continue
      print*,'        Error: coordinates section ended' 
      write(12,*)'       Error: coordinates section ended' 
      write(12,*)'       unexpectedly.'
      nerr=nerr+1
25    continue
      close(21)
      close(22)
      write(12,fsubend)
      end subroutine

c---------------------------------------------------------------------
      
      subroutine MBPP_transform_coords()

      use misc
      use defs
      implicit none

      ! read old coordinates from COORAT.IN file and writes new coordinates COORAT.OUT file 

      double precision matrix(1:3,1:3),tol,dummie(1:3)
      double precision inmat(1:3,1:3),outmat(1:3,1:3),outimat(1:3,1:3)
      double precision shift(1:3),shift_final(1:3),
     &                 invecs(1:3,1:3),outvecs(1:3,1:3)
      double precision,allocatable :: x_new(:,:),x_old(:,:)
      integer n,m,k,l,o,p,q,r,s,natom,natom_final,nuc(1:3),nsc(1:3)
      character line*70,line1*70,line_final*70,line2*100
      character transformation*20,infile*20,outfile*20
      logical,allocatable :: delete(:)
      logical domirror(1:3)
      parameter(tol=0.05D0)
      !character FMT1*1024,FMT2*1024,FMT3*1024,FMT4*1024,FMT5*1024
      
      write(12,fsubstart)trim(adjustl('MBPP_transform_coords'))
      
      ! read names of input and output file
      read(10,'(A100)',end=150) line2
      line2=adjustl(line2)
      if (line2(1:1).ne.'#') then
            read(line2,*) infile
      end if
      read(10,'(A100)',end=150) line2
      line2=adjustl(line2)
      if (line2(1:1).ne.'#') then
            read(line2,*) outfile
      end if
      write(12,'(A25,A)')'        I will read from ',infile
      write(12,'(A25,A)')'        and write to ',outfile

      ! choose transformation. Predefined are sc2fcc.
      ! For others specify 'userdef' and edit case 'userdef' below.
      !transformation='sc2fcc'   
      !transformation='fcc2sc'   
      read(10,'(A100)',end=150) line2
      line2=adjustl(line2)
      if (line2(1:1).ne.'#') then
            read(line2,*) transformation
      end if
      
      ! take also atoms of neighboring old unit cells into account, if
      ! the new cell lies partly outside the old one
      nuc(1)=3
      nuc(2)=3
      nuc(3)=3
      
      ! transform old coordinates into new ones: x_new = M * x_old

      SELECT CASE (transformation)
            CASE ('identity')
                matrix=0.0D0
                do n=1,3
                  matrix(n,n)=1.0D0
                end do
            CASE ('sc2fcc')
                matrix(1,1)=-1.0D0
                matrix(1,2)=1.0D0
                matrix(1,3)=1.0D0
                matrix(2,1)=1.0D0
                matrix(2,2)=-1.0D0
                matrix(2,3)=1.0D0
                matrix(3,1)=1.0D0
                matrix(3,2)=1.0D0
                matrix(3,3)=-1.0D0
                matrix=matrix*0.5D0
            CASE ('fcc2sc')
                matrix(1,1)=0.0D0
                matrix(1,2)=1.0D0
                matrix(1,3)=1.0D0
                matrix(2,1)=1.0D0
                matrix(2,2)=0.0D0
                matrix(2,3)=1.0D0
                matrix(3,1)=1.0D0
                matrix(3,2)=1.0D0
                matrix(3,3)=0.0D0
                !matrix=0.5D0*matrix
            CASE ('sc2hex')
                matrix(1,1)=2.0D0
                matrix(1,2)=2.0D0
                matrix(1,3)=-4.0D0
                matrix(2,1)=-4.0D0
                matrix(2,2)=2.0D0
                matrix(2,3)=2.0D0
                matrix(3,1)=1.0D0
                matrix(3,2)=1.0D0
                matrix(3,3)=1.0D0
                matrix=matrix/3.0D0
            CASE ('matrix')
                do n=1,3
                  read(10,*,end=150) matrix(n,1:3)
                end do
            CASE ('vector')
                do n=1,3
                  read(10,*,end=150) invecs(n,1:3)
                end do
                do n=1,3
                  read(10,*,end=150) outvecs(n,1:3)
                end do
                do n=1,3
                  inmat(1:3,n)=invecs(n,1:3)
                  outmat(1:3,n)=outvecs(n,1:3)
                end do
                call invert_matrix(outmat,outimat)
                matrix=matmul(outimat,inmat)
            CASE DEFAULT
                write(12,'(8x,"Error: unknown transformation.")') 
                write(12,'(8x,"Implemented transformations are:")') 
                write(12,'(8x,"sc2fcc,sc2hex,fcc2sc,matrix,")') 
                write(12,'(8x,"identity,vector")') 
                print*, '        Error: unknown transformation' 
                print*, '        (subroutine MBPP_transform_coords).' 
                nerr=nerr+1
      END SELECT
      
      ! read desired supercell dimensions
      read(10,'(A100)',end=150) line2
      line2=adjustl(line2)
      if (line2(1:1).ne.'#') then
            read(line2,*) nsc(1:3)
      end if
      
!      write(12,'(A35,3(I))')'        Your desired supercell dimension:',
!     &                    nsc
        WRITE(FMT2,*) nsc(1)
        WRITE(FMT3,*) nsc(2)
        WRITE(FMT4,*) nsc(3)
        WRITE(FMT1,*) '(A35,I',len(trim(FMT2)),',I',len(trim(FMT2)),',I'
     &   ,len(trim(FMT2)),')'
        write(12,FMT1)  nsc 
      do n=1,3
        matrix(n,1:3)=matrix(n,1:3)/dble(nsc(n))
      end do  
      write(12,1111)matrix
 1111 format(//8x,'Your transformation matrix is this:',
     &       //8x,3(f15.5),8x,3(f15.5),8x,3(f15.5))

      write(12,*)'        You chose this transformation: ',
     & transformation

      ! first shift
      shift(1)=0.0D0
      shift(2)=0.0D0
      shift(3)=0.0D0
      read(10,'(A100)',end=150) line2
      line2=adjustl(line2)
      if (line2(1:1).ne.'#') then
            read(line2,*) shift(1:3)
      end if
      write(12,'(A30,3(F10.6))')'        You chose as 1. shift:',shift

      ! then mirror
      domirror(1)=.False.
      domirror(2)=.False.
      domirror(3)=.False.
      read(10,'(A100)',end=150) line2
      line2=adjustl(line2)
      if (line2(1:1).ne.'#') then
            read(line2,*) domirror(1:3)
      end if
      write(12,'(A30,3(I2))')'        You chose to mirror: ',domirror
      
      ! then shift again
      shift_final(1)=0.0D0
      shift_final(2)=0.0D0
      shift_final(3)=0.0D0
      read(10,'(A100)',end=150) line2
      line2=adjustl(line2)
      if (line2(1:1).ne.'#') then
            read(line2,*) shift_final(1:3)
      end if
      write(12,'(A30,3(F10.6))')'        You chose as 2. shift:'
     &  ,shift_final

      open(40,file=infile,status='old')
      open(41,file=outfile,status='replace')

      ! read natom
 100  read(40,1000,end=125) line
      if (index(line,'natom').ge.1) then
            k=1
            do while (line(k:k+4).ne.'natom')
                  k=k+1
            end do
            do while (line(k:k).ne.'=')
                  k=k+1
            end do
            k=k+1
            do while (line(k:k).eq.' ')
                  k=k+1
            end do
            m=k
            do while (line(m:m).ne.' ')
                  m=m+1
            end do
            !open(12,file='temp',status='replace')
            !write(12,*) line(k:m)
            !rewind(12)
            !read(12,*) natom
            !close(12,status='delete')
            read(line(k:m),*) natom
            allocate(x_old(1:natom*nuc(1)*nuc(2)*nuc(3),1:3))
            allocate(x_new(1:natom*nuc(1)*nuc(2)*nuc(3),1:3))
            allocate(delete(1:natom*nuc(1)*nuc(2)*nuc(3)))
            ! read old coordinates for one species
            p=0
            do n=1,natom
            !read(40,*,end=200) x_old(n,1),x_old(n,2),x_old(n,3)
            read(40,'(A70)',end=200) line1
            line1=adjustl(line1)
            do while (line1(1:1).eq.'!') 
                  read(40,'(A70)',end=200) line1
                  line1=adjustl(line1)
            end do
            read(line1(1:70),*,end=200) x_old(n,1),x_old(n,2),x_old(n,3)
                  do q=0,nuc(1)-1
                    do r=0,nuc(2)-1
                      do s=0,nuc(3)-1
                        p=p+1
                        ! write old coordinates in units of new lattice
                        ! vectors
                        x_new(p,1)=(x_old(n,1)+float(q))*matrix(1,1)
     &                            +(x_old(n,2)+float(r))*matrix(1,2)
     &                            +(x_old(n,3)+float(s))*matrix(1,3)
                        x_new(p,2)=(x_old(n,1)+float(q))*matrix(2,1)
     &                            +(x_old(n,2)+float(r))*matrix(2,2)
     &                            +(x_old(n,3)+float(s))*matrix(2,3)
                        x_new(p,3)=(x_old(n,1)+float(q))*matrix(3,1)
     &                            +(x_old(n,2)+float(r))*matrix(3,2)
     &                            +(x_old(n,3)+float(s))*matrix(3,3)
                      end do  
                    end do  
                  end do
            end do
            ! shift new coordinates
            do n=1,natom*nuc(1)*nuc(2)*nuc(3)
                  do p=1,3
                        x_new(n,p)=x_new(n,p)+shift(p)
                  end do
            end do
            ! if new coordinates not between 0 and 1, shift them
            do p=1,natom*nuc(1)*nuc(2)*nuc(3)
              do o=1,3
                  do while (x_new(p,o).ge.1.0D0-tol) 
                        x_new(p,o)=x_new(p,o)-1.0D0
                  end do
                  do while (x_new(p,o).lt.0.0D0-tol) 
                        x_new(p,o)=x_new(p,o)+1.0D0
                  end do
              end do
            end do
            ! mirror new coordinates
            do p=1,3
                  if(domirror(p).eqv..True.) then
                        do n=1,natom*nuc(1)*nuc(2)*nuc(3)
                          x_new(n,p)=1.0D0-x_new(n,p)
                          if (x_new(n,p).ge.1.0D0-tol) then
                              x_new(n,p)=x_new(n,p)-1.0D0
                          end if
                        end do
                  end if
            end do
            ! shift new coordinates again
            do n=1,natom*nuc(1)*nuc(2)*nuc(3)
                  do p=1,3
                        x_new(n,p)=x_new(n,p)+shift_final(p)
                  end do
            end do
            ! if new coordinates not between 0 and 1, shift them
            do p=1,natom*nuc(1)*nuc(2)*nuc(3)
              do o=1,3
                  do while (x_new(p,o).ge.1.0D0-tol) 
                        x_new(p,o)=x_new(p,o)-1.0D0
                  end do
                  do while (x_new(p,o).lt.0.0D0-tol) 
                        x_new(p,o)=x_new(p,o)+1.0D0
                  end do
              end do
            end do
            ! sort new coordinates
            do n=1,natom*nuc(1)*nuc(2)*nuc(3)-1
                  do p=n+1,natom*nuc(1)*nuc(2)*nuc(3)
                        if(x_new(p,1).lt.x_new(n,1)) then
                              dummie(1:3)=x_new(p,1:3)
                              x_new(p,1:3)=x_new(n,1:3)
                              x_new(n,1:3)=dummie(1:3)
                        end if
                  end do
            end do
            do n=1,natom*nuc(1)*nuc(2)*nuc(3)-1
                  do p=n+1,natom*nuc(1)*nuc(2)*nuc(3)
                        if(x_new(p,2).lt.x_new(n,2).and.
     &                  abs(x_new(p,1)-x_new(n,1)).le.tol) then
                              dummie(1:3)=x_new(p,1:3)
                              x_new(p,1:3)=x_new(n,1:3)
                              x_new(n,1:3)=dummie(1:3)
                        end if
                  end do
            end do
            do n=1,natom*nuc(1)*nuc(2)*nuc(3)-1
                  do p=n+1,natom*nuc(1)*nuc(2)*nuc(3)
                        if(x_new(p,3).lt.x_new(n,3).and.
     &                  abs(x_new(p,1)-x_new(n,1)).le.tol.and.
     &                  abs(x_new(p,2)-x_new(n,2)).le.tol) then
                              dummie(1:3)=x_new(p,1:3)
                              x_new(p,1:3)=x_new(n,1:3)
                              x_new(n,1:3)=dummie(1:3)
                        end if
                  end do
            end do
            ! delete doubled coordinates
            delete=.False.
            do n=1,natom*nuc(1)*nuc(2)*nuc(3)-1
                  do p=n+1,natom*nuc(1)*nuc(2)*nuc(3)
                        if(abs(x_new(p,1)-x_new(n,1)).le.tol.and.
     &                     abs(x_new(p,2)-x_new(n,2)).le.tol.and.
     &                     abs(x_new(p,3)-x_new(n,3)).le.tol) then
                              delete(p)=.True.
                        end if
                  end do
            end do
            ! find final number of coordinates
            natom_final=0
            do n=1,natom*nuc(1)*nuc(2)*nuc(3)
                  if(delete(n).eqv..False.) natom_final=natom_final+1
            end do
            ! write to file
            line_final=line
!            write(41,'(A,1x,I5,3x,A)') line_final(1:k-2),natom_final,
!     &       line_final(k+3:70)
            write(41,'(A,1x,I5,3x,A)') line_final(1:k-1),natom_final,
     &       line_final(m+1:70)
            do n=1,natom*nuc(1)*nuc(2)*nuc(3)
                  if (delete(n).eqv..False.) then
                        write(41,2000) x_new(n,1),x_new(n,2),x_new(n,3)
                  end if
            end do
            deallocate(x_old,x_new,delete)
      end if
      goto 100

! 1000 format(A70)     
 1000 format(A)     
 2000 format(3(F15.10))     

125   close(40)
      close(41)
      
      write(12,fsubend)
      return

150   write(12,*)'      Error: File INPUT.COFIMA ended unexpectedly.'  
      close(40)
      close(41)
      nerr=nerr+1
      return
200   write(12,*)'      Error: File COORAT.IN ended unexpectedly.'  
      close(40)
      close(41)
      nerr=nerr+1
      return

      end subroutine

c---------------------------------------------------------------------

      subroutine coorat2gin()
      ! this subroutine transforms a COORAT file (MBPP input) to a gin
      ! file (GULP input file). 
      
      use defs
      implicit none
     
      integer i,j,natoms,natom
      character infilename*20,outfilename*20,units*4,line*100,nameat*2
      double precision xfrac(1:3)

      type(atom) , allocatable :: allatoms(:)

      write(12,fsubstart)trim(adjustl('coorat2gin'))
      read(10,'(A20)') infilename
      read(10,'(A20)') outfilename
      write(12,'(8x,"I will read from ",A)')infilename
      write(12,'(8x,"and write to ",A,".")')trim(adjustl(outfilename))
      open(21,file=infilename,status='old')
      open(22,file=outfilename,status='replace')
      write(22,*) 'frac'
      ! count atoms
      natoms=0
10    read(21,'(A100)',end=15) line
      line=adjustl(line)
      if (line(1:1).ne.'!') then
            if(index(line,'natom').le.0)natoms=natoms+1
      end if
      goto 10
15    continue
      rewind(21)
20    read(21,'(A100)',end=25) line
      !get number of atoms of each species
      line=adjustl(line)
      if (line(1:1).ne.'!') then
        if (index(line,'natom').gt.0) then
            i=index(line,'natom')+5
            do while (line(i:i).eq.'=') 
                  i=i+1
            end do
            do while (line(i:i).eq.' ') 
                  i=i+1
            end do
            j=i
            do while (line(j:j).ne.' ')
                  j=j+1
            end do
            read(line(i:j),*) natom
!            read(line(i:j),'(i)') natom
        end if
        !get name of each species
        if (index(line,'name').gt.0) then
            i=index(line,'name')+4
            do while (line(i:i).eq.'=') 
                  i=i+1
            end do
            do while (line(i:i).eq.' ') 
                  i=i+1
            end do
            nameat=line(i:i+1)
        end if
        ! write coordinates
        i=1
        do while (i.le.natom)
          read(21,'(A100)',end=23) line
          if(index(line,'natom').ge.1) goto 23
          line=adjustl(line)
          if(line(1:1).ne.'!') then
            read(line,*) xfrac(1:3)
            i=i+1
            do j=1,100
              if (nameat.eq.elements(j)) then
               write(22,'(A4,4x," core ",3(F15.10),A18)') 
     &               nameat,xfrac(1:3),' 1.0 1.0 0.0 1 1 1'
               write(22,'(A4,4x," shel ",3(F15.10),A18)') 
     &               nameat,xfrac(1:3),' 1.0 1.0 0.0 1 1 1'
              end if
            end do
          end if
        end do
      end if
      goto 20

23    continue
      print*,'        Error: coordinates section ended' 
      write(12,*)'       Error: coordinates section ended' 
      write(12,*)'       unexpectedly.'
      nerr=nerr+1
25    continue
      close(21)
      close(22)
      write(12,fsubend)
      end subroutine

c---------------------------------------------------------------------

      end module
