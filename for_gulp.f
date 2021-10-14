      module for_gulp
      implicit none

      contains

c---------------------------------------------------------------------

      subroutine trg2data()
      ! reads MD data from GULP trajectory output file
      
      use misc
      use defs
      implicit none
      
      integer i,j,natoms,nmirror,ndata
      double precision time,ekin,etot,temp,cell(1:3,1:3),alat(1:3)
      double precision angles(1:3),pol(1:3),polav(1:3),vol,totalcharge
      double precision  mirrorcoords(1:8,1:3),tolbound(1:3)
      character infilename*20,atomfilename*20,outfilename*20,line*100
      character ifpolstring*20,char1*20,char2*20
      logical ifpol
      type(atom), allocatable :: allatoms(:)
      logical isopen12
      !character FMT1*1024,FMT2*1024

      inquire(unit=12,opened=isopen12)  

      if(isopen12) then  
        write(12,fsubstart)trim(adjustl('trg2data'))
      else
        print fsubstart,trim(adjustl('trg2data'))
      end if
!      write(12,'(8x,"Reading input filename...")')
      read(10,'(A20)') infilename
      if(isopen12) then  
        write(12,'(8x,"Input file is ",A)')trim(infilename)
      else
        print '(8x,"Input file is ",A)',trim(infilename)
      end if
!      write(12,'(8x,"Reading res filename...")')
      read(10,'(A20)') atomfilename
      if(isopen12) then  
        write(12,'(8x,"Input file is ",A)')trim(atomfilename)
!        write(12,'(8x,"Reading output filename...")')
      else
        print '(8x,"Input file is ",A)',trim(atomfilename)
!        print '(8x,"Reading output filename...")'
      end if
      read(10,'(A20)') outfilename
!      write(12,'(8x,"Output file is ",A)')trim(outfilename)
!      write(12,'(8x,"Reading if pol. is to be calc....")')
      read(10,'(A20)') ifpolstring
      if(isopen12) then  
        write(12,'(8x," ",A)')trim(ifpolstring)
      else
        print '(8x," ",A)',trim(ifpolstring)
      end if
      tolbound=tolfrac      
7     read(10,'(A100)',end=8, err=100) line
      line=adjustl(line)
      if (line(1:1).ne.'#'.and.index(line,'tol').gt.0.and.
     &    index(line,'bound').gt.0) then
            read(line,*,err=100) char1,char2, tolbound(1:3)
      end if
8     continue
      if (index(ifpolstring,'pol').ge.1.and.index(ifpolstring,'yes')
     &    .ge.1) then
            ifpol=.true.
      else
            ifpol=.false.
      end if
      if (isopen12) then
        write(12,*)'       '
      else
        print*,'       '
      end if          
!      write(12,*)'       I will read MD data from file ',
!     &                   infilename
!      write(12,*)'       I will read the atom names from file ',
!     &                   atomfilename
!      write(12,*)'       I will write to file ',outfilename
!      if(ifpol) write(12,*)'       I will calculate polarizations.'
!      if(.not.ifpol)write(12,*)"       I won't calculate polarizations."

      open(23,file=infilename, status='old')
      rewind(23)
      read(23,'(A100)') line
      read(23,*) natoms
!      write(12,'(8x,"Found ",I," atoms.")') natoms
      WRITE(FMT2,*) natoms
      WRITE(FMT1,*) '(8x,"Found ",I',len(trim(FMT2)),
     &  ',"atoms")'
      if(isopen12) then
        write(12,FMT1) natoms
      else
        print FMT1, natoms
      end if
      allocate(allatoms(natoms))
      open(24,file=atomfilename,status='old')
      rewind(24)
      ! read order of atoms
      allatoms%name='        '
10    read(24,'(A100)',end=15) line
      line=adjustl(line)
      if (index(line,'frac').ge.1.or.index(line,'cart').ge.1) then
            i=1
            totalcharge=0.0D0
            do while (i.le.natoms)
                  read(24,'(A100)') line
                  line=adjustl(line)
                  if (line(1:1).ne.'#') then
                        read(line,*) allatoms(i)%name
                        read(line(len_trim(allatoms(i)%name)+1:100),*)
     &                  allatoms(i)%core,allatoms(i)%where(1:3),
     &                  allatoms(i)%charge  
                        totalcharge=totalcharge+allatoms(i)%charge
                        i=i+1
                  end if
            end do
            if(isopen12) then
             write(12,'(A23,F10.5)')'         Total charge:',totalcharge
            else
             print '(A23,F10.5)','         Total charge:',totalcharge
            end if
            goto 15
      end if
      goto 10

15    continue
      close(24)
      open(25,file=outfilename,status='replace')
      if(ifpol) then
            write(25,'(14(A15))') '# t (ps)', 'ekin (eV)','etot (eV)',
     &'T (K)',' a ','  b ', ' c (Angs)', 'alpha', 'beta', 'gamma (°)', 
     & 'vol (Angs³)','P_x', 'P_y', 'P_z (uC/cm²)'       
        write(12,'(8x,"Boundary tolerance: ",3(F8.4))') tolbound
      else
            write(25,'(11(A15))') '# t (ps)', 'ekin (eV)','etot (eV)',
     &'T (K)',' a ','  b ', ' c (Angs)', 'alpha', 'beta', 'gamma (°)',
     &'vol (Angs³)'       
      end if
      ndata=0
      polav=0.0D0
20    read(23,'(A100)',end=25) line
      ndata=ndata+1
      read(23,*) time, ekin, etot, temp
      read(23,'(A100)') line  ! now come the coordinates
      do i=1,natoms
            read(23,*) allatoms(i)%where(1:3)
      end do
      read(23,'(A100)') line   ! now come the velocities
      do i=1,natoms
        read(23,'(A100)') line
      end do
      read(23,'(A100)') line   
      ! if gulp version 4, now come the derivatives 
      if (index(line,'Derivatives').gt.0) then
        do i=1,natoms
          read(23,'(A100)') line
        end do
        read(23,'(A100)') line   
      end if  
      ! if gulp version 4, now come the site energies 
      if (index(line,'energies').gt.0) then
        do i=1,natoms
          read(23,'(A100)') line
        end do
        read(23,'(A100)') line   
      end if  
      ! now come the cell parameters
      do i=1,3
            read(23,*) cell(i,1:3)
      end do
      ! calculate the cell volume from the cell vectors
      vol=cell(1,1)*(cell(2,2)*cell(3,3)
     &          -cell(2,3)*cell(3,2))
     &          +cell(1,2)*(cell(2,3)*cell(3,1)
     &          -cell(2,1)*cell(3,3))
     &          +cell(1,3)*(cell(2,1)*cell(3,2)
     &          -cell(2,2)*cell(3,1))
      vol=abs(vol)

      ! calculate the cell edges and angles from the cell vectors
      do i=1,3
        alat(i)=sqrt(cell(i,1)**2+cell(i,2)**2+cell(i,3)**2)
      end do
      angles(1)=acos((cell(2,1)*cell(3,1)+cell(2,2)*cell(3,2)
     &       +cell(2,3)*cell(3,3))/(alat(2)*alat(3)))*180.0D0/Pi
      angles(2)=acos((cell(1,1)*cell(3,1)+cell(1,2)*cell(3,2)
     &       +cell(1,3)*cell(3,3))/(alat(1)*alat(3)))*180.0D0/Pi
      angles(3)=acos((cell(1,1)*cell(2,1)+cell(1,2)*cell(2,2)
     &       +cell(1,3)*cell(2,3))/(alat(1)*alat(2)))*180.0D0/Pi
      read(23,'(A100)') line   ! now comes the strain
      read(23,'(A100)') line   
      read(23,'(A100)') line   
      ! calculate the polarization
      if (ifpol) then
            pol=0.0D0
            !the following simply adds over coordinates times charges.
            !do i=1,natoms
            !  pol=pol+allatoms(i)%charge*allatoms(i)%where
            !end do
            !
            !the following takes care of atoms that sit on the cell
            !boundaries. 
            if(ndata.eq.1) then
                if(isopen12) then
                  write(12,'(8x,
     &             "atom multipl. & coords of first structure:")') 
                else
                  print '(8x,
     &             "atom multipl. & coords of first structure:")' 
                end if
            end if
            do i=1,natoms
              call mirror(allatoms(i)%where,cell,nmirror,mirrorcoords,
     &               tolbound)
              if(ndata.eq.1) then
                if(isopen12) then
!                  write(12,'(8x,A3,I)') allatoms(i)%name,nmirror
                  WRITE(FMT2,*) nmirror
                  WRITE(FMT1,*) '(8x,A3,I',len(trim(FMT2)),')'
                  write (12,FMT1)nmirror
                else 
                  WRITE(FMT2,*) nmirror
                  WRITE(FMT1,*) '(8x,A3,I',len(trim(FMT2)),')'
                  print FMT1,nmirror
                end if
              end if
              do j=1,nmirror
              if(ndata.eq.1) then
                if(isopen12) then
                   write(12,'(8x,3F12.5)') mirrorcoords(j,1:3)
                else
                   print '(8x,3F12.5)', mirrorcoords(j,1:3)
                end if
              end if
              pol(1:3)=pol(1:3)+mirrorcoords(j,1:3)
     &                   *allatoms(i)%charge/dble(nmirror)
              end do  
            end do
            !
            ! crudely correct for net charge
            do i=1,3
              pol=pol-totalcharge*0.5D0*
     &          (cell(1,i)+cell(2,i)+cell(3,i))
            end do
            if(isopen12) then
              write(12,*)'        Total dipole moment (e*Angs):',pol
            else
              print*,'        Total dipole moment (e*Angs):',pol
            end if 
            pol=pol*polfac_Angs_e/vol
            polav=polav+pol
            write(25,'(14(F15.6))') time,ekin,etot,temp,
     &           (alat(i),i=1,3),(angles(i),i=1,3),vol,pol(1:3)
      else
            write(25,'(11(F15.6))') time,ekin,etot,temp,
     &           (alat(i),i=1,3),(angles(i),i=1,3),vol
      end if
      goto 20

25    continue
      close(23)
      if (ifpol) then
            ! averqged pol. = sum(pol vaues)/number of pol values
            polav=polav/dble(ndata)
            ! print averaged pol. to file 
            write(25,'("# time-av. polariz. components in uC/cm^2:")')
            write(25,'("# ",3(F10.4))') polav(1:3)
      end if
      close(25)
      if(isopen12) then
        write(12,fsubend)
      else
        print fsubend
      end if
      return

100   print*,"Error: something wrong with INPUT.COFIMA."      
      if(isopen12) then
        write(12,'(8x,"Error: something wrong with INPUT.COFIMA.")')
      else
        print'(8x,"Error: something wrong with INPUT.COFIMA.")'
      end if
      nerr=nerr+1
      return
      end subroutine
              
c---------------------------------------------------------------------

      subroutine gulp_md_data()
      ! extracts temperature and pressure from gulp MD output file 

      use defs
      implicit none
      
      double precision eq_time
      double precision time,ekin,epot,etot,temp,tcs,press,alat(1:3),
     &            angle(1:3),vol,ekinav,epotav,etotav,tav,tcsav,pav,
     &            alatav(1:3),angleav(1:3),volav
      integer np,nt
      character line*150,time_char*13,infilename*20,outfilename*20
      character outfilename2*30
      
      write(12,fsubstart)trim(adjustl('gulp_md_data'))

      read(10,*) infilename
      read(10,*) outfilename
      outfilename2=trim(adjustl(outfilename))
      outfilename2(len_trim(adjustl(outfilename))+1:30)='.AV'
      open(26,file=infilename,status='old')
      open(27,file=outfilename,status='replace')
      open(28,file=outfilename2,status='replace')
      
      eq_time=0.0D0
      write(27,1001) 't (ps)','ekin','epot','etot','T','T c/s',
     &      'p (GPa)',
     &      ' a','b','c','alpha','beta','gamma','V'
      write(28,1001) 't (ps)','ekin av',
     &      'epotav','etotav','Tav','Tc/sav','pav','aav',
     &      'bav','cav','alphaav','betaav','gammaav','volav'
10    read(26,'(A150)',end=20) line
      if (index(line,'  ** Time : ').ge.1) then
            read(line(12:150),*) time
      end if
            
      if (index(line,'       Kinetic energy    (eV) =').ge.1) then
            read(line(32:64),*) ekin,ekinav
      end if
      
      if (index(line,'       Potential energy  (eV) =').ge.1) then
            read(line(32:64),*) epot,epotav
      end if
      
      if (index(line,'       Total energy      (eV) =').ge.1) then
            read(line(32:64),*) etot,etotav
      end if
      
      if (index(line,'       Temperature       (K)  =').ge.1) then
            read(line(32:64),*) temp,tav
      end if
      
      if (index(line,'       C/S temperature   (K)  =').ge.1) then
            read(line(32:64),*) tcs,tcsav
      end if
      
      if (index(line,'       Pressure         (GPa) =').ge.1) then
            read(line(32:64),*) press,pav
      end if

      if (index(line,'       Cell parameter : a (A) =').ge.1) then
            read(line(32:64),*) alat(1),alatav(1)
      end if

      if (index(line,'       Cell parameter : b (A) =').ge.1) then
            read(line(32:64),*) alat(2),alatav(2)
      end if

      if (index(line,'       Cell parameter : c (A) =').ge.1) then
            read(line(32:64),*) alat(3),alatav(3)
      end if

      if (index(line,'       Cell angle : alpha (o) =').ge.1) then
            read(line(32:64),*) angle(1),angleav(1)
      end if

      if (index(line,'       Cell angle : beta  (o) =').ge.1) then
            read(line(32:64),*) angle(2),angleav(2)
      end if

      if (index(line,'       Cell angle : gamma (o) =').ge.1) then
            read(line(32:64),*) angle(3),angleav(3)
      end if

      if (index(line,'       Cell volume :   (A**3) =').ge.1) then
            read(line(33:65),*) vol,volav
            write(27,1000) time,ekin,epot,etot,temp,tcs,press,alat,
     &            angle,vol
            write(28,1000) time,ekinav,epotav,etotav,tav,tcsav,pav,
     &                    alatav,angleav,volav
      end if
      
      if (index(line,'  Molecular dynamics production :').ge.1) then
            eq_time=time
            goto 20
      end if
      goto 10

20    continue

30    read(26,'(A150)',end=40) line
      if (index(line,'  ** Time : ').ge.1) then
            read(line(12:150),*) time
            time=time+eq_time
      end if
            
      if (index(line,'       Kinetic energy    (eV) =').ge.1) then
            read(line(32:64),*) ekin,ekinav
      end if
      
      if (index(line,'       Potential energy  (eV) =').ge.1) then
            read(line(32:64),*) epot,epotav
      end if
      
      if (index(line,'       Total energy      (eV) =').ge.1) then
            read(line(32:64),*) etot,etotav
      end if
      
      if (index(line,'       Temperature       (K)  =').ge.1) then
            read(line(32:64),*) temp,tav
      end if
      
      if (index(line,'       C/S temperature   (K)  =').ge.1) then
            read(line(32:64),*) tcs,tcsav
      end if
      
      if (index(line,'       Pressure         (GPa) =').ge.1) then
            read(line(32:64),*) press,pav
      end if

      if (index(line,'       Cell parameter : a (A) =').ge.1) then
            read(line(32:64),*) alat(1),alatav(1)
      end if

      if (index(line,'       Cell parameter : b (A) =').ge.1) then
            read(line(32:64),*) alat(2),alatav(2)
      end if

      if (index(line,'       Cell parameter : c (A) =').ge.1) then
            read(line(32:64),*) alat(3),alatav(3)
      end if

      if (index(line,'       Cell angle : alpha (o) =').ge.1) then
            read(line(32:64),*) angle(1),angleav(1)
      end if

      if (index(line,'       Cell angle : beta  (o) =').ge.1) then
            read(line(32:64),*) angle(2),angleav(2)
      end if

      if (index(line,'       Cell angle : gamma (o) =').ge.1) then
            read(line(32:64),*) angle(3),angleav(3)
      end if

      if (index(line,'       Cell volume :   (A**3) =').ge.1) then
            read(line(33:65),*) vol,volav
            write(27,1000) time,ekin,epot,etot,temp,tcs,press,alat,
     &            angle,vol
            write(28,1000) time,ekinav,epotav,etotav,tav,tcsav,pav,
     &            alatav,angleav,volav
      end if
      
      !if (index(line,'  ** Time : ').ge.1) then
      !  read(line(12:24),*) time
      !  time=time+eq_time
      !end if
      !if (index(line,'Pressure  ').ge.1) then
      !      read(line(33:65),*) press,pav
      !end if
      !if (index(line,'Temperature   ').ge.1) then
      !      read(line(33:65),*) temp,tav
      !      write(27,1000) time,press,pav,temp,tav
      !end if
      goto 30


40    continue
      close(26)
      close(27)
      close(28)

! 1000 format(14(F10.5))     
 1000 format(F8.5,F12.5,2(F15.5),F12.5,8(F10.5),F15.5)     
! 1001 format(14(A10))
 1001 format(A8,A12,2(A15),A12,8(A10),A15)     
      
      write(12,fsubend)
      end subroutine

c---------------------------------------------------------------------

      subroutine trg2coorat()
      ! reads GULP trajectory and writes MBPP COORAT files 
      use misc
      use defs
      implicit none
      
      integer i,j,natoms,nmirror,inc,iinc,ncores,nshels,icore,ishel
      integer iinclen
      double precision time,cell(1:3,1:3),alat(1:3),angles(1:3)
      character infilename*20,atomfilename*20,line*100
      character outfilecores*30,outfileshels*30
      character char1*20,char2*20
      type(atom), allocatable :: allatoms(:),cores(:),shels(:)

      write(12,fsubstart)trim(adjustl('trg2coorat'))
      read(10,'(A20)') infilename
      read(10,'(A20)') atomfilename
      read(10,*) inc
      write(12,*)'       '
      write(12,*)'        I will read MD data from file ',
     &                   infilename
      write(12,*)'        I will read the atom names from file ',
     &                   atomfilename
      outfilecores='COORAT.cores.'
      outfileshels='COORAT.shels.'
      write(12,*)'        I will write to files ',outfilecores,
     &           outfileshels 
      write(12,*)'        every ',inc,' steps.'

      open(23,file=infilename, status='old')
      read(23,'(A100)') line
      read(23,*) natoms
      allocate(allatoms(natoms))
      open(24,file=atomfilename,status='old')
      ! read order of atoms and number of cores and shels
      allatoms%name='        '
      ncores=0
      nshels=0
10    read(24,'(A100)',end=15) line
      line=adjustl(line)
      if (index(line,'frac').ge.1.or.index(line,'cart').ge.1) then
            i=1
            do while (i.le.natoms)
                  read(24,'(A100)') line
                  line=adjustl(line)
                  if (line(1:1).ne.'#') then
                      read(line,*) allatoms(i)%name
                      read(line(len_trim(allatoms(i)%name)+1:100),*)
     &                allatoms(i)%core,allatoms(i)%where(1:3)
                      if (allatoms(i)%core.eq.'core') ncores=ncores+1
                      if (allatoms(i)%core.eq.'shel') nshels=nshels+1
                      i=i+1
                  end if
            end do
            goto 15
      end if
      goto 10

15    continue
      close(24)
      iinc=0
      allocate(cores(ncores),shels(nshels))
20    read(23,'(A100)',end=25) line
      !read(23,'(A100)') line
      read(23,*) time
      read(23,'(A100)') line  ! now come the coordinates
      do i=1,natoms
            read(23,*) allatoms(i)%where(1:3)
      end do
      ! get cores and shels
            icore=0
            ishel=0
      do i=1,natoms
            if (allatoms(i)%core.eq.'core') then
                  icore=icore+1
                  cores(icore)=allatoms(i)
                  cores(icore)%core='core'
            end if  
            if (allatoms(i)%core.eq.'shel') then
                  ishel=ishel+1
                  shels(ishel)=allatoms(i)
                  shels(ishel)%core='shel'
            end if  
      end do
      read(23,'(A100)') line   ! now come the velocities
      do i=1,natoms
        read(23,'(A100)') line
      end do
      read(23,'(A100)') line   
      ! if gulp version 4, now come the derivatives 
      if (index(line,'Derivatives').gt.0) then
        do i=1,natoms
          read(23,'(A100)') line
        end do
        read(23,'(A100)') line   
      end if  
      ! if gulp version 4, now come the site energies 
      if (index(line,'energies').gt.0) then
        do i=1,natoms
          read(23,'(A100)') line
        end do
        read(23,'(A100)') line   
      end if  
      ! now come the cell parameters
      do i=1,3
            read(23,*) cell(i,1:3)
      end do
      read(23,'(A100)') line   ! now comes the strain
      read(23,'(A100)') line   
      read(23,'(A100)') line   
      if(mod(iinc,inc).eq.0) then
!            write(char1,'(I)') iinc
            write(char1,*) iinc
            char1=trim(adjustl(char1))
            iinclen=len_trim(adjustl(char1))
            write(outfilecores(14:13+iinclen),'(A)') 
     &       adjustl(char1(1:iinclen)) 
            write(outfileshels(14:13+iinclen),'(A)') 
     &       adjustl(char1(1:iinclen)) 
            call writecoorat(cores,ncores,cell,outfilecores) !write(25,*) iinc
            call writecoorat(shels,nshels,cell,outfileshels) !write(25,*) iinc
            ! calculate the cell edges and angles from the cell vectors
            do i=1,3
              alat(i)=sqrt(cell(i,1)**2+cell(i,2)**2+cell(i,3)**2)/Bohr
            end do
            angles(1)=acos((cell(2,1)*cell(3,1)+cell(2,2)*cell(3,2)
     &             +cell(2,3)*cell(3,3))/(alat(2)*alat(3)))*180.0D0/Pi
            angles(2)=acos((cell(1,1)*cell(3,1)+cell(1,2)*cell(3,2)
     &             +cell(1,3)*cell(3,3))/(alat(1)*alat(3)))*180.0D0/Pi
            angles(3)=acos((cell(1,1)*cell(2,1)+cell(1,2)*cell(2,2)
     &             +cell(1,3)*cell(2,3))/(alat(1)*alat(2)))*180.0D0/Pi
            ! write COORAT with commented cell parameter section
            open(25,file=outfilecores,status='old',access='append')
            write(25,'("!time in ps: ",F15.7)')time
            write(25,'("!unitcell edges and angles (Bohr,degree):")')
            write(25,'("!",6(f15.8))') alat(1:3),angles(1:3)
            close(25)
            open(25,file=outfileshels,status='old',access='append')
            write(25,'("!time in ps: ",F15.7)')time
            write(25,'("!unitcell edges and angles (Bohr,degree):")')
            write(25,'("!",6(f15.8))') alat(1:3),angles(1:3)
            close(25)
      end if
      iinc=iinc+1
      goto 20

25    continue      
      close(23)
      !close(25)
      write(12,fsubend)


      end subroutine

c---------------------------------------------------------------------

      subroutine res2config()
      ! this subroutine transforms gulp restart file to config file
      ! format, but the coordinates remain fractional -> run again with
      ! "fractional to absolute" command to get CONFIG file with absolute positions
      
      use defs
      implicit none

      character infilename*20,outfilename*20,line*100,name*10
      character xchar*10,ychar*10,zchar*10
      double precision a,b,c,alpha,beta,gamm,vecs(1:3,1:3),x,y,z
      double precision dumrvec(1:3),counter,denominator
      integer imcon,natoms
      integer i,i1
      
      !type atom
      !sequence
      !character name*8
      !logical written
      !double precision :: where(3)
      !end type atom

      type(atom) , allocatable :: allatoms(:)

      !---------------------------------------------------------------

      write(12,fsubstart)trim(adjustl('res2config'))
      read(10,'(A20)')infilename
      read(10,'(A20)')outfilename
      read(10,*)imcon
      write(12,*) '      Inputfile: ',infilename
      write(12,*) '      Output file: ',outfilename
      write(12,*) '      Reading cell parameters...'
      open(17,file=infilename,status='old')
      open(18,file=outfilename,status='replace')
      do i=1,8
        read(17,'(A100)')line
      end do
      ! read cell parameters
      read(17,*) a,b,c,alpha,beta,gamm
      vecs=0.0D0
      vecs(1,1)=a
      vecs(2,1)=b*cos(gamm*pi/180.0D0)
      vecs(2,2)=b*sin(gamm*Pi/180.0D0)
      vecs(3,1)=c*cos(beta*Pi/180.0D0)
      vecs(3,2)=c*(cos(alpha*pi/180.0D0)
     &       -cos(gamm*pi/180.0D0)*cos(beta*pi/180.0D0))
     &       /sin(gamm*pi/180.0D0)
      vecs(3,3)=sqrt(c**2-vecs(3,1)**2-vecs(3,2)**2)
      read(17,'(A100)')line
      ! determine number of atoms
      natoms=0
10    read(17,'(A100)',end=15) line 
      if (index(line,'core').ge.1) natoms=natoms+1
      if (index(line,'shel').ge.1) natoms=natoms+1
      if (index(line,'totalenergy').ge.1) goto 15
      goto 10
15    continue
      ! output
      write(18,*)' '
      write(18,'(I10,I10,I16)') 0,imcon,natoms
      do i=1,3
        write(18,'(3(F20.10))') vecs(i,1:3)
      end do
      ! read coordinates
      allocate(allatoms(natoms))
      rewind(17)
      do i=1,10
        read(17,'(A100)') line
      end do
      do i=1,natoms
        read(17,'(4(A10))') name,xchar,ychar,zchar
        if(xchar(3:3).eq.'/') then
            read(xchar(1:2),'(F2.0)') counter
            read(xchar(4:10),'(F7.0)') denominator
            x=counter/denominator
        else
            read(xchar,'(F10.7)') x
        end if
        if(ychar(3:3).eq.'/') then
            read(ychar(1:2),'(F2.0)') counter
            read(ychar(4:10),'(F7.0)') denominator
            y=counter/denominator
        else
            read(ychar,'(F10.7)') y
        end if
        if(zchar(3:3).eq.'/') then
            read(zchar(1:2),'(F2.0)') counter
            read(zchar(4:10),'(F7.0)') denominator
            z=counter/denominator
        else
            read(zchar,'(F10.7)') z
        end if
        ! if close to 1, shift by -1
        if (x.gt.1.0D0-tolfrac(1)) x=x-1.0d0
        if (y.gt.1.0D0-tolfrac(2)) y=y-1.0d0
        if (z.gt.1.0D0-tolfrac(3)) z=z-1.0d0
        if (index('0123456789 ',name(2:2)).ge.1) name(2:2)='_'
        ! remove spaces from atomic names and write to file 
        name(3:6)=name(7:10)
        name(7:10)='    '
        allatoms(i)%name=name(1:8)
        allatoms(i)%where(1)=x
        allatoms(i)%where(2)=y
        allatoms(i)%where(3)=z
      end do  
      ! sort coordinates
      ! sort x coordinates
      do i=1,natoms
        do i1=i+1,natoms
          if (allatoms(i1)%name.eq.allatoms(i)%name.and.
     &     allatoms(i1)%where(1).lt.allatoms(i)%where(1)-tolfrac(1))
     &    then
            dumrvec=allatoms(i)%where
            allatoms(i)%where=allatoms(i1)%where
            allatoms(i1)%where=dumrvec
          end if
        end do
      end do
      ! sort y coordinates
      do i=1,natoms
        do i1=i+1,natoms
          if (allatoms(i1)%name.eq.allatoms(i)%name.and.
     &     allatoms(i1)%where(2).lt.allatoms(i)%where(2)-tolfrac(2)
     &     .and.abs(allatoms(i1)%where(1)-allatoms(i)%where(1))
     &     .lt.tolfrac(1)) then
            dumrvec=allatoms(i)%where
            allatoms(i)%where=allatoms(i1)%where
            allatoms(i1)%where=dumrvec
          end if
        end do
      end do
      ! sort z coordinates
      do i=1,natoms
        do i1=i+1,natoms
          if (allatoms(i1)%name.eq.allatoms(i)%name.and.
     &     allatoms(i1)%where(3).lt.allatoms(i)%where(3)-tolfrac(3)
     &     .and.abs(allatoms(i1)%where(1)-allatoms(i)%where(1))
     &     .lt.tolfrac(1)
     &     .and.abs(allatoms(i1)%where(2)-allatoms(i)%where(2))
     &     .lt.tolfrac(2)) then
            dumrvec=allatoms(i)%where
            allatoms(i)%where=allatoms(i1)%where
            allatoms(i1)%where=dumrvec
          end if
        end do
      end do
      ! transform to absolute coordinates and print to file
      do i=1,natoms
        x=allatoms(i)%where(1)
        y=allatoms(i)%where(2)
        z=allatoms(i)%where(3)
        allatoms(i)%where(1:3)=x*vecs(1,1:3)+y*vecs(2,1:3)
     &          +z*vecs(3,1:3)
        write(18,'(A8,I12)') allatoms(i)%name,i
        write(18,'(F16.9,2(F20.9))') allatoms(i)%where(1:3)
      end do  
           
      close(17)
      close(18)
      write(12,*) '     '
      
      write(12,fsubend)
      end subroutine

c---------------------------------------------------------------------

      subroutine res2coorat()
      ! converts between gulp restart file and MBPP COORAT file 

      use defs
      implicit none
      character infilename*20,outfilename*30,outfilenamec*30,
     & outfilenames*30,line*100,name*10
      character xchar*10,ychar*10,zchar*10,coreshel*10
      double precision a,b,c,alpha,beta,gamm,vecs(1:3,1:3),x,y,z
      double precision dumrvec(1:3),counter,denominator
      integer natoms
      integer i,i1
      type(atom) , allocatable :: allatoms(:)
      type(element), allocatable :: elesc(:),eless(:)
      integer nelec,neles,numelec,numeles,j,ielec,ieles
      logical newelec,neweles

      write(12,fsubstart)trim(adjustl('res2coorat'))
      read(10,'(A20)')infilename
      read(10,'(A20)')outfilename
      outfilename=adjustl(outfilename)
      outfilenamec=outfilename
      outfilenamec(len_trim(outfilename)+1:30)='.cores'
      outfilenames=outfilename
      outfilenames(len_trim(outfilename)+1:30)='.shels'
      write(12,'(8x,"Inputfile: ",A)')infilename
      write(12,'(8x,"Outputfile: ",A)') outfilename
      write(12,*) '      Reading cell parameters...'
      open(17,file=infilename,status='old')
      open(18,file=outfilenamec,status='replace')
      open(19,file=outfilenames,status='replace')
5     read(17,'(A100)',end=6) line
      line=adjustl(line)
      if (index(line,'cell').ge.1.and.line(1:1).ne.'#') goto 7
      goto 5
      ! Stop if "cell" is not found
6     print*,"Error: no cell parameters found."
      write(12,'(8x,"Error: no cell parameters found.")')
      nerr=nerr+1
      close(17)
      close(18)
      close(19)
      return
7     continue
      ! read cell parameters
      read(17,*) a,b,c,alpha,beta,gamm
      vecs=0.0D0
      vecs(1,1)=a
      vecs(2,1)=b*cos(gamm*pi/180.0D0)
      vecs(2,2)=b*sin(gamm*Pi/180.0D0)
      vecs(3,1)=c*cos(beta*Pi/180.0D0)
      vecs(3,2)=c*(cos(alpha*pi/180.0D0)
     &       -cos(gamm*pi/180.0D0)*cos(beta*pi/180.0D0))
     &       /sin(gamm*pi/180.0D0)
      vecs(3,3)=sqrt(c**2-vecs(3,1)**2-vecs(3,2)**2)
      ! determine number of atoms
      natoms=0
8     read(17,'(A100)',end=9) line
      line=adjustl(line)
      if(index(line,'frac').ge.1.and.line(1:1).ne.'#') goto 10
      goto 8
      write(12,'(8x,"Reading atom info...")')
      print*,"reading atom info ..."
      ! Stop if "frac" is not found
9     print*,"Error: no fractional coordinates found."
      write(12,'(8x,"Error: no fractional coordiantes found.")')
      nerr=nerr+1
      close(17)
      close(18)
      close(19)
      return
10    read(17,'(A100)',end=15) line 
      line=adjustl(line)
      if(line(1:1).ne.'#') then
        if (index(line,'core').ge.1) natoms=natoms+1
        if (index(line,'shel').ge.1) natoms=natoms+1
        if (index(line,'totalenergy').ge.1) goto 15
        if(index(line,'core').le.0.and.index(line,'shel').le.0)goto 15
      end if
      goto 10
15    continue
      ! output
      write(12,'(8x,"There are ",I5," atoms in the file.")')natoms
      do i=1,3
        write(18,'("! cell vecs:",3(F20.10))') vecs(i,1:3)
        write(19,'("! cell vecs:",3(F20.10))') vecs(i,1:3)
      end do
      ! read coordinates
      allocate(allatoms(natoms))
      rewind(17)
20    read(17,'(A100)') line  
      line=adjustl(line)
      if (index(line,'frac').ge.1.and.line(1:1).ne.'#')goto 30
      goto 20
30    continue 
      do i=1,natoms
        read(17,'(A100)') line
        line=adjustl(line)
        do while (line(1:1).eq.'#') 
          read(17,'(A100)')line
          line=adjustl(line)
        end do 
        do j=1,len(line)
          if(line(j:j).eq.'/') line(j:j)='_'
        end do
        read(line,*) name,coreshel, xchar,ychar,zchar
        if(index(xchar,'_').ge.1) then
            read(xchar(1:index(xchar,'_')-1),*) counter
            read(xchar(index(xchar,'_')+1:10),*) denominator
            x=counter/denominator
        else
            read(xchar,'(F10.7)') x
        end if
        if(index(ychar,'_').ge.1) then
            read(ychar(1:index(ychar,'_')-1),*) counter
            read(ychar(index(ychar,'_')+1:10),*) denominator
            y=counter/denominator
        else
            read(ychar,'(F10.7)') y
        end if
        if(index(zchar,'_').ge.1) then
            read(zchar(1:index(zchar,'_')-1),*) counter
            read(zchar(4:10),'(F7.0)') denominator
            read(zchar(index(zchar,'_')+1:10),*) denominator
            z=counter/denominator
        else
            read(zchar,'(F10.7)') z
        end if
        ! if close to 1, shift by -1
        if (x.gt.1.0D0-tolfrac(1)) x=x-1.0d0
        if (y.gt.1.0D0-tolfrac(2)) y=y-1.0d0
        if (z.gt.1.0D0-tolfrac(3)) z=z-1.0d0
        !if (index('0123456789 ',name(2:2)).ge.1) name(2:2)='_'
        ! remove spaces from atomic names and write to file 
        name(3:6)=name(7:10)
        name(7:10)='    '
        allatoms(i)%name=name(1:8)
        allatoms(i)%core=trim(adjustl(coreshel))
        allatoms(i)%where(1)=x
        allatoms(i)%where(2)=y
        allatoms(i)%where(3)=z
      end do  
!      ! sort coordinates
!      ! sort x coordinates
!      do i=1,natoms
!        do i1=i+1,natoms
!          if (allatoms(i1)%name.eq.allatoms(i)%name.and.
!     &     allatoms(i1)%where(1).lt.allatoms(i)%where(1)-tolfrac(1))
!     &    then
!            dumrvec=allatoms(i)%where
!            allatoms(i)%where=allatoms(i1)%where
!            allatoms(i1)%where=dumrvec
!          end if
!        end do
!      end do
!      ! sort y coordinates
!      do i=1,natoms
!        do i1=i+1,natoms
!          if (allatoms(i1)%name.eq.allatoms(i)%name.and.
!     &     allatoms(i1)%where(2).lt.allatoms(i)%where(2)-tolfrac(2)
!     &     .and.abs(allatoms(i1)%where(1)-allatoms(i)%where(1))
!     &     .lt.tolfrac(1)) then
!            dumrvec=allatoms(i)%where
!            allatoms(i)%where=allatoms(i1)%where
!            allatoms(i1)%where=dumrvec
!          end if
!        end do
!      end do
!      ! sort z coordinates
!      do i=1,natoms
!        do i1=i+1,natoms
!          if (allatoms(i1)%name.eq.allatoms(i)%name.and.
!     &     allatoms(i1)%where(3).lt.allatoms(i)%where(3)-tolfrac(3)
!     &     .and.abs(allatoms(i1)%where(1)-allatoms(i)%where(1))
!     &     .lt.tolfrac(1)
!     &     .and.abs(allatoms(i1)%where(2)-allatoms(i)%where(2))
!     &     .lt.tolfrac(2)) then
!            dumrvec=allatoms(i)%where
!            allatoms(i)%where=allatoms(i1)%where
!            allatoms(i1)%where=dumrvec
!          end if
!        end do
!      end do
      !! transform to absolute coordinates 
      !do i=1,natoms
        !x=allatoms(i)%where(1)
        !y=allatoms(i)%where(2)
        !z=allatoms(i)%where(3)
        !allatoms(i)%where(1:3)=x*vecs(1,1:3)+y*vecs(2,1:3)
!     &  !        +z*vecs(3,1:3)
        !write(18,'(A8,I12)') allatoms(i)%name,i
        !write(18,'(F16.9,2(F20.9))') allatoms(i)%where(1:3)
      !end do  
      ! count elements and number of atoms of each element
      allatoms%written=.false.
      ! cores and shels
      nelec=0
      neles=0
      do i=1,natoms
            newelec=.true.
            neweles=.true.
            if (allatoms(i)%core(1:4).eq.'shel') newelec=.false.
            if (allatoms(i)%core(1:4).eq.'core') neweles=.false.
            do j=1,i-1
              if (allatoms(i)%name(1:8).eq.allatoms(j)%name(1:8).and.
     &           allatoms(i)%core(1:4).eq.'core'.and.
     &           allatoms(j)%core(1:4).eq.'core') newelec=.false.
              if (allatoms(i)%name(1:8).eq.allatoms(j)%name(1:8).and.
     &           allatoms(i)%core(1:4).eq.'shel'.and.
     &           allatoms(j)%core(1:4).eq.'shel') neweles=.false.
            end do
            if (newelec) nelec=nelec+1
            if (neweles) neles=neles+1
      end do
      !write(12,'(8x,"Found ",I4," core elements.")')nelec
      allocate(elesc(1:nelec))
!      ! shels
!      neles=0
!      do i=1,natoms
!            newele=.true.
!            if (allatoms(i)%core(1:4).eq.'core') newele=.false.
!            do j=1,i-1
!              if (allatoms(i)%name(1:8).eq.allatoms(j)%name(1:8).and.
!     &          allatoms(i)%core.eq.'shel'
!     &          .and.allatoms(j)%core.eq.'shel')
!     &         newele=.false.
!            end do
!            if (newele) neles=neles+1
!      end do
      !write(12,'(8x,"Found ",I4," shel elements.")')neles
      allocate(eless(1:neles))
      ! get names of elements
      ! cores
      ielec=1
      ieles=1
      do i=1,natoms
            newelec=.true.
            neweles=.true.
            if (allatoms(i)%core(1:4).eq.'shel') newelec=.false.
            if (allatoms(i)%core(1:4).eq.'core') neweles=.false.
            do j=1,i-1
              if (allatoms(i)%name(1:8).eq.allatoms(j)%name(1:8).and.
     &           allatoms(i)%core(1:4).eq.'core'.and.
     &           allatoms(j)%core(1:4).eq.'core') newelec=.false.
              if (allatoms(i)%name(1:8).eq.allatoms(j)%name(1:8).and.
     &           allatoms(i)%core(1:4).eq.'shel'.and.
     &           allatoms(j)%core(1:4).eq.'shel') neweles=.false.
!              if (allatoms(i)%name(1:8).eq.allatoms(j)%name(1:8).and.
!     &            allatoms(i)%core.eq.'core')
!     &         newele=.false.
            end do
            if (newelec) then
                  elesc(ielec)%name=allatoms(i)%name
                  ielec=ielec+1
            end if
            if (neweles) then
                  eless(ieles)%name=allatoms(i)%name
                  ieles=ieles+1
            end if
      end do
!      ! shels
!      iele=1
!      do i=1,natoms
!            newele=.true.
!            do j=1,i-1
!              if (allatoms(i)%name(1:8).eq.allatoms(j)%name(1:8).and.
!     &            allatoms(i)%core.eq.'shel')
!     &         newele=.false.
!            end do
!            if (newele) then
!                  eless(iele)%name=allatoms(i)%name
!                  iele=iele+1
!            end if
!      end do
      ! get numbers of elements
      ! cores
      do ielec=1,nelec
            numelec=0
            do j=1,natoms
               if(allatoms(j)%name.eq.elesc(ielec)%name) then
                  if(allatoms(j)%core.eq.'core')numelec=numelec+1
               end if
            end do
            elesc(ielec)%howmany=numelec
      end do
      ! shels
      do ieles=1,neles
            numeles=0
            do j=1,natoms
               if(allatoms(j)%name.eq.eless(ieles)%name) then
                  if(allatoms(j)%core.eq.'shel')numeles=numeles+1
               end if
            end do
            eless(ieles)%howmany=numeles
      end do
      write(12,'(8x,"Found ",I4," core elements.")')nelec
      write(12,'(8x,"Elements, number of atoms of each element:")')
      do i=1,nelec
            write(12,'(8x,A8,I6)') elesc(i)%name,elesc(i)%howmany
      end do
      write(12,'(8x,"Found ",I4," shel elements.")')neles
      write(12,'(8x,"Elements, number of atoms of each element:")')
      do i=1,neles
            write(12,'(8x,A8,I6)') eless(i)%name,eless(i)%howmany
      end do
      !write COORAT files
      !cores
      do ielec=1,nelec
            write(18,'(4x,"natom= ",I6," name= ",A8)')
     &       elesc(ielec)%howmany,elesc(ielec)%name
            do i=1,natoms
                  if(allatoms(i)%name.eq.elesc(ielec)%name) then
                        if(allatoms(i)%core.eq.'core')   
     &               write(18,'(3(F15.10))') allatoms(i)%where
                  end if
            end do
      end do
      ! shels
      do ieles=1,neles
            write(19,'(4x,"natom= ",I6," name= ",A8)')
     &       eless(ieles)%howmany,eless(ieles)%name
            do i=1,natoms
                  if(allatoms(i)%name.eq.eless(ieles)%name) then
                        if(allatoms(i)%core.eq.'shel')   
     &               write(19,'(3(F15.10))') allatoms(i)%where
                  end if
            end do
      end do
           
      close(17)
      close(18)
      close(19)
      write(12,*) '     '

      write(12,fsubend)
      end subroutine

c---------------------------------------------------------------------

      end module
