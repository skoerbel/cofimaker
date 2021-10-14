      module md
      implicit none

      contains

c---------------------------------------------------------------------

      subroutine mdana(trjfile,trjformat,trjatomfile,trjbndry)
      ! reads MD data from GULP trajectory output file
      
      use misc
      use defs
      implicit none
     
      double precision trjbndry(1:3)
      character (len=*) trjfile,trjatomfile,trjformat
      
      ! local variables
      integer i,j,natoms,nmirror,ndata
      ! for DL_POLY FILED file
      integer imol,nmolecules,natomsmol  
      integer ntrjdatalines,itimestep
      double precision timestep

      integer ispecies,nspecies,mult,nummols
      type(element), allocatable :: species(:)
      double precision time,ekin,etot,temp,cell(1:3,1:3),alat(1:3)
      double precision angles(1:3),pol(1:3),polav(1:3),vol,totalcharge
      double precision pol2(1:3)
      double precision  mirrorcoords(1:8,1:3)
      type(atom), allocatable :: allatoms(:)
      character line*100
      character char1*20,char2*20
      logical isopen12,MDOUTexists
      !character FMT1*1024,FMT2*1024

      ! say hello
      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
         write(12,fsubstart) trim(adjustl('mdana'))
         write(12,'(8x,"Input file is ",A)')trim(trjfile)
         write(12,'(8x,"Atom prop. are in ",A,".")')trim(trjatomfile)
         write(12,'(8x,"Format of traj. file:",A)')trim(trjformat)
         write(12,*)'       '
      else
         print fsubstart, trim(adjustl('mdana'))
         print '(8x,"Input file is ",A)', trim(trjfile)
         print '(8x,"Atom prop. are in ",A,".")',trim(trjatomfile)
         print '(8x,"Format of traj. file:",A)',trim(trjformat)
         print*, '       '
      end if

      ! read number of atoms 
      select case (trjformat)
      case ('trg')
        open(23,file=trjfile, status='old',err=230)
        read(23,'(A100)',end=230,err=230) line
        read(23,*,end=230,err=230) natoms
      case('hist','history','HIST','HISTORY')
        open(23,file=trjfile, status='old',err=231)
        read(23,'(A100)',end=231,err=231) line
        read(23,'(A100)',end=231,err=231) line
        read(line(1:10),*,err=231) ntrjdatalines
        if (ntrjdatalines.ge.1) print '(8x,"Will read pos.")'
        if (ntrjdatalines.ge.2) print '(8x,"Will read veloc. also.")'
        if (ntrjdatalines.ge.3) print '(8x,"Will read forces also.")'
        read(line(21:100),*,err=231) natoms
      case default
        if(isopen12) then
          write(12,'(8x,"Error: unknowm format of traj. file: ",A)'),
     &         trim(trjformat)
        end if
        print ferrmssg, "unknowm format of traj. file: ",
     &         trim(trjformat)
        return
      end select

      allocate(allatoms(natoms))
      if(isopen12) then
!            write(12,'(8x,"Found ",I," atoms.")') natoms
            WRITE(FMT2,*) natoms
            WRITE(FMT1,*) '(8x,"Found ",I',len(trim(FMT2)),
     &        '," atoms.")'
            write(12,FMT1) natoms 
      else
!            print '(8x,"Found ",I," atoms.")', natoms
            WRITE(FMT2,*) natoms
            WRITE(FMT1,*) '(8x,"Found ",I',len(trim(FMT2)),
     &        '," atoms.")'
            print FMT1, natoms 
      end if

!      select case(trjformat)
!      case('trg')
!            goto 9
!      case default
!            print '("Error: Unknown format of atom file. Quitting.")'
!            return
!      end select

!9      open(24,file=trjatomfile,status='old')
      allatoms%name='        '
      select case(trjformat)
      case('trg')
            open(24,file=trjatomfile,status='old')
            goto 10
      case ('hist','HIST','history','HISTORY')
            open(24,file='FIELD',status='old')
            goto 11
      end select

      ! if GULP trg: read order of atoms from GULP restart file
10    read(24,'(A100)',end=15,err=240) line
      line=adjustl(line)
      select case(trjformat)
      case('trg')
        if (index(line,'frac').ge.1.or.index(line,'cart').ge.1) then
          i=1
          totalcharge=0.0D0
          do while (i.le.natoms)
                read(24,'(A100)') line
                line=adjustl(line)
                if (line(1:1).ne.'#') then
                      read(line,*) allatoms(i)%name
                      read(line(len_trim(allatoms(i)%name)+1:100),*)
     &                allatoms(i)%core,allatoms(i)%where(1:3),
     &                allatoms(i)%charge  
                      totalcharge=totalcharge+allatoms(i)%charge
                      i=i+1
                end if
          end do
          if(isopen12) then
                write(12,'(A23,F10.5)')'         Total charge:',
     &            totalcharge
          else
                print '(A23,F10.5)','         Total charge:',
     &            totalcharge
          end if
          goto 15
        end if
      case default
      end select
      goto 10

11    continue

15    continue
      close(24)
      ! write header of output file
      print '(8x,"Will write MD output now.")'
      open(25,file="TRAJ.DAT",status='replace')
            write(25,'(14(A15))') '# t (ps)', 'ekin (eV)','etot (eV)',
     &'T (K)',' a ','  b ', ' c (Angs)', 'alpha', 'beta', 'gamma (°)', 
     & 'vol (Angs³)','P_x', 'P_y', 'P_z (uC/cm²)'      
      if(isopen12) then
            write(12,'(8x,"Boundary tolerance: ",3(F8.4))') trjbndry
      else
            print '(8x,"Boundary tolerance: ",3(F8.4))', trjbndry
      end if

      ! read trajectory and calculate polarization
      ndata=0
      polav=0.0D0
      pol2=0.0D0
20    read(23,'(A100)',end=25) line
      ndata=ndata+1
      select case(trjformat)
      case('trg')
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
      case('hist','history','HIST','HISTORY')
        read(line(9:18),*) itimestep
        read(line(49:60),*) timestep
        time=dble(itimestep)*timestep
        do i=1,3
          read(23,*) cell(i,1:3)
        end do
        do i=1,natoms
          read(23,'(A100)',end=25) line
          read(line(1:8),*) allatoms(i)%name
          allatoms(i)%core='core'
          if(index(allatoms(i)%name,'shel').gt.0) 
     &        allatoms(i)%core='shel'
          read(line(31:42),*) allatoms(i)%charge
          read(23,*) allatoms(i)%where(1:3)
          do j=1,ntrjdatalines
            read(23,'(A100)') line
          end do
        end do
      case default
      end select
      ! calculate the cell volume from the cell vectors
      vol=cell(1,1)*(cell(2,2)*cell(3,3)
     &              -cell(2,3)*cell(3,2))
     &   +cell(1,2)*(cell(2,3)*cell(3,1)
     &              -cell(2,1)*cell(3,3))
     &   +cell(1,3)*(cell(2,1)*cell(3,2)
     &              -cell(2,2)*cell(3,1))
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
      ! if GULP trg: now comes the strain. 
      if (trjformat.eq.'trg') then
        read(23,'(A100)') line   
        read(23,'(A100)') line   
        read(23,'(A100)') line   
      end if
      ! calculate the polarization
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
     &       "atom multipl. & coords of first structure:")') 
          else
                print '(8x,
     &       "atom multipl. & coords of first structure:")' 
          end if
      end if
      do i=1,natoms
        call mirror(allatoms(i)%where,cell,nmirror,mirrorcoords,
     &         trjbndry)
        if(ndata.eq.1) then
          if(isopen12) then
!            write(12,'(8x,A3,I)') allatoms(i)%name,nmirror
            WRITE(FMT2,*) nmirror
            WRITE(FMT1,*) '(8x,A3,I',len(trim(FMT2)),')'
            write(12,FMT1)  allatoms(i)%name,nmirror
          else
!            print '(8x,A3,I)', allatoms(i)%name,nmirror
            WRITE(FMT2,*) nmirror
            WRITE(FMT1,*) '(8x,A3,I',len(trim(FMT2)),')'
            print FMT1, allatoms(i)%name,nmirror
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
     &             *allatoms(i)%charge/dble(nmirror)
        end do  
      end do
      !
      ! crudely correct for net charge
      do i=1,3
        pol=pol-totalcharge*0.5D0*
     &    (cell(1,i)+cell(2,i)+cell(3,i))
      end do
!      if (isopen12) then
!        write(12,*)'        Total dipole moment (e*Angs):',pol
!      else  
!        print '("        Total dipole moment (e*Angs):",3(F10.4))'
!         ,pol
!      end if
      pol=pol*polfac_Angs_e/vol
      polav=polav+pol
      pol2=pol2+pol**2
      write(25,'(14(F15.6))') time,ekin,etot,temp,
     &     (alat(i),i=1,3),(angles(i),i=1,3),vol,pol(1:3)
      goto 20

25    continue
      close(23)
      ! averqged pol. = sum(pol vaues)/number of pol values
      polav=polav/dble(ndata)
      pol2=pol2/dble(ndata)
      ! print averaged pol. to file 
      write(25,'("# time-av. polarization in uC/cm^2:",3(F10.4))'),
     &        polav(1:3)
      write(25,'("# root of mean square             :",3(F10.4))'),
     &        sqrt(pol2(1:3))
      write(25,'("# variance                        :",3(F10.4))'),
     &        pol2(1:3)-polav(1:3)**2
      write(25,'("# standard deviation              :",3(F10.4))'),
     &        sqrt(pol2(1:3)-polav(1:3)**2)
      close(25)

      ! get MD data from other MD output files
      select case(trjformat)
      case ('trg') 
          ! if GULP, then get data from GULP MD output file
          INQUIRE(file="GOUT", exist=MDOUTexists)
          if (MDOUTexists) then
            print '(8x,"Will write MD output in GOUT to file
     &        MDOUT.DAT")'
            call gulpmdout()
          end if
      case('hist','history','HIST','HISTORY')
          ! if DL_POLY, then get data from DL_POLY MD output file
          INQUIRE(file="OUTPUT", exist=MDOUTexists)
          if (MDOUTexists) then
            print '(8x,"Will write MD output in OUTPUT to file
     &        MDOUT.DAT")'
            call dlpolymdout()
          end if
      end select

      if(isopen12) then
        write(12,fsubend)
      else
        print(fsubend)
      end if
      return

100   print*,"Error: something wrong with INPUT.COFIMA."   
      if(isopen12) then
        write(12,'(8x,"Error: something wrong with INPUT.COFIMA.")')
      else
        print '(8x,"Error: something wrong with INPUT.COFIMA.")'
      endif
      nerr=nerr+1
      return

230   print*,"Error: something wrong with the GULP trg file."   
      if(isopen12) then
        write(12,'(8x,"Error: something wrong with the GULP")')
        write(12,'(8x,"trg file.")')
      endif
      nerr=nerr+1
      return

231   print*,"Error: something wrong with the DL_POLY HISTORY file."   
      if(isopen12) then
        write(12,'(8x,"Error: something wrong with the DL_POLY")')
        write(12,'(8x,"HISTORY file.")')
      endif
      nerr=nerr+1
      return

240   print*,"Error: something wrong with the GULP restart file."   
      if(isopen12) then
        write(12,'(8x,"Error: something wrong with the GULP")')
        write(12,'(8x,"restart file.")')
      endif
      nerr=nerr+1
      return

241   print*,"Error: something wrong with the DL_POLY FIELD file."   
      if(isopen12) then
        write(12,'(8x,"Error: something wrong with the DL_POLY")')
        write(12,'(8x,"FIELD file.")')
      endif
      nerr=nerr+1
      return

      end subroutine 

c---------------------------------------------------------------------

      subroutine gulpmdout()
      ! extracts temperature and pressure etc from gulp MD output file 

      use defs
      implicit none
      
      ! variables
      
      !local variables
      double precision time,ekin,epot,etot,temp,tcs,press,alat(1:3),
     &            angle(1:3),vol,ekinav,epotav,etotav,tav,tcsav,pav,
     &            alatav(1:3),angleav(1:3),volav
      double precision ekin2,epot2,etot2,temp2,tcs2,press2,
     &            alat2(1:3),angle2(1:3),vol2
      integer i,ndata  ! ,np,nt
      character line*150,time_char*13      !,infilename*20,outfilename*20
      !character outfilename2*30
      double precision eq_time
      logical isopen12
      

      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
        write(12,fsubstart)trim(adjustl('gulpmdout'))
      else
        print fsubstart,trim(adjustl('gulpmdout'))
      end if

      !read(10,*) infilename
      !read(10,*) outfilename
      !outfilename2=trim(adjustl(outfilename))
      !outfilename2(len_trim(adjustl(outfilename))+1:30)='.AV'
      open(26,file="GOUT",status='old')
      open(27,file="MDOUT.DAT",status='replace')
      open(28,file="MDAV.DAT",status='replace')
      
      eq_time=0.0D0
      write(27,1001) 't (ps)','ekin','epot','etot','T','T c/s',
     &      'p (GPa)',
     &      ' a','b','c','alpha','beta','gamma','V'
      write(28,1001) 't (ps)','ekin av',
     &      'epotav','etotav','Tav','Tc/sav','pav','aav',
     &      'bav','cav','alphaav','betaav','gammaav','volav'
      ndata=0
      ekin2=0.0D0
      epot2=0.0D0
      etot2=0.0D0
      temp2=0.0D0
      tcs2=0.0D0
      press2=0.0D0
      do i=1,3
        alat2(i)=0.0D0
        angle2(i)=0.0D0
      end do
      vol2=0.0D0
10    read(26,'(A150)',end=20) line
      if (index(line,'  ** Time : ').ge.1) then
            read(line(12:150),*) time
            ndata=ndata+1
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
            ekin2=ekin2+ekin**2
            epot2=epot2+epot**2
            etot2=etot2+etot**2
            temp2=temp2+temp**2
            tcs2=tcs2+tcs**2
            press2=press2+press**2
            do i=1,3
              alat2(i)=alat2(i)+alat(i)**2
              angle2(i)=angle2(i)+angle(i)**2
            end do
            vol2=vol2+vol**2
      end if
      
      if (index(line,'  Molecular dynamics production :').ge.1) then
            eq_time=time
            ! write root of mean squares
            write(28,'("# root of mean square:")')
            write(28,2000) sqrt(ekin2/dble(ndata)),
     &                     sqrt(epot2/dble(ndata)),
     &                     sqrt(etot2/dble(ndata)),
     &                     sqrt(temp2/dble(ndata)),
     &                     sqrt(tcs2/dble(ndata)),
     &                     sqrt(press2/dble(ndata)),
     &                     sqrt(alat2/dble(ndata)),
     &                     sqrt(angle2/dble(ndata)),
     &                     sqrt(vol2/dble(ndata))
            ! write variance
            write(28,'("# variance:")')
            write(28,2000) ekin2/dble(ndata)-ekinav**2,
     &                     epot2/dble(ndata)-epotav**2,
     &                     etot2/dble(ndata)-etotav**2,
     &                     temp2/dble(ndata)-tav**2,
     &                     tcs2/dble(ndata)-tcsav**2,
     &                     press2/dble(ndata)-pav**2,
     &                     alat2/dble(ndata)-alatav**2,
     &                     angle2/dble(ndata)-angleav**2,
     &                     vol2/dble(ndata)-volav**2
            ! write standard deviations
            write(28,'("# standard deviations:")')
            write(28,2000) sqrt((ekin2/dble(ndata)-ekinav**2)),
     &                     sqrt((epot2/dble(ndata)-epotav**2)),
     &                     sqrt((etot2/dble(ndata)-etotav**2)),
     &                     sqrt((temp2/dble(ndata)-tav**2)),
     &                     sqrt((tcs2/dble(ndata)-tcsav**2)),
     &                     sqrt((press2/dble(ndata)-pav**2)),
     &                     sqrt((alat2/dble(ndata)-alatav**2)),
     &                     sqrt((angle2/dble(ndata)-angleav**2)),
     &                     sqrt((vol2/dble(ndata)-volav**2))
            goto 20
      end if
      goto 10

20    continue

      ndata=0
      ekin2=0.0D0
      epot2=0.0D0
      etot2=0.0D0
      temp2=0.0D0
      tcs2=0.0D0
      press2=0.0D0
      do i=1,3
        alat2(i)=0.0D0
        angle2(i)=0.0D0
      end do
      vol2=0.0D0
30    read(26,'(A150)',end=40) line
      if (index(line,'  ** Time : ').ge.1) then
            read(line(12:150),*) time
            ndata=ndata+1
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
            ekin2=ekin2+ekin**2
            epot2=epot2+epot**2
            etot2=etot2+etot**2
            temp2=temp2+temp**2
            tcs2=tcs2+tcs**2
            press2=press2+press**2
            do i=1,3
              alat2(i)=alat2(i)+alat(i)**2
              angle2(i)=angle2(i)+angle(i)**2
            end do
            vol2=vol2+vol**2
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
      ! write root of mean squares
      write(28,'("# root of mean square:")')
      write(28,2000) sqrt(ekin2/dble(ndata)),
     &               sqrt(epot2/dble(ndata)),
     &               sqrt(etot2/dble(ndata)),
     &               sqrt(temp2/dble(ndata)),
     &               sqrt(tcs2/dble(ndata)),
     &               sqrt(press2/dble(ndata)),
     &               sqrt(alat2/dble(ndata)),
     &               sqrt(angle2/dble(ndata)),
     &               sqrt(vol2/dble(ndata))
      ! write variance
      write(28,'("# variance:")')
      write(28,2000) ekin2/dble(ndata)-ekinav**2,
     &               epot2/dble(ndata)-epotav**2,
     &               etot2/dble(ndata)-etotav**2,
     &               temp2/dble(ndata)-tav**2,
     &               tcs2/dble(ndata)-tcsav**2,
     &               press2/dble(ndata)-pav**2,
     &               alat2/dble(ndata)-alatav**2,
     &               angle2/dble(ndata)-angleav**2,
     &               vol2/dble(ndata)-volav**2
      ! write standard deviations
      write(28,'("# standard deviation:")')
      write(28,2000) sqrt((ekin2/dble(ndata)-ekinav**2)),
     &               sqrt((epot2/dble(ndata)-epotav**2)),
     &               sqrt((etot2/dble(ndata)-etotav**2)),
     &               sqrt((temp2/dble(ndata)-tav**2)),
     &               sqrt((tcs2/dble(ndata)-tcsav**2)),
     &               sqrt((press2/dble(ndata)-pav**2)),
     &           sqrt(alat2(1:3)/dble(ndata)-alatav(1:3)**2),
     &          sqrt(angle2(1:3)/dble(ndata)-angleav(1:3)**2),
     &               sqrt((vol2/dble(ndata)-volav**2))
      close(26)
      close(27)
      close(28)

! 1000 format(14(F10.5))     
 1000 format(F10.5,F12.5,2(F15.5),F12.5,8(F10.5),F15.5)     
! 1001 format(14(A10))
 1001 format(A8,A12,2(A15),A12,8(A10),A15)     
 2000 format("#       ",2x,F12.5,2(F15.5),F12.5,8(F10.5),F15.5)     
      
      if(isopen12) then
        write(12,fsubend)
      else
        print fsubend
      end if
     
      end subroutine

c---------------------------------------------------------------------

      subroutine dlpolymdout()
      ! extracts temperature and pressure etc from DL_POLY MD output file 

      use defs
      implicit none
      
      ! variables

      !local variables
      double precision time,ekin,epot,etot,temp,tcs,press,alat(1:3),
     &            angle(1:3),vol,ekinav,epotav,etotav,tav,tcsav,pav,
     &            alatav(1:3),angleav(1:3),volav,cell(1:3,1:3)
      double precision ekin2,epot2,etot2,temp2,tcs2,press2,
     &            alat2(1:3),angle2(1:3),vol2
      integer ndata,i ! np,nt
      character line*150,time_char*13      !,infilename*20,outfilename*20
      double precision eq_time
      logical isopen12
      

      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
        write(12,fsubstart)trim(adjustl('dlpolymdout'))
      else
        print fsubstart,trim(adjustl('dlpolymdout'))
      end if

      !read(10,*) infilename
      !read(10,*) outfilename
      !outfilename2=trim(adjustl(outfilename))
      !outfilename2(len_trim(adjustl(outfilename))+1:30)='.AV'
      open(26,file="STATIS",status='old')
      open(27,file="MDOUT.DAT",status='replace')
      open(28,file="MDAV.DAT",status='replace')
      
      eq_time=0.0D0
      write(27,1001) 't (ps)','ekin','epot','etot','T','T c/s',
     &      'p (GPa)',
     &      ' a','b','c','alpha','beta','gamma','V'
      write(28,1001) 't (ps)','ekin av',
     &      'epotav','etotav','Tav','Tc/sav','pav','aav',
     &      'bav','cav','alphaav','betaav','gammaav','volav'
      ! first and second line are the title and units of the simulation
      read(26,'(A150)',end=20) line
      read(26,'(A150)',end=20) line
      ! print warning if energy is not in eV
      if(index(line,'electron Volts').le.0) then
            print '(8x,"Warning: energy unit is not eV.")'
      end if
      alatav=0.0D0
      alat2=0.0D0
      angleav=0.0D0
      angle2=0.0D0
      volav=0.0D0
      vol2=0.0D0
      pav=0.0D0
      press2=0.0D0
      ekinav=0.0D0
      ekin2=0.0D0
      epotav=0.0D0
      epot2=0.0D0
      etotav=0.0D0
      etot2=0.0D0
      tav=0.0D0
      temp2=0.0D0
      tcsav=0.0D0
      tcs2=0.0D0
      ndata=0
10    read(26,'(A150)',end=20) line 
      ndata=ndata+1
      ! first data line: time
      read(line(11:25),*) time
      
      read(26,'(A150)',end=20) line
      ! 2nd data line: etot,temp,epot 
      read(line(1:14),*) etot
      etotav=etotav+etot
      etot2=etot2+etot**2
      read(line(15:28),*) temp
      tav=tav+temp
      temp2=temp2+temp**2
      read(line(16:42),*) epot
      epotav=epotav+epot
      epot2=epot2+epot**2
      
      read(26,'(A150)',end=20) line
      read(26,'(A150)',end=20) line
      ! 3rd and 4th data line: not so interesting. Skip  
      
      read(26,'(A150)',end=20) line
      ! 5th data line: vol, temp_core/shel  
      read(line(43:70),*) vol,tcs
      volav=volav+vol
      vol2=vol2+vol**2
      tcsav=tcsav+tcs
      tcs2=tcs2+tcs**2
      
      read(26,'(A150)',end=20) line
      ! 6th data line: cell angles  
      read(line(29:70),*) angle(1:3)
      angleav=angleav+angle
      do i=1,3
        angle2(i)=angle2(i)+angle(i)**2
      end do
      
      read(26,'(A150)',end=20) line
      ! 7th data line: pressure  
      read(line(15:28),*) press
      ! DL_POLY gives pressure in kilo atm. We want GPa:
      press=press*1.0133D0/10.0D0
      pav=pav+press
      press2=press2+press**2
      
      do i=8,9
        read(26,'(A150)',end=20) line
      end do

      read(26,'(A150)',end=20) line
      ! 10th to 12th line: cell parameters
      read(line(29:150),*) cell(1,1:3)
      read(26,'(A150)',end=20) line
      read(line(1:150),*) cell(2,1:3),cell(3,1:2)
      read(26,'(A150)',end=20) line
      read(line(1:150),*) cell(3,3)
      do i=1,3
        alat(i)=sqrt(cell(i,1)**2+cell(i,2)**2+cell(i,3)**2)
        alatav(i)=alatav(i)+alat(i)
        alat2(i)=alat2(i)+alat(i)**2
      end do
            
!      end if
            write(27,1000) time,ekin,epot,etot,temp,tcs,press,alat,
     &            angle,vol
            write(28,1000) time,ekinav/dble(ndata),epotav/dble(ndata),
     &                    etotav/dble(ndata),tav/dble(ndata),
     &                    tcsav/dble(ndata),pav/dble(ndata),
     &                    alatav/dble(ndata),angleav/dble(ndata),
     &                    volav/dble(ndata)
!      end if
      
!      if (index(line,'  Molecular dynamics production :').ge.1) then
!            eq_time=time
!            goto 20
!      end if
      goto 10

20    continue

      ! write root of mean squares
      write(28,'("# root of mean square:")')
      write(28,2000) sqrt(ekin2/dble(ndata)),
     &               sqrt(epot2/dble(ndata)),
     &               sqrt(etot2/dble(ndata)),
     &               sqrt(temp2/dble(ndata)),
     &               sqrt(tcs2/dble(ndata)),
     &               sqrt(press2/dble(ndata)),
     &               sqrt(alat2/dble(ndata)),
     &               sqrt(angle2/dble(ndata)),
     &               sqrt(vol2/dble(ndata))
      ! write variance
      write(28,'("# variance:")')
      write(28,2000) (ekin2-ekinav**2/dble(ndata))/dble(ndata),
     &               (epot2-epotav**2/dble(ndata))/dble(ndata),
     &               (etot2-etotav**2/dble(ndata))/dble(ndata),
     &               (temp2-tav**2/dble(ndata))/dble(ndata),
     &               (tcs2-tcsav**2/dble(ndata))/dble(ndata),
     &               (press2-pav**2/dble(ndata))/dble(ndata),
     &               (alat2-alatav**2/dble(ndata))/dble(ndata),
     &               (angle2-angleav**2/dble(ndata))/dble(ndata),
     &               (vol2-volav**2/dble(ndata))/dble(ndata)
c      write(28,2000) ekin2/dble(ndata)-(ekinav/dble(ndata))**2,
c     &               epot2/dble(ndata)-epotav**2,
c     &               etot2/dble(ndata)-etotav**2,
c     &               temp2/dble(ndata)-tav**2,
c     &               tcs2/dble(ndata)-tcsav**2,
c     &               press2/dble(ndata)-pav**2,
c     &               alat2(1:3)/dble(ndata)-alatav(1:3)**2,
c     &               angle2(1:3)/dble(ndata)-angleav(1:3)**2,
c     &               vol2/dble(ndata)-volav**2
      ! write standard deviations
      write(28,'("# standard deviations:")')
      write(28,2000) sqrt((ekin2-ekinav**2/dble(ndata))/dble(ndata)),
     &               sqrt((epot2-epotav**2/dble(ndata))/dble(ndata)),
     &               sqrt((etot2-etotav**2/dble(ndata))/dble(ndata)),
     &               sqrt((temp2-tav**2/dble(ndata))/dble(ndata)),
     &               sqrt((tcs2-tcsav**2/dble(ndata))/dble(ndata)),
     &               sqrt((press2-pav**2/dble(ndata))/dble(ndata)),
     &               sqrt((alat2-alatav**2/dble(ndata))/dble(ndata)),
     &               sqrt((angle2-angleav**2/dble(ndata))/dble(ndata)),
     &               sqrt((vol2-volav**2/dble(ndata))/dble(ndata))

40    continue
      close(26)
      close(27)
      close(28)

! 1000 format(14(F10.5))     
 1000 format(F10.5,F12.5,2(F15.5),F12.5,8(F10.5),F15.5)     
! 1001 format(14(A10))
 1001 format(A8,A12,2(A15),A12,8(A10),A15)     
 2000 format("#       ",2x,F12.5,2(F15.5),F12.5,8(F10.5),F15.5)     
      
      if(isopen12) then
        write(12,fsubend)
      else
        print fsubend
      end if
     
      end subroutine

c---------------------------------------------------------------------

      subroutine lmpthermo(lmplog)
      use defs
      implicit none
      ! global
      character(len=*) lmplog
      ! local
      character line*256
      integer i
      logical isopen12

      ! write something to output file
      INQUIRE (unit=12, opened=isopen12)
      if(isopen12) then
        write(12,fsubstart) 'lmpthermo'
        write(12,'(8x,"writing MD data from LAMMPS logfile to",
     &   " MD.DAT.")')
      else
        print fsubstart, 'lmpthermo'
        print '(8x,"writing MD data from LAMMPS logfile to",
     &   " MD.DAT.")'
      end if

      open(26,file=lmplog,status="old",err=1000)
      rewind(26)
      open(27,file="MD.DAT",status="replace")
      ! read LAMMPS logfile and write MD data to "MD.DAT"
10    read(26,'(A256)',end=15,err=1000) line   
      line=adjustl(line)
      if(line(1:28).eq."Memory usage per processor =") then
         read(26,'(A256)',end=15,err=1000) line
         !write(27,'(A1,A255)') '#',line(1:255)
         write(27,'(A1,A)') '#',trim(line)
         read(26,'(A256)',end=15,err=1000) line
         do while(line(1:13).ne."Loop time of ")   
           !write(27,'(A256)') line(1:256)
           write(27,'(A)') trim(line)
           read(26,'(A256)',end=15,err=1000) line
         end do
      end if
      goto 10

15    continue      

      ! exit normally
      if(isopen12) then
        write(12,fsubend)
      else
        print fsubend
      end if
      close(26)
      close(27)
      return

      ! exit with an error
1000  nerr=nerr+1
      if(isopen12) then
        write(12,ferrmssg) 'Something wrong with LAMMPS logfile'
        write(12,fsubend)
      else
        print ferrmssg, 'Something wrong with LAMMPS logfile'
        print fsubend
      end if
      close(26)
      close(27)
      return

      end subroutine

c---------------------------------------------------------------------

      subroutine gulpnebana(infile)
      ! write structures and energies of GULP NEB images to files
      ! "NEBGEO.XSF" and "NEBENERGIES.DAT"
      use defs
      use misc
      implicit none
      
      ! global variables
      character(len=*) infile

      ! local variables
      integer i,i_rep,n_rep,n_at,nmin,nmax
      double precision aa,bb,cc,alph,beta,gamm,xxf,yyf,zzf,xx,yy,zz
      double precision avecs(1:3,1:3),energy,emin,emax
      double precision, allocatable :: initcoords(:,:),coords(:,:)
      double precision, allocatable :: reactcoord(:),neben(:)
      character line*120,ele*2,core*4
      character,allocatable :: elelist(:)*2
      logical isopen12

      ! write hello
      inquire(unit=12,opened=isopen12)
      if(isopen12) then
            write(12,fsubstart) "gulpnebana"
      else
            print fsubstart, "gulpnebana"
      end if

      ! open GULP NEB output file and read structures and energies 
      ! open files to write structures and energies 
      elelist=' '
      open(51,file=infile,status='old',err=100)
      open(52,file="NEBGEO.XSF",status='replace',err=110)
10    read(51,1000,end=11,err=100) line
      ! read number of atoms
      if (index(line,'Formula =')
     & .gt.0) then
            read(51,1000,err=100) line
            read(51,'(A38,I8)',err=100) line,n_at
            allocate(elelist(n_at))
            allocate(initcoords(n_at,1:3),coords(n_at,1:3))
      end if
      ! read number of images
      if (index(line,'Nudged Elastic Band')
     & .gt.0) then
            read(51,1000,err=100) line
            read(51,'(A22,I7)',err=100) line,n_rep
            write(52,*,err=110)  'ANIMSTEPS ',n_rep+2 
      end if
      allocate(reactcoord(n_rep+2),neben(n_rep+2))
      ! read and write energies of images
      if (index(line,'Energy of NEB configurations :')
     & .gt.0) then
            open(53,file='NEBENERGIES.DAT',status='replace')
            read(51,1000,err=100) line
            read(51,'(A23,F21.8)',err=100) line,energy
	    emin=energy
	    emax=energy
	    nmin=0
	    nmax=0
            write(53,1600,err=110)  0 ,energy 
            do i=1,n_rep
                  read(51,'(A23,F21.8)',err=100) line,energy
		  if (energy.lt.emin) then
		    emin=energy
		    nmin=i
		  end if
		  if (energy.gt.emax) then
		    emax=energy
		    nmax=i
		  end if
                  write(53,1600)  i ,energy 
            end do
            read(51,'(A23,F21.8)',err=100) line,energy
            write(53,1600,err=110)  n_rep+1 ,energy 
	    if (energy.lt.emin) then
	      emin=energy
	      nmin=n_rep+1
	    end if
	    if (energy.gt.emax) then
	      emax=energy
	      nmax=n_rep+1
	    end if
	    write(53,*) '# image with lowest energy:',nmin
	    write(53,*) '# image with highest energy:',nmax
	    write(53,'(A30,F9.6,A3)') '# energy difference:',
     &            emax-emin,' eV'
            close(53)
            goto 11
      end if
      goto 10

11    continue
      rewind(51)
12    read(51,1000,end=13,err=100) line
      ! get initial structure
      if (index(line,'Cartesian lattice vectors (Angstroms) :')
     & .gt.0) then
            read(51,1000,err=100) line
            do i=1,3
                  read(51,*,err=100) avecs(i,1:3)
            end do
            write(52,*,err=110)  'CRYSTAL'
            write(52,*,err=110)  'PRIMVEC ','1'
            write(52,1500,err=110) avecs(1,1:3)
            write(52,1500,err=110) avecs(2,1:3)
            write(52,1500,err=110) avecs(3,1:3)
      end if
      if (index(line,'Fractional coordinates of asymmetric unit :')
     & .gt.0) then
            do i=1,5
                  read(51,1000,err=100) line
            end do
            write(52,*,err=110) 'PRIMCOORD ','1'
            write(52,*,err=110) n_at,' 1'
            do i=1,n_at
            ! read fractional coordinates and calculates absolute
            ! ones
                 read(51,'(A8,A2,A8,F9.6,A3,F9.6,A3,F9.6,A20)',err=100) 
     &                line,ele,line,xxf,line,yyf,line,zzf,line
                 xx=xxf*avecs(1,1)+yyf*avecs(2,1)+zzf*avecs(3,1)     
                 yy=xxf*avecs(1,2)+yyf*avecs(2,2)+zzf*avecs(3,2)     
                 zz=xxf*avecs(1,3)+yyf*avecs(2,3)+zzf*avecs(3,3)     
                 write(52,3000,err=110) ele,xx,yy,zz
                 initcoords(i,1)=xx
                 initcoords(i,2)=yy
                 initcoords(i,3)=zz
            end do
      end if
      goto 12

13    continue
      rewind(51)
14    read(51,1000,end=15,err=100) line  
      ! look for structures of images
      if (index(line,'Replica coordinates in GULP input format')
     & .gt.0) then
            i_rep=1
            read(51,1000,err=100) line
            do while (i_rep.le.n_rep)
                  write(52,*,err=110)  'CRYSTAL'
                  write(52,*,err=100)  'PRIMVEC ',i_rep+1
                  do i=1,4
                        read(51,*,err=100) line  ! stuff,cell
                  end do
                  ! read cell parameters and calculate cell vectors
                  read(51,*,err=100) aa,bb,cc,alph,beta,gamm
                  avecs=0.0D0
                  avecs(1,1)=aa
                  avecs(2,1)=bb*cos(gamm*pi/180.0D0)
                  avecs(2,2)=bb*sin(gamm*pi/180.0D0)
                  avecs(3,1)=cc*cos(beta*pi/180.0D0)
                  ! avecs(3,2)=cc*cos(alph*pi/180.0D0) ! wrong !
	    	  avecs(3,2)=(cc*bb*cos(alph*Pi/180.0D0)
     &                  -avecs(2,1)*avecs(3,1))/avecs(2,2)
                  avecs(3,3)=sqrt(cc**2-avecs(3,1)**2-avecs(3,2)**2)
                  write(52,1500,err=110) avecs(1,1:3)
                  write(52,1500,err=110) avecs(2,1:3)
                  write(52,1500,err=110) avecs(3,1:3)
                  read(51,1000,err=100) line !fractional
                  write(52,*,err=110) 'PRIMCOORD ',i_rep+1
                  write(52,*,err=110) n_at,' 1'
                  do i=1,n_at
                        ! read fractional coordinates and calculates absolute
                        ! ones, write to output file
                        read(51,*,err=100) ele,core,xxf,yyf,zzf
                        elelist(i)=ele
                        xx=xxf*avecs(1,1)+yyf*avecs(2,1)+zzf*avecs(3,1)     
                        yy=xxf*avecs(1,2)+yyf*avecs(2,2)+zzf*avecs(3,2)     
                        zz=xxf*avecs(1,3)+yyf*avecs(2,3)+zzf*avecs(3,3)     
                        write(52,3000,err=110) ele,xx,yy,zz
                        coords(i,1)=xx
                        coords(i,2)=yy
                        coords(i,3)=zz
                  end do
                  ! get reaction coordinate
                  reactcoord(i_rep)=matnorm(coords-initcoords)
                  i_rep=i_rep+1
            end do
            goto 15
      end if
      goto 14

15    continue
      rewind(51)
16    read(51,1000,end=20,err=100) line
      ! get final structure
      if (index(line,'Final configuration :')
     & .gt.0) then
            read(51,1000,err=100) line
            write(52,*,err=110)  'CRYSTAL'
            write(52,*,err=110)  'PRIMVEC ',n_rep+2
            ! read cell parameters and calculate cell vectors
            read(51,'(A7,3(F11.5),3(F9.4))',err=100)line, aa,bb,cc,alph,
     &              beta,gamm
            avecs=0.0D0
            avecs(1,1)=aa
            avecs(2,1)=bb*cos(gamm*pi/180.0D0)
            avecs(2,2)=bb*sin(gamm*pi/180.0D0)
            avecs(3,1)=cc*cos(beta*pi/180.0D0)
            !avecs(3,2)=cc*cos(alph*pi/180.0D0) ! wrong !
	    avecs(3,2)=(cc*bb*cos(alph*Pi/180.0D0)
     &                  -avecs(2,1)*avecs(3,1))/avecs(2,2)
            avecs(3,3)=sqrt(cc**2-avecs(3,1)**2-avecs(3,2)**2)
            write(52,1500,err=110) avecs(1,1:3)
            write(52,1500,err=110) avecs(2,1:3)
            write(52,1500,err=110) avecs(3,1:3)
            read(51,1000,err=100) line
            write(52,*,err=110) 'PRIMCOORD ',n_rep+2
            write(52,*,err=110) n_at,' 1'
            do i=1,n_at
            ! read fractional coordinates and calculates absolute
            ! ones
                  read(51,'(A7,3(F13.6),A20)',err=100) line,xxf,yyf,zzf,
     &                   line
                  xx=xxf*avecs(1,1)+yyf*avecs(2,1)+zzf*avecs(3,1)     
                  yy=xxf*avecs(1,2)+yyf*avecs(2,2)+zzf*avecs(3,2)     
                  zz=xxf*avecs(1,3)+yyf*avecs(2,3)+zzf*avecs(3,3)     
                  write(52,3000,err=110) elelist(i),xx,yy,zz
            end do
      end if
      goto 16

20    continue
      close(52)

!     end normally
      if(isopen12) then
            write(12,fsubendext) "gulpnebana"
      else
            print fsubendext, "gulpnebana"
      end if
      return

! errors
! error when opening GULP NEB output file for reading
100   nerr=nerr+1
      close(51)
      if (isopen12) then
         write(12,ferrmssg) 
     &    "Something wrong with GULP NEB output file"
      else
         print ferrmssg,
     &    "Something wrong with GULP NEB output file"
      end if
      return
! error when opening output file for writing structures
110   nerr=nerr+1
      close(51)
      close(52)
      close(53)
      if (isopen12) then
         write(12,ferrmssg) 
     &    "Error during writing"
      else
         print ferrmssg,
     &    "Error during writing"
      end if
      return

! formats
 1000 format(A120)     
 1500 format(3(F12.6))  
 1600 format(I3,F15.6)     
 2000 format(6(F12.6))  
 3000 format(A3,3(F12.6))  

      end subroutine     

c---------------------------------------------------------------------

      subroutine nebana(infile,informat)
      ! write structures and energies of NEB images to files
      ! "NEBGEO.XSF" and "NEBENERGIES.DAT"
      use defs
      use misc
      use readcoords
      use writecoords
      implicit none
      
      ! global variables
      character(len=*) infile,informat

      ! local variables
      type(atom),allocatable :: atoms(:),coreatoms(:),shelatoms(:)
      type(element), allocatable :: species(:)
      integer i,i_rep,n_rep,n_at,nmin,nmax,ivec,nspecies,ncores,nshels,
     &  icore,ishel 
      double precision aa,bb,cc,alph,beta,gamm,xxf,yyf,zzf,xx,yy,zz
      double precision vecs(1:3,1:3),energy,emin,emax
      double precision, allocatable :: initcoords(:,:)
      double precision, allocatable :: reactcoord(:),neben(:),
     &   nebcoords(:,:,:),nebvecs(:,:,:)
      character filename*40,line*120,core*4
      logical isopen12,isdir
      !character FMT1*1024,FMT2*1024

      talk=.true.
      ! write hello
      inquire(unit=12,opened=isopen12)
      if(talk) then
        if(isopen12) then
            write(12,fsubstart) "nebana"
        else
            print fsubstart, "nebana"
        end if
      end if

      ! read NEB energies from NEB output file
      select case(informat)
      case('gulp','GULP','gin','GIN')
        print*,"Attempting to read output from GULP."
        ! open GULP NEB output file and read structures and energies 
        ! open files to write structures and energies 
        open(51,file=infile,status='old',err=100)
        rewind(51)
10      read(51,1000,end=11,err=100) line
        ! read number of atoms
        if (index(line,'Formula =')
     &   .gt.0) then
              read(51,1000,err=100) line
              read(51,'(A38,I8)',err=100) line,n_at
              print*,"Found ",n_at," atoms."
              if(.not.allocated(initcoords))
     &         allocate(initcoords(n_at,1:3))
              if(.not.allocated(atoms))allocate(atoms(n_at))
        end if
        ! read number of images
        if (index(line,'Nudged Elastic Band')
     &   .gt.0) then
           read(51,1000,err=100) line
           read(51,'(A22,I7)',err=100) line,n_rep
           print*,"Found ",n_rep," replicae."
           if(.not.allocated(reactcoord))
     &       allocate(reactcoord(n_rep+2),neben(n_rep+2))
           if(.not.allocated(nebcoords))
     &       allocate(nebcoords(n_rep+2,n_at,1:3),
     &       nebvecs(n_rep+2,1:3,1:3))
        end if
        ! read energies of images
        if (index(line,'Energy of NEB configurations :')
     &   .gt.0) then
              read(51,1000,err=100) line
              read(51,'(A23,F21.8)',err=100) line,neben(1)
              emin=neben(1)
              emax=neben(1)
              nmin=0
              nmax=0
              do i=1,n_rep
                  read(51,'(A23,F21.8)',err=100) line,neben(i+1)
          	  if (neben(i+1).lt.emin) then
          	    emin=neben(i+1)
          	    nmin=i
          	  end if
          	  if (neben(i+1).gt.emax) then
          	    emax=neben(i+1)
          	    nmax=i
          	  end if
              end do
              read(51,'(A23,F21.8)',err=100) line,neben(n_rep+2)
              if (neben(n_rep+2).lt.emin) then
                emin=neben(n_rep+2)
                nmin=n_rep+1
              end if
              if (neben(n_rep+2).gt.emax) then
                emax=neben(n_rep+2)
                nmax=n_rep+1
              end if
              goto 11
        end if
        goto 10

11      continue
        rewind(51)
        ! read structures from NEB output file
12      read(51,1000,end=13,err=100) line
        ! get initial structure
        if (index(line,'Cartesian lattice vectors (Angstroms) :')
     &   .gt.0) then
              read(51,1000,err=100) line
              do i=1,3
                  read(51,*,err=100) vecs(i,1:3)
                  nebvecs(1,i,1:3)=vecs(i,1:3)
              end do
        end if
        if (index(line,'Fractional coordinates of asymmetric unit :')
     &   .gt.0) then
            do i=1,5
                  read(51,1000,err=100) line
            end do
            do i=1,n_at
            ! read fractional coordinates and calculate absolute
            ! ones
               read(51,'(A8,A2,A8,F9.6,A3,F9.6,A3,F9.6,A20)',err=100) 
     &            line,atoms(i)%name(1:2),line,atoms(i)%where(1),line,
     &            atoms(i)%where(2),line,atoms(i)%where(3),line
               call frac2abs(atoms(i)%where,vecs,atoms(i)%abswhere)
               initcoords(i,1:3)=atoms(i)%abswhere(1:3)
               nebcoords(1,i,1:3)=atoms(i)%abswhere(1:3)
            end do
            reactcoord(1)=matnorm(nebcoords(1,1:n_at,1:3)
     &       -initcoords)
            ! get species
            call getspecies(atoms,species)
            nspecies=size(species)
        end if
        goto 12
13      continue
        rewind(51)
14      read(51,1000,end=15,err=100) line  
        ! look for structures of images
        if (index(line,'Replica coordinates in GULP input format')
     &   .gt.0) then
            i_rep=1
            read(51,1000,err=100) line
            do while (i_rep.le.n_rep)
               do i=1,4
                     read(51,*,err=100) line  ! stuff,cell
               end do
               ! read cell parameters and calculate cell vectors
               read(51,*,err=100) aa,bb,cc,alph,beta,gamm
               vecs=0.0D0
               vecs(1,1)=aa
               vecs(2,1)=bb*cos(gamm*pi/180.0D0)
               vecs(2,2)=bb*sin(gamm*pi/180.0D0)
               vecs(3,1)=cc*cos(beta*pi/180.0D0)
	       vecs(3,2)=(cc*bb*cos(alph*Pi/180.0D0)
     &               -vecs(2,1)*vecs(3,1))/vecs(2,2)
               vecs(3,3)=sqrt(cc**2-vecs(3,1)**2-vecs(3,2)**2)
               nebvecs(i_rep+1,1:3,1:3)=vecs(1:3,1:3)
               read(51,1000,err=100) line !fractional
               do i=1,n_at
                  ! read fractional coordinates and calculates absolute
                  ! ones, write to output file
                  read(51,*,err=100) atoms(i)%name,atoms(i)%core,
     &               atoms(i)%where(1:3)
                  call frac2abs(atoms(i)%where,vecs,atoms(i)%abswhere)
                  nebcoords(i_rep+1,i,1:3)=atoms(i)%abswhere(1:3)
                  ! shift by cell vectors if atoms have moved outside
                  ! cell edges
                  do ivec=1,3
                    do while (absvec(nebcoords(i_rep+1,i,1:3)
     &                 +vecs(ivec,1:3)-nebcoords(i_rep,i,1:3)).lt.
     &                  absvec(nebcoords(i_rep+1,i,1:3)
     &                 -nebcoords(i_rep,i,1:3)))
                       nebcoords(i_rep+1,i,1:3)=nebcoords(i_rep+1,i,1:3)
     &                    +vecs(ivec,1:3)   
                    end do
                    do while (absvec(nebcoords(i_rep+1,i,1:3)
     &                 -vecs(ivec,1:3)-nebcoords(i_rep,i,1:3)).lt.
     &                  absvec(nebcoords(i_rep+1,i,1:3)
     &                 -nebcoords(i_rep,i,1:3)))
                       nebcoords(i_rep+1,i,1:3)=nebcoords(i_rep+1,i,1:3)
     &                    -vecs(ivec,1:3)
                    end do
                  end do
               end do
               ! get reaction coordinate
               reactcoord(i_rep+1)=matnorm(nebcoords(i_rep+1,1:n_at,1:3)
     &       -initcoords)
               i_rep=i_rep+1
            end do
            goto 15
        end if
        goto 14

15      continue
        rewind(51)
16      read(51,1000,end=17,err=100) line
        ! get final structure
        if (index(line,'Final configuration :')
     &   .gt.0) then
            read(51,1000,err=100) line
            ! read cell parameters and calculate cell vectors
            read(51,'(A7,3(F11.5),3(F9.4))',err=100)line, aa,bb,cc,alph,
     &              beta,gamm
            vecs=0.0D0
            vecs(1,1)=aa
            vecs(2,1)=bb*cos(gamm*pi/180.0D0)
            vecs(2,2)=bb*sin(gamm*pi/180.0D0)
            vecs(3,1)=cc*cos(beta*pi/180.0D0)
	    vecs(3,2)=(cc*bb*cos(alph*Pi/180.0D0)
     &                  -vecs(2,1)*vecs(3,1))/vecs(2,2)
            vecs(3,3)=sqrt(cc**2-vecs(3,1)**2-vecs(3,2)**2)
            nebvecs(n_rep+2,1:3,1:3)=vecs(1:3,1:3)
            read(51,1000,err=100) line
            do i=1,n_at
            ! read fractional coordinates and calculates absolute
            ! ones
                  read(51,'(A7,3(F13.6),A20)',err=100) line,
     &                atoms(i)%where(1:3),line
                  call frac2abs(atoms(i)%where,vecs,atoms(i)%abswhere)
                  nebcoords(n_rep+2,i,1:3)=atoms(i)%abswhere(1:3)
                  ! shift by cell vectors if atoms have moved outside
                  ! cell edges
                  do ivec=1,3
                    do while (absvec(nebcoords(n_rep+2,i,1:3)
     &                 +vecs(ivec,1:3)-nebcoords(n_rep+1,i,1:3)).lt.
     &                  absvec(nebcoords(n_rep+2,i,1:3)
     &                 -nebcoords(n_rep+1,i,1:3)))
                       nebcoords(n_rep+2,i,1:3)=nebcoords(n_rep+2,i,1:3)
     &                    +vecs(ivec,1:3)   
                    end do
                    do while (absvec(nebcoords(n_rep+2,i,1:3)
     &                 -vecs(ivec,1:3)-nebcoords(n_rep+1,i,1:3)).lt.
     &                  absvec(nebcoords(n_rep+2,i,1:3)
     &                 -nebcoords(n_rep+1,i,1:3)))
                       nebcoords(n_rep+2,i,1:3)=nebcoords(n_rep+2,i,1:3)
     &                    -vecs(ivec,1:3)
                    end do
                  end do
            end do
            reactcoord(n_rep+2)=matnorm(nebcoords(n_rep+2,1:n_at,1:3)
     &       -initcoords)
        end if
        goto 16

17      continue
      case('out','OUT','mbpp','MBPP')
          ! determine number of replicae (images)
          n_rep=0
          filename="./img0/"
!          inquire(directory=filename,exist=isdir)
          inquire(file=filename,exist=isdir)
          do while (isdir)
             n_rep=n_rep+1
             filename=" "
             write(filename,'("./img",I0,"/")') n_rep
!             inquire(directory=filename,exist=isdir)
             inquire(file=filename,exist=isdir)
          end do   
          n_rep=n_rep-2
          ! determine initial structure
          i_rep=1 
          filename=" "
          write(filename,'("img0/OUT")') 
          ! get fractional coordinates and cell vectors
          call read_coords(filename,'mbpp',atoms,n_at,species,
     &         nspecies,vecs)
          if(.not.allocated(reactcoord))
     &      allocate(reactcoord(n_rep+2),neben(n_rep+2))
          if(.not.allocated(nebcoords))
     &      allocate(nebcoords(n_rep+2,n_at,1:3),
     &      nebvecs(n_rep+2,1:3,1:3))
          if(.not.allocated(initcoords))
     &         allocate(initcoords(n_at,1:3))
          ! get species
          call getspecies(atoms,species)
            nspecies=size(species)
          nebvecs(i_rep,1:3,1:3)=vecs(1:3,1:3)
          ! get absolute coordinates
          do i=1,n_at
             call frac2abs(atoms(i)%where,vecs,atoms(i)%abswhere)
             nebcoords(i_rep,i,1:3)=atoms(i)%abswhere(1:3)
          end do
          ! get coordinates of initial image
          initcoords(1:n_at,1:3)=nebcoords(i_rep,1:n_at,1:3)
          ! get energy
          call getenergy(filename,'mbpp',neben(i_rep))
          emin=neben(i_rep)
          emax=neben(i_rep)
          nmax=0
          nmin=0
          ! get reaction coordinate
          reactcoord(i_rep+1)=matnorm(nebcoords(i_rep,1:n_at,1:3)
     &       -initcoords)
          ! determine structures of the rest of the replicae 
          do i_rep=2,n_rep+1
             filename=" "
             write(filename,'("img",I0,"/",A)') i_rep-1,
     &           trim(adjustl(infile))
             ! get fractional coordinates and cell vectors
             call read_coords(filename,'mbpp',atoms,n_at,species,
     &            nspecies,vecs)
             nebvecs(i_rep,1:3,1:3)=vecs(1:3,1:3)
             ! get absolute coordinates
             do i=1,n_at
                call frac2abs(atoms(i)%where,vecs,atoms(i)%abswhere)
                nebcoords(i_rep,i,1:3)=atoms(i)%abswhere(1:3)
                ! shift by cell vectors if atoms have moved outside
                ! cell edges
                do ivec=1,3
                  do while (absvec(nebcoords(i_rep,i,1:3)
     &               +vecs(ivec,1:3)-nebcoords(i_rep-1,i,1:3)).lt.
     &                absvec(nebcoords(i_rep,i,1:3)
     &               -nebcoords(i_rep-1,i,1:3)))
                     nebcoords(i_rep,i,1:3)=nebcoords(i_rep,i,1:3)
     &                  +vecs(ivec,1:3)   
                  end do
                  do while (absvec(nebcoords(i_rep,i,1:3)
     &               -vecs(ivec,1:3)-nebcoords(i_rep-1,i,1:3)).lt.
     &                absvec(nebcoords(i_rep,i,1:3)
     &               -nebcoords(i_rep-1,i,1:3)))
                     nebcoords(i_rep,i,1:3)=nebcoords(i_rep,i,1:3)
     &                  -vecs(ivec,1:3)
                  end do
                end do
             end do
             ! get energy
             call getenergy(filename,'mbpp',neben(i_rep))
             if (neben(i_rep).lt.emin) then
               emin=neben(i_rep)
               nmin=i_rep-1
             end if
             if (neben(i_rep).gt.emax) then
               emax=neben(i_rep)
               nmax=n_rep-1
             end if
             ! get reaction coordinate
             reactcoord(i_rep)=matnorm(nebcoords(i_rep,1:n_at,1:3)
     &         -initcoords)
          end do
          ! determine last structure
          i_rep=n_rep+2
          filename=" "
          write(filename,'("img",I0,"/OUT")') i_rep-1
          ! get fractional coordinates and cell vectors
          call read_coords(filename,'mbpp',atoms,n_at,species,
     &         nspecies,vecs)
          nebvecs(i_rep,1:3,1:3)=vecs(1:3,1:3)
          ! get absolute coordinates
          do i=1,n_at
             call frac2abs(atoms(i)%where,vecs,atoms(i)%abswhere)
             nebcoords(i_rep,i,1:3)=atoms(i)%abswhere(1:3)
             ! shift by cell vectors if atoms have moved outside
             ! cell edges
             do ivec=1,3
               do while (absvec(nebcoords(i_rep,i,1:3)
     &            +vecs(ivec,1:3)-nebcoords(i_rep-1,i,1:3)).lt.
     &             absvec(nebcoords(i_rep,i,1:3)
     &            -nebcoords(i_rep-1,i,1:3)))
                  nebcoords(i_rep,i,1:3)=nebcoords(i_rep,i,1:3)
     &               +vecs(ivec,1:3)   
               end do
               do while (absvec(nebcoords(i_rep,i,1:3)
     &            -vecs(ivec,1:3)-nebcoords(i_rep-1,i,1:3)).lt.
     &             absvec(nebcoords(i_rep,i,1:3)
     &            -nebcoords(i_rep-1,i,1:3)))
                  nebcoords(i_rep,i,1:3)=nebcoords(i_rep,i,1:3)
     &               -vecs(ivec,1:3)
               end do
             end do
          end do
          ! get energy
          call getenergy(filename,'mbpp',neben(i_rep))
          if (neben(i_rep).lt.emin) then
            emin=neben(i_rep)
            nmin=i_rep-1
          end if
          if (neben(i_rep).gt.emax) then
            emax=neben(i_rep)
            nmax=n_rep-1
          end if
          ! get reaction coordinate
          reactcoord(i_rep)=matnorm(nebcoords(i_rep,1:n_at,1:3)
     &       -initcoords)
      case default
      end select

      ! write structures to animated xsf readable by xcrysden
      open(52,file="NEBGEO.XSF",status='replace',err=110)
      write(52,*,err=110)  'ANIMSTEPS ',n_rep+2 
      do i_rep=1,n_rep+2
        write(52,*,err=110)  'CRYSTAL'
        write(52,*,err=110)  'PRIMVEC ',i_rep
        write(52,1500,err=110) nebvecs(i_rep,1,1:3)
        write(52,1500,err=110) nebvecs(i_rep,2,1:3)
        write(52,1500,err=110) nebvecs(i_rep,3,1:3)
        write(52,*,err=110) 'PRIMCOORD ',i_rep
        write(52,*,err=110) n_at,' 1'
        do i=1,n_at
          write(52,3000,err=110) atoms(i)%name,
     &     nebcoords(i_rep,i,1:3)
        end do
      end do
      close(52)
      
      ! write each image structure to a single xsf file 
      do i_rep=1,n_rep+2
         filename=" "
         write(filename,'("GEO.IMG",I0,".xsf")') i_rep-1
         ! get fractional coordinates of the current image
         do i=1,n_at
           call abs2frac(nebcoords(i_rep,i,1:3),
     &        nebvecs(i_rep,1:3,1:3),atoms(i)%where)
         end do
         call write_coords(filename,"xsf",atoms,n_at,
     &      species,nspecies,nebvecs(i_rep,1:3,1:3))
      end do

      ! write each image structure to separate xsf files for cores and
      ! shels,
      ! first count numbers of cores and shels
      ncores=0
      nshels=0
      do i=1,n_at
        if(atoms(i)%core(1:4).eq."core") ncores=ncores+1
      end do
      nshels=n_at-ncores
      ! create lists of core and shel atoms
      allocate(coreatoms(ncores),shelatoms(nshels))
      icore=0
      ishel=0
      do i=1,n_at
        if(atoms(i)%core(1:4)=="core") then
              icore=icore+1
              coreatoms(icore)=atoms(i)
        else
              ishel=ishel+1
              shelatoms(ishel)=atoms(i)
        end if
      end do
      ! write
      do i_rep=1,n_rep+2
         ! get fractional coordinates of the current image
         icore=0
         ishel=0
         do i=1,n_at
           call abs2frac(nebcoords(i_rep,i,1:3),
     &        nebvecs(i_rep,1:3,1:3),atoms(i)%where)
           if(atoms(i)%core(1:4)=="core") then
                 icore=icore+1
                 coreatoms(icore)=atoms(i)
           else
                 ishel=ishel+1
                 shelatoms(ishel)=atoms(i)
           end if
         end do
         filename=" "
         write(filename,'("GEO.IMG",I0,".c.xsf")') i_rep-1
         call write_coords(filename,"xsf",coreatoms,ncores,
     &      species,nspecies,nebvecs(i_rep,1:3,1:3))
         filename=" "
         write(filename,'("GEO.IMG",I0,".s.xsf")') i_rep-1
         call write_coords(filename,"xsf",shelatoms,nshels,
     &      species,nspecies,nebvecs(i_rep,1:3,1:3))
      end do


      ! write energies
      open(53,file='NEBENERGIES.DAT',status='replace')
      do i_rep=1,n_rep+2
        write(53,1700,err=110) i_rep-1,neben(i_rep),
     &     reactcoord(i_rep)/dble(n_at) 
      end do
      write(53,*) '# image with lowest energy:',nmin
      write(53,*) '# image with highest energy:',nmax
      write(53,'(A30,F9.6,A3)') '# energy difference:',
     &      emax-emin,' eV'
      close(53)

!     end normally
      if(isopen12) then
            write(12,fsubendext) "nebana"
      else
            print fsubendext, "nebana"
      end if
      return

! errors
! error when opening GULP NEB output file for reading
100   nerr=nerr+1
      close(51)
      if (isopen12) then
         write(12,ferrmssg) 
     &    "Something wrong with NEB output file"
      else
         print ferrmssg,
     &    "Something wrong with NEB output file"
      end if
      return
! error when opening output file for writing structures
110   nerr=nerr+1
      close(51)
      close(52)
      close(53)
      if (isopen12) then
         write(12,ferrmssg) 
     &    "Error during writing"
      else
         print ferrmssg,
     &    "Error during writing"
      end if
      return

! formats
 1000 format(A120)     
 1500 format(3(F12.6))  
 1600 format(I3,F15.6)     
 1700 format(I3,F15.6,F18.12)    
 2000 format(6(F12.6))  
 3000 format(A3,3(F12.6))  

      end subroutine     

c---------------------------------------------------------------------

      subroutine nebchain(infile,informat,outformat,nimg)
      ! reads a chain of structure files and sets up an input file for a NEB calculation.
      use defs
      use misc
      use readcoords
      use writecoords
      implicit none
      
      ! global variables
      character(len=*) infile,informat,outformat
      integer nimg

      ! local variables
      type(atom),allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer iatom,iimg,natoms,nspecies 
      double precision aa,bb,cc,alph,beta,gamm,xxf,yyf,zzf,xx,yy,zz
      double precision vecs(1:3,1:3),cellpars(1:6)
      double precision, allocatable :: initcoords(:,:)
      double precision, allocatable :: nebcoords(:,:,:),nebvecs(:,:,:)
      character filename*40,line*120,core*4,outfile*40
      logical isopen12

      talk=.true.
      ! write hello
      inquire(unit=12,opened=isopen12)
      if(talk) then
        if(isopen12) then
            write(12,fsubstart) "nebchain"
        else
            print fsubstart, "nebchain"
        end if
      end if
 
      select case (outformat)
      case('gulp','GULP','gin','GIN','res','RES')
            ! read first structure 
            filename=adjustl(infile)
            write(filename(len_trim(filename)+1:len(filename)),
     &         '("0.",A)') adjustl(trim(informat))   
            print '(8x,"reading file",A)',filename
            call read_coords(filename,informat,atoms,natoms,species,
     &                       nspecies,vecs)
            ! write first structure to output file
            outfile="CHAIN.GIN"
            call write_coords(outfile,outformat,atoms,natoms,species,
     &                    nspecies,vecs)
            open(53,file=outfile,status="old",position="append")
            write(53,'("#",/,"#  Nudged elastic band data",/,"#",/,
     &"nebspring      0.000050")')
            ! write structures of intermediate images
            do iimg=1,nimg-1   
                  ! read image structure 
                  filename=adjustl(infile)
                  write(filename(len_trim(filename)+1:len(filename)),
     &             '(I0,".",A)') iimg,adjustl(informat)   
                  print '(8x,"reading file",A)',filename
                  call read_coords(filename,informat,atoms,natoms,
     &                 species,nspecies,vecs)
                  write(53,'("rcell    ",I0)') iimg
                  call vecs2cellpars(vecs,cellpars)
                  write(53,'(6(F10.6))') cellpars 
                  write(53,'("rfractional    ",I0)') iimg
                  do iatom=1,natoms
                    write(53,'(4(F10.6))')
     &               atoms(iatom)%where,atoms(iatom)%charge
                  end do
            end do
            ! read last structure 
            filename=adjustl(trim(infile))
            write(filename(len_trim(filename)+1:len(filename)),
     &       '(I0,".",A)') nimg,adjustl(trim(informat))   
            print '(8x,"reading file",A)',filename
            call read_coords(filename,informat,atoms,natoms,species,
     &                       nspecies,vecs)
            ! write last structure to output file
            write(53,'("fcell    ")') 
            call vecs2cellpars(vecs,cellpars)
            write(53,'(6(F10.6))') cellpars 
            write(53,'("ffractional    ")') 
            do iatom=1,natoms
              write(53,'(4(F10.6))')
     &         atoms(iatom)%where,atoms(iatom)%charge
            end do
            write(53,'("nebreplica    ",I0,/,"#")') nimg-1
            close(53)
      case default
            nerr=nerr+1
            if(isopen12) then
              write(12,ferrmssg) "unknown output format"
            else
              print ferrmssg, "unknown output format"
            end if
            return
      end select

!     end normally
      if(isopen12) then
            write(12,fsubendext) "nebchain"
      else
            print fsubendext, "nebchain"
      end if
      return

      end subroutine     

!--------------------------------------------------------------------------   

      subroutine create_neb_chain(file1,file2,format1,vecs1,vecs2,nimg)
      use defs
      use readcoords
      use writecoords
      implicit none
      type(atom), allocatable :: atoms1(:),atoms2(:),atoms(:)
      type(element), allocatable :: species1(:),species2(:)
      double precision vecs1(3,3),vecs2(3,3),vecs(3,3)
      integer, intent(in) :: nimg
      integer natoms1,natoms2,nspecies1,nspecies2,img,iatom
      character(len=*), intent(in) :: file1,file2,format1
      character*1024 filename

      print fsubstart,'create_neb_chain'
      print*,file1,file2,format1
      call read_coords(file1,format1,atoms1,natoms1,species1,             &
     &                 nspecies1, vecs1)
      call read_coords(file2,format1,atoms2,natoms2,species2,             &
     &                 nspecies2, vecs2)
      !
      ! begin sanity check
      !
      if (natoms1.ne.natoms2) call error_stop('files contain different n  &
     &umbers of atoms')
      if (nspecies1.ne.nspecies2) call error_stop('files contain differe  &
     &nt numbers of species')
      if(any(species1(:)%name.ne.species2(:)%name)) call error_stop('fil  &
     &es contain different species, or species are ordered differently')
      if (any(atoms1(:)%name.ne.atoms2(:)%name)) call error_stop('files   &
     &contain different atoms, or atoms are ordered differently')
      !
      ! end sanity check
      !
      allocate(atoms(natoms1))
      do img=2,nimg-1
        filename=trim(adjustl(trim(file1)))//"_IMG_" 
        write(filename(len_trim(filename)+1:),'(I0)') img
        !print*,filename  
        atoms=atoms1
        do iatom=1,natoms1
          atoms(iatom)%where(1:3)=(1.0d0-dble(img-1)/dble(nimg-1))        &
     &                 *atoms1(iatom)%where(1:3)                          &
     &                  +dble(img-1)/dble(nimg-1)                         &
     &                 *atoms2(iatom)%where(1:3)
          atoms(iatom)%abswhere=(1.0d0-dble(img-1)/dble(nimg-1))          &
     &   *atoms1(iatom)%abswhere+dble(img-1)/dble(nimg-1)                 &
     &   *atoms2(iatom)%abswhere
        end do   
        vecs=(1.0d0-dble(img-1)/dble(nimg-1))*vecs1                       &
     &       +dble(img-1)/dble(nimg-1)*vecs2
        call write_coords(filename,format1,atoms,natoms1,species1,        &
     &                    nspecies1,vecs)
      end do ! img
      deallocate(atoms)
      !
      print fsubendext,'create_neb_chain'
      !
      end subroutine create_neb_chain

c---------------------------------------------------------------------

      end module
