        module readinput
        implicit none

        contains

c---------------------------------------------------------------------

        subroutine read_10(strucfile,strucformat)
        use defs
        implicit none

        ! global
        character(len=*) strucfile,strucformat
        ! local
        character line*150
        logical strucfileread, strucformatread

        strucfileread=.false.
        strucformatread=.false.
        write(12,fsubstart) trim('read_10')
10      read(10,'(A150)',end=1000,err=1000) line  
        line=adjustl(line)
        if (line(1:1).eq.'#') goto 10
        if (line(1:1).eq.'!') goto 10
        if (index(line,'#').ge.1) line(index(line,'#'):len(line))=' '
        if (index(line,'!').ge.1) line(index(line,'!'):len(line))=' '
        if (line(1:2).eq."10") goto 20
        goto 10

20      read(10,'(A150)',end=1000) line  
        line=adjustl(line)
        if (line(1:1).eq.'#') goto 20
        if (line(1:1).eq.'!') goto 20
        if (index(line,'#').ge.1) line(index(line,'#'):len(line))=' '
        if (index(line,'!').ge.1) line(index(line,'!'):len(line))=' '
        if (index(line,"filename").ge.1) then
            read(line(index(line,"filename")+8:len(line)),*,err=1000) 
     &         strucfile
            strucfileread=.true.
        end if
        if (index(line,"format").ge.1) then
            read(line(index(line,"format")+6:len(line)),*,err=1000)
     &         strucformat
            strucformatread=.true.
        end if
        if(strucfileread.and.strucformatread) goto 30
        goto 20

30      write(12,fsubendext) 'read_10'
        return

1000    nerr=nerr+1
        !close(10)
        write(12,ferrmssg) "something wrong with structure section (10)"
        !close(12)
        return    

        end subroutine

c---------------------------------------------------------------------
        
        subroutine read_15(nsymm,symops,isymops,vecs)
        use defs
        use misc
        use symmetry
        implicit none

        ! global
        type(symop), allocatable :: symops(:),isymops(:)
        double precision, intent(in) :: vecs(1:3,1:3)
        integer nsymm
        ! local
        character line*150
        integer i,j
        double precision trala(1:3)

        write(12,fsubstart) trim('read_15')
        rewind(10)
10      read(10,'(A150)',end=1000,err=1000) line  
        line=adjustl(line)
        if (line(1:1).eq.'#') goto 10
        if (line(1:1).eq.'!') goto 10
        if (index(line,'#').gt.1) line(index(line,'#'):len(line))=' '
        if (index(line,'!').gt.1) line(index(line,'!'):len(line))=' '
        if (line(1:2).eq."15") goto 20
        goto 10

20      read(10,'(A150)',end=1000) line  
        line=adjustl(line)
        if (line(1:1).eq.'#') goto 20
        if (line(1:1).eq.'!') goto 20
        if (index(line,'#').gt.1) line(index(line,'#'):len(line))=' '
        if (index(line,'!').gt.1) line(index(line,'!'):len(line))=' '
        ! read nsymm
        if (index(line,"nsymm").ge.1) then
            read(line(index(line,"nsymm")+5:len(line)),*,err=1000) 
     &         nsymm
          allocate(symops(nsymm),isymops(nsymm))
          do i=1,nsymm
              read(10,'(A150)',err=1000,end=1000)line
              line=adjustl(line)
              do while (line(1:1).eq.'!'.or.line(1:1).eq.'#') 
                read(10,'(A150)',err=1000,end=1000)line
                line=adjustl(line)
              end do
              !read(line(1:len(line)),*,err=1000,end=1000) 
              read(line,*,err=1000,end=1000) 
     &        symops(i)%mat(1,1:3),
     &        symops(i)%mat(2,1:3),
     &        symops(i)%mat(3,1:3),
     &        trala(1:3) 
              symops(i)%trala(1:3)=trala(1)*vecs(1,1:3)
     &           +trala(2)*vecs(2,1:3)+trala(3)*vecs(3,1:3)
              call invert_matrix(symops(i)%mat,isymops(i)%mat) 
              isymops(i)%trala=
     &             -matmul(isymops(i)%mat,symops(i)%trala)
          end do
          ! check if symmetry operations form a group
          call checksymgroup(symops,isymops,nsymm)    
          goto 30
        end if
        goto 20

30      write (12,'(8x,"Symmetry operations: ")')
        write(12,'(8x,"Number of symmetry operations: ",I3)') nsymm
        write(12,'(8x,"Symmetry operations (matrix, translation:")')
        do i=1,nsymm
          write(12,'(8x,I5,9(F7.2),3x,3(F7.2))') i,symops(i)
        end do
        write(12,fsubendext) 'read_15'
        return

1000    nerr=nerr+1
        !close(10)
        write(12,ferrmssg) "something wrong with symmetry section (15)"
        !close(12)
        return    

        end subroutine

c---------------------------------------------------------------------
        
        subroutine read_20(emax,gmax,kstyle,kdens,kshift)
        use defs
        implicit none

        ! global
        double precision emax,gmax
        character kstyle*2
        integer kdens(1:3)
        double precision kshift(1:3) 
        ! local
        character line*150
        logical emaxread
        logical gmaxread
        logical kstyleread
        logical kdensread
        logical kshiftread

        emaxread=.false.
        gmaxread=.false.
        kstyleread=.false.
        kdensread=.false.
        kshiftread=.false.
        
        write(12,fsubstart) trim('read_20')

        ! find the line beginning with "20"
10      read(10,'(A150)',end=1000,err=1000) line  
        line=adjustl(line)
        if (line(1:1).eq.'#') goto 10
        if (line(1:1).eq.'!') goto 10
        if (index(line,'#').ge.1) line(index(line,'#'):len(line))=' '
        if (index(line,'!').ge.1) line(index(line,'!'):len(line))=' '
        if (line(1:2).eq."20") goto 20
        goto 10

20      read(10,'(A150)',end=1000) line  
        line=adjustl(line)
        if (line(1:1).eq.'#') goto 20
        if (line(1:1).eq.'!') goto 20
        if (index(line,'#').ge.1) line(index(line,'#'):len(line))=' '
        if (index(line,'!').ge.1) line(index(line,'!'):len(line))=' '
        ! read emax
        if (index(line,"emax").ge.1) then
            read(line(index(line,"emax")+4:len(line)),*,err=1000) 
     &         emax
            emaxread=.true.
            ! determine unit (default=Ryd):
            if(index(line,"eV").gt.0.or.index(line,"ev").gt.0) 
     &         emax=emax/Ryd
        end if
        ! read gmax
        if (index(line,"gmax").ge.1) then
            read(line(index(line,"gmax")+4:len(line)),*,err=1000)
     &         gmax
            gmaxread=.true.
            ! determine unit (default=1/bohr):
            if(index(line,"Angs").gt.0.or.index(line,"angs").gt.0) 
     &         gmax=gmax*bohr
        end if
        ! read kmesh parameters
        if (index(line,"kmesh").ge.1) then
            read(line(index(line,"style")+5:len(line)),*,err=1000)
     &         kstyle
            kstyleread=.true.
            read(line(index(line,"dens")+4:len(line)),*,err=1000)
     &         kdens(1:3)
            kdensread=.true.
            read(line(index(line,"shift")+5:len(line)),*,err=1000)
     &         kshift
            kshiftread=.true.
        end if
        if(emaxread.and.gmaxread.and.kstyleread.and.kdensread
     &      .and.kshiftread)
     &      goto 30
        goto 20

30      write (12,'(8x,"Basis data: ")')
        write(12,'(8x,"emax: ",F7.2," eV")') emax*Ryd
        write(12,'(8x,"gmax: ",F7.2," (1/bohr)")') gmax
        write(12,'(8x,"gmax: ",F7.2," (1/A)")') gmax/bohr
        write(12,'(8x,"kmesh style: ",A3)') trim(kstyle)
        write(12,'(8x,"# of k points along a, b, c: ",I4,I4,I4)') kdens
        write(12,'(8x,"shift of kmesh: ",F7.4)') kshift
        write(12,fsubendext) 'read_20'
        return

1000    nerr=nerr+1
        !close(10)
        write(12,ferrmssg) "something wrong with basis section (20)"
        !close(12)
        return    

        end subroutine

c---------------------------------------------------------------------
        
        subroutine read_25(species,pspots)
        use defs
        implicit none

        ! global
        type(element) :: species(:)
        type(pspot) pspots(:),pspotdum
        ! local
        character line*150
        character pseudostyle*256,pseudofile*256
        integer nspecies,ispecies,jspecies

        write(12,fsubstart) trim('read_25')
        write(12,'(8x,"Reading pseudopotential information")') 
        
        nspecies=size(species)
        pspots%name=' '
        ispecies=1
        !rewiind(10)
        ! search for the number 25 which announces the start of the
        ! pseudopotential information
10      read(10,'(A150)',end=1000,err=1000) line  
        line=adjustl(line)
        if (line(1:1).eq.'#') goto 10
        if (line(1:1).eq.'!') goto 10
        if (index(line,'#').gt.1) line(index(line,'#'):len(line))=' '
        if (index(line,'!').gt.1) line(index(line,'!'):len(line))=' '
        if (line(1:2).eq."25") goto 20
        goto 10

        ! read the formats and filenames of all pseudopotential files,
        ! then read the pseudopotentials from those files and store them
        ! in the variable "pspots"
20      read(10,'(A150)',end=1000) line  
        line=adjustl(line)
        if (line(1:1).eq.'#') goto 20
        if (line(1:1).eq.'!') goto 20
        if (index(line,'#').gt.1) line(index(line,'#'):len(line))=' '
        if (index(line,'!').gt.1) line(index(line,'!'):len(line))=' '
        if (index(line,"pseudo").ge.1) then
          read(line(1:2),*,err=1000) 
     &      pspots(ispecies)%name(1:2)
          read(line(index(line,"style")+5:len(line)),*,err=1000) 
     &      pseudostyle
          read(line(index(line,"file")+4:len(line)),*,err=1000) 
     &         pseudofile
          write (12,'(8x,A2," Pseudopotential style, file: ")')
     &         pspots(ispecies)%name
          write(12,*) '        ',trim(adjustl(pseudofile))
          select case(pseudostyle)
          case('fhi','FHI')
                  write(12,'(8x,"psp format: fhi")')
                  call read_ps_abinit_fhi(pseudofile,pspots(ispecies) )
                  if(nerr.gt.0) goto 1002
          case('mbpp','MBPP')
                  write(12,'(8x,"psp format: mbpp")')
                  call read_ps_mbpp(pseudofile,pspots(ispecies))
                  if(nerr.gt.0) goto 1002
          case default
                  goto 1001
          end select
          ispecies=ispecies+1
          if(ispecies==nspecies+1) goto 30
        end if
        goto 20

30      write (12,'(8x,"Pseudopotentials have been read.")') 
        ! sort pseudopotentials such that they have the same order as
        ! the elements in the "species" variable
        do ispecies=1,nspecies
          if (pspots(ispecies)%name(1:2).ne.species(ispecies)%name(1:2))
     &    then     
            do jspecies=1,nspecies
              if(pspots(jspecies)%name(1:2)
     &        .eq.species(ispecies)%name(1:2)) then
                pspotdum=pspots(jspecies)
                pspots(jspecies)=pspots(ispecies)
                pspots(ispecies)=pspotdum
              end if
            end do
          end if
        end do  
        write(12,fsubendext) 'read_25'
        return

1000    nerr=nerr+1
        close(10)
        write(12,ferrmssg) "something wrong with pseudo section (25)"
        !close(12)
        return    
1001    nerr=nerr+1
        close(10)
        write(12,ferrmssg) "Unknown pseudopotential format"
        return 
1002    nerr=nerr+1
        close(10)
        write(12,ferrmssg) "Could not read pseudopotentials."
        call error_stop('Could not read pseudopotential(s)')
        return 

        end subroutine read_25

c---------------------------------------------------------------------

        subroutine read_ps_abinit_fhi(pseudofile,pspots_i)
        ! reads pseudopotentials generated with fhi98pp as linked on
        ! abinit homepage
        use defs
        !use misc
        use integr
        implicit none
        type(pspot) pspots_i
        character (len=*) pseudofile
        integer l,nr,iline,ir,lmax_i,iorb
        integer idum
        double precision rdum,dens_i,occsum
        character line*1024    

        write(12,fsubstart) 'read_ps_fhi'
        open(51,file=pseudofile,status='old',form='formatted',
     &       err=1000)
        read(51,'(A1024)') line
        read(51,*) pspots_i%znuc,pspots_i%zion
        read(51,*) idum,pspots_i%ixc,lmax_i,pspots_i%lloc
        pspots_i%lmaxps=lmax_i
        do iline=1,15
          read(51,'(A1024)') line
        end do
        read(51,*) nr
        allocate(pspots_i%r(1:nr),pspots_i%wf(0:lmax_i,1:nr),
     &     pspots_i%vps(0:lmax_i,1:nr),pspots_i%pvdens(1:nr),
     &     pspots_i%pcdens(1:nr),pspots_i%occ(0:lmax_i))
        pspots_i%pvdens=0.0D0
        pspots_i%pcdens=0.0D0
        ! get atomic pseudo orbital occupations :
        pspots_i%occ=0.0D0
        ! 1. find the highest valence orbital
        iorb=0
        occsum=0.0D0
        do while (occsum.lt.pspots_i%znuc)
          iorb=iorb+1
          occsum=occsum+occ_list(iorb)
        end do
        write(12,'(8x,"Highest occupied valence orbital (n,l):",I2,I2)')
     &    n_list(iorb),l_list(iorb)  
        ! 2. go downwards and fill orbitals with valence electrons until
        ! all valence electrons (=zion) are placed.
        write(12,'(8x,"Occupations of valence orbitals (n,l,occ):")')
        occsum=0.0D0
        do while (occsum.lt.pspots_i%zion)
          l=l_list(iorb)
          if(l.gt.lmax_i) goto 1100
          pspots_i%occ(l)=min(occ_list(iorb),pspots_i%zion-occsum)
          occsum=occsum+pspots_i%occ(l)
          write(12,'(8x,I2,I2,F10.4)') n_list(iorb),l_list(iorb),
     &      pspots_i%occ(l)
          iorb=iorb-1  
        end do 
        ! read pseudo wavefunctions and potentials
        do l=0,lmax_i-1
          do ir=1,nr
            read(51,*) idum,pspots_i%r(ir),pspots_i%wf(l,ir),
     &        pspots_i%vps(l,ir)
          end do
          read(51,'(A1024)') line
        end do
        do ir=1,nr
          read(51,*) idum,pspots_i%r(ir),pspots_i%wf(lmax_i,ir),
     &    pspots_i%vps(lmax_i,ir)
        end do
        ! read partial core density, if present
        do ir=1,nr
          read(51,*,end=100) pspots_i%r(ir),pspots_i%pcdens(ir)
        end do
100     close(51) 
        ! calculate integrals over wavefunctions
!        do l=0,lmax_i
!          !dens_i=4.0D0*Pi*int_radial(pspots_i%r,pspots_i%wf(l,1:nr)**2)
!          dens_i=int_1d(pspots_i%r,pspots_i%wf(l,1:nr)**2)
!          write(12,'(8x,"integrated density of PSWF for l=",I1,":",
!     &F10.4)') l,dens_i
!        end do
        ! change units from Hartree to Ryd
        pspots_i%vps=pspots_i%vps*2.0d0
!        pspots_i%pcdens=pspots_i%pcdens*4.0D0*Pi/bohr**3
        pspots_i%pcdens=pspots_i%pcdens*4.0D0*Pi*pspots_i%r**2 !in the fhi file it lacks the factor r**2
        ! get atomic pseudo density. density=sum over (occupation * wf^2). The WF include
        ! already the factor r**2 that occurs in the radial integral.
        do l=0,lmax_i
          do ir=1,nr
            pspots_i%pvdens(ir)=pspots_i%pvdens(ir)+pspots_i%occ(l)
     *      *pspots_i%wf(l,ir)**2
          end do
        end do
        ! Now all WF's and densities already include the factor r**2, so
        ! that 1=int_0^00 WF(r)**2 dr, 1=int_0^00 dens(r) dr
        ! calculate integrals over wavefunctions**2
        do l=0,lmax_i
          !dens_i=4.0D0*Pi*int_radial(pspots_i%r,pspots_i%wf(l,1:nr)**2)
          !dens_i=int_1d(pspots_i%r,pspots_i%wf(l,1:nr)**2)
          dens_i=int_1d_taylorn(pspots_i%r,pspots_i%wf(l,1:nr)**2,6)
          write(12,'(8x,"integrated density of PSWF for l=",I1,":",
     &F10.4)') l,dens_i
        end do
        ! calculate integral over valence density
        ! try different integration schemes
        !write(12,'(8x,"integration scheme Riemann:")')
        !dens_i=int_1d(pspots_i%r,pspots_i%pvdens(1:nr))
        !write(12,'(8x,"Total integrated valence density:",F15.9)')dens_i
        !write(12,'(8x,"integration scheme Taylor 2:")')
        !dens_i=int_1d_taylorn(pspots_i%r,pspots_i%pvdens,2)
        !write(12,'(8x,"Total integrated valence density:",F15.9)')dens_i
        !write(12,'(8x,"integration scheme Taylor 6:")')
        dens_i=int_1d_taylorn(pspots_i%r,pspots_i%pvdens,6)
        write(12,'(8x,"Total integrated valence density:",F15.9)')dens_i
        ! calculate integral over partial core density
        !write(12,'(8x,"integration scheme Riemann:")')
        !dens_i=int_radial(pspots_i%r,pspots_i%pcdens(1:nr))
        !write(12,'(8x,"Integrated partial core density:",F15.9)')dens_i
        write(12,'(8x,"radial integration scheme Taylor 6:")')
        dens_i=int_1d_taylorn(pspots_i%r,pspots_i%pcdens(1:nr),6)
        write(12,'(8x,"Integrated partial core density:",F15.9)')dens_i
        !write(12,'(8x,"integration scheme Taylor 2:")')
        !dens_i=int_1d_taylorn(pspots_i%r,pspots_i%pcdens*pspots_i%r**2,
!     &            2)
!        write(12,'(8x,"Integrated partial core density:",F15.9)')dens_i
!        write(12,'(8x,"integration scheme Taylor 4:")')
!        dens_i=int_1d_taylorn(pspots_i%r,pspots_i%pcdens*pspots_i%r**2,
!     &            4)
!        write(12,'(8x,"Integrated partial core density:",F15.9)')dens_i
!        write(12,'(8x,"1d integration scheme Taylor 6:")')
!        dens_i=int_1d_taylorn(pspots_i%r,pspots_i%pcdens*pspots_i%r**2,
!     &            6)
!        write(12,'(8x,"Integrated partial core density:",F15.9)')dens_i
        write(12,fsubend) 
        return

1000    nerr=nerr+1
        close(51)
        write(12,ferrmssg) 'Pseudopotential could not be read' 
        call error_stop('Pseudopotential could not be read')
        return       
1100    nerr=nerr+1
        write(12,ferrmssg) 'Occupations could not be determined' 
        close(51)
        return       

        end subroutine

c---------------------------------------------------------------------

        subroutine read_ps_mbpp(pseudofile,pspots_i)
        ! reads pseudopotentials generated with Berkeley atomic program (to be read by MBPP)
        use defs
        use integr
        implicit none
        type(pspot) pspots_i
        character (len=*) pseudofile
        ! local variables
        integer l,nr,ir,lmax_i,iorb
        integer idum
        double precision rdum,dens_i,occsum
        character line*1024 
        
        character nameat*2,irel*3,nicore*4,icorr*2
        character*10 iray(1:6),ititle(1:7) 
        integer npotd,npotu,i
        double precision zval,a,b    

        write(12,fsubstart) 'read_ps_mbpp'
        ! read pseudopotential
        open(51,file=pseudofile,status='old',form='unformatted',
     &       err=1000)
        read(51) nameat,zval,icorr,irel,nicore,
     &    (iray(i),i=1,6),ititle,npotd,npotu,a,b,nr
        write(12,'(8x,"PSP file header:")') 
        write(12,'(8x,"element:",A2,", ionic charge:",F10.5,
     &    ", correlation:",A2,", relativist.:",A3,", partial core:",
     &    A4,",")')
     &    nameat,zval,icorr,irel,nicore
        write(12,'(8x,"title: ",6(A10))')
     &    (iray(i),i=1,6)
        write(12,'(8x,"pseudo occupations and cutoffs: ",7(A10))')ititle
        write(12,'(8x,"number of PSPs for down spin:",I2,
     &    ", for up-spin:",I2,", a=",F10.6,", b=",F10.6,
     &    ", number of radial gridpoints:",I6)') npotd,npotu,a,b,nr
        pspots_i%znuc=0.0D0
        do i=1,size(elements)
          if(nameat.eq.elements(i)) then
                  pspots_i%znuc=dble(i)
          end if
        end do
        ! determine ixc
        pspots_i%ixc=0
        if(icorr.eq."ca") pspots_i%ixc=7 ! Ceperley-Alder 
        if(icorr.eq."pb") pspots_i%ixc=11 ! Perdew-Burke-Ernzerhof 
        ! determine lmax
        lmax_i=npotd-1 
        pspots_i%lmaxps=lmax_i 
        ! like in MBPP, assume that lloc is 0
        pspots_i%lloc=0
        !determine zion:
        pspots_i%zion=zval
        allocate(pspots_i%r(1:nr),pspots_i%wf(0:lmax_i,1:nr),
     &     pspots_i%vps(0:lmax_i,1:nr),pspots_i%pvdens(1:nr),
     &     pspots_i%pcdens(1:nr),pspots_i%occ(0:lmax_i))
        pspots_i%pvdens=0.0D0
        pspots_i%pcdens=0.0D0
        pspots_i%wf(0:lmax_i,1:nr)=0.0D0
        ! get atomic pseudo orbital occupations :
        pspots_i%occ=0.0D0
        ! 1. find the highest valence orbital
        iorb=0
        occsum=0.0D0
        do while (occsum.lt.pspots_i%znuc)
          iorb=iorb+1
          occsum=occsum+occ_list(iorb)
        end do
        write(12,'(8x,"Highest occupied valence orbital (n,l):",I2,I2)')
     &    n_list(iorb),l_list(iorb)  
        ! 2. go downwards and fill orbitals with valence electrons until
        ! all valence electrons (=zion) are placed.
        write(12,'(8x,"Occupations of valence orbitals (n,l,occ):")')
        occsum=0.0D0
        do while (occsum.lt.pspots_i%zion)
          l=l_list(iorb)
          if(l.gt.lmax_i) goto 1100
          pspots_i%occ(l)=min(occ_list(iorb),pspots_i%zion-occsum)
          occsum=occsum+pspots_i%occ(l)
          write(12,'(8x,I2,I2,F10.4)') n_list(iorb),l_list(iorb),
     &      pspots_i%occ(l)
          iorb=iorb-1  
        end do 
        ! read pseudopotentials
        read(51) (pspots_i%r(ir),ir=1,nr)
        do l=0,npotd-1 !lmax_i-1
          read(51) idum,(pspots_i%vps(l,ir),ir=1,nr)
          if(idum.ne.l) then
                  write(12,ferrmssg) 'unexpected l in PP file'
                  nerr=nerr+1
                  close(51)
                  close(10)
                  close(12)
                  return
          end if
          ! pspots contain factor r, will be removed in the following.
          ! Also change from Ry to Hartree:
          do ir=2,nr
            pspots_i%vps(l,ir)=pspots_i%vps(l,ir)/(2.0d0*pspots_i%r(ir))
          end do
          if (pspots_i%r(1).gt.0.0d0) then
                  pspots_i%vps(l,1)=pspots_i%vps(l,1)/pspots_i%r(1)
          else
                  pspots_i%vps(l,1)=pspots_i%vps(l,2)
          end if
        end do
        ! there may be spin-orbit data in the file, ignore:
        do i=1,npotu
          read(51)
        end do 
        ! read core and valence density 
        read(51) (pspots_i%pcdens(ir),ir=1,nr)
        read(51) (pspots_i%pvdens(ir),ir=1,nr)
100     close(51) 
        ! calculate integral over valence and partial core density
        dens_i=int_1d_taylorn(pspots_i%r,pspots_i%pvdens,6)
        write(12,'(8x,"Total integrated valence density:",F15.9)')dens_i
        dens_i=int_1d_taylorn(pspots_i%r,pspots_i%pcdens(1:nr),6)
        write(12,'(8x,"Integrated partial core density:",F15.9)')dens_i
        write(12,fsubend) 
        return

1000    nerr=nerr+1
        close(51)
        write(12,ferrmssg) 'Pseudopotential could not be read' 
        call error_stop('Pseudopotential could not be read')
        return       
1100    nerr=nerr+1
        write(12,ferrmssg) 'Occupations could not be determined' 
        close(51)
        return       

        end subroutine

        end module
