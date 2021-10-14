      module convert
      implicit none

      contains

c---------------------------------------------------------------------

      subroutine conv(file1,file2,format1,
     &           format2,vecs)

      use defs
      use readcoords
      use writecoords
      implicit none

      ! variables
      character(len=*) :: file1
      character(len=*) :: file2
      character(len=*) :: format1       ! format of input file
      character(len=*) :: format2          ! format of outputfile

      double precision :: matrix(1:3,1:3)
      double precision vecs(1:3,1:3)
      type(atom), allocatable :: atoms(:)
      type(element), allocatable :: species(:)
      integer :: natoms,nspecies

      ! local variables
      !integer iatom,j,k,n,p,q,r,s,nuc(1:3),iallatom,nallatoms
      !integer ispecies
      !logical,allocatable :: delete(:)
      !type(atom), allocatable :: allatoms(:)
      
      talk=.true.
      if(talk)print fsubstart,"convert"
      
      ! read in  coordinates ...
      if(talk) then
          print '(" ")'
          print '(8x,"calling read_coords...")'
      end if
      call read_coords(file1,format1,atoms,natoms,species,
     &           nspecies,vecs)
      if(talk) print '(8x,"Coords read in.")'
      ! write coordinates
      call write_coords(file2,format2,atoms,natoms,
     &                  species,nspecies,vecs)

      if(talk)print fsubendext, 'convert'

      end subroutine

c---------------------------------------------------------------------

      end module
