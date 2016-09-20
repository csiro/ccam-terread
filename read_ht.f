! Conformal Cubic Atmospheric Model
    
! Copyright 2015 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
! This file is part of the Conformal Cubic Atmospheric Model (CCAM)
!
! CCAM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! CCAM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with CCAM.  If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------
      
      subroutine read_ht ( lat, ht, filepath10km )

! read topography data for specified latitude

      implicit none

! parameters
      integer      lui ! logical unit for input
      integer      nlong ! number of longitudes
      parameter (lui = 99)
      parameter (nlong = 4320) ! 360 * 12

! arguments
      integer      ht(nlong) ! topographic height in metres (output)
      real         lat ! latitude (input)
      character(len=*), intent(in) :: filepath10km

! variables
      character    buf(2,nlong) ! buffer
      integer      i ! subscript
      integer      ios ! i/o status
      integer(kind=8) ival ! int. value
      integer      j ! subscript
      logical      new ! true 1st call
      save new
      data new / .true. /

! open file
      if (new) then
         write(6,*)"open file='topo2' lui=",lui
         open (lui, file=trim(filepath10km)//'/'//'topo2',
     &   status='old', access='direct',
     &   form='unformatted', recl=nlong*2, iostat=ios)
         new = .false.
      end if

      i = nint((90.0 - lat) * 12.0 + 1.0)
      if ( i.lt.3 ) print *,'i=',i

      read (lui, rec=i, iostat=ios) buf

      do j = 1, nlong
        ival = 256 * ichar(buf(1,j)) + ichar(buf(2,j))
        if (ival .ge. 32768) then ! sign extend
           ival = ival - 65536
           !ival = ior( ival, X'FFFFFFFFFFFF0000' )
        end if
        ht(j) = ival
      end do

      end
