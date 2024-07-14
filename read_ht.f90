! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2024 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
      
subroutine read_ht( lat, ht, slope, ht_read, filepath10km )

use netcdf_m
      
! read topography data for specified latitude

implicit none

! parameters
integer, parameter :: lui=99 ! logical unit for input
integer, parameter :: nlong=4320 ! number of longitudes

! arguments
integer, dimension(nlong), intent(out) :: ht ! topographic height in metres (output)
real, dimension(nlong), intent(out) :: slope
real, dimension(nlong,3), intent(inout) :: ht_read
real, intent(in) :: lat ! latitude (input)
character(len=*), intent(in) :: filepath10km

! variables
character    buf(2,nlong) ! buffer
integer      i ! subscript
integer      ios ! i/o status
integer(kind=8) ival ! int. value
integer      j ! subscript
integer ierr
integer, save :: ncid, varid
real sx_e, sx_w, sy_n, sy_s, dx, dy
logical, save :: new=.true. ! true 1st call
logical, save :: ncfile=.false. ! is input a netcdf file?

! open file
if (new) then
   write(6,*)"open file='topo2' lui=",lui
   ierr = nf90_open(trim(filepath10km)//'/'//'topo2.nc',NF90_NOWRITE,ncid)
   if ( ierr==nf90_noerr ) then
     ncfile=.true.
     write(6,*) "Found netcdf version"
     ierr = nf90_inq_varid(ncid,"topo",varid)
   else
     ncfile=.false.
     open (lui, file=trim(filepath10km)//'/'//'topo2',status='old', access='direct', form='unformatted', recl=nlong*2, &
           iostat=ios)
   end if  
   new = .false.
end if

i = nint((90.0 - lat) * 12.0 + 1.0)
if ( i.lt.3 ) print *,'i=',i

! ht_read buffers terrain height with 1 being a higher latitude (i-1),
! 2 being the current latitude (i) and 3 being the lower latitude
! (i+1).  Only ht_read(:,1) needs to be read for the next loop as it moves
! from the south pole (i=2161) to the north pole (i=1)

if ( ncfile ) then
  if ( i==2161 ) then ! first call (south pole)
    ht_read(1:nlong,3) = 0.  
    ierr = nf90_get_var(ncid,varid,ht_read(1:nlong,1:2),start=(/1,2160/),count=(/nlong,2/))
  else if ( i>1 ) then
    ht_read(1:nlong,3) = ht_read(1:nlong,2)      
    ht_read(1:nlong,2) = ht_read(1:nlong,1)
    ierr = nf90_get_var(ncid,varid,ht_read(1:nlong,1),start=(/1,i-1/),count=(/nlong,1/))
  else ! last call (north pole)
    ht_read(1:nlong,3) = ht_read(1:nlong,2)
    ht_read(1:nlong,2) = ht_read(1:nlong,1)
    ht_read(1:nlong,1) = 0.
  end if
else    
  if ( i==2161 ) then ! first call (south pole)
    ht_read(1:nlong,3) = 0.    
    read (lui, rec=2161, iostat=ios) buf  
    do j = 1, nlong
      ival = 256 * ichar(buf(1,j)) + ichar(buf(2,j))
      if (ival .ge. 32768) then ! sign extend
        ival = ival - 65536
      end if
      ht_read(j,2) = ival
    end do
    read (lui, rec=2160, iostat=ios) buf  
    do j = 1, nlong
      ival = 256 * ichar(buf(1,j)) + ichar(buf(2,j))
      if (ival .ge. 32768) then ! sign extend
        ival = ival - 65536
      end if
      ht_read(j,1) = ival
    end do
  else if ( i>1 ) then
    ht_read(1:nlong,3) = ht_read(1:nlong,2)
    ht_read(1:nlong,2) = ht_read(1:nlong,1)
    read (lui, rec=i-1, iostat=ios) buf  
    do j = 1, nlong
      ival = 256 * ichar(buf(1,j)) + ichar(buf(2,j))
      if (ival .ge. 32768) then ! sign extend
        ival = ival - 65536
      end if
      ht_read(j,1) = ival
    end do
  else ! last call (north pole)
    ht_read(1:nlong,3) = ht_read(1:nlong,2)
    ht_read(1:nlong,2) = ht_read(1:nlong,1)
    ht_read(1:nlong,1) = 0.
  end if    
end if  

! return terrain height
ht(1:nlong) = ht_read(1:nlong,2)

! calculate slope
! as terread uses tiles, we estimate slope from the surrounding gradients
! which can be ommitted when near a boundary
if ( i==1 ) then ! north pole
  dx = 0.
  dy = 112000.*180./2160.  
  do j = 1,nlong
    !sx will always be zero at the pole
    sy_n = abs(ht_read(j,3) - ht_read(j,2))
    slope(j) = sy_n/dy
  end do  
else if ( i==2161 ) then ! south pole
  dx = 0.  
  dy = 112000.*180./2160.  
  do j = 1,nlong
    !sx will always be zero at the pole  
    sy_s = abs(ht_read(j,2) - ht_read(j,1))
    slope(j) = sy_s/dy
  end do    
else
  dx = cos(lat*3.1415927/180.)*112000.*360./real(nlong)
  dy = 112000.*180./2160.  
  sx_e = abs(ht_read(2,2) - ht_read(1,2))
  sx_w = abs(ht_read(1,2) - ht_read(nlong,2))
  sy_n = abs(ht_read(1,3) - ht_read(1,2))
  sy_s = abs(ht_read(1,2) - ht_read(1,1))
  slope(1) = 0.25*(sx_e/dx + sx_w/dx + sy_n/dy + sy_s/dy)
  do j = 2,nlong-1
    sx_e = abs(ht_read(j+1,2) - ht_read(j,2))
    sx_w = abs(ht_read(j,2) - ht_read(j-1,2))
    sy_n = abs(ht_read(j,3) - ht_read(j,2))
    sy_s = abs(ht_read(j,2) - ht_read(j,1))
    slope(j) = 0.25*(sx_e/dx + sx_w/dx + sy_n/dy + sy_s/dy)
  end do  
  sx_e = abs(ht_read(1,2) - ht_read(nlong,2))
  sx_w = abs(ht_read(nlong,2) - ht_read(nlong-1,2))
  sy_n = abs(ht_read(nlong,3) - ht_read(nlong,2))
  sy_s = abs(ht_read(nlong,2) - ht_read(nlong,1))
  slope(nlong) = 0.25*(sx_e/dx + sx_w/dx + sy_n/dy + sy_s/dy)
end if

end subroutine read_ht
