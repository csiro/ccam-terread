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
      
subroutine read10km(debug,do1km,il,cut10km,filepath10km)

use ccinterp
use rwork

implicit none

integer, parameter :: nx=4320
integer, parameter :: ny=2160

integer, intent(in) :: il
integer jl
integer izsmaxj, izsminj, j, i, lci, lcj, nface
integer, dimension(nx) :: izs
real dl
real lons,lats
real slope_l, aglat, aglon, alci, alcj, amask, zs
real, dimension(nx) :: slope_a
real, intent(in) :: cut10km
real, dimension(nx,3) :: ht_buffer
character(len=*), intent(in) :: filepath10km
logical, intent(in) :: debug, do1km

jl=6*il

lats=-90.
lons=0.
dl = 5./60.

!#######################################################################
do j=1,ny
!#######################################################################
  aglat = lats + (j-1)*dl  ! 0.+...

  call read_ht ( aglat, izs, slope_a, ht_buffer, filepath10km )

  izsmaxj=-99999
  izsminj=+99999
!#######################################################################
  do i=1,nx
!#######################################################################

    aglon=(lons)+real(i-1)*dl

!-----------------------------------------------------------------------
    izsmaxj=max(izsmaxj,izs(i))
    izsminj=min(izsminj,izs(i))

    call lltoijmod(aglon,aglat,alci,alcj,nface)  ! con-cubic
    lci = nint(alci)
    lcj = nint(alcj)
! convert to "double" (i,j) notation
    lcj=lcj+nface*il
    
    if (grid(lci,lcj).ge.20..or..not.do1km) then

      if (lci.gt.0.and.lci.le.il.and.lcj.gt.0.and.lcj.le.jl)then

        if (idsrtm(lci,lcj)==0.and.id250m(lci,lcj)==0.and.id1km(lci,lcj)==0) then ! exclude points which already have srtm data
                
! all other points
! (fixup to 5 min data coastline)
          if ( izs(i).lt.cut10km ) then
! ocean point
            amask = 0.
            zs=0.
            slope_l = 0.
          else  ! zs>=-1000
! land point
            amask = 1.
            zs=real(izs(i))
            slope_l = slope_a(i)
          end if  ! zs<-1000

          id10km(lci,lcj) = id10km(lci,lcj) + 1  
! accumulate topog. pnts
          zss(lci,lcj)=zss(lci,lcj)+zs
! accumulate lmask pnts
          almsk(lci,lcj) = almsk(lci,lcj) + amask
! accumulate number of  pnts in grid box
          inum(lci,lcj) = inum(lci,lcj) + 1
! find max/min topog. pnts in grid box
          tmax(lci,lcj)=max(tmax(lci,lcj),zs)
          tmin(lci,lcj)=min(tmin(lci,lcj),zs)
! sum of squares for sd calc.
          tsd(lci,lcj)=tsd(lci,lcj)+zs**2
! accumulate slope pnts
          slope(lci,lcj) = slope(lci,lcj) + slope_l
! sume of ssquares for sd calc
          slope_std(lci,lcj) = slope_std(lci,lcj) + slope_l**2
             
        end if ! idsrtm

      endif ! (lci.gt.0.and.lci.le.il.and.lcj.gt.0.and.lcj.le.jl)then
    endif !if (grid(lci,lcj).ge.20..or..not.do1km) then

!-----------------------------------------------------------------------

!#######################################################################
  enddo ! i
!#######################################################################
  if (mod(j,100).eq.0) then
    write(6,'("i,j,aglon,aglat,lci,lcj=",2i8,2f10.2,4i8)') &
               i,j,aglon,aglat,lci,lcj,izsmaxj,izsminj
  endif
!#######################################################################
enddo ! j
!#######################################################################

return
end subroutine read10km
