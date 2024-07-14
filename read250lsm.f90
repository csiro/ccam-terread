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
      
subroutine read250lsm(nx,ny,lons,lats,dl,debug,idia,jdia,il,ncid,varid,ncfile,mh,mv)

use ccinterp
use netcdf_m
use rwork
      
implicit none

integer, intent(in) :: il, idia, jdia
integer, intent(in) :: ncid, varid
integer, intent(in) :: nx, ny
integer, intent(in) :: mh, mv
integer ierr,ierrc
integer jl,zmax,zmin
integer ig, jg, lci, lcj, nface, n
integer, dimension(:,:), allocatable :: zz2d
real, intent(in) :: dl
real, intent(in) :: lons,lats
real fv,fh
real dtr
real alci, alcj, aglat, aglon
logical ok
logical, intent(in) :: ncfile
logical, intent(in) :: debug

zmin=10
zmax=-10

jl=6*il
dtr=3.14159265/180.

write(6,*)"read250lsm nx,ny,lons,lats=",nx,ny,lons,lats
write(6,*)"read250lsm debug,mh,mv=",debug,mh,mv
      
allocate( zz2d(nx,ny) )

! read data then close file
write(6,*) "ncfile=",ncfile
if ( ncfile ) then
  ierr=nf90_get_var(ncid,varid,zz2d,start=(/1,1/),count=(/nx,ny/))
  write(6,*)"varid,ierr=",varid,ierr
  ierrc = nf90_close(ncid)
  if ( ierr .ne. 0 ) return
  zmin=10
  zmax=-10
  do jg=1,ny
    do ig=1,nx
      zmax=max(zmax,zz2d(ig,jg))
      zmin=min(zmin,zz2d(ig,jg))
    enddo
  enddo
  write(6,*)"zz2d n,x=",zmin,zmax
  zmin=10
  zmax=-10
else
  return
end if  

!#######################################################################
write(6,*)"now loop through all points in this panel"
do jg=0,ny-1
!#######################################################################

  fv=real(mv)+(4800.-real(jg))/4800.

  aglat = (9.-fv)*10.

!#######################################################################
  do ig=0,nx-1
!#######################################################################

    n=(ig+1)+(jg)*nx
    fh=real(mh)+real(ig)/4800.
    aglon = (fh-18.)*10./cos(dtr*aglat)

!-----------------------------------------------------------------------
    ok=(aglat.gt.rlatn.and.aglat.lt.rlatx)
    ok=ok .and. (aglon.gt.rlonn.and.aglon.lt.rlonx)

    if ( ok ) then

!-----------------------------------------------------------------------

! compute model grid i,j
      call lltoijmod(aglon,aglat,alci,alcj,nface)  ! con-cubic/octagon
      lci = nint(alci)
      lcj = nint(alcj)
! convert to "double" (i,jg) notation
      lcj=lcj+nface*il

! check to make sure within grid dimensions
      if(lci.gt.0.and.lci.le.il.and.lcj.gt.0.and.lcj.le.jl) then

        if ( zz2d(ig+1,jg+1) .lt. 10 ) then

          if ( zz2d(ig+1,jg+1).eq.0) then
! accumulate lmask pnts
            lsmask(lci,lcj) = lsmask(lci,lcj) + 1
          end if

          id250lsm(lci,lcj) = id250lsm(lci,lcj) + 1  

          if ( debug ) then
            if ( lci.eq.idia .and. lcj.eq.jdia ) then
              write(6,'("lci,j,ig,jg,zz,glo,gla,lm, id250lsm=",4i5,4f9.3,i6)') &
                         lci,lcj,ig,jg,zz2d(ig+1,jg+1),aglon,aglat,lsmask(lci,lcj), id250lsm(lci,lcj)
            endif  ! selected points only
          endif  ! debug

          zmax=max(zmax,zz2d(ig+1,jg+1))
          zmin=min(zmin,zz2d(ig+1,jg+1))

        end if ! lt 10
             
      endif ! (lci.gt.0.and.lci.le.il.and.lcj.gt.0.and.lcj.le.jl)then

!-----------------------------------------------------------------------
    endif ! ( ok ) then
!-----------------------------------------------------------------------

!#######################################################################
  enddo ! i
!#######################################################################
  if ( mod(jg,100).eq.0 .and. ok ) then
    write(6,*)"ig,jg,lon,lat,lci,lcj=",ig,jg,aglon,aglat,lci,lcj
  endif
!#######################################################################
enddo ! j
write(6,*)"zmin,zmax=",zmin,zmax
!#######################################################################

deallocate( zz2d )
        
return
end subroutine read250lsm
