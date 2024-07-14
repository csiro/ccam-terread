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
      
subroutine read1km(file,nx,ny,lons,lats,debug,idia,jdia,il,cut1km)

use ccinterp
use netcdf_m
use rwork
      
implicit none

integer, intent(in) :: il, nx, ny, idia, jdia
integer jl, i, j, n
character(len=*), intent(in) :: file
integer(kind=2), dimension(7200) :: izs
integer, dimension(7200) :: i4zs
integer, dimension(7200,3) :: i4zs_read
integer, dimension(7200,3) :: i4zs_fix
integer nrecl
integer ncid, varid, ierr
integer izsmaxj, izsminj, lci, lcj, nface
real, intent(in) :: lons, lats
real slope_l, aglat, aglata, aglon
real alci, alcj, zs, amask
real dx, dy, sx_e, sx_w, sy_n, sy_s
real dl
real, dimension(7200) :: slope_a
real, intent(in) :: cut1km
logical, intent(in) :: debug
logical ok
logical ncfile

jl = 6*il

write(6,*)"read1km file,nx,ny=",file,nx,ny

dl = 1.0/120.0
nrecl = nx*2
      
ierr = nf90_open(trim(file)//'.nc',NF90_NOWRITE,ncid)
if ( ierr==NF90_NOERR ) then
  ncfile = .true.
  write(6,*) "Found netcdf version"
  ierr = nf90_inq_varid(ncid,"topo",varid)
else
  ncfile = .false.
  open(CONVERT='BIG_ENDIAN',unit=21,file=file,access='direct',form='unformatted',recl=nrecl)
end if  

!#######################################################################
do j=1,ny
!#######################################################################

  aglat=(lats)-real(j-1)*dl
  aglata=aglat

  ok = (aglat.gt.rlatn.and.aglat.lt.rlatx)

  ! i4zs_read is a buffer with 1 being the previous row (j-1), 2 being the current
  ! row (j) and 3 being the next row (j+1).
  
  if ( ncfile ) then
    if ( j==1 ) then ! first call
      i4zs_read(1:nx,1) = 0
      ierr = nf90_get_var(ncid,varid,i4zs_read(1:nx,2:3),start=(/1,1/),count=(/nx,2/))     
    else if ( j<ny ) then
      i4zs_read(1:nx,1) = i4zs_read(1:nx,2)
      i4zs_read(1:nx,2) = i4zs_read(1:nx,3)
      ierr = nf90_get_var(ncid,varid,i4zs_read(1:nx,3),start=(/1,j+1/),count=(/nx,1/))    
    else ! last call
      i4zs_read(1:nx,1) = i4zs_read(1:nx,2)
      i4zs_read(1:nx,2) = i4zs_read(1:nx,3)
      i4zs_read(1:nx,3) = 0
    end if  
  else    
    if ( j==1 ) then ! first call
      i4zs_read(1:nx,1) = 0  
      read(21,rec=1) izs(1:nx)
      i4zs_read(1:nx,2) = izs(1:nx)
      read(21,rec=2) izs(1:nx)
      i4zs_read(1:nx,3) = izs(1:nx)
    else if ( j<ny ) then
      i4zs_read(1:nx,1) = i4zs_read(1:nx,2)
      i4zs_read(1:nx,2) = i4zs_read(1:nx,3)
      read(21,rec=j+1) izs(1:nx)
      i4zs_read(1:nx,3) = izs(1:nx)
    else ! last call
      i4zs_read(1:nx,1) = i4zs_read(1:nx,2)
      i4zs_read(1:nx,2) = i4zs_read(1:nx,3)
      i4zs_read(1:nx,3) = 0
    end if    
  end if 
  i4zs(1:nx) = i4zs_read(1:nx,2)
  izsmaxj=-99999
  izsminj=+99999
  
  ! i4zs_fix removes missing values over water
  i4zs_fix(1:nx,1:3) = i4zs_read(1:nx,1:3)
  where ( i4zs_fix(1:nx,1:3)==-9999 )
    i4zs_fix(1:nx,1:3) = 0
  end where  
        
  ! estimate slope
  ! approximate slope using surrounding gradients can can be droped when near a boundary
  if ( j==1 ) then ! bottom row
    dx = cos(aglat*3.1415927/180.)*112000.*dl
    dy = 112000.*dl
    sx_e = abs( real(i4zs_fix(2,2)) - real(i4zs_fix(1,2)) )
    sy_n = abs( real(i4zs_fix(1,3)) - real(i4zs_fix(1,2)) )
    slope_a(1) = 0.5*(sx_e/dx+sy_n/dy)
    do i = 2,nx-1
      sx_e = abs( real(i4zs_fix(i+1,2)) - real(i4zs_fix(i,2)) )
      sx_w = abs( real(i4zs_fix(i,2)) - real(i4zs_fix(i-1,2)) )
      sy_n = abs( real(i4zs_fix(i,3)) - real(i4zs_fix(i,2)) )
      slope_a(i) = (1./3.)*(sx_e/dx+sx_w/dx+sy_n/dy)
    end do  
    sx_w = abs( real(i4zs_fix(nx,2)) - real(i4zs_fix(nx-1,2)) )
    sy_n = abs( real(i4zs_fix(nx,3)) - real(i4zs_fix(nx,2)) )
    slope_a(nx) = 0.5*(sx_w/dx+sy_n/dy)
  else if ( j<ny ) then
    dx = cos(aglat*3.1415927/180.)*112000.*dl
    dy = 112000.*dl
    sx_e = abs( real(i4zs_fix(2,2)) - real(i4zs_fix(1,2)) )
    sy_n = abs( real(i4zs_fix(1,3)) - real(i4zs_fix(1,2)) )
    sy_s = abs( real(i4zs_fix(1,2)) - real(i4zs_fix(1,1)) )
    slope_a(1) = (1./3.)*(sx_e/dx+sy_n/dy+sy_s/dy)
    do i = 2,nx-1
      sx_e = abs( real(i4zs_fix(i+1,2)) - real(i4zs_fix(i,2)) )
      sx_w = abs( real(i4zs_fix(i,2)) - real(i4zs_fix(i-1,2)) )
      sy_n = abs( real(i4zs_fix(i,3)) - real(i4zs_fix(i,2)) )
      sy_s = abs( real(i4zs_fix(i,2)) - real(i4zs_fix(i,1)) )
      slope_a(i) = 0.25*(sx_e/dx+sx_w/dx+sy_n/dy+sy_s/dy)
    end do  
    sx_w = abs( real(i4zs_fix(nx,2)) - real(i4zs_fix(nx-1,2)) )
    sy_n = abs( real(i4zs_fix(nx,3)) - real(i4zs_fix(nx,2)) )
    sy_s = abs( real(i4zs_fix(nx,2)) - real(i4zs_fix(nx,1)) )
    slope_a(nx) = (1./3.)*(sx_w/dx+sy_n/dy+sy_s/dy)
  else ! top row
    dx = cos(aglat*3.1415927/180.)*112000.*dl
    dy = 112000.*dl
    sx_e = abs( real(i4zs_fix(2,2)) - real(i4zs_fix(1,2)) )
    sy_s = abs( real(i4zs_fix(1,2)) - real(i4zs_fix(1,1)) )
    slope_a(1) = 0.5*(sx_e/dx+sy_s/dy)
    do i = 2,nx-1
      sx_e = abs( real(i4zs_fix(i+1,2)) - real(i4zs_fix(i,2)) )
      sx_w = abs( real(i4zs_fix(i,2)) - real(i4zs_fix(i-1,2)) )
      sy_s = abs( real(i4zs_fix(i,2)) - real(i4zs_fix(i,1)) )
      slope_a(i) = (1./3.)*(sx_e/dx+sx_w/dx+sy_s/dy)
    end do  
    sx_w = abs( real(i4zs_fix(nx,2)) - real(i4zs_fix(nx-1,2)) )
    sy_s = abs( real(i4zs_fix(nx,2)) - real(i4zs_fix(nx,1)) )
    slope_a(nx) = 0.5*(sx_w/dx+sy_s/dy)
  end if   

  if ( ok ) then

!#######################################################################
    do i=1,nx
!#######################################################################

      n=i+(ny+1-j)*nx
      aglon=lons+real(i-1)*dl
      aglat=aglata

!-----------------------------------------------------------------------
      izsmaxj=max(izsmaxj,i4zs(i))
      izsminj=min(izsminj,i4zs(i))

      ok=(aglat.gt.rlatn.and.aglat.lt.rlatx)
      ok=ok .and. (aglon.gt.rlonn.and.aglon.lt.rlonx)

      if ( ok ) then

!-----------------------------------------------------------------------

! compute model grid i,j
        call lltoijmod(aglon,aglat,alci,alcj,nface)  ! con-cubic
        lci = nint(alci)
        lcj = nint(alcj)
! convert to "double" (i,j) notation
        lcj=lcj+nface*il

! check to make sure within grid dimensions
        if (lci.gt.0.and.lci.le.il.and.lcj.gt.0.and.lcj.le.jl) then

          ! exclude points which already have srtm or 250m data
          if (idsrtm(lci,lcj)==0.and.id250m(lci,lcj)==0) then 
               
            if( i4zs(i).lt.cut1km ) then
! ocean point zs < 1
              amask = 0.
              zs=0.
              slope_l = 0.
            else  ! zs>=-1000
! land point zs >= 1
              amask = 1.
              zs=real(i4zs(i))
              slope_l = slope_a(i)
            end if  ! zs<-1000

            id1km(lci,lcj) = id1km(lci,lcj) + 1
! accumulate number of pnts in grid box
            inum(lci,lcj) = inum(lci,lcj) + 1
! accumulate topog. pnts
            zss(lci,lcj) = zss(lci,lcj) + zs
! accumulate lmask pnts
            almsk(lci,lcj) = almsk(lci,lcj) + amask
! find max/min topog. pnts in grid box
            tmax(lci,lcj)=max(tmax(lci,lcj),zs)
            tmin(lci,lcj)=min(tmin(lci,lcj),zs)
! sum of squares for sd calc.
            tsd(lci,lcj)=tsd(lci,lcj)+zs**2
! accumulate slope pnts
            slope(lci,lcj) = slope(lci,lcj) + slope_l
! sum of squares for slope_std calc
            slope_std(lci,lcj) = slope_std(lci,lcj) + slope_l**2

            if ( debug ) then
              if ( lci.eq.idia .and. lcj.eq.jdia ) then
                write(6,'("lci,lcj,i,izs(i),zs,zss(),almsk(), inum()=",4i6,3f10.3,i6)') &
                           lci,lcj,i,i4zs(i),zs,zss(lci,lcj),almsk(lci,lcj),inum(lci,lcj)
              endif  ! selected points only
            endif  ! debug

          end if !  idsrtm 
             
        end if ! (lci.gt.0.and.lci.le.il.and.lcj.gt.0.and.lcj.le.jl)then

!-----------------------------------------------------------------------
      endif ! ( ok ) then
!-----------------------------------------------------------------------

!#######################################################################

    end do ! i

!#######################################################################

    if ( mod(j,100).eq.0 ) then
      write(6,'("i,j,aglon,aglat,lci,lcj=",2i8,2f10.2,4i8)') &
                 i,j,aglon,aglat,lci,lcj,izsmaxj,izsminj
    endif

!#######################################################################

  end if! ( ok ) then

!#######################################################################
end do ! j
!#######################################################################

if ( ncfile ) then
  ierr = nf90_close(ncid)
else  
  close(21)
end if  

return
end subroutine read1km
!=======================================================================
