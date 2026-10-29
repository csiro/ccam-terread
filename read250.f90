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
      
subroutine read250(nx,ny,lons,lats,dl,debug,idia,jdia,il,ncid,varid,ncfile)

use ccinterp
use netcdf_m
use rwork
      
implicit none

integer, intent(in) :: il, idia, jdia
integer, intent(in) :: ncid, varid
integer, intent(in) :: nx, ny
integer ierr, i, ig, jg, n
integer jl, lci, lcj, nface
integer nrecl
integer(kind=2), dimension(:), allocatable :: zz
integer(kind=2), dimension(:), allocatable :: zz_fix
integer, dimension(:,:), allocatable :: zz2d
real, intent(in) :: dl
!real, dimension(:), allocatable :: slope_a
real aglat, aglon, alci, alcj, amask, zs
!real slope_l
!real sx_e, sx_w, sy_n, sy_s, dx, dy
real, intent(in) :: lons,lats
logical ok
logical patcha, patchb, patchc
logical patchd, patche, patchf
logical, intent(in) :: ncfile
logical, intent(in) :: debug

jl=6*il

write(6,*)"read250 nx,ny,lons,lats,dl=",nx,ny,lons,lats,dl,debug
      
allocate( zz(nx*ny), zz_fix(nx*ny) )
!allocate( slope_a(nx*ny) )

! read data then close file
if ( ncfile ) then
  allocate( zz2d(nx,ny) )  
  ierr = nf90_get_var(ncid,varid,zz2d(1:nx,1:ny),start=(/1,1/),count=(/nx,ny/))
  ierr = nf90_close(ncid)
  zz(1:nx*ny) = reshape( zz2d(1:nx,1:ny), (/ nx*ny /) )
  deallocate( zz2d )
else    
  read(13) (zz(i),i=1,nx*ny)
  close(13)
end if  

! fix missing values
zz_fix = zz
where ( zz<1 )
  zz_fix = 0
end where
      
!! calculate slope
!do jg = 1,ny
!  aglat=lats+real(jg-1)*dl  
!  dx = cos(aglat*3.1415927/180.)*112000.*dl
!  dy = 112000.*dl            
!  if ( jg==1 ) then
!    n = 1
!    sx_e = abs( real(zz_fix(n+1)) - real(zz_fix(n)) )
!    sy_n = abs( real(zz_fix(n+nx)) - real(zz_fix(n)) )
!    slope_a(n) = 0.5*(sx_e/dx+sy_n/dy)
!    do ig = 2,nx-1
!      n = ig
!      sx_e = abs( real(zz_fix(n+1)) - real(zz_fix(n)) )
!      sx_w = abs( real(zz_fix(n)) - real(zz_fix(n-1)) )
!      sy_n = abs( real(zz_fix(n+nx)) - real(zz_fix(n)) )
!      slope_a(n) = (1./3.)*(sx_e/dx+sx_w/dx+sy_n/dy)
!    end do    
!    n = nx
!    sx_w = abs( real(zz_fix(n)) - real(zz_fix(n-1)) )
!    sy_n = abs( real(zz_fix(n+nx)) - real(zz_fix(n)) )
!    slope_a(n) = 0.5*(sx_w/dx+sy_n/dy)
!  else if ( jg<ny ) then
!    n = 1 + (jg-1)*nx
!    sx_e = abs( real(zz_fix(n+1)) - real(zz_fix(n)) )
!    sy_n = abs( real(zz_fix(n+nx)) - real(zz_fix(n)) )
!    sy_s = abs( real(zz_fix(n)) - real(zz_fix(n-nx)) )
!    slope_a(n) = (1./3.)*(sx_e/dx+sy_n/dy+sy_s/dy)
!    do ig = 2,nx-1
!      n = ig + (jg-1)*nx  
!      sx_e = abs( real(zz_fix(n+1)) - real(zz_fix(n)) )
!      sx_w = abs( real(zz_fix(n)) - real(zz_fix(n-1)) )
!      sy_n = abs( real(zz_fix(n+nx)) - real(zz_fix(n)) )
!      sy_s = abs( real(zz_fix(n)) - real(zz_fix(n-nx)) )
!      slope_a(n) = 0.25*(sx_e/dx+sx_w/dx+sy_n/dy+sy_s/dy)
!    end do    
!    n = nx + (jg-1)*nx
!    sx_w = abs( real(zz_fix(n)) - real(zz_fix(n-1)) )
!    sy_n = abs( real(zz_fix(n+nx)) - real(zz_fix(n)) )
!    sy_s = abs( real(zz_fix(n)) - real(zz_fix(n-nx)) )
!    slope_a(n) = (1./3.)*(sx_w/dx+sy_n/dy+sy_s/dy)
!  else
!    n = 1 + (ny-1)*nx
!    sx_e = abs( real(zz_fix(n+1)) - real(zz_fix(n)) )
!    sy_s = abs( real(zz_fix(n)) - real(zz_fix(n-nx)) )
!    slope_a(n) = 0.5*(sx_e/dx+sy_s/dy)
!    do ig = 2,nx-1
!      n = ig + (ny-1)*nx  
!      sx_e = abs( real(zz_fix(n+1)) - real(zz_fix(n)) )
!      sx_w = abs( real(zz_fix(n)) - real(zz_fix(n-1)) )
!      sy_s = abs( real(zz_fix(n)) - real(zz_fix(n-nx)) )
!      slope_a(n) = (1./3.)*(sx_e/dx+sx_w/dx+sy_s/dy)
!    end do    
!    n = nx + (ny-1)*nx
!    sx_w = abs( real(zz_fix(n)) - real(zz_fix(n-1)) )
!    sy_s = abs( real(zz_fix(n)) - real(zz_fix(n-nx)) )
!    slope_a(n) = 0.5*(sx_w/dx+sy_s/dy)
!  end if
!end do

! MJT bug fix
patcha = ( lons==144.001 .and. lats==-27.9988 )
patchb = ( lons==126.001 .and. lats==-19.9988 )
patchc = ( lons==138.001 .and. lats==-19.9988 )      
patchd = ( lons==143.819 .and. lats==-43.7412 )
patche = ( lons==132.001 .and. lats==-15.9987 )            
patchf = ( lons==120.001 .and. lats==-27.9988 )

!#######################################################################
do jg=1,ny
!#######################################################################

  aglat=lats+real(jg-1)*dl

  ok = (aglat.gt.rlatn.and.aglat.lt.rlatx)

  if(ok)then

!#######################################################################
    do ig=1,nx
!#######################################################################

      n=ig+(jg-1)*nx
           
      ! MJT bug fix
      ! original
      !if (patcha.and.(jg>=526.and.jg<=527)) cycle
      !if (patchb.and.(jg>=334.and.jg<=335)) cycle
      !if (patchc.and.(jg>=706.and.jg<=707)) cycle
      !if (patchd.and.(jg>=1538.and.jg<=1539)) cycle
      !if (patche.and.(jg>=1.and.jg<=3)) cycle
      !if (patchf.and.(jg>=91.and.jg<=92)) cycle
      ! patch avoid slope problems
      if (patcha.and.(jg>=525.and.jg<=528)) cycle
      if (patchb.and.(jg>=333.and.jg<=336)) cycle
      if (patchc.and.(jg>=705.and.jg<=708)) cycle
      if (patchd.and.(jg>=1537.and.jg<=1540)) cycle
      if (patche.and.(jg>=1.and.jg<=4)) cycle
      if (patchf.and.(jg>=90.and.jg<=93)) cycle
	  
      aglat=lats+real(jg-1)*dl
      aglon=lons+real(ig-1)*dl

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
	   
          if (idsrtm(lci,lcj)==0) then ! exclude points which already have srtm data
               
            if( zz(n).lt.1 ) then
! ocean point if zz < 1
              amask = 0.
              zs=0.
              !slope_l = 0.               
            else
! land point zz >= 1
              amask = 1.
              zs=real(zz(n))
              !slope_l = slope_a(n)
            end if  ! zs<-1000

            id250m(lci,lcj) = id250m(lci,lcj) + 1  
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
!! accumulate slope pnts
!            slope(lci,lcj)=slope(lci,lcj)+slope_l
!! sum of squares for sd calc.
!            slope_std(lci,lcj)=slope_std(lci,lcj)+slope_l**2

            if ( debug ) then
              if ( lci.eq.idia .and. lcj.eq.jdia ) then
                write(6,'("lci,j,ig,jg,zz,glo,gla,zss,lm, inum()=",5i5,2f8.3,2f8.1,i6)') &
                           lci,lcj,ig,jg,zz(n),aglon,aglat,zss(lci,lcj),almsk(lci,lcj), inum(lci,lcj)
              endif  ! selected points only
            endif  ! debug
             
          end if ! idsrtm

        endif ! (lci.gt.0.and.lci.le.il.and.lcj.gt.0.and.lcj.le.jl)then

!-----------------------------------------------------------------------
      endif ! ( ok ) then
!-----------------------------------------------------------------------

!#######################################################################
    enddo ! i
  endif! (ok)then
!#######################################################################
  if ( mod(jg,100).eq.0 .and. ok ) then
    write(6,*)"ig,jg,lon,lat,lci,lcj=",ig,jg,aglon,aglat,lci,lcj
  endif
!#######################################################################
enddo ! j
!#######################################################################

deallocate( zz, zz_fix )
!deallocate( slope_a )
        
return
end subroutine read250
