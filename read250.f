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
      
      subroutine read250(nx,ny,lons,lats,dl,debug,idia,jdia,il,
     &                   ncid,varid,ncfile)

      use ccinterp
      use netcdf_m
      use rwork

      integer, intent(in) :: il
      integer, intent(in) :: ncid, varid
      integer ierr
      integer jl
      integer(kind=2), dimension(3600*2800) :: zz
      integer, dimension(nx,ny) :: zz2d
      real dl
      integer nrecl
      real lons,lats

      logical debug,ok,ow
      logical patcha, patchb, patchc
      logical patchd, patche, patchf
      logical, intent(in) :: ncfile

      data ow/.false./

      jl=6*il

      write(6,*)"read250 nx,ny,lons,lats,dl=",nx,ny,lons,lats,dl,debug

! read data then close file
      if ( ncfile ) then
        ierr = nf90_get_var(ncid,varid,zz2d(1:nx,1:ny),
     &    start=(/1,1/),count=(/nx,ny/))
        ierr = nf90_close(ncid)
        zz(1:nx*ny) = reshape( zz2d(1:nx,1:ny), (/ nx*ny /) )
      else    
        read(13) (zz(i),i=1,nx*ny)
        close(13)
      end if  
      
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
        if(ow)write(6,*)jg,aglat

        ok = (aglat.gt.rlatn.and.aglat.lt.rlatx)

        if(ok)then
          !write(6,*)"jg,lat,latn,latx,lons=",jg,aglat,rlatn,rlatx,lons
          if(ow)write(6,*)"zz=",(zz(n),n=1+(jg-1)*nx,nx+(jg-1)*nx)

!#######################################################################
         do ig=1,nx
!#######################################################################

           n=ig+(jg-1)*nx
           
           ! MJT bug fix
           if (patcha.and.(jg>=526.and.jg<=527)) cycle
           if (patchb.and.(jg>=334.and.jg<=335)) cycle
           if (patchc.and.(jg>=706.and.jp<=707)) cycle
           if (patchd.and.(jg>=1538.and.jg<=1539)) cycle
           if (patche.and.(jg>=1.and.jg<=3)) cycle
           if (patchf.and.(jg>=91.and.jg<=92)) cycle
	   
           aglat=lats+real(jg-1)*dl
           aglon=lons+real(ig-1)*dl
           if(ow)write(6,*)ig,aglon

!-----------------------------------------------------------------------
           ok=(aglat.gt.rlatn.and.aglat.lt.rlatx)
           ok=ok .and. (aglon.gt.rlonn.and.aglon.lt.rlonx)
           if ( ok ) then
            if(ow) write(6,*)"ig,lon,lonn,lonx=",ig,aglon,rlonn,rlonx
!-----------------------------------------------------------------------

! compute model grid i,j
           call lltoijmod(aglon,aglat,alci,alcj,nface)  ! con-cubic/octagon
           lci = nint(alci)
           lcj = nint(alcj)
! convert to "double" (i,jg) notation
           lcj=lcj+nface*il
           if(ow)write(6,*)ig,jg,aglon,aglat,lci,lcj

! check to make sure within grid dimensions
           if(lci.gt.0.and.lci.le.il.and.lcj.gt.0.and.lcj.le.jl) then
	   
           if (idsrtm(lci,lcj)==0) then ! exclude points which already have srtm data
               
             if( zz(n).lt.1 ) then
! ocean point if zz < 1
               amask = 0.
               zs=0.
             else
! land point zz >= 1
               amask = 1.
               zs=real(zz(n))
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

             if ( debug ) then
                if ( lci.eq.idia .and. lcj.eq.jdia ) then
                 write(6,'("lci,j,ig,jg,zz,glo,gla,zss,lm, inum()="
     &   ,5i5,2f8.3,2f8.1,i6)') lci,lcj,ig,jg,zz(n),aglon,aglat
     &          ,zss(lci,lcj),almsk(lci,lcj), inum(lci,lcj)
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

      return
      end ! read250
