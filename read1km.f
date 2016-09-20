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
      
      subroutine read1km(file,nx,ny,lons,lats,debug,idia,jdia,il,cut1km)

      use ccinterp
      use rwork

      integer, intent(in) :: il
      integer jl
      character(len=*), intent(in) :: file
      integer*2 izs(7200)
      real dl
      integer nrecl
      real lons,lats
      real, intent(in) :: cut1km
      parameter ( zmin=-100 )

      logical debug,ok

      !include 'newmpar.h'
      !include 'rwork.h' ! rmsk,inum,inumx,zss,almsk,tmax,tmin,tsd,rlatd,rlond,grid,id1km,rlonx,rlonn,rlatx,rlatn

      jl=6*il

      write(6,*)"read1km file,nx,ny=",file,nx,ny

      dl = 1.0/120.0

! on sol-as
!      nrecl = nx*2
!      open(21,file=file,access='direct',form='unformatted',recl=nrecl)

! on jjk-as
      nrecl = nx*2
      open(CONVERT='BIG_ENDIAN',unit=21,file=file,access='direct'
     &            ,form='unformatted',recl=nrecl)

!#######################################################################
      do j=1,ny
!#######################################################################

        aglat=(lats)-real(j-1)*dl
        !write(6,*)j,aglat
        aglata=aglat

        ok = (aglat.gt.rlatn.and.aglat.lt.rlatx)

        !if ( ok ) write(6,*)"lat,latn,latx,lons=",aglat,rlatn,rlatx,lons

         read(21,rec=j) izs(1:nx)
         !rmax=-1.e35
         !rmin=+1.e35
         !do i=1,nx
         !  rmax=max(rmax,real(izs(i)))
         !  rmin=min(rmin,real(izs(i)))
         !enddo !i=1,nx
         !if(rmax.gt.-999.)write(6,*)"j,rmax,rmin=",j,rmax,rmin
        izsmaxj=-99999
        izsminj=+99999

       if ( ok ) then

!#######################################################################
         do i=1,nx
!#######################################################################

           n=i+(ny+1-j)*nx
           aglon=lons+real(i-1)*dl
!!! fix up for fiji latitude
   !        if (aglon > 170. . and . aglon < 180. .and. 
   !  .         aglata > -20. .and. aglata < -10. ) then
   !           aglat=aglata + .06
   !        else
              aglat=aglata
   !        endif

!-----------------------------------------------------------------------
          izsmaxj=max(izsmaxj,izs(i))
          izsminj=min(izsminj,izs(i))

           ok=(aglat.gt.rlatn.and.aglat.lt.rlatx)
           ok=ok .and. (aglon.gt.rlonn.and.aglon.lt.rlonx)

           if ( ok ) then

           !write(6,*)"i,lon,lonn,lonx=",i,aglon,rlonn,rlonx
!-----------------------------------------------------------------------

! compute model grid i,j
           call lltoijmod(aglon,aglat,alci,alcj,nface)  ! con-cubic
           lci = nint(alci)
           lcj = nint(alcj)
! convert to "double" (i,j) notation
           lcj=lcj+nface*il
           !write(6,*)i,j,aglon,aglat,lci,lcj

!          if ( inum(lci,lcj) .eq. 0 ) then ! this would only take first input point

! check to make sure within grid dimensions
           if(lci.gt.0.and.lci.le.il.and.lcj.gt.0.and.lcj.le.jl)then
             !write(6,*)"i,j,izs(i)=",i,j,izs(i)

           ! exclude points which already have srtm or 250m data
           if (idsrtm(lci,lcj)==0.and.id250m(lci,lcj)==0) then 
               
             if( izs(i).lt.cut1km ) then
! ocean point zs < 1
               amask = 0.
               zs=0.
             else  ! zs>=-1000
! land point zs >= 1
               amask = 1.
               zs=real(izs(i))
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

             if ( debug ) then
                if ( lci.eq.idia .and. lcj.eq.jdia ) then
                  write(6,'("lci,lcj,i,izs(i),zs,zss(),almsk(), inum()="
     &           ,4i6,3f10.3,i6)') lci,lcj,i,izs(i),zs,zss(lci,lcj)
     &                       ,almsk(lci,lcj), inum(lci,lcj)
                endif  ! selected points only
             endif  ! debug

           end if !  idsrtm 
             
           endif ! (lci.gt.0.and.lci.le.il.and.lcj.gt.0.and.lcj.le.jl)then

!          endif ! ( inum(i,j) .eq. 0 ) then

!-----------------------------------------------------------------------
           endif ! ( ok ) then
!-----------------------------------------------------------------------

!#######################################################################

         enddo ! i

!#######################################################################

        if ( mod(j,100).eq.0 ) then
          write(6,'("i,j,aglon,aglat,lci,lcj=",2i8,2f10.2,4i8)')
     &               i,j,aglon,aglat,lci,lcj,izsmaxj,izsminj
        endif

!#######################################################################

        endif! ( ok ) then

!#######################################################################
      enddo ! j
!#######################################################################

      close(21)

      return
      end ! read1km
c=======================================================================
