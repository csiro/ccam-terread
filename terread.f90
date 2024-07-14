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
    
      program terread

      use ccinterp
      use netcdf_m
      use rwork
      
      implicit none
      
      include 'version.h'
      
      real rtd,dl,dlh,clons,clats,plon
      real ds,gridn,gridx
      real dlon,dlat
      real sbd,wbd,stl1,stl2
      real du,tanl,rnml
      integer nter,n,nproc
      integer idia,jdia,id,jd,numzer
      integer imin,imax,jmin,jmax,iq
      integer, dimension(2) :: iq2
      integer ier,ierr,idnc
      integer luout,h,v

      parameter (rtd=180.0/3.14159265 )
      parameter ( nter = 37 )

      logical olam

      logical debug, do1km, okx, oky, do250, dosrtm
      logical dosrtm1, do250lsm
      logical netout, topfilt
      logical ncfile
      character*1024 fileout
      character*9 formout

      integer nx,ny
      integer i,j
      integer il,jl,ilout
      integer ii,jj,kk,ll,nface,lci,lcj
      integer lcis,lcjs,incg,num250lsm
      integer ncid, varid
      integer il_old, nproc_old, work_old, helm_old, nxp_old, nyp_old
      integer il_test, nproc_test, work_test, helm_test, nxp_test, nyp_test
      integer il_new, nproc_new, work_new, helm_new
      integer il_max
      integer, dimension(2) :: iq_n, iq_s, iq_e, iq_w
      integer, dimension(2) :: ccdim,pxy
      integer, dimension(:), allocatable :: in,ie,is,iw
      integer*2, dimension(3601,3) :: idata1_read
      integer*2, dimension(1201,3) :: idata_read
      integer*2, dimension(3601) :: idata1
      integer*2, dimension(1201) :: idata
      integer, dimension(1201) :: i4data
      character*1 ns,ew
      character*1024 file
      character*1024 modislsm
      character*8 ausfile(nter)
      character*1 esg,nsg
      character*2 ncord
      character*3 ecord
      character*1024 fname
      character*1024 filepath10km
      character*1024 filepath1km
      character*1024 filepath250m
      character*1024 filepath250mlsm
      character*1024 filepathsrtm

      real lons,lats
      real lone,late
      real cutsrtm,cut10km,cut1km
      real rlong0,rlat0,schmidt,dst
      real zs,aglon,aglat,alci,alcj,atsd
      real slope_l, sx_e, sx_w, sy_n, sy_s, dx, dy
      real, dimension(2) :: lonlat
      real, dimension(:,:,:), allocatable :: rlld
      real, dimension(3601) :: slope1_a
      real, dimension(1201) :: slope_a
      
      data debug/.true./
      data do1km/.true./
      data do250/.true./
      data do250lsm/.false./
      data dosrtm/.false./
      data dosrtm1/.false./
      data netout/.false./
      data topfilt/.false./
      data idia/24/,jdia/72/
      data id/2/,jd/2/
      data il/48/
      data ds/0./
      data cutsrtm/0./,cut10km/0./,cut1km/0./
      data incg/0/,nproc/-1/
      data fileout/"top.nc"/
      data filepath10km/""/
      data filepath1km/""/
      data filepath250m/""/
      data filepath250mlsm/""/
      data filepathsrtm/""/

      namelist / topnml / ds, du, tanl, rnml, stl1, stl2, debug &
        ,luout, fileout, olam, wbd, sbd, dlon, dlat             &
        ,idia, jdia ,rlong0, rlat0, schmidt                     &
        ,do1km, do250, dosrtm, dosrtm1, id, jd, il, netout      &
        ,topfilt, cutsrtm, cut10km, cut1km, incg, filepath10km  &
        ,filepath1km, filepath250m, filepathsrtm                &
        ,do250lsm, filepath250mlsm, nproc

      ! Start banner
      write(6,*) "=============================================================================="
      write(6,*) "CCAM: Starting terread"
      write(6,*) "=============================================================================="
      write(6,*) version
      
      write(6,*) "Reading topnml namelist"
#ifdef unitnml
      open( unit=99, file='top.nml', status='old')
      read( 99, NML=topnml )
      close( 99 )
#else
      read ( 5,NML=topnml )
#endif
      write(6,*) "Finished reading namelist"
      write ( 6,topnml )

      
      ! optimise grid for requested processors
      if ( nproc > 0 ) then
        il_old = il       ! store original grid size
        call reducenproc(5,il_old,nproc,nproc_old,nxp_old,nyp_old)
        work_old = il_old**2*6/nproc_old ! estimated work with current grid size
        call helm_check(il_old,helm_old)
        nproc_new = nproc_old ! store expected nproc
        il_new = il_old       ! store expected grid size
        work_new = work_old   ! store expected work     
        helm_new = helm_old
        il_max = 2*il_new*ceiling(sqrt(real(nproc)/real(nproc_old))) ! estimate maximum grid size
        do il_test = il+1,il_max
          call reducenproc(5,il_test,nproc,nproc_test,nxp_test,nyp_test)
          if ( nproc_test <= nproc ) then
            work_new = il_new**2*6/nproc_new     ! work with the expected grid size
            work_test = il_test**2*6/nproc_test  ! work with the trial grid size
            if ( work_test < work_new ) then
              call helm_check(il_test,helm_test)
              if ( helm_test >= helm_new ) then
                il_new = il_test       ! update grid size if trial is successful
                nproc_new = nproc_test ! update nproc if trial is successful
                helm_new = helm_test
              end if  
            end if
          end if  
        end do
        write(6,*) "Optimise grid with:"
        write(6,*) "-> original il, new il       = ",il_old,il_new
        write(6,*) "-> original nproc, new nproc = ",nproc_old,nproc_new
        write(6,*) "-> original work, new work   = ",work_old,work_new
        write(6,*) "-> original helm, new helm   = ",helm_old,helm_new
        il = il_new
      end if    
      
      
      jl=6*il 
      call rworkalloc(il)

      if ( do250 ) then
        ausfile(1) = "sd52.ter"
        ausfile(2) = "sd53.ter"
        ausfile(3) = "sd54.ter"
        ausfile(4) = "se51.ter"
        ausfile(5) = "se52.ter"
        ausfile(6) = "se53.ter"
        ausfile(7) = "se54.ter"
        ausfile(8) = "se55.ter"
        ausfile(9) = "sf50.ter"
        ausfile(10) = "sf51.ter"
        ausfile(11) = "sf52.ter"
        ausfile(12) = "sf53.ter"
        ausfile(13) = "sf54.ter"
        ausfile(14) = "sf55.ter"
        ausfile(15) = "sg50.ter"
        ausfile(16) = "sg51.ter"
        ausfile(17) = "sg52.ter"
        ausfile(18) = "sg53.ter"
        ausfile(19) = "sg54.ter"
        ausfile(20) = "sg55.ter"
        ausfile(21) = "sg56.ter"
        ausfile(22) = "sh50.ter"
        ausfile(23) = "sh51.ter"
        ausfile(24) = "sh52.ter"
        ausfile(25) = "sh53.ter"
        ausfile(26) = "sh54.ter"
        ausfile(27) = "sh55.ter"
        ausfile(28) = "sh56.ter"
        ausfile(29) = "si50.ter"
        ausfile(30) = "si51.ter"
        ausfile(31) = "si53.ter"
        ausfile(32) = "si54.ter"
        ausfile(33) = "si55.ter"
        ausfile(34) = "si56.ter"
        ausfile(35) = "sj54.ter"
        ausfile(36) = "sj55.ter"
        ausfile(37) = "sk55.ter"
      end if

      ! set-up CC grid
      ccdim(1)=il
      ccdim(2)=jl
      lonlat(1)=rlong0
      lonlat(2)=rlat0
      allocate(rlld(ccdim(1),ccdim(2),2))
      allocate(in(il*jl),ie(il*jl),is(il*jl),iw(il*jl))
      call cgg2(rlld,grid,ccdim,lonlat,schmidt,dst,in,ie,is,iw)
      ds=dst
      rlond=rlld(:,:,1)
      rlatd=rlld(:,:,2)
      !call setxyz
      
      gridx=-1.
      rlatx=-999.
      rlonx=-999.
      gridn=+999.
      rlatn=+999.
      rlonn=+999.
      do j=1,jl
       do i=1,il
         if(rlond(i,j).gt.180.)rlond(i,j)=rlond(i,j)-360.
         gridx=max(gridx,grid(i,j))
         gridn=min(gridn,grid(i,j))
         if ( grid(i,j).lt.20. ) then    ! use 1km when grid spc < 20 km
           rlatx=max(rlatx,rlatd(i,j))
           rlonx=max(rlonx,rlond(i,j))
           rlatn=min(rlatn,rlatd(i,j))
           rlonn=min(rlonn,rlond(i,j))
         endif ! grid < 20 km
         if(rlond(i,j).gt.180.)rlond(i,j)=rlond(i,j)-360.	 
       enddo
      enddo
      write(6,*)"grid x,n=",gridx,gridn
      write(6,*)"rlon x,n=",rlonx,rlonn
      write(6,*)"rlat x,n=",rlatx,rlatn

!     nested model topography (output)
      if (.not.netout) then
        write(6,*) 'open',luout,fileout
        open(luout,file=fileout,form='formatted',status='unknown')
        call ncdf_setup('zsm.nc',idnc,3,il,jl,1,20020101,0000,rlong0,rlat0,schmidt)
      else
        call ncdf_setup(fileout,idnc,3,il,jl,1,20020101,0000,rlong0,rlat0,schmidt)
      end if

! netcdf file
      call nc2out(grid ,il,jl,1,1.,idnc,"grid","grid","km",0.,5000.)
      call nc2out(rlond,il,jl,1,1.,idnc,"lon","lon","degrees",0.,360.)
      call nc2out(rlatd,il,jl,1,1.,idnc,"lat","lat","degrees",-90.,90.)

! initialize min,max, and sd arrays
      do j=1,jl
       do i=1,il
        tmin(i,j)=99999.
        tmax(i,j)=-99999.
        tsd(i,j)=0.
        almsk(i,j)=0.
        zss(i,j)=0.
        inum(i,j)=0
        idsrtm1(i,j)=0
        idsrtm(i,j)=0
        id250m(i,j)=0
        id250lsm(i,j)=0
        id1km(i,j)=0
        id10km(i,j)=0
        slope(i,j)=0.
        slope_std(i,j)=0.
       enddo
      enddo
      
      if ( dosrtm1 ) then ! *.bil version
          
        write(6,*) "Process dosrtm1"
        
        do ii=-180,179
          write(6,*) "Searching ... ",ii
          if (ii.lt.0) then
            esg='w'
          else
            esg='e'
          end if
          write(ecord,'(I3.3)') abs(ii)
          do jj=-50,50
            if (jj.lt.0) then
              nsg='s'
            else
              nsg='n'
            end if
            write(ncord,'(I2.2)') abs(jj)
            fname=trim(filepathsrtm)//'/'//nsg//ncord//"_"//esg//ecord//'_1arc_v3.bil'
            
            open(13,file=fname,access='direct',form='unformatted', &
                 iostat=ierr,recl=7202,status='old')
            
            if (ierr.eq.0) then
            
              lons=real(ii)
              lone=lons+1.
              lats=real(jj)
              late=lats+1.
              
              okx = .false.
              okx = okx .or. (lons.ge.rlonn .and. lons.le.rlonx)
              okx = okx .or. (lone.ge.rlonn .and. lone.le.rlonx)
              okx = okx .or. (lons.le.rlonn .and. lone.ge.rlonx)
              oky = .false.
              oky = oky .or. (lats.ge.rlatn .and. lats.le.rlatx) ! start lat in area
              oky = oky .or. (late.ge.rlatn .and. late.le.rlatx) ! end lat in area
              oky = oky .or. (lats.le.rlatn .and. late.ge.rlatx) ! grid encompasses hr reg
              oky = oky .or. (lats.ge.rlatn .and. late.le.rlatx) ! grid within hr reg

              if ( okx .and. oky ) then
                write(6,*) "Reading ",trim(fname)
                write(6,*) "cutsrtm=",cutsrtm
                idata1_read = 0 ! buffer
                do kk = 1,3600 ! ignore extra overlap row
                  idata1_read(:,1) = idata1_read(:,2)  
                  idata1_read(:,2) = idata1_read(:,3)
                  idata1_read(:,3) = 0
                  if ( kk==1 ) then
                    read(13,rec=1) idata1_read(:,2)
                  end if  
                  if ( kk<3600 ) then
                    read(13,rec=kk+1) idata1_read(:,3)
                  end if  
                  idata1 = idata1_read(:,2)
                  aglat=late-real(kk-1)/3600.
                  ! calculate slope
                  dx = cos(aglat*3.1415927/180.)*112000./3600.
                  dy = 112000./3600.
                  if ( kk==1 ) then
                    sx_e = abs( real(idata1_read(2,2)) - real(idata1_read(1,2)) )  
                    sy_n = abs( real(idata1_read(1,3)) - real(idata1_read(1,2)) )
                    slope1_a(1) = 0.5*(sx_e/dx+sy_n/dy)
                    do ll = 2,3599
                      sx_e = abs( real(idata1_read(ll+1,2)) - real(idata1_read(ll,2)) )  
                      sx_w = abs( real(idata1_read(ll,2)) - real(idata1_read(ll-1,2)) )
                      sy_n = abs( real(idata1_read(ll,3)) - real(idata1_read(ll,2)) )
                      slope1_a(ll) = (1./3.)*(sx_e/dx+sx_w/dx+sy_n/dy)
                    end do
                    sx_w = abs( real(idata1_read(3600,2)) - real(idata1_read(3599,2)) )
                    sy_n = abs( real(idata1_read(3600,3)) - real(idata1_read(3600,2)) )
                    slope1_a(3600) = 0.5*(sx_w/dx+sy_n/dy)
                  else if ( kk<3600 ) then    
                    sx_e = abs( real(idata1_read(2,2)) - real(idata1_read(1,2)) )  
                    sy_n = abs( real(idata1_read(1,3)) - real(idata1_read(1,2)) )
                    sy_s = abs( real(idata1_read(1,2)) - real(idata1_read(1,1)) )
                    slope1_a(1) = (1./3.)*(sx_e/dx+sy_n/dy+sy_s/dy)
                    do ll = 2,3599
                      sx_e = abs( real(idata1_read(ll+1,2)) - real(idata1_read(ll,2)) )  
                      sx_w = abs( real(idata1_read(ll,2)) - real(idata1_read(ll-1,2)) )
                      sy_n = abs( real(idata1_read(ll,3)) - real(idata1_read(ll,2)) )
                      sy_s = abs( real(idata1_read(ll,2)) - real(idata1_read(ll,1)) )
                      slope1_a(ll) = 0.25*(sx_e/dx+sx_w/dx+sy_n/dy+sy_s/dy)
                    end do
                    sx_w = abs( real(idata1_read(3600,2)) - real(idata1_read(3599,2)) )
                    sy_n = abs( real(idata1_read(3600,3)) - real(idata1_read(3600,2)) )
                    sy_s = abs( real(idata1_read(3600,2)) - real(idata1_read(3600,1)) )
                    slope1_a(3600) = (1./3.)*(sx_w/dx+sy_n/dy+sy_s/dy)
                  else
                    sx_e = abs( real(idata1_read(2,2)) - real(idata1_read(1,2)) )  
                    sy_s = abs( real(idata1_read(1,2)) - real(idata1_read(1,1)) )
                    slope1_a(1) = 0.5*(sx_e/dx+sy_s/dy)
                    do ll = 2,3599
                      sx_e = abs( real(idata1_read(ll+1,2)) - real(idata1_read(ll,2)) )  
                      sx_w = abs( real(idata1_read(ll,2)) - real(idata1_read(ll-1,2)) )
                      sy_s = abs( real(idata1_read(ll,2)) - real(idata1_read(ll,1)) )
                      slope1_a(ll) = (1./3.)*(sx_e/dx+sx_w/dx+sy_s/dy)
                    end do
                    sx_w = abs( real(idata1_read(3600,2)) - real(idata1_read(3599,2)) )
                    sy_s = abs( real(idata1_read(3600,2)) - real(idata1_read(3600,1)) )
                    slope1_a(3600) = 0.5*(sx_w/dx+sy_s/dy)
                  end if
                  ! update cubic grid
                  do ll=1,3600
                    aglon=lons+real(ll-1)/3600.
                    call lltoijmod(aglon,aglat,alci,alcj,nface)
                    lcis = nint(alci)
                    lcjs = nint(alcj)
                    lcjs=lcj+nface*il
                    if (idata1(ll).gt.-32768) then
                      lcis = int(alci)
                      lcjs = int(alcj)
                      lcjs=lcjs+nface*il
                      
                      do lci=lcis,min(il,lcis+incg)
                      do lcj=lcjs,min(jl,lcjs+incg)
                    
                      zs=real(idata(ll))
                      slope_l=slope1_a(ll)
                      
                      if (zs.gt.cutsrtm) then
                        almsk(lci,lcj) = almsk(lci,lcj) + 1.
                      end if
                        
                      idsrtm1(lci,lcj) = idsrtm1(lci,lcj) + 1
                      inum(lci,lcj) = inum(lci,lcj) + 1 ! accumulate number of points in grid box
                      zss(lci,lcj) = zss(lci,lcj) + zs ! accumulate topog. pnts
                      tmax(lci,lcj)=max(tmax(lci,lcj),zs) ! find max/min topog. pnts in grid box
                      tmin(lci,lcj)=min(tmin(lci,lcj),zs)
                      tsd(lci,lcj)=tsd(lci,lcj)+zs**2 ! sum of squares for sd calc.
                      slope(lci,lcj)=slope(lci,lcj)+slope_l ! acculate slope pnts
                      slope_std(lci,lcj)=slope_std(lci,lcj)+slope_l**2 ! sum of squares for sd calc.
                      end do
                      end do
                      
                    end if ! idata1(ll).gt.-32768
                  end do ! ll
                end do   ! kk
              end if
              close(13)
            end if
          end do
        end do
          
      end if ! dosrtm1
      
      if ( dosrtm ) then ! *.hgt version (depreciated?)
        write(6,*) "Process dosrtm"
      
        do ii=-180,179
          write(6,*) "Searching ... ",ii
          if (ii.lt.0) then
            esg='W'
          else
            esg='E'
          end if
          write(ecord,'(I3.3)') abs(ii)
          do jj=-90,89
            if (jj.lt.0) then
              nsg='S'
            else
              nsg='N'
            end if
            write(ncord,'(I2.2)') abs(jj)
            fname=trim(filepathsrtm)//'/'//nsg//ncord//esg//ecord//'.hgt'
            ierr = nf90_open(trim(fname)//'.nc',NF90_NOWRITE,ncid)
            if ( ierr==nf90_noerr ) then
              ncfile = .true.
              write(6,*) "Found netcdf version"
              ierr = nf90_inq_varid(ncid,"topo",varid)
            else  
              ncfile = .false.  
              open (13,file=fname,access='direct',form='unformatted',convert='big_endian', &
                    iostat=ierr,recl=2402,status='old')
            end if  
	    
            if (ierr.eq.0) then
	    
              lons=real(ii)
              lone=lons+1.
              lats=real(jj)
              late=lats+1.

              okx = .false.
              okx = okx .or. (lons.ge.rlonn .and. lons.le.rlonx)
              okx = okx .or. (lone.ge.rlonn .and. lone.le.rlonx)
              okx = okx .or. (lons.le.rlonn .and. lone.ge.rlonx)
              oky = .false.
              oky = oky .or. (lats.ge.rlatn .and. lats.le.rlatx) ! start lat in area
              oky = oky .or. (late.ge.rlatn .and. late.le.rlatx) !  end  lat in area
              oky = oky .or. (lats.le.rlatn .and. late.ge.rlatx) ! grid encompasses hr reg
              oky = oky .or. (lats.ge.rlatn .and. late.le.rlatx) ! grid within hr reg
	    
              if (okx.and.oky) then
	            write(6,*) "Reading ",trim(fname)
                write(6,*) "cutsrtm=",cutsrtm
                idata1_read = 0
                do kk=1,1200 ! ignore extra overlap row
                  idata1_read(:,1) = idata1_read(:,2)  
                  idata1_read(:,2) = idata1_read(:,3)
                  idata1_read(:,3) = 0
                  if ( ncfile ) then
                    if ( kk==1 ) then  
                      ierr = nf90_get_var(ncid,varid,i4data,start=(/1,kk/),count=(/1201,1/))  
                      idata_read(:,2) = i4data
                    end if
                    if ( kk<1200 ) then
                      ierr = nf90_get_var(ncid,varid,i4data,start=(/1,kk+1/),count=(/1201,1/))  
                      idata_read(:,3) = i4data
                    end if    
                  else    
                    if ( kk==1 ) then  
                      read(13,rec=kk) idata_read(:,2)
                    end if
                    if ( kk<1200 ) then  
                      read(13,rec=kk+1) idata_read(:,3)
                    end if    
                  end if  
                  idata = idata_read(:,2)
                  aglat=late-real(kk-1)/1200.
                  ! calculate slope
                  dx = cos(aglat*3.1415927/180.)*112000./1200.
                  dy = 112000./1200.
                  if ( kk==1 ) then
                    sx_e = abs( real(idata_read(2,2)) - real(idata_read(1,2)) )  
                    sy_n = abs( real(idata_read(1,3)) - real(idata_read(1,2)) )
                    slope_a(1) = 0.5*(sx_e/dx+sy_n/dy)
                    do ll = 2,1199
                      sx_e = abs( real(idata_read(ll+1,2)) - real(idata_read(ll,2)) )  
                      sx_w = abs( real(idata_read(ll,2)) - real(idata_read(ll-1,2)) )
                      sy_n = abs( real(idata_read(ll,3)) - real(idata_read(ll,2)) )
                      slope_a(ll) = (1./3.)*(sx_e/dx+sx_w/dx+sy_n/dy)
                    end do
                    sx_w = abs( real(idata_read(1200,2)) - real(idata_read(1199,2)) )
                    sy_n = abs( real(idata_read(1200,3)) - real(idata_read(1200,2)) )
                    slope_a(1200) = 0.5*(sx_w/dx+sy_n/dy)
                  else if ( kk<1200 ) then    
                    sx_e = abs( real(idata_read(2,2)) - real(idata_read(1,2)) )  
                    sy_n = abs( real(idata_read(1,3)) - real(idata_read(1,2)) )
                    sy_s = abs( real(idata_read(1,2)) - real(idata_read(1,1)) )
                    slope_a(1) = (1./3.)*(sx_e/dx+sy_n/dy+sy_s/dy)
                    do ll = 2,1199
                      sx_e = abs( real(idata_read(ll+1,2)) - real(idata_read(ll,2)) )  
                      sx_w = abs( real(idata_read(ll,2)) - real(idata_read(ll-1,2)) )
                      sy_n = abs( real(idata_read(ll,3)) - real(idata_read(ll,2)) )
                      sy_s = abs( real(idata_read(ll,2)) - real(idata_read(ll,1)) )
                      slope_a(ll) = 0.25*(sx_e/dx+sx_w/dx+sy_n/dy+sy_s/dy)
                    end do
                    sx_w = abs( real(idata_read(1200,2)) - real(idata_read(1199,2)) )
                    sy_n = abs( real(idata_read(1200,3)) - real(idata_read(1200,2)) )
                    sy_s = abs( real(idata_read(1200,2)) - real(idata_read(1200,1)) )
                    slope_a(1200) = (1./3.)*(sx_w/dx+sy_n/dy+sy_s/dy)
                  else
                    sx_e = abs( real(idata_read(2,2)) - real(idata_read(1,2)) )  
                    sy_s = abs( real(idata_read(1,2)) - real(idata_read(1,1)) )
                    slope_a(1) = 0.5*(sx_e/dx+sy_s/dy)
                    do ll = 2,1199
                      sx_e = abs( real(idata_read(ll+1,2)) - real(idata_read(ll,2)) )  
                      sx_w = abs( real(idata_read(ll,2)) - real(idata_read(ll-1,2)) )
                      sy_s = abs( real(idata_read(ll,2)) - real(idata_read(ll,1)) )
                      slope_a(ll) = (1./3.)*(sx_e/dx+sx_w/dx+sy_s/dy)
                    end do
                    sx_w = abs( real(idata_read(1200,2)) - real(idata_read(1199,2)) )
                    sy_s = abs( real(idata_read(1200,2)) - real(idata_read(1200,1)) )
                    slope_a(1200) = 0.5*(sx_w/dx+sy_s/dy)
                  end if
                  ! update cubic grid
                  do ll=1,1200
                    aglon=lons+real(ll-1)/1200.
                    call lltoijmod(aglon,aglat,alci,alcj,nface)
                    if (idata(ll).gt.-32768) then
                      lcis = nint(alci)
                      lcjs = nint(alcj)
                      lcjs=lcjs+nface*il
                      
                      do lci=lcis,min(il,lcis+incg)
                      do lcj=lcjs,min(jl,lcjs+incg)
    
                      zs=real(idata(ll))
                      slope_l=slope_a(ll)
                      if( zs.gt.cutsrtm ) then
                        almsk(lci,lcj) = almsk(lci,lcj) + 1.
                      end if  ! zs<-1000
                      idsrtm(lci,lcj) = idsrtm(lci,lcj) + 1
                      inum(lci,lcj) = inum(lci,lcj) + 1 ! accumulate number of pnts in grid box
                      zss(lci,lcj) = zss(lci,lcj) + zs ! accumulate topog. pnts
                      tmax(lci,lcj)=max(tmax(lci,lcj),zs) ! find max/min topog. pnts in grid box
                      tmin(lci,lcj)=min(tmin(lci,lcj),zs)
                      tsd(lci,lcj)=tsd(lci,lcj)+zs**2 ! sum of squares for sd calc.
                      slope(lci,lcj)=slope(lci,lcj)+slope_l ! acculate slope pnts
                      slope_std(lci,lcj)=slope_std(lci,lcj)+slope_l**2 ! sum of squares for sd calc.
                      
                      end do
                      end do
                      
                    end if
                  end do
                end do
              end if
              if ( ncfile ) then
                ierr = nf90_close(ncid)
              else  
                close(13)
              end if  
            end if
          end do
        end do
      end if
      
!================================
      !if ( do250 ) then
      if ( do250 .and. ( rlatx .gt. rlatn ) ) then     ! MJT suggestion
!================================
      write(6,*)"Proccess Australian 250m files"

      do n = 1,37

        file=trim(filepath250m)//'/'//trim(ausfile(n))
        write(6,*) file
        ierr = nf90_open(trim(file)//'.nc',NF90_NOWRITE,ncid)
        if ( ierr==nf90_noerr ) then
          ncfile = .true.
          ierr = nf90_inq_varid(ncid,"topo",varid)
          ierr = nf90_get_att(ncid,nf90_global,"nx",nx)
          ierr = nf90_get_att(ncid,nf90_global,"ny",ny)
          ierr = nf90_get_att(ncid,nf90_global,"lons",lons)
          ierr = nf90_get_att(ncid,nf90_global,"lats",lats)
          ierr = nf90_get_att(ncid,nf90_global,"dl",dl)           
        else
          ncfile = .false.  
          open(13,file=file,form="unformatted")
          read(13) nx,ny,lons,lats,dl
        end if  
        lone=lons+(nx-1)*dl
        late=lats+(ny-1)*dl

        okx = .false.
        okx = okx .or. (lons.ge.rlonn .and. lons.le.rlonx)
        okx = okx .or. (lone.ge.rlonn .and. lone.le.rlonx)
        okx = okx .or. (lons.le.rlonn .and. lone.ge.rlonx)
        !okx = okx .or. (lons.ge.rlonn .and. lone.le.rlonx) ! MJT suggestion
        oky = .false.
        oky = oky .or. (lats.ge.rlatn .and. lats.le.rlatx) ! start lat in area
        oky = oky .or. (late.ge.rlatn .and. late.le.rlatx) !  end  lat in area
        oky = oky .or. (lats.le.rlatn .and. late.ge.rlatx) ! grid encompasses hr reg
        oky = oky .or. (lats.ge.rlatn .and. late.le.rlatx) ! grid within hr reg

        write(6,'("lons,lone=",2f6.1," lats,late=",2f6.1,2l4,2i6,a10)') lons,lone,lats,late,okx,oky,nx,ny,ausfile(n)

        if ( okx.and.oky ) then
          write(6,*)"call read250(nx,ny,lons,lats,dl,debug,idia,jdia)"
          call read250(nx,ny,lons,lats,dl,debug,idia,jdia,il,ncid,varid,ncfile)
        endif

      enddo ! n

      write(6,*)"Done with Australian 250m files"

!================================
      endif ! ( do250 ) then
!================================

!================================
      !if ( do250lsm ) then
      if ( do250lsm .and. ( rlatx .gt. rlatn ) ) then ! MJT suggrestion
!================================
      write(6,*)"Proccess modis 250m lsm files"

      do h = 0,37
      do v = 0,17

!modis_SRTM-LandWaterMask__LPDAAC__v5__250m-SIN-grid-tile-h31v08__2000-2009.nc
        write(modislsm,'("modis_SRTM-LandWaterMask__LPDAAC__v5__250m-SIN-grid-tile-h",i2.2,"v",i2.2,"__2000-2009.nc")') h,v
        file=trim(filepath250mlsm)//'/'//trim(modislsm)
        ierr = nf90_open(trim(file),NF90_NOWRITE,ncid)

!=================================================
        if ( ierr==nf90_noerr ) then

        write(6,*) "modislam=",trim(file)

        nx=4800
        ny=4800
        ncfile = .true.
        ierr = nf90_inq_varid(ncid,"lwd",varid)

        if ( ierr==nf90_noerr ) then

        lats = min((9.-v)*10.,(9.-(v+1))*10.)
        late = max((9.-v)*10.,(9.-(v+1))*10.)

        lons = min((h-18.)*10./cos(lats/rtd),(h-18.)*10./cos(late/rtd))
        lons = min(lons,((h+1)-18.)*10./cos(lats/rtd),((h+1)-18.)*10./cos(late/rtd))
        lone = max((h-18.)*10./cos(lats/rtd),(h-18.)*10./cos(late/rtd))
        lone = max(lone,((h+1)-18.)*10./cos(lats/rtd),((h+1)-18.)*10./cos(late/rtd))

        okx = .false.
        okx = okx .or. (lons.ge.rlonn .and. lons.le.rlonx)
        okx = okx .or. (lone.ge.rlonn .and. lone.le.rlonx)
        okx = okx .or. (lons.le.rlonn .and. lone.ge.rlonx)
        oky = .false.
        oky = oky .or. (lats.ge.rlatn .and. lats.le.rlatx) ! start lat in area
        oky = oky .or. (late.ge.rlatn .and. late.le.rlatx) !  end  lat in area
        oky = oky .or. (lats.le.rlatn .and. late.ge.rlatx) ! grid encompasses hr reg
        oky = oky .or. (lats.ge.rlatn .and. late.le.rlatx) ! grid within hr reg

        if ( okx.and.oky ) then
          write(6,*)"call read250lsm(nx,ny,lons,lats,dl,debug,idia,jdia,il,ncid,varid,ncfile,h,v)"
          call read250lsm(nx,ny,lons,lats,dl,debug,idia,jdia,il,ncid,varid,ncfile,h,v)
        endif

        end if  ! ierr var

        end if  ! ierr file

      enddo ! v
      enddo ! h
! end of loop over all modis panels

      write(6,*)"almsk"
      num250lsm=0 ! count number of grid points that have modis l/s mask data
      do j=1,jl
        do i=1,il
          if ( id250lsm(i,j) .gt. 0 ) then
             inumx(i,j)=nint(real(lsmask(i,j))/real(id250lsm(i,j)))
             num250lsm=num250lsm+1
          else
             inumx(i,j)=0
          endif
        enddo !i=1,il
      enddo !j=1,jl

      write(6,*)"Done with modis 250m files  num250lsm=",num250lsm

!================================
      endif ! ( do250lsm ) then
!================================

!================================
      write(6,*)"################## rlat x,n=",rlatx,rlatn
      do1km = do1km .and. ( rlatx .gt. rlatn )
      write(6,*)"###################### do1km=",do1km
      if ( do1km ) then
!================================
      write(6,*)"Processing DEM files"
      dl = 1.0/120.0
      dlh = 1.0/240.0
      nx=4800
      ny=6000
      do i = 1,9
       lons=-180.+real(i-1)*40.
       lone=lons+real(nx-1)*dl
       ew="W"
       if ( nint(lons).gt.0 ) ew="E"

       do j = 1,3
        lats=90.-real(j-1)*50.
        late=lats-real(ny-1)*dl 
        ns="N"
        if ( nint(lats).lt.0 ) ns="S"

        write(file,'(a1,i3.3,a1,i2.2,".DEM")') ew,abs(nint(lons)),ns,abs(nint(lats))

        okx = .false.
        okx = okx .or. (lons.gt.rlonn .and. lons.lt.rlonx)
        okx = okx .or. (lone.gt.rlonn .and. lone.lt.rlonx)
        okx = okx .or. (lons.lt.rlonn .and. lone.gt.rlonx)
        okx = okx .or. (lons.lt.-179. .and. rlonx.gt. 180.) ! fix for dateline data
        okx = okx .or. (lons.gt. 179. .and. rlonn.lt.-180.) ! fix for dateline data
        oky = .false.
        oky = oky .or. (lats.gt.rlatn .and. lats.lt.rlatx)
        oky = oky .or. (late.gt.rlatn .and. late.lt.rlatx)
        oky = oky .or. (lats.gt.rlatx .and. late.lt.rlatn)

        write(6,'("lons,lone=",2f6.1," lats,late=",2f6.1,a12,2l4)') lons,lone,lats,late,file,okx,oky
        file=trim(filepath1km)//'/'//trim(file)
        
        if ( okx.and.oky ) then
          write(6,*)"call read1km(file,nx,ny,lons,lats,debug,idia,jdia)"
! 9 sept 2004
          clons=lons+dlh
          clats=lats-dlh
          call read1km(file,nx,ny,clons,clats,debug,idia,jdia,il,cut1km)
! 9 sept 2004
        endif

       enddo ! j

      enddo ! i

      nx=7200
      ny=3600
      do i = 1,6
        lons=-180.+real(i-1)*60.
        lone=lons+real(nx-1)*dl
        ew="W"
        if ( nint(lons).gt.0 ) ew="E"
        lats=-60.
        late=lats-real(ny-1)*dl 
        ns="S"

        write(file,'(a1,i3.3,a1,i2.2,".DEM")') ew,abs(nint(lons)),ns,abs(nint(lats))

        okx = .false.
        okx = okx .or. (lons.gt.rlonn .and. lons.lt.rlonx)
        okx = okx .or. (lone.gt.rlonn .and. lone.lt.rlonx)
        okx = okx .or. (lons.lt.rlonn .and. lone.gt.rlonx)
        okx = okx .or. (lons.lt.-179. .and. rlonx.gt. 180.) ! fix for dateline data
        okx = okx .or. (lons.gt. 179. .and. rlonn.lt.-180.) ! fix for dateline data        
        oky = .false.
        oky = oky .or. (lats.gt.rlatn .and. lats.lt.rlatx)
        oky = oky .or. (late.gt.rlatn .and. late.lt.rlatx)
        oky = oky .or. (lats.gt.rlatx .and. late.lt.rlatn)

        write(6,'("lons,lone=",2f6.1," lats,late=",2f6.1,a12,2l4)') lons,lone,lats,late,file,okx,oky
        file=trim(filepath1km)//'/'//trim(file)

        if ( okx.and.oky ) then
          write(6,*)"call read1km(file,nx,ny,lons,lats,debug,idia,jdia)"
! 9 sept 2004
          clons=lons+dlh
          clats=lats-dlh
          call read1km(file,nx,ny,clons,clats,debug,idia,jdia,il,cut1km)
! 9 sept 2004
        endif

      enddo ! i

      write(6,*)"after do1km"

      write(6,*)"Done with DEM files"

!================================
      endif ! ( do1km ) then
!================================

!***********************************************************************

      write(6,*)"Now calling read10km"
      call read10km(debug,do1km,il,cut10km,filepath10km)
      write(6,*)"after read10km"

!================================

      ! PATCH !!!!
      Do i=1,il
        Do j=1,jl
          aglon=rlond(i,j)
          aglat=rlatd(i,j)
          If ((aglon.GT.135.).AND.(aglon.LT.140.)) then
            If ((aglat.GT.-32.).AND.(aglat.LT.-25.)) then
              almsk(i,j)=real(inum(i,j))
            End If
          End If
        End Do
      End Do

      write(6,*)"now compute grid average values"

      numzer=0
      imin=9999
      jmin=9999
      imax=0
      jmax=0

      do j=1,jl
       do i=1,il

         if ( rlond(i,j) .lt. 0. ) then
            plon=rlond(i,j)+360.
         else
            plon=rlond(i,j)
         endif
!mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
! to avoid using modis 250m lsm mask near dateline dues to missing data points in input dataset
! this will use the existing (srtm?) terrain data to define land/sea mask
         if ( plon .lt. 179. .or. plon .gt. 181. ) then
!mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
          if ( do250lsm .and. num250lsm .gt. 0 ) then
! set land-sea mask (rmsk) if using modis land-sea mask
           if ( id250lsm(i,j) .gt. 0 ) then
! here lsmask is accumlated ocean points
! so lsmask is zero over ocean
             almsk(i,j) = lsmask(i,j)/real(id250lsm(i,j))
! assume land point if less than half of pnts are water
             if ( almsk(i,j) .ge. 0.5 ) then
               rmsk(i,j) = 1.
             else 
               rmsk(i,j) = 0.
             endif 
           endif
           !write(6,*)i,j,inum(i,j),zss(i,j)
          endif !( do250lsm )
!mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
         else
          id250lsm(i,j)=0
         endif ! [[ rlond(i,j) .ge. 170. .and. rlond(i,j) .le. 180. ]]then
!mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm

          if ( inum(i,j).eq.0 ) then
! no input points found, assume sea level
! but fills in using neighboring points below
             numzer=numzer+1
             imin=min(i,imin)
             imax=max(i,imax)
             jmin=min(j,jmin)
             jmax=max(j,jmax)
             zss(i,j) = 0.
             rmsk(i,j)= 0.
             tsd(i,j)= 0.
             tmin(i,j)= 0.
             tmax(i,j)= 0.
             slope(i,j)=0.
             slope_std(i,j)=0.
          else
! now case where inum ne 0
             zss(i,j) = zss(i,j)/real(inum(i,j))
! compute sd
             tsd(i,j)=sqrt(abs(tsd(i,j)/real(inum(i,j))-zss(i,j)**2)) ! abs for rounding?
!            tsd(i,j)=sqrt(tsd(i,j)/real(inum(i,j))-zss(i,j)**2)

! compute slope
             slope(i,j) = slope(i,j)/real(inum(i,j))
! compute sd
             slope_std(i,j) = sqrt(abs(slope_std(i,j)/real(inum(i,j))-slope(i,j)**2))
             
!mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
! if no modis land-sea data, do following
             if ( id250lsm(i,j) .eq. 0 ) then

! compute fraction of gribbox with land
               almsk(i,j) = almsk(i,j)/real(inum(i,j))

! if less than half of gridbox is ocean, assume it is land pnt
               if ( almsk(i,j).le. 0.5 ) then
! ocean point
                 rmsk(i,j) = 0.
                 zss(i,j)=0. ! jjk added 18-8-2004
                 tsd(i,j)=0. !mjt 12-5-05
                 tmax(i,j)=0. !mjt 12-5-05
                 tmin(i,j)=0. !mjt 12-5-05
                 slope(i,j)=0.
                 slope_std(i,j)=0.
               else ! almsk
! land point
                 rmsk(i,j) = 1.
                 zss(i,j)=max(.1,zss(i,j))
               endif ! almsk

             endif !( id250lsm )
!mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm

          end if ! inum

       enddo ! i
      enddo ! j

! first fill 
      numzer = 999
      do while (numzer.gt.0)
!       this checks for non-data points, and inserts a neighbouring value
        print *,'**** numzer= ',numzer
        numzer=0
        do j=1,jl
         do i=1,il
          iq=i+(j-1)*il
          inumx(i,j)=inum(i,j)
          if(inum(i,j).eq.0)then
            iq2(:)=0
            iq_n(2) = (in(iq)-1)/il + 1
            iq_n(1) = in(iq) - (iq_n(2)-1)*il
            if(inum(iq_n(1),iq_n(2)).ne.0)then
              iq2(:)=iq_n(:)
            endif
            iq_e(2) = (ie(iq)-1)/il + 1
            iq_e(1) = ie(iq) - (iq_e(2)-1)*il
            if(inum(iq_e(1),iq_e(2)).ne.0)then
              iq2(:)=iq_e(:)
            endif
            iq_w(2) = (iw(iq)-1)/il + 1
            iq_w(1) = iw(iq) - (iq_w(2)-1)*il
            if(inum(iq_w(1),iq_w(2)).ne.0)then
              iq2(:)=iq_w(:)
            endif
            iq_s(2) = (is(iq)-1)/il + 1
            iq_s(1) = is(iq) - (iq_s(2)-1)*il
            if(inum(iq_s(1),iq_s(2)).ne.0)then
              iq2(:)=iq_s(:)
            endif
            if(iq2(1).ne.0)then
              zss(i,j)=zss(iq2(1),iq2(2))
              if ( .not. do250lsm ) rmsk(i,j)=rmsk(iq2(1),iq2(2))
              tsd(i,j)=tsd(iq2(1),iq2(2))
              tmin(i,j)=tmin(iq2(1),iq2(2))
              tmax(i,j)=tmax(iq2(1),iq2(2))
              almsk(i,j)=almsk(iq2(1),iq2(2))
              inumx(i,j)=1
              slope(i,j)=slope(iq2(1),iq2(2))
              slope_std(i,j)=slope_std(iq2(1),iq2(2))
            else
              numzer=numzer+1
            endif    !  (iq2.ne.0)
          endif      !  (inum(i,j).eq.0)
         enddo
        enddo
        do j=1,jl
         do i=1,il
          inum(i,j)=inumx(i,j)
         enddo ! i
        enddo ! j
      end do ! do while
! first fill 

! second fill to fix tsd
! now fill in tsd where only one data point exists over land
        numzer=0
        do j=1,jl
         do i=1,il
          if ( inum(i,j) .eq. 1 .and. zss(i,j) .gt. .1 ) then
            numzer=numzer+1
            inumx(i,j)=0
          else
            inumx(i,j)=1
          endif
         enddo ! i
        enddo ! j

      do while (numzer.gt.0)
!       this checks for non-data points, and inserts a neighbouring value
        print *,'**** numzer= ',numzer
        numzer=0
        do j=1,jl
         do i=1,il
          iq=i+(j-1)*il
          almsk(i,j)=inumx(i,j)
          atsd=0.
          ii=0
          if(inumx(i,j).eq.0)then
            iq2(:)=0
            iq_n(2) = (in(iq)-1)/il + 1
            iq_n(1) = in(iq) - (iq_n(2)-1)*il
            if(inumx(iq_n(1),iq_n(2)).ne.0)then
              atsd=atsd+tsd(iq_n(1),iq_n(2))
              iq2(:)=iq_n(:)
              ii=ii+1
            endif
            iq_e(2) = (ie(iq)-1)/il + 1
            iq_e(1) = ie(iq) - (iq_e(2)-1)*il
            if(inumx(iq_e(1),iq_e(2)).ne.0)then
              atsd=atsd+tsd(iq_e(1),iq_e(2))
              iq2(:)=iq_e(:)
              ii=ii+1
            endif
            iq_w(2) = (iw(iq)-1)/il + 1
            iq_w(1) = iw(iq) - (iq_w(2)-1)*il
            if(inumx(iq_w(1),iq_w(2)).ne.0)then
              atsd=atsd+tsd(iq_w(1),iq_w(2))
              iq2(:)=iq_w(:)
              ii=ii+1
            endif
            iq_s(2) = (is(iq)-1)/il + 1
            iq_s(1) = is(iq) - (iq_s(2)-1)*il
            if(inumx(iq_s(1),iq_s(2)).ne.0)then
              atsd=atsd+tsd(iq_s(1),iq_s(2))
              iq2(:)=iq_s(:)
              ii=ii+1
            endif
            if(ii.ne.0)then
              !tsd(i,j)=tsd(iq2(1),iq2(2))
              tsd(i,j)=atsd/real(ii)
              tmin(i,j)=tmin(iq2(1),iq2(2))
              tmax(i,j)=tmax(iq2(1),iq2(2))
              almsk(i,j)=1
            else
              numzer=numzer+1
            endif    !  (iq2.ne.0)
          endif      !  (inum(i,j).eq.0)
         enddo
        enddo
        do j=1,jl
         do i=1,il
          inumx(i,j)=nint(almsk(i,j))
         enddo ! i
        enddo ! j
      end do ! do while
!fix tsd
      
      if (topfilt) then
        write(6,*) "Filter orography"
        dum=zss
        do j=1,jl
          do i=1,il
            iq=i+(j-1)*il
            iq_n(2) = (in(iq)-1)/il + 1
            iq_n(1) = in(iq) - (iq_n(2)-1)*il
            iq_s(2) = (is(iq)-1)/il + 1
            iq_s(1) = is(iq) - (iq_s(2)-1)*il
            iq_e(2) = (ie(iq)-1)/il + 1
            iq_e(1) = ie(iq) - (iq_e(2)-1)*il
            iq_w(2) = (iw(iq)-1)/il + 1
            iq_w(1) = iw(iq) - (iq_w(2)-1)*il
            zss(i,j)=0.125*(dum(iq_n(1),iq_n(2))+dum(iq_s(1),iq_s(2))  &
                           +dum(iq_e(1),iq_e(2))+dum(iq_w(1),iq_w(2))) &
                           +0.5*dum(i,j)
          end do
        end do
      end if

      if (.not.netout) then
        write(luout,'(i3,i5,2f10.3,f6.3,f8.0,''  orog-mask-var'')') il,jl,rlong0,rlat0,schmidt,ds
      end if

      call nc2out(zss  ,il,jl,1,1.,idnc,"zs","zs","m",-100.,30000.)
      call nc2out(rmsk ,il,jl,1,1.,idnc,"lsm","lsm","none",-1.,1.)
      call nc2out(tsd  ,il,jl,1,1.,idnc,"tsd","tsd","m",0.,30000.)
      call nc2out(tmax ,il,jl,1,1.,idnc,"zmax","zmax","m",-100.,30000.)
      call nc2out(tmin ,il,jl,1,1.,idnc,"zmin","zmin","m",-100.,30000.)

      call nc2out(slope    ,il,jl,1,1.,idnc,"slope","slope","none",0.,10.)
      call nc2out(slope_std,il,jl,1,1.,idnc,"slope_std","slope_std","none",0.,10.)
      
!     write out g*zs(il,jl) to formatted file
      do j=1,jl
        do i=1,il
          zss(i,j)=9.80616*zss(i,j)
        end do ! i=1,il
      end do ! j=1,jl
      
      if (.not.netout) then
        ilout=min(il,30)
        write (formout,'(''('',i3,''f7.0)'')')ilout   !  i.e. (<il>f7.0)
        write(luout,formout) zss

!     write out land/sea mask to formatted file
        write (formout,'(''('',i3,''f4.1)'')')ilout   !  i.e. (<il>f4.1)
        write(luout,formout) rmsk

!     write out std.dev of top. to formatted file
        write (formout,'(''('',i3,''f6.0)'')')ilout   !  i.e. (<il>f6.0)
        write(luout,formout) tsd
      end if
      print *,'zss, rmsk and tsd written to unit',luout,' for il=',il

      do j=1,jl
        do i=1,il
          zss(i,j)=real(inum(i,j)) ! tmp array for printing
        end do ! i=1,il
      end do ! j=1,jl
      call nc2out(zss,il,jl,1,1.,idnc,"inum","inum","none",0.,2000.)

      do j=1,jl
        do i=1,il
          zss(i,j)=real(idsrtm1(i,j)) ! tmp array for printing
        end do ! i=1,il
      end do ! j=1,jl
      call nc2out(zss,il,jl,1,1.,idnc,"idsrtm1","idsrtm1","none",0.,2000.)

      do j=1,jl
        do i=1,il
          zss(i,j)=real(idsrtm(i,j)) ! tmp array for printing
        end do ! i=1,il
      end do ! j=1,jl
      call nc2out(zss,il,jl,1,1.,idnc,"idsrtm","idsrtm","none",0.,2000.)

      do j=1,jl
        do i=1,il
          zss(i,j)=real(id250m(i,j)) ! tmp array for printing
        end do ! i=1,il
      end do ! j=1,jl
      call nc2out(zss,il,jl,1,1.,idnc,"id250m","id250m","none",0.,2000.)

!mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
      if ( do250lsm ) then
      do j=1,jl
        do i=1,il
          zss(i,j)=(lsmask(i,j)) ! tmp array for printing
        end do ! i=1,il
      end do ! j=1,jl
      call nc2out(zss,il,jl,1,1.,idnc,"lsmask","lsmask","none",0.,100.)

      do j=1,jl
        do i=1,il
          zss(i,j)=real(id250lsm(i,j)) ! tmp array for printing
        end do ! i=1,il
      end do ! j=1,jl
      call nc2out(zss,il,jl,1,1.,idnc,"id250lsm","id250lsm","none",0.,2000.)
      endif
!mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm

      do j=1,jl
        do i=1,il
          zss(i,j)=real(id1km(i,j)) ! tmp array for printing
        end do ! i=1,il
      end do ! j=1,jl
      call nc2out(zss,il,jl,1,1.,idnc,"id1km","id1km","none",0.,2000.)

      do j=1,jl
        do i=1,il
          zss(i,j)=real(id10km(i,j)) ! tmp array for printing
        end do ! i=1,il
      end do ! j=1,jl
      call nc2out(zss,il,jl,1,1.,idnc,"id10km","id10km","none",0.,2000.)
      
      call ncclose(idnc)
      if (.not.netout) close(luout)
      
      deallocate(rlld)
      deallocate(in,ie,is,iw)
      call rworkdealloc

      ! Complete
      write(6,*) "CCAM: terread completed successfully"
      
      ! End banner
      write(6,*) "=============================================================================="
      write(6,*) "CCAM: Finished terread"
      write(6,*) "=============================================================================="
      
      stop
      end
    
      subroutine reducenproc(npanels,il_g,nproc,newnproc,nxp,nyp)
      implicit none
      integer, intent(in) :: il_g, nproc, npanels
      integer, intent(out) :: newnproc, nxp, nyp
      integer nproc_low, nxp_test, nyp_test
      nxp_test = 0
      nyp_test = 0
      ! try face decompositoin
      do nproc_low = nproc,1,-1
        call proctest_face(npanels,il_g,nproc_low,nxp_test,nyp_test)
        if ( nxp_test>0 ) exit
      end do
      newnproc = nproc_low
      nxp = nxp_test
      nyp = nyp_test
      return
      end subroutine reducenproc    
    
      subroutine proctest_face(npanels,il_g,nproc,nxp,nyp)
      implicit none
      integer, intent(in) :: il_g, nproc, npanels
      integer, intent(out) :: nxp, nyp
      integer jl_g
      if ( mod(nproc,6)/=0 .and. mod(6,nproc)/=0 ) then
        nxp = -1
        nyp = -1
      else
        jl_g = il_g + npanels*il_g                 ! size of grid along all panels (usually 6*il_g)
        nxp = max( 1, nint(sqrt(real(nproc)/6.)) ) ! number of processes in X direction
        nyp = nproc/nxp                            ! number of processes in Y direction
        ! search for valid process decomposition.  CCAM enforces the same grid size on each process
        do while ( (mod(il_g,max(nxp,1))/=0.or.mod(nproc/6,max(nxp,1))/=0.or.mod(jl_g,max(nyp,1))/=0) .and. nxp>0 )
          nxp = nxp - 1
          nyp = nproc/max(nxp,1)
        end do
      end if
      return
      end subroutine proctest_face    
    
      subroutine helm_check(il_g,g)
      implicit none
      integer, intent(in) :: il_g
      integer, intent(out) :: g
      integer gp
      g = 1
      gp = 2
      do while ( mod( il_g, gp )==0 )
        g = g + 1
        gp = 2*gp
      end do
      return
      end subroutine helm_check
!=======================================================================
