! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2026 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
! This subroutine is to extract (in memory) data from the CMIP aerosol dataset.
!

Subroutine getdata(dataout,grid,lsdata,rlld,sibdim,fname,month,year)

Use ccinterp
use netcdf_m

Implicit None

Integer, intent(in) :: month, year
Integer, dimension(2), intent(in) :: sibdim
integer, dimension(sibdim(1),sibdim(2)) :: countt
integer, dimension(4,2) :: arrsize
integer, dimension(4) :: ncsize
integer, dimension(4) :: nstart, ncount
integer, dimension(1) :: minpos
integer, dimension(:,:,:), allocatable :: lcmap
integer ncstatus,ncid
integer i,j,n,ix,ii,jj,fp,pos,ind,lci,lcj,nface
integer cmipmode, valident
integer iarchi
Real, dimension(sibdim(1),sibdim(2),19), intent(out) :: dataout
real, dimension(sibdim(1),sibdim(2)) :: datatmp
Real, dimension(sibdim(1),sibdim(2)), intent(in) :: grid,lsdata
Real, dimension(sibdim(1),sibdim(2),2), intent(in) :: rlld
Real, dimension(:,:), allocatable :: coverout,tmpout
Real, dimension(2,2) :: emlonlat
real, dimension(:), allocatable :: rlat,dis
real aglon,aglat,alci,alcj,ssum
real alon_test_a, alon_test_b, alat_test_a, alat_test_b, xp, yp
real aa00, aa10, aa01, aa11
character(len=*), dimension(13), intent(in) :: fname
character*160, dimension(2) :: varname
character*3 :: aname
character(len=5) :: mip_text
logical ltest

dataout = 0.

Write(6,*) 'Process CMIP aerosol datasets'
do j = 1,3 ! 1=Anth,2=Shipping,3=Biomass burning
  do n = 1,3 ! 1=SO2,BC,OC
      
    ! read emission array
    fp = (j-1)*3 + n
    ncstatus = nf_open(fname(fp+1),nf_nowrite,ncid)
    If ( ncstatus/=nf_noerr ) Then
      Write(6,*) "ERROR: Error opening NetCDF file ",trim(fname(fp+1))," (",ncstatus,")"
      call finishbanner
      Stop -1
    End If 
    write(6,*) "Processing ",trim(fname(fp+1))
    call getncdims(ncid,ncsize)
    Call getnclonlat(ncid,emlonlat)
    
    ! determine MIP era for data formap
    cmipmode = 5
    ncstatus = nf_get_att_text(ncid,nf_global,"mip_era",mip_text)
    if ( ncstatus == nf_noerr ) then
      if ( mip_text == "CMIP7" ) then
        cmipmode = 7
      else if ( mip_text == "CMIP6" ) then
        cmipmode = 6
      end if
    else
      ! assume CMIP5
      cmipmode = 5
    end if
    write(6,*) "cmipmode = ",cmipmode
    
    arrsize = 1
    arrsize(1:2,2) = ncsize(1:2)
    if ( allocated( coverout ) ) then
      deallocate( coverout, tmpout ) 
      deallocate( lcmap )
    end if
    allocate( coverout(arrsize(1,2),arrsize(2,2)), tmpout(arrsize(1,2),arrsize(2,2)) )
    allocate( lcmap(arrsize(1,2),arrsize(2,2),2) )
   
    select case( cmipmode )
      case(5,6)
        if ( ncsize(4) == 12 ) then  
          iarchi = month  
        else
          write(6,*) "ERROR: Expecting 12 month file for cmipmode = ",cmipmode
          stop
        end if  
      case(7)
        call findarchi(ncid,iarchi,month,year)
      case default
        write(6,*) "ERROR: Unable to determine cmipmode"
        stop
    end select
    arrsize(4,1) = iarchi

!    ! check for sector
!    cmipmode = 5
!    ncstatus = nf_inq_varid(ncid,'sector',valident)
!    if ( ncstatus==nf_noerr ) then
!      cmipmode = 6
!    end if

    aname='ERR'
    select case(n)
      case(1)
        aname = "SO2"
      case(2)
        aname = "BC"
      case(3)
        aname = "OC"
    end select
    
    ! no upper level ship emissions
    if ( j==2 ) then
      ix = 1
    else
      ix = 2
    end if
    
    do i = 1,ix  ! 1=Level1,2=Upper level

      coverout = 0.

      pos = (j-1)*6 + (n-1)*2 + i
      select case(pos)
        case(1,3,5) ! SO2,BC,OC Anth Level1
          select case(cmipmode)
            case(6,7)  
              nstart(1) = 1
              nstart(2) = 1
              nstart(4) = iarchi
              ncount(1) = arrsize(1,2)
              ncount(2) = arrsize(2,2)
              ncount(3) = 1
              ncount(4) = 1
              nstart(3) = 1 ! sector=0 (Agriculture)
              ncstatus = nf90_inq_varid(ncid,trim(aname)//'_em_anthro',valident)
              ncstatus = nf90_get_var(ncid,valident,tmpout,start=nstart,count=ncount)
              where ( tmpout<1.e20)
                coverout = coverout + tmpout
              end where  
              nstart(3) = 5 ! sector=4 (Residential/Commercial)
              ncstatus = nf90_inq_varid(ncid,trim(aname)//'_em_anthro',valident)
              ncstatus = nf90_get_var(ncid,valident,tmpout,start=nstart,count=ncount)
              where ( tmpout<1.e20)
                coverout = coverout + tmpout
              end where  
              nstart(3) = 4 ! sector=3 (Transport)
              ncstatus = nf90_inq_varid(ncid,trim(aname)//'_em_anthro',valident)
              ncstatus = nf90_get_var(ncid,valident,tmpout,start=nstart,count=ncount)
              where ( tmpout<1.e20)
                coverout = coverout + tmpout
              end where  
              nstart(3) = 8 ! sector=7 (Waste)
              ncstatus = nf90_inq_varid(ncid,trim(aname)//'_em_anthro',valident)
              ncstatus = nf90_get_var(ncid,valident,tmpout,start=nstart,count=ncount)
              where ( tmpout<1.e20)
                coverout = coverout + tmpout
              end where  
            case(5)  
              varname(1)='emiss_awb'
              varname(2)='kg m-2 s-1'
              Call getmeta(ncid,varname,tmpout,arrsize)
              where ( tmpout<1.e20)
                coverout = coverout + tmpout
              end where  
              varname(1)='emiss_dom'
              varname(2)='kg m-2 s-1'
              Call getmeta(ncid,varname,tmpout,arrsize)
              where ( tmpout<1.e20)
                coverout = coverout + tmpout
              end where  
              varname(1)='emiss_tra'
              varname(2)='kg m-2 s-1'
              Call getmeta(ncid,varname,tmpout,arrsize)
              where ( tmpout<1.e20)
                coverout = coverout + tmpout
              end where  
              varname(1)='emiss_wst'
              varname(2)='kg m-2 s-1'
              Call getmeta(ncid,varname,tmpout,arrsize)
              where ( tmpout<1.e20)
                coverout = coverout + tmpout
              end where  
          end select
          ind = (n-1)*2 + 1 ! 1=so2a1,3=bca1,5=oca1

        case(2,4,6) ! SO2,BC,OC Anth Upper level
          select case(cmipmode)
            case(6,7)  
              nstart(1) = 1
              nstart(2) = 1
              nstart(4) = iarchi
              ncount(1) = arrsize(1,2)
              ncount(2) = arrsize(2,2)
              ncount(3) = 1
              ncount(4) = 1
              nstart(3) = 2 ! sector=1 (Energy)
              ncstatus = nf90_inq_varid(ncid,trim(aname)//'_em_anthro',valident)
              ncstatus = nf90_get_var(ncid,valident,tmpout,start=nstart,count=ncount)
              where ( tmpout<1.e20)
                coverout = coverout + tmpout
              end where  
              nstart(3) = 3 ! sector=2 (Industry)
              ncstatus = nf90_inq_varid(ncid,trim(aname)//'_em_anthro',valident)
              ncstatus = nf90_get_var(ncid,valident,tmpout,start=nstart,count=ncount)
              where ( tmpout<1.e20)
               coverout = coverout + tmpout
              end where  
            case(5)  
              varname(1)='emiss_ene'
              varname(2)='kg m-2 s-1'
              Call getmeta(ncid,varname,tmpout,arrsize)
              where ( tmpout<1.e20)
                coverout = coverout + tmpout
              end where  
              varname(1)='emiss_ind'
              varname(2)='kg m-2 s-1'
              Call getmeta(ncid,varname,tmpout,arrsize)
              where ( tmpout<1.e20)
                coverout = coverout + tmpout
              end where  
          end select 
          ind = n*2 ! 2=so2a2,4=bca2,6=oca2

        case(7,9,11) ! SO2,BC,OC Ship Level1
          select case(cmipmode)
            case(6,7)  
              nstart(1) = 1
              nstart(2) = 1
              nstart(4) = iarchi
              ncount(1) = arrsize(1,2)
              ncount(2) = arrsize(2,2)
              ncount(3) = 1
              ncount(4) = 1
              nstart(3) = 8 ! sector=7 (Ship)
              ncstatus = nf90_inq_varid(ncid,trim(aname)//'_em_anthro',valident)
              ncstatus = nf90_get_var(ncid,valident,tmpout,start=nstart,count=ncount)
              where ( tmpout<1.e20)
                coverout = coverout + tmpout
              end where  
            case(5)  
              varname(1)='emiss_shp'
              varname(2)='kg m-2 s-1'
              Call getmeta(ncid,varname,tmpout,arrsize)
              where ( tmpout<1.e20)
                coverout = coverout + tmpout
              end where  
          end select
          ind = (n-1)*2 + 1 ! 1=so2a1,3=bca1,5=oca1

        case(13,15,17) ! SO2,BC,OC Biomass burning level1
          select case(cmipmode)
            case(6,7)
              ncstatus = nf90_inq_varid(ncid,trim(aname)//'_em_openburning',valident)
              if ( ncstatus == nf90_noerr ) then
                nstart(1) = 1
                nstart(2) = 1
                nstart(4) = iarchi
                ncount(1) = arrsize(1,2)
                ncount(2) = arrsize(2,2)
                ncount(3) = 1
                ncount(4) = 1
                nstart(3) = 3 ! sector=2 (Grassland)
                ncstatus = nf90_get_var(ncid,valident,tmpout,start=nstart,count=ncount)
              else
                ncstatus = nf90_inq_varid(ncid,'grassfire',valident)  
                if ( ncstatus /= nf90_noerr ) then
                  write(6,*) "ERROR: Need to add grassfire to emissions"
                  stop
                end if
                nstart(1) = 1
                nstart(2) = 1
                nstart(3) = iarchi
                ncount(1) = arrsize(1,2)
                ncount(2) = arrsize(2,2)
                ncount(3) = 1
                ncstatus = nf90_get_var(ncid,valident,tmpout,start=nstart(1:3),count=ncount(1:3))
              end if
              where ( tmpout<1.e20)
                coverout = coverout + tmpout
              end where  
            case(5)
              varname(1)='grassfire'
              varname(2)='kg m-2 s-1'
              Call getmeta(ncid,varname,tmpout,arrsize)
              where ( tmpout<1.e20)
                coverout = coverout + tmpout
              end where  
          end select
          ind = (n-1)*2 + 7 ! 7=so2b1,9=bcb1,11=ocb1       

        case(14,16,18) ! SO2,BC,OC Biomass burning upper level
          select case(cmipmode)
            case(6,7)
              ncstatus = nf90_inq_varid(ncid,trim(aname)//'_em_openburning',valident)
              if ( ncstatus == nf90_noerr ) then
                nstart(1) = 1
                nstart(2) = 1
                nstart(4) = iarchi ! = month
                ncount(1) = arrsize(1,2)
                ncount(2) = arrsize(2,2)
                ncount(3) = 1
                ncount(4) = 1
                nstart(3) = 2 ! sector=1 (Forest)
                ncstatus = nf90_get_var(ncid,valident,tmpout,start=nstart,count=ncount)
              else
                ncstatus = nf90_inq_varid(ncid,'forestfire',valident)
                if ( ncstatus /= nf90_noerr ) then
                  write(6,*) "ERROR: Need to add forestfire to emissions"
                  stop
                end if
                nstart(1) = 1
                nstart(2) = 1
                nstart(3) = iarchi ! = month
                ncount(1) = arrsize(1,2)
                ncount(2) = arrsize(2,2)
                ncount(3) = 1
                ncstatus = nf90_get_var(ncid,valident,tmpout,start=nstart(1:3),count=ncount(1:3))
              end if
              where ( tmpout<1.e20)
                coverout = coverout + tmpout
              end where  
            case(5)
              varname(1)='forestfire'
              varname(2)='kg m-2 s-1'
              Call getmeta(ncid,varname,tmpout,arrsize)
              where ( tmpout<1.e20)
                coverout = coverout + tmpout
              end where  
          end select
          ind = (n-1)*2 + 8 ! 8=so2b2,10=bcb2,12=ocb2

        case DEFAULT
          write(6,*) "ERROR: Internal error determining emission dataset"
          call finishbanner
          stop -1
      end select
      
      datatmp=0.
      countt=0
      
      ! bin tracer
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(arrsize,emlonlat,sibdim,lcmap) &
!$OMP   PRIVATE(jj,aglat,ii,aglon,alci,alcj,nface,lci,lcj)
      do jj=1,arrsize(2,2)
        aglat=(emlonlat(2,2)-emlonlat(2,1))*real(jj-1)/real(arrsize(2,2)-1)+emlonlat(2,1)
        do ii=1,arrsize(1,2)          
          aglon=(emlonlat(1,2)-emlonlat(1,1))*real(ii-1)/real(arrsize(1,2)-1)+emlonlat(1,1)
          call lltoijmod(aglon,aglat,alci,alcj,nface)
          lci = nint(alci)
          lcj = nint(alcj)
          lcj = lcj+nface*sibdim(1)
          lcmap(ii,jj,1) = lci
          lcmap(ii,jj,2) = lcj
        end do
      end do
!$OMP END PARALLEL DO
      do jj=1,arrsize(2,2)
        do ii=1,arrsize(1,2)
          lci = lcmap(ii,jj,1)
          lcj = lcmap(ii,jj,2)
          ! bin emission
          if ((nint(lsdata(lci,lcj))==1.and.j/=2).or. &
              (nint(lsdata(lci,lcj))==0.and.j==2)) then
            datatmp(lci,lcj)=datatmp(lci,lcj)+coverout(ii,jj)
          end if
          countt(lci,lcj)=countt(lci,lcj)+1
        end do
      end do
  
      ! fill missing values
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(sibdim,countt,lsdata,j,rlld,emlonlat,arrsize,datatmp,coverout) &
!$OMP   PRIVATE(lci,lcj,aglon,aglat,ii,jj)
      do lcj=1,sibdim(2)
        do lci=1,sibdim(1)
          if (countt(lci,lcj)==0) then
            if ((nint(lsdata(lci,lcj))==1.and.j/=2).or. &
                (nint(lsdata(lci,lcj))==0.and.j==2)) then
              aglon=rlld(lci,lcj,1)
              if (aglon<emlonlat(1,1)) aglon=aglon+360.
              if (aglon>emlonlat(1,1)+360.) aglon=aglon-360.
              aglat=rlld(lci,lcj,2)
              ii=nint((aglon-emlonlat(1,1))*real(arrsize(1,2)-1)/(emlonlat(1,2)-emlonlat(1,1)))+1
              if (ii>arrsize(1,2)) ii=ii-arrsize(1,2)
              jj=nint((aglat-emlonlat(2,1))*real(arrsize(2,2)-1)/(emlonlat(2,2)-emlonlat(2,1)))+1
              jj=min(max(jj,1),arrsize(2,2))
              datatmp(lci,lcj)=coverout(ii,jj)
            end if      
            countt(lci,lcj)=1
          end if
        end do
      end do
!$OMP END PARALLEL DO

      ! add to other emissions
      dataout(:,:,ind)=dataout(:,:,ind)+datatmp/real(countt)

    end do

    ! close nc file
    ncstatus=nf_close(ncid)
    
  end do
end do

! convert so2 to S
dataout(:,:,1)=0.5*dataout(:,:,1)
dataout(:,:,2)=0.5*dataout(:,:,2)
dataout(:,:,7)=0.5*dataout(:,:,7)
dataout(:,:,8)=0.5*dataout(:,:,8)

deallocate(coverout,tmpout)
deallocate(lcmap)

!--------------------------------------------------------------------
! volcanic emissions
Write(6,*) "Process volcanic emissions"
ncstatus=nf_open(fname(11),nf_nowrite,ncid)
If (ncstatus/=nf_noerr) Then
  Write(6,*) "ERROR: Error opening NetCDF file ",trim(fname(11))," (",ncstatus,")"
  call finishbanner
  Stop -1
End If
write(6,*) "Processing ",trim(fname(11))
Call getncdims(ncid,ncsize)
Call getnclonlat(ncid,emlonlat)
arrsize=1
arrsize(1:2,2)=ncsize(1:2)
arrsize(4,1)=1

allocate(coverout(arrsize(1,2),arrsize(2,2)))
allocate(lcmap(arrsize(1,2),arrsize(2,2),2))
coverout=0.
countt=0

varname(1)='field'
varname(2)='kg/yr'
Call getmeta(ncid,varname,coverout,arrsize)

! bin tracer
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(arrsize,emlonlat,sibdim,lcmap) &
!$OMP   PRIVATE(jj,aglat,ii,aglon,alci,alcj,nface,lci,lcj)
do jj=1,arrsize(2,2)
  aglat=(emlonlat(2,2)-emlonlat(2,1))*real(jj-1)/real(arrsize(2,2)-1)+emlonlat(2,1)
  do ii=1,arrsize(1,2)          
    aglon=(emlonlat(1,2)-emlonlat(1,1))*real(ii-1)/real(arrsize(1,2)-1)+emlonlat(1,1)
    call lltoijmod(aglon,aglat,alci,alcj,nface)
    lci = nint(alci)
    lcj = nint(alcj)
    lcj = lcj+nface*sibdim(1)
    lcmap(ii,jj,1) = lci
    lcmap(ii,jj,2) = lcj
  end do
end do
!$OMP END PARALLEL DO
do jj=1,arrsize(2,2)
  do ii=1,arrsize(1,2)
    lci = lcmap(ii,jj,1)
    lcj = lcmap(ii,jj,2)
    ! bin emission
    if (nint(lsdata(lci,lcj))==1) then
      dataout(lci,lcj,16)=dataout(lci,lcj,16)+coverout(ii,jj)
    end if
    countt(lci,lcj)=countt(lci,lcj)+1
  end do
end do
  
! fill missing values
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(sibdim,countt,lsdata,rlld,emlonlat,arrsize,dataout,coverout) &
!$OMP   PRIVATE(lci,lcj,aglon,aglat,ii,jj)
do lcj=1,sibdim(2)
  do lci=1,sibdim(1)
    if (countt(lci,lcj)==0) then
      if (nint(lsdata(lci,lcj))==1) then
        aglon=rlld(lci,lcj,1)
        aglat=rlld(lci,lcj,2)
        if (aglon<emlonlat(1,1)) aglon=aglon+360.
        if (aglon>emlonlat(1,1)+360.) aglon=aglon-360.
        ii=nint((aglon-emlonlat(1,1))*real(arrsize(1,2)-1)/(emlonlat(1,2)-emlonlat(1,1)))+1
        if (ii>arrsize(1,2)) ii=ii-arrsize(1,2)
        jj=nint((aglat-emlonlat(2,1))*real(arrsize(2,2)-1)/(emlonlat(2,2)-emlonlat(2,1)))+1
        jj=min(max(jj,1),arrsize(2,2))
        dataout(lci,lcj,16)=coverout(ii,jj)
      end if      
      countt(lci,lcj)=1
    end if
  end do
end do
!$OMP END PARALLEL DO

ncstatus=nf_close(ncid)
      
! add to other emissions
dataout(:,:,16)=dataout(:,:,16)/real(countt)

! Normalise to 1Tg/yr (note grid is in km)
ssum=sum(dataout(:,:,16))
if ( ssum>0. ) then
  dataout(:,:,16)=dataout(:,:,16)*1.E9/(86400.*365.25*ssum*grid*grid*1.E6)
end if

deallocate(coverout)
deallocate(lcmap)

!--------------------------------------------------------------------
! process DMS and natural organic emissions
Write(6,*) "Process DMS and natural organic emissions"

ncstatus=nf_open(fname(12),nf_nowrite,ncid)
If (ncstatus.NE.nf_noerr) Then
  Write(6,*) "ERROR: Error opening NetCDF file ",trim(fname(12))," (",ncstatus,")"
  call finishbanner
  Stop -1
End If
write(6,*) "Processing ",trim(fname(12))
Call getncdims(ncid,ncsize)
Call getnclonlat(ncid,emlonlat)

allocate(coverout(ncsize(1),ncsize(2)))
allocate(rlat(ncsize(2)),dis(ncsize(2)))
allocate(lcmap(ncsize(1),ncsize(2),2))
coverout=0.

arrsize=1
arrsize(2,2)=ncsize(2)
varname(1)='latitude'
varname(2)='degrees_north'
Call getmeta(ncid,varname,rlat,arrsize)

arrsize=1
arrsize(1:2,2)=ncsize(1:2)
arrsize(4,1)=month

do n=1,3

  select case(n)
    case(1)
      varname(1)='dmssea'
      varname(2)='conc'
      Call getmeta(ncid,varname,coverout,arrsize)
    case(2)
      varname(1)='dmsterr'
      varname(2)='kg m-2 s-1'
      Call getmeta(ncid,varname,coverout,arrsize)
    case(3)
      varname(1)='natorg'
      varname(2)='kg m-2 s-1'
      Call getmeta(ncid,varname,coverout,arrsize)
  end select

  countt=0

  ! bin tracer
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(arrsize,emlonlat,sibdim,lcmap,rlat) &
!$OMP   PRIVATE(jj,aglat,ii,aglon,alci,alcj,nface,lci,lcj)
  do jj=1,arrsize(2,2)
    aglat=rlat(jj)
    do ii=1,arrsize(1,2)          
      aglon=(emlonlat(1,2)-emlonlat(1,1))*real(ii-1)/real(arrsize(1,2)-1)+emlonlat(1,1)
      call lltoijmod(aglon,aglat,alci,alcj,nface)
      lci = nint(alci)
      lcj = nint(alcj)
      lcj = lcj+nface*sibdim(1)
      lcmap(ii,jj,1) = lci
      lcmap(ii,jj,2) = lcj
    end do
  end do
!$OMP END PARALLEL DO
  do jj=1,arrsize(2,2)
    do ii=1,arrsize(1,2)
      lci = lcmap(ii,jj,1)
      lcj = lcmap(ii,jj,2)
      ! bin emission
      select case(n)
        case(1)
          ltest=nint(lsdata(lci,lcj))==0
        case(2,3)
          ltest=nint(lsdata(lci,lcj))==1
      end select
      if (ltest) then
        dataout(lci,lcj,12+n)=dataout(lci,lcj,12+n)+coverout(ii,jj)
      end if
      countt(lci,lcj)=countt(lci,lcj)+1
    end do
  end do
  
  ! fill missing values
  do lcj=1,sibdim(2)
    do lci=1,sibdim(1)
      if (countt(lci,lcj)<=4) then
        select case(n)
          case(1)
            ltest=nint(lsdata(lci,lcj))==0
          case(2,3)
            ltest=nint(lsdata(lci,lcj))==1
        end select
        if (ltest) then
          aglon=rlld(lci,lcj,1)
          aglat=rlld(lci,lcj,2)
          if (aglon<emlonlat(1,1)) aglon=aglon+360.
          if (aglon>emlonlat(1,1)+360.) aglon=aglon-360.
          ii=int((aglon-emlonlat(1,1))*real(arrsize(1,2)-1)/(emlonlat(1,2)-emlonlat(1,1)))+1
          if (ii>arrsize(1,2)) ii=ii-arrsize(1,2)
          if (ii==arrsize(1,2)) ii=ii-1
          alon_test_a=(emlonlat(1,2)-emlonlat(1,1))*real(ii-1)/real(arrsize(1,2)-1)+emlonlat(1,1)
          alon_test_b=(emlonlat(1,2)-emlonlat(1,1))*real(ii)/real(arrsize(1,2)-1)+emlonlat(1,1)
          xp = (aglon-alon_test_a)/(alon_test_b-alon_test_a)
          xp = max( min( xp, 1. ), 0. )
          ! non-uniform spacing for rlat?
          dis=(rlat-aglat)**2
          minpos=minloc(dis)
          jj=minpos(1)
          if ( rlat(jj)>aglat ) jj=jj-1
          if (jj==arrsize(2,2)) jj=jj-1
          if (jj==0) jj=jj+1
          alat_test_a=rlat(jj)
          alat_test_b=rlat(jj+1)
          yp = (aglat-alat_test_a)/(alat_test_b-alat_test_a)
          yp = max( min( yp, 1. ), 0. )
          aa00 = coverout(ii,jj)
          aa10 = coverout(ii+1,jj) - coverout(ii,jj)
          aa01 = coverout(ii,jj+1) - coverout(ii,jj)
          aa11 = coverout(ii+1,jj+1) - coverout(ii+1,jj) - coverout(ii,jj+1) + coverout(ii,jj)
          dataout(lci,lcj,12+n)=aa00 + aa10*xp + aa01*yp + aa11*xp*yp
        else
          dataout(lci,lcj,12+n) = 0.  
        end if      
        countt(lci,lcj)=1
      end if
    end do
  end do

  dataout(:,:,12+n)=dataout(:,:,12+n)/real(countt)

end do

ncstatus=nf_close(ncid)
      
deallocate(coverout,rlat,dis)
deallocate(lcmap)

!--------------------------------------------------------------------
! process Sand, Silt and Clay fraction that can erode for dust emissions
Write(6,*) "Process dust emission datasets"

ncstatus=nf_open(fname(13),nf_nowrite,ncid)
If (ncstatus.NE.nf_noerr) Then
  Write(6,*) "ERROR: Error opening NetCDF file ",trim(fname(13))," (",ncstatus,")"
  call finishbanner
  Stop -1
End If
write(6,*) "Processing ",trim(fname(13))
Call getncdims(ncid,ncsize)
Call getnclonlat(ncid,emlonlat)

allocate(coverout(ncsize(1),ncsize(2)))
allocate(rlat(ncsize(2)),dis(ncsize(2)))
allocate(lcmap(ncsize(1),ncsize(2),2))
coverout=0.

arrsize=1
arrsize(2,2)=ncsize(2)

varname(1)='latitude'
varname(2)='degrees_north'
Call getmeta(ncid,varname,rlat,arrsize)

arrsize=1
arrsize(1:2,2)=ncsize(1:2)
arrsize(4,1)=1

do n=1,3

  select case(n)
    case(1)
      varname(1)='sand'
      varname(2)='none'
      Call getmeta(ncid,varname,coverout,arrsize)
    case(2)
      varname(1)='silt'
      varname(2)='none'
      Call getmeta(ncid,varname,coverout,arrsize)
    case(3)
      varname(1)='clay'
      varname(2)='none'
      Call getmeta(ncid,varname,coverout,arrsize)
  end select

  countt=0

  ! bin tracer
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(arrsize,emlonlat,sibdim,lcmap,rlat) &
!$OMP   PRIVATE(jj,aglat,ii,aglon,alci,alcj,nface,lci,lcj)
  do jj=1,arrsize(2,2)
    aglat=rlat(jj)
    do ii=1,arrsize(1,2)          
      aglon=(emlonlat(1,2)-emlonlat(1,1))*real(ii-1)/real(arrsize(1,2)-1)+emlonlat(1,1)
      call lltoijmod(aglon,aglat,alci,alcj,nface)
      lci = nint(alci)
      lcj = nint(alcj)
      lcj = lcj+nface*sibdim(1)
      lcmap(ii,jj,1) = lci
      lcmap(ii,jj,2) = lcj
    end do
  end do
!$OMP END PARALLEL DO
  do jj=1,arrsize(2,2)
    do ii=1,arrsize(1,2)
      lci = lcmap(ii,jj,1)
      lcj = lcmap(ii,jj,2)
      ! bin emission
      if (nint(lsdata(lci,lcj))==1) then
        dataout(lci,lcj,16+n)=dataout(lci,lcj,16+n)+coverout(ii,jj)
      end if
      countt(lci,lcj)=countt(lci,lcj)+1
    end do
  end do
  
  ! fill missing values
  do lcj=1,sibdim(2)
    do lci=1,sibdim(1)
      if (countt(lci,lcj)<=4) then
        if (nint(lsdata(lci,lcj))==1) then
          aglon=rlld(lci,lcj,1)
          aglat=rlld(lci,lcj,2)
          if (aglon<emlonlat(1,1)) aglon=aglon+360.
          if (aglon>emlonlat(1,1)+360.) aglon=aglon-360.
          ii=int((aglon-emlonlat(1,1))*real(arrsize(1,2)-1)/(emlonlat(1,2)-emlonlat(1,1)))+1
          if (ii>arrsize(1,2)) ii=ii-arrsize(1,2)
          if (ii==arrsize(1,2)) ii=ii-1
          alon_test_a=(emlonlat(1,2)-emlonlat(1,1))*real(ii-1)/real(arrsize(1,2)-1)+emlonlat(1,1)
          alon_test_b=(emlonlat(1,2)-emlonlat(1,1))*real(ii)/real(arrsize(1,2)-1)+emlonlat(1,1)
          xp = (aglon-alon_test_a)/(alon_test_b-alon_test_a)
          xp = max( min( xp, 1. ), 0. )
          dis=(rlat-aglat)**2
          minpos=minloc(dis)
          jj=minpos(1)
          if ( rlat(jj)>aglat ) jj=jj-1
          if (jj==arrsize(2,2)) jj=jj-1
          if (jj==0) jj=jj+1
          alat_test_a=rlat(jj)
          alat_test_b=rlat(jj+1)
          yp = (aglat-alat_test_a)/(alat_test_b-alat_test_a)
          yp = max( min( yp, 1. ), 0. )
          aa00 = coverout(ii,jj)
          aa10 = coverout(ii+1,jj) - coverout(ii,jj)
          aa01 = coverout(ii,jj+1) - coverout(ii,jj)
          aa11 = coverout(ii+1,jj+1) - coverout(ii+1,jj) - coverout(ii,jj+1) + coverout(ii,jj)
          dataout(lci,lcj,16+n)=aa00 + aa10*xp + aa01*yp + aa11*xp*yp
        else
          dataout(lci,lcj,16+n)=0.  
        end if      
        countt(lci,lcj)=1
      end if
    end do
  end do

  dataout(:,:,16+n)=dataout(:,:,16+n)/real(countt)

end do

ncstatus=nf_close(ncid)
      
deallocate(coverout,rlat,dis)
deallocate(lcmap)

Write(6,*) "Task complete"

Return
End
    
subroutine findarchi(ncid,iarchi,month,year)

use netcdf_m

implicit none

integer, intent(in) :: ncid, month, year
integer, intent(out) :: iarchi
integer ncstatus, ldid, maxarchi, idvtime
integer kdate_rsav, kdate_r
integer year_r, year_s, month_s
integer allleap
integer(kind=8) mtimer
real timer
logical ltest
character(len=80) datestring
character(len=80) calendarstring

ncstatus = nf90_inq_dimid(ncid,'time',ldid)
if ( ncstatus /= nf90_noerr ) then
  write(6,*) nf90_strerror(ncstatus)
  stop
end if
ncstatus = nf90_inquire_dimension(ncid,ldid,len=maxarchi)
if ( ncstatus /= nf90_noerr ) then
  write(6,*) nf90_strerror(ncstatus)
  stop
end if
ncstatus = nf90_inq_varid(ncid,'time',idvtime)
if ( ncstatus /= nf90_noerr ) then
  write(6,*) nf90_strerror(ncstatus)
  stop
end if
ncstatus = nf90_get_att(ncid,idvtime,'units',datestring)
if ( ncstatus /= nf90_noerr ) then
  write(6,*) nf90_strerror(ncstatus)
  stop
end if
ncstatus = nf90_get_att(ncid,idvtime,'calendar',calendarstring)
if ( ncstatus /= nf90_noerr ) then
  write(6,*) nf90_strerror(ncstatus)
  stop
end if

call processdatestring(datestring,kdate_rsav)
call processcalendarstring(calendarstring,allleap)

! fast read
iarchi = 1
kdate_r = kdate_rsav
ncstatus = nf90_get_var(ncid,idvtime,timer,start=(/iarchi/))
if ( ncstatus /= nf90_noerr ) then
  write(6,*) nf90_strerror(ncstatus)
  stop
end if
mtimer = nint(timer,8)*1440_8 ! units=days
call datefix(kdate_r,mtimer,allleap)
year_r = kdate_r/10000
year_s = year
month_s = month
!iarchi = max(year_s - year_r - 1,0)*12 ! assume 1 value per month
iarchi = 0
! search
ltest = .true.
do while ( ltest .and. iarchi<maxarchi )
  iarchi = iarchi + 1  
  kdate_r = kdate_rsav
  ncstatus = nf90_get_var(ncid,idvtime,timer,start=(/iarchi/))
  if ( ncstatus /= nf90_noerr ) then
    write(6,*) nf90_strerror(ncstatus)
    stop
  end if
  mtimer = nint(timer,8)*1440_8 ! units=days
  call datefix(kdate_r,mtimer,allleap)
  ltest = (kdate_r/100-year_s*100-month_s)<0
end do
if ( ltest ) then
  write(6,*) "ERROR: Search failed with ltest,iarchi = ",ltest,iarchi
  write(6,*) "kdate_r = ",kdate_r
  stop
end if

return
end subroutine findarchi

subroutine processdatestring(datestring,kdate_rsav)

implicit none

integer, intent(out) :: kdate_rsav
integer iposa, iposb, ierx
integer yyyy, mm, dd
character(len=*), intent(in) :: datestring

! process year
iposa = index(trim(datestring),'since')
iposa = iposa + 5 ! skip 'since'
iposb = index(trim(datestring(iposa:)),'-')
iposb = iposa + iposb - 2 ! remove '-'
read(datestring(iposa:iposb),FMT=*,iostat=ierx) yyyy
if ( ierx/=0 ) then
  write(6,*) "ERROR reading time units.  Expecting year but found ",datestring(iposa:iposb)
  stop
end if

! process month
iposa = iposb + 2 ! skip '-'
iposb = index(trim(datestring(iposa:)),'-')
iposb = iposa + iposb - 2 ! remove '-'
read(datestring(iposa:iposb),FMT=*,iostat=ierx) mm
if ( ierx/=0 ) then
  write(6,*) "ERROR reading time units.  Expecting month but found ",datestring(iposa:iposb)
  stop
end if

! process day
iposa = iposb + 2 ! skip '-'
iposb = index(trim(datestring(iposa:)),' ')
iposb = iposa + iposb - 2 ! remove ' '
if ( iposb<iposa ) then
  read(datestring(iposa:),FMT=*,iostat=ierx) dd
  if ( ierx/=0 ) then
    write(6,*) "ERROR reading time units.  Expecting day but found ",datestring(iposa:)
    stop
  end if
else
  read(datestring(iposa:iposb),FMT=*,iostat=ierx) dd
  if ( ierx/=0 ) then
    write(6,*) "ERROR reading time units.  Expecting day but found ",datestring(iposa:iposb)
    stop
  end if
end if

! final date and time
kdate_rsav = yyyy*10000 + mm*100 + dd

return
end subroutine processdatestring    
    
subroutine processcalendarstring(calendarstring,allleap)

implicit none

integer, intent(out) :: allleap
character(len=*), intent(in) :: calendarstring

allleap = -1
select case(calendarstring)
  case("")
    allleap = 1 ! standard
  case("365_day")
    allleap = 0 ! 365day
  case default
    write(6,*) "ERROR: Unknown calendar = ",trim(calendarstring)
    stop
end select

return
end subroutine processcalendarstring

subroutine datefix(kdate_r,mtimer_r,allleap)

implicit none

integer, intent(inout) :: kdate_r
integer(kind=8), intent(inout) :: mtimer_r
integer, intent(in) :: allleap
integer(kind=8), dimension(12) :: mdays
integer, dimension(12) :: mdays4
!integer leap_l
integer(kind=8) iyr,imo,iday
integer(kind=8) mtimerh,mtimerm
integer(kind=8) mdays_save
integer(kind=8), parameter :: minsday = 1440

iyr=int(kdate_r,8)/10000_8
imo=(int(kdate_r,8)-10000_8*iyr)/100_8
iday=int(kdate_r,8)-10000_8*iyr-100_8*imo

call calendar_function(mdays4,kdate_r,allleap)
mdays(:) = int( mdays4(:) )
do while ( mtimer_r>minsday*mdays(imo) )
  mtimer_r=mtimer_r-minsday*mdays(imo)
  imo=imo+1_8
  if ( imo>12_8 ) then
    imo=1_8
    iyr=iyr+1_8
    if ( allleap==1 ) then
      mdays(2)=28_8      
      if ( mod(iyr,4_8)==0   ) mdays(2)=29_8
      if ( mod(iyr,100_8)==0 ) mdays(2)=28_8
      if ( mod(iyr,400_8)==0 ) mdays(2)=29_8
    end if
  end if
end do
  
iday=iday+mtimer_r/minsday
mtimer_r=mod(mtimer_r,minsday)
  
! at this point mtimer_r has been reduced to fraction of a day
  
!mdays_save=mdays(imo)
!imo=imo+(iday-1_8)/mdays(imo)
!iday=mod(iday-1_8,mdays_save)+1_8
!
!iyr=iyr+(imo-1_8)/12_8
!imo=mod(imo-1_8,12_8)+1_8

kdate_r=int(iday+100_8*(imo+100_8*iyr),4)
!mtimer_r = 0.
  
return
end subroutine datefix  
    
subroutine calendar_function(mdays,kdate,leap)

integer, dimension(1:12), intent(out) :: mdays
integer, intent(in) :: kdate, leap
integer iyr, month

iyr = kdate/10000
month = (kdate-10000*iyr)/100
if ( leap==cal_365 ) then ! 365 day calendar
  mdays=(/31,28,31,30,31,30,31,31,30,31,30,31/)
else if ( leap==cal_leap ) then ! 365/366 day calendar
  mdays=(/31,28,31,30,31,30,31,31,30,31,30,31/)
  if (mod(iyr,4)==0) mdays(2)=29
  if (mod(iyr,100)==0) mdays(2)=28
  if (mod(iyr,400)==0) mdays(2)=29
else if ( leap==cal_360 ) then ! 360 day calendar
  mdays=(/30,30,30,30,30,30,30,30,30,30,30,30/)
else
  write(6,*) "ERROR: Unknown option for leap = ",leap
  stop -1
end if

return
end subroutine calendar_function    
