Program aeroemiss

! This code creates CCAM aerosol emission data using the CMIP5 datasets

Implicit None

Character*80, dimension(:,:), allocatable :: options
Character*160, dimension(13) :: fname
Character*160 topofile,so2_anth,so2_ship,so2_biom
Character*160 oc_anth,oc_ship,oc_biom
Character*160 bc_anth,bc_ship,bc_biom
Character*160 volcano,dmsfile,dustfile
Integer nopts,month

Namelist/aero/ topofile,month,so2_anth,so2_ship,so2_biom,oc_anth, &
               oc_ship,oc_biom,bc_anth,bc_ship,bc_biom,volcano,   &
               dmsfile,dustfile
                 

Write(6,*) 'AEROEMISS - CMIP5 aerosols to CC grid (MAR-13)'

! Read switches
nopts=1
Allocate (options(nopts,2))
options(:,1) = (/ '-o' /)
options(:,2) = ''

Call readswitch(options,nopts)
Call defaults(options,nopts)

! Read namelist
Write(6,*) 'Input &aero namelist'
Read(5,NML=aero)
Write(6,*) 'Namelist accepted'

! Generate veg data
fname(1)=topofile
fname(2)=so2_anth
fname(3)=bc_anth
fname(4)=oc_anth
fname(5)=so2_ship
fname(6)=bc_ship
fname(7)=oc_ship
fname(8)=so2_biom
fname(9)=bc_biom
fname(10)=oc_biom
fname(11)=volcano
fname(12)=dmsfile
fname(13)=dustfile

Call createaero(options,nopts,fname,month)

Deallocate(options)

Stop
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine displays the help message
!

Subroutine help()

Implicit None

Write(6,*)
Write(6,*) "Usage:"
Write(6,*) "  aeroemiss -o aero < aero.nml"
Write(6,*)
Write(6,*) "Options:"
Write(6,*) "  -o aero      Aerosol emissions output filename"
Write(6,*) "  aero.nml     Namelist file (see below)"
Write(6,*)
Write(6,*) "Namelist:"
Write(6,*) "  The namelist aero.nml specifies the CMIP5 files to"
Write(6,*) "  use for emission data.  The following example"
Write(6,*) "  illustrates the namelist syntax:"
Write(6,*)
Write(6,*) "  &aero"
Write(6,*) '    month    = 1'
Write(6,*) '    topofile = "topout"'
Write(6,*) '    so2_anth = "IPCC_emissions_RCP45_SO2_anthropogenic.nc"'
Write(6,*) '    so2_ship = "IPCC_emissions_RCP45_SO2_ships.nc"'
Write(6,*) '    so2_biom = "IPCC_emissions_RCP45_SO2_biomassburning.nc"'
Write(6,*) '    bc_anth  = "IPCC_emissions_RCP45_BC_anthropogenic.nc"'
Write(6,*) '    bc_ship  = "IPCC_emissions_RCP45_BC_ships.nc"'
Write(6,*) '    bc_biom  = "IPCC_emissions_RCP45_BC_biomassburning.nc"'
Write(6,*) '    oc_anth  = "IPCC_emissions_RCP45_OC_anthropogenic.nc"'
Write(6,*) '    oc_ship  = "IPCC_emissions_RCP45_OC_ships.nc"'
Write(6,*) '    oc_biom  = "IPCC_emissions_RCP45_OC_biomassburning.nc"'
Write(6,*) '    volcano  = "contineous_volc.nc"'
Write(6,*) '    dmsfile  = "dmsemiss.nc"'
Write(6,*) '    dustfile = "ginoux.nc"'
Write(6,*) "    /"
Write(6,*)
Write(6,*) "  where:"
Write(6,*) '    month         = Month for emissions'
Write(6,*) '    topofile      = Topography (input) file'
Write(6,*) '    so2_anth      = Anthropogenic SO2 emissions file'
Write(6,*) '    so2_ship      = Ships SO2 emissions file'
Write(6,*) '    so2_biom      = Biomass burning SO2 emissions file'
Write(6,*) '    bc_anth       = Anthropogenic BC emissions file'
Write(6,*) '    bc_ship       = Ships BC emissions file'
Write(6,*) '    bc_biom       = Biomass burning BC emissions file'
Write(6,*) '    oc_anth       = Anthropogenic OC emissions file'
Write(6,*) '    oc_ship       = Ships OC emissions file'
Write(6,*) '    oc_biom       = Biomass burning OC emissions file'
Write(6,*) '    volcano       = Volcanic emissions file'
Write(6,*) '    dmsfile       = DMS and natural organics emission file'
Write(6,*) '    dustfile      = Dust emission file'
Write(6,*)
Stop

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine determins the default values for the switches
!

Subroutine defaults(options,nopts)

Implicit None

Integer nopts
Character(len=*), dimension(nopts,2), intent(inout) :: options
Integer out
Integer locate

out=locate('-o',options(:,1),nopts)
If (options(out,2).EQ.'') then
  Write(6,*) "ERROR: Must specify output filename"
  stop
End if

Return
End


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine processes the sib data
!

Subroutine createaero(options,nopts,fname,month)

Use ccinterp

Implicit None

Integer, intent(in) :: nopts,month
Character(len=*), dimension(nopts,2), intent(in) :: options
Character*80, dimension(3) :: outputdesc
Character*160, dimension(13) :: fname
Character*80 returnoption,outfile
Character*45 header
real, dimension(:,:,:), allocatable :: rlld,aerosol
Real, dimension(:,:), allocatable :: gridout,lsdata,topdata
Real, dimension(2) :: lonlat
Real, dimension(3,2) :: alonlat
Real, dimension(1) :: alvl,atime
Real schmidt,dsx,ds
Integer, dimension(2) :: sibdim
Integer, dimension(4) :: dimnum,dimid,dimcount
Integer, dimension(0:4) :: ncidarr
Integer, dimension(6) :: adate
Integer, dimension(19) :: varid
Integer sibsize,tunit,i,j,k,ierr

outfile=returnoption('-o',options,nopts)

! Read topography file
tunit=1
call readtopography(tunit,fname(1),sibdim,lonlat,schmidt,dsx,header)

Write(6,*) "Dimension : ",sibdim
Write(6,*) "lon0,lat0 : ",lonlat
Write(6,*) "Schmidt   : ",schmidt
Allocate(gridout(sibdim(1),sibdim(2)),rlld(sibdim(1),sibdim(2),2))
Allocate(topdata(sibdim(1),sibdim(2)))
Allocate(lsdata(sibdim(1),sibdim(2)),aerosol(sibdim(1),sibdim(2),22))

call gettopols(tunit,fname(1),lsdata,sibdim)

! Determine lat/lon to CC mapping
Call ccgetgrid(rlld,gridout,sibdim,lonlat,schmidt,ds)

! Read CMIP5 aerosol data
Call getdata(aerosol,gridout,lsdata,rlld,sibdim,fname,month)

! Prep nc output
dimnum(1:2)=sibdim(1:2) ! CC grid dimensions
dimnum(3)=1 ! Turn off level
dimnum(4)=1 ! single month
adate=0 ! Turn off date
adate(2)=1 ! time units=months
Call ncinitcc(ncidarr,outfile,dimnum(1:3),dimid,adate)
call ncatt(ncidarr,'lat0',lonlat(2))
call ncatt(ncidarr,'lon0',lonlat(1))
call ncatt(ncidarr,'schmidt0',schmidt)
outputdesc=(/ 'so2a1', 'SO2 Anth lvl1 emission', 'kg m-2 s-1' /)
Call ncaddvargen(ncidarr,outputdesc,5,2,varid(1),1.,0.)
outputdesc=(/ 'so2a2', 'SO2 Anth lvl2 emission', 'kg m-2 s-1' /)
Call ncaddvargen(ncidarr,outputdesc,5,2,varid(2),1.,0.)
outputdesc=(/ 'bca1', 'BC Anth lvl1 emission', 'kg m-2 s-1' /)
Call ncaddvargen(ncidarr,outputdesc,5,2,varid(3),1.,0.)
outputdesc=(/ 'bca2', 'BC Anth lvl2 emission', 'kg m-2 s-1' /)
Call ncaddvargen(ncidarr,outputdesc,5,2,varid(4),1.,0.)
outputdesc=(/ 'oca1', 'OC Anth lvl1 emission', 'kg m-2 s-1' /)
Call ncaddvargen(ncidarr,outputdesc,5,2,varid(5),1.,0.)
outputdesc=(/ 'oca2', 'OC Anth lvl2 emission', 'kg m-2 s-1' /)
Call ncaddvargen(ncidarr,outputdesc,5,2,varid(6),1.,0.)
outputdesc=(/ 'so2b1', 'SO2 Bio lvl1 emission', 'kg m-2 s-1' /)
Call ncaddvargen(ncidarr,outputdesc,5,2,varid(7),1.,0.)
outputdesc=(/ 'so2b2', 'SO2 Bio lvl2 emission', 'kg m-2 s-1' /)
Call ncaddvargen(ncidarr,outputdesc,5,2,varid(8),1.,0.)
outputdesc=(/ 'bcb1', 'BC Bio lvl1 emission', 'kg m-2 s-1' /)
Call ncaddvargen(ncidarr,outputdesc,5,2,varid(9),1.,0.)
outputdesc=(/ 'bcb2', 'BC Bio lvl2 emission', 'kg m-2 s-1' /)
Call ncaddvargen(ncidarr,outputdesc,5,2,varid(10),1.,0.)
outputdesc=(/ 'ocb1', 'OC Bio lvl1 emission', 'kg m-2 s-1' /)
Call ncaddvargen(ncidarr,outputdesc,5,2,varid(11),1.,0.)
outputdesc=(/ 'ocb2', 'OC Bio lvl2 emission', 'kg m-2 s-1' /)
Call ncaddvargen(ncidarr,outputdesc,5,2,varid(12),1.,0.)
outputdesc=(/ 'dmso', 'DMS ocean emission', 'kg m-2 s-1' /)
Call ncaddvargen(ncidarr,outputdesc,5,2,varid(13),1.,0.)
outputdesc=(/ 'dmst', 'DMS land emission', 'kg m-2 s-1' /)
Call ncaddvargen(ncidarr,outputdesc,5,2,varid(14),1.,0.)
outputdesc=(/ 'ocna', 'Natural organic emission', 'kg m-2 s-1' /)
Call ncaddvargen(ncidarr,outputdesc,5,2,varid(15),1.,0.)
outputdesc=(/ 'vso2', 'Volcanic emissions', 'kg m-2 s-1' /)
Call ncaddvargen(ncidarr,outputdesc,5,2,varid(16),1.,0.)
outputdesc=(/ 'sandem', 'Sand fraction that can erode', 'none' /)
Call ncaddvargen(ncidarr,outputdesc,5,2,varid(17),1.,0.)
outputdesc=(/ 'siltem', 'Silt fraction that can erode', 'none' /)
Call ncaddvargen(ncidarr,outputdesc,5,2,varid(18),1.,0.)
outputdesc=(/ 'clayem', 'Clay fraction that can erode', 'none' /)
Call ncaddvargen(ncidarr,outputdesc,5,2,varid(19),1.,0.)
Call ncenddef(ncidarr)
alonlat(:,1)=(/ 1., real(sibdim(1)), 1. /)
alonlat(:,2)=(/ 1., real(sibdim(2)), 1. /)
alvl=1.
atime(1)=real(month)
Call nclonlatgen(ncidarr,dimid,alonlat,alvl,atime,dimnum)

Write(6,*) 'Write aerosol data'
dimcount=(/ sibdim(1), sibdim(2), 1, 1 /)
do i=1,19
  Call ncwritedatgen(ncidarr,aerosol(:,:,i),dimcount,varid(i))
end do
Call ncclose(ncidarr)

Deallocate(gridout,rlld,topdata,lsdata,aerosol)

Return
End

