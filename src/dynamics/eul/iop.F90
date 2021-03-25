module iop
!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: iop
! 
! !DESCRIPTION: 
! iop specific routines
!
! !USES:
!
  use cam_abortutils,   only: endrun
  use cam_logfile,      only: iulog
  use constituents,     only: readtrace, cnst_get_ind, pcnst, cnst_name
  use eul_control_mod,  only: eul_nsplit
  use netcdf,           only: NF90_NOERR,NF90_CLOSE,NF90_GET_VAR,NF90_INQUIRE_DIMENSION, &
                              NF90_INQ_DIMID, NF90_INQ_VARID, NF90_NOWRITE, NF90_OPEN, &
                              NF90_GET_ATT,NF90_GLOBAL,NF90_INQUIRE_ATTRIBUTE
  use phys_control,     only: phys_getopts
  use pmgrid,           only: beglat,endlat,plon,plev,plevp
  use prognostics,      only: n3,t3,q3,u3,v3,ps
  use scamMod,          only: use_camiop, ioptimeidx, have_ps, scm_backfill_iop_w_init, have_tsair, have_qref, &
                              tobs, have_t, tground, have_tg, qobs, have_q, have_cld,    &
                              have_clwp, divq, have_divq, vertdivq, have_vertdivq, divq3d, &
                              have_divq3d, dqfxcam, have_numliq, have_cldliq, have_cldice, &
                              have_numice, have_divu, have_divv, divt, have_divt, vertdivt, &
                              have_vertdivt, divt3d, have_divt3d, have_divu3d,  have_divv3d, &
                              have_ptend, ptend, wfld, uobs, have_u, uobs, vobs, have_v, &
                              vobs, have_prec, have_evap, have_q1, have_q2, have_lhflx, have_shflx, have_sflx, &
                              use_3dfrc, betacam, fixmascam, alphacam, doiopupdate, &
                              cldiceobs,  cldliqobs, cldobs, clwpobs, divu, &
                              divu3d, divv, divv3d, iopfile, lhflxobs, numiceobs, numliqobs, &
                              precobs, evapobs, q1obs, q2obs, scmlat, scmlon, shflxobs, tsair, qrefobs, have_omega, wfldh,qinitobs,sflxobs, &
                              have_wpqtpsfc_clasp, have_wpthlpsfc_clasp, have_qp2, have_rtp2_clasp, have_thlp2_clasp, have_rtpthlp_clasp, have_wp2_clasp, &
                              have_wp4_clasp, have_wp2thlp_clasp, have_wp2rtp_clasp, have_wprtp2_clasp, have_wpthlp2_clasp, &
                              have_wprtpthlp_clasp, have_wp3_clasp, have_upwp_clasp, have_uref_clasp, have_vref_clasp, have_vpwp_clasp, have_tsoil, have_lwupsrf, &
                              have_lwdnsrf, have_swupsrf, have_swdnsrf, &
                              lwdnsrfobs ,lwupsrfobs ,qp2obs ,swdnsrfobs ,swupsrfobs ,thlp2_clasp, rtp2_clasp, &
                              rtpthlp_clasp ,tsoilobs,upwp_clasp,vpwp_clasp,wp2_clasp ,wp2rtp_clasp, &
                              wp2thlp_clasp,wp3_clasp ,wp4_clasp ,wprtp2_clasp,wpqtpsfc_clasp,wpthlp2_clasp, &
                              wprtpthlp_clasp ,wpthlpsfc_clasp,uref_clasp,vref_clasp
  

  use shr_kind_mod,     only: r8 => shr_kind_r8, max_chars=>shr_kind_cl
  use shr_scam_mod,     only: shr_scam_GetCloseLatLon
  use spmd_utils,       only: masterproc
  use string_utils,     only: to_lower
  use time_manager,     only: timemgr_init, get_curr_date, get_curr_calday,&
                              get_nstep,is_first_step,get_start_date,timemgr_time_inc
  use wrap_nf,          only: wrap_inq_dimid,wrap_get_vara_realx
!
! !PUBLIC TYPES:
  implicit none


  private

  real(r8), allocatable,target :: dqfx3sav(:,:,:,:)       
  real(r8), allocatable,target :: t2sav(:,:,:)       
  real(r8), allocatable,target :: fusav(:,:,:)       
  real(r8), allocatable,target :: fvsav(:,:,:)       
  real(r8), allocatable,target :: divq3dsav(:,:,:,:)
  real(r8), allocatable,target :: divt3dsav(:,:,:)       
  real(r8), allocatable,target :: divu3dsav(:,:,:)       
  real(r8), allocatable,target :: divv3dsav(:,:,:)       
  real(r8), allocatable,target :: betasav(:)

  integer :: closelatidx,closelonidx,latid,lonid,levid,timeid

  real(r8):: closelat,closelon

!
! !PUBLIC MEMBER FUNCTIONS:
  public :: init_iop_fields
  public :: readiopdata         ! read iop boundary data
  public :: setiopupdate        ! find index in iopboundary data for current time
!  public :: scam_use_iop_srf
! !PUBLIC DATA:
  public betasav, &
         dqfx3sav, divq3dsav, divt3dsav,divu3dsav,divv3dsav,t2sav,fusav,fvsav

!
! !REVISION HISTORY:
! Created by John Truesdale
!
!EOP
!
! !PRIVATE MEMBER FUNCTIONS:
!----------------------------------------------------------------------- 

contains
   subroutine init_iop_fields()
!------------------------------------------------------------------------------
! Coupler for converting dynamics output variables into physics input variables
! also writes dynamics variables (on physics grid) to history file
!------------------------------------------------------------------------------
   implicit none
   character(len=*), parameter ::  sub = "init_iop_fields"
!-----------------------------------------------------------------------
   if (eul_nsplit>1) then
      call endrun('iop module cannot be used with eul_nsplit>1')
   endif
	        
   if(.not.allocated(betasav)) then
      allocate (betasav(beglat:endlat))
      betasav(:)=0._r8
   endif

   if(.not.allocated(dqfx3sav)) then
      allocate (dqfx3sav(plon,plev,pcnst,beglat:endlat))
      dqfx3sav(:,:,:,:)=0._r8
   endif
   if(.not.allocated(divq3dsav)) then
      allocate (divq3dsav(plon,plev,pcnst,beglat:endlat))
      divq3dsav(:,:,:,:)=0._r8
   endif
   if(.not.allocated(divt3dsav)) then
      allocate (divt3dsav(plon,plev,beglat:endlat))
      divt3dsav(:,:,:)=0._r8
   endif
   if(.not.allocated(divu3dsav)) then
      allocate (divu3dsav(plon,plev,beglat:endlat))
      divu3dsav(:,:,:)=0._r8
   endif
   if(.not.allocated(divv3dsav)) then
      allocate (divv3dsav(plon,plev,beglat:endlat))
      divv3dsav(:,:,:)=0._r8
   endif
   if(.not.allocated(t2sav)) then
      allocate (t2sav(plon,plev,beglat:endlat))  ! temp tendency
      t2sav(:,:,:)=0._r8
   endif
   if(.not.allocated(fusav)) then
      allocate (fusav(plon,plev,beglat:endlat))  ! U wind tendency
      fusav(:,:,:)=0._r8
   endif
   if(.not.allocated(fvsav)) then
      allocate (fvsav(plon,plev,beglat:endlat))  ! v wind tendency
      fvsav(:,:,:)=0._r8
   endif
  end subroutine init_iop_fields

subroutine readiopdata(timelevel)


!-----------------------------------------------------------------------
!     
!     Open and read netCDF file containing initial IOP  conditions
!     
!---------------------------Code history--------------------------------
!     
!     Written by J.  Truesdale    August, 1996, revised January, 1998
!     
!-----------------------------------------------------------------------
        use ppgrid,              only: begchunk, endchunk
	use phys_grid,           only: clat_p
	use commap,              only: latdeg, clat
        use getinterpnetcdfdata, only: getinterpncdata
        use shr_sys_mod,         only: shr_sys_flush
	use hycoef,              only: hyam, hybm
        use error_messages,      only: handle_ncerr
        use mo_chem_utls,        only : get_spc_ndx
        use mo_gas_phase_chemdr, only : map2chm
        use constituents,        only : sflxnam
!-----------------------------------------------------------------------
   implicit none
#if ( defined RS6000 )
   implicit automatic ( a-z )
#endif

   character(len=*), parameter ::  sub = "read_iop_data"

!------------------------------Input Arguments--------------------------
!     
integer, optional, intent(in) :: timelevel

!------------------------------Locals-----------------------------------
!     
   integer ntimelevel
   integer NCID, status
   integer time_dimID, lev_dimID,  lev_varID
   integer tsec_varID, bdate_varID,varid
   integer i,j
   integer nlev
   integer total_levs
   integer u_attlen

   integer bdate, ntime,nstep
   integer, allocatable :: tsec(:)
   integer k, m, n
   integer icldliq,icldice
   integer inumliq,inumice,idx

   logical have_srf              ! value at surface is available
   logical fill_ends             ! 
   logical swaplevs              ! 
   logical have_cnst(pcnst)
   real(r8) dummy
   real(r8) lat,xlat
   real(r8) srf(1)                  ! value at surface
   real(r8) pmid(plev)  ! pressure at model levels (time n)
   real(r8) pint(plevp) ! pressure at model interfaces (n  )
   real(r8) pdel(plev)  ! pdel(k)   = pint  (k+1)-pint  (k)
   real(r8) weight
   real(r8) tmpdata(1)
   real(r8) coldata(plev)
   real(r8), allocatable :: dplevs( : )
   integer strt4(4),cnt4(4),strt5(4),cnt5(4)
   character(len=16) :: lowername
   character(len=16) :: varname
   character(len=max_chars) :: units ! Units
   integer h2o_ndx

   nstep = get_nstep()
   fill_ends= .false.
   swaplevs = .false.

   if (present(timelevel)) then
      ntimelevel=timelevel
   else
      ntimelevel=n3
   end if

!     
!     Open IOP dataset
!     
  call handle_ncerr( nf90_open (iopfile, 0, ncid),&
       'readiopdata.F90', __LINE__)

!
!     if the dataset is a CAM generated dataset set use_camiop to true
!       CAM IOP datasets have a global attribute called CAM_GENERATED_IOP      
!
   if ( nf90_inquire_attribute( ncid, NF90_GLOBAL, 'CAM_GENERATED_FORCING', attnum=i )== NF90_NOERR ) then
      use_camiop = .true.
   else
      use_camiop = .false.
   endif

!=====================================================================
!     
!     Read time variables


   status = nf90_inq_dimid (ncid, 'time', time_dimID )
   if (status /= NF90_NOERR) status = nf90_inq_dimid (ncid, 'tsec', time_dimID )
   if (status /= NF90_NOERR) status = nf90_inq_dimid (ncid, 'time_offset', time_dimID )
   if (status /= NF90_NOERR) then
      if (masterproc) write(iulog,*) sub//':ERROR - readiopdata.F:Could not find dimension ID for time/tsec'
      status = NF90_CLOSE ( ncid )
      call endrun
   end if

   call handle_ncerr( nf90_inquire_dimension( ncid, time_dimID, len=ntime ),&
         'readiopdata.F90', __LINE__)

   allocate(tsec(ntime))

   status = nf90_inq_varid (ncid, 'tsec', tsec_varID )
   if (status /= NF90_NOERR) status = nf90_inq_varid (ncid, 'time_offset', tsec_varID )
      
   if (status /= NF90_NOERR) then
      if (masterproc) write(iulog,*) sub//':ERROR - readiopdata.F:Could not find variable ID for tsec'
      status = NF90_CLOSE ( ncid )
      call endrun
   else
      call handle_ncerr( nf90_get_var (ncid, tsec_varID, tsec),&
           'readiopdata.F90', __LINE__)
   end if

   status = nf90_inq_varid (ncid, 'nbdate', bdate_varID )
   if (status /= NF90_NOERR) status = nf90_inq_varid (ncid, 'bdate', bdate_varID )      
   if (status /= NF90_NOERR) then
      if (masterproc) write(iulog,*) sub//':ERROR - readiopdata.F:Could not find variable ID for bdate'
      status = NF90_CLOSE ( ncid )
      call endrun
   else
      call handle_ncerr( nf90_get_var (ncid, bdate_varID, bdate),&
           'readiopdata.F90', __LINE__)
   end if
!     
!======================================================
!     read level data
!     
   status = NF90_INQ_DIMID( ncid, 'lev', lev_dimID )
   if ( status .ne. nf90_noerr ) then
      if (masterproc) write(iulog,*) sub//':ERROR - readiopdata.F:Could not find variable dim ID  for lev'
      status = NF90_CLOSE ( ncid )
      return
   end if

   call handle_ncerr( nf90_inquire_dimension( ncid, lev_dimID, len=nlev ),&
         'readiopdata.F90', __LINE__)

   allocate(dplevs(nlev+1))

   status = NF90_INQ_VARID( ncid, 'lev', lev_varID )
   if ( status .ne. nf90_noerr ) then
      if (masterproc) write(iulog,*) sub//':ERROR - readiopdata.F:Could not find variable ID for lev'
      status = NF90_CLOSE ( ncid )
      return
   end if

   call handle_ncerr( nf90_get_var (ncid, lev_varID, dplevs(:nlev)),&
                    'readiopdata.F90', __LINE__)
   if (dplevs(1)-dplevs(nlev)>0._r8) then
      swaplevs=.true.
      dplevs(1:nlev)=dplevs(nlev:1:-1)
   end if
!
!CAM generated forcing already has pressure on millibars convert standard IOP if needed.
!
   call handle_ncerr(nf90_inquire_attribute(ncid, lev_varID, 'units', len=u_attlen),&
                    'readiopdata.F90', __LINE__)
   call handle_ncerr(nf90_get_att(ncid, lev_varID, 'units', units),&
                    'readiopdata.F90', __LINE__)
   units=trim(to_lower(units(1:u_attlen)))

   if ( units=='pa' .or. units=='pascal' .or. units=='pascals' ) then
!
!     convert pressure from Pascals to Millibars ( lev is expressed in pascals in iop datasets )
!
      do i=1,nlev
         dplevs( i ) = dplevs( i )/100._r8
      end do
   endif


   call shr_scam_GetCloseLatLon(ncid,scmlat,scmlon,closelat,closelon,closelatidx,closelonidx)

   lonid = 0
   latid = 0
   levid = 0
   timeid = 0

   call wrap_inq_dimid(ncid, 'lat', latid)
   call wrap_inq_dimid(ncid, 'lon', lonid)
   call wrap_inq_dimid(ncid, 'lev', levid)
   call wrap_inq_dimid(ncid, 'time', timeid)
 
   strt4(1) = closelonidx
   strt4(2) = closelatidx
   strt4(3) = iopTimeIdx
   strt4(4) = 1
   cnt4(1)  = 1
   cnt4(2)  = 1
   cnt4(3)  = 1
   cnt4(4)  = 1

   status = nf90_inq_varid( ncid, 'Ps', varid   )
   if ( status .ne. nf90_noerr ) status = nf90_inq_varid( ncid, 'p_srf_center', varid   )
   if ( status .ne. nf90_noerr ) then
      have_ps = .false.
      if (masterproc) write(iulog,*) sub//':Could not find variable Ps'
      if ( .not. scm_backfill_iop_w_init ) then
         status = NF90_CLOSE( ncid )
         return
      else
         if ( is_first_step() .and. masterproc) write(iulog,*) 'Using pressure value from Analysis Dataset'
      endif
   else
      status = nf90_get_var(ncid, varid, ps(1,1,ntimelevel), strt4)
      have_ps = .true.
      call handle_ncerr(nf90_inquire_attribute(ncid, varid, 'units', len=u_attlen),&
                    'readiopdata.F90', __LINE__)
      call handle_ncerr(nf90_get_att(ncid, varid, 'units', units),&
                    'readiopdata.F90', __LINE__)
      units=trim(to_lower(units(1:u_attlen)))
      ! convert ps to pascals if needed.
      if ( units=='mb'.or. units=='millibars' ) ps(1,1,ntimelevel)=ps(1,1,ntimelevel)*100._r8
   endif

!  If the IOP dataset has hyam,hybm,etc it is assumed to be a hybrid level
!  dataset.
	
   status =  nf90_inq_varid( ncid, 'hyam', varid   )
   if ( status == nf90_noerr ) then
      do i = 1, nlev
         dplevs( i ) = 1000.0_r8 * hyam( i ) + ps(1,1,ntimelevel) * hybm( i ) / 100.0_r8
      end do
   endif

!     add the surface pressure to the pressure level data, so that
!     surface boundary condition will be set properly,
!     making sure that it is the highest pressure in the array.
!

   total_levs = nlev+1
   dplevs(nlev+1) = ps(1,1,ntimelevel)/100.0_r8 ! ps is expressed in pascals
   do i= nlev, 1, -1
      if ( dplevs(i) > ps(1,1,ntimelevel)/100.0_r8) then
         total_levs = i
         dplevs(i) = ps(1,1,ntimelevel)/100.0_r8
      end if
   end do
   if (.not. use_camiop ) then
      nlev = total_levs
   endif
   if ( nlev == 1 ) then
      if (masterproc) write(iulog,*) sub//':Error - Readiopdata.F: Ps too low!'
      return
   endif

!=====================================================================


   status =  nf90_inq_varid( ncid, 'Tsair', varid   )
   if ( status .ne. nf90_noerr ) status =  nf90_inq_varid( ncid, 'T_srf', varid   )
   if ( status .ne. nf90_noerr ) then
      have_tsair = .false.
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,tsair)
      have_tsair = .true.
      call handle_ncerr(nf90_inquire_attribute(ncid, varid, 'units', len=u_attlen),&
           'readiopdata.F90', __LINE__)
      call handle_ncerr(nf90_get_att(ncid, varid, 'units', units),&
           'readiopdata.F90', __LINE__)
      units=trim(to_lower(units(1:u_attlen)))

      if ( units=='c' .or. units=='celsius' ) then
         !
         !     convert celsius to kelvin
         !
         tsair=tsair+273.15_r8
      endif
   endif

!
!      read in Tobs  For cam generated iop readin small t to avoid confusion
!      with capital T defined in cam
!

   tobs(:)= t3(1,:,1,ntimelevel)
   
   varname='t'
   status =  nf90_inq_varid( ncid, 't', varid   )
   if ( status .eq. nf90_noerr ) varname='t'
   status =  nf90_inq_varid( ncid, 'T', varid   )
   if ( status .eq. nf90_noerr ) varname='T'
   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx,trim(varname), have_tsair, &
        tsair(1), fill_ends, swaplevs, &
        dplevs, nlev,ps(1,1,ntimelevel),tobs, status )
   if ( status .ne. nf90_noerr ) then
      have_t = .false.
      if (masterproc) write(iulog,*) sub//':Could not find variable T'
      if ( .not. scm_backfill_iop_w_init ) then
         status = NF90_CLOSE( ncid )
         return
      else
         if (masterproc) write(iulog,*) sub//':Using value from Analysis Dataset'
      endif
!     
!     set T3 to Tobs on first time step
!     
   else
      have_t = .true.
      status =  nf90_inq_varid( ncid, trim(varname), varid   )
      call handle_ncerr(nf90_inquire_attribute(ncid, varid, 'units', len=u_attlen),&
           'readiopdata.F90', __LINE__)
      call handle_ncerr(nf90_get_att(ncid, varid, 'units', units),&
           'readiopdata.F90', __LINE__)
      units=trim(to_lower(units(1:u_attlen)))
      
      if ( units=='c' .or. units=='celsius' ) then
         !
         !     convert celsius to kelvin
         !
         tobs=tobs+273.15_r8
      endif
      
   endif

   status = nf90_inq_varid( ncid, 'Tg', varid   )
   if ( status .ne. nf90_noerr ) status =  nf90_inq_varid( ncid, 'T_skin', varid   )
   if (status .ne. nf90_noerr) then
      if (masterproc) write(iulog,*) sub//':Could not find variable Tg on IOP dataset'
      if ( have_tsair ) then
         if (masterproc) write(iulog,*) sub//':Using Tsair'
         tground = tsair     ! use surface value from T field
         have_Tg = .true.
      else
         have_Tg = .true.
         if (masterproc) write(iulog,*) sub//':Using T at lowest level from IOP dataset'
         tground = tobs(plev)
      endif
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,tground)
      have_Tg = .true.
      call handle_ncerr(nf90_inquire_attribute(ncid, varid, 'units', len=u_attlen),&
           'readiopdata.F90', __LINE__)
      call handle_ncerr(nf90_get_att(ncid, varid, 'units', units),&
           'readiopdata.F90', __LINE__)
      units=trim(to_lower(units(1:u_attlen)))
      
      if ( units=='c' .or. units=='celsius' ) then
         !
         !     convert celsius to kelvin
         !
         tground=tground+273.15_r8
      endif

   endif

   status = nf90_inq_varid( ncid, 'qsrf', varid   )
   if ( status .ne. nf90_noerr ) status =  nf90_inq_varid( ncid, 'q_srf', varid   )

   if ( status .ne. nf90_noerr ) then
      have_qref = .false.
   else
      status = nf90_get_var(ncid, varid, qrefobs(1), strt4)
      have_qref = .true.
      call handle_ncerr(nf90_inquire_attribute(ncid, varid, 'units', len=u_attlen),&
           'readiopdata.F90', __LINE__)
      call handle_ncerr(nf90_get_att(ncid, varid, 'units', units),&
           'readiopdata.F90', __LINE__)
      units=trim(to_lower(units(1:u_attlen)))
      
      if ( units=='g/kg' ) then
         !
         !     convert to kg/kg
         !
         qrefobs(1)=qrefobs(1)/1000._r8
      endif

   endif

   if (is_first_step()) then
      qinitobs(:,:)=q3(1,:,:,1,ntimelevel)
   end if

   qobs(:)= q3(1,:,1,1,ntimelevel)

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx,  'q', have_qref, &
      qrefobs(1), fill_ends, swaplevs, &
      dplevs, nlev,ps(1,1,ntimelevel), qobs, status )
   if ( status .ne. nf90_noerr ) then
      have_q = .false.
      if (masterproc) write(iulog,*) sub//':Could not find variable q'
      if ( .not. scm_backfill_iop_w_init ) then
         status = nf90_close( ncid )
         return
      else
         if (masterproc) write(iulog,*) sub//':Using values from Analysis Dataset'
      endif
   else
      have_q = .true.
      status =  nf90_inq_varid( ncid, 'q', varid   )
      call handle_ncerr(nf90_inquire_attribute(ncid, varid, 'units', len=u_attlen),&
           'readiopdata.F90', __LINE__)
      call handle_ncerr(nf90_get_att(ncid, varid, 'units', units),&
           'readiopdata.F90', __LINE__)
      units=trim(to_lower(units(1:u_attlen)))
      
      if ( units=='g/kg' ) then
         !
         !     convert to kg/kg
         !
         qobs=qobs/1000._r8
      endif

   endif

   cldobs = 0._r8
   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx,  'cld', .false., &
      dummy, fill_ends, swaplevs, dplevs, nlev,ps(1,1,ntimelevel), cldobs, status )
   if ( status .ne. nf90_noerr ) then
      have_cld = .false.
   else
      have_cld = .true.
   endif
   
   clwpobs = 0._r8
   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx,  'clwp', .false., &
      dummy, fill_ends, swaplevs, dplevs, nlev,ps(1,1,ntimelevel), clwpobs, status )
   if ( status .ne. nf90_noerr ) then
      have_clwp = .false.
   else
      have_clwp = .true.
   endif

!
!	read divq (horizontal advection)
!      
   status = nf90_inq_varid( ncid, 'divqsrf', varid   )
   if ( status .ne. nf90_noerr ) then
      have_srf = .false.
      srf=0._r8
   else
      status = nf90_get_var(ncid, varid, srf(1), strt4)
      have_srf = .true.
   endif

   divq(:,:)=0._r8
   varname='divq'
   status =  nf90_inq_varid( ncid, 'divq', varid   )
   if ( status .eq. nf90_noerr ) varname='divq'
   status =  nf90_inq_varid( ncid, 'q_adv_h', varid   )
   if ( status .eq. nf90_noerr ) varname='q_adv_h'

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, &
        trim(varname), have_srf, srf(1), fill_ends, swaplevs, &
        dplevs, nlev,ps(1,1,ntimelevel), divq(:,1), status )
   if ( status .ne. nf90_noerr ) then
      have_divq = .false.
   else
      have_divq = .true.
      status =  nf90_inq_varid( ncid, trim(varname), varid   )
      call handle_ncerr(nf90_inquire_attribute(ncid, varid, 'units', len=u_attlen),&
           'readiopdata.F90', __LINE__)
      call handle_ncerr(nf90_get_att(ncid, varid, 'units', units),&
           'readiopdata.F90', __LINE__)
      units=trim(to_lower(units(1:u_attlen)))
      select case (trim(units))
      case ('kg/kg/s') ! this is what we want
      case ('g/kg/s')
         divq(:,1)=divq(:,1)/1000._r8
      case ('kg/kg/h')
         divq(:,1)=divq(:,1)/3600._r8
      case ('kg/kg/hour')
         divq(:,1)=divq(:,1)/3600._r8
      case ('g/kg/hr')
         divq(:,1)=divq(:,1)/1000._r8/3600._r8
      case ('g/kg/hour')
         divq(:,1)=divq(:,1)/1000._r8/3600._r8
      case default
         call endrun(sub//':ERROR - no conversion for units of q horizontal advection term given:'//trim(units)//' need kg/kg/s')
      end select

   endif

!
!     read vertdivq if available
!
   status = nf90_inq_varid( ncid, 'vertdivqsrf', varid   )
   if ( status .ne. nf90_noerr ) then
      have_srf = .false.
      srf=0._r8
   else
      status = nf90_get_var(ncid, varid, srf(1), strt4)
      have_srf = .true.
   endif

   vertdivq=0._r8

   varname='vertdivq'
   status =  nf90_inq_varid( ncid, 'vertdivq', varid   )
   if ( status .eq. nf90_noerr ) varname='vertdivq'
   status =  nf90_inq_varid( ncid, 'q_adv_v', varid   )
   if ( status .eq. nf90_noerr ) varname='q_adv_v'
   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, trim(varname), &
        have_srf, srf(1), fill_ends, swaplevs, &
        dplevs, nlev,ps(1,1,ntimelevel), vertdivq(:,1), status )
   if ( status .ne. nf90_noerr ) then
      have_vertdivq = .false.
   else
      have_vertdivq = .true.
      status =  nf90_inq_varid( ncid, trim(varname), varid   )
      call handle_ncerr(nf90_inquire_attribute(ncid, varid, 'units', len=u_attlen),&
           'readiopdata.F90', __LINE__)
      call handle_ncerr(nf90_get_att(ncid, varid, 'units', units),&
           'readiopdata.F90', __LINE__)
      units=trim(to_lower(units(1:u_attlen)))
      select case (trim(units))
      case ('kg/kg/s') ! this is what we want
      case ('g/kg/s')
         vertdivq(:,1)=vertdivq(:,1)/1000._r8
      case ('kg/kg/hr')
         vertdivq(:,1)=vertdivq(:,1)/3600._r8
      case ('kg/kg/hour')
         vertdivq(:,1)=vertdivq(:,1)/3600._r8
      case ('g/kg/hr')
         vertdivq(:,1)=vertdivq(:,1)/1000._r8/3600._r8
      case ('g/kg/hour')
         vertdivq(:,1)=vertdivq(:,1)/1000._r8/3600._r8
      case default
         call endrun(sub//':ERROR - no conversion for units of q vertical advection term given:'//trim(units)//' need kg/kg/s')
      end select
   endif

!
!   add calls to get dynamics tendencies for all prognostic consts
!
   divq3d=0._r8
   srf=0._r8
   do m = 1, pcnst

      call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, trim(cnst_name(m))//'_dten', &
      .false., srf(1), fill_ends, swaplevs, &
      dplevs, nlev,ps(1,1,ntimelevel), divq3d(:,m), status )
      if ( status .ne. nf90_noerr ) then
         have_cnst(m) = .false.
         divq3d(1:,m)=0._r8
      else
         if (m==1) have_divq3d = .true.
         have_cnst(m) = .true.
      endif
      
      coldata = 0._r8
      call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, trim(cnst_name(m))//'_dqfx', &
      .false., srf(1), fill_ends, swaplevs, &
      dplevs, nlev,ps(1,1,ntimelevel), coldata, status )
      if ( STATUS .NE. NF90_NOERR ) then
         dqfxcam(1,:,m)=0._r8
      else
         dqfxcam(1,:,m)=coldata(:)
      endif

      tmpdata = 0._r8
      call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, trim(cnst_name(m))//'_alph', &
      .false., srf(1), fill_ends, swaplevs, &
      dplevs, nlev,ps(1,1,ntimelevel), tmpdata, status )
      if ( status .ne. nf90_noerr ) then
!         have_cnst(m) = .false.
         alphacam(m)=0._r8
      else
          alphacam(m)=tmpdata(1)
!         have_cnst(m) = .true.
      endif

   end do
!  read surface emissions
   h2o_ndx   = get_spc_ndx('H2O')
   do m = 1,pcnst
      n = map2chm(m)
      if ( n /= h2o_ndx .and. n > 0 ) then
         tmpdata = 0._r8
         status =  nf90_inq_varid( ncid, trim(sflxnam(m)), varid   )
         if ( status .ne. nf90_noerr ) then
            have_sflx(m) = .false.
         else
            status = nf90_get_var(ncid, varid, srf(1), strt4)
            have_sflx(m) = .true.
            sflxobs(m)=srf(1)
         endif
      endif
   end do

   numliqobs = 0._r8
   call cnst_get_ind('NUMLIQ', inumliq, abort=.false.)
   if ( inumliq > 0 ) then
      have_srf = .false.
      srf=0._r8
      call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'NUMLIQ', &
           have_srf, srf(1), fill_ends, swaplevs, &
           dplevs, nlev,ps(1,1,ntimelevel), numliqobs, status )
      if ( status .ne. nf90_noerr ) then
         have_numliq = .false.
      else
         have_numliq = .true.
         do i=1, PLEV
            q3(1,i,inumliq,1,ntimelevel)=numliqobs(i)
         end do
      endif
   else
         have_numliq = .false.
   end if

   have_srf = .false.
   srf=0._r8
   cldliqobs = 0._r8
   call cnst_get_ind('CLDLIQ', icldliq, abort=.false.)
   if ( icldliq > 0 ) then
      call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'CLDLIQ', &
           have_srf, srf(1), fill_ends, swaplevs, &
           dplevs, nlev,ps(1,1,ntimelevel), cldliqobs, status )
      if ( status .ne. nf90_noerr ) then
         have_cldliq = .false.
      else
         have_cldliq = .true.
         do i=1, PLEV
            q3(1,i,icldliq,1,ntimelevel)=cldliqobs(i)
         end do
      endif
   else
         have_cldliq = .false.
   endif

   have_srf = .false.
   srf=0._r8
   cldiceobs = 0._r8
   call cnst_get_ind('CLDICE', icldice, abort=.false.)
   if ( icldice > 0 ) then
      call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'CLDICE', &
           have_srf, srf(1), fill_ends, swaplevs, &
           dplevs, nlev,ps(1,1,ntimelevel), cldiceobs, status )
      if ( status .ne. nf90_noerr ) then
         have_cldice = .false.
      else
         have_cldice = .true.
         do i=1, PLEV
            q3(1,i,icldice,1,ntimelevel)=cldiceobs(i)
         end do
      endif
   else
      have_cldice = .false.
   endif

   have_srf = .false.
   srf=0._r8
   numiceobs = 0._r8
   call cnst_get_ind('NUMICE', inumice, abort=.false.)
   if ( inumice > 0 ) then
      have_srf = .false.
      call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'NUMICE', &
         have_srf, srf(1), fill_ends, swaplevs, &
         dplevs, nlev,ps(1,1,ntimelevel), numiceobs, status )
      if ( status .ne. nf90_noerr ) then
         have_numice = .false.
      else
         have_numice = .true.
         do i=1, PLEV
            q3(1,i,inumice,1,ntimelevel)=numiceobs(i)
         end do
      endif
   else
      have_numice = .false.
   end if

!
!	read divu (optional field)
!      
   status = nf90_inq_varid( ncid, 'divusrf', varid   )
   if ( status .ne. nf90_noerr ) then
      have_srf = .false.
      srf=0._r8
   else
      status = nf90_get_var(ncid, varid, srf(1), strt4)
      have_srf = .true.
   endif

   have_srf = .false.
   srf=0._r8
   divu = 0._r8
   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'divu', &
      have_srf, srf(1), fill_ends, swaplevs, &
      dplevs, nlev,ps(1,1,ntimelevel), divu, status )
   if ( status .ne. nf90_noerr ) then
      have_divu = .false.
   else
      have_divu = .true.
   endif
!
!	read divv (optional field)
!      
   status = nf90_inq_varid( ncid, 'divvsrf', varid   )
   if ( status .ne. nf90_noerr ) then
      have_srf = .false.
   else
      status = nf90_get_var(ncid, varid, srf(1), strt4)
      have_srf = .true.
   endif

   divv = 0._r8
   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'divv', &
      have_srf, srf(1), fill_ends, swaplevs, &
      dplevs, nlev,ps(1,1,ntimelevel), divv, status )
   if ( status .ne. nf90_noerr ) then
      have_divv = .false.
   else
      have_divv = .true.
   endif
!
!	read divt (optional field)
!      
   status = nf90_inq_varid( ncid, 'divtsrf', varid   )
   if ( status .ne. nf90_noerr ) then
      have_srf = .false.
   else
      status = nf90_get_var(ncid, varid, srf(1), strt4)
      have_srf = .true.
   endif

   divt=0._r8
   varname='divT'
   status =  nf90_inq_varid( ncid, 'divT', varid   )
   if ( status .eq. nf90_noerr ) varname='divT'
   status =  nf90_inq_varid( ncid, 'T_adv_h', varid   )
   if ( status .eq. nf90_noerr ) varname='T_adv_h'

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, &
      trim(varname), have_srf, srf(1), fill_ends, swaplevs, &
      dplevs, nlev,ps(1,1,ntimelevel), divt, status )
   if ( status .ne. nf90_noerr ) then
      have_divt = .false.
   else
      have_divt = .true.
      status =  nf90_inq_varid( ncid, trim(varname), varid   )
      call handle_ncerr(nf90_inquire_attribute(ncid, varid, 'units', len=u_attlen),&
           'readiopdata.F90', __LINE__)
      call handle_ncerr(nf90_get_att(ncid, varid, 'units', units),&
           'readiopdata.F90', __LINE__)
      units=trim(to_lower(units(1:u_attlen)))
      select case (trim(units))
      case ('k/s') ! this is what we want
      case ('k/hr')
         divt=divt/3600._r8
      case ('k/hour')
         divt=divt/3600._r8
      case ('c/s') 
         divt=divt+273.15_r8
      case ('c/hr')
         divt=divt+273.15_r8
         divt=divt/3600._r8
      case ('c/hour') 
         divt=divt+273.15_r8
         divt=divt/3600._r8
      case default
         call endrun(sub//':ERROR - no conversion for units of T horizontal advection term given:'//trim(units)//' need k/s')
      end select
   endif

!
!     read vertdivt if available
!
   status = nf90_inq_varid( ncid, 'vertdivTsrf', varid   )
   if ( status .ne. nf90_noerr ) then
      have_srf = .false.
   else
      status = nf90_get_var(ncid, varid, srf(1), strt4)
      have_srf = .true.
   endif
   
   vertdivt=0._r8
   varname='vertdivT'
   status =  nf90_inq_varid( ncid, 'vertdivT', varid   )
   if ( status .eq. nf90_noerr ) varname='vertdivT'
   status =  nf90_inq_varid( ncid, 'T_adv_v', varid   )
   if ( status .eq. nf90_noerr ) varname='T_adv_v'

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, trim(varname), &
      have_srf, srf(1), fill_ends, swaplevs, &
      dplevs, nlev,ps(1,1,ntimelevel), vertdivt, status )
   if ( status .ne. nf90_noerr ) then
      have_vertdivt = .false.
   else
      have_vertdivt = .true.
      status =  nf90_inq_varid( ncid, trim(varname), varid   )
      call handle_ncerr(nf90_inquire_attribute(ncid, varid, 'units', len=u_attlen),&
           'readiopdata.F90', __LINE__)
      call handle_ncerr(nf90_get_att(ncid, varid, 'units', units),&
           'readiopdata.F90', __LINE__)
      units=trim(to_lower(units(1:u_attlen)))
      select case (trim(units))
      case ('k/s') ! this is what we want
      case ('k/hr')
         vertdivt=vertdivt/3600._r8
      case ('k/hour')
         vertdivt=vertdivt/3600._r8
      case ('c/s') 
         vertdivt=vertdivt+273.15_r8
      case ('c/hr')
         vertdivt=vertdivt+273.15_r8
         vertdivt=vertdivt/3600._r8
      case ('c/hour')
         vertdivt=vertdivt+273.15_r8
         vertdivt=vertdivt/3600._r8
      case default
         call endrun(sub//':ERROR - no conversion for units of T vertical advection term given:'//trim(units)//' need k/s')
      end select
   endif
!
!	read divt3d (combined vertical/horizontal advection)
!      (optional field)

   status = nf90_inq_varid( ncid, 'divT3dsrf', varid   )
   if ( status .ne. nf90_noerr ) then
      have_srf = .false.
   else
      status = nf90_get_var(ncid, varid, srf(1), strt4)
      have_srf = .true.
   endif

   divT3d = 0._r8

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'divT3d', &
      have_srf, srf(1), fill_ends, swaplevs, &
      dplevs, nlev,ps(1,1,ntimelevel), divt3d, status )
   if ( status .ne. nf90_noerr ) then
      have_divt3d = .false.
   else
      have_divt3d = .true.
   endif

   divU3d = 0._r8

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'divU3d', &
      have_srf, srf(1), fill_ends, swaplevs, &
      dplevs, nlev,ps(1,1,ntimelevel), divu3d, status )
   if ( status .ne. nf90_noerr ) then
      have_divu3d = .false.
   else
      have_divu3d = .true.
   endif

   divV3d = 0._r8

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'divV3d', &
      have_srf, srf(1), fill_ends, swaplevs, &
      dplevs, nlev,ps(1,1,ntimelevel), divv3d, status )
   if ( status .ne. nf90_noerr ) then
      have_divv3d = .false.
   else
      have_divv3d = .true.
   endif
   
   status = nf90_inq_varid( ncid, 'Ptend', varid   )
   if ( status .ne. nf90_noerr ) status = nf90_inq_varid( ncid, 'omega_srf', varid   )
   if ( status .ne. nf90_noerr ) then
      have_ptend = .false.
      if (masterproc) write(iulog,*) sub//':Could not find variable Ptend. Setting to zero'
      ptend = 0.0_r8
   else
      status = nf90_get_var(ncid, varid, srf(1), strt4)
      have_ptend = .true.
      ptend= srf(1)
      call handle_ncerr(nf90_inquire_attribute(ncid, varid, 'units', len=u_attlen),&
           'readiopdata.F90', __LINE__)
      call handle_ncerr(nf90_get_att(ncid, varid, 'units', units),&
           'readiopdata.F90', __LINE__)
      units=trim(to_lower(units(1:u_attlen)))
      select case (trim(units))
      case ('pa/s') ! this is what we want
      case ('pascals/second') ! this is what we want
      case ('pa/hr')
         ptend=ptend/3600._r8
      case ('pa/hour')
         ptend=ptend/3600._r8
      case ('mb/s')
         ptend=ptend*100._r8
      case ('millibars/second')
         ptend=ptend*100._r8
      case ('mb/hr')
         ptend=ptend*100._r8/3600._r8
      case ('mb/hour')
         ptend=ptend*100._r8/3600._r8
      case ('millibars/hour')
         ptend=ptend*100._r8/3600._r8
      case default
         call endrun(sub//':ERROR - no conversion for units of pressure tendancy term given:'//trim(units)//' need pa/s')
      end select
   endif

   wfld=0._r8

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, &
      'omega', .true., ptend, fill_ends, swaplevs, &
      dplevs, nlev,ps(1,1,ntimelevel), wfld, status )
   if ( status .ne. nf90_noerr ) then
      have_omega = .false.
      if (masterproc) write(iulog,*) sub//':Could not find variable omega'
      if ( .not. scm_backfill_iop_w_init ) then
         status = nf90_close( ncid )
         return
      else
         if (masterproc) write(iulog,*) sub//'Using value from Analysis Dataset'
      endif
   else
      have_omega = .true.
      status =  nf90_inq_varid( ncid, 'omega', varid   )
      call handle_ncerr(nf90_inquire_attribute(ncid, varid, 'units', len=u_attlen),&
           'readiopdata.F90', __LINE__)
      call handle_ncerr(nf90_get_att(ncid, varid, 'units', units),&
           'readiopdata.F90', __LINE__)
      units=trim(to_lower(units(1:u_attlen)))
      select case (trim(units))
      case ('pa/s') ! this is what we want
      case ('pascals/second') ! this is what we want
      case ('pa/hr')
         wfld=wfld/3600._r8
      case ('pa/hour')
         wfld=wfld/3600._r8
      case ('mb/s')
         wfld=wfld*100._r8
      case ('millibars/second')
         wfld=wfld*100._r8
      case ('mb/hr')
         wfld=wfld*100._r8/3600._r8
      case ('mb/hour')
         wfld=wfld*100._r8/3600._r8
      case ('millibars/hour')
         wfld=wfld*100._r8/3600._r8
      case default
         call endrun(sub//':ERROR - no conversion for units of omega term given:'//trim(units)//' need pa/s')
      end select
   endif
   call plevs0(1    ,plon   ,plev    ,ps(1,1,ntimelevel)   ,pint,pmid ,pdel)
   call shr_sys_flush( iulog )
!
! Build interface vector for the specified omega profile
! (weighted average in pressure of specified level values)
!
   wfldh(:) = 0.0_r8

   do k=2,plev
      weight = (pint(k) - pmid(k-1))/(pmid(k) - pmid(k-1))
      wfldh(k) = (1.0_r8 - weight)*wfld(k-1) + weight*wfld(k)
   end do

   status = nf90_inq_varid( ncid, 'usrf', varid   )
   if ( status .ne. nf90_noerr ) status = nf90_inq_varid( ncid, 'u_srf', varid   )
   if ( status .ne. nf90_noerr ) then
      have_uref_clasp = .false.
      uref_clasp(1)=0._r8
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,uref_clasp)
      have_uref_clasp = .true.
   endif

   uobs=0._r8

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, &
      'u', have_uref_clasp, uref_clasp(1), fill_ends, swaplevs, &
      dplevs, nlev,ps(1,1,ntimelevel), uobs, status )
   if ( status .ne. nf90_noerr ) then
      have_u = .false.
   else
      have_u = .true.
      do i=1, PLEV
         u3(1,i,1,ntimelevel) = uobs(i)  !     set u to uobs at first time step
      end do
   endif

   status = nf90_inq_varid( ncid, 'vsrf', varid   )
   if ( status .ne. nf90_noerr ) status = nf90_inq_varid( ncid, 'v_srf', varid   )
   if ( status .ne. nf90_noerr ) then
      have_vref_clasp = .false.
      vref_clasp = 0._r8
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,vref_clasp)
      have_vref_clasp = .true.
   endif

   vobs=0._r8

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, &
      'v', have_vref_clasp, vref_clasp(1), fill_ends, swaplevs, &
      dplevs, nlev,ps(1,1,ntimelevel), vobs, status )
   if ( status .ne. nf90_noerr ) then
      have_v = .false.
   else
      have_v = .true.
      do i=1, PLEV
         v3(1,i,1,ntimelevel) = vobs(i)  !     set u to uobs at first time step
      end do
   endif
   call shr_sys_flush( iulog )

   status = nf90_inq_varid( ncid, 'Prec', varid   )
   if ( status .ne. nf90_noerr ) status = nf90_inq_varid( ncid, 'prec_srf', varid   )
   if ( status .ne. nf90_noerr ) then
      have_prec = .false.
      precobs = 0._r8
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,precobs)
      have_prec = .true.
      call handle_ncerr(nf90_inquire_attribute(ncid, varid, 'units', len=u_attlen),&
           'readiopdata.F90', __LINE__)
      call handle_ncerr(nf90_get_att(ncid, varid, 'units', units),&
           'readiopdata.F90', __LINE__)
      units=trim(to_lower(units(1:u_attlen)))
      select case (trim(units))
      case ('m/s') ! this is what we want
      case ('meters/second') ! this is what we want
      case ('m/hr') 
         precobs=precobs/3600._r8
      case ('m/hour') 
         precobs=precobs/3600._r8
      case ('mm/s') 
         precobs=precobs/1000._r8
      case ('mm/hr') 
         precobs=precobs/1000._r8/3600._r8
      case ('mm/hour') 
         precobs=precobs/1000._r8/3600._r8
      case default
         call endrun(sub//':ERROR - no conversion for units of precipitation given:'//trim(units)//' need m/s')
      end select
   endif

   q1obs = 0._r8
   varname='Q1'
   status = nf90_inq_varid( ncid, 'Q1', varid   )
   if ( status .eq. nf90_noerr ) varname='Q1'
   status = nf90_inq_varid( ncid, 'q1', varid   )
   if ( status .eq. nf90_noerr ) varname='q1'

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, trim(varname), &
      .false., dummy, fill_ends, swaplevs, & ! datasets don't contain Q1 at surface
      dplevs, nlev,ps(1,1,ntimelevel), q1obs, status )
   if ( status .ne. nf90_noerr ) then
      have_q1 = .false.
   else
      have_q1 = .true.
      status = nf90_inq_varid( ncid, trim(varname), varid   )
      call handle_ncerr(nf90_inquire_attribute(ncid, varid, 'units', len=u_attlen),&
           'readiopdata.F90', __LINE__)
      call handle_ncerr(nf90_get_att(ncid, varid, 'units', units),&
           'readiopdata.F90', __LINE__)
      units=trim(to_lower(units(1:u_attlen)))
      select case (trim(units))
      case ('k/s') ! this is what we want
      case ('k/hr') 
         q1obs=q1obs/3600._r8
      case ('k/hour') 
         q1obs=q1obs/3600._r8
      case default
         call endrun(sub//':ERROR - no conversion for units of Q1 given:'//trim(units)//' need k/s')
      end select
   endif

   q2obs = 0._r8

   varname='Q2'
   status = nf90_inq_varid( ncid, 'Q2', varid   )
   if ( status .eq. nf90_noerr ) varname='Q2'
   status = nf90_inq_varid( ncid, 'q2', varid   )
   if ( status .eq. nf90_noerr ) varname='q2'
   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, trim(varname), &
      .false., dummy, fill_ends, swaplevs, & ! datasets don't contain Q2 at surface
      dplevs, nlev,ps(1,1,ntimelevel), q2obs, status )
   if ( status .ne. nf90_noerr ) then
      have_q2 = .false.
   else
      have_q2 = .true.
      status = nf90_inq_varid( ncid, trim(varname), varid   )
      call handle_ncerr(nf90_inquire_attribute(ncid, varid, 'units', len=u_attlen),&
           'readiopdata.F90', __LINE__)
      call handle_ncerr(nf90_get_att(ncid, varid, 'units', units),&
           'readiopdata.F90', __LINE__)
      units=trim(to_lower(units(1:u_attlen)))
      select case (trim(units))
      case ('k/s') ! this is what we want
      case ('k/hr') 
         q2obs=q2obs/3600._r8
      case ('k/hour') 
         q2obs=q2obs/3600._r8
      case default
         call endrun(sub//':ERROR - no conversion for units of Q2 given:'//trim(units)//' need k/s')
      end select
   endif

   status = nf90_inq_varid( ncid, 'lhflx', varid   )
   if ( status .ne. nf90_noerr ) status = nf90_inq_varid( ncid, 'lh', varid   )
   if ( status .ne. nf90_noerr ) status = nf90_inq_varid( ncid, 'LH', varid   )
   
   if ( status .ne. nf90_noerr ) then
      have_lhflx = .false.
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,lhflxobs)
      have_lhflx = .true.
   endif


   status = nf90_inq_varid( ncid, 'shflx', varid   )
   if ( status .ne. nf90_noerr ) status = nf90_inq_varid( ncid, 'sh', varid   )
   if ( status .ne. nf90_noerr ) status = nf90_inq_varid( ncid, 'SH', varid   )
   
   if ( status .ne. nf90_noerr ) then
      have_shflx = .false.
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,shflxobs)
      have_shflx = .true.
   endif

   ! Read in new land fluxes

   ! read wpqtp_sfc (kg\\kg)m\\s)
   status = nf90_inq_varid( ncid, 'wpqtp_sfc', varid   )
   if ( status .ne. nf90_noerr ) then
      have_wpqtpsfc_clasp = .false.
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,wpqtpsfc_clasp)
      have_wpqtpsfc_clasp = .true.
   endif

   ! read wpthlp_sfc (mK\\s)
   status = nf90_inq_varid( ncid, 'wpthlp_sfc', varid   )
   if ( status .ne. nf90_noerr ) then
      have_wpthlpsfc_clasp = .false.
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,wpthlpsfc_clasp)
      have_wpthlpsfc_clasp = .true.
   endif

   ! read qp2 (kg2\\kg2)
   status = nf90_inq_varid( ncid, 'qp2', varid   )
   if ( status .ne. nf90_noerr ) then
      have_qp2 = .false.
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,qp2obs)
      have_qp2 = .true.
      ! MDF: use CLASP variables as well 
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,rtp2_clasp)
      have_rtp2_clasp = .true.
   endif

   ! read thlp2 (K2)
   status = nf90_inq_varid( ncid, 'thlp2', varid   )
   if ( status .ne. nf90_noerr ) then
      have_thlp2_clasp = .false.
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,thlp2_clasp)
      have_thlp2_clasp = .true.
      ! MDF: Write statment to log file
      write(iulog,*) 'Using CLASP thlp2 values.'
      ! End MDF addition
   endif

   ! read thlpqp (Kkg\\kg)
   status = nf90_inq_varid( ncid, 'thlpqp', varid   )
   if ( status .ne. nf90_noerr ) then
      have_rtpthlp_clasp = .false.
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,rtpthlp_clasp)
      have_rtpthlp_clasp = .true.
   endif

   ! read wp2 (Kkg\\kg)
   status = nf90_inq_varid( ncid, 'wp2', varid   )
   if ( status .ne. nf90_noerr ) then
      have_wp2_clasp = .false.
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,wp2_clasp)
      have_wp2_clasp = .true.
   endif

   ! read wp4 (Kkg\\kg)
   status = nf90_inq_varid( ncid, 'wp4', varid   )
   if ( status .ne. nf90_noerr ) then
      have_wp4_clasp = .false.
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,wp4_clasp)
      have_wp4_clasp = .true.
   endif

   ! read wp2thetap (Km2\\s2)
   status = nf90_inq_varid( ncid, 'wp2thetap', varid   )
   if ( status .ne. nf90_noerr ) then
      have_wp2thlp_clasp = .false.
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,wp2thlp_clasp)
      have_wp2thlp_clasp = .true.
      ! MDF: Write statment to log file
      write(iulog,*) 'Using CLASP wp2thlp values.'
      ! End MDF addition
   endif

   ! read wp2qp (kg\\kg)(m2\\s2)
   status = nf90_inq_varid( ncid, 'wp2qp', varid   )
   if ( status .ne. nf90_noerr ) then
      have_wp2rtp_clasp = .false.
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,wp2rtp_clasp)
      have_wp2rtp_clasp = .true.
   endif

   ! read wpqp2 (kg2\\kg2)(m\\s)
   status = nf90_inq_varid( ncid, 'wpqp2', varid   )
   if ( status .ne. nf90_noerr ) then
      have_wprtp2_clasp = .false.
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,wprtp2_clasp)
      have_wprtp2_clasp = .true.
      ! MDF: Write statment to log file
      write(iulog,*) 'Using CLASP wpqp2 values.'
      ! End MDF addition
   endif

   ! read wpthetap2 (K2\\K2)(m\\s)
   status = nf90_inq_varid( ncid, 'wpthetap2', varid   )
   if ( status .ne. nf90_noerr ) then
      have_wpthlp2_clasp = .false.
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,wpthlp2_clasp)
      have_wpthlp2_clasp = .true.
      ! MDF: Write statment to log file
      write(iulog,*) 'Using CLASP wpthlp2 values.'
      ! End MDF addition
   endif

   ! read wpthetapqp (Kmkg\\skg)
   status = nf90_inq_varid( ncid, 'wpthetapqp', varid   )
   if ( status .ne. nf90_noerr ) then
      have_wprtpthlp_clasp = .false.
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,wprtpthlp_clasp)
      have_wprtpthlp_clasp = .true.
   endif

   ! read wp3 (m3\\s3)
   status = nf90_inq_varid( ncid, 'wp3', varid   )
   if ( status .ne. nf90_noerr ) then
      have_wp3_clasp = .false.
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,wp3_clasp)
      have_wp3_clasp = .true.
   endif

   ! read upwp (m2\\s2)
   status = nf90_inq_varid( ncid, 'upwp', varid   )
   if ( status .ne. nf90_noerr ) then
      have_upwp_clasp = .false.
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,upwp_clasp)
      have_upwp_clasp = .true.
   endif

   ! read vpwp (m2\\s2)
   status = nf90_inq_varid( ncid, 'vpwp', varid   )
   if ( status .ne. nf90_noerr ) then
      have_vpwp_clasp = .false.
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,vpwp_clasp)
      have_vpwp_clasp = .true.
   endif

   
   ! read T_soil
   status = nf90_inq_varid( ncid, 'T_soil', varid   )
   if ( status .ne. nf90_noerr ) then
      have_tsoil = .false.
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,tsoilobs)
      have_tsoil = .true.
      if ( units=='c' .or. units=='celsius' ) then
         !
         !     convert celsius to kelvin
         !
         tsoilobs=tsoilobs+273.15_r8
      endif
   endif

   ! read evap
   status = nf90_inq_varid( ncid, 'evap_srf', varid   )
   if ( status .ne. nf90_noerr ) then
      have_evap = .false.
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,evapobs)
      have_evap = .true.
      status = nf90_inq_varid( ncid, 'evap_srf', varid   )
      call handle_ncerr(nf90_inquire_attribute(ncid, varid, 'units', len=u_attlen),&
           'readiopdata.F90', __LINE__)
      call handle_ncerr(nf90_get_att(ncid, varid, 'units', units),&
           'readiopdata.F90', __LINE__)
      units=trim(to_lower(units(1:u_attlen)))
      select case (trim(units))
      case ('m/s') ! this is what we want
      case ('meters/second') ! this is what we want
      case ('m/hr') 
         evapobs=evapobs/3600._r8
      case ('m/hour') 
         evapobs=evapobs/3600._r8
      case ('mm/s') 
         evapobs=evapobs/1000._r8
      case ('mm/hr') 
         evapobs=evapobs/1000._r8/3600._r8
      case ('mm/hour') 
         evapobs=evapobs/1000._r8/3600._r8
      case default
         call endrun(sub//':ERROR - no conversion for units of evaporation given:'//trim(units)//' need m/s')
      end select

   endif

   ! read lw_up_srf
   status = nf90_inq_varid( ncid, 'lw_up_srf', varid   )
   if ( status .ne. nf90_noerr ) then
      have_lwupsrf = .false.
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,lwupsrfobs)
      have_lwupsrf = .true.
   endif

   ! read lw_dn_srf
   status = nf90_inq_varid( ncid, 'lw_dn_srf', varid   )
   if ( status .ne. nf90_noerr ) then
      have_lwdnsrf = .false.
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,lwdnsrfobs)
      have_lwdnsrf = .true.
   endif

   ! read sw_up_srf
   status = nf90_inq_varid( ncid, 'sw_up_srf', varid   )
   if ( status .ne. nf90_noerr ) then
      have_swupsrf = .false.
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,swupsrfobs)
      have_swupsrf = .true.
   endif

   ! read sw_dn_srf
   status = nf90_inq_varid( ncid, 'sw_dn_srf', varid   )
   if ( status .ne. nf90_noerr ) then
      have_swdnsrf = .false.
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,swdnsrfobs)
      have_swdnsrf = .true.
   endif

   call shr_sys_flush( iulog )

!
!     fill in 3d forcing variables if we have both horizontal
!     and vertical components, but not the 3d
!
   if ( .not. have_cnst(1) .and. have_divq .and. have_vertdivq ) then
      do k=1,plev
         do m=1,pcnst
            divq3d(k,m) = divq(k,m) + vertdivq(k,m)
         enddo
      enddo
      have_divq3d = .true.
   endif

   if ( .not. have_divt3d .and. have_divt .and. have_vertdivt ) then
      if (masterproc) write(iulog,*) sub//'Don''t have divt3d - using divt and vertdivt'
      do k=1,plev
         divt3d(k) = divt(k) + vertdivt(k)
      enddo
      have_divt3d = .true.
   endif
!
!     make sure that use_3dfrc flag is set to true if we only have
!     3d forcing available
!
   if ( .not. have_divt .or. .not. have_divq ) then
      use_3dfrc = .true.
   endif
   call shr_sys_flush( iulog )

   status =  nf90_inq_varid( ncid, 'CLAT', varid   )
   if ( status == nf90_noerr ) then
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,clat)
      clat_p(1)=clat(1)
      latdeg(1) = clat(1)*45._r8/atan(1._r8)
   endif

   status =  nf90_inq_varid( ncid, 'beta', varid   )
   if ( status .ne. nf90_noerr ) then
      betacam = 0._r8
   else
      status = nf90_get_var(ncid, varid, srf(1), strt4)
      betacam=srf(1) 
   endif

   status =  nf90_inq_varid( ncid, 'fixmas', varid   )
   if ( status .ne. nf90_noerr ) then
      fixmascam=1.0_r8
   else
      status = nf90_get_var(ncid, varid, srf(1), strt4)
      fixmascam=srf(1)
   endif

   call shr_sys_flush( iulog )

   status = nf90_close( ncid )
   call shr_sys_flush( iulog )

   deallocate(dplevs,tsec)

   return
end subroutine readiopdata

subroutine setiopupdate

!-----------------------------------------------------------------------
!   
! Open and read netCDF file to extract time information
!
!---------------------------Code history--------------------------------
!
! Written by John Truesdale    August, 1996
! 
!-----------------------------------------------------------------------
  implicit none
#if ( defined RS6000 )
  implicit automatic (a-z)
#endif
   character(len=*), parameter ::  sub = "setiopupdate"

!------------------------------Locals-----------------------------------

   integer NCID,i
   integer tsec_varID, time_dimID
   integer, allocatable :: tsec(:)
   integer  ntime 
   integer bdate, bdate_varID
   integer STATUS
   integer next_date, next_sec, last_date, last_sec 
   integer :: ncsec,ncdate                      ! current time of day,date
   integer :: yr, mon, day                      ! year, month, and day component
   integer :: start_ymd,start_tod
   save tsec, ntime, bdate
   save last_date, last_sec 
!------------------------------------------------------------------------------

   if ( is_first_step() ) then
!     
!     Open  IOP dataset
!     
      STATUS = NF90_OPEN( iopfile, NF90_NOWRITE, NCID )
!     
!     Read time (tsec) variable 
!     
      STATUS = NF90_INQ_VARID( NCID, 'tsec', tsec_varID )
      if ( STATUS .NE. NF90_NOERR) STATUS = NF90_INQ_VARID( NCID, 'time_offset', tsec_varID )
      if ( STATUS .NE. NF90_NOERR .and. masterproc) write(iulog,*) &
         sub//':ERROR - setiopupdate.F:', &
         'Cant get variable ID for tsec'

      STATUS = NF90_INQ_VARID( NCID, 'bdate', bdate_varID )
      if ( STATUS .NE. NF90_NOERR ) then
         STATUS = NF90_INQ_VARID( NCID, 'basedate', bdate_varID )
         if ( STATUS .NE. NF90_NOERR .and. masterproc) write(iulog,*) &
              sub//':ERROR - setiopupdate.F:Cant get variable ID for bdate'
      endif

      STATUS = NF90_INQ_DIMID( NCID, 'time', time_dimID )
      if ( STATUS .NE. NF90_NOERR )  then
         STATUS = NF90_INQ_DIMID( NCID, 'tsec', time_dimID )
         if ( STATUS .NE. NF90_NOERR )  then
            write(iulog,* )'ERROR - setiopupdate.F:Could not find variable dim ID for time'
            STATUS = NF90_CLOSE ( NCID )
            return
         end if
      end if

      if ( STATUS .NE. NF90_NOERR .and. masterproc) write(iulog,*) &
           sub//':ERROR - setiopupdate.F:Cant get variable dim ID for time'

      STATUS = NF90_INQUIRE_DIMENSION( NCID, time_dimID, len=ntime )
      if ( STATUS .NE. NF90_NOERR ) then
         if (masterproc) write(iulog,*) sub//':ERROR - setiopupdate.F:Cant get time dimlen'
      endif

      if (.not.allocated(tsec)) allocate(tsec(ntime))

      STATUS = NF90_GET_VAR( NCID, tsec_varID, tsec )
      if ( STATUS .NE. NF90_NOERR )then
         if (masterproc) write(iulog,*) sub//':ERROR - setiopupdate.F:Cant get variable tsec'
      endif
      STATUS = NF90_GET_VAR( NCID, bdate_varID, bdate )
      if ( STATUS .NE. NF90_NOERR )then
         if (masterproc) write(iulog,*) sub//':ERROR - setiopupdate.F:Cant get variable bdate'
      endif
!     Close the netCDF file
      STATUS = NF90_CLOSE( NCID )
!     
!     determine the last date in the iop dataset
!     
      call timemgr_time_inc(bdate, 0, last_date, last_sec, inc_s=tsec(ntime))
!     
!     set the iop dataset index
!    
      iopTimeIdx=0
      do i=1,ntime           ! set the first ioptimeidx
         call timemgr_time_inc(bdate, 0, next_date, next_sec, inc_s=tsec(i))
         call get_start_date(yr,mon,day,start_tod)
         start_ymd = yr*10000 + mon*100 + day

         if ( start_ymd > next_date .or. (start_ymd == next_date &
            .and. start_tod >= next_sec)) then
            iopTimeIdx = i
         endif
      enddo

      call get_curr_date(yr,mon,day,ncsec)
      ncdate=yr*10000 + mon*100 + day

      if (iopTimeIdx == 0.or.iopTimeIdx >= ntime) then
         call timemgr_time_inc(bdate, 0, next_date, next_sec, inc_s=tsec(1))
         if (masterproc) then
            write(iulog,*) 'Error::setiopupdate: Current model time does not fall within IOP period'
            write(iulog,*) ' Current CAM Date is ',ncdate,' and ',ncsec,' seconds'
            write(iulog,*) ' IOP start is        ',next_date,' and ',next_sec,' seconds'
            write(iulog,*) ' IOP end is          ',last_date,' and ',last_sec,' seconds'
         end if
         call endrun
      endif

      doiopupdate = .true.

!------------------------------------------------------------------------------
!     Check if iop data needs to be updated and set doiopupdate accordingly
!------------------------------------------------------------------------------
   else                      ! endstep > 1

      call timemgr_time_inc(bdate, 0, next_date, next_sec, inc_s=tsec(iopTimeIdx+1))

      call get_curr_date(yr, mon, day, ncsec)
      ncdate = yr*10000 + mon*100 + day

      if ( ncdate > next_date .or. (ncdate == next_date &
         .and. ncsec >= next_sec)) then
         iopTimeIdx = iopTimeIdx + 1
         doiopupdate = .true.
#if DEBUG > 2
         if (masterproc) write(iulog,*) sub//'nstep = ',get_nstep()
         if (masterproc) write(iulog,*) sub//'ncdate=',ncdate,' ncsec=',ncsec
         if (masterproc) write(iulog,*) sub//'next_date=',next_date,' next_sec=',next_sec
         if (masterproc) write(iulog,*) sub//':******* do iop update'
#endif 
      else
         doiopupdate = .false.
      end if
   endif                     ! if (endstep == 0 )
!
!     make sure we're
!     not going past end of iop data
!
   if ( ncdate > last_date .or. (ncdate == last_date &
      .and. ncsec > last_sec))  then
      if ( .not. scm_backfill_iop_w_init ) then
         call endrun(sub//':ERROR - setiopupdate.c:Reached the end of the time varient dataset')
      else
         doiopupdate = .false.              
      end if
   endif

#if DEBUG > 1
   if (masterproc) write(iulog,*) sub//':iop time index = ' , ioptimeidx
#endif

   return

end subroutine setiopupdate

end module iop

