module ic_us_standard_atmosphere

  !-----------------------------------------------------------------------
  !
  ! Purpose: Set US standard atmosphere initial conditions based on input coordinates
  !
  !
  !-----------------------------------------------------------------------
  use cam_logfile,         only: iulog
  use shr_kind_mod,        only: r8 => shr_kind_r8
  use cam_abortutils,      only: endrun
  use spmd_utils,          only: masterproc
  use shr_sys_mod,         only: shr_sys_flush

  implicit none
  private

  ! Public interface
  public :: us_std_atm_set_ic

  integer, parameter  :: nreg = 15  ! number of regions
  real(r8), parameter :: hb(nreg) = & ! height a bottom of layer (m)
       (/0.0_r8, 1.1e4_r8, 2.0e4_r8, 3.2e4_r8, 4.7e4_r8, 5.1e4_r8, 7.1e4_r8, 8.6e4_r8, &
       9.1e4_r8, 1.1e5_r8, 1.2e5_r8, 1.5e5_r8, 2.0e5_r8, 3.0e5_r8, 7.e5_r8/)
  real(r8), parameter :: pb(nreg) = & ! standard pressure (Pa)
       (/101325._r8, 22632.1_r8, 5474.89_r8, 868.02_r8, 110.91_r8, 66.94_r8, 3.96_r8, 3.7e-1_r8,  &
       1.5e-1_r8, 7.1e-3_r8, 2.5e-3_r8, 4.5e-4_r8, 8.47e-5_r8, 8.77e-6_r8, 3.19e-8_r8/)
  real(r8), parameter :: tb(nreg) = & ! standard temperature (K)
       (/288.15_r8, 216.65_r8, 216.65_r8, 228.65_r8, 270.65_r8, 270.65_r8, 214.65_r8, 186.87_r8,  &
       186.87_r8, 240._r8, 360._r8, 634.39_r8, 854.56_r8, 976.01_r8, 1.e3_r8/)
  real(r8), parameter :: lb(nreg) = & ! temperature lapse rate (K/m)
       (/-0.0065_r8, 0.0_r8, 0.001_r8, 0.0028_r8, 0.0_r8, -0.0028_r8, -0.001852_r8, 0.0_r8,       &
       2.796e-3_r8, 0.012_r8, 9.15e-3_r8, 4.4e-3_r8, 1.21e-3_r8, 6.e-5_r8, 0.0_r8/)
  real(r8), parameter :: rg = 8.3144598_r8 ! universal gas constant (J/mol/K)
  real(r8), parameter :: g0 = 9.80665_r8   ! gravitational acceleration (m/s^2)
  real(r8), parameter :: mw = 0.0289644_r8 ! molar mass of dry air (kg/mol)
  real(r8), parameter :: c1 = g0*mw/rg
  
!==============================================================================
CONTAINS
!==============================================================================

  subroutine us_std_atm_set_ic(latvals, lonvals, U, V, T, PS, PHIS,           &
       Q, m_cnst, mask, verbose)
    use const_init,    only: cnst_init_default
    use constituents,  only: cnst_name
    use hycoef,        only: ps0,hyam,hybm
    use physconst,     only: tref, rair,gravit
    
    !-----------------------------------------------------------------------
    !
    ! Purpose: Set Held-Suarez initial values for dynamics state variables
    !
    !-----------------------------------------------------------------------

    ! Dummy arguments
    real(r8),           intent(in)    :: latvals(:) ! lat in degrees (ncol)
    real(r8),           intent(in)    :: lonvals(:) ! lon in degrees (ncol)
    real(r8), optional, intent(inout) :: U(:,:)     ! zonal velocity
    real(r8), optional, intent(inout) :: V(:,:)     ! meridional velocity
    real(r8), optional, intent(inout) :: T(:,:)     ! temperature
    real(r8), optional, intent(inout) :: PS(:)      ! surface pressure
    real(r8), optional, intent(in)    :: PHIS(:)    ! surface geopotential
    real(r8), optional, intent(inout) :: Q(:,:,:)   ! tracer (ncol, lev, m)
    integer,  optional, intent(in)    :: m_cnst(:)  ! tracer indices (reqd. if Q)
    logical,  optional, intent(in)    :: mask(:)    ! Only init where .true.
    logical,  optional, intent(in)    :: verbose    ! For internal use

    ! Local variables
    logical, allocatable              :: mask_use(:)
    logical                           :: verbose_use
    integer                           :: i, k, m
    integer                           :: ncol
    integer                           :: nlev
    integer                           :: ncnst
    character(len=*), parameter       :: subname = 'us_std_atm_set_ic'
    real(r8)                          :: ptmp(1), ztmp(1)

    ncol = size(latvals, 1)
    allocate(mask_use(ncol))
    if (present(mask)) then
      if (size(mask_use) /= size(mask)) then
        call endrun('cnst_init_default: input, mask, is wrong size')
      end if
      mask_use = mask
    else
      mask_use = .true.
    end if

    if (present(verbose)) then
      verbose_use = verbose
    else
      verbose_use = .true.
    end if

    nlev = -1
    if (present(U)) then
      nlev = size(U, 2)
      do k = 1, nlev
        where(mask_use)
          U(:,k) = 0.0_r8
        end where
      end do
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          U initialized by "',subname,'"'
      end if
    end if

    if (present(V)) then
      nlev = size(V, 2)
      do k = 1, nlev
        where(mask_use)
          V(:,k) = 0.0_r8
        end where
      end do
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          V initialized by "',subname,'"'
      end if
    end if

    if (present(T)) then
      if (.not.present(PHIS)) then
        call endrun('PHIS must be specified to initiallize T in ic_us_standard_atm')
      end if      
      nlev = size(T, 2)
      do k = 1, nlev
        do i=1,ncol
          if (mask_use(i)) then
            ! get surface pressure
            call std_atm_pres(PHIS(i:i)/gravit, ptmp(:))
            ! get pressure level
            ptmp = hyam(k)*ps0+hybm(k)*ptmp
            ! get height of pressure level            
            call std_atm_height(ztmp(:), ptmp(:))
            ! given height get temperature
            call std_atm_temp(ztmp(:), T(i:i,k))
            if (i==1) write(*,*) "k,ztmp",k,ztmp,ptmp,T(i,k)
          end if
        end do
      end do
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          T initialized by "',subname,'"'
      end if
    end if

    if (present(PS)) then
      if (.not.present(PHIS)) then
        call endrun('PHIS must be specified to initiallize PS in ic_us_standard_atm')
      end if
      
      do i=1,ncol
        if (mask_use(i)) then
          call std_atm_pres(PHIS(i:i)/gravit, PS(i:i))
        end if
      end do
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          PS initialized by "',subname,'"'
      end if
    end if

    if (present(Q)) then
      nlev = size(Q, 2)
      ncnst = size(m_cnst, 1)
      do m = 1, ncnst
        if (m_cnst(m) == 1) then
          ! No water vapor in Held-Suarez
          do k = 1, nlev
            where(mask_use)
              Q(:,k,m_cnst(m)) = 0.0_r8
            end where
          end do
          if(masterproc .and. verbose_use) then
            write(iulog,*) '          ', trim(cnst_name(m_cnst(m))), ' initialized by "',subname,'"'
          end if
        else
          call cnst_init_default(m_cnst(m), latvals, lonvals, Q(:,:,m_cnst(m)),&
               mask=mask_use, verbose=verbose_use, notfound=.false.)
        end if
      end do
    end if

    deallocate(mask_use)

  end subroutine us_std_atm_set_ic

  subroutine std_atm_pres(height, pstd)
    
    ! Use barometric formula for U.S. Standard Atmosphere to convert heights to pressures.
    ! This formula is valid up to 86 km.
    ! https://en.wikipedia.org/wiki/Barometric_formula
    
    ! arguments
    real(r8), intent(in)  :: height(:) ! height above sea level in meters
    real(r8), intent(out) :: pstd(:)   ! std pressure in Pa
    
    integer :: i, ii, k, nlev
    logical :: found_region
    character(len=*), parameter :: routine = 'ref_pres::std_atm_pres'
    !---------------------------------------------------------------------------
    
    nlev = size(height)
    do k = 1, nlev
      if (height(k)<0.0_r8) then
        ii=1
        found_region = .true.
      else
        ! find region containing height
        found_region = .false.
        find_region: do i = nreg, 1, -1
          if (height(k) >= hb(i)) then
            ii = i
            found_region = .true.
            exit find_region
          end if
        end do find_region
      end if
      
      if (.not. found_region) then
        write(iulog,*) routine, ': illegal height: ', height(k)
        call endrun(routine// ': illegal height < 0. ')
      end if
      
      if (lb(ii) /= 0._r8) then
        pstd(k) = pb(ii) * ( tb(ii) / (tb(ii) + lb(ii)*(height(k) - hb(ii)) ) )**(c1/lb(ii))
      else
        pstd(k) = pb(ii) * exp( -c1*(height(k) - hb(ii))/tb(ii) )
      end if
      
    end do
  end subroutine std_atm_pres


  subroutine std_atm_height(height, pstd)
    
    ! Use barometric formula for U.S. Standard Atmosphere to convert heights to pressures.
    ! This formula is valid up to 86 km.
    ! https://en.wikipedia.org/wiki/Barometric_formula
    
    ! arguments
    real(r8), intent(out)  :: height(:) ! height above sea level in meters
    real(r8), intent(in)   :: pstd(:)   ! std pressure in Pa
    
    integer :: i, ii, k, nlev
    logical :: found_region
    character(len=*), parameter :: routine = 'ref_pres::std_atm_height'
    !---------------------------------------------------------------------------
    
    nlev = size(height)
    do k = 1, nlev
      
      ! find region containing height
      found_region = .false.
      find_region: do i = 1,nreg
        if (pstd(k) > pb(i)) then
          ii = i
          found_region = .true.
          exit find_region
        end if
      end do find_region
      
      if (.not. found_region) then
        write(iulog,*) routine, ': illegal pressure: ', pstd(k)
        call endrun(routine// ': illegal pressure. ')
      end if
      
      if (lb(ii) /= 0._r8) then
        height(k)  =         (tb(ii)*(pstd(k)/pb(ii))**(-lb(ii)/c1) - tb(ii))/lb(ii)+hb(ii)
      else
        height(k) = -(tb(ii)/c1)*log(pstd(k)/pb(ii))+ hb(ii)        
      end if
    end do
  end subroutine std_atm_height

  subroutine std_atm_temp(height, temp)
    
    ! Use barometric formula for U.S. Standard Atmosphere to convert heights to pressures.
    ! This formula is valid up to 86 km.
    ! https://en.wikipedia.org/wiki/Barometric_formula
    
    ! arguments
    real(r8), intent(out)  :: temp(:)   ! temperature
    real(r8), intent(in)   :: height(:)   ! std pressure in Pa
    
    ! local vars

    
    integer :: i, ii, k, nlev
    logical :: found_region
    character(len=*), parameter :: routine = 'ref_pres::std_atm_temp'
    !---------------------------------------------------------------------------
    
    nlev = size(height)
    do k = 1, nlev
      if (height(k)<0.0_r8) then
        ii=1
        found_region = .true.
      else
        ! find region containing height
        found_region = .false.
        find_region: do i = nreg, 1, -1
          if (height(k) >= hb(i)) then
            ii = i
            found_region = .true.
            exit find_region
          end if
        end do find_region
      end if

      if (.not. found_region) then
        write(iulog,*) routine, ': illegal height: ', height(k)
        call endrun(routine// ': illegal height < 0. ')
      end if
      
      if (lb(ii) /= 0._r8) then
        temp(k) = tb(ii) + lb(ii)*height(k)
      else
        temp(k) = tb(ii)
      end if

    end do
  end subroutine std_atm_temp
end module ic_us_standard_atmosphere
