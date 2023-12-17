module module_NoahMP_hrldas_driver

  USE module_hrldas_netcdf_io
  USE module_sf_noahmp_groundwater
  USE module_sf_noahmpdrv, only: noahmp_init, noahmplsm, noahmp_urban, groundwater_init, noahmplsm_mosaic_hue, noahmp_mosaic_init
  USE module_sf_urban, only: urban_param_init, urban_var_init
  USE module_date_utilities
  USE noahmp_tables, ONLY:  LCZ_1_TABLE,LCZ_2_TABLE,LCZ_3_TABLE,LCZ_4_TABLE,LCZ_5_TABLE, &
                 LCZ_6_TABLE,LCZ_7_TABLE,LCZ_8_TABLE,LCZ_9_TABLE,LCZ_10_TABLE,LCZ_11_TABLE
#ifdef MPP_LAND
  use module_mpp_land, only: MPP_LAND_PAR_INI, mpp_land_init, getLocalXY, mpp_land_bcast_char
  use module_mpp_land, only: check_land, my_id , node_info
  use module_cpl_land, only: cpl_land_init
#endif

  IMPLICIT NONE

#ifdef MPP_LAND
  include "mpif.h"
#endif

#ifdef WRF_HYDRO
   REAL,    allocatable, DIMENSION(:,:)   :: infxsrt,sfcheadrt,soldrain,qtiledrain,ZWATBLE2D
   integer :: forc_typ, snow_assim, HRLDAS_ini_typ
   REAL,    allocatable, DIMENSION(:,:)   :: etpnd, greenfrac, prcp0
   real :: etpnd1
   character(len=19) :: forcDate
   character(len = 256):: GEO_STATIC_FLNM
   real, allocatable, dimension(:) :: zsoil
   integer :: kk
   integer  :: finemesh, finemesh_factor
#endif

  character(len=9), parameter :: version = "v20150506"
  integer :: LDASIN_VERSION

!------------------------------------------------------------------------
! Begin exact copy of declaration section from driver (substitute allocatable, remove intent)
!------------------------------------------------------------------------

! IN only (as defined in WRF)

  INTEGER                                 ::  ITIMESTEP ! timestep number
  INTEGER                                 ::  YR        ! 4-digit year
  REAL                                    ::  JULIAN    ! Julian day
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  COSZEN    ! cosine zenith angle
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  XLAT      ! latitude [rad]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  DZ8W      ! thickness of atmo layers [m]
  REAL                                    ::  DTBL      ! timestep [s]
  REAL,    ALLOCATABLE, DIMENSION(:)      ::  DZS       ! thickness of soil layers [m]
  INTEGER                                 ::  NSOIL     ! number of soil layers
  INTEGER                                 ::  NUM_SOIL_LAYERS     ! number of soil layers
  REAL                                    ::  DX        ! horizontal grid spacing [m]
  INTEGER, ALLOCATABLE, DIMENSION(:,:)    ::  IVGTYP    ! vegetation type
  INTEGER, ALLOCATABLE, DIMENSION(:,:)    ::  ISLTYP    ! soil type
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  VEGFRA    ! vegetation fraction []
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  TMN       ! deep soil temperature [K]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  XLAND     ! =2 ocean; =1 land/seaice
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  XICE      ! fraction of grid that is seaice
  REAL                                    ::  XICE_THRESHOLD! fraction of grid determining seaice
  INTEGER                                 ::  ISICE     ! land cover category for ice
  INTEGER                                 ::  ISURBAN   ! land cover category for urban
  INTEGER                                 ::  ISWATER   ! land cover category for water
  INTEGER                                 ::  ISLAKE    ! land cover category for lake
  INTEGER                                 ::  IDVEG     ! dynamic vegetation (1 -> off ; 2 -> on) with opt_crs = 1
  INTEGER                                 ::  IOPT_CRS  ! canopy stomatal resistance (1-> Ball-Berry; 2->Jarvis)
  INTEGER                                 ::  IOPT_BTR  ! soil moisture factor for stomatal resistance (1-> Noah; 2-> CLM; 3-> SSiB)
  INTEGER                                 ::  IOPT_RUN  ! runoff and groundwater (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS)
  INTEGER                                 ::  IOPT_SFC  ! surface layer drag coeff (CH & CM) (1->M-O; 2->Chen97)
  INTEGER                                 ::  IOPT_FRZ  ! supercooled liquid water (1-> NY06; 2->Koren99)
  INTEGER                                 ::  IOPT_INF  ! frozen soil permeability (1-> NY06; 2->Koren99)
  INTEGER                                 ::  IOPT_RAD  ! radiation transfer (1->gap=F(3D,cosz); 2->gap=0; 3->gap=1-Fveg)
  INTEGER                                 ::  IOPT_ALB  ! snow surface albedo (1->BATS; 2->CLASS)
  INTEGER                                 ::  IOPT_SNF  ! rainfall & snowfall (1-Jordan91; 2->BATS; 3->Noah)
  INTEGER                                 ::  IOPT_TBOT ! lower boundary of soil temperature (1->zero-flux; 2->Noah)
  INTEGER                                 ::  IOPT_STC  ! snow/soil temperature time scheme
  INTEGER                                 ::  IOPT_GLA  ! glacier option (1->phase change; 2->simple)
  INTEGER                                 ::  IOPT_RSF  ! surface resistance option (1->Zeng; 2->simple)
  INTEGER                                 ::  IZ0TLND   ! option of Chen adjustment of Czil (not used)
  INTEGER                                 ::  IOPT_SOIL ! soil configuration option
  INTEGER                                 ::  IOPT_PEDO ! soil pedotransfer function option
  INTEGER                                 ::  IOPT_CROP ! crop model option (0->none; 1->Liu et al.; 2->Gecros)
  INTEGER                                 ::  IOPT_IRR  ! irrigation scheme (0->none; >1 irrigation scheme ON)
  INTEGER                                 ::  IOPT_IRRM ! irrigation method (0->dynamic; 1-> sprinkler; 2-> micro; 3-> flood)
  INTEGER                                 ::  IOPT_INFDV!infiltration options for dynamic VIC (1->Philip; 2-> Green-Ampt;3->Smith-Parlange)
  INTEGER                                 ::  IOPT_TDRN ! drainage option (0->off; 1->simple scheme; 2->Hooghoudt's scheme)
  INTEGER                                 ::  IOPT_MOSAIC
  INTEGER                                 ::  IOPT_HUE
  REAL                                    ::  soiltstep ! soil time step (s) control namelist option (default=0: same as main NoahMP timstep)
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  T_PHY     ! 3D atmospheric temperature valid at mid-levels [K]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  QV_CURR   ! 3D water vapor mixing ratio [kg/kg_dry]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  U_PHY     ! 3D U wind component [m/s]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  V_PHY     ! 3D V wind component [m/s]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SWDOWN    ! solar down at surface [W m-2]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  GLW       ! longwave down at surface [W m-2]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  P8W       ! 3D pressure, valid at interface [Pa]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RAINBL    ! precipitation entering land model [mm] per time step
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RAINBL_tmp! precipitation forcingentering land model [mm/s]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SNOWBL    ! snow entering land model [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SR        ! frozen precip ratio entering land model [-]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RAINCV    ! convective precip forcing [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RAINNCV   ! non-convective precip forcing [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RAINSHV   ! shallow conv. precip forcing [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SNOWNCV   ! non-covective snow forcing (subset of rainncv) [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  GRAUPELNCV! non-convective graupel forcing (subset of rainncv) [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  HAILNCV   ! non-convective hail forcing (subset of rainncv) [mm]

! Spatially varying fields

  REAL, ALLOCATABLE, DIMENSION(:,:,:)     ::  bexp_3D    ! C-H B exponent
  REAL, ALLOCATABLE, DIMENSION(:,:,:)     ::  smcdry_3D  ! Soil Moisture Limit: Dry
  REAL, ALLOCATABLE, DIMENSION(:,:,:)     ::  smcwlt_3D  ! Soil Moisture Limit: Wilt
  REAL, ALLOCATABLE, DIMENSION(:,:,:)     ::  smcref_3D  ! Soil Moisture Limit: Reference
  REAL, ALLOCATABLE, DIMENSION(:,:,:)     ::  smcmax_3D  ! Soil Moisture Limit: Max
  REAL, ALLOCATABLE, DIMENSION(:,:,:)     ::  dksat_3D   ! Saturated Soil Conductivity
  REAL, ALLOCATABLE, DIMENSION(:,:,:)     ::  dwsat_3D   ! Saturated Soil Diffusivity
  REAL, ALLOCATABLE, DIMENSION(:,:,:)     ::  psisat_3D  ! Saturated Matric Potential
  REAL, ALLOCATABLE, DIMENSION(:,:,:)     ::  quartz_3D  ! Soil quartz content
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  refdk_2D   ! Reference Soil Conductivity
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  refkdt_2D  ! Soil Infiltration Parameter
  REAL, ALLOCATABLE, DIMENSION(:,:,:)     ::  soilcomp   ! Soil sand and clay content [fraction]
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  soilcl1    ! Soil texture class with depth
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  soilcl2    ! Soil texture class with depth
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  soilcl3    ! Soil texture class with depth
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  soilcl4    ! Soil texture class with depth
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  irr_frac_2D! irrigation Fraction
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  irr_har_2D ! number of days before harvest date to stop irrigation
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  irr_lai_2D ! Minimum lai to trigger irrigation
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  irr_mad_2D ! management allowable deficit (0-1)
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  filoss_2D  ! fraction of flood irrigation loss (0-1)
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  sprir_rate_2D! mm/h, sprinkler irrigation rate
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  micir_rate_2D! mm/h, micro irrigation rate
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  firtfac_2D ! flood application rate factor
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  ir_rain_2D ! maximum precipitation to stop irrigation trigger
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  bvic_2d    ! VIC model infiltration parameter [-] opt_run=6
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  axaj_2D    ! Tension water distribution inflection parameter [-] opt_run=7
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  bxaj_2D    ! Tension water distribution shape parameter [-] opt_run=7
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  xxaj_2D    ! Free water distribution shape parameter [-] opt_run=7
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  bdvic_2d   ! VIC model infiltration parameter [-] opt_run=8
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  gdvic_2d   ! Mean Capillary Drive (m) for infiltration models opt_run=8
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  bbvic_2d   ! DVIC heterogeniety parameter for infiltration [-] opt_run=8
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  KLAT_FAC   ! factor multiplier to hydraulic conductivity
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  TDSMC_FAC  ! factor multiplier to field capacity
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  TD_DC      ! drainage coefficient for simple
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  TD_DCOEF   ! drainge coefficient for Hooghoudt
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  TD_DDRAIN  ! depth of drain
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  TD_RADI    ! tile radius
  REAL, ALLOCATABLE, DIMENSION(:,:)       ::  TD_SPAC    ! tile spacing

! INOUT (with generic LSM equivalent) (as defined in WRF)

  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  TSK       ! surface radiative temperature [K]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  HFX       ! sensible heat flux [W m-2]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QFX       ! latent heat flux [kg s-1 m-2]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  LH        ! latent heat flux [W m-2]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  GRDFLX    ! ground/snow heat flux [W m-2]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SMSTAV    ! soil moisture avail. [not used]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SMSTOT    ! total soil water [mm][not used]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SFCRUNOFF ! accumulated surface runoff [m]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  UDRUNOFF  ! accumulated sub-surface runoff [m]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ALBEDO    ! total grid albedo []
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SNOWC     ! snow cover fraction []
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  SMOISEQ   ! volumetric soil moisture [m3/m3]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  SMOIS     ! volumetric soil moisture [m3/m3]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  SH2O      ! volumetric liquid soil moisture [m3/m3]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  TSLB      ! soil temperature [K]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SNOW      ! snow water equivalent [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SNOWH     ! physical snow depth [m]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  CANWAT    ! total canopy water + ice [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ACSNOM    ! accumulated snow melt leaving pack
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ACSNOW    ! accumulated snow on grid
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  EMISS     ! surface bulk emissivity
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QSFC      ! bulk surface specific humidity

! INOUT (with no Noah LSM equivalent) (as defined in WRF)

  INTEGER, ALLOCATABLE, DIMENSION(:,:)    ::  ISNOWXY   ! actual no. of snow layers
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  TVXY      ! vegetation leaf temperature
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  TGXY      ! bulk ground surface temperature
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  CANICEXY  ! canopy-intercepted ice (mm)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  CANLIQXY  ! canopy-intercepted liquid water (mm)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  EAHXY     ! canopy air vapor pressure (pa)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  TAHXY     ! canopy air temperature (k)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  CMXY      ! bulk momentum drag coefficient
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  CHXY      ! bulk sensible heat exchange coefficient
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  FWETXY    ! wetted or snowed fraction of the canopy (-)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SNEQVOXY  ! snow mass at last time step(mm h2o)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ALBOLDXY  ! snow albedo at last time step (-)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QSNOWXY   ! snowfall on the ground [mm/s]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QRAINXY   ! rainfall on the ground [mm/s]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  WSLAKEXY  ! lake water storage [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ZWTXY     ! water table depth [m]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  WAXY      ! water in the "aquifer" [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  WTXY      ! groundwater storage [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SMCWTDXY  ! groundwater storage [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  DEEPRECHXY! groundwater storage [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RECHXY    ! groundwater storage [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  TSNOXY    ! snow temperature [K]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  ZSNSOXY   ! snow layer depth [m]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  SNICEXY   ! snow layer ice [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  SNLIQXY   ! snow layer liquid water [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  LFMASSXY  ! leaf mass [g/m2]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RTMASSXY  ! mass of fine roots [g/m2]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  STMASSXY  ! stem mass [g/m2]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  WOODXY    ! mass of wood (incl. woody roots) [g/m2]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  GRAINXY   ! XING mass of grain!THREE
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  GDDXY     ! XINGgrowingdegressday
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  STBLCPXY  ! stable carbon in deep soil [g/m2]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  FASTCPXY  ! short-lived carbon, shallow soil [g/m2]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  LAI       ! leaf area index
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  LAI_tmp       ! leaf area index
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  XSAIXY    ! stem area index
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  TAUSSXY   ! snow age factor

  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  gecros_state   ! gecros crop model packed state vector

!irrigation
  REAL,    ALLOCATABLE, DIMENSION(:,:)    :: IRFRACT    ! irrigation fraction
  REAL,    ALLOCATABLE, DIMENSION(:,:)    :: SIFRACT    ! sprinkler irrigation fraction
  REAL,    ALLOCATABLE, DIMENSION(:,:)    :: MIFRACT    ! micro irrigation fraction
  REAL,    ALLOCATABLE, DIMENSION(:,:)    :: FIFRACT    ! flood irrigation fraction
  INTEGER, ALLOCATABLE, DIMENSION(:,:)    :: IRNUMSI    ! irrigation event number, Sprinkler
  INTEGER, ALLOCATABLE, DIMENSION(:,:)    :: IRNUMMI    ! irrigation event number, Micro
  INTEGER, ALLOCATABLE, DIMENSION(:,:)    :: IRNUMFI    ! irrigation event number, Flood
  REAL,    ALLOCATABLE, DIMENSION(:,:)    :: IRWATSI    ! irrigation water amount [m] to be applied, Sprinkler
  REAL,    ALLOCATABLE, DIMENSION(:,:)    :: IRWATMI    ! irrigation water amount [m] to be applied, Micro
  REAL,    ALLOCATABLE, DIMENSION(:,:)    :: IRWATFI    ! irrigation water amount [m] to be applied, Flood
  REAL,    ALLOCATABLE, DIMENSION(:,:)    :: IRELOSS    ! loss of irrigation water to evaporation,sprinkler [m/timestep]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    :: IRSIVOL    ! amount of irrigation by sprinkler (mm)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    :: IRMIVOL    ! amount of irrigation by micro (mm)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    :: IRFIVOL    ! amount of irrigation by micro (mm)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    :: IRRSPLH    ! latent heating from sprinkler evaporation (w/m2)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    :: LOCTIM     ! local time

! OUT (with no Noah LSM equivalent) (as defined in WRF)

  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  T2MVXY    ! 2m temperature of vegetation part
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  T2MBXY    ! 2m temperature of bare ground part
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  Q2MVXY    ! 2m mixing ratio of vegetation part
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  Q2MBXY    ! 2m mixing ratio of bare ground part
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  TRADXY    ! surface radiative temperature (k)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  NEEXY     ! net ecosys exchange (g/m2/s CO2)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  GPPXY     ! gross primary assimilation [g/m2/s C]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  NPPXY     ! net primary productivity [g/m2/s C]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  FVEGXY    ! Noah-MP vegetation fraction [-]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RUNSFXY   ! surface runoff [mm/s]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RUNSBXY   ! subsurface runoff [mm/s]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ECANXY    ! evaporation of intercepted water (mm/s)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  EDIRXY    ! soil surface evaporation rate (mm/s]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ETRANXY   ! transpiration rate (mm/s)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  FSAXY     ! total absorbed solar radiation (w/m2)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  FIRAXY    ! total net longwave rad (w/m2) [+ to atm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  APARXY    ! photosyn active energy by canopy (w/m2)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  PSNXY     ! total photosynthesis (umol co2/m2/s) [+]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SAVXY     ! solar rad absorbed by veg. (w/m2)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SAGXY     ! solar rad absorbed by ground (w/m2)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RSSUNXY   ! sunlit leaf stomatal resistance (s/m)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RSSHAXY   ! shaded leaf stomatal resistance (s/m)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  BGAPXY    ! between gap fraction
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  WGAPXY    ! within gap fraction
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  TGVXY     ! under canopy ground temperature[K]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  TGBXY     ! bare ground temperature [K]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  CHVXY     ! sensible heat exchange coefficient vegetated
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  CHBXY     ! sensible heat exchange coefficient bare-ground
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SHGXY     ! veg ground sen. heat [w/m2]   [+ to atm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SHCXY     ! canopy sen. heat [w/m2]   [+ to atm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SHBXY     ! bare sensible heat [w/m2]  [+ to atm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  EVGXY     ! veg ground evap. heat [w/m2]  [+ to atm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  EVBXY     ! bare soil evaporation [w/m2]  [+ to atm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  GHVXY     ! veg ground heat flux [w/m2]  [+ to soil]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  GHBXY     ! bare ground heat flux [w/m2] [+ to soil]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  IRGXY     ! veg ground net LW rad. [w/m2] [+ to atm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  IRCXY     ! canopy net LW rad. [w/m2] [+ to atm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  IRBXY     ! bare net longwave rad. [w/m2] [+ to atm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  TRXY      ! transpiration [w/m2]  [+ to atm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  EVCXY     ! canopy evaporation heat [w/m2]  [+ to atm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  CHLEAFXY  ! leaf exchange coefficient
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  CHUCXY    ! under canopy exchange coefficient
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  CHV2XY    ! veg 2m exchange coefficient
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  CHB2XY    ! bare 2m exchange coefficient
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RS        ! Total stomatal resistance (s/m)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  Z0        ! roughness length output to WRF
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ZNT       ! roughness length output to WRF
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QTDRAIN   ! tile drain discharge (mm)
! additional output variables
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  PAHXY     ! precipitation advected heat
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  PAHGXY    ! precipitation advected heat
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  PAHBXY    ! precipitation advected heat
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  PAHVXY    ! precipitation advected heat
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QINTSXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QINTRXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QDRIPSXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QDRIPRXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QTHROSXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QTHRORXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QSNSUBXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QMELTXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QSNFROXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QSUBCXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QFROCXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QEVACXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QDEWCXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QFRZCXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QMELTCXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QSNBOTXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  PONDINGXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  FPICEXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RAINLSM
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SNOWLSM
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  FORCTLSM
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  FORCQLSM
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  FORCPLSM
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  FORCZLSM
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  FORCWLSM
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ACC_SSOILXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ACC_QINSURXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ACC_QSEVAXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  EFLXBXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SOILENERGY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SNOWENERGY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  CANHSXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ACC_DWATERXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ACC_PRCPXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ACC_ECANXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ACC_ETRANXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ACC_EDIRXY
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  ACC_ETRANIXY

  INTEGER   ::  ids,ide, jds,jde, kds,kde,  &  ! d -> domain
   &            ims,ime, jms,jme, kms,kme,  &  ! m -> memory
   &            its,ite, jts,jte, kts,kte      ! t -> tile

!------------------------------------------------------------------------
! Needed for NoahMP init
!------------------------------------------------------------------------

  LOGICAL                                 ::  FNDSOILW    ! soil water present in input
  LOGICAL                                 ::  FNDSNOWH    ! snow depth present in input
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  CHSTARXY    ! for consistency with MP_init; delete later
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SEAICE      ! seaice fraction

!------------------------------------------------------------------------
! Needed for MMF_RUNOFF (IOPT_RUN = 5); not part of MP driver in WRF
!------------------------------------------------------------------------

  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  MSFTX
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  MSFTY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  EQZWT
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RIVERBEDXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RIVERCONDXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  PEXPXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  FDEPTHXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  AREAXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QRFSXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QSPRINGSXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QRFXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QSPRINGXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QSLATXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QLATXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RECHCLIM
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RIVERMASK
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  NONRIVERXY
  REAL                                    ::  WTDDT  = 30.0    ! frequency of groundwater call [minutes]
  INTEGER                                 ::  STEPWTD          ! step of groundwater call

!------------------------------------------------------------------------
! Needed for TILE DRAINAGE IF IOPT_TDRN = 1 OR 2
!------------------------------------------------------------------------
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  TD_FRACTION

!------------------------------------------------------------------------
! Needed for crop model (OPT_CROP=1)
!------------------------------------------------------------------------

  INTEGER, ALLOCATABLE, DIMENSION(:,:)   :: PGSXY
  INTEGER, ALLOCATABLE, DIMENSION(:,:)   :: CROPCAT
  REAL   , ALLOCATABLE, DIMENSION(:,:)   :: PLANTING
  REAL   , ALLOCATABLE, DIMENSION(:,:)   :: HARVEST
  REAL   , ALLOCATABLE, DIMENSION(:,:)   :: SEASON_GDD
  REAL   , ALLOCATABLE, DIMENSION(:,:,:) :: CROPTYPE

!------------------------------------------------------------------------
! Single- and Multi-layer Urban Models
!------------------------------------------------------------------------

  INTEGER                                 ::  num_urban_atmosphere ! atmospheric levels including ZLVL for BEP/BEM models

  REAL,    ALLOCATABLE                    ::  GMT       ! Hour of day (fractional) (needed for urban)
  INTEGER, ALLOCATABLE                    ::  JULDAY    ! Integer day (needed for urban)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  HRANG     ! hour angle (needed for urban)
  REAL,    ALLOCATABLE                    ::  DECLIN    ! declination (needed for urban)
  INTEGER                                 ::  num_roof_layers = 4
  INTEGER                                 ::  num_road_layers = 4
  INTEGER                                 ::  num_wall_layers = 4
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  cmr_sfcdif
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  chr_sfcdif
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  cmc_sfcdif
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  chc_sfcdif
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  cmgr_sfcdif
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  chgr_sfcdif
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  tr_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  tb_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  tg_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  tc_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  qc_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  uc_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  xxxr_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  xxxb_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  xxxg_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  xxxc_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  trl_urb3d
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  tbl_urb3d
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  tgl_urb3d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  sh_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  lh_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  g_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  rn_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ts_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  psim_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  psih_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  u10_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  v10_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  GZ1OZ0_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  AKMS_URB2D
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  th2_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  q2_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ust_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:)      ::  dzr
  REAL,    ALLOCATABLE, DIMENSION(:)      ::  dzb
  REAL,    ALLOCATABLE, DIMENSION(:)      ::  dzg
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  cmcr_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  tgr_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  tgrl_urb3d
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  smr_urb3d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  drelr_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  drelb_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  drelg_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  flxhumr_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  flxhumb_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  flxhumg_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  frc_urb2d
  INTEGER, ALLOCATABLE, DIMENSION(:,:)    ::  utype_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  chs
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  chs2
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  cqs2
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  trb_urb4d
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  tw1_urb4d
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  tw2_urb4d
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  tgb_urb4d
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  tlev_urb3d
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  qlev_urb3d
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  tw1lev_urb3d
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  tw2lev_urb3d
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  tglev_urb3d
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  tflev_urb3d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  sf_ac_urb3d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  lf_ac_urb3d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  cm_ac_urb3d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  sfvent_urb3d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  lfvent_urb3d
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  sfwin1_urb3d
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  sfwin2_urb3d
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  sfw1_urb3d
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  sfw2_urb3d
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  sfr_urb3d
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  sfg_urb3d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  lp_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  hi_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  lb_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  hgt_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  mh_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  stdh_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  lf_urb2d
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  theta_urban
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  u_urban
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  v_urban
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  dz_urban
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  rho_urban
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  p_urban
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ust
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  a_u_bep
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  a_v_bep
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  a_t_bep
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  a_q_bep
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  a_e_bep
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  b_u_bep
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  b_v_bep
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  b_t_bep
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  b_q_bep
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  b_e_bep
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  dlg_bep
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  dl_u_bep
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  sf_bep
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  vl_bep
  REAL                                    ::  height_urban
! new urban variables for green roof, PVP for BEP_BEM scheme=3, Zonato et al., 2021
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  EP_PV_URB3D
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QGR_URB3D
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  TGR_URB3D
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  DRAINGR_URB3D
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  T_PV_URB3D
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  TRV_URB4D
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  QR_URB4D
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  DRAIN_URB4D
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  SFRV_URB3D
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  LFRV_URB3D
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  DGR_URB3D
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  DG_URB3D
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  LFR_URB3D
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  LFG_URB3D
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SWDDIR    ! solar down at surface [W m-2]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SWDDIF

  !------------------------------------------------------------------------
  ! NOAHMP Mosaic Variables added -Aaron A.
  !------------------------------------------------------------------------
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: LANDUSEF
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: LANDUSEF2
  INTEGER, ALLOCATABLE,DIMENSION(:,:,:)  :: mosaic_cat_index

  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: TSK_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: HFX_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: QFX_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: LH_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: TMN_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: GRDFLX_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: SFCRUNOFF_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: UDRUNOFF_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: ALBEDO_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: SNOWC_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: CANWAT_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: SNOW_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: SNOWH_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: ACSNOM_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: ACSNOW_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: EMISS_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: QSFC_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: Z0_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: ZNT_mosaic

  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  ::  VEGFRA_mosaic    ! vegetation fraction []

  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: tvxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: tgxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: canicexy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: canliqxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: eahxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: tahxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: cmxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: chxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: fwetxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: sneqvoxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: alboldxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: qsnowxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: qrainxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: wslakexy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: zwtxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: waxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: wtxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: lfmassxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: rtmassxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: stmassxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: woodxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: grainxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: gddxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: PGSXY_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: smcwtdxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: stblcpxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: fastcpxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: xsaixy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: xlai_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: taussxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: smcwtdxy_mosiac
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: deeprechxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: rechxy_mosaic

  ! Irrigation Variables
  !irrigation
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    :: IRFRACT_mosaic    ! irrigation fraction
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    :: SIFRACT_mosaic    ! sprinkler irrigation fraction
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    :: MIFRACT_mosaic    ! micro irrigation fraction
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    :: FIFRACT_mosaic    ! flood irrigation fraction
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:)    :: IRNUMSI_mosaic    ! irrigation event number, Sprinkler
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:)    :: IRNUMMI_mosaic    ! irrigation event number, Micro
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:)    :: IRNUMFI_mosaic    ! irrigation event number, Flood
 
 REAL,    ALLOCATABLE, DIMENSION(:,:,:)    :: IRWATSI_mosaic    ! irrigation water amount [m] to be applied, Sprinkler
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    :: IRWATMI_mosaic    ! irrigation water amount [m] to be applied, Micro
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    :: IRWATFI_mosaic    ! irrigation water amount [m] to be applied, Flood
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    :: IRELOSS_mosaic    ! loss of irrigation water to evaporation,sprinkler [m/timestep]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    :: IRSIVOL_mosaic    ! amount of irrigation by sprinkler (mm)
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    :: IRMIVOL_mosaic    ! amount of irrigation by micro (mm)
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    :: IRFIVOL_mosaic    ! amount of irrigation by micro (mm)
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    :: IRRSPLH_mosaic    ! latent heating from sprinkler evaporation (w/m2)
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    :: LOCTIM_mosaic     ! local time

  !Out variables that need to be allocated for mosaic.
  ! OUT (with no Noah LSM equivalent) (as defined in WRF)

  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: t2mvxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: t2mbxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: chstarxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: q2mvxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: q2mbxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: tradxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: neexy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: gppxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: nppxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: fvegxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: runsfxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: runsbxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: ecanxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: edirxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: etranxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: fsaxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: firaxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: aparxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: psnxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: savxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: sagxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: rssunxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: rsshaxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: bgapxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: wgapxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: tgvxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: tgbxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: chvxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: chbxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: shgxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: shcxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: shbxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: evgxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: evbxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: ghvxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: ghbxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: irgxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: ircxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: irbxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: trxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: evcxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: chleafxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: chucxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: chv2xy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: chb2xy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: rs_mosaic

  ! Addtional outputs
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  PAHXY_mosaic     ! precipitation advected heat
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  PAHGXY_mosaic    ! precipitation advected heat
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  PAHBXY_mosaic    ! precipitation advected heat
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  PAHVXY_mosaic    ! precipitation advected heat
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  QINTSXY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  QINTRXY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  QDRIPSXY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  QDRIPRXY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  QTHROSXY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  QTHRORXY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  QSNSUBXY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  QMELTXY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  QSNFROXY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  QSUBCXY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  QFROCXY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  QEVACXY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  QDEWCXY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  QFRZCXY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  QMELTCXY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  QSNBOTXY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  PONDINGXY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  FPICEXY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  ACC_SSOILXY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  ACC_QINSURXY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  ACC_QSEVAXY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  EFLXBXY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  SOILENERGY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  SNOWENERGY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  CANHSXY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  ACC_DWATERXY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  ACC_PRCPXY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  ACC_ECANXY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  ACC_ETRANXY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  ACC_EDIRXY_mosaic
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  ACC_ETRANIXY_mosaic

  ! SNOW VARIABLES
  INTEGER,   ALLOCATABLE, DIMENSION(:,:,:)  :: isnowxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: zsnsoxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: tsnoxy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: snicexy_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: snliqxy_mosaic
  ! VARIABLES THAT ARE MULTI-DIMENSIONAL (these involve a careful initilization)
  REAL,   ALLOCATABLE, DIMENSION(:,:,:) :: TSLB_mosaic
  REAl,   ALLOCATABLE, DIMENSION(:,:,:) :: SMOIS_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:) :: SH2O_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:) :: SMOISEQ_mosaic

  !urban model paramters
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: TR_URB2D_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: TB_URB2D_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: TG_URB2D_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: TC_URB2D_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: QC_URB2D_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: UC_URB2D_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: SH_URB2D_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: LH_URB2D_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: G_URB2D_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: RN_URB2D_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: TS_URB2D_mosaic

  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: CMR_SFCDIF_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: CHR_SFCDIF_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: CMC_SFCDIF_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: CHC_SFCDIF_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: CMGR_SFCDIF_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: CHGR_SFCDIF_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: XXXR_URB2D_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: XXXB_URB2D_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: XXXG_URB2D_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: XXXC_URB2D_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: CMCR_URB2D_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: TGR_URB2D_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: DRELR_URB2D_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: DRELB_URB2D_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: DRELG_URB2D_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: FLXHUMR_URB2D_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: FLXHUMB_URB2D_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: FLXHUMG_URB2D_mosaic

  ! HUE Variables
  REAL,   ALLOCATABLE, DIMENSION(:,:)  :: RUNONSFXY
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: RUNONSFXY_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:)  :: DETENTION_STORAGEXY
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: DETENTION_STORAGEXY_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: VOL_FLUX_RUNONXY_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)  :: VOL_FLUX_SMXY_mosaic


  !urban 3d model

  REAL,   ALLOCATABLE, DIMENSION(:,:,:) :: TRL_URB3D_mosaic
  REAl,   ALLOCATABLE, DIMENSION(:,:,:) :: TBL_URB3D_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:) :: TGL_URB3D_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:) :: TGRL_URB3D_mosaic
  REAL,   ALLOCATABLE, DIMENSION(:,:,:) :: SMR_URB3D_mosaic
!------------------------------------------------------------------------
! END OF MOSAIC VARIABLES CALLS
!------------------------------------------------------------------------
!------------------------------------------------------------------------
! 2D variables not used in WRF - should be removed?
!------------------------------------------------------------------------

  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  XLONG       ! longitude
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  TERRAIN     ! terrain height
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  GVFMIN      ! annual minimum in vegetation fraction
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  GVFMAX      ! annual maximum in vegetation fraction

!------------------------------------------------------------------------
! End 2D variables not used in WRF
!------------------------------------------------------------------------

  CHARACTER(LEN=256) :: MMINSL  = 'STAS'  ! soil classification
  CHARACTER(LEN=256) :: LLANDUSE          ! (=USGS, using USGS landuse classification)

!------------------------------------------------------------------------
! Timing:
!------------------------------------------------------------------------

  INTEGER :: NTIME          ! timesteps
  integer :: clock_count_1 = 0
  integer :: clock_count_2 = 0
  integer :: clock_rate    = 0
  real    :: timing_sum    = 0.0

  integer :: sflx_count_sum
  integer :: count_before_sflx
  integer :: count_after_sflx

!---------------------------------------------------------------------
!  DECLARE/Initialize constants
!---------------------------------------------------------------------

    INTEGER                             :: I
    INTEGER                             :: J
    INTEGER                             :: SLOPETYP
    INTEGER                             :: YEARLEN
    INTEGER, PARAMETER                  :: NSNOW = 3    ! number of snow layers fixed to 3
    REAL, PARAMETER                     :: undefined_real = 9.9692099683868690E36 ! NetCDF float   FillValue
    INTEGER, PARAMETER                  :: undefined_int = -2147483647            ! NetCDF integer FillValue
    LOGICAL                             :: update_lai, update_veg
    INTEGER                             :: spinup_loop
    LOGICAL                             :: reset_spinup_date

!---------------------------------------------------------------------
!  File naming, parallel
!---------------------------------------------------------------------

  character(len=19)  :: olddate, newdate, startdate
  character          :: hgrid
  integer            :: igrid
  logical            :: lexist
  integer            :: imode
  integer            :: ixfull
  integer            :: jxfull
  integer            :: ixpar
  integer            :: jxpar
  integer            :: xstartpar
  integer            :: ystartpar
  integer            :: rank = 0
  CHARACTER(len=256) :: inflnm, outflnm, inflnm_template
  logical            :: restart_flag
  character(len=256) :: restart_flnm
  integer            :: ierr

!---------------------------------------------------------------------
! Attributes from LDASIN input file (or HRLDAS_SETUP_FILE, as the case may be)
!---------------------------------------------------------------------

  INTEGER           :: IX
  INTEGER           :: JX
  REAL              :: DY
  REAL              :: TRUELAT1
  REAL              :: TRUELAT2
  REAL              :: CEN_LON
  INTEGER           :: MAPPROJ
  REAL              :: LAT1
  REAL              :: LON1

#ifdef MPP_LAND
  integer ix_tmp, jx_tmp
#endif

!---------------------------------------------------------------------
!  NAMELIST start
!---------------------------------------------------------------------

  character(len=256) :: indir
  ! nsoil defined above
  integer            :: forcing_timestep
  integer            :: noah_timestep
  integer            :: start_year
  integer            :: start_month
  integer            :: start_day
  integer            :: start_hour
  integer            :: start_min
  character(len=256) :: outdir = "."
  character(len=256) :: restart_filename_requested = " "
  integer            :: restart_frequency_hours
  integer            :: output_timestep
  integer            :: spinup_loops

  !!Namelist modifications for NOAH-MP HUE added by Aaron Alexander
  integer            :: noahmp_mosaic_scheme
  character(len=256) :: geogrid_file_name_for_mosaic = " "
  integer            :: number_mosaic_catagories
  integer            :: number_land_use_catagories
  integer            :: noahmp_HUE_iopt
  !! End modifications for NOAH-MP HUE

  integer            :: sf_urban_physics = 0
  integer            :: use_wudapt_lcz   = 0 ! add for LCZ urban
  integer            :: num_urban_ndm    = 1
  integer            :: num_urban_ng     = 1
  integer            :: num_urban_nwr    = 1
  integer            :: num_urban_ngb    = 1
  integer            :: num_urban_nf     = 1
  integer            :: num_urban_nz     = 1
  integer            :: num_urban_nbui   = 1
  integer            :: num_urban_hi     = 15
  real               :: urban_atmosphere_thickness = 2.0
! new urban var for green roof and solar panel
  integer            :: num_urban_ngr    = 10 ! = ngr_u in bep_bem.F
  integer            :: urban_map_zgrd   = 1

  ! derived urban dimensions

  integer            :: urban_map_zrd
  integer            :: urban_map_zwd
  integer            :: urban_map_gd
  integer            :: urban_map_zd
  integer            :: urban_map_zdf
  integer            :: urban_map_bd
  integer            :: urban_map_wd
  integer            :: urban_map_gbd
  integer            :: urban_map_fbd

  character(len=256) :: forcing_name_T = "T2D"
  character(len=256) :: forcing_name_Q = "Q2D"
  character(len=256) :: forcing_name_U = "U2D"
  character(len=256) :: forcing_name_V = "V2D"
  character(len=256) :: forcing_name_P = "PSFC"
  character(len=256) :: forcing_name_LW = "LWDOWN"
  character(len=256) :: forcing_name_SW = "SWDOWN"
  character(len=256) :: forcing_name_PR = "RAINRATE"
  character(len=256) :: forcing_name_SN = ""

  integer            :: dynamic_veg_option
  integer            :: canopy_stomatal_resistance_option
  integer            :: btr_option
  integer            :: runoff_option
  integer            :: surface_drag_option
  integer            :: supercooled_water_option
  integer            :: frozen_soil_option
  integer            :: radiative_transfer_option
  integer            :: snow_albedo_option
  integer            :: pcp_partition_option
  integer            :: tbot_option
  integer            :: temp_time_scheme_option
  integer            :: glacier_option
  integer            :: surface_resistance_option

  integer            :: soil_data_option = 1
  integer            :: pedotransfer_option = 1
  integer            :: crop_option = 0
  integer            :: irrigation_option = 0
  integer            :: irrigation_method = 0
  integer            :: dvic_infiltration_option = 1
  integer            :: tile_drainage_option = 0
  real               :: soil_timestep_option = 0.0
  integer            :: noahmp_output = 0

  integer            :: split_output_count = 1
  logical            :: skip_first_output = .false.
  integer            :: khour
  integer            :: kday
  real               :: zlvl
  character(len=256) :: hrldas_setup_file = " "
  character(len=256) :: spatial_filename = " "
  character(len=256) :: external_veg_filename_template = " "
  character(len=256) :: external_lai_filename_template = " "
  character(len=256) :: agdata_flnm = " "
  character(len=256) :: tdinput_flnm = " "
  integer            :: xstart = 1
  integer            :: ystart = 1
  integer            ::   xend = 0
  integer            ::   yend = 0
  integer, PARAMETER    :: MAX_SOIL_LEVELS = 10   ! maximum soil levels in namelist
  REAL, DIMENSION(MAX_SOIL_LEVELS) :: soil_thick_input       ! depth to soil interfaces from namelist [m]

  namelist / NOAHLSM_OFFLINE /    &
#ifdef WRF_HYDRO
       finemesh,finemesh_factor,forc_typ, snow_assim , GEO_STATIC_FLNM, HRLDAS_ini_typ, &
#endif
       indir, nsoil, soil_thick_input, forcing_timestep, noah_timestep, &
       start_year, start_month, start_day, start_hour, start_min, &
       outdir, skip_first_output, &
       restart_filename_requested, restart_frequency_hours, output_timestep, &
       spinup_loops, &
       forcing_name_T,forcing_name_Q,forcing_name_U,forcing_name_V,forcing_name_P, &
       forcing_name_LW,forcing_name_SW,forcing_name_PR,forcing_name_SN, &

       dynamic_veg_option, canopy_stomatal_resistance_option, &
       btr_option, runoff_option, surface_drag_option, supercooled_water_option, &
       frozen_soil_option, radiative_transfer_option, snow_albedo_option, &
       pcp_partition_option, tbot_option, temp_time_scheme_option, &
       glacier_option, surface_resistance_option, &
       irrigation_option, irrigation_method, dvic_infiltration_option, &
       tile_drainage_option,soil_timestep_option,noahmp_output,&
       soil_data_option, pedotransfer_option, crop_option, &
       sf_urban_physics,use_wudapt_lcz,num_urban_hi,urban_atmosphere_thickness, &
       num_urban_ndm,num_urban_ng,num_urban_nwr ,num_urban_ngb , &
       num_urban_nf ,num_urban_nz,num_urban_nbui, &
       split_output_count, &
       khour, kday, zlvl, hrldas_setup_file, &
       spatial_filename, agdata_flnm, tdinput_flnm, &
       external_veg_filename_template, external_lai_filename_template, &
       xstart, xend, ystart, yend, &

       !!Added by Aaron Alexander
       noahmp_mosaic_scheme, geogrid_file_name_for_mosaic,number_mosaic_catagories,number_land_use_catagories, &
       noahmp_HUE_iopt
       !!End Mods by Aaron Alexander


  contains

  subroutine land_driver_ini(NTIME_out,wrfits,wrfite,wrfjts,wrfjte)

     USE module_bep_bem_helper, ONLY: nurbm ! C.He 03/31/2021: add for bep_bem treatment in WRFv4.3
     USE module_domain, ONLY : domain ! for groundwater_init WRF

     implicit  none
     integer:: NTIME_out

     REAL :: max_utype_urb2d
     TYPE (domain)          :: grid

    ! initilization for stand alone parallel code.
    integer, optional, intent(in) :: wrfits,wrfite,wrfjts,wrfjte
#ifdef MPP_LAND
    call  MPP_LAND_INIT()
#endif

! Initialize namelist variables to dummy values, so we can tell
! if they have not been set properly.

  nsoil                   = -999
  soil_thick_input        = -999
  dtbl                    = -999
  start_year              = -999
  start_month             = -999
  start_day               = -999
  start_hour              = -999
  start_min               = -999
  khour                   = -999
  kday                    = -999
  zlvl                    = -999
  forcing_timestep        = -999
  noah_timestep           = -999
  output_timestep         = -999
  spinup_loops            = 0
  restart_frequency_hours = -999

  !!NOAH-MP HUE added by Aaron Alexander
    noahmp_mosaic_scheme     = -999 !should true or false
    number_mosaic_catagories = -999 !should be greater than zero but less than 8
    number_land_use_catagories = -999 !should be set so we can read in the fractional land use data
    noahmp_HUE_iopt = -999
  !!NOAH-MP HUE end mods by Aaron Alexander

  open(30, file="namelist.hrldas", form="FORMATTED")
  read(30, NOAHLSM_OFFLINE, iostat=ierr)
  if (ierr /= 0) then
     write(*,'(/," ***** ERROR: Problem reading namelist NOAHLSM_OFFLINE",/)')
     rewind(30)
     read(30, NOAHLSM_OFFLINE)
     stop " ***** ERROR: Problem reading namelist NOAHLSM_OFFLINE"
  endif
  close(30)

  dtbl = real(noah_timestep)
  num_soil_layers = nsoil      ! because surface driver uses the long form
  IDVEG = dynamic_veg_option ! transfer from namelist to driver format
  IOPT_CRS = canopy_stomatal_resistance_option
  IOPT_BTR = btr_option
  IOPT_RUN = runoff_option
  IOPT_SFC = surface_drag_option
  IOPT_FRZ = supercooled_water_option
  IOPT_INF = frozen_soil_option
  IOPT_RAD = radiative_transfer_option
  IOPT_ALB = snow_albedo_option
  IOPT_SNF = pcp_partition_option
  IOPT_TBOT = tbot_option
  IOPT_STC = temp_time_scheme_option
  IOPT_GLA = glacier_option
  IOPT_RSF = surface_resistance_option
  IOPT_SOIL = soil_data_option
  IOPT_PEDO = pedotransfer_option
  IOPT_CROP = crop_option
  IOPT_IRR  = irrigation_option
  IOPT_IRRM = irrigation_method
  IOPT_INFDV= dvic_infiltration_option
  IOPT_TDRN = tile_drainage_option
  soiltstep = soil_timestep_option
  IOPT_MOSAIC = noahmp_mosaic_scheme
  IOPT_HUE = noahmp_HUE_iopt

!---------------------------------------------------------------------
!  NAMELIST end
!---------------------------------------------------------------------

!C.He: add for urban update in WRFv4.3
!-----------------------------------------------------------------------
! Urban physics set up. If the run-time option for use_wudapt_lcz = 0,
! then the number of urban classes is 3. Else, if the use_wudapt_lcz = 1,
! then the number increases to 11. The seemingly local variable
! assignment, "nurbm", is actually USE associated from the BEP BEM
! helper module.
!-----------------------------------------------------------------------
   IF ( use_wudapt_lcz .EQ. 0 ) THEN
      nurbm = 3
   ELSE IF ( use_wudapt_lcz .EQ. 1 ) THEN
      nurbm = 11
   END IF

!---------------------------------------------------------------------
!  NAMELIST check begin
!---------------------------------------------------------------------

  update_lai = .true.   ! default: use LAI if present in forcing file
  if (dynamic_veg_option == 1 .or. dynamic_veg_option == 2 .or. &
      dynamic_veg_option == 3 .or. dynamic_veg_option == 4 .or. &
      dynamic_veg_option == 5 .or. dynamic_veg_option == 6) &    ! remove dveg=10 and add dveg=1,3,4 into the update_lai flag false condition
    update_lai = .false.

  update_veg = .false.  ! default: don't use VEGFRA if present in forcing file
  if (dynamic_veg_option == 1 .or. dynamic_veg_option == 6 .or. dynamic_veg_option == 7) &
    update_veg = .true.

  if (nsoil < 0) then
     stop " ***** ERROR: NSOIL must be set in the namelist."
  endif

  if ((khour < 0) .and. (kday < 0)) then
     write(*, '(" ***** Namelist error: ************************************")')
     write(*, '(" ***** ")')
     write(*, '(" *****      Either KHOUR or KDAY must be defined.")')
     write(*, '(" ***** ")')
     stop
  else if (( khour < 0 ) .and. (kday > 0)) then
     khour = kday * 24
  else if ((khour > 0) .and. (kday > 0)) then
     write(*, '("Namelist warning:  KHOUR and KDAY both defined.")')
  else
     ! all is well.  KHOUR defined
  endif

  if (forcing_timestep < 0) then
        write(*, *)
        write(*, '(" ***** Namelist error: *****************************************")')
        write(*, '(" ***** ")')
        write(*, '(" *****       FORCING_TIMESTEP needs to be set greater than zero.")')
        write(*, '(" ***** ")')
        write(*, *)
        stop
  endif

  if (noah_timestep < 0) then
        write(*, *)
        write(*, '(" ***** Namelist error: *****************************************")')
        write(*, '(" ***** ")')
        write(*, '(" *****       NOAH_TIMESTEP needs to be set greater than zero.")')
        write(*, '(" *****                     900 seconds is recommended.       ")')
        write(*, '(" ***** ")')
        write(*, *)
        stop
  endif

  !
  ! Check that OUTPUT_TIMESTEP fits into NOAH_TIMESTEP:
  !
  if (output_timestep /= 0) then
     if (mod(output_timestep, noah_timestep) > 0) then
        write(*, *)
        write(*, '(" ***** Namelist error: *********************************************************")')
        write(*, '(" ***** ")')
        write(*, '(" *****       OUTPUT_TIMESTEP should set to an integer multiple of NOAH_TIMESTEP.")')
        write(*, '(" *****            OUTPUT_TIMESTEP = ", I12, " seconds")') output_timestep
        write(*, '(" *****            NOAH_TIMESTEP   = ", I12, " seconds")') noah_timestep
        write(*, '(" ***** ")')
        write(*, *)
        stop
     endif
  endif

  !
  ! Check that RESTART_FREQUENCY_HOURS fits into NOAH_TIMESTEP:
  !
  if (restart_frequency_hours /= 0) then
     if (mod(restart_frequency_hours*3600, noah_timestep) > 0) then
        write(*, *)
        write(*, '(" ***** Namelist error: ******************************************************")')
        write(*, '(" ***** ")')
        write(*, '(" *****       RESTART_FREQUENCY_HOURS (converted to seconds) should set to an ")')
        write(*, '(" *****       integer multiple of NOAH_TIMESTEP.")')
        write(*, '(" *****            RESTART_FREQUENCY_HOURS = ", I12, " hours:  ", I12, " seconds")') &
             restart_frequency_hours, restart_frequency_hours*3600
        write(*, '(" *****            NOAH_TIMESTEP           = ", I12, " seconds")') noah_timestep
        write(*, '(" ***** ")')
        write(*, *)
        stop
     endif
  endif

  if (dynamic_veg_option == 2 .or. dynamic_veg_option == 5 .or. dynamic_veg_option == 6) then
     if ( canopy_stomatal_resistance_option /= 1) then
        write(*, *)
        write(*, '(" ***** Namelist error: ******************************************************")')
        write(*, '(" ***** ")')
        write(*, '(" *****       CANOPY_STOMATAL_RESISTANCE_OPTION must be 1 when DYNAMIC_VEG_OPTION == 2/5/6")')
        write(*, *)
        stop
     endif
  endif

  if (soil_data_option == 4 .and. spatial_filename == " ") then
        write(*, *)
        write(*, '(" ***** Namelist error: ******************************************************")')
        write(*, '(" ***** ")')
        write(*, '(" *****       SPATIAL_FILENAME must be provided when SOIL_DATA_OPTION == 4")')
        write(*, *)
        stop
  endif

  if (sf_urban_physics == 2 .or. sf_urban_physics == 3) then
     if ( urban_atmosphere_thickness <= 0.0) then
        write(*, *)
        write(*, '(" ***** Namelist error: ******************************************************")')
        write(*, '(" ***** ")')
        write(*, '(" *****       When running BEP/BEM, URBAN_ATMOSPHERE_LEVELS must contain at least 3 levels")')
        write(*, *)
        stop
     endif
     num_urban_atmosphere = int(zlvl/urban_atmosphere_thickness)
     if (zlvl - num_urban_atmosphere*urban_atmosphere_thickness >= 0.5*urban_atmosphere_thickness)  &
            num_urban_atmosphere = num_urban_atmosphere + 1
     if ( num_urban_atmosphere <= 2) then
        write(*, *)
        write(*, '(" ***** Namelist error: ******************************************************")')
        write(*, '(" ***** ")')
        write(*, '(" *****       When running BEP/BEM, num_urban_atmosphere must contain at least 3 levels, ")')
        write(*, '(" *****        decrease URBAN_ATMOSPHERE_THICKNESS")')
        write(*, *)
        stop
     endif
  endif

  ! Checks for NOAH-MP HUE. By Aaron Alexander 19 May 2022
  if (IOPT_MOSAIC /= 0 .and. IOPT_MOSAIC /= 1) then
        write(*, *)
        write(*, '(" ***** Namelist error: ******************************************************")')
        write(*, '(" ***** ")')
        write(*, '(" *****       Mosaic option must be either 0 (off) or 1 (on)")')
        write(*, *)
        stop
  endif

  if (IOPT_MOSAIC == 1 ) then
     if (number_mosaic_catagories <2 .or. number_mosaic_catagories >8) then
        write(*, *)
        write(*, '(" ***** Namelist error: ******************************************************")')
        write(*, '(" ***** ")')
        write(*, '(" *****       Noahmp mosaic catagories must be larger than 2, but less than 8")')
        write(*, *)
        stop
     endif
  endif

  if (IOPT_MOSAIC == 1 ) then
     if (number_land_use_catagories == -999) then
        write(*, *)
        write(*, '(" ***** Namelist error: ******************************************************")')
        write(*, '(" ***** ")')
        write(*, '(" *****       HUE NOAH-MP requires you to set NUMBER_LAND_USE_CATAGORIES namlist variable")')
        write(*, *)
        stop
     endif
  endif

  if (IOPT_MOSAIC == 1 ) then
     if (sf_urban_physics >=2 ) then
        write(*, *)
        write(*, '(" ***** Namelist error: ******************************************************")')
        write(*, '(" ***** ")')
        write(*, '(" *****      HUE NOAH-MP is not currently set up to run with BEP or BEP-BEM")')
        write(*, *)
        stop
     endif
  endif

  if (IOPT_MOSAIC == 1 ) then
     if (tile_drainage_option > 0 ) then
        write(*, *)
        write(*, '(" ***** Namelist error: ******************************************************")')
        write(*, '(" ***** ")')
        write(*, '(" *****      HUE NOAH-MP is not currently set up to run with Tile Drainage Packages")')
        write(*, *)
        stop
     endif
  endif

  if (IOPT_HUE /= 0 .and. IOPT_HUE /= 1) then
        write(*, *)
        write(*, '(" ***** Namelist error: ******************************************************")')
        write(*, '(" ***** ")')
        write(*, '(" *****       HUE NOAH-MP option must be either 0 (off) or 1 (on)")')
        write(*, *)
        stop
  endif

  if (IOPT_HUE.eq.1) then
    if (number_land_use_catagories.le.40 ) then
       write(*, *)
       write(*, '(" ***** Namelist error: ******************************************************")')
       write(*, '(" ***** ")')
       write(*, '(" *****      HUE NOAH-MP must have more than 40 Land-types to be used. If not, no new physics will be implemented ")')
       write(*, *)
       stop
    endif
 endif

 if (IOPT_MOSAIC.eq.1) then
   if (iopt_crop.gt.0 ) then
      write(*, *)
      write(*, '(" ***** Namelist error: ******************************************************")')
      write(*, '(" ***** ")')
      write(*, '(" *****      Mosaic NOAH-MP does not support Geocros or Liu et al. 2016 crop model at this time ")')
      write(*, *)
      stop
   endif
endif

if (IOPT_MOSAIC.eq.1) then
  if (iopt_run.eq.5 ) then
     write(*, *)
     write(*, '(" ***** Namelist error: ******************************************************")')
     write(*, '(" ***** ")')
     write(*, '(" *****      Mosaic NOAH-MP does not support shallow groundwater model at this time  ")')
     write(*, *)
     stop
  endif
endif

if (IOPT_HUE.eq.1) then
  if (IOPT_MOSAIC.eq.0) then
     write(*, *)
     write(*, '(" ***** Namelist error: ******************************************************")')
     write(*, '(" ***** ")')
     write(*, '(" *****      Mosaic NOAH-MP must be used with HUE in spatially distributed format ")')
     write(*, *)
     stop
  endif
endif

!---------------------------------------------------------------------
!  NAMELIST check end
!---------------------------------------------------------------------

  ! derived urban dimensions

      if (sf_urban_physics > 0 ) then

         urban_map_zrd = num_urban_ndm * num_urban_nwr * num_urban_nz
         urban_map_zwd = num_urban_ndm * num_urban_nwr * num_urban_nz  * num_urban_nbui
         urban_map_gd  = num_urban_ndm * num_urban_ng
         urban_map_zd  = num_urban_ndm * num_urban_nz  * num_urban_nbui
         urban_map_zdf = num_urban_ndm * num_urban_nz
         urban_map_bd  = num_urban_nz  * num_urban_nbui
         urban_map_wd  = num_urban_ndm * num_urban_nz  * num_urban_nbui
         urban_map_gbd = num_urban_ndm * num_urban_ngb * num_urban_nbui
         urban_map_fbd = num_urban_ndm * (num_urban_nz - 1)  * &
                         num_urban_nf  * num_urban_nbui

! new urban var
        urban_map_zgrd = num_urban_ndm * num_urban_ngr * num_urban_nz

      end if

!----------------------------------------------------------------------
! Initialize gridded domain
!----------------------------------------------------------------------

#ifdef MPP_LAND
#ifdef WRF_HYDRO
    if(finemesh .ne. 0) then
         xstart = (wrfits-1)*finemesh_factor + 1
         xend = (wrfite-1)*finemesh_factor
         ystart = (wrfjts-1)*finemesh_factor + 1
         yend = (wrfjte-1)*finemesh_factor
         call CPL_LAND_INIT(xstart,xend, ystart,yend)
         ix_tmp = xend - xstart + 1
         jx_tmp = yend - ystart + 1
    else
#endif
       call read_dim(hrldas_setup_file,ix_tmp,jx_tmp)
       call MPP_LAND_PAR_INI(1,ix_tmp,jx_tmp,1)
       call getLocalXY(ix_tmp,jx_tmp,xstart,ystart,xend,yend)
#ifdef WRF_HYDRO
    endif
#endif
#endif

  call read_hrldas_hdrinfo(hrldas_setup_file, ix, jx, xstart, xend, ystart, yend, &
       iswater, islake, isurban, isice, llanduse, dx, dy, truelat1, truelat2, cen_lon, lat1, lon1, &
       igrid, mapproj)
  write(hgrid,'(I1)') igrid

  write(olddate,'(I4.4,"-",I2.2,"-",I2.2,"_",I2.2,":",I2.2,":",I2.2)') &
       start_year, start_month, start_day, start_hour, start_min, 0

  startdate = olddate

#ifdef MPP_LAND
   ix = ix_tmp
   jx = jx_tmp
#endif

#ifdef WRF_HYDRO
    forcdate = olddate
#endif

  ids = xstart
  ide = xend
  jds = ystart
  jde = yend
  kds = 1
  kde = 2
  its = xstart
  ite = xend
  jts = ystart
  jte = yend
  kts = 1
  kte = 2
  ims = xstart
  ime = xend
  jms = ystart
  jme = yend
  kms = 1
  kme = 2

  if (sf_urban_physics == 2 .or. sf_urban_physics == 3) then
    kde = num_urban_atmosphere
    kte = num_urban_atmosphere
    kme = num_urban_atmosphere
  endif

!---------------------------------------------------------------------
!  Allocate multi-dimension fields for subwindow calculation
!---------------------------------------------------------------------

  ixfull = xend-xstart+1
  jxfull = yend-ystart+1

  ixpar = ixfull
  jxpar = jxfull
  xstartpar = 1
  ystartpar = 1

  ALLOCATE ( COSZEN    (XSTART:XEND,YSTART:YEND) )    ! cosine zenith angle
  ALLOCATE ( XLAT      (XSTART:XEND,YSTART:YEND) )    ! latitude [rad]
  ALLOCATE ( DZ8W      (XSTART:XEND,KDS:KDE,YSTART:YEND) )  ! thickness of atmo layers [m]
  ALLOCATE ( DZS       (1:NSOIL)                   )  ! thickness of soil layers [m]
  ALLOCATE ( IVGTYP    (XSTART:XEND,YSTART:YEND) )    ! vegetation type
  ALLOCATE ( ISLTYP    (XSTART:XEND,YSTART:YEND) )    ! soil type
  ALLOCATE ( VEGFRA    (XSTART:XEND,YSTART:YEND) )    ! vegetation fraction []
  ALLOCATE ( TMN       (XSTART:XEND,YSTART:YEND) )    ! deep soil temperature [K]
  ALLOCATE ( XLAND     (XSTART:XEND,YSTART:YEND) )    ! =2 ocean; =1 land/seaice
  ALLOCATE ( XICE      (XSTART:XEND,YSTART:YEND) )    ! fraction of grid that is seaice
  ALLOCATE ( T_PHY     (XSTART:XEND,KDS:KDE,YSTART:YEND) )  ! 3D atmospheric temperature valid at mid-levels [K]
  ALLOCATE ( QV_CURR   (XSTART:XEND,KDS:KDE,YSTART:YEND) )  ! 3D water vapor mixing ratio [kg/kg_dry]
  ALLOCATE ( U_PHY     (XSTART:XEND,KDS:KDE,YSTART:YEND) )  ! 3D U wind component [m/s]
  ALLOCATE ( V_PHY     (XSTART:XEND,KDS:KDE,YSTART:YEND) )  ! 3D V wind component [m/s]
  ALLOCATE ( SWDOWN    (XSTART:XEND,YSTART:YEND) )    ! solar down at surface [W m-2]
  ALLOCATE ( SWDDIR    (XSTART:XEND,YSTART:YEND) )    ! solar down at surface [W m-2] for new urban solar panel
  ALLOCATE ( SWDDIF    (XSTART:XEND,YSTART:YEND) )    ! solar down at surface [W m-2] for new urban solar panel
  ALLOCATE ( GLW       (XSTART:XEND,YSTART:YEND) )    ! longwave down at surface [W m-2]
  ALLOCATE ( P8W       (XSTART:XEND,KDS:KDE,YSTART:YEND) )  ! 3D pressure, valid at interface [Pa]
  ALLOCATE ( RAINBL    (XSTART:XEND,YSTART:YEND) )    ! total precipitation entering land model [mm] per time step
  ALLOCATE ( SNOWBL    (XSTART:XEND,YSTART:YEND) )    ! snow entering land model [mm] per time step
  ALLOCATE ( RAINBL_tmp    (XSTART:XEND,YSTART:YEND) )    ! precipitation entering land model [mm]
  ALLOCATE ( SR        (XSTART:XEND,YSTART:YEND) )    ! frozen precip ratio entering land model [-]
  ALLOCATE ( RAINCV    (XSTART:XEND,YSTART:YEND) )    ! convective precip forcing [mm]
  ALLOCATE ( RAINNCV   (XSTART:XEND,YSTART:YEND) )    ! non-convective precip forcing [mm]
  ALLOCATE ( RAINSHV   (XSTART:XEND,YSTART:YEND) )    ! shallow conv. precip forcing [mm]
  ALLOCATE ( SNOWNCV   (XSTART:XEND,YSTART:YEND) )    ! non-covective snow forcing (subset of rainncv) [mm]
  ALLOCATE ( GRAUPELNCV(XSTART:XEND,YSTART:YEND) )    ! non-convective graupel forcing (subset of rainncv) [mm]
  ALLOCATE ( HAILNCV   (XSTART:XEND,YSTART:YEND) )    ! non-convective hail forcing (subset of rainncv) [mm]

  ALLOCATE ( bexp_3d    (XSTART:XEND,1:NSOIL,YSTART:YEND) )    ! C-H B exponent
  ALLOCATE ( smcdry_3D  (XSTART:XEND,1:NSOIL,YSTART:YEND) )    ! Soil Moisture Limit: Dry
  ALLOCATE ( smcwlt_3D  (XSTART:XEND,1:NSOIL,YSTART:YEND) )    ! Soil Moisture Limit: Wilt
  ALLOCATE ( smcref_3D  (XSTART:XEND,1:NSOIL,YSTART:YEND) )    ! Soil Moisture Limit: Reference
  ALLOCATE ( smcmax_3D  (XSTART:XEND,1:NSOIL,YSTART:YEND) )    ! Soil Moisture Limit: Max
  ALLOCATE ( dksat_3D   (XSTART:XEND,1:NSOIL,YSTART:YEND) )    ! Saturated Soil Conductivity
  ALLOCATE ( dwsat_3D   (XSTART:XEND,1:NSOIL,YSTART:YEND) )    ! Saturated Soil Diffusivity
  ALLOCATE ( psisat_3D  (XSTART:XEND,1:NSOIL,YSTART:YEND) )    ! Saturated Matric Potential
  ALLOCATE ( quartz_3D  (XSTART:XEND,1:NSOIL,YSTART:YEND) )    ! Soil quartz content
  ALLOCATE ( refdk_2D   (XSTART:XEND,YSTART:YEND) )            ! Reference Soil Conductivity
  ALLOCATE ( refkdt_2D  (XSTART:XEND,YSTART:YEND) )            ! Soil Infiltration Parameter
  ALLOCATE ( soilcomp   (XSTART:XEND,1:2*NSOIL,YSTART:YEND) )  ! Soil sand and clay content [fraction]
  ALLOCATE ( soilcl1    (XSTART:XEND,YSTART:YEND) )            ! Soil texture class with depth
  ALLOCATE ( soilcl2    (XSTART:XEND,YSTART:YEND) )            ! Soil texture class with depth
  ALLOCATE ( soilcl3    (XSTART:XEND,YSTART:YEND) )            ! Soil texture class with depth
  ALLOCATE ( soilcl4    (XSTART:XEND,YSTART:YEND) )            ! Soil texture class with depth
  ALLOCATE ( irr_frac_2D(XSTART:XEND,YSTART:YEND) )            ! irrigation Fraction
  ALLOCATE ( irr_har_2D (XSTART:XEND,YSTART:YEND) )            ! number of days before harvest date to stop irrigation
  ALLOCATE ( irr_lai_2D (XSTART:XEND,YSTART:YEND) )            ! Minimum lai to trigger irrigation
  ALLOCATE ( irr_mad_2D (XSTART:XEND,YSTART:YEND) )            ! management allowable deficit (0-1)
  ALLOCATE ( filoss_2D  (XSTART:XEND,YSTART:YEND) )            ! fraction of flood irrigation loss (0-1)
  ALLOCATE ( sprir_rate_2D(XSTART:XEND,YSTART:YEND) )          ! mm/h, sprinkler irrigation rate
  ALLOCATE ( micir_rate_2D(XSTART:XEND,YSTART:YEND) )          ! mm/h, micro irrigation rate
  ALLOCATE ( firtfac_2D (XSTART:XEND,YSTART:YEND) )            ! flood application rate factor
  ALLOCATE ( ir_rain_2D (XSTART:XEND,YSTART:YEND) )            ! maximum precipitation to stop irrigation trigger
  ALLOCATE ( bvic_2D    (XSTART:XEND,YSTART:YEND) )            ! VIC model infiltration parameter [-]
  ALLOCATE ( axaj_2D    (XSTART:XEND,YSTART:YEND) )            ! Tension water distribution inflection parameter [-]
  ALLOCATE ( bxaj_2D    (XSTART:XEND,YSTART:YEND) )            ! Tension water distribution shape parameter [-]
  ALLOCATE ( xxaj_2D    (XSTART:XEND,YSTART:YEND) )            ! Free water distribution shape parameter [-]
  ALLOCATE ( bdvic_2D   (XSTART:XEND,YSTART:YEND) )            ! DVIC model infiltration parameter [-]
  ALLOCATE ( gdvic_2D   (XSTART:XEND,YSTART:YEND) )            ! Mean Capillary Drive (m) for infiltration models
  ALLOCATE ( bbvic_2D   (XSTART:XEND,YSTART:YEND) )            ! DVIC heterogeniety parameter for infiltration [-]
  ALLOCATE ( KLAT_FAC   (XSTART:XEND,YSTART:YEND) )            ! factor multiplier to hydraulic conductivity
  ALLOCATE ( TDSMC_FAC  (XSTART:XEND,YSTART:YEND) )            ! factor multiplier to field capacity
  ALLOCATE ( TD_DC      (XSTART:XEND,YSTART:YEND) )            ! drainage coefficient for simple
  ALLOCATE ( TD_DCOEF   (XSTART:XEND,YSTART:YEND) )            ! drainge coefficient for Hooghoudt
  ALLOCATE ( TD_DDRAIN  (XSTART:XEND,YSTART:YEND) )            ! depth of drain
  ALLOCATE ( TD_RADI    (XSTART:XEND,YSTART:YEND) )            ! tile radius
  ALLOCATE ( TD_SPAC    (XSTART:XEND,YSTART:YEND) )            ! tile spacing

! INOUT (with generic LSM equivalent) (as defined in WRF)

  ALLOCATE ( TSK       (XSTART:XEND,YSTART:YEND) )  ! surface radiative temperature [K]
  ALLOCATE ( HFX       (XSTART:XEND,YSTART:YEND) )  ! sensible heat flux [W m-2]
  ALLOCATE ( QFX       (XSTART:XEND,YSTART:YEND) )  ! latent heat flux [kg s-1 m-2]
  ALLOCATE ( LH        (XSTART:XEND,YSTART:YEND) )  ! latent heat flux [W m-2]
  ALLOCATE ( GRDFLX    (XSTART:XEND,YSTART:YEND) )  ! ground/snow heat flux [W m-2]
  ALLOCATE ( SMSTAV    (XSTART:XEND,YSTART:YEND) )  ! soil moisture avail. [not used]
  ALLOCATE ( SMSTOT    (XSTART:XEND,YSTART:YEND) )  ! total soil water [mm][not used]
  ALLOCATE ( SFCRUNOFF (XSTART:XEND,YSTART:YEND) )  ! accumulated surface runoff [m]
  ALLOCATE ( UDRUNOFF  (XSTART:XEND,YSTART:YEND) )  ! accumulated sub-surface runoff [m]
  ALLOCATE ( ALBEDO    (XSTART:XEND,YSTART:YEND) )  ! total grid albedo []
  ALLOCATE ( SNOWC     (XSTART:XEND,YSTART:YEND) )  ! snow cover fraction []
  ALLOCATE ( SMOISEQ   (XSTART:XEND,1:NSOIL,YSTART:YEND) )     ! eq volumetric soil moisture [m3/m3]
  ALLOCATE ( SMOIS     (XSTART:XEND,1:NSOIL,YSTART:YEND) )     ! volumetric soil moisture [m3/m3]
  ALLOCATE ( SH2O      (XSTART:XEND,1:NSOIL,YSTART:YEND) )     ! volumetric liquid soil moisture [m3/m3]
  ALLOCATE ( TSLB      (XSTART:XEND,1:NSOIL,YSTART:YEND) )     ! soil temperature [K]
  ALLOCATE ( SNOW      (XSTART:XEND,YSTART:YEND) )  ! snow water equivalent [mm]
  ALLOCATE ( SNOWH     (XSTART:XEND,YSTART:YEND) )  ! physical snow depth [m]
  ALLOCATE ( CANWAT    (XSTART:XEND,YSTART:YEND) )  ! total canopy water + ice [mm]
  ALLOCATE ( ACSNOM    (XSTART:XEND,YSTART:YEND) )  ! accumulated snow melt leaving pack
  ALLOCATE ( ACSNOW    (XSTART:XEND,YSTART:YEND) )  ! accumulated snow on grid
  ALLOCATE ( EMISS     (XSTART:XEND,YSTART:YEND) )  ! surface bulk emissivity
  ALLOCATE ( QSFC      (XSTART:XEND,YSTART:YEND) )  ! bulk surface specific humidity

! INOUT (with no Noah LSM equivalent) (as defined in WRF)

  ALLOCATE ( ISNOWXY   (XSTART:XEND,YSTART:YEND) )  ! actual no. of snow layers
  ALLOCATE ( TVXY      (XSTART:XEND,YSTART:YEND) )  ! vegetation leaf temperature
  ALLOCATE ( TGXY      (XSTART:XEND,YSTART:YEND) )  ! bulk ground surface temperature
  ALLOCATE ( CANICEXY  (XSTART:XEND,YSTART:YEND) )  ! canopy-intercepted ice (mm)
  ALLOCATE ( CANLIQXY  (XSTART:XEND,YSTART:YEND) )  ! canopy-intercepted liquid water (mm)
  ALLOCATE ( EAHXY     (XSTART:XEND,YSTART:YEND) )  ! canopy air vapor pressure (pa)
  ALLOCATE ( TAHXY     (XSTART:XEND,YSTART:YEND) )  ! canopy air temperature (k)
  ALLOCATE ( CMXY      (XSTART:XEND,YSTART:YEND) )  ! bulk momentum drag coefficient
  ALLOCATE ( CHXY      (XSTART:XEND,YSTART:YEND) )  ! bulk sensible heat exchange coefficient
  ALLOCATE ( FWETXY    (XSTART:XEND,YSTART:YEND) )  ! wetted or snowed fraction of the canopy (-)
  ALLOCATE ( SNEQVOXY  (XSTART:XEND,YSTART:YEND) )  ! snow mass at last time step(mm h2o)
  ALLOCATE ( ALBOLDXY  (XSTART:XEND,YSTART:YEND) )  ! snow albedo at last time step (-)
  ALLOCATE ( QSNOWXY   (XSTART:XEND,YSTART:YEND) )  ! snowfall on the ground [mm/s]
  ALLOCATE ( QRAINXY   (XSTART:XEND,YSTART:YEND) )  ! rainfall on the ground [mm/s]
  ALLOCATE ( WSLAKEXY  (XSTART:XEND,YSTART:YEND) )  ! lake water storage [mm]
  ALLOCATE ( ZWTXY     (XSTART:XEND,YSTART:YEND) )  ! water table depth [m]
  ALLOCATE ( WAXY      (XSTART:XEND,YSTART:YEND) )  ! water in the "aquifer" [mm]
  ALLOCATE ( WTXY      (XSTART:XEND,YSTART:YEND) )  ! groundwater storage [mm]
  ALLOCATE ( SMCWTDXY  (XSTART:XEND,YSTART:YEND) )  ! soil moisture below the bottom of the column (m3m-3)
  ALLOCATE ( DEEPRECHXY(XSTART:XEND,YSTART:YEND) )  ! recharge to the water table when deep (m)
  ALLOCATE ( RECHXY    (XSTART:XEND,YSTART:YEND) )  ! recharge to the water table (diagnostic) (m)
  ALLOCATE ( TSNOXY    (XSTART:XEND,-NSNOW+1:0,    YSTART:YEND) )  ! snow temperature [K]
  ALLOCATE ( ZSNSOXY   (XSTART:XEND,-NSNOW+1:NSOIL,YSTART:YEND) )  ! snow layer depth [m]
  ALLOCATE ( SNICEXY   (XSTART:XEND,-NSNOW+1:0,    YSTART:YEND) )  ! snow layer ice [mm]
  ALLOCATE ( SNLIQXY   (XSTART:XEND,-NSNOW+1:0,    YSTART:YEND) )  ! snow layer liquid water [mm]
  ALLOCATE ( LFMASSXY  (XSTART:XEND,YSTART:YEND) )  ! leaf mass [g/m2]
  ALLOCATE ( RTMASSXY  (XSTART:XEND,YSTART:YEND) )  ! mass of fine roots [g/m2]
  ALLOCATE ( STMASSXY  (XSTART:XEND,YSTART:YEND) )  ! stem mass [g/m2]
  ALLOCATE ( WOODXY    (XSTART:XEND,YSTART:YEND) )  ! mass of wood (incl. woody roots) [g/m2]
  ALLOCATE ( GRAINXY   (XSTART:XEND,YSTART:YEND) )  ! mass of grain XING [g/m2]
  ALLOCATE ( GDDXY     (XSTART:XEND,YSTART:YEND) )  ! growing degree days XING FOUR
  ALLOCATE ( STBLCPXY  (XSTART:XEND,YSTART:YEND) )  ! stable carbon in deep soil [g/m2]
  ALLOCATE ( FASTCPXY  (XSTART:XEND,YSTART:YEND) )  ! short-lived carbon, shallow soil [g/m2]
  ALLOCATE ( LAI       (XSTART:XEND,YSTART:YEND) )  ! leaf area index
  ALLOCATE ( LAI_tmp   (XSTART:XEND,YSTART:YEND) )  ! leaf area index
  ALLOCATE ( XSAIXY    (XSTART:XEND,YSTART:YEND) )  ! stem area index
  ALLOCATE ( TAUSSXY   (XSTART:XEND,YSTART:YEND) )  ! snow age factor

  ALLOCATE ( gecros_state(XSTART:XEND,60,YSTART:YEND) )   ! gecros crop model packed state vector
! irrigation
  ALLOCATE ( IRFRACT   (XSTART:XEND,YSTART:YEND) )  ! irrigation fraction
  ALLOCATE ( SIFRACT   (XSTART:XEND,YSTART:YEND) )  ! sprinkler irrigation fraction
  ALLOCATE ( MIFRACT   (XSTART:XEND,YSTART:YEND) )  ! micro irrigation fraction
  ALLOCATE ( FIFRACT   (XSTART:XEND,YSTART:YEND) )  ! flood irrigation fraction
  ALLOCATE ( IRNUMSI   (XSTART:XEND,YSTART:YEND) )  ! irrigation event number, Sprinkler
  ALLOCATE ( IRNUMMI   (XSTART:XEND,YSTART:YEND) )  ! irrigation event number, Micro
  ALLOCATE ( IRNUMFI   (XSTART:XEND,YSTART:YEND) )  ! irrigation event number, Flood
  ALLOCATE ( IRWATSI   (XSTART:XEND,YSTART:YEND) )  ! irrigation water amount [m] to be applied, Sprinkler
  ALLOCATE ( IRWATMI   (XSTART:XEND,YSTART:YEND) )  ! irrigation water amount [m] to be applied, Micro
  ALLOCATE ( IRWATFI   (XSTART:XEND,YSTART:YEND) )  ! irrigation water amount [m] to be applied, Flood
  ALLOCATE ( IRELOSS   (XSTART:XEND,YSTART:YEND) )  ! loss of irrigation water to evaporation,sprinkler [mm]
  ALLOCATE ( IRSIVOL   (XSTART:XEND,YSTART:YEND) )  ! amount of irrigation by sprinkler (mm)
  ALLOCATE ( IRMIVOL   (XSTART:XEND,YSTART:YEND) )  ! amount of irrigation by micro (mm)
  ALLOCATE ( IRFIVOL   (XSTART:XEND,YSTART:YEND) )  ! amount of irrigation by micro (mm)
  ALLOCATE ( IRRSPLH   (XSTART:XEND,YSTART:YEND) )  ! latent heating from sprinkler evaporation (w/m2)
  ALLOCATE ( LOCTIM    (XSTART:XEND,YSTART:YEND) )  ! local time

! OUT (with no Noah LSM equivalent) (as defined in WRF)

  ALLOCATE ( T2MVXY    (XSTART:XEND,YSTART:YEND) )  ! 2m temperature of vegetation part
  ALLOCATE ( T2MBXY    (XSTART:XEND,YSTART:YEND) )  ! 2m temperature of bare ground part
  ALLOCATE ( Q2MVXY    (XSTART:XEND,YSTART:YEND) )  ! 2m mixing ratio of vegetation part
  ALLOCATE ( Q2MBXY    (XSTART:XEND,YSTART:YEND) )  ! 2m mixing ratio of bare ground part
  ALLOCATE ( TRADXY    (XSTART:XEND,YSTART:YEND) )  ! surface radiative temperature (k)
  ALLOCATE ( NEEXY     (XSTART:XEND,YSTART:YEND) )  ! net ecosys exchange (g/m2/s CO2)
  ALLOCATE ( GPPXY     (XSTART:XEND,YSTART:YEND) )  ! gross primary assimilation [g/m2/s C]
  ALLOCATE ( NPPXY     (XSTART:XEND,YSTART:YEND) )  ! net primary productivity [g/m2/s C]
  ALLOCATE ( FVEGXY    (XSTART:XEND,YSTART:YEND) )  ! Noah-MP vegetation fraction [-]
  ALLOCATE ( RUNSFXY   (XSTART:XEND,YSTART:YEND) )  ! surface runoff [mm/s]
  ALLOCATE ( RUNSBXY   (XSTART:XEND,YSTART:YEND) )  ! subsurface runoff [mm/s]
  ALLOCATE ( ECANXY    (XSTART:XEND,YSTART:YEND) )  ! evaporation of intercepted water (mm/s)
  ALLOCATE ( EDIRXY    (XSTART:XEND,YSTART:YEND) )  ! soil surface evaporation rate (mm/s]
  ALLOCATE ( ETRANXY   (XSTART:XEND,YSTART:YEND) )  ! transpiration rate (mm/s)
  ALLOCATE ( FSAXY     (XSTART:XEND,YSTART:YEND) )  ! total absorbed solar radiation (w/m2)
  ALLOCATE ( FIRAXY    (XSTART:XEND,YSTART:YEND) )  ! total net longwave rad (w/m2) [+ to atm]
  ALLOCATE ( APARXY    (XSTART:XEND,YSTART:YEND) )  ! photosyn active energy by canopy (w/m2)
  ALLOCATE ( PSNXY     (XSTART:XEND,YSTART:YEND) )  ! total photosynthesis (umol co2/m2/s) [+]
  ALLOCATE ( SAVXY     (XSTART:XEND,YSTART:YEND) )  ! solar rad absorbed by veg. (w/m2)
  ALLOCATE ( SAGXY     (XSTART:XEND,YSTART:YEND) )  ! solar rad absorbed by ground (w/m2)
  ALLOCATE ( RSSUNXY   (XSTART:XEND,YSTART:YEND) )  ! sunlit leaf stomatal resistance (s/m)
  ALLOCATE ( RSSHAXY   (XSTART:XEND,YSTART:YEND) )  ! shaded leaf stomatal resistance (s/m)
  ALLOCATE ( BGAPXY    (XSTART:XEND,YSTART:YEND) )  ! between gap fraction
  ALLOCATE ( WGAPXY    (XSTART:XEND,YSTART:YEND) )  ! within gap fraction
  ALLOCATE ( TGVXY     (XSTART:XEND,YSTART:YEND) )  ! under canopy ground temperature[K]
  ALLOCATE ( TGBXY     (XSTART:XEND,YSTART:YEND) )  ! bare ground temperature [K]
  ALLOCATE ( CHVXY     (XSTART:XEND,YSTART:YEND) )  ! sensible heat exchange coefficient vegetated
  ALLOCATE ( CHBXY     (XSTART:XEND,YSTART:YEND) )  ! sensible heat exchange coefficient bare-ground
  ALLOCATE ( SHGXY     (XSTART:XEND,YSTART:YEND) )  ! veg ground sen. heat [w/m2]   [+ to atm]
  ALLOCATE ( SHCXY     (XSTART:XEND,YSTART:YEND) )  ! canopy sen. heat [w/m2]   [+ to atm]
  ALLOCATE ( SHBXY     (XSTART:XEND,YSTART:YEND) )  ! bare sensible heat [w/m2]  [+ to atm]
  ALLOCATE ( EVGXY     (XSTART:XEND,YSTART:YEND) )  ! veg ground evap. heat [w/m2]  [+ to atm]
  ALLOCATE ( EVBXY     (XSTART:XEND,YSTART:YEND) )  ! bare soil evaporation [w/m2]  [+ to atm]
  ALLOCATE ( GHVXY     (XSTART:XEND,YSTART:YEND) )  ! veg ground heat flux [w/m2]  [+ to soil]
  ALLOCATE ( GHBXY     (XSTART:XEND,YSTART:YEND) )  ! bare ground heat flux [w/m2] [+ to soil]
  ALLOCATE ( IRGXY     (XSTART:XEND,YSTART:YEND) )  ! veg ground net LW rad. [w/m2] [+ to atm]
  ALLOCATE ( IRCXY     (XSTART:XEND,YSTART:YEND) )  ! canopy net LW rad. [w/m2] [+ to atm]
  ALLOCATE ( IRBXY     (XSTART:XEND,YSTART:YEND) )  ! bare net longwave rad. [w/m2] [+ to atm]
  ALLOCATE ( TRXY      (XSTART:XEND,YSTART:YEND) )  ! transpiration [w/m2]  [+ to atm]
  ALLOCATE ( EVCXY     (XSTART:XEND,YSTART:YEND) )  ! canopy evaporation heat [w/m2]  [+ to atm]
  ALLOCATE ( CHLEAFXY  (XSTART:XEND,YSTART:YEND) )  ! leaf exchange coefficient
  ALLOCATE ( CHUCXY    (XSTART:XEND,YSTART:YEND) )  ! under canopy exchange coefficient
  ALLOCATE ( CHV2XY    (XSTART:XEND,YSTART:YEND) )  ! veg 2m exchange coefficient
  ALLOCATE ( CHB2XY    (XSTART:XEND,YSTART:YEND) )  ! bare 2m exchange coefficient
  ALLOCATE ( RS        (XSTART:XEND,YSTART:YEND) )  ! Total stomatal resistance (s/m)
  ALLOCATE ( Z0        (XSTART:XEND,YSTART:YEND) )  ! roughness length output to WRF
  ALLOCATE ( ZNT       (XSTART:XEND,YSTART:YEND) )  ! roughness length output to WRF
  ALLOCATE ( QTDRAIN   (XSTART:XEND,YSTART:YEND) )  ! tile drainage (mm)
  ALLOCATE ( TD_FRACTION (XSTART:XEND,YSTART:YEND) )! tile drainage fraction
! additional output variables
  ALLOCATE ( PAHXY     (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( PAHGXY    (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( PAHBXY    (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( PAHVXY    (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( QINTSXY   (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( QINTRXY   (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( QDRIPSXY  (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( QDRIPRXY  (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( QTHROSXY  (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( QTHRORXY  (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( QSNSUBXY  (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( QSNFROXY  (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( QSUBCXY   (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( QFROCXY   (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( QEVACXY   (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( QDEWCXY   (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( QFRZCXY   (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( QMELTCXY  (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( QSNBOTXY  (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( QMELTXY   (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( PONDINGXY (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( FPICEXY   (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( RAINLSM   (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( SNOWLSM   (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( FORCTLSM  (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( FORCQLSM  (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( FORCPLSM  (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( FORCZLSM  (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( FORCWLSM  (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( ACC_SSOILXY (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( ACC_QINSURXY(XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( ACC_QSEVAXY (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( ACC_ETRANIXY(XSTART:XEND,1:NSOIL,YSTART:YEND) )
  ALLOCATE ( EFLXBXY   (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( SOILENERGY(XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( SNOWENERGY(XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( CANHSXY   (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( ACC_DWATERXY(XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( ACC_PRCPXY(XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( ACC_ECANXY(XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( ACC_ETRANXY (XSTART:XEND,YSTART:YEND) )
  ALLOCATE ( ACC_EDIRXY(XSTART:XEND,YSTART:YEND) )

  ALLOCATE ( XLONG     (XSTART:XEND,YSTART:YEND) )  ! longitude
  ALLOCATE ( TERRAIN   (XSTART:XEND,YSTART:YEND) )  ! terrain height
  ALLOCATE ( GVFMIN    (XSTART:XEND,YSTART:YEND) )  ! annual minimum in vegetation fraction
  ALLOCATE ( GVFMAX    (XSTART:XEND,YSTART:YEND) )  ! annual maximum in vegetation fraction

!------------------------------------------------------------------------
! Needed for MMF_RUNOFF (IOPT_RUN = 5); not part of MP driver in WRF
!------------------------------------------------------------------------

  ALLOCATE ( MSFTX       (XSTART:XEND,YSTART:YEND) )  !
  ALLOCATE ( MSFTY       (XSTART:XEND,YSTART:YEND) )  !
  ALLOCATE ( EQZWT       (XSTART:XEND,YSTART:YEND) )  !
  ALLOCATE ( RIVERBEDXY  (XSTART:XEND,YSTART:YEND) )  !
  ALLOCATE ( RIVERCONDXY (XSTART:XEND,YSTART:YEND) )  !
  ALLOCATE ( PEXPXY      (XSTART:XEND,YSTART:YEND) )  !
  ALLOCATE ( FDEPTHXY    (XSTART:XEND,YSTART:YEND) )  !
  ALLOCATE ( AREAXY      (XSTART:XEND,YSTART:YEND) )  !
  ALLOCATE ( QRFSXY      (XSTART:XEND,YSTART:YEND) )  !
  ALLOCATE ( QSPRINGSXY  (XSTART:XEND,YSTART:YEND) )  !
  ALLOCATE ( QRFXY       (XSTART:XEND,YSTART:YEND) )  !
  ALLOCATE ( QSPRINGXY   (XSTART:XEND,YSTART:YEND) )  !
  ALLOCATE ( QSLATXY     (XSTART:XEND,YSTART:YEND) )  !
  ALLOCATE ( QLATXY      (XSTART:XEND,YSTART:YEND) )  !
  ALLOCATE ( RECHCLIM    (XSTART:XEND,YSTART:YEND) )  !
  ALLOCATE ( RIVERMASK   (XSTART:XEND,YSTART:YEND) )  !
  ALLOCATE ( NONRIVERXY  (XSTART:XEND,YSTART:YEND) )  !

!------------------------------------------------------------------------
! Needed for crop model (OPT_CROP=1)
!------------------------------------------------------------------------

  ALLOCATE ( PGSXY       (XSTART:XEND,  YSTART:YEND) )
  ALLOCATE ( CROPCAT     (XSTART:XEND,  YSTART:YEND) )
  ALLOCATE ( PLANTING    (XSTART:XEND,  YSTART:YEND) )
  ALLOCATE ( HARVEST     (XSTART:XEND,  YSTART:YEND) )
  ALLOCATE ( SEASON_GDD  (XSTART:XEND,  YSTART:YEND) )
  ALLOCATE ( CROPTYPE    (XSTART:XEND,5,YSTART:YEND) )

!------------------------------------------------------------------------
! Single- and Multi-layer Urban Models
!------------------------------------------------------------------------

IF(SF_URBAN_PHYSICS > 0 )  THEN  ! any urban model

  ALLOCATE ( sh_urb2d       (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( lh_urb2d       (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( g_urb2d        (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( rn_urb2d       (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( ts_urb2d       (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( HRANG          (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( DECLIN                                                    )  !
  ALLOCATE ( GMT                                                       )  !
  ALLOCATE ( JULDAY                                                    )  !
  ALLOCATE ( frc_urb2d      (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( utype_urb2d    (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( lp_urb2d       (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( lb_urb2d       (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( hgt_urb2d      (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( ust            (XSTART:XEND,                 YSTART:YEND) )  !

!ENDIF

!IF(SF_URBAN_PHYSICS == 1 ) THEN  ! single layer urban model

  ALLOCATE ( cmr_sfcdif     (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( chr_sfcdif     (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( cmc_sfcdif     (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( chc_sfcdif     (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( cmgr_sfcdif    (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( chgr_sfcdif    (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( tr_urb2d       (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( tb_urb2d       (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( tg_urb2d       (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( tc_urb2d       (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( qc_urb2d       (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( uc_urb2d       (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( xxxr_urb2d     (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( xxxb_urb2d     (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( xxxg_urb2d     (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( xxxc_urb2d     (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( trl_urb3d      (XSTART:XEND, nsoil,          YSTART:YEND) )  !
  ALLOCATE ( tbl_urb3d      (XSTART:XEND, nsoil,          YSTART:YEND) )  !
  ALLOCATE ( tgl_urb3d      (XSTART:XEND, nsoil,          YSTART:YEND) )  !

  ALLOCATE ( psim_urb2d     (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( psih_urb2d     (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( u10_urb2d      (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( v10_urb2d      (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( GZ1OZ0_urb2d   (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( AKMS_URB2D     (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( th2_urb2d      (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( q2_urb2d       (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( ust_urb2d      (XSTART:XEND,                 YSTART:YEND) )  !

  ALLOCATE ( dzr            (             nsoil                      ) )  !
  ALLOCATE ( dzb            (             nsoil                      ) )  !
  ALLOCATE ( dzg            (             nsoil                      ) )  !
  ALLOCATE ( cmcr_urb2d     (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( tgr_urb2d      (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( tgrl_urb3d     (XSTART:XEND, nsoil,          YSTART:YEND) )  !
  ALLOCATE ( smr_urb3d      (XSTART:XEND, nsoil,          YSTART:YEND) )  !
  ALLOCATE ( drelr_urb2d    (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( drelb_urb2d    (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( drelg_urb2d    (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( flxhumr_urb2d  (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( flxhumb_urb2d  (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( flxhumg_urb2d  (XSTART:XEND,                 YSTART:YEND) )  !

  ALLOCATE ( chs            (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( chs2           (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( cqs2           (XSTART:XEND,                 YSTART:YEND) )  !

  ALLOCATE ( mh_urb2d       (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( stdh_urb2d     (XSTART:XEND,                 YSTART:YEND) )  !
  ALLOCATE ( lf_urb2d       (XSTART:XEND, 4,              YSTART:YEND) )  !

!ENDIF

!IF(SF_URBAN_PHYSICS == 2 .or. SF_URBAN_PHYSICS == 3) THEN  ! BEP or BEM urban models

  ALLOCATE ( trb_urb4d      (XSTART:XEND,urban_map_zrd,YSTART:YEND) )  !
  ALLOCATE ( tw1_urb4d      (XSTART:XEND,urban_map_zwd,YSTART:YEND) )  !
  ALLOCATE ( tw2_urb4d      (XSTART:XEND,urban_map_zwd,YSTART:YEND) )  !
  ALLOCATE ( tgb_urb4d      (XSTART:XEND,urban_map_gd ,YSTART:YEND) )  !
  ALLOCATE ( sfw1_urb3d     (XSTART:XEND,urban_map_zd ,YSTART:YEND) )  !
  ALLOCATE ( sfw2_urb3d     (XSTART:XEND,urban_map_zd ,YSTART:YEND) )  !
  ALLOCATE ( sfr_urb3d      (XSTART:XEND,urban_map_zdf,YSTART:YEND) )  !
  ALLOCATE ( sfg_urb3d      (XSTART:XEND,num_urban_ndm,YSTART:YEND) )  !

  ALLOCATE ( hi_urb2d       (XSTART:XEND, num_urban_hi,   YSTART:YEND) )  !

  ALLOCATE ( theta_urban    (XSTART:XEND,KDS:KDE,         YSTART:YEND) )  !
  ALLOCATE (     u_urban    (XSTART:XEND,KDS:KDE,         YSTART:YEND) )  !
  ALLOCATE (     v_urban    (XSTART:XEND,KDS:KDE,         YSTART:YEND) )  !
  ALLOCATE (    dz_urban    (XSTART:XEND,KDS:KDE,         YSTART:YEND) )  !
  ALLOCATE (   rho_urban    (XSTART:XEND,KDS:KDE,         YSTART:YEND) )  !
  ALLOCATE (     p_urban    (XSTART:XEND,KDS:KDE,         YSTART:YEND) )  !

  ALLOCATE ( a_u_bep        (XSTART:XEND,KDS:KDE,         YSTART:YEND) )  !
  ALLOCATE ( a_v_bep        (XSTART:XEND,KDS:KDE,         YSTART:YEND) )  !
  ALLOCATE ( a_t_bep        (XSTART:XEND,KDS:KDE,         YSTART:YEND) )  !
  ALLOCATE ( a_q_bep        (XSTART:XEND,KDS:KDE,         YSTART:YEND) )  !
  ALLOCATE ( a_e_bep        (XSTART:XEND,KDS:KDE,         YSTART:YEND) )  !
  ALLOCATE ( b_u_bep        (XSTART:XEND,KDS:KDE,         YSTART:YEND) )  !
  ALLOCATE ( b_v_bep        (XSTART:XEND,KDS:KDE,         YSTART:YEND) )  !
  ALLOCATE ( b_t_bep        (XSTART:XEND,KDS:KDE,         YSTART:YEND) )  !
  ALLOCATE ( b_q_bep        (XSTART:XEND,KDS:KDE,         YSTART:YEND) )  !
  ALLOCATE ( b_e_bep        (XSTART:XEND,KDS:KDE,         YSTART:YEND) )  !
  ALLOCATE ( dlg_bep        (XSTART:XEND,KDS:KDE,         YSTART:YEND) )  !
  ALLOCATE ( dl_u_bep       (XSTART:XEND,KDS:KDE,         YSTART:YEND) )  !
  ALLOCATE ( sf_bep         (XSTART:XEND,KDS:KDE,         YSTART:YEND) )  !
  ALLOCATE ( vl_bep         (XSTART:XEND,KDS:KDE,         YSTART:YEND) )  !

!ENDIF

!IF(SF_URBAN_PHYSICS == 3) THEN  ! BEM urban model

  ALLOCATE ( tlev_urb3d     (XSTART:XEND,urban_map_bd ,YSTART:YEND) )  !
  ALLOCATE ( qlev_urb3d     (XSTART:XEND,urban_map_bd ,YSTART:YEND) )  !
  ALLOCATE ( tw1lev_urb3d   (XSTART:XEND,urban_map_wd ,YSTART:YEND) )  !
  ALLOCATE ( tw2lev_urb3d   (XSTART:XEND,urban_map_wd ,YSTART:YEND) )  !
  ALLOCATE ( tglev_urb3d    (XSTART:XEND,urban_map_gbd,YSTART:YEND) )  !
  ALLOCATE ( tflev_urb3d    (XSTART:XEND,urban_map_fbd,YSTART:YEND) )  !
  ALLOCATE ( sf_ac_urb3d    (XSTART:XEND,              YSTART:YEND) )  !
  ALLOCATE ( lf_ac_urb3d    (XSTART:XEND,              YSTART:YEND) )  !
  ALLOCATE ( cm_ac_urb3d    (XSTART:XEND,              YSTART:YEND) )  !
  ALLOCATE ( sfvent_urb3d   (XSTART:XEND,              YSTART:YEND) )  !
  ALLOCATE ( lfvent_urb3d   (XSTART:XEND,              YSTART:YEND) )  !
  ALLOCATE ( sfwin1_urb3d   (XSTART:XEND,urban_map_wd ,YSTART:YEND) )  !
  ALLOCATE ( sfwin2_urb3d   (XSTART:XEND,urban_map_wd ,YSTART:YEND) )  !
! new urban variables greenroof & solar panel for BEM
  ALLOCATE ( ep_pv_urb3d    (XSTART:XEND,              YSTART:YEND) )  !
  ALLOCATE ( t_pv_urb3d     (XSTART:XEND,urban_map_zdf,YSTART:YEND) )  !
  ALLOCATE ( trv_urb4d      (XSTART:XEND,urban_map_zgrd,YSTART:YEND) ) !
  ALLOCATE ( qr_urb4d       (XSTART:XEND,urban_map_zgrd,YSTART:YEND) ) !
  ALLOCATE ( qgr_urb3d      (XSTART:XEND,              YSTART:YEND) )  !
  ALLOCATE ( tgr_urb3d      (XSTART:XEND,              YSTART:YEND) )  !
  ALLOCATE ( drain_urb4d    (XSTART:XEND,urban_map_zdf,YSTART:YEND) )  !
  ALLOCATE ( draingr_urb3d  (XSTART:XEND,              YSTART:YEND) )  !
  ALLOCATE ( sfrv_urb3d     (XSTART:XEND,urban_map_zdf,YSTART:YEND) )  !
  ALLOCATE ( lfrv_urb3d     (XSTART:XEND,urban_map_zdf,YSTART:YEND) )  !
  ALLOCATE ( dgr_urb3d      (XSTART:XEND,urban_map_zdf,YSTART:YEND) )  !
  ALLOCATE ( dg_urb3d       (XSTART:XEND,num_urban_ndm,YSTART:YEND) )  !
  ALLOCATE ( lfr_urb3d      (XSTART:XEND,urban_map_zdf,YSTART:YEND) )  !
  ALLOCATE ( lfg_urb3d      (XSTART:XEND,num_urban_ndm,YSTART:YEND) )  !

ENDIF

!------------------------------------------------------------------------

  ALLOCATE ( CHSTARXY  (XSTART:XEND,YSTART:YEND) )  ! for consistency with MP_init; delete later
  ALLOCATE ( SEAICE    (XSTART:XEND,YSTART:YEND) )  ! seaice fraction

!------------------------------------------------------------------------
!Added by Aaron Alexander for noahmp mosaic scheme
!------------------------------------------------------------------------

IF(IOPT_MOSAIC == 1) THEN
  ALLOCATE ( LANDUSEF         (XSTART:XEND,1:NUMBER_LAND_USE_CATAGORIES,YSTART:YEND) ) ! Fractioanl land-use data allocated. Written by Aaron A.
  ALLOCATE ( LANDUSEF2        (XSTART:XEND,1:NUMBER_LAND_USE_CATAGORIES,YSTART:YEND) ) !This is the landuse fraction that is normalized after the initlization (and will be used to average)
  ALLOCATE ( mosaic_cat_index (XSTART:XEND,1:NUMBER_LAND_USE_CATAGORIES,YSTART:YEND) ) !This is the actual catagory (for NLCD: 1 - 40). Used to call correct land use parameters
                                                                                       !

  ALLOCATE ( TSK_mosaic       (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! surface radiative temperature [K]
  ALLOCATE ( HFX_mosaic       (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! sensible heat flux [W m-2]
  ALLOCATE ( QFX_mosaic       (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! latent heat flux [kg s-1 m-2]
  ALLOCATE ( LH_mosaic        (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! latent heat flux [W m-2]
  ALLOCATE ( TMN_mosaic       (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! latent heat flux [W m-2]
  ALLOCATE ( GRDFLX_mosaic    (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! ground/snow heat flux [W m-2]
  ALLOCATE ( SFCRUNOFF_mosaic (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! accumulated surface runoff [m]
  ALLOCATE ( UDRUNOFF_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! accumulated sub-surface runoff [m]
  ALLOCATE ( ALBEDO_mosaic    (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! total grid albedo []
  ALLOCATE ( SNOWC_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! snow cover fraction []
  ALLOCATE ( CANWAT_mosaic    (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! total canopy water + ice [mm]
  ALLOCATE ( SNOW_mosaic      (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! snow water equivalent [mm]
  ALLOCATE ( SNOWH_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! physical snow depth [m]
  ALLOCATE ( ACSNOM_mosaic    (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! accumulated snow melt leaving pack
  ALLOCATE ( ACSNOW_mosaic    (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! accumulated snow on grid
  ALLOCATE ( EMISS_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! surface bulk emissivity
  ALLOCATE ( QSFC_mosaic      (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! bulk surface specific humidity
  ALLOCATE ( Z0_mosaic        (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! combined z0 sent to coupled model
  ALLOCATE ( ZNT_mosaic       (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! combined z0 sent to coupled model

  ALLOCATE ( tvxy_mosaic      (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! vegetation leaf temperature
  ALLOCATE ( tgxy_mosaic      (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! bulk ground surface temperature
  ALLOCATE ( canicexy_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! canopy-intercepted ice (mm)
  ALLOCATE ( canliqxy_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! canopy-intercepted liquid water (mm)
  ALLOCATE ( eahxy_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! canopy air vapor pressure (pa)
  ALLOCATE ( tahxy_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! canopy air temperature (k)
  ALLOCATE ( cmxy_mosaic      (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! bulk momentum drag coefficient
  ALLOCATE ( chxy_mosaic      (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! bulk sensible heat exchange coefficient
  ALLOCATE ( fwetxy_mosaic    (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! wetted or snowed fraction of the canopy (-)
  ALLOCATE ( sneqvoxy_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! snow mass at last time step(mm h2o)
  ALLOCATE ( alboldxy_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! snow albedo at last time step (-)
  ALLOCATE ( qsnowxy_mosaic   (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! snowfall on the ground [mm/s]
  ALLOCATE ( qrainxy_mosaic   (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! rainfall on the ground [mm/s]
  ALLOCATE ( wslakexy_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! lake water storage [mm]
  ALLOCATE ( zwtxy_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! water table depth [m]
  ALLOCATE ( waxy_mosaic      (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! water in the "aquifer" [mm]
  ALLOCATE ( wtxy_mosaic      (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! groundwater storage [mm]
  ALLOCATE ( lfmassxy_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! leaf mass [g/m2]
  ALLOCATE ( rtmassxy_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! mass of fine roots [g/m2]
  ALLOCATE ( stmassxy_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! stem mass [g/m2]
  ALLOCATE ( woodxy_mosaic    (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! mass of wood (incl. woody roots) [g/m2]
  ALLOCATE ( grainxy_mosaic   (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! mass of grain XING [g/m2]
  ALLOCATE ( gddxy_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! growing degree days XING (based on 10C)
  ALLOCATE ( pgsxy_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( smcwtdxy_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( stblcpxy_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! stable carbon in deep soil [g/m2]
  ALLOCATE ( fastcpxy_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! short-lived carbon, shallow soil [g/m2]
  ALLOCATE ( xsaixy_mosaic    (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! stem area index
  ALLOCATE ( xlai_mosaic      (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! leaf area index
  ALLOCATE ( taussxy_mosaic   (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( smcwtdxy_mosiac  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( deeprechxy_mosaic(XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( rechxy_mosaic    (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )

  ! irrigation
  ALLOCATE ( IRFRACT_mosaic   (XSTART:XEND,1:number_mosaic_catagories,YSTART:YEND) )  ! irrigation fraction
  ALLOCATE ( SIFRACT_mosaic   (XSTART:XEND,1:number_mosaic_catagories,YSTART:YEND) )  ! sprinkler irrigation fraction
  ALLOCATE ( MIFRACT_mosaic   (XSTART:XEND,1:number_mosaic_catagories,YSTART:YEND) )  ! micro irrigation fraction
  ALLOCATE ( FIFRACT_mosaic   (XSTART:XEND,1:number_mosaic_catagories,YSTART:YEND) )  ! flood irrigation fraction
  ALLOCATE ( IRNUMSI_mosaic   (XSTART:XEND,1:number_mosaic_catagories,YSTART:YEND) )  ! irrigation event number, Sprinkler
  ALLOCATE ( IRNUMMI_mosaic   (XSTART:XEND,1:number_mosaic_catagories,YSTART:YEND) )  ! irrigation event number, Micro
  ALLOCATE ( IRNUMFI_mosaic   (XSTART:XEND,1:number_mosaic_catagories,YSTART:YEND) )  ! irrigation event number, Flood
  ALLOCATE ( IRWATSI_mosaic   (XSTART:XEND,1:number_mosaic_catagories,YSTART:YEND) )  ! irrigation water amount [m] to be applied, Sprinkler
  ALLOCATE ( IRWATMI_mosaic   (XSTART:XEND,1:number_mosaic_catagories,YSTART:YEND) )  ! irrigation water amount [m] to be applied, Micro
  ALLOCATE ( IRWATFI_mosaic   (XSTART:XEND,1:number_mosaic_catagories,YSTART:YEND) )  ! irrigation water amount [m] to be applied, Flood
  ALLOCATE ( IRELOSS_mosaic   (XSTART:XEND,1:number_mosaic_catagories,YSTART:YEND) )  ! loss of irrigation water to evaporation,sprinkler [mm]
  ALLOCATE ( IRSIVOL_mosaic   (XSTART:XEND,1:number_mosaic_catagories,YSTART:YEND) )  ! amount of irrigation by sprinkler (mm)
  ALLOCATE ( IRMIVOL_mosaic   (XSTART:XEND,1:number_mosaic_catagories,YSTART:YEND) )  ! amount of irrigation by micro (mm)
  ALLOCATE ( IRFIVOL_mosaic   (XSTART:XEND,1:number_mosaic_catagories,YSTART:YEND) )  ! amount of irrigation by micro (mm)
  ALLOCATE ( IRRSPLH_mosaic   (XSTART:XEND,1:number_mosaic_catagories,YSTART:YEND) )  ! latent heating from sprinkler evaporation (w/m2)
  ALLOCATE ( LOCTIM_mosaic    (XSTART:XEND,1:number_mosaic_catagories,YSTART:YEND) )  ! local time

  ALLOCATE ( t2mvxy_mosaic    (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! 2m temperature of vegetation part
  ALLOCATE ( t2mbxy_mosaic    (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! 2m temperature of bare ground part
  ALLOCATE ( chstarxy_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! dummy mosaic [not used]
  ALLOCATE ( q2mvxy_mosaic    (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! 2m mixing ratio of vegetation part
  ALLOCATE ( q2mbxy_mosaic    (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! 2m mixing ratio of bare ground part
  ALLOCATE ( tradxy_mosaic    (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! surface radiative temperature (k)
  ALLOCATE ( neexy_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! net ecosys exchange (g/m2/s CO2)
  ALLOCATE ( gppxy_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! gross primary assimilation [g/m2/s C]
  ALLOCATE ( nppxy_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! net primary productivity [g/m2/s C]
  ALLOCATE ( fvegxy_mosaic    (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! Noah-MP vegetation fraction [-]
  ALLOCATE ( runsfxy_mosaic   (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! surface runoff [mm/s]
  ALLOCATE ( runsbxy_mosaic   (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! subsurface runoff [mm/s]
  ALLOCATE ( ecanxy_mosaic    (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! evaporation of intercepted water (mm/s)
  ALLOCATE ( edirxy_mosaic    (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! soil surface evaporation rate (mm/s]
  ALLOCATE ( etranxy_mosaic   (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! transpiration rate (mm/s)
  ALLOCATE ( fsaxy_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! total absorbed solar radiation (w/m2)
  ALLOCATE ( firaxy_mosaic    (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! total net longwave rad (w/m2) [+ to atm]
  ALLOCATE ( aparxy_mosaic    (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! photosyn active energy by canopy (w/m2)
  ALLOCATE ( psnxy_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! total photosynthesis (umol co2/m2/s) [+]
  ALLOCATE ( savxy_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! solar rad absorbed by veg. (w/m2)
  ALLOCATE ( sagxy_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! solar rad absorbed by ground (w/m2)
  ALLOCATE ( rssunxy_mosaic   (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! sunlit leaf stomatal resistance (s/m)
  ALLOCATE ( rsshaxy_mosaic   (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! shaded leaf stomatal resistance (s/m)
  ALLOCATE ( bgapxy_mosaic    (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! between gap fraction
  ALLOCATE ( wgapxy_mosaic    (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! within gap fraction
  ALLOCATE ( tgvxy_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! under canopy ground temperature[K]
  ALLOCATE ( tgbxy_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! bare ground temperature [K]
  ALLOCATE ( chvxy_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! sensible heat exchange coefficient vegetated
  ALLOCATE ( chbxy_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! sensible heat exchange coefficient bare-ground
  ALLOCATE ( shgxy_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! veg ground sen. heat [w/m2]   [+ to atm]
  ALLOCATE ( shcxy_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! canopy sen. heat [w/m2]   [+ to atm]
  ALLOCATE ( shbxy_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! bare sensible heat [w/m2]     [+ to atm]
  ALLOCATE ( evgxy_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! veg ground evap. heat [w/m2]  [+ to atm]
  ALLOCATE ( evbxy_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! bare soil evaporation [w/m2]  [+ to atm]
  ALLOCATE ( ghvxy_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! veg ground heat flux [w/m2]  [+ to soil]
  ALLOCATE ( ghbxy_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! bare ground heat flux [w/m2] [+ to soil]
  ALLOCATE ( irgxy_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! veg ground net LW rad. [w/m2] [+ to atm]
  ALLOCATE ( ircxy_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! canopy net LW rad. [w/m2] [+ to atm]
  ALLOCATE ( irbxy_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! bare net longwave rad. [w/m2] [+ to atm]
  ALLOCATE ( trxy_mosaic      (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! transpiration [w/m2]  [+ to atm]
  ALLOCATE ( evcxy_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! canopy evaporation heat [w/m2]  [+ to atm]
  ALLOCATE ( chleafxy_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! leaf exchange coefficient
  ALLOCATE ( chucxy_mosaic    (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! under canopy exchange coefficient
  ALLOCATE ( chv2xy_mosaic    (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! veg 2m exchange coefficient
  ALLOCATE ( chb2xy_mosaic    (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! bare 2m exchange coefficient
  ALLOCATE ( rs_mosaic        (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) ) ! Total stomatal resistance (s/m)
  ! The other outputs:
  ! additional output variables (unsure if these are needed)

  ALLOCATE ( PAHXY_mosaic     (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( PAHGXY_mosaic    (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( PAHBXY_mosaic    (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( PAHVXY_mosaic    (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( QINTSXY_mosaic   (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( QINTRXY_mosaic   (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( QDRIPSXY_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( QDRIPRXY_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( QTHROSXY_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( QTHRORXY_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( QSNSUBXY_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( QSNFROXY_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( QSUBCXY_mosaic   (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( QFROCXY_mosaic   (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( QEVACXY_mosaic   (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( QDEWCXY_mosaic   (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( QFRZCXY_mosaic   (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( QMELTCXY_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( QSNBOTXY_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( QMELTXY_mosaic   (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( PONDINGXY_mosaic (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( FPICEXY_mosaic   (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( ACC_SSOILXY_mosaic (XSTART:XEND,1:number_mosaic_catagories,YSTART:YEND) )
  ALLOCATE ( ACC_QINSURXY_mosaic(XSTART:XEND,1:number_mosaic_catagories,YSTART:YEND) )
  ALLOCATE ( ACC_QSEVAXY_mosaic (XSTART:XEND,1:number_mosaic_catagories,YSTART:YEND) )
  ALLOCATE ( ACC_ETRANIXY_mosaic(XSTART:XEND,1:NSOIL*number_mosaic_catagories  ,YSTART:YEND) )
  ALLOCATE ( EFLXBXY_mosaic   (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( SOILENERGY_mosaic(XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( SNOWENERGY_mosaic(XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( CANHSXY_mosaic   (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( ACC_DWATERXY_mosaic(XSTART:XEND,1:number_mosaic_catagories,YSTART:YEND) )
  ALLOCATE ( ACC_PRCPXY_mosaic(XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( ACC_ECANXY_mosaic(XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( ACC_ETRANXY_mosaic (XSTART:XEND,1:number_mosaic_catagories,YSTART:YEND) )
  ALLOCATE ( ACC_EDIRXY_mosaic(XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )

  ! Snow variables of interest
  ALLOCATE ( isnowxy_mosaic   (XSTART:XEND,1:number_mosaic_catagories,YSTART:YEND) )   ! actual no. of snow layers
  ALLOCATE ( zsnsoxy_mosaic   (XSTART:XEND,1:7*number_mosaic_catagories,YSTART:YEND) ) ! snow layer depth [m]
  ALLOCATE ( tsnoxy_mosaic    (XSTART:XEND,1:3*number_mosaic_catagories,YSTART:YEND) ) ! snow temperature [K]
  ALLOCATE ( snicexy_mosaic   (XSTART:XEND,1:3*number_mosaic_catagories,YSTART:YEND) ) ! snow layer ice [mm]
  ALLOCATE ( snliqxy_mosaic   (XSTART:XEND,1:3*number_mosaic_catagories,YSTART:YEND) ) ! snow layer liquid water [mm]

  ! Soil variables
  ALLOCATE ( TSLB_mosaic      (XSTART:XEND,1:NSOIL*number_mosaic_catagories, YSTART:YEND) ) ! soil temperature [K]
  ALLOCATE ( SMOIS_mosaic     (XSTART:XEND,1:NSOIL*number_mosaic_catagories, YSTART:YEND) ) ! volumetric soil moisture [m3/m3]
  ALLOCATE ( SH2O_mosaic      (XSTART:XEND,1:NSOIL*number_mosaic_catagories, YSTART:YEND) ) ! volumetric liquid soil moisture [m3/m3]
  ALLOCATE ( SMOISEQ_mosaic   (XSTART:XEND,1:NSOIL*number_mosaic_catagories, YSTART:YEND) )

  !urban variables
  ALLOCATE ( TR_URB2D_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( TB_URB2D_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( TG_URB2D_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( TC_URB2D_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( UC_URB2D_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( QC_URB2D_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( SH_URB2D_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( LH_URB2D_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( G_URB2D_mosaic   (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( RN_URB2D_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( TS_URB2D_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )

  ALLOCATE ( TRL_URB3D_mosaic  (XSTART:XEND,1:NSOIL*number_mosaic_catagories, YSTART:YEND) )
  ALLOCATE ( TBL_URB3D_mosaic  (XSTART:XEND,1:NSOIL*number_mosaic_catagories, YSTART:YEND) )
  ALLOCATE ( TGL_URB3D_mosaic  (XSTART:XEND,1:NSOIL*number_mosaic_catagories, YSTART:YEND) )

  ! Extra variables needed for the restart of the urban model
  ALLOCATE ( CMR_SFCDIF_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( CHR_SFCDIF_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( CMC_SFCDIF_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( CHC_SFCDIF_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( CMGR_SFCDIF_mosaic (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( CHGR_SFCDIF_mosaic (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( XXXR_URB2D_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( XXXB_URB2D_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( XXXG_URB2D_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( XXXC_URB2D_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( CMCR_URB2D_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( TGR_URB2D_mosaic   (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( TGRL_URB3D_mosaic  (XSTART:XEND,1:NSOIL*number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( SMR_URB3D_mosaic   (XSTART:XEND,1:NSOIL*number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( DRELR_URB2D_mosaic (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( DRELB_URB2D_mosaic (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( DRELG_URB2D_mosaic (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( FLXHUMR_URB2D_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( FLXHUMB_URB2D_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( FLXHUMG_URB2D_mosaic  (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )

  ALLOCATE ( RUNONSFXY (XSTART:XEND,  YSTART:YEND) )
  ALLOCATE ( RUNONSFXY_mosaic (XSTART:XEND,1:number_mosaic_catagories,YSTART:YEND) )
!IF (IOPT_HUE.eq.1) THEN

  !ALLOCATE ( RUNONSFXY (XSTART:XEND,  YSTART:YEND) )
  !ALLOCATE ( RUNONSFXY_mosaic (XSTART:XEND,1:number_mosaic_catagories,YSTART:YEND) )
  ALLOCATE ( DETENTION_STORAGEXY (XSTART:XEND,  YSTART:YEND) )
  ALLOCATE ( DETENTION_STORAGEXY_mosaic (XSTART:XEND,1:number_mosaic_catagories,YSTART:YEND) )
  ALLOCATE ( VOL_FLUX_RUNONXY_mosaic (XSTART:XEND,1:number_mosaic_catagories,  YSTART:YEND) )
  ALLOCATE ( VOL_FLUX_SMXY_mosaic (XSTART:XEND,1:NSOIL*number_mosaic_catagories,  YSTART:YEND) )

!end if

END IF
!------------------------------------------------------------------------
!Added by Aaron Alexander for noahmp mosaic scheme
!------------------------------------------------------------------------


#ifdef WRF_HYDRO
   allocate( greenfrac   (XSTART:XEND,YSTART:YEND))
   greenfrac = 0
#endif

  COSZEN     = undefined_real
  XLAT       = undefined_real
  DZ8W       = undefined_real
  DZS        = undefined_real
  IVGTYP     = undefined_int
  ISLTYP     = undefined_int
  SOILCL1    = undefined_real
  SOILCL2    = undefined_real
  SOILCL3    = undefined_real
  SOILCL4    = undefined_real
  SOILCOMP   = undefined_real
  VEGFRA     = undefined_real
  GVFMAX     = undefined_real
  TMN        = undefined_real
  XLAND      = undefined_real
  XICE       = undefined_real
  T_PHY      = undefined_real
  QV_CURR    = undefined_real
  U_PHY      = undefined_real
  V_PHY      = undefined_real
  SWDOWN     = undefined_real
  SWDDIR     = undefined_real
  SWDDIF     = undefined_real
  GLW        = undefined_real
  P8W        = undefined_real
  RAINBL     = undefined_real
  SNOWBL     = undefined_real
  RAINBL_tmp = undefined_real
  SR         = undefined_real
  RAINCV     = undefined_real
  RAINNCV    = undefined_real
  RAINSHV    = undefined_real
  SNOWNCV    = undefined_real
  GRAUPELNCV = undefined_real
  HAILNCV    = undefined_real
  TSK        = undefined_real
  QFX        = undefined_real
  SMSTAV     = undefined_real
  SMSTOT     = undefined_real
  SMOIS      = undefined_real
  SH2O       = undefined_real
  TSLB       = undefined_real
  SNOW       = undefined_real
  SNOWH      = undefined_real
  CANWAT     = undefined_real
  ACSNOM     = 0.0
  ACSNOW     = 0.0
  QSFC       = undefined_real
  SFCRUNOFF  = 0.0
  UDRUNOFF   = 0.0
  SMOISEQ    = undefined_real
  ALBEDO     = undefined_real
  ISNOWXY    = undefined_int
  TVXY       = undefined_real
  TGXY       = undefined_real
  CANICEXY   = undefined_real
  CANLIQXY   = undefined_real
  EAHXY      = undefined_real
  TAHXY      = undefined_real
  CMXY       = undefined_real
  CHXY       = undefined_real
  FWETXY     = undefined_real
  SNEQVOXY   = undefined_real
  ALBOLDXY   = undefined_real
  QSNOWXY    = undefined_real
  QRAINXY    = undefined_real
  WSLAKEXY   = undefined_real
  ZWTXY      = undefined_real
  WAXY       = undefined_real
  WTXY       = undefined_real
  TSNOXY     = undefined_real
  SNICEXY    = undefined_real
  SNLIQXY    = undefined_real
  LFMASSXY   = undefined_real
  RTMASSXY   = undefined_real
  STMASSXY   = undefined_real
  WOODXY     = undefined_real
  STBLCPXY   = undefined_real
  FASTCPXY   = undefined_real
  LAI        = undefined_real
  LAI_tmp    = undefined_real
  XSAIXY     = undefined_real
  TAUSSXY    = undefined_real
  XLONG      = undefined_real
  SEAICE     = undefined_real
  SMCWTDXY   = undefined_real
  DEEPRECHXY = 0.0
  RECHXY     = 0.0
  ZSNSOXY    = undefined_real
  GRDFLX     = undefined_real
  HFX        = undefined_real
  LH         = undefined_real
  EMISS      = undefined_real
  SNOWC      = undefined_real
  T2MVXY     = undefined_real
  T2MBXY     = undefined_real
  Q2MVXY     = undefined_real
  Q2MBXY     = undefined_real
  TRADXY     = undefined_real
  NEEXY      = undefined_real
  GPPXY      = undefined_real
  NPPXY      = undefined_real
  FVEGXY     = undefined_real
  RUNSFXY    = undefined_real
  RUNSBXY    = undefined_real
  ECANXY     = undefined_real
  EDIRXY     = undefined_real
  ETRANXY    = undefined_real
  FSAXY      = undefined_real
  FIRAXY     = undefined_real
  APARXY     = undefined_real
  PSNXY      = undefined_real
  SAVXY      = undefined_real
  FIRAXY     = undefined_real
  SAGXY      = undefined_real
  RSSUNXY    = undefined_real
  RSSHAXY    = undefined_real
  BGAPXY     = undefined_real
  WGAPXY     = undefined_real
  TGVXY      = undefined_real
  TGBXY      = undefined_real
  CHVXY      = undefined_real
  CHBXY      = undefined_real
  SHGXY      = undefined_real
  SHCXY      = undefined_real
  SHBXY      = undefined_real
  EVGXY      = undefined_real
  EVBXY      = undefined_real
  GHVXY      = undefined_real
  GHBXY      = undefined_real
  IRGXY      = undefined_real
  IRCXY      = undefined_real
  IRBXY      = undefined_real
  TRXY       = undefined_real
  EVCXY      = undefined_real
  CHLEAFXY   = undefined_real
  CHUCXY     = undefined_real
  CHV2XY     = undefined_real
  CHB2XY     = undefined_real
  RS         = undefined_real
! additional output
  PAHXY      = undefined_real
  PAHGXY     = undefined_real
  PAHBXY     = undefined_real
  PAHVXY     = undefined_real
  QINTSXY    = undefined_real
  QINTRXY    = undefined_real
  QDRIPSXY   = undefined_real
  QDRIPRXY   = undefined_real
  QTHROSXY   = undefined_real
  QTHRORXY   = undefined_real
  QSNSUBXY   = undefined_real
  QSNFROXY   = undefined_real
  QSUBCXY    = undefined_real
  QFROCXY    = undefined_real
  QEVACXY    = undefined_real
  QDEWCXY    = undefined_real
  QFRZCXY    = undefined_real
  QMELTCXY   = undefined_real
  QSNBOTXY   = undefined_real
  QMELTXY    = undefined_real
  PONDINGXY  = 0.0
  FPICEXY    = undefined_real
  RAINLSM    = undefined_real
  SNOWLSM    = undefined_real
  FORCTLSM   = undefined_real
  FORCQLSM   = undefined_real
  FORCPLSM   = undefined_real
  FORCZLSM   = undefined_real
  FORCWLSM   = undefined_real
  ACC_SSOILXY   = 0.0
  ACC_QINSURXY  = 0.0
  ACC_QSEVAXY   = 0.0
  ACC_ETRANIXY  = 0.0
  EFLXBXY    = undefined_real
  SOILENERGY = 0.0
  SNOWENERGY = 0.0
  CANHSXY    = undefined_real
  ACC_DWATERXY = 0.0
  ACC_PRCPXY   = 0.0
  ACC_ECANXY   = 0.0
  ACC_ETRANXY  = 0.0
  ACC_EDIRXY   = 0.0

  TERRAIN    = undefined_real
  GVFMIN     = undefined_real
  GVFMAX     = undefined_real
  MSFTX      = undefined_real
  MSFTY      = undefined_real
  EQZWT      = undefined_real
  RIVERBEDXY = undefined_real
  RIVERCONDXY= undefined_real
  PEXPXY     = undefined_real
  FDEPTHXY   = undefined_real
  AREAXY     = undefined_real
  QRFSXY     = undefined_real
  QSPRINGSXY = undefined_real
  QRFXY      = undefined_real
  QSPRINGXY  = undefined_real
  QSLATXY    = undefined_real
  QLATXY     = undefined_real
  CHSTARXY   = undefined_real
  Z0         = undefined_real
  ZNT        = undefined_real
  PGSXY      = undefined_int
  CROPCAT    = undefined_int
  PLANTING   = undefined_real
  HARVEST    = undefined_real
  SEASON_GDD = undefined_real
  CROPTYPE   = undefined_real
! tile drainage
  QTDRAIN    = undefined_real
  TD_FRACTION= undefined_real
! irrigation
  IRFRACT    = 0.0
  SIFRACT    = 0.0
  MIFRACT    = 0.0
  FIFRACT    = 0.0
  IRNUMSI    = 0
  IRNUMMI    = 0
  IRNUMFI    = 0
  IRWATSI    = 0.0
  IRWATMI    = 0.0
  IRWATFI    = 0.0
  IRELOSS    = 0.0
  IRSIVOL    = 0.0
  IRMIVOL    = 0.0
  IRFIVOL    = 0.0
  IRRSPLH    = 0.0
  LOCTIM     = undefined_real

  ! NOAH-MP HUE: Added by Aaron A.
  LANDUSEF   = undefined_real !added by Aaron A.

  ! Urban models
IF(SF_URBAN_PHYSICS > 0 )  THEN  ! any urban model
  HRANG      = undefined_real
  DECLIN     = undefined_real
  sh_urb2d   = undefined_real
  lh_urb2d   = undefined_real
  g_urb2d    = undefined_real
  rn_urb2d   = undefined_real
  ts_urb2d   = undefined_real
  GMT        = undefined_real
  JULDAY     = undefined_int
  frc_urb2d  = undefined_real
  utype_urb2d= undefined_int
  lp_urb2d   = undefined_real
  lb_urb2d   = undefined_real
  hgt_urb2d  = undefined_real
  ust        = undefined_real
  cmr_sfcdif = 0.
  chr_sfcdif = 0.
  cmc_sfcdif = 0.
  chc_sfcdif = 0.
  cmgr_sfcdif= 0.
  chgr_sfcdif= 0.
  tr_urb2d   = undefined_real
  tb_urb2d   = undefined_real
  tg_urb2d   = undefined_real
  tc_urb2d   = undefined_real
  qc_urb2d   = undefined_real
  uc_urb2d   = undefined_real
  xxxr_urb2d = undefined_real
  xxxb_urb2d = undefined_real
  xxxg_urb2d = undefined_real
  xxxc_urb2d = undefined_real
  trl_urb3d  = undefined_real
  tbl_urb3d  = undefined_real
  tgl_urb3d  = undefined_real
  psim_urb2d = undefined_real
  psih_urb2d = undefined_real
  u10_urb2d  = undefined_real
  v10_urb2d  = undefined_real
  GZ1OZ0_urb2d = undefined_real
  AKMS_URB2D = undefined_real
  th2_urb2d  = undefined_real
  q2_urb2d   = undefined_real
  ust_urb2d  = undefined_real
  dzr        = undefined_real
  dzb        = undefined_real
  dzg        = undefined_real
  cmcr_urb2d = undefined_real
  tgr_urb2d  = undefined_real
  tgrl_urb3d = undefined_real
  smr_urb3d  = undefined_real
  drelr_urb2d= undefined_real
  drelb_urb2d= undefined_real
  drelg_urb2d= undefined_real
  flxhumr_urb2d = undefined_real
  flxhumb_urb2d = undefined_real
  flxhumg_urb2d = undefined_real
  chs        = undefined_real
  chs2       = undefined_real
  cqs2       = undefined_real
  mh_urb2d   = undefined_real
  stdh_urb2d = undefined_real
  lf_urb2d   = undefined_real
  trb_urb4d  = undefined_real
  tw1_urb4d  = undefined_real
  tw2_urb4d  = undefined_real
  tgb_urb4d  = undefined_real
  sfw1_urb3d = undefined_real
  sfw2_urb3d = undefined_real
  sfr_urb3d  = undefined_real
  sfg_urb3d  = undefined_real
  hi_urb2d   = undefined_real
  theta_urban= undefined_real
  u_urban    = undefined_real
  v_urban    = undefined_real
  dz_urban   = undefined_real
  rho_urban  = undefined_real
  p_urban    = undefined_real
  a_u_bep    = undefined_real
  a_v_bep    = undefined_real
  a_t_bep    = undefined_real
  a_q_bep    = undefined_real
  a_e_bep    = undefined_real
  b_u_bep    = undefined_real
  b_v_bep    = undefined_real
  b_t_bep    = undefined_real
  b_q_bep    = undefined_real
  b_e_bep    = undefined_real
  dlg_bep    = undefined_real
  dl_u_bep   = undefined_real
  sf_bep     = undefined_real
  vl_bep     = undefined_real
  tlev_urb3d = undefined_real
  qlev_urb3d = undefined_real
  tw1lev_urb3d = undefined_real
  tw2lev_urb3d = undefined_real
  tglev_urb3d= undefined_real
  tflev_urb3d= undefined_real
  sf_ac_urb3d= undefined_real
  lf_ac_urb3d= undefined_real
  cm_ac_urb3d= undefined_real
  sfvent_urb3d = undefined_real
  lfvent_urb3d = undefined_real
  sfwin1_urb3d = undefined_real
  sfwin2_urb3d = undefined_real
  ep_pv_urb3d= undefined_real
  t_pv_urb3d = undefined_real
  trv_urb4d  = undefined_real
  qr_urb4d   = undefined_real
  qgr_urb3d  = undefined_real
  tgr_urb3d  = undefined_real
  drain_urb4d= undefined_real
  draingr_urb3d = undefined_real
  sfrv_urb3d = undefined_real
  lfrv_urb3d = undefined_real
  dgr_urb3d  = undefined_real
  dg_urb3d   = undefined_real
  lfr_urb3d  = undefined_real
  lfg_urb3d  = undefined_real

ENDIF

!IF (IOPT_HUE.eq.1) THEN
  RUNONSFXY = undefined_real
!END IF

  XLAND          = 1.0   ! water = 2.0, land = 1.0
  XICE           = 0.0   ! fraction of grid that is seaice
  XICE_THRESHOLD = 0.5   ! fraction of grid determining seaice (from WRF)

!----------------------------------------------------------------------
! Read Landuse Type and Soil Texture and Other Information
!----------------------------------------------------------------------
  ! Modified to add the NOAH-MP HUE
  CALL READLAND_HRLDAS(HRLDAS_SETUP_FILE, XSTART, XEND, YSTART, YEND,&
     ISWATER, ISLAKE, IVGTYP, ISLTYP, TERRAIN, TMN, XLAT, XLONG, XLAND, SEAICE,MSFTX,MSFTY, &
     IOPT_MOSAIC, geogrid_file_name_for_mosaic,number_land_use_catagories,LANDUSEF) !Added by Aaron A.

  WHERE(SEAICE > 0.0) XICE = 1.0

!------------------------------------------------------------------------
! For spatially-varying soil parameters, read in necessary extra fields
!------------------------------------------------------------------------

  if (soil_data_option == 2) then
    CALL READ_SOIL_TEXTURE(HRLDAS_SETUP_FILE, XSTART, XEND, YSTART, YEND, &
                      NSOIL,IVGTYP,SOILCL1,SOILCL2,SOILCL3,SOILCL4,ISICE,ISWATER)

    ISLTYP = nint(SOILCL1(:,:))

  end if

  if (soil_data_option == 3) then
    CALL READ_SOIL_COMPOSITION(HRLDAS_SETUP_FILE, XSTART, XEND, YSTART, YEND, &
                      NSOIL,IVGTYP,ISICE,ISWATER, &
		      SOILCOMP)
  end if

  if (soil_data_option == 4) then
    CALL READ_3D_SOIL(SPATIAL_FILENAME, XSTART, XEND, YSTART, YEND, &
                      NSOIL,BEXP_3D,SMCDRY_3D,SMCWLT_3D,SMCREF_3D,SMCMAX_3D,  &
		      DKSAT_3D,DWSAT_3D,PSISAT_3D,QUARTZ_3D,REFDK_2D,REFKDT_2D,&
          IRR_FRAC_2D,IRR_HAR_2D,IRR_LAI_2D,IRR_MAD_2D,FILOSS_2D,SPRIR_RATE_2D,&
          MICIR_RATE_2D,FIRTFAC_2D,IR_RAIN_2D,BVIC_2D,AXAJ_2D,BXAJ_2D,XXAJ_2D,&
          BDVIC_2D,GDVIC_2D,BBVIC_2D,&
          KLAT_FAC,TDSMC_FAC,TD_DC,TD_DCOEF,TD_DDRAIN,TD_RADI,TD_SPAC)
  end if

!----------------------------------------------------------------------
! For spatially-varying irrigation parameters, read in necessary extra fields
! only if IOPT_IRR >= 1
!----------------------------------------------------------------------
  if (irrigation_option >= 1) then
      CALL READ_AGRICULTURE_DATA(AGDATA_FLNM, XSTART, XEND, YSTART, YEND, &
                                 IRFRACT, SIFRACT, MIFRACT, FIFRACT)
  end if

!------------------------------------------------------------------------
! For IOPT_RUN = 5 (MMF groundwater), read in necessary extra fields
! This option is not tested for parallel use in the offline driver
!------------------------------------------------------------------------

  if (runoff_option == 5) then
    CALL READ_MMF_RUNOFF(HRLDAS_SETUP_FILE, XSTART, XEND, YSTART, YEND,&
                         FDEPTHXY,EQZWT,RECHCLIM,RIVERBEDXY)
  end if

!------------------------------------------------------------------------
! For OPT_CROP=1 (crop model), read in necessary extra fields
!------------------------------------------------------------------------

  CROPTYPE   = 0       ! make default 0% crops everywhere
  PLANTING   = 126     ! default planting date
  HARVEST    = 290     ! default harvest date
  SEASON_GDD = 1605    ! default total seasonal growing degree days
  if (crop_option == 1) then ! use crop_option=1 for reading crop extra fields
    CALL READ_CROP_INPUT(HRLDAS_SETUP_FILE, XSTART, XEND, YSTART, YEND,&
                         CROPTYPE,PLANTING,HARVEST,SEASON_GDD)
  end if

!------------------------------------------------------------------------
! For IOPT_TDRN = 1 or 2 READ TILE DRAIN MAP
!------------------------------------------------------------------------
  TD_FRACTION = 0.
  IF (IOPT_TDRN .GT. 0) THEN
     CALL READ_TILE_DRAIN_MAP(TDINPUT_FLNM,XSTART,XEND,YSTART,YEND,TD_FRACTION)
  ENDIF

!----------------------------------------------------------------------
! Initialize Model State
!----------------------------------------------------------------------

  SLOPETYP = 2
  DZS       =  SOIL_THICK_INPUT(1:NSOIL)

  ITIMESTEP = 1

  if (restart_filename_requested /= " ") then
     restart_flag = .TRUE.

     call find_restart_file(rank, trim(restart_filename_requested), startdate, khour, olddate, restart_flnm)

     call read_restart(trim(restart_flnm), xstart, xend, xstart, ixfull, jxfull, nsoil, olddate)

#ifdef MPP_LAND
     call mpp_land_bcast_char(19,OLDDATE(1:19))
#endif

     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SOIL_T"  , TSLB     )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SNOW_T"  , TSNOXY   )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SMC"     , SMOIS    )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SH2O"    , SH2O     )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ZSNSO"   , ZSNSOXY  )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SNICE"   , SNICEXY  )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SNLIQ"   , SNLIQXY  )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "QSNOW"   , QSNOWXY  )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "QRAIN"   , QRAINXY  )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "FWET"    , FWETXY   )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SNEQVO"  , SNEQVOXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "EAH"     , EAHXY    )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "TAH"     , TAHXY    )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ALBOLD"  , ALBOLDXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "CM"      , CMXY     )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "CH"      , CHXY     )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ISNOW"   , ISNOWXY  )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "CANLIQ"  , CANLIQXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "CANICE"  , CANICEXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SNEQV"   , SNOW     )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SNOWH"   , SNOWH    )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "TV"      , TVXY     )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "TG"      , TGXY     )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ZWT"     , ZWTXY    )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "WA"      , WAXY     )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "WT"      , WTXY     )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "WSLAKE"  , WSLAKEXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "LFMASS"  , LFMASSXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "RTMASS"  , RTMASSXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "STMASS"  , STMASSXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "WOOD"    , WOODXY   )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "GRAIN"   , GRAINXY  )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "GDD"     , GDDXY    )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "STBLCP"  , STBLCPXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "FASTCP"  , FASTCPXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "LAI"     , LAI      )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SAI"     , XSAIXY   )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "VEGFRA"  , VEGFRA   )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "GVFMIN"  , GVFMIN   )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "GVFMAX"  , GVFMAX   )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACMELT"  , ACSNOM   )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACSNOW"  , ACSNOW   )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "TAUSS"   , TAUSSXY  )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "QSFC"    , QSFC     )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SFCRUNOFF",SFCRUNOFF)
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "UDRUNOFF" ,UDRUNOFF )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "QTDRAIN"  ,QTDRAIN  )
! additional inout variable
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACC_SSOIL" , ACC_SSOILXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACC_QINSUR", ACC_QINSURXY)
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACC_QSEVA" , ACC_QSEVAXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACC_ETRANI", ACC_ETRANIXY)
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACC_DWATER", ACC_DWATERXY)
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACC_PRCP"  , ACC_PRCPXY  )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACC_ECAN"  , ACC_ECANXY  )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACC_ETRAN" , ACC_ETRANXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACC_EDIR"  , ACC_EDIRXY  )

     ! below for irrigation scheme
     IF (irrigation_option > 0) THEN
        call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "IRNUMSI" , IRNUMSI  )
        call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "IRNUMMI" , IRNUMMI  )
        call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "IRNUMFI" , IRNUMFI  )
        call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "IRWATSI" , IRWATSI  )
        call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "IRWATMI" , IRWATMI  )
        call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "IRWATFI" , IRWATFI  )
        call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "IRSIVOL" , IRSIVOL  )
        call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "IRMIVOL" , IRMIVOL  )
        call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "IRFIVOL" , IRFIVOL  )
     ENDIF
    ! below for opt_run = 5
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SMOISEQ"   , SMOISEQ    )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "AREAXY"    , AREAXY     )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SMCWTDXY"  , SMCWTDXY   )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "QRFXY"     , QRFXY      )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "DEEPRECHXY", DEEPRECHXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "QSPRINGXY" , QSPRINGXY  )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "QSLATXY"   , QSLATXY    )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "QRFSXY"    , QRFSXY     )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "QSPRINGSXY", QSPRINGSXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "RECHXY"    , RECHXY     )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "FDEPTHXY"   ,FDEPTHXY   )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "RIVERCONDXY",RIVERCONDXY)
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "RIVERBEDXY" ,RIVERBEDXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "EQZWT"      ,EQZWT      )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "PEXPXY"     ,PEXPXY     )

     IF(SF_URBAN_PHYSICS > 0 )  THEN  ! any urban model

       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,      "SH_URB2D"  ,     SH_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,      "LH_URB2D"  ,     LH_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,       "G_URB2D"  ,      G_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,      "RN_URB2D"  ,     RN_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,      "TS_URB2D"  ,     TS_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "FRC_URB2D"  ,    FRC_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,   "UTYPE_URB2D"  ,  UTYPE_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,      "LP_URB2D"  ,     LP_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,      "LB_URB2D"  ,     LB_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "HGT_URB2D"  ,    HGT_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,      "MH_URB2D"  ,     MH_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,    "STDH_URB2D"  ,   STDH_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,      "HI_URB2D"  ,     HI_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,      "LF_URB2D"  ,     LF_URB2D  )

     ENDIF

     IF(SF_URBAN_PHYSICS == 1 ) THEN  ! single layer urban model

       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "CMR_SFCDIF" ,    CMR_SFCDIF )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "CHR_SFCDIF" ,    CHR_SFCDIF )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "CMC_SFCDIF" ,    CMC_SFCDIF )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "CHC_SFCDIF" ,    CHC_SFCDIF )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,    "CMGR_SFCDIF" ,   CMGR_SFCDIF )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,    "CHGR_SFCDIF" ,   CHGR_SFCDIF )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,      "TR_URB2D"  ,     TR_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,      "TB_URB2D"  ,     TB_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,      "TG_URB2D"  ,     TG_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,      "TC_URB2D"  ,     TC_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,      "QC_URB2D"  ,     QC_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,      "UC_URB2D"  ,     UC_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,    "XXXR_URB2D"  ,   XXXR_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,    "XXXB_URB2D"  ,   XXXB_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,    "XXXG_URB2D"  ,   XXXG_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,    "XXXC_URB2D"  ,   XXXC_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "TRL_URB3D"  ,    TRL_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "TBL_URB3D"  ,    TBL_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "TGL_URB3D"  ,    TGL_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,    "CMCR_URB2D"  ,   CMCR_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "TGR_URB2D"  ,    TGR_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,    "TGRL_URB3D"  ,   TGRL_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "SMR_URB3D"  ,    SMR_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,   "DRELR_URB2D"  ,  DRELR_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,   "DRELB_URB2D"  ,  DRELB_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,   "DRELG_URB2D"  ,  DRELG_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "FLXHUMR_URB2D"  ,FLXHUMR_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "FLXHUMB_URB2D"  ,FLXHUMB_URB2D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "FLXHUMG_URB2D"  ,FLXHUMG_URB2D  )

     ENDIF

     IF(SF_URBAN_PHYSICS == 2 .or. SF_URBAN_PHYSICS == 3) THEN  ! BEP or BEM urban models

       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "TRB_URB4D"  ,    TRB_URB4D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "TW1_URB4D"  ,    TW1_URB4D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "TW2_URB4D"  ,    TW2_URB4D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "TGB_URB4D"  ,    TGB_URB4D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,    "SFW1_URB3D"  ,   SFW1_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,    "SFW2_URB3D"  ,   SFW2_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "SFR_URB3D"  ,    SFR_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "SFG_URB3D"  ,    SFG_URB3D  )

     ENDIF

     IF(SF_URBAN_PHYSICS == 3) THEN  ! BEM urban model

       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,    "TLEV_URB3D"  ,   TLEV_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,    "QLEV_URB3D"  ,   QLEV_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "TW1LEV_URB3D"  , TW1LEV_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "TW2LEV_URB3D"  , TW2LEV_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,   "TGLEV_URB3D"  ,  TGLEV_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,   "TFLEV_URB3D"  ,  TFLEV_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,   "SF_AC_URB3D"  ,  SF_AC_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,   "LF_AC_URB3D"  ,  LF_AC_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,   "CM_AC_URB3D"  ,  CM_AC_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "SFVENT_URB3D"  , SFVENT_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "LFVENT_URB3D"  , LFVENT_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "SFWIN1_URB3D"  , SFWIN1_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "SFWIN2_URB3D"  , SFWIN2_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,   "EP_PV_URB3D"  ,  EP_PV_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,    "T_PV_URB3D"  ,   T_PV_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "TRV_URB4D"  ,    TRV_URB4D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,      "QR_URB4D"  ,     QR_URB4D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "QGR_URB3D"  ,    QGR_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "TGR_URB3D"  ,    TGR_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,   "DRAIN_URB4D"  ,  DRAIN_URB4D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "DRAINGR_URB3D"  ,DRAINGR_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,    "SFRV_URB3D"  ,   SFRV_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,    "LFRV_URB3D"  ,   LFRV_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "DGR_URB3D"  ,    DGR_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,      "DG_URB3D"  ,     DG_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "LFR_URB3D"  ,    LFR_URB3D  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "LFG_URB3D"  ,    LFG_URB3D  )

     ENDIF

     IF (IOPT_MOSAIC.eq.1) THEN
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "TSLB_mosaic"  ,     TSLB_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,    "SMOIS_mosaic"  ,    SMOIS_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "SH2O_mosaic"  ,     SH2O_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,   "tsnoxy_mosaic"  ,   tsnoxy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "zsnsoxy_mosaic"  ,  zsnsoxy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "snicexy_mosaic"  ,  snicexy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "snliqxy_mosaic"  ,  snliqxy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "qsnowxy_mosaic"  ,  qsnowxy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "qrainxy_mosaic"  ,  qrainxy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,   "fwetxy_mosaic"  ,   fwetxy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "sneqvoxy_mosaic" , sneqvoxy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,    "eahxy_mosaic"  ,    eahxy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,    "tahxy_mosaic"  ,    tahxy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "alboldxy_mosaic"  , alboldxy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "cmxy_mosaic"  ,     cmxy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "chxy_mosaic"  ,     chxy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "isnowxy_mosaic"  ,  isnowxy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "canliqxy_mosaic" , canliqxy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "canicexy_mosaic" , canicexy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "SNOW_mosaic"  ,     SNOW_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,    "SNOWH_mosaic"  ,    SNOWH_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "tvxy_mosaic"  ,     tvxy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "tgxy_mosaic"  ,     tgxy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,    "zwtxy_mosaic"  ,    zwtxy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "waxy_mosaic"  ,     waxy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "wtxy_mosaic"  ,     wtxy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "lfmassxy_mosaic"  , lfmassxy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "rtmassxy_mosaic"  , rtmassxy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "stmassxy_mosaic"  , stmassxy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,   "woodxy_mosaic"  ,   woodxy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "grainxy_mosaic"  ,  grainxy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,    "gddxy_mosaic"  ,    gddxy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "stblcpxy_mosaic"  , stblcpxy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "fastcpxy_mosaic"  , fastcpxy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,   "xsaixy_mosaic"  ,   xsaixy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "xlai_mosaic"  ,     xlai_mosaic  )

       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,    "VEGFRA_mosaic" ,  VEGFRA_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,   "ACSNOM_mosaic"  ,   ACSNOM_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,   "ACSNOW_mosaic"  ,   ACSNOW_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,   "taussxy_mosaic" ,  taussxy_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,   "QSFC_mosaic"    ,     QSFC_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "SFCRUNOFF_mosaic",SFCRUNOFF_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "UDRUNOFF_mosaic" , UDRUNOFF_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACC_SSOILXY_mosaic" , ACC_SSOILXY_mosaic )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACC_QINSURXY_mosaic", ACC_QINSURXY_mosaic)
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACC_QSEVAXY_mosaic" , ACC_QSEVAXY_mosaic )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACC_ETRANIXY_mosaic", ACC_ETRANIXY_mosaic)
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACC_DWATERXY_mosaic", ACC_DWATERXY_mosaic)
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACC_PRCPXY_mosaic"  , ACC_PRCPXY_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACC_ECANXY_mosaic"  , ACC_ECANXY_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACC_ETRANXY_mosaic" , ACC_ETRANXY_mosaic )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACC_EDIRXY_mosaic"  , ACC_EDIRXY_mosaic  )

       
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,      "TSK_mosaic"  ,      TSK_mosaic  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull,   "CANWAT_mosaic"  ,   CANWAT_mosaic  )
       
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "smoiseq_mosaic"   , SMOISEQ_mosaic    )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "smcwtdxy_mosaic"  , smcwtdxy_mosaic   )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "deeprechxy_mosaic", deeprechxy_mosaic )



       !! Irrigation varaibles for re-writing
       IF (irrigation_option > 0) THEN
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "IRNUMSI_mosaic"  ,  IRNUMSI_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "IRNUMMI_mosaic"  ,  IRNUMMI_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "IRNUMFI_mosaic"  ,  IRNUMFI_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "IRWATSI_mosaic"  ,  IRWATSI_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "IRWATMI_mosaic"  ,  IRWATMI_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "IRWATFI_mosaic"  ,  IRWATFI_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "IRELOSS_mosaic"  ,  IRELOSS_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "IRSIVOL_mosaic"  ,  IRSIVOL_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "IRMIVOL_mosaic"  ,  IRMIVOL_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "IRFIVOL_mosaic"  ,  IRFIVOL_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,  "IRRSPLH_mosaic"  ,  IRRSPLH_mosaic  )
       ENDIF


       IF(SF_URBAN_PHYSICS > 0 )  THEN  ! any urban model

         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,      "SH_URB2D_mosaic"  ,     SH_URB2D_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,      "LH_URB2D_mosaic"  ,     LH_URB2D_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,       "G_URB2D_mosaic"  ,      G_URB2D_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,      "RN_URB2D_mosaic"  ,     RN_URB2D_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,      "TS_URB2D_mosaic"  ,     TS_URB2D_mosaic  )
       ENDIF

       IF(SF_URBAN_PHYSICS == 1 ) THEN  ! single layer urban model

         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,      "TR_URB2D_mosaic"  ,     TR_URB2D  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,      "TB_URB2D_mosaic"  ,     TB_URB2D  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,      "TG_URB2D_mosaic"  ,     TG_URB2D  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,      "TC_URB2D_mosaic"  ,     TC_URB2D  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,      "QC_URB2D_mosaic"  ,     QC_URB2D  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,      "UC_URB2D_mosaic"  ,     UC_URB2D  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "TRL_URB3D_mosaic"  ,    TRL_URB3D  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "TBL_URB3D_mosaic"  ,    TBL_URB3D  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "TGL_URB3D_mosaic"  ,    TGL_URB3D  )

         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "CMR_SFCDIF_mosaic"  ,    CMR_SFCDIF_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "CHR_SFCDIF_mosaic"  ,    CHR_SFCDIF_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "CMC_SFCDIF_mosaic"  ,    CMC_SFCDIF_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "CHC_SFCDIF_mosaic"  ,    CHC_SFCDIF_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "CMGR_SFCDIF_mosaic"  ,    CMGR_SFCDIF_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "CHGR_SFCDIF_mosaic"  ,    CHGR_SFCDIF_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "XXXR_URB2D_mosaic"  ,    XXXR_URB2D_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "XXXB_URB2D_mosaic"  ,    XXXB_URB2D_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "XXXG_URB2D_mosaic"  ,    XXXG_URB2D_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "XXXC_URB2D_mosaic"  ,    XXXC_URB2D_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "CMCR_URB2D_mosaic"  ,    CMCR_URB2D_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "TGR_URB2D_mosaic"   ,    TGR_URB2D_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "TGRL_URB3D_mosaic"  ,    TGRL_URB3D_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "SMR_URB3D_mosaic"  ,    SMR_URB3D_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "DRELR_URB2D_mosaic"  ,    DRELR_URB2D_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "DRELB_URB2D_mosaic"  ,    DRELB_URB2D_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "DRELG_URB2D_mosaic"  ,    DRELG_URB2D_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "FLXHUMR_URB2D_mosaic"  ,    FLXHUMR_URB2D_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "FLXHUMB_URB2D_mosaic"  ,    FLXHUMB_URB2D_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "FLXHUMG_URB2D_mosaic"  ,    FLXHUMG_URB2D_mosaic  )

       ENDIF

       !Hue variablesß
       IF(IOPT_HUE.eq.1) THEN
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "DETENTION_STORAGEXY_mosaic"  ,    DETENTION_STORAGEXY_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "RUNONSFXY_mosaic"  ,    RUNONSFXY_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "VOL_FLUX_SMXY_mosaic"  ,    VOL_FLUX_SMXY_mosaic  )
         call get_from_restart(xstart, xend, xstart, ixfull, jxfull,     "VOL_FLUX_RUNONXY_mosaic"  ,    VOL_FLUX_RUNONXY_mosaic  )
       ENDIF

     ENDIF
     STEPWTD = nint(WTDDT*60./DTBL)
     STEPWTD = max(STEPWTD,1)

! Must still call NOAHMP_INIT even in restart to set up parameter arrays (also done in WRF)

     CALL NOAHMP_INIT(    LLANDUSE,     SNOW,    SNOWH,   CANWAT,   ISLTYP,   IVGTYP, XLAT, &   ! call from WRF phys_init
                    TSLB,    SMOIS,     SH2O,      DZS, FNDSOILW, FNDSNOWH, &
                     TSK,  ISNOWXY,     TVXY,     TGXY, CANICEXY,      TMN,     XICE, &
                CANLIQXY,    EAHXY,    TAHXY,     CMXY,     CHXY,                     &
                  FWETXY, SNEQVOXY, ALBOLDXY,  QSNOWXY, QRAINXY,  WSLAKEXY,    ZWTXY,     WAXY, &
                    WTXY,   TSNOXY,  ZSNSOXY,  SNICEXY,  SNLIQXY, LFMASSXY, RTMASSXY, &
                STMASSXY,   WOODXY, STBLCPXY, FASTCPXY,   XSAIXY, LAI,                    &
                 GRAINXY,    GDDXY,                                                   &
                CROPTYPE,  CROPCAT,                                                   &
                irnumsi  ,irnummi  ,irnumfi  ,irwatsi,    &
                irwatmi  ,irwatfi  ,ireloss  ,irsivol,    &
                irmivol  ,irfivol  ,irrsplh  ,            &
                  T2MVXY,   T2MBXY, CHSTARXY,                                         &
                   NSOIL,  .true.,                                                   &
                  .true.,runoff_option, crop_option, irrigation_option, irrigation_method,   &
                  sf_urban_physics,                         &  ! urban scheme
                  ids,ide+1, jds,jde+1, kds,kde,                &  ! domain
                  ims,ime, jms,jme, kms,kme,                &  ! memory
                  its,ite, jts,jte, kts,kte                 &  ! tile
                     ,smoiseq  ,smcwtdxy ,rechxy, deeprechxy, qtdrain, areaxy ,dx, dy, msftx, msfty,&
                     wtddt    ,stepwtd  ,dtbl  ,qrfsxy ,qspringsxy  ,qslatxy,                  &
                     fdepthxy ,terrain       ,riverbedxy ,eqzwt ,rivercondxy ,pexpxy,              &
                     rechclim ,gecros_state                 &
                     )

         IF(SF_URBAN_PHYSICS > 0 ) THEN  !urban

                CALL urban_param_init(DZR,DZB,DZG,num_soil_layers,sf_urban_physics,use_wudapt_lcz)

                CALL urban_var_init(ISURBAN,TSK,TSLB,TMN,IVGTYP,                        & !urban
                      ims,ime,jms,jme,kms,kme,num_soil_layers,                          & !urban
                      LCZ_1_TABLE,LCZ_2_TABLE,LCZ_3_TABLE,LCZ_4_TABLE,LCZ_5_TABLE,      & !urban
		      LCZ_6_TABLE,LCZ_7_TABLE,LCZ_8_TABLE,LCZ_9_TABLE,LCZ_10_TABLE,     & !urban
		      LCZ_11_TABLE,                                                     & !urban
                      .true.,sf_urban_physics,                                          & !urban
                      XXXR_URB2D,    XXXB_URB2D,    XXXG_URB2D,XXXC_URB2D,              & !urban
                        TR_URB2D,      TB_URB2D,      TG_URB2D,  TC_URB2D, QC_URB2D,    & !urban
                       TRL_URB3D,     TBL_URB3D,     TGL_URB3D,                         & !urban
                        SH_URB2D,      LH_URB2D,       G_URB2D,  RN_URB2D, TS_URB2D,    & !urban
                   num_urban_ndm,  urban_map_zrd,  urban_map_zwd, urban_map_gd,         & !I multi-layer urban
                    urban_map_zd,  urban_map_zdf,   urban_map_bd, urban_map_wd,         & !I multi-layer urban
                   urban_map_gbd,  urban_map_fbd, urban_map_zgrd,                       & !I multi-layer urban
                    num_urban_hi,                                                       & !urban
                       TRB_URB4D,     TW1_URB4D,     TW2_URB4D, TGB_URB4D,              & !urban
                      TLEV_URB3D,    QLEV_URB3D,                                        & !urban
                    TW1LEV_URB3D,  TW2LEV_URB3D,                                        & !urban
                     TGLEV_URB3D,   TFLEV_URB3D,                                        & !urban
                     SF_AC_URB3D,   LF_AC_URB3D,   CM_AC_URB3D,                         & !urban
                    SFVENT_URB3D,  LFVENT_URB3D,                                        & !urban
                    SFWIN1_URB3D,  SFWIN2_URB3D,                                        & !urban
                      SFW1_URB3D,    SFW2_URB3D,     SFR_URB3D, SFG_URB3D,              & !urban
                     EP_PV_URB3D,    T_PV_URB3D,                                        & !GRZ
                       TRV_URB4D,      QR_URB4D,     QGR_URB3D, TGR_URB3D,              & !GRZ
                     DRAIN_URB4D, DRAINGR_URB3D,    SFRV_URB3D,                         & !GRZ
                      LFRV_URB3D,     DGR_URB3D,      DG_URB3D, LFR_URB3D, LFG_URB3D,   & !GRZ
                           SMOIS,                                                       & !GRZ
                        LP_URB2D,      HI_URB2D,      LB_URB2D,                         & !urban
                       HGT_URB2D,      MH_URB2D,    STDH_URB2D,                         & !urban
                        LF_URB2D,                                                       & !urban
                      CMCR_URB2D,     TGR_URB2D,    TGRL_URB3D, SMR_URB3D,              & !urban
                     DRELR_URB2D,   DRELB_URB2D,   DRELG_URB2D,                         & !urban
                   FLXHUMR_URB2D, FLXHUMB_URB2D, FLXHUMG_URB2D,                         & !urban
                         A_U_BEP,       A_V_BEP,       A_T_BEP,   A_Q_BEP,              & !multi-layer urban
                         A_E_BEP,       B_U_BEP,       B_V_BEP,                         & !multi-layer urban
                         B_T_BEP,       B_Q_BEP,       B_E_BEP,   DLG_BEP,              & !multi-layer urban
                        DL_U_BEP,        SF_BEP,        VL_BEP,                         & !multi-layer urban
                       FRC_URB2D,   UTYPE_URB2D, use_wudapt_lcz)                          !urban

               max_utype_urb2d = maxval(UTYPE_URB2D)*1.0
               IF (use_wudapt_lcz.eq.0 .and. max_utype_urb2d.gt.3.0) THEN  !new LCZ
                 CALL wrf_error_fatal &
                 ('USING 10 WUDAPT LCZ WITHOUT URBPARM_LCZ.TBL. SET USE_WUDAPT_LCZ=1')
               ENDIF
               IF (use_wudapt_lcz.eq.1 .and. max_utype_urb2d.le.3.0) THEN  ! new LCZ
                 CALL wrf_error_fatal &
                 ('USING URBPARM_LCZ.TBL WITH OLD 3 URBAN CLASSES. SET USE_WUDAPT_LCZ=0')
               ENDIF


          ENDIF

          IF (IOPT_MOSAIC.eq.1) then
          CALL NOAHMP_MOSAIC_INIT(XLAND, SNOW,    SNOWH,   CANWAT,   ISLTYP,   IVGTYP, XLAT,&   ! call from WRF phys_init
                        TSLB,    SMOIS,     SH2O,      DZS, FNDSOILW, FNDSNOWH,  &
                         TSK,  ISNOWXY,     TVXY,     TGXY, CANICEXY,      TMN,     XICE, &
                    CANLIQXY,    EAHXY,    TAHXY,     CMXY,     CHXY,                     &
                      FWETXY, SNEQVOXY, ALBOLDXY,  QSNOWXY, QRAINXY,  WSLAKEXY,    ZWTXY,     WAXY, &
                        WTXY,   TSNOXY,  ZSNSOXY,  SNICEXY,  SNLIQXY, LFMASSXY, RTMASSXY, &
                    STMASSXY,   WOODXY, STBLCPXY, FASTCPXY,   XSAIXY, LAI,                    &
                     GRAINXY,    GDDXY,                                                   &
                     CROPTYPE, CROPCAT,                                                   &
                     irnumsi, irnummi, irnumfi, irwatsi,                                  &
                     irwatmi, irwatfi, ireloss, irsivol,                                  &
                     irmivol, irfivol, irrsplh,                                           &
                      T2MVXY,   T2MBXY, CHSTARXY,                                         &
                       NSOIL,  .true.,                                                   &
                      .true.,runoff_option, crop_option, irrigation_option, irrigation_method,  &
                      sf_urban_physics, ISWATER, ISICE,                        &  
                      ISURBAN, 0, &                              ! urban scheme
                      ids,ide+1, jds,jde+1, kds,kde,                &  ! domain
                      ims,ime, jms,jme, kms,kme,                &  ! memory
                      its,ite, jts,jte, kts,kte,                 &  ! tile
                      smoiseq, smcwtdxy, rechxy, deeprechxy,     &
                      LANDUSEF, LANDUSEF2, number_land_use_catagories, IOPT_MOSAIC,                     &      ! Added by Aaron A. **IMPORTANT FOR THE FOR LOOPS OF FRACTIONAL 
                      mosaic_cat_index, number_mosaic_catagories,                                         &      ! Added by Aaron A. **IMPORTANT FOR THE FOR LOOPS OF FRACTIONAL LAND USE
                      IOPT_HUE,                                                      &
                      TSK_mosaic, TSLB_mosaic, SMOIS_mosaic, SH2O_mosaic,                   &      ! Added by Aaron A.
                      CANWAT_mosaic, SNOW_mosaic, SNOWH_mosaic,                             &      ! Added by Aaron A.
                      isnowxy_mosaic, tvxy_mosaic, tgxy_mosaic, canicexy_mosaic,            &      ! Added by Aaron A.
                      canliqxy_mosaic, eahxy_mosaic, tahxy_mosaic,              &      ! Added by Aaron A.
                      cmxy_mosaic, chxy_mosaic, fwetxy_mosaic, sneqvoxy_mosaic,             &      ! Added by Aaron A.
                      alboldxy_mosaic, qsnowxy_mosaic, qrainxy_mosaic, wslakexy_mosaic, zwtxy_mosaic,       &      ! Added by Aaron A.
                      waxy_mosaic, wtxy_mosaic, tsnoxy_mosaic, zsnsoxy_mosaic,              &      ! Added by Aaron A.
                      snicexy_mosaic, snliqxy_mosaic, lfmassxy_mosaic, rtmassxy_mosaic,     &      ! Added by Aaron A.
                      stmassxy_mosaic, woodxy_mosaic, stblcpxy_mosaic, fastcpxy_mosaic,     &      ! Added by Aaron A.
                      xsaixy_mosaic, xlai_mosaic, &

                      IRNUMSI_mosaic, IRNUMMI_mosaic, IRNUMFI_mosaic, IRWATSI_mosaic,       &
                      IRWATMI_mosaic, IRWATFI_mosaic, IRELOSS_mosaic, IRSIVOL_mosaic,       &
                      IRMIVOL_mosaic, IRFIVOL_mosaic, IRRSPLH_mosaic,                       &


                      smoiseq_mosaic,  smcwtdxy_mosaic, rechxy_mosaic, deeprechxy_mosaic,   &

                      TR_URB2D_mosaic, TB_URB2D_mosaic,                                     &      ! Added by Aaron A.
                      TG_URB2D_mosaic, TC_URB2D_mosaic, QC_URB2D_mosaic,                    &      ! Added by Aaron A.
                      TRL_URB3D_mosaic, TBL_URB3D_mosaic,                                   &      ! Added by Aaron A.
                      TGL_URB3D_mosaic, SH_URB2D_mosaic, LH_URB2D_mosaic,                   &      ! Added by Aaron A.
                      G_URB2D_mosaic, RN_URB2D_mosaic, TS_URB2D_mosaic,                     &      ! Added by Aaron A.
                      CMR_SFCDIF_mosaic, CHR_SFCDIF_mosaic, CMC_SFCDIF_mosaic,              &
                      CHC_SFCDIF_mosaic, CMGR_SFCDIF_mosaic, CHGR_SFCDIF_mosaic,            &
                      XXXR_URB2D_mosaic, XXXB_URB2D_mosaic, XXXG_URB2D_mosaic,              &
                      XXXC_URB2D_mosaic,                                                    &
                      CMCR_URB2D_mosaic, TGR_URB2D_mosaic,                                  &
                      TGRL_URB3D_mosaic, SMR_URB3D_mosaic,                                  &
                      DRELR_URB2D_mosaic, DRELB_URB2D_mosaic, DRELG_URB2D_mosaic,           &
                      FLXHUMR_URB2D_mosaic, FLXHUMB_URB2D_mosaic, FLXHUMG_URB2D_mosaic,      &
                      DETENTION_STORAGEXY_mosaic, &
                      Z0, ZNT_mosaic, Z0_mosaic)                                                            ! Added by Aaron A.

       ENDIF !end mosaic if else statement

  else
!----------------------------------------------------------------------
! First branch of the if else statement of the noah-mp mosaic scheme
! If there is no mosaic option happening, we just continue as normal
! Added by Aaron A. 19 May 2022
!-----------------------------------------------------------------------
if (IOPT_MOSAIC.eq.0) then

     restart_flag = .FALSE.

     SMOIS     =  undefined_real
     TSLB      =  undefined_real
     SH2O      =  undefined_real
     CANLIQXY  =  undefined_real
     TSK       =  undefined_real
     RAINBL_tmp    =  undefined_real
     SNOW      =  undefined_real
     SNOWH     =  undefined_real

#ifdef WRF_HYDRO
  if((forc_typ .gt. 2) .and. (forc_typ .ne. 6) ) HRLDAS_ini_typ = 0

  if(HRLDAS_ini_typ .eq. 1) then
     ! read initial parameters and conditions from the HRLDAS forcing data
     if(forc_typ .eq. 2) then
          inflnm = trim(indir)//"/"//&
          startdate(1:4)//startdate(6:7)//startdate(9:10)//startdate(12:13)//&
          startdate(15:16)//".LDASIN_DOMAIN"//hgrid
     else
          inflnm = trim(indir)//"/"//&
          startdate(1:4)//startdate(6:7)//startdate(9:10)//startdate(12:13)//&
          ".LDASIN_DOMAIN"//hgrid
     endif

#else

     inflnm = trim(indir)//"/"//&
          startdate(1:4)//startdate(6:7)//startdate(9:10)//startdate(12:13)//&
          ".LDASIN_DOMAIN"//hgrid
#endif

     CALL READINIT_HRLDAS(HRLDAS_SETUP_FILE, xstart, xend, ystart, yend,  &
          NSOIL, DZS, OLDDATE, LDASIN_VERSION, SMOIS,       &
          TSLB, CANWAT, TSK, SNOW, SNOWH, FNDSNOWH)

     VEGFRA    =  undefined_real
     LAI       =  undefined_real
     GVFMIN    =  undefined_real
     GVFMAX    =  undefined_real

     CALL READVEG_HRLDAS(HRLDAS_SETUP_FILE, xstart, xend, ystart, yend,  &
          OLDDATE, IVGTYP, VEGFRA, LAI, GVFMIN, GVFMAX)

#ifdef WRF_HYDRO
  else   !  HRLDAS_ini_typ

#ifdef MPP_LAND
    call HYDRO_HRLDAS_ini_mpp   &
#else
    call HYDRO_HRLDAS_ini   &
#endif
       (trim(hrldas_setup_file), xend-xstart+1,yend-ystart+1, &
       nsoil,SMOIS(:,1:NSOIL,:),TSLB(:,1:NSOIL,:),SH2O(:,1:NSOIL,:), CANWAT, TSK,SNOW,SNOWH,lai,VEGFRA,IVGTYP,FNDSNOWH)

       greenfrac = 0.0
#ifdef MPP_LAND
       call get_greenfrac_mpp &
#else
       call get_greenfrac &
#endif
        (trim(GEO_STATIC_FLNM),greenfrac, ix, jx, olddate, SHDMAX)
       !yw GVFMAX = maxval(greenfrac)
  endif
#endif

!     SNOW = SNOW * 1000.    ! Convert snow water equivalent to mm. MB: remove v3.7

     FNDSOILW = .FALSE.
     CALL NOAHMP_INIT(    LLANDUSE,     SNOW,    SNOWH,   CANWAT,   ISLTYP,   IVGTYP, XLAT, &   ! call from WRF phys_init
                    TSLB,    SMOIS,     SH2O,      DZS, FNDSOILW, FNDSNOWH, &
                     TSK,  ISNOWXY,     TVXY,     TGXY, CANICEXY,      TMN,     XICE, &
                CANLIQXY,    EAHXY,    TAHXY,     CMXY,     CHXY,                     &
                  FWETXY, SNEQVOXY, ALBOLDXY,  QSNOWXY, QRAINXY,  WSLAKEXY,    ZWTXY,     WAXY, &
                    WTXY,   TSNOXY,  ZSNSOXY,  SNICEXY,  SNLIQXY, LFMASSXY, RTMASSXY, &
                STMASSXY,   WOODXY, STBLCPXY, FASTCPXY,   XSAIXY, LAI,                    &
                 GRAINXY,    GDDXY,                                                   &
                CROPTYPE,  CROPCAT,                                                   &
                irnumsi  ,irnummi  ,irnumfi  ,irwatsi,    &
                irwatmi  ,irwatfi  ,ireloss  ,irsivol,    &
                irmivol  ,irfivol  ,irrsplh  ,            &
                  T2MVXY,   T2MBXY, CHSTARXY,                                         &
                   NSOIL,  .false.,                                                   &
                  .true.,runoff_option, crop_option, irrigation_option, irrigation_method,  &
                  sf_urban_physics,                         &  ! urban scheme
                  ids,ide+1, jds,jde+1, kds,kde,                &  ! domain
                  ims,ime, jms,jme, kms,kme,                &  ! memory
                  its,ite, jts,jte, kts,kte                 &  ! tile
                     ,smoiseq  ,smcwtdxy ,rechxy, deeprechxy, qtdrain, areaxy ,dx, dy, msftx, msfty,&
                     wtddt    ,stepwtd  ,dtbl  ,qrfsxy ,qspringsxy  ,qslatxy,                  &
                     fdepthxy ,terrain       ,riverbedxy ,eqzwt ,rivercondxy ,pexpxy,              &
                     rechclim ,gecros_state                 &
                     )

         if(iopt_run == 5) then
              call groundwater_init ( grid,                                     &
                 num_soil_layers, dzs, isltyp, ivgtyp, wtddt ,                  &
                 fdepthxy   , terrain , riverbedxy, eqzwt     ,                 &
                 rivercondxy, pexpxy  , areaxy    , zwtxy     ,                 &
                 smois      , sh2o    , smoiseq   , smcwtdxy  ,                 &
                 QLATXY     , qslatxy , QRFXY     , qrfsxy    ,                 &
                 deeprechxy , rechxy  , QSPRINGXY , qspringsxy,                 &
                 rechclim   ,                                                   &
                 ids,ide, jds,jde, kds,kde,                                     &
                 ims,ime, jms,jme, kms,kme,                                     &
                 ims,ime, jms,jme, kms,kme,                                     &
                 its,ite, jts,jte, kts,kte )
         endif

      TAUSSXY = 0.0   ! Need to be added to _INIT later

         IF(SF_URBAN_PHYSICS > 0 ) THEN  !urban

                 CALL urban_param_init(DZR,DZB,DZG,num_soil_layers,sf_urban_physics,use_wudapt_lcz)

                CALL urban_var_init(ISURBAN,TSK,TSLB,TMN,IVGTYP,                        & !urban
                      ims,ime,jms,jme,kms,kme,num_soil_layers,                          & !urban
                      LCZ_1_TABLE,LCZ_2_TABLE,LCZ_3_TABLE,LCZ_4_TABLE,LCZ_5_TABLE,      & !urban
                      LCZ_6_TABLE,LCZ_7_TABLE,LCZ_8_TABLE,LCZ_9_TABLE,LCZ_10_TABLE,     & !urban
                      LCZ_11_TABLE,                                                     & !urban
                      .false.,sf_urban_physics,                                         & !urban
                      XXXR_URB2D,    XXXB_URB2D,    XXXG_URB2D,XXXC_URB2D,              & !urban
                        TR_URB2D,      TB_URB2D,      TG_URB2D,  TC_URB2D, QC_URB2D,    & !urban
                       TRL_URB3D,     TBL_URB3D,     TGL_URB3D,                         & !urban
                        SH_URB2D,      LH_URB2D,       G_URB2D,  RN_URB2D, TS_URB2D,    & !urban
                   num_urban_ndm,  urban_map_zrd,  urban_map_zwd, urban_map_gd,         & !I multi-layer urban
                    urban_map_zd,  urban_map_zdf,   urban_map_bd, urban_map_wd,         & !I multi-layer urban
                   urban_map_gbd,  urban_map_fbd, urban_map_zgrd,                       & !I multi-layer urban
                    num_urban_hi,                                                       & !urban
                       TRB_URB4D,     TW1_URB4D,     TW2_URB4D, TGB_URB4D,              & !urban
                      TLEV_URB3D,    QLEV_URB3D,                                        & !urban
                    TW1LEV_URB3D,  TW2LEV_URB3D,                                        & !urban
                     TGLEV_URB3D,   TFLEV_URB3D,                                        & !urban
                     SF_AC_URB3D,   LF_AC_URB3D,   CM_AC_URB3D,                         & !urban
                    SFVENT_URB3D,  LFVENT_URB3D,                                        & !urban
                    SFWIN1_URB3D,  SFWIN2_URB3D,                                        & !urban
                      SFW1_URB3D,    SFW2_URB3D,     SFR_URB3D, SFG_URB3D,              & !urban
                     EP_PV_URB3D,    T_PV_URB3D,                                        & !GRZ
                       TRV_URB4D,      QR_URB4D,     QGR_URB3D, TGR_URB3D,              & !GRZ
                     DRAIN_URB4D, DRAINGR_URB3D,    SFRV_URB3D,                         & !GRZ
                      LFRV_URB3D,     DGR_URB3D,      DG_URB3D, LFR_URB3D, LFG_URB3D,   & !GRZ
                           SMOIS,                                                       & !GRZ
                        LP_URB2D,      HI_URB2D,      LB_URB2D,                         & !urban
                       HGT_URB2D,      MH_URB2D,    STDH_URB2D,                         & !urban
                        LF_URB2D,                                                       & !urban
                      CMCR_URB2D,     TGR_URB2D,    TGRL_URB3D, SMR_URB3D,              & !urban
                     DRELR_URB2D,   DRELB_URB2D,   DRELG_URB2D,                         & !urban
                   FLXHUMR_URB2D, FLXHUMB_URB2D, FLXHUMG_URB2D,                         & !urban
                         A_U_BEP,       A_V_BEP,       A_T_BEP,   A_Q_BEP,              & !multi-layer urban
                         A_E_BEP,       B_U_BEP,       B_V_BEP,                         & !multi-layer urban
                         B_T_BEP,       B_Q_BEP,       B_E_BEP,   DLG_BEP,              & !multi-layer urban
                        DL_U_BEP,        SF_BEP,        VL_BEP,                         & !multi-layer urban
                       FRC_URB2D,   UTYPE_URB2D, use_wudapt_lcz)                          !urban

               max_utype_urb2d = maxval(UTYPE_URB2D)*1.0
               IF (use_wudapt_lcz.eq.0 .and. max_utype_urb2d.gt.3.0) THEN  !new LCZ
                 CALL wrf_error_fatal &
                 ('USING 10 WUDAPT LCZ WITHOUT URBPARM_LCZ.TBL. SET USE_WUDAPT_LCZ=1')
               ENDIF
               IF (use_wudapt_lcz.eq.1 .and. max_utype_urb2d.le.3.0) THEN  ! new LCZ
                 CALL wrf_error_fatal &
                 ('USING URBPARM_LCZ.TBL WITH OLD 3 URBAN CLASSES. SET USE_WUDAPT_LCZ=0')
               ENDIF

          ENDIF

          IF(SF_URBAN_PHYSICS > 1 ) THEN  !urban
	    do i = 1, num_urban_atmosphere-1
              dz_urban(:,i,:) = urban_atmosphere_thickness  ! thickness of full levels
	    end do
	    dz_urban(:,num_urban_atmosphere,:) =                                   &
                  2*(zlvl - urban_atmosphere_thickness*(num_urban_atmosphere-1))
	    print*, dz_urban(1,:,1)
	  ENDIF
  !----------------------------------------------------------------------------
  ! HUE NOAH-MP Mosaic Scheme switch. This is the portion of the model
  ! that integrates multiple land-types.
  ! Added by Aaron A. 19 May 2022
  !----------------------------------------------------------------------------
  else if(IOPT_MOSAIC.eq.1) then

    restart_flag = .FALSE.

     SMOIS     =  undefined_real
     TSLB      =  undefined_real
     SH2O      =  undefined_real
     CANLIQXY  =  undefined_real
     TSK       =  undefined_real
     RAINBL_tmp    =  undefined_real
     SNOW      =  undefined_real
     SNOWH     =  undefined_real

#ifdef WRF_HYDRO
  if((forc_typ .gt. 2) .and. (forc_typ .ne. 6) ) HRLDAS_ini_typ = 0

  if(HRLDAS_ini_typ .eq. 1) then
     ! read initial parameters and conditions from the HRLDAS forcing data
     if(forc_typ .eq. 2) then
          inflnm = trim(indir)//"/"//&
          startdate(1:4)//startdate(6:7)//startdate(9:10)//startdate(12:13)//&
          startdate(15:16)//".LDASIN_DOMAIN"//hgrid
     else
          inflnm = trim(indir)//"/"//&
          startdate(1:4)//startdate(6:7)//startdate(9:10)//startdate(12:13)//&
          ".LDASIN_DOMAIN"//hgrid
     endif

#else

     inflnm = trim(indir)//"/"//&
          startdate(1:4)//startdate(6:7)//startdate(9:10)//startdate(12:13)//&
          ".LDASIN_DOMAIN"//hgrid
#endif

     CALL READINIT_HRLDAS(HRLDAS_SETUP_FILE, xstart, xend, ystart, yend,  &
          NSOIL, DZS, OLDDATE, LDASIN_VERSION, SMOIS,       &
          TSLB, CANWAT, TSK, SNOW, SNOWH, FNDSNOWH)

     VEGFRA    =  undefined_real
     LAI       =  undefined_real
     GVFMIN    =  undefined_real
     GVFMAX    =  undefined_real

     CALL READVEG_HRLDAS(HRLDAS_SETUP_FILE, xstart, xend, ystart, yend,  &
          OLDDATE, IVGTYP, VEGFRA, LAI, GVFMIN, GVFMAX)

#ifdef WRF_HYDRO
  else   !  HRLDAS_ini_typ

#ifdef MPP_LAND
    call HYDRO_HRLDAS_ini_mpp   &
#else
    call HYDRO_HRLDAS_ini   &
#endif
       (trim(hrldas_setup_file), xend-xstart+1,yend-ystart+1, &
       nsoil,SMOIS(:,1:NSOIL,:),TSLB(:,1:NSOIL,:),SH2O(:,1:NSOIL,:), CANWAT, TSK,SNOW,SNOWH,lai,VEGFRA,IVGTYP,FNDSNOWH)

       greenfrac = 0.0
#ifdef MPP_LAND
       call get_greenfrac_mpp &
#else
       call get_greenfrac &
#endif
        (trim(GEO_STATIC_FLNM),greenfrac, ix, jx, olddate, SHDMAX)
       !yw GVFMAX = maxval(greenfrac)
  endif
#endif

  !     SNOW = SNOW * 1000.    ! Convert snow water equivalent to mm. MB: remove v3.7

  FNDSOILW = .FALSE.
  CALL NOAHMP_INIT(    LLANDUSE,     SNOW,    SNOWH,   CANWAT,   ISLTYP,   IVGTYP, XLAT, &   ! call from WRF phys_init
                 TSLB,    SMOIS,     SH2O,      DZS, FNDSOILW, FNDSNOWH, &
                  TSK,  ISNOWXY,     TVXY,     TGXY, CANICEXY,      TMN,     XICE, &
             CANLIQXY,    EAHXY,    TAHXY,     CMXY,     CHXY,                     &
               FWETXY, SNEQVOXY, ALBOLDXY,  QSNOWXY, QRAINXY,  WSLAKEXY,    ZWTXY,     WAXY, &
                 WTXY,   TSNOXY,  ZSNSOXY,  SNICEXY,  SNLIQXY, LFMASSXY, RTMASSXY, &
             STMASSXY,   WOODXY, STBLCPXY, FASTCPXY,   XSAIXY, LAI,                    &
              GRAINXY,    GDDXY,                                                   &
             CROPTYPE,  CROPCAT,                                                   &
             irnumsi  ,irnummi  ,irnumfi  ,irwatsi,    &
             irwatmi  ,irwatfi  ,ireloss  ,irsivol,    &
             irmivol  ,irfivol  ,irrsplh  ,            &
               T2MVXY,   T2MBXY, CHSTARXY,                                         &
                NSOIL,  .false.,                                                   &
               .true.,runoff_option, crop_option, irrigation_option, irrigation_method,  &
               sf_urban_physics,                         &  ! urban scheme
               ids,ide+1, jds,jde+1, kds,kde,                &  ! domain
               ims,ime, jms,jme, kms,kme,                &  ! memory
               its,ite, jts,jte, kts,kte                 &  ! tile
                  ,smoiseq  ,smcwtdxy ,rechxy, deeprechxy, qtdrain, areaxy ,dx, dy, msftx, msfty,&
                  wtddt    ,stepwtd  ,dtbl  ,qrfsxy ,qspringsxy  ,qslatxy,                  &
                  fdepthxy ,terrain       ,riverbedxy ,eqzwt ,rivercondxy ,pexpxy,              &
                  rechclim ,gecros_state                 &
                  )

      if(iopt_run == 5) then
           call groundwater_init ( grid,                                     &
              num_soil_layers, dzs, isltyp, ivgtyp, wtddt ,                  &
              fdepthxy   , terrain , riverbedxy, eqzwt     ,                 &
              rivercondxy, pexpxy  , areaxy    , zwtxy     ,                 &
              smois      , sh2o    , smoiseq   , smcwtdxy  ,                 &
              QLATXY     , qslatxy , QRFXY     , qrfsxy    ,                 &
              deeprechxy , rechxy  , QSPRINGXY , qspringsxy,                 &
              rechclim   ,                                                   &
              ids,ide, jds,jde, kds,kde,                                     &
              ims,ime, jms,jme, kms,kme,                                     &
              ims,ime, jms,jme, kms,kme,                                     &
              its,ite, jts,jte, kts,kte )
      endif

   TAUSSXY = 0.0   ! Need to be added to _INIT later

      IF(SF_URBAN_PHYSICS > 0 ) THEN  !urban

              CALL urban_param_init(DZR,DZB,DZG,num_soil_layers,sf_urban_physics,use_wudapt_lcz)

             CALL urban_var_init(ISURBAN,TSK,TSLB,TMN,IVGTYP,                        & !urban
                   ims,ime,jms,jme,kms,kme,num_soil_layers,                          & !urban
                   LCZ_1_TABLE,LCZ_2_TABLE,LCZ_3_TABLE,LCZ_4_TABLE,LCZ_5_TABLE,      & !urban
                   LCZ_6_TABLE,LCZ_7_TABLE,LCZ_8_TABLE,LCZ_9_TABLE,LCZ_10_TABLE,     & !urban
                   LCZ_11_TABLE,                                                     & !urban
                   .false.,sf_urban_physics,                                         & !urban
                   XXXR_URB2D,    XXXB_URB2D,    XXXG_URB2D,XXXC_URB2D,              & !urban
                     TR_URB2D,      TB_URB2D,      TG_URB2D,  TC_URB2D, QC_URB2D,    & !urban
                    TRL_URB3D,     TBL_URB3D,     TGL_URB3D,                         & !urban
                     SH_URB2D,      LH_URB2D,       G_URB2D,  RN_URB2D, TS_URB2D,    & !urban
                num_urban_ndm,  urban_map_zrd,  urban_map_zwd, urban_map_gd,         & !I multi-layer urban
                 urban_map_zd,  urban_map_zdf,   urban_map_bd, urban_map_wd,         & !I multi-layer urban
                urban_map_gbd,  urban_map_fbd, urban_map_zgrd,                       & !I multi-layer urban
                 num_urban_hi,                                                       & !urban
                    TRB_URB4D,     TW1_URB4D,     TW2_URB4D, TGB_URB4D,              & !urban
                   TLEV_URB3D,    QLEV_URB3D,                                        & !urban
                 TW1LEV_URB3D,  TW2LEV_URB3D,                                        & !urban
                  TGLEV_URB3D,   TFLEV_URB3D,                                        & !urban
                  SF_AC_URB3D,   LF_AC_URB3D,   CM_AC_URB3D,                         & !urban
                 SFVENT_URB3D,  LFVENT_URB3D,                                        & !urban
                 SFWIN1_URB3D,  SFWIN2_URB3D,                                        & !urban
                   SFW1_URB3D,    SFW2_URB3D,     SFR_URB3D, SFG_URB3D,              & !urban
                  EP_PV_URB3D,    T_PV_URB3D,                                        & !GRZ
                    TRV_URB4D,      QR_URB4D,     QGR_URB3D, TGR_URB3D,              & !GRZ
                  DRAIN_URB4D, DRAINGR_URB3D,    SFRV_URB3D,                         & !GRZ
                   LFRV_URB3D,     DGR_URB3D,      DG_URB3D, LFR_URB3D, LFG_URB3D,   & !GRZ
                        SMOIS,                                                       & !GRZ
                     LP_URB2D,      HI_URB2D,      LB_URB2D,                         & !urban
                    HGT_URB2D,      MH_URB2D,    STDH_URB2D,                         & !urban
                     LF_URB2D,                                                       & !urban
                   CMCR_URB2D,     TGR_URB2D,    TGRL_URB3D, SMR_URB3D,              & !urban
                  DRELR_URB2D,   DRELB_URB2D,   DRELG_URB2D,                         & !urban
                FLXHUMR_URB2D, FLXHUMB_URB2D, FLXHUMG_URB2D,                         & !urban
                      A_U_BEP,       A_V_BEP,       A_T_BEP,   A_Q_BEP,              & !multi-layer urban
                      A_E_BEP,       B_U_BEP,       B_V_BEP,                         & !multi-layer urban
                      B_T_BEP,       B_Q_BEP,       B_E_BEP,   DLG_BEP,              & !multi-layer urban
                     DL_U_BEP,        SF_BEP,        VL_BEP,                         & !multi-layer urban
                    FRC_URB2D,   UTYPE_URB2D, use_wudapt_lcz)                          !urban

            max_utype_urb2d = maxval(UTYPE_URB2D)*1.0
            IF (use_wudapt_lcz.eq.0 .and. max_utype_urb2d.gt.3.0) THEN  !new LCZ
              CALL wrf_error_fatal &
              ('USING 10 WUDAPT LCZ WITHOUT URBPARM_LCZ.TBL. SET USE_WUDAPT_LCZ=1')
            ENDIF
            IF (use_wudapt_lcz.eq.1 .and. max_utype_urb2d.le.3.0) THEN  ! new LCZ
              CALL wrf_error_fatal &
              ('USING URBPARM_LCZ.TBL WITH OLD 3 URBAN CLASSES. SET USE_WUDAPT_LCZ=0')
            ENDIF

       ENDIF !end of if urban

       IF(SF_URBAN_PHYSICS > 1 ) THEN  !urban
         do i = 1, num_urban_atmosphere-1
           dz_urban(:,i,:) = urban_atmosphere_thickness  ! thickness of full levels
         end do
        dz_urban(:,num_urban_atmosphere,:) =                                   &
               2*(zlvl - urban_atmosphere_thickness*(num_urban_atmosphere-1))
    print*, dz_urban(1,:,1)
  ENDIF

      CALL NOAHMP_MOSAIC_INIT(XLAND, SNOW,    SNOWH,   CANWAT,   ISLTYP,   IVGTYP, XLAT, &   ! call from WRF phys_init
                    TSLB,    SMOIS,     SH2O,      DZS, FNDSOILW, FNDSNOWH,  &
                     TSK,  ISNOWXY,     TVXY,     TGXY, CANICEXY,      TMN,     XICE, &
                CANLIQXY,    EAHXY,    TAHXY,     CMXY,     CHXY,                     &
                  FWETXY, SNEQVOXY, ALBOLDXY,  QSNOWXY, QRAINXY,  WSLAKEXY,    ZWTXY,     WAXY, &
                    WTXY,   TSNOXY,  ZSNSOXY,  SNICEXY,  SNLIQXY, LFMASSXY, RTMASSXY, &
                STMASSXY,   WOODXY, STBLCPXY, FASTCPXY,   XSAIXY, LAI,                    &
                 GRAINXY,    GDDXY,                                                   &
                 CROPTYPE, CROPCAT,                                                   &
                 irnumsi, irnummi, irnumfi, irwatsi,                                  &
                 irwatmi, irwatfi, ireloss, irsivol,                                  &
                 irmivol, irfivol, irrsplh,                                           &
                  T2MVXY,   T2MBXY, CHSTARXY,                                         &
                   NSOIL,  .false.,                                                   &
                  .true.,runoff_option, crop_option, irrigation_option, irrigation_method,  &
                  sf_urban_physics,ISWATER, ISICE,                         &  ! urban scheme
                  ISURBAN, 0,                              &
                  ids,ide+1, jds,jde+1, kds,kde,                &  ! domain
                  ims,ime, jms,jme, kms,kme,                &  ! memory
                  its,ite, jts,jte, kts,kte,                 &  ! tile
                  smoiseq, smcwtdxy, rechxy, deeprechxy,                               &
                  LANDUSEF, landusef2, number_land_use_catagories, IOPT_MOSAIC,                     &      ! Added by Aaron A. **IMPORTANT FOR THE FOR LOOPS OF FRACTIONAL LAND USE
                  mosaic_cat_index, number_mosaic_catagories,                           &      ! Added by Aaron A. **IMPORTANT FOR THE FOR LOOPS OF FRACTIONAL LAND USE
                  IOPT_HUE,                                                      &
                  TSK_mosaic, TSLB_mosaic, SMOIS_mosaic, SH2O_mosaic,                   &      ! Added by Aaron A.
                  CANWAT_mosaic, SNOW_mosaic, SNOWH_mosaic,                             &      ! Added by Aaron A.
                  isnowxy_mosaic, tvxy_mosaic, tgxy_mosaic, canicexy_mosaic,            &      ! Added by Aaron A.
                  canliqxy_mosaic, eahxy_mosaic, tahxy_mosaic,              &      ! Added by Aaron A.
                  cmxy_mosaic, chxy_mosaic, fwetxy_mosaic, sneqvoxy_mosaic,             &      ! Added by Aaron A.
                  alboldxy_mosaic, qsnowxy_mosaic, qrainxy_mosaic, wslakexy_mosaic, zwtxy_mosaic,       &      ! Added by Aaron A.
                  waxy_mosaic, wtxy_mosaic, tsnoxy_mosaic, zsnsoxy_mosaic,              &      ! Added by Aaron A.
                  snicexy_mosaic, snliqxy_mosaic, lfmassxy_mosaic, rtmassxy_mosaic,     &      ! Added by Aaron A.
                  stmassxy_mosaic, woodxy_mosaic, stblcpxy_mosaic, fastcpxy_mosaic,     &      ! Added by Aaron A.
                  xsaixy_mosaic, xlai_mosaic,                                           &      ! Added by Aaron A.

                  IRNUMSI_mosaic, IRNUMMI_mosaic, IRNUMFI_mosaic, IRWATSI_mosaic,       &
                  IRWATMI_mosaic, IRWATFI_mosaic, IRELOSS_mosaic, IRSIVOL_mosaic,       &
                  IRMIVOL_mosaic, IRFIVOL_mosaic, IRRSPLH_mosaic,                       &

                  smoiseq_mosaic,  smcwtdxy_mosaic, rechxy_mosaic, deeprechxy_mosaic,   &

                  TR_URB2D_mosaic, TB_URB2D_mosaic,                                     &      ! Added by Aaron A.
                  TG_URB2D_mosaic, TC_URB2D_mosaic, QC_URB2D_mosaic,                    &      ! Added by Aaron A.
                  TRL_URB3D_mosaic, TBL_URB3D_mosaic,                                   &      ! Added by Aaron A.
                  TGL_URB3D_mosaic, SH_URB2D_mosaic, LH_URB2D_mosaic,                   &      ! Added by Aaron A.
                  G_URB2D_mosaic, RN_URB2D_mosaic, TS_URB2D_mosaic,                     &      ! Added by Aaron A.
                  CMR_SFCDIF_mosaic, CHR_SFCDIF_mosaic, CMC_SFCDIF_mosaic,              &
                  CHC_SFCDIF_mosaic, CMGR_SFCDIF_mosaic, CHGR_SFCDIF_mosaic,            &
                  XXXR_URB2D_mosaic, XXXB_URB2D_mosaic, XXXG_URB2D_mosaic,              &
                  XXXC_URB2D_mosaic,                                                    &
                  CMCR_URB2D_mosaic, TGR_URB2D_mosaic,                                  &
                  TGRL_URB3D_mosaic, SMR_URB3D_mosaic,                                  &
                  DRELR_URB2D_mosaic, DRELB_URB2D_mosaic, DRELG_URB2D_mosaic,           &
                  FLXHUMR_URB2D_mosaic, FLXHUMB_URB2D_mosaic, FLXHUMG_URB2D_mosaic,     &
                  DETENTION_STORAGEXY_mosaic,                                           &
                  Z0, ZNT_mosaic, Z0_mosaic )                                           ! Added by Aaron A.

   endif !end mosaic if else statement
  endif !end of the restart or non-restart call
  WRITE(*,*) "POST INITILIZATION"

  NTIME=(KHOUR)*3600./nint(dtbl)*(spinup_loops+1)
  spinup_loop = 0
  reset_spinup_date = .false.

  print*, "NTIME = ", NTIME , "KHOUR=",KHOUR,"dtbl = ", dtbl

  call system_clock(count=clock_count_1)   ! Start a timer


#ifdef WRF_HYDRO
   allocate( infxsrt   (ix,jx) )
   allocate( sfcheadrt (ix,jx) )
   allocate( soldrain  (ix,jx) )
   allocate( etpnd     (ix,jx) )
   allocate( prcp0     (ix,jx) )
   allocate( qtiledrain(ix,jx) )
   allocate( ZWATBLE2D (ix,jx) )

   prcp0     = 0
   sfcheadrt = 0.0
   infxsrt   = 0.0
   etpnd     = 0.0
   soldrain  = 0.0
   qtiledrain= 0.0
   ZWATBLE2D = 0.0

   allocate(zsoil (NSOIL))
   zsoil = 0

   zsoil(1) = -1*soil_thick_input(1)
   do kk = 2, NSOIL
      zsoil(kk) = zsoil(kk-1)-soil_thick_input(kk)
   end do

   print*, "zsoil/soil_thick_input = ", soil_thick_input(1:NSOIL)

   call hrldas_drv_HYDRO_ini(TSLB(:,1:NSOIL,:),SMOIS(:,1:NSOIL,:),SH2O(:,1:NSOIL,:), &
         infxsrt,sfcheadrt,soldrain,qtiledrain,ZWATBLE2D,ix,jx,NSOIL, SMOIS,real(noah_timestep),&
         olddate,zsoil(1:NSOIL))

   if(finemesh .ne. 0 ) then
       if(restart_flag) then
          NTIME_out =  10
       else
          NTIME_out =  1
       endif
       return
   endif
#endif

   NTIME_out = NTIME

end subroutine land_driver_ini

!===============================================================================
  subroutine land_driver_exe(itime)
     implicit  none
     integer :: itime          ! timestep loop

if(IOPT_MOSAIC.eq.0) then !this check looks and sees if we are going to be doing
  ! a mosaic tiling scheme call.
!---------------------------------------------------------------------------------
! Read the forcing data.
!---------------------------------------------------------------------------------

! For HRLDAS, we're assuming (for now) that each time period is in a
! separate file.  So we can open a new one right now.

     inflnm = trim(indir)//"/"//&
          olddate(1:4)//olddate(6:7)//olddate(9:10)//olddate(12:13)//&
          ".LDASIN_DOMAIN"//hgrid

     ! Build a filename template
     inflnm_template = trim(indir)//"/<date>.LDASIN_DOMAIN"//hgrid

#ifdef MPP_LAND
     call mpp_land_bcast_char(19,OLDDATE(1:19))
#endif

#ifdef WRF_HYDRO

#ifdef MPP_LAND
     call mpp_land_bcast_char(19,forcDate(1:19))
#endif
     if(finemesh .ne. 0) goto 991

     if(forc_typ .eq. 0) then
        CALL READFORC_HRLDAS(INFLNM_TEMPLATE, FORCING_TIMESTEP, OLDDATE,  &
             XSTART, XEND, YSTART, YEND,                                  &
       forcing_name_T,forcing_name_Q,forcing_name_U,forcing_name_V,forcing_name_P, &
       forcing_name_LW,forcing_name_SW,forcing_name_PR,forcing_name_SN, &
       T_PHY(:,1,:),QV_CURR(:,1,:),U_PHY(:,1,:),V_PHY(:,1,:),          &
       P8W(:,1,:), GLW, SWDOWN, RAINBL_tmp, SNOWBL, VEGFRA, update_veg, LAI, update_lai, reset_spinup_date, startdate)
        VEGFRA = VEGFRA * 100.
     else
        if(olddate == forcDate) then
           LAI_tmp = LAI
           CALL HYDRO_frocing_drv(trim(indir), forc_typ,snow_assim,olddate,                      &
               xstart, xend, ystart, yend,    &
               T_PHY(:,1,:),QV_CURR(:,1,:),U_PHY(:,1,:),V_PHY(:,1,:),P8W(:,1,:),    &
               GLW,SWDOWN,RAINBL_tmp,LAI,VEGFRA,SNOWH,ITIME,FORCING_TIMESTEP,prcp0)

           call geth_newdate(newdate, forcDate, FORCING_TIMESTEP)
           forcDate = newdate
        endif
     endif

#else
     CALL READFORC_HRLDAS(INFLNM_TEMPLATE, FORCING_TIMESTEP, OLDDATE,  &
          XSTART, XEND, YSTART, YEND,                                  &
       forcing_name_T,forcing_name_Q,forcing_name_U,forcing_name_V,forcing_name_P, &
       forcing_name_LW,forcing_name_SW,forcing_name_PR,forcing_name_SN, &
       T_PHY(:,1,:),QV_CURR(:,1,:),U_PHY(:,1,:),V_PHY(:,1,:),          &
       P8W(:,1,:), GLW, SWDOWN, RAINBL_tmp, SNOWBL, VEGFRA, update_veg, LAI, update_lai, reset_spinup_date, startdate)
       if(maxval(VEGFRA) <= 1 ) then
           VEGFRA = VEGFRA * 100.   ! added for input vegfra as a fraction (0~1)
       endif

#endif

991  continue

     where(XLAND > 1.5)   T_PHY(:,1,:) = 0.0  ! Prevent some overflow problems with ifort compiler [MB:20150812]
     where(XLAND > 1.5)   U_PHY(:,1,:) = 0.0
     where(XLAND > 1.5)   V_PHY(:,1,:) = 0.0
     where(XLAND > 1.5) QV_CURR(:,1,:) = 0.0
     where(XLAND > 1.5)     P8W(:,1,:) = 0.0
     where(XLAND > 1.5)     GLW        = 0.0
     where(XLAND > 1.5)  SWDOWN        = 0.0
     where(XLAND > 1.5) RAINBL_tmp     = 0.0
     where(XLAND > 1.5) SNOWBL         = 0.0

     QV_CURR(:,1,:) = QV_CURR(:,1,:)/(1.0 - QV_CURR(:,1,:))  ! Assuming input forcing are specific hum.;
                                                             ! WRF wants mixing ratio at driver level
     P8W(:,2,:)     = P8W(:,1,:)      ! WRF uses lowest two layers
     T_PHY(:,2,:)   = T_PHY(:,1,:)    ! Only pressure is needed in two layer but fill the rest
     U_PHY(:,2,:)   = U_PHY(:,1,:)    !
     V_PHY(:,2,:)   = V_PHY(:,1,:)    !
     QV_CURR(:,2,:) = QV_CURR(:,1,:)  !
     RAINBL = RAINBL_tmp * DTBL       ! RAINBL in WRF is [mm]
     SNOWBL = SNOWBL * DTBL           !
     SR         = 0.0                 ! Will only use component if opt_snf=4
     RAINCV     = 0.0
     RAINNCV    = RAINBL
     RAINSHV    = 0.0
     SNOWNCV    = SNOWBL
     GRAUPELNCV = 0.0
     HAILNCV    = 0.0
     DZ8W = 2*ZLVL                    ! 2* to be consistent with WRF model level

     SWDDIR = SWDOWN * 0.7  ! following noahmplsm ATM 70% direct radiation
     SWDDIF = SWDOWN * 0.3  ! following noahmplsm ATM 30% diffuse radiation

     IF(SF_URBAN_PHYSICS == 2 .or. SF_URBAN_PHYSICS == 3) THEN  ! BEP or BEM urban models
       where(XLAND < 1.5)     p_urban(:,kde,:) = p8w(:,1,:)
       where(XLAND < 1.5)   rho_urban(:,kde,:) = p8w(:,1,:)/287.04/t_phy(:,1,:)
       where(XLAND < 1.5) theta_urban(:,kde,:) = t_phy(:,1,:)*(100000.0/p8w(:,1,:))**(0.285714)
       where(XLAND < 1.5)     u_urban(:,kde,:) = u_phy(:,1,:)
       where(XLAND < 1.5)     v_urban(:,kde,:) = v_phy(:,1,:)
       height_urban = 0.5*(dz_urban(1,kde,1) + dz_urban(1,kde-1,1))

!      IF(ITIME == 0) THEN  ! BEP or BEM urban models
       DO I = KDE-1, KDS, -1
         where(XLAND < 1.5)     p_urban(:,i,:) = p8w(:,1,:)*exp(9.8*height_urban/287.04/t_phy(:,1,:))
         where(XLAND < 1.5)   rho_urban(:,i,:) = p_urban(:,i,:)/287.04/t_phy(:,1,:)
         where(XLAND < 1.5) theta_urban(:,i,:) = t_phy(:,1,:)*(100000.0/p_urban(:,i,:))**(0.285714)
         where(XLAND < 1.5)     u_urban(:,i,:) = u_urban(:,kde,:)*log(zlvl-height_urban)/log(zlvl)
         where(XLAND < 1.5)     v_urban(:,i,:) = v_urban(:,kde,:)*log(zlvl-height_urban)/log(zlvl)
         height_urban = height_urban + urban_atmosphere_thickness
       END DO
!      ENDIF
     ENDIF

!------------------------------------------------------------------------
! Noah-MP updates we can do before spatial loop.
!------------------------------------------------------------------------

   ! create a few fields that are IN in WRF - coszen, julian,yr

    DO J = YSTART,YEND
    DO I = XSTART,XEND
      IF(SF_URBAN_PHYSICS == 0) THEN
        CALL CALC_DECLIN(OLDDATE(1:19),XLAT(I,J), XLONG(I,J),COSZEN(I,J),JULIAN)
      ELSE
        CALL CALC_DECLIN(OLDDATE(1:19),XLAT(I,J), XLONG(I,J),COSZEN(I,J),JULIAN, &
                         HRANG(I,J), DECLIN, GMT, JULDAY)
      ENDIF
      IF(IOPT_IRR > 0 ) THEN
         CALL LOCAL_TIME(OLDDATE(1:19), XLONG(I,J), LOCTIM(I,J))
      ENDIF
    END DO
    END DO

    READ(OLDDATE(1:4),*)  YR
    YEARLEN = 365                      ! find length of year for phenology (also S Hemisphere)
    if (mod(YR,4) == 0) then
       YEARLEN = 366
       if (mod(YR,100) == 0) then
          YEARLEN = 365
          if (mod(YR,400) == 0) then
             YEARLEN = 366
          endif
       endif
    endif

    IF (ITIME == 1 .AND. .NOT. RESTART_FLAG ) THEN
      EAHXY = (P8W(:,1,:)*QV_CURR(:,1,:))/(0.622+QV_CURR(:,1,:)) ! Initial guess only.
      TAHXY = T_PHY(:,1,:)                                       ! Initial guess only.
      CHXY = 0.1
      CMXY = 0.1
    ENDIF

!------------------------------------------------------------------------
! Skip model call at t=1 since initial conditions are at start time; First model time is +1
!------------------------------------------------------------------------

   IF (ITIME > 0) THEN

!------------------------------------------------------------------------
! Call to Noah-MP driver same as surface_driver
!------------------------------------------------------------------------
     sflx_count_sum = 0 ! Timing

   ! Timing information for SFLX:
    call system_clock(count=count_before_sflx, count_rate=clock_rate)

         CALL noahmplsm(ITIMESTEP,       YR,   JULIAN,   COSZEN,  XLAT,XLONG, &
	           DZ8W,     DTBL,      DZS,     NUM_SOIL_LAYERS,         DX, &
		 IVGTYP,   ISLTYP,   VEGFRA,   GVFMAX,       TMN,             &
		  XLAND,     XICE,     XICE_THRESHOLD,   CROPCAT,             &
               PLANTING,  HARVEST,SEASON_GDD,                                 &
                  IDVEG, IOPT_CRS, IOPT_BTR, IOPT_RUN,  IOPT_SFC,   IOPT_FRZ, &
	       IOPT_INF, IOPT_RAD, IOPT_ALB, IOPT_SNF, IOPT_TBOT,   IOPT_STC, &
	       IOPT_GLA, IOPT_RSF,IOPT_SOIL,IOPT_PEDO,IOPT_CROP,    IOPT_IRR, &
               IOPT_IRRM,IOPT_INFDV,IOPT_TDRN,soiltstep,IOPT_MOSAIC, IOPT_HUE, IZ0TLND, sf_urban_physics,      &
	       SOILCOMP,  SOILCL1,  SOILCL2,   SOILCL3,  SOILCL4,             &
		  T_PHY,  QV_CURR,    U_PHY,    V_PHY,    SWDOWN,     SWDDIR, &
                 SWDDIF,      GLW,                                            &
		    P8W,   RAINBL,       SR,                                  &
                  IRFRACT, SIFRACT,   MIFRACT,  FIFRACT,                      & ! IN : Irrigation fractions
		    TSK,      HFX,      QFX,       LH,    GRDFLX,     SMSTAV, &
		 SMSTOT,SFCRUNOFF, UDRUNOFF,   ALBEDO,     SNOWC,      SMOIS, &
		   SH2O,     TSLB,     SNOW,    SNOWH,    CANWAT,     ACSNOM, &
		 ACSNOW,    EMISS,     QSFC,                                  &
 		     Z0,      ZNT,                                            & ! IN/OUT LSM eqv
                IRNUMSI,  IRNUMMI,  IRNUMFI,   IRWATSI,  IRWATMI,    IRWATFI, & ! IN/OUT Irrigation
                IRELOSS,  IRSIVOL,  IRMIVOL,   IRFIVOL,  IRRSPLH,   LLANDUSE, & ! IN/OUT Irrigation
		ISNOWXY,     TVXY,     TGXY, CANICEXY,  CANLIQXY,      EAHXY, &
		  TAHXY,     CMXY,     CHXY,   FWETXY,  SNEQVOXY,   ALBOLDXY, &
		QSNOWXY,  QRAINXY, WSLAKEXY,    ZWTXY, WAXY,WTXY,     TSNOXY, &
		ZSNSOXY,  SNICEXY,  SNLIQXY, LFMASSXY,  RTMASSXY,   STMASSXY, &
		 WOODXY, STBLCPXY, FASTCPXY,      LAI,    XSAIXY,    TAUSSXY, &
	        SMOISEQ, SMCWTDXY,DEEPRECHXY,  RECHXY,   GRAINXY,      GDDXY,PGSXY, & ! IN/OUT Noah MP only
                GECROS_STATE,                                                 & ! IN/OUT gecros model
                QTDRAIN,   TD_FRACTION,                                       & ! IN/OUT tile drainage
	         T2MVXY,   T2MBXY,   Q2MVXY,   Q2MBXY,                        &
                 TRADXY,    NEEXY,    GPPXY,    NPPXY,    FVEGXY,    RUNSFXY, &
	        RUNSBXY,   ECANXY,   EDIRXY,  ETRANXY,     FSAXY,     FIRAXY, &
                 APARXY,    PSNXY,    SAVXY,    SAGXY,   RSSUNXY,    RSSHAXY, &
                 BGAPXY,   WGAPXY,    TGVXY,    TGBXY,     CHVXY,      CHBXY, &
		  SHGXY,    SHCXY,    SHBXY,    EVGXY,     EVBXY,      GHVXY, &
		  GHBXY,    IRGXY,    IRCXY,    IRBXY,      TRXY,      EVCXY, &
	        CHLEAFXY,   CHUCXY,   CHV2XY,   CHB2XY,        RS,             &
                QINTSXY   ,QINTRXY   ,QDRIPSXY   ,&
                QDRIPRXY  ,QTHROSXY  ,QTHRORXY   ,&
                QSNSUBXY  ,QSNFROXY  ,QSUBCXY    ,&
                QFROCXY   ,QEVACXY   ,QDEWCXY    ,QFRZCXY   ,QMELTCXY   ,&
                QSNBOTXY  ,QMELTXY   ,PONDINGXY  ,PAHXY     ,PAHGXY, PAHVXY, PAHBXY,&
                FPICEXY,RAINLSM,SNOWLSM,FORCTLSM ,FORCQLSM,FORCPLSM,FORCZLSM,FORCWLSM,&
                ACC_SSOILXY, ACC_QINSURXY, ACC_QSEVAXY, ACC_ETRANIXY, EFLXBXY, &
                SOILENERGY, SNOWENERGY, CANHSXY, &
                ACC_DWATERXY, ACC_PRCPXY, ACC_ECANXY, ACC_ETRANXY, ACC_EDIRXY, &
!                BEXP_3D,SMCDRY_3D,SMCWLT_3D,SMCREF_3D,SMCMAX_3D,             &
!		 DKSAT_3D,DWSAT_3D,PSISAT_3D,QUARTZ_3D,                       &
!		 REFDK_2D,REFKDT_2D,                                          &
!                IRR_FRAC_2D,IRR_HAR_2D,IRR_LAI_2D,IRR_MAD_2D,FILOSS_2D,      &
!                SPRIR_RATE_2D,MICIR_RATE_2D,FIRTFAC_2D,IR_RAIN_2D,           &
!		 BVIC_2D,AXAJ_2D,BXAJ_2D,XXAJ_2D,BDVIC_2D,GDVIC_2D,BBVIC_2D,  &

#ifdef WRF_HYDRO
                 sfcheadrt,INFXSRT,soldrain,qtiledrain,ZWATBLE2D,             &    !O
#endif
                ids,ide, jds,jde, kds,kde,                      &
                ims,ime, jms,jme, kms,kme,                      &
                its,ite, jts,jte, kts,kte)!,        &
! variables below are optional
                 !MP_RAINC =  RAINCV, MP_RAINNC =    RAINNCV, MP_SHCV = RAINSHV,&
		!MP_SNOW  = SNOWNCV, MP_GRAUP  = GRAUPELNCV, MP_HAIL = HAILNCV )

          call system_clock(count=count_after_sflx, count_rate=clock_rate)
          sflx_count_sum = sflx_count_sum + ( count_after_sflx - count_before_sflx )

         IF(SF_URBAN_PHYSICS == 2 .or. SF_URBAN_PHYSICS == 3) THEN  ! BEP or BEM urban models
           dz8w = dz_urban
           u_phy = u_urban
           v_phy = v_urban
         ENDIF

         IF(SF_URBAN_PHYSICS > 0 ) THEN  !urban

	   call noahmp_urban (sf_urban_physics,     NUM_SOIL_LAYERS,     IVGTYP,ITIMESTEP,  & ! IN : Model configuration
                               DTBL,         COSZEN,           XLAT,                        & ! IN : Time/Space-related
                              T_PHY,        QV_CURR,          U_PHY,      V_PHY,   SWDOWN,  & ! IN : Forcing
                             SWDDIR,         SWDDIF,                                        &
		                GLW,            P8W,         RAINBL,       DZ8W,      ZNT,  & ! IN : Forcing
                                TSK,            HFX,            QFX,         LH,   GRDFLX,  & ! IN/OUT : LSM
		             ALBEDO,          EMISS,           QSFC,                        & ! IN/OUT : LSM
                            ids,ide,        jds,jde,        kds,kde,                        &
                            ims,ime,        jms,jme,        kms,kme,                        &
                            its,ite,        jts,jte,        kts,kte,                        &
                         cmr_sfcdif,     chr_sfcdif,     cmc_sfcdif,                        &
	                 chc_sfcdif,    cmgr_sfcdif,    chgr_sfcdif,                        &
                           tr_urb2d,       tb_urb2d,       tg_urb2d,                        & !H urban
	                   tc_urb2d,       qc_urb2d,       uc_urb2d,                        & !H urban
                         xxxr_urb2d,     xxxb_urb2d,     xxxg_urb2d, xxxc_urb2d,            & !H urban
                          trl_urb3d,      tbl_urb3d,      tgl_urb3d,                        & !H urban
                           sh_urb2d,       lh_urb2d,        g_urb2d,   rn_urb2d,  ts_urb2d, & !H urban
                         psim_urb2d,     psih_urb2d,      u10_urb2d,  v10_urb2d,            & !O urban
                       GZ1OZ0_urb2d,     AKMS_URB2D,                                        & !O urban
                          th2_urb2d,       q2_urb2d,      ust_urb2d,                        & !O urban
                             declin,          hrang,                                        & !I urban
                    num_roof_layers,num_wall_layers,num_road_layers,                        & !I urban
                                dzr,            dzb,            dzg,                        & !I urban
                         cmcr_urb2d,      tgr_urb2d,     tgrl_urb3d,  smr_urb3d,            & !H urban
                        drelr_urb2d,    drelb_urb2d,    drelg_urb2d,                        & !H urban
                      flxhumr_urb2d,  flxhumb_urb2d,  flxhumg_urb2d,                        & !H urban
                             julday,             yr,                                        & !H urban
                          frc_urb2d,    utype_urb2d,                                        & !I urban
                                chs,           chs2,           cqs2,                        & !H
                      num_urban_ndm,  urban_map_zrd,  urban_map_zwd, urban_map_gd,          & !I multi-layer urban
                       urban_map_zd,  urban_map_zdf,   urban_map_bd, urban_map_wd,          & !I multi-layer urban
                      urban_map_gbd,  urban_map_fbd, urban_map_zgrd,                        & !I multi-layer urban
                       num_urban_hi,                                                        & !I multi-layer urban
                          trb_urb4d,      tw1_urb4d,      tw2_urb4d,  tgb_urb4d,            & !H multi-layer urban
                         tlev_urb3d,     qlev_urb3d,                                        & !H multi-layer urban
                       tw1lev_urb3d,   tw2lev_urb3d,                                        & !H multi-layer urban
                        tglev_urb3d,    tflev_urb3d,                                        & !H multi-layer urban
                        sf_ac_urb3d,    lf_ac_urb3d,    cm_ac_urb3d,                        & !H multi-layer urban
                       sfvent_urb3d,   lfvent_urb3d,                                        & !H multi-layer urban
                       sfwin1_urb3d,   sfwin2_urb3d,                                        & !H multi-layer urban
                         sfw1_urb3d,     sfw2_urb3d,      sfr_urb3d,  sfg_urb3d,            & !H multi-layer urban
                        ep_pv_urb3d,     t_pv_urb3d,                                        & !RMS
                          trv_urb4d,       qr_urb4d,      qgr_urb3d,  tgr_urb3d,            & !RMS
                        drain_urb4d,  draingr_urb3d,     sfrv_urb3d, lfrv_urb3d,            & !RMS
                          dgr_urb3d,       dg_urb3d,      lfr_urb3d,  lfg_urb3d,            & !RMS
                           lp_urb2d,       hi_urb2d,       lb_urb2d,  hgt_urb2d,            & !H multi-layer urban
                           mh_urb2d,     stdh_urb2d,       lf_urb2d,                        & !SLUCM
                        theta_urban,      rho_urban,        p_urban,        ust,            & !I multi-layer urban
                                gmt,         julday,          XLONG,       XLAT,            & !I multi-layer urban
                            a_u_bep,        a_v_bep,        a_t_bep,    a_q_bep,            & !O multi-layer urban
                            a_e_bep,        b_u_bep,        b_v_bep,                        & !O multi-layer urban
                            b_t_bep,        b_q_bep,        b_e_bep,    dlg_bep,            & !O multi-layer urban
                           dl_u_bep,         sf_bep,         vl_bep)                          !O multi-layer urban

	 ENDIF

  IF(RUNOFF_OPTION.EQ.5.AND.MOD(ITIME,STEPWTD).EQ.0)THEN
           CALL wrf_message('calling WTABLE' )

!gmm update wtable from lateral flow and shed water to rivers
           CALL WTABLE_MMF_NOAHMP(                                        &
	       NUM_SOIL_LAYERS,  XLAND, XICE,       XICE_THRESHOLD, ISICE,    &
               ISLTYP,      SMOISEQ,    DZS,        WTDDT,                &
               FDEPTHXY,    AREAXY,     TERRAIN,    ISURBAN,    IVGTYP,   &
               RIVERCONDXY, RIVERBEDXY, EQZWT,      PEXPXY,               &
               SMOIS,       SH2O,       SMCWTDXY,   ZWTXY,                &
	       QLATXY, QRFXY, DEEPRECHXY, QSPRINGXY,                      &
               QSLATXY,     QRFSXY,     QSPRINGSXY, RECHXY,               &
               IDS,IDE, JDS,JDE, KDS,KDE,                                 &
               IMS,IME, JMS,JME, KMS,KME,                                 &
               ITS,ITE, JTS,JTE, KTS,KTE )

!         IF(SF_URBAN_PHYSICS.eq.1) THEN
!           DO j=jts,jte                                           !urban
!             DO i=its,ite                                         !urban
!             IF( IVGTYP(I,J) == ISURBAN    .or. IVGTYP(I,J) == LCZ_1 .or. IVGTYP(I,J) == LCZ_2 .or. &
!             IVGTYP(I,J) == LCZ_3      .or. IVGTYP(I,J) == LCZ_4 .or. IVGTYP(I,J) == LCZ_5 .or. &
!             IVGTYP(I,J) == LCZ_6      .or. IVGTYP(I,J) == LCZ_7 .or. IVGTYP(I,J) == LCZ_8 .or. &
!             IVGTYP(I,J) == LCZ_9      .or. IVGTYP(I,J) == LCZ_10 .or. IVGTYP(I,J) == LCZ_11 )THEN
!                 Q2(I,J)  = (FVEGXY(I,J)*Q2MVXY(I,J) + (1.-FVEGXY(I,J))*Q2MBXY(I,J))*(1.-FRC_URB2D(I,J)) +   &
!                             Q2_URB2D(I,J)*FRC_URB2D(I,J)
!                 T2(I,J)  = (FVEGXY(I,J)*T2MVXY(I,J) + (1.-FVEGXY(I,J))*T2MBXY(I,J))*(1.-FRC_URB2D(I,J)) +   &
!                             (TH2_URB2D(i,j)/((1.E5/PSFC(i,j))**RCP))*FRC_URB2D(I,J)
!                 TH2(I,J) = T2(i,j)*(1.E5/PSFC(i,j))**RCP
!                 U10(I,J)  = U10_URB2D(I,J)                       !urban
!                 V10(I,J)  = V10_URB2D(I,J)                       !urban
!                 PSIM(I,J) = PSIM_URB2D(I,J)                      !urban
!                 PSIH(I,J) = PSIH_URB2D(I,J)                      !urban
!                 GZ1OZ0(I,J) = GZ1OZ0_URB2D(I,J)                  !urban
!                 AKHS(I,J) = CHS(I,J)                             !urban
!                 AKMS(I,J) = AKMS_URB2D(I,J)                      !urban
!               END IF                                             !urban
!             ENDDO                                                !urban
!           ENDDO                                                  !urban
!         ENDIF

!         IF((SF_URBAN_PHYSICS.eq.2).OR.(SF_URBAN_PHYSICS.eq.3)) THEN
!           DO j=j_start(ij),j_end(ij)                             !urban
!             DO i=i_start(ij),i_end(ij)                           !urban
!             IF( IVGTYP(I,J) == ISURBAN    .or. IVGTYP(I,J) == LCZ_1 .or. IVGTYP(I,J) == LCZ_2 .or. &
!             IVGTYP(I,J) == LCZ_3      .or. IVGTYP(I,J) == LCZ_4 .or. IVGTYP(I,J) == LCZ_5 .or. &
!             IVGTYP(I,J) == LCZ_6      .or. IVGTYP(I,J) == LCZ_7 .or. IVGTYP(I,J) == LCZ_8 .or. &
!             IVGTYP(I,J) == LCZ_9      .or. IVGTYP(I,J) == LCZ_10 .or. IVGTYP(I,J) == LCZ_11 )THEN
!                T2(I,J)   = TH_PHY(i,1,j)/((1.E5/PSFC(I,J))**RCP) !urban
!                TH2(I,J) = TH_PHY(i,1,j) !urban
!                Q2(I,J)   = qv_curr(i,1,j)  !urban
!                U10(I,J)  = U_phy(I,1,J)                       !urban
!                V10(I,J)  = V_phy(I,1,J)                       !urban
!               END IF                                             !urban
!             ENDDO                                                !urban
!           ENDDO                                                  !urban
!         ENDIF

 ENDIF

!------------------------------------------------------------------------
! END of surface_driver consistent code
!------------------------------------------------------------------------

 ENDIF   ! SKIP FIRST TIMESTEP

  ELSEIF(IOPT_MOSAIC.eq.1) then !below this is the noah-mp mosiac scheme call Aaron A.

  inflnm = trim(indir)//"/"//&
       olddate(1:4)//olddate(6:7)//olddate(9:10)//olddate(12:13)//&
       ".LDASIN_DOMAIN"//hgrid

  ! Build a filename template
  inflnm_template = trim(indir)//"/<date>.LDASIN_DOMAIN"//hgrid

#ifdef MPP_LAND
  call mpp_land_bcast_char(19,OLDDATE(1:19))
#endif

#ifdef WRF_HYDRO

#ifdef MPP_LAND
  call mpp_land_bcast_char(19,forcDate(1:19))
#endif
  if(finemesh .ne. 0) goto 9955

  if(forc_typ .eq. 0) then
     CALL READFORC_HRLDAS(INFLNM_TEMPLATE, FORCING_TIMESTEP, OLDDATE,  &
          XSTART, XEND, YSTART, YEND,                                  &
    forcing_name_T,forcing_name_Q,forcing_name_U,forcing_name_V,forcing_name_P, &
    forcing_name_LW,forcing_name_SW,forcing_name_PR,forcing_name_SN, &
    T_PHY(:,1,:),QV_CURR(:,1,:),U_PHY(:,1,:),V_PHY(:,1,:),          &
    P8W(:,1,:), GLW, SWDOWN, RAINBL_tmp, SNOWBL, VEGFRA, update_veg, LAI, update_lai, reset_spinup_date, startdate)
     VEGFRA = VEGFRA * 100.
  else
     if(olddate == forcDate) then
        LAI_tmp = LAI
        CALL HYDRO_frocing_drv(trim(indir), forc_typ,snow_assim,olddate,                      &
            xstart, xend, ystart, yend,    &
            T_PHY(:,1,:),QV_CURR(:,1,:),U_PHY(:,1,:),V_PHY(:,1,:),P8W(:,1,:),    &
            GLW,SWDOWN,RAINBL_tmp,LAI,VEGFRA,SNOWH,ITIME,FORCING_TIMESTEP,prcp0)

        call geth_newdate(newdate, forcDate, FORCING_TIMESTEP)
        forcDate = newdate
     endif
  endif

#else
  CALL READFORC_HRLDAS(INFLNM_TEMPLATE, FORCING_TIMESTEP, OLDDATE,  &
       XSTART, XEND, YSTART, YEND,                                  &
    forcing_name_T,forcing_name_Q,forcing_name_U,forcing_name_V,forcing_name_P, &
    forcing_name_LW,forcing_name_SW,forcing_name_PR,forcing_name_SN, &
    T_PHY(:,1,:),QV_CURR(:,1,:),U_PHY(:,1,:),V_PHY(:,1,:),          &
    P8W(:,1,:), GLW, SWDOWN, RAINBL_tmp, SNOWBL, VEGFRA, update_veg, LAI, update_lai, reset_spinup_date, startdate)
    if(maxval(VEGFRA) <= 1 ) then
        VEGFRA = VEGFRA * 100.   ! added for input vegfra as a fraction (0~1)
    endif

#endif

9955  continue

  where(XLAND > 1.5)   T_PHY(:,1,:) = 0.0  ! Prevent some overflow problems with ifort compiler [MB:20150812]
  where(XLAND > 1.5)   U_PHY(:,1,:) = 0.0
  where(XLAND > 1.5)   V_PHY(:,1,:) = 0.0
  where(XLAND > 1.5) QV_CURR(:,1,:) = 0.0
  where(XLAND > 1.5)     P8W(:,1,:) = 0.0
  where(XLAND > 1.5)     GLW        = 0.0
  where(XLAND > 1.5)  SWDOWN        = 0.0
  where(XLAND > 1.5) RAINBL_tmp     = 0.0
  where(XLAND > 1.5) SNOWBL         = 0.0

  QV_CURR(:,1,:) = QV_CURR(:,1,:)/(1.0 - QV_CURR(:,1,:))  ! Assuming input forcing are specific hum.;
                                                          ! WRF wants mixing ratio at driver level
  P8W(:,2,:)     = P8W(:,1,:)      ! WRF uses lowest two layers
  T_PHY(:,2,:)   = T_PHY(:,1,:)    ! Only pressure is needed in two layer but fill the rest
  U_PHY(:,2,:)   = U_PHY(:,1,:)    !
  V_PHY(:,2,:)   = V_PHY(:,1,:)    !
  QV_CURR(:,2,:) = QV_CURR(:,1,:)  !
  RAINBL = RAINBL_tmp * DTBL       ! RAINBL in WRF is [mm]
  SNOWBL = SNOWBL * DTBL           !
  SR         = 0.0                 ! Will only use component if opt_snf=4
  RAINCV     = 0.0
  RAINNCV    = RAINBL
  RAINSHV    = 0.0
  SNOWNCV    = SNOWBL
  GRAUPELNCV = 0.0
  HAILNCV    = 0.0
  DZ8W = 2*ZLVL                    ! 2* to be consistent with WRF model level

  SWDDIR = SWDOWN * 0.7  ! following noahmplsm ATM 70% direct radiation
  SWDDIF = SWDOWN * 0.3  ! following noahmplsm ATM 30% diffuse radiation

  IF(SF_URBAN_PHYSICS == 2 .or. SF_URBAN_PHYSICS == 3) THEN  ! BEP or BEM urban models
    where(XLAND < 1.5)     p_urban(:,kde,:) = p8w(:,1,:)
    where(XLAND < 1.5)   rho_urban(:,kde,:) = p8w(:,1,:)/287.04/t_phy(:,1,:)
    where(XLAND < 1.5) theta_urban(:,kde,:) = t_phy(:,1,:)*(100000.0/p8w(:,1,:))**(0.285714)
    where(XLAND < 1.5)     u_urban(:,kde,:) = u_phy(:,1,:)
    where(XLAND < 1.5)     v_urban(:,kde,:) = v_phy(:,1,:)
    height_urban = 0.5*(dz_urban(1,kde,1) + dz_urban(1,kde-1,1))

!      IF(ITIME == 0) THEN  ! BEP or BEM urban models
    DO I = KDE-1, KDS, -1
      where(XLAND < 1.5)     p_urban(:,i,:) = p8w(:,1,:)*exp(9.8*height_urban/287.04/t_phy(:,1,:))
      where(XLAND < 1.5)   rho_urban(:,i,:) = p_urban(:,i,:)/287.04/t_phy(:,1,:)
      where(XLAND < 1.5) theta_urban(:,i,:) = t_phy(:,1,:)*(100000.0/p_urban(:,i,:))**(0.285714)
      where(XLAND < 1.5)     u_urban(:,i,:) = u_urban(:,kde,:)*log(zlvl-height_urban)/log(zlvl)
      where(XLAND < 1.5)     v_urban(:,i,:) = v_urban(:,kde,:)*log(zlvl-height_urban)/log(zlvl)
      height_urban = height_urban + urban_atmosphere_thickness
    END DO
!      ENDIF
  ENDIF

!------------------------------------------------------------------------
! Noah-MP updates we can do before spatial loop.
!------------------------------------------------------------------------

! create a few fields that are IN in WRF - coszen, julian,yr

 DO J = YSTART,YEND
 DO I = XSTART,XEND
   IF(SF_URBAN_PHYSICS == 0) THEN
     CALL CALC_DECLIN(OLDDATE(1:19),XLAT(I,J), XLONG(I,J),COSZEN(I,J),JULIAN)
   ELSE
     CALL CALC_DECLIN(OLDDATE(1:19),XLAT(I,J), XLONG(I,J),COSZEN(I,J),JULIAN, &
                      HRANG(I,J), DECLIN, GMT, JULDAY)
   ENDIF
   IF(IOPT_IRR > 0 ) THEN
      CALL LOCAL_TIME(OLDDATE(1:19), XLONG(I,J), LOCTIM(I,J))
   ENDIF
 END DO
 END DO

 READ(OLDDATE(1:4),*)  YR
 YEARLEN = 365                      ! find length of year for phenology (also S Hemisphere)
 if (mod(YR,4) == 0) then
    YEARLEN = 366
    if (mod(YR,100) == 0) then
       YEARLEN = 365
       if (mod(YR,400) == 0) then
          YEARLEN = 366
       endif
    endif
 endif

 IF (ITIME == 1 .AND. .NOT. RESTART_FLAG ) THEN
   EAHXY = (P8W(:,1,:)*QV_CURR(:,1,:))/(0.622+QV_CURR(:,1,:)) ! Initial guess only.
   TAHXY = T_PHY(:,1,:)                                       ! Initial guess only.
   CHXY = 0.1
   CMXY = 0.1
 ENDIF

!------------------------------------------------------------------------
! Skip model call at t=1 since initial conditions are at start time; First model time is +1
!------------------------------------------------------------------------

IF (ITIME > 0) THEN

!------------------------------------------------------------------------
! Call to Noah-MP driver same as surface_driver
!------------------------------------------------------------------------
  sflx_count_sum = 0 ! Timing

! Timing information for SFLX:
 call system_clock(count=count_before_sflx, count_rate=clock_rate)
 
 
 CALL noahmplsm_mosaic_hue(ITIMESTEP,        YR,   JULIAN,   COSZEN,XLAT,XLONG, & ! IN : Time/Space-related
                       DZ8W,       DTBL,       DZS,    NUM_SOIL_LAYERS,       DX,            & ! IN : Model configuration
                       IVGTYP,   ISLTYP,    VEGFRA,   GVFMAX,      TMN,            & ! IN : Vegetation/Soil characteristics
 		               XLAND,     XICE, XICE_THRESHOLD,  CROPCAT,                      & ! IN : Vegetation/Soil characteristics
                       PLANTING,  HARVEST,SEASON_GDD,                               &
                       IDVEG, IOPT_CRS,  IOPT_BTR, IOPT_RUN, IOPT_SFC, IOPT_FRZ,  & ! IN : User options
                       IOPT_INF, IOPT_RAD,  IOPT_ALB, IOPT_SNF,IOPT_TBOT, IOPT_STC,  & ! IN : User options
                       IOPT_GLA, IOPT_RSF, IOPT_SOIL,IOPT_PEDO,IOPT_CROP, IOPT_IRR,  & ! IN : User options
                       IOPT_IRRM, IOPT_INFDV, IOPT_TDRN, IOPT_MOSAIC, IOPT_HUE, soiltstep,  &
                       IZ0TLND, SF_URBAN_PHYSICS,                                    & ! IN : User options
                       SOILCOMP,  SOILCL1,  SOILCL2,   SOILCL3,  SOILCL4,            & ! IN : User options
                       T_PHY,     QV_CURR,     U_PHY,    V_PHY,   SWDOWN,      SWDDIR,  & ! IN : Forcing
                       SWDDIF,   GLW,                                         &
                       P8W, RAINBL,        SR,                                & ! IN : Forcing
                       IRFRACT, SIFRACT,   MIFRACT,  FIFRACT,                      & ! IN : Irrigation fractions
                       TSK,      HFX,      QFX,        LH,   GRDFLX,    SMSTAV, & ! IN/OUT LSM eqv
                       SMSTOT,SFCRUNOFF, UDRUNOFF,    ALBEDO,    SNOWC,     SMOIS, & ! IN/OUT LSM eqv
                       SH2O,     TSLB,     SNOW,     SNOWH,   CANWAT,    ACSNOM, & ! IN/OUT LSM eqv
                       ACSNOW,    EMISS,     QSFC,                                 & ! IN/OUT LSM eqv
                       Z0,      ZNT,                                           & ! IN/OUT LSM eqv
                       IRNUMSI,  IRNUMMI,  IRNUMFI,   IRWATSI,  IRWATMI,    IRWATFI, & ! IN/OUT Irrigation
                       IRELOSS,  IRSIVOL,  IRMIVOL,   IRFIVOL,  IRRSPLH,   LLANDUSE, & ! IN/OUT Irrigation
                       ISNOWXY,     TVXY,     TGXY,  CANICEXY, CANLIQXY,     EAHXY, & ! IN/OUT Noah MP only
                       TAHXY,     CMXY,     CHXY,    FWETXY, SNEQVOXY,  ALBOLDXY, & ! IN/OUT Noah MP only
                       QSNOWXY, QRAINXY, WSLAKEXY,    ZWTXY,      WAXY,     WTXY,    TSNOXY, & ! IN/OUT Noah MP only
                       ZSNSOXY,  SNICEXY,  SNLIQXY,  LFMASSXY, RTMASSXY,  STMASSXY, & ! IN/OUT Noah MP only
                       WOODXY, STBLCPXY, FASTCPXY,    LAI,   XSAIXY,   TAUSSXY, & ! IN/OUT Noah MP only
                       SMOISEQ, SMCWTDXY,DEEPRECHXY,   RECHXY,  GRAINXY,    GDDXY,PGSXY,  & ! IN/OUT Noah MP only
                       GECROS_STATE,                                                & ! IN/OUT gecros model
                       QTDRAIN, TD_FRACTION,                                  &
                       T2MVXY,   T2MBXY,    Q2MVXY,   Q2MBXY,                      & ! OUT Noah MP only
                       TRADXY,    NEEXY,    GPPXY,     NPPXY,   FVEGXY,   RUNSFXY, & ! OUT Noah MP only !added by Aaron A.
                       RUNSBXY,   ECANXY,   EDIRXY,   ETRANXY,    FSAXY,    FIRAXY, & ! OUT Noah MP only
                       APARXY,    PSNXY,    SAVXY,     SAGXY,  RSSUNXY,   RSSHAXY, & ! OUT Noah MP only
                       BGAPXY,   WGAPXY,    TGVXY,     TGBXY,    CHVXY,     CHBXY, & ! OUT Noah MP only
                       SHGXY,    SHCXY,    SHBXY,     EVGXY,    EVBXY,     GHVXY, & ! OUT Noah MP only
                       GHBXY,    IRGXY,    IRCXY,     IRBXY,     TRXY,     EVCXY, & ! OUT Noah MP only
                       CHLEAFXY,   CHUCXY,   CHV2XY,    CHB2XY, RS,                  & ! OUT Noah MP only
                       QINTSXY   ,QINTRXY   ,QDRIPSXY   ,&
                       QDRIPRXY  ,QTHROSXY  ,QTHRORXY   ,&
                       QSNSUBXY  ,QSNFROXY  ,QSUBCXY    ,&
                       QFROCXY   ,QEVACXY   ,QDEWCXY    ,QFRZCXY   ,QMELTCXY   ,& 
                       QSNBOTXY  ,QMELTXY   ,PONDINGXY  ,PAHXY     ,PAHGXY, PAHVXY, PAHBXY,&
                       FPICEXY,RAINLSM,SNOWLSM,FORCTLSM ,FORCQLSM,FORCPLSM,FORCZLSM,FORCWLSM,&
                       ACC_SSOILXY, ACC_QINSURXY, ACC_QSEVAXY, ACC_ETRANIXY, EFLXBXY, &
                       SOILENERGY, SNOWENERGY, CANHSXY, &
                       ACC_DWATERXY, ACC_PRCPXY, ACC_ECANXY, ACC_ETRANXY, ACC_EDIRXY, &
#ifdef WRF_HYDRO
                       sfcheadrt,INFXSRT,soldrain,qtiledrain,ZWATBLE2D,             &    !O
#endif
                       ids,ide,  jds,jde,  kds,kde,                    &
                       ims,ime,  jms,jme,  kms,kme,                    &
                       its,ite,  jts,jte,  kts,kte,                    &
                       LANDUSEF, landusef2, number_land_use_catagories,                    &      ! Added by Aaron A. **IMPORTANT FOR THE FOR LOOPS OF FRACTIONAL LAND USE
                       mosaic_cat_index, number_mosaic_catagories,                           &      ! Added by Aaron A. **IMPORTANT FOR THE FOR LOOPS OF FRACTIONAL LAND USE
                       TSK_mosaic, HFX_mosaic, QFX_mosaic, LH_mosaic,                        &      ! Added by Aaron A. IN/OUT LSM
                       GRDFLX_mosaic, SFCRUNOFF_mosaic, UDRUNOFF_mosaic,                     &      ! Added by Aaorn A. IN/OUT LSM
                       ALBEDO_mosaic, SNOWC_mosaic, TSLB_mosaic, SMOIS_mosaic,               &      ! Added by Aaron A. IN/OUT LSM
                       SH2O_mosaic,  CANWAT_mosaic, SNOW_mosaic, SNOWH_mosaic,               &      ! Added by Aaron A. IN/OUT LSM
                       ACSNOM_mosaic, ACSNOW_mosaic, EMISS_mosaic, QSFC_mosaic,              &      ! Added by Aaron A. IN/OUT LSM
                       Z0_mosaic, ZNT_mosaic, RS_mosaic,                                      &      ! Added by Aaron A. IN/OUT LSM
                       IRNUMSI_mosaic, IRNUMMI_mosaic, IRNUMFI_mosaic, IRWATSI_mosaic,       &      ! Added by Aaron A. IN/OUT Irrigation Scheme
                       IRWATMI_mosaic, IRWATFI_mosaic, IRELOSS_mosaic, IRSIVOL_mosaic,       &      ! Added by Aaron A. IN/OUT Irrigation Scheme
                       IRMIVOL_mosaic, IRFIVOL_mosaic, IRRSPLH_mosaic,                       &      ! Added by Aaron A. IN/OUT Irrigation Scheme
                       isnowxy_mosaic, tvxy_mosaic, tgxy_mosaic, canicexy_mosaic,            &      ! Added by Aaron A. IN/OUT NOAH MP ONLY
                       canliqxy_mosaic, eahxy_mosaic, tahxy_mosaic,                          &      ! Added by Aaron A. IN/OUT NOAH MP ONLY
                       cmxy_mosaic, chxy_mosaic, fwetxy_mosaic, sneqvoxy_mosaic,             &      ! Added by Aaron A. IN/OUT NOAH MP ONLY
                       alboldxy_mosaic, qsnowxy_mosaic, qrainxy_mosaic, wslakexy_mosaic, zwtxy_mosaic,       &      ! Added by Aaron A. IN/OUT NOAH MP ONLY
                       waxy_mosaic, wtxy_mosaic, tsnoxy_mosaic, zsnsoxy_mosaic,              &      ! Added by Aaron A. IN/OUT NOAH MP ONLY
                       snicexy_mosaic, snliqxy_mosaic, lfmassxy_mosaic, rtmassxy_mosaic,     &      ! Added by Aaron A. IN/OUT NOAH MP ONLY
                       stmassxy_mosaic, woodxy_mosaic, stblcpxy_mosaic, fastcpxy_mosaic,     &      ! Added by Aaron A. IN/OUT NOAH MP ONLY
                       xsaixy_mosaic, xlai_mosaic,                                           &      ! Added by Aaron A. IN/OUT NOAH MP ONLY
                       smoiseq_mosaic, smcwtdxy_mosaic,                                      &      ! Added by Aaron A. IN/OUT NOAH MP ONLY
                       deeprechxy_mosaic, rechxy_mosaic, taussxy_mosaic,                     &      ! Added by Aaron A. IN/OUT NOAH MP ONLY
                       ACC_SSOILXY_mosaic, ACC_QINSURXY_mosaic, ACC_QSEVAXY_mosaic,          &
                       ACC_ETRANIXY_mosaic,  &
                       ACC_DWATERXY_mosaic, ACC_PRCPXY_mosaic,               &
                       ACC_ECANXY_mosaic, ACC_EDIRXY_mosaic, ACC_ETRANXY_mosaic,             &
                       0,                                                              &
                       TR_URB2D_mosaic,TB_URB2D_mosaic,                                      & !H urban  Aaron A. Mosaic
                       TG_URB2D_mosaic,TC_URB2D_mosaic,                                      & !H urban  Aaron A. Mosaic
                       QC_URB2D_mosaic,UC_URB2D_mosaic,                                      & !H urban  Aaron A. Mosaic
                       TRL_URB3D_mosaic,TBL_URB3D_mosaic,                                    & !H urban  Aaron A. Mosaic
                       TGL_URB3D_mosaic,                                                     & !H urban  Aaron A. Mosaic
                       SH_URB2D_mosaic,LH_URB2D_mosaic,                                      & !H urban  Aaron A. Mosaic
                       G_URB2D_mosaic,RN_URB2D_mosaic,                                       & !H urban  Aaron A. Mosaic
                       TS_URB2D_mosaic,CMR_SFCDIF_mosaic, CHR_SFCDIF_mosaic,                 & !H urban  Aaron A. Mosaic
                       CMC_SFCDIF_mosaic, CHC_SFCDIF_mosaic, CMGR_SFCDIF_mosaic,             &
                       CHGR_SFCDIF_mosaic, XXXR_URB2D_mosaic, XXXB_URB2D_mosaic,  XXXC_URB2D_mosaic,             &
                       XXXG_URB2D_mosaic, CMCR_URB2D_mosaic, TGR_URB2D_mosaic,               &
                       TGRL_URB3D_mosaic, SMR_URB3D_mosaic, DRELR_URB2D_mosaic,              &
                       DRELB_URB2D_mosaic, DRELG_URB2D_mosaic, FLXHUMR_URB2D_mosaic,         &
                       FLXHUMB_URB2D_mosaic, FLXHUMG_URB2D_mosaic,                           &
                       cmr_sfcdif,     chr_sfcdif,     cmc_sfcdif,                           &
                       chc_sfcdif,    cmgr_sfcdif,    chgr_sfcdif,                           &
                       tr_urb2d,       tb_urb2d,       tg_urb2d,                             & !H urban
                       tc_urb2d,       qc_urb2d,       uc_urb2d,                             & !H urban
                       xxxr_urb2d,     xxxb_urb2d,     xxxg_urb2d, xxxc_urb2d,               & !H urban
                       trl_urb3d,      tbl_urb3d,      tgl_urb3d,                            & !H urban
                       sh_urb2d,       lh_urb2d,        g_urb2d,   rn_urb2d,  ts_urb2d,      & !H urban
                       psim_urb2d,     psih_urb2d,      u10_urb2d,  v10_urb2d,               & !O urban
                       GZ1OZ0_urb2d,     AKMS_URB2D,                                         & !O urban
                       th2_urb2d,       q2_urb2d,      ust_urb2d,                            & !O urban
                       declin,          hrang,                                           & !I urban
                       num_roof_layers,      num_wall_layers,      num_road_layers,          & !I urban
                       dzr,            dzb,            dzg,                                  & !I urban
                       cmcr_urb2d,      tgr_urb2d,     tgrl_urb3d,  smr_urb3d,               & !H urban
                       drelr_urb2d,    drelb_urb2d,    drelg_urb2d,                          & !H urban
                       flxhumr_urb2d,  flxhumb_urb2d,  flxhumg_urb2d,                        & !H urban
                       julday,             yr,                                            &
                       frc_urb2d,    utype_urb2d,                                            & !I urban
                       chs,           chs2,           cqs2,                                  & !H
                       lb_urb2d,      hgt_urb2d,  lp_urb2d,                                  & !H multi-layer urban
                       mh_urb2d,     stdh_urb2d,       lf_urb2d,                             & !SLUCM
                       theta_urban,      rho_urban,    p_urban,        ust,                  & !I multi-layer urban
                       gmt,                                                                  & !I multi-layer urban MODIFIED BY AARON A.
                       RUNONSFXY,RUNONSFXY_mosaic,DETENTION_STORAGEXY,DETENTION_STORAGEXY_mosaic, VOL_FLUX_RUNONXY_mosaic, VOL_FLUX_SMXY_mosaic )
                !MP_RAINC =  RAINCV, MP_RAINNC =    RAINNCV, MP_SHCV = RAINSHV,&
                !MP_SNOW  = SNOWNCV, MP_GRAUP  = GRAUPELNCV, MP_HAIL = HAILNCV )


    call system_clock(count=count_after_sflx, count_rate=clock_rate)
    sflx_count_sum = sflx_count_sum + ( count_after_sflx - count_before_sflx )

ENDIF ! End If skipping first time step

ENDIF ! End mosaic if else swap

! Output for history
     OUTPUT_FOR_HISTORY: if (output_timestep > 0) then
        if (mod(ITIME*noah_timestep, output_timestep) == 0 .and. .not. skip_first_output ) then

           call prepare_output_file (trim(outdir), version, &
                igrid, output_timestep, llanduse, split_output_count, hgrid,                &
                ixfull, jxfull, ixpar, jxpar, xstartpar, ystartpar,                         &
                iswater, mapproj, lat1, lon1, dx, dy, truelat1, truelat2, cen_lon,          &
                nsoil, nsnow, dzs, startdate, olddate, spinup_loop, spinup_loops, IVGTYP, ISLTYP,number_mosaic_catagories)

           DEFINE_MODE_LOOP : do imode = 1, 2

              call set_output_define_mode(imode)

              ! For 3D arrays, we need to know whether the Z dimension is snow layers, or soil layers.
        !! Added by Aaron A. Mosaic Capabilities built off of
        ! code infrastructure.
        ! MOS1 = number of mosaic catagories
        ! MOS2 = number of mosaic catagories * number of soil layers
        ! MOS3 = number of mosaic catagories * number of snow layers
        ! Properties - Assigned or predicted
              call add_to_output(IVGTYP     , "IVGTYP"  , "Dominant vegetation category"         , "category"              )
              call add_to_output(ISLTYP     , "ISLTYP"  , "Dominant soil category"               , "category"              )
              call add_to_output(FVEGXY     , "FVEG"    , "Green Vegetation Fraction"              , "-"                   )
              call add_to_output(LAI        , "LAI"     , "Leaf area index"                      , "-"                     )
              call add_to_output(XSAIXY     , "SAI"     , "Stem area index"                      , "-"                     )
        ! Forcing
              !call add_to_output(SWDOWN     , "SWFORC"  , "Shortwave forcing"                    , "W m{-2}"               )
              !call add_to_output(COSZEN     , "COSZ"    , "Cosine of zenith angle"                    , "W m{-2}"               )
              !call add_to_output(GLW        , "LWFORC"  , "Longwave forcing"                    , "W m{-2}"               )
              !call add_to_output(RAINBL_tmp , "RAINRATE", "Precipitation rate"                   , "mm s{-1}"        )
        ! Grid energy budget terms
              !call add_to_output(EMISS      , "EMISS"   , "Grid emissivity"                    , ""               )
              !call add_to_output(FSAXY      , "FSA"     , "Total absorbed SW radiation"          , "W m{-2}"               )
              !call add_to_output(FIRAXY     , "FIRA"    , "Total net LW radiation to atmosphere" , "W m{-2}"               )
              call add_to_output(GRDFLX     , "GRDFLX"  , "Heat flux into the soil"              , "W m{-2}"               )
              call add_to_output(HFX        , "HFX"     , "Total sensible heat to atmosphere"    , "W m{-2}"               )
              call add_to_output(LH         , "LH"      , "Total latent heat to atmosphere"      , "W m{-2}"               )
              call add_to_output(ECANXY     , "ECAN"    , "Canopy water evaporation rate"        , "kg m{-2} s{-1}"        )
              call add_to_output(ETRANXY    , "ETRAN"   , "Transpiration rate"                   , "kg m{-2} s{-1}"        )
              call add_to_output(EDIRXY     , "EDIR"    , "Direct from soil evaporation rate"    , "kg m{-2} s{-1}"        )
              call add_to_output(ALBEDO     , "ALBEDO"  , "Surface albedo"                       , "-"                     )
        ! Grid water budget terms - in addition to above
              call add_to_output(UDRUNOFF   , "UGDRNOFF", "Accumulated underground runoff"       , "mm"                    )
              call add_to_output(SFCRUNOFF  , "SFCRNOFF", "Accumulatetd surface runoff"          , "mm"                    )
              call add_to_output(CANLIQXY   , "CANLIQ"  , "Canopy liquid water content"          , "mm"                    )
              call add_to_output(CANICEXY   , "CANICE"  , "Canopy ice water content"             , "mm"                    )
              !call add_to_output(ZWTXY      , "ZWT"     , "Depth to water table"                 , "m"                     )
              !call add_to_output(WAXY       , "WA"      , "Water in aquifer"                     , "kg m{-2}"              )
              !call add_to_output(WTXY       , "WT"      , "Water in aquifer and saturated soil"  , "kg m{-2}"              )
              !call add_to_output(QTDRAIN    , "QTDRAIN" , "Accumulated tile drainage"            , "mm"                    )
        ! Additional needed to close the canopy energy budget
              !call add_to_output(SAVXY      , "SAV"     , "Solar radiative heat flux absorbed by vegetation", "W m{-2}"    )
              !call add_to_output(TRXY       , "TR"      , "Transpiration heat"                     , "W m{-2}"             )
              !call add_to_output(EVCXY      , "EVC"     , "Canopy evap heat"                       , "W m{-2}"             )
              !call add_to_output(IRCXY      , "IRC"     , "Canopy net LW rad"                      , "W m{-2}"             )
              !call add_to_output(SHCXY      , "SHC"     , "Canopy sensible heat"                   , "W m{-2}"             )
        ! Additional needed to close the under canopy ground energy budget
              !call add_to_output(IRGXY      , "IRG"     , "Ground net LW rad"                      , "W m{-2}"             )
              !call add_to_output(SHGXY      , "SHG"     , "Ground sensible heat"                   , "W m{-2}"             )
              !call add_to_output(EVGXY      , "EVG"     , "Ground evap heat"                       , "W m{-2}"             )
              !call add_to_output(GHVXY      , "GHV"     , "Ground heat flux + to soil vegetated"   , "W m{-2}"             )
        ! Needed to close the bare ground energy budget
              !call add_to_output(SAGXY      , "SAG"     , "Solar radiative heat flux absorbed by ground", "W m{-2}"        )
              !call add_to_output(IRBXY      , "IRB"     , "Net LW rad to atm bare"                 , "W m{-2}"             )
              !call add_to_output(SHBXY      , "SHB"     , "Sensible heat to atm bare"              , "W m{-2}"             )
              !call add_to_output(EVBXY      , "EVB"     , "Evaporation heat to atm bare"           , "W m{-2}"             )
              !call add_to_output(GHBXY      , "GHB"     , "Ground heat flux + to soil bare"        , "W m{-2}"             )
        ! Above-soil temperatures
              call add_to_output(TRADXY     , "TRAD"    , "Surface radiative temperature"        , "K"                     )
              call add_to_output(TGXY       , "TG"      , "Ground temperature"                   , "K"                     )
              call add_to_output(TVXY       , "TV"      , "Vegetation temperature"               , "K"                     )
              call add_to_output(TAHXY      , "TAH"     , "Canopy air temperature"               , "K"                     )
              call add_to_output(TGVXY      , "TGV"     , "Ground surface Temp vegetated"          , "K"                   )
              call add_to_output(TGBXY      , "TGB"     , "Ground surface Temp bare"               , "K"                   )
              !call add_to_output(T2MVXY     , "T2MV"    , "2m Air Temp vegetated"                  , "K"                   )
              !call add_to_output(T2MBXY     , "T2MB"    , "2m Air Temp bare"                       , "K"                   )
        ! Above-soil moisture
              !call add_to_output(Q2MVXY     , "Q2MV"    , "2m mixing ratio vegetated"              , "kg/kg"               )
              !call add_to_output(Q2MBXY     , "Q2MB"    , "2m mixing ratio bare"                   , "kg/kg"               )
              call add_to_output(EAHXY      , "EAH"     , "Canopy air vapor pressure"            , "Pa"                    )
              call add_to_output(FWETXY     , "FWET"    , "Wetted or snowed fraction of canopy"  , "fraction"              )
        ! Snow and soil - 3D terms
              !call add_to_output(ZSNSOXY(:,-nsnow+1:0,:),  "ZSNSO_SN" , "Snow layer depths from snow surface", "m", "SNOW")
              !call add_to_output(SNICEXY    , "SNICE"   , "Snow layer ice"                       , "mm"             , "SNOW")
              !call add_to_output(SNLIQXY    , "SNLIQ"   , "Snow layer liquid water"              , "mm"             , "SNOW")
              call add_to_output(TSLB       , "SOIL_T"  , "soil temperature"                     , "K"              , "SOIL")
              call add_to_output(SMOIS      , "SOIL_M"  , "volumetric soil moisture"             , "m{3} m{-3}"     , "SOIL")
              call add_to_output(SH2O       , "SOIL_W"  , "liquid volumetric soil moisture"      , "m3 m-3"         , "SOIL")
              call add_to_output(TSNOXY     , "SNOW_T"  , "snow temperature"                     , "K"              , "SNOW")
        ! Snow - 2D terms
              call add_to_output(SNOWH      , "SNOWH"   , "Snow depth"                           , "m"                     )
              call add_to_output(SNOW       , "SNEQV"   , "Snow water equivalent"                , "kg m{-2}"              )
              !call add_to_output(QSNOWXY    , "QSNOW"   , "Snowfall rate on the ground"          , "mm s{-1}"              )
              call add_to_output(QRAINXY    , "QRAIN"   , "Rainfall rate on the ground"          , "mm s{-1}"              )
              !call add_to_output(ISNOWXY    , "ISNOW"   , "Number of snow layers"                , "count"                 )
              call add_to_output(SNOWC      , "FSNO"    , "Snow-cover fraction on the ground"      , ""                    )
              !call add_to_output(ACSNOW     , "ACSNOW"  , "accumulated snow fall"                  , "mm"                  )
              !call add_to_output(ACSNOM     , "ACSNOM"  , "accumulated melting water out of snow bottom" , "mm"            )
        ! Exchange coefficients
              !call add_to_output(CMXY       , "CM"      , "Momentum drag coefficient"            , ""                      )
              !call add_to_output(CHXY       , "CH"      , "Sensible heat exchange coefficient"   , ""                      )
              !call add_to_output(CHVXY      , "CHV"     , "Exchange coefficient vegetated"         , "m s{-1}"             )
              !call add_to_output(CHBXY      , "CHB"     , "Exchange coefficient bare"              , "m s{-1}"             )
              !call add_to_output(CHLEAFXY   , "CHLEAF"  , "Exchange coefficient leaf"              , "m s{-1}"             )
              !call add_to_output(CHUCXY     , "CHUC"    , "Exchange coefficient bare"              , "m s{-1}"             )
              !call add_to_output(CHV2XY     , "CHV2"    , "Exchange coefficient 2-meter vegetated" , "m s{-1}"             )
              !call add_to_output(CHB2XY     , "CHB2"    , "Exchange coefficient 2-meter bare"      , "m s{-1}"             )
        ! Carbon allocation model
              !call add_to_output(LFMASSXY   , "LFMASS"  , "Leaf mass"                            , "g m{-2}"               )
              !call add_to_output(RTMASSXY   , "RTMASS"  , "Mass of fine roots"                   , "g m{-2}"               )
              !call add_to_output(STMASSXY   , "STMASS"  , "Stem mass"                            , "g m{-2}"               )
              !call add_to_output(WOODXY     , "WOOD"    , "Mass of wood and woody roots"         , "g m{-2}"               )
              !call add_to_output(GRAINXY    , "GRAIN"   , "Mass of grain "                       , "g m{-2}"               ) !XING!THREE
              !call add_to_output(GDDXY      , "GDD"     , "Growing degree days(10) "             , ""                      ) !XING
              !call add_to_output(STBLCPXY   , "STBLCP"  , "Stable carbon in deep soil"           , "g m{-2}"               )
              !call add_to_output(FASTCPXY   , "FASTCP"  , "Short-lived carbon in shallow soil"   , "g m{-2}"               )
              !call add_to_output(NEEXY      , "NEE"     , "Net ecosystem exchange"                 , "g m{-2} s{-1} CO2"   )
              !call add_to_output(GPPXY      , "GPP"     , "Net instantaneous assimilation"         , "g m{-2} s{-1} C"     )
              !call add_to_output(NPPXY      , "NPP"     , "Net primary productivity"               , "g m{-2} s{-1} C"     )
              !call add_to_output(PSNXY      , "PSN"     , "Total photosynthesis"                   , "umol CO@ m{-2} s{-1}")
              !call add_to_output(APARXY     , "APAR"    , "Photosynthesis active energy by canopy" , "W m{-2}"             )
        IF (IOPT_MOSAIC.eq.1) then

              call add_to_output(TSK_mosaic , "TRAD_mos", "Surface Radiative Temperature mosaic",       "K", "MOS1"             )
              call add_to_output(TSLB_mosaic , "SOIL_T_mos", "soil temperature mosaic",       "K", "MOS2"             )
              call add_to_output(SMOIS_mosaic , "SOIL_M_mos", "Volumetric Soil Moisture mosaic",       "m{3} m{-3}", "MOS2"             )
              call add_to_output(SH2O_mosaic , "SOIL_W_mos", "Volumetric Soil Moisture mosaic",       "m{3} m{-3}", "MOS2"             )
              call add_to_output(CANLIQXY_mosaic , "CANLIQ_mos", "Canopy Liquid Water Content mosaic",       "mm", "MOS1"             )
              call add_to_output(SNOW_mosaic , "SNEQV_mos", "Snow water equivalent mosaic",       "kg m{-2}", "MOS1"             )
              call add_to_output(SNOWH_mosaic , "SNOWH_mos", "Snow depth mosaic",       "m", "MOS1"             )
              call add_to_output(SNOWC_mosaic , "FSNO_mos", "snow-cover-fraction on ground mosaic",       "", "MOS1"             )
              call add_to_output(LANDUSEF2(:,1:number_mosaic_catagories,:), "FRAC_LND", "Mosaic Land Use Weights, add up to 1", "[-]", "MOS1")
              call add_to_output(mosaic_cat_index(:,1:number_mosaic_catagories,:), "MOSAIC_CATS", "Mosaic Catagories for Fractional Land-Use", "[-]", "MOS1")
              call add_to_output(LH_mosaic, "LH_mos", "things","w ^-2", "MOS1")
        END IF ! end if mosaic
        IF (IOPT_HUE.eq.1) THEN
              call add_to_output(RUNONSFXY  , "SFCRNON", "Accumulatetd surface runon grid averaged"          , "mm"                    )
              call add_to_output(DETENTION_STORAGEXY  , "DETENTION_STORAGE", "Depth of water stored in Green Roof"          , "mm"                    )
        ENDIF
        ! additional NoahMP output
          if (noahmp_output > 0) then
              ! additional water budget terms
              call add_to_output(QINTSXY    , "QINTS"   , "canopy interception (loading) rate for snowfall", "mm s{-1}"    )
              call add_to_output(QINTRXY    , "QINTR"   , "canopy interception rate for rain"              , "mm s{-1}"    )
              call add_to_output(QDRIPSXY   , "QDRIPS"  , "drip (unloading) rate for intercepted snow"     , "mm s{-1}"    )
              call add_to_output(QDRIPRXY   , "QDRIPR"  , "drip rate for canopy intercepted rain"          , "mm s{-1}"    )
              call add_to_output(QTHROSXY   , "QTHROS"  , "throughfall of snowfall"                        , "mm s{-1}"    )
              call add_to_output(QTHRORXY   , "QTHROR"  , "throughfall for rain"                           , "mm s{-1}"    )
              call add_to_output(QSNSUBXY   , "QSNSUB"  , "snow surface sublimation rate"                  , "mm s{-1}"    )
              call add_to_output(QSNFROXY   , "QSNFRO"  , "snow surface frost rate"                        , "mm s{-1}"    )
              call add_to_output(QSUBCXY    , "QSUBC"   , "canopy snow sublimation rate"                   , "mm s{-1}"    )
              call add_to_output(QFROCXY    , "QFROC"   , "canopy snow frost rate"                         , "mm s{-1}"    )
              call add_to_output(QEVACXY    , "QEVAC"   , "canopy snow evaporation rate"                   , "mm s{-1}"    )
              call add_to_output(QDEWCXY    , "QDEWC"   , "canopy snow dew rate"                           , "mm s{-1}"    )
              call add_to_output(QFRZCXY    , "QFRZC"   , "refreezing rate of canopy liquid water"         , "mm s{-1}"    )
              call add_to_output(QMELTCXY   , "QMELTC"  , "melting rate of canopy snow"                    , "mm s{-1}"    )
              call add_to_output(QSNBOTXY   , "QSNBOT"  , "total liquid water (melt+rain through) out of snow bottom", "mm s{-1}"    )
              call add_to_output(QMELTXY    , "QMELT"   , "snow melt due to phase change"                  , "mm s{-1}"    )
              call add_to_output(PONDINGXY  , "PONDING" , "total surface ponding per time step"            , "mm s{-1}"    )
              call add_to_output(FPICEXY    , "FPICE"   , "snow fraction in precipitation"                 , ""            )
              call add_to_output(ACC_QINSURXY,"ACC_QINSUR","accumuated water flux to soil within soil timestep","mm s{-1}" )
              call add_to_output(ACC_QSEVAXY, "ACC_QSEVA" ,"accumulated soil surface evap rate within soil timestep","mm s{-1}")
              call add_to_output(ACC_ETRANIXY,"ACC_ETRANI","accumualted transpiration rate within soil timestep", "mm s{-1}", "SOIL")
              call add_to_output(ACC_DWATERXY,"ACC_DWATER","accumulated snow,soil,canopy water change within soil timestep", "mm")
              call add_to_output(ACC_PRCPXY , "ACC_PRCP", "accumulated precipitation within soil timestep", "mm")
              call add_to_output(ACC_ECANXY , "ACC_ECAN", "accumulated net canopy evaporation within soil timestep", "mm")
              call add_to_output(ACC_ETRANXY, "ACC_ETRAN","accumulated transpiration within soil timestep", "mm")
              call add_to_output(ACC_EDIRXY , "ACC_EDIR", "accumulated net soil/snow evaporation within soil timestep", "mm")
              ! additional energy terms
              call add_to_output(PAHXY      , "PAH"     , "Precipitation advected heat flux"                       , "W m{-2}" )
              call add_to_output(PAHGXY     , "PAHG"    , "Precipitation advected heat flux to below-canopy ground", "W m{-2}" )
              call add_to_output(PAHBXY     , "PAHB"    , "Precipitation advected heat flux to bare ground"        , "W m{-2}" )
              call add_to_output(PAHVXY     , "PAHV"    , "Precipitation advected heat flux to canopy"             , "W m{-2}" )
              call add_to_output(ACC_SSOILXY, "ACC_SSOIL","accumulated heat flux into snow and soil layers within soil timestep","W m{-2}" )
              call add_to_output(EFLXBXY    , "EFLXB"   , "accumulated heat flux through soil bottom"              , "J m{-2}" )
              call add_to_output(SOILENERGY , "SOILENERGY","energy content in soil relative to 273.16"             , "KJ m{-2}")
              call add_to_output(SNOWENERGY , "SNOWENERGY","energy content in snow relative to 273.16"             , "KJ m{-2}")
              call add_to_output(CANHSXY    , "CANHS"   , "canopy heat storage change"                             , "W m{-2}" )
              ! additional forcing terms
              call add_to_output(RAINLSM    , "RAINLSM" , "lowest model liquid precipitation into LSM"             , "mm s{-1}")
              call add_to_output(SNOWLSM    , "SNOWLSM" , "lowest model snowfall into LSM"                         , "mm s{-1}")
              call add_to_output(FORCTLSM   , "FORCTLSM", "lowest model T into LSM"                                , "K"       )
              call add_to_output(FORCQLSM   , "FORCQLSM", "lowest model Q into LSM"                                , "kg kg{-1}")
              call add_to_output(FORCPLSM   , "FORCPLSM", "lowest model P into LSM"                                , "Pa"      )
              call add_to_output(FORCZLSM   , "FORCZLSM", "lowest model Z into LSM"                                , "m"       )
              call add_to_output(FORCWLSM   , "FORCWLSM", "lowest model wind speed into LSM"                       , "m s{-1}" )
          endif
        ! Irrigation
            IF (irrigation_option > 0) THEN
              call add_to_output(IRNUMSI    , "IRNUMSI" , "Sprinkler irrigation count"             , "-"                   )
              call add_to_output(IRNUMMI    , "IRNUMMI" , "Micro irrigation count"                 , "-"                   )
              call add_to_output(IRNUMFI    , "IRNUMFI" , "Flood irrigation count"                 , "-"                   )
              call add_to_output(IRELOSS    , "IRELOSS" , "Sprinkler Evaporation"                  , "mm"                  )
              call add_to_output(IRSIVOL    , "IRSIVOL" , "Sprinkler irrigation amount"            , "mm"                  )
              call add_to_output(IRMIVOL    , "IRMIVOL" , "Micro irrigation amount"                , "mm"                  )
              call add_to_output(IRFIVOL    , "IRFIVOL" , "Flood irrigation amount"                , "mm"                  )
              call add_to_output(IRRSPLH    , "IRRSPLH" , "Latent heating due to sprinkler"        , "W m{-2}"             )
            ENDIF
        ! MMF groundwater  model
	    IF(RUNOFF_OPTION == 5) THEN
              call add_to_output(SMCWTDXY   , "SMCWTD"   , "soil water content between bottom of the soil and water table", "m3 m{-3}")
              call add_to_output(RECHXY     , "RECH"     , "recharge to or from the water table when shallow"             , "m"       )
              call add_to_output(DEEPRECHXY , "DEEPRECH" , "recharge to or from the water table when deep"                , "m"       )
              call add_to_output(QRFSXY     , "QRFS"     , "accumulated groundwater baselow"                              , "mm"      )
              call add_to_output(QRFXY      , "QRF"      , "groundwater baseflow"                                         , "m"       )
              call add_to_output(QSPRINGSXY , "QSPRINGS" , "accumulated seeping water"                                    , "mm"      )
              call add_to_output(QSPRINGXY  , "QSPRING"  , "instantaneous seeping water"                                  , "m"       )
              call add_to_output(QSLATXY    , "QSLAT"    , "accumulated lateral flow"                                     , "mm"      )
              call add_to_output(QLATXY     , "QLAT"     , "instantaneous lateral flow"                                   , "m"       )
	    ENDIF
	! For now, no urban output variables included

           enddo DEFINE_MODE_LOOP

           call finalize_output_file(split_output_count,itime)

	endif

	if(skip_first_output) skip_first_output = .false.

     endif OUTPUT_FOR_HISTORY

     if (IVGTYP(xstart,ystart)==ISWATER) then
       write(*,'(" ***DATE=", A19)', advance="NO") olddate
     else
       write(*,'(" ***DATE=", A19, 6F10.5)', advance="NO") olddate, TSLB(xstart,1,ystart), LAI(xstart,ystart)
     endif

!------------------------------------------------------------------------
! Write Restart - timestamp equal to output will have same states
!------------------------------------------------------------------------

      if ( (restart_frequency_hours .gt. 0) .and. &
           (mod(ITIME, int(restart_frequency_hours*3600./nint(dtbl))) == 0)) then
       call lsm_restart()
      else
       if (restart_frequency_hours <= 0) then
          if ( (olddate( 9:10) == "01") .and. (olddate(12:13) == "00") .and. &
               (olddate(15:16) == "00") .and. (olddate(18:19) == "00") ) then
                call lsm_restart()  ! jlm - i moved all the restart code to a subroutine.
          endif
       endif
      endif

!------------------------------------------------------------------------
! Advance the time
!------------------------------------------------------------------------

     call geth_newdate(newdate, olddate, nint(dtbl))
     olddate = newdate

     if(itime > 0 .and. spinup_loops > 0 .and. mod(itime,ntime/(spinup_loops+1)) == 0) then
       spinup_loop = spinup_loop + 1
       call geth_newdate(olddate, startdate, nint(dtbl))
       reset_spinup_date = .true.
     end if

! update the timer
     call system_clock(count=clock_count_2, count_rate=clock_rate)
     timing_sum = timing_sum + float(clock_count_2-clock_count_1)/float(clock_rate)
     write(*,'("    Timing: ",f6.2," Cumulative:  ", f10.2, "  SFLX: ", f6.2 )') &
          float(clock_count_2-clock_count_1)/float(clock_rate), &
          timing_sum, real(sflx_count_sum) / real(clock_rate)
     clock_count_1 = clock_count_2

#ifdef WRF_HYDRO
      call hrldas_drv_HYDRO(TSLB(:,1:NSOIL,:),SMOIS(:,1:NSOIL,:),SH2O(:,1:NSOIL,:),infxsrt,&
                            sfcheadrt,soldrain,qtiledrain,ZWATBLE2D,ix,jx,NSOIL)
#endif

end subroutine land_driver_exe

!!===============================================================================
subroutine lsm_restart()
  implicit none

  print*, 'Write restart at '//olddate(1:13)

  call prepare_restart_file (trim(outdir), version, igrid, llanduse, olddate, startdate, &
       ixfull, jxfull, ixpar, jxpar, xstartpar, ystartpar,                               &
       nsoil, nsnow, num_urban_hi, dx, dy, truelat1, truelat2, mapproj, lat1, lon1,   &
       cen_lon, iswater, ivgtyp, number_mosaic_catagories)

  call add_to_restart(TSLB      , "SOIL_T", layers="SOIL")
  call add_to_restart(TSNOXY    , "SNOW_T", layers="SNOW")
  call add_to_restart(SMOIS     , "SMC"   , layers="SOIL")
  call add_to_restart(SH2O      , "SH2O"  , layers="SOIL")
  call add_to_restart(ZSNSOXY   , "ZSNSO" , layers="SOSN")
  call add_to_restart(SNICEXY   , "SNICE" , layers="SNOW")
  call add_to_restart(SNLIQXY   , "SNLIQ" , layers="SNOW")
  call add_to_restart(QSNOWXY   , "QSNOW" )
  call add_to_restart(QRAINXY   , "QRAIN" )
  call add_to_restart(FWETXY    , "FWET"  )
  call add_to_restart(SNEQVOXY  , "SNEQVO")
  call add_to_restart(EAHXY     , "EAH"   )
  call add_to_restart(TAHXY     , "TAH"   )
  call add_to_restart(ALBOLDXY  , "ALBOLD")
  call add_to_restart(CMXY      , "CM"    )
  call add_to_restart(CHXY      , "CH"    )
  call add_to_restart(ISNOWXY   , "ISNOW" )
  call add_to_restart(CANLIQXY  , "CANLIQ")
  call add_to_restart(CANICEXY  , "CANICE")
  call add_to_restart(SNOW      , "SNEQV" )
  call add_to_restart(SNOWH     , "SNOWH" )
  call add_to_restart(TVXY      , "TV"    )
  call add_to_restart(TGXY      , "TG"    )
  call add_to_restart(ZWTXY     , "ZWT"   )
  call add_to_restart(WAXY      , "WA"    )
  call add_to_restart(WTXY      , "WT"    )
  call add_to_restart(WSLAKEXY  , "WSLAKE")
  call add_to_restart(LFMASSXY  , "LFMASS")
  call add_to_restart(RTMASSXY  , "RTMASS")
  call add_to_restart(STMASSXY  , "STMASS")
  call add_to_restart(WOODXY    , "WOOD"  )
  call add_to_restart(GRAINXY   , "GRAIN" )
  call add_to_restart(GDDXY     , "GDD"   )
  call add_to_restart(STBLCPXY  , "STBLCP")
  call add_to_restart(FASTCPXY  , "FASTCP")
  call add_to_restart(LAI       , "LAI"   )
  call add_to_restart(XSAIXY    , "SAI"   )
  call add_to_restart(VEGFRA    , "VEGFRA")
  call add_to_restart(GVFMIN    , "GVFMIN")
  call add_to_restart(GVFMAX    , "GVFMAX")
  call add_to_restart(ACSNOM    , "ACMELT")
  call add_to_restart(ACSNOW    , "ACSNOW")
  call add_to_restart(TAUSSXY   , "TAUSS" )
  call add_to_restart(QSFC      , "QSFC"  )
  call add_to_restart(SFCRUNOFF , "SFCRUNOFF")
  call add_to_restart(UDRUNOFF  , "UDRUNOFF" )
  call add_to_restart(QTDRAIN   , "QTDRAIN"  )
  call add_to_restart(ACC_SSOILXY ,"ACC_SSOIL" )
  call add_to_restart(ACC_QINSURXY,"ACC_QINSUR")
  call add_to_restart(ACC_QSEVAXY ,"ACC_QSEVA" )
  call add_to_restart(ACC_ETRANIXY,"ACC_ETRANI", layers="SOIL")
  call add_to_restart(ACC_DWATERXY,"ACC_DWATER")
  call add_to_restart(ACC_PRCPXY  ,"ACC_PRCP"  )
  call add_to_restart(ACC_ECANXY  ,"ACC_ECAN"  )
  call add_to_restart(ACC_ETRANXY ,"ACC_ETRAN" )
  call add_to_restart(ACC_EDIRXY  ,"ACC_EDIR"  )
 ! irrigation scheme
  IF (irrigation_option > 0) THEN
     call add_to_restart(IRNUMSI   , "IRNUMSI")
     call add_to_restart(IRNUMMI   , "IRNUMMI")
     call add_to_restart(IRNUMFI   , "IRNUMFI")
     call add_to_restart(IRWATSI   , "IRWATSI")
     call add_to_restart(IRWATMI   , "IRWATMI")
     call add_to_restart(IRWATFI   , "IRWATFI")
     call add_to_restart(IRSIVOL   , "IRSIVOL")
     call add_to_restart(IRMIVOL   , "IRMIVOL")
     call add_to_restart(IRFIVOL   , "IRFIVOL")
  ENDIF
  ! below for opt_run = 5
  call add_to_restart(SMOISEQ   , "SMOISEQ"  , layers="SOIL"  )
  call add_to_restart(AREAXY    , "AREAXY"     )
  call add_to_restart(SMCWTDXY  , "SMCWTDXY"   )
  call add_to_restart(DEEPRECHXY, "DEEPRECHXY" )
  call add_to_restart(QSLATXY   , "QSLATXY"    )
  call add_to_restart(QRFSXY    , "QRFSXY"     )
  call add_to_restart(QSPRINGSXY, "QSPRINGSXY" )
  call add_to_restart(RECHXY    , "RECHXY"     )
  call add_to_restart(QRFXY     , "QRFXY"      )
  call add_to_restart(QSPRINGXY , "QSPRINGXY"  )
  call add_to_restart(FDEPTHXY , "FDEPTHXY"  )
  call add_to_restart(RIVERCONDXY , "RIVERCONDXY"  )
  call add_to_restart(RIVERBEDXY , "RIVERBEDXY"  )
  call add_to_restart(EQZWT , "EQZWT"  )
  call add_to_restart(PEXPXY , "PEXPXY"  )

  IF(SF_URBAN_PHYSICS > 0 )  THEN  ! any urban model

      call add_to_restart(     SH_URB2D  ,      "SH_URB2D" )
      call add_to_restart(     LH_URB2D  ,      "LH_URB2D" )
      call add_to_restart(      G_URB2D  ,       "G_URB2D" )
      call add_to_restart(     RN_URB2D  ,      "RN_URB2D" )
      call add_to_restart(     TS_URB2D  ,      "TS_URB2D" )
      call add_to_restart(    FRC_URB2D  ,     "FRC_URB2D" )
      call add_to_restart(  UTYPE_URB2D  ,   "UTYPE_URB2D" )
      call add_to_restart(     LP_URB2D  ,      "LP_URB2D" )
      call add_to_restart(     LB_URB2D  ,      "LB_URB2D" )
      call add_to_restart(    HGT_URB2D  ,     "HGT_URB2D" )
      call add_to_restart(     MH_URB2D  ,      "MH_URB2D" )
      call add_to_restart(   STDH_URB2D  ,    "STDH_URB2D" )
      call add_to_restart(     HI_URB2D  ,      "HI_URB2D", layers="URBN")
      call add_to_restart(     LF_URB2D  ,      "LF_URB2D", layers="URBN")

  ENDIF

  IF(SF_URBAN_PHYSICS == 1 ) THEN  ! single layer urban model

      call add_to_restart(    CMR_SFCDIF ,     "CMR_SFCDIF" )
      call add_to_restart(    CHR_SFCDIF ,     "CHR_SFCDIF" )
      call add_to_restart(    CMC_SFCDIF ,     "CMC_SFCDIF" )
      call add_to_restart(    CHC_SFCDIF ,     "CHC_SFCDIF" )
      call add_to_restart(   CMGR_SFCDIF ,    "CMGR_SFCDIF" )
      call add_to_restart(   CHGR_SFCDIF ,    "CHGR_SFCDIF" )
      call add_to_restart(     TR_URB2D  ,      "TR_URB2D" )
      call add_to_restart(     TB_URB2D  ,      "TB_URB2D" )
      call add_to_restart(     TG_URB2D  ,      "TG_URB2D" )
      call add_to_restart(     TC_URB2D  ,      "TC_URB2D" )
      call add_to_restart(     QC_URB2D  ,      "QC_URB2D" )
      call add_to_restart(     UC_URB2D  ,      "UC_URB2D" )
      call add_to_restart(   XXXR_URB2D  ,    "XXXR_URB2D" )
      call add_to_restart(   XXXB_URB2D  ,    "XXXB_URB2D" )
      call add_to_restart(   XXXG_URB2D  ,    "XXXG_URB2D" )
      call add_to_restart(   XXXC_URB2D  ,    "XXXC_URB2D" )
      call add_to_restart(    TRL_URB3D  ,     "TRL_URB3D", layers="SOIL" )
      call add_to_restart(    TBL_URB3D  ,     "TBL_URB3D", layers="SOIL" )
      call add_to_restart(    TGL_URB3D  ,     "TGL_URB3D", layers="SOIL" )
      call add_to_restart(   CMCR_URB2D  ,    "CMCR_URB2D" )
      call add_to_restart(    TGR_URB2D  ,     "TGR_URB2D" )
      call add_to_restart(   TGRL_URB3D  ,    "TGRL_URB3D", layers="SOIL" )
      call add_to_restart(    SMR_URB3D  ,     "SMR_URB3D", layers="SOIL" )
      call add_to_restart(  DRELR_URB2D  ,   "DRELR_URB2D" )
      call add_to_restart(  DRELB_URB2D  ,   "DRELB_URB2D" )
      call add_to_restart(  DRELG_URB2D  ,   "DRELG_URB2D" )
      call add_to_restart(FLXHUMR_URB2D  , "FLXHUMR_URB2D" )
      call add_to_restart(FLXHUMB_URB2D  , "FLXHUMB_URB2D" )
      call add_to_restart(FLXHUMG_URB2D  , "FLXHUMG_URB2D" )

  ENDIF

  IF(SF_URBAN_PHYSICS == 2 .or. SF_URBAN_PHYSICS == 3) THEN  ! BEP or BEM urban models

      call add_to_restart(    TRB_URB4D  ,     "TRB_URB4D", layers="URBN" )
      call add_to_restart(    TW1_URB4D  ,     "TW1_URB4D", layers="URBN" )
      call add_to_restart(    TW2_URB4D  ,     "TW2_URB4D", layers="URBN" )
      call add_to_restart(    TGB_URB4D  ,     "TGB_URB4D", layers="URBN" )
      call add_to_restart(   SFW1_URB3D  ,    "SFW1_URB3D", layers="URBN" )
      call add_to_restart(   SFW2_URB3D  ,    "SFW2_URB3D", layers="URBN" )
      call add_to_restart(    SFR_URB3D  ,     "SFR_URB3D", layers="URBN" )
      call add_to_restart(    SFG_URB3D  ,     "SFG_URB3D", layers="URBN" )

  ENDIF

  IF(SF_URBAN_PHYSICS == 3) THEN  ! BEM urban model

      call add_to_restart(   TLEV_URB3D  ,    "TLEV_URB3D", layers="URBN" )
      call add_to_restart(   QLEV_URB3D  ,    "QLEV_URB3D", layers="URBN" )
      call add_to_restart( TW1LEV_URB3D  ,  "TW1LEV_URB3D", layers="URBN" )
      call add_to_restart( TW2LEV_URB3D  ,  "TW2LEV_URB3D", layers="URBN" )
      call add_to_restart(  TGLEV_URB3D  ,   "TGLEV_URB3D", layers="URBN" )
      call add_to_restart(  TFLEV_URB3D  ,   "TFLEV_URB3D", layers="URBN" )
      call add_to_restart(  SF_AC_URB3D  ,   "SF_AC_URB3D" )
      call add_to_restart(  LF_AC_URB3D  ,   "LF_AC_URB3D" )
      call add_to_restart(  CM_AC_URB3D  ,   "CM_AC_URB3D" )
      call add_to_restart( SFVENT_URB3D  ,  "SFVENT_URB3D" )
      call add_to_restart( LFVENT_URB3D  ,  "LFVENT_URB3D" )
      call add_to_restart( SFWIN1_URB3D  ,  "SFWIN1_URB3D", layers="URBN" )
      call add_to_restart( SFWIN2_URB3D  ,  "SFWIN2_URB3D", layers="URBN" )
      call add_to_restart(  EP_PV_URB3D  ,   "EP_PV_URB3D" )
      call add_to_restart(   T_PV_URB3D  ,    "T_PV_URB3D", layers="URBN" )
      call add_to_restart(    TRV_URB4D  ,    "TRV_URB4D" , layers="URBN" )
      call add_to_restart(    QR_URB4D   ,     "QR_URB4D" , layers="URBN" )
      call add_to_restart(   QGR_URB3D   ,    "QGR_URB3D"  )
      call add_to_restart(   TGR_URB3D   ,    "TGR_URB3D"  )
      call add_to_restart(   DRAIN_URB4D ,   "DRAIN_URB4D", layers="URBN" )
      call add_to_restart( DRAINGR_URB3D , "DRAINGR_URB3D" )
      call add_to_restart(    SFRV_URB3D ,    "SFRV_URB3D", layers="URBN" )
      call add_to_restart(    LFRV_URB3D ,    "LFRV_URB3D", layers="URBN" )
      call add_to_restart(     DGR_URB3D ,     "DGR_URB3D", layers="URBN" )
      call add_to_restart(      DG_URB3D ,      "DG_URB3D", layers="URBN" )
      call add_to_restart(     LFR_URB3D ,     "LFR_URB3D", layers="URBN" )
      call add_to_restart(     LFG_URB3D ,     "LFG_URB3D", layers="URBN" )

  ENDIF
  ! MOS1 = only catagories of mosaic
  ! MOS2 = 4 soil catagoris for each mosiac
  ! MOS3 = number of soil layers times mosiac
  ! MOS4 = 7 possible layers (snow and soil) for each mosaic
  IF(IOPT_MOSAIC.eq.1) THEN
    call add_to_restart(        TSLB_mosaic,        "TSLB_mosaic" , layers="MOS2" )
    call add_to_restart(     tsnoxy_mosaic ,  "tsnoxy_mosaic"   , layers="MOS3")
    call add_to_restart(      SMOIS_mosaic ,       "SMOIS_mosaic", layers="MOS2" )
    call add_to_restart(        SH2O_mosaic ,       "SH2O_mosaic" , layers="MOS2" )
    call add_to_restart(   zsnsoxy_mosaic,   "zsnsoxy_mosaic"   , layers="MOS4" )
    call add_to_restart(   snicexy_mosaic,  "snicexy_mosaic"    , layers="MOS3")
    call add_to_restart(   snliqxy_mosaic,  "snliqxy_mosaic"   , layers="MOS3")
    call add_to_restart(    qsnowxy_mosaic , "qsnowxy_mosaic"   , layers="MOS1")
    call add_to_restart(    qrainxy_mosaic , "qrainxy_mosaic"    , layers="MOS1")
    call add_to_restart(        fwetxy_mosaic , "fwetxy_mosaic"  , layers="MOS1")
    call add_to_restart(   sneqvoxy_mosaic, "sneqvoxy_mosaic"   , layers="MOS1")
    call add_to_restart(        eahxy_mosaic , "eahxy_mosaic"   , layers="MOS1")
    call add_to_restart(     tahxy_mosaic , "tahxy_mosaic"   , layers="MOS1")
    call add_to_restart(   alboldxy_mosaic ,"alboldxy_mosaic"  , layers="MOS1" )
    call add_to_restart(        cmxy_mosaic,   "cmxy_mosaic"    , layers="MOS1")
    call add_to_restart(      chxy_mosaic ,    "chxy_mosaic"    , layers="MOS1")
    call add_to_restart(   isnowxy_mosaic,       "isnowxy_mosaic"  , layers="MOS1")
    call add_to_restart(   canliqxy_mosaic , "canliqxy_mosaic" , layers="MOS1")
    call add_to_restart(    canicexy_mosaic ,"canicexy_mosaic"  , layers="MOS1")
    call add_to_restart(         SNOW_mosaic ,      "SNOW_mosaic" , layers="MOS1" )
    call add_to_restart(       SNOWH_mosaic ,      "SNOWH_mosaic"  , layers="MOS1")
    call add_to_restart(           tvxy_mosaic,     "tvxy_mosaic"  , layers="MOS1")
    call add_to_restart(         tgxy_mosaic,       "tgxy_mosaic"  , layers="MOS1")
    call add_to_restart(      zwtxy_mosaic,   "zwtxy_mosaic"    , layers="MOS1")
    call add_to_restart(          waxy_mosaic, "waxy_mosaic"    , layers="MOS1")
    call add_to_restart(       wtxy_mosaic ,   "wtxy_mosaic"   , layers="MOS1")
    call add_to_restart(  wslakexy_mosaic , "wslakexy_mosaic"   , layers="MOS1") 
    call add_to_restart(  lfmassxy_mosaic, "lfmassxy_mosaic"   , layers="MOS1")
    call add_to_restart( rtmassxy_mosaic,  "rtmassxy_mosaic"   , layers="MOS1" )
    call add_to_restart(  stmassxy_mosaic, "stmassxy_mosaic"    , layers="MOS1")
    call add_to_restart(    woodxy_mosaic,   "woodxy_mosaic"    , layers="MOS1")
    call add_to_restart(   grainxy_mosaic,  "grainxy_mosaic"    , layers="MOS1")
    call add_to_restart(        gddxy_mosaic,  "gddxy_mosaic"   , layers="MOS1")
    call add_to_restart(  stblcpxy_mosaic, "stblcpxy_mosaic"    , layers="MOS1")
    call add_to_restart(  fastcpxy_mosaic , "fastcpxy_mosaic"   , layers="MOS1")
    call add_to_restart(    xsaixy_mosaic ,  "xsaixy_mosaic"    , layers="MOS1")
    call add_to_restart(       xlai_mosaic,    "xlai_mosaic"    , layers="MOS1")
    call add_to_restart(       VEGFRA_mosaic,    "VEGFRA_mosaic"    , layers="MOS1")    
    call add_to_restart(   ACSNOM_mosaic, "ACSNOM_mosaic"    , layers="MOS1" )
    call add_to_restart(   ACSNOW_mosaic, "ACSNOW_mosaic"    , layers="MOS1" )
    call add_to_restart(   taussxy_mosaic, "taussxy_mosaic"  , layers="MOS1" )
    call add_to_restart(  QSFC_mosaic,  "QSFC_mosaic"         , layers="MOS1")
    call add_to_restart( SFCRUNOFF_mosaic,  "SFCRUNOFF_mosaic", layers="MOS1")
    call add_to_restart(   UDRUNOFF_mosaic, "UDRUNOFF_mosaic" , layers="MOS1")
    call add_to_restart(  ACC_SSOILXY_mosaic, "ACC_SSOILXY_mosaic" , layers="MOS1" )
    call add_to_restart(  ACC_QINSURXY_mosaic, "ACC_QINSURXY_mosaic", layers="MOS1")
    call add_to_restart(  ACC_QSEVAXY_mosaic, "ACC_QSEVAXY_mosaic"  , layers="MOS1")
    call add_to_restart(  ACC_ETRANIXY_mosaic, "ACC_ETRANIXY_mosaic", layers="MOS2")
    call add_to_restart(  ACC_DWATERXY_mosaic, "ACC_DWATERXY_mosaic", layers="MOS1")
    call add_to_restart(  ACC_PRCPXY_mosaic,  "ACC_PRCPXY_mosaic" , layers="MOS1"   )
    call add_to_restart(  ACC_ECANXY_mosaic, "ACC_ECANXY_mosaic" , layers="MOS1"  )
    call add_to_restart(  ACC_ETRANXY_mosaic, "ACC_ETRANXY_mosaic" , layers="MOS1" )
    call add_to_restart(  ACC_EDIRXY_mosaic , "ACC_EDIRXY_mosaic" , layers="MOS1"  )
    ! unsure if we need this, but adding anyway
    call add_to_restart(    CANWAT_mosaic  ,      "CANWAT_mosaic" , layers="MOS1" )
    call add_to_restart(    TSK_mosaic ,       "TSK_mosaic", layers="MOS1"  )

    call add_to_restart(  SMOISEQ_mosaic,  "smoiseq_mosaic"   , layers="MOS2"    )
    call add_to_restart(  smcwtdxy_mosaic,  "smcwtdxy_mosaic"   , layers="MOS1"  )
    call add_to_restart(  deeprechxy_mosaic, "deeprechxy_mosaic" , layers="MOS1")

    IF (irrigation_option > 0) THEN
      call add_to_restart(    IRNUMSI_mosaic,  "IRNUMSI_mosaic" , layers="MOS1"  )
      call add_to_restart(    IRNUMMI_mosaic,  "IRNUMMI_mosaic" , layers="MOS1"   )
      call add_to_restart(    IRNUMFI_mosaic, "IRNUMFI_mosaic"  , layers="MOS1"  )
      call add_to_restart(    IRWATSI_mosaic, "IRWATSI_mosaic"  , layers="MOS1"  )
      call add_to_restart(    IRWATMI_mosaic , "IRWATMI_mosaic" , layers="MOS1"  )
      call add_to_restart(    IRWATFI_mosaic,   "IRWATFI_mosaic"  , layers="MOS1"  )
      call add_to_restart(   IRELOSS_mosaic , "IRELOSS_mosaic"   , layers="MOS1" )
      call add_to_restart(   IRSIVOL_mosaic, "IRSIVOL_mosaic"  , layers="MOS1"  )
      call add_to_restart(   IRMIVOL_mosaic ,  "IRMIVOL_mosaic" , layers="MOS1"  )
      call add_to_restart(    IRFIVOL_mosaic,  "IRFIVOL_mosaic" , layers="MOS1"  )
      call add_to_restart(     IRRSPLH_mosaic, "IRRSPLH_mosaic"  , layers="MOS1"  )
    ENDIF


    IF(SF_URBAN_PHYSICS > 0 )  THEN  ! any urban model

      call add_to_restart(        SH_URB2D_mosaic,  "SH_URB2D_mosaic" , layers="MOS1"   )
      call add_to_restart(        LH_URB2D_mosaic,  "LH_URB2D_mosaic"  , layers="MOS1"  )
      call add_to_restart(        G_URB2D_mosaic ,  "G_URB2D_mosaic" , layers="MOS1" )
      call add_to_restart(        RN_URB2D_mosaic,  "RN_URB2D_mosaic"  , layers="MOS1" )
      call add_to_restart(        TS_URB2D_mosaic,  "TS_URB2D_mosaic" , layers="MOS1"   )
      
    ENDIF

    IF(SF_URBAN_PHYSICS == 1 ) THEN  ! single layer urban model

      call add_to_restart(         TR_URB2D_mosaic,  "TR_URB2D_mosaic" , layers="MOS1"   )
      call add_to_restart(         TB_URB2D_mosaic , "TB_URB2D_mosaic" , layers="MOS1"  )
      call add_to_restart(         TG_URB2D_mosaic,  "TG_URB2D_mosaic"  , layers="MOS1" )
      call add_to_restart(         TC_URB2D_mosaic,  "TC_URB2D_mosaic"  , layers="MOS1"  )
      call add_to_restart(         QC_URB2D_mosaic,  "QC_URB2D_mosaic"  , layers="MOS1" )
      call add_to_restart(         UC_URB2D_mosaic,  "UC_URB2D_mosaic"  , layers="MOS1"  )
      call add_to_restart(         TRL_URB3D_mosaic, "TRL_URB3D_mosaic"  , layers="MOS2"  )
      call add_to_restart(         TBL_URB3D_mosaic,  "TBL_URB3D_mosaic"  , layers="MOS2"  )
      call add_to_restart(         TGL_URB3D_mosaic, "TGL_URB3D_mosaic" , layers="MOS2"  )

      call add_to_restart(    CMR_SFCDIF_mosaic ,     "CMR_SFCDIF_mosaic" ,layers="MOS1" )
      call add_to_restart(    CHR_SFCDIF_mosaic ,     "CHR_SFCDIF_mosaic" ,layers="MOS1" )
      call add_to_restart(    CMC_SFCDIF_mosaic ,     "CMC_SFCDIF_mosaic" ,layers="MOS1" )
      call add_to_restart(    CHC_SFCDIF_mosaic ,     "CHC_SFCDIF_mosaic",layers="MOS1"  )
      call add_to_restart(   CMGR_SFCDIF_mosaic ,    "CMGR_SFCDIF_mosaic",layers="MOS1"  )
      call add_to_restart(   CHGR_SFCDIF_mosaic ,    "CHGR_SFCDIF_mosaic" ,layers="MOS1" )
      call add_to_restart(   XXXR_URB2D_mosaic  ,    "XXXR_URB2D_mosaic" ,layers="MOS1" )
      call add_to_restart(   XXXB_URB2D_mosaic  ,    "XXXB_URB2D_mosaic" ,layers="MOS1" )
      call add_to_restart(   XXXG_URB2D_mosaic  ,    "XXXG_URB2D_mosaic" ,layers="MOS1" )
      call add_to_restart(   XXXC_URB2D_mosaic  ,    "XXXC_URB2D_mosaic" ,layers="MOS1" )
      call add_to_restart(   CMCR_URB2D_mosaic  ,    "CMCR_URB2D_mosaic" ,layers="MOS1" )
      call add_to_restart(    TGR_URB2D_mosaic  ,     "TGR_URB2D_mosaic" ,layers="MOS1" )
      call add_to_restart(   TGRL_URB3D_mosaic  ,    "TGRL_URB3D_mosaic", layers="MOS2" )
      call add_to_restart(    SMR_URB3D_mosaic  ,     "SMR_URB3D_mosaic", layers="MOS2" )
      call add_to_restart(  DRELR_URB2D_mosaic  ,   "DRELR_URB2D_mosaic",layers="MOS1"  )
      call add_to_restart(  DRELB_URB2D_mosaic  ,   "DRELB_URB2D_mosaic" ,layers="MOS1" )
      call add_to_restart(  DRELG_URB2D_mosaic  ,   "DRELG_URB2D_mosaic" ,layers="MOS1" )
      call add_to_restart(FLXHUMR_URB2D_mosaic  , "FLXHUMR_URB2D_mosaic" ,layers="MOS1" )
      call add_to_restart(FLXHUMB_URB2D_mosaic  , "FLXHUMB_URB2D_mosaic" ,layers="MOS1" )
      call add_to_restart(FLXHUMG_URB2D_mosaic  , "FLXHUMG_URB2D_mosaic" ,layers="MOS1" )


    ENDIF

    IF (IOPT_HUE.eq.1) THEN
      call add_to_restart(   RUNONSFXY_mosaic, "RUNONSFXY_mosaic" , layers="MOS1")
      call add_to_restart(   DETENTION_STORAGEXY_mosaic, "DETENTION_STORAGEXY_mosaic" , layers="MOS1")
      call add_to_restart( VOL_FLUX_RUNONXY_mosaic,  "VOL_FLUX_RUNONXY_mosaic", layers="MOS1")
      call add_to_restart(   VOL_FLUX_SMXY_mosaic, "VOL_FLUX_SMXY_mosaic" , layers="MOS2")
    END IF


  ENDIF


  call finalize_restart_file()

end subroutine lsm_restart

!------------------------------------------------------------------------------------

SUBROUTINE CALC_DECLIN ( NOWDATE, LATITUDE, LONGITUDE, COSZ, JULIAN, &
                         HRANG_OUT, DECLIN_OUT, GMT_OUT, JULDAY_OUT )

  USE MODULE_DATE_UTILITIES
!---------------------------------------------------------------------
   IMPLICIT NONE
!---------------------------------------------------------------------

! !ARGUMENTS:
   CHARACTER(LEN=19), INTENT(IN)  :: NOWDATE    ! YYYY-MM-DD_HH:MM:SS
   REAL,              INTENT(IN)  :: LATITUDE
   REAL,              INTENT(IN)  :: LONGITUDE
   REAL,              INTENT(OUT) :: COSZ
   REAL,              INTENT(OUT) :: JULIAN
   REAL,    OPTIONAL, INTENT(OUT) :: HRANG_OUT
   REAL,    OPTIONAL, INTENT(OUT) :: DECLIN_OUT
   REAL,    OPTIONAL, INTENT(OUT) :: GMT_OUT
   INTEGER, OPTIONAL, INTENT(OUT) :: JULDAY_OUT
   REAL                           :: OBECL
   REAL                           :: SINOB
   REAL                           :: SXLONG
   REAL                           :: ARG
   REAL                           :: TLOCTIM
   REAL                           :: HRANG
   REAL                           :: DECLIN
   REAL                           :: GMT
   INTEGER                        :: JULDAY
   INTEGER                        :: IDAY
   INTEGER                        :: IHOUR
   INTEGER                        :: IMINUTE
   INTEGER                        :: ISECOND

   REAL, PARAMETER :: DEGRAD = 3.14159265/180.
   REAL, PARAMETER :: DPD    = 360./365.

   CALL GETH_IDTS(NOWDATE(1:10), NOWDATE(1:4)//"-01-01", IDAY)
   READ(NOWDATE(12:13), *) IHOUR
   READ(NOWDATE(15:16), *) IMINUTE
   READ(NOWDATE(18:19), *) ISECOND
   GMT = REAL(IHOUR) + IMINUTE/60.0 + ISECOND/3600.0
   JULIAN = REAL(IDAY) + GMT/24.
   JULDAY = IDAY

!
! FOR SHORT WAVE RADIATION

   DECLIN=0.

!-----OBECL : OBLIQUITY = 23.5 DEGREE.

   OBECL=23.5*DEGRAD
   SINOB=SIN(OBECL)

!-----CALCULATE LONGITUDE OF THE SUN FROM VERNAL EQUINOX:

   IF(JULIAN.GE.80.)SXLONG=DPD*(JULIAN-80.)*DEGRAD
   IF(JULIAN.LT.80.)SXLONG=DPD*(JULIAN+285.)*DEGRAD
   ARG=SINOB*SIN(SXLONG)
   DECLIN=ASIN(ARG)

   TLOCTIM = REAL(IHOUR) + REAL(IMINUTE)/60.0 + REAL(ISECOND)/3600.0 + LONGITUDE/15.0 ! LOCAL TIME IN HOURS
   TLOCTIM = AMOD(TLOCTIM+24.0, 24.0)
   HRANG=15.*(TLOCTIM-12.)*DEGRAD
   COSZ=SIN(LATITUDE*DEGRAD)*SIN(DECLIN)+COS(LATITUDE*DEGRAD)*COS(DECLIN)*COS(HRANG)

   IF (PRESENT(HRANG_OUT ))  HRANG_OUT = HRANG
   IF (PRESENT(DECLIN_OUT)) DECLIN_OUT = DECLIN
   IF (PRESENT(GMT_OUT   ))    GMT_OUT = GMT
   IF (PRESENT(JULDAY_OUT)) JULDAY_OUT = JULDAY

 END SUBROUTINE CALC_DECLIN

!---------------------------------------------------------------------
 SUBROUTINE LOCAL_TIME(NOWDATE, LONGITUDE, TLOCTIM)
!---------------------------------------------------------------------
   IMPLICIT NONE
!---------------------------------------------------------------------
   CHARACTER(LEN=19), INTENT(IN)  :: NOWDATE    ! YYYY-MM-DD_HH:MM:SS
   REAL,              INTENT(IN)  :: LONGITUDE
   REAL,              INTENT(OUT) :: TLOCTIM
   INTEGER                        :: IHOUR
   INTEGER                        :: IMINUTE
   INTEGER                        :: ISECOND

   READ(NOWDATE(12:13), *) IHOUR
   READ(NOWDATE(15:16), *) IMINUTE
   READ(NOWDATE(18:19), *) ISECOND

   TLOCTIM = REAL(IHOUR) + REAL(IMINUTE)/60.0 + REAL(ISECOND)/3600.0 + LONGITUDE/15.0 ! LOCAL TIME IN HOURS
   TLOCTIM = AMOD(TLOCTIM+24.0, 24.0)

 END SUBROUTINE LOCAL_TIME


end module module_NoahMP_hrldas_driver

!subroutine wrf_message(msg)
!  implicit none
!  character(len=*), intent(in) :: msg
!  print*, msg
!end subroutine wrf_message

logical function wrf_dm_on_monitor() result(l)
  l = .TRUE.
  return
end function wrf_dm_on_monitor
