MODULE module_sf_noahmpdrv

!-------------------------------
#if ( WRF_CHEM == 1 )
  USE module_data_gocart_dust
#endif
!-------------------------------

!
CONTAINS
!
  SUBROUTINE noahmplsm(ITIMESTEP,        YR,   JULIAN,   COSZIN,XLAT,XLONG, & ! IN : Time/Space-related
                  DZ8W,       DT,       DZS,    NSOIL,       DX,            & ! IN : Model configuration
	        IVGTYP,   ISLTYP,    VEGFRA,   VEGMAX,      TMN,            & ! IN : Vegetation/Soil characteristics
		 XLAND,     XICE,XICE_THRES,  CROPCAT,                      & ! IN : Vegetation/Soil characteristics
	       PLANTING,  HARVEST,SEASON_GDD,                               &
                 IDVEG, IOPT_CRS,  IOPT_BTR, IOPT_RUN, IOPT_SFC, IOPT_FRZ,  & ! IN : User options
              IOPT_INF, IOPT_RAD,  IOPT_ALB, IOPT_SNF,IOPT_TBOT, IOPT_STC,  & ! IN : User options
              IOPT_GLA, IOPT_RSF, IOPT_SOIL,IOPT_PEDO,IOPT_CROP, IOPT_IRR,  & ! IN : User options
             IOPT_IRRM, IOPT_INFDV, IOPT_TDRN,soiltstep,IOPT_MOSAIC, IOPT_HUE, & ! IN : User options
              IZ0TLND, SF_URBAN_PHYSICS,                                    & ! IN : User options
	      SOILCOMP,  SOILCL1,  SOILCL2,   SOILCL3,  SOILCL4,            & ! IN : User options
                   T3D,     QV3D,     U_PHY,    V_PHY,   SWDOWN,     SWDDIR,&
                SWDDIF,      GLW,                                           & ! IN : Forcing
		 P8W3D,PRECIP_IN,        SR,                                & ! IN : Forcing
               IRFRACT,  SIFRACT,   MIFRACT,  FIFRACT,                      & ! IN : Noah MP only
                   TSK,      HFX,      QFX,        LH,   GRDFLX,    SMSTAV, & ! IN/OUT LSM eqv
                SMSTOT,SFCRUNOFF, UDRUNOFF,    ALBEDO,    SNOWC,     SMOIS, & ! IN/OUT LSM eqv
		  SH2O,     TSLB,     SNOW,     SNOWH,   CANWAT,    ACSNOM, & ! IN/OUT LSM eqv
		ACSNOW,    EMISS,     QSFC,                                 & ! IN/OUT LSM eqv
 		    Z0,      ZNT,                                           & ! IN/OUT LSM eqv
               IRNUMSI,  IRNUMMI,  IRNUMFI,   IRWATSI,  IRWATMI,   IRWATFI, & ! IN/OUT Noah MP only
               IRELOSS,  IRSIVOL,  IRMIVOL,   IRFIVOL,  IRRSPLH,  LLANDUSE, & ! IN/OUT Noah MP only
               ISNOWXY,     TVXY,     TGXY,  CANICEXY, CANLIQXY,     EAHXY, & ! IN/OUT Noah MP only
	         TAHXY,     CMXY,     CHXY,    FWETXY, SNEQVOXY,  ALBOLDXY, & ! IN/OUT Noah MP only
               QSNOWXY, QRAINXY,  WSLAKEXY, ZWTXY,  WAXY,  WTXY,    TSNOXY, & ! IN/OUT Noah MP only
	       ZSNSOXY,  SNICEXY,  SNLIQXY,  LFMASSXY, RTMASSXY,  STMASSXY, & ! IN/OUT Noah MP only
	        WOODXY, STBLCPXY, FASTCPXY,    XLAIXY,   XSAIXY,   TAUSSXY, & ! IN/OUT Noah MP only
	       SMOISEQ, SMCWTDXY,DEEPRECHXY,   RECHXY,  GRAINXY,    GDDXY,PGSXY,  & ! IN/OUT Noah MP only
               GECROS_STATE,                                                & ! IN/OUT gecros model
               QTDRAIN,   TD_FRACTION,                                      & ! IN/OUT tile drainage
	        T2MVXY,   T2MBXY,    Q2MVXY,   Q2MBXY,                      & ! OUT Noah MP only
	        TRADXY,    NEEXY,    GPPXY,     NPPXY,   FVEGXY,   RUNSFXY, & ! OUT Noah MP only
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
               QSNBOTXY  ,QMELTXY   ,PONDINGXY ,PAHXY      ,PAHGXY, PAHVXY, PAHBXY,&
               FPICEXY,RAINLSM,SNOWLSM,FORCTLSM ,FORCQLSM,FORCPLSM,FORCZLSM,FORCWLSM,&
               ACC_SSOILXY, ACC_QINSURXY, ACC_QSEVAXY, ACC_ETRANIXY, EFLXBXY, &
               SOILENERGY, SNOWENERGY, CANHSXY, &
               ACC_DWATERXY, ACC_PRCPXY, ACC_ECANXY, ACC_ETRANXY, ACC_EDIRXY, &
!                 BEXP_3D,SMCDRY_3D,SMCWLT_3D,SMCREF_3D,SMCMAX_3D,          & ! placeholders to activate 3D soil
!		 DKSAT_3D,DWSAT_3D,PSISAT_3D,QUARTZ_3D,                     &
!		 REFDK_2D,REFKDT_2D,                                        &
!                IRR_FRAC_2D,IRR_HAR_2D,IRR_LAI_2D,IRR_MAD_2D,FILOSS_2D,    &
!                SPRIR_RATE_2D,MICIR_RATE_2D,FIRTFAC_2D,IR_RAIN_2D,         &
!                BVIC_2D,AXAJ_2D,BXAJ_2D,XXAJ_2D,BDVIC_2D,GDVIC_2D,BBVIC_2D,&
!              KLAT_FAC,TDSMC_FAC,TD_DC,TD_DCOEF,TD_DDRAIN,TD_RADI,TD_SPAC, &
#ifdef WRF_HYDRO
               sfcheadrt,INFXSRT,soldrain,qtiledrain,ZWATBLE2D,             &
#endif
               ids,ide,  jds,jde,  kds,kde,                    &
               ims,ime,  jms,jme,  kms,kme,                    &
               its,ite,  jts,jte,  kts,kte,                    &
               MP_RAINC, MP_RAINNC, MP_SHCV, MP_SNOW, MP_GRAUP, MP_HAIL     )
!----------------------------------------------------------------
    USE MODULE_SF_NOAHMPLSM
!    USE MODULE_SF_NOAHMPLSM, only: noahmp_options, NOAHMP_SFLX, noahmp_parameters
    USE module_sf_noahmp_glacier
    USE NOAHMP_TABLES, ONLY: ISICE_TABLE, CO2_TABLE, O2_TABLE, DEFAULT_CROP_TABLE, ISCROP_TABLE, ISURBAN_TABLE, NATURAL_TABLE, &
                             LCZ_1_TABLE,LCZ_2_TABLE,LCZ_3_TABLE,LCZ_4_TABLE,LCZ_5_TABLE,LCZ_6_TABLE,LCZ_7_TABLE,LCZ_8_TABLE,  &
                             LCZ_9_TABLE,LCZ_10_TABLE,LCZ_11_TABLE

    USE module_sf_urban,    only: IRI_SCHEME
    USE module_ra_gfdleta,  only: cal_mon_day
!----------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------

! IN only

    INTEGER,                                         INTENT(IN   ) ::  ITIMESTEP ! timestep number
    INTEGER,                                         INTENT(IN   ) ::  YR        ! 4-digit year
    REAL,                                            INTENT(IN   ) ::  JULIAN    ! Julian day
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  COSZIN    ! cosine zenith angle
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  XLAT      ! latitude [rad]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  XLONG     ! latitude [rad]
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  DZ8W      ! thickness of atmo layers [m]
    REAL,                                            INTENT(IN   ) ::  DT        ! timestep [s]
    REAL,    DIMENSION(1:nsoil),                     INTENT(IN   ) ::  DZS       ! thickness of soil layers [m]
    INTEGER,                                         INTENT(IN   ) ::  NSOIL     ! number of soil layers
    REAL,                                            INTENT(IN   ) ::  DX        ! horizontal grid spacing [m]
    INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  IVGTYP    ! vegetation type
    INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  ISLTYP    ! soil type
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  VEGFRA    ! vegetation fraction []
    REAL,    DIMENSION( ims:ime ,         jms:jme ), INTENT(IN   ) ::  VEGMAX    ! annual max vegetation fraction []
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  TMN       ! deep soil temperature [K]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  XLAND     ! =2 ocean; =1 land/seaice
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  XICE      ! fraction of grid that is seaice
    REAL,                                            INTENT(IN   ) ::  XICE_THRES! fraction of grid determining seaice
    INTEGER,                                         INTENT(IN   ) ::  IDVEG     ! dynamic vegetation (1 -> off ; 2 -> on) with opt_crs = 1
    INTEGER,                                         INTENT(IN   ) ::  IOPT_CRS  ! canopy stomatal resistance (1-> Ball-Berry; 2->Jarvis)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_BTR  ! soil moisture factor for stomatal resistance (1-> Noah; 2-> CLM; 3-> SSiB)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_RUN  ! runoff and groundwater (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_SFC  ! surface layer drag coeff (CH & CM) (1->M-O; 2->Chen97)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_FRZ  ! supercooled liquid water (1-> NY06; 2->Koren99)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_INF  ! frozen soil permeability (1-> NY06; 2->Koren99)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_RAD  ! radiation transfer (1->gap=F(3D,cosz); 2->gap=0; 3->gap=1-Fveg)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_ALB  ! snow surface albedo (1->BATS; 2->CLASS)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_SNF  ! rainfall & snowfall (1-Jordan91; 2->BATS; 3->Noah)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_TBOT ! lower boundary of soil temperature (1->zero-flux; 2->Noah)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_STC  ! snow/soil temperature time scheme
    INTEGER,                                         INTENT(IN   ) ::  IOPT_GLA  ! glacier option (1->phase change; 2->simple)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_RSF  ! surface resistance (1->Sakaguchi/Zeng; 2->Seller; 3->mod Sellers; 4->1+snow)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_SOIL ! soil configuration option
    INTEGER,                                         INTENT(IN   ) ::  IOPT_PEDO ! soil pedotransfer function option
    INTEGER,                                         INTENT(IN   ) ::  IOPT_CROP ! crop model option (0->none; 1->Liu et al.; 2->Gecros)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_IRR  ! irrigation scheme (0->none; >1 irrigation scheme ON)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_IRRM ! irrigation method
    INTEGER,                                         INTENT(IN   ) ::  IOPT_INFDV! infiltration options for dynamic VIC infiltration (1->Philip; 2-> Green-Ampt;3->Smith-Parlange)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_TDRN ! tile drainage (0-> no tile drainage; 1-> simple tile drainage;2->Hooghoudt's)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_MOSAIC
    INTEGER,                                         INTENT(IN   ) ::  IOPT_HUE
    REAL,                                            INTENT(IN   ) ::  soiltstep ! soil timestep (s), default:0->same as main model timestep
    INTEGER,                                         INTENT(IN   ) ::  IZ0TLND   ! option of Chen adjustment of Czil (not used)
    INTEGER,                                         INTENT(IN   ) ::  sf_urban_physics   ! urban physics option
    REAL,    DIMENSION( ims:ime,       8, jms:jme ), INTENT(IN   ) ::  SOILCOMP  ! soil sand and clay percentage
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SOILCL1   ! soil texture in layer 1
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SOILCL2   ! soil texture in layer 2
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SOILCL3   ! soil texture in layer 3
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SOILCL4   ! soil texture in layer 4
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  T3D       ! 3D atmospheric temperature valid at mid-levels [K]
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  QV3D      ! 3D water vapor mixing ratio [kg/kg_dry]
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  U_PHY     ! 3D U wind component [m/s]
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  V_PHY     ! 3D V wind component [m/s]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SWDOWN    ! solar down at surface [W m-2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SWDDIF    ! solar down at surface [W m-2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SWDDIR    ! solar down at surface [W m-2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  GLW       ! longwave down at surface [W m-2]
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  P8W3D     ! 3D pressure, valid at interface [Pa]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  PRECIP_IN ! total input precipitation [mm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SR        ! frozen precipitation ratio [-]

!Optional Detailed Precipitation Partitioning Inputs
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ), OPTIONAL ::  MP_RAINC  ! convective precipitation entering land model [mm] ! MB/AN : v3.7
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ), OPTIONAL ::  MP_RAINNC ! large-scale precipitation entering land model [mm]! MB/AN : v3.7
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ), OPTIONAL ::  MP_SHCV   ! shallow conv precip entering land model [mm]      ! MB/AN : v3.7
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ), OPTIONAL ::  MP_SNOW   ! snow precipitation entering land model [mm]       ! MB/AN : v3.7
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ), OPTIONAL ::  MP_GRAUP  ! graupel precipitation entering land model [mm]    ! MB/AN : v3.7
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ), OPTIONAL ::  MP_HAIL   ! hail precipitation entering land model [mm]       ! MB/AN : v3.7

! Crop Model
    INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  CROPCAT   ! crop catagory
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  PLANTING  ! planting date
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  HARVEST   ! harvest date
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SEASON_GDD! growing season GDD
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  GRAINXY   ! mass of grain XING [g/m2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  GDDXY     ! growing degree days XING (based on 10C)
 INTEGER,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  PGSXY

! gecros model
    REAL,    DIMENSION( ims:ime,       60,jms:jme ), INTENT(INOUT) ::  gecros_state !  gecros crop

!Tile drain variables
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  QTDRAIN
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN)    ::  TD_FRACTION

#ifdef WRF_HYDRO
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  sfcheadrt,INFXSRT,soldrain,qtiledrain,ZWATBLE2D   ! for WRF-Hydro
#endif
! placeholders for 3D soil
!    REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(IN) ::  BEXP_3D   ! C-H B exponent
!    REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(IN) ::  SMCDRY_3D ! Soil Moisture Limit: Dry
!    REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(IN) ::  SMCWLT_3D ! Soil Moisture Limit: Wilt
!    REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(IN) ::  SMCREF_3D ! Soil Moisture Limit: Reference
!    REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(IN) ::  SMCMAX_3D ! Soil Moisture Limit: Max
!    REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(IN) ::  DKSAT_3D  ! Saturated Soil Conductivity
!    REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(IN) ::  DWSAT_3D  ! Saturated Soil Diffusivity
!    REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(IN) ::  PSISAT_3D ! Saturated Matric Potential
!    REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(IN) ::  QUARTZ_3D ! Soil quartz content
!    REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  REFDK_2D  ! Reference Soil Conductivity
!    REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  REFKDT_2D ! Soil Infiltration Parameter
!    REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  BVIC_2D   ! VIC model infiltration parameter [-] for opt_run=6
!    REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  AXAJ_2D   ! Xinanjiang: Tension water distribution inflection parameter [-] for opt_run=7
!    REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  BXAJ_2D   ! Xinanjiang: Tension water distribution shape parameter [-] for opt_run=7
!    REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  XXAJ_2D   ! Xinanjiang: Free water distribution shape parameter [-] for opt_run=7
!    REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  BDVIC_2D  ! VIC model infiltration parameter [-]
!    REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  GDVIC_2D  ! Mean Capillary Drive (m) for infiltration models
!    REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  BBVIC_2D  ! DVIC heterogeniety paramater [-]

! placeholders for 2D irrigation parameters
!    REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  IRR_FRAC_2D   ! irrigation Fraction
!    REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  IRR_HAR_2D    ! number of days before harvest date to stop irrigation
!    REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  IRR_LAI_2D    ! Minimum lai to trigger irrigation
!    REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  IRR_MAD_2D    ! management allowable deficit (0-1)
!    REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  FILOSS_2D     ! fraction of flood irrigation loss (0-1)
!    REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  SPRIR_RATE_2D ! mm/h, sprinkler irrigation rate
!    REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  MICIR_RATE_2D ! mm/h, micro irrigation rate
!    REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  FIRTFAC_2D    ! flood application rate factor
!    REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  IR_RAIN_2D    ! maximum precipitation to stop irrigation trigger

! placeholders for 2D tile drainage parameters
!    REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  KLAT_FAC   ! factor multiplier to hydraulic conductivity
!    REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  TDSMC_FAC  ! factor multiplier to field capacity
!    REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  TD_DC      ! drainage coefficient for simple
!    REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  TD_DCOEF   ! drainge coefficient for Hooghoudt
!    REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  TD_DDRAIN  ! depth of drain
!    REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  TD_RADI    ! tile radius
!    REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  TD_SPAC    ! tile spacing

! INOUT (with generic LSM equivalent)

    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  TSK       ! surface radiative temperature [K]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  HFX       ! sensible heat flux [W m-2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  QFX       ! latent heat flux [kg s-1 m-2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  LH        ! latent heat flux [W m-2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  GRDFLX    ! ground/snow heat flux [W m-2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SMSTAV    ! soil moisture avail. [not used]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SMSTOT    ! total soil water [mm][not used]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SFCRUNOFF ! accumulated surface runoff [m]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  UDRUNOFF  ! accumulated sub-surface runoff [m]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ALBEDO    ! total grid albedo []
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SNOWC     ! snow cover fraction []
    REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(INOUT) ::  SMOIS     ! volumetric soil moisture [m3/m3]
    REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(INOUT) ::  SH2O      ! volumetric liquid soil moisture [m3/m3]
    REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(INOUT) ::  TSLB      ! soil temperature [K]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SNOW      ! snow water equivalent [mm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SNOWH     ! physical snow depth [m]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  CANWAT    ! total canopy water + ice [mm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACSNOM    ! accumulated snow melt (mm)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACSNOW    ! accumulated snow on grid
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  EMISS     ! surface bulk emissivity
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  QSFC      ! bulk surface specific humidity
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  Z0        ! combined z0 sent to coupled model
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ZNT       ! combined z0 sent to coupled model
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  RS        ! Total stomatal resistance (s/m)

    INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ISNOWXY   ! actual no. of snow layers
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  TVXY      ! vegetation leaf temperature
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  TGXY      ! bulk ground surface temperature
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  CANICEXY  ! canopy-intercepted ice (mm)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  CANLIQXY  ! canopy-intercepted liquid water (mm)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  EAHXY     ! canopy air vapor pressure (pa)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  TAHXY     ! canopy air temperature (k)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  CMXY      ! bulk momentum drag coefficient
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  CHXY      ! bulk sensible heat exchange coefficient
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  FWETXY    ! wetted or snowed fraction of the canopy (-)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SNEQVOXY  ! snow mass at last time step(mm h2o)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ALBOLDXY  ! snow albedo at last time step (-)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  QSNOWXY   ! snowfall on the ground [mm/s]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  QRAINXY   ! rainfall on the ground [mm/s]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  WSLAKEXY  ! lake water storage [mm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ZWTXY     ! water table depth [m]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  WAXY      ! water in the "aquifer" [mm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  WTXY      ! groundwater storage [mm]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  TSNOXY    ! snow temperature [K]
    REAL,    DIMENSION( ims:ime,-2:NSOIL, jms:jme ), INTENT(INOUT) ::  ZSNSOXY   ! snow layer depth [m]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  SNICEXY   ! snow layer ice [mm]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  SNLIQXY   ! snow layer liquid water [mm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  LFMASSXY  ! leaf mass [g/m2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  RTMASSXY  ! mass of fine roots [g/m2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  STMASSXY  ! stem mass [g/m2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  WOODXY    ! mass of wood (incl. woody roots) [g/m2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  STBLCPXY  ! stable carbon in deep soil [g/m2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  FASTCPXY  ! short-lived carbon, shallow soil [g/m2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  XLAIXY    ! leaf area index
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  XSAIXY    ! stem area index
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  TAUSSXY   ! snow age factor
    REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(INOUT) ::  SMOISEQ   ! eq volumetric soil moisture [m3/m3]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SMCWTDXY  ! soil moisture content in the layer to the water table when deep
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  DEEPRECHXY ! recharge to the water table when deep
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  RECHXY    ! recharge to the water table (diagnostic)

! OUT (with no Noah LSM equivalent)

    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  T2MVXY    ! 2m temperature of vegetation part
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  T2MBXY    ! 2m temperature of bare ground part
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  Q2MVXY    ! 2m mixing ratio of vegetation part
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  Q2MBXY    ! 2m mixing ratio of bare ground part
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  TRADXY    ! surface radiative temperature (k)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  NEEXY     ! net ecosys exchange (g/m2/s CO2)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  GPPXY     ! gross primary assimilation [g/m2/s C]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  NPPXY     ! net primary productivity [g/m2/s C]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FVEGXY    ! Noah-MP vegetation fraction [-]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  RUNSFXY   ! surface runoff [mm] per soil timestep
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  RUNSBXY   ! subsurface runoff [mm] per soil timestep
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  ECANXY    ! evaporation of intercepted water (mm/s)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  EDIRXY    ! soil surface evaporation rate (mm/s]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  ETRANXY   ! transpiration rate (mm/s)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FSAXY     ! total absorbed solar radiation (w/m2)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FIRAXY    ! total net longwave rad (w/m2) [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  APARXY    ! photosyn active energy by canopy (w/m2)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  PSNXY     ! total photosynthesis (umol co2/m2/s) [+]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SAVXY     ! solar rad absorbed by veg. (w/m2)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SAGXY     ! solar rad absorbed by ground (w/m2)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  RSSUNXY   ! sunlit leaf stomatal resistance (s/m)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  RSSHAXY   ! shaded leaf stomatal resistance (s/m)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  BGAPXY    ! between gap fraction
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  WGAPXY    ! within gap fraction
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  TGVXY     ! under canopy ground temperature[K]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  TGBXY     ! bare ground temperature [K]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  CHVXY     ! sensible heat exchange coefficient vegetated
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  CHBXY     ! sensible heat exchange coefficient bare-ground
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SHGXY     ! veg ground sen. heat [w/m2]   [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SHCXY     ! canopy sen. heat [w/m2]   [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SHBXY     ! bare sensible heat [w/m2]     [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  EVGXY     ! veg ground evap. heat [w/m2]  [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  EVBXY     ! bare soil evaporation [w/m2]  [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  GHVXY     ! veg ground heat flux [w/m2]  [+ to soil]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  GHBXY     ! bare ground heat flux [w/m2] [+ to soil]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  IRGXY     ! veg ground net LW rad. [w/m2] [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  IRCXY     ! canopy net LW rad. [w/m2] [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  IRBXY     ! bare net longwave rad. [w/m2] [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  TRXY      ! transpiration [w/m2]  [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  EVCXY     ! canopy evaporation heat [w/m2]  [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  CHLEAFXY  ! leaf exchange coefficient
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  CHUCXY    ! under canopy exchange coefficient
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  CHV2XY    ! veg 2m exchange coefficient
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  CHB2XY    ! bare 2m exchange coefficient
! additional output variables
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  PAHXY     ! precipitation advected heat
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  PAHGXY    ! precipitation advected heat
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  PAHBXY    ! precipitation advected heat
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  PAHVXY    ! precipitation advected heat
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QINTSXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QINTRXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QDRIPSXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QDRIPRXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QTHROSXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QTHRORXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QSNSUBXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QSNFROXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QSUBCXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QFROCXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QEVACXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QDEWCXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QFRZCXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QMELTCXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QSNBOTXY  !total liquid water (snowmelt + rain through pack)out of snowpack bottom [mm/s]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QMELTXY   !snowmelt due to phase change (mm/s)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  PONDINGXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FPICEXY    !fraction of ice in precip
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  RAINLSM     !rain rate                   (mm/s)  AJN
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SNOWLSM     !liquid equivalent snow rate (mm/s)  AJN
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FORCTLSM
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FORCQLSM
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FORCPLSM
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FORCZLSM
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FORCWLSM
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_SSOILXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_QINSURXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_QSEVAXY
    REAL,    DIMENSION( ims:ime, 1:NSOIL, jms:jme ), INTENT(INOUT) ::  ACC_ETRANIXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  EFLXBXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SOILENERGY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SNOWENERGY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  CANHSXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_DWATERXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_PRCPXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_ECANXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_ETRANXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_EDIRXY

    INTEGER,  INTENT(IN   )   ::     ids,ide, jds,jde, kds,kde,  &  ! d -> domain
         &                           ims,ime, jms,jme, kms,kme,  &  ! m -> memory
         &                           its,ite, jts,jte, kts,kte      ! t -> tile

!2D inout irrigation variables
    REAL,    DIMENSION( ims:ime,          jms:jme ),  INTENT(IN)    :: IRFRACT    ! irrigation fraction
    REAL,    DIMENSION( ims:ime,          jms:jme ),  INTENT(IN)    :: SIFRACT    ! sprinkler irrigation fraction
    REAL,    DIMENSION( ims:ime,          jms:jme ),  INTENT(IN)    :: MIFRACT    ! micro irrigation fraction
    REAL,    DIMENSION( ims:ime,          jms:jme ),  INTENT(IN)    :: FIFRACT    ! flood irrigation fraction
    INTEGER, DIMENSION( ims:ime,          jms:jme ),  INTENT(INOUT) :: IRNUMSI    ! irrigation event number, Sprinkler
    INTEGER, DIMENSION( ims:ime,          jms:jme ),  INTENT(INOUT) :: IRNUMMI    ! irrigation event number, Micro
    INTEGER, DIMENSION( ims:ime,          jms:jme ),  INTENT(INOUT) :: IRNUMFI    ! irrigation event number, Flood
    REAL,    DIMENSION( ims:ime,          jms:jme ),  INTENT(INOUT) :: IRWATSI    ! irrigation water amount [m] to be applied, Sprinkler
    REAL,    DIMENSION( ims:ime,          jms:jme ),  INTENT(INOUT) :: IRWATMI    ! irrigation water amount [m] to be applied, Micro
    REAL,    DIMENSION( ims:ime,          jms:jme ),  INTENT(INOUT) :: IRWATFI    ! irrigation water amount [m] to be applied, Flood
    REAL,    DIMENSION( ims:ime,          jms:jme ),  INTENT(INOUT) :: IRELOSS    ! loss of irrigation water to evaporation,sprinkler [m/timestep]
    REAL,    DIMENSION( ims:ime,          jms:jme ),  INTENT(INOUT) :: IRSIVOL    ! amount of irrigation by sprinkler (mm)
    REAL,    DIMENSION( ims:ime,          jms:jme ),  INTENT(INOUT) :: IRMIVOL    ! amount of irrigation by micro (mm)
    REAL,    DIMENSION( ims:ime,          jms:jme ),  INTENT(INOUT) :: IRFIVOL    ! amount of irrigation by micro (mm)
    REAL,    DIMENSION( ims:ime,          jms:jme ),  INTENT(INOUT) :: IRRSPLH    ! latent heating from sprinkler evaporation (w/m2)
    CHARACTER(LEN=256),                               INTENT(IN)    :: LLANDUSE   ! landuse data name (USGS or MODIS_IGBP)

!ID local irrigation variables
    REAL                                                            :: IRRFRA     ! irrigation fraction
    REAL                                                            :: SIFAC      ! sprinkler irrigation fraction
    REAL                                                            :: MIFAC      ! micro irrigation fraction
    REAL                                                            :: FIFAC      ! flood irrigation fraction
    INTEGER                                                         :: IRCNTSI    ! irrigation event number, Sprinkler
    INTEGER                                                         :: IRCNTMI    ! irrigation event number, Micro
    INTEGER                                                         :: IRCNTFI    ! irrigation event number, Flood
    REAL                                                            :: IRAMTSI    ! irrigation water amount [m] to be applied, Sprinkler
    REAL                                                            :: IRAMTMI    ! irrigation water amount [m] to be applied, Micro
    REAL                                                            :: IRAMTFI    ! irrigation water amount [m] to be applied, Flood
    REAL                                                            :: IREVPLOS   ! loss of irrigation water to evaporation,sprinkler [m/timestep]
    REAL                                                            :: IRSIRATE   ! rate of irrigation by sprinkler [m/timestep]
    REAL                                                            :: IRMIRATE   ! rate of irrigation by micro [m/timestep]
    REAL                                                            :: IRFIRATE   ! rate of irrigation by micro [m/timestep]
    REAL                                                            :: FIRR       ! latent heating due to sprinkler evaporation (w m-2)
    REAL                                                            :: EIRR       ! evaporation due to sprinkler evaporation (mm/s)

! 1D equivalent of 2D/3D fields

! IN only

    REAL                                :: COSZ         ! cosine zenith angle
    REAL                                :: LAT          ! latitude [rad]
    REAL                                :: Z_ML         ! model height [m]
    INTEGER                             :: VEGTYP       ! vegetation type
    INTEGER,    DIMENSION(NSOIL)        :: SOILTYP      ! soil type
    INTEGER                             :: CROPTYPE     ! crop type
    REAL                                :: FVEG         ! vegetation fraction [-]
    REAL                                :: FVGMAX       ! annual max vegetation fraction []
    REAL                                :: TBOT         ! deep soil temperature [K]
    REAL                                :: T_ML         ! temperature valid at mid-levels [K]
    REAL                                :: Q_ML         ! water vapor mixing ratio [kg/kg_dry]
    REAL                                :: U_ML         ! U wind component [m/s]
    REAL                                :: V_ML         ! V wind component [m/s]
    REAL                                :: SWDN         ! solar down at surface [W m-2]
    REAL                                :: LWDN         ! longwave down at surface [W m-2]
    REAL                                :: P_ML         ! pressure, valid at interface [Pa]
    REAL                                :: PSFC         ! surface pressure [Pa]
    REAL                                :: PRCP         ! total precipitation entering  [mm/s]         ! MB/AN : v3.7
    REAL                                :: PRCPCONV     ! convective precipitation entering  [mm/s]    ! MB/AN : v3.7
    REAL                                :: PRCPNONC     ! non-convective precipitation entering [mm/s] ! MB/AN : v3.7
    REAL                                :: PRCPSHCV     ! shallow convective precip entering  [mm/s]   ! MB/AN : v3.7
    REAL                                :: PRCPSNOW     ! snow entering land model [mm/s]              ! MB/AN : v3.7
    REAL                                :: PRCPGRPL     ! graupel entering land model [mm/s]           ! MB/AN : v3.7
    REAL                                :: PRCPHAIL     ! hail entering land model [mm/s]              ! MB/AN : v3.7
    REAL                                :: PRCPOTHR     ! other precip, e.g. fog [mm/s]                ! MB/AN : v3.7

! INOUT (with generic LSM equivalent)

    REAL                                :: FSH          ! total sensible heat (w/m2) [+ to atm]
    REAL                                :: SSOIL        ! soil heat heat (w/m2)
    REAL                                :: SALB         ! surface albedo (-)
    REAL                                :: FSNO         ! snow cover fraction (-)
    REAL,   DIMENSION( 1:NSOIL)         :: SMCEQ        ! eq vol. soil moisture (m3/m3)
    REAL,   DIMENSION( 1:NSOIL)         :: SMC          ! vol. soil moisture (m3/m3)
    REAL,   DIMENSION( 1:NSOIL)         :: SMH2O        ! vol. soil liquid water (m3/m3)
    REAL,   DIMENSION(-2:NSOIL)         :: STC          ! snow/soil tmperatures
    REAL                                :: SWE          ! snow water equivalent (mm)
    REAL                                :: SNDPTH       ! snow depth (m)
    REAL                                :: EMISSI       ! net surface emissivity
    REAL                                :: QSFC1D       ! bulk surface specific humidity

! INOUT (with no Noah LSM equivalent)

    INTEGER                             :: ISNOW        ! actual no. of snow layers
    REAL                                :: TV           ! vegetation canopy temperature
    REAL                                :: TG           ! ground surface temperature
    REAL                                :: CANICE       ! canopy-intercepted ice (mm)
    REAL                                :: CANLIQ       ! canopy-intercepted liquid water (mm)
    REAL                                :: EAH          ! canopy air vapor pressure (pa)
    REAL                                :: TAH          ! canopy air temperature (k)
    REAL                                :: CM           ! momentum drag coefficient
    REAL                                :: CH           ! sensible heat exchange coefficient
    REAL                                :: FWET         ! wetted or snowed fraction of the canopy (-)
    REAL                                :: SNEQVO       ! snow mass at last time step(mm h2o)
    REAL                                :: ALBOLD       ! snow albedo at last time step (-)
    REAL                                :: QSNOW        ! snowfall on the ground [mm/s]
    REAL                                :: QRAIN        ! rainfall on the ground [mm/s]
    REAL                                :: WSLAKE       ! lake water storage [mm]
    REAL                                :: ZWT          ! water table depth [m]
    REAL                                :: WA           ! water in the "aquifer" [mm]
    REAL                                :: WT           ! groundwater storage [mm]
    REAL                                :: SMCWTD       ! soil moisture content in the layer to the water table when deep
    REAL                                :: DEEPRECH     ! recharge to the water table when deep
    REAL                                :: RECH         ! recharge to the water table (diagnostic)
    REAL, DIMENSION(-2:NSOIL)           :: ZSNSO        ! snow layer depth [m]
    REAL, DIMENSION(-2:              0) :: SNICE        ! snow layer ice [mm]
    REAL, DIMENSION(-2:              0) :: SNLIQ        ! snow layer liquid water [mm]
    REAL                                :: LFMASS       ! leaf mass [g/m2]
    REAL                                :: RTMASS       ! mass of fine roots [g/m2]
    REAL                                :: STMASS       ! stem mass [g/m2]
    REAL                                :: WOOD         ! mass of wood (incl. woody roots) [g/m2]
    REAL                                :: GRAIN        ! mass of grain XING [g/m2]
    REAL                                :: GDD          ! mass of grain XING[g/m2]
    INTEGER                             :: PGS          !stem respiration [g/m2/s]
    REAL                                :: STBLCP       ! stable carbon in deep soil [g/m2]
    REAL                                :: FASTCP       ! short-lived carbon, shallow soil [g/m2]
    REAL                                :: PLAI         ! leaf area index
    REAL                                :: PSAI         ! stem area index
    REAL                                :: TAUSS        ! non-dimensional snow age

! tile drainage
    REAL                                :: QTLDRN       ! tile drainage (mm)
    REAL                                :: TDFRACMP     ! tile drainage map
#ifdef WRF_HYDRO
    REAL                                :: WATBLED      ! water table depth for tile drainage
#endif

! OUT (with no Noah LSM equivalent)

    REAL                                :: Z0WRF        ! combined z0 sent to coupled model
    REAL                                :: T2MV         ! 2m temperature of vegetation part
    REAL                                :: T2MB         ! 2m temperature of bare ground part
    REAL                                :: Q2MV         ! 2m mixing ratio of vegetation part
    REAL                                :: Q2MB         ! 2m mixing ratio of bare ground part
    REAL                                :: TRAD         ! surface radiative temperature (k)
    REAL                                :: NEE          ! net ecosys exchange (g/m2/s CO2)
    REAL                                :: GPP          ! gross primary assimilation [g/m2/s C]
    REAL                                :: NPP          ! net primary productivity [g/m2/s C]
    REAL                                :: FVEGMP       ! greenness vegetation fraction [-]
    REAL                                :: RUNSF        ! surface runoff [mm] per soil timestep
    REAL                                :: RUNSB        ! subsurface runoff [mm] per soil timestep
    REAL                                :: ECAN         ! evaporation of intercepted water (mm/s)
    REAL                                :: ETRAN        ! transpiration rate (mm/s)
    REAL                                :: ESOIL        ! soil surface evaporation rate (mm/s]
    REAL                                :: FSA          ! total absorbed solar radiation (w/m2)
    REAL                                :: FIRA         ! total net longwave rad (w/m2) [+ to atm]
    REAL                                :: APAR         ! photosyn active energy by canopy (w/m2)
    REAL                                :: PSN          ! total photosynthesis (umol co2/m2/s) [+]
    REAL                                :: SAV          ! solar rad absorbed by veg. (w/m2)
    REAL                                :: SAG          ! solar rad absorbed by ground (w/m2)
    REAL                                :: RSSUN        ! sunlit leaf stomatal resistance (s/m)
    REAL                                :: RSSHA        ! shaded leaf stomatal resistance (s/m)
    REAL, DIMENSION(1:2)                :: ALBSND       ! snow albedo (direct)
    REAL, DIMENSION(1:2)                :: ALBSNI       ! snow albedo (diffuse)
    REAL                                :: RB           ! leaf boundary layer resistance (s/m)
    REAL                                :: LAISUN       ! sunlit leaf area index (m2/m2)
    REAL                                :: LAISHA       ! shaded leaf area index (m2/m2)
    REAL                                :: BGAP         ! between gap fraction
    REAL                                :: WGAP         ! within gap fraction
    REAL                                :: TGV          ! under canopy ground temperature[K]
    REAL                                :: TGB          ! bare ground temperature [K]
    REAL                                :: CHV          ! sensible heat exchange coefficient vegetated
    REAL                                :: CHB          ! sensible heat exchange coefficient bare-ground
    REAL                                :: IRC          ! canopy net LW rad. [w/m2] [+ to atm]
    REAL                                :: IRG          ! veg ground net LW rad. [w/m2] [+ to atm]
    REAL                                :: SHC          ! canopy sen. heat [w/m2]   [+ to atm]
    REAL                                :: SHG          ! veg ground sen. heat [w/m2]   [+ to atm]
    REAL                                :: EVG          ! veg ground evap. heat [w/m2]  [+ to atm]
    REAL                                :: GHV          ! veg ground heat flux [w/m2]  [+ to soil]
    REAL                                :: IRB          ! bare net longwave rad. [w/m2] [+ to atm]
    REAL                                :: SHB          ! bare sensible heat [w/m2]     [+ to atm]
    REAL                                :: EVB          ! bare evaporation heat [w/m2]  [+ to atm]
    REAL                                :: GHB          ! bare ground heat flux [w/m2] [+ to soil]
    REAL                                :: TR           ! transpiration [w/m2]  [+ to atm]
    REAL                                :: EVC          ! canopy evaporation heat [w/m2]  [+ to atm]
    REAL                                :: CHLEAF       ! leaf exchange coefficient
    REAL                                :: CHUC         ! under canopy exchange coefficient
    REAL                                :: CHV2         ! veg 2m exchange coefficient
    REAL                                :: CHB2         ! bare 2m exchange coefficient
    REAL                                :: QINTS
    REAL                                :: QINTR
    REAL                                :: QDRIPS
    REAL                                :: QDRIPR
    REAL                                :: QTHROS
    REAL                                :: QTHROR
    REAL                                :: QSNSUB
    REAL                                :: QSNFRO
    REAL                                :: QEVAC
    REAL                                :: QDEWC
    REAL                                :: QSUBC
    REAL                                :: QFROC
    REAL                                :: QFRZC
    REAL                                :: QMELTC
    REAL                                :: PAHV    !precipitation advected heat - vegetation net (W/m2)
    REAL                                :: PAHG    !precipitation advected heat - under canopy net (W/m2)
    REAL                                :: PAHB    !precipitation advected heat - bare ground net (W/m2)
    REAL                                :: PAH     !precipitation advected heat - total (W/m2)
    REAL                                :: RAININ  !rain rate                   (mm/s)
    REAL                                :: SNOWIN  !liquid equivalent snow rate (mm/s)
    REAL                                :: ACC_SSOIL
    REAL                                :: ACC_QINSUR
    REAL                                :: ACC_QSEVA
    REAL, DIMENSION( 1:NSOIL)           :: ACC_ETRANI       !transpiration rate (mm/s) [+]
    REAL                                :: EFLXB
    REAL                                :: XMF
    REAL, DIMENSION( -2:NSOIL )         :: HCPCT
    REAL                                :: DZSNSO
    REAL                                :: CANHS   ! canopy heat storage change (w/m2)
    REAL                                :: ACC_DWATER
    REAL                                :: ACC_PRCP
    REAL                                :: ACC_ECAN
    REAL                                :: ACC_ETRAN
    REAL                                :: ACC_EDIR

! Intermediate terms
    REAL                                :: FPICE        ! snow fraction of precip
    REAL                                :: FCEV         ! canopy evaporation heat (w/m2) [+ to atm]
    REAL                                :: FGEV         ! ground evaporation heat (w/m2) [+ to atm]
    REAL                                :: FCTR         ! transpiration heat flux (w/m2) [+ to atm]
    REAL                                :: QSNBOT       ! total liquid water (snowmelt + rain through pack)out of snowpack bottom [mm/s]
    REAL                                :: QMELT        ! snowmelt due to phase change (mm/s)
    REAL                                :: PONDING      ! snowmelt with no pack [mm]
    REAL                                :: PONDING1     ! snowmelt with no pack [mm]
    REAL                                :: PONDING2     ! snowmelt with no pack [mm]

! Local terms

    REAL, DIMENSION(1:60)               :: gecros1d     !  gecros crop
    REAL                                :: gecros_dd ,gecros_tbem,gecros_emb ,gecros_ema, &
                                           gecros_ds1,gecros_ds2 ,gecros_ds1x,gecros_ds2x

    REAL                                :: FSR          ! total reflected solar radiation (w/m2)
    REAL, DIMENSION(-2:0)               :: FICEOLD      ! snow layer ice fraction []
    REAL                                :: CO2PP        ! CO2 partial pressure [Pa]
    REAL                                :: O2PP         ! O2 partial pressure [Pa]
    REAL, DIMENSION(1:NSOIL)            :: ZSOIL        ! depth to soil interfaces [m]
    REAL                                :: FOLN         ! nitrogen saturation [%]

    REAL                                :: QC           ! cloud specific humidity for MYJ [not used]
    REAL                                :: PBLH         ! PBL height for MYJ [not used]
    REAL                                :: DZ8W1D       ! model level heights for MYJ [not used]

    INTEGER                             :: I
    INTEGER                             :: J
    INTEGER                             :: K
    INTEGER                             :: ICE
    INTEGER                             :: SLOPETYP
    LOGICAL                             :: IPRINT

    INTEGER                             :: SOILCOLOR          ! soil color index
    INTEGER                             :: IST          ! surface type 1-soil; 2-lake
    INTEGER                             :: YEARLEN
    REAL                                :: SOLAR_TIME
    INTEGER                             :: JMONTH, JDAY

    INTEGER, PARAMETER                  :: NSNOW = 3    ! number of snow layers fixed to 3
    REAL, PARAMETER                     :: undefined_value = -1.E36

    REAL, DIMENSION( 1:nsoil ) :: SAND
    REAL, DIMENSION( 1:nsoil ) :: CLAY
    REAL, DIMENSION( 1:nsoil ) :: ORGM
    ! added for complience with Noah-mp HUE scheme (doesnt change a thing, jsut stops an error) Aaron A.
    REAL,    DIMENSION( 1:nsoil)  ::  SMC_intermediate     ! volumetric soil moisture [m3/m3]
    REAL,    DIMENSION( 1:nsoil)  ::  SH2O_intermediate      ! volumetric liquid soil moisture [m3/m3]

    type(noahmp_parameters) :: parameters


! ----------------------------------------------------------------------
    CALL NOAHMP_OPTIONS(IDVEG  ,IOPT_CRS  ,IOPT_BTR  ,IOPT_RUN  ,IOPT_SFC  ,IOPT_FRZ , &
                IOPT_INF  ,IOPT_RAD  ,IOPT_ALB  ,IOPT_SNF  ,IOPT_TBOT, IOPT_STC  ,     &
		IOPT_RSF  ,IOPT_SOIL ,IOPT_PEDO ,IOPT_CROP ,IOPT_IRR , IOPT_IRRM ,     &
                IOPT_INFDV,IOPT_TDRN,IOPT_MOSAIC, IOPT_HUE )
    IPRINT    =  .false.                     ! debug printout

! for using soil update timestep difference from noahmp main timestep
    calculate_soil = .false.
    soil_update_steps = nint(soiltstep/DT)  ! 3600 = 1 hour
    soil_update_steps = max(soil_update_steps,1)
    if ( soil_update_steps == 1 ) then
      ACC_SSOILXY  = 0.0
      ACC_QINSURXY = 0.0
      ACC_QSEVAXY  = 0.0
      ACC_ETRANIXY = 0.0
      ACC_DWATERXY = 0.0
      ACC_PRCPXY   = 0.0
      ACC_ECANXY   = 0.0
      ACC_ETRANXY  = 0.0
      ACC_EDIRXY   = 0.0
    end if
    if ( soil_update_steps > 1 ) then
     if ( mod(itimestep,soil_update_steps) == 1 ) then
      ACC_SSOILXY  = 0.0
      ACC_QINSURXY = 0.0
      ACC_QSEVAXY  = 0.0
      ACC_ETRANIXY = 0.0
      ACC_DWATERXY = 0.0
      ACC_PRCPXY   = 0.0
      ACC_ECANXY   = 0.0
      ACC_ETRANXY  = 0.0
      ACC_EDIRXY   = 0.0
     end if
    end if

    if (mod(itimestep,soil_update_steps) == 0) calculate_soil = .true.
! end soil timestep

    YEARLEN = 365                            ! find length of year for phenology (also S Hemisphere)
    if (mod(YR,4) == 0) then
       YEARLEN = 366
       if (mod(YR,100) == 0) then
          YEARLEN = 365
          if (mod(YR,400) == 0) then
             YEARLEN = 366
          endif
       endif
    endif

    ZSOIL(1) = -DZS(1)                    ! depth to soil interfaces (<0) [m]
    DO K = 2, NSOIL
       ZSOIL(K) = -DZS(K) + ZSOIL(K-1)
    END DO

    JLOOP : DO J=jts,jte

       IF(ITIMESTEP == 1)THEN
          DO I=its,ite
             IF((XLAND(I,J)-1.5) >= 0.) THEN    ! Open water case
                IF(XICE(I,J) == 1. .AND. IPRINT) PRINT *,' sea-ice at water point, I=',I,'J=',J
                SMSTAV(I,J) = 1.0
                SMSTOT(I,J) = 1.0
                DO K = 1, NSOIL
                   SMOIS(I,K,J) = 1.0
                    TSLB(I,K,J) = 273.16
                ENDDO
             ELSE
                IF(XICE(I,J) == 1.) THEN        ! Sea-ice case
                   SMSTAV(I,J) = 1.0
                   SMSTOT(I,J) = 1.0
                   DO K = 1, NSOIL
                      SMOIS(I,K,J) = 1.0
                   ENDDO
                ENDIF
             ENDIF
          ENDDO
       ENDIF                                                               ! end of initialization over ocean


!-----------------------------------------------------------------------
   ILOOP : DO I = its, ite

    IF (XICE(I,J) >= XICE_THRES) THEN
       ICE = 1                            ! Sea-ice point

       SH2O  (i,1:NSOIL,j) = 1.0
       XLAIXY(i,j)         = 0.01

       CYCLE ILOOP ! Skip any processing at sea-ice points

    ELSE

       IF((XLAND(I,J)-1.5) >= 0.) CYCLE ILOOP   ! Open water case

!     2D to 1D

! IN only

       COSZ   = COSZIN  (I,J)                         ! cos zenith angle []
       LAT    = XLAT  (I,J)                           ! latitude [rad]
       Z_ML   = 0.5*DZ8W(I,1,J)                       ! DZ8W: thickness of full levels; ZLVL forcing height [m]
       VEGTYP = IVGTYP(I,J)                           ! vegetation type
       if(iopt_soil == 1) then
         SOILTYP= ISLTYP(I,J)                         ! soil type same in all layers
       elseif(iopt_soil == 2) then
         SOILTYP(1) = nint(SOILCL1(I,J))              ! soil type in layer1
         SOILTYP(2) = nint(SOILCL2(I,J))              ! soil type in layer2
         SOILTYP(3) = nint(SOILCL3(I,J))              ! soil type in layer3
         SOILTYP(4) = nint(SOILCL4(I,J))              ! soil type in layer4
       elseif(iopt_soil == 3) then
         SOILTYP= ISLTYP(I,J)                         ! to initialize with default
       end if
       FVEG   = VEGFRA(I,J)/100.                      ! vegetation fraction [0-1]
       FVGMAX = VEGMAX (I,J)/100.                     ! Vegetation fraction annual max [0-1]
       TBOT = TMN(I,J)                                ! Fixed deep soil temperature for land
       T_ML   = T3D(I,1,J)                            ! temperature defined at intermediate level [K]
       Q_ML   = QV3D(I,1,J)/(1.0+QV3D(I,1,J))         ! convert from mixing ratio to specific humidity [kg/kg]
       U_ML   = U_PHY(I,1,J)                          ! u-wind at interface [m/s]
       V_ML   = V_PHY(I,1,J)                          ! v-wind at interface [m/s]
       SWDN   = SWDOWN(I,J)                           ! shortwave down from SW scheme [W/m2]
       LWDN   = GLW(I,J)                              ! total longwave down from LW scheme [W/m2]
       P_ML   =(P8W3D(I,KTS+1,J)+P8W3D(I,KTS,J))*0.5  ! surface pressure defined at intermediate level [Pa]
	                                              !    consistent with temperature, mixing ratio
       PSFC   = P8W3D(I,1,J)                          ! surface pressure defined a full levels [Pa]
       PRCP   = PRECIP_IN (I,J) / DT                  ! timestep total precip rate (glacier) [mm/s]! MB: v3.7

       CROPTYPE = 0
       IF (IOPT_CROP > 0 .AND. VEGTYP == ISCROP_TABLE) CROPTYPE = DEFAULT_CROP_TABLE ! default croptype is generic dynamic vegetation crop
       IF (IOPT_CROP > 0 .AND. CROPCAT(I,J) > 0) THEN
         CROPTYPE = CROPCAT(I,J)                      ! crop type
	 VEGTYP = ISCROP_TABLE
         FVGMAX = 0.95
	 FVEG   = 0.95
       END IF

       IF (PRESENT(MP_RAINC) .AND. PRESENT(MP_RAINNC) .AND. PRESENT(MP_SHCV) .AND. &
           PRESENT(MP_SNOW)  .AND. PRESENT(MP_GRAUP)  .AND. PRESENT(MP_HAIL)   ) THEN

         PRCPCONV  = MP_RAINC (I,J)/DT                ! timestep convective precip rate [mm/s]     ! MB: v3.7
         PRCPNONC  = MP_RAINNC(I,J)/DT                ! timestep non-convective precip rate [mm/s] ! MB: v3.7
         PRCPSHCV  = MP_SHCV(I,J)  /DT                ! timestep shallow conv precip rate [mm/s]   ! MB: v3.7
         PRCPSNOW  = MP_SNOW(I,J)  /DT                ! timestep snow precip rate [mm/s]           ! MB: v3.7
         PRCPGRPL  = MP_GRAUP(I,J) /DT                ! timestep graupel precip rate [mm/s]        ! MB: v3.7
         PRCPHAIL  = MP_HAIL(I,J)  /DT                ! timestep hail precip rate [mm/s]           ! MB: v3.7

         PRCPOTHR  = PRCP - PRCPCONV - PRCPNONC - PRCPSHCV ! take care of other (fog) contained in rainbl
	 PRCPOTHR  = MAX(0.0,PRCPOTHR)
	 PRCPNONC  = PRCPNONC + PRCPOTHR
         PRCPSNOW  = PRCPSNOW + SR(I,J)  * PRCPOTHR
       ELSE
         PRCPCONV  = 0.
         PRCPNONC  = PRCP
         PRCPSHCV  = 0.
         PRCPSNOW  = SR(I,J) * PRCP
         PRCPGRPL  = 0.
         PRCPHAIL  = 0.
       ENDIF

! IN/OUT fields

       ISNOW                 = ISNOWXY (I,J)                ! snow layers []
       SMC  (      1:NSOIL)  = SMOIS   (I,      1:NSOIL,J)  ! soil total moisture [m3/m3]
       SMH2O(      1:NSOIL)  = SH2O    (I,      1:NSOIL,J)  ! soil liquid moisture [m3/m3]
       ! added by Aaron A.
       SMC_intermediate(1:NSOIL) = SMOIS(I, 1:NSOIL,J)     ! volumetric soil moisture [m3/m3] Aaron a. doesnt change a thing
       SH2O_intermediate(1:NSOIL) = SH2O(I, 1:NSOIL,J)      ! volumetric liquid soil moisture [m3/m3] Aaron a. Doesnt change a thing
       ! end added by Aaron A.
       STC  (-NSNOW+1:    0) = TSNOXY  (I,-NSNOW+1:    0,J) ! snow temperatures [K]
       STC  (      1:NSOIL)  = TSLB    (I,      1:NSOIL,J)  ! soil temperatures [K]
       SWE                   = SNOW    (I,J)                ! snow water equivalent [mm]
       SNDPTH                = SNOWH   (I,J)                ! snow depth [m]
       QSFC1D                = QSFC    (I,J)

! INOUT (with no Noah LSM equivalent)

       TV                    = TVXY    (I,J)                ! leaf temperature [K]
       TG                    = TGXY    (I,J)                ! ground temperature [K]
       CANLIQ                = CANLIQXY(I,J)                ! canopy liquid water [mm]
       CANICE                = CANICEXY(I,J)                ! canopy frozen water [mm]
       EAH                   = EAHXY   (I,J)                ! canopy vapor pressure [Pa]
       TAH                   = TAHXY   (I,J)                ! canopy temperature [K]
       CM                    = CMXY    (I,J)                ! avg. momentum exchange (MP only) [m/s]
       CH                    = CHXY    (I,J)                ! avg. heat exchange (MP only) [m/s]
       FWET                  = FWETXY  (I,J)                ! canopy fraction wet or snow
       SNEQVO                = SNEQVOXY(I,J)                ! SWE previous timestep
       ALBOLD                = ALBOLDXY(I,J)                ! albedo previous timestep, for snow aging
       QSNOW                 = QSNOWXY (I,J)                ! snow falling on ground
       QRAIN                 = QRAINXY (I,J)                ! rain falling on ground
       WSLAKE                = WSLAKEXY(I,J)                ! lake water storage (can be neg.) (mm)
       ZWT                   = ZWTXY   (I,J)                ! depth to water table [m]
       WA                    = WAXY    (I,J)                ! water storage in aquifer [mm]
       WT                    = WTXY    (I,J)                ! water in aquifer&saturated soil [mm]
       ZSNSO(-NSNOW+1:NSOIL) = ZSNSOXY (I,-NSNOW+1:NSOIL,J) ! depth to layer interface
       SNICE(-NSNOW+1:    0) = SNICEXY (I,-NSNOW+1:    0,J) ! snow layer ice content
       SNLIQ(-NSNOW+1:    0) = SNLIQXY (I,-NSNOW+1:    0,J) ! snow layer water content
       LFMASS                = LFMASSXY(I,J)                ! leaf mass
       RTMASS                = RTMASSXY(I,J)                ! root mass
       STMASS                = STMASSXY(I,J)                ! stem mass
       WOOD                  = WOODXY  (I,J)                ! mass of wood (incl. woody roots) [g/m2]
       STBLCP                = STBLCPXY(I,J)                ! stable carbon pool
       FASTCP                = FASTCPXY(I,J)                ! fast carbon pool
       PLAI                  = XLAIXY  (I,J)                ! leaf area index [-] (no snow effects)
       PSAI                  = XSAIXY  (I,J)                ! stem area index [-] (no snow effects)
       TAUSS                 = TAUSSXY (I,J)                ! non-dimensional snow age
       SMCEQ(       1:NSOIL) = SMOISEQ (I,       1:NSOIL,J)
       SMCWTD                = SMCWTDXY(I,J)
       RECH                  = 0.
       DEEPRECH              = 0.
       ACC_SSOIL             = ACC_SSOILXY (I,J)                 ! surface heat flux
       ACC_QSEVA             = ACC_QSEVAXY (I,J)
       ACC_QINSUR            = ACC_QINSURXY(I,J)
       ACC_ETRANI            = ACC_ETRANIXY(I,:,J)
       ACC_DWATER            = ACC_DWATERXY(I,J)
       ACC_PRCP              = ACC_PRCPXY  (I,J)
       ACC_ECAN              = ACC_ECANXY  (I,J)
       ACC_ETRAN             = ACC_ETRANXY (I,J)
       ACC_EDIR              = ACC_EDIRXY  (I,J)

! tile drainage
       QTLDRN                = 0.                           ! tile drainage (mm)
       TDFRACMP              = TD_FRACTION(I,J)             ! tile drainage map
#ifdef WRF_HYDRO
       WATBLED               = ZWATBLE2D (I,J)              ! water table depth for tile drainage
#endif

! irrigation vars
       IRRFRA                = IRFRACT(I,J)    ! irrigation fraction
       SIFAC                 = SIFRACT(I,J)    ! sprinkler irrigation fraction
       MIFAC                 = MIFRACT(I,J)    ! micro irrigation fraction
       FIFAC                 = FIFRACT(I,J)    ! flood irrigation fraction
       IRCNTSI               = IRNUMSI(I,J)    ! irrigation event number, Sprinkler
       IRCNTMI               = IRNUMMI(I,J)    ! irrigation event number, Micro
       IRCNTFI               = IRNUMFI(I,J)    ! irrigation event number, Flood
       IRAMTSI               = IRWATSI(I,J)    ! irrigation water amount [m] to be applied, Sprinkler
       IRAMTMI               = IRWATMI(I,J)    ! irrigation water amount [m] to be applied, Micro
       IRAMTFI               = IRWATFI(I,J)    ! irrigation water amount [m] to be applied, Flood
       IREVPLOS              = 0.0             ! loss of irrigation water to evaporation,sprinkler [m/timestep]
       IRSIRATE              = 0.0             ! rate of irrigation by sprinkler (mm)
       IRMIRATE              = 0.0             ! rate of irrigation by micro (mm)
       IRFIRATE              = 0.0             ! rate of irrigation by micro (mm)
       FIRR                  = 0.0             ! latent heating due to sprinkler evaporation (W m-2)
       EIRR                  = 0.0             ! evaporation from sprinkler (mm/s)

       if(iopt_crop == 2) then   ! gecros crop model

         gecros1d(1:60)      = gecros_state(I,1:60,J)       ! Gecros variables 2D -> local

         if(croptype == 1) then
           gecros_dd   =  2.5
           gecros_tbem =  2.0
           gecros_emb  = 10.2
           gecros_ema  = 40.0
           gecros_ds1  =  2.1 !BBCH 92
           gecros_ds2  =  2.0 !BBCH 90
           gecros_ds1x =  0.0
           gecros_ds2x = 10.0
         end if

         if(croptype == 2) then
           gecros_dd   =  5.0
           gecros_tbem =  8.0
           gecros_emb  = 15.0
           gecros_ema  =  6.0
           gecros_ds1  =  1.78  !BBCH 85
           gecros_ds2  =  1.63  !BBCH 80
           gecros_ds1x =  0.0
           gecros_ds2x = 14.0
         end if

       end if

       SLOPETYP     = 1                               ! set underground runoff slope term
       IST          = 1                               ! MP surface type: 1 = land; 2 = lake
       SOILCOLOR    = 4                               ! soil color: assuming a middle color category ?????????

       IF(any(SOILTYP == 14) .AND. XICE(I,J) == 0.) THEN
          IF(IPRINT) PRINT *, ' SOIL TYPE FOUND TO BE WATER AT A LAND-POINT'
          IF(IPRINT) PRINT *, i,j,'RESET SOIL in surfce.F'
          SOILTYP = 7
       ENDIF
         IF( IVGTYP(I,J) == ISURBAN_TABLE    .or. IVGTYP(I,J) == LCZ_1_TABLE .or. IVGTYP(I,J) == LCZ_2_TABLE .or. &
             IVGTYP(I,J) == LCZ_3_TABLE      .or. IVGTYP(I,J) == LCZ_4_TABLE .or. IVGTYP(I,J) == LCZ_5_TABLE .or. &
             IVGTYP(I,J) == LCZ_6_TABLE      .or. IVGTYP(I,J) == LCZ_7_TABLE .or. IVGTYP(I,J) == LCZ_8_TABLE .or. &
             IVGTYP(I,J) == LCZ_9_TABLE      .or. IVGTYP(I,J) == LCZ_10_TABLE .or. IVGTYP(I,J) == LCZ_11_TABLE ) THEN


         IF(SF_URBAN_PHYSICS == 0 ) THEN
           VEGTYP = ISURBAN_TABLE
         ELSE
           VEGTYP = NATURAL_TABLE  ! set urban vegetation type based on table natural
           FVGMAX = 0.96
         ENDIF

       ENDIF

! placeholders for 3D soil
!       parameters%bexp   = BEXP_3D  (I,1:NSOIL,J) ! C-H B exponent
!       parameters%smcdry = SMCDRY_3D(I,1:NSOIL,J) ! Soil Moisture Limit: Dry
!       parameters%smcwlt = SMCWLT_3D(I,1:NSOIL,J) ! Soil Moisture Limit: Wilt
!       parameters%smcref = SMCREF_3D(I,1:NSOIL,J) ! Soil Moisture Limit: Reference
!       parameters%smcmax = SMCMAX_3D(I,1:NSOIL,J) ! Soil Moisture Limit: Max
!       parameters%dksat  = DKSAT_3D (I,1:NSOIL,J) ! Saturated Soil Conductivity
!       parameters%dwsat  = DWSAT_3D (I,1:NSOIL,J) ! Saturated Soil Diffusivity
!       parameters%psisat = PSISAT_3D(I,1:NSOIL,J) ! Saturated Matric Potential
!       parameters%quartz = QUARTZ_3D(I,1:NSOIL,J) ! Soil quartz content
!       parameters%refdk  = REFDK_2D (I,J)         ! Reference Soil Conductivity
!       parameters%refkdt = REFKDT_2D(I,J)         ! Soil Infiltration Parameter

! placeholders for 2D additional runoff parameters
!       parameters%BVIC   = BVIC_2D(I,J)           ! VIC model infiltration parameter [-]
!       parameters%axaj   = AXAJ_2D(I,J)           ! Xinanjiang: Tension water distribution inflection parameter [-]
!       parameters%bxaj   = BXAJ_2D(I,J)           ! Xinanjiang: Tension water distribution shape parameter [-]
!       parameters%xxaj   = XXAJ_2D(I,J)           ! Xinanjiang: Free water distribution shape parameter [-]
!       parameters%BDVIC  = BDVIC_2D(I,J)          ! VIC model infiltration parameter [-]
!       parameters%GDVIC  = GDVIC_2D(I,J)          ! Mean Capillary Drive for infiltration models [m]
!       parameters%BBVIC  = BBVIC_2D(I,J)          ! DVIC heterogeniety parameter for infiltraton [-]

! placeholders for 2D irrigation params
!       parameters%IRR_FRAC   = IRR_FRAC_2D(I,J)   ! irrigation Fraction
!       parameters%IRR_HAR    = IRR_HAR_2D(I,J)    ! number of days before harvest date to stop irrigation
!       parameters%IRR_LAI    = IRR_LAI_2D(I,J)    ! Minimum lai to trigger irrigation
!       parameters%IRR_MAD    = IRR_MAD_2D(I,J)    ! management allowable deficit (0-1)
!       parameters%FILOSS     = FILOSS_2D(I,J)     ! fraction of flood irrigation loss (0-1)
!       parameters%SPRIR_RATE = SPRIR_RATE_2D(I,J) ! mm/h, sprinkler irrigation rate
!       parameters%MICIR_RATE = MICIR_RATE_2D(I,J) ! mm/h, micro irrigation rate
!       parameters%FIRTFAC    = FIRTFAC_2D(I,J)    ! flood application rate factor
!       parameters%IR_RAIN    = IR_RAIN_2D(I,J)    ! maximum precipitation to stop irrigation trigger

! placeholders for 2D tile drainage parameters
!       parameters%klat_fac  = KLAT_FAC (I,J)      ! factor multiplier to hydraulic conductivity
!       parameters%tdsmc_fac = TDSMC_FAC(I,J)      ! factor multiplier to field capacity
!       parameters%td_dc     = TD_DC    (I,J)      ! drainage coefficient for simple
!       parameters%td_dcoef  = TD_DCOEF (I,J)      ! drainge coefficient for Hooghoudt
!       parameters%td_ddrain = TD_DDRAIN(I,J)      ! depth of drain
!       parameters%td_radi   = TD_RADI  (I,J)      ! tile radius
!       parameters%td_spac   = TD_SPAC  (I,J)      ! tile spacing

       CALL TRANSFER_MP_PARAMETERS(VEGTYP,SOILTYP,SLOPETYP,SOILCOLOR,CROPTYPE,parameters)

       if(iopt_soil == 3 .and. .not. parameters%urban_flag) then

	sand = 0.01 * soilcomp(i,1:4,j)
	clay = 0.01 * soilcomp(i,5:8,j)
        orgm = 0.0

        if(opt_pedo == 1) call pedotransfer_sr2006(nsoil,sand,clay,orgm,parameters)

       end if

       GRAIN = GRAINXY (I,J)                ! mass of grain XING [g/m2]
       GDD   = GDDXY (I,J)                  ! growing degree days XING
       PGS   = PGSXY (I,J)                  ! growing degree days XING

       if(iopt_crop == 1 .and. croptype > 0) then
         parameters%PLTDAY = PLANTING(I,J)
	 parameters%HSDAY  = HARVEST (I,J)
	 parameters%GDDS1  = SEASON_GDD(I,J) / 1770.0 * parameters%GDDS1
	 parameters%GDDS2  = SEASON_GDD(I,J) / 1770.0 * parameters%GDDS2
	 parameters%GDDS3  = SEASON_GDD(I,J) / 1770.0 * parameters%GDDS3
	 parameters%GDDS4  = SEASON_GDD(I,J) / 1770.0 * parameters%GDDS4
	 parameters%GDDS5  = SEASON_GDD(I,J) / 1770.0 * parameters%GDDS5
       end if

       if(iopt_irr == 2) then
         parameters%PLTDAY = PLANTING(I,J)
         parameters%HSDAY  = HARVEST (I,J)
       end if

!=== hydrological processes for vegetation in urban model ===
!=== irrigate vegetaion only in urban area, MAY-SEP, 9-11pm

         IF( IVGTYP(I,J) == ISURBAN_TABLE    .or. IVGTYP(I,J) == LCZ_1_TABLE .or. IVGTYP(I,J) == LCZ_2_TABLE .or. &
               IVGTYP(I,J) == LCZ_3_TABLE      .or. IVGTYP(I,J) == LCZ_4_TABLE .or. IVGTYP(I,J) == LCZ_5_TABLE .or. &
               IVGTYP(I,J) == LCZ_6_TABLE      .or. IVGTYP(I,J) == LCZ_7_TABLE .or. IVGTYP(I,J) == LCZ_8_TABLE .or. &
               IVGTYP(I,J) == LCZ_9_TABLE      .or. IVGTYP(I,J) == LCZ_10_TABLE .or. IVGTYP(I,J) == LCZ_11_TABLE ) THEN

         IF(SF_URBAN_PHYSICS > 0 .AND. IRI_SCHEME == 1 ) THEN
	     SOLAR_TIME = (JULIAN - INT(JULIAN))*24 + XLONG(I,J)/15.0
	     IF(SOLAR_TIME < 0.) SOLAR_TIME = SOLAR_TIME + 24.
             CALL CAL_MON_DAY(INT(JULIAN),YR,JMONTH,JDAY)
             IF (SOLAR_TIME >= 21. .AND. SOLAR_TIME <= 23. .AND. JMONTH >= 5 .AND. JMONTH <= 9) THEN
                SMC(1) = max(SMC(1),parameters%SMCREF(1))
                SMC(2) = max(SMC(2),parameters%SMCREF(2))
             ENDIF
         ENDIF

       ENDIF

! Initialized local

       FICEOLD = 0.0
       FICEOLD(ISNOW+1:0) = SNICEXY(I,ISNOW+1:0,J) &  ! snow ice fraction
           /(SNICEXY(I,ISNOW+1:0,J)+SNLIQXY(I,ISNOW+1:0,J))
       CO2PP  = CO2_TABLE * P_ML                      ! partial pressure co2 [Pa]
       O2PP   = O2_TABLE  * P_ML                      ! partial pressure  o2 [Pa]
       FOLN   = 1.0                                   ! for now, set to nitrogen saturation
       QC     = undefined_value                       ! test dummy value
       PBLH   = undefined_value                       ! test dummy value ! PBL height
       DZ8W1D = DZ8W (I,1,J)                          ! thickness of atmospheric layers

       IF(VEGTYP == 25) FVEG = 0.0                  ! Set playa, lava, sand to bare
       IF(VEGTYP == 25) PLAI = 0.0
       IF(VEGTYP == 26) FVEG = 0.0                  ! hard coded for USGS
       IF(VEGTYP == 26) PLAI = 0.0
       IF(VEGTYP == 27) FVEG = 0.0
       IF(VEGTYP == 27) PLAI = 0.0

       IF ( VEGTYP == ISICE_TABLE ) THEN
         ICE = -1                           ! Land-ice point
         CALL NOAHMP_OPTIONS_GLACIER(IOPT_ALB  ,IOPT_SNF  ,IOPT_TBOT, IOPT_STC, IOPT_GLA )

         TBOT = MIN(TBOT,263.15)                      ! set deep temp to at most -10C
         CALL NOAHMP_GLACIER(     I,       J,    COSZ,   NSNOW,   NSOIL,      DT, & ! IN : Time/Space/Model-related
                               T_ML,    P_ML,    U_ML,    V_ML,    Q_ML,    SWDN, & ! IN : Forcing
                               PRCP,    LWDN,    TBOT,    Z_ML, FICEOLD,   ZSOIL, & ! IN : Forcing
                              QSNOW,  SNEQVO,  ALBOLD,      CM,      CH,   ISNOW, & ! IN/OUT :
                                SWE,     SMC,   ZSNSO,  SNDPTH,   SNICE,   SNLIQ, & ! IN/OUT :
                                 TG,     STC,   SMH2O,   TAUSS,  QSFC1D,          & ! IN/OUT :
                                FSA,     FSR,    FIRA,     FSH,    FGEV,   SSOIL, & ! OUT :
                               TRAD,   ESOIL,   RUNSF,   RUNSB,     SAG,    SALB, & ! OUT :
                              QSNBOT,PONDING,PONDING1,PONDING2,    T2MB,    Q2MB, & ! OUT :
			      EMISSI,  FPICE,    CHB2,   QMELT                    & ! OUT :
#ifdef WRF_HYDRO
                              , sfcheadrt(i,j)                                      &
#endif
                              )

         FSNO   = 1.0
         TV     = undefined_value     ! Output from standard Noah-MP undefined for glacier points
         TGB    = TG
         CANICE = undefined_value
         CANLIQ = undefined_value
         EAH    = undefined_value
         TAH    = undefined_value
         FWET   = undefined_value
         WSLAKE = undefined_value
!         ZWT    = undefined_value
         WA     = undefined_value
         WT     = undefined_value
         LFMASS = undefined_value
         RTMASS = undefined_value
         STMASS = undefined_value
         WOOD   = undefined_value
         QTLDRN = undefined_value
         GRAIN  = undefined_value
         GDD    = undefined_value
         STBLCP = undefined_value
         FASTCP = undefined_value
         PLAI   = undefined_value
         PSAI   = undefined_value
         T2MV   = undefined_value
         Q2MV   = undefined_value
         NEE    = undefined_value
         GPP    = undefined_value
         NPP    = undefined_value
         FVEGMP = 0.0
         ECAN   = 0.0
         ETRAN  = 0.0
         APAR   = undefined_value
         PSN    = undefined_value
         SAV    = 0.0
         RSSUN  = undefined_value
         RSSHA  = undefined_value
         RB     = undefined_value
         LAISUN = undefined_value
         LAISHA = undefined_value
         RS(I,J)= undefined_value
         BGAP   = undefined_value
         WGAP   = undefined_value
         TGV    = undefined_value
         CHV    = undefined_value
         CHB    = CH
         IRC    = 0.0
         IRG    = 0.0
         SHC    = 0.0
         SHG    = 0.0
         EVG    = 0.0
         GHV    = 0.0
         CANHS  = 0.0
         IRB    = FIRA
         SHB    = FSH
         EVB    = FGEV
         GHB    = SSOIL
         TR     = 0.0
         EVC    = 0.0
         PAH    = 0.0
         PAHG   = 0.0
         PAHB   = 0.0
         PAHV   = 0.0
         CHLEAF = undefined_value
         CHUC   = undefined_value
         CHV2   = undefined_value
         FCEV   = 0.0
         FCTR   = 0.0
         Z0WRF  = 0.002
         QFX(I,J) = ESOIL
         LH (I,J) = FGEV
         QINTS  = 0.0
         QINTR  = 0.0
         QDRIPS = 0.0
         QDRIPR = 0.0
         QTHROS = PRCP * FPICE
         QTHROR = PRCP * (1.0 - FPICE)
         QSNSUB = MAX( ESOIL, 0.)
         QSNFRO = ABS( MIN(ESOIL, 0.))
         QSUBC = 0.0
         QFROC = 0.0
         QFRZC = 0.0
         QMELTC = 0.0
         QEVAC = 0.0
         QDEWC = 0.0
         RAININ = PRCP * (1.0 - FPICE)
         SNOWIN = PRCP * FPICE
         CANICE = 0.0
         CANLIQ = 0.0
         QTLDRN = 0.0
         RUNSF  = RUNSF * dt
         RUNSB  = RUNSB * dt
    ELSE
         ICE=0                              ! Neither sea ice or land ice.
         CALL NOAHMP_SFLX (parameters, &
            I       , J       , LAT     , YEARLEN , JULIAN  , COSZ    , & ! IN : Time/Space-related
            DT      , DX      , DZ8W1D  , NSOIL   , ZSOIL   , NSNOW   , & ! IN : Model configuration
            FVEG    , FVGMAX  , VEGTYP  , ICE     , IST     , CROPTYPE, & ! IN : Vegetation/Soil characteristics
            SMCEQ   ,                                                   & ! IN : Vegetation/Soil characteristics
            T_ML    , P_ML    , PSFC    , U_ML    , V_ML    , Q_ML    , & ! IN : Forcing
            QC      , SWDN    , LWDN    ,                               & ! IN : Forcing
	    PRCPCONV, PRCPNONC, PRCPSHCV, PRCPSNOW, PRCPGRPL, PRCPHAIL, & ! IN : Forcing
            TBOT    , CO2PP   , O2PP    , FOLN    , FICEOLD , Z_ML    , & ! IN : Forcing
            IRRFRA  , SIFAC   , MIFAC   , FIFAC   , LLANDUSE,           & ! IN : Irrigation: fractions
            ALBOLD  , SNEQVO  ,                                         & ! IN/OUT :
            STC     , SMH2O   , SMC     , TAH     , EAH     , FWET    , & ! IN/OUT :
            CANLIQ  , CANICE  , TV      , TG      , QSFC1D  , QSNOW   , & ! IN/OUT :
            QRAIN   ,                                                   & ! IN/OUT :
            ISNOW   , ZSNSO   , SNDPTH  , SWE     , SNICE   , SNLIQ   , & ! IN/OUT :
            ZWT     , WA      , WT      , WSLAKE  , LFMASS  , RTMASS  , & ! IN/OUT :
            STMASS  , WOOD    , STBLCP  , FASTCP  , PLAI    , PSAI    , & ! IN/OUT :
            CM      , CH      , TAUSS   ,                               & ! IN/OUT :
            GRAIN   , GDD     , PGS     ,                               & ! IN/OUT
            SMCWTD  ,DEEPRECH , RECH    ,                               & ! IN/OUT :
            GECROS1D,                                                   & ! IN/OUT :
            QTLDRN  , TDFRACMP,                                         & ! IN/OUT : tile drainage
            Z0WRF   ,                                                   & ! OUT :
            IRCNTSI , IRCNTMI , IRCNTFI , IRAMTSI , IRAMTMI , IRAMTFI , & ! IN/OUT : Irrigation: vars
            IRSIRATE, IRMIRATE, IRFIRATE, FIRR    , EIRR    ,           & ! IN/OUT : Irrigation: vars
            FSA     , FSR     , FIRA    , FSH     , SSOIL   , FCEV    , & ! OUT :
            FGEV    , FCTR    , ECAN    , ETRAN   , ESOIL   , TRAD    , & ! OUT :
            TGB     , TGV     , T2MV    , T2MB    , Q2MV    , Q2MB    , & ! OUT :
            RUNSF   , RUNSB   , APAR    , PSN     , SAV     , SAG     , & ! OUT :
            FSNO    , NEE     , GPP     , NPP     , FVEGMP  , SALB    , & ! OUT :
            QSNBOT  , PONDING , PONDING1, PONDING2, RSSUN   , RSSHA   , & ! OUT :
            ALBSND  , ALBSNI  ,                                         & ! OUT :
            BGAP    , WGAP    , CHV     , CHB     , EMISSI  ,           & ! OUT :
            SHG     , SHC     , SHB     , EVG     , EVB     , GHV     , & ! OUT :
	    GHB     , IRG     , IRC     , IRB     , TR      , EVC     , & ! OUT :
	    CHLEAF  , CHUC    , CHV2    , CHB2    , FPICE   , PAHV    , & ! OUT :
            PAHG    , PAHB    , PAH     , LAISUN  , LAISHA  , RB      , & ! OUT :
            QINTS   , QINTR   , QDRIPS  , QDRIPR  , QTHROS  , QTHROR  , & ! OUT :
            QSNSUB  , QSNFRO  , QSUBC   , QFROC   , QFRZC   , QMELTC  , & ! OUT :
            QEVAC   , QDEWC   , QMELT   ,                               & ! OUT :
            RAININ  , SNOWIN  , ACC_SSOIL, ACC_QINSUR, ACC_QSEVA      , & ! OUT :
            ACC_ETRANI, HCPCT , EFLXB   , CANHS   ,                     & ! OUT :
            ACC_DWATER, ACC_PRCP, ACC_ECAN, ACC_ETRAN, ACC_EDIR         & ! INOUT
            ,SMC_INTERMEDIATE, SH2O_INTERMEDIATE                        & ! added by AARON A. changes nothing without mosaic
#ifdef WRF_HYDRO
            , sfcheadrt(i,j), WATBLED                                   &
#endif
            )            ! OUT :
            QFX(I,J) = ECAN + ESOIL + ETRAN + EIRR
            LH(I,J)  = FCEV + FGEV  + FCTR  + FIRR

   ENDIF ! glacial split ends

#ifdef WRF_HYDRO
!AD_CHANGE: Glacier cells can produce small negative subsurface runoff for mass balance.
!           This will crash channel routing, so only pass along positive runoff.
            !soldrain(i,j) = max(RUNSB*dt, 0.0)        !mm , underground runoff
            soldrain(i,j) = max(RUNSB, 0.0)        !mm , underground runoff
            INFXSRT(i,j) = RUNSF   !*dt        !mm , surface runoff
            qtiledrain (i,j) = QTLDRN   !*dt      !mm, tile drainage
#endif


! INPUT/OUTPUT

             TSK      (I,J)                = TRAD
             HFX      (I,J)                = FSH
             GRDFLX   (I,J)                = SSOIL
	     SMSTAV   (I,J)                = 0.0  ! [maintained as Noah consistency]
             SMSTOT   (I,J)                = 0.0  ! [maintained as Noah consistency]
             SFCRUNOFF(I,J)                = SFCRUNOFF(I,J) + RUNSF  !* DT
             UDRUNOFF (I,J)                = UDRUNOFF(I,J)  + RUNSB  !* DT
             QTDRAIN  (I,J)                = QTDRAIN (I,J)  + QTLDRN !* DT
             IF ( SALB > -999 ) THEN
                ALBEDO(I,J)                = SALB
             ENDIF
             SNOWC    (I,J)                = FSNO
             SMOIS    (I,      1:NSOIL,J)  = SMC   (      1:NSOIL)
             SH2O     (I,      1:NSOIL,J)  = SMH2O (      1:NSOIL)
             TSLB     (I,      1:NSOIL,J)  = STC   (      1:NSOIL)
             SNOW     (I,J)                = SWE
             SNOWH    (I,J)                = SNDPTH
             CANWAT   (I,J)                = CANLIQ + CANICE
             ACSNOW   (I,J)                = ACSNOW(I,J) +  PRECIP_IN(I,J) * FPICE
!             ACSNOM   (I,J)                = ACSNOM(I,J) + QSNBOT*DT + PONDING + PONDING1 + PONDING2
             ACSNOM   (I,J)                = ACSNOM(I,J) + QMELT*DT + PONDING + PONDING1 + PONDING2
             EMISS    (I,J)                = EMISSI
             QSFC     (I,J)                = QSFC1D

             ISNOWXY  (I,J)                = ISNOW
             TVXY     (I,J)                = TV
             TGXY     (I,J)                = TG
             CANLIQXY (I,J)                = CANLIQ
             CANICEXY (I,J)                = CANICE
             EAHXY    (I,J)                = EAH
             TAHXY    (I,J)                = TAH
             CMXY     (I,J)                = CM
             CHXY     (I,J)                = CH
             FWETXY   (I,J)                = FWET
             SNEQVOXY (I,J)                = SNEQVO
             ALBOLDXY (I,J)                = ALBOLD
             QSNOWXY  (I,J)                = QSNOW
             QRAINXY  (I,J)                = QRAIN
             WSLAKEXY (I,J)                = WSLAKE
             ZWTXY    (I,J)                = ZWT
             WAXY     (I,J)                = WA
             WTXY     (I,J)                = WT
             TSNOXY   (I,-NSNOW+1:    0,J) = STC   (-NSNOW+1:    0)
             ZSNSOXY  (I,-NSNOW+1:NSOIL,J) = ZSNSO (-NSNOW+1:NSOIL)
             SNICEXY  (I,-NSNOW+1:    0,J) = SNICE (-NSNOW+1:    0)
             SNLIQXY  (I,-NSNOW+1:    0,J) = SNLIQ (-NSNOW+1:    0)
             LFMASSXY (I,J)                = LFMASS
             RTMASSXY (I,J)                = RTMASS
             STMASSXY (I,J)                = STMASS
             WOODXY   (I,J)                = WOOD
             STBLCPXY (I,J)                = STBLCP
             FASTCPXY (I,J)                = FASTCP
             XLAIXY   (I,J)                = PLAI
             XSAIXY   (I,J)                = PSAI
             TAUSSXY  (I,J)                = TAUSS

! OUTPUT
             Z0       (I,J)                = Z0WRF
             ZNT      (I,J)                = Z0WRF
             T2MVXY   (I,J)                = T2MV
             T2MBXY   (I,J)                = T2MB
             Q2MVXY   (I,J)                = Q2MV/(1.0 - Q2MV)  ! specific humidity to mixing ratio
             Q2MBXY   (I,J)                = Q2MB/(1.0 - Q2MB)  ! consistent with registry def of Q2
             TRADXY   (I,J)                = TRAD
             NEEXY    (I,J)                = NEE
             GPPXY    (I,J)                = GPP
             NPPXY    (I,J)                = NPP
             FVEGXY   (I,J)                = FVEGMP
             RUNSFXY  (I,J)                = RUNSF
             RUNSBXY  (I,J)                = RUNSB
             ECANXY   (I,J)                = ECAN
             EDIRXY   (I,J)                = ESOIL
             ETRANXY  (I,J)                = ETRAN
             FSAXY    (I,J)                = FSA
             FIRAXY   (I,J)                = FIRA
             APARXY   (I,J)                = APAR
             PSNXY    (I,J)                = PSN
             SAVXY    (I,J)                = SAV
             SAGXY    (I,J)                = SAG
             RSSUNXY  (I,J)                = RSSUN
             RSSHAXY  (I,J)                = RSSHA
             LAISUN                        = MAX(LAISUN, 0.0)
             LAISHA                        = MAX(LAISHA, 0.0)
             RB                            = MAX(RB, 0.0)
! New Calculation of total Canopy/Stomatal Conductance Based on Bonan et al. (2011)
! -- Inverse of Canopy Resistance (below)
             IF(RSSUN .le. 0.0 .or. RSSHA .le. 0.0 .or. LAISUN .eq. 0.0 .or. LAISHA .eq. 0.0) THEN
                RS    (I,J)                = 0.0
             ELSE
                RS    (I,J)                = ((1.0/(RSSUN+RB)*LAISUN) + ((1.0/(RSSHA+RB))*LAISHA))
                RS    (I,J)                = 1.0/RS(I,J) !Resistance
             ENDIF
             BGAPXY   (I,J)                = BGAP
             WGAPXY   (I,J)                = WGAP
             TGVXY    (I,J)                = TGV
             TGBXY    (I,J)                = TGB
             CHVXY    (I,J)                = CHV
             CHBXY    (I,J)                = CHB
             IRCXY    (I,J)                = IRC
             IRGXY    (I,J)                = IRG
             SHCXY    (I,J)                = SHC
             SHGXY    (I,J)                = SHG
             EVGXY    (I,J)                = EVG
             GHVXY    (I,J)                = GHV
             IRBXY    (I,J)                = IRB
             SHBXY    (I,J)                = SHB
             EVBXY    (I,J)                = EVB
             GHBXY    (I,J)                = GHB
             canhsxy  (I,J)                = CANHS
             TRXY     (I,J)                = TR
             EVCXY    (I,J)                = EVC
             CHLEAFXY (I,J)                = CHLEAF
             CHUCXY   (I,J)                = CHUC
             CHV2XY   (I,J)                = CHV2
             CHB2XY   (I,J)                = CHB2
             PAHXY    (I,J)                = PAH
             PAHGXY   (I,J)                = PAHG
             PAHBXY   (I,J)                = PAHB
             PAHVXY   (I,J)                = PAHV
             QINTSXY  (I,J)                = QINTS
             QINTRXY  (I,J)                = QINTR
             QDRIPSXY (I,J)                = QDRIPS
             QDRIPRXY (I,J)                = QDRIPR
             QTHROSXY (I,J)                = QTHROS
             QTHRORXY (I,J)                = QTHROR
             QSNSUBXY (I,J)                = QSNSUB
             QSNFROXY (I,J)                = QSNFRO
             QSUBCXY  (I,J)                = QSUBC
             QFROCXY  (I,J)                = QFROC
             QEVACXY  (I,J)                = QEVAC
             QDEWCXY  (I,J)                = QDEWC
             QFRZCXY  (I,J)                = QFRZC
             QMELTCXY (I,J)                = QMELTC
             QSNBOTXY (I,J)                = QSNBOT
             QMELTXY  (I,J)                = QMELT
             PONDINGXY(I,J)                = PONDING + PONDING1 + PONDING2
             FPICEXY  (I,J)                = FPICE
             RAINLSM  (I,J)                = RAININ
             SNOWLSM  (I,J)                = SNOWIN
             FORCTLSM (I,J)                = T_ML
             FORCQLSM (I,J)                = Q_ML
             FORCPLSM (I,J)                = P_ML
             FORCZLSM (I,J)                = Z_ML
             FORCWLSM (I,J)                = SQRT(U_ML*U_ML + V_ML*V_ML)
             RECHXY   (I,J)                = RECHXY(I,J) + RECH*1.E3 !RECHARGE TO THE WATER TABLE
             DEEPRECHXY(I,J)               = DEEPRECHXY(I,J) + DEEPRECH
             SMCWTDXY(I,J)                 = SMCWTD
             ACC_SSOILXY (I,J)             = ACC_SSOIL
             ACC_QINSURXY(I,J)             = ACC_QINSUR
             ACC_QSEVAXY (I,J)             = ACC_QSEVA
             ACC_ETRANIXY(I,:,J)           = ACC_ETRANI
             ACC_DWATERXY(I,J)             = ACC_DWATER
             ACC_PRCPXY  (I,J)             = ACC_PRCP
             ACC_ECANXY  (I,J)             = ACC_ECAN
             ACC_ETRANXY (I,J)             = ACC_ETRAN
             ACC_EDIRXY  (I,J)             = ACC_EDIR
             EFLXBXY (I,J)                 = EFLXB
             SNOWENERGY(I,J)               = 0.0
             SOILENERGY(I,J)               = 0.0
             DO K = ISNOW+1, NSOIL
               IF(K == ISNOW+1) THEN
                 DZSNSO = - ZSNSO(K)
               ELSE
                 DZSNSO = ZSNSO(K-1) - ZSNSO(K)
               END IF
               IF(K >= 1) THEN
                 SOILENERGY(I,J) = SOILENERGY(I,J) + DZSNSO * HCPCT(K) * (STC(K)-273.16) * 0.001
               ELSE
                 SNOWENERGY(I,J) = SNOWENERGY(I,J) + DZSNSO * HCPCT(K) * (STC(K)-273.16) * 0.001
               END IF
             ENDDO

             GRAINXY  (I,J) = GRAIN !GRAIN XING
             GDDXY    (I,J) = GDD   !XING
	     PGSXY    (I,J) = PGS

             ! irrigation
             IRNUMSI(I,J)                  = IRCNTSI
             IRNUMMI(I,J)                  = IRCNTMI
             IRNUMFI(I,J)                  = IRCNTFI
             IRWATSI(I,J)                  = IRAMTSI
             IRWATMI(I,J)                  = IRAMTMI
             IRWATFI(I,J)                  = IRAMTFI
             IRSIVOL(I,J)                  = IRSIVOL(I,J)+(IRSIRATE*1000.0)
             IRMIVOL(I,J)                  = IRMIVOL(I,J)+(IRMIRATE*1000.0)
             IRFIVOL(I,J)                  = IRFIVOL(I,J)+(IRFIRATE*1000.0)
             IRELOSS(I,J)                  = IRELOSS(I,J)+(EIRR*DT) ! mm
             IRRSPLH(I,J)                  = IRRSPLH(I,J)+(FIRR*DT) ! Joules/m^2

             if(iopt_crop == 2) then   ! gecros crop model

               !*** Check for harvest
               if ((gecros1d(1) >= gecros_ds1).and.(gecros1d(42) < 0)) then
                 if (checkIfHarvest(gecros1d, DT, gecros_ds1, gecros_ds2, gecros_ds1x, &
                     gecros_ds2x) == 1) then

                   call gecros_reinit(gecros1d)
                 endif
               endif

               gecros_state (i,1:60,j)     = gecros1d(1:60)
             end if

          ENDIF                                                         ! endif of land-sea test

      ENDDO ILOOP                                                       ! of I loop
   ENDDO JLOOP                                                          ! of J loop

!------------------------------------------------------
  END SUBROUTINE noahmplsm
!------------------------------------------------------

SUBROUTINE TRANSFER_MP_PARAMETERS(VEGTYPE,SOILTYPE,SLOPETYPE,SOILCOLOR,CROPTYPE,parameters)

  USE NOAHMP_TABLES
  USE MODULE_SF_NOAHMPLSM

  implicit none

  INTEGER, INTENT(IN)    :: VEGTYPE
  INTEGER, INTENT(IN)    :: SOILTYPE(4)
  INTEGER, INTENT(IN)    :: SLOPETYPE
  INTEGER, INTENT(IN)    :: SOILCOLOR
  INTEGER, INTENT(IN)    :: CROPTYPE

  type (noahmp_parameters), intent(inout) :: parameters

  REAL    :: REFDK
  REAL    :: REFKDT
  REAL    :: FRZK
  REAL    :: FRZFACT
  INTEGER :: ISOIL

  parameters%ISWATER   =   ISWATER_TABLE
  parameters%ISBARREN  =  ISBARREN_TABLE
  parameters%ISICE     =     ISICE_TABLE
  parameters%ISCROP    =    ISCROP_TABLE
  parameters%EBLFOREST = EBLFOREST_TABLE

  parameters%URBAN_FLAG = .FALSE.
  IF( VEGTYPE == ISURBAN_TABLE    .or. VEGTYPE == LCZ_1_TABLE .or. VEGTYPE == LCZ_2_TABLE .or. &
             VEGTYPE == LCZ_3_TABLE      .or. VEGTYPE == LCZ_4_TABLE .or. VEGTYPE == LCZ_5_TABLE .or. &
             VEGTYPE == LCZ_6_TABLE      .or. VEGTYPE == LCZ_7_TABLE .or. VEGTYPE == LCZ_8_TABLE .or. &
             VEGTYPE == LCZ_9_TABLE      .or. VEGTYPE == LCZ_10_TABLE .or. VEGTYPE == LCZ_11_TABLE ) THEN
      parameters%URBAN_FLAG = .TRUE.
  ENDIF
  IF(VEGTYPE == 41 .or. VEGTYPE == 42 .or. VEGTYPE == 44 .or. VEGTYPE == 45) THEN
     parameters%URBAN_FLAG = .TRUE.
  ENDIF

!------------------------------------------------------------------------------------------!
! Transfer veg parameters
!------------------------------------------------------------------------------------------!
  parameters%CH2OP  =  CH2OP_TABLE(VEGTYPE)       !maximum intercepted h2o per unit lai+sai (mm)
  parameters%DLEAF  =  DLEAF_TABLE(VEGTYPE)       !characteristic leaf dimension (m)
  parameters%Z0MVT  =  Z0MVT_TABLE(VEGTYPE)       !momentum roughness length (m)
  parameters%HVT    =    HVT_TABLE(VEGTYPE)       !top of canopy (m)
  parameters%HVB    =    HVB_TABLE(VEGTYPE)       !bottom of canopy (m)
  parameters%DEN    =    DEN_TABLE(VEGTYPE)       !tree density (no. of trunks per m2)
  parameters%RC     =     RC_TABLE(VEGTYPE)       !tree crown radius (m)
  parameters%MFSNO  =  MFSNO_TABLE(VEGTYPE)       !snowmelt m parameter ()
  parameters%SCFFAC = SCFFAC_TABLE(VEGTYPE)       !snow cover factor (m) (originally hard-coded 2.5*z0 in SCF formulation)
  parameters%SAIM   =   SAIM_TABLE(VEGTYPE,:)     !monthly stem area index, one-sided
  parameters%LAIM   =   LAIM_TABLE(VEGTYPE,:)     !monthly leaf area index, one-sided
  parameters%SLA    =    SLA_TABLE(VEGTYPE)       !single-side leaf area per Kg [m2/kg]
  parameters%DILEFC = DILEFC_TABLE(VEGTYPE)       !coeficient for leaf stress death [1/s]
  parameters%DILEFW = DILEFW_TABLE(VEGTYPE)       !coeficient for leaf stress death [1/s]
  parameters%FRAGR  =  FRAGR_TABLE(VEGTYPE)       !fraction of growth respiration  !original was 0.3
  parameters%LTOVRC = LTOVRC_TABLE(VEGTYPE)       !leaf turnover [1/s]
  parameters%C3PSN  =  C3PSN_TABLE(VEGTYPE)       !photosynthetic pathway: 0. = c4, 1. = c3
  parameters%KC25   =   KC25_TABLE(VEGTYPE)       !co2 michaelis-menten constant at 25c (pa)
  parameters%AKC    =    AKC_TABLE(VEGTYPE)       !q10 for kc25
  parameters%KO25   =   KO25_TABLE(VEGTYPE)       !o2 michaelis-menten constant at 25c (pa)
  parameters%AKO    =    AKO_TABLE(VEGTYPE)       !q10 for ko25
  parameters%VCMX25 = VCMX25_TABLE(VEGTYPE)       !maximum rate of carboxylation at 25c (umol co2/m**2/s)
  parameters%AVCMX  =  AVCMX_TABLE(VEGTYPE)       !q10 for vcmx25
  parameters%BP     =     BP_TABLE(VEGTYPE)       !minimum leaf conductance (umol/m**2/s)
  parameters%MP     =     MP_TABLE(VEGTYPE)       !slope of conductance-to-photosynthesis relationship
  parameters%QE25   =   QE25_TABLE(VEGTYPE)       !quantum efficiency at 25c (umol co2 / umol photon)
  parameters%AQE    =    AQE_TABLE(VEGTYPE)       !q10 for qe25
  parameters%RMF25  =  RMF25_TABLE(VEGTYPE)       !leaf maintenance respiration at 25c (umol co2/m**2/s)
  parameters%RMS25  =  RMS25_TABLE(VEGTYPE)       !stem maintenance respiration at 25c (umol co2/kg bio/s)
  parameters%RMR25  =  RMR25_TABLE(VEGTYPE)       !root maintenance respiration at 25c (umol co2/kg bio/s)
  parameters%ARM    =    ARM_TABLE(VEGTYPE)       !q10 for maintenance respiration
  parameters%FOLNMX = FOLNMX_TABLE(VEGTYPE)       !foliage nitrogen concentration when f(n)=1 (%)
  parameters%TMIN   =   TMIN_TABLE(VEGTYPE)       !minimum temperature for photosynthesis (k)
  parameters%XL     =     XL_TABLE(VEGTYPE)       !leaf/stem orientation index
  parameters%RHOL   =   RHOL_TABLE(VEGTYPE,:)     !leaf reflectance: 1=vis, 2=nir
  parameters%RHOS   =   RHOS_TABLE(VEGTYPE,:)     !stem reflectance: 1=vis, 2=nir
  parameters%TAUL   =   TAUL_TABLE(VEGTYPE,:)     !leaf transmittance: 1=vis, 2=nir
  parameters%TAUS   =   TAUS_TABLE(VEGTYPE,:)     !stem transmittance: 1=vis, 2=nir
  parameters%MRP    =    MRP_TABLE(VEGTYPE)       !microbial respiration parameter (umol co2 /kg c/ s)
  parameters%CWPVT  =  CWPVT_TABLE(VEGTYPE)       !empirical canopy wind parameter
  parameters%WRRAT  =  WRRAT_TABLE(VEGTYPE)       !wood to non-wood ratio
  parameters%WDPOOL = WDPOOL_TABLE(VEGTYPE)       !wood pool (switch 1 or 0) depending on woody or not [-]
  parameters%TDLEF  =  TDLEF_TABLE(VEGTYPE)       !characteristic T for leaf freezing [K]
  parameters%NROOT  =  NROOT_TABLE(VEGTYPE)       !number of soil layers with root present
  parameters%RGL    =    RGL_TABLE(VEGTYPE)       !Parameter used in radiation stress function
  parameters%RSMIN  =     RS_TABLE(VEGTYPE)       !Minimum stomatal resistance [s m-1]
  parameters%HS     =     HS_TABLE(VEGTYPE)       !Parameter used in vapor pressure deficit function
  parameters%TOPT   =   TOPT_TABLE(VEGTYPE)       !Optimum transpiration air temperature [K]
  parameters%RSMAX  =  RSMAX_TABLE(VEGTYPE)       !Maximal stomatal resistance [s m-1]
!------------------------------------------------------------------------------------------!
! Transfer rad parameters
!------------------------------------------------------------------------------------------!
   parameters%ALBSAT    = ALBSAT_TABLE(SOILCOLOR,:)
   parameters%ALBDRY    = ALBDRY_TABLE(SOILCOLOR,:)
   parameters%ALBICE    = ALBICE_TABLE
   parameters%ALBLAK    = ALBLAK_TABLE
   parameters%OMEGAS    = OMEGAS_TABLE
   parameters%BETADS    = BETADS_TABLE
   parameters%BETAIS    = BETAIS_TABLE
   parameters%EG        = EG_TABLE

!------------------------------------------------------------------------------------------!
! Transfer crop parameters
!------------------------------------------------------------------------------------------!

  IF(CROPTYPE > 0) THEN
   parameters%PLTDAY    =    PLTDAY_TABLE(CROPTYPE)    ! Planting date
   parameters%HSDAY     =     HSDAY_TABLE(CROPTYPE)    ! Harvest date
   parameters%PLANTPOP  =  PLANTPOP_TABLE(CROPTYPE)    ! Plant density [per ha] - used?
   parameters%IRRI      =      IRRI_TABLE(CROPTYPE)    ! Irrigation strategy 0= non-irrigation 1=irrigation (no water-stress)
   parameters%GDDTBASE  =  GDDTBASE_TABLE(CROPTYPE)    ! Base temperature for GDD accumulation [C]
   parameters%GDDTCUT   =   GDDTCUT_TABLE(CROPTYPE)    ! Upper temperature for GDD accumulation [C]
   parameters%GDDS1     =     GDDS1_TABLE(CROPTYPE)    ! GDD from seeding to emergence
   parameters%GDDS2     =     GDDS2_TABLE(CROPTYPE)    ! GDD from seeding to initial vegetative
   parameters%GDDS3     =     GDDS3_TABLE(CROPTYPE)    ! GDD from seeding to post vegetative
   parameters%GDDS4     =     GDDS4_TABLE(CROPTYPE)    ! GDD from seeding to intial reproductive
   parameters%GDDS5     =     GDDS5_TABLE(CROPTYPE)    ! GDD from seeding to pysical maturity
   parameters%C3PSN     =     C3PSNI_TABLE(CROPTYPE)   ! parameters from stomata ! Zhe Zhang 2020-07-13
   parameters%KC25      =      KC25I_TABLE(CROPTYPE)
   parameters%AKC       =       AKCI_TABLE(CROPTYPE)
   parameters%KO25      =      KO25I_TABLE(CROPTYPE)
   parameters%AKO       =       AKOI_TABLE(CROPTYPE)
   parameters%AVCMX     =     AVCMXI_TABLE(CROPTYPE)
   parameters%VCMX25    =    VCMX25I_TABLE(CROPTYPE)
   parameters%BP        =        BPI_TABLE(CROPTYPE)
   parameters%MP        =        MPI_TABLE(CROPTYPE)
   parameters%FOLNMX    =    FOLNMXI_TABLE(CROPTYPE)
   parameters%QE25      =      QE25I_TABLE(CROPTYPE)   ! ends here
   parameters%C3C4      =      C3C4_TABLE(CROPTYPE)    ! photosynthetic pathway:  1. = c3 2. = c4
   parameters%AREF      =      AREF_TABLE(CROPTYPE)    ! reference maximum CO2 assimulation rate
   parameters%PSNRF     =     PSNRF_TABLE(CROPTYPE)    ! CO2 assimulation reduction factor(0-1) (caused by non-modeling part,e.g.pest,weeds)
   parameters%I2PAR     =     I2PAR_TABLE(CROPTYPE)    ! Fraction of incoming solar radiation to photosynthetically active radiation
   parameters%TASSIM0   =   TASSIM0_TABLE(CROPTYPE)    ! Minimum temperature for CO2 assimulation [C]
   parameters%TASSIM1   =   TASSIM1_TABLE(CROPTYPE)    ! CO2 assimulation linearly increasing until temperature reaches T1 [C]
   parameters%TASSIM2   =   TASSIM2_TABLE(CROPTYPE)    ! CO2 assmilation rate remain at Aref until temperature reaches T2 [C]
   parameters%K         =         K_TABLE(CROPTYPE)    ! light extinction coefficient
   parameters%EPSI      =      EPSI_TABLE(CROPTYPE)    ! initial light use efficiency
   parameters%Q10MR     =     Q10MR_TABLE(CROPTYPE)    ! q10 for maintainance respiration
   parameters%FOLN_MX   =   FOLN_MX_TABLE(CROPTYPE)    ! foliage nitrogen concentration when f(n)=1 (%)
   parameters%LEFREEZ   =   LEFREEZ_TABLE(CROPTYPE)    ! characteristic T for leaf freezing [K]
   parameters%DILE_FC   =   DILE_FC_TABLE(CROPTYPE,:)  ! coeficient for temperature leaf stress death [1/s]
   parameters%DILE_FW   =   DILE_FW_TABLE(CROPTYPE,:)  ! coeficient for water leaf stress death [1/s]
   parameters%FRA_GR    =    FRA_GR_TABLE(CROPTYPE)    ! fraction of growth respiration
   parameters%LF_OVRC   =   LF_OVRC_TABLE(CROPTYPE,:)  ! fraction of leaf turnover  [1/s]
   parameters%ST_OVRC   =   ST_OVRC_TABLE(CROPTYPE,:)  ! fraction of stem turnover  [1/s]
   parameters%RT_OVRC   =   RT_OVRC_TABLE(CROPTYPE,:)  ! fraction of root tunrover  [1/s]
   parameters%LFMR25    =    LFMR25_TABLE(CROPTYPE)    ! leaf maintenance respiration at 25C [umol CO2/m**2  /s]
   parameters%STMR25    =    STMR25_TABLE(CROPTYPE)    ! stem maintenance respiration at 25C [umol CO2/kg bio/s]
   parameters%RTMR25    =    RTMR25_TABLE(CROPTYPE)    ! root maintenance respiration at 25C [umol CO2/kg bio/s]
   parameters%GRAINMR25 = GRAINMR25_TABLE(CROPTYPE)    ! grain maintenance respiration at 25C [umol CO2/kg bio/s]
   parameters%LFPT      =      LFPT_TABLE(CROPTYPE,:)  ! fraction of carbohydrate flux to leaf
   parameters%STPT      =      STPT_TABLE(CROPTYPE,:)  ! fraction of carbohydrate flux to stem
   parameters%RTPT      =      RTPT_TABLE(CROPTYPE,:)  ! fraction of carbohydrate flux to root
   parameters%GRAINPT   =   GRAINPT_TABLE(CROPTYPE,:)  ! fraction of carbohydrate flux to grain
   parameters%LFCT      =      LFCT_TABLE(CROPTYPE,:)  ! fraction of translocation to grain ! Zhe Zhang 2020-07-13
   parameters%STCT      =      STCT_TABLE(CROPTYPE,:)  ! fraction of translocation to grain
   parameters%RTCT      =      RTCT_TABLE(CROPTYPE,:)  ! fraction of translocation to grain
   parameters%BIO2LAI   =   BIO2LAI_TABLE(CROPTYPE)    ! leaf are per living leaf biomass [m^2/kg]
  END IF

!------------------------------------------------------------------------------------------!
! Transfer global parameters
!------------------------------------------------------------------------------------------!
   parameters%CO2        =         CO2_TABLE
   parameters%O2         =          O2_TABLE
   parameters%TIMEAN     =      TIMEAN_TABLE
   parameters%FSATMX     =      FSATMX_TABLE
   parameters%Z0SNO      =       Z0SNO_TABLE
   parameters%SSI        =         SSI_TABLE
   parameters%SNOW_RET_FAC = SNOW_RET_FAC_TABLE
   parameters%SNOW_EMIS  =   SNOW_EMIS_TABLE
   parameters%SWEMX        =     SWEMX_TABLE
   parameters%TAU0         =      TAU0_TABLE
   parameters%GRAIN_GROWTH = GRAIN_GROWTH_TABLE
   parameters%EXTRA_GROWTH = EXTRA_GROWTH_TABLE
   parameters%DIRT_SOOT    =    DIRT_SOOT_TABLE
   parameters%BATS_COSZ    =    BATS_COSZ_TABLE
   parameters%BATS_VIS_NEW = BATS_VIS_NEW_TABLE
   parameters%BATS_NIR_NEW = BATS_NIR_NEW_TABLE
   parameters%BATS_VIS_AGE = BATS_VIS_AGE_TABLE
   parameters%BATS_NIR_AGE = BATS_NIR_AGE_TABLE
   parameters%BATS_VIS_DIR = BATS_VIS_DIR_TABLE
   parameters%BATS_NIR_DIR = BATS_NIR_DIR_TABLE
   parameters%RSURF_SNOW =  RSURF_SNOW_TABLE
   parameters%RSURF_EXP  =   RSURF_EXP_TABLE

! ----------------------------------------------------------------------
!  Transfer soil parameters
! ----------------------------------------------------------------------
! Adjust the parameter definitions associated with the soil IF the vegetation type
! is equal to permeable pavements. This soil type is based on the baseline soil
! type within the WRF SOIL.TBL, but has a higher saturated hydraulic conductivity

IF(VEGTYPE.EQ.44) THEN !Permeable Pavement Added by Aaron A.
    do isoil = 1, 2
      parameters%BEXP(isoil)   = BEXP_TABLE   (1)
      parameters%DKSAT(isoil)  = 1.0E-4
      parameters%DWSAT(isoil)  = DWSAT_TABLE  (1)
      parameters%PSISAT(isoil) = PSISAT_TABLE (1)
      parameters%QUARTZ(isoil) = QUARTZ_TABLE (1)
      parameters%SMCDRY(isoil) = SMCDRY_TABLE (1)
      parameters%SMCMAX(isoil) = 0.4
      parameters%SMCREF(isoil) = SMCREF_TABLE (1)
      parameters%SMCWLT(isoil) = SMCWLT_TABLE (1)
    end do
 ELSE IF(VEGTYPE.EQ.46) THEN ! Green Roof Added by Aaron A.
    do isoil = 1, size(soiltype)
      parameters%BEXP(isoil)   = BEXP_TABLE   (4)
      parameters%DKSAT(isoil)  = DKSAT_TABLE  (4)
      parameters%DWSAT(isoil)  = DWSAT_TABLE  (4)
      parameters%PSISAT(isoil) = PSISAT_TABLE (4)
      parameters%QUARTZ(isoil) = QUARTZ_TABLE (5)
      parameters%SMCDRY(isoil) = 0.1
      parameters%SMCMAX(isoil) = 0.4
      parameters%SMCREF(isoil) = 0.335
      parameters%SMCWLT(isoil) = 0.1

    end do
ELSE
    do isoil = 1, size(soiltype)
      parameters%BEXP(isoil)   = BEXP_TABLE   (SOILTYPE(isoil))
      parameters%DKSAT(isoil)  = DKSAT_TABLE  (SOILTYPE(isoil))
      parameters%DWSAT(isoil)  = DWSAT_TABLE  (SOILTYPE(isoil))
      parameters%PSISAT(isoil) = PSISAT_TABLE (SOILTYPE(isoil))
      parameters%QUARTZ(isoil) = QUARTZ_TABLE (SOILTYPE(isoil))
      parameters%SMCDRY(isoil) = SMCDRY_TABLE (SOILTYPE(isoil))
      parameters%SMCMAX(isoil) = SMCMAX_TABLE (SOILTYPE(isoil))
      parameters%SMCREF(isoil) = SMCREF_TABLE (SOILTYPE(isoil))
      parameters%SMCWLT(isoil) = SMCWLT_TABLE (SOILTYPE(isoil))
    end do
ENDIF


    parameters%F1     = F1_TABLE(SOILTYPE(1))
    parameters%REFDK  = REFDK_TABLE
    parameters%REFKDT = REFKDT_TABLE
    parameters%BVIC   = BVIC_TABLE(SOILTYPE(1))
    parameters%AXAJ   = AXAJ_TABLE(SOILTYPE(1))
    parameters%BXAJ   = BXAJ_TABLE(SOILTYPE(1))
    parameters%XXAJ   = XXAJ_TABLE(SOILTYPE(1))
    parameters%BDVIC  = BDVIC_TABLE(SOILTYPE(1))
    parameters%GDVIC  = GDVIC_TABLE(SOILTYPE(1))
    parameters%BBVIC  = BBVIC_TABLE(SOILTYPE(1))

!------------------------------------------------------------------------------------------!
! Transfer irrigation parameters
!------------------------------------------------------------------------------------------!
    parameters%IRR_FRAC   = IRR_FRAC_TABLE      ! irrigation Fraction
    parameters%IRR_HAR    = IRR_HAR_TABLE       ! number of days before harvest date to stop irrigation
    parameters%IRR_LAI    = IRR_LAI_TABLE       ! minimum lai to trigger irrigation
    parameters%IRR_MAD    = IRR_MAD_TABLE       ! management allowable deficit (0-1)
    parameters%FILOSS     = FILOSS_TABLE        ! fraction of flood irrigation loss (0-1)
    parameters%SPRIR_RATE = SPRIR_RATE_TABLE    ! mm/h, sprinkler irrigation rate
    parameters%MICIR_RATE = MICIR_RATE_TABLE    ! mm/h, micro irrigation rate
    parameters%FIRTFAC    = FIRTFAC_TABLE       ! flood application rate factor
    parameters%IR_RAIN    = IR_RAIN_TABLE       ! maximum precipitation to stop irrigation trigger

!------------------------------------------------------------------------------------------!
! Transfer tiledrain parameters
!------------------------------------------------------------------------------------------!
   parameters%KLAT_FAC        = KLAT_FAC_TABLE(SOILTYPE(1))
   parameters%TDSMC_FAC       = TDSMCFAC_TABLE(SOILTYPE(1))
   parameters%TD_DC           = TD_DC_TABLE(SOILTYPE(1))
   parameters%TD_DCOEF        = TD_DCOEF_TABLE(SOILTYPE(1))
   parameters%TD_RADI         = TD_RADI_TABLE(SOILTYPE(1))
   parameters%TD_SPAC         = TD_SPAC_TABLE(SOILTYPE(1))
   parameters%TD_DDRAIN       = TD_DDRAIN_TABLE(SOILTYPE(1))
   parameters%TD_DEPTH        = TD_DEPTH_TABLE(SOILTYPE(1))
   parameters%TD_ADEPTH       = TD_ADEPTH_TABLE(SOILTYPE(1))
   parameters%DRAIN_LAYER_OPT = DRAIN_LAYER_OPT_TABLE
   parameters%TD_D            = TD_D_TABLE(SOILTYPE(1))

! ----------------------------------------------------------------------
! Transfer GENPARM parameters
! ----------------------------------------------------------------------
    parameters%CSOIL  = CSOIL_TABLE
    parameters%ZBOT   = ZBOT_TABLE
    parameters%CZIL   = CZIL_TABLE

    FRZK   = FRZK_TABLE
    parameters%KDT    = parameters%REFKDT * parameters%DKSAT(1) / parameters%REFDK
    parameters%SLOPE  = SLOPE_TABLE(SLOPETYPE)

    IF((parameters%URBAN_FLAG).or.(VEGTYPE.EQ.41))THEN  ! Hardcoding some urban parameters for soil
       !parameters%SMCMAX = 0.45
       !parameters%SMCREF = 0.42
       !parameters%SMCWLT = 0.40
       !parameters%SMCDRY = 0.40
       parameters%CSOIL  = 3.E6
    ENDIF

! adjust FRZK parameter to actual soil type: FRZK * FRZFACT

    IF(SOILTYPE(1) /= 14) then
      FRZFACT = (parameters%SMCMAX(1) / parameters%SMCREF(1)) * (0.412 / 0.468)
      parameters%FRZX = FRZK * FRZFACT
    END IF

 END SUBROUTINE TRANSFER_MP_PARAMETERS

SUBROUTINE PEDOTRANSFER_SR2006(nsoil,sand,clay,orgm,parameters)

  use module_sf_noahmplsm
  use noahmp_tables

  implicit none

  integer,                    intent(in   ) :: nsoil     ! number of soil layers
  real, dimension( 1:nsoil ), intent(inout) :: sand
  real, dimension( 1:nsoil ), intent(inout) :: clay
  real, dimension( 1:nsoil ), intent(inout) :: orgm

  real, dimension( 1:nsoil ) :: theta_1500t
  real, dimension( 1:nsoil ) :: theta_1500
  real, dimension( 1:nsoil ) :: theta_33t
  real, dimension( 1:nsoil ) :: theta_33
  real, dimension( 1:nsoil ) :: theta_s33t
  real, dimension( 1:nsoil ) :: theta_s33
  real, dimension( 1:nsoil ) :: psi_et
  real, dimension( 1:nsoil ) :: psi_e

  type(noahmp_parameters), intent(inout) :: parameters
  integer :: k

  do k = 1,4
    if(sand(k) <= 0 .or. clay(k) <= 0) then
      sand(k) = 0.41
      clay(k) = 0.18
    end if
    if(orgm(k) <= 0 ) orgm(k) = 0.0
  end do

  theta_1500t =   sr2006_theta_1500t_a*sand       &
                + sr2006_theta_1500t_b*clay       &
                + sr2006_theta_1500t_c*orgm       &
                + sr2006_theta_1500t_d*sand*orgm  &
                + sr2006_theta_1500t_e*clay*orgm  &
                + sr2006_theta_1500t_f*sand*clay  &
                + sr2006_theta_1500t_g

  theta_1500  =   theta_1500t                      &
                + sr2006_theta_1500_a*theta_1500t  &
                + sr2006_theta_1500_b

  theta_33t   =   sr2006_theta_33t_a*sand       &
                + sr2006_theta_33t_b*clay       &
                + sr2006_theta_33t_c*orgm       &
                + sr2006_theta_33t_d*sand*orgm  &
                + sr2006_theta_33t_e*clay*orgm  &
                + sr2006_theta_33t_f*sand*clay  &
                + sr2006_theta_33t_g

  theta_33    =   theta_33t                              &
                + sr2006_theta_33_a*theta_33t*theta_33t  &
                + sr2006_theta_33_b*theta_33t            &
                + sr2006_theta_33_c

  theta_s33t  =   sr2006_theta_s33t_a*sand      &
                + sr2006_theta_s33t_b*clay      &
                + sr2006_theta_s33t_c*orgm      &
                + sr2006_theta_s33t_d*sand*orgm &
                + sr2006_theta_s33t_e*clay*orgm &
                + sr2006_theta_s33t_f*sand*clay &
                + sr2006_theta_s33t_g

  theta_s33   = theta_s33t                       &
                + sr2006_theta_s33_a*theta_s33t  &
                + sr2006_theta_s33_b

  psi_et      =   sr2006_psi_et_a*sand           &
                + sr2006_psi_et_b*clay           &
                + sr2006_psi_et_c*theta_s33      &
                + sr2006_psi_et_d*sand*theta_s33 &
                + sr2006_psi_et_e*clay*theta_s33 &
                + sr2006_psi_et_f*sand*clay      &
                + sr2006_psi_et_g

  psi_e       =   psi_et                        &
                + sr2006_psi_e_a*psi_et*psi_et  &
                + sr2006_psi_e_b*psi_et         &
                + sr2006_psi_e_c

  parameters%smcwlt = theta_1500
  parameters%smcref = theta_33
  parameters%smcmax =   theta_33    &
                      + theta_s33            &
                      + sr2006_smcmax_a*sand &
                      + sr2006_smcmax_b

  parameters%bexp   = 3.816712826 / (log(theta_33) - log(theta_1500) )
  parameters%psisat = psi_e
  parameters%dksat  = 1930.0 * (parameters%smcmax - theta_33) ** (3.0 - 1.0/parameters%bexp)
  parameters%quartz = sand

! Units conversion

  parameters%psisat = max(0.1,parameters%psisat)     ! arbitrarily impose a limit of 0.1kpa
  parameters%psisat = 0.101997 * parameters%psisat   ! convert kpa to m
  parameters%dksat  = parameters%dksat / 3600000.0   ! convert mm/h to m/s
  parameters%dwsat  = parameters%dksat * parameters%psisat *parameters%bexp / parameters%smcmax  ! units should be m*m/s
  parameters%smcdry = parameters%smcwlt

! Introducing somewhat arbitrary limits (based on SOILPARM) to prevent bad things

  parameters%smcmax = max(0.32 ,min(parameters%smcmax,             0.50 ))
  parameters%smcref = max(0.17 ,min(parameters%smcref,parameters%smcmax ))
  parameters%smcwlt = max(0.01 ,min(parameters%smcwlt,parameters%smcref ))
  parameters%smcdry = max(0.01 ,min(parameters%smcdry,parameters%smcref ))
  parameters%bexp   = max(2.50 ,min(parameters%bexp,               12.0 ))
  parameters%psisat = max(0.03 ,min(parameters%psisat,             1.00 ))
  parameters%dksat  = max(5.e-7,min(parameters%dksat,              1.e-5))
  parameters%dwsat  = max(1.e-6,min(parameters%dwsat,              3.e-5))
  parameters%quartz = max(0.05 ,min(parameters%quartz,             0.95 ))

 END SUBROUTINE PEDOTRANSFER_SR2006

  SUBROUTINE NOAHMP_INIT ( MMINLU, SNOW , SNOWH , CANWAT , ISLTYP ,   IVGTYP, XLAT, &
       TSLB , SMOIS , SH2O , DZS , FNDSOILW , FNDSNOWH ,             &
       TSK, isnowxy , tvxy     ,tgxy     ,canicexy ,         TMN,     XICE,   &
       canliqxy ,eahxy    ,tahxy    ,cmxy     ,chxy     ,                     &
       fwetxy   ,sneqvoxy ,alboldxy ,qsnowxy, qrainxy, wslakexy, zwtxy, waxy, &
       wtxy     ,tsnoxy   ,zsnsoxy  ,snicexy  ,snliqxy  ,lfmassxy ,rtmassxy , &
       stmassxy ,woodxy   ,stblcpxy ,fastcpxy ,xsaixy   ,lai      ,           &
       grainxy  ,gddxy    ,                                                   &
       croptype ,cropcat  ,                      &
       irnumsi  ,irnummi  ,irnumfi  ,irwatsi,    &
       irwatmi  ,irwatfi  ,ireloss  ,irsivol,    &
       irmivol  ,irfivol  ,irrsplh  ,            &
!jref:start
       t2mvxy   ,t2mbxy   ,chstarxy,             &
!jref:end
       NSOIL, restart,                 &
       allowed_to_read , iopt_run,  iopt_crop, iopt_irr, iopt_irrm,           &
       sf_urban_physics,                         &  ! urban scheme
       ids,ide, jds,jde, kds,kde,                &
       ims,ime, jms,jme, kms,kme,                &
       its,ite, jts,jte, kts,kte,                &
       smoiseq  ,smcwtdxy ,rechxy   ,deeprechxy, qtdrain, areaxy, dx, dy, msftx, msfty,& ! Optional groundwater
       wtddt    ,stepwtd  ,dt       ,qrfsxy     ,qspringsxy  , qslatxy    ,  &      ! Optional groundwater
       fdepthxy ,ht     ,riverbedxy ,eqzwt     ,rivercondxy ,pexpxy       ,  &      ! Optional groundwater
       rechclim,                                                             &      ! Optional groundwater
       gecros_state)                                                                ! Optional gecros crop

  USE NOAHMP_TABLES
  use module_sf_gecros, only: seednc,sla0,slnmin,ffat,flig,foac,fmin,npl,seedw,eg,fcrsh,seednc,lnci,cfv


  IMPLICIT NONE

! Initializing Canopy air temperature to 287 K seems dangerous to me [KWM].

    INTEGER, INTENT(IN   )    ::     ids,ide, jds,jde, kds,kde,  &
         &                           ims,ime, jms,jme, kms,kme,  &
         &                           its,ite, jts,jte, kts,kte

    INTEGER, INTENT(IN)       ::     NSOIL, iopt_run, iopt_crop, iopt_irr, iopt_irrm

    LOGICAL, INTENT(IN)       ::     restart,                    &
         &                           allowed_to_read
    INTEGER, INTENT(IN)       ::     sf_urban_physics                              ! urban, by yizhou

    REAL,    DIMENSION( NSOIL), INTENT(IN)    ::     DZS  ! Thickness of the soil layers [m]
    REAL,    INTENT(IN) , OPTIONAL ::     DX, DY
    REAL,    DIMENSION( ims:ime, jms:jme ) ,  INTENT(IN) , OPTIONAL :: MSFTX,MSFTY

    REAL,    DIMENSION( ims:ime, NSOIL, jms:jme ) ,    &
         &   INTENT(INOUT)    ::     SMOIS,                      &
         &                           SH2O,                       &
         &                           TSLB

    REAL,    DIMENSION( ims:ime, jms:jme ) ,                     &
         &   INTENT(INOUT)    ::     SNOW,                       &
         &                           SNOWH,                      &
         &                           CANWAT

    INTEGER, DIMENSION( ims:ime, jms:jme ),                      &
         &   INTENT(IN)       ::     ISLTYP,  &
                                     IVGTYP

    LOGICAL, INTENT(IN)       ::     FNDSOILW,                   &
         &                           FNDSNOWH

    REAL, DIMENSION(ims:ime,jms:jme), INTENT(IN) :: XLAT         !latitude
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(IN) :: TSK         !skin temperature (k)
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: TMN         !deep soil temperature (k)
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(IN) :: XICE         !sea ice fraction
    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: isnowxy     !actual no. of snow layers
    REAL, DIMENSION(ims:ime,-2:NSOIL,jms:jme), INTENT(INOUT) :: zsnsoxy  !snow layer depth [m]
    REAL, DIMENSION(ims:ime,-2:              0,jms:jme), INTENT(INOUT) :: tsnoxy   !snow temperature [K]
    REAL, DIMENSION(ims:ime,-2:              0,jms:jme), INTENT(INOUT) :: snicexy  !snow layer ice [mm]
    REAL, DIMENSION(ims:ime,-2:              0,jms:jme), INTENT(INOUT) :: snliqxy  !snow layer liquid water [mm]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: tvxy        !vegetation canopy temperature
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: tgxy        !ground surface temperature
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: canicexy    !canopy-intercepted ice (mm)
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: canliqxy    !canopy-intercepted liquid water (mm)
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: eahxy       !canopy air vapor pressure (pa)
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: tahxy       !canopy air temperature (k)
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: cmxy        !momentum drag coefficient
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: chxy        !sensible heat exchange coefficient
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: fwetxy      !wetted or snowed fraction of the canopy (-)
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: sneqvoxy    !snow mass at last time step(mm h2o)
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: alboldxy    !snow albedo at last time step (-)
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: qsnowxy     !snowfall on the ground [mm/s]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: qrainxy     !rainfall on the ground [mm/s]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: wslakexy    !lake water storage [mm]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: zwtxy       !water table depth [m]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: waxy        !water in the "aquifer" [mm]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: wtxy        !groundwater storage [mm]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: lfmassxy    !leaf mass [g/m2]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: rtmassxy    !mass of fine roots [g/m2]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: stmassxy    !stem mass [g/m2]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: woodxy      !mass of wood (incl. woody roots) [g/m2]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: grainxy     !mass of grain [g/m2] !XING
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: gddxy       !growing degree days !XING
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: stblcpxy    !stable carbon in deep soil [g/m2]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: fastcpxy    !short-lived carbon, shallow soil [g/m2]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: xsaixy      !stem area index
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: lai         !leaf area index
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: qtdrain     !tile drainage (mm)

    INTEGER, DIMENSION(ims:ime,  jms:jme), INTENT(OUT) :: cropcat
    REAL   , DIMENSION(ims:ime,5,jms:jme), INTENT(IN ) :: croptype

    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: irnumsi
    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: irnummi
    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: irnumfi
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: irwatsi
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: irwatmi
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: irwatfi
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: ireloss
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: irsivol
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: irmivol
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: irfivol
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: irrsplh

! IOPT_RUN = 5 option

    REAL, DIMENSION(ims:ime,1:nsoil,jms:jme), INTENT(INOUT) , OPTIONAL :: smoiseq !equilibrium soil moisture content [m3m-3]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) , OPTIONAL :: smcwtdxy    !deep soil moisture content [m3m-3]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) , OPTIONAL :: deeprechxy  !deep recharge [m]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) , OPTIONAL :: rechxy      !accumulated recharge [mm]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) , OPTIONAL :: qrfsxy      !accumulated flux from groundwater to rivers [mm]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) , OPTIONAL :: qspringsxy  !accumulated seeping water [mm]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) , OPTIONAL :: qslatxy     !accumulated lateral flow [mm]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) , OPTIONAL :: areaxy      !grid cell area [m2]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(IN) , OPTIONAL :: FDEPTHXY    !efolding depth for transmissivity (m)
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(IN) , OPTIONAL :: HT          !terrain height (m)
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) , OPTIONAL :: RIVERBEDXY  !riverbed depth (m)
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) , OPTIONAL :: EQZWT       !equilibrium water table depth (m)
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT), OPTIONAL :: RIVERCONDXY !river conductance
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT), OPTIONAL :: PEXPXY      !factor for river conductance
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(IN) , OPTIONAL :: rechclim

    REAL, DIMENSION(ims:ime,60,jms:jme), INTENT(INOUT),   OPTIONAL :: gecros_state                                     ! Optional gecros crop

    INTEGER,  INTENT(OUT) , OPTIONAL :: STEPWTD
    REAL, INTENT(IN) , OPTIONAL :: DT, WTDDT

!jref:start
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: t2mvxy        !2m temperature vegetation part (k)
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: t2mbxy        !2m temperature bare ground part (k)
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: chstarxy        !dummy
!jref:end


    REAL, DIMENSION(1:NSOIL)  :: ZSOIL      ! Depth of the soil layer bottom (m) from
    !                                                   the surface (negative)

    REAL                      :: BEXP, SMCMAX, PSISAT
    REAL                      :: FK, masslai,masssai

! gecros local variables
    REAL ::  hti,rdi,fpro,lncmin,fcar,cfo,clvi,crti,ygo,nlvi,laii,nrti,slnbi


    REAL, PARAMETER           :: BLIM  = 5.5
    REAL, PARAMETER           :: HLICE = 3.335E5
    REAL, PARAMETER           :: GRAV = 9.81
    REAL, PARAMETER           :: T0 = 273.15

    INTEGER                   :: errflag, i,j,itf,jtf,ns

    character(len=240) :: err_message
    character(len=4)  :: MMINSL
    character(len=*), intent(in) :: MMINLU
    MMINSL='STAS'

    call read_mp_veg_parameters(trim(MMINLU))
    call read_mp_soil_parameters()
    call read_mp_rad_parameters()
    call read_mp_global_parameters()
    call read_mp_crop_parameters()
    call read_tiledrain_parameters()
    call read_mp_optional_parameters()
    if(iopt_irr  >= 1) call read_mp_irrigation_parameters()

    IF( .NOT. restart ) THEN

       itf=min0(ite,ide-1)
       jtf=min0(jte,jde-1)

       !
       ! initialize physical snow height SNOWH
       !
       IF(.NOT.FNDSNOWH)THEN
          ! If no SNOWH do the following
          CALL wrf_message( 'SNOW HEIGHT NOT FOUND - VALUE DEFINED IN LSMINIT' )
          DO J = jts,jtf
             DO I = its,itf
                SNOWH(I,J)=SNOW(I,J)*0.005               ! SNOW in mm and SNOWH in m
             ENDDO
          ENDDO
       ENDIF


       ! Check if snow/snowh are consistent and cap SWE at 2000mm;
       !  the Noah-MP code does it internally but if we don't do it here, problems ensue
       DO J = jts,jtf
          DO I = its,itf
             IF ( SNOW(i,j) > 0. .AND. SNOWH(i,j) == 0. .OR. SNOWH(i,j) > 0. .AND. SNOW(i,j) == 0.) THEN
               WRITE(err_message,*)"problem with initial snow fields: snow/snowh>0 while snowh/snow=0 at i,j" &
                                     ,i,j,snow(i,j),snowh(i,j)
               CALL wrf_message(err_message)
             ENDIF
             IF ( SNOW( i,j ) > 2000.0 ) THEN
               SNOWH(I,J) = SNOWH(I,J) * 2000.0 / SNOW(I,J)      ! SNOW in mm and SNOWH in m
               SNOW (I,J) = 2000.0                               ! cap SNOW at 2000, maintain density
             ENDIF
          ENDDO
       ENDDO

       errflag = 0
       DO j = jts,jtf
          DO i = its,itf
             IF ( ISLTYP( i,j ) .LT. 1 ) THEN
                errflag = 1
                WRITE(err_message,*)"module_sf_noahlsm.F: lsminit: out of range ISLTYP ",i,j,ISLTYP( i,j )
                CALL wrf_message(err_message)
             ENDIF
          ENDDO
       ENDDO
       IF ( errflag .EQ. 1 ) THEN
          CALL wrf_error_fatal( "module_sf_noahlsm.F: lsminit: out of range value "// &
               "of ISLTYP. Is this field in the input?" )
       ENDIF
! GAC-->LATERALFLOW
! 20130219 - No longer need this - see module_data_gocart_dust
!#if ( WRF_CHEM == 1 )
!       !
!       ! need this parameter for dust parameterization in wrf/chem
!       !
!       do I=1,NSLTYPE
!          porosity(i)=maxsmc(i)
!       enddo
!#endif
! <--GAC

! initialize soil liquid water content SH2O

       DO J = jts , jtf
          DO I = its , itf
	    IF(IVGTYP(I,J)==ISICE_TABLE .AND. XICE(I,J) <= 0.0) THEN
              DO NS=1, NSOIL
	        SMOIS(I,NS,J) = 1.0                     ! glacier starts all frozen
	        SH2O(I,NS,J) = 0.0
	        TSLB(I,NS,J) = MIN(TSLB(I,NS,J),263.15) ! set glacier temp to at most -10C
              END DO
	        !TMN(I,J) = MIN(TMN(I,J),263.15)         ! set deep temp to at most -10C
		SNOW(I,J) = MAX(SNOW(I,J), 10.0)        ! set SWE to at least 10mm
                SNOWH(I,J)=SNOW(I,J)*0.01               ! SNOW in mm and SNOWH in m
	    ELSE

              BEXP   =   BEXP_TABLE(ISLTYP(I,J))
              SMCMAX = SMCMAX_TABLE(ISLTYP(I,J))
              PSISAT = PSISAT_TABLE(ISLTYP(I,J))

              DO NS=1, NSOIL
	        IF ( SMOIS(I,NS,J) > SMCMAX )  SMOIS(I,NS,J) = SMCMAX
              END DO
              IF ( ( BEXP > 0.0 ) .AND. ( SMCMAX > 0.0 ) .AND. ( PSISAT > 0.0 ) ) THEN
                DO NS=1, NSOIL
                   IF ( TSLB(I,NS,J) < 273.149 ) THEN    ! Use explicit as initial soil ice
                      FK=(( (HLICE/(GRAV*(-PSISAT))) *                              &
                           ((TSLB(I,NS,J)-T0)/TSLB(I,NS,J)) )**(-1/BEXP) )*SMCMAX
                      FK = MAX(FK, 0.02)
                      SH2O(I,NS,J) = MIN( FK, SMOIS(I,NS,J) )
                   ELSE
                      SH2O(I,NS,J)=SMOIS(I,NS,J)
                   ENDIF
                END DO
              ELSE
                DO NS=1, NSOIL
                   SH2O(I,NS,J)=SMOIS(I,NS,J)
                END DO
              ENDIF
            ENDIF
          ENDDO
       ENDDO
!  ENDIF


       DO J = jts,jtf
          DO I = its,itf
             qtdrain    (I,J) = 0.
             tvxy       (I,J) = TSK(I,J)
	       if(snow(i,j) > 0.0 .and. tsk(i,j) > 273.15) tvxy(I,J) = 273.15
             tgxy       (I,J) = TSK(I,J)
	       if(snow(i,j) > 0.0 .and. tsk(i,j) > 273.15) tgxy(I,J) = 273.15
             CANWAT     (I,J) = 0.0
             canliqxy   (I,J) = CANWAT(I,J)
             canicexy   (I,J) = 0.
             eahxy      (I,J) = 2000.
             tahxy      (I,J) = TSK(I,J)
	       if(snow(i,j) > 0.0 .and. tsk(i,j) > 273.15) tahxy(I,J) = 273.15
!             tahxy      (I,J) = 287.
!jref:start
             t2mvxy     (I,J) = TSK(I,J)
	       if(snow(i,j) > 0.0 .and. tsk(i,j) > 273.15) t2mvxy(I,J) = 273.15
             t2mbxy     (I,J) = TSK(I,J)
	       if(snow(i,j) > 0.0 .and. tsk(i,j) > 273.15) t2mbxy(I,J) = 273.15
             chstarxy     (I,J) = 0.1
!jref:end

             cmxy       (I,J) = 0.0
             chxy       (I,J) = 0.0
             fwetxy     (I,J) = 0.0
             sneqvoxy   (I,J) = 0.0
             alboldxy   (I,J) = 0.65
             qsnowxy    (I,J) = 0.0
             qrainxy    (I,J) = 0.0
             wslakexy   (I,J) = 0.0

             if(iopt_run.ne.5) then
                   waxy       (I,J) = 4900.                                       !???
                   wtxy       (I,J) = waxy(i,j)                                   !???
                   zwtxy      (I,J) = (25. + 2.0) - waxy(i,j)/1000/0.2            !???
             else
                   waxy       (I,J) = 0.
                   wtxy       (I,J) = 0.
                   areaxy     (I,J) = (DX * DY) / ( MSFTX(I,J) * MSFTY(I,J) )
             endif

           IF(IVGTYP(I,J) == ISBARREN_TABLE .OR. IVGTYP(I,J) == ISICE_TABLE .OR. &
	      ( SF_URBAN_PHYSICS == 0 .AND. IVGTYP(I,J) == ISURBAN_TABLE )  .OR. &
	      IVGTYP(I,J) == ISWATER_TABLE ) THEN

	     lai        (I,J) = 0.0
             xsaixy     (I,J) = 0.0
             lfmassxy   (I,J) = 0.0
             stmassxy   (I,J) = 0.0
             rtmassxy   (I,J) = 0.0
             woodxy     (I,J) = 0.0
             stblcpxy   (I,J) = 0.0
             fastcpxy   (I,J) = 0.0
             grainxy    (I,J) = 1E-10
             gddxy      (I,J) = 0
	     cropcat    (I,J) = 0

	   ELSE

	     lai        (I,J) = max(lai(i,j),0.05)             ! at least start with 0.05 for arbitrary initialization (v3.7)
             xsaixy     (I,J) = max(0.1*lai(I,J),0.05)         ! MB: arbitrarily initialize SAI using input LAI (v3.7)
             masslai = 1000. / max(SLA_TABLE(IVGTYP(I,J)),1.0) ! conversion from lai to mass  (v3.7)
             lfmassxy   (I,J) = lai(i,j)*masslai               ! use LAI to initialize (v3.7)
             masssai = 1000. / 3.0                             ! conversion from lai to mass (v3.7)
             stmassxy   (I,J) = xsaixy(i,j)*masssai            ! use SAI to initialize (v3.7)
             rtmassxy   (I,J) = 500.0                          ! these are all arbitrary and probably should be
             woodxy     (I,J) = 500.0                          ! in the table or read from initialization
             stblcpxy   (I,J) = 1000.0                         !
             fastcpxy   (I,J) = 1000.0                         !
             grainxy    (I,J) = 1E-10
             gddxy      (I,J) = 0

! Initialize crop for Liu crop model

	     if(iopt_crop == 1 ) then
	       cropcat    (i,j) = default_crop_table
               if(croptype(i,5,j) >= 0.5) then
                 rtmassxy(i,j) = 0.0
                 woodxy  (i,j) = 0.0

	         if(    croptype(i,1,j) > croptype(i,2,j) .and. &
		        croptype(i,1,j) > croptype(i,3,j) .and. &
		        croptype(i,1,j) > croptype(i,4,j) ) then   ! choose corn

		   cropcat (i,j) = 1
                   lfmassxy(i,j) =    lai(i,j)/0.015               ! Initialize lfmass Zhe Zhang 2020-07-13
                   stmassxy(i,j) = xsaixy(i,j)/0.003

	         elseif(croptype(i,2,j) > croptype(i,1,j) .and. &
		        croptype(i,2,j) > croptype(i,3,j) .and. &
		        croptype(i,2,j) > croptype(i,4,j) ) then   ! choose soybean

		   cropcat (i,j) = 2
                   lfmassxy(i,j) =    lai(i,j)/0.030               ! Initialize lfmass Zhe Zhang 2020-07-13
                   stmassxy(i,j) = xsaixy(i,j)/0.003

	         else

		   cropcat (i,j) = default_crop_table
                   lfmassxy(i,j) =    lai(i,j)/0.035
                   stmassxy(i,j) = xsaixy(i,j)/0.003

	         end if

	       end if
	     end if

! Initialize cropcat for gecros crop model

	     if(iopt_crop == 2) then
	       cropcat    (i,j) = 0
               if(croptype(i,5,j) >= 0.5) then
                  if(croptype(i,3,j) > 0.0)             cropcat(i,j) = 1 ! if any wheat, set to wheat
                  if(croptype(i,1,j) > croptype(i,3,j)) cropcat(i,j) = 2 ! change to maize
	       end if

               hti    = 0.01
               rdi    = 10.
               fpro   = 6.25*seednc
               lncmin = sla0*slnmin
               fcar   = 1.-fpro-ffat-flig-foac-fmin
               cfo    = 0.444*fcar+0.531*fpro+0.774*ffat+0.667*flig+0.368*foac
               clvi   = npl * seedw * cfo * eg * fcrsh
               crti   = npl * seedw * cfo * eg * (1.-fcrsh)
               ygo    = cfo/(1.275*fcar+1.887*fpro+3.189*ffat+2.231*flig+0.954* &
                        foac)*30./12.
               nlvi   = min(0.75 * npl * seedw * eg * seednc, lnci * clvi/cfv)
               laii   = clvi/cfv*sla0
               nrti   = npl * seedw * eg * seednc - nlvi
               slnbi  = nlvi/laii

               call gecros_init(xlat(i,j),hti,rdi,clvi,crti,nlvi,laii,nrti,slnbi,gecros_state(i,:,j))

             end if

! Noah-MP irrigation scheme !pvk
             if(iopt_irr >= 1 .and. iopt_irr <= 3) then
                if(iopt_irrm == 0 .or. iopt_irrm ==1) then       ! sprinkler
                   irnumsi(i,j) = 0
                   irwatsi(i,j) = 0.
                   ireloss(i,j) = 0.
                   irrsplh(i,j) = 0.
                else if (iopt_irrm == 0 .or. iopt_irrm ==2) then ! micro or drip
                   irnummi(i,j) = 0
                   irwatmi(i,j) = 0.
                   irmivol(i,j) = 0.
                else if (iopt_irrm == 0 .or. iopt_irrm ==3) then ! flood
                   irnumfi(i,j) = 0
                   irwatfi(i,j) = 0.
                   irfivol(i,j) = 0.
                end if
             end if

	   END IF

          enddo
       enddo


       ! Given the soil layer thicknesses (in DZS), initialize the soil layer
       ! depths from the surface.
       ZSOIL(1)         = -DZS(1)          ! negative
       DO NS=2, NSOIL
          ZSOIL(NS)       = ZSOIL(NS-1) - DZS(NS)
       END DO

       ! Initialize snow/soil layer arrays ZSNSOXY, TSNOXY, SNICEXY, SNLIQXY,
       ! and ISNOWXY
       CALL snow_init ( ims , ime , jms , jme , its , itf , jts , jtf , 3 , &
            &           NSOIL , zsoil , snow , tgxy , snowh ,     &
            &           zsnsoxy , tsnoxy , snicexy , snliqxy , isnowxy )

       !initialize arrays for groundwater dynamics iopt_run=5

       if(iopt_run.eq.5) then
          IF ( PRESENT(smoiseq)     .AND. &
            PRESENT(smcwtdxy)    .AND. &
            PRESENT(rechxy)      .AND. &
            PRESENT(deeprechxy)  .AND. &
            PRESENT(areaxy)      .AND. &
            PRESENT(dx)          .AND. &
            PRESENT(dy)          .AND. &
            PRESENT(msftx)       .AND. &
            PRESENT(msfty)       .AND. &
            PRESENT(wtddt)       .AND. &
            PRESENT(stepwtd)     .AND. &
            PRESENT(dt)          .AND. &
            PRESENT(qrfsxy)      .AND. &
            PRESENT(qspringsxy)  .AND. &
            PRESENT(qslatxy)     .AND. &
            PRESENT(fdepthxy)    .AND. &
            PRESENT(ht)          .AND. &
            PRESENT(riverbedxy)  .AND. &
            PRESENT(eqzwt)       .AND. &
            PRESENT(rivercondxy) .AND. &
            PRESENT(pexpxy)      .AND. &
            PRESENT(rechclim)    ) THEN

             STEPWTD = nint(WTDDT*60./DT)
             STEPWTD = max(STEPWTD,1)

          ELSE
             CALL wrf_error_fatal ('Not enough fields to use groundwater option in Noah-MP')
          END IF
       endif

    ENDIF

  END SUBROUTINE NOAHMP_INIT

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

  SUBROUTINE SNOW_INIT ( ims , ime , jms , jme , its , itf , jts , jtf ,                  &
       &                 NSNOW , NSOIL , ZSOIL , SWE , TGXY , SNODEP ,                    &
       &                 ZSNSOXY , TSNOXY , SNICEXY ,SNLIQXY , ISNOWXY )
!------------------------------------------------------------------------------------------
!   Initialize snow arrays for Noah-MP LSM, based in input SNOWDEP, NSNOW
!   ISNOWXY is an index array, indicating the index of the top snow layer.  Valid indices
!           for snow layers range from 0 (no snow) and -1 (shallow snow) to (-NSNOW)+1 (deep snow).
!   TSNOXY holds the temperature of the snow layer.  Snow layers are initialized with
!          temperature = ground temperature [?].  Snow-free levels in the array have value 0.0
!   SNICEXY is the frozen content of a snow layer.  Initial estimate based on SNODEP and SWE
!   SNLIQXY is the liquid content of a snow layer.  Initialized to 0.0
!   ZNSNOXY is the layer depth from the surface.
!------------------------------------------------------------------------------------------
    IMPLICIT NONE
!------------------------------------------------------------------------------------------
    INTEGER, INTENT(IN)                              :: ims, ime, jms, jme
    INTEGER, INTENT(IN)                              :: its, itf, jts, jtf
    INTEGER, INTENT(IN)                              :: NSNOW
    INTEGER, INTENT(IN)                              :: NSOIL
    REAL,    INTENT(IN), DIMENSION(ims:ime, jms:jme) :: SWE
    REAL,    INTENT(IN), DIMENSION(ims:ime, jms:jme) :: SNODEP
    REAL,    INTENT(IN), DIMENSION(ims:ime, jms:jme) :: TGXY
    REAL,    INTENT(IN), DIMENSION(1:NSOIL)          :: ZSOIL

    INTEGER, INTENT(OUT), DIMENSION(ims:ime, jms:jme)                :: ISNOWXY ! Top snow layer index
    REAL,    INTENT(OUT), DIMENSION(ims:ime, -NSNOW+1:NSOIL,jms:jme) :: ZSNSOXY ! Snow/soil layer depth from surface [m]
    REAL,    INTENT(OUT), DIMENSION(ims:ime, -NSNOW+1:    0,jms:jme) :: TSNOXY  ! Snow layer temperature [K]
    REAL,    INTENT(OUT), DIMENSION(ims:ime, -NSNOW+1:    0,jms:jme) :: SNICEXY ! Snow layer ice content [mm]
    REAL,    INTENT(OUT), DIMENSION(ims:ime, -NSNOW+1:    0,jms:jme) :: SNLIQXY ! snow layer liquid content [mm]

! Local variables:
!   DZSNO   holds the thicknesses of the various snow layers.
!   DZSNOSO holds the thicknesses of the various soil/snow layers.
    INTEGER                           :: I,J,IZ
    REAL,   DIMENSION(-NSNOW+1:    0) :: DZSNO
    REAL,   DIMENSION(-NSNOW+1:NSOIL) :: DZSNSO

!------------------------------------------------------------------------------------------

    DO J = jts , jtf
       DO I = its , itf
          IF ( SNODEP(I,J) < 0.025 ) THEN
             ISNOWXY(I,J) = 0
             DZSNO(-NSNOW+1:0) = 0.
          ELSE
             IF ( ( SNODEP(I,J) >= 0.025 ) .AND. ( SNODEP(I,J) <= 0.05 ) ) THEN
                ISNOWXY(I,J)    = -1
                DZSNO(0)  = SNODEP(I,J)
             ELSE IF ( ( SNODEP(I,J) > 0.05 ) .AND. ( SNODEP(I,J) <= 0.10 ) ) THEN
                ISNOWXY(I,J)    = -2
                DZSNO(-1) = SNODEP(I,J)/2.
                DZSNO( 0) = SNODEP(I,J)/2.
             ELSE IF ( (SNODEP(I,J) > 0.10 ) .AND. ( SNODEP(I,J) <= 0.25 ) ) THEN
                ISNOWXY(I,J)    = -2
                DZSNO(-1) = 0.05
                DZSNO( 0) = SNODEP(I,J) - DZSNO(-1)
             ELSE IF ( ( SNODEP(I,J) > 0.25 ) .AND. ( SNODEP(I,J) <= 0.45 ) ) THEN
                ISNOWXY(I,J)    = -3
                DZSNO(-2) = 0.05
                DZSNO(-1) = 0.5*(SNODEP(I,J)-DZSNO(-2))
                DZSNO( 0) = 0.5*(SNODEP(I,J)-DZSNO(-2))
             ELSE IF ( SNODEP(I,J) > 0.45 ) THEN
                ISNOWXY(I,J)     = -3
                DZSNO(-2) = 0.05
                DZSNO(-1) = 0.20
                DZSNO( 0) = SNODEP(I,J) - DZSNO(-1) - DZSNO(-2)
             ELSE
                CALL wrf_error_fatal("Problem with the logic assigning snow layers.")
             END IF
          END IF

          TSNOXY (I,-NSNOW+1:0,J) = 0.
          SNICEXY(I,-NSNOW+1:0,J) = 0.
          SNLIQXY(I,-NSNOW+1:0,J) = 0.
          DO IZ = ISNOWXY(I,J)+1 , 0
             TSNOXY(I,IZ,J)  = TGXY(I,J)  ! [k]
             SNLIQXY(I,IZ,J) = 0.00
             SNICEXY(I,IZ,J) = 1.00 * DZSNO(IZ) * (SWE(I,J)/SNODEP(I,J))  ! [kg/m3]
          END DO

          ! Assign local variable DZSNSO, the soil/snow layer thicknesses, for snow layers
          DO IZ = ISNOWXY(I,J)+1 , 0
             DZSNSO(IZ) = -DZSNO(IZ)
          END DO

          ! Assign local variable DZSNSO, the soil/snow layer thicknesses, for soil layers
          DZSNSO(1) = ZSOIL(1)
          DO IZ = 2 , NSOIL
             DZSNSO(IZ) = (ZSOIL(IZ) - ZSOIL(IZ-1))
          END DO

          ! Assign ZSNSOXY, the layer depths, for soil and snow layers
          ZSNSOXY(I,ISNOWXY(I,J)+1,J) = DZSNSO(ISNOWXY(I,J)+1)
          DO IZ = ISNOWXY(I,J)+2 , NSOIL
             ZSNSOXY(I,IZ,J) = ZSNSOXY(I,IZ-1,J) + DZSNSO(IZ)
          ENDDO

       END DO
    END DO

  END SUBROUTINE SNOW_INIT
! ==================================================================================================
! ----------------------------------------------------------------------
    SUBROUTINE GROUNDWATER_INIT (   &
            &            GRID, NSOIL , DZS, ISLTYP, IVGTYP, WTDDT , &
            &            FDEPTH, TOPO, RIVERBED, EQWTD, RIVERCOND, PEXP , AREA ,WTD ,  &
            &            SMOIS,SH2O, SMOISEQ, SMCWTDXY,  &
            &            QLATXY, QSLATXY, QRFXY, QRFSXY, &
            &            DEEPRECHXY, RECHXY , QSPRINGXY, QSPRINGSXY, &
            &            rechclim  ,                                   &
            &            ids,ide, jds,jde, kds,kde,                    &
            &            ims,ime, jms,jme, kms,kme,                    &
            &            ips,ipe, jps,jpe, kps,kpe,                    &
            &            its,ite, jts,jte, kts,kte                     )


  USE NOAHMP_TABLES, ONLY : BEXP_TABLE,SMCMAX_TABLE,PSISAT_TABLE,SMCWLT_TABLE,DWSAT_TABLE,DKSAT_TABLE, &
                                ISURBAN_TABLE, ISICE_TABLE ,ISWATER_TABLE
  USE module_sf_noahmp_groundwater, ONLY : LATERALFLOW
  USE module_domain, only: domain
#if (EM_CORE == 1)
#ifdef DM_PARALLEL
    USE module_dm        , ONLY : ntasks_x,ntasks_y,local_communicator,mytask,ntasks
    USE module_comm_dm , ONLY : halo_em_hydro_noahmp_sub
#endif
#endif

! ----------------------------------------------------------------------
  IMPLICIT NONE
! ----------------------------------------------------------------------

    INTEGER,  INTENT(IN   )   ::     ids,ide, jds,jde, kds,kde,  &
         &                           ims,ime, jms,jme, kms,kme,  &
         &                           ips,ipe, jps,jpe, kps,kpe,  &
         &                           its,ite, jts,jte, kts,kte
    TYPE(domain) , TARGET :: grid                             ! state
    INTEGER, INTENT(IN)                              :: NSOIL
    REAL,   INTENT(IN)                               ::     WTDDT
    REAL,    INTENT(IN), DIMENSION(1:NSOIL)          :: DZS
    INTEGER, INTENT(IN), DIMENSION(ims:ime, jms:jme) :: ISLTYP, IVGTYP
    REAL,    INTENT(IN), DIMENSION(ims:ime, jms:jme) :: FDEPTH, TOPO , AREA
    REAL,    INTENT(IN), DIMENSION(ims:ime, jms:jme) :: rechclim
    REAL,    INTENT(OUT), DIMENSION(ims:ime, jms:jme) :: RIVERCOND
    REAL,    INTENT(INOUT), DIMENSION(ims:ime, jms:jme) :: WTD, RIVERBED, EQWTD, PEXP
    REAL,     DIMENSION( ims:ime , 1:nsoil, jms:jme ), &
         &    INTENT(INOUT)   ::                          SMOIS, &
         &                                                 SH2O, &
         &                                                 SMOISEQ
    REAL,    INTENT(INOUT), DIMENSION(ims:ime, jms:jme) ::  &
                                                           SMCWTDXY, &
                                                           DEEPRECHXY, &
                                                           RECHXY, &
                                                           QSLATXY, &
                                                           QRFSXY, &
                                                           QSPRINGSXY, &
                                                           QLATXY, &
                                                           QRFXY, &
                                                           QSPRINGXY

! local
    INTEGER  :: I,J,K,ITER,itf,jtf, NITER, NCOUNT,NS
    REAL :: BEXP,SMCMAX,PSISAT,SMCWLT,DWSAT,DKSAT
    REAL :: FRLIQ,SMCEQDEEP
    REAL :: DELTAT,RCOND,TOTWATER
    REAL :: AA,BBB,CC,DD,DX,FUNC,DFUNC,DDZ,EXPON,SMC,FLUX
    REAL, DIMENSION(1:NSOIL) :: SMCEQ,ZSOIL
    REAL,      DIMENSION( ims:ime, jms:jme )    :: QLAT, QRF
    INTEGER,   DIMENSION( ims:ime, jms:jme )    :: LANDMASK !-1 for water (ice or no ice) and glacial areas, 1 for land where the LSM does its soil moisture calculations

       ! Given the soil layer thicknesses (in DZS), calculate the soil layer
       ! depths from the surface.
       ZSOIL(1)         = -DZS(1)          ! negative
       DO NS=2, NSOIL
          ZSOIL(NS)       = ZSOIL(NS-1) - DZS(NS)
       END DO


       itf=min0(ite,ide-1)
       jtf=min0(jte,jde-1)


    WHERE(IVGTYP.NE.ISWATER_TABLE.AND.IVGTYP.NE.ISICE_TABLE)
         LANDMASK=1
    ELSEWHERE
         LANDMASK=-1
    ENDWHERE

    PEXP = 1.0

    DELTAT=365.*24*3600. !1 year

!readjust the raw aggregated water table from hires, so that it is better compatible with topography

!use WTD here, to use the lateral communication routine
    WTD=EQWTD

    NCOUNT=0

 DO NITER=1,500

#if (EM_CORE == 1)
#ifdef DM_PARALLEL
#     include "HALO_EM_HYDRO_NOAHMP.inc"
#endif
#endif

!Calculate lateral flow

IF(NCOUNT.GT.0.OR.NITER.eq.1)THEN
    QLAT = 0.
    CALL LATERALFLOW(ISLTYP,WTD,QLAT,FDEPTH,TOPO,LANDMASK,DELTAT,AREA       &
                        ,ids,ide,jds,jde,kds,kde                      &
                        ,ims,ime,jms,jme,kms,kme                      &
                        ,its,ite,jts,jte,kts,kte                      )

    NCOUNT=0
    DO J=jts,jtf
       DO I=its,itf
          IF(LANDMASK(I,J).GT.0)THEN
            IF(QLAT(i,j).GT.1.e-2)THEN
                 NCOUNT=NCOUNT+1
                 WTD(i,j)=min(WTD(i,j)+0.25,0.)
            ENDIF
          ENDIF
        ENDDO
     ENDDO
ENDIF

 ENDDO

#if (EM_CORE == 1)
#ifdef DM_PARALLEL
#     include "HALO_EM_HYDRO_NOAHMP.inc"
#endif
#endif

EQWTD=WTD

!after adjusting, where qlat > 1cm/year now wtd is at the surface.
!it may still happen that qlat + rech > 0 and eqwtd-rbed <0. There the wtd can
!rise to the surface (poor drainage) but the et will then increase.


!now, calculate rcond:

    DO J=jts,jtf
       DO I=its,itf

        DDZ = EQWTD(I,J)- ( RIVERBED(I,J)-TOPO(I,J) )
!dont allow riverbed above water table
        IF(DDZ.LT.0.)then
               RIVERBED(I,J)=TOPO(I,J)+EQWTD(I,J)
               DDZ=0.
        ENDIF


        TOTWATER = AREA(I,J)*(QLAT(I,J)+RECHCLIM(I,J)*0.001)/DELTAT

        IF (TOTWATER.GT.0) THEN
              RIVERCOND(I,J) = TOTWATER / MAX(DDZ,0.05)
        ELSE
              RIVERCOND(I,J)=0.01
!and make riverbed  equal to eqwtd, otherwise qrf might be too big...
              RIVERBED(I,J)=TOPO(I,J)+EQWTD(I,J)
        ENDIF


       ENDDO
    ENDDO

!make riverbed to be height down from the surface instead of above sea level

    RIVERBED = min( RIVERBED-TOPO, 0.)

!now recompute lateral flow and flow to rivers to initialize deep soil moisture

    DELTAT = WTDDT * 60. !timestep in seconds for this calculation


!recalculate lateral flow

    QLAT = 0.
    CALL LATERALFLOW(ISLTYP,WTD,QLAT,FDEPTH,TOPO,LANDMASK,DELTAT,AREA       &
                        ,ids,ide,jds,jde,kds,kde                      &
                        ,ims,ime,jms,jme,kms,kme                      &
                        ,its,ite,jts,jte,kts,kte                      )

!compute flux from grounwater to rivers in the cell

    DO J=jts,jtf
       DO I=its,itf
          IF(LANDMASK(I,J).GT.0)THEN
             IF(WTD(I,J) .GT. RIVERBED(I,J) .AND.  EQWTD(I,J) .GT. RIVERBED(I,J)) THEN
               RCOND = RIVERCOND(I,J) * EXP(PEXP(I,J)*(WTD(I,J)-EQWTD(I,J)))
             ELSE
               RCOND = RIVERCOND(I,J)
             ENDIF
             QRF(I,J) = RCOND * (WTD(I,J)-RIVERBED(I,J)) * DELTAT/AREA(I,J)
!for now, dont allow it to go from river to groundwater
             QRF(I,J) = MAX(QRF(I,J),0.)
          ELSE
             QRF(I,J) = 0.
          ENDIF
       ENDDO
    ENDDO

!now compute eq. soil moisture, change soil moisture to be compatible with the water table and compute deep soil moisture

       DO J = jts,jtf
          DO I = its,itf
             BEXP   =   BEXP_TABLE(ISLTYP(I,J))
             SMCMAX = SMCMAX_TABLE(ISLTYP(I,J))
             SMCWLT = SMCWLT_TABLE(ISLTYP(I,J))
             IF(IVGTYP(I,J)==ISURBAN_TABLE)THEN
                 SMCMAX = 0.45
                 SMCWLT = 0.40
             ENDIF
             DWSAT  =   DWSAT_TABLE(ISLTYP(I,J))
             DKSAT  =   DKSAT_TABLE(ISLTYP(I,J))
             PSISAT = -PSISAT_TABLE(ISLTYP(I,J))
           IF ( ( BEXP > 0.0 ) .AND. ( smcmax > 0.0 ) .AND. ( -psisat > 0.0 ) ) THEN
             !initialize equilibrium soil moisture for water table diagnostic
                    CALL EQSMOISTURE(NSOIL ,  ZSOIL , SMCMAX , SMCWLT ,DWSAT, DKSAT  ,BEXP  , & !in
                                     SMCEQ                          )  !out

             SMOISEQ (I,1:NSOIL,J) = SMCEQ (1:NSOIL)


              !make sure that below the water table the layers are saturated and initialize the deep soil moisture
             IF(WTD(I,J) < ZSOIL(NSOIL)-DZS(NSOIL)) THEN

!initialize deep soil moisture so that the flux compensates qlat+qrf
!use Newton-Raphson method to find soil moisture

                         EXPON = 2. * BEXP + 3.
                         DDZ = ZSOIL(NSOIL) - WTD(I,J)
                         CC = PSISAT/DDZ
                         FLUX = (QLAT(I,J)-QRF(I,J))/DELTAT

                         SMC = 0.5 * SMCMAX

                         DO ITER = 1, 100
                           DD = (SMC+SMCMAX)/(2.*SMCMAX)
                           AA = -DKSAT * DD  ** EXPON
                           BBB = CC * ( (SMCMAX/SMC)**BEXP - 1. ) + 1.
                           FUNC =  AA * BBB - FLUX
                           DFUNC = -DKSAT * (EXPON/(2.*SMCMAX)) * DD ** (EXPON - 1.) * BBB &
                                   + AA * CC * (-BEXP) * SMCMAX ** BEXP * SMC ** (-BEXP-1.)

                           DX = FUNC/DFUNC
                           SMC = SMC - DX
                           IF ( ABS (DX) < 1.E-6)EXIT
                         ENDDO

                  SMCWTDXY(I,J) = MAX(SMC,1.E-4)

             ELSEIF(WTD(I,J) < ZSOIL(NSOIL))THEN
                  SMCEQDEEP = SMCMAX * ( PSISAT / ( PSISAT - DZS(NSOIL) ) ) ** (1./BEXP)
!                  SMCEQDEEP = MAX(SMCEQDEEP,SMCWLT)
                  SMCEQDEEP = MAX(SMCEQDEEP,1.E-4)
                  SMCWTDXY(I,J) = SMCMAX * ( WTD(I,J) -  (ZSOIL(NSOIL)-DZS(NSOIL))) + &
                                  SMCEQDEEP * (ZSOIL(NSOIL) - WTD(I,J))

             ELSE !water table within the resolved layers
                  SMCWTDXY(I,J) = SMCMAX
                  DO K=NSOIL,2,-1
                     IF(WTD(I,J) .GE. ZSOIL(K-1))THEN
                          FRLIQ = SH2O(I,K,J) / SMOIS(I,K,J)
                          SMOIS(I,K,J) = SMCMAX
                          SH2O(I,K,J) = SMCMAX * FRLIQ
                     ELSE
                          IF(SMOIS(I,K,J).LT.SMCEQ(K))THEN
                              WTD(I,J) = ZSOIL(K)
                          ELSE
                              WTD(I,J) = ( SMOIS(I,K,J)*DZS(K) - SMCEQ(K)*ZSOIL(K-1) + SMCMAX*ZSOIL(K) ) / &
                                         (SMCMAX - SMCEQ(K))
                          ENDIF
                          EXIT
                     ENDIF
                  ENDDO
             ENDIF
            ELSE
              SMOISEQ (I,1:NSOIL,J) = SMCMAX
              SMCWTDXY(I,J) = SMCMAX
              WTD(I,J) = 0.0
            ENDIF

!zero out some arrays

             QLATXY(I,J) = 0.0
             QSLATXY(I,J) = 0.0
             QRFXY(I,J) = 0.0
             QRFSXY(I,J) = 0.0
             DEEPRECHXY(I,J) = 0.0
             RECHXY(I,J) = 0.0
             QSPRINGXY(I,J) = 0.0
             QSPRINGSXY(I,J) = 0.0

          ENDDO
       ENDDO




    END  SUBROUTINE GROUNDWATER_INIT
! ==================================================================================================
! ----------------------------------------------------------------------
  SUBROUTINE EQSMOISTURE(NSOIL  ,  ZSOIL , SMCMAX , SMCWLT, DWSAT , DKSAT ,BEXP , & !in
                         SMCEQ                          )  !out
! ----------------------------------------------------------------------
  IMPLICIT NONE
! ----------------------------------------------------------------------
! input
  INTEGER,                         INTENT(IN) :: NSOIL !no. of soil layers
  REAL, DIMENSION(       1:NSOIL), INTENT(IN) :: ZSOIL !depth of soil layer-bottom [m]
  REAL,                            INTENT(IN) :: SMCMAX , SMCWLT, BEXP , DWSAT, DKSAT
!output
  REAL,  DIMENSION(      1:NSOIL), INTENT(OUT) :: SMCEQ  !equilibrium soil water  content [m3/m3]
!local
  INTEGER                                     :: K , ITER
  REAL                                        :: DDZ , SMC, FUNC, DFUNC , AA, BB , EXPON, DX

!gmmcompute equilibrium soil moisture content for the layer when wtd=zsoil(k)


   DO K=1,NSOIL

            IF ( K == 1 )THEN
                DDZ = -ZSOIL(K+1) * 0.5
            ELSEIF ( K < NSOIL ) THEN
                DDZ = ( ZSOIL(K-1) - ZSOIL(K+1) ) * 0.5
            ELSE
                DDZ = ZSOIL(K-1) - ZSOIL(K)
            ENDIF

!use Newton-Raphson method to find eq soil moisture

            EXPON = BEXP +1.
            AA = DWSAT/DDZ
            BB = DKSAT / SMCMAX ** EXPON

            SMC = 0.5 * SMCMAX

         DO ITER = 1, 100
            FUNC = (SMC - SMCMAX) * AA +  BB * SMC ** EXPON
            DFUNC = AA + BB * EXPON * SMC ** BEXP

            DX = FUNC/DFUNC
            SMC = SMC - DX
            IF ( ABS (DX) < 1.E-6)EXIT
         ENDDO

!             SMCEQ(K) = MIN(MAX(SMC,SMCWLT),SMCMAX*0.99)
             SMCEQ(K) = MIN(MAX(SMC,1.E-4),SMCMAX*0.99)
   ENDDO

END  SUBROUTINE EQSMOISTURE

! gecros initialization routines

SUBROUTINE gecros_init(xlat,hti,rdi,clvi,crti,nlvi,laii,nrti,slnbi,state_gecros)
implicit none
REAL, INTENT(IN)     :: HTI
REAL, INTENT(IN)     :: RDI
REAL, INTENT(IN)     :: CLVI
REAL, INTENT(IN)     :: CRTI
REAL, INTENT(IN)     :: NLVI
REAL, INTENT(IN)     :: LAII
REAL, INTENT(IN)     :: NRTI
REAL, INTENT(IN)     :: SLNBI
REAL, INTENT(IN)     :: XLAT
REAL, DIMENSION(1:60), INTENT(INOUT) :: STATE_GECROS

  !Inititalization of Gecros variables
  STATE_GECROS(1) = 0.      !DS
  STATE_GECROS(2) = 0.      !CTDURDI, HTI, CLVI, CRTI, NLVI, LAII, NRTI, SLNBI,
  STATE_GECROS(3) = 0.      !CVDU
  STATE_GECROS(4) = CLVI    !CLV
  STATE_GECROS(5) = 0.      !CLVD
  STATE_GECROS(6) = 0.      !CSST
  STATE_GECROS(7) = 0.      !CSO
  STATE_GECROS(8) = CRTI    !CSRT
  STATE_GECROS(9) =  0.     !CRTD
  STATE_GECROS(10) = 0.     !CLVDS
  STATE_GECROS(11) = NRTI   !NRT
  STATE_GECROS(12) = 0.     !NST
  STATE_GECROS(13) = NLVI   !NLV
  STATE_GECROS(14) = 0.     !NSO
  STATE_GECROS(15) = NLVI   !TNLV
  STATE_GECROS(16) = 0.     !NLVD
  STATE_GECROS(17) = 0.     !NRTD
  STATE_GECROS(18) = 0.     !CRVS
  STATE_GECROS(19) = 0.     !CRVR
  STATE_GECROS(20) = 0.     !NREOE
  STATE_GECROS(21) = 0.     !NREOF
  STATE_GECROS(22) = 0.     !DCDSR
  STATE_GECROS(23) = 0.     !DCDTR
  STATE_GECROS(24) = SLNBI  !SLNB
  STATE_GECROS(25) = LAII   !LAIC
  STATE_GECROS(26) = 0.     !RMUL
  STATE_GECROS(27) = 0.     !NDEMP
  STATE_GECROS(28) = 0.     !NSUPP
  STATE_GECROS(29) = 0.     !NFIXT
  STATE_GECROS(30) = 0.     !NFIXR
  STATE_GECROS(31) = 0.     !DCDTP
  STATE_GECROS(32) = 0.01   !HTI
  STATE_GECROS(33) = RDI    !RDI
  STATE_GECROS(34) = 0.     !TPCAN
  STATE_GECROS(35) = 0.     !TRESP
  STATE_GECROS(36) = 0.     !TNUPT
  STATE_GECROS(37) = 0.     !LITNT
  STATE_GECROS(38) = 0.     !daysSinceDS1
  STATE_GECROS(39) = 0.     !daysSinceDS2
  STATE_GECROS(40) = -1.    !drilled -1:false, 1:true
  STATE_GECROS(41) = -1.    !emerged -1:false, 1:true
  STATE_GECROS(42) = -1.    !harvested -1:false, 1:true
  STATE_GECROS(43) = 0.     !TTEM
  STATE_GECROS(44) = XLAT   !GLAT
  STATE_GECROS(45) = 0.     !WSO
  STATE_GECROS(46) = 0.     !WSTRAW
  STATE_GECROS(47) = 0.     !GrainNC
  STATE_GECROS(48) = 0.     !StrawNC
  STATE_GECROS(49) = 0.01   !GLAI
  STATE_GECROS(50) = 0.01   !TLAI
  STATE_GECROS(51) = HTI    !Fields 51-58 set for reinitialization
  STATE_GECROS(52) = RDI
  STATE_GECROS(53) = CLVI
  STATE_GECROS(54) = CRTI
  STATE_GECROS(55) = NRTI
  STATE_GECROS(56) = NLVI
  STATE_GECROS(57) = SLNBI
  STATE_GECROS(58) = LAII

END SUBROUTINE gecros_init

SUBROUTINE gecros_reinit(STATE_GECROS)
implicit none
REAL, DIMENSION(1:60), INTENT(INOUT) :: STATE_GECROS

  !Re-inititalization of Gecros variables after harvest
  STATE_GECROS(1) = 0.               !DS
  STATE_GECROS(2) = 0.               !CTDU
  STATE_GECROS(3) = 0.               !CVDU
  STATE_GECROS(4) = STATE_GECROS(53) !CLV
  STATE_GECROS(5) = 0.               !CLVD
  STATE_GECROS(6) = 0.               !CSST
  STATE_GECROS(7) = 0.               !CSO
  STATE_GECROS(8) = STATE_GECROS(54) !CRT
  STATE_GECROS(9) = 0.               !CRTD
  STATE_GECROS(10) = 0.              !CLVDS
  STATE_GECROS(11) = STATE_GECROS(55)!NRT
  STATE_GECROS(12) = 0.              !NST
  STATE_GECROS(13) = STATE_GECROS(56)!NLV
  STATE_GECROS(14) = 0.              !NSO
  STATE_GECROS(15) = STATE_GECROS(56)!TNLV
  STATE_GECROS(16) = 0.              !NLVD
  STATE_GECROS(17) = 0.              !NRTD
  STATE_GECROS(18) = 0.              !CRVS
  STATE_GECROS(19) = 0.              !CRVR
  STATE_GECROS(20) = 0.              !NREOE
  STATE_GECROS(21) = 0.              !NREOF
  STATE_GECROS(22) = 0.              !DCDSR
  STATE_GECROS(23) = 0.              !DCDTR
  STATE_GECROS(24) = STATE_GECROS(57)!SLNB
  STATE_GECROS(25) = STATE_GECROS(58)!LAIC
  STATE_GECROS(26) = 0.              !RMUL
  STATE_GECROS(27) = 0.              !NDEMP
  STATE_GECROS(28) = 0.              !NSUPP
  STATE_GECROS(29) = 0.              !NFIXT
  STATE_GECROS(30) = 0.              !NFIXR
  STATE_GECROS(31) = 0.              !DCDTP
  STATE_GECROS(32) = STATE_GECROS(51)!HT
  STATE_GECROS(33) = STATE_GECROS(52)!ROOTD
  STATE_GECROS(34) = 0.              !TPCAN
  STATE_GECROS(35) = 0.              !TRESP
  STATE_GECROS(36) = 0.              !TNUPT
  STATE_GECROS(37) = 0.              !LITNT
  STATE_GECROS(38) = 0.              !daysSinceDS1
  STATE_GECROS(39) = 0.              !daysSinceDS2
  STATE_GECROS(40) = -1.             !drilled -1:false, 1:true
  STATE_GECROS(41) = -1.             !emerged -1:false, 1:true
  STATE_GECROS(42) = 1.              !harvested -1:false, 1:true
  STATE_GECROS(43) = 0.              !TTEM
  STATE_GECROS(45) = 0.              !WSO
  STATE_GECROS(46) = 0.              !WSTRAW
  STATE_GECROS(47) = 0.              !GrainNC
  STATE_GECROS(48) = 0.              !StrawNC
  STATE_GECROS(49) = 0.01            !GLAI
  STATE_GECROS(50) = 0.01            !TLAI

END SUBROUTINE gecros_reinit

!***Function for HARVEST DATES:

!Determine if crop is to be harvested today
!function to be called once a day
!return codes: 0 - no, 1- yes
!requires two counters 'daysSinceDS2', 'daysSinceDS1' , zero-initialized to be maintained within caller
!STATE_GECROS(1) = current DS
!STATE_GECROS(38)=daysSinceDS1
!STATE_GECROS(39)=daysSinceDS2

function checkIfHarvest(STATE_GECROS, DT, harvestDS1, harvestDS2, harvestDS1ExtraDays, harvestDS2ExtraDays)
implicit none
real :: DT, harvestDS1, harvestDS2
real :: daysSinceDS1, daysSinceDS2
real :: harvestDS1ExtraDays, harvestDS2ExtraDays
integer :: checkIfHarvest
REAL, DIMENSION(1:60), INTENT(INOUT) :: STATE_GECROS


 !***check whether maturity (DS1) has been reached
 if (STATE_GECROS(1) >= harvestDS1) then

    if (STATE_GECROS(38) >= harvestDS1ExtraDays) then
        checkIfHarvest=1
 !if we are > DS1, but not over the limit, increase the counter of days
    else
        STATE_GECROS(38) = STATE_GECROS(38) + DT/86400.
    endif
 else

 !if maturity has not been reached, but we are close (> DS2)
 !check the number of days for which we have been > DS2
 !and harvest in case we are over the limit given for that stage
 !(in case that maturity will not be reached at all)

 checkIfHarvest=0
 if (STATE_GECROS(1) >= harvestDS2 ) then

       if (STATE_GECROS(39) >= harvestDS2ExtraDays) then
           checkIfHarvest=1
       else !if we are > DS2, but not over the limit, increase the counter of days
           STATE_GECROS(39) = STATE_GECROS(39) + DT/86400.
           checkIfHarvest=0
      endif
 endif
 endif
 return
end function checkIfHarvest

!------------------------------------------------------------------------------------------

  SUBROUTINE noahmp_urban(sf_urban_physics,   NSOIL,         IVGTYP,  ITIMESTEP,            & ! IN : Model configuration
                                 DT,     COSZ_URB2D,     XLAT_URB2D,                        & ! IN : Time/Space-related
                                T3D,           QV3D,          U_PHY,      V_PHY,   SWDOWN,  & ! IN : Forcing
                             SWDDIR,         SWDDIF,                                        &
		                GLW,          P8W3D,         RAINBL,       DZ8W,      ZNT,  & ! IN : Forcing
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
                         declin_urb,      omg_urb2d,                                        & !I urban
                    num_roof_layers,num_wall_layers,num_road_layers,                        & !I urban
                                dzr,            dzb,            dzg,                        & !I urban
                         cmcr_urb2d,      tgr_urb2d,     tgrl_urb3d,  smr_urb3d,            & !H urban
                        drelr_urb2d,    drelb_urb2d,    drelg_urb2d,                        & !H urban
                      flxhumr_urb2d,  flxhumb_urb2d,  flxhumg_urb2d,                        & !H urban
                             julian,          julyr,                                        & !H urban
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
                             th_phy,            rho,          p_phy,        ust,            & !I multi-layer urban
                                gmt,         julday,          xlong,       xlat,            & !I multi-layer urban
                            a_u_bep,        a_v_bep,        a_t_bep,    a_q_bep,            & !O multi-layer urban
                            a_e_bep,        b_u_bep,        b_v_bep,                        & !O multi-layer urban
                            b_t_bep,        b_q_bep,        b_e_bep,    dlg_bep,            & !O multi-layer urban
                           dl_u_bep,         sf_bep,         vl_bep                         & !O multi-layer urban
                 )

  USE module_sf_urban,    only: urban
  USE module_sf_bep,      only: bep
  USE module_sf_bep_bem,  only: bep_bem
  USE module_ra_gfdleta,  only: cal_mon_day
  USE NOAHMP_TABLES, ONLY: ISURBAN_TABLE, LCZ_1_TABLE, LCZ_2_TABLE, LCZ_3_TABLE
  USE module_model_constants, only: KARMAN, CP, XLV
!----------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------

    INTEGER,                                         INTENT(IN   ) ::  sf_urban_physics   ! urban physics option
    INTEGER,                                         INTENT(IN   ) ::  NSOIL     ! number of soil layers
    INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  IVGTYP    ! vegetation type
    INTEGER,                                         INTENT(IN   ) ::  ITIMESTEP ! timestep number
    REAL,                                            INTENT(IN   ) ::  DT        ! timestep [s]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  COSZ_URB2D
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  XLAT_URB2D
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  T3D       ! 3D atmospheric temperature valid at mid-levels [K]
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  QV3D      ! 3D water vapor mixing ratio [kg/kg_dry]
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  U_PHY     ! 3D U wind component [m/s]
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  V_PHY     ! 3D V wind component [m/s]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SWDOWN    ! solar down at surface [W m-2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SWDDIF    ! solar down at surface [W m-2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SWDDIR    ! solar down at surface [W m-2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  GLW       ! longwave down at surface [W m-2]
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  P8W3D     ! 3D pressure, valid at interface [Pa]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  RAINBL    ! total input precipitation [mm]
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  DZ8W      ! thickness of atmo layers [m]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ZNT       ! combined z0 sent to coupled model
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  TSK       ! surface radiative temperature [K]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  HFX       ! sensible heat flux [W m-2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  QFX       ! latent heat flux [kg s-1 m-2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  LH        ! latent heat flux [W m-2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  GRDFLX    ! ground/snow heat flux [W m-2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ALBEDO    ! total grid albedo []
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  EMISS     ! surface bulk emissivity
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  QSFC      ! bulk surface mixing ratio

    INTEGER,  INTENT(IN   )   ::     ids,ide, jds,jde, kds,kde,  &  ! d -> domain
         &                           ims,ime, jms,jme, kms,kme,  &  ! m -> memory
         &                           its,ite, jts,jte, kts,kte      ! t -> tile

! input variables surface_driver --> lsm

     INTEGER,                                                INTENT(IN   ) :: num_roof_layers
     INTEGER,                                                INTENT(IN   ) :: num_wall_layers
     INTEGER,                                                INTENT(IN   ) :: num_road_layers

     INTEGER,        DIMENSION( ims:ime, jms:jme ),          INTENT(IN   ) :: UTYPE_URB2D
     REAL,           DIMENSION( ims:ime, jms:jme ),          INTENT(IN   ) :: FRC_URB2D

     REAL, OPTIONAL, DIMENSION(1:num_roof_layers),           INTENT(IN   ) :: DZR
     REAL, OPTIONAL, DIMENSION(1:num_wall_layers),           INTENT(IN   ) :: DZB
     REAL, OPTIONAL, DIMENSION(1:num_road_layers),           INTENT(IN   ) :: DZG
     REAL, OPTIONAL,                                         INTENT(IN   ) :: DECLIN_URB
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),          INTENT(IN   ) :: OMG_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: TH_PHY
     REAL, OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: P_PHY
     REAL, OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: RHO

     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),          INTENT(INOUT) :: UST
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),          INTENT(INOUT) :: CHS, CHS2, CQS2

     INTEGER,  INTENT(IN   )   ::  julian, julyr                  !urban

! local variables lsm --> urban

     INTEGER :: UTYPE_URB ! urban type [urban=1, suburban=2, rural=3]
     REAL    :: TA_URB       ! potential temp at 1st atmospheric level [K]
     REAL    :: QA_URB       ! mixing ratio at 1st atmospheric level  [kg/kg]
     REAL    :: UA_URB       ! wind speed at 1st atmospheric level    [m/s]
     REAL    :: U1_URB       ! u at 1st atmospheric level             [m/s]
     REAL    :: V1_URB       ! v at 1st atmospheric level             [m/s]
     REAL    :: SSG_URB      ! downward total short wave radiation    [W/m/m]
     REAL    :: LLG_URB      ! downward long wave radiation           [W/m/m]
     REAL    :: RAIN_URB     ! precipitation                          [mm/h]
     REAL    :: RHOO_URB     ! air density                            [kg/m^3]
     REAL    :: ZA_URB       ! first atmospheric level                [m]
     REAL    :: DELT_URB     ! time step                              [s]
     REAL    :: SSGD_URB     ! downward direct short wave radiation   [W/m/m]
     REAL    :: SSGQ_URB     ! downward diffuse short wave radiation  [W/m/m]
     REAL    :: XLAT_URB     ! latitude                               [deg]
     REAL    :: COSZ_URB     ! cosz
     REAL    :: OMG_URB      ! hour angle
     REAL    :: ZNT_URB      ! roughness length                       [m]
     REAL    :: TR_URB
     REAL    :: TB_URB
     REAL    :: TG_URB
     REAL    :: TC_URB
     REAL    :: QC_URB
     REAL    :: UC_URB
     REAL    :: XXXR_URB
     REAL    :: XXXB_URB
     REAL    :: XXXG_URB
     REAL    :: XXXC_URB
     REAL, DIMENSION(1:num_roof_layers) :: TRL_URB  ! roof layer temp [K]
     REAL, DIMENSION(1:num_wall_layers) :: TBL_URB  ! wall layer temp [K]
     REAL, DIMENSION(1:num_road_layers) :: TGL_URB  ! road layer temp [K]
     LOGICAL  :: LSOLAR_URB

!===hydrological variable for single layer UCM===

     INTEGER :: jmonth, jday
     REAL    :: DRELR_URB
     REAL    :: DRELB_URB
     REAL    :: DRELG_URB
     REAL    :: FLXHUMR_URB
     REAL    :: FLXHUMB_URB
     REAL    :: FLXHUMG_URB
     REAL    :: CMCR_URB
     REAL    :: TGR_URB

     REAL, DIMENSION(1:num_roof_layers) :: SMR_URB  ! green roof layer moisture
     REAL, DIMENSION(1:num_roof_layers) :: TGRL_URB ! green roof layer temp [K]

     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                    INTENT(INOUT) :: DRELR_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                    INTENT(INOUT) :: DRELB_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                    INTENT(INOUT) :: DRELG_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                    INTENT(INOUT) :: FLXHUMR_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                    INTENT(INOUT) :: FLXHUMB_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                    INTENT(INOUT) :: FLXHUMG_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                    INTENT(INOUT) :: CMCR_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                    INTENT(INOUT) :: TGR_URB2D

     REAL, OPTIONAL, DIMENSION( ims:ime, 1:num_roof_layers, jms:jme ), INTENT(INOUT) :: TGRL_URB3D
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:num_roof_layers, jms:jme ), INTENT(INOUT) :: SMR_URB3D


! state variable surface_driver <--> lsm <--> urban

     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: TR_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: TB_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: TG_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: TC_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: QC_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: UC_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: XXXR_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: XXXB_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: XXXG_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: XXXC_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: SH_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: LH_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: G_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: RN_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: TS_URB2D

     REAL, OPTIONAL, DIMENSION( ims:ime, 1:num_roof_layers, jms:jme ), INTENT(INOUT) :: TRL_URB3D
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:num_wall_layers, jms:jme ), INTENT(INOUT) :: TBL_URB3D
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:num_road_layers, jms:jme ), INTENT(INOUT) :: TGL_URB3D

! output variable lsm --> surface_driver

     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT) :: PSIM_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT) :: PSIH_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT) :: GZ1OZ0_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT) :: U10_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT) :: V10_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT) :: TH2_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT) :: Q2_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT) :: AKMS_URB2D
     REAL,           DIMENSION( ims:ime, jms:jme ), INTENT(OUT) :: UST_URB2D


! output variables urban --> lsm

     REAL :: TS_URB           ! surface radiative temperature    [K]
     REAL :: QS_URB           ! surface humidity                 [-]
     REAL :: SH_URB           ! sensible heat flux               [W/m/m]
     REAL :: LH_URB           ! latent heat flux                 [W/m/m]
     REAL :: LH_KINEMATIC_URB ! latent heat flux, kinetic  [kg/m/m/s]
     REAL :: SW_URB           ! upward short wave radiation flux [W/m/m]
     REAL :: ALB_URB          ! time-varying albedo            [fraction]
     REAL :: LW_URB           ! upward long wave radiation flux  [W/m/m]
     REAL :: G_URB            ! heat flux into the ground        [W/m/m]
     REAL :: RN_URB           ! net radiation                    [W/m/m]
     REAL :: PSIM_URB         ! shear f for momentum             [-]
     REAL :: PSIH_URB         ! shear f for heat                 [-]
     REAL :: GZ1OZ0_URB       ! shear f for heat                 [-]
     REAL :: U10_URB          ! wind u component at 10 m         [m/s]
     REAL :: V10_URB          ! wind v component at 10 m         [m/s]
     REAL :: TH2_URB          ! potential temperature at 2 m     [K]
     REAL :: Q2_URB           ! humidity at 2 m                  [-]
     REAL :: CHS_URB
     REAL :: CHS2_URB
     REAL :: UST_URB

! NUDAPT Parameters urban --> lam

     REAL :: mh_urb
     REAL :: stdh_urb
     REAL :: lp_urb
     REAL :: hgt_urb
     REAL, DIMENSION(4) :: lf_urb

! Local variables

     INTEGER :: I,J,K
     REAL :: Q1

! Noah UA changes

     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                     INTENT(INOUT) :: CMR_SFCDIF
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                     INTENT(INOUT) :: CHR_SFCDIF
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                     INTENT(INOUT) :: CMGR_SFCDIF
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                     INTENT(INOUT) :: CHGR_SFCDIF
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                     INTENT(INOUT) :: CMC_SFCDIF
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                     INTENT(INOUT) :: CHC_SFCDIF

! Variables for multi-layer UCM

     REAL, OPTIONAL,                                                    INTENT(IN   ) :: GMT
     INTEGER, OPTIONAL,                                                 INTENT(IN   ) :: JULDAY
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                     INTENT(IN   ) :: XLAT, XLONG
     INTEGER,                                                           INTENT(IN   ) :: num_urban_ndm
     INTEGER,                                                           INTENT(IN   ) :: urban_map_zrd
     INTEGER,                                                           INTENT(IN   ) :: urban_map_zwd
     INTEGER,                                                           INTENT(IN   ) :: urban_map_gd
     INTEGER,                                                           INTENT(IN   ) :: urban_map_zd
     INTEGER,                                                           INTENT(IN   ) :: urban_map_zdf
     INTEGER,                                                           INTENT(IN   ) :: urban_map_bd
     INTEGER,                                                           INTENT(IN   ) :: urban_map_wd
     INTEGER,                                                           INTENT(IN   ) :: urban_map_gbd
     INTEGER,                                                           INTENT(IN   ) :: urban_map_fbd
     INTEGER,                                                           INTENT(IN   ) :: urban_map_zgrd
     INTEGER,                                                           INTENT(IN   ) :: NUM_URBAN_HI
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:num_urban_hi, jms:jme ),     INTENT(IN   ) :: hi_urb2d
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                     INTENT(IN   ) :: lp_urb2d
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                     INTENT(IN   ) :: lb_urb2d
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                     INTENT(IN   ) :: hgt_urb2d
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                     INTENT(IN   ) :: mh_urb2d
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                     INTENT(IN   ) :: stdh_urb2d
     REAL, OPTIONAL, DIMENSION( ims:ime, 4, jms:jme ),                  INTENT(IN   ) :: lf_urb2d

     REAL, OPTIONAL, DIMENSION( ims:ime, 1:urban_map_zrd, jms:jme ),    INTENT(INOUT) :: trb_urb4d
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:urban_map_zwd, jms:jme ),    INTENT(INOUT) :: tw1_urb4d
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:urban_map_zwd, jms:jme ),    INTENT(INOUT) :: tw2_urb4d
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:urban_map_gd , jms:jme ),    INTENT(INOUT) :: tgb_urb4d
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:urban_map_bd , jms:jme ),    INTENT(INOUT) :: tlev_urb3d
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:urban_map_bd , jms:jme ),    INTENT(INOUT) :: qlev_urb3d
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:urban_map_wd , jms:jme ),    INTENT(INOUT) :: tw1lev_urb3d
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:urban_map_wd , jms:jme ),    INTENT(INOUT) :: tw2lev_urb3d
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:urban_map_gbd, jms:jme ),    INTENT(INOUT) :: tglev_urb3d
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:urban_map_fbd, jms:jme ),    INTENT(INOUT) :: tflev_urb3d
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                     INTENT(INOUT) :: lf_ac_urb3d
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                     INTENT(INOUT) :: sf_ac_urb3d
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                     INTENT(INOUT) :: cm_ac_urb3d
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                     INTENT(INOUT) :: sfvent_urb3d
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                     INTENT(INOUT) :: lfvent_urb3d
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:urban_map_wd , jms:jme ),    INTENT(INOUT) :: sfwin1_urb3d
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:urban_map_wd , jms:jme ),    INTENT(INOUT) :: sfwin2_urb3d
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:urban_map_zd , jms:jme ),    INTENT(INOUT) :: sfw1_urb3d
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:urban_map_zd , jms:jme ),    INTENT(INOUT) :: sfw2_urb3d
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:urban_map_zdf, jms:jme ),    INTENT(INOUT) :: sfr_urb3d
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:num_urban_ndm, jms:jme ),    INTENT(INOUT) :: sfg_urb3d
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: ep_pv_urb3d !GRZ
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:urban_map_zdf, jms:jme ), INTENT(INOUT) :: t_pv_urb3d !GRZ
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:urban_map_zgrd, jms:jme ),INTENT(INOUT) :: trv_urb4d !GRZ
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:urban_map_zgrd, jms:jme ),INTENT(INOUT) :: qr_urb4d !GRZ
     REAL, OPTIONAL, DIMENSION( ims:ime,jms:jme ), INTENT(INOUT) :: qgr_urb3d  !GRZ
     REAL, OPTIONAL, DIMENSION( ims:ime,jms:jme ), INTENT(INOUT) :: tgr_urb3d  !GRZ
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:urban_map_zdf, jms:jme ),INTENT(INOUT) :: drain_urb4d !GRZ
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: draingr_urb3d !GRZ
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:urban_map_zdf, jms:jme ),INTENT(INOUT) :: sfrv_urb3d !GRZ
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:urban_map_zdf, jms:jme ),INTENT(INOUT) :: lfrv_urb3d !GRZ
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:urban_map_zdf, jms:jme ),INTENT(INOUT) :: dgr_urb3d !GRZ
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:num_urban_ndm, jms:jme ),INTENT(INOUT) :: dg_urb3d !GRZ
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:urban_map_zdf, jms:jme ),INTENT(INOUT) :: lfr_urb3d !GRZ
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:num_urban_ndm, jms:jme ),INTENT(INOUT) :: lfg_urb3d !GRZ
     REAL, OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ),            INTENT(INOUT) :: a_u_bep   !Implicit momemtum component X-direction
     REAL, OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ),            INTENT(INOUT) :: a_v_bep   !Implicit momemtum component Y-direction
     REAL, OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ),            INTENT(INOUT) :: a_t_bep   !Implicit component pot. temperature
     REAL, OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ),            INTENT(INOUT) :: a_q_bep   !Implicit momemtum component X-direction
     REAL, OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ),            INTENT(INOUT) :: a_e_bep   !Implicit component TKE
     REAL, OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ),            INTENT(INOUT) :: b_u_bep   !Explicit momentum component X-direction
     REAL, OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ),            INTENT(INOUT) :: b_v_bep   !Explicit momentum component Y-direction
     REAL, OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ),            INTENT(INOUT) :: b_t_bep   !Explicit component pot. temperature
     REAL, OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ),            INTENT(INOUT) :: b_q_bep   !Implicit momemtum component Y-direction
     REAL, OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ),            INTENT(INOUT) :: b_e_bep   !Explicit component TKE
     REAL, OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ),            INTENT(INOUT) :: vl_bep    !Fraction air volume in grid cell
     REAL, OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ),            INTENT(INOUT) :: dlg_bep   !Height above ground
     REAL, OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ),            INTENT(INOUT) :: sf_bep    !Fraction air at the face of grid cell
     REAL, OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ),            INTENT(INOUT) :: dl_u_bep  !Length scale

! Local variables for multi-layer UCM

     REAL,    DIMENSION( its:ite, jts:jte) :: HFX_RURAL,GRDFLX_RURAL          ! ,LH_RURAL,RN_RURAL
     REAL,    DIMENSION( its:ite, jts:jte) :: QFX_RURAL                       ! ,QSFC_RURAL,UMOM_RURAL,VMOM_RURAL
     REAL,    DIMENSION( its:ite, jts:jte) :: ALB_RURAL,EMISS_RURAL,TSK_RURAL ! ,UST_RURAL
     REAL,    DIMENSION( its:ite, jts:jte) :: HFX_URB,UMOM_URB,VMOM_URB
     REAL,    DIMENSION( its:ite, jts:jte) :: QFX_URB
     REAL,    DIMENSION( its:ite, jts:jte) :: EMISS_URB
     REAL,    DIMENSION( its:ite, jts:jte) :: RL_UP_URB
     REAL,    DIMENSION( its:ite, jts:jte) :: RS_ABS_URB
     REAL,    DIMENSION( its:ite, jts:jte) :: GRDFLX_URB

     REAL :: SIGMA_SB,RL_UP_RURAL,RL_UP_TOT,RS_ABS_TOT,UMOM,VMOM
     REAL :: r1,r2,r3
     REAL :: CMR_URB, CHR_URB, CMC_URB, CHC_URB, CMGR_URB, CHGR_URB
     REAL :: frc_urb,lb_urb
     REAL :: check

    character(len=80) :: message

    DO J=JTS,JTE
    DO I=ITS,ITE
      HFX_RURAL(I,J)                = HFX(I,J)
      QFX_RURAL(I,J)                = QFX(I,J)
      GRDFLX_RURAL(I,J)             = GRDFLX(I,J)
      EMISS_RURAL(I,J)              = EMISS(I,J)
      TSK_RURAL(I,J)                = TSK(I,J)
      ALB_RURAL(I,J)                = ALBEDO(I,J)
    END DO
    END DO

IF (SF_URBAN_PHYSICS == 1 ) THEN         ! Beginning of UCM CALL if block

!--------------------------------------
! URBAN CANOPY MODEL START
!--------------------------------------

JLOOP : DO J = jts, jte

ILOOP : DO I = its, ite


  IF( IVGTYP(I,J) == ISURBAN_TABLE .or. IVGTYP(I,J) == LCZ_1_TABLE .or. &
      IVGTYP(I,J) == LCZ_2_TABLE .or. IVGTYP(I,J) == LCZ_3_TABLE ) THEN

    UTYPE_URB = UTYPE_URB2D(I,J) !urban type (low, high or industrial)

    TA_URB    = T3D(I,1,J)                                ! [K]
    QA_URB    = QV3D(I,1,J)/(1.0+QV3D(I,1,J))             ! [kg/kg]
    UA_URB    = SQRT(U_PHY(I,1,J)**2.+V_PHY(I,1,J)**2.)
    U1_URB    = U_PHY(I,1,J)
    V1_URB    = V_PHY(I,1,J)
    IF(UA_URB < 1.) UA_URB=1.                             ! [m/s]
    SSG_URB   = SWDOWN(I,J)                               ! [W/m/m]
    SSGD_URB  = 0.8*SWDOWN(I,J)                           ! [W/m/m]
    SSGQ_URB  = SSG_URB-SSGD_URB                          ! [W/m/m]
    LLG_URB   = GLW(I,J)                                  ! [W/m/m]
    RAIN_URB  = RAINBL(I,J) / DT * 3600.0                 ! [mm/hr]
    RHOO_URB  = (P8W3D(I,KTS+1,J)+P8W3D(I,KTS,J))*0.5 / (287.04 * TA_URB * (1.0+ 0.61 * QA_URB)) ![kg/m/m/m]
    ZA_URB    = 0.5*DZ8W(I,1,J)                           ! [m]
    DELT_URB  = DT                                        ! [sec]
    XLAT_URB  = XLAT_URB2D(I,J)                           ! [deg]
    COSZ_URB  = COSZ_URB2D(I,J)
    OMG_URB   = OMG_URB2D(I,J)
    ZNT_URB   = ZNT(I,J)

    LSOLAR_URB = .FALSE.

    TR_URB = TR_URB2D(I,J)
    TB_URB = TB_URB2D(I,J)
    TG_URB = TG_URB2D(I,J)
    TC_URB = TC_URB2D(I,J)
    QC_URB = QC_URB2D(I,J)
    UC_URB = UC_URB2D(I,J)

    TGR_URB     = TGR_URB2D(I,J)
    CMCR_URB    = CMCR_URB2D(I,J)
    FLXHUMR_URB = FLXHUMR_URB2D(I,J)
    FLXHUMB_URB = FLXHUMB_URB2D(I,J)
    FLXHUMG_URB = FLXHUMG_URB2D(I,J)
    DRELR_URB   = DRELR_URB2D(I,J)
    DRELB_URB   = DRELB_URB2D(I,J)
    DRELG_URB   = DRELG_URB2D(I,J)

    DO K = 1,num_roof_layers
      TRL_URB(K) = TRL_URB3D(I,K,J)
      SMR_URB(K) = SMR_URB3D(I,K,J)
      TGRL_URB(K)= TGRL_URB3D(I,K,J)
    END DO

    DO K = 1,num_wall_layers
      TBL_URB(K) = TBL_URB3D(I,K,J)
    END DO

    DO K = 1,num_road_layers
      TGL_URB(K) = TGL_URB3D(I,K,J)
    END DO

    XXXR_URB = XXXR_URB2D(I,J)
    XXXB_URB = XXXB_URB2D(I,J)
    XXXG_URB = XXXG_URB2D(I,J)
    XXXC_URB = XXXC_URB2D(I,J)

! Limits to avoid dividing by small number
    IF (CHS(I,J) < 1.0E-02) THEN
      CHS(I,J)  = 1.0E-02
    ENDIF
    IF (CHS2(I,J) < 1.0E-02) THEN
      CHS2(I,J)  = 1.0E-02
    ENDIF
    IF (CQS2(I,J) < 1.0E-02) THEN
      CQS2(I,J)  = 1.0E-02
    ENDIF

    CHS_URB  = CHS(I,J)
    CHS2(I,J)= CQS2(I,J)
    CHS2_URB = CHS2(I,J)
    IF (PRESENT(CMR_SFCDIF)) THEN
      CMR_URB = CMR_SFCDIF(I,J)
      CHR_URB = CHR_SFCDIF(I,J)
      CMGR_URB = CMGR_SFCDIF(I,J)
      CHGR_URB = CHGR_SFCDIF(I,J)
      CMC_URB = CMC_SFCDIF(I,J)
      CHC_URB = CHC_SFCDIF(I,J)
    ENDIF

! NUDAPT for SLUCM

    MH_URB   = MH_URB2D(I,J)
    STDH_URB = STDH_URB2D(I,J)
    LP_URB   = LP_URB2D(I,J)
    HGT_URB  = HGT_URB2D(I,J)
    LF_URB   = 0.0
    DO K = 1,4
      LF_URB(K) = LF_URB2D(I,K,J)
    ENDDO
    FRC_URB  = FRC_URB2D(I,J)
    LB_URB   = LB_URB2D(I,J)
    CHECK    = 0

! Call urban
    CALL cal_mon_day(julian,julyr,jmonth,jday)
    CALL urban(LSOLAR_URB,                                                             & ! I
          num_roof_layers, num_wall_layers, num_road_layers,                           & ! C
                DZR,        DZB,        DZG, & ! C
          UTYPE_URB,     TA_URB,     QA_URB,     UA_URB,   U1_URB,  V1_URB, SSG_URB,   & ! I
           SSGD_URB,   SSGQ_URB,    LLG_URB,   RAIN_URB, RHOO_URB,                     & ! I
             ZA_URB, DECLIN_URB,   COSZ_URB,    OMG_URB,                               & ! I
           XLAT_URB,   DELT_URB,    ZNT_URB,                                           & ! I
            CHS_URB,   CHS2_URB,                                                       & ! I
             TR_URB,     TB_URB,     TG_URB,     TC_URB,   QC_URB,   UC_URB,           & ! H
            TRL_URB,    TBL_URB,    TGL_URB,                                           & ! H
           XXXR_URB,   XXXB_URB,   XXXG_URB,   XXXC_URB,                               & ! H
             TS_URB,     QS_URB,     SH_URB,     LH_URB, LH_KINEMATIC_URB,             & ! O
             SW_URB,    ALB_URB,     LW_URB,      G_URB,   RN_URB, PSIM_URB, PSIH_URB, & ! O
         GZ1OZ0_URB,                                                                   & !O
            CMR_URB,    CHR_URB,    CMC_URB,    CHC_URB,                               &
            U10_URB,    V10_URB,    TH2_URB,     Q2_URB,                               & ! O
            UST_URB,     mh_urb,   stdh_urb,     lf_urb,   lp_urb,                     & ! 0
            hgt_urb,    frc_urb,     lb_urb,      check, CMCR_URB,TGR_URB,             & ! H
           TGRL_URB,    SMR_URB,   CMGR_URB,   CHGR_URB,   jmonth,                     & ! H
          DRELR_URB,  DRELB_URB,                                                       & ! H
          DRELG_URB,FLXHUMR_URB,FLXHUMB_URB,FLXHUMG_URB )

    TS_URB2D(I,J) = TS_URB

    ALBEDO(I,J)   = FRC_URB2D(I,J) * ALB_URB + (1-FRC_URB2D(I,J)) * ALBEDO(I,J)        ![-]
    HFX(I,J)      = FRC_URB2D(I,J) * SH_URB  + (1-FRC_URB2D(I,J)) * HFX(I,J)           ![W/m/m]
    QFX(I,J)      = FRC_URB2D(I,J) * LH_KINEMATIC_URB &
                       + (1-FRC_URB2D(I,J))* QFX(I,J)                                  ![kg/m/m/s]
    LH(I,J)       = FRC_URB2D(I,J) * LH_URB  + (1-FRC_URB2D(I,J)) * LH(I,J)            ![W/m/m]
    GRDFLX(I,J)   = FRC_URB2D(I,J) * (G_URB) + (1-FRC_URB2D(I,J)) * GRDFLX(I,J)        ![W/m/m]
    TSK(I,J)      = FRC_URB2D(I,J) * TS_URB  + (1-FRC_URB2D(I,J)) * TSK(I,J)           ![K]
!    Q1            = QSFC(I,J)/(1.0+QSFC(I,J))
!    Q1            = FRC_URB2D(I,J) * QS_URB  + (1-FRC_URB2D(I,J)) * Q1                 ![-]

! Convert QSFC back to mixing ratio

!    QSFC(I,J)     = Q1/(1.0-Q1)
                   QSFC(I,J)= FRC_URB2D(I,J)*QS_URB+(1-FRC_URB2D(I,J))*QSFC(I,J)               !!   QSFC(I,J)=QSFC1D
    UST(I,J)      = FRC_URB2D(I,J) * UST_URB + (1-FRC_URB2D(I,J)) * UST(I,J)     ![m/s]

! Renew Urban State Variables

    TR_URB2D(I,J) = TR_URB
    TB_URB2D(I,J) = TB_URB
    TG_URB2D(I,J) = TG_URB
    TC_URB2D(I,J) = TC_URB
    QC_URB2D(I,J) = QC_URB
    UC_URB2D(I,J) = UC_URB

    TGR_URB2D(I,J)     = TGR_URB
    CMCR_URB2D(I,J)    = CMCR_URB
    FLXHUMR_URB2D(I,J) = FLXHUMR_URB
    FLXHUMB_URB2D(I,J) = FLXHUMB_URB
    FLXHUMG_URB2D(I,J) = FLXHUMG_URB
    DRELR_URB2D(I,J)   = DRELR_URB
    DRELB_URB2D(I,J)   = DRELB_URB
    DRELG_URB2D(I,J)   = DRELG_URB

    DO K = 1,num_roof_layers
      TRL_URB3D(I,K,J) = TRL_URB(K)
      SMR_URB3D(I,K,J) = SMR_URB(K)
      TGRL_URB3D(I,K,J)= TGRL_URB(K)
    END DO
    DO K = 1,num_wall_layers
      TBL_URB3D(I,K,J) = TBL_URB(K)
    END DO
    DO K = 1,num_road_layers
      TGL_URB3D(I,K,J) = TGL_URB(K)
    END DO

    XXXR_URB2D(I,J)    = XXXR_URB
    XXXB_URB2D(I,J)    = XXXB_URB
    XXXG_URB2D(I,J)    = XXXG_URB
    XXXC_URB2D(I,J)    = XXXC_URB

    SH_URB2D(I,J)      = SH_URB
    LH_URB2D(I,J)      = LH_URB
    G_URB2D(I,J)       = G_URB
    RN_URB2D(I,J)      = RN_URB
    PSIM_URB2D(I,J)    = PSIM_URB
    PSIH_URB2D(I,J)    = PSIH_URB
    GZ1OZ0_URB2D(I,J)  = GZ1OZ0_URB
    U10_URB2D(I,J)     = U10_URB
    V10_URB2D(I,J)     = V10_URB
    TH2_URB2D(I,J)     = TH2_URB
    Q2_URB2D(I,J)      = Q2_URB
    UST_URB2D(I,J)     = UST_URB
    AKMS_URB2D(I,J)    = KARMAN * UST_URB2D(I,J)/(GZ1OZ0_URB2D(I,J)-PSIM_URB2D(I,J))
    IF (PRESENT(CMR_SFCDIF)) THEN
      CMR_SFCDIF(I,J)  = CMR_URB
      CHR_SFCDIF(I,J)  = CHR_URB
      CMGR_SFCDIF(I,J) = CMGR_URB
      CHGR_SFCDIF(I,J) = CHGR_URB
      CMC_SFCDIF(I,J)  = CMC_URB
      CHC_SFCDIF(I,J)  = CHC_URB
    ENDIF

  ENDIF                                 ! urban land used type block

ENDDO ILOOP                             ! of I loop
ENDDO JLOOP                             ! of J loop

ENDIF                                   ! sf_urban_physics = 1 block

!--------------------------------------
! URBAN CANOPY MODEL END
!--------------------------------------

!--------------------------------------
! URBAN BEP and BEM MODEL BEGIN
!--------------------------------------

IF (SF_URBAN_PHYSICS == 2) THEN

DO J=JTS,JTE
DO I=ITS,ITE

  EMISS_URB(I,J)       = 0.
  RL_UP_URB(I,J)       = 0.
  RS_ABS_URB(I,J)      = 0.
  GRDFLX_URB(I,J)      = 0.
  B_Q_BEP(I,KTS:KTE,J) = 0.

END DO
END DO

  CALL BEP(frc_urb2d,  utype_urb2d, itimestep,       dz8w,         &
                  dt,        u_phy,     v_phy,                     &
              th_phy,          rho,     p_phy,     swdown,    glw, &
                 gmt,       julday,     xlong,       xlat,         &
          declin_urb,   cosz_urb2d, omg_urb2d,                     &
       num_urban_ndm, urban_map_zrd, urban_map_zwd, urban_map_gd,  &
        urban_map_zd, urban_map_zdf,  urban_map_bd, urban_map_wd,  &
       urban_map_gbd, urban_map_fbd,  num_urban_hi,                &
           trb_urb4d,    tw1_urb4d, tw2_urb4d,  tgb_urb4d,         &
          sfw1_urb3d,   sfw2_urb3d, sfr_urb3d,  sfg_urb3d,         &
            lp_urb2d,     hi_urb2d,  lb_urb2d,  hgt_urb2d,         &
             a_u_bep,      a_v_bep,   a_t_bep,                     &
             a_e_bep,      b_u_bep,   b_v_bep,                     &
             b_t_bep,      b_e_bep,   b_q_bep,    dlg_bep,         &
            dl_u_bep,       sf_bep,    vl_bep,                     &
           rl_up_urb,   rs_abs_urb, emiss_urb, grdflx_urb,         &
         ids,ide, jds,jde, kds,kde,                                &
         ims,ime, jms,jme, kms,kme,                                &
         its,ite, jts,jte, kts,kte )

ENDIF ! SF_URBAN_PHYSICS == 2

IF (SF_URBAN_PHYSICS == 3) THEN

DO J=JTS,JTE
DO I=ITS,ITE

  EMISS_URB(I,J)       = 0.
  RL_UP_URB(I,J)       = 0.
  RS_ABS_URB(I,J)      = 0.
  GRDFLX_URB(I,J)      = 0.
  B_Q_BEP(I,KTS:KTE,J) = 0.

END DO
END DO

  CALL BEP_BEM( frc_urb2d,  utype_urb2d,    itimestep,         dz8w,       &
                       dt,        u_phy,        v_phy,                     &
                   th_phy,          rho,        p_phy,       swdown,  glw, &
                      gmt,       julday,        xlong,         xlat,       &
               declin_urb,   cosz_urb2d,    omg_urb2d,                     &
            num_urban_ndm, urban_map_zrd, urban_map_zwd, urban_map_gd,     &
             urban_map_zd, urban_map_zdf,  urban_map_bd, urban_map_wd,     &
            urban_map_gbd, urban_map_fbd,  urban_map_zgrd,num_urban_hi,    &
                trb_urb4d,    tw1_urb4d,    tw2_urb4d,    tgb_urb4d,       &
               tlev_urb3d,   qlev_urb3d, tw1lev_urb3d, tw2lev_urb3d,       &
              tglev_urb3d,  tflev_urb3d,  sf_ac_urb3d,  lf_ac_urb3d,       &
              cm_ac_urb3d, sfvent_urb3d, lfvent_urb3d,                     &
             sfwin1_urb3d, sfwin2_urb3d,                                   &
               sfw1_urb3d,   sfw2_urb3d,    sfr_urb3d,    sfg_urb3d,       &
              ep_pv_urb3d,   t_pv_urb3d,                                   & !RMS
                trv_urb4d,     qr_urb4d,    qgr_urb3d,   tgr_urb3d,        & !RMS
              drain_urb4d,draingr_urb3d,   sfrv_urb3d,  lfrv_urb3d,        & !RMS
                dgr_urb3d,     dg_urb3d,    lfr_urb3d,   lfg_urb3d,        & !RMS
                   rainbl,       swddir,       swddif,                     &
                 lp_urb2d,     hi_urb2d,     lb_urb2d,    hgt_urb2d,       &
                  a_u_bep,      a_v_bep,      a_t_bep,                     &
                  a_e_bep,      b_u_bep,      b_v_bep,                     &
                  b_t_bep,      b_e_bep,      b_q_bep,      dlg_bep,       &
                 dl_u_bep,       sf_bep,       vl_bep,                     &
                rl_up_urb,   rs_abs_urb,    emiss_urb,   grdflx_urb, qv3d, &
             ids,ide, jds,jde, kds,kde,                                    &
             ims,ime, jms,jme, kms,kme,                                    &
             its,ite, jts,jte, kts,kte )

ENDIF ! SF_URBAN_PHYSICS == 3

IF((SF_URBAN_PHYSICS == 2).OR.(SF_URBAN_PHYSICS == 3))THEN

  sigma_sb=5.67e-08
  do j = jts, jte
  do i = its, ite
    UMOM_URB(I,J)     = 0.
    VMOM_URB(I,J)     = 0.
    HFX_URB(I,J)      = 0.
    QFX_URB(I,J)      = 0.

    do k=kts,kte
      a_u_bep(i,k,j) = a_u_bep(i,k,j)*frc_urb2d(i,j)
      a_v_bep(i,k,j) = a_v_bep(i,k,j)*frc_urb2d(i,j)
      a_t_bep(i,k,j) = a_t_bep(i,k,j)*frc_urb2d(i,j)
      a_q_bep(i,k,j) = 0.
      a_e_bep(i,k,j) = 0.
      b_u_bep(i,k,j) = b_u_bep(i,k,j)*frc_urb2d(i,j)
      b_v_bep(i,k,j) = b_v_bep(i,k,j)*frc_urb2d(i,j)
      b_t_bep(i,k,j) = b_t_bep(i,k,j)*frc_urb2d(i,j)
      b_q_bep(i,k,j) = b_q_bep(i,k,j)*frc_urb2d(i,j)
      b_e_bep(i,k,j) = b_e_bep(i,k,j)*frc_urb2d(i,j)
      HFX_URB(I,J)   = HFX_URB(I,J) + B_T_BEP(I,K,J)*RHO(I,K,J)*CP*DZ8W(I,K,J)*VL_BEP(I,K,J)
      QFX_URB(I,J)   = QFX_URB(I,J) + B_Q_BEP(I,K,J)*DZ8W(I,K,J)*VL_BEP(I,K,J)
      UMOM_URB(I,J)  = UMOM_URB(I,J)+ (A_U_BEP(I,K,J)*U_PHY(I,K,J)+B_U_BEP(I,K,J))*DZ8W(I,K,J)*VL_BEP(I,K,J)
      VMOM_URB(I,J)  = VMOM_URB(I,J)+ (A_V_BEP(I,K,J)*V_PHY(I,K,J)+B_V_BEP(I,K,J))*DZ8W(I,K,J)*VL_BEP(I,K,J)
      vl_bep(i,k,j)  = (1.-frc_urb2d(i,j)) + vl_bep(i,k,j)*frc_urb2d(i,j)
      sf_bep(i,k,j)  = (1.-frc_urb2d(i,j)) + sf_bep(i,k,j)*frc_urb2d(i,j)
    end do

    a_u_bep(i,1,j)   = (1.-frc_urb2d(i,j))*(-ust(I,J)*ust(I,J))/dz8w(i,1,j)/   &
                          ((u_phy(i,1,j)**2+v_phy(i,1,j)**2.)**.5)+a_u_bep(i,1,j)

    a_v_bep(i,1,j)   = (1.-frc_urb2d(i,j))*(-ust(I,J)*ust(I,J))/dz8w(i,1,j)/   &
                          ((u_phy(i,1,j)**2+v_phy(i,1,j)**2.)**.5)+a_v_bep(i,1,j)

    b_t_bep(i,1,j)   = (1.-frc_urb2d(i,j))*hfx_rural(i,j)/dz8w(i,1,j)/rho(i,1,j)/CP+ &
                           b_t_bep(i,1,j)

    b_q_bep(i,1,j)   = (1.-frc_urb2d(i,j))*qfx_rural(i,j)/dz8w(i,1,j)/rho(i,1,j)+b_q_bep(i,1,j)

    umom             = (1.-frc_urb2d(i,j))*ust(i,j)*ust(i,j)*u_phy(i,1,j)/               &
                         ((u_phy(i,1,j)**2+v_phy(i,1,j)**2.)**.5)+umom_urb(i,j)

    vmom             = (1.-frc_urb2d(i,j))*ust(i,j)*ust(i,j)*v_phy(i,1,j)/               &
                         ((u_phy(i,1,j)**2+v_phy(i,1,j)**2.)**.5)+vmom_urb(i,j)
    sf_bep(i,1,j)    = 1.

! using the emissivity and the total longwave upward radiation estimate the averaged skin temperature

  IF (FRC_URB2D(I,J).GT.0.) THEN
    rl_up_rural   = -emiss_rural(i,j)*sigma_sb*(tsk_rural(i,j)**4.)-(1.-emiss_rural(i,j))*glw(i,j)
    rl_up_tot     = (1.-frc_urb2d(i,j))*rl_up_rural     + frc_urb2d(i,j)*rl_up_urb(i,j)
    emiss(i,j)    = (1.-frc_urb2d(i,j))*emiss_rural(i,j)+ frc_urb2d(i,j)*emiss_urb(i,j)
    ts_urb2d(i,j) = (max(0.,(-rl_up_urb(i,j)-(1.-emiss_urb(i,j))*glw(i,j))/emiss_urb(i,j)/sigma_sb))**0.25
    tsk(i,j)      = (max(0., (-1.*rl_up_tot-(1.-emiss(i,j))*glw(i,j) )/emiss(i,j)/sigma_sb))**.25
    rs_abs_tot    = (1.-frc_urb2d(i,j))*swdown(i,j)*(1.-albedo(i,j))+frc_urb2d(i,j)*rs_abs_urb(i,j)

    if(swdown(i,j) > 0.)then
      albedo(i,j) = 1.-rs_abs_tot/swdown(i,j)
    else
      albedo(i,j) = alb_rural(i,j)
    endif

! rename *_urb to sh_urb2d,lh_urb2d,g_urb2d,rn_urb2d

    grdflx(i,j)   = (1.-frc_urb2d(i,j))*grdflx_rural(i,j)+ frc_urb2d(i,j)*grdflx_urb(i,j)
    qfx(i,j)      = (1.-frc_urb2d(i,j))*qfx_rural(i,j)   + qfx_urb(i,j)
    lh(i,j)       = qfx(i,j)*xlv
    hfx(i,j)      = hfx_urb(i,j)                         + (1-frc_urb2d(i,j))*hfx_rural(i,j)      ![W/m/m]
    sh_urb2d(i,j) = hfx_urb(i,j)/frc_urb2d(i,j)
    lh_urb2d(i,j) = qfx_urb(i,j)*xlv/frc_urb2d(i,j)
    g_urb2d(i,j)  = grdflx_urb(i,j)
    rn_urb2d(i,j) = rs_abs_urb(i,j)+emiss_urb(i,j)*glw(i,j)-rl_up_urb(i,j)
    ust(i,j)      = (umom**2.+vmom**2.)**.25

  ELSE

    sh_urb2d(i,j)    = 0.
    lh_urb2d(i,j)    = 0.
    g_urb2d(i,j)     = 0.
    rn_urb2d(i,j)    = 0.

  ENDIF

  enddo ! jloop
  enddo ! iloop

ENDIF ! SF_URBAN_PHYSICS == 2 or 3

!--------------------------------------
! URBAN BEP and BEM MODEL END
!--------------------------------------


END SUBROUTINE noahmp_urban

!------------------------------------------------------------------------------------------
!



!============================================================================================
!
! subroutine lsm_mosaic_init: initialization of mosaic state variables Generated by Aaron A.
! added on 24 May 2022
!
!============================================================================================
   SUBROUTINE NOAHMP_MOSAIC_INIT(XLAND, SNOW , SNOWH , CANWAT , ISLTYP ,   IVGTYP, XLAT, &
                           TSLB , SMOIS , SH2O , DZS , FNDSOILW , FNDSNOWH ,             &
                           TSK, isnowxy , tvxy     ,tgxy     ,canicexy ,         TMN,     XICE,   &
                           canliqxy ,eahxy    ,tahxy    ,cmxy     ,chxy     ,                     &
                           fwetxy   ,sneqvoxy ,alboldxy ,qsnowxy, qrainxy, wslakexy, zwtxy, waxy, &
                           wtxy     ,tsnoxy   ,zsnsoxy  ,snicexy  ,snliqxy  ,lfmassxy ,rtmassxy , &
                           stmassxy ,woodxy   ,stblcpxy ,fastcpxy ,xsaixy   ,lai      ,           &
                           grainxy  ,gddxy    ,                                                   &
                           croptype ,cropcat  ,                      &
                           irnumsi, irnummi, irnumfi, irwatsi,                                    &
                           irwatmi, irwatfi, ireloss, irsivol,                                    &
                           irmivol, irfivol, irrsplh,                                             &
                    !jref:start
                           t2mvxy   ,t2mbxy   ,chstarxy,            &
                    !jref:end
                           num_soil_layers, restart,                 &
                           allowed_to_read , iopt_run,  iopt_crop, iopt_irr, iopt_irrm,          &
                           sf_urban_physics, ISWATER, ISICE,                        &  ! urban scheme
                           ISURBAN, fractional_seaice,                              &
                           ids,ide, jds,jde, kds,kde,                &
                           ims,ime, jms,jme, kms,kme,                &
                           its,ite, jts,jte, kts,kte,                &
                           smoiseq, smcwtdxy, rechxy, deeprechxy,    &
                           LANDUSEF, LANDUSEF2, NLCAT, IOPT_MOSAIC,                     &      ! Added by Aaron A. **IMPORTANT FOR THE FOR LOOPS OF FRACTIONAL LAND USE
                           mosaic_cat_index, mosaic_cat,                                         &      ! Added by Aaron A. **IMPORTANT FOR THE FOR LOOPS OF FRACTIONAL LAND USE
                           IOPT_HUE,                                                      &
                           TSK_mosaic, TSLB_mosaic, SMOIS_mosaic, SH2O_mosaic,                   &      ! Added by Aaron A.
                           CANWAT_mosaic, SNOW_mosaic, SNOWH_mosaic,                             &      ! Added by Aaron A.
                           isnowxy_mosaic, tvxy_mosaic, tgxy_mosaic, canicexy_mosaic,            &      ! Added by Aaron A.
                           canliqxy_mosaic, eahxy_mosaic, tahxy_mosaic,              &      ! Added by Aaron A.
                           cmxy_mosaic, chxy_mosaic, fwetxy_mosaic, sneqvoxy_mosaic,             &      ! Added by Aaron A.
                           alboldxy_mosaic, qsnowxy_mosaic,qrainxy_mosaic, wslakexy_mosaic, zwtxy_mosaic,       &      ! Added by Aaron A.
                           waxy_mosaic, wtxy_mosaic, tsnoxy_mosaic, zsnsoxy_mosaic,              &      ! Added by Aaron A.
                           snicexy_mosaic, snliqxy_mosaic, lfmassxy_mosaic, rtmassxy_mosaic,     &      ! Added by Aaron A.
                           stmassxy_mosaic, woodxy_mosaic, stblcpxy_mosaic, fastcpxy_mosaic,     &      ! Added by Aaron A.
                           xsaixy_mosaic, lai_mosaic,               &      ! Added by Aaron A.

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
                           DETENTION_STORAGEXY_mosaic,                       & ! HUE Specific Variables
                           Z0, ZNT_mosaic, Z0_mosaic)                                                   ! Added by Aaron A. )                                                            ! Added by Aaron A.             ! Added by Aaron A.

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! This initilization functions to first generate an ordering of fractioanl land use based on each grid cell.
  ! Then, the data that has been initilzied in the normal NOAHMP INIT call are fed used to set each of the
  ! sub-landuse layer to the the same constant. This allows for continuity between the different calls of NOAHMP
  ! which may need all of the above variables to initilize.
  !
  ! In addition, the NOAHMP mosiac initilization will initilize the urban mosaic variables that are needed to inilize
  ! our runs.
  !
  ! Format of the code is as follows:
  ! -> Variable Declirations
  ! -> Land Use For loop that is nearly identical to the one that Dan Li authored for the NOAH LSM
  !     -> This uses a bubble sorting, in which the fractional land use is ordered through comparisoin
  !     -> of elemtns of land use that are directly next to each other.
  !     -> There are a few checks that occur to ensure that the LAND_SEA mask and the mosaic are
  !     -> Consistent. See below for more detailed information
  ! -> Setting filling the mosaic variables!
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !This is the original call for the NOAHMP INIT

  USE NOAHMP_TABLES


  INTEGER, INTENT(IN   )    ::     ids,ide, jds,jde, kds,kde,  &
         &                           ims,ime, jms,jme, kms,kme,  &
         &                           its,ite, jts,jte, kts,kte

    INTEGER, INTENT(IN)       ::     num_soil_layers, iopt_run, iopt_crop, ISWATER,ISURBAN, ISICE, fractional_seaice

    LOGICAL, INTENT(IN)       ::     restart,                    &
         &                           allowed_to_read
    INTEGER, INTENT(IN)       ::     sf_urban_physics                              ! urban, by yizhou

    REAL,    DIMENSION( num_soil_layers), INTENT(IN)    ::     DZS  ! Thickness of the soil layers [m]

    REAL,    DIMENSION( ims:ime, num_soil_layers, jms:jme ) ,    &
         &   INTENT(INOUT)    ::     SMOIS,                      &
         &                           SH2O,                       &
         &                           TSLB

    REAL,    DIMENSION( ims:ime, jms:jme ) ,                     &
         &   INTENT(INOUT)    ::     SNOW,                       &
         &                           SNOWH,                      &
         &                           CANWAT,                     &
         &                           XLAND

    INTEGER, DIMENSION( ims:ime, jms:jme ),                      &
         &   INTENT(IN)       ::     ISLTYP,  &
                                     IVGTYP

     LOGICAL, INTENT(IN)       ::     FNDSOILW,                   &
         &                           FNDSNOWH


    REAL, DIMENSION(ims:ime,jms:jme), INTENT(IN) :: XLAT         !latitude
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(IN) :: TSK         !skin temperature (k)
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: TMN         !deep soil temperature (k)
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(IN) :: XICE         !sea ice fraction
    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: isnowxy     !actual no. of snow layers
    REAL, DIMENSION(ims:ime,-2:num_soil_layers,jms:jme), INTENT(INOUT) :: zsnsoxy  !snow layer depth [m]
    REAL, DIMENSION(ims:ime,-2:              0,jms:jme), INTENT(INOUT) :: tsnoxy   !snow temperature [K]
    REAL, DIMENSION(ims:ime,-2:              0,jms:jme), INTENT(INOUT) :: snicexy  !snow layer ice [mm]
    REAL, DIMENSION(ims:ime,-2:              0,jms:jme), INTENT(INOUT) :: snliqxy  !snow layer liquid water [mm]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: tvxy        !vegetation canopy temperature
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: tgxy        !ground surface temperature
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: canicexy    !canopy-intercepted ice (mm)
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: canliqxy    !canopy-intercepted liquid water (mm)
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: eahxy       !canopy air vapor pressure (pa)
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: tahxy       !canopy air temperature (k)
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: cmxy        !momentum drag coefficient
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: chxy        !sensible heat exchange coefficient
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: fwetxy      !wetted or snowed fraction of the canopy (-)
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: sneqvoxy    !snow mass at last time step(mm h2o)
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: alboldxy    !snow albedo at last time step (-)
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: qsnowxy     !snowfall on the ground [mm/s]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: qrainxy     !rainfall on the ground [mm/s]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: wslakexy    !lake water storage [mm]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: zwtxy       !water table depth [m]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: waxy        !water in the "aquifer" [mm]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: wtxy        !groundwater storage [mm]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: lfmassxy    !leaf mass [g/m2]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: rtmassxy    !mass of fine roots [g/m2]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: stmassxy    !stem mass [g/m2]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: woodxy      !mass of wood (incl. woody roots) [g/m2]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: grainxy     !mass of grain [g/m2] !XING
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: gddxy       !growing degree days !XING
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: stblcpxy    !stable carbon in deep soil [g/m2]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: fastcpxy    !short-lived carbon, shallow soil [g/m2]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: xsaixy      !stem area index
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: lai         !leaf area index

!jref:start
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: t2mvxy        !2m temperature vegetation part (k)
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: t2mbxy        !2m temperature bare ground part (k)
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: chstarxy        !dummy
!jref:end

    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: irnumsi
    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: irnummi
    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: irnumfi
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: irwatsi
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: irwatmi
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: irwatfi
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: ireloss
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: irsivol
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: irmivol
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: irfivol
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: irrsplh

    INTEGER, DIMENSION(ims:ime,  jms:jme), INTENT(OUT) :: cropcat
    REAL   , DIMENSION(ims:ime,5,jms:jme), INTENT(IN ) :: croptype

    REAL, DIMENSION(ims:ime,1:num_soil_layers,jms:jme), INTENT(INOUT) , OPTIONAL :: smoiseq !equilibrium soil moisture content [m3m-3]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) , OPTIONAL :: smcwtdxy    !deep soil moisture content [m3m-3]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) , OPTIONAL :: deeprechxy  !deep recharge [m]
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) , OPTIONAL :: rechxy      !accumulated recharge [mm]
!added by Aaron A.
! From here on, we define the important variables that are needed to re-catagorize the landuse data sets

  INTEGER, INTENT(IN) :: NLCAT
  INTEGER, INTENT(IN) :: IOPT_MOSAIC                                       !This tells us that the mosaic scheme is active
  INTEGER, INTENT(IN) :: mosaic_cat                                               !This tells how many mosaic catagories that are of interet,
  REAL, DIMENSION( ims:ime, 1:NLCAT, jms:jme ) , INTENT(IN)::   LANDUSEF            !This is the original land-use fraction that was read in
  REAL, DIMENSION( ims:ime, 1:NLCAT, jms:jme ) , INTENT(INOUT)::   LANDUSEF2        !This is the land-use fraction that has been re-ordered

  INTEGER, DIMENSION( ims:ime, NLCAT, jms:jme ), INTENT(INOUT) :: mosaic_cat_index !This is the re-ordered mosaic catagory data

  !variables with dimensions 1:mosaic_cat
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: TSK_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: CANWAT_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: SNOW_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: SNOWH_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: tvxy_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: tgxy_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: canicexy_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: canliqxy_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: eahxy_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: tahxy_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: cmxy_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: chxy_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: fwetxy_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: sneqvoxy_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: alboldxy_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: qsnowxy_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: qrainxy_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: wslakexy_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: zwtxy_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: waxy_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: wtxy_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: lfmassxy_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: rtmassxy_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: stmassxy_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: woodxy_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: stblcpxy_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: fastcpxy_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: xsaixy_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: lai_mosaic


 !snow variables, which have dimensions of 7*mosaic cat and 3 * number of mosaic cats
  INTEGER, DIMENSION(ims:ime, 1:mosaic_cat, jms:jme), OPTIONAL, INTENT(INOUT) :: isnowxy_mosaic     !actual no. of snow layers
  REAL, DIMENSION(ims:ime, 1:7*mosaic_cat, jms:jme), OPTIONAL,  INTENT(INOUT) :: zsnsoxy_mosaic  !snow layer depth [m] **These have not been adjusted for indexing, Added in main driver
  REAL, DIMENSION(ims:ime, 1:3*mosaic_cat, jms:jme), OPTIONAL,  INTENT(INOUT) :: tsnoxy_mosaic   !snow temperature [K] **These have not been adjusted for indexing, Added in main driver
  REAL, DIMENSION(ims:ime, 1:3*mosaic_cat, jms:jme), OPTIONAL,  INTENT(INOUT) :: snicexy_mosaic  !snow layer ice [mm] **These have not been adjusted for indexing, Added in main driver
  REAL, DIMENSION(ims:ime, 1:3*mosaic_cat, jms:jme), OPTIONAL,  INTENT(INOUT) :: snliqxy_mosaic  !snow layer liquid water [mm] **These have not been adjusted for indexing, Added in main driver

 !variables that are for the soil layers
  REAL, DIMENSION( ims:ime, 1:num_soil_layers*mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT):: TSLB_mosaic
  REAL, DIMENSION( ims:ime, 1:num_soil_layers*mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT):: SMOIS_mosaic
  REAL, DIMENSION( ims:ime, 1:num_soil_layers*mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT):: SH2O_mosaic

  !values that are needed for the
  REAL, OPTIONAL, DIMENSION( ims:ime, 1:num_soil_layers*mosaic_cat, jms:jme ), INTENT(INOUT) :: TRL_URB3D_mosaic
  REAL, OPTIONAL, DIMENSION( ims:ime, 1:num_soil_layers*mosaic_cat, jms:jme ), INTENT(INOUT) :: TBL_URB3D_mosaic
  REAL, OPTIONAL, DIMENSION( ims:ime, 1:num_soil_layers*mosaic_cat, jms:jme ), INTENT(INOUT) :: TGL_URB3D_mosaic

  ! Needed for future integration of
  REAL, DIMENSION(ims:ime,1:num_soil_layers*mosaic_cat,jms:jme), INTENT(INOUT) , OPTIONAL :: smoiseq_mosaic ! equilibrium soil moisture content [m3m-3]
  REAL, DIMENSION(ims:ime,1:mosaic_cat, jms:jme), INTENT(INOUT) , OPTIONAL :: smcwtdxy_mosaic    !deep soil moisture content [m3m-3]
  REAL, DIMENSION(ims:ime,1:mosaic_cat, jms:jme), INTENT(INOUT) , OPTIONAL :: deeprechxy_mosaic  !deep recharge [m]
  REAL, DIMENSION(ims:ime,1:mosaic_cat, jms:jme), INTENT(INOUT) , OPTIONAL :: rechxy_mosaic      !accu

  ! Needed for irrigation model
  INTEGER, DIMENSION(ims:ime,1:mosaic_cat,jms:jme), INTENT(INOUT) :: IRNUMSI_mosaic
  INTEGER, DIMENSION(ims:ime,1:mosaic_cat,jms:jme), INTENT(INOUT) :: IRNUMMI_mosaic
  INTEGER, DIMENSION(ims:ime,1:mosaic_cat,jms:jme), INTENT(INOUT) :: IRNUMFI_mosaic
  REAL,    DIMENSION(ims:ime,1:mosaic_cat,jms:jme), INTENT(INOUT) :: IRWATSI_mosaic
  REAL,    DIMENSION(ims:ime,1:mosaic_cat,jms:jme), INTENT(INOUT) :: IRWATMI_mosaic
  REAL,    DIMENSION(ims:ime,1:mosaic_cat,jms:jme), INTENT(INOUT) :: IRWATFI_mosaic
  REAL,    DIMENSION(ims:ime,1:mosaic_cat,jms:jme), INTENT(INOUT) :: IRELOSS_mosaic
  REAL,    DIMENSION(ims:ime,1:mosaic_cat,jms:jme), INTENT(INOUT) :: IRSIVOL_mosaic
  REAL,    DIMENSION(ims:ime,1:mosaic_cat,jms:jme), INTENT(INOUT) :: IRMIVOL_mosaic
  REAL,    DIMENSION(ims:ime,1:mosaic_cat,jms:jme), INTENT(INOUT) :: IRFIVOL_mosaic
  REAL,    DIMENSION(ims:ime,1:mosaic_cat,jms:jme), INTENT(INOUT) :: IRRSPLH_mosaic


  !variables for 2d urban model
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: TR_URB2D_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: TB_URB2D_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: TG_URB2D_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: TC_URB2D_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: QC_URB2D_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: SH_URB2D_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: LH_URB2D_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: G_URB2D_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: RN_URB2D_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: TS_URB2D_mosaic

  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: CMR_SFCDIF_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: CHR_SFCDIF_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: CMC_SFCDIF_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: CHC_SFCDIF_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: CMGR_SFCDIF_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: CHGR_SFCDIF_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: XXXR_URB2D_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: XXXB_URB2D_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: XXXG_URB2D_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: XXXC_URB2D_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: CMCR_URB2D_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: TGR_URB2D_mosaic

  REAL, DIMENSION( ims:ime, 1:num_soil_layers*mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: TGRL_URB3D_mosaic
  REAL, DIMENSION( ims:ime, 1:num_soil_layers*mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: SMR_URB3D_mosaic

  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: DRELR_URB2D_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: DRELB_URB2D_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: DRELG_URB2D_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: FLXHUMR_URB2D_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: FLXHUMB_URB2D_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: FLXHUMG_URB2D_mosaic
  ! HUE Variables
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT) :: DETENTION_STORAGEXY_mosaic

  !variables for surface roughness
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme) , OPTIONAL, INTENT(INOUT) :: Z0_mosaic
  REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme) , OPTIONAL, INTENT(INOUT) :: ZNT_mosaic
  REAL, DIMENSION( ims:ime, jms:jme) , OPTIONAL, INTENT(INOUT) :: Z0



  INTEGER :: ij, i,j,mosaic_i,LastSwap,NumPairs,soil_k, Temp2,Temp5,Temp7, ICE,temp_index
  REAL :: Temp, Temp3,Temp4,Temp6,xice_threshold,x
  INTEGER, DIMENSION(7) :: HUE_noahmp_extracatagories ! added by Aaron A.
  INTEGER :: LC, CAT_interested, temp_master_index
  REAL :: temp_master_LUF
  LOGICAL :: IPRINT
  CHARACTER(len=256) :: message_text


  IPRINT=.false.


  if ( fractional_seaice == 0 ) then
     xice_threshold = 0.5
  else if ( fractional_seaice == 1 ) then
     xice_threshold = 0.02
  endif


  !===========================================================================
  ! CHOOSE THE TILES
  !===========================================================================

  itf=min0(ite,ide-1)
  jtf=min0(jte,jde-1)

  ! simple test

  DO i = its,itf
     DO j = jts,jtf
        IF ((xland(i,j).LT. 1.5 ) .AND. (IVGTYP(i,j) .EQ. ISWATER)) THEN
           PRINT*, 'BEFORE MOSAIC_INIT'
           CALL wrf_message("BEFORE MOSAIC_INIT")
           WRITE(message_text,fmt='(a,2I6,2F8.2,2I6)') 'I,J,xland,xice,mosaic_cat_index,ivgtyp = ', &
                 I,J,xland(i,j),xice(i,j),mosaic_cat_index(I,1,J),IVGTYP(i,j)
           CALL wrf_message(message_text)
        ENDIF
     ENDDO
  ENDDO

     DO i = its,itf                               !This loop populates the LANDUSEF2 and mosaic_cat
        DO j = jts,jtf
           DO mosaic_i=1,NLCAT
              LANDUSEF2(i,mosaic_i,j)=LANDUSEF(i,mosaic_i,j)
              mosaic_cat_index(i,mosaic_i,j)=mosaic_i

           ENDDO
        ENDDO
     ENDDO

     DO i = its,itf  !These do loops swap around the fractional land use
        DO j = jts,jtf

          NumPairs=NLCAT-1

          DO
               IF (NumPairs == 0) EXIT
                   LastSwap = 1
          DO  mosaic_i=1, NumPairs
            IF(LANDUSEF2(i,mosaic_i, j) < LANDUSEF2(i,mosaic_i+1, j)  ) THEN
               Temp = LANDUSEF2(i,mosaic_i, j)
               LANDUSEF2(i,mosaic_i, j)=LANDUSEF2(i,mosaic_i+1, j)
               LANDUSEF2(i,mosaic_i+1, j)=Temp
               LastSwap = mosaic_i

               Temp2 =  mosaic_cat_index(i,mosaic_i,j)
               mosaic_cat_index(i,mosaic_i,j)=mosaic_cat_index(i,mosaic_i+1,j)
               mosaic_cat_index(i,mosaic_i+1,j)=Temp2
            ENDIF
          ENDDO
               NumPairs = LastSwap - 1
          ENDDO

        ENDDO
      ENDDO
  !===========================================================================
  ! For non-seaice grids, eliminate the seaice-tiles
  !===========================================================================

     DO i = its,itf
        DO j = jts,jtf

         IF   (XLAND(I,J).LT.1.5)  THEN

             ICE = 0
                 IF( XICE(I,J).GE. XICE_THRESHOLD ) THEN
                   WRITE (message_text,fmt='(a,2I5)') 'sea-ice at point, I and J = ', i,j
                   CALL wrf_message(message_text)
                 ICE = 1
                 ENDIF

          IF (ICE == 1)   Then         ! sea-ice case , eliminate sea-ice if they are not the dominant ones

          IF (IVGTYP(i,j) == isice)  THEN    ! if this grid cell is dominanted by ice, then do nothing

          ELSE

                DO mosaic_i=2,mosaic_cat
                   IF (mosaic_cat_index(i,mosaic_i,j) == isice ) THEN
                       Temp4=LANDUSEF2(i,mosaic_i,j)
                       Temp5=mosaic_cat_index(i,mosaic_i,j)

                       LANDUSEF2(i,mosaic_i:NLCAT-1,j)=LANDUSEF2(i,mosaic_i+1:NLCAT,j)
                       mosaic_cat_index(i,mosaic_i:NLCAT-1,j)=mosaic_cat_index(i,mosaic_i+1:NLCAT,j)

                       LANDUSEF2(i,NLCAT,j)=Temp4
                       mosaic_cat_index(i,NLCAT,j)=Temp5
                   ENDIF
                 ENDDO

          ENDIF   ! for (IVGTYP(i,j) == isice )

          ELSEIF (ICE ==0)  THEN

          IF ((mosaic_cat_index(I,1,J) .EQ. ISWATER)) THEN

          ! xland < 1.5 but the dominant land use category based on our calculation is water

           IF (IVGTYP(i,j) .EQ. ISWATER) THEN

           ! xland < 1.5 but the dominant land use category based on the geogrid calculation is water, this must be wrong

              CALL wrf_message("IN MOSAIC_INIT")
              WRITE(message_text,fmt='(a,3I6,2F8.2)') 'I,J,IVGTYP,XLAND,XICE = ',I,J,IVGTYP(I,J),xland(i,j),xice(i,j)
              CALL wrf_message(message_text)
              CALL wrf_message("xland < 1.5 but the dominant land use category based on our calculation is water."//&
                   "In addition, the dominant land use category based on the geogrid calculation is water, this must be wrong")

           ENDIF  ! for (IVGTYP(i,j) .EQ. ISWATER)

           IF (IVGTYP(i,j) .NE. ISWATER) THEN

           ! xland < 1.5,   the dominant land use category based on our calculation is water, but based on the geogrid calculation is not water, which might be due to the inconsistence between land use data and land-sea mask

           Temp4=LANDUSEF2(i,1,j)
           Temp5=mosaic_cat_index(i,1,j)

           LANDUSEF2(i,1:NLCAT-1,j)=LANDUSEF2(i,2:NLCAT,j)
           mosaic_cat_index(i,1:NLCAT-1,j)=mosaic_cat_index(i,2:NLCAT,j)

           LANDUSEF2(i,NLCAT,j)=Temp4
           mosaic_cat_index(i,NLCAT,j)=Temp5

              CALL wrf_message("IN MOSAIC_INIT")
              WRITE(message_text,fmt='(a,3I6,2F8.2)') 'I,J,IVGTYP,XLAND,XICE = ',I,J,IVGTYP(I,J),xland(i,j),xice(i,j)
              CALL wrf_message(message_text)
              CALL wrf_message("xland < 1.5 but the dominant land use category based on our calculation is water."//&
                   "this is fine as long as we change our calculation so that the dominant land use category is"//&
                   "stwiched back to not water.")
              WRITE(message_text,fmt='(a,2I6)') 'land use category has been switched, before and after values are ', &
                   temp5,mosaic_cat_index(i,1,j)
              CALL wrf_message(message_text)
              WRITE(message_text,fmt='(a,2I6)') 'new dominant and second dominant cat are ', mosaic_cat_index(i,1,j),mosaic_cat_index(i,2,j)
              CALL wrf_message(message_text)

           ENDIF  ! for (IVGTYP(i,j) .NE. ISWATER)

           ELSE    !  for (mosaic_cat_index(I,1,J) .EQ. ISWATER)

                     DO mosaic_i=2,mosaic_cat
                    IF (mosaic_cat_index(i,mosaic_i,j) == iswater ) THEN
                       Temp4=LANDUSEF2(i,mosaic_i,j)
                       Temp5=mosaic_cat_index(i,mosaic_i,j)

                       LANDUSEF2(i,mosaic_i:NLCAT-1,j)=LANDUSEF2(i,mosaic_i+1:NLCAT,j)
                       mosaic_cat_index(i,mosaic_i:NLCAT-1,j)=mosaic_cat_index(i,mosaic_i+1:NLCAT,j)

                       LANDUSEF2(i,NLCAT,j)=Temp4
                       mosaic_cat_index(i,NLCAT,j)=Temp5
                    ENDIF
                  ENDDO

           ENDIF !  for (mosaic_cat_index(I,1,J) .EQ. ISWATER)

          ENDIF  !  for ICE == 1

      ELSE  ! FOR (XLAND(I,J).LT.1.5)

                 ICE = 0

                     IF( XICE(I,J).GE. XICE_THRESHOLD ) THEN
                       WRITE (message_text,fmt='(a,2I6)') 'sea-ice at water point, I and J = ', i,j
                       CALL wrf_message(message_text)
                       ICE = 1
                     ENDIF

           IF ((mosaic_cat_index(I,1,J) .NE. ISWATER)) THEN

                ! xland > 1.5 and the dominant land use category based on our calculation is not water

                 IF (IVGTYP(i,j) .NE. ISWATER) THEN

                 ! xland > 1.5 but the dominant land use category based on the geogrid calculation is not water, this must be wrong
                 CALL wrf_message("IN MOSAIC_INIT")
                 WRITE(message_text,fmt='(a,3I6,2F8.2)') 'I,J,IVGTYP,XLAND,XICE = ',I,J,IVGTYP(I,J),xland(i,j),xice(i,j)
                 CALL wrf_message(message_text)
                 CALL wrf_message("xland > 1.5 but the dominant land use category based on our calculation is not water."// &
                      "in addition, the dominant land use category based on the geogrid calculation is not water,"//  &
                      "this must be wrong.")
                 ENDIF  ! for (IVGTYP(i,j) .NE. ISWATER)

                 IF (IVGTYP(i,j) .EQ. ISWATER) THEN

                 ! xland > 1.5,   the dominant land use category based on our calculation is not water, but based on the geogrid calculation is water, which might be due to the inconsistence between land use data and land-sea mask

                 CALL wrf_message("IN MOSAIC_INIT")
                 WRITE(message_text,fmt='(a,3I6,2F8.2)') 'I,J,IVGTYP,XLAND,XICE = ',I,J,IVGTYP(I,J),xland(i,j),xice(i,j)
                 CALL wrf_message(message_text)
                 CALL wrf_message("xland > 1.5 but the dominant land use category based on our calculation is not water."// &
                      "however, the dominant land use category based on the geogrid calculation is water")
                 CALL wrf_message("This is fine. We do not need to do anyting because in the noaddrv, "//&
                      "we use xland as a criterion for whether using"// &
                      "mosaic or not when xland > 1.5, no mosaic will be used anyway")

                 ENDIF  ! for (IVGTYP(i,j) .NE. ISWATER)

           ENDIF !  for (mosaic_cat_index(I,1,J) .NE. ISWATER)

        ENDIF  ! FOR (XLAND(I,J).LT.1.5)

          ENDDO
      ENDDO

      !===========================================================================
      ! Swap the data so that 'paired' land-types are. We are going to
      ! order them in a specific way.
      ! Ordered so:
      ! Urban Tree over pavement (41)
      ! Urban Tree (47)
      ! Urban To Turfgrass (42)
      ! Urban Turfgrass (43)
      ! Urban to Permeable Pavement (45)
      ! Urban Permeable Pavement (44)
      ! Green Roof (46)
      ! EVERYTHING ELSE
      ! This assumes that you have done the work to add in extra grids
      ! Note we do not need to re-organize 46, because there
      ! is no paired land use
      !===========================================================================
      IF (iopt_HUE.eq.1) THEN ! this means that there are extra land-use values
        DO i = its,itf ! Loop over the i tile
          do j = jts, jtf ! loop over the j tile

            !set our indicies
            index_41 = -1
            index_47 = -1
            index_42 = -1
            index_43 = -1
            index_45 = -1
            index_44 = -1

            ! first find the locations of the indexes
            DO  mosaic_i=1, NLCAT
              ! first reorganizing 41 and 47
              IF (mosaic_cat_index(i,mosaic_i,j).eq.41) THEN
                index_41 = mosaic_i

              elseif (mosaic_cat_index(i,mosaic_i,j).eq.47) then
                index_47 = mosaic_i

              ENDIF
            ENDDO

            ! check to see if these are still negavitve. If they aren't, we can
            ! re-organize

            if ((index_41.gt.0).and.(index_47.gt.0) ) then !we can reorganize

              if ( index_41.lt.index_47 ) then !means that we are 41.....47, move together

                temp_master_index= mosaic_cat_index(i,index_47,j)
                mosaic_cat_index(i,index_41+2:index_47+1,j) = mosaic_cat_index(i,index_41+1:index_47,j)
                mosaic_cat_index(i,index_41+1,j) = temp_master_index

                temp_master_LUF = LANDUSEF2(i,index_47,j)
                LANDUSEF2(i,index_41+2:index_47+1,j) = LANDUSEF2(i,index_41+1:index_47,j)
                LANDUSEF2(i,index_41+1,j) = temp_master_LUF
              ELSE !This means we are ordered 47 .... 41, so we need to move and then swap

                temp_master_index = mosaic_cat_index(i,index_41,j)
                mosaic_cat_index(i,index_47+2:index_41+1,j) = mosaic_cat_index(i,index_47+1:index_41,j)
                mosaic_cat_index(i,index_47+1,j) = temp_master_index

                temp_master_LUF = LANDUSEF2(i,index_41,j)
                LANDUSEF2(i,index_47+2:index_41+1,j) = LANDUSEF2(i,index_41+1:index_47,j)
                LANDUSEF2(i,index_47+1,j) = temp_master_LUF

                ! now swap so it goes 41 47 instead of 47 41

                temp_master_index = mosaic_cat_index(i,index_47,j)
                mosaic_cat_index(i,index_47,j) = mosaic_cat_index(i,index_47+1,j)
                mosaic_cat_index(i,index_47+1,j) = temp_master_index

                temp_master_LUF = LANDUSEF2(i,index_47,j)
                LANDUSEF2(i,index_47,j) = LANDUSEF2(i,index_47+1,j)
                LANDUSEF2(i,index_47+1,j) = temp_master_LUF

              end if

            end if ! end the urban tree loop

      ! repeat process because index locations of 42 and 43 may have changed!

            ! first find the locations of the indexes for 42 and 43
            DO  mosaic_i=1, NLCAT
              ! first reorganizing 42 and 43
              IF (mosaic_cat_index(i,mosaic_i,j).eq.42) THEN
                index_42 = mosaic_i

              elseif (mosaic_cat_index(i,mosaic_i,j).eq.43) then
                index_43 = mosaic_i

              ENDIF
            ENDDO

            ! check to see if these are still negavitve. If they aren't, we can
            ! re-organize

            if ((index_42.gt.0).and.(index_43.gt.0) ) then !we can reorganize

              if ( index_42.lt.index_43 ) then !means that we are 42.....43, move together

                temp_master_index= mosaic_cat_index(i,index_43,j)
                mosaic_cat_index(i,index_42+2:index_43+1,j) = mosaic_cat_index(i,index_42+1:index_43,j)
                mosaic_cat_index(i,index_42+1,j) = temp_master_index

                temp_master_LUF = LANDUSEF2(i,index_43,j)
                LANDUSEF2(i,index_42+2:index_43+1,j) = LANDUSEF2(i,index_42+1:index_43,j)
                LANDUSEF2(i,index_42+1,j) = temp_master_LUF
              ELSE !This means we are ordered 43 .... 42, so we need to move and then swap

                temp_master_index = mosaic_cat_index(i,index_42,j)
                mosaic_cat_index(i,index_43+2:index_42+1,j) = mosaic_cat_index(i,index_43+1:index_42,j)
                mosaic_cat_index(i,index_43+1,j) = temp_master_index

                temp_master_LUF = LANDUSEF2(i,index_42,j)
                LANDUSEF2(i,index_43+2:index_42+1,j) = LANDUSEF2(i,index_43+1:index_42,j)
                LANDUSEF2(i,index_43+1,j) = temp_master_LUF

                ! now swap so it goes 42 43 instead of 43 42

                temp_master_index = mosaic_cat_index(i,index_43,j)
                mosaic_cat_index(i,index_43,j) = mosaic_cat_index(i,index_43+1,j)
                mosaic_cat_index(i,index_43+1,j) = temp_master_index

                temp_master_LUF = LANDUSEF2(i,index_43,j)
                LANDUSEF2(i,index_43,j) = LANDUSEF2(i,index_43+1,j)
                LANDUSEF2(i,index_43+1,j) = temp_master_LUF

              end if

            end if ! end urban disconnection to turfgrass

      ! repeat process because index locations of 45 and 44 may have changed!

            ! first find the locations of the indexes for 42 and 43
              DO  mosaic_i=1, NLCAT
              ! first reorganizing 42 and 43
                  IF (mosaic_cat_index(i,mosaic_i,j).eq.45) THEN
                    index_45 = mosaic_i

                  elseif (mosaic_cat_index(i,mosaic_i,j).eq.44) then
                    index_44 = mosaic_i

                  ENDIF
                ENDDO

              ! check to see if these are still negavitve. If they aren't, we can
              ! re-organize

              if ((index_45.gt.0).and.(index_44.gt.0) ) then !we can reorganize

                if ( index_45.lt.index_44) then !means that we are 45.....44, move together

                  temp_master_index= mosaic_cat_index(i,index_44,j)
                  mosaic_cat_index(i,index_45+2:index_44+1,j) = mosaic_cat_index(i,index_45+1:index_44,j)
                  mosaic_cat_index(i,index_45+1,j) = temp_master_index

                  temp_master_LUF = LANDUSEF2(i,index_44,j)
                  LANDUSEF2(i,index_45+2:index_44+1,j) = LANDUSEF2(i,index_45+1:index_44,j)
                  LANDUSEF2(i,index_45+1,j) = temp_master_LUF

                ELSE !This means we are ordered 44 .... 45, so we need to move and then swap

                  temp_master_index = mosaic_cat_index(i,index_45,j)
                  mosaic_cat_index(i,index_44+2:index_45+1,j) = mosaic_cat_index(i,index_44+1:index_45,j)
                  mosaic_cat_index(i,index_44+1,j) = temp_master_index

                  temp_master_LUF = LANDUSEF2(i,index_45,j)
                  LANDUSEF2(i,index_44+2:index_45+1,j) = LANDUSEF2(i,index_44+1:index_45,j)
                  LANDUSEF2(i,index_44+1,j) = temp_master_LUF

                  ! now swap so it goes 45 44 instead of 44 45

                  temp_master_index = mosaic_cat_index(i,index_44,j)
                  mosaic_cat_index(i,index_44,j) = mosaic_cat_index(i,index_44+1,j)
                  mosaic_cat_index(i,index_44+1,j) = temp_master_index

                  temp_master_LUF = LANDUSEF2(i,index_44,j)
                  LANDUSEF2(i,index_44,j) = LANDUSEF2(i,index_44+1,j)
                  LANDUSEF2(i,index_44+1,j) = temp_master_LUF

                end if

              end if ! end urban to permeable pavements
      ENDDO
      ENDDO
      ENDIF !end hue noahmp sort
      !===========================================================================
      ! normalize
      !===========================================================================

     DO i = its,itf
        DO j = jts,jtf

          Temp6=0

            DO mosaic_i=1,mosaic_cat
               Temp6=Temp6+LANDUSEF2(i,mosaic_i,j)
            ENDDO

            if (Temp6 .LT. 1e-5)  then

            Temp6 = 1e-5
            WRITE (message_text,fmt='(a,e8.1)') 'the total land surface fraction is less than ', temp6
            CALL wrf_message(message_text)
            WRITE (message_text,fmt='(a,2I6,4F8.2)') 'some landusef values at i,j are ', &
                 i,j,landusef2(i,1,j),landusef2(i,2,j),landusef2(i,3,j),landusef2(i,4,j)
            CALL wrf_message(message_text)
            WRITE (message_text,fmt='(a,2I6,3I6)') 'some mosaic cat values at i,j are ', &
                 i,j,mosaic_cat_index(i,1,j),mosaic_cat_index(i,2,j),mosaic_cat_index(i,3,j)
            CALL wrf_message(message_text)

            endif

            LANDUSEF2(i,1:mosaic_cat, j)=LANDUSEF2(i,1:mosaic_cat,j)*(1/Temp6)

          ENDDO
      ENDDO
WRITE(*,*) "POST NORMALIZATION"
  !===========================================================================
  ! initilize the variables
  !===========================================================================
IF(.not.restart)THEN !! only do these if this is a restart, otherwise,we have
  ! things already generated/read in
     DO i = its,itf
        DO j = jts,jtf
               
             DO mosaic_i=1,mosaic_cat

                !add qualifer for the mosaic croptpes
            ZNT_mosaic(i,mosaic_i,j)=Z0(i,j) !roughness following Dan Li in Noah model
            Z0_mosaic(i,mosaic_i,j)=Z0(i,j) !roughness following Dan LI in Noah Model

            TSK_mosaic(i,mosaic_i,j)=TSK(i,j)
            CANWAT_mosaic(i,mosaic_i,j)=CANWAT(i,j)
            SNOW_mosaic(i,mosaic_i,j)=SNOW(i,j)
            SNOWH_mosaic(i,mosaic_i,j)=SNOWH(i,j)

            tvxy_mosaic(i,mosaic_i,j)=tvxy(i,j)
            tgxy_mosaic(i,mosaic_i,j)=tgxy(i,j)

            canicexy_mosaic(i,mosaic_i,j)=canicexy(i,j)
            canliqxy_mosaic(i,mosaic_i,j)=canliqxy(i,j)

            eahxy_mosaic(i,mosaic_i,j)=eahxy(i,j)
            tahxy_mosaic(i,mosaic_i,j)=tahxy(i,j)
            cmxy_mosaic(i,mosaic_i,j)=cmxy(i,j)
            chxy_mosaic(i,mosaic_i,j)=chxy(i,j)

            fwetxy_mosaic(i,mosaic_i,j)=fwetxy(i,j)
            sneqvoxy_mosaic(i,mosaic_i,j)=sneqvoxy(i,j)
            alboldxy_mosaic(i,mosaic_i,j)=alboldxy(i,j)
            qsnowxy_mosaic(i,mosaic_i,j)=qsnowxy(i,j)
            qrainxy_mosaic(i,mosaic_i,j)=qrainxy(i,j)

            wslakexy_mosaic(i,mosaic_i,j)=wslakexy(i,j)
            zwtxy_mosaic(i,mosaic_i,j)=zwtxy(i,j)
            waxy_mosaic(i,mosaic_i,j)=waxy(i,j)
            wtxy_mosaic(i,mosaic_i,j)=wtxy(i,j)

            lfmassxy_mosaic(i,mosaic_i,j)=lfmassxy(i,j)
            rtmassxy_mosaic(i,mosaic_i,j)=rtmassxy(i,j)
            stmassxy_mosaic(i,mosaic_i,j)=stmassxy(i,j)

            woodxy_mosaic(i,mosaic_i,j)=woodxy(i,j)
            stblcpxy_mosaic(i,mosaic_i,j)=stblcpxy(i,j)
            fastcpxy_mosaic(i,mosaic_i,j)=fastcpxy(i,j)
            xsaixy_mosaic(i,mosaic_i,j)=xsaixy(i,j)
            lai_mosaic(i,mosaic_i,j)=lai(i,j)


              DO soil_k=1,num_soil_layers

              TSLB_mosaic(i,num_soil_layers*(mosaic_i-1)+soil_k,j)=TSLB(i,soil_k,j)
              SMOIS_mosaic(i,num_soil_layers*(mosaic_i-1)+soil_k,j)=SMOIS(i,soil_k,j)
              SH2O_mosaic(i,num_soil_layers*(mosaic_i-1)+soil_k,j)=SH2O(i,soil_k,j)

              ENDDO

              IF (IOPT_HUE.eq.1) THEN
                DETENTION_STORAGEXY_mosaic(i,mosaic_i,J) = 0.
              ENDIF
           !!snow things
           isnowxy_mosaic(i,mosaic_i,j)=isnowxy(i,j)

           DO snow_layer=1,7
            zsnsoxy_mosaic(i,7*(mosaic_i-1)+snow_layer,j)=zsnsoxy(i,snow_layer-3,j)
           ENDDO

          DO top_layer=1,3
            tsnoxy_mosaic(i,3*(mosaic_i-1)+top_layer,j)=tsnoxy(i,top_layer-3,j)
            snicexy_mosaic(i,3*(mosaic_i-1)+top_layer,j)=snicexy(i,top_layer-3,j)
            snliqxy_mosaic(i,3*(mosaic_i-1)+top_layer,j)=snliqxy(i,top_layer-3,j)
           ENDDO

           ! Noah-MP irrigation scheme !pvk
                        if(iopt_irr >= 1 .and. iopt_irr <= 3) then
                           if(iopt_irrm == 0 .or. iopt_irrm ==1) then       ! sprinkler
                              IRNUMSI_mosaic(i,mosaic_i,j) = 0
                              IRWATSI_mosaic(i,mosaic_i,j) = 0.
                              IRELOSS_mosaic(i,mosaic_i,j) = 0.
                              IRRSPLH_mosaic(i,mosaic_i,j) = 0.
                           else if (iopt_irrm == 0 .or. iopt_irrm ==2) then ! micro or drip
                              IRNUMMI_mosaic(i,mosaic_i,j) = 0
                              IRWATMI_mosaic(i,mosaic_i,j) = 0.
                              IRMIVOL_mosaic(i,mosaic_i,j) = 0.
                           else if (iopt_irrm == 0 .or. iopt_irrm ==3) then ! flood
                              IRNUMFI_mosaic(i,mosaic_i,j) = 0
                              IRWATFI_mosaic(i,mosaic_i,j) = 0.
                              IRFIVOL_mosaic(i,mosaic_i,j) = 0.
                           end if
                        end if


           TR_URB2D_mosaic(i,mosaic_i,j)=TSK(i,j)
           TB_URB2D_mosaic(i,mosaic_i,j)=TSK(i,j)
           TG_URB2D_mosaic(i,mosaic_i,j)=TSK(i,j)
           TC_URB2D_mosaic(i,mosaic_i,j)=TSK(i,j)
           TS_URB2D_mosaic(i,mosaic_i,j)=TSK(i,j)
           QC_URB2D_mosaic(i,mosaic_i,j)=0.01
           SH_URB2D_mosaic(i,mosaic_i,j)=0
           LH_URB2D_mosaic(i,mosaic_i,j)=0
           G_URB2D_mosaic(i,mosaic_i,j)=0
           RN_URB2D_mosaic(i,mosaic_i,j)=0

          TRL_URB3D_mosaic(I,4*(mosaic_i-1)+1,J)=TSLB(I,1,J)+0.
          TRL_URB3D_mosaic(I,4*(mosaic_i-1)+2,J)=0.5*(TSLB(I,1,J)+TSLB(I,2,J))
          TRL_URB3D_mosaic(I,4*(mosaic_i-1)+3,J)=TSLB(I,2,J)+0.
          TRL_URB3D_mosaic(I,4*(mosaic_i-1)+4,J)=TSLB(I,2,J)+(TSLB(I,3,J)-TSLB(I,2,J))*0.29

          TBL_URB3D_mosaic(I,4*(mosaic_i-1)+1,J)=TSLB(I,1,J)+0.
          TBL_URB3D_mosaic(I,4*(mosaic_i-1)+2,J)=0.5*(TSLB(I,1,J)+TSLB(I,2,J))
          TBL_URB3D_mosaic(I,4*(mosaic_i-1)+3,J)=TSLB(I,2,J)+0.
          TBL_URB3D_mosaic(I,4*(mosaic_i-1)+4,J)=TSLB(I,2,J)+(TSLB(I,3,J)-TSLB(I,2,J))*0.29

          TGL_URB3D_mosaic(I,4*(mosaic_i-1)+1,J)=TSLB(I,1,J)
          TGL_URB3D_mosaic(I,4*(mosaic_i-1)+2,J)=TSLB(I,2,J)
          TGL_URB3D_mosaic(I,4*(mosaic_i-1)+3,J)=TSLB(I,3,J)
          TGL_URB3D_mosaic(I,4*(mosaic_i-1)+4,J)=TSLB(I,4,J)

          TGRL_URB3D_mosaic(I,4*(mosaic_i-1)+1,J)=TSLB(I,1,J)+0.
          TGRL_URB3D_mosaic(I,4*(mosaic_i-1)+2,J)=0.5*(TSLB(I,1,J)+TSLB(I,2,J))
          TGRL_URB3D_mosaic(I,4*(mosaic_i-1)+3,J)=TSLB(I,2,J)+0.
          TGRL_URB3D_mosaic(I,4*(mosaic_i-1)+4,J)=TSLB(I,2,J)+(TSLB(I,3,J)-TSLB(I,2,J))*0.29

          SMR_URB3D_mosaic(I,4*(mosaic_i-1)+1,J)=0.2
          SMR_URB3D_mosaic(I,4*(mosaic_i-1)+2,J)=0.2
          SMR_URB3D_mosaic(I,4*(mosaic_i-1)+3,J)=0.2
          SMR_URB3D_mosaic(I,4*(mosaic_i-1)+4,J)=0.2

          CMR_SFCDIF_mosaic(I,mosaic_i,J)=0.
          CHR_SFCDIF_mosaic(I,mosaic_i,J)=0.
          CMC_SFCDIF_mosaic(I,mosaic_i,J)=0.
          CHC_SFCDIF_mosaic(I,mosaic_i,J)=0.
          CMGR_SFCDIF_mosaic(I,mosaic_i,J)=0.
          CHGR_SFCDIF_mosaic(I,mosaic_i,J)=0.

          XXXR_URB2D_mosaic(I,mosaic_i,J)=0.
          XXXB_URB2D_mosaic(I,mosaic_i,J)=0.
          XXXG_URB2D_mosaic(I,mosaic_i,J)=0.
          XXXC_URB2D_mosaic(I,mosaic_i,J)=0.
          CMCR_URB2D_mosaic(I,mosaic_i,J)=0.

          TGR_URB2D_mosaic(I,mosaic_i,J)=TSK(I,J)+0.

          DRELR_URB2D_mosaic(I,mosaic_i,J)=0.
          DRELB_URB2D_mosaic(I,mosaic_i,J)=0.
          DRELG_URB2D_mosaic(I,mosaic_i,J)=0.

          FLXHUMR_URB2D_mosaic(I,mosaic_i,J)=0.
          FLXHUMB_URB2D_mosaic(I,mosaic_i,J)=0.
          FLXHUMG_URB2D_mosaic(I,mosaic_i,J)=0.


            ENDDO
          ENDDO
      ENDDO

   ! simple test

       DO i = its,itf
        DO j = jts,jtf

           IF ((xland(i,j).LT. 1.5 ) .AND. (mosaic_cat_index(I,1,J) .EQ. ISWATER)) THEN
             CALL wrf_message("After MOSAIC_INIT")
             WRITE (message_text,fmt='(a,2I6,2F8.2,2I6)') 'weird xland,xice,mosaic_cat_index and ivgtyp at I,J = ', &
                i,j,xland(i,j),xice(i,j),mosaic_cat_index(I,1,J),IVGTYP(i,j)
             CALL wrf_message(message_text)
           ENDIF

        ENDDO
      ENDDO

 ENDIF      !  for not restart


!------------------------------------------------------------------------------
! The following code adds the mosaic for loop. It functions by taking a for loop
! over different land types, and then averaging them all together!
! Added into this newest version of Noah-MP on 25 May 2022.
! Aaron A.
!------------------------------------------------------------------------------

!--------------------------------
  END SUBROUTINE NOAHMP_MOSAIC_INIT
!--------------------------------

!--------------------------------
SUBROUTINE NOAHMPLSM_MOSAIC_HUE(ITIMESTEP,        YR,   JULIAN,   COSZIN,XLAT,XLONG, & ! IN : Time/Space-related
                  DZ8W,       DT,       DZS,    NSOIL,       DX,            & ! IN : Model configuration
                IVGTYP,   ISLTYP,    VEGFRA,   VEGMAX,      TMN,            & ! IN : Vegetation/Soil characteristics
                 XLAND,     XICE,XICE_THRES,  CROPCAT,                      & ! IN : Vegetation/Soil characteristics
               PLANTING,  HARVEST,SEASON_GDD,                               &
                 IDVEG, IOPT_CRS,  IOPT_BTR, IOPT_RUN, IOPT_SFC, IOPT_FRZ,  & ! IN : User options
              IOPT_INF, IOPT_RAD,  IOPT_ALB, IOPT_SNF,IOPT_TBOT, IOPT_STC,  & ! IN : User options
              IOPT_GLA, IOPT_RSF, IOPT_SOIL,IOPT_PEDO,IOPT_CROP, IOPT_IRR,  & ! IN : User options
        IOPT_IRRM, IOPT_INFDV, IOPT_TDRN, IOPT_MOSAIC, IOPT_HUE, soilstep,  & ! IN: User options (aaron a.)
              IZ0TLND, SF_URBAN_PHYSICS,                                    & ! IN : User options
              SOILCOMP,  SOILCL1,  SOILCL2,   SOILCL3,  SOILCL4,            & ! IN : User options
                   T3D,     QV3D,     U_PHY,    V_PHY,    SWDOWN,     SWDDIR,&
                SWDDIF,      GLW,                                           & ! IN : Forcing
                 P8W3D,PRECIP_IN,        SR,                                & ! IN : Forcing
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
                 WOODXY, STBLCPXY, FASTCPXY,   XLAIXY,   XSAIXY,   TAUSSXY, & ! IN/OUT Noah MP only
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
                          LANDUSEF, landusef2, NLCAT,                     &      ! Added by Aaron A. **IMPORTANT FOR THE FOR LOOPS OF FRACTIONAL LAND USE
                          mosaic_cat_index, mosaic_cat,                           &      ! Added by Aaron A. **IMPORTANT FOR THE FOR LOOPS OF FRACTIONAL LAND USE
                          TSK_mosaic, HFX_mosaic, QFX_mosaic, LH_mosaic,             &      ! Added by Aaron A. IN/OUT LSM
                          GRDFLX_mosaic, SFCRUNOFF_mosaic, UDRUNOFF_mosaic,                     &      ! Added by Aaorn A. IN/OUT LSM
                          ALBEDO_mosaic, SNOWC_mosaic, TSLB_mosaic, SMOIS_mosaic,               &      ! Added by Aaron A. IN/OUT LSM
                          SH2O_mosaic,  CANWAT_mosaic, SNOW_mosaic, SNOWH_mosaic,               &      ! Added by Aaron A. IN/OUT LSM
                          ACSNOM_mosaic, ACSNOW_mosaic, EMISS_mosaic, QSFC_mosaic,              &      ! Added by Aaron A. IN/OUT LSM
                          Z0_mosaic, ZNT_mosaic,  RS_mosaic,                                    &      ! Added by Aaron A. IN/OUT LSM
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
                          ACC_ETRANIXY_mosaic,                                                  &
                          ACC_DWATERXY_mosaic, ACC_PRCPXY_mosaic,                               &
                          ACC_ECANXY_mosaic, ACC_EDIRXY_mosaic, ACC_ETRANXY_mosaic,             &
                          use_wudapt_lcz,                                                       & !Aaron A. Mosiac
                          TR_URB2D_mosaic,TB_URB2D_mosaic,                                      & !H urban  Aaron A. Mosaic
                          TG_URB2D_mosaic,TC_URB2D_mosaic,                                      & !H urban  Aaron A. Mosaic
                          QC_URB2D_mosaic,UC_URB2D_mosaic,                                      & !H urban  Aaron A. Mosaic
                          TRL_URB3D_mosaic,TBL_URB3D_mosaic,                                    & !H urban  Aaron A. Mosaic
                          TGL_URB3D_mosaic,                                                     & !H urban  Aaron A. Mosaic
                          SH_URB2D_mosaic,LH_URB2D_mosaic,                                      & !H urban  Aaron A. Mosaic
                          G_URB2D_mosaic,RN_URB2D_mosaic,                                       & !H urban  Aaron A. Mosaic
                          TS_URB2D_mosaic,CMR_SFCDIF_mosaic, CHR_SFCDIF_mosaic,                 & !H urban  Aaron A. Mosaic
                          CMC_SFCDIF_mosaic, CHC_SFCDIF_mosaic, CMGR_SFCDIF_mosaic,             &
                          CHGR_SFCDIF_mosaic, XXXR_URB2D_mosaic, XXXB_URB2D_mosaic, XXXC_URB2D_mosaic,            &
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
                          declin_urb,      omg_urb2d,                                           & !I urban
                          num_roof_layers,      num_wall_layers,      num_road_layers,          & !I urban
                          dzr,            dzb,            dzg,                                  & !I urban
                          cmcr_urb2d,      tgr_urb2d,     tgrl_urb3d,  smr_urb3d,               & !H urban
                          drelr_urb2d,    drelb_urb2d,    drelg_urb2d,                          & !H urban
                          flxhumr_urb2d,  flxhumb_urb2d,  flxhumg_urb2d,                        & !H urban
                          julday,             julyr,                                            &
                          frc_urb2d,    utype_urb2d,                                            & !I urban
                          chs,           chs2,           cqs2,                                  & !H
                          lb_urb2d,  hgt_urb2d,        lp_urb2d,                                &
                          mh_urb2d,     stdh_urb2d,       lf_urb2d,                             & !SLUCM
                          th_phy,            rho,          p_phy,        ust,                   & !I multi-layer urban
                          gmt,                                                                  & !I multi-layer urban MODIFIED BY AARON A.
                          RUNONSFXY,RUNONSFXY_mosaic,DETENTION_STORAGEXY,DETENTION_STORAGEXY_mosaic, &
                          VOL_FLUX_RUNONXY_mosaic, VOL_FLUX_SMXY_mosaic,                        &
                          MP_RAINC, MP_RAINNC, MP_SHCV, MP_SNOW, MP_GRAUP, MP_HAIL )! HUE SPECIFIC VARIABLES


!-----------------------------------------------------------------------
! This was added by Aaron A.
! The goal is to mirror the sub-heterogeneity scheme (commonly called MOSAIC)
! that has been implemented within the WRF framework. This is the following outline
! of the code:
! Declare all model variables
! Declare loop over each of the number of land catagories
! Move data from the mosaic variables to single variables
! Loop over each of then land catagories (and the urban models)
! Create spatially averaged variables (still in loop)
! Finish land use loop
! Send the averaged intermediate variable to the named variables to be output
! **Note that we are going to be passing along BOTH the mosaic variables and the named variables
!-----------------------------------------------------------------------
! Synthesize both the normal calls and the urban calls (USE functions)
!-----------------------------------------------------------------------
        USE MODULE_SF_NOAHMPLSM
    !    USE MODULE_SF_NOAHMPLSM, only: noahmp_options, NOAHMP_SFLX, noahmp_parameters
        USE module_sf_noahmp_glacier
        USE NOAHMP_TABLES

        USE module_sf_urban
        USE module_ra_gfdleta,  only: cal_mon_day
        USE module_model_constants, only: KARMAN, CP, XLV
!-----------------------------------------------------------------------
        IMPLICIT NONE
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!       NOAH-MP DECLIRATIONS
!-----------------------------------------------------------------------

INTEGER,                                         INTENT(IN   ) ::  ITIMESTEP ! timestep number
INTEGER,                                         INTENT(IN   ) ::  YR        ! 4-digit year
REAL,                                            INTENT(IN   ) ::  JULIAN    ! Julian day
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  COSZIN    ! cosine zenith angle
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  XLAT      ! latitude [rad]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  XLONG     ! latitude [rad]
REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  DZ8W      ! thickness of atmo layers [m]
REAL,                                            INTENT(IN   ) ::  DT        ! timestep [s]
REAL,    DIMENSION(1:nsoil),                     INTENT(IN   ) ::  DZS       ! thickness of soil layers [m]
INTEGER,                                         INTENT(IN   ) ::  NSOIL     ! number of soil layers
REAL,                                            INTENT(IN   ) ::  DX        ! horizontal grid spacing [m]
INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  IVGTYP    ! vegetation type
INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  ISLTYP    ! soil type
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  VEGFRA    ! vegetation fraction []
REAL,    DIMENSION( ims:ime ,         jms:jme ), INTENT(IN   ) ::  VEGMAX    ! annual max vegetation fraction []
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  TMN       ! deep soil temperature [K]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  XLAND     ! =2 ocean; =1 land/seaice
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  XICE      ! fraction of grid that is seaice
REAL,                                            INTENT(IN   ) ::  XICE_THRES! fraction of grid determining seaice
INTEGER,                                         INTENT(IN   ) ::  IDVEG     ! dynamic vegetation (1 -> off ; 2 -> on) with opt_crs = 1
INTEGER,                                         INTENT(IN   ) ::  IOPT_CRS  ! canopy stomatal resistance (1-> Ball-Berry; 2->Jarvis)
INTEGER,                                         INTENT(IN   ) ::  IOPT_BTR  ! soil moisture factor for stomatal resistance (1-> Noah; 2-> CLM; 3-> SSiB)
INTEGER,                                         INTENT(IN   ) ::  IOPT_RUN  ! runoff and groundwater (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS)
INTEGER,                                         INTENT(IN   ) ::  IOPT_SFC  ! surface layer drag coeff (CH & CM) (1->M-O; 2->Chen97)
INTEGER,                                         INTENT(IN   ) ::  IOPT_FRZ  ! supercooled liquid water (1-> NY06; 2->Koren99)
INTEGER,                                         INTENT(IN   ) ::  IOPT_INF  ! frozen soil permeability (1-> NY06; 2->Koren99)
INTEGER,                                         INTENT(IN   ) ::  IOPT_RAD  ! radiation transfer (1->gap=F(3D,cosz); 2->gap=0; 3->gap=1-Fveg)
INTEGER,                                         INTENT(IN   ) ::  IOPT_ALB  ! snow surface albedo (1->BATS; 2->CLASS)
INTEGER,                                         INTENT(IN   ) ::  IOPT_SNF  ! rainfall & snowfall (1-Jordan91; 2->BATS; 3->Noah)
INTEGER,                                         INTENT(IN   ) ::  IOPT_TBOT ! lower boundary of soil temperature (1->zero-flux; 2->Noah)
INTEGER,                                         INTENT(IN   ) ::  IOPT_STC  ! snow/soil temperature time scheme
INTEGER,                                         INTENT(IN   ) ::  IOPT_GLA  ! glacier option (1->phase change; 2->simple)
INTEGER,                                         INTENT(IN   ) ::  IOPT_RSF  ! surface resistance (1->Sakaguchi/Zeng; 2->Seller; 3->mod Sellers; 4->1+snow)
INTEGER,                                         INTENT(IN   ) ::  IOPT_SOIL ! soil configuration option
INTEGER,                                         INTENT(IN   ) ::  IOPT_PEDO ! soil pedotransfer function option
INTEGER,                                         INTENT(IN   ) ::  IOPT_CROP ! crop model option (0->none; 1->Liu et al.; 2->Gecros)
INTEGER,                                         INTENT(IN   ) ::  IOPT_IRR  ! irrigation scheme (0->none; >1 irrigation scheme ON)
INTEGER,                                         INTENT(IN   ) ::  IOPT_IRRM ! irrigation method
INTEGER,                                         INTENT(IN   ) ::  IOPT_INFDV! infiltration options for dynamic VIC infiltration (1->Philip; 2-> Green-Ampt;3->Smith-Parlange)
INTEGER,                                         INTENT(IN   ) ::  IOPT_TDRN ! tile drainage (0-> no tile drainage; 1-> simple tile drainage;2->Hooghoudt's)
REAL,                                            INTENT(IN   ) ::  soilstep ! soil timestep (s), default:0->same as main model timestep
INTEGER,                                         INTENT(IN   ) ::  IZ0TLND   ! option of Chen adjustment of Czil (not used)
INTEGER,                                         INTENT(IN   ) ::  sf_urban_physics   ! urban physics option
REAL,    DIMENSION( ims:ime,       8, jms:jme ), INTENT(IN   ) ::  SOILCOMP  ! soil sand and clay percentage
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SOILCL1   ! soil texture in layer 1
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SOILCL2   ! soil texture in layer 2
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SOILCL3   ! soil texture in layer 3
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SOILCL4   ! soil texture in layer 4
REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  T3D       ! 3D atmospheric temperature valid at mid-levels [K]
REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  QV3D      ! 3D water vapor mixing ratio [kg/kg_dry]
REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  U_PHY     ! 3D U wind component [m/s]
REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  V_PHY     ! 3D V wind component [m/s]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SWDOWN    ! solar down at surface [W m-2]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SWDDIF    ! solar down at surface [W m-2]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SWDDIR    ! solar down at surface [W m-2]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  GLW       ! longwave down at surface [W m-2]
REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  P8W3D     ! 3D pressure, valid at interface [Pa]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  PRECIP_IN ! total input precipitation [mm]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SR        ! frozen precipitation ratio [-]

!Optional Detailed Precipitation Partitioning Inputs
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ), OPTIONAL ::  MP_RAINC  ! convective precipitation entering land model [mm] ! MB/AN : v3.7
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ), OPTIONAL ::  MP_RAINNC ! large-scale precipitation entering land model [mm]! MB/AN : v3.7
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ), OPTIONAL ::  MP_SHCV   ! shallow conv precip entering land model [mm]      ! MB/AN : v3.7
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ), OPTIONAL ::  MP_SNOW   ! snow precipitation entering land model [mm]       ! MB/AN : v3.7
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ), OPTIONAL ::  MP_GRAUP  ! graupel precipitation entering land model [mm]    ! MB/AN : v3.7
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ), OPTIONAL ::  MP_HAIL   ! hail precipitation entering land model [mm]       ! MB/AN : v3.7

! Crop Model
INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  CROPCAT   ! crop catagory
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  PLANTING  ! planting date
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  HARVEST   ! harvest date
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SEASON_GDD! growing season GDD
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  GRAINXY   ! mass of grain XING [g/m2]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  GDDXY     ! growing degree days XING (based on 10C)
INTEGER,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  PGSXY

! gecros model
REAL,    DIMENSION( ims:ime,       60,jms:jme ), INTENT(INOUT) ::  gecros_state !  gecros crop

!Tile drain variables
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  QTDRAIN
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN)    ::  TD_FRACTION

#ifdef WRF_HYDRO
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  sfcheadrt,INFXSRT,soldrain,qtiledrain,ZWATBLE2D   ! for WRF-Hydro
#endif

! INOUT (with generic LSM equivalent)

REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  TSK       ! surface radiative temperature [K]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  HFX       ! sensible heat flux [W m-2]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  QFX       ! latent heat flux [kg s-1 m-2]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  LH        ! latent heat flux [W m-2]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  GRDFLX    ! ground/snow heat flux [W m-2]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SMSTAV    ! soil moisture avail. [not used]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SMSTOT    ! total soil water [mm][not used]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SFCRUNOFF ! accumulated surface runoff [m]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  UDRUNOFF  ! accumulated sub-surface runoff [m]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ALBEDO    ! total grid albedo []
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SNOWC     ! snow cover fraction []
REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(INOUT) ::  SMOIS     ! volumetric soil moisture [m3/m3]
REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(INOUT) ::  SH2O      ! volumetric liquid soil moisture [m3/m3]
REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(INOUT) ::  TSLB      ! soil temperature [K]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SNOW      ! snow water equivalent [mm]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SNOWH     ! physical snow depth [m]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  CANWAT    ! total canopy water + ice [mm]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACSNOM    ! accumulated snow melt (mm)
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACSNOW    ! accumulated snow on grid
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  EMISS     ! surface bulk emissivity
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  QSFC      ! bulk surface specific humidity
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  Z0        ! combined z0 sent to coupled model
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ZNT       ! combined z0 sent to coupled model
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  RS        ! Total stomatal resistance (s/m)

INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ISNOWXY   ! actual no. of snow layers
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  TVXY      ! vegetation leaf temperature
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  TGXY      ! bulk ground surface temperature
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  CANICEXY  ! canopy-intercepted ice (mm)
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  CANLIQXY  ! canopy-intercepted liquid water (mm)
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  EAHXY     ! canopy air vapor pressure (pa)
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  TAHXY     ! canopy air temperature (k)
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  CMXY      ! bulk momentum drag coefficient
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  CHXY      ! bulk sensible heat exchange coefficient
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  FWETXY    ! wetted or snowed fraction of the canopy (-)
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SNEQVOXY  ! snow mass at last time step(mm h2o)
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ALBOLDXY  ! snow albedo at last time step (-)
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  QSNOWXY   ! snowfall on the ground [mm/s]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  QRAINXY   ! rainfall on the ground [mm/s]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  WSLAKEXY  ! lake water storage [mm]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ZWTXY     ! water table depth [m]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  WAXY      ! water in the "aquifer" [mm]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  WTXY      ! groundwater storage [mm]
REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  TSNOXY    ! snow temperature [K]
REAL,    DIMENSION( ims:ime,-2:NSOIL, jms:jme ), INTENT(INOUT) ::  ZSNSOXY   ! snow layer depth [m]
REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  SNICEXY   ! snow layer ice [mm]
REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  SNLIQXY   ! snow layer liquid water [mm]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  LFMASSXY  ! leaf mass [g/m2]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  RTMASSXY  ! mass of fine roots [g/m2]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  STMASSXY  ! stem mass [g/m2]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  WOODXY    ! mass of wood (incl. woody roots) [g/m2]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  STBLCPXY  ! stable carbon in deep soil [g/m2]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  FASTCPXY  ! short-lived carbon, shallow soil [g/m2]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  XLAIXY    ! leaf area index
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  XSAIXY    ! stem area index
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  TAUSSXY   ! snow age factor
REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(INOUT) ::  SMOISEQ   ! eq volumetric soil moisture [m3/m3]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SMCWTDXY  ! soil moisture content in the layer to the water table when deep
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  DEEPRECHXY ! recharge to the water table when deep
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  RECHXY    ! recharge to the water table (diagnostic)

! OUT (with no Noah LSM equivalent)

REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  T2MVXY    ! 2m temperature of vegetation part
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  T2MBXY    ! 2m temperature of bare ground part
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  Q2MVXY    ! 2m mixing ratio of vegetation part
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  Q2MBXY    ! 2m mixing ratio of bare ground part
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  TRADXY    ! surface radiative temperature (k)
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  NEEXY     ! net ecosys exchange (g/m2/s CO2)
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  GPPXY     ! gross primary assimilation [g/m2/s C]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  NPPXY     ! net primary productivity [g/m2/s C]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FVEGXY    ! Noah-MP vegetation fraction [-]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  RUNSFXY   ! surface runoff [mm] per soil timestep
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  RUNSBXY   ! subsurface runoff [mm] per soil timestep
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  ECANXY    ! evaporation of intercepted water (mm/s)
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  EDIRXY    ! soil surface evaporation rate (mm/s]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  ETRANXY   ! transpiration rate (mm/s)
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FSAXY     ! total absorbed solar radiation (w/m2)
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FIRAXY    ! total net longwave rad (w/m2) [+ to atm]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  APARXY    ! photosyn active energy by canopy (w/m2)
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  PSNXY     ! total photosynthesis (umol co2/m2/s) [+]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SAVXY     ! solar rad absorbed by veg. (w/m2)
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SAGXY     ! solar rad absorbed by ground (w/m2)
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  RSSUNXY   ! sunlit leaf stomatal resistance (s/m)
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  RSSHAXY   ! shaded leaf stomatal resistance (s/m)
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  BGAPXY    ! between gap fraction
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  WGAPXY    ! within gap fraction
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  TGVXY     ! under canopy ground temperature[K]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  TGBXY     ! bare ground temperature [K]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  CHVXY     ! sensible heat exchange coefficient vegetated
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  CHBXY     ! sensible heat exchange coefficient bare-ground
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SHGXY     ! veg ground sen. heat [w/m2]   [+ to atm]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SHCXY     ! canopy sen. heat [w/m2]   [+ to atm]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SHBXY     ! bare sensible heat [w/m2]     [+ to atm]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  EVGXY     ! veg ground evap. heat [w/m2]  [+ to atm]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  EVBXY     ! bare soil evaporation [w/m2]  [+ to atm]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  GHVXY     ! veg ground heat flux [w/m2]  [+ to soil]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  GHBXY     ! bare ground heat flux [w/m2] [+ to soil]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  IRGXY     ! veg ground net LW rad. [w/m2] [+ to atm]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  IRCXY     ! canopy net LW rad. [w/m2] [+ to atm]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  IRBXY     ! bare net longwave rad. [w/m2] [+ to atm]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  TRXY      ! transpiration [w/m2]  [+ to atm]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  EVCXY     ! canopy evaporation heat [w/m2]  [+ to atm]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  CHLEAFXY  ! leaf exchange coefficient
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  CHUCXY    ! under canopy exchange coefficient
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  CHV2XY    ! veg 2m exchange coefficient
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  CHB2XY    ! bare 2m exchange coefficient
! additional output variables
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  PAHXY     ! precipitation advected heat
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  PAHGXY    ! precipitation advected heat
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  PAHBXY    ! precipitation advected heat
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  PAHVXY    ! precipitation advected heat
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QINTSXY
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QINTRXY
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QDRIPSXY
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QDRIPRXY
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QTHROSXY
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QTHRORXY
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QSNSUBXY
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QSNFROXY
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QSUBCXY
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QFROCXY
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QEVACXY
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QDEWCXY
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QFRZCXY
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QMELTCXY
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QSNBOTXY  !total liquid water (snowmelt + rain through pack)out of snowpack bottom [mm/s]
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QMELTXY   !snowmelt due to phase change (mm/s)
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  PONDINGXY
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FPICEXY    !fraction of ice in precip
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  RAINLSM     !rain rate                   (mm/s)  AJN
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SNOWLSM     !liquid equivalent snow rate (mm/s)  AJN
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FORCTLSM
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FORCQLSM
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FORCPLSM
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FORCZLSM
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FORCWLSM
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_SSOILXY
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_QINSURXY
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_QSEVAXY
REAL,    DIMENSION( ims:ime, 1:NSOIL, jms:jme ), INTENT(INOUT) ::  ACC_ETRANIXY
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  EFLXBXY
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SOILENERGY
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SNOWENERGY
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  CANHSXY
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_DWATERXY
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_PRCPXY
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_ECANXY
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_ETRANXY
REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_EDIRXY

INTEGER,  INTENT(IN   )   ::     ids,ide, jds,jde, kds,kde,  &  ! d -> domain
     &                           ims,ime, jms,jme, kms,kme,  &  ! m -> memory
     &                           its,ite, jts,jte, kts,kte      ! t -> tile

!2D inout irrigation variables
REAL,    DIMENSION( ims:ime,          jms:jme ),  INTENT(IN)    :: IRFRACT    ! irrigation fraction
REAL,    DIMENSION( ims:ime,          jms:jme ),  INTENT(IN)    :: SIFRACT    ! sprinkler irrigation fraction
REAL,    DIMENSION( ims:ime,          jms:jme ),  INTENT(IN)    :: MIFRACT    ! micro irrigation fraction
REAL,    DIMENSION( ims:ime,          jms:jme ),  INTENT(IN)    :: FIFRACT    ! flood irrigation fraction
INTEGER, DIMENSION( ims:ime,          jms:jme ),  INTENT(INOUT) :: IRNUMSI    ! irrigation event number, Sprinkler
INTEGER, DIMENSION( ims:ime,          jms:jme ),  INTENT(INOUT) :: IRNUMMI    ! irrigation event number, Micro
INTEGER, DIMENSION( ims:ime,          jms:jme ),  INTENT(INOUT) :: IRNUMFI    ! irrigation event number, Flood
REAL,    DIMENSION( ims:ime,          jms:jme ),  INTENT(INOUT) :: IRWATSI    ! irrigation water amount [m] to be applied, Sprinkler
REAL,    DIMENSION( ims:ime,          jms:jme ),  INTENT(INOUT) :: IRWATMI    ! irrigation water amount [m] to be applied, Micro
REAL,    DIMENSION( ims:ime,          jms:jme ),  INTENT(INOUT) :: IRWATFI    ! irrigation water amount [m] to be applied, Flood
REAL,    DIMENSION( ims:ime,          jms:jme ),  INTENT(INOUT) :: IRELOSS    ! loss of irrigation water to evaporation,sprinkler [m/timestep]
REAL,    DIMENSION( ims:ime,          jms:jme ),  INTENT(INOUT) :: IRSIVOL    ! amount of irrigation by sprinkler (mm)
REAL,    DIMENSION( ims:ime,          jms:jme ),  INTENT(INOUT) :: IRMIVOL    ! amount of irrigation by micro (mm)
REAL,    DIMENSION( ims:ime,          jms:jme ),  INTENT(INOUT) :: IRFIVOL    ! amount of irrigation by micro (mm)
REAL,    DIMENSION( ims:ime,          jms:jme ),  INTENT(INOUT) :: IRRSPLH    ! latent heating from sprinkler evaporation (w/m2)
CHARACTER(LEN=256),                               INTENT(IN)    :: LLANDUSE   ! landuse data name (USGS or MODIS_IGBP)

!ID local irrigation variables
REAL                                                            :: IRRFRA     ! irrigation fraction
REAL                                                            :: SIFAC      ! sprinkler irrigation fraction
REAL                                                            :: MIFAC      ! micro irrigation fraction
REAL                                                            :: FIFAC      ! flood irrigation fraction
INTEGER                                                         :: IRCNTSI    ! irrigation event number, Sprinkler
INTEGER                                                         :: IRCNTMI    ! irrigation event number, Micro
INTEGER                                                         :: IRCNTFI    ! irrigation event number, Flood
REAL                                                            :: IRAMTSI    ! irrigation water amount [m] to be applied, Sprinkler
REAL                                                            :: IRAMTMI    ! irrigation water amount [m] to be applied, Micro
REAL                                                            :: IRAMTFI    ! irrigation water amount [m] to be applied, Flood
REAL                                                            :: IREVPLOS   ! loss of irrigation water to evaporation,sprinkler [m/timestep]
REAL                                                            :: IRSIRATE   ! rate of irrigation by sprinkler [m/timestep]
REAL                                                            :: IRMIRATE   ! rate of irrigation by micro [m/timestep]
REAL                                                            :: IRFIRATE   ! rate of irrigation by micro [m/timestep]
REAL                                                            :: FIRR       ! latent heating due to sprinkler evaporation (w m-2)
REAL                                                            :: EIRR       ! evaporation due to sprinkler evaporation (mm/s)

! 1D equivalent of 2D/3D fields

! IN only

REAL                                :: COSZ         ! cosine zenith angle
REAL                                :: LAT          ! latitude [rad]
REAL                                :: Z_ML         ! model height [m]
INTEGER                             :: VEGTYP       ! vegetation type
INTEGER,    DIMENSION(NSOIL)        :: SOILTYP      ! soil type
INTEGER                             :: CROPTYPE     ! crop type
REAL                                :: FVEG         ! vegetation fraction [-]
REAL                                :: FVGMAX       ! annual max vegetation fraction []
REAL                                :: TBOT         ! deep soil temperature [K]
REAL                                :: T_ML         ! temperature valid at mid-levels [K]
REAL                                :: Q_ML         ! water vapor mixing ratio [kg/kg_dry]
REAL                                :: U_ML         ! U wind component [m/s]
REAL                                :: V_ML         ! V wind component [m/s]
REAL                                :: SWDN         ! solar down at surface [W m-2]
REAL                                :: LWDN         ! longwave down at surface [W m-2]
REAL                                :: P_ML         ! pressure, valid at interface [Pa]
REAL                                :: PSFC         ! surface pressure [Pa]
REAL                                :: PRCP         ! total precipitation entering  [mm/s]         ! MB/AN : v3.7
REAL                                :: PRCPCONV     ! convective precipitation entering  [mm/s]    ! MB/AN : v3.7
REAL                                :: PRCPNONC     ! non-convective precipitation entering [mm/s] ! MB/AN : v3.7
REAL                                :: PRCPSHCV     ! shallow convective precip entering  [mm/s]   ! MB/AN : v3.7
REAL                                :: PRCPSNOW     ! snow entering land model [mm/s]              ! MB/AN : v3.7
REAL                                :: PRCPGRPL     ! graupel entering land model [mm/s]           ! MB/AN : v3.7
REAL                                :: PRCPHAIL     ! hail entering land model [mm/s]              ! MB/AN : v3.7
REAL                                :: PRCPOTHR     ! other precip, e.g. fog [mm/s]                ! MB/AN : v3.7

! INOUT (with generic LSM equivalent)

REAL                                :: FSH          ! total sensible heat (w/m2) [+ to atm]
REAL                                :: SSOIL        ! soil heat heat (w/m2)
REAL                                :: SALB         ! surface albedo (-)
REAL                                :: FSNO         ! snow cover fraction (-)
REAL,   DIMENSION( 1:NSOIL)         :: SMCEQ        ! eq vol. soil moisture (m3/m3)
REAL,   DIMENSION( 1:NSOIL)         :: SMC          ! vol. soil moisture (m3/m3)
REAL,   DIMENSION( 1:NSOIL)         :: SMH2O        ! vol. soil liquid water (m3/m3)
REAL,   DIMENSION(-2:NSOIL)         :: STC          ! snow/soil tmperatures
REAL                                :: SWE          ! snow water equivalent (mm)
REAL                                :: SNDPTH       ! snow depth (m)
REAL                                :: EMISSI       ! net surface emissivity
REAL                                :: QSFC1D       ! bulk surface specific humidity

! INOUT (with no Noah LSM equivalent)

INTEGER                             :: ISNOW        ! actual no. of snow layers
REAL                                :: TV           ! vegetation canopy temperature
REAL                                :: TG           ! ground surface temperature
REAL                                :: CANICE       ! canopy-intercepted ice (mm)
REAL                                :: CANLIQ       ! canopy-intercepted liquid water (mm)
REAL                                :: EAH          ! canopy air vapor pressure (pa)
REAL                                :: TAH          ! canopy air temperature (k)
REAL                                :: CM           ! momentum drag coefficient
REAL                                :: CH           ! sensible heat exchange coefficient
REAL                                :: FWET         ! wetted or snowed fraction of the canopy (-)
REAL                                :: SNEQVO       ! snow mass at last time step(mm h2o)
REAL                                :: ALBOLD       ! snow albedo at last time step (-)
REAL                                :: QSNOW        ! snowfall on the ground [mm/s]
REAL                                :: QRAIN        ! rainfall on the ground [mm/s]
REAL                                :: WSLAKE       ! lake water storage [mm]
REAL                                :: ZWT          ! water table depth [m]
REAL                                :: WA           ! water in the "aquifer" [mm]
REAL                                :: WT           ! groundwater storage [mm]
REAL                                :: SMCWTD       ! soil moisture content in the layer to the water table when deep
REAL                                :: DEEPRECH     ! recharge to the water table when deep
REAL                                :: RECH         ! recharge to the water table (diagnostic)
REAL, DIMENSION(-2:NSOIL)           :: ZSNSO        ! snow layer depth [m]
REAL, DIMENSION(-2:              0) :: SNICE        ! snow layer ice [mm]
REAL, DIMENSION(-2:              0) :: SNLIQ        ! snow layer liquid water [mm]
REAL                                :: LFMASS       ! leaf mass [g/m2]
REAL                                :: RTMASS       ! mass of fine roots [g/m2]
REAL                                :: STMASS       ! stem mass [g/m2]
REAL                                :: WOOD         ! mass of wood (incl. woody roots) [g/m2]
REAL                                :: GRAIN        ! mass of grain XING [g/m2]
REAL                                :: GDD          ! mass of grain XING[g/m2]
INTEGER                             :: PGS          !stem respiration [g/m2/s]
REAL                                :: STBLCP       ! stable carbon in deep soil [g/m2]
REAL                                :: FASTCP       ! short-lived carbon, shallow soil [g/m2]
REAL                                :: PLAI         ! leaf area index
REAL                                :: PSAI         ! stem area index
REAL                                :: TAUSS        ! non-dimensional snow age

! tile drainage
REAL                                :: QTLDRN       ! tile drainage (mm)
REAL                                :: TDFRACMP     ! tile drainage map
#ifdef WRF_HYDRO
REAL                                :: WATBLED      ! water table depth for tile drainage
#endif

! OUT (with no Noah LSM equivalent)

REAL                                :: Z0WRF        ! combined z0 sent to coupled model
REAL                                :: T2MV         ! 2m temperature of vegetation part
REAL                                :: T2MB         ! 2m temperature of bare ground part
REAL                                :: Q2MV         ! 2m mixing ratio of vegetation part
REAL                                :: Q2MB         ! 2m mixing ratio of bare ground part
REAL                                :: TRAD         ! surface radiative temperature (k)
REAL                                :: NEE          ! net ecosys exchange (g/m2/s CO2)
REAL                                :: GPP          ! gross primary assimilation [g/m2/s C]
REAL                                :: NPP          ! net primary productivity [g/m2/s C]
REAL                                :: FVEGMP       ! greenness vegetation fraction [-]
REAL                                :: RUNSF        ! surface runoff [mm] per soil timestep
REAL                                :: RUNSB        ! subsurface runoff [mm] per soil timestep
REAL                                :: ECAN         ! evaporation of intercepted water (mm/s)
REAL                                :: ETRAN        ! transpiration rate (mm/s)
REAL                                :: ESOIL        ! soil surface evaporation rate (mm/s]
REAL                                :: FSA          ! total absorbed solar radiation (w/m2)
REAL                                :: FIRA         ! total net longwave rad (w/m2) [+ to atm]
REAL                                :: APAR         ! photosyn active energy by canopy (w/m2)
REAL                                :: PSN          ! total photosynthesis (umol co2/m2/s) [+]
REAL                                :: SAV          ! solar rad absorbed by veg. (w/m2)
REAL                                :: SAG          ! solar rad absorbed by ground (w/m2)
REAL                                :: RSSUN        ! sunlit leaf stomatal resistance (s/m)
REAL                                :: RSSHA        ! shaded leaf stomatal resistance (s/m)
REAL, DIMENSION(1:2)                :: ALBSND       ! snow albedo (direct)
REAL, DIMENSION(1:2)                :: ALBSNI       ! snow albedo (diffuse)
REAL                                :: RB           ! leaf boundary layer resistance (s/m)
REAL                                :: LAISUN       ! sunlit leaf area index (m2/m2)
REAL                                :: LAISHA       ! shaded leaf area index (m2/m2)
REAL                                :: BGAP         ! between gap fraction
REAL                                :: WGAP         ! within gap fraction
REAL                                :: TGV          ! under canopy ground temperature[K]
REAL                                :: TGB          ! bare ground temperature [K]
REAL                                :: CHV          ! sensible heat exchange coefficient vegetated
REAL                                :: CHB          ! sensible heat exchange coefficient bare-ground
REAL                                :: IRC          ! canopy net LW rad. [w/m2] [+ to atm]
REAL                                :: IRG          ! veg ground net LW rad. [w/m2] [+ to atm]
REAL                                :: SHC          ! canopy sen. heat [w/m2]   [+ to atm]
REAL                                :: SHG          ! veg ground sen. heat [w/m2]   [+ to atm]
REAL                                :: EVG          ! veg ground evap. heat [w/m2]  [+ to atm]
REAL                                :: GHV          ! veg ground heat flux [w/m2]  [+ to soil]
REAL                                :: IRB          ! bare net longwave rad. [w/m2] [+ to atm]
REAL                                :: SHB          ! bare sensible heat [w/m2]     [+ to atm]
REAL                                :: EVB          ! bare evaporation heat [w/m2]  [+ to atm]
REAL                                :: GHB          ! bare ground heat flux [w/m2] [+ to soil]
REAL                                :: TR           ! transpiration [w/m2]  [+ to atm]
REAL                                :: EVC          ! canopy evaporation heat [w/m2]  [+ to atm]
REAL                                :: CHLEAF       ! leaf exchange coefficient
REAL                                :: CHUC         ! under canopy exchange coefficient
REAL                                :: CHV2         ! veg 2m exchange coefficient
REAL                                :: CHB2         ! bare 2m exchange coefficient
REAL                                :: QINTS
REAL                                :: QINTR
REAL                                :: QDRIPS
REAL                                :: QDRIPR
REAL                                :: QTHROS
REAL                                :: QTHROR
REAL                                :: QSNSUB
REAL                                :: QSNFRO
REAL                                :: QEVAC
REAL                                :: QDEWC
REAL                                :: QSUBC
REAL                                :: QFROC
REAL                                :: QFRZC
REAL                                :: QMELTC
REAL                                :: PAHV    !precipitation advected heat - vegetation net (W/m2)
REAL                                :: PAHG    !precipitation advected heat - under canopy net (W/m2)
REAL                                :: PAHB    !precipitation advected heat - bare ground net (W/m2)
REAL                                :: PAH     !precipitation advected heat - total (W/m2)
REAL                                :: RAININ  !rain rate                   (mm/s)
REAL                                :: SNOWIN  !liquid equivalent snow rate (mm/s)
REAL                                :: ACC_SSOIL
REAL                                :: ACC_QINSUR
REAL                                :: ACC_QSEVA
REAL, DIMENSION( 1:NSOIL)           :: ACC_ETRANI       !transpiration rate (mm/s) [+]
REAL                                :: EFLXB
REAL                                :: XMF
REAL, DIMENSION( -2:NSOIL )         :: HCPCT
REAL                                :: DZSNSO
REAL                                :: CANHS   ! canopy heat storage change (w/m2)
REAL                                :: ACC_DWATER
REAL                                :: ACC_PRCP
REAL                                :: ACC_ECAN
REAL                                :: ACC_ETRAN
REAL                                :: ACC_EDIR

! Intermediate terms
REAL                                :: FPICE        ! snow fraction of precip
REAL                                :: FCEV         ! canopy evaporation heat (w/m2) [+ to atm]
REAL                                :: FGEV         ! ground evaporation heat (w/m2) [+ to atm]
REAL                                :: FCTR         ! transpiration heat flux (w/m2) [+ to atm]
REAL                                :: QSNBOT       ! total liquid water (snowmelt + rain through pack)out of snowpack bottom [mm/s]
REAL                                :: QMELT        ! snowmelt due to phase change (mm/s)
REAL                                :: PONDING      ! snowmelt with no pack [mm]
REAL                                :: PONDING1     ! snowmelt with no pack [mm]
REAL                                :: PONDING2     ! snowmelt with no pack [mm]

! Local terms

REAL, DIMENSION(1:60)               :: gecros1d     !  gecros crop
REAL                                :: gecros_dd ,gecros_tbem,gecros_emb ,gecros_ema, &
                                       gecros_ds1,gecros_ds2 ,gecros_ds1x,gecros_ds2x

REAL                                :: FSR          ! total reflected solar radiation (w/m2)
REAL, DIMENSION(-2:0)               :: FICEOLD      ! snow layer ice fraction []
REAL                                :: CO2PP        ! CO2 partial pressure [Pa]
REAL                                :: O2PP         ! O2 partial pressure [Pa]
REAL, DIMENSION(1:NSOIL)            :: ZSOIL        ! depth to soil interfaces [m]
REAL                                :: FOLN         ! nitrogen saturation [%]

REAL                                :: QC           ! cloud specific humidity for MYJ [not used]
REAL                                :: PBLH         ! PBL height for MYJ [not used]
REAL                                :: DZ8W1D       ! model level heights for MYJ [not used]

INTEGER                             :: I
INTEGER                             :: J
INTEGER                             :: K
INTEGER                             :: ICE
INTEGER                             :: SLOPETYP
LOGICAL                             :: IPRINT

INTEGER                             :: SOILCOLOR          ! soil color index
INTEGER                             :: IST          ! surface type 1-soil; 2-lake
INTEGER                             :: YEARLEN
REAL                                :: SOLAR_TIME
INTEGER                             :: JMONTH, JDAY

INTEGER, PARAMETER                  :: NSNOW = 3    ! number of snow layers fixed to 3
REAL, PARAMETER                     :: undefined_value = -1.E36

REAL, DIMENSION( 1:nsoil ) :: SAND
REAL, DIMENSION( 1:nsoil ) :: CLAY
REAL, DIMENSION( 1:nsoil ) :: ORGM

type(noahmp_parameters) :: parameters


!-----------------------------------------------------------------------
!   DECLARIATIONS - URBAN
!-----------------------------------------------------------------------

! input variables surface_driver --> lsm
     INTEGER,                                                INTENT(IN   ) :: use_wudapt_lcz
     INTEGER,                                                INTENT(IN   ) :: num_roof_layers
     INTEGER,                                                INTENT(IN   ) :: num_wall_layers
     INTEGER,                                                INTENT(IN   ) :: num_road_layers
     INTEGER,        DIMENSION( ims:ime, jms:jme ),          INTENT(INOUT   ) :: UTYPE_URB2D !aaron A. 
     REAL,           DIMENSION( ims:ime, jms:jme ),          INTENT(INOUT   ) :: FRC_URB2D ! aarohn A. 

     REAL, OPTIONAL, DIMENSION(1:num_roof_layers),           INTENT(IN   ) :: DZR
     REAL, OPTIONAL, DIMENSION(1:num_wall_layers),           INTENT(IN   ) :: DZB
     REAL, OPTIONAL, DIMENSION(1:num_road_layers),           INTENT(IN   ) :: DZG
     REAL, OPTIONAL,                                         INTENT(IN   ) :: DECLIN_URB
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),          INTENT(IN   ) :: OMG_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: TH_PHY
     REAL, OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: P_PHY
     REAL, OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: RHO

     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),          INTENT(INOUT) :: UST
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),          INTENT(INOUT) :: CHS, CHS2, CQS2

     INTEGER,  INTENT(IN   )   ::  julyr                  !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                     INTENT(IN   ) :: lp_urb2d
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                     INTENT(IN   ) :: lb_urb2d
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                     INTENT(IN   ) :: hgt_urb2d
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                     INTENT(IN   ) :: mh_urb2d
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                     INTENT(IN   ) :: stdh_urb2d
     REAL, OPTIONAL, DIMENSION( ims:ime, 4, jms:jme ),                  INTENT(IN   ) :: lf_urb2d
! local variables lsm --> urban

     INTEGER :: UTYPE_URB ! urban type [urban=1, suburban=2, rural=3]
     REAL    :: TA_URB       ! potential temp at 1st atmospheric level [K]
     REAL    :: QA_URB       ! mixing ratio at 1st atmospheric level  [kg/kg]
     REAL    :: UA_URB       ! wind speed at 1st atmospheric level    [m/s]
     REAL    :: U1_URB       ! u at 1st atmospheric level             [m/s]
     REAL    :: V1_URB       ! v at 1st atmospheric level             [m/s]
     REAL    :: SSG_URB      ! downward total short wave radiation    [W/m/m]
     REAL    :: LLG_URB      ! downward long wave radiation           [W/m/m]
     REAL    :: RAIN_URB     ! precipitation                          [mm/h]
     REAL    :: RHOO_URB     ! air density                            [kg/m^3]
     REAL    :: ZA_URB       ! first atmospheric level                [m]
     REAL    :: DELT_URB     ! time step                              [s]
     REAL    :: SSGD_URB     ! downward direct short wave radiation   [W/m/m]
     REAL    :: SSGQ_URB     ! downward diffuse short wave radiation  [W/m/m]
     REAL    :: XLAT_URB     ! latitude                               [deg]
     REAL    :: COSZ_URB     ! cosz
     REAL    :: OMG_URB      ! hour angle
     REAL    :: ZNT_URB      ! roughness length                       [m]
     REAL    :: TR_URB
     REAL    :: TB_URB
     REAL    :: TG_URB
     REAL    :: TC_URB
     REAL    :: QC_URB
     REAL    :: UC_URB
     REAL    :: XXXR_URB
     REAL    :: XXXB_URB
     REAL    :: XXXG_URB
     REAL    :: XXXC_URB
     REAL, DIMENSION(1:num_roof_layers) :: TRL_URB  ! roof layer temp [K]
     REAL, DIMENSION(1:num_wall_layers) :: TBL_URB  ! wall layer temp [K]
     REAL, DIMENSION(1:num_road_layers) :: TGL_URB  ! road layer temp [K]
     LOGICAL  :: LSOLAR_URB

!===hydrological variable for single layer UCM===

     REAL    :: DRELR_URB
     REAL    :: DRELB_URB
     REAL    :: DRELG_URB
     REAL    :: FLXHUMR_URB
     REAL    :: FLXHUMB_URB
     REAL    :: FLXHUMG_URB
     REAL    :: CMCR_URB
     REAL    :: TGR_URB

     REAL, DIMENSION(1:num_roof_layers) :: SMR_URB  ! green roof layer moisture
     REAL, DIMENSION(1:num_roof_layers) :: TGRL_URB ! green roof layer temp [K]

     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                    INTENT(INOUT) :: DRELR_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                    INTENT(INOUT) :: DRELB_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                    INTENT(INOUT) :: DRELG_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                    INTENT(INOUT) :: FLXHUMR_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                    INTENT(INOUT) :: FLXHUMB_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                    INTENT(INOUT) :: FLXHUMG_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                    INTENT(INOUT) :: CMCR_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                    INTENT(INOUT) :: TGR_URB2D

     REAL, OPTIONAL, DIMENSION( ims:ime, 1:num_roof_layers, jms:jme ), INTENT(INOUT) :: TGRL_URB3D
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:num_roof_layers, jms:jme ), INTENT(INOUT) :: SMR_URB3D


! state variable surface_driver <--> lsm <--> urban

     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: TR_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: TB_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: TG_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: TC_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: QC_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: UC_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: XXXR_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: XXXB_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: XXXG_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: XXXC_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: SH_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: LH_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: G_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: RN_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: TS_URB2D

     REAL, OPTIONAL, DIMENSION( ims:ime, 1:num_roof_layers, jms:jme ), INTENT(INOUT) :: TRL_URB3D
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:num_wall_layers, jms:jme ), INTENT(INOUT) :: TBL_URB3D
     REAL, OPTIONAL, DIMENSION( ims:ime, 1:num_road_layers, jms:jme ), INTENT(INOUT) :: TGL_URB3D

! output variable lsm --> surface_driver

     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT) :: PSIM_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT) :: PSIH_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT) :: GZ1OZ0_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT) :: U10_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT) :: V10_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT) :: TH2_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT) :: Q2_URB2D
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT) :: AKMS_URB2D
     REAL,           DIMENSION( ims:ime, jms:jme ), INTENT(OUT) :: UST_URB2D


! output variables urban --> lsm

     REAL :: TS_URB           ! surface radiative temperature    [K]
     REAL :: QS_URB           ! surface humidity                 [-]
     REAL :: SH_URB           ! sensible heat flux               [W/m/m]
     REAL :: LH_URB           ! latent heat flux                 [W/m/m]
     REAL :: LH_KINEMATIC_URB ! latent heat flux, kinetic  [kg/m/m/s]
     REAL :: SW_URB           ! upward short wave radiation flux [W/m/m]
     REAL :: ALB_URB          ! time-varying albedo            [fraction]
     REAL :: LW_URB           ! upward long wave radiation flux  [W/m/m]
     REAL :: G_URB            ! heat flux into the ground        [W/m/m]
     REAL :: RN_URB           ! net radiation                    [W/m/m]
     REAL :: PSIM_URB         ! shear f for momentum             [-]
     REAL :: PSIH_URB         ! shear f for heat                 [-]
     REAL :: GZ1OZ0_URB       ! shear f for heat                 [-]
     REAL :: U10_URB          ! wind u component at 10 m         [m/s]
     REAL :: V10_URB          ! wind v component at 10 m         [m/s]
     REAL :: TH2_URB          ! potential temperature at 2 m     [K]
     REAL :: Q2_URB           ! humidity at 2 m                  [-]
     REAL :: CHS_URB
     REAL :: CHS2_URB
     REAL :: UST_URB

! NUDAPT Parameters urban --> lam

     REAL :: mh_urb
     REAL :: stdh_urb
     REAL :: lp_urb
     REAL :: hgt_urb
     REAL, DIMENSION(4) :: lf_urb

! Local variables

     REAL :: Q1

! Noah UA changes

     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                     INTENT(INOUT) :: CMR_SFCDIF
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                     INTENT(INOUT) :: CHR_SFCDIF
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                     INTENT(INOUT) :: CMGR_SFCDIF
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                     INTENT(INOUT) :: CHGR_SFCDIF
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                     INTENT(INOUT) :: CMC_SFCDIF
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ),                     INTENT(INOUT) :: CHC_SFCDIF

! Variables for multi-layer UCM

     REAL, OPTIONAL,                                                    INTENT(IN   ) :: GMT
     INTEGER, OPTIONAL,                                                 INTENT(IN   ) :: JULDAY

! Local variables for multi-layer UCM

     REAL,    DIMENSION( its:ite, jts:jte) :: HFX_RURAL,GRDFLX_RURAL          ! ,LH_RURAL,RN_RURAL
     REAL,    DIMENSION( its:ite, jts:jte) :: QFX_RURAL                       ! ,QSFC_RURAL,UMOM_RURAL,VMOM_RURAL
     REAL,    DIMENSION( its:ite, jts:jte) :: ALB_RURAL,EMISS_RURAL,TSK_RURAL ! ,UST_RURAL
     REAL,    DIMENSION( its:ite, jts:jte) :: HFX_URB,UMOM_URB,VMOM_URB
     REAL,    DIMENSION( its:ite, jts:jte) :: QFX_URB
     REAL,    DIMENSION( its:ite, jts:jte) :: EMISS_URB
     REAL,    DIMENSION( its:ite, jts:jte) :: RL_UP_URB
     REAL,    DIMENSION( its:ite, jts:jte) :: RS_ABS_URB
     REAL,    DIMENSION( its:ite, jts:jte) :: GRDFLX_URB

     REAL :: SIGMA_SB,RL_UP_RURAL,RL_UP_TOT,RS_ABS_TOT,UMOM,VMOM
     REAL :: r1,r2,r3
     REAL :: CMR_URB, CHR_URB, CMC_URB, CHC_URB, CMGR_URB, CHGR_URB
     REAL :: frc_urb,lb_urb
     REAL :: check

!-----------------------------------------------------------------------
!NOAH-MP MOSAIC Related Variables added to decliration (Aaron A. based on Dan Li)
!-----------------------------------------------------------------------
     INTEGER,           INTENT(IN) :: NLCAT
     INTEGER,           INTENT(IN) :: IOPT_MOSAIC                                       !This tells us that the mosaic scheme is active
     INTEGER,            INTENT(IN) :: mosaic_cat                                                 !This tells how many mosaic catagories that are of interet,
     INTEGER,           INTENT(IN) :: IOPT_HUE                                            !This passed to help to see if we need to use the new physics
     REAL, DIMENSION( ims:ime, 1:NLCAT, jms:jme ) , OPTIONAL, INTENT(IN)::   LANDUSEF            !This is the original land-use fraction that was read in
     REAL, DIMENSION( ims:ime, 1:NLCAT, jms:jme ) , OPTIONAL, INTENT(INOUT)::   LANDUSEF2        !This is the land-use fraction that has been re-ordered

     INTEGER, DIMENSION( ims:ime, NLCAT, jms:jme ), OPTIONAL, INTENT(INOUT) :: mosaic_cat_index !This is the re-ordered mosaic catagory data
     REAL,    DIMENSION( ims:ime,             jms:jme ), OPTIONAL, INTENT(OUT  ) ::  RUNONSFXY ! Accumulated Surface runon [mm] Added by Aaron A. [
     REAL,    DIMENSION( ims:ime,1:mosaic_cat,jms:jme ), OPTIONAL, INTENT(OUT  ) ::  RUNONSFXY_mosaic ! Accumulated Surface runon [mm] Added by Aaron A.
     REAL :: RUNONSRF
     !variables with dimensions 1:mosaic_cat
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: TSK_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: HFX_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: QFX_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: LH_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: GRDFLX_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: SFCRUNOFF_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: UDRUNOFF_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: ALBEDO_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: SNOWC_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: CANWAT_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: SNOW_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: SNOWH_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: ACSNOM_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: ACSNOW_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: EMISS_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: QSFC_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: Z0_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: ZNT_mosaic

      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: tvxy_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: tgxy_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: canicexy_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: canliqxy_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: eahxy_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: tahxy_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: cmxy_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: chxy_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: fwetxy_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: sneqvoxy_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: alboldxy_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: qsnowxy_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: qrainxy_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: wslakexy_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: zwtxy_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: waxy_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: wtxy_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: lfmassxy_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: rtmassxy_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: stmassxy_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: woodxy_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: stblcpxy_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: fastcpxy_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: xsaixy_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: xlai_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: taussxy_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: smcwtdxy_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: deeprechxy_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: rechxy_mosaic

     !Variables that are NOAH MP output
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: t2mvxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: t2mbxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: q2mvxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: q2mbxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: tradxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: neexy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: gppxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: nppxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: fvegxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: runsfxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: runsbxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: ecanxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: edirxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: etranxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: fsaxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: firaxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: aparxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: psnxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: savxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: sagxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: rssunxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: rsshaxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: bgapxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: wgapxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: tgvxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: tgbxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: chvxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: chbxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: shgxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: shcxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: shbxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: evgxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: evbxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: ghvxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: ghbxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: irgxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: ircxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: irbxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: trxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: evcxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: chleafxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: chucxy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: chv2xy_mosaic
      !REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: chb2xy_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: rs_mosaic

      ! These are additional variables
     ! REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  PAHXY_mosaic     ! precipitation advected heat
      !REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  PAHGXY_mosaic    ! precipitation advected heat
      !REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  PAHBXY_mosaic    ! precipitation advected heat
      !REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  PAHVXY_mosaic    ! precipitation advected heat
     ! REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  QINTSXY_mosaic
     ! REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  QINTRXY_mosaic
     ! REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  QDRIPSXY_mosaic
     ! REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  QDRIPRXY_mosaic
     ! REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  QTHROSXY_mosaic
     ! REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  QTHRORXY_mosaic
     ! REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  QSNSUBXY_mosaic
     ! REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  QSNFROXY_mosaic
     ! REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  QSUBCXY_mosaic
     ! REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  QFROCXY_mosaic
     ! REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  QEVACXY_mosaic
     ! REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  QDEWCXY_mosaic
     ! REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  QFRZCXY_mosaic
     ! REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  QMELTCXY_mosaic
     ! REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  QSNBOTXY_mosaic  !total liquid water (snowmelt + rain through pack)out of snowpack bottom [mm/s]
     ! REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  QMELTXY_mosaic   !snowmelt due to phase change (mm/s)
     ! REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  PONDINGXY_mosaic
     ! REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  FPICEXY_mosaic    !fraction of ice in precip
      REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  ACC_SSOILXY_mosaic
      REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  ACC_QINSURXY_mosaic
      REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  ACC_QSEVAXY_mosaic
      REAL,    DIMENSION( ims:ime, 1:NSOIL*mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  ACC_ETRANIXY_mosaic
     ! REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  EFLXBXY_mosaic
     ! REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  SOILENERGY_mosaic
     ! REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  SNOWENERGY_mosaic
     ! REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  CANHSXY_mosaic
     REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  ACC_DWATERXY_mosaic
     REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  ACC_PRCPXY_mosaic
     REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  ACC_ECANXY_mosaic
     REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  ACC_ETRANXY_mosaic
     REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT) ::  ACC_EDIRXY_mosaic

      ! Irrigation variables that are needed
      !2D inout irrigation variables
      INTEGER, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ),  INTENT(INOUT) :: IRNUMSI_mosaic    ! irrigation event number, Sprinkler
      INTEGER, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ),  INTENT(INOUT) :: IRNUMMI_mosaic    ! irrigation event number, Micro
      INTEGER, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ),  INTENT(INOUT) :: IRNUMFI_mosaic    ! irrigation event number, Flood
      REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ),  INTENT(INOUT) :: IRWATSI_mosaic    ! irrigation water amount [m] to be applied, Sprinkler
      REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ),  INTENT(INOUT) :: IRWATMI_mosaic    ! irrigation water amount [m] to be applied, Micro
      REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ),  INTENT(INOUT) :: IRWATFI_mosaic    ! irrigation water amount [m] to be applied, Flood
      REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ),  INTENT(INOUT) :: IRELOSS_mosaic    ! loss of irrigation water to evaporation,sprinkler [m/timestep]
      REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ),  INTENT(INOUT) :: IRSIVOL_mosaic    ! amount of irrigation by sprinkler (mm)
      REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ),  INTENT(INOUT) :: IRMIVOL_mosaic    ! amount of irrigation by micro (mm)
      REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ),  INTENT(INOUT) :: IRFIVOL_mosaic    ! amount of irrigation by micro (mm)
      REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ),  INTENT(INOUT) :: IRRSPLH_mosaic    ! latent heating from sprinkler evaporation (w/m2)


     !snow variables, which have dimensions of 7*mosaic cat and 3 * number of mosaic cats
      INTEGER, DIMENSION(ims:ime, 1:mosaic_cat, jms:jme), OPTIONAL, INTENT(INOUT) :: isnowxy_mosaic     !actual no. of snow layers
      REAL, DIMENSION(ims:ime, 1:(nsoil+3)*mosaic_cat, jms:jme), OPTIONAL,  INTENT(INOUT) :: zsnsoxy_mosaic  !snow layer depth [m] **These have not been adjusted for indexing, Added in main driver
      REAL, DIMENSION(ims:ime, 1:3*mosaic_cat, jms:jme), OPTIONAL,  INTENT(INOUT) :: tsnoxy_mosaic   !snow temperature [K] **These have not been adjusted for indexing, Added in main driver
      REAL, DIMENSION(ims:ime, 1:3*mosaic_cat, jms:jme), OPTIONAL,  INTENT(INOUT) :: snicexy_mosaic  !snow layer ice [mm] **These have not been adjusted for indexing, Added in main driver
      REAL, DIMENSION(ims:ime, 1:3*mosaic_cat, jms:jme), OPTIONAL,  INTENT(INOUT) :: snliqxy_mosaic  !snow layer liquid water [mm] **These have not been adjusted for indexing, Added in main driver

     !variables that are for the soil layers
      REAL, DIMENSION( ims:ime, 1:nsoil*mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT):: TSLB_mosaic
      REAL, DIMENSION( ims:ime, 1:nsoil*mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT):: SMOIS_mosaic
      REAL, DIMENSION( ims:ime, 1:nsoil*mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT):: SH2O_mosaic
      REAL, DIMENSION( ims:ime, 1:nsoil*mosaic_cat, jms:jme ), OPTIONAL, INTENT(INOUT):: SMOISEQ_mosaic

    !variables for 2d urban model
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: TR_URB2D_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: TB_URB2D_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: TG_URB2D_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: TC_URB2D_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: QC_URB2D_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: SH_URB2D_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: LH_URB2D_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: G_URB2D_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: RN_URB2D_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: TS_URB2D_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: UC_URB2D_mosaic
      !values that are needed for the
      REAL, OPTIONAL, DIMENSION( ims:ime, 1:nsoil*mosaic_cat, jms:jme ), INTENT(INOUT) :: TRL_URB3D_mosaic
      REAL, OPTIONAL, DIMENSION( ims:ime, 1:nsoil*mosaic_cat, jms:jme ), INTENT(INOUT) :: TBL_URB3D_mosaic
      REAL, OPTIONAL, DIMENSION( ims:ime, 1:nsoil*mosaic_cat, jms:jme ), INTENT(INOUT) :: TGL_URB3D_mosaic
      REAL, OPTIONAL, DIMENSION( ims:ime, 1:nsoil*mosaic_cat, jms:jme ), INTENT(INOUT) :: TGRL_URB3D_mosaic
      REAL, OPTIONAL, DIMENSION( ims:ime, 1:nsoil*mosaic_cat, jms:jme ), INTENT(INOUT) :: SMR_URB3D_mosaic

      ! Additional to be able to restart
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: CMR_SFCDIF_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: CHR_SFCDIF_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: CMC_SFCDIF_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: CHC_SFCDIF_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: CMGR_SFCDIF_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: CHGR_SFCDIF_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: XXXR_URB2D_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: XXXB_URB2D_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: XXXG_URB2D_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: XXXC_URB2D_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: CMCR_URB2D_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: TGR_URB2D_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: DRELR_URB2D_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: DRELB_URB2D_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: DRELG_URB2D_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: FLXHUMR_URB2D_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: FLXHUMB_URB2D_mosaic
      REAL, DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ) , OPTIONAL, INTENT(INOUT):: FLXHUMG_URB2D_mosaic

      ! HUE noahmp values
      REAL,    DIMENSION( ims:ime, 1:mosaic_cat, jms:jme ), INTENT(INOUT), OPTIONAL :: DETENTION_STORAGEXY_mosaic  ! detention storage for green roof
      REAL,    DIMENSION(ims:ime, jms:jme ), INTENT(INOUT), OPTIONAL :: DETENTION_STORAGEXY

      REAL,   DIMENSION( ims:ime, 1:NSOIL*mosaic_cat, jms:jme),  INTENT(INOUT), OPTIONAL :: VOL_FLUX_SMXY_mosaic !fractional soil moisture volume pulled from sharing
      REAL,   DIMENSION( ims:ime, 1:mosaic_cat       ,jms:jme),  INTENT(INOUT), OPTIONAL :: VOL_FLUX_RUNONXY_mosaic

      ! THESE ARE INTERMEDIATES, THEY WILL BE USED TO KEEP TRACK OF ALL MOSAIC VALUES
      REAL, DIMENSION( ims:ime, jms:jme ) :: TSK_mosaic_avg, HFX_mosaic_avg, QFX_mosaic_avg, LH_mosaic_avg,  & ! These are in and out LSM options
            GRDFLX_mosaic_avg, SFCRUNOFF_mosaic_sum, UDRUNOFF_mosaic_sum, ALBEDO_mosaic_avg,                 &
            SNOWC_mosaic_avg, CANWAT_mosaic_avg, SNOW_mosaic_avg, SNOWH_mosaic_avg, ACSNOM_mosaic_avg,       &
            ACSNOW_mosaic_avg, EMISS_mosaic_avg, QSFC_mosaic_avg,Z0_mosaic_avg, ZNT_mosaic_avg,FAREA_mosaic_avg

      REAL, DIMENSION( ims:ime, jms:jme ) :: tvxy_mosaic_avg, tgxy_mosaic_avg, canicexy_mosaic_avg,             &
            canliqxy_mosaic_avg, eahxy_mosaic_avg, tahxy_mosaic_avg, cmxy_mosaic_avg, chxy_mosaic_avg,          &
            fwetxy_mosaic_avg, sneqvoxy_mosaic_avg, alboldxy_mosaic_avg, qsnowxy_mosaic_avg, qrainxy_mosaic_avg,&
            wslakexy_mosaic_avg, zwtxy_mosaic_avg, waxy_mosaic_avg, wtxy_mosaic_avg, lfmassxy_mosaic_avg,       &
            rtmassxy_mosaic_avg, stmassxy_mosaic_avg, woodxy_mosaic_avg,  &
            stblcpxy_mosaic_avg, fastcpxy_mosaic_avg, xsaixy_mosaic_avg, xlai_mosaic_avg,     &
            taussxy_mosaic_avg, rechxy_mosaic_avg, deeprechxy_mosaic_avg, smcwtdxy_mosaic_avg

      REAL, DIMENSION( ims:ime, jms:jme ) :: t2mvxy_mosaic_avg, t2mbxy_mosaic_avg, q2mvxy_mosaic_avg,           &
            q2mbxy_mosaic_avg, tradxy_mosaic_avg, neexy_mosaic_avg, gppxy_mosaic_avg, nppxy_mosaic_avg,         &
            fvegxy_mosaic_avg, runsfxy_mosaic_avg, runsbxy_mosaic_avg, ecanxy_mosaic_avg, edirxy_mosaic_avg,    &
            etranxy_mosaic_avg, fsaxy_mosaic_avg, firaxy_mosaic_avg, aparxy_mosaic_avg, psnxy_mosaic_avg,       &
            savxy_mosaic_avg, sagxy_mosaic_avg, rssunxy_mosaic_avg, rsshaxy_mosaic_avg, bgapxy_mosaic_avg,      &
            wgapxy_mosaic_avg, tgvxy_mosaic_avg, tgbxy_mosaic_avg, chvxy_mosaic_avg, chbxy_mosaic_avg,          &
            shgxy_mosaic_avg, shcxy_mosaic_avg, shbxy_mosaic_avg, evgxy_mosaic_avg, evbxy_mosaic_avg,           &
            ghvxy_mosaic_avg, ghbxy_mosaic_avg,                                                                 &
            irgxy_mosaic_avg, ircxy_mosaic_avg, irbxy_mosaic_avg, trxy_mosaic_avg, evcxy_mosaic_avg,            &
            chleafxy_mosaic_avg, chucxy_mosaic_avg, chv2xy_mosaic_avg, chb2xy_mosaic_avg, rs_mosaic_avg

      ! irrigation intermediate variables
      REAL, DIMENSION( ims:ime, jms:jme) :: IRWATSI_mosaic_avg, IRWATMI_mosaic_avg, IRWATFI_mosaic_avg,         &
            IRELOSS_mosaic_avg, IRSIVOL_mosaic_avg, IRMIVOL_mosaic_avg,                     &
            IRFIVOL_mosaic_avg, IRRSPLH_mosaic_avg

      INTEGER, DIMENSION( ims:ime, jms:jme) :: IRNUMSI_mosaic_avg, IRNUMMI_mosaic_avg, IRNUMFI_mosaic_avg

      ! Extra variables averaged needed for extra outputs
      REAL, DIMENSION( ims:ime, jms:jme) :: QINTSXY_mosaic_avg, QINTRXY_mosaic_avg, QDRIPSXY_mosaic_avg,        &
            QDRIPRXY_mosaic_avg, QTHROSXY_mosaic_avg, QTHRORXY_mosaic_avg, QSNSUBXY_mosaic_avg,                 &
            QSNFROXY_mosaic_avg, QSUBCXY_mosaic_avg, QFROCXY_mosaic_avg, QEVACXY_mosaic_avg,                    &
            QDEWCXY_mosaic_avg, QFRZCXY_mosaic_avg, QMELTCXY_mosaic_avg, QSNBOTXY_mosaic_avg,                   &
            QMELTXY_mosaic_avg, PONDINGXY_mosaic_avg, PAHXY_mosaic_avg, PAHVXY_mosaic_avg,                      &
            PAHBXY_mosaic_avg, PAHGXY_mosaic_avg, FPICEXY_mosaic_avg

      ! Soil and ACC
      REAL, DIMENSION( ims:ime, jms:jme) :: ACC_SSOILXY_mosaic_avg, ACC_QINSURXY_mosaic_avg, ACC_QSEVAXY_mosaic_avg,  &
            EFLXBXY_mosaic_avg, SOILENERGY_mosaic_avg, SNOWENERGY_mosaic_avg, CANHSXY_mosaic_avg,                     &
            ACC_DWATERXY_mosaic_avg, ACC_PRCPXY_mosaic_avg, ACC_ECANXY_mosaic_avg, ACC_ETRANXY_mosaic_avg, ACC_EDIRXY_mosaic_avg

      REAL, DIMENSION (ims:ime, jms:jme) :: RUNONSFXY_mosaic_avg

      !REAL, DIMENSION (ims:ime, 1:(nsoil+3)*mosaic_cat, jms:jme) DZSNO_mosaic !thickness of the different layers for soil and snow

      INTEGER, DIMENSION( ims:ime, jms:jme ) :: isnowxy_mosaic_avg

      REAL, DIMENSION( ims:ime, 1:NSOIL, jms:jme )                     ::  TSLB_mosaic_avg,SMOIS_mosaic_avg,SH2O_mosaic_avg,SMOISEQ_mosaic_avg, ACC_ETRANIXY_mosaic_avg

      REAL, DIMENSION( ims:ime, 1:3, jms:jme) :: tsnoxy_mosaic_avg, snicexy_mosaic_avg, snliqxy_mosaic_avg

      REAL, DIMENSION( ims:ime, 1:7, jms:jme) :: zsnsoxy_mosaic_avg

     !HUE noah-mp variables
     INTEGER, DIMENSION( ims:ime, jms:jme ) ::    IVGTYP_dominant
     INTEGER ::  mosaic_i, URBAN_METHOD, zo_avg_option,LAYER
     REAL :: FAREA,FAREA2
     LOGICAL :: IPRINT_mosaic !may not be needed Aaron A.
     INTEGER, PARAMETER                   :: NSOIL_GR = 1

     REAL,    DIMENSION( 1:nsoil)  ::  SMC_intermediate     ! volumetric soil moisture [m3/m3]
     REAL,    DIMENSION( 1:nsoil)  ::  SH2O_intermediate      ! volumetric liquid soil moisture [m3/m3]
     REAL,    DIMENSION( 1:nsoil)  ::  BTRANI_dummy ! BTRANI_output
     REAL                          ::  DETENTION_STORAGE
     ! create intermediates for lateral transfers
     REAL                          :: FAREAXY
     REAL,   DIMENSION( 1:nsoil)   :: VOL_FLUX_SM
     REAL                          :: VOL_FLUX_RUNON
     INTEGER                       :: TRIGGER = 0
     !------------------------------------------------------------------
     !Begin code calls
     !------------------------------------------------------------------
      IPRINT_mosaic = .false.

         CALL NOAHMP_OPTIONS(IDVEG  ,IOPT_CRS  ,IOPT_BTR  ,IOPT_RUN  ,IOPT_SFC  ,IOPT_FRZ , &
                     IOPT_INF  ,IOPT_RAD  ,IOPT_ALB  ,IOPT_SNF  ,IOPT_TBOT, IOPT_STC  ,     &
     		IOPT_RSF  ,IOPT_SOIL ,IOPT_PEDO ,IOPT_CROP ,IOPT_IRR , IOPT_IRRM ,     &
                     IOPT_INFDV,IOPT_TDRN, IOPT_MOSAIC, IOPT_HUE )

         IPRINT    =  .false.                     ! debug printout

     ! for using soil update timestep difference from noahmp main timestep
         calculate_soil = .false.
         soil_update_steps = nint(soilstep/DT)  ! 3600 = 1 hour
         soil_update_steps = max(soil_update_steps,1)
         if ( soil_update_steps == 1 ) then
           ACC_SSOILXY  = 0.0
           ACC_QINSURXY = 0.0
           ACC_QSEVAXY  = 0.0
           ACC_ETRANIXY = 0.0
           ACC_DWATERXY = 0.0
           ACC_PRCPXY   = 0.0
           ACC_ECANXY   = 0.0
           ACC_ETRANXY  = 0.0
           ACC_EDIRXY   = 0.0
         end if
         if ( soil_update_steps > 1 ) then
          if ( mod(itimestep,soil_update_steps) == 1 ) then
           ACC_SSOILXY  = 0.0
           ACC_QINSURXY = 0.0
           ACC_QSEVAXY  = 0.0
           ACC_ETRANIXY = 0.0
           ACC_DWATERXY = 0.0
           ACC_PRCPXY   = 0.0
           ACC_ECANXY   = 0.0
           ACC_ETRANXY  = 0.0
           ACC_EDIRXY   = 0.0
          end if
         end if

         if (mod(itimestep,soil_update_steps) == 0) calculate_soil = .true.
     ! end soil timestep

         YEARLEN = 365                            ! find length of year for phenology (also S Hemisphere)
         if (mod(YR,4) == 0) then
            YEARLEN = 366
            if (mod(YR,100) == 0) then
               YEARLEN = 365
               if (mod(YR,400) == 0) then
                  YEARLEN = 366
               endif
            endif
         endif

         ZSOIL(1) = -DZS(1)                    ! depth to soil interfaces (<0) [m]
         DO K = 2, NSOIL
            ZSOIL(K) = -DZS(K) + ZSOIL(K-1)
         END DO

         JLOOP : DO J=jts,jte ! Begin DO loop for J Aaron A.
             IF(ITIMESTEP == 1)THEN
               DO I=its,ite
                  IF((XLAND(I,J)-1.5) >= 0.) THEN    ! Open water case
                     IF(XICE(I,J) == 1. .AND. IPRINT) PRINT *,' sea-ice at water point, I=',I,'J=',J
                     SMSTAV(I,J) = 1.0
                     SMSTOT(I,J) = 1.0
                     DO K = 1, NSOIL
                        SMOIS(I,K,J) = 1.0
                         TSLB(I,K,J) = 273.16
                     ENDDO
                  ELSE
                     IF(XICE(I,J) == 1.) THEN        ! Sea-ice case
                        SMSTAV(I,J) = 1.0
                        SMSTOT(I,J) = 1.0
                        DO K = 1, NSOIL
                           SMOIS(I,K,J) = 1.0
                        ENDDO
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF                                                               ! end of initialization over ocean

    !-----------------------------------------------------------------------
       ILOOP : DO I = its, ite

        IF (XICE(I,J) >= XICE_THRES) THEN !this is glacier points (NEED TO ADD END IF Statement )
           ICE = 1                            ! Sea-ice point

           SH2O  (i,1:NSOIL,j) = 1.0
           XLAIXY(i,j)         = 0.01

           CYCLE ILOOP ! Skip any processing at sea-ice points

        ELSE

          IF((XLAND(I,J)-1.5) >= 0.) CYCLE ILOOP   ! Open water case


               IVGTYP_dominant(I,J)=IVGTYP(I,J)                            ! SAVE THE DOMINANT Vegitation


               ! INITIALIZE THE AREA-AVERAGED variables (all average variables thtat I added)
               TSK_mosaic_avg(i,j) = 0.0
               HFX_mosaic_avg(i,j) = 0.0
               QFX_mosaic_avg(i,j) = 0.0
               LH_mosaic_avg(i,j)  = 0.0
               GRDFLX_mosaic_avg(i,j) = 0.0
               SFCRUNOFF_mosaic_sum(i,j) = 0.0
               UDRUNOFF_mosaic_sum(i,j) = 0.0
               ALBEDO_mosaic_avg(i,j) = 0.0
               SNOWC_mosaic_avg(i,j) = 0.0
               CANWAT_mosaic_avg(i,j) = 0.0
               SNOW_mosaic_avg(i,j) = 0.0
               ACSNOM_mosaic_avg(i,j) = 0.0
               ACSNOW_mosaic_avg(i,j) = 0.0
               EMISS_mosaic_avg(i,j) = 0.0
               QSFC_mosaic_avg(i,j) = 0.0
               Z0_mosaic_avg(i,j) = 0.0
               ZNT_mosaic_avg(i,j) = 0.0

               FAREA_mosaic_avg(i,j) = 0.0
               isnowxy_mosaic_avg(i,j) = 0
               tvxy_mosaic_avg(i,j) = 0.0
               tgxy_mosaic_avg(i,j) = 0.0
               canicexy_mosaic_avg(i,j) = 0.0
               canliqxy_mosaic_avg(i,j) = 0.0
               eahxy_mosaic_avg(i,j) = 0.0
               tahxy_mosaic_avg(i,j) = 0.0
               cmxy_mosaic_avg(i,j) = 0.0
               chxy_mosaic_avg(i,j) = 0.0
               fwetxy_mosaic_avg(i,j) = 0.0
               sneqvoxy_mosaic_avg(i,j) = 0.0
               alboldxy_mosaic_avg(i,j) = 0.0
               qsnowxy_mosaic_avg(i,j) = 0.0
               qrainxy_mosaic_avg(i,j) = 0.0
               wslakexy_mosaic_avg(i,j) = 0.0
               zwtxy_mosaic_avg(i,j) = 0.0
               waxy_mosaic_avg(i,j) = 0.0
               wtxy_mosaic_avg(i,j) = 0.0
               lfmassxy_mosaic_avg(i,j) = 0.0
               rtmassxy_mosaic_avg(i,j) = 0.0
               stmassxy_mosaic_avg(i,j) = 0.0
               woodxy_mosaic_avg(i,j) = 0.0
               stblcpxy_mosaic_avg(i,j) = 0.0
               fastcpxy_mosaic_avg(i,j) = 0.0
               xsaixy_mosaic_avg(i,j) = 0.0
               xlai_mosaic_avg(i,j) = 0.0
               taussxy_mosaic_avg(i,j) = 0.0
               smcwtdxy_mosaic_avg(i,j) = 0.0
               deeprechxy_mosaic_avg(i,j) = 0.0
               rechxy_mosaic_avg(i,j) = 0.0

               t2mvxy_mosaic_avg(i,j) = 0.0
               t2mbxy_mosaic_avg(i,j) = 0.0
               q2mvxy_mosaic_avg(i,j) = 0.0
               q2mbxy_mosaic_avg(i,j) = 0.0
               tradxy_mosaic_avg(i,j) = 0.0
               neexy_mosaic_avg(i,j) = 0.0
               gppxy_mosaic_avg(i,j) = 0.0
               nppxy_mosaic_avg(i,j) = 0.0
               fvegxy_mosaic_avg(i,j) = 0.0
               runsfxy_mosaic_avg(i,j) = 0.0
               runsbxy_mosaic_avg(i,j) = 0.0
               ecanxy_mosaic_avg(i,j) = 0.0
               edirxy_mosaic_avg(i,j) = 0.0
               etranxy_mosaic_avg(i,j) = 0.0
               fsaxy_mosaic_avg(i,j) = 0.0
               firaxy_mosaic_avg(i,j) = 0.0
               aparxy_mosaic_avg(i,j) = 0.0
               psnxy_mosaic_avg(i,j) = 0.0
               savxy_mosaic_avg(i,j) = 0.0
               sagxy_mosaic_avg(i,j) = 0.0
               rssunxy_mosaic_avg(i,j) = 0.0
               rsshaxy_mosaic_avg(i,j) = 0.0
               bgapxy_mosaic_avg(i,j) = 0.0
               wgapxy_mosaic_avg(i,j) = 0.0
               tgvxy_mosaic_avg(i,j) = 0.0
               tgbxy_mosaic_avg(i,j) = 0.0
               chvxy_mosaic_avg(i,j) = 0.0
               chbxy_mosaic_avg(i,j) = 0.0
               shgxy_mosaic_avg(i,j) = 0.0
               shcxy_mosaic_avg(i,j) = 0.0
               shbxy_mosaic_avg(i,j) = 0.0
               evgxy_mosaic_avg(i,j) = 0.0
               evbxy_mosaic_avg(i,j) = 0.0
               ghvxy_mosaic_avg(i,j) = 0.0
               ghbxy_mosaic_avg(i,j) = 0.0
               irgxy_mosaic_avg(i,j) = 0.0
               ircxy_mosaic_avg(i,j) = 0.0
               irbxy_mosaic_avg(i,j) = 0.0
               trxy_mosaic_avg(i,j) = 0.0
               evcxy_mosaic_avg(i,j) = 0.0
               chleafxy_mosaic_avg(i,j) = 0.0
               chucxy_mosaic_avg(i,j) = 0.0
               chv2xy_mosaic_avg(i,j) = 0.0
               chb2xy_mosaic_avg(i,j) = 0.0
               rs_mosaic_avg(i,j) = 0.0

               !Loops over the SOIL Layers
               DO LAYER=1,NSOIL

                   TSLB_mosaic_avg(i,LAYER,j) = 0.0
                   SMOIS_mosaic_avg(i,LAYER,j) = 0.0
                   SH2O_mosaic_avg(i,LAYER,j) = 0.0
                   smoiseq_mosaic_avg(i,LAYER,j) = 0.0
                   ACC_ETRANIXY_mosaic_avg(i,LAYER,j) = 0.0

               ENDDO

               DO LAYER=1,3

                   tsnoxy_mosaic_avg(i,LAYER,j) = 0.0
                   snicexy_mosaic_avg(i,LAYER,j) = 0.0
                   snliqxy_mosaic_avg(i,LAYER,j) = 0.0

               ENDDO

               DO LAYER=1,7

                   zsnsoxy_mosaic_avg(i,LAYER,j) = 0.0

               ENDDO

               !! additional varaibles that are going to be averaged

               QINTSXY_mosaic_avg(i,j) = 0.0
               QINTRXY_mosaic_avg(i,j) = 0.0
               QDRIPSXY_mosaic_avg(i,j) = 0.0
               QDRIPRXY_mosaic_avg(i,j) = 0.0
               QTHROSXY_mosaic_avg(i,j) = 0.0
               QTHRORXY_mosaic_avg(i,j) = 0.0
               QSNSUBXY_mosaic_avg(i,j) = 0.0
               QSNFROXY_mosaic_avg(i,j) = 0.0
               QSUBCXY_mosaic_avg(i,j) = 0.0
               QFROCXY_mosaic_avg(i,j) = 0.0
               QEVACXY_mosaic_avg(i,j) = 0.0
               QDEWCXY_mosaic_avg(i,j) = 0.0
               QFRZCXY_mosaic_avg(i,j) = 0.0
               QMELTCXY_mosaic_avg(i,j) = 0.0
               QSNBOTXY_mosaic_avg(i,j) = 0.0
               QMELTXY_mosaic_avg(i,j) = 0.0
               PONDINGXY_mosaic_avg(i,j) = 0.0
               PAHXY_mosaic_avg(i,j) = 0.0
               PAHVXY_mosaic_avg(i,j) = 0.0
               PAHBXY_mosaic_avg(i,j) = 0.0
               FPICEXY_mosaic_avg(i,j) = 0.0

               ! Soil and ACC averaged variables
               ACC_SSOILXY_mosaic_avg(i,j) = 0.0
               ACC_QINSURXY_mosaic_avg(i,j) = 0.0
               ACC_QSEVAXY_mosaic_avg(i,j) = 0.0
               EFLXBXY_mosaic_avg(i,j) = 0.0
               SOILENERGY_mosaic_avg(i,j) = 0.0
               SNOWENERGY_mosaic_avg(i,j) = 0.0
               CANHSXY_mosaic_avg(i,j) = 0.0
               ACC_DWATERXY_mosaic_avg(i,j) = 0.0
               ACC_PRCPXY_mosaic_avg(i,j) = 0.0
               ACC_ECANXY_mosaic_avg(i,j) = 0.0

               !! THESE are the irrigation values
               ! irrigation intermediate variables
               IRWATSI_mosaic_avg(i,j) = 0.0
               IRWATMI_mosaic_avg(i,j) = 0.0
               IRWATFI_mosaic_avg(i,j) = 0.0
               IRELOSS_mosaic_avg(i,j) = 0.0
               IRSIVOL_mosaic_avg(i,j) = 0.0
               IRSIVOL_mosaic_avg(i,j) = 0.0
               IRMIVOL_mosaic_avg(i,j) = 0.0
               IRFIVOL_mosaic_avg(i,j) = 0.0
               IRRSPLH_mosaic_avg(i,j) = 0.0

               IRNUMSI_mosaic_avg(i,j) = 0
               IRNUMMI_mosaic_avg(i,j) = 0
               IRNUMFI_mosaic_avg(i,j) = 0


               ! These are needed for the

               MOSAIC_LOOP : DO mosaic_i = mosaic_cat, 1, -1 !looping from the nth mosaic catagory to the largest mosaic catagory

               ! We need to move from our 3D data to the 2D data sets
               !These are the IN/Out values that have an LSM equivelant
               IVGTYP(I,J) = mosaic_cat_index(I,mosaic_i,J)
               TSK(I,J) = TSK_mosaic(I,mosaic_i,J)
               HFX(I,J) = HFX_mosaic(I,mosaic_i,J)
               QFX(I,J) = QFX_mosaic(I,mosaic_i,J)
               LH(I,J) = LH_mosaic(I,mosaic_i,J)
               GRDFLX(I,J) = GRDFLX_mosaic(I,mosaic_i,J) !do not need SMSTAV (not used)
               SFCRUNOFF(I,J) = SFCRUNOFF_mosaic(I,mosaic_i,J) !do not need SMSTOT (not used)
               UDRUNOFF(I,J) = UDRUNOFF_mosaic(I,mosaic_i,J)
               ALBEDO(I,J) = ALBEDO_mosaic(I,mosaic_i,J)
               SNOWC(I,J) = SNOWC_mosaic(I,mosaic_i,J)

               !These are the soil variables
               DO LAYER=1,NSOIL

                   SMOIS(I,LAYER,J) = SMOIS_mosaic(I,NSOIL*(mosaic_i-1)+LAYER,J)
                   SH2O(I,LAYER,J) = SH2O_mosaic(I,NSOIL*(mosaic_i-1)+LAYER,J)
                   TSLB(I,LAYER,J) = TSLB_mosaic(I,NSOIL*(mosaic_i-1)+LAYER,J)
                   SMOISEQ(I,LAYER,J) = SMOISEQ_mosaic(I,NSOIL*(mosaic_i-1)+LAYER,J)
                   ACC_ETRANIXY(I,LAYER,J) = ACC_ETRANIXY_mosaic(I,NSOIL*(mosaic_i-1)+LAYER,J)

               ENDDO

               SNOW(I,J) = SNOW_mosaic(I,mosaic_i,J)
               SNOWH(I,J) = SNOWH_mosaic(I,mosaic_i,J)
               CANWAT(I,J) = CANWAT_mosaic(I,mosaic_i,J)
               ACSNOM(I,J) = ACSNOM_mosaic(I,mosaic_i,J)
               ACSNOW(I,J) = ACSNOW_mosaic(I,mosaic_i,J)
               EMISS(I,J) = EMISS_mosaic(I,mosaic_i,J)
               QSFC(I,J) = QSFC_mosaic(I,mosaic_i,J)
               ZNT(I,J) = ZNT_mosaic(I,mosaic_i,J)
               Z0(I,J) = Z0_mosaic(I,mosaic_i,J)
               RS(I,J) = rs_mosaic(I,mosaic_i,J)
               !These are IN/OUT variables that do not have an LSM equivelant
               ISNOWXY(I,J) = isnowxy_mosaic(I,mosaic_i,J)
               TVXY(I,J) = tvxy_mosaic(I,mosaic_i,J)
               TGXY(I,J) = tgxy_mosaic(I,mosaic_i,J)
               CANICEXY(I,J) = canicexy_mosaic(I,mosaic_i,J)
               CANLIQXY(I,J) = canliqxy_mosaic(I,mosaic_i,J)
               EAHXY(I,J) = eahxy_mosaic(I,mosaic_i,J)
               TAHXY(I,J) = tahxy_mosaic(I,mosaic_i,J)
               CMXY(I,J) = cmxy_mosaic(I,mosaic_i,J)
               CHXY(I,J) = chxy_mosaic(I,mosaic_i,J)
               FWETXY(I,J) = fwetxy_mosaic(I,mosaic_i,J)
               SNEQVOXY(I,J) = sneqvoxy_mosaic(I,mosaic_i,J)
               ALBOLDXY(I,J) = alboldxy_mosaic(I,mosaic_i,J)
               QSNOWXY(I,J) = qsnowxy_mosaic(I,mosaic_i,J)
               QRAINXY(I,J) = qrainxy_mosaic(I,mosaic_i,J)
               WSLAKEXY(I,J) = wslakexy_mosaic(I,mosaic_i,J)
               ZWTXY(I,J) = zwtxy_mosaic(I,mosaic_i,J)
               WAXY(I,J) = waxy_mosaic(I,mosaic_i,J)
               WTXY(I,J) = wtxy_mosaic(I,mosaic_i,J)
               LFMASSXY(I,J) = lfmassxy_mosaic(I,mosaic_i,J)
               RTMASSXY(I,J) = rtmassxy_mosaic(I,mosaic_i,J)
               STMASSXY(I,J) = stmassxy_mosaic(I,mosaic_i,J)
               WOODXY(I,J) = woodxy_mosaic(I,mosaic_i,J)
               STBLCPXY(I,J) = stblcpxy_mosaic(I,mosaic_i,J)
               FASTCPXY(I,J) = fastcpxy_mosaic(I,mosaic_i,J)
               XLAIXY(I,J) = xlai_mosaic(I,mosaic_i,J)
               XSAIXY(I,J) = xsaixy_mosaic(I,mosaic_i,J)
               DEEPRECHXY(I,J) = deeprechxy_mosaic(I,mosaic_i,J)
               RECHXY(I,J) = rechxy_mosaic(I,mosaic_i,J)
               !SNOW VARIABLES

               DO LAYER=1,3
                   TSNOXY(I,LAYER-3,J) = tsnoxy_mosaic(I,3*(mosaic_i - 1)+LAYER,J)
                   SNICEXY(I,LAYER-3,J) = snicexy_mosaic(I,3*(mosaic_i - 1)+LAYER,J)
                   SNLIQXY(I,LAYER-3,J) = snliqxy_mosaic(I,3*(mosaic_i - 1)+LAYER,J)
               ENDDO

               DO LAYER=1,7
                   ZSNSOXY(I,LAYER-3,J) = zsnsoxy_mosaic(I,7*(mosaic_i - 1)+LAYER,J)
               ENDDO

               ! additional IN/OUT variables that have to be added
               ACC_SSOILXY(I,J) = ACC_SSOILXY_mosaic(I,mosaic_i,J)
               ACC_QINSURXY(I,J) = ACC_QINSURXY_mosaic(I,mosaic_i,J)
               ACC_QSEVAXY(I,J) = ACC_QSEVAXY_mosaic(I,mosaic_i,J)
               ACC_DWATERXY(I,J) = ACC_DWATERXY_mosaic(I,mosaic_i,J)
               ACC_PRCPXY(I,J) = ACC_PRCPXY_mosaic(I,mosaic_i,J)
               ACC_ECANXY(I,J) = ACC_ECANXY_mosaic(I,mosaic_i,J)
               ACC_ETRANXY(I,J) = ACC_ETRANXY_mosaic(I,mosaic_i,J)
               ACC_EDIRXY(I,J) = ACC_EDIRXY_mosaic(I,mosaic_i,J)

               !2D inout irrigation variables
               IRNUMSI(I,J) = IRNUMSI_mosaic(I,mosaic_i,J)
               IRNUMMI(I,J) = IRNUMMI_mosaic(I,mosaic_i,J)
               IRNUMFI(I,J) = IRNUMFI_mosaic(I,mosaic_i,J)
               IRWATSI(I,J) = IRWATSI_mosaic(I,mosaic_i,J)
               IRWATMI(I,J) = IRWATMI_mosaic(I,mosaic_i,J)
               IRWATFI(I,J) = IRWATFI_mosaic(I,mosaic_i,J)
               IRELOSS(I,J) = IRELOSS_mosaic(I,mosaic_i,J)
               IRSIVOL(I,J) = IRSIVOL_mosaic(I,mosaic_i,J)
               IRMIVOL(I,J) = IRMIVOL_mosaic(I,mosaic_i,J)
               IRFIVOL(I,J) = IRFIVOL_mosaic(I,mosaic_i,J)
               IRRSPLH(I,J) = IRRSPLH_mosaic(I,mosaic_i,J)

               ! do not need to add any of the OUTPUT values. These will be captured in the post processing steps!

               !Add in a (WE got here) Print statement

               !-----------------------------------------------------------------------
               ! insert the NOAHMP MODEL and the URBAN models
               !-----------------------------------------------------------------------
     ! add break down if this is just a mosaic scheme, or if it is a mosaic and HUE
     ! scheme. If there is no HUE option, we continue as normally

     !     2D to 1D

     ! IN only

 ! for using soil update timestep difference from noahmp main timestep
         calculate_soil = .false.
         soil_update_steps = nint(soilstep/DT)  ! 3600 = 1 hour
         soil_update_steps = max(soil_update_steps,1)
         if ( soil_update_steps == 1 ) then
           ACC_SSOILXY  = 0.0
           ACC_QINSURXY = 0.0
           ACC_QSEVAXY  = 0.0
           ACC_ETRANIXY = 0.0
           ACC_DWATERXY = 0.0
           ACC_PRCPXY   = 0.0
           ACC_ECANXY   = 0.0
           ACC_ETRANXY  = 0.0
           ACC_EDIRXY   = 0.0
         end if
         if ( soil_update_steps > 1 ) then
          if ( mod(itimestep,soil_update_steps) == 1 ) then
           ACC_SSOILXY  = 0.0
           ACC_QINSURXY = 0.0
           ACC_QSEVAXY  = 0.0
           ACC_ETRANIXY = 0.0
           ACC_DWATERXY = 0.0
           ACC_PRCPXY   = 0.0
           ACC_ECANXY   = 0.0
           ACC_ETRANXY  = 0.0
           ACC_EDIRXY   = 0.0
          end if
         end if

         if (mod(itimestep,soil_update_steps) == 0) calculate_soil = .true.
     ! end soil timestep


            COSZ   = COSZIN  (I,J)                         ! cos zenith angle []
            LAT    = XLAT  (I,J)                           ! latitude [rad]
            Z_ML   = 0.5*DZ8W(I,1,J)                       ! DZ8W: thickness of full levels; ZLVL forcing height [m]
            VEGTYP = IVGTYP(I,J)                           ! vegetation type
            if(iopt_soil == 1) then
              SOILTYP= ISLTYP(I,J)                         ! soil type same in all layers
            elseif(iopt_soil == 2) then
              SOILTYP(1) = nint(SOILCL1(I,J))              ! soil type in layer1
              SOILTYP(2) = nint(SOILCL2(I,J))              ! soil type in layer2
              SOILTYP(3) = nint(SOILCL3(I,J))              ! soil type in layer3
              SOILTYP(4) = nint(SOILCL4(I,J))              ! soil type in layer4
            elseif(iopt_soil == 3) then
              SOILTYP= ISLTYP(I,J)                         ! to initialize with default
            end if
            FVEG   = VEGFRA(I,J)/100.                      ! vegetation fraction [0-1]
            FVGMAX = VEGMAX (I,J)/100.                     ! Vegetation fraction annual max [0-1]
            TBOT = TMN(I,J)                                ! Fixed deep soil temperature for land
            T_ML   = T3D(I,1,J)                            ! temperature defined at intermediate level [K]
            Q_ML   = QV3D(I,1,J)/(1.0+QV3D(I,1,J))         ! convert from mixing ratio to specific humidity [kg/kg]
            U_ML   = U_PHY(I,1,J)                          ! u-wind at interface [m/s]
            V_ML   = V_PHY(I,1,J)                          ! v-wind at interface [m/s]
            SWDN   = SWDOWN(I,J)                           ! shortwave down from SW scheme [W/m2]
            LWDN   = GLW(I,J)                              ! total longwave down from LW scheme [W/m2]
            P_ML   =(P8W3D(I,KTS+1,J)+P8W3D(I,KTS,J))*0.5  ! surface pressure defined at intermediate level [Pa]
     	                                              !    consistent with temperature, mixing ratio
            PSFC   = P8W3D(I,1,J)                          ! surface pressure defined a full levels [Pa]
            PRCP   = PRECIP_IN (I,J) / DT                  ! timestep total precip rate (glacier) [mm/s]! MB: v3.7

            CROPTYPE = 0
            IF (IOPT_CROP > 0 .AND. VEGTYP == ISCROP_TABLE) CROPTYPE = DEFAULT_CROP_TABLE ! default croptype is generic dynamic vegetation crop
            IF (IOPT_CROP > 0 .AND. CROPCAT(I,J) > 0) THEN
              CROPTYPE = CROPCAT(I,J)                      ! crop type
     	 VEGTYP = ISCROP_TABLE
              FVGMAX = 0.95
     	 FVEG   = 0.95
            END IF

            IF (PRESENT(MP_RAINC) .AND. PRESENT(MP_RAINNC) .AND. PRESENT(MP_SHCV) .AND. &
                PRESENT(MP_SNOW)  .AND. PRESENT(MP_GRAUP)  .AND. PRESENT(MP_HAIL)   ) THEN

              PRCPCONV  = MP_RAINC (I,J)/DT                ! timestep convective precip rate [mm/s]     ! MB: v3.7
              PRCPNONC  = MP_RAINNC(I,J)/DT                ! timestep non-convective precip rate [mm/s] ! MB: v3.7
              PRCPSHCV  = MP_SHCV(I,J)  /DT                ! timestep shallow conv precip rate [mm/s]   ! MB: v3.7
              PRCPSNOW  = MP_SNOW(I,J)  /DT                ! timestep snow precip rate [mm/s]           ! MB: v3.7
              PRCPGRPL  = MP_GRAUP(I,J) /DT                ! timestep graupel precip rate [mm/s]        ! MB: v3.7
              PRCPHAIL  = MP_HAIL(I,J)  /DT                ! timestep hail precip rate [mm/s]           ! MB: v3.7

              PRCPOTHR  = PRCP - PRCPCONV - PRCPNONC - PRCPSHCV ! take care of other (fog) contained in rainbl
     	 PRCPOTHR  = MAX(0.0,PRCPOTHR)
     	 PRCPNONC  = PRCPNONC + PRCPOTHR
              PRCPSNOW  = PRCPSNOW + SR(I,J)  * PRCPOTHR
            ELSE
              PRCPCONV  = 0.
              PRCPNONC  = PRCP
              PRCPSHCV  = 0.
              PRCPSNOW  = SR(I,J) * PRCP
              PRCPGRPL  = 0.
              PRCPHAIL  = 0.
            ENDIF

     ! IN/OUT fields
            ! all HUE options that re needed (note that this could cause a memeory issue if 41 is right on the edge of the defined mosaic layers)
             IF (IOPT_HUE.eq.1) THEN

               IF (VEGTYP.EQ.41) THEN !Added by Aaron Alexander
                 ! This branch swaps the soil moisture for the canpoy of impermeable pavement
                 ! This assumes that the ordering is source (i.e. fully vegetated) is second. So its _pavement_ | _fully Vegetated_

                 SMC  (       1:NSOIL) = SMOIS(I, 1:NSOIL, J ) !  These are the soil moisture values of the pavement
                 SMH2O(       1:NSOIL) = SH2O( I, 1:NSOIL, J ) !  Soil water content that is under the pavement

                 DO LAYER=1,NSOIL
                   SMC_intermediate(LAYER) = SMOIS_mosaic(I+1,NSOIL*(mosaic_i)+LAYER,J)  !  This is the soil moisture of the fully vegetated square (one mosaic tile down)
                   SH2O_intermediate(LAYER) = SH2O_mosaic(I+1,NSOIL*(mosaic_i)+LAYER,J)  !  This is the soil moisture of the fully vegetated square (one mosaic tile down)

                 END DO

               ELSE ! All other land types that do not require the sharing of information across simulated soil areas
                 SMC  (      1:NSOIL)  = SMOIS   (I,      1:NSOIL,J)  ! soil total moisture [m3/m3]
                 SMH2O(      1:NSOIL)  = SH2O    (I,      1:NSOIL,J)  ! soil liquid moisture [m3/m3]

                 SMC_intermediate (1:NSOIL) = SMOIS(I,1:NSOIL,J)
                 SH2O_intermediate(1:NSOIL) = SH2O(I,1:NSOIL,J)

               ENDIF ! end soil moisture

               IF (VEGTYP.EQ.46) THEN
                 ! GREEN ROOF
                 DETENTION_STORAGE = DETENTION_STORAGEXY_mosaic(I,mosaic_i,J)
               ELSE
                 DETENTION_STORAGE = 0.
               END IF

               ! Turf grass flux runon
               IF( (VEGTYP.EQ.43).OR.(VEGTYP.EQ.44)) THEN
                 VOL_FLUX_RUNON = VOL_FLUX_RUNONXY_mosaic(I,mosaic_i-1,J) ! you want the previous one
               ELSE
                 VOL_FLUX_RUNON = 0.
               END IF

               ! soil moisture pulling flux from the previous call
               IF (VEGTYP.EQ.47) THEN
                 DO LAYER = 1,NSOIL
                    VOL_FLUX_SM(LAYER) = VOL_FLUX_SMXY_mosaic(I,NSOIL*(mosaic_i-2)+LAYER,J)
                 END DO
               ELSE
                 VOL_FLUX_SM(1:NSOIL) = 0.
               END IF

               !Scaling volume (FAREA of current )

               FAREAXY = LANDUSEF2(I,mosaic_i,J)

               !runonsf (TO KEEP TRACK OF EVERYTHING that is lateral transfer)
               RUNONSRF = RUNONSFXY_mosaic(I,mosaic_i,J)

             ELSE
               SMC  (      1:NSOIL)  = SMOIS   (I,      1:NSOIL,J)  ! soil total moisture [m3/m3]
               SMH2O(      1:NSOIL)  = SH2O    (I,      1:NSOIL,J)  ! soil liquid moisture [m3/m3]

               SMC_intermediate (1:NSOIL) = SMOIS(I,1:NSOIL,J)
               SH2O_intermediate(1:NSOIL) = SH2O(I,1:NSOIL,J)
             END IF !if else hue noahmp vs. just mosaic

            ISNOW                 = ISNOWXY (I,J)                ! snow layers []

            STC  (-NSNOW+1:    0) = TSNOXY  (I,-NSNOW+1:    0,J) ! snow temperatures [K]
            STC  (      1:NSOIL)  = TSLB    (I,      1:NSOIL,J)  ! soil temperatures [K]
            SWE                   = SNOW    (I,J)                ! snow water equivalent [mm]
            SNDPTH                = SNOWH   (I,J)                ! snow depth [m]
            QSFC1D                = QSFC    (I,J)

     ! INOUT (with no Noah LSM equivalent)

            TV                    = TVXY    (I,J)                ! leaf temperature [K]
            TG                    = TGXY    (I,J)                ! ground temperature [K]
            CANLIQ                = CANLIQXY(I,J)                ! canopy liquid water [mm]
            CANICE                = CANICEXY(I,J)                ! canopy frozen water [mm]
            EAH                   = EAHXY   (I,J)                ! canopy vapor pressure [Pa]
            TAH                   = TAHXY   (I,J)                ! canopy temperature [K]
            CM                    = CMXY    (I,J)                ! avg. momentum exchange (MP only) [m/s]
            CH                    = CHXY    (I,J)                ! avg. heat exchange (MP only) [m/s]
            FWET                  = FWETXY  (I,J)                ! canopy fraction wet or snow
            SNEQVO                = SNEQVOXY(I,J)                ! SWE previous timestep
            ALBOLD                = ALBOLDXY(I,J)                ! albedo previous timestep, for snow aging
            QSNOW                 = QSNOWXY (I,J)                ! snow falling on ground
            QRAIN                 = QRAINXY (I,J)                ! rain falling on ground
            WSLAKE                = WSLAKEXY(I,J)                ! lake water storage (can be neg.) (mm)
            ZWT                   = ZWTXY   (I,J)                ! depth to water table [m]
            WA                    = WAXY    (I,J)                ! water storage in aquifer [mm]
            WT                    = WTXY    (I,J)                ! water in aquifer&saturated soil [mm]
            ZSNSO(-NSNOW+1:NSOIL) = ZSNSOXY (I,-NSNOW+1:NSOIL,J) ! depth to layer interface
            SNICE(-NSNOW+1:    0) = SNICEXY (I,-NSNOW+1:    0,J) ! snow layer ice content
            SNLIQ(-NSNOW+1:    0) = SNLIQXY (I,-NSNOW+1:    0,J) ! snow layer water content
            LFMASS                = LFMASSXY(I,J)                ! leaf mass
            RTMASS                = RTMASSXY(I,J)                ! root mass
            STMASS                = STMASSXY(I,J)                ! stem mass
            WOOD                  = WOODXY  (I,J)                ! mass of wood (incl. woody roots) [g/m2]
            STBLCP                = STBLCPXY(I,J)                ! stable carbon pool
            FASTCP                = FASTCPXY(I,J)                ! fast carbon pool
            PLAI                  = XLAIXY  (I,J)                ! leaf area index [-] (no snow effects)
            PSAI                  = XSAIXY  (I,J)                ! stem area index [-] (no snow effects)
            TAUSS                 = TAUSSXY (I,J)                ! non-dimensional snow age
            SMCEQ(       1:NSOIL) = SMOISEQ (I,       1:NSOIL,J)
            SMCWTD                = SMCWTDXY(I,J)
            RECH                  = 0.
            DEEPRECH              = 0.
            ACC_SSOIL             = ACC_SSOILXY (I,J)                 ! surface heat flux
            ACC_QSEVA             = ACC_QSEVAXY (I,J)
            ACC_QINSUR            = ACC_QINSURXY(I,J)
            ACC_ETRANI            = ACC_ETRANIXY(I,:,J)
            ACC_DWATER            = ACC_DWATERXY(I,J)
            ACC_PRCP              = ACC_PRCPXY  (I,J)
            ACC_ECAN              = ACC_ECANXY  (I,J)
            ACC_ETRAN             = ACC_ETRANXY (I,J)
            ACC_EDIR              = ACC_EDIRXY  (I,J)

     ! tile drainage
            QTLDRN                = 0.                           ! tile drainage (mm)
            TDFRACMP              = TD_FRACTION(I,J)             ! tile drainage map
#ifdef WRF_HYDRO
            WATBLED               = ZWATBLE2D (I,J)              ! water table depth for tile drainage
#endif
     ! irrigation vars
            IRRFRA                = IRFRACT(I,J)    ! irrigation fraction
            SIFAC                 = SIFRACT(I,J)    ! sprinkler irrigation fraction
            MIFAC                 = MIFRACT(I,J)    ! micro irrigation fraction
            FIFAC                 = FIFRACT(I,J)    ! flood irrigation fraction
            IRCNTSI               = IRNUMSI(I,J)    ! irrigation event number, Sprinkler
            IRCNTMI               = IRNUMMI(I,J)    ! irrigation event number, Micro
            IRCNTFI               = IRNUMFI(I,J)    ! irrigation event number, Flood
            IRAMTSI               = IRWATSI(I,J)    ! irrigation water amount [m] to be applied, Sprinkler
            IRAMTMI               = IRWATMI(I,J)    ! irrigation water amount [m] to be applied, Micro
            IRAMTFI               = IRWATFI(I,J)    ! irrigation water amount [m] to be applied, Flood
            IREVPLOS              = 0.0             ! loss of irrigation water to evaporation,sprinkler [m/timestep]
            IRSIRATE              = 0.0             ! rate of irrigation by sprinkler (mm)
            IRMIRATE              = 0.0             ! rate of irrigation by micro (mm)
            IRFIRATE              = 0.0             ! rate of irrigation by micro (mm)
            FIRR                  = 0.0             ! latent heating due to sprinkler evaporation (W m-2)
            EIRR                  = 0.0             ! evaporation from sprinkler (mm/s)

            if(iopt_crop == 2) then   ! gecros crop model

              gecros1d(1:60)      = gecros_state(I,1:60,J)       ! Gecros variables 2D -> local

              if(croptype == 1) then
                gecros_dd   =  2.5
                gecros_tbem =  2.0
                gecros_emb  = 10.2
                gecros_ema  = 40.0
                gecros_ds1  =  2.1 !BBCH 92
                gecros_ds2  =  2.0 !BBCH 90
                gecros_ds1x =  0.0
                gecros_ds2x = 10.0
              end if

              if(croptype == 2) then
                gecros_dd   =  5.0
                gecros_tbem =  8.0
                gecros_emb  = 15.0
                gecros_ema  =  6.0
                gecros_ds1  =  1.78  !BBCH 85
                gecros_ds2  =  1.63  !BBCH 80
                gecros_ds1x =  0.0
                gecros_ds2x = 14.0
              end if

            end if

            SLOPETYP     = 1                               ! set underground runoff slope term
            IST          = 1                               ! MP surface type: 1 = land; 2 = lake
            SOILCOLOR    = 4                               ! soil color: assuming a middle color category ?????????

            IF(any(SOILTYP == 14) .AND. XICE(I,J) == 0.) THEN
               IF(IPRINT) PRINT *, ' SOIL TYPE FOUND TO BE WATER AT A LAND-POINT'
               IF(IPRINT) PRINT *, i,j,'RESET SOIL in surfce.F'
               SOILTYP = 7
            ENDIF
              IF( IVGTYP(I,J) == ISURBAN_TABLE    .or. IVGTYP(I,J) == LCZ_1_TABLE .or. IVGTYP(I,J) == LCZ_2_TABLE .or. &
                  IVGTYP(I,J) == LCZ_3_TABLE      .or. IVGTYP(I,J) == LCZ_4_TABLE .or. IVGTYP(I,J) == LCZ_5_TABLE .or. &
                  IVGTYP(I,J) == LCZ_6_TABLE      .or. IVGTYP(I,J) == LCZ_7_TABLE .or. IVGTYP(I,J) == LCZ_8_TABLE .or. &
                  IVGTYP(I,J) == LCZ_9_TABLE      .or. IVGTYP(I,J) == LCZ_10_TABLE .or. IVGTYP(I,J) == LCZ_11_TABLE ) THEN


              IF(SF_URBAN_PHYSICS == 0 ) THEN
                VEGTYP = ISURBAN_TABLE
              ELSE
                VEGTYP = NATURAL_TABLE  ! set urban vegetation type based on table natural
                FVGMAX = 0.96
              ENDIF

            ENDIF

            CALL TRANSFER_MP_PARAMETERS(VEGTYP,SOILTYP,SLOPETYP,SOILCOLOR,CROPTYPE,parameters)
            if(iopt_soil == 3 .and. .not. parameters%urban_flag) then
     	sand = 0.01 * soilcomp(i,1:4,j)
     	clay = 0.01 * soilcomp(i,5:8,j)
             orgm = 0.0
             if(opt_pedo == 1) call pedotransfer_sr2006(nsoil,sand,clay,orgm,parameters)
            end if
            GRAIN = GRAINXY (I,J)                ! mass of grain XING [g/m2]
            GDD   = GDDXY (I,J)                  ! growing degree days XING
            PGS   = PGSXY (I,J)                  ! growing degree days XING
            if(iopt_crop == 1 .and. croptype > 0) then
              parameters%PLTDAY = PLANTING(I,J)
     	 parameters%HSDAY  = HARVEST (I,J)
     	 parameters%GDDS1  = SEASON_GDD(I,J) / 1770.0 * parameters%GDDS1
     	 parameters%GDDS2  = SEASON_GDD(I,J) / 1770.0 * parameters%GDDS2
     	 parameters%GDDS3  = SEASON_GDD(I,J) / 1770.0 * parameters%GDDS3
     	 parameters%GDDS4  = SEASON_GDD(I,J) / 1770.0 * parameters%GDDS4
     	 parameters%GDDS5  = SEASON_GDD(I,J) / 1770.0 * parameters%GDDS5
            end if
            if(iopt_irr == 2) then
              parameters%PLTDAY = PLANTING(I,J)
              parameters%HSDAY  = HARVEST (I,J)
            end if

     !=== hydrological processes for vegetation in urban model ===
     !=== irrigate vegetaion only in urban area, MAY-SEP, 9-11pm
              IF( IVGTYP(I,J) == ISURBAN_TABLE    .or. IVGTYP(I,J) == LCZ_1_TABLE .or. IVGTYP(I,J) == LCZ_2_TABLE .or. &
                    IVGTYP(I,J) == LCZ_3_TABLE      .or. IVGTYP(I,J) == LCZ_4_TABLE .or. IVGTYP(I,J) == LCZ_5_TABLE .or. &
                    IVGTYP(I,J) == LCZ_6_TABLE      .or. IVGTYP(I,J) == LCZ_7_TABLE .or. IVGTYP(I,J) == LCZ_8_TABLE .or. &
                    IVGTYP(I,J) == LCZ_9_TABLE      .or. IVGTYP(I,J) == LCZ_10_TABLE .or. IVGTYP(I,J) == LCZ_11_TABLE ) THEN

              IF(SF_URBAN_PHYSICS > 0 .AND. IRI_SCHEME == 1 ) THEN
     	     SOLAR_TIME = (JULIAN - INT(JULIAN))*24 + XLONG(I,J)/15.0
     	     IF(SOLAR_TIME < 0.) SOLAR_TIME = SOLAR_TIME + 24.
                  CALL CAL_MON_DAY(INT(JULIAN),YR,JMONTH,JDAY)
                  IF (SOLAR_TIME >= 21. .AND. SOLAR_TIME <= 23. .AND. JMONTH >= 5 .AND. JMONTH <= 9) THEN
                     SMC(1) = max(SMC(1),parameters%SMCREF(1))
                     SMC(2) = max(SMC(2),parameters%SMCREF(2))
                  ENDIF
              ENDIF

            ENDIF
     ! Initialized local

            FICEOLD = 0.0
            FICEOLD(ISNOW+1:0) = SNICEXY(I,ISNOW+1:0,J) &  ! snow ice fraction
                /(SNICEXY(I,ISNOW+1:0,J)+SNLIQXY(I,ISNOW+1:0,J))
            CO2PP  = CO2_TABLE * P_ML                      ! partial pressure co2 [Pa]
            O2PP   = O2_TABLE  * P_ML                      ! partial pressure  o2 [Pa]
            FOLN   = 1.0                                   ! for now, set to nitrogen saturation
            QC     = undefined_value                       ! test dummy value
            PBLH   = undefined_value                       ! test dummy value ! PBL height
            DZ8W1D = DZ8W (I,1,J)                          ! thickness of atmospheric layers

            IF(VEGTYP == 25) FVEG = 0.0                  ! Set playa, lava, sand to bare
            IF(VEGTYP == 25) PLAI = 0.0
            IF(VEGTYP == 26) FVEG = 0.0                  ! hard coded for USGS
            IF(VEGTYP == 26) PLAI = 0.0
            IF(VEGTYP == 27) FVEG = 0.0
            IF(VEGTYP == 27) PLAI = 0.0

            IF ( VEGTYP == ISICE_TABLE ) THEN
              ICE = -1                           ! Land-ice point
              CALL NOAHMP_OPTIONS_GLACIER(IOPT_ALB  ,IOPT_SNF  ,IOPT_TBOT, IOPT_STC, IOPT_GLA )

              TBOT = MIN(TBOT,263.15)                      ! set deep temp to at most -10C
              CALL NOAHMP_GLACIER(     I,       J,    COSZ,   NSNOW,   NSOIL,      DT, & ! IN : Time/Space/Model-related
                                    T_ML,    P_ML,    U_ML,    V_ML,    Q_ML,    SWDN, & ! IN : Forcing
                                    PRCP,    LWDN,    TBOT,    Z_ML, FICEOLD,   ZSOIL, & ! IN : Forcing
                                   QSNOW,  SNEQVO,  ALBOLD,      CM,      CH,   ISNOW, & ! IN/OUT :
                                     SWE,     SMC,   ZSNSO,  SNDPTH,   SNICE,   SNLIQ, & ! IN/OUT :
                                      TG,     STC,   SMH2O,   TAUSS,  QSFC1D,          & ! IN/OUT :
                                     FSA,     FSR,    FIRA,     FSH,    FGEV,   SSOIL, & ! OUT :
                                    TRAD,   ESOIL,   RUNSF,   RUNSB,     SAG,    SALB, & ! OUT :
                                   QSNBOT,PONDING,PONDING1,PONDING2,    T2MB,    Q2MB, & ! OUT :
     			      EMISSI,  FPICE,    CHB2,   QMELT                    & ! OUT :
#ifdef WRF_HYDRO
                                   , sfcheadrt(i,j)                                      &
#endif
                                   )

              FSNO   = 1.0
              TV     = undefined_value     ! Output from standard Noah-MP undefined for glacier points
              TGB    = TG
              CANICE = undefined_value
              CANLIQ = undefined_value
              EAH    = undefined_value
              TAH    = undefined_value
              FWET   = undefined_value
              WSLAKE = undefined_value
     !         ZWT    = undefined_value
              WA     = undefined_value
              WT     = undefined_value
              LFMASS = undefined_value
              RTMASS = undefined_value
              STMASS = undefined_value
              WOOD   = undefined_value
              QTLDRN = undefined_value
              GRAIN  = undefined_value
              GDD    = undefined_value
              STBLCP = undefined_value
              FASTCP = undefined_value
              PLAI   = undefined_value
              PSAI   = undefined_value
              T2MV   = undefined_value
              Q2MV   = undefined_value
              NEE    = undefined_value
              GPP    = undefined_value
              NPP    = undefined_value
              FVEGMP = 0.0
              ECAN   = 0.0
              ETRAN  = 0.0
              APAR   = undefined_value
              PSN    = undefined_value
              SAV    = 0.0
              RSSUN  = undefined_value
              RSSHA  = undefined_value
              RB     = undefined_value
              LAISUN = undefined_value
              LAISHA = undefined_value
              RS(I,J)= undefined_value
              BGAP   = undefined_value
              WGAP   = undefined_value
              TGV    = undefined_value
              CHV    = undefined_value
              CHB    = CH
              IRC    = 0.0
              IRG    = 0.0
              SHC    = 0.0
              SHG    = 0.0
              EVG    = 0.0
              GHV    = 0.0
              CANHS  = 0.0
              IRB    = FIRA
              SHB    = FSH
              EVB    = FGEV
              GHB    = SSOIL
              TR     = 0.0
              EVC    = 0.0
              PAH    = 0.0
              PAHG   = 0.0
              PAHB   = 0.0
              PAHV   = 0.0
              CHLEAF = undefined_value
              CHUC   = undefined_value
              CHV2   = undefined_value
              FCEV   = 0.0
              FCTR   = 0.0
              Z0WRF  = 0.002
              QFX(I,J) = ESOIL
              LH (I,J) = FGEV
              QINTS  = 0.0
              QINTR  = 0.0
              QDRIPS = 0.0
              QDRIPR = 0.0
              QTHROS = PRCP * FPICE
              QTHROR = PRCP * (1.0 - FPICE)
              QSNSUB = MAX( ESOIL, 0.)
              QSNFRO = ABS( MIN(ESOIL, 0.))
              QSUBC = 0.0
              QFROC = 0.0
              QFRZC = 0.0
              QMELTC = 0.0
              QEVAC = 0.0
              QDEWC = 0.0
              RAININ = PRCP * (1.0 - FPICE)
              SNOWIN = PRCP * FPICE
              CANICE = 0.0
              CANLIQ = 0.0
              QTLDRN = 0.0
              RUNSF  = RUNSF * dt
              RUNSB  = RUNSB * dt
         ELSE
              ICE=0                              ! Neither sea ice or land ice.
              CALL NOAHMP_SFLX (parameters, &
                 I       , J       , LAT     , YEARLEN , JULIAN  , COSZ    , & ! IN : Time/Space-related
                 DT      , DX      , DZ8W1D  , NSOIL   , ZSOIL   , NSNOW   , & ! IN : Model configuration
                 FVEG    , FVGMAX  , VEGTYP  , ICE     , IST     , CROPTYPE, & ! IN : Vegetation/Soil characteristics
                 SMCEQ   ,                                                   & ! IN : Vegetation/Soil characteristics
                 T_ML    , P_ML    , PSFC    , U_ML    , V_ML    , Q_ML    , & ! IN : Forcing
                 QC      , SWDN    , LWDN    ,                               & ! IN : Forcing
     	    PRCPCONV, PRCPNONC, PRCPSHCV, PRCPSNOW, PRCPGRPL, PRCPHAIL, & ! IN : Forcing
                 TBOT    , CO2PP   , O2PP    , FOLN    , FICEOLD , Z_ML    , & ! IN : Forcing
                 IRRFRA  , SIFAC   , MIFAC   , FIFAC   , LLANDUSE,           & ! IN : Irrigation: fractions
                 ALBOLD  , SNEQVO  ,                                         & ! IN/OUT :
                 STC     , SMH2O   , SMC     , TAH     , EAH     , FWET    , & ! IN/OUT :
                 CANLIQ  , CANICE  , TV      , TG      , QSFC1D  , QSNOW   , & ! IN/OUT :
                 QRAIN   ,                                                   & ! IN/OUT :
                 ISNOW   , ZSNSO   , SNDPTH  , SWE     , SNICE   , SNLIQ   , & ! IN/OUT :
                 ZWT     , WA      , WT      , WSLAKE  , LFMASS  , RTMASS  , & ! IN/OUT :
                 STMASS  , WOOD    , STBLCP  , FASTCP  , PLAI    , PSAI    , & ! IN/OUT :
                 CM      , CH      , TAUSS   ,                               & ! IN/OUT :
                 GRAIN   , GDD     , PGS     ,                               & ! IN/OUT
                 SMCWTD  ,DEEPRECH , RECH    ,                               & ! IN/OUT :
                 GECROS1D,                                                   & ! IN/OUT :
                 QTLDRN  , TDFRACMP,                                         & ! IN/OUT : tile drainage
                 Z0WRF   ,                                                   & ! OUT :
                 IRCNTSI , IRCNTMI , IRCNTFI , IRAMTSI , IRAMTMI , IRAMTFI , & ! IN/OUT : Irrigation: vars
                 IRSIRATE, IRMIRATE, IRFIRATE, FIRR    , EIRR    ,           & ! IN/OUT : Irrigation: vars
                 FSA     , FSR     , FIRA    , FSH     , SSOIL   , FCEV    , & ! OUT :
                 FGEV    , FCTR    , ECAN    , ETRAN   , ESOIL   , TRAD    , & ! OUT :
                 TGB     , TGV     , T2MV    , T2MB    , Q2MV    , Q2MB    , & ! OUT :
                 RUNSF   , RUNSB   , APAR    , PSN     , SAV     , SAG     , & ! OUT :
                 FSNO    , NEE     , GPP     , NPP     , FVEGMP  , SALB    , & ! OUT :
                 QSNBOT  , PONDING , PONDING1, PONDING2, RSSUN   , RSSHA   , & ! OUT :
                 ALBSND  , ALBSNI  ,                                         & ! OUT :
                 BGAP    , WGAP    , CHV     , CHB     , EMISSI  ,           & ! OUT :
                 SHG     , SHC     , SHB     , EVG     , EVB     , GHV     , & ! OUT :
     	    GHB     , IRG     , IRC     , IRB     , TR      , EVC     , & ! OUT :
     	    CHLEAF  , CHUC    , CHV2    , CHB2    , FPICE   , PAHV    , & ! OUT :
                 PAHG    , PAHB    , PAH     , LAISUN  , LAISHA  , RB      , & ! OUT :
                 QINTS   , QINTR   , QDRIPS  , QDRIPR  , QTHROS  , QTHROR  , & ! OUT :
                 QSNSUB  , QSNFRO  , QSUBC   , QFROC   , QFRZC   , QMELTC  , & ! OUT :
                 QEVAC   , QDEWC   , QMELT   ,                               & ! OUT :
                 RAININ  , SNOWIN  , ACC_SSOIL, ACC_QINSUR, ACC_QSEVA      , & ! OUT :
                 ACC_ETRANI, HCPCT , EFLXB   , CANHS   ,                     & ! OUT :
                 ACC_DWATER, ACC_PRCP, ACC_ECAN, ACC_ETRAN, ACC_EDIR       , & ! INOUT
                 SMC_INTERMEDIATE, SH2O_INTERMEDIATE, BTRANI_dummy, RUNONSRF, NSOIL_GR,                  & ! addtiional variables for HUE added by Aaron A.
                DETENTION_STORAGE,FAREAXY,VOL_FLUX_SM,VOL_FLUX_RUNON &
#ifdef WRF_HYDRO
                 , sfcheadrt(i,j), WATBLED                                   &
#endif
 )            ! OUT :
                 QFX(I,J) = ECAN + ESOIL + ETRAN + EIRR
                 LH(I,J)  = FCEV + FGEV  + FCTR  + FIRR

        ENDIF ! glacial split ends
#ifdef WRF_HYDRO
     !AD_CHANGE: Glacier cells can produce small negative subsurface runoff for mass balance.
     !           This will crash channel routing, so only pass along positive runoff.
                 !soldrain(i,j) = max(RUNSB*dt, 0.0)        !mm , underground runoff
                 soldrain(i,j) = max(RUNSB, 0.0)        !mm , underground runoff
                 INFXSRT(i,j) = RUNSF   !*dt        !mm , surface runoff
                 qtiledrain (i,j) = QTLDRN   !*dt      !mm, tile drainage
#endif


     ! INPUT/OUTPUT
            IF (IOPT_HUE.eq.1) THEN

              ! Changes in the soil moisture with replacement
              IF (VEGTYP.EQ.41) THEN

                SMOIS(I, 1:NSOIL, J ) = SMC  (       1:NSOIL) !  These are the soil moisture values of the pavement
                SH2O( I, 1:NSOIL, J ) = SMH2O(       1:NSOIL) !  Soil water content that is under the pavement

                DO LAYER=1,NSOIL
                  SMOIS_mosaic(I+1,NSOIL*(mosaic_i)+LAYER,J) = SMC_intermediate (LAYER)   !  This is the soil moisture of the fully vegetated square
                  SH2O_mosaic(I+1,NSOIL*(mosaic_i)+LAYER,J)= SH2O_intermediate(LAYER)    !  This is the soil moisture of the fully vegetated square

                END DO
              ELSE
                SMOIS(I, 1:NSOIL, J ) = SMC(1:NSOIL)
                SH2O( I, 1:NSOIL, J ) = SMH2O(1:NSOIL)
              END IF

              ! the green roof detention storage
              DETENTION_STORAGEXY_mosaic(I,mosaic_i,J) = DETENTION_STORAGE
              ! runon flux
              VOL_FLUX_RUNONXY_mosaic(I,mosaic_i,J) = VOL_FLUX_RUNON
              ! soil moisutre flux
              DO LAYER=1,NSOIL
                VOL_FLUX_SMXY_mosaic(I,NSOIL*(mosaic_i-1)+LAYER,J) = VOL_FLUX_SM(LAYER)
              END DO

            ELSE
              SMOIS(I, 1:NSOIL, J ) = SMC(1:NSOIL)
              SH2O( I, 1:NSOIL, J ) = SMH2O(1:NSOIL)
            END IF




                  TSK      (I,J)                = TRAD
                  HFX      (I,J)                = FSH
                  GRDFLX   (I,J)                = SSOIL
     	     SMSTAV   (I,J)                = 0.0  ! [maintained as Noah consistency]
                  SMSTOT   (I,J)                = 0.0  ! [maintained as Noah consistency]
                  SFCRUNOFF(I,J)                = SFCRUNOFF(I,J) + RUNSF  !* DT
                  UDRUNOFF (I,J)                = UDRUNOFF(I,J)  + RUNSB  !* DT
                  QTDRAIN  (I,J)                = QTDRAIN (I,J)  + QTLDRN !* DT
                  IF ( SALB > -999 ) THEN
                     ALBEDO(I,J)                = SALB
                  ENDIF
                  SNOWC    (I,J)                = FSNO

                  SMOIS    (I,      1:NSOIL,J)  = SMC   (      1:NSOIL)
                  SH2O     (I,      1:NSOIL,J)  = SMH2O (      1:NSOIL)
                  TSLB     (I,      1:NSOIL,J)  = STC   (      1:NSOIL)
                  SNOW     (I,J)                = SWE
                  SNOWH    (I,J)                = SNDPTH
                  CANWAT   (I,J)                = CANLIQ + CANICE
                  ACSNOW   (I,J)                = ACSNOW(I,J) +  PRECIP_IN(I,J) * FPICE
     !             ACSNOM   (I,J)                = ACSNOM(I,J) + QSNBOT*DT + PONDING + PONDING1 + PONDING2
                  ACSNOM   (I,J)                = ACSNOM(I,J) + QMELT*DT + PONDING + PONDING1 + PONDING2
                  EMISS    (I,J)                = EMISSI
                  QSFC     (I,J)                = QSFC1D

                  ISNOWXY  (I,J)                = ISNOW
                  TVXY     (I,J)                = TV
                  TGXY     (I,J)                = TG
                  CANLIQXY (I,J)                = CANLIQ
                  CANICEXY (I,J)                = CANICE
                  EAHXY    (I,J)                = EAH
                  TAHXY    (I,J)                = TAH
                  CMXY     (I,J)                = CM
                  CHXY     (I,J)                = CH
                  FWETXY   (I,J)                = FWET
                  SNEQVOXY (I,J)                = SNEQVO
                  ALBOLDXY (I,J)                = ALBOLD
                  QSNOWXY  (I,J)                = QSNOW
                  QRAINXY  (I,J)                = QRAIN
                  WSLAKEXY (I,J)                = WSLAKE
                  ZWTXY    (I,J)                = ZWT
                  WAXY     (I,J)                = WA
                  WTXY     (I,J)                = WT
                  TSNOXY   (I,-NSNOW+1:    0,J) = STC   (-NSNOW+1:    0)
                  ZSNSOXY  (I,-NSNOW+1:NSOIL,J) = ZSNSO (-NSNOW+1:NSOIL)
                  SNICEXY  (I,-NSNOW+1:    0,J) = SNICE (-NSNOW+1:    0)
                  SNLIQXY  (I,-NSNOW+1:    0,J) = SNLIQ (-NSNOW+1:    0)
                  LFMASSXY (I,J)                = LFMASS
                  RTMASSXY (I,J)                = RTMASS
                  STMASSXY (I,J)                = STMASS
                  WOODXY   (I,J)                = WOOD
                  STBLCPXY (I,J)                = STBLCP
                  FASTCPXY (I,J)                = FASTCP
                  XLAIXY   (I,J)                = PLAI
                  XSAIXY   (I,J)                = PSAI
                  TAUSSXY  (I,J)                = TAUSS

     ! OUTPUT
                  Z0       (I,J)                = Z0WRF
                  ZNT      (I,J)                = Z0WRF
                  T2MVXY   (I,J)                = T2MV
                  T2MBXY   (I,J)                = T2MB
                  Q2MVXY   (I,J)                = Q2MV/(1.0 - Q2MV)  ! specific humidity to mixing ratio
                  Q2MBXY   (I,J)                = Q2MB/(1.0 - Q2MB)  ! consistent with registry def of Q2
                  TRADXY   (I,J)                = TRAD
                  NEEXY    (I,J)                = NEE
                  GPPXY    (I,J)                = GPP
                  NPPXY    (I,J)                = NPP
                  FVEGXY   (I,J)                = FVEGMP
                  RUNSFXY  (I,J)                = RUNSF
                  RUNSBXY  (I,J)                = RUNSB
                  ECANXY   (I,J)                = ECAN
                  EDIRXY   (I,J)                = ESOIL
                  ETRANXY  (I,J)                = ETRAN
                  FSAXY    (I,J)                = FSA
                  FIRAXY   (I,J)                = FIRA
                  APARXY   (I,J)                = APAR
                  PSNXY    (I,J)                = PSN
                  SAVXY    (I,J)                = SAV
                  SAGXY    (I,J)                = SAG
                  RSSUNXY  (I,J)                = RSSUN
                  RSSHAXY  (I,J)                = RSSHA
                  LAISUN                        = MAX(LAISUN, 0.0)
                  LAISHA                        = MAX(LAISHA, 0.0)
                  RB                            = MAX(RB, 0.0)
     ! New Calculation of total Canopy/Stomatal Conductance Based on Bonan et al. (2011)
     ! -- Inverse of Canopy Resistance (below)
                  IF(RSSUN .le. 0.0 .or. RSSHA .le. 0.0 .or. LAISUN .eq. 0.0 .or. LAISHA .eq. 0.0) THEN
                     RS    (I,J)                = 0.0
                  ELSE
                     RS    (I,J)                = ((1.0/(RSSUN+RB)*LAISUN) + ((1.0/(RSSHA+RB))*LAISHA))
                     RS    (I,J)                = 1.0/RS(I,J) !Resistance
                  ENDIF
                  BGAPXY   (I,J)                = BGAP
                  WGAPXY   (I,J)                = WGAP
                  TGVXY    (I,J)                = TGV
                  TGBXY    (I,J)                = TGB
                  CHVXY    (I,J)                = CHV
                  CHBXY    (I,J)                = CHB
                  IRCXY    (I,J)                = IRC
                  IRGXY    (I,J)                = IRG
                  SHCXY    (I,J)                = SHC
                  SHGXY    (I,J)                = SHG
                  EVGXY    (I,J)                = EVG
                  GHVXY    (I,J)                = GHV
                  IRBXY    (I,J)                = IRB
                  SHBXY    (I,J)                = SHB
                  EVBXY    (I,J)                = EVB
                  GHBXY    (I,J)                = GHB
                  canhsxy  (I,J)                = CANHS
                  TRXY     (I,J)                = TR
                  EVCXY    (I,J)                = EVC
                  CHLEAFXY (I,J)                = CHLEAF
                  CHUCXY   (I,J)                = CHUC
                  CHV2XY   (I,J)                = CHV2
                  CHB2XY   (I,J)                = CHB2
                  PAHXY    (I,J)                = PAH
                  PAHGXY   (I,J)                = PAHG
                  PAHBXY   (I,J)                = PAHB
                  PAHVXY   (I,J)                = PAHV
                  QINTSXY  (I,J)                = QINTS
                  QINTRXY  (I,J)                = QINTR
                  QDRIPSXY (I,J)                = QDRIPS
                  QDRIPRXY (I,J)                = QDRIPR
                  QTHROSXY (I,J)                = QTHROS
                  QTHRORXY (I,J)                = QTHROR
                  QSNSUBXY (I,J)                = QSNSUB
                  QSNFROXY (I,J)                = QSNFRO
                  QSUBCXY  (I,J)                = QSUBC
                  QFROCXY  (I,J)                = QFROC
                  QEVACXY  (I,J)                = QEVAC
                  QDEWCXY  (I,J)                = QDEWC
                  QFRZCXY  (I,J)                = QFRZC
                  QMELTCXY (I,J)                = QMELTC
                  QSNBOTXY (I,J)                = QSNBOT
                  QMELTXY  (I,J)                = QMELT
                  PONDINGXY(I,J)                = PONDING + PONDING1 + PONDING2
                  FPICEXY  (I,J)                = FPICE
                  RAINLSM  (I,J)                = RAININ
                  SNOWLSM  (I,J)                = SNOWIN
                  FORCTLSM (I,J)                = T_ML
                  FORCQLSM (I,J)                = Q_ML
                  FORCPLSM (I,J)                = P_ML
                  FORCZLSM (I,J)                = Z_ML
                  FORCWLSM (I,J)                = SQRT(U_ML*U_ML + V_ML*V_ML)
                  RECHXY   (I,J)                = RECHXY(I,J) + RECH*1.E3 !RECHARGE TO THE WATER TABLE
                  DEEPRECHXY(I,J)               = DEEPRECHXY(I,J) + DEEPRECH
                  SMCWTDXY(I,J)                 = SMCWTD
                  ACC_SSOILXY (I,J)             = ACC_SSOIL
                  ACC_QINSURXY(I,J)             = ACC_QINSUR
                  ACC_QSEVAXY (I,J)             = ACC_QSEVA
                  ACC_ETRANIXY(I,:,J)           = ACC_ETRANI
                  ACC_DWATERXY(I,J)             = ACC_DWATER
                  ACC_PRCPXY  (I,J)             = ACC_PRCP
                  ACC_ECANXY  (I,J)             = ACC_ECAN
                  ACC_ETRANXY (I,J)             = ACC_ETRAN
                  ACC_EDIRXY  (I,J)             = ACC_EDIR
                  EFLXBXY (I,J)                 = EFLXB
                  SNOWENERGY(I,J)               = 0.0
                  SOILENERGY(I,J)               = 0.0
                  DO K = ISNOW+1, NSOIL
                    IF(K == ISNOW+1) THEN
                      DZSNSO = - ZSNSO(K)
                    ELSE
                      DZSNSO = ZSNSO(K-1) - ZSNSO(K)
                    END IF
                    IF(K >= 1) THEN
                      SOILENERGY(I,J) = SOILENERGY(I,J) + DZSNSO * HCPCT(K) * (STC(K)-273.16) * 0.001
                    ELSE
                      SNOWENERGY(I,J) = SNOWENERGY(I,J) + DZSNSO * HCPCT(K) * (STC(K)-273.16) * 0.001
                    END IF
                  ENDDO

                  GRAINXY  (I,J) = GRAIN !GRAIN XING
                  GDDXY    (I,J) = GDD   !XING
     	     PGSXY    (I,J) = PGS

                  ! irrigation
                  IRNUMSI(I,J)                  = IRCNTSI
                  IRNUMMI(I,J)                  = IRCNTMI
                  IRNUMFI(I,J)                  = IRCNTFI
                  IRWATSI(I,J)                  = IRAMTSI
                  IRWATMI(I,J)                  = IRAMTMI
                  IRWATFI(I,J)                  = IRAMTFI
                  IRSIVOL(I,J)                  = IRSIVOL(I,J)+(IRSIRATE*1000.0)
                  IRMIVOL(I,J)                  = IRMIVOL(I,J)+(IRMIRATE*1000.0)
                  IRFIVOL(I,J)                  = IRFIVOL(I,J)+(IRFIRATE*1000.0)
                  IRELOSS(I,J)                  = IRELOSS(I,J)+(EIRR*DT) ! mm
                  IRRSPLH(I,J)                  = IRRSPLH(I,J)+(FIRR*DT) ! Joules/m^2

                  if(iopt_crop == 2) then   ! gecros crop model

                    !*** Check for harvest
                    if ((gecros1d(1) >= gecros_ds1).and.(gecros1d(42) < 0)) then
                      if (checkIfHarvest(gecros1d, DT, gecros_ds1, gecros_ds2, gecros_ds1x, &
                          gecros_ds2x) == 1) then

                        call gecros_reinit(gecros1d)
                      endif
                    endif

                    gecros_state (i,1:60,j)     = gecros1d(1:60)
                  end if

          !-----------------------------------------------------------------------
          !BEGIN URBAN MODEL
          !---------------------------------------------------------------------

          IF (SF_URBAN_PHYSICS == 1 ) THEN         ! Beginning of UCM CALL if block

          !--------------------------------------
          ! URBAN CANOPY MODEL START
          !--------------------------------------
           IF( IVGTYP(I,J) == ISURBAN_TABLE    .or. IVGTYP(I,J) == LCZ_1_TABLE .or. IVGTYP(I,J) == LCZ_2_TABLE .or. &
                     IVGTYP(I,J) == LCZ_3_TABLE      .or. IVGTYP(I,J) == LCZ_4_TABLE .or. IVGTYP(I,J) == LCZ_5_TABLE .or. &
                     IVGTYP(I,J) == LCZ_6_TABLE      .or. IVGTYP(I,J) == LCZ_7_TABLE .or. IVGTYP(I,J) == LCZ_8_TABLE .or. &
                     IVGTYP(I,J) == LCZ_9_TABLE      .or. IVGTYP(I,J) == LCZ_10_TABLE .or. IVGTYP(I,J) == LCZ_11_TABLE ) THEN
 
      
                !  UTYPE_URB = UTYPE_URB2D(I,J) !urban type (low, high or industrial)
                !  this need to be changed in the mosaic danli
                IF (use_wudapt_lcz == 1) THEN
                  IF(IVGTYP(I,J)==ISURBAN_TABLE) UTYPE_URB=5
                  IF(IVGTYP(I,J)==LCZ_1_TABLE) UTYPE_URB=1
                  IF(IVGTYP(I,J)==LCZ_2_TABLE) UTYPE_URB=2
                  IF(IVGTYP(I,J)==LCZ_3_TABLE) UTYPE_URB=3
                  IF(IVGTYP(I,J)==LCZ_4_TABLE) UTYPE_URB=4
                  IF(IVGTYP(I,J)==LCZ_5_TABLE) UTYPE_URB=5
                  IF(IVGTYP(I,J)==LCZ_6_TABLE) UTYPE_URB=6
                  IF(IVGTYP(I,J)==LCZ_7_TABLE) UTYPE_URB=7
                  IF(IVGTYP(I,J)==LCZ_8_TABLE) UTYPE_URB=8
                  IF(IVGTYP(I,J)==LCZ_9_TABLE) UTYPE_URB=9
                  IF(IVGTYP(I,J)==LCZ_10_TABLE) UTYPE_URB=10
                  IF(IVGTYP(I,J)==LCZ_11_TABLE) UTYPE_URB=11

                  IF(UTYPE_URB==1) FRC_URB2D(I,J)=1.
                  IF(UTYPE_URB==2) FRC_URB2D(I,J)=0.99
                  IF(UTYPE_URB==3) FRC_URB2D(I,J)=1.00
                  IF(UTYPE_URB==4) FRC_URB2D(I,J)=0.65
                  IF(UTYPE_URB==5) FRC_URB2D(I,J)=0.7
                  IF(UTYPE_URB==6) FRC_URB2D(I,J)=0.65
                  IF(UTYPE_URB==7) FRC_URB2D(I,J)=0.3
                  IF(UTYPE_URB==8) FRC_URB2D(I,J)=0.85
                  IF(UTYPE_URB==9) FRC_URB2D(I,J)=0.3
                  IF(UTYPE_URB==10) FRC_URB2D(I,J)=0.55
                  IF(UTYPE_URB==11) FRC_URB2D(I,J)=1.
                ELSE
                  IF(IVGTYP(I,J)==ISURBAN_TABLE) UTYPE_URB=2
                  IF(IVGTYP(I,J)==LCZ_1_TABLE) UTYPE_URB=1 ! LOW_DENSITY_RESIDENTIAL
                  IF(IVGTYP(I,J)==LCZ_2_TABLE) UTYPE_URB=2 ! HIGH_DENSITY_RESIDENTIAL
                  IF(IVGTYP(I,J)==LCZ_3_TABLE) UTYPE_URB=3 ! HIGH_INTENSITY_INDUSTRIAL

                  IF(UTYPE_URB==1) FRC_URB2D(I,J)=0.5
                  IF(UTYPE_URB==2) FRC_URB2D(I,J)=0.9
                  IF(UTYPE_URB==3) FRC_URB2D(I,J)=0.95
                END IF



              ! Begin addition by Aaron A. on 27 May 2022.
              ! This reads in the mosaic values to hopefully
              ! cut down on interations

              DRELR_URB2D(I,J) = DRELR_URB2D_mosaic(I,mosaic_i,J)
              DRELB_URB2D(I,J) = DRELB_URB2D_mosaic(I,mosaic_i,J)
              DRELG_URB2D(I,J) = DRELG_URB2D_mosaic(I,mosaic_i,J)
              FLXHUMR_URB2D(I,J) = FLXHUMR_URB2D_mosaic(I,mosaic_i,J)
              FLXHUMB_URB2D(I,J) = FLXHUMB_URB2D_mosaic(I,mosaic_i,J)
              FLXHUMG_URB2D(I,J) = FLXHUMG_URB2D_mosaic(I,mosaic_i,J)

              CMCR_URB2D(I,J) = CMCR_URB2D_mosaic(I,mosaic_i,J)
              TGR_URB2D(I,J) = TGR_URB2D_mosaic(I,mosaic_i,J)

              DO LAYER=1,NSOIL

                  TGRL_URB3D(I,LAYER,J) = TGRL_URB3D_mosaic(I,NSOIL*(mosaic_i-1)+LAYER,J)
                  SMR_URB3D(I,LAYER,J) = SMR_URB3D_mosaic(I,NSOIL*(mosaic_i-1)+LAYER,J)

                  TRL_URB3D(I,LAYER,J) = TRL_URB3D_mosaic(I,NSOIL*(mosaic_i-1)+LAYER,J)
                  TBL_URB3D(I,LAYER,J) = TBL_URB3D_mosaic(I,NSOIL*(mosaic_i-1)+LAYER,J)
                  TGL_URB3D(I,LAYER,J) = TGL_URB3D_mosaic(I,NSOIL*(mosaic_i-1)+LAYER,J)
              ENDDO


         ! state variable surface_driver <--> lsm <--> urban

              TR_URB2D(I,J) = TR_URB2D_mosaic(I,mosaic_i,J)
              TB_URB2D(I,J) = TB_URB2D_mosaic(I,mosaic_i,J)
              TG_URB2D(I,J) = TG_URB2D_mosaic(I,mosaic_i,J)
              TC_URB2D(I,J) = TC_URB2D_mosaic(I,mosaic_i,J)
              QC_URB2D(I,J) = QC_URB2D_mosaic(I,mosaic_i,J)
              UC_URB2D(I,J) = UC_URB2D_mosaic(I,mosaic_i,J)
              XXXR_URB2D(I,J) = XXXR_URB2D_mosaic(I,mosaic_i,J)
              XXXB_URB2D(I,J) = XXXB_URB2D_mosaic(I,mosaic_i,J)
              XXXG_URB2D(I,J) = XXXG_URB2D_mosaic(I,mosaic_i,J)
              XXXC_URB2D(I,J) = XXXC_URB2D_mosaic(I,mosaic_i,J)
              SH_URB2D(I,J) = SH_URB2D_mosaic(I,mosaic_i,J)
              LH_URB2D(I,J) = LH_URB2D_mosaic(I,mosaic_i,J)
              G_URB2D(I,J) = G_URB2D_mosaic(I,mosaic_i,J)
              RN_URB2D(I,J) = RN_URB2D_mosaic(I,mosaic_i,J)
              TS_URB2D(I,J) = TS_URB2D_mosaic(I,mosaic_i,J)



              !UTYPE_URB = UTYPE_URB2D(I,J) !urban type (low, high or industrial)

              TA_URB    = T3D(I,1,J)                                ! [K]
              QA_URB    = QV3D(I,1,J)/(1.0+QV3D(I,1,J))             ! [kg/kg]
              UA_URB    = SQRT(U_PHY(I,1,J)**2.+V_PHY(I,1,J)**2.)
              U1_URB    = U_PHY(I,1,J)
              V1_URB    = V_PHY(I,1,J)
              IF(UA_URB < 1.) UA_URB=1.                             ! [m/s]
              SSG_URB   = SWDOWN(I,J)                               ! [W/m/m]
              SSGD_URB  = 0.8*SWDOWN(I,J)                           ! [W/m/m]
              SSGQ_URB  = SSG_URB-SSGD_URB                          ! [W/m/m]
              LLG_URB   = GLW(I,J)                                  ! [W/m/m]
              RAIN_URB  = PRECIP_IN(I,J) / DT * 3600.0                 ! [mm/hr]
              RHOO_URB  = (P8W3D(I,KTS+1,J)+P8W3D(I,KTS,J))*0.5 / (287.04 * TA_URB * (1.0+ 0.61 * QA_URB)) ![kg/m/m/m]
              ZA_URB    = 0.5*DZ8W(I,1,J)                           ! [m]
              DELT_URB  = DT                                        ! [sec]
              XLAT_URB  = XLAT(I,J)                           ! [deg]
              COSZ_URB  = COSZIN(I,J)
              OMG_URB   = OMG_URB2D(I,J)
              ZNT_URB   = ZNT(I,J)

              LSOLAR_URB = .FALSE.

              TR_URB = TR_URB2D(I,J)
              TB_URB = TB_URB2D(I,J)
              TG_URB = TG_URB2D(I,J)
              TC_URB = TC_URB2D(I,J)
              QC_URB = QC_URB2D(I,J)
              UC_URB = UC_URB2D(I,J)

              TGR_URB     = TGR_URB2D(I,J)
              CMCR_URB    = CMCR_URB2D(I,J)
              FLXHUMR_URB = FLXHUMR_URB2D(I,J)
              FLXHUMB_URB = FLXHUMB_URB2D(I,J)
              FLXHUMG_URB = FLXHUMG_URB2D(I,J)
              DRELR_URB   = DRELR_URB2D(I,J)
              DRELB_URB   = DRELB_URB2D(I,J)
              DRELG_URB   = DRELG_URB2D(I,J)

              DO K = 1,num_roof_layers
                TRL_URB(K) = TRL_URB3D(I,K,J)
                SMR_URB(K) = SMR_URB3D(I,K,J)
                TGRL_URB(K)= TGRL_URB3D(I,K,J)
              END DO

              DO K = 1,num_wall_layers
                TBL_URB(K) = TBL_URB3D(I,K,J)
              END DO

              DO K = 1,num_road_layers
                TGL_URB(K) = TGL_URB3D(I,K,J)
              END DO

              XXXR_URB = XXXR_URB2D(I,J)
              XXXB_URB = XXXB_URB2D(I,J)
              XXXG_URB = XXXG_URB2D(I,J)
              XXXC_URB = XXXC_URB2D(I,J)

          ! Limits to avoid dividing by small number
              IF (CHS(I,J) < 1.0E-02) THEN
                CHS(I,J)  = 1.0E-02
              ENDIF
              IF (CHS2(I,J) < 1.0E-02) THEN
                CHS2(I,J)  = 1.0E-02
              ENDIF
              IF (CQS2(I,J) < 1.0E-02) THEN
                CQS2(I,J)  = 1.0E-02
              ENDIF

              CHS_URB  = CHS(I,J)
              CHS2(I,J)= CQS2(I,J)
              CHS2_URB = CHS2(I,J)
              IF (PRESENT(CMR_SFCDIF)) THEN
                CMR_URB = CMR_SFCDIF(I,J)
                CHR_URB = CHR_SFCDIF(I,J)
                CMGR_URB = CMGR_SFCDIF(I,J)
                CHGR_URB = CHGR_SFCDIF(I,J)
                CMC_URB = CMC_SFCDIF(I,J)
                CHC_URB = CHC_SFCDIF(I,J)
              ENDIF

          ! NUDAPT for SLUCM

              MH_URB   = MH_URB2D(I,J)
              STDH_URB = STDH_URB2D(I,J)
              LP_URB   = LP_URB2D(I,J)
              HGT_URB  = HGT_URB2D(I,J)
              LF_URB   = 0.0
              DO K = 1,4
                LF_URB(K) = LF_URB2D(I,K,J)
              ENDDO
              FRC_URB  = FRC_URB2D(I,J)
              LB_URB   = LB_URB2D(I,J)
              CHECK    = 0
              IF (I.EQ.73.AND.J.EQ.125)THEN
                CHECK = 1
              END IF

          ! Call urban
              CALL cal_mon_day(INT(julian),julyr,jmonth,jday)
              CALL urban(LSOLAR_URB,                                                             & ! I
                    num_roof_layers, num_wall_layers, num_road_layers,                           & ! C
                          DZR,        DZB,        DZG, & ! C
                    UTYPE_URB,     TA_URB,     QA_URB,     UA_URB,   U1_URB,  V1_URB, SSG_URB,   & ! I
                     SSGD_URB,   SSGQ_URB,    LLG_URB,   RAIN_URB, RHOO_URB,                     & ! I
                       ZA_URB, DECLIN_URB,   COSZ_URB,    OMG_URB,                               & ! I
                     XLAT_URB,   DELT_URB,    ZNT_URB,                                           & ! I
                      CHS_URB,   CHS2_URB,                                                       & ! I
                       TR_URB,     TB_URB,     TG_URB,     TC_URB,   QC_URB,   UC_URB,           & ! H
                      TRL_URB,    TBL_URB,    TGL_URB,                                           & ! H
                     XXXR_URB,   XXXB_URB,   XXXG_URB,   XXXC_URB,                               & ! H
                       TS_URB,     QS_URB,     SH_URB,     LH_URB, LH_KINEMATIC_URB,             & ! O
                       SW_URB,    ALB_URB,     LW_URB,      G_URB,   RN_URB, PSIM_URB, PSIH_URB, & ! O
                   GZ1OZ0_URB,                                                                   & !O
                      CMR_URB,    CHR_URB,    CMC_URB,    CHC_URB,                               &
                      U10_URB,    V10_URB,    TH2_URB,     Q2_URB,                               & ! O
                      UST_URB,     mh_urb,   stdh_urb,     lf_urb,   lp_urb,                     & ! 0
                      hgt_urb,    frc_urb,     lb_urb,      check, CMCR_URB,TGR_URB,             & ! H
                     TGRL_URB,    SMR_URB,   CMGR_URB,   CHGR_URB,   jmonth,                     & ! H
                    DRELR_URB,  DRELB_URB,                                                       & ! H
                    DRELG_URB,FLXHUMR_URB,FLXHUMB_URB,FLXHUMG_URB )

              TS_URB2D(I,J) = TS_URB

              ALBEDO(I,J)   = FRC_URB2D(I,J) * ALB_URB + (1-FRC_URB2D(I,J)) * ALBEDO(I,J)        ![-]
              HFX(I,J)      = FRC_URB2D(I,J) * SH_URB  + (1-FRC_URB2D(I,J)) * HFX(I,J)           ![W/m/m]
              QFX(I,J)      = FRC_URB2D(I,J) * LH_KINEMATIC_URB &
                                 + (1-FRC_URB2D(I,J))* QFX(I,J)                                  ![kg/m/m/s]
              LH(I,J)       = FRC_URB2D(I,J) * LH_URB  + (1-FRC_URB2D(I,J)) * LH(I,J)            ![W/m/m]
              GRDFLX(I,J)   = FRC_URB2D(I,J) * (G_URB) + (1-FRC_URB2D(I,J)) * GRDFLX(I,J)        ![W/m/m]
              TSK(I,J)      = FRC_URB2D(I,J) * TS_URB  + (1-FRC_URB2D(I,J)) * TSK(I,J)           ![K]
          !    Q1            = QSFC(I,J)/(1.0+QSFC(I,J))
          !    Q1            = FRC_URB2D(I,J) * QS_URB  + (1-FRC_URB2D(I,J)) * Q1                 ![-]

          ! Convert QSFC back to mixing ratio

          !    QSFC(I,J)     = Q1/(1.0-Q1)
                             QSFC(I,J)= FRC_URB2D(I,J)*QS_URB+(1-FRC_URB2D(I,J))*QSFC(I,J)               !!   QSFC(I,J)=QSFC1D
              UST(I,J)      = FRC_URB2D(I,J) * UST_URB + (1-FRC_URB2D(I,J)) * UST(I,J)     ![m/s]
              ZNT(I,J)= EXP(FRC_URB2D(I,J)*ALOG(ZNT_URB)+(1-FRC_URB2D(I,J))* ALOG(ZNT(I,J)))   ! ADD BY DAN
          ! Renew Urban State Variables

              TR_URB2D(I,J) = TR_URB
              TB_URB2D(I,J) = TB_URB
              TG_URB2D(I,J) = TG_URB
              TC_URB2D(I,J) = TC_URB
              QC_URB2D(I,J) = QC_URB
              UC_URB2D(I,J) = UC_URB

              TGR_URB2D(I,J)     = TGR_URB
              CMCR_URB2D(I,J)    = CMCR_URB
              FLXHUMR_URB2D(I,J) = FLXHUMR_URB
              FLXHUMB_URB2D(I,J) = FLXHUMB_URB
              FLXHUMG_URB2D(I,J) = FLXHUMG_URB
              DRELR_URB2D(I,J)   = DRELR_URB
              DRELB_URB2D(I,J)   = DRELB_URB
              DRELG_URB2D(I,J)   = DRELG_URB

              DO K = 1,num_roof_layers
                TRL_URB3D(I,K,J) = TRL_URB(K)
                SMR_URB3D(I,K,J) = SMR_URB(K)
                TGRL_URB3D(I,K,J)= TGRL_URB(K)
              END DO
              DO K = 1,num_wall_layers
                TBL_URB3D(I,K,J) = TBL_URB(K)
              END DO
              DO K = 1,num_road_layers
                TGL_URB3D(I,K,J) = TGL_URB(K)
              END DO

              XXXR_URB2D(I,J)    = XXXR_URB
              XXXB_URB2D(I,J)    = XXXB_URB
              XXXG_URB2D(I,J)    = XXXG_URB
              XXXC_URB2D(I,J)    = XXXC_URB

              SH_URB2D(I,J)      = SH_URB
              LH_URB2D(I,J)      = LH_URB
              G_URB2D(I,J)       = G_URB
              RN_URB2D(I,J)      = RN_URB
              PSIM_URB2D(I,J)    = PSIM_URB
              PSIH_URB2D(I,J)    = PSIH_URB
              GZ1OZ0_URB2D(I,J)  = GZ1OZ0_URB
              U10_URB2D(I,J)     = U10_URB
              V10_URB2D(I,J)     = V10_URB
              TH2_URB2D(I,J)     = TH2_URB
              Q2_URB2D(I,J)      = Q2_URB
              UST_URB2D(I,J)     = UST_URB
              AKMS_URB2D(I,J)    = KARMAN * UST_URB2D(I,J)/(GZ1OZ0_URB2D(I,J)-PSIM_URB2D(I,J))
              IF (PRESENT(CMR_SFCDIF)) THEN
                CMR_SFCDIF(I,J)  = CMR_URB
                CHR_SFCDIF(I,J)  = CHR_URB
                CMGR_SFCDIF(I,J) = CMGR_URB
                CHGR_SFCDIF(I,J) = CHGR_URB
                CMC_SFCDIF(I,J)  = CMC_URB
                CHC_SFCDIF(I,J)  = CHC_URB
              ENDIF

            ! We also need to re-new all of the mosaic variables:
            DRELR_URB2D_mosaic(I,mosaic_i,J) = DRELR_URB2D(I,J)
            DRELB_URB2D_mosaic(I,mosaic_i,J) = DRELB_URB2D(I,J)
            DRELG_URB2D_mosaic(I,mosaic_i,J) = DRELG_URB2D(I,J)
            FLXHUMR_URB2D_mosaic(I,mosaic_i,J) = FLXHUMR_URB2D(I,J)
            FLXHUMB_URB2D_mosaic(I,mosaic_i,J) = FLXHUMR_URB2D(I,J)
            FLXHUMG_URB2D_mosaic(I,mosaic_i,J) = FLXHUMR_URB2D(I,J)

            CMCR_URB2D_mosaic(I,mosaic_i,J) = CMCR_URB2D(I,J)
            TGR_URB2D_mosaic(I,mosaic_i,J) = TGR_URB2D(I,J)

            DO LAYER=1,NSOIL

                TGRL_URB3D_mosaic(I,NSOIL*(mosaic_i-1)+LAYER,J) = TGRL_URB3D(I,LAYER,J)
                SMR_URB3D_mosaic(I,NSOIL*(mosaic_i-1)+LAYER,J) = SMR_URB3D(I,LAYER,J)

                TRL_URB3D_mosaic(I,NSOIL*(mosaic_i-1)+LAYER,J) = TRL_URB3D(I,LAYER,J)
                TBL_URB3D_mosaic(I,NSOIL*(mosaic_i-1)+LAYER,J) = TBL_URB3D(I,LAYER,J)
                TGL_URB3D_mosaic(I,NSOIL*(mosaic_i-1)+LAYER,J) = TGL_URB3D(I,LAYER,J)
            ENDDO


       ! state variable surface_driver <--> lsm <--> urban

            TR_URB2D_mosaic(I,mosaic_i,J) = TR_URB2D(I,J)
            TB_URB2D_mosaic(I,mosaic_i,J) = TB_URB2D(I,J)
            TG_URB2D_mosaic(I,mosaic_i,J) = TB_URB2D(I,J)
            TC_URB2D_mosaic(I,mosaic_i,J) = TC_URB2D(I,J)
            QC_URB2D_mosaic(I,mosaic_i,J) = QC_URB2D(I,J)
            UC_URB2D_mosaic(I,mosaic_i,J) = UC_URB2D(I,J)
            XXXR_URB2D_mosaic(I,mosaic_i,J) = XXXR_URB2D(I,J)
            XXXB_URB2D_mosaic(I,mosaic_i,J) = XXXB_URB2D(I,J)
            XXXG_URB2D_mosaic(I,mosaic_i,J) = XXXG_URB2D(I,J)
            XXXC_URB2D_mosaic(I,mosaic_i,J) = XXXC_URB2D(I,J)

            SH_URB2D_mosaic(I,mosaic_i,J) = SH_URB2D(I,J)
            LH_URB2D_mosaic(I,mosaic_i,J) = LH_URB2D(I,J)
            G_URB2D_mosaic(I,mosaic_i,J) = G_URB2D(I,J)
            RN_URB2D_mosaic(I,mosaic_i,J) = RN_URB2D(I,J)
            TS_URB2D_mosaic(I,mosaic_i,J) = TS_URB2D(I,J)

            ENDIF                                 ! urban land used type block



          ENDIF                                   ! sf_urban_physics = 1 block

          !--------------------------------------
          ! URBAN CANOPY MODEL END
          !--------------------------------------

          !-----------------------------------------------------------------------
          ! Done with the NOAH-UCM MOSAIC based off of DAN LI Aaron A.
          !-----------------------------------------------------------------------
              !we are now moving back from 2D data to 3D data

              !IN/OUT with generic LSM equivielants
              TSK_mosaic(I,mosaic_i,J) = TSK(I,J)
              HFX_mosaic(I,mosaic_i,J) = HFX(I,J)
              QFX_mosaic(I,mosaic_i,J) = QFX(I,J)
              LH_mosaic(I,mosaic_i,J) = LH(I,J)
              GRDFLX_mosaic(I,mosaic_i,J) = GRDFLX(I,J)
              SFCRUNOFF_mosaic(I,mosaic_i,J) = SFCRUNOFF(I,J)
              UDRUNOFF_mosaic(I,mosaic_i,J) = UDRUNOFF(I,J)
              ALBEDO_mosaic(I,mosaic_i,J) = ALBEDO(I,J)
              SNOWC_mosaic(I,mosaic_i,J) = SNOWC(I,J)
              SNOW_mosaic(I,mosaic_i,J) = SNOW(I,J)
              SNOWH_mosaic(I,mosaic_i,J) = SNOWH(I,J)
              CANWAT_mosaic(I,mosaic_i,J) = CANWAT(I,J)
              ACSNOM_mosaic(I,mosaic_i,J) = ACSNOM(I,J)
              ACSNOW_mosaic(I,mosaic_i,J) = ACSNOW(I,J)
              EMISS_mosaic(I,mosaic_i,J) = EMISS(I,J)
              QSFC_mosaic(I,mosaic_i,J) = QSFC(I,J)
              Z0_mosaic(I,mosaic_i,J) = Z0(I,J)
              ZNT_mosaic(I,mosaic_i,J) = ZNT(I,J)
              rs_mosaic(I,mosaic_i,J) = RS(I,J)

              !These are the soil variables
                      DO LAYER=1,NSOIL

                          SMOIS_mosaic(I,NSOIL*(mosaic_i-1)+LAYER,J) = SMOIS(I,LAYER,J)
                          SH2O_mosaic(I,NSOIL*(mosaic_i-1)+LAYER,J) = SH2O(I,LAYER,J)
                          TSLB_mosaic(I,NSOIL*(mosaic_i - 1) + LAYER,J) = TSLB(I,LAYER,J)
                          SMOISEQ_mosaic(I,NSOIL*(mosaic_i - 1) + LAYER,J) = SMOISEQ(I,LAYER,J)
                          ACC_ETRANIXY_mosaic(I,NSOIL*(mosaic_i - 1)+LAYER,J) = ACC_ETRANIXY(I,LAYER,J)
                      ENDDO


              !IN/OUT with no GENERIC LSM

              isnowxy_mosaic(I,mosaic_i,J) = ISNOWXY(I,J)
              tvxy_mosaic(I,mosaic_i,J) = TVXY(I,J)
              tgxy_mosaic(I,mosaic_i,J) = TGXY(I,J)
              canicexy_mosaic(I,mosaic_i,J) = CANICEXY(I,J)
              canliqxy_mosaic(I,mosaic_i,J) = CANLIQXY(I,J)
              eahxy_mosaic(I,mosaic_i,J) = EAHXY(I,J)
              tahxy_mosaic(I,mosaic_i,J) = TAHXY(I,J)
              cmxy_mosaic(I,mosaic_i,J) = CMXY(I,J)
              chxy_mosaic(I,mosaic_i,J) = CHXY(I,J)
              fwetxy_mosaic(I,mosaic_i,J) = FWETXY(I,J)
              sneqvoxy_mosaic(I,mosaic_i,J) = SNEQVOXY(I,J)
              alboldxy_mosaic(I,mosaic_i,J) = ALBOLDXY(I,J)
              qsnowxy_mosaic(I,mosaic_i,J) = QSNOWXY(I,J)
              qrainxy_mosaic(I,mosaic_i,J) = QRAINXY(I,J)
              wslakexy_mosaic(I,mosaic_i,J) = WSLAKEXY(I,J)
              zwtxy_mosaic(I,mosaic_i,J) = ZWTXY(I,J)
              waxy_mosaic(I,mosaic_i,J) = WAXY(I,J)
              wtxy_mosaic(I,mosaic_i,J) = WTXY(I,J)
              lfmassxy_mosaic(I,mosaic_i,J) = LFMASSXY(I,J)
              rtmassxy_mosaic(I,mosaic_i,J) = RTMASSXY(I,J)
              stmassxy_mosaic(I,mosaic_i,J) = STMASSXY(I,J)
              woodxy_mosaic(I,mosaic_i,J) = WOODXY(I,J)
              stblcpxy_mosaic(I,mosaic_i,J) = STBLCPXY(I,J)
              fastcpxy_mosaic(I,mosaic_i,J) = FASTCPXY(I,J)
              xlai_mosaic(I,mosaic_i,J) = XLAIXY(I,J)
              xsaixy_mosaic(I,mosaic_i,J) = XSAIXY(I,J)
              taussxy_mosaic(I,mosaic_i,J) = TAUSSXY(I,J)
              smcwtdxy_mosaic(I,mosaic_i,J) = SMCWTDXY(I,J)
              deeprechxy_mosaic(I,mosaic_i,J) = DEEPRECHXY(I,J)
              rechxy_mosaic(I,mosaic_i,J) = RECHXY(I,J)

              !OUT ONLY

              ! t2mvxy_mosaic(I,mosaic_i,J) = T2MVXY(I,J)
              ! t2mbxy_mosaic(I,mosaic_i,J) = T2MBXY(I,J)
              ! q2mvxy_mosaic(I,mosaic_i,J) = Q2MVXY(I,J)
              ! q2mbxy_mosaic(I,mosaic_i,J) = Q2MBXY(I,J)
              ! tradxy_mosaic(I,mosaic_i,J) = TRADXY(I,J)
              ! neexy_mosaic(I,mosaic_i,J) = NEEXY(I,J)
              ! gppxy_mosaic(I,mosaic_i,J) = GPPXY(I,J)
              ! nppxy_mosaic(I,mosaic_i,J) = NPPXY(I,J)
              ! fvegxy_mosaic(I,mosaic_i,J) = FVEGXY(I,J)
              ! runsfxy_mosaic(I,mosaic_i,J) = RUNSFXY(I,J)
              ! runsbxy_mosaic(I,mosaic_i,J) = RUNSBXY(I,J)
              ! ecanxy_mosaic(I,mosaic_i,J) = ECANXY(I,J)
              ! edirxy_mosaic(I,mosaic_i,J) = EDIRXY(I,J)
              ! etranxy_mosaic(I,mosaic_i,J) = ETRANXY(I,J)
              ! fsaxy_mosaic(I,mosaic_i,J) = FSAXY(I,J)
              ! firaxy_mosaic(I,mosaic_i,J) = FIRAXY(I,J)
              ! aparxy_mosaic(I,mosaic_i,J) = APARXY(I,J)
              ! psnxy_mosaic(I,mosaic_i,J) = PSNXY(I,J)
              ! savxy_mosaic(I,mosaic_i,J) = SAVXY(I,J)
              ! sagxy_mosaic(I,mosaic_i,J) = SAGXY(I,J)
              ! rssunxy_mosaic(I,mosaic_i,J) = RSSUNXY(I,J)
              ! rsshaxy_mosaic(I,mosaic_i,J) = RSSHAXY(I,J)
              ! bgapxy_mosaic(I,mosaic_i,J) = BGAPXY(I,J)
              ! wgapxy_mosaic(I,mosaic_i,J) = WGAPXY(I,J)
              ! tgvxy_mosaic(I,mosaic_i,J) = TGVXY(I,J)
              ! tgbxy_mosaic(I,mosaic_i,J) = TGBXY(I,J)
              ! chvxy_mosaic(I,mosaic_i,J) = CHVXY(I,J)
              ! chbxy_mosaic(I,mosaic_i,J) = CHBXY(I,J)
              ! shgxy_mosaic(I,mosaic_i,J) = SHGXY(I,J)
              ! shcxy_mosaic(I,mosaic_i,J) = SHCXY(I,J)
              ! shbxy_mosaic(I,mosaic_i,J) = SHBXY(I,J)
              ! evgxy_mosaic(I,mosaic_i,J) = EVGXY(I,J)
              ! evbxy_mosaic(I,mosaic_i,J) = EVBXY(I,J)
              ! ghvxy_mosaic(I,mosaic_i,J) = GHVXY(I,J)
              ! ghbxy_mosaic(I,mosaic_i,J) = GHBXY(I,J)
              ! irgxy_mosaic(I,mosaic_i,J) = IRGXY(I,J)
              ! ircxy_mosaic(I,mosaic_i,J) = IRCXY(I,J)
              ! irbxy_mosaic(I,mosaic_i,J) = IRBXY(I,J)
              ! trxy_mosaic(I,mosaic_i,J) = TRXY(I,J)
              ! evcxy_mosaic(I,mosaic_i,J) = EVCXY(I,J)
              ! chleafxy_mosaic(I,mosaic_i,J) = CHLEAFXY(I,J)
              ! chucxy_mosaic(I,mosaic_i,J) = CHUCXY(I,J)
              ! chv2xy_mosaic(I,mosaic_i,J) = CHV2XY(I,J)
              ! chb2xy_mosaic(I,mosaic_i,J) = CHB2XY(I,J)

              !!! NEED TO ADD NEW VARIABLES

              ! PAHXY_mosaic(I,mosaic_i,J) = PAHXY(I,J)
              ! PAHGXY_mosaic(I,mosaic_i,J) = PAHGXY(I,J)
              ! PAHBXY_mosaic(I,mosaic_i,J) = PAHBXY(I,J)    ! precipitation advected heat
              ! PAHVXY_mosaic(I,mosaic_i,J) = PAHVXY(I,J)    ! precipitation advected heat
              ! QINTSXY_mosaic(I,mosaic_i,J) = QINTSXY(I,J)
              ! QINTRXY_mosaic(I,mosaic_i,J) = QINTRXY(I,J)
              ! QDRIPSXY_mosaic(I,mosaic_i,J) = QDRIPSXY(I,J)
              ! QDRIPRXY_mosaic(I,mosaic_i,J) = QDRIPRXY(I,J)
              ! QTHROSXY_mosaic(I,mosaic_i,J) = QTHROSXY(I,J)
              ! QTHRORXY_mosaic(I,mosaic_i,J) = QTHRORXY(I,J)
              ! QSNSUBXY_mosaic(I,mosaic_i,J) = QSNSUBXY(I,J)
              ! QSNFROXY_mosaic(I,mosaic_i,J) = QSNFROXY(I,J)
              ! QSUBCXY_mosaic(I,mosaic_i,J) = QSUBCXY(I,J)
              ! QFROCXY_mosaic(I,mosaic_i,J) = QFROCXY(I,J)
              ! QEVACXY_mosaic(I,mosaic_i,J) = QEVACXY(I,J)
              ! QDEWCXY_mosaic(I,mosaic_i,J) = QDEWCXY(I,J)
              ! QFRZCXY_mosaic(I,mosaic_i,J) = QFRZCXY(I,J)
              ! QMELTCXY_mosaic(I,mosaic_i,J) = QMELTCXY(I,J)
              ! QSNBOTXY_mosaic(I,mosaic_i,J) = QSNBOTXY(I,J)
              ! QMELTXY_mosaic(I,mosaic_i,J) = QMELTXY(I,J)
              ! PONDINGXY_mosaic(I,mosaic_i,J) = PONDINGXY(I,J)
              ! FPICEXY_mosaic(I,mosaic_i,J) = FPICEXY(I,J)
              ACC_SSOILXY_mosaic(I,mosaic_i,J) = ACC_SSOILXY(I,J)
              ACC_QINSURXY_mosaic(I,mosaic_i,J) = ACC_QINSURXY(I,J)
              ACC_QSEVAXY_mosaic(I,mosaic_i,J) = ACC_QSEVAXY(I,J)
              ! EFLXBXY_mosaic(I,mosaic_i,J) = EFLXBXY(I,J)
              ! SOILENERGY_mosaic(I,mosaic_i,J) = SOILENERGY(I,J)
              ! SNOWENERGY_mosaic(I,mosaic_i,J) = SNOWENERGY(I,J)
              ! CANHSXY_mosaic(I,mosaic_i,J) = CANHSXY(I,J)
              ACC_DWATERXY_mosaic(I,mosaic_i,J) = ACC_DWATERXY(I,J)
              ACC_PRCPXY_mosaic(I,mosaic_i,J) = ACC_PRCPXY(I,J)
              ACC_ECANXY_mosaic(I,mosaic_i,J) = ACC_ECANXY(I,J)
              ACC_ETRANXY_mosaic(I,mosaic_i,J) = ACC_ETRANXY(I,J)
              ACC_EDIRXY_mosaic(I,mosaic_i,J) = ACC_EDIRXY(I,J)

              ! Irrigation variables that are needed
              !2D inout irrigation variables
              IRNUMSI_mosaic(I,mosaic_i,J) = IRNUMSI(I,J)  ! irrigation event number, Sprinkler
              IRNUMMI_mosaic(I,mosaic_i,J) = IRNUMMI(I,J)   ! irrigation event number, Micro
              IRNUMFI_mosaic(I,mosaic_i,J) = IRNUMFI(I,J)   ! irrigation event number, Flood
              IRWATSI_mosaic(I,mosaic_i,J) = IRWATSI(I,J)    ! irrigation water amount [m] to be applied, Sprinkler
              IRWATMI_mosaic(I,mosaic_i,J) = IRWATMI(I,J)    ! irrigation water amount [m] to be applied, Micro
              IRWATFI_mosaic(I,mosaic_i,J) = IRWATFI(I,J)    ! irrigation water amount [m] to be applied, Flood
              IRELOSS_mosaic(I,mosaic_i,J) = IRELOSS(I,J)    ! loss of irrigation water to evaporation,sprinkler [m/timestep]
              IRSIVOL_mosaic(I,mosaic_i,J) = IRSIVOL(I,J)    ! amount of irrigation by sprinkler (mm)
              IRMIVOL_mosaic(I,mosaic_i,J) = IRMIVOL(I,J)    ! amount of irrigation by micro (mm)
              IRFIVOL_mosaic(I,mosaic_i,J) = IRFIVOL(I,J)    ! amount of irrigation by micro (mm)
              IRRSPLH_mosaic(I,mosaic_i,J) = IRRSPLH(I,J)    ! latent heating from sprinkler evaporation (w/m2)




              !SNOW VARIABLES

                  DO LAYER=1,3
                      tsnoxy_mosaic(I,3*(mosaic_i - 1) + LAYER,J) = TSNOXY(I,LAYER-3,J)
                      snicexy_mosaic(I,3*(mosaic_i - 1) + LAYER,J) = SNICEXY(I,LAYER-3,J)
                      snliqxy_mosaic(I,3*(mosaic_i - 1) + LAYER,J) = SNLIQXY(I,LAYER-3,J)
                  ENDDO

                  DO LAYER=1,7
                      zsnsoxy_mosaic(I,7*(mosaic_i - 1) + LAYER, J) = ZSNSOXY(I,LAYER-3,J)
                  ENDDO

              !HUE THINGS
              RUNONSFXY_mosaic(I,mosaic_i,J) =  RUNONSFXY_mosaic(I,mosaic_i,J) + RUNONSRF * DT
              DETENTION_STORAGEXY_mosaic(I,mosaic_i,J) = DETENTION_STORAGE
                IF (VEGTYP.eq.46) THEN
                  DETENTION_STORAGEXY(I,J) = DETENTION_STORAGE
                else
                  DETENTION_STORAGEXY(I,J) = 0.0
                ENDIF

              !-------------------------------------------------------------------
              !We now do Grid Averaging!!!!
              !-------------------------------------------------------------------

              FAREA = landusef2(I,mosaic_i,J) !This logical option is going to be used multiply
              FAREA2 = landusef2(I,mosaic_i,J)

              IF( IVGTYP(I,J) == ISURBAN_TABLE .or. IVGTYP(I,J) == LCZ_1_TABLE .or. IVGTYP(I,J) == LCZ_2_TABLE .or. &
                IVGTYP(I,J) == LCZ_3_TABLE .or. IVGTYP(I,J) == 42 .or. IVGTYP(I,J) == 44 .or. &
                IVGTYP(I,J) == 45.or. IVGTYP(I,J) == ISBARREN_TABLE ) THEN
                FAREA2 = 0.0 !This is a fix to change how the urban areas/barren areas
                ! We are only going to be outputting the average of the vegetated area.
               ENDIF


              TSK_mosaic_avg(I,J) = TSK_mosaic_avg(I,J) + (EMISS_mosaic(I,mosaic_i,J)*TSK_mosaic(I,mosaic_i,J)**4)*FAREA  ! Conserve the longwave radiation
              !TSK_mosaic_avg(I,J) = TSK_mosaic_avg(I,J) + TSK_mosaic(I,mosaic_i,J)*FAREA
              FAREA_mosaic_avg(I,J) = FAREA_mosaic_avg(I,J) + FAREA
              HFX_mosaic_avg(I,J) = HFX_mosaic_avg(I,J) + HFX_mosaic(I,mosaic_i,J)*FAREA
              QFX_mosaic_avg(I,J) = QFX_mosaic_avg(I,J) + QFX_mosaic(I,mosaic_i,J)*FAREA
              LH_mosaic_avg(I,J) =  LH_mosaic_avg(I,J)  + LH_mosaic(I,mosaic_i,J)*FAREA
              GRDFLX_mosaic_avg(I,J) = GRDFLX_mosaic_avg(I,J) + GRDFLX_mosaic(I,mosaic_i,J)*FAREA
              SFCRUNOFF_mosaic_sum(I,J) = SFCRUNOFF_mosaic_sum(I,J) + SFCRUNOFF_mosaic(I,mosaic_i,J)*FAREA
              UDRUNOFF_mosaic_sum(I,J) = UDRUNOFF_mosaic_sum(I,J) + UDRUNOFF_mosaic(I,mosaic_i,J)*FAREA
              ALBEDO_mosaic_avg(I,J) = ALBEDO_mosaic_avg(I,J) + ALBEDO_mosaic(I,mosaic_i,J)*FAREA
              SNOWC_mosaic_avg(I,J) = SNOWC_mosaic_avg(I,J) + SNOWC_mosaic(I,mosaic_i,J)*FAREA
              CANWAT_mosaic_avg(I,J) = CANWAT_mosaic_avg(I,J) + CANWAT_mosaic(I,mosaic_i,J)*FAREA
              SNOW_mosaic_avg(I,J) = SNOW_mosaic_avg(I,J) + SNOW_mosaic(I,mosaic_i,J)*FAREA
              SNOWH_mosaic_avg(I,J) = SNOWH_mosaic_avg(I,J) + SNOWH_mosaic(I,mosaic_i,J)*FAREA
              ACSNOM_mosaic_avg(I,J) = ACSNOM_mosaic_avg(I,J) + ACSNOM_mosaic(I,mosaic_i,J)*FAREA
              ACSNOW_mosaic_avg(I,J) = ACSNOW_mosaic_avg(I,J) +  ACSNOW_mosaic(I,mosaic_i,J)*FAREA
              EMISS_mosaic_avg(I,J) = EMISS_mosaic_avg(I,J) + EMISS_mosaic(I,mosaic_i,J)*FAREA
              QSFC_mosaic_avg(I,J) = QSFC_mosaic_avg(I,J) + QSFC_mosaic(I,mosaic_i,J)*FAREA
              Z0_mosaic_avg(I,J) = Z0_mosaic_avg(i,j)+ALOG(Z0_mosaic(i,mosaic_i,j))*FAREA
              ZNT_mosaic_avg(I,J) = ZNT_mosaic_avg(i,j)+ALOG(ZNT_mosaic(i,mosaic_i,j))*FAREA
              rs_mosaic_avg(I,J) = rs_mosaic_avg(I,J) + rs_mosaic(I,mosaic_i,J)*FAREA

              DO LAYER=1,NSOIL

                  TSLB_mosaic_avg(I,LAYER,J) = TSLB_mosaic_avg(I,LAYER,J) + TSLB_mosaic(I,NSOIL*(mosaic_i-1)+LAYER,J)*FAREA
                  SMOIS_mosaic_avg(I,LAYER,J) = SMOIS_mosaic_avg(I,LAYER,J) + SMOIS_mosaic(I,NSOIL*(mosaic_i-1)+LAYER,J)*FAREA
                  SH2O_mosaic_avg(I,LAYER,J) = SH2O_mosaic_avg(I,LAYER,J) + SH2O_mosaic(I,NSOIL*(mosaic_i-1)+LAYER,J)*FAREA
                  SMOISEQ_mosaic_avg(I,LAYER,J) = SMOISEQ_mosaic_avg(I,LAYER,J) + SMOISEQ_mosaic(I,NSOIL*(mosaic_i-1)+LAYER,J)*FAREA
                  ACC_ETRANIXY_mosaic_avg(I,LAYER,J) = ACC_ETRANIXY_mosaic_avg(I,LAYER,J) + ACC_ETRANIXY_mosaic(I,NSOIL*(mosaic_i - 1)+LAYER,J)*FAREA2
              ENDDO

              !IN OUT LSM NO LSM EQUIVELANTS
              isnowxy_mosaic_avg(I,J) = isnowxy_mosaic_avg(I,J) + isnowxy_mosaic(I,mosaic_i,J) !not normalized by the area because it is an integer
              tvxy_mosaic_avg(I,J) = tvxy_mosaic_avg(I,J) + tvxy_mosaic(I,mosaic_i,J)*FAREA2
              tgxy_mosaic_avg(I,J) = tgxy_mosaic_avg(I,J) + tgxy_mosaic(I,mosaic_i,J)*FAREA
              canicexy_mosaic_avg(I,J) = canicexy_mosaic_avg(I,J) + canicexy_mosaic(I,mosaic_i,J)*FAREA
              canliqxy_mosaic_avg(I,J) = canliqxy_mosaic_avg(I,J) + canliqxy_mosaic(I,mosaic_i,J)*FAREA
              eahxy_mosaic_avg(I,J) = eahxy_mosaic_avg(I,J) + eahxy_mosaic(I,mosaic_i,J)*FAREA2
              tahxy_mosaic_avg(I,J) = tahxy_mosaic_avg(I,J) + tahxy_mosaic(I,mosaic_i,J)*FAREA2
              cmxy_mosaic_avg(I,J) = cmxy_mosaic_avg(I,J) + cmxy_mosaic(I,mosaic_i,J)*FAREA
              chxy_mosaic_avg(I,J) = chxy_mosaic_avg(I,J) + chxy_mosaic(I,mosaic_i,J)*FAREA
              fwetxy_mosaic_avg(I,J) = fwetxy_mosaic_avg(I,J) + fwetxy_mosaic(I,mosaic_i,J)*FAREA2
              sneqvoxy_mosaic_avg(I,J) = sneqvoxy_mosaic_avg(I,J) + sneqvoxy_mosaic(I,mosaic_i,J)*FAREA
              alboldxy_mosaic_avg(I,J) = alboldxy_mosaic_avg(I,J) + alboldxy_mosaic(I,mosaic_i,J)*FAREA
              qsnowxy_mosaic_avg(I,J) = qsnowxy_mosaic_avg(I,J) + qsnowxy_mosaic(I,mosaic_i,J)*FAREA
              qrainxy_mosaic_avg(I,J) = qrainxy_mosaic_avg(I,J) + qrainxy_mosaic(I,mosaic_i,J)*FAREA
              wslakexy_mosaic_avg(I,J) = wslakexy_mosaic_avg(I,J) + wslakexy_mosaic(I,mosaic_i,J)*FAREA
              zwtxy_mosaic_avg(I,J) = zwtxy_mosaic_avg(I,J) + zwtxy_mosaic(I,mosaic_i,J)*FAREA
              waxy_mosaic_avg(I,J) = waxy_mosaic_avg(I,J) + waxy_mosaic(I,mosaic_i,J)*FAREA
              wtxy_mosaic_avg(I,J) = wtxy_mosaic_avg(I,J) + wtxy_mosaic(I,mosaic_i,J)*FAREA
              lfmassxy_mosaic_avg(I,J) = lfmassxy_mosaic_avg(I,J) + lfmassxy_mosaic(I,mosaic_i,J)*FAREA2
              rtmassxy_mosaic_avg(I,J) = rtmassxy_mosaic_avg(I,J) + rtmassxy_mosaic(I,mosaic_i,J)*FAREA2
              stmassxy_mosaic_avg(I,J) = stmassxy_mosaic_avg(I,J) + stmassxy_mosaic(I,mosaic_i,J)*FAREA2
              woodxy_mosaic_avg(I,J) = woodxy_mosaic_avg(I,J) + woodxy_mosaic(I,mosaic_i,J)*FAREA2
              stblcpxy_mosaic_avg(I,J) = stblcpxy_mosaic_avg(I,J) + stblcpxy_mosaic(I,mosaic_i,J)*FAREA2
              fastcpxy_mosaic_avg(I,J) = fastcpxy_mosaic_avg(I,J) + fastcpxy_mosaic(I,mosaic_i,J)*FAREA2
              xsaixy_mosaic_avg(I,J) = xsaixy_mosaic_avg(I,J) + xsaixy_mosaic(I,mosaic_i,J)*FAREA2
              xlai_mosaic_avg(I,J) = xlai_mosaic_avg(I,J) + xlai_mosaic(I,mosaic_i,J)*FAREA2
              taussxy_mosaic_avg(I,J) = taussxy_mosaic_avg(I,J) + taussxy_mosaic(I,mosaic_i,J)*FAREA
              rechxy_mosaic_avg(I,J) = rechxy_mosaic_avg(I,J) + rechxy_mosaic(I,mosaic_i,J)*FAREA
              deeprechxy_mosaic_avg(I,J) = deeprechxy_mosaic_avg(I,J) + deeprechxy_mosaic(I,mosaic_i,J)*FAREA
              smcwtdxy_mosaic_avg(I,J) = smcwtdxy_mosaic_avg(I,J) + smcwtdxy_mosaic(I,mosaic_i,J)*FAREA

              !OUT LSM EQUIVELANTS

              t2mvxy_mosaic_avg(I,J) = t2mvxy_mosaic_avg(I,J) + t2mvxy(I,J)*FAREA2
              t2mbxy_mosaic_avg(I,J) = t2mbxy_mosaic_avg(I,J) + t2mbxy(I,J)*FAREA
              q2mvxy_mosaic_avg(I,J) = q2mvxy_mosaic_avg(I,J) + q2mvxy(I,J)*FAREA2
              q2mbxy_mosaic_avg(I,J) = q2mbxy_mosaic_avg(I,J) + q2mbxy(I,J)*FAREA
              tradxy_mosaic_avg(I,J) = tradxy_mosaic_avg(I,J) + tradxy(I,J)*FAREA
              neexy_mosaic_avg(I,J) = neexy_mosaic_avg(I,J) + neexy(I,J)*FAREA2
              gppxy_mosaic_avg(I,J) = gppxy_mosaic_avg(I,J) + gppxy(I,J)*FAREA2
              nppxy_mosaic_avg(I,J) = nppxy_mosaic_avg(I,J) + nppxy(I,J)*FAREA2
              fvegxy_mosaic_avg(I,J) = fvegxy_mosaic_avg(I,J) + fvegxy(I,J)*FAREA2
              runsfxy_mosaic_avg(I,J) = runsfxy_mosaic_avg(I,J) + runsfxy(I,J)*FAREA
              runsbxy_mosaic_avg(I,J) = runsbxy_mosaic_avg(I,J) + runsbxy(I,J)*FAREA
              ecanxy_mosaic_avg(I,J) = ecanxy_mosaic_avg(I,J) + ecanxy(I,J)*FAREA2
              edirxy_mosaic_avg(I,J) = edirxy_mosaic_avg(I,J) + edirxy(I,J)*FAREA2
              etranxy_mosaic_avg(I,J) = etranxy_mosaic_avg(I,J) + etranxy(I,J)*FAREA2
              fsaxy_mosaic_avg(I,J) = fsaxy_mosaic_avg(I,J) + fsaxy(I,J)*FAREA
              firaxy_mosaic_avg(I,J) = firaxy_mosaic_avg(I,J) + firaxy(I,J)*FAREA
              aparxy_mosaic_avg(I,J) = aparxy_mosaic_avg(I,J) + aparxy(I,J)*FAREA2
              psnxy_mosaic_avg(I,J) = psnxy_mosaic_avg(I,J) + psnxy(I,J)*FAREA2
              savxy_mosaic_avg(I,J) = savxy_mosaic_avg(I,J) + savxy(I,J)*FAREA2
              sagxy_mosaic_avg(I,J) = sagxy_mosaic_avg(I,J) + sagxy(I,J)*FAREA
              rssunxy_mosaic_avg(I,J) = rssunxy_mosaic_avg(I,J) + rssunxy(I,J)*FAREA2
              rsshaxy_mosaic_avg(I,J) = rsshaxy_mosaic_avg(I,J) + rsshaxy(I,J)*FAREA2
              bgapxy_mosaic_avg(I,J) = bgapxy_mosaic_avg(I,J) + bgapxy(I,J)*FAREA2
              wgapxy_mosaic_avg(I,J) = wgapxy_mosaic_avg(I,J) + wgapxy(I,J)*FAREA2
              tgvxy_mosaic_avg(I,J) = tgvxy_mosaic_avg(I,J) + tgvxy(I,J)*FAREA2
              tgbxy_mosaic_avg(I,J) = tgbxy_mosaic_avg(I,J) + tgbxy(I,J)*FAREA
              chvxy_mosaic_avg(I,J) = chvxy_mosaic_avg(I,J) + chvxy(I,J)*FAREA2
              chbxy_mosaic_avg(I,J) = chbxy_mosaic_avg(I,J) + chbxy(I,J)*FAREA
              shgxy_mosaic_avg(I,J) = shgxy_mosaic_avg(I,J) + shgxy(I,J)*FAREA2
              shcxy_mosaic_avg(I,J) = shcxy_mosaic_avg(I,J) + shcxy(I,J)*FAREA2
              shbxy_mosaic_avg(I,J) = shbxy_mosaic_avg(I,J) + shbxy(I,J)*FAREA
              evgxy_mosaic_avg(I,J) = evgxy_mosaic_avg(I,J) + evgxy(I,J)*FAREA2
              evbxy_mosaic_avg(I,J) = evbxy_mosaic_avg(I,J) + evbxy(I,J)*FAREA
              ghvxy_mosaic_avg(I,J) = ghvxy_mosaic_avg(I,J) + ghvxy(I,J)*FAREA2
              ghbxy_mosaic_avg(I,J) = ghbxy_mosaic_avg(I,J) + ghbxy(I,J)*FAREA
              irgxy_mosaic_avg(I,J) = irgxy_mosaic_avg(I,J) + irgxy(I,J)*FAREA2
              ircxy_mosaic_avg(I,J) = ircxy_mosaic_avg(I,J) + ircxy(I,J)*FAREA2
              irbxy_mosaic_avg(I,J) = irbxy_mosaic_avg(I,J) + irbxy(I,J)*FAREA
              trxy_mosaic_avg(I,J) = trxy_mosaic_avg(I,J) + trxy(I,J)*FAREA2
              evcxy_mosaic_avg(I,J) = evcxy_mosaic_avg(I,J) + evcxy(I,J)*FAREA2
              chleafxy_mosaic_avg(I,J) = chleafxy_mosaic_avg(I,J) + chleafxy(I,J)*FAREA2
              chucxy_mosaic_avg(I,J) = chucxy_mosaic_avg(I,J) + chucxy(I,J)*FAREA2
              chv2xy_mosaic_avg(I,J) = chv2xy_mosaic_avg(I,J) + chv2xy(I,J)*FAREA2
              chb2xy_mosaic_avg(I,J) = chb2xy_mosaic_avg(I,J) + chb2xy(I,J)*FAREA

              ! Extra Variables for the averaging
              ! irrigation intermediate variables
              IRWATSI_mosaic_avg(I,J) =  IRWATSI_mosaic_avg(I,J) + IRWATSI_mosaic(I,mosaic_i,J)*FAREA
              IRWATMI_mosaic_avg(I,J) =  IRWATMI_mosaic_avg(I,J) + IRWATMI_mosaic(I,mosaic_i,J)*FAREA
              IRWATFI_mosaic_avg(I,J) =  IRWATFI_mosaic_avg(I,J) + IRWATFI_mosaic(I,mosaic_i,J)*FAREA
              IRELOSS_mosaic_avg(I,J) =  IRELOSS_mosaic_avg(I,J) + IRELOSS_mosaic(I,mosaic_i,J)*FAREA
              IRSIVOL_mosaic_avg(I,J) =  IRSIVOL_mosaic_avg(I,J) + IRSIVOL_mosaic(I,mosaic_i,J)*FAREA
              IRMIVOL_mosaic_avg(I,J) =  IRMIVOL_mosaic_avg(I,J) + IRMIVOL_mosaic(I,mosaic_i,J)*FAREA
              IRFIVOL_mosaic_avg(I,J) =  IRFIVOL_mosaic_avg(I,J) + IRFIVOL_mosaic(I,mosaic_i,J)*FAREA
              IRRSPLH_mosaic_avg(I,J) =  IRRSPLH_mosaic_avg(I,J) + IRRSPLH_mosaic(I,mosaic_i,J)*FAREA

              IRNUMSI_mosaic_avg(I,J) = IRNUMSI_mosaic_avg(I,J) + IRNUMSI_mosaic(I,mosaic_i,J)
              IRNUMMI_mosaic_avg(I,J) = IRNUMMI_mosaic_avg(I,J) + IRNUMMI_mosaic(I,mosaic_i,J)
              IRNUMFI_mosaic_avg(I,J) = IRNUMFI_mosaic_avg(I,J) + IRNUMFI_mosaic(I,mosaic_i,J)

               DO LAYER=1,3
                      tsnoxy_mosaic_avg(I, LAYER,J) = tsnoxy_mosaic_avg(I, LAYER,J) +  tsnoxy_mosaic(I, 3*(mosaic_i - 1) + LAYER,J)*FAREA
                      snicexy_mosaic_avg(I,LAYER,J) = snicexy_mosaic_avg(I,LAYER,J) +  snicexy_mosaic(I,3*(mosaic_i - 1) + LAYER,J)*FAREA
                      snliqxy_mosaic_avg(I,LAYER,J) = SNLIQXY_mosaic_avg(I,LAYER,J) + snliqxy_mosaic(I,3*(mosaic_i - 1) + LAYER,J)*FAREA
               ENDDO

               DO LAYER=1,7
                      zsnsoxy_mosaic_avg(I,LAYER,J) = zsnsoxy_mosaic_avg(I,LAYER,J) + zsnsoxy_mosaic(I,7*(mosaic_i - 1) + LAYER,J)*FAREA
               ENDDO

              ! Extra variables averaged needed for extra outputs
              QINTSXY_mosaic_avg(I,J) = QINTSXY_mosaic_avg(I,J) + QINTSXY(I,J)*FAREA
              QINTRXY_mosaic_avg(I,J) = QINTRXY_mosaic_avg(I,J) + QINTRXY(I,J)*FAREA
              QDRIPSXY_mosaic_avg(I,J) = QDRIPSXY_mosaic_avg(I,J) + QDRIPSXY(I,J)*FAREA
              QDRIPRXY_mosaic_avg(I,J) = QDRIPRXY_mosaic_avg(I,J) + QDRIPRXY(I,J)*FAREA
              QTHROSXY_mosaic_avg(I,J) = QTHROSXY_mosaic_avg(I,J) + QTHRORXY(I,J)*FAREA
              QTHRORXY_mosaic_avg(I,J) = QTHRORXY_mosaic_avg(I,J) + QTHRORXY(I,J)*FAREA
              QSNSUBXY_mosaic_avg(I,J) = QSNSUBXY_mosaic_avg(I,J) + QSNSUBXY(I,J)*FAREA
              QSNFROXY_mosaic_avg(I,J) = QSNFROXY_mosaic_avg(I,J) + QSNFROXY(I,J)*FAREA
              QSUBCXY_mosaic_avg(I,J) = QSUBCXY_mosaic_avg(I,J) + QSUBCXY(I,J)*FAREA
              QFROCXY_mosaic_avg(I,J) = QFROCXY_mosaic_avg(I,J) + QFROCXY(I,J)*FAREA
              QEVACXY_mosaic_avg(I,J) = QEVACXY_mosaic_avg(I,J) + QEVACXY(I,J)*FAREA
              QDEWCXY_mosaic_avg(I,J) = QDEWCXY_mosaic_avg(I,J) + QDEWCXY(I,J)*FAREA
              QFRZCXY_mosaic_avg(I,J) = QFRZCXY_mosaic_avg(I,J) + QFRZCXY(I,J)*FAREA
              QMELTCXY_mosaic_avg(I,J) = QMELTCXY_mosaic_avg(I,J) + QMELTCXY(I,J)*FAREA
              QSNBOTXY_mosaic_avg(I,J) = QSNBOTXY_mosaic_avg(I,J) + QSNBOTXY(I,J)*FAREA
              QMELTXY_mosaic_avg(I,J) = QMELTXY_mosaic_avg(I,J) + QMELTXY(I,J)*FAREA
              PONDINGXY_mosaic_avg(I,J) = PONDINGXY_mosaic_avg(I,J) + PONDINGXY(I,J)*FAREA
              PAHXY_mosaic_avg(I,J) = PAHXY_mosaic_avg(I,J) + PAHXY(I,J)*FAREA
              PAHVXY_mosaic_avg(I,J) = PAHVXY_mosaic_avg(I,J) + PAHVXY(I,J)*FAREA
              PAHBXY_mosaic_avg(I,J) = PAHBXY_mosaic_avg(I,J) + PAHBXY(I,J)*FAREA
              PAHGXY_mosaic_avg(I,J) = PAHGXY_mosaic_avg(I,J) + PAHGXY(I,J)*FAREA
              FPICEXY_mosaic_avg(I,J) = FPICEXY_mosaic_avg(I,J) + FPICEXY(I,J)*FAREA

              ! Soil and ACC
              ACC_SSOILXY_mosaic_avg(I,J) = ACC_SSOILXY_mosaic_avg(I,J) + ACC_SSOILXY_mosaic(I,mosaic_i,J)*FAREA
              ACC_QINSURXY_mosaic_avg(I,J) = ACC_QINSURXY_mosaic_avg(I,J) + ACC_QINSURXY_mosaic(I,mosaic_i,J)*FAREA
              ACC_QSEVAXY_mosaic_avg(I,J) = ACC_QSEVAXY_mosaic_avg(I,J) + ACC_QSEVAXY_mosaic(I,mosaic_i,J)*FAREA
              EFLXBXY_mosaic_avg(I,J) = EFLXBXY_mosaic_avg(I,J) + EFLXBXY(I,J)*FAREA
              SOILENERGY_mosaic_avg(I,J) = SOILENERGY_mosaic_avg(I,J) + SOILENERGY(I,J)*FAREA
              SNOWENERGY_mosaic_avg(I,J) = SNOWENERGY_mosaic_avg(I,J) + SNOWENERGY(I,J)*FAREA
              CANHSXY_mosaic_avg(I,J) = CANHSXY_mosaic_avg(I,J) + CANHSXY(I,J)*FAREA
              ACC_DWATERXY_mosaic_avg(I,J) = ACC_DWATERXY_mosaic_avg(I,J) + ACC_DWATERXY_mosaic(I,mosaic_i,J)*FAREA
              ACC_PRCPXY_mosaic_avg(I,J) = ACC_PRCPXY_mosaic_avg(I,J) + ACC_PRCPXY_mosaic(I,mosaic_i,J)*FAREA
              ACC_ECANXY_mosaic_avg(I,J) = ACC_ECANXY_mosaic_avg(I,J) + ACC_ECANXY_mosaic(I,mosaic_i,J)*FAREA
              ACC_EDIRXY_mosaic_avg(I,J) = ACC_EDIRXY_mosaic_avg(I,J) + ACC_EDIRXY_mosaic(I,mosaic_i,J)*FAREA
              ACC_ETRANXY_mosaic_avg(I,J) = ACC_ETRANXY_mosaic_avg(I,J) + ACC_ETRANXY_mosaic(I,mosaic_i,J)*FAREA
              !HUE THINGS
              RUNONSFXY_mosaic_avg(I,J) = RUNONSFXY_mosaic_avg(I,J) + RUNONSFXY_mosaic(I,mosaic_i,J)*FAREA

             END DO MOSAIC_LOOP  !Loop over each of the mosaic variables
              !-------------------------------------------------------------------
              !WE now send the Mosaic Values to 2D variables so they can be called in other routines
              !-------------------------------------------------------------------

              !Final Step, we adjust outputs that have only outputs
              !for vegetation. We do this because some variables are not
              !defined in other areas, and

              FAREA2 = 0.

              DO mosaic_i= 1,mosaic_cat
                IF( (mosaic_cat_index(I,mosaic_i,J) /= ISURBAN_TABLE) .and.( mosaic_cat_index(I,mosaic_i,J) /= LCZ_1_TABLE) .and. &
                   (mosaic_cat_index(I,mosaic_i,J) /= LCZ_2_TABLE) .and. (mosaic_cat_index(I,mosaic_i,J) /= LCZ_3_TABLE) .and. (mosaic_cat_index(I,mosaic_i,J) /= 42) .and. &
                   (mosaic_cat_index(I,mosaic_i,J)/= 44) .and. (IVGTYP(I,J) /= 45).and. (mosaic_cat_index(I,mosaic_i,J) /= ISBARREN_TABLE) ) THEN

                  FAREA2 = FAREA2 + landusef2(I,mosaic_i,J)
                ENDIF
              END DO
              
              IF (FAREA2 == 0.) THEN 
                  WRITE(*,*) "Grid point I = ",I, "Grid Point J = ", J ," Has no vegetated area across mosaic values. Veg outputs to -9999"
                  TRIGGER = 1
              endif
              ! FAREA is now the total area that was used
              ! to average our vegetation canopy. We divide the
              ! variables that we just added by this to scale to be correct
              ! if FAREA2 was multiplied by in the previous lines.
              ! if not, we just pass the average value (which is correct)
              ! into its final place to be written out

              IVGTYP(I,J) = IVGTYP_dominant(I,J) !dominant vegetation catagory
              TSK(I,J) = (TSK_mosaic_avg(I,J)/EMISS_mosaic_avg(I,J))**(0.25) 
              HFX(I,J) = HFX_mosaic_avg(I,J)
              QFX(I,J) = QFX_mosaic_avg(I,J)
              LH(I,J) = LH_mosaic_avg(I,J)
              GRDFLX(I,J) = GRDFLX_mosaic_avg(I,J)
              SFCRUNOFF(I,J) = SFCRUNOFF_mosaic_sum(I,J)
              UDRUNOFF(I,J) = UDRUNOFF_mosaic_sum(I,J)
              ALBEDO(I,J) = ALBEDO_mosaic_avg(I,J)
              SNOWC(I,J) = SNOWC_mosaic_avg(I,J)
              SNOW(I,J) = SNOW_mosaic_avg(I,J)
              SNOWH(I,J) = SNOWH_mosaic_avg(I,J)
              CANWAT(I,J) = CANWAT_mosaic_avg(I,J)
              ACSNOM(I,J) = ACSNOM_mosaic_avg(I,J)
              ACSNOW(I,J) = ACSNOW_mosaic_avg(I,J)
              EMISS(I,J) = EMISS_mosaic_avg(I,J)
              QSFC(I,J) = QSFC_mosaic_avg(I,J)
              Z0(I,J) =  EXP(Z0_mosaic_avg(i,j)/FAREA_mosaic_avg(i,j))
              ZNT(I,J) =  EXP(ZNT_mosaic_avg(i,j)/FAREA_mosaic_avg(i,j))
              RS(I,J) = rs_mosaic_avg(I,J)

              
              DO LAYER=1,NSOIL

                  TSLB(I,LAYER,J) = TSLB_mosaic_avg(I,LAYER,J)
                  SMOIS(I,LAYER,J) = SMOIS_mosaic_avg(I,LAYER,J)
                  SH2O(I,LAYER,J) = SH2O_mosaic_avg(I,LAYER,J)
                  SMOISEQ(I,LAYER,J) = SMOISEQ_mosaic_avg(I,LAYER,J)
                  IF (TRIGGER == 0) THEN
                      ACC_ETRANIXY(I,LAYER,J) = ACC_ETRANIXY_mosaic_avg(I,LAYER,J)/FAREA2
                  else
                       ACC_ETRANIXY(I,LAYER,J) = -9999.
                  endif
                   
              ENDDO
              

              ISNOWXY(I,J) = isnowxy_mosaic_avg(I,J)/mosaic_cat
              IF(SNOWH(I,J).ne.0.0.and.QSNOWXY(I,J).eq.0.0) THEN
              DO LAYER= 1,mosaic_cat
                 IF ((SNOWH_MOSAIC(I,LAYER,J).ne.0.0).and.(QSNOWXY_mosaic(I,LAYER,J).eq.0.0)) THEN
                     SNOWH_mosaic(I,1:mosaic_cat,J) = 0.0
                 ENDIF
              ENDDO
              SNOWH(I,J) = 0.0
              ELSE
              SNOWH(I,J) = SNOWH_mosaic_avg(I,J)
              ENDIF
              
              IF (TRIGGER == 0) THEN
                  TVXY(I,J) = tvxy_mosaic_avg(I,J)/FAREA2
                  EAHXY(I,J) = eahxy_mosaic_avg(I,J)/FAREA2
                  TAHXY(I,J) = tahxy_mosaic_avg(I,J)/FAREA2
                  LFMASSXY(I,J) = lfmassxy_mosaic_avg(I,J)/FAREA2
                  RTMASSXY(I,J) = rtmassxy_mosaic_avg(I,J)/FAREA2
                  STMASSXY(I,J) = stmassxy_mosaic_avg(I,J)/FAREA2
                  WOODXY(I,J) = woodxy_mosaic_avg(I,J)/FAREA2
                  STBLCPXY(I,J) = stblcpxy_mosaic_avg(I,J)/FAREA2
                  FASTCPXY(I,J) = fastcpxy_mosaic_avg(I,J)/FAREA2
                  XLAIXY(I,J) = xlai_mosaic_avg(I,J)/FAREA2
                  XSAIXY(I,J) = xsaixy_mosaic_avg(I,J)/FAREA2
                  T2MVXY(I,J) = t2mvxy_mosaic_avg(I,J)/FAREA2
                  Q2MVXY(I,J) = q2mvxy_mosaic_avg(I,J)/FAREA2
                  NEEXY(I,J) = neexy_mosaic_avg(I,J)/FAREA2
                  GPPXY(I,J) = gppxy_mosaic_avg(I,J)/FAREA2
                  NPPXY(I,J) = nppxy_mosaic_avg(I,J)/FAREA2
                  FVEGXY(I,J) = fvegxy_mosaic_avg(I,J)/FAREA2
                  ECANXY(I,J) = ecanxy_mosaic_avg(I,J)/FAREA2
                  EDIRXY(I,J) = edirxy_mosaic_avg(I,J)/FAREA2
                  ETRANXY(I,J) = etranxy_mosaic_avg(I,J)/FAREA2
                  APARXY(I,J) = aparxy_mosaic_avg(I,J)/FAREA2
                  PSNXY(I,J) = psnxy_mosaic_avg(I,J)/FAREA2
                  SAVXY(I,J) = savxy_mosaic_avg(I,J)/FAREA2
                  RSSUNXY(I,J) = rssunxy_mosaic_avg(I,J)/FAREA2
                  RSSHAXY(I,J) = rsshaxy_mosaic_avg(I,J)/FAREA2
                  BGAPXY(I,J) = bgapxy_mosaic_avg(I,J)/FAREA2
                  WGAPXY(I,J) = wgapxy_mosaic_avg(I,J)/FAREA2
                  TGVXY(I,J) = tgvxy_mosaic_avg(I,J)/FAREA2
                  CHVXY(I,J) = chvxy_mosaic_avg(I,J)/FAREA2
                  SHBXY(I,J) = shbxy_mosaic_avg(I,J)/FAREA2
                  SHCXY(I,J) = shcxy_mosaic_avg(I,J)/FAREA2
                  EVGXY(I,J) = evgxy_mosaic_avg(I,J)/FAREA2
                  EVBXY(I,J) = evbxy_mosaic_avg(I,J)/FAREA2
                  IRGXY(I,J) = irgxy_mosaic_avg(I,J)/FAREA2
                  IRCXY(I,J) = ircxy_mosaic_avg(I,J)/FAREA2
                  TRXY(I,J) = trxy_mosaic_avg(I,J)/FAREA2
                  EVCXY(I,J) = evcxy_mosaic_avg(I,J)/FAREA2
                  CHLEAFXY(I,J) = chleafxy_mosaic_avg(I,J)/FAREA2
                  CHUCXY(I,J) = chucxy_mosaic_avg(I,J)/FAREA2
                  CHV2XY(I,J) = chv2xy_mosaic_avg(I,J)/FAREA2
              ELSE
                  TVXY(I,J) = -9999.
                  EAHXY(I,J) = -9999.
                  TAHXY(I,J) = -9999.
                  LFMASSXY(I,J) = -9999.
                  RTMASSXY(I,J) = -9999.
                  STMASSXY(I,J) = -9999.
                  WOODXY(I,J) = -9999.
                  STBLCPXY(I,J) = -9999.
                  FASTCPXY(I,J) = -9999.
                  XLAIXY(I,J) = -9999.
                  XSAIXY(I,J) = -9999.
                  T2MVXY(I,J) = -9999.
                  Q2MVXY(I,J) = -9999.
                  NEEXY(I,J) = -9999.
                  GPPXY(I,J) = -9999.
                  NPPXY(I,J) = -9999.
                  FVEGXY(I,J) = -9999.
                  ECANXY(I,J) = -9999.
                  EDIRXY(I,J) = -9999.
                  ETRANXY(I,J) = -9999.
                  APARXY(I,J) = -9999.
                  PSNXY(I,J) = -9999.
                  SAVXY(I,J) = -9999.
                  RSSUNXY(I,J) = -9999.
                  RSSHAXY(I,J) = -9999.
                  BGAPXY(I,J) = -9999.
                  WGAPXY(I,J) = -9999.
                  TGVXY(I,J) = -9999.
                  CHVXY(I,J) = -9999.
                  SHBXY(I,J) = -9999.
                  SHCXY(I,J) = -9999.
                  EVGXY(I,J) = -9999.
                  EVBXY(I,J) = -9999.
                  IRGXY(I,J) = -9999.
                  IRCXY(I,J) = -9999.2
                  TRXY(I,J) = -9999.
                  EVCXY(I,J) = -9999.
                  CHLEAFXY(I,J) = -9999.
                  CHUCXY(I,J) = -9999.
                  CHV2XY(I,J) = -9999.
                  
                  ! Set trigger back to 0
                  TRIGGER = 0

              ENDIF
              TGXY(I,J) = tgxy_mosaic_avg(I,J)
              CANICEXY(I,J) = canicexy_mosaic_avg(I,J)
              CANLIQXY(I,J) = canliqxy_mosaic_avg(I,J)
              
              
              CMXY(I,J) = cmxy_mosaic_avg(I,J)
              CHXY(I,J) = chxy_mosaic_avg(I,J)
              FWETXY(I,J) = fwetxy_mosaic_avg(I,J)
              SNEQVOXY(I,J) = sneqvoxy_mosaic_avg(I,J)
              ALBOLDXY(I,J) = alboldxy_mosaic_avg(I,J)
              QSNOWXY(I,J) = qsnowxy_mosaic_avg(I,J)
              QRAINXY(I,J) = qrainxy_mosaic_avg(I,J)
              WSLAKEXY(I,J) = wslakexy_mosaic_avg(I,J)
              ZWTXY(I,J) = zwtxy_mosaic_avg(I,J)
              WAXY(I,J) = waxy_mosaic_avg(I,J)
              WTXY(I,J) = wtxy_mosaic_avg(I,J)
              
              TAUSSXY(I,J) = taussxy_mosaic_avg(I,J)
              SMCWTDXY(I,J) = smcwtdxy_mosaic_avg(I,J)
              DEEPRECHXY(I,J) = deeprechxy_mosaic_avg(I,J)
              RECHXY(I,J) = rechxy_mosaic_avg(I,J)

              !out variables only
              
              T2MBXY(I,J) = t2mbxy_mosaic_avg(I,J)
              
              Q2MBXY(I,J) = q2mbxy_mosaic_avg(I,J)
              TRADXY(I,J) = tradxy_mosaic_avg(I,J)
              
              RUNSFXY(I,J) = runsfxy_mosaic_avg(I,J)
              RUNSBXY(I,J) = runsbxy_mosaic_avg(I,J)
              
              FSAXY(I,J) = fsaxy_mosaic_avg(I,J)
              FIRAXY(I,J) = firaxy_mosaic_avg(I,J)
              SAGXY(I,J) = sagxy_mosaic_avg(I,J)
              TGBXY(I,J) = tgbxy_mosaic_avg(I,J)
              CHBXY(I,J) = chbxy_mosaic_avg(I,J)
              SHGXY(I,J) = shbxy_mosaic_avg(I,J)
              
              GHVXY(I,J) = ghvxy_mosaic_avg(I,J)
              GHBXY(I,J) = ghbxy_mosaic_avg(I,J)
             
              IRBXY(I,J) = irbxy_mosaic_avg(I,J)
              CHB2XY(I,J) = chb2xy_mosaic_avg(I,J)


              ! irrigation intermediate variables
              IRWATSI(I,J) =  IRWATSI_mosaic_avg(I,J)
              IRWATMI(I,J) =  IRWATMI_mosaic_avg(I,J)
              IRWATFI(I,J) =  IRWATFI_mosaic_avg(I,J)
              IRELOSS(I,J) =  IRELOSS_mosaic_avg(I,J)
              IRSIVOL(I,J) =  IRSIVOL_mosaic_avg(I,J)
              IRMIVOL(I,J) =  IRMIVOL_mosaic_avg(I,J)
              IRFIVOL(I,J) =  IRFIVOL_mosaic_avg(I,J)
              IRRSPLH(I,J) =  IRRSPLH_mosaic_avg(I,J)

              IRNUMSI(I,J) = IRNUMSI_mosaic_avg(I,J)/mosaic_cat
              IRNUMMI(I,J) = IRNUMMI_mosaic_avg(I,J)/mosaic_cat
              IRNUMFI(I,J) = IRNUMFI_mosaic_avg(I,J)/mosaic_cat

              ! Extra variables averaged needed for extra outputs
              QINTSXY(I,J) = QINTSXY_mosaic_avg(I,J)
              QINTRXY(I,J) = QINTRXY_mosaic_avg(I,J)
              QDRIPSXY(I,J) = QDRIPSXY_mosaic_avg(I,J)
              QDRIPRXY(I,J) = QDRIPRXY_mosaic_avg(I,J)
              QTHROSXY(I,J) = QTHROSXY_mosaic_avg(I,J)
              QTHRORXY(I,J) = QTHRORXY_mosaic_avg(I,J)
              QSNSUBXY(I,J) = QSNSUBXY_mosaic_avg(I,J)
              QSNFROXY(I,J) = QSNFROXY_mosaic_avg(I,J)
              QSUBCXY(I,J) = QSUBCXY_mosaic_avg(I,J)
              QFROCXY(I,J) = QFROCXY_mosaic_avg(I,J)
              QEVACXY(I,J) = QEVACXY_mosaic_avg(I,J)
              QDEWCXY(I,J) = QDEWCXY_mosaic_avg(I,J)
              QFRZCXY(I,J) = QFRZCXY_mosaic_avg(I,J)
              QMELTCXY(I,J) = QMELTCXY_mosaic_avg(I,J)
              QSNBOTXY(I,J) = QSNBOTXY_mosaic_avg(I,J)
              QMELTXY(I,J) = QMELTXY_mosaic_avg(I,J)
              PONDINGXY(I,J) = PONDINGXY_mosaic_avg(I,J)
              PAHXY(I,J) = PAHXY_mosaic_avg(I,J)
              PAHVXY(I,J) = PAHVXY_mosaic_avg(I,J)
              PAHBXY(I,J) = PAHBXY_mosaic_avg(I,J)
              PAHGXY(I,J) = PAHGXY_mosaic_avg(I,J)
              FPICEXY_mosaic_avg(I,J) = FPICEXY_mosaic_avg(I,J)

              ACC_SSOILXY(I,J) = ACC_SSOILXY_mosaic_avg(I,J)
              ACC_QINSURXY(I,J) = ACC_QINSURXY_mosaic_avg(I,J)
              ACC_QSEVAXY(I,J) = ACC_QSEVAXY_mosaic_avg(I,J)
              EFLXBXY(I,J) = EFLXBXY_mosaic_avg(I,J)
              SOILENERGY(I,J) = SOILENERGY_mosaic_avg(I,J)
              SNOWENERGY(I,J) = SNOWENERGY_mosaic_avg(I,J)
              CANHSXY(I,J) = CANHSXY_mosaic_avg(I,J)
              ACC_DWATERXY(I,J) = ACC_DWATERXY_mosaic_avg(I,J)
              ACC_PRCPXY(I,J) = ACC_PRCPXY_mosaic_avg(I,J)
              ACC_ECANXY(I,J) = ACC_ECANXY_mosaic_avg(I,J)
              ACC_EDIRXY(I,J) = ACC_EDIRXY_mosaic_avg(I,J)
              ACC_ETRANXY(I,J) = ACC_ETRANXY_mosaic_avg(I,J)

               DO LAYER=1,3
                      TSNOXY(I, LAYER-3,J) = tsnoxy_mosaic_avg(I, LAYER,J)
                      snicexy(I,LAYER-3,J) = snicexy_mosaic_avg(I,LAYER,J)
                      snliqxy(I,LAYER-3,J) = SNLIQXY_mosaic_avg(I,LAYER,J)
               ENDDO

               DO LAYER=1,7
                      zsnsoxy(I,LAYER-3,J) = zsnsoxy_mosaic_avg(I,LAYER,J)
               ENDDO
              !HUE THINGS
              RUNONSFXY(I,J) = RUNONSFXY_mosaic_avg(I,J)


            ENDIF ! End of the Soil-Water_ice If else statement

           ENDDO ILOOP !END of the I LOOP
           ENDDO JLOOP !END of the J LOOP

         END SUBROUTINE NOAHMPLSM_MOSAIC_HUE


END MODULE module_sf_noahmpdrv
