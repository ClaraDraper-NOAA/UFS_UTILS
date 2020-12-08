MODULE M_Snow_Analysis

USE NETCDF
USE M_DA
!USE MPI
Use, Intrinsic :: IEEE_ARITHMETIC       

CONTAINS

 subroutine Snow_Analysis_OI(SNOW_OI_TYPE, MAX_TASKS, MYRANK, NPROCS, IDIM, JDIM, IY, IM, ID, IH, &  
                                num_assim_steps, dT_Asssim,  & 
                                LENSFC, IVEGSRC, PERCENT_OBS_WITHHELD, &
                                L_horz , h_ver, obs_tolerance, &
                                obs_srch_rad, bkgst_srch_rad, max_num_nearStn, max_num_nearIMS, &                                
                                ims_max_ele, num_subgrd_ims_cels, &
                                assim_SnowPack_obs, assim_SnowCov_obs, ims_correlated, stn_var, &
                                GHCND_SNOWDEPTH_PATH, IMS_SNOWCOVER_PATH, &
                                IMS_INDEXES_PATH, SFC_FORECAST_PREFIX, &
                                SNOANL)
                                                        
        !----------------------------------------------------------------------
        ! Input arguments: 
        ! IDIM * JDIM = LENSFC: number of grid cells in tile = xdim * ydim   
        ! IY, IM, ID, IH = year, month, day, hour of current model step   
        ! MYRANK: rank/id of the MPI process
        ! ...
        !
        ! Inputs, read from file:
        ! RLA, RLO: lat lon information for the tile
        ! SNOFCS(LENSFC): forecast  snowdepth, or snow water equivalent forecat, as requested
        !
        ! Outputs:
        ! SNOANL: snow analysis (of SWE or SND, see below)
        ! 
        ! Draper - changes to snow variables to clarify names, removed unnecesary conversions
        !          SWE - snow water equivalent 
        !          SND - snow depth 
        !          SNO - generic variable to refer to either (determined by value of 
        !                snow_OI_type
        !                   
        !----------------------------------------------------------------------
        IMPLICIT NONE
        !
        include 'mpif.h'
        
        integer, parameter :: dp = kind(1.d0)

        INTEGER, intent(in)    :: SNOW_OI_TYPE, MAX_TASKS, MYRANK, NPROCS, IDIM, JDIM, &
                                  IY, IM, ID, IH, LENSFC, IVEGSRC, num_assim_steps 
        LOGICAL, intent(in)    :: assim_SnowPack_obs, assim_SnowCov_obs, ims_correlated
        CHARACTER(len=4), intent(in)   :: stn_var ! should probably be called control_var
        CHARACTER(LEN=*), Intent(In)   :: GHCND_SNOWDEPTH_PATH, IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH
        CHARACTER(LEN=*), Intent(In)   :: SFC_FORECAST_PREFIX
!  GHCND_SNOWDEPTH_PATH="/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/GHCND/",
!  IMS_SNOWCOVER_PATH="/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS/",
!  IMS_INDEXES_PATH="/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS_INDEXES/",
!  SFC_FORECAST_PREFIX = "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/fnbgsio/snow."
        REAL, intent(In)    :: PERCENT_OBS_WITHHELD, dT_Asssim 
        Real, intent(In)    :: L_horz , h_ver, obs_tolerance
        Real, intent(In)    :: obs_srch_rad, bkgst_srch_rad, ims_max_ele        
        INTEGER, intent(in) :: max_num_nearStn, max_num_nearIMS, num_subgrd_ims_cels 
        REAL, intent(Out)   :: SNOANL(LENSFC) 
        CHARACTER(LEN=5)    :: TILE_NUM
        Character(LEN=3)    :: rank_str
        INTEGER                     :: IERR     
        REAL        :: RLA(LENSFC), RLO(LENSFC), RLO_Tile(LENSFC), OROG(LENSFC)  !, OROG_UF(LENSFC)
        REAL                :: SNDFCS(LENSFC), SWEFCS(LENSFC), SNDANL(LENSFC), SWEANL(LENSFC)
        REAL                :: SNOFCS(LENSFC), VETFCS(LENSFC), SNUP_Array(LENSFC)
        INTEGER             :: LANDMASK(LENSFC)
        CHARACTER(len=250)   :: dim_name, ghcnd_inp_file, ims_inp_file, ims_inp_file_indices
        CHARACTER(len=5)     :: y_str, m_str, d_Str, h_str, fvs_tile
        REAL, ALLOCATABLE    :: SNOOBS_stn(:), SNOFCS_at_stn(:)               
        REAL, ALLOCATABLE    :: Lat_stn(:), Lon_stn(:), OROGFCS_at_stn(:)  
        REAL                 :: lat_min, lat_max, lon_min, lon_max      
        Real                 :: SNCOV_IMS(LENSFC)  ! ims resampled at each grid
        Real                 :: SNO_IMS_at_Grid(LENSFC)

        INTEGER :: num_stn, Num_Ims, num_Eval !num_subgrd_ims_cels !Num_Ims_Lat, Num_Ims_Lon
        !Real    :: bkgst_srch_rad   ! radius_of_influence for selecting state at observation point
        INTEGER :: jndx, zndx, ncol, nrow
        Integer, Allocatable   :: index_back_at_nearStn(:), index_back_at_nearIMS(:) !loc_near_Obs(:), 
        Integer                :: num_loc, num_loc_1, num_loc_2
        Real, Parameter         :: Stdev_back_depth = 30., Stdev_Obsv_depth = 40., Stdev_Obsv_ims = 80. ! mm 
        real                    :: stdev_obsv, stdev_back
        Integer                         :: ims_assm_hour
        !Real                           :: obs_tolerance, ims_max_ele

        Real(dp), Allocatable    :: B_cov_mat(:,:), b_cov_vect(:)
        Real(dp), Allocatable    :: O_cov_mat(:,:), W_wght_vect(:)
        Real, Allocatable        :: back_at_Obs(:), obs_Array(:), Lat_Obs(:), Lon_Obs(:), orogfcs_Obs(:)
        REAL                     :: innov_at_Grid(LENSFC)
        Real, Allocatable        :: obs_Innov(:), OmB_innov_at_stn(:)

        CHARACTER(len=250)       :: forc_inp_path, da_out_file  !anl_inp_path,
        CHARACTER(LEN=3)         :: RANKCH 

        REAL, ALLOCATABLE  :: SNOFCS_atEvalPts(:), innov_atEvalPts(:), SNOANL_atEvalPts(:)  !evalution points 
        REAL, ALLOCATABLE  :: Lat_atEvalPts(:), Lon_atEvalPts(:), Obs_atEvalPts(:)     !evalution points
        Integer, ALLOCATABLE    :: index_back_atEval(:)     ! background locations at eval points 
        Integer, ALLOCATABLE    :: index_back_atObs(:)   ! the location of background corresponding obs

        Real               :: snodens, SNODENS_DIST(LENSFC)
        LOGICAL            :: assim_snpack_stn, assim_SWE    !use swe as model state var? (instead of snow depth)
        LOGICAL            :: assim_sncov, assim_sncov_thisGridCell    !assimilate sncov, 

        Integer            :: veg_type_landice  ! 10.21.20: no assmn over land ice
    ! for mpi par
        INTEGER            :: Np_ext, Np_til, p_tN, p_tRank, N_sA, N_sA_Ext, mp_start, mp_end
        INTEGER            :: send_proc, rec_stat(MPI_STATUS_SIZE), dest_Aoffset, pindex
        INTEGER            :: mpiReal_size, rsize
        REAL               :: tmp
        INTEGER            :: istep, IY_loc, IM_loc, ID_loc, IH_loc
        REAL               :: IH_real
        INTEGER            :: num_print

!=============================================================================================
! 1. initialise vars,set-up processors, and read lat/lon from orog files.
!=============================================================================================
        assim_snpack_stn = assim_SnowPack_obs ! if assimilating GHCND, doesn't have SWE obs. need to convert to SWE
        assim_sncov =  assim_SnowCov_obs  

        If (SNOW_OI_TYPE == 2) then 
                assim_SWE = .True.      ! note: if this is set true, need to establish stdev_obsv and stdev_back
                                        !       and need to code readers for appropriate data sets
                stdev_obsv = 0.1 * stdev_obsv_depth   ! assumes snodens = 0.1; need to refined
                stdev_back = 0.1 * stdev_back_depth
        ELSEIF (SNOW_OI_TYPE == 1) then 
                assim_SWE = .False.                 
                stdev_obsv = stdev_obsv_depth
                stdev_back = stdev_back_depth    
        ENDIF
        !obs_srch_rad = 250. ! radius of observation search
        ims_assm_hour = 18 
        ! noah models specific? Needed to ID glaciers.
        if (IVEGSRC == 2) then   ! sib
                veg_type_landice=13
        else
                veg_type_landice=15
        endif

!total number of processors used = Multiple of 6: any extra processors sent to end of subroutine
        IF (myrank ==0) PRINT*,"Total num proc ", NPROCS, " Num tiles /Max.tasks: ", MAX_TASKS

        Np_ext = MOD(NPROCS, MAX_TASKS)  ! extra/inactive procs
        if (MYRANK >  NPROCS - Np_ext - 1) goto 999
        Np_til = NPROCS / MAX_TASKS  ! num proc. per tile 
        p_tN = MOD(MYRANK, MAX_TASKS)  ! tile for proc.
        p_tRank = MYRANK / MAX_TASKS  ! proc. rank within tile
        N_sA = LENSFC / Np_til  ! sub array length per proc
        N_sA_Ext = LENSFC - N_sA * Np_til ! extra grid cells
        if(p_tRank == 0) then 
                mp_start = 1
        else
                mp_start = p_tRank * N_sA + N_sA_Ext + 1   ! start index of subarray for proc
        endif
        mp_end = (p_tRank + 1) * N_sA + N_sA_Ext                ! end index of subarray for proc        
        If(myrank == 0 )PRINT*,"sub array length ", N_sA, " extra sub array: ", N_sA_Ext
! if (p_tN /= 4 ) goto 999

! must have at least one kind of obs: depth or area. If no obs files: go to end 
        if (GHCND_SNOWDEPTH_PATH(1:8).eq.'        ') assim_snpack_stn = .false.
        if ((IMS_SNOWCOVER_PATH(1:8).eq.'        ') .OR. (IMS_INDEXES_PATH(1:8).eq.'        ')) then
                assim_sncov = .false.
        end if
        if( (.not. assim_snpack_stn) .and. (.not. assim_sncov) ) then
                print*, "Observation paths don't exist!, skipping the OI DA"
                goto 999
        end if

! READ THE OROGRAPHY AND GRID POINT LAT/LONS FOR THE CUBED-SPHERE TILE p_tN
        CALL READ_LAT_LON_OROG_atRank(p_tN, RLA,RLO,OROG,TILE_NUM,IDIM,JDIM,LENSFC) !OROG_UF,
        if (print_deb) PRINT*,"Snow anl on ", MYRANK, " Tile group: ", p_tN, " Tile: ", TILE_NUM
        
        RLO_Tile = RLO ! copy so that RLO used later is not modified
        Where(RLO_Tile > 180) RLO_Tile = RLO_Tile - 360
        lat_min = MAX(MINVAL(RLA) - 1., -90.) ! CSD - why is the limit 1?
        lat_max = MIN(MAXVAL(RLA) + 1., 90.)
        lon_min = MAX(MINVAL(RLO_Tile) - 1., -180.)
        lon_max = MIN(MAXVAL(RLO_Tile) + 1., 180.)      
        if (p_tN==3)  then  
                lon_min = 125.       ! lon_min is left, lon_max is right, not necessary min/max value
                lon_max = -145.
        endif   
        if ((p_tRank==0) .and. print_deb) then
                print*, TILE_NUM, " min/max lat/lon ", lat_min, lat_max, lon_min, lon_max
        endif

! If multiple time steps are simulated in a given time period (window) given by num_assim_steps * dT_Assim
! Note: the DA outputs from the last time step are returned        
        IH_real = IH; IY_loc = IY; IM_loc=IM; ID_loc=ID; IH_loc=IH; ! these variables change inside loop below
    Do istep = 1, num_assim_steps        
                
        write(y_str, "(I4)") IY_loc
        write(m_str, "(I0.2)") IM_loc
        write(d_str, "(I0.2)") ID_loc
        write(h_str, "(I0.2)") IH_loc
        write(fvs_tile, "(I3)") IDIM

! controls calling of obs operator for stn data. If remains 0 will not be called.
        num_stn = 0 
        num_Eval = 0

!=============================================================================================
! 2. Read model forecast here, as need VETFCS and snow density for IMS snow depth calc. (later, separate read routines) 
!=============================================================================================

       ! READ THE INPUT SURFACE DATA ON THE CUBED-SPHERE TILE p_tN. A
       ! Also get vegtype (VETFCS) to identify glacier.    

        if (SFC_FORECAST_PREFIX(1:8).eq.'        ') then
                ! FNBGSI = "./fnbgsi." // RANKCH
                WRITE(RANKCH, '(I3.3)') (p_tN+1)
                forc_inp_path = "./fnbgsi." // RANKCH
        else
                forc_inp_path = TRIM(SFC_FORECAST_PREFIX)//TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"."//TRIM(h_str)// &
                                     "0000.sfc_data."//TILE_NUM//".nc"
        end if
                                     
        Call READ_Forecast_Data_atPath(forc_inp_path, veg_type_landice, LENSFC, SWEFCS, SNDFCS, &
                                       VETFCS, LANDMASK)

        if (assim_SWE) then 
           SNOFCS=SWEFCS  
        else 
           SNOFCS=SNDFCS 
        endif
        ! average snow density of snow over non-glacier land.
        SNODENS_DIST = SWEFCS/SNDFCS
        snodens = SUM(SWEFCS/SNDFCS, Mask = (LANDMASK==1 .and. SNDFCS>0.01 )) &
                         / COUNT (LANDMASK==1 .and. SNDFCS> 0.01) 
        ! for grid cells with no valid density, fill in the average snodens
        Where(.not.(LANDMASK==1 .and. SNDFCS>0.01 )) SNODENS_DIST = snodens
        If (p_tRank==0)  print*, "Tile ", p_tN, ' mean snow density', snodens
        tmp = SUM(SWEFCS,  Mask = (LANDMASK==1 .and. SNDFCS>0.01 )) &
                         / COUNT (LANDMASK==1 .and. SNDFCS> 0.01)
        If (p_tRank==0)  print*, "Tile ", p_tN,  ' mean SWE', tmp
        tmp = SUM(SNDFCS,  Mask = (LANDMASK==1 .and. SNDFCS>0.01 )) &
                         / COUNT (LANDMASK==1 .and. SNDFCS> 0.01)
        If (p_tRank==0)  print*, "Tile ", p_tN,  ' mean SND', tmp
        
!=============================================================================================
! 3. Read observations
!=============================================================================================

! 3a. Read station obs (of snow depth or SWE)

        if (assim_snpack_stn) then 
             ghcnd_inp_file = TRIM(GHCND_SNOWDEPTH_PATH)//"GHCND.SNWD."// &
                                         TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//TRIM(h_str)//".nc"
             dim_name = "Site_Id"    
             call Observation_Read_GHCND_Tile_excNaN(p_tN, ghcnd_inp_file, dim_name, &
                        lat_min, lat_max, lon_min, lon_max, & 
                        num_stn, SNOOBS_stn,              &
                        Lat_stn, Lon_stn, MYRANK) 

            if ((p_tRank==0) .and. print_deb) then
                    print*, "Tile ", p_tN, " num. Stn obs ", num_stn
            endif
            if ((p_tRank==0) .and. (p_tN==2) .and. print_deb) then
                    PRINT*, "Stn SND from rank: ", MYRANK
                    PRINT*, SNOOBS_stn
                    PRINT*, "Lat at Stn from rank: ", MYRANK
                    PRINT*, Lat_stn
                    PRINT*, "Lon at Stn from rank: ", MYRANK
                    PRINT*, Lon_stn
            endif
            if (myrank==0) PRINT*,'Finished reading station data'

       endif

! CSD beyond this point, there should be no specific mention of the station data source

! 3b. Read remotely sensed snow cover, and convert to  snow depth or SWE. 

       if (assim_sncov) then
            ! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS/"
            ims_inp_file = TRIM(IMS_SNOWCOVER_PATH)//"IMS.SNCOV."// &
                                       TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//TRIM(h_str)//".nc"                      !
            ! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS_INDEXES/"
            ims_inp_file_indices = TRIM(IMS_INDEXES_PATH)//"C"//TRIM(ADJUSTL(fvs_tile))//&
                                                       ".IMS.Indices."//TRIM(TILE_NUM)//".nc"                       
            Call Observation_Read_IMS_Full(ims_inp_file, ims_inp_file_indices, &
                                                       MYRANK, JDIM, IDIM, num_subgrd_ims_cels, SNCOV_IMS)
            if((p_tRank==0) .and. (p_tN==1) .and. print_deb) then
                    PRINT*, "SNCOV from rank: ", MYRANK
                    PRINT*, SNCOV_IMS
            endif 
            if (myrank==0) PRINT*,'Finished reading SNCOV, converting to snow depth' 

            call CalcSWEFromSnowCover(SNCOV_IMS, VETFCS, LENSFC, SNO_IMS_at_Grid, SNUP_Array)

            if (.not. assim_SWE)  SNO_IMS_at_Grid = SNO_IMS_at_Grid/SNODENS_DIST ! convert SWE to SND
            
            if ((p_tN==4) .and. (p_tRank==0) .and. print_deb) then
                    PRINT*, "SNCOV obs at each grid cell from rank: ", MYRANK
                    PRINT*, SNO_IMS_at_Grid
            endif
            if (myrank==0) PRINT*,'Finished converting SNCOV observations'
        
        endif ! read_IMS 

!=============================================================================================
! 4. Get H(x): Read forecast snow fields from restart files, then interpolate to obs location.
!=============================================================================================

! 4a. read the forecast file on model grid : this was done earlier, as need veg type and 
!      snow density for IMS snow depth conversion

! 4b. get H(x) for station obs
        ! Get model states at obs points
        if (num_stn > 0) then ! skip if not reading in station data / no obs were available
            ALLOCATE(SNOFCS_at_stn(num_stn))
            ALLOCATE(OROGFCS_at_stn(num_stn)) 
            ALLOCATE(index_back_atObs(num_stn)) 
            ALLOCATE(OmB_innov_at_stn(num_stn)) 
             ! using PERCENT_OBS_WITHHELD % of stn locations for evaluation
            num_Eval = floor(0.01 * PERCENT_OBS_WITHHELD * num_stn)  
            if (num_Eval > 0) then 
                ALLOCATE(index_back_atEval(num_Eval)) 
                ALLOCATE(Obs_atEvalPts(num_Eval)) 
                ALLOCATE(SNOFCS_atEvalPts(num_Eval)) 
                ALLOCATE(Lat_atEvalPts(num_Eval))
                ALLOCATE(Lon_atEvalPts(num_Eval)) 
                ALLOCATE(innov_atEvalPts(num_Eval))
                ALLOCATE(SNOANL_atEvalPts(num_Eval))    
                if(p_tRank == 0) then 
                        PRINT*, "Tile ", p_tN+1, " ", num_Eval, ' points for evaluation excluded from DA'       
                endif  
            endif 

    ! CSD todo: separate out eval call and obs operator call. model info (SNOOBS_stn shouldn't be in the obs operator)
    !           for JEDI, probably also want to add snow depth derived from snow cover here.
            Call Observation_Operator_Parallel(Myrank, MAX_TASKS, p_tN, p_tRank, Np_til, & 
                                RLA, RLO, OROG, Lat_stn, Lon_stn,   &
                                LENSFC, num_stn, num_Eval, bkgst_srch_rad, SNOFCS, SNOOBS_stn,  &
                                SNOFCS_at_stn, OROGFCS_at_stn, index_back_atObs, index_back_atEval, &
                                Obs_atEvalPts, SNOFCS_atEvalPts, Lat_atEvalPts, Lon_atEvalPts)
            if ((p_tN==4) .and. (p_tRank==0) .and. print_deb) then 
                    PRINT*, "Background Indices at eval points"
                    PRINT*, index_back_atEval       
                    PRINT*, "Obs at Eval Points" 
                    PRINT*, Obs_atEvalPts   
                    PRINT*, "Forecast at Eval Points"
                    PRINT*, SNOFCS_atEvalPts                             
                    PRINT*, "Lat at Eval Points"
                    PRINT*, Lat_atEvalPts
                    PRINT*, "Lon at Eval Points"
                    PRINT*, Lon_atEvalPts
            endif

            OmB_innov_at_stn = SNOOBS_stn - SNOFCS_at_stn

            if ((p_tN==4) .and. (p_tRank==0) .and. print_deb) then
                    PRINT*, "station Lat range from rank: ", MYRANK, MINVAL(Lat_stn), " ", MAXVAL(Lat_stn)
                    PRINT*, "station Lon range from rank: ", MYRANK, MINVAL(Lon_stn), " ", MAXVAL(Lon_stn)
                !     PRINT*, "Lat at Obs Points"
                !     PRINT*, Lat_stn
                !     PRINT*, "Lon at Obs Points"
                !     PRINT*, Lon_stn
                    PRINT*, "Model elevation at station locations from rank: ", MYRANK
                    PRINT*, OROGFCS_at_stn
                    PRINT*, "Background Indices at obs points"
                    PRINT*, index_back_atObs
                    PRINT*, "Background ", stn_var, " at station locations from rank: ", MYRANK
                    PRINT*, SNOFCS_at_stn  
                !     PRINT*, "Obs at obs stns" 
                !     PRINT*, SNOOBS_stn   
                    PRINT*, "O - B (innovation at obs points)"
                    PRINT*, OmB_innov_at_stn 
            endif
            if (myrank==0) PRINT*,'Finished observation operator for station data'         
        endif ! num_stn > 0

!=============================================================================================
! 5.  Obs-based obs QC goes here
!=============================================================================================

!CSDCSD - todo. Add QC here.

! QC steps:
! if station elevation >1500, discard stn_obs
! is IMS observations > threshold, discard IMS obs * 
! min/max limits on station obs
! temperature check on all obs 

! if abs (model - obs ) elevation > ??,  discard obs 

! QC obs for land cover (below QC's update cell, but not obs) 
! gross error check * 
! screen stn obs for IMS snow cover 
! screen snow cover-derived snow depth (IMS) if model has snow  * 

!=============================================================================================
! 6. Perform the DA update, by looping over each grid cell 
!=============================================================================================
 
        !obs_srch_rad = 250. ! radius of observation search
        if (myrank==0) PRINT*,'Starting DA loop'
        Do jndx = mp_start, mp_end     !LENSFC/2, LENSFC ! 442369, 442370   !
                num_loc_1 = 0
                num_loc_2 = 0
                assim_sncov_thisGridCell = .FALSE.
                if (print_deb) print*, "proc ", myrank, " grid: ", jndx
                if(num_stn>0) then 
                ! currently: find station obs in radius, do gross error check, and limit to 50 obs
                ! QC: gross error check is done in this call.
                        call nearest_Observations_Locations(RLA(jndx), RLO(jndx),    &
                                        Lat_stn, Lon_stn,  num_stn, obs_srch_rad, max_num_nearStn,   &
                                        stdev_back, stdev_obsv, obs_tolerance,                 &
                                        SNOFCS_at_stn, SNOOBS_stn,                                              &
                                        index_back_at_nearStn,  num_loc_1) !,      &LENSFC,
                        if (print_deb) print*, "number of stn sndpth obs ", num_loc_1
                endif         
                ! legacy code to collect IMS obs from surrounding grid cells.                 
                ! if ( assim_sncov .AND. (IH_loc == ims_assm_hour) ) then
                !         Call nearest_IMS_Locations(LENSFC, max_num_nearIMS, obs_srch_rad, &
                !                         ims_max_ele, RLA(jndx), RLO(jndx), RLA, RLO, OROG, & 
                !                         SNOFCS, SWEFCS, SNUP_Array, &
                !                         SNO_IMS_at_Grid, num_loc_2, index_back_at_nearIMS)
                !         if (num_loc_2 > 0) assim_sncov_thisGridCell = .TRUE.
                !         if (print_deb) print*, " number of IMS obs: ", num_loc_2                                                                       
                ! endif
                ! QC: -for this grid cell, test IMS above ims_max_ele 
                !     -if model and IMS have 100% snow cover, delete the IMS snow depth value
                !     
                if( assim_sncov .AND. (IH_loc == ims_assm_hour) .AND. &
                    (.NOT. IEEE_IS_NAN(SNOFCS(jndx))) .AND. &
                    (.NOT. IEEE_IS_NAN(SNO_IMS_at_Grid(jndx))) .AND. &
                    (OROG(jndx) <= ims_max_ele) .AND. &
                    ( .not.( (SWEFCS(jndx) >= SNUP_Array(jndx)) .AND. & 
                       SNO_IMS_at_Grid(jndx) >= SNUP_Array(jndx)) ) ) then
                        num_loc_2 = 1 
                        assim_sncov_thisGridCell = .TRUE.                
                endif
                num_loc = num_loc_1 + num_loc_2                
                ! if assim_sncov=false >> num_loc_1=num_loc
                ! QC: only update this grid cell if it is land.
                if((num_loc > 0) .and. ( LANDMASK(jndx) == 1 )) then 
                        ! get background states
                        Allocate(back_at_Obs(num_loc))
                        Allocate(obs_Array(num_loc))
                        Allocate(Lat_Obs(num_loc))
                        Allocate(Lon_Obs(num_loc))
                        Allocate(orogfcs_Obs(num_loc))                        
                        if(num_loc_1 > 0) then
                                Do zndx = 1, num_loc_1     
                                        back_at_Obs(zndx) = SNOFCS_at_stn(index_back_at_nearStn(zndx))
                                        obs_Array(zndx) = SNOOBS_stn(index_back_at_nearStn(zndx))
                                        Lat_Obs(zndx) = Lat_stn(index_back_at_nearStn(zndx))
                                        Lon_Obs(zndx) = Lon_stn(index_back_at_nearStn(zndx))
                                        orogfcs_Obs(zndx) = OROGFCS_at_stn(index_back_at_nearStn(zndx)) 
                                End Do
                        End if
                        !multi-IMS
                        ! if(assim_sncov_thisGridCell) then                                
                        !         Do zndx = 1, num_loc_2
                        !                 back_at_Obs(num_loc_1+zndx) = SNOFCS(index_back_at_nearIMS(zndx))
                        !                 obs_Array(num_loc_1+zndx) = SNO_IMS_at_Grid(index_back_at_nearIMS(zndx))
                        !                 Lat_Obs(num_loc_1+zndx) = RLA(index_back_at_nearIMS(zndx))   !Lat_IMS_atGrid(jndx)
                        !                 Lon_Obs(num_loc_1+zndx) = RLO(index_back_at_nearIMS(zndx))   !Lon_IMS_atGrid(jndx)
                        !                 orogfcs_Obs(num_loc_1+zndx) = OROG(index_back_at_nearIMS(zndx))  !Ele_IMS(jndx)
                        !         End Do
                        ! endif

                        ! Append IMS-derived snow depth to the obs array 
                        if(assim_sncov_thisGridCell) then                                
                                back_at_Obs(num_loc) = SNOFCS(jndx)
                                obs_Array(num_loc) = SNO_IMS_at_Grid(jndx)
                                Lat_Obs(num_loc) = RLA(jndx)   
                                Lon_Obs(num_loc) = RLO(jndx) 
                                orogfcs_Obs(num_loc) = OROG(jndx)
                        endif
                        ! compute covariances
                        Allocate(B_cov_mat(num_loc, num_loc))
                        Allocate(b_cov_vect(num_loc))
                        Allocate(O_cov_mat(num_loc, num_loc))
                        Allocate(W_wght_vect(num_loc))   
                        ! call compute_covariances_multiIMS(RLA(jndx), RLO(jndx), OROG(jndx), SNOFCS(jndx),    &
                        !         Lat_Obs, Lon_Obs, orogfcs_Obs, num_loc, num_loc_1, num_loc_2,                &
                        !         Stdev_back, Stdev_Obsv, Stdev_Obsv_ims,      &
                        !         L_horz, h_ver,                                   &   !L_horz in Km, h_ver in m
                        !         assim_sncov_thisGridCell, ims_correlated,            &
                        !         B_cov_mat, b_cov_vect, O_cov_mat, W_wght_vect)
                        call compute_covariances(RLA(jndx), RLO(jndx), OROG(jndx), SNOFCS(jndx),    &
                                Lat_Obs, Lon_Obs, orogfcs_Obs, num_loc,   &
                                Stdev_back, Stdev_Obsv, Stdev_Obsv_ims,      &
                                L_horz, h_ver,                                   &   !L_horz in Km, h_ver in m
                                assim_sncov_thisGridCell,                          &
                                B_cov_mat, b_cov_vect, O_cov_mat, W_wght_vect)
                        Allocate(obs_Innov(num_loc))
                        call Snow_DA_OI(back_at_Obs, obs_Array, num_loc, W_wght_vect,            &
                                SNOFCS(jndx), innov_at_Grid(jndx), SNOANL(jndx), obs_Innov)
                        if ((p_tN==4) .and. (p_tRank==0) .and. print_deb) then  !
                                print*, "proc ", myrank, "loop ", jndx, "num depth obs ", num_loc_1, "total obs", num_loc
                                PRINT*, " background at obs pts: "
                                PRINT*, back_at_Obs     
                                PRINT*, "Observed"
                                PRINT*,  obs_Array
                                PRINT*, "Obs innovation: "
                                PRINT*, obs_Innov
                                PRINT*, "Weight vector: "
                                PRINT*, W_wght_vect     
                                print*, "innov: ", innov_at_Grid(jndx), "forec: ", SNOFCS(jndx), " anl: ", SNOANL(jndx)
                        endif           
                        !free mem
                        DEALLOCATE(back_at_Obs, obs_Array)
                        DEALLOCATE(Lat_Obs, Lon_Obs, orogfcs_Obs, obs_Innov)
                        DEALLOCATE(B_cov_mat, b_cov_vect, O_cov_mat, W_wght_vect)

                else
                        SNOANL(jndx) = SNOFCS(jndx) ! if not land, or no obs avail keep forecast
                endif
                if (allocated(index_back_at_nearStn))  Deallocate(index_back_at_nearStn) 
                if (allocated(index_back_at_nearIMS))  Deallocate(index_back_at_nearIMS)                 
        End do
        CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
        if (myrank==0) PRINT*, ' Finished DA loops'
        
!=============================================================================================
! 7. Write debug outputs
!=============================================================================================

! ToDO: Better way to handle this? ! CSD - I'll fix this later.
! Real data type size corresponding to mpi
        rsize = SIZEOF(snodens)
        Call MPI_TYPE_SIZE(MPI_REAL, mpiReal_size, IERR) 
        If (rsize == 4 ) then 
                mpiReal_size = MPI_REAL4
        elseif (rsize == 8 ) then 
                mpiReal_size = MPI_REAL8
        elseif (rsize == 16 ) then 
                mpiReal_size = MPI_REAL16
        else
                PRINT*," Possible mismatch between Fortran Real ", rsize," and Mpi Real ", mpiReal_size
                Stop
        endif
        ! send analyses arrays to 'tile-level root' proc.               
        if (MYRANK > (MAX_TASKS - 1) ) then
                call MPI_SEND(SNOANL(mp_start:mp_end), N_sA, mpiReal_size, p_tN,   &
                                          MYRANK, MPI_COMM_WORLD, IERR) 
                call MPI_SEND(innov_at_Grid(mp_start:mp_end), N_sA, mpiReal_size, p_tN,   &
                                          MYRANK*100, MPI_COMM_WORLD, IERR)
        else    !if(p_tRank == 0) then  
                Do pindex =  1, (Np_til - 1)   ! sender proc index within tile group
                        dest_Aoffset = pindex * N_sA + N_sA_Ext + 1   ! dest array offset
                        send_proc = MYRANK +  pindex * MAX_TASKS
                        call MPI_RECV(SNOANL(dest_Aoffset:dest_Aoffset+N_sA-1), N_sA, mpiReal_size, send_proc,      &
                                          send_proc, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                        call MPI_RECV(innov_at_Grid(dest_Aoffset:dest_Aoffset+N_sA-1), N_sA, mpiReal_size, send_proc,      &
                                          send_proc*100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                enddo
        endif
        if (myrank==0) PRINT*,'Finished Data copy'

        if (MYRANK > MAX_TASKS - 1 ) goto 998   ! if(p_tRank /= 0 ) goto 998

        ! avoid -ve anl
        Where(SNOANL < 0.) SNOANL = 0.
        if (print_deb) then
                PRINT*, "Innovation SWE/snwd from rank: ", MYRANK
            PRINT*, innov_at_Grid       
            PRINT*, "Analysis SWE/ snwd  from rank: ", MYRANK
            PRINT*, SNOANL
        endif
        ! d = swe/snd
        ! this is for writing outputs at observation and evaluation points
        if (assim_SWE)  then 
                SWEANL = SNOANL 
                SNDANL = SNOANL / SNODENS_DIST
        else
                SWEANL = SNOANL * SNODENS_DIST
                SNDANL = SNOANL
        endif

        ! !Compute updated snocov 
        ! !Call update_snow_cover_fraction(LENSFC, SNOANL, VETFCS, anl_fSCA)

        ! copy values at eval points
        innov_atEvalPts = IEEE_VALUE(innov_atEvalPts, IEEE_QUIET_NAN)
        SNOANL_atEvalPts = IEEE_VALUE(SNOANL_atEvalPts, IEEE_QUIET_NAN)
        !SNOANL_Cur_atEvalPts = IEEE_VALUE(SNOANL_Cur_atEvalPts, IEEE_QUIET_NAN)
        Do jndx = 1, num_Eval
                if ((index_back_atEval(jndx) > 0) .and. &
                    (LANDMASK(jndx)  == 1 )) then 
                        innov_atEvalPts(jndx) = innov_at_Grid(index_back_atEval(jndx))
                        SNOANL_atEvalPts(jndx) = SNOANL(index_back_atEval(jndx))
                        !SNOANL_Cur_atEvalPts(jndx) = SNOANL_Cur(index_back_atEval(jndx))
                endif
        End do
        ! write outputs 
        Write(rank_str, '(I3.3)') (MYRANK+1)
        ! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/Analysis/"
        da_out_file = "./SNOANLOI."// &  !
                                  TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18_tile"//rank_str//".nc"  !  
        call Write_DA_Outputs(da_out_file, IDIM, JDIM, LENSFC, MYRANK, &
                              SWEFCS, SWEANL, SNDFCS, SNDANL, &  !
                              num_stn, Lat_stn, Lon_stn, SNOOBS_stn, SNOFCS_at_stn, OmB_innov_at_stn,  &                                
                              innov_at_Grid, SNCOV_IMS, &
                              num_Eval, Lat_atEvalPts, Lon_atEvalPts, Obs_atEvalPts, & 
                              SNOFCS_atEvalPts, innov_atEvalPts, SNOANL_atEvalPts) !, SNOANL_Cur_atEvalPts)  !, anl_fSCA) !updated snocov

998 CONTINUE
        ! clean up
             if (allocated(SNOOBS_stn))      DEALLOCATE(SNOOBS_stn)
             if (allocated(SNOFCS_at_stn))   DEALLOCATE(SNOFCS_at_stn)
             if (allocated(OmB_innov_at_stn))   DEALLOCATE(OmB_innov_at_stn)
             if (allocated(Lat_stn))         DEALLOCATE(Lat_stn) 
             if (allocated(Lon_stn))         DEALLOCATE(Lon_stn) 
             if (allocated(OROGFCS_at_stn))  DEALLOCATE(OROGFCS_at_stn) 
             if (allocated(Obs_atEvalPts))   DEALLOCATE(Obs_atEvalPts)
             if (allocated(SNOFCS_atEvalPts)) DEALLOCATE(SNOFCS_atEvalPts)
             if (allocated(innov_atEvalPts))  DEALLOCATE(innov_atEvalPts)
             if (allocated(SNOANL_atEvalPts)) DEALLOCATE(SNOANL_atEvalPts)
             if (allocated(index_back_atObs)) DEALLOCATE(index_back_atObs)
             if (allocated(Lat_atEvalPts))   DEALLOCATE(Lat_atEvalPts) 
             if (allocated(Lon_atEvalPts))   DEALLOCATE(Lon_atEvalPts)
        
        Call UPDATEtime(IY_loc, IM_loc, ID_loc, IH_real, dT_Asssim)
        IH_loc = INT(IH_real)    ! (for the obs. available currently) this should be 18
        if (myrank==0) PRINT*,'Finished OI DA at datetime: ', y_str, m_str, d_str, h_str
        
    End do

999 CONTINUE
        PRINT*,'Finished OI DA ON RANK: ', MYRANK
        CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

        !STOP

        RETURN

 END subroutine Snow_Analysis_OI


 ! SWE threshold for 100% snow cover
 subroutine get_SWE_Threshold(VETFCS_in, SNUP)
        
        IMPLICIT NONE
        !
        Real, Intent(In)        :: VETFCS_in
        Real, Intent(Out)       :: SNUP
        
        INTEGER            :: VETFCS
        REAL               :: snupx(30)

        !This is for the IGBP veg classification scheme.
        ! SWE at which snow cover reaches 100%
        snupx = (/0.080, 0.080, 0.080, 0.080, 0.080, 0.020,     &
                0.020, 0.060, 0.040, 0.020, 0.010, 0.020,                       &
                0.020, 0.020, 0.013, 0.013, 0.010, 0.020,                       &
                0.020, 0.020, 0.000, 0.000, 0.000, 0.000,                       &
                0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)

        VETFCS = INT(VETFCS_in)
        If (VETFCS==0) VETFCS = 7  !vtype_tile[vtype_tile==0] = 7
        
        SNUP = snupx(VETFCS) * 1000.0 ! mm
           
        RETURN
            
 END SUBROUTINE get_SWE_Threshold

 Subroutine update_snow_cover_fraction(LENSFC, SNOANL, VETFCS_in, anl_fSCA)

        IMPLICIT NONE
        !
        include 'mpif.h'
        
        integer, parameter :: dp = kind(1.d0)

        INTEGER, intent(in) :: LENSFC
        REAL, intent(In)   :: SNOANL(LENSFC), VETFCS_in(LENSFC)
        REAL, intent(Out)   ::  anl_fSCA(LENSFC)
        INTEGER                :: VETFCS(LENSFC)

        REAL               :: snupx(30), SNEQV(LENSFC), SNUP, SALP, RSNOW
        Integer                    :: indx, vtype_int

        !This is for the IGBP veg classification scheme.
        snupx = (/0.080, 0.080, 0.080, 0.080, 0.080, 0.020,     &
                        0.020, 0.060, 0.040, 0.020, 0.010, 0.020,                       &
                        0.020, 0.020, 0.013, 0.013, 0.010, 0.020,                       &
                        0.020, 0.020, 0.000, 0.000, 0.000, 0.000,                       &
                        0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)

        SNEQV = 0.001 * SNOANL   ! units mm->m
        SALP = -4.0
        
        VETFCS = INT(VETFCS_in)
        Where(VETFCS==0) VETFCS = 7  !vtype_tile[vtype_tile==0] = 7
        
        Do indx=1, LENSFC
                SNUP = snupx(VETFCS(indx))
                if (SNUP == 0.) then
                        print*, " 0.0 snup value, check vegclasses", vtype_int
                        Stop
                endif

                IF (SNEQV(indx) .LT. SNUP) THEN
                        RSNOW = SNEQV(indx)/SNUP
                        anl_fSCA(indx) = 1. - (EXP(SALP*RSNOW) - RSNOW*EXP(SALP))
                ELSE
                        anl_fSCA(indx) = 1.0
                ENDIF

                if (SNEQV(indx) < 0.00001)  anl_fSCA(indx) = 0.0        

        End do
        
        RETURN

 End Subroutine update_snow_cover_fraction

 ! the following code based on write_data() in read_write_data.f90
 Subroutine Write_DA_Outputs(output_file, idim, jdim, lensfc, myrank,   &
                                 snoforc, snoanl, snwdforc, snwdanal,   &
                                 num_stn, Lat_atObs, Lon_atObs, OBS_stn, FCS_at_stn, OmB_innov_at_stn,  & 
                                 inovatgrid, SNCOV_IMS, &
                                 num_Eval, Lat_atEvalPts, Lon_atEvalPts, Obs_atEvalPts, & 
            SNOFCS_atEvalPts, innov_atEvalPts, SNOANL_atEvalPts) !, SNOANL_Cur_atEvalPts)  !, anl_fSCA) !updated snocov
        !------------------------------------------------------------------
        ! Write DA ouputs: 
        ! forecast SWE
        ! analysis SWE
        ! analysis Snow Depth
        ! innovation at grid
        !------------------------------------------------------------------
        implicit none

        CHARACTER(LEN=*), Intent(In)      :: output_file
        integer, intent(in)         :: idim, jdim, lensfc, num_stn
        real, intent(in)            :: snoforc(lensfc), snoanl(lensfc), snwdforc(lensfc), snwdanal(lensfc)
        Real, intent(in)            :: OBS_stn(num_stn), FCS_at_stn(num_stn), OmB_innov_at_stn(num_stn)
        Real, intent(in)            :: Lat_atObs(num_stn), Lon_atObs(num_stn)
        Real, intent(in)            :: inovatgrid(lensfc), SNCOV_IMS(lensfc)  !, anl_fSCA(lensfc)
        integer, intent(in)         :: num_Eval
        real, intent(in)    :: Lat_atEvalPts(num_Eval), Lon_atEvalPts(num_Eval), Obs_atEvalPts(num_Eval)
        real, intent(in)    :: SNOFCS_atEvalPts(num_Eval), SNOANL_atEvalPts(num_Eval) 
        real, intent(in)    :: innov_atEvalPts(num_Eval)  !, SNOANL_Cur_atEvalPts(num_Eval)

        integer                     :: fsize=65536, inital=0
        integer                     :: header_buffer_val = 16384
        integer                     :: dims_3d(3), dims_strt(3), dims_end(3)
        integer                     :: error, i, ncid
        integer                     :: dim_x, dim_y, dim_time, dim_eval, dim_stn
        integer                     :: id_x, id_y, id_time
        integer       :: id_swe_forc, id_swe, id_snwdf, id_snwd, id_innov, id_imscov   !, id_anlscov
        integer       :: id_latstn, id_lonstn, id_obsstn, id_forcstn, id_innovstn
        integer       :: id_lateval, id_loneval, id_obseval, id_forceval, id_anleval, id_innoveval   !, id_cur_anleval, , id_anlscov
        
        integer                     :: myrank

        real(kind=4)                :: times
        real(kind=4), allocatable   :: x_data(:), y_data(:)
        real(kind=8), allocatable   :: dum2d(:,:)

        include "mpif.h"

        ! print*
        ! print*,"Process ", myrank, "writing output data to: ",trim(output_file)

        !--- create the file
        error = NF90_CREATE(output_file, IOR(NF90_NETCDF4,NF90_CLASSIC_MODEL), ncid, initialsize=inital, chunksize=fsize)
        call netcdf_err(error, 'CREATING FILE='//trim(output_file) )

        !--- define dimensions
        error = nf90_def_dim(ncid, 'xaxis_1', idim, dim_x)
        call netcdf_err(error, 'DEFINING XAXIS DIMENSION' )
        error = nf90_def_dim(ncid, 'yaxis_1', jdim, dim_y)
        call netcdf_err(error, 'DEFINING YAXIS DIMENSION' )
        error = nf90_def_dim(ncid, 'Time', 1, dim_time)
        call netcdf_err(error, 'DEFINING TIME DIMENSION' )
        ! obs points
        error = nf90_def_dim(ncid, 'obs_points', num_stn, dim_stn)
        call netcdf_err(error, 'DEFINING obs_points' )
        ! eval obs points
        if (num_Eval>0) then
                error = nf90_def_dim(ncid, 'eval_points', num_Eval, dim_eval)
                call netcdf_err(error, 'DEFINING eval_points' )
        endif

        !--- define fields
        error = nf90_def_var(ncid, 'xaxis_1', NF90_FLOAT, dim_x, id_x)
        call netcdf_err(error, 'DEFINING XAXIS_1 FIELD' )
        error = nf90_put_att(ncid, id_x, "long_name", "xaxis_1")
        call netcdf_err(error, 'DEFINING XAXIS_1 LONG NAME' )
        error = nf90_put_att(ncid, id_x, "units", "none")
        call netcdf_err(error, 'DEFINING XAXIS_1 UNITS' )
        error = nf90_put_att(ncid, id_x, "cartesian_axis", "X")
        call netcdf_err(error, 'WRITING XAXIS_1 FIELD' )

        error = nf90_def_var(ncid, 'yaxis_1', NF90_FLOAT, dim_y, id_y)
        call netcdf_err(error, 'DEFINING YAXIS_1 FIELD' )
        error = nf90_put_att(ncid, id_y, "long_name", "yaxis_1")
        call netcdf_err(error, 'DEFINING YAXIS_1 LONG NAME' )
        error = nf90_put_att(ncid, id_y, "units", "none")
        call netcdf_err(error, 'DEFINING YAXIS_1 UNITS' )
        error = nf90_put_att(ncid, id_y, "cartesian_axis", "Y")
        call netcdf_err(error, 'WRITING YAXIS_1 FIELD' )

        error = nf90_def_var(ncid, 'Time', NF90_FLOAT, dim_time, id_time)
        call netcdf_err(error, 'DEFINING TIME FIELD' )
        error = nf90_put_att(ncid, id_time, "long_name", "Time")
        call netcdf_err(error, 'DEFINING TIME LONG NAME' )
        error = nf90_put_att(ncid, id_time, "units", "time level")
        call netcdf_err(error, 'DEFINING TIME UNITS' )
        error = nf90_put_att(ncid, id_time, "cartesian_axis", "T")
        call netcdf_err(error, 'WRITING TIME FIELD' )

        dims_3d(1) = dim_x
        dims_3d(2) = dim_y
        dims_3d(3) = dim_time

        error = nf90_def_var(ncid, 'SWE_Forecast', NF90_DOUBLE, dims_3d, id_swe_forc)
        call netcdf_err(error, 'DEFINING SWE_Forecast' )
        error = nf90_put_att(ncid, id_swe_forc, "long_name", "Forecast Snow Water Equivalent")
        call netcdf_err(error, 'DEFINING SWE Forecast LONG NAME' )
        error = nf90_put_att(ncid, id_swe_forc, "units", "mm")
        call netcdf_err(error, 'DEFINING SWE Forecast UNITS' )

        error = nf90_def_var(ncid, 'SWE_Analysis', NF90_DOUBLE, dims_3d, id_swe)
        call netcdf_err(error, 'DEFINING SWE_Analysis' )
        error = nf90_put_att(ncid, id_swe, "long_name", "Analysis Snow Water Equivalent")
        call netcdf_err(error, 'DEFINING SWE LONG NAME' )
        error = nf90_put_att(ncid, id_swe, "units", "mm")
        call netcdf_err(error, 'DEFINING SWE UNITS' )

        error = nf90_def_var(ncid, 'SND_Forecast', NF90_DOUBLE, dims_3d, id_snwdf)
        call netcdf_err(error, 'DEFINING SNF Forecast' )
        error = nf90_put_att(ncid, id_snwdf, "long_name", "Forecast Snow Depth")
        call netcdf_err(error, 'DEFINING SND Forecast LONG NAME' )
        error = nf90_put_att(ncid, id_snwdf, "units", "mm")
        call netcdf_err(error, 'DEFINING SND Forecast UNITS' )

        error = nf90_def_var(ncid, 'SND_Analysis', NF90_DOUBLE, dims_3d, id_snwd)
        call netcdf_err(error, 'DEFINING SND Analyis' )
        error = nf90_put_att(ncid, id_snwd, "long_name", "Analysis Snow Depth")
        call netcdf_err(error, 'DEFINING SND Analysis LONG NAME' )
        error = nf90_put_att(ncid, id_snwd, "units", "mm")
        call netcdf_err(error, 'DEFINING SND Analysis UNITS' )

        error = nf90_def_var(ncid, 'DA_Innovation', NF90_DOUBLE, dims_3d, id_innov)
        call netcdf_err(error, 'DEFINING DA_Innovation' )
        error = nf90_put_att(ncid, id_innov, "long_name", "DA Innovation at model grid")
        call netcdf_err(error, 'DEFINING DA_Innovation LONG NAME' )
        error = nf90_put_att(ncid, id_innov, "units", "mm")
        call netcdf_err(error, 'DEFINING DA_Innovation UNITS' )

        error = nf90_def_var(ncid, 'imsfSCA', NF90_DOUBLE, dims_3d, id_imscov)
        call netcdf_err(error, 'DEFINING imsfSCA' )
        error = nf90_put_att(ncid, id_imscov, "long_name", "IMS fractional Snow Covered Area")
        call netcdf_err(error, 'DEFINING imsfSCA LONG NAME' )
        error = nf90_put_att(ncid, id_imscov, "units", "-")
        call netcdf_err(error, 'DEFINING imsfSCA UNITS' )

        ! obs points
        error = nf90_def_var(ncid, 'latitude@MetaData', NF90_DOUBLE, dim_stn, id_latstn)
        call netcdf_err(error, 'DEFINING  latitude@MetaData' )
        error = nf90_put_att(ncid, id_latstn, "long_name", "Latitude at Observation Points")
        call netcdf_err(error, 'DEFINING  latitude@MetaData LONG NAME' )
        error = nf90_put_att(ncid, id_latstn, "units", "deg")
        call netcdf_err(error, 'DEFINING  latitude@MetaData UNITS' )

        error = nf90_def_var(ncid, 'longitude@MetaData', NF90_DOUBLE, dim_stn, id_lonstn)
        call netcdf_err(error, 'DEFINING longitude@MetaData' )
        error = nf90_put_att(ncid, id_lonstn, "long_name", "Longitude at Observation Points")
        call netcdf_err(error, 'DEFINING longitude@MetaData LONG NAME' )
        error = nf90_put_att(ncid, id_lonstn, "units", "deg")
        call netcdf_err(error, 'DEFINING longitude@MetaData UNITS' )
        
        error = nf90_def_var(ncid, 'snwdph@ObsValue', NF90_DOUBLE, dim_stn, id_obsstn)
        call netcdf_err(error, 'DEFINING snwdph@ObsValue' )
        error = nf90_put_att(ncid, id_obsstn, "long_name", "Observed at Observation Points")
        call netcdf_err(error, 'DEFINING snwdph@ObsValue LONG NAME' )
        error = nf90_put_att(ncid, id_obsstn, "units", "mm")
        call netcdf_err(error, 'DEFINING snwdph@ObsValue UNITS' )
        
        error = nf90_def_var(ncid, 'snwdph@hofx', NF90_DOUBLE, dim_stn, id_forcstn) 
        call netcdf_err(error, 'DEFINING snwdph@hofx' )
        error = nf90_put_att(ncid, id_forcstn, "long_name", "Forecast at Observation Points")
        call netcdf_err(error, 'DEFINING snwdph@hofx LONG NAME' )
        error = nf90_put_att(ncid, id_forcstn, "units", "mm")
        call netcdf_err(error, 'DEFINING snwdph@hofx UNITS' )

        error = nf90_def_var(ncid, 'snwdph@omb', NF90_DOUBLE, dim_stn, id_innovstn)
        call netcdf_err(error, 'DEFINING snwdph@omb' )
        error = nf90_put_att(ncid, id_innovstn, "long_name", "Innovation(O-B) at Observation Points")
        call netcdf_err(error, 'DEFINING snwdph@omb LONG NAME' )
        error = nf90_put_att(ncid, id_innovstn, "units", "mm")
        call netcdf_err(error, 'DEFINING snwdph@omb UNITS' )

        ! error = nf90_def_var(ncid, 'anlfSCA', NF90_DOUBLE, dims_3d, id_anlscov)
        ! call netcdf_err(error, 'DEFINING anlfSCA' )
        ! error = nf90_put_att(ncid, id_anlscov, "long_name", "Analysis fractional Snow Covered Area")
        ! call netcdf_err(error, 'DEFINING anlfSCA LONG NAME' )
        ! error = nf90_put_att(ncid, id_anlscov, "units", "-")
        ! call netcdf_err(error, 'DEFINING anlfSCA UNITS' )

        ! eval points
        if (num_Eval>0) then 
                error = nf90_def_var(ncid, 'LatEvalPoints', NF90_DOUBLE, dim_eval, id_lateval)
                call netcdf_err(error, 'DEFINING LatEvalPoints' )
                error = nf90_put_att(ncid, id_lateval, "long_name", "Latitude at Evaluation Points")
                call netcdf_err(error, 'DEFINING LatEvalPoints LONG NAME' )
                error = nf90_put_att(ncid, id_lateval, "units", "deg")
                call netcdf_err(error, 'DEFINING LatEvalPoints UNITS' )

                error = nf90_def_var(ncid, 'LonEvalPoints', NF90_DOUBLE, dim_eval, id_loneval)
                call netcdf_err(error, 'DEFINING LonEvalPoints' )
                error = nf90_put_att(ncid, id_loneval, "long_name", "Longitude at Evaluation Points")
                call netcdf_err(error, 'DEFINING LonEvalPoints LONG NAME' )
                error = nf90_put_att(ncid, id_loneval, "units", "deg")
                call netcdf_err(error, 'DEFINING LonEvalPoints UNITS' )
                
                error = nf90_def_var(ncid, 'Obs_atEvalPts', NF90_DOUBLE, dim_eval, id_obseval)
                call netcdf_err(error, 'DEFINING Obs_atEvalPts' )
                error = nf90_put_att(ncid, id_obseval, "long_name", "Observed at Evaluation Points")
                call netcdf_err(error, 'DEFINING Obs_atEvalPts LONG NAME' )
                error = nf90_put_att(ncid, id_obseval, "units", "mm")
                call netcdf_err(error, 'DEFINING Obs_atEvalPts UNITS' )
                
                error = nf90_def_var(ncid, 'SNOFCS_atEvalPts', NF90_DOUBLE, dim_eval, id_forceval)
                call netcdf_err(error, 'DEFINING SNOFCS_atEvalPts' )
                error = nf90_put_att(ncid, id_forceval, "long_name", "Forecast at Evaluation Points")
                call netcdf_err(error, 'DEFINING SNOFCS_atEvalPts LONG NAME' )
                error = nf90_put_att(ncid, id_forceval, "units", "mm")
                call netcdf_err(error, 'DEFINING SNOFCS_atEvalPts UNITS' )

                error = nf90_def_var(ncid, 'SNOANL_atEvalPts', NF90_DOUBLE, dim_eval, id_anleval)
                call netcdf_err(error, 'DEFINING SNOANL_atEvalPts' )
                error = nf90_put_att(ncid, id_anleval, "long_name", "Analysis at Evaluation Points")
                call netcdf_err(error, 'DEFINING SNOANL_atEvalPts LONG NAME' )
                error = nf90_put_att(ncid, id_anleval, "units", "mm")
                call netcdf_err(error, 'DEFINING SNOANL_atEvalPts UNITS' )

                ! error = nf90_def_var(ncid, 'SNOANL_Cur_atEvalPts', NF90_DOUBLE, dim_eval, id_cur_anleval)
                ! call netcdf_err(error, 'DEFINING SNOANL_Cur_atEvalPts' )
                ! error = nf90_put_att(ncid, id_cur_anleval, "long_name", "Current Analysis at Evaluation Points")
                ! call netcdf_err(error, 'DEFINING SNOANL_Cur_atEvalPts LONG NAME' )
                ! error = nf90_put_att(ncid, id_cur_anleval, "units", "mm")
                ! call netcdf_err(error, 'DEFINING SNOANL_Cur_atEvalPts UNITS' )
                
                error = nf90_def_var(ncid, 'Innov_atEvalPts', NF90_DOUBLE, dim_eval, id_innoveval)
                call netcdf_err(error, 'DEFINING Innov_atEvalPts' )
                error = nf90_put_att(ncid, id_innoveval, "long_name", "Innovation at Evaluation Points")
                call netcdf_err(error, 'DEFINING Innov_atEvalPts LONG NAME' )
                error = nf90_put_att(ncid, id_innoveval, "units", "mm")
                call netcdf_err(error, 'DEFINING Innov_atEvalPts UNITS' )

                ! error = nf90_def_var(ncid, 'anlfSCA', NF90_DOUBLE, dims_3d, id_anlscov)
                ! call netcdf_err(error, 'DEFINING anlfSCA' )
                ! error = nf90_put_att(ncid, id_anlscov, "long_name", "Analysis fractional Snow Covered Area")
                ! call netcdf_err(error, 'DEFINING anlfSCA LONG NAME' )
                ! error = nf90_put_att(ncid, id_anlscov, "units", "-")
                ! call netcdf_err(error, 'DEFINING anlfSCA UNITS' )
        endif

        error = nf90_enddef(ncid, header_buffer_val,4,0,4)
        call netcdf_err(error, 'DEFINING HEADER' )

        allocate(x_data(idim))
        do i = 1, idim
        x_data(i) = float(i)
        enddo
        allocate(y_data(jdim))
        do i = 1, jdim
        y_data(i) = float(i)
        enddo
        times = 1.0

        error = nf90_put_var( ncid, id_x, x_data)
        call netcdf_err(error, 'WRITING XAXIS RECORD' )
        error = nf90_put_var( ncid, id_y, y_data)
        call netcdf_err(error, 'WRITING YAXIS RECORD' )
        error = nf90_put_var( ncid, id_time, times)
        call netcdf_err(error, 'WRITING TIME RECORD' )

        allocate(dum2d(idim,jdim))
        dims_strt(1:3) = 1
        dims_end(1) = idim
        dims_end(2) = jdim
        dims_end(3) = 1
        
        dum2d = reshape(snoforc, (/idim,jdim/))
        error = nf90_put_var( ncid, id_swe_forc, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING SWE Forecast RECORD' )

        dum2d = reshape(snoanl, (/idim,jdim/))
        error = nf90_put_var( ncid, id_swe, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING SWE Analysis RECORD' ) 

        dum2d = reshape(snwdforc, (/idim,jdim/))
        error = nf90_put_var( ncid, id_snwdf, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING SND Forecast RECORD' )

        dum2d = reshape(snwdanal, (/idim,jdim/))
        error = nf90_put_var( ncid, id_snwd, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING SND analysis RECORD' )

        dum2d = reshape(inovatgrid, (/idim,jdim/))
        error = nf90_put_var( ncid, id_innov, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING innovation RECORD' )

        dum2d = reshape(SNCOV_IMS, (/idim,jdim/))
        error = nf90_put_var( ncid, id_imscov, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING imsfSCA RECORD' )

        ! dum2d = reshape(anl_fSCA, (/idim,jdim/))
        ! error = nf90_put_var( ncid, id_anlscov, dum2d, dims_strt, dims_end)
        ! call netcdf_err(error, 'WRITING anlfSCA RECORD' )
        
        ! obs points (obs, hofx, omb) 
        error = nf90_put_var( ncid, id_latstn, Lat_atObs)
        call netcdf_err(error, 'WRITING Lat_atObsPts RECORD' )

        error = nf90_put_var( ncid, id_lonstn, Lon_atObs)
        call netcdf_err(error, 'WRITING Lon_atObsPts RECORD' )

        error = nf90_put_var( ncid, id_obsstn, OBS_stn)
        call netcdf_err(error, 'WRITING Obs_atObsPts RECORD' )

        error = nf90_put_var( ncid, id_forcstn, FCS_at_stn)
        call netcdf_err(error, 'WRITING SNOFCS_atObsPts RECORD' )

        error = nf90_put_var( ncid, id_innovstn, OmB_innov_at_stn)
        call netcdf_err(error, 'WRITING innov_atObsPts RECORD' )

        ! eval points
        if (num_Eval>0) then 
                error = nf90_put_var( ncid, id_lateval, Lat_atEvalPts)
                call netcdf_err(error, 'WRITING Lat_atEvalPts RECORD' )

                error = nf90_put_var( ncid, id_loneval, Lon_atEvalPts)
                call netcdf_err(error, 'WRITING Lon_atEvalPts RECORD' )

                error = nf90_put_var( ncid, id_obseval, Obs_atEvalPts)
                call netcdf_err(error, 'WRITING Obs_atEvalPts RECORD' )

                error = nf90_put_var( ncid, id_forceval, SNOFCS_atEvalPts)
                call netcdf_err(error, 'WRITING SNOFCS_atEvalPts RECORD' )

                error = nf90_put_var( ncid, id_anleval, SNOANL_atEvalPts)
                call netcdf_err(error, 'WRITING SNOANL_atEvalPts RECORD' )
                ! error = nf90_put_var( ncid, id_cur_anleval, SNOANL_Cur_atEvalPts)
                ! call netcdf_err(error, 'WRITING SNOANL_Cur_atEvalPts RECORD' )
                error = nf90_put_var( ncid, id_innoveval, innov_atEvalPts)
                call netcdf_err(error, 'WRITING innov_atEvalPts RECORD' )
        endif

        deallocate(x_data, y_data)
        deallocate(dum2d)

        error = nf90_close(ncid)
    
 End subroutine Write_DA_Outputs

 SUBROUTINE Observation_Read_atEval(ghcnd_inp_file, &   !dim_name,                      &
                    NDIM, Lat_GHCND, Lon_GHCND, MYRANK)
        
        IMPLICIT NONE
    
        include 'mpif.h'
        !Open netCDF for a snotel and read the SWE, SnowDepth,..., Lat, Lon, at a given datetime
        !ToDO: Can you use variable length char array ?
        CHARACTER(LEN=*), Intent(In)      :: ghcnd_inp_file     !, dim_name
        INTEGER                :: ERROR, NCID
        INTEGER                :: MYRANK
        INTEGER                :: ID_DIM, ID_VAR
        INTEGER, Intent(Out)   :: NDIM
    
        REAL, ALLOCATABLE, Intent(Out)     :: Lat_GHCND(:), Lon_GHCND(:) !, Ele_GHCND(:)
    
        ERROR=NF90_OPEN(TRIM(ghcnd_inp_file),NF90_NOWRITE,NCID)
        CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(ghcnd_inp_file) )
    
        ERROR=NF90_INQ_DIMID(NCID, 'eval_points', ID_DIM)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Dimension' )
    
        ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NDIM)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Dimension' )
    
        ALLOCATE(Lat_GHCND(NDIM))
        ALLOCATE(Lon_GHCND(NDIM))
        !ALLOCATE(Ele_GHCND(NDIM))
    
        ERROR=NF90_INQ_VARID(NCID, 'LatEvalPoints', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lat ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, Lat_GHCND)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lat RECORD' )
    
        ERROR=NF90_INQ_VARID(NCID, 'LonEvalPoints', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lon ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, Lon_GHCND)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lon RECORD' )
    
        ERROR = NF90_CLOSE(NCID)
                  
        RETURN
        
 End SUBROUTINE Observation_Read_atEval

 SUBROUTINE Find_Nearest_GridIndices_Parallel(Myrank, MAX_TASKS, p_tN, p_tRank, Np_til, &
                                            num_src, num_tar, RLA_cg, RLO_cg, RLA, RLO, index_ens_atGrid)
                            
        IMPLICIT NONE
        !
        !USE intrinsic::ieee_arithmetic
        include "mpif.h"
        
        Integer, Intent(In)          :: Myrank, MAX_TASKS, p_tN, p_tRank, Np_til, num_src, num_tar
        Real, Intent(In)                :: RLA_cg(num_src), RLO_cg(num_src), RLA(num_tar), RLO(num_tar)
        Integer, Intent(Out)         :: index_ens_atGrid(num_tar)
        Real                        :: RLA_cg_rad(num_src), RLO_cg_rad(num_src)
        Real                        :: RLA_rad(num_tar), RLO_rad(num_tar)
        
        INTEGER                     :: indx, min_indx
        Real                        :: distArr(num_src), haversinArr(num_src)
        Real                        :: d_latArr(num_src), d_lonArr(num_src)
        Real(16), Parameter         :: PI_16 = 4 * atan (1.0_16)        
        Real(16), Parameter         :: pi_div_180 = PI_16/180.0
        Real, Parameter                 :: earth_rad = 6371.
        
        ! for mpi par
        INTEGER            :: N_sA, N_sA_Ext, mp_start, mp_end 
        INTEGER            :: send_proc, rec_proc, rec_stat(MPI_STATUS_SIZE), dest_Aoffset, pindex
        INTEGER            :: mpiInt_size, isize, IERR

        !Np_til ! num proc. per tile p_tRank ! proc. rank within tile !p_tN  ! tile for proc.
        N_sA = num_tar / Np_til  ! sub array length per proc
        N_sA_Ext = num_tar - N_sA * Np_til ! extra grid cells
        if(p_tRank == 0) then 
                mp_start = 1
        else
                mp_start = p_tRank * N_sA + N_sA_Ext + 1   ! start index of subarray for proc
        endif
        mp_end = (p_tRank + 1) * N_sA + N_sA_Ext                ! end index of subarray for proc
                
        index_ens_atGrid = -1  
        ! at each target point compute its distance from source RLA/RLO pairs and find the position of the minimum      
        RLA_rad =  pi_div_180 * RLA
        RLO_rad =  pi_div_180 * RLO
        RLA_cg_rad =  pi_div_180 * RLA_cg
        RLO_cg_rad =  pi_div_180 * RLO_cg       
        ! https://en.wikipedia.org/wiki/Haversine_formula
        Do indx = mp_start, mp_end    ! num_tar 
                d_latArr = (RLA_rad(indx) - RLA_cg_rad) / 2.
                d_lonArr = (RLO_rad(indx) - RLO_cg_rad) / 2.
                haversinArr = sin(d_latArr)**2 + cos(RLA_rad(indx)) * cos(RLA_cg_rad) * sin(d_lonArr)**2
                WHERE(haversinArr > 1) haversinArr = 1.   ! ensure numerical errors don't make h>1
                distArr = 2 * earth_rad * asin(sqrt(haversinArr))               
                min_indx = MINLOC(distArr, dim = 1)  !, MASK=ieee_is_nan(distArr))
                index_ens_atGrid(indx) = min_indx
        end do

        isize = SIZEOF(N_sA) 
        Call MPI_TYPE_SIZE(MPI_INTEGER, mpiInt_size, IERR) 
        If (isize == 2 ) then 
                mpiInt_size = MPI_INTEGER2
        elseif (isize == 4 ) then 
                mpiInt_size = MPI_INTEGER4
        elseif (isize == 8 ) then 
                mpiInt_size = MPI_INTEGER8
        else
                PRINT*," Possible mismatch between Fortran Int ", isize," and Mpi Int ", mpiInt_size
                Stop
        endif
    
        if (MYRANK > (MAX_TASKS - 1) ) then
                call MPI_SEND(index_ens_atGrid(mp_start:mp_end), N_sA, mpiInt_size, p_tN,   &
                                                MYRANK*1000, MPI_COMM_WORLD, IERR)
        else !if (MYRANK == p_tN ) then  
                Do pindex =  1, (Np_til - 1)   ! sender proc index within tile group
                        dest_Aoffset = pindex * N_sA + N_sA_Ext + 1   ! dest array offset
                        send_proc = MYRANK +  pindex * MAX_TASKS
                        call MPI_RECV(index_ens_atGrid(dest_Aoffset:dest_Aoffset+N_sA-1), N_sA, mpiInt_size, send_proc, &
                                                send_proc*1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                enddo
        endif
    !ToDO: better way to do this?
        ! now share the whole array
        if (MYRANK < MAX_TASKS ) then   !if (MYRANK == p_tN ) then      
                Do pindex =  1, (Np_til - 1)   ! receiving proc index within tile group
                        rec_proc = MYRANK +  pindex * MAX_TASKS
                        call MPI_SEND(index_ens_atGrid, num_tar, mpiInt_size, rec_proc, MYRANK*100, MPI_COMM_WORLD, IERR)
                enddo
        else 
                call MPI_RECV(index_ens_atGrid, num_tar, mpiInt_size, p_tN, p_tN*100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
        endif
             
    RETURN
        
 END SUBROUTINE Find_Nearest_GridIndices_Parallel
     
 SUBROUTINE read_ensemble_forcing(MYRANK, MAX_TASKS, TILE_NUM, ens_size, LENSFC, IDIM_In, JDIM_In, & 
                                        index_ens_atGrid,    &
                                    ens_inp_path, y_str, m_str, d_str, h_str, hprev_str, &
                                    assim_SWE, save_Ens, SNOFCS_Inp_Ens, da_out_file)    ! assim_SWE, 
    
        IMPLICIT NONE

        include "mpif.h"

        INTEGER, INTENT(IN)       :: MYRANK, MAX_TASKS, ens_size, LENSFC        
        CHARACTER(LEN=5), INTENT(In)      :: TILE_NUM
        INTEGER, Intent(In)               :: IDIM_In, JDIM_In
        Integer, Intent(In)                :: index_ens_atGrid(LENSFC)
        !REAL, INTENT(In)                  :: RLA_cg(LENSFC/4), RLO_cg(LENSFC/4), RLA(LENSFC), RLO(LENSFC)
        CHARACTER(LEN=*), Intent(In)      :: ens_inp_path, y_str, m_str, d_str, h_str, hprev_str        
        LOGICAL, Intent(In)                       :: assim_SWE, save_Ens
        REAL, INTENT(OUT)         :: SNOFCS_Inp_Ens(ens_size, LENSFC) !, SNOCOV_Ens(ens_size, LENSFC)  !VEGFCS(LENSFC), 
        !REAL, INTENT(OUT)        :: FMM(LENSFC), FHH(LENSFC), SRFLAG(LENSFC)
        CHARACTER(LEN=*)          :: da_out_file

        CHARACTER(LEN=50)         :: FNBGSI
        CHARACTER(LEN=3)          :: ens_str
        CHARACTER(LEN=250)        :: ens_inp_path_full

        INTEGER                   :: ERROR, NCID
        INTEGER                   :: IDIM, JDIM, ID_DIM, ens_indx       
        INTEGER                   :: i, ix, jy, dim_x, dim_y, dim_ens, id_x, id_y, id_ens, id_forc
        INTEGER                   :: ID_VAR

        REAL(KIND=8)              :: DUMMY(IDIM_In/2, JDIM_In/2), DUMMY2D(IDIM_In, JDIM_In)
        REAL(KIND=8)              :: SNOFCS_Inp_Ens_In(LENSFC/4)
        real(kind=4), allocatable   :: x_data(:), y_data(:), ens_data(:)
        integer                     :: fsize=65536, inital=0
        integer                     :: header_buffer_val = 16384
        integer                     :: dims_3d(3), dims_strt(3), dims_end(3)
        
        Do ens_indx = 1, ens_size
                WRITE(ens_str, '(I3.3)') (ens_indx)
                ens_inp_path_full = TRIM(ens_inp_path)//"mem"//ens_str//"/RESTART/"// &
                                                TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"." // &
                                                TRIM(h_str)// "0000.sfc_data."//TILE_NUM//".nc"

                ERROR=NF90_OPEN(TRIM(ens_inp_path_full),NF90_NOWRITE,NCID)
                CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(ens_inp_path_full) )

                ERROR=NF90_INQ_DIMID(NCID, 'xaxis_1', ID_DIM)
                CALL NETCDF_ERR(ERROR, 'READING xaxis_1' )
                ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=IDIM)
                CALL NETCDF_ERR(ERROR, 'READING xaxis_1' )
                ERROR=NF90_INQ_DIMID(NCID, 'yaxis_1', ID_DIM)
                CALL NETCDF_ERR(ERROR, 'READING yaxis_1' )
                ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=JDIM)
                CALL NETCDF_ERR(ERROR, 'READING yaxis_1' )

                IF ((IDIM*JDIM) /= (LENSFC/4)) THEN
                PRINT*,'FATAL ERROR: DIMENSIONS WRONG.'
                CALL MPI_ABORT(MPI_COMM_WORLD, 88)
                ENDIF
                if(assim_SWE) then
                        ERROR=NF90_INQ_VARID(NCID, "sheleg", ID_VAR)
                        CALL NETCDF_ERR(ERROR, 'READING sheleg ID' )
                        ERROR=NF90_GET_VAR(NCID, ID_VAR, DUMMY)
                        CALL NETCDF_ERR(ERROR, 'READING sheleg' )
                else
                        ERROR=NF90_INQ_VARID(NCID, "snwdph", ID_VAR)
                        CALL NETCDF_ERR(ERROR, 'READING snwdph ID' )
                        ERROR=NF90_GET_VAR(NCID, ID_VAR, DUMMY)
                        CALL NETCDF_ERR(ERROR, 'READING snwdph' )
                endif
                ! Do ix = 1, IDIM_In
                !       Do jy = 1, JDIM_In
                !               DUMMY2D(ix, jy) = DUMMY((ix+1)/2, (jy+1)/2)
                !       End do
                ! End do
                ! SNOFCS_Inp_Ens(ens_indx, :) = RESHAPE(DUMMY2D, (/LENSFC/))
                SNOFCS_Inp_Ens_In = RESHAPE(DUMMY, (/(LENSFC/4)/))
                Do ix = 1, LENSFC
                        SNOFCS_Inp_Ens(ens_indx, ix) = SNOFCS_Inp_Ens_In(index_ens_atGrid(ix))
                End do          
                ! ERROR=NF90_INQ_VARID(NCID, "sncovr", ID_VAR)
                ! CALL NETCDF_ERR(ERROR, 'READING sncovr ID' )
                ! ERROR=NF90_GET_VAR(NCID, ID_VAR, DUMMY)
                ! CALL NETCDF_ERR(ERROR, 'READING sncovr' )
                ! Do ix = 1, IDIM_In
                !       Do jy = 1, JDIM_In
                !               DUMMY2D(ix, jy) = DUMMY((ix+1)/2, (jy+1)/2)
                !       End do
                ! End do
                ! SNOCOV_Ens(ens_indx, :) = RESHAPE(DUMMY2D, (/LENSFC/))

                ERROR = NF90_CLOSE(NCID)
        End do
        
        If (save_Ens .and. (Myrank < MAX_TASKS)) then
        
                print*,"Process ", myrank, "writing ens snow forc data to: ",trim(da_out_file)
                
                !--- create the file
                error = NF90_CREATE(da_out_file, IOR(NF90_NETCDF4,NF90_CLASSIC_MODEL), &
                                                                        ncid, initialsize=inital, chunksize=fsize)
                call netcdf_err(error, 'CREATING FILE='//trim(da_out_file) )

                !--- define dimensions
                error = nf90_def_dim(ncid, 'xaxis_1', IDIM_In, dim_x)
                call netcdf_err(error, 'DEFINING XAXIS DIMENSION' )
                error = nf90_def_dim(ncid, 'yaxis_1', JDIM_In, dim_y)
                call netcdf_err(error, 'DEFINING YAXIS DIMENSION' )
                error = nf90_def_dim(ncid, 'Ens_mem', ens_size, dim_ens)
                call netcdf_err(error, 'DEFINING ENS DIMENSION' )

                !--- define fields
                error = nf90_def_var(ncid, 'xaxis_1', NF90_FLOAT, dim_x, id_x)
                call netcdf_err(error, 'DEFINING XAXIS_1 FIELD' )
                error = nf90_put_att(ncid, id_x, "long_name", "xaxis_1")
                call netcdf_err(error, 'DEFINING XAXIS_1 LONG NAME' )
                error = nf90_put_att(ncid, id_x, "units", "none")
                call netcdf_err(error, 'DEFINING XAXIS_1 UNITS' )
                error = nf90_put_att(ncid, id_x, "cartesian_axis", "X")
                call netcdf_err(error, 'WRITING XAXIS_1 FIELD' )

                error = nf90_def_var(ncid, 'yaxis_1', NF90_FLOAT, dim_y, id_y)
                call netcdf_err(error, 'DEFINING YAXIS_1 FIELD' )
                error = nf90_put_att(ncid, id_y, "long_name", "yaxis_1")
                call netcdf_err(error, 'DEFINING YAXIS_1 LONG NAME' )
                error = nf90_put_att(ncid, id_y, "units", "none")
                call netcdf_err(error, 'DEFINING YAXIS_1 UNITS' )
                error = nf90_put_att(ncid, id_y, "cartesian_axis", "Y")
                call netcdf_err(error, 'WRITING YAXIS_1 FIELD' )

                error = nf90_def_var(ncid, 'Ens_mem', NF90_FLOAT, dim_ens, id_ens)
                call netcdf_err(error, 'DEFINING ENS FIELD' )
                error = nf90_put_att(ncid, id_ens, "long_name", "Ensemble members")
                call netcdf_err(error, 'DEFINING ENS LONG NAME' )
                error = nf90_put_att(ncid, id_ens, "units", "-")
                call netcdf_err(error, 'DEFINING ENS UNITS' )
                error = nf90_put_att(ncid, id_ens, "cartesian_axis", "ENS")
                call netcdf_err(error, 'WRITING ENS FIELD' )

                dims_3d(1) = dim_x
                dims_3d(2) = dim_y
                dims_3d(3) = dim_ens

                error = nf90_def_var(ncid, 'SND_Forecast', NF90_DOUBLE, dims_3d, id_forc)
                call netcdf_err(error, 'DEFINING SND_Forecast' )
                error = nf90_put_att(ncid, id_forc, "long_name", "Forecast Snow Depth")
                call netcdf_err(error, 'DEFINING SND Forecast LONG NAME' )
                error = nf90_put_att(ncid, id_forc, "units", "mm")
                call netcdf_err(error, 'DEFINING SND Forecast UNITS' )

                error = nf90_enddef(ncid, header_buffer_val,4,0,4)
                call netcdf_err(error, 'DEFINING HEADER' )

                allocate(x_data(IDIM_In))
                do i = 1, IDIM_In
                x_data(i) = float(i)
                enddo
                allocate(y_data(JDIM_In))
                do i = 1, JDIM_In
                y_data(i) = float(i)
                enddo
                allocate(ens_data(ens_size))
                do i = 1, ens_size
                        ens_data(i) = float(i)
                enddo

                error = nf90_put_var( ncid, id_x, x_data)
                call netcdf_err(error, 'WRITING XAXIS RECORD' )
                error = nf90_put_var( ncid, id_y, y_data)
                call netcdf_err(error, 'WRITING YAXIS RECORD' )
                error = nf90_put_var( ncid, id_ens, ens_data)
                call netcdf_err(error, 'WRITING Ensmem RECORD' )

                Do ens_indx = 1, ens_size
                        dims_strt(1:2) = 1
                        dims_strt(3) = ens_indx
                        dims_end(1) = IDIM_In
                        dims_end(2) = JDIM_In
                        dims_end(3) = 1
                        
                        DUMMY2D = reshape(SNOFCS_Inp_Ens(ens_indx, :), (/IDIM_In,JDIM_In/))
                        error = nf90_put_var( ncid, id_forc, DUMMY2D, dims_strt, dims_end)
                        call netcdf_err(error, 'WRITING SWE Forecast RECORD' )
                End do
                error = nf90_close(ncid)
                deallocate(x_data, y_data, ens_data)
        Endif

        !deallocate(DUMMY, DUMMY2D)
        
 END SUBROUTINE read_ensemble_forcing
    
 SUBROUTINE READ_Forecast_Data(MYRANK, LENSFC, veg_type_landice, SWEFCS, SNDFCS, VETFCS, LANDMASK)
    
        IMPLICIT NONE

        include "mpif.h"

        INTEGER, INTENT(IN)       :: MYRANK, LENSFC, veg_type_landice
        REAL, INTENT(OUT)         :: SWEFCS(LENSFC),SNDFCS(LENSFC),VETFCS(LENSFC)
        INTEGER, INTENT(OUT)      :: LANDMASK(LENSFC) 

        CHARACTER(LEN=50)         :: FNBGSI
        CHARACTER(LEN=3)          :: RANKCH

        INTEGER                   :: ERROR, NCID
        INTEGER                   :: IDIM, JDIM, ID_DIM
        INTEGER                   :: ID_VAR, i

        REAL(KIND=8), ALLOCATABLE :: DUMMY(:,:)
        REAL                      :: SLMASK(LENSFC)

        !CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, ERROR)

        WRITE(RANKCH, '(I3.3)') (MYRANK+1)

        FNBGSI = "./fnbgsi." // RANKCH
        if (print_deb) PRINT*, "READ INPUT SFC DATA FROM: "//TRIM(FNBGSI)

        ERROR=NF90_OPEN(TRIM(FNBGSI),NF90_NOWRITE,NCID)
        CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(FNBGSI) )

        ERROR=NF90_INQ_DIMID(NCID, 'xaxis_1', ID_DIM)
        CALL NETCDF_ERR(ERROR, 'READING xaxis_1' )
        ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=IDIM)
        CALL NETCDF_ERR(ERROR, 'READING xaxis_1' )

        ERROR=NF90_INQ_DIMID(NCID, 'yaxis_1', ID_DIM)
        CALL NETCDF_ERR(ERROR, 'READING yaxis_1' )
        ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=JDIM)
        CALL NETCDF_ERR(ERROR, 'READING yaxis_1' )

        IF ((IDIM*JDIM) /= LENSFC) THEN
        PRINT*,'FATAL ERROR: DIMENSIONS WRONG.'
        CALL MPI_ABORT(MPI_COMM_WORLD, 88)
        ENDIF

        ALLOCATE(DUMMY(IDIM,JDIM))

        ERROR=NF90_INQ_VARID(NCID, "sheleg", ID_VAR)
        CALL NETCDF_ERR(ERROR, 'READING sheleg ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING sheleg' )
        SWEFCS = RESHAPE(DUMMY, (/LENSFC/)) 

        ERROR=NF90_INQ_VARID(NCID, "snwdph", ID_VAR)
        CALL NETCDF_ERR(ERROR, 'READING snwdph ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING snwdph' )
        SNDFCS = RESHAPE(DUMMY, (/LENSFC/))

        ERROR=NF90_INQ_VARID(NCID, "vtype", ID_VAR)
        CALL NETCDF_ERR(ERROR, 'READING vtype ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING vtype' )
        VETFCS = RESHAPE(DUMMY, (/LENSFC/))    

        ERROR=NF90_INQ_VARID(NCID, "slmsk", ID_VAR)
        CALL NETCDF_ERR(ERROR, 'READING slmsk ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING slmsk' )
        SLMASK = RESHAPE(DUMMY, (/LENSFC/))    

        do i = 1, LENSFC 
           ! if land, but not land ice, set mask to 1.
           if ( (NINT(SLMASK(i)) == 1 ) .and.   & 
                ( NINT(VETFCS(i)) /=  veg_type_landice  )) then 
                LANDMASK(i) = 1 
           else 
                LANDMASK(i) = 0
           endif 
        enddo

        ! slmsk is 0 - ocean, 1 - land, 2 -seaice 
        ! convert to integer  0 not land or glacier, 1 - non-glacier covered land
        
        DEALLOCATE(DUMMY)

        ERROR = NF90_CLOSE(NCID)
    
 END SUBROUTINE READ_Forecast_Data

 SUBROUTINE READ_Forecast_Data_atPath(forc_inp_path, veg_type_landice, LENSFC, SWEFCS, SNDFCS, VETFCS, LANDMASK) !VEGFCS, !SRFLAG)
    
        IMPLICIT NONE

        include "mpif.h"
        
        CHARACTER(LEN=*), Intent(In)      :: forc_inp_path
        INTEGER, INTENT(IN)               :: LENSFC, veg_type_landice
        REAL, INTENT(OUT)                 :: SWEFCS(LENSFC), SNDFCS(LENSFC), VETFCS(LENSFC)  !VEGFCS(LENSFC), 
        INTEGER, INTENT(OUT)              :: LANDMASK(LENSFC) 
        !REAL, INTENT(OUT)        :: FMM(LENSFC), FHH(LENSFC), SRFLAG(LENSFC)

        CHARACTER(LEN=50)         :: FNBGSI
        CHARACTER(LEN=3)          :: RANKCH

        INTEGER                   :: ERROR, NCID, i
        INTEGER                   :: IDIM, JDIM, ID_DIM
        INTEGER                   :: ID_VAR

        REAL(KIND=8), ALLOCATABLE :: DUMMY(:,:)
        REAL                      :: SLMASK(LENSFC)

        !CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, ERROR)
        ! WRITE(RANKCH, '(I3.3)') (MYRANK+1)
        ! FNBGSI = "./fnbgsi." // RANKCH
        ! if (print_deb) PRINT*, "READ INPUT SFC DATA FROM: "//TRIM(FNBGSI)

        !CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, ERROR)
        ! if (print_deb) PRINT*, "READ INPUT SFC DATA FROM: "//TRIM(forc_inp_path)

        ERROR=NF90_OPEN(TRIM(forc_inp_path), NF90_NOWRITE,NCID)
        CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(forc_inp_path) )
        ! ERROR=NF90_OPEN(TRIM(FNBGSI),NF90_NOWRITE,NCID)
        ! CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(FNBGSI) )

        ERROR=NF90_INQ_DIMID(NCID, 'xaxis_1', ID_DIM)
        CALL NETCDF_ERR(ERROR, 'READING xaxis_1' )
        ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=IDIM)
        CALL NETCDF_ERR(ERROR, 'READING xaxis_1' )

        ERROR=NF90_INQ_DIMID(NCID, 'yaxis_1', ID_DIM)
        CALL NETCDF_ERR(ERROR, 'READING yaxis_1' )
        ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=JDIM)
        CALL NETCDF_ERR(ERROR, 'READING yaxis_1' )

        IF ((IDIM*JDIM) /= LENSFC) THEN
        PRINT*,'FATAL ERROR: DIMENSIONS WRONG.'
        CALL MPI_ABORT(MPI_COMM_WORLD, 88)
        ENDIF

        ALLOCATE(DUMMY(IDIM,JDIM))

        ERROR=NF90_INQ_VARID(NCID, "sheleg", ID_VAR)
        CALL NETCDF_ERR(ERROR, 'READING sheleg ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING sheleg' )
        SWEFCS = RESHAPE(DUMMY, (/LENSFC/))

        ERROR=NF90_INQ_VARID(NCID, "snwdph", ID_VAR)
        CALL NETCDF_ERR(ERROR, 'READING snwdph ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING snwdph' )
        SNDFCS = RESHAPE(DUMMY, (/LENSFC/))

        ERROR=NF90_INQ_VARID(NCID, "vtype", ID_VAR)
        CALL NETCDF_ERR(ERROR, 'READING vtype ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING vtype' )
        VETFCS = RESHAPE(DUMMY, (/LENSFC/))    

        ERROR=NF90_INQ_VARID(NCID, "slmsk", ID_VAR)
        CALL NETCDF_ERR(ERROR, 'READING slmsk ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING slmsk' )
        SLMASK = RESHAPE(DUMMY, (/LENSFC/))    

        do i = 1, LENSFC 
           ! if land, but not land ice, set mask to 1.
           if ( (NINT(SLMASK(i)) == 1 ) .and.   & 
                ( NINT(VETFCS(i)) /=  veg_type_landice  )) then 
                LANDMASK(i) = 1 
           else 
                LANDMASK(i) = 0
           endif 
        enddo

        ! slmsk is 0 - ocean, 1 - land, 2 -seaice 
        ! convert to integer  0 not land or glacier, 1 - non-glacier covered land
        DEALLOCATE(DUMMY)

        ERROR = NF90_CLOSE(NCID)
    
 END SUBROUTINE READ_Forecast_Data_atPath

 SUBROUTINE READ_Analysis_Data(anl_inp_path, MYRANK, LENSFC, SWEFCS, SNDFCS) 
    
        IMPLICIT NONE

        include "mpif.h"
        
        CHARACTER(LEN=*), Intent(In)      :: anl_inp_path
        INTEGER, INTENT(IN)       :: MYRANK, LENSFC
        REAL, INTENT(OUT)         :: SWEFCS(LENSFC), SNDFCS(LENSFC)

        INTEGER                   :: ERROR, NCID
        INTEGER                   :: IDIM, JDIM, ID_DIM
        INTEGER                   :: ID_VAR

        REAL(KIND=8), ALLOCATABLE :: DUMMY(:,:)

        !CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, ERROR)
        if (print_deb) PRINT*, "READ INPUT SFC DATA FROM: "//TRIM(anl_inp_path)

        ERROR=NF90_OPEN(TRIM(anl_inp_path),NF90_NOWRITE,NCID)
        CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(anl_inp_path) )

        ERROR=NF90_INQ_DIMID(NCID, 'xaxis_1', ID_DIM)
        CALL NETCDF_ERR(ERROR, 'READING xaxis_1' )
        ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=IDIM)
        CALL NETCDF_ERR(ERROR, 'READING xaxis_1' )

        ERROR=NF90_INQ_DIMID(NCID, 'yaxis_1', ID_DIM)
        CALL NETCDF_ERR(ERROR, 'READING yaxis_1' )
        ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=JDIM)
        CALL NETCDF_ERR(ERROR, 'READING yaxis_1' )

        IF ((IDIM*JDIM) /= LENSFC) THEN
        PRINT*,'FATAL ERROR: DIMENSIONS WRONG.'
        CALL MPI_ABORT(MPI_COMM_WORLD, 88)
        ENDIF

        ALLOCATE(DUMMY(IDIM, JDIM))

        ERROR=NF90_INQ_VARID(NCID, "sheleg", ID_VAR)
        CALL NETCDF_ERR(ERROR, 'READING sheleg ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING sheleg' )
        SWEFCS = RESHAPE(DUMMY, (/LENSFC/))

        ERROR=NF90_INQ_VARID(NCID, "snwdph", ID_VAR)
        CALL NETCDF_ERR(ERROR, 'READING snwdph ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING snwdph' )
        SNDFCS = RESHAPE(DUMMY, (/LENSFC/)) 
 
        
        DEALLOCATE(DUMMY)

        ERROR = NF90_CLOSE(NCID)
    
 END SUBROUTINE READ_Analysis_Data
    
 SUBROUTINE READ_LAT_LON_CoarseRes(inp_path,RLA,RLO,IDIM,JDIM,IJDIM)
    
     IMPLICIT NONE
    
     include "mpif.h"
    
     INTEGER, INTENT(IN)    :: IDIM, JDIM, IJDIM
         
         CHARACTER(LEN=*), Intent(In)      :: inp_path
    
     REAL, INTENT(OUT)      :: RLA(IJDIM),RLO(IJDIM)
    
     INTEGER                :: ERROR, NCID
     INTEGER                :: I, II, J, JJ
     INTEGER                :: ID_DIM, ID_VAR, NX, NY
    
     REAL, ALLOCATABLE      :: DUMMY(:,:), GEOLAT(:,:), GEOLON(:,:)
    
     !CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, ERROR)
    
     if (print_deb) then
        PRINT*, "READ FV3 GRID INFO FROM: "//TRIM(inp_path)
     endif
    
     ERROR=NF90_OPEN(TRIM(inp_path),NF90_NOWRITE,NCID)
     CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(inp_path) )
    
     ERROR=NF90_INQ_DIMID(NCID, 'nx', ID_DIM)
     CALL NETCDF_ERR(ERROR, 'ERROR READING NX ID' )
    
     ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NX)
     CALL NETCDF_ERR(ERROR, 'ERROR READING NX' )
    
     ERROR=NF90_INQ_DIMID(NCID, 'ny', ID_DIM)
     CALL NETCDF_ERR(ERROR, 'ERROR READING NY ID' )
    
     ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NY)
     CALL NETCDF_ERR(ERROR, 'ERROR READING NY' )
    
     IF ((NX/2) /= IDIM .OR. (NY/2) /= JDIM) THEN
       PRINT*,'FATAL ERROR: DIMENSIONS IN FILE: ',(NX/2),(NY/2)
       PRINT*,'DO NOT MATCH GRID DIMENSIONS: ',IDIM,JDIM
       CALL MPI_ABORT(MPI_COMM_WORLD, 130)
     ENDIF
    
     ALLOCATE(GEOLON(NX+1,NY+1))
     ALLOCATE(GEOLAT(NX+1,NY+1))
    
     ERROR=NF90_INQ_VARID(NCID, 'x', ID_VAR)
     CALL NETCDF_ERR(ERROR, 'ERROR READING X ID' )
     ERROR=NF90_GET_VAR(NCID, ID_VAR, GEOLON)
     CALL NETCDF_ERR(ERROR, 'ERROR READING X RECORD' )
    
     ERROR=NF90_INQ_VARID(NCID, 'y', ID_VAR)
     CALL NETCDF_ERR(ERROR, 'ERROR READING Y ID' )
     ERROR=NF90_GET_VAR(NCID, ID_VAR, GEOLAT)
     CALL NETCDF_ERR(ERROR, 'ERROR READING Y RECORD' )
    
     ALLOCATE(DUMMY(IDIM,JDIM))
    
     DO J = 1, JDIM
       DO I = 1, IDIM
         II = 2*I
         JJ = 2*J
         DUMMY(I,J) = GEOLON(II,JJ)
       ENDDO
     ENDDO
    
     RLO = RESHAPE(DUMMY, (/IJDIM/))
    
     DEALLOCATE(GEOLON)
    
     DO J = 1, JDIM
       DO I = 1, IDIM
         II = 2*I
         JJ = 2*J
         DUMMY(I,J) = GEOLAT(II,JJ)
       ENDDO
     ENDDO
    
     RLA = RESHAPE(DUMMY, (/IJDIM/))
    
     DEALLOCATE(GEOLAT, DUMMY)
    
     ERROR = NF90_CLOSE(NCID)
    
 END SUBROUTINE READ_LAT_LON_CoarseRes

 SUBROUTINE READ_LAT_LON_OROG_atRank(MYRANK, RLA,RLO,OROG,TILE_NUM,IDIM,JDIM,IJDIM) !OROG_UF,
    
    !--------------------------------------------------------------
    ! READ LATITUDE, LONGITUDE, FILTERED OROGRAPHY, AND
    ! UNFILTERED OROGRAPHY FOR THE CUBED-SPHERE TILE FROM
    ! THE "GRID" FILE.
    !--------------------------------------------------------------
    
     IMPLICIT NONE
    
     include "mpif.h"
    
     INTEGER, INTENT(IN)    :: IDIM, JDIM, IJDIM
    
     CHARACTER(LEN=5), INTENT(OUT) :: TILE_NUM
    
     REAL, INTENT(OUT)      :: RLA(IJDIM),RLO(IJDIM)
     REAL, INTENT(OUT)      :: OROG(IJDIM)   !,OROG_UF(IJDIM)
    
     CHARACTER(LEN=50)      :: FNOROG, FNGRID
     CHARACTER(LEN=3)       :: RANKCH
    
     INTEGER                :: ERROR, NCID, NCID_OROG
     INTEGER                :: I, II, J, JJ, MYRANK
     INTEGER                :: ID_DIM, ID_VAR, NX, NY
    
     REAL, ALLOCATABLE         :: DUMMY(:,:), GEOLAT(:,:), GEOLON(:,:)
     REAL(KIND=4), ALLOCATABLE :: DUMMY4(:,:)
    
     !CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, ERROR)
    
     WRITE(RANKCH, '(I3.3)') (MYRANK+1)
    
     FNGRID = "./fngrid." // RANKCH
    
     if (print_deb) then
        PRINT*, "READ FV3 GRID INFO FROM: "//TRIM(FNGRID)
     endif
    
     ERROR=NF90_OPEN(TRIM(FNGRID),NF90_NOWRITE,NCID)
     CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(FNGRID) )
    
     ERROR=NF90_INQ_DIMID(NCID, 'nx', ID_DIM)
     CALL NETCDF_ERR(ERROR, 'ERROR READING NX ID' )
    
     ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NX)
     CALL NETCDF_ERR(ERROR, 'ERROR READING NX' )
    
     ERROR=NF90_INQ_DIMID(NCID, 'ny', ID_DIM)
     CALL NETCDF_ERR(ERROR, 'ERROR READING NY ID' )
    
     ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NY)
     CALL NETCDF_ERR(ERROR, 'ERROR READING NY' )
    
     IF ((NX/2) /= IDIM .OR. (NY/2) /= JDIM) THEN
       PRINT*,'FATAL ERROR: DIMENSIONS IN FILE: ',(NX/2),(NY/2)
       PRINT*,'DO NOT MATCH GRID DIMENSIONS: ',IDIM,JDIM
       CALL MPI_ABORT(MPI_COMM_WORLD, 130)
     ENDIF
    
     ALLOCATE(GEOLON(NX+1,NY+1))
     ALLOCATE(GEOLAT(NX+1,NY+1))
    
     ERROR=NF90_INQ_VARID(NCID, 'x', ID_VAR)
     CALL NETCDF_ERR(ERROR, 'ERROR READING X ID' )
     ERROR=NF90_GET_VAR(NCID, ID_VAR, GEOLON)
     CALL NETCDF_ERR(ERROR, 'ERROR READING X RECORD' )
    
     ERROR=NF90_INQ_VARID(NCID, 'y', ID_VAR)
     CALL NETCDF_ERR(ERROR, 'ERROR READING Y ID' )
     ERROR=NF90_GET_VAR(NCID, ID_VAR, GEOLAT)
     CALL NETCDF_ERR(ERROR, 'ERROR READING Y RECORD' )
    
     ALLOCATE(DUMMY(IDIM,JDIM))
    
     DO J = 1, JDIM
       DO I = 1, IDIM
         II = 2*I
         JJ = 2*J
         DUMMY(I,J) = GEOLON(II,JJ)
       ENDDO
     ENDDO
    
     RLO = RESHAPE(DUMMY, (/IJDIM/))
    
     DEALLOCATE(GEOLON)
    
     DO J = 1, JDIM
       DO I = 1, IDIM
         II = 2*I
         JJ = 2*J
         DUMMY(I,J) = GEOLAT(II,JJ)
       ENDDO
     ENDDO
    
     RLA = RESHAPE(DUMMY, (/IJDIM/))
    
     DEALLOCATE(GEOLAT, DUMMY)
    
     ERROR=NF90_INQ_VARID(NCID, 'tile', ID_VAR)
     CALL NETCDF_ERR(ERROR, 'ERROR READING TILE ID' )
     ERROR=NF90_GET_VAR(NCID, ID_VAR, TILE_NUM)
     CALL NETCDF_ERR(ERROR, 'ERROR READING TILE RECORD' )
    
     ERROR = NF90_CLOSE(NCID)
    
     FNOROG = "./fnorog." // RANKCH
    
     if (print_deb) PRINT*, "READ FV3 OROG INFO FROM: "//TRIM(FNOROG)
    
     ERROR=NF90_OPEN(TRIM(FNOROG),NF90_NOWRITE,NCID_OROG)
     CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(FNOROG) )
    
     ALLOCATE(DUMMY4(IDIM,JDIM))
    
!      ERROR=NF90_INQ_VARID(NCID_OROG, 'orog_raw', ID_VAR)
!      CALL NETCDF_ERR(ERROR, 'ERROR READING orog_raw ID' )
!      ERROR=NF90_GET_VAR(NCID_OROG, ID_VAR, DUMMY4)
!      CALL NETCDF_ERR(ERROR, 'ERROR READING orog_raw RECORD' )
!      OROG_UF = RESHAPE(DUMMY4, (/IJDIM/))
    
     ERROR=NF90_INQ_VARID(NCID_OROG, 'orog_filt', ID_VAR)
     CALL NETCDF_ERR(ERROR, 'ERROR READING orog_filt ID' )
     ERROR=NF90_GET_VAR(NCID_OROG, ID_VAR, DUMMY4)
     CALL NETCDF_ERR(ERROR, 'ERROR READING orog_filt RECORD' )
     OROG = RESHAPE(DUMMY4, (/IJDIM/))
    
     DEALLOCATE(DUMMY4)
    
     ERROR = NF90_CLOSE(NCID_OROG)
    
 END SUBROUTINE READ_LAT_LON_OROG_atRank

 ! 10.14.20: subroutine copied from UEBFortran at
 ! https://github.com/dtarb/UEBFortran/blob/master/Model/snowxv.f90
 !                 Update time for each time step
 SUBROUTINE UPDATEtime(YEAR,MONTH,DAY,HOUR,DT)
        
        IMPLICIT NONE

        INTEGER   :: YEAR, MONTH, DAY, DMON(12), DM, I       ! 30/03/2004 ITB 
        !INTEGER   :: LYEAR  ! 30/03/2004 ITB  
        Real      :: hour, dt  ! DGT Dec 10, 2004.  Fixing ITB errors 
 
        DATA (DMON(I),I=1,12)/31,28,31,30,31,30,31,31,30,31,30,31/
        HOUR=HOUR+DT
        DM=DMON(MONTH)
!  check for leap years 
        if(month .eq. 2)dm=lyear(year)
10   continue
        IF(HOUR.GE.24.0) THEN
          HOUR=HOUR-24.0
          DAY=DAY+1
          go to 10
        ENDIF
20   continue
        IF(DAY.GT.DM) THEN
          DAY=day - dm
          MONTH=MONTH+1
          IF(MONTH.GT.12) THEN
                MONTH=1
                YEAR=YEAR+1
                  DM=DMON(MONTH)
                if(month .eq. 2)dm=lyear(year)
                endif
                go to 20
        ENDIF
        RETURN

 END SUBROUTINE UPDATEtime

 ! 10.14.20: subroutine copied from UEBFortran at
 ! https://github.com/dtarb/UEBFortran/blob/master/Model/snowxv.f90
!    function to return number of days in February checking for leap years
 function lyear(year) result(Alyear)  
        
        IMPLICIT NONE

        Integer  :: year, Alyear
        IF(MOD(YEAR,4).GT.0 .or. &
                (mod(year,100) .eq.0 .and. mod(year,400) .ne. 0)) THEN
        ! Leap years are every 4 years 
        ! - except for years that are multiples of centuries (e.g. 1800, 1900)
        ! - except again that when the century is divisible by 4 (e.g. 1600, 2000)
        !   then it is a leap year 
                Alyear=28
        ELSE
                Alyear=29
        ENDIF
          
 end function lyear

! subroutine Snow_Analysis_OI_multiTime(SNOW_OI_TYPE, MAX_TASKS, MYRANK, NPROCS, IDIM, JDIM, IY, IM, ID, IH, DELTSFC, num_steps, &  
! !                                                         LENSFC, IVEGSRC, PERCENT_OBS_WITHHELD, &
! !                                                         L_horz , h_ver, obs_tolerance, ims_max_ele, , &
! !                                                         GHCND_SNOWDEPTH_PATH, IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH, &
! !                                                         CURRENT_ANALYSIS_PATH,   &
! !                                                         SNOANL)
                                                        
! !         !----------------------------------------------------------------------
!         ! Input arguments: 
!         ! IDIM * JDIM = LENSFC: number of grid cells in tile = xdim * ydim   
!         ! IY, IM, ID, IH = year, month, day, hour of current model step   
!         ! MYRANK: rank/id of the MPI process
!         ! ...
!         !
!         ! Inputs, read from file:
!         ! RLA, RLO: lat lon information for the tile
!         ! SNOFCS(LENSFC): forecast  snowdepth, or snow water equivalent forecat, as requested
!         !
!         ! Outputs:
!         ! SNOANL: snow analysis (of SWE or SND, see below)
!         ! 
!         ! Draper - changes to snow variables to clarify names, removed unnecesary conversions
!         !          SWE - snow water equivalent 
!         !          SND - snow depth 
!         !          SNO - generic variable to refer to either (determined by value of 
!         !                snow_OI_type
!         !                   
!         !----------------------------------------------------------------------
!         IMPLICIT NONE
!         !
!         include 'mpif.h'
        
!         integer, parameter :: dp = kind(1.d0)

!         INTEGER, intent(in) :: SNOW_OI_TYPE, MAX_TASKS, MYRANK, NPROCS, IDIM, JDIM, IY, IM, ID, IH, LENSFC, IVEGSRC     
!         CHARACTER(LEN=*), Intent(In)   :: GHCND_SNOWDEPTH_PATH, IMS_SNOWCOVER_PATH, & 
!                                                                           IMS_INDEXES_PATH              !, CURRENT_ANALYSIS_PATH
! !  GHCND_SNOWDEPTH_PATH="/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/GHCND/",
! !  IMS_SNOWCOVER_PATH="/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS/",
! !  IMS_INDEXES_PATH="/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS_INDEXES/",
!         REAL, intent(In)    :: PERCENT_OBS_WITHHELD
!         Real, intent(In)    :: L_horz , h_ver, obs_tolerance, ims_max_ele
!         INTEGER, intent(in) :: max_numobs_atgrid 
!         REAL, intent(Out)   :: SNOANL(LENSFC) 
!         CHARACTER(LEN=5)    :: TILE_NUM
!         Character(LEN=3)    :: rank_str
!         INTEGER                     :: IERR     
!         REAL        :: RLA(LENSFC), RLO(LENSFC), RLO_Tile(LENSFC), OROG(LENSFC)  !, OROG_UF(LENSFC)
!         REAL                :: SNDFCS(LENSFC), SWEFCS(LENSFC), SNDANL(LENSFC), SWEANL(LENSFC)
!         REAL                :: SNOFCS(LENSFC), VETFCS(LENSFC), SNUP
!         INTEGER             :: LANDMASK(LENSFC)
!         CHARACTER(len=250)   :: dim_name, ghcnd_inp_file, ims_inp_file, ims_inp_file_indices
!         CHARACTER(len=5)     :: y_str, m_str, d_Str, h_str, fvs_tile
!         REAL, ALLOCATABLE    :: SNOOBS_stn(:), SNOFCS_at_stn(:)               
!         REAL, ALLOCATABLE    :: Lat_stn(:), Lon_stn(:), OROGFCS_at_stn(:)  
!         REAL                 :: lat_min, lat_max, lon_min, lon_max      
!         Real                 :: SNCOV_IMS(LENSFC)  ! ims resampled at each grid
!         Real                 :: SNO_IMS_at_Grid(LENSFC)

!         INTEGER :: num_stn, Num_Ims, num_Eval, !num_subgrd_ims_cels !Num_Ims_Lat, Num_Ims_Lon
!         Real    :: bkgst_srch_rad   ! radius_of_influence for selecting state at observation point
!         INTEGER :: jndx, zndx, ncol, nrow
!         Integer, Allocatable   :: loc_near_Obs(:), index_back_atIMS(:)
!         Integer                :: num_loc, num_loc_1, num_loc_2
!         Real, Parameter         :: Stdev_back_depth = 30., Stdev_Obsv_depth = 40., Stdev_Obsv_ims = 60. !80. ! mm 
!         real                    :: stdev_obsv, stdev_back
!         Integer                         :: ims_assm_hour
!         !Real                           :: obs_tolerance, ims_max_ele

!         Real(dp), Allocatable    :: B_cov_mat(:,:), b_cov_vect(:)
!         Real(dp), Allocatable    :: O_cov_mat(:,:), W_wght_vect(:)
!         Real, Allocatable   :: back_at_Obs(:), obs_Array(:), Lat_Obs(:), Lon_Obs(:), orogfcs_Obs(:)
!         REAL                :: innov_at_Grid(LENSFC)
!         Real, Allocatable   :: obs_Innov(:), OmB_innov_at_stn(:)

!         CHARACTER(len=250)       :: da_out_file  !anl_inp_path, 

!         REAL, ALLOCATABLE  :: SNOFCS_atEvalPts(:), innov_atEvalPts(:), SNOANL_atEvalPts(:)  !evalution points 
!         REAL, ALLOCATABLE  :: Lat_atEvalPts(:), Lon_atEvalPts(:), Obs_atEvalPts(:)     !evalution points
!         Integer, ALLOCATABLE    :: index_back_atEval(:)     ! background locations at eval points 
!         Integer, ALLOCATABLE    :: index_back_atObs(:)   ! the location of background corresponding obs

!         Real                       :: snodens, SNODENS_DIST(LENSFC)
!         LOGICAL                    :: assim_SWE, assim_GHCND   !assimilate swe? (instead of snow depth)
!         LOGICAL                    :: assim_IMS, assim_IMS_thisGridCell    !assimilate sncov, 

!         Integer                    :: veg_type_landice  ! 10.21.20: no assmn over land ice
!     ! for mpi par
!         INTEGER            :: Np_ext, Np_til, p_tN, p_tRank, N_sA, N_sA_Ext, mp_start, mp_end
!         INTEGER            :: send_proc, rec_stat(MPI_STATUS_SIZE), dest_Aoffset, pindex
!         INTEGER            :: mpiReal_size, rsize
!         CHARACTER(len=4)   :: stn_var ! should probably be called control_var
!         REAL :: tmp

! !=============================================================================================
! ! 1. initialise vars,set-up processors, and read lat/lon from orog files.
! !=============================================================================================

!         If (SNOW_OI_TYPE == 2) then 
!                 assim_SWE = .True.  ! note: if this is set true, need to establish stdev_obsv and stdev_back
!                                     !       and need to code readers for appropriate data sets
!                 stn_var = 'SWE' 
!                 assim_IMS = .True.  ! Draper, didn't test this, but not important
!                 assim_GHCND = .False. ! Doesn't have SWE obs.
!         ELSEIF (SNOW_OI_TYPE == 1) then 
!                 assim_SWE = .False. 
!                 stn_var = 'SND' 
!                 assim_GHCND = .True.
!                 stdev_obsv = stdev_obsv_depth
!                 stdev_back = stdev_back_depth
!                 assim_IMS = .True.  
!         ENDIF
!         num_stn = 0 ! controls calling of obs operator for stn data. If remains 0 will not be called.
!         num_Eval = 0

!         ! noah models specific? Needed to ID glaciers.
!         if (IVEGSRC == 2) then   ! sib
!                 veg_type_landice=13
!         else
!                 veg_type_landice=15
!         endif

! ! CSD  - do we need some checks in here to make sure we have enough processors? 
! ! perhaps also write out any extra processors
!         IF (myrank ==0) PRINT*,"Total num proc ", NPROCS, " Num tiles /Max.tasks: ", MAX_TASKS

!         Np_ext = MOD(NPROCS, MAX_TASKS)  ! extra/inactive procs
!         if (MYRANK >  NPROCS - Np_ext - 1) goto 999
!         Np_til = NPROCS / MAX_TASKS  ! num proc. per tile 
!         p_tN = MOD(MYRANK, MAX_TASKS)  ! tile for proc.
!         p_tRank = MYRANK / MAX_TASKS  ! proc. rank within tile
!         N_sA = LENSFC / Np_til  ! sub array length per proc
!         N_sA_Ext = LENSFC - N_sA * Np_til ! extra grid cells
!         if(p_tRank == 0) then 
!                 mp_start = 1
!         else
!                 mp_start = p_tRank * N_sA + N_sA_Ext + 1   ! start index of subarray for proc
!         endif
!         mp_end = (p_tRank + 1) * N_sA + N_sA_Ext                ! end index of subarray for proc        
!         If(myrank == 0 )PRINT*,"sub array length ", N_sA, " extra sub array: ", N_sA_Ext
! if (p_tN /= 4 ) goto 999
!         if (GHCND_SNOWDEPTH_PATH(1:8).eq.'        ') assim_GHCND = .false.
!         if (IMS_SNOWCOVER_PATH(1:8).eq.'        ') assim_IMS = .false.
!         if( (.not. assim_GHCND) .and. (.not. assim_IMS) ) then
!                 print*, "Observation paths don't exist!, skipping the OI DA"
!                 goto 999
!         end if
   
! ! ToDO: use file location the way the rest of the program accesses
!         ! 4.9.20 for now hard code the file location 
!         write(y_str, "(I4)") IY
!         write(m_str, "(I0.2)") IM
!         write(d_str, "(I0.2)") ID
!         write(h_str, "(I0.2)") IH
!         write(fvs_tile, "(I3)") IDIM    

! ! CSD - can we delete orog_uf? 
!     ! READ THE OROGRAPHY AND GRID POINT LAT/LONS FOR THE CUBED-SPHERE TILE p_tN
!         CALL READ_LAT_LON_OROG_atRank(p_tN, RLA,RLO,OROG,TILE_NUM,IDIM,JDIM,LENSFC) !OROG_UF,
!         PRINT*,"Snow anl on ", MYRANK, " Tile group: ", p_tN, " Tile: ", TILE_NUM
        
!         RLO_Tile = RLO ! CSD, why copy?
!         Where(RLO_Tile > 180) RLO_Tile = RLO_Tile - 360
!         lat_min = MAX(MINVAL(RLA) - 1., -90.) ! CSD - why is the limit 1?
!         lat_max = MIN(MAXVAL(RLA) + 1., 90.)
!         lon_min = MAX(MINVAL(RLO_Tile) - 1., -180.)
!         lon_max = MIN(MAXVAL(RLO_Tile) + 1., 180.)      
!         if (p_tN==3)  then  
!                 lon_min = 125.       ! lon_min is left, lon_max is right, not necessary min/max value
!                 lon_max = -145.
!         endif   
!         if ((p_tRank==0) .and. print_deb) then
!                 print*, TILE_NUM, " min/max lat/lon ", lat_min, lat_max, lon_min, lon_max
!         endif

! !=============================================================================================
! ! 3a. Read model forecast here, as need VETFCS and snow density for IMS snow depth calc. (later, separate read routines) 
! !=============================================================================================

!        ! READ THE INPUT SURFACE DATA ON THE CUBED-SPHERE TILE p_tN. A
!        ! Also get vegtype (VETFCS) to identify glacier. 
        
!         Call READ_Forecast_Data(p_tN, LENSFC, veg_type_landice, SWEFCS, SNDFCS, & 
!                                  VETFCS, LANDMASK)

!         if (assim_SWE) then 
!            SNOFCS=SWEFCS  
!         else 
!            SNOFCS=SNDFCS 
!         endif
!         ! average snow density of snow over non-glacier land.
!         SNODENS_DIST = SWEFCS/SNDFCS
!         snodens = SUM(SWEFCS/SNDFCS, Mask = (LANDMASK==1 .and. SNDFCS>0.01 )) &
!                          / COUNT (LANDMASK==1 .and. SNDFCS> 0.01) 
!         ! for grid cells with no valid density, fill in the average snodens
!         Where((LANDMASK==1 .and. SNDFCS>0.01 )) SNODENS_DIST = snodens
!         If (p_tRank==0)  print*, "Tile ", p_tN, ' mean snow density', snodens
!         tmp = SUM(SWEFCS,  Mask = (LANDMASK==1 .and. SNDFCS>0.01 )) &
!                          / COUNT (LANDMASK==1 .and. SNDFCS> 0.01)
!         If (p_tRank==0)  print*, "Tile ", p_tN,  ' mean SWE', tmp
!         tmp = SUM(SNDFCS,  Mask = (LANDMASK==1 .and. SNDFCS>0.01 )) &
!                          / COUNT (LANDMASK==1 .and. SNDFCS> 0.01)
!         If (p_tRank==0)  print*, "Tile ", p_tN,  ' mean SND', tmp
        
! !=============================================================================================
! ! 2. Read observations
! !=============================================================================================

! ! 2a. Read station obs (of snow depth or SWE)

!         if (assim_GHCND) then 
!              ghcnd_inp_file = TRIM(GHCND_SNOWDEPTH_PATH)//"GHCND.SNWD."// &
!                                          TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//TRIM(h_str)//".nc"
!              dim_name = "Site_Id"    
!              call Observation_Read_GHCND_Tile_excNaN(p_tN, ghcnd_inp_file, dim_name, &
!                         lat_min, lat_max, lon_min, lon_max, & 
!                         num_stn, SNOOBS_stn,              &
!                         Lat_stn, Lon_stn, MYRANK) 

!             if ((p_tRank==0) .and. print_deb) then
!                     print*, "Tile ", p_tN, " num. GHCND obs ", num_stn
!             endif
!             if ((p_tRank==0) .and. (p_tN==2) .and. print_deb) then
!                     PRINT*, "GHCND SND from rank: ", MYRANK
!                     PRINT*, SNOOBS_stn
!                     PRINT*, "Lat at GHCND from rank: ", MYRANK
!                     PRINT*, Lat_stn
!                     PRINT*, "Lon at GHCND from rank: ", MYRANK
!                     PRINT*, Lon_stn
!                     PRINT*,'Finished reading GHCND'
!             endif
!             if (myrank==0) PRINT*,'Finished reading GHCND'

!        endif

! ! beyond this point, there should be no specific mention of the station data source

! ! 2b. Read remotely sensed snow cover, and convert to  snow depth or SWE. 

!        if (assim_IMS) then
!             ! (max) number of IMS subcells within a tile grid cell
!             If (IDIM == 96) then   
!                 num_subgrd_ims_cels = 627               
!                 bkgst_srch_rad = 240.     !Km distance from gridcell to search for corresponding background state
!             elseif (IDIM == 128) then
!                num_subgrd_ims_cels = 627               
!                bkgst_srch_rad = 240.          ! CSD  - this should be in a namelist. What is num_sub?
!             elseif (IDIM == 768) then
!                 num_subgrd_ims_cels = 30
!                 bkgst_srch_rad = 27.                              !Km  
!             else
!                 PRINT*,'Error, tile res. not known ', IDIM
!                 stop   
!             endif
              
!             ! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS/"
!             ims_inp_file = TRIM(IMS_SNOWCOVER_PATH)//"IMS.SNCOV."// &
!                                        TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//TRIM(h_str)//".nc"                      !
!             ! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS_INDEXES/"
!             ims_inp_file_indices = TRIM(IMS_INDEXES_PATH)//"C"//TRIM(ADJUSTL(fvs_tile))//&
!                                                        ".IMS.Indices."//TRIM(TILE_NUM)//".nc"                       
!             Call Observation_Read_IMS_Full(ims_inp_file, ims_inp_file_indices, &
!                                                        MYRANK, JDIM, IDIM, num_subgrd_ims_cels, SNCOV_IMS)
!             if((p_tRank==0) .and. (p_tN==1) .and. print_deb) then
!                     PRINT*, "IMS SNCOV from rank: ", MYRANK
!                     PRINT*, SNCOV_IMS
!                     PRINT*,'Finished reading IMS'
!             endif 
!             if (myrank==0) PRINT*,'Finished reading IMS, converting to snow depth' 

!             call CalcSWEFromSnowCover(SNCOV_IMS, VETFCS, LENSFC, SNO_IMS_at_Grid)
! ! CSDCSD - here  filter SNO_IMS_at_Grid, to remove obs where both model and IMS show snow
!             if (.not. assim_SWE)  SNO_IMS_at_Grid = SNO_IMS_at_Grid/SNODENS_DIST !snodens ! convert SWE to SND

!             if ((p_tN==2) .and. (p_tRank==0) .and. print_deb) then
!                     PRINT*, "IMS obs at each grid cell from rank: ", MYRANK
!                     PRINT*, SNO_IMS_at_Grid
!                     PRINT*,'Finished comverting IMS observations'
!             endif
!             if (myrank==0) PRINT*,'Finished converting IMS observations'

!         endif ! read_IMS 

! !=============================================================================================
! ! 3. Get H(x): Read forecast snow fields from restart files, then interpolate to obs location.
! !=============================================================================================

! ! 3a. read the forecast file on model grid : this was done earlier, as need veg type and 
! !      snow density for IMS snow depth conversion

! ! 3b. get H(x) for station obs
!         ! Get model states at obs points
!         if (num_stn > 0) then ! skip if not reading in station data / no obs were available
!             ALLOCATE(SNOFCS_at_stn(num_stn))
!             ALLOCATE(OROGFCS_at_stn(num_stn)) 
!             ALLOCATE(index_back_atObs(num_stn)) 
!             ALLOCATE(OmB_innov_at_stn(num_stn)) 
!              ! using PERCENT_OBS_WITHHELD % of stn locations for evaluation
! ! CSD - todo, add code to skip all eval sets if num_Eval is 0.
!             num_Eval = floor(0.01 * PERCENT_OBS_WITHHELD * num_stn)  
!             if (num_Eval > 0) then 
!                 ALLOCATE(index_back_atEval(num_Eval)) 
!                 ALLOCATE(Obs_atEvalPts(num_Eval)) 
!                 ALLOCATE(SNOFCS_atEvalPts(num_Eval)) 
!                 ALLOCATE(Lat_atEvalPts(num_Eval))
!                 ALLOCATE(Lon_atEvalPts(num_Eval)) 
!                 ALLOCATE(innov_atEvalPts(num_Eval))
!                 ALLOCATE(SNOANL_atEvalPts(num_Eval))    
!                 if(p_tRank == 0) then 
!                         PRINT*, "Tile ", p_tN+1, " ", num_Eval, ' points for evaluation excluded from DA'       
!                 endif  
!             endif 

!     ! CSD todo: separate out eval call and obs operator call. model info (SNOOBS_stn shouldn't be in the obs operator)

! ! CSDCSD - could either add IMS in here, or take advantage of the obs being located on model grid cells to simplify 
!             Call Observation_Operator_Parallel(Myrank, MAX_TASKS, p_tN, p_tRank, Np_til, & 
!                                 RLA, RLO, OROG, Lat_stn, Lon_stn,   &
!                                 LENSFC, num_stn, num_Eval, bkgst_srch_rad, SNOFCS, SNOOBS_stn,  &
!                                 SNOFCS_at_stn, OROGFCS_at_stn, index_back_atObs, index_back_atEval, &
!                                 Obs_atEvalPts, SNOFCS_atEvalPts, Lat_atEvalPts, Lon_atEvalPts)
            
!             OmB_innov_at_stn = SNOOBS_stn - SNOFCS_at_stn
!             if ((p_tN==2) .and. (p_tRank==0) .and. print_deb) then 
!                     PRINT*, "Background Indices at eval points"
!                     PRINT*, index_back_atEval       
!                     PRINT*, "Obs at Eval Points" 
!                     PRINT*, Obs_atEvalPts   
!                     PRINT*, "Forecast at Eval Points"
!                     PRINT*, SNOFCS_atEvalPts
!                     PRINT*, "O - B (innovation at obs points)"
!                     PRINT*, OmB_innov_at_stn          
!                     PRINT*, "Lat at Eval Points"
!                     PRINT*, Lat_atEvalPts
!                     PRINT*, "Lon at Eval Points"
!                     PRINT*, Lon_atEvalPts
!             endif

!             if ((p_tN==2) .and. (p_tRank==0) .and. print_deb) then
!                     PRINT*, "station Lat range from rank: ", MYRANK, MINVAL(Lat_stn), " ", MAXVAL(Lat_stn)
!                     PRINT*, "station Lon range from rank: ", MYRANK, MINVAL(Lon_stn), " ", MAXVAL(Lon_stn)
!                     PRINT*, "Model elevation at station locations from rank: ", MYRANK
!                     PRINT*, OROGFCS_at_stn
!                     PRINT*, "Background ", stn_var, " at station locations from rank: ", MYRANK
!                     PRINT*, SNOFCS_at_stn  
!             endif
!             if (myrank==0) PRINT*,'Finished observation operator for station data'         
!         endif ! num_stn > 0

! !=============================================================================================
! ! 4.  Obs-based obs QC goes here
! !=============================================================================================

! !CSDCSD - todo. Add  QC here.

! ! if model? elevation >1500, discard stn_obs
! ! if abs (model - obs ) elevation > ? m discard obs 
! ! QC obs (or model background might be easier?) with the LANDMASK
! ! min/max limits on obs? 

! ! CSDCSD - probably here, construct an obs vector that contains both IMS and station obs
! ! this then replaces the obs_stn input below (you'll have to be able to track the obs source for R too)

! !=============================================================================================
! ! 5. Perform the DA update, by looping over each grid cell 
! !=============================================================================================

!         ims_assm_hour = 18  
!         obs_srch_rad = 250. ! radius of observation search

!         if (myrank==0) PRINT*,'Starting DA loop'
!         Do jndx = mp_start, mp_end     
!                 num_loc_1 = 0
!                 num_loc_2 = 0
!                 assim_IMS_thisGridCell = .FALSE.
!                 call debug_print("loop ", float(jndx))
!                 if(num_stn>0) then 
!                 ! currently: find station obs in radius, do gross error check, and limit to 50 obs
!                         call nearest_Observations_Locations(RLA(jndx), RLO(jndx),    &
!                                         Lat_stn, Lon_stn,  num_stn, obs_srch_rad, max_numobs_atgrid,   &
!                                         stdev_back, stdev_obsv, obs_tolerance,                 &
!                                         SNOFCS_at_stn, SNOOBS_stn,                                              &
!                                         loc_near_Obs,  num_loc_1) !,      &LENSFC,
!                         call debug_print("number of stn sndpth obs ", float(num_loc))
!                 endif                         
!                 if ( assim_IMS .AND. (IH == ims_assm_hour) .AND. & 
!                      (OROG(jndx) <= ims_max_ele) .AND. &
!                      (.NOT. IEEE_IS_NAN(SNO_IMS_at_Grid(jndx))) ) then
!                         Call get_SWE_Threshold(VETFCS(jndx), SNUP)
!                         If ((SWEFCS(jndx) >= SNUP) .AND. &    !no IMS assmn: if model & ims both 100% snow covered 
!                             (SNO_IMS_at_Grid(jndx) >= SNUP) ) then              
!                                 assim_IMS_thisGridCell = .False.
!                                 num_loc_2 = 0
!                         else                      
!                                 assim_IMS_thisGridCell = .TRUE.
!                                 num_loc_2 = 1
!                         endif
!                 endif
!                 ! CSDCSD - this is only assimilating the IMS obs at the current grid cell. 
!                 ! want to use all within radius ( I bet this is why you were adding snow where 
!                 ! wgere IMS doesn't have it.
!                 !check IMS is assimilated
!                 ! if (myrank==1) print*, "loop ", jndx
!                 ! if ( assim_IMS .AND. (IH == ims_assm_hour) .AND. & 
!                 !      (OROG(jndx) <= ims_max_ele) .AND. &
!                 !      (.NOT. IEEE_IS_NAN(SNO_IMS_at_Grid(jndx))) ) then

!                 !         Call get_SWE_Threshold(VETFCS(jndx), SNUP)
!                 !         If ((SWEFCS(jndx) >= SNUP) .AND. &    !no IMS assmn: if model & ims both 100% snow covered 
!                 !             (SNO_IMS_at_Grid(jndx) >= SNUP) ) then              
!                 !                 assim_IMS_thisGridCell = .False.
!                 !                 num_loc_2 = 0
!                 !         else
!                 !                 Call nearest_IMS_Locations(LENSFC, 10, obs_srch_rad, &
!                 !                         ims_max_ele, RLA(jndx), RLO(jndx), RLA, RLO, OROG, & 
!                 !                         SNOFCS, SNO_IMS_at_Grid, num_loc_2, index_back_atIMS)
!                 !                 if (num_loc_2 > 0) then 
!                 !                         assim_IMS_thisGridCell = .TRUE.
!                 !                         if ((print_deb) .AND. &
!                 !                             (myrank==1) ) print*, " number of IMS obs: ", num_loc_2                        
!                 !                 end if
!                 !         end if                        
!                 ! endif     
!                 num_loc = num_loc_1 + num_loc_2
!                 ! if assim_IMS=false >> num_loc_1=num_loc
!                 if((num_loc > 0) .and. ( LANDMASK(jndx) == 1 )) then ! .and. (SNCOV_IMS(jndx) > 0.)) then    
!                         ! get background states
!                         Allocate(back_at_Obs(num_loc))
!                         Allocate(obs_Array(num_loc))
!                         Allocate(Lat_Obs(num_loc))
!                         Allocate(Lon_Obs(num_loc))
!                         Allocate(orogfcs_Obs(num_loc))
!                         ! snow depth/swe 
!                         if(num_loc_1 > 0) then
!                                 Do zndx = 1, num_loc_1     
!                                         back_at_Obs(zndx) = SNOFCS_at_stn(loc_near_Obs(zndx))
!                                         obs_Array(zndx) = SNOOBS_stn(loc_near_Obs(zndx))
!                                         Lat_Obs(zndx) = Lat_stn(loc_near_Obs(zndx))
!                                         Lon_Obs(zndx) = Lon_stn(loc_near_Obs(zndx))
!                                         orogfcs_Obs(zndx) = OROGFCS_at_stn(loc_near_Obs(zndx)) 
!                                 End Do
!                         End if
!                         !ims
!                         if(assim_IMS_thisGridCell) then
!                                 back_at_Obs(num_loc) = SNOFCS(jndx)
! 				obs_Array(num_loc) = SNO_IMS_at_Grid(jndx)
! 				Lat_Obs(num_loc) = RLA(jndx)   !Lat_IMS_atGrid(jndx)
! 				Lon_Obs(num_loc) = RLO(jndx)   !Lon_IMS_atGrid(jndx)
!                                 orogfcs_Obs(num_loc) = OROG(jndx)  !Ele_IMS(jndx)
!                                 !multi-IMS
!                                 ! Do zndx = 1, num_loc_2
!                                 !         back_at_Obs(num_loc_1+zndx) = SNOFCS(index_back_atIMS(zndx))
!                                 !         obs_Array(num_loc_1+zndx) = SNO_IMS_at_Grid(index_back_atIMS(zndx))
!                                 !         Lat_Obs(num_loc_1+zndx) = RLA(index_back_atIMS(zndx))   !Lat_IMS_atGrid(jndx)
!                                 !         Lon_Obs(num_loc_1+zndx) = RLO(index_back_atIMS(zndx))   !Lon_IMS_atGrid(jndx)
!                                 !         orogfcs_Obs(num_loc_1+zndx) = OROG(index_back_atIMS(zndx))  !Ele_IMS(jndx)
!                                 ! End Do
!                         endif
!                         ! compute covariances
!                         Allocate(B_cov_mat(num_loc, num_loc))
!                         Allocate(b_cov_vect(num_loc))
!                         Allocate(O_cov_mat(num_loc, num_loc))
!                         Allocate(W_wght_vect(num_loc))   
!                         call compute_covariances(RLA(jndx), RLO(jndx), OROG(jndx), SNOFCS(jndx),    &
!                                 Lat_Obs, Lon_Obs, orogfcs_Obs, num_loc,    			 &
!                                 Stdev_back, Stdev_Obsv, Stdev_Obsv_ims,      &
!                                 L_horz, h_ver,                                   &   !L_horz in Km, h_ver in m
!                                 assim_IMS_thisGridCell,                          &
!                                 B_cov_mat, b_cov_vect, O_cov_mat, W_wght_vect)

!                         ! CSDCSD - will need to change R spec in this function, now that have > 1 ims obs.
!                         ! call compute_covariances_multiIMS(RLA(jndx), RLO(jndx), OROG(jndx), SNOFCS(jndx),    &
!                         !         Lat_Obs, Lon_Obs, orogfcs_Obs, num_loc, num_loc_1, num_loc_2,                &
!                         !         Stdev_back, Stdev_Obsv, Stdev_Obsv_ims,      &
!                         !         L_horz, h_ver,                                   &   !L_horz in Km, h_ver in m
!                         !         assim_IMS_thisGridCell,                          &
!                         !         B_cov_mat, b_cov_vect, O_cov_mat, W_wght_vect)
!                         ! call OI DA
!                         Allocate(obs_Innov(num_loc))
!                         call Snow_DA_OI(back_at_Obs, obs_Array, num_loc, W_wght_vect,            &
!                                 SNOFCS(jndx), innov_at_Grid(jndx), SNOANL(jndx), obs_Innov)
!                         if ((p_tN==1) .and. (p_tRank==0) .and. print_deb) then  
!                                 print*, "proc ", myrank, "loop ", jndx, "num depth obs ", num_loc_1, "total obs", num_loc
!                                 PRINT*, " background at obs pts: "
!                                 PRINT*, back_at_Obs     
!                                 PRINT*, "Observed"
!                                 PRINT*,  obs_Array
!                                 PRINT*, "Obs innovation: "
!                                 PRINT*, obs_Innov
!                                 PRINT*, "Weight vector: "
!                                 PRINT*, W_wght_vect     
!                                 print*, "innov: ", innov_at_Grid(jndx), "forec: ", SNOFCS(jndx), " anl: ", SNOANL(jndx)
!                         endif           
!                         !free mem
!                         DEALLOCATE(back_at_Obs, obs_Array)
!                         DEALLOCATE(Lat_Obs, Lon_Obs, orogfcs_Obs, obs_Innov)
!                         DEALLOCATE(B_cov_mat, b_cov_vect, O_cov_mat, W_wght_vect)
!                         ! QCC by ims--use ims as snow mask
!                         if((SNCOV_IMS(jndx) >= 0.5) .and. (SNOANL(jndx) < 50.)) then
!                                 SNOANL(jndx) = 50.
!                         ! elseif((SNCOV_IMS(jndx) < 0.5) .and. (SNCOV_IMS(jndx) > 0.1) .and. (anl_at_Grid(jndx) >= 50)) then
!                         !       anl_at_Grid(jndx) = 50.
!                         elseif((SNCOV_IMS(jndx) < 0.1) .and. (SNOANL(jndx) >= 0.)) then
!                                         SNOANL(jndx) = 0.
!                         endif
!                         if ((p_tN==1) .and. (p_tRank==0) .and. print_deb) then                                          
!                                 PRINT*, "analyis at grid after mask: ", SNOANL(jndx)
!                         endif
!                 else
!                         SNOANL(jndx) = SNOFCS(jndx) ! if not land, or no obs avail keep forecast
!                 endif
!                 if (allocated(loc_near_Obs))  Deallocate(loc_near_Obs) 
!                 if (allocated(index_back_atIMS))  Deallocate(index_back_atIMS)                 
!         End do
!         CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
!         if (myrank==0) PRINT*, ' Finished DA loops'
        

! ! ToDO: Better way to handle this? ! CSD - I'll fix this later.
! ! Real data type size corresponding to mpi
!         rsize = SIZEOF(snodens)
!         Call MPI_TYPE_SIZE(MPI_REAL, mpiReal_size, IERR) 
!         If (rsize == 4 ) then 
!                 mpiReal_size = MPI_REAL4
!         elseif (rsize == 8 ) then 
!                 mpiReal_size = MPI_REAL8
!         elseif (rsize == 16 ) then 
!                 mpiReal_size = MPI_REAL16
!         else
!                 PRINT*," Possible mismatch between Fortran Real ", rsize," and Mpi Real ", mpiReal_size
!                 Stop
!         endif
!         ! send analyses arrays to 'tile-level root' proc.               
!         if (MYRANK > (MAX_TASKS - 1) ) then
!                 call MPI_SEND(SNOANL(mp_start:mp_end), N_sA, mpiReal_size, p_tN,   &
!                                           MYRANK, MPI_COMM_WORLD, IERR) 
!                 call MPI_SEND(innov_at_Grid(mp_start:mp_end), N_sA, mpiReal_size, p_tN,   &
!                                           MYRANK*100, MPI_COMM_WORLD, IERR)
!         else    !if(p_tRank == 0) then  
!                 Do pindex =  1, (Np_til - 1)   ! sender proc index within tile group
!                         dest_Aoffset = pindex * N_sA + N_sA_Ext + 1   ! dest array offset
!                         send_proc = MYRANK +  pindex * MAX_TASKS
!                         call MPI_RECV(SNOANL(dest_Aoffset:dest_Aoffset+N_sA-1), N_sA, mpiReal_size, send_proc,      &
!                                           send_proc, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
!                         call MPI_RECV(innov_at_Grid(dest_Aoffset:dest_Aoffset+N_sA-1), N_sA, mpiReal_size, send_proc,      &
!                                           send_proc*100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
!                 enddo
!         endif
!         if (myrank==0) PRINT*,'Finished Data copy'

!         if (MYRANK > MAX_TASKS - 1 ) goto 998   ! if(p_tRank /= 0 ) goto 998

!         ! avoid -ve anl
!         Where(SNOANL < 0.) SNOANL = 0.
!         if (print_deb) then
!                 PRINT*, "Innovation SWE/snwd from rank: ", MYRANK
!             PRINT*, innov_at_Grid       
!             PRINT*, "Analysis SWE/ snwd  from rank: ", MYRANK
!             PRINT*, SNOANL
!         endif
!         ! d = swe/snd
!         ! this is for writing outputs at observation and evaluation points
!         if (assim_SWE)  then 
!                 SWEANL = SNOANL 
!                 SNDANL = SNOANL / SNODENS_DIST
!         else
!                 SWEANL = SNOANL * SNODENS_DIST
!                 SNDANL = SNOANL
!         endif

!         ! !Compute updated snocov 
!         ! !Call update_snow_cover_fraction(LENSFC, SNOANL, VETFCS, anl_fSCA)

!         ! copy values at eval points
!         innov_atEvalPts = IEEE_VALUE(innov_atEvalPts, IEEE_QUIET_NAN)
!         SNOANL_atEvalPts = IEEE_VALUE(SNOANL_atEvalPts, IEEE_QUIET_NAN)
!         !SNOANL_Cur_atEvalPts = IEEE_VALUE(SNOANL_Cur_atEvalPts, IEEE_QUIET_NAN)
!         Do jndx = 1, num_Eval
!                 if ((index_back_atEval(jndx) > 0) .and. &
!                     (LANDMASK(jndx)  == 1 )) then 
!                         innov_atEvalPts(jndx) = innov_at_Grid(index_back_atEval(jndx))
!                         SNOANL_atEvalPts(jndx) = SNOANL(index_back_atEval(jndx))
!                         !SNOANL_Cur_atEvalPts(jndx) = SNOANL_Cur(index_back_atEval(jndx))
!                 endif
!         End do
!         ! Do jndx = 1, Num_Ghcnd
!         !       if ((index_back_atObs(jndx) > 0) .and. &
!         !           (NINT(VETFCS(jndx)) /= veg_type_landice)) then  !10.21.20: exclude if land-ice      
!         !                       SNOANL_atEvalPts(jndx) = anl_at_Grid(index_back_atObs(jndx))
!         !                       SNOANL_Cur_atEvalPts(jndx) = SNOANL_Cur(index_back_atObs(jndx))
!         !       endif
!         ! End do

!         ! write outputs 
!         Write(rank_str, '(I3.3)') (MYRANK+1)
!         ! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/Analysis/"
!         da_out_file = "./SNOANLOI."// &  !
!                                   TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18_tile"//rank_str//".nc"  !  
!         call Write_DA_Outputs(da_out_file, IDIM, JDIM, LENSFC, MYRANK, &
!                               SWEFCS, SWEANL, SNDFCS, SNDANL, &  !
!                               num_stn, Lat_stn, Lon_stn, SNOOBS_stn, SNOFCS_at_stn, OmB_innov_at_stn,  &                                
!                               innov_at_Grid, SNCOV_IMS, &
!                               num_Eval, Lat_atEvalPts, Lon_atEvalPts, Obs_atEvalPts, & 
!                               SNOFCS_atEvalPts, innov_atEvalPts, SNOANL_atEvalPts) !, SNOANL_Cur_atEvalPts)  !, anl_fSCA) !updated snocov

! 998 CONTINUE
!         ! clean up
!              if (allocated(SNOOBS_stn))      DEALLOCATE(SNOOBS_stn)
!              if (allocated(SNOFCS_at_stn))   DEALLOCATE(SNOFCS_at_stn)
!              if (allocated(OmB_innov_at_stn))   DEALLOCATE(OmB_innov_at_stn)
!              if (allocated(Lat_stn))         DEALLOCATE(Lat_stn) 
!              if (allocated(Lon_stn))         DEALLOCATE(Lon_stn) 
!              if (allocated(OROGFCS_at_stn))  DEALLOCATE(OROGFCS_at_stn) 
!              if (allocated(Obs_atEvalPts))   DEALLOCATE(Obs_atEvalPts)
!              if (allocated(SNOFCS_atEvalPts)) DEALLOCATE(SNOFCS_atEvalPts)
!              if (allocated(innov_atEvalPts))  DEALLOCATE(innov_atEvalPts)
!              if (allocated(SNOANL_atEvalPts)) DEALLOCATE(SNOANL_atEvalPts)
!              if (allocated(index_back_atObs)) DEALLOCATE(index_back_atObs)
!              if (allocated(Lat_atEvalPts))   DEALLOCATE(Lat_atEvalPts) 
!              if (allocated(Lon_atEvalPts))   DEALLOCATE(Lon_atEvalPts) 
! 999 CONTINUE
!         PRINT*,'Finished OI DA ON RANK: ', MYRANK
!         CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

!         !STOP

!         RETURN

!  END subroutine Snow_Analysis_OI_multiTime
    

 END MODULE M_Snow_Analysis
 
