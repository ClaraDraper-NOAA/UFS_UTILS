MODULE M_Snow_Analysis

USE NETCDF
USE M_DA
!USE MPI
Use, Intrinsic :: IEEE_ARITHMETIC	

CONTAINS

subroutine Snow_Analysis_EnSRF(SNOW_IO_TYPE, MAX_TASKS, MYRANK, NPROCS, IDIM, JDIM, IY, IM, ID, IH, &
				   DELTSFC, LENSFC, IVEGSRC, GHCND_SNOWDEPTH_PATH, IMS_SNOWCOVER_PATH, & 
				   IMS_INDEXES_PATH, CURRENT_ANALYSIS_PATH, ENKFGDAS_TOP_DIR, SNOANL)
						   
   !----------------------------------------------------------------------
   ! Input: forecast/background states for a single tile by a single MPI process
   ! reads observations: snow depth / SWE and cover
   ! calls observation operator for a given tile 
   ! does OI update per tile
   ! returns updated/Analysis states back to the caller, surface drive (sfcdrv) program
   !
   ! RLA, RLO: lat lon information for the tile
   ! SNOFCS(LENSFC), SWEFCS(LENSFC): snowdepth and snow water equivalent forecat (background)
   ! (background/forecast in arrays of LENSFC)
   ! compute snow depnsity (SNWDEN) and fraction of snow (if needed) from SWE/SNWD 
   !
   ! IDIM * JDIM = LENSFC: number of grid cells in tile = xdim * ydim   
   ! IY, IM, ID, IH = year, month, day, hour of current model step   
   ! MYRANK: rank/id of the MPI process
   !
   ! Outputs:
   ! SNOANL:  SWE and snowdepthafter DA
   !                   
   !----------------------------------------------------------------------
   IMPLICIT NONE
   !
   include 'mpif.h'
   
   integer, parameter :: dp = kind(1.d0)

   INTEGER, intent(in) :: SNOW_IO_TYPE, MAX_TASKS, MYRANK, NPROCS, IDIM, JDIM, IY, IM, ID, IH, LENSFC, IVEGSRC
   REAL, intent(in)    :: DELTSFC
   CHARACTER(LEN=*), Intent(In)   :: GHCND_SNOWDEPTH_PATH, IMS_SNOWCOVER_PATH, & 
								IMS_INDEXES_PATH, CURRENT_ANALYSIS_PATH, ENKFGDAS_TOP_DIR
!  GHCND_SNOWDEPTH_PATH="/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/GHCND/",
!  IMS_SNOWCOVER_PATH="/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS/",
!  IMS_INDEXES_PATH="/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS_INDEXES/",
!  CURRENT_ANALYSIS_PATH = "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/Cur_Analysis/"
   REAL, intent(Out)   :: SNOANL(LENSFC)   !, SWDANL(LENSFC), anl_fSCA(LENSFC)
   REAL                :: SWDANL(LENSFC)
   CHARACTER(LEN=5)    :: TILE_NUM
   Character(LEN=3)    :: rank_str
   INTEGER			    :: IERR	
   REAL       :: RLA(LENSFC), RLO(LENSFC), RLO_Tile(LENSFC), OROG(LENSFC), OROG_UF(LENSFC)
   REAL     :: SNOFCS(LENSFC), SWEFCS(LENSFC), SNWDEN(LENSFC), VETFCS(LENSFC)	
   CHARACTER(len=250)   :: dim_name, ghcnd_inp_file, ims_inp_file, ims_inp_file_indices
   CHARACTER(len=250) 	 :: anl_inp_path, da_out_file
   LOGICAL 		         :: save_Ens  ! save ens inputs for inspection
   CHARACTER(len=5)     :: y_str, m_str, d_Str, h_str, hprev_str, fvs_tile
   REAL, ALLOCATABLE    :: SNOFCS_atGHCND(:), SNOObs_atGHCND(:)		!, SWDFCS_atGHCND(:)
   REAL, ALLOCATABLE    :: Lat_GHCND(:), Lon_GHCND(:), Ele_GHCND(:)  
   REAL                 :: lat_min, lat_max, lon_min, lon_max  	
   Real		     :: SNCOV_IMS(LENSFC)  ! ims resampled at each grid
   Real		     :: SNWD_IMS_at_Grid(LENSFC)
   INTEGER :: Num_Ghcnd, Num_Ims, num_Eval, num_sub !Num_Ims_Lat, Num_Ims_Lon
   Real	:: max_distance   ! radius_of_influence for selecting state at observation point
   INTEGER :: jndx, zndx, ncol, nrow
   Integer, Allocatable   :: loc_nearest_Obs(:)
   Integer				   :: num_loc, num_loc_1, max_num_loc
   Real, Parameter		:: Stdev_back = 30., Stdev_Obs_depth = 40., Stdev_Obs_ims = 80. ! mm 
   Integer				:: ims_assm_hour
   Real				:: obs_tolerance, ims_max_ele
   Real                :: L_horz, h_ver
   Real, Allocatable   :: obs_Array(:), Lat_Obs(:), Lon_Obs(:), Ele_Obs(:)
   Real               :: SNOANL_Cur(LENSFC), SWDANL_Cur(LENSFC)
   REAL, ALLOCATABLE  :: SNOANL_Cur_atEvalPts(:)  !evalution points 
   REAL, ALLOCATABLE  :: SNOFCS_atEvalPts(:), innov_atEvalPts(:), SNOANL_atEvalPts(:)  !evalution points 
   REAL, ALLOCATABLE  :: Lat_atEvalPts(:), Lon_atEvalPts(:), Obs_atEvalPts(:)     !evalution points
   Integer, ALLOCATABLE 	:: index_back_atEval(:)     ! background locations at eval points 
   Integer, ALLOCATABLE	:: index_back_atObs(:)   ! the location of background corresponding obs
! ens
   Real, ALLOCATABLE    :: SNOFCS_Inp_Ens(:,:), SNOFCS_atObs_ens(:,:)  !, back_at_Obs_ens(:,:)
   INTEGER              :: ens_size
   REAL, Allocatable    :: anl_at_Grid_ens(:,:)  !, innov_at_Grid_ensM(:),
   REAL                 :: innov_at_Grid_ensM(LENSFC)
   Real, Allocatable    :: obs_Innov_ensM(:)
   REAL                 :: RLA_cg(LENSFC/4), RLO_cg(LENSFC/4)

   Real			   :: ims_threshold      ! threshold for converting IMS fSCA to binary 1, 0 values	
   Real			   :: snwden_val
   LOGICAL 		   :: assim_SWE, assim_GHCND, assim_IMS, assim_IMS_thisGridCell, rcov_localize
   CHARACTER(len=5)   :: fvs_tile_coarse
   CHARACTER(len=250) 	 :: ens_inp_path, gdenkf_grid_inp_path
   Integer              :: index_ens_atGrid(LENSFC)

   Integer 		   :: veg_type_landice  ! 10.21.20: no assmn over land ice
   ! for mpi par
   INTEGER            :: Np_ext, Np_til, p_tN, p_tRank, N_sA, N_sA_Ext, mp_start, mp_end
   INTEGER            :: send_proc, rec_stat(MPI_STATUS_SIZE), dest_Aoffset, pindex
   INTEGER            :: mpiReal_size, rsize
   
   assim_SWE = .False.  ! note: if this is set true, may need to adjust background and obs stddev above
   If (SNOW_IO_TYPE == 2) assim_SWE = .True.
   assim_GHCND = .True. 
   assim_IMS = .True.     !!assimilate sncov; 
   assim_IMS_thisGridCell = .FALSE.    ! if assimilating ims, skip this grid cell for this time step
   rcov_localize = .False.  ! localize R Covariance 

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
   mp_end = (p_tRank + 1) * N_sA + N_sA_Ext 		! end index of subarray for proc	
   If(myrank == 0 )PRINT*,"sub array length ", N_sA, " extra sub array: ", N_sA_Ext
! if (p_tN /= 2 ) goto 999

   if( (GHCND_SNOWDEPTH_PATH(1:8).eq.'        ').and. & 
	    (IMS_SNOWCOVER_PATH(1:8).eq.'        ') ) then
		print*, "Observation paths don't exist!, skipping the OI DA"
		goto 999
	end if

! ToDO: use file location the way the rest of the program accesses
   ! 4.9.20 for now hard code the file location 
   write(y_str, "(I4)") IY
   write(m_str, "(I0.2)") IM
   write(d_str, "(I0.2)") ID
   write(h_str, "(I0.2)") IH
   write(hprev_str, "(I0.2)") (IH - INT(DELTSFC))
   write(fvs_tile, "(I3)") IDIM
   Write(rank_str, '(I3.3)') (MYRANK+1)

   ! READ THE OROGRAPHY AND GRID POINT LAT/LONS FOR THE CUBED-SPHERE TILE.
   CALL READ_LAT_LON_OROG_atRank(p_tN, RLA,RLO,OROG,OROG_UF,TILE_NUM,IDIM,JDIM,LENSFC)
   PRINT*,"Snow anl on ", MYRANK, " Tile group: ", p_tN, " Tile: ", TILE_NUM
  ! READ THE INPUT SURFACE DATA ON THE CUBED-SPHERE TILE.
   Call READ_Forecast_Data(p_tN, LENSFC, SWEFCS, SNOFCS, VETFCS)  !VEGFCS, 
! 10.22.20 Cur_Analysis: exiting SNODEP based analysis
   !"/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/Cur_Analysis/"
   anl_inp_path = TRIM(CURRENT_ANALYSIS_PATH)//"snow."// &
				TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"." // &
				TRIM(h_str)// "0000.sfcanl_data."//TILE_NUM//".nc"
   Call READ_Analysis_Data(anl_inp_path, p_tN, LENSFC, SNOANL_Cur, SWDANL_Cur)
   
   write(fvs_tile_coarse, "(I3)") IDIM/2
   !/scratch1/NCEPDEV/global/glopara/fix/fix_fv3/C768/C768_grid.tileX.nc
   gdenkf_grid_inp_path = "/scratch1/NCEPDEV/global/glopara/fix/fix_fv3/C"// & 
	   TRIM(ADJUSTL(fvs_tile_coarse))//"/C"//TRIM(ADJUSTL(fvs_tile_coarse))// &
	   "_grid."//TILE_NUM//".nc"
   Call READ_LAT_LON_CoarseRes(gdenkf_grid_inp_path,RLA_cg,RLO_cg,IDIM/2,JDIM/2,LENSFC/4)
   ! ensemble grid indices from coarse (eg C384) to high res (C768)
   Call Find_Nearest_GridIndices_Parallel(MYRANK, MAX_TASKS, p_tN, p_tRank, Np_til, &
							   LENSFC/4, LENSFC, RLA_cg, RLO_cg, RLA, RLO, index_ens_atGrid)
   ! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/Analysis
   da_out_file = "./SNOENS."// &  !
				 TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18_tile"//rank_str//".nc"  !
   ! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/ROTDIRS/dec15/"
   ens_inp_path = TRIM(ENKFGDAS_TOP_DIR)//"enkfgdas."//TRIM(y_str)//TRIM(m_str)// &
				  TRIM(d_str)//"/"//TRIM(hprev_str)//"/"
   ens_size = 20					
   Allocate(SNOFCS_Inp_Ens(ens_size, LENSFC))
   !Allocate(innov_at_Grid_ensM(LENSFC))
   Allocate(anl_at_Grid_ens(LENSFC, ens_size+1))	
   save_Ens = .False.
   Call read_ensemble_forcing(MYRANK, MAX_TASKS, TILE_NUM, ens_size, LENSFC, IDIM, JDIM, & 
							   index_ens_atGrid,   &  !RLA_cg, RLO_cg, RLA, RLO,    &
							   ens_inp_path, y_str, m_str, d_str, h_str, hprev_str, &
							   assim_SWE, save_Ens, SNOFCS_Inp_Ens, da_out_file) 
   ! Stop   

   ! data_dir = /scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/GHCND/GHCND.SNWD.2019100118.nc
   RLO_Tile = RLO
   Where(RLO_Tile > 180) RLO_Tile = RLO_Tile - 360
   lat_min = MAX(MINVAL(RLA) - 1., -90.)
   lat_max = MIN(MAXVAL(RLA) + 1., 90.)
   lon_min = MAX(MINVAL(RLO_Tile) - 1., -180.)
   lon_max = MIN(MAXVAL(RLO_Tile) + 1., 180.)	
   ! lon_min = MAX(MINVAL(RLO) - 1., 0.)
   ! lon_max = MIN(MAXVAL(RLO) + 1., 360.)	
   if (p_tN==3)  then  
	   lon_min = 125.       ! lon_min is left, lon_max is right, not necessary min/max value
	   lon_max = -145.
   endif
   if ((p_tRank==0) ) then  !.and. print_deb
	   print*, "Tile ", p_tN, " min/max lat/lon ", lat_min, lat_max, lon_min, lon_max
   endif
   ! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/GHCND/
   ghcnd_inp_file = TRIM(GHCND_SNOWDEPTH_PATH)//"GHCND.SNWD."// &
						TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18.nc"
   dim_name = "Site_Id"	
   call Observation_Read_GHCND_Tile_excNaN(p_tN, ghcnd_inp_file, dim_name, &
				   lat_min, lat_max, lon_min, lon_max, & 
				   Num_Ghcnd, SNOObs_atGHCND,		&
				   Lat_GHCND, Lon_GHCND, MYRANK)  !Ele_GHCND,		&
   ! Do jndx=0, 5
   ! 	if ((p_tN == jndx) .and. (p_tRank==0)) then
   ! 		print*, "Tile ", p_tN, " num. GHCND obs ", Num_Ghcnd
   ! 		PRINT*, "GHCND SNWD from rank: ", MYRANK
   ! 		PRINT*, SNOObs_atGHCND
   ! 	endif
   ! End do
   ! Stop
   ! call Observation_Read_GHCND_Tile(ghcnd_inp_file, dim_name, &
   ! 				lat_min, lat_max, lon_min, lon_max, & 
   ! 				Num_Ghcnd, SNOObs_atGHCND,		&
   ! 				Lat_GHCND, Lon_GHCND, MYRANK)  !Ele_GHCND,		&
   ! call Observation_Read_GHCND(ghcnd_inp_file, dim_name, Num_Ghcnd, SNOObs_atGHCND,		&
   ! 				Lat_GHCND, Lon_GHCND, MYRANK)  !Ele_GHCND,		&
   
   if ((p_tRank==0) .and. print_deb) then
	   print*, "Tile ", p_tN, " num. GHCND obs ", Num_Ghcnd
   endif
   if ((p_tRank==0) .and. (p_tN==2) .and. print_deb) then
	   PRINT*, "GHCND SNWD from rank: ", MYRANK
	   PRINT*, SNOObs_atGHCND
	   PRINT*, "Lat at GHCND from rank: ", MYRANK
	   PRINT*, Lat_GHCND
	   PRINT*, "Lon at GHCND from rank: ", MYRANK
	   PRINT*, Lon_GHCND
	   PRINT*,'Finished reading GHCND'
   endif
   if (myrank==0) PRINT*,'Finished reading GHCND'
  
   ! (max) number of IMS subcells within a tile grid cell
   If (IDIM == 96) then          
	   num_sub = 627               
	   max_distance = 240.		!Km radius of influence: distance from gridcell to search for observations
   elseif (IDIM == 128) then
	   num_sub = 627               
	   max_distance = 240.		!Km 
   elseif (IDIM == 192) then
	   !num_sub = 30 
	   PRINT*,'Error, tile type not known '
	   stop
   elseif (IDIM == 384) then
	   !num_sub = 30  
	   PRINT*,'Error, tile type not known '
	   stop      
   elseif (IDIM == 768) then
	   num_sub = 30
	   max_distance = 27.        			!Km  
   else
	   PRINT*,'Error, tile type not known '
	   stop   
   endif
   ! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS/"
	ims_inp_file = TRIM(IMS_SNOWCOVER_PATH)//"IMS.SNCOV."// &
				   TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18.nc"                      !
	! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS_INDEXES/"
	ims_inp_file_indices = TRIM(IMS_INDEXES_PATH)//"C"//TRIM(ADJUSTL(fvs_tile))//&
						   ".IMS.Indices."//TRIM(TILE_NUM)//".nc"				
   Call Observation_Read_IMS_Full(ims_inp_file, ims_inp_file_indices, &
								  MYRANK, JDIM, IDIM, num_sub, SNCOV_IMS)
   if((p_tRank==0) .and. (p_tN==2) .and. print_deb) then
	   PRINT*, "IMS SNCOV from rank: ", MYRANK
	   PRINT*, SNCOV_IMS
	   PRINT*,'Finished reading IMS'
   endif 
   if (myrank==0) PRINT*,'Finished reading IMS'
   
   ! 4.8.20: compute snow density ratio from forecast swe and snwdepth
   SNWDEN = SNOFCS / SWEFCS
   snwden_val = SUM(SNWDEN, Mask = (.not. IEEE_IS_NAN(SNWDEN))) &
					/ COUNT (.not. IEEE_IS_NAN(SNWDEN))
   SNWDEN = snwden_val  !10.   ! ratio of sndepth to swe when swe /=0 : SNOFCS / SWEFCS
   if(p_tRank == 0) then
	   print*, "Process ", MYRANK, " average snwd/swe ratio: ", snwden_val
   endif
   ! Get model states at obs points
   ALLOCATE(SNOFCS_atGHCND(Num_Ghcnd))
   ALLOCATE(Ele_GHCND(Num_Ghcnd)) 
   ALLOCATE(index_back_atObs(Num_Ghcnd)) 
   !ALLOCATE(index_back_atGHCND(Num_Ghcnd)) 
   num_Eval = floor(0.05 * Num_Ghcnd)      ! using 5% of ghcnd locations for evaluation
   ALLOCATE(index_back_atEval(num_Eval)) 
   ALLOCATE(Obs_atEvalPts(num_Eval)) 
   ALLOCATE(SNOFCS_atEvalPts(num_Eval)) 
   ALLOCATE(Lat_atEvalPts(num_Eval))
   ALLOCATE(Lon_atEvalPts(num_Eval)) 
   ALLOCATE(innov_atEvalPts(num_Eval))
   ALLOCATE(SNOANL_atEvalPts(num_Eval)) 
   ALLOCATE(SNOANL_Cur_atEvalPts(num_Eval)) 	 	
   if(p_tRank == 0) then 
	   PRINT*, "Tile ", p_tN+1, " ", num_Eval, ' evaluation points to be excluded from DA'	
   endif	
   Call Observation_Operator_Parallel(Myrank, MAX_TASKS, p_tN, p_tRank, Np_til, & 
					   RLA, RLO, OROG, Lat_GHCND, Lon_GHCND,   &
					   LENSFC, Num_Ghcnd, num_Eval, max_distance, SNOFCS, SNOObs_atGHCND,  &
					   SNOFCS_atGHCND, Ele_GHCND, index_back_atObs, index_back_atEval, Obs_atEvalPts,    &
					   SNOFCS_atEvalPts, Lat_atEvalPts, Lon_atEvalPts)
   ! Call Observation_Operator_Snofcs_Parallel(Myrank, MAX_TASKS, p_tN, p_tRank, Np_til, & 
   !                     RLA, RLO, OROG, Lat_GHCND, Lon_GHCND,   &
   ! 					LENSFC, Num_Ghcnd, num_Eval, max_distance, SNOFCS, SNOObs_atGHCND,  &
   ! 					Ele_GHCND, index_back_atObs, index_back_atEval, Obs_atEvalPts,    &
   ! 					SNOFCS_atEvalPts, Lat_atEvalPts, Lon_atEvalPts)
   ALLOCATE(SNOFCS_atObs_ens(ens_size, Num_Ghcnd))
   SNOFCS_atObs_ens = IEEE_VALUE(SNOFCS_atObs_ens, IEEE_QUIET_NAN)
   Do jndx = 1, Num_Ghcnd
	   if (index_back_atObs(jndx) > 0) then
		   SNOFCS_atObs_ens(:, jndx) = SNOFCS_Inp_Ens(:, index_back_atObs(jndx))
	   endif
   end do

   if ((p_tN==2) .and. (p_tRank==0) .and. print_deb) then 
	   PRINT*, "Background Indices at eval points"
	   PRINT*, index_back_atEval
	   PRINT*, "Background Indices at obs points"
	   PRINT*, index_back_atObs	
	   PRINT*, "Obs at Eval Points" 
	   PRINT*, Obs_atEvalPts	
	   PRINT*, "Forecast at Eval Points"
	   PRINT*, SNOFCS_atEvalPts	
	   PRINT*, "Lat at Eval Points"
	   PRINT*, Lat_atEvalPts
	   PRINT*, "Lon at Eval Points"
	   PRINT*, Lon_atEvalPts
   endif
   !Stop
   ! call Observation_Operator(RLA, RLO, OROG, Lat_GHCND, Lon_GHCND,   &
   ! 					LENSFC, Num_Ghcnd, max_distance,            &
   ! 					SNOFCS,                                  &
   ! 					SNOFCS_atGHCND, Ele_GHCND, index_back_atGHCND)
   if ((p_tN==2) .and. (p_tRank==0) .and. print_deb) then
	   PRINT*, "GHCND Lat range from rank: ", MYRANK, MINVAL(Lat_GHCND), " ", MAXVAL(Lat_GHCND)
	   PRINT*, "GHCND Lon range from rank: ", MYRANK, MINVAL(Lon_GHCND), " ", MAXVAL(Lon_GHCND)
	   PRINT*, "Background at GHCND locations from rank: ", MYRANK
	   PRINT*, SNOFCS_atGHCND   !SNOFCS_atObs_ens   !	
	   PRINT*, "Elevation at GHCND locations from rank: ", MYRANK
	   PRINT*, Ele_GHCND
	   PRINT*,'Finished observation operator for GHCND'	
   endif
   if (myrank==0) PRINT*,'Finished observation operator for GHCND'		
   ims_threshold = 0.5  ! threshold for converting IMS fSCA to binary 1, 0 values
   ! Call Observation_Operator_IMS_fSCA_Threshold(SNCOV_IMS, SNOFCS, SNWDEN, assim_SWE,	   &
   ! 						LENSFC, ims_threshold, SNWD_IMS_at_Grid) 
   call Observation_Operator_IMS_fSCA(SNCOV_IMS, SNWDEN, VETFCS, assim_SWE, LENSFC, SNWD_IMS_at_Grid) 
   if ((p_tN==2) .and. (p_tRank==0) .and. print_deb) then
	   PRINT*, "IMS obs at each grid cell from rank: ", MYRANK
	   PRINT*, SNWD_IMS_at_Grid
	   PRINT*,'Finished observation operator for IMS'
   endif
   if (myrank==0) PRINT*,'Finished observation operator for IMS'
   
   L_horz = 55.  !120.  !
   h_ver = 800.  !1200.	!
   obs_tolerance = 5.0
   max_num_loc = 50   !100	!
   ims_max_ele = 1500.
   ims_assm_hour = 18
   max_distance = 250.   !Km 120.			!
   !if (myrank==4) then 
   if (myrank==0) PRINT*,'Starting DA loop'

   ! 10.21.20: no assimilation if land-ice
	if (IVEGSRC == 2) then   ! sib
		veg_type_landice=13
	else
		veg_type_landice=15
	endif
   !Do ncol=1, IDIM
   !Do nrow = 1, Num_Snotel !JDIM/2
	   Do jndx = mp_start, mp_end     !1, LENSFC		!jndx = (ncol-1)*JDIM + nrow		! 
		   num_loc_1 = 0
		   assim_IMS_thisGridCell = .FALSE.
		   call debug_print("loop ", float(jndx))
		   ! GHCND			 
		   if(assim_GHCND) then
			   ! 8.19.20: we are assuming here if deterministic forecast exists/not null
			   ! at a point, then ensemble members also have valid value
			   call nearest_Observations_Locations(RLA(jndx), RLO(jndx),    &
					   Lat_GHCND, Lon_GHCND,  Num_Ghcnd, max_distance, max_num_loc,   &
					   Stdev_back, Stdev_Obs_depth, obs_tolerance,                 &
					   SNOFCS_atGHCND, SNOObs_atGHCND,						 &
					   loc_nearest_Obs,  num_loc_1) !,      &LENSFC,
			   call debug_print("number of GHCND sndpth obs ", float(num_loc))
		   endif
		   num_loc = num_loc_1
		   !check IMS is assimilated
		   if (assim_IMS) then
			   if((.NOT. IEEE_IS_NAN(SNWD_IMS_at_Grid(jndx))) .AND. &
				  (OROG(jndx) <= ims_max_ele) .AND. &
				  (IH == ims_assm_hour)) then
				   num_loc = num_loc + 1
				   assim_IMS_thisGridCell = .TRUE.
			   endif
		   endif
		   ! if assim_IMS=false >> num_loc_1=num_loc
		   ! 10.21.20: no assimilation if land-ice	.NE. veg_type_landice	
			if((num_loc > 0) .and. (NINT(VETFCS(jndx)) /= veg_type_landice)) then ! .and. (SNCOV_IMS(jndx) > 0.)) then    
			   Allocate(obs_Array(num_loc))
			   Allocate(Lat_Obs(num_loc))
			   Allocate(Lon_Obs(num_loc))
			   Allocate(Ele_Obs(num_loc))
			   ! ghcnd
			   if(num_loc_1 > 0) then
				   Do zndx = 1, num_loc_1     
					   obs_Array(zndx) = SNOObs_atGHCND(loc_nearest_Obs(zndx))
					   Lat_Obs(zndx) = Lat_GHCND(loc_nearest_Obs(zndx))
					   Lon_Obs(zndx) = Lon_GHCND(loc_nearest_Obs(zndx))
					   Ele_Obs(zndx) = Ele_GHCND(loc_nearest_Obs(zndx))
				   End Do
			   End if
			   !ims
			   if(assim_IMS_thisGridCell) then
				   obs_Array(num_loc) = SNWD_IMS_at_Grid(jndx)
				   Lat_Obs(num_loc) = RLA(jndx)   !Lat_IMS_atGrid(jndx)
				   Lon_Obs(num_loc) = RLO(jndx)   !Lon_IMS_atGrid(jndx)
				   Ele_Obs(num_loc) = OROG(jndx)  !Ele_IMS(jndx)
			   endif	
			   ! EnKF
			   Allocate(obs_Innov_ensM(num_loc))	
			   Call snow_DA_EnSRF(RLA(jndx), RLO(jndx), OROG(jndx),     &
			   Lat_Obs, Lon_Obs, Ele_Obs,     				&
			   L_horz, h_ver,                                      &   !L_horz in Km, h_ver in m
			   assim_IMS_thisGridCell, rcov_localize,                   &
			   jndx, ens_size, LENSFC, SNOFCS_Inp_Ens,          &
			   num_loc_1, num_loc, loc_nearest_Obs, SNOFCS_atObs_ens,	 &
			   Stdev_Obs_depth, Stdev_Obs_ims,                  &
			   obs_Array,                          &
			   obs_Innov_ensM, innov_at_Grid_ensM(jndx), anl_at_Grid_ens(jndx,:))
			   if ((p_tN==4) .and. (p_tRank==0) .and. print_deb) then	
				   print*, "proc ", myrank, "loop ", jndx, &
				   "num depth obs ", num_loc_1, "total obs", num_loc, " ens size ", ens_size
				   PRINT*, "Ens background at obs pts: "
				   PRINT*, SNOFCS_Inp_Ens(:, jndx)	
				   PRINT*, "Observed"
				   PRINT*,  obs_Array
				   PRINT*, "Ens mean Obs innovation: "
				   PRINT*, obs_Innov_ensM
				   PRINT*, "Ens mean innovation at grid: "
				   PRINT*, innov_at_Grid_ensM(jndx)
				   PRINT*, "Ens analysis at grid: "
				   PRINT*, anl_at_Grid_ens(jndx, :)
			   endif					
			   DEALLOCATE(obs_Array, obs_Innov_ensM)
			   DEALLOCATE(Lat_Obs, Lon_Obs, Ele_Obs)
			   ! QCC by ims--use ims as snow mask
			   Do zndx = 1, ens_size+1  
				   if((SNCOV_IMS(jndx) >= 0.5) .and. (anl_at_Grid_ens(jndx, zndx) < 50.)) then
					   anl_at_Grid_ens(jndx, zndx) = 50.
				   ! elseif((SNCOV_IMS(jndx) < 0.5) .and. (SNCOV_IMS(jndx) > 0.1) .and. (anl_at_Grid(jndx) >= 50)) then
				   ! 	anl_at_Grid(jndx) = 50.
				   elseif((SNCOV_IMS(jndx) < 0.1) .and. (anl_at_Grid_ens(jndx, zndx) >= 0.)) then
					   anl_at_Grid_ens(jndx, zndx) = 0.
				   endif
			   End do
			   if ((p_tN==4) .and. (p_tRank==0) .and. print_deb) then						
				   PRINT*, "Ens analysis at grid after mask: "
				   PRINT*, anl_at_Grid_ens(jndx, :)
			   endif
		   else
			   anl_at_Grid_ens(jndx, :) = SNOFCS_Inp_Ens(jndx, :)
		   endif
		   if (assim_GHCND) Deallocate(loc_nearest_Obs) 
	   End do
   !End do
   if (myrank==0) PRINT*, 'Finished DA loops'

! ToDO: Better way to handle this?
! Real data type size corresponding to mpi
   rsize = SIZEOF(snwden_val) 
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
   ! ens analysis innov_at_Grid_ens(jndx), anl_at_Grid_ens(jndx)
   if (MYRANK > (MAX_TASKS - 1) ) then
	   call MPI_SEND(anl_at_Grid_ens(mp_start:mp_end,:),N_sA*(ens_size+1), mpiReal_size,    &
					 p_tN, MYRANK*1000, MPI_COMM_WORLD, IERR) 
	   call MPI_SEND(innov_at_Grid_ensM(mp_start:mp_end), N_sA, mpiReal_size,    &
					 p_tN, MYRANK*10000, MPI_COMM_WORLD, IERR)
   else    !if(p_tRank == 0) then  
	   Do pindex =  1, (Np_til - 1)   ! sender proc index within tile group
		   dest_Aoffset = pindex * N_sA + N_sA_Ext + 1   ! dest array offset
		   send_proc = MYRANK +  pindex * MAX_TASKS
		   call MPI_RECV(anl_at_Grid_ens(dest_Aoffset:dest_Aoffset+N_sA-1, :), N_sA*(ens_size+1),      &
		   mpiReal_size, send_proc, send_proc*1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
		   call MPI_RECV(innov_at_Grid_ensM(dest_Aoffset:dest_Aoffset+N_sA-1), N_sA,       &
		   mpiReal_size, send_proc, send_proc*10000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
	   enddo
   endif
   if (myrank==0) PRINT*,'Finished Data copy'

   if (MYRANK > MAX_TASKS - 1 ) goto 998   ! if(p_tRank /= 0 ) goto 998

   ! avoid -ve anl
   Where(anl_at_Grid_ens < 0.) anl_at_Grid_ens = 0.
   ! swe and snwd
   SNOANL = anl_at_Grid_ens(:,ens_size+1) / SNWDEN		!snwden_val 
   SWDANL = anl_at_Grid_ens(:,ens_size+1)	
!    if(assim_SWE) then 
!       SNOANL_Cur = SNOANL_Cur 
!    else
! 		SNOANL_Cur = SWDANL_Cur
!    end if
   if (print_deb) then
	   PRINT*, "Innovation SWE/snwd from rank: ", MYRANK
	   PRINT*, innov_at_Grid_ensM	
	   PRINT*, "Analysis SWE/ snwd  from rank: ", MYRANK
	   PRINT*, anl_at_Grid_ens	
   endif

   !Compute updated snocov	
   !Call update_snow_cover_fraction(LENSFC, SNOANL, VETFCS, anl_fSCA)

   ! copy values at eval points
   innov_atEvalPts = IEEE_VALUE(innov_atEvalPts, IEEE_QUIET_NAN)
   SNOANL_atEvalPts = IEEE_VALUE(SNOANL_atEvalPts, IEEE_QUIET_NAN)
   SNOANL_Cur_atEvalPts = IEEE_VALUE(SNOANL_Cur_atEvalPts, IEEE_QUIET_NAN)
   Do jndx = 1, num_Eval
	   if ((index_back_atEval(jndx) > 0) .and. &
	       (NINT(VETFCS(jndx)) /= veg_type_landice)) then  !10.21.20: exclude if land-ice
		    	innov_atEvalPts(jndx) = innov_at_Grid_ensM(index_back_atEval(jndx))
		   		SNOANL_atEvalPts(jndx) = anl_at_Grid_ens(index_back_atEval(jndx), ens_size+1)
		   		SNOANL_Cur_atEvalPts(jndx) = SWDANL_Cur(index_back_atEval(jndx)) !SNOANL_Cur
	   endif
   End do

   ! write outputs	
   Write(rank_str, '(I3.3)') (MYRANK+1)
   ! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/Analysis"
   da_out_file = "./SNOANLEnSRF."// &  !
				 TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18_tile"//rank_str//".nc"  !   
   call Write_DA_Outputs(da_out_file, IDIM, JDIM, LENSFC, MYRANK, &
						 SWEFCS, SNOANL, SNOFCS, SWDANL, &  !
						 innov_at_Grid_ensM, SNCOV_IMS, &
				 !   Num_Ghcnd, Lat_GHCND, Lon_GHCND, SNOObs_atGHCND, &
			   ! SNOFCS_atGHCND, SNOFCS_atGHCND, SNOFCS_atGHCND)
			num_Eval, Lat_atEvalPts, Lon_atEvalPts, Obs_atEvalPts, & 
		SNOFCS_atEvalPts, innov_atEvalPts, SNOANL_atEvalPts, SNOANL_Cur_atEvalPts)  !, anl_fSCA) !updated snocov

998 CONTINUE
   DEALLOCATE(SNOObs_atGHCND, SNOFCS_atGHCND, Lat_GHCND, Lon_GHCND, Ele_GHCND)  !, index_back_atGHCND)
   DEALLOCATE(Obs_atEvalPts, SNOFCS_atEvalPts, innov_atEvalPts, SNOANL_atEvalPts)
   DEALLOCATE(SNOANL_Cur_atEvalPts) 
   DEALLOCATE(index_back_atObs, index_back_atEval, Lat_atEvalPts, Lon_atEvalPts) !, Ele_atEvalPts)
   DEAllocate(SNOFCS_Inp_Ens, SNOFCS_atObs_ens)
   DEAllocate(anl_at_Grid_ens)  !innov_at_Grid_ens, 
   !DEALLOCATE(SNCOV_IMS, Lat_IMS, Lon_IMS, SNOFCS_atIMS, SWDFCS_atIMS)
999 CONTINUE
   PRINT*,'Finished EnSRF DA ON RANK: ', MYRANK
   CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

   !STOP

   RETURN

END subroutine Snow_Analysis_EnSRF

subroutine Snow_Analysis_EnKF(SNOW_IO_TYPE, MAX_TASKS, MYRANK, NPROCS, IDIM, JDIM, IY, IM, ID, IH, &
	               DELTSFC, LENSFC, IVEGSRC, GHCND_SNOWDEPTH_PATH, IMS_SNOWCOVER_PATH, & 
	            IMS_INDEXES_PATH, CURRENT_ANALYSIS_PATH, ENKFGDAS_TOP_DIR, SNOANL)
							
	!----------------------------------------------------------------------
	! Input: forecast/background states for a single tile by a single MPI process
	! reads observations: snow depth / SWE and cover
	! calls observation operator for a given tile 
	! does OI update per tile
	! returns updated/Analysis states back to the caller, surface drive (sfcdrv) program
	!
	! RLA, RLO: lat lon information for the tile
	! SNOFCS(LENSFC), SWEFCS(LENSFC): snowdepth and snow water equivalent forecat (background)
	! (background/forecast in arrays of LENSFC)
	! compute snow depnsity (SNWDEN) and fraction of snow (if needed) from SWE/SNWD 
	!
	! IDIM * JDIM = LENSFC: number of grid cells in tile = xdim * ydim   
	! IY, IM, ID, IH = year, month, day, hour of current model step   
	! MYRANK: rank/id of the MPI process
	!
	! Outputs:
	! SNOANL:  SWE and snowdepthafter DA
	!                   
	!----------------------------------------------------------------------
	IMPLICIT NONE
	!
	include 'mpif.h'
	
	integer, parameter :: dp = kind(1.d0)

	INTEGER, intent(in) :: SNOW_IO_TYPE, MAX_TASKS, MYRANK, NPROCS, IDIM, JDIM, IY, IM, ID, IH, LENSFC, IVEGSRC
	REAL, intent(in)    :: DELTSFC
	CHARACTER(LEN=*), Intent(In)   :: GHCND_SNOWDEPTH_PATH, IMS_SNOWCOVER_PATH, & 
							IMS_INDEXES_PATH, CURRENT_ANALYSIS_PATH, ENKFGDAS_TOP_DIR
!  GHCND_SNOWDEPTH_PATH="/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/GHCND/",
!  IMS_SNOWCOVER_PATH="/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS/",
!  IMS_INDEXES_PATH="/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS_INDEXES/",
!  CURRENT_ANALYSIS_PATH = "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/Cur_Analysis/"
	REAL, intent(Out)   :: SNOANL(LENSFC)   !, SWDANL(LENSFC), anl_fSCA(LENSFC)
	REAL                :: SWDANL(LENSFC)
	CHARACTER(LEN=5)    :: TILE_NUM
	Character(LEN=3)    :: rank_str
	INTEGER			    :: IERR	
	REAL       :: RLA(LENSFC), RLO(LENSFC), RLO_Tile(LENSFC), OROG(LENSFC), OROG_UF(LENSFC)
	REAL     :: SNOFCS(LENSFC), SWEFCS(LENSFC), SNWDEN(LENSFC), VETFCS(LENSFC)	
	CHARACTER(len=250)   :: dim_name, ghcnd_inp_file, ims_inp_file, ims_inp_file_indices
	CHARACTER(len=250) 	 :: anl_inp_path, da_out_file
	LOGICAL 		     :: save_Ens  ! save ens inputs for inspection
	CHARACTER(len=5)     :: y_str, m_str, d_Str, h_str, hprev_str, fvs_tile
	REAL, ALLOCATABLE    :: SNOFCS_atGHCND(:), SNOObs_atGHCND(:)		!, SWDFCS_atGHCND(:)
	REAL, ALLOCATABLE    :: Lat_GHCND(:), Lon_GHCND(:), Ele_GHCND(:)  
	REAL                 :: lat_min, lat_max, lon_min, lon_max  	
	Real		     :: SNCOV_IMS(LENSFC)  ! ims resampled at each grid
	Real		     :: SNWD_IMS_at_Grid(LENSFC)
	INTEGER :: Num_Ghcnd, Num_Ims, num_Eval, num_sub !Num_Ims_Lat, Num_Ims_Lon
	Real	:: max_distance   ! radius_of_influence for selecting state at observation point
	INTEGER :: jndx, zndx, ncol, nrow
	Integer, Allocatable   :: loc_nearest_Obs(:)
	Integer				   :: num_loc, num_loc_1, max_num_loc
	Real, Parameter		:: Stdev_back = 30., Stdev_Obs_depth = 40., Stdev_Obs_ims = 80. ! mm 
	Integer				:: ims_assm_hour
	Real				:: obs_tolerance, ims_max_ele
	Real                :: L_horz, h_ver
	Real, Allocatable   :: obs_Array(:), Lat_Obs(:), Lon_Obs(:), Ele_Obs(:)
	Real               :: SNOANL_Cur(LENSFC), SWDANL_Cur(LENSFC)
	REAL, ALLOCATABLE  :: SNOANL_Cur_atEvalPts(:)  !evalution points 
	REAL, ALLOCATABLE  :: SNOFCS_atEvalPts(:), innov_atEvalPts(:), SNOANL_atEvalPts(:)  !evalution points 
	REAL, ALLOCATABLE  :: Lat_atEvalPts(:), Lon_atEvalPts(:), Obs_atEvalPts(:)     !evalution points
	Integer, ALLOCATABLE 	:: index_back_atEval(:)     ! background locations at eval points 
	Integer, ALLOCATABLE	:: index_back_atObs(:)   ! the location of background corresponding obs
! ens
	Real, ALLOCATABLE    :: SNOFCS_Inp_Ens(:,:), SNOFCS_atObs_ens(:,:)  !, back_at_Obs_ens(:,:)
	INTEGER              :: ens_size
	REAL, Allocatable    :: innov_at_Grid_ens(:,:), anl_at_Grid_ens(:,:)
	Real, Allocatable    :: obs_Innov_ens(:,:)
	REAL                 :: RLA_cg(LENSFC/4), RLO_cg(LENSFC/4)

	Real			   :: ims_threshold      ! threshold for converting IMS fSCA to binary 1, 0 values	
	Real			   :: snwden_val
	LOGICAL 		   :: assim_SWE, assim_GHCND, assim_IMS, assim_IMS_thisGridCell    !assimilate sncov
	CHARACTER(len=5)   :: fvs_tile_coarse
	CHARACTER(len=250) 	 :: ens_inp_path, gdenkf_grid_inp_path
	Integer              :: index_ens_atGrid(LENSFC)

	Integer 		   :: veg_type_landice ! 10.21.20: no assmn over land ice
	! for mpi par
	INTEGER            :: Np_ext, Np_til, p_tN, p_tRank, N_sA, N_sA_Ext, mp_start, mp_end
	INTEGER            :: send_proc, rec_stat(MPI_STATUS_SIZE), dest_Aoffset, pindex
	INTEGER            :: mpiReal_size, rsize
	
	assim_SWE = .False.  ! note: if this is set true, may need to adjust background and obs stddev above
	If (SNOW_IO_TYPE == 2) assim_SWE = .True.
	assim_GHCND = .True.
	assim_IMS = .True.     !!assimilate sncov; 
	assim_IMS_thisGridCell = .FALSE.    ! if assimilating ims, skip this grid cell for this time step

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
	mp_end = (p_tRank + 1) * N_sA + N_sA_Ext 		! end index of subarray for proc	
	If(myrank == 0 )PRINT*,"sub array length ", N_sA, " extra sub array: ", N_sA_Ext
! if (p_tN /= 2 ) goto 999

	if( (GHCND_SNOWDEPTH_PATH(1:8).eq.'        ').and. & 
	    (IMS_SNOWCOVER_PATH(1:8).eq.'        ') ) then
		print*, "Observation paths don't exist!, skipping the OI DA"
		goto 999
	end if

! ToDO: use file location the way the rest of the program accesses
	! 4.9.20 for now hard code the file location 
	write(y_str, "(I4)") IY
	write(m_str, "(I0.2)") IM
	write(d_str, "(I0.2)") ID
	write(h_str, "(I0.2)") IH
	write(hprev_str, "(I0.2)") (IH - INT(DELTSFC))
	write(fvs_tile, "(I3)") IDIM
	Write(rank_str, '(I3.3)') (MYRANK+1)

	! READ THE OROGRAPHY AND GRID POINT LAT/LONS FOR THE CUBED-SPHERE TILE.
	CALL READ_LAT_LON_OROG_atRank(p_tN, RLA,RLO,OROG,OROG_UF,TILE_NUM,IDIM,JDIM,LENSFC)
	PRINT*,"Snow anl on ", MYRANK, " Tile group: ", p_tN, " Tile: ", TILE_NUM

   ! READ THE INPUT SURFACE DATA ON THE CUBED-SPHERE TILE.
	Call READ_Forecast_Data(p_tN, LENSFC, SWEFCS, SNOFCS, VETFCS)  !VEGFCS, 
! 10.22.20 Cur_Analysis: exiting SNODEP based analysis
	!"/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/Cur_Analysis/"
	anl_inp_path = TRIM(CURRENT_ANALYSIS_PATH)//"snow."// &
				TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"." // &
				TRIM(h_str)// "0000.sfcanl_data."//TILE_NUM//".nc"
	Call READ_Analysis_Data(anl_inp_path, p_tN, LENSFC, SNOANL_Cur, SWDANL_Cur)

	write(fvs_tile_coarse, "(I3)") IDIM/2
	!/scratch1/NCEPDEV/global/glopara/fix/fix_fv3/C768/C768_grid.tileX.nc
	gdenkf_grid_inp_path = "/scratch1/NCEPDEV/global/glopara/fix/fix_fv3/C"// & 
		TRIM(ADJUSTL(fvs_tile_coarse))//"/C"//TRIM(ADJUSTL(fvs_tile_coarse))// &
		"_grid."//TILE_NUM//".nc"
	Call READ_LAT_LON_CoarseRes(gdenkf_grid_inp_path,RLA_cg,RLO_cg,IDIM/2,JDIM/2,LENSFC/4)
	! ensemble grid indices from coarse (eg C384) to high res (C768)
	Call Find_Nearest_GridIndices_Parallel(MYRANK, MAX_TASKS, p_tN, p_tRank, Np_til, &
		                        LENSFC/4, LENSFC, RLA_cg, RLO_cg, RLA, RLO, index_ens_atGrid)
	! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/Analysis"
	da_out_file = "./SNOENS."// &  !
				  TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18_tile"//rank_str//".nc"  !
	! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/ROTDIRS/dec15/"
   ens_inp_path = TRIM(ENKFGDAS_TOP_DIR)//"enkfgdas."//TRIM(y_str)//TRIM(m_str)// &
				  TRIM(d_str)//"/"//TRIM(hprev_str)//"/"
	ens_size = 20					
	Allocate(SNOFCS_Inp_Ens(ens_size, LENSFC))
	Allocate(innov_at_Grid_ens(LENSFC, ens_size+1))
	Allocate(anl_at_Grid_ens(LENSFC, ens_size+1))	
	save_Ens = .False.
	Call read_ensemble_forcing(MYRANK, MAX_TASKS, TILE_NUM, ens_size, LENSFC, IDIM, JDIM, & 
								index_ens_atGrid,   &  !RLA_cg, RLO_cg, RLA, RLO,    &
	                            ens_inp_path, y_str, m_str, d_str, h_str, hprev_str, &
								assim_SWE, save_Ens, SNOFCS_Inp_Ens, da_out_file) 
	! Stop   

	! data_dir = /scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/GHCND/GHCND.SNWD.2019100118.nc
	RLO_Tile = RLO
	Where(RLO_Tile > 180) RLO_Tile = RLO_Tile - 360
	lat_min = MAX(MINVAL(RLA) - 1., -90.)
	lat_max = MIN(MAXVAL(RLA) + 1., 90.)
	lon_min = MAX(MINVAL(RLO_Tile) - 1., -180.)
	lon_max = MIN(MAXVAL(RLO_Tile) + 1., 180.)	
	! lon_min = MAX(MINVAL(RLO) - 1., 0.)
	! lon_max = MIN(MAXVAL(RLO) + 1., 360.)	
	if (p_tN==3)  then  
		lon_min = 125.       ! lon_min is left, lon_max is right, not necessary min/max value
		lon_max = -145.
	endif
	if ((p_tRank==0) ) then  !.and. print_deb
		print*, "Tile ", p_tN, " min/max lat/lon ", lat_min, lat_max, lon_min, lon_max
	endif
	! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/GHCND/
	ghcnd_inp_file = TRIM(GHCND_SNOWDEPTH_PATH)//"GHCND.SNWD."// &
					 TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18.nc"
	dim_name = "Site_Id"	
	call Observation_Read_GHCND_Tile_excNaN(p_tN, ghcnd_inp_file, dim_name, &
					lat_min, lat_max, lon_min, lon_max, & 
					Num_Ghcnd, SNOObs_atGHCND,		&
					Lat_GHCND, Lon_GHCND, MYRANK)  !Ele_GHCND,		&
	! Do jndx=0, 5
	! 	if ((p_tN == jndx) .and. (p_tRank==0)) then
	! 		print*, "Tile ", p_tN, " num. GHCND obs ", Num_Ghcnd
	! 		PRINT*, "GHCND SNWD from rank: ", MYRANK
	! 		PRINT*, SNOObs_atGHCND
	! 	endif
	! End do
	! Stop
	! call Observation_Read_GHCND_Tile(ghcnd_inp_file, dim_name, &
	! 				lat_min, lat_max, lon_min, lon_max, & 
	! 				Num_Ghcnd, SNOObs_atGHCND,		&
	! 				Lat_GHCND, Lon_GHCND, MYRANK)  !Ele_GHCND,		&
	! call Observation_Read_GHCND(ghcnd_inp_file, dim_name, Num_Ghcnd, SNOObs_atGHCND,		&
	! 				Lat_GHCND, Lon_GHCND, MYRANK)  !Ele_GHCND,		&
	
	if ((p_tRank==0) .and. print_deb) then
		print*, "Tile ", p_tN, " num. GHCND obs ", Num_Ghcnd
	endif
	if ((p_tRank==0) .and. (p_tN==2) .and. print_deb) then
		PRINT*, "GHCND SNWD from rank: ", MYRANK
		PRINT*, SNOObs_atGHCND
		PRINT*, "Lat at GHCND from rank: ", MYRANK
		PRINT*, Lat_GHCND
		PRINT*, "Lon at GHCND from rank: ", MYRANK
		PRINT*, Lon_GHCND
		PRINT*,'Finished reading GHCND'
	endif
	if (myrank==0) PRINT*,'Finished reading GHCND'
   
	! (max) number of IMS subcells within a tile grid cell
	If (IDIM == 96) then          
		num_sub = 627               
		max_distance = 240.		!Km radius of influence: distance from gridcell to search for observations
	elseif (IDIM == 128) then
	   num_sub = 627               
	   max_distance = 240.		!Km 
	elseif (IDIM == 192) then
		!num_sub = 30 
		PRINT*,'Error, tile type not known '
		stop
	elseif (IDIM == 384) then
		!num_sub = 30  
		PRINT*,'Error, tile type not known '
		stop      
	elseif (IDIM == 768) then
		num_sub = 30
		max_distance = 27.        			!Km  
	else
		PRINT*,'Error, tile type not known '
		stop   
	endif
	! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS/"
	ims_inp_file = TRIM(IMS_SNOWCOVER_PATH)//"IMS.SNCOV."// &
				   TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18.nc"                      !
	! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS_INDEXES/"
	ims_inp_file_indices = TRIM(IMS_INDEXES_PATH)//"C"//TRIM(ADJUSTL(fvs_tile))//&
						   ".IMS.Indices."//TRIM(TILE_NUM)//".nc"				
	Call Observation_Read_IMS_Full(ims_inp_file, ims_inp_file_indices, &
				                   MYRANK, JDIM, IDIM, num_sub, SNCOV_IMS)
	if((p_tRank==0) .and. (p_tN==2) .and. print_deb) then
		PRINT*, "IMS SNCOV from rank: ", MYRANK
		PRINT*, SNCOV_IMS
		PRINT*,'Finished reading IMS'
	endif 
	if (myrank==0) PRINT*,'Finished reading IMS'
	
	! 4.8.20: compute snow density ratio from forecast swe and snwdepth
	SNWDEN = SNOFCS / SWEFCS
	snwden_val = SUM(SNWDEN, Mask = (.not. IEEE_IS_NAN(SNWDEN))) &
	                 / COUNT (.not. IEEE_IS_NAN(SNWDEN))
	SNWDEN = snwden_val  !10.   ! ratio of sndepth to swe when swe /=0 : SNOFCS / SWEFCS
	if(p_tRank == 0) then
		print*, "Process ", MYRANK, " average snwd/swe ratio: ", snwden_val
	endif
	! Get model states at obs points
	ALLOCATE(SNOFCS_atGHCND(Num_Ghcnd))
	ALLOCATE(Ele_GHCND(Num_Ghcnd)) 
	ALLOCATE(index_back_atObs(Num_Ghcnd)) 
	!ALLOCATE(index_back_atGHCND(Num_Ghcnd)) 
	num_Eval = floor(0.05 * Num_Ghcnd)      ! using 5% of ghcnd locations for evaluation
	ALLOCATE(index_back_atEval(num_Eval)) 
	ALLOCATE(Obs_atEvalPts(num_Eval)) 
	ALLOCATE(SNOFCS_atEvalPts(num_Eval)) 
	ALLOCATE(Lat_atEvalPts(num_Eval))
	ALLOCATE(Lon_atEvalPts(num_Eval)) 
	ALLOCATE(innov_atEvalPts(num_Eval))
	ALLOCATE(SNOANL_atEvalPts(num_Eval))
	ALLOCATE(SNOANL_Cur_atEvalPts(num_Eval))  	 	
	if(p_tRank == 0) then 
		PRINT*, "Tile ", p_tN+1, " ", num_Eval, ' evaluation points to be excluded from DA'	
	endif	
	Call Observation_Operator_Parallel(Myrank, MAX_TASKS, p_tN, p_tRank, Np_til, & 
	                    RLA, RLO, OROG, Lat_GHCND, Lon_GHCND,   &
						LENSFC, Num_Ghcnd, num_Eval, max_distance, SNOFCS, SNOObs_atGHCND,  &
						SNOFCS_atGHCND, Ele_GHCND, index_back_atObs, index_back_atEval, Obs_atEvalPts,    &
						SNOFCS_atEvalPts, Lat_atEvalPts, Lon_atEvalPts)
	! Call Observation_Operator_Snofcs_Parallel(Myrank, MAX_TASKS, p_tN, p_tRank, Np_til, & 
	!                     RLA, RLO, OROG, Lat_GHCND, Lon_GHCND,   &
	! 					LENSFC, Num_Ghcnd, num_Eval, max_distance, SNOFCS, SNOObs_atGHCND,  &
	! 					Ele_GHCND, index_back_atObs, index_back_atEval, Obs_atEvalPts,    &
	! 					SNOFCS_atEvalPts, Lat_atEvalPts, Lon_atEvalPts)
    ALLOCATE(SNOFCS_atObs_ens(ens_size, Num_Ghcnd))
	SNOFCS_atObs_ens = IEEE_VALUE(SNOFCS_atObs_ens, IEEE_QUIET_NAN)
	Do jndx = 1, Num_Ghcnd
		if (index_back_atObs(jndx) > 0) then
			SNOFCS_atObs_ens(:, jndx) = SNOFCS_Inp_Ens(:, index_back_atObs(jndx))
		endif
	end do

	if ((p_tN==2) .and. (p_tRank==0) .and. print_deb) then 
		PRINT*, "Background Indices at eval points"
		PRINT*, index_back_atEval
		PRINT*, "Background Indices at obs points"
		PRINT*, index_back_atObs	
		PRINT*, "Obs at Eval Points" 
		PRINT*, Obs_atEvalPts	
		PRINT*, "Forecast at Eval Points"
		PRINT*, SNOFCS_atEvalPts	
		PRINT*, "Lat at Eval Points"
		PRINT*, Lat_atEvalPts
		PRINT*, "Lon at Eval Points"
		PRINT*, Lon_atEvalPts
	endif
	!Stop
	! call Observation_Operator(RLA, RLO, OROG, Lat_GHCND, Lon_GHCND,   &
	! 					LENSFC, Num_Ghcnd, max_distance,            &
	! 					SNOFCS,                                  &
	! 					SNOFCS_atGHCND, Ele_GHCND, index_back_atGHCND)
	if ((p_tN==2) .and. (p_tRank==0) .and. print_deb) then
		PRINT*, "GHCND Lat range from rank: ", MYRANK, MINVAL(Lat_GHCND), " ", MAXVAL(Lat_GHCND)
		PRINT*, "GHCND Lon range from rank: ", MYRANK, MINVAL(Lon_GHCND), " ", MAXVAL(Lon_GHCND)
		PRINT*, "Background at GHCND locations from rank: ", MYRANK
		PRINT*, SNOFCS_atGHCND   !SNOFCS_atObs_ens   !	
		PRINT*, "Elevation at GHCND locations from rank: ", MYRANK
		PRINT*, Ele_GHCND
		PRINT*,'Finished observation operator for GHCND'	
	endif
	if (myrank==0) PRINT*,'Finished observation operator for GHCND'		
	ims_threshold = 0.5  ! threshold for converting IMS fSCA to binary 1, 0 values
	! Call Observation_Operator_IMS_fSCA_Threshold(SNCOV_IMS, SNOFCS, SNWDEN, assim_SWE,	   &
	! 						LENSFC, ims_threshold, SNWD_IMS_at_Grid) 
	call Observation_Operator_IMS_fSCA(SNCOV_IMS, SNWDEN, VETFCS, assim_SWE, LENSFC, SNWD_IMS_at_Grid) 
	if ((p_tN==2) .and. (p_tRank==0) .and. print_deb) then
		PRINT*, "IMS obs at each grid cell from rank: ", MYRANK
		PRINT*, SNWD_IMS_at_Grid
		PRINT*,'Finished observation operator for IMS'
	endif
	if (myrank==0) PRINT*,'Finished observation operator for IMS'
	
	L_horz = 55.  !120.  !
	h_ver = 800.  !1200.	!
	obs_tolerance = 5.0
	max_num_loc = 50   !100	!
	ims_max_ele = 1500.
	ims_assm_hour = 18
	max_distance = 250.   !Km 120.			!
	!if (myrank==4) then 
	if (myrank==0) PRINT*,'Starting DA loop'

	! 10.21.20: no assimilation if land-ice
	if (IVEGSRC == 2) then   ! sib
		veg_type_landice=13
	else
		veg_type_landice=15
	endif
	!Do ncol=1, IDIM
	!Do nrow = 1, Num_Snotel !JDIM/2
	    Do jndx = mp_start, mp_end     !1, LENSFC		!jndx = (ncol-1)*JDIM + nrow		! 
			num_loc_1 = 0
			assim_IMS_thisGridCell = .FALSE.
			call debug_print("loop ", float(jndx))
			! GHCND			 
			if(assim_GHCND) then
				! 8.19.20: we are assuming here if deterministic forecast exists/not null
				! at a point, then ensemble members also have valid value
				call nearest_Observations_Locations(RLA(jndx), RLO(jndx),    &
						Lat_GHCND, Lon_GHCND,  Num_Ghcnd, max_distance, max_num_loc,   &
						Stdev_back, Stdev_Obs_depth, obs_tolerance,                 &
						SNOFCS_atGHCND, SNOObs_atGHCND,						 &
						loc_nearest_Obs,  num_loc_1) !,      &LENSFC,
				call debug_print("number of GHCND sndpth obs ", float(num_loc))
			endif
			num_loc = num_loc_1
			!check IMS is assimilated
			if (assim_IMS) then
				if((.NOT. IEEE_IS_NAN(SNWD_IMS_at_Grid(jndx))) .AND. &
				   (OROG(jndx) <= ims_max_ele) .AND. &
				   (IH == ims_assm_hour)) then
					num_loc = num_loc + 1
					assim_IMS_thisGridCell = .TRUE.
				endif
			endif
            ! if assim_IMS=false >> num_loc_1=num_loc
			! 10.21.20: no assimilation if land-ice	.NE. veg_type_landice	
			if((num_loc > 0) .and. (NINT(VETFCS(jndx)) /= veg_type_landice)) then ! .and. (SNCOV_IMS(jndx) > 0.)) then        
				Allocate(obs_Array(num_loc))
				Allocate(Lat_Obs(num_loc))
				Allocate(Lon_Obs(num_loc))
				Allocate(Ele_Obs(num_loc))
				! ghcnd
				if(num_loc_1 > 0) then
					Do zndx = 1, num_loc_1     
						obs_Array(zndx) = SNOObs_atGHCND(loc_nearest_Obs(zndx))
						Lat_Obs(zndx) = Lat_GHCND(loc_nearest_Obs(zndx))
						Lon_Obs(zndx) = Lon_GHCND(loc_nearest_Obs(zndx))
						Ele_Obs(zndx) = Ele_GHCND(loc_nearest_Obs(zndx))
					End Do
				End if
				!ims
				if(assim_IMS_thisGridCell) then
					obs_Array(num_loc) = SNWD_IMS_at_Grid(jndx)
					Lat_Obs(num_loc) = RLA(jndx)   !Lat_IMS_atGrid(jndx)
					Lon_Obs(num_loc) = RLO(jndx)   !Lon_IMS_atGrid(jndx)
					Ele_Obs(num_loc) = OROG(jndx)  !Ele_IMS(jndx)
				endif	
				! EnKF
				Allocate(obs_Innov_ens(num_loc, ens_size))				
				Call snow_DA_EnKF_LocCov(RLA(jndx), RLO(jndx), OROG(jndx),     &
						Lat_Obs, Lon_Obs, Ele_Obs,     				&
						L_horz, h_ver,                                      &   !L_horz in Km, h_ver in m
					assim_IMS_thisGridCell,                    &
					jndx, ens_size, LENSFC, SNOFCS_Inp_Ens,          &
					num_loc_1, num_loc, loc_nearest_Obs, SNOFCS_atObs_ens,	 &
					Stdev_Obs_depth, Stdev_Obs_ims,                  &
					obs_Array,                          &
					obs_Innov_ens, innov_at_Grid_ens(jndx,:), anl_at_Grid_ens(jndx,:))
				! call snow_DA_EnKF(assim_IMS_thisGridCell,                    &
				! 	jndx, ens_size, LENSFC, SNOFCS_Inp_Ens,          &
				! 	num_loc_1, num_loc, loc_nearest_Obs, SNOFCS_atObs_ens,	 &
				! 	Stdev_Obs_depth, Stdev_Obs_ims,                  &
				! 	obs_Array,                          &
				! 	obs_Innov_ens, innov_at_Grid_ens(jndx,:), anl_at_Grid_ens(jndx,:))
				if ((p_tN==4) .and. (p_tRank==0) .and. print_deb) then	
					print*, "proc ", myrank, "loop ", jndx, &
					"num depth obs ", num_loc_1, "total obs", num_loc, " ens size ", ens_size
					PRINT*, "Ens background at obs pts: "
					PRINT*, SNOFCS_Inp_Ens(:, jndx)	
					PRINT*, "Observed"
					PRINT*,  obs_Array
					PRINT*, "Ens Obs innovation: "
					PRINT*, obs_Innov_ens
					PRINT*, "Ens innovation at grid: "
					PRINT*, innov_at_Grid_ens(jndx, :)
					PRINT*, "Ens analysis at grid: "
					PRINT*, anl_at_Grid_ens(jndx, :)
				endif					
				DEALLOCATE(obs_Array, obs_Innov_ens)
				DEALLOCATE(Lat_Obs, Lon_Obs, Ele_Obs)
				! QCC by ims--use ims as snow mask
				Do zndx = 1, ens_size+1  
					if((SNCOV_IMS(jndx) >= 0.5) .and. (anl_at_Grid_ens(jndx, zndx) < 50.)) then
						anl_at_Grid_ens(jndx, zndx) = 50.
					! elseif((SNCOV_IMS(jndx) < 0.5) .and. (SNCOV_IMS(jndx) > 0.1) .and. (anl_at_Grid(jndx) >= 50)) then
					! 	anl_at_Grid(jndx) = 50.
					elseif((SNCOV_IMS(jndx) < 0.1) .and. (anl_at_Grid_ens(jndx, zndx) >= 0.)) then
						anl_at_Grid_ens(jndx, zndx) = 0.
					endif
				End do
				if ((p_tN==4) .and. (p_tRank==0) .and. print_deb) then						
					PRINT*, "Ens analysis at grid after mask: "
					PRINT*, anl_at_Grid_ens(jndx, :)
				endif
			else
				anl_at_Grid_ens(jndx, :) = SNOFCS_Inp_Ens(jndx, :)
			endif
			if (assim_GHCND) Deallocate(loc_nearest_Obs) 
		End do
	!End do
	if (myrank==0) PRINT*, 'Finished DA loops'

! ToDO: Better way to handle this?
! Real data type size corresponding to mpi
	rsize = SIZEOF(snwden_val) 
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
	! ens analysis innov_at_Grid_ens(jndx), anl_at_Grid_ens(jndx)
	if (MYRANK > (MAX_TASKS - 1) ) then
		call MPI_SEND(anl_at_Grid_ens(mp_start:mp_end,:),N_sA*(ens_size+1), mpiReal_size,    &
					  p_tN, MYRANK*1000, MPI_COMM_WORLD, IERR) 
		call MPI_SEND(innov_at_Grid_ens(mp_start:mp_end,:),N_sA*(ens_size+1), mpiReal_size,    &
					  p_tN, MYRANK*10000, MPI_COMM_WORLD, IERR)
	else    !if(p_tRank == 0) then  
		Do pindex =  1, (Np_til - 1)   ! sender proc index within tile group
			dest_Aoffset = pindex * N_sA + N_sA_Ext + 1   ! dest array offset
			send_proc = MYRANK +  pindex * MAX_TASKS
			call MPI_RECV(anl_at_Grid_ens(dest_Aoffset:dest_Aoffset+N_sA-1, :), N_sA*(ens_size+1),      &
			mpiReal_size, send_proc, send_proc*1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
			call MPI_RECV(innov_at_Grid_ens(dest_Aoffset:dest_Aoffset+N_sA-1, :), N_sA*(ens_size+1),       &
			mpiReal_size, send_proc, send_proc*10000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
		enddo
	endif
	if (myrank==0) PRINT*,'Finished Data copy'

	if (MYRANK > MAX_TASKS - 1 ) goto 998   ! if(p_tRank /= 0 ) goto 998

	! avoid -ve anl
	Where(anl_at_Grid_ens < 0.) anl_at_Grid_ens = 0.
	! swe and snwd
	SNOANL = anl_at_Grid_ens(:,ens_size+1) / SNWDEN		!snwden_val 
	SWDANL = anl_at_Grid_ens(:,ens_size+1)	
!    if(assim_SWE) then 
!       SNOANL_Cur = SNOANL_Cur 
!    else
! 		SNOANL_Cur = SWDANL_Cur
!    end if
	if (print_deb) then
		PRINT*, "Innovation SWE/snwd from rank: ", MYRANK
	    PRINT*, innov_at_Grid_ens	
	    PRINT*, "Analysis SWE/ snwd  from rank: ", MYRANK
	    PRINT*, anl_at_Grid_ens	
	endif

	!Compute updated snocov	
	!Call update_snow_cover_fraction(LENSFC, SNOANL, VETFCS, anl_fSCA)

	! copy values at eval points
	innov_atEvalPts = IEEE_VALUE(innov_atEvalPts, IEEE_QUIET_NAN)
	SNOANL_atEvalPts = IEEE_VALUE(SNOANL_atEvalPts, IEEE_QUIET_NAN)
	SNOANL_Cur_atEvalPts = IEEE_VALUE(SNOANL_Cur_atEvalPts, IEEE_QUIET_NAN)
	Do jndx = 1, num_Eval
		if ((index_back_atEval(jndx) > 0) .and. &
		    (NINT(VETFCS(jndx)) /= veg_type_landice)) then  !10.21.20: exclude if land-ice
				innov_atEvalPts(jndx) = innov_at_Grid_ens(index_back_atEval(jndx), ens_size+1)
				SNOANL_atEvalPts(jndx) = anl_at_Grid_ens(index_back_atEval(jndx), ens_size+1)
				SNOANL_Cur_atEvalPts(jndx) = SWDANL_Cur(index_back_atEval(jndx)) !SNOANL_Cur
		endif
	End do

	! write outputs	
	Write(rank_str, '(I3.3)') (MYRANK+1)
	! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/Analysis
	da_out_file = "./SNOANLEnKF."// &  !
				  TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18_tile"//rank_str//".nc"  !	
	call Write_DA_Outputs(da_out_file, IDIM, JDIM, LENSFC, MYRANK, &
	                      SWEFCS, SNOANL, SNOFCS, SWDANL, &  !
						  innov_at_Grid_ens(:,ens_size+1), SNCOV_IMS, &
			      !   Num_Ghcnd, Lat_GHCND, Lon_GHCND, SNOObs_atGHCND, &
			    !   SNOFCS_atGHCND, SNOFCS_atGHCND, SNOFCS_atGHCND)
			num_Eval, Lat_atEvalPts, Lon_atEvalPts, Obs_atEvalPts, & 
		SNOFCS_atEvalPts, innov_atEvalPts, SNOANL_atEvalPts, SNOANL_Cur_atEvalPts)  !, anl_fSCA) !updated snocov

998 CONTINUE
	DEALLOCATE(SNOObs_atGHCND, SNOFCS_atGHCND, Lat_GHCND, Lon_GHCND, Ele_GHCND)  !, index_back_atGHCND)
	DEALLOCATE(Obs_atEvalPts, SNOFCS_atEvalPts, innov_atEvalPts, SNOANL_atEvalPts)
	DEALLOCATE(SNOANL_Cur_atEvalPts)
	DEALLOCATE(index_back_atObs, index_back_atEval, Lat_atEvalPts, Lon_atEvalPts) !, Ele_atEvalPts)
	DEAllocate(SNOFCS_Inp_Ens, SNOFCS_atObs_ens)
	DEAllocate(innov_at_Grid_ens, anl_at_Grid_ens)
	!DEALLOCATE(SNCOV_IMS, Lat_IMS, Lon_IMS, SNOFCS_atIMS, SWDFCS_atIMS)
999 CONTINUE
    PRINT*,'Finished EnKF DA ON RANK: ', MYRANK
	CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

	!STOP

	RETURN

 END subroutine Snow_Analysis_EnKF 

 subroutine Snow_Analysis_OI_multiTime(SNOW_IO_TYPE, MAX_TASKS, MYRANK, NPROCS, IDIM, JDIM, & 
	                       LENSFC, IVEGSRC, GHCND_SNOWDEPTH_PATH, IMS_SNOWCOVER_PATH, & 
	                    IMS_INDEXES_PATH, CURRENT_ANALYSIS_PATH, SNOANL) !IY, IM, ID, IH, 
							
	!----------------------------------------------------------------------
	! Input: forecast/background states for a single tile by a single MPI process
	! reads observations: snow depth / SWE and cover
	! calls observation operator for a given tile 
	! does OI update per tile
	! returns updated/Analysis states back to the caller, surface drive (sfcdrv) program
	!
	! RLA, RLO: lat lon information for the tile
	! SNOFCS(LENSFC), SWDFCS(LENSFC): snowdepth and snow water equivalent forecat (background)
	! (background/forecast in arrays of LENSFC)
	! compute snow depnsity (SNWDEN) and fraction of snow (if needed) from SWE/SNWD 
	!
	! IDIM * JDIM = LENSFC: number of grid cells in tile = xdim * ydim   
	! IY, IM, ID, IH = year, month, day, hour of current model step   
	! MYRANK: rank/id of the MPI process
	!
	! Outputs:
	! SNOANL, SWDANL:  SWE and snowdepthafter DA
	!                   
	!----------------------------------------------------------------------
	IMPLICIT NONE
	!
	include 'mpif.h'
	
	integer, parameter :: dp = kind(1.d0)

	INTEGER, intent(in) :: SNOW_IO_TYPE, MAX_TASKS, MYRANK, NPROCS, IDIM, JDIM, LENSFC, IVEGSRC
	CHARACTER(LEN=*), Intent(In)   :: GHCND_SNOWDEPTH_PATH, IMS_SNOWCOVER_PATH, & 
									  IMS_INDEXES_PATH, CURRENT_ANALYSIS_PATH
!  GHCND_SNOWDEPTH_PATH="/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/GHCND/",
!  IMS_SNOWCOVER_PATH="/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS/",
!  IMS_INDEXES_PATH="/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS_INDEXES/",
!  CURRENT_ANALYSIS_PATH = "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/Cur_Analysis/"
	INTEGER             :: IY, IM, ID, IH
	INTEGER			    :: dTIncrement
	REAL                :: DELTSFC, IH_real
	REAL, intent(Out)   :: SNOANL(LENSFC)   !, SWDANL(LENSFC), anl_fSCA(LENSFC)
	CHARACTER(LEN=5)    :: TILE_NUM
	Character(LEN=3)    :: rank_str
	INTEGER			    :: IERR	
	REAL        :: RLA(LENSFC), RLO(LENSFC), RLO_Tile(LENSFC), OROG(LENSFC), OROG_UF(LENSFC)
	REAL                :: SNOFCS_Inp(LENSFC), SWDFCS(LENSFC)
	REAL                :: SNOFCS(LENSFC), SNWDEN(LENSFC), VETFCS(LENSFC), SWDANL(LENSFC)	
	CHARACTER(len=250)   :: dim_name, ghcnd_inp_file, ims_inp_file, ims_inp_file_indices
	CHARACTER(len=250) 	 :: forc_inp_path, anl_inp_path, da_out_file
	CHARACTER(len=5)     :: y_str, m_str, d_Str, h_str, fvs_tile
	REAL, ALLOCATABLE    :: SNWD_GHCND(:), SNOFCS_atGHCND(:)		!, SWDFCS_atGHCND(:)
	REAL, ALLOCATABLE    :: Lat_GHCND(:), Lon_GHCND(:), Ele_GHCND(:)  
	REAL                 :: lat_min, lat_max, lon_min, lon_max  	
	Real		     :: SNCOV_IMS(LENSFC)  ! ims resampled at each grid
	Real		     :: SNWD_IMS_at_Grid(LENSFC)

	INTEGER :: Num_Ghcnd, Num_Ims, num_Eval, num_sub !Num_Ims_Lat, Num_Ims_Lon
	Real	:: max_distance   ! radius_of_influence for selecting state at observation point
	INTEGER :: jndx, zndx, ncol, nrow
	Integer, Allocatable   :: Loc_backSt_atObs(:)
	Integer				   :: num_loc, num_loc_1, max_num_loc
	Real, Parameter		:: Stdev_back = 30., Stdev_Obs_depth = 40., Stdev_Obs_ims = 80. ! mm 
	Integer				:: ims_assm_hour
	Real				:: obs_tolerance, ims_max_ele
	Real                :: L_horz, h_ver
	Real(dp), Allocatable 	 :: B_cov_mat(:,:), b_cov_vect(:)
	Real(dp), Allocatable 	 :: O_cov_mat(:,:), W_wght_vect(:)
	Real, Allocatable   :: back_at_Obs(:), obs_Array(:), Lat_Obs(:), Lon_Obs(:), Ele_Obs(:)
	REAL                :: innov_at_Grid(LENSFC), anl_at_Grid(LENSFC)
	Real, Allocatable  :: obs_Innov(:)
	Real               :: SNOANL_Cur(LENSFC), SWDANL_Cur(LENSFC)
	REAL, ALLOCATABLE  :: SNOANL_Cur_atEvalPts(:)  !evalution points 
	REAL, ALLOCATABLE  :: SNOFCS_atEvalPts(:), innov_atEvalPts(:), SNOANL_atEvalPts(:)  !evalution points 
	REAL, ALLOCATABLE  :: Lat_atEvalPts(:), Lon_atEvalPts(:), Obs_atEvalPts(:)     !evalution points
	Integer, ALLOCATABLE 	:: index_back_atEval(:)     ! background locations at eval points 
	Integer, ALLOCATABLE	:: index_back_atObs(:)   ! the location of background corresponding obs

	Real			   :: ims_threshold      ! threshold for converting IMS fSCA to binary 1, 0 values	
	Real			   :: snwden_val
	LOGICAL 		   :: assim_SWE, assim_GHCND   !assimilate swe? (instead of snow depth)
	LOGICAL 		   :: assim_IMS, assim_IMS_thisGridCell    !assimilate sncov

	Integer 		   :: veg_type_landice ! 10.21.20: no assmn over land ice
    ! for mpi par
	INTEGER            :: Np_ext, Np_til, p_tN, p_tRank, N_sA, N_sA_Ext, mp_start, mp_end
	INTEGER            :: send_proc, rec_stat(MPI_STATUS_SIZE), dest_Aoffset, pindex
	INTEGER            :: mpiReal_size, rsize

	assim_SWE = .False.  ! note: if this is set true, may need to adjust background and obs stddev above
	If (SNOW_IO_TYPE == 2) assim_SWE = .True.
	assim_GHCND = .True.
	assim_IMS = .True.     !!assimilate sncov; 
	assim_IMS_thisGridCell = .FALSE.    ! if assimilating ims, skip this grid cell for this time step

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
	mp_end = (p_tRank + 1) * N_sA + N_sA_Ext 		! end index of subarray for proc	
	If(myrank == 0 ) PRINT*,"sub array length ", N_sA, " extra sub array: ", N_sA_Ext
! if (p_tN /= 2 ) goto 999

	if( (GHCND_SNOWDEPTH_PATH(1:8).eq.'        ').and. & 
	    (IMS_SNOWCOVER_PATH(1:8).eq.'        ') ) then
		print*, "Observation paths don't exist!, skipping the OI DA"
		goto 999
	end if

	! READ THE INPUT SURFACE DATA ON THE CUBED-SPHERE TILE.
	Call READ_Forecast_Data(p_tN, LENSFC, SNOFCS_Inp, SWDFCS, VETFCS)  !VEGFCS,
	! READ THE OROGRAPHY AND GRID POINT LAT/LONS FOR THE CUBED-SPHERE TILE.
	CALL READ_LAT_LON_OROG_atRank(p_tN, RLA,RLO,OROG,OROG_UF,TILE_NUM,IDIM,JDIM,LENSFC)
!!PRINT*,"Snow anl on ", MYRANK, " Tile group: ", p_tN, " Tile: ", TILE_NUM
	RLO_Tile = RLO
	Where(RLO_Tile > 180) RLO_Tile = RLO_Tile - 360
	lat_min = MAX(MINVAL(RLA) - 1., -90.)
	lat_max = MIN(MAXVAL(RLA) + 1., 90.)
	lon_min = MAX(MINVAL(RLO_Tile) - 1., -180.)
	lon_max = MIN(MAXVAL(RLO_Tile) + 1., 180.)	
	! lon_min = MAX(MINVAL(RLO) - 1., 0.)
	! lon_max = MIN(MAXVAL(RLO) + 1., 360.)	
	if (p_tN==3)  then  
		lon_min = 125.       ! lon_min is left, lon_max is right, not necessary min/max value
		lon_max = -145.
	endif	
	if ((p_tRank==0) .and. print_deb) then
		print*, "Tile ", p_tN, " min/max lat/lon ", lat_min, lat_max, lon_min, lon_max
	endif

	IY = 2019; IM = 11; ID = 29; IH = 18; IH_real = 18.; 
	DELTSFC  = 24.0          !* dTIncrement
	
Do dTIncrement = 1, 63

	Call UPDATEtime(IY, IM, ID, IH_real, DELTSFC)
	! IH = INT(IH_real) this should be always 18
	
    ! ToDO: use file location the way the rest of the program accesses
	! 4.9.20 for now hard code the file location
	write(y_str, "(I4)") IY
	write(m_str, "(I0.2)") IM
	write(d_str, "(I0.2)") ID 
	write(h_str, "(I0.2)") IH
    write(fvs_tile, "(I3)") IDIM	
	!write(hprev_str, "(I0.2)") (IH - INT(DELTSFC))

	! forc_inp_path = "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/ROTDIRS/dec15/"// &
	!                "gdas."//TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"/atmos/"//TRIM(h_str)// &
	! 			    !"/mem"//ens_str//"/RESTART/"// &
	forc_inp_path = "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/fnbgsio/"// &
					"snow."//TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"." // &
					TRIM(h_str)// "0000.sfc_data."//TILE_NUM//".nc"
	Call READ_Forecast_Data_atPath(forc_inp_path, p_tN, LENSFC, SNOFCS_Inp, SWDFCS) !, VETFCS) !VEGFCS, !SRFLAG)
! 10.22.20 Cur_Analysis: exiting SNODEP based analysis
	!"/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/Cur_Analysis/"
	anl_inp_path = TRIM(CURRENT_ANALYSIS_PATH)//"snow."// &
				TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"." // &
				TRIM(h_str)// "0000.sfcanl_data."//TILE_NUM//".nc"
	Call READ_Analysis_Data(anl_inp_path, p_tN, LENSFC, SNOANL_Cur, SWDANL_Cur)	

	! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/GHCND/
	ghcnd_inp_file = TRIM(GHCND_SNOWDEPTH_PATH)//"GHCND.SNWD."// &
					 TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18.nc"
	dim_name = "Site_Id"	
	call Observation_Read_GHCND_Tile_excNaN(p_tN, ghcnd_inp_file, dim_name, &
					lat_min, lat_max, lon_min, lon_max, & 
					Num_Ghcnd, SNWD_GHCND,              & !SNOObs_atGHCND,		&
					Lat_GHCND, Lon_GHCND, MYRANK)  !Ele_GHCND,		&
	! Do jndx=0, 5
	! 	if ((p_tN == jndx) .and. (p_tRank==0)) then
	! 		print*, "Tile ", p_tN, " num. GHCND obs ", Num_Ghcnd
	! 		PRINT*, "GHCND SNWD from rank: ", MYRANK
	! 		PRINT*, SNOObs_atGHCND
	! 	endif
	! End do
	! Stop
	if ((p_tRank==0) .and. print_deb) then
		print*, "Tile ", p_tN, " num. GHCND obs ", Num_Ghcnd
	endif
	if ((p_tRank==0) .and. (p_tN==2) .and. print_deb) then
		PRINT*, "GHCND SNWD from rank: ", MYRANK
		PRINT*, SNWD_GHCND
		PRINT*, "Lat at GHCND from rank: ", MYRANK
		PRINT*, Lat_GHCND
		PRINT*, "Lon at GHCND from rank: ", MYRANK
		PRINT*, Lon_GHCND
		PRINT*,'Finished reading GHCND'
	endif
	if (myrank==0) PRINT*,'Finished reading GHCND'
   
	! (max) number of IMS subcells within a tile grid cell
	If (IDIM == 96) then          
		num_sub = 627               
		max_distance = 240.		!Km radius of influence: distance from gridcell to search for observations
	elseif (IDIM == 128) then
	   num_sub = 627               
	   max_distance = 240.		!Km 
	elseif (IDIM == 192) then
		!num_sub = 30 
		PRINT*,'Error, tile type not known '
		stop
	elseif (IDIM == 384) then
		!num_sub = 30  
		PRINT*,'Error, tile type not known '
		stop      
	elseif (IDIM == 768) then
		num_sub = 30
		max_distance = 27.        			!Km  
	else
		PRINT*,'Error, tile type not known '
		stop   
	endif
	! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS/"
	ims_inp_file = TRIM(IMS_SNOWCOVER_PATH)//"IMS.SNCOV."// &
				   TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18.nc"                      !
	! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS_INDEXES/"
	ims_inp_file_indices = TRIM(IMS_INDEXES_PATH)//"C"//TRIM(ADJUSTL(fvs_tile))//&
						   ".IMS.Indices."//TRIM(TILE_NUM)//".nc"				
	Call Observation_Read_IMS_Full(ims_inp_file, ims_inp_file_indices, &
				                   MYRANK, JDIM, IDIM, num_sub, SNCOV_IMS)
	if((p_tRank==0) .and. (p_tN==2) .and. print_deb) then
		PRINT*, "IMS SNCOV from rank: ", MYRANK
		PRINT*, SNCOV_IMS
		PRINT*,'Finished reading IMS'
	endif 
	if (myrank==0) PRINT*,'Finished reading IMS'
	
	! 4.8.20: compute snow density ratio from forecast swe and snwdepth
	SNWDEN = SWDFCS / SNOFCS_Inp
	snwden_val = SUM(SNWDEN, Mask = (.not. IEEE_IS_NAN(SNWDEN))) &
	                 / COUNT (.not. IEEE_IS_NAN(SNWDEN))
	SNWDEN = snwden_val  !10.   ! ratio of sndepth to swe when swe /=0 : SWDFCS / SNOFCS
	if(p_tRank == 0 .and. print_deb) then
		print*, "Process ", MYRANK, " average snwd/swe ratio: ", snwden_val
	endif
	if (assim_SWE) then 
		SNWD_GHCND = SNWD_GHCND / SNWDEN		!snwden_val
		SNOFCS = SNOFCS_Inp
	else
		SNOFCS = SNOFCS_Inp * SNWDEN    !snwden_val  !SWDFCS			! assimilate depth
		! SNWD_GHCND = SNWD_GHCND 
	endif

	! Get model states at obs points
	ALLOCATE(SNOFCS_atGHCND(Num_Ghcnd))
	ALLOCATE(Ele_GHCND(Num_Ghcnd)) 
	ALLOCATE(index_back_atObs(Num_Ghcnd)) 
	num_Eval = floor(0.05 * Num_Ghcnd)      ! using 5% of ghcnd locations for evaluation
	ALLOCATE(index_back_atEval(num_Eval)) 
	ALLOCATE(Obs_atEvalPts(num_Eval)) 
	ALLOCATE(SNOFCS_atEvalPts(num_Eval)) 
	ALLOCATE(Lat_atEvalPts(num_Eval))
	ALLOCATE(Lon_atEvalPts(num_Eval)) 
	ALLOCATE(innov_atEvalPts(num_Eval))
	ALLOCATE(SNOANL_atEvalPts(num_Eval))
	ALLOCATE(SNOANL_Cur_atEvalPts(num_Eval))  	 	
	if(p_tRank == 0) then 
		PRINT*, "Tile ", p_tN+1, " ", num_Eval, ' points for evaluation excluded from DA'	
	endif	

	Call Observation_Operator_Parallel(Myrank, MAX_TASKS, p_tN, p_tRank, Np_til, & 
	                    RLA, RLO, OROG, Lat_GHCND, Lon_GHCND,   &
						LENSFC, Num_Ghcnd, num_Eval, max_distance, SNOFCS, SNWD_GHCND,  &
						SNOFCS_atGHCND, Ele_GHCND, index_back_atObs, index_back_atEval, &
						Obs_atEvalPts, SNOFCS_atEvalPts, Lat_atEvalPts, Lon_atEvalPts)
	if ((p_tN==2) .and. (p_tRank==0) .and. print_deb) then 
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
	!Stop
	! call Observation_Operator(RLA, RLO, OROG, Lat_GHCND, Lon_GHCND,   &
	! 					LENSFC, Num_Ghcnd, max_distance,            &
	! 					SNOFCS,                                  &
	! 					SNOFCS_atGHCND, Ele_GHCND, index_back_atGHCND)
	if ((p_tN==2) .and. (p_tRank==0) .and. print_deb) then
		PRINT*, "GHCND Lat range from rank: ", MYRANK, MINVAL(Lat_GHCND), " ", MAXVAL(Lat_GHCND)
		PRINT*, "GHCND Lon range from rank: ", MYRANK, MINVAL(Lon_GHCND), " ", MAXVAL(Lon_GHCND)
		PRINT*, "Background SWE at GHCND locations from rank: ", MYRANK
		PRINT*, SNOFCS_atGHCND	
		PRINT*, "Elevation at GHCND locations from rank: ", MYRANK
		PRINT*, Ele_GHCND
		PRINT*,'Finished observation operator for GHCND'	
	endif
	if (myrank==0) PRINT*,'Finished observation operator for GHCND'		
	ims_threshold = 0.5  ! threshold for converting IMS fSCA to binary 1, 0 values
	! Call Observation_Operator_IMS_fSCA_Threshold(SNCOV_IMS, SNOFCS, SNWDEN, assim_SWE,	   &
	! 						LENSFC, ims_threshold, SNWD_IMS_at_Grid) 
	call Observation_Operator_IMS_fSCA(SNCOV_IMS, SNWDEN, VETFCS, assim_SWE, LENSFC, SNWD_IMS_at_Grid) 
	if ((p_tN==2) .and. (p_tRank==0) .and. print_deb) then
		PRINT*, "IMS obs at each grid cell from rank: ", MYRANK
		PRINT*, SNWD_IMS_at_Grid
		PRINT*,'Finished observation operator for IMS'
	endif
	if (myrank==0) PRINT*,'Finished observation operator for IMS'
	
	L_horz = 55.  !120.  !
	h_ver = 800.  !1200.	!
	obs_tolerance = 5.0
	max_num_loc = 50   !100	!
	ims_max_ele = 1500.
	ims_assm_hour = 18
	max_distance = 250.   !120.			!Km
	!if (myrank==4) then 
	if (myrank==0) PRINT*,'Starting DA loop'

	! 10.21.20: no assimilation if land-ice
	if (IVEGSRC == 2) then   ! sib
		veg_type_landice=13
	else
		veg_type_landice=15
	endif
	!Do ncol=1, IDIM
	!Do nrow = 1, Num_Snotel !JDIM/2
	    Do jndx = mp_start, mp_end     !1, LENSFC		!jndx = (ncol-1)*JDIM + nrow		! 
			num_loc_1 = 0
			assim_IMS_thisGridCell = .FALSE.
			call debug_print("loop ", float(jndx))
			! GHCND			 
			if(assim_GHCND) then
				call nearest_Observations_Locations(RLA(jndx), RLO(jndx),    &
						Lat_GHCND, Lon_GHCND,  Num_Ghcnd, max_distance, max_num_loc,   &
						Stdev_back, Stdev_Obs_depth, obs_tolerance,                 &
						SNOFCS_atGHCND, SNWD_GHCND,						 &
						Loc_backSt_atObs,  num_loc_1) !,      &LENSFC,
				call debug_print("number of GHCND sndpth obs ", float(num_loc))
			endif
			num_loc = num_loc_1
			!check IMS is assimilated
			if (assim_IMS) then
				if((.NOT. IEEE_IS_NAN(SNWD_IMS_at_Grid(jndx))) .AND. &
				   (OROG(jndx) <= ims_max_ele) .AND. &
				   (IH == ims_assm_hour)) then
					num_loc = num_loc + 1
					assim_IMS_thisGridCell = .TRUE.
				endif
			endif 
			! if assim_IMS=false >> num_loc_1=num_loc
			! 10.21.20: no assimilation if land-ice	.NE. veg_type_landice	
			if((num_loc > 0) .and. (NINT(VETFCS(jndx)) /= veg_type_landice)) then ! .and. (SNCOV_IMS(jndx) > 0.)) then      
				! get background states
				Allocate(back_at_Obs(num_loc))
				Allocate(obs_Array(num_loc))
				Allocate(Lat_Obs(num_loc))
				Allocate(Lon_Obs(num_loc))
				Allocate(Ele_Obs(num_loc))
				! ghcnd
				if(num_loc_1 > 0) then
					Do zndx = 1, num_loc_1     
						back_at_Obs(zndx) = SNOFCS_atGHCND(Loc_backSt_atObs(zndx))
						obs_Array(zndx) = SNWD_GHCND(Loc_backSt_atObs(zndx))
						Lat_Obs(zndx) = Lat_GHCND(Loc_backSt_atObs(zndx))
						Lon_Obs(zndx) = Lon_GHCND(Loc_backSt_atObs(zndx))
						Ele_Obs(zndx) = Ele_GHCND(Loc_backSt_atObs(zndx))
					End Do
				End if
				!ims
				if(assim_IMS_thisGridCell) then
					back_at_Obs(num_loc) = SNOFCS(jndx)
					obs_Array(num_loc) = SNWD_IMS_at_Grid(jndx)
					Lat_Obs(num_loc) = RLA(jndx)   !Lat_IMS_atGrid(jndx)
					Lon_Obs(num_loc) = RLO(jndx)   !Lon_IMS_atGrid(jndx)
					Ele_Obs(num_loc) = OROG(jndx)  !Ele_IMS(jndx)
				endif
				! compute covariances
				Allocate(B_cov_mat(num_loc, num_loc))
				Allocate(b_cov_vect(num_loc))
				Allocate(O_cov_mat(num_loc, num_loc))
				Allocate(W_wght_vect(num_loc))
				call compute_covariances(RLA(jndx), RLO(jndx), OROG(jndx), SNOFCS(jndx),    &
					Lat_Obs, Lon_Obs, Ele_Obs, num_loc,    			 &
					Stdev_back, Stdev_Obs_depth, Stdev_Obs_ims,      &
					L_horz, h_ver,                                   &   !L_horz in Km, h_ver in m
					assim_IMS_thisGridCell,                          &
					B_cov_mat, b_cov_vect, O_cov_mat, W_wght_vect)
				! call OI DA
				Allocate(obs_Innov(num_loc))
				call Snow_DA_OI(back_at_Obs, obs_Array, num_loc, W_wght_vect,            &
					SNOFCS(jndx), innov_at_Grid(jndx), anl_at_Grid(jndx), obs_Innov)
				if ((p_tN==2) .and. (p_tRank==0) .and. print_deb) then	
					print*, "proc ", myrank, "loop ", jndx, "num depth obs ", num_loc_1, "total obs", num_loc
					PRINT*, " background at obs pts: "
					PRINT*, back_at_Obs	
					PRINT*, "Observed"
					PRINT*,  obs_Array
					PRINT*, "Obs innovation: "
					PRINT*, obs_Innov
					PRINT*, "Weight vector: "
					PRINT*, W_wght_vect	
					print*, "innov: ", innov_at_Grid(jndx), "forec: ", SNOFCS(jndx), " anl: ", anl_at_Grid(jndx)
				endif		
				!free mem
				DEALLOCATE(back_at_Obs, obs_Array)
				DEALLOCATE(Lat_Obs, Lon_Obs, Ele_Obs, obs_Innov)
				DEALLOCATE(B_cov_mat, b_cov_vect, O_cov_mat, W_wght_vect)
				! QCC by ims--use ims as snow mask
				if((SNCOV_IMS(jndx) >= 0.5) .and. (anl_at_Grid(jndx) < 50.)) then
					anl_at_Grid(jndx) = 50.
				! elseif((SNCOV_IMS(jndx) < 0.5) .and. (SNCOV_IMS(jndx) > 0.1) .and. (anl_at_Grid(jndx) >= 50)) then
				! 	anl_at_Grid(jndx) = 50.
				elseif((SNCOV_IMS(jndx) < 0.1) .and. (anl_at_Grid(jndx) >= 0.)) then
						anl_at_Grid(jndx) = 0.
				endif
				if ((p_tN==2) .and. (p_tRank==0) .and. print_deb) then						
					PRINT*, "analyis at grid after mask: ", anl_at_Grid(jndx)
				endif
			else
				anl_at_Grid(jndx) = SNOFCS(jndx)
			endif
			if (assim_GHCND) Deallocate(Loc_backSt_atObs) 
		End do
	!End do
	if (myrank==0) PRINT*, 'Finished DA loops'

! ToDO: Better way to handle this?
! Real data type size corresponding to mpi
	rsize = SIZEOF(snwden_val) 
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
		call MPI_SEND(anl_at_Grid(mp_start:mp_end), N_sA, mpiReal_size, p_tN,   &
					  MYRANK, MPI_COMM_WORLD, IERR) 
		call MPI_SEND(innov_at_Grid(mp_start:mp_end), N_sA, mpiReal_size, p_tN,   &
					  MYRANK*100, MPI_COMM_WORLD, IERR)
	else    !if(p_tRank == 0) then  
		Do pindex =  1, (Np_til - 1)   ! sender proc index within tile group
			dest_Aoffset = pindex * N_sA + N_sA_Ext + 1   ! dest array offset
			send_proc = MYRANK +  pindex * MAX_TASKS
			call MPI_RECV(anl_at_Grid(dest_Aoffset:dest_Aoffset+N_sA-1), N_sA, mpiReal_size, send_proc,      &
					  send_proc, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
			call MPI_RECV(innov_at_Grid(dest_Aoffset:dest_Aoffset+N_sA-1), N_sA, mpiReal_size, send_proc,      &
					  send_proc*100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
		enddo
	endif
	if (myrank==0) PRINT*,'Finished Data copy'

	if (MYRANK > MAX_TASKS - 1 ) goto 998   ! if(p_tRank /= 0 ) goto 998

	! avoid -ve anl
	Where(anl_at_Grid < 0.) anl_at_Grid = 0.
	! swe and snwd
	if(assim_SWE) then
		SNOANL = anl_at_Grid 
		SWDANL = anl_at_Grid * SNWDEN		!snwden_val 
		!SNOANL_Cur = SNOANL_Cur    
	else
		SNOANL = anl_at_Grid / SNWDEN		!snwden_val 
		SWDANL = anl_at_Grid
		SNOANL_Cur = SWDANL_Cur
		!innov_at_Grid = innov_at_Grid / SNWDEN
	endif	
	if (print_deb) then
		PRINT*, "Innovation SWE/snwd from rank: ", MYRANK
	    PRINT*, innov_at_Grid	
	    PRINT*, "Analysis SWE/ snwd  from rank: ", MYRANK
	    PRINT*, anl_at_Grid	
	endif

	!Compute updated snocov	
	!Call update_snow_cover_fraction(LENSFC, SNOANL, VETFCS, anl_fSCA)

	! copy values at eval points
	innov_atEvalPts = IEEE_VALUE(innov_atEvalPts, IEEE_QUIET_NAN)
	SNOANL_atEvalPts = IEEE_VALUE(SNOANL_atEvalPts, IEEE_QUIET_NAN)
	SNOANL_Cur_atEvalPts = IEEE_VALUE(SNOANL_Cur_atEvalPts, IEEE_QUIET_NAN)
	Do jndx = 1, num_Eval
		if ((index_back_atEval(jndx) > 0) .and. &
		    (NINT(VETFCS(jndx)) /= veg_type_landice)) then  !10.21.20: exclude if land-ice	
				innov_atEvalPts(jndx) = innov_at_Grid(index_back_atEval(jndx))
				SNOANL_atEvalPts(jndx) = anl_at_Grid(index_back_atEval(jndx))
				SNOANL_Cur_atEvalPts(jndx) = SNOANL_Cur(index_back_atEval(jndx))
		endif
	End do

	! write outputs	
	Write(rank_str, '(I3.3)') (MYRANK+1)
	! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/Analysis
	da_out_file = "./SNOANLOI."// &  !
				  TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18_tile"//rank_str//".nc"  !	
	call Write_DA_Outputs(da_out_file, IDIM, JDIM, LENSFC, MYRANK, &
	                      SNOFCS_Inp, SNOANL, SWDFCS, SWDANL, &  !
			      innov_at_Grid, SNCOV_IMS, &
			      !   Num_Ghcnd, Lat_GHCND, Lon_GHCND, SNWD_GHCND, &
			      !   SNOFCS_atGHCND, SNOFCS_atGHCND, SNOFCS_atGHCND)
			      num_Eval, Lat_atEvalPts, Lon_atEvalPts, Obs_atEvalPts, & 
			SNOFCS_atEvalPts, innov_atEvalPts, SNOANL_atEvalPts, SNOANL_Cur_atEvalPts)  !, anl_fSCA) !updated snocov

998 CONTINUE
	DEALLOCATE(SNWD_GHCND, SNOFCS_atGHCND, Lat_GHCND, Lon_GHCND, Ele_GHCND)  !, index_back_atGHCND)
	DEALLOCATE(Obs_atEvalPts, SNOFCS_atEvalPts, innov_atEvalPts, SNOANL_atEvalPts)
	DEALLOCATE(SNOANL_Cur_atEvalPts) 
	DEALLOCATE(index_back_atObs, index_back_atEval, Lat_atEvalPts, Lon_atEvalPts) !, Ele_atEvalPts)
	!DEALLOCATE(SNCOV_IMS, Lat_IMS, Lon_IMS, SNOFCS_atIMS, SWDFCS_atIMS)
	If(MYRANK==0) PRINT*,'Finished OI DA at time: ', y_str, m_str, d_str
End do

999 CONTINUE
    !PRINT*,'Finished OI DA ON RANK: ', MYRANK
	CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

	STOP

	RETURN

 END subroutine Snow_Analysis_OI_multiTime

 subroutine Snow_Analysis_OI(SNOW_IO_TYPE, MAX_TASKS, MYRANK, NPROCS, IDIM, JDIM, IY, IM, ID, IH, &  
							LENSFC, IVEGSRC, GHCND_SNOWDEPTH_PATH, IMS_SNOWCOVER_PATH, & 
							IMS_INDEXES_PATH, CURRENT_ANALYSIS_PATH, SNOANL)
							
	!----------------------------------------------------------------------
	! Input: forecast/background states for a single tile by a single MPI process
	! reads observations: snow depth / SWE and cover
	! calls observation operator for a given tile 
	! does OI update per tile
	! returns updated/Analysis states back to the caller, surface drive (sfcdrv) program
	!
	! RLA, RLO: lat lon information for the tile
	! SNOFCS(LENSFC), SWDFCS(LENSFC): snowdepth and snow water equivalent forecat (background)
	! (background/forecast in arrays of LENSFC)
	! compute snow depnsity (SNWDEN) and fraction of snow (if needed) from SWE/SNWD 
	!
	! IDIM * JDIM = LENSFC: number of grid cells in tile = xdim * ydim   
	! IY, IM, ID, IH = year, month, day, hour of current model step   
	! MYRANK: rank/id of the MPI process
	!
	! Outputs:
	! SNOANL, SWDANL:  SWE and snowdepthafter DA
	!                   
	!----------------------------------------------------------------------
	IMPLICIT NONE
	!
	include 'mpif.h'
	
	integer, parameter :: dp = kind(1.d0)

	INTEGER, intent(in) :: SNOW_IO_TYPE, MAX_TASKS, MYRANK, NPROCS, IDIM, JDIM, IY, IM, ID, IH, LENSFC, IVEGSRC	
	CHARACTER(LEN=*), Intent(In)   :: GHCND_SNOWDEPTH_PATH, IMS_SNOWCOVER_PATH, & 
									  IMS_INDEXES_PATH, CURRENT_ANALYSIS_PATH
!  GHCND_SNOWDEPTH_PATH="/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/GHCND/",
!  IMS_SNOWCOVER_PATH="/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS/",
!  IMS_INDEXES_PATH="/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS_INDEXES/",
!  CURRENT_ANALYSIS_PATH = "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/Cur_Analysis/"
	REAL, intent(Out)   :: SNOANL(LENSFC)   !, SWDANL(LENSFC), anl_fSCA(LENSFC)
	CHARACTER(LEN=5)    :: TILE_NUM
	Character(LEN=3)    :: rank_str
	INTEGER			    :: IERR	
	REAL        :: RLA(LENSFC), RLO(LENSFC), RLO_Tile(LENSFC), OROG(LENSFC), OROG_UF(LENSFC)
	REAL                :: SNOFCS_Inp(LENSFC), SWDFCS(LENSFC)
	REAL                :: SNOFCS(LENSFC), SNWDEN(LENSFC), VETFCS(LENSFC), SWDANL(LENSFC)	
	CHARACTER(len=250)   :: dim_name, ghcnd_inp_file, ims_inp_file, ims_inp_file_indices
	CHARACTER(len=250) 	 :: anl_inp_path, da_out_file
	CHARACTER(len=5)     :: y_str, m_str, d_Str, h_str, fvs_tile
	REAL, ALLOCATABLE    :: SNWD_GHCND(:), SNOFCS_atGHCND(:)		!, SWDFCS_atGHCND(:)
	REAL, ALLOCATABLE    :: Lat_GHCND(:), Lon_GHCND(:), Ele_GHCND(:)  
	REAL                 :: lat_min, lat_max, lon_min, lon_max  	
	Real		     :: SNCOV_IMS(LENSFC)  ! ims resampled at each grid
	Real		     :: SNWD_IMS_at_Grid(LENSFC)

	INTEGER :: Num_Ghcnd, Num_Ims, num_Eval, num_sub !Num_Ims_Lat, Num_Ims_Lon
	Real	:: max_distance   ! radius_of_influence for selecting state at observation point
	INTEGER :: jndx, zndx, ncol, nrow
	Integer, Allocatable   :: Loc_backSt_atObs(:)
	Integer				   :: num_loc, num_loc_1, max_num_loc
	Real, Parameter		:: Stdev_back = 30., Stdev_Obs_depth = 40., Stdev_Obs_ims = 80. ! mm 
	Integer				:: ims_assm_hour
	Real				:: obs_tolerance, ims_max_ele
	Real                :: L_horz, h_ver
	Real(dp), Allocatable 	 :: B_cov_mat(:,:), b_cov_vect(:)
	Real(dp), Allocatable 	 :: O_cov_mat(:,:), W_wght_vect(:)
	Real, Allocatable   :: back_at_Obs(:), obs_Array(:), Lat_Obs(:), Lon_Obs(:), Ele_Obs(:)
	REAL                :: innov_at_Grid(LENSFC), anl_at_Grid(LENSFC)
	Real, Allocatable  :: obs_Innov(:)
	Real               :: SNOANL_Cur(LENSFC), SWDANL_Cur(LENSFC)
	REAL, ALLOCATABLE  :: SNOANL_Cur_atEvalPts(:)  !evalution points 
	REAL, ALLOCATABLE  :: SNOFCS_atEvalPts(:), innov_atEvalPts(:), SNOANL_atEvalPts(:)  !evalution points 
	REAL, ALLOCATABLE  :: Lat_atEvalPts(:), Lon_atEvalPts(:), Obs_atEvalPts(:)     !evalution points
	Integer, ALLOCATABLE 	:: index_back_atEval(:)     ! background locations at eval points 
	Integer, ALLOCATABLE	:: index_back_atObs(:)   ! the location of background corresponding obs

	Real			   :: ims_threshold      ! threshold for converting IMS fSCA to binary 1, 0 values	
	Real			   :: snwden_val
	LOGICAL 		   :: assim_SWE, assim_GHCND   !assimilate swe? (instead of snow depth)
	LOGICAL 		   :: assim_IMS, assim_IMS_thisGridCell    !assimilate sncov

	Integer 		   :: veg_type_landice  ! 10.21.20: no assmn over land ice
    ! for mpi par
	INTEGER            :: Np_ext, Np_til, p_tN, p_tRank, N_sA, N_sA_Ext, mp_start, mp_end
	INTEGER            :: send_proc, rec_stat(MPI_STATUS_SIZE), dest_Aoffset, pindex
	INTEGER            :: mpiReal_size, rsize

	assim_SWE = .False.  ! note: if this is set true, may need to adjust background and obs stddev above
	If (SNOW_IO_TYPE == 2) assim_SWE = .True.
	assim_GHCND = .True.
	assim_IMS = .True.     !!assimilate sncov; 
	assim_IMS_thisGridCell = .FALSE.    ! if assimilating ims, skip this grid cell for this time step

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
	mp_end = (p_tRank + 1) * N_sA + N_sA_Ext 		! end index of subarray for proc	
	If(myrank == 0 )PRINT*,"sub array length ", N_sA, " extra sub array: ", N_sA_Ext
! if (p_tN /= 2 ) goto 999

	if( (GHCND_SNOWDEPTH_PATH(1:8).eq.'        ').and. & 
	    (IMS_SNOWCOVER_PATH(1:8).eq.'        ') ) then
		print*, "Observation paths don't exist!, skipping the OI DA"
		goto 999
	end if
   
! ToDO: use file location the way the rest of the program accesses
	! 4.9.20 for now hard code the file location 
	write(y_str, "(I4)") IY
	write(m_str, "(I0.2)") IM
	write(d_str, "(I0.2)") ID
	write(h_str, "(I0.2)") IH
	write(fvs_tile, "(I3)") IDIM	
    ! READ THE OROGRAPHY AND GRID POINT LAT/LONS FOR THE CUBED-SPHERE TILE.
	CALL READ_LAT_LON_OROG_atRank(p_tN, RLA,RLO,OROG,OROG_UF,TILE_NUM,IDIM,JDIM,LENSFC)
	PRINT*,"Snow anl on ", MYRANK, " Tile group: ", p_tN, " Tile: ", TILE_NUM
   ! READ THE INPUT SURFACE DATA ON THE CUBED-SPHERE TILE.
	Call READ_Forecast_Data(p_tN, LENSFC, SNOFCS_Inp, SWDFCS, VETFCS)  !VEGFCS, 	
! 10.22.20 Cur_Analysis: existing SNODEP based analysis
	!"/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/Cur_Analysis/"
	anl_inp_path = TRIM(CURRENT_ANALYSIS_PATH)//"snow."// &
				   TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"." // &
				   TRIM(h_str)// "0000.sfcanl_data."//TILE_NUM//".nc"
	Call READ_Analysis_Data(anl_inp_path, p_tN, LENSFC, SNOANL_Cur, SWDANL_Cur)
	
	RLO_Tile = RLO
	Where(RLO_Tile > 180) RLO_Tile = RLO_Tile - 360
	lat_min = MAX(MINVAL(RLA) - 1., -90.)
	lat_max = MIN(MAXVAL(RLA) + 1., 90.)
	lon_min = MAX(MINVAL(RLO_Tile) - 1., -180.)
	lon_max = MIN(MAXVAL(RLO_Tile) + 1., 180.)	
	! lon_min = MAX(MINVAL(RLO) - 1., 0.)
	! lon_max = MIN(MAXVAL(RLO) + 1., 360.)	
	if (p_tN==3)  then  
		lon_min = 125.       ! lon_min is left, lon_max is right, not necessary min/max value
		lon_max = -145.
	endif	
	if ((p_tRank==0) .and. print_deb) then
		print*, "Tile ", p_tN, " min/max lat/lon ", lat_min, lat_max, lon_min, lon_max
	endif
	! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/GHCND/
	ghcnd_inp_file = TRIM(GHCND_SNOWDEPTH_PATH)//"GHCND.SNWD."// &
					 TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18.nc"
	dim_name = "Site_Id"	
	call Observation_Read_GHCND_Tile_excNaN(p_tN, ghcnd_inp_file, dim_name, &
					lat_min, lat_max, lon_min, lon_max, & 
					Num_Ghcnd, SNWD_GHCND,              & !SNOObs_atGHCND,		&
					Lat_GHCND, Lon_GHCND, MYRANK)  !Ele_GHCND,		&
	! Do jndx=0, 5
	! 	if ((p_tN == jndx) .and. (p_tRank==0)) then
	! 		print*, "Tile ", p_tN, " num. GHCND obs ", Num_Ghcnd
	! 		PRINT*, "GHCND SNWD from rank: ", MYRANK
	! 		PRINT*, SNOObs_atGHCND
	! 	endif
	! End do
	! Stop
	! call Observation_Read_GHCND_Tile(ghcnd_inp_file, dim_name, &
	! 				lat_min, lat_max, lon_min, lon_max, & 
	! 				Num_Ghcnd, SNWD_GHCND,		&
	! 				Lat_GHCND, Lon_GHCND, MYRANK)  !Ele_GHCND,		&
	! call Observation_Read_GHCND(ghcnd_inp_file, dim_name, Num_Ghcnd, SNWD_GHCND,		&
	! 				Lat_GHCND, Lon_GHCND, MYRANK)  !Ele_GHCND,		&
	if ((p_tRank==0) .and. print_deb) then
		print*, "Tile ", p_tN, " num. GHCND obs ", Num_Ghcnd
	endif
	if ((p_tRank==0) .and. (p_tN==2) .and. print_deb) then
		PRINT*, "GHCND SNWD from rank: ", MYRANK
		PRINT*, SNWD_GHCND
		PRINT*, "Lat at GHCND from rank: ", MYRANK
		PRINT*, Lat_GHCND
		PRINT*, "Lon at GHCND from rank: ", MYRANK
		PRINT*, Lon_GHCND
		PRINT*,'Finished reading GHCND'
	endif
	if (myrank==0) PRINT*,'Finished reading GHCND'
   
	! (max) number of IMS subcells within a tile grid cell
	If (IDIM == 96) then          
		num_sub = 627               
		max_distance = 240.		!Km radius of influence: distance from gridcell to search for observations
	elseif (IDIM == 128) then
	   num_sub = 627               
	   max_distance = 240.		!Km 
	elseif (IDIM == 192) then
		!num_sub = 30 
		PRINT*,'Error, tile type not known '
		stop
	elseif (IDIM == 384) then
		!num_sub = 30  
		PRINT*,'Error, tile type not known '
		stop      
	elseif (IDIM == 768) then
		num_sub = 30
		max_distance = 27.        			!Km  
	else
		PRINT*,'Error, tile type not known '
		stop   
	endif
	  
	! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS/"
	ims_inp_file = TRIM(IMS_SNOWCOVER_PATH)//"IMS.SNCOV."// &
				   TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18.nc"                      !
	! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS_INDEXES/"
	ims_inp_file_indices = TRIM(IMS_INDEXES_PATH)//"C"//TRIM(ADJUSTL(fvs_tile))//&
						   ".IMS.Indices."//TRIM(TILE_NUM)//".nc"			
	Call Observation_Read_IMS_Full(ims_inp_file, ims_inp_file_indices, &
				                   MYRANK, JDIM, IDIM, num_sub, SNCOV_IMS)
	if((p_tRank==0) .and. (p_tN==2) .and. print_deb) then
		PRINT*, "IMS SNCOV from rank: ", MYRANK
		PRINT*, SNCOV_IMS
		PRINT*,'Finished reading IMS'
	endif 
	if (myrank==0) PRINT*,'Finished reading IMS'
	
	! 4.8.20: compute snow density ratio from forecast swe and snwdepth
	SNWDEN = SWDFCS / SNOFCS_Inp
	snwden_val = SUM(SNWDEN, Mask = (.not. IEEE_IS_NAN(SNWDEN))) &
	                 / COUNT (.not. IEEE_IS_NAN(SNWDEN))
	SNWDEN = snwden_val  !10.   ! ratio of sndepth to swe when swe /=0 : SWDFCS / SNOFCS
	if(p_tRank == 0) then
		print*, "Process ", MYRANK, " average snwd/swe ratio: ", snwden_val
	endif
	if (assim_SWE) then 
		SNWD_GHCND = SNWD_GHCND / SNWDEN		!snwden_val
		SNOFCS = SNOFCS_Inp
	else
		SNOFCS = SNOFCS_Inp * SNWDEN    !snwden_val  !SWDFCS			! assimilate depth
		! SNWD_GHCND = SNWD_GHCND 
	endif

	! Get model states at obs points
	ALLOCATE(SNOFCS_atGHCND(Num_Ghcnd))
	ALLOCATE(Ele_GHCND(Num_Ghcnd)) 
	ALLOCATE(index_back_atObs(Num_Ghcnd)) 
    num_Eval = floor(0.05 * Num_Ghcnd)      ! using 5% of ghcnd locations for evaluation
	ALLOCATE(index_back_atEval(num_Eval)) 
	ALLOCATE(Obs_atEvalPts(num_Eval)) 
	ALLOCATE(SNOFCS_atEvalPts(num_Eval)) 
	ALLOCATE(Lat_atEvalPts(num_Eval))
	ALLOCATE(Lon_atEvalPts(num_Eval)) 
	ALLOCATE(innov_atEvalPts(num_Eval))
	ALLOCATE(SNOANL_atEvalPts(num_Eval)) 	
	ALLOCATE(SNOANL_Cur_atEvalPts(num_Eval)) 
	! ALLOCATE(SNOANL_atEvalPts(Num_Ghcnd)) 	
	! ALLOCATE(SNOANL_Cur_atEvalPts(Num_Ghcnd)) 		
	if(p_tRank == 0) then 
		PRINT*, "Tile ", p_tN+1, " ", num_Eval, ' points for evaluation excluded from DA'	
	endif	

	Call Observation_Operator_Parallel(Myrank, MAX_TASKS, p_tN, p_tRank, Np_til, & 
	                    RLA, RLO, OROG, Lat_GHCND, Lon_GHCND,   &
						LENSFC, Num_Ghcnd, num_Eval, max_distance, SNOFCS, SNWD_GHCND,  &
						SNOFCS_atGHCND, Ele_GHCND, index_back_atObs, index_back_atEval, &
						Obs_atEvalPts, SNOFCS_atEvalPts, Lat_atEvalPts, Lon_atEvalPts)
	if ((p_tN==2) .and. (p_tRank==0) .and. print_deb) then 
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
	!Stop
	! call Observation_Operator(RLA, RLO, OROG, Lat_GHCND, Lon_GHCND,   &
	! 					LENSFC, Num_Ghcnd, max_distance,            &
	! 					SNOFCS,                                  &
	! 					SNOFCS_atGHCND, Ele_GHCND, index_back_atGHCND)
	if ((p_tN==2) .and. (p_tRank==0) .and. print_deb) then
		PRINT*, "GHCND Lat range from rank: ", MYRANK, MINVAL(Lat_GHCND), " ", MAXVAL(Lat_GHCND)
		PRINT*, "GHCND Lon range from rank: ", MYRANK, MINVAL(Lon_GHCND), " ", MAXVAL(Lon_GHCND)
		PRINT*, "Background SWE at GHCND locations from rank: ", MYRANK
		PRINT*, SNOFCS_atGHCND	
		PRINT*, "Elevation at GHCND locations from rank: ", MYRANK
		PRINT*, Ele_GHCND
		PRINT*,'Finished observation operator for GHCND'	
	endif
	if (myrank==0) PRINT*,'Finished observation operator for GHCND'		
	ims_threshold = 0.5  ! threshold for converting IMS fSCA to binary 1, 0 values
	! Call Observation_Operator_IMS_fSCA_Threshold(SNCOV_IMS, SNOFCS, SNWDEN, assim_SWE,	   &
	! 						LENSFC, ims_threshold, SNWD_IMS_at_Grid) 
	call Observation_Operator_IMS_fSCA(SNCOV_IMS, SNWDEN, VETFCS, assim_SWE, LENSFC, SNWD_IMS_at_Grid) 
	if ((p_tN==2) .and. (p_tRank==0) .and. print_deb) then
		PRINT*, "IMS obs at each grid cell from rank: ", MYRANK
		PRINT*, SNWD_IMS_at_Grid
		PRINT*,'Finished observation operator for IMS'
	endif
	if (myrank==0) PRINT*,'Finished observation operator for IMS'
	
	L_horz = 55.  !120.  !
	h_ver = 800.  !1200.	!
	obs_tolerance = 5.0
	max_num_loc = 50   !100	!
	ims_max_ele = 1500.
	ims_assm_hour = 18
	max_distance = 250.   !120.			!Km
	!if (myrank==4) then 
	if (myrank==0) PRINT*,'Starting DA loop'

	! 10.21.20: no assimilation if land-ice
	if (IVEGSRC == 2) then   ! sib
		veg_type_landice=13
	else
		veg_type_landice=15
	endif
	!Do ncol=1, IDIM
	!Do nrow = 1, Num_Snotel !JDIM/2
	    Do jndx = mp_start, mp_end     !1, LENSFC		!jndx = (ncol-1)*JDIM + nrow		! 
			num_loc_1 = 0
			assim_IMS_thisGridCell = .FALSE.
			call debug_print("loop ", float(jndx))
			! GHCND			 
			if(assim_GHCND) then
				call nearest_Observations_Locations(RLA(jndx), RLO(jndx),    &
						Lat_GHCND, Lon_GHCND,  Num_Ghcnd, max_distance, max_num_loc,   &
						Stdev_back, Stdev_Obs_depth, obs_tolerance,                 &
						SNOFCS_atGHCND, SNWD_GHCND,						 &
						Loc_backSt_atObs,  num_loc_1) !,      &LENSFC,
				call debug_print("number of GHCND sndpth obs ", float(num_loc))
			endif
			num_loc = num_loc_1
			!check IMS is assimilated
			if (assim_IMS) then
				if((.NOT. IEEE_IS_NAN(SNWD_IMS_at_Grid(jndx))) .AND. &
				   (OROG(jndx) <= ims_max_ele) .AND. &
				   (IH == ims_assm_hour)) then
					num_loc = num_loc + 1
					assim_IMS_thisGridCell = .TRUE.
				endif
			endif
			! if assim_IMS=false >> num_loc_1=num_loc
			! 10.21.20: no assimilation if land-ice	.NE. veg_type_landice	
			if((num_loc > 0) .and. (NINT(VETFCS(jndx)) /= veg_type_landice)) then ! .and. (SNCOV_IMS(jndx) > 0.)) then    
				! get background states
				Allocate(back_at_Obs(num_loc))
				Allocate(obs_Array(num_loc))
				Allocate(Lat_Obs(num_loc))
				Allocate(Lon_Obs(num_loc))
				Allocate(Ele_Obs(num_loc))
				! ghcnd
				if(num_loc_1 > 0) then
					Do zndx = 1, num_loc_1     
						back_at_Obs(zndx) = SNOFCS_atGHCND(Loc_backSt_atObs(zndx))
						obs_Array(zndx) = SNWD_GHCND(Loc_backSt_atObs(zndx))
						Lat_Obs(zndx) = Lat_GHCND(Loc_backSt_atObs(zndx))
						Lon_Obs(zndx) = Lon_GHCND(Loc_backSt_atObs(zndx))
						Ele_Obs(zndx) = Ele_GHCND(Loc_backSt_atObs(zndx))
					End Do
				End if
				!ims
				if(assim_IMS_thisGridCell) then
					back_at_Obs(num_loc) = SNOFCS(jndx)
					obs_Array(num_loc) = SNWD_IMS_at_Grid(jndx)
					Lat_Obs(num_loc) = RLA(jndx)   !Lat_IMS_atGrid(jndx)
					Lon_Obs(num_loc) = RLO(jndx)   !Lon_IMS_atGrid(jndx)
					Ele_Obs(num_loc) = OROG(jndx)  !Ele_IMS(jndx)
				endif
				! compute covariances
				Allocate(B_cov_mat(num_loc, num_loc))
				Allocate(b_cov_vect(num_loc))
				Allocate(O_cov_mat(num_loc, num_loc))
				Allocate(W_wght_vect(num_loc))
				call compute_covariances(RLA(jndx), RLO(jndx), OROG(jndx), SNOFCS(jndx),    &
					Lat_Obs, Lon_Obs, Ele_Obs, num_loc,    			 &
					Stdev_back, Stdev_Obs_depth, Stdev_Obs_ims,      &
					L_horz, h_ver,                                   &   !L_horz in Km, h_ver in m
					assim_IMS_thisGridCell,                          &
					B_cov_mat, b_cov_vect, O_cov_mat, W_wght_vect)
				! call OI DA
				Allocate(obs_Innov(num_loc))
				call Snow_DA_OI(back_at_Obs, obs_Array, num_loc, W_wght_vect,            &
					SNOFCS(jndx), innov_at_Grid(jndx), anl_at_Grid(jndx), obs_Innov)
				if ((p_tN==2) .and. (p_tRank==0) .and. print_deb) then	
					print*, "proc ", myrank, "loop ", jndx, "num depth obs ", num_loc_1, "total obs", num_loc
					PRINT*, " background at obs pts: "
					PRINT*, back_at_Obs	
					PRINT*, "Observed"
					PRINT*,  obs_Array
					PRINT*, "Obs innovation: "
					PRINT*, obs_Innov
					PRINT*, "Weight vector: "
					PRINT*, W_wght_vect	
					print*, "innov: ", innov_at_Grid(jndx), "forec: ", SNOFCS(jndx), " anl: ", anl_at_Grid(jndx)
				endif		
				!free mem
				DEALLOCATE(back_at_Obs, obs_Array)
				DEALLOCATE(Lat_Obs, Lon_Obs, Ele_Obs, obs_Innov)
				DEALLOCATE(B_cov_mat, b_cov_vect, O_cov_mat, W_wght_vect)
				! QCC by ims--use ims as snow mask
				if((SNCOV_IMS(jndx) >= 0.5) .and. (anl_at_Grid(jndx) < 50.)) then
					anl_at_Grid(jndx) = 50.
				! elseif((SNCOV_IMS(jndx) < 0.5) .and. (SNCOV_IMS(jndx) > 0.1) .and. (anl_at_Grid(jndx) >= 50)) then
				! 	anl_at_Grid(jndx) = 50.
				elseif((SNCOV_IMS(jndx) < 0.1) .and. (anl_at_Grid(jndx) >= 0.)) then
						anl_at_Grid(jndx) = 0.
				endif
				if ((p_tN==2) .and. (p_tRank==0) .and. print_deb) then						
					PRINT*, "analyis at grid after mask: ", anl_at_Grid(jndx)
				endif
			else
				anl_at_Grid(jndx) = SNOFCS(jndx)
			endif
			if (assim_GHCND) Deallocate(Loc_backSt_atObs) 
		End do
	!End do
	if (myrank==0) PRINT*, 'Finished DA loops'

! ToDO: Better way to handle this?
! Real data type size corresponding to mpi
	rsize = SIZEOF(snwden_val) 
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
		call MPI_SEND(anl_at_Grid(mp_start:mp_end), N_sA, mpiReal_size, p_tN,   &
					  MYRANK, MPI_COMM_WORLD, IERR) 
		call MPI_SEND(innov_at_Grid(mp_start:mp_end), N_sA, mpiReal_size, p_tN,   &
					  MYRANK*100, MPI_COMM_WORLD, IERR)
	else    !if(p_tRank == 0) then  
		Do pindex =  1, (Np_til - 1)   ! sender proc index within tile group
			dest_Aoffset = pindex * N_sA + N_sA_Ext + 1   ! dest array offset
			send_proc = MYRANK +  pindex * MAX_TASKS
			call MPI_RECV(anl_at_Grid(dest_Aoffset:dest_Aoffset+N_sA-1), N_sA, mpiReal_size, send_proc,      &
					  send_proc, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
			call MPI_RECV(innov_at_Grid(dest_Aoffset:dest_Aoffset+N_sA-1), N_sA, mpiReal_size, send_proc,      &
					  send_proc*100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
		enddo
	endif
	if (myrank==0) PRINT*,'Finished Data copy'

	if (MYRANK > MAX_TASKS - 1 ) goto 998   ! if(p_tRank /= 0 ) goto 998

	! avoid -ve anl
	Where(anl_at_Grid < 0.) anl_at_Grid = 0.
	! swe and snwd
	if(assim_SWE) then
		SNOANL = anl_at_Grid 
		SWDANL = anl_at_Grid * SNWDEN		!snwden_val    
		!SNOANL_Cur = SNOANL_Cur 
	else
		SNOANL = anl_at_Grid / SNWDEN		!snwden_val 
		SWDANL = anl_at_Grid
		SNOANL_Cur = SWDANL_Cur
		!innov_at_Grid = innov_at_Grid / SNWDEN
	endif	
	if (print_deb) then
		PRINT*, "Innovation SWE/snwd from rank: ", MYRANK
	    PRINT*, innov_at_Grid	
	    PRINT*, "Analysis SWE/ snwd  from rank: ", MYRANK
	    PRINT*, anl_at_Grid	
	endif

	!Compute updated snocov	
	!Call update_snow_cover_fraction(LENSFC, SNOANL, VETFCS, anl_fSCA)

	! copy values at eval points
	innov_atEvalPts = IEEE_VALUE(innov_atEvalPts, IEEE_QUIET_NAN)
	SNOANL_atEvalPts = IEEE_VALUE(SNOANL_atEvalPts, IEEE_QUIET_NAN)
	SNOANL_Cur_atEvalPts = IEEE_VALUE(SNOANL_Cur_atEvalPts, IEEE_QUIET_NAN)
	Do jndx = 1, num_Eval
		if ((index_back_atEval(jndx) > 0) .and. &
		    (NINT(VETFCS(jndx)) /= veg_type_landice)) then  !10.21.20: exclude if land-ice	
				innov_atEvalPts(jndx) = innov_at_Grid(index_back_atEval(jndx))
				SNOANL_atEvalPts(jndx) = anl_at_Grid(index_back_atEval(jndx))
				SNOANL_Cur_atEvalPts(jndx) = SNOANL_Cur(index_back_atEval(jndx))
		endif
	End do
	! Do jndx = 1, Num_Ghcnd
	! 	if ((index_back_atObs(jndx) > 0) .and. &
	! 	    (NINT(VETFCS(jndx)) /= veg_type_landice)) then  !10.21.20: exclude if land-ice	
	! 			SNOANL_atEvalPts(jndx) = anl_at_Grid(index_back_atObs(jndx))
	! 			SNOANL_Cur_atEvalPts(jndx) = SNOANL_Cur(index_back_atObs(jndx))
	! 	endif
	! End do

	! write outputs	
	Write(rank_str, '(I3.3)') (MYRANK+1)
	! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/Analysis/"
	da_out_file = "./SNOANLOI."// &  !
				  TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18_tile"//rank_str//".nc"  !	
	call Write_DA_Outputs(da_out_file, IDIM, JDIM, LENSFC, MYRANK, &
	                      SNOFCS_Inp, SNOANL, SWDFCS, SWDANL, &  !
			            innov_at_Grid, SNCOV_IMS, &
			    !       Num_Ghcnd, Lat_GHCND, Lon_GHCND, SNWD_GHCND, &
			    ! SNOFCS_atGHCND, SNOFCS_atGHCND, SNOANL_atEvalPts, SNOANL_Cur_atEvalPts)
			    num_Eval, Lat_atEvalPts, Lon_atEvalPts, Obs_atEvalPts, & 
			SNOFCS_atEvalPts, innov_atEvalPts, SNOANL_atEvalPts, SNOANL_Cur_atEvalPts)  !, anl_fSCA) !updated snocov

998 CONTINUE
	DEALLOCATE(SNWD_GHCND, SNOFCS_atGHCND, Lat_GHCND, Lon_GHCND, Ele_GHCND)  !, index_back_atGHCND)
	DEALLOCATE(Obs_atEvalPts, SNOFCS_atEvalPts, innov_atEvalPts, SNOANL_atEvalPts)
	DEALLOCATE(SNOANL_Cur_atEvalPts) 
	DEALLOCATE(index_back_atObs, index_back_atEval, Lat_atEvalPts, Lon_atEvalPts) !, Ele_atEvalPts)
	!DEALLOCATE(SNCOV_IMS, Lat_IMS, Lon_IMS, SNOFCS_atIMS, SWDFCS_atIMS)
999 CONTINUE
    PRINT*,'Finished OI DA ON RANK: ', MYRANK
	CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

	!STOP

	RETURN

 END subroutine Snow_Analysis_OI

 subroutine map_outputs_toObs(MAX_TASKS, MYRANK, NPROCS, IDIM, JDIM, IY, IM, ID, IH, & 
	                          LENSFC, CURRENT_ANALYSIS_PATH)
							
	IMPLICIT NONE
	!
	include 'mpif.h'
	
	integer, parameter :: dp = kind(1.d0)

	INTEGER, intent(in) :: MAX_TASKS, MYRANK, NPROCS, IDIM, JDIM, IY, IM, ID, IH, LENSFC
	CHARACTER(LEN=*), Intent(In)   :: CURRENT_ANALYSIS_PATH
	CHARACTER(LEN=5)    :: TILE_NUM
	Character(LEN=3)    :: rank_str
	INTEGER			    :: IERR	
	REAL        :: RLA(LENSFC), RLO(LENSFC), OROG(LENSFC), OROG_UF(LENSFC)
	REAL                 :: SNOANL(LENSFC), SWDANL(LENSFC)	
	CHARACTER(len=250)   :: dim_name, anl_inp_path, da_out_file
	CHARACTER(len=5)     :: y_str, m_str, d_Str, h_str, fvs_tile
	REAL, ALLOCATABLE    :: Lat_GHCND(:), Lon_GHCND(:), Lon_Obs(:)
	INTEGER              :: num_Eval, jndx
	Integer, Allocatable   :: index_back_atEval(:)
	REAL, ALLOCATABLE  :: SNOANL_atEvalPts(:)  !evalution points	
    ! for mpi par
	INTEGER            :: Np_ext, Np_til, p_tN, p_tRank, N_sA, N_sA_Ext, mp_start, mp_end

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
	mp_end = (p_tRank + 1) * N_sA + N_sA_Ext 		! end index of subarray for proc	
	If (myrank == 0 ) PRINT*,"sub array length ", N_sA, " extra sub array: ", N_sA_Ext
	! if (p_tN /= 2 ) goto 999
	
    ! READ THE OROGRAPHY AND GRID POINT LAT/LONS FOR THE CUBED-SPHERE TILE.
	CALL READ_LAT_LON_OROG_atRank(p_tN, RLA,RLO,OROG,OROG_UF,TILE_NUM,IDIM,JDIM,LENSFC) 	
	PRINT*,"Snow anl on ", MYRANK, " Tile group: ", p_tN, " Tile: ", TILE_NUM
 
	write(y_str, "(I4)") IY
	write(m_str, "(I0.2)") IM
	write(d_str, "(I0.2)") ID
	write(h_str, "(I0.2)") IH
	write(fvs_tile, "(I3)") IDIM
	!"/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/Cur_Analysis/"
	anl_inp_path = TRIM(CURRENT_ANALYSIS_PATH)//"snow."// &
				TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"." // &
				TRIM(h_str)// "0000.sfcanl_data."//TILE_NUM//".nc"
	Call READ_Analysis_Data(anl_inp_path, p_tN, LENSFC, SNOANL, SWDANL)
	
	Write(rank_str, '(I3.3)') (p_tN+1)	
	! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/Analysis2/
	da_out_file = "./SNOANL."// &  !
				  TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18_tile"//rank_str//".nc"						
	Call Observation_Read_atEval(da_out_file, num_Eval, Lat_GHCND, Lon_GHCND, p_tN)					
	if ((p_tRank==0) ) then
		print*, "Tile ", p_tN, " num. Eval ", num_Eval
		! PRINT*, "Lat at Eval from rank: ", MYRANK
		! PRINT*, Lat_GHCND
		! PRINT*, "Lon at Eval from rank: ", MYRANK
		! PRINT*, Lon_GHCND
	endif
	ALLOCATE(index_back_atEval(num_Eval)) 
	ALLOCATE(SNOANL_atEvalPts(num_Eval)) 
	ALLOCATE(Lon_Obs(num_Eval)) 

    Lon_Obs = Lon_GHCND
    Where(Lon_Obs < 0) Lon_Obs = 360. + Lon_Obs
	! grid indices from coarse (eg C384) to high res (C768)
	Call Find_Nearest_GridIndices_Parallel(MYRANK, MAX_TASKS, p_tN, p_tRank, Np_til, &
					LENSFC, num_Eval, RLA, RLO, Lat_GHCND, Lon_Obs, index_back_atEval)	
	if ( (p_tRank==0) ) then 
		PRINT*, "Background Indices at eval points rank ",MYRANK
		PRINT*, index_back_atEval	
	endif

	if (MYRANK > MAX_TASKS - 1 ) goto 998   ! if(p_tRank /= 0 ) goto 998

	SNOANL_atEvalPts = IEEE_VALUE(SNOANL_atEvalPts, IEEE_QUIET_NAN)
	Do jndx = 1, num_Eval
		if (index_back_atEval(jndx) > 0) then
			SNOANL_atEvalPts(jndx) = SWDANL(index_back_atEval(jndx))
		endif
	End do
	dim_name = 'SNOANL_Cur_atEvalPts'
	call Append_DA_Outputs(da_out_file, MYRANK, num_Eval, dim_name, SNOANL_atEvalPts)  

998 CONTINUE
	DEALLOCATE(Lat_GHCND, Lon_GHCND, Lon_Obs, SNOANL_atEvalPts, index_back_atEval) 
999 CONTINUE
    PRINT*,'Finished ON RANK: ', MYRANK
	CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

	STOP

	RETURN

End subroutine map_outputs_toObs

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
	Integer 		   :: indx, vtype_int

	!This is for the IGBP veg classification scheme.
	snupx = (/0.080, 0.080, 0.080, 0.080, 0.080, 0.020, 	&
			0.020, 0.060, 0.040, 0.020, 0.010, 0.020,			&
			0.020, 0.020, 0.013, 0.013, 0.010, 0.020,			&
			0.020, 0.020, 0.000, 0.000, 0.000, 0.000,			&
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
                                 snoforc, snoanl, snwdforc, snwdanal, inovatgrid, SNCOV_IMS, &
                                 num_Eval, Lat_atEvalPts, Lon_atEvalPts, Obs_atEvalPts, & 
            SNOFCS_atEvalPts, innov_atEvalPts, SNOANL_atEvalPts, SNOANL_Cur_atEvalPts)  !, anl_fSCA) !updated snocov
	!------------------------------------------------------------------
	! Write DA ouputs: 
	! forecast SWE
	! analysis SWE
	! analysis Snow Depth
	! innovation at grid
	!------------------------------------------------------------------
	implicit none

	CHARACTER(LEN=*), Intent(In)      :: output_file
	integer, intent(in)         :: idim, jdim, lensfc
	real, intent(in)            :: snoforc(lensfc), snoanl(lensfc), snwdforc(lensfc), snwdanal(lensfc)
	Real, intent(in)            :: inovatgrid(lensfc), SNCOV_IMS(lensfc)  !, anl_fSCA(lensfc)
	integer, intent(in)         :: num_Eval
	real, intent(in)    :: Lat_atEvalPts(num_Eval), Lon_atEvalPts(num_Eval), Obs_atEvalPts(num_Eval)
	real, intent(in)    :: SNOFCS_atEvalPts(num_Eval), SNOANL_atEvalPts(num_Eval) 
	real, intent(in)    :: innov_atEvalPts(num_Eval), SNOANL_Cur_atEvalPts(num_Eval)

	integer                     :: fsize=65536, inital=0
	integer                     :: header_buffer_val = 16384
	integer                     :: dims_3d(3), dims_strt(3), dims_end(3)
	integer                     :: error, i, ncid
	integer                     :: dim_x, dim_y, dim_time, dim_eval
	integer                     :: id_x, id_y, id_time
	integer                     :: id_swe_forc, id_swe, id_snwdf, id_snwd, id_innov, id_imscov   !, id_anlscov
	integer       :: id_lateval, id_loneval, id_obseval, id_forceval
	integer       :: id_anleval, id_cur_anleval, id_innoveval   !, id_anlscov

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
	error = nf90_def_dim(ncid, 'eval_points', num_Eval, dim_eval)
	call netcdf_err(error, 'DEFINING eval_points' )

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

	error = nf90_def_var(ncid, 'SNWD_Forecast', NF90_DOUBLE, dims_3d, id_snwdf)
	call netcdf_err(error, 'DEFINING SNWDPH Forecast' )
	error = nf90_put_att(ncid, id_snwdf, "long_name", "Forecast Snow Depth")
	call netcdf_err(error, 'DEFINING SNWDPH Forecast LONG NAME' )
	error = nf90_put_att(ncid, id_snwdf, "units", "mm")
	call netcdf_err(error, 'DEFINING SNWDPH Forecast UNITS' )

	error = nf90_def_var(ncid, 'SNWD_Analysis', NF90_DOUBLE, dims_3d, id_snwd)
	call netcdf_err(error, 'DEFINING SNWDPH Analyis' )
	error = nf90_put_att(ncid, id_snwd, "long_name", "Analysis Snow Depth")
	call netcdf_err(error, 'DEFINING SNWDPH Analysis LONG NAME' )
	error = nf90_put_att(ncid, id_snwd, "units", "mm")
	call netcdf_err(error, 'DEFINING SNWDPH Analysis UNITS' )

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

	! eval points
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

	error = nf90_def_var(ncid, 'SNOANL_Cur_atEvalPts', NF90_DOUBLE, dim_eval, id_cur_anleval)
	call netcdf_err(error, 'DEFINING SNOANL_Cur_atEvalPts' )
	error = nf90_put_att(ncid, id_cur_anleval, "long_name", "Current Analysis at Evaluation Points")
	call netcdf_err(error, 'DEFINING SNOANL_Cur_atEvalPts LONG NAME' )
	error = nf90_put_att(ncid, id_cur_anleval, "units", "mm")
	call netcdf_err(error, 'DEFINING SNOANL_Cur_atEvalPts UNITS' )
	
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
	call netcdf_err(error, 'WRITING SWE Ananlysis RECORD' ) 

	dum2d = reshape(snwdforc, (/idim,jdim/))
	error = nf90_put_var( ncid, id_snwdf, dum2d, dims_strt, dims_end)
	call netcdf_err(error, 'WRITING SNWDPH Forecast RECORD' )

	dum2d = reshape(snwdanal, (/idim,jdim/))
	error = nf90_put_var( ncid, id_snwd, dum2d, dims_strt, dims_end)
	call netcdf_err(error, 'WRITING SNWDPH analysis RECORD' )

	dum2d = reshape(inovatgrid, (/idim,jdim/))
	error = nf90_put_var( ncid, id_innov, dum2d, dims_strt, dims_end)
	call netcdf_err(error, 'WRITING innovation RECORD' )

	dum2d = reshape(SNCOV_IMS, (/idim,jdim/))
	error = nf90_put_var( ncid, id_imscov, dum2d, dims_strt, dims_end)
	call netcdf_err(error, 'WRITING imsfSCA RECORD' )

	! dum2d = reshape(anl_fSCA, (/idim,jdim/))
	! error = nf90_put_var( ncid, id_anlscov, dum2d, dims_strt, dims_end)
	! call netcdf_err(error, 'WRITING anlfSCA RECORD' )
	
	! eval points
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

	error = nf90_put_var( ncid, id_cur_anleval, SNOANL_Cur_atEvalPts)
	call netcdf_err(error, 'WRITING SNOANL_Cur_atEvalPts RECORD' )

	error = nf90_put_var( ncid, id_innoveval, innov_atEvalPts)
	call netcdf_err(error, 'WRITING innov_atEvalPts RECORD' )

	deallocate(x_data, y_data)
	deallocate(dum2d)

	error = nf90_close(ncid)
    
 End subroutine Write_DA_Outputs


 Subroutine Append_DA_Outputs(output_file, myrank,   &
                                 num_Eval, var_name, SNOANL_atEvalPts)  
	!------------------------------------------------------------------
	! Append var_name to existing nc file: 
	!------------------------------------------------------------------
	implicit none

	CHARACTER(LEN=*), Intent(In)      :: output_file, var_name
	integer, intent(in)         :: num_Eval
	real, intent(in)            :: SNOANL_atEvalPts(num_Eval)

	integer                     :: fsize=65536, inital=0
	integer                     :: header_buffer_val = 16384
	integer                     :: error, i, ncid, dim_eval, id_anleval
	integer                     :: myrank

	include "mpif.h"

	print*, "Process ", myrank, "writing output data to: ",trim(output_file)

	! error = NF90_CREATE(output_file, IOR(NF90_NETCDF4,NF90_CLASSIC_MODEL), ncid, initialsize=inital, chunksize=fsize)
	error = NF90_OPEN(trim(output_file), NF90_WRITE, NCID)
	call netcdf_err(error, 'Opening FILE='//trim(output_file) )
	
	! obs points
	! error = nf90_def_dim(ncid, 'eval_points', num_Eval, dim_eval)
    error = NF90_INQ_DIMID(ncid, 'eval_points', dim_eval)
    call netcdf_err(error, 'ERROR READING Dimension' )
	
	error=NF90_REDEF(ncid) 
	call netcdf_err(error, 'ERROR Redefining header' )

	error = nf90_def_var(ncid, trim(var_name), NF90_DOUBLE, dim_eval, id_anleval)
	call netcdf_err(error, 'DEFINING '//trim(var_name) )
	error = nf90_put_att(ncid, id_anleval, "long_name", trim(var_name))
	call netcdf_err(error, 'DEFINING '//trim(var_name)//' LONG NAME' )
	error = nf90_put_att(ncid, id_anleval, "units", "mm")
	call netcdf_err(error, 'DEFINING '//trim(var_name)//' UNITS' )

	error = nf90_enddef(ncid, header_buffer_val,4,0,4)
	call netcdf_err(error, 'DEFINING HEADER' )

	
	! eval points
	error = nf90_put_var( ncid, id_anleval, SNOANL_atEvalPts)
	call netcdf_err(error, 'WRITING SNOANL_Cur_atEvalPts RECORD' )

	error = nf90_close(ncid)
    
 End subroutine Append_DA_Outputs

 SUBROUTINE Observation_Read_atEval(ghcnd_inp_file, &   !dim_name,			&
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
    
        ! REAL, ALLOCATABLE, Intent(Out)    :: SNWD_GHCND(:)
        REAL, ALLOCATABLE, Intent(Out)	   :: Lat_GHCND(:), Lon_GHCND(:) !, Ele_GHCND(:)
    
        ERROR=NF90_OPEN(TRIM(ghcnd_inp_file),NF90_NOWRITE,NCID)
        CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(ghcnd_inp_file) )
    
        ERROR=NF90_INQ_DIMID(NCID, 'eval_points', ID_DIM)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Dimension' )
    
        ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NDIM)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Dimension' )
    
        ! ALLOCATE(SNWD_GHCND(NDIM))
        ALLOCATE(Lat_GHCND(NDIM))
        ALLOCATE(Lon_GHCND(NDIM))
        !ALLOCATE(Ele_GHCND(NDIM))
    
        ! ERROR=NF90_INQ_VARID(NCID, 'SNWD', ID_VAR)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING SNWD ID' )
        ! ERROR=NF90_GET_VAR(NCID, ID_VAR, SNWD_GHCND)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING SNWD RECORD' )
    
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
	Real, Intent(In) 	        :: RLA_cg(num_src), RLO_cg(num_src), RLA(num_tar), RLO(num_tar)
	Integer, Intent(Out)         :: index_ens_atGrid(num_tar)
	Real 	                    :: RLA_cg_rad(num_src), RLO_cg_rad(num_src)
	Real 	                    :: RLA_rad(num_tar), RLO_rad(num_tar)
	
	INTEGER                     :: indx, min_indx
	Real                        :: distArr(num_src), haversinArr(num_src)
	Real 	                    :: d_latArr(num_src), d_lonArr(num_src)
	Real(16), Parameter         :: PI_16 = 4 * atan (1.0_16)	
	Real(16), Parameter         :: pi_div_180 = PI_16/180.0
	Real, Parameter		        :: earth_rad = 6371.
	
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
	mp_end = (p_tRank + 1) * N_sA + N_sA_Ext 		! end index of subarray for proc
		
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
	LOGICAL, Intent(In)		          :: assim_SWE, save_Ens
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
		! 	Do jy = 1, JDIM_In
		! 		DUMMY2D(ix, jy) = DUMMY((ix+1)/2, (jy+1)/2)
		! 	End do
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
		! 	Do jy = 1, JDIM_In
		! 		DUMMY2D(ix, jy) = DUMMY((ix+1)/2, (jy+1)/2)
		! 	End do
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

		error = nf90_def_var(ncid, 'SNWD_Forecast', NF90_DOUBLE, dims_3d, id_forc)
		call netcdf_err(error, 'DEFINING SNWD_Forecast' )
		error = nf90_put_att(ncid, id_forc, "long_name", "Forecast Snow Depth")
		call netcdf_err(error, 'DEFINING SNWD Forecast LONG NAME' )
		error = nf90_put_att(ncid, id_forc, "units", "mm")
		call netcdf_err(error, 'DEFINING SNWD Forecast UNITS' )

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
    
 SUBROUTINE READ_Forecast_Data(MYRANK, LENSFC, SNOFCS, SWDFCS, VETFCS) !VEGFCS, !SRFLAG)
    
	IMPLICIT NONE

	include "mpif.h"

	INTEGER, INTENT(IN)       :: MYRANK, LENSFC
	REAL, INTENT(OUT)         :: SNOFCS(LENSFC), SWDFCS(LENSFC), VETFCS(LENSFC)  !VEGFCS(LENSFC), 
	!REAL, INTENT(OUT)        :: FMM(LENSFC), FHH(LENSFC), SRFLAG(LENSFC)

	CHARACTER(LEN=50)         :: FNBGSI
	CHARACTER(LEN=3)          :: RANKCH

	INTEGER                   :: ERROR, NCID
	INTEGER                   :: IDIM, JDIM, ID_DIM
	INTEGER                   :: ID_VAR

	REAL(KIND=8), ALLOCATABLE :: DUMMY(:,:)

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
	SNOFCS = RESHAPE(DUMMY, (/LENSFC/))

	ERROR=NF90_INQ_VARID(NCID, "snwdph", ID_VAR)
	CALL NETCDF_ERR(ERROR, 'READING snwdph ID' )
	ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
	CALL NETCDF_ERR(ERROR, 'READING snwdph' )
	SWDFCS = RESHAPE(DUMMY, (/LENSFC/))
	!  ERROR=NF90_INQ_VARID(NCID, "vfrac", ID_VAR)
	!  CALL NETCDF_ERR(ERROR, 'READING vfrac ID' )
	!  ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
	!  CALL NETCDF_ERR(ERROR, 'READING vfrac' )
	!  VEGFCS = RESHAPE(DUMMY, (/LENSFC/))  
	ERROR=NF90_INQ_VARID(NCID, "vtype", ID_VAR)
	CALL NETCDF_ERR(ERROR, 'READING vtype ID' )
	ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
	CALL NETCDF_ERR(ERROR, 'READING vtype' )
	VETFCS = RESHAPE(DUMMY, (/LENSFC/))    
	!  ERROR=NF90_INQ_VARID(NCID, "srflag", ID_VAR)
	!  CALL NETCDF_ERR(ERROR, 'READING srflag ID' )
	!  ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
	!  CALL NETCDF_ERR(ERROR, 'READING srflag' )
	!  SRFLAG = RESHAPE(DUMMY, (/LENSFC/))  
	
	DEALLOCATE(DUMMY)

	ERROR = NF90_CLOSE(NCID)
    
 END SUBROUTINE READ_Forecast_Data

 SUBROUTINE READ_Forecast_Data_atPath(forc_inp_path, MYRANK, LENSFC, SNOFCS, SWDFCS) !, VETFCS) !VEGFCS, !SRFLAG)
    
	IMPLICIT NONE

	include "mpif.h"
	
	CHARACTER(LEN=*), Intent(In)      :: forc_inp_path
	INTEGER, INTENT(IN)       :: MYRANK, LENSFC
	REAL, INTENT(OUT)         :: SNOFCS(LENSFC), SWDFCS(LENSFC) !, VETFCS(LENSFC)  !VEGFCS(LENSFC), 
	!REAL, INTENT(OUT)        :: FMM(LENSFC), FHH(LENSFC), SRFLAG(LENSFC)

	CHARACTER(LEN=50)         :: FNBGSI
	CHARACTER(LEN=3)          :: RANKCH

	INTEGER                   :: ERROR, NCID
	INTEGER                   :: IDIM, JDIM, ID_DIM
	INTEGER                   :: ID_VAR

	REAL(KIND=8), ALLOCATABLE :: DUMMY(:,:)

	!CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, ERROR)
	! WRITE(RANKCH, '(I3.3)') (MYRANK+1)
	! FNBGSI = "./fnbgsi." // RANKCH
	! if (print_deb) PRINT*, "READ INPUT SFC DATA FROM: "//TRIM(FNBGSI)

	!CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, ERROR)
	if (print_deb) PRINT*, "READ INPUT SFC DATA FROM: "//TRIM(forc_inp_path)

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
	SNOFCS = RESHAPE(DUMMY, (/LENSFC/))

	ERROR=NF90_INQ_VARID(NCID, "snwdph", ID_VAR)
	CALL NETCDF_ERR(ERROR, 'READING snwdph ID' )
	ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
	CALL NETCDF_ERR(ERROR, 'READING snwdph' )
	SWDFCS = RESHAPE(DUMMY, (/LENSFC/))
	
	DEALLOCATE(DUMMY)

	ERROR = NF90_CLOSE(NCID)
    
 END SUBROUTINE READ_Forecast_Data_atPath

 SUBROUTINE READ_Analysis_Data(anl_inp_path, MYRANK, LENSFC, SNOFCS, SWDFCS) 
    
	IMPLICIT NONE

	include "mpif.h"
	
	CHARACTER(LEN=*), Intent(In)      :: anl_inp_path
	INTEGER, INTENT(IN)       :: MYRANK, LENSFC
	REAL, INTENT(OUT)         :: SNOFCS(LENSFC), SWDFCS(LENSFC) 

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
	SNOFCS = RESHAPE(DUMMY, (/LENSFC/))

	ERROR=NF90_INQ_VARID(NCID, "snwdph", ID_VAR)
	CALL NETCDF_ERR(ERROR, 'READING snwdph ID' )
	ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
	CALL NETCDF_ERR(ERROR, 'READING snwdph' )
	SWDFCS = RESHAPE(DUMMY, (/LENSFC/)) 
 
	
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

 SUBROUTINE READ_LAT_LON_OROG_atRank(MYRANK, RLA,RLO,OROG,OROG_UF,TILE_NUM,IDIM,JDIM,IJDIM)
    
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
     REAL, INTENT(OUT)      :: OROG(IJDIM),OROG_UF(IJDIM)
    
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
    
     ERROR=NF90_INQ_VARID(NCID_OROG, 'orog_raw', ID_VAR)
     CALL NETCDF_ERR(ERROR, 'ERROR READING orog_raw ID' )
     ERROR=NF90_GET_VAR(NCID_OROG, ID_VAR, DUMMY4)
     CALL NETCDF_ERR(ERROR, 'ERROR READING orog_raw RECORD' )
     OROG_UF = RESHAPE(DUMMY4, (/IJDIM/))
    
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
    

 END MODULE M_Snow_Analysis
 
