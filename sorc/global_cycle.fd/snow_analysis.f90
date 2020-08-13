MODULE M_Snow_Analysis

USE NETCDF
USE M_DA_OI
!USE MPI
Use, Intrinsic :: IEEE_ARITHMETIC	

Logical, parameter		  :: print_deb = .False.  ! change this to /= 1 to turn off print

CONTAINS
   
subroutine Snow_Analysis(MAX_TASKS, MYRANK, NPROCS, IDIM, JDIM, IY, IM, ID, IH, LENSFC, SNOANL)
							
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

	INTEGER, intent(in) :: MAX_TASKS, MYRANK, NPROCS, IDIM, JDIM, IY, IM, ID, IH, LENSFC
	REAL, intent(Out)   :: SNOANL(LENSFC)   !, SWDANL(LENSFC), anl_fSCA(LENSFC)
	CHARACTER(LEN=5)    :: TILE_NUM
	Character(LEN=3)    :: rank_str
	INTEGER			    :: IERR	
	REAL        :: RLA(LENSFC), RLO(LENSFC), RLO_Tile(LENSFC), OROG(LENSFC), OROG_UF(LENSFC)
	REAL                :: SNOFCS_Inp(LENSFC), SWDFCS(LENSFC)
	REAL                :: SNOFCS(LENSFC), SNWDEN(LENSFC), VETFCS(LENSFC), SWDANL(LENSFC)	
	CHARACTER(len=250)   :: dim_name, ghcnd_inp_file, ims_inp_file, ims_inp_file_indices
	CHARACTER(len=250) 	 :: da_out_file
	CHARACTER(len=5)     :: y_str, m_str, d_Str, fvs_tile
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
	REAL, ALLOCATABLE  :: SNOFCS_atEvalPts(:), innov_atEvalPts(:), SNOANL_atEvalPts(:)  !evalution points 
	REAL, ALLOCATABLE  :: Lat_atEvalPts(:), Lon_atEvalPts(:), Obs_atEvalPts(:)     !evalution points
	Integer, ALLOCATABLE 	:: index_back_atEval(:)     ! background locations at eval points 
	
	Real			   :: ims_threshold      ! threshold for converting IMS fSCA to binary 1, 0 values	
	Real			   :: snwden_val
	LOGICAL 		   :: assim_SWE, assim_GHCND   !assimilate swe? (instead of snow depth)
	LOGICAL 		   :: assim_IMS, assim_IMS_thisGridCell    !assimilate sncov
    ! for mpi par
	INTEGER            :: Np_ext, Np_til, p_tN, p_tRank, N_sA, N_sA_Ext, mp_start, mp_end
	INTEGER            :: send_proc, rec_stat(MPI_STATUS_SIZE), dest_Aoffset, pindex
	INTEGER            :: mpiReal_size, rsize

	assim_SWE = .False.  ! note: if this is set true, may need to adjust background and obs stddev above
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
	
    ! READ THE OROGRAPHY AND GRID POINT LAT/LONS FOR THE CUBED-SPHERE TILE.
	CALL READ_LAT_LON_OROG_atRank(p_tN, RLA,RLO,OROG,OROG_UF,TILE_NUM,IDIM,JDIM,LENSFC)
   ! READ THE INPUT SURFACE DATA ON THE CUBED-SPHERE TILE.
	Call READ_Forecast_Data(p_tN, LENSFC, SNOFCS_Inp, SWDFCS, VETFCS)  !VEGFCS, 
	PRINT*,"Snow anl on ", MYRANK, " Tile group: ", p_tN, " Tile: ", TILE_NUM

! ToDO: use file location the way the rest of the program accesses
	! 4.9.20 for now hard code the file location 
	write(y_str, "(I4)") IY
	write(m_str, "(I0.2)") IM
	write(d_str, "(I0.2)") ID
	write(fvs_tile, "(I3)") IDIM
	! data_dir = /scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/GHCND/GHCND.SNWD.2019100118.nc
	RLO_Tile = RLO
	Where(RLO_Tile > 180) RLO_Tile = RLO_Tile - 360
	lat_min = MAX(MINVAL(RLA) - 1., -90.)
	lat_max = MIN(MAXVAL(RLA) + 1., 90.)
	lon_min = MAX(MINVAL(RLO_Tile) - 1., -180.)
	lon_max = MIN(MAXVAL(RLO_Tile) + 1., 180.)		
	if ((p_tRank==0) .and. print_deb) then
		print*, "Tile ", p_tN, " min/max lat/lon ", lat_min, lat_max, lon_min, lon_max
	endif
	ghcnd_inp_file = "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/GHCND/GHCND.SNWD."// &
						TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18.nc"
	dim_name = "Site_Id"	
	call Observation_Read_GHCND_Tile(ghcnd_inp_file, dim_name, &
					lat_min, lat_max, lon_min, lon_max, & 
					Num_Ghcnd, SNWD_GHCND,		&
					Lat_GHCND, Lon_GHCND, MYRANK)  !Ele_GHCND,		&
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
	ims_inp_file = "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS/IMS.SNCOV."// &
						TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18.nc"                      !
	ims_inp_file_indices = "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS/C"// &
						TRIM(ADJUSTL(fvs_tile))//".IMS.Indices."//TRIM(TILE_NUM)//".nc"			
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
	!ALLOCATE(index_back_atGHCND(Num_Ghcnd)) 
	num_Eval = floor(0.05 * Num_Ghcnd)      ! using 5% of ghcnd locations for evaluation
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
	Call Observation_Operator_Parallel(Myrank, MAX_TASKS, p_tN, p_tRank, Np_til, & 
	                    RLA, RLO, OROG, Lat_GHCND, Lon_GHCND,   &
						LENSFC, Num_Ghcnd, num_Eval, max_distance, SNOFCS, SNWD_GHCND,  &
						SNOFCS_atGHCND, Ele_GHCND, index_back_atEval, Obs_atEvalPts,    &
						SNOFCS_atEvalPts, Lat_atEvalPts, Lon_atEvalPts)
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
			if(num_loc > 0) then     !.and. (SNCOV_IMS(jndx) > 0.)) then    
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
	else
		SNOANL = anl_at_Grid / SNWDEN		!snwden_val 
		SWDANL = anl_at_Grid
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
	Do jndx = 1, num_Eval
		if (index_back_atEval(jndx) > 0) then
			innov_atEvalPts(jndx) = innov_at_Grid(index_back_atEval(jndx))
			SNOANL_atEvalPts(jndx) = anl_at_Grid(index_back_atEval(jndx))
		endif
	End do

	! write outputs
!ToDO: Standard output locations	
	Write(rank_str, '(I3.3)') (MYRANK+1)
	da_out_file = "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/Analysis/SNOANL."// &  !
						TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18_tile"//rank_str//".nc"  !
	
	call Write_DA_Outputs(da_out_file, IDIM, JDIM, LENSFC, MYRANK, &
	                      SNOFCS_Inp, SNOANL, SWDFCS, SWDANL, &  !
			      innov_at_Grid, SNCOV_IMS, &
			      !   Num_Ghcnd, Lat_GHCND, Lon_GHCND, SNWD_GHCND, &
			      !   SNOFCS_atGHCND, SNOFCS_atGHCND, SNOFCS_atGHCND)
			      num_Eval, Lat_atEvalPts, Lon_atEvalPts, Obs_atEvalPts, & 
			      SNOFCS_atEvalPts, innov_atEvalPts, SNOANL_atEvalPts)  !, anl_fSCA) !updated snocov

998 CONTINUE
	DEALLOCATE(SNWD_GHCND, SNOFCS_atGHCND, Lat_GHCND, Lon_GHCND, Ele_GHCND)  !, index_back_atGHCND)
	DEALLOCATE(Obs_atEvalPts, SNOFCS_atEvalPts, innov_atEvalPts, SNOANL_atEvalPts)
	DEALLOCATE(index_back_atEval, Lat_atEvalPts, Lon_atEvalPts) !, Ele_atEvalPts)
	!DEALLOCATE(SNCOV_IMS, Lat_IMS, Lon_IMS, SNOFCS_atIMS, SWDFCS_atIMS)
999 CONTINUE
    PRINT*,'Finished OI DA ON RANK: ', MYRANK
	CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

	!STOP

	RETURN

 END subroutine Snow_Analysis

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
 
 ! Gets obs snow depth from IMS based on exponential/log 'depletion curve' 
SUBROUTINE Observation_Operator_IMS_fSCA(SNCOV_IMS, SNWDEN, VETFCS_in, assim_SWE, LENSFC, 		&
							              SNWD_IMS_at_Grid) !,      &
	
	IMPLICIT NONE
	!
	Real, Intent(In) 	:: SNCOV_IMS(LENSFC), SNWDEN(LENSFC), VETFCS_in(LENSFC)
	Logical, Intent(In)	:: assim_SWE
	INTEGER, Intent(In)	:: LENSFC
	Real, Intent(Out) 	:: SNWD_IMS_at_Grid(LENSFC)
	
	INTEGER             :: VETFCS(LENSFC)
	REAL               :: snupx(30), SNUP, SALP, RSNOW
	Integer 		   :: indx, vtype_int

	!Fill background values to nan (to differentiate those that don't have value)
	SNWD_IMS_at_Grid = IEEE_VALUE(SNWD_IMS_at_Grid, IEEE_QUIET_NAN)

	! call debug_print("Here ", 1.)

	!This is for the IGBP veg classification scheme.
	snupx = (/0.080, 0.080, 0.080, 0.080, 0.080, 0.020, 	&
			0.020, 0.060, 0.040, 0.020, 0.010, 0.020,			&
			0.020, 0.020, 0.013, 0.013, 0.010, 0.020,			&
			0.020, 0.020, 0.000, 0.000, 0.000, 0.000,			&
			0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)

	SALP = -4.0
	VETFCS = INT(VETFCS_in)
	Where(VETFCS==0) VETFCS = 7  !vtype_tile[vtype_tile==0] = 7
	
	Do indx = 1, LENSFC  
		if (.NOT. IEEE_IS_NAN(SNCOV_IMS(indx))) then
			SNUP = snupx(VETFCS(indx))
			if (SNUP == 0.) then
				print*, " 0.0 snup value, check vegclasses", vtype_int
				Stop
			endif

			if (SNCOV_IMS(indx) >= 1.0) then
				SNWD_IMS_at_Grid(indx) = SNUP * 1000.0  ! units mm
			elseif (SNCOV_IMS(indx) < 0.001) then
				SNWD_IMS_at_Grid(indx) = 0.0  ! units mm
			else
				RSNOW = LOG(1. - SNCOV_IMS(indx)) / SALP
				if (RSNOW > 1.) RSNOW = 1.
				SNWD_IMS_at_Grid(indx) = RSNOW * SNUP * 1000.   ! units mm
			endif

		endif
	end do	
	! print*, "IMS Sndpth at model grids"
	! print*, SNWD_IMS_at_Grid
	! print*
	! assim_SWE = True >> swe assimilated; SNOFCS  = SWDFCS / SNWDEN
	if (assim_SWE) SNWD_IMS_at_Grid = SNWD_IMS_at_Grid / SNWDEN
	! print*, "IMS SWE at model grids"
	! print*, SNWD_IMS_at_Grid

	RETURN

END SUBROUTINE Observation_Operator_IMS_fSCA

 ! Gets obs snow depth from IMS based on threshold fSCA 
SUBROUTINE Observation_Operator_IMS_fSCA_Threshold(SNCOV_IMS, SNOFCS, SNWDEN, assim_SWE,	   &
							LENSFC, ims_threshold, 		&
							SNWD_IMS_at_Grid) !,      &

    IMPLICIT NONE
	!
	Real, Intent(In) 	:: SNCOV_IMS(LENSFC), SNOFCS(LENSFC), SNWDEN(LENSFC)
	Logical				:: assim_SWE
	INTEGER 			:: LENSFC
	Real, Intent(In)    :: ims_threshold

	Real, Intent(Out) 	:: SNWD_IMS_at_Grid(LENSFC)
	
	INTEGER :: indx
	Real, Parameter		:: SWE_Tolerance = 0.001    ! smallest swe value

	!Fill background values to nan (to differentiate those that don't have value)
	SNWD_IMS_at_Grid = IEEE_VALUE(SNWD_IMS_at_Grid, IEEE_QUIET_NAN)

    ! call debug_print("Here ", 1.)

	Do indx = 1, LENSFC  
		if (.NOT. IEEE_IS_NAN(SNCOV_IMS(indx))) then
			if (SNCOV_IMS(indx) >= ims_threshold) then
				! ims snow, model no snow => obs=50 mm
				if (SNOFCS(indx) < SWE_Tolerance) SNWD_IMS_at_Grid(indx) = 50.
				! ims snow, model snow => no assimilation
			else  ! IMS fSCA < thresh => Ims obs = 0
				SNWD_IMS_at_Grid(indx) = 0.
			endif   ! all others nan	
		endif
	end do
	
	! print*, "IMS Sndpth at model grids"
	! print*, SNWD_IMS_at_Grid
	! print*
	! assim_SWE = True >> swe assimilated; SNOFCS  = SWDFCS / SNWDEN
	if (assim_SWE) SNWD_IMS_at_Grid = SNWD_IMS_at_Grid / SNWDEN

	! print*, "IMS SWE at model grids"
	! print*, SNWD_IMS_at_Grid

	RETURN

END SUBROUTINE Observation_Operator_IMS_fSCA_Threshold

subroutine resample_to_model_tiles_intrp(data_grid_ims, data_grid_ims_ind, &
                                            nlat_ims, nlon_ims, n_lat, n_lon, num_sub, &  !myrank, &
											grid_dat)
											
    Use, Intrinsic :: IEEE_ARITHMETIC

    Implicit None

    Integer, Intent(In)     :: nlat_ims, nlon_ims, n_lat, n_lon, num_sub 
    Integer, Intent(In)     :: data_grid_ims(nlon_ims, nlat_ims), data_grid_ims_ind(num_sub, n_lon, n_lat) 
    Real, Intent(Out)       :: grid_dat(n_lon, n_lat)

	Integer   :: jc, jy, ix, num_loc_counter
	Integer   :: lonlatcoord_ims, loncoord_ims, latcoord_ims
    
	grid_dat = IEEE_VALUE(grid_dat, IEEE_QUIET_NAN)

    Do jy=1, n_lat
    !print*, "process: ", myrank, "loop ", indx

        Do ix=1, n_lon
            
            num_loc_counter = data_grid_ims_ind(1, ix, jy)
            if (num_loc_counter < 1) then 
                !print*, "no matching values!"
                cycle
            end if
            !print*, "jy ", jy, " ix ", ix
            grid_dat(ix, jy) = 0.
			Do jc = 2, num_loc_counter+1
				lonlatcoord_ims = data_grid_ims_ind(jc, ix, jy) - 1 
				latcoord_ims = lonlatcoord_ims / nlon_ims + 1
				loncoord_ims = mod(lonlatcoord_ims, nlon_ims) + 1
				if(latcoord_ims > nlat_ims) then
					latcoord_ims = nlat_ims
					print*, "Warning! lat coordinate outside domain boundary"
				endif
				if(loncoord_ims > nlon_ims) then
					loncoord_ims = nlon_ims
					print*, "Warning! lon coordinate outside domain boundary"
				endif
                grid_dat(ix, jy) =  grid_dat(ix, jy) + data_grid_ims(loncoord_ims, latcoord_ims)              
            End do

            grid_dat(ix, jy) =  grid_dat(ix, jy) / num_loc_counter ! first location, num obs

        End do

    End do

    return !grid_dat

End subroutine resample_to_model_tiles_intrp

 !This reads the whole IMS file and uses a-priori prepared indices to sample those wihin the grid cel
 SUBROUTINE Observation_Read_IMS_Full(inp_file, inp_file_indices, &
				MYRANK, n_lat, n_lon, num_sub, &
				SNCOV_IMS)
				! Ele		&
				
	IMPLICIT NONE

	include 'mpif.h'		  

	!ToDO: Can you use variable length char array ?
	CHARACTER(LEN=*), Intent(In)   :: inp_file, inp_file_indices !, dim_name
	INTEGER, Intent(In)            :: MYRANK, n_lat, n_lon, num_sub
	! ToDO: ims snow cover is of 'byte' type (Chcek the right one)	
	Real, Intent(Out)       :: SNCOV_IMS(n_lat * n_lon) 	

	INTEGER, ALLOCATABLE    :: SNCOV_IMS_2D_full(:,:)    !SNCOV_IMS_1D(:), 
	Integer                 :: data_grid_ims_ind(num_sub, n_lon, n_lat) 
	Real                    :: grid_dat(n_lon, n_lat)
	
	INTEGER                :: ERROR, NCID, ID_DIM, ID_VAR, DIM_LEN, DIM_LEN_lat, DIM_LEN_lon

	ERROR=NF90_OPEN(TRIM(inp_file),NF90_NOWRITE,NCID)
	CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(inp_file) )

	ERROR=NF90_INQ_DIMID(NCID, 'lat', ID_DIM)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Dimension lat' )
	
	ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=DIM_LEN_lat)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Dimension Lat' )
	
	ERROR=NF90_INQ_DIMID(NCID, 'lon', ID_DIM)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Dimension lon' )

	ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=DIM_LEN_lon)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Dimension Lon' )

    ALLOCATE(SNCOV_IMS_2D_full(DIM_LEN_lon, DIM_LEN_lat))	
	! print*, "initial IMS array size (lon, lat)= ", DIM_LEN_lon, " ",DIM_LEN_lat
	!ALLOCATE(Ele(DIM_LEN_lat, DIM_LEN_lon))
	ERROR=NF90_INQ_VARID(NCID, 'Band1', ID_VAR)
	CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV ID' )
	ERROR=NF90_GET_VAR(NCID, ID_VAR, SNCOV_IMS_2D_full, start = (/ 1, 1 /), &
							count = (/ DIM_LEN_lon, DIM_LEN_lat/))
	CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV RECORD' )
	! ERROR=NF90_INQ_VARID(NCID, 'lat', ID_VAR)
	! CALL NETCDF_ERR(ERROR, 'ERROR READING Lat var ID' )
	! ERROR=NF90_GET_VAR(NCID, ID_VAR, Lat_IMS_1D)
	! CALL NETCDF_ERR(ERROR, 'ERROR READING Lat RECORD' )
	! ERROR=NF90_INQ_VARID(NCID, 'lon', ID_VAR)
	! CALL NETCDF_ERR(ERROR, 'ERROR READING Lon var ID' )
	! ERROR=NF90_GET_VAR(NCID, ID_VAR, Lon_IMS_1D)
	! CALL NETCDF_ERR(ERROR, 'ERROR READING Lon RECORD' )

	! need to read corresponding elevation values 
	! ERROR=NF90_INQ_VARID(NCID, 'Elevation', ID_VAR)
	! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation ID' )
	! ERROR=NF90_GET_VAR(NCID, ID_VAR, Ele)
	! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation RECORD' )
	
	ERROR = NF90_CLOSE(NCID)

	ERROR=NF90_OPEN(TRIM(inp_file_indices),NF90_NOWRITE, NCID)
	CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(inp_file_indices) )

	ERROR=NF90_INQ_VARID(NCID, 'IMS_Indices', ID_VAR)
	CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV Indices ID' )
	ERROR=NF90_GET_VAR(NCID, ID_VAR, data_grid_ims_ind, start = (/ 1, 1, 1 /), &
							count = (/ num_sub, n_lon, n_lat/))
	CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV Indices' )

	ERROR = NF90_CLOSE(NCID)

	! print*, "IMS 2D"
	! print*, SNCOV_IMS_2D_full
	! print*
	
	! DIM_LEN = DIM_LEN_lat * DIM_LEN_lon
	! ALLOCATE(SNCOV_IMS_1D(DIM_LEN))
	! SNCOV_IMS_1D = Reshape(SNCOV_IMS_2D_full, (/DIM_LEN/))

	! print*, "IMS 1D"
	! print*, SNCOV_IMS_1D
	! print*
    Where(SNCOV_IMS_2D_full /= 4) SNCOV_IMS_2D_full = 0
	Where(SNCOV_IMS_2D_full == 4) SNCOV_IMS_2D_full = 1
	

	! print*, "IMS binary"
	! print*, SNCOV_IMS_1D
	! print*

	! print*, "IMS Indices 3D"
	! print*, data_grid_ims_ind
	! print*

	call resample_to_model_tiles_intrp(SNCOV_IMS_2D_full, data_grid_ims_ind, &
	                                   DIM_LEN_lat, DIM_LEN_lon, n_lat, n_lon, num_sub, &  !myrank, &
		                               grid_dat)

	SNCOV_IMS = Reshape(grid_dat, (/n_lat * n_lon/))

	DEALLOCATE(SNCOV_IMS_2D_full)
			  
	RETURN
	
 End SUBROUTINE Observation_Read_IMS_Full

 SUBROUTINE Observation_Read_IMS_tile(inp_file, MYRANK, n_lat, n_lon, &
									  SNCOV_IMS)
									  		
	IMPLICIT NONE

	include 'mpif.h'		  

	!ToDO: Can you use variable length char array ?
	CHARACTER(LEN=*), Intent(In)   :: inp_file
	INTEGER, Intent(In)            :: MYRANK, n_lat, n_lon
	! ToDO: ims snow cover is of 'byte' type (Chcek the right one)	
	Real, Intent(Out)       :: SNCOV_IMS(n_lat * n_lon) 
	Real                    :: SNCOV_IMS_2D_full( n_lon, n_lat)
	
	INTEGER                :: ERROR, NCID, ID_DIM, ID_VAR, DIM_LEN, DIM_LEN_lat, DIM_LEN_lon

	ERROR=NF90_OPEN(TRIM(inp_file),NF90_NOWRITE,NCID)
	CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(inp_file) )

	!ALLOCATE(Ele(DIM_LEN_lat, DIM_LEN_lon))
	ERROR=NF90_INQ_VARID(NCID, 'imsfSCA', ID_VAR)
	CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV ID' )
	ERROR=NF90_GET_VAR(NCID, ID_VAR, SNCOV_IMS_2D_full, start = (/ 1, 1 /), &
							count = (/n_lon, n_lat/))
	CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV RECORD' )

	! print*, "IMS 2D"
	! print*, SNCOV_IMS_2D_full
	! print*
	
	SNCOV_IMS = Reshape(SNCOV_IMS_2D_full, (/n_lat * n_lon/))
	! print*, "IMS 1D"
	! print*, SNCOV_IMS
	! print*
			  
	RETURN
	
 End SUBROUTINE Observation_Read_IMS_tile

 SUBROUTINE Observation_Read_IMS(inp_file, MYRANK,  &
	            lat_min, lat_max, lon_min, lon_max, &
				DIM_LEN,       & !_lat, DIM_LEN_lon,  &
				SNCOV_IMS, Lat_IMS, Lon_IMS)
				! Ele		&
				
	IMPLICIT NONE

	include 'mpif.h'		  

	!ToDO: Can you use variable length char array ?
	CHARACTER(LEN=*), Intent(In)   :: inp_file !, dim_name
	INTEGER, Intent(In)            :: MYRANK
    REAL, Intent(In)         	   :: lat_min, lat_max, lon_min, lon_max    
	! ToDO: ims snow cover is of 'byte' type (Chcek the right one)
	INTEGER, Intent(Out)    :: DIM_LEN 	
	INTEGER, ALLOCATABLE, Intent(Out)    :: SNCOV_IMS(:) 	
	REAL, ALLOCATABLE, Intent(Out)	   :: Lat_IMS(:), Lon_IMS(:)	!, Ele(:,:)

	INTEGER, ALLOCATABLE    :: SNCOV_IMS_2D(:,:), SNCOV_IMS_2D_full(:,:)
	REAL, ALLOCATABLE	   :: Lat_IMS_1D(:), Lon_IMS_1D(:)	
	REAL, ALLOCATABLE	   :: Lat_IMS_2D(:,:), Lon_IMS_2D(:,:)
	REAL, ALLOCATABLE	   :: Lat_minmax_Diff(:), Lon_minmax_Diff(:)	
	
	INTEGER                :: ERROR, NCID, ID_DIM, ID_VAR, DIM_LEN_lat, DIM_LEN_lon
	INTEGER                :: indx, jndx, minlat_indx, maxlat_indx, minlon_indx, maxlon_indx
	INTEGER                :: icounter, jcounter, iincr, jincr

	ERROR=NF90_OPEN(TRIM(inp_file),NF90_NOWRITE,NCID)
	CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(inp_file) )

	ERROR=NF90_INQ_DIMID(NCID, 'lat', ID_DIM)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Dimension lat' )
	
	ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=DIM_LEN_lat)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Dimension Lat' )
	
	ERROR=NF90_INQ_DIMID(NCID, 'lon', ID_DIM)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Dimension lon' )

	ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=DIM_LEN_lon)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Dimension Lon' )

    ALLOCATE(SNCOV_IMS_2D_full(DIM_LEN_lon, DIM_LEN_lat))	
	ALLOCATE(Lat_IMS_1D(DIM_LEN_lat))
	ALLOCATE(Lon_IMS_1D(DIM_LEN_lon))

	! print*, "initial IMS array size (lon, lat)= ", DIM_LEN_lon, " ",DIM_LEN_lat
	!ALLOCATE(Ele(DIM_LEN_lat, DIM_LEN_lon))

	ERROR=NF90_INQ_VARID(NCID, 'Band1', ID_VAR)
	CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV ID' )
	ERROR=NF90_GET_VAR(NCID, ID_VAR, SNCOV_IMS_2D_full, start = (/ 1, 1 /), &
							count = (/ DIM_LEN_lon, DIM_LEN_lat/))
	CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV RECORD' )
	
	ERROR=NF90_INQ_VARID(NCID, 'lat', ID_VAR)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Lat var ID' )
	ERROR=NF90_GET_VAR(NCID, ID_VAR, Lat_IMS_1D)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Lat RECORD' )

	ERROR=NF90_INQ_VARID(NCID, 'lon', ID_VAR)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Lon var ID' )
	ERROR=NF90_GET_VAR(NCID, ID_VAR, Lon_IMS_1D)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Lon RECORD' )

	! find indices within grid cell boundaries
	ALLOCATE(Lat_minmax_Diff(DIM_LEN_lat))
	ALLOCATE(Lon_minmax_Diff(DIM_LEN_lon))
	! find index with min abs value difference
	Lat_minmax_Diff = abs(Lat_IMS_1D - lat_min)
	minlat_indx = MINLOC(Lat_minmax_Diff, dim = 1)
	Lat_minmax_Diff = abs(Lat_IMS_1D - lat_max)
	maxlat_indx = MINLOC(Lat_minmax_Diff, dim = 1)
	Lon_minmax_Diff = abs(Lon_IMS_1D - lon_min)
	minlon_indx = MINLOC(Lon_minmax_Diff, dim =1 )
	Lon_minmax_Diff = abs(Lon_IMS_1D - lon_max)
	maxlon_indx = MINLOC(Lon_minmax_Diff, dim =1 )
	
	DIM_LEN_lat = 1 + abs(maxlat_indx - minlat_indx)
	DIM_LEN_lon = 1 + abs(maxlon_indx - minlon_indx)
	! print*, "New (Tile-specific) IMS array size = ", DIM_LEN_lon, " ",DIM_LEN_lat
	
	DIM_LEN = DIM_LEN_lat * DIM_LEN_lon
	ALLOCATE(SNCOV_IMS(DIM_LEN))
	ALLOCATE(Lat_IMS(DIM_LEN))
	ALLOCATE(Lon_IMS(DIM_LEN))
	ALLOCATE(Lat_IMS_2D(DIM_LEN_lon, DIM_LEN_lat))
	ALLOCATE(Lon_IMS_2D(DIM_LEN_lon, DIM_LEN_lat))
	ALLOCATE(SNCOV_IMS_2D(DIM_LEN_lon, DIM_LEN_lat))	
	
	! print*," Dims of SNCOV_IMS_2D_full ", size(SNCOV_IMS_2D_full,1), size(SNCOV_IMS_2D_full,2)
	iincr = 1; jincr = 1;
	if (maxlat_indx < minlat_indx) jincr = -1
	if (maxlon_indx < minlon_indx) iincr = -1
	! print*, "lat subarray indices = ", minlat_indx, " ", maxlat_indx
	! print*, "lon subarray indices = ", minlon_indx, " ", maxlon_indx
	! print*, "jincr, iincr = ", jincr, " ", iincr
	jcounter = 1
	Do jndx=minlat_indx, maxlat_indx, jincr
		!print*, "jcounter ", jcounter
		icounter = 1
		Do indx=minlon_indx, maxlon_indx, iincr
			!print*, "icounter ", icounter
			Lat_IMS_2D(icounter, jcounter) = Lat_IMS_1D(jndx)
			Lon_IMS_2D(icounter, jcounter) = Lon_IMS_1D(indx)
			SNCOV_IMS_2D(icounter, jcounter) = SNCOV_IMS_2D_full(indx, jndx)
			icounter = icounter + 1
		End do
		jcounter = jcounter + 1
	End do
	Lat_IMS = Reshape(Lat_IMS_2D, (/DIM_LEN/))
	Lon_IMS = Reshape(Lon_IMS_2D, (/DIM_LEN/))
	SNCOV_IMS = Reshape(SNCOV_IMS_2D, (/DIM_LEN/))

	! need to read corresponding elevation values 
	! ERROR=NF90_INQ_VARID(NCID, 'Elevation', ID_VAR)
	! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation ID' )
	! ERROR=NF90_GET_VAR(NCID, ID_VAR, Ele)
	! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation RECORD' )
	
	ERROR = NF90_CLOSE(NCID)

	DEALLOCATE(SNCOV_IMS_2D, SNCOV_IMS_2D_full)
	DEALLOCATE(Lat_IMS_1D, Lon_IMS_1D, Lat_IMS_2D, Lon_IMS_2D)
	DEALLOCATE(Lat_minmax_Diff, Lon_minmax_Diff)	
			  
	RETURN
	
 End SUBROUTINE Observation_Read_IMS

 ! Get model states at obs points
 ! Warning: This assumes all distance coordinates are valid; 
 ! do quality control of coordinates beforehand
 SUBROUTINE Observation_Operator(RLA, RLO, OROG, Lat_Obs, Lon_Obs,   &
						LENSFC, num_Obs, max_distance, 		&
						SNWD_back,  				&
						SNWD_atObs, Ele_atObs, index_back_atObs) !,      &
						!intp_mode) 
						!  SWE_atObs  SWE_back

    IMPLICIT NONE
	!
	!USE intrinsic::ieee_arithmetic
	Real, Intent(In) 	:: RLA(LENSFC), RLO(LENSFC), OROG(LENSFC)
	Real, Intent(In) 	:: Lat_Obs(num_Obs), Lon_Obs(num_Obs)  ! don't want to alter these
	INTEGER :: LENSFC, num_Obs
	Real	:: max_distance   ! radius_of_influence
	Real, Intent(In) 	:: SNWD_back(LENSFC)

	Real, Intent(Out) 	:: SNWD_atObs(num_Obs), Ele_atObs(num_Obs)
	Integer, Intent(Out) 	:: index_back_atObs(num_Obs)   ! the location of the corresponding obs
	
	Real 	::  Lon_Obs_2(num_Obs)		!RLO_2(LENSFC), 	
	Real 	:: RLA_rad(LENSFC), RLO_rad(LENSFC)
	Real 	:: Lat_Obs_rad(num_Obs), Lon_Obs_rad(num_Obs)	
	INTEGER :: indx, jndx, zndx, min_indx
	Real    :: distArr(LENSFC), haversinArr(LENSFC)
	Real 	:: d_latArr(LENSFC), d_lonArr(LENSFC)
	Real(16), Parameter :: PI_16 = 4 * atan (1.0_16)	
	Real(16), Parameter :: pi_div_180 = PI_16/180.0
	Real, Parameter		:: earth_rad = 6371.
	! PRINT*, "PI: ", PI_16
	! PRINT*, "PI / 180: ", pi_div_180

	!Fill background values to nan (to differentiate those htat don't have value)
	SNWD_atObs = IEEE_VALUE(SNWD_atObs, IEEE_QUIET_NAN)	
	Ele_atObs = IEEE_VALUE(Ele_atObs, IEEE_QUIET_NAN)	
	index_back_atObs = -1   ! when corresponding value doesn't exit
	
	!if intp_mode == 'near'		! [bilinear, customInterpol])

	! RLO from 0 to 360 (no -ve lon)
	Lon_Obs_2 = Lon_Obs
	Where(Lon_Obs_2 < 0) Lon_Obs_2 = 360. + Lon_Obs_2
	! Do zndx = 1, num_Obs 
	! 	if (Lon_Obs_2(zndx) < 0) Lon_Obs_2(zndx) = 360. + Lon_Obs_2(zndx)
	! end do

	! at each obs point compute its distance from RLA/RLO pairs 
	! then find the position of the minimum

	! shortest distance over sphere using great circle distance 	
	RLA_rad =  pi_div_180 * RLA
	RLO_rad =  pi_div_180 * RLO
	Lat_Obs_rad =  pi_div_180 * Lat_Obs
	Lon_Obs_rad =  pi_div_180 * Lon_Obs_2   
	
    ! https://en.wikipedia.org/wiki/Haversine_formula
    ! https://www.geeksforgeeks.org/program-distance-two-points-earth/
	! Distance, d = R * arccos[(sin(lat1) * sin(lat2)) + cos(lat1) * cos(lat2) * cos(long2 – long1)]
	! dist = 2 * R * asin { sqrt [sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2]}
	Do indx = 1, num_Obs 
		d_latArr = (Lat_Obs_rad(indx) - RLA_rad) / 2.
		d_lonArr = (Lon_Obs_rad(indx) - RLO_rad) / 2.
		haversinArr = sin(d_latArr)**2 + cos(Lat_Obs_rad(indx)) * cos(RLA_rad) * sin(d_lonArr)**2
		WHERE(haversinArr > 1) haversinArr = 1.   ! ensure numerical errors don't make h>1
		
		distArr = 2 * earth_rad * asin(sqrt(haversinArr))		
		!distArr = (Lat_Obs(indx) - RLA)**2 + (Lon_Obs_2(indx) - RLO)**2 
		min_indx = MINLOC(distArr, dim = 1)  !, MASK=ieee_is_nan(distArr))

		if(distArr(min_indx) < max_distance) then
			SNWD_atObs(indx) = SNWD_back(min_indx) 
			Ele_atObs(indx) = OROG(min_indx)
			index_back_atObs(indx) = min_indx
		! else
			! Print*, " Warning! distance greater than ",max_distance," km ", distArr(min_indx)
		endif
	end do
	
	
	RETURN
	
 END SUBROUTINE Observation_Operator

 SUBROUTINE Observation_Operator_Parallel(Myrank, MAX_TASKS, p_tN, p_tRank, Np_til, & 
	                    RLA, RLO, OROG, Lat_Obs, Lon_Obs,               &
						LENSFC, num_Obs, num_Eval, max_distance, SNOFCS_back, SNWD_GHCND,  &
						SNOFCS_atObs, Ele_atObs, index_back_atEval, Obs_atEvalPts,      &
						SNOFCS_atEvalPts, Lat_atEvalPts, Lon_atEvalPts)
						
    IMPLICIT NONE
	!
	!USE intrinsic::ieee_arithmetic
	include "mpif.h"

	Real, Intent(In) 	:: RLA(LENSFC), RLO(LENSFC), OROG(LENSFC)
	Real, Intent(In) 	:: Lat_Obs(num_Obs), Lon_Obs(num_Obs)  ! don't want to alter these
	INTEGER             :: Myrank, MAX_TASKS, p_tN, p_tRank, Np_til, LENSFC, num_Obs, num_Eval
	Real	            :: max_distance   ! radius_of_influence
	Real, Intent(In) 	:: SNOFCS_back(LENSFC)
	Real, Intent(InOut) 	:: SNWD_GHCND(num_Obs) 
	Real, Intent(Out) 	    :: SNOFCS_atObs(num_Obs), Ele_atObs(num_Obs), Obs_atEvalPts(num_Eval)
	Integer, Intent(Out) 	:: index_back_atEval(num_Eval)   ! the location of evaluation points
	Real, Intent(Out) 	    :: SNOFCS_atEvalPts(num_Eval), Lat_atEvalPts(num_Eval), Lon_atEvalPts(num_Eval)
	
	Integer	:: index_back_atObs(num_Obs)   ! the location of background corresponding obs
	Real 	:: Lon_Obs_2(num_Obs)		!RLO_2(LENSFC), 	
	Real 	:: RLA_rad(LENSFC), RLO_rad(LENSFC)
	Real 	:: Lat_Obs_rad(num_Obs), Lon_Obs_rad(num_Obs)	
	INTEGER :: indx, jndx, jzndx, zndx, min_indx
	Real    :: distArr(LENSFC), haversinArr(LENSFC)
	Real 	:: d_latArr(LENSFC), d_lonArr(LENSFC)
	Real(16), Parameter :: PI_16 = 4 * atan (1.0_16)	
	Real(16), Parameter :: pi_div_180 = PI_16/180.0
	Real, Parameter		:: earth_rad = 6371.
	
	! for mpi par
	INTEGER            :: N_sA, N_sA_Ext, mp_start, mp_end 
	INTEGER            :: send_proc, rec_proc, rec_stat(MPI_STATUS_SIZE), dest_Aoffset, pindex
	INTEGER            :: mpiReal_size, rsize, mpiInt_size, isize, IERR
	Real               :: rand_nextVal  ! randomly select evalution points to exclude from DA
	Integer            :: rand_evalPoint(num_Eval)    ! randomly select evalution points to exclude from DA

	!Np_til ! num proc. per tile p_tRank ! proc. rank within tile !p_tN  ! tile for proc.
	N_sA = num_Obs / Np_til  ! sub array length per proc
	N_sA_Ext = num_Obs - N_sA * Np_til ! extra grid cells
	if(p_tRank == 0) then 
		mp_start = 1
	else
		mp_start = p_tRank * N_sA + N_sA_Ext + 1   ! start index of subarray for proc
	endif
	mp_end = (p_tRank + 1) * N_sA + N_sA_Ext 		! end index of subarray for proc

	!Fill background values to nan (to differentiate those htat don't have value)
	SNOFCS_atObs = IEEE_VALUE(SNOFCS_atObs, IEEE_QUIET_NAN)	
	Ele_atObs = IEEE_VALUE(Ele_atObs, IEEE_QUIET_NAN)	
	index_back_atObs = -1   ! when corresponding value doesn't exit	
	rand_evalPoint = -1
	index_back_atEval = -1

	! RLO from 0 to 360 (no -ve lon)
	Lon_Obs_2 = Lon_Obs
	Where(Lon_Obs_2 < 0) Lon_Obs_2 = 360. + Lon_Obs_2
	! Do zndx = 1, num_Obs 
	!     !Lon_Obs[Lon_Obs<0]= 360.0 + Lon_Obs[Lon_Obs<0]
	! 	if (Lon_Obs(zndx) < 0) Lon_Obs_2(zndx) = 360. + Lon_Obs(zndx)
	! end do
	! at each obs point compute its distance from RLA/RLO pairs 
	! then find the position of the minimum

	! shortest distance over sphere using great circle distance 	
	RLA_rad =  pi_div_180 * RLA
	RLO_rad =  pi_div_180 * RLO
	Lat_Obs_rad =  pi_div_180 * Lat_Obs
	Lon_Obs_rad =  pi_div_180 * Lon_Obs_2   	
    ! https://en.wikipedia.org/wiki/Haversine_formula
    ! https://www.geeksforgeeks.org/program-distance-two-points-earth/
	! Distance, d = R * arccos[(sin(lat1) * sin(lat2)) + cos(lat1) * cos(lat2) * cos(long2 – long1)]
	! dist = 2 * R * asin { sqrt [sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2]}
	Do indx = mp_start, mp_end   !1, num_Obs 
		d_latArr = (Lat_Obs_rad(indx) - RLA_rad) / 2.
		d_lonArr = (Lon_Obs_rad(indx) - RLO_rad) / 2.
		haversinArr = sin(d_latArr)**2 + cos(Lat_Obs_rad(indx)) * cos(RLA_rad) * sin(d_lonArr)**2
		WHERE(haversinArr > 1) haversinArr = 1.   ! ensure numerical errors don't make h>1
		
		distArr = 2 * earth_rad * asin(sqrt(haversinArr))		
		!distArr = (Lat_Obs(indx) - RLA)**2 + (Lon_Obs_2(indx) - RLO)**2 
		min_indx = MINLOC(distArr, dim = 1)  !, MASK=ieee_is_nan(distArr))

		if(distArr(min_indx) < max_distance) then
			SNOFCS_atObs(indx) = SNOFCS_back(min_indx) 
			Ele_atObs(indx) = OROG(min_indx)
			index_back_atObs(indx) = min_indx
		! else
			! Print*, " Warning! distance greater than ",max_distance," km ", distArr(min_indx)
		endif
	end do

! ToDO: Better way to handle this?
! Real data type size corresponding to mpi
	rsize = SIZEOF(max_distance) 
	Call MPI_TYPE_SIZE(MPI_REAL, mpiReal_size, IERR) 
	If (rsize == 4 ) then 
		mpiReal_size = MPI_REAL4
	elseif (rsize == 8 ) then 
		mpiReal_size = MPI_REAL8
	elseif (rsize == 16 ) then 
		mpiReal_size = MPI_REAL16
	else
		PRINT*," Possible mismatch between Fortran Real ", rsize," and Mpi Real", mpiReal_size
		Stop
	endif
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
		call MPI_SEND(SNOFCS_atObs(mp_start:mp_end), N_sA, mpiReal_size, p_tN,   &
					  MYRANK, MPI_COMM_WORLD, IERR) 
		call MPI_SEND(Ele_atObs(mp_start:mp_end), N_sA, mpiReal_size, p_tN,   &
					  MYRANK*100, MPI_COMM_WORLD, IERR)
		call MPI_SEND(index_back_atObs(mp_start:mp_end), N_sA, mpiInt_size, p_tN,   &
					  MYRANK*1000, MPI_COMM_WORLD, IERR)
	else !if (MYRANK == p_tN ) then  
		Do pindex =  1, (Np_til - 1)   ! sender proc index within tile group
			dest_Aoffset = pindex * N_sA + N_sA_Ext + 1   ! dest array offset
			send_proc = MYRANK +  pindex * MAX_TASKS
			call MPI_RECV(SNOFCS_atObs(dest_Aoffset:dest_Aoffset+N_sA-1), N_sA, mpiReal_size, send_proc,  &
					  send_proc, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
			call MPI_RECV(Ele_atObs(dest_Aoffset:dest_Aoffset+N_sA-1), N_sA, mpiReal_size, send_proc,   &
					  send_proc*100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
			call MPI_RECV(index_back_atObs(dest_Aoffset:dest_Aoffset+N_sA-1), N_sA, mpiInt_size, send_proc, &
					  send_proc*1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
		enddo
	endif
!ToDO: better way to do this?
	! now share the whole array
	if (MYRANK < MAX_TASKS ) then   !if (MYRANK == p_tN ) then 	
		! Print*, "Started selecting obs points"
		! Select obs points to exclude from DA 
		Call random_number(rand_nextVal)
		rand_evalPoint(1) = floor(rand_nextVal * num_Obs) + 1
		index_back_atEval(1) = index_back_atObs(rand_evalPoint(1))
		jndx = 2
		jzndx = 1
		Do  While (jndx <= num_Eval)        !jndx = 2, num_Eval
			Call random_number(rand_nextVal)
			! Print*, rand_nextVal
			rand_evalPoint(jndx) = floor(rand_nextVal * num_Obs) + 1
			if((rand_evalPoint(jndx) /= rand_evalPoint(jndx-1)) .AND.      &
			   (.NOT. IEEE_IS_NAN(SNWD_GHCND(rand_evalPoint(jndx))/SNWD_GHCND(rand_evalPoint(jndx)))) )  then
				index_back_atEval(jndx) = index_back_atObs(rand_evalPoint(jndx))
				jndx = jndx + 1
			else
				cycle
			endif
			jzndx = jzndx + 1
			If (jzndx >= 2*num_Eval) exit 
		Enddo 
		! Print*, "Finished selecting obs points"	
		Do pindex =  1, (Np_til - 1)   ! receiving proc index within tile group
			rec_proc = MYRANK +  pindex * MAX_TASKS
			call MPI_SEND(SNOFCS_atObs, num_Obs, mpiReal_size, rec_proc, MYRANK, MPI_COMM_WORLD, IERR) 
		    call MPI_SEND(Ele_atObs, num_Obs, mpiReal_size, rec_proc, MYRANK*100, MPI_COMM_WORLD, IERR)
			call MPI_SEND(index_back_atObs, num_Obs, mpiInt_size, rec_proc, MYRANK*1000, MPI_COMM_WORLD, IERR)
			call MPI_SEND(index_back_atEval, num_Eval, mpiInt_size, rec_proc, MYRANK*10000, MPI_COMM_WORLD, IERR)
			call MPI_SEND(rand_evalPoint, num_Eval, mpiInt_size, rec_proc, MYRANK*100000, MPI_COMM_WORLD, IERR)
		enddo
	else 
		call MPI_RECV(SNOFCS_atObs, num_Obs, mpiReal_size, p_tN, p_tN, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
		call MPI_RECV(Ele_atObs, num_Obs, mpiReal_size, p_tN, p_tN*100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
		call MPI_RECV(index_back_atObs, num_Obs, mpiInt_size, p_tN, p_tN*1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
		call MPI_RECV(index_back_atEval, num_Eval, mpiInt_size, p_tN, p_tN*10000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
		call MPI_RECV(rand_evalPoint, num_Eval, mpiInt_size, p_tN, p_tN*100000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
	endif
	Obs_atEvalPts = IEEE_VALUE(rand_nextVal, IEEE_QUIET_NAN)
	SNOFCS_atEvalPts = IEEE_VALUE(rand_nextVal, IEEE_QUIET_NAN)
	Lat_atEvalPts = IEEE_VALUE(rand_nextVal, IEEE_QUIET_NAN)
	Lon_atEvalPts = IEEE_VALUE(rand_nextVal, IEEE_QUIET_NAN)
	Do  jndx = 1, num_Eval
		if (index_back_atEval(jndx) > 0) then
			Obs_atEvalPts(jndx) = SNWD_GHCND(rand_evalPoint(jndx))
			SNOFCS_atEvalPts(jndx) = SNOFCS_atObs(rand_evalPoint(jndx))
			Lat_atEvalPts(jndx) = Lat_Obs(rand_evalPoint(jndx)) 
			Lon_atEvalPts(jndx) = Lon_Obs(rand_evalPoint(jndx))
			SNWD_GHCND(rand_evalPoint(jndx)) = IEEE_VALUE(rand_nextVal, IEEE_QUIET_NAN) ! exclude point from DA	
		endif	
	Enddo 	
	! if (Myrank == 0) then 
	! 	Print*, "Started selecting obs points"
	! 	Call random_number(rand_nextVal)
	! 	Print*, rand_nextVal
	! 	rand_evalPoint(1) = floor(rand_nextVal * num_Obs) + 1
	! 	index_back_atEval(1) = index_back_atObs(rand_evalPoint(1))
	! 	jndx = 2
	! 	Do  While (jndx <= num_Eval)        !jndx = 2, num_Eval
	! 		Call random_number(rand_nextVal)
	! 		Print*, rand_nextVal
	! 		rand_evalPoint(jndx) = floor(rand_nextVal * num_Obs) + 1
	! 		if( (rand_evalPoint(jndx) /= rand_evalPoint(jndx-1)) .AND. &
	! 			(.NOT. IEEE_IS_NAN(SNOFCS_atObs(rand_evalPoint(jndx)))) .AND. &
	! 			(.NOT. IEEE_IS_NAN(SNWD_GHCND(rand_evalPoint(jndx))/SNWD_GHCND(rand_evalPoint(jndx)))) ) then
	! 			index_back_atEval(jndx) = index_back_atObs(rand_evalPoint(jndx))
	! 			jndx = jndx + 1
	! 		else
	! 			cycle
	! 		endif
	! 		!If (jndx >= num_Eval) exit 
	! 	Enddo 
	! 	Print*, "Finished selecting obs points"	
	! Endif	
	! Call MPI_BCAST(index_back_atEval, num_Eval, mpiInt_size, 0, MPI_COMM_WORLD, IERR)
	! Call MPI_BCAST(rand_evalPoint, num_Eval, mpiInt_size, 0, MPI_COMM_WORLD, IERR)	
	
	RETURN
	
 END SUBROUTINE Observation_Operator_Parallel
 
 SUBROUTINE Observation_Read_SNOTEL(snotel_inp_file,  &
				dim_name,			&
				NDIM, &
				SWE_SNOTEL,      &
				SNWD_SNOTEL,		&
				Lat_SNOTEL,      &
				Lon_SNOTEL,		&
                !  Ele_SNOTEL		&
				MYRANK)
	IMPLICIT NONE

	include 'mpif.h'

	!Open netCDF for a snotel and read the SWE, SnowDepth,..., Lat, Lon, at a given datetime
	!ToDO: Can you use variable length char array ?
	CHARACTER(LEN=*), Intent(In)      :: snotel_inp_file, dim_name

	INTEGER                :: ERROR, NCID
	INTEGER                :: MYRANK
	INTEGER                :: ID_DIM, ID_VAR 
	INTEGER, Intent(out)   :: NDIM

	REAL, ALLOCATABLE, Intent(Out)    :: SWE_SNOTEL(:), SNWD_SNOTEL(:)
	REAL, ALLOCATABLE, Intent(Out)	   :: Lat_SNOTEL(:), Lon_SNOTEL(:) !, Ele_SNOTEL(:)
	REAL, ALLOCATABLE	   :: SWE_Ratio(:) !, Ele_SNOTEL(:)
	Integer                :: num_invalid

	ERROR=NF90_OPEN(TRIM(snotel_inp_file),NF90_NOWRITE,NCID)
	CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(snotel_inp_file) )

	ERROR=NF90_INQ_DIMID(NCID, TRIM(dim_name), ID_DIM)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Dimension' )

	ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NDIM)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Dimension' )


	ALLOCATE(SWE_SNOTEL(NDIM))
	ALLOCATE(SNWD_SNOTEL(NDIM))
	ALLOCATE(Lat_SNOTEL(NDIM))
	ALLOCATE(Lon_SNOTEL(NDIM))
	ALLOCATE(SWE_Ratio(NDIM))	
	!ALLOCATE(Ele_SNOTEL(NDIM))

	ERROR=NF90_INQ_VARID(NCID, 'SWE', ID_VAR)
	CALL NETCDF_ERR(ERROR, 'ERROR READING SWE ID' )
	ERROR=NF90_GET_VAR(NCID, ID_VAR, SWE_SNOTEL)
	CALL NETCDF_ERR(ERROR, 'ERROR READING SWE RECORD' )

	ERROR=NF90_INQ_VARID(NCID, 'SNWD', ID_VAR)
	CALL NETCDF_ERR(ERROR, 'ERROR READING SNWD ID' )
	ERROR=NF90_GET_VAR(NCID, ID_VAR, SNWD_SNOTEL)
	CALL NETCDF_ERR(ERROR, 'ERROR READING SNWD RECORD' )

	ERROR=NF90_INQ_VARID(NCID, 'lat', ID_VAR)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Lat ID' )
	ERROR=NF90_GET_VAR(NCID, ID_VAR, Lat_SNOTEL)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Lat RECORD' )

	ERROR=NF90_INQ_VARID(NCID, 'lon', ID_VAR)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Lon ID' )
	ERROR=NF90_GET_VAR(NCID, ID_VAR, Lon_SNOTEL)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Lon RECORD' )

	SWE_Ratio = SWE_SNOTEL / SWE_SNOTEL
	num_invalid = COUNT (IEEE_IS_NAN(SWE_Ratio))
	call debug_print("number of invalid values ", float(num_invalid))

	! need to read corresponding elevation values 
	! ERROR=NF90_INQ_VARID(NCID, 'Elevation', ID_VAR)
	! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation ID' )
	! ERROR=NF90_GET_VAR(NCID, ID_VAR, Ele_SNOTEL)
	! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation RECORD' )

	ERROR = NF90_CLOSE(NCID)

	DEALLOCATE(SWE_Ratio)
			  
	RETURN
	
 End SUBROUTINE Observation_Read_SNOTEL
!
SUBROUTINE Observation_Read_GHCND_Tile(ghcnd_inp_file, dim_name,			&
				lat_min, lat_max, lon_min, lon_max, &
				NDIM, 			&
				SNWD_GHCND,		&
				Lat_GHCND,      &
				Lon_GHCND,		&
				!Ele_GHCND,		&
				MYRANK)
	
	IMPLICIT NONE

	include 'mpif.h'
	!Open netCDF and read the SWE, SnowDepth,..., Lat, Lon, at a given datetime
	!ToDO: Can you use variable length char array ?
	CHARACTER(LEN=*), Intent(In)      :: ghcnd_inp_file, dim_name
    REAL, Intent(In)       :: lat_min, lat_max, lon_min, lon_max 
	INTEGER, Intent(Out)   :: NDIM
	REAL, ALLOCATABLE, Intent(Out)    :: SNWD_GHCND(:)
	REAL, ALLOCATABLE, Intent(Out)	  :: Lat_GHCND(:), Lon_GHCND(:) !, Ele_GHCND(:)

	INTEGER                :: MYRANK, ERROR, NCID, ID_DIM, ID_VAR, NDIM_In
	REAL, ALLOCATABLE      :: SNWD_GHCND_In(:), Lat_GHCND_In(:), Lon_GHCND_In(:) !, Ele_GHCND(:)
	INTEGER, ALLOCATABLE   :: index_Array(:)
	INTEGER                :: jndx, jcounter

	ERROR=NF90_OPEN(TRIM(ghcnd_inp_file),NF90_NOWRITE,NCID)
	CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(ghcnd_inp_file) )

	ERROR=NF90_INQ_DIMID(NCID, TRIM(dim_name), ID_DIM)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Dimension' )

	ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NDIM_In)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Dimension' )

	ALLOCATE(SNWD_GHCND_In(NDIM_In))
	ALLOCATE(Lat_GHCND_In(NDIM_In))
	ALLOCATE(Lon_GHCND_In(NDIM_In))
	!ALLOCATE(Ele_GHCND(NDIM_In))

	ERROR=NF90_INQ_VARID(NCID, 'SNWD', ID_VAR)
	CALL NETCDF_ERR(ERROR, 'ERROR READING SNWD ID' )
	ERROR=NF90_GET_VAR(NCID, ID_VAR, SNWD_GHCND_In)
	CALL NETCDF_ERR(ERROR, 'ERROR READING SNWD RECORD' )

	ERROR=NF90_INQ_VARID(NCID, 'lat', ID_VAR)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Lat ID' )
	ERROR=NF90_GET_VAR(NCID, ID_VAR, Lat_GHCND_In)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Lat RECORD' )

	ERROR=NF90_INQ_VARID(NCID, 'lon', ID_VAR)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Lon ID' )
	ERROR=NF90_GET_VAR(NCID, ID_VAR, Lon_GHCND_In)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Lon RECORD' )
	! need to read corresponding elevation values 
	! ERROR=NF90_INQ_VARID(NCID, 'Elevation', ID_VAR)
	! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation ID' )
	! ERROR=NF90_GET_VAR(NCID, ID_VAR, Ele_GHCND_In)
	! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation RECORD' )
	ALLOCATE(index_Array(NDIM_In))
	NDIM = 0
	jcounter = 1
	Do jndx = 1, NDIM_In
		! Print*, "jndx = ", jndx
		If((Lat_GHCND_In(jndx) > lat_min) .and. (Lat_GHCND_In(jndx) < lat_max) .and. &
			(Lon_GHCND_In(jndx) > lon_min) .and. (Lon_GHCND_In(jndx) < lon_max)) then
				NDIM = NDIM + 1
				index_Array(jcounter) = jndx
				jcounter = jcounter + 1
		Endif
	End do
	ALLOCATE(SNWD_GHCND(NDIM))
	ALLOCATE(Lat_GHCND(NDIM))
	ALLOCATE(Lon_GHCND(NDIM))
	!ALLOCATE(Ele_GHCND(NDIM))
	If(NDIM > 0) then
		Do jndx = 1, NDIM
			SNWD_GHCND(jndx) = SNWD_GHCND_In(index_Array(jndx))
			Lat_GHCND(jndx) = Lat_GHCND_In(index_Array(jndx))
			Lon_GHCND(jndx) = Lon_GHCND_In(index_Array(jndx))
		End do
	Endif
	DEALLOCATE(index_Array)
	! jcounter = 1
	! If(NDIM > 0) then
	! 	Do jndx = 1, NDIM_In
	! 		If((Lat_GHCND_In(jndx) > lat_min) .and. (Lat_GHCND_In(jndx) < lat_max)) then
	! 			If((Lon_GHCND_In(jndx) > lon_min) .and. (Lon_GHCND_In(jndx) < lon_max)) then
	! 				SNWD_GHCND(jcounter) = SNWD_GHCND_In(jndx)
	! 				Lat_GHCND(jcounter) = Lat_GHCND_In(jndx)
	! 				Lon_GHCND(jcounter) = Lon_GHCND_In(jndx)
	! 				jcounter = jcounter + 1
	! 			Endif
	! 		Endif
	! 	End do
	! Endif	

	ERROR = NF90_CLOSE(NCID)
			  
	RETURN
	
 End SUBROUTINE Observation_Read_GHCND_Tile

 SUBROUTINE Observation_Read_GHCND(ghcnd_inp_file, dim_name,			&
				NDIM, 			&
				SNWD_GHCND,		&
				Lat_GHCND,      &
				Lon_GHCND,		&
				!Ele_GHCND,		&
				MYRANK)
	
	IMPLICIT NONE

	include 'mpif.h'
	!Open netCDF for a snotel and read the SWE, SnowDepth,..., Lat, Lon, at a given datetime
	!ToDO: Can you use variable length char array ?
	CHARACTER(LEN=*), Intent(In)      :: ghcnd_inp_file, dim_name
	INTEGER                :: ERROR, NCID
	INTEGER                :: MYRANK
	INTEGER                :: ID_DIM, ID_VAR
	INTEGER, Intent(Out)   :: NDIM

	REAL, ALLOCATABLE, Intent(Out)    :: SNWD_GHCND(:)
	REAL, ALLOCATABLE, Intent(Out)	   :: Lat_GHCND(:), Lon_GHCND(:) !, Ele_GHCND(:)

	ERROR=NF90_OPEN(TRIM(ghcnd_inp_file),NF90_NOWRITE,NCID)
	CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(ghcnd_inp_file) )

	ERROR=NF90_INQ_DIMID(NCID, TRIM(dim_name), ID_DIM)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Dimension' )

	ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NDIM)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Dimension' )

	ALLOCATE(SNWD_GHCND(NDIM))
	ALLOCATE(Lat_GHCND(NDIM))
	ALLOCATE(Lon_GHCND(NDIM))
	!ALLOCATE(Ele_GHCND(NDIM))

	ERROR=NF90_INQ_VARID(NCID, 'SNWD', ID_VAR)
	CALL NETCDF_ERR(ERROR, 'ERROR READING SNWD ID' )
	ERROR=NF90_GET_VAR(NCID, ID_VAR, SNWD_GHCND)
	CALL NETCDF_ERR(ERROR, 'ERROR READING SNWD RECORD' )

	ERROR=NF90_INQ_VARID(NCID, 'lat', ID_VAR)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Lat ID' )
	ERROR=NF90_GET_VAR(NCID, ID_VAR, Lat_GHCND)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Lat RECORD' )

	ERROR=NF90_INQ_VARID(NCID, 'lon', ID_VAR)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Lon ID' )
	ERROR=NF90_GET_VAR(NCID, ID_VAR, Lon_GHCND)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Lon RECORD' )

	! need to read corresponding elevation values 
	! ERROR=NF90_INQ_VARID(NCID, 'Elevation', ID_VAR)
	! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation ID' )
	! ERROR=NF90_GET_VAR(NCID, ID_VAR, Ele_GHCND)
	! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation RECORD' )

	ERROR = NF90_CLOSE(NCID)
			  
	RETURN
	
 End SUBROUTINE Observation_Read_GHCND

 ! the following code based on write_data() in read_write_data.f90
 Subroutine Write_DA_Outputs(output_file, idim, jdim, lensfc, myrank,   &
							 snoforc, snoanl, snwdforc, snwdanal, inovatgrid, SNCOV_IMS, &
							 num_Eval, Lat_atEvalPts, Lon_atEvalPts, Obs_atEvalPts, & 
						     SNOFCS_atEvalPts, innov_atEvalPts, SNOANL_atEvalPts)  !, anl_fSCA) !updated snocov
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
	real, intent(in)    :: SNOFCS_atEvalPts(num_Eval), SNOANL_atEvalPts(num_Eval), innov_atEvalPts(num_Eval)

	integer                     :: fsize=65536, inital=0
	integer                     :: header_buffer_val = 16384
	integer                     :: dims_3d(3), dims_strt(3), dims_end(3)
	integer                     :: error, i, ncid
	integer                     :: dim_x, dim_y, dim_time, dim_eval
	integer                     :: id_x, id_y, id_time
	integer                     :: id_swe_forc, id_swe, id_snwdf, id_snwd, id_innov, id_imscov   !, id_anlscov
	integer       :: id_lateval, id_loneval, id_obseval, id_forceval, id_anleval, id_innoveval   !, id_anlscov

	integer                     :: myrank

	real(kind=4)                :: times
	real(kind=4), allocatable   :: x_data(:), y_data(:)
	real(kind=8), allocatable   :: dum2d(:,:)

	include "mpif.h"

	print*
	print*,"Process ", myrank, "writing output data to: ",trim(output_file)

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
	error = nf90_put_att(ncid, id_forceval, "units", "deg")
	call netcdf_err(error, 'DEFINING SNOFCS_atEvalPts UNITS' )

	error = nf90_def_var(ncid, 'SNOANL_atEvalPts', NF90_DOUBLE, dim_eval, id_anleval)
	call netcdf_err(error, 'DEFINING SNOANL_atEvalPts' )
	error = nf90_put_att(ncid, id_anleval, "long_name", "Analysis at Evaluation Points")
	call netcdf_err(error, 'DEFINING SNOANL_atEvalPts LONG NAME' )
	error = nf90_put_att(ncid, id_anleval, "units", "mm")
	call netcdf_err(error, 'DEFINING SNOANL_atEvalPts UNITS' )
	
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
	
	error = nf90_put_var( ncid, id_innoveval, innov_atEvalPts)
	call netcdf_err(error, 'WRITING innov_atEvalPts RECORD' )

	deallocate(x_data, y_data)
	deallocate(dum2d)

	error = nf90_close(ncid)

 End subroutine Write_DA_Outputs

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

 ! copied from read_write_data.f90---may be better to 'include'
 SUBROUTINE debug_print(STRING, num_val )

	!--------------------------------------------------------------
	! prints ERROR  MESSAGE
	!--------------------------------------------------------------
	
	 IMPLICIT NONE
	
	 CHARACTER(LEN=*), INTENT(IN) :: STRING
	 real, Intent(in)			  :: num_val
	 CHARACTER(LEN=20)			  :: numval_Str

	 write(numval_Str, "(F18.3)"),  num_val
	 	
	 IF(print_deb) PRINT*, TRIM(STRING), " ", numval_Str
	
	 RETURN
 END SUBROUTINE debug_print

 SUBROUTINE NETCDF_ERR( ERR, STRING )

	!--------------------------------------------------------------
	! IF AT NETCDF CALL RETURNS AN ERROR, PRINT OUT A MESSAGE
	! AND STOP PROCESSING.
	!--------------------------------------------------------------
	
	 IMPLICIT NONE
	
	 include 'mpif.h'
	
	 INTEGER, INTENT(IN) :: ERR
	 CHARACTER(LEN=*), INTENT(IN) :: STRING
	 CHARACTER(LEN=80) :: ERRMSG
	
	 IF( ERR == NF90_NOERR )RETURN
	 ERRMSG = NF90_STRERROR(ERR)
	 PRINT*,''
	 PRINT*,'FATAL ERROR: ', TRIM(STRING), ': ', TRIM(ERRMSG)
	 PRINT*,'STOP.'
	 CALL MPI_ABORT(MPI_COMM_WORLD, 999)
	
	 RETURN
 END SUBROUTINE NETCDF_ERR


 END MODULE M_Snow_Analysis
 
