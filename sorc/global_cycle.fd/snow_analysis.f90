MODULE M_Snow_Analysis

USE NETCDF
USE M_DA_OI
Use, Intrinsic :: IEEE_ARITHMETIC	

Logical, parameter		  :: print_deb = .False.  ! change this to /= 1 to turn off print

CONTAINS
   
subroutine Snow_Analysis(RLA, RLO, OROG, MYRANK, IDIM, JDIM, &
							IY, IM, ID, IH, &
							TILE_NUM,         &
							LENSFC,	 		&
							SNOFCS_Inp, SWDFCS,  &
							SNWDEN,	VETFCS,		 &
							SNOANL, SWDANL, anl_fSCA)
							
	!----------------------------------------------------------------------
	! Input: forecast/background states for a single tile by a single MPI process
	! reads observations: snow depth / SWE and cover
	! calls observation operator for a given tile 
	! does OI update per tile
	! returns updated/Analysis states back to the caller, surface drive (sfcdrv) program
	!
	! RLA, RLO: lat lon information for the tile, INPUT to this subroutine
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

	INTEGER, intent(in) :: IY, IM, ID, IH, IDIM, JDIM
	INTEGER, intent(In) :: MYRANK, LENSFC
	CHARACTER(LEN=5), INTENT(In) :: TILE_NUM
	Character(LEN=3)       :: rank_str
	INTEGER			 :: IERR
	
	REAL, intent(In)    :: RLA(LENSFC), RLO(LENSFC), OROG(LENSFC)
	REAL, intent(In)    :: SNOFCS_Inp(LENSFC), SWDFCS(LENSFC)
	REAL                :: SNOFCS(LENSFC), SNWDEN(LENSFC), VETFCS(LENSFC)
	REAL, intent(Out)   :: SNOANL(LENSFC), SWDANL(LENSFC), anl_fSCA(LENSFC)

	CHARACTER(len=250)   :: snotel_inp_file, dim_name, ghcnd_inp_file, ims_inp_file, ims_inp_file_indices
	CHARACTER(len=250) 	 :: da_out_file
	CHARACTER(len=5)     :: y_str, m_str, d_Str, fvs_tile
	REAL, ALLOCATABLE    :: SWE_SNOTEL(:), SNWD_SNOTEL(:), SNWDensity_SNOTEL(:)
	REAL, ALLOCATABLE    :: Lat_SNOTEL(:), Lon_SNOTEL(:), Ele_SNOTEL(:)
	REAL, ALLOCATABLE    ::  SNOFCS_atSNOTEL(:) !, SWDFCS_atSNOTEL(:) 
	REAL, ALLOCATABLE    :: SNWD_GHCND(:), SNOFCS_atGHCND(:)		!, SWDFCS_atGHCND(:)
	REAL, ALLOCATABLE	 :: Lat_GHCND(:), Lon_GHCND(:), Ele_GHCND(:)
	Integer, ALLOCATABLE 	:: index_back_atSNOTEL(:), index_back_atGHCND(:)   ! the location of the corresponding obs
    REAL, ALLOCATABLE    ::  SNOANL_atEvalPts(:), innov_atEvalPts(:)            !analysis at evalution points 
	
	!--ims resampled to each grid REAL            	 :: lat_min, lat_max, lon_min, lon_max 
	Real				 :: SNCOV_IMS(LENSFC)  !, SNOFCS_atIMS(:), SWDFCS_atIMS(:)
	!--ims resampled to each grid Integer              :: IMS_Loc_Array(LENSFC)
	Real				 :: SNWD_IMS_at_Grid(LENSFC)
	!--ims resampled to each grid Real				   :: Lat_IMS_atGrid(LENSFC), Lon_IMS_atGrid(LENSFC)
	!--ims resampled to each grid REAL, ALLOCATABLE	   :: Lat_IMS(:), Lon_IMS(:)  !, Ele_Ims(:)

	INTEGER :: Num_Snotel, Num_Ghcnd, Num_Ims, num_sub !Num_Ims_Lat, Num_Ims_Lon
	Real	:: max_distance   ! radius_of_influence for selecting state at observation point

	INTEGER :: jndx, zndx, ncol, nrow
	Integer, Allocatable   :: Loc_backSt_atObs(:), Loc_backSt_atObs2(:)
	Integer				   :: num_loc, num_loc_1, num_loc_2, max_num_loc
	Real, Parameter		:: Stdev_back = 30., Stdev_Obs_depth = 40., Stdev_Obs_ims = 80. ! mm 
	Integer				:: ims_assm_hour
	Real				:: obs_tolerance, ims_max_ele
	Real                :: L_horz, h_ver
	Real(dp), Allocatable 			:: B_cov_mat(:,:), b_cov_vect(:)
	Real(dp), Allocatable 		    :: O_cov_mat(:,:), W_wght_vect(:)
	Real, Allocatable   :: back_at_Obs(:), obs_Array(:), Lat_Obs(:), Lon_Obs(:), Ele_Obs(:)
	REAL                :: innov_at_Grid(LENSFC), anl_at_Grid(LENSFC)
	Real, Allocatable  :: obs_Innov(:)
	Real			   :: ims_threshold  ! threshold for converting IMS fSCA to binary 1, 0 values	
	Real			   :: snwden_val, snwden_val_SNOTEL
	LOGICAL 		   :: assim_SWE, assim_SNOTEL, assim_GHCND   !assimilate swe? (instead of snow depth)
	LOGICAL 		   :: assim_IMS, assim_IMS_thisGridCell    !assimilate sncov
	
	! Real                 :: max_value
    ! INTEGER(SIZEOF_SIZE_T)   :: i_isize
    ! i_isize = SIZEOF(max_value)
	! print*, "size of real: ", i_isize
	! print*, "size of int: ", SIZEOF_SIZE_T
	
	assim_SWE = .False.  ! note: if this is set true, may need to adjust background and obs stddev above
	assim_SNOTEL = .False.
	assim_GHCND = .True.
	assim_IMS = .True.     !!assimilate sncov; 
	assim_IMS_thisGridCell = .FALSE.    ! if assimilating ims, skip this grid cell for this time step

	! integer :: A_Test(5)
	! A_Test = (/1, -6, 3, 4, 5/)

	!NUM_THREADS = NUM_PARTHDS()  ! ? no include?

	!CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
	
	!if (MYRANK /= 4 ) goto 999

	PRINT*,"STARTING Snow analysis on RANK ",MYRANK, " Tile number: ", TILE_NUM
	PRINT*
	! first call observation readers 
	! observatoin in arrays of LENOBS

	! 4.9.20 for now hard code the file location
	! ToDO: use file location the way the rest of the program accesses 
	! SNOTEL.SWE.SNWD.2020033018.nc
	! data_dir = /scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/SNOTEL/
	write(y_str, "(I4)") IY
	write(m_str, "(I0.2)") IM
	write(d_str, "(I0.2)") ID
	write(fvs_tile, "(I3)") IDIM
	snotel_inp_file = "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/SNOTEL/SNOTEL.SWE.SNWD."// &
						TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18.nc"
	dim_name = "Site_Id"
	
	!ALLOCATE(Ele_SNOTEL(Num_Snotel))
		! if (myrank==0) then 
	call Observation_Read_SNOTEL(snotel_inp_file,  &
					dim_name,			&
					Num_Snotel,      &
					SWE_SNOTEL,      &
					SNWD_SNOTEL,		&
					Lat_SNOTEL,      &
					Lon_SNOTEL,		&
					!Ele_SNOTEL,		&
					MYRANK)

	! if (myrank==4) then
	! 	PRINT*, "SNOTEL SWE from rank: ", MYRANK
	! 	PRINT*, SWE_SNOTEL	
	! endif
	PRINT*,'Finished reading SNOTEL ON RANK: ', MYRANK	
	ALLOCATE(SNWDensity_SNOTEL(Num_Snotel))
	SNWDensity_SNOTEL = SNWD_SNOTEL / SWE_SNOTEL
	snwden_val_SNOTEL = SUM(SNWDensity_SNOTEL, Mask = (.not. IEEE_IS_NAN(SNWDensity_SNOTEL))) &
	                 / COUNT (.not. IEEE_IS_NAN(SNWDensity_SNOTEL))
	print*, "Process ", MYRANK, " SNOTEL average swe/snwd ratio: ", snwden_val_SNOTEL
	DEALLOCATE(SNWDensity_SNOTEL)

	! snow depth from GHCND 
	! GHCND.SNWD.2019100118.nc
	ghcnd_inp_file = "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/GHCND/GHCND.SNWD."// &
						TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18.nc"
	!dim_name = "Site_Id"
	! ALLOCATE(Ele_GHCND(Num_Ghcnd))
	call Observation_Read_GHCND(ghcnd_inp_file,  &
					dim_name,			&
					Num_Ghcnd, 		&
					SNWD_GHCND,		&
					Lat_GHCND,      &
					Lon_GHCND,		&
					!Ele_GHCND,		&
					MYRANK)

	!	PRINT*, "GHCND SNWD from rank: ", MYRANK
	!	PRINT*, SNWD_GHCND
	PRINT*,'Finished reading GHCND ON RANK: ', MYRANK	

	! IMS.SNCOV.2020022318.nc
	! ims_inp_file = "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS/IMS.SNCOV."// &
	! 					TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18.nc"
	! lat_min = MINVAL(RLA) - 1.; lat_max = MAXVAL(RLA) + 1.;
	! lon_min = MINVAL(RLO) - 1.; lon_max = MAXVAL(RLO) + 1.;	
	! if (lon_min > 180) lon_min = lon_min - 360
	! if (lon_max > 180) lon_max = lon_max - 360
	! print*, "min/max lat ", lat_min,  lat_max
	! print*, "min/max lon ", lon_min,  lon_max
	! call Observation_Read_IMS(ims_inp_file, MYRANK,  &
	! 					lat_min, lat_max, lon_min, lon_max, &
	! 					Num_Ims,       & !_lat, DIM_LEN_lon,  &
	! 					SNCOV_IMS, Lat_IMS, Lon_IMS)
	! ims_inp_file = "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS/IMS.SNCOV."// &
	! 					TRIM(y_str)//TRIM(m_str)//TRIM(d_str)// &
	! 					'.180000.sfcanl_data.'//TRIM(TILE_NUM)//".nc"
	! Call Observation_Read_IMS_tile(ims_inp_file, MYRANK, JDIM, IDIM, &
	! 								  SNCOV_IMS)	
	! if (myrank==4) then
	! 	! Print*, "RLA"
	! 	! Print*, RLA
	! 	PRINT*, "IMS SNCOV from rank: ", MYRANK
	! 	PRINT*, SNCOV_IMS
	! endif ! if myrnak==4

	num_sub = 30   !627  ! (max) number of ims subcells within a tile grid cell
	ims_inp_file = "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS/IMS.SNCOV."// &
						TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18.nc"                      !
	ims_inp_file_indices = "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS/C"// &
						TRIM(ADJUSTL(fvs_tile))//".IMS.Indices."//TRIM(TILE_NUM)//".nc"			
	Call Observation_Read_IMS_Full(ims_inp_file, ims_inp_file_indices, &
				MYRANK, JDIM, IDIM, num_sub, &
				SNCOV_IMS)
	! if (myrank==4) then
	! 	! Print*, "RLA"
	! 	! Print*, RLA
	! 	PRINT*, "IMS SNCOV from rank: ", MYRANK
	! 	PRINT*, SNCOV_IMS
	! endif ! if myrnak==4
	PRINT*,'Finished reading IMS ON RANK: ', MYRANK
	! goto 999
	
	! broadcast observed arrays to the other processes
	!  CALL MPI_BCAST(SWE_SNOTEL, Num_Snotel, MPI_REAL, 0, MPI_COMM_WORLD, IERR)
	!  CALL MPI_BCAST(SNWD_SNOTEL, Num_Snotel, MPI_REAL, 0, MPI_COMM_WORLD, IERR)
	!  CALL MPI_BCAST(SNWD_GHCND, Num_Ghcnd, MPI_REAL, 0, MPI_COMM_WORLD, IERR)
	!  CALL MPI_BCAST(SNCOV_IMS, Num_Ims_Lat * Num_Ims_Lon, MPI_INT, 0, MPI_COMM_WORLD, IERR)	
	!CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
	
	! 4.8.20: compute snow density from forecast swe and snwdepth
	SNWDEN = SWDFCS / SNOFCS_Inp
	snwden_val = SUM(SNWDEN, Mask = (.not. IEEE_IS_NAN(SNWDEN))) &
	                 / COUNT (.not. IEEE_IS_NAN(SNWDEN))
	SNWDEN = snwden_val  !10.   ! ratio of sndepth to swe when swe /=0 : SWDFCS / SNOFCS
	print*, "Process ", MYRANK, " average swe/snwd ratio: ", snwden_val
	if (assim_SWE) then 
		SNWD_SNOTEL = SWE_SNOTEL
		SNWD_GHCND = SNWD_GHCND / SNWDEN		!snwden_val
		SNOFCS = SNOFCS_Inp
	else
		SNOFCS = SNOFCS_Inp * SNWDEN    !snwden_val  !SWDFCS			! assimilate depth
		! SNWD_SNOTEL = SNWD_SNOTEL
		! SNWD_GHCND = SNWD_GHCND 
	endif
	! Get model states at obs points
	ALLOCATE(SNOFCS_atSNOTEL(Num_Snotel))
	!ALLOCATE(SWDFCS_atSNOTEL(Num_Snotel))
	ALLOCATE(Ele_SNOTEL(Num_Snotel))
	ALLOCATE(SNOFCS_atGHCND(Num_Ghcnd))
	!ALLOCATE(SWDFCS_atGHCND(Num_Ghcnd)) 
	ALLOCATE(Ele_GHCND(Num_Ghcnd)) 
	ALLOCATE(index_back_atSNOTEL(Num_Snotel))
	ALLOCATE(index_back_atGHCND(Num_Ghcnd)) 
	! 6.22.20 for evluation at snotel points
	ALLOCATE(SNOANL_atEvalPts(Num_Snotel)) 
	ALLOCATE(innov_atEvalPts(Num_Snotel))

	!ALLOCATE(SNOFCS_atIMS(Num_Ims)) !_Lat * Num_Ims_Lon))
	!ALLOCATE(SWDFCS_atIMS(Num_Ims)) !_Lat * Num_Ims_Lon))
	max_distance = 27.  !240.			!Km 
	call Observation_Operator(RLA, RLO, OROG, Lat_SNOTEL, Lon_SNOTEL,   &
						LENSFC, Num_Snotel, max_distance, 		&
						SNOFCS, 				&
						SNOFCS_atSNOTEL, Ele_SNOTEL, index_back_atSNOTEL)
							!,                 &						intp_mode)
	!CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
	if (myrank==4) then
		! PRINT*, "Background SWE  from rank: ", MYRANK
		! PRINT*, SNOFCS	
		! PRINT*
		PRINT*, "Latitude  from rank: ", MYRANK, " range: ", MINVAL(RLA), " ", MAXVAL(RLA)
		!PRINT*, RLA	
		PRINT*
		PRINT*, "Longitude from rank: ", MYRANK, " range: ", MINVAL(RLO), " ", MAXVAL(RLO)
		!PRINT*, RLO	
		PRINT*
		PRINT*, "Latitude at SNOTEL locations from rank: ", MYRANK, " range: ", &
										MINVAL(Lat_SNOTEL), " ", MAXVAL(Lat_SNOTEL)
		!PRINT*, Lat_SNOTEL	
		PRINT*
		PRINT*, "Longitude at SNOTEL locations from rank: ", MYRANK, " range: ", &
										MINVAL(Lon_SNOTEL), " ", MAXVAL(Lon_SNOTEL)
		!PRINT*, Lon_SNOTEL	
		! PRINT*
		! PRINT*, "Background SWE at SNOTEL locations from rank: ", MYRANK
		! PRINT*, SNOFCS_atSNOTEL	
		! PRINT*
		! PRINT*, "OROG at model grid cells rank: ", MYRANK
		! PRINT*, OROG
		! PRINT*
		! PRINT*, "Elevation at SNOTEL locations from rank: ", MYRANK
		! PRINT*, Ele_SNOTEL
	endif
	!CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

	PRINT*,'Finished observation operator for SNOTEL ON RANK: ', MYRANK
	
	max_distance = 27.  !120.			!Km 
	call Observation_Operator(RLA, RLO, OROG, Lat_GHCND, Lon_GHCND,   &
						LENSFC, Num_Ghcnd, max_distance,            &
						SNOFCS,                                  &
						SNOFCS_atGHCND, Ele_GHCND, index_back_atGHCND)
	! 						!,                 &						intp_mode) 
	! CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
	if (myrank==4) then
		PRINT*, "Latitude at GHCND locations from rank: ", MYRANK, " range: ", &
										MINVAL(Lat_GHCND), " ", MAXVAL(Lat_GHCND)
		!PRINT*, Lat_SNOTEL	
		PRINT*
		PRINT*, "Longitude at GHCND locations from rank: ", MYRANK, " range: ", &
										MINVAL(Lon_GHCND), " ", MAXVAL(Lon_GHCND)
		!PRINT*, Lon_SNOTEL	
		! PRINT*
		! PRINT*, "Background SWE at GHCND locations from rank: ", MYRANK
		! PRINT*, SNOFCS_atGHCND	
		! PRINT*
		! PRINT*, "Elevation at GHCND locations from rank: ", MYRANK
		! PRINT*, Ele_GHCND
	endif
	! CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
	PRINT*,'Finished observation operator for GHCND ON RANK: ', MYRANK	

	! max_distance = 60.			!Km , ims_threshold,
	! call Observation_Operator_IMS(RLA, RLO, Lat_IMS, Lon_IMS,   &
	! 						SNCOV_IMS, SNOFCS, SNWDEN, assim_SWE,	   &
	! 						LENSFC, Num_Ims, max_distance, 		&
	! 						IMS_Loc_Array, SNWD_IMS_at_Grid, Lat_IMS_atGrid, Lon_IMS_atGrid) 

	ims_threshold = 0.5  ! threshold for converting IMS fSCA to binary 1, 0 values
	Call Observation_Operator_IMS_fSCA(SNCOV_IMS, SNOFCS, SNWDEN, assim_SWE,	   &
							LENSFC, ims_threshold, 		&
							SNWD_IMS_at_Grid) !,      &
	if (myrank==4) then
		! PRINT*, "Latitude at IMS locations from rank: ", MYRANK, " range: ", &
		! 								MINVAL(Lat_IMS), " ", MAXVAL(Lat_IMS)
		! !PRINT*, Lat_IMS	
		! PRINT*
		! PRINT*, "Longitude at IMS locations from rank: ", MYRANK, " range: ", &
		! 								MINVAL(Lon_IMS), " ", MAXVAL(Lon_IMS)
		!PRINT*, Lon_IMS	
		! PRINT*
		! PRINT*, "IMS locations corresponding to each grid cell from rank: ", MYRANK
		! PRINT*, IMS_Loc_Array	
		! PRINT*
		PRINT*, "IMS obs at each grid cell from rank: ", MYRANK
		!PRINT*, SNWD_IMS_at_Grid
	endif
	! CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
	PRINT*,'Finished observation operator for IMS ON RANK: ', MYRANK

	! ens loc: /scratch2/BMC/gsienkf/Tseganeh.Gichamo/ROTDIRS/dec15/enkfgdas.20191215/18/mem001/RESTART

	! 4.13.20 use IEEE_IS_NAN() in assimilation DA to exclude points too far from obs
	
	L_horz = 55.  !120.  !
	h_ver = 800.  !1200.	!
	obs_tolerance = 5.0
	max_num_loc = 50   !100	!
	ims_max_ele = 1500.
	ims_assm_hour = 18
	!if (myrank==4) then ! for now do da for tile 5
	!Do ncol=1, IDIM
	!Do nrow = 1, Num_Snotel !JDIM/2
	    Do jndx = 1, LENSFC		!jndx = index_back_atSNOTEL(nrow)   ! 	jndx = (ncol-1)*JDIM + nrow		! 
			num_loc_1 = 0
			num_loc_2 = 0
			assim_IMS_thisGridCell = .FALSE.
			call debug_print("loop ", float(jndx))
			max_distance = 250.   !240.			!Km 
			if(assim_SNOTEL) then  
				call nearest_Observations_Locations(RLA(jndx), RLO(jndx),    &
						Lat_SNOTEL, Lon_SNOTEL, Num_Snotel, max_distance, max_num_loc,    &
						Stdev_back, Stdev_Obs_depth, obs_tolerance,                 &
						SNOFCS_atSNOTEL, SNWD_SNOTEL,						 &
						Loc_backSt_atObs,  num_loc_1) !,      &LENSFC,
				call debug_print("number of SNOTEL sndpth obs ", float(num_loc_1))
			endif

			! GHCND
			max_distance = 250.   !120.			!Km 
			if(assim_GHCND) then
				call nearest_Observations_Locations(RLA(jndx), RLO(jndx),    &
						Lat_GHCND, Lon_GHCND,  Num_Ghcnd, max_distance, max_num_loc,   &
						Stdev_back, Stdev_Obs_depth, obs_tolerance,                 &
						SNOFCS_atGHCND, SNWD_GHCND,						 &
						Loc_backSt_atObs2,  num_loc_2) !,      &LENSFC,
				call debug_print("number of GHCND sndpth obs ", float(num_loc_2))
			endif

			num_loc = num_loc_1 + num_loc_2

			!check IMS 
			if(assim_IMS) then
				if((.NOT. IEEE_IS_NAN(SNWD_IMS_at_Grid(jndx))) .AND. &
				   (OROG(jndx) <= ims_max_ele) .AND. &
				   (IH == ims_assm_hour)) then
					num_loc = num_loc_1 + num_loc_2 + 1
					assim_IMS_thisGridCell = .TRUE.
				endif
			endif
            
			if(num_loc > 0) then     !.and. (SNCOV_IMS(jndx) > 0.)) then    
				! get background states
				Allocate(back_at_Obs(num_loc))
				Allocate(obs_Array(num_loc))
				Allocate(Lat_Obs(num_loc))
				Allocate(Lon_Obs(num_loc))
				Allocate(Ele_Obs(num_loc))
				if(num_loc_1 > 0) then
					Do zndx = 1, num_loc_1    ! is assim_IMS=false>>num_loc_1+num_loc2=num_loc
						back_at_Obs(zndx) = SNOFCS_atSNOTEL(Loc_backSt_atObs(zndx))
						obs_Array(zndx) = SNWD_SNOTEL(Loc_backSt_atObs(zndx))
						Lat_Obs(zndx) = Lat_SNOTEL(Loc_backSt_atObs(zndx))
						Lon_Obs(zndx) = Lon_SNOTEL(Loc_backSt_atObs(zndx))
						Ele_Obs(zndx) = Ele_SNOTEL(Loc_backSt_atObs(zndx))
					End Do
				End if
		!TODO: what if there is no num_loc_2?????
				! ghcnd
				if(num_loc_2 > 0) then
					Do zndx = num_loc_1 + 1, num_loc_1 + num_loc_2     ! is assim_IMS=false>>num_loc_1+num_loc2=num_loc
						back_at_Obs(zndx) = SNOFCS_atGHCND(Loc_backSt_atObs2(zndx-num_loc_1))
						obs_Array(zndx) = SNWD_GHCND(Loc_backSt_atObs2(zndx-num_loc_1))
						Lat_Obs(zndx) = Lat_GHCND(Loc_backSt_atObs2(zndx-num_loc_1))
						Lon_Obs(zndx) = Lon_GHCND(Loc_backSt_atObs2(zndx-num_loc_1))
						Ele_Obs(zndx) = Ele_GHCND(Loc_backSt_atObs2(zndx-num_loc_1))
					End Do
				End if
				!ims obs
				if (assim_IMS_thisGridCell) then
					back_at_Obs(num_loc) = SNOFCS(jndx)
					obs_Array(num_loc) = SNWD_IMS_at_Grid(jndx)
					Lat_Obs(num_loc) = RLA(jndx) !Lat_IMS_atGrid(jndx)
					Lon_Obs(num_loc) = RLO(jndx)  !Lon_IMS_atGrid(jndx)
					Ele_Obs(num_loc) = OROG(jndx)
				endif
				! compute covariances
				Allocate(B_cov_mat(num_loc, num_loc))
				Allocate(b_cov_vect(num_loc))
				Allocate(O_cov_mat(num_loc, num_loc))
				Allocate(W_wght_vect(num_loc))
				call compute_covariances(RLA(jndx), RLO(jndx), OROG(jndx), SNOFCS(jndx),    &
					Lat_Obs, Lon_Obs, Ele_Obs, num_loc,    			 &
					Stdev_back, Stdev_Obs_depth, Stdev_Obs_ims,         &
					L_horz, h_ver,                                      &   !L_horz in Km, h_ver in m
					assim_IMS_thisGridCell,                                          &
					B_cov_mat, b_cov_vect, O_cov_mat, W_wght_vect)
				! call OI DA
				Allocate(obs_Innov(num_loc))
				call Snow_DA_OI(back_at_Obs, obs_Array, num_loc,  &
					W_wght_vect,            &
					SNOFCS(jndx), innov_at_Grid(jndx), anl_at_Grid(jndx), obs_Innov)
				if (myrank==4) then	
					print*, "loop ", jndx, "  num snotel obs ", num_loc_1, "num ghcnd obs ", num_loc_2, &
					        "total obs", num_loc
					PRINT*, " background at obs pts: "
					PRINT*, back_at_Obs	
					PRINT*, "Observed"
					PRINT*,  obs_Array
					PRINT*, "Obs innovation: "
					PRINT*, obs_Innov
					PRINT*, "Weight vector: "
					PRINT*, W_wght_vect	
					PRINT*, "innovation at grid:"
					PRINT*,  innov_at_Grid(jndx)
					PRINT*, "forecast:"
					PRINT*, SNOFCS(jndx)
					PRINT*, "analyis at grid: "
					PRINT*, anl_at_Grid(jndx)
				endif		
				!free mem
				DEALLOCATE(back_at_Obs, obs_Array)
				DEALLOCATE(Lat_Obs, Lon_Obs, Ele_Obs, obs_Innov)
				DEALLOCATE(B_cov_mat, b_cov_vect, O_cov_mat, W_wght_vect)
				! QCC by ims--use ims as snow mask
				if((SNCOV_IMS(jndx) >= 0.5) .and. (anl_at_Grid(jndx) < 50.)) then
					anl_at_Grid(jndx) = 50.
				elseif((SNCOV_IMS(jndx) < 0.1) .and. (anl_at_Grid(jndx) >= 0.)) then
						anl_at_Grid(jndx) = 0.
				endif
				if (myrank==4) then						
					PRINT*, "analyis at grid after mask: "
					PRINT*, anl_at_Grid(jndx)
				endif
				! if((SNCOV_IMS(jndx) >= 0.5) .and. (anl_at_Grid(jndx) < 50.)) then
				! 	anl_at_Grid(jndx) = 50.
				! elseif((SNCOV_IMS(jndx) < 0.5) .and. (SNCOV_IMS(jndx) > 0.1) .and. (anl_at_Grid(jndx) >= 50)) then
				! 	anl_at_Grid(jndx) = 50.
				! elseif((SNCOV_IMS(jndx) < 0.1) .and. (anl_at_Grid(jndx) >= 0.)) then
				! 		anl_at_Grid(jndx) = 0.
				! endif
			else
				anl_at_Grid(jndx) = SNOFCS(jndx)
			endif
			
			if(assim_SNOTEL) then  
				Deallocate(Loc_backSt_atObs) 
			endif
			if(assim_GHCND) then 
				Deallocate(Loc_backSt_atObs2) 
			endif

		End do
	!End do

		! avoid -ve anl
		Where(anl_at_Grid < 0.) anl_at_Grid = 0.
		!print*, "I am here 3, ", anl_at_Grid

		! swe and snwd
		if(assim_SWE) then
			SNOANL = anl_at_Grid 
			SWDANL = anl_at_Grid * SNWDEN		!snwden_val   ! snwden_val_SNOTEL !  
		else
			SNOANL = anl_at_Grid / SNWDEN		!snwden_val !  snwden_val_SNOTEL  ! 
			SWDANL = anl_at_Grid
			!innov_at_Grid = innov_at_Grid / SNWDEN
			!SNOFCS_atSNOTEL = SNOFCS_atSNOTEL / snwden_val
		endif
		
		! PRINT*
		! PRINT*, "Innovation SWE  from rank: ", MYRANK
		! PRINT*, innov_at_Grid	
		! PRINT*
		! PRINT*, "Analysis SWE  from rank: ", MYRANK
		! PRINT*, anl_at_Grid	
		! PRINT*
	!end if

	!Compute updated snocov	
	!Call update_snow_cover_fraction(LENSFC, SNOANL, VETFCS, anl_fSCA)

	! copy values at eval points
	SNOANL_atEvalPts = IEEE_VALUE(SNOANL_atEvalPts, IEEE_QUIET_NAN)
	Do jndx = 1, Num_Snotel
		if (index_back_atSNOTEL(jndx) > 0) then 
			innov_atEvalPts(jndx) = innov_at_Grid(index_back_atSNOTEL(jndx))
			if(assim_SWE) then
				SNOANL_atEvalPts(jndx) = SNOANL(index_back_atSNOTEL(jndx))
			else
				SNOANL_atEvalPts(jndx) = SWDANL(index_back_atSNOTEL(jndx))
			endif
		endif
	End do

	! write outputs
!ToDO: Standard output locations
	
	Write(rank_str, '(I3.3)') (MYRANK+1)
	da_out_file = "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/Analysis/SNOANL."// &  !
						TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18_tile"//rank_str//".nc"  !
	
	call Write_DA_Outputs(da_out_file, IDIM, JDIM, LENSFC, MYRANK, SNOFCS, SNOANL, SWDANL, &
						  innov_at_Grid, SNCOV_IMS, &
						  Num_Snotel, Lat_SNOTEL, Lon_SNOTEL, SNWD_SNOTEL, & 
						  SNOFCS_atSNOTEL, SNOANL_atEvalPts, innov_atEvalPts)  !, anl_fSCA) !updated snocov

	DEALLOCATE(SWE_SNOTEL, SNWD_SNOTEL, Lat_SNOTEL, Lon_SNOTEL, Ele_SNOTEL)
	DEALLOCATE(SNOFCS_atSNOTEL)   !, SWDFCS_atSNOTEL)
	DEALLOCATE(SNWD_GHCND, SNOFCS_atGHCND)  !, SWDFCS_atGHCND) 
	DEALLOCATE(Lat_GHCND, Lon_GHCND)
	DEALLOCATE(index_back_atSNOTEL, index_back_atGHCND)
	DEALLOCATE(SNOANL_atEvalPts, innov_atEvalPts)
	!DEALLOCATE(SNCOV_IMS, Lat_IMS, Lon_IMS) !SNOFCS_atIMS, SWDFCS_atIMS)
	!DEALLOCATE(Ele_SNOTEL, Ele_GHCND)
999 CONTINUE
	CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
	PRINT*,'Finished OI DA ON RANK: ', MYRANK

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
 
 ! Get model states at obs points
 ! Warning: This assumes all distance coordinates are valid; 
 ! do quality control of coordinates beforehand
 SUBROUTINE Observation_Operator(RLA, RLO, OROG, Lat_Obs, Lon_Obs,   &
						LENSFC, num_Obs, max_distance, 		&
						SNWD_back,  				&
						SNWD_atObs, Ele_atObs, index_back_atObs) !,      &
						!intp_mode) 
						!  SWE_atObs  SWE_back

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
	Do zndx = 1, num_Obs 
	    !Lon_Obs[Lon_Obs<0]= 360.0 + Lon_Obs[Lon_Obs<0]
		if (Lon_Obs(zndx) < 0) Lon_Obs_2(zndx) = 360. + Lon_Obs(zndx)
	end do
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
		! Do jndx = 1, LENSFC 
		! 	if (haversinArr(jndx) > 1) haversinArr(jndx) = 1 
		! end do
		! Do jndx = 1, LENSFC 
		! 	if (haversinArr(jndx) < 0) haversinArr(jndx) = 0 ! ensure <0
		! end do
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

 ! for each model grid cell, select the corresponding ims grid
 SUBROUTINE Observation_Operator_IMS(RLA, RLO, Lat_Obs, Lon_Obs,   &
							SNCOV_IMS, SNOFCS, SNWDEN, assim_SWE,	   &
							LENSFC, num_Obs, max_distance, 		&
							IMS_Loc_Array, SNWD_IMS_at_Grid, Lat_IMS_atGrid, Lon_IMS_atGrid) !,      &
		!intp_mode) 
		!  SWE_atObs  SWE_back

	!USE intrinsic::ieee_arithmetic
	Real, Intent(In) 	:: RLA(LENSFC), RLO(LENSFC)		!, OROG(LENSFC)
	Real, Intent(In) 	:: Lat_Obs(num_Obs), Lon_Obs(num_Obs)  ! don't want to alter these
	INTEGER, Intent(In) 	:: SNCOV_IMS(num_Obs)
	Real, Intent(In) 	:: SNOFCS(LENSFC), SNWDEN(LENSFC)
	Logical				:: assim_SWE
	INTEGER :: LENSFC, num_Obs
	Real	:: max_distance   ! radius_of_influence

	INTEGER, Intent(Out) 	:: IMS_Loc_Array(LENSFC)
	Real, Intent(Out) 	:: SNWD_IMS_at_Grid(LENSFC), Lat_IMS_atGrid(LENSFC), Lon_IMS_atGrid(LENSFC)
	
	Real 	:: Lon_Obs_2(num_Obs)		!RLO_2(LENSFC), 	
	Real 	:: RLA_rad(LENSFC), RLO_rad(LENSFC)
	Real 	:: Lat_Obs_rad(num_Obs), Lon_Obs_rad(num_Obs)	
	INTEGER :: indx, jndx, zndx, min_indx
	Real 	:: d_latArr(num_Obs), d_lonArr(num_Obs)
	Real    :: distArr(num_Obs), haversinArr(num_Obs)
	Real(16), Parameter :: PI_16 = 4 * atan (1.0_16)	
	Real(16), Parameter :: pi_div_180 = PI_16/180.0
	Real, Parameter		:: earth_rad = 6371.
	Real, Parameter		:: SWE_Tolerance = 0.000001    ! smallest obs swe value
	! PRINT*, "PI: ", PI_16
	! PRINT*, "PI / 180: ", pi_div_180
	
	IMS_Loc_Array = 0
	!Fill background values to nan (to differentiate those htat don't have value)
	SNWD_IMS_at_Grid = IEEE_VALUE(SNWD_IMS_at_Grid, IEEE_QUIET_NAN)
	Lat_IMS_atGrid = IEEE_VALUE(Lat_IMS_atGrid, IEEE_QUIET_NAN)
	Lon_IMS_atGrid = IEEE_VALUE(Lon_IMS_atGrid, IEEE_QUIET_NAN)	
	
    call debug_print("Here ", 1.)
	! RLO from 0 to 360 (no -ve lon)
	! Do zndx = 1, num_Obs 
	! 	if (Lon_Obs(zndx) < 0) Lon_Obs_2(zndx) = 360. + Lon_Obs(zndx)
	! end do
	Lon_Obs_2 = Lon_Obs
	WHERE(Lon_Obs < 0) Lon_Obs_2 = 360. + Lon_Obs
	call debug_print("Here ", 2.)
	! at each obs point compute its distance from RLA/RLO pairs 
	! then find the position of the minimum

	! shortest distance over sphere using great circle distance 	
	RLA_rad =  pi_div_180 * RLA
	RLO_rad =  pi_div_180 * RLO
	Lat_Obs_rad =  pi_div_180 * Lat_Obs
	Lon_Obs_rad =  pi_div_180 * Lon_Obs_2  
    call debug_print("Here ", 3.)
	! https://en.wikipedia.org/wiki/Haversine_formula
	! https://www.geeksforgeeks.org/program-distance-two-points-earth/
	! Distance, d = R * arccos[(sin(lat1) * sin(lat2)) + cos(lat1) * cos(lat2) * cos(long2 – long1)]
	! dist = 2 * R * asin { sqrt [sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2]}
	Do indx = 1, LENSFC  
		!print*,"loop ", indx
		d_latArr = (Lat_Obs_rad - RLA_rad(indx)) / 2.
		d_lonArr = (Lon_Obs_rad - RLO_rad(indx)) / 2.
		haversinArr = sin(d_latArr)**2 + cos(Lat_Obs_rad) * cos(RLA_rad(indx)) * sin(d_lonArr)**2

		WHERE(haversinArr > 1) haversinArr = 1. ! ensure numerical errors don't make h>1
		! Do jndx = 1, num_Obs 
		! 	if (haversinArr(jndx) > 1) haversinArr(jndx) = 1 ! ensure numerical errors don't make h>1
		! end do
		! Do jndx = 1, num_Obs
		! 	if (haversinArr(jndx) < 0) haversinArr(jndx) = 0 ! ensure <0
		! end do
		distArr = 2 * earth_rad * asin(sqrt(haversinArr))		
		!distArr = (Lat_Obs(indx) - RLA)**2 + (Lon_Obs_2(indx) - RLO)**2 
		min_indx = MINLOC(distArr, dim = 1)  !, MASK=ieee_is_nan(distArr))
		
		if(distArr(min_indx) < max_distance) then
			IMS_Loc_Array(indx) = min_indx 							
			Lat_IMS_atGrid(indx) = Lat_Obs(min_indx)
			Lon_IMS_atGrid(indx) = Lon_Obs(min_indx)
			! if land but no snow ims_snow_depth=0
			! if land covered by snow and model snow=0, ims_snow_depth = 50 mm
			! all other cases NAN
			if(SNCOV_IMS(min_indx) == 2) SNWD_IMS_at_Grid(indx) = 0. 	! land no-snow 
			if((SNCOV_IMS(min_indx) == 4) .AND. (SNOFCS(indx) < SWE_Tolerance)) SNWD_IMS_at_Grid(indx) = 50.	!SNCOV_IMS(min_indx)
		! else
		!  	Print*, " Warning distance greater than ",max_distance," km ", distArr(min_indx)
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

END SUBROUTINE Observation_Operator_IMS

SUBROUTINE Observation_Operator_IMS_fSCA(SNCOV_IMS, SNOFCS, SNWDEN, assim_SWE,	   &
							LENSFC, ims_threshold, 		&
							SNWD_IMS_at_Grid) !,      &

	Real, Intent(In) 	:: SNCOV_IMS(LENSFC), SNOFCS(LENSFC), SNWDEN(LENSFC)
	Logical				:: assim_SWE
	INTEGER 			:: LENSFC
	Real, Intent(In)    :: ims_threshold

	Real, Intent(Out) 	:: SNWD_IMS_at_Grid(LENSFC)
	
	INTEGER :: indx
	Real, Parameter		:: SWE_Tolerance = 0.001    ! smallest swe value

	!Fill background values to nan (to differentiate those that don't have value)
	SNWD_IMS_at_Grid = IEEE_VALUE(SNWD_IMS_at_Grid, IEEE_QUIET_NAN)

    call debug_print("Here ", 1.)

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

END SUBROUTINE Observation_Operator_IMS_fSCA

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

	! if (maxlat_indx > minlat_indx) then
	! 	DIM_LEN_lat = 1 + maxlat_indx - minlat_indx
	! else
	! 	DIM_LEN_lat = 1 + minlat_indx - maxlat_indx
	! if (maxlon_indx > minlon_indx) then
	! 	DIM_LEN_lon = 1 + maxlon_indx - minlon_indx
	! else
	! 	DIM_LEN_lon = 1 + minlon_indx - maxlon_indx

	! DIM_LEN = DIM_LEN_lat * DIM_LEN_lon
	! ALLOCATE(SNCOV_IMS(DIM_LEN))
	! ALLOCATE(Lat_IMS(DIM_LEN))
	! ALLOCATE(Lon_IMS(DIM_LEN))
	! ALLOCATE(LatLon_IMS_2D(DIM_LEN_lon, DIM_LEN_lat))
	! ALLOCATE(SNCOV_IMS_2D(DIM_LEN_lon, DIM_LEN_lat))	
	! !Do jndx=1, DIM_LEN_lat
	! 	Do indx=1, DIM_LEN_lon
	! 		LatLon_IMS_2D(indx, :) = Lat_IMS_1D(minlat_indx:maxlat_indx)
	! 	End do
	! !End do
	! Lat_IMS = Reshape(LatLon_IMS_2D, (/DIM_LEN/))
	! Do jndx=1, DIM_LEN_lat
	! 	!Do indx=1, DIM_LEN_lon
	! 		LatLon_IMS_2D(:, jndx) = Lon_IMS_1D(minlon_indx:maxlon_indx)
	! 	!End do
	! End do
	! Lon_IMS = Reshape(LatLon_IMS_2D, (/DIM_LEN/))

	! SNCOV_IMS_2D = SNCOV_IMS_2D_full(minlon_indx:maxlon_indx, minlat_indx:maxlat_indx)
	! SNCOV_IMS = Reshape(SNCOV_IMS_2D, (/DIM_LEN/))
			  
	RETURN
	
 End SUBROUTINE Observation_Read_IMS
 
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
SUBROUTINE Observation_Read_GHCND(ghcnd_inp_file,  &
				dim_name,			&
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
							 snoforc, snoanl, snwdanal, inovatgrid, SNCOV_IMS, &
							 Num_Snotel, Lat_SNOTEL, Lon_SNOTEL, Obs_SNOTEL, &
						     SNOFCS_atSNOTEL, SNOANL_atEvalPts, innov_atEvalPts) !, anl_fSCA)

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
	real, intent(in)            :: snoforc(lensfc), snoanl(lensfc), snwdanal(lensfc)
	Real, intent(in)            :: inovatgrid(lensfc), SNCOV_IMS(lensfc)  !, anl_fSCA(lensfc)
	integer, intent(in)         :: Num_Snotel
	real, intent(in)    :: Lat_SNOTEL(Num_Snotel), Lon_SNOTEL(Num_Snotel), Obs_SNOTEL(Num_Snotel)
	real, intent(in)    :: SNOFCS_atSNOTEL(Num_Snotel), SNOANL_atEvalPts(Num_Snotel), innov_atEvalPts(Num_Snotel)

	integer                     :: fsize=65536, inital=0
	integer                     :: header_buffer_val = 16384
	integer                     :: dims_3d(3), dims_strt(3), dims_end(3)
	integer                     :: error, i, ncid
	integer                     :: dim_x, dim_y, dim_time, dim_eval
	integer                     :: id_x, id_y, id_time
	integer                     :: id_swe_forc, id_swe, id_snwd, id_innov, id_imscov   !, id_anlscov
	integer       :: id_latsnotel, id_lonsnotel, id_obssnotel, id_forcsnotel, id_anlsnotel, id_innovsnotel   !, id_anlscov

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
	error = nf90_def_dim(ncid, 'eval_points', Num_Snotel, dim_eval)
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

	error = nf90_def_var(ncid, 'SNWD_Analysis', NF90_DOUBLE, dims_3d, id_snwd)
	call netcdf_err(error, 'DEFINING SNWDPH' )
	error = nf90_put_att(ncid, id_snwd, "long_name", "Analysis Snow Depth")
	call netcdf_err(error, 'DEFINING SNWDPH LONG NAME' )
	error = nf90_put_att(ncid, id_snwd, "units", "mm")
	call netcdf_err(error, 'DEFINING SNWDPH UNITS' )

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
	error = nf90_def_var(ncid, 'LatSNOTEL', NF90_DOUBLE, dim_eval, id_latsnotel)
	call netcdf_err(error, 'DEFINING LatSNOTEL' )
	error = nf90_put_att(ncid, id_latsnotel, "long_name", "Latitude at SNOTEL Station")
	call netcdf_err(error, 'DEFINING LatSNOTEL LONG NAME' )
	error = nf90_put_att(ncid, id_latsnotel, "units", "deg")
	call netcdf_err(error, 'DEFINING LatSNOTEL UNITS' )

	error = nf90_def_var(ncid, 'LonSNOTEL', NF90_DOUBLE, dim_eval, id_lonsnotel)
	call netcdf_err(error, 'DEFINING LonSNOTEL' )
	error = nf90_put_att(ncid, id_lonsnotel, "long_name", "Longitude at SNOTEL Station")
	call netcdf_err(error, 'DEFINING LonSNOTEL LONG NAME' )
	error = nf90_put_att(ncid, id_lonsnotel, "units", "deg")
	call netcdf_err(error, 'DEFINING LonSNOTEL UNITS' )
	
	error = nf90_def_var(ncid, 'Obs_SNOTEL', NF90_DOUBLE, dim_eval, id_obssnotel)
	call netcdf_err(error, 'DEFINING Obs_SNOTEL' )
	error = nf90_put_att(ncid, id_obssnotel, "long_name", "Observed at SNOTEL Station")
	call netcdf_err(error, 'DEFINING Obs_SNOTEL LONG NAME' )
	error = nf90_put_att(ncid, id_obssnotel, "units", "mm")
	call netcdf_err(error, 'DEFINING Obs_SNOTEL UNITS' )
	
	error = nf90_def_var(ncid, 'SNOFCS_atSNOTEL', NF90_DOUBLE, dim_eval, id_forcsnotel)
	call netcdf_err(error, 'DEFINING SNOFCS_atSNOTEL' )
	error = nf90_put_att(ncid, id_forcsnotel, "long_name", "Forecast at SNOTEL Station")
	call netcdf_err(error, 'DEFINING SNOFCS_atSNOTEL LONG NAME' )
	error = nf90_put_att(ncid, id_forcsnotel, "units", "deg")
	call netcdf_err(error, 'DEFINING SNOFCS_atSNOTEL UNITS' )

	error = nf90_def_var(ncid, 'SNOANL_atSNOTEL', NF90_DOUBLE, dim_eval, id_anlsnotel)
	call netcdf_err(error, 'DEFINING SNOANL_atSNOTEL' )
	error = nf90_put_att(ncid, id_anlsnotel, "long_name", "Analysis at SNOTEL Station")
	call netcdf_err(error, 'DEFINING SNOANL_atSNOTEL LONG NAME' )
	error = nf90_put_att(ncid, id_anlsnotel, "units", "mm")
	call netcdf_err(error, 'DEFINING SNOANL_atSNOTEL UNITS' )
	
	error = nf90_def_var(ncid, 'Innov_atSNOTEL', NF90_DOUBLE, dim_eval, id_innovsnotel)
	call netcdf_err(error, 'DEFINING Innov_atSNOTEL' )
	error = nf90_put_att(ncid, id_innovsnotel, "long_name", "Innovation at SNOTEL Station")
	call netcdf_err(error, 'DEFINING Innov_atSNOTEL LONG NAME' )
	error = nf90_put_att(ncid, id_innovsnotel, "units", "mm")
	call netcdf_err(error, 'DEFINING Innov_atSNOTEL UNITS' )

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
	error = nf90_put_var( ncid, id_latsnotel, Lat_SNOTEL)
	call netcdf_err(error, 'WRITING Lat_SNOTEL RECORD' )

	error = nf90_put_var( ncid, id_lonsnotel, Lon_SNOTEL)
	call netcdf_err(error, 'WRITING Lon_SNOTEL RECORD' )

	error = nf90_put_var( ncid, id_obssnotel, Obs_SNOTEL)
	call netcdf_err(error, 'WRITING Obs_SNOTEL RECORD' )

	error = nf90_put_var( ncid, id_forcsnotel, SNOFCS_atSNOTEL)
	call netcdf_err(error, 'WRITING SNOFCS_atSNOTEL RECORD' )

	error = nf90_put_var( ncid, id_anlsnotel, SNOANL_atEvalPts)
	call netcdf_err(error, 'WRITING SNOANL_atEvalPts RECORD' )
	
	error = nf90_put_var( ncid, id_innovsnotel, innov_atEvalPts)
	call netcdf_err(error, 'WRITING innov_atEvalPts RECORD' )

	deallocate(x_data, y_data)
	deallocate(dum2d)

	error = nf90_close(ncid)

 End subroutine Write_DA_Outputs

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
 
