MODULE M_DA_OI

    USE NETCDF
    Use, Intrinsic :: IEEE_ARITHMETIC	
    
    CONTAINS
       
    subroutine Select_Nearest_Observations(RLA, RLO, MYRANK, IDIM, JDIM, &
                                IY, IM, ID, IH, &
                                TILE_NUM,         &
                                LENSFC,	 		&
                                SNOFCS, SWDFCS,  &
                                SNOANL, SWDANL, SNWDEN)
                                
    !----------------------------------------------------------------------
    ! Input: forecast/background states for a single tile by a single MPI process
    ! reads observations: snow depth / SWE and cover
    ! calls observation operator for a given tile 
    !                   
    !----------------------------------------------------------------------
     IMPLICIT NONE
    !
     include 'mpif.h'
    
     INTEGER, intent(in) :: IY, IM, ID, IH, IDIM, JDIM
     INTEGER, intent(In) :: MYRANK, LENSFC
     CHARACTER(LEN=5), INTENT(In) :: TILE_NUM
     INTEGER			 :: IERR
     
     REAL                :: RLA(LENSFC), RLO(LENSFC)
     REAL                :: SNOFCS(LENSFC), SWDFCS(LENSFC), SNWDEN(LENSFC)
     REAL                :: SNOANL(LENSFC), SWDANL(LENSFC)
    
     CHARACTER(len=250) :: snotel_inp_file, dim_name, ghcnd_inp_file, ims_inp_file
     CHARACTER(len=5) :: y_str, m_str, d_Str
     REAL, ALLOCATABLE    :: SWE_SNOTEL(:), SNWD_SNOTEL(:)
     REAL, ALLOCATABLE    :: Lat_SNOTEL(:), Lon_SNOTEL(:), SNOFCS_atSNOTEL(:) !, Ele_SNOTEL(:)
     REAL, ALLOCATABLE    :: SNWD_GHCND(:), SNOFCS_atGHCND(:), SWDFCS_atGHCND(:)
     REAL, ALLOCATABLE	  :: Lat_GHCND(:), Lon_GHCND(:) !, Ele_GHCND(:)
     INTEGER, ALLOCATABLE    :: SNCOV_IMS(:,:), SNOFCS_atIMS(:,:), SWDFCS_atIMS(:,:)
     REAL, ALLOCATABLE	   :: Lat_IMS(:), Lon_IMS(:)	!, Ele(:,:)
     INTEGER :: Num_Snotel, Num_Ghcnd, Num_Ims_Lat, Num_Ims_Lon
     Real	:: max_distance   ! radius_of_influence for observations
    
     !NUM_THREADS = NUM_PARTHDS()  ! ? no include?
    
     CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    
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
     snotel_inp_file = "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/SNOTEL/SNOTEL.SWE.SNWD."// &
                         TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18.nc"
     dim_name = "Site_Id"
     
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
    
     if (myrank==0) then
        PRINT*, "SNOTEL SWE from rank: ", MYRANK
        PRINT*, SWE_SNOTEL	
     endif
     PRINT*,'Finished reading SNOTEL ON RANK: ', MYRANK			
     ! snow depth from GHCND 
     ! GHCND.SNWD.2019100118.nc
    ghcnd_inp_file = "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/GHCND/GHCND.SNWD."// &
                        TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18.nc"
    !dim_name = "Site_Id"
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
    ims_inp_file = "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS/IMS.SNCOV."// &
                        TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18.nc"	
    call Observation_Read_IMS(ims_inp_file,  &
                    Num_Ims_Lat, Num_Ims_Lon,	&
                    SNCOV_IMS,		&
                    Lat_IMS,      &
                    Lon_IMS,		&
                    !Ele,		&
                    MYRANK)
    !	
    PRINT*,'Finished reading IMS ON RANK: ', MYRANK	
    
    ! if (myrank==4) then
    ! 	PRINT*, "IMS SNCOV from rank: ", MYRANK
    ! 	PRINT*, SNCOV_IMS
    ! endif
    ! endif ! if myrnak==0
     
    !  CALL MPI_BCAST(SWE_SNOTEL, Num_Snotel, MPI_REAL, 0, MPI_COMM_WORLD, IERR)
    !  CALL MPI_BCAST(SNWD_SNOTEL, Num_Snotel, MPI_REAL, 0, MPI_COMM_WORLD, IERR)
    !  CALL MPI_BCAST(SNWD_GHCND, Num_Ghcnd, MPI_REAL, 0, MPI_COMM_WORLD, IERR)
    !  CALL MPI_BCAST(SNCOV_IMS, Num_Ims_Lat * Num_Ims_Lon, MPI_INT, 0, MPI_COMM_WORLD, IERR)
     
     CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    ! broadcast observed arrays to the other processes
    
    ! 4.8.20: compute snow density from forecast swe and snwdepth
     SNWDEN = SWDFCS / SNOFCS
    ! Get model states at obs points
     ALLOCATE(SNOFCS_atSNOTEL(Num_Snotel))
     ALLOCATE(SNOFCS_atGHCND(Num_Ghcnd))
     ALLOCATE(SWDFCS_atGHCND(Num_Ghcnd)) 
     ALLOCATE(SNOFCS_atIMS(Num_Ims_Lat * Num_Ims_Lon))
     ALLOCATE(SWDFCS_atIMS(Num_Ims_Lat * Num_Ims_Lon))
     max_distance = 120			!Km 
     call Observation_Operator(RLA, RLO, Lat_SNOTEL, Lon_SNOTEL,   &
                        SNOFCS, SNOFCS_atSNOTEL,            &
                        LENSFC, Num_Snotel, max_distance)
                            !,                 &						intp_mode)
    call Observation_Operator(RLA, RLO, Lat_GHCND, Lon_GHCND,   &
                        SNOFCS, SNOFCS_atGHCND,            &
                        LENSFC, Num_Ghcnd, max_distance)
                            !,                 &						intp_mode) 
    call Observation_Operator(RLA, RLO, Lat_GHCND, Lon_GHCND,   &
                        SWDFCS, SWDFCS_atGHCND,            &
                        LENSFC, Num_Ghcnd, max_distance)
                            !,                 &						intp_mode) 
    call Observation_Operator(RLA, RLO, Lat_IMS, Lon_IMS,    &
                        SNOFCS, SNOFCS_atIMS,            &
                        LENSFC, Num_Ims_Lat * Num_Ims_Lon, max_distance)
                            !,                 &						intp_mode) 
    call Observation_Operator(RLA, RLO, Lat_IMS, Lon_IMS,   &
                        SWDFCS, SWDFCS_atIMS,            &
                        LENSFC, Num_Ims_Lat * Num_Ims_Lon, max_distance)
                            !,                 &						intp_mode) 
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    if (myrank==4) then
        PRINT*, "Background SWE  from rank: ", MYRANK
        PRINT*, SNOFCS	
        PRINT*
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
        PRINT*
        PRINT*, "Background SWE at SNOTEL locations from rank: ", MYRANK
        PRINT*, SNOFCS_atSNOTEL	
    endif
    
    ! 4.13.20 use IEEE_IS_NAN() in assimilation DA to exclude points too far from obs
    
    !  CALL Snow_DA_OI(LUGB,IDIM,JDIM,LENSFC, DELTSFC,  &
    !              IY,IM,ID,IH,FH,IALB,                  &
    !              SNOFCS,  SWDFCS,      &
    !              OROG, OROG_UF,MYRANK)
    
     DEALLOCATE(SWE_SNOTEL, SNWD_SNOTEL, Lat_SNOTEL, Lon_SNOTEL)
     DEALLOCATE(SNOFCS_atSNOTEL, SNWD_GHCND, SNOFCS_atGHCND, SWDFCS_atGHCND) 
     DEALLOCATE(SNCOV_IMS, SNOFCS_atIMS, SWDFCS_atIMS)
     DEALLOCATE(Lat_GHCND, Lon_GHCND, Lat_IMS, Lon_IMS)
    
     PRINT*,'Finished observation operator ON RANK: ', MYRANK
     CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    
     STOP !RETURN
    
     END subroutine Snow_Analysis
     
     ! Get model states at obs points
     ! Warning: This assumes all distance coordinates are valid; 
     ! do quality control of coordinates beforehand
     SUBROUTINE Observation_Operator(RLA, RLO, Lat_Obs, Lon_Obs,   &
                            SNWD_back,  				&
                            SNWD_atObs, 		&	
                            LENSFC, num_Obs, max_distance) !,      &
                            !intp_mode) 
                            !  SWE_atObs  SWE_back
    
        !USE intrinsic::ieee_arithmetic
        INTEGER :: LENSFC, num_Obs
        INTEGER :: indx, min_indx
        Real    :: distArr(LENSFC), haversinArr(LENSFC)
        Real 	:: d_latArr(LENSFC), d_lonArr(LENSFC)
        real	:: max_distance   ! radius_of_influence
        !CHARACTER(len=*), intent(In):: intp_mode
        Real, Intent(In) 	:: RLA(LENSFC), RLO(LENSFC), Lon_Obs(num_Obs)  ! don't want to alter these
        Real 	:: RLO_2(LENSFC), Lat_Obs(num_Obs), Lon_Obs_2(num_Obs)
        Real    :: SNWD_back(LENSFC), SNWD_atObs(num_Obs)		
        Real 	:: RLA_rad(LENSFC), RLO_rad(LENSFC)
        Real 	:: Lat_Obs_rad(num_Obs), Lon_Obs_rad(num_Obs)	
        Real(16), Parameter :: PI_16 = 4 * atan (1.0_16)	
        Real(16), Parameter :: pi_div_180 = PI_16/180.0
        Real, Parameter		:: earth_rad = 6371.
        ! PRINT*, "PI: ", PI_16
        ! PRINT*, "PI / 180: ", pi_div_180
    
        !Fill background values to nan (to differentiate those htat don't have value)
        SNWD_atObs = IEEE_VALUE(SNWD_atObs, IEEE_QUIET_NAN)	
        
        !if intp_mode == 'near'		! [bilinear, customInterpol])
    
        ! RLO from 0 to 360 (no -ve lon)
        ! Do indx = 1, num_Obs 
            ! Lon_Obs[Lon_Obs<0]= 360.0 + Lon_Obs[Lon_Obs<0]
        ! 	if (Lon_Obs(indx) < 0) Lon_Obs_2(indx) = 360. + Lon_Obs(indx)
        ! end do
        Do jndx = 1, LENSFC 
            if (RLO(jndx) > 180) RLO_2(jndx) = RLO(jndx) - 360.0 
        end do
        Do zndx = 1, num_Obs 
            if (Lon_Obs(zndx) > 180) Lon_Obs_2(zndx) = Lon_Obs(zndx) - 360.0
        end do
        ! at each obs point compute its distance from RLA/RLO pairs 
        ! then find the position of the minimum
    
        ! shortest distance over sphere using great circle distance 	
        RLA_rad =  pi_div_180 * RLA
        RLO_rad =  pi_div_180 * RLO_2
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
            Do jndx = 1, LENSFC 
                if (haversinArr(jndx) > 1) haversinArr(jndx) = 1 ! ensure numerical errors don't make h>1
            end do
            distArr = 2 * earth_rad * asin(sqrt(haversinArr))		
            !distArr = (Lat_Obs(indx) - RLA)**2 + (Lon_Obs_2(indx) - RLO)**2 
            min_indx = MINLOC(distArr, dim = 1)  !, MASK=ieee_is_nan(distArr))
            if(distArr(min_indx) < max_distance) then
                SNWD_atObs(indx) = SNWD_back(min_indx) 
                !SWE_atObs[indx] = SWE_back[min_indx]
            endif
        end do
        
        
        RETURN
        
     END SUBROUTINE Observation_Operator
     !
     SUBROUTINE Observation_Read_IMS(inp_file,  &
                    DIM_LEN_lat, DIM_LEN_lon,  &
                    SNCOV_IMS,		&
                    Lat,      &
                    Lon,		&
                    ! Ele		&
                    MYRANK)
        IMPLICIT NONE
    
        include 'mpif.h'		  
    
        !Open netCDF for a snotel and read the SWE, SnowDepth,..., Lat, Lon, at a given datetime
        !ToDO: Can you use variable length char array ?
        CHARACTER(LEN=*), Intent(In)      :: inp_file !, dim_name
    
        INTEGER                :: ERROR, NCID
        INTEGER                :: MYRANK
        INTEGER                :: ID_DIM, ID_VAR
        INTEGER, Intent(Out)   :: DIM_LEN_lat, DIM_LEN_lon
        
        ! ToDO: ims snow cover is of 'byte' type (Chcek the right one)
        INTEGER, ALLOCATABLE, Intent(Out)    :: SNCOV_IMS(:,:) 
        REAL, ALLOCATABLE, Intent(Out)	   :: Lat(:), Lon(:)	!, Ele(:,:)
    
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
    
    
        ALLOCATE(SNCOV_IMS(DIM_LEN_lon, DIM_LEN_lat))
        ALLOCATE(Lat(DIM_LEN_lat))
        ALLOCATE(Lon(DIM_LEN_lon))
        !ALLOCATE(Ele(DIM_LEN_lat, DIM_LEN_lon))
    
        ERROR=NF90_INQ_VARID(NCID, 'Band1', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, SNCOV_IMS, start = (/ 1, 1 /), &
                                count = (/ DIM_LEN_lon, DIM_LEN_lat/))
        CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV RECORD' )
    
        ERROR=NF90_INQ_VARID(NCID, 'lat', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lat var ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, Lat)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lat RECORD' )
    
        ERROR=NF90_INQ_VARID(NCID, 'lon', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lon var ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, Lon)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lon RECORD' )
    
        ! need to read corresponding elevation values 
        ! ERROR=NF90_INQ_VARID(NCID, 'Elevation', ID_VAR)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation ID' )
        ! ERROR=NF90_GET_VAR(NCID, ID_VAR, Ele)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation RECORD' )
    
        ERROR = NF90_CLOSE(NCID)
                  
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
    
        ! need to read corresponding elevation values 
        ! ERROR=NF90_INQ_VARID(NCID, 'Elevation', ID_VAR)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation ID' )
        ! ERROR=NF90_GET_VAR(NCID, ID_VAR, Ele_SNOTEL)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation RECORD' )
    
        ERROR = NF90_CLOSE(NCID)
                  
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
     
     SUBROUTINE Snow_DA_OI(LUGB,IDIM,JDIM,LENSFC, DELTSFC,  &
                            IY,IM,ID,IH,FH,                  &
                            MYRANK,				&
                            PRNT, 				&
                            TILE_NUM,                 &
                            SNOFCS,  SWDFCS,      &
                            OROG, RLA, RLO)
                            
    !--------------------------------------------------------------------------------
    ! 
    !  Snow vars:
    ! 
    !   RLA, RLO
    !	SNOFCS(LENSFC), SWDFCS(LENSFC), ALBFCS(LENSFC,4), ABSFCS(LENSFC) 
    !	OROG(LENSFC), OROG_UF(LENSFC), ZORFCS(LENSFC)
    !
    !     SNO  .. Liquid-equivalent snow depth  => 
    !     SWD  .. Actual snow depth
    !	  ABS  .. Maximum snow albedo
    !     ALB  .. Snow-free albedo
    !
    !     OROG .. Orography
    !     OROG_UF .. Orography unfiltered
    !     RLA,
    !	  RLO
    !     ZOR  .. Surface roughness length!
    !     VET  .. Vegetation type
    !     VEG  .. Vegetation cover
    !
    !     SLI  .. LAND/SEA/SEA-ICE mask. (1/0/2 respectively)
    !     VMN  .. Vegetation cover minimum
    !     VMX  .. Vegetation cover maximum
    !     SLP  .. Slope type
    ! ------------------------------------------------------------------------------------
     !USE READ_WRITE_DATA
    
     IMPLICIT NONE
    
     include 'mpif.h'
    
     INTEGER, INTENT(IN) :: IDIM, JDIM, LENSFC
     INTEGER, INTENT(IN) :: LUGB, IY, IM, ID, IH
     INTEGER, INTENT(IN) ::  MYRANK
     LOGICAL, INTENT(IN) :: PRNT 
     REAL, INTENT(IN)    :: FH, DELTSFC
    
     CHARACTER(LEN=5)    :: TILE_NUM
    
     INTEGER             :: I, IERR
    
     REAL                :: OROG(LENSFC)
     REAL                :: SNOFCS(LENSFC), SWDFCS(LENSFC) 
     REAL                :: RLA(LENSFC), RLO(LENSFC)
    
     
    
     !TYPE(Observation_Str)     :: Obs_Str
    
    
     RETURN
    
     END SUBROUTINE Snow_DA_OI
    
     ! copied from read_write_data.f90---may be better to 'include'
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
    
     END MODULE M_DA_OI
     
    