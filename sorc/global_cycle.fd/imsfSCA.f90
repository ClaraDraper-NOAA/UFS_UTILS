Program imsfSCA

    IMPLICIT NONE

    USE NETCDF

    Use, Intrinsic :: IEEE_ARITHMETIC	

    include 'mpif.h'
	
	integer, parameter :: dp = kind(1.d0)

	INTEGER   :: IY, IM, ID, IH, IDIM, JDIM, MYRANK, LENSFC

	INTEGER			 :: IERR
	
	CHARACTER(len=250)   :: ims_inp_file, ims_inp_file_indices
	CHARACTER(len=250) 	 :: da_out_file
	CHARACTER(len=5)     :: y_str, m_str, d_Str, fvs_tile
    Character(LEN=3)       :: rank_str

	!--ims resampled to each grid REAL            	 :: lat_min, lat_max, lon_min, lon_max 
    Real, Allocatable		:: SNCOV_IMS(:)  !, SNOFCS_atIMS(:), SWDFCS_atIMS(:)
    Integer :: num_sub

    PRINT*,"STARTING Snow analysis on RANK "
    IY = 2019
    IM = 11
    ID = 15
    IH = 18
    IDIM = 768
    JDIM = 768 
    MYRANK = 4
    LENSFC = IDIM * JDIM
    Allocate(SNCOV_IMS(LENSFC))

    write(y_str, "(I4)") IY
	write(m_str, "(I0.2)") IM
	write(d_str, "(I0.2)") ID
	write(fvs_tile, "(I3)") IDIM


	num_sub = 30   !627  ! (max) number of ims subcells within a tile grid cell
	ims_inp_file = "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS/IMS.SNCOV."// &
						TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18.nc"
	ims_inp_file_indices = "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS/C"// &
						TRIM(ADJUSTL(fvs_tile))//".IMS.Indices."//TRIM(TILE_NUM)//".nc"			
	Call Observation_Read_IMS_Full(ims_inp_file, ims_inp_file_indices, &
				JDIM, IDIM, num_sub, &
				SNCOV_IMS)

	PRINT*,'Finished reading IMS ON RANK: ', MYRANK
	
	
	Write(rank_str, '(I3.3)') (MYRANK+1)
	da_out_file = "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/Analysis/SNOANL.SWE.SNWD."// &
						TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18_tile"//rank_str//".nc"
	
	call Write_DA_Outputs(da_out_file, IDIM, JDIM, LENSFC, SNCOV_IMS) !updated snocov

    
    Deallocate(SNCOV_IMS)

	PRINT*,'Finished'

	STOP

END Program imsfSCA


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
				n_lat, n_lon, num_sub, &
				SNCOV_IMS)
				! Ele		&
				
	IMPLICIT NONE

	include 'mpif.h'		  

	!ToDO: Can you use variable length char array ?
	CHARACTER(LEN=*), Intent(In)   :: inp_file, inp_file_indices !, dim_name
	INTEGER, Intent(In)            :: n_lat, n_lon, num_sub
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

 

 ! the following code based on write_data() in read_write_data.f90
 Subroutine Write_DA_Outputs(output_file, idim, jdim, lensfc, SNCOV_IMS)

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

	Real, intent(in)            ::  SNCOV_IMS(lensfc)

	integer                     :: fsize=65536, inital=0
	integer                     :: header_buffer_val = 16384
	integer                     :: dims_3d(3), dims_strt(3), dims_end(3)
	integer                     :: error, i, ncid
	integer                     :: dim_x, dim_y, dim_time
	integer                     :: id_x, id_y, id_time
	integer                     ::  id_imscov


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

	error = nf90_def_var(ncid, 'imsfSCA', NF90_DOUBLE, dims_3d, id_imscov)
	call netcdf_err(error, 'DEFINING imsfSCA' )
	error = nf90_put_att(ncid, id_imscov, "long_name", "IMS fractional Snow Covered Area")
	call netcdf_err(error, 'DEFINING imsfSCA LONG NAME' )
	error = nf90_put_att(ncid, id_imscov, "units", "-")
	call netcdf_err(error, 'DEFINING imsfSCA UNITS' )

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
	
	dum2d = reshape(SNCOV_IMS, (/idim,jdim/))
	error = nf90_put_var( ncid, id_imscov, dum2d, dims_strt, dims_end)
	call netcdf_err(error, 'WRITING imsfSCA RECORD' )

	deallocate(x_data, y_data)
	deallocate(dum2d)

	error = nf90_close(ncid)

 End subroutine Write_DA_Outputs

 