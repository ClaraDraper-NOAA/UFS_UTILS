! Program ImsGridMap

!     Implicit None
!     include 'mpif.h'
!     use MPI 
!     Use NETCDF

!     Integer :: num_procs, ierr, myrank

!     Integer, Intent(In)     :: ntot_ims, n_lon, n_lat, ntot, myrank 
!     Real, Intent(In)        :: lat_grid_ims(ntot_ims), lon_grid_ims(ntot_ims), lat_grid(ntot), lon_grid(ntot)
!     Integer, Intent(Out)    :: grid_dat(ntot_ims, 10), data_grid_ims(ntot_ims)
    
!     Integer   :: i, j, indx, num_loc_counter, data_grid_ims(ntot_ims)
!     Real      :: lat_grid_1(ntot), lon_grid_1(ntot), maxdist(ntot), dist(ntot_ims), outside_dist(ntot_ims)

!     CALL MPI_INIT(ierr)
!     CALL MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
!     CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)


!     file_in = 'snow.20200101.180000.sfcanl_data.tile' + str(til + 1) + '.nc'
!     file_out = "IMS.Indices.tile" + str(til + 1) + '.nc'
!     cmdString = "copy " + file_in + " " + file_out
!     callSubprocess_ps(cmdString)

!     lat_min = np.amin(lat_grid[til]); lat_max = np.amax(lat_grid[til])
!     lon_min = np.amin(lon_grid[til]); lon_max = np.amax(lon_grid[til])

!     # find index  with min abs value difference
!     Lat_minmax_Diff = np.abs(lat_grid_ims_1 - lat_min)
!     minlat_indx = np.nanargmin(Lat_minmax_Diff)
!     Lat_minmax_Diff = np.abs(lat_grid_ims_1 - lat_max)
!     maxlat_indx = np.nanargmin(Lat_minmax_Diff)
!     Lon_minmax_Diff = np.abs(lon_grid_ims_1 - lon_min)
!     minlon_indx = np.nanargmin(Lon_minmax_Diff)
!     Lon_minmax_Diff = np.abs(lon_grid_ims_1 - lon_max)
!     maxlon_indx =  np.nanargmin(Lon_minmax_Diff)
!     print("lat index ", minlat_indx," ", maxlat_indx)
!     print("lon index ", minlon_indx, " ", maxlon_indx)
!     DIM_LEN_lat = 1 + abs(maxlat_indx - minlat_indx)
!     DIM_LEN_lon = 1 + abs(maxlon_indx - minlon_indx)

!     lat_grid_ims = np.full((DIM_LEN_lat, DIM_LEN_lon), np.nan)
!     lon_grid_ims = np.full((DIM_LEN_lat, DIM_LEN_lon), np.nan)
!     #data_grid_ims_til = np.full((DIM_LEN_lat, DIM_LEN_lon), np.nan)
!     if (minlat_indx > maxlat_indx or minlon_indx > maxlon_indx):
!         print("order of indices wrong")
!         exit(1)
!     for ilon in range(DIM_LEN_lon):
!         lat_grid_ims[:, ilon] = lat_grid_ims_1[minlat_indx:maxlat_indx+1]
!     for ilat in range(DIM_LEN_lat):
!         lon_grid_ims[ilat, :] = lon_grid_ims_1[minlon_indx:maxlon_indx+1]

!     call resample_to_model_tiles(ntot_ims, lat_grid_ims, lon_grid_ims, &
!                                  lat_grid, lon_grid, n_lon, n_lat, ntot, myrank, &
!                                  grid_dat)


! End program ImsGridMap


subroutine select_imscells_within_tile_grid(lat_grid_ims, lon_grid_ims, &
                                   lat_grid, lon_grid, n_lat_out, n_lon_out, num_sub, &
                                   n_lat, n_lon, ntot_ims, &  !myrank, &
                                   grid_dat)
    
    Use, Intrinsic :: IEEE_ARITHMETIC
    
    Implicit None
    
    Integer, Intent(In)     :: ntot_ims, n_lat, n_lon, num_sub, n_lat_out, n_lon_out    !n_lon, n_lat,     !, myrank 
    Real, Intent(In)        :: lat_grid_ims(ntot_ims), lon_grid_ims(ntot_ims) 
    Real, Intent(In)        :: lat_grid(n_lat, n_lon), lon_grid(n_lat, n_lon)
    Integer, Intent(Out)    :: grid_dat(n_lat_out, n_lon_out, num_sub)
    
    Integer   :: j, iy, jx, num_loc_counter
    Real      :: min_y, max_y, min_x, max_x

    ! #lon_grid[lon_grid > 180.] = lon_grid[lon_grid > 180.] - 360.
    print*, "total number ims grid cells = ", ntot_ims

    Do iy=2, n_lat-1, 2
        
        !print*, "process: ", myrank, "loop ", indx

        Do jx=2, n_lon-1, 2
            print*, "iy ", iy/2, " jx ", jx/2

            min_y = MIN(lat_grid(iy-1, jx-1), lat_grid(iy-1, jx), lat_grid(iy-1, jx+1), &
                        lat_grid(iy, jx-1), lat_grid(iy, jx), lat_grid(iy, jx+1), &
                        lat_grid(iy+1, jx-1), lat_grid(iy+1, jx), lat_grid(iy+1, jx+1))
            max_y = MAX(lat_grid(iy-1, jx-1), lat_grid(iy-1, jx), lat_grid(iy-1, jx+1), &
                        lat_grid(iy, jx-1), lat_grid(iy, jx), lat_grid(iy, jx+1), &
                        lat_grid(iy+1, jx-1), lat_grid(iy+1, jx), lat_grid(iy+1, jx+1))
            min_x = MIN(lon_grid(iy-1, jx-1), lon_grid(iy-1, jx), lon_grid(iy-1, jx+1), &
                        lon_grid(iy, jx-1), lon_grid(iy, jx), lon_grid(iy, jx+1), &
                        lon_grid(iy+1, jx-1), lon_grid(iy+1, jx), lon_grid(iy+1, jx+1))
            max_x = MAX(lon_grid(iy-1, jx-1), lon_grid(iy-1, jx), lon_grid(iy-1, jx+1), &
                        lon_grid(iy, jx-1), lon_grid(iy, jx), lon_grid(iy, jx+1), &
                        lon_grid(iy+1, jx-1), lon_grid(iy+1, jx), lon_grid(iy+1, jx+1))

            if(min_y > max_y) then
                print*, "Error! min y value > max y" 
                stop
            endif
            if(min_x > max_x) then 
                print*, "Error! min x value > max x"
                stop
            endif
            num_loc_counter = 0   
            Do j = 1, ntot_ims
                if((lat_grid_ims(j) .ge. min_y) .AND. (lat_grid_ims(j) .le. max_y) .AND. &
                    (lon_grid_ims(j) .ge. min_x) .AND. (lon_grid_ims(j) .le. max_x)) then  
                    num_loc_counter = num_loc_counter + 1
                    grid_dat(iy/2, jx/2, num_loc_counter+1) =  j             
                endif
                if (num_loc_counter > num_sub-2) then 
                    print*, "warning more than ",num_sub-2," grids found!"
                    exit
                end if
            End do
            
            grid_dat(iy/2, jx/2, 1) =  num_loc_counter ! first location, num obs
        End do

    End do

    return !grid_dat

End subroutine select_imscells_within_tile_grid

subroutine select_imscells_within_tile_grid_dist(lat_grid_ims, lon_grid_ims, &
                                    lat_grid, lon_grid, n_lat_out, n_lon_out, num_sub, &
                                    n_lat, n_lon, ntot_ims, &  !myrank, &
                                    grid_dat)       


    Use, Intrinsic :: IEEE_ARITHMETIC

    Implicit None

    Integer, Intent(In)     :: ntot_ims, n_lat, n_lon, num_sub, n_lat_out, n_lon_out    !n_lon, n_lat,     !, myrank 
    Real, Intent(In)        :: lat_grid_ims(ntot_ims), lon_grid_ims(ntot_ims) 
    Real, Intent(In)        :: lat_grid(n_lat, n_lon), lon_grid(n_lat, n_lon)
    Integer, Intent(Out)    :: grid_dat(n_lat_out, n_lon_out, num_sub)

    Integer   :: i, j, iy, jx, num_loc_counter, data_grid_ims(ntot_ims)
    Real      :: diag_1, diag_2, diag_g, dist_g(ntot_ims)  !min_y, max_y, min_x, max_x, temp_y, temp_x, 
    !Real      :: lat_grid_1(ntot), lon_grid_1(ntot), maxdist(ntot), dist(ntot_ims) !, outside_dist(ntot_ims)

    ! #lon_grid[lon_grid > 180.] = lon_grid[lon_grid > 180.] - 360.

    ! #grid_dat = np.full([n_lon * n_lat], np.nan)
    ! # lat_grid_ims[np.isnan(data_grid_ims)] = np.nan
    ! # lon_grid_ims[np.isnan(data_grid_ims)] = np.nan

    ! lat_grid_1 = IEEE_VALUE(lat_grid_1, IEEE_QUIET_NAN)
    ! lon_grid_1 = IEEE_VALUE(lon_grid_1, IEEE_QUIET_NAN)
    !grid_dat = IEEE_VALUE(grid_dat, IEEE_QUIET_NAN)

    ! lat_grid_1(1:ntot-1) = lat_grid(2:ntot)
    ! lat_grid_1(ntot) = lat_grid(ntot-1)
    ! lon_grid_1(1:ntot-1) = lon_grid(2:ntot)
    ! lon_grid_1(ntot) = lon_grid(ntot-1)

    Do i=1, ntot_ims
        data_grid_ims(i) = i   !indices of array 
    End do

    !maxdist = (lat_grid_1 - lat_grid) ** 2 + (lon_grid_1 - lon_grid) ** 2

    Do iy=2, n_lat-1, 2

        !print*, "process: ", myrank, "loop ", indx

        Do jx=2, n_lon-1, 2
            print*, "iy ", iy/2, " jx ", jx/2
            !dist = (lat_grid(indx) - lat_grid_ims) ** 2 + (lon_grid(indx) - lon_grid_ims) ** 2
            ! min_y = MIN(lat_grid(iy-1, jx-1), lat_grid(iy-1, jx), lat_grid(iy-1, jx+1))
            ! max_y = MAX(lat_grid(iy+1, jx-1), lat_grid(iy+1, jx), lat_grid(iy+1, jx+1))
            ! min_x = MIN(lon_grid(iy-1, jx-1), lon_grid(iy, jx-1), lon_grid(iy+1, jx-1))
            ! max_x = MAX(lon_grid(iy-1, jx+1), lon_grid(iy, jx+1), lon_grid(iy+1, jx+1))

            ! if(min_y > max_y) then
            !     temp_y = min_y
            !     min_y = max_y
            !     max_y = temp_y
            !     !stop
            ! end if
            ! if(min_x > max_x) then
            !     temp_x = min_x
            !     min_x = max_x
            !     max_x = temp_x
            !     !stop
            ! end if
            ! if(min_y > max_y) print*, "warning! min y value > max y" 
            ! if(min_x > max_x) print*, "Warning! min x value > max x"
            diag_1 = (lat_grid(iy-1, jx-1) - lat_grid(iy+1, jx+1)) ** 2 + &
                     (lon_grid(iy-1, jx-1) - lon_grid(iy+1, jx+1)) ** 2
            diag_2 = (lat_grid(iy+1, jx-1) - lat_grid(iy-1, jx+1)) ** 2 + &
                     (lon_grid(iy+1, jx-1) - lon_grid(iy-1, jx+1)) ** 2
            diag_g = MAX(diag_1, diag_2) / 4.

            dist_g = (lat_grid(iy, jx) - lat_grid_ims) ** 2 + (lon_grid(iy, jx) - lon_grid_ims) ** 2
            !dist_g =sqrt(dist_g)
            num_loc_counter = 0
            Do j = 1, ntot_ims
                if(dist_g(j) <= diag_g) then  
                    num_loc_counter = num_loc_counter + 1
                    grid_dat(iy/2, jx/2, num_loc_counter+1) =  data_grid_ims(j)             
                endif
                if (num_loc_counter > num_sub-2) then 
                    print*, "warning more than ",num_sub-2," grids found!"
                    exit
                end if
            End do
            ! Do j = 1, ntot_ims
            !     if((lat_grid_ims(j) > min_y) .AND. (lat_grid_ims(j) < max_y) .AND. &
            !         (lon_grid_ims(j) > min_x) .AND. (lon_grid_ims(j) < max_x)) then  
            !         num_loc_counter = num_loc_counter + 1
            !         grid_dat(iy/2, jx/2, num_loc_counter+1) =  data_grid_ims(j)             
            !     endif
            !     if (num_loc_counter > 18) then 
            !         print*, "warning more than 18 grids found!"
            !         exit
            !     end if
            ! End do
            ! Do j = 1, ntot_ims
            !     if((dist(j) < maxdist(indx)) .AND. (.NOT. IEEE_IS_NAN(dist(j)))) then  
            !         num_loc_counter = num_loc_counter + 1
            !         grid_dat(indx, num_loc_counter+1) =  data_grid_ims(j)             
            !     endif
            !     if (num_loc_counter > 8) then 
            !         print*, "warning more than 8 grids found!"
            !         exit
            !     end if
            ! End do
            grid_dat(iy/2, jx/2, 1) =  num_loc_counter ! first location, num obs
        End do

    End do

    return !grid_dat

End subroutine select_imscells_within_tile_grid_dist

subroutine resample_to_model_tiles_intrp(data_grid_ims, data_grid_ims_ind, &
                                            ntot_ims, n_lat, n_lon, num_sub, &  !myrank, &
                                            grid_dat)


    Use, Intrinsic :: IEEE_ARITHMETIC

    Implicit None

    Integer, Intent(In)     :: ntot_ims, n_lat, n_lon, num_sub
    Integer, Intent(In)     :: data_grid_ims(ntot_ims), data_grid_ims_ind(n_lat, n_lon, num_sub) 
    Real, Intent(Out)       :: grid_dat(n_lat, n_lon)

    Integer   :: jc, iy, jx, num_loc_counter
    
    grid_dat = IEEE_VALUE(grid_dat, IEEE_QUIET_NAN)

    Do iy=1, n_lat
    !print*, "process: ", myrank, "loop ", indx

        Do jx=1, n_lon
            
            num_loc_counter = data_grid_ims_ind(iy, jx, 1)
            if (num_loc_counter < 1) then 
                print*, "no matching values!"
                cycle
            end if

            print*, "iy ", iy, " jx ", jx
            grid_dat(iy, jx) = 0.
            Do jc = 2, num_loc_counter+1
                grid_dat(iy, jx) =  grid_dat(iy, jx) + data_grid_ims(data_grid_ims_ind(iy, jx, jc))              
            End do

            grid_dat(iy, jx) =  grid_dat(iy, jx) / num_loc_counter ! first location, num obs

        End do

    End do

    return !grid_dat

End subroutine resample_to_model_tiles_intrp

        
