
integer(2) function compar(a1, a2) 

    IMPLICIT NONE
    Real :: a1, a2

    if (a1 .LT. a2) then 
        compar = -1
    else if (a1 .GT. a2) then
        compar = 1
    else
        compar = 0
    endif

end function compar

MODULE M_DA_OI

    USE NETCDF
    Use, Intrinsic :: IEEE_ARITHMETIC
    use, Intrinsic :: iso_c_binding
    use IFPORT
    
    Logical, parameter :: print_debug = .False.
    
    Integer(2), External :: compar
    
    CONTAINS
       
    subroutine nearest_Observations_Locations(RLA_jndx, RLO_jndx,    &
        Lat_Obs, Lon_Obs, num_Obs, max_distance, max_num_loc,  				&
        Stdev_back, Stdev_Obs_depth, obs_tolerance,                 &
        SNOFCS_atObs, OBS_atOBs,                                 &
        Loc_backSt_atObs,  num_loc) !,      &LENSFC,
        !intp_mode) 
        !  SWE_atObs  SWE_back
        IMPLICIT NONE

        Real, Intent(In) 	:: RLA_jndx, RLO_jndx  ! don't want to alter these
        Real, Intent(In) 	:: Lat_Obs(num_Obs), Lon_Obs(num_Obs)
        Integer, Intent(In) :: num_Obs, max_num_loc
        Real, Intent(In)	:: max_distance   ! radius_of_influence
        Real, Intent(In) 	:: Stdev_back, Stdev_Obs_depth, obs_tolerance
        Real, Intent(In)    :: SNOFCS_atObs(num_Obs)
        Real, Intent(In)    :: OBS_atOBs(num_Obs)
        Integer, Allocatable, Intent(Out)    :: Loc_backSt_atObs(:)
        !Real, Allocatable, Intent(Out)    :: dist_atObs(:)
        Integer, Intent(Out) :: num_loc
        
        Integer :: indx, jndx, zndx, num_loc_counter        
        Real    :: distArr(num_Obs), haversinArr(num_Obs)
        Real 	:: d_latArr(num_Obs), d_lonArr(num_Obs)
        Real    :: Lon_Obs_2(num_Obs)   !RLO_2_jndx, 
        Real 	:: RLA_rad_jndx, RLO_rad_jndx
        Real 	:: Lat_Obs_rad(num_Obs), Lon_Obs_rad(num_Obs)	
        Real(16), Parameter :: PI_16 = 4 * atan (1.0_16)	
        Real(16), Parameter :: pi_div_180 = PI_16/180.0
        Real, Parameter		:: earth_rad = 6371.
        Real                :: innov_criteria
        Real, Allocatable   :: dist_atObs(:)
        Real, Allocatable   :: dist_atObs_dummy(:)
        Integer, Allocatable   :: Loc_backSt_atObs_dummy(:)
        Real                 :: max_value

        INTEGER(SIZEOF_SIZE_T)   :: i_len, i_isize  !C_SIZE_T
        i_isize = SIZEOF(max_value)
        !print*, "size of real: ", i_isize

        !Integer              :: loc_indices(max_num_loc)
        
        ! PRINT*, "PI: ", PI_16
        ! PRINT*, "PI / 180: ", pi_div_180	
    
        !Do jndx = 1, LENSFC 
        !if (RLO_jndx > 180) RLO_2_jndx = RLO_jndx - 360.0 
        !end do
        !print*, " here 1"
        Lon_Obs_2 = Lon_Obs
        Where(Lon_Obs_2 < 0) Lon_Obs_2 = 360. + Lon_Obs_2
        ! Do zndx = 1, num_Obs 
        !     if (Lon_Obs(zndx) < 0) Lon_Obs_2(zndx) = 360. + Lon_Obs(zndx)
        ! end do

        ! at each obs point compute its distance from RLA/RLO pairs 
        ! then find the position of the minimum
        !print*, " here 2"
        ! shortest distance over sphere using great circle distance 	
        RLA_rad_jndx =  pi_div_180 * RLA_jndx
        RLO_rad_jndx =  pi_div_180 * RLO_jndx
        Lat_Obs_rad =  pi_div_180 * Lat_Obs
        Lon_Obs_rad =  pi_div_180 * Lon_Obs_2   
        
        ! https://en.wikipedia.org/wiki/Haversine_formula
        ! https://www.geeksforgeeks.org/program-distance-two-points-earth/
        ! Distance, d = R * arccos[(sin(lat1) * sin(lat2)) + cos(lat1) * cos(lat2) * cos(long2 – long1)]
        ! dist = 2 * R * asin { sqrt [sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2]}
        !Do jndx = 1, LENSFC 
        d_latArr = (Lat_Obs_rad - RLA_rad_jndx) / 2.
        d_lonArr = (Lon_Obs_rad - RLO_rad_jndx) / 2.
        haversinArr = sin(d_latArr)**2 + cos(Lat_Obs_rad) * cos(RLA_rad_jndx) * sin(d_lonArr)**2
        Where (haversinArr > 1) haversinArr = 1.
        ! Do indx = 1, num_Obs 
        !     if (haversinArr(indx) > 1) haversinArr(indx) = 1 ! ensure numerical errors don't make h>1
        ! end do
        Where (haversinArr < 0) haversinArr = 0.
        ! Do indx = 1, num_Obs 
		! 	if (haversinArr(indx) < 0) haversinArr(indx) = 0 ! ensure <0
		! end do
        distArr = 2 * earth_rad * asin(sqrt(haversinArr))
        
        innov_criteria =  obs_tolerance * sqrt(Stdev_back**2 + Stdev_Obs_depth**2)
        ! 4.15.20: can you do the following without loop?       
        num_loc_counter = 0
        Do indx = 1, num_Obs	
            if((distArr(indx) < max_distance) .AND. &
               (.NOT. IEEE_IS_NAN(SNOFCS_atObs(indx))) .AND. &
               (.NOT. IEEE_IS_NAN(OBS_atOBs(indx)/OBS_atOBs(indx))) .AND. &    ! this excludes 0's ?
               (abs(OBS_atOBs(indx) - SNOFCS_atObs(indx)) < innov_criteria)) then                
                num_loc_counter = num_loc_counter + 1
            endif
        End do
        num_loc = num_loc_counter
        Allocate(Loc_backSt_atObs_dummy(num_loc))
        Allocate(dist_atObs_dummy(num_loc))
        jndx = 1
        Do indx = 1, num_Obs	
            if((distArr(indx) < max_distance) .AND. &
               (.NOT. IEEE_IS_NAN(SNOFCS_atObs(indx))) .AND. &
               (.NOT. IEEE_IS_NAN(OBS_atOBs(indx)/OBS_atOBs(indx))) .AND. &
               (abs(OBS_atOBs(indx) - SNOFCS_atObs(indx)) < innov_criteria)) then 
                Loc_backSt_atObs_dummy(jndx) = indx
                dist_atObs_dummy(jndx) = distArr(indx)
                jndx = jndx  + 1
            endif
        End do
        
        ! if num of obs > 50, choose the 50 closest obs
        if (num_loc > max_num_loc) then 
            Allocate(dist_atObs(num_loc))
            dist_atObs = dist_atObs_dummy
            i_len = num_loc
            Call QSORT (dist_atObs_dummy, i_len, i_isize, compar)
            ! print*, "unsorted dist array ",  dist_atObs
            ! print*
            ! print*, "Sorted dist array ",  dist_atObs_dummy
            max_value = dist_atObs_dummy(max_num_loc)
            !loc_indices = findloc(dist_atObs <= max_value) 
            Allocate(Loc_backSt_atObs(max_num_loc)) 
            indx = 1     
            Do jndx = 1, num_loc
                if(dist_atObs(jndx) <= max_value) then
                    Loc_backSt_atObs(indx) = Loc_backSt_atObs_dummy(jndx)
                    indx = indx + 1
                endif
            End do
            ! print*, "unsorted locations ",  Loc_backSt_atObs_dummy
            ! print*
            ! print*, "Sorted locations ",  Loc_backSt_atObs
            num_loc = max_num_loc			
            Deallocate(dist_atObs)
        else
            Allocate(Loc_backSt_atObs(num_loc)) 
            Loc_backSt_atObs = Loc_backSt_atObs_dummy
        endif

        Deallocate(dist_atObs_dummy, Loc_backSt_atObs_dummy)
        
        RETURN

    End subroutine nearest_Observations_Locations

    subroutine select_max_num_obs(dist_atObs, SNOFCS_atObs, OBS_atOBs, &
        num_loc, max_num_loc,  				&
        out_SNOFCS_atObs, out_OBS_atOBs) 

        IMPLICIT NONE
    
        !Integer(2), External :: compar
        !USE intrinsic::ieee_arithmetic
        
        Real, Intent(In)    :: SNOFCS_atObs(num_loc), OBS_atOBs(num_loc), dist_atObs(num_loc)
        Integer, Intent(In) :: num_loc, max_num_loc
        Real, Intent(Out)    :: out_SNOFCS_atObs(max_num_loc), out_OBS_atOBs(max_num_loc)
        Real                 :: dist_atObs_dummy(num_loc)
        Real                 :: max_value
        !Integer              :: loc_indices(max_num_loc)
        Integer :: indx, jndx

        INTEGER(SIZEOF_SIZE_T)   :: i_len, i_isize    !C_SIZE_T
        i_isize = SIZEOF(max_value)
        i_len = num_loc
        
        dist_atObs_dummy = dist_atObs
        Call QSORT (dist_atObs_dummy, i_len, i_isize, compar)
        max_value = dist_atObs_dummy(max_num_loc)
        !loc_indices = findloc(dist_atObs <= max_value) 
        indx = 1       
        Do jndx = 1, max_num_loc
            if(dist_atObs(jndx) <= max_value) then
                out_SNOFCS_atObs(indx) = SNOFCS_atObs(jndx)
                out_OBS_atOBs(indx) = OBS_atOBs(jndx)
                indx = indx + 1
            endif
        End do
        
    End subroutine select_max_num_obs

    subroutine compute_covariances(RLA_jndx, RLO_jndx, Orog_jndx,  SNOforc_jndx,  &
        Lat_Obs, Lon_Obs, Ele_Obs, num_Obs,    				&
        Stdev_back, Stdev_Obs_depth, Stdev_Obs_ims,         &
        L_horz, h_ver,                                      &   !L_horz in Km, h_ver in m
        assim_IMS,                                          &
        B_cov_mat, b_cov_vect, O_cov_mat, W_wght_vect) !,      &LENSFC,
        
        IMPLICIT NONE
        !USE intrinsic::ieee_arithmetic
        integer, parameter :: dp = kind(1.d0)

        Real, Intent(In) 	:: RLA_jndx, RLO_jndx, Orog_jndx, SNOforc_jndx
        Real, Intent(In) 	:: Stdev_back, Stdev_Obs_depth, Stdev_Obs_ims 
        Real, Intent(In) 	:: Lon_Obs(num_Obs), Lat_Obs(num_Obs), Ele_Obs(num_Obs)
        Real, Intent(In) 	:: L_horz, h_ver  !L_horz in Km, h_ver in m
        Integer, Intent(In) :: num_Obs
        LOGICAL, Intent(In) :: assim_IMS

        Real(dp), Intent(Out)    :: B_cov_mat(num_obs,num_obs), b_cov_vect(num_obs)
        Real(dp), Intent(Out)    :: O_cov_mat(num_obs,num_obs), W_wght_vect(num_obs) !6.10.20: W_wght_vect(num_obs, 1)
        
        Real(dp)    :: W_wght_vect_intr(1, num_obs)	
        
        Integer :: indx, jndx, zndx    
        Real    :: rjk_distArr(num_Obs, num_Obs), zjk_distArr(num_Obs, num_Obs)    
        Real    :: l_distArr(num_Obs), h_distArr(num_Obs), haversinArr(num_Obs)
        Real(dp)    :: Innov_cov_mat(num_obs,num_obs), Innov_cov_mat_inv(num_obs,num_obs)
        Real 	:: d_latArr(num_Obs), d_lonArr(num_Obs)
        Real    ::  Lon_Obs_2(num_Obs)      !RLO_2_jndx,
        Real 	:: RLA_rad_jndx, RLO_rad_jndx
        Real 	:: Lat_Obs_rad(num_Obs), Lon_Obs_rad(num_Obs)	
        Real(16), Parameter :: PI_16 = 4 * atan (1.0_16)	
        Real(16), Parameter :: pi_div_180 = PI_16/180.0
        Real, Parameter		:: earth_rad = 6371.
        Real, Parameter		:: snforc_tol = 0.001
        Real, Parameter		:: Bcov_scale_factor = 0.1
        ! PRINT*, "PI: ", PI_16
        ! PRINT*, "PI / 180: ", pi_div_180	
        
        !Lon between -180 and 180 for some inputs
        Lon_Obs_2 = Lon_Obs
        Where(Lon_Obs_2 < 0) Lon_Obs_2 = 360. + Lon_Obs_2
        ! Do zndx = 1, num_Obs 
        !     if (Lon_Obs(zndx) < 0) Lon_Obs_2(zndx) = 360. + Lon_Obs(zndx)
        ! end do   
        ! deg to rad	
        RLA_rad_jndx =  pi_div_180 * RLA_jndx
        RLO_rad_jndx =  pi_div_180 * RLO_jndx
        Lat_Obs_rad =  pi_div_180 * Lat_Obs
        Lon_Obs_rad =  pi_div_180 * Lon_Obs_2   

        !1. Observation covariance 
        ! O = stdev_o*stdev_o * I , I = Identitity matrix
        O_cov_mat = 0.
        Do indx = 1, num_Obs
            O_cov_mat(indx, indx) = Stdev_Obs_depth * Stdev_Obs_depth
        end do
        if (assim_IMS) O_cov_mat(num_Obs, num_Obs) = Stdev_Obs_ims * Stdev_Obs_ims

        if (print_debug) then
            print*, "Obs cov"
            print*, O_cov_mat
        endif

        !2. Background states covariance (at observation space)
        ! B = stddev_back*stdev_back * Corr(j, k)
        ! Corr(j, k) = (1+rjk/L)exp(-rjk/L)exp(-(zjk/h)^2)
        ! rjk = horizontal distance between j and k
        ! zjk = vertical distance between j and k
        ! L = horizontal correlation length (note in Km)
        ! h = vertical correlation length   (in m)
        
        ! https://en.wikipedia.org/wiki/Haversine_formula
        ! dist = 2 * R * asin { sqrt [sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2]}
! 4.16.20 ToDO: This is a symmetric matrix: can revise later to doing only half of the computations
        Do jndx = 1, num_Obs 
            d_latArr = (Lat_Obs_rad(jndx) - Lat_Obs_rad) / 2.
            d_lonArr = (Lon_Obs_rad(jndx) - Lon_Obs_rad) / 2.
            haversinArr = sin(d_latArr)**2 + cos(Lat_Obs_rad) * cos(Lat_Obs_rad(jndx)) * sin(d_lonArr)**2
            Where (haversinArr > 1) haversinArr = 1
            ! Do indx = 1, num_Obs 
            !     if (haversinArr(indx) > 1) haversinArr(indx) = 1 ! ensure numerical errors don't make h>1
            ! end do
            rjk_distArr(jndx,:) = 2 * earth_rad * asin(sqrt(haversinArr))	! rjk, k = 1, Num obs for a given j
            zjk_distArr(jndx,:) = Ele_Obs(jndx) - Ele_Obs       ! zjk, k = 1, Num obs for a given j
        End do
        !Corr_j_k = (1+rjk/L)exp(-rjk/L)exp(-(zjk/h)^2)
        if (print_debug) then
            print*, "Dist for Back corr at obs pts"
            print*, rjk_distArr
            print*, "Vertical dist for Back corr at obs pts"
            print*, zjk_distArr
        endif
        B_cov_mat = (1. + rjk_distArr/L_horz) * exp(-1. * rjk_distArr/L_horz) !L_horz in Km, h_ver in m
        B_cov_mat = B_cov_mat * exp(-1. * (zjk_distArr/h_ver)**2)
        if (print_debug) then
            print*, "Backround corr at obs pts"
            print*, B_cov_mat
        endif   
        !B_cov_mat(num_Obs, num_Obs) = 1.
        B_cov_mat = B_cov_mat * Stdev_back * stdev_back   
        if (print_debug) then
            print*, "Backround cov at obs pts"
            print*, B_cov_mat
        endif   
        
        !3. b background covariance between model grid and obs points
        ! similar to the above (B_cov_mat) except just for one point againt N obs points 
        d_latArr = (RLA_rad_jndx - Lat_Obs_rad) / 2.
        d_lonArr = (RLO_rad_jndx - Lon_Obs_rad) / 2.
        haversinArr = sin(d_latArr)**2 + cos(Lat_Obs_rad) * cos(RLA_rad_jndx) * sin(d_lonArr)**2
        Where (haversinArr > 1) haversinArr = 1.
        ! Do indx = 1, num_Obs 
        !     if (haversinArr(indx) > 1) haversinArr(indx) = 1 ! ensure numerical errors don't make h>1
        ! end do
        l_distArr = 2 * earth_rad * asin(sqrt(haversinArr))	! rjk, k = 1, Num obs for a given j
        h_distArr = Orog_jndx - Ele_Obs       ! zjk, k = 1, Num obs for a given j
        if (print_debug) then
            print*, "Horz Dist for Back corr at obs pts and model grid"
            print*, l_distArr
            print*, "Vertical dist for Back corr at obs pts and model grid"
            print*, h_distArr
        endif
        !Corr_j_k = (1+rjk/L)exp(-rjk/L)exp(-(zjk/h)^2)
        b_cov_vect =  (1. + l_distArr/L_horz) * exp(-1. * l_distArr/L_horz)  !L_horz in Km, h_ver in m
        b_cov_vect = b_cov_vect * exp(-1. * (h_distArr/h_ver)**2)
        if (print_debug) then
            print*, "b corr between model grid and obs pts"
            print*, b_cov_vect
        endif
        b_cov_vect = b_cov_vect * stdev_back * stdev_back 
        if (print_debug) then
            print*, "b cov between model grid and obs pts"
            print*, b_cov_vect
        endif

        ! 4. Weight vector
        ! W = (B_cov_mat + Obs_cov_mat)^-1 b
        Innov_cov_mat = B_cov_mat + O_cov_mat       
        Innov_cov_mat_inv = inv(Innov_cov_mat)        
        !W_wght_vect = matmul(Innov_cov_mat_inv, RESHAPE(b_cov_vect,(/num_obs,1/)))
        W_wght_vect_intr = matmul(RESHAPE(b_cov_vect,(/1, num_obs/)), Innov_cov_mat_inv) ! [1,m]x[m,m]=[1,m]
        W_wght_vect = RESHAPE(W_wght_vect_intr,(/num_obs/))
        if (print_debug) then
            print*, "Innov cov"
            print*, Innov_cov_mat
            print*, "inverse of Innov cov"
            print*, Innov_cov_mat_inv
            print*, "weights vector"
            print*, W_wght_vect
        endif
        
        RETURN

    END SUBROUTINE compute_covariances   

    SUBROUTINE Snow_DA_OI(back_at_Obs, obs_Array, num_Obs,  &
        W_wght_vect,            &
        back_at_Grid, innov_at_Grid, anl_at_Grid, obs_Innov)
    
        IMPLICIT NONE

        include 'mpif.h'

        integer, parameter :: dp = kind(1.d0)

        REAL, INTENT(In)    :: back_at_Obs(num_Obs), obs_Array(num_Obs)
        Integer, Intent(In)  :: num_Obs
        Real(dp), Intent(In)    :: W_wght_vect(num_Obs)    !, weighted_Innov
        REAL, INTENT(In)    :: back_at_Grid
        REAL, INTENT(Out)    :: innov_at_Grid, anl_at_Grid
        Real, INTENT(Out)  :: obs_Innov(num_Obs)
        Real               :: innov_at_Grid_interm(1,1)

        obs_Innov = obs_Array - back_at_Obs
        innov_at_Grid_interm = matmul(RESHAPE(W_wght_vect,(/1,num_Obs/)), &
                                    RESHAPE(obs_Innov,(/num_Obs,1/))) ! weighted innovation
        
        innov_at_Grid = innov_at_Grid_interm(1,1)
        anl_at_Grid = back_at_Grid + innov_at_Grid

        RETURN

    END SUBROUTINE Snow_DA_OI

    SUBROUTINE Obs_Quality_Check(obs_Array, num_Obs,  &
        Obs_QC, num_obs_qc)
    
        IMPLICIT NONE

        include 'mpif.h'

        integer, parameter :: dp = kind(1.d0)

        REAL, INTENT(In)    :: obs_Array(num_Obs)
        Integer, Intent(In)  :: num_Obs
        REAL, Allocatable, INTENT(Out)    :: Obs_QC(:)
        Integer, INTENT(Out)  :: num_obs_qc

        

        RETURN

    END SUBROUTINE Obs_Quality_Check

    SUBROUTINE Anl_Quality_Check(obs_Array, num_Obs,  &
        Obs_QC, num_obs_qc)
    
        IMPLICIT NONE

        include 'mpif.h'

        integer, parameter :: dp = kind(1.d0)

        REAL, INTENT(In)    :: obs_Array(num_Obs)
        Integer, Intent(In)  :: num_Obs
        REAL, Allocatable, INTENT(Out)    :: Obs_QC(:)
        Integer, INTENT(Out)  :: num_obs_qc

        

        RETURN

    END SUBROUTINE Anl_Quality_Check
    
    !http://fortranwiki.org/fortran/show/Matrix+inversion
    ! ! Returns the inverse of a matrix calculated by finding the LU
    ! decomposition.  Depends on LAPACK.
    function inv(A) result(Ainv)
        
        IMPLICIT NONE
        
        integer, parameter :: dp = kind(1.d0)

        real(dp), dimension(:,:), intent(in) :: A
        real(dp), dimension(size(A,1),size(A,2)) :: Ainv
    
        real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
        integer, dimension(size(A,1)) :: ipiv   ! pivot indices
        integer :: n, info
    
        ! External procedures defined in LAPACK
        external DGETRF
        external DGETRI
    
        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = size(A,1)
    
        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call DGETRF(n, n, Ainv, n, ipiv, info)
    
        if (info /= 0) then
        stop 'Matrix is numerically singular!'
        end if
    
        ! DGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.
        call DGETRI(n, Ainv, n, ipiv, work, n, info)
    
        if (info /= 0) then
        stop 'Matrix inversion failed!'
        end if
    end function inv

END MODULE M_DA_OI
     
    