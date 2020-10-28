
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

MODULE M_DA

    USE NETCDF
    Use, Intrinsic :: IEEE_ARITHMETIC
    use, Intrinsic :: iso_c_binding
    use IFPORT
    
    Logical, parameter :: print_debug = .False.

    Logical, parameter :: print_deb = .False.  ! change this to /= 1 to turn off print
    
    Integer(2), External :: compar
    
    CONTAINS
       
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
    
    Subroutine snow_DA_EnSRF(RLA_jndx, RLO_jndx, Orog_jndx,   &
        Lat_Obs, Lon_Obs, Ele_Obs,     				&
        L_horz, h_ver,                                      &   !L_horz in Km, h_ver in m
        assim_IMS, rcov_localize,                     &
        jindx, ens_size, LENSFC, SNOFCS_Inp_Ens,          &
        num_Obs_1, num_Obs, loc_nearest_Obs, SNOFCS_atObs_ens,	 &
        Stdev_Obs_depth, Stdev_Obs_ims,                  &
        obs_Array,                          &
        obs_Innov, innov_atGrid_ensM, anl_at_Grid_ens)
        
        IMPLICIT NONE
        !USE intrinsic::ieee_arithmetic
        integer, parameter :: dp = kind(1.d0)

        Real, Intent(In) 	:: RLA_jndx, RLO_jndx, Orog_jndx
        Real, Intent(In) 	:: Lat_Obs(num_Obs), Lon_Obs(num_Obs), Ele_Obs(num_Obs)
        Real, Intent(In) 	:: L_horz, h_ver  !L_horz in Km, h_ver in m

        LOGICAL, Intent(In) :: assim_IMS, rcov_localize
        Integer, Intent(In) :: jindx, ens_size, LENSFC 
        Real, Intent(In) 	:: SNOFCS_Inp_Ens(ens_size, LENSFC)
        Integer, Intent(In) :: num_Obs_1, num_Obs, loc_nearest_Obs(num_Obs_1)        
        Real, Intent(In) 	:: SNOFCS_atObs_ens(ens_size, num_Obs_1)
        Real, Intent(In) 	:: Stdev_Obs_depth, Stdev_Obs_ims
        REAL, INTENT(In)    :: obs_Array(num_Obs)

        Real, INTENT(Out)   :: obs_Innov(num_Obs)
        Real, INTENT(Out)   :: innov_atGrid_ensM
        Real, INTENT(Out)   :: anl_at_Grid_ens(ens_size+1)          

        Real                :: X_state(1, ens_size), Xh_state_atObs(num_Obs, ens_size)
        Real(dp)            :: X_ave_State, Xh_ave_State(num_obs)
        Real(dp)            :: X_ens_Anomaly(1, ens_size), Xh_ens_Anomaly(num_Obs, ens_size)
        Real(dp)            :: Pxz_cov(1, num_obs), Pzz_cov(num_Obs, num_Obs)
        Real(dp)            :: Pzz_s(num_Obs, num_Obs), Pzz_s_i(num_Obs, num_Obs)
        Real(dp)            :: R_cov(num_obs, num_obs), Rs_cov(num_obs, num_obs)
        Real(dp)            :: K_gain(1, num_obs), K_anom_gain(1, num_obs)   !, Pzz_cov_inv(num_obs, num_obs)
        Real(dp)            :: innov_ens_Anomaly(1, ens_size)
        Real(dp)            :: innov_atGrid_ensM_interm(1, 1) 
        Real(dp)            :: Xu_ens_Anomaly_upd(1, ens_size)
        Real(dp)            :: anl_at_Grid_ens_interm(1, ens_size)

        !Integer             :: indx, zndx 
        Integer :: indx, jndx, zndx    
        ! Real    :: rjk_distArr(num_Obs, num_Obs), zjk_distArr(num_Obs, num_Obs)    
        Real    :: l_distArr(num_Obs), h_distArr(num_Obs), haversinArr(num_Obs)
        Real(dp)    :: R_cov_loc(num_obs)  !,num_obs)
        
        Real 	:: d_latArr(num_Obs), d_lonArr(num_Obs)
        Real    ::  Lon_Obs_2(num_Obs)      !RLO_2_jndx,
        Real 	:: RLA_rad_jndx, RLO_rad_jndx
        Real 	:: Lat_Obs_rad(num_Obs), Lon_Obs_rad(num_Obs)	
        Real(16), Parameter :: PI_16 = 4 * atan (1.0_16)	
        Real(16), Parameter :: pi_div_180 = PI_16/180.0
        Real, Parameter		:: earth_rad = 6371.
        
        !Lon between -180 and 180 for some inputs
        Lon_Obs_2 = Lon_Obs
        Where(Lon_Obs_2 < 0) Lon_Obs_2 = 360. + Lon_Obs_2
	
        RLA_rad_jndx =  pi_div_180 * RLA_jndx
        RLO_rad_jndx =  pi_div_180 * RLO_jndx
        Lat_Obs_rad =  pi_div_180 * Lat_Obs
        Lon_Obs_rad =  pi_div_180 * Lon_Obs_2          
        ! https://en.wikipedia.org/wiki/Haversine_formula
        d_latArr = (RLA_rad_jndx - Lat_Obs_rad) / 2.
        d_lonArr = (RLO_rad_jndx - Lon_Obs_rad) / 2.
        haversinArr = sin(d_latArr)**2 + cos(Lat_Obs_rad) * cos(RLA_rad_jndx) * sin(d_lonArr)**2
        Where (haversinArr > 1) haversinArr = 1.
        l_distArr = 2 * earth_rad * asin(sqrt(haversinArr))	! rjk, k = 1, Num obs for a given j
        h_distArr = Orog_jndx - Ele_Obs       ! zjk, k = 1, Num obs for a given j
        if (print_debug) then
            print*, "Horz Dist for Back corr at obs pts and model grid"
            print*, l_distArr
            print*, "Vertical dist for Back corr at obs pts and model grid"
            print*, h_distArr
        endif      
        ! factor for localization of R covariance 
        R_cov_loc = exp(1. * (l_distArr/L_horz)**2) !L_horz in Km, h_ver in m (1. + rjk_distArr/L_horz) * 
        R_cov_loc = R_cov_loc * exp(1. * (h_distArr/h_ver)**2)
        if (print_debug) then
            print*, "Obs cov localization factor"
            print*, R_cov_loc
        endif   
        
        ! Background states X and Xh (Dim [E] [m, E], E = ensemble size, m = num obs)
        X_state(1, :) = SNOFCS_Inp_Ens(:, jindx)
        Do zndx = 1, num_Obs_1
            Xh_state_atObs(zndx, :) = SNOFCS_atObs_ens(:, loc_nearest_Obs(zndx))
        End do
        if(assim_IMS) Xh_state_atObs(num_Obs, :) = SNOFCS_Inp_Ens(:, jindx)
        !SUM(SNWDEN, Mask = (.not. IEEE_IS_NAN(SNWDEN))) / COUNT (.not. IEEE_IS_NAN(SNWDEN))
        X_ave_State = SUM(X_state(1,:)) / ens_size
        Do zndx = 1, num_Obs
            Xh_ave_State(zndx) = SUM(Xh_state_atObs(zndx, :)) / ens_size
        End do  
        ! ens anomaly X' and Xh'    
        X_ens_Anomaly(1,:) = X_state(1,:) - X_ave_State
        Do zndx = 1, num_Obs
            Xh_ens_Anomaly(zndx, :) = Xh_state_atObs(zndx, :) - Xh_ave_State(zndx)
        End do
        ! State-obs xCov: Pxz = X'*Xh'_T [1,m] (_T = transpose)
        Pxz_cov = matmul(X_ens_Anomaly, TRANSPOSE(Xh_ens_Anomaly)) ![1,E][E,m] = [1,m]
        Pxz_cov = Pxz_cov / (ens_size - 1)
        ! Observation covariance 
        ! R = stdev_o*stdev_o * I , I = Identitity matrix
        R_cov = 0.
        Do indx = 1, num_Obs
            R_cov(indx, indx) = Stdev_Obs_depth * Stdev_Obs_depth
        end do
        if (assim_IMS) R_cov(num_Obs, num_Obs) = Stdev_Obs_ims * Stdev_Obs_ims
        if (print_debug) then
            print*, "Obs cov before localization"
            print*, R_cov
        endif
        ! R cov localization 
        ! R_cov = R_cov * R_cov_loc
        if(rcov_localize) then 
            Do indx = 1, num_Obs
                R_cov(indx, indx) = R_cov(indx, indx) * R_cov_loc(indx)
            end do
            if (print_debug) then
                print*, "Obs cov after localization"
                print*, R_cov
            endif
        endif  
        !sqrt of R
        Rs_cov = matsqrt(R_cov)
        ! State at obs Cov: Pzz_cov = Xh'*Xh'_T [m,m] 
        Pzz_cov = matmul(Xh_ens_Anomaly, TRANSPOSE(Xh_ens_Anomaly)) ![m,E][E,m] = [m,m]
        Pzz_cov = Pzz_cov / (ens_size - 1)        
        ! Innovvation Cov: Pzz = Pzz + R = Xh'*Xh'_T + R [m,m] 
        Pzz_cov = Pzz_cov + R_cov
        ! sqrt pzz
        Pzz_s = matsqrt(Pzz_cov)
        ! inv
        Pzz_s_i = inv(Pzz_s)                       
        ! Kalman gain K = Pxz*Pzz_i [1,m][m,m] = [1,m] !𝐾=  𝑃_𝑥𝑧 𝑃_𝑧𝑧^(−1)
        K_gain = matmul(Pxz_cov, inv(Pzz_cov))    
        ! 𝑲^′=  𝑃_𝑥𝑧 ∗ (𝑃_(𝐳𝐳_𝒔_𝒊) )^𝑻 ∗ (𝑃_(𝐳𝑧_𝒔)+𝑹_𝒔 )^(−𝟏)
        K_anom_gain = matmul(Pxz_cov, TRANSPOSE(Pzz_s_i))
        K_anom_gain = matmul(K_anom_gain, inv(Pzz_s + Rs_cov))       
        if (print_debug) then
            print*, "Innov cov"
            print*, Pzz_cov
            print*, "inverse of Innov cov"
            print*, inv(Pzz_cov)
            print*, "Ens mean Gain"
            print*, K_gain
            print*, "Ens anomaly Gain"
            print*, K_anom_gain
        endif
        ! Innovation at obs d = Z - Xh [m,E]
        Do zndx = 1, num_Obs
            obs_Innov(zndx) = obs_Array(zndx) - Xh_ave_State(zndx)
        End do
        ! Ens mean Innovation at model dX = K*d [1,m][m,1] = [1,1]
        innov_atGrid_ensM_interm = matmul(K_gain, RESHAPE(obs_Innov, (/num_Obs, 1/))) 
        innov_atGrid_ensM = innov_atGrid_ensM_interm(1, 1) 
        ! Ens mean Update 
        anl_at_Grid_ens(ens_size + 1) = X_ave_State + innov_atGrid_ensM
        ! Ens anomaly Innovation = 𝑲′ * 𝑿_𝒉′ [1,m][m,E] = [1,E]
        innov_ens_Anomaly = matmul(K_anom_gain, Xh_ens_Anomaly)
        ! Ens anomaly Update (1,E) 𝑿′𝒂= 𝑿′ − 𝑲′*𝑿_𝒉′ 
        Xu_ens_Anomaly_upd = X_ens_Anomaly - innov_ens_Anomaly
        ! non-ens_mean, ensemble members update
        anl_at_Grid_ens_interm(1,:) = Xu_ens_Anomaly_upd(1, :) + anl_at_Grid_ens(ens_size + 1)
        anl_at_Grid_ens(1:ens_size) = RESHAPE(anl_at_Grid_ens_interm,(/ens_size/))
        if (print_debug) then
            print*, "Obs Innov"
            print*, obs_Innov
            print*, "Innov at grid ens mean"
            print*, innov_atGrid_ensM
            print*, "Innov at grid ens anomaly"
            print*, innov_ens_Anomaly
            print*, "Analysis"
            print*, anl_at_Grid_ens
        endif

        RETURN

    End Subroutine snow_DA_EnSRF

    Subroutine snow_DA_EnKF_LocCov(RLA_jndx, RLO_jndx, Orog_jndx,   &
        Lat_Obs, Lon_Obs, Ele_Obs,     				&
        L_horz, h_ver,                                      &   !L_horz in Km, h_ver in m
        assim_IMS,                    &
        jindx, ens_size, LENSFC, SNOFCS_Inp_Ens,          &
        num_Obs_1, num_Obs, loc_nearest_Obs, SNOFCS_atObs_ens,	 &
        Stdev_Obs_depth, Stdev_Obs_ims,                  &
        obs_Array,                          &
        obs_Innov_ens, innov_at_Grid_ens, anl_at_Grid_ens)
        
        IMPLICIT NONE
        !USE intrinsic::ieee_arithmetic
        integer, parameter :: dp = kind(1.d0)

        Real, Intent(In) 	:: RLA_jndx, RLO_jndx, Orog_jndx
        Real, Intent(In) 	:: Lat_Obs(num_Obs), Lon_Obs(num_Obs), Ele_Obs(num_Obs)
        Real, Intent(In) 	:: L_horz, h_ver  !L_horz in Km, h_ver in m

        LOGICAL, Intent(In) :: assim_IMS
        Integer, Intent(In) :: jindx, ens_size, LENSFC 
        Real, Intent(In) 	:: SNOFCS_Inp_Ens(ens_size, LENSFC)
        Integer, Intent(In) :: num_Obs_1, num_Obs, loc_nearest_Obs(num_Obs_1)        
        Real, Intent(In) 	:: SNOFCS_atObs_ens(ens_size, num_Obs_1)
        Real, Intent(In) 	:: Stdev_Obs_depth, Stdev_Obs_ims
        REAL, INTENT(In)    :: obs_Array(num_Obs)

        Real, INTENT(Out)   :: obs_Innov_ens(num_Obs, ens_size)
        Real, INTENT(Out)   :: innov_at_Grid_ens(ens_size+1)
        Real, INTENT(Out)   :: anl_at_Grid_ens(ens_size+1)          

        Real                :: X_state(1, ens_size), Xh_state_atObs(num_Obs, ens_size)
        Real(dp)            :: X_ave_State, Xh_ave_State(num_obs)
        Real(dp)            :: X_ens_Anomaly(1, ens_size), Xh_ens_Anomaly(num_Obs, ens_size)
        Real(dp)            :: Pxz_cov(1, num_obs), Pzz_h_cov(num_Obs, num_Obs)
        Real(dp)            :: Pzz_cov(num_Obs, num_Obs), R_cov(num_obs, num_obs)
        Real(dp)            :: Pzz_cov_inv(num_obs, num_obs), K_gain(1, num_obs)
        Real(dp)            :: innov_at_Grid_ens_interm(1, ens_size)
        Real(dp)            :: anl_at_Grid_ens_interm(1, ens_size) 

        !Integer             :: indx, zndx 
        Integer :: indx, jndx, zndx    
        ! Real    :: rjk_distArr(num_Obs, num_Obs), zjk_distArr(num_Obs, num_Obs)    
        Real    :: l_distArr(num_Obs), h_distArr(num_Obs), haversinArr(num_Obs)
        Real(dp)    :: R_cov_loc(num_obs)  !,num_obs)
        
        Real 	:: d_latArr(num_Obs), d_lonArr(num_Obs)
        Real    ::  Lon_Obs_2(num_Obs)      !RLO_2_jndx,
        Real 	:: RLA_rad_jndx, RLO_rad_jndx
        Real 	:: Lat_Obs_rad(num_Obs), Lon_Obs_rad(num_Obs)	
        Real(16), Parameter :: PI_16 = 4 * atan (1.0_16)	
        Real(16), Parameter :: pi_div_180 = PI_16/180.0
        Real, Parameter		:: earth_rad = 6371.
        
        !Lon between -180 and 180 for some inputs
        Lon_Obs_2 = Lon_Obs
        Where(Lon_Obs_2 < 0) Lon_Obs_2 = 360. + Lon_Obs_2
	
        RLA_rad_jndx =  pi_div_180 * RLA_jndx
        RLO_rad_jndx =  pi_div_180 * RLO_jndx
        Lat_Obs_rad =  pi_div_180 * Lat_Obs
        Lon_Obs_rad =  pi_div_180 * Lon_Obs_2          
        ! https://en.wikipedia.org/wiki/Haversine_formula
        ! Do jndx = 1, num_Obs 
        !     d_latArr = (Lat_Obs_rad(jndx) - Lat_Obs_rad) / 2.
        !     d_lonArr = (Lon_Obs_rad(jndx) - Lon_Obs_rad) / 2.
        !     haversinArr = sin(d_latArr)**2 + cos(Lat_Obs_rad) * cos(Lat_Obs_rad(jndx)) * sin(d_lonArr)**2
        !     Where (haversinArr > 1) haversinArr = 1
        !     rjk_distArr(jndx,:) = 2 * earth_rad * asin(sqrt(haversinArr))	! rjk, k = 1, Num obs for a given j
        !     zjk_distArr(jndx,:) = Ele_Obs(jndx) - Ele_Obs       ! zjk, k = 1, Num obs for a given j
        ! End do
        !Corr_j_k = (1+rjk/L)exp(-rjk/L)exp(-(zjk/h)^2)
        ! if (print_debug) then
        !     print*, "Dist for Back corr at obs pts"
        !     print*, rjk_distArr
        !     print*, "Vertical dist for Back corr at obs pts"
        !     print*, zjk_distArr
        ! endif
        ! factor for localization of R covariance 
        ! R_cov_loc = exp(1. * (rjk_distArr/L_horz)**2) !L_horz in Km, h_ver in m (1. + rjk_distArr/L_horz) * 
        ! R_cov_loc = R_cov_loc * exp(1. * (zjk_distArr/h_ver)**2)
        d_latArr = (RLA_rad_jndx - Lat_Obs_rad) / 2.
        d_lonArr = (RLO_rad_jndx - Lon_Obs_rad) / 2.
        haversinArr = sin(d_latArr)**2 + cos(Lat_Obs_rad) * cos(RLA_rad_jndx) * sin(d_lonArr)**2
        Where (haversinArr > 1) haversinArr = 1.
        l_distArr = 2 * earth_rad * asin(sqrt(haversinArr))	! rjk, k = 1, Num obs for a given j
        h_distArr = Orog_jndx - Ele_Obs       ! zjk, k = 1, Num obs for a given j
        if (print_debug) then
            print*, "Horz Dist for Back corr at obs pts and model grid"
            print*, l_distArr
            print*, "Vertical dist for Back corr at obs pts and model grid"
            print*, h_distArr
        endif      
        ! factor for localization of R covariance 
        R_cov_loc = exp(1. * (l_distArr/L_horz)**2) !L_horz in Km, h_ver in m (1. + rjk_distArr/L_horz) * 
        R_cov_loc = R_cov_loc * exp(1. * (h_distArr/h_ver)**2)
        if (print_debug) then
            print*, "Obs cov localization factor"
            print*, R_cov_loc
        endif   
        
        ! Background states X and Xh (Dim [E] [m, E], E = ensemble size, m = num obs)
        X_state(1, :) = SNOFCS_Inp_Ens(:, jindx)
        Do zndx = 1, num_Obs_1
            Xh_state_atObs(zndx, :) = SNOFCS_atObs_ens(:, loc_nearest_Obs(zndx))
        End do
        if(assim_IMS) Xh_state_atObs(num_Obs, :) = SNOFCS_Inp_Ens(:, jindx)
        !SUM(SNWDEN, Mask = (.not. IEEE_IS_NAN(SNWDEN))) / COUNT (.not. IEEE_IS_NAN(SNWDEN))
        X_ave_State = SUM(X_state(1,:)) / ens_size
        Do zndx = 1, num_Obs
            Xh_ave_State(zndx) = SUM(Xh_state_atObs(zndx, :)) / ens_size
        End do  
        ! ens anomaly X' and Xh'    
        X_ens_Anomaly(1,:) = X_state(1,:) - X_ave_State
        Do zndx = 1, num_Obs
            Xh_ens_Anomaly(zndx, :) = Xh_state_atObs(zndx, :) - Xh_ave_State(zndx)
        End do
        ! State-obs xCov: Pxz = X'*Xh'_T [1,m] (_T = transpose)
        Pxz_cov = matmul(X_ens_Anomaly, TRANSPOSE(Xh_ens_Anomaly)) ![1,E][E,m] = [1,m]
        Pxz_cov = Pxz_cov / (ens_size - 1)
        ! State at obs Cov: Pzz_h = Xh'*Xh'_T [m,m] 
        Pzz_h_cov = matmul(Xh_ens_Anomaly, TRANSPOSE(Xh_ens_Anomaly)) ![m,E][E,m] = [m,m]
        Pzz_h_cov = Pzz_h_cov / (ens_size - 1)
        ! Observation covariance 
        ! R = stdev_o*stdev_o * I , I = Identitity matrix
        R_cov = 0.
        Do indx = 1, num_Obs
            R_cov(indx, indx) = Stdev_Obs_depth * Stdev_Obs_depth
        end do
        if (assim_IMS) R_cov(num_Obs, num_Obs) = Stdev_Obs_ims * Stdev_Obs_ims
        if (print_debug) then
            print*, "Obs cov before localization"
            print*, R_cov
        endif
        ! R cov localization 
        ! R_cov = R_cov * R_cov_loc
        Do indx = 1, num_Obs
            R_cov(indx, indx) = R_cov(indx, indx) * R_cov_loc(indx)
        end do
        if (print_debug) then
            print*, "Obs cov after localization"
            print*, R_cov
        endif
        ! Innovvation Cov: Pzz = Pzz_h + R = Xh'*Xh'_T + R [m,m] 
        Pzz_cov = Pzz_h_cov + R_cov
        Pzz_cov_inv = inv(Pzz_cov)                       
        ! Kalman gain K = Pxz*Pzz_i [1,m][m,m] = [1,m]
        K_gain = matmul(Pxz_cov, Pzz_cov_inv)        
        if (print_debug) then
            print*, "Innov cov"
            print*, Pzz_cov
            print*, "inverse of Innov cov"
            print*, Pzz_cov_inv
            print*, "Gain vector"
            print*, K_gain
        endif
        ! Innovation at obs d = Z - Xh [m,E]
        Do zndx = 1, num_Obs
            obs_Innov_ens(zndx, :) = obs_Array(zndx) - Xh_state_atObs(zndx, :)
        End do
        ! Innovation at model dX = K*d [1,m][m,E] = [1,E]
        innov_at_Grid_ens_interm = matmul(K_gain, obs_Innov_ens)
        ! Update (1,E)
        anl_at_Grid_ens_interm = X_state + innov_at_Grid_ens_interm
        !
        innov_at_Grid_ens(1:ens_size) = RESHAPE(innov_at_Grid_ens_interm,(/ens_size/))
        anl_at_Grid_ens(1:ens_size) = RESHAPE(anl_at_Grid_ens_interm,(/ens_size/))
        ! ens mean
        anl_at_Grid_ens(ens_size + 1) = SUM(anl_at_Grid_ens(1:ens_size)) / ens_size
        innov_at_Grid_ens(ens_size + 1) = SUM(innov_at_Grid_ens(1:ens_size)) / ens_size
        if (print_debug) then
            print*, "Obs Innov"
            print*, obs_Innov_ens
            print*, "Innov at grid"
            print*, innov_at_Grid_ens
            print*, "Analysis"
            print*, anl_at_Grid_ens
        endif

        RETURN

    End Subroutine snow_DA_EnKF_LocCov

    Subroutine snow_DA_EnKF(assim_IMS,                    &
        jindx, ens_size, LENSFC, SNOFCS_Inp_Ens,          &
        num_Obs_1, num_Obs, loc_nearest_Obs, SNOFCS_atObs_ens,	 &
        Stdev_Obs_depth, Stdev_Obs_ims,                  &
        obs_Array,                          &
        obs_Innov_ens, innov_at_Grid_ens, anl_at_Grid_ens)
        
        IMPLICIT NONE
        !USE intrinsic::ieee_arithmetic
        integer, parameter :: dp = kind(1.d0)

        LOGICAL, Intent(In) :: assim_IMS
        Integer, Intent(In) :: jindx, ens_size, LENSFC 
        Real, Intent(In) 	:: SNOFCS_Inp_Ens(ens_size, LENSFC)
        Integer, Intent(In) :: num_Obs_1, num_Obs, loc_nearest_Obs(num_Obs_1)        
        Real, Intent(In) 	:: SNOFCS_atObs_ens(ens_size, num_Obs_1)
        Real, Intent(In) 	:: Stdev_Obs_depth, Stdev_Obs_ims
        REAL, INTENT(In)    :: obs_Array(num_Obs)

        Real, INTENT(Out)   :: obs_Innov_ens(num_Obs, ens_size)
        Real, INTENT(Out)   :: innov_at_Grid_ens(ens_size+1)
        Real, INTENT(Out)   :: anl_at_Grid_ens(ens_size+1)          

        Real                :: X_state(1, ens_size), Xh_state_atObs(num_Obs, ens_size)
        Real(dp)            :: X_ave_State, Xh_ave_State(num_obs)
        Real(dp)            :: X_ens_Anomaly(1, ens_size), Xh_ens_Anomaly(num_Obs, ens_size)
        Real(dp)            :: Pxz_cov(1, num_obs), Pzz_h_cov(num_Obs, num_Obs)
        Real(dp)            :: Pzz_cov(num_Obs, num_Obs), R_cov(num_obs, num_obs)
        Real(dp)            :: Pzz_cov_inv(num_obs, num_obs), K_gain(1, num_obs)
        Real(dp)            :: innov_at_Grid_ens_interm(1, ens_size)
        Real(dp)            :: anl_at_Grid_ens_interm(1, ens_size) 

        Integer             :: indx, zndx 
        
        ! Background states X and Xh (Dim [E] [m, E], E = ensemble size, m = num obs)
        X_state(1, :) = SNOFCS_Inp_Ens(:, jindx)
        Do zndx = 1, num_Obs_1
            Xh_state_atObs(zndx, :) = SNOFCS_atObs_ens(:, loc_nearest_Obs(zndx))
        End do
        if(assim_IMS) Xh_state_atObs(num_Obs, :) = SNOFCS_Inp_Ens(:, jindx)
        !SUM(SNWDEN, Mask = (.not. IEEE_IS_NAN(SNWDEN))) / COUNT (.not. IEEE_IS_NAN(SNWDEN))
        X_ave_State = SUM(X_state(1,:)) / ens_size
        Do zndx = 1, num_Obs
            Xh_ave_State(zndx) = SUM(Xh_state_atObs(zndx, :)) / ens_size
        End do  
        ! ens anomaly X' and Xh'    
        X_ens_Anomaly(1,:) = X_state(1,:) - X_ave_State
        Do zndx = 1, num_Obs
            Xh_ens_Anomaly(zndx, :) = Xh_state_atObs(zndx, :) - Xh_ave_State(zndx)
        End do
        ! State-obs xCov: Pxz = X'*Xh'_T [1,m] (_T = transpose)
        Pxz_cov = matmul(X_ens_Anomaly, TRANSPOSE(Xh_ens_Anomaly)) ![1,E][E,m] = [1,m]
        Pxz_cov = Pxz_cov / (ens_size - 1)
        ! State at obs Cov: Pzz_h = Xh'*Xh'_T [m,m] 
        Pzz_h_cov = matmul(Xh_ens_Anomaly, TRANSPOSE(Xh_ens_Anomaly)) ![m,E][E,m] = [m,m]
        Pzz_h_cov = Pzz_h_cov / (ens_size - 1)
        ! Observation covariance 
        ! R = stdev_o*stdev_o * I , I = Identitity matrix
        R_cov = 0.
        Do indx = 1, num_Obs
            R_cov(indx, indx) = Stdev_Obs_depth * Stdev_Obs_depth
        end do
        if (assim_IMS) R_cov(num_Obs, num_Obs) = Stdev_Obs_ims * Stdev_Obs_ims
        if (print_debug) then
            print*, "Obs cov"
            print*, R_cov
        endif
        ! Innovvation Cov: Pzz = Pzz_h + R = Xh'*Xh'_T + R [m,m] 
        Pzz_cov = Pzz_h_cov + R_cov
        Pzz_cov_inv = inv(Pzz_cov)                       
        ! Kalman gain K = Pxz*Pzz_i [1,m][m,m] = [1,m]
        K_gain = matmul(Pxz_cov, Pzz_cov_inv)        
        if (print_debug) then
            print*, "Innov cov"
            print*, Pzz_cov
            print*, "inverse of Innov cov"
            print*, Pzz_cov_inv
            print*, "Gain vector"
            print*, K_gain
        endif
        ! Innovation at obs d = Z - Xh [m,E]
        Do zndx = 1, num_Obs
            obs_Innov_ens(zndx, :) = obs_Array(zndx) - Xh_state_atObs(zndx, :)
        End do
        ! Innovation at model dX = K*d [1,m][m,E] = [1,E]
        innov_at_Grid_ens_interm = matmul(K_gain, obs_Innov_ens)
        ! Update (1,E)
        anl_at_Grid_ens_interm = X_state + innov_at_Grid_ens_interm
        !
        innov_at_Grid_ens(1:ens_size) = RESHAPE(innov_at_Grid_ens_interm,(/ens_size/))
        anl_at_Grid_ens(1:ens_size) = RESHAPE(anl_at_Grid_ens_interm,(/ens_size/))
        ! ens mean
        anl_at_Grid_ens(ens_size + 1) = SUM(anl_at_Grid_ens(1:ens_size)) / ens_size
        innov_at_Grid_ens(ens_size + 1) = SUM(innov_at_Grid_ens(1:ens_size)) / ens_size
        if (print_debug) then
            print*, "Obs Innov"
            print*, obs_Innov_ens
            print*, "Innov at grid"
            print*, innov_at_Grid_ens
            print*, "Analysis"
            print*, anl_at_Grid_ens
        endif

        RETURN

    End Subroutine snow_DA_EnKF

    ! ! LAPACK.
    function matsqrt(A) result(Asqrt)
        
        IMPLICIT NONE
        
        integer, parameter :: dp = kind(1.d0)

        real(dp), dimension(:,:), intent(in) :: A
        real(dp), dimension(size(A,1),size(A,2)) :: Asqrt
        integer :: N, info
     
        ! External procedures defined in LAPACK
        external DPOTRF
    
        ! Backup
        Asqrt = A
        N = size(A,1)
    
        call DPOTRF ('L', N, Asqrt, N, info)
    
        if (info /= 0) then
            stop 'Matrix factorization failed!'
        end if

    end function matsqrt

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
                            SNOFCS_atObs, Ele_atObs, index_back_atObs, index_back_atEval, Obs_atEvalPts,      &
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
        Integer, Intent(Out)	:: index_back_atObs(num_Obs)   ! the location of background corresponding obs
        Real, Intent(Out) 	    :: SNOFCS_atEvalPts(num_Eval), Lat_atEvalPts(num_Eval), Lon_atEvalPts(num_Eval)
        
        !Integer	:: index_back_atObs(num_Obs)   ! the location of background corresponding obs
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
            if(num_Eval > 0) then
                Call random_number(rand_nextVal)
                rand_evalPoint(1) = floor(rand_nextVal * num_Obs) + 1
                index_back_atEval(1) = index_back_atObs(rand_evalPoint(1))
                if(.NOT. IEEE_IS_NAN(SNWD_GHCND(rand_evalPoint(1))/   & 
                                     SNWD_GHCND(rand_evalPoint(1)))) then
                    index_back_atEval(1) = index_back_atObs(rand_evalPoint(1))
                else   ! try one more time
                    Call random_number(rand_nextVal)
                    rand_evalPoint(1) = floor(rand_nextVal * num_Obs) + 1
                    index_back_atEval(1) = index_back_atObs(rand_evalPoint(1))
                endif
            end if
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
        
        RETURN
        
     END SUBROUTINE Observation_Operator_Parallel

    ! SUBROUTINE Observation_Operator_Snofcs_Parallel(Myrank, MAX_TASKS, p_tN, p_tRank, Np_til, & 
    !                         RLA, RLO, OROG, Lat_Obs, Lon_Obs,               &
    !                         LENSFC, num_Obs, num_Eval, max_distance, SNOFCS_back, SNWD_GHCND,  &
    !                         Ele_atObs, index_back_atObs, index_back_atEval, Obs_atEvalPts,      &
    !                         SNOFCS_atEvalPts, Lat_atEvalPts, Lon_atEvalPts)
                            
    !     IMPLICIT NONE
    !     !
    !     !USE intrinsic::ieee_arithmetic
    !     include "mpif.h"
    
    !     Real, Intent(In) 	:: RLA(LENSFC), RLO(LENSFC), OROG(LENSFC)
    !     Real, Intent(In) 	:: Lat_Obs(num_Obs), Lon_Obs(num_Obs)  ! don't want to alter these
    !     INTEGER             :: Myrank, MAX_TASKS, p_tN, p_tRank, Np_til, LENSFC, num_Obs, num_Eval
    !     Real	            :: max_distance   ! radius_of_influence
    !     Real, Intent(In) 	:: SNOFCS_back(LENSFC)
    !     Real, Intent(InOut) 	:: SNWD_GHCND(num_Obs) 
    !     Real, Intent(Out) 	    :: Ele_atObs(num_Obs), Obs_atEvalPts(num_Eval)
    !     Integer, Intent(Out) 	:: index_back_atEval(num_Eval)   ! the location of evaluation points
    !     Integer, Intent(Out)	:: index_back_atObs(num_Obs)   ! the location of background corresponding obs
    !     Real, Intent(Out) 	    :: SNOFCS_atEvalPts(num_Eval)
    !     Real, Intent(Out) 	    :: Lat_atEvalPts(num_Eval), Lon_atEvalPts(num_Eval)
        
    !     !Integer	:: index_back_atObs(num_Obs)   ! the location of background corresponding obs
    !     Real 	:: Lon_Obs_2(num_Obs)		!RLO_2(LENSFC), 	
    !     Real 	:: RLA_rad(LENSFC), RLO_rad(LENSFC)
    !     Real 	:: Lat_Obs_rad(num_Obs), Lon_Obs_rad(num_Obs)	
    !     INTEGER :: indx, jndx, jzndx, zndx, min_indx
    !     Real    :: distArr(LENSFC), haversinArr(LENSFC)
    !     Real 	:: d_latArr(LENSFC), d_lonArr(LENSFC)
    !     Real(16), Parameter :: PI_16 = 4 * atan (1.0_16)	
    !     Real(16), Parameter :: pi_div_180 = PI_16/180.0
    !     Real, Parameter		:: earth_rad = 6371.
        
    !     ! for mpi par
    !     INTEGER            :: N_sA, N_sA_Ext, mp_start, mp_end 
    !     INTEGER            :: send_proc, rec_proc, rec_stat(MPI_STATUS_SIZE), dest_Aoffset, pindex
    !     INTEGER            :: mpiReal_size, rsize, mpiInt_size, isize, IERR
    !     Real               :: rand_nextVal  ! randomly select evalution points to exclude from DA
    !     Integer            :: rand_evalPoint(num_Eval)    ! randomly select evalution points to exclude from DA
    
    !     !Np_til ! num proc. per tile p_tRank ! proc. rank within tile !p_tN  ! tile for proc.
    !     N_sA = num_Obs / Np_til  ! sub array length per proc
    !     N_sA_Ext = num_Obs - N_sA * Np_til ! extra grid cells
    !     if(p_tRank == 0) then 
    !         mp_start = 1
    !     else
    !         mp_start = p_tRank * N_sA + N_sA_Ext + 1   ! start index of subarray for proc
    !     endif
    !     mp_end = (p_tRank + 1) * N_sA + N_sA_Ext 		! end index of subarray for proc
    
    !     !Fill values to nan (to differentiate those htat don't have val
    !     Ele_atObs = IEEE_VALUE(Ele_atObs, IEEE_QUIET_NAN)	
    !     index_back_atObs = -1   ! when corresponding value doesn't exit	
    !     rand_evalPoint = -1
    !     index_back_atEval = -1
    
    !     ! RLO from 0 to 360 (no -ve lon)
    !     Lon_Obs_2 = Lon_Obs
    !     Where(Lon_Obs_2 < 0) Lon_Obs_2 = 360. + Lon_Obs_2
    
    !     ! at each obs point compute its distance from RLA/RLO pairs 
    !     ! then find the position of the minimum    
    !     ! shortest distance over sphere using great circle distance 	
    !     RLA_rad =  pi_div_180 * RLA
    !     RLO_rad =  pi_div_180 * RLO
    !     Lat_Obs_rad =  pi_div_180 * Lat_Obs
    !     Lon_Obs_rad =  pi_div_180 * Lon_Obs_2   	
    !     ! https://en.wikipedia.org/wiki/Haversine_formula
    !     ! https://www.geeksforgeeks.org/program-distance-two-points-earth/
    !     ! Distance, d = R * arccos[(sin(lat1) * sin(lat2)) + cos(lat1) * cos(lat2) * cos(long2 – long1)]
    !     ! dist = 2 * R * asin { sqrt [sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2]}
    !     Do indx = mp_start, mp_end   !1, num_Obs 
    !         d_latArr = (Lat_Obs_rad(indx) - RLA_rad) / 2.
    !         d_lonArr = (Lon_Obs_rad(indx) - RLO_rad) / 2.
    !         haversinArr = sin(d_latArr)**2 + cos(Lat_Obs_rad(indx)) * cos(RLA_rad) * sin(d_lonArr)**2
    !         WHERE(haversinArr > 1) haversinArr = 1.   ! ensure numerical errors don't make h>1
            
    !         distArr = 2 * earth_rad * asin(sqrt(haversinArr))		
    !         min_indx = MINLOC(distArr, dim = 1)  !, MASK=ieee_is_nan(distArr))
    
    !         if(distArr(min_indx) < max_distance) then
    !             Ele_atObs(indx) = OROG(min_indx)
    !             index_back_atObs(indx) = min_indx
    !         ! else 
    !             ! Print*, " Warning! distance greater than ",max_distance," km ", distArr(min_indx)
    !         endif
    !     end do
    
    ! ! ToDO: Better way to handle this?
    ! ! Real data type size corresponding to mpi
    !     rsize = SIZEOF(max_distance) 
    !     Call MPI_TYPE_SIZE(MPI_REAL, mpiReal_size, IERR) 
    !     If (rsize == 4 ) then 
    !         mpiReal_size = MPI_REAL4
    !     elseif (rsize == 8 ) then 
    !         mpiReal_size = MPI_REAL8
    !     elseif (rsize == 16 ) then 
    !         mpiReal_size = MPI_REAL16
    !     else
    !         PRINT*," Possible mismatch between Fortran Real ", rsize," and Mpi Real", mpiReal_size
    !         Stop
    !     endif
    !     isize = SIZEOF(N_sA) 
    !     Call MPI_TYPE_SIZE(MPI_INTEGER, mpiInt_size, IERR) 
    !     If (isize == 2 ) then 
    !         mpiInt_size = MPI_INTEGER2
    !     elseif (isize == 4 ) then 
    !         mpiInt_size = MPI_INTEGER4
    !     elseif (isize == 8 ) then 
    !         mpiInt_size = MPI_INTEGER8
    !     else
    !         PRINT*," Possible mismatch between Fortran Int ", isize," and Mpi Int ", mpiInt_size
    !         Stop
    !     endif
    
    !     if (MYRANK > (MAX_TASKS - 1) ) then
    !         call MPI_SEND(Ele_atObs(mp_start:mp_end), N_sA, mpiReal_size, p_tN,   &
    !                       MYRANK*100, MPI_COMM_WORLD, IERR)
    !         call MPI_SEND(index_back_atObs(mp_start:mp_end), N_sA, mpiInt_size, p_tN,   &
    !                       MYRANK*1000, MPI_COMM_WORLD, IERR)
    !     else !if (MYRANK == p_tN ) then  
    !         Do pindex =  1, (Np_til - 1)   ! sender proc index within tile group
    !             dest_Aoffset = pindex * N_sA + N_sA_Ext + 1   ! dest array offset
    !             send_proc = MYRANK +  pindex * MAX_TASKS
    !             call MPI_RECV(Ele_atObs(dest_Aoffset:dest_Aoffset+N_sA-1), N_sA, mpiReal_size, send_proc,   &
    !                       send_proc*100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
    !             call MPI_RECV(index_back_atObs(dest_Aoffset:dest_Aoffset+N_sA-1), N_sA, mpiInt_size, send_proc, &
    !                       send_proc*1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
    !         enddo
    !     endif
    ! !ToDO: better way to do this?
    !     ! now share the whole array
    !     if (MYRANK < MAX_TASKS ) then   !if (MYRANK == p_tN ) then 	
    !         ! Print*, "Started selecting obs points"
    !         ! Select obs points to exclude from DA 
    !         Call random_number(rand_nextVal)
    !         rand_evalPoint(1) = floor(rand_nextVal * num_Obs) + 1
    !         index_back_atEval(1) = index_back_atObs(rand_evalPoint(1))
    !         jndx = 2
    !         jzndx = 1
    !         Do  While (jndx <= num_Eval)        !jndx = 2, num_Eval
    !             Call random_number(rand_nextVal)
    !             ! Print*, rand_nextVal
    !             rand_evalPoint(jndx) = floor(rand_nextVal * num_Obs) + 1
    !             if((rand_evalPoint(jndx) /= rand_evalPoint(jndx-1)) .AND.      &
    !                (.NOT. IEEE_IS_NAN(SNWD_GHCND(rand_evalPoint(jndx))/SNWD_GHCND(rand_evalPoint(jndx)))) )  then
    !                 index_back_atEval(jndx) = index_back_atObs(rand_evalPoint(jndx))
    !                 jndx = jndx + 1
    !             else
    !                 cycle
    !             endif
    !             jzndx = jzndx + 1
    !             If (jzndx >= 2*num_Eval) exit 
    !         Enddo 
    !         ! Print*, "Finished selecting obs points"	
    !         Do pindex =  1, (Np_til - 1)   ! receiving proc index within tile group
    !             rec_proc = MYRANK +  pindex * MAX_TASKS
    !             call MPI_SEND(Ele_atObs, num_Obs, mpiReal_size, rec_proc, MYRANK*100, MPI_COMM_WORLD, IERR)
    !             call MPI_SEND(index_back_atObs, num_Obs, mpiInt_size, rec_proc, MYRANK*1000, MPI_COMM_WORLD, IERR)
    !             call MPI_SEND(index_back_atEval, num_Eval, mpiInt_size, rec_proc, MYRANK*10000, MPI_COMM_WORLD, IERR)
    !             call MPI_SEND(rand_evalPoint, num_Eval, mpiInt_size, rec_proc, MYRANK*100000, MPI_COMM_WORLD, IERR)
    !         enddo
    !     else 
    !         call MPI_RECV(Ele_atObs, num_Obs, mpiReal_size, p_tN, p_tN*100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
    !         call MPI_RECV(index_back_atObs, num_Obs, mpiInt_size, p_tN, p_tN*1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
    !         call MPI_RECV(index_back_atEval, num_Eval, mpiInt_size, p_tN, p_tN*10000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
    !         call MPI_RECV(rand_evalPoint, num_Eval, mpiInt_size, p_tN, p_tN*100000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
    !     endif
    !     Obs_atEvalPts = IEEE_VALUE(rand_nextVal, IEEE_QUIET_NAN)
    !     SNOFCS_atEvalPts = IEEE_VALUE(rand_nextVal, IEEE_QUIET_NAN)
    !     Lat_atEvalPts = IEEE_VALUE(rand_nextVal, IEEE_QUIET_NAN)
    !     Lon_atEvalPts = IEEE_VALUE(rand_nextVal, IEEE_QUIET_NAN)
    !     Do  jndx = 1, num_Eval
    !         if (index_back_atEval(jndx) > 0) then
    !             Obs_atEvalPts(jndx) = SNWD_GHCND(rand_evalPoint(jndx))
    !             SNOFCS_atEvalPts(jndx) = SNOFCS_atObs(rand_evalPoint(jndx))
    !             Lat_atEvalPts(jndx) = Lat_Obs(rand_evalPoint(jndx)) 
    !             Lon_atEvalPts(jndx) = Lon_Obs(rand_evalPoint(jndx))
    !             SNWD_GHCND(rand_evalPoint(jndx)) = IEEE_VALUE(rand_nextVal, IEEE_QUIET_NAN) ! exclude point from DA	
    !         endif	
    !     Enddo 	
        
    !     RETURN
        
    !  END SUBROUTINE Observation_Operator_Snofcs_Parallel
    
    SUBROUTINE Observation_Read_GHCND_Tile_excNaN(p_tN, ghcnd_inp_file, dim_name,	&
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
        INTEGER, Intent(In)               :: p_tN
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
        If(p_tN == 3) then
            Do jndx = 1, NDIM_In
                ! Print*, "jndx = ", jndx
                If((Lat_GHCND_In(jndx) > lat_min) .and. (Lat_GHCND_In(jndx) < lat_max) .and. &
                   (((Lon_GHCND_In(jndx) >= lon_min) .and. (Lon_GHCND_In(jndx) <= 180.)) .or. &
                    ((Lon_GHCND_In(jndx) >= -180.) .and. (Lon_GHCND_In(jndx) <= lon_max))) .and. &
                    (SNWD_GHCND_In(jndx) >= 0 )) then
                        !(.NOT. IEEE_IS_NAN(SNWD_GHCND_In(jndx)))) then
                        NDIM = NDIM + 1
                        index_Array(jcounter) = jndx
                        jcounter = jcounter + 1
                Endif
            End do
        Else
            Do jndx = 1, NDIM_In
                ! Print*, "jndx = ", jndx
                If((Lat_GHCND_In(jndx) > lat_min) .and. (Lat_GHCND_In(jndx) < lat_max) .and. &
                    (Lon_GHCND_In(jndx) > lon_min) .and. (Lon_GHCND_In(jndx) < lon_max) .and. &
                    (SNWD_GHCND_In(jndx) >= 0 )) then
                        !(.NOT. IEEE_IS_NAN(SNWD_GHCND_In(jndx)))) then
                        NDIM = NDIM + 1
                        index_Array(jcounter) = jndx
                        jcounter = jcounter + 1
                Endif
            End do
        Endif
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
        
     End SUBROUTINE Observation_Read_GHCND_Tile_excNaN

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

END MODULE M_DA
     
    