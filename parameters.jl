
module Bounds


    ##? Unscaled Initial POints/ Guesses
    x0_us = [60.0;  60.0;   60.0]
    z0_us = [0.0;   60.0;   60.0;   60.0]
    u0_us = [10.0;  0.0;    0.0;    0.0]

    # des0_us = [100.0]

    ##? Lower and Upper scales
    ls_x = [30.0;   30.0;   30.0]
    us_x = [100.0;  100.0;  100.0]

    ls_z = [0.0;    30.0;   30.0;   30.0]
    us_z = [50.0;  100.0;  100.0;  100.0]

    ls_u = [0.0;    0.0;    0.0;    0.0]
    us_u = [50.0;  50.0;  50.0;  10000000]

    # ls_des = [75.0]
    # us_des = [125.0]

    #region-> Generic Code
    Nx = size(x0_us, 1)
    Nz = size(z0_us, 1)
    Nu = size(u0_us, 1)
    # Ndes = size(des0_us, 1)

    ##Scaled Points
    x0 =    (x0_us - ls_x)      ./ (us_x - ls_x)
    z0 =    (z0_us - ls_z)      ./ (us_z - ls_z)
    u0 =    (u0_us - ls_u)      ./ (us_u - ls_u)
    # des0 =  (des0_us - ls_des)  ./ (us_des - ls_des)
    #endregion

    export  x0,      z0,     u0,    #des0, 
            Nx,      Nz,     Nu,    #Ndes,
            ls_x,    ls_z,   ls_u,  #ls_des,
            us_x,    us_z,   us_u  #us_des   

end

module Model_parameters

    ## District Heat Parameters 
    ρ_dh = 1000         
    Cp_dh = 4.18        #kJ/(kg-K)
    q_dh = 10            #m3/hr

    ##Waste Heat Parameters
    ρ_wh = 1000
    Cp_wh = 4.18
    q_wh = 10

    T_dh_ret = 30
    T_dh_minSup = 60

    V_tes = 100
    V_whb = 10
    V_phb = 10

    export  ρ_dh, Cp_dh, q_dh, 
            ρ_wh, Cp_wh, q_wh, 
            T_dh_ret, T_dh_minSup, 
            V_tes, 
            V_whb, V_phb

end