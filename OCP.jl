using Plots
using JuMP
using Ipopt


##* Reading Model Parameters from Module
include("parameters.jl")
      using Main.Bounds     
            x_guess = Bounds.x0;    z_guess = Bounds.z0;    u_guess = Bounds.u0;    #des_guess = Bounds.des0
            Nx = Bounds.Nx;         Nz = Bounds.Nz;         Nu = Bounds.Nu;         #Ndes = Bounds.Ndes 
            ls_x = Bounds.ls_x;     ls_z = Bounds.ls_z;     ls_u = Bounds.ls_u;     #ls_des = Bounds.ls_des 
            us_x = Bounds.us_x;     us_z = Bounds.us_z;     us_u = Bounds.us_u;     #us_des = Bounds.us_des

      using Main.Model_parameters
            ρ_dh = Model_parameters.ρ_dh;             Cp_dh = Model_parameters.Cp_dh;                 q_dh = Model_parameters.q_dh 
            ρ_wh = Model_parameters.ρ_wh;             Cp_wh = Model_parameters.Cp_wh;                 q_wh = Model_parameters.q_wh
            T_dh_ret= Model_parameters.T_dh_ret;      T_dh_minSup = Model_parameters.T_dh_minSup     
            V_tes = Model_parameters.V_tes           
            V_whb = Model_parameters.V_whb;           V_phb = Model_parameters.V_phb     
      
function Collocation_Matrix()
      #Radau
      t1 = 0.155051
      t2 = 0.644949
      t3 = 1.0

      M1 = [
            t1 1 / 2 * t1^2 1 / 3 * t1^3
            t2 1 / 2 * t2^2 1 / 3 * t2^3
            t3 1 / 2 * t3^2 1 / 3 * t3^3
            ]
      M2 = [
            1 t1 t1^2
            1 t2 t2^2
            1 t3 t3^2
            ]

      M = M1 * inv(M2)
      return M
end

##* Function to solve OCP
# function Build_OCP(Q_whb, Tf, (ns,np) )

#region-> Value of Arguments for Debugging
      # Q_whb = vcat(1.0*ones(10,1), ones(10,1), 1.0*ones(10,1)) *1.2539999996092727e6
      Q_whb = vcat(1.0*ones(10,1), 1.2*ones(10,1), 0.8*ones(10,1)) *1.254e6
      Tf = 30.0
      (ns, np) = (1,1)
#endregion

                  #region-> Setting Initial guesses and Dimensions
                  # dx0_us  = 0*x_guess
                  # alg0_us = 0 * z0
                  # Nx = size(x0, 1)
                  # Nz = size(z0, 1)
                  # Nu = size(u0, 1)
                  #endregion

      ##? Set OCP Parameters here
      Solve_OCP         = true
      Display_Plots     = true
      T0 = 0.0
      dt = 1.0
      NFE = convert(Int32, (Tf - T0)/dt)
      # NFE = 30
      # dt = (Tf - T0) / NFE
      NCP = 3

##* Defining Solver
      model1 = Model(with_optimizer(Ipopt.Optimizer, print_level=5 , max_iter = 1500))

      #region -> Variables and Objective
            ## Declare Variables
            # @variable(model1, des[1:Ndes])                  #scaled volume
            @variable(model1, x0[1:Nx])                     #initial guess - scaled

            @variable(model1, x[1:Nx, 1:NFE, 1:NCP])
            @variable(model1, dx_us[1:Nx, 1:NFE, 1:NCP])

            @variable(model1, z[1:Nz, 1:NFE, 1:NCP])
            # @variable(model1, alg_us[1:Nz, 1:NFE, 1:NCP])

            @variable(model1, u[1:Nu, 1:NFE])

                  #region-> #? Set Variable Bounds AND Initial Guesses (scaled)
                  for nx in 1:Nx, nz in 1:Nz, nu in 1:Nu, nfe in 1:NFE, ncp in 1:NCP     
                        # set_lower_bound(des[ndes], 0)
                        # set_upper_bound(des[ndes], 1)

                        # set_lower_bound(x0[nx], 0)
                        # set_upper_bound(x0[nx], 1)

                        set_lower_bound(x[nx, nfe, ncp], 0)
                        set_upper_bound(x[nx, nfe, ncp], 1)

                        set_lower_bound(z[nz, nfe, ncp], 0)
                        #set_upper_bound(z[nz, nfe, ncp], 999)

                        #set_lower_bound(dx_us[nx, nfe, ncp], 0)
                        #set_upper_bound(dx_us[nx, nfe, ncp], 999)

                        #set_lower_bound(alg_us[nz, nfe, ncp], 0)
                        #set_upper_bound(alg_us[nz, nfe, ncp], 999)

                        set_lower_bound(u[nu, nfe], 0)
                        set_upper_bound(u[nu, nfe], 1)
                  end

                        #Initial Guesses (Scaled)
                  for nx in 1:Nx, nz in 1:Nz, nu in 1:Nu, nfe in 1:NFE, ncp in 1:NCP
                        # set_start_value(des[ndes],                des_guess[ndes])
                        # set_start_value(x0[nx],                   x_guess[nx])
                        set_start_value(x[nx, nfe, ncp],          x_guess[nx])
                        set_start_value(z[nz, nfe, ncp],          z_guess[nz])
                        set_start_value(dx_us[nx, nfe, ncp],      0.0)
                        # set_start_value(alg_us[nz, nfe, ncp],     0)
                        set_start_value(u[nu, nfe],               u_guess[nu])

                  end
                  #endregion

                  #region-> Expressions for Unscaling Variables (makes it easier to write DAE Equation) #?add variables and change indices 
                  @NLexpressions(model1, begin                             
                        T_tes[nfe in 1:NFE, ncp in 1:NCP],        x[1, nfe, ncp]            *  (us_x[1]     - ls_x[1])   + ls_x[1]
                        T_phb[nfe in 1:NFE, ncp in 1:NCP],        x[2, nfe, ncp]            *  (us_x[2]     - ls_x[2])   + ls_x[2]
                        T_whb[nfe in 1:NFE, ncp in 1:NCP],        x[3, nfe, ncp]            *  (us_x[3]     - ls_x[3])   + ls_x[3]


                        q_bypass[nfe in 1:NFE, ncp in 1:NCP],    z[1, nfe, ncp]            *  (us_z[1]     - ls_z[1])   + ls_z[1]
                        T_a[nfe in 1:NFE, ncp in 1:NCP],         z[2, nfe, ncp]            *  (us_z[2]     - ls_z[2])   + ls_z[2]
                        T_b[nfe in 1:NFE, ncp in 1:NCP],         z[3, nfe, ncp]            *  (us_z[3]     - ls_z[3])   + ls_z[3]
                        T_c[nfe in 1:NFE, ncp in 1:NCP],         z[4, nfe, ncp]            *  (us_z[4]     - ls_z[4])   + ls_z[4]

                        q_whb[nfe in 1:NFE],                     u[1, nfe]                 *  (us_u[1]     - ls_u[1])   + ls_u[1] 
                        q_a[nfe in 1:NFE],                       u[2, nfe]                 *  (us_u[2]     - ls_u[2])   + ls_u[2] 
                        q_b[nfe in 1:NFE],                       u[3, nfe]                 *  (us_u[3]     - ls_u[3])   + ls_u[3] 
                        Q_phb[nfe in 1:NFE],                     u[4, nfe]                 *  (us_u[4]     - ls_u[4])   + ls_u[4] 

                  end)
                  #endregion
                  
                  # fix(z_ADMM, 60.0)
                  if np == 1
                        for nx in 1:Nx
                              fix(x0[nx], x_guess[nx], force = true)
                        end
                  end

                  #!For degugging
                  for nfe in 1:NFE, ncp in 1:NCP
                        # fix(u[1,nfe], 0.2, force=true)      #q_whb at
                        # fix(u[2,nfe], 0.0, force=true)      #q_a
                        # fix(u[3,nfe], 0.0, force=true)      #q_b
                        fix(z[1,nfe,ncp], 0.0, force = true)


                  end
                  for nfe in 1:Int(NFE/3)
                        # fix(u[3,nfe], 0.0, force=true)      #q_b
                  end


            # Objective
            @NLobjective(model1, Min, sum( u[4,nfe]^2       + 1e0*u[2,nfe]*u[3,nfe]    + (1e-2)*u[2,nfe]^2  + (1e-2)*u[3,nfe]^2    for nfe in 1:NFE) )  #Duty in KJ

      #endregion

##* Defining Constraints

      @NLconstraints(model1, begin
            #?Defining the model ODEs in each line
            Constr_ODE1[nfe in 1:NFE, ncp in 1:NCP], dx_us[1, nfe, ncp]      ==  (q_a[nfe]*(T_a[nfe]        - T_tes[nfe, ncp]) + q_b[nfe]*(T_b[nfe]   - T_tes[nfe, ncp])  ) /V_tes            #DegC/hr
            Constr_ODE2[nfe in 1:NFE, ncp in 1:NCP], dx_us[2, nfe, ncp]      ==  (q_dh *(T_c[nfe, ncp]      - T_phb[nfe, ncp]) )/V_phb                + Q_phb[nfe]/( V_phb*ρ_dh*Cp_dh)
            Constr_ODE3[nfe in 1:NFE, ncp in 1:NCP], dx_us[2, nfe, ncp]      ==  (q_whb[nfe]*(T_a[nfe, ncp] - T_whb[nfe, ncp])) /V_whb                + Q_whb[nfe]/( V_whb*ρ_dh*Cp_dh)
            #In case of more states - pattern
            #Constr_ODE999[nfe=1:NFE, ncp=1:NCP], dx[999,nfe,ncp] ==
      end)

      @NLconstraints(model1, begin
            #?Defining Model Algebraic Equations in each line
            Constr_Alg1[nfe in 1:NFE, ncp in 1:NCP], q_bypass[nfe, ncp]       ==   q_dh - q_a[nfe] + q_b[nfe] - q_whb[nfe] 
            Constr_Alg2[nfe in 1:NFE, ncp in 1:NCP], T_a[nfe, ncp]            ==  (q_dh*T_dh_ret + q_b[nfe]*T_tes[nfe,ncp])/(q_dh + q_b[nfe]) 
            Constr_Alg3[nfe in 1:NFE, ncp in 1:NCP], T_b[nfe, ncp]            ==  (q_whb[nfe]*T_whb[nfe,ncp] + q_a[nfe]*T_tes[nfe,ncp])/(q_whb[nfe] + q_a[nfe]) 
            Constr_Alg4[nfe in 1:NFE, ncp in 1:NCP], T_c[nfe, ncp]            ==  (q_bypass[nfe,ncp]*T_a[nfe,ncp] + (q_dh - q_bypass[nfe,ncp])*T_b[nfe,ncp])/q_dh 
            #In case of more states - pattern
            #Constr_Alg999[nfe=1:NFE, ncp=1:NCP], alg[999,nfe,ncp] ==
      end)

      @NLconstraints(model1, begin
            #?Defining any Inequality Constraints in each line
            [nfe in 1:NFE, ncp in 1:NCP], T_dh_minSup <= T_phb[nfe, ncp] <= 70.0
            
            #In case of more states - pattern
            #Constr_Ineq999[nfe=1:NFE, ncp=1:NCP], alg[999,nfe,ncp] ==
      end)

                  #region-> std code -> Collocation Equations
                        collMat = Collocation_Matrix()
                        @NLconstraints(model1, begin
                              #Collocation Equation for Differential Equations (scaled form)
                              #t = 0
                              Constr_Coll_Diff0[nx in 1:Nx,  nfe = 1, ncp in 1:NCP],        x[nx, nfe, ncp]    == x0[nx]               + dt * sum(collMat[ncp, i] * dx_us[nx, nfe, i] for i = 1:NCP) /(us_x[nx] - ls_x[nx])
                              #t = 1 ... (N-1)
                              Constr_Coll_Diff[nx in 1:Nx,   nfe in 2:NFE, ncp = 1:NCP],    x[nx, nfe, ncp]    == x[nx, nfe-1, NCP]    + dt * sum(collMat[ncp, i] * dx_us[nx, nfe, i] for i = 1:NCP) /(us_x[nx] - ls_x[nx])
                        end)
                        
                  #endregion-> std code     

##* Solve the model
      # if Solve_OCP == true
            optimize!(model1)
            JuMP.termination_status(model1)
            JuMP.solve_time(model1::Model)
            star_Obj = JuMP.objective_value(model1)
            
            star_dx_us = JuMP.value.(dx_us[:,:,NCP])
            star_x0 = JuMP.value.(x0)
            star_u = JuMP.value.(u)

            #Getting Varables directly from expressions - Unscaled Units
            star_x0_us = star_x0            .*  (us_x - ls_x) + ls_x

            star_T_tes = JuMP.value.(T_tes[:, NCP])
            star_T_phb  = JuMP.value.(T_phb[:, NCP])
            star_T_whb  = JuMP.value.(T_whb[:, NCP])
                              star_T_tes = cat(star_x0_us[1], star_T_tes, dims = 1)  
                              star_T_phb = cat(star_x0_us[2], star_T_phb, dims = 1) 
                              star_T_whb = cat(star_x0_us[3], star_T_whb, dims = 1) 


            star_q_bypass = JuMP.value.(q_bypass[:,NCP])
            star_T_a      = JuMP.value.(T_a[:,NCP])
            star_T_b      = JuMP.value.(T_b[:, NCP])
            star_T_c      = JuMP.value.(T_c[:, NCP])

            
            star_q_whb  = JuMP.value.(q_whb[:])
            star_q_a    = JuMP.value.(q_a[:])  
            star_q_b    = JuMP.value.(q_b[:]) 
            star_Q_phb  = JuMP.value.(Q_phb[:])
      # end
      
##* Plot Solution
      # if Display_Plots == true
                  star_cost = sum(star_Q_phb[nfe] for nfe in 1:NFE)
                  t_plot = collect(T0:dt:Tf)    #Returns NFE+1 dimensional vector
                  
                  #choose backend for plots
                  plotly()
                  # gr()

                        #Differential States
                        p11 = plot(t_plot, star_T_tes,                  label = "T_tes")
                        p11 = plot!(t_plot, star_T_phb,                 label = "T_phb")
                        p11 = plot!(t_plot, star_T_whb,                 label = "T_whb",                                title = "Differential States, @Obj $star_cost")

                        #Manipulated Variables
                        p12 = plot(t_plot[1:end-1], star_q_whb,         label = "q_whb",        linetype = :steppost)
                        p12 = plot!(t_plot[1:end-1], star_q_a,          label = "q_a",          linetype = :steppost)
                        p12 = plot!(t_plot[1:end-1], star_q_b,          label = "q_b",          linetype = :steppost,   title = "Inputs")
                        
                        p13 = plot(t_plot[1:end-1], star_Q_phb,         label ="Q_phb (kJ/hr)", linetype = :steppost,   title = "Inputs")

                        #Algebraic States
                        p14 = plot(t_plot[2:end], star_q_bypass,        label = "q_bypass",                             title = "Alg states")
                        
                        p15 = plot(t_plot[2:end], star_T_a,             label = "T_a")  
                        p15 = plot!(t_plot[2:end], star_T_b,            label = "T_b")  
                        p15 = plot!(t_plot[2:end], star_T_c,            label = "T_c",                                  title = "Alg states") 


      # end
##* Display Plots
p11
p12
p13
p14
p15

println("Total cost is $star_cost")
# return model1
# end



