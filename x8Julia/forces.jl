using MAT

function forces(t, y, P, u, wind)
    P = matread("/Users/ethanclement/Documents/Internship DGA/x8Julia/x8_param.mat")
    # y : state vector
    # P : .mat file containing the parameters of the X8 Skywalker
    # u : control inputs elevator, aileron and throttle (no rudder)
    P["gravity"] = 9.81
    P["rho"] = 1.2250

    
    # pos = y[1:3]
    # Unpack state vector
    Theta = y[4:6]
    vel = y[7:9]
    rate = y[10:12]

    phi = Theta[1]
    theta = Theta[2]
    psi = Theta[3]
    p = rate[1] + wind[4]
    q = rate[2] + wind[5]
    r = rate[3] + wind[6]
    
    # Control inputs
    elevator = u[1]
    aileron = u[2]
    rudder = u[3]
    throttle = u[4]

    # Relative velocity
    wind_b = wind[1:3]
    vel_r = vel - wind_b
    u_r = vel_r[1]
    v_r = vel_r[2]
    w_r = vel_r[3]

    # Compute airspeed Va, angle-of-attack alpha, and side-slip beta
    Va = sqrt(u_r^2 + v_r^2 + w_r^2)
    Va = Va == 0 ? 1e-5 : Va  # Avoid division by zero

    alpha = atan(w_r, u_r) # Angle of attack
    beta = asin(v_r / Va)  # Side-slip angle

    # Compute gravitational force in body frame
    fg_N = [0.0, 0.0, P["mass"] * P["gravity"]]  # Gravity in NED frame
    fg_b = Rzyx(phi, theta, psi)' * fg_N   # Convert to body frame by multiplying NED frame to rotation matrix Rzyx

    # Longitudinal mode
    C_L_alpha = P["C_L_0"] + P["C_L_alpha"] * alpha # Lift coefficient 
    f_lift_s = 0.5 * P["rho"] * Va^2 * P["S_wing"] * (C_L_alpha + P["C_L_q"] * P["c"] / (2 * Va) * q + P["C_L_delta_e"] * elevator) # Lift equation

    C_D_alpha = P["C_D_0"] + P["C_D_alpha1"] * alpha + P["C_D_alpha2"] * alpha^2 # Drag coefficient (angle of attack)
    C_D_beta = P["C_D_beta1"] * beta + P["C_D_beta2"] * beta^2 # Drag coefficient (side-slip angle)
    f_drag_s = 0.5 * P["rho"] * Va^2 * P["S_wing"] * (C_D_alpha + C_D_beta + P["C_D_q"] * P["c"] / (2 * Va) * q + P["C_D_delta_e"] * elevator^2) # Drag equation

    # pitch moment
    m_a = P["C_m_0"] + P["C_m_alpha"] * alpha 
    m = 0.5 * P["rho"] * Va^2 * P["S_wing"] * P["c"] * (m_a + P["C_m_q"] * P["c"] / (2 * Va) * q + P["C_m_delta_e"] * elevator) # Pitching moment

    # Lateral mode
    f_y = 0.5 * P["rho"] * Va^2 * P["S_wing"] * (P["C_Y_0"] + P["C_Y_beta"] * beta + P["C_Y_p"] * P["b"] / (2 * Va) * p + P["C_Y_r"] * P["b"] / (2 * Va) * r + P["C_Y_delta_a"] * aileron + P["C_Y_delta_r"] * rudder) # Side force
    l = 0.5 * P["rho"] * Va^2 * P["b"] * P["S_wing"] * (P["C_l_0"] + P["C_l_beta"] * beta + P["C_l_p"] * P["b"] / (2 * Va) * p + P["C_l_r"] * P["b"] / (2 * Va) * r + P["C_l_delta_a"] * aileron + P["C_l_delta_r"] * rudder) # Roll moment
    n = 0.5 * P["rho"] * Va^2 * P["b"] * P["S_wing"] * (P["C_n_0"] + P["C_n_beta"] * beta + P["C_n_p"] * P["b"] / (2 * Va) * p + P["C_n_r"] * P["b"] / (2 * Va) * r + P["C_n_delta_a"] * aileron + P["C_n_delta_r"] * rudder) # Yaw moment

    # Aero forces and torques
    F_aero = Rzyx(0, alpha, beta)' * [-f_drag_s, f_y, -f_lift_s] # Aerodynamic forces in the body frame
    T_aero = [l, m, n] # Aerodynamic moments 

    # Propulsion force
    Vd = Va + throttle * (P["k_motor"] - Va) # Dynamic pressure velocity equation
    F_prop = [0.5 * P["rho"] * P["S_prop"] * P["C_prop"] * Vd * (Vd - Va), 0, 0] # Propulsion forces
    T_prop = [-P["k_T_P"] * (P["k_Omega"] * throttle)^2, 0, 0] # Propulsion moments 

    # Sum forces and torques
    Force = F_prop + fg_b + F_aero
    Torque = T_aero + T_prop

    return vcat(Force, Torque)
end
