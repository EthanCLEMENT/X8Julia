function dynamics(t, y, P, tau)
    #= Function to compute the dynamics of the X8-Skywalker =#
    
    # Unpack state vector
    pos = y[1:3]
    Theta = y[4:6]
    vel = y[7:9]
    Omega = y[10:12]

    # Calculate Coriolis and centripetal matrix C_rb
    C_rb = [zeros(3, 3)              -P["mass"] * Smtrx(vel) - P["mass"] * Smtrx(Omega) * Smtrx(P["r_cg"]);
            -P["mass"] * Smtrx(vel) + P["mass"] * Smtrx(P["r_cg"]) * Smtrx(Omega)  -Smtrx(P["I_cg"] * Omega)]

    # Calculate generalized acceleration ny_dot
    ny_dot = P["M_rb"] \ (tau - C_rb * vcat(vel, Omega))
    vel_dot = ny_dot[1:3]
    Omega_dot = ny_dot[4:6]

    # Calculate linear and angular velocity
    pos_dot = Rzyx(Theta[1], Theta[2], Theta[3]) * vel
    Theta_dot = TransformationMatrix(Theta) * Omega

    # Return the state derivative
    xdot = vcat(pos_dot, Theta_dot, vel_dot, Omega_dot)
    return xdot
end
