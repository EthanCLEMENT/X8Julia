using LinearAlgebra
using DifferentialEquations
using Plots
using DelimitedFiles
using MAT
include("Smtrx.jl")
include("Rxyz.jl")
include("Smtrx.jl")
include("dynamics.jl")
include("forces.jl")
include("TransformationMatrix.jl")

mutable struct ControllerState
    i_phi::Float64
    i_h::Float64
    i_V::Float64
    i_theta::Float64
    i_chi::Float64
end

function simX8()
    # Load parameters from 'x8_param.mat'
    P = matread("x8Julia/x8_param.mat")

    # Inertia matrix: assuming symmetry with respect to xz-plane -> Jxy=Jyz=0
    P["I_cg"] = [P["Jx"] 0 -P["Jxz"];
                 0 P["Jy"] 0;
                 -P["Jxz"] 0 P["Jz"]]

    # Mass matrix
    P["M_rb"] = [I(3) * P["mass"] -P["mass"] * Smtrx(P["r_cg"]);
                 P["mass"] * Smtrx(P["r_cg"]) P["I_cg"]]

    P["rho"] = 1.2250
    P["gravity"] = 9.81

    tend = 400.0
    tspan = (0.0, tend)

    y0 = [0.0, 0.0, -200.0,   # Position
          0.0, 0.0, 0.0,      # Euler angles
          18.0, 0.0, 0.0,     # Velocity
          0.0, 0.0, 0.0]      # Rates

    wind = zeros(6)  # Wind velocity components in body frame

    do_trim = false

    if do_trim
        y_trim, u_trim = findTrim(y0, P)
    else
        u_trim = [0.0370, 0.0000, 0.0, 0.1219]
        y_trim = [0.0000, -0.0000, -200.0000, 0.0000, 0.0308, 0.0000,
                  17.9914, 0.0000, 0.5551, 0.0000, 0.0000, 0.0000]
    end

    # Set controller gains
    P["kp_h"]     = -0.025
    P["ki_h"]     = 0.0

    P["kp_theta"] = 0.1
    P["kd_theta"] = -0.01
    P["ki_theta"] = 0.0

    P["kp_V"]     = -0.05
    P["ki_V"]     = -0.01

    P["kp_phi"]   = -0.5
    P["ki_phi"]   = 0.0
    P["kd_phi"]   = 0.0

    P["kp_chi"]   = -0.05
    P["ki_chi"]   = 0.0

    # Initialize controller state
    ctrl_state = ControllerState(0.0, 0.0, 0.0, 0.0, 0.0)

    # Define the ODE function
    function ODEfunc!(dy, y, p, t)
        # Compute the control input
        ref_vals = ref(t, y)
        u = controller(t, y, P, u_trim, ref_vals, ctrl_state)

        # Compute the forces and moments
        F = forces(t, y, P, u, wind)

        # Compute the dynamics
        dy[:] = dynamics(t, y, P, F)
    end

    # Set up the problem
    prob = ODEProblem(ODEfunc!, y_trim, tspan)

    # Solve the ODE
    sol = solve(prob, Tsit5(), saveat=0.1)

    # Extract time and state variables
    t = sol.t
    y = sol.u

    # Plotting
    D = [yi[3] for yi in y]  # Altitude (D)
    phi = [rad2deg(yi[4]) for yi in y]
    theta = [rad2deg(yi[5]) for yi in y]
    psi = [rad2deg(yi[6]) for yi in y]
    u_vel = [yi[7] for yi in y]
    v_vel = [yi[8] for yi in y]
    w_vel = [yi[9] for yi in y]
    p_rate = [rad2deg(yi[10]) for yi in y]
    q_rate = [rad2deg(yi[11]) for yi in y]
    r_rate = [rad2deg(yi[12]) for yi in y]

    # Create plots
    plot(t, D, label="D", xlabel="Time (s)", ylabel="Altitude (m)", title="Altitude vs Time")
    plot(t, phi, label="φ", xlabel="Time (s)", ylabel="Angle (deg)", title="Euler Angles vs Time")
    plot!(t, theta, label="θ")
    plot!(t, psi, label="ψ")

    # Plotting u_vel, v_vel, w_vel
    plot(t, u_vel, label="u_vel", xlabel="Time (s)", ylabel="Velocity (m/s)", title="Velocities vs Time")
    plot!(t, v_vel, label="v_vel")
    plot!(t, w_vel, label="w_vel")

    # Plotting p_rate, q_rate, r_rate
    plot(t, p_rate, label="p_rate", xlabel="Time (s)", ylabel="Rate (deg/s)", title="Rates vs Time")
    plot!(t, q_rate, label="q_rate")
    plot!(t, r_rate, label="r_rate")


end

function ref(t, y)
    V = 18.0
    if t < 20
        chi = 0.0
    elseif t < 35
        chi = (t - 20) * 1 * pi / 180
    else
        chi = 15 * pi / 180
    end

    climb_rate = 0.2
    if t < 75
        h = 200.0
    elseif t < 325
        h = 200.0 + (t - 75) * climb_rate
    else
        h = 250.0
    end

    return [V; chi; h]
end

function controller(t, y, P, u_trim, ref, ctrl_state)
    pos = y[1:3]
    Theta = y[4:6]
    vel = y[7:9]
    Omega = y[10:12]

    V_ref = ref[1]
    chi_ref = ref[2]
    h_ref = ref[3]

    # Compute navigation frame velocity
    v_n = Rzyx(Theta[1], Theta[2], Theta[3]) * vel
    chi = atan(v_n[2], v_n[1])

    # Update integrator states
    ctrl_state.i_chi += chi - chi_ref
    phi_ref = P["kp_chi"] * (chi - chi_ref) + P["ki_chi"] * ctrl_state.i_chi

    ctrl_state.i_h += -pos[3] - h_ref
    theta_ref = P["kp_h"] * (-pos[3] - h_ref) + P["ki_h"] * ctrl_state.i_h

    ctrl_state.i_theta += Theta[2] - theta_ref
    delta_e = P["kp_theta"] * (Theta[2] - theta_ref) + P["ki_theta"] * ctrl_state.i_theta - P["kd_theta"] * Omega[2]

    ctrl_state.i_phi += Theta[1] - phi_ref
    delta_a = P["kp_phi"] * (Theta[1] - phi_ref) + P["ki_phi"] * ctrl_state.i_phi - P["kd_phi"] * Omega[1]

    delta_r = 0.0  # Rudder input (assuming no sideslip control)

    ctrl_state.i_V += norm(vel) - V_ref
    delta_t = P["kp_V"] * (norm(vel) - V_ref) + P["ki_V"] * ctrl_state.i_V

    u = [delta_e; delta_a; delta_r; delta_t] .+ u_trim

    # Saturate control inputs
    u = clamp.(u, -1.0, 1.0)
    u[4] = max(0.0, u[4])  # Throttle saturation between 0 and 1

    return u
end

simX8()
