% CR3BP Orbital Propagator
% Purdue Astrodynamics
% Jack McKinley

clc
clear
close all

%% Initializations

% Characteristic Properties
charM = 6.046804e24; % Mass of Earth-Luna system
charL = 3.844e5; % Distance between Earth and Luna
charT = 3.751903e5; % Time
charMu = .012150585; % Mu of the Earth-Luna system
charV = charL / charT; % Velocity of the moon

% Celestial Bodies
mass1 = "Earth";
mass2 = "Luna";

% Initial States
spctInit = transpose([1.1849919972792, 0, 0, 0, -.367025154279506, -.231116501800945]); % Halo Orbit
% spctInit = transpose([.846915127, 0, 0, 0, -.0782405216, 0]); % Lyapunov
% spctInit = transpose([.979865884449209, 0, .02486, 0, -1.8106932305339, 0]); % Marginally stable orbit
% spctInit = transpose([]); % Template

% Positions
m1Pos = (-charMu) * charL;
m2Pos = (1 - charMu) * charL;

%% Calculations

% ODE45 Components
t_final = 4 * pi;
tspan = linspace(0, t_final, 1e4);

tols = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[t, stateout] = ode45(@threeBody_prop, tspan, spctInit, tols);

% Parse Stateout Vector into Positions
x_pos = stateout(:,1);
y_pos = stateout(:,2);
z_pos = stateout(:,3);

% Parse Stateout Vector into Velocities
x_vel = stateout(:,4);
y_vel = stateout(:,5);
z_vel = stateout(:,6);

% Colinear Lagrange Points
L1 = colinear_Lagrange(charMu, 0.8);
L2 = colinear_Lagrange(charMu, 1);
L3 = colinear_Lagrange(charMu, -1);

L1_dim = L1 * charL;
L2_dim = L2 * charL;
L3_dim = L3 * charL;

% Equilateral Lagrange Points
L4_x = 1/2 - charMu;
L4_y = sqrt(3) / 2;
L5_x = 1/2 - charMu;
L5_y = -sqrt(3) / 2;

L4_xdim = L4_x * charL;
L4_ydim = L4_y * charL;
L5_xdim = L5_x * charL;
L5_ydim = L5_y * charL;

% Position Dimensionalization
x_posDim = x_pos * charL;
y_posDim = y_pos * charL;
z_posDim = z_pos * charL;

% Velocity Dimensionalization
x_velDim = x_vel * charL;
y_velDim = y_vel * charL;
z_velDim = z_vel * charL;

%% Plots

% Plotting Earth, Luna and Lagrange Points
fig = figure();
plot3(m1Pos, 0, 0, "Marker", ".", "MarkerSize", 20, "MarkerEdgeColor", "green")
hold on
plot3(m2Pos, 0, 0, "Marker", ".", "MarkerSize", 20, "MarkerEdgeColor", "#999999")
plot3(L1_dim, 0, 0, "Marker", "v", "MarkerSize", 5, "MarkerFaceColor", "b", "MarkerEdgeColor", "b")
plot3(L2_dim, 0, 0, "Marker", "^", "MarkerSize", 5, "MarkerFaceColor", "b", "MarkerEdgeColor", "b")
plot3(L3_dim, 0, 0, "Marker", "p", "MarkerSize", 5, "MarkerFaceColor", "b", "MarkerEdgeColor", "b")
plot3(L4_xdim, L4_ydim, 0, "Marker", "d", "MarkerSize", 5, "MarkerFaceColor", "b", "MarkerEdgeColor", "b")
plot3(L5_xdim, L5_ydim, 0, "Marker", "s", "MarkerSize", 5, "MarkerFaceColor", "b", "MarkerEdgeColor", "b")

% Plotting spacecraft motion
plot3(x_posDim, y_posDim, z_posDim, "Color", "black")
grid on
camlight;
lighting gouraud
axis equal
title("Earth-Luna CR3BP")
xlabel("X-Position [km]")
ylabel("Y-Position [km]")
zlabel("Z-Position [km]")
legend("Earth", "Luna", "L_1", "L_2", "L_3", "L_4", "L_5")

%% Outputs

fprintf("Lagrange Point L1 is located at %.4f times the radius between Earth and Luna \n", L1)
fprintf("Lagrange Point L2 is located at %.4f times the radius between Earth and Luna \n", L2)
fprintf("Lagrange Point L3 is located at %.4f times the radius between Earth and Luna \n", L3)

%% Functions

% Calculates the state derivative via ODE45
function state_dot = threeBody_prop(t, state)
    %First parameter will always be "t", second parameter will always be "state"
    charMu = .012150585; % Mu of the Earth-Luna system
    
    % Positions (from state)
    x = state(1);
    y = state(2);
    z = state(3);
    
    % Velocities (from state)
    x_dot = state(4);
    y_dot = state(5);
    z_dot = state(6);

    % Radii
    r1 = [x + charMu; y; z];
    r2 = [x - 1 + charMu; y; z];
    r1Mag = norm(r1);
    r2Mag = norm(r2);
    
    xVel_dot = (2 * y_dot) + x - ((1 - charMu) * (x + charMu) / r1Mag^3) - ...
        (charMu * (x - 1 + charMu) / r2Mag^3);
    yVel_dot = (-2 * x_dot) + y - ((1 - charMu) * (y) / r1Mag^3) - ...
        (charMu * (y) / r2Mag^3);
    zVel_dot = - ((1 - charMu) * (z) / r1Mag^3) - (charMu * (z) / r2Mag^3);
    
    state_dot = [x_dot; y_dot; z_dot; xVel_dot; yVel_dot; zVel_dot];
end

% Calculates the Colinear Lagrange Points
function lagrange = colinear_Lagrange(charMu, x_n)
    % Initializations
    max_iterations = 100;
    precision = 10;
    tolerance = 1e-10;

    % Builing Colinear Lagrange Point Formula (from EoMs)
    syms x
    EoM_x = x - ((1 - charMu) * (x + charMu)) / abs(power((x + charMu), 3)) - ...
        (charMu * (x - 1 + charMu)) / abs(power((x - 1 + charMu), 3)); 
    diff_EoM = diff(EoM_x);

    % Iteratively Calculating Lagrange Points (Newton-Raphson Method)
    for i = 1:max_iterations
        x_n1 = vpa(x_n - (subs(EoM_x, x, x_n) / subs(diff_EoM, x, x_n)), precision);

        if (tolerance > abs(x_n1 - x_n))
            break
        end
        
        x_n = x_n1;
    end
    
    % Checking for Max Iterations
    if i == max_iterations
        disp('Maximum number of iterations reached. Solution may not have converged.');
    end

    lagrange = x_n1;
end

