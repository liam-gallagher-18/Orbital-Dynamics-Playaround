clc; clear; close all;

%% Nine-Neptune-Sun 3BP
% Constants
G = 4 * pi^2;             % Gravitational constant
M_sun = 1;                % Mass of the Sun
M_neptune = 5.15e-5;      % Neptune mass (in Solar masses)
M_9 = 3.0e-5;             % Planet Nine estimated mass (in Solar masses) (~10 Earth masses)

% Neptune
a_nept = 30.069;  
e_nept = 0.0086;
r0_nept = a_nept*(1 + e_nept);                         % Start at aphelion
v0_nept = sqrt(G*M_sun*(2/r0_nept - 1/a_nept));        % Orbital velocity

% Planet Nine
a_9 = 410;                                              % Estimated semi-major axis (~380 AU, ranges from 300 to 520)
e_9 = 0.22;                                             % Estimated eccentricity (ranges from 0.15 to 0.4)
r0_9 = a_9*(1 - e_9);                                   % Start at perihelion (closer to Neptune)
v0_9 = sqrt(G*M_sun*(2/r0_9 - 1/a_9));

% State vector: [xN, yN, vxN, vyN, x9, y9, vx9, vy9]
y0 = [r0_nept; 0; 0; v0_nept; r0_9; 0; 0; v0_9];

% Neptune orbit is ~165 years. Planet Nine orbit ~15,000 years.
tspan = [0 150000]; % Years - ADJUST AS NEEDED... over long periods of time you can really see the drift between 2BP and 3BP (Figure 2)

% Solve ODE
options = odeset('RelTol', 1e-8);
[t, sol] = ode45(@(t, y) n_body_system(t, y, G, M_sun, M_neptune, M_9), tspan, y0, options);

%% Neptune-Sun 2BP
mu = 4 * pi^2; % Gravitational parameter for the Sun

state_nept_2BP = [r0_nept; 0; 0; v0_nept];
[tN, solN] = ode45(@(t2, y2) gravity_ode(t2, y2, mu), tspan, state_nept_2BP, options);

%% Plots
figure('Color', 'w', 'Position', [100, 100, 900, 400]);

% Subplot 1: Solar System
subplot(1,2,1);
plot(0, 0, 'yo', 'MarkerSize', 10, 'MarkerFaceColor', 'y'); hold on; % Sun
plot(sol(:,1), sol(:,2), 'b', 'LineWidth', 1);       % Neptune
plot(sol(:,5), sol(:,6), 'Color', 'm', 'LineWidth', 2); % Planet Nine
title('Solar System View (AU)');
axis equal; grid on; legend('Sun', 'Neptune', 'Planet Nine');
% Subplot 2: Neptune Focus
subplot(1,2,2);
hold on
plot(0, 0, 'yo', 'MarkerSize', 10, 'MarkerFaceColor', 'y'); 
plot(sol(:,1), sol(:,2), 'b', 'LineWidth', 1.5); hold on;
plot(solN(:,1), solN(:,2), 'r:');
title('Focus on Neptune Orbit');
xlabel('AU'); ylabel('AU');
xlim([-29.9 -29.7]); ylim([-0.1 0.1]);
grid on;
hold off

% Calculate the distance between the Neptune's expected position without Nine and its expected position with Nine
dist_difference = sqrt((sol(:,1) - solN(:,1)).^2 + (sol(:,2) - solN(:,2)).^2);
t_plot = t(1:end-5); % Remove last few points for plotting
dist_plot = 149597871*(dist_difference(1:end-5)); % Remove last few points and convert to km

figure('Color', 'w');

plot(t_plot, dist_plot, 'r', 'LineWidth', 2);
title('Neptune Displacement due to Planet Nine');
xlabel('Years');
ylabel('Displacement (km)');
grid on;

%% 3BP ODE Fn
function dydt = n_body_system(~, y, G, Ms, Mn, M9)
    % Extract positions
    r_nept = y(1:2); % Neptune [x; y]
    r_9 = y(5:6); % Planet Nine [x; y]
    
    % Distances
    dist_sNept = norm(r_nept);       % Sun to Neptune
    dist_s9 = norm(r_9);             % Sun to Nine
    dist_n9 = norm(r_nept - r_9);    % Neptune to Nine
    
    % Accelerations (Gravity: -G*M*r / |r|^3)
    % Accel on Neptune = (Pull from Sun) + (Pull from Planet Nine)
    a_n = -G*Ms*r_nept/dist_sNept^3 + G*M9*(r_9 - r_nept)/dist_n9^3;
    
    % Accel on Nine = (Pull from Sun) + (Pull from Neptune)
    a_9 = -G*Ms*r_9/dist_s9^3 + G*Mn*(r_nept - r_9)/dist_n9^3;
    
    dydt = [y(3); y(4); a_n(1); a_n(2); y(7); y(8); a_9(1); a_9(2)];
end

%% 2BP ODE Fn
function dydt2 = gravity_ode(~, state_nept, mu)
    x_nept = state_nept(1);
    y_nept = state_nept(2);
    vx_nept = state_nept(3);
    vy_nept = state_nept(4);
    
    r_nept = sqrt(x_nept^2 + y_nept^2);
    
    % Derivatives: [vx; vy; ax; ay]
    dydt2 = [vx_nept; 
            vy_nept; 
            -mu * x_nept / r_nept^3; 
            -mu * y_nept / r_nept^3];
end