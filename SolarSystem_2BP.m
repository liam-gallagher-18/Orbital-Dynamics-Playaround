clc; clear; close all;

%% Solar System Playaround %%
% Using separate 2BP propagation for each individual planet with the Sun
% (effects of the planets on each other not considered)

% Constants (Units: AU and Years)
mu = 4 * pi^2; % Gravitational parameter for the Sun
tspan = [0 1]; % Simulate for 1 year

% Planetary Data (Semi-major axis and Eccentricity)
% Mercury
a_merc = 0.387; 
e_merc = 0.2056;
% Venus
a_ven = 0.723;  
e_ven = 0.0067;
% Earth
a_earth = 1;  
e_earth = 0.017;

% Initial Conditions, Use: r_p = a(1-e), v_p = sqrt( mu/a * (1+e)/(1-e) )

% Mercury Initial [x; y; vx; vy]
r0_merc = a_merc * (1 - e_merc);
v0_merc = sqrt(mu/a_merc * (1 + e_merc)/(1 - e_merc));
state0_merc = [r0_merc; 0; 0; v0_merc];

% Venus Initial [x; y; vx; vy]
r0_ven = a_ven * (1 - e_ven);
v0_ven = sqrt(mu/a_ven * (1 + e_ven)/(1 - e_ven));
state0_ven = [r0_ven; 0; 0; v0_ven];

% Earth Initial [x; y; vx; vy]
r0_earth = a_earth * (1 - e_earth);
v0_earth = sqrt(mu/a_earth * (1 + e_earth)/(1 - e_earth));
state0_earth = [r0_earth; 0; 0; v0_earth];

% Numerical Integration
options = odeset('RelTol', 1e-8);
[tM, solM] = ode45(@(t, y) gravity_ode(t, y, mu), tspan, state0_merc, options);
[tV, solV] = ode45(@(t, y) gravity_ode(t, y, mu), tspan, state0_ven, options);
[tE, solE] = ode45(@(t, y) gravity_ode(t, y, mu), tspan, state0_earth, options);

% Plotting and Animation
figure('Color', 'k'); 
hold on; 
axis equal; 
grid on;
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
title('Orbits of Planets', 'Color', 'w');
xlabel('AU'); 
ylabel('AU');

% Draw the Sun
plot(0, 0, 'yo', 'MarkerSize', 10, 'MarkerFaceColor', 'y'); 

% Static Orbits
plot(solM(:,1), solM(:,2), 'r--');
plot(solV(:,1), solV(:,2), 'w--');
plot(solE(:,1), solE(:,2), 'g--');

% Moving handles for animation
hM = plot(solM(1,1), solM(1,2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6); % Mercury
hV = plot(solV(1,1), solV(1,2), 'co', 'MarkerFaceColor', 'c', 'MarkerSize', 8); % Venus
hE = plot(solE(1,1), solE(1,2), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8); % Earth
% legend('Sun', 'Mercury Orbit', 'Venus Orbit', 'Mercury', 'Venus', 'TextColor', 'w');

% Animation Loop
for i = 1:length(tM)
    % Update Mercury position
    set(hM, 'XData', solM(i,1), 'YData', solM(i,2));
    
    % For other planets, find the position at the same time index
    % (Interpolate because ode45 uses variable steps)
    posV = interp1(tV, solV(:,1:2), tM(i));
    set(hV, 'XData', posV(1), 'YData', posV(2));

    posE = interp1(tE, solE(:,1:2), tM(i));
    set(hE, 'XData', posE(1), 'YData', posE(2));
    
    drawnow;
    pause(0.02); % Slow down animation to see it
end

%% ODE Fn
function dydt = gravity_ode(~, state, mu)
    x = state(1);
    y = state(2);
    vx = state(3);
    vy = state(4);
    
    r = sqrt(x^2 + y^2);
    
    % Derivatives: [vx; vy; ax; ay]
    dydt = [vx; 
            vy; 
            -mu * x / r^3; 
            -mu * y / r^3];
end