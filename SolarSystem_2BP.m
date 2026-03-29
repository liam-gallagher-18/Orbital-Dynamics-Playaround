clc; clear; close all;

%% Solar System Playaround %%
% Using separate 2BP propagation for each individual planet with the Sun
% (effects of the planets on each other not considered)

% Constants (Units: AU and Years)
mu = 4 * pi^2; % Gravitational parameter for the Sun
tspan = [0 165]; % Simulate for 1 year

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
% Mars
a_mars = 1.524;  
e_mars = 0.0934;
% Jupiter
a_jupiter = 5.203;  
e_jupiter = 0.0484;
% Saturn
a_saturn = 9.537;  
e_saturn = 0.0542;
% Uranus
a_uranus = 19.191;  
e_uranus = 0.0472;
% Neptune
a_neptune = 30.069;  
e_neptune = 0.0086;

% Initial Conditions, Use: r_p = a*(1-e), v_p = sqrt(mu / a*(1+e) / (1-e))
% Mercury Initial [x; y; vx; vy]
r0_merc = a_merc*(1 - e_merc);
v0_merc = sqrt(mu/a_merc*(1 + e_merc)/(1 - e_merc));
state0_merc = [r0_merc; 0; 0; v0_merc];

% Venus Initial [x; y; vx; vy]
r0_ven = a_ven * (1 - e_ven);
v0_ven = sqrt(mu/a_ven * (1 + e_ven)/(1 - e_ven));
state0_ven = [r0_ven; 0; 0; v0_ven];

% Earth Initial [x; y; vx; vy]
r0_earth = a_earth * (1 - e_earth);
v0_earth = sqrt(mu/a_earth * (1 + e_earth)/(1 - e_earth));
state0_earth = [r0_earth; 0; 0; v0_earth];

% Mars Initial [x; y; vx; vy]
r0_mars = a_mars * (1 - e_mars);
v0_mars = sqrt(mu/a_mars * (1 + e_mars)/(1 - e_mars));
state0_mars = [r0_mars; 0; 0; v0_mars];

% Jupiter Initial [x; y; vx; vy]
r0_jupiter = a_jupiter * (1 - e_jupiter);
v0_jupiter = sqrt(mu/a_jupiter * (1 + e_jupiter)/(1 - e_jupiter));
state0_jupiter = [r0_jupiter; 0; 0; v0_jupiter];

% Saturn Initial [x; y; vx; vy]
r0_saturn = a_saturn * (1 - e_saturn);
v0_saturn = sqrt(mu/a_saturn * (1 + e_saturn)/(1 - e_saturn));
state0_saturn = [r0_saturn; 0; 0; v0_saturn];

% Uranus Initial [x; y; vx; vy]
r0_uranus = a_uranus * (1 - e_uranus);
v0_uranus = sqrt(mu/a_uranus * (1 + e_uranus)/(1 - e_uranus));
state0_uranus = [r0_uranus; 0; 0; v0_uranus];

% Neptune Initial [x; y; vx; vy]
r0_neptune = a_neptune * (1 - e_neptune);
v0_neptune = sqrt(mu/a_neptune * (1 + e_neptune)/(1 - e_neptune));
state0_neptune = [r0_neptune; 0; 0; v0_neptune];


% Numerical Integration
options = odeset('RelTol', 1e-8);
[tM, solM] = ode45(@(t, y) gravity_ode(t, y, mu), tspan, state0_merc, options);
[tV, solV] = ode45(@(t, y) gravity_ode(t, y, mu), tspan, state0_ven, options);
[tE, solE] = ode45(@(t, y) gravity_ode(t, y, mu), tspan, state0_earth, options);
[tMa, solMa] = ode45(@(t, y) gravity_ode(t, y, mu), tspan, state0_mars, options);
[tJ, solJ] = ode45(@(t, y) gravity_ode(t, y, mu), tspan, state0_jupiter, options);
[tS, solS] = ode45(@(t, y) gravity_ode(t, y, mu), tspan, state0_saturn, options);
[tU, solU] = ode45(@(t, y) gravity_ode(t, y, mu), tspan, state0_uranus, options);
[tN, solN] = ode45(@(t, y) gravity_ode(t, y, mu), tspan, state0_neptune, options);

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
plot(solM(:,1), solM(:,2), 'w--');
plot(solV(:,1), solV(:,2), 'y--');
plot(solE(:,1), solE(:,2), 'g--');
plot(solMa(:,1), solMa(:,2), 'r--');
plot(solJ(:,1), solJ(:,2), '--', 'Color', '#FF8C00', 'LineWidth', 1);
plot(solS(:,1), solS(:,2), 'm--');
plot(solU(:,1), solU(:,2), 'c--');
plot(solN(:,1), solN(:,2), 'b--');

% Moving handles for animation
hM = plot(solM(1,1), solM(1,2), 'wo', 'MarkerFaceColor', 'w', 'MarkerSize', 2); % Mercury
hV = plot(solV(1,1), solV(1,2), 'yo', 'MarkerFaceColor', 'y', 'MarkerSize', 4); % Venus
hE = plot(solE(1,1), solE(1,2), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 4); % Earth
hMa = plot(solMa(1,1), solMa(1,2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 3); % Mars
hJ = plot(solJ(1,1), solJ(1,2), 'o', 'MarkerFaceColor', '#FF8C00', 'MarkerSize', 8); % Jupiter
hS = plot(solS(1,1), solS(1,2), 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 8); % Saturn
hU = plot(solU(1,1), solU(1,2), 'co', 'MarkerFaceColor', 'c', 'MarkerSize', 6); % Uranus
hN = plot(solN(1,1), solN(1,2), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6); % Neptune
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

    posMa = interp1(tMa, solMa(:,1:2), tM(i));
    set(hMa, 'XData', posMa(1), 'YData', posMa(2));

    posJ = interp1(tJ, solJ(:,1:2), tM(i));
    set(hJ, 'XData', posJ(1), 'YData', posJ(2));

    posS = interp1(tS, solS(:,1:2), tM(i));
    set(hS, 'XData', posS(1), 'YData', posS(2));

    posU = interp1(tU, solU(:,1:2), tM(i));
    set(hU, 'XData', posU(1), 'YData', posU(2));

    posN = interp1(tN, solN(:,1:2), tM(i));
    set(hN, 'XData', posN(1), 'YData', posN(2));
    
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