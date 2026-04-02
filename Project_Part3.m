clear;
clc;
close all;

%% Number of Orbits Around Mars
N = 5;

%% Mars Parameters
R_mars = 3396000;                  % Radius of Mars [m]
G_mars = 6.674e-11;
M_mars = 6.41693e23;
mu_mars = G_mars*M_mars;           % m^3/s^2
sol = 88775;                       % [s]
omega_mars = 2*pi / sol;           % rad/s

%% Combined Maneuver Calculation
% Givens
a_0 = 33000000;             % [m]
r_0 = a_0;
a_f = 250000 + R_mars;      % [m]
r_f = a_f;
a_trans = (r_0 + r_f)/2;    % [m]

v_0         = sqrt(mu_mars/r_0);                            % [m/s]
v_trans_0   = sqrt((2*mu_mars/r_0) - (mu_mars/a_trans));
v_f         = sqrt(mu_mars/r_f);
v_trans_f   = sqrt((2*mu_mars/r_f) - (mu_mars/a_trans));

Delta_i     = 90;       % [deg]

% Define Delta_v function
Delta_v_tot = @(s) ...
    sqrt(v_0^2 + v_trans_0^2 - 2*v_0*v_trans_0*cosd(s*Delta_i)) + ...
    sqrt(v_f^2 + v_trans_f^2 - 2*v_f*v_trans_f*cosd((1-s)*Delta_i));

% Minimize Delta_v_tot
[s_opt, dv_min] = fminbnd(Delta_v_tot, 0, 1);

fprintf('Optimal s = %.6f\n', s_opt);
fprintf('Minimum total Delta v = %.6f m/s\n', dv_min);

%% Plot Transfer
% Initial Orbit
i_0 = 0;

% Transfer Orbit (Elliptical)
i_trans = s_opt * Delta_i; 
e_trans = (r_0 - r_f) / (r_0 + r_f); 
p_trans = a_trans * (1 - e_trans^2);

% Final Orbit
i_f = 90;

% Plot
figure('Color','w'); 
hold on; grid on; axis equal;
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('3D Combined Inclination + Altitude Change Maneuver');
view(3);

% Draw Mars (Scaled to km)
[Xm,Ym,Zm] = sphere(60);
surf(R_mars*Xm/1000, R_mars*Ym/1000, R_mars*Zm/1000, ...
    'FaceColor','r','EdgeColor','none','FaceAlpha',0.6, 'DisplayName', 'Mars');

% Trajectories
% Initial
theta = linspace(-pi, pi, 200);
x_1 = r_0 * cos(theta) / 1000;
y_1 = r_0 * sin(theta) / 1000;
z_1 = zeros(size(theta));
plot3(x_1, y_1, z_1, 'b', 'LineWidth', 1.5, 'DisplayName', 'Initial Orbit');

% Transfer
nu_trans = linspace(pi, 2*pi, 100); 
r_vals = p_trans ./ (1 + e_trans * cos(nu_trans));

x_pq = r_vals .* cos(nu_trans);
y_pq = r_vals .* sin(nu_trans);
z_pq = zeros(size(nu_trans));

points_trans = zeros(3, length(nu_trans));
Rot_Trans = [1, 0, 0; 
             0, cosd(i_trans), -sind(i_trans); 
             0, sind(i_trans),  cosd(i_trans)];

for k1 = 1:length(nu_trans)
    pt = [x_pq(k1); y_pq(k1); z_pq(k1)];
    points_trans(:,k1) = Rot_Trans * pt;
end

plot3(points_trans(1,:)/1000, points_trans(2,:)/1000, points_trans(3,:)/1000, ...
      'm', 'LineWidth', 2, 'DisplayName', 'Transfer Orbit');

% Final
theta = linspace(0, 2*pi, 200);
x_2_peri = r_f * cos(theta);
y_2_peri = r_f * sin(theta);
z_2_peri = zeros(size(theta));

points_final = zeros(3, length(theta));
Rot_Final = [1, 0, 0; 
             0, cosd(i_f), -sind(i_f); 
             0, sind(i_f),  cosd(i_f)];

for k1 = 1:length(theta)
    pt = [x_2_peri(k1); y_2_peri(k1); z_2_peri(k1)];
    points_final(:,k1) = Rot_Final * pt;
end

plot3(points_final(1,:)/1000, points_final(2,:)/1000, points_final(3,:)/1000, ...
      'k', 'LineWidth', 1.5, 'DisplayName', 'Final Orbit');

legend show;
hold off;

%% Propagate Orbit Around Mars
% Initial time
t0 = 0;
theta0 = 0;

% Orbital period
n = sqrt(mu_mars/a_f^3);          % rad/s
T = 2*pi*sqrt(a_f^3/mu_mars);
tf = N*T;                         % Total seconds (changes based on # or orbit at top of code)

% Time vector FOR ORBITS AROUND MARS, NOT TRANSFER
t = linspace(t0, tf, N*500);       
tspan = [t0 tf];
nu_0 = 0;

% Loop over time
lat = zeros(1,length(t));
lon = zeros(1,length(t));

for k1 = 1:length(t)
    dt = t(k1); % Time since start [s]
    
    % Propagate Mars' rotation
    theta = theta0 + omega_mars*dt;         % radians

    % Propagate true anomaly
    nu = nu_0 + rad2deg(n*dt);              % degrees

    % COEs to ECI
    r = a_f;                                  
    rPeri = [ r*cosd(nu);
              r*sind(nu);
              0 ];

    C_NP = [1, 0, 0;
            0, cosd(i_f), -sind(i_f);
            0, sind(i_f),  cosd(i_f)];

    rECI = C_NP*rPeri;

    % ECI to ECEF
    C_EE = [ cos(theta),  sin(theta), 0;
            -sin(theta),  cos(theta), 0;
             0,           0,          1 ];

    rECEF = C_EE*rECI;

    % ECEF to LLA
    x = rECEF(1);
    y = rECEF(2);
    z = rECEF(3);
    lat(k1) = asind(z/norm(rECEF));
    lon(k1) = rad2deg(atan2(y, x));
end

% Clean up horizontal wrap lines for plot
lon_diff = diff(lon);
jump_indices = find(abs(lon_diff) > 300); 
lon(jump_indices) = NaN;
lat(jump_indices) = NaN;

%% When Will Orbit Repeat?
% Satellite rotation time is given by T
% Mars rotation time is one sol
orbits_per_sol = sol / T;

% Does this ratio ever equal a whole number?

repeats = false;
for j = 1:1:N
    future_orbits = orbits_per_sol*j;
    if abs(future_orbits - round(future_orbits)) < 10^(-9)
        fprintf('Orbit repeats at orbit #%.0f!', N)
        repeats = true;
        break
    end
end

if ~repeats
    fprintf('Orbit does not repeat within %.0f orbits.', N);
    disp(' ')
end

% How wide is the LineWidth?
ax_gt = gca; 
originalUnits = get(ax_gt, 'Units');
set(ax_gt, 'Units', 'points');
pos = get(ax_gt, 'Position');
axesWidthInPoints = pos(3);
percentage = (1.5 / axesWidthInPoints) * 100;
set(ax_gt, 'Units', originalUnits);
fprintf('A LineWidth of 1.5 covers %.4f%% of the graph width.\n', percentage);
%% Plot Orbit Path
% Set up 2BP propagation
r0 = [a_f; 0; 0];
% Initial velocity for 90 deg inclination
v_circ = sqrt(mu_mars/a_f);
v0_peri = [0; v_circ; 0]; % Velocity in perifocal frame

% Perifocal to ECI
C_NP = [1, 0, 0;
        0, cosd(i_f), -sind(i_f);
        0, sind(i_f),  cosd(i_f)];

r_prop = C_NP * r0;
v_prop = C_NP * v0_peri;
state0 = [r_prop; v_prop];

% Solve 2BP
opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t_ode, R] = ode45(@(t,y) twoBodyProp(t, y, mu_mars), tspan, state0, opts);

% ECI to ECEF
r_final = R(:,1:3);
rECEF_3D = zeros(size(r_final));

for k1 = 1:length(t_ode)
    theta = theta0 + omega_mars * t_ode(k1);
    C_EE = [ cos(theta)  sin(theta)  0;
            -sin(theta)  cos(theta)  0;
             0           0           1 ];
    rECEF_3D(k1,:) = (C_EE * r_final(k1,:)')';
end

% Plot
figure('Color','w');
plot3(rECEF_3D(:,1)/1000, rECEF_3D(:,2)/1000, rECEF_3D(:,3)/1000, 'g', 'LineWidth', 1.0); 
hold on
surf(R_mars*Xm/1000, R_mars*Ym/1000, R_mars*Zm/1000, ...
    'FaceColor','r','EdgeColor','none','FaceAlpha',0.8);
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
grid on; axis equal;
title('Mars Orbit (Fixed Frame)');
hold off;

%% Plot Ground Track
figure('Color','w');
ax = gca;
xlim(ax, [-180 180]);
ylim(ax, [-90 90]);
hold(ax, 'on');

plot(lon,lat,'linewidth',1.5)

I = imread('mars_map_labeled.jpg');
h = imagesc([-180 180], [90 -90], I);
uistack(h,'bottom')
set(ax, 'YDir', 'normal'); 
axis(ax, 'normal');        
set(ax, 'Layer', 'top');
grid on;
box on;
xticks(-180:60:180);
yticks(-90:30:90);
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
title(['Ground Track (' num2str(N) ' Orbits)'], 'Interpreter', 'none', 'FontSize', 12);
hold off

%% ANIMATION 1: The Transfer Maneuver
% Create a new figure for the animation
fig1 = figure('Color','k', 'Name', 'Transfer Orbit', 'Position',[100 100 1260 946], 'Renderer','opengl');
axis equal; grid on; box on;
hold on;
view(3);
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w');
xlim([-40000 40000]); ylim([-40000 40000]); zlim([-40000 40000]);

% Draw Mars
surf(R_mars*Xm/1000, R_mars*Ym/1000, R_mars*Zm/1000, ...
'FaceColor','r', 'EdgeColor','none', 'FaceAlpha', 0.9);
lighting gouraud; light('Position', [1 1 0], 'Style', 'infinite');

% Trajectories
% Initial Orbit
P1 = [x_1; y_1; z_1];
% Transfer Orbit
P2 = points_trans/1000;
% Final Orbit
P3 = points_final/1000;
% Combine into one continuous path for the animation loop
Full_Path = [P1, P2, P3];

% Create graphics objects
hTrail = animatedline('Color', 'c', 'LineWidth', 1.5);
hSat   = plot3(0,0,0, 'o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);

title('Animation: Transfer Maneuver', 'Color', 'w');
vid1 = VideoWriter('C:\Users\13143\OneDrive - University of Dayton\Documents\Grad School\Orbital Mechanics\MarsTransfer.mp4', 'MPEG-4');
% vid1.FrameRate = 30;
open(vid1);

% Animation Loop
for k1 = 1:4:length(Full_Path)
    frame = getframe(fig1);
    writeVideo(vid1, frame);

    x = Full_Path(1,k1);
    y = Full_Path(2,k1);
    z = Full_Path(3,k1);
    
    % Update Satellite Position
    set(hSat, 'XData', x, 'YData', y, 'ZData', z);
    
    % Add to Trail
    addpoints(hTrail, x, y, z);
    
    % Dynamic Title
    if k1 < length(P1)
        txt = 'Phase 1: Initial Orbit';
    elseif k1 < length(P1) + length(P2)
        txt = 'Phase 2: Transfer Orbit';
    else
        txt = 'Phase 3: Final Orbit';
    end
    title(txt, 'Color', 'w', 'FontSize', 14);

    drawnow limitrate;
end

close(vid1);

%% ANIMATION 2: Long-Term Orbit & Ground Drift
fig2 = figure('Color','k', 'Name', 'Orbit Propagation', 'Position',[100 100 1260 946], 'Renderer','opengl');
axis equal; grid on; box on;
hold on;
view(3);
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w');

% Draw Mars with Map
try
    I_tex = flipud(I); 
    surf(R_mars*Xm/1000, R_mars*Ym/1000, R_mars*Zm/1000, ...
        I_tex, 'FaceColor', 'texturemap', 'EdgeColor', 'none', 'FaceAlpha', 1);
catch
    surf(R_mars*Xm/1000, R_mars*Ym/1000, R_mars*Zm/1000, ...
    'FaceColor',[0.8 0.4 0.2], 'EdgeColor','none', 'FaceAlpha', 0.8);
end

% Set Axis Limits
limit_val = (R_mars + 1000000)/1000;
xlim([-limit_val limit_val]);
ylim([-limit_val limit_val]);
zlim([-limit_val limit_val]);

% Prepare Data
rECEF_animation = rECEF_3D / 1000; % Convert to km

% Create Graphics Objects
hTrail2 = animatedline('Color', 'g', 'LineWidth', 1);
hSat2   = plot3(rECEF_animation(1,1), rECEF_animation(1,2), rECEF_animation(1,3), ...
    'o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);

title('Animation: 3D Orbit Path', 'Color', 'w');

% Video Setup
vid2 = VideoWriter('C:\Users\13143\OneDrive - University of Dayton\Documents\Grad School\Orbital Mechanics\MarsOrbit.mp4', 'MPEG-4');
% vid2.FrameRate = 30;
open(vid2);

target_frames = 144000;
step_size = ceil(length(rECEF_animation) / target_frames);

% Animation Loop
for k2 = 1:step_size:length(rECEF_animation)
    frame = getframe(fig2);
    writeVideo(vid2, frame);
    
    x = rECEF_animation(k2,1);
    y = rECEF_animation(k2,2);
    z = rECEF_animation(k2,3);
    
    set(hSat2, 'XData', x, 'YData', y, 'ZData', z);
    addpoints(hTrail2, x, y, z);

    drawnow; 
end

close(vid2);

%% 2BP Function
function xdot = twoBodyProp(~, r, mu)
    r_fn = r(1:3);
    v_fn = r(4:6);
    r_scalar_fn = norm(r_fn);
    a_fn = -mu * r_fn / r_scalar_fn^3;
    xdot = [v_fn; a_fn];
end