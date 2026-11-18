%% AERSP 309 Project 2
% Compare the orbit with and without the J2 effect.
% All functions used by this script are included at the bottom.
%
% Use km for distance, km/s for velocity, and seconds for time.
% Calculations use radians; printed results and plots use degrees.
% Saving the figures requires MATLAB R2020a or later.

clear; clc; close all;

%% Constants and initial conditions
mu = 3.986e5;          % km^3/s^2
RE = 6371;            % km
J2 = 1.081874e-3;      % J2 coefficient given in the project

r0 = [8207; 7114; 5253];
v0 = [-3.0347; 2.7144; 3.2901];
state0 = [r0; v0];

% Save the plots in a folder called figures next to this script.
scriptFolder = fileparts(mfilename('fullpath'));
figureFolder = fullfile(scriptFolder, 'figures');

if ~exist(figureFolder, 'dir')
    mkdir(figureFolder);
end

%% Part 1: Find the initial orbital elements
% The function returns [a, e, i, RAAN, omega, nu].
coe0 = rv2coe(r0, v0, mu);

a0 = coe0(1);
e0 = coe0(2);
i0 = coe0(3);

T = 2*pi*sqrt(a0^3/mu);

fprintf('INITIAL ORBITAL ELEMENTS\n');
fprintf('a     = %.10f km\n', a0);
fprintf('e     = %.12f\n', e0);
fprintf('i     = %.10f deg\n', rad2deg(i0));
fprintf('RAAN  = %.10f deg\n', rad2deg(coe0(4)));
fprintf('omega = %.10f deg\n', rad2deg(coe0(5)));
fprintf('nu    = %.10f deg\n', rad2deg(coe0(6)));
fprintf('Period = %.10f s = %.10f h\n\n', T, T/3600);

%% Part 2: Simulate one orbit with point-mass gravity
options = odeset('RelTol', 1e-8, ...
                 'AbsTol', 1e-8, ...
                 'MaxStep', T/200);

[tPoint, statePoint] = ode45( ...
    @(t,y) equationsOfMotion(t,y,mu,RE,0), ...
    linspace(0,T,2001), state0, options);

coePoint = elementHistory(statePoint, mu);
plotPoint = elementsForPlot(coePoint);

%% Part 3: Add J2 and simulate five orbital periods
% Use the initial two-body period T to set the simulation length.
% With J2, the orbit does not close exactly after T.

[tJ2, stateJ2] = ode45( ...
    @(t,y) equationsOfMotion(t,y,mu,RE,J2), ...
    linspace(0,5*T,10001), state0, options);

coeJ2 = elementHistory(stateJ2, mu);
plotJ2 = elementsForPlot(coeJ2);

% Keep a separate copy of the results for the first period.
tJ2One = tJ2(1:2001);
stateJ2One = stateJ2(1:2001,:);
coeJ2One = coeJ2(1:2001,:);

%% Optional: Compare the RAAN and perigee precession rates
% p = a*(1-e^2), so (RE/p)^2 = (RE/a)^2/(1-e^2)^2.
% 5*cos(i)^2 - 1 = 4 - 5*sin(i)^2.
%
% Estimate the average precession rates using the initial orbital elements.
% These estimates leave out the short-period variations.

p0 = a0*(1-e0^2);
n0 = sqrt(mu/a0^3);

RAANrate = -1.5*n0*J2*(RE/p0)^2*cos(i0);
omegaRate = 0.75*n0*J2*(RE/p0)^2*(5*cos(i0)^2-1);

RAANpredicted = rad2deg(coe0(4) + RAANrate*tJ2);
omegaPredicted = rad2deg(coe0(5) + omegaRate*tJ2);

%% Check the results and print a summary
pointDeviation = max( ...
    abs(plotPoint(:,1:5) - plotPoint(1,1:5)), [], 1);

fprintf('POINT-MASS MAXIMUM DEVIATIONS FROM INITIAL VALUES\n');
fprintf('a:     %.6g km\n', pointDeviation(1));
fprintf('e:     %.6g\n', pointDeviation(2));
fprintf('i:     %.6g deg\n', pointDeviation(3));
fprintf('RAAN:  %.6g deg\n', pointDeviation(4));
fprintf('omega: %.6g deg\n', pointDeviation(5));

fprintf('True-anomaly advance: %.10f deg\n', ...
    plotPoint(end,6)-plotPoint(1,6));

closureError = norm(statePoint(end,1:3)-r0.');
fprintf('One-orbit position closure error: %.6g km\n\n', ...
    closureError);

% Check that energy and angular momentum stay constant without J2.
rPoint = vecnorm(statePoint(:,1:3),2,2);
energyPoint = sum(statePoint(:,4:6).^2,2)/2 - mu./rPoint;
hPoint = cross(statePoint(:,1:3),statePoint(:,4:6),2);

pointEnergyError = max(abs( ...
    (energyPoint-energyPoint(1))/energyPoint(1)));

pointHError = max(vecnorm(hPoint-hPoint(1,:),2,2)) ...
    /norm(hPoint(1,:));

% With J2, check the full energy and the z-component of angular momentum.
rJ2 = vecnorm(stateJ2(:,1:3),2,2);
zJ2 = stateJ2(:,3);

potentialJ2 = -mu./rJ2 ...
    + mu*J2*RE^2./(2*rJ2.^3).*(3*(zJ2./rJ2).^2-1);

energyJ2 = sum(stateJ2(:,4:6).^2,2)/2 + potentialJ2;
hJ2 = cross(stateJ2(:,1:3),stateJ2(:,4:6),2);

J2EnergyError = max(abs( ...
    (energyJ2-energyJ2(1))/energyJ2(1)));

J2HzError = max(abs(hJ2(:,3)-hJ2(1,3))) ...
    /abs(hJ2(1,3));

assert(pointEnergyError < 1e-7, ...
    'Point-mass energy conservation check failed.');
assert(pointHError < 1e-7, ...
    'Point-mass angular momentum check failed.');
assert(J2EnergyError < 1e-7, ...
    'J2 energy conservation check failed.');
assert(J2HzError < 1e-7, ...
    'J2 z-angular-momentum check failed.');

fprintf('CONSERVATION CHECKS\n');
fprintf('Point-mass relative energy error: %.6g\n', ...
    pointEnergyError);
fprintf('Point-mass relative angular momentum error: %.6g\n', ...
    pointHError);
fprintf('J2 relative total energy error: %.6g\n', ...
    J2EnergyError);
fprintf('J2 relative h_z error: %.6g\n\n', J2HzError);

J2range = max(plotJ2,[],1)-min(plotJ2,[],1);
J2change = plotJ2(end,:)-plotJ2(1,:);

elementNames = {'a (km)'; 'e'; 'i (deg)'; 'RAAN (deg)'; ...
                'omega (deg)'; 'nu (deg, unwrapped)'};

J2summary = table(elementNames,J2range.',J2change.', ...
    'VariableNames', {'Element','PeakToPeak','FinalMinusInitial'});

fprintf('J2 RESULTS OVER FIVE NOMINAL PERIODS\n');
disp(J2summary);

fprintf('Predicted RAAN rate: %.10f deg/day\n', ...
    rad2deg(RAANrate)*86400);
fprintf('Predicted omega rate: %.10f deg/day\n', ...
    rad2deg(omegaRate)*86400);

fprintf('Predicted five-period RAAN change: %.10f deg\n', ...
    rad2deg(RAANrate)*5*T);
fprintf('Predicted five-period omega change: %.10f deg\n\n', ...
    rad2deg(omegaRate)*5*T);

%% Figure 1: Point-mass trajectory
fig1 = orbitPlot(statePoint,RE, ...
    'Point-mass Earth: one orbit');

saveFigure(fig1,figureFolder,'01_point_mass_orbit');

%% Figure 2: Point-mass orbital elements
fig2 = figure('Color','w','Visible','off', ...
    'Position',[100,100,1100,850]);

axPoint = elementPlot(tPoint,plotPoint, ...
    'Point-mass Earth: orbital elements over one orbit');

% Choose axis limits that keep small numerical errors in perspective.
% The actual deviations are printed in the Command Window.
for k = 1:5
    center = plotPoint(1,k);
    halfWidth = max(abs(center)*1e-3,1e-4);
    ylim(axPoint(k),center+[-halfWidth,halfWidth]);
end

saveFigure(fig2,figureFolder,'02_point_mass_elements');

%% Figure 3: Point-mass true anomaly
fig3 = figure('Color','w','Visible','off', ...
    'Position',[100,100,1000,600]);

plot(tPoint/3600,plotPoint(:,6),'LineWidth',1.6);
grid on;
xlabel('Time (h)');
ylabel('Unwrapped true anomaly (deg)');
title('Point-mass Earth: true anomaly over one orbit');

saveFigure(fig3,figureFolder,'03_point_mass_true_anomaly');

%% Figure 4: J2 trajectory
fig4 = orbitPlot(stateJ2,RE, ...
    'Point-mass + J2 Earth: five nominal orbits');

saveFigure(fig4,figureFolder,'04_J2_five_orbits');

%% Figure 5: J2 orbital elements
fig5 = figure('Color','w','Visible','off', ...
    'Position',[100,100,1100,850]);

axJ2 = elementPlot(tJ2,plotJ2, ...
    'Point-mass + J2 Earth: five nominal orbits');

hold(axJ2(4),'on');
plot(axJ2(4),tJ2/3600,RAANpredicted,'--','LineWidth',1.4);
legend(axJ2(4),'Numerical osculating','First-order secular', ...
    'Location','best');

hold(axJ2(5),'on');
plot(axJ2(5),tJ2/3600,omegaPredicted,'--','LineWidth',1.4);
legend(axJ2(5),'Numerical osculating','First-order secular', ...
    'Location','best');

saveFigure(fig5,figureFolder,'05_J2_elements');

% Show all five figures once they have been saved.
set([fig1,fig2,fig3,fig4,fig5],'Visible','on');

fprintf('All conservation checks passed.\n');
fprintf('Five PNG figures saved in:\n%s\n',figureFolder);

%% Functions used by the script

function coe = rv2coe(r,v,mu)
% Convert ECI position and velocity to [a, e, i, RAAN, omega, nu].
% Return angles in radians and use NaN for any undefined angles.

    r = r(:);
    v = v(:);

    R = norm(r);
    V = norm(v);

    if R == 0
        error('Position must be nonzero.');
    end

    h = cross(r,v);
    H = norm(h);

    if H <= 1e-12*R*max(V,1)
        error('Classical elements are undefined for radial motion.');
    end

    node = cross([0;0;1],h);
    N = norm(node);

    eccentricityVector = cross(v,h)/mu-r/R;
    e = norm(eccentricityVector);

    energy = V^2/2-mu/R;

    if abs(energy) <= 10*eps(mu/R)
        a = Inf;
    else
        a = -mu/(2*energy);
    end

    inclination = atan2(hypot(h(1),h(2)),h(3));

    RAAN = NaN;
    omega = NaN;
    nu = NaN;

    if N/H > 1e-12
        RAAN = mod(atan2(node(2),node(1)),2*pi);

        if e > 1e-12
            sinOmega = dot(cross(node,eccentricityVector),h) ...
                /(N*e*H);
            cosOmega = dot(node,eccentricityVector)/(N*e);

            omega = mod(atan2(sinOmega,cosOmega),2*pi);
        end
    end

    if e > 1e-12
        sinNu = dot(cross(eccentricityVector,r),h)/(e*R*H);
        cosNu = dot(eccentricityVector,r)/(e*R);

        nu = mod(atan2(sinNu,cosNu),2*pi);
    end

    coe = [a,e,inclination,RAAN,omega,nu];
end

function dydt = equationsOfMotion(~,y,mu,RE,J2)
% Calculate the position and velocity derivatives in the ECI frame.
% Set J2 to zero to use only point-mass gravity.

    r = y(1:3);
    v = y(4:6);
    R = norm(r);

    acceleration = -mu*r/R^3;

    if J2 ~= 0
        q = (r(3)/R)^2;
        factor = 1.5*J2*mu*RE^2/R^5;

        J2acceleration = factor*[ ...
            r(1)*(5*q-1);
            r(2)*(5*q-1);
            r(3)*(5*q-3)];

        acceleration = acceleration+J2acceleration;
    end

    dydt = [v;acceleration];
end

function coe = elementHistory(states,mu)
% Find all six orbital elements at each output time.

    numberOfStates = size(states,1);
    coe = zeros(numberOfStates,6);

    for k = 1:numberOfStates
        coe(k,:) = rv2coe(states(k,1:3),states(k,4:6),mu);
    end
end

function values = elementsForPlot(coe)
% Convert angles to degrees for plotting.
% Unwrap RAAN, omega, and nu so they do not jump at 360 degrees.

    values = coe;
    values(:,3) = rad2deg(coe(:,3));
    values(:,4:6) = rad2deg(unwrap(coe(:,4:6),[],1));
end

function fig = orbitPlot(states,RE,heading)
% Plot the spacecraft trajectory in ECI coordinates.
% The Earth sphere shows the reference radius; J2 is handled in the equations.

    fig = figure('Color','w','Visible','off', ...
        'Position',[100,100,1000,800]);

    [sx,sy,sz] = sphere(60);

    surf(RE*sx,RE*sy,RE*sz, ...
        'FaceColor',[0.4,0.65,0.85], ...
        'EdgeColor','none', ...
        'FaceAlpha',0.65, ...
        'DisplayName','Earth (radius reference)');

    hold on;

    plot3(states(:,1),states(:,2),states(:,3), ...
        'Color',[0.75,0.15,0.12], ...
        'LineWidth',1.6, ...
        'DisplayName','Trajectory');

    plot3(states(1,1),states(1,2),states(1,3),'ko', ...
        'MarkerFaceColor','k', ...
        'DisplayName','Initial state');

    axis equal;
    grid on;
    view(40,25);

    xlabel('ECI x (km)');
    ylabel('ECI y (km)');
    zlabel('ECI z (km)');
    title(heading);
    legend('Location','best');
end

function ax = elementPlot(t,values,heading)
% Plot the six orbital elements against time.

    labels = {'a (km)', ...
              'e', ...
              'i (deg)', ...
              'RAAN (deg)', ...
              'Argument of perigee (deg)', ...
              'True anomaly (deg, unwrapped)'};

    tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
    ax = gobjects(6,1);

    for k = 1:6
        ax(k) = nexttile;
        plot(t/3600,values(:,k),'LineWidth',1.3);
        grid on;
        xlabel('Time (h)');
        ylabel(labels{k});
        ax(k).YAxis.Exponent = 0;
    end

    sgtitle(heading);
end

function saveFigure(fig,folder,name)
% Save each figure as a PNG using the filename passed to this function.

    exportgraphics(fig,fullfile(folder,[name,'.png']), ...
        'Resolution',160);
end