% -------------------------------------------------------------------
% Patrick Heng
% 14 Nov 2025 - 6 Dec 2025 (Transient Solver)
% Script to plot dynamic results of the PFR boundary value problem 
% using a non-linear finite volume scheme.
% Uses the CUT loop to solve the coupled transport problems:
%   C - Concentration equation
%   U - U momentum equation
%   T - Temperature equation
% which is iteratively solved until convergence
%
% Run this AFTER running 'PFR_FV_test_CUT_transient.m'
% -------------------------------------------------------------------

close all; clear variables; clc
% Load dynamic data
dyn_data = load('PFR_dynamic_data.mat');

% Load variables
ZZ = dyn_data.Z;
TT = dyn_data.T;
m = dyn_data.flue_flow;
CA = dyn_data.CA_copy;
T = dyn_data.T_copy;
Ta = dyn_data.Ta_copy;

% Colors to plot with
colororder([0,0,1; 1,0,0; 0,0.7,0; 1,0.7,0; 0.8,0,1]);

% ----- CA line plot -----
figure(1)
delta = 33;     % Density of lines
strt = 1;       % Start index for axial position

plot3(ZZ(:,strt:delta:end),TT(:,strt:delta:end),CA(strt:delta:end,:), ...
    linewidth=2)

% Pretty plot parameters
view(90,0)
ylabel('$t$ (s)',interpreter='latex'); 
zlabel('$C_A\ (\mathrm{mol/m^3})$',interpreter='latex');
grid on; box on
fontname('Serif'); fontsize(14,'points');

z = ZZ(1,strt:delta:end);
% Legend labels
labels = cell(length(z),1);     
for i = 1:length(z)
    labels{i} = sprintf('%.2f', z(i));
end
legend(labels,interpreter='latex', location='southeast', ...
    orientation='horizontal')


% ----- T line plot -----
figure(2)
% Colors to plot with
colororder([0,0,1; 1,0,0; 0,0.7,0; 1,0.7,0; 0.8,0,1]);

plot3(ZZ(:,strt:delta:end),TT(:,strt:delta:end),T(strt:delta:end,:), ...
    linewidth=2)

% Pretty plot parameters
view(90,0)
ylabel('$t$ (s)',interpreter='latex'); 
zlabel('$T\ (\mathrm{K})$',interpreter='latex');
grid on; box on
fontname('Serif'); fontsize(14,'points');

% Legend labels
labels = cell(length(z),1);
for i = 1:length(z)
    labels{i} = sprintf('%.2f', z(i));
end
legend(labels,interpreter='latex', location='southeast', ...
    orientation='horizontal')




% ----- Control action (mass flow) plot -----
figure(3)

% Colors to plot with
colororder([0,0,1; 1,0,0; 0,0.7,0; 1,0.7,0; 0.8,0,1]);
plot(TT(:,1),m,linewidth=2)

% Pretty plot parameters
xlabel('$t$ (s)',interpreter='latex'); 
ylabel('$\dot{m}\ (\mathrm{kg/s})$',interpreter='latex');
legend('Flue Gas Flow Rate',interpreter='latex', location='southeast', ...
    orientation='horizontal')
grid on; box on
fontname('Serif'); fontsize(14,'points');


% ----- CA surface plot -----
figure(4)
surf(ZZ,TT,CA',facealpha=1,linestyle='none')

% Pretty plot parameters
cb = colorbar;
view(0,90)
xlabel('$z$',interpreter='latex'); 
ylabel('$t$ (s)',interpreter='latex');
ylabel(cb,'$C_A\ (\mathrm{mol/m^3})$',interpreter='latex')
clim([0,17])

grid on; box on
fontname('Serif'); fontsize(14,'points');

% ----- T surface plot -----
figure(5)
surf(ZZ,TT,T',facealpha=1,linestyle='none')

% Pretty plot parameters
cb = colorbar;
view(0,90)
xlabel('$z$',interpreter='latex'); 
ylabel('$t$ (s)',interpreter='latex');
ylabel(cb,'$T$ (K)',interpreter='latex')
clim([940,1040])

grid on; box on
fontname('Serif'); fontsize(14,'points');
