clc; clear all;

numel = 2;

% WALL TEMPS [K] (arbitaryly choosen)
T_iw = [650 500];
T_ow = [600 450];

% ENGINE PROPERTIES (arbitaryly choosen)
% ===Design Choosen===
% inner radius [m]
r_iw = [1.9375/39.37 1.8/39.37];
% Wall thicknesses [m]
t_iw = [1e-3 0.6e-3];
t_ow = [2e-3 2.1e-3];
% Channel wall thicknesses
w = [1e-3 1.2e-3]; % rib width
d = 2e-3;    % channel width
% Number of channels
N = 40;

% ===Design Calculated===
% outer radius
r_ow = r_iw+t_iw+d+t_ow;
% channel width
h = (2*pi*r_iw/N)-w;
% ===Area calcs=== [m^2]
r_chamber = r_iw(1); % [m] chamber radius
A_chamber = pi*r_chamber^2;
A_throat = pi*(0.75/39.37)^2;
A_channels = N*h*d;
%A_wall = 2*pi*r_iw.*t_iw - A_channels;
A_wall = pi*(r_ow.^2 - r_iw.^2) - A_channels;
A_cc = A_chamber-A_throat;

% Chamber pressure (arbitaryly choosen)
P_chamber = 400*6895; % [psi->Pa] 

% Pressure difference [Pa] channel(higher) to chamber (30  psi)
deltaP = 206843;

% Enviroment
P0 = 101325;  % [Pa]


% Material properties [Inconel 718]
[E_iw,E_ow,E_avg,alpha_iw,alpha_ow,alpha_avg,nu_iw,nu_ow,nu_avg] = MatProperties_Inconel718(T_iw,T_ow);


% Stress Calcs
[stress_total,stress_T_ow, stress_T_iw, stress_T_s, stress_T_c,stress_P_b, stress_P_c, stress_P_s, stress_P_t, stress_P_hoop, stress_P_a] = Stress_v1_BL(T_iw,T_ow,E_iw,E_ow,E_avg,alpha_iw,alpha_ow,alpha_avg,nu_iw,nu_ow,nu_avg,t_iw,t_ow,h,w,d,r_ow,r_iw,deltaP,P_chamber,A_wall,A_cc,numel)
