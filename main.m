clc; clear; close all;

%% Inputs

% CEA Inputs
Pc = 413; % psia
Pe = 13.7; % psia
O_F = 1.1;
T_inlet = 300; % K

k_w = 20; %W/m K
res = 0; % thermal resistance coating
Thrust = 1600; % lbf
%Pe = 13; % psi 
C_star_eff = 0.94;
C_F_eff = 0.9;
L_star = 35; % in
h_ch = 0.001/0.0254; % in
w_ch_min = 0.001/0.0254; % in
%num_ch = 50;
w_rib = 0.001/0.0254; % in
D_c = 3.875; % in
t_ins = 0.001/0.0254; % in
t_out = 0.001/0.0254; % in

% CEA Outputs
MW_g = [21.834   21.836   21.929   22.040];
gamma = [1.1590   1.1591   1.1698   1.2082];
mu_g = [0.96585  0.96454  0.91604  0.65672]*1e-4;
Cp_g = [3.3244   3.3184   2.9049   2.1909]*1000;
k_g = [5.6321   5.6144   4.5272   2.0893]*1e-1;
Pr_g = [0.5701   0.5701   0.5878   0.6887];
C_star = 1627.5;
C_F = 1.4862;
T_tot = 2863.55;

%% Solution

num_nodes = 50;
angle_conv = 40;
[A,D,M,P,T,w_ch,h_ch,D_h,w_rib,num_ch,t_ins,t_out,step,pos,D_t,A_t,Pc,Pe,mdot,MW_g,gamma,mu_g,Cp_g,k_g,Pr_g,id_th,id_c] = geometry(Thrust,Pc,Pe,C_star,C_star_eff,C_F,C_F_eff,L_star,angle_conv,h_ch,w_rib,w_ch_min,MW_g,gamma,mu_g,Cp_g,k_g,Pr_g,t_ins,t_out,T_tot,num_nodes);
mdot_f = mdot*1/(1+O_F); %kg/s

ratio = 0.75;

% D_t = D_t*0.0254; % m
% A_t = A_t*0.0254^2; % m^2

[T_c,T_sat_c,P_c,q,T_chg,T_rhg,T_ci,T_ri,T_cb,T_rb,T_rt,T_ro,T_co,T_ct] = heatTransfer2D(Pc,M,P,A,D,D_h,w_ch,h_ch,w_rib,num_ch,t_ins,t_out,step,pos,gamma,mu_g,Cp_g,k_g,Pr_g,C_star,T,k_w,D_t,A_t,T_inlet,ratio,mdot_f,C_star_eff,res);


% ===STRESS===
% Temperature Dependent Material properties
% Is averaging these temps the right approach???
% Inner Wall
T_iw = (T_rhg+T_chg)/2;
T_iw = double(T_iw);
[row_iw,E_iw,nu_iw,alpha_iw,k_iw,Cp_iw,Yield_iw] = MatProperties_Steel174ph(T_iw);
% Outer Wall
T_ow = (T_ro+T_co)/2;
T_ow = double(T_ow);
[row_ow,E_ow,nu_ow,alpha_ow,k_ow,Cp_ow,Yield_ow] = MatProperties_Steel174ph(T_ow);

% TEMPORARY VARIABLE NAME CONVERSION (Will update code)
t_iw = t_ins;
t_ow = t_out;
N = num_ch;
h = h_ch;
w = w_rib;
d = w_ch_min;
r_iw = D_c/2;
r_ow = r_iw+t_iw+d+t_ow;
deltaP = 206843;  % Pressure difference [Pa] channel(higher) to chamber (30  psi)
P_chamber = Pc*6895; % [psi to Pa]
A_channels = N*h*d;
A_wall = pi*(r_ow.^2 - r_iw.^2) - A_channels;
A_chamber = pi*r_iw.^2;
A_throat = pi*(D_t/2).^2;
A_cc = A_chamber-A_throat;

% Conley's function call
%[stress_total,stress_T_ow, stress_T_iw, stress_T_s, stress_T_c,stress_P_b, stress_P_c, stress_P_s, stress_P_t, stress_P_hoop, stress_P_a,yield] = stress(T_ow,T_iw,E_iw,E_ow,alpha_iw,alpha_ow,nu_iw,nu_ow,P,P_c,A,D,w_ch,w_rib,t_ins,h_ch,t_out,num_ch);
% Brandon's function call
[stressTotal,stressTow, stressTiw, stressTs, stressTc, stressPb, stressPs, stressPt, stressP_hoopow, stressP_hoopiw, stressPa,vonMises] = stress(T_iw,T_ow,E_iw,E_ow,alpha_iw,alpha_ow,nu_iw,nu_ow,t_iw,t_ow,h,w,d,r_ow,r_iw,deltaP,P_chamber,A_wall,A_cc)

% Total Inner wall stress only
stressTotaliw = stressTiw + stressP_hoopiw;
% ===========


P_c = P_c/6894.76; % psi

Q = q(1:end-1).*D(1:end-1).*step*pi;
% Loss = sum(Q) % W

figure (1)
plot(pos,D/2)
title('contour')
axis equal

figure(2)
clf
hold on
plot(pos,T_c,'LineWidth',1,'Color','b')
plot(pos,T_sat_c,'--','LineWidth',1,'Color','r')
xline(0, '--', 'Exit','HandleVisibility','off')
xline(pos(id_th), '--', 'Throat','HandleVisibility','off')
xline(pos(id_c), '--', 'Chamber','HandleVisibility','off')
hold off
title('Coolant Temperature')
xlabel('Axial Distance (m)')
ylabel('Temperature (K)')
legend('T_c','T_{saturation}')
grid on

figure(3)
clf
hold on
plot(pos,P_c,'LineWidth',1)
xline(0, '--', 'Exit','HandleVisibility','off')
xline(pos(id_th), '--', 'Throat','HandleVisibility','off')
xline(pos(id_c), '--', 'Chamber','HandleVisibility','off')
hold off
title('Coolant Pressure')
xlabel('Axial Distance (m)')
ylabel('Pressure (psia)')
grid on

figure(4)
clf
hold on
plot(pos,q,'LineWidth',1)
xline(0, '--', 'Exit','HandleVisibility','off')
xline(pos(id_th), '--', 'Throat','HandleVisibility','off')
xline(pos(id_c), '--', 'Chamber','HandleVisibility','off')
hold off
title('Heat Flux')
xlabel('Axial Distance (m)')
ylabel('Heat Flux (W/m^2)')
grid on


figure(5)
clf
hold on
xline(pos(id_th), '--', 'Throat','HandleVisibility','off')
xline(pos(id_c), '--', 'Chamber','HandleVisibility','off')
plot(pos,T_chg,'LineWidth',1)
plot(pos,T_rhg,'LineWidth',1)
xline(0, '--', 'Exit','HandleVisibility','off')
hold off
title('Hot Wall Temperature')
xlabel('Axial Distance (m)')
ylabel('Temperature (K)')
legend('Coolant Section','Rib Section')
grid on

figure(6)
clf
hold on
xline(pos(id_th), '--', 'Throat','HandleVisibility','off')
xline(pos(id_c), '--', 'Chamber','HandleVisibility','off')
plot(pos,stressTotaliw/6895000,'LineWidth',1)
plot(pos,Yield_iw/6895000,'LineWidth',1)
xline(0, '--', 'Exit','HandleVisibility','off')
hold off
title('Stress')
xlabel('Axial Distance (m)')
ylabel('Stress (ksi)')
legend('Inner Wall Stress','Inner Wall Yield')
grid on