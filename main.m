clc; clear; close all;

%% Inputs

% CEA Inputs
Pc = 413; % psia
Pe = 13; % psia
O_F = 1.1;
T_inlet = 300; % K

k_w = 20; %W/m K
res = 0; % thermal resistance coating
Thrust = 1500; % lbf
Pe = 13; % psi
C_star_eff = 1;
C_F_eff = 1;
L_star = 35; % in
h_ch = 0.001/0.0254; % in
w_ch_min = 0.001/0.0254; % in
num_ch = 50;
w_rib = 0.001/0.0254; % in
D_c = 3.875; % in
t_ins = 0.001/0.0254; % in
t_out = 0.0625; % in

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

num_nodes = 100;
angle_conv = 40;
[A,D,M,P,T,w_ch,h_ch,D_h,w_rib,t_ins,t_out,step,pos,D_t,A_t,Pc,Pe,mdot,MW_g,gamma,mu_g,Cp_g,k_g,Pr_g] = geometry(Thrust,Pc,Pe,C_star,C_star_eff,C_F,C_F_eff,L_star,angle_conv,h_ch,w_rib,w_ch_min,MW_g,gamma,mu_g,Cp_g,k_g,Pr_g,t_ins,t_out,T_tot,num_nodes);
mdot_f = mdot*1/(1+O_F); %kg/s

ratio = 0.75;

% D_t = D_t*0.0254; % m
% A_t = A_t*0.0254^2; % m^2

[T_c,T_sat_c,P_c,q,T_wc,T_wr] = heatTransfer2D(Pc,M,P,A,D,D_h,w_ch,h_ch,w_rib,num_ch,t_ins,t_out,step,pos,gamma,mu_g,Cp_g,k_g,Pr_g,C_star,T,k_w,D_t,A_t,T_inlet,ratio,mdot_f,C_star_eff,res);

P_c = P_c/6894.76; % psi

Q = q(1:end-1).*D(1:end-1).*step*pi;
Loss = sum(Q) % W

figure (1)
plot(pos,D/2)
title('contour')
axis equal
figure(2)
hold on
plot(pos,T_c)
plot(pos,T_sat_c)
title('T_c')
legend('coolant','saturation')
figure(3)
plot(pos,P_c)
title('P_c')
figure(4)
plot(pos,q)
title('q')
figure(5)
hold on
plot(pos,T_wc)
plot(pos,T_wr)
hold off
title('T_w')
legend('channel','rib')