clc; clear; close all;

%% Inputs

Pc = 400; % psi
O_F = 1.1;
gamma = 1.21;
Cp_g = 3346; %J/kgK

C_star = 1627.4; % m/s
mu_g = 0.74167e-4; %Pa-s
T_cns = 2861; % K
k_w = 20; %W/m K
T_inlet = 300; % K

res = 0; % thermal resistance coating
Thrust = 6672; % N
Pe = 13; % psi
C_star_eff = 0.94;
C_F = 1.4;
C_F_eff = 0.9;
L_star = 35; % in
w_ch = 0.001/0.0254; % in
num_ch = 60;
l_rib = 0.001/0.0254; % in
D_c = 3.875; % in
t_ins = 0.001; % m
t_out = 0.0625*0.0254; % m

%% Solution

% [A,D,M,P,l_ch,w_ch,D_h,step,pos] = geometryAutomatic(Thrust,Pc,O_F,C_star,P_a,conv_angle,t_ins,t_out,L_star,gamma);
% 

% [C_star,C_F,gamma,MW_g,Cp_g,mu_g,T_cns] = CEAproperties(Pc,Pe,O_F,eth_ratio)

[A,D,M,P,l_ch,w_ch,D_h,step,pos,D_t,A_t,mdot] = geometry(Thrust,Pc,Pe,C_star,C_star_eff,C_F,C_F_eff,L_star,w_ch,num_ch,l_rib,D_c,gamma);

mdot_f = mdot*1/(1+O_F);
A = A*0.0254^2; % m^2
D = D*0.0254; % m
D_h = D_h*0.0254; % m
l_ch = l_ch*0.0254; % m
l_rib = l_rib*0.0254; % m
w_ch = w_ch*0.0254; % m
step = step*0.0254; % m
pos = pos*0.0254; % m
P = P*6894.76; % Pa
Pc = Pc*6894.76; % Pa

ratio = 0.75;

D_t = D_t*0.0254; % m
A_t = A_t*0.0254^2; % m^2

[T_c,T_sat_c,P_c,q,T_wc,T_wr] = heatTransfer2D(Pc,M,P,A,D,D_h,l_ch,w_ch,l_rib,num_ch,t_ins,t_out,step,pos,gamma,mu_g,Cp_g,C_star,T_cns,k_w,D_t,A_t,T_inlet,ratio,mdot_f,C_star_eff,res);

P_c = P_c/6894.76; % psi

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