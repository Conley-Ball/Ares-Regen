clc; clear; close all;

%% Inputs

Pc = 400; % psi
gamma = 1.15;
MW_g = 21.83; % g/mol
Cp_g = 3.346; %kJ/kgK

C_star = 1627.4; % m/s
mu_g = 0.74167; %Pa-s
T_cns = 2861; % K
k_w = 20; %W/m K
T_inlet = 300; % K


mdot_f = 1.5;
Thrust = 1500;
Pe = 13;
C_star_eff = 0.94;
C_F = 1.4;
C_F_eff = 0.9;
L_star = 35;
w_ch = 0.0984;
num_ch = 42;
l_rib = 0.02;
D_c = 3.875;
t_ins = 0.001;
t_out = 0.0625*0.0254;

%% Solution

% [A,D,M,P,l_ch,w_ch,D_h,step,pos] = geometryAutomatic(Thrust,Pc,O_F,C_star,P_a,conv_angle,t_ins,t_out,L_star,gamma);
% 

% [C_star,C_F,gamma,MW_g,Cp_g,mu_g,T_cns] = CEAproperties(Pc,Pe,O_F,eth_ratio)

[A,D,M,P,l_ch,w_ch,D_h,step,pos,D_t,A_t] = geometry(Thrust,Pc,Pe,C_star,C_star_eff,C_F,C_F_eff,L_star,w_ch,num_ch,l_rib,D_c,gamma);

A = A*0.0254^2;
D = D*0.0254;
D_h = D_h*0.0254;
l_ch = l_ch*0.0254;
l_rib = l_rib*0.0254;
w_ch = w_ch*0.0254;
step = step*0.0254;
pos = pos*0.0254;
P = P*6894.76;
Pc = Pc*6894.76;

ratio = 0.75;
num = 40;

D_t = D_t*0.0254;
A_t = A_t*0.0254^2;

[T_c,P_cool,q] = heatTransfer2D(Pc,M,P,A,D,D_h,l_ch,w_ch,l_rib,num,t_ins,t_out,step,pos,gamma,MW_g,mu_g,Cp_g,C_star,T_cns,k_w,D_t,A_t,T_inlet,ratio,mdot_f);