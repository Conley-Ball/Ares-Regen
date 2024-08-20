clc; clear; close all;

%% Inputs

Pc = 400; % psi
gamma = 1.15;
MW_g = 21.83; % g/mol
Cp_g = 3.346; %kJ/kgK
C_star = 1627.4; % m/s
mu_g = 0.74167; %Pa-s
T_aw = 2861; % K
k_w = 20; %W/m K
T_inlet = 550; % K
mdot_f = 1.5;

% Geometry

%Chamber dimensions
l_c = 6.5; %in
th_in = 0.022638; %in
th_out = 0.0625; %in
d_c_in = 3.875; %in
l_c_ch = 0.2; %in
w_c_ch = 0.0984;%in
nodes_c = 20;
%Converging dimensions
l_conv = 2.335; %in
angle_conv = 35; %deg
nodes_conv = 40;
%Diversing dimensions
l_div = 2.626; %in 
d_t = 1.7938; %in
d_e = 3.724; %in
l_div_ch = 0.16; %in
w_div_ch = 0.0984; %in
nodes_div = 20;
%Additional Paramenters
th_l_ch = 0.0689; %in
th_w_ch = 0.041; %in

%% Solution

% [A,D,M,P,l_ch,w_ch,D_h,step,pos] = geometryAutomatic(Thrust,Pc,O_F,C_star,P_a,conv_angle,t_ins,t_out,L_star,gamma);
% 
% [A,D,M,P,l_ch,w_ch,D_h,step,pos] = geometryManual(Pc,conv_angle,t_ins,t_out,A_c,L_c,A_th,L_conv,A_e,L_div,l_c_ch,w_c_ch,l_th_ch,w_th_ch,l_e_ch,w_e_ch,gamma);

[M,P,D_h,A,D,l_ch,w_ch,step,pos] = geometry(Pc,gamma,l_c,th_in,th_out,d_c_in,l_c_ch,w_c_ch,nodes_c,l_conv,angle_conv,nodes_conv,l_div,d_t,d_e,l_div_ch,w_div_ch,nodes_div,th_l_ch,th_w_ch);

A = A*0.0254^2;
% D = D*0.0254;
% D_h = D_h*0.0254;
% l_ch = l_ch*0.0254;
% w_ch = w_ch*0.0254;
step = step*0.0254;
pos = pos*0.0254;
P = P*6894.76;
Pc = Pc*6894.76;

ratio = 0.75;
num = 40;

[T_c,P_c,q] = heatTransfer2D(Pc,M,P,A,D,D_h,l_ch,w_ch,num,step,pos,gamma,MW_g,mu_g,Cp_g,C_star,T_aw,k_w,T_inlet,ratio,mdot_f);

