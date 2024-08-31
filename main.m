
function [v_m_stress, Yield_iw, T_chg] = main(w_ch_min, w_rib, t_ins, h_ch)

%% Inputs

% w_ch_min = 0.000635;
% w_rib = 0.000381;
% t_ins = 0.00042;
% h_ch= 0.000635;

%CEA Inputs
Pc = 413; % psia
Pe = 13.7; % psia
O_F = 1.1;
T_inlet = 300; % K

res = 0; % thermal resistance coating
Thrust = 1600; % lbf
%Pe = 13; % psi 
C_star_eff = 0.94;
C_F_eff = 0.90;
L_star = 28; % in
h_ch = h_ch/0.0254; % in
w_ch_min = w_ch_min/0.0254; % in
%num_ch = 50;
w_rib = w_rib/0.0254; % in
%D_c = 3.875; % in
t_ins = t_ins/0.0254; % in
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
T_thr = 2657.58*C_star_eff^2;


%% Solution

num_nodes = 300;
angle_conv = 40;
[A,D,M,P,T,w_ch,h_ch,D_h,w_rib,num_ch,t_ins,t_out,step,pos,D_t,A_t,Pc,Pe,mdot,MW_g,gamma,mu_g,Cp_g,k_g,Pr_g,id_th,id_c,l_div] = geometry(Thrust,Pc,Pe,C_star,C_star_eff,C_F,C_F_eff,L_star,angle_conv,h_ch,w_rib,w_ch_min,MW_g,gamma,mu_g,Cp_g,k_g,Pr_g,t_ins,t_out,T_thr,num_nodes);
mdot_f = mdot*1/(1+O_F); %kg/s

ratio = 0.75;

% D_t = D_t*0.0254; % m
% A_t = A_t*0.0254^2; % m^2
tic
[T_c,T_sat_c,P_c,q,T_chg,T_rhg,T_ci,T_ri,T_cb,T_rb,T_rt,T_ro,T_co,T_ct] = heatTransfer2D(Pc,M,A,D,D_h,w_ch,h_ch,w_rib,num_ch,t_ins,t_out,step,pos,gamma,mu_g,Cp_g,Pr_g,C_star,T,D_t,A_t,T_inlet,ratio,mdot_f,C_star_eff,res);
toc

% ===STRESS===
% Temperature Dependent Material properties
% Inner Wall
% T_ci = eval(T_ci);
[row_iw,E_iw,nu_iw,alpha_iw,k_iw,Cp_iw,Yield_iw] = MatProperties_Steel174ph(T_ci);
% Outer Wall
% T_co = eval(T_co);
[row_ow,E_ow,nu_ow,alpha_ow,k_ow,Cp_ow,Yield_ow] = MatProperties_Steel174ph(T_co);

% Stress function call
[v_m_stress] = stress_new(P_c,P,w_ch,t_ins,w_rib,D,t_out,pos,alpha_iw,E_iw,nu_iw,h_ch,T_ci,T_co,D_t,num_ch,l_div);
% Total Inner wall stress only
%stressTotaliw = stressTiw + stressP_hoop;
% ===========


P_c = P_c/6894.76; % psi

Q = q(1:end-1).*D(1:end-1).*step*pi;
% Loss = sum(Q) % W

% figure (1)
% plot(pos,D/2)
% title('contour')
% axis equal
% 
% figure(2)
% clf
% hold on
% plot(pos,T_c,'LineWidth',1,'Color','b')
% plot(pos,T_sat_c,'--','LineWidth',1,'Color','r')
% xline(0, '--', 'Exit','HandleVisibility','off')
% xline(pos(id_th), '--', 'Throat','HandleVisibility','off')
% xline(pos(id_c), '--', 'Chamber','HandleVisibility','off')
% hold off
% title('Coolant Temperature')
% xlabel('Axial Distance (m)')
% ylabel('Temperature (K)')
% legend('T_c','T_{saturation}')
% grid on
% 
% figure(3)
% clf
% hold on
% plot(pos,P_c,'LineWidth',1)
% xline(0, '--', 'Exit','HandleVisibility','off')
% xline(pos(id_th), '--', 'Throat','HandleVisibility','off')
% xline(pos(id_c), '--', 'Chamber','HandleVisibility','off')
% hold off
% title('Coolant Pressure')
% xlabel('Axial Distance (m)')
% ylabel('Pressure (psia)')
% grid on
% 
% figure(4)
% clf
% hold on
% plot(pos,q,'LineWidth',1)
% xline(0, '--', 'Exit','HandleVisibility','off')
% xline(pos(id_th), '--', 'Throat','HandleVisibility','off')
% xline(pos(id_c), '--', 'Chamber','HandleVisibility','off')
% hold off
% title('Heat Flux')
% xlabel('Axial Distance (m)')
% ylabel('Heat Flux (W/m^2)')
% grid on
% 
% 
% figure(5)
% clf
% hold on
% xline(pos(id_th), '--', 'Throat','HandleVisibility','off')
% xline(pos(id_c), '--', 'Chamber','HandleVisibility','off')
% plot(pos,T_chg,'LineWidth',1)
% plot(pos,T_rhg,'LineWidth',1)
% xline(0, '--', 'Exit','HandleVisibility','off')
% hold off
% title('Hot Wall Temperature')
% xlabel('Axial Distance (m)')
% ylabel('Temperature (K)')
% legend('Coolant Section','Rib Section')
% grid on
% 
% figure(6)
% clf
% hold on
% xline(pos(id_th), '--', 'Throat','HandleVisibility','off')
% xline(pos(id_c), '--', 'Chamber','HandleVisibility','off')
% plot(pos,v_m_stress/6895000,'LineWidth',1)
% plot(pos,Yield_iw/6895000,'LineWidth',1)
% xline(0, '--', 'Exit','HandleVisibility','off')
% hold off
% title('Stress')
% xlabel('Axial Distance (m)')
% ylabel('Stress (ksi)')
% legend('vonMises','Inner Wall Yield')
% grid on
