clc; clear; close all;

% ===INPUTS===

material =              'aluminum'; % Select from: 'steel' 'aluminum' 'inconel'
Pc =                    313.7000; % psia
%% 
O_F =                   0.85;
Thrust =                2000; % lbf
fos =                   600;


roughness =             19e-6; % m
stiffness =             0.2;
L_star =                25; % in
angle_conv =            40; % deg
pitch =                 0; % deg
throat_only =           false;
max_angle =             45; % deg
Pe =                    13.7; % psia
Cd =                    0.7;

res =                   0.00005/1; % thermal resistance coating t/k
% res =                   0;
C_star_eff =            0.94;
C_F_eff =               0.99;
T_inlet =               300; % K
T_CEA_f =               370; % K
T_initial =             300; % K
eth_ratio =             0.75;

h_ch =                  [0.001 0.001 0.001 0.00175]/0.0254;  % higher raises Q, lower lowers Q
w_ch_min =              0.001/0.0254;                      % higher raises Q, lower lowers Q
w_rib =                 [0.001 0.001 0.001 0.001]/0.0254;                      
t_ins =                 [0.002 0.001 0.001 0.001]/0.0254;  % higher lowers Q, lower raises Q (and lowers FOS)
t_out =                 [0.002 0.0015 0.0015 0.0015]/0.0254; % in
fillet =                0.00025; % m radius



% ===CEA===

[mdot,A_t,D_t,A_e,D_e,C_star,C_F,gamma,MW_g,Cp_g,mu_g,k_g,T_thr,Pr_g] = sizing(Pc,Pe,O_F,eth_ratio,T_inlet,Thrust,C_star_eff,C_F_eff,T_CEA_f);
fprintf('CEA Finished\n')

% ===GEOMETRY===

num_nodes = 1500; % Station Resolution
[A,D,M,P,T,w_ch,D_h,w_rib,num_ch,t_ins,t_out,step,pos,D_t,A_t,Pc,Pe,mdot,A_e,D_e,MW_g,gamma,mu_g,Cp_g,k_g,Pr_g,id_th,id_c,l_div,h_ch,psi,phi,dtheta] = geometry(Thrust,Pc,Pe,mdot,A_t,D_t,A_e,D_e,C_star,C_star_eff,C_F,C_F_eff,L_star,angle_conv,w_rib,w_ch_min,MW_g,gamma,mu_g,Cp_g,k_g,Pr_g,t_ins,t_out,T_thr,num_nodes,h_ch,pitch,max_angle,throat_only);
mdot_f = mdot*1/(1+O_F); %kg/s
fprintf('Geometry Finished\n')

% ===HEAT TRANSFER===

fprintf('Heat Transfer       \n')
ch_resolution = 1e-5; % m
if res > 0
[T_c,T_sat_c,P_c,q,T_chg,T_rhg,T_ci,T_ri,T_cb,T_rb,T_rt,T_ro,T_co,T_ct,h_chg,h_rhg,h_ch,t_ins,q_crit,q_crit2,T_cw,T_rw,T_aw,h_c] = heatTransfer2DTC(Pc,M,A,D,D_h,w_ch,w_rib,num_ch,t_ins,t_out,step,pos,gamma,mu_g,Cp_g,Pr_g,C_star,T,D_t,A_t,T_inlet,eth_ratio,mdot_f,C_star_eff,res,num_ch,l_div,fos,P,fillet,stiffness,material,roughness,ch_resolution,h_ch,phi);
else
[T_c,T_sat_c,P_c,q,T_chg,T_rhg,T_ci,T_ri,T_cb,T_rb,T_rt,T_ro,T_co,T_ct,h_chg,h_rhg,h_ch,t_ins,q_crit,q_crit2] = heatTransfer2D(Pc,M,A,D,D_h,w_ch,w_rib,num_ch,t_ins,t_out,step,pos,gamma,mu_g,Cp_g,Pr_g,C_star,T,D_t,A_t,T_inlet,eth_ratio,mdot_f,C_star_eff,res,num_ch,l_div,fos,P,fillet,stiffness,material,roughness,ch_resolution,h_ch,phi);
end
fprintf('\b\b\b\b\b\b\bFinished\n')

% ===STRESS===

% Temperature Dependent Material properties
[E_iw,nu_iw,alpha_iw,k_iw,Cp_iw,Yield_iw] = materialProperties(T_ci,material);
[E_ow,nu_ow,alpha_ow,k_ow,Cp_ow,Yield_ow] = materialProperties(T_co,material);
[E_hydro,nu_hydro,alpha_hydro,k_hydro,Cp_hydro,Yield_hydro] = materialProperties(T_initial*ones(size(pos)),material);

if res == 0
    T_cw = T_chg;
end
[sigma_s_i,sigma_s_o,sigma_c_r,sigma_t_r,sigma_s_r,sigma_r,sigma_a,sigma_s_i_hydro,sigma_b_i_hydro,sigma_t_r_hydro,sigma_total_inner,sigma_total_outer,sigma_h_i_hydro,sigma_h_o_hydro,sigma_h_i,sigma_phi_i] = stress(P_c,P,w_ch,t_ins,w_rib,D,t_out,pos,alpha_iw,alpha_ow,E_iw,E_ow,nu_iw,nu_ow,h_ch,T_ci,T_co,T_cw,T_cb,T_ct,D_t,num_ch,l_div);
% Total Inner wall stress only
%stressTotaliw = stressTiw + stressP_hoop;
fprintf('Stress Finished\n')

% ===INJECTOR===

rho_water = py.CoolProp.CoolProp.PropsSI('D','T', T_c(end), 'P', P_c(end), 'water');
rho_ethanol = py.CoolProp.CoolProp.PropsSI('D','T', T_c(end), 'P', P_c(end), 'ethanol');
rho_c = ((1-eth_ratio)/rho_water+eth_ratio/rho_ethanol)^-1;
A_fuel = mdot_f/(Cd*(2*rho_c*0.2*Pc)^0.5);
A_ox = (mdot-mdot_f)/(0.75*(2*1141*0.2*Pc)^0.5);

% ==SOLIDWORKS==

curveHelper(D, t_ins, h_ch, t_out, pos, dtheta, psi, w_rib, w_ch, num_ch)
fprintf('Geometry Exported\n')

% ===GRAPHS===

f1 = figure(1);
hold on
plot(pos,D/2,'Color','k')
plot((pos+(t_ins).*sin(-psi)),(D/2+(t_ins).*cos(psi)),'Color','k')
plot((pos+(t_ins+h_ch).*sin(-psi)),(D/2+(t_ins+h_ch).*cos(psi)),'Color','k')
plot((pos+(t_ins+h_ch+t_out).*sin(-psi)),(D/2+(t_ins+h_ch+t_out).*cos(psi)),'Color','k')
plot(pos,-D/2,'Color','k')
plot((pos+(t_ins).*sin(-psi)),-(D/2+(t_ins).*cos(psi)),'Color','k')
plot((pos+(t_ins+h_ch).*sin(-psi)),-(D/2+(t_ins+h_ch).*cos(psi)),'Color','k')
plot((pos+(t_ins+h_ch+t_out).*sin(-psi)),-(D/2+(t_ins+h_ch+t_out).*cos(psi)),'Color','k')
hold off
axis equal
title('Contour')
grid on
% f1.Position = [200,200,600,400];

f2 = figure(2);
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
% f2.Position = [200,200,600,400];

f3 = figure(3);
clf
hold on
plot(pos,P_c/6894.76,'LineWidth',1)
xline(0, '--', 'Exit','HandleVisibility','off')
xline(pos(id_th), '--', 'Throat','HandleVisibility','off')
xline(pos(id_c), '--', 'Chamber','HandleVisibility','off')
hold off
title('Coolant Pressure')
xlabel('Axial Distance (m)')
ylabel('Pressure (psia)')
grid on
% f3.Position = [200,200,600,400];

f4 = figure(4);
clf
hold on
plot(pos,q,'LineWidth',1)
plot(pos,q_crit,'--','LineWidth',1)
plot(pos,q_crit2,'--','LineWidth',1)
xline(0, '--', 'Exit','HandleVisibility','off')
xline(pos(id_th), '--', 'Throat','HandleVisibility','off')
xline(pos(id_c), '--', 'Chamber','HandleVisibility','off')
hold off
title('Heat Flux')
xlabel('Axial Distance (m)')
ylabel('Heat Flux (W/m^2)')
legend('Heat Flux (W/m^2)','Critical Heat Flux(W/m^2)','Rosenhow Critical Heat Flux(W/m^2)')
grid on
% f4.Position = [200,200,500,400];


f5 = figure(5);
clf
hold on
xline(pos(id_th), '--', 'Throat','HandleVisibility','off')
xline(pos(id_c), '--', 'Chamber','HandleVisibility','off')
plot(pos,T_chg,'LineWidth',1)
plot(pos,T_cb,'LineWidth',1)
plot(pos,T_ci,'LineWidth',1)
plot(pos,T_ct,'LineWidth',1)
plot(pos,T_co,'LineWidth',1)
if res > 0
plot(pos,T_cw,'LineWidth',1)
plot(pos,T_rw,'LineWidth',1)
end
xline(0, '--', 'Exit','HandleVisibility','off')
hold off
title('Hot Wall Temperature')
xlabel('Axial Distance (m)')
ylabel('Temperature (K)')
legend('T_{w}','T_{cb}','T_{ci}','T_{ct}','T_{co}','T_{cw}','T_{rw}')
grid on
% f5.Position = [200,200,500,400];

f6 = figure(6);
clf
hold on
xline(pos(id_th), '--', 'Throat','HandleVisibility','off')
xline(pos(id_c), '--', 'Chamber','HandleVisibility','off')
% plot(pos,T_chg,'LineWidth',1)
plot(pos,h_chg,'LineWidth',1)
plot(pos,h_rhg,'LineWidth',1)
xline(0, '--', 'Exit','HandleVisibility','off')
hold off
title('Heat Transfer Coefficients')
xlabel('Axial Distance (m)')
ylabel('heat')
legend('h_{chg}','h_{rhg}')
grid on
% f6.Position = [200,200,500,400];

f7 = figure(7);
clf
hold on
xline(pos(id_th), '--', 'Throat','HandleVisibility','off')
xline(pos(id_c), '--', 'Chamber','HandleVisibility','off')
plot(pos,sigma_total_inner/6895000,'LineWidth',1)
plot(pos,Yield_iw/6895000,'--','LineWidth',1)
plot(pos,sigma_total_outer/6895000,'LineWidth',1)
plot(pos,Yield_ow/6895000,'--','LineWidth',1)
xline(0, '--', 'Exit','HandleVisibility','off')
hold off
title('Stress')
xlabel('Axial Distance (m)')
ylabel('Stress (ksi)')
legend('Total inner wall stress','Inner Wall Yield','Total outer wall stress','Outer Wall Yield')
grid on
% f7.Position = [200,200,500,400];

f8 = figure(8);
clf
hold on
xline(pos(id_th), '--', 'Throat','HandleVisibility','off')
xline(pos(id_c), '--', 'Chamber','HandleVisibility','off')
plot(pos,w_ch,'--','LineWidth',1.5)
plot(pos,h_ch,'--','LineWidth',1.5)
plot(pos,w_rib,'--','LineWidth',1.5)
plot(pos,t_ins,'--','LineWidth',1.5)
plot(pos,t_out,'--','LineWidth',1.5)
xline(0, '--', 'Exit','HandleVisibility','off')
hold off
title('Feature Sizes')
xlabel('Axial Distance (m)')
ylabel('Feature Size (m)')
ylim([0,0.004])
legend('Channel Width','Channel Height','Rib Width','Inner Wall Thickness','Outer Wall Thickness')
grid on
% f8.Position = [200,200,600,400];

f10 = figure(10);
clf
hold on
xline(pos(id_th), '--', 'Throat','HandleVisibility','off')
xline(pos(id_c), '--', 'Chamber','HandleVisibility','off')
plot(pos,sigma_h_i_hydro/6895000,'LineWidth',1,'DisplayName','Hydrostatic Inner Hoop')
plot(pos,sigma_h_o_hydro/6895000,'LineWidth',1,'DisplayName','Hydrostatic Outer Hoop')
plot(pos,sigma_s_i_hydro/6895000,'LineWidth',1,'DisplayName','Hydrostatic Inner Shear')
plot(pos,sigma_b_i_hydro/6895000,'LineWidth',1,'DisplayName','Hydrostatic Inner Bending')
plot(pos,sigma_t_r_hydro/6895000,'LineWidth',1,'DisplayName','Hydrostatic Rib Tensile')
plot(pos,Yield_hydro/6895000, '--','LineWidth',1,'DisplayName','Hydrostatic Yield')
xline(0, '--', 'Exit','HandleVisibility','off')
hold off
title('Hydrostatic Stress')
xlabel('Axial Distance (m)')
ylabel('Stress (ksi)')
legend()
grid on


f11 = figure(11);
clf
hold on
xline(pos(id_th), '--', 'Throat','HandleVisibility','off')
xline(pos(id_c), '--', 'Chamber','HandleVisibility','off')
plot(pos,sigma_total_inner/6895000,'LineWidth',1)
plot(pos,sigma_h_i/6895000,'LineWidth',1)
plot(pos,sigma_phi_i/6895000,'LineWidth',1)
plot(pos,Yield_iw/6895000,'--','LineWidth',1)
xline(0, '--', 'Exit','HandleVisibility','off')
hold off
title('Inner wall Stress')
xlabel('Axial Distance (m)')
ylabel('Stress (ksi)')
legend('Total inner wall stress','Inner wall hoop','Inner wall thermal','Yield inner wall')
grid on



%ANSYS Data
Data_h_chg = [pos(:), h_chg(:)];
Data_T_aw = [pos(:), T_aw(:)];
Data_h_c = [pos(:), h_c(:)];
Data_T_c = [pos(:), T_c(:)];
Data_P = [pos(:), P(:)];
Data_P_c = [pos(:), P_c(:)];