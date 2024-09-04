function [v_m_stress, Yield_iw, T_chg, P_c] = main_function(w_ch_min, w_rib, t_ins, h_ch_th, h_ch_c, h_ch_e, Pc, Pe, O_F, T_inlet, res, Thrust, C_star_eff, C_F_eff, L_star, t_out, ratio, AR, C_star, C_F, gamma, MW_g, Cp_g, mu_g, k_g, T_thr, Pr_g)



% res = 0.00005/20; % thermal resistance coating t/k



% h_ch_c = 1.5*h_ch_th; % in
% h_ch_e = 1.5*h_ch_th; % in







%% Solution

num_nodes = 300;
angle_conv = 40;
[A,D,M,P,T,w_ch,h_ch,D_h,w_rib,num_ch,t_ins,t_out,step,pos,D_t,A_t,Pc,Pe,mdot,MW_g,gamma,mu_g,Cp_g,k_g,Pr_g,id_th,id_c,l_div] = geometry(Thrust,Pc,Pe,C_star,C_star_eff,C_F,C_F_eff,L_star,angle_conv,h_ch_th,h_ch_c,h_ch_e,w_rib,w_ch_min,MW_g,gamma,mu_g,Cp_g,k_g,Pr_g,t_ins,t_out,T_thr,num_nodes);
mdot_f = mdot*1/(1+O_F); %kg/s

% D_t = D_t*0.0254; % m
% A_t = A_t*0.0254^2; % m^2
tic
[T_c,T_sat_c,P_c,q,T_chg,T_rhg,T_ci,T_ri,T_cb,T_rb,T_rt,T_ro,T_co,T_ct,h_chg,h_rhg] = heatTransfer2D(Pc,M,A,D,D_h,w_ch,h_ch,w_rib,num_ch,t_ins,t_out,step,pos,gamma,mu_g,Cp_g,Pr_g,C_star,T,D_t,A_t,T_inlet,ratio,mdot_f,C_star_eff,res);
toc

% ===STRESS===
% Temperature Dependent Material properties
% Inner Wall
[row_iw,E_iw,nu_iw,alpha_iw,k_iw,Cp_iw,Yield_iw] = MatProperties_AlSi10Mg(T_ci);
% Outer Wall
[row_ow,E_ow,nu_ow,alpha_ow,k_ow,Cp_ow,Yield_ow] = MatProperties_AlSi10Mg(T_co);

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
% plot(pos,T_cb,'LineWidth',1)
% plot(pos,T_ci,'LineWidth',1)
% plot(pos,T_ct,'LineWidth',1)
% plot(pos,T_co,'LineWidth',1)
% xline(0, '--', 'Exit','HandleVisibility','off')
% hold off
% title('Hot Wall Temperature')
% xlabel('Axial Distance (m)')
% ylabel('Temperature (K)')
% legend('T_{w}','T_{cb}','T_{ci}','T_{ct}','T_{co}')
% grid on
% 
% figure(6)
% clf
% hold on
% xline(pos(id_th), '--', 'Throat','HandleVisibility','off')
% xline(pos(id_c), '--', 'Chamber','HandleVisibility','off')
% % plot(pos,T_chg,'LineWidth',1)
% plot(pos,h_chg,'LineWidth',1)
% plot(pos,h_rhg,'LineWidth',1)
% xline(0, '--', 'Exit','HandleVisibility','off')
% hold off
% title('Heat Transfer Coefficients')
% xlabel('Axial Distance (m)')
% ylabel('heat')
% legend('h_{chg}','h_{rhg}')
% grid on
% 
% figure(7)
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
end