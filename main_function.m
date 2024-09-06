function [v_m_stress, Yield_iw, T_chg, P_c] = main_function(w_ch_min, w_rib, t_ins, h_ch, Pc, Pe, O_F, T_inlet, res, Thrust, C_star_eff, C_F_eff, L_star, t_out, ratio, AR, C_star, C_F, gamma, MW_g, Cp_g, mu_g, k_g, T_thr, Pr_g)
%% Solution
num_nodes = 300;
angle_conv = 40;
[A,D,M,P,T,w_ch,h_ch,D_h,w_rib,num_ch,t_ins,t_out,step,pos,D_t,A_t,Pc,Pe,mdot,MW_g,gamma,mu_g,Cp_g,k_g,Pr_g,id_th,id_c,l_div] = geometry(Thrust,Pc,Pe,C_star,C_star_eff,C_F,C_F_eff,L_star,angle_conv,h_ch,w_rib,w_ch_min,MW_g,gamma,mu_g,Cp_g,k_g,Pr_g,t_ins,t_out,T_thr,num_nodes);
mdot_f = mdot*1/(1+O_F); %kg/s


[T_c,T_sat_c,P_c,q,T_chg,T_rhg,T_ci,T_ri,T_cb,T_rb,T_rt,T_ro,T_co,T_ct,h_chg,h_rhg] = heatTransfer2D(Pc,M,A,D,D_h,w_ch,h_ch,w_rib,num_ch,t_ins,t_out,step,pos,gamma,mu_g,Cp_g,Pr_g,C_star,T,D_t,A_t,T_inlet,ratio,mdot_f,C_star_eff,res);

[row_iw,E_iw,nu_iw,alpha_iw,k_iw,Cp_iw,Yield_iw] = MatProperties_AlSi10Mg(T_ci);

[row_ow,E_ow,nu_ow,alpha_ow,k_ow,Cp_ow,Yield_ow] = MatProperties_AlSi10Mg(T_co);

[v_m_stress] = stress_new(P_c,P,w_ch,t_ins,w_rib,D,t_out,pos,alpha_iw,E_iw,nu_iw,h_ch,T_ci,T_co,D_t,num_ch,l_div);

P_c = P_c/6894.76; % psi

Q = q(1:end-1).*D(1:end-1).*step*pi;
% Loss = sum(Q) % W
end