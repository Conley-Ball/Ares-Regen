function [sigma_s_i,sigma_s_o,sigma_c_r,sigma_t_r,sigma_s_r,sigma_r,sigma_a,sigma_s_i_hydro,sigma_b_i_hydro,sigma_t_r_hydro,sigma_total_inner,sigma_total_outer,sigma_h_i_hydro,sigma_h_o_hydro] = stress(P_c,P,w_ch,t_ins,w_rib,D,t_out,pos,alpha_iw,alpha_ow,E_iw,E_ow,nu_iw,nu_ow,h_ch,T_ci,T_co,T_chg,T_cb,T_ct,D_t,num_ch,l_div)
sigma_h_i = zeros(1,length(pos));
sigma_b_i = zeros(1,length(pos));
sigma_phi_i = zeros(1,length(pos));
sigma_total_inner = zeros(1,length(pos));
sigma_h_o = zeros(1,length(pos));
sigma_b_o = zeros(1,length(pos));
sigma_phi_o = zeros(1,length(pos));
sigma_total_outer = zeros(1,length(pos));
sigma_s_i = zeros(1,length(pos));
sigma_s_o = zeros(1,length(pos));
sigma_c_r = zeros(1,length(pos));
sigma_t_r = zeros(1,length(pos));
sigma_s_r = zeros(1,length(pos));
sigma_r = zeros(1,length(pos));
sigma_a = zeros(1,length(pos));
sigma_h_i_hydro = zeros(1,length(pos));
sigma_h_o_hydro = zeros(1,length(pos));
sigma_s_i_hydro = zeros(1,length(pos));
sigma_b_i_hydro = zeros(1,length(pos));
sigma_t_r_hydro = zeros(1,length(pos));
for i = 1:length(pos)
q1 = P_c(i) - P(i);
q2 = P_c(i) - 62052.8156; %pa
r0(i) = D(i)/2+t_ins(i)+h_ch(i)+t_out(i); 
%% Inner wall
% Hoop stress in inner wall
sigma_h_i(i) = (q1*(D(i)/2+t_ins(i)/2))/t_ins(i);
% Maximum bending stress in the inner wall (4)
sigma_b_i(i) = q1/2*(w_ch(i).^2/t_ins(i).^2); % channel hydrostat if chamber pressure is zero
% Circumferential thermal stress in the inner wall (B3)
sigma_phi_i(i) = (alpha_iw(i)*(T_chg(i)- T_cb(i))*E_iw(i))/(2*(1-nu_iw(i)));
% Total inner
sigma_total_inner(i) = sigma_h_i(i)+sigma_b_i(i)+sigma_phi_i(i);

%% Outer wall
% Hoop stress in outer wall
sigma_h_o(i) = (q2*(D(i)/2+t_ins(i)+h_ch(i)+t_out(i)/2))/t_out(i);
% Maximum bending in the outer wall (4)
sigma_b_o(i) = q2/2*(w_ch(i).^2/t_out(i)^2); % channel hydrostat if chamber pressure is zero
% Circumferential thermal stress in the outer wall (B3)
sigma_phi_o(i) = (alpha_ow(i)*abs(T_ct(i)-T_co(i))*E_ow(i))/(1-nu_ow(i));
% Total outer
sigma_total_outer(i) = sigma_h_o(i)+sigma_b_o(i)+sigma_phi_o(i);

%% Rib
% Maximum shear stress in the inner wall is at the rib wall joint (3)
sigma_s_i(i) = q1*w_ch(i)/2*t_ins(i); % channel hydrostat if chamber pressure is zero
% Maximum shear stress in the otuer wall is at the rib wall joint (3)
sigma_s_o(i) = q2*w_ch(i)/2*t_out(i);


% Compressive or tensile stress in the ribs due to pressure difference between channel and chamber (5),(6)
if P(i)>P_c(i)
    sigma_c_r(i) = q1/w_rib(i)*(w_ch(i)+w_rib(i));
    sigma_t_r(i) = 0;
else
    sigma_c_r(i) = 0;
    sigma_t_r(i) = q1/w_rib(i)*w_ch(i); % channel hydrostat if chamber pressure is zero
end

% Largest shear stress if shear strain in ribs is negligible (B4)
sigma_s_r(i) = ((w_ch(i)+w_rib(i)).*alpha_iw(i).*E_iw(i).*(T_ci(i)-T_co(i)))./(5*w_rib(i).*(1-nu_iw(i))).*(t_out(i).*t_ins(i))/(t_out(i)+t_ins(i))^.2;
% Compressive loading in ribs due to imposed tempurature difference between inner and outer walls (B4)
sigma_r(i) = (alpha_iw(i).*E_iw(i).*(T_ci(i)-T_co(i)))./(((1-nu_iw(i)).*r0(i).*w_rib(i)./(w_ch(i)+w_rib(i))).*((t_ins(i)+t_out(i))./(t_ins(i).*t_out(i))));

%% Axial stress                                                                      
if pos(i)>l_div
    sigma_a(i) = ((P(end)*pi.*(D(i)^2/4-D_t^2/4))./ (pi .* ((D(i)/2+t_ins(i)+h_ch(i)+t_out(i))^2 - D(i)^2/4) - (num_ch .* h_ch(i) .* w_ch(i)) ));
else
    sigma_a(i) = 0;
end

%% Channel hydrostat is max shear, bending in wall, and tensile in ribs when chamber pressure is zero and q is max
q_max = P_c(1);
sigma_h_i_hydro(i) = (q_max*(D(i)/2+t_ins(i)/2))/t_ins(i);
sigma_h_o_hydro(i) = (q_max*(D(i)/2+t_ins(i)+h_ch(i)+t_out(i)/2))/t_out(i);
sigma_s_i_hydro(i) = q_max*w_ch(i)/2*t_ins(i); %max shear
sigma_b_i_hydro(i) = q_max/2*(w_ch(i).^2/t_ins(i).^2); %max bending
sigma_t_r_hydro(i) = q_max/w_rib(i)*w_ch(i); %max tensile in ribs
end
end
