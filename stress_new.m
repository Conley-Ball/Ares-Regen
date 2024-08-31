function [v_m_stress] = stress_new(P_c,P,w_ch,t_ins,w_rib,D,t_out,pos,alpha_iw,E_iw,nu_iw,h_ch,T_ci,T_co,D_t,num_ch,l_div)

%Changes I made to main function (for Conley to update)
%Replaced stress function with this one
%set l_div as an output for the geometry function for axial stress if statement
%changed plot to v_m_stress
%commented out "stressTotaliw = stressTiw + stressP_hoop" line as its not no longer needed

for i = 1:length(pos)
q = P_c(i) - P(i);
r0 = D(i)/2+t_ins+h_ch+t_out; 
ri = D(i)/2;
r_bar = (r0+ri)/2;

sigma_s_i(i) = q*w_ch(i)/2*t_ins; %3 max shear stress in the inner wall (rib-wall joint)
sigma_b_i(i) = q/2*(w_ch(i).^2/t_ins); %4 max bending stress
%sigma_cr(i) = q/w_rib*(w_ch(i)+w_rib); %5 P is greater than p_c, compressive
sigma_t_r(i) = q/w_rib*w_ch(i); %6 p_c is greater than P, tensile (not being used in von mises?)
sigma_phi(i) = q*r_bar/(t_ins+t_out); %7 %circumferential stress inner wall
%need to add outer wall still but not super important as its not the max stress case
sigma_phi_o(i) = alpha_iw(i).*E_iw(i).*(T_ci(i)-T_co(i))./(1-nu_iw(i)).*(t_out./(t_ins+t_out)); %B4 average circumferential thermal stress
sigma_s_r(i) = ((w_ch(i)+w_rib).*alpha_iw(i).*E_iw(i).*(T_ci(i)-T_co(i)))./(5*w_rib.*(1-nu_iw(i))).*(t_out.*t_ins)/(t_out+t_ins)^.2; %B4 max thermal shear stress
sigma_r(i) = (alpha_iw(i).*E_iw(i).*(T_ci(i)-T_co(i)))./((1-nu_iw(i)).*r0.*w_rib./(w_ch(i)+w_rib)).*((t_ins+t_out)./(t_ins.*t_out)); %B4 thermal compressive
            %sigma_phi_hydro(i) = q/2*((d+h_ch+t_in+t_out)/(t_in+t_out));%B5
            %sigma_t_r_hydro(i) = q*(w_ch+w_rib)/w_rib; %B5
            %sigma_s_i_hydro(i) = q/2*(w_ch/t_ins); %B5
            %sigma_b_i_hydro(i) = q/2*(w_ch/t_ins)^2; %B6
sigma_total_hoop(i) = sigma_phi(i)+sigma_phi_o(i);
sigma_total_shear(i) = sigma_s_i(i)+sigma_s_r(i);


% axial stress was taken from Adams definiton and is not from any paper
if pos(i)>l_div
    sigma_a(i) = ((P(end)*pi.*(D(i)^2/4-D_t^2/4))./ (pi .* ((D(i)/2+t_ins+h_ch+t_out)^2 - D(i)^2/4) - (num_ch .* h_ch .* w_ch(i)) ));
else
    sigma_a(i) = 0;
end
% general 3d von-mises stress taken from wiki where axial, hoop, and bending are sigma11, sigma22, sigma33
% https://en.wikipedia.org/wiki/Von_Mises_yield_criterion
v_m_stress(i) = sqrt((1/2*((sigma_a(i)-sigma_total_hoop(i))^2+(sigma_total_hoop(i)-sigma_b_i(i))^2+(sigma_b_i(i)-sigma_a(i))^2))+(3*sigma_total_shear(i)^2));
%v_m_stress(i) = ((sigma_total_hoop(i))^2+sigma_a(i)^2 - (sigma_a(i)*sigma_total_hoop(i)))^0.5; %from adam
end
end
