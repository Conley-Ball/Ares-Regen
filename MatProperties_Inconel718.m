function [E_iw,E_ow,E_avg,alpha_iw,alpha_ow,alpha_avg,nu_iw,nu_ow,nu_avg] = MatProperties_Inconel718(T_iw,T_ow)

% Material properties [Inconel 718]
% Equations derived from best fit curves on experimental data using Logger Pro
% _iw: inner wall, _ow: outer wall, _avg: average

% Elastic Modulus [Pa] 
% Source: https://www.hightempmetals.com/techdata/hitempInconel718data.php
E_iw = 2.540E+11 + (-4.861E+08)*T_iw + (2.036E+06)*T_iw.^2 + (-4121)*T_iw.^3 + (3.765)*T_iw.^4 + (-0.001282)*T_iw.^5;
E_ow = 2.540E+11 + (-4.861E+08)*T_ow + (2.036E+06)*T_ow.^2 + (-4121)*T_ow.^3 + (3.765)*T_ow.^4 + (-0.001282)*T_ow.^5;
E_avg = (E_iw+E_ow)/2;

% Thermal diffusivity [m^2/s]
% Source: "Thermophysical properties of Inconel 718 alloy A Sh Agazhanov, D A Samoshkin and Yu M Kozlovskii"
alpha_iw = 0.004469 + (-1.682E-05)*T_iw + (5.685E-08)*T_iw.^2 + (-7.305E-11)*T_iw.^3 + (4.190E-14)*T_iw.^4 + (-8.867E-18)*T_iw.^5;
alpha_ow = 0.004469 + (-1.682E-05)*T_ow + (5.685E-08)*T_ow.^2 + (-7.305E-11)*T_ow.^3 + (4.190E-14)*T_ow.^4 + (-8.867E-18)*T_ow.^5;
alpha_avg = (alpha_iw+alpha_ow)/2;

% Poisons Ratio
% Source: https://www.hightempmetals.com/techdata/hitempInconel718data.php
nu_iw = (-1836)*exp(-(T_iw - 646.0).^2/(-8.935E+04)^2) + 1836; % Gausian
nu_ow = (-1836)*exp(-(T_ow - 646.0).^2/(-8.935E+04)^2) + 1836;
nu_avg = (nu_iw+nu_ow)/2;

% % Thermal Conductivity [W/mK]
%k = (-4.607E-06)*T^2 + (0.02220)*T + (3.468);
% 
% Linear thermal expansion coefficient [1/K]
% a = 7.216E+06 + (3.680E+04)*T + (-73.32)*T^2 + (0.07016)*T^3 + (-2.215E-05)*T^4;

end