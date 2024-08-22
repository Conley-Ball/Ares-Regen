function [E_iw,E_ow,E_avg,alpha_iw,alpha_ow,alpha_avg,nu_iw,nu_ow,nu_avg,k] = MatProperties_AlSi10Mg(T_iw,T_ow,T)
% Material Properties for Additively Manufactured AlSi10Mg
% Source: Ansys Material Database

%===PROPERTIES===
% Density [kg/m^3]
row_table = [295	2670;
843	1710];

% Elastic modulous [Pa]
E_table = [298	76600000000;
323	76100000000;
373	74300000000;
423	72700000000;
473	70600000000;
523	68900000000;
573	67000000000];

% Poison's Ratio
nu_table = [298	0.33;
323	0.33;
373	0.33;
423	0.33;
473	0.33;
523	0.33;
573	0.33;];

% Thermal Expansion [1/K]
alpha_table = [373	2.06E-05;
423	2.36E-05;
473	2.47E-05;
523	2.58E-05;
573	3.04E-05;
623	3.29E-05;
673	2.71E-05;
723	2.44E-05];

% Thermal Conductivity [W/mK]
k_table = [295	110;
323	111;
373	112;
423	114;
473	113;
523	109;
573	116;
673	116;
723	115;
773	109;
803	109;
843	101;
893	54;
913	48;
973	51];

% Specific heat [J/kg]
Cp_table = [295	91
423	1010;
548	1025;
698	1136;
823	1326];


% ===Interpolation===
% Density [kg/m^3]
row = interp1(row_table(:,1), row_table(:,2), T, 'spline');

% Elastic modulous [Pa]
E = interp1(E_table(:,1), E_table(:,2), T, 'spline');
E_iw = interp1(E_table(:,1), E_table(:,2), T_iw, 'spline');
E_ow = interp1(E_table(:,1), E_table(:,2), T_ow, 'spline');
E_avg = (E_iw+E_ow)/2;

% Poison's Ratio
nu = interp1(nu_table(:,1), nu_table(:,2), T, 'spline');
nu_iw = interp1(nu_table(:,1), nu_table(:,2), T_iw, 'spline');
nu_ow = interp1(nu_table(:,1), nu_table(:,2), T_ow, 'spline');
nu_avg = (nu_iw+nu_ow)/2;

% Thermal Expansion [1/K]
alpha = interp1(alpha_table(:,1), alpha_table(:,2), T, 'spline');
alpha_iw = interp1(alpha_table(:,1), alpha_table(:,2), T_iw, 'spline');
alpha_ow = interp1(alpha_table(:,1), alpha_table(:,2), T_ow, 'spline');
alpha_avg = (nu_iw+nu_ow)/2;

% Thermal Conductivity [W/mK]
k = interp1(k_table(:,1), k_table(:,2), T, 'spline');

% Specific heat [J/kg]
Cp = interp1(Cp_table(:,1), Cp_table(:,2), T, 'spline');

end