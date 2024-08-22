function [E_iw,E_ow,E_avg,alpha_iw,alpha_ow,alpha_avg,nu_iw,nu_ow,nu_avg,k] = MatProperties_Inconel718(T_iw,T_ow,T)
% Material Properties for Additively Manufactured Inconel 718
% Source: Ansys Material Database

%===PROPERTIES===
% Density [kg/m^3]
row_table = [293	8220;
500	8121;
700	8048;
1000 7961;
1200 7875;
1400 7787;
1533 7733;
1617 7579;
1723 7488;
1800 7488;
2100 7341];

% Elastic modulous [Pa]
E_table = [294	165000000000;
810	    152000000000;
1088	110000000000;
1255	55000000000;
1366	34000000000];

% Poison's Ratio
nu_table = [294	0.3;
810	    0.28;
1088	0.323;
1255	0.368;
1366	0.4];

% Thermal Expansion [1/K]
alpha_table = [700	1.44E-05;
811	    1.49E-05;
922	    1.54E-05;
1144	1.75E-05;
1366	1.83E-05;
2100	1.83E-05];

% Thermal Conductivity [W/mK]
k_table = [295	11.9;
506	    13.7;
721	    16.9;
930	    21.7;
1139	25.6;
1352	22.9;
1562	19.1;
1773	17.7];

% Specific heat [J/kg]
Cp_table = [293	421;
373	    442;
473	    453;
573	    472;
673	    481;
773	    502;
873	    527;
973	    562;
1073	606;
1123	628;
1173	636;
1273	647;
1373	651;
1773	652];


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

