function [E_iw,E_ow,E_avg,alpha_iw,alpha_ow,alpha_avg,nu_iw,nu_ow,nu_avg,k] = MatProperties_AlSi10Mg(T_iw,T_ow,T)
% Material Properties for Additively Manufactured Steel174ph
% Source: Ansys Material Database

%===PROPERTIES===
% Density [kg/m^3] (a
row_table = [295	7790];

% Elastic modulous [Pa]
E_table = [298	204000000000;
494	    195000000000;
580	    187000000000;
650	    182000000000;
728	    176000000000;
798	    168000000000;
885	    153000000000;
957	    142000000000;
1067	129000000000;
1162	117000000000];

% Poison's Ratio
nu_table = [298	0.291;
494	    0.295;
580	    0.296;
650	    0.305;
728	    0.316;
798	    0.309;
885	    0.322;
957	    0.332;
1067	0.348;
1162	0.361];

% Thermal Expansion [1/K]
alpha_table = [473	1.1E-05;
573	1.14E-05;
673	1.18E-05;
773	1.2E-05];

% Thermal Conductivity [W/mK]
k_table = [294	15.2;
461	18;
627	19.9;
794	20.6;
961	28.1];

% Specific heat [J/kg]
Cp_table = [350	475;
375	488.4;
400	497.3;
425	505.4;
450	514.7;
460	522;
475	522.8;
500	533.5;
525	547.7;
550	554.8;
600	563.6;
625	572.6;
650	580.5;
675	588.5;
700	597.8;
725	608;
750	634.3;
775	659.4;
795	681.2;
800	687.5;
825	724.6;
850	740.9;
875	757.6;
900	767.2;
925	778.9];


% ===Interpolation===
% Density [kg/m^3]
row = row_table(1,2);

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
