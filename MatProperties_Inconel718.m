function [row,E,nu,alpha,k,Cp,Yield] = MatProperties_Inconel718(T)
% Material Properties for Additively Manufactured Inconel 718
% Source: Ansys Material Database
%{ There are memory and performance benefits to using "griddedInterpolant" objects over the "interp" functions. 
% griddedInterpolant offers substantial performance improvements for repeated queries of the interpolant object, 
% whereas the interp functions perform a new calculation each time they are called. Also, griddedInterpolant 
% stores the sample points in a memory-efficient format (as a compact grid) and is multithreaded to take 
% advantage of multicore computer processors.
%}

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

% Yield Strength [Pa]
Yield_table = [294	648000000;
811	    558000000;
1089	338000000;
1255	90000000;
1366	34000000];


% ===Interpolation===
% Density [kg/m^3]
row_interp = griddedInterpolant(row_table(:,1), row_table(:,2));
row = row_interp(T);

% Elastic modulous [Pa]
E_interp = griddedInterpolant(E_table(:,1), E_table(:,2));
E = E_interp(T);

% Poison's Ratio
nu_interp = griddedInterpolant(nu_table(:,1), nu_table(:,2));
nu = nu_interp(T);

% Thermal Expansion [1/K]
alpha_interp = griddedInterpolant(alpha_table(:,1), alpha_table(:,2));
alpha = alpha_interp(T);

% Thermal Conductivity [W/mK]
k_interp = griddedInterpolant(k_table(:,1), k_table(:,2));
k = k_interp(T);

% Specific heat [J/kg]
Cp_interp = griddedInterpolant(Cp_table(:,1), Cp_table(:,2));
Cp = Cp_interp(T);

% Yield Strength [Pa]
Yield_interp = griddedInterpolant(Yield_table(:,1), Yield_table(:,2));
Yield = Yield_interp(T);


end

