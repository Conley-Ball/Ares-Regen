function [row,E,nu,alpha,k,Cp,Yield] = MatProperties_AlSi10Mg(T)
% Material Properties for Additively Manufactured AlSi10Mg
% Source: Ansys Material Database
%{ There are memory and performance benefits to using "griddedInterpolant" objects over the "interp" functions. 
% griddedInterpolant offers substantial performance improvements for repeated queries of the interpolant object, 
% whereas the interp functions perform a new calculation each time they are called. Also, griddedInterpolant 
% stores the sample points in a memory-efficient format (as a compact grid) and is multithreaded to take 
% advantage of multicore computer processors.
%}

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

% Yield Strength [Pa]
Yield_table = [298	251000000;
373	232000000;
423	221000000;
473	197000000;
523	148000000];


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