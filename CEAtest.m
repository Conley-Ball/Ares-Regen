clc; clear; close all;

% RKT1=CEA('problem','rocket','equilibrium','fac','ma,kg/m^2',3,'o/f',2,'case','CEAM-RKT1','p(psi)',413,'pi/p',30,'reactants','fuel','C2H5OH(L)','wt%',75,'fuel','H2O(L)','wt%',25,'oxid','O2(L)','output','transport','mks','end','screen');


Pc = 313;
Pe = 13;
OF = 0.9;
eth_frac = 0.75;
C_star_eff = 0.94;
C_F_eff = 0.99;
Thrust = 2000;

gamma = sizing(Pc, Pe, OF, eth_frac,Thrust,C_star_eff,C_F_eff);


