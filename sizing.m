function [mdot,A_t,D_t,A_e,D_e,C_star,C_F,gamma,MW_g,Cp_g,mu_g,k_g,T_thr,Pr_g] = sizing(Pc,Pe,OF,eth_ratio,T_inlet,Thrust,C_star_eff,C_F_eff)
    % Initial Guess
    mdot = 1;
    A_e = 1;
    A_t = 1;
    error = 1;
    % CEA Convergence
    while error > 1e-2
        RKT1=CEA('problem','rocket','equilibrium','fac','ma,kg/s',mdot,'o/f',OF,'p(psi)',Pc,'pi/p',Pc/Pe,'supsonic(ae/at)',A_e/A_t,'reactants','fuel','C2H5OH(L)','wt%',eth_ratio*100,'t(k)',T_inlet,'fuel','H2O(L)','wt%',(1-eth_ratio)*100,'t(k)',T_inlet,'oxid','O2(L)','wt%',100,'t(k)',90.0,'output','transport','mks','end');
        C_star = RKT1.output.eql.cstar(1);
        C_F = RKT1.output.eql.cf(end);
        gamma = RKT1.output.eql.gamma(end);
        mdot=(Thrust*4.448)/(C_star*C_star_eff*C_F*C_F_eff); %kg/s
        A_t = (C_star*C_star_eff*mdot)/(Pc*6895); %m^2
        D_t = sqrt(4*A_t/pi); %m
        M_e = sqrt((2/(gamma-1))*((Pc/Pe)^((gamma-1)/gamma)-1)); %exit mach number using isentropic relations
        error = A_e;
        A_e = A_t*((gamma+1)/2)^(-(gamma+1)/(2*gamma-2))*((1+((gamma-1)/2)*M_e^2)^((gamma+1)/(2*gamma-2)))/M_e; %exit area using isentropic relations (m^2)
        error = abs(error-A_e);
        D_e = sqrt(4*A_e/pi); %m
    end
    % Output
    gamma = RKT1.output.eql.gamma(1:4)';
    MW_g = RKT1.output.eql.mw(1:4)';
    Cp_g = RKT1.output.eql.cp_tran.froz(1:4)'*1000;
    mu_g = RKT1.output.eql.viscosity(1:4)'*1e-6;
    k_g = RKT1.output.eql.conduct.froz(1:4)';
    T_thr = RKT1.output.eql.temperature(3)'*C_star_eff^2;
    Pr_g = RKT1.output.eql.prandtl.froz(1:4)';
end