function F = solveFunction(T_w,w,b,t,t_r,l,epsilon,eta,T_aw,M,T,A,gamma,Cp_g,mu_g,Pr_g,Pc,C_star,C_star_eff,D_t,R,A_t,res,h_c,T_c,kT_chg,kT_rhg,kT_ci,kT_ri)
    sigma_c = ( (0.5*(T_w(1)/T)*(1+(gamma-1)/2*M^2)+0.5)^0.68*(1+(gamma-1)/2*M^2)^0.12 )^-1; % []
    sigma_r = ( (0.5*(T_w(2)/T)*(1+(gamma-1)/2*M^2)+0.5)^0.68*(1+(gamma-1)/2*M^2)^0.12 )^-1; % []
    
    h_chg = (0.026/(D_t^0.2))*(mu_g^0.2*Cp_g*Pr_g^-0.6)*(Pc/(C_star*C_star_eff))^0.8*(D_t/R)^0.1*(A_t/A)^0.9*sigma_c; % [W/m^2-K]
    h_rhg = (0.026/(D_t^0.2))*(mu_g^0.2*Cp_g*Pr_g^-0.6)*(Pc/(C_star*C_star_eff))^0.8*(D_t/R)^0.1*(A_t/A)^0.9*sigma_r; % [W/m^2-K]
    

    T_ci = T_w(1)-h_chg*t_r/kT_chg*(T_aw-T_w(1)); % K
    T_ri = T_w(2)-h_rhg*t_r/kT_rhg*(T_aw-T_w(2)); % K
    T_cb = (kT_ci/(t-t_r)*T_ci+h_c*T_c)/(kT_ci/(t-t_r)+h_c); % K
    T_rb = (kT_ri/(t-t_r)*T_ri+epsilon*eta*h_c*T_c)/(kT_ri/(t-t_r)+epsilon*eta*h_c); % K

    Q_cb = b*h_c*(T_cb-T_c); % W/m
    Q_rb = w*epsilon*eta*h_c*(T_rb-T_c); % W/m
    Q_cr = t*kT_ci/l*(T_ci-T_ri); % W/m
    Q_chg = b*(T_aw-T_w(1))/(1/h_chg+res); % W/m
    Q_rhg = w*(T_aw-T_w(2))/(1/h_rhg+res); % W/m
    
    F(1) = -Q_chg+2*Q_cr+Q_cb;
    F(2) = -Q_rhg-2*Q_cr+Q_rb;

end