function [T_c,T_sat_c,P_c,q,T_chg,T_rhg,T_ci,T_ri,T_cb,T_rb,T_rt,T_ro,T_co,T_ct] = heatTransfer2D(Pc,M,P,A,D,D_h,w_ch,h_ch,w_rib,num,t_ins,t_out,step,pos,gamma,mu_g,Cp_g,k_g,Pr_g,C_star,T,k_w,D_t,A_t,T_inlet,ratio,mdot_f,C_star_eff,res)
    T_c = zeros(1,length(pos));
    P_c = zeros(1,length(pos));
    q = zeros(1,length(pos));
    T = T*C_star_eff^2;

    h_ch = h_ch*ones(1,length(w_ch));
    step = [step 1];

    T_c(1) = T_inlet;
    pgs = 0;
    T = T*C_star_eff^2;
    P_c_inj = 0;
    P_c_prev = Pc*1.5;
    tol = 6000;
    
    while abs(P_c_inj-Pc) > tol
        P_c_prev = P_c(1);
        plot(pos,P_c)
        P_c = zeros(1,length(pos));
        P_c(1) = P_c_prev+Pc-P_c_inj;
        for i=1:length(pos)
            %% Coolant Properties
            [rho_c,mu_c,k_c,Cp_c,T_sat_c(i),v_c,h_lg,rho_c_liquid,rho_c_vapor] = coolant(T_c(i),P_c(i),ratio,mdot_f,w_ch(i),h_ch(i),num);
            %% Coolant Side Heat Transfer
    
            % [~,~,~,k,~] = materialPropertiesSteel174PH(T);
            Re_c = D_h(i)*v_c*rho_c/mu_c; % [m] * [m/s] * [kg/m^3] / [Pa-s] = []
            Pr_c = mu_c*Cp_c/k_c; % [kg/m-s] * [J/kg-K] / [W/m-K] = []
            Nu_c = 0.023*Re_c^0.8*Pr_c^0.4; % []


            h_c = Nu_c*k_c/D_h(i); % [W/m^2-K]
            if T_sat_c(i) < T_c(i)
                warning('Coolant is boiling. This may cause an error. Proceed?')
                pause
            end
            
            %% Gas Side Heat Transfer
            
            w = w_rib; % m
            b = w_ch(i); % m
            h = h_ch(i); % m
            t = t_ins; % m
            t_r = t/2; % m
            d = t_out; % m
            d_r = d/2; % m
            l_i = (b+w)/2; % m
            l_o = l_i; % m
            l = l_i; % m
    
    
            
            syms T_chgs
            syms T_rhgs
    
            % Pr_g = mu_g(i)*Cp_g(i)/k_g(i);
            % Pr_g = 4*gamma(i)/(9*gamma(i)-5); % []
            r = Pr_g(i)^0.33; % []
            R = 0.5258*D_t/2; % m
    
            epsilon = ((2*k_w)/(w*h_c))^0.5; % []
            m = ((2*h_c)/(k_w*w))^0.5; % 1/m
            H = (d_r/k_w+(w*l)/(2*d*k_w)+(w*d_r)/(b*k_w)+w/(h_c*b))^-1; % [W/m^2-K]
            beta = (H/(k_w*m)*sinh(m*h)+cosh(m*h))^-1; % []
            eta = (cosh(m*h)-beta)/sinh(m*h); % []
    
            
            T_aw = T(i)*((1+r*(gamma(i)-1)/2*M(i)^2)/(1+(gamma(i)-1)/2*M(i)^2)); % K
            
            sigma_c = ( (0.5*(T_chgs/T(i))*(1+(gamma(i)-1)/2*M(i)^2)+0.5)^0.68*(1+(gamma(i)-1)/2*M(i)^2)^0.2 )^-1; % []
            sigma_r = ( (0.5*(T_rhgs/T(i))*(1+(gamma(i)-1)/2*M(i)^2)+0.5)^0.68*(1+(gamma(i)-1)/2*M(i)^2)^0.2 )^-1; % []
            
            h_chg = (0.026/(D_t^0.2))*(mu_g(i)^0.2*Cp_g(i)*Pr_g(i)^-0.6)*(Pc/C_star)^0.8*(D_t/R)^0.1*(A_t/A(i))^0.9*sigma_c; % [W/m^2-K]
            h_rhg = (0.026/(D_t^0.2))*(mu_g(i)^0.2*Cp_g(i)*Pr_g(i)^-0.6)*(Pc/C_star)^0.8*(D_t/R)^0.1*(A_t/A(i))^0.9*sigma_r; % [W/m^2-K]
             
            T_ci(i) = T_chgs-h_chg*t_r/k_w*(T_aw-T_chgs); % K
            T_ri(i) = T_rhgs-h_rhg*t_r/k_w*(T_aw-T_rhgs); % K
            T_cb(i) = (k_w/(t-t_r)*T_ci(i)+h_c*T_c(i))/(k_w/(t-t_r)+h_c); % K
            T_rb(i) = (k_w/(t-t_r)*T_ri(i)+epsilon*eta*h_c*T_c(i))/(k_w/(t-t_r)+epsilon*eta*h_c); % K
    
            Q_cb = b*h_c*(T_cb(i)-T_c(i)); % W/m
            Q_rb = w*epsilon*eta*h_c*(T_rb(i)-T_c(i)); % W/m
            Q_cr = t*k_w/l*(T_ci(i)-T_ri(i)); % W/m
            Q_chg = b*(T_aw-T_chgs)/(1/h_chg+res); % W/m
            Q_rhg = w*(T_aw-T_rhgs)/(1/h_rhg+res); % W/m
            
            eq1 = sym(matlabFunction(-Q_chg+2*Q_cr+Q_cb));
            eq2 = sym(matlabFunction(-Q_rhg-2*Q_cr+Q_rb));
    
            [T_chg(i),T_rhg(i)] = vpasolve(eq1,eq2,T_chgs,T_rhgs);
    
            sigma_c = ( (0.5*(T_chg(i)/T(i))*(1+(gamma(i)-1)/2*M(i)^2)+0.5)^0.68*(1+(gamma(i)-1)/2*M(i)^2)^0.2 )^-1; % []
            sigma_r = ( (0.5*(T_rhg(i)/T(i))*(1+(gamma(i)-1)/2*M(i)^2)+0.5)^0.68*(1+(gamma(i)-1)/2*M(i)^2)^0.2 )^-1; % []
            
            h_chg = (0.026/(D_t^0.2))*(mu_g(i)^0.2*Cp_g(i)*Pr_g(i)^-0.6)*(Pc/C_star)^0.8*(D_t/R)^0.1*(A_t/A(i))^0.9*sigma_c; % [W/m^2-K]
            h_rhg = (0.026/(D_t^0.2))*(mu_g(i)^0.2*Cp_g(i)*Pr_g(i)^-0.6)*(Pc/C_star)^0.8*(D_t/R)^0.1*(A_t/A(i))^0.9*sigma_r; % [W/m^2-K]
            
            T_ci(i) = T_chg(i)-h_chg*t_r/k_w*(T_aw-T_chg(i)); % K
            T_ri(i) = T_rhg(i)-h_rhg*t_r/k_w*(T_aw-T_rhg(i)); % K
            T_cb(i) = (k_w/(t-t_r)*T_ci(i)+h_c*T_c(i))/(k_w/(t-t_r)+h_c); % K
            T_rb(i) = (k_w/(t-t_r)*T_ri(i)+epsilon*eta*h_c*T_c(i))/(k_w/(t-t_r)+epsilon*eta*h_c); % K

            T_rt(i) = beta*(T_rb(i)-T_c(i)) + T_c(i); % K
            T_ro(i) = T_rt(i) - H*d_r/k_w*(T_rt(i)-T_c(i)); % K
            T_co(i) = T_ro(i) - H*w*l/(4*k_w*d)*(T_rt(i)-T_c(i)); % K
            T_ct(i) = T_co(i) - H*w*d_r/(k_w*b)*(T_rt(i)-T_c(i)); % K

            Q_cb = b*h_c*(T_cb(i)-T_c(i)); % W/m
            Q_rb = w*epsilon*eta*h_c*(T_rb(i)-T_c(i)); % W/m
            Q_cr = t*k_w/l*(T_ci(i)-T_ri(i)); % W/m
            Q_chg = b*(T_aw-T_chg(i))/(1/h_chg+res); % W/m
            Q_rhg = w*(T_aw-T_rhg(i))/(1/h_rhg+res); % W/m
            
            Q_tot = Q_chg + Q_rhg; % W/m
            q(i) = Q_tot/(b+w); % W/m^2
            if q < 0
                error('Negative Heat Flux')
            end
    
            epsilon_w = 2e-6; % m
            f = (-1.8*log10((epsilon_w/3.7/D_h(i))^1.11+6.9/Re_c))^-2; % []
            if i+1 <= length(pos)
                if D_h(i+1) < D_h(i)
                    K_L = ((D_h(i)/D_h(i+1))^2-1)^2;
                elseif D_h(i+1) > D_h(i)
                    K_L = 0.5-0.167*(D_h(i+1)/D_h(i))-0.125*(D_h(i+1)/D_h(i))^2-0.208*(D_h(i+1)/D_h(i))^3;
                else
                    K_L = 0;
                end
                P_c(i+1) = P_c(i) - rho_c*v_c^2/2*(K_L+f*step(i)/D_h(i)); % Pa
                T_c(i+1) = T_c(i) + pi*D(i)*step(i)*q(i)/(mdot_f*Cp_c); % K
                if P_c(i+1) < 0
                    error('Negative Coolant Pressure')
                end
            end
            pgs = pgs + 1;
            progress = pgs/length(pos)
    
        end
        P_c_inj = P_c(end);
    end

end