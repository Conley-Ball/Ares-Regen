function [T_c,T_sat_c,P_c,q,T_chg,T_rhg,T_ci,T_ri,T_cb,T_rb,T_rt,T_ro,T_co,T_ct] = heatTransfer2D(Pc,M,A,D,D_h,w_ch,h_ch,w_rib,num,t_ins,t_out,step,pos,gamma,mu_g,Cp_g,Pr_g,C_star,T,D_t,A_t,T_inlet,ratio,mdot_f,C_star_eff,res)
    T_c = zeros(1,length(pos));
    P_c = zeros(1,length(pos));

    T_sat_c = zeros(1,length(pos));
    h_c = zeros(1,length(pos));
    h_chg = zeros(1,length(pos));
    h_rhg = zeros(1,length(pos));

    T_ci = zeros(1,length(pos));
    T_ri = zeros(1,length(pos));
    T_cb = zeros(1,length(pos));
    T_rb = zeros(1,length(pos));
    T_rt = zeros(1,length(pos));
    T_ro = zeros(1,length(pos));
    T_co = zeros(1,length(pos));
    T_ct = zeros(1,length(pos));

    T_chg = zeros(1,length(pos));
    T_rhg = zeros(1,length(pos));
    q = zeros(1,length(pos));

    step = [step 1];

    T_c(1) = T_inlet;
    pgs = 0;
    P_c_inj = -Pc;
    tol = 6894.76;
    
    while abs(P_c_inj-Pc) > tol
        P_c_prev = P_c(1);
        P_c = zeros(1,length(pos));
        P_c(1) = P_c_prev+Pc-P_c_inj;
        for i=1:length(pos)
            %% Coolant Properties
            [rho_c,mu_c,k_c,Cp_c,T_sat_c(i),v_c] = coolant(T_c(i),P_c(i),ratio,mdot_f,w_ch(i),h_ch(i),num);
            
            %% Coolant Side Heat Transfer
    
            Re_c = D_h(i)*v_c*rho_c/mu_c;
            Pr_c = mu_c*Cp_c/k_c;
            Nu_c = 0.023*Re_c^0.8*Pr_c^0.4;

            h_c(i) = Nu_c*k_c/D_h(i); % [W/m^2-K];
            if T_sat_c(i) < T_c(i)
                warning('Coolant is boiling. This may cause an error. Proceed?')
                pause
            end
            
            %% Gas Side Heat Transfer
            
            % Geometry
            w = w_rib; % m
            b = w_ch(i); % m
            h = h_ch(i); % m
            t = t_ins; % m
            t_r = t/2; % m
            d = t_out; % m
            d_r = d/2; % m
            l = (b+w)/2; % m
            R = 0.5258*D_t/2; % m

            % Conductivity

            if i > 1
                kT_chg = kSteel174ph(T_chg(i-1));
                kT_rhg = kSteel174ph(T_rhg(i-1));
                kT_ci = kSteel174ph(T_ci(i-1));
                kT_ri = kSteel174ph(T_ri(i-1));
                kT_cb = kSteel174ph(T_cb(i-1));
                kT_rb = kSteel174ph(T_rb(i-1));
                kT_rt = kSteel174ph(T_rt(i-1));
            else
                kT_chg = kSteel174ph(300);
                kT_rhg = kSteel174ph(300);
                kT_ci = kSteel174ph(300);
                kT_ri = kSteel174ph(300);
                kT_cb = kSteel174ph(300);
                kT_rb = kSteel174ph(300);
                kT_rt = kSteel174ph(300);
            end
            error = 300;
            while abs(error) > 1
            
                r = Pr_g(i)^0.33; % []
                epsilon = ((2*kT_rb)/(w*h_c(i)))^0.5; % []
                m = ((2*h_c(i))/(kT_rb*w))^0.5; % 1/m
                H = (d_r/kT_rb+(w*l)/(2*d*kT_rb)+(w*d_r)/(b*kT_rb)+w/(h_c(i)*b))^-1; % [W/m^2-K]
                beta = (H/(kT_rb*m)*sinh(m*h)+cosh(m*h))^-1; % []
                eta = (cosh(m*h)-beta)/sinh(m*h); % []
                
                T_aw = T(i)*((1+r*(gamma(i)-1)/2*M(i)^2)/(1+(gamma(i)-1)/2*M(i)^2)); % K
    
                eq = @(T_w) solveFunction(T_w,w,b,t,t_r,l,epsilon,eta,T_aw,M(i),T(i),A(i),gamma(i),Cp_g,mu_g,Pr_g,Pc,C_star,C_star_eff,D_t,R,A_t,res,h_c(i),T_c(i),kT_chg,kT_rhg,kT_ci,kT_ri);
                
                if i > 1
                    sol = fsolve(eq,[T_chg(i-1),T_rhg(i-1)],optimset('Display', 'off', 'TolX',1e-4));
                else
                    sol = fsolve(eq,[0,0],optimset('Display', 'off', 'TolX',1e-4));
                end
                T_chg_old = T_chg(i);
                T_chg(i) = sol(1);
                T_rhg(i) = sol(2);
    
                error = T_chg(i) - T_chg_old;
    
                sigma_c = ( (0.5*(T_chg(i)/T(i))*(1+(gamma(i)-1)/2*M(i)^2)+0.5)^0.68*(1+(gamma(i)-1)/2*M(i)^2)^0.12 )^-1; % []
                sigma_r = ( (0.5*(T_rhg(i)/T(i))*(1+(gamma(i)-1)/2*M(i)^2)+0.5)^0.68*(1+(gamma(i)-1)/2*M(i)^2)^0.12 )^-1; % []
                
                h_chg(i) = (0.026/(D_t^0.2))*(mu_g(end)^0.2*Cp_g(end)*Pr_g(end)^-0.6)*(Pc/(C_star*C_star_eff))^0.8*(D_t/R)^0.1*(A_t/A(i))^0.9*sigma_c; % [W/m^2-K]
                h_rhg(i) = (0.026/(D_t^0.2))*(mu_g(end)^0.2*Cp_g(end)*Pr_g(end)^-0.6)*(Pc/(C_star*C_star_eff))^0.8*(D_t/R)^0.1*(A_t/A(i))^0.9*sigma_r; % [W/m^2-K]
                
                T_ci(i) = T_chg(i)-h_chg(i)*t_r/kT_chg*(T_aw-T_chg(i)); % K
                T_ri(i) = T_rhg(i)-h_rhg(i)*t_r/kT_rhg*(T_aw-T_rhg(i)); % K
                T_cb(i) = (kT_ci/(t-t_r)*T_ci(i)+h_c(i)*T_c(i))/(kT_ci/(t-t_r)+h_c(i)); % K
                T_rb(i) = (kT_ri/(t-t_r)*T_ri(i)+epsilon*eta*h_c(i)*T_c(i))/(kT_ri/(t-t_r)+epsilon*eta*h_c(i)); % K
    
                T_rt(i) = beta*(T_rb(i)-T_c(i)) + T_c(i); % K
                T_ro(i) = T_rt(i) - H*d_r/kT_rt*(T_rt(i)-T_c(i)); % K
                T_co(i) = T_ro(i) - H*w*l/(4*kT_rt*d)*(T_rt(i)-T_c(i)); % K
                T_ct(i) = T_co(i) - H*w*d_r/(kT_rt*b)*(T_rt(i)-T_c(i)); % K
    
                kT_chg = kSteel174ph(T_chg(i));
                kT_rhg = kSteel174ph(T_rhg(i));
                kT_ci = kSteel174ph(T_ci(i));
                kT_ri = kSteel174ph(T_ri(i));
                kT_cb = kSteel174ph(T_cb(i));
                kT_rb = kSteel174ph(T_rb(i));
                kT_rt = kSteel174ph(T_rt(i));

            end

            % Q_chg = b*(T_aw-T_chg(i))/(1/h_chg(i)+res); % W/m
            % Q_rhg = w*(T_aw-T_rhg(i))/(1/h_rhg(i)+res); % W/m
            Q_chg = (T_aw-T_chg(i))/(1/(b*h_chg(i))+res/b); % W/m
            Q_rhg = (T_aw-T_rhg(i))/(1/(w*h_rhg(i))+res/w); % W/m
            
            Q_tot = Q_chg + Q_rhg; % W/m
            q(i) = Q_tot/(b+w); % W/m^2
            if q < 0
                error('Negative Heat Flux')
            end
    
            epsilon_w = 2e-6; % m
            f = (-1.8*log10((epsilon_w/3.7/D_h(i))^1.11+6.9/Re_c))^-2; % []
            if i+1 <= length(pos)
                if D_h(i+1) > D_h(i)
                    K_L = ((D_h(i)/D_h(i+1))^2-1)^2;
                elseif D_h(i+1) < D_h(i)
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
            progress = pgs/length(pos);
            clc
            fprintf('Progress: %.2f%%\n',progress*100)
    
        end
        P_c_inj = P_c(end);
        pgs = 0;
    end

end