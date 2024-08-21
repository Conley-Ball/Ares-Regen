function [T_c,P_cool,q] = heatTransfer2D(Pc,M,P,A,D,D_h,l_ch,w_ch,l_rib,num,t_ins,t_out,step,pos,gamma,MW_g,mu_g,Cp_g,C_star,T_cns,k_w,D_t,A_t,T_inlet,ratio,mdot_f)
    T_c = zeros(1,length(pos));
    P_c = zeros(1,length(pos));
    q = zeros(1,length(pos));

    T_c(1) = T_inlet;
    P_c(1) = Pc;
    
    for i=1:length(pos)
        %% Coolant Properties
        [rho_c,mu_c,k_c,Cp_c,st_c,T_sat_c,v_c,h_lg,rho_c_liquid,rho_c_vapor] = coolant(T_c(i),P_c(i),ratio,mdot_f,l_ch(i),w_ch(i),num);
        %% Coolant Side Heat Transfer

        Re_c = D_h(i)*v_c*rho_c/mu_c;
        Pr_c = mu_c*Cp_c/k_c;
        Nu_c = 0.023*Re_c^0.8*Pr_c^0.4;

        h_c = Nu_c*k_c/D_h(i);
        if T_sat_c < T_c(i)
            warning('boiling')
        end
        
        %% Gas Side Heat Transfer
        
        w = l_rib(i);
        b = l_ch(i);
        h = w_ch(i);
        t = t_ins;
        t_r = t/2;
        d = t_out;
        d_r = d/2;
        l_i = (b+w)/2;
        l_o = l_i;
        l = l_i;



        syms T_chg
        syms T_rhg

        % Pr_g = mu_g(i)*Cp_g(i)/k_g(i);
        Pr_g = 4*gamma(i)/(9*gamma(i)-5);
        r = Pr_g^0.33;
        R = 0.5258*D_t/2;

        epsilon = ((2*k_w)/(w*h_c))^0.5;
        m = ((2*h_c)/(k_w*w))^0.5;
        H = (d_r/k_w+(w*l)/(2*d*k_w)+(w*d_r)/(b*k_w)+w/(h_c*b))^-1;
        beta = (H/(k_w*m)*sinh(m*h)+cosh(m*h))^-1;
        eta = (cosh(m*h)-beta)/sinh(m*h);

        T_aw = T_cns*((1+r*(gamma(i)-1)/2*M(i)^2)/(1+(gamma(i)-1)/2*M(i)^2))*0.9;
        sigma_c = ( (0.5*(T_chg/T_cns(i))*(1+(gamma(i)-1)/2*M(i)^2)+0.5)^0.68*(1+(gamma(i)-1)/2*M(i)^2)^0.2 )^-1;
        sigma_r = ( (0.5*(T_rhg/T_cns(i))*(1+(gamma(i)-1)/2*M(i)^2)+0.5)^0.68*(1+(gamma(i)-1)/2*M(i)^2)^0.2 )^-1;
        % sigma_c = simplify(sigma_c);
        % sigma_r = simplify(sigma_r);
        
        h_chg = (0.026/(D_t^0.2))*(mu_g(i)^0.2*Cp_g(i)*Pr_g^-0.6)*(Pc*9.81/C_star)^0.8*(D_t/R)^0.1*(A_t/A(i))^0.9*sigma_c;
        h_rhg = (0.026/(D_t^0.2))*(mu_g(i)^0.2*Cp_g(i)*Pr_g^-0.6)*(Pc*9.81/C_star)^0.8*(D_t/R)^0.1*(A_t/A(i))^0.9*sigma_r;
        % h_chg = simplify(h_chg);
        % h_rhg = simplify(h_rhg);
        
        T_ci = T_chg-h_chg*t_r/k_w*(T_aw-T_chg);
        T_ri = T_rhg-h_rhg*t_r/k_w*(T_aw-T_rhg);
        T_cb = (k_w/(t-t_r)*T_ci+h_c*T_c(i))/(k_w/(t-t_r)+h_c);
        T_rb = (k_w/(t-t_r)*T_ri+epsilon*eta*h_c*T_c(i))/(k_w/(t-t_r)+epsilon*eta*h_c);
        
        
        % T_ci = simplify(T_ci);
        % T_ri = simplify(T_ri);
        % T_cb = simplify(T_cb);
        % T_rb = simplify(T_rb);

        Q_cb = b*h_c*(T_cb-T_c(i));
        Q_rb = w*epsilon*eta*h_c*(T_rb-T_c(i));
        Q_cr = t*k_w/l*(T_ci-T_ri);
        Q_chg = b*h_chg*(T_aw-T_chg);
        Q_rhg = w*h_rhg*(T_aw-T_rhg);

        % Q_cb = simplify(Q_cb);
        % Q_rb = simplify(Q_rb);
        % Q_cr = simplify(Q_cr);
        % Q_chg = simplify(Q_chg);
        % Q_rhg = simplify(Q_rhg);
        
        eq1 = sym(matlabFunction(0==-Q_chg+2*Q_cr+Q_cb));
        eq2 = sym(matlabFunction(0==-Q_rhg-2*Q_cr+Q_rb));

        vpasolve(eq1,eq2,T_chg,T_rhg)


        
        % hold on
        % rectangle('Position',[0,0,b+w,t+h+d])
        % rectangle('Position',[0,t,b,t+h])
        % hold off




        syms x;
        syms y;
        a1 = 2*x+3*y;
        a2 = 7*y-1;
        fun1 = matlabFunction(a1)
        fun2 = matlabFunction(a2)
        solve(fun1,fun2,x,y)




    end

end