%Rocket Project Ares 2024-2025
function [A,D,M,P,T,w_ch,h_ch,D_h,w_rib,t_ins,t_out,step,pos,D_t,A_t,Pc,Pe,mdot,MW_g,gamma,mu_g,Cp_g,k_g,Pr_g] = geometry(Thrust,Pc,Pe,C_star,C_star_eff,C_F,C_F_eff,L_star,angle_conv,h_ch,w_rib,w_ch_min,MW_g,gamma,mu_g,Cp_g,k_g,Pr_g,t_ins,t_out,T_tot,num_nodes)


    %% Geometry inputs


    %% Conversion to SI units
    
    Pc = Pc*6894.76; % pa
    Pe = Pe*6894.76; % pa
    t_ins = t_ins*0.0254; % m
    t_out = t_out*0.0254; % m
    h_ch = h_ch*0.0254; % m
    w_rib=w_rib*0.0254; % m
    w_ch_min=w_ch_min*0.0254; % m
    L_star=L_star*0.0254; % m
    Thrust=Thrust*4.44822; % N
    
    %% General Calcs
    %determine mdot using thrust, C*, CF
    mdot=Thrust/(C_star*C_star_eff*C_F*C_F_eff); %kg/s
    A_t = (C_star*C_star_eff*mdot)/Pc; %m^2
    D_t = sqrt(4*A_t/pi); %m
    %isentropic relations
    M_e = sqrt((2/(gamma(3)-1))*((Pc/Pe)^((gamma(3)-1)/gamma(3))-1)); %exit mach number using isentropic relations
    A_e = A_t*((gamma(3)+1)/2)^(-(gamma(3)+1)/(2*gamma(3)-2))*((1+((gamma(3)-1)/2)*M_e^2)^((gamma(3)+1)/(2*gamma(3)-2)))/M_e; %exit area using isentropic relations (m^2)
    D_e = sqrt(4*A_e/pi); %m
    
    %% Inital and final Rao nozzle
    
    %Screenshot for inital and final angles was taken from this paper
    %http://www.aspirespace.org.uk/downloads/Thrust%20optimised%20parabolic%20nozzle.pdf
    % and data was retrieved using https://automeris.io/
    %rao1 is for calculating initial angle, rao2 is final angle
    x_rao1 = [3.5491553337381196, 4.141624925980432, 4.879808892675213, 5.7495633367608825, 6.774338768275435, ...
         7.981765407112365, 9.404398155068069, 11.080594348242615, 13.055547955946532, 15.382508110411303, ...
         18.124214821568884, 21.354590586956895, 25.160733506306656, 29.6452656396067, 34.92909992558685, ...
         41.15470026288546, 48.489922652923525, 57.13254097022197, 67.31557938497336, 79.31359521181543, ...
         92.95695238881558];
    y_rao1 = [20.870097189228723, 21.832350192216822, 22.79030841045174, 23.63098477351592, 24.478980711103766, ...
         25.273225752568504, 25.981142327706802, 26.65016462867063, 27.295412256006063, 27.89160380240665, ...
         28.399078982351654, 28.935117695057162, 29.481278357426717, 30.04706189464332, 30.470412067126137, ...
         30.95457926718375, 31.541731415134752, 32.03808964536178, 32.51019749011783, 32.94023695653978, ...
         33.3902302207864];
    x_rao2 = [3.5904857429958956, 4.070428171054432, 5.719223314634882, 6.7385910816982015, ...
         7.939646219119776, 9.354771838877735, 11.022122868235531, 12.986654791257335, ...
         15.301335748427004, 18.028574675266487, 21.241903985740148, 25.02796216932407, ...
         29.488829756956964, 34.74478163869922, 40.937529941692276, 48.2340449036045, ...
         56.83105676079198, 66.96036003206942, 78.89506321335335, 92.71136447839058];
    y_rao2 = [14.610636017385616, 13.785853882713283, 12.227179902624336, 11.551761549093506, ...
         11.10107751906277, 10.530562062409295, 10.143993247575622, 9.761766942691338, ...
         9.46120661676666, 9.166282992159424, 8.876707492752566, 8.537921923493954, ...
         8.32444630434437, 8.086041139339187, 7.797795933292498, 7.602551481553931, ...
         7.383324445109793, 7.212617116207298, 7.0450458492968195, 6.9231916986699105];
    % Interpolation and extrapolation to any angle outside the region. 
    % Linear is a good approximation as the data covers most of the region, so
    % any low expansion ratio will still be accurate
    rao_i = interp1(x_rao1, y_rao1, A_e/A_t, 'linear', 'extrap');
    rao_f = interp1(x_rao2, y_rao2, A_e/A_t, 'linear', 'extrap');
    
    %% Optimize contraction ratio using empiracal data
    
    % Optimal contraction ratio and chamber diameter are found using the same
    % approach as the Rao angles, using the throat diameter/contraction ratio figure from Huzel
    % and Huang (fig 4-9)
    x_con = [0.4877466235754126, 0.6082422818553651, 0.8793641367127586, ...
         1.2409627479820822, 1.5472041635006222, 2.3002770856155674, ...
         2.8161301189578865, 3.4591676476892683, 4.618407539324998, ...
         6.824311241770623, 8.50469344365885];
    y_con = [7.701508441865165, 6.949877424291748, 5.865059345567369, ...
         5.03645148284939, 4.590384852410534, 3.849097108400973, ...
         3.534808108721339, 3.269344931244534, 2.9198504441387643, ...
         2.5177748988955715, 2.34145726849852];
    CR = interp1(x_con,y_con,D_t*39.3701, 'linear');
    D_c=sqrt(CR*D_t^2); %m
    
    %% Diverging
    
    % Equations for the diverging, converging and chamber geometry were used
    % from this paper unless otherwise specified
    %http://www.aspirespace.org.uk/downloads/Thrust%20optimised%20parabolic%20nozzle.pdf
    %length of diverging section is defined as 80% of a half angle 15 degree conical nozzle
    l_div = (.8*(sqrt(A_e/A_t)-1)*D_t/2)/tand(15);
    %Connection between radius at throat and diverging section is defined by a curve from
    %-90 degrees to the inital angle-90 degrees.
    %X and Y coordinates of the curve are defined in the Rao nozzle sizing
    %paper
    t3 = linspace(-90, rao_i - 90, .0225*num_nodes);
    x3 = 0.382 * D_t / 2 .* cosd(t3);
    y3 = 0.382 * D_t / 2 * sind(t3) + 0.382 * D_t / 2 + D_t/2;
    %N is defined as the point at the start of the parabolic arc. It is
    %the end point of the connecting curve
    n_x = x3(end);
    n_y = y3(end); 
    %E is the last point of the arc, and is at the nozzle exit
    e_x = l_div;
    e_y = D_e/2;
    %m1 and m2 are the slopes of the tangential lines of N and E
    m_1 = tand(rao_i);
    m_2 = tand(rao_f);
    %c_1 and c_2 are the intercepts 
    c_1 = n_y-m_1*n_x;
    c_2 = e_y-m_2*e_x;
    %q is the intersection of the tangential lines from the inital and final
    %angles
    q_x = (c_2-c_1)/(m_1-m_2);
    q_y = (m_1*c_2-m_2*c_1)/(m_1-m_2);
    %The bell contour is a Bezier curve defined by the x4 and y4 equations
    % found in the reference paper attached at the start of the diverging geometry definition
    t4 = linspace(0,1,.2725*num_nodes);
    x4 = (1-t4).^2*n_x+2*(1-t4).*t4*q_x+t4.^2*e_x;
    y4 = (1-t4).^2*n_y+2*(1-t4).*t4*q_y+t4.^2*e_y;

    x_div = flip([-x3 -x4(2:end)] + x4(end));
    y_div = flip([y3 y4(2:end)]);
    
    %% Converging
    
    %chamber to converging section is defined as a curve from 90 degrees to
    %90-the specified converging angle
    t0 = linspace(90,90-angle_conv,.05*num_nodes);
    x0 = (0.5258*D_t/2).*cosd(t0);
    y0 = D_c/2-0.5258*D_t/2+(0.5258*D_t/2).*sind(t0);
    %x2 and y2 represent the curve between the converging section and the
    %throat. Again, it is defined by a curve specified by the throat diameter times a scalar.
    t2 = linspace(-90-angle_conv,-90, .14*num_nodes);
    x2 = 1.5*D_t/2*cosd(t2);
    y2 = 1.5*D_t/2*sind(t2)+1.5*D_t/2+D_t/2;
    %x1 and y1 is the straight section characterizing the converging section,
    %which connects the two curves
    x1_right_point = (y0(end)-y2(1))/tand(angle_conv)+x0(end);
    x1 = linspace(x0(end),x1_right_point, .07*num_nodes);
    y1 = linspace(y0(end), y2(1), .07*num_nodes);
    %curve 2 is offset so sections can be combined
    x2 = x2 + (x1(end)+abs(x2(1)));
    %end of x0,x1, y0,y1 are removed due to overlapping x values when used to
    %calcualte volume
    x_conv = flip([-x0(1:end), -x1(2:end-1), -x2(1:end-1)]+x2(end)+x_div(end));
    y_conv = flip([y0(1:end), y1(2:end-1), y2(1:end-1)]);
    
    %% Chamber
    
    %converging section is splined together
    f_conv = @(var) interp1(x_conv, y_conv, var, 'spline');
    %volume of converging section is integrated to obtain volume
    %https://www.rit.edu/academicsuccesscenter/sites/rit.edu.academicsuccesscenter/files/documents/math-handouts/C8_VolumesbyIntegration_BP_9_22_14.pdf
    vol_conv = pi*integral(@(x_conv) f_conv(x_conv).^2, x_conv(1), x_conv(end));
    %volume calculated using L* definition
    vol_c = L_star * A_t;
    a_c = pi/4*D_c^2;
    %volume of combustion chamber adjusted as L* includes the converging section
    vol_cc = vol_c - vol_conv;
    %finally the length of the the chamber section can be computed and linspace
    %can be used to create nodes
    l_c = vol_cc /(a_c);
    xc=linspace(x_conv(end),x_conv(end)+l_c,floor(.45*num_nodes));
    xc = xc(2:end);
    yc=D_c/2*ones(1,floor(.45*num_nodes));
    yc = yc(2:end);
    
    %% Full Geometry
    
    %the positon of the nodes is an array of all sections combined together.
    %Sections are adjusted to line up, flip is used to reverse position
    %definition so that pos starts at zero, and *-1 is used to flip positon so
    %that the exit is positioned at x=0
    pos = [x_div x_conv xc]; %m
    %radius is just y position of each section
    r = [y_div y_conv yc]; %m
    D = 2*r; %m
    A = pi*r.^2; %m
    %step is distance between nodes
    step=diff(pos);
    
    %% Channels
    
    %circumference of the outer side of the inner wall is calculated
    circ = 2*pi*(r+t_ins);
    %throat circumference is calculated
    circ_t = 2*pi*D_t/2;
    %number of channels is the minimum that can fit the throat with a minimum
    %rib thickness and minimum channel width (channel geometry method may be updated in the future)
    num_ch = floor(circ_t/(w_rib+w_ch_min));
    %actual channel width can be found now based on the number of channels
    w_ch = circ/num_ch-w_rib;
    %hydraulic diameter definition
    D_h = 4*h_ch*w_ch./(2*h_ch+2*w_ch);
    
    %% Flow calcs
    
    %gamma is linearly interpolated between exit, throat, and chamber values.
    %position vector is split into diverging, converging, and chamber sections
    %to determine length of gamma vectors
    pos_div = length(x_div);
    pos_conv = length(x_conv);
    pos_c = length(xc);

    %gamma is combined
    MW_g = [linspace(MW_g(4),MW_g(3), pos_div),linspace(MW_g(3),MW_g(2),pos_conv),linspace(MW_g(2),MW_g(1),pos_c)];
    gamma = [linspace(gamma(4),gamma(3), pos_div),linspace(gamma(3),gamma(2),pos_conv),linspace(gamma(2),gamma(1),pos_c)];
    mu_g = [linspace(mu_g(4),mu_g(3), pos_div),linspace(mu_g(3),mu_g(2),pos_conv),linspace(mu_g(2),mu_g(1),pos_c)];
    Cp_g = [linspace(Cp_g(4),Cp_g(3), pos_div),linspace(Cp_g(3),Cp_g(2),pos_conv),linspace(Cp_g(2),Cp_g(1),pos_c)];
    k_g = [linspace(k_g(4),k_g(3), pos_div),linspace(k_g(3),k_g(2),pos_conv),linspace(k_g(2),k_g(1),pos_c)];
    Pr_g = [linspace(Pr_g(4),Pr_g(3), pos_div),linspace(Pr_g(3),Pr_g(2),pos_conv),linspace(Pr_g(2),Pr_g(1),pos_c)];

    
    %A is seperated into subsonic and supersonic ssections
    A_sup = A(1:pos_div); %m^2
    A_sub = A(pos_div+1:end); %m^2
    %flow isentropic returns the mach number given an area and a value of gamma
    for i=1:length(A_sub)
        M_sub(i) = flowisentropic(gamma(i),A_sub(i)/A_t,'sub');
    end
    %calculation is repeated for the supersonic section of the nozzle
    for i=1:length(A_sup)
        if A_sup(i)/A_t < 1
            A_sup(i) = A_t;
        end
        M_sup(i) = flowisentropic(gamma(i),A_sup(i)/A_t,'sup');
    end
    %Mach number arrays are concatenated together
    M = [M_sup M_sub];
    %static pressure is computed using the local mach number
    P = Pc * (1+(gamma-1)/2.*M.^2).^(-gamma/(gamma-1)); %pa
    T = (1+(gamma-1)/2.*M.^2).^-1*T_tot;



%% Test

% figure()
% scatter(pos,r)
% hold on
% plot(pos,-1*r)
% axis equal
% title('Engine Contour')
% hold off
% figure()
% plot(pos,M);
% title('Mach number')
% figure()
% plot(pos,P)
% title('Static Pressure (pa)')
% figure()
% plot(pos,T)
% title('Static Temperature (K)')

end