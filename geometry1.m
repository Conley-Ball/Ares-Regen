function [A,D,M,P,l_ch,w_ch,D_h,step,pos,D_t,A_t,mdoti] = geometry(Thrust,Pc,Pe,C_star,C_star_eff,C_F,C_F_eff,L_star,w_ch,num_ch,l_rib,D_c,gamma,th_in,num_nodes)
%% Conversions
Pc = Pc*6894.76;
Pe = Pe*6894.76;
th_in = th_in*0.0254;
w_ch = w_ch*0.0254;
l_rib=l_rib*0.000508;
D_c=D_c*0.0254;
L_star=L_star*0.0254;
Thrust=Thrust*4.44822;
%% General Calcs
mdot=Thrust/(C_star*C_star_eff*C_F*C_F_eff);
mdoti=mdot;
A_t = (C_star*C_star_eff*mdot)/Pc;
D_t = sqrt(4*A_t/pi); %in
m_e = sqrt((2/(gamma-1))*((Pc/Pe)^((gamma-1)/gamma)-1));
a_e = A_t*((gamma+1)/2)^(-(gamma+1)/(2*gamma-2))*((1+((gamma-1)/2)*m_e^2)^((gamma+1)/(2*gamma-2)))/m_e;
d_e = sqrt(4*a_e/pi);
%% Inital and final Rao nozzles
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
% Perform the interpolation with extrapolation
rao_i = interp1(x_rao1, y_rao1, a_e/A_t, 'linear', 'extrap');
rao_f = interp1(x_rao2, y_rao2, a_e/A_t, 'linear', 'extrap');
%% Diverging
l_div = (.8*(sqrt(a_e/A_t)-1)*D_t/2)/tand(15);
t3 = linspace(-90, rao_i - 90, .04*num_nodes);
x3 = 0.382 * D_t / 2 .* cosd(t3);
y3 = 0.382 * D_t / 2 * sind(t3) + 0.382 * D_t / 2 + D_t/2;
n_x = x3(end);
n_y = y3(end); 
e_x = l_div;
e_y = d_e/2;
m_1 = tand(rao_i);
m_2 = tand(rao_f);
c_1 = n_y-m_1*n_x;
c_2 = e_y-m_2*e_x;
q_x = (c_2-c_1)/(m_1-m_2);
q_y = (m_1*c_2-m_2*c_1)/(m_1-m_2);
t4 = linspace(0,1,.2*num_nodes);
x4 = (1-t4).^2*n_x+2*(1-t4).*t4*q_x+t4.^2*e_x;
y4 = (1-t4).^2*n_y+2*(1-t4).*t4*q_y+t4.^2*e_y;
%% Converging
t0 = linspace(90,45,.05*num_nodes);
x0 = (0.5258*D_t/2).*cosd(t0);
y0 = D_c/2-0.5258*D_t/2+(0.5258*D_t/2).*sind(t0);
t2 = linspace(-135,-90, .18*num_nodes);
x2 = 1.5*D_t/2*cosd(t2);
y2 = 1.5*D_t/2*sind(t2)+1.5*D_t/2+D_t/2;
x1_right_point = y0(end) + x0(end) - y2(1);
x1 = linspace(x0(end),x1_right_point, .03*num_nodes);
y1 = linspace(y0(end), y2(1), .03*num_nodes);
x2 = x2 + (x1(end)+abs(x2(1)));
%% Chamber
x_conv = [x0(1:end-1), x1(1:end-1), x2];
y_conv = [y0(1:end-1), y1(1:end-1), y2];
f_conv = @(var) interp1(x_conv, y_conv, var, 'spline');
vol_conv = pi*integral(@(x_conv) f_conv(x_conv).^2, x_conv(1), x_conv(end));
vol_c = L_star * A_t;
a_c = pi/4*D_c^2;
vol_cc = vol_c - vol_conv;
l_c = vol_cc /(a_c);
xc=linspace(-l_c,0,.5*num_nodes);
yc=D_c/2*ones(1,.5*num_nodes);
%% Full Geometry
pos = flip((cat(2,xc- x2(end),x0- x2(end),x1- x2(end), x2- x2(end), x3, x4))-x4(end))*-1;
r = flip(cat(2, yc, y0, y1, y2, y3, y4));    
D = 2*r;
A = pi*r.^2;

step=diff(pos);
A_sup = A(1:length(y3)+length(y4));
A_sub = A(length(y3)+length(y4):end);
for i=1:length(A_sub)
    if A_sub(i)/A_t < 1
        A_sub(i) = A_t;
    end
    M_sub(i) = flowisentropic(gamma,A_sub(i)/A_t,'sub');
end
for i=1:length(A_sup)
    if A_sup(i)/A_t < 1
        A_sup(i) = A_t;
    end
    M_sup(i) = flowisentropic(gamma,A_sup(i)/A_t,'sup');
end
M = [M_sup M_sub];
M = M(1:end-1);

P = Pc * (1+(gamma-1)/2*M.^2).^(-gamma/(gamma-1));


%% Channels
th_in=th_in*ones(1,length(pos));
circ = 2*pi*(r+th_in);
l_ch = circ/num_ch-l_rib;
D_h = 4*w_ch*l_ch./(2*l_ch+2*w_ch);
end
