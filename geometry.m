%% Rocket Project 2024-2025
function [A,D,M,P,l_ch,w_ch,D_h,step,pos,D_t,A_t] = geometry(Thrust,Pc,Pe,C_star,C_star_eff,C_F,C_F_eff,L_star,w_ch,num_ch,l_rib,D_c,gamma)


l_rao_bell = 0.8;
resolution = 100;
rao_i = 32.5;
rao_f = 12;
%% General Calcs
mdot=Thrust/(C_star*C_star_eff*C_F*C_F_eff);%kg/s
C_star=C_star*39.3701; %m/2 to in/s
mdot=mdot*0.0057101471301634; % kg/s to slinches/s
A_t = (C_star*C_star_eff*mdot)/Pc;
D_t = sqrt(4*A_t/pi); %in
M_e = sqrt((2/(gamma-1))*((Pc/Pe)^((gamma-1)/gamma)-1));
A_e = A_t*((gamma+1)/2)^(-(gamma+1)/(2*gamma-2))*((1+((gamma-1)/2)*M_e^2)^((gamma+1)/(2*gamma-2)))/M_e;
D_e = sqrt(4*A_e/pi);
%% Diverging
l_div = (l_rao_bell*(sqrt(A_e/A_t)-1)*D_t/2)/tand(15);
t3 = linspace(-90, rao_i - 90, resolution);
x3 = 0.382 * D_t / 2 .* cosd(t3);
y3 = 0.382 * D_t / 2 * sind(t3) + 0.382 * D_t / 2 + D_t/2;
n_x = x3(end);
n_y = y3(end); 
e_x = l_div;
e_y = D_e/2;
m_1 = tand(rao_i);
m_2 = tand(rao_f);
c_1 = n_y-m_1*n_x;
c_2 = e_y-m_2*e_x;
q_x = (c_2-c_1)/(m_1-m_2);
q_y = (m_1*c_2-m_2*c_1)/(m_1-m_2);
t4 = linspace(0,1,resolution);
x4 = (1-t4).^2*n_x+2*(1-t4).*t4*q_x+t4.^2*e_x;
y4 = (1-t4).^2*n_y+2*(1-t4).*t4*q_y+t4.^2*e_y;
%% Converging
t0 = linspace(90,45,resolution);
x0 = (0.5258*D_t/2).*cosd(t0);
y0 = D_c/2-0.5258*D_t/2+(0.5258*D_t/2).*sind(t0);
t2 = linspace(-135,-90, resolution);
x2 = 1.5*D_t/2*cosd(t2);
y2 = 1.5*D_t/2*sind(t2)+1.5*D_t/2+D_t/2;
x1_right_point = y0(end) + x0(end) - y2(1);
x1 = linspace(x0(end),x1_right_point, (resolution/10)+2);
y1 = linspace(y0(end), y2(1), (resolution/10)+2);
x2 = x2 + (x1(end)+abs(x2(1)));
%% Chamber
f0 = @(var) interp1(x0, y0, var, 'spline');
vol0 = pi*integral(@(x0) f0(x0).^2, x0(1), x0(end));
f1 = @(var)interp1(x1, y1, var, 'spline');
vol1 = pi*integral(@(x1)f1(x1).^2, x1(1), x1(end));
f2 = @(var2)interp1(x2, y2, var2, 'spline');
vol2 = pi*integral(@(x2)f2(x2).^2, x2(1), x2(end));
vol_c = L_star * A_t;
A_c = pi/4*D_c^2;
vol_cc = vol_c - (vol0 + vol1 + vol2);
L_c = vol_cc /(A_c);
xc=linspace(-L_c,0,resolution);
yc=D_c/2*ones(1,resolution);
%% Full Geometry
pos = flip((cat(2,xc- x2(end),x0- x2(end),x1- x2(end), x2- x2(end), x3, x4))-x4(end))*-1;

r = flip(cat(2, yc, y0, y1, y2, y3, y4));  


D = 2*r;
A = pi*r.^2;

A_sup = A(1:length(y3)+length(y4));
A_sub = A(length(y3)+length(y4):end);
for i=1:length(A_sub)
    M_sub(i) = flowisentropic(gamma,A_sub(i)/A_t,'sub');
end
for i=1:length(A_sup)
    M_sup(i) = flowisentropic(gamma,A_sup(i)/A_t,'sup');
end
M = [M_sup M_sub];
M = M(1,end-1);

P = Pc * (1+(gamma-1)/2*M.^2).^(-gamma/(gamma-1));

circ = 2*pi*r;
step = diff(pos);
%% Channels
l_ch = circ/num_ch-l_rib;
D_h = 4*w_ch*l_ch./(2*l_ch+2*w_ch);
%% Gas
end