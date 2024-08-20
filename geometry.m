%% Rocket Project 2024-2025
function [mx,p_stat,d_h,a,d,l_ch,w_ch,step,pos] = geometry(pc,gamma,l_c,th_in,th_out,d_c_in,l_c_ch,w_c_ch,nodes_c,l_conv,angle_conv,nodes_conv,l_div,d_t,d_e,l_div_ch,w_div_ch,nodes_div,th_l_ch,th_w_ch)

%%
th_l_ch = th_l_ch*0.0254; %m
th_w_ch = th_w_ch*0.0254; %m
%% Diverging calcs
l_div_ch = l_div_ch*0.0254; %m
w_div_ch = w_div_ch*0.0254; %m
a_t = pi*d_t^2/4; %in^2
a_e = d_e^2*pi/4; %in^2
step_div = l_div/(nodes_div-1); %in
pos_div = 0:step_div:l_div; %in
slope_div = (a_t-a_e)./(l_div);
slope_div_l = (th_l_ch-l_div_ch)/(l_div); %m/in
slope_div_w = (th_w_ch-w_div_ch)/(l_div); %m/in
a_div = slope_div.*pos_div+a_e; %in^2
l_div_ch_array = slope_div_l.*pos_div+l_div_ch; %m
w_div_ch_array = slope_div_w.*pos_div+w_div_ch; %m
d_h_div = 4.*l_div_ch_array.*w_div_ch_array./(2*l_div_ch_array+2*w_div_ch_array); %m
%Flow calcs
for i=1:length(pc)
    for j=1:length(pos_div)
        mx_div(i,j) = flowisentropic(gamma(i), a_div(1,j)/a_t,'sup');
        p_stat_div(i,j) = pc(i)*(1+(gamma(i)-1)/2*mx_div(i,j)^2)^(-gamma(i)/(gamma(i)-1)); %psi
    end
end
%% Chamber calcs
l_c_ch = l_c_ch*0.0254; %m
w_c_ch = w_c_ch*0.0254; %m
step_c = l_c/(nodes_c-1); %in
pos_c = 0:step_c:l_c; %in
d_c_out = d_c_in + 2*(th_in + w_c_ch/0.0254 + th_out); %in
a_c = pi*d_c_in^2/4*ones(1,nodes_c); %in^2
d_h_c = 4*l_c_ch*w_c_ch/(2*l_c_ch+2*w_c_ch)*ones(1,nodes_c); %m
l_c_ch_array = l_c_ch*ones(1,nodes_c); %m
w_c_ch_array = w_c_ch*ones(1,nodes_c); %m
for i=1:length(pc)
    for j=1:length(pos_div)
        mx_c(i,j) = flowisentropic(gamma(i), a_c(j)/a_t, 'sub');
        p_stat_c(i,j) = pc(i)*(1+(gamma(i)-1)/2*mx_c(i,j)^2)^(-gamma(i)/(gamma(i)-1)); %psi
    end
end
%% Converging calcs
angle_conv=tand(angle_conv);
step_conv=l_conv/(nodes_conv-1); %in
pos_conv=0:step_conv:l_conv; %in
r1 = 0;
r2= 0; 
tol = 1e-7;
error = 0.1;
while error>tol
    r1 = r1+error;
    xa = (angle_conv*r1)/((1+angle_conv)^0.5);
    ya = -(r1^2-xa^2)^0.5+d_t/2+r1;
    ya1 = angle_conv*(xa-l_conv/2)+(d_c_in+d_t)/4;
    error = ya-ya1;
end
tol=1e-5;
error=0.1;
while error>tol
    r2 = r2+error;
    xb = (l_conv+angle_conv^2*l_conv-angle_conv*r2*(angle_conv^2+1)^0.5)/(1+angle_conv^2);
    yb = (r2^2-(xb-l_conv).^2).^0.5+d_c_in/2-r2;
    yb1 = angle_conv*(xb-l_conv/2)+(d_c_in+d_t)/4;
    error = yb1-yb;
end
pos_conv_1 = -(r1^2-(pos_conv).^2).^0.5+d_t/2+r1; %in
pos_conv_2 = angle_conv*(pos_conv-l_conv/2)+(d_c_in+d_t)/4;
pos_conv_3 = (r2^2-(pos_conv-l_conv).^2).^0.5+d_c_in/2-r2; %in
slope_conv_l = (l_c_ch-th_l_ch)/(l_conv);
slope_conv_w = (w_c_ch-th_w_ch)/(l_conv);
d_conv = 2*cat(2,pos_conv_1(1:floor(xa/step_conv)),pos_conv_2(floor(xa/step_conv)+1:floor(xb/step_conv)),pos_conv_3(floor(xb/step_conv)+1:end)); %in
a_conv = pi*d_conv.^2/4; %in^2
l_conv_ch_array = slope_conv_l.*pos_conv+th_l_ch; %m
w_conv_ch_array = slope_conv_w.*pos_conv+th_w_ch; %m
d_h_conv = 4.*l_conv_ch_array.*w_conv_ch_array./(2*l_conv_ch_array+2*w_conv_ch_array); %m
for i=1:length(pc)
    for j=1:length(pos_conv)
        mx_conv(i,j) = flowisentropic(gamma(i), a_conv(1,j)/a_t, 'sub');
        p_stat_conv(i,j) = pc(i)*(1+(gamma(i)-1)/2*mx_conv(i,j)^2)^(-gamma(i)/(gamma(i)-1)); %psi
    end
end
%% Full Geometry
mx = [mx_div,mx_conv,mx_c];
p_stat = [p_stat_div,p_stat_conv,p_stat_c]; %psi
d_h = [d_h_div, d_h_conv, d_h_c]; %m
a = [a_div,a_conv,a_c]; %in^2
d = (4*a/pi).^0.5; %m
l_ch = [l_div_ch_array,l_conv_ch_array,l_c_ch_array]; %m
w_ch = [w_div_ch_array,w_conv_ch_array,w_c_ch_array]; %m
step = [step_div*ones(1,nodes_div), step_conv*ones(1,nodes_conv), step_c*ones(1,nodes_c)]; %in
pos = [pos_div,pos_conv+pos_div(nodes_div),pos_div(nodes_div)+pos_conv(nodes_conv)+pos_c]; %in
plot(pos,d)
axis equal