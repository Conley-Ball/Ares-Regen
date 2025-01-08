function [Parameter,Geometry,Gas] = runGeometry(Parameter,Geometry,Gas)
% RUNGEOMETRY Generate and/or process engine geometry.
% 
%   [Parameter,Geometry,Gas] = runGeometry(Parameter,Geometry,Gas)
%   calculated geometry on a set of input parameters in struct
%   Parameters, geometric features in struct Geometry, and gas properties
%   in struct Gas.

%% Engine Sizing
%  Function constructs engine geometry from inputs.

if Parameter.design == true

    % Diameters are calculated.
    Geometry.throat_diameter = 2*sqrt(Geometry.throat_area/pi); % m
    Geometry.exit_diameter = 2*sqrt(Geometry.exit_area/pi); % m
    
    % Radii are calculated.
    throat_radius = Geometry.throat_diameter/2; % m
    exit_radius = Geometry.exit_diameter/2; % m
    chamber_radius = Geometry.chamber_diameter/2; % m

    %% Diverging Section
    
    % Initial and final angles for a Rao nozzle are calculated using
    % relations derived from the graph [1].
    initial_angle = -29/(Geometry.exit_area/Geometry.throat_area-1)^0.2+44.9; % deg
    final_angle=15/(Geometry.exit_area/Geometry.throat_area-1)^0.5+5.4; % deg
    
    % Diverging nozzle section is split into 2 sections: parabolic and 
    % circular (near throat).
    % The parabolic section length is found using 80% of a 15 degree
    % conical design as per Rao nozzle design [1].
    diverging_length_1 = 0.8*((Geometry.exit_area/Geometry.throat_area)^0.5-1)*throat_radius/tand(15); % m
    % The circular section length is calculated.
    diverging_length_2 = 0.682*throat_radius*sind(initial_angle); % m
    % Total diverging section length is the sum.
    Geometry.diverging_length = diverging_length_1+diverging_length_2; % m
    
    % Converging nozzle section is split into 3 sections: circular (near
    % throat), conical, and circular (near chamber).
    % The circular section length is calculated.
    converging_length_1 = 1.5*throat_radius*sind(Geometry.conv_angle); % m
    % The conical section length is calculated as a line with the given 
    % converging angle, and intersection with the other sections.
    converging_length_2 = (chamber_radius-throat_radius-(0.5+1.5)*throat_radius*(1-cosd(Geometry.conv_angle)))/tand(Geometry.conv_angle); % m
    % The circular section length is calculated.
    converging_length_3 = 0.5*throat_radius*sind(Geometry.conv_angle); % m
    % Total converging section length is the sum.
    Geometry.converging_length = converging_length_1+converging_length_2+converging_length_3; % m
    
    % Parabolic section geometry depends on a quadratic Bezier curve 
    % defined by point Q, the intersection of two lines with the Rao
    % angles. The point is calculated by a simple intersection.
    Q_x = ((throat_radius+0.682*throat_radius*(1-cosd(initial_angle))-tand(-initial_angle)*diverging_length_1)-exit_radius)/(tand(-final_angle)-tand(-initial_angle)); % m
    Q_r = tand(-final_angle)*Q_x+exit_radius; % m
    
    % The parabolic curve is defined as parametric with variable t. The
    % number of nodes for this section is calculated to ensure even
    % distribution across the whole engine.
    t_1=linspace(0,1,ceil(0.5*Parameter.resolution*diverging_length_1/(Geometry.converging_length+Geometry.diverging_length)));
    % The x and r functions are calculated.
    x_1 = (1-t_1).^2.*diverging_length_1+2.*(1-t_1).*t_1.*Q_x+t_1.^2*0; % m
    r_1 = (1-t_1).^2.*(throat_radius+0.682.*throat_radius.*(1-cosd(initial_angle)))+2.*(1-t_1).*t_1.*Q_r+t_1.^2.*exit_radius; % m
    
    % The circular curve is defined as parametric with variable t.
    t_2 = linspace(initial_angle,0,ceil(0.5*Parameter.resolution*diverging_length_2/(Geometry.converging_length+Geometry.diverging_length)));
    % The x and r functions are calculated.
    x_2 = Geometry.diverging_length-0.682*throat_radius*sind(t_2); % m
    r_2 = throat_radius+0.682*throat_radius*(1-cosd(t_2)); % m
    
    % The two sections are combined to form the complete diverging section 
    % of the nozzle. The parabolic section was calculated right-to-left and
    % is thus flipped.
    x_diverging = [flip(x_1) x_2(2:end)]; % m
    r_diverging = [flip(r_1) r_2(2:end)]; % m
    
    %% Converging Section
    
    % The circular curve is defined as parametric with variable t.
    t_1 = linspace(0,Geometry.conv_angle,ceil(0.5*Parameter.resolution*converging_length_1/(Geometry.converging_length+Geometry.diverging_length)));
    % The x and r functions are calculated.
    x_1 = Geometry.diverging_length+1.5*throat_radius*sind(t_1); % m
    r_1 = throat_radius+1.5*throat_radius*(1-cosd(t_1)); % m
    
    % The linear section is calculated as a simple y=mx+b
    x_2 = linspace(Geometry.diverging_length+converging_length_1,Geometry.diverging_length+converging_length_1+converging_length_2,ceil(0.5*Parameter.resolution*converging_length_2/(Geometry.converging_length+Geometry.diverging_length))); % m
    r_2 = throat_radius+1.5*throat_radius*(1-cosd(Geometry.conv_angle))+tand(Geometry.conv_angle)*(x_2-converging_length_1-Geometry.diverging_length); % m
    
    % The circular curve is defined as parametric with variable t.
    t_3 = linspace(Geometry.conv_angle,0,ceil(0.5*Parameter.resolution*converging_length_3/(Geometry.converging_length+Geometry.diverging_length)));
    % The x and r functions are calculated.
    x_3 = Geometry.diverging_length+Geometry.converging_length-0.5*throat_radius*sind(t_3); % m
    r_3 = chamber_radius-0.5*throat_radius*(1-cosd(t_3)); % m
    
    % The three sections are combined to form the complete converging 
    % section of the nozzle.
    x_converging = [x_1 x_2(2:end) x_3(2:end)]; % m
    r_converging = [r_1 r_2(2:end) r_3(2:end)]; % m
    
    %% Chamber Section

    % The volume of the converging section is calculated by use of a line
    % integral
    converging_volume = sum(pi*r_converging.^2*Geometry.converging_length/(length(r_converging)-1)); % m^3
    % The volume of the chamber is calculated using the characteristic
    % length and throat area. Since combustion is considered in the
    % converging section, it is subtracted.
    chamber_volume = Parameter.L_star*Geometry.throat_area-converging_volume; % m^3
    % The length of the chamber section is calculated.
    Geometry.chamber_length = chamber_volume/(pi*chamber_radius^2); % m
    
    % The x and r values for the chamber are defined with half of the total
    % number of nodes.
    x_chamber = linspace(Geometry.diverging_length+Geometry.converging_length,Geometry.diverging_length+Geometry.converging_length+Geometry.chamber_length,0.5*Parameter.resolution); % m
    r_chamber = chamber_radius*ones(1,0.5*Parameter.resolution); % m
    
    % The axial coordinates are combined, excluding overlapping values.
    Geometry.axial_distance = [x_diverging x_converging(2:end) x_chamber(2:end)]; % m
    % The radial coordinates are combined similarly.
    r = [r_diverging r_converging(2:end) r_chamber(2:end)]; % m
    % Diameter is calculated.
    Geometry.diameter = 2*r; % m

    % Various station indices are calculated for plotting: exit, throat,
    % converging, injector.
    Geometry.stations = [1,length(x_diverging),length(x_diverging)+length(x_converging)-1,length(r)]; % m


%% Pre-Defined Geometry
%  Function interprets raw data from the imported geometry.

else

    % Diameters are calculated.
    Geometry.throat_diameter = 2*sqrt(Geometry.throat_area/pi); % m
    Geometry.exit_diameter = 2*sqrt(Geometry.exit_area/pi); % m
    Geometry.chamber_diameter = Geometry.diameter(end); % m
    
    % Radii are calculated.
    throat_radius = Geometry.throat_diameter/2; % m
    exit_radius = Geometry.exit_diameter/2; % m
    chamber_radius = Geometry.chamber_diameter/2; % m
    r = Geometry.diameter/2; % m

    % The throat is located as the minimum radius.
    [~,throat_index] = min(Geometry.diameter); % m
    % The converging section is located by ignoring the chamber and 
    % finding the length.
    converging_index = length(Geometry.diameter(Geometry.diameter~=Geometry.chamber_diameter)); % m
    % Various station indices are calculated for plotting: exit, throat,
    % converging, injector.
    Geometry.stations = [1,throat_index,converging_index+1,length(Geometry.diameter)]; % m

    % Section lengths are calculated.
    Geometry.diverging_length = Geometry.axial_distance(throat_index); % m
    Geometry.converging_length = Geometry.axial_distance(converging_index+1)-Geometry.diverging_length; % m
    Geometry.chamber_length = Geometry.axial_distance(end)-Geometry.converging_length-Geometry.diverging_length; % m

    % Chamber volume is calculated.
    converging_volume = sum(pi*r(throat_index+1:converging_index+1).^2*Geometry.converging_length/(converging_index-throat_index+1)); % m^3
    chamber_volume = Geometry.chamber_length*pi*chamber_radius^2; % m^3
    
    % Characterisic length is calculated.
    Parameter.L_star = (converging_volume+chamber_volume)/Geometry.throat_area; % m
    
end

% Engine area is calculated.
Geometry.area = pi*r.^2; % m^2

% Number of channels is calculated from the circumference at the throat,
% and the minimum rib width.
throat_circumference = 2*pi*(throat_radius+Geometry.inner_wall_thickness(2)); % m
Geometry.channel_number = floor(throat_circumference/(Geometry.rib_width(2)+Geometry.minimum_channel_width));

%% Channel Geometries

% Variable X is defined with the axial distances at each station for use in
% interpolation.
X = [0,Geometry.diverging_length,Geometry.diverging_length+Geometry.converging_length,Geometry.diverging_length+Geometry.converging_length+Geometry.chamber_length]; % m
% Channel geometries are defined as a linear interpolation between each
% station.
Geometry.inner_wall_thickness = interp1(X,Geometry.inner_wall_thickness,Geometry.axial_distance); % m
Geometry.outer_wall_thickness = interp1(X,Geometry.outer_wall_thickness,Geometry.axial_distance); % m
Geometry.channel_height = interp1(X,Geometry.channel_height,Geometry.axial_distance); % m
Geometry.rib_width = interp1(X,Geometry.rib_width,Geometry.axial_distance); % m

% Circumference is calculated.
Geometry.circumference = 2*pi*(r+Geometry.inner_wall_thickness); % m
% The channel width is calculated as from the rib width an circumference
Geometry.channel_width = Geometry.circumference/Geometry.channel_number-Geometry.rib_width; % m

Geometry.channel_area = Geometry.channel_width.*Geometry.channel_height + (2*Geometry.fillet)^2 - pi*(Geometry.fillet)^2;
channel_perimeter = 2*(Geometry.channel_width+Geometry.channel_height-4*Geometry.fillet) + 2*pi*Geometry.fillet;
Geometry.hydraulic_diameter = 4*Geometry.channel_area./channel_perimeter;

%% Gas Parameters

% Gas properties at each station are linearly interpolated in the same
% method as the channel geometries.
Gas.specific_heat_ratio = interp1(X,Gas.specific_heat_ratio,Geometry.axial_distance);
Gas.specific_heat_capacity = interp1(X,Gas.specific_heat_capacity,Geometry.axial_distance); % J/kg-K
Gas.viscosity = interp1(X,Gas.viscosity,Geometry.axial_distance); % Pa-s
Gas.Prandtl_number = interp1(X,Gas.Prandtl_number,Geometry.axial_distance);

%% Channel Angles

% Axial step is calculated.
dx = diff(Geometry.axial_distance); % m
% Since the final value is linear, the value can be copied.
Geometry.dx = [dx dx(end)]; % m

% Radial step is calculated.
dr = diff(r); % m
dr = [dr dr(end)]; % m

% Angle alpha represents the pitch angle for helical channels.
Geometry.alpha = Geometry.pitch*pi/180; % rad

% Theta step is calculated in cylindrical coordinates.
dtheta = tan(Geometry.alpha)./r.*(Geometry.dx.^2+dr.^2).^0.5; % rad

% Angle psi represents the angle of the engine contour in the r-x plane.
Geometry.psi = atan(dr./Geometry.dx); % rad

% Angle phi represents the angle of a channel with respect to the main
% coordinate system axis.
Geometry.phi = atan2((dr.^2+(r.*dtheta).^2).^0.5,Geometry.dx); % rad

% The maximum such angle is calculated for additive manufacturing
% considerations.
Geometry.maximum_angle = max(Geometry.phi); % rad

end

% References:
% [1] Newlands, Rick. The Thrust Optimised Parabolic nozzle