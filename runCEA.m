function [Parameter,Geometry,Gas] = runCEA(Parameter,Geometry)
% RUNCEA Run CEAM.
% 
%   [Parameter,Geometry,Gas] = runCEA(Parameter,Geometry)
%   runs CEAM on a set of input parameters in struct
%   Parameters and geometric features in struct Geometry.
%
%   Function relies on CEA.p and thermo_lib.mat to function.

%% Pre-Defined Geometry
%  Function runs CEA once if engine geometry is pre-defined.

if Parameter.design == false
    % CEAM is run using CEA.p, outputs are stored in struct RKT1.
    RKT1=CEA('problem','rocket','equilibrium','fac','ma,kg/s',Parameter.mass_flow,'o/f',Parameter.O_F_ratio,'p(psi)',Parameter.chamber_pressure/6894.76,'pi/p',Parameter.chamber_pressure/Parameter.ambient_pressure,'supsonic(ae/at)',Geometry.exit_area/Geometry.throat_area,'reactants','fuel','C2H5OH(L)','wt%',Parameter.ethanol_ratio*100,'t(k)',Parameter.T,'fuel','H2O(L)','wt%',(1-Parameter.ethanol_ratio)*100,'t(k)',Parameter.T,'oxid','O2(L)','wt%',100,'t(k)',90.0,'output','transport','mks','end');
    
    % Characteristic velocity and thrust coefficient are taken from CEA
    % output
    Parameter.characteristic_velocity = RKT1.output.eql.cstar(1); % m/s
    Parameter.thrust_coefficient = RKT1.output.eql.cf(end);

    % Thurst is calculated from mass flow, characteristic veloicty, and
    % thrust coefficient, along with the corresponding efficiencies.
    Parameter.thrust = Parameter.mass_flow*Parameter.characteristic_velocity*Parameter.C_star_efficiency*Parameter.thrust_coefficient*Parameter.C_F_efficiency; % N

%% Engine Sizing
%  Function runs CEA in a loop to properly size the engine.

else
    % A while loop is set up to converge CEA inputs, since mass flow and
    % nozzle properties are not initially known.
    mass_flow_previous = 0; % kg/s
    
    % Variables are initialized with arbitrary values
    Parameter.mass_flow = 1; % kg/s
    Geometry.exit_area = 1; % m^2
    Geometry.throat_area = 1; % m^2

    % While loop convergence relies on CEA values matching those of the
    % previous iteration.
    while abs(mass_flow_previous-Parameter.mass_flow) > 1e-6
        % CEAM is run using CEA.p, outputs are stored in struct RKT1.
        RKT1=CEA('problem','rocket','equilibrium','fac','ma,kg/s',Parameter.mass_flow,'o/f',Parameter.O_F_ratio,'p(psi)',Parameter.chamber_pressure/6894.76,'pi/p',Parameter.chamber_pressure/Parameter.ambient_pressure,'supsonic(ae/at)',Geometry.exit_area/Geometry.throat_area,'reactants','fuel','C2H5OH(L)','wt%',Parameter.ethanol_ratio*100,'t(k)',Parameter.T,'fuel','H2O(L)','wt%',(1-Parameter.ethanol_ratio)*100,'t(k)',Parameter.T,'oxid','O2(L)','wt%',100,'t(k)',90.0,'output','transport','mks','end');
        
        % Characteristic velocity and thrust coefficient are taken from CEA
        % output
        Parameter.characteristic_velocity = RKT1.output.eql.cstar(1); % m/s
        Parameter.thrust_coefficient = RKT1.output.eql.cf(end);

        % Mass flow convergence is updated
        mass_flow_previous = Parameter.mass_flow; % kg/s

        % Mass flow rate is calculated from thrust, characterisic velocity,
        % thrust coefficient, along with the corresponding efficiencies.
        Parameter.mass_flow=(Parameter.thrust)/(Parameter.characteristic_velocity*Parameter.C_star_efficiency*Parameter.thrust_coefficient*Parameter.C_F_efficiency); % kg/s
        
        % Throat area is calculated from the same quantities along with
        % chamber pressure.
        Geometry.throat_area = (Parameter.characteristic_velocity*Parameter.C_star_efficiency*Parameter.mass_flow)/(Parameter.chamber_pressure); % m^2
        
        % Exit area is calculated using area ratio taken from CEA.
        Geometry.exit_area = Geometry.throat_area*RKT1.output.eql.aeat(4); % m^2
    end
end

%% Mass Flow
%  Applies O/F ratio

Parameter.fuel_mass_flow = Parameter.mass_flow*1/(1+Parameter.O_F_ratio);
Parameter.oxidizer_mass_flow = Parameter.mass_flow*(Parameter.O_F_ratio/(1+Parameter.O_F_ratio));

%% Gas Properties
%  Properties of the gas in the engine are recorded in struct Gas.

% Gas properties are taken from CEA at 4 stations: injector, converging,
% throat, and exit.
Gas.specific_heat_ratio = flip(RKT1.output.eql.gamma(1:4)'); 
Gas.specific_heat_capacity = flip(RKT1.output.eql.cp_tran.froz(1:4)')*1000; % J/kg-K
Gas.viscosity = flip(RKT1.output.eql.viscosity(1:4)')*1e-6; % Pa-s
Gas.Prandtl_number = flip(RKT1.output.eql.prandtl.froz(1:4)');

% Gas temperature at the throat is taken for future calculations.
Gas.throat_temperature = RKT1.output.eql.temperature(3)'; % K

end