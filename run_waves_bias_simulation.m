% WAVE-INDUCED BIAS SIMULATION
% 
% This is the main script, it sets up and runs the simulation
% 

% 2024, Andrey Shcherbina (shcher@uw.edu)

clear

%% set up the wave field -- monochromatic here
disp('Monochromatic waves')
WS.a = 1;   % m, wave amplitude    
WS.t0 = 5;  % s, wave period
WS.ph = 0;  % rad, wave phase 
WS.fr = 1/WS.t0; % Hz, wave frequency
WS.NPERIODS = 5; % simulation length   
WS.direction = 1; % wave direction (i.e., sign of omega)

% Create the wave functions with these parameters
wave_functions;

%% Set up the platform parameters

% What is the platform's motion response? 
% Could be:
%   0 - platform is fixed in space (bottom-mounted)
%   L - fully-lagrangian platform (float)
%   AUV - self-propelled AUV/ASV (specify propulsion speed U_TTW below)
% (other motions can also be considered - but not here)

PL.MOTION = 'L';
PL.U_TTW = 0; % through-the-water speed (use with PL.MOTION = 'AUV')


% What is the platform tilt response?
% Could be '0', 'G', 'M':
%   0 - orientation is fixed
%   G - orient along the local gravity vector, i.e., normal to the isopotential surfaces - flat plates would do that (aks "hydrostatic").
%   M - orient along the material line that is vertical at rest (disregarding tilting due to Stokes shear) - spars would do that (aka "inertial"). 
PL.TILT = 'M'; 

% Consider some simplifications that we may want for testing:
PL.LINEARIZED = false; % whether to use linearized velocity estimation
PL.SWEEP = true; % whether to include the "sweeping" motion
PL.FRAME_ROTATION = true; % whether to include the frame rotation


PL.x0 = 0; % initial platform position (x)
% PL.z0 = 0; % initial platform position (z) will be defined later, based on the ADCP configuration 

%% set up the ADCP (and set the platform depth)
% Choose one of the pre-defined configuration scenarios - "DOWN","UP", or "UPDOWN"
ADCP.ORIENTATION = 'UPDOWN';

switch ADCP.ORIENTATION
    case 'UP',
        PL.z0 = -1/WS.k; % initial platform position (z)
        % Upward ADCP parameters
        ADCP.bins = [0:20]'/10*abs(PL.z0);% positive = upward
        ADCP.N_BEAMS = 5; 
        ADCP.BEAM_ANGLE = 25;
        fn = sprintf('wm-up-%s%c',PL.MOTION,PL.TILT); % create a file name
    case 'DOWN',
        PL.z0 = 0; % initial platform position (z)
        % Downward ADCP parameters
        ADCP.bins = -[0:20]'/10/WS.k; % positive = upward
        ADCP.N_BEAMS = 5; s
        ADCP.BEAM_ANGLE = 25;
        fn = sprintf('wm-down-%s%c',PL.MOTION,PL.TILT); % create a file name
    case 'UPDOWN',
        PL.z0 = -1/WS.k; %  initial platform position (z)
        % Up- & downward ADCP parameters
        ADCP.bins = [-10:1:10]'/10*abs(PL.z0); % positive = upward
        ADCP.N_BEAMS = 5; 
        ADCP.BEAM_ANGLE = 25;
        fn = sprintf('wm-updown-%s%c',PL.MOTION,PL.TILT);
end

ADCP = calculate_ADCP_beams(ADCP); % calculate the beam projection coefficients
nbins = length(ADCP.bins);


%% run the simulator
simulate_wave_bias;

%% save and plot the results
% save(fn,'-struct','PL');

plot_wave_bias_components(PL);




