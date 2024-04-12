% Time stepping the model platform behavior in the wave field
% 
%   The following structures need to be defined:
%   WS - wave field
%   PL - platform response
%   ADCP - ADCP configuration
%
% Here, the platform trajectory is first-order (linear), so we don't need to
% integrate. 
% We don't need to account for Doppler shift due to the Stokes, either.
% (But we do need to consider Doppler shift due to the TTW velocity.)
% This is self-consistent to the second order (I think). 
%
% This code works a bit differently from the real ADCP sampling (for
% efficiency). The main difference is in how the UP/DOWN beams and cells
% are handled. Here, I am allowing negative ranges. When this happens, beam
% direction remains the same (upward). 
% Because of that, the same beam (#1) always points towards +X and +Z
% Therefore, U = a*(B1-B2)


% store the wave parameters and the ADCP setup
PL.ADCP = ADCP;
PL.WS = WS;



% Need to Doppler-shift the period so that the instrument samples full number of cycles (in its own Lagrangian frame)
% This is only necessary for TTW propulsion velocity, *not* the Stokes 

% From Longuet-Higgins (1986), Lagrangian period is 
% T_L = L/(c-U) = 2*pi/k/(omega/k - U) = 2*pi/(omega - k*U) 

if strcmp(PL.MOTION,'AUV')
    Td = 2*pi/( min(abs(WS.omega)) - min(abs(WS.k))*PL.U_TTW );
    fprintf('Doppler shifting of wave period: %.1f->%.1f\n',WS.t0,Td);
else
    Td = max(WS.t0);
end
Td = abs(Td);

%% save the output variable descriptions
PL.notes = { ...
    'xz0 - platform trajectory without the waves'; ...
    'xz - platform trajectory with the waves'; ...
    'xz_eu - nominal velocity measurement locations (without the motion/tilt/beam effects)';...
    'bxz - nominal xz positions of each cell of each beam (without the motion/tilt!), * relative to the instrument* ';...
    'bxz0 - nominal xz positions of each cell of each beam (without the motion/tilt!), * absolute, accounting for platform propulsion *';...
    'bxzi - tilted and advected xz positions of the ADCP bins, for each beam';...
    'uw0 - platform velocity';...
    'uw_eu - true Eulerian velocity at the nominal measurement locations (xz_eu)';... 
    'uw - true Eulerian velocity at each cell of each beam (bxzi)';...
    'vr_abs - radial beam velocities at each cell of each beam, absolute';...
    'vr - radial beam velocities at each cell of each beam, relative to the instrument -- this is the only thing a real ADCP would measure';...
    };
  
%% Establish time axis
dt = min(WS.t0)/80;
t = (0:dt:(WS.NPERIODS*Td-dt))';

nt = length(t);

PL.t = t(:); % time in seconds


    
    
    %% Handle platform motion
    % Here, I initialize PL.xz0 based on the steady platform advection, if any (the waves are added on top of that).
    % Stokes drift does not need to be included (to the first order).
    % TODO: Add provisions for profiling motions!
    
    PL.xz0 = PL.x0 + 1i*PL.z0 + 0*t; % start with no motion

    switch PL.MOTION
        case '0', % fixed in space
            PL.xz = PL.xz0;
            PL.uw0 = 0*t;
        case 'L',
            % the platform follows the orbital motions
            PL.xz = XZ(PL.xz0,t);
            PL.uw0 = UW(PL.xz,t);
        case 'AUV',
            % self-propelled AUV/ASV (through the water speed is PL.U_TTW)
            % use the "local gravity" vector as a proxy of local surface slope
            PL.xz0 = PL.xz0 + PL.U_TTW*t;
            PL.xz = XZ(PL.xz0,t);
            % fix the orbitals:
            gp = 1./(1-PL.U_TTW/WS.Cph);% gamma_p, orbital stretch due to platform propulsion
            PL.xz = PL.xz-PL.xz0;
            PL.xz = real(PL.xz)*gp+i*imag(PL.xz)+PL.xz0;
            PL.uw0 = UW(PL.xz0,t) + UStokes(PL.z0)+PL.U_TTW; 
    end
   
%% Establish the effective gravity and material line vectors
% Note that we use the "no-waves" positions here (first order!) 
PL.Gv = Gv(real(PL.xz0),imag(PL.xz0),t); % save gravity
PL.Mv = Mv(real(PL.xz0),imag(PL.xz0),t); % save vertical
PL.Acceleration = Axz(PL.xz0,t);


%% Handle the platform tilt
switch PL.TILT,
    % In the first-order approximation,
    % the platform can follow either the vertical or the horizontal material line
    % (there's probably frequency-dependent phase lag in realty, too)
    % so there are two possible platform axes:
    case 'G',
        PL.axis = PL.Gv;
    case 'M',
        PL.axis = PL.Mv;
    otherwise
        PL.axis = 1i+t*0; % stay vertical
end

%% Establish ADCP beams and cells

nbins = size(ADCP.bins,1);
bxz = shiftdim(ADCP.bins*(ADCP.bv),-1); % Nominal xz positions of each cell of each beam (without the motion/tilt!), * relative to the instrument*
bxz0  = bxz + PL.xz0;% Nominal xz positions of each cell of each beam (without the motion/tilt!), * absolute, accounting for platform propulsion *

bv = normalize(shiftdim(ADCP.bv,-1)); % nominal beam vectors

if PL.FRAME_ROTATION
    bvi = bv.*PL.axis/1i; % tilted beam vectors
else
    bvi = bv; % no beam vector rotation
end

if PL.SWEEP,
    bxzi = bxz.*PL.axis/1i + PL.xz; % tilted and advected ADCP bins, for each beam
else
    bxzi = bxz + PL.xz;     % no sweep -- just advected ADCP bins, for each beam
end

% bxzi0 = mean(bxzi,3); % this is the nominal velocity estimate location (midpoints between the beams)

% save coordinates
PL.xz_eu = PL.xz0+i*ADCP.bins'; % Nominal velocity measurement locations (without the motion/tilt)
PL.bxz = bxz;   % Nominal xz positions of each cell of each beam (without the motion/tilt), * relative to the instrument* !
PL.bxz0 = bxz0; % Nominal xz positions of each cell of each beam (without the motion/tilt), * in space * !
PL.bxzi = bxzi; % Tilted and advected xz positions of the ADCP bins, for each beam




%% Handle velocity sampling
% Start with the "true" velocities at the beam cells
% We should use the first-order formula U(...), W(...) to obtain
% Eulerian velocities
% Note also, that U (Eulerian) will be different from dx/dt (Lagrangian)

if PL.LINEARIZED,
    % linearized velocity field (for testing)
    uwi =   UW_linearized(PL.bxzi,PL.bxz0,t); % true velocities at each cell
else
    uwi =   UW(PL.bxzi,t); % true velocities at each cell
end

% % Perhaps blank the measurements that happen to be above the surface?
% % Where *is* the surface?
% zs = Eta(real(bxzi),t);
% out = imag(PL.bxzi)>zs;
% uwi(out) = iNaN;

PL.uw = uwi; % save true vel. @ the bins (its average will already be somewhat biased re. Eulerian!)
PL.vr_abs = real(PL.uw.*conj(bvi)); % radial absolute velocities

uwr = uwi - PL.uw0; % relative velocities -- remove platform motion
vr = real(uwr.*conj(bvi)); % radial relative velocities (i.e., beam velocities)
PL.vr = vr;

% Eulerian ("true") velocities at the nominal cell locations, for the reference 
PL.uw_eu = UW(PL.xz_eu,t);
