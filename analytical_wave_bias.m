function WB = analytical_wave_bias(z1,PL,PROC)
% ANALYTICAL_WAVE_BIAS - compute analytical expression for the wave-induced
% bias components 

WS = PL.WS; % recover the wave parameters
ADCP = PL.ADCP; % recover the ADCP parameters


% shorthand notation
k = WS.k;
a = WS.a;
omega = WS.omega;

z0 = PL.z0;
r = z1-z0; % range

gp = 1./(1-PL.U_TTW/WS.Cph);% gamma_p, orbital stretch due to platform propulsion
SPF = 0.5*(gp+1); % self-propulsion factor


US0 = SPF.*omega.*k.*a.^2;% surface Stokes drift
USp = SPF.*omega.*k.*a.^2.*exp(2*k.*z0); % platform Stokes drift


% tilt response factor?
switch PL.TILT
    case '0',
        gamma = 0;
    case 'M',
        gamma = 1;
    case 'G',
        gamma = -1;
end



if PROC.ADCP
    % ADCP Beam Geometry response functions
    th = ADCP.BEAM_ANGLE*pi/180;
    Ru = sin(th+k.*r*tan(th))/sin(th);
    Rw = cos(th+k.*r*tan(th))/cos(th);
    Rt = 2*sin(2*th+k.*r*tan(th))/sin(2*th); % =Ru+Rw
else
    % Perfect ("Vector") sampler
    Ru = 1;
    Rw = 1;
    Rt = 1;
end

% wave-induced platform motion bias, relative
WB.motion_bias_relative = US0.*exp(k.*z0).*(Ru.*exp(k.*z1) - exp(k.*z0)); 
        
% wave-induced platform tilt "sweeping" bias, relative
WB.sweep_bias_relative = 0.5*gamma* k .* r .* US0.*exp(k.*(z0+z1)).*Rt; % "sweeping" tilt response

% wave-induced platform tilt frame rotation bias, relative
WB.rotation_bias_relative  = 0.5*gamma.*WB.motion_bias_relative;

%% net bias - depends on the configuration
WB.net_bias_relative = 0;
switch PL.MOTION
        case 'L',
            % Lagrangian platform add motion bias
        WB.net_bias_relative = WB.net_bias_relative + WB.motion_bias_relative;
    case 'AUV',
        % Self-propelled platform add motion bias and TTW
        WB.net_bias_relative = WB.net_bias_relative + WB.motion_bias_relative-PL.U_TTW;
end

if PL.SWEEP,
    % add the platform tilt "sweeping" bias, relative
    WB.net_bias_relative = WB.net_bias_relative + WB.sweep_bias_relative;
end

if PL.FRAME_ROTATION,
    % add the platform tilt frame rotation bias, relative
    WB.net_bias_relative = WB.net_bias_relative + WB.rotation_bias_relative;
end


% absolute bias (includes the platform's own Stokes drift)
WB.net_bias_absolute = WB.net_bias_relative  + USp;


% include UStokes just in case
WB.UStokes = SPF.*omega.*k.*a.^2.*exp(2*k.*z1); 
