% WAVE-INDUCED BIAS SIMULATION
% This subroutine defines a collection of functions that describe deep-water surface
% wave dynamics as functions of (x,z,t).
% It is expected that the wave field descriptor structure (WS) is already populated.
% The functions can handle a spectrum of waves if supplied with several
% frequency/amplitude pairs.
% 
% NB: Velocities are Eulerian, displacements - Lagrangian.
% NB: This is first-order approximation, but velocities are correct to
% O[(ak)^4]. 
%

% 2024, Andrey Shcherbina (shcher@uw.edu)
% Based, in part, on Eric D'Asaro's WaveEquations.pdf, following Phillips (1980).


% shift the wavenumber to the 4th dimension (for spectral calculations)
WS.a = shiftdim(WS.a(:),-3);
WS.ph = shiftdim(WS.ph(:),-3);

% calculate some constants
g = 9.81; % m/s^2
rho = 1e3; % kg/m^3
WS.omega =  WS.direction*2*pi*shiftdim(WS.fr(:),-3);
WS.k =   WS.omega.^2/g;
WS.Cph = WS.omega./WS.k;

% Stokes drift profile:
UStokes = @(z) WS.omega.*WS.k.*WS.a.^2.*exp(2*abs(WS.k).*z); % We always compute Stokes, but we never add it to the velocity (these are first-order equations!) Will need to add lated, if desired


% Wave phase:
PH = @(x,t) WS.k.*x - WS.omega.*t +  WS.ph;


% Free surface elevation:
Eta = @(x,t) sum( WS.a.*sin(PH(x,t)),4);


% Particle trajectories (first order, no Stokes)
Xsp = @(x,z,t) x +      WS.a .* exp(abs(WS.k) .* z) .* cos(PH( x ,t));
X = @(x,z,t) x +  sum(       WS.a .* exp(abs(WS.k) .* z) .* cos(PH( x ,t)),   4);

Zsp = @(x,z,t) z +  WS.a.* exp(abs(WS.k).*z) .*sin(PH( x ,t)); % spectral
Z = @(x,z,t) z + sum( WS.a.* exp(abs(WS.k).*z) .*sin(PH( x ,t)), 4);

% Particle trajectory gradients (i.e., material line deformation)
Xx = @(x,z,t) 1+sum( WS.a .* WS.k .* exp(abs(WS.k) .* z) .* -sin(PH(x,t)), 4);
Xz = @(x,z,t)   sum( WS.a .* WS.k .* exp(abs(WS.k) .* z) .* cos(PH(x,t)), 4); 

Zx = @(x,z,t) sum( WS.a .* WS.k .* exp(abs(WS.k) .* z) .* cos(PH(x,t)),   4);% isopotential slope (dz/dx)
Zz = @(x,z,t) 1+sum( WS.a .* WS.k .* exp(abs(WS.k) .* z) .* sin(PH(x,t)),   4); % vertical stretch


% Hydrostatic pressure (in dbar)
P = @(x,z,t) (-z + sum( WS.a.* exp(abs(WS.k).*z) .*sin(PH(x,t)), 4))*rho*g/1e4;

      
% Flow potential:
Q = @(x,z,t) sum( WS.a .* WS.omega./WS.k .* exp(abs(WS.k) .* z) .* -cos(PH(x,t)), 4); % NB: z is up!

% Velocities:
Usp = @(x,z,t)  WS.a .* WS.omega .* exp(abs(WS.k) .* z) .* sin(PH(x,t)); % spectral...
U = @(x,z,t) sum( WS.a .* WS.omega .* exp(abs(WS.k) .* z) .* sin(PH(x,t)), 4); 

Wsp = @(x,z,t)  WS.a .* WS.omega .* exp(abs(WS.k) .* z) .* -cos(PH(x,t)); % spectral...
W = @(x,z,t) sum( WS.a .* WS.omega .* exp(abs(WS.k) .* z) .* -cos(PH(x,t)), 4); 

% Velocity gradients
Ux = @(x,z,t) sum( WS.a .* WS.omega .* WS.k.* exp(abs(WS.k) .* z) .* cos(PH(x,t)), 4); 
Uz = @(x,z,t) sum( WS.a .* WS.omega .* WS.k.* exp(abs(WS.k) .* z) .* sin(PH(x,t)), 4); 

Wx = @(x,z,t) sum( WS.a .* WS.omega .* WS.k.* exp(abs(WS.k) .* z) .* sin(PH(x,t)), 4); 
Wz = @(x,z,t) sum( WS.a .* WS.omega .* WS.k.* exp(abs(WS.k) .* z) .* -cos(PH(x,t)), 4); 


% Lagrangian accelerations
Ax = @(x,z,t) sum( WS.a .* WS.omega .* WS.omega .* exp(abs(WS.k).*z) .* -cos(PH(x,t)), 4); % NB: z is up!
Az = @(x,z,t) sum( WS.a .* WS.omega .* WS.omega .* exp(abs(WS.k).*z) .* -sin(PH(x,t)), 4); % NB: z is up!



normalize = @(x) x./abs(x); % vector normalization function

% Acceleration vector (i.e. effective gravity), in complex notation (Re=horizontal, Im=vertical)
Gv = @(x,z,t) normalize(Ax(x,z,t)/g+1i*Az(x,z,t)/g+1i);

% Material "vertical" vector
Mv = @(x,z,t) normalize(Xz(x,z,t)+1i*Zz(x,z,t));


%% Complex shorthands for the functions
XZ = @(xz,t) X(real(xz),imag(xz),t) + 1i*Z(real(xz),imag(xz),t);
UW = @(xz,t) U(real(xz),imag(xz),t) + 1i*W(real(xz),imag(xz),t);
Axz = @(xz,t) Ax(real(xz),imag(xz),t) + 1i*Az(real(xz),imag(xz),t);

%% Linearized velocity
UWx = @(xz,t) Ux(real(xz),imag(xz),t) + 1i*Wx(real(xz),imag(xz),t);
UWz = @(xz,t) Uz(real(xz),imag(xz),t) + 1i*Wz(real(xz),imag(xz),t);

UW_linearized = @(xz,xz0,t) UW(xz0,t) + ...
                            UWx(xz0,t).*real(xz-xz0)+...
                            UWz(xz0,t).*imag(xz-xz0) ;
