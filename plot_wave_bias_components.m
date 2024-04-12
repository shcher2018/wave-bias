function plot_wave_bias_components(PL)
% PLOT_WAVE_BIAS_* - plot the results of wave-induced bias simulation.
% This version displays the individual components of the bias
% Call:
%   plot_wave_bias_components(file) or plot_wave_bias_components(PL)
% where PL is the structure created by make_waves.

if ~isstruct(PL)
    PL = load(PL);
end
z1 = imag(PL.xz_eu(1,:));

% set up the processing accordingly
PROC.ADCP = true; % "ADCP" processing
PROC.EARTH_COORDINATES = false;
PROC.U_REFERENCE = false;
PROC.W_REFERENCE = false;
PROC.TILT_BIN_MAPPING = 0;

% infer the "measured" velocities
uw = infer_velocity(PL,PROC); % this will take care of bin-mapping, referencing, rotation, etc. (if required by parameters in PROC)

% wave-average to get the bias
guw = mean(uw,1);

% Eulerian mean (should be zero, so a sanity check)
euw = squeeze(mean(PL.uw_eu));

% analytical expression
WB = analytical_wave_bias(z1,PL,PROC);



%%
figure
clf

plot(WB.UStokes,z1,'r--','linewidth',1);
hold on;
grid on
plot(WB.net_bias_relative,z1,'b','linewidth',1);

plot(WB.motion_bias_relative,z1,'m-','linewidth',0.5);
plot(WB.sweep_bias_relative,z1,'m--','linewidth',0.5);
plot(WB.rotation_bias_relative,z1,'m:','linewidth',0.5);
plot(real(guw),z1,'g+');

plot(xlim,PL.z0*[1 1],'k--')
legend({'Stokes','Analytical bias','... motion','... sweeping','... rotation','Simulation','Platform depth'},'location','SE');
xlabel('U_{bias} (m/s)')
ylabel('z (m)');
title ('Wave-induced bias, relative , ADCP processing')
