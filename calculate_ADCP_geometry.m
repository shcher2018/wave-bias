function ADCP = calculate_ADCP_beams(ADCP)
% CALCULATE_ADCP_BEAMS - calculate the ADCP beam vectors and velocity projection
%   coefficients based on the beam angle and number of beams. 
%   
%   ADCP = calculate_ADCP_beams(ADCP)
%   

% 2024, Andrey Shcherbina (shcher@uw.edu)


% Since the simulation is 2D, we only need the beams in X-Z plane.
% Even if an ADCP doesn't have a vertical beam (e.g., N_BEAMS==4), I'll
% include it anyway since it will be used in "vector" calculations. 

switch ADCP.N_BEAMS,
    case 3,
        % beam vectors NB: these point *away*
        % Beam 1 points forward,
        % Beam 3 is straight up
        ADCP.bv = [1 -0.5 0]*tan(ADCP.BEAM_ANGLE/180*pi) + 1i; % not yet normalized!
        
        % beam transformation matrix (for 3-beam instrument)
        % |u|   | a11 a12 |   | vr1 |
        % | | = |         | * |     |
        % |w|   | a21 a22 |   | vr2 |
        
        % the actual Aquadopp's transformation matrix is
        % 1.5774 -0.7891 -0.7891
        % 0.0000 -1.3662 1.3662
        % 0.3677 0.3677 0.3677
        % so we're basically add the last two columns together (assuming
        % the two beams pointing "back" measure the same thing):
        ADCP.a11 =   1.5774677271977;
        ADCP.a12 = -1.4680160396636;
        ADCP.a21 =  0.367792637277844;
        ADCP.a22 =  0.684547121589507;
    case {4,5}
        % beam vectors NB: these point *away*
        % Beam 1 points forward,
        % Beam 3 is straight up
        ADCP.bv = [1 -1 0]*tan(ADCP.BEAM_ANGLE/180*pi) + 1i; % not yet normalized!
        
        % beam transformation matrix (for 2-beam instrument)
        % |u|   | a11 a12 |   | vr1 |   | a -a |   | vr1 |
        % | | = |         | * |     | = |      | * |     |
        % |w|   | a21 a22 |   | vr2 |   | b  b |   | vr1 |
        
        a = 1/(2*sin(ADCP.BEAM_ANGLE/180*pi));
        b = 1/(2*cos(ADCP.BEAM_ANGLE/180*pi));
        ADCP.a11 = a;
        ADCP.a12 = -a;
        ADCP.a21 = b;
        ADCP.a22 = b;
        
end
