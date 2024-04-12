function uw = infer_velocity(PL, PROC)
% INFER_VELOCITY - Infer the velocities from the "virtual platform observations" created by "make_waves"
%   This function replicates the ADCP velocity processing (including
%   referencing, if required).
%   It can also perform the "Vector" processing that avoids the issues
%   arising from the ADCP beam geometry.
%   Call as
%       uw = infer_velocity(PL, PROC)
%   where
%       PL is the "platform" structure created by make_waves, and
%       PROC is the processing parameter structure:
%   PROC.ADCP - if true, do the ADCP beam processing; otherwise, do the idealized sampler ("vector") processing
%   PROC.TILT_BIN_MAPPING - if 0, no bin mapping
%                           if 1, apply normal (RDI-style) bin mapping
%                           if 2, apply bin interpolation (Ott 2002, doi:10.1175/1520-0426(2002)019<1738:AIITCO>2.0.CO;2)
%                           This option has no effect if PROC.ADCP = 0
%   PROC.EARTH_COORDINATES - if true, return velocities rotated to Earth coordinates; otherwise, in "ship" coordinates
%   PROC.U_REFERENCE       - if true, reference U using platform's velocity (assumed known); otherwise, U is relative
%   PROC.W_REFERENCE       - if true, reference W using platform's velocity (assumed known); otherwise, W is relative 

% 2024, Andrey Shcherbina (shcher@uw.edu)


if PROC.ADCP
    %% ADCP beam processing

    % start with the relative beam velocities (pick the two slanted beams)
    vr1 = PL.vr(:,:,1);
    vr2 = PL.vr(:,:,2);
    nbin = size(PL.xz_eu,2);

    switch PROC.TILT_BIN_MAPPING
        case 1,
            %% Bin (cell) mapping
            % We need to figure out which bin pairs should be considered in claculations at each nominal level
            % (we could do this the other way around, i.e. where does each bin go, but it would be less straightforward later on).
            % This is saved as the "bin index (ib)", with the dimesions of [time, nominal range, beam].

            % First, add a row of NaNs,to use for out-of-range bins
            vr1(:,end+1) = NaN;
            vr2(:,end+1) = NaN;

            % Because of the tilting, other cells may be better aligned with the nominal levels.
            zb = sq(imag(PL.bxz(:,:,1))); % nominal levels, relative to the platform

            % ensure that the cells are equally spaced, so we could be more efficient
            h = diff(zb);
            if max(abs(h-h(1)))/abs(h(1))>1e-3,
                error('Bin mapping is not (yet) implemented for uneven bin spacing')
            end
            h = h(1);
            za = imag(PL.bxzi(:,:,1:2)-PL.xz); % actual relative height of each bin of each beam
            dza = za - zb; % deviation of the actual bin relative to the nominal level

            % modify the bin index: if a bin goes up, the index should go down
            ib = ib - (dza/h);

            ib = round(ib);
            ib(ib<1 | ib>nbin) = nbin+1; % refer out-of-range bins to the row of NaNs

            i1 = colsub2ind(ib(:,:,1));
            i2 = colsub2ind(ib(:,:,2));

            % apply cell sorting
            vr1i = vr1(i1);
            vr2i = vr2(i2);
        case 2,
            % implement the "improved" bin mapping described in  Ott 2002 (doi:10.1175/1520-0426(2002)019<1738:AIITCO>2.0.CO;2)
            % We need to interpolate the beam velocities onto the nominal levels.

            zb = sq(imag(PL.bxz(:,:,1))); % nominal levels, relative to the platform
            za = imag(PL.bxzi(:,:,1:2)-PL.xz); % actual relative height of each bin of each beam

            % this can probably be vectorized...
            vr1i = vr1*iNaN;
            vr2i = vr2*iNaN;
            for k=1:size(za,1)
                vr1i(k,:) = interp1(za(k,:,1),vr1(k,:),zb,'*pchip');
                vr2i(k,:) = interp1(za(k,:,2),vr2(k,:),zb,'*pchip');
            end
        otherwise,
            % no bin mapping
            vr1i = vr1;
            vr2i = vr2;
    end


%% Infer the UV velocity from the beam velocities using the beam transformation coefficients

% velocity inference
uw =    vr1i*PL.ADCP.a11 + vr2i*PL.ADCP.a12+...
    1i*(vr1i*PL.ADCP.a21 + vr2i*PL.ADCP.a22);


% Do we need to rotate to Earth coordinates?
if PROC.EARTH_COORDINATES
    % Rotate to Earth coordinates, still relative to platform
    uw = uw.*PL.axis/1i;
end


% do we know platform's horizontal velocity for referencing?
if PROC.U_REFERENCE
    % yes, reference...
    uw = uw + real(PL.uw0);
end


% do we know platform's vertical velocity for referencing?
if PROC.W_REFERENCE
    % yes, reference...
    uw = uw + 1i*imag(PL.uw0);
end

else
    %% vector processing
    % use the true velocity at the sampling cells corresponding to the
    % middle beam; ignore the slanted beams.
    % There is no cell mapping.
    uw = PL.uw(:,:,3) - PL.uw0; % relative velocity profile

    % do we know platform's horizontal velocity for referencing?
    if PROC.U_REFERENCE
        % yes, reference...
        uw = uw + real(PL.uw0);
    end


    % do we know platform's vertical velocity for referencing?
    if PROC.W_REFERENCE
        % yes, reference...
        uw = uw + 1i*imag(PL.uw0);
    end

    % coordinate system rotation?
    % This is different from the "ADCP" estimator: There, after the beam calculations,
    % velocities are in the instrument frame of reference and may need to be
    % rotated to Earth frame.
    % Here, we are _already_ in Earth frame, but we may need to rotate to instrument coordinates.

    if ~PROC.EARTH_COORDINATES
        theta = angle(PL.axis/1i);
        uw = uw.*exp(-i*theta);
    end

end
