function design = designMatrix_BH(cfg0, phi)
% [design] = designMatrix_BH(cfg, phi)
%    Returns hypothetical channel responses given a presented orientation, cf. Brouwer & Heeger.
%
%    phi         Array of length N, where N is the number of trials, that specifies the presented
%                orientation on each trial. Orientation is specified in degrees and is expected to
%                have a range of 0-180.
%    cfg         Configuration struct that can possess the following fields:
%                .numC                               The number of hypothetical channels C to use. The
%                                                    channels are equally distributed across the circle,
%                                                    starting at 0.
%                .tuningCurve                        The tuning curve according to which the responses
%                                                    are calculated.
%                .tuningCurve = 'vonMises'           Von Mises curve. Kappa: concentration parameter.
%                .tuningCurve = 'halfRectCos'        Half-wave rectified cosine. Kappa: power.
%                .tuningCurve = [function_handle]    User-specified function that can take a matrix as input,
%                                                    containing angles in radians with a range of 0-pi.
%                .kappa                              Parameter(s) to be passed on to the tuning curve.
%                .offset                             The orientation of the first channel. (default = 0)
%           
%    design      The design matrix, of size C x N, containing the hypothetical channel responses.
%
%    See also designMatrix_dummy designMatrix_discriminant

%    Created by Pim Mostert, 2016
phi = phi(:);

numN = length(phi);

switch cfg0.tuningCurve
    case 'vonMises'
        cfg0.tuningCurve = @(X) vonMises_curve(X*2, cfg0.kappa);
    case 'halfRectCos'
        cfg0.tuningCurve = @(X) abs(cos(X)).^cfg0.kappa;       
end    

% Apply offset
if ~isfield(cfg0, 'offset')
    cfg0.offset = 0;
end

phi = phi - cfg0.offset;

% Create matrix
design = (0:(cfg0.numC-1))'*ones(1, numN) * (180/cfg0.numC);
design = design - ones(cfg0.numC, 1)*phi';
design = design * (pi/180);
design = cfg0.tuningCurve(design);

end


