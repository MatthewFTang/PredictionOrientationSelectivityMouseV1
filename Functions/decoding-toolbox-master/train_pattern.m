function decoder = train_pattern(cfg0, X, Y)
% [decoder] = train_pattern(cfg, X, Y)
%    Linearly estimates the feature pattern corresponding to activity in a latent component as
%    prescribed in X. Several patterns may be trained indepedently, corresponding to several
%    latent components.
%
%    X           Vector or matrix of size C x N, where C is the number of components and N is
%                the number of trials, that contains the expected/prescribed component activity
%                in the training data.
%    Y           Matrix of size F x N, where F is the number of features, that contains the
%                training data.
%    cfg         Configuration struct that can possess the following fields:
%                .discardNan = 'yes' or 'no'      Whether trials with NaN in either X or Y should be
%                                                 removed prior to training. Default is 'no'.
%                .demean = 'yes' or 'no'          Whether to demean the data first (per feature, over
%                                                 trials). Default = 'yes';.
%
%    decoder     The (set of) decoder(s), that may be passed on to an appropriate decoding function,
%                e.g. decode_beamformer. It may contain a field .pattern of size C x F
%
%    See also DECODE_PATTERN.

%    Created by Pim Mostert, 2017

if ~isfield(cfg0, 'discardNan')
    cfg0.discardNan = 0;
end
if ~isfield(cfg0, 'demean')
    cfg0.demean = 'yes';
end    

decoder = [];

% Tranpose
Y = Y';
X = X';

if cfg0.discardNan
    X = X(~any(isnan(Y), 2), :);
    Y = Y(~any(isnan(Y), 2), :);
end

numF = size(Y, 2);
numC = size(X, 2);
numN = size(Y, 1);

% Demean dataand design matrix
if strcmp(cfg0.demean, 'yes')
    mY = mean(Y, 1);
    Y = Y - repmat(mY, [numN, 1]);
    
    decoder.mY = mY';
end

% Estimate patterns
decoder.W = zeros(numC, numF);

for ic = 1:numC
    % Estimate leadfield for current channel    
    decoder.W(ic, :) = (X(:, ic)'*X(:, ic))\(X(:, ic)'*Y);
end

end
