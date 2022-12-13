function decoder = train_beamformer(cfg0, X, Y)
% [decoder] = train_beamformer(cfg, X, Y)
%    Trains a linear decoder "beamformer style" to optimally recover the latent components as
%    prescribed in X. Several decoders may be trained indepedently, corresponding to several
%    latent components.
%
%    X           Vector or matrix of size C x N, where C is the number of components and N is
%                the number of trials, that contains the expected/prescribed component activity
%                in the training data.
%    Y           Matrix of size F x N, where F is the number of features, that contains the
%                training data.
%    cfg         Configuration struct that can possess the following fields:
%                .gamma = [scalar]                Shrinkage regularization parameter, with range [0 1].
%                                                 No default given.
%                .discardNan = 'yes' or 'no'      Whether trials with NaN in either X or Y should be
%                                                 removed prior to training. Default is 'no'.
%                .returnPattern = 'yes' or 'no'   Whether the spatial patterns of the components should
%                                                 be returned. Default = 'no';
%                .demean = 'yes' or 'no'          Whether to demean the data first (per feature, over
%                                                 trials). Default = 'yes';.
%
%    decoder     The (set of) decoder(s), that may be passed on to an appropriate decoding function,
%                e.g. decode_beamformer. It may contain a field .pattern of size C x F
%
%    See also DECODE_BEAMFORMER.

%    Created by Pim Mostert, 2016

if ~isfield(cfg0, 'discardNan')
    cfg0.discardNan = 0;
end
if ~isfield(cfg0, 'returnPattern')
    cfg0.returnPattern = 'no';
end
if ~isfield(cfg0, 'demean')
    cfg0.demean = 'yes';
end

decoder = [];

% Tranpose
Y = Y'; % makes the data trials x channels
X = X';

if cfg0.discardNan
    X = X(~any(isnan(Y), 2), :);
    Y = Y(~any(isnan(Y), 2), :);
end

numF = size(Y, 2); % N channels
numC = size(X, 2); % N orientation channels
numN = size(Y, 1); % N trials

% Demean dataand design matrix
if strcmp(cfg0.demean, 'yes')
    mY = mean(Y, 1);
    Y = Y - repmat(mY, [numN, 1]);
    
    decoder.mY = mY';
end

% Estimate filters
if strcmp(cfg0.returnPattern, 'yes')
    decoder.pattern = zeros(numF, numC);
end
decoder.W = zeros(numC, numF);

for ic = 1:numC
    % Estimate leadfield for current channel
    l = (X(:, ic)'*X(:, ic))\(X(:, ic)'*Y); %returns 1 x 64
    
    if strcmp(cfg0.returnPattern, 'yes')
        decoder.pattern(:, ic) = l';
    end
    
    % Estimate noise
    N = Y - X(:, ic)*l;
    
    % Estimate noise covariance
    S = cov(N);
    
    % Regularize
    S = (1-cfg0.gamma)*S + cfg0.gamma*eye(numF)*trace(S)/numF;
    
    % Calculate filter
    decoder.W(ic, :) = l/S;
    decoder.W(ic, :) = decoder.W(ic, :)/(decoder.W(ic, :)*l');
end

end
