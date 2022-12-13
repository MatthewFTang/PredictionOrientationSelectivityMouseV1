function decoder = train_LDA(cfg0, X, Y)
% [decoder] = train_LDA(cfg, X, Y)
%    Trains a linear discriminant analysis decoder.
%
%    X           Array of length N, where N is the number of trials, that specifies the class
%                label (0 or 1) for each trial.
%    Y           Matrix of size F x N, where F is the number of features, that contains the
%                training data.
%    cfg         Configuration struct that can possess the following fields:
%                .gamma = [scalar]                Shrinkage regularization parameter, with range [0 1]. 
%                                                 No default given.
%                .discardNan = 'yes' or 'no'      Whether trials with NaN in either X or Y should be
%                                                 removed prior to training. Default is 'no'.
%
%    decoder     The trained decoder, that may be passed to an appropriate decoding function.
%
%    See also DECODE_LDA

%    Created by Pim Mostert, 2016

decoder = [];

%% Pre-process cfg-struct
if ~isfield(cfg0, 'discardNan')
    cfg0.discardNan = 'no';
end
if ~isfield(cfg0, 'gamma')
    warning(sprintf('No regularization (cfg.gamma) specified!\nIf this is intended, then please specifyc cfg0.gamma = 0'));
end

%% Pre-process data
X = X(:);
Y = Y';

if strcmp(cfg0.discardNan, 'yes')
    iNan = isnan(X) | any(isnan(Y), 2);
    
    X = X(~iNan);
    Y = Y(~iNan, :);
end

numF = size(Y, 2);

%% Calculate decoder
% Calculate group means
m0 = mean(Y(X==0, :), 1);
m1 = mean(Y(X==1, :), 1);

% Mean difference
d = m1 - m0;

% Calculate covariances
S0 = cov(Y(X==0, :), 1);
S1 = cov(Y(X==1, :), 1);

% Pool covariance
S = 0.5*(S0 + S1);

% Regularize pooled covariance matrix
if isfield(cfg0, 'gamma')
    S = (1-cfg0.gamma)*S + cfg0.gamma*eye(numF)*trace(S)/numF;
end

% Calculate weights
W = d/S;

% Normalize weights
W = W/(W*d');

%% Include overall (unweighted) mean for demeaning during decoding
decoder.mY = 0.5*(m0 + m1);

%% Return
decoder.W = W;

end
