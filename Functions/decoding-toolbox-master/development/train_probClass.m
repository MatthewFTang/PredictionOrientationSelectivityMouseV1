function decoder = train_probClass(cfg0, X, Y)
% [decoder] = train_probClass(cfg, X, Y)
%    Trains a multiclass probabilistic classifier. It (if I'm not mistaking) is identical to logistic
%    regression under the assumptions that 1) the data within each class are Gaussian distributed
%    and 2) the covariances are identical across classes.
%
%    X           Array of length N, where N is the number of trials, that contains the class
%                labels. Must be numeric, though can have arbitrary values, and can contain
%                multiple classes.
%    Y           Matrix of size F x N, where F is the number of features, that contains the
%                training data.
%    cfg         Configuration struct that can possess the following fields:
%                .gamma = [scalar]                Shrinkage regularization parameter, with range [0 1]. 
%                                                 No default given.
%                .discardNan = 'yes' or 'no'      Whether trials with NaN in either X or Y should be
%                                                 removed prior to training. Default is 'no'.
%                .demean = 'yes' or 'no'          Whether to demean the data first (per feature, over
%                                                 trials). Default = 'yes';.
%
%    decoder     The trained decoder, that may be passed to an appropriate decoding function.
%
%    See also DECODE_PROBCLASS

%    Created by Pim Mostert, 2017

decoder = [];

%% Pre-process cfg-struct
if ~isfield(cfg0, 'discardNan')
    cfg0.discardNan = 'no';
end
if ~isfield(cfg0, 'gamma')
    error(sprintf('No regularization (cfg.gamma) specified!\nIf this is intended, then please specifyc cfg0.gamma = 0'));
end
if ~isfield(cfg0, 'demean')
    cfg0.demean = 'yes';
end    

%% Pre-process data
X = X(:);
Y = Y';

if strcmp(cfg0.discardNan, 'yes')
    iNan = isnan(X) | any(isnan(Y), 2);
    
    X = X(~iNan);
    Y = Y(~iNan, :);
end

numN = size(Y, 1);
numF = size(Y, 2);

% Demean data
if strcmp(cfg0.demean, 'yes')
    mY = mean(Y, 1);
    Y = Y - repmat(mY, [numN, 1]);
    
    decoder.mY = mY';
end

% Extract class labels
decoder.classLabels = unique(X);
numC = length(decoder.classLabels);

%% Calculate decoder
% Calculate means and shared covariance matrix
m = zeros(numF, numC);
S = zeros(numF, numF);

for ic = 1:numC
    m(:, ic) = mean(Y(X==decoder.classLabels(ic), :), 1)';
    
    S = S + cov(Y(X==decoder.classLabels(ic), :));
end
S = S/numC;

% Regularize
S = (1-cfg0.gamma)*S + cfg0.gamma*eye(numF)*trace(S)/numF;

% Calculate weights and bias
decoder.W = m'/S;
decoder.b = -0.5*sum(m' .* decoder.W, 2);

end
