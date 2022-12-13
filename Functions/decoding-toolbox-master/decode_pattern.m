function Xhat = decode_pattern(cfg0, decoder, Y)
% [Xhat] = decode_pattern(cfg, decoder, Y)
%    Estimate the activity of latent components using a linear decoder, obtained from an
%    appropriate training function. Several components may be estimated independently.
%
%    decoder     The linear decoder obtained from e.g. train_pattern.
%    Y           Matrix of size F x N, where F is the number of features and the N the number of trials,
%                that contains the data that is to be decoded.
%    cfg         Configuration struct that can possess the following fields:
%                .demean                          Whether the data should be demeaned (per feature,
%                                                 over trials) prior to decoding. The mean can be
%                                                 specified in the following ways:
%                        = 'trainData'            The nanmean of the training data (default).
%                        = 'testData'             The mean of the testing data.
%                        = [F x 1] vector         Manually specified mean.
%                        = 'no'                   No demeaning.
%                .covariance                      Whether the pattern should be multiplied by the
%                                                 inverse of a covariance matrix.
%                        = 'testData'             The nancov of the testing data. Specify
%                                                 cfg.gamma for regularization.
%                        = [F x F] vector         Manually specified covariance matrix, where F is
%                                                 the number of features (e.g. sensors).
%                        = 'no'                   (default)
%
%    Xhat        Vector or matrix of size C x N, where C is the number of components, containing
%                the decoded data.
%
%    See also TRAIN_PATTERN.

%    Created by Pim Mostert, 2017

if ~isfield(cfg0, 'demean')
    cfg0.demean = 'trainData';
end
if ~isfield(cfg0, 'covariance')
    cfg0.covariance = 'no';
end

% Convert to matrix
numN = size(Y, 2);
numF = size(Y, 1);

% Demean
if strcmp(cfg0.demean, 'trainData')
    if ~isfield(decoder, 'mY')
        error('No mean found in decoder');
    end
    
    Y = Y - repmat(decoder.mY, [1, numN]);
elseif strcmp(cfg0.demean, 'testData')
    m = nanmean(Y, 2);
    
    Y = Y - repmat(m, [1, numN]);
elseif isnumeric(cfg0.demean)
    Y = Y - repmat(cfg0.demean, [1, numN]);
else
    error('Demeaning configuration ''%s'' is unknown', cfg0.demean);
end

% Calculate filter
if strcmp(cfg0.covariance, 'testData')
    S = nancov(Y');
    
    % Regularize
    if isfield(cfg0, 'gamma')
        S = (1-cfg0.gamma)*S + cfg0.gamma*eye(numF)*trace(S)/numF;
    end
    
    % Calculate filter
    W = decoder.W/S;
elseif strcmp(cfg0.covariance, 'no');
    W = decoder.W;
elseif isnumeric(cfg0.covariance)
    W = decoder.W/cfg0.covariance;
end

% Decode
Xhat = W*Y;

end
