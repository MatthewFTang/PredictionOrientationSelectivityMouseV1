function Xhat = decode_beamformer(cfg0, decoder, Y)
% [Xhat] = decode_beamformer(cfg, decoder, Y)
%    Estimate the activity of latent components using a linear decoder, obtained from an
%    appropriate training function. Several components may be estimated independently.
%
%    decoder     The linear decoder obtained from e.g. train_beamformer.
%    Y           Matrix of size F x N, where F is the number of features and the N the number of trials,
%                that contains the data that is to be decoded.
%    cfg         Configuration struct that can possess the following fields:
%                .demean                          Whether the data should be demeaned (per feature,
%                                                 over trials) prior to decoding. The mean can be
%                                                 specified in the following ways:
%                        = 'trainData'            The mean of the training data (default).
%                        = 'testData'             The mean of the testing data.
%                        = [F x 1] vector         Manually specified mean, where F is the number of
%                                                 features (e.g. sensors).
%                        = 'no'                   No demeaning.
%
%    Xhat        Vector or matrix of size C x N, where C is the number of components, containing
%                the decoded data.
%
%    See also TRAIN_BEAMFORMER.

%    Created by Pim Mostert, 2016

if ~isfield(cfg0, 'demean')
    cfg0.demean = 'trainData';
end

% Convert to matrix
numN = size(Y, 2);

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
elseif strcmp(cfg0.demean, 'no')
else
    error('Demeaning configuration ''%s'' is unknown', cfg0.demean);
end
        
% Decode
Xhat = decoder.W*Y;

end
