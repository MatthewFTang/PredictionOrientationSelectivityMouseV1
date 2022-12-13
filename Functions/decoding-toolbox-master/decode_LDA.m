function Xhat = decode_LDA(cfg0, decoder, Y)
% [Xhat] = decode_LDA(cfg, decoder, Y)
%    Decodes data using a linear discriminant analysis decoder as obtained from
%    an appropriate training function, e.g. train_LDA.
%
%    decoder     The decoder as obtained from an appropriate training function.
%
%    Y           Matrix of size F x N, where N is the number of trials and F the number of
%                features, that contains the data to be decoded.
%
%    cfg         Configuration struct that can possess the following fields:
%                .demean = 'yes' or 'no'          Whether the data should be demeaned (per feature,
%                                                 over trials) prior to decoding, based on the mean
%                                                 of the training data. Default = 'yes'.
%
%    Xhat        Vector of length N that contains the decoded data.
%
%    See also TRAIN_LDA

%    Created by Pim Mostert, 2016

%% Pre-process cfg-struct
if ~isfield(cfg0, 'demean')
    cfg0.demean = 'yes';
end

%% Pre-process data
numN = size(Y, 2);

% Demean
if strcmp(cfg0.demean, 'yes')
    Y = Y - repmat(decoder.mY', [1, numN]);
end

%% Decode
Xhat = decoder.W*Y;

end
