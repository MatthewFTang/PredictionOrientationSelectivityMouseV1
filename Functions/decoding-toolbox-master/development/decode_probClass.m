function pPost = decode_probClass(cfg0, decoder, Y)
% pPost = decode_beamformer(cfg, decoder, Y)
%    Apply probabilistic multiclass decoder, obtained from an appropriate training function. It 
%    returns the posterior probabilities of a the trial belonging to a class, given the data. At 
%    this moment the inclusion of (non-uniform) prior probabilities is however not yet implemented.
%
%    decoder     The probabilistic multiclass decoder obtained from e.g. train_probClass.
%    Y           Matrix of size F x N, where F is the number of features and the N the number of trials,
%                that contains the data that is to be classified.
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
%    pPost       Matrix of size C x N, where C is the number of classes, that contains the posterior
%                probabilities of the trial belonging to each of the classes, given the data.
%
%    See also TRAIN_PROBCLASS.

%    Created by Pim Mostert, 2017

if ~isfield(cfg0, 'demean')
    cfg0.demean = 'trainData';
end

% Useful variables
numN = size(Y, 2);
numC = size(decoder.W, 1);

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
pPost = exp(decoder.W*Y + repmat(decoder.b, [1, numN]));
pPost = pPost ./ repmat(sum(pPost, 1), [numC, 1]);

end
