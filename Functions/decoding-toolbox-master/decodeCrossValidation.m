function Xhat = decodeCrossValidation(cfg0, X, Y)
% [Xhat] = decodeCrossValidation(cfg, X, Y)
%    Implements k-fold cross-validation, in which a subset of trials is left out
%    in each iteration as testing data, while training on the remaining trials.
%
%    cfg         Configuration struct that can possess the following fields:
%                .trainfun = [function_name]      The training function that is used for training.
%                .traincfg = [struct]             The configuration struct that will be passed on to
%                                                 the training function. Default = [];
%                .decodefun = [function_name]     The decoding function that is used for decoding.
%                .decodecfg = [struct]            The configuration struct that will be passed on to
%                                                 the decoding function. Default = [].
%                .folds = [cell_array]            A cell-array of length k, where k is the number of folds,
%                                                 in which each cell contains a vector with the trial numbers
%                                                 belonging to that particular fold.
%                .feedback = 'yes' or 'no'        Whether the function should provide feedback on its progress.
%                                                 Default = 'no'.            
%
%    X           Matrix of arbitrary dimensions, but of which the last dimension is N, that contains
%                the training information. In each fold, a selection of this matrix (along the last
%                dimension) is sent to the training function.
%
%    Y           Matrix of arbitrary dimensions, but of which the last dimension corresponds to the
%                number of trials N, that contains the data. In each fold, a selection of this matrix
%                (along the last dimension) is sent to the training and decoding function.
%
%    Xhat        Matrix of dimensions as output by the decoding functiong, plus an additional dimension
%                of length N, that contains the decoded data.
%
%    See also CREATEFOLDS

%    Created by Pim Mostert, 2016

tStart = tic;

if ~isfield(cfg0, 'traincfg')
    cfg0.traincfg = [];
end
if ~isfield(cfg0, 'decodecfg')
    cfg0.decodecfg = [];
end
if ~isfield(cfg0, 'feedback')
    cfg0.feedback = 'no';
end

dimsY = size(Y);

numN = dimsY(end);
numFold = length(cfg0.folds);

%% Reshape data to allow for arbitrary dimensionality
Y = reshape(Y, [prod(dimsY(1:(end-1))), numN]);
if isvector(X)
    X = X(:)';
    dimsX = size(X);
else
    dimsX = size(X);
    X = reshape(X, [prod(dimsX(1:(end-1))), numN]);
end

%% Do first fold manually, to determine output size of decoder
iFold = 1;

tFold = tic;

index_train = cell2mat(cfg0.folds((1:numFold) ~= iFold)');
index_decode = cfg0.folds{iFold};

% Select training data
Y_train = reshape(Y(:, index_train), [dimsY(1:(end-1)), length(index_train)]);
X_train = reshape(X(:, index_train), [dimsX(1:(end-1)), length(index_train)]);

% Train decoder
decoder = feval(cfg0.trainfun, cfg0.traincfg, X_train, Y_train);

% Select data to be decoded
Y_decode = reshape(Y(:, index_decode), [dimsY(1:(end-1)), length(index_decode)]);

% Decode data
Xhat_curFold = feval(cfg0.decodefun, cfg0.decodecfg, decoder, Y_decode);

% Feedback
if strcmp(cfg0.feedback, 'yes')
    fprintf('%s: finished fold %g/%g - it took %.2f s\n', mfilename, iFold, numFold, toc(tFold));
end

%% Allocate memory for results and do rest of folds
dimsOut = size(Xhat_curFold);
dimsOut = dimsOut(1:(end-1));

Xhat = zeros([prod(dimsOut), numN]);
Xhat(:, index_decode) = reshape(Xhat_curFold, [prod(dimsOut), length(index_decode)]);

for iFold = 2:numFold
    tFold = tic;

    index_train = cell2mat(cfg0.folds((1:numFold) ~= iFold)');
    index_decode = cfg0.folds{iFold};

    % Select training data
    Y_train = reshape(Y(:, index_train), [dimsY(1:(end-1)), length(index_train)]);
    X_train = reshape(X(:, index_train), [dimsX(1:(end-1)), length(index_train)]);

    % Train decoder
    decoder = feval(cfg0.trainfun, cfg0.traincfg, X_train, Y_train);

    % Select data to be decoded
    Y_decode = reshape(Y(:, index_decode), [dimsY(1:(end-1)), length(index_decode)]);

    % Decode data
    Xhat_curFold = feval(cfg0.decodefun, cfg0.decodecfg, decoder, Y_decode);
    Xhat(:, index_decode) = reshape(Xhat_curFold, [prod(dimsOut), length(index_decode)]);

    % Feedback
    if strcmp(cfg0.feedback, 'yes')
        fprintf('%s: finished fold %g/%g - it took %.2f s\n', mfilename, iFold, numFold, toc(tFold));
    end
end

%% Return
Xhat = reshape(Xhat, [dimsOut, numN]);

if strcmp(cfg0.feedback, 'yes')
    fprintf('%s - all finished - it took %.2f s\n', mfilename, toc(tStart));
end

end

    
