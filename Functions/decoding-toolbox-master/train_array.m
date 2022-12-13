function decoder = train_array(cfg0, X, Y)
% [decoders] = train_array(cfg, X, Y)
%   HELP NEEDS TO BE UPDATED
%    Trains an array of decoders along a specified dimension in the data, e.g. time.
%
%    X           Cell-vector of length N, where N is the number of trials, that contains the training
%                information (e.g. labels). 
%
%    Y           Cell-vector of length N that contains the training data. Each cell should now contain
%                not just a vector containing the features, but instead a (higher-dimensional) matrix 
%                that also contains an additional dimension, e.g. [sensors x time].
%
%    cfg         Configuration struct that can possess the fields below. It may also be a struct array,
%                in which case the corresponding element is passed within each iteration across the
%                specified dimension.
%                .trainfun = [function_name]      The training function that is used for training each decoder.
%                .traincfg = [struct]             The configuration struct that will be passed on to
%                                                 the training function. Default = [];
%                .dim = [scalar]                  Dimension of the data in each cell along which the array
%                                                 of decoders is constructed. Default = 2.
%                .feedback = 'yes' or 'no'        Whether the function should provide feedback on its progress.
%                                                 Default = 'no'.            
%
%    decoders    The array of trained decoders, that may be passed to an appropriate decoding function, e.g.
%                decode_multiple or decode_multipleGeneralization.
%
%    See also DECODE_ARRAY and DECODE_ARRAYGENERALIZATION

%    Created by Pim Mostert, 2016

if ~isfield(cfg0, 'traincfg')
    cfg0.traincfg = [];
end
if ~isfield(cfg0, 'feedback');
    cfg0.feedback = 'no';
end

dims = size(Y);

numN = dims(end);
numDec = dims(end-1);

dimsSub = dims(1:(end-2));

%% Reshape data to allow for arbitrary dimensionality
Y = reshape(Y, [prod(dimsSub), numDec, numN]);

%% Iterate over decoders
tDec = tic;

decoder = cell(numDec, 1);
for iDec = 1:numDec
    % Select data for current decoder
    curY = reshape(Y(:, iDec, :), [dimsSub, numN]);

    % Train decoder
    if (length(cfg0.traincfg) > 1)
        curTraincfg = cfg0.traincfg(iDec);
    else
        curTraincfg = cfg0.traincfg;
    end
    
    decoder{iDec} = feval(cfg0.trainfun, curTraincfg, X, curY);

    % Feedback    
    if (toc(tDec) > 2) && strcmp(cfg0.feedback, 'yes')
        fprintf('%s: finished training decoder %g/%g\n', mfilename, iDec, numDec);
        tDec = tic;
    end    
end

end

