function Xhat = decode_arrayGeneralization(cfg0, decoder, Y)
% [Xhat] = train_arrayGeneralization(cfg, decoders, Y)
%    Decodes the data along a specified dimension, e.g. time, using an array of decoders. Additionally,
%    each decoder is also applied to all other points along this dimension. For example, this function
%    can implement the temporal generalization method (King & Dehaene, 2014).
%
%    decoders    Cell-vector of length D, that contains an array of decoders, where D is the number
%                of decoders. This array may be obtained from an appropriate training function 
%                such as train_array.
%
%    Y           Matrix of arbitrary dimension that contains the data to be decoded, though the last 
%                two dimensions should be T and N, respectively, where T is the dimension along which all
%                decoders should be applied iteratively and N is the number of trials. For example, 
%                [sensors x time x trials] if the array of decoders should be applied over time. Note that
%                although D and T are likely to correspond to the same quantity, e.g. time, they do not
%                necessarily have the same size, for instance when training on one task and generalizing
%                to another task. If Y is 2-dimensional, then it is assumed to correspond to one single
%                trial, i.e. [sensors x time x 1].
%
%    cfg         Configuration struct that can possess the following fields:
%                .decodefun = [function_name]     The decoding function to which each of the decoders 
%                                                 is passed on.
%                .decodecfg = [struct]            The configuration struct that will be passed on to
%                                                 the decoding function. Default = [];
%                .feedback = 'yes' or 'no'        Whether the function should provide feedback on its progress.
%                                                 Default = 'no'.            
%                .oneTrial = 'yes' or 'no'        <yet to be described>
%
%    Xhat        Matrix of dimensions as output by the decoding functiong, plus the additional dimensions
%                of D, D and N, that contains the decoded data for each decoder, applied to all other
%                points along the dimension of interast, for each trial.
%
%    See also TRAIN_ARRAY and DECODE_ARRAY

%    Created by Pim Mostert, 2016

tStart = tic;

if ~isfield(cfg0, 'decodecfg')
    cfg0.decodecfg = [];
end
if ~isfield(cfg0, 'feedback')
    cfg0.feedback = 'no';
end    
if ~isfield(cfg0, 'oneTrial')
    cfg0.oneTrial = 'no';
end

dims = size(Y);

if (length(size(Y)) == 2)
    dims = [dims 1];
end

numN = dims(end);
numD = length(decoder);
numT = dims(end-1);

dimsSub = dims(1:(end-2));

%% Reshape data to allow for arbitrary dimensionality
Y = reshape(Y, [prod(dimsSub), numT, numN]);

%% Do first decoder manually, to obtain output size of decoder
Xhat_curDec = feval(cfg0.decodefun, cfg0.decodecfg, decoder{1}, reshape(Y, [dimsSub, numT*numN]));
dimsOut = size(Xhat_curDec);
dimsOut = dimsOut(1:(end-1));

%% Allocate memory for output and iterate over remaining decoders
Xhat = zeros([prod(dimsOut), numD, numT, numN]);
Xhat(:, 1, :, :) = reshape(Xhat_curDec, [prod(dimsOut), numT, numN]);

tDec = tic;
for iDec = 2:numD
    Xhat_curDec = nan(prod(dimsOut), numT, numN);
    for it = 1:numT
        Xhat_curDec(:, it, :) = feval(cfg0.decodefun, cfg0.decodecfg, decoder{iDec}, squeeze(Y(:, it, :)));
    end
    Xhat(:, iDec, :, :) = reshape(Xhat_curDec, [prod(dimsOut), numT, numN]);

    if (toc(tDec) > 2) && strcmp(cfg0.feedback, 'yes')
        fprintf('%s - finished decoder %g/%g\n', mfilename, iDec, numD);
        tDec = tic;
    end
end    
    
%% Reshape output
Xhat = reshape(Xhat, [dimsOut, numD, numT, numN]);

if strcmp(cfg0.feedback, 'yes')
    fprintf('%s - all finished - it took %.2f s\n', mfilename, toc(tStart));
end

end        
        
    
    
    
    
    