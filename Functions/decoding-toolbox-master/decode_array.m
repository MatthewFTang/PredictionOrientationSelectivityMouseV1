function Xhat = decode_array(cfg0, decoder, Y)
% [Xhat] = decode_array(cfg, decoders, Y)
%    Decodes the data along a specified dimension, e.g. time, using an array of decoders.
%
%    decoders    Cell-vector of length D, that contains an array of decoders, where D is the number
%                of decoders. This array may be obtained from an appropriate training function 
%                such as train_array.
%
%    Y           Matrix of arbitrary dimension that contains the data to be decoded, though the last 
%                two dimensions should be D and N, respectively, where N is the number of trials.
%                For example, [sensors x time x trials] if the array of decoders was trained along time.
%
%    cfg         Configuration struct that can possess the following fields:
%                .decodefun = [function_name]     The decoding function to which each of the decoders 
%                                                 is passed on.
%                .decodecfg = [struct]            The configuration struct that will be passed on to
%                                                 the decoding function. Default = [];
%                .feedback = 'yes' or 'no'        Whether the function should provide feedback on its progress.
%                                                 Default = 'no'.            
%
%    Xhat        Matrix of dimensions as output by the decoding functiong, plus the additional dimensions
%                of D and N, that contains the decoded data for each decoder and trial.
%
%    See also TRAIN_ARRAY and DECODE_ARRAYGENERALIZATION

%    Created by Pim Mostert, 2016

tStart = tic;

if ~isfield(cfg0, 'decodecfg')
    cfg0.decodecfg = [];
end
if ~isfield(cfg0, 'feedback')
    cfg0.feedback = 'no';
end    

dims = size(Y);

numN = dims(end);
numDec = dims(end-1);

dimsSub = dims(1:(end-2));

%% Reshape data to allow for arbitrary dimensionality
Y = reshape(Y, [prod(dimsSub), numDec, numN]);

%% Do first decoder manually, to obtain output size of decoder
Xhat_curDec = feval(cfg0.decodefun, cfg0.decodecfg, decoder{1}, squeeze(Y(:, 1, :)));
dimsOut = size(Xhat_curDec);
dimsOut = dimsOut(1:(end-1));

%% Allocate memory for output and iterate over remaining decoders
Xhat = zeros([prod(dimsOut), numDec, numN]);
Xhat(:, 1, :) = reshape(Xhat_curDec, [prod(dimsOut), numN]);

tDec = tic;
for iDec = 2:numDec
    Xhat_curDec = feval(cfg0.decodefun, cfg0.decodecfg, decoder{iDec}, squeeze(Y(:, iDec, :)));
    Xhat(:, iDec, :) = reshape(Xhat_curDec, [prod(dimsOut), numN]);

    if (toc(tDec) > 2) && strcmp(cfg0.feedback, 'yes')
        fprintf('%s - finished decoder %g/%g\n', mfilename, iDec, numDec);
        tDec = tic;
    end
end    
    
%% Reshape output
Xhat = reshape(Xhat, [dimsOut, numDec, numN]);

if strcmp(cfg0.feedback, 'yes')
    fprintf('%s - all finished - it took %.2f s\n', mfilename, toc(tStart));
end

end        
        
    
    
    
    
    