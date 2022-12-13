function gamma = shrinkageGamma(X, memEff, feedback)
% [gamma] = shrinkageGamma(X[, memEff, feedback])
%    Returns the optimal shrinkage parameter for data X,
%    as described in Blankertz et al. (2011) Neuroimage.
%
%    X           Matrix of size F x N, where F is the number of features and N the number of trials,
%                that contains the data for which the optimal shrinkage parameter is to be calculated.
%
%    memEff      Optional argument. Specificy memEff=1 to reduce the memory requirement of the function,
%                for instance when X is large. If X is small and/or enough memory is available, then
%                memEff=0 is probably faster, although I haven't thoroughly verified this. Default = 0.
%
%    feedback    Optional argument, 'yes' or 'no'. Whether or not to provide feedback on the progress
%                of the calculation. This only has an effect when memEff=1. Default = 'yes'.
%
%    gamma       The calculated optimal shrinkage parameter that may be passed on to an appropriate
%                training function.
%
%    See also TRAIN_LDA and TRAIN_BEAMFORMER.

%    Created by Pim Mostert, 2016

if nargin < 2
    memEff = 0;

end
if nargin < 3
    feedback = 'yes';
end

numF = size(X, 1);
numN = size(X, 2);

if ~memEff
    m = mean(X, 2);
    S = cov(X');

    nu = trace(S)/numF;

    z = zeros(numF, numF, numN);
    for n = 1:numN
        z(:, :, n) = (X(:, n) - m)*(X(:, n) - m)';
    end

    gamma = ...
        (numN/((numN-1)^2)) * ...
        sum(sum(var(z, [], 3))) / ...
        sum(sum((S - nu*eye(numF)).^2)) ...
    ;
else
    % Demean X
    X = X - mean(X, 2)*ones(1, numN);

    sumVarDiag = 0;

    % Do diagonal elements first to determine nu
    diagS = zeros(numF, 1);
    for iF = 1:numF
        s = X(iF, :).^2;

        diagS(iF) = sum(s)/(numN-1);
        s = s - mean(s);
        sumVarDiag = sumVarDiag + (s*s')/(numN-1);
    end

    nu = mean(diagS);

    diagS = diagS - nu;
    sumSdiag = diagS'*diagS;

    % Do off-diagonal elements, but only one triangle
    sumS = 0;
    sumVar = 0;

    tStart = tic;
    tCur = tStart;
    for iF1 = 2:numF
        for iF2 = 1:(iF1-1)
            s = X(iF1, :) .* X(iF2, :);

            sumS = sumS + (sum(s)/(numN-1)).^2;
            s = s - mean(s);
            sumVar = sumVar + (s*s')/(numN-1); 
            
            if (toc(tCur) > 2) && strcmp(feedback, 'yes')
                pDone = (((iF1-2)*(iF1-1)/2) + iF2) / ((numF-1)*numF/2);
                fprintf('%s - finished calculating off-diagonal elements %.2f%% - time remaining: %.2f s\n', mfilename, pDone*100, toc(tStart) * ((1-pDone)/pDone))
                tCur = tic;
            end
        end
    end
    sumS = sumS*2 + sumSdiag;
    sumVar = sumVar*2 + sumVarDiag;

    gamma = (numN/((numN-1)^2))*sumVar/sumS;
end
            
end            
            
            
            
           
        
        
        
        
        