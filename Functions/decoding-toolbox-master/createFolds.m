function folds = createFolds(cfg0, G)
% [folds] = createFolds(cfg, G)
%    Divides trials into folds for cross-validation. The function attempts 
%    to divide the trials as evenly as possible over the number of folds, 
%    while trying to balance the different conditions as much as possible.
%
%    cfg         Configuration struct that can possess the following fields:
%                .nFold = [scalar]                Number of folds the trials should be divided over.
%
%    G           Vector of length N, where N is the number of trials, that specifies
%                the condition to which that trial belongs, as identified by a unique number.
%
%    folds       A cell-array of length cfg.nFold that contains in each cell the indices
%                of the trials belonging to that particular fold.

%    Created by Pim Mostert, 2016

G = G(:);

CONDS = unique(G);
N_conds = length(CONDS);

folds = cell(cfg0.nFold, 1);
for iCond = 1:N_conds
    % Find indices
    index = find(G == CONDS(iCond));
    nIndex = length(index);
    
    % Shuffle
    index = index(randperm(nIndex));

    % Distribute across folds
    groupNumber = floor((0:(nIndex-1))*(cfg0.nFold/nIndex))+1;
    
    for iFold = 1:cfg0.nFold
        folds{iFold} = [folds{iFold}, index(groupNumber==iFold)'];
    end
end

end