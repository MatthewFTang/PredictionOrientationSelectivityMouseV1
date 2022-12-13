function design = designMatrix_dummy(cfg0, group)
% [design] = designMatrix_dummy(cfg, G)
%    Returns a design matrix using dummy coding for the groups.
%
%    G           Array of length N, where N is the number of trials, that contains the group each trial
%                belongs to, as identified by a unique number.
%    cfg         Configuration struct that can possess the following fields:
%                .intercept = 'yes' or 'no'             Whether to include an intercept or not.
%                                                       Default = 'no';
%
%    design      The design matrix of size G x N, where G is the number of G, containing the dummy codes.
%
%    See also designMatrix_BH designMatrix_discriminant

%    Created by Pim Mostert, 2016
if ~isfield(cfg0, 'intercept')
    cfg0.intercept = 'no';
end

group = group(:);

GROUP = unique(group);

numG = length(GROUP);
numN = length(group);

design = zeros(numG, numN);

for iGroup = 1:numG
    design(iGroup, group==GROUP(iGroup)) = 1;
end

% Include intercept
if strcmp(cfg0.intercept, 'yes')
    design(1, :) = 1;
end

end


