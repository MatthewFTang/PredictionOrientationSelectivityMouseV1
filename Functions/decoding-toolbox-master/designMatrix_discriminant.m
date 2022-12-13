function design = designMatrix_discriminant(cfg0, group)
% [design] = designMatrix_dummy(cfg, G)
%    Returns a design matrix using discriminant coding for two groups.
%
%    G           Array of length N, where N is the number of trials, that contains the group each trial
%                belongs to, as identified by a unique number. G must consist entirely of two unique values.
%    cfg         Configuration struct that can possess the following fields:
%                .intercept = 'yes' or 'no'             Whether to include an intercept or not.
%                                                       Default = 'no';
%
%    design      The design matrix of size G x N, where G is the number of G, containing the discriminant codes.
%
%    See also designMatrix_BH designMatrix_dummy

%    Created by Pim Mostert, 2016
if ~isfield(cfg0, 'intercept')
    cfg0.intercept = 'no';
end

group = group(:);
[GROUP, ~, design] = unique(group);

if length(GROUP) ~= 2
    error('G must contain exaclty two unique values');
end

design = (design - 1)*2 - 1;

% Include intercept
if strcmp(cfg0.intercept, 'yes')
    design = [ones(length(group), 1), design];
end

design = design';

end


