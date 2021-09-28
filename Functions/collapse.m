
function data = collapse (data, factors2collapse, averageTrials)

% COLLAPSE  Collapses factor in a multi-dimensional array of condition-average data or trial data.
%   
%   [data] = COLLAPSE(data,2) collapses second dimension of data
%   [data] = COLLAPSE(data,[2 3]) collapses second and third dimensions of data (in that order)
%
%   inputs
%       data :              n-dimensional numeric matrix, where each dimension represents one factor/condition
%                           OR n-dimensional cell array with each cell representing one cell in design matrix, 
%                               with trials in format elec x time x trial
%
%   edits
%       27/06/2016: edited to maintain non-collapsed singleton dims
%       15/02/2017: added functionality for trial-based data
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   %  Cooper Smout                         %
%   %  c.smout@uq.edu.au                    %
%   %  Queensland Brain Institute           %
%   %  University of Queensland, Australia  %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% settings
if nargin<3
    averageTrials = 0;
end

% orig size
origNdims = ndims(data);

% check no dimensions larger than actual data
if any(factors2collapse>origNdims)
    factors2collapse = factors2collapse(factors2collapse<=origNdims);
    warning('Input dimension/s exceed matrix dimensions')
end

% determine data type
if isnumeric(data) || islogical(data) % condition average data

    % collapse variables
    for FACTOR = factors2collapse
        data = nanmean(data,FACTOR);
    end

    % permute to remove collapsed singletons
    conds = 1:origNdims;
    permuteOrder = find(~ismember(conds,factors2collapse));
    if length(permuteOrder)>1
        for FACTOR = factors2collapse
            if FACTOR~=origNdims % last dimension will collapse automatically
                permuteOrder = [permuteOrder FACTOR];
            end
        end
        data = permute(data,permuteOrder);
    else % only one dimension remaining, put all else at end
        data = permute(data,[permuteOrder factors2collapse]);
    end

elseif iscell(data) % trial data (elec x time x trial)

    warning('CHECK FUNCTION WORKS --- CODED LATE AT NIGHT!!!')

    % check size
    origSize = size(data);
    if origNdims>10
        error('ONLY CODED TO ACCEPT UP TO 10 FACTORS')
    end
    
    % check for singletons
    if ~isequal(size(data),size(squeeze(data)))
        error('NOT CODED FOR SINGLETONS (SEE SQUEEZE AT END)')
    end
    
    % check format
    sz = size(data{1});
    disp([num2str(sz(1)) ' electrode/s found'])
    disp([num2str(sz(2)) ' timepoint/s found'])
    
    % loop collapsed factors
    for FACTOR = factors2collapse 
        
        % initiate new data structure
        newData = {};
        
        % loop dimensions
        for d1 = 1:size(data,1) % condition 1
            for d2 = 1:size(data,2) % condition 2
                for d3 = 1:size(data,3) % condition 3
                    for d4 = 1:size(data,4) % condition 4
                        for d5 = 1:size(data,5) % condition 5
                            for d6 = 1:size(data,6) % condition 6
                                for d7 = 1:size(data,7) % condition 7
                                    for d8 = 1:size(data,8) % condition 8
                                        for d9 = 1:size(data,9) % condition 9
                                            for d10 = 1:size(data,10) % condition 10

                                                % create list of new cell IDs
                                                list = {'d1','d2','d3','d4','d5','d6','d7','d8','d9','d10'};
                                                list{FACTOR} = '1';
                                                keptString = strjoin(list,',');

                                                % create list of old cell IDs
                                                list{FACTOR} = ':';
                                                origString = strjoin(list,',');

                                                % concatenate old cells into new data cell
                                                evalstr = ['newData{' keptString '} = cat(3,data{' origString '});'];
                                                eval(evalstr)
                                                
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
        % update data
        data = newData;
        
    end
    
    % remove collapsed singleton cell dimensions
    data = squeeze(data);
    
    % average across trials
    if averageTrials
        newData = nan([sz(1:2) size(data)]);
        for d1 = 1:size(data,1) % condition 1
            for d2 = 1:size(data,2) % condition 2
                for d3 = 1:size(data,3) % condition 3
                    for d4 = 1:size(data,4) % condition 4
                        for d5 = 1:size(data,5) % condition 5
                            for d6 = 1:size(data,6) % condition 6
                                for d7 = 1:size(data,7) % condition 7
                                    for d8 = 1:size(data,8) % condition 8
                                        for d9 = 1:size(data,9) % condition 9
                                            for d10 = 1:size(data,10) % condition 10
                                                newData(:,:,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10) = mean(data{d1,d2,d3,d4,d5,d6,d7,d8,d9,d10},3);
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
        % rearrange to place electrodes and timepoints at end of array?
        data = newData; %permute(newData,[3:ndims(newData) 1 2]);
        
    end
    
else
    error('UNKNOWN DATA TYPE')
end
end



%% EMBEDDED FUNCTION TO PREVENT FUNCTION CONFLICTS WITH EEGLAB (CLEANLINE PLUGIN)
function output = strjoin(input, separator)
%STRJOIN Concatenate an array into a single string.
%
%     S = strjoin(C)
%     S = strjoin(C, separator)
%
% Description
%
% S = strjoin(C) takes an array C and returns a string S which concatenates
% array elements with comma. C can be a cell array of strings, a character
% array, a numeric array, or a logical array. If C is a matrix, it is first
% flattened to get an array and concateneted. S = strjoin(C, separator) also
% specifies separator for string concatenation. The default separator is comma.
%
% Examples
%
%     >> str = strjoin({'this','is','a','cell','array'})
%     str =
%     this,is,a,cell,array
%
%     >> str = strjoin([1,2,2],'_')
%     str =
%     1_2_2
%
%     >> str = strjoin({1,2,2,'string'},'\t')
%     str =
%     1 2 2 string
%

  if nargin < 2, separator = ','; end
  assert(ischar(separator), 'Invalid separator input: %s', class(separator));
  separator = strrep(separator, '%', '%%');

  output = '';
  if ~isempty(input)
    if ischar(input)
      input = cellstr(input);
    end
    if isnumeric(input) || islogical(input)
      output = [repmat(sprintf(['%.15g', separator], input(1:end-1)), ...
                       1, ~isscalar(input)), ...
                sprintf('%.15g', input(end))];
    elseif iscellstr(input)
      output = [repmat(sprintf(['%s', separator], input{1:end-1}), ...
                       1, ~isscalar(input)), ...
                sprintf('%s', input{end})];
    elseif iscell(input)
      output = strjoin(cellfun(@(x)strjoin(x, separator), input, ...
                               'UniformOutput', false), ...
                       separator);
    else
      error('strjoin:invalidInput', 'Unsupported input: %s', class(input));
    end
  end
end

