function [y] = filtfast(x,dim,i,ftype,fsize)

%author: V Wyart
%edits: N Myers (added Butterworth filter)
if nargin < 5
    error('Missing input arguments!');
end

if dim > ndims(x)
    error('Filter dimension exceeds input dimensionality!');
end
if isempty(i)
    i = 1:size(x,dim);
end

fungaussian = @(x,m,s)exp(-0.5*((x-m)./s).^2);
switch ftype
    case 'boxcar'
        fsize = floor(fsize/2)*2+1;
        w = ones(1,fsize);
        w = w/sum(w);
        b = 1;
    case 'gaussian'
        w = fungaussian(-ceil(2*fsize):+ceil(2*fsize),0,fsize);
        w = w/sum(w);
        fsize = length(w);
        b = 1;
    case 'butter'
        %Wn    = [4 8]/data.fsample*2;
        [b w] = butter(1,fsize);
    otherwise
        error('Unknown filter type!');
end

dims = [dim,setdiff(1:ndims(x),dim)];

y = permute(x,dims);
y = filter(w,b,y,[],1); ysize = size(y);
if strcmp(ftype,'butter')
    y = y(min(i+floor(fsize(1)*250/2),ysize(1)),:);
else
    y = y(min(i+floor(fsize/2),ysize(1)),:);
end
y = reshape(y,[length(i),ysize(2:end)]);
y = ipermute(y,dims);

end