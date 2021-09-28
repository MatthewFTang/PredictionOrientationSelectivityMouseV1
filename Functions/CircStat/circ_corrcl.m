function [rho pval] = circ_corrcl(alpha, x)
%
% [rho pval ts] = circ_corrcc(alpha, x)
%   Correlation coefficient between one circular and one or more linear random
%   variables (1 column per variable).
%
%   Input:
%     alpha   sample of angles in radians
%     x       sample of linear random variable
%
%   Output:
%     rho     correlation coefficient
%     pval    p-value
%
% References:
%     Biostatistical Analysis, J. H. Zar, p. 651
%
% PHB 6/7/2008
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html


if size(alpha,2) > size(alpha,1)
	alpha = alpha';
end

if size(alpha,2) > size(alpha,1)
	x = x';
end

if size(alpha,1)~=size(x,1)
  error('Input dimensions do not match.')
end

n = size(alpha,1);

% compute correlation coefficent for sin and cos independently

usefast = true;
if usefast
    %disp('Fast!')
    rxs = corrcoef([x,sin(alpha)]); rxs = rxs(1:(end-1),end);
    rxc = corrcoef([x,cos(alpha)]); rxc = rxc(1:(end-1),end);
    rcs = corrcoef([sin(alpha),cos(alpha)]); rcs = rcs(1,2);
else
    rxs = corr(x,sin(alpha));
    rxc = corr(x,cos(alpha));
    rcs = corr(sin(alpha),cos(alpha));
end
% compute angular-linear correlation (equ. 27.47)
rho = sqrt((rxc.^2 + rxs.^2 - 2*rxc.*rxs*rcs)/(1-rcs^2));

% compute pvalue
pval = 1 - chi2cdf(n*rho.^2,2);

