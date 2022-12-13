function [y]= fitCic180(params,x)


% y = cg(x,params(1),params(2),params(3),params(4),params(5),180);
% (x,b,a,po,tw,c,P) 
b=params(1);
a=params(2);
po=params(3);
tw=params(4);
P=180;

y=zeros(size(x));
for i=1:length(x(:))
    y(i)=b;
    for j=-4:4
        y(i)=y(i)+a*exp(-(x(i)-po-j*P)^2 / (2*tw^2));
    end
end
