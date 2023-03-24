%ACF.M
%y = acf(x) computes the first-order autocorrelation of the vector x
%y = acf(x,J) is a J-by-1 vector containing the autocorrelations of x up to order J. 
%Element j of y is the autocorrelation of order j of x.
function y = acf(x,J)

if nargin<2; J=1; end

for j=1:J
Y = [x(1+j:end) x(1:end-j)];
ACY = corrcoef(Y);
y(j,1) = ACY(2,1);
end