function   table_tex = tabletex(x,n,digits,first_column);
%table_tex = tabletex(x,n,digits) transforms the numerical table x into a LaTeX-READABLE string table named table_tex. 
% Inputs: x a matrix of numbers
%n is either the number of digits of precision  or the number of decimals.
%The default value for n is 2. 
%digits takes  the value 1 or 0. The default value for digits is 1. 
%If digits == 1, then n refers to number of decimals. 
%If digits == 0, then n is the number of digits of precision.
%table_ex = tabletex(x,n,digits,first_column) 
%adds a column on the left side. the input first_column is a cell with strints. Example 
%first_column={'$y_t$';'$x_t'}
%(c) S. Schmitt-Grohe and M. Uribe, 1998. 

if nargin < 2
n=2;
end

if nargin<3
digits = 1; %this means that the default is the number of digits after the period. 
end

[r,c]=size(x);

amper = [];
w = [];

for j=1:r
amper = [amper;'&  '];
w = [w;'\\'];
end


table_tex=[];

if digits==1
    %choose this to control the number of digits after the decimal point. 
    s = ['%6.' num2str(n) 'f'];
    %
else
    s = n;
end


for j=1:c

    table_tex =  [table_tex amper num2str(x(:,j),s)];

end

table_tex = [table_tex w];
 
if nargin==4
w = table_tex;
table_tex = cell(0);
for j=1:size(x,1)
table_tex{j}=[first_column{j}   w(j,:)];
end
table_tex = char(table_tex);
end

%clc
%table_tex
