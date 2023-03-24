%bias_sequential_one_run.m
%estimate the small-sample bias in the contribution of world shocks mediated by commodity prices to the variance of NIPA variables  (Y, C, I, and TBY)  country-by-country and one NIPA variable at the time (sequential estimation). 
%Output: matrix b of order 138x4 where b(i,j) represents the bias in country i and variable j. 
%© Stephanie Schmitt-Grohé and Martín Uribe, June 2016 

clear all
warning('off') 
%This program estimates possibly hundred of thousands of SVARs with artificial data. 
%Every time an SVAR fails to have all roots inside the unit circle, Matlab
%will generates a warning message. To avoid this, warning is set to off. 

load est_sequential_one_p.mat vF vETA TT country_name
%produced by running est_sequential_one_p.m
%vF is a 4x4x138x4 matrix containing the matrices F of the SVAR in the paper. 
%The first two dimensions of vF is the matrix F for one SVAR. 
%The third dimension is the country i, for i=1:38;  
%And the fourth dimension is the NIPA variable included in the SVAR (since
%the estimation is sequential). For example vF(:, :, 1,1) is the matrix F
%for country i=1 and NIPA variable 1, which is output. 
%vETA is a 4x4x138x4 matrix containing the lower Cholesky decomposition of the matrix G*SIGMA*G' of the  SVAR given in the paper . 
%TT is a 138x1 vector containing the sample size for the domestic block for
%each country. 

nv = 2; %total number of variables in the SVAR, 1 world price plus 1 NIPA. 
np = 1; %total number of prices in the SVAR, 1. 
Tmontecarlo = 1e3; %number of draws for the Monte Carlo experiment. 
b=zeros(numel(country_name), 4,4); 
for k = 1:4%one regression for each of 4 possible world price shocks 
for i=1:numel(country_name);
for j=1:4 %NIPA variables
hx = vF(:,:,i,j); %matrix F for country i, for SVAR including NIPA variable j
ETA = vETA(:,:,i,k, j);
b(i,k,j) =  bias_est(hx,ETA,TT(i),55,np, Tmontecarlo);
end %for j=1:4 %NIPA variables
disp(i)
disp(squeeze(b(i,k,:))')
end %for i=1:numel(country_name);
end %for k = 1:4%one regression for each of 4 possible world price shocks 

save bias_sequential_one_p_run.mat

load bias_sequential_one_p_run.mat
%produced by running the present  program (bias_sequential_one_p_run.m)
load est_sequential_one_p.mat result
%produced  with est_sequential_one_p.m

%Median variance shares for 1 world price at a time
table = squeeze([nanmedian(result-b,1)])

%Best Single Price for each NIPA variable. 
share=result-b;
bestshare = squeeze(max(share, [],2));%this takes the max across the 4 different world prices considered. 
table = [table;  nanmedian(bestshare,1)];

%Best Single World Price for Output only
%for the list of countries by shock, see the script one_p_best_y_table.m

[max_share, max_shareix]=max(squeeze(share(:,:, 1)),[],2);
bestshareY=zeros(length(country_name), 4);  
for i=1:length(country_name)
bestshareY(i, 1:4) = share(i, max_shareix(i), :); 
end

table = [table;  nanmedian(bestshareY,1)];


first_column = {'One World Price, $p^a$';'One World Price, $p^f$';'One World Price, $p^m$';'One World Price, $r$';'Best Single World Price'; 'Best Single World Price for $y$'}

tabletex(table,2,1,first_column)

