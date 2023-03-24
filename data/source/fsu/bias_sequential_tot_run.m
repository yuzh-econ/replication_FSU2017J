%bias_sequential_tot_run.m
%estimate the small-sample bias in the contribution of world shocks mediated by the terms of trade  to the variance of NIPA variables  (Y, C, I, and TBY)  country-by-country and one NIPA variable at the time (sequential estimation). 
%Output: matrix b of order 138x4 where b(i,j) represents the bias in country i and variable j. 
%© Martín Uribe, June 2016 

clear all
warning('off') 
%This program estimates possibly hundred of thousands of SVARs with artificial data. 
%Every time an SVAR fails to have all roots inside the unit circle, Matlab
%will generates a warning message. To avoid this, warning is set to off. 

load est_sequential_tot.mat vF vETA TT country_name
%produced by running est_sequential_tot.m
%vF is a 2x2x138x4 matrix containing the matrices F of the SVAR in the paper. 
%The first two dimensions of vF is the matrix F for one SVAR. 
%The third dimension is the country i, for i=1:38;  
%And the fourth dimension is the NIPA variable included in the SVAR (since
%the estimation is sequential). For example vF(:, :, 1,1) is the matrix F
%for country i=1 and NIPA variable 1, which is output. 
%vETA is a 2x2x138x4 matrix containing the lower Cholesky decomposition of the matrix G*SIGMA*G' of the  SVAR given in the paper . 
%TT is a 138x1 vector containing the sample size for the domestic block for
%each country. 

nv = size(vF,1); %total number of variables in the SVAR, here 1 TOTplus 1 NIPA. 
np = nv-1; %total number of prices in the SVAR, 1 (TOT). 
Tmontecarlo = 1e3; %number of draws for the Monte Carlo experiment. 

for i=1:numel(country_name);
for j=1:4
hx = vF(:,:,i,j); %matrix F for country i, for SVAR including NIPA variable j
ETA = vETA(:,:,i,j);
b(i,j) =  bias_est(hx,ETA,TT(i),TT(i),np, Tmontecarlo);
end
disp(i);
disp(b(i,:))
end

save bias_sequential_tot_run.mat

%to produce a table after running this program, simply run the following lines (no need to run the above lines uagain):

first_column ={'Noncorrected Estimate';'Small-Sample Bias';'Corrected Estimate';'MAD of Corrected Estimate'}

load bias_sequential_tot_run.mat
%produced by running the present  program (bias_sequential_tot_run.m)
load est_sequential_tot.mat result
%produced  with est_sequential_tot.m

table = [nanmedian(result) 
nanmedian(b)
nanmedian(result-b)
mad(result-b)
]

tabletex(table,2,1,first_column)