%bias_joint_r_run.m
%estimate the small-sample bias in the contribution of world shocks mediated by commodity prices and the world interest rate  to the variance of NIPA variables country-by-country. 
%© Martín Uribe, June 2016 

clear all
warning('off') %this program estimates possibly hundred of thousands of SVARs with artificial data. Every time  an SVAR fails to have all unit roots inside the unit circle, matlab generates a warning message. 

load est_joint_r.mat vF vETA TT country_name
%produced  with est_joint_r.m

nv = 8; %number of variables in the SVAR
np = 4; %number of world prices

 Tmontecarlo = 1e3;  %#of simulaitons

for i=1:numel(country_name);
hx = vF(:,:,i);
ETA = vETA(:,:,i);
b(i,1:nv-np) =  bias_est(hx,ETA,TT(i),55,np, Tmontecarlo);
disp(i);
disp(b(i,:));
end

save bias_joint_r_run.mat

%to produce the a table after running this program, simply run the folloiwng lines (no need to run the above lines again):

first_column ={'Noncorrected Estimate';'Small-Sample Bias';'Corrected Estimate';'MAD of Corrected Estimate'}

load bias_joint_r_run.mat
%produced by running the present  program (bias_joint_r_run.m)
load est_joint_r.mat result
%produced  with est_joint_r.m


table = [nanmedian(result) 
nanmedian(b)
nanmedian(result-b)
mad(result-b)
]

%tabletex(table,2,1,first_column)