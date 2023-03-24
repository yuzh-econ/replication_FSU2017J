%bias_post04.m
%estimate the small-sample bias in the contribution of world shocks mediated by commodity prices and the world interest rate  to the variance of output  country-by-country. 
%© Martín Uribe, June 2016 

clear all
warning('off') %this program estimates possibly hundred of thousands of SVARs with artificial data. Every time  an SVAR fails to have all unit roots inside the unit circle, matlab generates a warning message. 

load est_post04.mat vF vETA TT country_name tt
%produced  with est_joint_r.m in 
%z:\joint\isom16\quarterly

nv = 5; %number of variables in the SVAR
np = 4; %number of world prices

 Tmontecarlo = 1e3;  %#of simulaitons

for i=1:numel(country_name);
if ismember(i,tt);
hx = vF(:,:,i);
ETA = vETA(:,:,i);
b(i,1:nv-np) =  bias_est(hx,ETA,TT(i),TT(i),np, Tmontecarlo);
disp(i);
disp(b(i,:));
end %if ismember(i,tt)
end

country_name_q100  = country_name(tt);

save bias_post04.mat

%to produce the a table after running this program, simply run the folloiwng lines (no need to run the above lines again):

first_column ={'Noncorrected Estimate';'Small-Sample Bias';'Corrected Estimate';'MAD of Corrected Estimate'}

load bias_post04.mat
%produced by running the present  program (bias_post04.m)
load est_post04.mat result
%produced  with est_post04.m

table = [nanmedian(result(tt)) 
nanmedian(b(tt))
nanmedian(result(tt)-b(tt))
mad(result(tt)-b(tt))
]

tabletex(table,2,1,first_column)