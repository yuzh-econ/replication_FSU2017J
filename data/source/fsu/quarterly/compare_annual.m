%compare_annual.m 
%Produces the bias-corrected contribution of world shocks to the variance of output using annual data for the set of (38) countries for which quarterly output data is available for 100 quarters or more
%© Stephanie Schmitt-Grohé and Martín Uribe, August 2016 

clear all

load ..\138cou\bias_sequential_r_run
%(exact location of file  must be adjusted by user) 
%produced by running
%bias_sequential_r_run.m  in 
% z:\joint\isom16\138cou

load ..\138cou\est_sequential_r result
%(exact location of file  must be adjusted by user)
%produced by running
%est_sequential_r.m  in 
% z:\joint\isom16\138cou
%Annual data 138 countries
nc138 =  nanmedian(result(:,1));
b138 = nanmedian(b(:,1));
result-b;
r138 = nanmedian(ans(:,1)); 
mad138 = mad(ans(:,1));
table138 = [nc138;b138;r138;mad138]

load bias_joint_r_run country_name_q100
%(exact location of file  must be adjusted by user)
%produced by running
%bias_joint_r_run.m  in 
%z:\joint\isom16\quarterly


n = length(country_name_q100);

for i=1:n
[~,I(i,1)] = intersect(country_name,country_name_q100(i));
end

%Annual data 38 countries
nca38 =  nanmedian(result(I,1));
ba38 = nanmedian(b(I,1));
result-b;
ra38 = nanmedian(ans(I,1)); 
mada38 = mad(ans(I,1));
tablea38 = [nca38;ba38;ra38;mada38]

clear result b tt
load bias_joint_r_run b
load est_joint_r result tt


%Quarterly data 38 countries
I = tt;
ncq38 =  nanmedian(result(I,1));
bq38 = nanmedian(b(I,1));
result-b;
rq38 = nanmedian(ans(I,1)); 
madq38 = mad(ans(I,1));
tableq38 = [ncq38;bq38;rq38;madq38]

table = [tableq38 tablea38 table138]



first_column ={'Noncorrected Estimate';'Small-Sample Bias';'Corrected Estimate';'MAD of Corrected Estimate'}

%tabletex(table,2,1,first_column)