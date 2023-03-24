%compare_corr.m computes correlations and standard deviations of HP(1600) real commodity prices over three samples: 1960:Q1-2015:4, 1960:Q1-2003:Q4, and 2004:Q1-2015:Q4.
%© Stephanie Schmitt-Grohé and Martín Uribe, August 2016 

clear all 
format compact

load vdet_hp1600 BP country_name readme b raw_data 
%produced by running vdet.m in z:\joint\isom16\quarterly

nc = length(country_name);  %number of countries

%Agriculture
i=2; 
A = log(raw_data(:,i,1));

%Fuels
i=3; 
F = log(raw_data(:,i,1));

%Metals
i=4; 
M = log(raw_data(:,i,1));

%R
i=5; 
R = log(1+raw_data(:,i,1)/100);

%HP Filter (1600)
w = 1600; %smoothing parameter value
a = hpfilter(A,w);
f = hpfilter(F,w);
m = hpfilter(M,w);
r = hpfilter(R,w);

t2003Q4 = find(raw_data(:,1,1)==2003.75);

corr([a f m ]);
all_sample = [ans(2,1) ans(3,1) ans(2,3) std([a f m])]';

a1 = a(1:t2003Q4); 
f1 = f(1:t2003Q4); 
m1 = m(1:t2003Q4); 
corr([a1 f1 m1]);
pre04 = [ans(2,1) ans(3,1) ans(2,3) std([a1 f1 m1])]';

a2 = a(1+t2003Q4:end); 
f2 = f(1+t2003Q4:end); 
m2 = m(1+t2003Q4:end); 
corr([a2 f2 m2]);
post04 = [ans(2,1) ans(3,1) ans(2,3) std([a2 f2 m2])]';

table = [all_sample pre04 post04]

first_column = {'$\rho(p^a,p^f)$', '$\rho(p^a,p^m)$', '$\rho(p^f,p^m)$', '$\sigma(p^a)$',  '$\sigma(p^f)$',  '$\sigma(p^m)$'}

tabletex(table,2,1,first_column)