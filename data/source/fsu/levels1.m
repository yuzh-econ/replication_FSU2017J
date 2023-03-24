%levels1.m
%plot the level and cyclical components of real commodity prices and second moments of the cyclical components of commodity prices and the real world interest rae.  
%© Martín Uribe, June 2016 

clear all , clf
orient tall
format compact

load vdet_hp100  BP country_name readme b raw_data 
%the string variable readme explains the content of vdet_hp100.mat
%produced by running vdet.m in z:\joint\isom16\138cou

%Agriculture
i=3; 
A = (raw_data(:,i,1));

%Fuels
i=4; 
F = (raw_data(:,i,1));

%Metals
i=5; 
M = (raw_data(:,i,1));



%Metals
i=6; 
R = (raw_data(:,i,1));

%date
t = raw_data(:,1,1);

subplot(3,2,1)
h=plot(t,A/A(1),'-')
title('Price Level, Agricultural Commodities')
grid 
xlim([t(1) t(end)+1])
set(gca, 'XTick', [1960 1970 1980 1990 2000 2010 ])
ylim([0.2 1.6])

subplot(3,2,3), 
plot(t,M/M(1),'-');
title('Price Level, Metals')
grid
xlim([t(1) t(end)+1])
set(gca, 'XTick', [1960 1970 1980 1990 2000 2010 ])
ylim([0.2 1.6])

subplot(3,2,5)
plot(t,F/F(1),'-')
title('Price Level, Fuels')
xlim([t(1) t(end)+1])
grid
xlim([t(1) t(end)+1])
set(gca, 'XTick', [1960 1970 1980 1990 2000 2010 ])

[a,at] = hpfilter(log(A),100);
[f,ft] = hpfilter(log(F),100);
[m,mt] = hpfilter(log(M),100);
[r,rt] = hpfilter(log((1+R/100)),100);

subplot(3,2,2)
plot(t,a,'-')
title('Cyclical Component, Agricultural Commodities')
grid
xlim([t(1) t(end)+1])
set(gca, 'XTick', [1960 1970 1980 1990 2000 2010 ])
ylim(0.4*[-1 1])


subplot(3,2,4)
plot(t,m,'-')
title('Cyclical Component, Metals')
grid
xlim([t(1) t(end)+1])
set(gca, 'XTick', [1960 1970 1980 1990 2000 2010 ])
ylim(0.4*[-1 1])

subplot(3,2,6)
plot(t,f,'-')
title('Cyclical Component, Fuels')
grid
xlim([t(1) t(end)+1])
set(gca, 'XTick', [1960 1970 1980 1990 2000 2010 ])
axis('tight')

shg

x = [a  m f r];
table = [std(x);
[acf(x(:,1)) acf(x(:,2)) acf(x(:,3)) acf(x(:,4))];
corr(x)];

for i=1:size(BP)
stdy(i,1) = std(BP{i}(:,7));
end

table = [table;std(x)/mean(stdy)]

col1 = {'Standard Deviation, $\sigma(p)$'; 'Serial Correlation, $\rho(p)$';  'Correlation with Agri., $\rho(p^a,p)$'; 'Correlation with Metals, $\rho(p^m,p)$'; 'Correlation with Fuels, $\rho(p^f,p)$';  'Correlation with Interest Rate, $\rho(r,p)$'; 'Relative Std.Dev, $\sigma(p)/\sigma(GDP)$'};

%tabletex(table,2,1,col1)