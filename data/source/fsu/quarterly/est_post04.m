%est_post04.m 
%re-estimates the matrices
% A and Sigma_mu defining the foreign bloc  and the matrices B, C, D, and Sigma_e defining the domestic bloc of the SVAR system 
%on the post 2003 sample. 
% the foreign bloc consists of: pa, pf, pm, and r
%quarterly data
%38 countries, to be included a country must have at least 100 consecutive observations. 
% The SVAR system is the one studied in the paper ``World Shocks, World Prices,  And Business Cycles.''
%The main output  of this program is the 86x1 variable `result'  containg  estimates of the share of the variances of otuput  
%explained by world shocks mediated by commodity prices and the interest rate  for the 86  countries in the quarterly  panel. 
%© Martín Uribe, June 2016 

clear all 
format compact
clc

load vdet_hp1600 BP country_name readme b raw_data 
%produced by running vdet.m in z:\joint\isom16\quarterly

nc = length(country_name);  %number of countries

% set the sample to start in 2004:Q1
for i=1:nc
tfirst = find(BP{i}(:,1)==2004)
BP{i} = BP{i}(tfirst:end,:)
end


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

%HP Filter (1600) detrend on the entire sample 
w = 1600; %smoothing parameter value
a = hpfilter(A,w);
f = hpfilter(F,w);
m = hpfilter(M,w);
r = hpfilter(R,w);

%set the new sample to start in 2004Q1
tlast1 = find(raw_data(:,1)==2003.75)
a = a(tlast1+1:end);
f = f(tlast1+1:end);
m = m(tlast1+1:end);
r = r(tlast1+1:end);


%re-estimate matrices A and Sigma_mu of the foreign bloc
lagg([a f m r]);
p = ans(:,1:4);
p1 = ans(:,5:end);
T = size(p,1);
cons = ones(T,1);
X = [p1 cons];
b = X\p;
A = b(1:end-1,:)';
mu = p-X*b;
Sigma_mu = cov(mu);

nv = 5; %# of variables in SVAR
np = 4; %# of world prices 
gx = eye(nv);

nc = length(country_name); %# of countries
vA = zeros(np,np,nc);
vSigma_mu = zeros(np,np,nc);
vB = zeros(nv-np,np,nc);
vC = zeros(nv-np,nv-np,nc);
vD = zeros(nv-np,np,nc);
vSigma_e = zeros(nv-np,nv-np,nc);
vF = zeros(nv,nv,nc);
vETA = zeros(nv,nv,nc);

load est_joint_r tt
%produced by running est_joint_r.m 
%tt contains the index of countries with at least 100 quarterly
%observations. (There are 86 countries, but only 38 have at least 100
%quarters of data.)


for i = 1:nc; %country
if ismember(i,tt)
BP{i}(:,[2 3 4 5 6]); %data for domestic block
lagg(ans);
D = ans(:,1:nv);
D1 = ans(:,nv+1:end);
T = size(D,1);
TT(i,1) = T+1; %sample size domestic block
cons = ones(T,1);
X = [D1 D(:,1:np) cons];
Y = D(:,np+1:nv);
b = X\Y;
e = Y-X*b;
b = b(1:end-1,:)';
B = b(:,1:np);
C = b(:,np+1:nv);
D = b(:,nv+1:nv+np);

Sigma_e = cov(e);

F = [A zeros(np,nv-np);
        D*A+B C];
G = [eye(np) zeros(np,nv-np); 
         D eye(nv-np)];
Sigma = [Sigma_mu zeros(np,nv-np);
                  zeros(nv-np,np) Sigma_e];
%variance decomposition 
ETA = chol(G*Sigma*G','lower');
 Vyr = variance_decomposition(gx,F,ETA);
result(i,1:nv-np) = sum(Vyr(1:np,np+1:end));
disp([ num2str(round(result(i,:)*100)/100) '  ' country_name{i}])

vA(:,:,i) = A;
vSigma_mu(:,:,i) = Sigma_mu;
vB(:,:,i) = B;
vC(:,:,i) = C;
vD(:,:,i) = D;
vSigma_e(:,:,i) = Sigma_e;
vF(:,:,i) = F;
vETA(:,:,i) = ETA;
end %if ismember(i,tt)
end

disp([ num2str(round(nanmedian(result(tt,:))*100)/100) '  median'])

disp([ num2str(round(nanmean(result(tt,1))*100)/100) '  mean'])

save est_post04.mat 