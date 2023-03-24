%est_joint_r.m 
%estimates the matrices F, G, and Sigma defining the SVAR system of owrld prices and output at quarterly frequency   for 86  countries and cross-coutnry averages for countries with 100 quarters of data or more. The SVAR system is the one studied in the paper ``World Shocks, World Prices,  And Business Cycles.''
%The main output  of this program is the 86x1 variable `result'  containg  estimates of the share of the variances of otuput  explained by world shocks mediated by commodity prices and the interest rate  for the 86  countries in thequarterly  panel. 
%© Martín Uribe, June 2016 

clear all 
format compact
clc

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
for i = 1:nc; %country
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

end

tt= find(TT>=100);
disp([ num2str(round(nanmedian(result(tt,:))*100)/100) '  median'])

disp([ num2str(round(nanmean(result(tt,1))*100)/100) '  mean'])

save est_joint_r.mat 