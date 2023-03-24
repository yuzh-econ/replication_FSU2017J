%est_sequential_tot.m 
%estimates the matrices F, G, and Sigma defining the SVAR system of TOT and NIPA variables (Y,C,I,TBY) one at the time for 138 countries. The SVAR system is the one studied in the paper ``World Shocks, World Prices,  And Business Cycles.''
%The main output  of this program is the 138x4 variable `result'  containing  estimates of the share of the variances of the NIPA variables explained by world shocks mediated by the terms of trade for the 138 countries in the panel. 
%© Martín Uribe, June 2016 

clear all 
warning('off')
format compact
clc

load vdet_tot_hp100 BP country_name readme b raw_data 
%produced by running vdet.m in z:\joint\isom16\138cou



nv = 2; %# of variables in SVAR
np = 1; %size of foreign bloc
gx = eye(nv);

for i = 1:length(country_name); %country

%Foreign bloc (Terms of Trade, tot)
tot = BP{i}(:,2); %Hp filtered log of TOT
lagg([tot]);
p = ans(:,1:1);
p1 = ans(:,2:end);
T = size(p,1);
cons = ones(T,1);
X = [p1 cons];
b = X\p;
A = b(1:end-1,:)';
mu = p-X*b;
Sigma_mu = cov(mu);
R2p = 1-var(mu)./var(p);

for j = 1:4; %NIPA variables
BP{i}(:,[2 6+j]); %data for domestic block
lagg(ans);
D = ans(:,1:nv);
D1 = ans(:,nv+1:end);
T = size(D,1)
cons = ones(T,1);
X = [D1 D(:,1:np) cons];
Y = D(:,end);
b = X\Y;
e = Y-X*b;
b = b(1:end-1,:)';
B = b(1:np);
C = b(np+1:nv);
D = b(nv+1:nv+np);
Sigma_eps = cov(e);
F = [A zeros(np,nv-np);
        D*A+B C];
G = [eye(np) zeros(np,nv-np); 
         D eye(nv-np)];
Sigma = [Sigma_mu zeros(np,nv-np);
                  zeros(nv-np,np) Sigma_eps];
%variance decomposition 
ETA = chol(G*Sigma*G','lower');
 [Vyr,Vxr,Vy,Vx]=variance_decomposition(gx,F,ETA);
result(i,j) = sum(Vyr(1:np,end));

vF(1:nv,1:nv,i,j) = F;
vETA(1:nv,1:nv,i,j) = ETA;

end % for j = 1:4; %NIPA variables

TT(i,1) = T+1; %sample size domestic block. Why add 1? Because we took lags before and lost 1 obs. 

disp([ num2str(round(result(i,:)*100)/100) '  ' country_name{i}])

end % for i = 1:length(country_name)

disp([ num2str(round(nanmedian(result)*100)/100) '  median'])

disp([ num2str(round(mad(result,1)*100)/100) '  mad'])

disp([ num2str(round(nanmean(result)*100)/100) '  mean'])


save  est_sequential_tot.mat vF vETA TT country_name result R2p Sigma_mu A