%est_sequential_pc.m 
%est_sequential.m 
%estimates the matrices F, G, and Sigma defining the SVAR system of owrld prices and NIPA variables (Y,C,I,TBY) one at the time for 138 countries. The SVAR system is the one studied in the paper ``World Shocks, World Prices,  And Business Cycles.'' This program considers the case in which the foreign bloc is the first principal compoenent of [pa, pf, pm, r],.
%The main output  of this program is the 138x4 variable `result'  containg  estimates of the share of the variances of the NIPA variables explained by world shocks mediated by the first principal component of commodity prices and the interest rate  for the 138 countries in the panel. 
%© Martín Uribe, June 2016 
clear all 
warning('off')
format compact
clc

load vdet_hp100 BP country_name readme b raw_data 
%produced by running vdet.m in z:\joint\isom16\138cou

%Agriculture
i=3; 
A = log(raw_data(:,i,1));

%Fuels
i=4; 
F = log(raw_data(:,i,1));

%Metals
i=5; 
M = log(raw_data(:,i,1));

%Interest Rate
i=6; 
R = log(1+raw_data(:,i,1)/100);

%HP Filter (100)
a = hpfilter(A,100);
f = hpfilter(F,100);
m = hpfilter(M,100);
r = hpfilter(R,100);

p = [a f m r];

%Principal Component
[loadings,pc] = princomp(p);

pc1 = pc(:,1); %main pc
tpc = raw_data(:,1,1);  %date of principal component

lagg(pc1);
p = ans(:,1);
p1 = ans(:,2);
T = size(p,1);
cons = ones(T,1);
X = [p1 cons];

b = X\p;
A = b(1:end-1,:)';
mu = p-X*b;
Sigma_mu = cov(mu);
R2p = 1-var(mu)./var(p);

nv = 2; %# of variables in SVAR
npc = 1; %# of principal components included
gx = eye(nv);
for i = 1:length(country_name); %country
for j = 1:4; %j=[1 2 3 4 ]==>[Y C I TBY]
d = BP{i}(:,[6+j]); %data for domestic block
t1 = BP{i}(1,1);
t2 = BP{i}(end,1);
itpc1 = find(tpc==t1);
itpc2 = find(tpc==t2);
 [pc1(itpc1:itpc2) d]; 
lagg(ans);
D = ans(:,1:nv);
D1 = ans(:,nv+1:2*nv);
T = size(D,1);
cons = ones(T,1);
X = [D1 D(:,1:npc) cons];
Y = D(:,end);
b = X\Y;
e = Y-X*b;
b = b(1:end-1,:)';
B = b(1:npc);
C = b(npc+1:nv);
D = b(nv+1:nv+npc);

Sigma_eps = cov(e);

F = [A zeros(npc,nv-npc);
        D*A+B C];
G = [eye(npc) zeros(npc,nv-npc); 
         D eye(nv-npc)];
Sigma = [Sigma_mu zeros(npc,nv-npc);
                  zeros(nv-npc,npc) Sigma_eps];
%variance decomposition 
ETA = chol(G*Sigma*G','lower');
 [Vyr,Vxr,Vy,Vx]=variance_decomposition(gx,F,ETA);
result(i,j) = sum(Vyr(1:npc,end));

vF(1:nv,1:nv,i,j) = F;
vETA(1:nv,1:nv,i,j) = ETA;
end

TT(i,1) = T+1; %sample size domestic block. Why add 1? Because we took lags before and lost 1 obs. 

disp([ num2str(round(result(i,:)*100)/100) '  ' country_name{i}])

end

disp([ num2str(round(median(result)*100)/100) '  median'])

disp([ num2str(round(mean(result)*100)/100) '  mean'])

save  est_sequential_pc.mat