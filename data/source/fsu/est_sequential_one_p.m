%est_sequential_one_p.m 
clear all 
format compact

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

price= [a f m r]; 

priceix = [3 4 5 6];%column of the price, agr = 3, fuels = 4, metals = 5, interest rate =6 

for k=1:4 %one price at a time 

lagg(price(:,k));

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
np = 1; %# of world  prices 
 
pix = priceix(k); 

gx = eye(nv);

for i = 1:length(country_name); %country
for j = 1:4; %NIPA variables
%BP{i}(:,[3 4 5 6+j]); %data for domestic block
BP{i}(:,[pix 6+j]); %data for domestic block

lagg(ans);
D = ans(:,1:nv);
D1 = ans(:,nv+1:end);
T = size(D,1);
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
result(i,k, j) = sum(Vyr(1:np,end));%i is country, k is the world price shock; j is the NIPA variable


vF(1:nv,1:nv,i,k,j) = F;%first dimension in nv, second dimension is nv, third i is country i, fourth is the world price shock; and fifth is the  NIPA variable (index j)
vETA(1:nv,1:nv,i,k,j) = ETA;

end % for j = 1:4; %NIPA variables

TT(i,1) = T+1; %sample size domestic block. Why add 1? Because we took lags before and lost 1 obs. 

end %for i = 1:length(country_name); %country
end %for k=1:4 %one price at a time 

save est_sequential_one_p.mat vF vETA TT country_name result