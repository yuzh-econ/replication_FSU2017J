function [sB,stdB,Shx,SETA] = bias_est(hx,ETA,Tshort, Tlong,np, Tmontecarlo);
%[sB,stdB,Shx,SETA] = bias_est(hx,ETA,Tshort, Tlong,np) returns a measure of the small sample bias. The true system is assumed to be
%x_t = hx x_t-1 + ETA eps_t
%Elements 1:np of x_t, denoted p_t, are assumed to follow an autonomous AR(1) process and to be estimated with time series of length Tlong.
%The remaining elements of x_t, denoted y_t,  follows a process that depends upon past values of p_t and y_t. This  block is estimated with Tshort observations. 
%Input Tmontecarlo indicates the number of Monte Carlo simulations to be
%exectued. The default is Tmontecarlo = 1e3; 
%sB is a mesure of the bias in the estimated  contribution of p_t to the variance of y_t. The true contributiion is denoted VS. 
%The matrices Shx and SETA are the average finite-sample estimates of hx and ETA using artificial data stemming from the true model. 
%The bias is estimated by a Monte Carlo method as follows:
%(1) Produce an artificial data set of a desired length (say 250) from the true model. 
%(2) Use the last Tlong artificial observations to estimate the p_t block of the SVAR.
%(3) Use the last Tshort observations to estimate the y_t block. 
%(4) steps (2) and (3) yield an estimate of hx and ETA. Use these estimates to compute the share of the variance of y_t explained by p_t. 
%(5) Repeat (2)-(4) a desired number of times (say 1e3 times) and compute averages of the resulting estimates of hx ETA and VS, and denote them  Shx, SETA, and SVS. Then sB = SVS-VS. stdB is the standard deviation of the estimated biases.  
%� Mart�n Uribe, June 2016 

if nargin<6
    Tmontecarlo = 1e3; 
end


nv = size(hx,1); %# variables in SVAR
gx = eye(nv);

trueVD = variance_decomposition(gx,hx,ETA);
trueVS = sum(trueVD(1:np,np+1:end),1);

vsVS = [];
Shx = hx*0;
SETA = ETA*0;
for i=1:Tmontecarlo
Y = simu_1st(gx, hx, ETA, 250);%Simulate artifical times series. The order of Y is: rows are observations (250 years), the columns are the variables of the vector x_t. 
Dshort = Y(end-Tshort+1:end,:);
Dlong = Y(end-Tlong+1:end,:);
%Use the long sample, Tlong (which are 55 obs) to estimate the commodity price block 
%and the short sample (Tshort, which is different for each country) to estimate the domestic block. 

lagg(Dlong(:,1:np));
p = ans(:,1:np);
p1 = ans(:,np+1:end);
T = size(p,1);
cons = ones(T,1);
X = [p1 cons];
b = X\p;
A = b(1:end-1,:)';
mu = p-X*b;
Sigma_mu = cov(mu);

lagg(Dshort);
D = ans(:,1:nv);
D1 = ans(:,nv+1:end);
T = size(D,1);
cons = ones(T,1);
X = [D1 D(:,1:np) cons];
Y = D(:,np+1:end);
b = X\Y;
e = Y-X*b;
b = b(1:end-1,:)';
B = b(:,1:np);
C = b(:,np+1:nv);
D = b(:,nv+1:nv+np);
Sigma_eps = cov(e);

shx = [A zeros(np,nv-np);
        D*A+B C];
G = [eye(np) zeros(np,nv-np); 
         D eye(nv-np)];
Sigma = [Sigma_mu zeros(np,nv-np);
                  zeros(nv-np,np) Sigma_eps];
%variance decomposition 
sETA = chol(G*Sigma*G','lower');

%To filter a possible small sample bias in the estimation of the price VAR, unpercentage the following two lines:
%shx(1:np,1:np) = hxp;
%sETA(1:np,1:np) = ETAp;

V = variance_decomposition(gx,shx,sETA);

vsVS = [vsVS;sum(V(1:np,np+1:end),1)];

Shx = (i-1)/i*Shx + shx/i;
SETA = (i-1)/i*SETA + sETA/i;
end

sVS = nanmean(vsVS);
sB = sVS-trueVS;
stdB = nanstd(sVS);

if size(sVS)~=size(trueVS)
error('sVS not computed correctly')
end