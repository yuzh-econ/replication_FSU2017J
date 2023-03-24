clear all
load est_joint_r.mat 
load est_post04.mat vA vSigma_mu tt
A  = vA(:,:,tt(1));
Sigma_mu = vSigma_mu(:,:,tt(1));

for i = 1:nc; %country
if ismember(i,tt)

B = vB(:,:,i);
C = vC(:,:,i);
D = vD(:,:,i);
Sigma_e = vSigma_e(:,:,i);
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

save est_fin.mat 