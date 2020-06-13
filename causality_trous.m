function [cb cbt ifail rr pp]=causality_trous(X,x,type,par,th)
% Input:    X       : matrix (nvar x m x n) of driver data;
%           x       : matrix n x nt  target variables
%                       n    = number of samples
%           type    : kernel function 'p' polinomial 'g' gaussian ;
%           par     : parameter of the kernel function;
%           m       : order of the model   % default m=1;
% Output:   cb      : cb(i,j) = i->j
%           ifail   :  0 no error 
%                      1 cholesky algoritm fail
%                      2 complex eigensystem
%                      3 error in polypower
%Last versione (11/01/2011)
[nvar m n]=size(X);
[n nt]=size(x);
f=1.e-6;
%th=0.05;
Xr=reshape(X,nvar*m,n);
rr=0;
pp=0;
[VV D ifail V1 D1 xnorm]=filtro(Xr,type,par,f,true);
if ifail>0
    return
end
VT=VV*D.^0.5;
polycall=true;
kk=0;
for driver=1:nvar
    XX=X;
    XX(driver,:,:)=[];
    Xr=reshape(XX,(nvar-1)*m,n);
    [V D ifail L2 VN2 D2]=filtro(Xr,type,par,f,polycall);
    polycall=false;
    [VN ifail]=vnorma(VT,V,VV);
    if ifail>0
        return
    end
    xv=x-V*V'*x;
    %for k=1:length(target)
        [rr ppt]=corr(xv,VN);
    %[rr ppt]=corr(xv(:,target(k)),VN);
   
rn=rr.^2; %questa ti fa tornare alla versione normalizzata della gc
%questa invece è l'unnormalized
% [rrt pp]=corr(x,VN);
% rn=rrt.^2;


cbt=sum(rn);
thb=th/length(rr);
indpr=find(ppt>thb);
rn(indpr)=0;
cb(driver,1:nt)=sum(rn,2);
    end

