function [V D ifail V1 D1 xnorm]=filtro(X,type,p,f,polycall)
if nargin<5
    polycall=true;
end
if nargin<4
    f=1.e-6;
end
[m n]=size(X);
if type=='g'
    fast=false;
    rmax=n;
else
    rmax=min(n,nchoosek(p+m,m));
    fast=nchoosek(p+m,m)<rmax+1;
end
if fast
    %polynomial kernel (p+m)!/(p!m!)-1< rmax'
    [L ifail]=Leval(X,m,p,polycall);
    if ifail>0
        V=0;
        D=0;
        return
    end
    [VN D1]=eig(L'*L);
    V1=L*VN;
    ifail=0;
else
    [L P ifa]=cholesky(X,type,p,rmax,f);
%     size(L)
    if ifa
        ifail=1;
        V=0;
        D=0;
        return
    end
    ifail=0;
    [VN D1]=eig(L'*L);
    V1=L(P,:)*VN;
end
xnorm=repmat(sqrt(dot(V1,V1)),n,1);
V1=V1./xnorm;
[s ind]=sort(abs(diag(D1)),'descend');
r=sum(s>f*s(1));
ind=ind(1:r);
V=V1(:,ind);
D=D1(ind,ind);
