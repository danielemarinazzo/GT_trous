function [VN ifail]=vnorma(VT,V,VV)
n=size(V,1)-1;
fast=size(VV,2)<n-2;
ifail=0;
if fast
    A=VV'*VT;
    B=(VV'*V)*(V'*VT);
    KKN=A'*A-B*A-A'*B'+B*B';
    [VVN D]=eig(KKN);
    VVN=VV*VVN;
    [VN d ifail]=chop(VVN,diag(D),1.e-8);
else
    K=VT*VT';
    P=V*V';
    KT=(K-P*K-K*P+P*K*P);
    [VVN D]=eig(KT);
    [VN d ifail]=chop(VVN,diag(D),1.e-8);
end
