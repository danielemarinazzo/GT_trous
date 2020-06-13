function [V d ifail]=chop(VV,dd,val);
if nargin <2
    val=1.e-8;
end
ifail=0;
f=1.e-8;
[n nv]=size(VV);
v=[VV' dd]; 
ind=find(abs(imag(v))<val);
v(ind)=real(v(ind));
d=v(:,n+1);
ind=find(d>f*max(d) & d >f);
V=v(ind,(1:n))';
if sum(sum(abs(imag(v(ind,:)))))>0;
    ifail=2;
end
