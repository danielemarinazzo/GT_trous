function [cc lag]=func_corr(x,y,delta)
n=length(x);
lag=-delta:delta;
for i=0:delta
    cc(i+delta+1)=corr(x(1:n-i),y(1+i:n));
end
for i=1:delta
    cc(delta-i+1)=corr(y(1:n-i),x(1+i:n));
end