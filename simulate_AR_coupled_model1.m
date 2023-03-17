function X=simulate_AR_coupled_model1(N,a)
% a is the coupling
settleTime=1000;
N = N + settleTime;
X = 0.5*randn(2,N);

for i=9:N
    X(1,i) = X(1,i) + 0.3.*X(1,i-1);
    X(2,i) = X(2,i) + 0.1.*X(2,i-1) + a.*X(1,i-8);
end
X = X(:,settleTime+1:end);