function [w c W C]=trousBspline(x,J,m)
%calcola la trasformata trous usando la trasformata di haar
%m ordine modello autoregressivo
N=length(x);
y=zeros(1,2*N);y(N+1:2*N)=x;
w=zeros(J,2*N);
c=zeros(J,2*N);
h=[1 4 6 4 1]/16;
%Aggiungo artificialmente i punti prima, in modo da riempire anche i primi
%coefficienti
for t=N+1:2*N
    c(1,t)=h*y(t-4:t)';
    w(1,t)=y(t)-c(1,t);
    for j=2:J
        a=0;for f=1:5;a=a+h(f)*c(j-1,t-2^(j-1)*(f-1));end
        c(j,t)=a;
        w(j,t)=c(j-1,t)-c(j,t);
    end
end

%ora estraggo il data set
W=zeros(J,m,2*N);
C=zeros(J,m,2*N);
%C=zeros(m,2*N);
for t=N+1:2*N
   for j=1:J
       for k=1:m
           W(j,k,t)=w(j,t-(k-1)*2^j);
           C(j,k,t)=c(j,t-(k-1)*2^j);
       end
   end
%    for k=1:m
%        W(J+1,k,t)=c(J,t-(k-1)*2^J);
%    end
end



w=w(:,N+1:2*N);
c=c(:,N+1:2*N);
W=W(:,:,N+1:2*N);
C=C(:,:,N+1:2*N);
