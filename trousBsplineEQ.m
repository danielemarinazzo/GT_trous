function [w c]=trousBsplineEQ(x,J)
%calcola la trasformata trous usando la trasformata di haar
%m ordine modello autoregressivo
N=length(x);
w=zeros(J,N);
c=zeros(J,N);
h=[1 4 6 4 1]/16;
%Aggiungo artificialmente i punti prima, in modo da riempire anche i primi
%coefficienti
for t=3:N-2
    c(1,t)=h*x(t-2:t+2);
    w(1,t)=x(t)-c(1,t);
end
cc=3;
    for j=2:J
        cc=cc+2^j;
        for t=cc:N-cc+1
        a=0;for f=-2:2;a=a+h(f+3)*c(j-1,t+2^(j-1)*f);end
        c(j,t)=a;
        w(j,t)=c(j-1,t)-c(j,t);
    end
end




w=w(:,cc:N-cc+1);
c=c(:,cc:N-cc+1);

