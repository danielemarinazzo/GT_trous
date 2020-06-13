function [cb co]=causality_scale(x,J,type,par,m)
% x data matrix samples x variables
[n m1]=size(x);
th=0.05/(m1*(m1-1)/2);
for i=1:m1
    [w c W C]=trousBspline(squeeze(x(:,i)),J,1);
    for j=i+1:m1
        [w1 c1 W1 C1]=trousBspline(x(:,j),J,1);
        for h=1:J
            cb(i,j,h)=causality_trials(w(h,1:n-1)',x(1:n-1,j),x(2:n,j),type,par,m,1,th);
            cb(j,i,h)=causality_trials(w1(h,1:n-1)',x(1:n-1,i),x(2:n,i),type,par,m,1,th);
        end
        
        for h=1:J
            [co(i,j,h) p]=corr(w(h,:)',w1(h,:)');
            if p>th
                co(i,j,h)=0;
            end
            co(j,i,h)=co(i,j,h);
        end
    end
end

end

