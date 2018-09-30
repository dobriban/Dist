%Experiments of ARE for MP models, randomly distributed case, actual
%regression, not just compare the trace
n = 10000;
p = 100;
gamma = p/n;

for j = 1:10

rng(j);
J=j;
savefigs=1;    closefigs=1;
X = normrnd(0,1,[n,p])*eye(p);
ep = normrnd(0,1,[n,1]);

%intialize ARE for each number of machine
ARE1 = zeros(1,1/gamma);

for i = 1:1/gamma
 a = p*ones(i,1);
    for k = 1:(n-i*p)
    j = randi(i);
    a(j) = a(j)+1;
    end
 a = tril(ones(i))*a;
t = 1/(ep(1:a(1))'*X(1:a(1),:)*inv(X(1:a(1),:)'*X(1:a(1),:))*inv(X(1:a(1),:)'*X(1:a(1),:))*X(1:a(1),:)'*ep(1:a(1)));
for s = 1:i-1
    t = t+1/(ep(a(s)+1:a(s+1))'*X(a(s)+1:a(s+1),:)*inv(X(a(s)+1:a(s+1),:)'*X(a(s)+1:a(s+1),:))*inv(X(a(s)+1:a(s+1),:)'*X(a(s)+1:a(s+1),:))*X(a(s)+1:a(s+1),:)'*ep(a(s)+1:a(s+1)));
end
ARE1(i) = (ep'*X*inv(X'*X)*inv(X'*X)*X'*ep)*t;
end
     
m = 1:1:1/gamma;
plot(m,ARE1,'LineWidth',3)
hold on 
ARE2 = (1-m*gamma)./(1-gamma);
plot(m,ARE2,':','LineWidth',3)
hold off
xlabel('Number of Machines');
ylabel('ARE');
set(gca,'fontsize',20)
legend('Numerical','Theoretical');
str = sprintf( 'p/n=%.3f',p/n);
title(str);

  if savefigs==1
        filename = ...
            sprintf( './Actual-n=%d-p=%d-gamma=%.3f-rng(%d).png',...
            n,p,gamma,J);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        if closefigs==1
            close(gcf)
        end
  end
end