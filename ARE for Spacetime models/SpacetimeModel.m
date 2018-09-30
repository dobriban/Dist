%Experiments of ARE for spacetime models, tridiagonal covariance, sample
%randomly distributed i.e. different gamma_i and G_i.
cd('C:\Dropbox\Projects\Distributed\Experiments\\ARE for Spacetime models')
%%
n = 10000;
p_arr = [10,20,50,100];
L = length(p_arr);
for l=1:length(p_arr)
    
p = p_arr(l);
gamma = p/n;

%create a test covariance matrix
rng(2);
savefigs=1;    closefigs=1;
Sigma1 = eye(p);
%Sigma1 = diag(1000*rand(p,1));
%Sigma2 = diag(v);
%Y = normrnd(0,1,[n,n]);
%Sigma2 = (Y'*Y)/n;
%Sigma2 = diag(eig(Sigma2));
Sigma2 = full(gallery('tridiag',n,-2,5,-2));

%generate a random matrix with prescribed covariance matrix
X = sqrtm(Sigma2)*normrnd(0,1,[n,p])*sqrtm(Sigma1);

%intialize ARE for each number of machine
ARE1 = zeros(1,1/gamma);

for i = 1:1/gamma
 a = p*ones(i,1);
    for k = 1:(n-i*p)
    j = randi(i);
    a(j) = a(j)+1;
    end
 a = tril(ones(i))*a;
t = 1/trace(inv(X(1:a(1),:)'*X(1:a(1),:)));
for s = 1:i-1
    t = t+1/trace(inv(X(a(s)+1:a(s+1),:)'*X(a(s)+1:a(s+1),:)));
end
ARE1(i) = trace(inv(X'*X))*t;
end

m = 1:1:1/gamma;
plot(m,ARE1,'LineWidth',3)
hold on 
ARE2 = (1/gamma-m)/(1/gamma-1);
plot(m,ARE2,':','LineWidth',3)
hold off
xlabel('Number of Machines');
ylabel('ARE');
set(gca,'fontsize',20);
xlim([min(m),max(m)]);
legend('Numerical','Optimal');
str=sprintf('p/n=%.3f' ,p/n);
title(str);

  if savefigs==1
        filename = ...
            sprintf( './Spacetime-randomly-distributed-tridiagonal-covariance-n=%d-p=%d-gamma=%.3f.png',...
            n,p,gamma);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        if closefigs==1
            close(gcf)
        end
  end
end

