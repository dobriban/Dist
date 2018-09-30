%Experiments of ARE for MP models, randomly distributed case.
cd('C:\Dropbox\Projects\Distributed\Experiments\ARE for Marchenko-Pastur models')
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
v = rand(p,1)+1;
Sigma = diag(v);
%Sigma = eye(p);
%Y = normrnd(10,10,[p,p]);
%Sigma = Y'*Y+100*eye(p);

%generate a random matrix with prescribed covariance matrix
X = normrnd(0,1,[n,p])*sqrtm(Sigma);
%X = rand(n,p)*sqrtm(Sigma);

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
xlim([min(m),max(m)]);
set(gca,'fontsize',20)
legend('Numerical','Theoretical');
str = sprintf( 'p/n=%.3f',p/n);
title(str);

  if savefigs==1
        filename = ...
            sprintf( './Randomly-distributed-diagonal-covariance-n=%d-p=%d-gamma=%.3f.png',...
            n,p,gamma);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        if closefigs==1
            close(gcf)
        end
  end
end