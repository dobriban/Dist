%Experiments of ARE for Elliptical models, uniformly distributed diagonal entries.
cd('C:\Dropbox\Projects\Distributed\Experiments\ARE for Elliptical models')
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
v = rand(n,1);
Sigma1 = eye(p);
%Sigma1 = diag(1000*rand(p,1));
Sigma2 = diag(v);
%Y = normrnd(0,1,[n,n]);
%Sigma2 = Y'*Y+eye(n);
%Sigma2 = full(gallery('tridiag',n,-2,5,-2));

%generate a random matrix with prescribed covariance matrix
X = sqrtm(Sigma2)*normrnd(0,1,[n,p])*sqrtm(Sigma1);

%intialize ARE for each number of machine
ARE1 = zeros(1,1/gamma);

ARE1(1)=1;
for i = 2:1/gamma
    t = 1/trace(inv((X(1:floor(n/i),:))'*X(1:floor(n/i),:)));
    for j = 1:i-2
        X_j = X(j*floor(n/i)+1:(j+1)*floor(n/i),:);
        t = t+1/trace(inv(X_j'*X_j));
    end
    t = t+1/trace(inv((X((i-1)*floor(n/i)+1:n,:))'*X((i-1)*floor(n/i)+1:n,:)));
    ARE1(i)=trace(inv(X'*X))*t;
end

ARE2 = zeros(1,1/gamma);
myfun = @(x,b) log(1+x)/x-b;  % parameterized function
b = 1-gamma;                    % parameter
fun = @(x) myfun(x,b); % function of x alone
phi = fzero(fun,[0.00000000001 10000000000]);

for i = 1:(1/gamma-1)
   
myfun = @(x,b) log(1+x)/x-b;  % parameterized function
b = 1-i*gamma;                    % parameter
fun = @(x) myfun(x,b); % function of x alone
x = fzero(fun,[0.00000000001 10000000000]);

ARE2(i) = i*phi/x;
end
ARE2(1/gamma)=0;


m = 1:1:1/gamma;
plot(m,ARE1,'LineWidth',3)
hold on 
plot(m,ARE2,'--','LineWidth',3)
ARE3 = (1/gamma-m)/(1/gamma-1);
plot(m,ARE3,':','LineWidth',3)
hold off
xlabel('Number of Machines');
ylabel('ARE');
set(gca,'fontsize',20);
xlim([min(m),max(m)]);
legend('Numerical', 'Theoretical','Optimal');
str=sprintf('p/n=%.3f' ,p/n);
title(str);

  if savefigs==1
        filename = ...
            sprintf( './Elliptical-uniformly-distributed-covariance-n=%d-p=%d-gamma=%.3f.png',...
            n,p,gamma);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        if closefigs==1
            close(gcf)
        end
  end
end
