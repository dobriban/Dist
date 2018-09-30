%Experiments of ARE for Elliptical models, the worst case.
cd('C:\Dropbox\Projects\Distributed\Experiments\\ARE for Elliptical models')
%%
n = 10000;
epsilon = 1;

alpha_arr = [100,1000,10000];

for L = 1:length(alpha_arr)
alpha = alpha_arr(L);

%c-greater-than-2gamma
c_arr = [0.05, 0.1, 0.5, 0.8];
%create a test covariance matrix
for L1 = 1:length(c_arr)
c = c_arr(L1);  
p = 100;
gamma = p/n;

%create a test covariance matrix
rng(2);
savefigs=1;    closefigs=1;

v = binornd(1,c,[n,1]);
for i = 1:n
    if v(i)==0
        v(i)=v(i)+1;
    else
        v(i)=alpha*v(i);
    end
end

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
myfun = @(x,b) (1-c)/(1+epsilon*x)+c/(1+alpha*epsilon*x)-b;  % parameterized function
b = 1-gamma;                    % parameter
fun = @(x) myfun(x,b); % function of x alone
phi = fzero(fun,[0.00000000001 10000000000]);

for i = 1:(1/gamma-1)
   
myfun = @(x,b) (1-c)/(1+epsilon*x)+c/(1+alpha*epsilon*x)-b;  % parameterized function
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
xlim([min(m),max(m)]);
set(gca,'fontsize',20);
legend('Numerical', 'Theoretical','Optimal');
str1=sprintf('p/n=%.3f, ',p/n);
str2=sprintf('alpha=%d, ',alpha);
str3=sprintf('c=%.3f',c);
s=[str1 str2 str3];
title(s);

  if savefigs==1
        filename = ...
            sprintf( './Elliptical-WorstCase-c-greater-than-2gamma-n=%d-p=%d-gamma=%.3f-alpha=%d-c=%.3f.png',...
            n,p,gamma,alpha,c);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        if closefigs==1
            close(gcf)
        end
  end
end

%0-less-than-c-less-than-gamma
c_arr = [0.001, 0.005];
%create a test covariance matrix
for L2 = 1:length(c_arr)
c = c_arr(L2);  
p = 100;
gamma = p/n;

%create a test covariance matrix
rng(2);
savefigs=1;    closefigs=1;

v = binornd(1,c,[n,1]);
for i = 1:n
    if v(i)==0
        v(i)=v(i)+1;
    else
        v(i)=alpha*v(i);
    end
end

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
myfun = @(x,b) (1-c)/(1+epsilon*x)+c/(1+alpha*epsilon*x)-b;  % parameterized function
b = 1-gamma;                    % parameter
fun = @(x) myfun(x,b); % function of x alone
phi = fzero(fun,[0.00000000001 10000000000]);

for i = 1:(1/gamma-1)
   
myfun = @(x,b) (1-c)/(1+epsilon*x)+c/(1+alpha*epsilon*x)-b;  % parameterized function
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
xlim([min(m),max(m)]);
set(gca,'fontsize',20);
legend('Numerical', 'Theoretical','Optimal');
str1=sprintf('p/n=%.3f, ',p/n);
str2=sprintf('alpha=%d, ',alpha);
str3=sprintf('c=%.3f',c);
s=[str1 str2 str3];
title(s);

  if savefigs==1
        filename = ...
            sprintf( './Elliptical-WorstCase-0-less-than-c-less-than-gamma-n=%d-p=%d-gamma=%.3f-alpha=%d-c=%.3f.png',...
            n,p,gamma,alpha,c);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        if closefigs==1
            close(gcf)
        end
  end
end

%gamma-less-than-c-less-than-2gamma
c = 0.015;
%create a test covariance matrix
p = 100;
gamma = p/n;

%create a test covariance matrix
rng(2);
savefigs=1;    closefigs=1;

v = binornd(1,c,[n,1]);
for i = 1:n
    if v(i)==0
        v(i)=v(i)+1;
    else
        v(i)=alpha*v(i);
    end
end

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
myfun = @(x,b) (1-c)/(1+epsilon*x)+c/(1+alpha*epsilon*x)-b;  % parameterized function
b = 1-gamma;                    % parameter
fun = @(x) myfun(x,b); % function of x alone
phi = fzero(fun,[0.00000000001 10000000000]);

for i = 1:(1/gamma-1)
   
myfun = @(x,b) (1-c)/(1+epsilon*x)+c/(1+alpha*epsilon*x)-b;  % parameterized function
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
xlim([min(m),max(m)]);
set(gca,'fontsize',20);
legend('Numerical', 'Theoretical','Optimal');
str1=sprintf('p/n=%.3f, ',p/n);
str2=sprintf('alpha=%d, ',alpha);
str3=sprintf('c=%.3f',c);
s=[str1 str2 str3];
title(s);

  if savefigs==1
        filename = ...
            sprintf( './Elliptical-WorstCase-gamma-less-than-c-less-than-2gamma-n=%d-p=%d-gamma=%.3f-alpha=%d-c=%.3f.png',...
            n,p,gamma,alpha,c);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        if closefigs==1
            close(gcf)
        end
  end




%c=gamma
c = 0.01;
%create a test covariance matrix
p = 100;
gamma = p/n;

%create a test covariance matrix
rng(2);
savefigs=1;    closefigs=1;

v = binornd(1,c,[n,1]);
for i = 1:n
    if v(i)==0
        v(i)=v(i)+1;
    else
        v(i)=alpha*v(i);
    end
end

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
myfun = @(x,b) (1-c)/(1+epsilon*x)+c/(1+alpha*epsilon*x)-b;  % parameterized function
b = 1-gamma;                    % parameter
fun = @(x) myfun(x,b); % function of x alone
phi = fzero(fun,[0.00000000001 10000000000]);

for i = 1:(1/gamma-1)
   
myfun = @(x,b) (1-c)/(1+epsilon*x)+c/(1+alpha*epsilon*x)-b;  % parameterized function
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
xlim([min(m),max(m)]);
set(gca,'fontsize',20);
legend('Numerical', 'Theoretical','Optimal');
str1=sprintf('p/n=%.3f, ',p/n);
str2=sprintf('alpha=%d, ',alpha);
str3=sprintf('c=%.3f',c);
s=[str1 str2 str3];
title(s);

  if savefigs==1
        filename = ...
            sprintf( './Elliptical-WorstCase-c=gamma-n=%d-p=%d-gamma=%.3f-alpha=%d-c=%.3f.png',...
            n,p,gamma,alpha,c);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        if closefigs==1
            close(gcf)
        end
  end



%c=2gamma
c = 0.02;
%create a test covariance matrix
p = 100;
gamma = p/n;

%create a test covariance matrix
rng(2);
savefigs=1;    closefigs=1;

v = binornd(1,c,[n,1]);
for i = 1:n
    if v(i)==0
        v(i)=v(i)+1;
    else
        v(i)=alpha*v(i);
    end
end

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
myfun = @(x,b) (1-c)/(1+epsilon*x)+c/(1+alpha*epsilon*x)-b;  % parameterized function
b = 1-gamma;                    % parameter
fun = @(x) myfun(x,b); % function of x alone
phi = fzero(fun,[0.00000000001 10000000000]);

for i = 1:(1/gamma-1)
   
myfun = @(x,b) (1-c)/(1+epsilon*x)+c/(1+alpha*epsilon*x)-b;  % parameterized function
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
xlim([min(m),max(m)]);
set(gca,'fontsize',20);
legend('Numerical', 'Theoretical','Optimal');
str1=sprintf('p/n=%.3f, ',p/n);
str2=sprintf('alpha=%d, ',alpha);
str3=sprintf('c=%.3f',c);
s=[str1 str2 str3];
title(s);

  if savefigs==1
        filename = ...
            sprintf( './Elliptical-WorstCase-c=2gamma-n=%d-p=%d-gamma=%.3f-alpha=%d-c=%.3f.png',...
            n,p,gamma,alpha,c);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        if closefigs==1
            close(gcf)
        end
  end

end