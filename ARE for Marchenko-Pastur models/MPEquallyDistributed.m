%Experiments of ARE for MP models, equally distributed case.
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
%Y = normrnd(10,10,[p,p]);
%Sigma = Y'*Y+100*eye(p);

%generate a random matrix with prescribed covariance matrix
X = normrnd(0,1,[n,p])*sqrtm(Sigma);

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
            sprintf( './MP-Equally-distributed-diagonal-covariance-n=%d-p=%d-gamma=%.3f.png',...
            n,p,gamma);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        if closefigs==1
            close(gcf)
        end
  end
end
