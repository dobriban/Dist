n = 10000;
p = 100;
gamma = p/n;
alpha = 10000;

c = 1.5*gamma;
g_t = [1,alpha];
savefigs=1;    closefigs=1;

for k=1:length(g_t)
    
    rng(2);

ARE = zeros(1,1/gamma);

myfun = @(x,b) (1-c)/(1+x)+c/(1+alpha*x)-b;  % parameterized function
b = 1-gamma;                    % parameter
fun = @(x) myfun(x,b); % function of x alone
phi = fzero(fun,[0.00000000001 10000000000]);

for i = 1:(1/gamma-1)
   
myfun = @(x,b) (1-c)/(1+x)+c/(1+alpha*x)-b;  % parameterized function
b = 1-i*gamma;                    % parameter
fun = @(x) myfun(x,b); % function of x alone
x = fzero(fun,[0.00000000001 10000000000]);

ARE(i) = i*phi/x;
end
ARE(1/gamma)=0;



AIE = zeros(1,1/gamma);

for i = 1:(1/gamma-1)
   
myfun = @(x,b) (1-c)/(1+x)+c/(1+alpha*x)-b;  % parameterized function
b = 1-i*gamma;                    % parameter
fun = @(x) myfun(x,b); % function of x alone
x = fzero(fun,[0.00000000001 10000000000]);

AIE(i) = (i*i*(1-gamma))/(i*gamma+(1+(alpha-1)*c)*(i-1)*x+i*i*(1-2*gamma));
end
AIE(1/gamma)=0;

AOE = zeros(1,1/gamma);


myfun = @(x,b) (1-c)/(1+x)+c/(1+alpha*x)-b;  % parameterized function
b = 1-gamma;                    % parameter
fun = @(x) myfun(x,b); % function of x alone
psi = fzero(fun,[0.00000000001 10000000000]);

for i = 1:(1/gamma-1)
   
myfun = @(x,b) (1-c)/(1+x)+c/(1+alpha*x)-b;  % parameterized function
b = 1-i*gamma;                    % parameter
fun = @(x) myfun(x,b); % function of x alone
x = fzero(fun,[0.00000000001 10000000000]);

AOE(i) = (i+i*g_t(k)*psi)/(i+g_t(k)*x);
end
AOE(1/gamma)=0;

m = 1:1:1/gamma;
plot(m,ARE,'-','linewidth',4,'color',rand(1,3));
hold on 
plot(m,AIE,'--','linewidth',4,'color',rand(1,3));
plot(m,AOE,'-.','linewidth',4,'color',rand(1,3));
hold off
xlabel('Number of Machines');
ylabel('Efficiency');
set(gca,'fontsize',20)
legend('Estimation','Train','Test');
str = sprintf( 'g_t=%d,c=%.3f,p/n=%.3f',g_t(k),c,gamma);
title(str);

    if savefigs==1
        filename = sprintf( './Elliptical-g_t=%d-c=%.3f-gamma=%.3f.png',g_t(k),c,gamma);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        if closefigs==1
            close(gcf)
        end
    end

end
