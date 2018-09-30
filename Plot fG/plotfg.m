%Plots of f,G
cd('C:\Dropbox\Projects\Distributed\Experiments\Plot fG')


%rng(2)
x = [0:0.1:100];
y = 1./(1+x);
plot(x,y,'LineWidth',3)
xlabel('x');
ylabel('\eta(x)');
set(gca,'fontsize',20);
xlim([min(x),max(x)]);

savefigs = 1;
  if savefigs==1
        filename = sprintf( './eta-plot.png');
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
  end


%rng(2)
x = [0:0.001:1];
y = x./(1-x);
plot(x,y,'LineWidth',3)
xlabel('\gamma');
ylabel('f(\gamma)');
set(gca,'fontsize',20);
xlim([min(x),max(x)]);

savefigs = 1;
  if savefigs==1
        filename = sprintf( './f-plot.png');
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
  end


