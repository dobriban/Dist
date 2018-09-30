%Experiments of ARE for MP models.
cd('C:\Dropbox\Projects\Distributed\Experiments\Plot Eff')
%%
gamma_array = [0.01,0.02,0.05,0.1,0.2,0.3];
a = {'-','--','-.',':'};
savefigs=1;    closefigs=1;

for i=1:length(gamma_array)
    figure, hold on
    rng(2);
    
    gamma  = gamma_array(i);
    AE  = @(k) (1-k.*gamma)/(1-gamma);
    [x,y]= fplot(AE,[1,1/gamma]);
    h1=plot(x,y,'linewidth',4,'color',rand(1,3));
    set(h1,'LineStyle',a{1});
    
    OE  = @(k) 1./(1+(k-1).*gamma^2./(1-k.*gamma));
    [x,y]= fplot(OE,[1,1/gamma]);
    h2=plot(x,y,'linewidth',4,'color',rand(1,3));
    set(h2,'LineStyle',a{2});
       
    xlabel('Number of Machines');
    ylabel('Efficiency');
    xlim([1,1/gamma]);
    set(gca,'fontsize',20)
    legend('Estimation Error','Test Error');
    str = sprintf( 'p/n=%.3f',gamma);
    title(str);
    
    if savefigs==1
        filename = sprintf( './E-MP-gamma=%.3f.png',gamma);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        if closefigs==1
            close(gcf)
        end
    end
end
