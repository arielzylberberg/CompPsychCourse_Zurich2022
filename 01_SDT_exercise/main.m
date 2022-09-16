addpath(genpath('../functions_addtopath/'));

%% parameters
sigma = 10;
mu = [25,35];

%% with simulations
ntrials = 100000;
s = randn(ntrials,2);
s = s * sigma + mu;

csim = mean(s(:,2)>s(:,1));

%% without simulations
sigma_dif = sqrt(2)*sigma;
mu_dif = diff(mu);
c = 1 - normcdf(0,mu_dif,sigma_dif);

%% plot
p = publish_plot(1,1);
set(gcf,'Position',[224  246  730  283])
x = linspace(-60,60,1000);
y = normpdf(x,mu_dif,sigma_dif);
plot(x,y,'k');
hold all
ha = area(x(x<0),y(x<0));
set(ha,'EdgeColor','none','FaceColor',0.5*[1,1,1]);
set(p.h_ax(1),'ycolor','none');
xlabel('Difference in perceived line length [cm]')
p.format('FontSize',25,'LineWidthPlot',2);

%% probability correct vs delta length

MU_DIF = linspace(0,50,100);
C = 1-normcdf(0,MU_DIF,sigma_dif);

p = publish_plot(1,1);
plot(MU_DIF, C);
ylabel('Probability correct');
xlabel('\Delta line length [cm]');
p.format('presentation')