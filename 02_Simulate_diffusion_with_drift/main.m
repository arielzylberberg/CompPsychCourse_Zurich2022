addpath(genpath('../functions_addtopath/'));

%%

seed = 104046;
rng(seed,'twister');

%% some figure defaults
set(0,'defaultAxesFontSize',18);
set(0, 'DefaultLineLineWidth', 1);

%%  parameters
NumTrials = 5000; % number of simulated trials
NumTimeSteps = 1000; % number of time steps
DeltaT = 0.001; % time step

sigma = 1; % standard deviation of the DV after 1 sec of accumulation
mu = 5;    % mean of the DV after 1 sec of accumulation

%%  momentary evidence
mean_per_step = mu * DeltaT;
stdev_per_step = sigma * sqrt(DeltaT);
e = mean_per_step + stdev_per_step * randn(NumTrials, NumTimeSteps);

%% accumulated evidence
dv = cumsum(e,2);

%% FIGURES

%% plot the DV one trial
figure();
t = [1:NumTimeSteps] * DeltaT;
plot(t, dv(1,:));
xlabel('Time [s]');
ylabel('Decision variable');

%% plot the DV for a few trials
figure();
t = [1:NumTimeSteps] * DeltaT;
plot(t, dv(1:10,:));
xlabel('Time [s]');
ylabel('Decision variable');

%% mean and variance of the DV as a function of time
figure()
subplot(1,2,1);
plot(t, mean(dv))
xlabel('Time [s]')
ylabel('E[DV]')

subplot(1,2,2);
plot(t, std(dv).^2)
xlabel('Time [s]')
ylabel('VAR [DV]')

%% histogram of the state of the DV at different times
tind = findclose(t,[0.1,0.5,1]);
n = length(tind);
p = publish_plot(n,1);
set(gcf,'Position',[505  229  390  468]);
edges = linspace(-9,9,200);
for i=1:n
    p.next();
    histogram(dv(:,tind(i)),edges);
    de = edges(2)-edges(1);
    pr = normpdf(edges,mu*t(tind(i)),sigma*sqrt(t(tind(i))));
    pr = pr*de*NumTrials;
    hold all
    plot(edges,pr,'r-','LineWidth',1);
%     str = ['t = ',num2str(t(tind(i))),' s'];
%     ht(i) = p.text_draw(i,str);
    ylabel('Trial count');
    xlabel('Accumulated evidence (a.u.)');
end

p.format('FontSize',15);
% set(ht,'FontSize',15);
same_ylim(p.h_ax);
% p.unlabel_center_plots();

%% with bounds
% set simulation parameters
% recalculate the dv, with a lower mu to have more errors
mu = 1;    % mean of the DV after 1 sec of accumulation
sigma = 1; 
NumTrials = 5000; % number of simulated trials
NumTimeSteps = 4000; % number of time steps
DeltaT = 0.001; % time step
t = [1:NumTimeSteps] * DeltaT;

%%
mean_per_step = mu * DeltaT;
stdev_per_step = sigma * sqrt(DeltaT);
e = mean_per_step + stdev_per_step * randn(NumTrials, NumTimeSteps);
dv = cumsum(e,2);

Bup = 1;
[choice, DecT, DVnan] = bound_cross_flat(t, dv, Bup);

%% plots
p = publish_plot(1,1);
plot(t, DVnan(1:100,:)');
xlabel('Time [s]');
ylabel('Accumulated evidence');
p.format();

p = publish_plot(1,1);
histogram(DVnan(:,500),100);
xlabel('Decision variable (a.u.)');
ylabel('# trials');
p.format();

p = publish_plot(1,1);
[tt,xx,ss] = curva_media(DecT, choice, [],0);
terrorbar(tt,xx,ss);
xlabel('Accuracy');
ylabel('Decime time [s]');
xlim([-0.1,1.1]);
p.format();

p = publish_plot(1,1);
histogram(DecT,100);
ylabel('# trials');
xlabel('Decision Time [s]');
p.format();


