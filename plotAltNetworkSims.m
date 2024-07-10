% script for simulating networks with other connectivity and
% nonlinearities, and plotting the results shown in ED Fig 5

%% 1d simulations

clear

p = [];             %parameter structure
p.n = 8;            %number of neurons
p.tau = 1;          %single neuron time constant
p.T = 500;          %simulation duration
p.J0 = -30;         %inhibition strength
p.I = 10;           %external input 
p.j1v = 32:.1:37;   %array of excitation strengths
p.ntrial = 15;      %number of sampled angles for the dynamics
p.flag_lin = 0;     % 1: threshold-linear transfers function; 0:logexp

% estimate circular variance for varying excitation strength
v = zeros(size(p.j1v));
for i = 1:length(p.j1v)
    % disp(i/length(p.j1v))
    p.J1 = p.j1v(i);
    % initial (psi0) and final (psi1) angles
    [psi0, psi1] = run_test(p); 
    delta_angle = psi1 - psi0;
    tmp = sum(sin(delta_angle)).^2 + sum(cos(delta_angle)).^2;
    v(i) = 1 - sqrt(tmp)/length(delta_angle); %circular variance
end
save alternateNetworkSimulations/scan_1d_500timesteps_circvar_v2.mat p v

% change inhibition strength
p.J0 = -20; 

% estimate circular variance for varying excitation strength
v = zeros(size(p.j1v));
for i = 1:length(p.j1v)
    % disp(i/length(p.j1v))
    p.J1 = p.j1v(i);
    [psi0, psi1] = run_test(p);
    delta_angle = psi1 - psi0;
    tmp = sum(sin(delta_angle)).^2 + sum(cos(delta_angle)).^2;
    v(i) = 1 - sqrt(tmp)/length(delta_angle); %circular_variance
end
save alternateNetworkSimulations/scan_1d_500timesteps_circvar_JI20_v2.mat p v

%% 2d simulations

p = [];                 %parameter structure
p.n = 4^2;              %number of neurons
p.tau = 1;              %single neuron time constant
p.I = 10;               %external input 
p.T = 500;              %simulation duration
p.J0 = -20;             %inhibition strength
p.ntrial = 6^2;         %number of sampled angles for the dynamics
p.flag_lin = 1;         % 1: threshold-linear transfers function; 0:logexp
p.j1v = 15.9:.01:16.1;  %array of excitation strengths

v = zeros(size(p.j1v));
m = v;
for i = 1:length(p.j1v)
    % disp(i/length(p.j1v))
    p.J1 = p.j1v(i);
    % initial (psi0_1,psi0_2) and final (psi1_1,psi1_2) angles
    [psi0_1, psi0_2, psi1_1, psi1_2] = run_test_2d(p);
    delta_angle_1 = psi1_1 - psi0_1;
    delta_angle_2 = psi1_2 - psi0_2;
    tmp = (sum(sin(delta_angle_1)).^2 + sum(cos(delta_angle_1)).^2+...
        sum(sin(delta_angle_2)).^2 + sum(cos(delta_angle_2)).^2)/2;
    v(i) = 1 - sqrt(tmp)/length(delta_angle_1); %variance
end

save alternateNetworkSimulations/scan_2d_500timesteps_circvar_v2.mat p v
%% plotting 1d simulations
clear
fs = 12; %fontsize

figure(1)
clf
subplot(2,3,[1 2 3])
load alternateNetworkSimulations/scan_1d_500timesteps_circvar_JI20_v2.mat
scatter(p.j1v,v,30,[0 0 0])
load alternateNetworkSimulations/scan_1d_500timesteps_circvar_v2.mat
hold on
scatter(p.j1v,v,80,[0 0 0],'filled')
scatter(p.j1v(15),v(15),150,[1 0 0],'filled')
scatter(p.j1v(28),v(28),150,[0 0 1],'filled')
scatter(p.j1v(40),v(40),150,[0 1 0],'filled')
set(gca,'Box','off','TickDir','out','FontSize',fs)
xlim([33 36])
xlabel('J_E')
ylabel('circular variance')

%panels for specific values of excitation strength
p.ntrial = 50; % #of initial angles
p.J1 = p.j1v(15); 
[psi0, psi1] = run_test(p); % initial (psi0) and final (psi1) angles
subplot(2,3,4)
scatter(psi0/2/pi*360,psi1/2/pi*360,40,[1 0 0],'filled')
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-180 -90 0 90 180])
set(gca,'Box','off','TickDir','out','FontSize',fs)
axis square
xlabel('initial')
ylabel('final')

p.J1 = p.j1v(28);
[psi0, psi1] = run_test(p);
subplot(2,3,5)
scatter(psi0/2/pi*360,psi1/2/pi*360,40,[0 0 1],'filled')
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-180 -90 0 90 180])
set(gca,'Box','off','TickDir','out','FontSize',fs)
axis square
xlabel('initial')
ylabel('final')

p.J1 = p.j1v(40);
[psi0, psi1] = run_test(p);
subplot(2,3,6)
scatter(psi0/2/pi*360,psi1/2/pi*360,40,[0 1 0],'filled')
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-180 -90 0 90 180])
set(gca,'Box','off','TickDir','out','FontSize',fs)
axis square
xlabel('initial')
ylabel('final')

fig = gcf;
fig.Theme = "light"; 

%% plotting 2d simulations
clear
fs = 12; %font size
load alternateNetworkSimulations/scan_2d_500timesteps_circvar_v2.mat

figure(2)
clf
subplot(2,3,[1 2 3])
scatter(p.j1v(1:1:end),v(1:1:end),80,[0 0 0],'filled')
hold on
scatter(p.j1v(3),v(3),150,[1 0 0],'filled')
scatter(p.j1v(13),v(13),150,[0 0 1],'filled')
scatter(p.j1v(18),v(18),150,[0 1 0],'filled')
set(gca,'Box','off','TickDir','out','FontSize',fs)
xlim([p.j1v(1) p.j1v(end)])
xlim([15.9 16.1])
xlabel('J_E')
ylabel('circular variance')

%panels for specific values of excitation strength
p.ntrial = 10^2; % #of initial angles
p.J1 = p.j1v(3);
[psi1, psi2, d1, d2] = run_test_2d(p);
subplot(2,3,4)
scatter(psi1/2/pi*360,psi2/2/pi*360,50,[0 0 0],'LineWidth',2)
hold on
scatter(d1/2/pi*360,...
    d2/2/pi*360,20,[1 0 0],'filled')
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-180 -90 0 90 180])
set(gca,'Box','off','TickDir','out','FontSize',fs)
axis square
xlabel('angle 1')
ylabel('angle 2')

p.J1 = p.j1v(13);
[psi1, psi2, d1, d2] = run_test_2d(p);
subplot(2,3,5)
scatter(psi1/2/pi*360,psi2/2/pi*360,50,[0 0 0],'LineWidth',2)
hold on
scatter(d1/2/pi*360,...
    d2/2/pi*360,20,[0 0 1],'filled')
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-180 -90 0 90 180])
set(gca,'Box','off','TickDir','out','FontSize',fs)
axis square
xlabel('angle 1')
ylabel('angle 2')

p.J1 = p.j1v(18);
[psi1, psi2, d1, d2] = run_test_2d(p);
subplot(2,3,6)
scatter(psi1/2/pi*360,psi2/2/pi*360,50,[0 0 0],'LineWidth',2)
hold on
scatter(d1/2/pi*360,...
    d2/2/pi*360,20,[0 1 0],'filled')
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-180 -90 0 90 180])
set(gca,'Box','off','TickDir','out','FontSize',fs)
axis square
xlabel('angle 1')
ylabel('angle 2')

fig = gcf;
fig.Theme = "light"; 

%% functions
function [initial_angle, final_angle] = run_test(p)
    if (p.flag_lin)
        p.io = @(c,p) linthr(c,p);
    else
        p.io = @(c,p) logexp(c,p);
    end

    %weight matix
    p.th = linspace(-pi,pi,p.n+1); 
    p.th(end) = [];
    p.W = (p.J1*exp(2*cos(p.th-p.th'))/besseli(0,2)/2/pi+p.J0)/p.n; 

    final_angle = zeros(p.ntrial, 1); 
    initial_angle = zeros(p.ntrial, 1);
    for k = 1:p.ntrial
        %dynamics
        initial_angle(k) = 2*pi*(k-1)/p.ntrial-pi;
        p.incond = 5*cos(p.th-initial_angle(k));
        opts = odeset('RelTol',1e-3,'AbsTol',1e-6);
        [~,y] = ode45(@diff_ring, [0,p.T], p.incond, opts, p);
        z = exp(1i*(p.th))*y'/p.n;
        psi = angle(z);
        final_angle(k) = psi(end);
    end
end

function [psi1, psi2, d1, d2] = run_test_2d(p)
    if (p.flag_lin)
        p.io = @(c,p) linthr(c,p);
    else
        p.io = @(c,p) logexp(c,p);
    end

    %weight matrix
    p.th = linspace(-pi,pi,sqrt(p.n)+1); 
    p.th(end) = [];
    [th1,th2]=meshgrid(p.th);
    p.th1 = th1(:)'; p.th2 = th2(:)';
    tmp = ( cos(p.th1-p.th1') + cos(p.th2-p.th2') )/2;
    p.W = (p.J1*tmp+p.J0)/p.n;
    
    psi = linspace(-pi,pi,sqrt(p.ntrial)+1); psi(end) = [];
    [psi1,psi2]=meshgrid(psi);
    psi1 = psi1(:); psi2 = psi2(:);
    d1 = zeros(p.ntrial,1); d2 = d1; 
    for k = 1:p.ntrial 
        %dynamics
        p.incond = (cos(p.th1-psi1(k))+cos(p.th2-psi2(k)))/2;
        [~,y] = ode45(@diff_ring, [0,p.T], p.incond, [], p);
        z1 = exp(1i*(p.th1))*y'/p.n;
        z2 = exp(1i*(p.th2))*y'/p.n;
        ph1 = angle(z1);
        ph2 = angle(z2);
        d1(k) = ph1(end);
        d2(k) = ph2(end);    
    end
end

function dhdt = diff_ring(~,h,p)
dhdt = (-h + p.W*p.io(h,p) + p.I)/p.tau;
end

function r = logexp(c,~)
r = log(1+exp(c));
r(isinf(r))=c(isinf(r));
end

function r = linthr(c,~)
r = c; r(c<0) = 0;
end

