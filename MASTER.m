% MASTER SCRIPT Inference Simulation for HMMs with Preferential Attachment
%   MASTER() Generates simulation results for the 2017 ICASSP paper: 
%   Nonparametric Learning for Hidden Markov Models with Preferential 
%   Dynamics
%
%   MIT License
% 
%   Copyright (c) 2017 Asher A. Hensley Research
%   $Revision: 0.1 $  $Date: 2016/09/06 19:23:54 $
%   
%   Permission is hereby granted, free of charge, to any person obtaining a 
%   copy of this software and associated documentation files (the 
%   "Software"), to deal in the Software without restriction, including 
%   without limitation the rights to use, copy, modify, merge, publish, 
%   distribute, sublicense, and/or sell copies of the Software, and to 
%   permit persons to whom the Software is furnished to do so, subject to 
%   the following conditions:
% 
%   The above copyright notice and this permission notice shall be included 
%   in all copies or substantial portions of the Software.
% 
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
%   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
%   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
%   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
%   CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
%   TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
%   SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%Clean Up
close all
clear
clc

%Init Hyperparameters
a = 1;
b = 1;
c = 1;
d = 1;

%Set Random Seed
seed = 88;
rng(seed);

%Generate Sample Realization
nsamp = 2000;
alpha = 0.75;
[ztrue,n] = makeStateSequence(alpha,nsamp);
lambda = gamrnd(c,1/d,1,length(n));
x = makeObservations(ztrue,lambda);

%Initialize Markov Chain
a0 = 1;
[zi,stats] = initMarkovChain(x,c,d,a0);

%Run Gibbs Sampler
ngibbs = 2000;
nchains = 1;
S = cell(1,nchains);
for kk = 1:nchains
    disp(['Running Chain #' num2str(kk) ' of ' num2str(nchains)])
    S{kk} = runGibbsSampler(x,zi,stats.lambda,a0,ngibbs,a,b,c,d);
end

%Init Results Figure
figure,
set(gcf,'outerposition',[201,110,1006,642],...
    'paperpositionmode','auto')

%Plot Measurements
subplot(2,3,1)
plot(x)
grid on
set(gca,'fontsize',12)
xlabel('Time')
ylabel('Measurement')
title('(A)')

%Plot Number of States
subplot(2,3,2)
muMax = zeros(1,ngibbs);
for kk = 1:nchains
    muMax = muMax+max(S{kk}.z);
end
muMax = muMax/nchains;
plot(muMax,'linewidth',1.5),hold on
ax = axis;
plot(ax(1:2),max(ztrue)*[1,1],'r','linewidth',1.5)
grid on
set(gca,'fontsize',12)
xlabel('Gibbs Sample No.')
ylabel('Number of States')
legend('Inferred','True')
title('(B)')

%Plot Tails
subplot(2,3,3)
muTail = zeros(1,ngibbs);
for kk = 1:nchains
    muTail = muTail+S{kk}.tail;
end
muTail = muTail/nchains;
plot(muTail,'linewidth',1.5),hold on
ax = axis;
plot(ax(1:2),4.55*[1,1],'r','linewidth',1.5)
ylim([0,10])
grid on
set(gca,'fontsize',12)
xlabel('Gibbs Sample No.')
ylabel('% Observations > 2\sigma')
legend('Measured','Theoretical')
title('(C)')

%Plot Learned Sigmas
subplot(2,3,4)
nburn = 1000;
muLambda = zeros(1,nsamp);
for kk = 1:nchains
    muLambda = muLambda+mean(S{kk}.lambda(:,nburn+1:end),2)';
end
muLambda = muLambda/nchains;
plot(muLambda,'linewidth',2),hold on
plot(lambda(ztrue),'r--','linewidth',2)
grid on
set(gca,'fontsize',12)
xlabel('Time')
ylabel('Precision')
legend('Inferred','True')
title('(D)')

%Plot States
clear h
subplot(2,3,5)
muState = zeros(1,nsamp);
for kk = 1:nchains
    plot(S{kk}.z(:,nburn+1:end),'b')
    hold on
    muState = muState+mean(S{kk}.z(:,nburn+1:end),2)';
end
h(1) = plot(S{kk}.z(:,nburn+1));
h(2) = plot(ztrue,'r','linewidth',2);
grid on
set(gca,'fontsize',12)
xlabel('Time')
ylabel('State Number')
legend(h,'Sampled','True','location','northwest')
title('(E)')

%Plot Alpha
subplot(2,3,6)
muAlpha = zeros(1,ngibbs);
for kk = 1:nchains
    muAlpha = muAlpha+S{kk}.alpha;
end
plot(muAlpha/nchains,'b','linewidth',2)
hold on
ax = axis;
plot(ax(1:2),0.75*[1,1],'r','linewidth',2)
grid on
set(gca,'fontsize',12,'box','on')
xlabel('Gibbs Sample No.')
ylabel('Alpha')
legend('Sampled','True')
title('(F)')