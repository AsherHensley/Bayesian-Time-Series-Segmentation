function S = runGibbsSampler(x,z0,lambda0,alpha0,ngibbs,a,b,c,d)
%RUNGIBBSSAMPLER Run Gibbs Sampling Inference Procedure
%   S = runGibbsSampler(x,z0,lambda0,alpha0,ngibbs,a,b,c,d) applies
%   Gibbs sampling inference procedure and returns sampling structure S. 
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

%Initialize
nsamp = length(x);
n0 = hist(z0,unique(z0));

%Setup MCMC Struct
gibbs.a = a;
gibbs.b = b;
gibbs.c = c;
gibbs.d = d;
gibbs.z = z0;
gibbs.n = n0;
gibbs.alpha = alpha0;
gibbs.lambda = lambda0;

%Samples
S.z = zeros(nsamp,ngibbs);
S.lambda = zeros(nsamp,ngibbs);
S.alpha = nan(1,ngibbs);
S.tail = nan(1,ngibbs);

%Run Sampler
hw = waitbar(0,'Sample 0');
dT = 0;
for ii = 1:ngibbs
    
    %Sample State Sequence
    tic
    
    %First Point
    t = 1;
    gibbs = sampleFirstPoint(x,t,gibbs);
    
    %Points 2:nsamp-1
    t = 2;
    while t<nsamp
        
        %Select Sampling Case
        leftEdge = gibbs.z(t)==gibbs.z(t-1);
        rightEdge = gibbs.z(t)==gibbs.z(t+1);
        type = [num2str(leftEdge),num2str(rightEdge)];
        
        %Sample Current Assignment
        switch type
            case '00'
                gibbs = sampleDoubleBoundary(x,t,gibbs);
                t = t+1;
                
            case '01'
                gibbs = sampleLeftBoundary(x,t,gibbs);
                t = t+1;
                
            case '10'
                gibbs = sampleRightBoundary(x,t,gibbs);
                t = t+1;
                
            case '11'
                t = t+1;
                continue
                
        end
    end
    
    %Last Point
    gibbs.z(end) = gibbs.z(end-1);

    %Sample Alpha
    gibbs.alpha = sampleAlphaGibbs(gibbs.alpha,gibbs.n,gibbs.a,gibbs.b);
    
    %Sample Parameters
    nstates = length(gibbs.n);
    gibbs.lambda = zeros(1,nstates);
    for kk = 1:nstates
        mask = gibbs.z==kk;
        gibbs.lambda(kk) = sampleGammaPosterior(x(mask),0,gibbs.c,gibbs.d);
    end
    
    %Save Samples
    S.z(:,ii) = gibbs.z;
    S.lambda(:,ii) = gibbs.lambda(gibbs.z);
    S.alpha(ii) = gibbs.alpha;
    S.tail(ii) = 100*sum(abs(x)>2./sqrt(gibbs.lambda(gibbs.z)))/nsamp;
    
    %Update Waitbar
    dt = toc;
    dT = dT+dt;
    waitbar(ii/ngibbs,hw,['Sample ' num2str(ii) ' of '...
        num2str(ngibbs) ' Complete (Time Remaining: ',...
        num2str((ngibbs-ii)*dT/ii/60) ' minutes)'])
end
delete(hw)