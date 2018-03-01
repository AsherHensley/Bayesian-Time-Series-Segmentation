function [z,stats] = initMarkovChain(x,c,d,alpha)
%MAKESTATESEQUENCE Initialize Markov Chain 
%   [z,stats] = initMarkovChain(x,c,d,alpha) returns initialized latent
%   state sequence z and stats structure given measurements x and
%   hyperparameters c, d, and alpha.
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

%Setup
k = 1;
N = length(x);
n = 1;
z = [1,zeros(1,N-1)];
lambda = sampleGammaPosterior(x(1),0,c,d);

stats.lambda = [lambda,zeros(1,N-1)];
stats.nstates = [1,zeros(1,N-1)];

%Run
for t = 2:N

    %Setup Urn
    prior = [n(k)./(n(k)+alpha),alpha./(n(k)+alpha)];
    L0 = normpdf(x(t),0,1/sqrt(lambda));
    L1 = stupdf(x(t),0,c/d,2*c);
    urn = [L0,L1].*prior;
    urn = urn./sum(urn);
    
    %State Transition
    U = rand<urn(1);
    
    %Update
    switch U
        case 1  %Same State
            z(t) = k;
            n(k) = n(k)+1;
            lambda = sampleGammaPosterior(x(z==k),0,c,d);
        case 0  %New state
            k = k+1;
            z(t) = k;
            n = [n,1];
            lambda = sampleGammaPosterior(x(t),0,c,d);
    end
    
    %Update Stats
    stats.nstates(t) = length(n);
    stats.lambda(t) = lambda;

end

%Likelihood
L = length(n);
stats.LL = L*log(alpha)+sum(log(beta(n,alpha+1)))...
    +sum(log(normpdf(x,0,1./sqrt(stats.lambda))));
