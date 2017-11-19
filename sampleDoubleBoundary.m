function gibbs = sampleDoubleBoundary(x,t,gibbs)
%SAMPLEDOUBLEBOUNDARY Sample Double Boundary Point of Latent State Sequence
%   gibbs = sampleDoubleBoundary(x,t,gibbs) returns updated gibbs struct 
%   given measurements x, time stamp t, and current gibbs struct.
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

%Unpack
j = gibbs.z(t);
n = gibbs.n;
A = gibbs.alpha;
L = gibbs.lambda;
c = gibbs.c;
d = gibbs.d;

%Compute Weights
moveLeft = n(j-1)/(n(j-1)+A+1)*normpdf(x(t),0,sqrt(1/L(j-1)));
moveRight = n(j+1)/(n(j+1)+A+1)*normpdf(x(t),0,sqrt(1/L(j+1)));
splitEvt = A/(1+A)*stupdf(x(t),0,c/d,2*c);

%Normalize
w = [moveLeft,moveRight,splitEvt];
w = w/sum(w);

%Draw Sample
cdf = cumsum(w);
U = rand;

%Update
if U<cdf(1)
    
    %Assign to Left State
    gibbs.z(t) = j-1;
    gibbs.z(t+1:end) = gibbs.z(t+1:end)-1;
    gibbs.lambda(j) = [];
    
elseif U>=cdf(1) && U<cdf(2)
    
    %Asign to Right State
    gibbs.z(t) = j;
    gibbs.z(t+1:end) = gibbs.z(t+1:end)-1;
    gibbs.lambda(j) = [];
    
else
    
    %Split
    gibbs.z(t) = j;
    
    %Sample Precision for New State
    gibbs.lambda(j) = sampleGammaPosterior(x(t),0,c,d);
    
end

%Update n
gibbs.n = hist(gibbs.z,unique(gibbs.z));