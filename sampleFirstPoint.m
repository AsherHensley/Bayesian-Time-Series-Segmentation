function gibbs = sampleFirstPoint(x,t,gibbs)
%SAMPLEFIRSTPOINT Sample First Point of Latent State Sequence
%   gibbs = sampleFirstPoint(x,t,gibbs) returns updated gibbs struct given
%   measurements x, time stamp t, and current gibbs struct.
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

%Error Chk
if j~=1
    error('Indexing Error at First Point')
end

%Determine Edge Case
rightEdge = gibbs.z(t+1)~=1;

%Handle Cases
if rightEdge
    
    %Weights
    moveRight = n(j+1)/(n(j+1)+A+1)*normpdf(x(t),0,sqrt(1/L(j+1)));
    splitEvt = A/(1+A)*stupdf(x(t),0,c/d,2*c);
    
    %Normalize
    w = [moveRight,splitEvt];
    w = w/sum(w);

    %Update
    U = rand;
    if U<w(1)
        
        %Assign to Right State
        gibbs.z(t) = j+1;
        gibbs.z = gibbs.z-1;
        gibbs.lambda(j) = [];
        
    else
        
        %Split
        gibbs.z(t) = j;
        
        %Sample Precision for New State
        gibbs.lambda(j) = sampleGammaPosterior(x(t),0,c,d);
        
    end
    
else
    
    %Weights
    noChange = (n(j)-1)/(n(j)+A)*normpdf(x(t),0,sqrt(1/L(j)));
    splitEvt = A/(1+A)*stupdf(x(t),0,c/d,2*c);
    
    %Normalize
    w = [noChange,splitEvt];
    w = w/sum(w);
    
    %Update
    U = rand;
    if U<w(1)
        
        %No Change
        gibbs.z(t) = j;
        
    else
        
        %Split
        gibbs.z(t) = j;
        gibbs.z(t+1:end) = gibbs.z(t+1:end)+1;
        
        %Sample Precision for New State
        newLambda = sampleGammaPosterior(x(t),0,c,d);
        gibbs.lambda = [newLambda,gibbs.lambda];
        
    end
    
end

%Update n
gibbs.n = hist(gibbs.z,unique(gibbs.z));