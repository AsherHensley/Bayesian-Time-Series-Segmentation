function u = sampleGammaPosterior(x,mu,a0,b0)
%SAMPLEGAMMAPOSTERIOR Draw Sample from Gamma Posterior
%   u = sampleGammaPosterior(x,mu,a0,b0) returns sample from gamma
%   posterior u given data x, mean mu, and hyperparameters a0 and b0. 
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

%Compute Posterior Parameters
n = length(x);
a = a0+n/2;
b = b0+0.5*sum((x-mu).^2);

%Draw Sample
u = gamrnd(a,1/b);