function alpha = sampleAlphaGibbs(alpha0,n,a,b)
%SAMPLEALPHAGIBBS Run Gibbs Sampling Inference Procedure for alpha
%   alpha = sampleAlphaGibbs(alpha0,n,a,b) returns sample for alpha
%   posterior
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

ngibbs = 1000;
N = length(n);
alpha = zeros(1,N);
alpha(1) = alpha0;
for kk = 2:ngibbs
    w = -log(betarnd(alpha(kk-1)+1,n));
    alpha(kk) = gamrnd(a+N,1/(b+sum(w)));
end
alpha = alpha(end);