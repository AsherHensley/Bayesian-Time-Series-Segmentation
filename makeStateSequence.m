function [z,n] = makeStateSequence(alpha,nsamp)
%MAKESTATESEQUENCE Sample Yule-Simon State Sequence
%   [z,n] = makeStateSequence(alpha,nsamp) returns state sequence vector z
%   with vector of self transitions n given hyperparameter alpha and
%   sequence duration nsamp
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
z = zeros(1,nsamp);
z(1) = 1;
n = 1;

%Run Forward Process
for kk = 2:nsamp
    
    %Draw Sample
    nk = n(z(kk-1));
    p = nk/(nk+alpha);
    U = rand;
    
    %Update
    if U<=p
        z(kk) = z(kk-1);
        n(z(kk)) = n(z(kk))+1;
    else
        z(kk) = z(kk-1)+1;
        n = [n,1];
    end
    
end