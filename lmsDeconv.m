function [B,MSE] = lmsDeconv(f,d,M)
%LSQFilt - Demonstration routine for Least-Squares FIR filter design
%[B,MSE] = LSQFilt(f,d,M)
%   f - rowvector of data samples   - length N
%   d - row vector of desired values - length N
%   M - filter order
% Returns:
%   B - vector of optimal filter coefficients
%   MSE - minimized value of the mean-square-error
%
% Note: This routine is for tutorial purposes only. The Levinson method for
%   toeplitz matrix inversion would be used in practical methods.
%
% Based on:
% Author: D. Rowell
% Revised: 10/29/07
% https://ocw.mit.edu/courses/mechanical-engineering/2-161-signal-processing-continuous-and-discrete-fall-2008/lecture-notes/lecture_24.pdf
%-----------------------------------------------------------------------­
    pkg load signal;

N = length(f);
% Compute the correlation coefficients.
% Note that matlab defines the cross-correlaton backwards!! and
% we need to reverse the order of the subscripts.
%
    phiff=xcorr(f);
    phifd=xcorr(d,f);
%
% Extract out the central regions (low-lag values) and form
% the autocorrelation matrix.
%
    rff=phiff(N:N+M-1);
    P=phifd(N:N+M-1);
%
% Compute the optimal filter coefficients
%
    B = toepsolve(rff, P, M, 0);
%
% and the residual mean-square-error
%
    phidd=xcorr(d);
    MSE=phidd(N) - P'*B;
%
endfunction
