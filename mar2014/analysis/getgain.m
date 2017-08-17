% [f,gain] = getgain(h,Fs)
% DEA 3 Mar 2011
% bugs fixed 10 Mar 2013
%
% Returns filter gain given a kernel - it's just the one-sided amplitudes of the
% Fourier transform
%
% % INPUT VARIABLES 
% h                             % time series
% Fs                            % Sampling freq, in units^{-1}
%
% % INTERNAL VARIABLES
% T = 1/Fs;                     % Sample time
% L = length(h);                % Length of signal
% t = (0:L-1)*T;                % Time vector
%
% % OUTPUT VARIABLES
% f                             % Freqs corresponding to the F coefs
% gain                          % amplitudes of Y (single-sided)

function [f,gain] = getgain(h,Fs)

if nargin == 1
    Fs = 1;
end

L = length(h);                % Length of signal

L2f = floor(L/2);

Y = fft(h);
f = Fs/2*linspace(0,1,L2f+1);

% compute single-sided amplitude spectrum and phases
gain = abs((Y(1:L2f+1)));
%phis = angle(Y(1:L2f+1));
