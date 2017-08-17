% [Y,f,amps,phis] = myfft(h,Fs)
% DEA 3 Mar 2011
% bugs fixed 10 Mar 2013
%
% streamlined FFT.
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
% Y                             % Fourier coefficients (complex)
% f                             % Freqs corresponding to the F coefs
% amps                          % amplitudes of Y (single-sided)
% phis                          % phases (single-sided)


function [Y,f,amps,phis] = myfft(h,Fs)

if nargin == 1
    Fs = 1;
end

L = length(h);                % Length of signal

L2f = floor(L/2);

Y = fft(h);
f = Fs/2*linspace(0,1,L2f+1);

% compute single-sided amplitude spectrum and phases
amps = 2*abs((Y(1:L2f+1)));

% there is no need to double the mean (zero frequency) value
amps(1) = 0.5*amps(1);

% if the record length is even, then there is no need to double the Nyquist frequency value
if ~mod(L,2)
  amps(end) = 0.5*amps(end);
end

%phis = phase(Y(1:L2f+1));
phis = angle(Y(1:L2f+1));

% Comment: March 2013
% nb that the one-sided amplitude spectrum is twice one side of the amplitudes.
% the one-sided power spectrum is twice one side of the power.
% thus the one-sided power spectrum is NOT the square of the one-sided amplitude spectrum.
% the one-sided power spectrum is 0.5 times the one-sided amplitude spectrum squared.
% The discrete Parseval relationship (see 6.39 in Carl's notes) is 
% N-1             N/2
% sum x_m^2 = 1/N*sum |xhat(s_n)|^2
% m=0             k=-N/2
% Thus in order to verify that Parseval's relation holds, try 
% 1/20*sum(abs(Y).^2)-sum(h.^2)
