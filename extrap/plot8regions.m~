% Remember: different values of T come out of the ode solver!
% To do: 
% 1. find reasonable comparison to values equilibrating on
% their own. Is this useful?
% 2. Is the RMS metric a reasonable one? what do the fractions look like?

%load ../mar2014/setup/cadjs/MD99_2334Kadj.mat
%load ../mar2014/setup/cadjs/MD07_3076adj.mat
load ../mar2014/setup/cadjs/TR163_31badj.mat

% Difference in time to get Heavisides (need to look back at my
% notes to find why this is)
TSTEP = 200;
Ti = T(1):TSTEP:T(end);
addpath ../mar2014/utils
Csd = DGradient(C,T,1,'2ndOrder');
Csd(1,:) = 0; % require that inl times have 0 conc
Csd2 = DGradient(Csd,T,1,'2ndOrder');
Csd2(1,:) = 0;% require that inl times have 0 conc
Csd2i = interp1(T,Csd2,Ti);
nc = cumsum(Csd2);

% normalize by total at each time
Cn = bsxfun(@rdivide,nc,sum(nc,2)+.01);
Cnn = bsxfun(@rdivide,Cn,sum(nc,2)+.01);
% subtract off the final fraction:
Cns = bsxfun(@minus,Cn,Cn(end,:));
% RMS time series
D = sqrt(mean(Cns.^2,2));

figure(1); plot(D)

figure(2); plot(1-sum(nc,2))