

for ii = 1:8
    load(['ttdsummary' num2str(ii)])
    Call(:,:,ii) = C;
end

%%
ind = 71000;
c8 = squeeze(Call(:,ind,2:8));
c8n = bsxfun(@rdivide,c8,Call(:,ind,1)+.1);
c8nn = bsxfun(@rdivide,bsxfun(@rdivide,c8,Call(:,ind,1)),Call(:,ind,1));
figure(1)
clf
subplot(2,1,1)
plot(T,c8)
subplot(2,1,2)
plot(T,c8n)
%figure(2)
%plot(T,c8nn)

%%

[u s v] = svd(c8);

plot(T,u(:,1)*s(1,1)*v(:,1)')
plot(T,u(:,2)*s(2,2)*v(:,2)')

%%

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