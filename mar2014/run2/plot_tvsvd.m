% plot_tvsvd.m

clear
L = 5; % number of patterns to plot
load rs
load v230r
%% plot
close all
set(gcf,'color','w','position',[1 1 1024 1201])

% isolate the temporal filter part
% Vk = tf*tp (temporal filter times toeplitz)
%  tf = v230\toeplitz([u(:,1)',zeros(2806*100,1)],[u(1,1);zeros(100,1)])

% Assume that the temporal filter is nearly stationary

for ii = 1:L
    pii = (ii-1)*3 + 1; % subplot index
    subplot(L,3,pii:(pii+1))
    plotsurf(uv(:,ii));
    shading flat
    axis equal
    axis tight
    colorbar
    % compute the temporal filter
    % weights: w computed this way:
     w = v230r'*uv(:,ii);
    % is the same as vv(:,ii)*s(ii,ii).
    s1 = reshape(vv(:,ii),101,230);
    s1 = reshape(w,101,230);
    tf=s1*s1';
    % plot the temporal filter gain
    subplot(L,3,pii+2)
    [f,gain] = getgain(tf(50,:),1/200);
    plot(f,gain)
%    set(gca,'xticklabel',num2cell(1./get(gca,'xtick')))
    set(gca,'xticklabel',{'Inf','2000','1000','667','500',400'})
    ylabel('Gain')
    xlabel('Period (years)')
end

