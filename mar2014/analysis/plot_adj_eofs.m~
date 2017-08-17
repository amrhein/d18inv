clear
figure(1), clf, set(gcf,'color','w')
figure(2), clf, set(gcf,'color','w')
load /home/dan/Dropbox/2012-2013/Data/core_order_new

load /home/dan/Dropbox/2012-2013/MSfinal/4weddell/gfulls.mat

NS = 2806;
[r c] = size(gfulls);
nc = c/NS; % number of cores

%% all
% for ii = 1:nc
%
%     gg = gfulls(:,(1+NS*(ii-1)):NS*ii);
%     %[ug, sg, vg] = svd(full(gg),'econ');
%     [vy,iy] = max(sum(gg,2)); % time of max var
%     v = gg(iy,:)';
%     vg = v/norm(v);
%     ug = gg*vg; % PC-like quantity
%     pv = norm(ug)/norm(gg(:)); % percent var accounted for
%     ug = ug/norm(ug);
%
%     figure(1)
%     subplot(4,2,ii)
%     plot(ug(:,1))
% %    pv1 = sg(1,1).^2/sum(diag(sg).^2);
%     title(num2str(pv*100,3))
%
%     figure(2)
%     subplot(4,2,ii)
%     bubbsurf(vg,0.001,300)
%     title(core_order_new(ii))
% end

%% just ss cores
SC = 400;

figure(1), clf
figure(2), clf
iind = [6 8];
for jj = 1:2
    ii = iind(jj);
    
    gg = gfulls(:,(1+NS*(ii-1)):NS*ii);
    %[ug, sg, vg] = svd(full(gg),'econ');
    [vy,iy] = max(sum(gg,2)); % time of max var
    v = gg(iy,:)';
    vg = v/norm(v);
    ug = gg*vg; % project 'EOF' onto the field to obtain a PC-like quantity
    pv = norm(ug)/norm(gg(:)); % percent var accounted for
    ug = ug/norm(ug);
    
    figure(1)
    hold on
%    subplot(1,2,jj)
    if jj==1
        plot([0:200:5000],ug)
    else
        plot([0:200:5000],ug,'r')
    end
    xlim([0,5000])
%    title(num2str(pv*100,3))
    
    figure(2)
    hold on
    if jj == 1
        bubbsurf(-vg,0.001,SC)
    else
        bubbsurf_nomask(vg,0.001,SC)
    end
    %title(core_order_new(ii))
end


figure(2)
axis off
scatter(62,-90,SC*0.5,'r')
text(75,-90,core_order_new(8))
scatter(2,-90,SC*0.5,'b')
text(15,-90,core_order_new(6))
export_fig('-pdf','Figs/adj_EOFs')

figure(1)
xlabel('Time (yrs)')
ylabel('Normalized units')
set(gcf,'position',[377   642   560   187])
legend('NA3146', 'EP3210')

export_fig('-pdf','Figs/adj_PCs')

