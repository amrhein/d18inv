load solout230

subplot(3,1,1)
plot(diag(s),'k')
xlabel('Singular value index')
ylabel('Singular values')
axis tight

set(gcf,'color','w','position',[974 1355 358 419])

export_fig('-pdf','Figs/singvals230')
