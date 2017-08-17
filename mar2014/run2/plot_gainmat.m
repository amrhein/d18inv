% plot_gainmat.m

addpath ../analysis
close all
set(gcf,'color','w','position',[1822 297 540 400])

load gainmat
SCALE = 500;
bubbsurf(-gainmat(1,:),10e-4,SCALE);

hold on

 lr = 1;
 scatter(2,-90,SCALE*1,'b')
 text(15,-90,[num2str(lr,1) ])
 scatter(71,-90,SCALE*lr/2,'b')
 text(82,-90,[num2str(lr/2,1)])
 scatter(135,-90,SCALE*0.1,'b')
 text(143,-90,[num2str(0.1,1) ])
 scatter(195,-90,SCALE*0.01,'b')
 text(202,-90,[num2str(0.01,1)])
 axis off
