%%
load('../23july/4weddell/Gbpusv23July.mat','vs')
%%

eofs = [];
% iset = [1 105 210 315 420];
% liset = length(iset);
% 
% M = 856;
% 
% for ii = 1:liset
%     vst = vs(:,iset(ii));
%     xrt = reshape(vst,2806,[]);
%     [uxrt sxrt vxrt] = svd(full(xrt));
%     eofs(:,ii) = sparse(uxrt(:,1));
% end

% load Gmusv
liset = 9;
% eofs = vm(:,1:liset);

load ../23july/4weddell/gfulls23july.mat

Gm = reshape(sum(gfulls23July),2806,9);
eofs = Gm;
%%
M = 9;
eofw = nan(liset,M);
for ii = 1:M
    vst = vs(:,ii);
    xr = reshape(vst,2806,[]);
    %xrst = sum(xr,2);
    eofw(:,ii) = var(xr'*eofs)/sum(var(xr'));
    if any(eofw(:)>1), keyboard; end
end
    
save eofwgm eofw
    

plot(eofw')               
hold on
plot(diag(s)/max(diag(s)),'k')

