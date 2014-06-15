function plotcomp( C, s, chanlocs, times )
nComp = size(C,2);
s = detrend(s')';
s = s - repmat(mean(s(:,times<0),2),[1,size(s,2)]);
s = s ./ repmat(std(s(:,times<0),0,2),[1,size(s,2)]);
for iComp = 1:nComp
    subplot(2,nComp,iComp), topoplot(C(:,iComp),chanlocs);
end
subplot(2,nComp,nComp+1:2*nComp), plot(times,s');
end