function plotcomp( C, z, chanlocs, times )
nComp = size(C,2);
for iComp = 1:nComp
    subplot(2,nComp,iComp), topoplot(C(:,iComp),chanlocs);
end
subplot(2,nComp,nComp+1:2*nComp), plot(times,z');
end

