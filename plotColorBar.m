figure,
cbar_handle = colorbar('SouthOutside');
caxis([0 1]);
colormap('gray');
set(gcf,'Color','None');
export_fig('colorbar.pdf', cbar_handle);