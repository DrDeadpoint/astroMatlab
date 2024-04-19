function [cmap,lineCols] = colmap(ax,cmap,minVal,maxVal,allVals)
caxis(ax,[minVal maxVal])
colormap(ax,cmap)

ndim = ndims(allVals);
colRat = (allVals - minVal)./(maxVal - minVal);
colRat(colRat == 0) = 0.001;
colRat(colRat > 1) = 1;
lc1 = cmap(ceil(colRat*length(cmap)),1);
lc2 = cmap(ceil(colRat*length(cmap)),2);
lc3 = cmap(ceil(colRat*length(cmap)),3);
lc1 = reshape(lc1,size(allVals));
lc2 = reshape(lc2,size(allVals));
lc3 = reshape(lc3,size(allVals));
lineCols = cat(ndim+1,lc1,lc2,lc3);
end
