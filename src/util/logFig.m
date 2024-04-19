function logFig(fh,fname_noExt,description,parentFunc)
% logFig(fh,fname_noExt,description,parentFunc)
%
% save figures to a specified location and write a log about it
% to get parentFunc, use: P = mfilename('fullpath');
% if contains(fname_noExt,'.')
%    error('filename may not contain a period') 
% end
savepath = 'C:\Users\Alex\Code\Astrodynamics\figures\'; %laptop
set(fh,'PaperPositionMode','auto')
saveas(fh,[savepath fname_noExt '.fig']);
% exportgraphics(fh,[savepath fname_noExt '.png'],'ContentType','vector');
export_fig([savepath fname_noExt], '-png', '-transparent', '-r300', fh) %png are basically twice the size
% print(fh,[savepath fname_noExt],'-dpng','-r0',...
%                                 'InvertHardcopy','off',...
%                                 'PaperUnits','normalized',...
%                                 'PaperPosition',[0.05 0.05 0.95 0.95]); %grab current axes and scale by those values?
%print(fh,[savepath fname_noExt],'-dpdf','-bestfit',...
%                                'InvertHardcopy','off',...
%                                'PaperUnits','normalized');

if isfile([savepath 'figlog.mat'])
   load figlog.mat figCell
else
    figCell = cell(1,3);
    figCell{1,1} = 'Filenames';
    figCell{1,2} = 'Descriptions';
    figCell{1,3} = 'Parent Functions';
end
[r,~] = size(figCell);
figLine = r+1;
for i = 1:r
   if strcmp(figCell(i,1),fname_noExt)
       if strcmp(figCell(i,2),description)
           if strcmp(figCell(i,3),parentFunc)
               figLine = i;
               break
           end
       end
   end              
end
figCell{figLine,1} = fname_noExt;
figCell{figLine,2} = description;
figCell{figLine,3} = parentFunc;
save([savepath 'figlog.mat'],'figCell')
end