function con = con_timeofflight(segment,desTOF)
whichSeg = segment.seg_index;
F = @(segs,nodes) conF(segs,nodes,desTOF);
partials = @(segs,nodes,desPartials) calc_partials(segs,nodes,desPartials);
conname = ['Time of Flight: seg ' num2str(whichSeg)];
con = c_constraint(conname,F,partials,whichSeg,[]);
end
%% constraint evaluation
function F = conF(segs,nodes,desTOF)
    DT = getDT(segs,nodes);
    F = DT - desTOF; 
end

%% partials
function partials = calc_partials(segs,nodes,desPartials)
%     partials = nan(1,length(desPartials));
%     for i = 1:length(desPartials)
%         desiredPartialCell = desPartials(i,:);
%         desiredPartial = desiredPartialCell{3};
%         if ~strcmp(desiredPartialCell{1},'seg')
%             error('Incorrect partial desired for timeofflight')
%         end
%         thisSegInd = segs(desiredPartialCell{2}.SEG_INDEX);
%         thisSeg = nan;
%         for j = 1:length(segs)
%             if segs{j}.SEG_INDEX == thisSegInd
%                 thisSeg = segs{j};
%                 break
%             end
%         end
%         switch desiredPartial %string of which partial to take
%             case 'dt'
%                 thisPartial = 1*sign(thisSeg.DT);
%             case {'thrust','SOCalpha0',...
%                     'SOCalphadot','SOCbeta0','SOCbetadot'}
%                 thisPartial = 0;
%             otherwise
%                 error(['Partial "' desiredPartial...
%                     '" does not exist for Continuity constraint'])
%         end
%         partials(i) = thisPartial;
%     end    
end
%% utility
function DT = getDT(segs,nodes)
    DT = segs{1}.getDuration().value;
end