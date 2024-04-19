% Old Constraint:
function F = con_inclination(allSTs,whichSeg,whichNode,whichPrimary,desInclination)
% i think this only works in the P1P2 rotating frame right now
ST = allSTs{whichSeg};
switch whichNode
    case 1 %t0 of segment
        seg_state = ST.STATE(1:6,1);
    case 2
        seg_state = ST.STATE(1:6,end);
    otherwise
        error('bad node')
end
% dot product between relative position and relative velociy vectors
% F = (nodePos - primaryPos) * (nodeVel - primaryVel)
mu = ST.SYSTEM_MODEL.muP1P2;
switch whichPrimary
    case 1
        primaryLocation = [-mu;0;0];
    case 2
        primaryLocation = [1-mu;0;0];
    otherwise
        error('Not a valid primary')
end
r = seg_state(1:3); v = seg_state(4:6);
r_rel = r - primaryLocation;
r_rel = r_rel./norm(r_rel); v = v./norm(v);
h = cross(r_rel,v); %angular momentum of orbit relative to primary
inc = acos(dot(h,[0;0;1])); %angle between ang mom and z-axis
F = inc - desInclination; %don't change this line unless you know what you're doing!!!
end
% function con = con_primaryApse(whichSeg,whichNode,whichPrimary)
% warning('Primary Apse currently only works when in the P1-P2 rotating frame')
% F = @(allSTs) conF(allSTs,whichSeg,whichNode,whichPrimary);
% DFrow = @(allSTs,x0key) partials(allSTs,x0key,whichSeg,whichNode,whichPrimary);
% conname = ['Apse Constraint: seg ' num2str(whichSeg)...
%     ', node ' num2str(whichNode) ', primary ' num2str(whichPrimary)];
% con = CONSTRAINT(conname,F,DFrow);
% end
% %% constraint evaluation
% function F = conF(allSTs,whichSeg,whichNode,whichPrimary)
%     [r,v,primaryLocation,~] = getrv(allSTs,whichSeg,whichNode,whichPrimary);
%     r_rel = r - primaryLocation;
%     % if i want to use this for for than the P1P2 rotating frame, need to
%     % account for motion of the primaries as well (location and velocity)
%     F = dot(r_rel,v); 
% end
% 
% %% partials
% function DFrow = partials(allSTs,x0key,whichSeg,whichNode,whichPrimary)
%     DFrow = nan(1,length(x0key));
%     [r,v,primaryLocation,ST] = getrv(allSTs,whichSeg,whichNode,whichPrimary);
%     r_rel = r - primaryLocation;
%     for i = 1:length(x0key)
%         [segString,freeVar] = strtok(x0key{i},' '); % of the form: 'seg# freeVar'
%         segNum = str2double(segString(4:end));
%         if segNum ~= whichSeg %not relevant
%             thisPartial = 0;
%         else 
%             desiredPartial = freeVar(2:end);
%             thisPartial = nan;
%             switch whichNode
%                 case 1 %with respect to initial node
%                     switch desiredPartial %string of which partial to take
%                         case 'x'
%                             thisPartial = v(1);
%                         case 'y'
%                             thisPartial = v(2);
%                         case 'z'
%                             thisPartial = v(3);
%                         case 'xd'
%                             thisPartial = r_rel(1);
%                         case 'yd'
%                             thisPartial = r(2);
%                         case 'zd'   
%                             thisPartial = r(3);
%                         case {'mass','t0','dt','SOCalpha0',...
%                                 'SOCalphadot','SOCbeta0','SOCbetadot'}
%                             thisPartial = 0;
%                         otherwise
%                             error('Partial does not exist')
%                     end
%                 case 2 %with respect to final node
%                     STM = ST.STM;
%                     switch desiredPartial
%                         case {'x','y','z','xd','yd','zd'}
%                             STMcol = find(strcmp(desiredPartial,{'x','y','z','xd','yd','zd'}));
%                             thisPartial = dot(STM(4:6,STMcol),r_rel) + dot(STM(1:3,STMcol),v);
%                         case 'mass'
%             
%                         case 't0'
%             
%                         case 'dt'
%             
%                         case 'uCoef'
%                         
%                         otherwise
%                             error('Partial does not exist')
%                     end
%             end
%         end
%         DFrow(i) = thisPartial;
%     end    
% end
% %% utility
% function [r,v,primaryLocation,segST] = getrv(allSTs,whichSeg,whichNode,whichPrimary)
%     segST = allSTs{whichSeg};
%     switch whichNode
%         case 1 %t0 of segment
%             seg_state = segST.STATE(1:6,1);
%         case 2
%             seg_state = segST.STATE(1:6,end);
%         otherwise
%             error('bad node')
%     end
%     % dot product between relative position and relative velociy vectors
%     % F = (nodePos - primaryPos) * (nodeVel - primaryVel)
%     mu = segST.SYSTEM_MODEL.muP1P2;
%     switch whichPrimary
%         case 1
%             primaryLocation = [-mu;0;0];
%         case 2
%             primaryLocation = [1-mu;0;0];
%         otherwise
%             error('Not a valid primary')
%     end
%     r = seg_state(1:3); v = seg_state(4:6);
% end

