function F = con_declination(allSTs,whichSeg,whichNode,desDeclination)
ST = allSTs{whichSeg};
switch whichNode
    case 1 %t0 of segment
        seg_vel = ST.STATE(4:6,1);
    case 2
        seg_vel = ST.STATE(4:6,end);
    otherwise
        error('bad node')
end
% returns the difference between current inclination and desired
z = [0;0;1];
seg_vel = seg_vel./norm(seg_vel);
decl = acos(dot(z,seg_vel));
F = decl - desDeclination;
end