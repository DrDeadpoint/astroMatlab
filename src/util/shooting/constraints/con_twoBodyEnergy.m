function F = con_twoBodyEnergy(allSTs,whichSeg,whichNode,whichPrimary,des2Benergy)
ST_in = allSTs{whichSeg};
switch whichPrimary
    case 1
        mu = ST_in.SYSTEM_MODEL.P1DATA.MU;
        desFrame = 'P1centinert';
    case 2
        mu = ST_in.SYSTEM_MODEL.P2DATA.MU;
        desFrame = 'P2centinert';
    otherwise
        error('Not a valid primary')
end
% returns the difference between current two body energy and desired two body energy
ST_bodyCenteredInertial = CONVframe(ST_in,desFrame);
switch whichNode
    case 1 %t0 of segment
        seg_state = ST_bodyCenteredInertial.STATE(1:6,1);
    case 2
        seg_state = ST_bodyCenteredInertial.STATE(1:6,end);
    otherwise
        error('bad node')
end
r = seg_state(1:3); v = seg_state(4:6);
if strcmp('nd',ST_bodyCenteredInertial.STATE_UNIT)
    lstar = ST_bodyCenteredInertial.SYSTEM_MODEL.lstarP1P2;
    tstar = ST_bodyCenteredInertial.SYSTEM_MODEL.tstarP1P2;
    r = r*lstar;
    v = v*lstar/tstar;
end
F = twoBodyEnergy(r,v,mu)/des2Benergy - 1; %don't change this line unless you know what you're doing!!!
end