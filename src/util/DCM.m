function rot = DCM(angle,XYZ,isPassive)
% rot = DCM(angle,XYZ,isPassive)
%   this function produces a rotation matrix about a single axis
%       use om_conv to work with Euler Angles
% angle is in radians
% XYZ is 'x', 'y', or 'z'
% isPassive is a logical stating whether you want a Passive rotation (set it to false for row vector format)
%       (passive rotates the reference frame, not the vector)
% rot is the direction cosines matrix used to rotate a coordinate frame
% You should use Passive when creating a transformation matrix from euler angles (ie XYZ to xyz)

if strcmpi(XYZ,'x') || XYZ == 1 || strcmp(XYZ,'1')
    %rotation about x axis
    rot = [1, 0, 0; 0, cos(angle), sin(angle); 0, -sin(angle), cos(angle)];
elseif strcmpi(XYZ,'y') || XYZ == 2 || strcmp(XYZ,'2')
    %rotation about y axis
    rot = [cos(angle), 0, -sin(angle); 0, 1, 0; sin(angle), 0, cos(angle)];
elseif strcmpi(XYZ,'z') || XYZ == 3 || strcmp(XYZ,'3')
    %rotation about z axis
    rot = [cos(angle), sin(angle), 0; -sin(angle), cos(angle), 0; 0, 0, 1];
else
    error('no good rotation input')
end
if ~isPassive %rotate vector, active
    rot = rot'; %also row vector format
end

end