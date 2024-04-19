function h = plotStates(ax,stateVec,varargin)
% h = plotStates(ax,stateVec)
% takes the top [3xn] matrix of a stateVec and plots it on ax
[m,n] = size(stateVec);
if m < 3
    error('stateVec must have at least 3 rows')
end
% numPoints = 100;
% if n > 2
%     newStates = interparc(numPoints,stateVec(1,:),stateVec(2,:),stateVec(3,:));
%     newStates = newStates';
% else
    newStates = stateVec(1:3,:);
% end
% newStates = stateVec(1:3,:);
hh = plot3(ax,newStates(1,:),newStates(2,:),newStates(3,:),varargin{:});
if nargout == 1
    h = hh;
end
end