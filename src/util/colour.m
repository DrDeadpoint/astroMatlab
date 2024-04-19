function rgb = colour(col_string)
% rgb = colour(col_string)
% gets a color based on input string
% if numeric: rgbymconpwbg is the order
% possible colours:
%   'k' black
%   'b' blue
%   'r' red
%   'g' green
%   'y' yellow
%   'm' magenta
%   'w' white
%   'c' cyan
%   'gray'
%   'o' orange
%   'n' navy
%   'plum'
%   'p' purple
%   'teal'
%   'tan'

colList = 'wrgbymcotp';
if isnumeric(col_string)
   colNum = mod(col_string,length(colList))+1;
   col_string = colList(colNum);
end
switch col_string
    case {'k','black'}
        rgb = [0 0 0];
    case {'b','blue'}
        rgb = [0 0 153];
    case {'r','red'}
        rgb = [153 0 0];
    case {'g', 'green'}
        rgb = [0 125 0];
    case {'y','yellow'}
        rgb = [179 151 8];
        rgb = [219 180 8];
    case {'m','magenta'}
        rgb = [199 21 133];
    case {'w','white'}
        rgb = [255 255 255];
    case {'c','cyan'}
        rgb = [0 220 220];
    case {'gray'}
        rgb = [96 96 96];
    case {'o','orange'}
        rgb = [255 128 0];
    case {'n','navy'}
        rgb = [0 0 102];
    case {'plum'}
        rgb = [102 0 51];
    case {'p','purple'}
        rgb = [102 0 204];
    case {'teal'}
        rgb = [24 198 122];
    case {'t','tan'}
        rgb = [245 222 179];
    case {'test'}
        ffigure;
        hold on
        for i = 1:length(colList)
            plot([i-1 i],[0 1],'-','Color',colour(i),'LineWidth',3);
            colNum = mod(i,length(colList))+1;
            col_string = colList(colNum);
            text(i-0.5,0.5,col_string);
        end
        rgb = 0;
    otherwise
        error([col_string ' is an unknown color string'])
end
rgb = rgb/255;
end