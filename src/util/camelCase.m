function str = camelCase(str)
% Capitalizes only the first letter of each word
str=lower(str);
idx=regexp([' ' str],'(?<=\s+)\S','start')-1;
str(idx)=upper(str(idx));