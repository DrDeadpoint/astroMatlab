function [bool,errCode] = checkEquality(input1,input2)
%CHECKCLASSEQUALITY 
bool = true;
errCode = '';
if isa(input2,class(input1))
    switch class(input1)
        case {'double', 'char', 'logical'}
            bool = isequal(input1,input2);
            if ~bool
                [r,c] = size(input1);
                errInd = 0;
                if r == 1 || c==1 %vector or scalar
                    for i = 1:length(input1)
                        if ~isequal(input1(1), input2(1))
                            errInd = i;
                        end
                    end
                end
                errCode = ['"' class(input1) '" mismatch at index ' num2str(errInd)];
            end
        case 'cell'
            [r,c] = size(input1);
            for i = 1:r
                for j = 1:c
                    [bool,errCode] = checkEquality(input1{i,j},input2{i,j});
                    if ~bool
                        errCode = ['cell mismatch at index (' num2str(i) ',' num2str(j) '): ' errCode];
                        break
                    end
                end
            end
        case 'function_handle'
            funinfo1 = functions(input1);
            funinfo2 = functions(input2);
            bool = strcmp(funinfo1.function, funinfo2.function);
            if ~bool
                errCode = 'function_handle strings not the same';
            else
                [bool,errCode] = checkEquality(funinfo1.workspace{1}, funinfo2.workspace{1});
                if ~bool
                    errCode = ['function_handle mismatch: ' errCode];
                end
            end
        otherwise %structure or another class
            fnames = fieldnames(input1);
            fnames2 = fieldnames(input2);
            if ~isequal(fnames,fnames2)
                bool = false;
                errCode = 'Field mismatch';
            else
                for i = 1:length(fnames)
                    [bool,errCode] = checkEquality(input1.(fnames{i}), input2.(fnames{i}));
                    if ~bool
                        errCode = ['In field "' fnames{i} '": ' errCode];
                        break
                    end
                end
            end
    end
else
    bool = false;
    errCode = ['Class of input2 (' class(input2) ') does not match class of input1 (' class(input1) ')'];
end
end

