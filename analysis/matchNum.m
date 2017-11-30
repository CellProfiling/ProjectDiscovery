function [result] = matchNum(code,merge_level) 

if nargin<2
    merge_level = 0;
end

dictclasses = [101,102,111,112,113,121,122,123,131,201,202,203,204,211,212,213,214,215,221,222,231,232,233,301,302,303,801,901,902];


if merge_level==1
    dictclasses = floor([101,102,111,112,113,121,122,123,131,201,202,203,204,211,212,213,214,215,221,222,231,232,233,301,302,303,801,901,902]./10);
elseif merge_level==2
    dictclasses = floor([101,102,111,112,113,121,122,123,131,201,202,203,204,211,212,213,214,215,221,222,231,232,233,301,302,303,801,901,902]./100);
end
u_classes = unique(dictclasses);


if iscell(code)
    number = zeros(length(code),1);
    for i = 1:length(code)
        if ischar(code{i})|| isnan(code{i}) || code{i}==0
            if code{i}==0
                wait = 1;
            end
            continue
        end
        if merge_level==1
            code{i} = floor(code{i}./10);
        elseif merge_level==2
            code{i} = floor(code{i}./100);
        end
        %[~,number(i)] = find(dictclasses==code{i});
        [~,number(i)] = find(u_classes==code{i});
    end
elseif code==0
    number = 0;
    %     found = 0;

else
    number = zeros(length(code),1);
    for i = 1:length(code)
        if ~isnan(code(i))
            if merge_level==1
                code{i} = floor(code(i)./10);
            elseif merge_level==2
                code{i} = floor(code(i)./100);
            end
        %[~,number(i)] = find(dictclasses==code(i));
        [~,number(i)] = find(u_classes==code(i));
        else
            number(i) = nan;
        end
    end
    number(isnan(number)) = [];
end

%result = false(1,length(dictclasses));
result = false(1,length(u_classes));
if number ~= 0
result(number) = 1;
end