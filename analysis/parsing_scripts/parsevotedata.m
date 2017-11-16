function [vc,vp,p] = parsevotedata(inputstr,dictclasses,merge_level,nplayers)
%This function takes a string from MMOS summary output and parses it
%
%INPUT:
%inputstr = a string in MMOS summary format for the column 'votes'
%dictclasses = a 1xn vector containing the MMOS codes for the various
%classes
%
%OUTPUT:
%vc = 1xn vector containing vote counts for each of the classes
%vp = 1xn vector containing vote percentages for each class
%p = 1xn probability of statistical significance for each class
%
%Written by: Devin Sullivan

if nargin<3
    merge_level = 0;
end

numclasses = length(dictclasses);

%allocate the outputs
vc = zeros(1,numclasses);
vp = zeros(1,numclasses);
p = ones(1,numclasses);


%split input string
%first replace beginning and end so they get split properly
inputstr2 = strrep(inputstr,'[{','');
inputstr2 = strrep(inputstr2,'}]','');
votedclasses = strsplit(inputstr2,'},{');

%go through each vote
for i = 1:length(votedclasses)
    r = 0;
    %split apart vote
    currvote = strrep(votedclasses{i},'"','');
    voteparts = strsplit(currvote,',');
    
    voteparts_split = cellfun(@(x) strsplit(x,':'),voteparts,'UniformOutput',false);
    
    %we need to update r first or we will not have a place to put stuff
    for j = 1:length(voteparts)
        
        %look for the different outputs
        if strcmpi(voteparts_split{j}{1},'r')
            %r output
            r = floor(str2num(voteparts_split{j}{2}).*10^(-merge_level));
            break
        end
    end
    
    currloc = find(dictclasses==r);
    
    for k = 1:length(voteparts)
        
        if strcmpi(voteparts_split{k}{1},'vc')
            vc(currloc) = vc(currloc)+str2num(voteparts_split{k}{2});
        elseif strcmpi(voteparts_split{k}{1},'vp')
            vp(currloc) = vp(currloc)+str2num(voteparts_split{k}{2});
        elseif strcmpi(voteparts_split{k}{1},'p')
            if strcmpi(voteparts_split{k}{2},'null')
                %p(currloc) = NaN;
                check = 1;
            else
                p(currloc) = min(p(currloc),str2num(voteparts_split{k}{2}));
            end
        elseif strcmpi(voteparts_split{k}{1},'r')
        else
            warning('unrecognized output. Make sure the script is up to date')
        end
    end
    
end


maxchoose = 5;
p_recalc = ones(size(p));
for i = 1:numclasses
    p_recalc(i) = calcOdds(vc(i),nplayers,numclasses,maxchoose);

end
p = p_recalc;