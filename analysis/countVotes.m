function [u_classes,votes,votes_tot] = countVotes(votes,dictclasses,merge_level)

if ~iscell(votes)
    votes = {votes};
end

currnums = cell(1,length(votes));

for i = 1:length(votes)
    currvotes = votes{i};
    currnums{i} = cellfun(@str2num,strsplit(strrep(strrep(currvotes,'[',''),']',''),','));
end

allnums = cell2mat(currnums);
allnums = floor(allnums.*10^-merge_level);
[u_classes,ia,ic] = unique(allnums);

%force vector to ensure that each number is a bin. 
u_classes = [0,u_classes];
[votes] = hist(allnums,u_classes);
peripheryind = u_classes==3;
if merge_level==2 && sum(peripheryind)~=0
    cytoind = u_classes==2;
    votes(cytoind) = votes(cytoind)+votes(peripheryind);
    votes(peripheryind) = [];
    u_classes(peripheryind) = [];
end
if merge_level==2
    dictclasses(dictclasses==3) = [];
end



u_classes(1) = [];
votes(1) = [];
votes_tot = zeros(size(dictclasses));
for i = 1:length(u_classes)
    currind = (u_classes(i)==dictclasses);
    votes_tot(currind) = votes(i);
end