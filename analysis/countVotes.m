function [u_classes,vote_nums,votes_tot,player_vote] = countVotes(votes,dictclasses,merge_level)

if ~iscell(votes)
    votes = {votes};
end

numclasses = length(dictclasses);
numplayers = length(votes);

currnums = cell(1,numplayers);
player_vote = zeros(numplayers,numclasses);
for i = 1:length(votes)
    currvotes = votes{i};
    currnums{i} = cellfun(@str2num,strsplit(strrep(strrep(currvotes,'[',''),']',''),','));
    currnums{i} = floor(currnums{i}.*10^-merge_level);
    player_vote(i,:) = ismember(dictclasses,currnums{i});
end


allnums = cell2mat(currnums);
% allnums = floor(allnums.*10^-merge_level);
[u_classes,ia,ic] = unique(allnums);
votes_tot = sum(player_vote);


% %force vector to ensure that each number is a bin. 
u_classes = [0,u_classes];
[vote_nums] = hist(allnums,u_classes);
% peripheryind = u_classes==3;
% if merge_level==2 && sum(peripheryind)~=0
%     cytoind = u_classes==2;
%     vote_nums(cytoind) = vote_nums(cytoind)+vote_nums(peripheryind);
%     vote_nums(peripheryind) = [];
%     u_classes(peripheryind) = [];
% end
% if merge_level==2
%     dictclasses(dictclasses==3) = [];
% end

u_classes(1) = [];
vote_nums(1) = [];
% votes_tot = zeros(size(dictclasses));
% for i = 1:length(u_classes)
%     currind = (u_classes(i)==dictclasses);
%     votes_tot(currind) = vote_nums(i);
% end