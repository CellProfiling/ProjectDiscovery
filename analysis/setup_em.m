function setup_em(inpath,intresultspath)
%Expected inputs: inpath,intresultspath

if nargin<1
inpath = '../data/detailedData/results-detailed.tab';
end
if nargin<2
intresultspath = '../results_archive/intermediate/parsedDetailedData_v3_mergelevel0.mat';
end

%parameters 
merge_level = 0;
maj_vote = 0.5;


load(intresultspath,'hpasol_tot');
w_0 = sum(hpasol_tot)./size(hpasol_tot,1);

%grab the data
%tq_data = tdfread(inpath,'\t');
fid = fopen(inpath);
%This is an outdated format Updated 16,03,16
%C = textscan(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s');
% headers = fgetl(fid);
% formatspec = [repmat('%q\t',1,length(headers)-1),'%q'];
%added because of bad_parsing
% headers = [headers(1:end-1),'updated_at2',headers(end)];
% formatspec = [formatspec,'\t%q'];
formatspec = '%d\t%s\t%d\t%f\t%q';

%C = textscan(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%q\t%s');
C = textscan(fid,formatspec);
fclose(fid);

% taskID = C{1};
originalCode = C{2};
playerID = C{3};
playerScore = C{4};
player_result = C{5};

%setup dictionaries
[dictclasses,dictnames,dict_hash] = getDictionaries(merge_level);
num_classes = length(dictclasses);
%load HPA data - this shouldn't actually be needed for the base algorithm
% [v14_hash,hpahash,hparaw] = loadHPAdata(if_imgs_path,dictnames,dict_hash);

[uplayerID,pia,pic] = unique(playerID);
num_players = length(uplayerID);
[uniqueCodes,ia,ic] = unique(originalCode);

% for z = 1:length(originalCode)
%     
% end

num_tasks = length(uniqueCodes);

cat_votes_tot = zeros(1,num_classes);
num_votes_pertask = zeros(num_tasks,1);
players_pertask = zeros(num_tasks,1);
task_votes = cell(num_tasks,1);
prop_votes = zeros(num_tasks,num_classes);
seen_task = sparse(num_tasks,num_players);
task_player_nmj = cell(num_classes,1);
for i = 1:num_classes
    task_player_nmj{i} = sparse(num_tasks,num_players);
end
Tinit = sparse(num_tasks,num_classes);
pihat_k_num = zeros(num_classes,num_classes,num_players);

%find out which players have seen which tasks 
for m = 1:num_tasks
    tic
%     code = uniqueCodes(i);
    currinds = ic==m;
    curr_results = player_result(currinds);
    curr_players = pic(currinds);
    seen_task(m,curr_players) = 1;
    num_votes_pertask(m) = length(curr_results);
    players_pertask(m) = length(unique(curr_players));
    %convert votes into actual counts
    [~,~,votes_sample,task_votes] = countVotes(curr_results,dictclasses,merge_level);
    prop_votes(m,:) = votes_sample./num_votes_pertask(m);
    Tinit(m,:) = prop_votes(m,:)>maj_vote;
    %assign votes to players
    for k = 1:num_votes_pertask(m)
%         pihat_k_num = 0;
        for i = 1:num_classes
            %if it's currently a solution, add to confusion matrix
            if Tinit(m,i)
                %%%Something needs to be fixed here
                pihat_k_num(i,:,curr_players(k)) = pihat_k_num(i,:,curr_players(k))+task_votes(k,:);
                %%%
            end
            
            %if it doesn't have a vote, skip it
            if task_votes(k,i)==0
                continue
            end
            %if it does, put it in the correct cellmat
            task_player_nmj{i}(m,curr_players(k)) = task_votes(k,i);%+vote_task(curr_players(k),:);
        end
    end

    cat_votes_tot = cat_votes_tot+votes_sample;
   toc
   holdup=m
    
end

try
save('setup_em_intermediate.mat')
catch
    whoops
end

pi_hat = pihat_k_num./repmat(sum(pihat_k_num,1),num_classes,1,1);
p_hati = zeros(1,num_classes);
for i = 1:num_classes
    p_hati(i) = nnz(Tinit(:,i))./num_tasks;
end

middle_prod = zeros(num_tasks,num_classes);
for m = 1:num_tasks
    inner_prod = ones(num_players,num_classes);
    for i = 1:num_classes
        for j = 1:num_classes
            inner_prod(:,i) = inner_prod(:,i).*squeeze(pi_hat(i,j,:)).^(task_player_nmj{j}(m,:)');
        end
    end
    middle_prod(m) = prod(inner_prod);
end
outer_prod = prod(p_hati.*middle_prod);

min_votes = min(num_votes_pertask);
% num_votes_tot = sum(num_votes);
tot_votes_used = min_votes*num_tasks;
mu_0 = cat_votes_tot./sum(num_votes_pertask);
%task_votes_trim = cellfun(@(x) x(1:min_votes,:),task_votes,'UniformOutput',false);
currcol = 26;
task_votes_trim = cellfun(@(x) x(1:min_votes,currcol),task_votes,'UniformOutput',false);
pd_em(task_votes_trim,mu_0(currcol),w_0(currcol),prop_votes(:,currcol)>0.5)

