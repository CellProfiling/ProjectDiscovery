function score_Players(inpath,if_imgs_path,merge_level,outputdir)
%if_images_path,cell_line,merge_level,parsed_outpath)
fstart = tic;
if nargin<1 || isempty(inpath)
    inpath = '../data/detailedData/results-detailed.tab';
end

if nargin<2 || isempty(if_imgs_path)
    if_imgs_path = '../data/hpa_results/IF_images_13062016.csv';
end

if nargin<3 || isempty(merge_level)
    merge_level = 0;
end

if nargin<4 || isempty(outputdir)
    outputdir = '../results/intermediate/';
end

%create output path 
outpath = [outputdir,filesep,'playerTrueAcc.mat'];
%check if it's already done 
if exist(outpath,'file')
    disp('Already computed score_Players, skipping  recompute');
    return
end



% if nargin<3 || isempty(if_images_path)
%     if_images_path = '../hpa_results/IF_images_13062016.csv';
% end

% if nargin<4 || isempty(cell_line)
%     cell_line = '';
% end

if ~isdir(outputdir)
    mkdir(outputdir)
end


[dictclasses,dictnames] = getDictionaries(merge_level);
nclasses = length(dictclasses);



%%%parameters
min_tasks = 10;


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

clear C

%double check that we got the right stuff
av_player_score = mean(playerScore);

% uniqueCodes = unique(originalCode);
[unique_pIDs,~,IDinds] = unique(playerID);
% pID_hash = java.util.HashMap;
% for i = 1:length(playerID)
%     if ~pID_hash.containsKey(playerID(i))
%         pID_hash.put(playerID(i),java.util.HashSet);
%     end
%     indset = pID_hash.get(playerID(i));
%     indset.add(i);
%    
% end

numPlayers = length(unique_pIDs);
% num_tasks = length(uniqueCodes);
% numvotes = cell(num_tasks,1);
% numplayers_tot = zeros(num_tasks,1);
% u_numplayers = zeros(num_tasks,1);
% tot_results = cell(num_tasks,1);
player_meanHamming  = nan(numPlayers,1);
avplayer_precision = nan(numPlayers,1);
avplayer_recall = nan(numPlayers,1);
avplayer_f1score = nan(numPlayers,1);
totplayer_precision = nan(numPlayers,1);
totplayer_recall = nan(numPlayers,1);
totplayer_f1score = nan(numPlayers,1);
naive_meanHamming  = nan(numPlayers,1);
player_minscore = nan(numPlayers,1);
player_maxscore = nan(numPlayers,1);
player_meanscore = nan(numPlayers,1);
% score_Hamming = cell(numPlayers,1);
player_numtasks  = nan(numPlayers,1);
playerID_mHamming_hash = java.util.HashMap;
badplayers = false(numPlayers,1);
invalid_submissions = false(length(IDinds),1);
%if ~exist(hpahashfile,'file')

[dictclasses,dictnames,dict_hash] = getDictionaries(merge_level);
nclasses = length(dictclasses);


[v14_hash,hpahash,hparaw] = loadHPAdata(if_imgs_path,dictnames,dict_hash);
clear dict_hash
%     save(hpahashfile,'v14_hash','hpahash','hparaw','dictclasses','dict_names')
%else
%    load(hpahashfile)
%end
toc(fstart)
loopstart = tic;
votes_naive = zeros(1,nclasses);
votes_naive(2) = 1;
for j = 1:numPlayers
    %     exp_cellline_hash = java.util.HashMap;
    %     vote_inds = strcmp(uniqueCodes{i},originalCode);
    currinds = IDinds==j;
    player_codes = originalCode(currinds);
    player_numtasks(j) = length(player_codes);
    curr_results = player_result(currinds);
    score_Hamming = nan(player_numtasks(j),1);
    naive_Hamming = nan(player_numtasks(j),1);
    all_hparesult = zeros(player_numtasks(j),nclasses);
    loc_truepos = zeros(1,nclasses);
    loc_trueneg = zeros(1,nclasses);
    loc_falsepos = zeros(1,nclasses);
    loc_falseneg = zeros(1,nclasses);
    for i = 1:player_numtasks(j)
        %         curr_playerIDs = playerID(vote_inds);
        %         curr_scores = playerScore(vote_inds);
        
        %track number of votes and number of players on each task
        %         numplayers_tot(i) = sum(vote_inds);
        %     numplayers_tot(i) = length(curr_playerIDs);
        %         u_numplayers(i) = length(unique(curr_playerIDs));
        %         tot_results{i} = player_result(vote_inds);
        %         [classes,numvotes{i}] = countVotes(tot_results{i});
        
        %         ptot(i) = calcOdds(numvotes{i}, numplayers_tot(i),nclasses,maxchoose)
        [~,~,votes_tot] = countVotes(curr_results{i},dictclasses,merge_level);
        
        currexpname = player_codes{i};%[if_imgs_data{2}{i},'_',if_imgs_data{3}{i},'_',if_imgs_data{4}{i},'_'];
        %     exp_cellline_hash.put(currexpname,ptot(i));
        hparesult = hpahash.get(currexpname);
        
        if isempty(hparesult)
            disp(['Annotation not found!',currexpname,'\n'])
            continue
        end
        
        if size(hparesult,1)~=size(votes_tot,1)
            hparesult = hparesult';
        end
        all_hparesult(i,:) = hparesult;
        [~,~,score_Hamming(i),~,~,...
            truepos,trueneg,falsepos,falseneg] = assess_task(hparesult,votes_tot);
        %calc truepos/falseneg etc
        %calc the true pos
        loc_truepos = loc_truepos+truepos;
        %calc the true neg
        loc_trueneg = loc_trueneg+trueneg;
        %calc the false positive
        loc_falsepos = double(loc_falsepos)+falsepos;
        %calc the false neg
        loc_falseneg = double(loc_falseneg)+falseneg;

%         score_Hamming(i) = sum(and(hparesult,votes_tot))/sum(or(hparesult,votes_tot));
        naive_Hamming(i) = sum(and(hparesult,votes_naive))/sum(or(hparesult,votes_naive));
%         loc_truepos = loc_truepos+and(hparesult,gamerresult);
%         loc_falsepos = loc_truepos+and(hparesult,gamerresult);
        %         [~,~,score_Hamming(i),~,~,...
        %             ~,~,~,~] = assess_task(hparesult,votes_tot);
        %if exphash.containsKey(well_name)
    end
    player_meanHamming(j) = mean(score_Hamming(~isnan(score_Hamming)));
    perclass_precision = loc_truepos./(loc_truepos+loc_falsepos);
    avplayer_precision(j) = mean(perclass_precision(~isnan(perclass_precision)));
    totplayer_precision(j) = sum(loc_truepos)./(sum(loc_truepos)+sum(loc_falsepos));
    perclass_recall = loc_truepos./(loc_truepos+loc_falseneg);
    avplayer_recall(j) = mean(perclass_recall(~isnan(perclass_recall)));
    totplayer_recall(j) = sum(loc_truepos)./(sum(loc_truepos)+sum(loc_falseneg));
    avplayer_f1score(j) = mean(2.*(avplayer_precision(j).*avplayer_recall(j))./(avplayer_precision(j)+avplayer_recall(j)));
    totplayer_f1score(j) = mean(2.*(totplayer_precision(j).*totplayer_recall(j))./(totplayer_precision(j)+totplayer_recall(j)));
    naive_meanHamming(j) = mean(naive_Hamming(~isnan(score_Hamming)));
    player_minscore(j) = min(playerScore(currinds));
    if player_minscore>0.6
        what= 1;
    end
    player_maxscore(j) = max(playerScore(currinds));
    player_meanscore(j) = mean(playerScore(currinds));

    if (naive_meanHamming(j)>player_meanHamming(j)) || (player_numtasks(j)<min_tasks)
    %if naive_meanHamming(j)>player_meanHamming(j)
%         if player_numtasks(j)>20
            badplayers(j) = true;
            invalid_submissions(currinds) = true;
%         end

    end
    playerID_mHamming_hash.put(unique_pIDs(j),player_meanHamming(j));
    
    %Uncomment this to save temporary results. I/O is expensive, but so is
    %re-running. 
    %     if mod(j,1000)==0
    %         %         toc(fstart)
    %         toc(loopstart)
    %         savestart = tic;
    %         save('playerTrueAcc.mat')
    %         toc(savestart)
    %         break
    %     end
    
end
looptime = toc(loopstart)
%     curr_code = uniqueCodes{i};
try
%     outpath = [outputdir,filesep,'playerTrueAcc.mat'];
    save(outpath)
    niceplot(player_numtasks,player_meanHamming,'num tasks (per-player)','mean Hamming score','','semilogx','.')
    player_performance(outpath);
catch
    error('problem saving or ploting player score. Check that plotly dependency is fulfilled for plotting.');
end


end


