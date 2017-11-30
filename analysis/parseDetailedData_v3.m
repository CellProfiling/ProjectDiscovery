function parseDetailedData_v3(inpath,if_imgs_path, outpath, player_datafile, merge_level,tuning_cutoffs,cell_line)
%if_images_path,cell_line,merge_level,parsed_outpath)

%add all necessary files
addpath(genpath('./'))
%set the random seed
rng(3)

%%%Parameters 
score_cut = 0;
manual_cut = 0.01;
min_p_exp = 50;
z_099 = 2.58;%z_995=1.96;
p = 0.5;
me = 0.05;
min_vote_perc = 0.05;%minimum % of votes as decimal for ll_ratio test, 
%this avoids noise classification for very rare classes
% ll_dof = 1;%degrees of freedom for each per-class log-likelihood ratio test
% ll_alpha = 0.01;%significance level for ll_ratio test
maxchoose = 5;
all_PD_vote_path = './all_PD_votes.mat';

%%%

if nargin<1 || isempty(inpath)
    %inpath = '../detailedData/results-detailed.tab';
    inpath = '../detailedData/classifications.tab';
end

if nargin<2 || isempty(if_imgs_path)
    if_imgs_path = '../hpa_results/IF_images_13062016.csv';
end

if nargin <3 || isempty(outpath)
    outpath = ['./pdDetailedIntermediate',filesep];
end

if nargin<4 || isempty(player_datafile)
    player_datafile = './playerTrueAcc.mat';
end

if nargin<5 || isempty(merge_level)
    merge_level = 0;
end

if nargin<6 || isempty(tuning_cutoffs)
    tuning_cutoffs = true;
end


if nargin<7 || isempty(cell_line)
    cell_line = [];
end

% if nargin<3 || isempty(if_images_path)
%     if_images_path = '../hpa_results/IF_images_13062016.csv';
% end

% if nargin<4 || isempty(cell_line)
%     cell_line = '';
% end
whole_prog = tic;

if ~isdir(outpath)
    mkdir(outpath)
end


[dictclasses,dictnames,dict_hash] = getDictionaries(merge_level);
nclasses = length(dictnames);

% tic
[v14_hash,hpahash,hparaw,hpaheaders] = loadHPAdata(if_imgs_path,dictnames,dict_hash);
% toc

%grab the data
%tq_data = tdfread(inpath,'\t');
fid = fopen(inpath);
headerformat = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s';
headers = textscan(fid,headerformat,1);
%This is an outdated format Updated 16,03,16
%C = textscan(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s');
% headers = fgetl(fid);
% formatspec = [repmat('%q\t',1,length(headers)-1),'%q'];
%added because of bad_parsing
% headers = [headers(1:end-1),'updated_at2',headers(end)];
% formatspec = [formatspec,'\t%q'];
%formatspec = '%d\t%s\t%d\t%f\t%q';
% formatspec = '%d\t%s\t%d\t%f\t%q\t%s\t%s\t%{yyyy-mm-dd HH:MM:SS}D';
formatspec = '%d%s%d%f%q%s%s%s%s';

%C = textscan(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%q\t%s');
C = textscan(fid,formatspec,'Delimiter',{'\t',' '});
fclose(fid);

% taskID = C{1};
originalCode = C{2};
playerID = C{3};
playerScore = C{4};
playerResult = C{5};
remark = C{6};
is_remark = double(strcmp(remark,'{"remark":true}'));
clear remark
%circumstances = C{7};
createdAt_date = C{8};
createdAt_time = C{9};

clear C

%%Make per-cell line filtering
if ~isempty(cell_line)
    [exp_cellline_hash,~] = getCellLines(if_imgs_path);
    [is_cell_line] = get_selected_cellline(originalCode,exp_cellline_hash,cell_line);
    
    originalCode(~is_cell_line) = [];
    playerID(~is_cell_line) = [];
    playerScore(~is_cell_line) = [];
    playerResult(~is_cell_line) = [];
    is_remark(~is_cell_line) = [];
    %circumstances(~is_cell_line) = [];
    createdAt_date(~is_cell_line) = [];
    createdAt_time(~is_cell_line) = [];
end
%%

num_submissions = size(originalCode,1);

% load(player_datafile,'invalid_submissions','badplayers');
%
% originalCode(invalid_submissions) = [];
% playerID(invalid_submissions) = [];
% playerScore(invalid_submissions) = [];
% playerResult(invalid_submissions) = [];
% createdAt(invalid_submissions) = [];



%double check that we got the right stuff
av_player_score = mean(playerScore);

%get the number and list of unique tasks
[uniqueCodes,~,ic_code] = unique(originalCode);
num_tasks = length(uniqueCodes);

%calculate the number of samples needed for a good proportion for cutoff
%tuning
% z_099 = 2.58;%z_995=1.96;
% p = 0.5;
N =num_tasks;
% me = 0.05;
num_proportions = ((z_099^2*p*(1-p))/me^2)/(1+(z_099^2*p*(1-p))/(me^2*N));

%initialize stuff
numvotes = cell(num_tasks,1);
numplayers_tot = zeros(num_tasks,1);
u_numplayers = zeros(num_tasks,1);
tot_results = cell(num_tasks,1);
tot_remark_perc = zeros(num_tasks,1);

p_tot = nan(num_tasks,nclasses);
rr_tot = nan(num_tasks,nclasses);
h_tot = nan(num_tasks,nclasses);
% gamersol_tot = zeros(num_tasks,nclasses);
hpasol_tot = zeros(num_tasks,nclasses);
votes_tot = zeros(num_tasks,nclasses);


%initialize stuff
matchmat = zeros(1,nclasses);
tot_tasksv14 = 0;
nummatch_tot = 0;

%confusion matrices
locconfmat = zeros(nclasses,nclasses);
singlelocconfmat = zeros(nclasses,nclasses);
locconfmatv14 = zeros(nclasses,nclasses);
singlelocconfmatv14 = zeros(nclasses,nclasses);


%tot_resultcounts
tot_hparesult = zeros(1,nclasses);
tot_gamerresult = zeros(1,nclasses);
tot_hparesultv14 = zeros(1,nclasses);
tot_gamerresultv14 = zeros(1,nclasses);

%hamming scores
score_Hamming = nan(num_tasks,1);
score_Hammingv14 = nan(num_tasks,1);

%fp,fn,tp,tn
loc_truepos = zeros(1,nclasses);
loc_trueneg = zeros(1,nclasses);
loc_falsepos = zeros(1,nclasses);
loc_falseneg = zeros(1,nclasses);
loc_trueposv14 = zeros(1,nclasses);
loc_truenegv14 = zeros(1,nclasses);
loc_falseposv14 = zeros(1,nclasses);
loc_falsenegv14 = zeros(1,nclasses);

valid_players = playerScore>score_cut;


try
    load(['./pdDetailedIntermediate/parsedDetailedData_v7_mergelevel',num2str(merge_level),'_cell',cell_line,'.mat'],'p_tot','hpasol_tot');
catch
    %create a hashmap of original Codes.
    originalCode_hash = java.util.HashMap(num_tasks);
    disp('computing big map of sets')
    tic;
    for i = 1:length(originalCode)
        % for i = 1:1000
        if ~originalCode_hash.containsKey(originalCode{i})
            originalCode_hash.put(originalCode{i},java.util.HashSet);
        end
        originalCode_hash.get(originalCode{i}).add(i);
    end
    toc
    
    %get the total number of votes for each class in the project. This is
    %used for calculating the probability of a class being voted 
%     if ~exist(all_PD_vote_path,'file')
%         [voted_classes,tot_votes_per_class,votes_alldata] = countVotes(playerResult,dictclasses,merge_level);
%     else
%         load(all_PD_vote_path);
%     end
%     vote_prop_PD = votes_alldata./length(playerResult);
    hpa_ans = load(['./pdDetailedIntermediate/parsedDetailedData_v3_mergelevel',num2str(merge_level),'.mat'],'hpasol_tot');
    hpa_odds = sum(hpa_ans.hpasol_tot)./size(hpasol_tot,1);
    clear hpa_ans
    
    for i = 1:num_tasks
        %         for i = 1:100000
        curr_votes = originalCode_hash.get(uniqueCodes{i});
        if isempty(curr_votes)
            continue
        end
        currexp_name = uniqueCodes{i};
        vote_inds = cell2mat(cell(originalCode_hash.get(uniqueCodes{i}).toArray));%ic_code==i;
        vote_bool = false(num_submissions,1);
        vote_bool(vote_inds) = true;
        curr_playerIDs = playerID(vote_bool);
        %curr_scores = playerScore(vote_bool);
        
        vote_bool_trim = and(vote_bool,valid_players);
        
        %track number of votes and number of players on each task
        numplayers_tot(i) = sum(vote_bool_trim);
        %     numplayers_tot(i) = length(curr_playerIDs);
        u_numplayers(i) = length(unique(curr_playerIDs));
        tot_results{i} = playerResult(vote_bool_trim);
        
        %deal with total numer of 'remarks' - denotes strange sample
        tot_remark_perc(i) = sum(is_remark(vote_bool_trim))./numplayers_tot(i);
        clc
        
        [~,numvotes{i},votes_tot(i,:)] = countVotes(tot_results{i},dictclasses,merge_level);
        
        nchoose = sum(numvotes{i});
        %     nchoose = maxchoose*numplayers_tot(i);
        %%%Have to add a shift of 0.5 (or any number <1) to include the current
        %%%vote count!
%         p_tot(i,:) = hygecdf(votes_tot(i,:)-0.5,nclasses*numplayers_tot(i),numplayers_tot(i),nchoose,'upper');
        %%%Testing what happens when we use "hidden" or "blank" balls for
        %%%non-votes 20,10,2017
        p_tot(i,:) = hygecdf(votes_tot(i,:)-0.5,(nclasses+maxchoose-1)*numplayers_tot(i),numplayers_tot(i),maxchoose*numplayers_tot(i),'upper');
        gamer_odds = votes_tot(i,:)./numplayers_tot(i);
        rr_tot(i,:) = gamer_odds./hpa_odds;
        h_tot(i,:) = 2*(asin(gamer_odds)-asin(hpa_odds));
        
%         p_tot(i,:) = hygecdf(votes_tot(i,:)-0.5,nclasses*numplayers_tot(i),numplayers_tot(i),nchoose,'upper');
%         votes_tot(i,votes_tot(i,:)==0) = -realmin;
%         p_tot(i,:) = binocdf(votes_tot(i,:),nchoose,prob_class,'upper');
        %perform a log-ratio test
        %first compute likelihood 
%         null_likelihood = binopdf(votes_tot(i,:),numplayers_tot(i),vote_prop_PD);
%         alt_likelihood = binopdf(round(vote_prop_PD.*numplayers_tot(i)),numplayers_tot(i),vote_prop_PD);
%         [h,p_tot(i,:),~,~] = lratiotest(log(null_likelihood),log(alt_likelihood),ll_dof,ll_alpha);
%         [h,p_tot(i,:),~,~] = lratiotest(log(alt_likelihood),log(null_likelihood),ll_dof,ll_alpha);
        %Lastly, make sure we have a minimum of min_vote_perc (5%)
%         curr_vote_prop = (votes_tot(i,:)./numplayers_tot(i));
%         curr_vote_prop(curr_vote_prop<min_vote_perc) = 0;
        %double check that we only take over-represented things (one sided)
%         p_tot(i,:) = p_tot(i,:)+(curr_vote_prop<=vote_prop_PD)+realmin;
%         h_tot(i,:) = double( h.*(curr_vote_prop>vote_prop_PD));
        
        
        hparesult = hpahash.get(currexp_name)';
        
        if isempty(hparesult)
            disp(['Annotation not found!',uniqueCodes{i},'\n'])
            continue
        end
        hpasol_tot(i,:) = hparesult;
        
    end
end

%perform b-h correction for multiple hypothesis - doesn't actually make
%much of a difference 
% h_adj = nan(size(p_tot));
% adj_p = nan(size(p_tot));
% crit_p = nan(nclasses,1);
% adj_ci_cvrg = nan(nclasses,1);
% for i = 1:size(p_tot,2)
%     [h_adj(:,i), crit_p(i), adj_ci_cvrg(i), adj_p(:,i)]=fdr_bh(p_tot(:,i),manual_cut,'pdep','no');
% end


%compute the percentage of PD answers for each label for varying cutoffs of
%p.
PD_num = zeros(size(p_tot,2),min_p_exp);
for i = 1:size(p_tot,2)
    for j = 2:min_p_exp
        PD_num(i,j) = sum(p_tot(:,i)<10^-j);
    end
end
PD_perc = PD_num./repmat(num_tasks,size(PD_num,1),size(PD_num,2));

%grab a random sub-set of HPA annotated images based on the required
%sample-size calculation.
train_inds = randperm(num_tasks,ceil(num_proportions));
hpa_test = hpasol_tot(train_inds,:);
%estimate the HPA class distribution
hpa_perc = sum(hpa_test)./size(hpa_test,1);
%give rare classes a value so they don't get the maximium cutoff
% hpa_perc(hpa_perc==0) = mean(hpa_perc(hpa_perc~=0));

if tuning_cutoffs
    sig_cuts = zeros(1,size(hpa_perc,2));
    for i = 1:size(PD_perc,1)
        diff_prop = abs(PD_perc(i,:)-repmat(hpa_perc(i),1,size(PD_perc,2)));
        min_cut = find(diff_prop==min(diff_prop));
        if length(min_cut)>1
            min_cut = min_cut(1);
        end
        sig_cuts(i) = min_cut;%find(diff_prop==min(diff_prop));
        ideal_y(i,1) = PD_perc(i,sig_cuts(i));
        ideal_y(i,2) = hpa_perc(i);
    end
    sig_cuts = 10.^(-sig_cuts);
    sig_cuts(hpa_perc==0) = manual_cut;
else
    sig_cuts = repmat(manual_cut,1,size(hpa_perc,2));
end

gamersol_tot = p_tot<repmat(sig_cuts,num_tasks,1);
tot_gamerresult = sum(gamersol_tot);
tot_hparesult = sum(hpasol_tot);
tot_tasks = 0;

for i = 1:num_tasks
    %if the task was used for training cutoffs, skip it so as not to bias the results. 
    if any(i==train_inds)
        continue
    end
    %     ptot(i) = calcOdds(numvotes{i}, numplayers_tot(i),nclasses,maxchoose)
    %     gamersol_tot(i,:) = p_tot(i,:)<sig_cuts;%(10.^-sig_cuts);
    %     gamersol_tot(i,:) = p_tot(i,:)<sig_cuts;%(10.^-sig_cuts);
    gamerresult = double(gamersol_tot(i,:));
    
    %get the hpasolution from our hashmap
    
    %     hparesult = hpahash.get(currexp_name)';
    hparesult = hpasol_tot(i,:);
    %
        if isempty(hparesult)
            disp(['Annotation not found!',uniqueCodes{i},'\n'])
            continue
        end
    %     hpasol_tot(i,:) = hparesult;
    
    %track the total number of classifications of each category
    %     tot_hparesult = tot_hparesult+double(hparesult>0);
    %     tot_gamerresult = tot_gamerresult+gamerresult;
    %     toc
    %     tic
    [match_vec,nummatch,score_Hamming(i),singlelocconfpt,locconfpt,...
        truepos,trueneg,falsepos,falseneg] = assess_task(hparesult,gamerresult);
    %     toc
    %     tic
    matchmat = matchmat+match_vec;
    nummatch_tot = nummatch_tot+nummatch;
    singlelocconfmat = singlelocconfmat+singlelocconfpt;
    locconfmat  = locconfmat+locconfpt;
    
    %calc truepos/falseneg etc
    %calc the true pos
    loc_truepos = loc_truepos+truepos;
    %calc the true neg
    loc_trueneg = loc_trueneg+trueneg;
    %calc the false positive
    loc_falsepos = double(loc_falsepos)+falsepos;
    %calc the false neg
    loc_falseneg = double(loc_falseneg)+falseneg;
    
    
    
    
    %track those tasks that are only v14
    if v14_hash.contains(uniqueCodes{i})
        %tot class counts
        tot_hparesultv14 = tot_hparesultv14+hparesult;
        tot_gamerresultv14 = tot_gamerresultv14+gamerresult;
        %hamming score
        score_Hammingv14(i) = score_Hamming(i);
        %truepos etc
        loc_trueposv14 = loc_trueposv14+truepos;
        loc_truenegv14 = loc_truenegv14+trueneg;
        loc_falseposv14 = double(loc_falseposv14)+falsepos;
        loc_falsenegv14 = double(loc_falsenegv14)+falseneg;
        %heatmaps
        singlelocconfmatv14 = singlelocconfmatv14+singlelocconfpt;
        locconfmatv14 = locconfmatv14+locconfpt;
        tot_tasksv14 = tot_tasksv14+1;
    end
    tot_tasks = tot_tasks+1;
    %     toc
end

precision = loc_truepos./(loc_truepos+loc_falsepos);
recall = loc_truepos./(loc_truepos+loc_falseneg);
f1_score = 2.*(precision.*recall)./(precision+recall);
try
    if tuning_cutoffs
        outname = [outpath,'parsedDetailedData_v4_mergelevel',num2str(merge_level),'_cell',cell_line,'.mat'];
        hmp_outname = ['../pd_results_2017/heatmap_tune_predictions_cell',cell_line,'_v4.csv'];
        hmtp_outname = ['../pd_results_2017/heatmap_tune_truepos_cell',cell_line,'_v4.csv'];
        ml_outname = ['cell_',cell_line,'merge_level_v4_',num2str(merge_level),'.csv'];
        cutoffs_outname = ['../pd_results_2017/cutoffs_cell',cell_line,'_v4.csv'];

        fid = fopen(cutoffs_outname,'w');
        formatstr = '%s,%d\n';
        for i = 1:length(dictnames)
            fprintf(fid,formatstr,dictnames{i},sig_cuts(i));
        end
        fclose(fid);
        
        

    else
        outname = [outpath,'parsedDetailedData_v4_notuning_mergelevel',num2str(merge_level),'_cell',cell_line,'.mat'];
        hmp_outname = ['../pd_results_2017/heatmap_predictions_cell',cell_line,'_v4.csv'];
        hmtp_outname = ['../pd_results_2017/heatmap_truepos_cell',cell_line,'_v4.csv'];
        ml_outname = ['cell_',cell_line,'merge_level_notuning_v4_',num2str(merge_level),'.csv'];
    end
    save('-v7.3',outname,...
        'score_Hamming','loc_truepos','loc_trueneg','loc_falsepos','loc_falseneg',...
        'tot_hparesult','tot_gamerresult','tot_hparesultv14','tot_gamerresultv14',...
        'matchmat','singlelocconfmat','locconfmat','uniqueCodes','v14_hash',...
        'score_Hammingv14','loc_trueposv14','loc_truenegv14','loc_falseposv14','loc_falsenegv14',...
        'singlelocconfmatv14','locconfmatv14','precision','recall','f1_score','dictnames','hparaw',...
        'gamersol_tot','p_tot','hpasol_tot','sig_cuts','votes_tot','cell_line','rr_tot','h_tot');
    
    fid = fopen(hmp_outname,'w');
    formatstr = ['%s,',repmat('%d,',1,size(gamersol_tot,2)-1),'%d\n'];
    for i = 1:size(gamersol_tot,1)
        fprintf(fid,formatstr,uniqueCodes{i},gamersol_tot(i,:));
    end
    fclose(fid);
    
    fid = fopen(hmtp_outname,'w');
    formatstr = ['%s,',repmat('%d,',1,size(hpasol_tot,2)-1),'%d\n'];
    for i = 1:size(hpasol_tot,1)
        fprintf(fid,formatstr,uniqueCodes{i},hpasol_tot(i,:));
    end
    fclose(fid);
    
    fid = fopen(ml_outname,'w');
    fprintf(fid,'%s,%s,%s,%s,%s\n','Name','precision','recall','f1_score','count');
    for i = 1:length(dictnames)
        fprintf(fid,'%s,%f,%f,%f,%f\n',dictnames{i},precision(i),recall(i),f1_score(i),tot_hparesult(i));
    end    
    
    toc(whole_prog)
    
    
    
    
catch
    wait = 1;
end





end

function fdr = bhfdr(p)
% Adapted from the built-in mafdr package to work for matrices where
%columns are treated as independent tests of an experiment
% Compute the fdr using BH Linear step-up procedure
%
%Adapted by: Devin P Sullivan 30,10,2017
%Note: This is similar to the DATAMATRIX mode of the original but doesn't
%require a DATAMATRIX object


m = size(p,1);
[p_ord, idx] = sort(p);
fdr_ord =  p_ord .* repmat((m./(1:m))',1,size(p,2));
% Running min in reverse order (in-place)
%Add for loop for each column
fdr = nan(size(p));
for i = 1:size(p,2)
    bioinfoprivate.cumminmaxmex(fdr_ord(:,i),'min','reverse');
    fdr(idx(:,i),i) = fdr_ord;
end
end


