function compare_precision_recall(reannotation_path,pd_datapath_tuned,pd_datapath_notuning,dnn_pr_datapath,player_path)

fid = fopen(reannotation_path);

rawdata = textscan(fid,'%s%f','Whitespace','','TreatAsEmpty',' ','Delimiter',{':'});

labels = rawdata{1};
numbers = rawdata{2};

tot_loc = find(strcmpi(labels,'Total'));

tot_loc_labels = labels(tot_loc(1)+1:tot_loc(1)+8);
tot_loc_stats = numbers(tot_loc(1)+1:tot_loc(1)+8);

% expert_precision = tot_loc_stats(strcmpi(tot_loc_labels,'precision'));
% expert_recall = tot_loc_stats(strcmpi(tot_loc_labels,'recall'));
expert_precision = 0.74;
expert_recall = 0.6909524;


load(player_path,'badplayers','player_numtasks','avplayer_precision','avplayer_recall','avplayer_f1score')

%load the detailed data with tuning 
%load('./pdDetailedIntermediate/parsedDetailedData20170724_mergelevel0.mat','precision','recall');
load(pd_datapath_tuned,'precision','recall','f1_score');
consensus_precision = mean(precision(~isnan(recall)));
consensus_recall = mean(recall(~isnan(recall)));
clear precision 
clear recall

%load the detailed data without tuning
load(pd_datapath_notuning,'precision','recall');
notune_precision = mean(precision(~isnan(recall)));
notune_recall = mean(recall(~isnan(recall)));
clear precision 
clear recall

%Dnn results - taken from Accuracies_and_costs.ods
% dnn_pr_path = '../data/avg_p_r_LocCAT_allcells.txt';
fid = fopen(dnn_pr_datapath,'rb');
headers = textscan(fid, '%s',1, 'delimiter', '\n');
headers = strsplit(headers{1}{1},',');
num_columns = length(headers);

dnn_data = textscan(fid,'%s%q%q','Delimiter',',');
fclose(fid);

p_column = strcmpi(headers,'avg_p');
r_column = strcmpi(headers,'avg_r');

dnn_precision_list = str2double(dnn_data{p_column});
dnn_recall_list = str2double(dnn_data{r_column});
% remove_inds = dnn_precision==0&dnn_recall==0;
dnn_ind = strcmpi(dnn_data{1},'DNN');
dnn_precision = dnn_precision_list(dnn_ind);
dnn_recall = dnn_recall_list(dnn_ind);


TLv14_ind = strcmpi(dnn_data{1},'TLDNN_v14');
TLv14_pseudogamer_ind = strcmpi(dnn_data{1},'TLDNN_v14_FG');
TL_Precision = dnn_precision_list(TLv14_ind);%w/o pseudo gamer floored (Updated 10.10.2017)
TL_pseudoPrecision = dnn_precision_list(TLv14_pseudogamer_ind);%w/pseudo gamer
TL_Recall = dnn_recall_list(TLv14_ind);%w/o pseudo gamer floored (Updated 10.10.2017)
TL_pseudoRecall = dnn_recall_list(TLv14_pseudogamer_ind);%w/ pseudo gamer


contour_plotly2

wait = 1;

%hold on
%scatter(expert_precision,expert_recall,'MarkerSize',100)