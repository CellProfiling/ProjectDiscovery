function masterPDanalysis_v2
%This is the wrapper script for the data analysis related to Project
%Discovery.
%
%Written By: Devin P Sullivan

%move to the analysis directory and perform analysis 
workingdir = which('masterPDanalysis_v2');
workingdir = strsplit(workingdir,'masterPDanalysis_v2');
cd(workingdir{1})
addpath(genpath('./'));


%Define paths
rawdatapath = '../data/20170307-tasks2.tab';
intresultspath = '../results/intermediate';
outputpath = '../results';
if ~isdir(outputpath)
    mkdir(outputpath)
end
if_images_path = '../data/hpa_results/IF_images_13062016.csv';
pd_data_path = '../data/detailedData/classifications.tab';
pd_results_path = '../data/detailedData/results-detailed.tab';
pd_per_round_path = '../data/mmos_data/20170307-tasks3.tab';
hpa_data_path = '../data/hpa_results/IF_images_13062016.csv';

%%Parsing parameters 
%merge_level - used for building the heirarchical trees with 0 being
%the leaf node
merge_level = 0;%no merging  
% merge_level = 1;%merge at 10^1: (101,102), (111,112) etc
% merge_level = 2;%merge at 10^2: (101-113), (201-214) etc
%cell_line - can be used to analyze a specific cell line
cell_line = [];%when empty, all cell lines will be analyzed 
%player_datafile - was used to prune gamers with low accuracy..
player_datafile = './playerTrueAcc.mat';%Depricated?
%tuning_cutoffs - boolean flag indicating either to use cutoff tuning or
%the standard p-value cutoff (manual_cut = 0.01); 
tuning_cutoffs = true;
%version - string controlling intermediate results name
tuning_version = 'v4';
notuning_version = 'v4_notuning';

%intermediate results path
intresults_file_tuning = [intresultspath,filesep,...
    'parsedDetailedData_',tuning_version,'_mergelevel',num2str(merge_level),'_cell',cell_line,'.mat'];
%     'parsedDetailedData_',tuning_version,'_mergelevel',num2str(merge_level),'_cell',cell_line,'.mat'];
intresults_file_notuning = [intresultspath,filesep,...
    'parsedDetailedData_',notuning_version,'_mergelevel',num2str(merge_level),'_cell',cell_line,'.mat'];
intresults_file_perday = [intresultspath,filesep,...
    'parsedTQdata20170307-tasks3_mergelevel0.mat'];

if ~exist(intresults_file_tuning,'file')
    parseDetailedData_v3(pd_data_path,hpa_data_path, intresultspath, player_datafile, merge_level,tuning_cutoffs,cell_line);
end

if ~exist(intresults_file_notuning,'file')
    tuning_cutoffs = false;
    parseDetailedData_v3(pd_data_path,hpa_data_path, intresultspath, player_datafile, merge_level,tuning_cutoffs,cell_line);
end

if ~exist(intresults_file_perday,'file')
    intresults_file_perday = parseTQdata_merge(pd_per_round_path,intresultspath,if_images_path,cell_line,merge_level);
end

%Generate per-player F1 plot and intermediate result
score_Players(pd_results_path,if_images_path,merge_level,intresultspath)

%Generate combined-round outputs
% generate_PDoutput_detailed(intresults_file_tuning,outputpath,merge_level)

%%Calculate precision and recall metrics
computePD_precisionRecall(intresults_file_tuning,outputpath)

%%Compute over/underrepresentation of labels
%paired bar graph
label_frequencyPD(intresults_file_notuning,intresults_file_tuning,outputpath)
%close all

%%Calculate multi-label heatmap
probHeatmapPD(intresults_file_tuning, if_images_path,outputpath)

%%Calculate per-day accuracy plot
per_day_stats(intresults_file_perday,outputpath)

%%Calculate player-rating graph

%%Create output for heirarchical tree plotting (done in python)

%%Generate precision-recall plot (Figure 6), final generation done in
%%python
reannotation_path = '../data/reannotation_results.txt';
pd_datapath_tuned = '../results/intermediate/parsedDetailedData_v4_mergelevel0_cell.mat';
pd_datapath_notuning = '../results/intermediate/parsedDetailedData_v4_notuning_mergelevel0_cell.mat';
dnn_pr_datapath = '../data/avg_p_r_LocCAT_allcells.txt';
player_path = '../results/intermediate/playerTrueAcc.mat';
compare_precision_recall(reannotation_path,pd_datapath_tuned,pd_datapath_notuning,dnn_pr_datapath,player_path)



%%%Calculate general citizen science metrics
%Citizen Science Metrics Paths
ccp_tutorialdata = '../data/ccp_data/tutorialstartfinish.xlsx';
mmos_stats = '../data/mmos_data/mmos_stats.txt';
pd_metricfile = '../data/PD_citizenscience_metrics.txt';
cox_file = '../data/Cox_et_al_Raw_Data_Summary.xlsx';
makeCitizenScienceMetrics(intresults_file_tuning,outputpath,pd_metricfile,ccp_tutorialdata,mmos_stats)

%compare PD to other citizen science projects
compareCitizenScienceEfforts(cox_file,[outputpath,filesep,pd_metricfile],outputpath)
