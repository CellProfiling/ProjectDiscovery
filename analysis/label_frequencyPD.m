function label_frequencyPD(datapath_notuning,datapath_tuning,outpath,backgroundmode,version)
%This function creates the paired bar plot and statistical analysis of
%frequency of annotation for Project Discovery 
%
%Written By: Devin P. Sullivan 

addpath(genpath('~/Documents/jkitchin-matlab-cmu-1bc9774'))


%binknnpath = './binary_knn/binary_knn_predictions.mat';
% binknnpath = './binary_knn_predictions.mat';

if nargin<1 || isempty(datapath_notuning)
    %datapath = './pdAnalysisIntermediate/parsedTQdata201607121468334823219_tasks.mat';
    datapath_notuning = './pdDetailedIntermediate/parsedDetailedData_v3_notuning_mergelevel0.mat';
end
if nargin<2 || isempty(outpath)
    datapath_tuning = './pdDetailedIntermediate/parsedDetailedData_v3_mergelevel0.mat';
end
if nargin<3 || isempty(outpath)
    outpath = '../pd_results';
end
if nargin<4 || isempty(backgroundmode)
    backgroundmode = 'w';
end
if nargin<5 || isempty(version)
    warning('version not specified, defaulting to "tmp".');
    version = 'tmp';
end
if isnumeric(version)
    version = num2str(version);
end

%checks if the figure will be displayed on 'w' or black ('k') background
if strcmpi(backgroundmode,'w')
    text_color = 'k';
else
    text_color = 'w';
end

%load the necessary data
%load(datapath,'singlelocconfmat','locconfmat','locconfmatv14','dictnames','overannotation','underannotation','overannotationv14','underannotationv14');

load(datapath_tuning,'tot_gamerresult','tot_gamerresultv14')
tot_gamerresult_prop = tot_gamerresult;
tot_gamerresultv14_prop = tot_gamerresultv14;

load(datapath_notuning,'tot_hparesult','tot_gamerresult','tot_hparesultv14','tot_gamerresultv14','dictnames','score_Hammingv14');

num_tasks = length(score_Hammingv14);
num_tasksv14 = sum(~isnan(score_Hammingv14));
clear uniqueCodes 
%create percentages
%locconfmat_perc = locconfmat./sum(locconfmat(:));
%Sum out the EVE axis
% loc_percHPA = (sum(locconfmat,2)'+underannotation)./(repmat(sum(locconfmat(:)),1,size(locconfmat,1))+underannotation);
loc_percHPA = tot_hparesult./num_tasks;
% loc_percHPAv14 = (sum(locconfmatv14,2)'+underannotationv14)./(repmat(sum(locconfmatv14(:)),1,size(locconfmatv14,1))+underannotationv14);
loc_percHPAv14 = tot_hparesultv14./num_tasksv14;

%sum out the HPA axis
loc_percEVE = tot_gamerresult./num_tasks;
loc_percEVEv14 = tot_gamerresultv14/num_tasksv14;

%make the perc for proportional tuning 
loc_percEVE_prop = tot_gamerresult_prop./num_tasks;
loc_percEVEv14_prop = tot_gamerresultv14_prop./num_tasksv14;

%plot stuff
save_tot = [outpath,filesep,'freq_plotv',version,'_tune.png'];
save_v14 = [outpath,filesep,'freq_plotv',version,'v14_tune.png'];
% plot_freq(loc_percHPA(1:end-3),loc_percEVE(1:end-3),...
%     loc_percEVE_tune(1:end-3),perc_binknn(1:end-3),dictnames(1:end-3),save_tot,text_color)
% plot_freq(loc_percHPAv14(1:end-3),loc_percEVEv14(1:end-3),...
%     loc_percEVEv14_tune(1:end-3),dictnames(1:end-3),save_v14,text_color)
plot_freq(loc_percHPA(1:end-3),loc_percEVE(1:end-3),...
    loc_percEVE_prop(1:end-3),dictnames(1:end-3),save_tot,text_color)
plot_freq(loc_percHPAv14(1:end-3),loc_percEVEv14(1:end-3),...
    loc_percEVEv14_prop(1:end-3),dictnames(1:end-3),save_v14,text_color)



%compute statistics 


end

function plot_freq(freqHPA,freqEVE,freqEVE_prop,dictnames,savename,text_color)

figure('Position', [1000 114 2000 700],'PaperPositionMode','auto');
%h = bar([freqHPA',freqEVE',freqEVE_tune',freqEVE_prop']);
datamat = [freqHPA',freqEVE',freqEVE_prop'];
h = bar(datamat);

%h = bar([freqHPA',freqEVE']);
% ylim([0,max([freqHPA,freqEVE,freqEVE_prop])+0.1]);
ylim([0,round(max(datamat(:))+0.05,2)]);
% ylim([0,max([freqHPA,freqEVE,freqEVE_tune,freqEVE_prop])+0.1]);
set(gca,'ytick',0:0.1:1.0)
% h(1).FaceColor = [0,0,0.8];
% h(2).FaceColor = [0.8,0,0];
% h(3).FaceColor = [0,0.8,0];
% h(4).FaceColor = [0.8,0,0.8];
cmu_colors = @cmu.colors;
h(1).FaceColor = cmu_colors('gold (metallic)');
h(2).FaceColor = cmu_colors('cadmium red');
% h(3).FaceColor = cmu_colors('fern green');
%h(4).FaceColor = cmu_colors('celestial blue');
h(3).FaceColor = cmu_colors('celestial blue');

% h2 = legend('HPA','EVE','EVE_{proportional}')
%h2 = legend('HPA','EVE','EVE_{tuned}','EVE_{proportional}')
%h2 = legend('HPA','EVE')
% set(h2,'box','off','TextColor',text_color);
ylabel('Relative frequency')
%set(gca,'FontSize',28,'FontName','Ariel','YTick',[],'box','off','color','none','Xtick',[0:0.05:0.5]);
set(gca,'FontSize',28,'FontName','Arial','box','off','color','none',...
    'XTickLabel',dictnames,'XTick',1:length(dictnames),...
    'XTickLabelRotation',-45,'XColor',text_color,'YColor',text_color,'TickLength',[0,0],'YGrid','on');
export_fig(savename,'-transparent')
% saveas(gcf,[outpath,filesep,'freq_plot4.png'])

end

