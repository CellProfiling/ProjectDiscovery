function probHeatmapPD(pd_datapath, if_images_path,outputdir,backgroundmode,max_exp)
%This program
%%%
%%%

%%Info for eventual parser. Need to change function calls elsewhere first.
% (varargin)
%
% %set up parser
% p = inputParser;
%
% %set defaults
% default_datapath = './pdAnalysisIntermediate/parsedTQdata201611031478199274957_tasks.mat';
% default_if_images_path = '../hpa_results/IF_images_13062016.csv';
% default_outputdir = '../pd_results_2017/';
% default_backgroundmode = 'w';
% default_max_exp = 10;
%
% %set variable values
% addParameter(p,'pd_datapath',default_datapath);
% addParameter(p,'if_images_path',default_if_images_path);
% addParameter(p,'outputdir',default_outputdir);
% addParameter(p,'backgroundmode',default_backgroundmode);
% addParameter(p,'max_exp',default_max_exp);

addpath(genpath('/Users/devinsullivan/Documents/CellCycle/heatmaps/'))


if nargin<1 || isempty(pd_datapath)
    pd_datapath = './pdAnalysisIntermediate_2017/parsedTQdata20170307-tasks_mergelevel0.mat';
end

if nargin<2 || isempty(if_images_path)
    if_images_path = '../hpa_results/IF_images_13062016.csv';
end

if nargin<3 || isempty(outputdir)
    outputdir = '../pd_results_2017/';
end

if nargin<4 || isempty(backgroundmode)
    backgroundmode = 'w';
end
if nargin<5 || isempty(max_exp)
    max_exp = 10;
end

if strcmpi(backgroundmode,'w')
    text_color = 'k';
    nan_color = [1,1,1];
else
    text_color = 'w';
    nan_color = [0,0,0];
end

psig_cut = 0.01;
frac_min = 0.1;

%get colocalization probabilities from the HPA data
[coloc_nums,coloc_probs,class_probs,if_categories] = getColocProb(if_images_path);
[coloc_nums,coloc_probsv14,class_probs,if_categories] = getColocProb(if_images_path,true);

if_categories = cellfun(@(x) strrep(x,'Cytoskeleton (',''),if_categories,'UniformOutput',false);
if_categories = cellfun(@(x) strrep(x,')',''),if_categories,'UniformOutput',false);
if_categories{strcmpi('Cytoplasm',if_categories)} = 'Cytosol';
if_categories{strcmpi('Microtubule end',if_categories)} = 'Microtubule ends';
if_categories{strcmpi('Nucleoli (Fibrillar center',if_categories)} = 'Nucleoli fibrillar center';
if_categories{strcmpi('Focal Adhesions',if_categories)} = 'Focal adhesion sites';

%load counts of hpa-gamer data
% load(pd_datapath,'locconfmat','dictnames','tot_tasks');
load(pd_datapath,'locconfmat','dictnames',...
    'tot_tasks','loc_truepos','loc_falseneg','tot_gamerresult',...
    'tot_hparesult','tot_tasksv14','loc_trueposv14','loc_falsenegv14','tot_gamerresultv14',...
    'locconfmatv14','tot_hparesultv14');
totHPA_per_class = loc_truepos+loc_falseneg;
totHPA_per_classv14 = loc_trueposv14+loc_falsenegv14;
% tot_tasks = num_tasks;

%expected_mat = coloc_probs.*tot_tasks;
nclasses = size(locconfmat,1);
% dictnames_matched = dictnames;

% num_locs = length(if_categories);
% pval_mat = zeros(size(locconfmat));
expected_mat_ordered = zeros(size(locconfmat));
bfpval_heatmap = nan(size(locconfmat));
expected_mat_orderedv14 = zeros(size(locconfmat));
bfpval_heatmapv14 = nan(size(locconfmat));

colocPerc = zeros(size(locconfmat));
curr_coloc = zeros(size(locconfmat));
colocPercv14 = zeros(size(locconfmat));
curr_colocv14 = zeros(size(locconfmat));

%organize
for i = 1:nclasses
    curr_if_cat1 = strcmpi(dictnames{i},if_categories);
    if sum(curr_if_cat1)==0
        warning(['no if category found for ',dictnames{i}]);
        %         continue
    end
    
    for j = 1:nclasses
        
        curr_if_cat2 = strcmpi(dictnames{j},if_categories);
        %get the colocalization percentage
        colocPerc(i,j) = locconfmat(i,j)./totHPA_per_class(i);
        colocPercv14(i,j) = locconfmatv14(i,j)./totHPA_per_classv14(i);
        %check that this is an HPA-defined category, otherwise there is no
        %point in making the pvalue matrix for this.
        if sum(curr_if_cat2)==0 || sum(curr_if_cat1)==0
            continue
        end
        curr_coloc(i,j) = coloc_probs(curr_if_cat1,curr_if_cat2);
        curr_colocv14(i,j) = coloc_probsv14(curr_if_cat1,curr_if_cat2);
        expected_mat_ordered(i,j) = curr_coloc(i,j).*tot_hparesult(i);
        expected_mat_orderedv14(i,j) = curr_colocv14(i,j).*tot_hparesultv14(i);
        
        [if_categories(curr_if_cat1),if_categories(curr_if_cat2)]
        [expected_mat_ordered(i,j),locconfmat(i,j)]
        
        bfpval_heatmap(i,j) = min([1,(1-binocdf(locconfmat(i,j),tot_hparesult(i),curr_coloc(i,j)))*nclasses]);
        bfpval_heatmapv14(i,j) = min([1,(1-binocdf(locconfmatv14(i,j),tot_hparesultv14(i),curr_colocv14(i,j)))*nclasses]);
        
        
        
    end
    
    
end

%remove non-significant
bfpval_heatmap(bfpval_heatmap>psig_cut) = nan;
bfpval_heatmapv14(bfpval_heatmapv14>psig_cut) = nan;
% pval_heatmap(:,pval_heatmap>psig_cut) = nan;
%set non-hpa categories to nan;
bfpval_heatmap(totHPA_per_class==0,:) = nan;
bfpval_heatmap_log10 = log10(bfpval_heatmap);

outname = 'perc_hpagamer_heatmap_';
if ~isempty(strfind(pd_datapath,'notuning'))
    outname = [outname,'notuning_'];
end
bf_savename = [outputdir,filesep,outname,'v4.png'];
abs_bfpval_heatmap_log10 = abs(bfpval_heatmap_log10);
abs_bfpval_heatmap_log10(abs_bfpval_heatmap_log10>max_exp) = max_exp;
abs_bfpval_heatmap_log10(isnan(abs_bfpval_heatmap_log10)) = 0;
save_heatmap(abs_bfpval_heatmap_log10',dictnames,[bf_savename(1:end-3),'csv']);
plot_heatmap(bfpval_heatmap_log10,dictnames,nan_color,bf_savename,[1,11]);

bfpval_heatmapv14(totHPA_per_classv14==0,:) = nan;
bfpval_heatmapv14_log10 = log10(bfpval_heatmapv14);
abs_bfpval_heatmapv14_log10 = abs(bfpval_heatmapv14_log10);
abs_bfpval_heatmapv14_log10(abs_bfpval_heatmapv14_log10>max_exp) = max_exp;
abs_bfpval_heatmapv14_log10(isnan(abs_bfpval_heatmapv14_log10)) = 0;
outname = 'binomial_heatmap_';
if ~isempty(strfind(pd_datapath,'notuning'))
    outname = [outname,'notuning_'];
end
bfv14_savename = [outputdir,filesep,outname,'v14_v4.png'];
save_heatmap(abs_bfpval_heatmapv14_log10',dictnames,[bfv14_savename(1:end-3),'csv']);
plot_heatmap(bfpval_heatmapv14_log10,dictnames,nan_color,bfv14_savename,[1,11]);
% pval_heatmap(:,totHPA_per_class==0) = nan;


locconfper = locconfmat./repmat(tot_gamerresult,size(locconfmat,1),1);
locconfper(locconfper<frac_min) = NaN;
outname = 'perc_hpagamer_heatmap_';
if ~isempty(strfind(pd_datapath,'notuning'))
    outname = [outname,'notuning_'];
end
locconf_savename = [outputdir,filesep,outname,'v4.png'];
% abs_locconfper =

plot_heatmap(locconfper,dictnames,nan_color,locconf_savename,[0,1]);

locconfperv14 = locconfmatv14./repmat(tot_gamerresultv14,size(locconfmatv14,1),1);
locconfperv14(locconfperv14<frac_min) = NaN;
locconfv14_savename = [outputdir,filesep,outname,'v14_v4.png'];
plot_heatmap(locconfperv14,dictnames,nan_color,locconfv14_savename,[0.001,1.0])

end


function plot_heatmap(table_data,dictnames,nan_color,savename,min_maxvals)

if min(min_maxvals)<1
    fp_val = '%0.2f';
else
    fp_val = '%0.f';
end

%Display the pvalue heatmap
figure('Position', [100, 100, 1800, 1040],'PaperPositionMode','auto');
heatmap(abs(table_data), dictnames, dictnames, fp_val,'TickAngle',45,'ShowAllTicks',true, 'Colormap', @parula, ...
    'ColorBar', true,'ColorLevels',20,'NaNColor', nan_color,'UseLogColormap', false,'MinColorValue', ...
    min_maxvals(1),'MaxColorValue',min_maxvals(2),'FontSize',20,'GridLines','-');
set(gca,'FontSize',20,'FontName','Ariel','box','off','color','none',...
    'XColor',[0,0,0.8],'YColor',[0.8,0,0]);
export_fig(savename,'-transparent')
% saveas(gcf,[outputdir,filesep,'binomial_heatmap.tif']);
close all

end



function save_heatmap(table_data,dictnames,savename)

dispnames = dictnames;%translate_dictnames(dictnames);

numclasses = length(dictnames);
fid = fopen(savename,'w');
fprintf(fid,'X,');
for i = 1:numclasses-1;fprintf(fid,'%s,',dispnames{i});end
fprintf(fid,'%s\n',dispnames{end});
for i = 1:numclasses
    fprintf(fid,['%s,',repmat('%d,',1,numclasses-1),'%d\n'],dispnames{i},table_data(i,:));
end
fclose(fid);

end
