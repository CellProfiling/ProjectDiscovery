function make_remark_list(intfilepath,outpath)

if nargin<2 || isempty(outpath)
   outpath = '../results/high_remarks.csv'; 
end
    
addpath(genpath('~/Documents/jkitchin-matlab-cmu-1bc9774'))

load(intfilepath,'origIDremarksort','remarknumsort','remarkind','solUCstrcell','tot_hparesult','tot_gamerresult');

numtasks = length(origIDremarksort);
solUC_remarksort = solUCstrcell(remarkind);

[dictclasses,dictnames,dict_hash] = getDictionaries(0);
numclasses = length(dictnames);

high_remarks = remarknumsort>10;

remarknum_high = remarknumsort(high_remarks);
origID_high = origIDremarksort(high_remarks);
solUC_remark_high = solUC_remarksort(high_remarks);

num_high_remarks = sum(high_remarks);

fid = fopen(outpath,'w');
fprintf(fid,'Plate:well,FOV,Number of remarks,Annotated classes\n');

formatstr = '%s,%s,%d,%s\n';

annotated_classes = zeros(num_high_remarks,numclasses);

for i = 1:num_high_remarks
    curr_id_split = strsplit(origID_high{i},'_');
    
    curr_solUC = solUC_remark_high{i};
    
    for j = 1:length(curr_solUC)
        annotated_classes(i,:) = annotated_classes(i,:)+matchNum(str2num(curr_solUC{j}));
    end
    
    fprintf(fid,formatstr,strjoin(curr_id_split(1:2),':'),curr_id_split{3},...
        remarknum_high(i),strjoin(dictnames(logical(annotated_classes(i,:))),';'));
    
    
end

fclose(fid);

cmu_colors = @cmu.colors;
figure('Position', [1000 114 2000 700],'PaperPositionMode','auto');
remark_freq = sum(annotated_classes)./num_high_remarks;
hpa_freq = tot_hparesult./numtasks;
gamer_freq = tot_gamerresult./numtasks;

h = bar([remark_freq',gamer_freq',hpa_freq']);

h(1).FaceColor = cmu_colors('army green');
h(2).FaceColor = cmu_colors('celestial blue');
% h(3).FaceColor = cmu_colors('fern green');
%h(4).FaceColor = cmu_colors('celestial blue');
h(3).FaceColor = cmu_colors('gold (metallic)');

set(gca,'FontSize',28,'FontName','Arial','box','off','color','none',...
    'XTickLabel',dictnames,'XTick',1:length(dictnames),...
    'XTickLabelRotation',-45,'XColor',text_color,'YColor',text_color,'TickLength',[0,0],'YGrid','on');
legend({'Remarked','PD','HPA'})
export_fig(savename,'-transparent')
