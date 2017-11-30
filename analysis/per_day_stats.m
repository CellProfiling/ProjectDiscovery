function per_day_stats(parsedTQfile,outputdir)

addpath(genpath('./altmany-export_fig-f0af704'))

if nargin<2 || isempty(outputdir)
    outputdir = '../pd_results2017';
end

pthresh = 0.05;

%load the data
load(parsedTQfile,'runningdates','solUCnum',...
    'numclass_perfmatch','pmat','score_Hammingv14');

% [closevotes,perf_ind,entropy_onT,entropy_offT,mse] = analyze_Votes(vcmat,pmat,votecountnum,solUCnum,originalCode,'tmp_');

if isempty(runningdates)
    return
end

%convert HPA result to number
hparesults = zeros(size(pmat));
for i = 1:length(solUCnum)
    if all(isnan(pmat(i,:)))
        continue
    end
    hparesults(i,:) = matchNum(solUCnum{i});
end

%calculate truepos, falsepos, falseneg for F1 score later
pmat_binary = pmat<pthresh;
truepos = hparesults&(pmat_binary);
falsepos = (hparesults==0)&(pmat_binary);
falseneg = (hparesults==1)&(~(pmat_binary));

%get each day of the data
days = unique(runningdates);
numdays = max(days);
days = 1:max(days);
movavg_days = min(10,numdays/10);

percperf_perday = zeros(numdays,1);

f1_perday = zeros(numdays,1);
tot_correct = zeros(numdays,1);
% perf_perc = zeros(numdays,1);
% meanent_on = zeros(numdays,1);
% meanent_off = zeros(numdays,1);
numcurrday = zeros(numdays,1);
numcurrdayv14 = zeros(numdays,1);
totnum = 0;
%loop through each day and perform the analysis
for i = 1:numdays
    
    currdayinds = (runningdates==i);
    numcurrday(i) = sum(currdayinds);
    currday_solUCnum = solUCnum(currdayinds);
    
    
    %grab the info for those days
    numcorrect_day = numclass_perfmatch(currdayinds);
    totnum = totnum + length(numcorrect_day);
    
    
    %what the percentage per-day 
    num_correct = sum(numcorrect_day~=0);
    tot_correct(i) = sum(num_correct);
%     num_wrong = sum(numcorrect_day==0);
    percperf_perday(i) = num_correct./size(numcorrect_day,1);%(num_correct+num_wrong);
    
    %calculate f1 scores
    curr_tp = sum(truepos(currdayinds,:));
    curr_fp = sum(falsepos(currdayinds,:));
    curr_fn = sum(falseneg(currdayinds,:));
    curr_precision = sum(curr_tp)./(sum(curr_tp)+sum(curr_fp));
    curr_recall = sum(curr_tp)./(sum(curr_tp)+sum(curr_fn));
    f1_perday(i) = 2*curr_precision*curr_recall/(curr_precision+curr_recall);
    
    
    %curr_f1 = score_f1(currdayinds);
    %curr_f1 = curr_f1(~isnan(curr_f1));
    %score_f1_perday(i) = mean(curr_f1);
    %std_f1_perday(i) = std(curr_f1);
    %curr_f1v14 = score_f1v14(currdayinds);
    %curr_f1v14 = curr_f1v14(~isnan(curr_f1v14));
    %score_f1_perdayv14(i) = mean(curr_f1v14);
    %std_f1_perdayv14(i) = std(curr_f1v14);

    curr_hammingv14 = score_Hammingv14(currdayinds);
    curr_hammingv14 = curr_hammingv14(~isnan(curr_hammingv14));
    numcurrdayv14(i) = length(curr_hammingv14);
    
    %Get the current day info
%     pmat_day = pmat(currdayinds,:);
%     vcmat_day = vcmat(currdayinds,:);
%     vpmat_day = vpmat(currdayinds,:);
%     vc_day = votecountnum(currdayinds);
%     origID_day = originalCode(currdayinds);
   
       
%     [closevotes,perf_ind,ent_on,ent_off,mse(i)] =...
%         analyze_Votes(vcmat_day,pmat_day,vc_day,currday_solUCnum,origID_day,[pwd,filesep,'junedata',filesep,'day',num2str(i),filesep,'day',num2str(i)]);
%     
%     meanent_on(i) = mean(ent_on);
%     meanent_off(i) = mean(ent_off);
    
    %perf_perc(i) = sum(perf_ind)/size(perf_ind,1);
    
    
    
    
    
end

xlabelstr = 'Day';
ylabelstr = '% perfect match';
titlestr = '';
plottype = 'plot';
plotstyle = '.-';
newfigure = 1;
plotcolor = [0.3,0.3,0.3];
H = niceplot(days,percperf_perday,xlabelstr,ylabelstr,titlestr,...
    plottype,plotstyle,newfigure,plotcolor)


% set(gcf, 'Color', 'k')
% set(gca, 'Color', 'k')
% set(gca, 'XColor','w')
% set(gca, 'YColor','w')
% cb = colorbar;
% cb.Color = [1,1,1];
% set(gcf, 'InvertHardCopy', 'off');
H.PaperPositionMode = 'auto';
% export_fig([outputdir,filesep,'full_overlap_perday.png'],'-transparent')
%saveas(H,[outputdir,filesep,'full_overlap_perday',date,'.tif'])

%%f1 per-day
%output f1 score per-day
savename_f1 = [outputdir,filesep,'f1_score_perday.png'];
plot_f1perday(days,f1_perday,movavg_days,savename_f1);

%%Images analyzed per-day
%output images per-day 
savename_imgs = [outputdir,filesep,'participation_perday.png'];
plot_activity(days,numcurrday,movavg_days,savename_imgs)

%output images per-day 
savename_imgsv14 = [outputdir,filesep,'participation_perdayv14.png'];
plot_activity(days,numcurrdayv14,movavg_days,savename_imgsv14)
%%

end



function plot_f1perday(days,f1_daily,movavg_days,savename_f1)

close all

%remove any nans
f1_daily(isnan(f1_daily)) = 0;

%calulate moving average 
mov_av_f1 = tsmovavg(f1_daily,'s',movavg_days,1);

% score_f1_perday(isnan(score_f1_perday)) = 0;
% std_f1_perday(isnan(std_f1_perday)) = 0;
% score_f1_perdayv14(isnan(score_f1_perdayv14)) = 0;
% std_f1_perdayv14(isnan(std_f1_perdayv14)) = 0;


%sem = std_daily./sqrt(numcurrday);
%sem(isnan(sem)) = 0;

ylabelstr = 'f1 score';
xlabelstr = 'Day';
H = niceplot(days,f1_daily,xlabelstr,ylabelstr,'','plot','-',1,[0.5,0.5,0.5])
set(gca,'FontSize',28,'FontName','Arial','box','off','color','none');%,...
    %'XColor',[0.3,0.3,0.3],'YColor',[1,0.5,0.5]);

hold on
plot(mov_av_f1,'k','LineWidth',2)

export_fig(savename_f1,'-transparent')
hold off

end

function plot_activity(days,numcurrday,movavg_days,savename_imgs)

close all

%calculate the moving average 
mov_av_numcurrday = tsmovavg(numcurrday,'s',movavg_days,1);

%yyaxis right
xlabelstr = 'Day';
ylabelstr = 'Images analyzed';
niceplot(days,numcurrday,xlabelstr,ylabelstr,'','plot','-',1,[0.5,0.5,0.5])
hold on
plot(days,mov_av_numcurrday,'k','LineWidth',2)

set(gca,'FontSize',28,'FontName','Arial','box','off','color','none');%,...
    %'XColor',[0.3,0.3,0.3],'YColor',[0.5,1,0.5]);
%export_fig([outputdir,filesep,'f1_score_perday.png'],'-transparent')
export_fig(savename_imgs,'-transparent')
end


