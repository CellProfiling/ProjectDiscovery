function [precision_perclass,recall_perclass] = computePD_precisionRecall(datapath,outpath)

if nargin<1 || isempty(datapath)
    datapath = './pdAnalysisIntermediate/parsedTQdata201607121468334823219_tasks.mat';
end
if nargin<2 || isempty(outpath)
    outpath = '../pd_results';
end

if ~isdir(outpath)
    mkdir(outpath);
end


load(datapath,'loc_truepos','loc_falsepos','loc_falseneg','dictnames','loc_trueposv14','loc_falseposv14','loc_falsenegv14');

%calculate precision
precision_perclass = loc_truepos./(loc_truepos+loc_falsepos);
precision_perclassv14 = loc_trueposv14./(loc_trueposv14+loc_falseposv14);
%calculate recall
recall_perclass = loc_truepos./(loc_truepos+loc_falseneg);
recall_perclassv14 = loc_trueposv14./(loc_trueposv14+loc_falsenegv14);
%calculate hamming score - This isn't right because we are now double
%counting false numbers (e.g. not all "0" hamming scores are created equal)
hamming_perclass = loc_truepos./(loc_truepos+loc_falsepos+loc_falseneg);
%hamming_perclassv14 = loc_trueposv14./(loc_trueposv14+loc_falseposv14+loc_falsenegv14);


%create the file 
%headers = [{'metric'},dictnames];
% headers = {'location','precision','recall','Hamming score'};
headers = {'location','precision','recall'};


savepath = [outpath,filesep,'precision_recall2_all.txt'];
fid = fopen(savepath,'w');
fprintf(fid,'%s,',headers{1:end-1});
fprintf(fid,'%s\n',headers{end});
for i = 1:length(precision_perclass)
%     fprintf(fid,'%s,%3.2f,%3.2f,%3.2f\n',dictnames{i},precision_perclassv14(i),recall_perclassv14(i),hamming_perclassv14(i));
    fprintf(fid,'%s,%3.2f,%3.2f\n',dictnames{i},precision_perclass(i),recall_perclass(i));
%     fprintf(fid,'%3.2f\nrecall,',precision_perclass(end));
%     fprintf(fid,'%3.2f,',recall_perclass(1:end-1));
%     fprintf(fid,'%3.2f\nHamming score,',recall_perclass(end));
%     fprintf(fid,'%3.2f,',hamming_perclass(1:end-1));
%     fprintf(fid,'%3.2f',hamming_perclass(end));
end
% mean_precisionv14 = mean(precision_perclassv14(~isnan(recall_perclassv14)));
% mean_recallv14 = mean(recall_perclassv14(~isnan(recall_perclassv14)));
% mean_hammingv14 = mean(hamming_perclassv14(~isnan(recall_perclassv14)));
% fprintf(fid,'%s,%3.2f,%3.2f,%3.2f\n','Average',mean_precisionv14,mean_recallv14,mean_hammingv14);
mean_precision = mean(precision_perclass(~isnan(recall_perclass)));
mean_recall = mean(recall_perclass(~isnan(recall_perclass)));
mean_hamming = mean(hamming_perclass(~isnan(recall_perclass)));
fprintf(fid,'%s,%3.2f,%3.2f\n','Average',mean_precision,mean_recall);


fclose(fid);

