function [match_vec,nummatch,score_Hamming,singlelocconfmat,locconfmat,...
    loc_truepos,loc_trueneg,loc_falsepos,loc_falseneg] = assess_task(hparesult,gamerresult)

%%Variables
singlelocconfmat = zeros(length(hparesult),length(gamerresult));
locconfmat = zeros(length(hparesult),length(gamerresult));
wrong = 0;
numclass_perfmatch = 0;
numclass_perfcov = 0;

hparesult = logical(hparesult);
gamerresult = logical(gamerresult);

if size(hparesult,1)~=size(gamerresult,1)
    hparesult = hparesult';
end
match_vec=and(hparesult,gamerresult);
nummatch = sum(match_vec);
% numsol_0=numsol_0+sum(match_vec);
% matchmat = matchmat+match_vec;

%calc the true pos
loc_truepos = match_vec;
%calc the true neg
loc_trueneg = (~match_vec.*~hparesult.*~gamerresult);

%calc the Hamming score for this task
score_Hamming = sum(and(hparesult,gamerresult))/sum(or(hparesult,gamerresult));

%%Check if any are completely wrong
% if nummatch==0
if score_Hamming == 0
    wrong = 1;
    %     elseif nummatch(i)==length(resultsnum{i})
elseif nummatch==sum(gamerresult)
    if score_Hamming==1
        numclass_perfmatch = nummatch;
    end
    numclass_perfcov = nummatch;
end

%calc mismatches
diff_result = hparesult-gamerresult;

%add a point to a blanck single location confusion matrix
if sum(hparesult)==1 && sum(gamerresult)==1
    singlelocconfmat(hparesult,gamerresult) = singlelocconfmat(hparesult,gamerresult)+1;
end
%add points to the multi-label heatmap of common confusion/co-annotation
locconfmat(hparesult,gamerresult) = locconfmat(hparesult,gamerresult)+1;

%calc false positive and false negative
loc_falsepos = double(diff_result<0);
loc_falseneg = double(diff_result>0);

%Calculate over- and under-annotation
%Currently Unused
%overanntotation
% if max(diff_result)<0
%     overannotation = double(diff_result<0);
% end
%underannotation
% if min(diff_result)>0
%     underannotation = double(diff_result>0);
% end



