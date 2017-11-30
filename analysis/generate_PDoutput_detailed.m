function generate_PDoutput_detailed(parsedDetailedfile,outputdir,merge_level)

if nargin<2 || isempty(outputdir)
    outputdir = '../pd_results/';
end

if ~isdir(outputdir)
    mkdir(outputdir)
end

if nargin<3 || isempty(merge_level)
    merge_level = 0;
end


%load the data
load(parsedDetailedfile,'uniqueCodes','hpasol_tot','dictnames',...
    'v14_hash','gamersol_tot','tot_gamerresult','tot_hparesult','tot_gamerresultv14','tot_hparesultv14');
num_origIDs = length(uniqueCodes);
%print out everything
% [closevotes,perf_ind,entropy_onT,entropy_offT,mse] = analyze_Votes(vcmat,pmat,votecountnum,solUCnum,originalCode,'PD_output');

%combine repeats
all_PDmatch = ones(num_origIDs,1);
maj_PDmatch = ones(num_origIDs,1);
% all_PDrounds = cell(num_origIDs,1);
% maj_PDrounds = cell(num_origIDs,1);
% any_PDrounds = cell(num_origIDs,1);
% all_PDrounds_names = cell(num_origIDs,1);
% maj_PDrounds_names = cell(num_origIDs,1);
% any_PDrounds_names = cell(num_origIDs,1);
solUC_names = cell(num_origIDs,1);
all_perf = nan(num_origIDs,1);
maj_perf = nan(num_origIDs,1);
any_perf = nan(num_origIDs,1);
uorig_v14bool = zeros(num_origIDs,1);
all_PDmatch_correct = zeros(num_origIDs,1);
all_PDmatch_score_Hamming = nan(num_origIDs,1);

numclasses = length(dictnames);
any_freq = zeros(1,numclasses);
all_freq = zeros(1,numclasses);
maj_freq = zeros(1,numclasses);
hpa_freq = zeros(1,numclasses);
totEVE_freq = tot_gamerresult;
all_PDmatch_freq = zeros(1,numclasses);
any_count = 0;
all_count = 0;
maj_count = 0;
hpa_count = 0;
EVE_count = 0;
% numrounds_task = zeros(num_origIDs,1);
for i = 1:num_origIDs
%     numrounds_task(i) = sum(currinds);
    
    if v14_hash.contains([uniqueCodes{i}])
        uorig_v14bool(i) = true;
    end
    
    %     curr_vcmat = vcmat(currinds,:);
    %     curr_pmat = pmat(currinds,:);
    %     curr_vcnum = votecountnum(currinds);
    curr_resultsnum = gamersol_tot(i,:);
    curr_resultsnum_tmp = curr_resultsnum;%cellfun(@(x) unique(floor(cell2mat(x).*10^-merge_level)),curr_resultsnum,'UniformOutput',false);
    
%         [unique_con,~,con_ind] = unique(curr_resultsnum_tmp);
    [~,~,con_ind] = unique(curr_resultsnum_tmp);
    hparesult = hpasol_tot(i,:);%matchNum(u_solUCnum{i},merge_level);
    %     curr_hpasol = u_solUCnum(i);%u_origID{1};
%     sol_cell = {cellfun(@(x) floor(x.*10^-merge_level),curr_resultsnum{1})};

%     numrounds_task(i)=sum(cellfun(@isempty,curr_resultsnum)==0);
    %     all_PDmatch(i) = 1;
%     nummatches = 1;
%     curr_freq = matchNum(curr_resultsnum{1},merge_level);
%     for j = 2:numrounds_task(i)
%         sol_cell = [sol_cell,{cellfun(@(x) floor(x.*10^-merge_level),curr_resultsnum{j})}];
%         %         curr_res = sol_cell{j};%floor([curr_resultsnum{j}{:}].*10^-merge_level);
%         %         prev_res = sol_cell{j-1};%floor([curr_resultsnum{j-1}{:}].*10^-merge_level);
%         %         if (length(curr_res)~=length(prev_res)) || ~all(curr_res==prev_res)
%         %             if ~isnan(curr_resultsnum{j}{:})
%         %                 all_PDmatch(i) = 0;
%         %             else
%         %                 all_PDmatch(i) = -1;
%         %             end
%         %         else
%         %             nummatches = nummatches+1;
%         %         end
%         curr_freq = curr_freq+matchNum(curr_resultsnum{j},merge_level);
%         %         totEVE_freq = totEVE_freq+matchNum(curr_resultsnum{j},merge_level);
%         EVE_count = EVE_count+1;
%     end
%     totEVE_freq = totEVE_freq+curr_freq;
    %     if nummatches<(length(curr_resultsnum)/2)
    nummatches = sum(con_ind==mode(con_ind));
    if nummatches<(length(curr_resultsnum)/2)
        maj_PDmatch(i) = 0;
        all_PDmatch(i) = 0;
    elseif nummatches<length(curr_resultsnum)
        all_PDmatch(i) = 0;
    else
        all_PDmatchresult = matchNum(curr_resultsnum{1},merge_level);
        all_PDmatch_freq = all_PDmatch_freq+all_PDmatchresult;
        all_PDmatch_correct(i) = all(hparesult==all_PDmatchresult);
        all_PDmatch_score_Hamming(i) = sum(and(hparesult,all_PDmatchresult))/sum(or(hparesult,all_PDmatchresult));
    end
    
    [any_labels_mat,f_ind_lab] = unique(cell2mat(sol_cell));
    %any_labels = sol_cell(f_ind_lab);
    sol_mat = cell2mat(sol_cell);
    any_labels = unique(sol_mat);
    allresult = curr_freq==size(curr_resultsnum,1);
    always_sol = dictclasses(allresult);%{};
    majresult = curr_freq>=(size(curr_resultsnum,1)/2);
    majority_sol = dictclasses(majresult);%{};
    %     for k = 1:length(any_labels)
    %
    %         %do the total solution
    % %         if all(~strcmpi(any_lables{k},u_solUCnum{i}))
    % %             any_perf(i) = 0;
    % %         end
    %
    %         %do the always solutions
    %         if sum(sol_mat==any_labels{k})==length(curr_resultsnum)
    %             always_sol = [always_sol,any_labels(k)];
    % %             if any(~cellfun(@(x) x==any_lables{k},u_solUCnum{i}))
    % %                 all_perf(i) = 0;
    % %             end
    %         end
    %
    %         %do the majority solutions
    %         if sum(sol_mat==any_labels{k})>=(length(curr_resultsnum)/2)
    %             majority_sol = [majority_sol,any_labels(k)];
    % %             if any(~strcmpi(any_lables{k},u_solUCnum{i}))
    % %                 maj_perf(i) = 0;
    % %             end
    %         end
    %
    %     end
    
    %     hparesult = matchNum(solUCnum{i});
    anyresult = curr_freq>0;%zeros(size(dictclasses));
    %     for l = 1:length(any_labels)
    %         anyresult = anyresult+any_labels{l}==dictclasses;
    %     end
    %     anyresult = matchNum(any_labels,merge_level);
    %allresult = matchNum(always_sol,merge_level);
    %     allresult = zeros(size(dictclasses));
    %     for l = 1:length(always_sol)
    %         allresult = allresult+always_sol{l}==dictclasses;
    %     end
    %     allresult = matchNum(always_sol,merge_level);
    %     majresult = matchNum(majority_sol,merge_level);
    %     majresult = zeros(size(dictclasses));
    %     for l = 1:length(majority_sol)
    %         majresult = majresult+majority_sol{l}==dictclasses;
    %     end
    %     majresult = matchNum(majority_sol,merge_level);
    
    any_perf(i) = all(hparesult==anyresult);
    any_freq = any_freq+anyresult;
    any_count = any_count+double(sum(anyresult)>0);
    all_perf(i) = all(hparesult==allresult);
    all_freq = all_freq+allresult;
    all_count = all_count+double(sum(allresult)>0);
    maj_perf(i) = all(hparesult==majresult);
    maj_freq = maj_freq+majresult;
    maj_count = maj_count+double(sum(majresult)>0);
    hpa_freq = hpa_freq+hparesult;
    hpa_count = hpa_count+double(sum(hparesult)>0);
    

    %     any_perf(i) = (length(any_lables)==length(u_solUCnum{i})) && (all(cell2mat(any_lables)==cell2mat(u_solUCnum{i})));
    %     all_perf(i) = (length(always_sol)==length(u_solUCnum{i})) && all(cell2mat(always_sol)==cell2mat(u_solUCnum{i}));
    %     maj_perf(i) = (length(majority_sol)==length(u_solUCnum{i})) && all(cell2mat(majority_sol)==cell2mat(u_solUCnum{i}));
    %     any_perf(i) = cell2mat(any_lables
    
    any_PDrounds{i} = any_labels;
    all_PDrounds{i} = always_sol;
    maj_PDrounds{i} = majority_sol;
    any_PDrounds_names{i} = dictnames(anyresult);%matchNum(any_PDrounds{i},merge_level));
    all_PDrounds_names{i} = dictnames(allresult);%matchNum(all_PDrounds{i},merge_level));
    maj_PDrounds_names{i} = dictnames(majresult);%matchNum(maj_PDrounds{i},merge_level));
    solUC_names{i} = dictnames(matchNum(u_solUCnum{i},merge_level));
    %     curr_code = originalCode(currinds);
    
end

%take only the perfect ones
make_output_dat_solonly(all_perf,uniqueCodes,solUC_names,all_PDrounds_names,[outputdir,filesep,'all_roundsPD.txt'])
make_output_dat_solonly(any_perf,uniqueCodes,solUC_names,any_PDrounds_names,[outputdir,filesep,'tot_roundsPD.txt'])
make_output_dat_solonly(maj_perf,uniqueCodes,solUC_names,maj_PDrounds_names,[outputdir,filesep,'maj_roundsPD.txt'])

%%make outputs for v14 data only
uorig_v14bool = boolean(uorig_v14bool);
all_perfv14 = all_perf(uorig_v14bool);
any_perfv14 = any_perf(uorig_v14bool);
maj_perfv14 = maj_perf(uorig_v14bool);
uniqueCodesv14 = uniqueCodes(uorig_v14bool);
solUC_namesv14 = solUC_names(uorig_v14bool);
all_PDrounds_namesv14 = all_PDrounds_names(uorig_v14bool);
any_PDrounds_namesv14 = any_PDrounds_names(uorig_v14bool);
maj_PDrounds_namesv14 = maj_PDrounds_names(uorig_v14bool);
make_output_dat_solonly(all_perfv14,uniqueCodesv14,solUC_namesv14,all_PDrounds_namesv14,[outputdir,filesep,'all_roundsv14PD.txt'])
make_output_dat_solonly(any_perfv14,uniqueCodesv14,solUC_namesv14,any_PDrounds_namesv14,[outputdir,filesep,'tot_roundsv14PD.txt'])
make_output_dat_solonly(maj_perfv14,uniqueCodesv14,solUC_namesv14,maj_PDrounds_namesv14,[outputdir,filesep,'maj_roundsv14PD.txt'])

%%make self-consistency metrics
perc_allmatch = sum(all_PDmatch)./length(all_PDmatch);
perc_allmatch_v14 = sum(all_PDmatch.*uorig_v14bool)./sum(uorig_v14bool);
perc_majmatch = sum(maj_PDmatch)./length(maj_PDmatch);
perc_majmatch_v14 = sum(maj_PDmatch.*uorig_v14bool)./sum(uorig_v14bool);
metricnames = {'all rounds self-match','all rounds self-match v14',...
    'majority rounds self-match','majority round self-match v14',...
    '"all rounds (and)" labels hpa-match','"all rounds (and)" labels hpa-match v14'...
    '"any rounds (or)" labels hpa-match','"any rounds (or)" labels hpa-match v14'...
    '"majority rounds" labels hpa-match','"majority rounds" labels hpa-match v14'};
consistency_outname = [outputdir,filesep,'consistency_PD.txt'];
fid = fopen(consistency_outname,'w');
formatSpec = '%s,%f\n';

fprintf(fid,formatSpec,metricnames{1},perc_allmatch);
fprintf(fid,formatSpec,metricnames{2},perc_allmatch_v14);
fprintf(fid,formatSpec,metricnames{3},perc_majmatch);
fprintf(fid,formatSpec,metricnames{4},perc_majmatch_v14);
fprintf(fid,formatSpec,metricnames{5},sum(all_perf)./length(all_perf));
fprintf(fid,formatSpec,metricnames{6},sum(all_perfv14)./length(all_perfv14));
fprintf(fid,formatSpec,metricnames{7},sum(any_perf)./length(any_perf));
fprintf(fid,formatSpec,metricnames{8},sum(any_perfv14)./length(any_perfv14));
fprintf(fid,formatSpec,metricnames{9},sum(maj_perf)./length(maj_perf));
fprintf(fid,formatSpec,metricnames{10},sum(maj_perfv14)./length(maj_perfv14));

fclose(fid)
%%

%[closevotes,perf_ind,entropy_onT,entropy_offT,mse] = analyze_Votes(vcmat,pmat,votecountnum,solUCnum,originalCode,'PD_output');
%take all of them
%make_output_dat(perf_ind,origIDs,vctot,entropy_onT,entropy_offT,entropy_ratioT,HPA_solnames,gamersolnames,tot_outname);
%[closevotes,perf_ind,entropy_onT,entropy_offT,mse] = analyze_Votes(vcmat,pmat,votecountnum,solUCnum,originalCode,'PD_output');
