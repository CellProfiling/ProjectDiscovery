function parsed_outpath = parseTQdata_merge(inpath,outpath,if_images_path,cell_line,merge_level,parsed_outpath)
%This function parses the data from the EVE online players in Project
%Discovery as delivered from MMOS and generates metrics seen in the
%publication: Sullivan et al. 2016
%
%INPUTS:
%inpath - string to the most recent data export from MMOS. This should be a
%tab delimited file.
%
%OUTPUTS:
%matfile of intermediate results parsed from the MMOS data download
%
%
%
%Written By: Devin P Sullivan May 2015

sig_cut = 0.01;
if nargin<1 || isempty(inpath)
    inpath = '../data/mmos_data/201607121468334823219_tasks.tab';%D.Sullivan 29,06,2018%20170307-tasks3.tab';
end

if nargin <2 || isempty(outpath)
    outpath = '../results/intermediate_TQ/';
end
if nargin<3 || isempty(if_images_path)
    if_images_path = '../data/hpa_results/IF_images_13062016.csv';
end

if nargin<4 || isempty(cell_line)
    cell_line = '';
end

if ~isdir(outpath)
    mkdir(outpath)
end

%Merge level -- parameter that determines what power of 10 to merge classes
%at from their codes. 
if nargin< 5 || isempty(merge_level)
    merge_level=0;%no merging
    % merge_level = 1;%merging 101-102 (10), 111-113 (11) etc
    % merge_level = 2;%merging 101-131 (1), 201-233 (2) etc
end

[path,inpath_base,filetype] = fileparts(inpath);
if nargin<6 || isempty(parsed_outpath)
    parsed_outpath = [outpath,filesep,'parsedTQdata',inpath_base,cell_line,'_mergelevel',merge_level,'.mat'];
end


if exist(parsed_outpath,'file')
    disp('parseTQdata_merge file found. Skipping recompute');
    return
end

dictclasses = [101,102,111,112,113,121,122,123,131,201,202,203,204,211,212,213,214,215,221,222,231,232,233,301,302,303,801,901,902];
[dictclasses,ia,merged_inds] = unique(floor(dictclasses./10^merge_level));
%     dictnames = {'Nucleus','Nucleoplasm','Nuc. bodies (few)','Nuc. bodies (many)','Nuc. speckles','Nucleoli',...
%         'Nucleoli (rim)','Nucleoli (fib center)','Nuclear membrane','Cytoplasm','Aggresome','Mitochondria','R&R',...
%         'Int. filaments','Tubule ends','Microtubules','Actin','Ck. bridge','MTOC','Centrosome',...
%         'ER','Golgi','Vesicles','Junctions','Focal adhesions','Plasma Membrane','CCD','Negative','Unspecific'};
dictnames = {'Nucleus','Nucleoplasm','Nuclear bodies','Nuclear bodies (many)','Nuclear speckles','Nucleoli',...
    'Nucleoli (rim)','Nucleoli (fib center)','Nuclear membrane','Cytoplasm','Aggresome','Mitochondria','R&R',...
    'Intermediate filaments','Microtubule end','Microtubules','Actin filaments','Cytokinetic bridge',...
    'Microtubule organizing center','Centrosome',...
    'Endoplasmic reticulum','Golgi apparatus','Vesicles','Cell Junctions',...
    'Focal Adhesions','Plasma membrane','CCD','Negative','Unspecific'};

dictclasses_ml1 = [10,11,12,13,20,21,22,23,30,80,90];
dictnames_ml1 = {'Nucleus smooth','Nuclear bodies/speckles','Nucleoli_tot','Nuc membrane',...
    'Cytoplasmic structures','Cytoskeleton','MTOC','Secretory',...
    'Periphery_2','CCV_2','Unidentifiable_2'};
dictclasses_ml2 = [1,2,3,8,9];
dictnames_ml2 = {'Nuclear',...
    'Cytoplasmic',...
    'Periphery','CCV','Unidentifiable'};

dictnames2 = cell(length(dictclasses),1);
for i = 1:length(dictclasses)
    dictnames2{i} = strjoin(dictnames(merged_inds==i),',');
end
dictnames = dictnames2;
clear dictnames2


%grab the data
%tq_data = tdfread(inpath,'\t');
fid = fopen(inpath);
%This is an outdated format Updated 16,03,16
%C = textscan(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s');
headers = fgetl(fid);
headers = strsplit(headers,'\t');
formatspec = [repmat('%q\t',1,length(headers)-1),'%q'];
%added because of bad_parsing
headers = [headers(1:end-1),'updated_at2',headers(end)];
formatspec = [formatspec,'\t%q'];

%C = textscan(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%q\t%s');
C = textscan(fid,formatspec);
fclose(fid);



for i = 1:length(C)
%     eval([C{i}{1},'= C{',num2str(i),'}(2:end);'])
    eval([headers{i},'= C{',num2str(i),'};'])
end


for i = 1:length(updatedAt)
    updatedAt{i} = [updatedAt{i},'  ',updated_at2{i}];
end

clear updated_at2

%grab only the specified cell type. If none are specified, then skip this. 
[exp_cellline_hash,~] = getCellLines(if_images_path);
[is_cell_line] = get_selected_cellline(originalCode,exp_cellline_hash,cell_line);

if sum(is_cell_line)<length(originalCode)
    for i = 1:length(C)
        eval([C{i}{1},'= ',C{i}{1},'(is_cell_line);'])
        
    end
end


%rank the list based on reports
remarknum = cellfun(@(x) str2num(x),remarkCount);
[remarknumsort,remarkind] = sort(remarknum,'descend');
remarkorigsort = remarkCount(remarkind);
origIDremarksort = originalCode(remarkind);
outputremark = [{'origID','remarkCount'};origIDremarksort,remarkorigsort];
reID = fopen(['remarkSort',inpath_base,'.dat'],'w');
reformatspec = '%s\t%s\n';
for i = 1:size(outputremark,1)
    fprintf(reID,reformatspec,outputremark{i,:});
end
fclose(reID);

%rank the list based on total classifications
votecountnum = cellfun(@(x) str2num(x),classificationCount);
[vcsort,vcind] = sort(votecountnum,'descend');
origIDvcsort = originalCode(vcind);
vcorigsort = classificationCount(vcind);
outputvc = [{'origID','voteCount'};origIDvcsort,vcorigsort];
vcID = fopen(['voteSort.dat',inpath_base,'.dat'],'w');
reformatspec = '%s\t%s\n';
for i = 1:size(outputvc,1)
    fprintf(vcID,reformatspec,outputvc{i,:});
end
fclose(vcID);

%format the date-time into reasonable dates counting from the start-date
runningdates = convertDate(updatedAt);

%get the perfect match accurracy
percacc = sum(strcmp(result,solutionUnconfirmed))/size(solutionUnconfirmed,1);

%get the classifications split apart
solUCstrcell = cellfun(@(x) strsplit(x(2:end-1),','),solutionUnconfirmed,'UniformOutput',false);
resultstrcell = cellfun(@(x) strsplit(x(2:end-1),','),result,'UniformOutput',false);

%make hashmap for v14 
%read IF_images
fid=fopen(if_images_path);

headers = textscan(fid, '%s',1, 'delimiter', '\n');

headers = strsplit(headers{1}{1},',');
numcolumns = length(headers);

linestr = [repmat('%q ',1,numcolumns-1),'%q'];

if_imgs_data = textscan(fid,linestr,'Delimiter',',','HeaderLines',1);
fclose(fid);

if_imgs_lastversion = if_imgs_data{16};
% if_imgs_lastversion = cellfun(@(x) str2num(x),if_imgs_data{16},'UniformOutput',false);
hpa_filename = cellfun(@(x) strsplit(x,'/'),if_imgs_data{1},'UniformOutput',false);
hpa_code = cellfun(@(x) x(end),hpa_filename);
is_v14 = cell2mat(cellfun(@(x) strcmp(x,'14'),if_imgs_lastversion,'UniformOutput',false));
v14_hash = java.util.HashSet;
for i = 1:length(is_v14)
    if is_v14(i)
        v14_hash.add(hpa_code{i});
    end
end



%pre allocate
%constants
num_tasks = size(solUCstrcell,1);
numclasses = length(dictclasses);
%per task
solUCnum = cell(num_tasks,1);
resultsnum = cell(num_tasks,1);
nummatch = zeros(num_tasks,1);
numclass_perfcov = zeros(num_tasks,1);
numclass_perfmatch = zeros(num_tasks,1);
score_Hamming = nan(num_tasks,1);
score_Hammingv14 = nan(num_tasks,1);

%totals
match = 0;
wrong = 0;
cytospam = 0;
nextclosest = 0;
numempty = 0;
numsol_0 = 0;

%per class
tot_gamerresult = zeros(1,numclasses);
tot_gamerresultv14 = zeros(1,numclasses);
%gamerresult = zeros(1,numclasses);
tot_hparesult = zeros(1,numclasses);
tot_hparesultv14 = zeros(1,numclasses);
%hparesult = zeros(1,numclasses);
matchmat = zeros(1,numclasses);
%isclose = zeros(1,numclasses);
%notclose = zeros(1,numclasses);
overannotation = zeros(1,numclasses);
underannotation = zeros(1,numclasses);
overannotationv14 = zeros(1,numclasses);
underannotationv14 = zeros(1,numclasses);

%confusion matrices
locconfmat = zeros(numclasses,numclasses);
singlelocconfmat = zeros(numclasses,numclasses);
locconfmatv14 = zeros(numclasses,numclasses);
singlelocconfmatv14 = zeros(numclasses,numclasses);

%matrices of answers
pmat = zeros(num_tasks,numclasses);
vcmat = zeros(num_tasks,numclasses);
vpmat = zeros(num_tasks,numclasses);

%fp,fn,tp,tn
loc_falsepos = zeros(1,numclasses);
loc_falseneg = zeros(1,numclasses);
loc_truepos = zeros(1,numclasses);
loc_trueneg = zeros(1,numclasses);
loc_falseposv14 = zeros(1,numclasses);
loc_falsenegv14 = zeros(1,numclasses);
loc_trueposv14 = zeros(1,numclasses);
loc_truenegv14 = zeros(1,numclasses);

classificationCount = cellfun(@str2num,classificationCount);
tot_tasks = 0;
tot_tasksv14 = 0; 
num_hpa_tasks = zeros(1,numclasses);
%parse the data
% for i = 1:100
for i = 1:num_tasks
    
    %parse votes
    [vcmat(i,:),vpmat(i,:),pmat(i,:)] = parsevotedata(votes{i},dictclasses,merge_level,classificationCount(i));
    
    %convert HPA and gamer solutions
    solUCnum{i} = cellfun(@(x) str2double(x),solUCstrcell{i},'UniformOutput',false);
    numsolUC = length(solUCnum{i});
    resultsnum{i} = cellfun(@(x) str2double(x),resultstrcell{i},'UniformOutput',false);
    
    unmatchedHPA = solUCnum{i};
    unmatchedGamer = resultsnum{i};
    %get the number of answers per class from the HPA
    hparesult = matchNum(solUCnum{i},merge_level);
    %         for m = 1:numsolUC
    %             if isempty(solUCnum{i}{m})
    %                 continue
    %             end
    %             hpanum = matchNum(solUCnum{i}{m});
    %             if hpanum==0
    %                 continue
    %             end
    %             hparesult(hpanum) = hparesult(hpanum)+1;
    %         end
    
    %look through gamer solutions for a match to HPA
    [gamerresult] = pmat(i,:)<sig_cut;%matchNum(resultsnum{i},merge_level);
    
    if sum(hparesult)==0
        warning('no HPA solution')
        continue
    end
    if sum(gamerresult)==0
        warning('no gamer solution')
        %             loc_falseneg = hparesult;
        continue
    end
    
    
    match_vec=and(hparesult,gamerresult);
    nummatch(i) = sum(match_vec);
    numsol_0=numsol_0+sum(match_vec);
    matchmat = matchmat+match_vec;
    
    %         locconfmat(match_vec,match_vec) = locconfmat(match_vec,match_vec)+1;
%     if sum(hparesult)==1
%         singlelocconfmat(match_vec,match_vec) = singlelocconfmat(match_vec,match_vec)+1;
%     end
    
    loc_truepos = loc_truepos+match_vec;
    loc_trueneg = loc_trueneg+(~match_vec.*~hparesult.*~gamerresult);
    
    
    
    %2016,08,15 - adding Hamming Score so we can be comparible with
    %Casper's model
    score_Hamming(i) = sum(and(hparesult,gamerresult))/sum(or(hparesult,gamerresult));
        
    %find the things that are unmatched
    %         unmatchedGamer(cellfun(@(x) isempty(x),unmatchedGamer)) = [];
    %         unmatchedHPA(cellfun(@(x) isempty(x),unmatchedHPA)) = [];
    
    %remove the ones that have a bad HPA number ==0
    %         unmatchedHPA(cellfun(@(x) x==0,unmatchedHPA)) = [];
    diff_result = hparesult-gamerresult;
    if sum(hparesult)==1 && sum(gamerresult)==1
        singlelocconfmat(hparesult,gamerresult) = singlelocconfmat(hparesult,gamerresult)+1;
    end
    if sum(match_vec)>0
        checkthis = 1;
    end
    locconfmat(hparesult,gamerresult) = locconfmat(hparesult,gamerresult)+1;
    
    loc_falsepos = double(loc_falsepos)+double(diff_result<0);
    loc_falseneg = double(loc_falseneg)+double(diff_result>0);

%     mean(sum(loc_truepos)./sum([loc_truepos,loc_falseneg,loc_falsepos]))
%     xxx = score_Hamming(i)
    
    %if sum(diff_result)==0
    if all(diff_result==0)
        
    elseif max(diff_result)<=0
        overannotation = overannotation+double(diff_result<0);
    elseif min(diff_result)>=0
        underannotation = underannotation+double(diff_result>0);
    else

%         loc_falsepos = double(loc_falsepos)+double(diff_result<0);
%         loc_falseneg = double(loc_falseneg)+double(diff_result>0);
    end
    
    
    
    if nummatch(i)==0
        wrong = wrong+1;
%     elseif nummatch(i)==length(resultsnum{i})
    elseif nummatch(i)==sum(gamerresult)
        if score_Hamming(i)==1
            numclass_perfmatch(i) = nummatch(i);
        end
        numclass_perfcov(i) = nummatch(i);
    end
    
    %track the total number of classifications of each category
    tot_hparesult = tot_hparesult+hparesult;
    tot_gamerresult = tot_gamerresult+gamerresult;
    
    if v14_hash.contains(originalCode{i})
        score_Hammingv14(i) = sum(and(hparesult,gamerresult))/sum(or(hparesult,gamerresult));
        loc_falseposv14 = double(loc_falseposv14)+double(diff_result<0);
        loc_falsenegv14 = double(loc_falsenegv14)+double(diff_result>0);
        loc_trueposv14 = loc_trueposv14+match_vec;
        loc_truenegv14 = loc_truenegv14+(~match_vec.*~hparesult.*~gamerresult);
        tot_hparesultv14 = tot_hparesultv14+hparesult;
        tot_gamerresultv14 = tot_gamerresultv14+gamerresult;
        if sum(hparesult)==1 && sum(gamerresult)==1
            singlelocconfmatv14(hparesult,gamerresult) = singlelocconfmatv14(hparesult,gamerresult)+1;
        end
        locconfmatv14(hparesult,gamerresult) = locconfmatv14(hparesult,gamerresult)+1;

        if all(diff_result==0)
        elseif max(diff_result)<=0
            overannotationv14 = overannotationv14+double(diff_result<0);
        elseif min(diff_result)>=0
            underannotationv14 = underannotationv14+double(diff_result>0);
        else
            
            %         loc_falsepos = double(loc_falsepos)+double(diff_result<0);
            %         loc_falseneg = double(loc_falseneg)+double(diff_result>0);
        end
        
        tot_tasksv14 = tot_tasksv14+1;
        
        
    end
        
    
    tot_tasks = tot_tasks+1;
    
end


save(parsed_outpath)


