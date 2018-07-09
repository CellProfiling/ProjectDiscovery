function compute_time_vs_score(intresultspath,mmos_results)

if nargin<1
    intresultspath = '../results_archive/intermediate/playerTrueAcc.mat';
end

if nargin<2
    mmos_results = '../data/detailedData/classifications.tab';
end

fid = fopen(mmos_results);


%formatspec = '%d\t%s\t%d\t%f\t%q\t%s\t%q\t%q';
formatspec = '%d%s%d%f%q%s%q%q';

%C = textscan(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%q\t%s');
C = textscan(fid,formatspec,'Delimiter','\t','HeaderLines',1);
fclose(fid);

% taskID = C{1};
% originalCode = C{2};
playerID = C{3};
% playerScore = C{4};
% player_result = C{5};
player_time = C{7};


clear C

% player_time = cellfun(strsplit,player_time,'UniformOutput',false);


load(intresultspath,'unique_pIDs','player_numtasks','unique_pIDs','avplayer_recall','avplayer_precision','avplayer_f1score');


player_hash = java.util.HashMap;

for i = 1:length(player_time)
    
    if mod(i,10000)==0
        disp(i)
    end
    
    time_num = str2double(player_time{i}(6:end-1));
    
    if player_hash.containsKey(playerID(i))
        tmp = player_hash.get(playerID(i));
        player_hash.put(playerID(i),[tmp;time_num]);
    else
        player_hash.put(playerID(i),time_num);
    end
    
end

save('time_v_score_tmp.mat');

median_time = zeros(length(unique_pIDs),1);

for i = 1:length(unique_pIDs)
    
    if mod(i,10000)==0
        disp(i)
    end
    
    try 
        median_time(i) = median(player_hash.get(unique_pIDs(i)));
    catch
    end
    
end

H = niceplot(median_time,avplayer_f1score,'log10 player time (seconds)',...
    'F1 score','','semilogx','.')
plot(median_time,avplayer_f1score)
