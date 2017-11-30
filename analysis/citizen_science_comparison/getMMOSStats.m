function [number_of_seconds,number_of_players_tot,...
    mean_classifications_per_player,median_classifications_per_player,...
    players_with_no_tasks] = getMMOSStats(mmos_stats)

fid = fopen(mmos_stats);
mmos_data = textscan(fid,'%q %f','Delimiter',',');
fclose(fid);

number_of_seconds = mmos_data{2}(1);
number_of_players_tot = mmos_data{2}(2);
mean_classifications_per_player = mmos_data{2}(3);
median_classifications_per_player = mmos_data{2}(4);
players_with_no_tasks = mmos_data{2}(5);