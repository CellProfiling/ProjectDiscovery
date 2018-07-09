function loccat_heatmap(loccat_table_path,if_images_path)

fid = fopen(loccat_table_path);


dnn_data = xlsread('./binom_transfer_pseudo_gamer_v14.xls');
dnn_numbers = csvread('./binom_transver_pseudo_gamer_v14_counts.txt',1);