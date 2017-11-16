function [exp_cellline_hash,all_celllines] = getCellLines(if_images_path)

if nargin<1 || isempty(if_images_path)
    if_images_path = '../hpa_results/IF_images_13062016.csv';
end


%for the probability of figure 5 - cytoscape. What co-localizations are
%significantly enriched. 

%read IF_images
fid=fopen(if_images_path);

headers = textscan(fid, '%s',1, 'delimiter', '\n');

headers = strsplit(headers{1}{1},',');
numcolumns = length(headers);

linestr = [repmat('%q ',1,numcolumns-1),'%q'];

if_imgs_data = textscan(fid,linestr,'Delimiter',',','HeaderLines',1);

if_imgs_cell_line = if_imgs_data{12};
all_celllines = unique(if_imgs_cell_line);

num_FOV = length(if_imgs_data{1});

exp_cellline_hash = java.util.HashMap;
for i = 1:num_FOV
    currexpname = [if_imgs_data{2}{i},'_',if_imgs_data{3}{i},'_',if_imgs_data{4}{i},'_'];
    exp_cellline_hash.put(currexpname,if_imgs_cell_line{i});
    
    %if exphash.containsKey(well_name)
    
end
