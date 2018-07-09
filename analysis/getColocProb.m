function [coloc_nums,coloc_probs,class_probs,if_categories,class_nums] = getColocProb(if_images_path,v14_flag)

if nargin<1 || isempty(if_images_path)
    if_images_path = '../hpa_results/IF_images_13062016.csv';
end

if nargin<2 || isempty(v14_flag)
   v14_flag = false; 
end


%for the probability of figure 5 - cytoscape. What co-localizations are
%significantly enriched. 

%read IF_images
fid=fopen(if_images_path);

headers = textscan(fid, '%s',1, 'delimiter', '\n');

headers = strsplit(headers{1}{1},',');
loc_col = strcmpi(headers,'locations');
version_col = strcmpi(headers,'latest_version');
numcolumns = length(headers);

linestr = [repmat('%q ',1,numcolumns-1),'%q'];

if_imgs_data = textscan(fid,linestr,'Delimiter',',');
% if_imgs_data = textscan(fid,linestr,'Delimiter',',','HeaderLines',1);

if_imgs_locs = if_imgs_data{loc_col};
if v14_flag
    is_v14 = cell2mat(cellfun(@(x) strcmp(x,'14'),if_imgs_data{version_col},'UniformOutput',false));
end


if_loc_strings = strjoin(if_imgs_locs,',');
if_categories = unique(strsplit(if_loc_strings,','));
num_categories = length(if_categories);
num_if_exps = size(if_imgs_locs,1);

if_locmat = zeros(num_if_exps,num_categories);

%build locmat
tot_count = 0;
for i = 1:num_if_exps
    if v14_flag && ~is_v14(i)
        continue
    end
   curr_locs = strsplit(if_imgs_locs{i},',');
   if isempty(curr_locs)
       continue;
   end
   tot_count = tot_count+1;
   for j = 1:num_categories
       if any(strcmp(if_categories{j},curr_locs))
           if_locmat(i,j) = 1;
       end
   end
end

if_locmat = logical(if_locmat);
num_annotations = sum(if_locmat,2);
colocexps = if_locmat(num_annotations>1,:);
has_data = num_annotations>0;
coloc_nums = zeros(num_categories,num_categories);
coloc_probs = zeros(num_categories,num_categories);
for i = 1:num_categories
    currexp = if_locmat(if_locmat(:,i),:);
    coloc_nums(i,:) = sum(currexp,1);
%     coloc_probs(i,:) = coloc_nums(i,:)./size(currexp,1);
    coloc_probs(i,:) = coloc_nums(i,:)./size(currexp,1);
    
end

if v14_flag
    class_probs = diag(coloc_nums)./sum(is_v14.*has_data);%sum(num_if_exps);
    class_nums = sum(is_v14.*has_data);
else
    class_probs = diag(coloc_nums)./sum(has_data);
    class_nums = sum(has_data);
end
