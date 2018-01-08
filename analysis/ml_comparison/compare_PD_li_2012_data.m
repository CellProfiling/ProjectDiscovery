function compare_PD_li_2012_data(li_datapath,PD_datapath)

if nargin<1 || isempty(li_datapath)
    li_datapath = '/Users/devinsullivan/Documents/PD_paper/Comparisons/Li_et_al_2012_svm/li_et_al_labels.csv';
end

if nargin<2 || isempty(PD_datapath)
    PD_datapath = '/Users/devinsullivan/Documents/ProjectDiscovery/results_archive/intermediate/parsedDetailedData_v3_mergelevel0.mat';
end


fid = fopen(li_datapath);
headers = textscan(fid,'%s',1);
headers = strsplit(headers{1}{1},',');
class_names_li = headers(2:end);
formattot = ['%s',repmat('%d',1,length(headers)-2),'%d'];
data = textscan(fid,formattot,'Delimiter',',');
fclose(fid);

%get the list of experiments
well_id = data{1};
%get list of results HPAv4
hpav4_labelmat = cell2mat(data(2:end));

%load PD data
load(PD_datapath,'dictnames','gamersol_tot','uniqueCodes','hpasol_tot');

%translate the Li solution key to our solution key
hpav4_translated = translate_li_classes(hpav4_labelmat,class_names_li,dictnames);

%create well id hashset
well_hash = java.util.HashMap;
for i = 1:length(well_id)
    well_hash.put(well_id{i},hpav4_translated(i,:));
end
%Get PD data for those experiments 


%convert per-image to per-well ids
img_parts = cellfun(@(x) strsplit(x,'_'),uniqueCodes,'UniformOutput',false);
uC_wells = cellfun(@(x) strjoin(x(1:2),'_'),img_parts,'UniformOutput',false);

num_imgs = length(uniqueCodes);

keep_exp = false(num_imgs,1);
for i = 1:num_imgs
    if ~well_hash.containsKey(uC_wells{i})
        continue
    else
        if sum(gamersol_tot(i,:))==0
            continue
        end
        keep_exp(i) = true;
    end

end

%extract Li data from all data
gamersol_lidata = gamersol_tot(keep_exp,:);
hpasol_lidata = hpasol_tot(keep_exp,:);
uC_wells_lidata = uC_wells(keep_exp);


%check the solutions 
tp_gamer = zeros(1,size(gamersol_lidata,2));
fp_gamer = zeros(1,size(gamersol_lidata,2));
fn_gamer = zeros(1,size(gamersol_lidata,2));
for i = 1:length(uC_wells_lidata)
    curr_sol = well_hash.get(uC_wells_lidata{i});
    tp_gamer = tp_gamer+curr_sol.*gamersol_lidata(i,:);
    fp_gamer = fp_gamer+double(~curr_sol).*gamersol_lidata(i,:);
    fn_gamer = fn_gamer+curr_sol.*double(~gamersol_lidata(i,:));    
end

