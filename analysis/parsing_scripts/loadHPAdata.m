function [v14_hash,hpasol_hash,if_imgs_data,headers] = loadHPAdata(if_images_path,dictnames,dict_translator)
fid=fopen(if_images_path);

headers = textscan(fid, '%s',1, 'delimiter', '\n');

headers = strsplit(headers{1}{1},',');
numcolumns = length(headers);
loc_header_ind = strcmpi(headers,'locations');
unspec_header_ind = strcmpi(headers,'unspecific');

linestr = [repmat('%q ',1,numcolumns-1),'%q'];

if_imgs_data = textscan(fid,linestr,'Delimiter',',','HeaderLines',1);
fclose(fid);

dict_map = java.util.HashMap;
for i = 1:length(dictnames)
    dict_map.put(dictnames{i},i);
end

% dict_map.put('Nucleoli (Fibrillar center)',dict_map.get('Nucleoli (fib center)'));
dict_translator.put('Nucleoli (Fibrillar center)',dict_translator.get('Nucleoli (fib center)'));
dict_translator.put('Cytoskeleton (Intermediate filaments)',dict_translator.get('Intermediate filaments'));
dict_translator.put('Cytoskeleton (Microtubule end)',dict_translator.get('Microtubule end'));
dict_translator.put('Cytoskeleton (Microtubules)',dict_translator.get('Microtubules'));
dict_translator.put('Cytoskeleton (Actin filaments)',dict_translator.get('Actin filaments'));
dict_translator.put('Cytoskeleton (Cytokinetic bridge)',dict_translator.get('Cytokinetic bridge'));
dict_translator.put('NULL',dict_map.get('Negative'));
% dict_map.put('Cytoskeleton (Actin filaments)',dict_map.get('Nucleoli (fib center)'));

if_imgs_lastversion = if_imgs_data{16};
% if_imgs_lastversion = cellfun(@(x) str2num(x),if_imgs_data{16},'UniformOutput',false);
hpa_filename = cellfun(@(x) strsplit(x,'/'),if_imgs_data{1},'UniformOutput',false);
hpa_code = cellfun(@(x) x(end),hpa_filename);
is_v14 = cell2mat(cellfun(@(x) strcmp(x,'14'),if_imgs_lastversion,'UniformOutput',false));
v14_hash = java.util.HashSet;
hpasol_hash = java.util.HashMap;
hpa_soltot = zeros(size(dictnames));
hpa_solperc_prev = zeros(size(dictnames));
perc_calc = false;
hpa_solperc_diff = zeros(length(is_v14),1);
for i = 1:length(is_v14)
    
    currnames = strsplit(if_imgs_data{loc_header_ind}{i},',');
    currnames_num = zeros(size(dictnames));
    is_unspecific = logical(str2num(if_imgs_data{unspec_header_ind}{i}));
    if is_unspecific
        currnames_num(dict_map.get('Unspecific')) = currnames_num(dict_map.get('Unspecific'))+1;
    else
%         if strcmpi(if_imgs_data{6}{i},'NULL')
%             holdup = 1;
%         end
%         
        for j = 1:length(currnames)
            
                    if strfind(currnames{j},'Nucleoplasm')
                        holdup = 1;
            %             tmpname = strsplit(currnames{j},'(');
            %             currnames{j} = tmpname{end}(1:end-1);
                    end
            %         if strcmp(currnames{j},'Nucleoli (Fibrillar center)')
            %             checkme = 1;
            % %             currnames{j} = 'Nucleoli (fib center)';
            %         end
            %         switch currnames
            %             case 'Cytoskeleton (Intermediate filaments)'
            %                 currnames = 'Intermediate filaments';
            %             case 'Cytoskeleton (Microtubule end)'
            %                 currnames = 'Microtubule end';
            %             case 'Cytoskeleton (Microtubules)'
            %                 currnames = 'Microtubules';
            %             case 'Cytoskeleton (Actin filaments)'
            %                 currnames = 'Actin filaments';
            %             case 'Cytoskeleton (Cytokinetic bridge)'
            %                 currnames = 'Cytokinetic bridge';
            %             otherwise
            %         end
            %         currnames_num = currnames_num+strcmp(dictnames,currnames{j});

            currnames_num(dict_map.get(dict_translator.get(currnames{j}))) = currnames_num(dict_map.get(dict_translator.get(currnames{j})))+1;
        end
    end
    if is_v14(i)
        v14_hash.add(hpa_code{i});
    end
    hpasol_hash.put(hpa_code{i},currnames_num>0);
        
end

