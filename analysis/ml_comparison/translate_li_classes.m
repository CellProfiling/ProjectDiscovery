function li_translated = translate_li_classes(li_datamat,class_names_li,dictnames)


li_translated = zeros(size(li_datamat,1),length(dictnames));
for i = 1:length(class_names_li)
    curr_name = 0;
    switch class_names_li{i}
        case 'loc_Centrosome'
            curr_name = strcmpi(dictnames,'Centrosome');
        case 'loc_Cytoplasm'
            curr_name = strcmpi(dictnames,'Cytosol');
        case 'loc_Cytoskeleton_actin_filaments'
            curr_name = strcmpi(dictnames,'Actin filaments');
        case 'loc_Cytoskeleton_intermediate_filaments'
            curr_name = strcmpi(dictnames,'Intermediate filaments');
        case 'loc_Cytoskeleton_microtubules'
            curr_name = strcmpi(dictnames,'Microtubules');
        case 'loc_Endoplasmic_reticulum'
            curr_name = strcmpi(dictnames,'Endoplasmic reticulum');
        case 'loc_Golgi'
            curr_name = strcmpi(dictnames,'Golgi apparatus');
        case 'loc_Mitochondria'
            curr_name = strcmpi(dictnames,'Mitochondria');
        case 'loc_Nuclei'
            curr_name = strcmpi(dictnames,'Nucleus');
        case 'loc_Nuclei_but_not_nucleoli'
            curr_name = strcmpi(dictnames,'Nucleoplasm');
        case 'loc_Nucleoli'
            curr_name = strcmpi(dictnames,'Nucleoli');
        case 'loc_Plasma_membrane'
            curr_name = strcmpi(dictnames,'Plasma membrane');
        case 'loc_Vesicles'
            curr_name = strcmpi(dictnames,'Vesicles');
        otherwise
            error('unknown class')
    end
    li_translated(:,curr_name) = li_datamat(:,i);
end