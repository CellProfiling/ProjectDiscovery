function [is_cell_line] = get_selected_cellline(originalCode,exp_cellline_hash,cell_line)



if nargin<3 || isempty(cell_line)
    is_cell_line = ones(size(originalCode));
    return
end


is_cell_line = false(length(originalCode),1);

cell_line_hash = java.util.HashSet;

if ~iscell(cell_line)
    cell_line = {cell_line};
end

for j = 1:length(cell_line)
    cell_line_hash.add(cell_line{j});
end

for j = 1:length(originalCode)
    exp_cell_line = exp_cellline_hash.get(originalCode{j});
    is_cell_line(j) = cell_line_hash.contains(exp_cell_line);
end