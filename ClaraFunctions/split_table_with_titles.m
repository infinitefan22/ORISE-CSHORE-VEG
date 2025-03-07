% function result = eliminate_empty_rows(tablename)
% tablename = tablename(~cellfun('isempty',tablename)) ;
% result = tablename ;
% end

function res = split_table_with_titles(titles, data)
n= 1 %indicating title rows
for i = 1:height(data)
    if data{i,1}
