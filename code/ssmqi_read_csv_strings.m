function data = ssmqi_read_csv_strings(path)
% Read an SSMQI CSV without delimiter or date-type auto-detection surprises.

options = detectImportOptions(path,'Delimiter',',','VariableNamingRule','preserve');
options = setvartype(options,options.VariableNames,'string');
data = readtable(path,options);
end
