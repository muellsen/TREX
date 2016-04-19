function TSM = readTSMData(drugClass)
%% Download the raw data file, if necessary.

url_base = 'http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/MUTATIONLISTS/NP_TSM/';
data_path = strcat(drugClass, '_TSM.txt');
data_url = strcat(url_base, drugClass);

if ~(exist(data_path, 'file') == 2)
    raw_data = urlread(data_url);
    fid = fopen(data_path, 'w');
    fwrite(fid, raw_data);
    fclose(fid);
end

%% Parse the data file

% Keep only the mutation positions (first column).
data = readtable(data_path, 'Delimiter', '\t', 'ReadVariableNames', false);
TSM = data{:,1};

end