function [X,Y] = readGeneData(drugClass)
%% Return already prepared data, if possible.

out_x_path = strcat(drugClass, '_X.txt');
out_y_path = strcat(drugClass, '_Y.txt');

if exist(out_x_path, 'file') && exist(out_y_path, 'file')
    X = readtable(out_x_path);
    Y = readtable(out_y_path);
    return
end

%% Download the raw data file, if necessary.

url_base = 'http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/DATA/';
data_path = strcat(drugClass, '_DATA.txt');
data_url = strcat(url_base, data_path);

if ~(exist(data_path, 'file') == 2)
    raw_data = urlread(data_url);
    fid = fopen(data_path, 'w');
    fwrite(fid, raw_data);
    fclose(fid);
end

%% Parse the data file

w = warning('off', 'MATLAB:table:ModifiedVarnames');
data = readtable(data_path, 'Delimiter', '\t', 'TreatAsEmpty', {'', 'NA'});
warning(w);

% Standard mutations
muts = [char('A'-1+(1:26)) 'i' 'd'];

% Columns 1:3 are labels
% Columns 4:k-1 are drug resistance measurements
% Columns k:end are mutations at different positions
k = find(strcmp(data.Properties.VariableNames,'P1'), 1);
drug_cols = 4:(k-1);
pos_cols = k:size(data,2);

drug_names = data.Properties.VariableNames(drug_cols);
pos_names = data.Properties.VariableNames(pos_cols);

% Remove samples with any non-standard mutations or error flags.
bad_muts = cellfun(@isempty, regexp(data{:,pos_cols}, '^(\.|\-|[A-Zid]+)$'));
bad_rows = any(bad_muts, 2);
data(bad_rows,:) = [];

%% Construct the design matrix X and the response vector Y

% X(i,j) = i-th patient sample, j-th mutation+position combination (0 or 1)
% Y(i,j) = i-th patient sample, j-th drug (drug restance measurement)

n = size(data, 1);
X = zeros(n, length(muts), length(pos_cols));
for i = 1:n
    for k = 1:length(pos_cols)
        j = ismember(muts, data{i,pos_cols(k)}{1});
        X(i,j,k) = 1;
    end
end
X = reshape(X, n, []);

X_names = cell(length(muts), length(pos_cols));
for j = 1:length(muts)
    for k = 1:length(pos_cols)
        X_names(j,k) = strcat(pos_names(k), '_', muts(j));
    end
end
X_names = reshape(X_names, 1, []);

% Remove any mutation+position combos that never appear in the data.
empty_cols = find(sum(X) == 0);
X(:,empty_cols) = [];
X_names(:,empty_cols) = [];

% Package X and Y as tables.
X = array2table(X, 'VariableNames', X_names);
Y = data(:,drug_cols);

%% Write X and Y to disk

writetable(X, out_x_path);
writetable(Y, out_y_path);

end