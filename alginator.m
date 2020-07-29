
[~,~,GNRD] = xlsread('../structures.xlsx', 'GNRD');

fileID = fopen('GNRD.pdb', 'w');

[rows, ~] = size(GNRD);
num_atoms = rows - 1;

for i = 1:num_atoms
    
    if ~isnan(GNRD{i+1,4})
        fprintf(fileID, ['ATOM    ', '%d  ', GNRD{i+1,2}, '   1    ', num2str(GNRD{i+1,4}), '  ',...
            num2str(GNRD{i+1,5}), '  ', num2str(GNRD{i+1,6}), '  1.00  0.00\n'], i);
    end
end

fclose(fileID);