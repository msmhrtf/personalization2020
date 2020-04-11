
cipic = load('C:\Users\jtissi16\Documents\MATLAB\TH\CIPIC\CIPIC_hrtf_database\standard_hrir_database\subject_021\hrir_final.mat');

matrix_l = cipic.hrir_l;
matrix_r = cipic.hrir_r;

mat_l = zeros(25,50);
mat_r = zeros(25,50);
for i = 1:25
    mat_l = squeeze(matrix_l(i,:,:));
    mat_r = squeeze(matrix_r(i,:,:));
    
    idx = num2str(i);
    if numel(idx) == 1 % add a zero if below 10
        file_name_l = strcat('L','0',idx,'.txt');
        file_name_r = strcat('R','0',idx,'.txt');
        dlmwrite(file_name_l, mat_l,'delimiter',' ', 'precision', 16);
        dlmwrite(file_name_r, mat_r,'delimiter',' ', 'precision', 16);
    else
        file_name_l = strcat('L',idx,'.txt');
        file_name_r = strcat('R',idx,'.txt');
        dlmwrite(file_name_l, mat_l,'delimiter',' ', 'precision', 16);
        dlmwrite(file_name_r, mat_r,'delimiter',' ', 'precision', 16);
    end
    
    
    
end

