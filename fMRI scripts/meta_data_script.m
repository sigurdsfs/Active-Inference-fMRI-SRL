%% fMRI analysis

clear all

%% Display images
%data_dir = 'C:/Users/Bruger/OneDrive - Aarhus universitet/5th Semester/Bachelor Thesis/Data/raw_srl/';


data_path = 'C:\Users\Bruger\Documents\raw_srl';
ID = importdata('C:/Users/Bruger/OneDrive - Aarhus universitet/5th Semester/Bachelor Thesis/Data/raw_srl/names.txt');

path = append(data_path,'\TNU_PRSSI_',ID{1},'/scandata/SRL_fMRI.nii');

path2 ='C:/Users/Bruger/OneDrive - Aarhus universitet/5th Semester/Bachelor Thesis/Data/raw_srl/TNU_PRSSI_A0502/scandata/SRL_fMRI.nii,1';


spm_image('Display', path)

%% See structure of the subject
func_images = spm_vol(path);

func_images_4d = spm_read_vols(func_images);

%% Get TR
for n = 1:size(ID,1)
    path_temp = append(data_path,'\TNU_PRSSI_',ID{n},'/scandata/SRL_fMRI.nii');
    temp_spm_vol = spm_vol(path_temp);

    TR(n) = temp_spm_vol(1).private.timing.tspace;
end