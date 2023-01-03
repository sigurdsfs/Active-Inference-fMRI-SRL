% Contrast
clear all
ID = importdata('C:/Users/Bruger/OneDrive - Aarhus universitet/5th Semester/Bachelor Thesis/Data/raw_srl/names.txt');


for n = 1:1 %size(ID,1)
    matlabbatch{1}.spm.stats.con.spmmat = {append('C:\Users\Bruger\Documents\raw_srl\TNU_PRSSI_',ID{n},'\SPM.mat')};
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Cue';
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'Prediction';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [0 0 0 1];
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'Outcome';
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [0 0 0 0 1];
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'Cue Dopamine';
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = [0 1 1];
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'Outcome dopamine';
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = [0 0 0 0 0 1 1];
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = 'All dopamine';
    matlabbatch{1}.spm.stats.con.consess{6}.tcon.weights = [0 1 1 0 0 1 1];
    matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{7}.tcon.name = 'Tonic dopamine';
    matlabbatch{1}.spm.stats.con.consess{7}.tcon.weights = [0 1 0 0 0 1 0];
    matlabbatch{1}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{8}.tcon.name = 'Phasic Dopamine';
    matlabbatch{1}.spm.stats.con.consess{8}.tcon.weights = [0 0 1 0 0 0 1];
    matlabbatch{1}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{9}.tcon.name = 'Cue tonic dopamine';
    matlabbatch{1}.spm.stats.con.consess{9}.tcon.weights = [0 1 0];
    matlabbatch{1}.spm.stats.con.consess{9}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{10}.tcon.name = 'Cue phasic dopamine';
    matlabbatch{1}.spm.stats.con.consess{10}.tcon.weights = [0 0 1];
    matlabbatch{1}.spm.stats.con.consess{10}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{11}.tcon.name = 'Outcome tonic dopamine';
    matlabbatch{1}.spm.stats.con.consess{11}.tcon.weights = [0 0 0 0 0 1 0];
    matlabbatch{1}.spm.stats.con.consess{11}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{12}.tcon.name = 'Outcome phasic dopamine';
    matlabbatch{1}.spm.stats.con.consess{12}.tcon.weights = [0 0 0 0 0 0 1];
    matlabbatch{1}.spm.stats.con.consess{12}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.delete = 1;

    % Run
    spm('defaults', 'FMRI');
    spm_jobman('run', matlabbatch);

    disp('###########################')
    fprintf("%03d done with \n",n)
    disp('##########################')


    clear matlabbatch
end