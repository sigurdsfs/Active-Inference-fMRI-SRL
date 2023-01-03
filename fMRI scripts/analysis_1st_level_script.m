%%
clear all
%%
ID = importdata('C:/Users/Bruger/OneDrive - Aarhus universitet/5th Semester/Bachelor Thesis/Data/raw_srl/names.txt');
dir = 'C:\Users\Bruger\Documents\raw_srl\';

for n = 1:size(ID,1)

    subject_dir = append(dir, 'TNU_PRSSI_', ID{n});
    
    matlabbatch{1}.spm.stats.fmri_spec.dir = {subject_dir};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2.5;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    %% Scans
    scans = {};
    for i = 1:294
        scans = [scans ; append(subject_dir, '\scandata\swSRL_fMRI.nii,',num2str(i))];
    end
    %scans = cellstr(scans);

    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = scans;
    %%
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {append('C:\Users\Bruger\Documents\raw_srl\Paradigm_files\paradigm_file_valid_invalid_trials_',ID{n},'.mat')}; %find matching paradigm file
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {append('C:\Users\Bruger\Documents\raw_srl\TNU_PRSSI_',ID{n},'\scandata\rp_SRL_fMRI.txt')};
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Cue';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 0 0 0 0 0 0];
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Decision';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 0 0 1 0 0 0];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'Outcome';
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 0 0 0 1 0 0];
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'Cue_tonic';
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0 1];
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'cue_phasic';
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [0 0 1];
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'outcome_tonic';
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [0 0 0 0 0 1];
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = 'outcome_phasic';
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [0 0 0 0 0 0 1];
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.delete = 0;

    % Run
    spm('defaults', 'FMRI');
    spm_jobman('run', matlabbatch);
    
    disp('###########################')
    fprintf("%03d done with \n",n)
    disp('##########################')

    %clear matlabbatch subject_dir scans
end