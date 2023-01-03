clear all
%% Setup
ID = importdata('C:/Users/Bruger/OneDrive - Aarhus universitet/5th Semester/Bachelor Thesis/Data/raw_srl/names.txt');

%create folder for results
resdir = "C:\Users\Bruger\Documents\raw_srl\second level\";
if isfolder(resdir) == 1
    rmdir(resdir, "s")
end

%%

for con = 1:12
    
    resdir_con = append(resdir,sprintf('con_%04d',con));

%     if isfolder(resdir_con) == 1
%         try
%             rmdir(resdir_con, "s"); %remove old files
%         catch
%             warning('Problem deleting files in folder.');
%         end
%     end
    %%
    mkdir(resdir_con);
    
    resdir_con = cellstr(resdir_con);% put the name in a cell so that it is accepted by spm_jobman
    
    %%
    files = {};
    for i = 1:size(ID,1)
       files = [files ; append('C:\Users\Bruger\Documents\raw_srl\TNU_PRSSI_',ID{i},'\',sprintf('con_%04d.nii,1',con))];
    
    end
    
    
    matlabbatch{1}.spm.stats.factorial_design.dir = resdir_con;
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = files;
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov.files = {'C:\Users\Bruger\Documents\raw_srl\posterior_list_alpha_beta_omega.txt'};
    matlabbatch{1}.spm.stats.factorial_design.multi_cov.iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.multi_cov.iCC = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Main effect';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 0 0];
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'alpha';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 1 0];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'beta';
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 0 1];
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.delete = 1;

    % Run
    spm('defaults', 'FMRI');
    spm_jobman('run', matlabbatch);
end