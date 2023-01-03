clear all 
%%
data_path = 'C:\Users\Bruger\Documents\raw_srl';
nt = 294;
%%
ID = importdata('C:/Users/Bruger/OneDrive - Aarhus universitet/5th Semester/Bachelor Thesis/Data/raw_srl/names.txt');


%% Check if all files are preprocessed. 
n_bad = 0;
for n = 1:size(ID,1)
    if isfile(append(data_path,'\TNU_PRSSI_',ID{n},'\scandata\swSRL_fMRI.nii'))
        sprintf('%s ok', ID{n})
    else
        sprintf('%s not ok', ID{n})
        n_bad = n_bad + 1;
    end
end

sprintf("%d are not preprocessed", n_bad)
%% Check for number of participants who should be excluded based on too much rotation or movement. 

for n = 1:size(ID,1)
    subj(n).ID = ID{n};
    subj(n).excluded = 0;
end

p_excluded = 0;

fd_max_trans = zeros(size(ID,1),1);
fd_max_rotat = zeros(size(ID,1),1);

fd_mean_trans = zeros(size(ID,1),1);
fd_mean_rotat = zeros(size(ID,1),1);

for subji = 1:size(ID,1)                                % For loop on number od subjects                                              % Exclude subjects not present in other directory

    % Get rest folder

    rp_rest = load(append(data_path,'\TNU_PRSSI_',ID{subji},'/scandata/rp_SRL_fMRI.txt'));              % input the text file generated at the Realignment stage
   
    rp_diff_trans = diff(rp_rest(:,1:3));                      % The first 3 parameters tell the displacement in x,y, and z direction in mm. Here the vector difference operator is used to get the derivative of vector. (framewise difference)
    rp_diff_rotat = diff(rp_rest(:,4:6)*180/pi);               % The last  3 parameters tell the rotation values pitch, yaw and roll in radians (here we convert them to degrees)
   
    multp = rp_diff_trans*rp_diff_trans';                      % diagonals will represent dx^2 + dy^2 + dz^2 
    fd_trans = sqrt(diag(multp));                              % sqrt (dx^2 + dy^2 + dz^2)
   
    multp = rp_diff_rotat*rp_diff_rotat';                      % diagonals will represent dpitch^2 + dyaw^2 + droll^2
    fd_rotat = sqrt(diag(multp));                              % sqrt (dpitch^2 + dyaw^2 + droll^2)
   
    fd_max_trans(subji,1) = max(fd_trans);                     % Get the maximum traslational framewize displacement
    fd_max_rotat(subji,1) = max(fd_rotat);                     % Get the maximum rotational framewize displacement
   
    fd_mean_trans(subji,1) = mean(fd_trans);                   % Get the mean traslational framewize displacement
    fd_mean_rotat(subji,1) = mean(fd_rotat);                   % Get the mean rotational framewize displacement
    
    if fd_max_trans(subji)>2 || fd_max_rotat(subji)>2 || fd_mean_trans(subji)>0.3 || fd_mean_rotat(subji)>0.3
        subj(subji).excluded = 1;
        p_excluded = p_excluded+1
    
    end
end
sprintf("%d should be excluded due to too much movement", p_excluded)

%% Setup data for paradigm file creation
all_data = {};
%----- LOOP through all ID's -----
for n = 1:size(ID,1)
    path = append(data_path,"\TNU_PRSSI_",ID(n),"\behavior\srl\PRSSI_",ID(n),".mat");
    data = importdata(path);
    
    %------ Create a green selected? col (0 = no, 1 = yes)--------------
    for i = 1:size(data.alldata,1)
        if data.alldata(i,1) == 0 & data.alldata(i,4) == 3
            data.alldata(i,16) = 1;
        
        elseif data.alldata(i,1) == 1 & data.alldata(i,4) == 2 
            data.alldata(i,16) = 1;
        
        else 
            data.alldata(i,16) = 0;
        end 
    end 

    % -------Create a valid column--------------- (0 = no, 1 = yes)
    for i = 1:size(data.alldata,1)
        data.alldata_columns(1,17) = {'Valid'};
        if data.alldata(i,8) < 0
            data.alldata(i,17) = 0;
        else 
            data.alldata(i,17) = 1;
        end 
    end
    data.alldata_columns{16} = 'Select_Green';
    
    
    % ------ Turn into structure and then table----------- 
    data_struct = struct(data.alldata_columns{1}, data.alldata(:,1), ...
    data.alldata_columns{2}, data.alldata(:,2), ...
    data.alldata_columns{3}, data.alldata(:,3) - data.startScan.GetSecs, ... % T-reaction
    data.alldata_columns{4}, data.alldata(:,4), ...
    data.alldata_columns{5}, data.alldata(:,5), ...
    data.alldata_columns{6}, data.alldata(:,6), ...
    data.alldata_columns{7}, data.alldata(:,7) - data.startScan.GetSecs, ... % T-cue
    data.alldata_columns{8}, data.alldata(:,8) - data.startScan.GetSecs, ... % T-decision
    data.alldata_columns{9}, data.alldata(:,9) - data.startScan.GetSecs, ... % T-target
    data.alldata_columns{10}, data.alldata(:,10) - data.startScan.GetSecs, ... % T-keypress (BUT THIS DOESNT MAKE SENSE).... so we use the other 3 cue,decision & target.
    data.alldata_columns{11}, data.alldata(:,11), ...
    data.alldata_columns{12}, data.alldata(:,12), ...
    data.alldata_columns{13}, data.alldata(:,13), ...
    data.alldata_columns{14}, data.alldata(:,14), ...
    data.alldata_columns{15}, data.alldata(:,15), ...
    data.alldata_columns{16}, data.alldata(:,16), ...
    data.alldata_columns{17}, data.alldata(:,17), ...
    'ID', repmat(ID{n},160,1));

    
    data_table = struct2table(data_struct);

    % ----add all subject tables together----
    if n == 1
        data_combined = data_table;
    else 
        data_combined = [data_combined ; data_table];
    end 
    
    data_combined.t_reaction = data_combined.t_dec - data_combined.t_cue;
    
    all_data = [all_data ; {data}];
    

end

%data_combined.ID = categorical(data_combined.ID)

save("C:\Users\Bruger\Documents\raw_srl\Database_20112022\data_combined", 'data_combined')



%% Paradigmn creation
%Loop through all participants
ID = importdata('C:/Users/Bruger/OneDrive - Aarhus universitet/5th Semester/Bachelor Thesis/Data/raw_srl/names.txt');


pardir = "C:\Users\Bruger\Documents\raw_srl\Paradigm_files";
data_combined = load("C:\Users\Bruger\Documents\raw_srl\Database_20112022\data_combined").data_combined;

MDP_fit_list = load("C:\Users\Bruger\OneDrive - Aarhus universitet\5th Semester\Bachelor Thesis\Active Inference Models\Two factor - 5 Time Steps\MDP_fit_list.mat").MDP_fit_list;

for n = 1:size(ID,1)
    fprintf("%03d started\n",n)
    %---------Onsets-------------------
    %Cue
    Onset_cue_valid = table2array(data_combined(data_combined.Valid == 1 & strcmp(data_combined.ID, ID(n)) == 1,"t_cue"));
    Onset_cue_invalid = table2array(data_combined(data_combined.Valid == 0 & strcmp(data_combined.ID, ID(n)) == 1,"t_cue"));
    
    %Decision
    Onset_decision_valid = table2array(data_combined(data_combined.Valid == 1 & strcmp(data_combined.ID, ID(n)) == 1,"t_dec"));
    
    %Onset_decision_invalid = data_combined(data_combined.Valid == 0,"t_dec"); %They
    %didn't make a decision so they therefore cannot have a timing.
    
    %Target
    Onset_target_valid = table2array(data_combined(data_combined.Valid == 1 & strcmp(data_combined.ID, ID(n)) == 1,"t_target"));
    Onset_target_invalid = table2array(data_combined(data_combined.Valid == 0 & strcmp(data_combined.ID, ID(n)) == 1,"t_target"));
    
    
    %----------Duration----------------
    % Cue
    dur_cue_valid = table2array(data_combined(data_combined.Valid == 1 & strcmp(data_combined.ID, ID(n)) == 1,"t_reaction")); % How long did it take them to make a decision
    dur_cue_invalid = repmat(1.5, 1, size(Onset_cue_invalid,1));
    % Decision
    dur_decision_valid = Onset_target_valid - Onset_decision_valid;   %repmat(0, 1, size(Onset_decision_valid,1));
    % Target
    dur_target_valid = repmat(.5, 1, size(Onset_target_valid,1));
    dur_target_invalid = repmat(.5, 1, size(Onset_target_invalid,1));
    
    %--------Pmod---------------------
    %valid
    tonic_dopamine_cue_valid = [];
    phasic_dopamine_cue_valid = [];
    
    tonic_dopmaine_outcome_valid = [];
    phasic_dopamine_outcome_valid = [];
    
    %invalid
    tonic_dopamine_cue_invalid = [];
    phasic_dopamine_cue_invalid = [];
    
    tonic_dopmaine_outcome_invalid = [];
    phasic_dopamine_outcome_invalid = [];

    for i = 1:size(MDP_fit_list{1},2) %loop through all trials

        if any(MDP_fit_list{n}(i).s(2,:) ~= 1)
            
            %Find indexes for the predicted dopamine responses. 
            cue_index = 1:(find(MDP_fit_list{n}(i).s(2,:) ~= 1)-1)*16-1; % Index from (start:decison) message passes
            outcome_index = (find(MDP_fit_list{n}(i).s(2,:) ~= 1)-1)*16 : find(MDP_fit_list{n}(i).s(2,:) ~= 1)*16; %(deicion:decision +16) message passes
    
            %append to list
            tonic_dopamine_cue_valid =  [tonic_dopamine_cue_valid ; mean(MDP_fit_list{n}(i).wn(cue_index))];
            phasic_dopamine_cue_valid = [phasic_dopamine_cue_valid ; mean(MDP_fit_list{n}(i).dn(cue_index))];
    
            tonic_dopmaine_outcome_valid = [tonic_dopmaine_outcome_valid ; mean(MDP_fit_list{n}(i).wn(outcome_index))];
            phasic_dopamine_outcome_valid = [phasic_dopamine_outcome_valid ; mean(MDP_fit_list{n}(i).dn(outcome_index))];
       
        else 
            %Find indexes for the predicted dopamine responses.
            cue_index = 1:47; % Index from (start:decison) message passes
            outcome_index = 48:64; %(deicion:decision +16) message passes
    
            %append to list
            tonic_dopamine_cue_invalid =  [tonic_dopamine_cue_invalid ; mean(MDP_fit_list{n}(i).wn(cue_index))];
            phasic_dopamine_cue_invalid = [phasic_dopamine_cue_invalid ; mean(MDP_fit_list{n}(i).dn(cue_index))];
    
            tonic_dopmaine_outcome_invalid = [tonic_dopmaine_outcome_invalid ; mean(MDP_fit_list{n}(i).wn(outcome_index))];
            phasic_dopamine_outcome_invalid = [phasic_dopamine_outcome_invalid ; mean(MDP_fit_list{n}(i).dn(outcome_index))];
    
        end 
     
    
    end

    %Put it into proper pmod structure
    pmod=struct('name',{''},'param',{''},'poly',{''});
    
    % Names
    pmod(1).name{1}= 'Tonic Dopamine'; % Cue
    pmod(1).name{2} = 'Phasic Dopamine';
    pmod(3).name{1}= 'Tonic Dopamine'; % Outcome
    pmod(3).name{2} = 'Phasic Dopamine';
    
    % Value   
    pmod(1).param{1} = zscore(tonic_dopamine_cue_valid); %Cue
    pmod(1).param{2} = zscore(phasic_dopamine_cue_valid);
    pmod(3).param{1} = zscore(tonic_dopmaine_outcome_valid); %Outcome 
    pmod(3).param{2} = zscore(phasic_dopamine_outcome_valid);
    
    % Order of fit
    pmod(1).poly{1}= 1;
    pmod(1).poly{2}= 1;
    pmod(3).poly{1}= 1;
    pmod(3).poly{2}= 1;
    
    
    if size(data_combined(data_combined.Valid == 0 & strcmp(data_combined.ID, ID(n)) == 1,:),1) == 1 %Check if invalid trails exist
        pmod(4).name{1}= 'Tonic Dopamine'; % Cue
        pmod(4).name{2} = 'Phasic Dopamine';
        pmod(5).name{1}= 'Tonic Dopamine'; % Outcome
        pmod(5).name{2} = 'Phasic Dopamine';
        
        % Value   
        pmod(4).param{1} = zscore(tonic_dopamine_cue_invalid); %Cue
        pmod(4).param{2} = zscore(phasic_dopamine_cue_invalid);
        pmod(5).param{1} = zscore(tonic_dopmaine_outcome_invalid); %Outcome 
        pmod(5).param{2} = zscore(phasic_dopamine_outcome_invalid);
        
        % Order of fit
        pmod(4).poly{1}= 1;
        pmod(4).poly{2}= 1;
        pmod(5).poly{1}= 1;
        pmod(5).poly{2}= 1;
            

    end


    if size(data_combined(data_combined.Valid == 0 & strcmp(data_combined.ID, ID(n)) == 1,:),1) == 0 % If invalid trials does not exist then we do not include
       names = {'ValidTrialCue' 'ValidTrialDecision' 'ValidTrialsTarget'};
       onsets = {Onset_cue_valid Onset_decision_valid Onset_target_valid};
       durations = {dur_cue_valid' dur_decision_valid' dur_target_valid'};
       pmod  = pmod;

    else % Or else don't.
       names = {'ValidTrialCue' 'ValidTrialDecision' 'ValidTrialsTarget'  'InvalidTrialCue' 'InvalidTrialsTarget'};
       onsets =   {Onset_cue_valid Onset_decision_valid Onset_target_valid Onset_cue_invalid Onset_target_invalid};
       durations = {dur_cue_valid' dur_decision_valid' dur_target_valid'  dur_cue_invalid' dur_target_invalid' };
       pmod = pmod;
    end
    
    
    disp("Done")
    
    paradigm_file_name=['paradigm_file_valid_invalid_trials_' ID{n}]; %number in consequtive order.
    cd(pardir)
    save(paradigm_file_name, 'durations', 'onsets', 'names', 'pmod');
    
    
    clear names onsets duration pmod Onset_cue_valid Onset_cue_invalid Onset_decision_valid Onset_target_valid Onset_target_invalid dur_cue_valid dur_cue_invalid dur_decision_valid dur_target_valid dur_target_invalid
end

%%

wrong = [];
for i = 1:160
    if all(mdp_list{3}(i).s(2,:) == 1)
        wrong = [wrong ; i];
    end
end
wrong
%%
wrong = [];
for n = 1:size(ID,1)
    for i = 1:160
       if all(MDP_fit_list{n}(i).s(2,:) == 1)
         wrong = [wrong ; i];
       end
    end
end
wrong
%% Testing the dopamine finding procedure.

%valid
tonic_dopamine_cue_valid = [];
phasic_dopamine_cue_valid = [];

tonic_dopmaine_outcome_valid = [];
phasic_dopamine_outcome_valid = [];

%invalid
tonic_dopamine_cue_invalid = [];
phasic_dopamine_cue_invalid = [];

tonic_dopmaine_outcome_invalid = [];
phasic_dopamine_outcome_invalid = [];

for i = 1:size(MDP_fit_list{1},2) %loop through all trials

    if any(MDP_fit_list{n}(i).s(2,:) ~= 1)
        
        %Find indexes for the predicted dopamine responses. 
        cue_index = 1:(find(MDP_fit_list{n}(i).s(2,:) ~= 1)-1)*16-1; % Index from (start:decison) message passes
        outcome_index = (find(MDP_fit_list{n}(i).s(2,:) ~= 1)-1)*16 : find(MDP_fit_list{n}(i).s(2,:) ~= 1)*16; %(deicion:decision +16) message passes

        %append to list
        tonic_dopamine_cue_valid =  [tonic_dopamine_cue_valid ; mean(MDP_fit_list{n}(i).wn(cue_index))];
        phasic_dopamine_cue_valid = [phasic_dopamine_cue_valid ; mean(MDP_fit_list{n}(i).dn(cue_index))];

        tonic_dopmaine_outcome_valid = [tonic_dopmaine_outcome_valid ; mean(MDP_fit_list{n}(i).wn(outcome_index))];
        phasic_dopamine_outcome_valid = [phasic_dopamine_outcome_valid ; mean(MDP_fit_list{n}(i).dn(outcome_index))];
   
    else 
        %Find indexes for the predicted dopamine responses.
        cue_index = 1:47; % Index from (start:decison) message passes
        outcome_index = 48:64; %(deicion:decision +16) message passes

        %append to list
        tonic_dopmaine_cue_invalid =  [tonic_dopamine_cue_invalid ; mean(MDP_fit_list{n}(i).wn(cue_index))];
        phasic_dopamine_cue_invalid = [phasic_dopamine_cue_invalid ; mean(MDP_fit_list{n}(i).dn(cue_index))];

        tonic_dopmaine_outcome_invalid = [tonic_dopmaine_outcome_invalid ; mean(MDP_fit_list{n}(i).wn(outcome_index))];
        phasic_dopamine_outcome_invalid = [phasic_dopamine_outcome_invalid ; mean(MDP_fit_list{n}(i).dn(outcome_index))];

    end 
 

end




%%

cue_index = 1:(find(MDP_fit_list{1}(1).s(2,:) ~= 1)-1)*16
outcome_index = (find(MDP_fit_list{1}(1).s(2,:) ~= 1)-1)*16+1 : find(MDP_fit_list{1}(1).s(2,:) ~= 1)*16
%%

tonic_dopamine_cue = mean(MDP_fit_list{1}(1).wn(cue_index))
phasic_dopamine_cue = mean(MDP_fit_list{1}(1).dn(cue_index))

mean(MDP_fit_list{1}(1).wn(outcome_index))
mean(MDP_fit_list{1}(1).dn(outcome_index))

%% Random investigations of number columns
i = 48;
format long g
all_data{i}.alldata(:,7:10) - all_data{i}.startScan.GetSecs

%%
i = 1
format long g
all_data{i}.alldata(:,3) - all_data{i}.startScan.GetSecs


%% Onset times 
format long g
data.alldata(:,7:10) - data.startScan.GetSecs

%% reaction times 
format long g
data.alldata(:,3) - data.startScan.GetSecs 

%% reaction times 2 ( THIS ONE WORKS!)
format long g
(data.alldata(:,8) - data.startScan.GetSecs) - (data.alldata(:,7) - data.startScan.GetSecs)
%%
max(data_combined.t_dec - data_combined.t_cue)
%%
format long
data.alldata(:,7:10)
%%
for n = 1:size(ID,1)
    all_data{n}.reward
end