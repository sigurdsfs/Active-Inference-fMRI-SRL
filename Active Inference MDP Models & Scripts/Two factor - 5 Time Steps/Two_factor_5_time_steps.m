%% Model with two hidden state factors (5 Time steps) 
% slower or faster responses are modelled differently. 
clear all 
close all 
%%
%---- Start---
D{1} = [1 0]'; % Green Better or Yellow Better 
 
D{2} = [1 0 0]'; %start, Choose green , choose yellow

d{1} = [.5 .5]';
d{2} = [1 0 0]'; 
%Likelihood

Ns = [length(D{1}) length(D{2})];

pWin = .8;


A{1}(:,:,1) = [1 1;  % null 
               0 0;  % loss 
               0 0]; % Win

A{1}(:,:,2) = [0      0;     % Null        
               1-pWin pWin;  % Loss
               pWin 1-pWin]; % Win

A{1}(:,:,3) = [0      0;     % Null
               pWin 1-pWin;  % Loss
               1-pWin pWin]; % Win

for i = 1:Ns(2) 

    A{2}(i,:,i) = [1 1];

end

a{1} = A{1} *200;
a{2} = A{2} *200;

% Transistion
B{1}(:,:,1) = [1 0;  % 'Green Better' Context
               0 1]; % 'Yellow Better' Context
 

B{2}(:,:,1) = [1 1 1; %Start 
               0 0 0; % Green 
               0 0 0];% Yellow


B{2}(:,:,2) = [0 0 0; %Start 
               1 1 1; %Green 
               0 0 0];%Yellow

B{2}(:,:,3) = [0 0 0; %Start
               0 0 0; %Green
               1 1 1];%Yellow 


% Preference
No = [size(A{1},1) size(A{2},1)];

T = 4;
C{1} = zeros(No(1),T); %Wins/Losses

C{2} = zeros(No(2),T); %Observed

rs = 6;
la = 1;

C{1}(1,2:T) = 0; %Null
C{1}(2,2:T) = -la; %loss
C{1}(3,2:T) = rs; %Win


%shallow Policies 
Np = 3;
Nf = 2;

%U = ones(1,Np,Nf);
%U(:,:,1) = [1 1 1];

%U(:,:,2) = [1 2 3];


Np = (T-1)*2+1; % Number of policies
Nf = 2; % Number of state factors

%V = ones(T-1,Np,Nf);

V(:,:,1) = ones(T-1,Np);

V(:,:,2) = [1 2 1 1 3 1 1; 
            1 1 2 1 1 3 1;
            1 1 1 2 1 1 3];
              
              


% Eta: learning rate (0-1) controlling the magnitude of concentration parameter
% updates after each trial (if learning is enabled).

    eta = 0.5; % By default we here set this to 0.5, but try changing its value  
               % to see how it affects model behavior

% Omega: forgetting rate (0-1) controlling the reduction in concentration parameter
% magnitudes after each trial (if learning is enabled). This controls the
% degree to which newer experience can 'over-write' what has been learned
% from older experiences. It is adaptive in environments where the true
% parameters in the generative process (priors, likelihoods, etc.) can
% change over time. A low value for omega can be seen as a prior that the
% world is volatile and that contingencies change over time.

    omega = .75; % By default we here set this to 1 (indicating no forgetting, 
               % but try changing its value to see how it affects model behavior. 
               % Values below 1 indicate greater rates of forgetting.
               
% Beta: Expected precision of expected free energy (G) over policies (a 
% positive value, with higher values indicating lower expected precision).
% Lower values increase the influence of habits (E) and otherwise make
% policy selection less deteriministic. For our example simulations we will
% simply set this to its default value of 1:

     beta = 1; % By default this is set to 1, but try increasing its value 
               % to lower precision and see how it affects model behavior

% Alpha: An 'inverse temperature' or 'action precision' parameter that 
% controls how much randomness there is when selecting actions (e.g., how 
% often the agent might choose not to take the hint, even if the model 
% assigned the highest probability to that action. This is a positive 
% number, where higher values indicate less randomness. Here we set this to 
% a high value:

    alpha = 1;  % Any positive number. 1 is very low, 32 is fairly high; 
                 % an extremely high value can be used to specify
                 % deterministic action (e.g., 512)

% ERP: This parameter controls the degree of belief resetting at each 
% time point in a trial when simulating neural responses. A value of 1
% indicates no resetting, in which priors smoothly carry over. Higher
% values indicate degree of loss in prior confidence at each time step.

    erp = 1; % By default we here set this to 1, but try increasing its value  
             % to see how it affects simulated neural (and behavioral) responses
                          
% tau: Time constant for evidence accumulation. This parameter controls the
% magnitude of updates at each iteration of gradient descent. Larger values 
% of tau will lead to smaller updates and slower convergence time, 
% but will also promote greater stability in posterior beliefs. 

    tau = 12; % Here we set this to 12 to simulate smooth physiological responses,   
              % but try adjusting its value to see how it affects simulated
              % neural (and behavioral) responses



%%
mdp.T = T;                    % Number of time steps
mdp.V = V;                    % allowable (deep) policies

%mdp.U = U;                   % We could have instead used shallow 
                                  % policies (specifying U instead of V).

mdp.A = A;                    % state-outcome mapping
mdp.B = B;                    % transition probabilities
mdp.C = C;                    % preferred states
mdp.D = D;                    % priors over initial states

mdp.d = d;                    % enable learning priors over initial states
mdp.a = a;

mdp.eta = eta;                % learning rate
mdp.omega = omega;            % forgetting rate
mdp.alpha = alpha;            % action precision
mdp.beta = beta;              % expected precision of expected free energy over policies
mdp.erp = erp;                % degree of belief resetting at each timestep
mdp.tau = tau;                % time constant for evidence accumulation

label.factor{1}   = 'contexts';   label.name{1}    = {'Green-better','Yellow-better'};
label.factor{2}   = 'choice states';     label.name{2}    = {'start','choose Green','choose yellow'};
label.modality{1} = 'lose/win';  label.outcome{1} = {'null','lose','win'};
label.modality{2} = 'observed action';  label.outcome{2} = {'start','choose Green','choose yellow'};
label.action{2} = {'start','Green','Yellow'};
mdp.label = label;

mdp = spm_MDP_check(mdp);
 
%spm_MDP_factor_graph(mdp)
%% Single trial
MDP = spm_MDP_VB_X_tutorial(mdp);

spm_figure('GetWin', "Trial Figure")
spm_MDP_VB_trial(MDP)

spm_figure('GetWin', 'Trial Bio')
spm_MDP_VB_LFP(MDP)

%% Multi trial

N = 64;

MDP = mdp;

[MDP(1:N)] = deal(MDP);


MDP = spm_MDP_VB_X(MDP);


spm_figure('GetWin','Figure 2'); clf
spm_MDP_VB_game_tutorial_sigurd_5_time_steps(MDP)

%% Multi trial reversal learning
N = 64;

MDP = mdp;

[MDP(1:N)] = deal(MDP);

for i = 1:N/4
    MDP(i).D{1}   = [0 1]'; % Start in the 'Green-better' context for 
                            % early trials
end

for i = N/2+1:N
    MDP(i).D{1}   = [1 0]'; % Switch to 'Yellow-better' context for 
                            % the remainder of the trials
end
 

MDP = spm_MDP_VB_X(MDP);


spm_figure('GetWin','Figure 2'); clf
spm_MDP_VB_game_tutorial_sigurd_5_time_steps(MDP)
%%
spm_MDP_VB_trial(MDP(32))
%MDP_fit = spm_MDP_VB(MDP);
%spm_MDP_DEM()
%% Simulated data used for parameter recovery
%spm_figure("GetWin","Multiple trials plot")
%spm_MDP_VB_game(MDP);

% Invert
mdp.la_true = la;
mdp.rs_true = rs;



DCM.MDP = mdp;

DCM.field = {'eta','beta'};

DCM.U = {MDP.o};
DCM.Y = {MDP.u};

%DCM.X = {MDP.X};

DCM = spm_dcm_mdp_sigurd(DCM);

%spm_DEM(DCM)

subplot(2,2,3)
xticklabels(DCM.field),xlabel('Parameter')
subplot(2,2,4)
xticklabels(DCM.field),xlabel('Parameter')

% Check deviation of prior and posterior means & posterior covariance
%==========================================================================

%--------------------------------------------------------------------------
% re-transform values and compare prior with posterior estimates
%--------------------------------------------------------------------------

field = fieldnames(DCM.M.pE);
for i = 1:length(field)
    if strcmp(field{i},'eta')
        prior(i) = 1/(1+exp(-DCM.M.pE.(field{i})));
        posterior(i) = 1/(1+exp(-DCM.Ep.(field{i}))); 
    elseif strcmp(field{i},'omega')
        prior(i) = 1/(1+exp(-DCM.M.pE.(field{i})));
        posterior(i) = 1/(1+exp(-DCM.Ep.(field{i})));
    else
        prior(i) = exp(DCM.M.pE.(field{i}));
        posterior(i) = exp(DCM.Ep.(field{i}));
    end
end

figure, set(gcf,'color','white')
subplot(2,1,1),hold on
title('Means')
bar(prior,'FaceColor',[.5,.5,.5]),bar(posterior,0.5,'k')
xlim([0,length(prior)+1]),set(gca, 'XTick', 1:length(prior)),set(gca, 'XTickLabel', DCM.field)
legend({'Prior','Posterior'})
hold off
subplot(2,1,2)
imagesc(DCM.Cp),caxis([0 1]),colorbar
title('(Co-)variance')
set(gca, 'XTick', 1:length(prior)),set(gca, 'XTickLabel', DCM.field)
set(gca, 'YTick', 1:length(prior)),set(gca, 'YTickLabel', DCM.field)

%% Empirical States and Observations
format long g

mdp_list = {};
MDP_fit_list = {};

%Data 
datadir = ("C:/Users/Bruger/OneDrive - Aarhus universitet/5th Semester/Bachelor Thesis/Data/raw_srl/");
ID = importdata(append(datadir, "/names.txt")); % Unique ID's

%All behavioral data
data_combined2 = load("C:\Users\Bruger\Documents\raw_srl\Database_20112022\data_combined.mat");
data_combined2 = data_combined2.data_combined;

%create reaction time
data_combined2.t_reaction = data_combined2.t_dec - data_combined2.t_cue;

for n = 1:size(ID,1)
    
    data_subset = data_combined2(strcmp(data_combined2.ID, ID(n)) == 1, :);
    data_subset.t_reaction = data_subset.t_dec - data_subset.t_cue;
    % =========================================================================
    %                                   MDP PART
    % =========================================================================
    MDP_multi = mdp;
    N_trial = size(data_subset,1); % Check number of trials
    
    [ MDP_multi(1:N_trial) ] = deal(MDP_multi);
    
    %--------------------------------------------------------------------------
    % Set the emperical values for states, observations & actions.
    %--------------------------------------------------------------------------

    for i = 1:N_trial
    
    
        % ===== Specify context (first hidden state factor)============
        if data_subset.blockProb(i) == 0.8
            MDP_multi(i).D{1} = [1 0]';
        else
            MDP_multi(i).D{1} = [0 1]';
        end
    
        %============= s matrix =======================================
        % ---------Better Context-----------
        if data_subset.blockProb(i) == 0.8 % Green best 
            MDP_multi(i).s(:,:) = [1 1 1 1;
                                   0 0 0 0];
        else % Yellow best
            MDP_multi(i).s(:,:) = [2 2 2 2;
                                   0 0 0 0];
        end
    
    
        % ------Start green yellow state factor--------
        % If the agent were in chose green state.
        if data_subset.Select_Green(i) == 1 
            if data_subset.t_reaction(i) < 0 
                MDP_multi(i).s(2,:) = [1 1 1 1]; % never chose
            elseif data_subset.t_reaction(i) >= 0 && data_subset.t_reaction(i) < 0.5
                MDP_multi(i).s(2,:) = [1 2 1 1];
            elseif data_subset.t_reaction(i) >= 0.5 && data_subset.t_reaction(i) < 1
                MDP_multi(i).s(2,:) = [1 1 2 1];
            elseif data_subset.t_reaction(i) >= 1 && data_subset.t_reaction(i) < 1.6
                MDP_multi(i).s(2,:) = [1 1 1 2];

            end
    
        % If the agent chose yellow
        elseif data_subset.Select_Green(i) == 0 
            if data_subset.t_reaction(i) < 0 
                MDP_multi(i).s(2,:) = [1 1 1 1]; % never chose
            elseif data_subset.t_reaction(i) >= 0 && data_subset.t_reaction(i) < 0.5
                MDP_multi(i).s(2,:) = [1 3 1 1];
            elseif data_subset.t_reaction(i) >= 0.5 && data_subset.t_reaction(i) < 1
                MDP_multi(i).s(2,:) = [1 1 3 1];
            elseif data_subset.t_reaction(i) >= 1 && data_subset.t_reaction(i) < 1.6
                MDP_multi(i).s(2,:) = [1 1 1 3];
            end
    
        end
          
    
     %=========== o matrix =================
    
       %---------Win/loss----------
    
        if data_subset.correctness(i) == 1 %Win
            if data_subset.t_reaction(i) < 0 
                MDP_multi(i).o = [1 1 1 1]; % never chose
            elseif data_subset.t_reaction(i) >= 0 && data_subset.t_reaction(i) < 0.5
                MDP_multi(i).o = [1 3 1 1];
            elseif data_subset.t_reaction(i) >= 0.5 && data_subset.t_reaction(i) < 1
                MDP_multi(i).o = [1 1 3 1];
            elseif data_subset.t_reaction(i) >= 1 && data_subset.t_reaction(i) < 1.6
                MDP_multi(i).o = [1 1 1 3];
            end
       
        elseif data_subset.correctness(i) == 0
            if data_subset.t_reaction(i) < 0 
                MDP_multi(i).o = [1 1 1 1]; % never chose
            elseif data_subset.t_reaction(i) >= 0 && data_subset.t_reaction(i) < 0.5
                MDP_multi(i).o = [1 2 1 1];
            elseif data_subset.t_reaction(i) >= 0.5 && data_subset.t_reaction(i) < 1
                MDP_multi(i).o = [1 1 2 1];
            elseif data_subset.t_reaction(i) >= 1 && data_subset.t_reaction(i) < 1.6
                MDP_multi(i).o = [1 1 1 2];
            end
        end
          %-------Observing Actions-----
            MDP_multi(i).o(2,:) = MDP_multi(i).s(2,:); 
    
    %=========== u matrix ================
        % If the agente chose green
        if data_subset.Select_Green(i) == 1 
            if data_subset.t_reaction(i) < 0 
                MDP_multi(i).u = [1 1 1;
                                  1 1 1]; % never chose
    
            elseif data_subset.t_reaction(i) >= 0 && data_subset.t_reaction(i) < 0.5
                MDP_multi(i).u = [1 1 1;
                                  2 1 1];
            elseif data_subset.t_reaction(i) >= 0.5 && data_subset.t_reaction(i) < 1
                MDP_multi(i).u = [1 1 1;
                                  1 2 1];
            elseif data_subset.t_reaction(i) >= 1 && data_subset.t_reaction(i) < 1.6
                MDP_multi(i).u = [1 1 1;
                                  1 1 2];

            end
    
        % If the agent chose yellow
        elseif data_subset.Select_Green(i) == 0 
            if data_subset.t_reaction(i) < 0 
                MDP_multi(i).u = [1 1 1;
                                  1 1 1]; % never chose
            elseif data_subset.t_reaction(i) >= 0 && data_subset.t_reaction(i) < 0.5
                MDP_multi(i).u = [1 1 1;
                                  3 1 1];
            elseif data_subset.t_reaction(i) >= 0.5 && data_subset.t_reaction(i) < 1
                MDP_multi(i).u = [1 1 1;
                                  1 3 1];
            elseif data_subset.t_reaction(i) >= 1 && data_subset.t_reaction(i) < 1.6
                MDP_multi(i).u = [1 1 1;
                                  1 1 3];
            end
        end
    end
    
    %------------ Excecute MDP --------------------------------
    % Setup mdp structure with actions, states & observations.
    mdp_list = [mdp_list; {MDP_multi} ];
    

    %temp_fit = spm_MDP_VB_X_tutorial(MDP_multi);
    
    % Simulate behavior with actions, states & observations predefined.
    %MDP_fit_list = [MDP_fit_list; {temp_fit} ];
    

    %Display end time
    disp(datetime("now"))
    fprintf("%03d done\n",n) %Counts number of participants done.

    % Reset workspace
    clear MDP_multi temp_fit;
    clear N_trial;
    clear data_subset;


end 

%Display end time
disp(datetime("now"))

dir = "C:\Users\Bruger\OneDrive - Aarhus universitet\5th Semester\Bachelor Thesis\Active Inference Models\Two factor - 5 Time Steps\";
save(append(dir,'mdp_list.mat'), "mdp_list") % Save the mdp structure
save(append(dir,'MDP_fit_list.mat'), "MDP_fit_list") % Save the trial data
%%
format long g
%data_combined2.t_dec - data_combined2.t_cue
max(data_combined2.t_reaction)


%%
mdp_fit = spm_MDP_VB_X(mdp_list{1});
%%
spm_figure('GetWin','Figure 2'); clf
spm_MDP_VB_game_tutorial_sigurd_5_time_steps(mdp_fit)

%% alpha only model

%Data 
datadir = ("C:/Users/Bruger/OneDrive - Aarhus universitet/5th Semester/Bachelor Thesis/Data/raw_srl/");
ID = importdata(append(datadir, "/names.txt")); % Unique ID's

%save Dir
dir = "C:\Users\Bruger\OneDrive - Aarhus universitet\5th Semester\Bachelor Thesis\Active Inference Models\Two factor - 5 Time Steps\";
MDP_list = load("C:\Users\Bruger\OneDrive - Aarhus universitet\5th Semester\Bachelor Thesis\Active Inference Models\Two factor - 5 Time Steps\mdp_list.mat");
MDP_list = MDP_list.mdp_list;


GCM_1_alpha = {};
F_1_alpha_params = [];
for i = 1:size(ID,1)

    mdp.la_true = 1;
    mdp.rs_true = 4;
    
    DCM.MDP = mdp;
    DCM.field = {"alpha"};
    
    DCM.Y = {MDP_list{i}.u};
    DCM.U =  {MDP_list{i}.o};
    DCM.S = {MDP_list{i}.s};
    
    %DCM_fit = spm_dcm_mdp_sigurd(DCM);
    DCM_fit = Estimate_parameters_sigurd_5_time_step(DCM);
    subplot(2,2,3)
    xticklabels(DCM_fit.field),xlabel('Parameter')
    subplot(2,2,4)
    xticklabels(DCM_fit.field),xlabel('Parameter')
    


    field = fieldnames(DCM_fit.M.pE);
    for i = 1:length(field)
        if strcmp(field{i},'eta')
            DCM_fit.prior(i) = 1/(1+exp(-DCM_fit.M.pE.(field{i})));
            DCM_fit.posterior(i) = 1/(1+exp(-DCM_fit.Ep.(field{i})));
        elseif strcmp(field{i},'omega')
            DCM_fit.prior(i) = 1/(1+exp(-DCM_fit.M.pE.(field{i})));
            DCM_fit.posterior(i) = 1/(1+exp(-DCM_fit.Ep.(field{i})));
        else
            DCM_fit.prior(i) = exp(DCM_fit.M.pE.(field{i}));
            DCM_fit.posterior(i) = exp(DCM_fit.Ep.(field{i}));
        end
    end



    F_1_alpha_params = [F_1_alpha_params DCM_fit.F];% Get free energies for each participant's model

    %Save each participants output.
    GCM_1_alpha = [GCM_1_alpha;{DCM_fit}];
    
    disp(datetime("now"))
    fprintf("%03d done\n",i) %Counts number of participants done.dw
    clear DCM
end 
disp(datetime("now"))

save(append(dir,'GCM_alpha_list.mat'), 'GCM_1_alpha') % Save the output from the DCM for each participant. 
save(append(dir,'F_alpha_params.mat'), 'F_1_alpha_params') % Save the output from the DCM for each participant.

%% beta only model

%Data 
datadir = ("C:/Users/Bruger/OneDrive - Aarhus universitet/5th Semester/Bachelor Thesis/Data/raw_srl/");
ID = importdata(append(datadir, "/names.txt")); % Unique ID's

%save Dir
dir = "C:\Users\Bruger\OneDrive - Aarhus universitet\5th Semester\Bachelor Thesis\Active Inference Models\Two factor - 5 Time Steps\";
MDP_list = load("C:\Users\Bruger\OneDrive - Aarhus universitet\5th Semester\Bachelor Thesis\Active Inference Models\Two factor - 5 Time Steps\mdp_list.mat");
MDP_list = MDP_list.mdp_list;


GCM_1b = {};
F_1b_params = [];
for i = 1:size(ID,1)

    mdp.la_true = 1;
    mdp.rs_true = 4;
    
    DCM.MDP = mdp;
    DCM.field = {"beta"};
    
    DCM.Y = {MDP_list{i}.u};
    DCM.U =  {MDP_list{i}.o};
    DCM.S = {MDP_list{i}.s};
    
    %DCM_fit = spm_dcm_mdp_sigurd(DCM);
    DCM_fit = Estimate_parameters_sigurd_5_time_step(DCM);
    subplot(2,2,3)
    xticklabels(DCM_fit.field),xlabel('Parameter')
    subplot(2,2,4)
    xticklabels(DCM_fit.field),xlabel('Parameter')
    


    field = fieldnames(DCM_fit.M.pE);
    for i = 1:length(field)
        if strcmp(field{i},'eta')
            DCM_fit.prior(i) = 1/(1+exp(-DCM_fit.M.pE.(field{i})));
            DCM_fit.posterior(i) = 1/(1+exp(-DCM_fit.Ep.(field{i})));
        elseif strcmp(field{i},'omega')
            DCM_fit.prior(i) = 1/(1+exp(-DCM_fit.M.pE.(field{i})));
            DCM_fit.posterior(i) = 1/(1+exp(-DCM_fit.Ep.(field{i})));
        else
            DCM_fit.prior(i) = exp(DCM_fit.M.pE.(field{i}));
            DCM_fit.posterior(i) = exp(DCM_fit.Ep.(field{i}));
        end
    end



    F_1b_params = [F_1b_params DCM_fit.F];% Get free energies for each participant's model

    %Save each participants output.
    GCM_1b = [GCM_1b;{DCM_fit}];
    
    disp(datetime("now"))
    fprintf("%03d done\n",i) %Counts number of participants done.dw
    clear DCM
end 
disp(datetime("now"))

save(append(dir,'GCM_beta_list.mat'), 'GCM_1b') % Save the output from the DCM for each participant. 
save(append(dir,'F_beta_params.mat'), 'F_1b_params') % Save the output from the DCM for each participant. 

%% Random stuff
GCM_1b = load("C:\Users\Bruger\OneDrive - Aarhus universitet\5th Semester\Bachelor Thesis\Active Inference Models\Two factor - 5 Time Steps\DCM_beta_list.mat");
GCM_1b = GCM_1b.GCM_1b;

temp = [];
for i = 1:size(GCM_1b,1)
  temp =  [temp;GCM_1b{i}.posterior];
end
%bar(GCM_1b{1}.posterior)
x= 1:48;

figure(2);
scatter(x,temp)



%% Beta Eta Model

%Data 
datadir = ("C:/Users/Bruger/OneDrive - Aarhus universitet/5th Semester/Bachelor Thesis/Data/raw_srl/");
ID = importdata(append(datadir, "/names.txt")); % Unique ID's

%save Dir
dir = "C:\Users\Bruger\OneDrive - Aarhus universitet\5th Semester\Bachelor Thesis\Active Inference Models\Two factor - 5 Time Steps\";
MDP_list = load("C:\Users\Bruger\OneDrive - Aarhus universitet\5th Semester\Bachelor Thesis\Active Inference Models\Two factor - 5 Time Steps\mdp_list.mat");
MDP_list = MDP_list.mdp_list;


GCM_2_beta_eta = {};
F_2_beta_eta_params = [];
for i = 1:size(ID,1)
    disp(datetime("now"))
    mdp.la_true = 1;
    mdp.rs_true = 4;
    
    DCM.MDP = mdp;
    DCM.field = {"beta", "eta"};
    
    DCM.Y = {MDP_list{i}.u};
    DCM.U =  {MDP_list{i}.o};
    DCM.S = {MDP_list{i}.s};
    
    %DCM_fit = spm_dcm_mdp_sigurd(DCM);
    DCM_fit = Estimate_parameters_sigurd_5_time_step(DCM);
    subplot(2,2,3)
    xticklabels(DCM_fit.field),xlabel('Parameter')
    subplot(2,2,4)
    xticklabels(DCM_fit.field),xlabel('Parameter')
    


    field = fieldnames(DCM_fit.M.pE);
    for i = 1:length(field)
        if strcmp(field{i},'eta')
            DCM_fit.prior(i) = 1/(1+exp(-DCM_fit.M.pE.(field{i})));
            DCM_fit.posterior(i) = 1/(1+exp(-DCM_fit.Ep.(field{i})));
        elseif strcmp(field{i},'omega')
            DCM_fit.prior(i) = 1/(1+exp(-DCM_fit.M.pE.(field{i})));
            DCM_fit.posterior(i) = 1/(1+exp(-DCM_fit.Ep.(field{i})));
        else
            DCM_fit.prior(i) = exp(DCM_fit.M.pE.(field{i}));
            DCM_fit.posterior(i) = exp(DCM_fit.Ep.(field{i}));
        end
    end



    F_2_beta_eta_params = [F_2_beta_eta_params DCM_fit.F];% Get free energies for each participant's model

    %Save each participants output.
    GCM_2_beta_eta = [GCM_2_beta_eta;{DCM_fit}];
    
    disp(datetime("now"))
    fprintf("%03d done\n",i) %Counts number of participants done.dw

    clear DCM
end 
disp(datetime("now"))

save(append(dir,'GCM_2_beta_eta_list.mat'), 'GCM_2_beta_eta') % Save the output from the DCM for each participant. 
save(append(dir,'F_2_beta_eta_params.mat'), 'F_2_beta_eta_params') % Save the output from the DCM for each participant. 
%% Alpha Beta Model

%Data 
datadir = ("C:/Users/Bruger/OneDrive - Aarhus universitet/5th Semester/Bachelor Thesis/Data/raw_srl/");
ID = importdata(append(datadir, "/names.txt")); % Unique ID's

%save Dir
dir = "C:\Users\Bruger\OneDrive - Aarhus universitet\5th Semester\Bachelor Thesis\Active Inference Models\Two factor - 5 Time Steps\";
MDP_list = load("C:\Users\Bruger\OneDrive - Aarhus universitet\5th Semester\Bachelor Thesis\Active Inference Models\Two factor - 5 Time Steps\mdp_list.mat");
MDP_list = MDP_list.mdp_list;


%GCM_2_alpha_beta = {};
%F_2_alpha_beta_params = [];
for i = 41:size(ID,1)
    disp(datetime("now"))
    mdp.la_true = 1;
    mdp.rs_true = 4;
    
    DCM.MDP = mdp;
    DCM.field = {"alpha", "beta"};
    
    DCM.Y = {MDP_list{i}.u};
    DCM.U =  {MDP_list{i}.o};
    DCM.S = {MDP_list{i}.s};
    
    %DCM_fit = spm_dcm_mdp_sigurd(DCM);
    DCM_fit = Estimate_parameters_sigurd_5_time_step(DCM);
    subplot(2,2,3)
    xticklabels(DCM_fit.field),xlabel('Parameter')
    subplot(2,2,4)
    xticklabels(DCM_fit.field),xlabel('Parameter')
    


    field = fieldnames(DCM_fit.M.pE);
    for i = 1:length(field)
        if strcmp(field{i},'eta')
            DCM_fit.prior(i) = 1/(1+exp(-DCM_fit.M.pE.(field{i})));
            DCM_fit.posterior(i) = 1/(1+exp(-DCM_fit.Ep.(field{i})));
        elseif strcmp(field{i},'omega')
            DCM_fit.prior(i) = 1/(1+exp(-DCM_fit.M.pE.(field{i})));
            DCM_fit.posterior(i) = 1/(1+exp(-DCM_fit.Ep.(field{i})));
        else
            DCM_fit.prior(i) = exp(DCM_fit.M.pE.(field{i}));
            DCM_fit.posterior(i) = exp(DCM_fit.Ep.(field{i}));
        end
    end



    F_2_alpha_beta_params = [F_2_alpha_beta_params DCM_fit.F];% Get free energies for each participant's model

    %Save each participants output.
    GCM_2_alpha_beta = [GCM_2_alpha_beta;{DCM_fit}];
    
    disp(datetime("now"))
    fprintf("%03d done\n",i) %Counts number of participants done.dw

    clear DCM
end 
disp(datetime("now"))

save(append(dir,'GCM_2_alpha_beta_list.mat'), 'GCM_2_alpha_beta') % Save the output from the DCM for each participant. 
save(append(dir,'F_2_alpha_beta_params.mat'), 'F_2_alpha_beta_params') % Save the output from the DCM for each participant. 



%% Alpha Beta Eta 

%Data 
datadir = ("C:/Users/Bruger/OneDrive - Aarhus universitet/5th Semester/Bachelor Thesis/Data/raw_srl/");
ID = importdata(append(datadir, "/names.txt")); % Unique ID's

%save Dir
dir = "C:\Users\Bruger\OneDrive - Aarhus universitet\5th Semester\Bachelor Thesis\Active Inference Models\Two factor - 5 Time Steps\";
MDP_list = load("C:\Users\Bruger\OneDrive - Aarhus universitet\5th Semester\Bachelor Thesis\Active Inference Models\Two factor - 5 Time Steps\mdp_list.mat");
MDP_list = MDP_list.mdp_list;


%GCM_2_alpha_beta_eta = {};
%F_2_alpha_beta_eta_params = [];
for i = 48:size(ID,1)
    disp(datetime("now"))
    mdp.la_true = 1;
    mdp.rs_true = 4;
    
    DCM.MDP = mdp;
    DCM.field = {"alpha","beta", "eta"};
    disp(DCM.field)
    DCM.Y = {MDP_list{i}.u};
    DCM.U =  {MDP_list{i}.o};
    DCM.S = {MDP_list{i}.s};
    
    %DCM_fit = spm_dcm_mdp_sigurd(DCM);
    DCM_fit = Estimate_parameters_sigurd_5_time_step(DCM);
    subplot(2,2,3)
    xticklabels(DCM_fit.field),xlabel('Parameter')
    subplot(2,2,4)
    xticklabels(DCM_fit.field),xlabel('Parameter')
    


    field = fieldnames(DCM_fit.M.pE);
    for i = 1:length(field)
        if strcmp(field{i},'eta')
            DCM_fit.prior(i) = 1/(1+exp(-DCM_fit.M.pE.(field{i})));
            DCM_fit.posterior(i) = 1/(1+exp(-DCM_fit.Ep.(field{i})));
        elseif strcmp(field{i},'omega')
            DCM_fit.prior(i) = 1/(1+exp(-DCM_fit.M.pE.(field{i})));
            DCM_fit.posterior(i) = 1/(1+exp(-DCM_fit.Ep.(field{i})));
        else
            DCM_fit.prior(i) = exp(DCM_fit.M.pE.(field{i}));
            DCM_fit.posterior(i) = exp(DCM_fit.Ep.(field{i}));
        end
    end



    F_2_alpha_beta_eta_params = [F_2_alpha_beta_eta_params DCM_fit.F];% Get free energies for each participant's model

    %Save each participants output.
    GCM_2_alpha_beta_eta = [GCM_2_alpha_beta_eta;{DCM_fit}];
    
    disp(datetime("now"))
    fprintf("%03d done\n",i) %Counts number of participants done.dw

    clear DCM
end 
disp(datetime("now"))

save(append(dir,'GCM_3_alpha_beta_eta_list.mat'), 'GCM_2_alpha_beta_eta') % Save the output from the DCM for each participant. 
save(append(dir,'F_3_alpha_beta_eta_params.mat'), 'F_2_alpha_beta_eta_params') % Save the output from the DCM for each participant. 

%% Alpha Beta Omega 

%Data 
datadir = ("C:/Users/Bruger/OneDrive - Aarhus universitet/5th Semester/Bachelor Thesis/Data/raw_srl/");
ID = importdata(append(datadir, "/names.txt")); % Unique ID's

%save Dir
dir = "C:\Users\Bruger\OneDrive - Aarhus universitet\5th Semester\Bachelor Thesis\Active Inference Models\Two factor - 5 Time Steps\";
MDP_list = load("C:\Users\Bruger\OneDrive - Aarhus universitet\5th Semester\Bachelor Thesis\Active Inference Models\Two factor - 5 Time Steps\mdp_list.mat");
MDP_list = MDP_list.mdp_list;

GCM_2_alpha_beta_omega = {};
F_2_alpha_beta_omega_params = [];
for i = 1:size(ID,1)
    disp(datetime("now"))
    mdp.la_true = 1;
    mdp.rs_true = 4;
    
    DCM.MDP = mdp;
    DCM.field = {"alpha","beta","omega"};
    disp(DCM.field)
    DCM.Y = {MDP_list{i}.u};
    DCM.U =  {MDP_list{i}.o};
    DCM.S = {MDP_list{i}.s};
    
    %DCM_fit = spm_dcm_mdp_sigurd(DCM);
    DCM_fit = Estimate_parameters_sigurd_5_time_step(DCM);
    subplot(2,2,3)
    xticklabels(DCM_fit.field),xlabel('Parameter')
    subplot(2,2,4)
    xticklabels(DCM_fit.field),xlabel('Parameter')
    


    field = fieldnames(DCM_fit.M.pE);
    for i = 1:length(field)
        if strcmp(field{i},'eta')
            DCM_fit.prior(i) = 1/(1+exp(-DCM_fit.M.pE.(field{i})));
            DCM_fit.posterior(i) = 1/(1+exp(-DCM_fit.Ep.(field{i})));
        elseif strcmp(field{i},'omega')
            DCM_fit.prior(i) = 1/(1+exp(-DCM_fit.M.pE.(field{i})));
            DCM_fit.posterior(i) = 1/(1+exp(-DCM_fit.Ep.(field{i})));
        else
            DCM_fit.prior(i) = exp(DCM_fit.M.pE.(field{i}));
            DCM_fit.posterior(i) = exp(DCM_fit.Ep.(field{i}));
        end
    end



    F_2_alpha_beta_omega_params = [F_2_alpha_beta_omega_params DCM_fit.F];% Get free energies for each participant's model

    %Save each participants output.
    GCM_2_alpha_beta_omega = [GCM_2_alpha_beta_omega;{DCM_fit}];
    
    disp(datetime("now"))
    fprintf("%03d done\n",i) %Counts number of participants done.dw

    clear DCM
end 
disp(datetime("now"))

save(append(dir,'GCM_3_alpha_beta_omega_list.mat'), 'GCM_2_alpha_beta_omega') % Save the output from the DCM for each participant. 
save(append(dir,'F_3_alpha_beta_omega_params.mat'), 'F_2_alpha_beta_omega_params') % Save the output from the DCM for each participant.


%% Alpha Beta Eta Omega 

%Data 
datadir = ("C:/Users/Bruger/OneDrive - Aarhus universitet/5th Semester/Bachelor Thesis/Data/raw_srl/");
ID = importdata(append(datadir, "/names.txt")); % Unique ID's

%save Dir
dir = "C:\Users\Bruger\OneDrive - Aarhus universitet\5th Semester\Bachelor Thesis\Active Inference Models\Two factor - 5 Time Steps\";
MDP_list = load("C:\Users\Bruger\OneDrive - Aarhus universitet\5th Semester\Bachelor Thesis\Active Inference Models\Two factor - 5 Time Steps\mdp_list.mat");
MDP_list = MDP_list.mdp_list;

GCM_2_alpha_beta_eta_omega_1_37 = {};
F_2_alpha_beta_eta_omega_params_1_37 = [];
for i = 1:36 %size(ID,1)
    disp(datetime("now"))
    mdp.la_true = 1;
    mdp.rs_true = 4;
    
    DCM.MDP = mdp;
    DCM.field = {"alpha","beta","eta","omega"};
    disp(DCM.field)
    DCM.Y = {MDP_list{i}.u};
    DCM.U =  {MDP_list{i}.o};
    DCM.S = {MDP_list{i}.s};
    
    %DCM_fit = spm_dcm_mdp_sigurd(DCM);
    DCM_fit = Estimate_parameters_sigurd_5_time_step(DCM);
    subplot(2,2,3)
    xticklabels(DCM_fit.field),xlabel('Parameter')
    subplot(2,2,4)
    xticklabels(DCM_fit.field),xlabel('Parameter')
    


    field = fieldnames(DCM_fit.M.pE);
    for i = 1:length(field)
        if strcmp(field{i},'eta')
            DCM_fit.prior(i) = 1/(1+exp(-DCM_fit.M.pE.(field{i})));
            DCM_fit.posterior(i) = 1/(1+exp(-DCM_fit.Ep.(field{i})));
        elseif strcmp(field{i},'omega')
            DCM_fit.prior(i) = 1/(1+exp(-DCM_fit.M.pE.(field{i})));
            DCM_fit.posterior(i) = 1/(1+exp(-DCM_fit.Ep.(field{i})));
        else
            DCM_fit.prior(i) = exp(DCM_fit.M.pE.(field{i}));
            DCM_fit.posterior(i) = exp(DCM_fit.Ep.(field{i}));
        end
    end



    F_2_alpha_beta_eta_omega_params_1_37 = [F_2_alpha_beta_eta_omega_params_1_37 DCM_fit.F];% Get free energies for each participant's model

    %Save each participants output.
    GCM_2_alpha_beta_eta_omega_1_37 = [GCM_2_alpha_beta_eta_omega_1_37;{DCM_fit}];
    
    disp(datetime("now"))
    fprintf("%03d done\n",i) %Counts number of participants done.dw

    clear DCM
end 
disp(datetime("now"))

save(append(dir,'GCM_4_alpha_beta_eta_omega_list.mat'), 'GCM_2_alpha_beta_eta_omega') % Save the output from the DCM for each participant. 
save(append(dir,'F_4_alpha_beta_eta_omega_params.mat'), 'F_2_alpha_beta_eta_omega_params') % Save the output from the DCM for each participant.

%%
F_4_alpha_beta_eta_omega_params = [F_2_alpha_beta_eta_omega_params_1_37  F_2_alpha_beta_eta_omega_params ];
GCM_4_alpha_beta_eta_omega_list = [GCM_2_alpha_beta_eta_omega_1_37 ; GCM_2_alpha_beta_eta_omega];

save(append(dir,'GCM_4_alpha_beta_eta_omega_list.mat'), 'GCM_4_alpha_beta_eta_omega_list') % Save the output from the DCM for each participant. 
save(append(dir,'F_4_alpha_beta_eta_omega_params.mat'), 'F_4_alpha_beta_eta_omega_params')

%% MODEL COMPARISON

%--------------------------------------------------------------------------
% Load in data from all models to avoid running them again at a later time
% point
%--------------------------------------------------------------------------
data_path = "C:\Users\Bruger\OneDrive - Aarhus universitet\5th Semester\Bachelor Thesis\Active Inference Models\Two factor - 5 Time Steps\";

%Load in F params files for all participants
F_alpha_params = load(append(data_path, 'F_alpha_params.mat')).F_1_alpha_params;
F_beta_params = load(append(data_path, 'F_beta_params.mat')).F_1b_params;
F_2_beta_eta_params = load(append(data_path, 'F_2_beta_eta_params.mat')).F_2_beta_eta_params;
F_2_alpha_beta_params = load(append(data_path, 'F_2_alpha_beta_params.mat')).F_2_alpha_beta_params;
F_3_alpha_beta_eta_params = load(append(data_path, 'F_3_alpha_beta_eta_params.mat')).F_2_alpha_beta_eta_params;
F_3_alpha_beta_omega_params = load(append(data_path, 'F_3_alpha_beta_omega_params.mat')).F_2_alpha_beta_omega_param_new;
F_4_alpha_beta_eta_omega_params = load(append(data_path, 'F_4_alpha_beta_eta_omega_params.mat')).F_4_alpha_beta_eta_omega_params;

% Load in DCM files for all participants
GCM_alpha_list = load(append(data_path, 'GCM_alpha_list.mat')).GCM_1_alpha; % Alpha
GCM__beta_list = load(append(data_path, 'GCM_beta_list.mat')).GCM_1b;  % Beta
GCM_2_alpha_beta_list = load(append(data_path, 'GCM_2_alpha_beta_list.mat')).GCM_2_alpha_beta; %Alpha Beta
GCM_2_beta_eta_list = load(append(data_path, 'GCM_2_beta_eta_list.mat')).GCM_2_beta_eta;     %Beta eta 
GCM_3_alpha_beta_eta_list = load(append(data_path, 'GCM_3_alpha_beta_eta_list.mat')).GCM_2_alpha_beta_eta; %Alpha beta eta 
GCM_3_alpha_beta_omega_list = load(append(data_path, 'GCM_3_alpha_beta_omega_list.mat')).GCM_3_alpha_beta_omega; % Alpha beta omega
GCM_4_alpha_beta_eta_omega_list = load(append(data_path, 'GCM_4_alpha_beta_eta_omega_list.mat')).GCM_4_alpha_beta_eta_omega_list;
%% Run model comparison
% Model comparison
[alpha,exp_r,xp,pxp,bor]  = spm_BMS([F_alpha_params' F_beta_params' F_2_beta_eta_params' F_2_alpha_beta_params' F_3_alpha_beta_eta_params' F_3_alpha_beta_omega_params' F_4_alpha_beta_eta_omega_params']);

disp(' ');
disp(' ');
disp('Protected exceedance probability (pxp):');
disp(pxp);
disp(' ')
disp('expectation of the posterior p(r|y)');
disp(exp_r);
disp(' ');

% ----Summary of output:-----------
% All of the metrics supports the model with one parameter (alpha, action precision).   

%% Test model recoverability

% The model which only fits alpha has the highest exceedance probability &
% protected exceedance probability. Exceedance probability is the
% probability that model k is more commonly expressed than any other
% model in the model space

% We will therefore now test the parameter recoverability of said model
% using the posterior estimates from our subject level model fitting. 

sim_params = [];
true_params = [];
MDP_best_list = {};
%------------------Setup simulation-------------------------------------
% We will simulate in an exact match of the 2-armed bandit as the
% participants were presented to. All other model parameters than alpha
% will be equal to our initial setup. While subject level estimates of
% alpha will be used in the simulation. 

%Data 
datadir = ("C:/Users/Bruger/OneDrive - Aarhus universitet/5th Semester/Bachelor Thesis/Data/raw_srl/");
ID = importdata(append(datadir, "/names.txt")); % Unique ID's

for n = 1:size(ID,1)
    %Save true parameters
    true_params = [true_params ; GCM_3_alpha_beta_omega_list{n}.posterior];
    
    [MDP_best(1:80)] = deal(mdp);
    %Alpha
    [MDP_best(1:80).alpha] = deal(GCM_3_alpha_beta_omega_list{n}.posterior(1));
    [MDP_best(1:80).alpha_true] = deal(GCM_3_alpha_beta_omega_list{n}.posterior(1)); 
    %Beta 
    [MDP_best(1:80).beta] = deal(GCM_3_alpha_beta_omega_list{n}.posterior(2));
    [MDP_best(1:80).beta_true] = deal(GCM_3_alpha_beta_omega_list{n}.posterior(2)); 
    % Eta 
    [MDP_best(1:80).omega] = deal(GCM_3_alpha_beta_omega_list{n}.posterior(3));
    [MDP_best(1:80).omega_true] = deal(GCM_3_alpha_beta_omega_list{n}.posterior(3)); 


    for i = 1:40
     MDP_best(i).D{1} = [1 0]';
    end
    
    for i = 41:55
        MDP_best(i).D{1} = [0 1]';
    end 
    
    for i = 56:80
        MDP_best(i).D{1} = [1 0]';
    end
    
%     for i = 81:105
%         MDP_best(i).D{1} = [0 1]';
%     end
%     
%     for i = 106:120
%         MDP_best(i).D{1} = [1 0]';
%     end
%     for i = 121:160
%         MDP_best(i).D{1} = [0 1]';
%     end
%     
%     for i = 81:105
%         MDP_best(i).D{1} = [0 1]';
%     end
%     
%     for i = 106:120
%         MDP_best(i).D{1} = [1 0]';
%     end
%     for i = 121:160
%         MDP_best(i).D{1} = [0 1]';
%     end
    % Simulate
    MDP_best_sim = spm_MDP_VB_X(MDP_best);
    
    % Save the simulated data for all subjects for plotting later
    MDP_best_list = [MDP_best_list ; {MDP_best_sim}];

    %-------- plot-------------
    
    %  and to show posterior beliefs and behavior:
    %spm_figure('GetWin','Figure 1'); clf    % display behavior
    %spm_MDP_VB_trial(MDP_best_sim(1));
    
    % Visualize simulated neural responses
    %spm_figure('GetWin','Figure 2'); clf    % display behavior
    %spm_MDP_VB_game_tutorial_sigurd(MDP_best_sim)
    
    %spm_figure('GetWin','Figure 3'); clf    % display behavior
    %spm_MDP_VB_game_tutorial_sigurd(MDP_fit_list{1})
    
    
    % ========================================================================= 
    %              Model fitting to attempt recovery of alpha
    %==========================================================================
    
    mdp.la_true = 1;
    mdp.rs_true = 4;
    
    DCM.MDP = mdp;
    DCM.field = {"alpha","beta", "omega"};
    
    DCM.Y = {MDP_best_sim.u};
    DCM.U =  {MDP_best_sim.o};
    DCM.S = {MDP_best_sim.s};
    
    DCM_fit = Estimate_parameters_sigurd_5_time_step(DCM);
    
    subplot(2,2,3)
    xticklabels(DCM_fit.field),xlabel('Parameter')
    subplot(2,2,4)
    xticklabels(DCM_fit.field),xlabel('Parameter')
    

    
    field = fieldnames(DCM_fit.M.pE);
    for i = 1:length(field)
        if strcmp(field{i},'eta')
            DCM_fit.prior(i) = 1/(1+exp(-DCM_fit.M.pE.(field{i})));
            DCM_fit.posterior(i) = 1/(1+exp(-DCM_fit.Ep.(field{i})));
        elseif strcmp(field{i},'omega')
            DCM_fit.prior(i) = 1/(1+exp(-DCM_fit.M.pE.(field{i})));
            DCM_fit.posterior(i) = 1/(1+exp(-DCM_fit.Ep.(field{i})));
        else
            DCM_fit.prior(i) = exp(DCM_fit.M.pE.(field{i}));
            DCM_fit.posterior(i) = exp(DCM_fit.Ep.(field{i}));
        end
    end
    
    %Save all the simulated parameters
    sim_params = [sim_params ; DCM_fit.posterior];
    clear MDP_best DCM_fit DCM
end

save("C:\Users\Bruger\OneDrive - Aarhus universitet\5th Semester\Bachelor Thesis\Active Inference Models\Two factor - 5 Time Steps\Simulated_param.mat",'sim_params')
%% Posterior estimates
posterior_list = [];
for i = 1:48
  posterior_list = [posterior_list ; GCM_3_alpha_beta_omega_list{i}.posterior];
end

%save('C:\Users\Bruger\Documents\raw_srl\posterior_list_alpha_beta_omega.txt',"posterior_list" ,-ascii, -double)

dlmwrite('C:\Users\Bruger\Documents\raw_srl\posterior_list_alpha_beta_omega.txt', posterior_list , 'delimiter','\t','newline','pc','precision',13);

mean(posterior_list(:,:),1)
std(posterior_list(:,:),1)
%%
% Conditional Covariance
posterior_list_cov = [];
for i = 1:48
  posterior_list_cov = [posterior_list_cov ; GCM_2_alpha_beta_list{i}.Cp(1,2)];
end

mean(posterior_list_cov)
std(posterior_list_cov)

%% Correlation of model recovery
% Get correlations and significance
[Correlations_alpha, Significance_alpha] = corrcoef(true_params(:,1), sim_params(:,1));

[Correlations_beta, Significance_beta] = corrcoef(true_params(:,2), sim_params(:,2));

[Correlations_omega, Significance_omega] = corrcoef(true_params(:,3), sim_params(:,3));
%%
%table 
disp(' ');
disp('3-parameter alpha beta omega model:');
disp(' ');
fprintf('Alpha recoverability: r = %.2g\n',Correlations_alpha(1,2));
fprintf('Correlation significance: p = %.2g\n',Significance_alpha(1,2));

disp(' ');
disp('3-parameter alpha beta omega model:');
disp(' ');
fprintf('Beta recoverability: r = %.2g\n',Correlations_beta(1,2));
fprintf('Correlation significance: p = %.2g\n',Significance_beta(1,2));


disp(' ');
disp('3-parameter alpha beta omega model:');
disp(' ');
fprintf('Omega recoverability: r = %.2g\n',Correlations_omega(1,2));
fprintf('Correlation significance: p = %.2g\n',Significance_omega(1,2));
%%

%scatter plot
figure
scatter(true_params(:,1),sim_params(:,1),'filled')
lsline
title('Recoverability: Alpha (two-parameter model)')
xlabel('True (Generative) Alpha') 
ylabel('Estimated Alpha')
text(2.5, 8, ['r = ' num2str(Correlations_alpha(1,2))])
text(2.5, 7, ['p = ' num2str(Significance_alpha(1,2))])
%% Get predicted neural responses for winning model using fitted model parameters and empirical states, observation & actions.

datadir = ("C:/Users/Bruger/OneDrive - Aarhus universitet/5th Semester/Bachelor Thesis/Data/raw_srl/");
ID = importdata(append(datadir, "/names.txt")); % Unique ID's

GCM_3_alpha_beta_omega_list = load('C:\Users\Bruger\OneDrive - Aarhus universitet\5th Semester\Bachelor Thesis\Active Inference Models\Two factor - 5 Time Steps\GCM_3_alpha_beta_omega_list.mat').GCM_3_alpha_beta_omega; % Alpha beta omega
mdp_list = load("C:\Users\Bruger\OneDrive - Aarhus universitet\5th Semester\Bachelor Thesis\Active Inference Models\Two factor - 5 Time Steps\mdp_list.mat").mdp_list;


MDP_fit_list = {};
for n = 1:size(ID,1)
    fprintf("%03d done\n",n) %Counts number of participants done
    
    for t = 1:160
        mdp_list{n}(t).alpha = GCM_3_alpha_beta_omega_list{n}.posterior(1);
        mdp_list{n}(t).beta = GCM_3_alpha_beta_omega_list{n}.posterior(2);
        mdp_list{n}(t).omega = GCM_3_alpha_beta_omega_list{n}.posterior(3);
    end
    
    temp_fit = spm_MDP_VB_X(mdp_list{n});
    
    % Simulate behavior with actions, states & observations predefined.
    MDP_fit_list = [MDP_fit_list; {temp_fit}];
    
end

%Save data
dir = "C:\Users\Bruger\OneDrive - Aarhus universitet\5th Semester\Bachelor Thesis\Active Inference Models\Two factor - 5 Time Steps\";
save(append(dir,'MDP_fit_list.mat'), "MDP_fit_list") % Save the trial data