%% Bachelor Thesis Cognitive Science 2022 Aarhus Univeristy.
% Author: Sigurd Fyhn Sørensen

clear all
close all
rng('shuffle') % This sets a random seed for every run

% ===== Experiment description ===========

% The general setup of the experiment is fairly simple. 
% There is a yellow and a green card. One of the cards will always 
% have a 80% probability of giving points while the other card will then
% have a 20% probability of giving points. 
% 
% Which card has the 80%/20% changes through the experiment.
% 
% The two cards are presented side by side, it is ambiguous which side
% which color appears on. 
%
% Participants accumulated points by picking a winning card which in the 
% end of the experiment would be converted into a money bonus on top of 
% the baseline pay. 


%% ==================== MODEL ==============================================
% First thing first we have two context scenarios. 
%   1) Green Machine is best
%   2) Yellow Machine is best.

% The task phase two phases (in this version, could createa preceeding 
% perceptual time step where they think) 
%   1) Choose green.
%   2) Choose yellow card.


T = 2; % (1) Start state/observing stimuli, (2) Select green or yellow.

%==========================================================================
%  Priors about initial states: D and d
% =========================================================================

%--------------------------------------------------------------------------
% Specify prior probabilities about initial states in the generative 
% process (D)
% Note: By default, these will also be the priors for the generative model
%--------------------------------------------------------------------------

% ----First hidden state factor is which machine is best (CONTEXT).----- 
% Capital letter refers to the generative process whereas lowercase
% represent generative model. 

D{1} = [1 1]';  % {'green better','yellow better'}
d{1} = [.25 .25]';  % {'green better','yellow better'}

% -----Second hidden state factor (BEHAVIOR)----------
% 100% certainity that I am gonna start in state1
D{2} = [1 0 0]'; % {'start','choose-green','choose-yellow'}

d{2} = [1 0 0]'; % {'start','choose-green','choose-yellow'}


%==========================================================================
%  State-outcome mappings and beliefs: A and a
%==========================================================================
%--------------------------------------------------------------------------
% Specify the probabilities of outcomes given each state in the generative 
% process (A)
% This includes one matrix per outcome modality
% Note: By default, these will also be the beliefs in the generative model
%--------------------------------------------------------------------------

% Our likelihood matrix has Outcomes on as 
%   Rows: Observations 
%   Columns: Context (which machine is best/first state factor)
%   Third dimension: Corresponds to behavior/second state factor)
%   If we have several more state factor we would add additional dimensions


% ---------Outcome Modality 1----------------
%  First we initialize the basic A{1} for all behaviors

Ns = [length(D{1}) length(D{2})];

for i = 1:Ns(2) 

    A{1}(:,:,i) = [1 1; % Null
                   0 0; % Loss
                   0 0];% Win
end


pWin = 0.8; % 80% Chance if you select the right color max(P(win|color)).
% Choose green 
A{1}(:,:,2) = [0      0;     % Null        
               1-pWin pWin;  % Loss
               pWin 1-pWin]; % Win
% Choose yellow
A{1}(:,:,3) = [0         0  ;  % Null        
               pWin   1-pWin;  % Loss
               1-pWin  pWin ]; % Win

% ----------Outcome modality 2--------------- (observering own behavior)
% Identity mapping between behavior states and observed behaviors.
% To ensure the agent knows that behaviors were carried out as planned.
% Here, each row corresponds to each behavior state.

for i = 1:Ns(2) 

    A{2}(i,:,i) = [1 1];

end


a{1} = A{1}*200;
a{2} = A{2}*200;


%--------------------------------------------------------------------------
% Specify prior beliefs about state-outcome mappings in the generative model 
% (a)
% Note: This will simulate learning state-outcome mappings if specified.
%--------------------------------------------------------------------------
     
a{1}(:,:,2) =  [0  0;  % Null        
               .5 .5;  % Loss
               .5 .5]; % Win
     
     
a{1}(:,:,3) = [0  0;  % Null        
              .5 .5;  % Loss
              .5 .5]; % win


%==========================================================================
% Controlled transitions and transition beliefs : B{:,:,u} and b(:,:,u)
%==========================================================================

%--------------------------------------------------------------------------
% Next, we have to specify the probabilistic transitions between hidden states
% under each action (sometimes called 'control states'). 
%--------------------------------------------------------------------------

%------------ State Factor 1 ---------- (Context)
%  Columns: States at time t
%  Rows: States at time t+1 (also refered to as s' in some litterature)
%  Third Dimension: Is an action.


% Context state never changes within a trial.
B{1}(:,:,1) = [1 0; % Green better
               0 1]; %Yellow better


%--------- State Factor 2------------
%  Columns: States at time t
%  Rows: States at time t+1 (also refered to as s' in some litterature)
%  Third Dimension: Is an action.

% Move to the Start state from any other state
 B{2}(:,:,1) = [1 1 1;  % Start State
                0 0 0;  % Choose Green Card
                0 0 0]; % Choose Yellow Card
           
% Move to the Choose Greenn state from any other state
B{2}(:,:,2) = [0 0 0;  % Start State
               1 1 1;  % Choose Green Card
               0 0 0]; % Choose Yellow Card

% Move to the Choose Yellow state from any other state
B{2}(:,:,3) = [0 0 0;  % Start State
               0 0 0;  % Choose Green Card
               1 1 1]; % Choose Yellow Card



%==========================================================================
% Preferred outcomes: C and c
%==========================================================================
%--------------------------------------------------------------------------
% Next, we have to specify the 'prior preferences', encoded here as log
% probabilities. 
%--------------------------------------------------------------------------

% One matrix per outcome modality. 
%   ROWS: Each row is an observation, and each
%   Columns: Is a time point. 
% Negative values indicate lower preference,
%   positive values indicate a high preference. Stronger preferences promote
%   risky choices and reduced information-seeking.

% We can start by setting a 0 preference for all outcomes:

No = [size(A{1},1) size(A{2},1)]; % number of outcomes in 
                                               % each outcome modality
C{1}      = zeros(No(1),T); % Wins/Losses
C{2}      = zeros(No(2),T); % Observed Behaviors

% Then we can specify a 'loss aversion' magnitude (la) at time points 2 
% and 3, and a 'reward seeking' (or 'risk-seeking') magnitude (rs). Here,
% rs is divided by 2 at the third time point to encode a smaller win ($2
% instead of $4) if taking the hint before choosing a slot machine.

la = 1; % By default we set this to 1, but try changing its value to 
        % see how it affects model behavior

rs = 4; % We set this value at the top of the script. 
          % By default we set it to 4, but try changing its value to 
          % see how it affects model behavior (higher values will promote
          % risk-seeking, as described in the main text)

C{1}(:,:) =    [0  0;  % Null
                0  -la;  % Loss
                0  rs];  % win

C{1}(:,2:T)
%==========================================================================
% Allowable policies: U or V. 
%==========================================================================

%--------------------------------------------------------------------------
% Each policy is a sequence of actions over time that the agent can 
% consider. 
%--------------------------------------------------------------------------

% ========== Shalow Policies ======
Np = 3;
Nf = 2;

U = ones(1,Np, Nf);

U(:,:,1) = [1 1 1]; % Context state is not controllable
U(:,:,2) = [1 2 3]; % All three actions in B{2} are allowed

% ========== DEEP POLICIES  =======

%Np = 3; % Number of policies
%Nf = 2; % Number of state factors

%V         = ones(T-1,Np,Nf);

%V(:,:,1) = [1 1 1;
%            1 1 1]; % Context state is not controllable

%V(:,:,2) = [1 1 1;
%            1 2 3];

%==========================================================================
% Habits: E and e. 
%==========================================================================

%--------------------------------------------------------------------------
% Optional: a columns vector with one entry per policy, indicating the 
% prior probability of choosing that policy (i.e., independent of other 
% beliefs). 
%--------------------------------------------------------------------------


% We will not equip our agent with habits in our example simulations, 
% but this could be specified as a follows if one wanted to include a
% strong habit to choose the 4th policy:

% E = [.1 .1 .1 .6 .1]';

% To incorporate habit learning, where policies become more likely after 
% each time they are chosen, one can also specify concentration parameters
% by specifying e. For example:

% e = [1 1 1 1 1]';


%==========================================================================
% Additional optional parameters. 
%==========================================================================
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

    omega = 1; % By default we here set this to 1 (indicating no forgetting, 
               % but try changing its value to see how it affects model behavior. 
               % Values below 1 indicate greater rates of forgetting.
               
% Beta: Expected precision of expected free energy (G) over policies (a 
% positive value, with higher values indicating lower expected precision).
% Lower values increase the influence of habits (E) and otherwise make
% policy selection less deteriministic. For our example simulations we will
% simply set this to its default value of 1:

     beta = 4; % By default this is set to 1, but try increasing its value 
               % to lower precision and see how it affects model behavior

% Alpha: An 'inverse temperature' or 'action precision' parameter that 
% controls how much randomness there is when selecting actions (e.g., how 
% often the agent might choose not to take the hint, even if the model 
% assigned the highest probability to that action. This is a positive 
% number, where higher values indicate less randomness. Here we set this to 
% a high value:

    alpha = 16;  % Any positive number. 1 is very low, 32 is fairly high; 
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
 
% 2. Define MDP Structure
%==========================================================================
%==========================================================================

mdp.T = T;                    % Number of time steps

%mdp.V = V;                    % allowable (deep) policies
mdp.U = U;                   % We could have instead used shallow 
                                  % policies (specifying U instead of V).

mdp.A = A;                    % state-outcome mapping
mdp.B = B;                    % transition probabilities
mdp.C = C;                    % preferred states
mdp.D = D;                    % priors over initial states

mdp.d = d;                    % enable learning priors over initial states


mdp.eta = eta;                % learning rate
mdp.omega = omega;            % forgetting rate
mdp.alpha = alpha;            % action precision
mdp.beta = beta;              % expected precision of expected free energy over policies
mdp.erp = erp;                % degree of belief resetting at each timestep
mdp.tau = tau;                % time constant for evidence accumulation

% Note, here we are not including habits:

    % mdp.E = E;

% or learning other parameters:
    mdp.a = a;                    
    % mdp.b = b;
    % mdp.c = c;
    % mdp.e = e;         

% or specifying true states or outcomes:

    % mdp.s = s;
    % mdp.o = o;
    
% or specifying other optional parameters (described above):

    % mdp.chi = chi;    % confidence threshold for ceasing evidence
                        % accumulation in lower levels of hierarchical models
    % mdp.zeta = zeta;  % occams window for ceasing to consider implausible
                        % policies
      
% We can add labels to states, outcomes, and actions for subsequent plotting:

label.factor{1}   = 'Contexts';          label.name{1}    = {'Green-better','Yellow-better'};
label.factor{2}   = 'Choice states';     label.name{2}    = {'start','choose green','choose yellow'};
label.modality{1} = 'win/lose';          label.outcome{1} = {'null','lose','win'};
label.modality{2} = 'observed action';   label.outcome{2} = {'start','choose green','choose yellow'};
label.action{2} = {'Green','Yellow'};
mdp.label = label;



%--------------------------------------------------------------------------
% Use a script to check if all matrix-dimensions are correct:
%--------------------------------------------------------------------------
mdp = spm_MDP_check(mdp);

clear a A B C D d Nf No Np Ns T




%% Single trial simulation

MDP = spm_MDP_VB_X_tutorial(mdp);


% We can then use standard plotting routines to visualize simulated neural 
% responses

spm_figure('GetWin','Figure 1'); clf    % display behavior
spm_MDP_VB_LFP(MDP); 

%  and to show posterior beliefs and behavior:

spm_figure('GetWin','Figure 2'); clf    % display behavior

spm_MDP_VB_trial(MDP);

%spm_MDP_plot(MDP)

%% 4. Multi-trial simulations
% Next, we can expand the mdp structure to include multiple trials
N = 64; % number of trials

MDP = mdp;

[MDP(1:N)] = deal(MDP);


%MDP_fit = spm_MDP_VB_XX(MDP);
MDP_fit = spm_MDP_VB_X_tutorial(MDP);

% We can again visualize simulated neural responses
spm_figure('GetWin','Figure 2'); clf
spm_MDP_VB_game(MDP_fit)

spm_figure('GetWin','Figure 3'); clf    % display behavior
spm_MDP_VB_game_tutorial_sigurd(MDP_fit); 
%% Multi trials (reversal learning)
N =64; % number of trials (must be multiple of 8)

MDP_multi2 = mdp;

[MDP_multi2(1:N)] = deal(MDP_multi2);

for i = 1:N
    if i < N/4
        MDP_multi2(i).D{1}   = [0 1]'; % 1/4 Left best
    else 
        MDP_multi2(i).D{1}   = [1 0]'; % 1/2 right best
    end
end
  
MDP_multi2_fit = spm_MDP_VB_X_tutorial(MDP_multi2);

% We can then use standard plotting routines to visualize simulated neural 
% responses


%  and to show posterior beliefs and behavior:
spm_figure('GetWin','Figure 1'); clf    % display behavior
spm_MDP_VB_trial(MDP_multi2_fit(1));

% Visualize simulated neural responses
spm_figure('GetWin','Figure 2'); clf    % display behavior
spm_MDP_VB_game_tutorial_sigurd(MDP_multi2_fit); 
spm_MDP_VB_game_tutorial()
%spm_figure('GetWin','Figure 3'); clf    % display behavior
%spm_MDP_VB_LFP(MDP_multi2_fit(1)); 
%% Empirical States and Observations

mdp_list = {};
MDP_fit_list = {};

datadir = ("C:/Users/Bruger/OneDrive - Aarhus universitet/5th Semester/Bachelor Thesis/Data/raw_srl/");
ID = importdata(append(datadir, "/names.txt"));

%----- LOOP through all ID's -----
for n = 1:size(ID,1)
    path = append(datadir,"TNU_PRSSI_",ID(n),"/behavior/srl/PRSSI_",ID(n),".mat");
    
    
    data = importdata(path);
    
    % Create a green selected? col (0 = no, 1 = yes)
    for i = 1:size(data.alldata,1)
        if data.alldata(i,1) == 0 & data.alldata(i,4) == 3
            data.alldata(i,16) = 1;
        
        elseif data.alldata(i,1) == 1 & data.alldata(i,4) == 2 
            data.alldata(i,16) = 1;
        
        else 
            data.alldata(i,16) = 0;
        end 
    end 
    data.alldata_columns{16} = 'Select_Green';
    
    % =========================================================================
    %                                   MDP PART
    % =========================================================================
    Mdp_multi = mdp;
    
    N_trial = size(data.alldata,1);
    
    [ Mdp_multi(1:N_trial) ] = deal(Mdp_multi);
    
    %--------------------------------------------------------------------------
    % Set the emperical values for states, observations & actions.
    %--------------------------------------------------------------------------
    for i = 1:size(data.alldata,1)
        
     % Specify context (first hidden state factor)
        if data.alldata(i,13) == 0.8 % Green best
            Mdp_multi(i).D{1} = [1 0]';
        else 
            Mdp_multi(i).D{1} = [0 1]';
        end
     %============= s matrix ==============
        % -----------Context---------
        if data.alldata(i,13) == 0.8 % Green best
            Mdp_multi(i).s(:,:) = [1 1;
                                   0 0];
        elseif data.alldata(i,13) == 0.2 %Yellow best
           Mdp_multi(i).s(:,:) = [2 2;
                                  0 0];
        end 
    
        % ---------Action-----------
        if data.alldata(i,16) == 1
            Mdp_multi(i).s(2,:) = [1 2]; % Choose green
        elseif data.alldata(i,16) == 0
            Mdp_multi(i).s(2,:) = [1 3]; % Choose yellow
        end
        
     %=========== o matrix =================
    
       %---------Win/loss----------
        if data.alldata(i,5) == 1   %Win
            Mdp_multi(i).o = [1 3]; 
        else                        %Loss
            Mdp_multi(i).o = [1 2];
        end 
    
      %-------Observing Actions-----
        Mdp_multi(i).o(2,:) = Mdp_multi(i).s(2,:); 
    
    
     %=========== u matrix ================
        if data.alldata(i,16) == 1 %Choose green
           Mdp_multi(i).u = [1;
                             2];
        elseif data.alldata(i,16) == 0 %Choose yellow
           Mdp_multi(i).u = [1;
                             3];
        end 
    end 
    

    %------------ Excecute MDP --------------------------------
    temp_fit = spm_MDP_VB_X_tutorial(Mdp_multi);

    % Setup mdp structure with actions, states & observations.
    mdp_list = [mdp_list; {Mdp_multi} ];
    
    % Simulate behavior with actions, states & observations predefined.
    MDP_fit_list = [MDP_fit_list; {temp_fit} ];
    


    fprintf("%03d done\n",n) %Counts number of participants done.
    
    % Reset workspace
    clear Mdp_multi;
    clear temp_fit;
    clear path;
    clear data;
end

%Display end time
disp(datetime("now"))

save('Active Inference Models/Two Factor - 2 Time Steps/mdp_list.mat', "mdp_list") % Save the mdp structure
save("Active Inference Models/Two Factor - 2 Time Steps/MDP_fit_list.mat", "MDP_fit_list") % Save the trial data

%% beta only model
GCM_1b = {};
F_1b_params = [];
for i =1:size(ID,1)

    mdp.la_true = 1;
    mdp.rs_true = 4;
    
    DCM.MDP = mdp;
    DCM.field = {"beta"};
    
    DCM.Y = {MDP_fit_list{i}.u};
    DCM.U =  {MDP_fit_list{i}.o};
    DCM.S = {MDP_fit_list{i}.s};
    
    DCM_fit = Estimate_parameters_sigurd_experiment(DCM);
    
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

    fprintf("%03d done\n",i) %Counts number of participants done.dw
    clear DCM
end 
disp(datetime("now"))

save('Active Inference Models/Two Factor - 2 Time Steps/DCM_beta_list.mat', 'GCM_1b') % Save the output from the DCM for each participant. 
save('Active Inference Models/Two Factor - 2 Time Steps/F_beta_params.mat', 'F_1b_params') % Save the output from the DCM for each participant. 

%%
GCM_1a = {};
F_1a_params = [];
for i =1:size(ID,1)

    mdp.la_true = 1;
    mdp.rs_true = 4;
    
    DCM.MDP = mdp;
    DCM.field = {"alpha"};
    
    DCM.Y = {MDP_fit_list{i}.u};
    DCM.U =  {MDP_fit_list{i}.o};
    DCM.S = {MDP_fit_list{i}.s};
    
    DCM_fit = Estimate_parameters_sigurd_experiment(DCM);
    
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



    F_1a_params = [F_1a_params DCM_fit.F];% Get free energies for each participant's model

    %Save each participants output.
    GCM_1a = [GCM_1a;{DCM_fit}];

    fprintf("%03d done\n",i) %Counts number of participants done.dw
    clear DCM
end 
disp(datetime("now"))

save('Active Inference Models/Two Factor - 2 Time Steps/DCM_alpha_list.mat', 'GCM_1a') % Save the output from the DCM for each participant. 
save('Active Inference Models/Two Factor - 2 Time Steps/F_alpha_params.mat', 'F_1a_params') % Save the output from the DCM for each participant. 




%% 2 parameter model
GCM_2 = {};
F_2_params = [];
for i =1:size(ID,1)

    mdp.la_true = 1;
    mdp.rs_true = 4;
    
    DCM.MDP = mdp;
    DCM.field = {"beta","alpha"};
    
    DCM.Y = {MDP_fit_list{i}.u};
    DCM.U =  {MDP_fit_list{i}.o};
    DCM.S = {MDP_fit_list{i}.s};
    
    DCM_fit = Estimate_parameters_sigurd_experiment(DCM);
    
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



    F_2_params = [F_2_params DCM_fit.F];% Get free energies for each participant's model

    %Save each participants output.
    GCM_2 = [GCM_2;{DCM_fit}];

    fprintf("%03d done\n",i) %Counts number of participants done.dw
    clear DCM
end 
disp(datetime("now"))


save('Active Inference Models/Two Factor - 2 Time Steps/DCM_2_list.mat', 'GCM_2') % Save the output from the DCM for each participant. 
save('Active Inference Models/Two Factor - 2 Time Steps/F_2_params.mat', 'F_2_params') % Save the output from the DCM for each participant. 

%load('Active Inference Models/Two Factor - 2 Time Steps/DCM_list.mat')

%% Random test
for n= 1:size(GCM,1)
    GCM{n}.M.pE
end

for n = 1:size(mdp_list,1)
    mdp_list{n}.D{1}
end

%% Random Plot
for n = 1:size(GCM,1)
    field = fieldnames(GCM{n}.M.pE);
    for i = 1:length(field)
        if strcmp(field{i},'eta')
            GCM{n}.prior(i) = 1/(1+exp(-GCM{n}.M.pE.(field{i})));
            GCM{n}.posterior(i) = 1/(1+exp(-GCM{n}.Ep.(field{i}))); 
        elseif strcmp(field{i},'omega')
            GCM{n}.prior(i) = 1/(1+exp(-GCM{n}.M.pE.(field{i})));
            GCM{n}.posterior(i) = 1/(1+exp(-GCM{n}.Ep.(field{i})));
        else
            GCM{n}.prior(i) = exp(GCM{n}.M.pE.(field{i}));
            GCM{n}.posterior(i) = exp(GCM{n}.Ep.(field{i}));
        end
    end

    
end 

%temp_post = [GCM{1}.posterior ; GCM{2}.posterior]


%% 3 parameter model

GCM_3 = {};
F_3_params = [];
for i =1:size(ID,1)

    mdp.la_true = 1;
    mdp.rs_true = 4;
    
    DCM.MDP = mdp;
    DCM.field = {"beta","alpha", "eta"};
    
    DCM.Y = {MDP_fit_list{i}.u};
    DCM.U =  {MDP_fit_list{i}.o};
    DCM.S = {MDP_fit_list{i}.s};
    
    DCM_fit = Estimate_parameters_sigurd_experiment(DCM);
    
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



    F_3_params = [F_3_params DCM_fit.F];% Get free energies for each participant's model

    %Save each participants output.
    GCM_3 = [GCM_3;{DCM_fit}];

    fprintf("%03d done\n",i) %Counts number of participants done.dw
    clear DCM
end 
disp(datetime("now"))


save('Active Inference Models/Two Factor - 2 Time Steps/DCM_3_list.mat', 'GCM_3') % Save the output from the DCM for each participant. 
save('Active Inference Models/Two Factor - 2 Time Steps/F_3_params.mat', 'F_3_params') % Save the F (variational free energy) from the DCM for each participant. 


%% MODEL COMPARISON

%--------------------------------------------------------------------------
% Load in data from all models to avoid running them again at a later time
% point
%--------------------------------------------------------------------------
data_path = 'C:\Users\Bruger\OneDrive - Aarhus universitet\5th Semester\Bachelor Thesis\Active Inference Models\Two Factor - 2 Time Steps\';
%Load in F params files for all participants
F_2_params = load(append(data_path, 'F_2_params.mat'));
F_3_params = load(append(data_path, 'F_3_params.mat'));
F_1a_params = load(append(data_path, 'F_alpha_params.mat'));
F_1b_params = load(append(data_path, 'F_beta_params.mat'));

% Load in DCM files for all participants
GCM_2 = load(append(data_path, 'DCM_2_list.mat'));
GCM_3 = load(append(data_path, 'DCM_3_list.mat'));
GCM_1a = load(append(data_path, 'DCM_alpha_list.mat'));
GCM_1b = load(append(data_path, 'DCM_beta_list.mat'));


% Model comparison
[alpha,exp_r,xp,pxp,bor]  = spm_BMS([F_1a_params' F_1b_params' F_2_params' F_3_params']);

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

sim_params_alpha = [];
true_params_alpha = [];
MDP_best_alpha_list = {};
%------------------Setup simulation-------------------------------------
% We will simulate in an exact match of the 2-armed bandit as the
% participants were presented to. All other model parameters than alpha
% will be equal to our initial setup. While subject level estimates of
% alpha will be used in the simulation. 

for n = 1:size(ID,1)
    %Save true parameters
    true_params_alpha = [true_params_alpha ; GCM_1a{n}.posterior];
    
    [MDP_best_alpha(1:120)] = deal(mdp);
    
    [MDP_best_alpha(1:120).alpha] = deal(GCM_1a{n}.posterior);
    
    [MDP_best_alpha(1:120).alpha_true] = deal(GCM_1a{n}.posterior); 
    
    for i = 1:40
     MDP_best_alpha(i).D{1} = [1 0]';
    end
    
    for i = 41:55
        MDP_best_alpha(i).D{1} = [0 1]';
    end 
    
    for i = 56:80
        MDP_best_alpha(i).D{1} = [1 0]';
    end
    
    for i = 81:105
        MDP_best_alpha(i).D{1} = [0 1]';
    end
    
    for i = 106:120
        MDP_best_alpha(i).D{1} = [1 0]';
    end
    for i = 121:160
        MDP_best_alpha(i).D{1} = [0 1]';
    end
    % Simulate
    MDP_best_alpha_sim = spm_MDP_VB_X_tutorial(MDP_best_alpha);
    
    % Save the simulated data for all subjects for plotting later
    MDP_best_alpha_list = [MDP_best_alpha_list ; {MDP_best_alpha_sim}];

    %-------- plot-------------
    
    %  and to show posterior beliefs and behavior:
    %spm_figure('GetWin','Figure 1'); clf    % display behavior
    %spm_MDP_VB_trial(MDP_best_alpha_sim(1));
    
    % Visualize simulated neural responses
    %spm_figure('GetWin','Figure 2'); clf    % display behavior
    %spm_MDP_VB_game_tutorial_sigurd(MDP_best_alpha_sim)
    
    %spm_figure('GetWin','Figure 3'); clf    % display behavior
    %spm_MDP_VB_game_tutorial_sigurd(MDP_fit_list{1})
    
    
    % ========================================================================= 
    %              Model fitting to attempt recovery of alpha
    %==========================================================================
    
    mdp.la_true = 1;
    mdp.rs_true = 4;
    
    DCM.MDP = mdp;
    DCM.field = {"alpha"};
    
    DCM.Y = {MDP_best_alpha_sim.u};
    DCM.U =  {MDP_best_alpha_sim.o};
    DCM.S = {MDP_best_alpha_sim.s};
    
    DCM_fit_alpha = Estimate_parameters_sigurd_experiment(DCM);
    
    subplot(2,2,3)
    xticklabels(DCM_fit.field),xlabel('Parameter')
    subplot(2,2,4)
    xticklabels(DCM_fit.field),xlabel('Parameter')
    

    
    field = fieldnames(DCM_fit_alpha.M.pE);
    for i = 1:length(field)
        if strcmp(field{i},'eta')
            DCM_fit_alpha.prior(i) = 1/(1+exp(-DCM_fit_alpha.M.pE.(field{i})));
            DCM_fit_alpha.posterior(i) = 1/(1+exp(-DCM_fit_alpha.Ep.(field{i})));
        elseif strcmp(field{i},'omega')
            DCM_fit_alpha.prior(i) = 1/(1+exp(-DCM_fit_alpha.M.pE.(field{i})));
            DCM_fit_alpha.posterior(i) = 1/(1+exp(-DCM_fit_alpha.Ep.(field{i})));
        else
            DCM_fit_alpha.prior(i) = exp(DCM_fit_alpha.M.pE.(field{i}));
            DCM_fit_alpha.posterior(i) = exp(DCM_fit_alpha.Ep.(field{i}));
        end
    end
    
    %Save all the simulated parameters
    sim_params_alpha = [sim_params_alpha ; DCM_fit_alpha.posterior];
    
end

%% Correlation of model recovery
% Get correlations and significance
[Correlations_alpha_1, Significance_alpha_1] = corrcoef(true_params_alpha, sim_params_alpha);

%table 
disp(' ');
disp('1-parameter alpha model:');
disp(' ');
fprintf('Alpha recoverability: r = %.2g\n',Correlations_alpha_1(1,2));
fprintf('Correlation significance: p = %.2g\n',Significance_alpha_1(1,2));

%scatter plot
figure
scatter(true_params_alpha,sim_params_alpha,'filled')
lsline
title('Recoverability: Alpha (two-parameter model)')
xlabel('True (Generative) Alpha') 
ylabel('Estimated Alpha')
text(2.5, 8, ['r = ' num2str(Correlations_alpha_1(1,2))])
text(2.5, 7, ['p = ' num2str(Significance_alpha_1(1,2))])
%%

%%
%model1 = temp_post
for i = 1:size(GCM_3,1)
    GCM_3{i}.posterior
end
%%

figure, set(gcf,'color','white')
subplot(2,1,1),hold on
title('Means')
bar(prior,'FaceColor',[.5,.5,.5]),bar(posterior,0.5,'k')
xlim([0,length(prior)+1]),set(gca, 'XTick', 1:length(prior)),set(gca, 'XTickLabel', GCM{1}.field)
legend({'Prior','Posterior'})
hold off
subplot(2,1,2)
imagesc(GCM{1}.Cp),caxis([0 1]),colorbar
title('(Co-)variance')
set(gca, 'XTick', 1:length(prior)),set(gca, 'XTickLabel', GCM{1}.field)
set(gca, 'YTick', 1:length(prior)),set(gca, 'YTickLabel', GCM{1}.field)


%%







field = fieldnames(GCM{1}.M.pE);
for i = 1:length(field)
    if strcmp(field{i},'eta')
        prior(i) = 1/(1+exp(-GCM{1}.M.pE.(field{i})));
        posterior(i) = 1/(1+exp(-GCM{1}.Ep.(field{i}))); 
    elseif strcmp(field{i},'omega')
        prior(i) = 1/(1+exp(-GCM{1}.M.pE.(field{i})));
        posterior(i) = 1/(1+exp(-GCM{1}.Ep.(field{i})));
    else
        prior(i) = exp(GCM{1}.M.pE.(field{i}));
        posterior(i) = exp(GCM{1}.Ep.(field{i}));
    end
end

figure, set(gcf,'color','white')
subplot(2,1,1),hold on
title('Means')
bar(prior,'FaceColor',[.5,.5,.5]),bar(posterior,0.5,'k')
xlim([0,length(prior)+1]),set(gca, 'XTick', 1:length(prior)),set(gca, 'XTickLabel', GCM{1}.field)
legend({'Prior','Posterior'})
hold off
subplot(2,1,2)
imagesc(GCM{1}.Cp),caxis([0 1]),colorbar
title('(Co-)variance')
set(gca, 'XTick', 1:length(prior)),set(gca, 'XTickLabel', GCM{1}.field)
set(gca, 'YTick', 1:length(prior)),set(gca, 'YTickLabel', GCM{1}.field)

%%
i = 6;

figure, set(gcf,'color','white')
subplot(2,1,1),hold on
title('Means')
bar(GCM{i}.prior,'FaceColor',[.5,.5,.5]),bar(GCM{i}.posterior,0.5,'k')
xlim([0,length(GCM{i}.prior)+1]),set(gca, 'XTick', 1:length(GCM{i}.prior)),set(gca, 'XTickLabel', GCM{1}.field)
legend({'Prior','Posterior'})
hold off
subplot(2,1,2)
imagesc(GCM{i}.Cp),caxis([0 1]),colorbar
title('(Co-)variance')
set(gca, 'XTick', 1:length(GCM{i}.prior)),set(gca, 'XTickLabel', GCM{i}.field)
set(gca, 'YTick', 1:length(GCM{i}.prior)),set(gca, 'YTickLabel', GCM{i}.field)

%%
% Check deviation of prior and posterior means & posterior covariance
%==========================================================================

%--------------------------------------------------------------------------
% re-transform values and compare prior with posterior estimates
%--------------------------------------------------------------------------

field = fieldnames(DCM_fit.M.pE);
for i = 1:length(field)
    if strcmp(field{i},'eta')
        prior(i) = 1/(1+exp(-DCM_fit.M.pE.(field{i})));
        posterior(i) = 1/(1+exp(-DCM_fit.Ep.(field{i}))); 
    elseif strcmp(field{i},'omega')
        prior(i) = 1/(1+exp(-DCM_fit.M.pE.(field{i})));
        posterior(i) = 1/(1+exp(-DCM_fit.Ep.(field{i})));
    else
        prior(i) = exp(DCM_fit.M.pE.(field{i}));
        posterior(i) = exp(DCM_fit.Ep.(field{i}));
    end
end

figure, set(gcf,'color','white')
subplot(2,1,1),hold on
title('Means')
bar(prior,'FaceColor',[.5,.5,.5]),bar(posterior,0.5,'k')
xlim([0,length(prior)+1]),set(gca, 'XTick', 1:length(prior)),set(gca, 'XTickLabel', DCM_fit.field)
legend({'Prior','Posterior'})
hold off
subplot(2,1,2)
imagesc(DCM_fit.Cp),caxis([0 1]),colorbar
title('(Co-)variance')
set(gca, 'XTick', 1:length(prior)),set(gca, 'XTickLabel', DCM_fit.field)
set(gca, 'YTick', 1:length(prior)),set(gca, 'YTickLabel', DCM_fit.field)



%% Invert model (FITTING WITH SIMULATED) 
 
N =64; % number of trials (must be multiple of 8)

MDP = mdp;

[MDP(1:N)] = deal(MDP);

for i = 1:N
    if i < N/4
        MDP(i).D{1}  = [1 0]'; % 1/4 Left best
    else 
        MDP(i).D{1}  = [0 1]'; % 1/2 right best
    end
end
  
alpha = 8;
eta = 0.7;

[MDP(1:N).alpha] = deal(alpha);
[MDP(1:N).eta] = deal(eta);

MDP = spm_MDP_VB_X_tutorial(MDP);
%%

mdp.alpha_true = alpha;   % Carries over true la value for use in estimation script
mdp.eta_true = eta;   % Carries over true rs value for use in estimation script
mdp.rs_true = 4;
mdp.la_true = 1;


DCM.MDP   = mdp;                  % MDP model that will be estimated
DCM.field = {'beta', "alpha"};       % parameter (field) names to optimise

% Note: If you wanted to fit other parameters, you can simply add their
% field names, such as:

 % DCM.field = {'alpha','eta'}; % if you wanted to fit learning rate
 
% This requires that those parameters are also included in the possible
% parameters specified in the Estimate_parameters script.

% Next we add the true observations and actions of a (simulated)
% participant

DCM.U     = {MDP.o};              % include the observations made by (real 
                                  % or simulated) participants
                                  
DCM.Y     = {MDP.u};              % include the actions made by (real or 
                                  % simulated) participants
%%
DCM       = Estimate_parameters_sigurd(DCM); % Run the parameter estimation function


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


%%
u = [1 1; % Context, 1) Left Best 2) Right Best
     2 3]; % Behhavior, 1) Start, 2) Hint, 3) Choose-Left, 4) Choose-Right  

% Empirical Observations of Outcome Modalities
% Columns: Time Point (1:T)
% Rows: Outcome Modality (A-tensors)

o = [1 2 1; % 1) No Hint, 2) Left-Hint, 3) Right-Hint
     1 1 3; % 1) Null, 2) Loss, 3) Win
     1 2 3]; % Oberserved, 1) Start, 2) Choose Hint, 3) Choose-Left 4) Choose-Right 

mdp.o = o;
mdp.u = u;

MDP = spm_MDP_VB_X_tutorial(mdp);

