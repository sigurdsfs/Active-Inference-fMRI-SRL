%% Start 
clear all 
close all
rng('shuffle')
%% Model with three hidden state factors
clear all 
close all
rng('shuffle')
T = 2;

D{1} = [1 0]'; %Green better or Yellow better

D{2} = [1 1]'; % Yellow;L/Green;R or Green;L/Yellow;R 
 
D{3} = [1 0 0]'; %start, left , right

d{1} = [.25, .25]';
d{2} = [.25 .25]';
d{3} = [.25 .25 .25]';

Ns = [length(D{1}) length(D{2}) length(D{3})];

for i = 1:Ns(2)
    for j = 1:Ns(3)
        A{1}(:,:,i,j) = [1 1; %null
                         0 0; %Loss
                         0 0]; %win
    end 
end 
pWin = .8;

%When you choose LEFT & Yellow is left
%        YL L    G_best  Y_best  
A{1}(:,:,1,2) = [  0       0 ;   %Null
                  pWin 1-pWin;   %Loss
                 1-pWin  pWin];  %Win

%When you choose LEFT & Green is left.
%        GL L    G_best  Y_best  
A{1}(:,:,2,2) = [  0       0 ;   %Null
                 1-pWin  pWin;   %Loss
                 pWin  1-pWin];  %Win

%When you choose right & Green is right 
%        GR R    G_best  Y_best  
A{1}(:,:,1,3) = [  0       0 ;   %Null
                 1-pWin  pWin;   %Loss
                 pWin  1-pWin];  %Win

%When you choose Right and yellow is right 
%       YR R    G_best  Y_best  
A{1}(:,:,2,3) = [  0        0 ;  %Null
                  pWin  1-pWin;  %Loss
                 1-pWin   pWin]; %Win



for j = 1:Ns(3)
    for i = 1:Ns(2) 

        A{2}(i,:,i,j) = [1 1];

    end
end
a{1} = A{1} * 200;
a{2} = A{2} * 200;

% Hidden state factor 1 (green or yellow better)
B{1} = eye(2,2);

% Hidden state factor 2 (Green left/right
B{2} = eye(2,2);

% Hidden state factor 3 (Start, choose left, choose right)
for i = 1:Ns(3)
    B{3}(i,:,i) = [1 1 1]; 
end


No = [size(A{1},1) size(A{2},1)];

C{1}      = zeros(No(1),T); % Wins/Losses

rs = 4;
la = 4;
C{1}(:,:) =    [0  0;  % Null
                0  -la;  % Loss
                0  rs];  % win


C{2}      = zeros(No(2),T); % Observed Behaviors


Np = 2;
Nf = 3;

U = ones(1,Np, Nf);

U(:,:,1) = [1 1]; % Context state is not controllable
U(:,:,2) = [1 1]; % Green Left or Right is not controlable
U(:,:,3) = [2 3];

eta = 0.5;
omega = 1;
beta = 1;
alpha = 32;
erp = 1;
tau = 32;


mdp.T = T;                    % Number of time steps
%mdp.V = V;                    % allowable (deep) policies

mdp.U = U;                   % We could have instead used shallow 
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


mdp = spm_MDP_check(mdp);


MDP = spm_MDP_VB_X_tutorial(mdp);


%Plot trial
spm_figure('GetWin','trial'); clf
spm_MDP_VB_trial(MDP);


%%


%% Three factor multiple trials
% Next, we can expand the mdp structure to include multiple trials
N = 64; % number of trials

MDP = mdp;

[MDP(1:N)] = deal(MDP);


%MDP_fit = spm_MDP_VB_XX(MDP);
MDP_fit = spm_MDP_VB_X_tutorial(MDP);


spm_figure('GetWin','Figure 3'); clf    % display behavior
spm_MDP_VB_game(MDP_fit); 


DEM.M = mdp;
DEM.U = {MDP_fit.o};
DEM.Y = {MDP_fit.u};

spm_DEM(DEM)
%%
spm_figure('GetWin', 'Single trial, trial = 1'); clf
spm_MDP_VB_trial(MDP_fit(1))


spm_figure('GetWin', 'Single trial, trial = 5'); clf
spm_MDP_VB_trial(MDP_fit(5))

spm_figure('GetWin', 'Single trial, trial = 20'); clf
spm_MDP_VB_trial(MDP_fit(20))


spm_figure('GetWin', 'Single trial, trial = 45'); clf
spm_MDP_VB_trial(MDP_fit(45))
% We can again visualize simulated neural responses



