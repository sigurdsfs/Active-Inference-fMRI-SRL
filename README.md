# Active-Inference-fMRI-SRL
\- "MDP Models" Contains 3 proposed task speicif active inference MDPs.

\- "fMRI" Contains all the scripts used for preprocessing, 1st- and 2nd-level of analysis.

## The Project in brief
In this fMRI project, we intend to use computational active inference models in shape of discrete POMDPs to identify brain regions activated by our task. In this project, 48 healthy male volunteers performed a reward-based associative learning task with winning probabilities of a card (a green or a yellow fractal) changing between 0.8 and 0.2 over 160 trials divided in six blocks. Model parameters will be recovered by using the SPM dynamic expectation maximization (DEM) library to achieve model inversion in relation to the empirically collected states and observation of the participants. Recovered parameter and state trajectories will inform the design matrices for our fMRI analysis. More specifically we will be investigating the brain regions found to be related with learning rate (eta:η), action precision (alpha:α) and expected precision of expected free energy (EFE; beta:β), risk seeking (RS) as well as trial wise phasic dopamine predicted by the neural process theory of active inference.   

## Motivation:
Model based and model free computational modelling combined with fMRI analysis has seen extensive application in the reward learning and decision making literature (O’Doherty et al., 2007). Using an active inference model has the advantage that our cost function is EFE. Since the EFE term contains both an epistemic and pragmatic value, it is ideal for investigating exploit-explore tasks such as the n-armed-bandit task. We hope that extending model-based fMRI analysis to include active-inference models will yield additional insight into decision-making and learning to an extent which traditional reinforcement learning models could not achieve.
