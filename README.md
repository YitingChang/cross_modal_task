# cross_modal_task
Yi-Ting Chang, Eric A. Finkel, Duo Xu, Daniel H. O’Connor\
Johns Hopkins University\
\
Flexible responses to sensory stimuli based on changing rules are critical for adapting to a dynamic environment. However, it remains unclear how the brain encodes rule information and uses this information to guide behavioral responses to sensory stimuli. Here, we made single-unit recordings while head-fixed mice performed a cross-modal sensory selection task in which they switched between two rules in different blocks of trials: licking in response to tactile stimuli applied to a whisker while rejecting visual stimuli, or licking to visual stimuli while rejecting the tactile stimuli. Along a cortical sensorimotor processing stream including the primary (S1) and secondary (S2) somatosensory areas, and the medial (MM) and anterolateral (ALM) motor areas, the single-trial activity of individual neurons distinguished between the two rules both prior to and in response to the tactile stimulus. Variable rule-dependent responses to identical stimuli could in principle occur via appropriate configuration of pre-stimulus preparatory states of a neural population, which would shape the subsequent response. We hypothesized that neural populations in S1, S2, MM and ALM would show preparatory activity states that were set in a rule-dependent manner to cause processing of sensory information according to the current rule. This hypothesis was supported for the motor cortical areas by findings that (1) the current task rule could be decoded from pre-stimulus population activity in ALM and MM; (2) neural subspaces containing the population activity differed between the two rules both prior to the stimulus and during the stimulus-evoked response; and (3) optogenetic disruption of pre-stimulus states within ALM and MM impaired task performance. Our findings indicate that flexible selection of an appropriate action in response to a sensory input can occur via configuration of preparatory states in the motor cortex.

**Hightlights**

- Task rules are reflected in preparatory activity in sensory and motor cortices.

- Neural subspaces for processing tactile signals depend on the current task rule.

- Motor cortical activity tracks rule switches and is required for flexible rule-guided behavior.

Preprint: [Rule-based modulation of a sensorimotor transformation across cortical areas](https://doi.org/10.1101/2023.08.21.554194) 

## Data
Our data was collected, processed, and stored in the following formats.     
### Raw data
o	Electrophysiology, behavior, and stimulus signals were recorded via the Intan system (.rhd)\
o	The task structure was recorded via the Bcontrol system (.mat)
### Preprocessed data
o	Neuron data and non-neuron data were extracted and combined in MSessionExplorer, a trial-based system for data organization. \
o	Kilosort 1 was used for spike sorting.\
o	Session information including mouse name, date of birth, brain area was added to MSessionExplorer.\
o	Code: [CM_preprocessing_silicon_probe.m](CM_preprocessing_silicon_probe.m) 
### Processed data
o	A table contains processed data for each session.\
\
o	**firingRate**: Firing rates were organized in a multidimensional array (neuron, block, stimulus, decision, time, trial) with size(N, B, S, D, T, K). Unit: Hz. Options include bin size and smooth window for Gaussian filter.
  -- The analyses for sensory-evoked activity use 10 ms bin with Gaussian smooth (50 ms window)\
  -- The analyses for initial activity use 100 ms bin without Gaussian smooth.\
\
o	**trialNum**: Total trial numbers were organized in a multidimensional array (neuron, block, stimulus, decision, trial) with size (N, B, S, D, K).\
\
o	**earlyTrans**: Boolean array (1: early transition). For each block, the early transition is defined as a window from the first trial of that block to the false alarm (included). It is organized in a multidimensional array (neuron, block, stimulus, decision, trial) with size(N,B,S,D,K).\
\
o	**lateTrans**: Boolean array (1: late transition). For each block, the late transition is defined as a window from the first false alarm to the first hit(not included). It is organized in a multidimensional array (neuron, block, stimulus, decision, trial) with size(N,B,S,D,K).\
\
o	**beforeCue**: Boolean array (1: before a transition cue). It is organized in a multidimensional array (neuron, block, stimulus, decision, trial) with size(N,B,S,D,K).\
\
o	**afterCue**: Boolean array (1: after a transition cue). It is organized in a multidimensional array (neuron, block, stimulus, decision, trial) with size(N,B,S,D,K).\
\
o	**behavPerformance**: Table contains (1) the fractions of trial outcomes in respond-to-touch and respond-to-light blocks, (2) the fractions of correct trials in respond-to-touch blocks, respond-to-light blocks, a session.\
\
o	**firstHit**: The trial number of a first hit after block switch.\  
\
o	Code: [CM_data_array.m](CM_data_array.m)

## Analysis
Scripts used to analyze the data and to generate figures in our preprint can be found in the figure folders.

## Contact
Data will be available in [DANDI](https://www.dandiarchive.org/) soon.

If you have any questions or want to request our datasets, please contact:\
Yi-Ting Chang (ytchang[at]jhmi.edu)



