# cross_modal_task
This repository contains methods to analyze single-neuron activity, neural population activity, and animal behaviors in a cross-modal sensory selection task. The findings are shown in our preprint: [Rule-based modulation of a sensorimotor transformation across cortical areas](https://doi.org/10.1101/2023.08.21.554194).

## Summary
Yi-Ting Chang, Eric A. Finkel, Duo Xu, Daniel H. O’Connor\
Johns Hopkins University

Flexible responses to sensory stimuli based on changing rules are critical for adapting to a dynamic environment. However, it remains unclear how the brain encodes and uses rule information to guide behavior. Here, we made single-unit recordings while head-fixed mice performed a cross-modal sensory selection task where they switched between two rules: licking in response to tactile stimuli while rejecting visual stimuli, or vice versa. Along a cortical sensorimotor processing stream including the primary (S1) and secondary (S2) somatosensory areas, and the medial (MM) and anterolateral (ALM) motor areas, single-neuron activity distinguished between the two rules both prior to and in response to the tactile stimulus. We hypothesized that neural populations in these areas would show rule-dependent preparatory states, which would shape the subsequent sensory processing and behavior. This hypothesis was supported for the motor cortical areas (MM and ALM) by findings that (1) the current task rule could be decoded from pre-stimulus population activity; (2) neural subspaces containing the population activity differed between the two rules; and (3) optogenetic disruption of pre-stimulus states impaired task performance. Our findings indicate that flexible action selection in response to sensory input can occur via configuration of preparatory states in the motor cortex. 

**Hightlights**

- Task rules are reflected in preparatory activity in sensory and motor cortices.

- Neural subspaces for processing tactile signals depend on the current task rule.

- Motor cortical activity tracks rule switches and is required for flexible rule-guided behavior.

## Electrophysiology Data 
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

## Data of Inhibition experiments
### Raw data
o	Behavior, stimulus and laser signals were recorded via the Intan system (.rhd)\
o	The task structure was recorded via the Bcontrol system (.mat)
### Preprocessed data
o	Behavior, stimulus, and laser signals were extracted and combined in MSessionExplorer, a trial-based system for data organization.\ 
o	Session information including mouse name, date of birth, brain area was added to MSessionExplorer.\
o	Code: [photoinhibition_preprocessing.m](photoinhibition_preprocessing.m)
### Processed data
o	A table contains processed data for each session.\
\
o	**lick**: a 4D logical array (b,s,l,o) for lick info in different conditions (block, stimulus, decision, opto-inhibition)\
\
b: respond-to-touch (1) and respond-to-light (2) blocks\
s: tactile (1), visual (2), no (3-5) stimuli\
l: right (1), left (2), no (3) licks\
o: no (1), pre- (2), post- (3) optoinhibition \
For example, data(1,1,1,2) includes all tactile stimulus trials in respond-to-touch blocks and optoinhibition is before stimulus onsets. \
One means right lick is detected whereas zero does not for that trial. \
Note: For trials without stimulus, s=3 for all catch trials, s=4 for short catch trials (0.8 sec), s=5 for long catch trials (2 sec); o=1 for ITI, o=3 for laser trials.  \
Response window: 0~2 second from stimulus onsets.\
\
o	**isPassed**: whether the baseline behavioral performance was passed the criteria (see Method). \
\
o	**Baseline_performance**: hit rates, block correction rates, overall correction rates, laser catch rates\
\
o	**trialNum_tBlock**: a 4D array (b,s,l,o) for trial number info in respond-to-touch blocks.\
\
o	**trialNum_vBlock**: a 4D array (b,s,l,o) for trial number info in respond-to-light blocks.\
\
o	**transition_array**: a 4D logical array (b,s,l,o) for transition info in different conditions.\
\
o	Code: [photoinhibition_processing.m](photoinhibition_processing.m)

## Analysis
Scripts used to analyze the data and to generate figures in our preprint can be found in the figure folders.

## Data availability
The preprocessed data converted to the [Neurodata Without Borders (NWB)](https://www.nwb.org/) format is available at [DANDI archive](https://doi.org/10.48324/dandi.000232/0.240510.2038).\
The processed data is available at [Zenodo](https://zenodo.org/records/11176244). 

## Contact
If you have any questions or want to request our datasets, please contact:\
Yi-Ting Chang (ytchang[at]jhmi.edu)



