## Figure 2 Decoding of task rule in individual neurons
### a.	Example units
  -- Input:  YT071_190721_se.mat (preprocessed data) unit 19\
             EF0151_190506_se.mat (preprocessed data) unit 108\
  -- Code: [CM_RasterPlot_baseline.m]
### b.	Heatmap of normalized baseline activity
  -- Input: AllTrials_100bin_noSmooth\ data_array.mat (processed data)\
  -- Code: [CM_zScore.m]
  -- Normalization: (x-mean(baseline))/std(baseline); baseline (-1-0 s from stimulus onset)
### c.	Distribution of respond-to-touch and respond-to-light block selectivity before stimulus onset
### d.	Distribution of tactile and visual stimulus selectivity before stimulus onset
### e.	Distribution of tactile and visual stimulus selectivity after stimulus onset
### f.	Relationship between block-type selectivity before stimulus onset and tHit-tCR selectivity after stimulus onset
  -- Figure 2c-f
  -- Input: AllTrials_100bin_noSmooth\ data_array.mat (processed data) 
  -- Code: [CM_ROC.m]
