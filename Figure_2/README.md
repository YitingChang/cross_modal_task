## Figure 2. Single-unit activity before stimulus delivery reflects task rules.
### a.	Example units
  -- Input:  YT071_190721_se.mat (preprocessed data) unit 19, EF0151_190506_se.mat (preprocessed data) unit 108\
  -- Code: [CM_RasterPlot_prestim.m](CM_RasterPlot_prestim.m)
### b.	Heatmap of normalized baseline activity
  -- Input: data_array_100msBin.mat (processed data)\
  -- Code: [CM_zScore_prestim.m](CM_zScore_prestim.m)\
  -- Normalization: (x-mean(baseline))/std(baseline); baseline (-1-0 s from stimulus onset)
### c.	Distribution of respond-to-touch and respond-to-light block selectivity before stimulus onset
### d.	Distribution of tactile and visual stimulus selectivity before stimulus onset
### e.	Distribution of tactile and visual stimulus selectivity after stimulus onset
### f.	Relationship between block-type selectivity before stimulus onset and tHit-tCR selectivity after stimulus onset
  -- Figure 2c-f\
  -- Input: data_array_100msBin.mat (processed data)\
  -- Output: Unit_AUC_BP_100msBin_1000nBoot_-1to2p5.mat (processed data)\
  --         Unit_AUC_SP_100msBin_1000nBoot_-p3top3.mat (processed data)\
  -- Code: [CM_ROC_block_stimulus.m](CM_ROC_block_stimulus.m)
