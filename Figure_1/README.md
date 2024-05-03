
## Figure 1. Touch-evoked activity of individual neurons is modulated by task rules.
### a.	Schematic of the cross-modal selection task.
### b.	Example behavioral session.
  -- Input: YT084_200124_se.mat (preprocessed data)\
  -- Code: [CM_behavioral_session.m](CM_behavioral_session.m)
### c.	Task design and trial outcomes.
### d.	Fractions of trial outcomes.
  -- Input: data_array_10msBin.mat (processed data)\
  -- Code: [CM_behavioral_performance.m](CM_behavioral_performance.m)
### e.	Reconstructed locations of silicon probes.
### f.	Example unit.
  -- Input: YT084_200201_se.mat (preprocessed data) unit # 127\
  -- Code: [CM_RasterPlot_response.m](CM_RasterPlot_response.m)
### g.	Normalized population activity
  -- Input: data_array_10msBin.mat (processed data)\
  -- Code: [CM_zScore.m](CM_zScore.m)\
  -- Normalization: (x-mean(baseline))/std(baseline); baseline (-1-0 s from stimulus onset)
### h.	Distribution of the tHit and tCR selectivity.
  -- Input: data_array_10msBin.mat (processed data)\
  -- Output: Unit_AUC_tCorrect_1000nBoot.mat (processed data)\
  -- Code: [CM_ROC_tHit_tCR.m](CM_ROC_tHit_tCR.m)
