# Analysis-Scripts for ''Single-neuron mechanisms of neural adaptation in the human temporal lobe''

## Data & dependencies
1. Sorted Spiking and behavioral data are stored in `sessions.mat` (~5GB), which can be downloaded here: 
2. Segmented intracranial EEG traces are stored in `ieeg_linked_mastoids_256Hz.mat` (~5GB) and can be downloaded here: 
3. Some scripts in this repository depend on scripts that can be found here: [https://github.com/mormannlab/common_analysis_scripts](https://github.com/mormannlab/common_analysis_scripts)

## Behavioral Results & Figure 1B
1. adjust paths in `BehaviorAggregate.m` to your situation
2. running `BehaviorAggregate.m` will generate the file `reactiontimes_primed_control_category.mat`.
2. running `BehaviorAnalysesAndPlots.m` will output stats in the MATALB prompt and generate `Figure1B.png`.

## iEEG ERPs
