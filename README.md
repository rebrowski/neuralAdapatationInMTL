# Analysis-Scripts for ''Single-neuron mechanisms of neural adaptation in the human temporal lobe''

## Data & dependencies
1. Sorted Spiking and behavioral data are stored in `sessions.mat` (~5GB), which can be downloaded here: 
2. Segmented intracranial EEG traces are stored in `ieeg_linked_mastoids_256Hz.mat` (~5GB) and can be downloaded here: 
3. Some scripts in this repository depend on scripts that can be found here: [https://github.com/mormannlab/common_analysis_scripts](https://github.com/mormannlab/common_analysis_scripts)

## Paths and startup.m
It is assumed that there is directory to the data, one for the scripts in this repo, and one for the abovementioned `common_analysis_scripts`. To adjust these paths to your situation, edit `startup.m` and run it once at startup of matlab.

## Behavioral Results & Figure 1B
1. running `BehaviorAggregate.m` will generate the file `reactiontimes_primed_control_category.mat`.
2. running `BehaviorAnalysesAndPlots.m` will output stats in the MATALB prompt and generate `Figure1B.png`.

## iEEG ERPs
1. running `ERPsPlotfigure1C.m` will generate the figure `Figure1B.png`
2. running `ERPsLatencyAnalysesSITable1.m` will print output to terminal that was used for Supporting Information Table 1. 
