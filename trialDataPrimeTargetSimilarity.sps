* Encoding: UTF-8.

PRESERVE.
SET DECIMAL DOT.

GET DATA  /TYPE=TXT
  /FILE=
    "/Users/treber/scripts/neuralAdapatationInMTL/trialDataPrimeTargetSimilarity.csv"
  /ENCODING='UTF8'
  /DELIMITERS=","
  /QUALIFIER='"'
  /ARRANGEMENT=DELIMITED
  /FIRSTCASE=2
  /DATATYPEMIN PERCENTAGE=95.0
  /VARIABLES=
  unitNo AUTO
  zscores AUTO
  condition AUTO
  regionname AUTO
  prevStimName AUTO
  stimName AUTO
  similarityToPrevStim AUTO
  similarityCondition AUTO
  /MAP.
RESTORE.
CACHE.
EXECUTE.
DATASET NAME DataSet1 WINDOW=FRONT.

DATASET ACTIVATE DataSet1.
DATASET DECLARE agg.
AGGREGATE
  /OUTFILE='agg'
  /BREAK=unitNo similarityCondition
  /zscores_mean=MEAN(zscores) 
  /regionname_first=FIRST(regionname) 
  /similarityToPrevStim_mean=MEAN(similarityToPrevStim)
  /N_BREAK=N.

DATASET ACTIVATE agg.


SORT CASES BY unitNo similarityCondition.
CASESTOVARS
  /ID=unitNo
  /INDEX=similarityCondition
  /GROUPBY=VARIABLE.

DATASET ACTIVATE agg.

SORT CASES  BY regionname_first.
SPLIT FILE SEPARATE BY regionname_first.


T-TEST PAIRS=zscores_mean.PrimedHighSimilarity WITH zscores_mean.PrimedLowSimilarity (PAIRED)
  /ES DISPLAY(TRUE) STANDARDIZER(SD)
  /CRITERIA=CI(.9500)
  /MISSING=ANALYSIS.

SPLIT FILE OFF. 

