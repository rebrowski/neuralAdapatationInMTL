* Encoding: UTF-8.

PRESERVE.
SET DECIMAL DOT.

GET DATA  /TYPE=TXT
  /FILE="/Users/treber/scripts/neuralAdapatationInMTL/firingRatesByRankCondition.csv" 
  /ENCODING='UTF8'
  /DELIMITERS=","
  /QUALIFIER='"'
  /ARRANGEMENT=DELIMITED
  /FIRSTCASE=2
  /DATATYPEMIN PERCENTAGE=95.0
  /VARIABLES=
  fr AUTO
  zfr AUTO
  rank AUTO
  condition AUTO
  region AUTO
  subject AUTO
  sessid AUTO
  unitId AUTO
  /MAP.
RESTORE.
CACHE.
EXECUTE.
DATASET NAME DataSet1 WINDOW=FRONT.

SORT CASES BY unitId condition subject sessid rank region.
CASESTOVARS
  /ID=unitId region subject sessid
  /INDEX=rank condition
  /GROUPBY=VARIABLE.

DATASET ACTIVATE DataSet1.

GLM fr.1.control fr.1.primed fr.2.control fr.2.primed fr.3.control fr.3.primed 
    fr.4.control fr.4.primed BY region
  /WSFACTOR=rank 4 Polynomial condition 2 Polynomial 
  /MEASURE=fr
  /METHOD=SSTYPE(3)
  /PLOT=PROFILE(rank*condition*region) TYPE=LINE ERRORBAR=SE(2) MEANREFERENCE=NO YAXIS=AUTO
  /PRINT=DESCRIPTIVE ETASQ 
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=rank condition rank*condition
  /DESIGN=region.


GLM zfr.1.control zfr.1.primed zfr.2.control zfr.2.primed zfr.3.control zfr.3.primed zfr.4.control 
    zfr.4.primed BY region
  /WSFACTOR=rank 4 Polynomial condition 2 Polynomial 
  /MEASURE=zfr 
  /METHOD=SSTYPE(3)
  /PLOT=PROFILE(rank*region) TYPE=LINE ERRORBAR=SE(2) MEANREFERENCE=NO YAXIS=AUTO
  /PRINT=DESCRIPTIVE ETASQ 
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=rank condition rank*condition
  /DESIGN=region.
