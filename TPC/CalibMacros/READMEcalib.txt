In this README you can find the list of calibration components used for calibration using reconstructed information.  


Base classes:
AliTPCAnalysisTaskcalib - analysis task owner of the calibration components using reconstructed events and tracks

AliTPCcalibBase  - base class for other components
AliTPCcalibCalib - not very nice name - component to reaply calibration 
                   according the setup in the configuration macro


See components below:
=============================================================================
Standard CPass0,CPass1
=============================================================================
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
AliTPCcalibTime
Input- reconstructed tracks and events
OCDB output:  TPC/Calib/TimeDrift array of graphs
	      +
              in past used as one of input for the TPC/Calib/Corrections 
              
Details:
Kalman fit: drift, drfit velocity gradient t0 shift+ alignment 
            output: part of TPC/CalibTimeDrift
Histogramming:
	a.) TPC/ITS, TPC/Vertex, TPC/TOF, TPC/TRD residual histograms in 4 D
             -cretaion of distortion maps for run by run alignment
            OCDB output: part of TPC/CalibTimeDrift
            additional output: distortion map in meanITSVertex.root
        b.) TPCAside vertex, TPCCside vertex, Global vertex 
            combination histogramming 
             -implemented as a robust method - not used 
              for extraction of calibration because of worse precission, and 
              problems starting strting at high luminosity running.
        c.) TPCdrift laser histogramming 
             -OCDB output: part of TPC/CalibTimeDrift
             -as laser is not 100 % stable, calibration is not used in 
              the reconstruction unless all other methods mentioned
              above failed



xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
AliTPCcalibTimeGain:
Input- reconstructed tracks and events
OCDB output:  TPC/Calib/TimeGain - array of graphs
Details: the calibration of the time dependence of the TPC gain due to pressure and temperature changes etc. In addition signal attenuation calibration 




xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
AliTPCcalibGainMult: 
Input- reconstructed tracks and events
OCDB output:  TPC/Calib/TimeGain - array of graphs
Details: gain/tail as function of multiplicity, gain per chamber
         + experimental code for the  dEdx calibration presented 
           one year ago (disabled) 


=============================================================================
New in CPass0,CPass1
=============================================================================

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
AliTPCcalibAlign         - 
Input - reconstructed events/tracks
Histograms: THnSparse -  residuals  at sector  boundaries,
            THn - 4D cluster/track residuals
Status: code cleaning needed 
Output - in past used as one of input for the TPC/Calib/Corrections 
          automatic procedure missing
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
AliTPCcalibLaser         - 
Input  - reconstructed events/tracks
Actions: drift velocity fit; histogramming of the cluster residuals with 
         respect to survey 
Output - drift velocity (integrated in AliTPCcalibTime)
         in past used as one of input for the TPC/Calib/Corrections 


xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx     AliTPCcalibTracks
Input  - tracks and clusters
Output - TPC/CalibClusterParam
Action: getting resolution and shape parameterizatio and space point bias calibration 
Histogramming : the local residuals between the tracklet and cluster

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx     AliTPCcalibSummary
dumping of the information into the calibration tree used later for calibration visualization

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
Components to be replaced once the high pt trees are created by default:
AliTPCcalibTrigger.cxxAliTPCcalibV0.cxx 

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
Obsolete:
AliTPCcalibTrigger.cxx AliTPCcalibAlignment.cxx 

                 

# calibration visualization
-rw-rw-r-- 1 miranov miranov  40882 Oct  6 23:50 AliTPCcalibSummary.cxx

-rw-rw-r-- 1 miranov miranov   7292 Oct  6 23:50 AliTPCAnalysisTaskcalib.cxx
-rw-rw-r-- 1 miranov miranov  84811 Oct  6 23:50 AliTPCcalibTime.cxx
-rw-rw-r-- 1 miranov miranov  24779 Oct  6 23:50 AliTPCcalibTimeGain.cxx

-rw-rw-r-- 1 miranov miranov 127790 Oct  6 23:50 AliTPCcalibAlign.cxx
-rw-rw-r-- 1 miranov miranov   4571 Oct  6 23:50 AliTPCcalibAlignment.cxx
-rw-rw-r-- 1 miranov miranov  13159 Oct  6 23:50 AliTPCcalibBase.cxx
-rw-rw-r-- 1 miranov miranov  16622 Oct  6 23:50 AliTPCcalibCalib.cxx
-rw-rw-r-- 1 miranov miranov  57212 Oct  6 23:50 AliTPCcalibCosmic.cxx
-rw-rw-r-- 1 miranov miranov  73173 Oct  6 23:50 AliTPCcalibDB.cxx
-rw-rw-r-- 1 miranov miranov 109458 Oct  6 23:50 AliTPCcalibDButil.cxx
-rw-rw-r-- 1 miranov miranov  75162 Oct  6 23:50 AliTPCcalibGainMult.cxx
-rw-rw-r-- 1 miranov miranov 171739 Oct  6 23:50 AliTPCcalibLaser.cxx
-rw-rw-r-- 1 miranov miranov   6754 Oct  6 23:50 AliTPCcalibTracksCuts.cxx
-rw-rw-r-- 1 miranov miranov  50184 Oct  6 23:50 AliTPCcalibTracksGain.cxx
-rw-rw-r-- 1 miranov miranov  11188 Oct  6 23:50 AliTPCcalibTrigger.cxx
-rw-rw-r-- 1 miranov miranov  31505 Oct  6 23:50 AliTPCcalibV0.cxx

