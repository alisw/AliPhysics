// $Id$

/*! 

\page README_rec Reconstruction
 
The reconstruction is a multistage process, driven by the AliMUONReconstructor class :

- We read the RAW data, convert them (convert them back for simulated data) to 
AliMUONDigit. Performed by AliMUONDigitMaker. 
- We calibrate the digits, by subtracting pedestals and multiplying by gains. Performed
 by AliMUONDigitCalibrator. Calibrated digits might be saved (back) to TreeD.
- We make clusters by associating digits, producing AliMUONRawCluster objects that end up
  in TreeR of MUON.RecPoints#.root file(s). Steered by AliMUONClusterReconstructor. @ref AliMUONVClusterFinder "More..."

Part of the reconstruction, but not really steered by AliMUONReconstructor, is the
 last stage, the tracking, i.e. associating clusters to form tracks. Steered by AliMUONTracker.
 Producing AliMUONTrack objects that end up in TreeT of MUON.Tracks#.root and AliESDMuonTrack objects
that end up in the ESD (Event Summary Data), AliESDs.root file. @ref AliMUONVTrackReconstructor "More..."


\section rec_s1 How to tune muon track reconstruction

Several options and adjustable parameters are available for both Kalman and Original
tracking algorithms (hard coded for the moment in AliMUONVTrackReconstructor.cxx):
- *fgkSigmaToCutForTracking* : quality cut used to select new clusters to be
  attached to the track candidate and to select good tracks.
- *fgkMakeTrackCandidatesFast* : if this flag is set to 'true', the track candidates
  are made assuming linear propagation between stations 4 and 5.
- *fgkTrackAllTracks* : according to the value of this flag, in case that several
  new clusters pass the quality cut, either we consider all the possibilities
  (duplicating tracks) or we attach only the best cluster.
- *fgkRecoverTracks* : if this flag is set to 'true', we try to recover the tracks
  lost during the tracking by removing the worst of the 2 clusters attached in the
  previous station.
- *fgkImproveTracks* : if this flag is set to 'true', we try to improve the quality
  of the tracks at the end of the tracking by removing clusters that do not pass
  new quality cut (within the limit that we must keep at least one cluster per
  the station).
- *fgkSigmaToCutForImprovement* : quality cut used when we try to improve the
  quality of the tracks.

This chapter is defined in the READMEraw.txt file.

*/
