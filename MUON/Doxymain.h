/*! \mainpage MUON code documentation 

This is the documentation for the MUON simulation and reconstruction code.
It is a mix of general concepts and code implementation details...

\section ssimulation Simulation

The simulation encompasses the following tasks :

- generation of MC particles (the kinematics of the event ends up in the TreeK
 of Kinematics#.root)                              
  
- tracking particles through the detector using 
the Virtual Monte Carlo, producing AliMUONHit objects, that end up in 
 the TreeH of MUON.Hits#.root file(s). This part is steered by AliMUON and its child
AliMUONv1 classes.

- converting MC hits into AliMUONVDigit, called SDigits, that end up in the TreeS
 of the MUON.SDigits#.root file(s). A S(ummable)Digit is a pad with its associated
charge, but no noise or electronics response function applied. Steered by AliMUONSDigitizerV2 class.

- converting SDigits into Digits, by applying electronics calibrations. Actually, we de-calibrate
 the digits at this stage, by adding a pedestal and dividing by a gain, more or less. Steered
 by AliMUONDigitizerV3 class. Digits end up in TreeD of MUON.Digits#.root file(s). In addition,
 for the trigger, we create AliMUONLocalTrigger, AliMUONRegionalTrigger and AliMUONGlobalTrigger objects 
at this stage, that ends up in TreeD as well.

- convert the Digits into RAW data, in a format that should be exactly the same as real data from the
DAQ. Performed by AliMUONRawWriter.

From there on, the reconstruction can proceed, in the very same way for real or simulated data,
 as long as they are in RAW format.

\section sreconstruction Reconstruction

The reconstruction is a multistage process, driven by the AliMUONReconstructor class :

- we read the RAW data, convert them (convert them back for simulated data) to 
AliMUONDigit. Performed by AliMUONDigitMaker. 
- we calibrate the digits, by subtracting pedestals and multiplying by gains. Performed
 by AliMUONDigitCalibrator. Calibrated digits might be saved (back) to TreeD.
- we make clusters by associating digits, producing AliMUONRawCluster objects that end up
  in TreeR of MUON.RecPoints#.root file(s). Steered by AliMUONClusterReconstructor. @ref AliMUONVClusterFinder "More..."

Part of the reconstruction, but not really steered by AliMUONReconstructor, is the
 last stage, the tracking, i.e. associating clusters to form tracks. Steered by AliMUONTracker.
 Producing AliMUONTrack objects that end up in TreeT of MUON.Tracks#.root and AliESDMuonTrack objects
that end up in the ESD (Event Summary Data), AliESDs.root file. @ref AliMUONVTrackReconstructor "More..."

\section sdataaccess Data access

All the simulation and reconstruction use containers (called stores in MUON jargon) to hold
 the data we're dealing with : hits, (s)digits, trigger, clusters, tracks and trigger tracks.
 All those stores share some commonalities, in particular with respect to how they are read/written from/to
 TTree. @ref AliMUONVStore "More..."

See more in:
- \ref README_main 
- \ref README_raw 
- \ref README_mapping 
- \ref README_calib 
- \ref README_geometry 
- \ref README_trigger 
- \ref README_shuttle 
- \ref README_evaluation 
- \ref README_fast

*/
