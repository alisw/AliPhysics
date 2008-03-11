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


\section rec_s1 How to tune the muon reconstruction

Several options and adjustable parameters allow to tune the entire reconstruction.
These can be changed by adding the following lines in the reconstruction macro (@ref runReconstruction.C) :

<pre>
  AliMUONRecoParam *muonRecoParam = AliMUONRecoParam::GetLow(High)FluxParam();
  muonRecoParam->Use...();
  muonRecoParam->Set...();
  ...
  AliRecoParam::Instance()->RegisterRecoParam(muonRecoParam);
</pre>

Have a look at @ref AliMUONRecoParam for the complete list
of options/parameters with their purpose.


\section rec_s2 How to convert MUON digit/cluster/track into ESD pad/cluster/track and vice versa

The class @ref AliMUONESDInterface is doing the job.
There are 2 ways of using this class:
1) using the static methods to convert the objects one by one (and possibly put them
   into the provided store).
2) loading an entire ESDEvent and using the getters and/or the iterators
   to access the corresponding MUON objects.
Note: You can change (via static method) the type of the store this class is using.

This chapter is defined in the READMErec.txt file.

*/
