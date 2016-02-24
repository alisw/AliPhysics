// $Id$

/*! 

\page README_rec MUON Reconstruction
 
The reconstruction is a multistage process, driven by the AliMUONTracker and AliMUONReconstructor classes
via the AliReconstruction class, which is divided into three parts:
- the digitization of the electronic response,
- the clustering of the digits to locate the crossing point of the muon with the chamber,
- the tracking to reconstruct the trajectory of the muon in the spectrometer from which we can extract the kinematics.

All the adjustable options and parameters used to tune the different part of the reconstruction are handled by the class AliMUONRecoParam.


\section rec_s1 Digitization

- We read the RAW data, convert them (convert them back for simulated data) to digit (object inheriting from AliMUONVDigit
stored into containers inheriting from AliMUONVDigitStore). This conversion is performed by the class AliMUONDigitMaker.
- We calibrate the digits, via AliMUONDigitCalibrator, by subtracting pedestals and multiplying by a constant gain. All the calibration parameters
(pedestals, HV) are read from the OCDB and stored into AliMUONCalibrationData objects.
- We create the status of the digit (e.g. pedestal higher than maximum or HV switched off), using AliMUONPadStatusMaker.
- We create the status map for each digit, i.e the global status (good/bad) of that digit and of its neighbords, using AliMUONPadStatusMapMaker.
- Calibrated digits might be saved (back) to TreeD in MUON.Digits.root file.


\section rec_s2 Clustering

- We convert the digits having a positive charge into pads (AliMUONPad objects), which also contain information about the digit geometrical
position.
- We loop over pads in the bending and non-bending planes of the DE to form groups of contiguous pads. We then merge the overlapping groups
of pads from both cathodes to build the pre-clusters that are the objects to be clusterized.
- We unfold each pre-cluster in order to extract the number and the position of individual clusters merged in it (complex pre-clusters are
made of a superimposition of signals from muon, from physical background (e.g. hadrons) and from electronic noise).
- We finally determine the MC label: take the one of the simulated track that contribute the most to the total charge of the 2 (bending and the
non bending) pads located below the cluster position. This is possible only if we perform the reconstruction from simulated digits (which contain
the list of MC track contributions). We set it to -1 when reconstructing from raw data or in case of failure.

Several versions of pre-clustering are available, all inheriting from AliMUONVClusterFinder, with different ways to loop over pads to form
pre-clusters:
- AliMUONPreClusterFinder
- AliMUONPreClusterFinderV2
- AliMUONPreClusterFinderV3

Several version of clustering are available, all inheriting from AliMUONVClusterFinder, with different degrees of complexity:
- AliMUONClusterFinderCOG simply compute the Center Of Gravity of the charge distribution in the pre-cluster.
- AliMUONClusterFinderSimpleFit simply fit the charge distribution with a single 2D Mathieson function.
- AliMUONClusterFinderMLEM uses the Maximum Likelihood Expectation Minimization algorithm.
This is a recursive procedure which determines the number and the approximate position of clusters into the pre-cluster that are needed
to reproduce the whole charge distribution. It assumes that the charge distribution of each single cluster follow a 2D Mathieson function.
If the estimated number of clusters is too high (>3), the pre-cluster is split into several groups of 1-2 or 3 clusters selected having
the minimum total coupling to all the other clusters into the pre-cluster. Each group of clusters is then fitted with a sum of 2D Mathieson
functions to extract their exact position.
- AliMUONClusterFinderPeakCOG is a simplified version of the MLEM clusterizer, without splitting and computing the Center Of Gravity of the
local charge distribution to extract the position of every clusters found in the pre-cluster.
- AliMUONClusterFinderPeakFit is another simplified version of the MLEM clusterizer again without splitting. The pre-cluster is fitted with
a sum of 2D Mathieson if it contains less than 3 clusters or we switch to the above COG method.

The cluster recontruction is driven by the class AliMUONSimpleClusterServer, inheriting from AliMUONVClusterServer.
It can be performed either before or during the tracking. In the first case, all the chambers are fully clusterized and the clusters (objects
inheriting from AliMUONVCluster stored into containers inheriting from AliMUONVClusterStore) are saved to TreeR in Muon.RecPoints.root file.
We use the class AliMUONLegacyClusterServer (also inheriting from AliMUONVClusterServer) read back the TreeR and provide clusters to the tracking.
In the second case, we clusterize the chambers only in the region where we are looking for new clusters to be attached to the track candidates.
This makes the clustering faster but the clusters cannot be saved to the TreeR.


\section rec_s3 Tracking

The MUON code provides two different algorithms to reconstruct the muon trajectory. In both cases the general tracking procedure is the same,
the only difference being the way the track parameters are computed from the cluster positions. The "Original" algorithm perform a fit of the
track parameters using the MINUIT package of Root, while the "Kalman" algorithm compute them using analytical formulae. The classes driving
the tracking are AliMUONTrackReconstructor and AliMUONTrackReconstructorK for the "Original" and the "Kalman" algorithms respectively,
both inheriting from AliMUONVTrackReconstructor. The reconstructed muon tracks are objects of the class AliMUONTrack.

The general tracking procedure is as follow:
- Build primary track candidates using clusters on station 4 and 5: Make all combination of clusters between the two chambers of station 5(4).
For each combination compute the local position and impact parameter of the tracklet at vertex and estimate its bending momentum given the averaged
magnetic field inside the dipole and assuming that the track is coming from the vertex. Also compute the corresponding error and covariances of
these parameters. Then select pairs for which the estimated bending momentum and the non-bending impact parameter at vertex are within given limits
taking into account the errors. Extrapolate the primary track candidates to station 4(5), look for at least one compatible cluster to validate them
and recompute the track parameters and covariances.
- Remove the identical track candidates (i.e. the ones sharing exactly the same clusters), and the ones whose bending momentum and non-bending
impact parameter at vertex are out of given limits taking into account the errors.
- Propagate the track to stations 3, 2 then 1. At each station, ask the "ClusterServer" to provide clusters in the region of interest defined in
the reconstruction parameters. Select the one(s) compatible with the track and recompute the track parameters and covariances. Remove the track if
no good cluster is found or if its re-computed bending momentum and non-bending impact parameter at vertex are out of given limits taking into
account the errors.
- Remove the connected tracks (i.e. the ones sharing one cluster or more in stations 3, 4 or 5) keeping the one with the largest number of cluster
or the one with the lowest chi2 in case of equality. Then recompute the track parameters and covariances at each attached cluster (using the
so-called Smoother algorithm in the case of the "Kalman" tracking).
- Find the MC label from the label of each attached cluster (if available): more than 50% of clusters must share the same label, including 1 before
and 1 after the dipole. Set it to -1 when reconstructing real data or in case of failure.
- The reconstructed tracks are finally matched with the trigger tracks (reconstructed from the local response of the trigger) to identify the
muon(s) that made the trigger.

The new clusters to be attached to the track are selected according to their local chi2 (i.e. their transverse position relatively to the track,
normalized by the convolution of the cluster resolution with the resolution of the track extrapolated to the cluster location).
If several compatible clusters are found on the same chamber, the track candidate is duplicated to consider all the possibilities.

The last part of the tracking is the extrapolation of the reconstructed tracks to the vertex of the collision. The vertex position is measured
by the SPD (the Silicon Pixel layers of the ITS Detector). In order to be able to perform any kind of muon analysis, we need to compute the track
parameters assuming the muon has been produced in the initial collision as well as the track parameters in the vertex plane. The first set of
parameters is obtained by correcting for energy loss and multiple Coulomb scattering in the front absorber (we force the track to come from the
exact vertex position (x,y,z) by using the Branson correction), while the second one is obtained by correcting for energy loss only.

The final results of the reconstruction - from which we will perform the physical analyses, compute detector efficiencies and perform calibration
checks - are stored in objects of the class AliESDMuonTrack and saved in AliESD.root file. Three kinds of track can be saved: a tracker track
matched with a trigger track, a tracker track alone and a trigger track alone (unused data members are set to default in the last two cases).
The complete list of MUON data saved into ESD is given in section @ref rec_s5.


\section rec_s4 How to tune the muon reconstruction

Several options and adjustable parameters allow to tune the entire reconstruction. They are stored in the OCDB in the directory MUON/Calib/RecoParam.
However, it is possible to customize the parameters by adding the following lines in the reconstruction macro (runReconstruction.C):
\verbatim
  AliMUONRecoParam *muonRecoParam = AliMUONRecoParam::Get...Param();
  muonRecoParam->Use...();
  muonRecoParam->Set...();
  ...
  MuonRec->SetRecoParam("MUON",muonRecoParam);
\endverbatim

Three sets of default parameters are available:
- <code>GetLowFluxParam()</code>: parameters for p-p collisions
- <code>GetHighFluxParam()</code>: parameters for Pb-Pb collisions
- <code>GetCosmicParam()</code>: parameters for cosmic runs
- <code>GetCalibrationParam()</code>: parameters for cosmic runs

The latter is a dummy set which allows to avoid any reconstruction in case a software trigger event is taken.
Software triggers are sent to trigger electronics during physics run in order to read the scalers: no action from the MUON tracker is required during such events whose reconstruction has to be skipped.

Every option/parameter can be set one by one. Here is the complete list of available setters:
- <code>SetClusteringMode("mode")</code>: set the clustering (pre-clustering) mode: NOCLUSTERING, PRECLUSTER, PRECLUSTERV2, PRECLUSTERV3, COG,
  SIMPLEFIT, SIMPLEFITV3, MLEM:DRAW, MLEM, MLEMV2, MLEMV3.
- <code>SetTrackingMode("mode")</code>: Set the tracking mode: ORIGINAL, KALMAN.
- <code>CombineClusterTrackReco(flag)</code>: switch on/off the combined cluster/track reconstruction
- <code>SaveFullClusterInESD(flag, % of event)</code>: save all cluster info (including pads) in ESD, for the given percentage of events
  (100% by default)
- <code>SelectOnTrackSlope(flag)</code>: switch to select tracks on their slope instead of impact parameter at vertex and/or bending momentum.
- <code>SetMinBendingMomentum(value)</code>: set the minimum acceptable value (GeV/c) of track momentum in bending plane
- <code>SetMaxBendingMomentum(value)</code>: set the maximum acceptable value (GeV/c) of track momentum in bending plane
- <code>SetMaxNonBendingSlope(value)</code>: set the maximum value of the track slope in non bending plane (used when selecting on track slope).
- <code>SetMaxBendingSlope(value)</code>: set the maximum value of the track slope in non bending plane (used when selecting on track slope).
- <code>SetNonBendingVertexDispersion(value)</code>: set the vertex dispersion (cm) in non bending plane (used for the original tracking and to
  select track on their non-bending impact parameter at vertex).
- <code>SetBendingVertexDispersion(value)</code>: set the vertex dispersion (cm) in bending plane (used for the original tracking, to compute the
  error on the estimated bending momentum at the very begining and to select track on their bending impact parameter at vertex (used when B=0)).
- <code>SetMaxNonBendingDistanceToTrack(value)</code>: set the maximum distance to the track to search for compatible cluster(s) in non bending
  direction. This value is convoluted with both the track and the cluster resolutions to define the region of interest.
- <code>SetMaxBendingDistanceToTrack(value)</code>: set the maximum distance to the track to search for compatible cluster(s) in bending direction
  This value is convoluted with both the track and the cluster resolutions to define the region of interest.
- <code>SetSigmaCutForTracking(value)</code>: set the cut in sigma to apply on cluster (local chi2) and track (global chi2) during tracking
- <code>ImproveTracks(flag, sigma cut)</code>: recompute the local chi2 of each cluster with the final track parameters and removed the ones that
  do not pass a new quality cut. The track is removed if we do not end with at least one good cluster per requested station and two clusters in
  station 4 and 5 together whatever they are requested or not.
- <code>ImproveTracks(flag)</code>: same as above using the default quality cut
- <code>SetSigmaCutForTrigger(value)</code>: set the cut in sigma to apply on track during trigger hit pattern search
- <code>SetStripCutForTrigger(value)</code>: set the cut in strips to apply on trigger track during trigger chamber efficiency
- <code>SetMaxStripAreaForTrigger(value)</code>: set the maximum search area in strips to apply on trigger track during trigger chamber efficiency
- <code>SetMaxNormChi2MatchTrigger(value)</code>: set the maximum normalized chi2 for tracker/trigger track matching
- <code>TrackAllTracks(flag)</code>: consider all the clusters passing the sigma cut (duplicate the track) or only the best one
- <code>RecoverTracks(flag)</code>: during the tracking procedure, if no cluster is found in station 1 or 2, we try it again after having removed
  (if possible with respect to the condition to keep at least 1 cluster per requested station) the worst cluster attached in the previous station
  (assuming it was a cluster from background).
- <code>MakeTrackCandidatesFast(flag)</code>: make the primary track candidates formed by cluster on stations 4 and 5 assuming there is no
  magnetic field in that region to speed up the reconstruction.
- <code>MakeMoreTrackCandidates(Bool_t flag)</code>: make the primary track candidate using 1 cluster on station 4 and 1 cluster on station 5
  instead of starting from 2 clusters in the same station.
- <code>ComplementTracks(Bool_t flag)</code>: look for potentially missing cluster to be attached to the track (a track may contain up to 2
  clusters per chamber do to the superimposition of DE, while the tracking procedure is done in such a way that only 1 can be attached).
- <code>RemoveConnectedTracksInSt12(Bool_t flag)</code>: extend the definition of connected tracks to be removed at the end of the tracking
  procedure to the ones sharing one cluster on more in any station, including stations 1 and 2.
- <code>UseSmoother(Bool_t flag)</code>: use or not the smoother to recompute the track parameters at each attached cluster
  (used for Kalman tracking only)
- <code>UseChamber(Int_t iCh, Bool_t flag)</code>: set the chambers to be used (disable the clustering if the chamber is not used).
- <code>RequestStation(Int_t iSt, Bool_t flag)</code>: impose/release the condition "at least 1 cluster per station" for that station.
- <code>BypassSt45(Bool_t st4, Bool_t st5)</code>: make the primary track candidate from the trigger track instead of using stations 4 and/or 5.
- <code>SetHVSt12Limits(float low, float high)</code>: Set Low and High threshold for St12 HV
- <code>SetHVSt345Limits(float low, float high)</code>: Set Low and High threshold for St345 HV
- <code>SetPedMeanLimits(float low, float high)</code>: Set Low and High threshold for pedestal mean
- <code>SetPedSigmaLimits(float low, float high)</code>: Set Low and High threshold for pedestal sigma
- <code>SetPadGoodnessMask(UInt_t mask)</code>: Set the goodness mask (see AliMUONPadStatusMapMaker)
- <code>ChargeSigmaCut(Double_t value)</code>: Number of sigma cut we must apply when cutting on adc-ped
- <code>SetDefaultNonBendingReso(Int_t iCh, Double_t val)</code>: Set the default non bending resolution of chamber iCh
- <code>SetDefaultBendingReso(Int_t iCh, Double_t val)</code>: Set the default bending resolution of chamber iCh
- <code>SetMaxTriggerTracks(Int_t val)</code>: Set the maximum number of trigger tracks above which the tracking is cancelled
- <code>SetMaxTrackCandidates(Int_t val)</code>: Set the maximum number of track candidates above which the tracking abort
- <code>SetManuOccupancyLimits(float low, float high)</code>: Set the limits for the acceptable manu occupancy
- <code>SetBuspatchOccupancyLimits(float low, float high)</code>: Set the limits for the acceptable bp occupancy
- <code>SetDEOccupancyLimits(float low, float high)</code>: Set the limits for the acceptable DE occupancy

We can use the method Print("FULL") to printout all the parameters and options set in the class AliMUONRecoParam.

RecoParams can be put into OCDB using the MakeMUONSingleRecoParam.C or MakeMUONRecoParamArray.C macros.
The first stores only one (default) RecoParam.
The latter allows to store either:
 - LowFlux (default)
 - Calibration

for real data with bunch crossing or
 - Cosmic (default)
 - Calibration

for cosmic runs.

\section rec_s5 ESD content

Three kinds of track can be saved in ESD: a tracker track matched with a trigger track, a tracker track alone and a trigger track alone (unused
data members are set to default values in the last two cases). These tracks are stored in objects of the class AliESDMuonTrack. Two methods can be
used to know the content of an ESD track:
- <code>ContainTrackerData()</code>: Return kTRUE if the track contain tracker data
- <code>ContainTriggerData()</code>: Return kTRUE if the track contain trigger data

The AliESDMuonTrack objects contain:
- Tracker track parameters (x, theta_x, y, theta_y, 1/p_yz) at vertex (x=x_vtx; y=y_vtx)
- Tracker track parameters in the vertex plane
- Tracker track parameters at first cluster
- Tracker track parameter covariances at first cluster
- Tracker track global informations (track ID, chi2, number of clusters, cluster map, MC label if any)
- Array of Ids of associated clusters (clusters are stored in a separate TClonesArray in ESD)
- Trigger track informations (local trigger decision, strip pattern, hit pattern, ...)
- Chi2 of tracker/trigger track matching

The AliESDMuonCluster objects contain:
- Cluster ID providing information about the location of the cluster (chamber ID and DE ID)
- Cluster position (x,y,z)
- Cluster resolution (sigma_x,sigma_y)
- Charge
- Chi2
- MC label if any
- Array of IDs of associated pads for a given fraction of events (pads are stored in a separate TClonesArray in ESD)

The AliESDMuonPad objects contain:
- Digit ID providing information about the location of the digit (DE ID, Manu ID, Manu channel and cathode)
- Raw charge (ADC value)
- Calibrated charge
- One saturation bit and one calibration bit to say whether it is saturated/calibrated or not


\section rec_s6 Conversion between MUON/ESD objects

Every conversion between MUON objects (AliMUONVDigit/AliMUONVCluster/AliMUONTrack) and ESD objects
(AliESDMuonPad/AliESDMuonCluster/AliESDMuonTrack) is done by the class AliMUONESDInterface.

WARNING: some of these conversions require input from outside, namely the magnetic field map, the geometry, the reconstruction parameters
and/or the mapping segmentation. In particular:
- The conversion of ESDPads to MUON digits requires the mapping segmentation.
- The conversion of a MUON track to an ESD track requires the magnetic field and the geometry to extrapolate the track parameters at vertex
and compute the correction of multiple scattering and energy loss in the front absorber.
- While converting an ESD track to a MUON track, the track is refitted by using the cluster position stored in ESD in order to recover the
track parameters at each cluster. This refitting needs both the magnetic field and the reconstruction parameters initially used to
reconstruct the tracks to be correct. The reconstruction parameters can be passed to the interface by using the static method
AliMUONESDTrack::ResetTracker(const AliMUONRecoParam* recoParam, Bool_t info). The refitting can however be disconnected by user (using flag
in the fonction parameters). In that case, none of these external inputs is necessary anymore but only the track parameters at first cluster,
which is then copied directly from the ESD, is meaningful.

There are 2 ways of using this class:

1) Using the static methods to convert the objects one by one (and possibly put them into the provided store):
- Get track parameters at vertex, at DCA, ...:
\verbatim
  ...
  AliESDMuonTrack* esdTrack = esd->GetMuonTrack(iTrack);
  AliMUONTrackParam param;
  AliMUONESDInterface::GetParamAtVertex(*esdTrack, param);
\endverbatim
  
- Convert an AliMUONVDigit to an AliESDMuonPad:
\verbatim
  ...
  AliMUONVDigit *digit = ...;
  AliESDMuonPad esdPad;
  AliMUONESDInterface::MUONToESD(*digit, esdPad);
\endverbatim
  
- Convert an AliMUONLocalTrigger to a ghost AliESDMuonTrack (containing only trigger informations):
\verbatim
  ...
  AliMUONLocalTrigger* locTrg = ...;
  AliMUONTriggerTrack* triggerTrack = ...;
  AliESDMuonTrack esdTrack;
  AliMUONESDInterface::MUONToESD(*locTrg, esdTrack, trackId, triggerTrack);
\endverbatim

- Convert an AliESDMuonTrack to an AliMUONTrack (the parameters at each clusters are recomputed or not according to the flag "refit".
if not, only the parameters at first cluster are relevant):
\verbatim
  ...
  AliESDMuonTrack* esdTrack = esd->GetMuonTrack(iTrack);
  AliMUONTrack track;
  AliMUONESDInterface::ESDToMUON(*esdTrack, track, refit);
\endverbatim
  
- Add an AliESDMuonTrack (converted into AliMUONTrack object) into an AliMUONVTrackStore (same remark as above about the flag "refit"):
\verbatim
  ...
  AliESDMuonTrack* esdTrack = esd->GetMuonTrack(iTrack);
  AliMUONVTrackStore *trackStore = AliMUONESDInteface::NewTrackStore();
  AliMUONTrack* trackInStore = AliMUONESDInterface::Add(*esdTrack, *trackStore, refit);
\endverbatim

2) Loading an entire ESDEvent and using the finders and/or the iterators to access the corresponding MUON objects:
- First load the ESD event:
\verbatim
  AliMUONESDInterface esdInterface;
  esdInterface.LoadEvent(*esd);
\endverbatim
  
- Get the track store:
\verbatim
  AliMUONVTrackStore* trackStore = esdInterface.GetTracks();
\endverbatim
  
- Access the number of digits in a particular cluster:
\verbatim
  Int_t nDigits = esdInterface.GetNDigitsInCluster(clusterId);
\endverbatim
  
- Find a particular digit using its ID:
\verbatim
  AliMUONVDigit *digit = esdInterface.FindDigit(digitId);
\endverbatim

- Find a particular cluster in a given track using their IDs:
\verbatim
  AliMUONVCluster* cluster = esdInterface.FindCluster(trackId, clusterId);
\endverbatim
  
- Iterate over all clusters of a particular track using an iterator:
\verbatim
  TIterator* nextCluster = esdInterface.CreateClusterIterator(trackId);
  while ((cluster = static_cast<AliMUONVCluster*>(nextCluster()))) {...}
\endverbatim

Note: You can change (via static method) the type of the store this class is using:
\verbatim
  AliMUONESDInterface::UseTrackStore("name");
  AliMUONESDInterface::UseClusterStore("name");
  AliMUONESDInterface::UseDigitStore("name");
  AliMUONESDInterface::UseTriggerStore("name");
\endverbatim


\section rec_s7 ESD cluster/track refitting

We can re-clusterize and re-track the clusters/tracks stored into the ESD by using the class AliMUONRefitter. This class gets the MUON objects
to be refitted from an instance of AliMUONESDInterface (see section @ref rec_s6), then uses the reconstruction framework to refit them. The new
reconstruction parameters are still set via the class AliMUONRecoParam (see section @ref rec_s4) and passed to refitter through its constructor.
The old reconstruction parameters, the mapping, the magnetic field and/or the geometry may also be needed to convert the ESD objects to MUON ones
and/or to refit them. The initial data are not changed. Results are stored into new MUON objects. The aim of the refitting is to be able to
study effects of changing the reconstruction parameter, the calibration parameters or the alignment without re-running the entire reconstruction.

To use this class we first have to connect it to the ESD interface containing MUON objects:
\verbatim
  AliMUONRefitter refitter;
  refitter.Connect(&esdInterface);
\endverbatim

We can then:
- Re-clusterize the ESD clusters using the attached ESD pads (several new clusters can be reconstructed per ESD cluster):
\verbatim
  AliMUONVClusterStore* clusterStore = refitter.ReClusterize(iTrack, iCluster);
  AliMUONVClusterStore* clusterStore = refitter.ReClusterize(clusterId);
\endverbatim
  
- Re-fit the ESD tracks using the attached ESD clusters:
\verbatim
  AliMUONTrack* track = refitter.RetrackFromClusters(iTrack);
  AliMUONVTrackStore* trackStore = refitter.ReconstructFromClusters();
\endverbatim
  
- Reconstruct the ESD tracks from ESD pads (i.e. re-clusterize the attached clusters). Consider all the combination of clusters and return only
  the best one:
\verbatim
  AliMUONTrack* track = refitter.RetrackFromDigits(iTrack);
  AliMUONVTrackStore* trackStore = refitter.ReconstructFromDigits();
\endverbatim

The macro MUONRefit.C is an example of using this class. The results are stored in a new AliESDs.root file.


This chapter is defined in the READMErec.txt file.

*/
