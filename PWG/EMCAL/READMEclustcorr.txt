/*! \page READMEclustcorr EMCal/DCal cluster corrections
The AliAODCaloCluster and AliESDCaloCluster objects can store different level of corrections to the cluster energy. The "bare" energy is obtained using the method `cluster->E()`. This energy usually:

- implements all the basic energy/time calibrations
- implements bad channel removal
- does not include non-linearity correction
- does not include exotic cell removal
- does not include any analysis-specific corrections (such as the "hadronic correction")

For some dataset one needs to re-run the clusterizer from the cells. This is dataset-specific. Analyzers that are unfamiliar with a specific dataset should communicate with EMCal/DCal detector experts to determine whether the basic corrections and the bad channel map were already available and used in a specific ESD/AOD production.

In the following a workflow is suggested. This should be good for most analysis but it is not the only possible one. A lot of code is in fact duplicated and the same result can be obtained using different pieces of the framework.

# EMCal Tender

<strong>Class names</strong>: AliTender, AliEmcalTenderTask (PWG/EMCAL) and AliEMCALTenderSupply (TENDER/TenderSupplies).

<strong>Add task macro</strong>: PWG/EMCAL/macros/AddTaskEMCALTender.C

~~~{.cxx}
AliAnalysisTaskSE *AddTaskEMCALTender(
  Bool_t distBC         = kTRUE,   //distance to bad channel
  Bool_t recalibClus    = kTRUE,   //recalibrate cluster energy
  Bool_t recalcClusPos  = kTRUE,   //recalculate cluster position
  Bool_t nonLinearCorr  = kTRUE,   //apply non-linearity
  Bool_t remExoticCell  = kTRUE,   //remove exotic cells
  Bool_t remExoticClus  = kTRUE,   //remove exotic clusters
  Bool_t fidRegion      = kFALSE,  //apply fiducial cuts
  Bool_t calibEnergy    = kTRUE,   //calibrate energy
  Bool_t calibTime      = kTRUE,   //calibrate timing
  Bool_t remBC          = kTRUE,   //remove bad channels
  UInt_t nonLinFunct    = AliEMCALRecoUtils::kBeamTestCorrected,
  Bool_t reclusterize   = kTRUE,   //reclusterize
  Float_t seedthresh    = 0.100,   //seed threshold
  Float_t cellthresh    = 0.050,   //cell threshold
  UInt_t clusterizer    = AliEMCALRecParam::kClusterizerv2,
  Bool_t trackMatch     = kTRUE,   //track matching
  Bool_t updateCellOnly = kFALSE,  //only change if you run your own clusterizer task
  Float_t timeMin       = 100e-9,  //minimum time of physical signal in a cell/digit (s)
  Float_t timeMax       = 900e-9,  //maximum time of physical signal in a cell/digit (s)
  Float_t timeCut       = 900e-9,  //maximum time difference between the digits inside EMC cluster (s)
  const char *pass      = 0,       //string defining pass (use none if figured out from path)
  Bool_t  remapMcAod    = kFALSE,  //switch on the remaping for the MC labels in AOD productions,
  TString cdbStorage    = "raw://" // "local://"
)
~~~

<strong>Suggested parameter list</strong>:
~~~
kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kTRUE,kTRUE,kTRUE,0,kFALSE,0.1,0.05,0,kFALSE,kTRUE
~~~

The EMCal tender is used to do **cell level** corrections. This means that, unless you are going to run a clusterizer after the tender, you probably do not need the tender. However in most cases using the tender + clusterizer might be necessary. The tender is used to apply the following corrections:

- apply latest available energy calibration
- remove bad (hot) channels

The tender will calibrate the cells in the standard branch in the ESD/AOD event (respectively "EMCALCells" and "emcalCells"). This means that original cell information in the event **will be overwritten**.

# EMCal clusterizer

<strong>Class name</strong>: AliAnalysisTaskEMCALClusterizerFast (PWG/EMCAL).

<strong>Add task macro</strong>: PWG/EMCAL/macros/AddTaskClusterizerFast.C

~~~{.cxx}
AliAnalysisTaskEMCALClusterizeFast* AddTaskClusterizerFast(
  const char* taskname  = "ClusterizerFast",
  const char* cellsName = "",
  const char* clusName  = "",
  UInt_t clusterizer    = AliEMCALRecParam::kClusterizerv2,
  Double_t cellE        = 0.05,
  Double_t seedE        = 0.1,
  const Float_t timeMin = -1,      //minimum time of physical signal in a cell/digit (s)
  const Float_t timeMax = +1,      //maximum time of physical signal in a cell/digit (s)
  const Float_t timeCut =  1,      //maximum time difference between the digits inside EMC cluster (s)
  Bool_t remExoticCell  = kTRUE,
  Bool_t calcDistToBC   = kFALSE,
  UInt_t inputCellType  = AliAnalysisTaskEMCALClusterizeFast::kFEEData)
~~~

<strong>Suggested parameter list</strong>:
~~~
"ClusterizerFast", "", "", kClusterizerType, 0.05,0.1, kEMCtimeMin, kEMCtimeMax, kEMCtimeCut, kFALSE, kFALSE, AliAnalysisTaskEMCALClusterizeFast::kFEEData
~~~

The clusterizer type and the time cuts are to be chosen appropriately for each dataset. Usually the v1 clusterizer is used for pp and the v2 clusterizer is used for PbPb. Sometimes for pp reference runs with the same collision energy as PbPb the v2 clusterizer is employed, but this is analysis dependent. EMCal detector experts are to be contacted for the time cuts.

The clusterizer will use the cells from the standard branch in the ESD/AOD event (respectively "EMCALCells" and "emcalCells") and rewrite the cluster collection ("CaloCluster" or "caloCluster").

At this point the energy of the cluster will be available through `cluster->E()` where cluster is the pointer to the AliAODCaloCluster or AliESDCaloCluster object.

# Cluster "maker"
<strong>Class name</strong>: AliEmcalClusterMaker (PWG/EMCAL).

<strong>Add task macro</strong>: PWG/EMCAL/macros/AddTaskEmcalClusterMaker.C

~~~{.cxx}
AliEmcalClusterMaker* AddTaskEmcalClusterMaker(const UInt_t nonLinFunct   = AliEMCALRecoUtils::kBeamTestCorrected,
                                               const Bool_t remExClus     = kTRUE,
                                               const char *nClusters      = 0,
                                               const char *outClusName    = "EmcCaloClusters",
                                               const Double_t emin        = 0.3,
                                               const Bool_t   histo       = kFALSE,
                                               const char *outputname     = "AnalysisResults.root"
)
~~~

<strong>Suggested parameter list</strong>: AliEMCALRecoUtils::kBeamTestCorrectedv3, kTRUE, 0, "", 0.

Non-linearity correction and "exotic" cluster removal are performed in the cluster "maker". Non-linearity correction to the cluster energy is necessary because the response of the calorimeter is not linear for very low momentum particles or very high momentum (shower leakage). "Exotic" cluster are energetic clusters where most energy deposition is concentrated in one single cell. This clusters are not reproduced in MC simulations and are believed to arise from neutrons showering directly into the APD. This clusters need to be flagged, so that they can be easily rejected during the analysis.

The energy of the cluster **after** the non-linearity correction can be retrieved using the method `cluster->GetNonLinCorrEnergy()`. The "exotic" flag can be retrieved using `cluster->GetIsExotic()`. "Exotic" clusters can be easily rejected if clusters are accessed using an AliClusterContainer object. "Exotic" cluster removal is switched on by default in AliClusterContainer, however **it is necessary to run the cluster maker to flag "exotic" cluster beforehand**.

# Track-cluster matching
<strong>Class name</strong>: AliEmcalClusTrackMatcherTask (PWG/EMCAL).

<strong>Add task macro</strong>: PWG/EMCAL/macros/AddTaskEmcalClusTrackMatcher.C

~~~{.cxx}
AliEmcalClusTrackMatcherTask* AddTaskEmcalClusTrackMatcher(const char *nTracks          = "usedefault",
                                                           const char *nClusters        = "usedefault",
                                                           const Double_t maxDist       = 0.1,
                                                           const Bool_t attachEmcalPart = kFALSE,
                                                           const Bool_t updateClusters  = kTRUE,
                                                           const Bool_t updateTracks    = kTRUE,
                                                           const Bool_t createHisto     = kFALSE)
~~~

<strong>Suggested parameter list</strong>: use default.

Tracks and clusters are matched using a simple geometrical algorithm. Multiple tracks can be matched to a single clusters; however only one cluster can be matched to a track. The default configuration of the task is such that it will attempt track propagation to the EMCal surface (440 cm) if the track is not already propagated. This means that the OCDB has to be loaded beforehand (e.g. using the AliEmcalSetup task) in order to have access to the ALICE magnetic field and geometry. This should usually work in both AOD and ESD events.

The number of tracks matched to a cluster can be retrieved using `cluster->GetNTracksMatched()`. Unfortunately the method to access the tracks matched to a cluster depend on the data format. For ESD clusters (AliESDCaloClusters):
~~~{.cxx}
Int_t iTrack = cluster->GetTrackMatchedIndex(i);
~~~
will return the position of the track in the array. The integer i is a number from 0 to `cluster->GetNTracksMatched() -1`. A pointer to the track object can be retrieved using:
~~~{.cxx}
AliVTrack* track = static_cast<AliVTrack*>(GetParticleContainer(0)->GetParticle(iTrack));
~~~
(assuming that the task is derived from AliAnalysisTaskEmcal or AliAnalysisTaskEmcalJet).

For AOD clusters (AliAODCaloClusters) the method:
~~~{.cxx}
AliVTrack* track = static_cast<AliVTrack*>(cluster->GetTrackMatched(i));
~~~
will directly return a pointer to the matched track.

To get the cluster matched to a track one can use (both ESD and AOD):
~~~{.cxx}
Int_t iCluster = track->GetEMCALcluster();
~~~
This will return the index of the cluster. To get a pointer to the cluster object:
~~~{.cxx}
AliVCluster *cluster = GetClusterContainer(0)->GetCluster(iCluster);
~~~
(again assuming that the task is derived from AliAnalysisTaskEmcal or AliAnalysisTaskEmcalJet).


# Charged particle correction (aka "hadronic correction")

<strong>Class name</strong>: AliHadCorrTask (PWG/EMCAL).

<strong>Add task macro</strong>: PWG/EMCAL/macros/AddTaskHadCorr.C
~~~{.cxx}
AliHadCorrTask* AddTaskHadCorr(
  const char *nTracks        = "EmcalTracks",
  const char *nClusters      = "EmcalClusters",
  const char *outClusName    = "CaloClustersCorr",
  const Double_t hadcorr     = 1,
  const Double_t minPt       = 0.15,
  const Double_t phiMatch    = 0.050,
  const Double_t etaMatch    = 0.025,
  const Double_t Eexcl       = 0,
  const Bool_t trackClus     = kTRUE,
  const Bool_t   histo       = kFALSE,
  const char *outputname     = "AnalysisResults.root"
)
~~~

<strong>Suggested parameter list</strong>:

~~~
"usedefault", "usedefault", "",  2.0, 0.15, 0.030, 0.015, 0, kTRUE, kTRUE
~~~

Charged particles deposit some energy in the calorimeter. Most of the charged particle are hadrons, such as pions, kaons and protons. The hadronic response of the calorimeter has been studied in some details. Most of the high energetic particles ( > 1 GeV) only release a small amount of energy. These are usually called "minimum ionizing particles" (MIP). Occasionally hadrons may interact strongly with the nuclei of the material in the calorimeter and start a hadronic shower. In this case the energy deposition is much higher. High momentum muons are also MIP, but they never shower in the calorimeter. Finally electrons do shower in the calorimeter, in a way that is quite similar to a shower initiated by a photon.

All charged particles are measured in the tracking detectors (ITS+TPC+TOF). To avoid double counting their contribution to the jet energy flow, they need to be subtracted cluster-by-cluster. This is done by matching EMCal/DCal clusters with charged tracks and then subtracting a certain fraction of the matched tracks from the cluster energy. The most common choice is to subtract 100% of the momentum of the sum of the matched tracks.

The energy of the cluster **after** the hadronic correction can be obtained using the method `cluster->GetHadCorrEnergy()`.
*/
