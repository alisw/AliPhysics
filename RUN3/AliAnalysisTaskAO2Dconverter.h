/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskAO2Dconverter_H
#define AliAnalysisTaskAO2Dconverter_H

#include "AliAnalysisFilter.h"
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"

#include <TString.h>

#include "TClass.h"

#include <Rtypes.h>

class AliESDEvent;

class AliAnalysisTaskAO2Dconverter : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskAO2Dconverter() = default;
  AliAnalysisTaskAO2Dconverter(const char *name);
  virtual ~AliAnalysisTaskAO2Dconverter();

  AliAnalysisTaskAO2Dconverter(const AliAnalysisTaskAO2Dconverter &) = default;
  AliAnalysisTaskAO2Dconverter &operator=(const AliAnalysisTaskAO2Dconverter &) = delete;

  void SetUseEventCuts(Bool_t useEventCuts=kTRUE) { fUseEventCuts = useEventCuts;}
  Bool_t GetUseEventCuts() const {return fUseEventCuts;}

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  void SetNumberOfEventsPerCluster(int n) { fNumberOfEventsPerCluster = n; }

  static AliAnalysisTaskAO2Dconverter* AddTask(TString suffix = "");
  enum TreeIndex { // Index of the output trees
    kEvents = 0,
    kTracks,
    kCalo,
    kCaloTrigger,
    kMuon,
    kMuonCls,
    kZdc,
    kRun2V0,
    kV0s,
    kCascades,
    kTOF,
    kKinematics,
    kMCvtx,
    kRange,
    kLabels,
    kBC,
    kTrees
  };
  enum TaskModes { // Flag for the task operation mode
    kStandard = 0,
    kMC
  };
  enum MCGeneratorID { // Generator type
    kAliGenEventHeader = 0,
    kAliGenCocktailEventHeader,
    kAliGenDPMjetEventHeader,
    kAliGenEpos3EventHeader,
    kAliGenEposEventHeader,
    kAliGenEventHeaderTunedPbPb,
    kAliGenGeVSimEventHeader,
    kAliGenHepMCEventHeader,
    kAliGenHerwigEventHeader,
    kAliGenHijingEventHeader,
    kAliGenPythiaEventHeader,
    kAliGenToyEventHeader,
    kGenerators
  };
  enum TrackTypeEnum : uint8_t {
    GlobalTrack,
    ITSStandalone,
    MFTStandalone,
    Run2Tracklet
  }; // corresponds to O2/Core/Framework/include/Framework/DataTypes.h
  static const TClass* Generator[kGenerators]; // Generators

  TTree* CreateTree(TreeIndex t);
  void PostTree(TreeIndex t);
  void EnableTree(TreeIndex t) { fTreeStatus[t] = kTRUE; };
  void DisableTree(TreeIndex t) { fTreeStatus[t] = kFALSE; };
  static const TString TreeName[kTrees];  //! Names of the TTree containers
  static const TString TreeTitle[kTrees]; //! Titles of the TTree containers

  void Prune(TString p) { fPruneList = p; }; // Setter of the pruning list
  void SetMCMode() { fTaskMode = kMC; };     // Setter of the MC running mode

  AliAnalysisFilter fTrackFilter; // Standard track filter object
private:
  Bool_t fUseEventCuts = kFALSE;         //! Use or not event cuts
  AliEventCuts fEventCuts;      //! Standard event cuts
  AliESDEvent *fESD = nullptr;  //! input event
  TList *fOutputList = nullptr; //! output list
  
  Int_t fEventCount = 0; //! event count

  // Output TTree
  TTree* fTree[kTrees] = { nullptr }; //! Array with all the output trees
  void Prune();                       // Function to perform tree pruning
  void FillTree(TreeIndex t);         // Function to fill the trees (only the active ones)

  // Task configuration variables
  TString fPruneList = "";                // Names of the branches that will not be saved to output file
  Bool_t fTreeStatus[kTrees] = { kTRUE }; // Status of the trees i.e. kTRUE (enabled) or kFALSE (disabled)
  int fNumberOfEventsPerCluster = 1000;   // Maximum basket size of the trees

  TaskModes fTaskMode = kStandard; // Running mode of the task. Useful to set for e.g. MC mode

  // Data structures

  struct {
    // Start indices and numbers of elements for data in the other trees matching this vertex.
    // Needed for random access of collision-related data, allowing skipping data discarded by the user
    Int_t     fStart[kTrees]    = {0}; /// Start entry indices for data in the other trees matching this vertex
    Int_t     fNentries[kTrees] = {0}; /// Numbers of entries for data in the other trees matching this vertex
    // Event data
    Int_t fBCsID = 0u;       /// Index to BC table
    // Primary vertex position
    Float_t  fPosX = -999.f;       /// Primary vertex x coordinate
    Float_t  fPosY = -999.f;       /// Primary vertex y coordinate
    Float_t  fPosZ = -999.f;       /// Primary vertex z coordinate
    // Primary vertex covariance matrix
    Float_t  fCovXX = 999.f;    /// cov[0]
    Float_t  fCovXY = 0.f;      /// cov[1]
    Float_t  fCovXZ = 0.f;      /// cov[2]
    Float_t  fCovYY = 999.f;    /// cov[3]
    Float_t  fCovYZ = 0.f;      /// cov[4]
    Float_t  fCovZZ = 999.f;    /// cov[5]
    // Quality parameters
    Float_t  fChi2;             /// Chi2 of the vertex
    UInt_t   fN;                /// Number of contributors

    // The calculation of event time certainly will be modified in Run3
    // The prototype below can be switched on request
    Float_t fCollisionTime = -999.f;    /// Event time (t0) obtained with different methods (best, T0, T0-TOF, ...)
    Float_t fCollisionTimeRes = -999.f; /// Resolution on the event time (t0) obtained with different methods (best, T0, T0-TOF, ...)
    UChar_t fCollisionTimeMask = 0u;    /// Mask with the method used to compute the event time (0x1=T0-TOF,0x2=T0A,0x3=TOC) for each momentum bins

  } vtx; //! structure to keep the primary vertex (avoid name conflicts)

  struct {
    int fRunNumber = -1;         /// Run number
    ULong64_t fGlobalBC = 0u;    /// Unique bunch crossing id. Contains period, orbit and bunch crossing numbers
    ULong64_t fTriggerMask = 0u; /// Trigger class mask
  } bc; //! structure to keep trigger-related info
  
  struct {
    // Track data

    Int_t   fCollisionsID;    /// The index of the collision vertex in the TF, to which the track is attached
    
    uint8_t fTrackType;       // Type of track: global, ITS standalone, tracklet, ...
    
    // In case we need connection to TOF clusters, activate next lines
    // Int_t   fTOFclsIndex;     /// The index of the associated TOF cluster
    // Int_t   fNTOFcls;         /// The number of TOF clusters
    
    

    // Coordinate system parameters
    Float_t fX = -999.f;     /// X coordinate for the point of parametrisation
    Float_t fAlpha = -999.f; /// Local <--> global coor.system rotation angle

    // Track parameters
    Float_t fY = -999.f;          /// fP[0] local Y-coordinate of a track (cm)
    Float_t fZ = -999.f;          /// fP[1] local Z-coordinate of a track (cm)
    Float_t fSnp = -999.f;        /// fP[2] local sine of the track momentum azimuthal angle
    Float_t fTgl = -999.f;        /// fP[3] tangent of the track momentum dip angle
    Float_t fSigned1Pt = -999.f;  /// fP[4] 1/pt (1/(GeV/c))

    // Covariance matrix
    Float_t fCYY = -999.f;       /// fC[0]
    Float_t fCZY = -999.f;       /// fC[1]
    Float_t fCZZ = -999.f;       /// fC[2]
    Float_t fCSnpY = -999.f;     /// fC[3]
    Float_t fCSnpZ = -999.f;     /// fC[4]
    Float_t fCSnpSnp = -999.f;   /// fC[5]
    Float_t fCTglY = -999.f;     /// fC[6]
    Float_t fCTglZ = -999.f;     /// fC[7]
    Float_t fCTglSnp = -999.f;   /// fC[8]
    Float_t fCTglTgl = -999.f;   /// fC[9]
    Float_t fC1PtY = -999.f;     /// fC[10]
    Float_t fC1PtZ = -999.f;     /// fC[11]
    Float_t fC1PtSnp = -999.f;   /// fC[12]
    Float_t fC1PtTgl = -999.f;   /// fC[13]
    Float_t fC1Pt21Pt2 = -999.f; /// fC[14]

    // Additional track parameters
    Float_t fTPCinnerP = -999.f; /// Full momentum at the inner wall of TPC for dE/dx PID

    // Track quality parameters
    ULong64_t fFlags = 0u;       /// Reconstruction status flags

    // Clusters and tracklets
    UChar_t fITSClusterMap = 0u;   /// ITS map of clusters, one bit per a layer
    UChar_t fTPCNClsFindable = 0u; /// number of clusters that could be assigned in the TPC
    Char_t fTPCNClsFindableMinusFound = 0;       /// difference between foundable and found clusters
    Char_t fTPCNClsFindableMinusCrossedRows = 0; ///  difference between foundable clsuters and crossed rows
    UChar_t fTPCNClsShared = 0u;   /// Number of shared clusters
    UChar_t fTRDNTracklets = 0u;   /// number of TRD tracklets used for tracking/PID (TRD/TOF pattern)

    // Chi2
    Float_t fITSChi2NCl = -999.f; /// chi2/Ncl ITS
    Float_t fTPCChi2NCl = -999.f; /// chi2/Ncl TPC
    Float_t fTRDChi2 = -999.f;    /// chi2 TRD match (?)
    Float_t fTOFChi2 = -999.f;    /// chi2 TOF match (?)

    // PID
    Float_t fTPCSignal = -999.f; /// dE/dX TPC
    Float_t fTRDSignal = -999.f; /// dE/dX TRD
    Float_t fTOFSignal = -999.f; /// TOFsignal
    Float_t fLength = -999.f;    /// Int.Lenght @ TOF
  } tracks;                      //! structure to keep track information

  struct {
    // MC information on the event
    Short_t fGeneratorsID = 0u; /// Generator ID used for the MC
    Float_t fX = -999.f;  /// Primary vertex x coordinate from MC
    Float_t fY = -999.f;  /// Primary vertex y coordinate from MC
    Float_t fZ = -999.f;  /// Primary vertex z coordinate from MC
    Float_t fT = -999.f;  /// Time of the collision from MC
    Float_t fWeight = -999.f;  /// Weight from MC
    Int_t fNProduced = 0;  /// Number of stable or undecayed particles from MC
  } mcvtx;  //! MC vertices

  struct {
    // Range of labels:
    // for track i the labels are located between fRange[i] and fRange[i+1]
    // so we have Ntracks+1 entries in the table (tree) and fRange[0]=0
    UInt_t fRange = 0; // Current upper limit of the range
  } range; //! Range of labels for each track

  struct {
    // Track labels
    Int_t fLabel = -1;           /// Track label
    // Int_t fTOFLabel[3] = { -1 }; /// Label of the track matched to TOF
  } labels; //! Track labels
  
  struct {
    // MC particle

    Int_t   fCollisionsID;    /// The index of the MC collision vertex

    // MC information (modified version of TParticle
    Int_t fPdgCode    = -99999; /// PDG code of the particle
    Int_t fStatusCode = -99999; /// generation status code
    Int_t fMother[2]   = { 0 }; /// Indices of the mother particles
    Int_t fDaughter[2] = { 0 }; /// Indices of the daughter particles
    Float_t fWeight    = 1;     /// particle weight from the generator or ML

    Float_t fPx = -999.f; /// x component of momentum
    Float_t fPy = -999.f; /// y component of momentum
    Float_t fPz = -999.f; /// z component of momentum
    Float_t fE  = -999.f; /// Energy (covers the case of resonances, no need for calculated mass)

    Float_t fVx = -999.f; /// x of production vertex
    Float_t fVy = -999.f; /// y of production vertex
    Float_t fVz = -999.f; /// z of production vertex
    Float_t fVt = -999.f; /// t of production vertex
    // We do not use the polarisation so far
  } mcparticle;  //! MC particles from the kinematics tree

  // To test the compilation uncoment the line below
  // #define USE_TOF_CLUST 1
#ifdef USE_TOF_CLUST
  struct {
    // TOF clusters
    // PH: Do we store the TOF information per track?
    Int_t fTOFChannel = -1;        /// Index of the matched channel
    Short_t fTOFncls = -1;         /// Number of matchable clusters of the track
    Float_t fDx = -999.f;          /// Residual along x
    Float_t fDz = -999.f;          /// Residual along z
    Float_t fToT = -999.f;         /// ToT
    Float_t fLengthRatio = -999.f; /// Ratio of the integrated track length @ TOF to the cluster with respect to the matched cluster
  } tofClusters;                   //! structure to keep TOF clusters
#endif

  struct {
    // Calorimeter data (EMCAL & PHOS)

    Int_t fBCsID = 0u;       /// Index to BC table

    Short_t fCellNumber = -1;     /// Cell absolute Id. number
    Float_t fAmplitude = -999.f;  /// Cell amplitude (= energy!)
    Float_t fTime = -999.f;       /// Cell time
    Char_t fCellType = -1;        /// EMCAL: High Gain: 0 / Low Gain: 1 / TRU: 2 / LEDmon 3 (see DataFromatsEMCAL/Constants.h)
    Char_t fCaloType = -1;            /// Cell type (-1 is undefined, 0 is PHOS, 1 is EMCAL)
  } calo;                         //! structure to keep EMCAL info
  
  struct {
    // Calorimeter trigger data (EMCAL & PHOS)
    Int_t fBCsID = 0u;        /// Index to BC table
    Short_t fFastOrAbsID = - 1;   /// FastOR absolute ID
    Float_t fL0Amplitude = -1.f;  /// L0 amplitude (ADC) := Peak Amplitude
    Float_t fL0Time = -1.f;       /// L0 time
    Int_t fL1TimeSum = -1;        /// L1 amplitude (ADC) := Integral over L0 time samples
    Char_t fNL0Times = -1;        /// Number of L0 times
    Int_t fTriggerBits = 0;       /// Online trigger bits
    Char_t fCaloType = -1;            /// Calorimeter type (-1 is undefined, 0 is PHOS, 1 is EMCAL)
  } calotrigger;                  //! structure to keep calo trigger info

  struct {
    // MUON track data

    Int_t fBCsID = 0u;            /// Index to BC table
    // In case we need connection to muon clusters, activate next lines
    // Int_t   fClusterIndex;        /// The index of the associated MUON clusters
    // Int_t   fNclusters;           /// The number of MUON clusters

    /// Parameters at vertex
    Float_t fInverseBendingMomentum; ///< Inverse bending momentum (GeV/c ** -1) times the charge 
    Float_t fThetaX;                 ///< Angle of track at vertex in X direction (rad)
    Float_t fThetaY;                 ///< Angle of track at vertex in Y direction (rad)
    Float_t fZMu;                    ///< Z coordinate (cm)
    Float_t fBendingCoor;            ///< bending coordinate (cm)
    Float_t fNonBendingCoor;         ///< non bending coordinate (cm)

    /// Reduced covariance matrix of UNCORRECTED track parameters, ordered as follow:      <pre>
    /// [0] =  <X,X>
    /// [1] =<X,ThetaX>  [2] =<ThetaX,ThetaX>
    /// [3] =  <X,Y>     [4] =  <Y,ThetaX>     [5] =  <Y,Y>
    /// [6] =<X,ThetaY>  [7] =<ThetaX,ThetaY>  [8] =<Y,ThetaY>  [9] =<ThetaY,ThetaY>
    /// [10]=<X,InvP_yz> [11]=<ThetaX,InvP_yz> [12]=<Y,InvP_yz> [13]=<ThetaY,InvP_yz> [14]=<InvP_yz,InvP_yz>  </pre>
    Float_t fCovariances[15]; ///< \brief reduced covariance matrix of parameters AT FIRST CHAMBER

    /// Global tracking info
    Float_t fChi2;                ///< chi2 in the MUON track fit
    Float_t fChi2MatchTrigger;    ///< chi2 of trigger/track matching
  } muons;                        //! structure to keep muons information

  struct {
    // Muon cluster data
    
    Int_t   fMuonsID; /// The index of the muon track to which the clusters are attached
    Float_t fX;         ///< cluster X position
    Float_t fY;         ///< cluster Y position
    Float_t fZ;         ///< cluster Z position
    Float_t fErrX;      ///< transverse position errors
    Float_t fErrY;      ///< transverse position errors
    Float_t fCharge;    ///< cluster charge
    Float_t fChi2;      ///< cluster chi2
  } mucls;              //! structure to keep muon clusters information

  struct {
    // ZDC: it is not clear what is the minimal set of information (PH)

    Int_t fBCsID = 0u;       /// Index to BC table

    Float_t   fZEM1Energy;   	     ///< E in ZEM1
    Float_t   fZEM2Energy;	     ///< E in ZEM2

    Float_t   fZNCTowerEnergy[5];    ///< E in 5 ZNC sectors - high gain chain
    Float_t   fZNATowerEnergy[5];    ///< E in 5 ZNA sectors - high gain chain
    Float_t   fZPCTowerEnergy[5];    ///< E in 5 ZPC sectors - high gain chain
    Float_t   fZPATowerEnergy[5];    ///< E in 5 ZPA sectors - high gain chain

    Float_t   fZNCTowerEnergyLR[5];  ///< E in 5 ZNC sectors - low gain chain
    Float_t   fZNATowerEnergyLR[5];  ///< E in 5 ZNA sectors - low gain chain
    Float_t   fZPCTowerEnergyLR[5];  ///< E in 5 ZPC sectors - low gain chain
    Float_t   fZPATowerEnergyLR[5];  ///< E in 5 ZPA sectors - low gain chain

    Float_t   fZDCTDCCorrected[32][4]; /// ZDC TDC data in ns corrected 4 phase shift

    UChar_t   fFired;                  /// Bits: 0 - ZNA, 1 - ZNC, 2 - ZPA, 3 - ZPC, 4 - ZEM1, 5 - ZEM2
  } zdc;                             //! structure to keep ZDC information

  struct {
    /// Run 2 VZERO Legacy table 

    Int_t fBCsID = 0u;       /// Index to BC table

    Float_t fAdc[64];          ///  adc for each channel
    Float_t fTime[64];         ///  time for each channel
    Float_t fWidth[64];        ///  time width for each channel
    ULong64_t fBBFlag;         ///  BB Flags from Online V0 Electronics
    ULong64_t fBGFlag;         ///  BG Flags from Online V0 Electronics
  } vzero;                     //! structure to keep VZERO information

  struct {
    /// V0s (Ks, Lambda)

    Int_t fPosTrackID; // Positive track ID
    Int_t fNegTrackID; // Negative track ID
  } v0s;               //! structure to keep v0sinformation

  struct {
    /// Cascades

    Int_t fV0sID; // V0 ID
    Int_t fTracksID; // Bachelor track ID
  } cascs;             //! structure to keep cascades information

  /// Offsets to convert the IDs within one collision to global IDs
  Int_t fOffsetMuTrackID = 0; ///! Offset of MUON track IDs (used in the clusters)
  Int_t fOffsetTrackID = 0;   ///! Offset of track IDs (used in V0s)
  Int_t fOffsetV0ID = 0;      ///! Offset of track IDs (used in cascades)
  Int_t fOffsetLabel = 0;      ///! Offset of track IDs (used in cascades)

  ClassDef(AliAnalysisTaskAO2Dconverter, 7);
};

#endif
