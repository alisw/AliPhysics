// Analysis task for creating a reduced tree containing event, track and resonance candidate information 
// Creation date: 2012/06/21 
// Authors: Ionut-Cristian Arsene (i.c.arsene@gsi.de,i.c.arsene@cern.ch)                                 

#ifndef ALIANALYSISTASKREDUCEDTREEMAKER_H
#define ALIANALYSISTASKREDUCEDTREEMAKER_H 

#include <TList.h>
#include <AliAnalysisTaskSE.h>
#include <AliTimeRangeCut.h>
#include <AliSignalMC.h>

class AliAnalysisCuts;
class TTree;
class TFile;
class TBits;
class AliESDtrack;
class AliAODTrack;
class AliESDv0;
class AliAODv0;
class AliESDv0Cuts;
class AliESDv0KineCuts;
class AliAODv0KineCuts;
class AliKFVertex;
class AliReducedBaseEvent;
class AliReducedPairInfo;
class AliAnalysisUtils;
class AliFlowTrackCuts;
class AliMCEvent;

//_________________________________________________________________________
class AliAnalysisTaskReducedTreeMaker : public AliAnalysisTaskSE {
  
public:   
   enum ETreeWritingOptions {
      kBaseEventsWithBaseTracks=0,    // write basic event info and basic track info
      kBaseEventsWithFullTracks,      // write basic event info and full track info
      kFullEventsWithBaseTracks,      // write full event info and base track info
      kFullEventsWithFullTracks,       // write full event info and full track info
      kNMaxTrackFilters=32
   };
   
   enum ETreeMCwritingOptions {
      kBaseTrack=0,
      kFullTrack,
      kMaxMCsignals=32
   };
   
public:
  AliAnalysisTaskReducedTreeMaker();
  AliAnalysisTaskReducedTreeMaker(const char *name, Bool_t writeTree=kTRUE);
  virtual ~AliAnalysisTaskReducedTreeMaker(){  }

  virtual void UserExec(Option_t *option);
  virtual void UserCreateOutputObjects();
  virtual void FinishTaskOutput();
  
  void UsePhysicsSelection(Bool_t phy=kTRUE) {fSelectPhysics=phy;}
  void SetTriggerMask(UInt_t mask) {fTriggerMask=mask;}
  UInt_t GetTriggerMask() const { return fTriggerMask; }
  
  
  // Configure AliAnalysisUtils (used for p-Pb run)
  void SetUseAnalysisUtils(Bool_t flag=kTRUE) {fUseAnalysisUtils=flag;}
  void ConfigAnalysisUtils(Int_t minVtxContrib=1, Double_t maxVtxZ=10.0, Bool_t cutOnSPDVtxZ=kTRUE) {
    fMinVtxContributors = minVtxContrib; fMaxVtxZ=maxVtxZ; fCutOnSPDVtxZ=cutOnSPDVtxZ;
  }
  void SetRejectPileup(Bool_t flag=kTRUE) {fRejectPileup = flag;}
  
  // Cuts for selection of event to be written to tree
  void SetEventFilter(AliAnalysisCuts * const filter) {fEventFilter=filter;}
  void SetTimeRangeReject(Bool_t reject=kTRUE) {fTimeRangeReject = reject;}
  // Cuts for selecting tracks included in the tree
  void AddTrackFilter(AliAnalysisCuts * const filter, Bool_t option=kFALSE, Int_t minSel=-1, Int_t maxSel=-1);
  // Cuts for calorimeter clusters included in the tree
  void AddCaloClusterFilter(AliAnalysisCuts * const filter);
    
  // Cuts for selecting V0s
  void SetK0sPionCuts(AliAnalysisCuts * const filter) {fK0sPionCuts=filter;}
  void SetLambdaProtonCuts(AliAnalysisCuts * const filter) {fLambdaProtonCuts=filter;}
  void SetLambdaPionCuts(AliAnalysisCuts * const filter) {fLambdaPionCuts=filter;}
  void SetGammaElectronCuts(AliAnalysisCuts* const filter) {fGammaElectronCuts=filter;}
  void SetK0sCuts(AliESDv0Cuts* const cuts) {fK0sCuts = cuts;}
  void SetLambdaCuts(AliESDv0Cuts* const cuts) {fLambdaCuts = cuts;}
  void SetGammaConvCuts(AliESDv0KineCuts* const cuts) {fGammaConvCuts = cuts;}
  void SetV0OpenCuts(AliESDv0KineCuts* const cuts) {fV0OpenCuts = cuts;}
  void SetV0StrongCuts(AliESDv0KineCuts* const cuts) {fV0StrongCuts = cuts;}
  void SetV0CutsAOD(AliAODv0KineCuts* const cuts) {fV0CutsAOD = cuts;}
  void SetK0sMassRange(Double_t min=0.4, Double_t max=0.6) {fK0sMassRange[0]=min; fK0sMassRange[1]=max;}
  void SetLambdaMassRange(Double_t min=1.08, Double_t max=1.15) {fLambdaMassRange[0]=min; fLambdaMassRange[1]=max;}
  void SetGammaConvMassRange(Double_t min=0.0, Double_t max=0.1) {fGammaMassRange[0]=min; fGammaMassRange[1]=max;}

  // TStrings with active or inactive branches
  void SetTreeActiveBranch(TString b)   {fActiveBranches+=b+";";}
  void SetTreeInactiveBranch(TString b) {fInactiveBranches+=b+";";}
  
  // Select the type of information to be written
  void SetTreeWritingOption(Int_t option)         {fTreeWritingOption = option;}
  // Suppress writing the tree to disk
  void SetWriteTree(Bool_t option=kTRUE)  {fWriteTree = option;}
  Bool_t WriteTree() const {return fWriteTree;}
  
  // Toggle on/off information branches
  void SetFillTrackInfo(Bool_t flag=kTRUE)        {fFillTrackInfo = flag;}
  void SetFillV0Info(Bool_t flag=kTRUE)           {fFillV0Info = flag;}
  void SetFillGammaConversions(Bool_t flag=kTRUE) {fFillGammaConversions = flag;}
  void SetFillK0s(Bool_t flag=kTRUE)              {fFillK0s = flag;}
  void SetFillLambda(Bool_t flag=kTRUE)           {fFillLambda = flag;}
  void SetFillALambda(Bool_t flag=kTRUE)             {fFillALambda = flag;}
  void SetFillCaloClusterInfo(Bool_t flag=kTRUE)   {fFillCaloClusterInfo = flag;}
  void SetFillFMDInfo(Bool_t flag=kTRUE)               {fFillFMDInfo = flag;}
  //void SetFillBayesianPIDInfo(Bool_t flag=kTRUE)  {fFillBayesianPIDInfo = flag;}
  void SetFillEventPlaneInfo(Bool_t flag=kTRUE, Float_t tpcEtaGap=1.0)    {fFillEventPlaneInfo = flag; fEventPlaneTPCetaGap=tpcEtaGap;}
  void SetFlowTrackFilter(AliAnalysisCuts* const filter)  {fFlowTrackFilter=filter;}
  void SetFillMCInfo(Bool_t flag=kTRUE)               {fFillMCInfo = flag;}
  void AddMCsignal(AliSignalMC* mc, Int_t wOpt=kFullTrack) {
     if(fMCsignals.GetEntries()>=kMaxMCsignals) return; 
     fMCsignals.Add(mc); 
     fMCsignalsWritingOptions[fMCsignals.GetEntries()-1] = wOpt;}
  void SetFillHFInfo(Bool_t flag=kTRUE)               {fFillHFInfo = flag;}
  void SetFillTRDMatchedTracks(Bool_t flag1=kTRUE, Bool_t flag2=kFALSE)   {fFillTRDMatchedTracks = flag1; fFillAllTRDMatchedTracks=flag2;}
  Float_t GetInvPtDevFromBC(Int_t b, Int_t c); // calculates the sagitta value from the online tracks
  void SetWriteUnbiasedEvents(Double_t scaleDown=1.0) {fScaleDownEvents=scaleDown;}

    
 private:

  AliAnalysisUtils* fAnalysisUtils;      // Analysis Utils instance (used to select p-Pb events)
  Bool_t            fUseAnalysisUtils;   // Enable using AnalysisUtils
  Int_t             fMinVtxContributors; // Minimum vtx. contributors used in AliAnalysisUtils
  Double_t          fMaxVtxZ;            // Max. value of the cut on the |z|
  Bool_t            fCutOnSPDVtxZ;       // Use cut on the SPD vtx.
  
  Bool_t fSelectPhysics;             // Whether to use physics selection
  UInt_t fTriggerMask;               // Event trigger mask
  Bool_t fRejectPileup;              // pileup rejection wanted
  
  Int_t     fTreeWritingOption;                 // one of the options described by ETreeWritingOptions
  Bool_t    fWriteTree;                         // if kFALSE don't write the tree, use task only to produce on the fly reduced events
  Double_t  fScaleDownEvents;                      // allow writing events which do not fulfill the minimum number of tracks criteria with scale down factor (default is zero)
  Bool_t    fWriteSecondTrackArray;       // write second array only if full+base tracks requested
  std::vector<Bool_t> fWriteBaseTrack;  // specifier if tracks for certain track filter are reduced or base tracks
  std::vector<Int_t>  fMinSelectedTracks;     // array of min required selected tracks for each track filter
  std::vector<Int_t>  fMaxSelectedTracks;     // array of max required selected tracks for each track filter
  std::vector<Int_t>  fNSelectedFullTracks;       // array of number of full (ReducedTrack) selected tracks in the current event for each track filter + for V0 prongs
  std::vector<Int_t>  fNSelectedBaseTracks;       // array of number of selected base tracks in the current event for each track filter + for V0 prongs
  
  TList* fEventsList;      				// List of event statistics histogram
  TH2I*  fEventsHistogram;            // event statistics histogram
  TH2I*  fTRDEventsHistogram;         // TRD event statistics histogram
  TH2I*  fEMCalEventsHistogram;       // EMCal event statistics histogram
  TList* fCentEventsList; 						// Cent event statistics List of histograms for different estimators
  TH2I*  fTracksHistogram;      			// track statistics histogram
  TH2I*  fMCSignalsHistogram;    			// MC tracks statistics histogram

  Bool_t fFillTrackInfo;            // fill track information
  Bool_t fFillV0Info;               // fill the V0 information
  Bool_t fFillGammaConversions;     // fill gamma conversions
  Bool_t fFillK0s;                  // fill the K0s V0s
  Bool_t fFillLambda;               // fill the lambda V0s
  Bool_t fFillALambda;              // fill the anti-lambda V0s
  Bool_t fFillCaloClusterInfo;      // fill the calorimeter clusters
  Bool_t fFillFMDInfo;              // fill the FMD info
  Bool_t fFillEventPlaneInfo;       // Write event plane information
  Float_t fEventPlaneTPCetaGap;     //  Eta gap between the two TPC sub events
  Bool_t fFillMCInfo;               // Write MC truth information
  Bool_t fFillHFInfo;               // Write HF Process information
  TList   fMCsignals;               // list of AliSignalMC objects to select which particles from the Kinematics stack will be written in the tree
  Int_t   fMCsignalsWritingOptions[kMaxMCsignals];   // writing options for each of the MC signals (either base or full track)     
  Bool_t  fFillTRDMatchedTracks;     // Write global tracks with matched TRD tracks
  Bool_t  fFillAllTRDMatchedTracks;  // if true, fill all global tracks matched in TRD; if false, fill just those global tracks which were selected with the track of V0 filters
  UChar_t fTRDtrglayerMaskEl;       // layer mask for tracklet requirements


  AliAnalysisCuts *fEventFilter;      // event filter
  TList            fTrackFilter;      // filter for the hadrons to be correlated with the dielectrons
  AliAnalysisCuts *fFlowTrackFilter;  // filter for the barrel tracks to be used for the Q-vector
  TList            fClusterFilter;    // filter for the calorimeter clusters

  AliESDv0Cuts *fK0sCuts;            // v0 standard filter for K0s->pi+pi-
  AliESDv0Cuts *fLambdaCuts;         // v0 standard filter for Lambda0->p + pi
  AliESDv0KineCuts *fGammaConvCuts;  // v0 standard filter for gamma conversions
  AliAnalysisCuts *fK0sPionCuts;     // filter for pions from K0s
  AliAnalysisCuts *fLambdaProtonCuts;   // filter for protons from Lambda
  AliAnalysisCuts *fLambdaPionCuts;     // filter for pions from Lambda
  AliAnalysisCuts *fGammaElectronCuts;  // filter for electrons from gamma conversions
  AliESDv0KineCuts *fV0OpenCuts;        // v0 strong filter for tagged V0s
  AliESDv0KineCuts *fV0StrongCuts;      // v0 strong filter for tagged V0s
  AliAODv0KineCuts *fV0CutsAOD;        // filter for tagged V0s in AODs
  
  TH2D* fFMDhist;  // FMD map from AliForwardMult

  Double_t fK0sMassRange[2];        // mass range for allowed K0s pairs
  Double_t fLambdaMassRange[2];     // mass range for allowed Lambda pairs
  Double_t fGammaMassRange[2];      // mass range for allowed Gamma conversion pairs

  TString fActiveBranches;          // list of active output tree branches
  TString fInactiveBranches;        // list of inactive output tree branches

  TFile *fTreeFile;                  //! output file containing the tree
  TTree *fTree;                      //! Reduced event tree
  
  Int_t fNevents;

  AliReducedBaseEvent *fReducedEvent;     //! reduced event wise information
  TBits* fUsedVars;                       // used variables for the AliDielectronVarManager
  
  AliTimeRangeCut     fTimeRangeCut;      //! time range selection based on OADB
  Bool_t              fTimeRangeReject;   //  do not accept events if these are marked by fTimeRangeCut
  
  void FillEventInfo();                     // fill reduced event information
  void FillTrackInfo();                     // fill reduced track information
  void FillMCTruthInfo();                   // fill MC truth particles
  void FillV0PairInfo();                    // fill V0 reduced pair information
  AliReducedPairInfo* FillV0PairInfo(AliESDv0* v0, Int_t id, AliESDtrack* legPos, AliESDtrack* legNeg, AliKFVertex* vtxKF, Bool_t chargesAreCorrect);
  void FillV0PairInfoAOD();                    // fill V0 reduced pair information
  AliReducedPairInfo* FillV0PairInfoAOD(AliAODv0* v0, Int_t id, AliAODTrack* legPos, AliAODTrack* legNeg, AliKFVertex* vtxKF);
  UChar_t EncodeTPCClusterMap(AliVParticle* track, Bool_t isAOD);
  void FillCaloClusters();
  void FillFMDInfo(Bool_t isAOD);
  Int_t GetSPDTrackletMultiplicity(AliVEvent* event, Float_t lowEta, Float_t highEta);
  
  UInt_t MatchMCsignals(Int_t iparticle);     // check all MC signals for MC particle i and write decisions in a bit map
  Bool_t CheckMCtruthWriteFormat(UInt_t bitMap);   // make a decision on which format (base/full) to be used for the MC track
  Bool_t CheckParticleSource(AliMCEvent* event, Int_t ipart, AliSignalMC* mcSignal);
  Bool_t CheckPDGcode(AliMCEvent* event, Int_t ipart, AliSignalMC* mcSignal);
  
  void   FillStatisticsHistograms(Bool_t isSelected, UInt_t physSel, UChar_t trdTrigMap, UInt_t emcalTrigMap, Double_t xbin, Double_t* percentiles, Int_t nEstimators);
  void   FillTrackStatisticsHistogram();
  
  Bool_t  IsTrackSelected(AliVParticle* track, Double_t* values, std::vector<Bool_t>& filterDecision);
  Bool_t  IsSelectedTrackRequestedBaseTrack(std::vector<Bool_t> filterDecision, Bool_t usedForV0Or);
  Bool_t  IsClusterSelected(AliVCluster* cluster, std::vector<Bool_t>& filterDecision);
  void    SetTrackFilterQualityFlags(AliReducedBaseTrack* track, std::vector<Bool_t> filterDecision);
  
  Double_t Rapidity(Double_t r, Double_t z);
  Double_t Radius(Double_t eta, Double_t z);
  
  AliAnalysisTaskReducedTreeMaker(const AliAnalysisTaskReducedTreeMaker &c);
  AliAnalysisTaskReducedTreeMaker& operator= (const AliAnalysisTaskReducedTreeMaker &c);

  ClassDef(AliAnalysisTaskReducedTreeMaker, 18); //Analysis Task for creating a reduced event information tree
};
#endif
