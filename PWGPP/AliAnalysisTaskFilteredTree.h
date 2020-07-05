#ifndef ALIDNDPTTRACKDUMPTASK_H
#define ALIDNDPTTRACKDUMPTASK_H


//------------------------------------------------------------------------------
/*
   Class to process/filter reconstruction information from ESD, ESD friends, MC and provide them for later reprocessing
   Filtering schema - low pt part is downscaled - to have flat pt spectra of selected topologies (tracks and V0s)
   Downscaling schema is controlled by downscaling factors
   Usage:
     1.) Filtering on Lego train
     2.) expert QA for tracking (resolution efficiency)
     3.) pt resolution studies using V0s
     4.) dEdx calibration using V0s
     5.) pt resolution and dEdx studies using cosmic
     +
     6.) Info used for later raw data OFFLINE triggering  (highPt, V0, laser, cosmic, high dEdx)

   Exported trees (with full objects and derived variables):
   1.) "highPt"     - filtered trees with esd tracks, derived variables(propagated tracks), optional MC info +optional space points
   2.) "V0s" -      - filtered trees with selected V0s (rough KF chi2 cut), KF particle and corresponding esd tracks + optional space points
   3.) "Laser"      - dump laser tracks with space points if exists
   4.) "CosmicTree" - cosmic track candidate (random or triggered) + esdTracks(up/down)+ optional points
   5.) "dEdx"       - tree with high dEdx tpc tracks
*/
class AliESDEvent;
class AliMCEvent;
class AliFilteredTreeEventCuts;
class AliFilteredTreeAcceptanceCuts;
class AliESDtrackCuts;
class AliMagFMaps;
class AliESDEvent; 
class AliMCEvent; 
class AliKFParticle; 
class AliESDv0; 
class AliExternalTrackParam;
class AliESDtrack;
class AliESDfriendTrack;
class AliESDVertex;
class AliStack;
class TList;
class TObjArray;
class TTree;
class TTreeSRedirector;
class TParticle;
class TH3D;
class AliESDtools;
#include <string>

#include "AliTriggerAnalysis.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskFilteredTree : public AliAnalysisTaskSE {
 public:

  enum EAnalysisMode { kInvalidAnalysisMode=-1,
                      kTPCITSAnalysisMode=0,
                      kTPCAnalysisMode=1 };

  AliAnalysisTaskFilteredTree(const char *name = "AliAnalysisTaskFilteredTree");
  virtual ~AliAnalysisTaskFilteredTree();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  virtual Bool_t Notify();
  virtual void   FinishTaskOutput();
  void SetUseMCInfo(Bool_t info)           { fUseMCInfo = info; }
  Bool_t IsUseMCInfo() const               { return (fMC)?kTRUE:kFALSE; }
  void SetUseESDfriends(Bool_t friends)    { fUseESDfriends = friends; }
  Bool_t IsUseESDfriends() const              { return fUseESDfriends; }
  
  // Process events
  void ProcessAll(AliESDEvent *const esdEvent=0, AliMCEvent *const mcEvent=0, AliESDfriend *const esdFriend=0);
  void Process(AliESDEvent *const esdEvent=0, AliMCEvent *const mcEvent=0, AliESDfriend *const esdFriend=0);
  void ProcessV0(AliESDEvent *const esdEvent=0, AliMCEvent *const mcEvent=0, AliESDfriend *const esdFriend=0);
  void ProcessdEdx(AliESDEvent *const esdEvent=0, AliMCEvent *const mcEvent=0, AliESDfriend *const esdFriend=0);
  void ProcessLaser(AliESDEvent *const esdEvent=0, AliMCEvent *const mcEvent=0, AliESDfriend *const esdFriend=0);
  void ProcessMCEff(AliESDEvent *const esdEvent=0, AliMCEvent *const mcEvent=0, AliESDfriend *const esdFriend=0);
  void ProcessCosmics(AliESDEvent *const esdEvent=0, AliESDfriend* esdFriend=0); 
  void ProcessMC();  //TODO  - not yet finished/tested/enabled

  void ProcessITSTPCmatchOut(AliESDEvent *const esdEvent=0,  AliESDfriend *const esdFriend=0);
  void ProcessTrackMatch(AliESDEvent *const esdEvent=0,  AliESDfriend *const esdFriend=0);

  void SetEventCuts(AliFilteredTreeEventCuts* const cuts)              { fFilteredTreeEventCuts = cuts; }
  void SetAcceptanceCuts(AliFilteredTreeAcceptanceCuts* const cuts)    { fFilteredTreeAcceptanceCuts = cuts; }
  void SetRecAcceptanceCuts(AliFilteredTreeAcceptanceCuts* const cuts) { fFilteredTreeRecAcceptanceCuts = cuts; }
  void SetTrackCuts(AliESDtrackCuts* const cuts)                { fEsdTrackCuts = cuts; }
  void SetTrigger(const AliTriggerAnalysis::Trigger trigger)    { fTrigger = trigger; }
  void SetAnalysisMode(EAnalysisMode mode) { fAnalysisMode = mode; }

  AliFilteredTreeEventCuts* GetEventCuts() const                       { return fFilteredTreeEventCuts; }
  AliFilteredTreeAcceptanceCuts* GetAcceptanceCuts() const             { return fFilteredTreeAcceptanceCuts; }
  AliFilteredTreeAcceptanceCuts* GetRecAcceptanceCuts() const          { return fFilteredTreeRecAcceptanceCuts; }  
  AliESDtrackCuts* GetTrackCuts() const                         { return fEsdTrackCuts; }
  AliTriggerAnalysis::Trigger GetTrigger() const                { return fTrigger; }
  EAnalysisMode GetAnalysisMode() const          { return fAnalysisMode; }

  TString GetCentralityEstimator() const {return fCentralityEstimator; }
  void SetCentralityEstimator(TString centEst="V0M") { fCentralityEstimator = centEst; }

  Bool_t IsFromConversion(Int_t label, AliStack *const stack);
  Bool_t IsFromMaterial(Int_t label, AliStack *const stack);
  Bool_t IsFromStrangeness(Int_t label, AliStack *const stack);
  TParticle *GetMother(TParticle *const particle, AliStack *const stack);

  Bool_t ConstrainTPCInner(AliExternalTrackParam *const tpcInnerC, const AliESDVertex* vtx, Double_t b[3]);
  Bool_t ConstrainTrackInner(AliExternalTrackParam *const trackInnerC, const AliESDVertex* vtx, Double_t mass, Double_t b[3]);

  // v0s selection
  Int_t  GetKFParticle(AliESDv0 *const v0, AliESDEvent * const event, AliKFParticle & kfparticle);
  Bool_t IsV0Downscaled(AliESDv0 *const v0);
  Int_t  V0DownscaledMask(AliESDv0 *const v0);
  Bool_t IsHighDeDxParticle(AliESDtrack * const track);

  void SetLowPtTrackDownscaligF(Double_t fact) { fLowPtTrackDownscaligF = fact; }
  void SetLowPtV0DownscaligF(Double_t fact)    { fLowPtV0DownscaligF = fact; }
  void SetFriendDownscaling(Double_t fact)    { fFriendDownscaling = fact; }
  
  void   SetProcessCosmics(Bool_t flag) { fProcessCosmics = flag; }
  Bool_t GetProcessCosmics() { return fProcessCosmics; }
  //
  void   SetProcessProcessITSTPCmatchOut(Bool_t flag) { fProcessITSTPCmatchOut = flag; }
  Bool_t GetProcessProcessITSTPCmatchOut() { return fProcessITSTPCmatchOut; }

  
  void SetProcessAll(Bool_t proc) { fProcessAll = proc; }
  static Int_t GetMCTrueTrackMult(AliMCEvent *const mcEvent, AliFilteredTreeEventCuts *const evtCuts, AliFilteredTreeAcceptanceCuts *const accCuts);

  void SetFillTrees(Bool_t filltree) { fFillTree = filltree ;}
  Bool_t GetFillTrees() { return fFillTree ;}

  void FillHistograms(AliESDtrack* const ptrack, AliExternalTrackParam* const ptpcInnerC, Double_t centralityF, Double_t chi2TPCInnerC);
  Int_t   GetNearestTrack(const AliExternalTrackParam * trackMatch, Int_t indexSkip, AliESDEvent*event, Int_t trackType, Int_t paramType,  AliExternalTrackParam & paramNearest);
  static void SetDefaultAliasesV0(TTree *treeV0);
  static void  SetDefaultAliasesV0PID(TTree *treeV0, Int_t pidHash);
  static void SetDefaultAliasesHighPt(TTree *tree);
  static void SetDefaultAliasesEvents(TTree *treeEvent);
  Int_t GetMCInfoTrack(Int_t label,   std::map<std::string,float> &trackInfoF, std::map<std::string,TObject*> &trackInfoO);  //TODO- test before enabling
  Int_t GetMCInfoKink(Int_t label,    std::map<std::string,float> &kinkInfoF, std::map<std::string,TObject*> &kinkInfoO);  // TODO
  static Int_t GetMCTrackDiff(const TParticle &particle, const AliExternalTrackParam &param, TClonesArray &trackRefArray, TVectorF &mcDiff); //TODO test before enabling
  /// sqrt s - mass dependent downsampling trigger (pt spectra as parameterized in https://iopscience.iop.org/article/10.1088/2399-6528/aab00f/pdf)
  static Double_t TsalisCharged(Double_t pt, Double_t mass, Double_t sqrts);
  static Int_t    DownsampleTsalisCharged(Double_t pt, Double_t factorPt, Double_t factor1Pt,  Double_t sqrts=5020, Double_t mass=0.2);
  Int_t  PIDSelection(AliESDtrack *track, TParticle *particle = nullptr);
 private:
  AliESDEvent *fESD;    //! ESD event
  AliMCEvent *fMC;      //! MC event
  AliESDfriend *fESDfriend; //! ESDfriend event
  AliESDtools *fESDtool;      /// tools to calculate derived variables from the ESD
  TList* fOutput;       //! list send on output slot 0
  TIterator *fPitList;  //! iterator over the output objetcs
  Bool_t fUseMCInfo;        // use MC information
  Bool_t fUseESDfriends;    // use esd friends
  Bool_t fReducePileUp;     // downscale the information for the pile-up TPC tracks
  Bool_t fFillTree;         // do not fill trees

  AliFilteredTreeEventCuts      *fFilteredTreeEventCuts;      // event cuts
  AliFilteredTreeAcceptanceCuts *fFilteredTreeAcceptanceCuts; // acceptance cuts  
  AliFilteredTreeAcceptanceCuts *fFilteredTreeRecAcceptanceCuts; // additional recontruction acceptance cuts (not used for MC truth)
  AliESDtrackCuts *fEsdTrackCuts;          // esd track cuts
  AliTriggerAnalysis::Trigger fTrigger;    // trigger settings
  EAnalysisMode fAnalysisMode;   // analysis mode TPC only, TPC + ITS

  TTreeSRedirector* fTreeSRedirector;      //! temp tree to dump output

  TString fCentralityEstimator;     // use centrality can be "VOM" (default), "FMD", "TRK", "TKL", "CL0", "CL1", "V0MvsFMD", "TKLvsV0M", "ZEMvsZDC"

  Double_t fLowPtTrackDownscaligF; // low pT track downscaling factor
  Double_t fLowPtV0DownscaligF;    // low pT V0 downscaling factor
  Double_t fFriendDownscaling;     // friend info downscaling )absolute value used), Modes>=1 downscaling in respect to the amount of tracks, Mode<=-1 (downscaling in respect to the data volume)
  Double_t fSqrtS;                 // sqrt(s) used for downsampling to approximate spectra function
  Double_t fChargedEffectiveMass;           // mass used for downsampling to approximate spectra function (pion,Kaon,prootn)
  Double_t fV0EffectiveMass;           // mass used for downsampling to approximate spectra function for V0 (K0s,Lambda)
  Double_t fProcessAll; // Calculate all track properties including MC
  
  Bool_t fProcessCosmics; // look for cosmic pairs from random trigger
  Bool_t fProcessITSTPCmatchOut;  // switch to process ITS/TPC standalone tracks

  TTree* fHighPtTree;       //! list send on output slot 0
  TTree* fV0Tree;           //! list send on output slot 0
  TTree* fdEdxTree;         //! list send on output slot 0
  TTree* fLaserTree;        //! list send on output slot 0
  TTree* fMCEffTree;        //! list send on output slot 0
  TTree* fCosmicPairsTree;  //! list send on output slot 0
  TH1F * fSelectedTracksMask;   //! histogram of the selected tracks
  TH1F * fSelectedPIDMask;   //! histogram of the selected tracks
  TH1F * fSelectedV0Mask;   //! histogram of the selected tracks

  TH3D* fPtResPhiPtTPC;    //! sigma(pt)/pt vs Phi vs Pt for prim. TPC tracks
  TH3D* fPtResPhiPtTPCc;   //! sigma(pt)/pt vs Phi vs Pt for prim. TPC contrained to vertex tracks
  TH3D* fPtResPhiPtTPCITS; //! sigma(pt)/pt vs Phi vs Pt for prim. TPC+ITS tracks

  TH3D* fPtResEtaPtTPC;    //! sigma(pt)/pt vs Eta vs Pt for prim. TPC tracks
  TH3D* fPtResEtaPtTPCc;   //! sigma(pt)/pt vs Eta vs Pt for prim. TPC contrained to vertex tracks
  TH3D* fPtResEtaPtTPCITS; //! sigma(pt)/pt vs Eta vs Pt for prim. TPC+ITS tracks

  TH3D* fPtResCentPtTPC;    //! sigma(pt)/pt vs Cent vs Pt for prim. TPC tracks
  TH3D* fPtResCentPtTPCc;   //! sigma(pt)/pt vs Cent vs Pt for prim. TPC contrained to vertex tracks
  TH3D* fPtResCentPtTPCITS; //! sigma(pt)/pt vs Cent vs Pt for prim. TPC+ITS tracks
  TObjString fCurrentFileName; // cached value of current file name
  AliESDtrack* fDummyTrack; //! dummy track for tree init

  AliAnalysisTaskFilteredTree(const AliAnalysisTaskFilteredTree&); // not implemented
  AliAnalysisTaskFilteredTree& operator=(const AliAnalysisTaskFilteredTree&); // not implemented
  ClassDef(AliAnalysisTaskFilteredTree, 1); // example of analysis
};

#endif
