#ifndef ALIDNDPTTRACKDUMPTASK_H
#define ALIDNDPTTRACKDUMPTASK_H

//------------------------------------------------------------------------------
// Task to dump track information 
// TPC constrained and TPC+ITS combined 
// for outliers analysis.
// 
// Author: J.Otwinowski 19/06/2011 
//------------------------------------------------------------------------------

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
  static void SetDefaultAliasesHighPt(TTree *treeV0);
 private:

  AliESDEvent *fESD;    //! ESD event
  AliMCEvent *fMC;      //! MC event
  AliESDfriend *fESDfriend; //! ESDfriend event
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
  Double_t fProcessAll; // Calculate all track properties including MC
  
  Bool_t fProcessCosmics; // look for cosmic pairs from random trigger
  Bool_t fProcessITSTPCmatchOut;  // swittch to process ITS/TPC standalone tracks

  TTree* fHighPtTree;       //! list send on output slot 0
  TTree* fV0Tree;           //! list send on output slot 0
  TTree* fdEdxTree;         //! list send on output slot 0
  TTree* fLaserTree;        //! list send on output slot 0
  TTree* fMCEffTree;        //! list send on output slot 0
  TTree* fCosmicPairsTree;  //! list send on output slot 0

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
