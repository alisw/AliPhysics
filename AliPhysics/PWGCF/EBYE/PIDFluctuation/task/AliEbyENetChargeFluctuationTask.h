#ifndef AliEbyENetChargeFluctuationTask_cxx
#define AliEbyENetChargeFluctuationTask_cxx

//=========================================================================//
//                                                                         //
//             AliEbyE Analysis for Net-Charge     Fluctuation             //
//              Author:   Deepika Rathee  || Satyajit Jena                 //
//                        drathee@cern.ch || sjena@cern.ch                 //
//                                                                         //
//=========================================================================//

#include "THnBase.h"
#include "THn.h"
#include "TH1F.h"
#include "TF1.h"
#include "TProfile2D.h"
#include "TRandom3.h"

class TList;

class AliESDtrack;
class AliMCEvent;
class AliStack;
class AliPIDResponse;
class AliESDtrackCuts;
class AliInputEventHandler;
class AliESDInputHandler;
class AliAODInputHandler;
class AliAODEvent;
class AliAODTrack;
class AliAODMCParticle;
class AliMCParticle;
class TClonesArray;
class AliHelperPID;

#include "AliAnalysisTaskSE.h"

class AliEbyENetChargeFluctuationTask: public AliAnalysisTaskSE {

 public:

  AliEbyENetChargeFluctuationTask( const char *name = "NetChargeFluctuation");
  virtual ~AliEbyENetChargeFluctuationTask();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  void SetVertexDiamond(Double_t vx, Double_t vy, Double_t vz) {
    fVxMax = vx;fVyMax = vy; fVzMax = vz;
  }
  void SetKinematicsCuts(Double_t ptl, Double_t pth, Double_t eta, Double_t rap) {
    fPtMin = ptl; fPtMax = pth; fEtaMin = -eta; fEtaMax = eta; fRapMin = -rap; fRapMax = rap;
  }
  void SetTrackFilterBit(Int_t bit) {fAODtrackCutBit = bit; }
  void SetCentralityEstimator(const char* cent)  { fCentralityEstimator = cent;}
  void SetSystemType(Int_t i) { fSystemType = i; }
  void SetEventSelectionBit(UInt_t val) { fSelectBit = val; }
  void SetPhi(Double_t phil) {fPhiMax = phil;}
  void SetIsMC() {fIsMC = kTRUE;}
  void SetNSubSamples(Int_t i) {fNSubSamples = i;}
  void Debug() {fDebug = kTRUE;}
  void DoBasicQA() {fNeedQa = kTRUE;}
  void SetIsAOD() {fIsAOD = kTRUE;}
  void SetAnal(Int_t i, Int_t j);// {fAnalType = i; };
  void SetHelperPID(AliHelperPID* pid){ fHelperPID = pid; }
  void SetDca(Double_t dcaxy,Double_t dcaz) { fDcaXy = dcaxy; fDcaZ = dcaz;}
 
  void SetAnalysisCutObject(AliESDtrackCuts *const trackCuts) {
    fESDtrackCuts = trackCuts;}
  

 private:

  Int_t GetPDG(Int_t i) {
    if (i == 0) return 0;  else if (i == 1) return  211;  
    else if (i == 2) return  321; else if (i == 3) return  2212;  
    else return 0; }
  Double_t *CreateLogAxis(Int_t nbins, Double_t xmin, Double_t xmax);
  
  /*---- Setup events ----*/

  Int_t SetupEvent();
  Int_t SetupESD();
  Int_t SetupAOD();
  Int_t SetupMC();
  Int_t SetupEventCR(AliESDInputHandler *esdHandler, AliAODInputHandler *aodHandler, AliMCEvent *mcEvent);
  void  ResetCurrentEvent();
  /* ---- Check Events/Track/AOD/ESD ---- */
  
  Bool_t ChargedTrack(AliVTrack* track);
  Bool_t TriggeredEvents();
  Bool_t RejectedEvent();
  Bool_t ParticleRapidity(AliVParticle *particle, Double_t &yP, Int_t gCurPid);
  Bool_t TrackRapidity(AliVTrack *track, Double_t &yP, Int_t gCurPid);
  Bool_t IsFindableInTPC(Int_t label);
  Bool_t AcceptEvent(AliAODEvent *event, Int_t cent) const; //! accept eventc
  Bool_t IsEventStats(Int_t *aEventCuts);
  Bool_t AcceptTrackL(AliVTrack *track) const;
  Bool_t AcceptTrackLDCA(AliVTrack *track) const;
  Bool_t AcceptTrackLMC(AliVParticle *particle, Int_t idxMC) const;

  /* ------ Containers ------ */
  
  void InitPhy();  
  void CreateQA();
  void CreateBasicQA();
  void CreateSourceHistos(const Char_t *title, Bool_t isMC);
  void CreateBasicHistos(const Char_t *title, Bool_t isMC, Bool_t isPer);
  void CreateRatioHistos(const Char_t *title,Bool_t isMC, Bool_t isPer);
  void CreateGroupHistos(const Char_t *name, const Char_t *title, Int_t nSample,Bool_t isPer);
  void FillSourceHistos(Bool_t isMC);
  void FillBasicHistos(Bool_t isMC, Bool_t isPer);
  void FillRatioHistos(Bool_t isMC,Bool_t isPer);
  void FillGroupHistos(const Char_t *name, Int_t iSub, Bool_t isMC,Bool_t isPer);
  
 
  void FillCC(Int_t i);
  void FillRCC(Int_t i);
  void FillRecCE(Int_t i);
  void FillRecDE(Int_t i);
  void FillRecDED(Int_t i);
  void FillGenCE(Int_t i);
  void FillQAThnRec(AliVTrack *track, Int_t gPid, Double_t rap);
  void FillQAThnMc(AliVParticle *particle, Int_t gPid, Double_t rap);

  /* ------ Get-Calculate-Fill ------ */

  void ExecAA();
  void ExecpA();
  void Execpp();
  void CalEC();
  void CalED();
  void CreateCE(); 
  void CreateDEM(); 
  void CreateDED(); 
  void CalculateCE(Int_t gPid);
  void CalculateDE(Int_t gPid);
  void CalculateDED(Int_t gPid);
  
  /* Global Members Private */

  AliInputEventHandler *fInputEventHandler; // Input Handler
  AliESDEvent          *fESD;               // Current ESD events
  AliAODEvent          *fAOD;               // Current AOD events
  AliMCEvent           *fMCEvent;           // Current MC Event
  AliStack             *fStack;             // Stak tree
  AliAODInputHandler   *fAODHandler;        // AOD handler
  AliESDInputHandler   *fESDHandler;        // ESD handler
  AliStack             *fMCStack;           // MC Stack
  TClonesArray         *fArrayMC;           // AOD MC stack
  AliESDtrackCuts      *fESDtrackCuts;      // ESD Track Cuts

  TList     *fQaList;                       // Tree of QA
  TList     *fPhyList;                      // For Analysis
  TList     *fDcaList;                      // dca for Both Data and MC 
  TList     *fEffList;                      // Correction List

  Int_t      fSystemType;                   // Translated to Int_t 0:"pp", 1:"pA", 2:"AA"
  TString    fCentralityEstimator;          // "V0M","TRK","TKL","ZDC","FMD"  
 
  Double_t   fVxMax;                        // X vertex  Range
  Double_t   fVyMax;                        // Y vertex Range
  Double_t   fVzMax;                        // Z vertex Range
  Double_t   fPhiMin;                       // Phi Minimum
  Double_t   fPhiMax;                       // Phi Maximum
  Double_t   fPtMin;                        // pT Minimum
  Double_t   fPtMax;                        // Pt Maximum
  Double_t   fEtaMin;                       // Eta Minimum
  Double_t   fEtaMax;                       // Eta Maximum
  Double_t   fRapMin;                       // Rapidity Minimum
  Double_t   fRapMax;                       // Rapidity Maximum
  Double_t   fDcaXy;                        // DCA Xy
  Double_t   fDcaZ;                         // DCA Z
  Double_t   fCentralityBin;                // Centrality bin of current event within max centrality bin
  Double_t   fCentralityPercentile;         // Centrality percentile of current event
  Double_t   fNp[4][2];                     // Array of particle/anti-particle counts
  Double_t   fMCNp[4][2];                   // Array of MC particle/anti-particle counts
  Double_t   fRedFactp[9][2];               // Array of particle/anti-particle reduced factorial
  Double_t   fCurGen[8];                    // Current Gen Track
  Double_t   fCurGenD[8];                   // Current info on DGen
  Double_t   fCurRecD[7];                   // Current DRec
  Double_t   fCurRec[5];                    // Current Rec info
  Double_t   fCurCont[6];                   // Current Contamination

  Float_t    fMinTrackLengthMC;             // Min track length for MC tracks
  UInt_t     fSelectBit;                    // Trigger Bit
  Int_t      fAODtrackCutBit;               // AOD BITs
  Int_t      fNSubSamples;                  // N subsamples
  Int_t      fSubSampleIdx;                 // Subsample idx for current event
  Int_t      fOrder;                        // Max order of higher order distributions
  Int_t      fNTriggers;                    // N triggers used
  Int_t      fHEventStatMax;                // Max N cuts to be included in HEventSta
  Int_t      fNCentralityBins;              // N centrality bins used
  Int_t      fCentralityBinMax;             // Max Cent
  Int_t      fNTracks;                      // Number of Tracks of Current Events
  Int_t      fNbwcBin;                      // FIXME

  Bool_t     fIsMC;                         // Is MC event - Auto set by Add Task
  Bool_t     fIsRatio;                      // Is Ratio
  Bool_t     fIsAOD;                        // analysis mode: 0 = ESDs  | 1 = AODs
  Bool_t     fIsSub;                        // analysis mode SS  
  Bool_t     fIsBS;                         // analysis mode BS
  Bool_t     fIsPer;                        // analysis mode PER
  Bool_t     fIsEff;                        // analysis mode Eff
  Bool_t     fDebug;                        // Debug
  Bool_t     fIsQa;                         // Check for QA
  Bool_t     fNeedQa;                       // All QA When Needed
  Bool_t     fIsPhy;                        // Check for Phy
  Bool_t     fIsDca;                        // Check for Dca
  Bool_t     fIsNu;                         // Check for Nu
  Bool_t     fIsTen;                        // 10% bin 
  Bool_t     fIs3D;                         // 3D Mapping

  TRandom3  *fRan;                          // Radom Number BS
  TRandom3  *fRanIdx;                       // Random Number SS

  AliHelperPID *fHelperPID;                 // Customised HelperPID class

  AliEbyENetChargeFluctuationTask(const AliEbyENetChargeFluctuationTask&);
  AliEbyENetChargeFluctuationTask& operator = (const AliEbyENetChargeFluctuationTask&);
  ClassDef(AliEbyENetChargeFluctuationTask, 1);

};

#endif

 
