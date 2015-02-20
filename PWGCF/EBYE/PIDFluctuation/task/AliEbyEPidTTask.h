#ifndef AliEbyEPidTTask_cxx
#define AliEbyEPidTTask_cxx

//=========================================================================//
//                                                                         //
//             AliEbyE Analysis for Particle Ratio Fluctuation             //
//              Author:   Deepika Rathee  || Satyajit Jena                 //
//                        drathee@cern.ch || sjena@cern.ch                 //
//                                                                         //
//=========================================================================//

#include "TH1F.h"
#include "TF1.h"

class TList;
class TTree;

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

class AliEbyEPidTTask: public AliAnalysisTaskSE {
 public:
  AliEbyEPidTTask( const char *name = "HigherMomentAnalysis");
  virtual ~AliEbyEPidTTask();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetAODtrackCutBit(Int_t bit) {fAODtrackCutBit = bit; }
  void SetIsMC() {fIsMC = kTRUE;}
  void SetIsAOD() {fIsAOD = kTRUE;}
  void RunQA() {fIsQa = kTRUE;}
  void Debug() {fDebug = kTRUE;}
  void SetHelperPID(AliHelperPID* pid){fHelperPID = pid;}
 void SetAnalysisCutObject(AliESDtrackCuts *const trackCuts) {
   fESDtrackCuts = trackCuts;}
 void SetDca(Double_t dcaxy,Double_t dcaz) { fDcaXy = dcaxy; fDcaZ = dcaz;}
 void SetVertexDiamond(Double_t vx, Double_t vy, Double_t vz) {
   fVxMax = vx;fVyMax = vy; fVzMax = vz;
  }
 void SetKinematicsCuts(Double_t ptl, Double_t pth, Double_t eta) {
   fPtMin = ptl; fPtMax = pth; fEtaMin = -eta; fEtaMax = eta; 
  }
 
 static const Int_t kTrack = 12000;

 private:

  Bool_t RejectedEvent();
  Bool_t TriggeredEvents();
  Bool_t AcceptTrackLMC(AliVParticle *particle) const;
  Bool_t AcceptTrackL(AliVTrack *track) const;
  void ExecEvents();
  void ResetCurrentEvent(); 
   
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
     
  TList        *fThnList;       //!
  Int_t        fAODtrackCutBit; //
  AliHelperPID *fHelperPID;     //
  TH1D         *fEventCounter;  //
  TTree        *fPidCont;       //!
  
  Double_t   fVxMax;                        // X vertex  Range
  Double_t   fVyMax;                        // Y vertex Range
  Double_t   fVzMax;                        // Z vertex Range
  Double_t   fPtMin;                        // pT Minimum
  Double_t   fPtMax;                        // Pt Maximum
  Double_t   fEtaMin;                       // Eta Minimum
  Double_t   fEtaMax;                       // Eta Maximum
  Double_t   fDcaXy;                        // DCA Xy
  Double_t   fDcaZ;                         // DCA Z

  Bool_t     fIsMC;                         // Is MC event - Auto set by Add Task
  Bool_t     fIsAOD;                        // analysis mode: 0 = ESDs  | 1 = AODs
  Bool_t     fDebug;                        // Debug
  Bool_t     fIsQa;                         // Check for QA
  
  Int_t   fRunNumber;           //
  Int_t   fNumberOfTracks;      //
  Int_t   fNumberOfTracksM;     //
  Int_t   fNTracks;            // Number of Tracks of Current Events
  Float_t fCentrality[6];       //
  Float_t fVtx[3];              // 
  Int_t   fTrigMask[5];         //
   
  Float_t fTrackPt[kTrack];     //
  Float_t fTrackEta[kTrack];    //
  Float_t fTrackPhi[kTrack];    //
  Float_t fTrackPtM[kTrack];    //
  Float_t fTrackEtaM[kTrack];   //
  Float_t fTrackPhiM[kTrack];   //
  Float_t fTrackDxy[kTrack];    //
  Float_t fTrackDz[kTrack];     //
  Int_t   fTrackPid[kTrack];    //
  Int_t   fTrackPidM[kTrack];   //
 

  //________________________________
  AliEbyEPidTTask(const AliEbyEPidTTask&);
  AliEbyEPidTTask& operator = (const AliEbyEPidTTask&);
  ClassDef(AliEbyEPidTTask, 1);

};

#endif

 
