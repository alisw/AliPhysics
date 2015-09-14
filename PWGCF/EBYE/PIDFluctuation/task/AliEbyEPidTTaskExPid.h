#ifndef AliEbyEPidTTaskExPid_cxx
#define AliEbyEPidTTaskExPid_cxx

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
class TH2F;

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


#include "AliAnalysisTaskSE.h"

class AliEbyEPidTTaskExPid: public AliAnalysisTaskSE {
 public:
  AliEbyEPidTTaskExPid( const char *name = "HigherMomentAnalysis");
  virtual ~AliEbyEPidTTaskExPid();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetAODtrackCutBit(Int_t bit) {fAODtrackCutBit = bit; }
  void SetIsMC() {fIsMC = kTRUE;}
  void SetIsAOD() {fIsAOD = kTRUE;}
  
  void SetAnalysisCutObject(AliESDtrackCuts *const trackCuts) {
   fESDtrackCuts = trackCuts;}
  void SetKinematicsCuts(Double_t ptl, Double_t pth, Double_t eta) {
    fPtMin = ptl; fPtMax = pth; fEtaMin = -eta; fEtaMax = eta; 
  }
  
  static const Int_t kTrack = 15000;
  static const Int_t   fgkPidEx[4];       // Histogram N bins for eta

 private:
 
  Bool_t AcceptTrackLMC(AliVParticle *particle) const;
  Bool_t AcceptTrackL(AliVTrack *track) const;
  TClonesArray         *fArrayMC;           // AOD MC stack
  AliESDtrackCuts      *fESDtrackCuts;      // ESD Track Cuts
  AliMCEvent           *fMCEvent;           // Current MC Event
  AliStack             *fMCStack;             // Stak tree   
  AliPIDResponse   *fPIDResponse;     //! PID response object
  Int_t        fAODtrackCutBit; //
  TTree        *fPidCont;       //!
  
  Double_t   fPtMin;                        // pT Minimum
  Double_t   fPtMax;                        // Pt Maximum
  Double_t   fEtaMin;                       // Eta Minimum
  Double_t   fEtaMax;                       // Eta Maximum

  Bool_t     fIsMC;                         // Is MC event - Auto set by Add Task
  Bool_t     fIsAOD;                        // analysis mode: 0 = ESDs  | 1 = AODs
  
  Int_t   fRunNumber;           //
  Int_t   fNumberOfTracks;      //
  Int_t   fNumberOfTracksM;     //
  Int_t   fNTracks;            // Number of Tracks of Current Events

  Float_t fCentrality[6];       //
  Float_t fVtx[3];              // 

  Float_t fPParam[kTrack][15];//
  Float_t fTParam[kTrack][11];//

  Int_t   fPidStat[kTrack];     //
  Int_t   fTrackLabel[kTrack]; //
  Int_t   fTrackLabelM[kTrack]; //
  Float_t fTrackPtM[kTrack];    //
  Float_t fTrackEtaM[kTrack];   //
  Float_t fTrackPhiM[kTrack];   //
  Int_t   fTrackPidM[kTrack];   //

  //________________________________
  AliEbyEPidTTaskExPid(const AliEbyEPidTTaskExPid&);
  AliEbyEPidTTaskExPid& operator = (const AliEbyEPidTTaskExPid&);
  ClassDef(AliEbyEPidTTaskExPid, 1);
};

#endif

 
