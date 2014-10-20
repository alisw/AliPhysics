#ifndef AliEbyEPidTTaskMC_cxx
#define AliEbyEPidTTaskMC_cxx

//=========================================================================//
//                                                                         //
//             AliEbyE Analysis for Particle Ratio Fluctuation             //
//              Author:   Deepika Rathee  || Satyajit Jena                 //
//                        drathee@cern.ch || sjena@cern.ch                 //
//                                                                         //
//=========================================================================//

class TH1D;
class TH2F;
class TH3F;
class TString;
class AliAODTrack;
class AliAODMCParticle;
class TList;
class TTree;
class AliHelperPID;

#include "AliAnalysisTaskSE.h"


class AliEbyEPidTTaskMC: public AliAnalysisTaskSE {
 public:
  AliEbyEPidTTaskMC( const char *name = "HigherMomentAnalysis");
  virtual ~AliEbyEPidTTaskMC();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetAODtrackCutBit(Int_t bit) {fAODtrackCutBit = bit; }
  void SetCentralityEstimator(const char* cent) { fCentralityEstimator = cent;}
  void SetAnalysisType(Bool_t ismc) {isMC = ismc;}
  void SetHelperPID(AliHelperPID* pid){fHelperPID = pid;}

  static const Int_t kTrack = 90000;

 private:
 
  Bool_t       AcceptTrack(AliAODTrack *track) const; 
  Bool_t       AcceptMCTrack(AliAODMCParticle *track) const; 
  TList        *fThnList;
  Bool_t       isMC;                   
  TString      fCentralityEstimator;   
  Int_t        fAODtrackCutBit;
  AliHelperPID *fHelperPID;
  TH1D         *fEventCounter;
  TTree        *fEventTree;
  
  Int_t   fRunNumber;
  Int_t   fNumberOfTracks;
  Int_t   fNumberOfTracksM;
  Int_t   fFilterBit;
  Float_t fCentPercentile;
  Float_t fVertexX;
  Float_t fVertexY;
  Float_t fVertexZ;
  
  Float_t fTrackPt[kTrack];
  Float_t fTrackEta[kTrack];
  Float_t fTrackPhi[kTrack];

  Float_t fTrackPtM[kTrack];
  Float_t fTrackEtaM[kTrack];
  Float_t fTrackPhiM[kTrack];

  Int_t   fTrackCharge[kTrack];
  Int_t   fTrackPid[kTrack];

  Int_t   fTrackChargeM[kTrack];
  Int_t   fTrackPidM[kTrack];
 

  //________________________________
  AliEbyEPidTTaskMC(const AliEbyEPidTTaskMC&);
  AliEbyEPidTTaskMC& operator = (const AliEbyEPidTTaskMC&);
  ClassDef(AliEbyEPidTTaskMC, 1);

};

#endif

 
