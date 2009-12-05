#ifndef AliAnalysisTaskQASym_cxx
#define AliAnalysisTaskQASym_cxx
 

class TH1F;
class TH2F;
class TList;
class TNtuple;

class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;


#include "AliAnalysisTaskSE.h"
#include "TFile.h"
#include "TNtuple.h"

class AliAnalysisTaskQASym : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskQASym(const char *name = "AliAnalysisTaskQASym");
  virtual ~AliAnalysisTaskQASym() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  
  
  virtual void   SetCuts(AliESDtrackCuts* cuts)
     {fCuts = cuts;}

  virtual void   SetFieldOn(Bool_t b = kTRUE){fFieldOn = b;} 
  
 private:

 
  Bool_t      fFieldOn;

  TList       *fHists;          // List of histos

  //old
  TH1F        *fHistRECpt;      // pt 
  TH1F        *fEta;            // eta
  TH1F        *fEtaPt;          // eta over pt 
  TH1F        *fQPt;            // charge over pt 
  TH1F        *fDca;            // distance of closest approach
  TH1F        *fqPtRec[7];      // charge over pt divided for ITS layer cases
  TH1F        *fqRec;           // reconstrcuted charge
  TH1F        *fsigmaPt;        // sigma_pT
  TH2F        *fDcaSigmaPos[7]; // distance of closest approach against sigma_pT
  TH2F        *fDcaSigmaNeg[7]; // distance of closest approach against sigma_pT


  //positive und negative particles
  TH1F        *fRecPtPos;      // pt of pos partice
  TH1F        *fRecPtNeg;      // pt of neg particle
  TH1F        *fRecPhiPos;     // phi of pos. particle
  TH1F        *fRecPhiNeg;     // phi of neg. particle
  TH1F        *fRecEtaPos;     // eta of neg. particle
  TH1F        *fRecEtaNeg;     // eta of neg. particle
  TH1F        *fRecEtaPtPos;   // eta over pt of neg. particle
  TH1F        *fRecEtaPtNeg;   // eta over pt of neg. particle
  TH1F        *fRecDcaPos;     // distance of closest approach of neg. particle
  TH1F        *fRecDcaNeg;     // distance of closest of neg. particle
  TH1F        *fRecDcaNegInv;  // invers dca of neg. particle
  TH1F        *fRecDPos;       // impact parameter of neg. particle
  TH1F        *fRecDNeg;       // impact parameter of neg. particle




  // two sides of TPC -> Eta/Theta
  TH1F        *fRecQPtPosEta; 
  TH1F        *fRecQPtNegEta;
  TH1F        *fRecPtPosEta;  
  TH1F        *fRecPtNegEta; 
  TH1F        *fRecPhiPosEta;
  TH1F        *fRecPhiNegEta;
  TH1F        *fRecDcaPosEta;
  TH1F        *fRecDcaNegEta;
  TH1F        *fRecDPosEta;
  TH1F        *fRecDNegEta;


 
  // sectors of TPC (with pt>xGeV?), TODO: extent to TPC standalone tracks
  TH1F        *fRecPtTpcSector[18];
  TH1F        *fRecEtaTpcSector[18];
  TH1F        *fRecQPtTpcSector[18];
  TH1F        *fRecEtaPtTpcSector[18];
  TH1F        *fSignedDcaTpcSector[18];


  // 7 different case of hit in ITS ladders
  TH1F        *fSignDcaPos[7];
  TH1F        *fSignDcaNeg[7];
  TH1F        *fSignDcaNegInv[7];
  TH1F        *fPtSigmaPos[7];
  TH1F        *fPtSigmaNeg[7];
  TH1F        *fRecPtPosLadder[7];
  TH1F        *fRecPtNegLadder[7];
  TH1F        *fRecPhiPosLadder[7];
  TH1F        *fRecPhiNegLadder[7];
  TH1F        *fRecEtaPosLadder[7];
  TH1F        *fRecEtaNegLadder[7];


  // 2D: all measures as function of z of vertex
  TH2F        *fRecPtPosVz;
  TH2F        *fRecPtNegVz;
  TH2F        *fRecEtaPosVz;
  TH2F        *fRecEtaNegVz;
  TH2F        *fRecPhiPosVz;
  TH2F        *fRecPhiNegVz;
  TH2F        *fSignedDcaPosVz;
  TH2F        *fSignedDcaNegVz;
  TH2F        *fRecQPtPosEtaVz;
  TH2F        *fRecQPtNegEtaVz;
  TH2F        *fRecEtaPtPosVz;
  TH2F        *fRecEtaPtNegVz;


  //high
  TH1F * fDeltaPhiAll;
  TH2F * fDeltaPhiLeading;
  TH1F * fDiffDcaD;

  //sim
  TH1F * fPhiRec;
  TH1F * fThetaRec;
  TH1F * fNumber;
  TH1F * fVx;
  TH1F * fVy;
  TH1F * fVz;
  TNtuple * test;



  Double_t  sdca_tr;
  Float_t xy, z;

  AliESDtrackCuts* fCuts;                      // List of cuts
  AliAnalysisTaskQASym(const AliAnalysisTaskQASym&); // not implemented
  AliAnalysisTaskQASym& operator=(const AliAnalysisTaskQASym&); // not implemented
  
  ClassDef(AliAnalysisTaskQASym, 1); // Basic QA exploiting symmetries
};

#endif
