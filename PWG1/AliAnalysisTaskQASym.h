#ifndef ALIANALYSISTASKQASYM_H
#define ALIANALYSISTASKQASYM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id:$ */
 
//------------------------------
// Analysis task for quality-assurance of central tracking
// mainly based on fundamental symmetries 
//
// eva.sicking@cern.ch

class TH1F;
class TH2F;
class TList;
class TNtuple;

class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;


#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskQASym : public AliAnalysisTaskSE {
 public:
    AliAnalysisTaskQASym();
    AliAnalysisTaskQASym(const char *name);
    virtual ~AliAnalysisTaskQASym() {}
  
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    
    
  
    virtual void   SetCuts(AliESDtrackCuts* cuts)
	{fCuts = cuts;}
    
    virtual void   SetFieldOn(Bool_t b = kTRUE){fFieldOn = b;} 
  
 private:

 
  Bool_t      fFieldOn;         // field flag

  TList       *fHists;          // List of histos

  //old
  TH1F        *fHistRECpt;      // pt 
  TH1F        *fEta;            // eta
  TH2F        *fEtaPhi;         // eta-phi
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
  TH1F        *fRecQPtPosEta;  // Q x Pt for pos. eta
  TH1F        *fRecQPtNegEta;  // Q x Pt for neg. eta
  TH1F        *fRecPtPosEta;   //     Pt for pos. eta
  TH1F        *fRecPtNegEta;   //     Pt for neg. eta
  TH1F        *fRecPhiPosEta;  // phi for pos. eta
  TH1F        *fRecPhiNegEta;  // phi for neg. eta 
  TH1F        *fRecDcaPosEta;  // dca for pos. eta 
  TH1F        *fRecDcaNegEta;  // dca for neg. eta
  TH1F        *fRecDPosEta;    // d   for pos. eta
  TH1F        *fRecDNegEta;    // d   for neg. eta

  // sectors of TPC (with pt>xGeV?), TODO: extent to TPC standalone tracks
  TH1F        *fRecPtTpcSector[18];      // pt per sector
  TH1F        *fRecEtaTpcSector[18];     // eta per sector
  TH1F        *fRecQPtTpcSector[18];     // Qxpt per sector
  TH1F        *fRecEtaPtTpcSector[18];   // eta x pt per sector
  TH1F        *fSignedDcaTpcSector[18];  // dca per sector

  // 7 different case of hit in ITS ladders
  TH1F        *fSignDcaPos[7];           // dca for pos. charge       
  TH1F        *fSignDcaNeg[7];           // dca for neg. charge
  TH1F        *fSignDcaNegInv[7];        // dca for neg. charge
  TH1F        *fPtSigmaPos[7];           // sigma pt for pos. charge
  TH1F        *fPtSigmaNeg[7];           // sigma pt for neg. charge
  TH1F        *fRecPtPosLadder[7];       // pt for pos. charge
  TH1F        *fRecPtNegLadder[7];       // pt for neg. charge
  TH1F        *fRecPhiPosLadder[7];      // phi for pos. charge
  TH1F        *fRecPhiNegLadder[7];      // phi for neg. charge
  TH1F        *fRecEtaPosLadder[7];      // eta for pos. charge
  TH1F        *fRecEtaNegLadder[7];      // eta for neg. charge

  // 2D: all measures as function of z of vertex
  TH2F        *fRecPtPosVz;              // pt for pos. charge
  TH2F        *fRecPtNegVz;              // pt for neg. charge
  TH2F        *fRecEtaPosVz;             // eta for pos. charge
  TH2F        *fRecEtaNegVz;             // eta for neg. charge
  TH2F        *fRecPhiPosVz;             // phi for pos. charge
  TH2F        *fRecPhiNegVz;             // phi for neg. charge
  TH2F        *fSignedDcaPosVz;          // dca for pos. charge
  TH2F        *fSignedDcaNegVz;          // dca for neg. charge
  TH2F        *fRecQPtPosEtaVz;          // qxpt for pos. charge
  TH2F        *fRecQPtNegEtaVz;          // qxpt for neg. charge
  TH2F        *fRecEtaPtPosVz;           // etaxpt for pos. charge
  TH2F        *fRecEtaPtNegVz;           // etaxpt for neg. charge

  //high
  TH1F * fDeltaPhiAll;                   // dphi
  TH2F * fDeltaPhiLeading;               // dphi rel. to leading
  TH1F * fDiffDcaD;                      // delta dca

  //sim
  TH1F * fPhiRec;                        // phi
  TH1F * fThetaRec;                      // theta
  TH1F * fNumber;                        // n
  TH1F * fVx;                            // vx
  TH1F * fVy;                            // vy
  TH1F * fVz;                            // vz

  AliESDtrackCuts* fCuts;                // List of cuts
  AliAnalysisTaskQASym(const AliAnalysisTaskQASym&); // not implemented
  AliAnalysisTaskQASym& operator=(const AliAnalysisTaskQASym&); // not implemented
  
  ClassDef(AliAnalysisTaskQASym, 1); // Basic QA exploiting symmetries
};

#endif
