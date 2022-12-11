/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskPsi2Spolarization_H
#define AliAnalysisTaskPsi2Spolarization_H
#include "AliAnalysisTaskSE.h"



class AliAnalysisTaskPsi2Spolarization : public AliAnalysisTaskSE  
{
 public:
  // two class constructors
  AliAnalysisTaskPsi2Spolarization();
  AliAnalysisTaskPsi2Spolarization(const char *name);
  // class destructor
  virtual                 ~AliAnalysisTaskPsi2Spolarization();
  // called once at beginning of runtime

  virtual void            NotifyRun();
  virtual void            UserCreateOutputObjects();
  // called for each event
  virtual void            UserExec(Option_t* option);
  // called at end of analysis
  virtual void            Terminate(Option_t* option);
  AliMuonTrackCuts* GetTrackCuts() { return fMuonTrackCuts; }
    
    void SetBeamEnergy(Double_t en) {fBeamEnergy=en;}
 
 private:
    
  AliMuonTrackCuts* fMuonTrackCuts;
  Double_t fBeamEnergy;   // Energy of the beam (required for the CS angle)
  AliAODEvent*            fAOD;           //!<! input event
  TList*                  fOutputList;    //!<! output list
  TH1F*                   fHistPtDimu;        //!<! dummy histogram
  TH1F*                   fHistRap;        //!<! dummy histogram
  TH1F*                   fHistInvMass;
  TH1F*                   fVertexZ;
  TH1F*                   fHistCosThetaCSDimu;
  TH1F*                   fHistCosThetaHEDimu;
  TH1F*                   fHistPhiCSDimu;
  TH1F*                   fHistPhiHEDimu;
  TTree*                  tree;
  TTree*                  tree1;
  //AliAnalysisTaskPsi2Spolarization(const AliAnalysisTaskPsi2Spolarization&); // not implemented
  Float_t     Rap(Float_t e, Float_t pz) const;
  //AliAnalysisTaskPsi2Spolarization& operator=(const AliAnalysisTaskPsi2Spolarization&); // not implemented
    
  Double_t CostCS(AliAODTrack* , AliAODTrack* );
  Double_t CostHE(AliAODTrack* , AliAODTrack* );
  Double_t PhiCS(AliAODTrack* , AliAODTrack* );
  Double_t PhiHE(AliAODTrack* , AliAODTrack* );
    
  ClassDef(AliAnalysisTaskPsi2Spolarization, 1);
};

#endif
