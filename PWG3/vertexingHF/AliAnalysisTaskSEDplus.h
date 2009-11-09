#ifndef ALIANALYSISTASKDPLUS_H
#define ALIANALYSISTASKDPLUS_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskSEDplus
// AliAnalysisTaskSE for the D+ candidates Invariant Mass Histogram and 
//comparison of heavy-flavour decay candidates
// to MC truth (kinematics stored in the AOD)
// Renu Bala, bala@to.infn.it
//*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TH1F.h>

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisVertexingHF.h"

class AliAnalysisTaskSEDplus : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSEDplus();
  AliAnalysisTaskSEDplus(const char *name, Bool_t fillNtuple=kFALSE);
  virtual ~AliAnalysisTaskSEDplus();
  void SetReadMC(Bool_t readMC=kTRUE){fReadMC=readMC;}
  void SetDoLikeSign(Bool_t dols=kTRUE){fDoLS=dols;}
  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  void LSAnalysis(TClonesArray *arrayOppositeSign,TClonesArray *arrayLikeSign,AliAODEvent *aod,AliAODVertex *vtx1, Int_t nDplusOS);
    
 private:

  AliAnalysisTaskSEDplus(const AliAnalysisTaskSEDplus &source);
  AliAnalysisTaskSEDplus& operator=(const AliAnalysisTaskSEDplus& source); 
  TList   *fOutput; //! list send on output slot 0
  TH1F    *fHistNEvents; //!hist. for No. of events
  TNtuple *fNtupleDplus; //! output ntuple
  TH1F *fHistLS;     //! LS tripl normalized using all charged tracks in the event
  TH1F *fHistLStrip; //! LS tripl normalized using OS vs LS number of triplets ratio
  TH1F *fHistLStripcut;//! LS tripl normalized using OS vs LS number of triplets passing dplus cuts ratio
  TH1F *fHistOS; //! OS triplets
  TH1F *fHistOSbkg;//! only Background (need MC truth)
  TH1F *fHistLSnoweight;//! LS triplets wo normalization
  TH1F *fHistLSsinglecut;//! LS tripl normalized using tracks with Pt>0.4GeV
  
  Bool_t fFillNtuple;   // flag for filling ntuple
  Bool_t fReadMC;       //flag for access to MC
  Bool_t fDoLS;        //flag for switching on/off Like Sign analysis
  AliAnalysisVertexingHF *fVHF;  // Vertexer heavy flavour (used to pass the cuts)
  
  ClassDef(AliAnalysisTaskSEDplus,4); // AliAnalysisTaskSE for the MC association of heavy-flavour decay candidates
};

#endif

