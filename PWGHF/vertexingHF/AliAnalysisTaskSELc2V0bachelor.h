#ifndef ALIANALYSISTASKSELC2V0BACHELOR_H
#define ALIANALYSISTASKSELC2V0BACHELOR_H
/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */ 

#include <TH2F.h>
#include "TROOT.h"
#include "TSystem.h"

#include "AliAnalysisTaskSE.h"
#include "AliAODEvent.h"
#include "AliPID.h"
#include "AliAODTrack.h"
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"
#include "AliRDHFCutsLctoV0.h"
#include "AliNormalizationCounter.h"

class AliAnalysisTaskSELc2V0bachelor : public AliAnalysisTaskSE 
{
  
 public:
  
  AliAnalysisTaskSELc2V0bachelor();
  AliAnalysisTaskSELc2V0bachelor(const Char_t* name, AliRDHFCutsLctoV0* cutsA, AliRDHFCutsLctoV0* cutsB,
				 Bool_t useOnTheFly=kFALSE);
  virtual ~AliAnalysisTaskSELc2V0bachelor();

  // Implementation of interface methods  
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
 
  // histos
  void FillLc2pK0Sspectrum(AliAODRecoCascadeHF *part, Int_t isLc,
			   Int_t &nSelectedProd, AliRDHFCutsLctoV0 *cutsProd,
			   Int_t &nSelectedAnal, AliRDHFCutsLctoV0 *cutsAnal);

  void DefineHistograms();
  Int_t CheckOrigin(TClonesArray* arrayMC, AliAODMCParticle *mcPartCandidate) const;

  void MakeAnalysisForLc2prK0S(AliAODVertex * vtx1,
			       TClonesArray *arrayLctopK0s,
			       TClonesArray *mcArray,
			       Int_t &nSelectedProd, AliRDHFCutsLctoV0 *cutsProd,
			       Int_t &nSelectedAnal, AliRDHFCutsLctoV0 *cutsAnal);
 
  Int_t MatchToMC(AliAODRecoCascadeHF *lc2bacV0,
		  Int_t *pdgDgLc2bacV0, Int_t *pdgDgV0,
		  TClonesArray *mcArray);

  void SearchLcDaughter(TClonesArray *arrayMC);

  void DefineK0SHistos();

  // set MC usage
  void SetMC(Bool_t theMCon) {fUseMCInfo = theMCon;}
  Bool_t GetMC() const {return fUseMCInfo;}

  void FillArmPodDistribution(AliAODv0 *vZero, TString histoTitle, TList *histoList);

  void SetK0sAnalysis(Bool_t a) {fIsK0sAnalysis=a;}
  Bool_t GetK0sAnalysis() const {return fIsK0sAnalysis;}

  void SetUseOnTheFlyV0(Bool_t a) { fUseOnTheFlyV0=a; }
  Bool_t GetUseOnTheFlyV0() { return fUseOnTheFlyV0; }

 private:
  
  AliAnalysisTaskSELc2V0bachelor(const AliAnalysisTaskSELc2V0bachelor &source);
  AliAnalysisTaskSELc2V0bachelor& operator=(const AliAnalysisTaskSELc2V0bachelor& source); 
  
  Bool_t fUseMCInfo;          // Use MC info
  TList *fOutput;             // User output1 // general histos
  TList *fOutputAll;          // User output2 // histos without pid and cut on V0
  TList *fOutputPIDBach;      // User output3 // histos with PId on Bachelor

  // define the histograms
  TH1I *fCEvents;                    // Histogram to check selected events
  AliPIDResponse *fPIDResponse;      //! PID response object
  Bool_t fIsK0sAnalysis;             // switch between Lpi and K0sp
  AliNormalizationCounter *fCounter; // AliNormalizationCounter on output slot 4
  AliRDHFCutsLctoV0 *fProdCuts;      // Cuts - sent to output slot 5
  AliRDHFCutsLctoV0 *fAnalCuts;      // Cuts - sent to output slot 5
  TList *fListCuts;                  // list of cuts
  Bool_t fUseOnTheFlyV0;             // flag to analyze also on-the-fly V0 candidates
  Bool_t fIsEventSelected;           // flag for event selected

  ClassDef(AliAnalysisTaskSELc2V0bachelor,2); // class for Lc->p K0
};

#endif

