#ifndef ALIANALYSISTASKSEDSTARJETS_H
#define ALIANALYSISTASKSEDSTARJETS_H
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

//-----------------------------------------------------------------------
// Author : A. Grelli, UTRECHT
//-----------------------------------------------------------------------


#include <TH2F.h>
#include "AliAnalysisTaskSE.h"
#include "AliAODEvent.h"
#include "AliRDHFCutsDStartoKpipi.h"

class TParticle ;
class TClonesArray ;
class AliAODMCParticle;


class AliAnalysisTaskSEDStarJets : public AliAnalysisTaskSE 
{
  
 public:
  
  AliAnalysisTaskSEDStarJets();
  AliAnalysisTaskSEDStarJets(const Char_t* name,AliRDHFCutsDStartoKpipi* cuts);
  virtual ~AliAnalysisTaskSEDStarJets();
  
  virtual void     UserCreateOutputObjects();
  virtual void     UserExec(Option_t *option);
  virtual void     Terminate(Option_t *);
  virtual void     Init();
  virtual void     LocalInit() {Init();}

  //side band background eval
  void     SideBandBackground(Double_t finvM, Double_t finvMDStar,  Double_t dStarMomBkg, Double_t fejet, Double_t ejet);
  
  // inizializations
  Bool_t   DefineHistoFroAnalysis();  
  //MC FF
  double   FillMCFF(AliAODMCParticle* mcPart, TClonesArray* mcArray, Int_t mcLabel);
  // correction for UA1 cone algorithm
  void     SetChargeFractionCorrection(Int_t chargeFrCorr) {fchargeFrCorr =  chargeFrCorr;}
  Int_t    GetChargeFractionCorrection() const {return fchargeFrCorr;}

  // set MC usage
  void    SetMC(Bool_t theMCon) {fUseMCInfo = theMCon;}
  Bool_t  GetMC() const {return fUseMCInfo;}
  
 private:
  
  AliAnalysisTaskSEDStarJets(const AliAnalysisTaskSEDStarJets &source);
  AliAnalysisTaskSEDStarJets& operator=(const AliAnalysisTaskSEDStarJets& source); 

  Int_t  fEvents;                //  n. of events
  Int_t  fchargeFrCorr;          //  Charge fraction correction UA1 algorithm
  Bool_t fUseMCInfo;             //  Use MC info
  Bool_t fRequireNormalization;  //  normalization 
  
  TList *fOutput;                  //! user output
  AliRDHFCutsDStartoKpipi *fCuts;  // Cuts 

  // define the histograms 

  TH1F *ftrigger;        //!
  TH1F *fPtPion;         //!
  TH1F *fInvMass;        //!
  TH1F *fRECOPtDStar;    //!
  TH1F *fRECOPtBkg;      //!
  TH1F *fDStar;          //!
  TH1F *fDiff;           //!
  TH1F *fDiffSideBand;   //!
  TH1F *fDStarMass;      //!
  TH1F *fPhi;            //!
  TH1F *fPhiBkg;         //!
  TH1F *fTrueDiff;       //!
  TH1F *fResZ;           //!
  TH1F *fResZBkg;        //!  
  TH1F *fEjet;           //!
  TH1F *fPhijet;         //!
  TH1F *fEtaJet;         //!
  TH1F *theMCFF;         //!
  TH1F *fDphiD0Dstar;    //!
  TH1F *fPtJet;          //!

  ClassDef(AliAnalysisTaskSEDStarJets,3); // class for charm-jet correlations
};

#endif
