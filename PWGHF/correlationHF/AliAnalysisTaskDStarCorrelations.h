#ifndef AliAnalysisTaskDStarCorrelations_H
#define AliAnalysisTaskDStarCorrelations_H

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

//-----------------------------------------------------------------------
//          
//
//						   Author S.Bjelogrlic
//                         Utrecht University 
//                      sandro.bjelogrlic@cern.ch
//
//-----------------------------------------------------------------------

#include <TH2F.h>
#include <TH1D.h>
#include <TH3D.h>
#include "AliAnalysisTaskSE.h"
#include "AliAODEvent.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliEventPoolManager.h"
#include "AliAODRecoCascadeHF.h"
#include "AliHFAssociatedTrackCuts.h"
#include "AliNormalizationCounter.h"
#include "AliHFCorrelator.h"
#include <THnSparse.h>
#include "AliAnalysisUtils.h"
#include "AliVertexingHFUtils.h"
class TParticle ;
class TClonesArray ;
class AliAODMCParticle;
class AliAODEvent;
class AliVParticle;
class TObjArray;
class AliEventPoolManager;
class AliESDEvent;



class AliAnalysisTaskDStarCorrelations : public AliAnalysisTaskSE
{

 public :
  
  enum CollSyst {pp,pA,AA};
  enum DEffVariable{kNone,kMult,kCentr,kRapidity,kEta};
  
  AliAnalysisTaskDStarCorrelations();
  AliAnalysisTaskDStarCorrelations(const Char_t* name,AliRDHFCutsDStartoKpipi* cuts, AliHFAssociatedTrackCuts *AsscCuts, AliAnalysisTaskDStarCorrelations::CollSyst syst,Bool_t mode);
  virtual ~AliAnalysisTaskDStarCorrelations();
  
  
  
  // Implementation of interface methods  
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  
  void DefineThNSparseForAnalysis();
  void DefineHistoForAnalysis();
  void EnlargeDZeroMassWindow();
    Bool_t IsDDaughter(AliAODMCParticle* d, AliAODMCParticle* track) const ;
  
  // checker for event mixing
  void EventMixingChecks(AliAODEvent * AOD); 
  // setters
  void SetMonteCarlo(Bool_t k) {fmontecarlo = k;}
  void SetUseMixing (Bool_t j) {fmixing = j;}
  void SetCorrelator(Int_t l) {fselect = l;} // select 1 for hadrons, 2 for Kaons, 3 for Kzeros
  void SetUseDisplacement(Int_t m) {fDisplacement=m;} // select 0 for no displ, 1 for abs displ, 2 for d0/sigma_d0
  void SetCollSys(CollSyst system){fSystem=system;} // select between pp (kFALSE) or PbPb (kTRUE)
  void SetEfficiencyVariable(DEffVariable var){fEfficiencyVariable = var;} // set the efficiency variable to use
  void SetLevelOfDebug(Int_t debug){fDebugLevel=debug;} // set debug level
  void SetUseReconstruction(Bool_t reco){fReco = reco;}
  void SetDMesonSigmas(Float_t DStarWin, Float_t D0Win, Float_t SBmin, Float_t SBmax){
    fDMesonSigmas[0] = DStarWin;
    fDMesonSigmas[1] = D0Win;
    fDMesonSigmas[2] = SBmin;
    fDMesonSigmas[3] = SBmax;
    
  }

  void SetUseEfficiencyCorrection(Bool_t correction){fUseEfficiencyCorrection = correction;} // setter for using the single track efficiency correction
  void SetUseDmesonEfficiencyCorrection(Bool_t correction){fUseDmesonEfficiencyCorrection = correction;} // setter for using the single track efficiency correction
  void SetUseHadronicChannelAtKineLevel (Bool_t use){fUseHadronicChannelAtKineLevel = use;}
    
  void SetDim(){fDim = 4;
    fDMesonSigmas = new Float_t[4];}
  void SetDeffMapvsPt(TH1D * map){fDeffMapvsPt = map;}
  void SetDeffMapvsPtvsMult(TH2D * map){fDeffMapvsPtvsMult = (TH2D*)map;}
  void SetDeffMapvsPtvsMultvsEta(TH2D * map){fDeffMapvsPtvsEta = map;}
  void SetNofPhiBins(Int_t nbins){fPhiBins = nbins;}
    void SetMaxDStarEta(Double_t eta){fMaxEtaDStar = eta;}
  
  

private:
  
  AliAnalysisTaskDStarCorrelations(const AliAnalysisTaskDStarCorrelations &source);
  AliAnalysisTaskDStarCorrelations& operator=(const AliAnalysisTaskDStarCorrelations& source);
  
  TObject* fhandler; //! Analysis Handler
  TClonesArray* fmcArray; //mcarray
  AliNormalizationCounter *fCounter; // counter
  AliHFCorrelator * fCorrelator; // object for correlations
  
  
  
  Int_t fselect; // select what to correlate with a D* 1-chargedtracks,2-chargedkaons,3-k0s
  Bool_t fmontecarlo;//switch for MC
  Bool_t fmixing;// switch for event mixing
  Bool_t fFullmode;
  CollSyst fSystem; // pp, pPb or PbPb
  DEffVariable  fEfficiencyVariable; // set second variable to study efficiency (mult, centr, y, eta)
  Bool_t fReco; // use reconstruction or MC truth
  Bool_t fUseEfficiencyCorrection; // boolean variable to use or not the efficiency correction
  Bool_t fUseDmesonEfficiencyCorrection; // boolean flag for the use of Dmeson efficiency correction
  Bool_t fUseCentrality;// boolean to switch in between centrality or multiplicity
  Bool_t fUseHadronicChannelAtKineLevel; //
 
  Int_t fPhiBins;
  Int_t fEvents; //! number of event
  Int_t fDebugLevel; //! debug level
  Int_t fDisplacement; // set 0 for no displacement cut, 1 for absolute d0, 2 for d0/sigma_d0
  Int_t fDim;//
  Int_t fNofPtBins;
     Double_t fMaxEtaDStar;
  Float_t *fDMesonSigmas;//[fDim]
  Float_t * fD0Window;  //[fNofPtBins]
   
  
  
  
  TList *fOutput;                  //! user output data
  TList *fOutputMC;                //! outpu for MC
  AliRDHFCutsDStartoKpipi *fCuts;  // Cuts D*
  AliHFAssociatedTrackCuts *fCuts2; // cuts for associated
  AliAnalysisUtils *fUtils;
  AliAODTracklets * fTracklets; // AliAODtracklets
  
  TH1D * fDeffMapvsPt; // histo for Deff mappin
  TH2D * fDeffMapvsPtvsMult; // histo for Deff mappin
  TH2D * fDeffMapvsPtvsEta; // histo for Deff mappin
  
  ClassDef(AliAnalysisTaskDStarCorrelations,6); // class for D meson correlations
  
};

#endif
