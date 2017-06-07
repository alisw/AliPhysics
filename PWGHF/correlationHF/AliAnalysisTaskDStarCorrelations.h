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

/* $Id: AliAnalysisTaskDStarCorrelations.h 65139 2013-11-25 14:47:45Z fprino $ */

//-----------------------------------------------------------------------
//          
//
//						   Author S.Bjelogrlic
//                         Utrecht University 
//                      sandro.bjelogrlic@cern.ch                       
//                      Mandeep Kour->CorrelationVsMultiplicity
//                      mandeep.kour@cern.ch
//-----------------------------------------------------------------------
#include <TROOT.h>
#include <TSystem.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TProfile.h>
#include <TArrayD.h>
#include <TRandom.h>
#include <THnSparse.h>
#include "AliAnalysisTaskSE.h"
#include "AliAODEvent.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliEventPoolManager.h"
#include "AliAODRecoCascadeHF.h"
#include "AliHFAssociatedTrackCuts.h"
#include "AliNormalizationCounter.h"
#include "AliHFCorrelator.h"
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
  enum BkgMethod{kDZeroSB, kDStarSB};
  
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
    void SetCorrelator(Int_t l) {fselect = l;} // select 1 for hadrons, 2 for Kaons, 3 for Kzeros
    void SetMonteCarlo(Bool_t k) {fmontecarlo = k;}
    void SetUseMixing (Bool_t j) {fmixing = j;}
    void SetUseMult (Bool_t j) {fmult = j;}
    void SetUseFullMode (Bool_t j) {fFullmode = j;}
    void SetCollSys(CollSyst system){fSystem=system;} // select between pp (kFALSE) or PbPb (kTRUE)
    void SetEfficiencyVariable(DEffVariable var){fEfficiencyVariable = var;} // set the efficiency variable to use
    void SetBkgEstimationMethod(BkgMethod var){fBkgMethod = var;} // set the efficiency variable to use
    void SetUseReconstruction(Bool_t reco){fReco = reco;}
    void SetUseEfficiencyCorrection(Bool_t correction){fUseEfficiencyCorrection = correction;} // setter for using the single track efficiency correction
    void SetUseDmesonEfficiencyCorrection(Bool_t correction){fUseDmesonEfficiencyCorrection = correction;} // setter for using the single track efficiency correction
    void SetUseCentrality (Bool_t j) {fUseCentrality = j;} // switch for centrality (kTRUE)/multiplicity(kFALSE)
    void SetUseHadronicChannelAtKineLevel (Bool_t use){fUseHadronicChannelAtKineLevel = use;}
    void SetUseRemoveMoreThanOneCDmesonCandidate (Bool_t use){fRemoveMoreThanOneDmesonCandidate = use;}
    void SetLimitAcceptanceForMC (Bool_t use){fLimitAcceptanceForMC = use;}
    void SetAODMismatchProtection(Int_t opt=1) {fAODProtection=opt;}   
    void SetUseSmallSizePlots(Bool_t smsize){fUseSmallSizePlots=smsize;}
    
    void SetNofPhiBins(Int_t nbins){fPhiBins = nbins;} // number of delta phi bins
    void SetLevelOfDebug(Int_t debug){fDebugLevel=debug;} // set debug level
    void SetUseDisplacement(Int_t m) {fDisplacement=m;} // select 0 for no displ, 1 for abs displ, 2 for d0/sigma_d0
  
    
    void SetDim(){fDim = 4;
        fDMesonSigmas = new Float_t[4];} // standard definedt = cannot be changed from outside
    
    void SetMaxDStarEta(Double_t eta){fMaxEtaDStar = eta;} // maximum eta D*

  
  void SetDMesonSigmas(Float_t DStarWin, Float_t D0Win, Float_t SBmin, Float_t SBmax){
    fDMesonSigmas[0] = DStarWin;
    fDMesonSigmas[1] = D0Win;
    fDMesonSigmas[2] = SBmin;
    fDMesonSigmas[3] = SBmax;
    
  }
void SetUseMultarray(Int_t nMultBin, Double_t *MultBinLimits,Double_t minMult,Double_t maxMult)
  { fmultarray=new  Double_t[nMultBin+1];
    fnmultBins=nMultBin;
    for(Int_t i=0;i<fnmultBins+1;i++)
      {
    fmultarray[i]=MultBinLimits[i];
      }
    
    fminMult= minMult;
    fmaxMult= maxMult;
}

void SetMultiplVsZProfileLHC10b(TProfile* hprof){
    if(fMultEstimatorAvg[0]) delete fMultEstimatorAvg[0];
    fMultEstimatorAvg[0]=new TProfile(*hprof);
  }
    void SetMultiplVsZProfileLHC10c(TProfile* hprof){
    if(fMultEstimatorAvg[1]) delete fMultEstimatorAvg[1];
    fMultEstimatorAvg[1]=new TProfile(*hprof);
  }
    void SetMultiplVsZProfileLHC10d(TProfile* hprof){
    if(fMultEstimatorAvg[2]) delete fMultEstimatorAvg[2];
    fMultEstimatorAvg[2]=new TProfile(*hprof);
    }
void SetMultiplVsZProfileLHC10e(TProfile* hprof){
    if(fMultEstimatorAvg[3]) delete fMultEstimatorAvg[3];
    fMultEstimatorAvg[3]=new TProfile(*hprof);
    }

  enum { kNtrk10=0, kNtrk10to16=1, kVZERO=2, kNtrk03=3, kNtrk05=4, kVZEROA=5, kVZEROEq=6, kVZEROAEq=7 };
  void SetMultiplicityEstimator(Int_t value){ fMultiplicityEstimator=value; }
  Int_t GetMultiplicityEstimator(){ return fMultiplicityEstimator; }
  TProfile* GetEstimatorHistogram(const AliVEvent *event);
  void SetReferenceMultiplcity(Double_t rmu){fRefMult=rmu;}
  enum { kEta10=0, kEta10to16=1, kEtaVZERO=2, kEta03=3, kEta05=5, kEtaVZEROA=5};
 
  void SetDeffMapvsPt(TH1D * map){fDeffMapvsPt = map;}
  void SetDeffMapvsPtvsMult(TH2D * map){fDeffMapvsPtvsMult = (TH2D*)map;}
  void SetDeffMapvsPtvsMultvsEta(TH2D * map){fDeffMapvsPtvsEta = map;}
  
  void SetUseMCEventType(Bool_t k){fMCEventType = k;}  
  

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
  Bool_t fmult;// switch for multiplicity analysis
  Int_t fnmultBins;
  Double_t* fmultarray;//[fnmultBins]
  Double_t fminMult;
  Double_t fmaxMult;
  Bool_t fFullmode;
  CollSyst fSystem; // pp, pPb or PbPb
  DEffVariable  fEfficiencyVariable; // set second variable to study efficiency (mult, centr, y, eta)
  BkgMethod fBkgMethod; // bkg estimation method (dstar or dzero sidebands)
  Bool_t fReco; // use reconstruction or MC truth
  Bool_t fUseEfficiencyCorrection; // boolean variable to use or not the efficiency correction
  Bool_t fUseDmesonEfficiencyCorrection; // boolean flag for the use of Dmeson efficiency correction
  Bool_t fUseCentrality;// boolean to switch in between centrality or multiplicity
  Bool_t fUseHadronicChannelAtKineLevel; //
  Bool_t fRemoveMoreThanOneDmesonCandidate; // flag to remove a second, 3rd etc candidate if there is any - useful in PbPb
  Bool_t fLimitAcceptanceForMC; // flag to remove a second, 3rd etc candidate if there is any - useful in PbPb
  Bool_t fUseSmallSizePlots; //flag to reduce number o bins in THnSparse (for merging issues)

  Int_t fMultiplicityEstimator; // Definition of the multiplicity estimator: kNtrk10=0, kNtrk10to16=1, kVZERO=2
  Int_t fDoVZER0ParamVertexCorr; 

  Int_t fPhiBins;
  Int_t fEvents; //! number of event
  Int_t fDebugLevel; //! debug level
  Int_t fDisplacement; // set 0 for no displacement cut, 1 for absolute d0, 2 for d0/sigma_d0
  Int_t fDim;//
  Int_t fNofPtBins;
  Double_t fMaxEtaDStar;
  Float_t *fDMesonSigmas;//[fDim]
  Float_t * fD0Window;  //[fNofPtBins]
   
  Bool_t fMCEventType; // Use MC event type 
  
  Double_t fRefMult;   // refrence multiplcity (period b)
  Int_t fAODProtection;            // flag to activate protection against AOD-dAOD mismatch.
  
  TList *fOutput;                  //! user output data
  TList *fDmesonOutput; //!output related to d meson
  TList *fTracksOutput; //!output related to tracks
  TList *fEMOutput; //! output with EM checks
  TList *fCorrelationOutput; // ! output with correlation sparses
  TList *fOutputMC;                //! outpu for MC
  AliRDHFCutsDStartoKpipi *fCuts;  // Cuts D*
  AliHFAssociatedTrackCuts *fAssocCuts; // cuts for associated
  AliAnalysisUtils *fUtils;
  AliAODTracklets * fTracklets; // AliAODtracklets
  
  TH1D * fDeffMapvsPt; // histo for Deff mappin
  TH2D * fDeffMapvsPtvsMult; // histo for Deff mappin
  TH2D * fDeffMapvsPtvsEta; // histo for Deff mappin
  TProfile* fMultEstimatorAvg[4]; 
  ClassDef(AliAnalysisTaskDStarCorrelations,11); // class for D meson correlations
  
};

#endif
