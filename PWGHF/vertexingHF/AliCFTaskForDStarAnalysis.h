#ifndef ALICFTASKFORDSTARANALYSIS_H
#define ALICFTASKFORDSTARANALYSIS_H
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
// Class for D* corrections --  

#include "AliAnalysisTaskSE.h"

class TH1I;
class TParticle ;
class TFile ;
class TClonesArray ;
class AliCFManager;
class AliAODRecoDecay;
class AliAODRecoDecayHF2Prong;
class AliAODMCParticle;
class THnSparse;

class AliCFTaskForDStarAnalysis : public AliAnalysisTaskSE {
  public:

  enum {
    kStepGenerated       = 0,
    kStepAcceptance      = 1,
    kStepVertex          = 2,
    kStepRefit           = 3,
    kStepReconstructed   = 4,
    kStepRecoAcceptance  = 5,
    kStepRecoITSClusters = 6,
    kStepRecoCuts        = 7
  };

  AliCFTaskForDStarAnalysis();
  AliCFTaskForDStarAnalysis(const Char_t* name);
  AliCFTaskForDStarAnalysis& operator= (const AliCFTaskForDStarAnalysis& c);
  AliCFTaskForDStarAnalysis(const AliCFTaskForDStarAnalysis& c);
  virtual ~AliCFTaskForDStarAnalysis();

  // ANALYSIS FRAMEWORK STUFF to loop on data and fill output objects
  void     UserCreateOutputObjects();
  void     UserExec(Option_t *option);
  void     Terminate(Option_t *);

 // UNFOLDING
  void     SetCorrelationMatrix(THnSparse* h) {fCorrelation=h;}
  void     SetAcceptanceUnf(Bool_t AcceptanceUnf) {fAcceptanceUnf = AcceptanceUnf;}
  Bool_t   GetAcceptanceUnf() const {return fAcceptanceUnf;}

  // CORRECTION FRAMEWORK
  void           SetCFManager(AliCFManager* io) {fCFManager = io;}   // global correction manager
  AliCFManager * GetCFManager()                 {return fCFManager;} // get corr manager

  Bool_t   GetDStarMCParticle(AliAODMCParticle* mcPart, TClonesArray* mcArray, Double_t* vectorMC)const;
  Bool_t   EvaluateIfD0toKpi(AliAODMCParticle* neutralDaugh, TClonesArray* mcArray, Double_t* VectorD0)const;
  // for the D0 
  void     SetMinITSClusters(Int_t minITSClusters) {fMinITSClusters = minITSClusters;}
  Int_t    GetMinITSClusters() const {return fMinITSClusters;}
  // for the soft pion
  void     SetMinITSClustersSoft(Int_t minITSClustersSoft) {fMinITSClustersSoft = minITSClustersSoft;}
  Int_t    GetMinITSClustersSoft() const {return fMinITSClustersSoft;}

 protected:
 
  AliCFManager   *fCFManager;   //  pointer to the CF manager
  TH1I *fHistEventsProcessed;   //! simple histo for monitoring the number of events processed
  THnSparse* fCorrelation;      //  response matrix for unfolding
  Int_t fCountRecoDStarSel;     //  Reco particle found that satisfy cuts in D* selection
  Int_t fEvents;                //  n. of events
  Int_t fMinITSClusters;        //  min n. of ITS clusters for RecoDecay
  Int_t fMinITSClustersSoft;    //  min n. of ITS clusters for RecoDecay soft pion
  Bool_t fAcceptanceUnf;        //  flag for unfolding before or after cuts.
  
  ClassDef(AliCFTaskForDStarAnalysis,3); // class for D* corrections
};

#endif
