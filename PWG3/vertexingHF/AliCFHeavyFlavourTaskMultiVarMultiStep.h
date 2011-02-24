#ifndef ALICFHEAVYFLAVOURTASKMULTIVARMULTISTEP_H
#define ALICFHEAVYFLAVOURTASKMULTIVARMULTISTEP_H
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
// Class for HF corrections as a function of many variables and step 
// Author : C. Zampolli, CERN
// Base class for HF Unfolding - agrelli@uu.nl
//-----------------------------------------------------------------------


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
class AliRDHFCutsD0toKpi;

class AliCFHeavyFlavourTaskMultiVarMultiStep : public AliAnalysisTaskSE {
  public:

  enum {
    kStepGeneratedLimAcc = 0,
    kStepGenerated       = 1,
    kStepAcceptance      = 2,
    kStepVertex          = 3,
    kStepRefit           = 4,
    kStepReconstructed   = 5,
    kStepRecoAcceptance  = 6,
    kStepRecoITSClusters = 7,
    kStepRecoPPR         = 8,
    kStepRecoPID         = 9
  };

  AliCFHeavyFlavourTaskMultiVarMultiStep();
  AliCFHeavyFlavourTaskMultiVarMultiStep(const Char_t* name, AliRDHFCutsD0toKpi* cuts);
  AliCFHeavyFlavourTaskMultiVarMultiStep& operator= (const AliCFHeavyFlavourTaskMultiVarMultiStep& c);
  AliCFHeavyFlavourTaskMultiVarMultiStep(const AliCFHeavyFlavourTaskMultiVarMultiStep& c);
  virtual ~AliCFHeavyFlavourTaskMultiVarMultiStep();

  // ANALYSIS FRAMEWORK STUFF to loop on data and fill output objects
  void     UserCreateOutputObjects();
  void     UserExec(Option_t *option);
  void     Init();
  void     LocalInit() {Init();}
  void     Terminate(Option_t *);

 // UNFOLDING
  void     SetCorrelationMatrix(THnSparse* h) {fCorrelation=h;}
  void     SetAcceptanceUnf(Bool_t AcceptanceUnf) {fAcceptanceUnf = AcceptanceUnf;}
  Bool_t   GetAcceptanceUnf() const {return fAcceptanceUnf;}

  
  // CORRECTION FRAMEWORK RELATED FUNCTIONS
  void           SetCFManager(AliCFManager* io) {fCFManager = io;}   // global correction manager
  AliCFManager * GetCFManager()                 {return fCFManager;} // get corr manager

  void     SetPDG(Int_t code) {fPDG = code; }     // defines the PDG code of searched HF
  Double_t CosThetaStar(AliAODMCParticle* mcPart, AliAODMCParticle* mcPartDaughter0, AliAODMCParticle* mcPartDaughter1) const;  // returns cos(ThetaStar) of the D0 decay
  Double_t CT(AliAODMCParticle* mcPart, AliAODMCParticle* mcPartDaughter0, AliAODMCParticle* mcPartDaughter1) const;            // returns cT of the D0 decay
  void     SetFillFromGenerated(Bool_t flag) {fFillFromGenerated = flag;}
  Bool_t   GetFillFromGenerated() const {return fFillFromGenerated;}
  Bool_t   GetGeneratedValuesFromMCParticle(AliAODMCParticle* mcPart, TClonesArray* mcArray, Double_t* vectorMC) const;
  void     SetMinITSClusters(Int_t minITSClusters) {fMinITSClusters = minITSClusters;}
  Int_t    GetMinITSClusters() const {return fMinITSClusters;}
  Int_t    CheckOrigin(AliAODMCParticle* mcPart, TClonesArray* mcArray) const;

  void SetKeepD0fromB(Bool_t keepD0fromB){fKeepD0fromB=keepD0fromB;}
  Bool_t GetKeepD0fromB(){return fKeepD0fromB;}
  void SetKeepD0fromBOnly(Bool_t keepD0fromBOnly){ fKeepD0fromBOnly=keepD0fromBOnly; }
  Bool_t GetKeepD0fromBOnly(){ return fKeepD0fromBOnly;}
  void SetUseWeight(Bool_t useWeight){fUseWeight=useWeight;}
  Bool_t GetUseWeight() const {return fUseWeight;}
  Double_t GetWeight(Float_t pt);
  Double_t dNdptFit(Float_t pt, Double_t* par);
  void SetSign(Char_t isSign) {fSign = isSign;}
  Char_t GetSign() {return fSign;}

 protected:
  Int_t           fPDG;         //  PDG code of searched V0's
  AliCFManager   *fCFManager;   //  pointer to the CF manager
  TH1I *fHistEventsProcessed;   //! simple histo for monitoring the number of events processed
  THnSparse* fCorrelation;      //  response matrix for unfolding
  Int_t fCountMC;               //  MC particle found
  Int_t fCountAcc;              //  MC particle found that satisfy acceptance cuts
  Int_t fCountVertex;       //  Reco particle found that satisfy vertex constrained
  Int_t fCountRefit;        //  Reco particle found that satisfy kTPCrefit and kITSrefit
  Int_t fCountReco;             //  Reco particle found that satisfy cuts
  Int_t fCountRecoAcc;          //  Reco particle found that satisfy cuts in requested acceptance
  Int_t fCountRecoITSClusters;  //  Reco particle found that satisfy cuts in n. of ITS clusters
  Int_t fCountRecoPPR;          //  Reco particle found that satisfy cuts in PPR
  Int_t fCountRecoPID;          //  Reco particle found that satisfy cuts in PPR, and PID condition
  Int_t fEvents;                //  n. of events
  Bool_t fFillFromGenerated;    //  flag to indicate whether data container should be filled 
                                //  with generated values also for reconstructed particles
  Int_t fMinITSClusters;        //  min n. of ITS clusters for RecoDecay
  Bool_t fAcceptanceUnf;        //  flag for unfolding before or after cuts.

  Bool_t fKeepD0fromB;          // flag to consider also D0 coming from B
  Bool_t fKeepD0fromBOnly;      // flag to consider _only_ D0 coming from B
  AliRDHFCutsD0toKpi* fCuts;    // cuts
  Bool_t fUseWeight;            // flag to decide whether to use weights != 1 when filling the container or not
  Double_t fWeight;             // weight used to fill the container
  Char_t fSign;                 // flag to decide wheter to keep D0 only (0), D0bar only (1), or both D0 and D0bar (2)
  
  ClassDef(AliCFHeavyFlavourTaskMultiVarMultiStep,7); // class for HF corrections as a function of many variables
};

#endif
