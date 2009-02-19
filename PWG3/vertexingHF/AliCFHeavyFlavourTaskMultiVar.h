#ifndef ALICFHEAVYFLAVOURTASKMULTIVAR_H
#define ALICFHEAVYFLAVOURTASKMULTIVAR_H
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
// Class for HF corrections as a function of many variables
// Author : C. Zampolli, CERN
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

class AliCFHeavyFlavourTaskMultiVar : public AliAnalysisTaskSE {
  public:

  enum {
    kStepGenerated       = 0,
    kStepReconstructed   = 1
  };

  AliCFHeavyFlavourTaskMultiVar();
  AliCFHeavyFlavourTaskMultiVar(const Char_t* name);
  AliCFHeavyFlavourTaskMultiVar& operator= (const AliCFHeavyFlavourTaskMultiVar& c);
  AliCFHeavyFlavourTaskMultiVar(const AliCFHeavyFlavourTaskMultiVar& c);
  virtual ~AliCFHeavyFlavourTaskMultiVar();

  // ANALYSIS FRAMEWORK STUFF to loop on data and fill output objects
  void     UserCreateOutputObjects();
  void     UserExec(Option_t *option);
  void     Terminate(Option_t *);
  
  // CORRECTION FRAMEWORK RELATED FUNCTIONS
  void           SetCFManager(AliCFManager* io) {fCFManager = io;}   // global correction manager
  AliCFManager * GetCFManager()                 {return fCFManager;} // get corr manager

  void     SetPDG(Int_t code) {fPDG = code; }     // defines the PDG code of searched HF
  Int_t    IsMcVtx(AliAODRecoDecayHF2Prong* vtx); // checks if the AliAODv0 can be associated, returns mother label
  Int_t    GetVtxLabel(Int_t* labels);            // returns label of vertex given the daughter labels
  Double_t CosThetaStar(AliAODMCParticle* mcPart, AliAODMCParticle* mcPartDaughter0, AliAODMCParticle* mcPartDaughter1) const;  // returns cos(ThetaStar) of the D0 decay
  Double_t CT(AliAODMCParticle* mcPart, AliAODMCParticle* mcPartDaughter0, AliAODMCParticle* mcPartDaughter1) const;            // returns cT of the D0 decay
  void     SetFillFromGenerated(Bool_t flag) {fFillFromGenerated = flag;}
  Bool_t   GetFillFromGenerated() const {return fFillFromGenerated;}
  Bool_t   GetGeneratedValuesFromMCParticle(AliAODMCParticle* mcPart, TClonesArray* mcArray, Double_t* vectorMC) const;

 protected:
  Int_t           fPDG;         //  PDG code of searched V0's
  AliCFManager   *fCFManager  ; //  pointer to the CF manager
  TH1I *fHistEventsProcessed;   //! simple histo for monitoring the number of events processed
  Int_t fCountMC;               //  MC particle found
  Int_t fEvents;                //  n. of events
  Bool_t fFillFromGenerated;    //  flag to indicate whether data container should be filled 
                                //  with generated values also for reconstructed particles
  Int_t fSkipped;               //  n. of reco decays skipped after checking on the decay (MUST be D0->K+pi)

  
  ClassDef(AliCFHeavyFlavourTaskMultiVar,0); // class for HF corrections as a function of many variables
};

#endif
