#ifndef ALIEMCALJETUTILITYCONSTSUBTRACTOR_H
#define ALIEMCALJETUTILITYCONSTSUBTRACTOR_H

#include <TNamed.h>

#include "AliEmcalJetUtility.h"
#include "AliFJWrapper.h"

class TClonesArray;
class AliEmcalJetTask;
class AliEmcalJet;
class AliFJWrapper;
class AliRhoParameter;

class AliEmcalJetUtilityConstSubtractor : public AliEmcalJetUtility
{
 public:

  AliEmcalJetUtilityConstSubtractor();
  AliEmcalJetUtilityConstSubtractor(const char* name);
  AliEmcalJetUtilityConstSubtractor(const AliEmcalJetUtilityConstSubtractor &jet);
  AliEmcalJetUtilityConstSubtractor& operator=(const AliEmcalJetUtilityConstSubtractor &jet);
  ~AliEmcalJetUtilityConstSubtractor() {;}

  void                   SetRhoName(const char *n)           { fRhoName      = n         ; }
  void                   SetRhomName(const char *n)          { fRhomName     = n         ; }
  void                   SetUseExternalBkg(Bool_t b)         { fUseExternalBkg   = b     ; }

  void                   SetJetsSubName(const char *n)       { fJetsSubName      = n     ; }
  void                   SetParticlesSubName(const char *n)  { fParticlesSubName = n     ; }
  void                   SetAlpha(const Double_t a)            { fAlpha            = a     ; }
  void                   SetMaxDelR(const Double_t r)          { fMaxDelR          = r     ; }

  void Init();
  void InitEvent(AliFJWrapper& fjw);
  void Prepare(AliFJWrapper& fjw);
  void ProcessJet(AliEmcalJet* jet, Int_t ij, AliFJWrapper& fjw);
  void Terminate(AliFJWrapper& fjw);

 protected:

  TString                fJetsSubName;                        // name of subtracted jet collection
  TString                fParticlesSubName;                   // name of subtracted particle collection
  Bool_t                 fUseExternalBkg;                     // use external background for generic subtractor
  TString                fRhoName;                            // name of rho
  TString                fRhomName;                           // name of rhom
  Double_t               fRho;                                // pT background density
  Double_t               fRhom;                               // mT background density
  Double_t               fAlpha;                              // pT weight exponent applied in const sub
  Double_t               fMaxDelR;                            // Max distance between ghost and constituent pair in subtraction

  TClonesArray          *fJetsSub;                            //!subtracted jet collection
  TClonesArray          *fParticlesSub;                       //!subtracted particle collection
  AliRhoParameter       *fRhoParam;                           //!event rho
  AliRhoParameter       *fRhomParam;                          //!event rhom

  ClassDef(AliEmcalJetUtilityConstSubtractor, 2) // Emcal jet utility that implements the constituent subtractor form the fastjet contrib
};
#endif
